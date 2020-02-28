//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForKH.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSKSH(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSKSH_0_95(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSKSH_95_190(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSKSH_190_285(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSH_285_380(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSH_380_474(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSH_474_568(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSH_568_662(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSH_662_756(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSKSH_0_95(      CMemBlock2D<double>* primBuffer,
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
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_xxxxxx_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx); 

                auto tg_xxxxxx_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 1); 

                auto tg_xxxxxx_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 2); 

                auto tg_xxxxxx_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 3); 

                auto tg_xxxxxx_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 4); 

                auto tg_xxxxxx_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 5); 

                auto tg_xxxxxx_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 6); 

                auto tg_xxxxxx_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 7); 

                auto tg_xxxxxx_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 8); 

                auto tg_xxxxxx_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 9); 

                auto tg_xxxxxx_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 10); 

                auto tg_xxxxxx_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 11); 

                auto tg_xxxxxx_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 12); 

                auto tg_xxxxxx_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 13); 

                auto tg_xxxxxx_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 14); 

                auto tg_xxxxxx_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 15); 

                auto tg_xxxxxx_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 16); 

                auto tg_xxxxxx_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 17); 

                auto tg_xxxxxx_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 18); 

                auto tg_xxxxxx_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 19); 

                auto tg_xxxxxx_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 20); 

                auto tg_xxxxxy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 21); 

                auto tg_xxxxxy_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 22); 

                auto tg_xxxxxy_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 23); 

                auto tg_xxxxxy_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 24); 

                auto tg_xxxxxy_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 25); 

                auto tg_xxxxxy_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 26); 

                auto tg_xxxxxy_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 27); 

                auto tg_xxxxxy_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 28); 

                auto tg_xxxxxy_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 29); 

                auto tg_xxxxxy_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 30); 

                auto tg_xxxxxy_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 31); 

                auto tg_xxxxxy_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 32); 

                auto tg_xxxxxy_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 33); 

                auto tg_xxxxxy_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 34); 

                auto tg_xxxxxy_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 35); 

                auto tg_xxxxxy_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 36); 

                auto tg_xxxxxy_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 37); 

                auto tg_xxxxxy_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 38); 

                auto tg_xxxxxy_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 39); 

                auto tg_xxxxxy_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 40); 

                auto tg_xxxxxy_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 41); 

                auto tg_xxxxxz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 42); 

                auto tg_xxxxxz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 43); 

                auto tg_xxxxxz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 44); 

                auto tg_xxxxxz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 45); 

                auto tg_xxxxxz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 46); 

                auto tg_xxxxxz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 47); 

                auto tg_xxxxxz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 48); 

                auto tg_xxxxxz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 49); 

                auto tg_xxxxxz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 50); 

                auto tg_xxxxxz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 51); 

                auto tg_xxxxxz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 52); 

                auto tg_xxxxxz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 53); 

                auto tg_xxxxxz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 54); 

                auto tg_xxxxxz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 55); 

                auto tg_xxxxxz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 56); 

                auto tg_xxxxxz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 57); 

                auto tg_xxxxxz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 58); 

                auto tg_xxxxxz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 59); 

                auto tg_xxxxxz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 60); 

                auto tg_xxxxxz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 61); 

                auto tg_xxxxxz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 62); 

                auto tg_xxxxyy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 63); 

                auto tg_xxxxyy_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 64); 

                auto tg_xxxxyy_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 65); 

                auto tg_xxxxyy_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 66); 

                auto tg_xxxxyy_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 67); 

                auto tg_xxxxyy_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 68); 

                auto tg_xxxxyy_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 69); 

                auto tg_xxxxyy_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 70); 

                auto tg_xxxxyy_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 71); 

                auto tg_xxxxyy_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 72); 

                auto tg_xxxxyy_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 73); 

                auto tg_xxxxyy_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 74); 

                auto tg_xxxxyy_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 75); 

                auto tg_xxxxyy_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 76); 

                auto tg_xxxxyy_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 77); 

                auto tg_xxxxyy_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 78); 

                auto tg_xxxxyy_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 79); 

                auto tg_xxxxyy_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 80); 

                auto tg_xxxxyy_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 81); 

                auto tg_xxxxyy_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 82); 

                auto tg_xxxxyy_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 83); 

                auto tg_xxxxyz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 84); 

                auto tg_xxxxyz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 85); 

                auto tg_xxxxyz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 86); 

                auto tg_xxxxyz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 87); 

                auto tg_xxxxyz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 88); 

                auto tg_xxxxyz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 89); 

                auto tg_xxxxyz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 90); 

                auto tg_xxxxyz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 91); 

                auto tg_xxxxyz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 92); 

                auto tg_xxxxyz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 93); 

                auto tg_xxxxyz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 94); 

                auto tg_xxxxxx_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx); 

                auto tg_xxxxxx_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 1); 

                auto tg_xxxxxx_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 2); 

                auto tg_xxxxxx_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 3); 

                auto tg_xxxxxx_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 4); 

                auto tg_xxxxxx_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 5); 

                auto tg_xxxxxx_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 6); 

                auto tg_xxxxxx_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 7); 

                auto tg_xxxxxx_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 8); 

                auto tg_xxxxxx_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 9); 

                auto tg_xxxxxx_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 10); 

                auto tg_xxxxxx_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 11); 

                auto tg_xxxxxx_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 12); 

                auto tg_xxxxxx_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 13); 

                auto tg_xxxxxx_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 14); 

                auto tg_xxxxxx_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 15); 

                auto tg_xxxxxx_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 16); 

                auto tg_xxxxxx_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 17); 

                auto tg_xxxxxx_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 18); 

                auto tg_xxxxxx_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 19); 

                auto tg_xxxxxx_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 20); 

                auto tg_xxxxxy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 21); 

                auto tg_xxxxxy_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 22); 

                auto tg_xxxxxy_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 23); 

                auto tg_xxxxxy_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 24); 

                auto tg_xxxxxy_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 25); 

                auto tg_xxxxxy_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 26); 

                auto tg_xxxxxy_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 27); 

                auto tg_xxxxxy_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 28); 

                auto tg_xxxxxy_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 29); 

                auto tg_xxxxxy_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 30); 

                auto tg_xxxxxy_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 31); 

                auto tg_xxxxxy_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 32); 

                auto tg_xxxxxy_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 33); 

                auto tg_xxxxxy_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 34); 

                auto tg_xxxxxy_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 35); 

                auto tg_xxxxxy_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 36); 

                auto tg_xxxxxy_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 37); 

                auto tg_xxxxxy_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 38); 

                auto tg_xxxxxy_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 39); 

                auto tg_xxxxxy_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 40); 

                auto tg_xxxxxy_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 41); 

                auto tg_xxxxxz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 42); 

                auto tg_xxxxxz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 43); 

                auto tg_xxxxxz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 44); 

                auto tg_xxxxxz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 45); 

                auto tg_xxxxxz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 46); 

                auto tg_xxxxxz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 47); 

                auto tg_xxxxxz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 48); 

                auto tg_xxxxxz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 49); 

                auto tg_xxxxxz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 50); 

                auto tg_xxxxxz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 51); 

                auto tg_xxxxxz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 52); 

                auto tg_xxxxxz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 53); 

                auto tg_xxxxxz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 54); 

                auto tg_xxxxxz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 55); 

                auto tg_xxxxxz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 56); 

                auto tg_xxxxxz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 57); 

                auto tg_xxxxxz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 58); 

                auto tg_xxxxxz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 59); 

                auto tg_xxxxxz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 60); 

                auto tg_xxxxxz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 61); 

                auto tg_xxxxxz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 62); 

                auto tg_xxxxyy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 63); 

                auto tg_xxxxyy_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 64); 

                auto tg_xxxxyy_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 65); 

                auto tg_xxxxyy_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 66); 

                auto tg_xxxxyy_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 67); 

                auto tg_xxxxyy_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 68); 

                auto tg_xxxxyy_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 69); 

                auto tg_xxxxyy_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 70); 

                auto tg_xxxxyy_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 71); 

                auto tg_xxxxyy_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 72); 

                auto tg_xxxxyy_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 73); 

                auto tg_xxxxyy_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 74); 

                auto tg_xxxxyy_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 75); 

                auto tg_xxxxyy_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 76); 

                auto tg_xxxxyy_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 77); 

                auto tg_xxxxyy_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 78); 

                auto tg_xxxxyy_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 79); 

                auto tg_xxxxyy_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 80); 

                auto tg_xxxxyy_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 81); 

                auto tg_xxxxyy_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 82); 

                auto tg_xxxxyy_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 83); 

                auto tg_xxxxyz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 84); 

                auto tg_xxxxyz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 85); 

                auto tg_xxxxyz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 86); 

                auto tg_xxxxyz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 87); 

                auto tg_xxxxyz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 88); 

                auto tg_xxxxyz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 89); 

                auto tg_xxxxyz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 90); 

                auto tg_xxxxyz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 91); 

                auto tg_xxxxyz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 92); 

                auto tg_xxxxyz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 93); 

                auto tg_xxxxyz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 94); 

                auto tg_xxxxx_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx); 

                auto tg_xxxxx_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 1); 

                auto tg_xxxxx_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 2); 

                auto tg_xxxxx_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 3); 

                auto tg_xxxxx_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 4); 

                auto tg_xxxxx_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 5); 

                auto tg_xxxxx_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 6); 

                auto tg_xxxxx_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 7); 

                auto tg_xxxxx_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 8); 

                auto tg_xxxxx_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 9); 

                auto tg_xxxxx_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 10); 

                auto tg_xxxxx_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 11); 

                auto tg_xxxxx_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 12); 

                auto tg_xxxxx_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 13); 

                auto tg_xxxxx_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 14); 

                auto tg_xxxxx_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 15); 

                auto tg_xxxxx_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 16); 

                auto tg_xxxxx_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 17); 

                auto tg_xxxxx_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 18); 

                auto tg_xxxxx_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 19); 

                auto tg_xxxxx_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 20); 

                auto tg_xxxxy_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 21); 

                auto tg_xxxxy_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 22); 

                auto tg_xxxxy_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 23); 

                auto tg_xxxxy_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 24); 

                auto tg_xxxxy_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 25); 

                auto tg_xxxxy_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 26); 

                auto tg_xxxxy_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 27); 

                auto tg_xxxxy_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 28); 

                auto tg_xxxxy_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 29); 

                auto tg_xxxxy_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 30); 

                auto tg_xxxxy_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 31); 

                auto tg_xxxxy_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 32); 

                auto tg_xxxxy_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 33); 

                auto tg_xxxxy_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 34); 

                auto tg_xxxxy_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 35); 

                auto tg_xxxxy_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 36); 

                auto tg_xxxxy_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 37); 

                auto tg_xxxxy_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 38); 

                auto tg_xxxxy_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 39); 

                auto tg_xxxxy_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 40); 

                auto tg_xxxxy_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 41); 

                auto tg_xxxxz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 42); 

                auto tg_xxxxz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 43); 

                auto tg_xxxxz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 44); 

                auto tg_xxxxz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 45); 

                auto tg_xxxxz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 46); 

                auto tg_xxxxz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 47); 

                auto tg_xxxxz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 48); 

                auto tg_xxxxz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 49); 

                auto tg_xxxxz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 50); 

                auto tg_xxxxz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 51); 

                auto tg_xxxxz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 52); 

                auto tg_xxxxz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 53); 

                auto tg_xxxxz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 54); 

                auto tg_xxxxz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 55); 

                auto tg_xxxxz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 56); 

                auto tg_xxxxz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 57); 

                auto tg_xxxxz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 58); 

                auto tg_xxxxz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 59); 

                auto tg_xxxxz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 60); 

                auto tg_xxxxz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 61); 

                auto tg_xxxxz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 62); 

                auto tg_xxxyy_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 63); 

                auto tg_xxxyy_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 64); 

                auto tg_xxxyy_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 65); 

                auto tg_xxxyy_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 66); 

                auto tg_xxxyy_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 67); 

                auto tg_xxxyy_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 68); 

                auto tg_xxxyy_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 69); 

                auto tg_xxxyy_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 70); 

                auto tg_xxxyy_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 71); 

                auto tg_xxxyy_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 72); 

                auto tg_xxxyy_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 73); 

                auto tg_xxxyy_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 74); 

                auto tg_xxxyy_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 75); 

                auto tg_xxxyy_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 76); 

                auto tg_xxxyy_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 77); 

                auto tg_xxxyy_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 78); 

                auto tg_xxxyy_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 79); 

                auto tg_xxxyy_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 80); 

                auto tg_xxxyy_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 81); 

                auto tg_xxxyy_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 82); 

                auto tg_xxxyy_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 83); 

                auto tg_xxxyz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 84); 

                auto tg_xxxyz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 85); 

                auto tg_xxxyz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 86); 

                auto tg_xxxyz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 87); 

                auto tg_xxxyz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 88); 

                auto tg_xxxyz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 89); 

                auto tg_xxxyz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 90); 

                auto tg_xxxyz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 91); 

                auto tg_xxxyz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 92); 

                auto tg_xxxyz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 93); 

                auto tg_xxxyz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 94); 

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

                auto tg_xxxxxx_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx); 

                auto tg_xxxxxx_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 1); 

                auto tg_xxxxxx_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 2); 

                auto tg_xxxxxx_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 3); 

                auto tg_xxxxxx_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 4); 

                auto tg_xxxxxx_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 5); 

                auto tg_xxxxxx_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 6); 

                auto tg_xxxxxx_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 7); 

                auto tg_xxxxxx_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 8); 

                auto tg_xxxxxx_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 9); 

                auto tg_xxxxxx_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 10); 

                auto tg_xxxxxx_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 11); 

                auto tg_xxxxxx_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 12); 

                auto tg_xxxxxx_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 13); 

                auto tg_xxxxxx_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 14); 

                auto tg_xxxxxy_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 15); 

                auto tg_xxxxxy_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 16); 

                auto tg_xxxxxy_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 17); 

                auto tg_xxxxxy_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 18); 

                auto tg_xxxxxy_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 19); 

                auto tg_xxxxxy_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 20); 

                auto tg_xxxxxy_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 21); 

                auto tg_xxxxxy_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 22); 

                auto tg_xxxxxy_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 23); 

                auto tg_xxxxxy_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 24); 

                auto tg_xxxxxy_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 25); 

                auto tg_xxxxxy_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 26); 

                auto tg_xxxxxy_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 27); 

                auto tg_xxxxxy_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 28); 

                auto tg_xxxxxy_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 29); 

                auto tg_xxxxxz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 30); 

                auto tg_xxxxxz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 31); 

                auto tg_xxxxxz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 32); 

                auto tg_xxxxxz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 33); 

                auto tg_xxxxxz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 34); 

                auto tg_xxxxxz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 35); 

                auto tg_xxxxxz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 36); 

                auto tg_xxxxxz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 37); 

                auto tg_xxxxxz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 38); 

                auto tg_xxxxxz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 39); 

                auto tg_xxxxxz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 40); 

                auto tg_xxxxxz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 41); 

                auto tg_xxxxxz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 42); 

                auto tg_xxxxxz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 43); 

                auto tg_xxxxxz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 44); 

                auto tg_xxxxyy_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 45); 

                auto tg_xxxxyy_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 46); 

                auto tg_xxxxyy_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 47); 

                auto tg_xxxxyy_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 48); 

                auto tg_xxxxyy_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 49); 

                auto tg_xxxxyy_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 50); 

                auto tg_xxxxyy_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 51); 

                auto tg_xxxxyy_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 52); 

                auto tg_xxxxyy_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 53); 

                auto tg_xxxxyy_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 54); 

                auto tg_xxxxyy_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 55); 

                auto tg_xxxxyy_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 56); 

                auto tg_xxxxyy_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 57); 

                auto tg_xxxxyy_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 58); 

                auto tg_xxxxyy_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 59); 

                auto tg_xxxxyz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 60); 

                auto tg_xxxxyz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 61); 

                auto tg_xxxxyz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 62); 

                auto tg_xxxxyz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 63); 

                auto tg_xxxxyz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 64); 

                auto tg_xxxxyz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 65); 

                auto tg_xxxxyz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 66); 

                auto tg_xxxxyz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 67); 

                auto tg_xxxxyz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 68); 

                auto tg_xxxxyz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 69); 

                auto tg_xxxxyz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 70); 

                // set up pointers to integrals

                auto tg_xxxxxxx_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx); 

                auto tg_xxxxxxx_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 1); 

                auto tg_xxxxxxx_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 2); 

                auto tg_xxxxxxx_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 3); 

                auto tg_xxxxxxx_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 4); 

                auto tg_xxxxxxx_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 5); 

                auto tg_xxxxxxx_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 6); 

                auto tg_xxxxxxx_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 7); 

                auto tg_xxxxxxx_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 8); 

                auto tg_xxxxxxx_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 9); 

                auto tg_xxxxxxx_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 10); 

                auto tg_xxxxxxx_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 11); 

                auto tg_xxxxxxx_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 12); 

                auto tg_xxxxxxx_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 13); 

                auto tg_xxxxxxx_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 14); 

                auto tg_xxxxxxx_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 15); 

                auto tg_xxxxxxx_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 16); 

                auto tg_xxxxxxx_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 17); 

                auto tg_xxxxxxx_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 18); 

                auto tg_xxxxxxx_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 19); 

                auto tg_xxxxxxx_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 20); 

                auto tg_xxxxxxy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 21); 

                auto tg_xxxxxxy_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 22); 

                auto tg_xxxxxxy_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 23); 

                auto tg_xxxxxxy_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 24); 

                auto tg_xxxxxxy_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 25); 

                auto tg_xxxxxxy_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 26); 

                auto tg_xxxxxxy_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 27); 

                auto tg_xxxxxxy_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 28); 

                auto tg_xxxxxxy_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 29); 

                auto tg_xxxxxxy_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 30); 

                auto tg_xxxxxxy_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 31); 

                auto tg_xxxxxxy_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 32); 

                auto tg_xxxxxxy_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 33); 

                auto tg_xxxxxxy_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 34); 

                auto tg_xxxxxxy_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 35); 

                auto tg_xxxxxxy_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 36); 

                auto tg_xxxxxxy_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 37); 

                auto tg_xxxxxxy_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 38); 

                auto tg_xxxxxxy_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 39); 

                auto tg_xxxxxxy_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 40); 

                auto tg_xxxxxxy_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 41); 

                auto tg_xxxxxxz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 42); 

                auto tg_xxxxxxz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 43); 

                auto tg_xxxxxxz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 44); 

                auto tg_xxxxxxz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 45); 

                auto tg_xxxxxxz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 46); 

                auto tg_xxxxxxz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 47); 

                auto tg_xxxxxxz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 48); 

                auto tg_xxxxxxz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 49); 

                auto tg_xxxxxxz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 50); 

                auto tg_xxxxxxz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 51); 

                auto tg_xxxxxxz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 52); 

                auto tg_xxxxxxz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 53); 

                auto tg_xxxxxxz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 54); 

                auto tg_xxxxxxz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 55); 

                auto tg_xxxxxxz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 56); 

                auto tg_xxxxxxz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 57); 

                auto tg_xxxxxxz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 58); 

                auto tg_xxxxxxz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 59); 

                auto tg_xxxxxxz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 60); 

                auto tg_xxxxxxz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 61); 

                auto tg_xxxxxxz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 62); 

                auto tg_xxxxxyy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 63); 

                auto tg_xxxxxyy_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 64); 

                auto tg_xxxxxyy_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 65); 

                auto tg_xxxxxyy_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 66); 

                auto tg_xxxxxyy_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 67); 

                auto tg_xxxxxyy_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 68); 

                auto tg_xxxxxyy_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 69); 

                auto tg_xxxxxyy_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 70); 

                auto tg_xxxxxyy_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 71); 

                auto tg_xxxxxyy_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 72); 

                auto tg_xxxxxyy_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 73); 

                auto tg_xxxxxyy_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 74); 

                auto tg_xxxxxyy_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 75); 

                auto tg_xxxxxyy_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 76); 

                auto tg_xxxxxyy_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 77); 

                auto tg_xxxxxyy_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 78); 

                auto tg_xxxxxyy_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 79); 

                auto tg_xxxxxyy_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 80); 

                auto tg_xxxxxyy_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 81); 

                auto tg_xxxxxyy_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 82); 

                auto tg_xxxxxyy_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 83); 

                auto tg_xxxxxyz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 84); 

                auto tg_xxxxxyz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 85); 

                auto tg_xxxxxyz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 86); 

                auto tg_xxxxxyz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 87); 

                auto tg_xxxxxyz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 88); 

                auto tg_xxxxxyz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 89); 

                auto tg_xxxxxyz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 90); 

                auto tg_xxxxxyz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 91); 

                auto tg_xxxxxyz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 92); 

                auto tg_xxxxxyz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 93); 

                auto tg_xxxxxyz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 94); 

                // Batch of Integrals (0,95)

                #pragma omp simd aligned(fxn, fza, tg_xxxxx_xxxxx_0, tg_xxxxx_xxxxx_1, tg_xxxxx_xxxxy_0, \
                                         tg_xxxxx_xxxxy_1, tg_xxxxx_xxxxz_0, tg_xxxxx_xxxxz_1, tg_xxxxx_xxxyy_0, \
                                         tg_xxxxx_xxxyy_1, tg_xxxxx_xxxyz_0, tg_xxxxx_xxxyz_1, tg_xxxxx_xxxzz_0, \
                                         tg_xxxxx_xxxzz_1, tg_xxxxx_xxyyy_0, tg_xxxxx_xxyyy_1, tg_xxxxx_xxyyz_0, \
                                         tg_xxxxx_xxyyz_1, tg_xxxxx_xxyzz_0, tg_xxxxx_xxyzz_1, tg_xxxxx_xxzzz_0, \
                                         tg_xxxxx_xxzzz_1, tg_xxxxx_xyyyy_0, tg_xxxxx_xyyyy_1, tg_xxxxx_xyyyz_0, \
                                         tg_xxxxx_xyyyz_1, tg_xxxxx_xyyzz_0, tg_xxxxx_xyyzz_1, tg_xxxxx_xyzzz_0, \
                                         tg_xxxxx_xyzzz_1, tg_xxxxx_xzzzz_0, tg_xxxxx_xzzzz_1, tg_xxxxx_yyyyy_0, \
                                         tg_xxxxx_yyyyy_1, tg_xxxxx_yyyyz_0, tg_xxxxx_yyyyz_1, tg_xxxxx_yyyzz_0, \
                                         tg_xxxxx_yyyzz_1, tg_xxxxx_yyzzz_0, tg_xxxxx_yyzzz_1, tg_xxxxx_yzzzz_0, \
                                         tg_xxxxx_yzzzz_1, tg_xxxxx_zzzzz_0, tg_xxxxx_zzzzz_1, tg_xxxxxx_xxxx_1, \
                                         tg_xxxxxx_xxxxx_0, tg_xxxxxx_xxxxx_1, tg_xxxxxx_xxxxy_0, tg_xxxxxx_xxxxy_1, \
                                         tg_xxxxxx_xxxxz_0, tg_xxxxxx_xxxxz_1, tg_xxxxxx_xxxy_1, tg_xxxxxx_xxxyy_0, \
                                         tg_xxxxxx_xxxyy_1, tg_xxxxxx_xxxyz_0, tg_xxxxxx_xxxyz_1, tg_xxxxxx_xxxz_1, \
                                         tg_xxxxxx_xxxzz_0, tg_xxxxxx_xxxzz_1, tg_xxxxxx_xxyy_1, tg_xxxxxx_xxyyy_0, \
                                         tg_xxxxxx_xxyyy_1, tg_xxxxxx_xxyyz_0, tg_xxxxxx_xxyyz_1, tg_xxxxxx_xxyz_1, \
                                         tg_xxxxxx_xxyzz_0, tg_xxxxxx_xxyzz_1, tg_xxxxxx_xxzz_1, tg_xxxxxx_xxzzz_0, \
                                         tg_xxxxxx_xxzzz_1, tg_xxxxxx_xyyy_1, tg_xxxxxx_xyyyy_0, tg_xxxxxx_xyyyy_1, \
                                         tg_xxxxxx_xyyyz_0, tg_xxxxxx_xyyyz_1, tg_xxxxxx_xyyz_1, tg_xxxxxx_xyyzz_0, \
                                         tg_xxxxxx_xyyzz_1, tg_xxxxxx_xyzz_1, tg_xxxxxx_xyzzz_0, tg_xxxxxx_xyzzz_1, \
                                         tg_xxxxxx_xzzz_1, tg_xxxxxx_xzzzz_0, tg_xxxxxx_xzzzz_1, tg_xxxxxx_yyyy_1, \
                                         tg_xxxxxx_yyyyy_0, tg_xxxxxx_yyyyy_1, tg_xxxxxx_yyyyz_0, tg_xxxxxx_yyyyz_1, \
                                         tg_xxxxxx_yyyz_1, tg_xxxxxx_yyyzz_0, tg_xxxxxx_yyyzz_1, tg_xxxxxx_yyzz_1, \
                                         tg_xxxxxx_yyzzz_0, tg_xxxxxx_yyzzz_1, tg_xxxxxx_yzzz_1, tg_xxxxxx_yzzzz_0, \
                                         tg_xxxxxx_yzzzz_1, tg_xxxxxx_zzzz_1, tg_xxxxxx_zzzzz_0, tg_xxxxxx_zzzzz_1, \
                                         tg_xxxxxxx_xxxxx_0, tg_xxxxxxx_xxxxy_0, tg_xxxxxxx_xxxxz_0, tg_xxxxxxx_xxxyy_0, \
                                         tg_xxxxxxx_xxxyz_0, tg_xxxxxxx_xxxzz_0, tg_xxxxxxx_xxyyy_0, tg_xxxxxxx_xxyyz_0, \
                                         tg_xxxxxxx_xxyzz_0, tg_xxxxxxx_xxzzz_0, tg_xxxxxxx_xyyyy_0, tg_xxxxxxx_xyyyz_0, \
                                         tg_xxxxxxx_xyyzz_0, tg_xxxxxxx_xyzzz_0, tg_xxxxxxx_xzzzz_0, tg_xxxxxxx_yyyyy_0, \
                                         tg_xxxxxxx_yyyyz_0, tg_xxxxxxx_yyyzz_0, tg_xxxxxxx_yyzzz_0, tg_xxxxxxx_yzzzz_0, \
                                         tg_xxxxxxx_zzzzz_0, tg_xxxxxxy_xxxxx_0, tg_xxxxxxy_xxxxy_0, tg_xxxxxxy_xxxxz_0, \
                                         tg_xxxxxxy_xxxyy_0, tg_xxxxxxy_xxxyz_0, tg_xxxxxxy_xxxzz_0, tg_xxxxxxy_xxyyy_0, \
                                         tg_xxxxxxy_xxyyz_0, tg_xxxxxxy_xxyzz_0, tg_xxxxxxy_xxzzz_0, tg_xxxxxxy_xyyyy_0, \
                                         tg_xxxxxxy_xyyyz_0, tg_xxxxxxy_xyyzz_0, tg_xxxxxxy_xyzzz_0, tg_xxxxxxy_xzzzz_0, \
                                         tg_xxxxxxy_yyyyy_0, tg_xxxxxxy_yyyyz_0, tg_xxxxxxy_yyyzz_0, tg_xxxxxxy_yyzzz_0, \
                                         tg_xxxxxxy_yzzzz_0, tg_xxxxxxy_zzzzz_0, tg_xxxxxxz_xxxxx_0, tg_xxxxxxz_xxxxy_0, \
                                         tg_xxxxxxz_xxxxz_0, tg_xxxxxxz_xxxyy_0, tg_xxxxxxz_xxxyz_0, tg_xxxxxxz_xxxzz_0, \
                                         tg_xxxxxxz_xxyyy_0, tg_xxxxxxz_xxyyz_0, tg_xxxxxxz_xxyzz_0, tg_xxxxxxz_xxzzz_0, \
                                         tg_xxxxxxz_xyyyy_0, tg_xxxxxxz_xyyyz_0, tg_xxxxxxz_xyyzz_0, tg_xxxxxxz_xyzzz_0, \
                                         tg_xxxxxxz_xzzzz_0, tg_xxxxxxz_yyyyy_0, tg_xxxxxxz_yyyyz_0, tg_xxxxxxz_yyyzz_0, \
                                         tg_xxxxxxz_yyzzz_0, tg_xxxxxxz_yzzzz_0, tg_xxxxxxz_zzzzz_0, tg_xxxxxy_xxxx_1, \
                                         tg_xxxxxy_xxxxx_0, tg_xxxxxy_xxxxx_1, tg_xxxxxy_xxxxy_0, tg_xxxxxy_xxxxy_1, \
                                         tg_xxxxxy_xxxxz_0, tg_xxxxxy_xxxxz_1, tg_xxxxxy_xxxy_1, tg_xxxxxy_xxxyy_0, \
                                         tg_xxxxxy_xxxyy_1, tg_xxxxxy_xxxyz_0, tg_xxxxxy_xxxyz_1, tg_xxxxxy_xxxz_1, \
                                         tg_xxxxxy_xxxzz_0, tg_xxxxxy_xxxzz_1, tg_xxxxxy_xxyy_1, tg_xxxxxy_xxyyy_0, \
                                         tg_xxxxxy_xxyyy_1, tg_xxxxxy_xxyyz_0, tg_xxxxxy_xxyyz_1, tg_xxxxxy_xxyz_1, \
                                         tg_xxxxxy_xxyzz_0, tg_xxxxxy_xxyzz_1, tg_xxxxxy_xxzz_1, tg_xxxxxy_xxzzz_0, \
                                         tg_xxxxxy_xxzzz_1, tg_xxxxxy_xyyy_1, tg_xxxxxy_xyyyy_0, tg_xxxxxy_xyyyy_1, \
                                         tg_xxxxxy_xyyyz_0, tg_xxxxxy_xyyyz_1, tg_xxxxxy_xyyz_1, tg_xxxxxy_xyyzz_0, \
                                         tg_xxxxxy_xyyzz_1, tg_xxxxxy_xyzz_1, tg_xxxxxy_xyzzz_0, tg_xxxxxy_xyzzz_1, \
                                         tg_xxxxxy_xzzz_1, tg_xxxxxy_xzzzz_0, tg_xxxxxy_xzzzz_1, tg_xxxxxy_yyyy_1, \
                                         tg_xxxxxy_yyyyy_0, tg_xxxxxy_yyyyy_1, tg_xxxxxy_yyyyz_0, tg_xxxxxy_yyyyz_1, \
                                         tg_xxxxxy_yyyz_1, tg_xxxxxy_yyyzz_0, tg_xxxxxy_yyyzz_1, tg_xxxxxy_yyzz_1, \
                                         tg_xxxxxy_yyzzz_0, tg_xxxxxy_yyzzz_1, tg_xxxxxy_yzzz_1, tg_xxxxxy_yzzzz_0, \
                                         tg_xxxxxy_yzzzz_1, tg_xxxxxy_zzzz_1, tg_xxxxxy_zzzzz_0, tg_xxxxxy_zzzzz_1, \
                                         tg_xxxxxyy_xxxxx_0, tg_xxxxxyy_xxxxy_0, tg_xxxxxyy_xxxxz_0, tg_xxxxxyy_xxxyy_0, \
                                         tg_xxxxxyy_xxxyz_0, tg_xxxxxyy_xxxzz_0, tg_xxxxxyy_xxyyy_0, tg_xxxxxyy_xxyyz_0, \
                                         tg_xxxxxyy_xxyzz_0, tg_xxxxxyy_xxzzz_0, tg_xxxxxyy_xyyyy_0, tg_xxxxxyy_xyyyz_0, \
                                         tg_xxxxxyy_xyyzz_0, tg_xxxxxyy_xyzzz_0, tg_xxxxxyy_xzzzz_0, tg_xxxxxyy_yyyyy_0, \
                                         tg_xxxxxyy_yyyyz_0, tg_xxxxxyy_yyyzz_0, tg_xxxxxyy_yyzzz_0, tg_xxxxxyy_yzzzz_0, \
                                         tg_xxxxxyy_zzzzz_0, tg_xxxxxyz_xxxxx_0, tg_xxxxxyz_xxxxy_0, tg_xxxxxyz_xxxxz_0, \
                                         tg_xxxxxyz_xxxyy_0, tg_xxxxxyz_xxxyz_0, tg_xxxxxyz_xxxzz_0, tg_xxxxxyz_xxyyy_0, \
                                         tg_xxxxxyz_xxyyz_0, tg_xxxxxyz_xxyzz_0, tg_xxxxxyz_xxzzz_0, tg_xxxxxyz_xyyyy_0, \
                                         tg_xxxxxz_xxxx_1, tg_xxxxxz_xxxxx_0, tg_xxxxxz_xxxxx_1, tg_xxxxxz_xxxxy_0, \
                                         tg_xxxxxz_xxxxy_1, tg_xxxxxz_xxxxz_0, tg_xxxxxz_xxxxz_1, tg_xxxxxz_xxxy_1, \
                                         tg_xxxxxz_xxxyy_0, tg_xxxxxz_xxxyy_1, tg_xxxxxz_xxxyz_0, tg_xxxxxz_xxxyz_1, \
                                         tg_xxxxxz_xxxz_1, tg_xxxxxz_xxxzz_0, tg_xxxxxz_xxxzz_1, tg_xxxxxz_xxyy_1, \
                                         tg_xxxxxz_xxyyy_0, tg_xxxxxz_xxyyy_1, tg_xxxxxz_xxyyz_0, tg_xxxxxz_xxyyz_1, \
                                         tg_xxxxxz_xxyz_1, tg_xxxxxz_xxyzz_0, tg_xxxxxz_xxyzz_1, tg_xxxxxz_xxzz_1, \
                                         tg_xxxxxz_xxzzz_0, tg_xxxxxz_xxzzz_1, tg_xxxxxz_xyyy_1, tg_xxxxxz_xyyyy_0, \
                                         tg_xxxxxz_xyyyy_1, tg_xxxxxz_xyyyz_0, tg_xxxxxz_xyyyz_1, tg_xxxxxz_xyyz_1, \
                                         tg_xxxxxz_xyyzz_0, tg_xxxxxz_xyyzz_1, tg_xxxxxz_xyzz_1, tg_xxxxxz_xyzzz_0, \
                                         tg_xxxxxz_xyzzz_1, tg_xxxxxz_xzzz_1, tg_xxxxxz_xzzzz_0, tg_xxxxxz_xzzzz_1, \
                                         tg_xxxxxz_yyyy_1, tg_xxxxxz_yyyyy_0, tg_xxxxxz_yyyyy_1, tg_xxxxxz_yyyyz_0, \
                                         tg_xxxxxz_yyyyz_1, tg_xxxxxz_yyyz_1, tg_xxxxxz_yyyzz_0, tg_xxxxxz_yyyzz_1, \
                                         tg_xxxxxz_yyzz_1, tg_xxxxxz_yyzzz_0, tg_xxxxxz_yyzzz_1, tg_xxxxxz_yzzz_1, \
                                         tg_xxxxxz_yzzzz_0, tg_xxxxxz_yzzzz_1, tg_xxxxxz_zzzz_1, tg_xxxxxz_zzzzz_0, \
                                         tg_xxxxxz_zzzzz_1, tg_xxxxy_xxxxx_0, tg_xxxxy_xxxxx_1, tg_xxxxy_xxxxy_0, \
                                         tg_xxxxy_xxxxy_1, tg_xxxxy_xxxxz_0, tg_xxxxy_xxxxz_1, tg_xxxxy_xxxyy_0, \
                                         tg_xxxxy_xxxyy_1, tg_xxxxy_xxxyz_0, tg_xxxxy_xxxyz_1, tg_xxxxy_xxxzz_0, \
                                         tg_xxxxy_xxxzz_1, tg_xxxxy_xxyyy_0, tg_xxxxy_xxyyy_1, tg_xxxxy_xxyyz_0, \
                                         tg_xxxxy_xxyyz_1, tg_xxxxy_xxyzz_0, tg_xxxxy_xxyzz_1, tg_xxxxy_xxzzz_0, \
                                         tg_xxxxy_xxzzz_1, tg_xxxxy_xyyyy_0, tg_xxxxy_xyyyy_1, tg_xxxxy_xyyyz_0, \
                                         tg_xxxxy_xyyyz_1, tg_xxxxy_xyyzz_0, tg_xxxxy_xyyzz_1, tg_xxxxy_xyzzz_0, \
                                         tg_xxxxy_xyzzz_1, tg_xxxxy_xzzzz_0, tg_xxxxy_xzzzz_1, tg_xxxxy_yyyyy_0, \
                                         tg_xxxxy_yyyyy_1, tg_xxxxy_yyyyz_0, tg_xxxxy_yyyyz_1, tg_xxxxy_yyyzz_0, \
                                         tg_xxxxy_yyyzz_1, tg_xxxxy_yyzzz_0, tg_xxxxy_yyzzz_1, tg_xxxxy_yzzzz_0, \
                                         tg_xxxxy_yzzzz_1, tg_xxxxy_zzzzz_0, tg_xxxxy_zzzzz_1, tg_xxxxyy_xxxx_1, \
                                         tg_xxxxyy_xxxxx_0, tg_xxxxyy_xxxxx_1, tg_xxxxyy_xxxxy_0, tg_xxxxyy_xxxxy_1, \
                                         tg_xxxxyy_xxxxz_0, tg_xxxxyy_xxxxz_1, tg_xxxxyy_xxxy_1, tg_xxxxyy_xxxyy_0, \
                                         tg_xxxxyy_xxxyy_1, tg_xxxxyy_xxxyz_0, tg_xxxxyy_xxxyz_1, tg_xxxxyy_xxxz_1, \
                                         tg_xxxxyy_xxxzz_0, tg_xxxxyy_xxxzz_1, tg_xxxxyy_xxyy_1, tg_xxxxyy_xxyyy_0, \
                                         tg_xxxxyy_xxyyy_1, tg_xxxxyy_xxyyz_0, tg_xxxxyy_xxyyz_1, tg_xxxxyy_xxyz_1, \
                                         tg_xxxxyy_xxyzz_0, tg_xxxxyy_xxyzz_1, tg_xxxxyy_xxzz_1, tg_xxxxyy_xxzzz_0, \
                                         tg_xxxxyy_xxzzz_1, tg_xxxxyy_xyyy_1, tg_xxxxyy_xyyyy_0, tg_xxxxyy_xyyyy_1, \
                                         tg_xxxxyy_xyyyz_0, tg_xxxxyy_xyyyz_1, tg_xxxxyy_xyyz_1, tg_xxxxyy_xyyzz_0, \
                                         tg_xxxxyy_xyyzz_1, tg_xxxxyy_xyzz_1, tg_xxxxyy_xyzzz_0, tg_xxxxyy_xyzzz_1, \
                                         tg_xxxxyy_xzzz_1, tg_xxxxyy_xzzzz_0, tg_xxxxyy_xzzzz_1, tg_xxxxyy_yyyy_1, \
                                         tg_xxxxyy_yyyyy_0, tg_xxxxyy_yyyyy_1, tg_xxxxyy_yyyyz_0, tg_xxxxyy_yyyyz_1, \
                                         tg_xxxxyy_yyyz_1, tg_xxxxyy_yyyzz_0, tg_xxxxyy_yyyzz_1, tg_xxxxyy_yyzz_1, \
                                         tg_xxxxyy_yyzzz_0, tg_xxxxyy_yyzzz_1, tg_xxxxyy_yzzz_1, tg_xxxxyy_yzzzz_0, \
                                         tg_xxxxyy_yzzzz_1, tg_xxxxyy_zzzz_1, tg_xxxxyy_zzzzz_0, tg_xxxxyy_zzzzz_1, \
                                         tg_xxxxyz_xxxx_1, tg_xxxxyz_xxxxx_0, tg_xxxxyz_xxxxx_1, tg_xxxxyz_xxxxy_0, \
                                         tg_xxxxyz_xxxxy_1, tg_xxxxyz_xxxxz_0, tg_xxxxyz_xxxxz_1, tg_xxxxyz_xxxy_1, \
                                         tg_xxxxyz_xxxyy_0, tg_xxxxyz_xxxyy_1, tg_xxxxyz_xxxyz_0, tg_xxxxyz_xxxyz_1, \
                                         tg_xxxxyz_xxxz_1, tg_xxxxyz_xxxzz_0, tg_xxxxyz_xxxzz_1, tg_xxxxyz_xxyy_1, \
                                         tg_xxxxyz_xxyyy_0, tg_xxxxyz_xxyyy_1, tg_xxxxyz_xxyyz_0, tg_xxxxyz_xxyyz_1, \
                                         tg_xxxxyz_xxyz_1, tg_xxxxyz_xxyzz_0, tg_xxxxyz_xxyzz_1, tg_xxxxyz_xxzz_1, \
                                         tg_xxxxyz_xxzzz_0, tg_xxxxyz_xxzzz_1, tg_xxxxyz_xyyy_1, tg_xxxxyz_xyyyy_0, \
                                         tg_xxxxyz_xyyyy_1, tg_xxxxyz_xyyz_1, tg_xxxxyz_xyzz_1, tg_xxxxyz_xzzz_1, \
                                         tg_xxxxyz_yyyy_1, tg_xxxxz_xxxxx_0, tg_xxxxz_xxxxx_1, tg_xxxxz_xxxxy_0, \
                                         tg_xxxxz_xxxxy_1, tg_xxxxz_xxxxz_0, tg_xxxxz_xxxxz_1, tg_xxxxz_xxxyy_0, \
                                         tg_xxxxz_xxxyy_1, tg_xxxxz_xxxyz_0, tg_xxxxz_xxxyz_1, tg_xxxxz_xxxzz_0, \
                                         tg_xxxxz_xxxzz_1, tg_xxxxz_xxyyy_0, tg_xxxxz_xxyyy_1, tg_xxxxz_xxyyz_0, \
                                         tg_xxxxz_xxyyz_1, tg_xxxxz_xxyzz_0, tg_xxxxz_xxyzz_1, tg_xxxxz_xxzzz_0, \
                                         tg_xxxxz_xxzzz_1, tg_xxxxz_xyyyy_0, tg_xxxxz_xyyyy_1, tg_xxxxz_xyyyz_0, \
                                         tg_xxxxz_xyyyz_1, tg_xxxxz_xyyzz_0, tg_xxxxz_xyyzz_1, tg_xxxxz_xyzzz_0, \
                                         tg_xxxxz_xyzzz_1, tg_xxxxz_xzzzz_0, tg_xxxxz_xzzzz_1, tg_xxxxz_yyyyy_0, \
                                         tg_xxxxz_yyyyy_1, tg_xxxxz_yyyyz_0, tg_xxxxz_yyyyz_1, tg_xxxxz_yyyzz_0, \
                                         tg_xxxxz_yyyzz_1, tg_xxxxz_yyzzz_0, tg_xxxxz_yyzzz_1, tg_xxxxz_yzzzz_0, \
                                         tg_xxxxz_yzzzz_1, tg_xxxxz_zzzzz_0, tg_xxxxz_zzzzz_1, tg_xxxyy_xxxxx_0, \
                                         tg_xxxyy_xxxxx_1, tg_xxxyy_xxxxy_0, tg_xxxyy_xxxxy_1, tg_xxxyy_xxxxz_0, \
                                         tg_xxxyy_xxxxz_1, tg_xxxyy_xxxyy_0, tg_xxxyy_xxxyy_1, tg_xxxyy_xxxyz_0, \
                                         tg_xxxyy_xxxyz_1, tg_xxxyy_xxxzz_0, tg_xxxyy_xxxzz_1, tg_xxxyy_xxyyy_0, \
                                         tg_xxxyy_xxyyy_1, tg_xxxyy_xxyyz_0, tg_xxxyy_xxyyz_1, tg_xxxyy_xxyzz_0, \
                                         tg_xxxyy_xxyzz_1, tg_xxxyy_xxzzz_0, tg_xxxyy_xxzzz_1, tg_xxxyy_xyyyy_0, \
                                         tg_xxxyy_xyyyy_1, tg_xxxyy_xyyyz_0, tg_xxxyy_xyyyz_1, tg_xxxyy_xyyzz_0, \
                                         tg_xxxyy_xyyzz_1, tg_xxxyy_xyzzz_0, tg_xxxyy_xyzzz_1, tg_xxxyy_xzzzz_0, \
                                         tg_xxxyy_xzzzz_1, tg_xxxyy_yyyyy_0, tg_xxxyy_yyyyy_1, tg_xxxyy_yyyyz_0, \
                                         tg_xxxyy_yyyyz_1, tg_xxxyy_yyyzz_0, tg_xxxyy_yyyzz_1, tg_xxxyy_yyzzz_0, \
                                         tg_xxxyy_yyzzz_1, tg_xxxyy_yzzzz_0, tg_xxxyy_yzzzz_1, tg_xxxyy_zzzzz_0, \
                                         tg_xxxyy_zzzzz_1, tg_xxxyz_xxxxx_0, tg_xxxyz_xxxxx_1, tg_xxxyz_xxxxy_0, \
                                         tg_xxxyz_xxxxy_1, tg_xxxyz_xxxxz_0, tg_xxxyz_xxxxz_1, tg_xxxyz_xxxyy_0, \
                                         tg_xxxyz_xxxyy_1, tg_xxxyz_xxxyz_0, tg_xxxyz_xxxyz_1, tg_xxxyz_xxxzz_0, \
                                         tg_xxxyz_xxxzz_1, tg_xxxyz_xxyyy_0, tg_xxxyz_xxyyy_1, tg_xxxyz_xxyyz_0, \
                                         tg_xxxyz_xxyyz_1, tg_xxxyz_xxyzz_0, tg_xxxyz_xxyzz_1, tg_xxxyz_xxzzz_0, \
                                         tg_xxxyz_xxzzz_1, tg_xxxyz_xyyyy_0, tg_xxxyz_xyyyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxxx_xxxxx_0[j] = pb_x * tg_xxxxxx_xxxxx_0[j] + fr * tg_xxxxxx_xxxxx_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxx_0[j] - tg_xxxxx_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxx_xxxx_1[j];

                    tg_xxxxxxx_xxxxy_0[j] = pb_x * tg_xxxxxx_xxxxy_0[j] + fr * tg_xxxxxx_xxxxy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxy_0[j] - tg_xxxxx_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxx_xxxy_1[j];

                    tg_xxxxxxx_xxxxz_0[j] = pb_x * tg_xxxxxx_xxxxz_0[j] + fr * tg_xxxxxx_xxxxz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxz_0[j] - tg_xxxxx_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxx_xxxz_1[j];

                    tg_xxxxxxx_xxxyy_0[j] = pb_x * tg_xxxxxx_xxxyy_0[j] + fr * tg_xxxxxx_xxxyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxyy_0[j] - tg_xxxxx_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxyy_1[j];

                    tg_xxxxxxx_xxxyz_0[j] = pb_x * tg_xxxxxx_xxxyz_0[j] + fr * tg_xxxxxx_xxxyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxyz_0[j] - tg_xxxxx_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxyz_1[j];

                    tg_xxxxxxx_xxxzz_0[j] = pb_x * tg_xxxxxx_xxxzz_0[j] + fr * tg_xxxxxx_xxxzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxzz_0[j] - tg_xxxxx_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxzz_1[j];

                    tg_xxxxxxx_xxyyy_0[j] = pb_x * tg_xxxxxx_xxyyy_0[j] + fr * tg_xxxxxx_xxyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyyy_0[j] - tg_xxxxx_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyyy_1[j];

                    tg_xxxxxxx_xxyyz_0[j] = pb_x * tg_xxxxxx_xxyyz_0[j] + fr * tg_xxxxxx_xxyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyyz_0[j] - tg_xxxxx_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyyz_1[j];

                    tg_xxxxxxx_xxyzz_0[j] = pb_x * tg_xxxxxx_xxyzz_0[j] + fr * tg_xxxxxx_xxyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyzz_0[j] - tg_xxxxx_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyzz_1[j];

                    tg_xxxxxxx_xxzzz_0[j] = pb_x * tg_xxxxxx_xxzzz_0[j] + fr * tg_xxxxxx_xxzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxzzz_0[j] - tg_xxxxx_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xzzz_1[j];

                    tg_xxxxxxx_xyyyy_0[j] = pb_x * tg_xxxxxx_xyyyy_0[j] + fr * tg_xxxxxx_xyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyyy_0[j] - tg_xxxxx_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyyy_1[j];

                    tg_xxxxxxx_xyyyz_0[j] = pb_x * tg_xxxxxx_xyyyz_0[j] + fr * tg_xxxxxx_xyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyyz_0[j] - tg_xxxxx_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyyz_1[j];

                    tg_xxxxxxx_xyyzz_0[j] = pb_x * tg_xxxxxx_xyyzz_0[j] + fr * tg_xxxxxx_xyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyzz_0[j] - tg_xxxxx_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyzz_1[j];

                    tg_xxxxxxx_xyzzz_0[j] = pb_x * tg_xxxxxx_xyzzz_0[j] + fr * tg_xxxxxx_xyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyzzz_0[j] - tg_xxxxx_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yzzz_1[j];

                    tg_xxxxxxx_xzzzz_0[j] = pb_x * tg_xxxxxx_xzzzz_0[j] + fr * tg_xxxxxx_xzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xzzzz_0[j] - tg_xxxxx_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_zzzz_1[j];

                    tg_xxxxxxx_yyyyy_0[j] = pb_x * tg_xxxxxx_yyyyy_0[j] + fr * tg_xxxxxx_yyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyyy_0[j] - tg_xxxxx_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxx_yyyyz_0[j] = pb_x * tg_xxxxxx_yyyyz_0[j] + fr * tg_xxxxxx_yyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyyz_0[j] - tg_xxxxx_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxx_yyyzz_0[j] = pb_x * tg_xxxxxx_yyyzz_0[j] + fr * tg_xxxxxx_yyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyzz_0[j] - tg_xxxxx_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxx_yyzzz_0[j] = pb_x * tg_xxxxxx_yyzzz_0[j] + fr * tg_xxxxxx_yyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyzzz_0[j] - tg_xxxxx_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxx_yzzzz_0[j] = pb_x * tg_xxxxxx_yzzzz_0[j] + fr * tg_xxxxxx_yzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yzzzz_0[j] - tg_xxxxx_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxx_zzzzz_0[j] = pb_x * tg_xxxxxx_zzzzz_0[j] + fr * tg_xxxxxx_zzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_zzzzz_0[j] - tg_xxxxx_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_xxxxx_0[j] = pb_x * tg_xxxxxy_xxxxx_0[j] + fr * tg_xxxxxy_xxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxx_0[j] - tg_xxxxy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxy_xxxx_1[j];

                    tg_xxxxxxy_xxxxy_0[j] = pb_x * tg_xxxxxy_xxxxy_0[j] + fr * tg_xxxxxy_xxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxy_0[j] - tg_xxxxy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxy_xxxy_1[j];

                    tg_xxxxxxy_xxxxz_0[j] = pb_x * tg_xxxxxy_xxxxz_0[j] + fr * tg_xxxxxy_xxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxz_0[j] - tg_xxxxy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxy_xxxz_1[j];

                    tg_xxxxxxy_xxxyy_0[j] = pb_x * tg_xxxxxy_xxxyy_0[j] + fr * tg_xxxxxy_xxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxyy_0[j] - tg_xxxxy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxyy_1[j];

                    tg_xxxxxxy_xxxyz_0[j] = pb_x * tg_xxxxxy_xxxyz_0[j] + fr * tg_xxxxxy_xxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxyz_0[j] - tg_xxxxy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxyz_1[j];

                    tg_xxxxxxy_xxxzz_0[j] = pb_x * tg_xxxxxy_xxxzz_0[j] + fr * tg_xxxxxy_xxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxzz_0[j] - tg_xxxxy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxzz_1[j];

                    tg_xxxxxxy_xxyyy_0[j] = pb_x * tg_xxxxxy_xxyyy_0[j] + fr * tg_xxxxxy_xxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyyy_0[j] - tg_xxxxy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyyy_1[j];

                    tg_xxxxxxy_xxyyz_0[j] = pb_x * tg_xxxxxy_xxyyz_0[j] + fr * tg_xxxxxy_xxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyyz_0[j] - tg_xxxxy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyyz_1[j];

                    tg_xxxxxxy_xxyzz_0[j] = pb_x * tg_xxxxxy_xxyzz_0[j] + fr * tg_xxxxxy_xxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyzz_0[j] - tg_xxxxy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyzz_1[j];

                    tg_xxxxxxy_xxzzz_0[j] = pb_x * tg_xxxxxy_xxzzz_0[j] + fr * tg_xxxxxy_xxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxzzz_0[j] - tg_xxxxy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xzzz_1[j];

                    tg_xxxxxxy_xyyyy_0[j] = pb_x * tg_xxxxxy_xyyyy_0[j] + fr * tg_xxxxxy_xyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyyy_0[j] - tg_xxxxy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyyy_1[j];

                    tg_xxxxxxy_xyyyz_0[j] = pb_x * tg_xxxxxy_xyyyz_0[j] + fr * tg_xxxxxy_xyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyyz_0[j] - tg_xxxxy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyyz_1[j];

                    tg_xxxxxxy_xyyzz_0[j] = pb_x * tg_xxxxxy_xyyzz_0[j] + fr * tg_xxxxxy_xyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyzz_0[j] - tg_xxxxy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyzz_1[j];

                    tg_xxxxxxy_xyzzz_0[j] = pb_x * tg_xxxxxy_xyzzz_0[j] + fr * tg_xxxxxy_xyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyzzz_0[j] - tg_xxxxy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yzzz_1[j];

                    tg_xxxxxxy_xzzzz_0[j] = pb_x * tg_xxxxxy_xzzzz_0[j] + fr * tg_xxxxxy_xzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xzzzz_0[j] - tg_xxxxy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_zzzz_1[j];

                    tg_xxxxxxy_yyyyy_0[j] = pb_x * tg_xxxxxy_yyyyy_0[j] + fr * tg_xxxxxy_yyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyyy_0[j] - tg_xxxxy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxy_yyyyz_0[j] = pb_x * tg_xxxxxy_yyyyz_0[j] + fr * tg_xxxxxy_yyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyyz_0[j] - tg_xxxxy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxy_yyyzz_0[j] = pb_x * tg_xxxxxy_yyyzz_0[j] + fr * tg_xxxxxy_yyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyzz_0[j] - tg_xxxxy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxy_yyzzz_0[j] = pb_x * tg_xxxxxy_yyzzz_0[j] + fr * tg_xxxxxy_yyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyzzz_0[j] - tg_xxxxy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_yzzzz_0[j] = pb_x * tg_xxxxxy_yzzzz_0[j] + fr * tg_xxxxxy_yzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yzzzz_0[j] - tg_xxxxy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_zzzzz_0[j] = pb_x * tg_xxxxxy_zzzzz_0[j] + fr * tg_xxxxxy_zzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_zzzzz_0[j] - tg_xxxxy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_xxxxx_0[j] = pb_x * tg_xxxxxz_xxxxx_0[j] + fr * tg_xxxxxz_xxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxx_0[j] - tg_xxxxz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxz_xxxx_1[j];

                    tg_xxxxxxz_xxxxy_0[j] = pb_x * tg_xxxxxz_xxxxy_0[j] + fr * tg_xxxxxz_xxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxy_0[j] - tg_xxxxz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxz_xxxy_1[j];

                    tg_xxxxxxz_xxxxz_0[j] = pb_x * tg_xxxxxz_xxxxz_0[j] + fr * tg_xxxxxz_xxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxz_0[j] - tg_xxxxz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxz_xxxz_1[j];

                    tg_xxxxxxz_xxxyy_0[j] = pb_x * tg_xxxxxz_xxxyy_0[j] + fr * tg_xxxxxz_xxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxyy_0[j] - tg_xxxxz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxyy_1[j];

                    tg_xxxxxxz_xxxyz_0[j] = pb_x * tg_xxxxxz_xxxyz_0[j] + fr * tg_xxxxxz_xxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxyz_0[j] - tg_xxxxz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxyz_1[j];

                    tg_xxxxxxz_xxxzz_0[j] = pb_x * tg_xxxxxz_xxxzz_0[j] + fr * tg_xxxxxz_xxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxzz_0[j] - tg_xxxxz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxzz_1[j];

                    tg_xxxxxxz_xxyyy_0[j] = pb_x * tg_xxxxxz_xxyyy_0[j] + fr * tg_xxxxxz_xxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyyy_0[j] - tg_xxxxz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyyy_1[j];

                    tg_xxxxxxz_xxyyz_0[j] = pb_x * tg_xxxxxz_xxyyz_0[j] + fr * tg_xxxxxz_xxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyyz_0[j] - tg_xxxxz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyyz_1[j];

                    tg_xxxxxxz_xxyzz_0[j] = pb_x * tg_xxxxxz_xxyzz_0[j] + fr * tg_xxxxxz_xxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyzz_0[j] - tg_xxxxz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyzz_1[j];

                    tg_xxxxxxz_xxzzz_0[j] = pb_x * tg_xxxxxz_xxzzz_0[j] + fr * tg_xxxxxz_xxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxzzz_0[j] - tg_xxxxz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xzzz_1[j];

                    tg_xxxxxxz_xyyyy_0[j] = pb_x * tg_xxxxxz_xyyyy_0[j] + fr * tg_xxxxxz_xyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyyy_0[j] - tg_xxxxz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyyy_1[j];

                    tg_xxxxxxz_xyyyz_0[j] = pb_x * tg_xxxxxz_xyyyz_0[j] + fr * tg_xxxxxz_xyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyyz_0[j] - tg_xxxxz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyyz_1[j];

                    tg_xxxxxxz_xyyzz_0[j] = pb_x * tg_xxxxxz_xyyzz_0[j] + fr * tg_xxxxxz_xyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyzz_0[j] - tg_xxxxz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyzz_1[j];

                    tg_xxxxxxz_xyzzz_0[j] = pb_x * tg_xxxxxz_xyzzz_0[j] + fr * tg_xxxxxz_xyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyzzz_0[j] - tg_xxxxz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yzzz_1[j];

                    tg_xxxxxxz_xzzzz_0[j] = pb_x * tg_xxxxxz_xzzzz_0[j] + fr * tg_xxxxxz_xzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xzzzz_0[j] - tg_xxxxz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_zzzz_1[j];

                    tg_xxxxxxz_yyyyy_0[j] = pb_x * tg_xxxxxz_yyyyy_0[j] + fr * tg_xxxxxz_yyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyyy_0[j] - tg_xxxxz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxz_yyyyz_0[j] = pb_x * tg_xxxxxz_yyyyz_0[j] + fr * tg_xxxxxz_yyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyyz_0[j] - tg_xxxxz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxz_yyyzz_0[j] = pb_x * tg_xxxxxz_yyyzz_0[j] + fr * tg_xxxxxz_yyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyzz_0[j] - tg_xxxxz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxz_yyzzz_0[j] = pb_x * tg_xxxxxz_yyzzz_0[j] + fr * tg_xxxxxz_yyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyzzz_0[j] - tg_xxxxz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_yzzzz_0[j] = pb_x * tg_xxxxxz_yzzzz_0[j] + fr * tg_xxxxxz_yzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yzzzz_0[j] - tg_xxxxz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_zzzzz_0[j] = pb_x * tg_xxxxxz_zzzzz_0[j] + fr * tg_xxxxxz_zzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_zzzzz_0[j] - tg_xxxxz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_xxxxx_0[j] = pb_x * tg_xxxxyy_xxxxx_0[j] + fr * tg_xxxxyy_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxx_0[j] - tg_xxxyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyy_xxxx_1[j];

                    tg_xxxxxyy_xxxxy_0[j] = pb_x * tg_xxxxyy_xxxxy_0[j] + fr * tg_xxxxyy_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxy_0[j] - tg_xxxyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyy_xxxy_1[j];

                    tg_xxxxxyy_xxxxz_0[j] = pb_x * tg_xxxxyy_xxxxz_0[j] + fr * tg_xxxxyy_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxz_0[j] - tg_xxxyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyy_xxxz_1[j];

                    tg_xxxxxyy_xxxyy_0[j] = pb_x * tg_xxxxyy_xxxyy_0[j] + fr * tg_xxxxyy_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxyy_0[j] - tg_xxxyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxyy_1[j];

                    tg_xxxxxyy_xxxyz_0[j] = pb_x * tg_xxxxyy_xxxyz_0[j] + fr * tg_xxxxyy_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxyz_0[j] - tg_xxxyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxyz_1[j];

                    tg_xxxxxyy_xxxzz_0[j] = pb_x * tg_xxxxyy_xxxzz_0[j] + fr * tg_xxxxyy_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxzz_0[j] - tg_xxxyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxzz_1[j];

                    tg_xxxxxyy_xxyyy_0[j] = pb_x * tg_xxxxyy_xxyyy_0[j] + fr * tg_xxxxyy_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyyy_0[j] - tg_xxxyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyyy_1[j];

                    tg_xxxxxyy_xxyyz_0[j] = pb_x * tg_xxxxyy_xxyyz_0[j] + fr * tg_xxxxyy_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyyz_0[j] - tg_xxxyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyyz_1[j];

                    tg_xxxxxyy_xxyzz_0[j] = pb_x * tg_xxxxyy_xxyzz_0[j] + fr * tg_xxxxyy_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyzz_0[j] - tg_xxxyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyzz_1[j];

                    tg_xxxxxyy_xxzzz_0[j] = pb_x * tg_xxxxyy_xxzzz_0[j] + fr * tg_xxxxyy_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxzzz_0[j] - tg_xxxyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xzzz_1[j];

                    tg_xxxxxyy_xyyyy_0[j] = pb_x * tg_xxxxyy_xyyyy_0[j] + fr * tg_xxxxyy_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyyy_0[j] - tg_xxxyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyyy_1[j];

                    tg_xxxxxyy_xyyyz_0[j] = pb_x * tg_xxxxyy_xyyyz_0[j] + fr * tg_xxxxyy_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyyz_0[j] - tg_xxxyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyyz_1[j];

                    tg_xxxxxyy_xyyzz_0[j] = pb_x * tg_xxxxyy_xyyzz_0[j] + fr * tg_xxxxyy_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyzz_0[j] - tg_xxxyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyzz_1[j];

                    tg_xxxxxyy_xyzzz_0[j] = pb_x * tg_xxxxyy_xyzzz_0[j] + fr * tg_xxxxyy_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyzzz_0[j] - tg_xxxyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yzzz_1[j];

                    tg_xxxxxyy_xzzzz_0[j] = pb_x * tg_xxxxyy_xzzzz_0[j] + fr * tg_xxxxyy_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xzzzz_0[j] - tg_xxxyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_zzzz_1[j];

                    tg_xxxxxyy_yyyyy_0[j] = pb_x * tg_xxxxyy_yyyyy_0[j] + fr * tg_xxxxyy_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyyy_0[j] - tg_xxxyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxyy_yyyyz_0[j] = pb_x * tg_xxxxyy_yyyyz_0[j] + fr * tg_xxxxyy_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyyz_0[j] - tg_xxxyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxyy_yyyzz_0[j] = pb_x * tg_xxxxyy_yyyzz_0[j] + fr * tg_xxxxyy_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyzz_0[j] - tg_xxxyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxyy_yyzzz_0[j] = pb_x * tg_xxxxyy_yyzzz_0[j] + fr * tg_xxxxyy_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyzzz_0[j] - tg_xxxyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_yzzzz_0[j] = pb_x * tg_xxxxyy_yzzzz_0[j] + fr * tg_xxxxyy_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yzzzz_0[j] - tg_xxxyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_zzzzz_0[j] = pb_x * tg_xxxxyy_zzzzz_0[j] + fr * tg_xxxxyy_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_zzzzz_0[j] - tg_xxxyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_xxxxx_0[j] = pb_x * tg_xxxxyz_xxxxx_0[j] + fr * tg_xxxxyz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxx_0[j] - tg_xxxyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyz_xxxx_1[j];

                    tg_xxxxxyz_xxxxy_0[j] = pb_x * tg_xxxxyz_xxxxy_0[j] + fr * tg_xxxxyz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxy_0[j] - tg_xxxyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyz_xxxy_1[j];

                    tg_xxxxxyz_xxxxz_0[j] = pb_x * tg_xxxxyz_xxxxz_0[j] + fr * tg_xxxxyz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxz_0[j] - tg_xxxyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyz_xxxz_1[j];

                    tg_xxxxxyz_xxxyy_0[j] = pb_x * tg_xxxxyz_xxxyy_0[j] + fr * tg_xxxxyz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxyy_0[j] - tg_xxxyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxyy_1[j];

                    tg_xxxxxyz_xxxyz_0[j] = pb_x * tg_xxxxyz_xxxyz_0[j] + fr * tg_xxxxyz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxyz_0[j] - tg_xxxyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxyz_1[j];

                    tg_xxxxxyz_xxxzz_0[j] = pb_x * tg_xxxxyz_xxxzz_0[j] + fr * tg_xxxxyz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxzz_0[j] - tg_xxxyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxzz_1[j];

                    tg_xxxxxyz_xxyyy_0[j] = pb_x * tg_xxxxyz_xxyyy_0[j] + fr * tg_xxxxyz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyyy_0[j] - tg_xxxyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyyy_1[j];

                    tg_xxxxxyz_xxyyz_0[j] = pb_x * tg_xxxxyz_xxyyz_0[j] + fr * tg_xxxxyz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyyz_0[j] - tg_xxxyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyyz_1[j];

                    tg_xxxxxyz_xxyzz_0[j] = pb_x * tg_xxxxyz_xxyzz_0[j] + fr * tg_xxxxyz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyzz_0[j] - tg_xxxyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyzz_1[j];

                    tg_xxxxxyz_xxzzz_0[j] = pb_x * tg_xxxxyz_xxzzz_0[j] + fr * tg_xxxxyz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxzzz_0[j] - tg_xxxyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xzzz_1[j];

                    tg_xxxxxyz_xyyyy_0[j] = pb_x * tg_xxxxyz_xyyyy_0[j] + fr * tg_xxxxyz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyyy_0[j] - tg_xxxyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSH_95_190(      CMemBlock2D<double>* primBuffer,
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
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_xxxxyz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 95); 

                auto tg_xxxxyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 96); 

                auto tg_xxxxyz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 97); 

                auto tg_xxxxyz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 98); 

                auto tg_xxxxyz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 99); 

                auto tg_xxxxyz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 100); 

                auto tg_xxxxyz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 101); 

                auto tg_xxxxyz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 102); 

                auto tg_xxxxyz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 103); 

                auto tg_xxxxyz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 104); 

                auto tg_xxxxzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 105); 

                auto tg_xxxxzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 106); 

                auto tg_xxxxzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 107); 

                auto tg_xxxxzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 108); 

                auto tg_xxxxzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 109); 

                auto tg_xxxxzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 110); 

                auto tg_xxxxzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 111); 

                auto tg_xxxxzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 112); 

                auto tg_xxxxzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 113); 

                auto tg_xxxxzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 114); 

                auto tg_xxxxzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 115); 

                auto tg_xxxxzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 116); 

                auto tg_xxxxzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 117); 

                auto tg_xxxxzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 118); 

                auto tg_xxxxzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 119); 

                auto tg_xxxxzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 120); 

                auto tg_xxxxzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 121); 

                auto tg_xxxxzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 122); 

                auto tg_xxxxzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 123); 

                auto tg_xxxxzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 124); 

                auto tg_xxxxzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 125); 

                auto tg_xxxyyy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 126); 

                auto tg_xxxyyy_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 127); 

                auto tg_xxxyyy_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 128); 

                auto tg_xxxyyy_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 129); 

                auto tg_xxxyyy_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 130); 

                auto tg_xxxyyy_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 131); 

                auto tg_xxxyyy_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 132); 

                auto tg_xxxyyy_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 133); 

                auto tg_xxxyyy_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 134); 

                auto tg_xxxyyy_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 135); 

                auto tg_xxxyyy_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 136); 

                auto tg_xxxyyy_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 137); 

                auto tg_xxxyyy_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 138); 

                auto tg_xxxyyy_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 139); 

                auto tg_xxxyyy_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 140); 

                auto tg_xxxyyy_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 141); 

                auto tg_xxxyyy_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 142); 

                auto tg_xxxyyy_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 143); 

                auto tg_xxxyyy_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 144); 

                auto tg_xxxyyy_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 145); 

                auto tg_xxxyyy_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 146); 

                auto tg_xxxyyz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 147); 

                auto tg_xxxyyz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 148); 

                auto tg_xxxyyz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 149); 

                auto tg_xxxyyz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 150); 

                auto tg_xxxyyz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 151); 

                auto tg_xxxyyz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 152); 

                auto tg_xxxyyz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 153); 

                auto tg_xxxyyz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 154); 

                auto tg_xxxyyz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 155); 

                auto tg_xxxyyz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 156); 

                auto tg_xxxyyz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 157); 

                auto tg_xxxyyz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 158); 

                auto tg_xxxyyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 159); 

                auto tg_xxxyyz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 160); 

                auto tg_xxxyyz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 161); 

                auto tg_xxxyyz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 162); 

                auto tg_xxxyyz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 163); 

                auto tg_xxxyyz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 164); 

                auto tg_xxxyyz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 165); 

                auto tg_xxxyyz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 166); 

                auto tg_xxxyyz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 167); 

                auto tg_xxxyzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 168); 

                auto tg_xxxyzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 169); 

                auto tg_xxxyzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 170); 

                auto tg_xxxyzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 171); 

                auto tg_xxxyzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 172); 

                auto tg_xxxyzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 173); 

                auto tg_xxxyzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 174); 

                auto tg_xxxyzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 175); 

                auto tg_xxxyzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 176); 

                auto tg_xxxyzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 177); 

                auto tg_xxxyzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 178); 

                auto tg_xxxyzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 179); 

                auto tg_xxxyzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 180); 

                auto tg_xxxyzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 181); 

                auto tg_xxxyzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 182); 

                auto tg_xxxyzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 183); 

                auto tg_xxxyzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 184); 

                auto tg_xxxyzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 185); 

                auto tg_xxxyzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 186); 

                auto tg_xxxyzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 187); 

                auto tg_xxxyzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 188); 

                auto tg_xxxzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 189); 

                auto tg_xxxxyz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 95); 

                auto tg_xxxxyz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 96); 

                auto tg_xxxxyz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 97); 

                auto tg_xxxxyz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 98); 

                auto tg_xxxxyz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 99); 

                auto tg_xxxxyz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 100); 

                auto tg_xxxxyz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 101); 

                auto tg_xxxxyz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 102); 

                auto tg_xxxxyz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 103); 

                auto tg_xxxxyz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 104); 

                auto tg_xxxxzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 105); 

                auto tg_xxxxzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 106); 

                auto tg_xxxxzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 107); 

                auto tg_xxxxzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 108); 

                auto tg_xxxxzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 109); 

                auto tg_xxxxzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 110); 

                auto tg_xxxxzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 111); 

                auto tg_xxxxzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 112); 

                auto tg_xxxxzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 113); 

                auto tg_xxxxzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 114); 

                auto tg_xxxxzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 115); 

                auto tg_xxxxzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 116); 

                auto tg_xxxxzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 117); 

                auto tg_xxxxzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 118); 

                auto tg_xxxxzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 119); 

                auto tg_xxxxzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 120); 

                auto tg_xxxxzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 121); 

                auto tg_xxxxzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 122); 

                auto tg_xxxxzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 123); 

                auto tg_xxxxzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 124); 

                auto tg_xxxxzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 125); 

                auto tg_xxxyyy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 126); 

                auto tg_xxxyyy_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 127); 

                auto tg_xxxyyy_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 128); 

                auto tg_xxxyyy_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 129); 

                auto tg_xxxyyy_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 130); 

                auto tg_xxxyyy_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 131); 

                auto tg_xxxyyy_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 132); 

                auto tg_xxxyyy_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 133); 

                auto tg_xxxyyy_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 134); 

                auto tg_xxxyyy_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 135); 

                auto tg_xxxyyy_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 136); 

                auto tg_xxxyyy_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 137); 

                auto tg_xxxyyy_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 138); 

                auto tg_xxxyyy_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 139); 

                auto tg_xxxyyy_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 140); 

                auto tg_xxxyyy_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 141); 

                auto tg_xxxyyy_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 142); 

                auto tg_xxxyyy_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 143); 

                auto tg_xxxyyy_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 144); 

                auto tg_xxxyyy_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 145); 

                auto tg_xxxyyy_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 146); 

                auto tg_xxxyyz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 147); 

                auto tg_xxxyyz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 148); 

                auto tg_xxxyyz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 149); 

                auto tg_xxxyyz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 150); 

                auto tg_xxxyyz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 151); 

                auto tg_xxxyyz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 152); 

                auto tg_xxxyyz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 153); 

                auto tg_xxxyyz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 154); 

                auto tg_xxxyyz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 155); 

                auto tg_xxxyyz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 156); 

                auto tg_xxxyyz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 157); 

                auto tg_xxxyyz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 158); 

                auto tg_xxxyyz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 159); 

                auto tg_xxxyyz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 160); 

                auto tg_xxxyyz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 161); 

                auto tg_xxxyyz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 162); 

                auto tg_xxxyyz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 163); 

                auto tg_xxxyyz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 164); 

                auto tg_xxxyyz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 165); 

                auto tg_xxxyyz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 166); 

                auto tg_xxxyyz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 167); 

                auto tg_xxxyzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 168); 

                auto tg_xxxyzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 169); 

                auto tg_xxxyzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 170); 

                auto tg_xxxyzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 171); 

                auto tg_xxxyzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 172); 

                auto tg_xxxyzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 173); 

                auto tg_xxxyzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 174); 

                auto tg_xxxyzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 175); 

                auto tg_xxxyzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 176); 

                auto tg_xxxyzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 177); 

                auto tg_xxxyzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 178); 

                auto tg_xxxyzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 179); 

                auto tg_xxxyzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 180); 

                auto tg_xxxyzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 181); 

                auto tg_xxxyzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 182); 

                auto tg_xxxyzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 183); 

                auto tg_xxxyzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 184); 

                auto tg_xxxyzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 185); 

                auto tg_xxxyzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 186); 

                auto tg_xxxyzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 187); 

                auto tg_xxxyzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 188); 

                auto tg_xxxzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 189); 

                auto tg_xxxyz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 95); 

                auto tg_xxxyz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 96); 

                auto tg_xxxyz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 97); 

                auto tg_xxxyz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 98); 

                auto tg_xxxyz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 99); 

                auto tg_xxxyz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 100); 

                auto tg_xxxyz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 101); 

                auto tg_xxxyz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 102); 

                auto tg_xxxyz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 103); 

                auto tg_xxxyz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 104); 

                auto tg_xxxzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 105); 

                auto tg_xxxzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 106); 

                auto tg_xxxzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 107); 

                auto tg_xxxzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 108); 

                auto tg_xxxzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 109); 

                auto tg_xxxzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 110); 

                auto tg_xxxzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 111); 

                auto tg_xxxzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 112); 

                auto tg_xxxzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 113); 

                auto tg_xxxzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 114); 

                auto tg_xxxzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 115); 

                auto tg_xxxzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 116); 

                auto tg_xxxzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 117); 

                auto tg_xxxzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 118); 

                auto tg_xxxzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 119); 

                auto tg_xxxzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 120); 

                auto tg_xxxzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 121); 

                auto tg_xxxzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 122); 

                auto tg_xxxzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 123); 

                auto tg_xxxzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 124); 

                auto tg_xxxzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 125); 

                auto tg_xxyyy_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 126); 

                auto tg_xxyyy_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 127); 

                auto tg_xxyyy_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 128); 

                auto tg_xxyyy_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 129); 

                auto tg_xxyyy_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 130); 

                auto tg_xxyyy_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 131); 

                auto tg_xxyyy_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 132); 

                auto tg_xxyyy_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 133); 

                auto tg_xxyyy_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 134); 

                auto tg_xxyyy_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 135); 

                auto tg_xxyyy_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 136); 

                auto tg_xxyyy_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 137); 

                auto tg_xxyyy_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 138); 

                auto tg_xxyyy_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 139); 

                auto tg_xxyyy_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 140); 

                auto tg_xxyyy_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 141); 

                auto tg_xxyyy_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 142); 

                auto tg_xxyyy_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 143); 

                auto tg_xxyyy_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 144); 

                auto tg_xxyyy_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 145); 

                auto tg_xxyyy_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 146); 

                auto tg_xxyyz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 147); 

                auto tg_xxyyz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 148); 

                auto tg_xxyyz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 149); 

                auto tg_xxyyz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 150); 

                auto tg_xxyyz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 151); 

                auto tg_xxyyz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 152); 

                auto tg_xxyyz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 153); 

                auto tg_xxyyz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 154); 

                auto tg_xxyyz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 155); 

                auto tg_xxyyz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 156); 

                auto tg_xxyyz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 157); 

                auto tg_xxyyz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 158); 

                auto tg_xxyyz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 159); 

                auto tg_xxyyz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 160); 

                auto tg_xxyyz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 161); 

                auto tg_xxyyz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 162); 

                auto tg_xxyyz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 163); 

                auto tg_xxyyz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 164); 

                auto tg_xxyyz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 165); 

                auto tg_xxyyz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 166); 

                auto tg_xxyyz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 167); 

                auto tg_xxyzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 168); 

                auto tg_xxyzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 169); 

                auto tg_xxyzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 170); 

                auto tg_xxyzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 171); 

                auto tg_xxyzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 172); 

                auto tg_xxyzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 173); 

                auto tg_xxyzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 174); 

                auto tg_xxyzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 175); 

                auto tg_xxyzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 176); 

                auto tg_xxyzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 177); 

                auto tg_xxyzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 178); 

                auto tg_xxyzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 179); 

                auto tg_xxyzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 180); 

                auto tg_xxyzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 181); 

                auto tg_xxyzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 182); 

                auto tg_xxyzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 183); 

                auto tg_xxyzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 184); 

                auto tg_xxyzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 185); 

                auto tg_xxyzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 186); 

                auto tg_xxyzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 187); 

                auto tg_xxyzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 188); 

                auto tg_xxzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 189); 

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

                auto tg_xxxxyz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 71); 

                auto tg_xxxxyz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 72); 

                auto tg_xxxxyz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 73); 

                auto tg_xxxxyz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 74); 

                auto tg_xxxxzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 75); 

                auto tg_xxxxzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 76); 

                auto tg_xxxxzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 77); 

                auto tg_xxxxzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 78); 

                auto tg_xxxxzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 79); 

                auto tg_xxxxzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 80); 

                auto tg_xxxxzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 81); 

                auto tg_xxxxzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 82); 

                auto tg_xxxxzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 83); 

                auto tg_xxxxzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 84); 

                auto tg_xxxxzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 85); 

                auto tg_xxxxzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 86); 

                auto tg_xxxxzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 87); 

                auto tg_xxxxzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 88); 

                auto tg_xxxxzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 89); 

                auto tg_xxxyyy_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 90); 

                auto tg_xxxyyy_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 91); 

                auto tg_xxxyyy_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 92); 

                auto tg_xxxyyy_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 93); 

                auto tg_xxxyyy_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 94); 

                auto tg_xxxyyy_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 95); 

                auto tg_xxxyyy_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 96); 

                auto tg_xxxyyy_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 97); 

                auto tg_xxxyyy_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 98); 

                auto tg_xxxyyy_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 99); 

                auto tg_xxxyyy_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 100); 

                auto tg_xxxyyy_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 101); 

                auto tg_xxxyyy_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 102); 

                auto tg_xxxyyy_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 103); 

                auto tg_xxxyyy_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 104); 

                auto tg_xxxyyz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 105); 

                auto tg_xxxyyz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 106); 

                auto tg_xxxyyz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 107); 

                auto tg_xxxyyz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 108); 

                auto tg_xxxyyz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 109); 

                auto tg_xxxyyz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 110); 

                auto tg_xxxyyz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 111); 

                auto tg_xxxyyz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 112); 

                auto tg_xxxyyz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 113); 

                auto tg_xxxyyz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 114); 

                auto tg_xxxyyz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 115); 

                auto tg_xxxyyz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 116); 

                auto tg_xxxyyz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 117); 

                auto tg_xxxyyz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 118); 

                auto tg_xxxyyz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 119); 

                auto tg_xxxyzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 120); 

                auto tg_xxxyzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 121); 

                auto tg_xxxyzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 122); 

                auto tg_xxxyzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 123); 

                auto tg_xxxyzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 124); 

                auto tg_xxxyzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 125); 

                auto tg_xxxyzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 126); 

                auto tg_xxxyzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 127); 

                auto tg_xxxyzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 128); 

                auto tg_xxxyzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 129); 

                auto tg_xxxyzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 130); 

                auto tg_xxxyzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 131); 

                auto tg_xxxyzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 132); 

                auto tg_xxxyzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 133); 

                auto tg_xxxyzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 134); 

                auto tg_xxxzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 135); 

                // set up pointers to integrals

                auto tg_xxxxxyz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 95); 

                auto tg_xxxxxyz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 96); 

                auto tg_xxxxxyz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 97); 

                auto tg_xxxxxyz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 98); 

                auto tg_xxxxxyz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 99); 

                auto tg_xxxxxyz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 100); 

                auto tg_xxxxxyz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 101); 

                auto tg_xxxxxyz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 102); 

                auto tg_xxxxxyz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 103); 

                auto tg_xxxxxyz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 104); 

                auto tg_xxxxxzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 105); 

                auto tg_xxxxxzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 106); 

                auto tg_xxxxxzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 107); 

                auto tg_xxxxxzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 108); 

                auto tg_xxxxxzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 109); 

                auto tg_xxxxxzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 110); 

                auto tg_xxxxxzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 111); 

                auto tg_xxxxxzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 112); 

                auto tg_xxxxxzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 113); 

                auto tg_xxxxxzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 114); 

                auto tg_xxxxxzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 115); 

                auto tg_xxxxxzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 116); 

                auto tg_xxxxxzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 117); 

                auto tg_xxxxxzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 118); 

                auto tg_xxxxxzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 119); 

                auto tg_xxxxxzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 120); 

                auto tg_xxxxxzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 121); 

                auto tg_xxxxxzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 122); 

                auto tg_xxxxxzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 123); 

                auto tg_xxxxxzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 124); 

                auto tg_xxxxxzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 125); 

                auto tg_xxxxyyy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 126); 

                auto tg_xxxxyyy_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 127); 

                auto tg_xxxxyyy_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 128); 

                auto tg_xxxxyyy_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 129); 

                auto tg_xxxxyyy_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 130); 

                auto tg_xxxxyyy_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 131); 

                auto tg_xxxxyyy_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 132); 

                auto tg_xxxxyyy_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 133); 

                auto tg_xxxxyyy_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 134); 

                auto tg_xxxxyyy_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 135); 

                auto tg_xxxxyyy_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 136); 

                auto tg_xxxxyyy_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 137); 

                auto tg_xxxxyyy_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 138); 

                auto tg_xxxxyyy_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 139); 

                auto tg_xxxxyyy_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 140); 

                auto tg_xxxxyyy_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 141); 

                auto tg_xxxxyyy_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 142); 

                auto tg_xxxxyyy_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 143); 

                auto tg_xxxxyyy_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 144); 

                auto tg_xxxxyyy_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 145); 

                auto tg_xxxxyyy_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 146); 

                auto tg_xxxxyyz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 147); 

                auto tg_xxxxyyz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 148); 

                auto tg_xxxxyyz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 149); 

                auto tg_xxxxyyz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 150); 

                auto tg_xxxxyyz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 151); 

                auto tg_xxxxyyz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 152); 

                auto tg_xxxxyyz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 153); 

                auto tg_xxxxyyz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 154); 

                auto tg_xxxxyyz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 155); 

                auto tg_xxxxyyz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 156); 

                auto tg_xxxxyyz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 157); 

                auto tg_xxxxyyz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 158); 

                auto tg_xxxxyyz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 159); 

                auto tg_xxxxyyz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 160); 

                auto tg_xxxxyyz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 161); 

                auto tg_xxxxyyz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 162); 

                auto tg_xxxxyyz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 163); 

                auto tg_xxxxyyz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 164); 

                auto tg_xxxxyyz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 165); 

                auto tg_xxxxyyz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 166); 

                auto tg_xxxxyyz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 167); 

                auto tg_xxxxyzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 168); 

                auto tg_xxxxyzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 169); 

                auto tg_xxxxyzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 170); 

                auto tg_xxxxyzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 171); 

                auto tg_xxxxyzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 172); 

                auto tg_xxxxyzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 173); 

                auto tg_xxxxyzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 174); 

                auto tg_xxxxyzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 175); 

                auto tg_xxxxyzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 176); 

                auto tg_xxxxyzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 177); 

                auto tg_xxxxyzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 178); 

                auto tg_xxxxyzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 179); 

                auto tg_xxxxyzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 180); 

                auto tg_xxxxyzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 181); 

                auto tg_xxxxyzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 182); 

                auto tg_xxxxyzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 183); 

                auto tg_xxxxyzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 184); 

                auto tg_xxxxyzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 185); 

                auto tg_xxxxyzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 186); 

                auto tg_xxxxyzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 187); 

                auto tg_xxxxyzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 188); 

                auto tg_xxxxzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 189); 

                // Batch of Integrals (95,190)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxyz_xyyyz_0, tg_xxxxxyz_xyyzz_0, tg_xxxxxyz_xyzzz_0, \
                                         tg_xxxxxyz_xzzzz_0, tg_xxxxxyz_yyyyy_0, tg_xxxxxyz_yyyyz_0, tg_xxxxxyz_yyyzz_0, \
                                         tg_xxxxxyz_yyzzz_0, tg_xxxxxyz_yzzzz_0, tg_xxxxxyz_zzzzz_0, tg_xxxxxzz_xxxxx_0, \
                                         tg_xxxxxzz_xxxxy_0, tg_xxxxxzz_xxxxz_0, tg_xxxxxzz_xxxyy_0, tg_xxxxxzz_xxxyz_0, \
                                         tg_xxxxxzz_xxxzz_0, tg_xxxxxzz_xxyyy_0, tg_xxxxxzz_xxyyz_0, tg_xxxxxzz_xxyzz_0, \
                                         tg_xxxxxzz_xxzzz_0, tg_xxxxxzz_xyyyy_0, tg_xxxxxzz_xyyyz_0, tg_xxxxxzz_xyyzz_0, \
                                         tg_xxxxxzz_xyzzz_0, tg_xxxxxzz_xzzzz_0, tg_xxxxxzz_yyyyy_0, tg_xxxxxzz_yyyyz_0, \
                                         tg_xxxxxzz_yyyzz_0, tg_xxxxxzz_yyzzz_0, tg_xxxxxzz_yzzzz_0, tg_xxxxxzz_zzzzz_0, \
                                         tg_xxxxyyy_xxxxx_0, tg_xxxxyyy_xxxxy_0, tg_xxxxyyy_xxxxz_0, tg_xxxxyyy_xxxyy_0, \
                                         tg_xxxxyyy_xxxyz_0, tg_xxxxyyy_xxxzz_0, tg_xxxxyyy_xxyyy_0, tg_xxxxyyy_xxyyz_0, \
                                         tg_xxxxyyy_xxyzz_0, tg_xxxxyyy_xxzzz_0, tg_xxxxyyy_xyyyy_0, tg_xxxxyyy_xyyyz_0, \
                                         tg_xxxxyyy_xyyzz_0, tg_xxxxyyy_xyzzz_0, tg_xxxxyyy_xzzzz_0, tg_xxxxyyy_yyyyy_0, \
                                         tg_xxxxyyy_yyyyz_0, tg_xxxxyyy_yyyzz_0, tg_xxxxyyy_yyzzz_0, tg_xxxxyyy_yzzzz_0, \
                                         tg_xxxxyyy_zzzzz_0, tg_xxxxyyz_xxxxx_0, tg_xxxxyyz_xxxxy_0, tg_xxxxyyz_xxxxz_0, \
                                         tg_xxxxyyz_xxxyy_0, tg_xxxxyyz_xxxyz_0, tg_xxxxyyz_xxxzz_0, tg_xxxxyyz_xxyyy_0, \
                                         tg_xxxxyyz_xxyyz_0, tg_xxxxyyz_xxyzz_0, tg_xxxxyyz_xxzzz_0, tg_xxxxyyz_xyyyy_0, \
                                         tg_xxxxyyz_xyyyz_0, tg_xxxxyyz_xyyzz_0, tg_xxxxyyz_xyzzz_0, tg_xxxxyyz_xzzzz_0, \
                                         tg_xxxxyyz_yyyyy_0, tg_xxxxyyz_yyyyz_0, tg_xxxxyyz_yyyzz_0, tg_xxxxyyz_yyzzz_0, \
                                         tg_xxxxyyz_yzzzz_0, tg_xxxxyyz_zzzzz_0, tg_xxxxyz_xyyyz_0, tg_xxxxyz_xyyyz_1, \
                                         tg_xxxxyz_xyyzz_0, tg_xxxxyz_xyyzz_1, tg_xxxxyz_xyzzz_0, tg_xxxxyz_xyzzz_1, \
                                         tg_xxxxyz_xzzzz_0, tg_xxxxyz_xzzzz_1, tg_xxxxyz_yyyyy_0, tg_xxxxyz_yyyyy_1, \
                                         tg_xxxxyz_yyyyz_0, tg_xxxxyz_yyyyz_1, tg_xxxxyz_yyyz_1, tg_xxxxyz_yyyzz_0, \
                                         tg_xxxxyz_yyyzz_1, tg_xxxxyz_yyzz_1, tg_xxxxyz_yyzzz_0, tg_xxxxyz_yyzzz_1, \
                                         tg_xxxxyz_yzzz_1, tg_xxxxyz_yzzzz_0, tg_xxxxyz_yzzzz_1, tg_xxxxyz_zzzz_1, \
                                         tg_xxxxyz_zzzzz_0, tg_xxxxyz_zzzzz_1, tg_xxxxyzz_xxxxx_0, tg_xxxxyzz_xxxxy_0, \
                                         tg_xxxxyzz_xxxxz_0, tg_xxxxyzz_xxxyy_0, tg_xxxxyzz_xxxyz_0, tg_xxxxyzz_xxxzz_0, \
                                         tg_xxxxyzz_xxyyy_0, tg_xxxxyzz_xxyyz_0, tg_xxxxyzz_xxyzz_0, tg_xxxxyzz_xxzzz_0, \
                                         tg_xxxxyzz_xyyyy_0, tg_xxxxyzz_xyyyz_0, tg_xxxxyzz_xyyzz_0, tg_xxxxyzz_xyzzz_0, \
                                         tg_xxxxyzz_xzzzz_0, tg_xxxxyzz_yyyyy_0, tg_xxxxyzz_yyyyz_0, tg_xxxxyzz_yyyzz_0, \
                                         tg_xxxxyzz_yyzzz_0, tg_xxxxyzz_yzzzz_0, tg_xxxxyzz_zzzzz_0, tg_xxxxzz_xxxx_1, \
                                         tg_xxxxzz_xxxxx_0, tg_xxxxzz_xxxxx_1, tg_xxxxzz_xxxxy_0, tg_xxxxzz_xxxxy_1, \
                                         tg_xxxxzz_xxxxz_0, tg_xxxxzz_xxxxz_1, tg_xxxxzz_xxxy_1, tg_xxxxzz_xxxyy_0, \
                                         tg_xxxxzz_xxxyy_1, tg_xxxxzz_xxxyz_0, tg_xxxxzz_xxxyz_1, tg_xxxxzz_xxxz_1, \
                                         tg_xxxxzz_xxxzz_0, tg_xxxxzz_xxxzz_1, tg_xxxxzz_xxyy_1, tg_xxxxzz_xxyyy_0, \
                                         tg_xxxxzz_xxyyy_1, tg_xxxxzz_xxyyz_0, tg_xxxxzz_xxyyz_1, tg_xxxxzz_xxyz_1, \
                                         tg_xxxxzz_xxyzz_0, tg_xxxxzz_xxyzz_1, tg_xxxxzz_xxzz_1, tg_xxxxzz_xxzzz_0, \
                                         tg_xxxxzz_xxzzz_1, tg_xxxxzz_xyyy_1, tg_xxxxzz_xyyyy_0, tg_xxxxzz_xyyyy_1, \
                                         tg_xxxxzz_xyyyz_0, tg_xxxxzz_xyyyz_1, tg_xxxxzz_xyyz_1, tg_xxxxzz_xyyzz_0, \
                                         tg_xxxxzz_xyyzz_1, tg_xxxxzz_xyzz_1, tg_xxxxzz_xyzzz_0, tg_xxxxzz_xyzzz_1, \
                                         tg_xxxxzz_xzzz_1, tg_xxxxzz_xzzzz_0, tg_xxxxzz_xzzzz_1, tg_xxxxzz_yyyy_1, \
                                         tg_xxxxzz_yyyyy_0, tg_xxxxzz_yyyyy_1, tg_xxxxzz_yyyyz_0, tg_xxxxzz_yyyyz_1, \
                                         tg_xxxxzz_yyyz_1, tg_xxxxzz_yyyzz_0, tg_xxxxzz_yyyzz_1, tg_xxxxzz_yyzz_1, \
                                         tg_xxxxzz_yyzzz_0, tg_xxxxzz_yyzzz_1, tg_xxxxzz_yzzz_1, tg_xxxxzz_yzzzz_0, \
                                         tg_xxxxzz_yzzzz_1, tg_xxxxzz_zzzz_1, tg_xxxxzz_zzzzz_0, tg_xxxxzz_zzzzz_1, \
                                         tg_xxxxzzz_xxxxx_0, tg_xxxyyy_xxxx_1, tg_xxxyyy_xxxxx_0, tg_xxxyyy_xxxxx_1, \
                                         tg_xxxyyy_xxxxy_0, tg_xxxyyy_xxxxy_1, tg_xxxyyy_xxxxz_0, tg_xxxyyy_xxxxz_1, \
                                         tg_xxxyyy_xxxy_1, tg_xxxyyy_xxxyy_0, tg_xxxyyy_xxxyy_1, tg_xxxyyy_xxxyz_0, \
                                         tg_xxxyyy_xxxyz_1, tg_xxxyyy_xxxz_1, tg_xxxyyy_xxxzz_0, tg_xxxyyy_xxxzz_1, \
                                         tg_xxxyyy_xxyy_1, tg_xxxyyy_xxyyy_0, tg_xxxyyy_xxyyy_1, tg_xxxyyy_xxyyz_0, \
                                         tg_xxxyyy_xxyyz_1, tg_xxxyyy_xxyz_1, tg_xxxyyy_xxyzz_0, tg_xxxyyy_xxyzz_1, \
                                         tg_xxxyyy_xxzz_1, tg_xxxyyy_xxzzz_0, tg_xxxyyy_xxzzz_1, tg_xxxyyy_xyyy_1, \
                                         tg_xxxyyy_xyyyy_0, tg_xxxyyy_xyyyy_1, tg_xxxyyy_xyyyz_0, tg_xxxyyy_xyyyz_1, \
                                         tg_xxxyyy_xyyz_1, tg_xxxyyy_xyyzz_0, tg_xxxyyy_xyyzz_1, tg_xxxyyy_xyzz_1, \
                                         tg_xxxyyy_xyzzz_0, tg_xxxyyy_xyzzz_1, tg_xxxyyy_xzzz_1, tg_xxxyyy_xzzzz_0, \
                                         tg_xxxyyy_xzzzz_1, tg_xxxyyy_yyyy_1, tg_xxxyyy_yyyyy_0, tg_xxxyyy_yyyyy_1, \
                                         tg_xxxyyy_yyyyz_0, tg_xxxyyy_yyyyz_1, tg_xxxyyy_yyyz_1, tg_xxxyyy_yyyzz_0, \
                                         tg_xxxyyy_yyyzz_1, tg_xxxyyy_yyzz_1, tg_xxxyyy_yyzzz_0, tg_xxxyyy_yyzzz_1, \
                                         tg_xxxyyy_yzzz_1, tg_xxxyyy_yzzzz_0, tg_xxxyyy_yzzzz_1, tg_xxxyyy_zzzz_1, \
                                         tg_xxxyyy_zzzzz_0, tg_xxxyyy_zzzzz_1, tg_xxxyyz_xxxx_1, tg_xxxyyz_xxxxx_0, \
                                         tg_xxxyyz_xxxxx_1, tg_xxxyyz_xxxxy_0, tg_xxxyyz_xxxxy_1, tg_xxxyyz_xxxxz_0, \
                                         tg_xxxyyz_xxxxz_1, tg_xxxyyz_xxxy_1, tg_xxxyyz_xxxyy_0, tg_xxxyyz_xxxyy_1, \
                                         tg_xxxyyz_xxxyz_0, tg_xxxyyz_xxxyz_1, tg_xxxyyz_xxxz_1, tg_xxxyyz_xxxzz_0, \
                                         tg_xxxyyz_xxxzz_1, tg_xxxyyz_xxyy_1, tg_xxxyyz_xxyyy_0, tg_xxxyyz_xxyyy_1, \
                                         tg_xxxyyz_xxyyz_0, tg_xxxyyz_xxyyz_1, tg_xxxyyz_xxyz_1, tg_xxxyyz_xxyzz_0, \
                                         tg_xxxyyz_xxyzz_1, tg_xxxyyz_xxzz_1, tg_xxxyyz_xxzzz_0, tg_xxxyyz_xxzzz_1, \
                                         tg_xxxyyz_xyyy_1, tg_xxxyyz_xyyyy_0, tg_xxxyyz_xyyyy_1, tg_xxxyyz_xyyyz_0, \
                                         tg_xxxyyz_xyyyz_1, tg_xxxyyz_xyyz_1, tg_xxxyyz_xyyzz_0, tg_xxxyyz_xyyzz_1, \
                                         tg_xxxyyz_xyzz_1, tg_xxxyyz_xyzzz_0, tg_xxxyyz_xyzzz_1, tg_xxxyyz_xzzz_1, \
                                         tg_xxxyyz_xzzzz_0, tg_xxxyyz_xzzzz_1, tg_xxxyyz_yyyy_1, tg_xxxyyz_yyyyy_0, \
                                         tg_xxxyyz_yyyyy_1, tg_xxxyyz_yyyyz_0, tg_xxxyyz_yyyyz_1, tg_xxxyyz_yyyz_1, \
                                         tg_xxxyyz_yyyzz_0, tg_xxxyyz_yyyzz_1, tg_xxxyyz_yyzz_1, tg_xxxyyz_yyzzz_0, \
                                         tg_xxxyyz_yyzzz_1, tg_xxxyyz_yzzz_1, tg_xxxyyz_yzzzz_0, tg_xxxyyz_yzzzz_1, \
                                         tg_xxxyyz_zzzz_1, tg_xxxyyz_zzzzz_0, tg_xxxyyz_zzzzz_1, tg_xxxyz_xyyyz_0, \
                                         tg_xxxyz_xyyyz_1, tg_xxxyz_xyyzz_0, tg_xxxyz_xyyzz_1, tg_xxxyz_xyzzz_0, \
                                         tg_xxxyz_xyzzz_1, tg_xxxyz_xzzzz_0, tg_xxxyz_xzzzz_1, tg_xxxyz_yyyyy_0, \
                                         tg_xxxyz_yyyyy_1, tg_xxxyz_yyyyz_0, tg_xxxyz_yyyyz_1, tg_xxxyz_yyyzz_0, \
                                         tg_xxxyz_yyyzz_1, tg_xxxyz_yyzzz_0, tg_xxxyz_yyzzz_1, tg_xxxyz_yzzzz_0, \
                                         tg_xxxyz_yzzzz_1, tg_xxxyz_zzzzz_0, tg_xxxyz_zzzzz_1, tg_xxxyzz_xxxx_1, \
                                         tg_xxxyzz_xxxxx_0, tg_xxxyzz_xxxxx_1, tg_xxxyzz_xxxxy_0, tg_xxxyzz_xxxxy_1, \
                                         tg_xxxyzz_xxxxz_0, tg_xxxyzz_xxxxz_1, tg_xxxyzz_xxxy_1, tg_xxxyzz_xxxyy_0, \
                                         tg_xxxyzz_xxxyy_1, tg_xxxyzz_xxxyz_0, tg_xxxyzz_xxxyz_1, tg_xxxyzz_xxxz_1, \
                                         tg_xxxyzz_xxxzz_0, tg_xxxyzz_xxxzz_1, tg_xxxyzz_xxyy_1, tg_xxxyzz_xxyyy_0, \
                                         tg_xxxyzz_xxyyy_1, tg_xxxyzz_xxyyz_0, tg_xxxyzz_xxyyz_1, tg_xxxyzz_xxyz_1, \
                                         tg_xxxyzz_xxyzz_0, tg_xxxyzz_xxyzz_1, tg_xxxyzz_xxzz_1, tg_xxxyzz_xxzzz_0, \
                                         tg_xxxyzz_xxzzz_1, tg_xxxyzz_xyyy_1, tg_xxxyzz_xyyyy_0, tg_xxxyzz_xyyyy_1, \
                                         tg_xxxyzz_xyyyz_0, tg_xxxyzz_xyyyz_1, tg_xxxyzz_xyyz_1, tg_xxxyzz_xyyzz_0, \
                                         tg_xxxyzz_xyyzz_1, tg_xxxyzz_xyzz_1, tg_xxxyzz_xyzzz_0, tg_xxxyzz_xyzzz_1, \
                                         tg_xxxyzz_xzzz_1, tg_xxxyzz_xzzzz_0, tg_xxxyzz_xzzzz_1, tg_xxxyzz_yyyy_1, \
                                         tg_xxxyzz_yyyyy_0, tg_xxxyzz_yyyyy_1, tg_xxxyzz_yyyyz_0, tg_xxxyzz_yyyyz_1, \
                                         tg_xxxyzz_yyyz_1, tg_xxxyzz_yyyzz_0, tg_xxxyzz_yyyzz_1, tg_xxxyzz_yyzz_1, \
                                         tg_xxxyzz_yyzzz_0, tg_xxxyzz_yyzzz_1, tg_xxxyzz_yzzz_1, tg_xxxyzz_yzzzz_0, \
                                         tg_xxxyzz_yzzzz_1, tg_xxxyzz_zzzz_1, tg_xxxyzz_zzzzz_0, tg_xxxyzz_zzzzz_1, \
                                         tg_xxxzz_xxxxx_0, tg_xxxzz_xxxxx_1, tg_xxxzz_xxxxy_0, tg_xxxzz_xxxxy_1, \
                                         tg_xxxzz_xxxxz_0, tg_xxxzz_xxxxz_1, tg_xxxzz_xxxyy_0, tg_xxxzz_xxxyy_1, \
                                         tg_xxxzz_xxxyz_0, tg_xxxzz_xxxyz_1, tg_xxxzz_xxxzz_0, tg_xxxzz_xxxzz_1, \
                                         tg_xxxzz_xxyyy_0, tg_xxxzz_xxyyy_1, tg_xxxzz_xxyyz_0, tg_xxxzz_xxyyz_1, \
                                         tg_xxxzz_xxyzz_0, tg_xxxzz_xxyzz_1, tg_xxxzz_xxzzz_0, tg_xxxzz_xxzzz_1, \
                                         tg_xxxzz_xyyyy_0, tg_xxxzz_xyyyy_1, tg_xxxzz_xyyyz_0, tg_xxxzz_xyyyz_1, \
                                         tg_xxxzz_xyyzz_0, tg_xxxzz_xyyzz_1, tg_xxxzz_xyzzz_0, tg_xxxzz_xyzzz_1, \
                                         tg_xxxzz_xzzzz_0, tg_xxxzz_xzzzz_1, tg_xxxzz_yyyyy_0, tg_xxxzz_yyyyy_1, \
                                         tg_xxxzz_yyyyz_0, tg_xxxzz_yyyyz_1, tg_xxxzz_yyyzz_0, tg_xxxzz_yyyzz_1, \
                                         tg_xxxzz_yyzzz_0, tg_xxxzz_yyzzz_1, tg_xxxzz_yzzzz_0, tg_xxxzz_yzzzz_1, \
                                         tg_xxxzz_zzzzz_0, tg_xxxzz_zzzzz_1, tg_xxxzzz_xxxx_1, tg_xxxzzz_xxxxx_0, \
                                         tg_xxxzzz_xxxxx_1, tg_xxyyy_xxxxx_0, tg_xxyyy_xxxxx_1, tg_xxyyy_xxxxy_0, \
                                         tg_xxyyy_xxxxy_1, tg_xxyyy_xxxxz_0, tg_xxyyy_xxxxz_1, tg_xxyyy_xxxyy_0, \
                                         tg_xxyyy_xxxyy_1, tg_xxyyy_xxxyz_0, tg_xxyyy_xxxyz_1, tg_xxyyy_xxxzz_0, \
                                         tg_xxyyy_xxxzz_1, tg_xxyyy_xxyyy_0, tg_xxyyy_xxyyy_1, tg_xxyyy_xxyyz_0, \
                                         tg_xxyyy_xxyyz_1, tg_xxyyy_xxyzz_0, tg_xxyyy_xxyzz_1, tg_xxyyy_xxzzz_0, \
                                         tg_xxyyy_xxzzz_1, tg_xxyyy_xyyyy_0, tg_xxyyy_xyyyy_1, tg_xxyyy_xyyyz_0, \
                                         tg_xxyyy_xyyyz_1, tg_xxyyy_xyyzz_0, tg_xxyyy_xyyzz_1, tg_xxyyy_xyzzz_0, \
                                         tg_xxyyy_xyzzz_1, tg_xxyyy_xzzzz_0, tg_xxyyy_xzzzz_1, tg_xxyyy_yyyyy_0, \
                                         tg_xxyyy_yyyyy_1, tg_xxyyy_yyyyz_0, tg_xxyyy_yyyyz_1, tg_xxyyy_yyyzz_0, \
                                         tg_xxyyy_yyyzz_1, tg_xxyyy_yyzzz_0, tg_xxyyy_yyzzz_1, tg_xxyyy_yzzzz_0, \
                                         tg_xxyyy_yzzzz_1, tg_xxyyy_zzzzz_0, tg_xxyyy_zzzzz_1, tg_xxyyz_xxxxx_0, \
                                         tg_xxyyz_xxxxx_1, tg_xxyyz_xxxxy_0, tg_xxyyz_xxxxy_1, tg_xxyyz_xxxxz_0, \
                                         tg_xxyyz_xxxxz_1, tg_xxyyz_xxxyy_0, tg_xxyyz_xxxyy_1, tg_xxyyz_xxxyz_0, \
                                         tg_xxyyz_xxxyz_1, tg_xxyyz_xxxzz_0, tg_xxyyz_xxxzz_1, tg_xxyyz_xxyyy_0, \
                                         tg_xxyyz_xxyyy_1, tg_xxyyz_xxyyz_0, tg_xxyyz_xxyyz_1, tg_xxyyz_xxyzz_0, \
                                         tg_xxyyz_xxyzz_1, tg_xxyyz_xxzzz_0, tg_xxyyz_xxzzz_1, tg_xxyyz_xyyyy_0, \
                                         tg_xxyyz_xyyyy_1, tg_xxyyz_xyyyz_0, tg_xxyyz_xyyyz_1, tg_xxyyz_xyyzz_0, \
                                         tg_xxyyz_xyyzz_1, tg_xxyyz_xyzzz_0, tg_xxyyz_xyzzz_1, tg_xxyyz_xzzzz_0, \
                                         tg_xxyyz_xzzzz_1, tg_xxyyz_yyyyy_0, tg_xxyyz_yyyyy_1, tg_xxyyz_yyyyz_0, \
                                         tg_xxyyz_yyyyz_1, tg_xxyyz_yyyzz_0, tg_xxyyz_yyyzz_1, tg_xxyyz_yyzzz_0, \
                                         tg_xxyyz_yyzzz_1, tg_xxyyz_yzzzz_0, tg_xxyyz_yzzzz_1, tg_xxyyz_zzzzz_0, \
                                         tg_xxyyz_zzzzz_1, tg_xxyzz_xxxxx_0, tg_xxyzz_xxxxx_1, tg_xxyzz_xxxxy_0, \
                                         tg_xxyzz_xxxxy_1, tg_xxyzz_xxxxz_0, tg_xxyzz_xxxxz_1, tg_xxyzz_xxxyy_0, \
                                         tg_xxyzz_xxxyy_1, tg_xxyzz_xxxyz_0, tg_xxyzz_xxxyz_1, tg_xxyzz_xxxzz_0, \
                                         tg_xxyzz_xxxzz_1, tg_xxyzz_xxyyy_0, tg_xxyzz_xxyyy_1, tg_xxyzz_xxyyz_0, \
                                         tg_xxyzz_xxyyz_1, tg_xxyzz_xxyzz_0, tg_xxyzz_xxyzz_1, tg_xxyzz_xxzzz_0, \
                                         tg_xxyzz_xxzzz_1, tg_xxyzz_xyyyy_0, tg_xxyzz_xyyyy_1, tg_xxyzz_xyyyz_0, \
                                         tg_xxyzz_xyyyz_1, tg_xxyzz_xyyzz_0, tg_xxyzz_xyyzz_1, tg_xxyzz_xyzzz_0, \
                                         tg_xxyzz_xyzzz_1, tg_xxyzz_xzzzz_0, tg_xxyzz_xzzzz_1, tg_xxyzz_yyyyy_0, \
                                         tg_xxyzz_yyyyy_1, tg_xxyzz_yyyyz_0, tg_xxyzz_yyyyz_1, tg_xxyzz_yyyzz_0, \
                                         tg_xxyzz_yyyzz_1, tg_xxyzz_yyzzz_0, tg_xxyzz_yyzzz_1, tg_xxyzz_yzzzz_0, \
                                         tg_xxyzz_yzzzz_1, tg_xxyzz_zzzzz_0, tg_xxyzz_zzzzz_1, tg_xxzzz_xxxxx_0, \
                                         tg_xxzzz_xxxxx_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxyz_xyyyz_0[j] = pb_x * tg_xxxxyz_xyyyz_0[j] + fr * tg_xxxxyz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyyz_0[j] - tg_xxxyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyyz_1[j];

                    tg_xxxxxyz_xyyzz_0[j] = pb_x * tg_xxxxyz_xyyzz_0[j] + fr * tg_xxxxyz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyzz_0[j] - tg_xxxyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyzz_1[j];

                    tg_xxxxxyz_xyzzz_0[j] = pb_x * tg_xxxxyz_xyzzz_0[j] + fr * tg_xxxxyz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyzzz_0[j] - tg_xxxyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yzzz_1[j];

                    tg_xxxxxyz_xzzzz_0[j] = pb_x * tg_xxxxyz_xzzzz_0[j] + fr * tg_xxxxyz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xzzzz_0[j] - tg_xxxyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_zzzz_1[j];

                    tg_xxxxxyz_yyyyy_0[j] = pb_x * tg_xxxxyz_yyyyy_0[j] + fr * tg_xxxxyz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyyy_0[j] - tg_xxxyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxyz_yyyyz_0[j] = pb_x * tg_xxxxyz_yyyyz_0[j] + fr * tg_xxxxyz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyyz_0[j] - tg_xxxyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxyz_yyyzz_0[j] = pb_x * tg_xxxxyz_yyyzz_0[j] + fr * tg_xxxxyz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyzz_0[j] - tg_xxxyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxyz_yyzzz_0[j] = pb_x * tg_xxxxyz_yyzzz_0[j] + fr * tg_xxxxyz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyzzz_0[j] - tg_xxxyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_yzzzz_0[j] = pb_x * tg_xxxxyz_yzzzz_0[j] + fr * tg_xxxxyz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yzzzz_0[j] - tg_xxxyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_zzzzz_0[j] = pb_x * tg_xxxxyz_zzzzz_0[j] + fr * tg_xxxxyz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_zzzzz_0[j] - tg_xxxyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_xxxxx_0[j] = pb_x * tg_xxxxzz_xxxxx_0[j] + fr * tg_xxxxzz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxx_0[j] - tg_xxxzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxzz_xxxx_1[j];

                    tg_xxxxxzz_xxxxy_0[j] = pb_x * tg_xxxxzz_xxxxy_0[j] + fr * tg_xxxxzz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxy_0[j] - tg_xxxzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzz_xxxy_1[j];

                    tg_xxxxxzz_xxxxz_0[j] = pb_x * tg_xxxxzz_xxxxz_0[j] + fr * tg_xxxxzz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxz_0[j] - tg_xxxzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzz_xxxz_1[j];

                    tg_xxxxxzz_xxxyy_0[j] = pb_x * tg_xxxxzz_xxxyy_0[j] + fr * tg_xxxxzz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxyy_0[j] - tg_xxxzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxyy_1[j];

                    tg_xxxxxzz_xxxyz_0[j] = pb_x * tg_xxxxzz_xxxyz_0[j] + fr * tg_xxxxzz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxyz_0[j] - tg_xxxzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxyz_1[j];

                    tg_xxxxxzz_xxxzz_0[j] = pb_x * tg_xxxxzz_xxxzz_0[j] + fr * tg_xxxxzz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxzz_0[j] - tg_xxxzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxzz_1[j];

                    tg_xxxxxzz_xxyyy_0[j] = pb_x * tg_xxxxzz_xxyyy_0[j] + fr * tg_xxxxzz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyyy_0[j] - tg_xxxzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyyy_1[j];

                    tg_xxxxxzz_xxyyz_0[j] = pb_x * tg_xxxxzz_xxyyz_0[j] + fr * tg_xxxxzz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyyz_0[j] - tg_xxxzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyyz_1[j];

                    tg_xxxxxzz_xxyzz_0[j] = pb_x * tg_xxxxzz_xxyzz_0[j] + fr * tg_xxxxzz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyzz_0[j] - tg_xxxzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyzz_1[j];

                    tg_xxxxxzz_xxzzz_0[j] = pb_x * tg_xxxxzz_xxzzz_0[j] + fr * tg_xxxxzz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxzzz_0[j] - tg_xxxzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xzzz_1[j];

                    tg_xxxxxzz_xyyyy_0[j] = pb_x * tg_xxxxzz_xyyyy_0[j] + fr * tg_xxxxzz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyyy_0[j] - tg_xxxzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyyy_1[j];

                    tg_xxxxxzz_xyyyz_0[j] = pb_x * tg_xxxxzz_xyyyz_0[j] + fr * tg_xxxxzz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyyz_0[j] - tg_xxxzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyyz_1[j];

                    tg_xxxxxzz_xyyzz_0[j] = pb_x * tg_xxxxzz_xyyzz_0[j] + fr * tg_xxxxzz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyzz_0[j] - tg_xxxzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyzz_1[j];

                    tg_xxxxxzz_xyzzz_0[j] = pb_x * tg_xxxxzz_xyzzz_0[j] + fr * tg_xxxxzz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyzzz_0[j] - tg_xxxzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yzzz_1[j];

                    tg_xxxxxzz_xzzzz_0[j] = pb_x * tg_xxxxzz_xzzzz_0[j] + fr * tg_xxxxzz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xzzzz_0[j] - tg_xxxzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_zzzz_1[j];

                    tg_xxxxxzz_yyyyy_0[j] = pb_x * tg_xxxxzz_yyyyy_0[j] + fr * tg_xxxxzz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyyy_0[j] - tg_xxxzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxzz_yyyyz_0[j] = pb_x * tg_xxxxzz_yyyyz_0[j] + fr * tg_xxxxzz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyyz_0[j] - tg_xxxzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxzz_yyyzz_0[j] = pb_x * tg_xxxxzz_yyyzz_0[j] + fr * tg_xxxxzz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyzz_0[j] - tg_xxxzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxzz_yyzzz_0[j] = pb_x * tg_xxxxzz_yyzzz_0[j] + fr * tg_xxxxzz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyzzz_0[j] - tg_xxxzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_yzzzz_0[j] = pb_x * tg_xxxxzz_yzzzz_0[j] + fr * tg_xxxxzz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yzzzz_0[j] - tg_xxxzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_zzzzz_0[j] = pb_x * tg_xxxxzz_zzzzz_0[j] + fr * tg_xxxxzz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_zzzzz_0[j] - tg_xxxzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_xxxxx_0[j] = pb_x * tg_xxxyyy_xxxxx_0[j] + fr * tg_xxxyyy_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxx_0[j] - tg_xxyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyy_xxxx_1[j];

                    tg_xxxxyyy_xxxxy_0[j] = pb_x * tg_xxxyyy_xxxxy_0[j] + fr * tg_xxxyyy_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxy_0[j] - tg_xxyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyy_xxxy_1[j];

                    tg_xxxxyyy_xxxxz_0[j] = pb_x * tg_xxxyyy_xxxxz_0[j] + fr * tg_xxxyyy_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxz_0[j] - tg_xxyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyy_xxxz_1[j];

                    tg_xxxxyyy_xxxyy_0[j] = pb_x * tg_xxxyyy_xxxyy_0[j] + fr * tg_xxxyyy_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxyy_0[j] - tg_xxyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxyy_1[j];

                    tg_xxxxyyy_xxxyz_0[j] = pb_x * tg_xxxyyy_xxxyz_0[j] + fr * tg_xxxyyy_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxyz_0[j] - tg_xxyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxyz_1[j];

                    tg_xxxxyyy_xxxzz_0[j] = pb_x * tg_xxxyyy_xxxzz_0[j] + fr * tg_xxxyyy_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxzz_0[j] - tg_xxyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxzz_1[j];

                    tg_xxxxyyy_xxyyy_0[j] = pb_x * tg_xxxyyy_xxyyy_0[j] + fr * tg_xxxyyy_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyyy_0[j] - tg_xxyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyyy_1[j];

                    tg_xxxxyyy_xxyyz_0[j] = pb_x * tg_xxxyyy_xxyyz_0[j] + fr * tg_xxxyyy_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyyz_0[j] - tg_xxyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyyz_1[j];

                    tg_xxxxyyy_xxyzz_0[j] = pb_x * tg_xxxyyy_xxyzz_0[j] + fr * tg_xxxyyy_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyzz_0[j] - tg_xxyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyzz_1[j];

                    tg_xxxxyyy_xxzzz_0[j] = pb_x * tg_xxxyyy_xxzzz_0[j] + fr * tg_xxxyyy_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxzzz_0[j] - tg_xxyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xzzz_1[j];

                    tg_xxxxyyy_xyyyy_0[j] = pb_x * tg_xxxyyy_xyyyy_0[j] + fr * tg_xxxyyy_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyyy_0[j] - tg_xxyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyyy_1[j];

                    tg_xxxxyyy_xyyyz_0[j] = pb_x * tg_xxxyyy_xyyyz_0[j] + fr * tg_xxxyyy_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyyz_0[j] - tg_xxyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyyz_1[j];

                    tg_xxxxyyy_xyyzz_0[j] = pb_x * tg_xxxyyy_xyyzz_0[j] + fr * tg_xxxyyy_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyzz_0[j] - tg_xxyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyzz_1[j];

                    tg_xxxxyyy_xyzzz_0[j] = pb_x * tg_xxxyyy_xyzzz_0[j] + fr * tg_xxxyyy_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyzzz_0[j] - tg_xxyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yzzz_1[j];

                    tg_xxxxyyy_xzzzz_0[j] = pb_x * tg_xxxyyy_xzzzz_0[j] + fr * tg_xxxyyy_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xzzzz_0[j] - tg_xxyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_zzzz_1[j];

                    tg_xxxxyyy_yyyyy_0[j] = pb_x * tg_xxxyyy_yyyyy_0[j] + fr * tg_xxxyyy_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyyy_0[j] - tg_xxyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyyy_yyyyz_0[j] = pb_x * tg_xxxyyy_yyyyz_0[j] + fr * tg_xxxyyy_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyyz_0[j] - tg_xxyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyyy_yyyzz_0[j] = pb_x * tg_xxxyyy_yyyzz_0[j] + fr * tg_xxxyyy_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyzz_0[j] - tg_xxyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyyy_yyzzz_0[j] = pb_x * tg_xxxyyy_yyzzz_0[j] + fr * tg_xxxyyy_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyzzz_0[j] - tg_xxyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_yzzzz_0[j] = pb_x * tg_xxxyyy_yzzzz_0[j] + fr * tg_xxxyyy_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yzzzz_0[j] - tg_xxyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_zzzzz_0[j] = pb_x * tg_xxxyyy_zzzzz_0[j] + fr * tg_xxxyyy_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_zzzzz_0[j] - tg_xxyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_xxxxx_0[j] = pb_x * tg_xxxyyz_xxxxx_0[j] + fr * tg_xxxyyz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxx_0[j] - tg_xxyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyz_xxxx_1[j];

                    tg_xxxxyyz_xxxxy_0[j] = pb_x * tg_xxxyyz_xxxxy_0[j] + fr * tg_xxxyyz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxy_0[j] - tg_xxyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyz_xxxy_1[j];

                    tg_xxxxyyz_xxxxz_0[j] = pb_x * tg_xxxyyz_xxxxz_0[j] + fr * tg_xxxyyz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxz_0[j] - tg_xxyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyz_xxxz_1[j];

                    tg_xxxxyyz_xxxyy_0[j] = pb_x * tg_xxxyyz_xxxyy_0[j] + fr * tg_xxxyyz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxyy_0[j] - tg_xxyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxyy_1[j];

                    tg_xxxxyyz_xxxyz_0[j] = pb_x * tg_xxxyyz_xxxyz_0[j] + fr * tg_xxxyyz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxyz_0[j] - tg_xxyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxyz_1[j];

                    tg_xxxxyyz_xxxzz_0[j] = pb_x * tg_xxxyyz_xxxzz_0[j] + fr * tg_xxxyyz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxzz_0[j] - tg_xxyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxzz_1[j];

                    tg_xxxxyyz_xxyyy_0[j] = pb_x * tg_xxxyyz_xxyyy_0[j] + fr * tg_xxxyyz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyyy_0[j] - tg_xxyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyyy_1[j];

                    tg_xxxxyyz_xxyyz_0[j] = pb_x * tg_xxxyyz_xxyyz_0[j] + fr * tg_xxxyyz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyyz_0[j] - tg_xxyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyyz_1[j];

                    tg_xxxxyyz_xxyzz_0[j] = pb_x * tg_xxxyyz_xxyzz_0[j] + fr * tg_xxxyyz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyzz_0[j] - tg_xxyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyzz_1[j];

                    tg_xxxxyyz_xxzzz_0[j] = pb_x * tg_xxxyyz_xxzzz_0[j] + fr * tg_xxxyyz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxzzz_0[j] - tg_xxyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xzzz_1[j];

                    tg_xxxxyyz_xyyyy_0[j] = pb_x * tg_xxxyyz_xyyyy_0[j] + fr * tg_xxxyyz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyyy_0[j] - tg_xxyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyyy_1[j];

                    tg_xxxxyyz_xyyyz_0[j] = pb_x * tg_xxxyyz_xyyyz_0[j] + fr * tg_xxxyyz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyyz_0[j] - tg_xxyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyyz_1[j];

                    tg_xxxxyyz_xyyzz_0[j] = pb_x * tg_xxxyyz_xyyzz_0[j] + fr * tg_xxxyyz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyzz_0[j] - tg_xxyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyzz_1[j];

                    tg_xxxxyyz_xyzzz_0[j] = pb_x * tg_xxxyyz_xyzzz_0[j] + fr * tg_xxxyyz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyzzz_0[j] - tg_xxyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yzzz_1[j];

                    tg_xxxxyyz_xzzzz_0[j] = pb_x * tg_xxxyyz_xzzzz_0[j] + fr * tg_xxxyyz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xzzzz_0[j] - tg_xxyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_zzzz_1[j];

                    tg_xxxxyyz_yyyyy_0[j] = pb_x * tg_xxxyyz_yyyyy_0[j] + fr * tg_xxxyyz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyyy_0[j] - tg_xxyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyyz_yyyyz_0[j] = pb_x * tg_xxxyyz_yyyyz_0[j] + fr * tg_xxxyyz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyyz_0[j] - tg_xxyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyyz_yyyzz_0[j] = pb_x * tg_xxxyyz_yyyzz_0[j] + fr * tg_xxxyyz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyzz_0[j] - tg_xxyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyyz_yyzzz_0[j] = pb_x * tg_xxxyyz_yyzzz_0[j] + fr * tg_xxxyyz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyzzz_0[j] - tg_xxyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_yzzzz_0[j] = pb_x * tg_xxxyyz_yzzzz_0[j] + fr * tg_xxxyyz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yzzzz_0[j] - tg_xxyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_zzzzz_0[j] = pb_x * tg_xxxyyz_zzzzz_0[j] + fr * tg_xxxyyz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_zzzzz_0[j] - tg_xxyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_xxxxx_0[j] = pb_x * tg_xxxyzz_xxxxx_0[j] + fr * tg_xxxyzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxx_0[j] - tg_xxyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyzz_xxxx_1[j];

                    tg_xxxxyzz_xxxxy_0[j] = pb_x * tg_xxxyzz_xxxxy_0[j] + fr * tg_xxxyzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxy_0[j] - tg_xxyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzz_xxxy_1[j];

                    tg_xxxxyzz_xxxxz_0[j] = pb_x * tg_xxxyzz_xxxxz_0[j] + fr * tg_xxxyzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxz_0[j] - tg_xxyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzz_xxxz_1[j];

                    tg_xxxxyzz_xxxyy_0[j] = pb_x * tg_xxxyzz_xxxyy_0[j] + fr * tg_xxxyzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxyy_0[j] - tg_xxyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxyy_1[j];

                    tg_xxxxyzz_xxxyz_0[j] = pb_x * tg_xxxyzz_xxxyz_0[j] + fr * tg_xxxyzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxyz_0[j] - tg_xxyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxyz_1[j];

                    tg_xxxxyzz_xxxzz_0[j] = pb_x * tg_xxxyzz_xxxzz_0[j] + fr * tg_xxxyzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxzz_0[j] - tg_xxyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxzz_1[j];

                    tg_xxxxyzz_xxyyy_0[j] = pb_x * tg_xxxyzz_xxyyy_0[j] + fr * tg_xxxyzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyyy_0[j] - tg_xxyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyyy_1[j];

                    tg_xxxxyzz_xxyyz_0[j] = pb_x * tg_xxxyzz_xxyyz_0[j] + fr * tg_xxxyzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyyz_0[j] - tg_xxyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyyz_1[j];

                    tg_xxxxyzz_xxyzz_0[j] = pb_x * tg_xxxyzz_xxyzz_0[j] + fr * tg_xxxyzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyzz_0[j] - tg_xxyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyzz_1[j];

                    tg_xxxxyzz_xxzzz_0[j] = pb_x * tg_xxxyzz_xxzzz_0[j] + fr * tg_xxxyzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxzzz_0[j] - tg_xxyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xzzz_1[j];

                    tg_xxxxyzz_xyyyy_0[j] = pb_x * tg_xxxyzz_xyyyy_0[j] + fr * tg_xxxyzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyyy_0[j] - tg_xxyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyyy_1[j];

                    tg_xxxxyzz_xyyyz_0[j] = pb_x * tg_xxxyzz_xyyyz_0[j] + fr * tg_xxxyzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyyz_0[j] - tg_xxyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyyz_1[j];

                    tg_xxxxyzz_xyyzz_0[j] = pb_x * tg_xxxyzz_xyyzz_0[j] + fr * tg_xxxyzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyzz_0[j] - tg_xxyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyzz_1[j];

                    tg_xxxxyzz_xyzzz_0[j] = pb_x * tg_xxxyzz_xyzzz_0[j] + fr * tg_xxxyzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyzzz_0[j] - tg_xxyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yzzz_1[j];

                    tg_xxxxyzz_xzzzz_0[j] = pb_x * tg_xxxyzz_xzzzz_0[j] + fr * tg_xxxyzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xzzzz_0[j] - tg_xxyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_zzzz_1[j];

                    tg_xxxxyzz_yyyyy_0[j] = pb_x * tg_xxxyzz_yyyyy_0[j] + fr * tg_xxxyzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyyy_0[j] - tg_xxyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyzz_yyyyz_0[j] = pb_x * tg_xxxyzz_yyyyz_0[j] + fr * tg_xxxyzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyyz_0[j] - tg_xxyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyzz_yyyzz_0[j] = pb_x * tg_xxxyzz_yyyzz_0[j] + fr * tg_xxxyzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyzz_0[j] - tg_xxyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyzz_yyzzz_0[j] = pb_x * tg_xxxyzz_yyzzz_0[j] + fr * tg_xxxyzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyzzz_0[j] - tg_xxyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_yzzzz_0[j] = pb_x * tg_xxxyzz_yzzzz_0[j] + fr * tg_xxxyzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yzzzz_0[j] - tg_xxyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_zzzzz_0[j] = pb_x * tg_xxxyzz_zzzzz_0[j] + fr * tg_xxxyzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_zzzzz_0[j] - tg_xxyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_xxxxx_0[j] = pb_x * tg_xxxzzz_xxxxx_0[j] + fr * tg_xxxzzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxx_0[j] - tg_xxzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxzzz_xxxx_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSH_190_285(      CMemBlock2D<double>* primBuffer,
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
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_xxxzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 190); 

                auto tg_xxxzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 191); 

                auto tg_xxxzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 192); 

                auto tg_xxxzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 193); 

                auto tg_xxxzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 194); 

                auto tg_xxxzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 195); 

                auto tg_xxxzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 196); 

                auto tg_xxxzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 197); 

                auto tg_xxxzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 198); 

                auto tg_xxxzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 199); 

                auto tg_xxxzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 200); 

                auto tg_xxxzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 201); 

                auto tg_xxxzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 202); 

                auto tg_xxxzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 203); 

                auto tg_xxxzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 204); 

                auto tg_xxxzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 205); 

                auto tg_xxxzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 206); 

                auto tg_xxxzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 207); 

                auto tg_xxxzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 208); 

                auto tg_xxxzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 209); 

                auto tg_xxyyyy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 210); 

                auto tg_xxyyyy_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 211); 

                auto tg_xxyyyy_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 212); 

                auto tg_xxyyyy_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 213); 

                auto tg_xxyyyy_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 214); 

                auto tg_xxyyyy_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 215); 

                auto tg_xxyyyy_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 216); 

                auto tg_xxyyyy_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 217); 

                auto tg_xxyyyy_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 218); 

                auto tg_xxyyyy_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 219); 

                auto tg_xxyyyy_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 220); 

                auto tg_xxyyyy_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 221); 

                auto tg_xxyyyy_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 222); 

                auto tg_xxyyyy_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 223); 

                auto tg_xxyyyy_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 224); 

                auto tg_xxyyyy_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 225); 

                auto tg_xxyyyy_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 226); 

                auto tg_xxyyyy_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 227); 

                auto tg_xxyyyy_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 228); 

                auto tg_xxyyyy_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 229); 

                auto tg_xxyyyy_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 230); 

                auto tg_xxyyyz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 231); 

                auto tg_xxyyyz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 232); 

                auto tg_xxyyyz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 233); 

                auto tg_xxyyyz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 234); 

                auto tg_xxyyyz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 235); 

                auto tg_xxyyyz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 236); 

                auto tg_xxyyyz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 237); 

                auto tg_xxyyyz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 238); 

                auto tg_xxyyyz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 239); 

                auto tg_xxyyyz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 240); 

                auto tg_xxyyyz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 241); 

                auto tg_xxyyyz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 242); 

                auto tg_xxyyyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 243); 

                auto tg_xxyyyz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 244); 

                auto tg_xxyyyz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 245); 

                auto tg_xxyyyz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 246); 

                auto tg_xxyyyz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 247); 

                auto tg_xxyyyz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 248); 

                auto tg_xxyyyz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 249); 

                auto tg_xxyyyz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 250); 

                auto tg_xxyyyz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 251); 

                auto tg_xxyyzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 252); 

                auto tg_xxyyzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 253); 

                auto tg_xxyyzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 254); 

                auto tg_xxyyzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 255); 

                auto tg_xxyyzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 256); 

                auto tg_xxyyzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 257); 

                auto tg_xxyyzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 258); 

                auto tg_xxyyzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 259); 

                auto tg_xxyyzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 260); 

                auto tg_xxyyzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 261); 

                auto tg_xxyyzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 262); 

                auto tg_xxyyzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 263); 

                auto tg_xxyyzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 264); 

                auto tg_xxyyzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 265); 

                auto tg_xxyyzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 266); 

                auto tg_xxyyzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 267); 

                auto tg_xxyyzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 268); 

                auto tg_xxyyzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 269); 

                auto tg_xxyyzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 270); 

                auto tg_xxyyzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 271); 

                auto tg_xxyyzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 272); 

                auto tg_xxyzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 273); 

                auto tg_xxyzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 274); 

                auto tg_xxyzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 275); 

                auto tg_xxyzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 276); 

                auto tg_xxyzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 277); 

                auto tg_xxyzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 278); 

                auto tg_xxyzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 279); 

                auto tg_xxyzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 280); 

                auto tg_xxyzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 281); 

                auto tg_xxyzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 282); 

                auto tg_xxyzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 283); 

                auto tg_xxyzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 284); 

                auto tg_xxxzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 190); 

                auto tg_xxxzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 191); 

                auto tg_xxxzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 192); 

                auto tg_xxxzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 193); 

                auto tg_xxxzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 194); 

                auto tg_xxxzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 195); 

                auto tg_xxxzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 196); 

                auto tg_xxxzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 197); 

                auto tg_xxxzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 198); 

                auto tg_xxxzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 199); 

                auto tg_xxxzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 200); 

                auto tg_xxxzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 201); 

                auto tg_xxxzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 202); 

                auto tg_xxxzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 203); 

                auto tg_xxxzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 204); 

                auto tg_xxxzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 205); 

                auto tg_xxxzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 206); 

                auto tg_xxxzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 207); 

                auto tg_xxxzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 208); 

                auto tg_xxxzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 209); 

                auto tg_xxyyyy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 210); 

                auto tg_xxyyyy_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 211); 

                auto tg_xxyyyy_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 212); 

                auto tg_xxyyyy_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 213); 

                auto tg_xxyyyy_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 214); 

                auto tg_xxyyyy_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 215); 

                auto tg_xxyyyy_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 216); 

                auto tg_xxyyyy_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 217); 

                auto tg_xxyyyy_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 218); 

                auto tg_xxyyyy_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 219); 

                auto tg_xxyyyy_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 220); 

                auto tg_xxyyyy_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 221); 

                auto tg_xxyyyy_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 222); 

                auto tg_xxyyyy_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 223); 

                auto tg_xxyyyy_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 224); 

                auto tg_xxyyyy_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 225); 

                auto tg_xxyyyy_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 226); 

                auto tg_xxyyyy_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 227); 

                auto tg_xxyyyy_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 228); 

                auto tg_xxyyyy_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 229); 

                auto tg_xxyyyy_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 230); 

                auto tg_xxyyyz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 231); 

                auto tg_xxyyyz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 232); 

                auto tg_xxyyyz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 233); 

                auto tg_xxyyyz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 234); 

                auto tg_xxyyyz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 235); 

                auto tg_xxyyyz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 236); 

                auto tg_xxyyyz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 237); 

                auto tg_xxyyyz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 238); 

                auto tg_xxyyyz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 239); 

                auto tg_xxyyyz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 240); 

                auto tg_xxyyyz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 241); 

                auto tg_xxyyyz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 242); 

                auto tg_xxyyyz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 243); 

                auto tg_xxyyyz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 244); 

                auto tg_xxyyyz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 245); 

                auto tg_xxyyyz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 246); 

                auto tg_xxyyyz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 247); 

                auto tg_xxyyyz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 248); 

                auto tg_xxyyyz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 249); 

                auto tg_xxyyyz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 250); 

                auto tg_xxyyyz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 251); 

                auto tg_xxyyzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 252); 

                auto tg_xxyyzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 253); 

                auto tg_xxyyzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 254); 

                auto tg_xxyyzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 255); 

                auto tg_xxyyzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 256); 

                auto tg_xxyyzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 257); 

                auto tg_xxyyzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 258); 

                auto tg_xxyyzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 259); 

                auto tg_xxyyzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 260); 

                auto tg_xxyyzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 261); 

                auto tg_xxyyzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 262); 

                auto tg_xxyyzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 263); 

                auto tg_xxyyzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 264); 

                auto tg_xxyyzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 265); 

                auto tg_xxyyzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 266); 

                auto tg_xxyyzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 267); 

                auto tg_xxyyzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 268); 

                auto tg_xxyyzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 269); 

                auto tg_xxyyzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 270); 

                auto tg_xxyyzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 271); 

                auto tg_xxyyzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 272); 

                auto tg_xxyzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 273); 

                auto tg_xxyzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 274); 

                auto tg_xxyzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 275); 

                auto tg_xxyzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 276); 

                auto tg_xxyzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 277); 

                auto tg_xxyzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 278); 

                auto tg_xxyzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 279); 

                auto tg_xxyzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 280); 

                auto tg_xxyzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 281); 

                auto tg_xxyzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 282); 

                auto tg_xxyzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 283); 

                auto tg_xxyzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 284); 

                auto tg_xxzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 190); 

                auto tg_xxzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 191); 

                auto tg_xxzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 192); 

                auto tg_xxzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 193); 

                auto tg_xxzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 194); 

                auto tg_xxzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 195); 

                auto tg_xxzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 196); 

                auto tg_xxzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 197); 

                auto tg_xxzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 198); 

                auto tg_xxzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 199); 

                auto tg_xxzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 200); 

                auto tg_xxzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 201); 

                auto tg_xxzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 202); 

                auto tg_xxzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 203); 

                auto tg_xxzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 204); 

                auto tg_xxzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 205); 

                auto tg_xxzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 206); 

                auto tg_xxzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 207); 

                auto tg_xxzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 208); 

                auto tg_xxzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 209); 

                auto tg_xyyyy_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 210); 

                auto tg_xyyyy_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 211); 

                auto tg_xyyyy_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 212); 

                auto tg_xyyyy_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 213); 

                auto tg_xyyyy_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 214); 

                auto tg_xyyyy_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 215); 

                auto tg_xyyyy_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 216); 

                auto tg_xyyyy_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 217); 

                auto tg_xyyyy_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 218); 

                auto tg_xyyyy_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 219); 

                auto tg_xyyyy_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 220); 

                auto tg_xyyyy_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 221); 

                auto tg_xyyyy_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 222); 

                auto tg_xyyyy_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 223); 

                auto tg_xyyyy_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 224); 

                auto tg_xyyyy_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 225); 

                auto tg_xyyyy_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 226); 

                auto tg_xyyyy_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 227); 

                auto tg_xyyyy_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 228); 

                auto tg_xyyyy_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 229); 

                auto tg_xyyyy_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 230); 

                auto tg_xyyyz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 231); 

                auto tg_xyyyz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 232); 

                auto tg_xyyyz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 233); 

                auto tg_xyyyz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 234); 

                auto tg_xyyyz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 235); 

                auto tg_xyyyz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 236); 

                auto tg_xyyyz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 237); 

                auto tg_xyyyz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 238); 

                auto tg_xyyyz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 239); 

                auto tg_xyyyz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 240); 

                auto tg_xyyyz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 241); 

                auto tg_xyyyz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 242); 

                auto tg_xyyyz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 243); 

                auto tg_xyyyz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 244); 

                auto tg_xyyyz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 245); 

                auto tg_xyyyz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 246); 

                auto tg_xyyyz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 247); 

                auto tg_xyyyz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 248); 

                auto tg_xyyyz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 249); 

                auto tg_xyyyz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 250); 

                auto tg_xyyyz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 251); 

                auto tg_xyyzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 252); 

                auto tg_xyyzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 253); 

                auto tg_xyyzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 254); 

                auto tg_xyyzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 255); 

                auto tg_xyyzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 256); 

                auto tg_xyyzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 257); 

                auto tg_xyyzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 258); 

                auto tg_xyyzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 259); 

                auto tg_xyyzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 260); 

                auto tg_xyyzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 261); 

                auto tg_xyyzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 262); 

                auto tg_xyyzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 263); 

                auto tg_xyyzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 264); 

                auto tg_xyyzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 265); 

                auto tg_xyyzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 266); 

                auto tg_xyyzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 267); 

                auto tg_xyyzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 268); 

                auto tg_xyyzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 269); 

                auto tg_xyyzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 270); 

                auto tg_xyyzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 271); 

                auto tg_xyyzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 272); 

                auto tg_xyzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 273); 

                auto tg_xyzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 274); 

                auto tg_xyzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 275); 

                auto tg_xyzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 276); 

                auto tg_xyzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 277); 

                auto tg_xyzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 278); 

                auto tg_xyzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 279); 

                auto tg_xyzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 280); 

                auto tg_xyzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 281); 

                auto tg_xyzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 282); 

                auto tg_xyzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 283); 

                auto tg_xyzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 284); 

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

                auto tg_xxxzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 136); 

                auto tg_xxxzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 137); 

                auto tg_xxxzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 138); 

                auto tg_xxxzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 139); 

                auto tg_xxxzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 140); 

                auto tg_xxxzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 141); 

                auto tg_xxxzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 142); 

                auto tg_xxxzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 143); 

                auto tg_xxxzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 144); 

                auto tg_xxxzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 145); 

                auto tg_xxxzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 146); 

                auto tg_xxxzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 147); 

                auto tg_xxxzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 148); 

                auto tg_xxxzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 149); 

                auto tg_xxyyyy_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 150); 

                auto tg_xxyyyy_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 151); 

                auto tg_xxyyyy_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 152); 

                auto tg_xxyyyy_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 153); 

                auto tg_xxyyyy_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 154); 

                auto tg_xxyyyy_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 155); 

                auto tg_xxyyyy_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 156); 

                auto tg_xxyyyy_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 157); 

                auto tg_xxyyyy_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 158); 

                auto tg_xxyyyy_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 159); 

                auto tg_xxyyyy_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 160); 

                auto tg_xxyyyy_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 161); 

                auto tg_xxyyyy_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 162); 

                auto tg_xxyyyy_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 163); 

                auto tg_xxyyyy_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 164); 

                auto tg_xxyyyz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 165); 

                auto tg_xxyyyz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 166); 

                auto tg_xxyyyz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 167); 

                auto tg_xxyyyz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 168); 

                auto tg_xxyyyz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 169); 

                auto tg_xxyyyz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 170); 

                auto tg_xxyyyz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 171); 

                auto tg_xxyyyz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 172); 

                auto tg_xxyyyz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 173); 

                auto tg_xxyyyz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 174); 

                auto tg_xxyyyz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 175); 

                auto tg_xxyyyz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 176); 

                auto tg_xxyyyz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 177); 

                auto tg_xxyyyz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 178); 

                auto tg_xxyyyz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 179); 

                auto tg_xxyyzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 180); 

                auto tg_xxyyzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 181); 

                auto tg_xxyyzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 182); 

                auto tg_xxyyzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 183); 

                auto tg_xxyyzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 184); 

                auto tg_xxyyzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 185); 

                auto tg_xxyyzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 186); 

                auto tg_xxyyzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 187); 

                auto tg_xxyyzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 188); 

                auto tg_xxyyzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 189); 

                auto tg_xxyyzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 190); 

                auto tg_xxyyzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 191); 

                auto tg_xxyyzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 192); 

                auto tg_xxyyzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 193); 

                auto tg_xxyyzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 194); 

                auto tg_xxyzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 195); 

                auto tg_xxyzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 196); 

                auto tg_xxyzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 197); 

                auto tg_xxyzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 198); 

                auto tg_xxyzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 199); 

                auto tg_xxyzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 200); 

                auto tg_xxyzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 201); 

                auto tg_xxyzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 202); 

                auto tg_xxyzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 203); 

                auto tg_xxyzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 204); 

                auto tg_xxyzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 205); 

                auto tg_xxyzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 206); 

                // set up pointers to integrals

                auto tg_xxxxzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 190); 

                auto tg_xxxxzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 191); 

                auto tg_xxxxzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 192); 

                auto tg_xxxxzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 193); 

                auto tg_xxxxzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 194); 

                auto tg_xxxxzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 195); 

                auto tg_xxxxzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 196); 

                auto tg_xxxxzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 197); 

                auto tg_xxxxzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 198); 

                auto tg_xxxxzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 199); 

                auto tg_xxxxzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 200); 

                auto tg_xxxxzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 201); 

                auto tg_xxxxzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 202); 

                auto tg_xxxxzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 203); 

                auto tg_xxxxzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 204); 

                auto tg_xxxxzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 205); 

                auto tg_xxxxzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 206); 

                auto tg_xxxxzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 207); 

                auto tg_xxxxzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 208); 

                auto tg_xxxxzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 209); 

                auto tg_xxxyyyy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 210); 

                auto tg_xxxyyyy_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 211); 

                auto tg_xxxyyyy_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 212); 

                auto tg_xxxyyyy_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 213); 

                auto tg_xxxyyyy_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 214); 

                auto tg_xxxyyyy_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 215); 

                auto tg_xxxyyyy_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 216); 

                auto tg_xxxyyyy_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 217); 

                auto tg_xxxyyyy_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 218); 

                auto tg_xxxyyyy_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 219); 

                auto tg_xxxyyyy_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 220); 

                auto tg_xxxyyyy_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 221); 

                auto tg_xxxyyyy_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 222); 

                auto tg_xxxyyyy_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 223); 

                auto tg_xxxyyyy_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 224); 

                auto tg_xxxyyyy_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 225); 

                auto tg_xxxyyyy_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 226); 

                auto tg_xxxyyyy_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 227); 

                auto tg_xxxyyyy_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 228); 

                auto tg_xxxyyyy_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 229); 

                auto tg_xxxyyyy_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 230); 

                auto tg_xxxyyyz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 231); 

                auto tg_xxxyyyz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 232); 

                auto tg_xxxyyyz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 233); 

                auto tg_xxxyyyz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 234); 

                auto tg_xxxyyyz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 235); 

                auto tg_xxxyyyz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 236); 

                auto tg_xxxyyyz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 237); 

                auto tg_xxxyyyz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 238); 

                auto tg_xxxyyyz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 239); 

                auto tg_xxxyyyz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 240); 

                auto tg_xxxyyyz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 241); 

                auto tg_xxxyyyz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 242); 

                auto tg_xxxyyyz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 243); 

                auto tg_xxxyyyz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 244); 

                auto tg_xxxyyyz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 245); 

                auto tg_xxxyyyz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 246); 

                auto tg_xxxyyyz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 247); 

                auto tg_xxxyyyz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 248); 

                auto tg_xxxyyyz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 249); 

                auto tg_xxxyyyz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 250); 

                auto tg_xxxyyyz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 251); 

                auto tg_xxxyyzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 252); 

                auto tg_xxxyyzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 253); 

                auto tg_xxxyyzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 254); 

                auto tg_xxxyyzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 255); 

                auto tg_xxxyyzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 256); 

                auto tg_xxxyyzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 257); 

                auto tg_xxxyyzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 258); 

                auto tg_xxxyyzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 259); 

                auto tg_xxxyyzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 260); 

                auto tg_xxxyyzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 261); 

                auto tg_xxxyyzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 262); 

                auto tg_xxxyyzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 263); 

                auto tg_xxxyyzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 264); 

                auto tg_xxxyyzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 265); 

                auto tg_xxxyyzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 266); 

                auto tg_xxxyyzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 267); 

                auto tg_xxxyyzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 268); 

                auto tg_xxxyyzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 269); 

                auto tg_xxxyyzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 270); 

                auto tg_xxxyyzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 271); 

                auto tg_xxxyyzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 272); 

                auto tg_xxxyzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 273); 

                auto tg_xxxyzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 274); 

                auto tg_xxxyzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 275); 

                auto tg_xxxyzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 276); 

                auto tg_xxxyzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 277); 

                auto tg_xxxyzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 278); 

                auto tg_xxxyzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 279); 

                auto tg_xxxyzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 280); 

                auto tg_xxxyzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 281); 

                auto tg_xxxyzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 282); 

                auto tg_xxxyzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 283); 

                auto tg_xxxyzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 284); 

                // Batch of Integrals (190,285)

                #pragma omp simd aligned(fxn, fza, tg_xxxxzzz_xxxxy_0, tg_xxxxzzz_xxxxz_0, tg_xxxxzzz_xxxyy_0, \
                                         tg_xxxxzzz_xxxyz_0, tg_xxxxzzz_xxxzz_0, tg_xxxxzzz_xxyyy_0, tg_xxxxzzz_xxyyz_0, \
                                         tg_xxxxzzz_xxyzz_0, tg_xxxxzzz_xxzzz_0, tg_xxxxzzz_xyyyy_0, tg_xxxxzzz_xyyyz_0, \
                                         tg_xxxxzzz_xyyzz_0, tg_xxxxzzz_xyzzz_0, tg_xxxxzzz_xzzzz_0, tg_xxxxzzz_yyyyy_0, \
                                         tg_xxxxzzz_yyyyz_0, tg_xxxxzzz_yyyzz_0, tg_xxxxzzz_yyzzz_0, tg_xxxxzzz_yzzzz_0, \
                                         tg_xxxxzzz_zzzzz_0, tg_xxxyyyy_xxxxx_0, tg_xxxyyyy_xxxxy_0, tg_xxxyyyy_xxxxz_0, \
                                         tg_xxxyyyy_xxxyy_0, tg_xxxyyyy_xxxyz_0, tg_xxxyyyy_xxxzz_0, tg_xxxyyyy_xxyyy_0, \
                                         tg_xxxyyyy_xxyyz_0, tg_xxxyyyy_xxyzz_0, tg_xxxyyyy_xxzzz_0, tg_xxxyyyy_xyyyy_0, \
                                         tg_xxxyyyy_xyyyz_0, tg_xxxyyyy_xyyzz_0, tg_xxxyyyy_xyzzz_0, tg_xxxyyyy_xzzzz_0, \
                                         tg_xxxyyyy_yyyyy_0, tg_xxxyyyy_yyyyz_0, tg_xxxyyyy_yyyzz_0, tg_xxxyyyy_yyzzz_0, \
                                         tg_xxxyyyy_yzzzz_0, tg_xxxyyyy_zzzzz_0, tg_xxxyyyz_xxxxx_0, tg_xxxyyyz_xxxxy_0, \
                                         tg_xxxyyyz_xxxxz_0, tg_xxxyyyz_xxxyy_0, tg_xxxyyyz_xxxyz_0, tg_xxxyyyz_xxxzz_0, \
                                         tg_xxxyyyz_xxyyy_0, tg_xxxyyyz_xxyyz_0, tg_xxxyyyz_xxyzz_0, tg_xxxyyyz_xxzzz_0, \
                                         tg_xxxyyyz_xyyyy_0, tg_xxxyyyz_xyyyz_0, tg_xxxyyyz_xyyzz_0, tg_xxxyyyz_xyzzz_0, \
                                         tg_xxxyyyz_xzzzz_0, tg_xxxyyyz_yyyyy_0, tg_xxxyyyz_yyyyz_0, tg_xxxyyyz_yyyzz_0, \
                                         tg_xxxyyyz_yyzzz_0, tg_xxxyyyz_yzzzz_0, tg_xxxyyyz_zzzzz_0, tg_xxxyyzz_xxxxx_0, \
                                         tg_xxxyyzz_xxxxy_0, tg_xxxyyzz_xxxxz_0, tg_xxxyyzz_xxxyy_0, tg_xxxyyzz_xxxyz_0, \
                                         tg_xxxyyzz_xxxzz_0, tg_xxxyyzz_xxyyy_0, tg_xxxyyzz_xxyyz_0, tg_xxxyyzz_xxyzz_0, \
                                         tg_xxxyyzz_xxzzz_0, tg_xxxyyzz_xyyyy_0, tg_xxxyyzz_xyyyz_0, tg_xxxyyzz_xyyzz_0, \
                                         tg_xxxyyzz_xyzzz_0, tg_xxxyyzz_xzzzz_0, tg_xxxyyzz_yyyyy_0, tg_xxxyyzz_yyyyz_0, \
                                         tg_xxxyyzz_yyyzz_0, tg_xxxyyzz_yyzzz_0, tg_xxxyyzz_yzzzz_0, tg_xxxyyzz_zzzzz_0, \
                                         tg_xxxyzzz_xxxxx_0, tg_xxxyzzz_xxxxy_0, tg_xxxyzzz_xxxxz_0, tg_xxxyzzz_xxxyy_0, \
                                         tg_xxxyzzz_xxxyz_0, tg_xxxyzzz_xxxzz_0, tg_xxxyzzz_xxyyy_0, tg_xxxyzzz_xxyyz_0, \
                                         tg_xxxyzzz_xxyzz_0, tg_xxxyzzz_xxzzz_0, tg_xxxyzzz_xyyyy_0, tg_xxxyzzz_xyyyz_0, \
                                         tg_xxxzzz_xxxxy_0, tg_xxxzzz_xxxxy_1, tg_xxxzzz_xxxxz_0, tg_xxxzzz_xxxxz_1, \
                                         tg_xxxzzz_xxxy_1, tg_xxxzzz_xxxyy_0, tg_xxxzzz_xxxyy_1, tg_xxxzzz_xxxyz_0, \
                                         tg_xxxzzz_xxxyz_1, tg_xxxzzz_xxxz_1, tg_xxxzzz_xxxzz_0, tg_xxxzzz_xxxzz_1, \
                                         tg_xxxzzz_xxyy_1, tg_xxxzzz_xxyyy_0, tg_xxxzzz_xxyyy_1, tg_xxxzzz_xxyyz_0, \
                                         tg_xxxzzz_xxyyz_1, tg_xxxzzz_xxyz_1, tg_xxxzzz_xxyzz_0, tg_xxxzzz_xxyzz_1, \
                                         tg_xxxzzz_xxzz_1, tg_xxxzzz_xxzzz_0, tg_xxxzzz_xxzzz_1, tg_xxxzzz_xyyy_1, \
                                         tg_xxxzzz_xyyyy_0, tg_xxxzzz_xyyyy_1, tg_xxxzzz_xyyyz_0, tg_xxxzzz_xyyyz_1, \
                                         tg_xxxzzz_xyyz_1, tg_xxxzzz_xyyzz_0, tg_xxxzzz_xyyzz_1, tg_xxxzzz_xyzz_1, \
                                         tg_xxxzzz_xyzzz_0, tg_xxxzzz_xyzzz_1, tg_xxxzzz_xzzz_1, tg_xxxzzz_xzzzz_0, \
                                         tg_xxxzzz_xzzzz_1, tg_xxxzzz_yyyy_1, tg_xxxzzz_yyyyy_0, tg_xxxzzz_yyyyy_1, \
                                         tg_xxxzzz_yyyyz_0, tg_xxxzzz_yyyyz_1, tg_xxxzzz_yyyz_1, tg_xxxzzz_yyyzz_0, \
                                         tg_xxxzzz_yyyzz_1, tg_xxxzzz_yyzz_1, tg_xxxzzz_yyzzz_0, tg_xxxzzz_yyzzz_1, \
                                         tg_xxxzzz_yzzz_1, tg_xxxzzz_yzzzz_0, tg_xxxzzz_yzzzz_1, tg_xxxzzz_zzzz_1, \
                                         tg_xxxzzz_zzzzz_0, tg_xxxzzz_zzzzz_1, tg_xxyyyy_xxxx_1, tg_xxyyyy_xxxxx_0, \
                                         tg_xxyyyy_xxxxx_1, tg_xxyyyy_xxxxy_0, tg_xxyyyy_xxxxy_1, tg_xxyyyy_xxxxz_0, \
                                         tg_xxyyyy_xxxxz_1, tg_xxyyyy_xxxy_1, tg_xxyyyy_xxxyy_0, tg_xxyyyy_xxxyy_1, \
                                         tg_xxyyyy_xxxyz_0, tg_xxyyyy_xxxyz_1, tg_xxyyyy_xxxz_1, tg_xxyyyy_xxxzz_0, \
                                         tg_xxyyyy_xxxzz_1, tg_xxyyyy_xxyy_1, tg_xxyyyy_xxyyy_0, tg_xxyyyy_xxyyy_1, \
                                         tg_xxyyyy_xxyyz_0, tg_xxyyyy_xxyyz_1, tg_xxyyyy_xxyz_1, tg_xxyyyy_xxyzz_0, \
                                         tg_xxyyyy_xxyzz_1, tg_xxyyyy_xxzz_1, tg_xxyyyy_xxzzz_0, tg_xxyyyy_xxzzz_1, \
                                         tg_xxyyyy_xyyy_1, tg_xxyyyy_xyyyy_0, tg_xxyyyy_xyyyy_1, tg_xxyyyy_xyyyz_0, \
                                         tg_xxyyyy_xyyyz_1, tg_xxyyyy_xyyz_1, tg_xxyyyy_xyyzz_0, tg_xxyyyy_xyyzz_1, \
                                         tg_xxyyyy_xyzz_1, tg_xxyyyy_xyzzz_0, tg_xxyyyy_xyzzz_1, tg_xxyyyy_xzzz_1, \
                                         tg_xxyyyy_xzzzz_0, tg_xxyyyy_xzzzz_1, tg_xxyyyy_yyyy_1, tg_xxyyyy_yyyyy_0, \
                                         tg_xxyyyy_yyyyy_1, tg_xxyyyy_yyyyz_0, tg_xxyyyy_yyyyz_1, tg_xxyyyy_yyyz_1, \
                                         tg_xxyyyy_yyyzz_0, tg_xxyyyy_yyyzz_1, tg_xxyyyy_yyzz_1, tg_xxyyyy_yyzzz_0, \
                                         tg_xxyyyy_yyzzz_1, tg_xxyyyy_yzzz_1, tg_xxyyyy_yzzzz_0, tg_xxyyyy_yzzzz_1, \
                                         tg_xxyyyy_zzzz_1, tg_xxyyyy_zzzzz_0, tg_xxyyyy_zzzzz_1, tg_xxyyyz_xxxx_1, \
                                         tg_xxyyyz_xxxxx_0, tg_xxyyyz_xxxxx_1, tg_xxyyyz_xxxxy_0, tg_xxyyyz_xxxxy_1, \
                                         tg_xxyyyz_xxxxz_0, tg_xxyyyz_xxxxz_1, tg_xxyyyz_xxxy_1, tg_xxyyyz_xxxyy_0, \
                                         tg_xxyyyz_xxxyy_1, tg_xxyyyz_xxxyz_0, tg_xxyyyz_xxxyz_1, tg_xxyyyz_xxxz_1, \
                                         tg_xxyyyz_xxxzz_0, tg_xxyyyz_xxxzz_1, tg_xxyyyz_xxyy_1, tg_xxyyyz_xxyyy_0, \
                                         tg_xxyyyz_xxyyy_1, tg_xxyyyz_xxyyz_0, tg_xxyyyz_xxyyz_1, tg_xxyyyz_xxyz_1, \
                                         tg_xxyyyz_xxyzz_0, tg_xxyyyz_xxyzz_1, tg_xxyyyz_xxzz_1, tg_xxyyyz_xxzzz_0, \
                                         tg_xxyyyz_xxzzz_1, tg_xxyyyz_xyyy_1, tg_xxyyyz_xyyyy_0, tg_xxyyyz_xyyyy_1, \
                                         tg_xxyyyz_xyyyz_0, tg_xxyyyz_xyyyz_1, tg_xxyyyz_xyyz_1, tg_xxyyyz_xyyzz_0, \
                                         tg_xxyyyz_xyyzz_1, tg_xxyyyz_xyzz_1, tg_xxyyyz_xyzzz_0, tg_xxyyyz_xyzzz_1, \
                                         tg_xxyyyz_xzzz_1, tg_xxyyyz_xzzzz_0, tg_xxyyyz_xzzzz_1, tg_xxyyyz_yyyy_1, \
                                         tg_xxyyyz_yyyyy_0, tg_xxyyyz_yyyyy_1, tg_xxyyyz_yyyyz_0, tg_xxyyyz_yyyyz_1, \
                                         tg_xxyyyz_yyyz_1, tg_xxyyyz_yyyzz_0, tg_xxyyyz_yyyzz_1, tg_xxyyyz_yyzz_1, \
                                         tg_xxyyyz_yyzzz_0, tg_xxyyyz_yyzzz_1, tg_xxyyyz_yzzz_1, tg_xxyyyz_yzzzz_0, \
                                         tg_xxyyyz_yzzzz_1, tg_xxyyyz_zzzz_1, tg_xxyyyz_zzzzz_0, tg_xxyyyz_zzzzz_1, \
                                         tg_xxyyzz_xxxx_1, tg_xxyyzz_xxxxx_0, tg_xxyyzz_xxxxx_1, tg_xxyyzz_xxxxy_0, \
                                         tg_xxyyzz_xxxxy_1, tg_xxyyzz_xxxxz_0, tg_xxyyzz_xxxxz_1, tg_xxyyzz_xxxy_1, \
                                         tg_xxyyzz_xxxyy_0, tg_xxyyzz_xxxyy_1, tg_xxyyzz_xxxyz_0, tg_xxyyzz_xxxyz_1, \
                                         tg_xxyyzz_xxxz_1, tg_xxyyzz_xxxzz_0, tg_xxyyzz_xxxzz_1, tg_xxyyzz_xxyy_1, \
                                         tg_xxyyzz_xxyyy_0, tg_xxyyzz_xxyyy_1, tg_xxyyzz_xxyyz_0, tg_xxyyzz_xxyyz_1, \
                                         tg_xxyyzz_xxyz_1, tg_xxyyzz_xxyzz_0, tg_xxyyzz_xxyzz_1, tg_xxyyzz_xxzz_1, \
                                         tg_xxyyzz_xxzzz_0, tg_xxyyzz_xxzzz_1, tg_xxyyzz_xyyy_1, tg_xxyyzz_xyyyy_0, \
                                         tg_xxyyzz_xyyyy_1, tg_xxyyzz_xyyyz_0, tg_xxyyzz_xyyyz_1, tg_xxyyzz_xyyz_1, \
                                         tg_xxyyzz_xyyzz_0, tg_xxyyzz_xyyzz_1, tg_xxyyzz_xyzz_1, tg_xxyyzz_xyzzz_0, \
                                         tg_xxyyzz_xyzzz_1, tg_xxyyzz_xzzz_1, tg_xxyyzz_xzzzz_0, tg_xxyyzz_xzzzz_1, \
                                         tg_xxyyzz_yyyy_1, tg_xxyyzz_yyyyy_0, tg_xxyyzz_yyyyy_1, tg_xxyyzz_yyyyz_0, \
                                         tg_xxyyzz_yyyyz_1, tg_xxyyzz_yyyz_1, tg_xxyyzz_yyyzz_0, tg_xxyyzz_yyyzz_1, \
                                         tg_xxyyzz_yyzz_1, tg_xxyyzz_yyzzz_0, tg_xxyyzz_yyzzz_1, tg_xxyyzz_yzzz_1, \
                                         tg_xxyyzz_yzzzz_0, tg_xxyyzz_yzzzz_1, tg_xxyyzz_zzzz_1, tg_xxyyzz_zzzzz_0, \
                                         tg_xxyyzz_zzzzz_1, tg_xxyzzz_xxxx_1, tg_xxyzzz_xxxxx_0, tg_xxyzzz_xxxxx_1, \
                                         tg_xxyzzz_xxxxy_0, tg_xxyzzz_xxxxy_1, tg_xxyzzz_xxxxz_0, tg_xxyzzz_xxxxz_1, \
                                         tg_xxyzzz_xxxy_1, tg_xxyzzz_xxxyy_0, tg_xxyzzz_xxxyy_1, tg_xxyzzz_xxxyz_0, \
                                         tg_xxyzzz_xxxyz_1, tg_xxyzzz_xxxz_1, tg_xxyzzz_xxxzz_0, tg_xxyzzz_xxxzz_1, \
                                         tg_xxyzzz_xxyy_1, tg_xxyzzz_xxyyy_0, tg_xxyzzz_xxyyy_1, tg_xxyzzz_xxyyz_0, \
                                         tg_xxyzzz_xxyyz_1, tg_xxyzzz_xxyz_1, tg_xxyzzz_xxyzz_0, tg_xxyzzz_xxyzz_1, \
                                         tg_xxyzzz_xxzz_1, tg_xxyzzz_xxzzz_0, tg_xxyzzz_xxzzz_1, tg_xxyzzz_xyyy_1, \
                                         tg_xxyzzz_xyyyy_0, tg_xxyzzz_xyyyy_1, tg_xxyzzz_xyyyz_0, tg_xxyzzz_xyyyz_1, \
                                         tg_xxyzzz_xyyz_1, tg_xxyzzz_xyzz_1, tg_xxyzzz_xzzz_1, tg_xxyzzz_yyyy_1, \
                                         tg_xxyzzz_yyyz_1, tg_xxzzz_xxxxy_0, tg_xxzzz_xxxxy_1, tg_xxzzz_xxxxz_0, \
                                         tg_xxzzz_xxxxz_1, tg_xxzzz_xxxyy_0, tg_xxzzz_xxxyy_1, tg_xxzzz_xxxyz_0, \
                                         tg_xxzzz_xxxyz_1, tg_xxzzz_xxxzz_0, tg_xxzzz_xxxzz_1, tg_xxzzz_xxyyy_0, \
                                         tg_xxzzz_xxyyy_1, tg_xxzzz_xxyyz_0, tg_xxzzz_xxyyz_1, tg_xxzzz_xxyzz_0, \
                                         tg_xxzzz_xxyzz_1, tg_xxzzz_xxzzz_0, tg_xxzzz_xxzzz_1, tg_xxzzz_xyyyy_0, \
                                         tg_xxzzz_xyyyy_1, tg_xxzzz_xyyyz_0, tg_xxzzz_xyyyz_1, tg_xxzzz_xyyzz_0, \
                                         tg_xxzzz_xyyzz_1, tg_xxzzz_xyzzz_0, tg_xxzzz_xyzzz_1, tg_xxzzz_xzzzz_0, \
                                         tg_xxzzz_xzzzz_1, tg_xxzzz_yyyyy_0, tg_xxzzz_yyyyy_1, tg_xxzzz_yyyyz_0, \
                                         tg_xxzzz_yyyyz_1, tg_xxzzz_yyyzz_0, tg_xxzzz_yyyzz_1, tg_xxzzz_yyzzz_0, \
                                         tg_xxzzz_yyzzz_1, tg_xxzzz_yzzzz_0, tg_xxzzz_yzzzz_1, tg_xxzzz_zzzzz_0, \
                                         tg_xxzzz_zzzzz_1, tg_xyyyy_xxxxx_0, tg_xyyyy_xxxxx_1, tg_xyyyy_xxxxy_0, \
                                         tg_xyyyy_xxxxy_1, tg_xyyyy_xxxxz_0, tg_xyyyy_xxxxz_1, tg_xyyyy_xxxyy_0, \
                                         tg_xyyyy_xxxyy_1, tg_xyyyy_xxxyz_0, tg_xyyyy_xxxyz_1, tg_xyyyy_xxxzz_0, \
                                         tg_xyyyy_xxxzz_1, tg_xyyyy_xxyyy_0, tg_xyyyy_xxyyy_1, tg_xyyyy_xxyyz_0, \
                                         tg_xyyyy_xxyyz_1, tg_xyyyy_xxyzz_0, tg_xyyyy_xxyzz_1, tg_xyyyy_xxzzz_0, \
                                         tg_xyyyy_xxzzz_1, tg_xyyyy_xyyyy_0, tg_xyyyy_xyyyy_1, tg_xyyyy_xyyyz_0, \
                                         tg_xyyyy_xyyyz_1, tg_xyyyy_xyyzz_0, tg_xyyyy_xyyzz_1, tg_xyyyy_xyzzz_0, \
                                         tg_xyyyy_xyzzz_1, tg_xyyyy_xzzzz_0, tg_xyyyy_xzzzz_1, tg_xyyyy_yyyyy_0, \
                                         tg_xyyyy_yyyyy_1, tg_xyyyy_yyyyz_0, tg_xyyyy_yyyyz_1, tg_xyyyy_yyyzz_0, \
                                         tg_xyyyy_yyyzz_1, tg_xyyyy_yyzzz_0, tg_xyyyy_yyzzz_1, tg_xyyyy_yzzzz_0, \
                                         tg_xyyyy_yzzzz_1, tg_xyyyy_zzzzz_0, tg_xyyyy_zzzzz_1, tg_xyyyz_xxxxx_0, \
                                         tg_xyyyz_xxxxx_1, tg_xyyyz_xxxxy_0, tg_xyyyz_xxxxy_1, tg_xyyyz_xxxxz_0, \
                                         tg_xyyyz_xxxxz_1, tg_xyyyz_xxxyy_0, tg_xyyyz_xxxyy_1, tg_xyyyz_xxxyz_0, \
                                         tg_xyyyz_xxxyz_1, tg_xyyyz_xxxzz_0, tg_xyyyz_xxxzz_1, tg_xyyyz_xxyyy_0, \
                                         tg_xyyyz_xxyyy_1, tg_xyyyz_xxyyz_0, tg_xyyyz_xxyyz_1, tg_xyyyz_xxyzz_0, \
                                         tg_xyyyz_xxyzz_1, tg_xyyyz_xxzzz_0, tg_xyyyz_xxzzz_1, tg_xyyyz_xyyyy_0, \
                                         tg_xyyyz_xyyyy_1, tg_xyyyz_xyyyz_0, tg_xyyyz_xyyyz_1, tg_xyyyz_xyyzz_0, \
                                         tg_xyyyz_xyyzz_1, tg_xyyyz_xyzzz_0, tg_xyyyz_xyzzz_1, tg_xyyyz_xzzzz_0, \
                                         tg_xyyyz_xzzzz_1, tg_xyyyz_yyyyy_0, tg_xyyyz_yyyyy_1, tg_xyyyz_yyyyz_0, \
                                         tg_xyyyz_yyyyz_1, tg_xyyyz_yyyzz_0, tg_xyyyz_yyyzz_1, tg_xyyyz_yyzzz_0, \
                                         tg_xyyyz_yyzzz_1, tg_xyyyz_yzzzz_0, tg_xyyyz_yzzzz_1, tg_xyyyz_zzzzz_0, \
                                         tg_xyyyz_zzzzz_1, tg_xyyzz_xxxxx_0, tg_xyyzz_xxxxx_1, tg_xyyzz_xxxxy_0, \
                                         tg_xyyzz_xxxxy_1, tg_xyyzz_xxxxz_0, tg_xyyzz_xxxxz_1, tg_xyyzz_xxxyy_0, \
                                         tg_xyyzz_xxxyy_1, tg_xyyzz_xxxyz_0, tg_xyyzz_xxxyz_1, tg_xyyzz_xxxzz_0, \
                                         tg_xyyzz_xxxzz_1, tg_xyyzz_xxyyy_0, tg_xyyzz_xxyyy_1, tg_xyyzz_xxyyz_0, \
                                         tg_xyyzz_xxyyz_1, tg_xyyzz_xxyzz_0, tg_xyyzz_xxyzz_1, tg_xyyzz_xxzzz_0, \
                                         tg_xyyzz_xxzzz_1, tg_xyyzz_xyyyy_0, tg_xyyzz_xyyyy_1, tg_xyyzz_xyyyz_0, \
                                         tg_xyyzz_xyyyz_1, tg_xyyzz_xyyzz_0, tg_xyyzz_xyyzz_1, tg_xyyzz_xyzzz_0, \
                                         tg_xyyzz_xyzzz_1, tg_xyyzz_xzzzz_0, tg_xyyzz_xzzzz_1, tg_xyyzz_yyyyy_0, \
                                         tg_xyyzz_yyyyy_1, tg_xyyzz_yyyyz_0, tg_xyyzz_yyyyz_1, tg_xyyzz_yyyzz_0, \
                                         tg_xyyzz_yyyzz_1, tg_xyyzz_yyzzz_0, tg_xyyzz_yyzzz_1, tg_xyyzz_yzzzz_0, \
                                         tg_xyyzz_yzzzz_1, tg_xyyzz_zzzzz_0, tg_xyyzz_zzzzz_1, tg_xyzzz_xxxxx_0, \
                                         tg_xyzzz_xxxxx_1, tg_xyzzz_xxxxy_0, tg_xyzzz_xxxxy_1, tg_xyzzz_xxxxz_0, \
                                         tg_xyzzz_xxxxz_1, tg_xyzzz_xxxyy_0, tg_xyzzz_xxxyy_1, tg_xyzzz_xxxyz_0, \
                                         tg_xyzzz_xxxyz_1, tg_xyzzz_xxxzz_0, tg_xyzzz_xxxzz_1, tg_xyzzz_xxyyy_0, \
                                         tg_xyzzz_xxyyy_1, tg_xyzzz_xxyyz_0, tg_xyzzz_xxyyz_1, tg_xyzzz_xxyzz_0, \
                                         tg_xyzzz_xxyzz_1, tg_xyzzz_xxzzz_0, tg_xyzzz_xxzzz_1, tg_xyzzz_xyyyy_0, \
                                         tg_xyzzz_xyyyy_1, tg_xyzzz_xyyyz_0, tg_xyzzz_xyyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxzzz_xxxxy_0[j] = pb_x * tg_xxxzzz_xxxxy_0[j] + fr * tg_xxxzzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxy_0[j] - tg_xxzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzz_xxxy_1[j];

                    tg_xxxxzzz_xxxxz_0[j] = pb_x * tg_xxxzzz_xxxxz_0[j] + fr * tg_xxxzzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxz_0[j] - tg_xxzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzz_xxxz_1[j];

                    tg_xxxxzzz_xxxyy_0[j] = pb_x * tg_xxxzzz_xxxyy_0[j] + fr * tg_xxxzzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxyy_0[j] - tg_xxzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxyy_1[j];

                    tg_xxxxzzz_xxxyz_0[j] = pb_x * tg_xxxzzz_xxxyz_0[j] + fr * tg_xxxzzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxyz_0[j] - tg_xxzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxyz_1[j];

                    tg_xxxxzzz_xxxzz_0[j] = pb_x * tg_xxxzzz_xxxzz_0[j] + fr * tg_xxxzzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxzz_0[j] - tg_xxzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxzz_1[j];

                    tg_xxxxzzz_xxyyy_0[j] = pb_x * tg_xxxzzz_xxyyy_0[j] + fr * tg_xxxzzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyyy_0[j] - tg_xxzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyyy_1[j];

                    tg_xxxxzzz_xxyyz_0[j] = pb_x * tg_xxxzzz_xxyyz_0[j] + fr * tg_xxxzzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyyz_0[j] - tg_xxzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyyz_1[j];

                    tg_xxxxzzz_xxyzz_0[j] = pb_x * tg_xxxzzz_xxyzz_0[j] + fr * tg_xxxzzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyzz_0[j] - tg_xxzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyzz_1[j];

                    tg_xxxxzzz_xxzzz_0[j] = pb_x * tg_xxxzzz_xxzzz_0[j] + fr * tg_xxxzzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxzzz_0[j] - tg_xxzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xzzz_1[j];

                    tg_xxxxzzz_xyyyy_0[j] = pb_x * tg_xxxzzz_xyyyy_0[j] + fr * tg_xxxzzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyyy_0[j] - tg_xxzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyyy_1[j];

                    tg_xxxxzzz_xyyyz_0[j] = pb_x * tg_xxxzzz_xyyyz_0[j] + fr * tg_xxxzzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyyz_0[j] - tg_xxzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyyz_1[j];

                    tg_xxxxzzz_xyyzz_0[j] = pb_x * tg_xxxzzz_xyyzz_0[j] + fr * tg_xxxzzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyzz_0[j] - tg_xxzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyzz_1[j];

                    tg_xxxxzzz_xyzzz_0[j] = pb_x * tg_xxxzzz_xyzzz_0[j] + fr * tg_xxxzzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyzzz_0[j] - tg_xxzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yzzz_1[j];

                    tg_xxxxzzz_xzzzz_0[j] = pb_x * tg_xxxzzz_xzzzz_0[j] + fr * tg_xxxzzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xzzzz_0[j] - tg_xxzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_zzzz_1[j];

                    tg_xxxxzzz_yyyyy_0[j] = pb_x * tg_xxxzzz_yyyyy_0[j] + fr * tg_xxxzzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyyy_0[j] - tg_xxzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxzzz_yyyyz_0[j] = pb_x * tg_xxxzzz_yyyyz_0[j] + fr * tg_xxxzzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyyz_0[j] - tg_xxzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxzzz_yyyzz_0[j] = pb_x * tg_xxxzzz_yyyzz_0[j] + fr * tg_xxxzzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyzz_0[j] - tg_xxzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxzzz_yyzzz_0[j] = pb_x * tg_xxxzzz_yyzzz_0[j] + fr * tg_xxxzzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyzzz_0[j] - tg_xxzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_yzzzz_0[j] = pb_x * tg_xxxzzz_yzzzz_0[j] + fr * tg_xxxzzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yzzzz_0[j] - tg_xxzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_zzzzz_0[j] = pb_x * tg_xxxzzz_zzzzz_0[j] + fr * tg_xxxzzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_zzzzz_0[j] - tg_xxzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_xxxxx_0[j] = pb_x * tg_xxyyyy_xxxxx_0[j] + fr * tg_xxyyyy_xxxxx_1[j] + fl1_fx * (tg_xyyyy_xxxxx_0[j] - tg_xyyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyy_xxxx_1[j];

                    tg_xxxyyyy_xxxxy_0[j] = pb_x * tg_xxyyyy_xxxxy_0[j] + fr * tg_xxyyyy_xxxxy_1[j] + fl1_fx * (tg_xyyyy_xxxxy_0[j] - tg_xyyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyy_xxxy_1[j];

                    tg_xxxyyyy_xxxxz_0[j] = pb_x * tg_xxyyyy_xxxxz_0[j] + fr * tg_xxyyyy_xxxxz_1[j] + fl1_fx * (tg_xyyyy_xxxxz_0[j] - tg_xyyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyy_xxxz_1[j];

                    tg_xxxyyyy_xxxyy_0[j] = pb_x * tg_xxyyyy_xxxyy_0[j] + fr * tg_xxyyyy_xxxyy_1[j] + fl1_fx * (tg_xyyyy_xxxyy_0[j] - tg_xyyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxyy_1[j];

                    tg_xxxyyyy_xxxyz_0[j] = pb_x * tg_xxyyyy_xxxyz_0[j] + fr * tg_xxyyyy_xxxyz_1[j] + fl1_fx * (tg_xyyyy_xxxyz_0[j] - tg_xyyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxyz_1[j];

                    tg_xxxyyyy_xxxzz_0[j] = pb_x * tg_xxyyyy_xxxzz_0[j] + fr * tg_xxyyyy_xxxzz_1[j] + fl1_fx * (tg_xyyyy_xxxzz_0[j] - tg_xyyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxzz_1[j];

                    tg_xxxyyyy_xxyyy_0[j] = pb_x * tg_xxyyyy_xxyyy_0[j] + fr * tg_xxyyyy_xxyyy_1[j] + fl1_fx * (tg_xyyyy_xxyyy_0[j] - tg_xyyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyyy_1[j];

                    tg_xxxyyyy_xxyyz_0[j] = pb_x * tg_xxyyyy_xxyyz_0[j] + fr * tg_xxyyyy_xxyyz_1[j] + fl1_fx * (tg_xyyyy_xxyyz_0[j] - tg_xyyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyyz_1[j];

                    tg_xxxyyyy_xxyzz_0[j] = pb_x * tg_xxyyyy_xxyzz_0[j] + fr * tg_xxyyyy_xxyzz_1[j] + fl1_fx * (tg_xyyyy_xxyzz_0[j] - tg_xyyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyzz_1[j];

                    tg_xxxyyyy_xxzzz_0[j] = pb_x * tg_xxyyyy_xxzzz_0[j] + fr * tg_xxyyyy_xxzzz_1[j] + fl1_fx * (tg_xyyyy_xxzzz_0[j] - tg_xyyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xzzz_1[j];

                    tg_xxxyyyy_xyyyy_0[j] = pb_x * tg_xxyyyy_xyyyy_0[j] + fr * tg_xxyyyy_xyyyy_1[j] + fl1_fx * (tg_xyyyy_xyyyy_0[j] - tg_xyyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyyy_1[j];

                    tg_xxxyyyy_xyyyz_0[j] = pb_x * tg_xxyyyy_xyyyz_0[j] + fr * tg_xxyyyy_xyyyz_1[j] + fl1_fx * (tg_xyyyy_xyyyz_0[j] - tg_xyyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyyz_1[j];

                    tg_xxxyyyy_xyyzz_0[j] = pb_x * tg_xxyyyy_xyyzz_0[j] + fr * tg_xxyyyy_xyyzz_1[j] + fl1_fx * (tg_xyyyy_xyyzz_0[j] - tg_xyyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyzz_1[j];

                    tg_xxxyyyy_xyzzz_0[j] = pb_x * tg_xxyyyy_xyzzz_0[j] + fr * tg_xxyyyy_xyzzz_1[j] + fl1_fx * (tg_xyyyy_xyzzz_0[j] - tg_xyyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yzzz_1[j];

                    tg_xxxyyyy_xzzzz_0[j] = pb_x * tg_xxyyyy_xzzzz_0[j] + fr * tg_xxyyyy_xzzzz_1[j] + fl1_fx * (tg_xyyyy_xzzzz_0[j] - tg_xyyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_zzzz_1[j];

                    tg_xxxyyyy_yyyyy_0[j] = pb_x * tg_xxyyyy_yyyyy_0[j] + fr * tg_xxyyyy_yyyyy_1[j] + fl1_fx * (tg_xyyyy_yyyyy_0[j] - tg_xyyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyyy_yyyyz_0[j] = pb_x * tg_xxyyyy_yyyyz_0[j] + fr * tg_xxyyyy_yyyyz_1[j] + fl1_fx * (tg_xyyyy_yyyyz_0[j] - tg_xyyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyyy_yyyzz_0[j] = pb_x * tg_xxyyyy_yyyzz_0[j] + fr * tg_xxyyyy_yyyzz_1[j] + fl1_fx * (tg_xyyyy_yyyzz_0[j] - tg_xyyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyyy_yyzzz_0[j] = pb_x * tg_xxyyyy_yyzzz_0[j] + fr * tg_xxyyyy_yyzzz_1[j] + fl1_fx * (tg_xyyyy_yyzzz_0[j] - tg_xyyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_yzzzz_0[j] = pb_x * tg_xxyyyy_yzzzz_0[j] + fr * tg_xxyyyy_yzzzz_1[j] + fl1_fx * (tg_xyyyy_yzzzz_0[j] - tg_xyyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_zzzzz_0[j] = pb_x * tg_xxyyyy_zzzzz_0[j] + fr * tg_xxyyyy_zzzzz_1[j] + fl1_fx * (tg_xyyyy_zzzzz_0[j] - tg_xyyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_xxxxx_0[j] = pb_x * tg_xxyyyz_xxxxx_0[j] + fr * tg_xxyyyz_xxxxx_1[j] + fl1_fx * (tg_xyyyz_xxxxx_0[j] - tg_xyyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyz_xxxx_1[j];

                    tg_xxxyyyz_xxxxy_0[j] = pb_x * tg_xxyyyz_xxxxy_0[j] + fr * tg_xxyyyz_xxxxy_1[j] + fl1_fx * (tg_xyyyz_xxxxy_0[j] - tg_xyyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyz_xxxy_1[j];

                    tg_xxxyyyz_xxxxz_0[j] = pb_x * tg_xxyyyz_xxxxz_0[j] + fr * tg_xxyyyz_xxxxz_1[j] + fl1_fx * (tg_xyyyz_xxxxz_0[j] - tg_xyyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyz_xxxz_1[j];

                    tg_xxxyyyz_xxxyy_0[j] = pb_x * tg_xxyyyz_xxxyy_0[j] + fr * tg_xxyyyz_xxxyy_1[j] + fl1_fx * (tg_xyyyz_xxxyy_0[j] - tg_xyyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxyy_1[j];

                    tg_xxxyyyz_xxxyz_0[j] = pb_x * tg_xxyyyz_xxxyz_0[j] + fr * tg_xxyyyz_xxxyz_1[j] + fl1_fx * (tg_xyyyz_xxxyz_0[j] - tg_xyyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxyz_1[j];

                    tg_xxxyyyz_xxxzz_0[j] = pb_x * tg_xxyyyz_xxxzz_0[j] + fr * tg_xxyyyz_xxxzz_1[j] + fl1_fx * (tg_xyyyz_xxxzz_0[j] - tg_xyyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxzz_1[j];

                    tg_xxxyyyz_xxyyy_0[j] = pb_x * tg_xxyyyz_xxyyy_0[j] + fr * tg_xxyyyz_xxyyy_1[j] + fl1_fx * (tg_xyyyz_xxyyy_0[j] - tg_xyyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyyy_1[j];

                    tg_xxxyyyz_xxyyz_0[j] = pb_x * tg_xxyyyz_xxyyz_0[j] + fr * tg_xxyyyz_xxyyz_1[j] + fl1_fx * (tg_xyyyz_xxyyz_0[j] - tg_xyyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyyz_1[j];

                    tg_xxxyyyz_xxyzz_0[j] = pb_x * tg_xxyyyz_xxyzz_0[j] + fr * tg_xxyyyz_xxyzz_1[j] + fl1_fx * (tg_xyyyz_xxyzz_0[j] - tg_xyyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyzz_1[j];

                    tg_xxxyyyz_xxzzz_0[j] = pb_x * tg_xxyyyz_xxzzz_0[j] + fr * tg_xxyyyz_xxzzz_1[j] + fl1_fx * (tg_xyyyz_xxzzz_0[j] - tg_xyyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xzzz_1[j];

                    tg_xxxyyyz_xyyyy_0[j] = pb_x * tg_xxyyyz_xyyyy_0[j] + fr * tg_xxyyyz_xyyyy_1[j] + fl1_fx * (tg_xyyyz_xyyyy_0[j] - tg_xyyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyyy_1[j];

                    tg_xxxyyyz_xyyyz_0[j] = pb_x * tg_xxyyyz_xyyyz_0[j] + fr * tg_xxyyyz_xyyyz_1[j] + fl1_fx * (tg_xyyyz_xyyyz_0[j] - tg_xyyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyyz_1[j];

                    tg_xxxyyyz_xyyzz_0[j] = pb_x * tg_xxyyyz_xyyzz_0[j] + fr * tg_xxyyyz_xyyzz_1[j] + fl1_fx * (tg_xyyyz_xyyzz_0[j] - tg_xyyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyzz_1[j];

                    tg_xxxyyyz_xyzzz_0[j] = pb_x * tg_xxyyyz_xyzzz_0[j] + fr * tg_xxyyyz_xyzzz_1[j] + fl1_fx * (tg_xyyyz_xyzzz_0[j] - tg_xyyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yzzz_1[j];

                    tg_xxxyyyz_xzzzz_0[j] = pb_x * tg_xxyyyz_xzzzz_0[j] + fr * tg_xxyyyz_xzzzz_1[j] + fl1_fx * (tg_xyyyz_xzzzz_0[j] - tg_xyyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_zzzz_1[j];

                    tg_xxxyyyz_yyyyy_0[j] = pb_x * tg_xxyyyz_yyyyy_0[j] + fr * tg_xxyyyz_yyyyy_1[j] + fl1_fx * (tg_xyyyz_yyyyy_0[j] - tg_xyyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyyz_yyyyz_0[j] = pb_x * tg_xxyyyz_yyyyz_0[j] + fr * tg_xxyyyz_yyyyz_1[j] + fl1_fx * (tg_xyyyz_yyyyz_0[j] - tg_xyyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyyz_yyyzz_0[j] = pb_x * tg_xxyyyz_yyyzz_0[j] + fr * tg_xxyyyz_yyyzz_1[j] + fl1_fx * (tg_xyyyz_yyyzz_0[j] - tg_xyyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyyz_yyzzz_0[j] = pb_x * tg_xxyyyz_yyzzz_0[j] + fr * tg_xxyyyz_yyzzz_1[j] + fl1_fx * (tg_xyyyz_yyzzz_0[j] - tg_xyyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_yzzzz_0[j] = pb_x * tg_xxyyyz_yzzzz_0[j] + fr * tg_xxyyyz_yzzzz_1[j] + fl1_fx * (tg_xyyyz_yzzzz_0[j] - tg_xyyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_zzzzz_0[j] = pb_x * tg_xxyyyz_zzzzz_0[j] + fr * tg_xxyyyz_zzzzz_1[j] + fl1_fx * (tg_xyyyz_zzzzz_0[j] - tg_xyyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_xxxxx_0[j] = pb_x * tg_xxyyzz_xxxxx_0[j] + fr * tg_xxyyzz_xxxxx_1[j] + fl1_fx * (tg_xyyzz_xxxxx_0[j] - tg_xyyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyzz_xxxx_1[j];

                    tg_xxxyyzz_xxxxy_0[j] = pb_x * tg_xxyyzz_xxxxy_0[j] + fr * tg_xxyyzz_xxxxy_1[j] + fl1_fx * (tg_xyyzz_xxxxy_0[j] - tg_xyyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzz_xxxy_1[j];

                    tg_xxxyyzz_xxxxz_0[j] = pb_x * tg_xxyyzz_xxxxz_0[j] + fr * tg_xxyyzz_xxxxz_1[j] + fl1_fx * (tg_xyyzz_xxxxz_0[j] - tg_xyyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzz_xxxz_1[j];

                    tg_xxxyyzz_xxxyy_0[j] = pb_x * tg_xxyyzz_xxxyy_0[j] + fr * tg_xxyyzz_xxxyy_1[j] + fl1_fx * (tg_xyyzz_xxxyy_0[j] - tg_xyyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxyy_1[j];

                    tg_xxxyyzz_xxxyz_0[j] = pb_x * tg_xxyyzz_xxxyz_0[j] + fr * tg_xxyyzz_xxxyz_1[j] + fl1_fx * (tg_xyyzz_xxxyz_0[j] - tg_xyyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxyz_1[j];

                    tg_xxxyyzz_xxxzz_0[j] = pb_x * tg_xxyyzz_xxxzz_0[j] + fr * tg_xxyyzz_xxxzz_1[j] + fl1_fx * (tg_xyyzz_xxxzz_0[j] - tg_xyyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxzz_1[j];

                    tg_xxxyyzz_xxyyy_0[j] = pb_x * tg_xxyyzz_xxyyy_0[j] + fr * tg_xxyyzz_xxyyy_1[j] + fl1_fx * (tg_xyyzz_xxyyy_0[j] - tg_xyyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyyy_1[j];

                    tg_xxxyyzz_xxyyz_0[j] = pb_x * tg_xxyyzz_xxyyz_0[j] + fr * tg_xxyyzz_xxyyz_1[j] + fl1_fx * (tg_xyyzz_xxyyz_0[j] - tg_xyyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyyz_1[j];

                    tg_xxxyyzz_xxyzz_0[j] = pb_x * tg_xxyyzz_xxyzz_0[j] + fr * tg_xxyyzz_xxyzz_1[j] + fl1_fx * (tg_xyyzz_xxyzz_0[j] - tg_xyyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyzz_1[j];

                    tg_xxxyyzz_xxzzz_0[j] = pb_x * tg_xxyyzz_xxzzz_0[j] + fr * tg_xxyyzz_xxzzz_1[j] + fl1_fx * (tg_xyyzz_xxzzz_0[j] - tg_xyyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xzzz_1[j];

                    tg_xxxyyzz_xyyyy_0[j] = pb_x * tg_xxyyzz_xyyyy_0[j] + fr * tg_xxyyzz_xyyyy_1[j] + fl1_fx * (tg_xyyzz_xyyyy_0[j] - tg_xyyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyyy_1[j];

                    tg_xxxyyzz_xyyyz_0[j] = pb_x * tg_xxyyzz_xyyyz_0[j] + fr * tg_xxyyzz_xyyyz_1[j] + fl1_fx * (tg_xyyzz_xyyyz_0[j] - tg_xyyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyyz_1[j];

                    tg_xxxyyzz_xyyzz_0[j] = pb_x * tg_xxyyzz_xyyzz_0[j] + fr * tg_xxyyzz_xyyzz_1[j] + fl1_fx * (tg_xyyzz_xyyzz_0[j] - tg_xyyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyzz_1[j];

                    tg_xxxyyzz_xyzzz_0[j] = pb_x * tg_xxyyzz_xyzzz_0[j] + fr * tg_xxyyzz_xyzzz_1[j] + fl1_fx * (tg_xyyzz_xyzzz_0[j] - tg_xyyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yzzz_1[j];

                    tg_xxxyyzz_xzzzz_0[j] = pb_x * tg_xxyyzz_xzzzz_0[j] + fr * tg_xxyyzz_xzzzz_1[j] + fl1_fx * (tg_xyyzz_xzzzz_0[j] - tg_xyyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_zzzz_1[j];

                    tg_xxxyyzz_yyyyy_0[j] = pb_x * tg_xxyyzz_yyyyy_0[j] + fr * tg_xxyyzz_yyyyy_1[j] + fl1_fx * (tg_xyyzz_yyyyy_0[j] - tg_xyyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyzz_yyyyz_0[j] = pb_x * tg_xxyyzz_yyyyz_0[j] + fr * tg_xxyyzz_yyyyz_1[j] + fl1_fx * (tg_xyyzz_yyyyz_0[j] - tg_xyyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyzz_yyyzz_0[j] = pb_x * tg_xxyyzz_yyyzz_0[j] + fr * tg_xxyyzz_yyyzz_1[j] + fl1_fx * (tg_xyyzz_yyyzz_0[j] - tg_xyyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyzz_yyzzz_0[j] = pb_x * tg_xxyyzz_yyzzz_0[j] + fr * tg_xxyyzz_yyzzz_1[j] + fl1_fx * (tg_xyyzz_yyzzz_0[j] - tg_xyyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_yzzzz_0[j] = pb_x * tg_xxyyzz_yzzzz_0[j] + fr * tg_xxyyzz_yzzzz_1[j] + fl1_fx * (tg_xyyzz_yzzzz_0[j] - tg_xyyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_zzzzz_0[j] = pb_x * tg_xxyyzz_zzzzz_0[j] + fr * tg_xxyyzz_zzzzz_1[j] + fl1_fx * (tg_xyyzz_zzzzz_0[j] - tg_xyyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_xxxxx_0[j] = pb_x * tg_xxyzzz_xxxxx_0[j] + fr * tg_xxyzzz_xxxxx_1[j] + fl1_fx * (tg_xyzzz_xxxxx_0[j] - tg_xyzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyzzz_xxxx_1[j];

                    tg_xxxyzzz_xxxxy_0[j] = pb_x * tg_xxyzzz_xxxxy_0[j] + fr * tg_xxyzzz_xxxxy_1[j] + fl1_fx * (tg_xyzzz_xxxxy_0[j] - tg_xyzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzz_xxxy_1[j];

                    tg_xxxyzzz_xxxxz_0[j] = pb_x * tg_xxyzzz_xxxxz_0[j] + fr * tg_xxyzzz_xxxxz_1[j] + fl1_fx * (tg_xyzzz_xxxxz_0[j] - tg_xyzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzz_xxxz_1[j];

                    tg_xxxyzzz_xxxyy_0[j] = pb_x * tg_xxyzzz_xxxyy_0[j] + fr * tg_xxyzzz_xxxyy_1[j] + fl1_fx * (tg_xyzzz_xxxyy_0[j] - tg_xyzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxyy_1[j];

                    tg_xxxyzzz_xxxyz_0[j] = pb_x * tg_xxyzzz_xxxyz_0[j] + fr * tg_xxyzzz_xxxyz_1[j] + fl1_fx * (tg_xyzzz_xxxyz_0[j] - tg_xyzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxyz_1[j];

                    tg_xxxyzzz_xxxzz_0[j] = pb_x * tg_xxyzzz_xxxzz_0[j] + fr * tg_xxyzzz_xxxzz_1[j] + fl1_fx * (tg_xyzzz_xxxzz_0[j] - tg_xyzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxzz_1[j];

                    tg_xxxyzzz_xxyyy_0[j] = pb_x * tg_xxyzzz_xxyyy_0[j] + fr * tg_xxyzzz_xxyyy_1[j] + fl1_fx * (tg_xyzzz_xxyyy_0[j] - tg_xyzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyyy_1[j];

                    tg_xxxyzzz_xxyyz_0[j] = pb_x * tg_xxyzzz_xxyyz_0[j] + fr * tg_xxyzzz_xxyyz_1[j] + fl1_fx * (tg_xyzzz_xxyyz_0[j] - tg_xyzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyyz_1[j];

                    tg_xxxyzzz_xxyzz_0[j] = pb_x * tg_xxyzzz_xxyzz_0[j] + fr * tg_xxyzzz_xxyzz_1[j] + fl1_fx * (tg_xyzzz_xxyzz_0[j] - tg_xyzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyzz_1[j];

                    tg_xxxyzzz_xxzzz_0[j] = pb_x * tg_xxyzzz_xxzzz_0[j] + fr * tg_xxyzzz_xxzzz_1[j] + fl1_fx * (tg_xyzzz_xxzzz_0[j] - tg_xyzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xzzz_1[j];

                    tg_xxxyzzz_xyyyy_0[j] = pb_x * tg_xxyzzz_xyyyy_0[j] + fr * tg_xxyzzz_xyyyy_1[j] + fl1_fx * (tg_xyzzz_xyyyy_0[j] - tg_xyzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyyy_1[j];

                    tg_xxxyzzz_xyyyz_0[j] = pb_x * tg_xxyzzz_xyyyz_0[j] + fr * tg_xxyzzz_xyyyz_1[j] + fl1_fx * (tg_xyzzz_xyyyz_0[j] - tg_xyzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSH_285_380(      CMemBlock2D<double>* primBuffer,
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
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_xxyzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 285); 

                auto tg_xxyzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 286); 

                auto tg_xxyzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 287); 

                auto tg_xxyzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 288); 

                auto tg_xxyzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 289); 

                auto tg_xxyzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 290); 

                auto tg_xxyzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 291); 

                auto tg_xxyzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 292); 

                auto tg_xxyzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 293); 

                auto tg_xxzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 294); 

                auto tg_xxzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 295); 

                auto tg_xxzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 296); 

                auto tg_xxzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 297); 

                auto tg_xxzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 298); 

                auto tg_xxzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 299); 

                auto tg_xxzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 300); 

                auto tg_xxzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 301); 

                auto tg_xxzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 302); 

                auto tg_xxzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 303); 

                auto tg_xxzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 304); 

                auto tg_xxzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 305); 

                auto tg_xxzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 306); 

                auto tg_xxzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 307); 

                auto tg_xxzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 308); 

                auto tg_xxzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 309); 

                auto tg_xxzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 310); 

                auto tg_xxzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 311); 

                auto tg_xxzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 312); 

                auto tg_xxzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 313); 

                auto tg_xxzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 314); 

                auto tg_xyyyyy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 315); 

                auto tg_xyyyyy_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 316); 

                auto tg_xyyyyy_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 317); 

                auto tg_xyyyyy_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 318); 

                auto tg_xyyyyy_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 319); 

                auto tg_xyyyyy_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 320); 

                auto tg_xyyyyy_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 321); 

                auto tg_xyyyyy_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 322); 

                auto tg_xyyyyy_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 323); 

                auto tg_xyyyyy_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 324); 

                auto tg_xyyyyy_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 325); 

                auto tg_xyyyyy_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 326); 

                auto tg_xyyyyy_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 327); 

                auto tg_xyyyyy_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 328); 

                auto tg_xyyyyy_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 329); 

                auto tg_xyyyyy_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 330); 

                auto tg_xyyyyy_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 331); 

                auto tg_xyyyyy_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 332); 

                auto tg_xyyyyy_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 333); 

                auto tg_xyyyyy_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 334); 

                auto tg_xyyyyy_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 335); 

                auto tg_xyyyyz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 336); 

                auto tg_xyyyyz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 337); 

                auto tg_xyyyyz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 338); 

                auto tg_xyyyyz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 339); 

                auto tg_xyyyyz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 340); 

                auto tg_xyyyyz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 341); 

                auto tg_xyyyyz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 342); 

                auto tg_xyyyyz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 343); 

                auto tg_xyyyyz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 344); 

                auto tg_xyyyyz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 345); 

                auto tg_xyyyyz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 346); 

                auto tg_xyyyyz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 347); 

                auto tg_xyyyyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 348); 

                auto tg_xyyyyz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 349); 

                auto tg_xyyyyz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 350); 

                auto tg_xyyyyz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 351); 

                auto tg_xyyyyz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 352); 

                auto tg_xyyyyz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 353); 

                auto tg_xyyyyz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 354); 

                auto tg_xyyyyz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 355); 

                auto tg_xyyyyz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 356); 

                auto tg_xyyyzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 357); 

                auto tg_xyyyzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 358); 

                auto tg_xyyyzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 359); 

                auto tg_xyyyzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 360); 

                auto tg_xyyyzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 361); 

                auto tg_xyyyzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 362); 

                auto tg_xyyyzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 363); 

                auto tg_xyyyzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 364); 

                auto tg_xyyyzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 365); 

                auto tg_xyyyzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 366); 

                auto tg_xyyyzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 367); 

                auto tg_xyyyzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 368); 

                auto tg_xyyyzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 369); 

                auto tg_xyyyzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 370); 

                auto tg_xyyyzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 371); 

                auto tg_xyyyzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 372); 

                auto tg_xyyyzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 373); 

                auto tg_xyyyzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 374); 

                auto tg_xyyyzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 375); 

                auto tg_xyyyzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 376); 

                auto tg_xyyyzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 377); 

                auto tg_xyyzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 378); 

                auto tg_xyyzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 379); 

                auto tg_xxyzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 285); 

                auto tg_xxyzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 286); 

                auto tg_xxyzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 287); 

                auto tg_xxyzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 288); 

                auto tg_xxyzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 289); 

                auto tg_xxyzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 290); 

                auto tg_xxyzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 291); 

                auto tg_xxyzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 292); 

                auto tg_xxyzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 293); 

                auto tg_xxzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 294); 

                auto tg_xxzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 295); 

                auto tg_xxzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 296); 

                auto tg_xxzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 297); 

                auto tg_xxzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 298); 

                auto tg_xxzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 299); 

                auto tg_xxzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 300); 

                auto tg_xxzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 301); 

                auto tg_xxzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 302); 

                auto tg_xxzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 303); 

                auto tg_xxzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 304); 

                auto tg_xxzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 305); 

                auto tg_xxzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 306); 

                auto tg_xxzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 307); 

                auto tg_xxzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 308); 

                auto tg_xxzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 309); 

                auto tg_xxzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 310); 

                auto tg_xxzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 311); 

                auto tg_xxzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 312); 

                auto tg_xxzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 313); 

                auto tg_xxzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 314); 

                auto tg_xyyyyy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 315); 

                auto tg_xyyyyy_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 316); 

                auto tg_xyyyyy_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 317); 

                auto tg_xyyyyy_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 318); 

                auto tg_xyyyyy_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 319); 

                auto tg_xyyyyy_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 320); 

                auto tg_xyyyyy_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 321); 

                auto tg_xyyyyy_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 322); 

                auto tg_xyyyyy_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 323); 

                auto tg_xyyyyy_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 324); 

                auto tg_xyyyyy_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 325); 

                auto tg_xyyyyy_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 326); 

                auto tg_xyyyyy_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 327); 

                auto tg_xyyyyy_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 328); 

                auto tg_xyyyyy_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 329); 

                auto tg_xyyyyy_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 330); 

                auto tg_xyyyyy_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 331); 

                auto tg_xyyyyy_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 332); 

                auto tg_xyyyyy_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 333); 

                auto tg_xyyyyy_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 334); 

                auto tg_xyyyyy_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 335); 

                auto tg_xyyyyz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 336); 

                auto tg_xyyyyz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 337); 

                auto tg_xyyyyz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 338); 

                auto tg_xyyyyz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 339); 

                auto tg_xyyyyz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 340); 

                auto tg_xyyyyz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 341); 

                auto tg_xyyyyz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 342); 

                auto tg_xyyyyz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 343); 

                auto tg_xyyyyz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 344); 

                auto tg_xyyyyz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 345); 

                auto tg_xyyyyz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 346); 

                auto tg_xyyyyz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 347); 

                auto tg_xyyyyz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 348); 

                auto tg_xyyyyz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 349); 

                auto tg_xyyyyz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 350); 

                auto tg_xyyyyz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 351); 

                auto tg_xyyyyz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 352); 

                auto tg_xyyyyz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 353); 

                auto tg_xyyyyz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 354); 

                auto tg_xyyyyz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 355); 

                auto tg_xyyyyz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 356); 

                auto tg_xyyyzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 357); 

                auto tg_xyyyzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 358); 

                auto tg_xyyyzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 359); 

                auto tg_xyyyzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 360); 

                auto tg_xyyyzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 361); 

                auto tg_xyyyzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 362); 

                auto tg_xyyyzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 363); 

                auto tg_xyyyzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 364); 

                auto tg_xyyyzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 365); 

                auto tg_xyyyzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 366); 

                auto tg_xyyyzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 367); 

                auto tg_xyyyzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 368); 

                auto tg_xyyyzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 369); 

                auto tg_xyyyzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 370); 

                auto tg_xyyyzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 371); 

                auto tg_xyyyzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 372); 

                auto tg_xyyyzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 373); 

                auto tg_xyyyzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 374); 

                auto tg_xyyyzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 375); 

                auto tg_xyyyzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 376); 

                auto tg_xyyyzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 377); 

                auto tg_xyyzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 378); 

                auto tg_xyyzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 379); 

                auto tg_xyzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 285); 

                auto tg_xyzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 286); 

                auto tg_xyzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 287); 

                auto tg_xyzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 288); 

                auto tg_xyzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 289); 

                auto tg_xyzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 290); 

                auto tg_xyzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 291); 

                auto tg_xyzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 292); 

                auto tg_xyzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 293); 

                auto tg_xzzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 294); 

                auto tg_xzzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 295); 

                auto tg_xzzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 296); 

                auto tg_xzzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 297); 

                auto tg_xzzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 298); 

                auto tg_xzzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 299); 

                auto tg_xzzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 300); 

                auto tg_xzzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 301); 

                auto tg_xzzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 302); 

                auto tg_xzzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 303); 

                auto tg_xzzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 304); 

                auto tg_xzzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 305); 

                auto tg_xzzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 306); 

                auto tg_xzzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 307); 

                auto tg_xzzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 308); 

                auto tg_xzzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 309); 

                auto tg_xzzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 310); 

                auto tg_xzzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 311); 

                auto tg_xzzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 312); 

                auto tg_xzzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 313); 

                auto tg_xzzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 314); 

                auto tg_yyyyy_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 315); 

                auto tg_yyyyy_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 316); 

                auto tg_yyyyy_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 317); 

                auto tg_yyyyy_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 318); 

                auto tg_yyyyy_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 319); 

                auto tg_yyyyy_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 320); 

                auto tg_yyyyy_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 321); 

                auto tg_yyyyy_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 322); 

                auto tg_yyyyy_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 323); 

                auto tg_yyyyy_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 324); 

                auto tg_yyyyy_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 325); 

                auto tg_yyyyy_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 326); 

                auto tg_yyyyy_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 327); 

                auto tg_yyyyy_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 328); 

                auto tg_yyyyy_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 329); 

                auto tg_yyyyy_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 330); 

                auto tg_yyyyy_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 331); 

                auto tg_yyyyy_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 332); 

                auto tg_yyyyy_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 333); 

                auto tg_yyyyy_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 334); 

                auto tg_yyyyy_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 335); 

                auto tg_yyyyz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 336); 

                auto tg_yyyyz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 337); 

                auto tg_yyyyz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 338); 

                auto tg_yyyyz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 339); 

                auto tg_yyyyz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 340); 

                auto tg_yyyyz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 341); 

                auto tg_yyyyz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 342); 

                auto tg_yyyyz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 343); 

                auto tg_yyyyz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 344); 

                auto tg_yyyyz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 345); 

                auto tg_yyyyz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 346); 

                auto tg_yyyyz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 347); 

                auto tg_yyyyz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 348); 

                auto tg_yyyyz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 349); 

                auto tg_yyyyz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 350); 

                auto tg_yyyyz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 351); 

                auto tg_yyyyz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 352); 

                auto tg_yyyyz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 353); 

                auto tg_yyyyz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 354); 

                auto tg_yyyyz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 355); 

                auto tg_yyyyz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 356); 

                auto tg_yyyzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 357); 

                auto tg_yyyzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 358); 

                auto tg_yyyzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 359); 

                auto tg_yyyzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 360); 

                auto tg_yyyzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 361); 

                auto tg_yyyzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 362); 

                auto tg_yyyzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 363); 

                auto tg_yyyzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 364); 

                auto tg_yyyzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 365); 

                auto tg_yyyzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 366); 

                auto tg_yyyzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 367); 

                auto tg_yyyzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 368); 

                auto tg_yyyzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 369); 

                auto tg_yyyzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 370); 

                auto tg_yyyzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 371); 

                auto tg_yyyzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 372); 

                auto tg_yyyzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 373); 

                auto tg_yyyzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 374); 

                auto tg_yyyzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 375); 

                auto tg_yyyzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 376); 

                auto tg_yyyzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 377); 

                auto tg_yyzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 378); 

                auto tg_yyzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 379); 

                auto tg_xyzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 285); 

                auto tg_xyzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 286); 

                auto tg_xyzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 287); 

                auto tg_xyzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 288); 

                auto tg_xyzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 289); 

                auto tg_xyzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 290); 

                auto tg_xyzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 291); 

                auto tg_xyzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 292); 

                auto tg_xyzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 293); 

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

                auto tg_yyyzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 371); 

                auto tg_yyyzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 372); 

                auto tg_yyyzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 373); 

                auto tg_yyyzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 374); 

                auto tg_yyyzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 375); 

                auto tg_yyyzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 376); 

                auto tg_yyyzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 377); 

                auto tg_yyzzz_xxxxx_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 378); 

                auto tg_yyzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 379); 

                auto tg_xxyzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 207); 

                auto tg_xxyzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 208); 

                auto tg_xxyzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 209); 

                auto tg_xxzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 210); 

                auto tg_xxzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 211); 

                auto tg_xxzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 212); 

                auto tg_xxzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 213); 

                auto tg_xxzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 214); 

                auto tg_xxzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 215); 

                auto tg_xxzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 216); 

                auto tg_xxzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 217); 

                auto tg_xxzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 218); 

                auto tg_xxzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 219); 

                auto tg_xxzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 220); 

                auto tg_xxzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 221); 

                auto tg_xxzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 222); 

                auto tg_xxzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 223); 

                auto tg_xxzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 224); 

                auto tg_xyyyyy_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 225); 

                auto tg_xyyyyy_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 226); 

                auto tg_xyyyyy_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 227); 

                auto tg_xyyyyy_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 228); 

                auto tg_xyyyyy_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 229); 

                auto tg_xyyyyy_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 230); 

                auto tg_xyyyyy_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 231); 

                auto tg_xyyyyy_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 232); 

                auto tg_xyyyyy_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 233); 

                auto tg_xyyyyy_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 234); 

                auto tg_xyyyyy_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 235); 

                auto tg_xyyyyy_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 236); 

                auto tg_xyyyyy_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 237); 

                auto tg_xyyyyy_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 238); 

                auto tg_xyyyyy_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 239); 

                auto tg_xyyyyz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 240); 

                auto tg_xyyyyz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 241); 

                auto tg_xyyyyz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 242); 

                auto tg_xyyyyz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 243); 

                auto tg_xyyyyz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 244); 

                auto tg_xyyyyz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 245); 

                auto tg_xyyyyz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 246); 

                auto tg_xyyyyz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 247); 

                auto tg_xyyyyz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 248); 

                auto tg_xyyyyz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 249); 

                auto tg_xyyyyz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 250); 

                auto tg_xyyyyz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 251); 

                auto tg_xyyyyz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 252); 

                auto tg_xyyyyz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 253); 

                auto tg_xyyyyz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 254); 

                auto tg_xyyyzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 255); 

                auto tg_xyyyzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 256); 

                auto tg_xyyyzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 257); 

                auto tg_xyyyzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 258); 

                auto tg_xyyyzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 259); 

                auto tg_xyyyzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 260); 

                auto tg_xyyyzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 261); 

                auto tg_xyyyzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 262); 

                auto tg_xyyyzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 263); 

                auto tg_xyyyzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 264); 

                auto tg_xyyyzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 265); 

                auto tg_xyyyzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 266); 

                auto tg_xyyyzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 267); 

                auto tg_xyyyzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 268); 

                auto tg_xyyyzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 269); 

                auto tg_xyyzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 270); 

                auto tg_xyyzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 271); 

                // set up pointers to integrals

                auto tg_xxxyzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 285); 

                auto tg_xxxyzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 286); 

                auto tg_xxxyzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 287); 

                auto tg_xxxyzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 288); 

                auto tg_xxxyzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 289); 

                auto tg_xxxyzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 290); 

                auto tg_xxxyzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 291); 

                auto tg_xxxyzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 292); 

                auto tg_xxxyzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 293); 

                auto tg_xxxzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 294); 

                auto tg_xxxzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 295); 

                auto tg_xxxzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 296); 

                auto tg_xxxzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 297); 

                auto tg_xxxzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 298); 

                auto tg_xxxzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 299); 

                auto tg_xxxzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 300); 

                auto tg_xxxzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 301); 

                auto tg_xxxzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 302); 

                auto tg_xxxzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 303); 

                auto tg_xxxzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 304); 

                auto tg_xxxzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 305); 

                auto tg_xxxzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 306); 

                auto tg_xxxzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 307); 

                auto tg_xxxzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 308); 

                auto tg_xxxzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 309); 

                auto tg_xxxzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 310); 

                auto tg_xxxzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 311); 

                auto tg_xxxzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 312); 

                auto tg_xxxzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 313); 

                auto tg_xxxzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 314); 

                auto tg_xxyyyyy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 315); 

                auto tg_xxyyyyy_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 316); 

                auto tg_xxyyyyy_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 317); 

                auto tg_xxyyyyy_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 318); 

                auto tg_xxyyyyy_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 319); 

                auto tg_xxyyyyy_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 320); 

                auto tg_xxyyyyy_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 321); 

                auto tg_xxyyyyy_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 322); 

                auto tg_xxyyyyy_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 323); 

                auto tg_xxyyyyy_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 324); 

                auto tg_xxyyyyy_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 325); 

                auto tg_xxyyyyy_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 326); 

                auto tg_xxyyyyy_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 327); 

                auto tg_xxyyyyy_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 328); 

                auto tg_xxyyyyy_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 329); 

                auto tg_xxyyyyy_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 330); 

                auto tg_xxyyyyy_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 331); 

                auto tg_xxyyyyy_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 332); 

                auto tg_xxyyyyy_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 333); 

                auto tg_xxyyyyy_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 334); 

                auto tg_xxyyyyy_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 335); 

                auto tg_xxyyyyz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 336); 

                auto tg_xxyyyyz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 337); 

                auto tg_xxyyyyz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 338); 

                auto tg_xxyyyyz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 339); 

                auto tg_xxyyyyz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 340); 

                auto tg_xxyyyyz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 341); 

                auto tg_xxyyyyz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 342); 

                auto tg_xxyyyyz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 343); 

                auto tg_xxyyyyz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 344); 

                auto tg_xxyyyyz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 345); 

                auto tg_xxyyyyz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 346); 

                auto tg_xxyyyyz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 347); 

                auto tg_xxyyyyz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 348); 

                auto tg_xxyyyyz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 349); 

                auto tg_xxyyyyz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 350); 

                auto tg_xxyyyyz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 351); 

                auto tg_xxyyyyz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 352); 

                auto tg_xxyyyyz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 353); 

                auto tg_xxyyyyz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 354); 

                auto tg_xxyyyyz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 355); 

                auto tg_xxyyyyz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 356); 

                auto tg_xxyyyzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 357); 

                auto tg_xxyyyzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 358); 

                auto tg_xxyyyzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 359); 

                auto tg_xxyyyzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 360); 

                auto tg_xxyyyzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 361); 

                auto tg_xxyyyzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 362); 

                auto tg_xxyyyzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 363); 

                auto tg_xxyyyzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 364); 

                auto tg_xxyyyzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 365); 

                auto tg_xxyyyzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 366); 

                auto tg_xxyyyzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 367); 

                auto tg_xxyyyzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 368); 

                auto tg_xxyyyzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 369); 

                auto tg_xxyyyzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 370); 

                auto tg_xxyyyzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 371); 

                auto tg_xxyyyzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 372); 

                auto tg_xxyyyzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 373); 

                auto tg_xxyyyzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 374); 

                auto tg_xxyyyzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 375); 

                auto tg_xxyyyzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 376); 

                auto tg_xxyyyzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 377); 

                auto tg_xxyyzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 378); 

                auto tg_xxyyzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 379); 

                // Batch of Integrals (285,380)

                #pragma omp simd aligned(fxn, fza, tg_xxxyzzz_xyyzz_0, tg_xxxyzzz_xyzzz_0, tg_xxxyzzz_xzzzz_0, \
                                         tg_xxxyzzz_yyyyy_0, tg_xxxyzzz_yyyyz_0, tg_xxxyzzz_yyyzz_0, tg_xxxyzzz_yyzzz_0, \
                                         tg_xxxyzzz_yzzzz_0, tg_xxxyzzz_zzzzz_0, tg_xxxzzzz_xxxxx_0, tg_xxxzzzz_xxxxy_0, \
                                         tg_xxxzzzz_xxxxz_0, tg_xxxzzzz_xxxyy_0, tg_xxxzzzz_xxxyz_0, tg_xxxzzzz_xxxzz_0, \
                                         tg_xxxzzzz_xxyyy_0, tg_xxxzzzz_xxyyz_0, tg_xxxzzzz_xxyzz_0, tg_xxxzzzz_xxzzz_0, \
                                         tg_xxxzzzz_xyyyy_0, tg_xxxzzzz_xyyyz_0, tg_xxxzzzz_xyyzz_0, tg_xxxzzzz_xyzzz_0, \
                                         tg_xxxzzzz_xzzzz_0, tg_xxxzzzz_yyyyy_0, tg_xxxzzzz_yyyyz_0, tg_xxxzzzz_yyyzz_0, \
                                         tg_xxxzzzz_yyzzz_0, tg_xxxzzzz_yzzzz_0, tg_xxxzzzz_zzzzz_0, tg_xxyyyyy_xxxxx_0, \
                                         tg_xxyyyyy_xxxxy_0, tg_xxyyyyy_xxxxz_0, tg_xxyyyyy_xxxyy_0, tg_xxyyyyy_xxxyz_0, \
                                         tg_xxyyyyy_xxxzz_0, tg_xxyyyyy_xxyyy_0, tg_xxyyyyy_xxyyz_0, tg_xxyyyyy_xxyzz_0, \
                                         tg_xxyyyyy_xxzzz_0, tg_xxyyyyy_xyyyy_0, tg_xxyyyyy_xyyyz_0, tg_xxyyyyy_xyyzz_0, \
                                         tg_xxyyyyy_xyzzz_0, tg_xxyyyyy_xzzzz_0, tg_xxyyyyy_yyyyy_0, tg_xxyyyyy_yyyyz_0, \
                                         tg_xxyyyyy_yyyzz_0, tg_xxyyyyy_yyzzz_0, tg_xxyyyyy_yzzzz_0, tg_xxyyyyy_zzzzz_0, \
                                         tg_xxyyyyz_xxxxx_0, tg_xxyyyyz_xxxxy_0, tg_xxyyyyz_xxxxz_0, tg_xxyyyyz_xxxyy_0, \
                                         tg_xxyyyyz_xxxyz_0, tg_xxyyyyz_xxxzz_0, tg_xxyyyyz_xxyyy_0, tg_xxyyyyz_xxyyz_0, \
                                         tg_xxyyyyz_xxyzz_0, tg_xxyyyyz_xxzzz_0, tg_xxyyyyz_xyyyy_0, tg_xxyyyyz_xyyyz_0, \
                                         tg_xxyyyyz_xyyzz_0, tg_xxyyyyz_xyzzz_0, tg_xxyyyyz_xzzzz_0, tg_xxyyyyz_yyyyy_0, \
                                         tg_xxyyyyz_yyyyz_0, tg_xxyyyyz_yyyzz_0, tg_xxyyyyz_yyzzz_0, tg_xxyyyyz_yzzzz_0, \
                                         tg_xxyyyyz_zzzzz_0, tg_xxyyyzz_xxxxx_0, tg_xxyyyzz_xxxxy_0, tg_xxyyyzz_xxxxz_0, \
                                         tg_xxyyyzz_xxxyy_0, tg_xxyyyzz_xxxyz_0, tg_xxyyyzz_xxxzz_0, tg_xxyyyzz_xxyyy_0, \
                                         tg_xxyyyzz_xxyyz_0, tg_xxyyyzz_xxyzz_0, tg_xxyyyzz_xxzzz_0, tg_xxyyyzz_xyyyy_0, \
                                         tg_xxyyyzz_xyyyz_0, tg_xxyyyzz_xyyzz_0, tg_xxyyyzz_xyzzz_0, tg_xxyyyzz_xzzzz_0, \
                                         tg_xxyyyzz_yyyyy_0, tg_xxyyyzz_yyyyz_0, tg_xxyyyzz_yyyzz_0, tg_xxyyyzz_yyzzz_0, \
                                         tg_xxyyyzz_yzzzz_0, tg_xxyyyzz_zzzzz_0, tg_xxyyzzz_xxxxx_0, tg_xxyyzzz_xxxxy_0, \
                                         tg_xxyzzz_xyyzz_0, tg_xxyzzz_xyyzz_1, tg_xxyzzz_xyzzz_0, tg_xxyzzz_xyzzz_1, \
                                         tg_xxyzzz_xzzzz_0, tg_xxyzzz_xzzzz_1, tg_xxyzzz_yyyyy_0, tg_xxyzzz_yyyyy_1, \
                                         tg_xxyzzz_yyyyz_0, tg_xxyzzz_yyyyz_1, tg_xxyzzz_yyyzz_0, tg_xxyzzz_yyyzz_1, \
                                         tg_xxyzzz_yyzz_1, tg_xxyzzz_yyzzz_0, tg_xxyzzz_yyzzz_1, tg_xxyzzz_yzzz_1, \
                                         tg_xxyzzz_yzzzz_0, tg_xxyzzz_yzzzz_1, tg_xxyzzz_zzzz_1, tg_xxyzzz_zzzzz_0, \
                                         tg_xxyzzz_zzzzz_1, tg_xxzzzz_xxxx_1, tg_xxzzzz_xxxxx_0, tg_xxzzzz_xxxxx_1, \
                                         tg_xxzzzz_xxxxy_0, tg_xxzzzz_xxxxy_1, tg_xxzzzz_xxxxz_0, tg_xxzzzz_xxxxz_1, \
                                         tg_xxzzzz_xxxy_1, tg_xxzzzz_xxxyy_0, tg_xxzzzz_xxxyy_1, tg_xxzzzz_xxxyz_0, \
                                         tg_xxzzzz_xxxyz_1, tg_xxzzzz_xxxz_1, tg_xxzzzz_xxxzz_0, tg_xxzzzz_xxxzz_1, \
                                         tg_xxzzzz_xxyy_1, tg_xxzzzz_xxyyy_0, tg_xxzzzz_xxyyy_1, tg_xxzzzz_xxyyz_0, \
                                         tg_xxzzzz_xxyyz_1, tg_xxzzzz_xxyz_1, tg_xxzzzz_xxyzz_0, tg_xxzzzz_xxyzz_1, \
                                         tg_xxzzzz_xxzz_1, tg_xxzzzz_xxzzz_0, tg_xxzzzz_xxzzz_1, tg_xxzzzz_xyyy_1, \
                                         tg_xxzzzz_xyyyy_0, tg_xxzzzz_xyyyy_1, tg_xxzzzz_xyyyz_0, tg_xxzzzz_xyyyz_1, \
                                         tg_xxzzzz_xyyz_1, tg_xxzzzz_xyyzz_0, tg_xxzzzz_xyyzz_1, tg_xxzzzz_xyzz_1, \
                                         tg_xxzzzz_xyzzz_0, tg_xxzzzz_xyzzz_1, tg_xxzzzz_xzzz_1, tg_xxzzzz_xzzzz_0, \
                                         tg_xxzzzz_xzzzz_1, tg_xxzzzz_yyyy_1, tg_xxzzzz_yyyyy_0, tg_xxzzzz_yyyyy_1, \
                                         tg_xxzzzz_yyyyz_0, tg_xxzzzz_yyyyz_1, tg_xxzzzz_yyyz_1, tg_xxzzzz_yyyzz_0, \
                                         tg_xxzzzz_yyyzz_1, tg_xxzzzz_yyzz_1, tg_xxzzzz_yyzzz_0, tg_xxzzzz_yyzzz_1, \
                                         tg_xxzzzz_yzzz_1, tg_xxzzzz_yzzzz_0, tg_xxzzzz_yzzzz_1, tg_xxzzzz_zzzz_1, \
                                         tg_xxzzzz_zzzzz_0, tg_xxzzzz_zzzzz_1, tg_xyyyyy_xxxx_1, tg_xyyyyy_xxxxx_0, \
                                         tg_xyyyyy_xxxxx_1, tg_xyyyyy_xxxxy_0, tg_xyyyyy_xxxxy_1, tg_xyyyyy_xxxxz_0, \
                                         tg_xyyyyy_xxxxz_1, tg_xyyyyy_xxxy_1, tg_xyyyyy_xxxyy_0, tg_xyyyyy_xxxyy_1, \
                                         tg_xyyyyy_xxxyz_0, tg_xyyyyy_xxxyz_1, tg_xyyyyy_xxxz_1, tg_xyyyyy_xxxzz_0, \
                                         tg_xyyyyy_xxxzz_1, tg_xyyyyy_xxyy_1, tg_xyyyyy_xxyyy_0, tg_xyyyyy_xxyyy_1, \
                                         tg_xyyyyy_xxyyz_0, tg_xyyyyy_xxyyz_1, tg_xyyyyy_xxyz_1, tg_xyyyyy_xxyzz_0, \
                                         tg_xyyyyy_xxyzz_1, tg_xyyyyy_xxzz_1, tg_xyyyyy_xxzzz_0, tg_xyyyyy_xxzzz_1, \
                                         tg_xyyyyy_xyyy_1, tg_xyyyyy_xyyyy_0, tg_xyyyyy_xyyyy_1, tg_xyyyyy_xyyyz_0, \
                                         tg_xyyyyy_xyyyz_1, tg_xyyyyy_xyyz_1, tg_xyyyyy_xyyzz_0, tg_xyyyyy_xyyzz_1, \
                                         tg_xyyyyy_xyzz_1, tg_xyyyyy_xyzzz_0, tg_xyyyyy_xyzzz_1, tg_xyyyyy_xzzz_1, \
                                         tg_xyyyyy_xzzzz_0, tg_xyyyyy_xzzzz_1, tg_xyyyyy_yyyy_1, tg_xyyyyy_yyyyy_0, \
                                         tg_xyyyyy_yyyyy_1, tg_xyyyyy_yyyyz_0, tg_xyyyyy_yyyyz_1, tg_xyyyyy_yyyz_1, \
                                         tg_xyyyyy_yyyzz_0, tg_xyyyyy_yyyzz_1, tg_xyyyyy_yyzz_1, tg_xyyyyy_yyzzz_0, \
                                         tg_xyyyyy_yyzzz_1, tg_xyyyyy_yzzz_1, tg_xyyyyy_yzzzz_0, tg_xyyyyy_yzzzz_1, \
                                         tg_xyyyyy_zzzz_1, tg_xyyyyy_zzzzz_0, tg_xyyyyy_zzzzz_1, tg_xyyyyz_xxxx_1, \
                                         tg_xyyyyz_xxxxx_0, tg_xyyyyz_xxxxx_1, tg_xyyyyz_xxxxy_0, tg_xyyyyz_xxxxy_1, \
                                         tg_xyyyyz_xxxxz_0, tg_xyyyyz_xxxxz_1, tg_xyyyyz_xxxy_1, tg_xyyyyz_xxxyy_0, \
                                         tg_xyyyyz_xxxyy_1, tg_xyyyyz_xxxyz_0, tg_xyyyyz_xxxyz_1, tg_xyyyyz_xxxz_1, \
                                         tg_xyyyyz_xxxzz_0, tg_xyyyyz_xxxzz_1, tg_xyyyyz_xxyy_1, tg_xyyyyz_xxyyy_0, \
                                         tg_xyyyyz_xxyyy_1, tg_xyyyyz_xxyyz_0, tg_xyyyyz_xxyyz_1, tg_xyyyyz_xxyz_1, \
                                         tg_xyyyyz_xxyzz_0, tg_xyyyyz_xxyzz_1, tg_xyyyyz_xxzz_1, tg_xyyyyz_xxzzz_0, \
                                         tg_xyyyyz_xxzzz_1, tg_xyyyyz_xyyy_1, tg_xyyyyz_xyyyy_0, tg_xyyyyz_xyyyy_1, \
                                         tg_xyyyyz_xyyyz_0, tg_xyyyyz_xyyyz_1, tg_xyyyyz_xyyz_1, tg_xyyyyz_xyyzz_0, \
                                         tg_xyyyyz_xyyzz_1, tg_xyyyyz_xyzz_1, tg_xyyyyz_xyzzz_0, tg_xyyyyz_xyzzz_1, \
                                         tg_xyyyyz_xzzz_1, tg_xyyyyz_xzzzz_0, tg_xyyyyz_xzzzz_1, tg_xyyyyz_yyyy_1, \
                                         tg_xyyyyz_yyyyy_0, tg_xyyyyz_yyyyy_1, tg_xyyyyz_yyyyz_0, tg_xyyyyz_yyyyz_1, \
                                         tg_xyyyyz_yyyz_1, tg_xyyyyz_yyyzz_0, tg_xyyyyz_yyyzz_1, tg_xyyyyz_yyzz_1, \
                                         tg_xyyyyz_yyzzz_0, tg_xyyyyz_yyzzz_1, tg_xyyyyz_yzzz_1, tg_xyyyyz_yzzzz_0, \
                                         tg_xyyyyz_yzzzz_1, tg_xyyyyz_zzzz_1, tg_xyyyyz_zzzzz_0, tg_xyyyyz_zzzzz_1, \
                                         tg_xyyyzz_xxxx_1, tg_xyyyzz_xxxxx_0, tg_xyyyzz_xxxxx_1, tg_xyyyzz_xxxxy_0, \
                                         tg_xyyyzz_xxxxy_1, tg_xyyyzz_xxxxz_0, tg_xyyyzz_xxxxz_1, tg_xyyyzz_xxxy_1, \
                                         tg_xyyyzz_xxxyy_0, tg_xyyyzz_xxxyy_1, tg_xyyyzz_xxxyz_0, tg_xyyyzz_xxxyz_1, \
                                         tg_xyyyzz_xxxz_1, tg_xyyyzz_xxxzz_0, tg_xyyyzz_xxxzz_1, tg_xyyyzz_xxyy_1, \
                                         tg_xyyyzz_xxyyy_0, tg_xyyyzz_xxyyy_1, tg_xyyyzz_xxyyz_0, tg_xyyyzz_xxyyz_1, \
                                         tg_xyyyzz_xxyz_1, tg_xyyyzz_xxyzz_0, tg_xyyyzz_xxyzz_1, tg_xyyyzz_xxzz_1, \
                                         tg_xyyyzz_xxzzz_0, tg_xyyyzz_xxzzz_1, tg_xyyyzz_xyyy_1, tg_xyyyzz_xyyyy_0, \
                                         tg_xyyyzz_xyyyy_1, tg_xyyyzz_xyyyz_0, tg_xyyyzz_xyyyz_1, tg_xyyyzz_xyyz_1, \
                                         tg_xyyyzz_xyyzz_0, tg_xyyyzz_xyyzz_1, tg_xyyyzz_xyzz_1, tg_xyyyzz_xyzzz_0, \
                                         tg_xyyyzz_xyzzz_1, tg_xyyyzz_xzzz_1, tg_xyyyzz_xzzzz_0, tg_xyyyzz_xzzzz_1, \
                                         tg_xyyyzz_yyyy_1, tg_xyyyzz_yyyyy_0, tg_xyyyzz_yyyyy_1, tg_xyyyzz_yyyyz_0, \
                                         tg_xyyyzz_yyyyz_1, tg_xyyyzz_yyyz_1, tg_xyyyzz_yyyzz_0, tg_xyyyzz_yyyzz_1, \
                                         tg_xyyyzz_yyzz_1, tg_xyyyzz_yyzzz_0, tg_xyyyzz_yyzzz_1, tg_xyyyzz_yzzz_1, \
                                         tg_xyyyzz_yzzzz_0, tg_xyyyzz_yzzzz_1, tg_xyyyzz_zzzz_1, tg_xyyyzz_zzzzz_0, \
                                         tg_xyyyzz_zzzzz_1, tg_xyyzzz_xxxx_1, tg_xyyzzz_xxxxx_0, tg_xyyzzz_xxxxx_1, \
                                         tg_xyyzzz_xxxxy_0, tg_xyyzzz_xxxxy_1, tg_xyyzzz_xxxy_1, tg_xyzzz_xyyzz_0, \
                                         tg_xyzzz_xyyzz_1, tg_xyzzz_xyzzz_0, tg_xyzzz_xyzzz_1, tg_xyzzz_xzzzz_0, \
                                         tg_xyzzz_xzzzz_1, tg_xyzzz_yyyyy_0, tg_xyzzz_yyyyy_1, tg_xyzzz_yyyyz_0, \
                                         tg_xyzzz_yyyyz_1, tg_xyzzz_yyyzz_0, tg_xyzzz_yyyzz_1, tg_xyzzz_yyzzz_0, \
                                         tg_xyzzz_yyzzz_1, tg_xyzzz_yzzzz_0, tg_xyzzz_yzzzz_1, tg_xyzzz_zzzzz_0, \
                                         tg_xyzzz_zzzzz_1, tg_xzzzz_xxxxx_0, tg_xzzzz_xxxxx_1, tg_xzzzz_xxxxy_0, \
                                         tg_xzzzz_xxxxy_1, tg_xzzzz_xxxxz_0, tg_xzzzz_xxxxz_1, tg_xzzzz_xxxyy_0, \
                                         tg_xzzzz_xxxyy_1, tg_xzzzz_xxxyz_0, tg_xzzzz_xxxyz_1, tg_xzzzz_xxxzz_0, \
                                         tg_xzzzz_xxxzz_1, tg_xzzzz_xxyyy_0, tg_xzzzz_xxyyy_1, tg_xzzzz_xxyyz_0, \
                                         tg_xzzzz_xxyyz_1, tg_xzzzz_xxyzz_0, tg_xzzzz_xxyzz_1, tg_xzzzz_xxzzz_0, \
                                         tg_xzzzz_xxzzz_1, tg_xzzzz_xyyyy_0, tg_xzzzz_xyyyy_1, tg_xzzzz_xyyyz_0, \
                                         tg_xzzzz_xyyyz_1, tg_xzzzz_xyyzz_0, tg_xzzzz_xyyzz_1, tg_xzzzz_xyzzz_0, \
                                         tg_xzzzz_xyzzz_1, tg_xzzzz_xzzzz_0, tg_xzzzz_xzzzz_1, tg_xzzzz_yyyyy_0, \
                                         tg_xzzzz_yyyyy_1, tg_xzzzz_yyyyz_0, tg_xzzzz_yyyyz_1, tg_xzzzz_yyyzz_0, \
                                         tg_xzzzz_yyyzz_1, tg_xzzzz_yyzzz_0, tg_xzzzz_yyzzz_1, tg_xzzzz_yzzzz_0, \
                                         tg_xzzzz_yzzzz_1, tg_xzzzz_zzzzz_0, tg_xzzzz_zzzzz_1, tg_yyyyy_xxxxx_0, \
                                         tg_yyyyy_xxxxx_1, tg_yyyyy_xxxxy_0, tg_yyyyy_xxxxy_1, tg_yyyyy_xxxxz_0, \
                                         tg_yyyyy_xxxxz_1, tg_yyyyy_xxxyy_0, tg_yyyyy_xxxyy_1, tg_yyyyy_xxxyz_0, \
                                         tg_yyyyy_xxxyz_1, tg_yyyyy_xxxzz_0, tg_yyyyy_xxxzz_1, tg_yyyyy_xxyyy_0, \
                                         tg_yyyyy_xxyyy_1, tg_yyyyy_xxyyz_0, tg_yyyyy_xxyyz_1, tg_yyyyy_xxyzz_0, \
                                         tg_yyyyy_xxyzz_1, tg_yyyyy_xxzzz_0, tg_yyyyy_xxzzz_1, tg_yyyyy_xyyyy_0, \
                                         tg_yyyyy_xyyyy_1, tg_yyyyy_xyyyz_0, tg_yyyyy_xyyyz_1, tg_yyyyy_xyyzz_0, \
                                         tg_yyyyy_xyyzz_1, tg_yyyyy_xyzzz_0, tg_yyyyy_xyzzz_1, tg_yyyyy_xzzzz_0, \
                                         tg_yyyyy_xzzzz_1, tg_yyyyy_yyyyy_0, tg_yyyyy_yyyyy_1, tg_yyyyy_yyyyz_0, \
                                         tg_yyyyy_yyyyz_1, tg_yyyyy_yyyzz_0, tg_yyyyy_yyyzz_1, tg_yyyyy_yyzzz_0, \
                                         tg_yyyyy_yyzzz_1, tg_yyyyy_yzzzz_0, tg_yyyyy_yzzzz_1, tg_yyyyy_zzzzz_0, \
                                         tg_yyyyy_zzzzz_1, tg_yyyyz_xxxxx_0, tg_yyyyz_xxxxx_1, tg_yyyyz_xxxxy_0, \
                                         tg_yyyyz_xxxxy_1, tg_yyyyz_xxxxz_0, tg_yyyyz_xxxxz_1, tg_yyyyz_xxxyy_0, \
                                         tg_yyyyz_xxxyy_1, tg_yyyyz_xxxyz_0, tg_yyyyz_xxxyz_1, tg_yyyyz_xxxzz_0, \
                                         tg_yyyyz_xxxzz_1, tg_yyyyz_xxyyy_0, tg_yyyyz_xxyyy_1, tg_yyyyz_xxyyz_0, \
                                         tg_yyyyz_xxyyz_1, tg_yyyyz_xxyzz_0, tg_yyyyz_xxyzz_1, tg_yyyyz_xxzzz_0, \
                                         tg_yyyyz_xxzzz_1, tg_yyyyz_xyyyy_0, tg_yyyyz_xyyyy_1, tg_yyyyz_xyyyz_0, \
                                         tg_yyyyz_xyyyz_1, tg_yyyyz_xyyzz_0, tg_yyyyz_xyyzz_1, tg_yyyyz_xyzzz_0, \
                                         tg_yyyyz_xyzzz_1, tg_yyyyz_xzzzz_0, tg_yyyyz_xzzzz_1, tg_yyyyz_yyyyy_0, \
                                         tg_yyyyz_yyyyy_1, tg_yyyyz_yyyyz_0, tg_yyyyz_yyyyz_1, tg_yyyyz_yyyzz_0, \
                                         tg_yyyyz_yyyzz_1, tg_yyyyz_yyzzz_0, tg_yyyyz_yyzzz_1, tg_yyyyz_yzzzz_0, \
                                         tg_yyyyz_yzzzz_1, tg_yyyyz_zzzzz_0, tg_yyyyz_zzzzz_1, tg_yyyzz_xxxxx_0, \
                                         tg_yyyzz_xxxxx_1, tg_yyyzz_xxxxy_0, tg_yyyzz_xxxxy_1, tg_yyyzz_xxxxz_0, \
                                         tg_yyyzz_xxxxz_1, tg_yyyzz_xxxyy_0, tg_yyyzz_xxxyy_1, tg_yyyzz_xxxyz_0, \
                                         tg_yyyzz_xxxyz_1, tg_yyyzz_xxxzz_0, tg_yyyzz_xxxzz_1, tg_yyyzz_xxyyy_0, \
                                         tg_yyyzz_xxyyy_1, tg_yyyzz_xxyyz_0, tg_yyyzz_xxyyz_1, tg_yyyzz_xxyzz_0, \
                                         tg_yyyzz_xxyzz_1, tg_yyyzz_xxzzz_0, tg_yyyzz_xxzzz_1, tg_yyyzz_xyyyy_0, \
                                         tg_yyyzz_xyyyy_1, tg_yyyzz_xyyyz_0, tg_yyyzz_xyyyz_1, tg_yyyzz_xyyzz_0, \
                                         tg_yyyzz_xyyzz_1, tg_yyyzz_xyzzz_0, tg_yyyzz_xyzzz_1, tg_yyyzz_xzzzz_0, \
                                         tg_yyyzz_xzzzz_1, tg_yyyzz_yyyyy_0, tg_yyyzz_yyyyy_1, tg_yyyzz_yyyyz_0, \
                                         tg_yyyzz_yyyyz_1, tg_yyyzz_yyyzz_0, tg_yyyzz_yyyzz_1, tg_yyyzz_yyzzz_0, \
                                         tg_yyyzz_yyzzz_1, tg_yyyzz_yzzzz_0, tg_yyyzz_yzzzz_1, tg_yyyzz_zzzzz_0, \
                                         tg_yyyzz_zzzzz_1, tg_yyzzz_xxxxx_0, tg_yyzzz_xxxxx_1, tg_yyzzz_xxxxy_0, \
                                         tg_yyzzz_xxxxy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyzzz_xyyzz_0[j] = pb_x * tg_xxyzzz_xyyzz_0[j] + fr * tg_xxyzzz_xyyzz_1[j] + fl1_fx * (tg_xyzzz_xyyzz_0[j] - tg_xyzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyzz_1[j];

                    tg_xxxyzzz_xyzzz_0[j] = pb_x * tg_xxyzzz_xyzzz_0[j] + fr * tg_xxyzzz_xyzzz_1[j] + fl1_fx * (tg_xyzzz_xyzzz_0[j] - tg_xyzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yzzz_1[j];

                    tg_xxxyzzz_xzzzz_0[j] = pb_x * tg_xxyzzz_xzzzz_0[j] + fr * tg_xxyzzz_xzzzz_1[j] + fl1_fx * (tg_xyzzz_xzzzz_0[j] - tg_xyzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_zzzz_1[j];

                    tg_xxxyzzz_yyyyy_0[j] = pb_x * tg_xxyzzz_yyyyy_0[j] + fr * tg_xxyzzz_yyyyy_1[j] + fl1_fx * (tg_xyzzz_yyyyy_0[j] - tg_xyzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyzzz_yyyyz_0[j] = pb_x * tg_xxyzzz_yyyyz_0[j] + fr * tg_xxyzzz_yyyyz_1[j] + fl1_fx * (tg_xyzzz_yyyyz_0[j] - tg_xyzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyzzz_yyyzz_0[j] = pb_x * tg_xxyzzz_yyyzz_0[j] + fr * tg_xxyzzz_yyyzz_1[j] + fl1_fx * (tg_xyzzz_yyyzz_0[j] - tg_xyzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyzzz_yyzzz_0[j] = pb_x * tg_xxyzzz_yyzzz_0[j] + fr * tg_xxyzzz_yyzzz_1[j] + fl1_fx * (tg_xyzzz_yyzzz_0[j] - tg_xyzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_yzzzz_0[j] = pb_x * tg_xxyzzz_yzzzz_0[j] + fr * tg_xxyzzz_yzzzz_1[j] + fl1_fx * (tg_xyzzz_yzzzz_0[j] - tg_xyzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_zzzzz_0[j] = pb_x * tg_xxyzzz_zzzzz_0[j] + fr * tg_xxyzzz_zzzzz_1[j] + fl1_fx * (tg_xyzzz_zzzzz_0[j] - tg_xyzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_xxxxx_0[j] = pb_x * tg_xxzzzz_xxxxx_0[j] + fr * tg_xxzzzz_xxxxx_1[j] + fl1_fx * (tg_xzzzz_xxxxx_0[j] - tg_xzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzzzz_xxxx_1[j];

                    tg_xxxzzzz_xxxxy_0[j] = pb_x * tg_xxzzzz_xxxxy_0[j] + fr * tg_xxzzzz_xxxxy_1[j] + fl1_fx * (tg_xzzzz_xxxxy_0[j] - tg_xzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzz_xxxy_1[j];

                    tg_xxxzzzz_xxxxz_0[j] = pb_x * tg_xxzzzz_xxxxz_0[j] + fr * tg_xxzzzz_xxxxz_1[j] + fl1_fx * (tg_xzzzz_xxxxz_0[j] - tg_xzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzz_xxxz_1[j];

                    tg_xxxzzzz_xxxyy_0[j] = pb_x * tg_xxzzzz_xxxyy_0[j] + fr * tg_xxzzzz_xxxyy_1[j] + fl1_fx * (tg_xzzzz_xxxyy_0[j] - tg_xzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxyy_1[j];

                    tg_xxxzzzz_xxxyz_0[j] = pb_x * tg_xxzzzz_xxxyz_0[j] + fr * tg_xxzzzz_xxxyz_1[j] + fl1_fx * (tg_xzzzz_xxxyz_0[j] - tg_xzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxyz_1[j];

                    tg_xxxzzzz_xxxzz_0[j] = pb_x * tg_xxzzzz_xxxzz_0[j] + fr * tg_xxzzzz_xxxzz_1[j] + fl1_fx * (tg_xzzzz_xxxzz_0[j] - tg_xzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxzz_1[j];

                    tg_xxxzzzz_xxyyy_0[j] = pb_x * tg_xxzzzz_xxyyy_0[j] + fr * tg_xxzzzz_xxyyy_1[j] + fl1_fx * (tg_xzzzz_xxyyy_0[j] - tg_xzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyyy_1[j];

                    tg_xxxzzzz_xxyyz_0[j] = pb_x * tg_xxzzzz_xxyyz_0[j] + fr * tg_xxzzzz_xxyyz_1[j] + fl1_fx * (tg_xzzzz_xxyyz_0[j] - tg_xzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyyz_1[j];

                    tg_xxxzzzz_xxyzz_0[j] = pb_x * tg_xxzzzz_xxyzz_0[j] + fr * tg_xxzzzz_xxyzz_1[j] + fl1_fx * (tg_xzzzz_xxyzz_0[j] - tg_xzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyzz_1[j];

                    tg_xxxzzzz_xxzzz_0[j] = pb_x * tg_xxzzzz_xxzzz_0[j] + fr * tg_xxzzzz_xxzzz_1[j] + fl1_fx * (tg_xzzzz_xxzzz_0[j] - tg_xzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xzzz_1[j];

                    tg_xxxzzzz_xyyyy_0[j] = pb_x * tg_xxzzzz_xyyyy_0[j] + fr * tg_xxzzzz_xyyyy_1[j] + fl1_fx * (tg_xzzzz_xyyyy_0[j] - tg_xzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyyy_1[j];

                    tg_xxxzzzz_xyyyz_0[j] = pb_x * tg_xxzzzz_xyyyz_0[j] + fr * tg_xxzzzz_xyyyz_1[j] + fl1_fx * (tg_xzzzz_xyyyz_0[j] - tg_xzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyyz_1[j];

                    tg_xxxzzzz_xyyzz_0[j] = pb_x * tg_xxzzzz_xyyzz_0[j] + fr * tg_xxzzzz_xyyzz_1[j] + fl1_fx * (tg_xzzzz_xyyzz_0[j] - tg_xzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyzz_1[j];

                    tg_xxxzzzz_xyzzz_0[j] = pb_x * tg_xxzzzz_xyzzz_0[j] + fr * tg_xxzzzz_xyzzz_1[j] + fl1_fx * (tg_xzzzz_xyzzz_0[j] - tg_xzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yzzz_1[j];

                    tg_xxxzzzz_xzzzz_0[j] = pb_x * tg_xxzzzz_xzzzz_0[j] + fr * tg_xxzzzz_xzzzz_1[j] + fl1_fx * (tg_xzzzz_xzzzz_0[j] - tg_xzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_zzzz_1[j];

                    tg_xxxzzzz_yyyyy_0[j] = pb_x * tg_xxzzzz_yyyyy_0[j] + fr * tg_xxzzzz_yyyyy_1[j] + fl1_fx * (tg_xzzzz_yyyyy_0[j] - tg_xzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxzzzz_yyyyz_0[j] = pb_x * tg_xxzzzz_yyyyz_0[j] + fr * tg_xxzzzz_yyyyz_1[j] + fl1_fx * (tg_xzzzz_yyyyz_0[j] - tg_xzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxzzzz_yyyzz_0[j] = pb_x * tg_xxzzzz_yyyzz_0[j] + fr * tg_xxzzzz_yyyzz_1[j] + fl1_fx * (tg_xzzzz_yyyzz_0[j] - tg_xzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxzzzz_yyzzz_0[j] = pb_x * tg_xxzzzz_yyzzz_0[j] + fr * tg_xxzzzz_yyzzz_1[j] + fl1_fx * (tg_xzzzz_yyzzz_0[j] - tg_xzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_yzzzz_0[j] = pb_x * tg_xxzzzz_yzzzz_0[j] + fr * tg_xxzzzz_yzzzz_1[j] + fl1_fx * (tg_xzzzz_yzzzz_0[j] - tg_xzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_zzzzz_0[j] = pb_x * tg_xxzzzz_zzzzz_0[j] + fr * tg_xxzzzz_zzzzz_1[j] + fl1_fx * (tg_xzzzz_zzzzz_0[j] - tg_xzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_xxxxx_0[j] = pb_x * tg_xyyyyy_xxxxx_0[j] + fr * tg_xyyyyy_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxx_0[j] - tg_yyyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyy_xxxx_1[j];

                    tg_xxyyyyy_xxxxy_0[j] = pb_x * tg_xyyyyy_xxxxy_0[j] + fr * tg_xyyyyy_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxy_0[j] - tg_yyyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyy_xxxy_1[j];

                    tg_xxyyyyy_xxxxz_0[j] = pb_x * tg_xyyyyy_xxxxz_0[j] + fr * tg_xyyyyy_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxz_0[j] - tg_yyyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyy_xxxz_1[j];

                    tg_xxyyyyy_xxxyy_0[j] = pb_x * tg_xyyyyy_xxxyy_0[j] + fr * tg_xyyyyy_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxyy_0[j] - tg_yyyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxyy_1[j];

                    tg_xxyyyyy_xxxyz_0[j] = pb_x * tg_xyyyyy_xxxyz_0[j] + fr * tg_xyyyyy_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxyz_0[j] - tg_yyyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxyz_1[j];

                    tg_xxyyyyy_xxxzz_0[j] = pb_x * tg_xyyyyy_xxxzz_0[j] + fr * tg_xyyyyy_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxzz_0[j] - tg_yyyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxzz_1[j];

                    tg_xxyyyyy_xxyyy_0[j] = pb_x * tg_xyyyyy_xxyyy_0[j] + fr * tg_xyyyyy_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyyy_0[j] - tg_yyyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyyy_1[j];

                    tg_xxyyyyy_xxyyz_0[j] = pb_x * tg_xyyyyy_xxyyz_0[j] + fr * tg_xyyyyy_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyyz_0[j] - tg_yyyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyyz_1[j];

                    tg_xxyyyyy_xxyzz_0[j] = pb_x * tg_xyyyyy_xxyzz_0[j] + fr * tg_xyyyyy_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyzz_0[j] - tg_yyyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyzz_1[j];

                    tg_xxyyyyy_xxzzz_0[j] = pb_x * tg_xyyyyy_xxzzz_0[j] + fr * tg_xyyyyy_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxzzz_0[j] - tg_yyyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xzzz_1[j];

                    tg_xxyyyyy_xyyyy_0[j] = pb_x * tg_xyyyyy_xyyyy_0[j] + fr * tg_xyyyyy_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyyy_0[j] - tg_yyyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyyy_1[j];

                    tg_xxyyyyy_xyyyz_0[j] = pb_x * tg_xyyyyy_xyyyz_0[j] + fr * tg_xyyyyy_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyyz_0[j] - tg_yyyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyyz_1[j];

                    tg_xxyyyyy_xyyzz_0[j] = pb_x * tg_xyyyyy_xyyzz_0[j] + fr * tg_xyyyyy_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyzz_0[j] - tg_yyyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyzz_1[j];

                    tg_xxyyyyy_xyzzz_0[j] = pb_x * tg_xyyyyy_xyzzz_0[j] + fr * tg_xyyyyy_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyzzz_0[j] - tg_yyyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yzzz_1[j];

                    tg_xxyyyyy_xzzzz_0[j] = pb_x * tg_xyyyyy_xzzzz_0[j] + fr * tg_xyyyyy_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xzzzz_0[j] - tg_yyyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_zzzz_1[j];

                    tg_xxyyyyy_yyyyy_0[j] = pb_x * tg_xyyyyy_yyyyy_0[j] + fr * tg_xyyyyy_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyyy_0[j] - tg_yyyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyyy_yyyyz_0[j] = pb_x * tg_xyyyyy_yyyyz_0[j] + fr * tg_xyyyyy_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyyz_0[j] - tg_yyyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyyy_yyyzz_0[j] = pb_x * tg_xyyyyy_yyyzz_0[j] + fr * tg_xyyyyy_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyzz_0[j] - tg_yyyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyyy_yyzzz_0[j] = pb_x * tg_xyyyyy_yyzzz_0[j] + fr * tg_xyyyyy_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyzzz_0[j] - tg_yyyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_yzzzz_0[j] = pb_x * tg_xyyyyy_yzzzz_0[j] + fr * tg_xyyyyy_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yzzzz_0[j] - tg_yyyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_zzzzz_0[j] = pb_x * tg_xyyyyy_zzzzz_0[j] + fr * tg_xyyyyy_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_zzzzz_0[j] - tg_yyyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_xxxxx_0[j] = pb_x * tg_xyyyyz_xxxxx_0[j] + fr * tg_xyyyyz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxx_0[j] - tg_yyyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyz_xxxx_1[j];

                    tg_xxyyyyz_xxxxy_0[j] = pb_x * tg_xyyyyz_xxxxy_0[j] + fr * tg_xyyyyz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxy_0[j] - tg_yyyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyz_xxxy_1[j];

                    tg_xxyyyyz_xxxxz_0[j] = pb_x * tg_xyyyyz_xxxxz_0[j] + fr * tg_xyyyyz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxz_0[j] - tg_yyyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyz_xxxz_1[j];

                    tg_xxyyyyz_xxxyy_0[j] = pb_x * tg_xyyyyz_xxxyy_0[j] + fr * tg_xyyyyz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxyy_0[j] - tg_yyyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxyy_1[j];

                    tg_xxyyyyz_xxxyz_0[j] = pb_x * tg_xyyyyz_xxxyz_0[j] + fr * tg_xyyyyz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxyz_0[j] - tg_yyyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxyz_1[j];

                    tg_xxyyyyz_xxxzz_0[j] = pb_x * tg_xyyyyz_xxxzz_0[j] + fr * tg_xyyyyz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxzz_0[j] - tg_yyyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxzz_1[j];

                    tg_xxyyyyz_xxyyy_0[j] = pb_x * tg_xyyyyz_xxyyy_0[j] + fr * tg_xyyyyz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyyy_0[j] - tg_yyyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyyy_1[j];

                    tg_xxyyyyz_xxyyz_0[j] = pb_x * tg_xyyyyz_xxyyz_0[j] + fr * tg_xyyyyz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyyz_0[j] - tg_yyyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyyz_1[j];

                    tg_xxyyyyz_xxyzz_0[j] = pb_x * tg_xyyyyz_xxyzz_0[j] + fr * tg_xyyyyz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyzz_0[j] - tg_yyyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyzz_1[j];

                    tg_xxyyyyz_xxzzz_0[j] = pb_x * tg_xyyyyz_xxzzz_0[j] + fr * tg_xyyyyz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxzzz_0[j] - tg_yyyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xzzz_1[j];

                    tg_xxyyyyz_xyyyy_0[j] = pb_x * tg_xyyyyz_xyyyy_0[j] + fr * tg_xyyyyz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyyy_0[j] - tg_yyyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyyy_1[j];

                    tg_xxyyyyz_xyyyz_0[j] = pb_x * tg_xyyyyz_xyyyz_0[j] + fr * tg_xyyyyz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyyz_0[j] - tg_yyyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyyz_1[j];

                    tg_xxyyyyz_xyyzz_0[j] = pb_x * tg_xyyyyz_xyyzz_0[j] + fr * tg_xyyyyz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyzz_0[j] - tg_yyyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyzz_1[j];

                    tg_xxyyyyz_xyzzz_0[j] = pb_x * tg_xyyyyz_xyzzz_0[j] + fr * tg_xyyyyz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyzzz_0[j] - tg_yyyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yzzz_1[j];

                    tg_xxyyyyz_xzzzz_0[j] = pb_x * tg_xyyyyz_xzzzz_0[j] + fr * tg_xyyyyz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xzzzz_0[j] - tg_yyyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_zzzz_1[j];

                    tg_xxyyyyz_yyyyy_0[j] = pb_x * tg_xyyyyz_yyyyy_0[j] + fr * tg_xyyyyz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyyy_0[j] - tg_yyyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyyz_yyyyz_0[j] = pb_x * tg_xyyyyz_yyyyz_0[j] + fr * tg_xyyyyz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyyz_0[j] - tg_yyyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyyz_yyyzz_0[j] = pb_x * tg_xyyyyz_yyyzz_0[j] + fr * tg_xyyyyz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyzz_0[j] - tg_yyyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyyz_yyzzz_0[j] = pb_x * tg_xyyyyz_yyzzz_0[j] + fr * tg_xyyyyz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyzzz_0[j] - tg_yyyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_yzzzz_0[j] = pb_x * tg_xyyyyz_yzzzz_0[j] + fr * tg_xyyyyz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yzzzz_0[j] - tg_yyyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_zzzzz_0[j] = pb_x * tg_xyyyyz_zzzzz_0[j] + fr * tg_xyyyyz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_zzzzz_0[j] - tg_yyyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_xxxxx_0[j] = pb_x * tg_xyyyzz_xxxxx_0[j] + fr * tg_xyyyzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxx_0[j] - tg_yyyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyzz_xxxx_1[j];

                    tg_xxyyyzz_xxxxy_0[j] = pb_x * tg_xyyyzz_xxxxy_0[j] + fr * tg_xyyyzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxy_0[j] - tg_yyyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzz_xxxy_1[j];

                    tg_xxyyyzz_xxxxz_0[j] = pb_x * tg_xyyyzz_xxxxz_0[j] + fr * tg_xyyyzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxz_0[j] - tg_yyyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzz_xxxz_1[j];

                    tg_xxyyyzz_xxxyy_0[j] = pb_x * tg_xyyyzz_xxxyy_0[j] + fr * tg_xyyyzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxyy_0[j] - tg_yyyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxyy_1[j];

                    tg_xxyyyzz_xxxyz_0[j] = pb_x * tg_xyyyzz_xxxyz_0[j] + fr * tg_xyyyzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxyz_0[j] - tg_yyyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxyz_1[j];

                    tg_xxyyyzz_xxxzz_0[j] = pb_x * tg_xyyyzz_xxxzz_0[j] + fr * tg_xyyyzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxzz_0[j] - tg_yyyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxzz_1[j];

                    tg_xxyyyzz_xxyyy_0[j] = pb_x * tg_xyyyzz_xxyyy_0[j] + fr * tg_xyyyzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyyy_0[j] - tg_yyyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyyy_1[j];

                    tg_xxyyyzz_xxyyz_0[j] = pb_x * tg_xyyyzz_xxyyz_0[j] + fr * tg_xyyyzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyyz_0[j] - tg_yyyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyyz_1[j];

                    tg_xxyyyzz_xxyzz_0[j] = pb_x * tg_xyyyzz_xxyzz_0[j] + fr * tg_xyyyzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyzz_0[j] - tg_yyyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyzz_1[j];

                    tg_xxyyyzz_xxzzz_0[j] = pb_x * tg_xyyyzz_xxzzz_0[j] + fr * tg_xyyyzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxzzz_0[j] - tg_yyyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xzzz_1[j];

                    tg_xxyyyzz_xyyyy_0[j] = pb_x * tg_xyyyzz_xyyyy_0[j] + fr * tg_xyyyzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyyy_0[j] - tg_yyyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyyy_1[j];

                    tg_xxyyyzz_xyyyz_0[j] = pb_x * tg_xyyyzz_xyyyz_0[j] + fr * tg_xyyyzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyyz_0[j] - tg_yyyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyyz_1[j];

                    tg_xxyyyzz_xyyzz_0[j] = pb_x * tg_xyyyzz_xyyzz_0[j] + fr * tg_xyyyzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyzz_0[j] - tg_yyyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyzz_1[j];

                    tg_xxyyyzz_xyzzz_0[j] = pb_x * tg_xyyyzz_xyzzz_0[j] + fr * tg_xyyyzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyzzz_0[j] - tg_yyyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yzzz_1[j];

                    tg_xxyyyzz_xzzzz_0[j] = pb_x * tg_xyyyzz_xzzzz_0[j] + fr * tg_xyyyzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xzzzz_0[j] - tg_yyyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_zzzz_1[j];

                    tg_xxyyyzz_yyyyy_0[j] = pb_x * tg_xyyyzz_yyyyy_0[j] + fr * tg_xyyyzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyyy_0[j] - tg_yyyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyzz_yyyyz_0[j] = pb_x * tg_xyyyzz_yyyyz_0[j] + fr * tg_xyyyzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyyz_0[j] - tg_yyyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyzz_yyyzz_0[j] = pb_x * tg_xyyyzz_yyyzz_0[j] + fr * tg_xyyyzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyzz_0[j] - tg_yyyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyzz_yyzzz_0[j] = pb_x * tg_xyyyzz_yyzzz_0[j] + fr * tg_xyyyzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyzzz_0[j] - tg_yyyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_yzzzz_0[j] = pb_x * tg_xyyyzz_yzzzz_0[j] + fr * tg_xyyyzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yzzzz_0[j] - tg_yyyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_zzzzz_0[j] = pb_x * tg_xyyyzz_zzzzz_0[j] + fr * tg_xyyyzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_zzzzz_0[j] - tg_yyyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_xxxxx_0[j] = pb_x * tg_xyyzzz_xxxxx_0[j] + fr * tg_xyyzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxx_0[j] - tg_yyzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyzzz_xxxx_1[j];

                    tg_xxyyzzz_xxxxy_0[j] = pb_x * tg_xyyzzz_xxxxy_0[j] + fr * tg_xyyzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxy_0[j] - tg_yyzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzz_xxxy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSH_380_474(      CMemBlock2D<double>* primBuffer,
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

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_xyyzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 380); 

                auto tg_xyyzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 381); 

                auto tg_xyyzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 382); 

                auto tg_xyyzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 383); 

                auto tg_xyyzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 384); 

                auto tg_xyyzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 385); 

                auto tg_xyyzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 386); 

                auto tg_xyyzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 387); 

                auto tg_xyyzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 388); 

                auto tg_xyyzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 389); 

                auto tg_xyyzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 390); 

                auto tg_xyyzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 391); 

                auto tg_xyyzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 392); 

                auto tg_xyyzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 393); 

                auto tg_xyyzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 394); 

                auto tg_xyyzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 395); 

                auto tg_xyyzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 396); 

                auto tg_xyyzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 397); 

                auto tg_xyyzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 398); 

                auto tg_xyzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 399); 

                auto tg_xyzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 400); 

                auto tg_xyzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 401); 

                auto tg_xyzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 402); 

                auto tg_xyzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 403); 

                auto tg_xyzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 404); 

                auto tg_xyzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 405); 

                auto tg_xyzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 406); 

                auto tg_xyzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 407); 

                auto tg_xyzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 408); 

                auto tg_xyzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 409); 

                auto tg_xyzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 410); 

                auto tg_xyzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 411); 

                auto tg_xyzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 412); 

                auto tg_xyzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 413); 

                auto tg_xyzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 414); 

                auto tg_xyzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 415); 

                auto tg_xyzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 416); 

                auto tg_xyzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 417); 

                auto tg_xyzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 418); 

                auto tg_xyzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 419); 

                auto tg_xzzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 420); 

                auto tg_xzzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 421); 

                auto tg_xzzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 422); 

                auto tg_xzzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 423); 

                auto tg_xzzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 424); 

                auto tg_xzzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 425); 

                auto tg_xzzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 426); 

                auto tg_xzzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 427); 

                auto tg_xzzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 428); 

                auto tg_xzzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 429); 

                auto tg_xzzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 430); 

                auto tg_xzzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 431); 

                auto tg_xzzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 432); 

                auto tg_xzzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 433); 

                auto tg_xzzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 434); 

                auto tg_xzzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 435); 

                auto tg_xzzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 436); 

                auto tg_xzzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 437); 

                auto tg_xzzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 438); 

                auto tg_xzzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 439); 

                auto tg_xzzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 440); 

                auto tg_yyyyyy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 441); 

                auto tg_yyyyyy_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 442); 

                auto tg_yyyyyy_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 443); 

                auto tg_yyyyyy_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 444); 

                auto tg_yyyyyy_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 445); 

                auto tg_yyyyyy_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 446); 

                auto tg_yyyyyy_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 447); 

                auto tg_yyyyyy_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 448); 

                auto tg_yyyyyy_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 449); 

                auto tg_yyyyyy_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 450); 

                auto tg_yyyyyy_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 451); 

                auto tg_yyyyyy_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 452); 

                auto tg_yyyyyy_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 453); 

                auto tg_yyyyyy_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 454); 

                auto tg_yyyyyy_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 455); 

                auto tg_yyyyyy_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 456); 

                auto tg_yyyyyy_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 457); 

                auto tg_yyyyyy_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 458); 

                auto tg_yyyyyy_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 459); 

                auto tg_yyyyyy_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 460); 

                auto tg_yyyyyy_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 461); 

                auto tg_yyyyyz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 462); 

                auto tg_yyyyyz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 463); 

                auto tg_yyyyyz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 464); 

                auto tg_yyyyyz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 465); 

                auto tg_yyyyyz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 466); 

                auto tg_yyyyyz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 467); 

                auto tg_yyyyyz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 468); 

                auto tg_yyyyyz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 469); 

                auto tg_yyyyyz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 470); 

                auto tg_yyyyyz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 471); 

                auto tg_yyyyyz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 472); 

                auto tg_yyyyyz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 473); 

                auto tg_xyyzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 380); 

                auto tg_xyyzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 381); 

                auto tg_xyyzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 382); 

                auto tg_xyyzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 383); 

                auto tg_xyyzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 384); 

                auto tg_xyyzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 385); 

                auto tg_xyyzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 386); 

                auto tg_xyyzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 387); 

                auto tg_xyyzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 388); 

                auto tg_xyyzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 389); 

                auto tg_xyyzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 390); 

                auto tg_xyyzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 391); 

                auto tg_xyyzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 392); 

                auto tg_xyyzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 393); 

                auto tg_xyyzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 394); 

                auto tg_xyyzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 395); 

                auto tg_xyyzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 396); 

                auto tg_xyyzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 397); 

                auto tg_xyyzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 398); 

                auto tg_xyzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 399); 

                auto tg_xyzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 400); 

                auto tg_xyzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 401); 

                auto tg_xyzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 402); 

                auto tg_xyzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 403); 

                auto tg_xyzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 404); 

                auto tg_xyzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 405); 

                auto tg_xyzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 406); 

                auto tg_xyzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 407); 

                auto tg_xyzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 408); 

                auto tg_xyzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 409); 

                auto tg_xyzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 410); 

                auto tg_xyzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 411); 

                auto tg_xyzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 412); 

                auto tg_xyzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 413); 

                auto tg_xyzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 414); 

                auto tg_xyzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 415); 

                auto tg_xyzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 416); 

                auto tg_xyzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 417); 

                auto tg_xyzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 418); 

                auto tg_xyzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 419); 

                auto tg_xzzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 420); 

                auto tg_xzzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 421); 

                auto tg_xzzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 422); 

                auto tg_xzzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 423); 

                auto tg_xzzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 424); 

                auto tg_xzzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 425); 

                auto tg_xzzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 426); 

                auto tg_xzzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 427); 

                auto tg_xzzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 428); 

                auto tg_xzzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 429); 

                auto tg_xzzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 430); 

                auto tg_xzzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 431); 

                auto tg_xzzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 432); 

                auto tg_xzzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 433); 

                auto tg_xzzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 434); 

                auto tg_xzzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 435); 

                auto tg_xzzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 436); 

                auto tg_xzzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 437); 

                auto tg_xzzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 438); 

                auto tg_xzzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 439); 

                auto tg_xzzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 440); 

                auto tg_yyyyyy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 441); 

                auto tg_yyyyyy_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 442); 

                auto tg_yyyyyy_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 443); 

                auto tg_yyyyyy_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 444); 

                auto tg_yyyyyy_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 445); 

                auto tg_yyyyyy_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 446); 

                auto tg_yyyyyy_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 447); 

                auto tg_yyyyyy_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 448); 

                auto tg_yyyyyy_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 449); 

                auto tg_yyyyyy_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 450); 

                auto tg_yyyyyy_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 451); 

                auto tg_yyyyyy_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 452); 

                auto tg_yyyyyy_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 453); 

                auto tg_yyyyyy_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 454); 

                auto tg_yyyyyy_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 455); 

                auto tg_yyyyyy_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 456); 

                auto tg_yyyyyy_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 457); 

                auto tg_yyyyyy_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 458); 

                auto tg_yyyyyy_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 459); 

                auto tg_yyyyyy_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 460); 

                auto tg_yyyyyy_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 461); 

                auto tg_yyyyyz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 462); 

                auto tg_yyyyyz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 463); 

                auto tg_yyyyyz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 464); 

                auto tg_yyyyyz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 465); 

                auto tg_yyyyyz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 466); 

                auto tg_yyyyyz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 467); 

                auto tg_yyyyyz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 468); 

                auto tg_yyyyyz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 469); 

                auto tg_yyyyyz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 470); 

                auto tg_yyyyyz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 471); 

                auto tg_yyyyyz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 472); 

                auto tg_yyyyyz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 473); 

                auto tg_yyzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 380); 

                auto tg_yyzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 381); 

                auto tg_yyzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 382); 

                auto tg_yyzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 383); 

                auto tg_yyzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 384); 

                auto tg_yyzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 385); 

                auto tg_yyzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 386); 

                auto tg_yyzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 387); 

                auto tg_yyzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 388); 

                auto tg_yyzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 389); 

                auto tg_yyzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 390); 

                auto tg_yyzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 391); 

                auto tg_yyzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 392); 

                auto tg_yyzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 393); 

                auto tg_yyzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 394); 

                auto tg_yyzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 395); 

                auto tg_yyzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 396); 

                auto tg_yyzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 397); 

                auto tg_yyzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 398); 

                auto tg_yzzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 399); 

                auto tg_yzzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 400); 

                auto tg_yzzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 401); 

                auto tg_yzzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 402); 

                auto tg_yzzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 403); 

                auto tg_yzzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 404); 

                auto tg_yzzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 405); 

                auto tg_yzzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 406); 

                auto tg_yzzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 407); 

                auto tg_yzzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 408); 

                auto tg_yzzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 409); 

                auto tg_yzzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 410); 

                auto tg_yzzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 411); 

                auto tg_yzzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 412); 

                auto tg_yzzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 413); 

                auto tg_yzzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 414); 

                auto tg_yzzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 415); 

                auto tg_yzzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 416); 

                auto tg_yzzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 417); 

                auto tg_yzzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 418); 

                auto tg_yzzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 419); 

                auto tg_zzzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 420); 

                auto tg_zzzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 421); 

                auto tg_zzzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 422); 

                auto tg_zzzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 423); 

                auto tg_zzzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 424); 

                auto tg_zzzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 425); 

                auto tg_zzzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 426); 

                auto tg_zzzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 427); 

                auto tg_zzzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 428); 

                auto tg_zzzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 429); 

                auto tg_zzzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 430); 

                auto tg_zzzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 431); 

                auto tg_zzzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 432); 

                auto tg_zzzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 433); 

                auto tg_zzzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 434); 

                auto tg_zzzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 435); 

                auto tg_zzzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 436); 

                auto tg_zzzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 437); 

                auto tg_zzzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 438); 

                auto tg_zzzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 439); 

                auto tg_zzzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 440); 

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

                auto tg_xyyzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 272); 

                auto tg_xyyzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 273); 

                auto tg_xyyzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 274); 

                auto tg_xyyzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 275); 

                auto tg_xyyzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 276); 

                auto tg_xyyzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 277); 

                auto tg_xyyzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 278); 

                auto tg_xyyzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 279); 

                auto tg_xyyzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 280); 

                auto tg_xyyzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 281); 

                auto tg_xyyzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 282); 

                auto tg_xyyzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 283); 

                auto tg_xyyzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 284); 

                auto tg_xyzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 285); 

                auto tg_xyzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 286); 

                auto tg_xyzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 287); 

                auto tg_xyzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 288); 

                auto tg_xyzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 289); 

                auto tg_xyzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 290); 

                auto tg_xyzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 291); 

                auto tg_xyzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 292); 

                auto tg_xyzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 293); 

                auto tg_xyzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 294); 

                auto tg_xyzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 295); 

                auto tg_xyzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 296); 

                auto tg_xyzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 297); 

                auto tg_xyzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 298); 

                auto tg_xyzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 299); 

                auto tg_xzzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 300); 

                auto tg_xzzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 301); 

                auto tg_xzzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 302); 

                auto tg_xzzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 303); 

                auto tg_xzzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 304); 

                auto tg_xzzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 305); 

                auto tg_xzzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 306); 

                auto tg_xzzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 307); 

                auto tg_xzzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 308); 

                auto tg_xzzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 309); 

                auto tg_xzzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 310); 

                auto tg_xzzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 311); 

                auto tg_xzzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 312); 

                auto tg_xzzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 313); 

                auto tg_xzzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 314); 

                auto tg_yyyyyy_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 315); 

                auto tg_yyyyyy_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 316); 

                auto tg_yyyyyy_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 317); 

                auto tg_yyyyyy_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 318); 

                auto tg_yyyyyy_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 319); 

                auto tg_yyyyyy_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 320); 

                auto tg_yyyyyy_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 321); 

                auto tg_yyyyyy_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 322); 

                auto tg_yyyyyy_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 323); 

                auto tg_yyyyyy_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 324); 

                auto tg_yyyyyy_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 325); 

                auto tg_yyyyyy_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 326); 

                auto tg_yyyyyy_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 327); 

                auto tg_yyyyyy_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 328); 

                auto tg_yyyyyy_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 329); 

                auto tg_yyyyyz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 330); 

                auto tg_yyyyyz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 331); 

                auto tg_yyyyyz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 332); 

                auto tg_yyyyyz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 333); 

                auto tg_yyyyyz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 334); 

                auto tg_yyyyyz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 335); 

                auto tg_yyyyyz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 336); 

                auto tg_yyyyyz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 337); 

                auto tg_yyyyyz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 338); 

                auto tg_yyyyyz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 339); 

                auto tg_yyyyyz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 340); 

                auto tg_yyyyyz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 341); 

                // set up pointers to integrals

                auto tg_xxyyzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 380); 

                auto tg_xxyyzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 381); 

                auto tg_xxyyzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 382); 

                auto tg_xxyyzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 383); 

                auto tg_xxyyzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 384); 

                auto tg_xxyyzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 385); 

                auto tg_xxyyzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 386); 

                auto tg_xxyyzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 387); 

                auto tg_xxyyzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 388); 

                auto tg_xxyyzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 389); 

                auto tg_xxyyzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 390); 

                auto tg_xxyyzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 391); 

                auto tg_xxyyzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 392); 

                auto tg_xxyyzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 393); 

                auto tg_xxyyzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 394); 

                auto tg_xxyyzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 395); 

                auto tg_xxyyzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 396); 

                auto tg_xxyyzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 397); 

                auto tg_xxyyzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 398); 

                auto tg_xxyzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 399); 

                auto tg_xxyzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 400); 

                auto tg_xxyzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 401); 

                auto tg_xxyzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 402); 

                auto tg_xxyzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 403); 

                auto tg_xxyzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 404); 

                auto tg_xxyzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 405); 

                auto tg_xxyzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 406); 

                auto tg_xxyzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 407); 

                auto tg_xxyzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 408); 

                auto tg_xxyzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 409); 

                auto tg_xxyzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 410); 

                auto tg_xxyzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 411); 

                auto tg_xxyzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 412); 

                auto tg_xxyzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 413); 

                auto tg_xxyzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 414); 

                auto tg_xxyzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 415); 

                auto tg_xxyzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 416); 

                auto tg_xxyzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 417); 

                auto tg_xxyzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 418); 

                auto tg_xxyzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 419); 

                auto tg_xxzzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 420); 

                auto tg_xxzzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 421); 

                auto tg_xxzzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 422); 

                auto tg_xxzzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 423); 

                auto tg_xxzzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 424); 

                auto tg_xxzzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 425); 

                auto tg_xxzzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 426); 

                auto tg_xxzzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 427); 

                auto tg_xxzzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 428); 

                auto tg_xxzzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 429); 

                auto tg_xxzzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 430); 

                auto tg_xxzzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 431); 

                auto tg_xxzzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 432); 

                auto tg_xxzzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 433); 

                auto tg_xxzzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 434); 

                auto tg_xxzzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 435); 

                auto tg_xxzzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 436); 

                auto tg_xxzzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 437); 

                auto tg_xxzzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 438); 

                auto tg_xxzzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 439); 

                auto tg_xxzzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 440); 

                auto tg_xyyyyyy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 441); 

                auto tg_xyyyyyy_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 442); 

                auto tg_xyyyyyy_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 443); 

                auto tg_xyyyyyy_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 444); 

                auto tg_xyyyyyy_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 445); 

                auto tg_xyyyyyy_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 446); 

                auto tg_xyyyyyy_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 447); 

                auto tg_xyyyyyy_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 448); 

                auto tg_xyyyyyy_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 449); 

                auto tg_xyyyyyy_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 450); 

                auto tg_xyyyyyy_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 451); 

                auto tg_xyyyyyy_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 452); 

                auto tg_xyyyyyy_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 453); 

                auto tg_xyyyyyy_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 454); 

                auto tg_xyyyyyy_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 455); 

                auto tg_xyyyyyy_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 456); 

                auto tg_xyyyyyy_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 457); 

                auto tg_xyyyyyy_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 458); 

                auto tg_xyyyyyy_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 459); 

                auto tg_xyyyyyy_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 460); 

                auto tg_xyyyyyy_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 461); 

                auto tg_xyyyyyz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 462); 

                auto tg_xyyyyyz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 463); 

                auto tg_xyyyyyz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 464); 

                auto tg_xyyyyyz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 465); 

                auto tg_xyyyyyz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 466); 

                auto tg_xyyyyyz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 467); 

                auto tg_xyyyyyz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 468); 

                auto tg_xyyyyyz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 469); 

                auto tg_xyyyyyz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 470); 

                auto tg_xyyyyyz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 471); 

                auto tg_xyyyyyz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 472); 

                auto tg_xyyyyyz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 473); 

                // Batch of Integrals (380,474)

                #pragma omp simd aligned(fxn, fza, tg_xxyyzzz_xxxxz_0, tg_xxyyzzz_xxxyy_0, tg_xxyyzzz_xxxyz_0, \
                                         tg_xxyyzzz_xxxzz_0, tg_xxyyzzz_xxyyy_0, tg_xxyyzzz_xxyyz_0, tg_xxyyzzz_xxyzz_0, \
                                         tg_xxyyzzz_xxzzz_0, tg_xxyyzzz_xyyyy_0, tg_xxyyzzz_xyyyz_0, tg_xxyyzzz_xyyzz_0, \
                                         tg_xxyyzzz_xyzzz_0, tg_xxyyzzz_xzzzz_0, tg_xxyyzzz_yyyyy_0, tg_xxyyzzz_yyyyz_0, \
                                         tg_xxyyzzz_yyyzz_0, tg_xxyyzzz_yyzzz_0, tg_xxyyzzz_yzzzz_0, tg_xxyyzzz_zzzzz_0, \
                                         tg_xxyzzzz_xxxxx_0, tg_xxyzzzz_xxxxy_0, tg_xxyzzzz_xxxxz_0, tg_xxyzzzz_xxxyy_0, \
                                         tg_xxyzzzz_xxxyz_0, tg_xxyzzzz_xxxzz_0, tg_xxyzzzz_xxyyy_0, tg_xxyzzzz_xxyyz_0, \
                                         tg_xxyzzzz_xxyzz_0, tg_xxyzzzz_xxzzz_0, tg_xxyzzzz_xyyyy_0, tg_xxyzzzz_xyyyz_0, \
                                         tg_xxyzzzz_xyyzz_0, tg_xxyzzzz_xyzzz_0, tg_xxyzzzz_xzzzz_0, tg_xxyzzzz_yyyyy_0, \
                                         tg_xxyzzzz_yyyyz_0, tg_xxyzzzz_yyyzz_0, tg_xxyzzzz_yyzzz_0, tg_xxyzzzz_yzzzz_0, \
                                         tg_xxyzzzz_zzzzz_0, tg_xxzzzzz_xxxxx_0, tg_xxzzzzz_xxxxy_0, tg_xxzzzzz_xxxxz_0, \
                                         tg_xxzzzzz_xxxyy_0, tg_xxzzzzz_xxxyz_0, tg_xxzzzzz_xxxzz_0, tg_xxzzzzz_xxyyy_0, \
                                         tg_xxzzzzz_xxyyz_0, tg_xxzzzzz_xxyzz_0, tg_xxzzzzz_xxzzz_0, tg_xxzzzzz_xyyyy_0, \
                                         tg_xxzzzzz_xyyyz_0, tg_xxzzzzz_xyyzz_0, tg_xxzzzzz_xyzzz_0, tg_xxzzzzz_xzzzz_0, \
                                         tg_xxzzzzz_yyyyy_0, tg_xxzzzzz_yyyyz_0, tg_xxzzzzz_yyyzz_0, tg_xxzzzzz_yyzzz_0, \
                                         tg_xxzzzzz_yzzzz_0, tg_xxzzzzz_zzzzz_0, tg_xyyyyyy_xxxxx_0, tg_xyyyyyy_xxxxy_0, \
                                         tg_xyyyyyy_xxxxz_0, tg_xyyyyyy_xxxyy_0, tg_xyyyyyy_xxxyz_0, tg_xyyyyyy_xxxzz_0, \
                                         tg_xyyyyyy_xxyyy_0, tg_xyyyyyy_xxyyz_0, tg_xyyyyyy_xxyzz_0, tg_xyyyyyy_xxzzz_0, \
                                         tg_xyyyyyy_xyyyy_0, tg_xyyyyyy_xyyyz_0, tg_xyyyyyy_xyyzz_0, tg_xyyyyyy_xyzzz_0, \
                                         tg_xyyyyyy_xzzzz_0, tg_xyyyyyy_yyyyy_0, tg_xyyyyyy_yyyyz_0, tg_xyyyyyy_yyyzz_0, \
                                         tg_xyyyyyy_yyzzz_0, tg_xyyyyyy_yzzzz_0, tg_xyyyyyy_zzzzz_0, tg_xyyyyyz_xxxxx_0, \
                                         tg_xyyyyyz_xxxxy_0, tg_xyyyyyz_xxxxz_0, tg_xyyyyyz_xxxyy_0, tg_xyyyyyz_xxxyz_0, \
                                         tg_xyyyyyz_xxxzz_0, tg_xyyyyyz_xxyyy_0, tg_xyyyyyz_xxyyz_0, tg_xyyyyyz_xxyzz_0, \
                                         tg_xyyyyyz_xxzzz_0, tg_xyyyyyz_xyyyy_0, tg_xyyyyyz_xyyyz_0, tg_xyyzzz_xxxxz_0, \
                                         tg_xyyzzz_xxxxz_1, tg_xyyzzz_xxxyy_0, tg_xyyzzz_xxxyy_1, tg_xyyzzz_xxxyz_0, \
                                         tg_xyyzzz_xxxyz_1, tg_xyyzzz_xxxz_1, tg_xyyzzz_xxxzz_0, tg_xyyzzz_xxxzz_1, \
                                         tg_xyyzzz_xxyy_1, tg_xyyzzz_xxyyy_0, tg_xyyzzz_xxyyy_1, tg_xyyzzz_xxyyz_0, \
                                         tg_xyyzzz_xxyyz_1, tg_xyyzzz_xxyz_1, tg_xyyzzz_xxyzz_0, tg_xyyzzz_xxyzz_1, \
                                         tg_xyyzzz_xxzz_1, tg_xyyzzz_xxzzz_0, tg_xyyzzz_xxzzz_1, tg_xyyzzz_xyyy_1, \
                                         tg_xyyzzz_xyyyy_0, tg_xyyzzz_xyyyy_1, tg_xyyzzz_xyyyz_0, tg_xyyzzz_xyyyz_1, \
                                         tg_xyyzzz_xyyz_1, tg_xyyzzz_xyyzz_0, tg_xyyzzz_xyyzz_1, tg_xyyzzz_xyzz_1, \
                                         tg_xyyzzz_xyzzz_0, tg_xyyzzz_xyzzz_1, tg_xyyzzz_xzzz_1, tg_xyyzzz_xzzzz_0, \
                                         tg_xyyzzz_xzzzz_1, tg_xyyzzz_yyyy_1, tg_xyyzzz_yyyyy_0, tg_xyyzzz_yyyyy_1, \
                                         tg_xyyzzz_yyyyz_0, tg_xyyzzz_yyyyz_1, tg_xyyzzz_yyyz_1, tg_xyyzzz_yyyzz_0, \
                                         tg_xyyzzz_yyyzz_1, tg_xyyzzz_yyzz_1, tg_xyyzzz_yyzzz_0, tg_xyyzzz_yyzzz_1, \
                                         tg_xyyzzz_yzzz_1, tg_xyyzzz_yzzzz_0, tg_xyyzzz_yzzzz_1, tg_xyyzzz_zzzz_1, \
                                         tg_xyyzzz_zzzzz_0, tg_xyyzzz_zzzzz_1, tg_xyzzzz_xxxx_1, tg_xyzzzz_xxxxx_0, \
                                         tg_xyzzzz_xxxxx_1, tg_xyzzzz_xxxxy_0, tg_xyzzzz_xxxxy_1, tg_xyzzzz_xxxxz_0, \
                                         tg_xyzzzz_xxxxz_1, tg_xyzzzz_xxxy_1, tg_xyzzzz_xxxyy_0, tg_xyzzzz_xxxyy_1, \
                                         tg_xyzzzz_xxxyz_0, tg_xyzzzz_xxxyz_1, tg_xyzzzz_xxxz_1, tg_xyzzzz_xxxzz_0, \
                                         tg_xyzzzz_xxxzz_1, tg_xyzzzz_xxyy_1, tg_xyzzzz_xxyyy_0, tg_xyzzzz_xxyyy_1, \
                                         tg_xyzzzz_xxyyz_0, tg_xyzzzz_xxyyz_1, tg_xyzzzz_xxyz_1, tg_xyzzzz_xxyzz_0, \
                                         tg_xyzzzz_xxyzz_1, tg_xyzzzz_xxzz_1, tg_xyzzzz_xxzzz_0, tg_xyzzzz_xxzzz_1, \
                                         tg_xyzzzz_xyyy_1, tg_xyzzzz_xyyyy_0, tg_xyzzzz_xyyyy_1, tg_xyzzzz_xyyyz_0, \
                                         tg_xyzzzz_xyyyz_1, tg_xyzzzz_xyyz_1, tg_xyzzzz_xyyzz_0, tg_xyzzzz_xyyzz_1, \
                                         tg_xyzzzz_xyzz_1, tg_xyzzzz_xyzzz_0, tg_xyzzzz_xyzzz_1, tg_xyzzzz_xzzz_1, \
                                         tg_xyzzzz_xzzzz_0, tg_xyzzzz_xzzzz_1, tg_xyzzzz_yyyy_1, tg_xyzzzz_yyyyy_0, \
                                         tg_xyzzzz_yyyyy_1, tg_xyzzzz_yyyyz_0, tg_xyzzzz_yyyyz_1, tg_xyzzzz_yyyz_1, \
                                         tg_xyzzzz_yyyzz_0, tg_xyzzzz_yyyzz_1, tg_xyzzzz_yyzz_1, tg_xyzzzz_yyzzz_0, \
                                         tg_xyzzzz_yyzzz_1, tg_xyzzzz_yzzz_1, tg_xyzzzz_yzzzz_0, tg_xyzzzz_yzzzz_1, \
                                         tg_xyzzzz_zzzz_1, tg_xyzzzz_zzzzz_0, tg_xyzzzz_zzzzz_1, tg_xzzzzz_xxxx_1, \
                                         tg_xzzzzz_xxxxx_0, tg_xzzzzz_xxxxx_1, tg_xzzzzz_xxxxy_0, tg_xzzzzz_xxxxy_1, \
                                         tg_xzzzzz_xxxxz_0, tg_xzzzzz_xxxxz_1, tg_xzzzzz_xxxy_1, tg_xzzzzz_xxxyy_0, \
                                         tg_xzzzzz_xxxyy_1, tg_xzzzzz_xxxyz_0, tg_xzzzzz_xxxyz_1, tg_xzzzzz_xxxz_1, \
                                         tg_xzzzzz_xxxzz_0, tg_xzzzzz_xxxzz_1, tg_xzzzzz_xxyy_1, tg_xzzzzz_xxyyy_0, \
                                         tg_xzzzzz_xxyyy_1, tg_xzzzzz_xxyyz_0, tg_xzzzzz_xxyyz_1, tg_xzzzzz_xxyz_1, \
                                         tg_xzzzzz_xxyzz_0, tg_xzzzzz_xxyzz_1, tg_xzzzzz_xxzz_1, tg_xzzzzz_xxzzz_0, \
                                         tg_xzzzzz_xxzzz_1, tg_xzzzzz_xyyy_1, tg_xzzzzz_xyyyy_0, tg_xzzzzz_xyyyy_1, \
                                         tg_xzzzzz_xyyyz_0, tg_xzzzzz_xyyyz_1, tg_xzzzzz_xyyz_1, tg_xzzzzz_xyyzz_0, \
                                         tg_xzzzzz_xyyzz_1, tg_xzzzzz_xyzz_1, tg_xzzzzz_xyzzz_0, tg_xzzzzz_xyzzz_1, \
                                         tg_xzzzzz_xzzz_1, tg_xzzzzz_xzzzz_0, tg_xzzzzz_xzzzz_1, tg_xzzzzz_yyyy_1, \
                                         tg_xzzzzz_yyyyy_0, tg_xzzzzz_yyyyy_1, tg_xzzzzz_yyyyz_0, tg_xzzzzz_yyyyz_1, \
                                         tg_xzzzzz_yyyz_1, tg_xzzzzz_yyyzz_0, tg_xzzzzz_yyyzz_1, tg_xzzzzz_yyzz_1, \
                                         tg_xzzzzz_yyzzz_0, tg_xzzzzz_yyzzz_1, tg_xzzzzz_yzzz_1, tg_xzzzzz_yzzzz_0, \
                                         tg_xzzzzz_yzzzz_1, tg_xzzzzz_zzzz_1, tg_xzzzzz_zzzzz_0, tg_xzzzzz_zzzzz_1, \
                                         tg_yyyyyy_xxxx_1, tg_yyyyyy_xxxxx_0, tg_yyyyyy_xxxxx_1, tg_yyyyyy_xxxxy_0, \
                                         tg_yyyyyy_xxxxy_1, tg_yyyyyy_xxxxz_0, tg_yyyyyy_xxxxz_1, tg_yyyyyy_xxxy_1, \
                                         tg_yyyyyy_xxxyy_0, tg_yyyyyy_xxxyy_1, tg_yyyyyy_xxxyz_0, tg_yyyyyy_xxxyz_1, \
                                         tg_yyyyyy_xxxz_1, tg_yyyyyy_xxxzz_0, tg_yyyyyy_xxxzz_1, tg_yyyyyy_xxyy_1, \
                                         tg_yyyyyy_xxyyy_0, tg_yyyyyy_xxyyy_1, tg_yyyyyy_xxyyz_0, tg_yyyyyy_xxyyz_1, \
                                         tg_yyyyyy_xxyz_1, tg_yyyyyy_xxyzz_0, tg_yyyyyy_xxyzz_1, tg_yyyyyy_xxzz_1, \
                                         tg_yyyyyy_xxzzz_0, tg_yyyyyy_xxzzz_1, tg_yyyyyy_xyyy_1, tg_yyyyyy_xyyyy_0, \
                                         tg_yyyyyy_xyyyy_1, tg_yyyyyy_xyyyz_0, tg_yyyyyy_xyyyz_1, tg_yyyyyy_xyyz_1, \
                                         tg_yyyyyy_xyyzz_0, tg_yyyyyy_xyyzz_1, tg_yyyyyy_xyzz_1, tg_yyyyyy_xyzzz_0, \
                                         tg_yyyyyy_xyzzz_1, tg_yyyyyy_xzzz_1, tg_yyyyyy_xzzzz_0, tg_yyyyyy_xzzzz_1, \
                                         tg_yyyyyy_yyyy_1, tg_yyyyyy_yyyyy_0, tg_yyyyyy_yyyyy_1, tg_yyyyyy_yyyyz_0, \
                                         tg_yyyyyy_yyyyz_1, tg_yyyyyy_yyyz_1, tg_yyyyyy_yyyzz_0, tg_yyyyyy_yyyzz_1, \
                                         tg_yyyyyy_yyzz_1, tg_yyyyyy_yyzzz_0, tg_yyyyyy_yyzzz_1, tg_yyyyyy_yzzz_1, \
                                         tg_yyyyyy_yzzzz_0, tg_yyyyyy_yzzzz_1, tg_yyyyyy_zzzz_1, tg_yyyyyy_zzzzz_0, \
                                         tg_yyyyyy_zzzzz_1, tg_yyyyyz_xxxx_1, tg_yyyyyz_xxxxx_0, tg_yyyyyz_xxxxx_1, \
                                         tg_yyyyyz_xxxxy_0, tg_yyyyyz_xxxxy_1, tg_yyyyyz_xxxxz_0, tg_yyyyyz_xxxxz_1, \
                                         tg_yyyyyz_xxxy_1, tg_yyyyyz_xxxyy_0, tg_yyyyyz_xxxyy_1, tg_yyyyyz_xxxyz_0, \
                                         tg_yyyyyz_xxxyz_1, tg_yyyyyz_xxxz_1, tg_yyyyyz_xxxzz_0, tg_yyyyyz_xxxzz_1, \
                                         tg_yyyyyz_xxyy_1, tg_yyyyyz_xxyyy_0, tg_yyyyyz_xxyyy_1, tg_yyyyyz_xxyyz_0, \
                                         tg_yyyyyz_xxyyz_1, tg_yyyyyz_xxyz_1, tg_yyyyyz_xxyzz_0, tg_yyyyyz_xxyzz_1, \
                                         tg_yyyyyz_xxzz_1, tg_yyyyyz_xxzzz_0, tg_yyyyyz_xxzzz_1, tg_yyyyyz_xyyy_1, \
                                         tg_yyyyyz_xyyyy_0, tg_yyyyyz_xyyyy_1, tg_yyyyyz_xyyyz_0, tg_yyyyyz_xyyyz_1, \
                                         tg_yyyyyz_xyyz_1, tg_yyyyyz_xyzz_1, tg_yyyyyz_xzzz_1, tg_yyyyyz_yyyy_1, \
                                         tg_yyyyyz_yyyz_1, tg_yyzzz_xxxxz_0, tg_yyzzz_xxxxz_1, tg_yyzzz_xxxyy_0, \
                                         tg_yyzzz_xxxyy_1, tg_yyzzz_xxxyz_0, tg_yyzzz_xxxyz_1, tg_yyzzz_xxxzz_0, \
                                         tg_yyzzz_xxxzz_1, tg_yyzzz_xxyyy_0, tg_yyzzz_xxyyy_1, tg_yyzzz_xxyyz_0, \
                                         tg_yyzzz_xxyyz_1, tg_yyzzz_xxyzz_0, tg_yyzzz_xxyzz_1, tg_yyzzz_xxzzz_0, \
                                         tg_yyzzz_xxzzz_1, tg_yyzzz_xyyyy_0, tg_yyzzz_xyyyy_1, tg_yyzzz_xyyyz_0, \
                                         tg_yyzzz_xyyyz_1, tg_yyzzz_xyyzz_0, tg_yyzzz_xyyzz_1, tg_yyzzz_xyzzz_0, \
                                         tg_yyzzz_xyzzz_1, tg_yyzzz_xzzzz_0, tg_yyzzz_xzzzz_1, tg_yyzzz_yyyyy_0, \
                                         tg_yyzzz_yyyyy_1, tg_yyzzz_yyyyz_0, tg_yyzzz_yyyyz_1, tg_yyzzz_yyyzz_0, \
                                         tg_yyzzz_yyyzz_1, tg_yyzzz_yyzzz_0, tg_yyzzz_yyzzz_1, tg_yyzzz_yzzzz_0, \
                                         tg_yyzzz_yzzzz_1, tg_yyzzz_zzzzz_0, tg_yyzzz_zzzzz_1, tg_yzzzz_xxxxx_0, \
                                         tg_yzzzz_xxxxx_1, tg_yzzzz_xxxxy_0, tg_yzzzz_xxxxy_1, tg_yzzzz_xxxxz_0, \
                                         tg_yzzzz_xxxxz_1, tg_yzzzz_xxxyy_0, tg_yzzzz_xxxyy_1, tg_yzzzz_xxxyz_0, \
                                         tg_yzzzz_xxxyz_1, tg_yzzzz_xxxzz_0, tg_yzzzz_xxxzz_1, tg_yzzzz_xxyyy_0, \
                                         tg_yzzzz_xxyyy_1, tg_yzzzz_xxyyz_0, tg_yzzzz_xxyyz_1, tg_yzzzz_xxyzz_0, \
                                         tg_yzzzz_xxyzz_1, tg_yzzzz_xxzzz_0, tg_yzzzz_xxzzz_1, tg_yzzzz_xyyyy_0, \
                                         tg_yzzzz_xyyyy_1, tg_yzzzz_xyyyz_0, tg_yzzzz_xyyyz_1, tg_yzzzz_xyyzz_0, \
                                         tg_yzzzz_xyyzz_1, tg_yzzzz_xyzzz_0, tg_yzzzz_xyzzz_1, tg_yzzzz_xzzzz_0, \
                                         tg_yzzzz_xzzzz_1, tg_yzzzz_yyyyy_0, tg_yzzzz_yyyyy_1, tg_yzzzz_yyyyz_0, \
                                         tg_yzzzz_yyyyz_1, tg_yzzzz_yyyzz_0, tg_yzzzz_yyyzz_1, tg_yzzzz_yyzzz_0, \
                                         tg_yzzzz_yyzzz_1, tg_yzzzz_yzzzz_0, tg_yzzzz_yzzzz_1, tg_yzzzz_zzzzz_0, \
                                         tg_yzzzz_zzzzz_1, tg_zzzzz_xxxxx_0, tg_zzzzz_xxxxx_1, tg_zzzzz_xxxxy_0, \
                                         tg_zzzzz_xxxxy_1, tg_zzzzz_xxxxz_0, tg_zzzzz_xxxxz_1, tg_zzzzz_xxxyy_0, \
                                         tg_zzzzz_xxxyy_1, tg_zzzzz_xxxyz_0, tg_zzzzz_xxxyz_1, tg_zzzzz_xxxzz_0, \
                                         tg_zzzzz_xxxzz_1, tg_zzzzz_xxyyy_0, tg_zzzzz_xxyyy_1, tg_zzzzz_xxyyz_0, \
                                         tg_zzzzz_xxyyz_1, tg_zzzzz_xxyzz_0, tg_zzzzz_xxyzz_1, tg_zzzzz_xxzzz_0, \
                                         tg_zzzzz_xxzzz_1, tg_zzzzz_xyyyy_0, tg_zzzzz_xyyyy_1, tg_zzzzz_xyyyz_0, \
                                         tg_zzzzz_xyyyz_1, tg_zzzzz_xyyzz_0, tg_zzzzz_xyyzz_1, tg_zzzzz_xyzzz_0, \
                                         tg_zzzzz_xyzzz_1, tg_zzzzz_xzzzz_0, tg_zzzzz_xzzzz_1, tg_zzzzz_yyyyy_0, \
                                         tg_zzzzz_yyyyy_1, tg_zzzzz_yyyyz_0, tg_zzzzz_yyyyz_1, tg_zzzzz_yyyzz_0, \
                                         tg_zzzzz_yyyzz_1, tg_zzzzz_yyzzz_0, tg_zzzzz_yyzzz_1, tg_zzzzz_yzzzz_0, \
                                         tg_zzzzz_yzzzz_1, tg_zzzzz_zzzzz_0, tg_zzzzz_zzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyzzz_xxxxz_0[j] = pb_x * tg_xyyzzz_xxxxz_0[j] + fr * tg_xyyzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxz_0[j] - tg_yyzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzz_xxxz_1[j];

                    tg_xxyyzzz_xxxyy_0[j] = pb_x * tg_xyyzzz_xxxyy_0[j] + fr * tg_xyyzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxyy_0[j] - tg_yyzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxyy_1[j];

                    tg_xxyyzzz_xxxyz_0[j] = pb_x * tg_xyyzzz_xxxyz_0[j] + fr * tg_xyyzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxyz_0[j] - tg_yyzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxyz_1[j];

                    tg_xxyyzzz_xxxzz_0[j] = pb_x * tg_xyyzzz_xxxzz_0[j] + fr * tg_xyyzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxzz_0[j] - tg_yyzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxzz_1[j];

                    tg_xxyyzzz_xxyyy_0[j] = pb_x * tg_xyyzzz_xxyyy_0[j] + fr * tg_xyyzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyyy_0[j] - tg_yyzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyyy_1[j];

                    tg_xxyyzzz_xxyyz_0[j] = pb_x * tg_xyyzzz_xxyyz_0[j] + fr * tg_xyyzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyyz_0[j] - tg_yyzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyyz_1[j];

                    tg_xxyyzzz_xxyzz_0[j] = pb_x * tg_xyyzzz_xxyzz_0[j] + fr * tg_xyyzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyzz_0[j] - tg_yyzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyzz_1[j];

                    tg_xxyyzzz_xxzzz_0[j] = pb_x * tg_xyyzzz_xxzzz_0[j] + fr * tg_xyyzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxzzz_0[j] - tg_yyzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xzzz_1[j];

                    tg_xxyyzzz_xyyyy_0[j] = pb_x * tg_xyyzzz_xyyyy_0[j] + fr * tg_xyyzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyyy_0[j] - tg_yyzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyyy_1[j];

                    tg_xxyyzzz_xyyyz_0[j] = pb_x * tg_xyyzzz_xyyyz_0[j] + fr * tg_xyyzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyyz_0[j] - tg_yyzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyyz_1[j];

                    tg_xxyyzzz_xyyzz_0[j] = pb_x * tg_xyyzzz_xyyzz_0[j] + fr * tg_xyyzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyzz_0[j] - tg_yyzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyzz_1[j];

                    tg_xxyyzzz_xyzzz_0[j] = pb_x * tg_xyyzzz_xyzzz_0[j] + fr * tg_xyyzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyzzz_0[j] - tg_yyzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yzzz_1[j];

                    tg_xxyyzzz_xzzzz_0[j] = pb_x * tg_xyyzzz_xzzzz_0[j] + fr * tg_xyyzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xzzzz_0[j] - tg_yyzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_zzzz_1[j];

                    tg_xxyyzzz_yyyyy_0[j] = pb_x * tg_xyyzzz_yyyyy_0[j] + fr * tg_xyyzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyyy_0[j] - tg_yyzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyzzz_yyyyz_0[j] = pb_x * tg_xyyzzz_yyyyz_0[j] + fr * tg_xyyzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyyz_0[j] - tg_yyzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyzzz_yyyzz_0[j] = pb_x * tg_xyyzzz_yyyzz_0[j] + fr * tg_xyyzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyzz_0[j] - tg_yyzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyzzz_yyzzz_0[j] = pb_x * tg_xyyzzz_yyzzz_0[j] + fr * tg_xyyzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyzzz_0[j] - tg_yyzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_yzzzz_0[j] = pb_x * tg_xyyzzz_yzzzz_0[j] + fr * tg_xyyzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yzzzz_0[j] - tg_yyzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_zzzzz_0[j] = pb_x * tg_xyyzzz_zzzzz_0[j] + fr * tg_xyyzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_zzzzz_0[j] - tg_yyzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_xxxxx_0[j] = pb_x * tg_xyzzzz_xxxxx_0[j] + fr * tg_xyzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxx_0[j] - tg_yzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzzzz_xxxx_1[j];

                    tg_xxyzzzz_xxxxy_0[j] = pb_x * tg_xyzzzz_xxxxy_0[j] + fr * tg_xyzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxy_0[j] - tg_yzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzz_xxxy_1[j];

                    tg_xxyzzzz_xxxxz_0[j] = pb_x * tg_xyzzzz_xxxxz_0[j] + fr * tg_xyzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxz_0[j] - tg_yzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzz_xxxz_1[j];

                    tg_xxyzzzz_xxxyy_0[j] = pb_x * tg_xyzzzz_xxxyy_0[j] + fr * tg_xyzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxyy_0[j] - tg_yzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxyy_1[j];

                    tg_xxyzzzz_xxxyz_0[j] = pb_x * tg_xyzzzz_xxxyz_0[j] + fr * tg_xyzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxyz_0[j] - tg_yzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxyz_1[j];

                    tg_xxyzzzz_xxxzz_0[j] = pb_x * tg_xyzzzz_xxxzz_0[j] + fr * tg_xyzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxzz_0[j] - tg_yzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxzz_1[j];

                    tg_xxyzzzz_xxyyy_0[j] = pb_x * tg_xyzzzz_xxyyy_0[j] + fr * tg_xyzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyyy_0[j] - tg_yzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyyy_1[j];

                    tg_xxyzzzz_xxyyz_0[j] = pb_x * tg_xyzzzz_xxyyz_0[j] + fr * tg_xyzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyyz_0[j] - tg_yzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyyz_1[j];

                    tg_xxyzzzz_xxyzz_0[j] = pb_x * tg_xyzzzz_xxyzz_0[j] + fr * tg_xyzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyzz_0[j] - tg_yzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyzz_1[j];

                    tg_xxyzzzz_xxzzz_0[j] = pb_x * tg_xyzzzz_xxzzz_0[j] + fr * tg_xyzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxzzz_0[j] - tg_yzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xzzz_1[j];

                    tg_xxyzzzz_xyyyy_0[j] = pb_x * tg_xyzzzz_xyyyy_0[j] + fr * tg_xyzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyyy_0[j] - tg_yzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyyy_1[j];

                    tg_xxyzzzz_xyyyz_0[j] = pb_x * tg_xyzzzz_xyyyz_0[j] + fr * tg_xyzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyyz_0[j] - tg_yzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyyz_1[j];

                    tg_xxyzzzz_xyyzz_0[j] = pb_x * tg_xyzzzz_xyyzz_0[j] + fr * tg_xyzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyzz_0[j] - tg_yzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyzz_1[j];

                    tg_xxyzzzz_xyzzz_0[j] = pb_x * tg_xyzzzz_xyzzz_0[j] + fr * tg_xyzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyzzz_0[j] - tg_yzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yzzz_1[j];

                    tg_xxyzzzz_xzzzz_0[j] = pb_x * tg_xyzzzz_xzzzz_0[j] + fr * tg_xyzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xzzzz_0[j] - tg_yzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_zzzz_1[j];

                    tg_xxyzzzz_yyyyy_0[j] = pb_x * tg_xyzzzz_yyyyy_0[j] + fr * tg_xyzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyyy_0[j] - tg_yzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyzzzz_yyyyz_0[j] = pb_x * tg_xyzzzz_yyyyz_0[j] + fr * tg_xyzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyyz_0[j] - tg_yzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyzzzz_yyyzz_0[j] = pb_x * tg_xyzzzz_yyyzz_0[j] + fr * tg_xyzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyzz_0[j] - tg_yzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyzzzz_yyzzz_0[j] = pb_x * tg_xyzzzz_yyzzz_0[j] + fr * tg_xyzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyzzz_0[j] - tg_yzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_yzzzz_0[j] = pb_x * tg_xyzzzz_yzzzz_0[j] + fr * tg_xyzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yzzzz_0[j] - tg_yzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_zzzzz_0[j] = pb_x * tg_xyzzzz_zzzzz_0[j] + fr * tg_xyzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_zzzzz_0[j] - tg_yzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_xxxxx_0[j] = pb_x * tg_xzzzzz_xxxxx_0[j] + fr * tg_xzzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxx_0[j] - tg_zzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzzzz_xxxx_1[j];

                    tg_xxzzzzz_xxxxy_0[j] = pb_x * tg_xzzzzz_xxxxy_0[j] + fr * tg_xzzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxy_0[j] - tg_zzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzz_xxxy_1[j];

                    tg_xxzzzzz_xxxxz_0[j] = pb_x * tg_xzzzzz_xxxxz_0[j] + fr * tg_xzzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxz_0[j] - tg_zzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzz_xxxz_1[j];

                    tg_xxzzzzz_xxxyy_0[j] = pb_x * tg_xzzzzz_xxxyy_0[j] + fr * tg_xzzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyy_0[j] - tg_zzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxyy_1[j];

                    tg_xxzzzzz_xxxyz_0[j] = pb_x * tg_xzzzzz_xxxyz_0[j] + fr * tg_xzzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyz_0[j] - tg_zzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxyz_1[j];

                    tg_xxzzzzz_xxxzz_0[j] = pb_x * tg_xzzzzz_xxxzz_0[j] + fr * tg_xzzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxzz_0[j] - tg_zzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxzz_1[j];

                    tg_xxzzzzz_xxyyy_0[j] = pb_x * tg_xzzzzz_xxyyy_0[j] + fr * tg_xzzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyy_0[j] - tg_zzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyyy_1[j];

                    tg_xxzzzzz_xxyyz_0[j] = pb_x * tg_xzzzzz_xxyyz_0[j] + fr * tg_xzzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyz_0[j] - tg_zzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyyz_1[j];

                    tg_xxzzzzz_xxyzz_0[j] = pb_x * tg_xzzzzz_xxyzz_0[j] + fr * tg_xzzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyzz_0[j] - tg_zzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyzz_1[j];

                    tg_xxzzzzz_xxzzz_0[j] = pb_x * tg_xzzzzz_xxzzz_0[j] + fr * tg_xzzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxzzz_0[j] - tg_zzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xzzz_1[j];

                    tg_xxzzzzz_xyyyy_0[j] = pb_x * tg_xzzzzz_xyyyy_0[j] + fr * tg_xzzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyy_0[j] - tg_zzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyyy_1[j];

                    tg_xxzzzzz_xyyyz_0[j] = pb_x * tg_xzzzzz_xyyyz_0[j] + fr * tg_xzzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyz_0[j] - tg_zzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyyz_1[j];

                    tg_xxzzzzz_xyyzz_0[j] = pb_x * tg_xzzzzz_xyyzz_0[j] + fr * tg_xzzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyzz_0[j] - tg_zzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyzz_1[j];

                    tg_xxzzzzz_xyzzz_0[j] = pb_x * tg_xzzzzz_xyzzz_0[j] + fr * tg_xzzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyzzz_0[j] - tg_zzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yzzz_1[j];

                    tg_xxzzzzz_xzzzz_0[j] = pb_x * tg_xzzzzz_xzzzz_0[j] + fr * tg_xzzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzzzz_0[j] - tg_zzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_zzzz_1[j];

                    tg_xxzzzzz_yyyyy_0[j] = pb_x * tg_xzzzzz_yyyyy_0[j] + fr * tg_xzzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyy_0[j] - tg_zzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxzzzzz_yyyyz_0[j] = pb_x * tg_xzzzzz_yyyyz_0[j] + fr * tg_xzzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyz_0[j] - tg_zzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxzzzzz_yyyzz_0[j] = pb_x * tg_xzzzzz_yyyzz_0[j] + fr * tg_xzzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyzz_0[j] - tg_zzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxzzzzz_yyzzz_0[j] = pb_x * tg_xzzzzz_yyzzz_0[j] + fr * tg_xzzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyzzz_0[j] - tg_zzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_yzzzz_0[j] = pb_x * tg_xzzzzz_yzzzz_0[j] + fr * tg_xzzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzzzz_0[j] - tg_zzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_zzzzz_0[j] = pb_x * tg_xzzzzz_zzzzz_0[j] + fr * tg_xzzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzzzz_0[j] - tg_zzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xyyyyyy_xxxxx_0[j] = pb_x * tg_yyyyyy_xxxxx_0[j] + fr * tg_yyyyyy_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyyy_xxxx_1[j];

                    tg_xyyyyyy_xxxxy_0[j] = pb_x * tg_yyyyyy_xxxxy_0[j] + fr * tg_yyyyyy_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyyy_xxxy_1[j];

                    tg_xyyyyyy_xxxxz_0[j] = pb_x * tg_yyyyyy_xxxxz_0[j] + fr * tg_yyyyyy_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyyy_xxxz_1[j];

                    tg_xyyyyyy_xxxyy_0[j] = pb_x * tg_yyyyyy_xxxyy_0[j] + fr * tg_yyyyyy_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxyy_1[j];

                    tg_xyyyyyy_xxxyz_0[j] = pb_x * tg_yyyyyy_xxxyz_0[j] + fr * tg_yyyyyy_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxyz_1[j];

                    tg_xyyyyyy_xxxzz_0[j] = pb_x * tg_yyyyyy_xxxzz_0[j] + fr * tg_yyyyyy_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxzz_1[j];

                    tg_xyyyyyy_xxyyy_0[j] = pb_x * tg_yyyyyy_xxyyy_0[j] + fr * tg_yyyyyy_xxyyy_1[j] + fl1_fxn * tg_yyyyyy_xyyy_1[j];

                    tg_xyyyyyy_xxyyz_0[j] = pb_x * tg_yyyyyy_xxyyz_0[j] + fr * tg_yyyyyy_xxyyz_1[j] + fl1_fxn * tg_yyyyyy_xyyz_1[j];

                    tg_xyyyyyy_xxyzz_0[j] = pb_x * tg_yyyyyy_xxyzz_0[j] + fr * tg_yyyyyy_xxyzz_1[j] + fl1_fxn * tg_yyyyyy_xyzz_1[j];

                    tg_xyyyyyy_xxzzz_0[j] = pb_x * tg_yyyyyy_xxzzz_0[j] + fr * tg_yyyyyy_xxzzz_1[j] + fl1_fxn * tg_yyyyyy_xzzz_1[j];

                    tg_xyyyyyy_xyyyy_0[j] = pb_x * tg_yyyyyy_xyyyy_0[j] + fr * tg_yyyyyy_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyyy_1[j];

                    tg_xyyyyyy_xyyyz_0[j] = pb_x * tg_yyyyyy_xyyyz_0[j] + fr * tg_yyyyyy_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyyz_1[j];

                    tg_xyyyyyy_xyyzz_0[j] = pb_x * tg_yyyyyy_xyyzz_0[j] + fr * tg_yyyyyy_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyzz_1[j];

                    tg_xyyyyyy_xyzzz_0[j] = pb_x * tg_yyyyyy_xyzzz_0[j] + fr * tg_yyyyyy_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yzzz_1[j];

                    tg_xyyyyyy_xzzzz_0[j] = pb_x * tg_yyyyyy_xzzzz_0[j] + fr * tg_yyyyyy_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_zzzz_1[j];

                    tg_xyyyyyy_yyyyy_0[j] = pb_x * tg_yyyyyy_yyyyy_0[j] + fr * tg_yyyyyy_yyyyy_1[j];

                    tg_xyyyyyy_yyyyz_0[j] = pb_x * tg_yyyyyy_yyyyz_0[j] + fr * tg_yyyyyy_yyyyz_1[j];

                    tg_xyyyyyy_yyyzz_0[j] = pb_x * tg_yyyyyy_yyyzz_0[j] + fr * tg_yyyyyy_yyyzz_1[j];

                    tg_xyyyyyy_yyzzz_0[j] = pb_x * tg_yyyyyy_yyzzz_0[j] + fr * tg_yyyyyy_yyzzz_1[j];

                    tg_xyyyyyy_yzzzz_0[j] = pb_x * tg_yyyyyy_yzzzz_0[j] + fr * tg_yyyyyy_yzzzz_1[j];

                    tg_xyyyyyy_zzzzz_0[j] = pb_x * tg_yyyyyy_zzzzz_0[j] + fr * tg_yyyyyy_zzzzz_1[j];

                    tg_xyyyyyz_xxxxx_0[j] = pb_x * tg_yyyyyz_xxxxx_0[j] + fr * tg_yyyyyz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyyz_xxxx_1[j];

                    tg_xyyyyyz_xxxxy_0[j] = pb_x * tg_yyyyyz_xxxxy_0[j] + fr * tg_yyyyyz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyyz_xxxy_1[j];

                    tg_xyyyyyz_xxxxz_0[j] = pb_x * tg_yyyyyz_xxxxz_0[j] + fr * tg_yyyyyz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyyz_xxxz_1[j];

                    tg_xyyyyyz_xxxyy_0[j] = pb_x * tg_yyyyyz_xxxyy_0[j] + fr * tg_yyyyyz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxyy_1[j];

                    tg_xyyyyyz_xxxyz_0[j] = pb_x * tg_yyyyyz_xxxyz_0[j] + fr * tg_yyyyyz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxyz_1[j];

                    tg_xyyyyyz_xxxzz_0[j] = pb_x * tg_yyyyyz_xxxzz_0[j] + fr * tg_yyyyyz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxzz_1[j];

                    tg_xyyyyyz_xxyyy_0[j] = pb_x * tg_yyyyyz_xxyyy_0[j] + fr * tg_yyyyyz_xxyyy_1[j] + fl1_fxn * tg_yyyyyz_xyyy_1[j];

                    tg_xyyyyyz_xxyyz_0[j] = pb_x * tg_yyyyyz_xxyyz_0[j] + fr * tg_yyyyyz_xxyyz_1[j] + fl1_fxn * tg_yyyyyz_xyyz_1[j];

                    tg_xyyyyyz_xxyzz_0[j] = pb_x * tg_yyyyyz_xxyzz_0[j] + fr * tg_yyyyyz_xxyzz_1[j] + fl1_fxn * tg_yyyyyz_xyzz_1[j];

                    tg_xyyyyyz_xxzzz_0[j] = pb_x * tg_yyyyyz_xxzzz_0[j] + fr * tg_yyyyyz_xxzzz_1[j] + fl1_fxn * tg_yyyyyz_xzzz_1[j];

                    tg_xyyyyyz_xyyyy_0[j] = pb_x * tg_yyyyyz_xyyyy_0[j] + fr * tg_yyyyyz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyyy_1[j];

                    tg_xyyyyyz_xyyyz_0[j] = pb_x * tg_yyyyyz_xyyyz_0[j] + fr * tg_yyyyyz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSH_474_568(      CMemBlock2D<double>* primBuffer,
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

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyyyyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 474); 

                auto tg_yyyyyz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 475); 

                auto tg_yyyyyz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 476); 

                auto tg_yyyyyz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 477); 

                auto tg_yyyyyz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 478); 

                auto tg_yyyyyz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 479); 

                auto tg_yyyyyz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 480); 

                auto tg_yyyyyz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 481); 

                auto tg_yyyyyz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 482); 

                auto tg_yyyyzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 483); 

                auto tg_yyyyzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 484); 

                auto tg_yyyyzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 485); 

                auto tg_yyyyzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 486); 

                auto tg_yyyyzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 487); 

                auto tg_yyyyzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 488); 

                auto tg_yyyyzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 489); 

                auto tg_yyyyzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 490); 

                auto tg_yyyyzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 491); 

                auto tg_yyyyzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 492); 

                auto tg_yyyyzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 493); 

                auto tg_yyyyzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 494); 

                auto tg_yyyyzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 495); 

                auto tg_yyyyzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 496); 

                auto tg_yyyyzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 497); 

                auto tg_yyyyzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 498); 

                auto tg_yyyyzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 499); 

                auto tg_yyyyzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 500); 

                auto tg_yyyyzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 501); 

                auto tg_yyyyzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 502); 

                auto tg_yyyyzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 503); 

                auto tg_yyyzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 504); 

                auto tg_yyyzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 505); 

                auto tg_yyyzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 506); 

                auto tg_yyyzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 507); 

                auto tg_yyyzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 508); 

                auto tg_yyyzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 509); 

                auto tg_yyyzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 510); 

                auto tg_yyyzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 511); 

                auto tg_yyyzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 512); 

                auto tg_yyyzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 513); 

                auto tg_yyyzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 514); 

                auto tg_yyyzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 515); 

                auto tg_yyyzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 516); 

                auto tg_yyyzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 517); 

                auto tg_yyyzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 518); 

                auto tg_yyyzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 519); 

                auto tg_yyyzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 520); 

                auto tg_yyyzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 521); 

                auto tg_yyyzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 522); 

                auto tg_yyyzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 523); 

                auto tg_yyyzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 524); 

                auto tg_yyzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 525); 

                auto tg_yyzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 526); 

                auto tg_yyzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 527); 

                auto tg_yyzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 528); 

                auto tg_yyzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 529); 

                auto tg_yyzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 530); 

                auto tg_yyzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 531); 

                auto tg_yyzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 532); 

                auto tg_yyzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 533); 

                auto tg_yyzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 534); 

                auto tg_yyzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 535); 

                auto tg_yyzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 536); 

                auto tg_yyzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 537); 

                auto tg_yyzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 538); 

                auto tg_yyzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 539); 

                auto tg_yyzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 540); 

                auto tg_yyzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 541); 

                auto tg_yyzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 542); 

                auto tg_yyzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 543); 

                auto tg_yyzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 544); 

                auto tg_yyzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 545); 

                auto tg_yzzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 546); 

                auto tg_yzzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 547); 

                auto tg_yzzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 548); 

                auto tg_yzzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 549); 

                auto tg_yzzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 550); 

                auto tg_yzzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 551); 

                auto tg_yzzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 552); 

                auto tg_yzzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 553); 

                auto tg_yzzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 554); 

                auto tg_yzzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 555); 

                auto tg_yzzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 556); 

                auto tg_yzzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 557); 

                auto tg_yzzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 558); 

                auto tg_yzzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 559); 

                auto tg_yzzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 560); 

                auto tg_yzzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 561); 

                auto tg_yzzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 562); 

                auto tg_yzzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 563); 

                auto tg_yzzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 564); 

                auto tg_yzzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 565); 

                auto tg_yzzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 566); 

                auto tg_zzzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 567); 

                auto tg_yyyyyz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 474); 

                auto tg_yyyyyz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 475); 

                auto tg_yyyyyz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 476); 

                auto tg_yyyyyz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 477); 

                auto tg_yyyyyz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 478); 

                auto tg_yyyyyz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 479); 

                auto tg_yyyyyz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 480); 

                auto tg_yyyyyz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 481); 

                auto tg_yyyyyz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 482); 

                auto tg_yyyyzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 483); 

                auto tg_yyyyzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 484); 

                auto tg_yyyyzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 485); 

                auto tg_yyyyzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 486); 

                auto tg_yyyyzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 487); 

                auto tg_yyyyzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 488); 

                auto tg_yyyyzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 489); 

                auto tg_yyyyzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 490); 

                auto tg_yyyyzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 491); 

                auto tg_yyyyzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 492); 

                auto tg_yyyyzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 493); 

                auto tg_yyyyzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 494); 

                auto tg_yyyyzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 495); 

                auto tg_yyyyzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 496); 

                auto tg_yyyyzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 497); 

                auto tg_yyyyzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 498); 

                auto tg_yyyyzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 499); 

                auto tg_yyyyzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 500); 

                auto tg_yyyyzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 501); 

                auto tg_yyyyzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 502); 

                auto tg_yyyyzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 503); 

                auto tg_yyyzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 504); 

                auto tg_yyyzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 505); 

                auto tg_yyyzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 506); 

                auto tg_yyyzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 507); 

                auto tg_yyyzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 508); 

                auto tg_yyyzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 509); 

                auto tg_yyyzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 510); 

                auto tg_yyyzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 511); 

                auto tg_yyyzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 512); 

                auto tg_yyyzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 513); 

                auto tg_yyyzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 514); 

                auto tg_yyyzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 515); 

                auto tg_yyyzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 516); 

                auto tg_yyyzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 517); 

                auto tg_yyyzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 518); 

                auto tg_yyyzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 519); 

                auto tg_yyyzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 520); 

                auto tg_yyyzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 521); 

                auto tg_yyyzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 522); 

                auto tg_yyyzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 523); 

                auto tg_yyyzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 524); 

                auto tg_yyzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 525); 

                auto tg_yyzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 526); 

                auto tg_yyzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 527); 

                auto tg_yyzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 528); 

                auto tg_yyzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 529); 

                auto tg_yyzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 530); 

                auto tg_yyzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 531); 

                auto tg_yyzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 532); 

                auto tg_yyzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 533); 

                auto tg_yyzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 534); 

                auto tg_yyzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 535); 

                auto tg_yyzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 536); 

                auto tg_yyzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 537); 

                auto tg_yyzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 538); 

                auto tg_yyzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 539); 

                auto tg_yyzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 540); 

                auto tg_yyzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 541); 

                auto tg_yyzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 542); 

                auto tg_yyzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 543); 

                auto tg_yyzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 544); 

                auto tg_yyzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 545); 

                auto tg_yzzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 546); 

                auto tg_yzzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 547); 

                auto tg_yzzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 548); 

                auto tg_yzzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 549); 

                auto tg_yzzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 550); 

                auto tg_yzzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 551); 

                auto tg_yzzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 552); 

                auto tg_yzzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 553); 

                auto tg_yzzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 554); 

                auto tg_yzzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 555); 

                auto tg_yzzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 556); 

                auto tg_yzzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 557); 

                auto tg_yzzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 558); 

                auto tg_yzzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 559); 

                auto tg_yzzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 560); 

                auto tg_yzzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 561); 

                auto tg_yzzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 562); 

                auto tg_yzzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 563); 

                auto tg_yzzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 564); 

                auto tg_yzzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 565); 

                auto tg_yzzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 566); 

                auto tg_zzzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 567); 

                auto tg_yyyyyz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 342); 

                auto tg_yyyyyz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 343); 

                auto tg_yyyyyz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 344); 

                auto tg_yyyyzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 345); 

                auto tg_yyyyzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 346); 

                auto tg_yyyyzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 347); 

                auto tg_yyyyzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 348); 

                auto tg_yyyyzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 349); 

                auto tg_yyyyzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 350); 

                auto tg_yyyyzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 351); 

                auto tg_yyyyzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 352); 

                auto tg_yyyyzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 353); 

                auto tg_yyyyzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 354); 

                auto tg_yyyyzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 355); 

                auto tg_yyyyzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 356); 

                auto tg_yyyyzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 357); 

                auto tg_yyyyzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 358); 

                auto tg_yyyyzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 359); 

                auto tg_yyyzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 360); 

                auto tg_yyyzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 361); 

                auto tg_yyyzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 362); 

                auto tg_yyyzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 363); 

                auto tg_yyyzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 364); 

                auto tg_yyyzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 365); 

                auto tg_yyyzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 366); 

                auto tg_yyyzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 367); 

                auto tg_yyyzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 368); 

                auto tg_yyyzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 369); 

                auto tg_yyyzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 370); 

                auto tg_yyyzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 371); 

                auto tg_yyyzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 372); 

                auto tg_yyyzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 373); 

                auto tg_yyyzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 374); 

                auto tg_yyzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 375); 

                auto tg_yyzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 376); 

                auto tg_yyzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 377); 

                auto tg_yyzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 378); 

                auto tg_yyzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 379); 

                auto tg_yyzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 380); 

                auto tg_yyzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 381); 

                auto tg_yyzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 382); 

                auto tg_yyzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 383); 

                auto tg_yyzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 384); 

                auto tg_yyzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 385); 

                auto tg_yyzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 386); 

                auto tg_yyzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 387); 

                auto tg_yyzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 388); 

                auto tg_yyzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 389); 

                auto tg_yzzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 390); 

                auto tg_yzzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 391); 

                auto tg_yzzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 392); 

                auto tg_yzzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 393); 

                auto tg_yzzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 394); 

                auto tg_yzzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 395); 

                auto tg_yzzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 396); 

                auto tg_yzzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 397); 

                auto tg_yzzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 398); 

                auto tg_yzzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 399); 

                auto tg_yzzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 400); 

                auto tg_yzzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 401); 

                auto tg_yzzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 402); 

                auto tg_yzzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 403); 

                auto tg_yzzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 404); 

                auto tg_zzzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 405); 

                // set up pointers to integrals

                auto tg_xyyyyyz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 474); 

                auto tg_xyyyyyz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 475); 

                auto tg_xyyyyyz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 476); 

                auto tg_xyyyyyz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 477); 

                auto tg_xyyyyyz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 478); 

                auto tg_xyyyyyz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 479); 

                auto tg_xyyyyyz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 480); 

                auto tg_xyyyyyz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 481); 

                auto tg_xyyyyyz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 482); 

                auto tg_xyyyyzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 483); 

                auto tg_xyyyyzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 484); 

                auto tg_xyyyyzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 485); 

                auto tg_xyyyyzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 486); 

                auto tg_xyyyyzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 487); 

                auto tg_xyyyyzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 488); 

                auto tg_xyyyyzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 489); 

                auto tg_xyyyyzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 490); 

                auto tg_xyyyyzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 491); 

                auto tg_xyyyyzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 492); 

                auto tg_xyyyyzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 493); 

                auto tg_xyyyyzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 494); 

                auto tg_xyyyyzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 495); 

                auto tg_xyyyyzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 496); 

                auto tg_xyyyyzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 497); 

                auto tg_xyyyyzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 498); 

                auto tg_xyyyyzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 499); 

                auto tg_xyyyyzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 500); 

                auto tg_xyyyyzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 501); 

                auto tg_xyyyyzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 502); 

                auto tg_xyyyyzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 503); 

                auto tg_xyyyzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 504); 

                auto tg_xyyyzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 505); 

                auto tg_xyyyzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 506); 

                auto tg_xyyyzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 507); 

                auto tg_xyyyzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 508); 

                auto tg_xyyyzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 509); 

                auto tg_xyyyzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 510); 

                auto tg_xyyyzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 511); 

                auto tg_xyyyzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 512); 

                auto tg_xyyyzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 513); 

                auto tg_xyyyzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 514); 

                auto tg_xyyyzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 515); 

                auto tg_xyyyzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 516); 

                auto tg_xyyyzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 517); 

                auto tg_xyyyzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 518); 

                auto tg_xyyyzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 519); 

                auto tg_xyyyzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 520); 

                auto tg_xyyyzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 521); 

                auto tg_xyyyzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 522); 

                auto tg_xyyyzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 523); 

                auto tg_xyyyzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 524); 

                auto tg_xyyzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 525); 

                auto tg_xyyzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 526); 

                auto tg_xyyzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 527); 

                auto tg_xyyzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 528); 

                auto tg_xyyzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 529); 

                auto tg_xyyzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 530); 

                auto tg_xyyzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 531); 

                auto tg_xyyzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 532); 

                auto tg_xyyzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 533); 

                auto tg_xyyzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 534); 

                auto tg_xyyzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 535); 

                auto tg_xyyzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 536); 

                auto tg_xyyzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 537); 

                auto tg_xyyzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 538); 

                auto tg_xyyzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 539); 

                auto tg_xyyzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 540); 

                auto tg_xyyzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 541); 

                auto tg_xyyzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 542); 

                auto tg_xyyzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 543); 

                auto tg_xyyzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 544); 

                auto tg_xyyzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 545); 

                auto tg_xyzzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 546); 

                auto tg_xyzzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 547); 

                auto tg_xyzzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 548); 

                auto tg_xyzzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 549); 

                auto tg_xyzzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 550); 

                auto tg_xyzzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 551); 

                auto tg_xyzzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 552); 

                auto tg_xyzzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 553); 

                auto tg_xyzzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 554); 

                auto tg_xyzzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 555); 

                auto tg_xyzzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 556); 

                auto tg_xyzzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 557); 

                auto tg_xyzzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 558); 

                auto tg_xyzzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 559); 

                auto tg_xyzzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 560); 

                auto tg_xyzzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 561); 

                auto tg_xyzzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 562); 

                auto tg_xyzzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 563); 

                auto tg_xyzzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 564); 

                auto tg_xyzzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 565); 

                auto tg_xyzzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 566); 

                auto tg_xzzzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 567); 

                // Batch of Integrals (474,568)

                #pragma omp simd aligned(fxn, tg_xyyyyyz_xyyzz_0, tg_xyyyyyz_xyzzz_0, tg_xyyyyyz_xzzzz_0, \
                                         tg_xyyyyyz_yyyyy_0, tg_xyyyyyz_yyyyz_0, tg_xyyyyyz_yyyzz_0, tg_xyyyyyz_yyzzz_0, \
                                         tg_xyyyyyz_yzzzz_0, tg_xyyyyyz_zzzzz_0, tg_xyyyyzz_xxxxx_0, tg_xyyyyzz_xxxxy_0, \
                                         tg_xyyyyzz_xxxxz_0, tg_xyyyyzz_xxxyy_0, tg_xyyyyzz_xxxyz_0, tg_xyyyyzz_xxxzz_0, \
                                         tg_xyyyyzz_xxyyy_0, tg_xyyyyzz_xxyyz_0, tg_xyyyyzz_xxyzz_0, tg_xyyyyzz_xxzzz_0, \
                                         tg_xyyyyzz_xyyyy_0, tg_xyyyyzz_xyyyz_0, tg_xyyyyzz_xyyzz_0, tg_xyyyyzz_xyzzz_0, \
                                         tg_xyyyyzz_xzzzz_0, tg_xyyyyzz_yyyyy_0, tg_xyyyyzz_yyyyz_0, tg_xyyyyzz_yyyzz_0, \
                                         tg_xyyyyzz_yyzzz_0, tg_xyyyyzz_yzzzz_0, tg_xyyyyzz_zzzzz_0, tg_xyyyzzz_xxxxx_0, \
                                         tg_xyyyzzz_xxxxy_0, tg_xyyyzzz_xxxxz_0, tg_xyyyzzz_xxxyy_0, tg_xyyyzzz_xxxyz_0, \
                                         tg_xyyyzzz_xxxzz_0, tg_xyyyzzz_xxyyy_0, tg_xyyyzzz_xxyyz_0, tg_xyyyzzz_xxyzz_0, \
                                         tg_xyyyzzz_xxzzz_0, tg_xyyyzzz_xyyyy_0, tg_xyyyzzz_xyyyz_0, tg_xyyyzzz_xyyzz_0, \
                                         tg_xyyyzzz_xyzzz_0, tg_xyyyzzz_xzzzz_0, tg_xyyyzzz_yyyyy_0, tg_xyyyzzz_yyyyz_0, \
                                         tg_xyyyzzz_yyyzz_0, tg_xyyyzzz_yyzzz_0, tg_xyyyzzz_yzzzz_0, tg_xyyyzzz_zzzzz_0, \
                                         tg_xyyzzzz_xxxxx_0, tg_xyyzzzz_xxxxy_0, tg_xyyzzzz_xxxxz_0, tg_xyyzzzz_xxxyy_0, \
                                         tg_xyyzzzz_xxxyz_0, tg_xyyzzzz_xxxzz_0, tg_xyyzzzz_xxyyy_0, tg_xyyzzzz_xxyyz_0, \
                                         tg_xyyzzzz_xxyzz_0, tg_xyyzzzz_xxzzz_0, tg_xyyzzzz_xyyyy_0, tg_xyyzzzz_xyyyz_0, \
                                         tg_xyyzzzz_xyyzz_0, tg_xyyzzzz_xyzzz_0, tg_xyyzzzz_xzzzz_0, tg_xyyzzzz_yyyyy_0, \
                                         tg_xyyzzzz_yyyyz_0, tg_xyyzzzz_yyyzz_0, tg_xyyzzzz_yyzzz_0, tg_xyyzzzz_yzzzz_0, \
                                         tg_xyyzzzz_zzzzz_0, tg_xyzzzzz_xxxxx_0, tg_xyzzzzz_xxxxy_0, tg_xyzzzzz_xxxxz_0, \
                                         tg_xyzzzzz_xxxyy_0, tg_xyzzzzz_xxxyz_0, tg_xyzzzzz_xxxzz_0, tg_xyzzzzz_xxyyy_0, \
                                         tg_xyzzzzz_xxyyz_0, tg_xyzzzzz_xxyzz_0, tg_xyzzzzz_xxzzz_0, tg_xyzzzzz_xyyyy_0, \
                                         tg_xyzzzzz_xyyyz_0, tg_xyzzzzz_xyyzz_0, tg_xyzzzzz_xyzzz_0, tg_xyzzzzz_xzzzz_0, \
                                         tg_xyzzzzz_yyyyy_0, tg_xyzzzzz_yyyyz_0, tg_xyzzzzz_yyyzz_0, tg_xyzzzzz_yyzzz_0, \
                                         tg_xyzzzzz_yzzzz_0, tg_xyzzzzz_zzzzz_0, tg_xzzzzzz_xxxxx_0, tg_yyyyyz_xyyzz_0, \
                                         tg_yyyyyz_xyyzz_1, tg_yyyyyz_xyzzz_0, tg_yyyyyz_xyzzz_1, tg_yyyyyz_xzzzz_0, \
                                         tg_yyyyyz_xzzzz_1, tg_yyyyyz_yyyyy_0, tg_yyyyyz_yyyyy_1, tg_yyyyyz_yyyyz_0, \
                                         tg_yyyyyz_yyyyz_1, tg_yyyyyz_yyyzz_0, tg_yyyyyz_yyyzz_1, tg_yyyyyz_yyzz_1, \
                                         tg_yyyyyz_yyzzz_0, tg_yyyyyz_yyzzz_1, tg_yyyyyz_yzzz_1, tg_yyyyyz_yzzzz_0, \
                                         tg_yyyyyz_yzzzz_1, tg_yyyyyz_zzzz_1, tg_yyyyyz_zzzzz_0, tg_yyyyyz_zzzzz_1, \
                                         tg_yyyyzz_xxxx_1, tg_yyyyzz_xxxxx_0, tg_yyyyzz_xxxxx_1, tg_yyyyzz_xxxxy_0, \
                                         tg_yyyyzz_xxxxy_1, tg_yyyyzz_xxxxz_0, tg_yyyyzz_xxxxz_1, tg_yyyyzz_xxxy_1, \
                                         tg_yyyyzz_xxxyy_0, tg_yyyyzz_xxxyy_1, tg_yyyyzz_xxxyz_0, tg_yyyyzz_xxxyz_1, \
                                         tg_yyyyzz_xxxz_1, tg_yyyyzz_xxxzz_0, tg_yyyyzz_xxxzz_1, tg_yyyyzz_xxyy_1, \
                                         tg_yyyyzz_xxyyy_0, tg_yyyyzz_xxyyy_1, tg_yyyyzz_xxyyz_0, tg_yyyyzz_xxyyz_1, \
                                         tg_yyyyzz_xxyz_1, tg_yyyyzz_xxyzz_0, tg_yyyyzz_xxyzz_1, tg_yyyyzz_xxzz_1, \
                                         tg_yyyyzz_xxzzz_0, tg_yyyyzz_xxzzz_1, tg_yyyyzz_xyyy_1, tg_yyyyzz_xyyyy_0, \
                                         tg_yyyyzz_xyyyy_1, tg_yyyyzz_xyyyz_0, tg_yyyyzz_xyyyz_1, tg_yyyyzz_xyyz_1, \
                                         tg_yyyyzz_xyyzz_0, tg_yyyyzz_xyyzz_1, tg_yyyyzz_xyzz_1, tg_yyyyzz_xyzzz_0, \
                                         tg_yyyyzz_xyzzz_1, tg_yyyyzz_xzzz_1, tg_yyyyzz_xzzzz_0, tg_yyyyzz_xzzzz_1, \
                                         tg_yyyyzz_yyyy_1, tg_yyyyzz_yyyyy_0, tg_yyyyzz_yyyyy_1, tg_yyyyzz_yyyyz_0, \
                                         tg_yyyyzz_yyyyz_1, tg_yyyyzz_yyyz_1, tg_yyyyzz_yyyzz_0, tg_yyyyzz_yyyzz_1, \
                                         tg_yyyyzz_yyzz_1, tg_yyyyzz_yyzzz_0, tg_yyyyzz_yyzzz_1, tg_yyyyzz_yzzz_1, \
                                         tg_yyyyzz_yzzzz_0, tg_yyyyzz_yzzzz_1, tg_yyyyzz_zzzz_1, tg_yyyyzz_zzzzz_0, \
                                         tg_yyyyzz_zzzzz_1, tg_yyyzzz_xxxx_1, tg_yyyzzz_xxxxx_0, tg_yyyzzz_xxxxx_1, \
                                         tg_yyyzzz_xxxxy_0, tg_yyyzzz_xxxxy_1, tg_yyyzzz_xxxxz_0, tg_yyyzzz_xxxxz_1, \
                                         tg_yyyzzz_xxxy_1, tg_yyyzzz_xxxyy_0, tg_yyyzzz_xxxyy_1, tg_yyyzzz_xxxyz_0, \
                                         tg_yyyzzz_xxxyz_1, tg_yyyzzz_xxxz_1, tg_yyyzzz_xxxzz_0, tg_yyyzzz_xxxzz_1, \
                                         tg_yyyzzz_xxyy_1, tg_yyyzzz_xxyyy_0, tg_yyyzzz_xxyyy_1, tg_yyyzzz_xxyyz_0, \
                                         tg_yyyzzz_xxyyz_1, tg_yyyzzz_xxyz_1, tg_yyyzzz_xxyzz_0, tg_yyyzzz_xxyzz_1, \
                                         tg_yyyzzz_xxzz_1, tg_yyyzzz_xxzzz_0, tg_yyyzzz_xxzzz_1, tg_yyyzzz_xyyy_1, \
                                         tg_yyyzzz_xyyyy_0, tg_yyyzzz_xyyyy_1, tg_yyyzzz_xyyyz_0, tg_yyyzzz_xyyyz_1, \
                                         tg_yyyzzz_xyyz_1, tg_yyyzzz_xyyzz_0, tg_yyyzzz_xyyzz_1, tg_yyyzzz_xyzz_1, \
                                         tg_yyyzzz_xyzzz_0, tg_yyyzzz_xyzzz_1, tg_yyyzzz_xzzz_1, tg_yyyzzz_xzzzz_0, \
                                         tg_yyyzzz_xzzzz_1, tg_yyyzzz_yyyy_1, tg_yyyzzz_yyyyy_0, tg_yyyzzz_yyyyy_1, \
                                         tg_yyyzzz_yyyyz_0, tg_yyyzzz_yyyyz_1, tg_yyyzzz_yyyz_1, tg_yyyzzz_yyyzz_0, \
                                         tg_yyyzzz_yyyzz_1, tg_yyyzzz_yyzz_1, tg_yyyzzz_yyzzz_0, tg_yyyzzz_yyzzz_1, \
                                         tg_yyyzzz_yzzz_1, tg_yyyzzz_yzzzz_0, tg_yyyzzz_yzzzz_1, tg_yyyzzz_zzzz_1, \
                                         tg_yyyzzz_zzzzz_0, tg_yyyzzz_zzzzz_1, tg_yyzzzz_xxxx_1, tg_yyzzzz_xxxxx_0, \
                                         tg_yyzzzz_xxxxx_1, tg_yyzzzz_xxxxy_0, tg_yyzzzz_xxxxy_1, tg_yyzzzz_xxxxz_0, \
                                         tg_yyzzzz_xxxxz_1, tg_yyzzzz_xxxy_1, tg_yyzzzz_xxxyy_0, tg_yyzzzz_xxxyy_1, \
                                         tg_yyzzzz_xxxyz_0, tg_yyzzzz_xxxyz_1, tg_yyzzzz_xxxz_1, tg_yyzzzz_xxxzz_0, \
                                         tg_yyzzzz_xxxzz_1, tg_yyzzzz_xxyy_1, tg_yyzzzz_xxyyy_0, tg_yyzzzz_xxyyy_1, \
                                         tg_yyzzzz_xxyyz_0, tg_yyzzzz_xxyyz_1, tg_yyzzzz_xxyz_1, tg_yyzzzz_xxyzz_0, \
                                         tg_yyzzzz_xxyzz_1, tg_yyzzzz_xxzz_1, tg_yyzzzz_xxzzz_0, tg_yyzzzz_xxzzz_1, \
                                         tg_yyzzzz_xyyy_1, tg_yyzzzz_xyyyy_0, tg_yyzzzz_xyyyy_1, tg_yyzzzz_xyyyz_0, \
                                         tg_yyzzzz_xyyyz_1, tg_yyzzzz_xyyz_1, tg_yyzzzz_xyyzz_0, tg_yyzzzz_xyyzz_1, \
                                         tg_yyzzzz_xyzz_1, tg_yyzzzz_xyzzz_0, tg_yyzzzz_xyzzz_1, tg_yyzzzz_xzzz_1, \
                                         tg_yyzzzz_xzzzz_0, tg_yyzzzz_xzzzz_1, tg_yyzzzz_yyyy_1, tg_yyzzzz_yyyyy_0, \
                                         tg_yyzzzz_yyyyy_1, tg_yyzzzz_yyyyz_0, tg_yyzzzz_yyyyz_1, tg_yyzzzz_yyyz_1, \
                                         tg_yyzzzz_yyyzz_0, tg_yyzzzz_yyyzz_1, tg_yyzzzz_yyzz_1, tg_yyzzzz_yyzzz_0, \
                                         tg_yyzzzz_yyzzz_1, tg_yyzzzz_yzzz_1, tg_yyzzzz_yzzzz_0, tg_yyzzzz_yzzzz_1, \
                                         tg_yyzzzz_zzzz_1, tg_yyzzzz_zzzzz_0, tg_yyzzzz_zzzzz_1, tg_yzzzzz_xxxx_1, \
                                         tg_yzzzzz_xxxxx_0, tg_yzzzzz_xxxxx_1, tg_yzzzzz_xxxxy_0, tg_yzzzzz_xxxxy_1, \
                                         tg_yzzzzz_xxxxz_0, tg_yzzzzz_xxxxz_1, tg_yzzzzz_xxxy_1, tg_yzzzzz_xxxyy_0, \
                                         tg_yzzzzz_xxxyy_1, tg_yzzzzz_xxxyz_0, tg_yzzzzz_xxxyz_1, tg_yzzzzz_xxxz_1, \
                                         tg_yzzzzz_xxxzz_0, tg_yzzzzz_xxxzz_1, tg_yzzzzz_xxyy_1, tg_yzzzzz_xxyyy_0, \
                                         tg_yzzzzz_xxyyy_1, tg_yzzzzz_xxyyz_0, tg_yzzzzz_xxyyz_1, tg_yzzzzz_xxyz_1, \
                                         tg_yzzzzz_xxyzz_0, tg_yzzzzz_xxyzz_1, tg_yzzzzz_xxzz_1, tg_yzzzzz_xxzzz_0, \
                                         tg_yzzzzz_xxzzz_1, tg_yzzzzz_xyyy_1, tg_yzzzzz_xyyyy_0, tg_yzzzzz_xyyyy_1, \
                                         tg_yzzzzz_xyyyz_0, tg_yzzzzz_xyyyz_1, tg_yzzzzz_xyyz_1, tg_yzzzzz_xyyzz_0, \
                                         tg_yzzzzz_xyyzz_1, tg_yzzzzz_xyzz_1, tg_yzzzzz_xyzzz_0, tg_yzzzzz_xyzzz_1, \
                                         tg_yzzzzz_xzzz_1, tg_yzzzzz_xzzzz_0, tg_yzzzzz_xzzzz_1, tg_yzzzzz_yyyy_1, \
                                         tg_yzzzzz_yyyyy_0, tg_yzzzzz_yyyyy_1, tg_yzzzzz_yyyyz_0, tg_yzzzzz_yyyyz_1, \
                                         tg_yzzzzz_yyyz_1, tg_yzzzzz_yyyzz_0, tg_yzzzzz_yyyzz_1, tg_yzzzzz_yyzz_1, \
                                         tg_yzzzzz_yyzzz_0, tg_yzzzzz_yyzzz_1, tg_yzzzzz_yzzz_1, tg_yzzzzz_yzzzz_0, \
                                         tg_yzzzzz_yzzzz_1, tg_yzzzzz_zzzz_1, tg_yzzzzz_zzzzz_0, tg_yzzzzz_zzzzz_1, \
                                         tg_zzzzzz_xxxx_1, tg_zzzzzz_xxxxx_0, tg_zzzzzz_xxxxx_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyyyyz_xyyzz_0[j] = pb_x * tg_yyyyyz_xyyzz_0[j] + fr * tg_yyyyyz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyzz_1[j];

                    tg_xyyyyyz_xyzzz_0[j] = pb_x * tg_yyyyyz_xyzzz_0[j] + fr * tg_yyyyyz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yzzz_1[j];

                    tg_xyyyyyz_xzzzz_0[j] = pb_x * tg_yyyyyz_xzzzz_0[j] + fr * tg_yyyyyz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_zzzz_1[j];

                    tg_xyyyyyz_yyyyy_0[j] = pb_x * tg_yyyyyz_yyyyy_0[j] + fr * tg_yyyyyz_yyyyy_1[j];

                    tg_xyyyyyz_yyyyz_0[j] = pb_x * tg_yyyyyz_yyyyz_0[j] + fr * tg_yyyyyz_yyyyz_1[j];

                    tg_xyyyyyz_yyyzz_0[j] = pb_x * tg_yyyyyz_yyyzz_0[j] + fr * tg_yyyyyz_yyyzz_1[j];

                    tg_xyyyyyz_yyzzz_0[j] = pb_x * tg_yyyyyz_yyzzz_0[j] + fr * tg_yyyyyz_yyzzz_1[j];

                    tg_xyyyyyz_yzzzz_0[j] = pb_x * tg_yyyyyz_yzzzz_0[j] + fr * tg_yyyyyz_yzzzz_1[j];

                    tg_xyyyyyz_zzzzz_0[j] = pb_x * tg_yyyyyz_zzzzz_0[j] + fr * tg_yyyyyz_zzzzz_1[j];

                    tg_xyyyyzz_xxxxx_0[j] = pb_x * tg_yyyyzz_xxxxx_0[j] + fr * tg_yyyyzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyzz_xxxx_1[j];

                    tg_xyyyyzz_xxxxy_0[j] = pb_x * tg_yyyyzz_xxxxy_0[j] + fr * tg_yyyyzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyzz_xxxy_1[j];

                    tg_xyyyyzz_xxxxz_0[j] = pb_x * tg_yyyyzz_xxxxz_0[j] + fr * tg_yyyyzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyzz_xxxz_1[j];

                    tg_xyyyyzz_xxxyy_0[j] = pb_x * tg_yyyyzz_xxxyy_0[j] + fr * tg_yyyyzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxyy_1[j];

                    tg_xyyyyzz_xxxyz_0[j] = pb_x * tg_yyyyzz_xxxyz_0[j] + fr * tg_yyyyzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxyz_1[j];

                    tg_xyyyyzz_xxxzz_0[j] = pb_x * tg_yyyyzz_xxxzz_0[j] + fr * tg_yyyyzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxzz_1[j];

                    tg_xyyyyzz_xxyyy_0[j] = pb_x * tg_yyyyzz_xxyyy_0[j] + fr * tg_yyyyzz_xxyyy_1[j] + fl1_fxn * tg_yyyyzz_xyyy_1[j];

                    tg_xyyyyzz_xxyyz_0[j] = pb_x * tg_yyyyzz_xxyyz_0[j] + fr * tg_yyyyzz_xxyyz_1[j] + fl1_fxn * tg_yyyyzz_xyyz_1[j];

                    tg_xyyyyzz_xxyzz_0[j] = pb_x * tg_yyyyzz_xxyzz_0[j] + fr * tg_yyyyzz_xxyzz_1[j] + fl1_fxn * tg_yyyyzz_xyzz_1[j];

                    tg_xyyyyzz_xxzzz_0[j] = pb_x * tg_yyyyzz_xxzzz_0[j] + fr * tg_yyyyzz_xxzzz_1[j] + fl1_fxn * tg_yyyyzz_xzzz_1[j];

                    tg_xyyyyzz_xyyyy_0[j] = pb_x * tg_yyyyzz_xyyyy_0[j] + fr * tg_yyyyzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyyy_1[j];

                    tg_xyyyyzz_xyyyz_0[j] = pb_x * tg_yyyyzz_xyyyz_0[j] + fr * tg_yyyyzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyyz_1[j];

                    tg_xyyyyzz_xyyzz_0[j] = pb_x * tg_yyyyzz_xyyzz_0[j] + fr * tg_yyyyzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyzz_1[j];

                    tg_xyyyyzz_xyzzz_0[j] = pb_x * tg_yyyyzz_xyzzz_0[j] + fr * tg_yyyyzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yzzz_1[j];

                    tg_xyyyyzz_xzzzz_0[j] = pb_x * tg_yyyyzz_xzzzz_0[j] + fr * tg_yyyyzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_zzzz_1[j];

                    tg_xyyyyzz_yyyyy_0[j] = pb_x * tg_yyyyzz_yyyyy_0[j] + fr * tg_yyyyzz_yyyyy_1[j];

                    tg_xyyyyzz_yyyyz_0[j] = pb_x * tg_yyyyzz_yyyyz_0[j] + fr * tg_yyyyzz_yyyyz_1[j];

                    tg_xyyyyzz_yyyzz_0[j] = pb_x * tg_yyyyzz_yyyzz_0[j] + fr * tg_yyyyzz_yyyzz_1[j];

                    tg_xyyyyzz_yyzzz_0[j] = pb_x * tg_yyyyzz_yyzzz_0[j] + fr * tg_yyyyzz_yyzzz_1[j];

                    tg_xyyyyzz_yzzzz_0[j] = pb_x * tg_yyyyzz_yzzzz_0[j] + fr * tg_yyyyzz_yzzzz_1[j];

                    tg_xyyyyzz_zzzzz_0[j] = pb_x * tg_yyyyzz_zzzzz_0[j] + fr * tg_yyyyzz_zzzzz_1[j];

                    tg_xyyyzzz_xxxxx_0[j] = pb_x * tg_yyyzzz_xxxxx_0[j] + fr * tg_yyyzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyzzz_xxxx_1[j];

                    tg_xyyyzzz_xxxxy_0[j] = pb_x * tg_yyyzzz_xxxxy_0[j] + fr * tg_yyyzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyzzz_xxxy_1[j];

                    tg_xyyyzzz_xxxxz_0[j] = pb_x * tg_yyyzzz_xxxxz_0[j] + fr * tg_yyyzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyzzz_xxxz_1[j];

                    tg_xyyyzzz_xxxyy_0[j] = pb_x * tg_yyyzzz_xxxyy_0[j] + fr * tg_yyyzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxyy_1[j];

                    tg_xyyyzzz_xxxyz_0[j] = pb_x * tg_yyyzzz_xxxyz_0[j] + fr * tg_yyyzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxyz_1[j];

                    tg_xyyyzzz_xxxzz_0[j] = pb_x * tg_yyyzzz_xxxzz_0[j] + fr * tg_yyyzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxzz_1[j];

                    tg_xyyyzzz_xxyyy_0[j] = pb_x * tg_yyyzzz_xxyyy_0[j] + fr * tg_yyyzzz_xxyyy_1[j] + fl1_fxn * tg_yyyzzz_xyyy_1[j];

                    tg_xyyyzzz_xxyyz_0[j] = pb_x * tg_yyyzzz_xxyyz_0[j] + fr * tg_yyyzzz_xxyyz_1[j] + fl1_fxn * tg_yyyzzz_xyyz_1[j];

                    tg_xyyyzzz_xxyzz_0[j] = pb_x * tg_yyyzzz_xxyzz_0[j] + fr * tg_yyyzzz_xxyzz_1[j] + fl1_fxn * tg_yyyzzz_xyzz_1[j];

                    tg_xyyyzzz_xxzzz_0[j] = pb_x * tg_yyyzzz_xxzzz_0[j] + fr * tg_yyyzzz_xxzzz_1[j] + fl1_fxn * tg_yyyzzz_xzzz_1[j];

                    tg_xyyyzzz_xyyyy_0[j] = pb_x * tg_yyyzzz_xyyyy_0[j] + fr * tg_yyyzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyyy_1[j];

                    tg_xyyyzzz_xyyyz_0[j] = pb_x * tg_yyyzzz_xyyyz_0[j] + fr * tg_yyyzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyyz_1[j];

                    tg_xyyyzzz_xyyzz_0[j] = pb_x * tg_yyyzzz_xyyzz_0[j] + fr * tg_yyyzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyzz_1[j];

                    tg_xyyyzzz_xyzzz_0[j] = pb_x * tg_yyyzzz_xyzzz_0[j] + fr * tg_yyyzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yzzz_1[j];

                    tg_xyyyzzz_xzzzz_0[j] = pb_x * tg_yyyzzz_xzzzz_0[j] + fr * tg_yyyzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_zzzz_1[j];

                    tg_xyyyzzz_yyyyy_0[j] = pb_x * tg_yyyzzz_yyyyy_0[j] + fr * tg_yyyzzz_yyyyy_1[j];

                    tg_xyyyzzz_yyyyz_0[j] = pb_x * tg_yyyzzz_yyyyz_0[j] + fr * tg_yyyzzz_yyyyz_1[j];

                    tg_xyyyzzz_yyyzz_0[j] = pb_x * tg_yyyzzz_yyyzz_0[j] + fr * tg_yyyzzz_yyyzz_1[j];

                    tg_xyyyzzz_yyzzz_0[j] = pb_x * tg_yyyzzz_yyzzz_0[j] + fr * tg_yyyzzz_yyzzz_1[j];

                    tg_xyyyzzz_yzzzz_0[j] = pb_x * tg_yyyzzz_yzzzz_0[j] + fr * tg_yyyzzz_yzzzz_1[j];

                    tg_xyyyzzz_zzzzz_0[j] = pb_x * tg_yyyzzz_zzzzz_0[j] + fr * tg_yyyzzz_zzzzz_1[j];

                    tg_xyyzzzz_xxxxx_0[j] = pb_x * tg_yyzzzz_xxxxx_0[j] + fr * tg_yyzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyzzzz_xxxx_1[j];

                    tg_xyyzzzz_xxxxy_0[j] = pb_x * tg_yyzzzz_xxxxy_0[j] + fr * tg_yyzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyzzzz_xxxy_1[j];

                    tg_xyyzzzz_xxxxz_0[j] = pb_x * tg_yyzzzz_xxxxz_0[j] + fr * tg_yyzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyzzzz_xxxz_1[j];

                    tg_xyyzzzz_xxxyy_0[j] = pb_x * tg_yyzzzz_xxxyy_0[j] + fr * tg_yyzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxyy_1[j];

                    tg_xyyzzzz_xxxyz_0[j] = pb_x * tg_yyzzzz_xxxyz_0[j] + fr * tg_yyzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxyz_1[j];

                    tg_xyyzzzz_xxxzz_0[j] = pb_x * tg_yyzzzz_xxxzz_0[j] + fr * tg_yyzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxzz_1[j];

                    tg_xyyzzzz_xxyyy_0[j] = pb_x * tg_yyzzzz_xxyyy_0[j] + fr * tg_yyzzzz_xxyyy_1[j] + fl1_fxn * tg_yyzzzz_xyyy_1[j];

                    tg_xyyzzzz_xxyyz_0[j] = pb_x * tg_yyzzzz_xxyyz_0[j] + fr * tg_yyzzzz_xxyyz_1[j] + fl1_fxn * tg_yyzzzz_xyyz_1[j];

                    tg_xyyzzzz_xxyzz_0[j] = pb_x * tg_yyzzzz_xxyzz_0[j] + fr * tg_yyzzzz_xxyzz_1[j] + fl1_fxn * tg_yyzzzz_xyzz_1[j];

                    tg_xyyzzzz_xxzzz_0[j] = pb_x * tg_yyzzzz_xxzzz_0[j] + fr * tg_yyzzzz_xxzzz_1[j] + fl1_fxn * tg_yyzzzz_xzzz_1[j];

                    tg_xyyzzzz_xyyyy_0[j] = pb_x * tg_yyzzzz_xyyyy_0[j] + fr * tg_yyzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyyy_1[j];

                    tg_xyyzzzz_xyyyz_0[j] = pb_x * tg_yyzzzz_xyyyz_0[j] + fr * tg_yyzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyyz_1[j];

                    tg_xyyzzzz_xyyzz_0[j] = pb_x * tg_yyzzzz_xyyzz_0[j] + fr * tg_yyzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyzz_1[j];

                    tg_xyyzzzz_xyzzz_0[j] = pb_x * tg_yyzzzz_xyzzz_0[j] + fr * tg_yyzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yzzz_1[j];

                    tg_xyyzzzz_xzzzz_0[j] = pb_x * tg_yyzzzz_xzzzz_0[j] + fr * tg_yyzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_zzzz_1[j];

                    tg_xyyzzzz_yyyyy_0[j] = pb_x * tg_yyzzzz_yyyyy_0[j] + fr * tg_yyzzzz_yyyyy_1[j];

                    tg_xyyzzzz_yyyyz_0[j] = pb_x * tg_yyzzzz_yyyyz_0[j] + fr * tg_yyzzzz_yyyyz_1[j];

                    tg_xyyzzzz_yyyzz_0[j] = pb_x * tg_yyzzzz_yyyzz_0[j] + fr * tg_yyzzzz_yyyzz_1[j];

                    tg_xyyzzzz_yyzzz_0[j] = pb_x * tg_yyzzzz_yyzzz_0[j] + fr * tg_yyzzzz_yyzzz_1[j];

                    tg_xyyzzzz_yzzzz_0[j] = pb_x * tg_yyzzzz_yzzzz_0[j] + fr * tg_yyzzzz_yzzzz_1[j];

                    tg_xyyzzzz_zzzzz_0[j] = pb_x * tg_yyzzzz_zzzzz_0[j] + fr * tg_yyzzzz_zzzzz_1[j];

                    tg_xyzzzzz_xxxxx_0[j] = pb_x * tg_yzzzzz_xxxxx_0[j] + fr * tg_yzzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yzzzzz_xxxx_1[j];

                    tg_xyzzzzz_xxxxy_0[j] = pb_x * tg_yzzzzz_xxxxy_0[j] + fr * tg_yzzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yzzzzz_xxxy_1[j];

                    tg_xyzzzzz_xxxxz_0[j] = pb_x * tg_yzzzzz_xxxxz_0[j] + fr * tg_yzzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yzzzzz_xxxz_1[j];

                    tg_xyzzzzz_xxxyy_0[j] = pb_x * tg_yzzzzz_xxxyy_0[j] + fr * tg_yzzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxyy_1[j];

                    tg_xyzzzzz_xxxyz_0[j] = pb_x * tg_yzzzzz_xxxyz_0[j] + fr * tg_yzzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxyz_1[j];

                    tg_xyzzzzz_xxxzz_0[j] = pb_x * tg_yzzzzz_xxxzz_0[j] + fr * tg_yzzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxzz_1[j];

                    tg_xyzzzzz_xxyyy_0[j] = pb_x * tg_yzzzzz_xxyyy_0[j] + fr * tg_yzzzzz_xxyyy_1[j] + fl1_fxn * tg_yzzzzz_xyyy_1[j];

                    tg_xyzzzzz_xxyyz_0[j] = pb_x * tg_yzzzzz_xxyyz_0[j] + fr * tg_yzzzzz_xxyyz_1[j] + fl1_fxn * tg_yzzzzz_xyyz_1[j];

                    tg_xyzzzzz_xxyzz_0[j] = pb_x * tg_yzzzzz_xxyzz_0[j] + fr * tg_yzzzzz_xxyzz_1[j] + fl1_fxn * tg_yzzzzz_xyzz_1[j];

                    tg_xyzzzzz_xxzzz_0[j] = pb_x * tg_yzzzzz_xxzzz_0[j] + fr * tg_yzzzzz_xxzzz_1[j] + fl1_fxn * tg_yzzzzz_xzzz_1[j];

                    tg_xyzzzzz_xyyyy_0[j] = pb_x * tg_yzzzzz_xyyyy_0[j] + fr * tg_yzzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyyy_1[j];

                    tg_xyzzzzz_xyyyz_0[j] = pb_x * tg_yzzzzz_xyyyz_0[j] + fr * tg_yzzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyyz_1[j];

                    tg_xyzzzzz_xyyzz_0[j] = pb_x * tg_yzzzzz_xyyzz_0[j] + fr * tg_yzzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyzz_1[j];

                    tg_xyzzzzz_xyzzz_0[j] = pb_x * tg_yzzzzz_xyzzz_0[j] + fr * tg_yzzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yzzz_1[j];

                    tg_xyzzzzz_xzzzz_0[j] = pb_x * tg_yzzzzz_xzzzz_0[j] + fr * tg_yzzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_zzzz_1[j];

                    tg_xyzzzzz_yyyyy_0[j] = pb_x * tg_yzzzzz_yyyyy_0[j] + fr * tg_yzzzzz_yyyyy_1[j];

                    tg_xyzzzzz_yyyyz_0[j] = pb_x * tg_yzzzzz_yyyyz_0[j] + fr * tg_yzzzzz_yyyyz_1[j];

                    tg_xyzzzzz_yyyzz_0[j] = pb_x * tg_yzzzzz_yyyzz_0[j] + fr * tg_yzzzzz_yyyzz_1[j];

                    tg_xyzzzzz_yyzzz_0[j] = pb_x * tg_yzzzzz_yyzzz_0[j] + fr * tg_yzzzzz_yyzzz_1[j];

                    tg_xyzzzzz_yzzzz_0[j] = pb_x * tg_yzzzzz_yzzzz_0[j] + fr * tg_yzzzzz_yzzzz_1[j];

                    tg_xyzzzzz_zzzzz_0[j] = pb_x * tg_yzzzzz_zzzzz_0[j] + fr * tg_yzzzzz_zzzzz_1[j];

                    tg_xzzzzzz_xxxxx_0[j] = pb_x * tg_zzzzzz_xxxxx_0[j] + fr * tg_zzzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_zzzzzz_xxxx_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSH_568_662(      CMemBlock2D<double>* primBuffer,
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

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_yyyyyy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 441); 

                auto tg_yyyyyy_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 442); 

                auto tg_yyyyyy_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 443); 

                auto tg_yyyyyy_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 444); 

                auto tg_yyyyyy_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 445); 

                auto tg_yyyyyy_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 446); 

                auto tg_yyyyyy_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 447); 

                auto tg_yyyyyy_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 448); 

                auto tg_yyyyyy_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 449); 

                auto tg_yyyyyy_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 450); 

                auto tg_yyyyyy_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 451); 

                auto tg_yyyyyy_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 452); 

                auto tg_yyyyyy_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 453); 

                auto tg_yyyyyy_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 454); 

                auto tg_yyyyyy_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 455); 

                auto tg_yyyyyy_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 456); 

                auto tg_yyyyyy_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 457); 

                auto tg_yyyyyy_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 458); 

                auto tg_yyyyyy_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 459); 

                auto tg_yyyyyy_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 460); 

                auto tg_yyyyyy_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 461); 

                auto tg_yyyyyz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 462); 

                auto tg_yyyyyz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 463); 

                auto tg_yyyyyz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 464); 

                auto tg_yyyyyz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 465); 

                auto tg_yyyyyz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 466); 

                auto tg_yyyyyz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 467); 

                auto tg_yyyyyz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 468); 

                auto tg_yyyyyz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 469); 

                auto tg_yyyyyz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 470); 

                auto tg_yyyyyz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 471); 

                auto tg_yyyyyz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 472); 

                auto tg_yyyyyz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 473); 

                auto tg_yyyyyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 474); 

                auto tg_yyyyyz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 475); 

                auto tg_yyyyyz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 476); 

                auto tg_yyyyyz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 477); 

                auto tg_yyyyyz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 478); 

                auto tg_yyyyyz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 479); 

                auto tg_yyyyyz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 480); 

                auto tg_yyyyyz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 481); 

                auto tg_yyyyyz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 482); 

                auto tg_yyyyzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 483); 

                auto tg_yyyyzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 484); 

                auto tg_yyyyzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 485); 

                auto tg_yyyyzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 486); 

                auto tg_yyyyzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 487); 

                auto tg_yyyyzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 488); 

                auto tg_yyyyzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 489); 

                auto tg_yyyyzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 490); 

                auto tg_yyyyzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 491); 

                auto tg_yyyyzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 492); 

                auto tg_yyyyzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 493); 

                auto tg_yyyyzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 494); 

                auto tg_yyyyzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 495); 

                auto tg_yyyyzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 496); 

                auto tg_yyyyzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 497); 

                auto tg_yyyyzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 498); 

                auto tg_yyyyzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 499); 

                auto tg_yyyyzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 500); 

                auto tg_yyyyzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 501); 

                auto tg_yyyyzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 502); 

                auto tg_yyyyzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 503); 

                auto tg_yyyzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 504); 

                auto tg_yyyzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 505); 

                auto tg_yyyzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 506); 

                auto tg_yyyzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 507); 

                auto tg_yyyzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 508); 

                auto tg_yyyzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 509); 

                auto tg_yyyzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 510); 

                auto tg_yyyzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 511); 

                auto tg_yyyzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 512); 

                auto tg_yyyzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 513); 

                auto tg_yyyzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 514); 

                auto tg_zzzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 568); 

                auto tg_zzzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 569); 

                auto tg_zzzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 570); 

                auto tg_zzzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 571); 

                auto tg_zzzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 572); 

                auto tg_zzzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 573); 

                auto tg_zzzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 574); 

                auto tg_zzzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 575); 

                auto tg_zzzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 576); 

                auto tg_zzzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 577); 

                auto tg_zzzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 578); 

                auto tg_zzzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 579); 

                auto tg_zzzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 580); 

                auto tg_zzzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 581); 

                auto tg_zzzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 582); 

                auto tg_zzzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 583); 

                auto tg_zzzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 584); 

                auto tg_zzzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 585); 

                auto tg_zzzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 586); 

                auto tg_zzzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 587); 

                auto tg_yyyyyy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 441); 

                auto tg_yyyyyy_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 442); 

                auto tg_yyyyyy_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 443); 

                auto tg_yyyyyy_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 444); 

                auto tg_yyyyyy_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 445); 

                auto tg_yyyyyy_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 446); 

                auto tg_yyyyyy_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 447); 

                auto tg_yyyyyy_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 448); 

                auto tg_yyyyyy_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 449); 

                auto tg_yyyyyy_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 450); 

                auto tg_yyyyyy_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 451); 

                auto tg_yyyyyy_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 452); 

                auto tg_yyyyyy_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 453); 

                auto tg_yyyyyy_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 454); 

                auto tg_yyyyyy_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 455); 

                auto tg_yyyyyy_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 456); 

                auto tg_yyyyyy_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 457); 

                auto tg_yyyyyy_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 458); 

                auto tg_yyyyyy_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 459); 

                auto tg_yyyyyy_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 460); 

                auto tg_yyyyyy_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 461); 

                auto tg_yyyyyz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 462); 

                auto tg_yyyyyz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 463); 

                auto tg_yyyyyz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 464); 

                auto tg_yyyyyz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 465); 

                auto tg_yyyyyz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 466); 

                auto tg_yyyyyz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 467); 

                auto tg_yyyyyz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 468); 

                auto tg_yyyyyz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 469); 

                auto tg_yyyyyz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 470); 

                auto tg_yyyyyz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 471); 

                auto tg_yyyyyz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 472); 

                auto tg_yyyyyz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 473); 

                auto tg_yyyyyz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 474); 

                auto tg_yyyyyz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 475); 

                auto tg_yyyyyz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 476); 

                auto tg_yyyyyz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 477); 

                auto tg_yyyyyz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 478); 

                auto tg_yyyyyz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 479); 

                auto tg_yyyyyz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 480); 

                auto tg_yyyyyz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 481); 

                auto tg_yyyyyz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 482); 

                auto tg_yyyyzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 483); 

                auto tg_yyyyzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 484); 

                auto tg_yyyyzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 485); 

                auto tg_yyyyzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 486); 

                auto tg_yyyyzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 487); 

                auto tg_yyyyzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 488); 

                auto tg_yyyyzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 489); 

                auto tg_yyyyzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 490); 

                auto tg_yyyyzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 491); 

                auto tg_yyyyzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 492); 

                auto tg_yyyyzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 493); 

                auto tg_yyyyzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 494); 

                auto tg_yyyyzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 495); 

                auto tg_yyyyzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 496); 

                auto tg_yyyyzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 497); 

                auto tg_yyyyzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 498); 

                auto tg_yyyyzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 499); 

                auto tg_yyyyzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 500); 

                auto tg_yyyyzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 501); 

                auto tg_yyyyzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 502); 

                auto tg_yyyyzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 503); 

                auto tg_yyyzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 504); 

                auto tg_yyyzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 505); 

                auto tg_yyyzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 506); 

                auto tg_yyyzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 507); 

                auto tg_yyyzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 508); 

                auto tg_yyyzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 509); 

                auto tg_yyyzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 510); 

                auto tg_yyyzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 511); 

                auto tg_yyyzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 512); 

                auto tg_yyyzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 513); 

                auto tg_yyyzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 514); 

                auto tg_zzzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 568); 

                auto tg_zzzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 569); 

                auto tg_zzzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 570); 

                auto tg_zzzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 571); 

                auto tg_zzzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 572); 

                auto tg_zzzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 573); 

                auto tg_zzzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 574); 

                auto tg_zzzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 575); 

                auto tg_zzzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 576); 

                auto tg_zzzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 577); 

                auto tg_zzzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 578); 

                auto tg_zzzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 579); 

                auto tg_zzzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 580); 

                auto tg_zzzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 581); 

                auto tg_zzzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 582); 

                auto tg_zzzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 583); 

                auto tg_zzzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 584); 

                auto tg_zzzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 585); 

                auto tg_zzzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 586); 

                auto tg_zzzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 587); 

                auto tg_yyyyy_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 315); 

                auto tg_yyyyy_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 316); 

                auto tg_yyyyy_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 317); 

                auto tg_yyyyy_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 318); 

                auto tg_yyyyy_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 319); 

                auto tg_yyyyy_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 320); 

                auto tg_yyyyy_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 321); 

                auto tg_yyyyy_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 322); 

                auto tg_yyyyy_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 323); 

                auto tg_yyyyy_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 324); 

                auto tg_yyyyy_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 325); 

                auto tg_yyyyy_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 326); 

                auto tg_yyyyy_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 327); 

                auto tg_yyyyy_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 328); 

                auto tg_yyyyy_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 329); 

                auto tg_yyyyy_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 330); 

                auto tg_yyyyy_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 331); 

                auto tg_yyyyy_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 332); 

                auto tg_yyyyy_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 333); 

                auto tg_yyyyy_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 334); 

                auto tg_yyyyy_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 335); 

                auto tg_yyyyz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 336); 

                auto tg_yyyyz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 337); 

                auto tg_yyyyz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 338); 

                auto tg_yyyyz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 339); 

                auto tg_yyyyz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 340); 

                auto tg_yyyyz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 341); 

                auto tg_yyyyz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 342); 

                auto tg_yyyyz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 343); 

                auto tg_yyyyz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 344); 

                auto tg_yyyyz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 345); 

                auto tg_yyyyz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 346); 

                auto tg_yyyyz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 347); 

                auto tg_yyyyz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 348); 

                auto tg_yyyyz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 349); 

                auto tg_yyyyz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 350); 

                auto tg_yyyyz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 351); 

                auto tg_yyyyz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 352); 

                auto tg_yyyyz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 353); 

                auto tg_yyyyz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 354); 

                auto tg_yyyyz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 355); 

                auto tg_yyyyz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 356); 

                auto tg_yyyzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 357); 

                auto tg_yyyzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 358); 

                auto tg_yyyzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 359); 

                auto tg_yyyzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 360); 

                auto tg_yyyzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 361); 

                auto tg_yyyzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 362); 

                auto tg_yyyzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 363); 

                auto tg_yyyzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 364); 

                auto tg_yyyzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 365); 

                auto tg_yyyzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 366); 

                auto tg_yyyzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 367); 

                auto tg_yyyzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 368); 

                auto tg_yyyzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 369); 

                auto tg_yyyzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 370); 

                auto tg_yyyzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 371); 

                auto tg_yyyzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 372); 

                auto tg_yyyzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 373); 

                auto tg_yyyzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 374); 

                auto tg_yyyzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 375); 

                auto tg_yyyzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 376); 

                auto tg_yyyzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 377); 

                auto tg_yyzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 378); 

                auto tg_yyzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 379); 

                auto tg_yyzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 380); 

                auto tg_yyzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 381); 

                auto tg_yyzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 382); 

                auto tg_yyzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 383); 

                auto tg_yyzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 384); 

                auto tg_yyzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 385); 

                auto tg_yyzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 386); 

                auto tg_yyzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 387); 

                auto tg_yyzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 388); 

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

                auto tg_yyzzz_xyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 388); 

                auto tg_yyyyyy_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 315); 

                auto tg_yyyyyy_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 316); 

                auto tg_yyyyyy_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 317); 

                auto tg_yyyyyy_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 318); 

                auto tg_yyyyyy_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 319); 

                auto tg_yyyyyy_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 320); 

                auto tg_yyyyyy_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 321); 

                auto tg_yyyyyy_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 322); 

                auto tg_yyyyyy_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 323); 

                auto tg_yyyyyy_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 324); 

                auto tg_yyyyyy_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 325); 

                auto tg_yyyyyy_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 326); 

                auto tg_yyyyyy_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 327); 

                auto tg_yyyyyy_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 328); 

                auto tg_yyyyyy_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 329); 

                auto tg_yyyyyz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 330); 

                auto tg_yyyyyz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 331); 

                auto tg_yyyyyz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 332); 

                auto tg_yyyyyz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 333); 

                auto tg_yyyyyz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 334); 

                auto tg_yyyyyz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 335); 

                auto tg_yyyyyz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 336); 

                auto tg_yyyyyz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 337); 

                auto tg_yyyyyz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 338); 

                auto tg_yyyyyz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 339); 

                auto tg_yyyyyz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 340); 

                auto tg_yyyyyz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 341); 

                auto tg_yyyyyz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 342); 

                auto tg_yyyyyz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 343); 

                auto tg_yyyyyz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 344); 

                auto tg_yyyyzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 345); 

                auto tg_yyyyzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 346); 

                auto tg_yyyyzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 347); 

                auto tg_yyyyzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 348); 

                auto tg_yyyyzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 349); 

                auto tg_yyyyzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 350); 

                auto tg_yyyyzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 351); 

                auto tg_yyyyzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 352); 

                auto tg_yyyyzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 353); 

                auto tg_yyyyzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 354); 

                auto tg_yyyyzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 355); 

                auto tg_yyyyzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 356); 

                auto tg_yyyyzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 357); 

                auto tg_yyyyzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 358); 

                auto tg_yyyyzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 359); 

                auto tg_yyyzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 360); 

                auto tg_yyyzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 361); 

                auto tg_yyyzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 362); 

                auto tg_yyyzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 363); 

                auto tg_yyyzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 364); 

                auto tg_yyyzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 365); 

                auto tg_yyyzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 366); 

                auto tg_zzzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 406); 

                auto tg_zzzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 407); 

                auto tg_zzzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 408); 

                auto tg_zzzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 409); 

                auto tg_zzzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 410); 

                auto tg_zzzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 411); 

                auto tg_zzzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 412); 

                auto tg_zzzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 413); 

                auto tg_zzzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 414); 

                auto tg_zzzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 415); 

                auto tg_zzzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 416); 

                auto tg_zzzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 417); 

                auto tg_zzzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 418); 

                auto tg_zzzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 419); 

                // set up pointers to integrals

                auto tg_xzzzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 568); 

                auto tg_xzzzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 569); 

                auto tg_xzzzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 570); 

                auto tg_xzzzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 571); 

                auto tg_xzzzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 572); 

                auto tg_xzzzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 573); 

                auto tg_xzzzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 574); 

                auto tg_xzzzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 575); 

                auto tg_xzzzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 576); 

                auto tg_xzzzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 577); 

                auto tg_xzzzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 578); 

                auto tg_xzzzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 579); 

                auto tg_xzzzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 580); 

                auto tg_xzzzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 581); 

                auto tg_xzzzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 582); 

                auto tg_xzzzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 583); 

                auto tg_xzzzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 584); 

                auto tg_xzzzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 585); 

                auto tg_xzzzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 586); 

                auto tg_xzzzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 587); 

                auto tg_yyyyyyy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 588); 

                auto tg_yyyyyyy_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 589); 

                auto tg_yyyyyyy_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 590); 

                auto tg_yyyyyyy_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 591); 

                auto tg_yyyyyyy_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 592); 

                auto tg_yyyyyyy_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 593); 

                auto tg_yyyyyyy_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 594); 

                auto tg_yyyyyyy_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 595); 

                auto tg_yyyyyyy_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 596); 

                auto tg_yyyyyyy_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 597); 

                auto tg_yyyyyyy_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 598); 

                auto tg_yyyyyyy_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 599); 

                auto tg_yyyyyyy_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 600); 

                auto tg_yyyyyyy_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 601); 

                auto tg_yyyyyyy_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 602); 

                auto tg_yyyyyyy_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 603); 

                auto tg_yyyyyyy_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 604); 

                auto tg_yyyyyyy_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 605); 

                auto tg_yyyyyyy_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 606); 

                auto tg_yyyyyyy_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 607); 

                auto tg_yyyyyyy_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 608); 

                auto tg_yyyyyyz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 609); 

                auto tg_yyyyyyz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 610); 

                auto tg_yyyyyyz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 611); 

                auto tg_yyyyyyz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 612); 

                auto tg_yyyyyyz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 613); 

                auto tg_yyyyyyz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 614); 

                auto tg_yyyyyyz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 615); 

                auto tg_yyyyyyz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 616); 

                auto tg_yyyyyyz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 617); 

                auto tg_yyyyyyz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 618); 

                auto tg_yyyyyyz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 619); 

                auto tg_yyyyyyz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 620); 

                auto tg_yyyyyyz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 621); 

                auto tg_yyyyyyz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 622); 

                auto tg_yyyyyyz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 623); 

                auto tg_yyyyyyz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 624); 

                auto tg_yyyyyyz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 625); 

                auto tg_yyyyyyz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 626); 

                auto tg_yyyyyyz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 627); 

                auto tg_yyyyyyz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 628); 

                auto tg_yyyyyyz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 629); 

                auto tg_yyyyyzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 630); 

                auto tg_yyyyyzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 631); 

                auto tg_yyyyyzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 632); 

                auto tg_yyyyyzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 633); 

                auto tg_yyyyyzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 634); 

                auto tg_yyyyyzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 635); 

                auto tg_yyyyyzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 636); 

                auto tg_yyyyyzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 637); 

                auto tg_yyyyyzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 638); 

                auto tg_yyyyyzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 639); 

                auto tg_yyyyyzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 640); 

                auto tg_yyyyyzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 641); 

                auto tg_yyyyyzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 642); 

                auto tg_yyyyyzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 643); 

                auto tg_yyyyyzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 644); 

                auto tg_yyyyyzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 645); 

                auto tg_yyyyyzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 646); 

                auto tg_yyyyyzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 647); 

                auto tg_yyyyyzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 648); 

                auto tg_yyyyyzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 649); 

                auto tg_yyyyyzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 650); 

                auto tg_yyyyzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 651); 

                auto tg_yyyyzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 652); 

                auto tg_yyyyzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 653); 

                auto tg_yyyyzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 654); 

                auto tg_yyyyzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 655); 

                auto tg_yyyyzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 656); 

                auto tg_yyyyzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 657); 

                auto tg_yyyyzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 658); 

                auto tg_yyyyzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 659); 

                auto tg_yyyyzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 660); 

                auto tg_yyyyzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 661); 

                // Batch of Integrals (568,662)

                #pragma omp simd aligned(fxn, fza, tg_xzzzzzz_xxxxy_0, tg_xzzzzzz_xxxxz_0, tg_xzzzzzz_xxxyy_0, \
                                         tg_xzzzzzz_xxxyz_0, tg_xzzzzzz_xxxzz_0, tg_xzzzzzz_xxyyy_0, tg_xzzzzzz_xxyyz_0, \
                                         tg_xzzzzzz_xxyzz_0, tg_xzzzzzz_xxzzz_0, tg_xzzzzzz_xyyyy_0, tg_xzzzzzz_xyyyz_0, \
                                         tg_xzzzzzz_xyyzz_0, tg_xzzzzzz_xyzzz_0, tg_xzzzzzz_xzzzz_0, tg_xzzzzzz_yyyyy_0, \
                                         tg_xzzzzzz_yyyyz_0, tg_xzzzzzz_yyyzz_0, tg_xzzzzzz_yyzzz_0, tg_xzzzzzz_yzzzz_0, \
                                         tg_xzzzzzz_zzzzz_0, tg_yyyyy_xxxxx_0, tg_yyyyy_xxxxx_1, tg_yyyyy_xxxxy_0, \
                                         tg_yyyyy_xxxxy_1, tg_yyyyy_xxxxz_0, tg_yyyyy_xxxxz_1, tg_yyyyy_xxxyy_0, \
                                         tg_yyyyy_xxxyy_1, tg_yyyyy_xxxyz_0, tg_yyyyy_xxxyz_1, tg_yyyyy_xxxzz_0, \
                                         tg_yyyyy_xxxzz_1, tg_yyyyy_xxyyy_0, tg_yyyyy_xxyyy_1, tg_yyyyy_xxyyz_0, \
                                         tg_yyyyy_xxyyz_1, tg_yyyyy_xxyzz_0, tg_yyyyy_xxyzz_1, tg_yyyyy_xxzzz_0, \
                                         tg_yyyyy_xxzzz_1, tg_yyyyy_xyyyy_0, tg_yyyyy_xyyyy_1, tg_yyyyy_xyyyz_0, \
                                         tg_yyyyy_xyyyz_1, tg_yyyyy_xyyzz_0, tg_yyyyy_xyyzz_1, tg_yyyyy_xyzzz_0, \
                                         tg_yyyyy_xyzzz_1, tg_yyyyy_xzzzz_0, tg_yyyyy_xzzzz_1, tg_yyyyy_yyyyy_0, \
                                         tg_yyyyy_yyyyy_1, tg_yyyyy_yyyyz_0, tg_yyyyy_yyyyz_1, tg_yyyyy_yyyzz_0, \
                                         tg_yyyyy_yyyzz_1, tg_yyyyy_yyzzz_0, tg_yyyyy_yyzzz_1, tg_yyyyy_yzzzz_0, \
                                         tg_yyyyy_yzzzz_1, tg_yyyyy_zzzzz_0, tg_yyyyy_zzzzz_1, tg_yyyyyy_xxxx_1, \
                                         tg_yyyyyy_xxxxx_0, tg_yyyyyy_xxxxx_1, tg_yyyyyy_xxxxy_0, tg_yyyyyy_xxxxy_1, \
                                         tg_yyyyyy_xxxxz_0, tg_yyyyyy_xxxxz_1, tg_yyyyyy_xxxy_1, tg_yyyyyy_xxxyy_0, \
                                         tg_yyyyyy_xxxyy_1, tg_yyyyyy_xxxyz_0, tg_yyyyyy_xxxyz_1, tg_yyyyyy_xxxz_1, \
                                         tg_yyyyyy_xxxzz_0, tg_yyyyyy_xxxzz_1, tg_yyyyyy_xxyy_1, tg_yyyyyy_xxyyy_0, \
                                         tg_yyyyyy_xxyyy_1, tg_yyyyyy_xxyyz_0, tg_yyyyyy_xxyyz_1, tg_yyyyyy_xxyz_1, \
                                         tg_yyyyyy_xxyzz_0, tg_yyyyyy_xxyzz_1, tg_yyyyyy_xxzz_1, tg_yyyyyy_xxzzz_0, \
                                         tg_yyyyyy_xxzzz_1, tg_yyyyyy_xyyy_1, tg_yyyyyy_xyyyy_0, tg_yyyyyy_xyyyy_1, \
                                         tg_yyyyyy_xyyyz_0, tg_yyyyyy_xyyyz_1, tg_yyyyyy_xyyz_1, tg_yyyyyy_xyyzz_0, \
                                         tg_yyyyyy_xyyzz_1, tg_yyyyyy_xyzz_1, tg_yyyyyy_xyzzz_0, tg_yyyyyy_xyzzz_1, \
                                         tg_yyyyyy_xzzz_1, tg_yyyyyy_xzzzz_0, tg_yyyyyy_xzzzz_1, tg_yyyyyy_yyyy_1, \
                                         tg_yyyyyy_yyyyy_0, tg_yyyyyy_yyyyy_1, tg_yyyyyy_yyyyz_0, tg_yyyyyy_yyyyz_1, \
                                         tg_yyyyyy_yyyz_1, tg_yyyyyy_yyyzz_0, tg_yyyyyy_yyyzz_1, tg_yyyyyy_yyzz_1, \
                                         tg_yyyyyy_yyzzz_0, tg_yyyyyy_yyzzz_1, tg_yyyyyy_yzzz_1, tg_yyyyyy_yzzzz_0, \
                                         tg_yyyyyy_yzzzz_1, tg_yyyyyy_zzzz_1, tg_yyyyyy_zzzzz_0, tg_yyyyyy_zzzzz_1, \
                                         tg_yyyyyyy_xxxxx_0, tg_yyyyyyy_xxxxy_0, tg_yyyyyyy_xxxxz_0, tg_yyyyyyy_xxxyy_0, \
                                         tg_yyyyyyy_xxxyz_0, tg_yyyyyyy_xxxzz_0, tg_yyyyyyy_xxyyy_0, tg_yyyyyyy_xxyyz_0, \
                                         tg_yyyyyyy_xxyzz_0, tg_yyyyyyy_xxzzz_0, tg_yyyyyyy_xyyyy_0, tg_yyyyyyy_xyyyz_0, \
                                         tg_yyyyyyy_xyyzz_0, tg_yyyyyyy_xyzzz_0, tg_yyyyyyy_xzzzz_0, tg_yyyyyyy_yyyyy_0, \
                                         tg_yyyyyyy_yyyyz_0, tg_yyyyyyy_yyyzz_0, tg_yyyyyyy_yyzzz_0, tg_yyyyyyy_yzzzz_0, \
                                         tg_yyyyyyy_zzzzz_0, tg_yyyyyyz_xxxxx_0, tg_yyyyyyz_xxxxy_0, tg_yyyyyyz_xxxxz_0, \
                                         tg_yyyyyyz_xxxyy_0, tg_yyyyyyz_xxxyz_0, tg_yyyyyyz_xxxzz_0, tg_yyyyyyz_xxyyy_0, \
                                         tg_yyyyyyz_xxyyz_0, tg_yyyyyyz_xxyzz_0, tg_yyyyyyz_xxzzz_0, tg_yyyyyyz_xyyyy_0, \
                                         tg_yyyyyyz_xyyyz_0, tg_yyyyyyz_xyyzz_0, tg_yyyyyyz_xyzzz_0, tg_yyyyyyz_xzzzz_0, \
                                         tg_yyyyyyz_yyyyy_0, tg_yyyyyyz_yyyyz_0, tg_yyyyyyz_yyyzz_0, tg_yyyyyyz_yyzzz_0, \
                                         tg_yyyyyyz_yzzzz_0, tg_yyyyyyz_zzzzz_0, tg_yyyyyz_xxxx_1, tg_yyyyyz_xxxxx_0, \
                                         tg_yyyyyz_xxxxx_1, tg_yyyyyz_xxxxy_0, tg_yyyyyz_xxxxy_1, tg_yyyyyz_xxxxz_0, \
                                         tg_yyyyyz_xxxxz_1, tg_yyyyyz_xxxy_1, tg_yyyyyz_xxxyy_0, tg_yyyyyz_xxxyy_1, \
                                         tg_yyyyyz_xxxyz_0, tg_yyyyyz_xxxyz_1, tg_yyyyyz_xxxz_1, tg_yyyyyz_xxxzz_0, \
                                         tg_yyyyyz_xxxzz_1, tg_yyyyyz_xxyy_1, tg_yyyyyz_xxyyy_0, tg_yyyyyz_xxyyy_1, \
                                         tg_yyyyyz_xxyyz_0, tg_yyyyyz_xxyyz_1, tg_yyyyyz_xxyz_1, tg_yyyyyz_xxyzz_0, \
                                         tg_yyyyyz_xxyzz_1, tg_yyyyyz_xxzz_1, tg_yyyyyz_xxzzz_0, tg_yyyyyz_xxzzz_1, \
                                         tg_yyyyyz_xyyy_1, tg_yyyyyz_xyyyy_0, tg_yyyyyz_xyyyy_1, tg_yyyyyz_xyyyz_0, \
                                         tg_yyyyyz_xyyyz_1, tg_yyyyyz_xyyz_1, tg_yyyyyz_xyyzz_0, tg_yyyyyz_xyyzz_1, \
                                         tg_yyyyyz_xyzz_1, tg_yyyyyz_xyzzz_0, tg_yyyyyz_xyzzz_1, tg_yyyyyz_xzzz_1, \
                                         tg_yyyyyz_xzzzz_0, tg_yyyyyz_xzzzz_1, tg_yyyyyz_yyyy_1, tg_yyyyyz_yyyyy_0, \
                                         tg_yyyyyz_yyyyy_1, tg_yyyyyz_yyyyz_0, tg_yyyyyz_yyyyz_1, tg_yyyyyz_yyyz_1, \
                                         tg_yyyyyz_yyyzz_0, tg_yyyyyz_yyyzz_1, tg_yyyyyz_yyzz_1, tg_yyyyyz_yyzzz_0, \
                                         tg_yyyyyz_yyzzz_1, tg_yyyyyz_yzzz_1, tg_yyyyyz_yzzzz_0, tg_yyyyyz_yzzzz_1, \
                                         tg_yyyyyz_zzzz_1, tg_yyyyyz_zzzzz_0, tg_yyyyyz_zzzzz_1, tg_yyyyyzz_xxxxx_0, \
                                         tg_yyyyyzz_xxxxy_0, tg_yyyyyzz_xxxxz_0, tg_yyyyyzz_xxxyy_0, tg_yyyyyzz_xxxyz_0, \
                                         tg_yyyyyzz_xxxzz_0, tg_yyyyyzz_xxyyy_0, tg_yyyyyzz_xxyyz_0, tg_yyyyyzz_xxyzz_0, \
                                         tg_yyyyyzz_xxzzz_0, tg_yyyyyzz_xyyyy_0, tg_yyyyyzz_xyyyz_0, tg_yyyyyzz_xyyzz_0, \
                                         tg_yyyyyzz_xyzzz_0, tg_yyyyyzz_xzzzz_0, tg_yyyyyzz_yyyyy_0, tg_yyyyyzz_yyyyz_0, \
                                         tg_yyyyyzz_yyyzz_0, tg_yyyyyzz_yyzzz_0, tg_yyyyyzz_yzzzz_0, tg_yyyyyzz_zzzzz_0, \
                                         tg_yyyyz_xxxxx_0, tg_yyyyz_xxxxx_1, tg_yyyyz_xxxxy_0, tg_yyyyz_xxxxy_1, \
                                         tg_yyyyz_xxxxz_0, tg_yyyyz_xxxxz_1, tg_yyyyz_xxxyy_0, tg_yyyyz_xxxyy_1, \
                                         tg_yyyyz_xxxyz_0, tg_yyyyz_xxxyz_1, tg_yyyyz_xxxzz_0, tg_yyyyz_xxxzz_1, \
                                         tg_yyyyz_xxyyy_0, tg_yyyyz_xxyyy_1, tg_yyyyz_xxyyz_0, tg_yyyyz_xxyyz_1, \
                                         tg_yyyyz_xxyzz_0, tg_yyyyz_xxyzz_1, tg_yyyyz_xxzzz_0, tg_yyyyz_xxzzz_1, \
                                         tg_yyyyz_xyyyy_0, tg_yyyyz_xyyyy_1, tg_yyyyz_xyyyz_0, tg_yyyyz_xyyyz_1, \
                                         tg_yyyyz_xyyzz_0, tg_yyyyz_xyyzz_1, tg_yyyyz_xyzzz_0, tg_yyyyz_xyzzz_1, \
                                         tg_yyyyz_xzzzz_0, tg_yyyyz_xzzzz_1, tg_yyyyz_yyyyy_0, tg_yyyyz_yyyyy_1, \
                                         tg_yyyyz_yyyyz_0, tg_yyyyz_yyyyz_1, tg_yyyyz_yyyzz_0, tg_yyyyz_yyyzz_1, \
                                         tg_yyyyz_yyzzz_0, tg_yyyyz_yyzzz_1, tg_yyyyz_yzzzz_0, tg_yyyyz_yzzzz_1, \
                                         tg_yyyyz_zzzzz_0, tg_yyyyz_zzzzz_1, tg_yyyyzz_xxxx_1, tg_yyyyzz_xxxxx_0, \
                                         tg_yyyyzz_xxxxx_1, tg_yyyyzz_xxxxy_0, tg_yyyyzz_xxxxy_1, tg_yyyyzz_xxxxz_0, \
                                         tg_yyyyzz_xxxxz_1, tg_yyyyzz_xxxy_1, tg_yyyyzz_xxxyy_0, tg_yyyyzz_xxxyy_1, \
                                         tg_yyyyzz_xxxyz_0, tg_yyyyzz_xxxyz_1, tg_yyyyzz_xxxz_1, tg_yyyyzz_xxxzz_0, \
                                         tg_yyyyzz_xxxzz_1, tg_yyyyzz_xxyy_1, tg_yyyyzz_xxyyy_0, tg_yyyyzz_xxyyy_1, \
                                         tg_yyyyzz_xxyyz_0, tg_yyyyzz_xxyyz_1, tg_yyyyzz_xxyz_1, tg_yyyyzz_xxyzz_0, \
                                         tg_yyyyzz_xxyzz_1, tg_yyyyzz_xxzz_1, tg_yyyyzz_xxzzz_0, tg_yyyyzz_xxzzz_1, \
                                         tg_yyyyzz_xyyy_1, tg_yyyyzz_xyyyy_0, tg_yyyyzz_xyyyy_1, tg_yyyyzz_xyyyz_0, \
                                         tg_yyyyzz_xyyyz_1, tg_yyyyzz_xyyz_1, tg_yyyyzz_xyyzz_0, tg_yyyyzz_xyyzz_1, \
                                         tg_yyyyzz_xyzz_1, tg_yyyyzz_xyzzz_0, tg_yyyyzz_xyzzz_1, tg_yyyyzz_xzzz_1, \
                                         tg_yyyyzz_xzzzz_0, tg_yyyyzz_xzzzz_1, tg_yyyyzz_yyyy_1, tg_yyyyzz_yyyyy_0, \
                                         tg_yyyyzz_yyyyy_1, tg_yyyyzz_yyyyz_0, tg_yyyyzz_yyyyz_1, tg_yyyyzz_yyyz_1, \
                                         tg_yyyyzz_yyyzz_0, tg_yyyyzz_yyyzz_1, tg_yyyyzz_yyzz_1, tg_yyyyzz_yyzzz_0, \
                                         tg_yyyyzz_yyzzz_1, tg_yyyyzz_yzzz_1, tg_yyyyzz_yzzzz_0, tg_yyyyzz_yzzzz_1, \
                                         tg_yyyyzz_zzzz_1, tg_yyyyzz_zzzzz_0, tg_yyyyzz_zzzzz_1, tg_yyyyzzz_xxxxx_0, \
                                         tg_yyyyzzz_xxxxy_0, tg_yyyyzzz_xxxxz_0, tg_yyyyzzz_xxxyy_0, tg_yyyyzzz_xxxyz_0, \
                                         tg_yyyyzzz_xxxzz_0, tg_yyyyzzz_xxyyy_0, tg_yyyyzzz_xxyyz_0, tg_yyyyzzz_xxyzz_0, \
                                         tg_yyyyzzz_xxzzz_0, tg_yyyyzzz_xyyyy_0, tg_yyyzz_xxxxx_0, tg_yyyzz_xxxxx_1, \
                                         tg_yyyzz_xxxxy_0, tg_yyyzz_xxxxy_1, tg_yyyzz_xxxxz_0, tg_yyyzz_xxxxz_1, \
                                         tg_yyyzz_xxxyy_0, tg_yyyzz_xxxyy_1, tg_yyyzz_xxxyz_0, tg_yyyzz_xxxyz_1, \
                                         tg_yyyzz_xxxzz_0, tg_yyyzz_xxxzz_1, tg_yyyzz_xxyyy_0, tg_yyyzz_xxyyy_1, \
                                         tg_yyyzz_xxyyz_0, tg_yyyzz_xxyyz_1, tg_yyyzz_xxyzz_0, tg_yyyzz_xxyzz_1, \
                                         tg_yyyzz_xxzzz_0, tg_yyyzz_xxzzz_1, tg_yyyzz_xyyyy_0, tg_yyyzz_xyyyy_1, \
                                         tg_yyyzz_xyyyz_0, tg_yyyzz_xyyyz_1, tg_yyyzz_xyyzz_0, tg_yyyzz_xyyzz_1, \
                                         tg_yyyzz_xyzzz_0, tg_yyyzz_xyzzz_1, tg_yyyzz_xzzzz_0, tg_yyyzz_xzzzz_1, \
                                         tg_yyyzz_yyyyy_0, tg_yyyzz_yyyyy_1, tg_yyyzz_yyyyz_0, tg_yyyzz_yyyyz_1, \
                                         tg_yyyzz_yyyzz_0, tg_yyyzz_yyyzz_1, tg_yyyzz_yyzzz_0, tg_yyyzz_yyzzz_1, \
                                         tg_yyyzz_yzzzz_0, tg_yyyzz_yzzzz_1, tg_yyyzz_zzzzz_0, tg_yyyzz_zzzzz_1, \
                                         tg_yyyzzz_xxxx_1, tg_yyyzzz_xxxxx_0, tg_yyyzzz_xxxxx_1, tg_yyyzzz_xxxxy_0, \
                                         tg_yyyzzz_xxxxy_1, tg_yyyzzz_xxxxz_0, tg_yyyzzz_xxxxz_1, tg_yyyzzz_xxxy_1, \
                                         tg_yyyzzz_xxxyy_0, tg_yyyzzz_xxxyy_1, tg_yyyzzz_xxxyz_0, tg_yyyzzz_xxxyz_1, \
                                         tg_yyyzzz_xxxz_1, tg_yyyzzz_xxxzz_0, tg_yyyzzz_xxxzz_1, tg_yyyzzz_xxyy_1, \
                                         tg_yyyzzz_xxyyy_0, tg_yyyzzz_xxyyy_1, tg_yyyzzz_xxyyz_0, tg_yyyzzz_xxyyz_1, \
                                         tg_yyyzzz_xxyz_1, tg_yyyzzz_xxyzz_0, tg_yyyzzz_xxyzz_1, tg_yyyzzz_xxzz_1, \
                                         tg_yyyzzz_xxzzz_0, tg_yyyzzz_xxzzz_1, tg_yyyzzz_xyyy_1, tg_yyyzzz_xyyyy_0, \
                                         tg_yyyzzz_xyyyy_1, tg_yyzzz_xxxxx_0, tg_yyzzz_xxxxx_1, tg_yyzzz_xxxxy_0, \
                                         tg_yyzzz_xxxxy_1, tg_yyzzz_xxxxz_0, tg_yyzzz_xxxxz_1, tg_yyzzz_xxxyy_0, \
                                         tg_yyzzz_xxxyy_1, tg_yyzzz_xxxyz_0, tg_yyzzz_xxxyz_1, tg_yyzzz_xxxzz_0, \
                                         tg_yyzzz_xxxzz_1, tg_yyzzz_xxyyy_0, tg_yyzzz_xxyyy_1, tg_yyzzz_xxyyz_0, \
                                         tg_yyzzz_xxyyz_1, tg_yyzzz_xxyzz_0, tg_yyzzz_xxyzz_1, tg_yyzzz_xxzzz_0, \
                                         tg_yyzzz_xxzzz_1, tg_yyzzz_xyyyy_0, tg_yyzzz_xyyyy_1, tg_zzzzzz_xxxxy_0, \
                                         tg_zzzzzz_xxxxy_1, tg_zzzzzz_xxxxz_0, tg_zzzzzz_xxxxz_1, tg_zzzzzz_xxxy_1, \
                                         tg_zzzzzz_xxxyy_0, tg_zzzzzz_xxxyy_1, tg_zzzzzz_xxxyz_0, tg_zzzzzz_xxxyz_1, \
                                         tg_zzzzzz_xxxz_1, tg_zzzzzz_xxxzz_0, tg_zzzzzz_xxxzz_1, tg_zzzzzz_xxyy_1, \
                                         tg_zzzzzz_xxyyy_0, tg_zzzzzz_xxyyy_1, tg_zzzzzz_xxyyz_0, tg_zzzzzz_xxyyz_1, \
                                         tg_zzzzzz_xxyz_1, tg_zzzzzz_xxyzz_0, tg_zzzzzz_xxyzz_1, tg_zzzzzz_xxzz_1, \
                                         tg_zzzzzz_xxzzz_0, tg_zzzzzz_xxzzz_1, tg_zzzzzz_xyyy_1, tg_zzzzzz_xyyyy_0, \
                                         tg_zzzzzz_xyyyy_1, tg_zzzzzz_xyyyz_0, tg_zzzzzz_xyyyz_1, tg_zzzzzz_xyyz_1, \
                                         tg_zzzzzz_xyyzz_0, tg_zzzzzz_xyyzz_1, tg_zzzzzz_xyzz_1, tg_zzzzzz_xyzzz_0, \
                                         tg_zzzzzz_xyzzz_1, tg_zzzzzz_xzzz_1, tg_zzzzzz_xzzzz_0, tg_zzzzzz_xzzzz_1, \
                                         tg_zzzzzz_yyyy_1, tg_zzzzzz_yyyyy_0, tg_zzzzzz_yyyyy_1, tg_zzzzzz_yyyyz_0, \
                                         tg_zzzzzz_yyyyz_1, tg_zzzzzz_yyyz_1, tg_zzzzzz_yyyzz_0, tg_zzzzzz_yyyzz_1, \
                                         tg_zzzzzz_yyzz_1, tg_zzzzzz_yyzzz_0, tg_zzzzzz_yyzzz_1, tg_zzzzzz_yzzz_1, \
                                         tg_zzzzzz_yzzzz_0, tg_zzzzzz_yzzzz_1, tg_zzzzzz_zzzz_1, tg_zzzzzz_zzzzz_0, \
                                         tg_zzzzzz_zzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xzzzzzz_xxxxy_0[j] = pb_x * tg_zzzzzz_xxxxy_0[j] + fr * tg_zzzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xxxy_1[j];

                    tg_xzzzzzz_xxxxz_0[j] = pb_x * tg_zzzzzz_xxxxz_0[j] + fr * tg_zzzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xxxz_1[j];

                    tg_xzzzzzz_xxxyy_0[j] = pb_x * tg_zzzzzz_xxxyy_0[j] + fr * tg_zzzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxyy_1[j];

                    tg_xzzzzzz_xxxyz_0[j] = pb_x * tg_zzzzzz_xxxyz_0[j] + fr * tg_zzzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxyz_1[j];

                    tg_xzzzzzz_xxxzz_0[j] = pb_x * tg_zzzzzz_xxxzz_0[j] + fr * tg_zzzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxzz_1[j];

                    tg_xzzzzzz_xxyyy_0[j] = pb_x * tg_zzzzzz_xxyyy_0[j] + fr * tg_zzzzzz_xxyyy_1[j] + fl1_fxn * tg_zzzzzz_xyyy_1[j];

                    tg_xzzzzzz_xxyyz_0[j] = pb_x * tg_zzzzzz_xxyyz_0[j] + fr * tg_zzzzzz_xxyyz_1[j] + fl1_fxn * tg_zzzzzz_xyyz_1[j];

                    tg_xzzzzzz_xxyzz_0[j] = pb_x * tg_zzzzzz_xxyzz_0[j] + fr * tg_zzzzzz_xxyzz_1[j] + fl1_fxn * tg_zzzzzz_xyzz_1[j];

                    tg_xzzzzzz_xxzzz_0[j] = pb_x * tg_zzzzzz_xxzzz_0[j] + fr * tg_zzzzzz_xxzzz_1[j] + fl1_fxn * tg_zzzzzz_xzzz_1[j];

                    tg_xzzzzzz_xyyyy_0[j] = pb_x * tg_zzzzzz_xyyyy_0[j] + fr * tg_zzzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyyy_1[j];

                    tg_xzzzzzz_xyyyz_0[j] = pb_x * tg_zzzzzz_xyyyz_0[j] + fr * tg_zzzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyyz_1[j];

                    tg_xzzzzzz_xyyzz_0[j] = pb_x * tg_zzzzzz_xyyzz_0[j] + fr * tg_zzzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyzz_1[j];

                    tg_xzzzzzz_xyzzz_0[j] = pb_x * tg_zzzzzz_xyzzz_0[j] + fr * tg_zzzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yzzz_1[j];

                    tg_xzzzzzz_xzzzz_0[j] = pb_x * tg_zzzzzz_xzzzz_0[j] + fr * tg_zzzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zzzz_1[j];

                    tg_xzzzzzz_yyyyy_0[j] = pb_x * tg_zzzzzz_yyyyy_0[j] + fr * tg_zzzzzz_yyyyy_1[j];

                    tg_xzzzzzz_yyyyz_0[j] = pb_x * tg_zzzzzz_yyyyz_0[j] + fr * tg_zzzzzz_yyyyz_1[j];

                    tg_xzzzzzz_yyyzz_0[j] = pb_x * tg_zzzzzz_yyyzz_0[j] + fr * tg_zzzzzz_yyyzz_1[j];

                    tg_xzzzzzz_yyzzz_0[j] = pb_x * tg_zzzzzz_yyzzz_0[j] + fr * tg_zzzzzz_yyzzz_1[j];

                    tg_xzzzzzz_yzzzz_0[j] = pb_x * tg_zzzzzz_yzzzz_0[j] + fr * tg_zzzzzz_yzzzz_1[j];

                    tg_xzzzzzz_zzzzz_0[j] = pb_x * tg_zzzzzz_zzzzz_0[j] + fr * tg_zzzzzz_zzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyyyy_xxxxx_0[j] = pb_y * tg_yyyyyy_xxxxx_0[j] + fr * tg_yyyyyy_xxxxx_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxx_0[j] - tg_yyyyy_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyyy_xxxxy_0[j] = pb_y * tg_yyyyyy_xxxxy_0[j] + fr * tg_yyyyyy_xxxxy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxy_0[j] - tg_yyyyy_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxxx_1[j];

                    tg_yyyyyyy_xxxxz_0[j] = pb_y * tg_yyyyyy_xxxxz_0[j] + fr * tg_yyyyyy_xxxxz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxz_0[j] - tg_yyyyy_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyyy_xxxyy_0[j] = pb_y * tg_yyyyyy_xxxyy_0[j] + fr * tg_yyyyyy_xxxyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxyy_0[j] - tg_yyyyy_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xxxy_1[j];

                    tg_yyyyyyy_xxxyz_0[j] = pb_y * tg_yyyyyy_xxxyz_0[j] + fr * tg_yyyyyy_xxxyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxyz_0[j] - tg_yyyyy_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxxz_1[j];

                    tg_yyyyyyy_xxxzz_0[j] = pb_y * tg_yyyyyy_xxxzz_0[j] + fr * tg_yyyyyy_xxxzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxzz_0[j] - tg_yyyyy_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyyy_xxyyy_0[j] = pb_y * tg_yyyyyy_xxyyy_0[j] + fr * tg_yyyyyy_xxyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyyy_0[j] - tg_yyyyy_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_xxyy_1[j];

                    tg_yyyyyyy_xxyyz_0[j] = pb_y * tg_yyyyyy_xxyyz_0[j] + fr * tg_yyyyyy_xxyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyyz_0[j] - tg_yyyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xxyz_1[j];

                    tg_yyyyyyy_xxyzz_0[j] = pb_y * tg_yyyyyy_xxyzz_0[j] + fr * tg_yyyyyy_xxyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyzz_0[j] - tg_yyyyy_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxzz_1[j];

                    tg_yyyyyyy_xxzzz_0[j] = pb_y * tg_yyyyyy_xxzzz_0[j] + fr * tg_yyyyyy_xxzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxzzz_0[j] - tg_yyyyy_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyyy_xyyyy_0[j] = pb_y * tg_yyyyyy_xyyyy_0[j] + fr * tg_yyyyyy_xyyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyyy_0[j] - tg_yyyyy_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyy_xyyy_1[j];

                    tg_yyyyyyy_xyyyz_0[j] = pb_y * tg_yyyyyy_xyyyz_0[j] + fr * tg_yyyyyy_xyyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyyz_0[j] - tg_yyyyy_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_xyyz_1[j];

                    tg_yyyyyyy_xyyzz_0[j] = pb_y * tg_yyyyyy_xyyzz_0[j] + fr * tg_yyyyyy_xyyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyzz_0[j] - tg_yyyyy_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xyzz_1[j];

                    tg_yyyyyyy_xyzzz_0[j] = pb_y * tg_yyyyyy_xyzzz_0[j] + fr * tg_yyyyyy_xyzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyzzz_0[j] - tg_yyyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xzzz_1[j];

                    tg_yyyyyyy_xzzzz_0[j] = pb_y * tg_yyyyyy_xzzzz_0[j] + fr * tg_yyyyyy_xzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xzzzz_0[j] - tg_yyyyy_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyyy_yyyyy_0[j] = pb_y * tg_yyyyyy_yyyyy_0[j] + fr * tg_yyyyyy_yyyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyyy_0[j] - tg_yyyyy_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyy_yyyy_1[j];

                    tg_yyyyyyy_yyyyz_0[j] = pb_y * tg_yyyyyy_yyyyz_0[j] + fr * tg_yyyyyy_yyyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyyz_0[j] - tg_yyyyy_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyy_yyyz_1[j];

                    tg_yyyyyyy_yyyzz_0[j] = pb_y * tg_yyyyyy_yyyzz_0[j] + fr * tg_yyyyyy_yyyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyzz_0[j] - tg_yyyyy_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_yyzz_1[j];

                    tg_yyyyyyy_yyzzz_0[j] = pb_y * tg_yyyyyy_yyzzz_0[j] + fr * tg_yyyyyy_yyzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyzzz_0[j] - tg_yyyyy_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_yzzz_1[j];

                    tg_yyyyyyy_yzzzz_0[j] = pb_y * tg_yyyyyy_yzzzz_0[j] + fr * tg_yyyyyy_yzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yzzzz_0[j] - tg_yyyyy_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_zzzz_1[j];

                    tg_yyyyyyy_zzzzz_0[j] = pb_y * tg_yyyyyy_zzzzz_0[j] + fr * tg_yyyyyy_zzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_zzzzz_0[j] - tg_yyyyy_zzzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxxx_0[j] = pb_y * tg_yyyyyz_xxxxx_0[j] + fr * tg_yyyyyz_xxxxx_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxx_0[j] - tg_yyyyz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxxy_0[j] = pb_y * tg_yyyyyz_xxxxy_0[j] + fr * tg_yyyyyz_xxxxy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxy_0[j] - tg_yyyyz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxxx_1[j];

                    tg_yyyyyyz_xxxxz_0[j] = pb_y * tg_yyyyyz_xxxxz_0[j] + fr * tg_yyyyyz_xxxxz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxz_0[j] - tg_yyyyz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxyy_0[j] = pb_y * tg_yyyyyz_xxxyy_0[j] + fr * tg_yyyyyz_xxxyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxyy_0[j] - tg_yyyyz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xxxy_1[j];

                    tg_yyyyyyz_xxxyz_0[j] = pb_y * tg_yyyyyz_xxxyz_0[j] + fr * tg_yyyyyz_xxxyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxyz_0[j] - tg_yyyyz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxxz_1[j];

                    tg_yyyyyyz_xxxzz_0[j] = pb_y * tg_yyyyyz_xxxzz_0[j] + fr * tg_yyyyyz_xxxzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxzz_0[j] - tg_yyyyz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxyyy_0[j] = pb_y * tg_yyyyyz_xxyyy_0[j] + fr * tg_yyyyyz_xxyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyyy_0[j] - tg_yyyyz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_xxyy_1[j];

                    tg_yyyyyyz_xxyyz_0[j] = pb_y * tg_yyyyyz_xxyyz_0[j] + fr * tg_yyyyyz_xxyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyyz_0[j] - tg_yyyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xxyz_1[j];

                    tg_yyyyyyz_xxyzz_0[j] = pb_y * tg_yyyyyz_xxyzz_0[j] + fr * tg_yyyyyz_xxyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyzz_0[j] - tg_yyyyz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxzz_1[j];

                    tg_yyyyyyz_xxzzz_0[j] = pb_y * tg_yyyyyz_xxzzz_0[j] + fr * tg_yyyyyz_xxzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxzzz_0[j] - tg_yyyyz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xyyyy_0[j] = pb_y * tg_yyyyyz_xyyyy_0[j] + fr * tg_yyyyyz_xyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyyy_0[j] - tg_yyyyz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyz_xyyy_1[j];

                    tg_yyyyyyz_xyyyz_0[j] = pb_y * tg_yyyyyz_xyyyz_0[j] + fr * tg_yyyyyz_xyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyyz_0[j] - tg_yyyyz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_xyyz_1[j];

                    tg_yyyyyyz_xyyzz_0[j] = pb_y * tg_yyyyyz_xyyzz_0[j] + fr * tg_yyyyyz_xyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyzz_0[j] - tg_yyyyz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xyzz_1[j];

                    tg_yyyyyyz_xyzzz_0[j] = pb_y * tg_yyyyyz_xyzzz_0[j] + fr * tg_yyyyyz_xyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyzzz_0[j] - tg_yyyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xzzz_1[j];

                    tg_yyyyyyz_xzzzz_0[j] = pb_y * tg_yyyyyz_xzzzz_0[j] + fr * tg_yyyyyz_xzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xzzzz_0[j] - tg_yyyyz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_yyyyy_0[j] = pb_y * tg_yyyyyz_yyyyy_0[j] + fr * tg_yyyyyz_yyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyyy_0[j] - tg_yyyyz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyz_yyyy_1[j];

                    tg_yyyyyyz_yyyyz_0[j] = pb_y * tg_yyyyyz_yyyyz_0[j] + fr * tg_yyyyyz_yyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyyz_0[j] - tg_yyyyz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyz_yyyz_1[j];

                    tg_yyyyyyz_yyyzz_0[j] = pb_y * tg_yyyyyz_yyyzz_0[j] + fr * tg_yyyyyz_yyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyzz_0[j] - tg_yyyyz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_yyzz_1[j];

                    tg_yyyyyyz_yyzzz_0[j] = pb_y * tg_yyyyyz_yyzzz_0[j] + fr * tg_yyyyyz_yyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyzzz_0[j] - tg_yyyyz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_yzzz_1[j];

                    tg_yyyyyyz_yzzzz_0[j] = pb_y * tg_yyyyyz_yzzzz_0[j] + fr * tg_yyyyyz_yzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yzzzz_0[j] - tg_yyyyz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_zzzz_1[j];

                    tg_yyyyyyz_zzzzz_0[j] = pb_y * tg_yyyyyz_zzzzz_0[j] + fr * tg_yyyyyz_zzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_zzzzz_0[j] - tg_yyyyz_zzzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxxx_0[j] = pb_y * tg_yyyyzz_xxxxx_0[j] + fr * tg_yyyyzz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxx_0[j] - tg_yyyzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxxy_0[j] = pb_y * tg_yyyyzz_xxxxy_0[j] + fr * tg_yyyyzz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxy_0[j] - tg_yyyzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxxx_1[j];

                    tg_yyyyyzz_xxxxz_0[j] = pb_y * tg_yyyyzz_xxxxz_0[j] + fr * tg_yyyyzz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxz_0[j] - tg_yyyzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxyy_0[j] = pb_y * tg_yyyyzz_xxxyy_0[j] + fr * tg_yyyyzz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxyy_0[j] - tg_yyyzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xxxy_1[j];

                    tg_yyyyyzz_xxxyz_0[j] = pb_y * tg_yyyyzz_xxxyz_0[j] + fr * tg_yyyyzz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxyz_0[j] - tg_yyyzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxxz_1[j];

                    tg_yyyyyzz_xxxzz_0[j] = pb_y * tg_yyyyzz_xxxzz_0[j] + fr * tg_yyyyzz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxzz_0[j] - tg_yyyzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxyyy_0[j] = pb_y * tg_yyyyzz_xxyyy_0[j] + fr * tg_yyyyzz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyyy_0[j] - tg_yyyzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_xxyy_1[j];

                    tg_yyyyyzz_xxyyz_0[j] = pb_y * tg_yyyyzz_xxyyz_0[j] + fr * tg_yyyyzz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyyz_0[j] - tg_yyyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xxyz_1[j];

                    tg_yyyyyzz_xxyzz_0[j] = pb_y * tg_yyyyzz_xxyzz_0[j] + fr * tg_yyyyzz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyzz_0[j] - tg_yyyzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxzz_1[j];

                    tg_yyyyyzz_xxzzz_0[j] = pb_y * tg_yyyyzz_xxzzz_0[j] + fr * tg_yyyyzz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxzzz_0[j] - tg_yyyzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xyyyy_0[j] = pb_y * tg_yyyyzz_xyyyy_0[j] + fr * tg_yyyyzz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyyy_0[j] - tg_yyyzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzz_xyyy_1[j];

                    tg_yyyyyzz_xyyyz_0[j] = pb_y * tg_yyyyzz_xyyyz_0[j] + fr * tg_yyyyzz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyyz_0[j] - tg_yyyzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_xyyz_1[j];

                    tg_yyyyyzz_xyyzz_0[j] = pb_y * tg_yyyyzz_xyyzz_0[j] + fr * tg_yyyyzz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyzz_0[j] - tg_yyyzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xyzz_1[j];

                    tg_yyyyyzz_xyzzz_0[j] = pb_y * tg_yyyyzz_xyzzz_0[j] + fr * tg_yyyyzz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyzzz_0[j] - tg_yyyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xzzz_1[j];

                    tg_yyyyyzz_xzzzz_0[j] = pb_y * tg_yyyyzz_xzzzz_0[j] + fr * tg_yyyyzz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xzzzz_0[j] - tg_yyyzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_yyyyy_0[j] = pb_y * tg_yyyyzz_yyyyy_0[j] + fr * tg_yyyyzz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyyy_0[j] - tg_yyyzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyzz_yyyy_1[j];

                    tg_yyyyyzz_yyyyz_0[j] = pb_y * tg_yyyyzz_yyyyz_0[j] + fr * tg_yyyyzz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyyz_0[j] - tg_yyyzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzz_yyyz_1[j];

                    tg_yyyyyzz_yyyzz_0[j] = pb_y * tg_yyyyzz_yyyzz_0[j] + fr * tg_yyyyzz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyzz_0[j] - tg_yyyzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_yyzz_1[j];

                    tg_yyyyyzz_yyzzz_0[j] = pb_y * tg_yyyyzz_yyzzz_0[j] + fr * tg_yyyyzz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyzzz_0[j] - tg_yyyzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_yzzz_1[j];

                    tg_yyyyyzz_yzzzz_0[j] = pb_y * tg_yyyyzz_yzzzz_0[j] + fr * tg_yyyyzz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yzzzz_0[j] - tg_yyyzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_zzzz_1[j];

                    tg_yyyyyzz_zzzzz_0[j] = pb_y * tg_yyyyzz_zzzzz_0[j] + fr * tg_yyyyzz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_zzzzz_0[j] - tg_yyyzz_zzzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxxx_0[j] = pb_y * tg_yyyzzz_xxxxx_0[j] + fr * tg_yyyzzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxx_0[j] - tg_yyzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxxy_0[j] = pb_y * tg_yyyzzz_xxxxy_0[j] + fr * tg_yyyzzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxy_0[j] - tg_yyzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxxx_1[j];

                    tg_yyyyzzz_xxxxz_0[j] = pb_y * tg_yyyzzz_xxxxz_0[j] + fr * tg_yyyzzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxz_0[j] - tg_yyzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxyy_0[j] = pb_y * tg_yyyzzz_xxxyy_0[j] + fr * tg_yyyzzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxyy_0[j] - tg_yyzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xxxy_1[j];

                    tg_yyyyzzz_xxxyz_0[j] = pb_y * tg_yyyzzz_xxxyz_0[j] + fr * tg_yyyzzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxyz_0[j] - tg_yyzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxxz_1[j];

                    tg_yyyyzzz_xxxzz_0[j] = pb_y * tg_yyyzzz_xxxzz_0[j] + fr * tg_yyyzzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxzz_0[j] - tg_yyzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxyyy_0[j] = pb_y * tg_yyyzzz_xxyyy_0[j] + fr * tg_yyyzzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyyy_0[j] - tg_yyzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_xxyy_1[j];

                    tg_yyyyzzz_xxyyz_0[j] = pb_y * tg_yyyzzz_xxyyz_0[j] + fr * tg_yyyzzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyyz_0[j] - tg_yyzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xxyz_1[j];

                    tg_yyyyzzz_xxyzz_0[j] = pb_y * tg_yyyzzz_xxyzz_0[j] + fr * tg_yyyzzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyzz_0[j] - tg_yyzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxzz_1[j];

                    tg_yyyyzzz_xxzzz_0[j] = pb_y * tg_yyyzzz_xxzzz_0[j] + fr * tg_yyyzzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxzzz_0[j] - tg_yyzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xyyyy_0[j] = pb_y * tg_yyyzzz_xyyyy_0[j] + fr * tg_yyyzzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyyy_0[j] - tg_yyzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzz_xyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSH_662_756(      CMemBlock2D<double>* primBuffer,
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
                                             {7, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
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

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_yyyzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 515); 

                auto tg_yyyzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 516); 

                auto tg_yyyzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 517); 

                auto tg_yyyzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 518); 

                auto tg_yyyzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 519); 

                auto tg_yyyzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 520); 

                auto tg_yyyzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 521); 

                auto tg_yyyzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 522); 

                auto tg_yyyzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 523); 

                auto tg_yyyzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 524); 

                auto tg_yyzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 525); 

                auto tg_yyzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 526); 

                auto tg_yyzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 527); 

                auto tg_yyzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 528); 

                auto tg_yyzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 529); 

                auto tg_yyzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 530); 

                auto tg_yyzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 531); 

                auto tg_yyzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 532); 

                auto tg_yyzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 533); 

                auto tg_yyzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 534); 

                auto tg_yyzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 535); 

                auto tg_yyzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 536); 

                auto tg_yyzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 537); 

                auto tg_yyzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 538); 

                auto tg_yyzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 539); 

                auto tg_yyzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 540); 

                auto tg_yyzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 541); 

                auto tg_yyzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 542); 

                auto tg_yyzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 543); 

                auto tg_yyzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 544); 

                auto tg_yyzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 545); 

                auto tg_yzzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 546); 

                auto tg_yzzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 547); 

                auto tg_yzzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 548); 

                auto tg_yzzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 549); 

                auto tg_yzzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 550); 

                auto tg_yzzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 551); 

                auto tg_yzzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 552); 

                auto tg_yzzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 553); 

                auto tg_yzzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 554); 

                auto tg_yzzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 555); 

                auto tg_yzzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 556); 

                auto tg_yzzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 557); 

                auto tg_yzzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 558); 

                auto tg_yzzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 559); 

                auto tg_yzzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 560); 

                auto tg_yzzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 561); 

                auto tg_yzzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 562); 

                auto tg_yzzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 563); 

                auto tg_yzzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 564); 

                auto tg_yzzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 565); 

                auto tg_yzzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 566); 

                auto tg_zzzzzz_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 567); 

                auto tg_zzzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 568); 

                auto tg_zzzzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 569); 

                auto tg_zzzzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 570); 

                auto tg_zzzzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 571); 

                auto tg_zzzzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 572); 

                auto tg_zzzzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 573); 

                auto tg_zzzzzz_xxyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 574); 

                auto tg_zzzzzz_xxyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 575); 

                auto tg_zzzzzz_xxzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 576); 

                auto tg_zzzzzz_xyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 577); 

                auto tg_zzzzzz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 578); 

                auto tg_zzzzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 579); 

                auto tg_zzzzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 580); 

                auto tg_zzzzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 581); 

                auto tg_zzzzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 582); 

                auto tg_zzzzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 583); 

                auto tg_zzzzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 584); 

                auto tg_zzzzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 585); 

                auto tg_zzzzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 586); 

                auto tg_zzzzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 587); 

                auto tg_yyyzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 515); 

                auto tg_yyyzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 516); 

                auto tg_yyyzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 517); 

                auto tg_yyyzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 518); 

                auto tg_yyyzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 519); 

                auto tg_yyyzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 520); 

                auto tg_yyyzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 521); 

                auto tg_yyyzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 522); 

                auto tg_yyyzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 523); 

                auto tg_yyyzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 524); 

                auto tg_yyzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 525); 

                auto tg_yyzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 526); 

                auto tg_yyzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 527); 

                auto tg_yyzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 528); 

                auto tg_yyzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 529); 

                auto tg_yyzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 530); 

                auto tg_yyzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 531); 

                auto tg_yyzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 532); 

                auto tg_yyzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 533); 

                auto tg_yyzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 534); 

                auto tg_yyzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 535); 

                auto tg_yyzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 536); 

                auto tg_yyzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 537); 

                auto tg_yyzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 538); 

                auto tg_yyzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 539); 

                auto tg_yyzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 540); 

                auto tg_yyzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 541); 

                auto tg_yyzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 542); 

                auto tg_yyzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 543); 

                auto tg_yyzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 544); 

                auto tg_yyzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 545); 

                auto tg_yzzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 546); 

                auto tg_yzzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 547); 

                auto tg_yzzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 548); 

                auto tg_yzzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 549); 

                auto tg_yzzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 550); 

                auto tg_yzzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 551); 

                auto tg_yzzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 552); 

                auto tg_yzzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 553); 

                auto tg_yzzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 554); 

                auto tg_yzzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 555); 

                auto tg_yzzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 556); 

                auto tg_yzzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 557); 

                auto tg_yzzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 558); 

                auto tg_yzzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 559); 

                auto tg_yzzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 560); 

                auto tg_yzzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 561); 

                auto tg_yzzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 562); 

                auto tg_yzzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 563); 

                auto tg_yzzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 564); 

                auto tg_yzzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 565); 

                auto tg_yzzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 566); 

                auto tg_zzzzzz_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 567); 

                auto tg_zzzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 568); 

                auto tg_zzzzzz_xxxxz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 569); 

                auto tg_zzzzzz_xxxyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 570); 

                auto tg_zzzzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 571); 

                auto tg_zzzzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 572); 

                auto tg_zzzzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 573); 

                auto tg_zzzzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 574); 

                auto tg_zzzzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 575); 

                auto tg_zzzzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 576); 

                auto tg_zzzzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 577); 

                auto tg_zzzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 578); 

                auto tg_zzzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 579); 

                auto tg_zzzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 580); 

                auto tg_zzzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 581); 

                auto tg_zzzzzz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 582); 

                auto tg_zzzzzz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 583); 

                auto tg_zzzzzz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 584); 

                auto tg_zzzzzz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 585); 

                auto tg_zzzzzz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 586); 

                auto tg_zzzzzz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 587); 

                auto tg_yyzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 389); 

                auto tg_yyzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 390); 

                auto tg_yyzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 391); 

                auto tg_yyzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 392); 

                auto tg_yyzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 393); 

                auto tg_yyzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 394); 

                auto tg_yyzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 395); 

                auto tg_yyzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 396); 

                auto tg_yyzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 397); 

                auto tg_yyzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 398); 

                auto tg_yzzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 399); 

                auto tg_yzzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 400); 

                auto tg_yzzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 401); 

                auto tg_yzzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 402); 

                auto tg_yzzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 403); 

                auto tg_yzzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 404); 

                auto tg_yzzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 405); 

                auto tg_yzzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 406); 

                auto tg_yzzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 407); 

                auto tg_yzzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 408); 

                auto tg_yzzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 409); 

                auto tg_yzzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 410); 

                auto tg_yzzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 411); 

                auto tg_yzzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 412); 

                auto tg_yzzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 413); 

                auto tg_yzzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 414); 

                auto tg_yzzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 415); 

                auto tg_yzzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 416); 

                auto tg_yzzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 417); 

                auto tg_yzzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 418); 

                auto tg_yzzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 419); 

                auto tg_zzzzz_xxxxx_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 420); 

                auto tg_zzzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 421); 

                auto tg_zzzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 422); 

                auto tg_zzzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 423); 

                auto tg_zzzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 424); 

                auto tg_zzzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 425); 

                auto tg_zzzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 426); 

                auto tg_zzzzz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 427); 

                auto tg_zzzzz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 428); 

                auto tg_zzzzz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 429); 

                auto tg_zzzzz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 430); 

                auto tg_zzzzz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 431); 

                auto tg_zzzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 432); 

                auto tg_zzzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 433); 

                auto tg_zzzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 434); 

                auto tg_zzzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 435); 

                auto tg_zzzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 436); 

                auto tg_zzzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 437); 

                auto tg_zzzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 438); 

                auto tg_zzzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 439); 

                auto tg_zzzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 440); 

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

                auto tg_yyyzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 367); 

                auto tg_yyyzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 368); 

                auto tg_yyyzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 369); 

                auto tg_yyyzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 370); 

                auto tg_yyyzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 371); 

                auto tg_yyyzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 372); 

                auto tg_yyyzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 373); 

                auto tg_yyyzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 374); 

                auto tg_yyzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 375); 

                auto tg_yyzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 376); 

                auto tg_yyzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 377); 

                auto tg_yyzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 378); 

                auto tg_yyzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 379); 

                auto tg_yyzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 380); 

                auto tg_yyzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 381); 

                auto tg_yyzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 382); 

                auto tg_yyzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 383); 

                auto tg_yyzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 384); 

                auto tg_yyzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 385); 

                auto tg_yyzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 386); 

                auto tg_yyzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 387); 

                auto tg_yyzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 388); 

                auto tg_yyzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 389); 

                auto tg_yzzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 390); 

                auto tg_yzzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 391); 

                auto tg_yzzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 392); 

                auto tg_yzzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 393); 

                auto tg_yzzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 394); 

                auto tg_yzzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 395); 

                auto tg_yzzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 396); 

                auto tg_yzzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 397); 

                auto tg_yzzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 398); 

                auto tg_yzzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 399); 

                auto tg_yzzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 400); 

                auto tg_yzzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 401); 

                auto tg_yzzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 402); 

                auto tg_yzzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 403); 

                auto tg_yzzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 404); 

                auto tg_zzzzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 405); 

                auto tg_zzzzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 406); 

                auto tg_zzzzzz_xxxz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 407); 

                auto tg_zzzzzz_xxyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 408); 

                auto tg_zzzzzz_xxyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 409); 

                auto tg_zzzzzz_xxzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 410); 

                auto tg_zzzzzz_xyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 411); 

                auto tg_zzzzzz_xyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 412); 

                auto tg_zzzzzz_xyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 413); 

                auto tg_zzzzzz_xzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 414); 

                auto tg_zzzzzz_yyyy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 415); 

                auto tg_zzzzzz_yyyz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 416); 

                auto tg_zzzzzz_yyzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 417); 

                auto tg_zzzzzz_yzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 418); 

                auto tg_zzzzzz_zzzz_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 419); 

                // set up pointers to integrals

                auto tg_yyyyzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 662); 

                auto tg_yyyyzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 663); 

                auto tg_yyyyzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 664); 

                auto tg_yyyyzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 665); 

                auto tg_yyyyzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 666); 

                auto tg_yyyyzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 667); 

                auto tg_yyyyzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 668); 

                auto tg_yyyyzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 669); 

                auto tg_yyyyzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 670); 

                auto tg_yyyyzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 671); 

                auto tg_yyyzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 672); 

                auto tg_yyyzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 673); 

                auto tg_yyyzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 674); 

                auto tg_yyyzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 675); 

                auto tg_yyyzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 676); 

                auto tg_yyyzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 677); 

                auto tg_yyyzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 678); 

                auto tg_yyyzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 679); 

                auto tg_yyyzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 680); 

                auto tg_yyyzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 681); 

                auto tg_yyyzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 682); 

                auto tg_yyyzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 683); 

                auto tg_yyyzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 684); 

                auto tg_yyyzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 685); 

                auto tg_yyyzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 686); 

                auto tg_yyyzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 687); 

                auto tg_yyyzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 688); 

                auto tg_yyyzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 689); 

                auto tg_yyyzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 690); 

                auto tg_yyyzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 691); 

                auto tg_yyyzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 692); 

                auto tg_yyzzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 693); 

                auto tg_yyzzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 694); 

                auto tg_yyzzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 695); 

                auto tg_yyzzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 696); 

                auto tg_yyzzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 697); 

                auto tg_yyzzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 698); 

                auto tg_yyzzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 699); 

                auto tg_yyzzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 700); 

                auto tg_yyzzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 701); 

                auto tg_yyzzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 702); 

                auto tg_yyzzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 703); 

                auto tg_yyzzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 704); 

                auto tg_yyzzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 705); 

                auto tg_yyzzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 706); 

                auto tg_yyzzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 707); 

                auto tg_yyzzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 708); 

                auto tg_yyzzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 709); 

                auto tg_yyzzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 710); 

                auto tg_yyzzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 711); 

                auto tg_yyzzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 712); 

                auto tg_yyzzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 713); 

                auto tg_yzzzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 714); 

                auto tg_yzzzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 715); 

                auto tg_yzzzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 716); 

                auto tg_yzzzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 717); 

                auto tg_yzzzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 718); 

                auto tg_yzzzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 719); 

                auto tg_yzzzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 720); 

                auto tg_yzzzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 721); 

                auto tg_yzzzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 722); 

                auto tg_yzzzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 723); 

                auto tg_yzzzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 724); 

                auto tg_yzzzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 725); 

                auto tg_yzzzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 726); 

                auto tg_yzzzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 727); 

                auto tg_yzzzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 728); 

                auto tg_yzzzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 729); 

                auto tg_yzzzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 730); 

                auto tg_yzzzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 731); 

                auto tg_yzzzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 732); 

                auto tg_yzzzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 733); 

                auto tg_yzzzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 734); 

                auto tg_zzzzzzz_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 735); 

                auto tg_zzzzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 736); 

                auto tg_zzzzzzz_xxxxz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 737); 

                auto tg_zzzzzzz_xxxyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 738); 

                auto tg_zzzzzzz_xxxyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 739); 

                auto tg_zzzzzzz_xxxzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 740); 

                auto tg_zzzzzzz_xxyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 741); 

                auto tg_zzzzzzz_xxyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 742); 

                auto tg_zzzzzzz_xxyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 743); 

                auto tg_zzzzzzz_xxzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 744); 

                auto tg_zzzzzzz_xyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 745); 

                auto tg_zzzzzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 746); 

                auto tg_zzzzzzz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 747); 

                auto tg_zzzzzzz_xyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 748); 

                auto tg_zzzzzzz_xzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 749); 

                auto tg_zzzzzzz_yyyyy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 750); 

                auto tg_zzzzzzz_yyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 751); 

                auto tg_zzzzzzz_yyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 752); 

                auto tg_zzzzzzz_yyzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 753); 

                auto tg_zzzzzzz_yzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 754); 

                auto tg_zzzzzzz_zzzzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 755); 

                // Batch of Integrals (662,756)

                #pragma omp simd aligned(fxn, fza, tg_yyyyzzz_xyyyz_0, tg_yyyyzzz_xyyzz_0, tg_yyyyzzz_xyzzz_0, \
                                         tg_yyyyzzz_xzzzz_0, tg_yyyyzzz_yyyyy_0, tg_yyyyzzz_yyyyz_0, tg_yyyyzzz_yyyzz_0, \
                                         tg_yyyyzzz_yyzzz_0, tg_yyyyzzz_yzzzz_0, tg_yyyyzzz_zzzzz_0, tg_yyyzzz_xyyyz_0, \
                                         tg_yyyzzz_xyyyz_1, tg_yyyzzz_xyyz_1, tg_yyyzzz_xyyzz_0, tg_yyyzzz_xyyzz_1, \
                                         tg_yyyzzz_xyzz_1, tg_yyyzzz_xyzzz_0, tg_yyyzzz_xyzzz_1, tg_yyyzzz_xzzz_1, \
                                         tg_yyyzzz_xzzzz_0, tg_yyyzzz_xzzzz_1, tg_yyyzzz_yyyy_1, tg_yyyzzz_yyyyy_0, \
                                         tg_yyyzzz_yyyyy_1, tg_yyyzzz_yyyyz_0, tg_yyyzzz_yyyyz_1, tg_yyyzzz_yyyz_1, \
                                         tg_yyyzzz_yyyzz_0, tg_yyyzzz_yyyzz_1, tg_yyyzzz_yyzz_1, tg_yyyzzz_yyzzz_0, \
                                         tg_yyyzzz_yyzzz_1, tg_yyyzzz_yzzz_1, tg_yyyzzz_yzzzz_0, tg_yyyzzz_yzzzz_1, \
                                         tg_yyyzzz_zzzz_1, tg_yyyzzz_zzzzz_0, tg_yyyzzz_zzzzz_1, tg_yyyzzzz_xxxxx_0, \
                                         tg_yyyzzzz_xxxxy_0, tg_yyyzzzz_xxxxz_0, tg_yyyzzzz_xxxyy_0, tg_yyyzzzz_xxxyz_0, \
                                         tg_yyyzzzz_xxxzz_0, tg_yyyzzzz_xxyyy_0, tg_yyyzzzz_xxyyz_0, tg_yyyzzzz_xxyzz_0, \
                                         tg_yyyzzzz_xxzzz_0, tg_yyyzzzz_xyyyy_0, tg_yyyzzzz_xyyyz_0, tg_yyyzzzz_xyyzz_0, \
                                         tg_yyyzzzz_xyzzz_0, tg_yyyzzzz_xzzzz_0, tg_yyyzzzz_yyyyy_0, tg_yyyzzzz_yyyyz_0, \
                                         tg_yyyzzzz_yyyzz_0, tg_yyyzzzz_yyzzz_0, tg_yyyzzzz_yzzzz_0, tg_yyyzzzz_zzzzz_0, \
                                         tg_yyzzz_xyyyz_0, tg_yyzzz_xyyyz_1, tg_yyzzz_xyyzz_0, tg_yyzzz_xyyzz_1, \
                                         tg_yyzzz_xyzzz_0, tg_yyzzz_xyzzz_1, tg_yyzzz_xzzzz_0, tg_yyzzz_xzzzz_1, \
                                         tg_yyzzz_yyyyy_0, tg_yyzzz_yyyyy_1, tg_yyzzz_yyyyz_0, tg_yyzzz_yyyyz_1, \
                                         tg_yyzzz_yyyzz_0, tg_yyzzz_yyyzz_1, tg_yyzzz_yyzzz_0, tg_yyzzz_yyzzz_1, \
                                         tg_yyzzz_yzzzz_0, tg_yyzzz_yzzzz_1, tg_yyzzz_zzzzz_0, tg_yyzzz_zzzzz_1, \
                                         tg_yyzzzz_xxxx_1, tg_yyzzzz_xxxxx_0, tg_yyzzzz_xxxxx_1, tg_yyzzzz_xxxxy_0, \
                                         tg_yyzzzz_xxxxy_1, tg_yyzzzz_xxxxz_0, tg_yyzzzz_xxxxz_1, tg_yyzzzz_xxxy_1, \
                                         tg_yyzzzz_xxxyy_0, tg_yyzzzz_xxxyy_1, tg_yyzzzz_xxxyz_0, tg_yyzzzz_xxxyz_1, \
                                         tg_yyzzzz_xxxz_1, tg_yyzzzz_xxxzz_0, tg_yyzzzz_xxxzz_1, tg_yyzzzz_xxyy_1, \
                                         tg_yyzzzz_xxyyy_0, tg_yyzzzz_xxyyy_1, tg_yyzzzz_xxyyz_0, tg_yyzzzz_xxyyz_1, \
                                         tg_yyzzzz_xxyz_1, tg_yyzzzz_xxyzz_0, tg_yyzzzz_xxyzz_1, tg_yyzzzz_xxzz_1, \
                                         tg_yyzzzz_xxzzz_0, tg_yyzzzz_xxzzz_1, tg_yyzzzz_xyyy_1, tg_yyzzzz_xyyyy_0, \
                                         tg_yyzzzz_xyyyy_1, tg_yyzzzz_xyyyz_0, tg_yyzzzz_xyyyz_1, tg_yyzzzz_xyyz_1, \
                                         tg_yyzzzz_xyyzz_0, tg_yyzzzz_xyyzz_1, tg_yyzzzz_xyzz_1, tg_yyzzzz_xyzzz_0, \
                                         tg_yyzzzz_xyzzz_1, tg_yyzzzz_xzzz_1, tg_yyzzzz_xzzzz_0, tg_yyzzzz_xzzzz_1, \
                                         tg_yyzzzz_yyyy_1, tg_yyzzzz_yyyyy_0, tg_yyzzzz_yyyyy_1, tg_yyzzzz_yyyyz_0, \
                                         tg_yyzzzz_yyyyz_1, tg_yyzzzz_yyyz_1, tg_yyzzzz_yyyzz_0, tg_yyzzzz_yyyzz_1, \
                                         tg_yyzzzz_yyzz_1, tg_yyzzzz_yyzzz_0, tg_yyzzzz_yyzzz_1, tg_yyzzzz_yzzz_1, \
                                         tg_yyzzzz_yzzzz_0, tg_yyzzzz_yzzzz_1, tg_yyzzzz_zzzz_1, tg_yyzzzz_zzzzz_0, \
                                         tg_yyzzzz_zzzzz_1, tg_yyzzzzz_xxxxx_0, tg_yyzzzzz_xxxxy_0, tg_yyzzzzz_xxxxz_0, \
                                         tg_yyzzzzz_xxxyy_0, tg_yyzzzzz_xxxyz_0, tg_yyzzzzz_xxxzz_0, tg_yyzzzzz_xxyyy_0, \
                                         tg_yyzzzzz_xxyyz_0, tg_yyzzzzz_xxyzz_0, tg_yyzzzzz_xxzzz_0, tg_yyzzzzz_xyyyy_0, \
                                         tg_yyzzzzz_xyyyz_0, tg_yyzzzzz_xyyzz_0, tg_yyzzzzz_xyzzz_0, tg_yyzzzzz_xzzzz_0, \
                                         tg_yyzzzzz_yyyyy_0, tg_yyzzzzz_yyyyz_0, tg_yyzzzzz_yyyzz_0, tg_yyzzzzz_yyzzz_0, \
                                         tg_yyzzzzz_yzzzz_0, tg_yyzzzzz_zzzzz_0, tg_yzzzz_xxxxx_0, tg_yzzzz_xxxxx_1, \
                                         tg_yzzzz_xxxxy_0, tg_yzzzz_xxxxy_1, tg_yzzzz_xxxxz_0, tg_yzzzz_xxxxz_1, \
                                         tg_yzzzz_xxxyy_0, tg_yzzzz_xxxyy_1, tg_yzzzz_xxxyz_0, tg_yzzzz_xxxyz_1, \
                                         tg_yzzzz_xxxzz_0, tg_yzzzz_xxxzz_1, tg_yzzzz_xxyyy_0, tg_yzzzz_xxyyy_1, \
                                         tg_yzzzz_xxyyz_0, tg_yzzzz_xxyyz_1, tg_yzzzz_xxyzz_0, tg_yzzzz_xxyzz_1, \
                                         tg_yzzzz_xxzzz_0, tg_yzzzz_xxzzz_1, tg_yzzzz_xyyyy_0, tg_yzzzz_xyyyy_1, \
                                         tg_yzzzz_xyyyz_0, tg_yzzzz_xyyyz_1, tg_yzzzz_xyyzz_0, tg_yzzzz_xyyzz_1, \
                                         tg_yzzzz_xyzzz_0, tg_yzzzz_xyzzz_1, tg_yzzzz_xzzzz_0, tg_yzzzz_xzzzz_1, \
                                         tg_yzzzz_yyyyy_0, tg_yzzzz_yyyyy_1, tg_yzzzz_yyyyz_0, tg_yzzzz_yyyyz_1, \
                                         tg_yzzzz_yyyzz_0, tg_yzzzz_yyyzz_1, tg_yzzzz_yyzzz_0, tg_yzzzz_yyzzz_1, \
                                         tg_yzzzz_yzzzz_0, tg_yzzzz_yzzzz_1, tg_yzzzz_zzzzz_0, tg_yzzzz_zzzzz_1, \
                                         tg_yzzzzz_xxxx_1, tg_yzzzzz_xxxxx_0, tg_yzzzzz_xxxxx_1, tg_yzzzzz_xxxxy_0, \
                                         tg_yzzzzz_xxxxy_1, tg_yzzzzz_xxxxz_0, tg_yzzzzz_xxxxz_1, tg_yzzzzz_xxxy_1, \
                                         tg_yzzzzz_xxxyy_0, tg_yzzzzz_xxxyy_1, tg_yzzzzz_xxxyz_0, tg_yzzzzz_xxxyz_1, \
                                         tg_yzzzzz_xxxz_1, tg_yzzzzz_xxxzz_0, tg_yzzzzz_xxxzz_1, tg_yzzzzz_xxyy_1, \
                                         tg_yzzzzz_xxyyy_0, tg_yzzzzz_xxyyy_1, tg_yzzzzz_xxyyz_0, tg_yzzzzz_xxyyz_1, \
                                         tg_yzzzzz_xxyz_1, tg_yzzzzz_xxyzz_0, tg_yzzzzz_xxyzz_1, tg_yzzzzz_xxzz_1, \
                                         tg_yzzzzz_xxzzz_0, tg_yzzzzz_xxzzz_1, tg_yzzzzz_xyyy_1, tg_yzzzzz_xyyyy_0, \
                                         tg_yzzzzz_xyyyy_1, tg_yzzzzz_xyyyz_0, tg_yzzzzz_xyyyz_1, tg_yzzzzz_xyyz_1, \
                                         tg_yzzzzz_xyyzz_0, tg_yzzzzz_xyyzz_1, tg_yzzzzz_xyzz_1, tg_yzzzzz_xyzzz_0, \
                                         tg_yzzzzz_xyzzz_1, tg_yzzzzz_xzzz_1, tg_yzzzzz_xzzzz_0, tg_yzzzzz_xzzzz_1, \
                                         tg_yzzzzz_yyyy_1, tg_yzzzzz_yyyyy_0, tg_yzzzzz_yyyyy_1, tg_yzzzzz_yyyyz_0, \
                                         tg_yzzzzz_yyyyz_1, tg_yzzzzz_yyyz_1, tg_yzzzzz_yyyzz_0, tg_yzzzzz_yyyzz_1, \
                                         tg_yzzzzz_yyzz_1, tg_yzzzzz_yyzzz_0, tg_yzzzzz_yyzzz_1, tg_yzzzzz_yzzz_1, \
                                         tg_yzzzzz_yzzzz_0, tg_yzzzzz_yzzzz_1, tg_yzzzzz_zzzz_1, tg_yzzzzz_zzzzz_0, \
                                         tg_yzzzzz_zzzzz_1, tg_yzzzzzz_xxxxx_0, tg_yzzzzzz_xxxxy_0, tg_yzzzzzz_xxxxz_0, \
                                         tg_yzzzzzz_xxxyy_0, tg_yzzzzzz_xxxyz_0, tg_yzzzzzz_xxxzz_0, tg_yzzzzzz_xxyyy_0, \
                                         tg_yzzzzzz_xxyyz_0, tg_yzzzzzz_xxyzz_0, tg_yzzzzzz_xxzzz_0, tg_yzzzzzz_xyyyy_0, \
                                         tg_yzzzzzz_xyyyz_0, tg_yzzzzzz_xyyzz_0, tg_yzzzzzz_xyzzz_0, tg_yzzzzzz_xzzzz_0, \
                                         tg_yzzzzzz_yyyyy_0, tg_yzzzzzz_yyyyz_0, tg_yzzzzzz_yyyzz_0, tg_yzzzzzz_yyzzz_0, \
                                         tg_yzzzzzz_yzzzz_0, tg_yzzzzzz_zzzzz_0, tg_zzzzz_xxxxx_0, tg_zzzzz_xxxxx_1, \
                                         tg_zzzzz_xxxxy_0, tg_zzzzz_xxxxy_1, tg_zzzzz_xxxxz_0, tg_zzzzz_xxxxz_1, \
                                         tg_zzzzz_xxxyy_0, tg_zzzzz_xxxyy_1, tg_zzzzz_xxxyz_0, tg_zzzzz_xxxyz_1, \
                                         tg_zzzzz_xxxzz_0, tg_zzzzz_xxxzz_1, tg_zzzzz_xxyyy_0, tg_zzzzz_xxyyy_1, \
                                         tg_zzzzz_xxyyz_0, tg_zzzzz_xxyyz_1, tg_zzzzz_xxyzz_0, tg_zzzzz_xxyzz_1, \
                                         tg_zzzzz_xxzzz_0, tg_zzzzz_xxzzz_1, tg_zzzzz_xyyyy_0, tg_zzzzz_xyyyy_1, \
                                         tg_zzzzz_xyyyz_0, tg_zzzzz_xyyyz_1, tg_zzzzz_xyyzz_0, tg_zzzzz_xyyzz_1, \
                                         tg_zzzzz_xyzzz_0, tg_zzzzz_xyzzz_1, tg_zzzzz_xzzzz_0, tg_zzzzz_xzzzz_1, \
                                         tg_zzzzz_yyyyy_0, tg_zzzzz_yyyyy_1, tg_zzzzz_yyyyz_0, tg_zzzzz_yyyyz_1, \
                                         tg_zzzzz_yyyzz_0, tg_zzzzz_yyyzz_1, tg_zzzzz_yyzzz_0, tg_zzzzz_yyzzz_1, \
                                         tg_zzzzz_yzzzz_0, tg_zzzzz_yzzzz_1, tg_zzzzz_zzzzz_0, tg_zzzzz_zzzzz_1, \
                                         tg_zzzzzz_xxxx_1, tg_zzzzzz_xxxxx_0, tg_zzzzzz_xxxxx_1, tg_zzzzzz_xxxxy_0, \
                                         tg_zzzzzz_xxxxy_1, tg_zzzzzz_xxxxz_0, tg_zzzzzz_xxxxz_1, tg_zzzzzz_xxxy_1, \
                                         tg_zzzzzz_xxxyy_0, tg_zzzzzz_xxxyy_1, tg_zzzzzz_xxxyz_0, tg_zzzzzz_xxxyz_1, \
                                         tg_zzzzzz_xxxz_1, tg_zzzzzz_xxxzz_0, tg_zzzzzz_xxxzz_1, tg_zzzzzz_xxyy_1, \
                                         tg_zzzzzz_xxyyy_0, tg_zzzzzz_xxyyy_1, tg_zzzzzz_xxyyz_0, tg_zzzzzz_xxyyz_1, \
                                         tg_zzzzzz_xxyz_1, tg_zzzzzz_xxyzz_0, tg_zzzzzz_xxyzz_1, tg_zzzzzz_xxzz_1, \
                                         tg_zzzzzz_xxzzz_0, tg_zzzzzz_xxzzz_1, tg_zzzzzz_xyyy_1, tg_zzzzzz_xyyyy_0, \
                                         tg_zzzzzz_xyyyy_1, tg_zzzzzz_xyyyz_0, tg_zzzzzz_xyyyz_1, tg_zzzzzz_xyyz_1, \
                                         tg_zzzzzz_xyyzz_0, tg_zzzzzz_xyyzz_1, tg_zzzzzz_xyzz_1, tg_zzzzzz_xyzzz_0, \
                                         tg_zzzzzz_xyzzz_1, tg_zzzzzz_xzzz_1, tg_zzzzzz_xzzzz_0, tg_zzzzzz_xzzzz_1, \
                                         tg_zzzzzz_yyyy_1, tg_zzzzzz_yyyyy_0, tg_zzzzzz_yyyyy_1, tg_zzzzzz_yyyyz_0, \
                                         tg_zzzzzz_yyyyz_1, tg_zzzzzz_yyyz_1, tg_zzzzzz_yyyzz_0, tg_zzzzzz_yyyzz_1, \
                                         tg_zzzzzz_yyzz_1, tg_zzzzzz_yyzzz_0, tg_zzzzzz_yyzzz_1, tg_zzzzzz_yzzz_1, \
                                         tg_zzzzzz_yzzzz_0, tg_zzzzzz_yzzzz_1, tg_zzzzzz_zzzz_1, tg_zzzzzz_zzzzz_0, \
                                         tg_zzzzzz_zzzzz_1, tg_zzzzzzz_xxxxx_0, tg_zzzzzzz_xxxxy_0, tg_zzzzzzz_xxxxz_0, \
                                         tg_zzzzzzz_xxxyy_0, tg_zzzzzzz_xxxyz_0, tg_zzzzzzz_xxxzz_0, tg_zzzzzzz_xxyyy_0, \
                                         tg_zzzzzzz_xxyyz_0, tg_zzzzzzz_xxyzz_0, tg_zzzzzzz_xxzzz_0, tg_zzzzzzz_xyyyy_0, \
                                         tg_zzzzzzz_xyyyz_0, tg_zzzzzzz_xyyzz_0, tg_zzzzzzz_xyzzz_0, tg_zzzzzzz_xzzzz_0, \
                                         tg_zzzzzzz_yyyyy_0, tg_zzzzzzz_yyyyz_0, tg_zzzzzzz_yyyzz_0, tg_zzzzzzz_yyzzz_0, \
                                         tg_zzzzzzz_yzzzz_0, tg_zzzzzzz_zzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyzzz_xyyyz_0[j] = pb_y * tg_yyyzzz_xyyyz_0[j] + fr * tg_yyyzzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyyz_0[j] - tg_yyzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_xyyz_1[j];

                    tg_yyyyzzz_xyyzz_0[j] = pb_y * tg_yyyzzz_xyyzz_0[j] + fr * tg_yyyzzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyzz_0[j] - tg_yyzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xyzz_1[j];

                    tg_yyyyzzz_xyzzz_0[j] = pb_y * tg_yyyzzz_xyzzz_0[j] + fr * tg_yyyzzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyzzz_0[j] - tg_yyzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xzzz_1[j];

                    tg_yyyyzzz_xzzzz_0[j] = pb_y * tg_yyyzzz_xzzzz_0[j] + fr * tg_yyyzzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xzzzz_0[j] - tg_yyzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_yyyyy_0[j] = pb_y * tg_yyyzzz_yyyyy_0[j] + fr * tg_yyyzzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyyy_0[j] - tg_yyzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyzzz_yyyy_1[j];

                    tg_yyyyzzz_yyyyz_0[j] = pb_y * tg_yyyzzz_yyyyz_0[j] + fr * tg_yyyzzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyyz_0[j] - tg_yyzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzz_yyyz_1[j];

                    tg_yyyyzzz_yyyzz_0[j] = pb_y * tg_yyyzzz_yyyzz_0[j] + fr * tg_yyyzzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyzz_0[j] - tg_yyzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_yyzz_1[j];

                    tg_yyyyzzz_yyzzz_0[j] = pb_y * tg_yyyzzz_yyzzz_0[j] + fr * tg_yyyzzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyzzz_0[j] - tg_yyzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_yzzz_1[j];

                    tg_yyyyzzz_yzzzz_0[j] = pb_y * tg_yyyzzz_yzzzz_0[j] + fr * tg_yyyzzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yzzzz_0[j] - tg_yyzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_zzzz_1[j];

                    tg_yyyyzzz_zzzzz_0[j] = pb_y * tg_yyyzzz_zzzzz_0[j] + fr * tg_yyyzzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_zzzzz_0[j] - tg_yyzzz_zzzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxxx_0[j] = pb_y * tg_yyzzzz_xxxxx_0[j] + fr * tg_yyzzzz_xxxxx_1[j] + fl1_fx * (tg_yzzzz_xxxxx_0[j] - tg_yzzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxxy_0[j] = pb_y * tg_yyzzzz_xxxxy_0[j] + fr * tg_yyzzzz_xxxxy_1[j] + fl1_fx * (tg_yzzzz_xxxxy_0[j] - tg_yzzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxxx_1[j];

                    tg_yyyzzzz_xxxxz_0[j] = pb_y * tg_yyzzzz_xxxxz_0[j] + fr * tg_yyzzzz_xxxxz_1[j] + fl1_fx * (tg_yzzzz_xxxxz_0[j] - tg_yzzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxyy_0[j] = pb_y * tg_yyzzzz_xxxyy_0[j] + fr * tg_yyzzzz_xxxyy_1[j] + fl1_fx * (tg_yzzzz_xxxyy_0[j] - tg_yzzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xxxy_1[j];

                    tg_yyyzzzz_xxxyz_0[j] = pb_y * tg_yyzzzz_xxxyz_0[j] + fr * tg_yyzzzz_xxxyz_1[j] + fl1_fx * (tg_yzzzz_xxxyz_0[j] - tg_yzzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxxz_1[j];

                    tg_yyyzzzz_xxxzz_0[j] = pb_y * tg_yyzzzz_xxxzz_0[j] + fr * tg_yyzzzz_xxxzz_1[j] + fl1_fx * (tg_yzzzz_xxxzz_0[j] - tg_yzzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxyyy_0[j] = pb_y * tg_yyzzzz_xxyyy_0[j] + fr * tg_yyzzzz_xxyyy_1[j] + fl1_fx * (tg_yzzzz_xxyyy_0[j] - tg_yzzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_xxyy_1[j];

                    tg_yyyzzzz_xxyyz_0[j] = pb_y * tg_yyzzzz_xxyyz_0[j] + fr * tg_yyzzzz_xxyyz_1[j] + fl1_fx * (tg_yzzzz_xxyyz_0[j] - tg_yzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xxyz_1[j];

                    tg_yyyzzzz_xxyzz_0[j] = pb_y * tg_yyzzzz_xxyzz_0[j] + fr * tg_yyzzzz_xxyzz_1[j] + fl1_fx * (tg_yzzzz_xxyzz_0[j] - tg_yzzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxzz_1[j];

                    tg_yyyzzzz_xxzzz_0[j] = pb_y * tg_yyzzzz_xxzzz_0[j] + fr * tg_yyzzzz_xxzzz_1[j] + fl1_fx * (tg_yzzzz_xxzzz_0[j] - tg_yzzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xyyyy_0[j] = pb_y * tg_yyzzzz_xyyyy_0[j] + fr * tg_yyzzzz_xyyyy_1[j] + fl1_fx * (tg_yzzzz_xyyyy_0[j] - tg_yzzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzz_xyyy_1[j];

                    tg_yyyzzzz_xyyyz_0[j] = pb_y * tg_yyzzzz_xyyyz_0[j] + fr * tg_yyzzzz_xyyyz_1[j] + fl1_fx * (tg_yzzzz_xyyyz_0[j] - tg_yzzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_xyyz_1[j];

                    tg_yyyzzzz_xyyzz_0[j] = pb_y * tg_yyzzzz_xyyzz_0[j] + fr * tg_yyzzzz_xyyzz_1[j] + fl1_fx * (tg_yzzzz_xyyzz_0[j] - tg_yzzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xyzz_1[j];

                    tg_yyyzzzz_xyzzz_0[j] = pb_y * tg_yyzzzz_xyzzz_0[j] + fr * tg_yyzzzz_xyzzz_1[j] + fl1_fx * (tg_yzzzz_xyzzz_0[j] - tg_yzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xzzz_1[j];

                    tg_yyyzzzz_xzzzz_0[j] = pb_y * tg_yyzzzz_xzzzz_0[j] + fr * tg_yyzzzz_xzzzz_1[j] + fl1_fx * (tg_yzzzz_xzzzz_0[j] - tg_yzzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_yyyyy_0[j] = pb_y * tg_yyzzzz_yyyyy_0[j] + fr * tg_yyzzzz_yyyyy_1[j] + fl1_fx * (tg_yzzzz_yyyyy_0[j] - tg_yzzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzzzz_yyyy_1[j];

                    tg_yyyzzzz_yyyyz_0[j] = pb_y * tg_yyzzzz_yyyyz_0[j] + fr * tg_yyzzzz_yyyyz_1[j] + fl1_fx * (tg_yzzzz_yyyyz_0[j] - tg_yzzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzz_yyyz_1[j];

                    tg_yyyzzzz_yyyzz_0[j] = pb_y * tg_yyzzzz_yyyzz_0[j] + fr * tg_yyzzzz_yyyzz_1[j] + fl1_fx * (tg_yzzzz_yyyzz_0[j] - tg_yzzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_yyzz_1[j];

                    tg_yyyzzzz_yyzzz_0[j] = pb_y * tg_yyzzzz_yyzzz_0[j] + fr * tg_yyzzzz_yyzzz_1[j] + fl1_fx * (tg_yzzzz_yyzzz_0[j] - tg_yzzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_yzzz_1[j];

                    tg_yyyzzzz_yzzzz_0[j] = pb_y * tg_yyzzzz_yzzzz_0[j] + fr * tg_yyzzzz_yzzzz_1[j] + fl1_fx * (tg_yzzzz_yzzzz_0[j] - tg_yzzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_zzzz_1[j];

                    tg_yyyzzzz_zzzzz_0[j] = pb_y * tg_yyzzzz_zzzzz_0[j] + fr * tg_yyzzzz_zzzzz_1[j] + fl1_fx * (tg_yzzzz_zzzzz_0[j] - tg_yzzzz_zzzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxxx_0[j] = pb_y * tg_yzzzzz_xxxxx_0[j] + fr * tg_yzzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxx_0[j] - tg_zzzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxxy_0[j] = pb_y * tg_yzzzzz_xxxxy_0[j] + fr * tg_yzzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxy_0[j] - tg_zzzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxxx_1[j];

                    tg_yyzzzzz_xxxxz_0[j] = pb_y * tg_yzzzzz_xxxxz_0[j] + fr * tg_yzzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxz_0[j] - tg_zzzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxyy_0[j] = pb_y * tg_yzzzzz_xxxyy_0[j] + fr * tg_yzzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyy_0[j] - tg_zzzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xxxy_1[j];

                    tg_yyzzzzz_xxxyz_0[j] = pb_y * tg_yzzzzz_xxxyz_0[j] + fr * tg_yzzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyz_0[j] - tg_zzzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxxz_1[j];

                    tg_yyzzzzz_xxxzz_0[j] = pb_y * tg_yzzzzz_xxxzz_0[j] + fr * tg_yzzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxzz_0[j] - tg_zzzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxyyy_0[j] = pb_y * tg_yzzzzz_xxyyy_0[j] + fr * tg_yzzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyy_0[j] - tg_zzzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_xxyy_1[j];

                    tg_yyzzzzz_xxyyz_0[j] = pb_y * tg_yzzzzz_xxyyz_0[j] + fr * tg_yzzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyz_0[j] - tg_zzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xxyz_1[j];

                    tg_yyzzzzz_xxyzz_0[j] = pb_y * tg_yzzzzz_xxyzz_0[j] + fr * tg_yzzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyzz_0[j] - tg_zzzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxzz_1[j];

                    tg_yyzzzzz_xxzzz_0[j] = pb_y * tg_yzzzzz_xxzzz_0[j] + fr * tg_yzzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxzzz_0[j] - tg_zzzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xyyyy_0[j] = pb_y * tg_yzzzzz_xyyyy_0[j] + fr * tg_yzzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyy_0[j] - tg_zzzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzz_xyyy_1[j];

                    tg_yyzzzzz_xyyyz_0[j] = pb_y * tg_yzzzzz_xyyyz_0[j] + fr * tg_yzzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyz_0[j] - tg_zzzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_xyyz_1[j];

                    tg_yyzzzzz_xyyzz_0[j] = pb_y * tg_yzzzzz_xyyzz_0[j] + fr * tg_yzzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyzz_0[j] - tg_zzzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xyzz_1[j];

                    tg_yyzzzzz_xyzzz_0[j] = pb_y * tg_yzzzzz_xyzzz_0[j] + fr * tg_yzzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyzzz_0[j] - tg_zzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xzzz_1[j];

                    tg_yyzzzzz_xzzzz_0[j] = pb_y * tg_yzzzzz_xzzzz_0[j] + fr * tg_yzzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzzzz_0[j] - tg_zzzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_yyyyy_0[j] = pb_y * tg_yzzzzz_yyyyy_0[j] + fr * tg_yzzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyy_0[j] - tg_zzzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzzzz_yyyy_1[j];

                    tg_yyzzzzz_yyyyz_0[j] = pb_y * tg_yzzzzz_yyyyz_0[j] + fr * tg_yzzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyz_0[j] - tg_zzzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzz_yyyz_1[j];

                    tg_yyzzzzz_yyyzz_0[j] = pb_y * tg_yzzzzz_yyyzz_0[j] + fr * tg_yzzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyzz_0[j] - tg_zzzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_yyzz_1[j];

                    tg_yyzzzzz_yyzzz_0[j] = pb_y * tg_yzzzzz_yyzzz_0[j] + fr * tg_yzzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyzzz_0[j] - tg_zzzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_yzzz_1[j];

                    tg_yyzzzzz_yzzzz_0[j] = pb_y * tg_yzzzzz_yzzzz_0[j] + fr * tg_yzzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzzzz_0[j] - tg_zzzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_zzzz_1[j];

                    tg_yyzzzzz_zzzzz_0[j] = pb_y * tg_yzzzzz_zzzzz_0[j] + fr * tg_yzzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzzzz_0[j] - tg_zzzzz_zzzzz_1[j] * fl1_fza);

                    tg_yzzzzzz_xxxxx_0[j] = pb_y * tg_zzzzzz_xxxxx_0[j] + fr * tg_zzzzzz_xxxxx_1[j];

                    tg_yzzzzzz_xxxxy_0[j] = pb_y * tg_zzzzzz_xxxxy_0[j] + fr * tg_zzzzzz_xxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxxx_1[j];

                    tg_yzzzzzz_xxxxz_0[j] = pb_y * tg_zzzzzz_xxxxz_0[j] + fr * tg_zzzzzz_xxxxz_1[j];

                    tg_yzzzzzz_xxxyy_0[j] = pb_y * tg_zzzzzz_xxxyy_0[j] + fr * tg_zzzzzz_xxxyy_1[j] + fl1_fxn * tg_zzzzzz_xxxy_1[j];

                    tg_yzzzzzz_xxxyz_0[j] = pb_y * tg_zzzzzz_xxxyz_0[j] + fr * tg_zzzzzz_xxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxxz_1[j];

                    tg_yzzzzzz_xxxzz_0[j] = pb_y * tg_zzzzzz_xxxzz_0[j] + fr * tg_zzzzzz_xxxzz_1[j];

                    tg_yzzzzzz_xxyyy_0[j] = pb_y * tg_zzzzzz_xxyyy_0[j] + fr * tg_zzzzzz_xxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxyy_1[j];

                    tg_yzzzzzz_xxyyz_0[j] = pb_y * tg_zzzzzz_xxyyz_0[j] + fr * tg_zzzzzz_xxyyz_1[j] + fl1_fxn * tg_zzzzzz_xxyz_1[j];

                    tg_yzzzzzz_xxyzz_0[j] = pb_y * tg_zzzzzz_xxyzz_0[j] + fr * tg_zzzzzz_xxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxzz_1[j];

                    tg_yzzzzzz_xxzzz_0[j] = pb_y * tg_zzzzzz_xxzzz_0[j] + fr * tg_zzzzzz_xxzzz_1[j];

                    tg_yzzzzzz_xyyyy_0[j] = pb_y * tg_zzzzzz_xyyyy_0[j] + fr * tg_zzzzzz_xyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xyyy_1[j];

                    tg_yzzzzzz_xyyyz_0[j] = pb_y * tg_zzzzzz_xyyyz_0[j] + fr * tg_zzzzzz_xyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xyyz_1[j];

                    tg_yzzzzzz_xyyzz_0[j] = pb_y * tg_zzzzzz_xyyzz_0[j] + fr * tg_zzzzzz_xyyzz_1[j] + fl1_fxn * tg_zzzzzz_xyzz_1[j];

                    tg_yzzzzzz_xyzzz_0[j] = pb_y * tg_zzzzzz_xyzzz_0[j] + fr * tg_zzzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xzzz_1[j];

                    tg_yzzzzzz_xzzzz_0[j] = pb_y * tg_zzzzzz_xzzzz_0[j] + fr * tg_zzzzzz_xzzzz_1[j];

                    tg_yzzzzzz_yyyyy_0[j] = pb_y * tg_zzzzzz_yyyyy_0[j] + fr * tg_zzzzzz_yyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzzzz_yyyy_1[j];

                    tg_yzzzzzz_yyyyz_0[j] = pb_y * tg_zzzzzz_yyyyz_0[j] + fr * tg_zzzzzz_yyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_yyyz_1[j];

                    tg_yzzzzzz_yyyzz_0[j] = pb_y * tg_zzzzzz_yyyzz_0[j] + fr * tg_zzzzzz_yyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_yyzz_1[j];

                    tg_yzzzzzz_yyzzz_0[j] = pb_y * tg_zzzzzz_yyzzz_0[j] + fr * tg_zzzzzz_yyzzz_1[j] + fl1_fxn * tg_zzzzzz_yzzz_1[j];

                    tg_yzzzzzz_yzzzz_0[j] = pb_y * tg_zzzzzz_yzzzz_0[j] + fr * tg_zzzzzz_yzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zzzz_1[j];

                    tg_yzzzzzz_zzzzz_0[j] = pb_y * tg_zzzzzz_zzzzz_0[j] + fr * tg_zzzzzz_zzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzzz_xxxxx_0[j] = pb_z * tg_zzzzzz_xxxxx_0[j] + fr * tg_zzzzzz_xxxxx_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxx_0[j] - tg_zzzzz_xxxxx_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxxy_0[j] = pb_z * tg_zzzzzz_xxxxy_0[j] + fr * tg_zzzzzz_xxxxy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxy_0[j] - tg_zzzzz_xxxxy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxxz_0[j] = pb_z * tg_zzzzzz_xxxxz_0[j] + fr * tg_zzzzzz_xxxxz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxz_0[j] - tg_zzzzz_xxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxxx_1[j];

                    tg_zzzzzzz_xxxyy_0[j] = pb_z * tg_zzzzzz_xxxyy_0[j] + fr * tg_zzzzzz_xxxyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxyy_0[j] - tg_zzzzz_xxxyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxyz_0[j] = pb_z * tg_zzzzzz_xxxyz_0[j] + fr * tg_zzzzzz_xxxyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxyz_0[j] - tg_zzzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxxy_1[j];

                    tg_zzzzzzz_xxxzz_0[j] = pb_z * tg_zzzzzz_xxxzz_0[j] + fr * tg_zzzzzz_xxxzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxzz_0[j] - tg_zzzzz_xxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xxxz_1[j];

                    tg_zzzzzzz_xxyyy_0[j] = pb_z * tg_zzzzzz_xxyyy_0[j] + fr * tg_zzzzzz_xxyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyyy_0[j] - tg_zzzzz_xxyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxyyz_0[j] = pb_z * tg_zzzzzz_xxyyz_0[j] + fr * tg_zzzzzz_xxyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyyz_0[j] - tg_zzzzz_xxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxyy_1[j];

                    tg_zzzzzzz_xxyzz_0[j] = pb_z * tg_zzzzzz_xxyzz_0[j] + fr * tg_zzzzzz_xxyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyzz_0[j] - tg_zzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xxyz_1[j];

                    tg_zzzzzzz_xxzzz_0[j] = pb_z * tg_zzzzzz_xxzzz_0[j] + fr * tg_zzzzzz_xxzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxzzz_0[j] - tg_zzzzz_xxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_xxzz_1[j];

                    tg_zzzzzzz_xyyyy_0[j] = pb_z * tg_zzzzzz_xyyyy_0[j] + fr * tg_zzzzzz_xyyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyyy_0[j] - tg_zzzzz_xyyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xyyyz_0[j] = pb_z * tg_zzzzzz_xyyyz_0[j] + fr * tg_zzzzzz_xyyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyyz_0[j] - tg_zzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xyyy_1[j];

                    tg_zzzzzzz_xyyzz_0[j] = pb_z * tg_zzzzzz_xyyzz_0[j] + fr * tg_zzzzzz_xyyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyzz_0[j] - tg_zzzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xyyz_1[j];

                    tg_zzzzzzz_xyzzz_0[j] = pb_z * tg_zzzzzz_xyzzz_0[j] + fr * tg_zzzzzz_xyzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyzzz_0[j] - tg_zzzzz_xyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_xyzz_1[j];

                    tg_zzzzzzz_xzzzz_0[j] = pb_z * tg_zzzzzz_xzzzz_0[j] + fr * tg_zzzzzz_xzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xzzzz_0[j] - tg_zzzzz_xzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzz_xzzz_1[j];

                    tg_zzzzzzz_yyyyy_0[j] = pb_z * tg_zzzzzz_yyyyy_0[j] + fr * tg_zzzzzz_yyyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyyy_0[j] - tg_zzzzz_yyyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_yyyyz_0[j] = pb_z * tg_zzzzzz_yyyyz_0[j] + fr * tg_zzzzzz_yyyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyyz_0[j] - tg_zzzzz_yyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_yyyy_1[j];

                    tg_zzzzzzz_yyyzz_0[j] = pb_z * tg_zzzzzz_yyyzz_0[j] + fr * tg_zzzzzz_yyyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyzz_0[j] - tg_zzzzz_yyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_yyyz_1[j];

                    tg_zzzzzzz_yyzzz_0[j] = pb_z * tg_zzzzzz_yyzzz_0[j] + fr * tg_zzzzzz_yyzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyzzz_0[j] - tg_zzzzz_yyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_yyzz_1[j];

                    tg_zzzzzzz_yzzzz_0[j] = pb_z * tg_zzzzzz_yzzzz_0[j] + fr * tg_zzzzzz_yzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yzzzz_0[j] - tg_zzzzz_yzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzz_yzzz_1[j];

                    tg_zzzzzzz_zzzzz_0[j] = pb_z * tg_zzzzzz_zzzzz_0[j] + fr * tg_zzzzzz_zzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_zzzzz_0[j] - tg_zzzzz_zzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzzzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

