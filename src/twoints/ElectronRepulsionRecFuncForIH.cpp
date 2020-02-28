//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForIH.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSISH(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSISH_0_98(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSISH_98_196(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSISH_196_294(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISH_294_392(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISH_392_490(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSISH_490_588(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSISH_0_98(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxxyz_xyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 95); 

                auto tg_xxxyz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 96); 

                auto tg_xxxyz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 97); 

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

                auto tg_xxxyz_xyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 95); 

                auto tg_xxxyz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 96); 

                auto tg_xxxyz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 97); 

                auto tg_xxxx_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx); 

                auto tg_xxxx_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 1); 

                auto tg_xxxx_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 2); 

                auto tg_xxxx_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 3); 

                auto tg_xxxx_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 4); 

                auto tg_xxxx_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 5); 

                auto tg_xxxx_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 6); 

                auto tg_xxxx_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 7); 

                auto tg_xxxx_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 8); 

                auto tg_xxxx_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 9); 

                auto tg_xxxx_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 10); 

                auto tg_xxxx_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 11); 

                auto tg_xxxx_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 12); 

                auto tg_xxxx_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 13); 

                auto tg_xxxx_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 14); 

                auto tg_xxxx_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 15); 

                auto tg_xxxx_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 16); 

                auto tg_xxxx_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 17); 

                auto tg_xxxx_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 18); 

                auto tg_xxxx_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 19); 

                auto tg_xxxx_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 20); 

                auto tg_xxxy_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 21); 

                auto tg_xxxy_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 22); 

                auto tg_xxxy_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 23); 

                auto tg_xxxy_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 24); 

                auto tg_xxxy_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 25); 

                auto tg_xxxy_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 26); 

                auto tg_xxxy_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 27); 

                auto tg_xxxy_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 28); 

                auto tg_xxxy_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 29); 

                auto tg_xxxy_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 30); 

                auto tg_xxxy_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 31); 

                auto tg_xxxy_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 32); 

                auto tg_xxxy_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 33); 

                auto tg_xxxy_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 34); 

                auto tg_xxxy_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 35); 

                auto tg_xxxy_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 36); 

                auto tg_xxxy_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 37); 

                auto tg_xxxy_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 38); 

                auto tg_xxxy_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 39); 

                auto tg_xxxy_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 40); 

                auto tg_xxxy_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 41); 

                auto tg_xxxz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 42); 

                auto tg_xxxz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 43); 

                auto tg_xxxz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 44); 

                auto tg_xxxz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 45); 

                auto tg_xxxz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 46); 

                auto tg_xxxz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 47); 

                auto tg_xxxz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 48); 

                auto tg_xxxz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 49); 

                auto tg_xxxz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 50); 

                auto tg_xxxz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 51); 

                auto tg_xxxz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 52); 

                auto tg_xxxz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 53); 

                auto tg_xxxz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 54); 

                auto tg_xxxz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 55); 

                auto tg_xxxz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 56); 

                auto tg_xxxz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 57); 

                auto tg_xxxz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 58); 

                auto tg_xxxz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 59); 

                auto tg_xxxz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 60); 

                auto tg_xxxz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 61); 

                auto tg_xxxz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 62); 

                auto tg_xxyy_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 63); 

                auto tg_xxyy_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 64); 

                auto tg_xxyy_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 65); 

                auto tg_xxyy_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 66); 

                auto tg_xxyy_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 67); 

                auto tg_xxyy_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 68); 

                auto tg_xxyy_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 69); 

                auto tg_xxyy_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 70); 

                auto tg_xxyy_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 71); 

                auto tg_xxyy_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 72); 

                auto tg_xxyy_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 73); 

                auto tg_xxyy_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 74); 

                auto tg_xxyy_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 75); 

                auto tg_xxyy_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 76); 

                auto tg_xxyy_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 77); 

                auto tg_xxyy_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 78); 

                auto tg_xxyy_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 79); 

                auto tg_xxyy_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 80); 

                auto tg_xxyy_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 81); 

                auto tg_xxyy_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 82); 

                auto tg_xxyy_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 83); 

                auto tg_xxyz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 84); 

                auto tg_xxyz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 85); 

                auto tg_xxyz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 86); 

                auto tg_xxyz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 87); 

                auto tg_xxyz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 88); 

                auto tg_xxyz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 89); 

                auto tg_xxyz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 90); 

                auto tg_xxyz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 91); 

                auto tg_xxyz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 92); 

                auto tg_xxyz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 93); 

                auto tg_xxyz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 94); 

                auto tg_xxyz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 95); 

                auto tg_xxyz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 96); 

                auto tg_xxyz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 97); 

                auto tg_xxxx_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx); 

                auto tg_xxxx_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 1); 

                auto tg_xxxx_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 2); 

                auto tg_xxxx_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 3); 

                auto tg_xxxx_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 4); 

                auto tg_xxxx_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 5); 

                auto tg_xxxx_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 6); 

                auto tg_xxxx_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 7); 

                auto tg_xxxx_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 8); 

                auto tg_xxxx_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 9); 

                auto tg_xxxx_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 10); 

                auto tg_xxxx_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 11); 

                auto tg_xxxx_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 12); 

                auto tg_xxxx_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 13); 

                auto tg_xxxx_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 14); 

                auto tg_xxxx_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 15); 

                auto tg_xxxx_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 16); 

                auto tg_xxxx_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 17); 

                auto tg_xxxx_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 18); 

                auto tg_xxxx_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 19); 

                auto tg_xxxx_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 20); 

                auto tg_xxxy_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 21); 

                auto tg_xxxy_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 22); 

                auto tg_xxxy_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 23); 

                auto tg_xxxy_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 24); 

                auto tg_xxxy_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 25); 

                auto tg_xxxy_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 26); 

                auto tg_xxxy_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 27); 

                auto tg_xxxy_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 28); 

                auto tg_xxxy_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 29); 

                auto tg_xxxy_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 30); 

                auto tg_xxxy_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 31); 

                auto tg_xxxy_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 32); 

                auto tg_xxxy_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 33); 

                auto tg_xxxy_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 34); 

                auto tg_xxxy_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 35); 

                auto tg_xxxy_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 36); 

                auto tg_xxxy_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 37); 

                auto tg_xxxy_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 38); 

                auto tg_xxxy_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 39); 

                auto tg_xxxy_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 40); 

                auto tg_xxxy_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 41); 

                auto tg_xxxz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 42); 

                auto tg_xxxz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 43); 

                auto tg_xxxz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 44); 

                auto tg_xxxz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 45); 

                auto tg_xxxz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 46); 

                auto tg_xxxz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 47); 

                auto tg_xxxz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 48); 

                auto tg_xxxz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 49); 

                auto tg_xxxz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 50); 

                auto tg_xxxz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 51); 

                auto tg_xxxz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 52); 

                auto tg_xxxz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 53); 

                auto tg_xxxz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 54); 

                auto tg_xxxz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 55); 

                auto tg_xxxz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 56); 

                auto tg_xxxz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 57); 

                auto tg_xxxz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 58); 

                auto tg_xxxz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 59); 

                auto tg_xxxz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 60); 

                auto tg_xxxz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 61); 

                auto tg_xxxz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 62); 

                auto tg_xxyy_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 63); 

                auto tg_xxyy_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 64); 

                auto tg_xxyy_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 65); 

                auto tg_xxyy_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 66); 

                auto tg_xxyy_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 67); 

                auto tg_xxyy_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 68); 

                auto tg_xxyy_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 69); 

                auto tg_xxyy_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 70); 

                auto tg_xxyy_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 71); 

                auto tg_xxyy_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 72); 

                auto tg_xxyy_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 73); 

                auto tg_xxyy_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 74); 

                auto tg_xxyy_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 75); 

                auto tg_xxyy_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 76); 

                auto tg_xxyy_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 77); 

                auto tg_xxyy_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 78); 

                auto tg_xxyy_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 79); 

                auto tg_xxyy_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 80); 

                auto tg_xxyy_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 81); 

                auto tg_xxyy_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 82); 

                auto tg_xxyy_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 83); 

                auto tg_xxyz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 84); 

                auto tg_xxyz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 85); 

                auto tg_xxyz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 86); 

                auto tg_xxyz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 87); 

                auto tg_xxyz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 88); 

                auto tg_xxyz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 89); 

                auto tg_xxyz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 90); 

                auto tg_xxyz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 91); 

                auto tg_xxyz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 92); 

                auto tg_xxyz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 93); 

                auto tg_xxyz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 94); 

                auto tg_xxyz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 95); 

                auto tg_xxyz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 96); 

                auto tg_xxyz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 97); 

                auto tg_xxxxx_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx); 

                auto tg_xxxxx_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 1); 

                auto tg_xxxxx_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 2); 

                auto tg_xxxxx_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 3); 

                auto tg_xxxxx_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 4); 

                auto tg_xxxxx_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 5); 

                auto tg_xxxxx_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 6); 

                auto tg_xxxxx_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 7); 

                auto tg_xxxxx_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 8); 

                auto tg_xxxxx_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 9); 

                auto tg_xxxxx_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 10); 

                auto tg_xxxxx_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 11); 

                auto tg_xxxxx_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 12); 

                auto tg_xxxxx_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 13); 

                auto tg_xxxxx_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 14); 

                auto tg_xxxxy_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 15); 

                auto tg_xxxxy_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 16); 

                auto tg_xxxxy_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 17); 

                auto tg_xxxxy_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 18); 

                auto tg_xxxxy_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 19); 

                auto tg_xxxxy_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 20); 

                auto tg_xxxxy_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 21); 

                auto tg_xxxxy_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 22); 

                auto tg_xxxxy_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 23); 

                auto tg_xxxxy_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 24); 

                auto tg_xxxxy_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 25); 

                auto tg_xxxxy_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 26); 

                auto tg_xxxxy_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 27); 

                auto tg_xxxxy_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 28); 

                auto tg_xxxxy_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 29); 

                auto tg_xxxxz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 30); 

                auto tg_xxxxz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 31); 

                auto tg_xxxxz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 32); 

                auto tg_xxxxz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 33); 

                auto tg_xxxxz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 34); 

                auto tg_xxxxz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 35); 

                auto tg_xxxxz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 36); 

                auto tg_xxxxz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 37); 

                auto tg_xxxxz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 38); 

                auto tg_xxxxz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 39); 

                auto tg_xxxxz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 40); 

                auto tg_xxxxz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 41); 

                auto tg_xxxxz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 42); 

                auto tg_xxxxz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 43); 

                auto tg_xxxxz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 44); 

                auto tg_xxxyy_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 45); 

                auto tg_xxxyy_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 46); 

                auto tg_xxxyy_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 47); 

                auto tg_xxxyy_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 48); 

                auto tg_xxxyy_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 49); 

                auto tg_xxxyy_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 50); 

                auto tg_xxxyy_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 51); 

                auto tg_xxxyy_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 52); 

                auto tg_xxxyy_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 53); 

                auto tg_xxxyy_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 54); 

                auto tg_xxxyy_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 55); 

                auto tg_xxxyy_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 56); 

                auto tg_xxxyy_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 57); 

                auto tg_xxxyy_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 58); 

                auto tg_xxxyy_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 59); 

                auto tg_xxxyz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 60); 

                auto tg_xxxyz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 61); 

                auto tg_xxxyz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 62); 

                auto tg_xxxyz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 63); 

                auto tg_xxxyz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 64); 

                auto tg_xxxyz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 65); 

                auto tg_xxxyz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 66); 

                auto tg_xxxyz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 67); 

                auto tg_xxxyz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 68); 

                auto tg_xxxyz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 69); 

                auto tg_xxxyz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 70); 

                auto tg_xxxyz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 71); 

                auto tg_xxxyz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 72); 

                auto tg_xxxyz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 73); 

                // set up pointers to integrals

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

                auto tg_xxxxyz_xyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 95); 

                auto tg_xxxxyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 96); 

                auto tg_xxxxyz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 97); 

                // Batch of Integrals (0,98)

                #pragma omp simd aligned(fxn, fza, tg_xxxx_xxxxx_0, tg_xxxx_xxxxx_1, tg_xxxx_xxxxy_0, \
                                         tg_xxxx_xxxxy_1, tg_xxxx_xxxxz_0, tg_xxxx_xxxxz_1, tg_xxxx_xxxyy_0, tg_xxxx_xxxyy_1, \
                                         tg_xxxx_xxxyz_0, tg_xxxx_xxxyz_1, tg_xxxx_xxxzz_0, tg_xxxx_xxxzz_1, tg_xxxx_xxyyy_0, \
                                         tg_xxxx_xxyyy_1, tg_xxxx_xxyyz_0, tg_xxxx_xxyyz_1, tg_xxxx_xxyzz_0, tg_xxxx_xxyzz_1, \
                                         tg_xxxx_xxzzz_0, tg_xxxx_xxzzz_1, tg_xxxx_xyyyy_0, tg_xxxx_xyyyy_1, tg_xxxx_xyyyz_0, \
                                         tg_xxxx_xyyyz_1, tg_xxxx_xyyzz_0, tg_xxxx_xyyzz_1, tg_xxxx_xyzzz_0, tg_xxxx_xyzzz_1, \
                                         tg_xxxx_xzzzz_0, tg_xxxx_xzzzz_1, tg_xxxx_yyyyy_0, tg_xxxx_yyyyy_1, tg_xxxx_yyyyz_0, \
                                         tg_xxxx_yyyyz_1, tg_xxxx_yyyzz_0, tg_xxxx_yyyzz_1, tg_xxxx_yyzzz_0, tg_xxxx_yyzzz_1, \
                                         tg_xxxx_yzzzz_0, tg_xxxx_yzzzz_1, tg_xxxx_zzzzz_0, tg_xxxx_zzzzz_1, tg_xxxxx_xxxx_1, \
                                         tg_xxxxx_xxxxx_0, tg_xxxxx_xxxxx_1, tg_xxxxx_xxxxy_0, tg_xxxxx_xxxxy_1, \
                                         tg_xxxxx_xxxxz_0, tg_xxxxx_xxxxz_1, tg_xxxxx_xxxy_1, tg_xxxxx_xxxyy_0, \
                                         tg_xxxxx_xxxyy_1, tg_xxxxx_xxxyz_0, tg_xxxxx_xxxyz_1, tg_xxxxx_xxxz_1, \
                                         tg_xxxxx_xxxzz_0, tg_xxxxx_xxxzz_1, tg_xxxxx_xxyy_1, tg_xxxxx_xxyyy_0, \
                                         tg_xxxxx_xxyyy_1, tg_xxxxx_xxyyz_0, tg_xxxxx_xxyyz_1, tg_xxxxx_xxyz_1, \
                                         tg_xxxxx_xxyzz_0, tg_xxxxx_xxyzz_1, tg_xxxxx_xxzz_1, tg_xxxxx_xxzzz_0, \
                                         tg_xxxxx_xxzzz_1, tg_xxxxx_xyyy_1, tg_xxxxx_xyyyy_0, tg_xxxxx_xyyyy_1, \
                                         tg_xxxxx_xyyyz_0, tg_xxxxx_xyyyz_1, tg_xxxxx_xyyz_1, tg_xxxxx_xyyzz_0, \
                                         tg_xxxxx_xyyzz_1, tg_xxxxx_xyzz_1, tg_xxxxx_xyzzz_0, tg_xxxxx_xyzzz_1, \
                                         tg_xxxxx_xzzz_1, tg_xxxxx_xzzzz_0, tg_xxxxx_xzzzz_1, tg_xxxxx_yyyy_1, \
                                         tg_xxxxx_yyyyy_0, tg_xxxxx_yyyyy_1, tg_xxxxx_yyyyz_0, tg_xxxxx_yyyyz_1, \
                                         tg_xxxxx_yyyz_1, tg_xxxxx_yyyzz_0, tg_xxxxx_yyyzz_1, tg_xxxxx_yyzz_1, \
                                         tg_xxxxx_yyzzz_0, tg_xxxxx_yyzzz_1, tg_xxxxx_yzzz_1, tg_xxxxx_yzzzz_0, \
                                         tg_xxxxx_yzzzz_1, tg_xxxxx_zzzz_1, tg_xxxxx_zzzzz_0, tg_xxxxx_zzzzz_1, \
                                         tg_xxxxxx_xxxxx_0, tg_xxxxxx_xxxxy_0, tg_xxxxxx_xxxxz_0, tg_xxxxxx_xxxyy_0, \
                                         tg_xxxxxx_xxxyz_0, tg_xxxxxx_xxxzz_0, tg_xxxxxx_xxyyy_0, tg_xxxxxx_xxyyz_0, \
                                         tg_xxxxxx_xxyzz_0, tg_xxxxxx_xxzzz_0, tg_xxxxxx_xyyyy_0, tg_xxxxxx_xyyyz_0, \
                                         tg_xxxxxx_xyyzz_0, tg_xxxxxx_xyzzz_0, tg_xxxxxx_xzzzz_0, tg_xxxxxx_yyyyy_0, \
                                         tg_xxxxxx_yyyyz_0, tg_xxxxxx_yyyzz_0, tg_xxxxxx_yyzzz_0, tg_xxxxxx_yzzzz_0, \
                                         tg_xxxxxx_zzzzz_0, tg_xxxxxy_xxxxx_0, tg_xxxxxy_xxxxy_0, tg_xxxxxy_xxxxz_0, \
                                         tg_xxxxxy_xxxyy_0, tg_xxxxxy_xxxyz_0, tg_xxxxxy_xxxzz_0, tg_xxxxxy_xxyyy_0, \
                                         tg_xxxxxy_xxyyz_0, tg_xxxxxy_xxyzz_0, tg_xxxxxy_xxzzz_0, tg_xxxxxy_xyyyy_0, \
                                         tg_xxxxxy_xyyyz_0, tg_xxxxxy_xyyzz_0, tg_xxxxxy_xyzzz_0, tg_xxxxxy_xzzzz_0, \
                                         tg_xxxxxy_yyyyy_0, tg_xxxxxy_yyyyz_0, tg_xxxxxy_yyyzz_0, tg_xxxxxy_yyzzz_0, \
                                         tg_xxxxxy_yzzzz_0, tg_xxxxxy_zzzzz_0, tg_xxxxxz_xxxxx_0, tg_xxxxxz_xxxxy_0, \
                                         tg_xxxxxz_xxxxz_0, tg_xxxxxz_xxxyy_0, tg_xxxxxz_xxxyz_0, tg_xxxxxz_xxxzz_0, \
                                         tg_xxxxxz_xxyyy_0, tg_xxxxxz_xxyyz_0, tg_xxxxxz_xxyzz_0, tg_xxxxxz_xxzzz_0, \
                                         tg_xxxxxz_xyyyy_0, tg_xxxxxz_xyyyz_0, tg_xxxxxz_xyyzz_0, tg_xxxxxz_xyzzz_0, \
                                         tg_xxxxxz_xzzzz_0, tg_xxxxxz_yyyyy_0, tg_xxxxxz_yyyyz_0, tg_xxxxxz_yyyzz_0, \
                                         tg_xxxxxz_yyzzz_0, tg_xxxxxz_yzzzz_0, tg_xxxxxz_zzzzz_0, tg_xxxxy_xxxx_1, \
                                         tg_xxxxy_xxxxx_0, tg_xxxxy_xxxxx_1, tg_xxxxy_xxxxy_0, tg_xxxxy_xxxxy_1, \
                                         tg_xxxxy_xxxxz_0, tg_xxxxy_xxxxz_1, tg_xxxxy_xxxy_1, tg_xxxxy_xxxyy_0, \
                                         tg_xxxxy_xxxyy_1, tg_xxxxy_xxxyz_0, tg_xxxxy_xxxyz_1, tg_xxxxy_xxxz_1, \
                                         tg_xxxxy_xxxzz_0, tg_xxxxy_xxxzz_1, tg_xxxxy_xxyy_1, tg_xxxxy_xxyyy_0, \
                                         tg_xxxxy_xxyyy_1, tg_xxxxy_xxyyz_0, tg_xxxxy_xxyyz_1, tg_xxxxy_xxyz_1, \
                                         tg_xxxxy_xxyzz_0, tg_xxxxy_xxyzz_1, tg_xxxxy_xxzz_1, tg_xxxxy_xxzzz_0, \
                                         tg_xxxxy_xxzzz_1, tg_xxxxy_xyyy_1, tg_xxxxy_xyyyy_0, tg_xxxxy_xyyyy_1, \
                                         tg_xxxxy_xyyyz_0, tg_xxxxy_xyyyz_1, tg_xxxxy_xyyz_1, tg_xxxxy_xyyzz_0, \
                                         tg_xxxxy_xyyzz_1, tg_xxxxy_xyzz_1, tg_xxxxy_xyzzz_0, tg_xxxxy_xyzzz_1, \
                                         tg_xxxxy_xzzz_1, tg_xxxxy_xzzzz_0, tg_xxxxy_xzzzz_1, tg_xxxxy_yyyy_1, \
                                         tg_xxxxy_yyyyy_0, tg_xxxxy_yyyyy_1, tg_xxxxy_yyyyz_0, tg_xxxxy_yyyyz_1, \
                                         tg_xxxxy_yyyz_1, tg_xxxxy_yyyzz_0, tg_xxxxy_yyyzz_1, tg_xxxxy_yyzz_1, \
                                         tg_xxxxy_yyzzz_0, tg_xxxxy_yyzzz_1, tg_xxxxy_yzzz_1, tg_xxxxy_yzzzz_0, \
                                         tg_xxxxy_yzzzz_1, tg_xxxxy_zzzz_1, tg_xxxxy_zzzzz_0, tg_xxxxy_zzzzz_1, \
                                         tg_xxxxyy_xxxxx_0, tg_xxxxyy_xxxxy_0, tg_xxxxyy_xxxxz_0, tg_xxxxyy_xxxyy_0, \
                                         tg_xxxxyy_xxxyz_0, tg_xxxxyy_xxxzz_0, tg_xxxxyy_xxyyy_0, tg_xxxxyy_xxyyz_0, \
                                         tg_xxxxyy_xxyzz_0, tg_xxxxyy_xxzzz_0, tg_xxxxyy_xyyyy_0, tg_xxxxyy_xyyyz_0, \
                                         tg_xxxxyy_xyyzz_0, tg_xxxxyy_xyzzz_0, tg_xxxxyy_xzzzz_0, tg_xxxxyy_yyyyy_0, \
                                         tg_xxxxyy_yyyyz_0, tg_xxxxyy_yyyzz_0, tg_xxxxyy_yyzzz_0, tg_xxxxyy_yzzzz_0, \
                                         tg_xxxxyy_zzzzz_0, tg_xxxxyz_xxxxx_0, tg_xxxxyz_xxxxy_0, tg_xxxxyz_xxxxz_0, \
                                         tg_xxxxyz_xxxyy_0, tg_xxxxyz_xxxyz_0, tg_xxxxyz_xxxzz_0, tg_xxxxyz_xxyyy_0, \
                                         tg_xxxxyz_xxyyz_0, tg_xxxxyz_xxyzz_0, tg_xxxxyz_xxzzz_0, tg_xxxxyz_xyyyy_0, \
                                         tg_xxxxyz_xyyyz_0, tg_xxxxyz_xyyzz_0, tg_xxxxyz_xyzzz_0, tg_xxxxz_xxxx_1, \
                                         tg_xxxxz_xxxxx_0, tg_xxxxz_xxxxx_1, tg_xxxxz_xxxxy_0, tg_xxxxz_xxxxy_1, \
                                         tg_xxxxz_xxxxz_0, tg_xxxxz_xxxxz_1, tg_xxxxz_xxxy_1, tg_xxxxz_xxxyy_0, \
                                         tg_xxxxz_xxxyy_1, tg_xxxxz_xxxyz_0, tg_xxxxz_xxxyz_1, tg_xxxxz_xxxz_1, \
                                         tg_xxxxz_xxxzz_0, tg_xxxxz_xxxzz_1, tg_xxxxz_xxyy_1, tg_xxxxz_xxyyy_0, \
                                         tg_xxxxz_xxyyy_1, tg_xxxxz_xxyyz_0, tg_xxxxz_xxyyz_1, tg_xxxxz_xxyz_1, \
                                         tg_xxxxz_xxyzz_0, tg_xxxxz_xxyzz_1, tg_xxxxz_xxzz_1, tg_xxxxz_xxzzz_0, \
                                         tg_xxxxz_xxzzz_1, tg_xxxxz_xyyy_1, tg_xxxxz_xyyyy_0, tg_xxxxz_xyyyy_1, \
                                         tg_xxxxz_xyyyz_0, tg_xxxxz_xyyyz_1, tg_xxxxz_xyyz_1, tg_xxxxz_xyyzz_0, \
                                         tg_xxxxz_xyyzz_1, tg_xxxxz_xyzz_1, tg_xxxxz_xyzzz_0, tg_xxxxz_xyzzz_1, \
                                         tg_xxxxz_xzzz_1, tg_xxxxz_xzzzz_0, tg_xxxxz_xzzzz_1, tg_xxxxz_yyyy_1, \
                                         tg_xxxxz_yyyyy_0, tg_xxxxz_yyyyy_1, tg_xxxxz_yyyyz_0, tg_xxxxz_yyyyz_1, \
                                         tg_xxxxz_yyyz_1, tg_xxxxz_yyyzz_0, tg_xxxxz_yyyzz_1, tg_xxxxz_yyzz_1, \
                                         tg_xxxxz_yyzzz_0, tg_xxxxz_yyzzz_1, tg_xxxxz_yzzz_1, tg_xxxxz_yzzzz_0, \
                                         tg_xxxxz_yzzzz_1, tg_xxxxz_zzzz_1, tg_xxxxz_zzzzz_0, tg_xxxxz_zzzzz_1, \
                                         tg_xxxy_xxxxx_0, tg_xxxy_xxxxx_1, tg_xxxy_xxxxy_0, tg_xxxy_xxxxy_1, tg_xxxy_xxxxz_0, \
                                         tg_xxxy_xxxxz_1, tg_xxxy_xxxyy_0, tg_xxxy_xxxyy_1, tg_xxxy_xxxyz_0, tg_xxxy_xxxyz_1, \
                                         tg_xxxy_xxxzz_0, tg_xxxy_xxxzz_1, tg_xxxy_xxyyy_0, tg_xxxy_xxyyy_1, tg_xxxy_xxyyz_0, \
                                         tg_xxxy_xxyyz_1, tg_xxxy_xxyzz_0, tg_xxxy_xxyzz_1, tg_xxxy_xxzzz_0, tg_xxxy_xxzzz_1, \
                                         tg_xxxy_xyyyy_0, tg_xxxy_xyyyy_1, tg_xxxy_xyyyz_0, tg_xxxy_xyyyz_1, tg_xxxy_xyyzz_0, \
                                         tg_xxxy_xyyzz_1, tg_xxxy_xyzzz_0, tg_xxxy_xyzzz_1, tg_xxxy_xzzzz_0, tg_xxxy_xzzzz_1, \
                                         tg_xxxy_yyyyy_0, tg_xxxy_yyyyy_1, tg_xxxy_yyyyz_0, tg_xxxy_yyyyz_1, tg_xxxy_yyyzz_0, \
                                         tg_xxxy_yyyzz_1, tg_xxxy_yyzzz_0, tg_xxxy_yyzzz_1, tg_xxxy_yzzzz_0, tg_xxxy_yzzzz_1, \
                                         tg_xxxy_zzzzz_0, tg_xxxy_zzzzz_1, tg_xxxyy_xxxx_1, tg_xxxyy_xxxxx_0, \
                                         tg_xxxyy_xxxxx_1, tg_xxxyy_xxxxy_0, tg_xxxyy_xxxxy_1, tg_xxxyy_xxxxz_0, \
                                         tg_xxxyy_xxxxz_1, tg_xxxyy_xxxy_1, tg_xxxyy_xxxyy_0, tg_xxxyy_xxxyy_1, \
                                         tg_xxxyy_xxxyz_0, tg_xxxyy_xxxyz_1, tg_xxxyy_xxxz_1, tg_xxxyy_xxxzz_0, \
                                         tg_xxxyy_xxxzz_1, tg_xxxyy_xxyy_1, tg_xxxyy_xxyyy_0, tg_xxxyy_xxyyy_1, \
                                         tg_xxxyy_xxyyz_0, tg_xxxyy_xxyyz_1, tg_xxxyy_xxyz_1, tg_xxxyy_xxyzz_0, \
                                         tg_xxxyy_xxyzz_1, tg_xxxyy_xxzz_1, tg_xxxyy_xxzzz_0, tg_xxxyy_xxzzz_1, \
                                         tg_xxxyy_xyyy_1, tg_xxxyy_xyyyy_0, tg_xxxyy_xyyyy_1, tg_xxxyy_xyyyz_0, \
                                         tg_xxxyy_xyyyz_1, tg_xxxyy_xyyz_1, tg_xxxyy_xyyzz_0, tg_xxxyy_xyyzz_1, \
                                         tg_xxxyy_xyzz_1, tg_xxxyy_xyzzz_0, tg_xxxyy_xyzzz_1, tg_xxxyy_xzzz_1, \
                                         tg_xxxyy_xzzzz_0, tg_xxxyy_xzzzz_1, tg_xxxyy_yyyy_1, tg_xxxyy_yyyyy_0, \
                                         tg_xxxyy_yyyyy_1, tg_xxxyy_yyyyz_0, tg_xxxyy_yyyyz_1, tg_xxxyy_yyyz_1, \
                                         tg_xxxyy_yyyzz_0, tg_xxxyy_yyyzz_1, tg_xxxyy_yyzz_1, tg_xxxyy_yyzzz_0, \
                                         tg_xxxyy_yyzzz_1, tg_xxxyy_yzzz_1, tg_xxxyy_yzzzz_0, tg_xxxyy_yzzzz_1, \
                                         tg_xxxyy_zzzz_1, tg_xxxyy_zzzzz_0, tg_xxxyy_zzzzz_1, tg_xxxyz_xxxx_1, \
                                         tg_xxxyz_xxxxx_0, tg_xxxyz_xxxxx_1, tg_xxxyz_xxxxy_0, tg_xxxyz_xxxxy_1, \
                                         tg_xxxyz_xxxxz_0, tg_xxxyz_xxxxz_1, tg_xxxyz_xxxy_1, tg_xxxyz_xxxyy_0, \
                                         tg_xxxyz_xxxyy_1, tg_xxxyz_xxxyz_0, tg_xxxyz_xxxyz_1, tg_xxxyz_xxxz_1, \
                                         tg_xxxyz_xxxzz_0, tg_xxxyz_xxxzz_1, tg_xxxyz_xxyy_1, tg_xxxyz_xxyyy_0, \
                                         tg_xxxyz_xxyyy_1, tg_xxxyz_xxyyz_0, tg_xxxyz_xxyyz_1, tg_xxxyz_xxyz_1, \
                                         tg_xxxyz_xxyzz_0, tg_xxxyz_xxyzz_1, tg_xxxyz_xxzz_1, tg_xxxyz_xxzzz_0, \
                                         tg_xxxyz_xxzzz_1, tg_xxxyz_xyyy_1, tg_xxxyz_xyyyy_0, tg_xxxyz_xyyyy_1, \
                                         tg_xxxyz_xyyyz_0, tg_xxxyz_xyyyz_1, tg_xxxyz_xyyz_1, tg_xxxyz_xyyzz_0, \
                                         tg_xxxyz_xyyzz_1, tg_xxxyz_xyzz_1, tg_xxxyz_xyzzz_0, tg_xxxyz_xyzzz_1, \
                                         tg_xxxyz_xzzz_1, tg_xxxyz_yyyy_1, tg_xxxyz_yyyz_1, tg_xxxyz_yyzz_1, tg_xxxyz_yzzz_1, \
                                         tg_xxxz_xxxxx_0, tg_xxxz_xxxxx_1, tg_xxxz_xxxxy_0, tg_xxxz_xxxxy_1, tg_xxxz_xxxxz_0, \
                                         tg_xxxz_xxxxz_1, tg_xxxz_xxxyy_0, tg_xxxz_xxxyy_1, tg_xxxz_xxxyz_0, tg_xxxz_xxxyz_1, \
                                         tg_xxxz_xxxzz_0, tg_xxxz_xxxzz_1, tg_xxxz_xxyyy_0, tg_xxxz_xxyyy_1, tg_xxxz_xxyyz_0, \
                                         tg_xxxz_xxyyz_1, tg_xxxz_xxyzz_0, tg_xxxz_xxyzz_1, tg_xxxz_xxzzz_0, tg_xxxz_xxzzz_1, \
                                         tg_xxxz_xyyyy_0, tg_xxxz_xyyyy_1, tg_xxxz_xyyyz_0, tg_xxxz_xyyyz_1, tg_xxxz_xyyzz_0, \
                                         tg_xxxz_xyyzz_1, tg_xxxz_xyzzz_0, tg_xxxz_xyzzz_1, tg_xxxz_xzzzz_0, tg_xxxz_xzzzz_1, \
                                         tg_xxxz_yyyyy_0, tg_xxxz_yyyyy_1, tg_xxxz_yyyyz_0, tg_xxxz_yyyyz_1, tg_xxxz_yyyzz_0, \
                                         tg_xxxz_yyyzz_1, tg_xxxz_yyzzz_0, tg_xxxz_yyzzz_1, tg_xxxz_yzzzz_0, tg_xxxz_yzzzz_1, \
                                         tg_xxxz_zzzzz_0, tg_xxxz_zzzzz_1, tg_xxyy_xxxxx_0, tg_xxyy_xxxxx_1, tg_xxyy_xxxxy_0, \
                                         tg_xxyy_xxxxy_1, tg_xxyy_xxxxz_0, tg_xxyy_xxxxz_1, tg_xxyy_xxxyy_0, tg_xxyy_xxxyy_1, \
                                         tg_xxyy_xxxyz_0, tg_xxyy_xxxyz_1, tg_xxyy_xxxzz_0, tg_xxyy_xxxzz_1, tg_xxyy_xxyyy_0, \
                                         tg_xxyy_xxyyy_1, tg_xxyy_xxyyz_0, tg_xxyy_xxyyz_1, tg_xxyy_xxyzz_0, tg_xxyy_xxyzz_1, \
                                         tg_xxyy_xxzzz_0, tg_xxyy_xxzzz_1, tg_xxyy_xyyyy_0, tg_xxyy_xyyyy_1, tg_xxyy_xyyyz_0, \
                                         tg_xxyy_xyyyz_1, tg_xxyy_xyyzz_0, tg_xxyy_xyyzz_1, tg_xxyy_xyzzz_0, tg_xxyy_xyzzz_1, \
                                         tg_xxyy_xzzzz_0, tg_xxyy_xzzzz_1, tg_xxyy_yyyyy_0, tg_xxyy_yyyyy_1, tg_xxyy_yyyyz_0, \
                                         tg_xxyy_yyyyz_1, tg_xxyy_yyyzz_0, tg_xxyy_yyyzz_1, tg_xxyy_yyzzz_0, tg_xxyy_yyzzz_1, \
                                         tg_xxyy_yzzzz_0, tg_xxyy_yzzzz_1, tg_xxyy_zzzzz_0, tg_xxyy_zzzzz_1, tg_xxyz_xxxxx_0, \
                                         tg_xxyz_xxxxx_1, tg_xxyz_xxxxy_0, tg_xxyz_xxxxy_1, tg_xxyz_xxxxz_0, tg_xxyz_xxxxz_1, \
                                         tg_xxyz_xxxyy_0, tg_xxyz_xxxyy_1, tg_xxyz_xxxyz_0, tg_xxyz_xxxyz_1, tg_xxyz_xxxzz_0, \
                                         tg_xxyz_xxxzz_1, tg_xxyz_xxyyy_0, tg_xxyz_xxyyy_1, tg_xxyz_xxyyz_0, tg_xxyz_xxyyz_1, \
                                         tg_xxyz_xxyzz_0, tg_xxyz_xxyzz_1, tg_xxyz_xxzzz_0, tg_xxyz_xxzzz_1, tg_xxyz_xyyyy_0, \
                                         tg_xxyz_xyyyy_1, tg_xxyz_xyyyz_0, tg_xxyz_xyyyz_1, tg_xxyz_xyyzz_0, tg_xxyz_xyyzz_1, \
                                         tg_xxyz_xyzzz_0, tg_xxyz_xyzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxx_xxxxx_0[j] = pb_x * tg_xxxxx_xxxxx_0[j] + fr * tg_xxxxx_xxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxx_0[j] - tg_xxxx_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxx_xxxx_1[j];

                    tg_xxxxxx_xxxxy_0[j] = pb_x * tg_xxxxx_xxxxy_0[j] + fr * tg_xxxxx_xxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxy_0[j] - tg_xxxx_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxx_xxxy_1[j];

                    tg_xxxxxx_xxxxz_0[j] = pb_x * tg_xxxxx_xxxxz_0[j] + fr * tg_xxxxx_xxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxxz_0[j] - tg_xxxx_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxx_xxxz_1[j];

                    tg_xxxxxx_xxxyy_0[j] = pb_x * tg_xxxxx_xxxyy_0[j] + fr * tg_xxxxx_xxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxyy_0[j] - tg_xxxx_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxx_xxyy_1[j];

                    tg_xxxxxx_xxxyz_0[j] = pb_x * tg_xxxxx_xxxyz_0[j] + fr * tg_xxxxx_xxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxyz_0[j] - tg_xxxx_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxx_xxyz_1[j];

                    tg_xxxxxx_xxxzz_0[j] = pb_x * tg_xxxxx_xxxzz_0[j] + fr * tg_xxxxx_xxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxxzz_0[j] - tg_xxxx_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxx_xxzz_1[j];

                    tg_xxxxxx_xxyyy_0[j] = pb_x * tg_xxxxx_xxyyy_0[j] + fr * tg_xxxxx_xxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxyyy_0[j] - tg_xxxx_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xyyy_1[j];

                    tg_xxxxxx_xxyyz_0[j] = pb_x * tg_xxxxx_xxyyz_0[j] + fr * tg_xxxxx_xxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxyyz_0[j] - tg_xxxx_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xyyz_1[j];

                    tg_xxxxxx_xxyzz_0[j] = pb_x * tg_xxxxx_xxyzz_0[j] + fr * tg_xxxxx_xxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxyzz_0[j] - tg_xxxx_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xyzz_1[j];

                    tg_xxxxxx_xxzzz_0[j] = pb_x * tg_xxxxx_xxzzz_0[j] + fr * tg_xxxxx_xxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xxzzz_0[j] - tg_xxxx_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxx_xzzz_1[j];

                    tg_xxxxxx_xyyyy_0[j] = pb_x * tg_xxxxx_xyyyy_0[j] + fr * tg_xxxxx_xyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyyyy_0[j] - tg_xxxx_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yyyy_1[j];

                    tg_xxxxxx_xyyyz_0[j] = pb_x * tg_xxxxx_xyyyz_0[j] + fr * tg_xxxxx_xyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyyyz_0[j] - tg_xxxx_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yyyz_1[j];

                    tg_xxxxxx_xyyzz_0[j] = pb_x * tg_xxxxx_xyyzz_0[j] + fr * tg_xxxxx_xyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyyzz_0[j] - tg_xxxx_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yyzz_1[j];

                    tg_xxxxxx_xyzzz_0[j] = pb_x * tg_xxxxx_xyzzz_0[j] + fr * tg_xxxxx_xyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xyzzz_0[j] - tg_xxxx_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_yzzz_1[j];

                    tg_xxxxxx_xzzzz_0[j] = pb_x * tg_xxxxx_xzzzz_0[j] + fr * tg_xxxxx_xzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_xzzzz_0[j] - tg_xxxx_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxx_zzzz_1[j];

                    tg_xxxxxx_yyyyy_0[j] = pb_x * tg_xxxxx_yyyyy_0[j] + fr * tg_xxxxx_yyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyyyy_0[j] - tg_xxxx_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxx_yyyyz_0[j] = pb_x * tg_xxxxx_yyyyz_0[j] + fr * tg_xxxxx_yyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyyyz_0[j] - tg_xxxx_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxx_yyyzz_0[j] = pb_x * tg_xxxxx_yyyzz_0[j] + fr * tg_xxxxx_yyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyyzz_0[j] - tg_xxxx_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxx_yyzzz_0[j] = pb_x * tg_xxxxx_yyzzz_0[j] + fr * tg_xxxxx_yyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yyzzz_0[j] - tg_xxxx_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxx_yzzzz_0[j] = pb_x * tg_xxxxx_yzzzz_0[j] + fr * tg_xxxxx_yzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_yzzzz_0[j] - tg_xxxx_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxx_zzzzz_0[j] = pb_x * tg_xxxxx_zzzzz_0[j] + fr * tg_xxxxx_zzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxx_zzzzz_0[j] - tg_xxxx_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxy_xxxxx_0[j] = pb_x * tg_xxxxy_xxxxx_0[j] + fr * tg_xxxxy_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxx_0[j] - tg_xxxy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxy_xxxx_1[j];

                    tg_xxxxxy_xxxxy_0[j] = pb_x * tg_xxxxy_xxxxy_0[j] + fr * tg_xxxxy_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxy_0[j] - tg_xxxy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxy_xxxy_1[j];

                    tg_xxxxxy_xxxxz_0[j] = pb_x * tg_xxxxy_xxxxz_0[j] + fr * tg_xxxxy_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxxz_0[j] - tg_xxxy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxy_xxxz_1[j];

                    tg_xxxxxy_xxxyy_0[j] = pb_x * tg_xxxxy_xxxyy_0[j] + fr * tg_xxxxy_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxyy_0[j] - tg_xxxy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxy_xxyy_1[j];

                    tg_xxxxxy_xxxyz_0[j] = pb_x * tg_xxxxy_xxxyz_0[j] + fr * tg_xxxxy_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxyz_0[j] - tg_xxxy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxy_xxyz_1[j];

                    tg_xxxxxy_xxxzz_0[j] = pb_x * tg_xxxxy_xxxzz_0[j] + fr * tg_xxxxy_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxxzz_0[j] - tg_xxxy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxy_xxzz_1[j];

                    tg_xxxxxy_xxyyy_0[j] = pb_x * tg_xxxxy_xxyyy_0[j] + fr * tg_xxxxy_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxyyy_0[j] - tg_xxxy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xyyy_1[j];

                    tg_xxxxxy_xxyyz_0[j] = pb_x * tg_xxxxy_xxyyz_0[j] + fr * tg_xxxxy_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxyyz_0[j] - tg_xxxy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xyyz_1[j];

                    tg_xxxxxy_xxyzz_0[j] = pb_x * tg_xxxxy_xxyzz_0[j] + fr * tg_xxxxy_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxyzz_0[j] - tg_xxxy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xyzz_1[j];

                    tg_xxxxxy_xxzzz_0[j] = pb_x * tg_xxxxy_xxzzz_0[j] + fr * tg_xxxxy_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xxzzz_0[j] - tg_xxxy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxy_xzzz_1[j];

                    tg_xxxxxy_xyyyy_0[j] = pb_x * tg_xxxxy_xyyyy_0[j] + fr * tg_xxxxy_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyyyy_0[j] - tg_xxxy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yyyy_1[j];

                    tg_xxxxxy_xyyyz_0[j] = pb_x * tg_xxxxy_xyyyz_0[j] + fr * tg_xxxxy_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyyyz_0[j] - tg_xxxy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yyyz_1[j];

                    tg_xxxxxy_xyyzz_0[j] = pb_x * tg_xxxxy_xyyzz_0[j] + fr * tg_xxxxy_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyyzz_0[j] - tg_xxxy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yyzz_1[j];

                    tg_xxxxxy_xyzzz_0[j] = pb_x * tg_xxxxy_xyzzz_0[j] + fr * tg_xxxxy_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xyzzz_0[j] - tg_xxxy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_yzzz_1[j];

                    tg_xxxxxy_xzzzz_0[j] = pb_x * tg_xxxxy_xzzzz_0[j] + fr * tg_xxxxy_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_xzzzz_0[j] - tg_xxxy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxy_zzzz_1[j];

                    tg_xxxxxy_yyyyy_0[j] = pb_x * tg_xxxxy_yyyyy_0[j] + fr * tg_xxxxy_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyyyy_0[j] - tg_xxxy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxy_yyyyz_0[j] = pb_x * tg_xxxxy_yyyyz_0[j] + fr * tg_xxxxy_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyyyz_0[j] - tg_xxxy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxy_yyyzz_0[j] = pb_x * tg_xxxxy_yyyzz_0[j] + fr * tg_xxxxy_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyyzz_0[j] - tg_xxxy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxy_yyzzz_0[j] = pb_x * tg_xxxxy_yyzzz_0[j] + fr * tg_xxxxy_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yyzzz_0[j] - tg_xxxy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxy_yzzzz_0[j] = pb_x * tg_xxxxy_yzzzz_0[j] + fr * tg_xxxxy_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_yzzzz_0[j] - tg_xxxy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxy_zzzzz_0[j] = pb_x * tg_xxxxy_zzzzz_0[j] + fr * tg_xxxxy_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxy_zzzzz_0[j] - tg_xxxy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxz_xxxxx_0[j] = pb_x * tg_xxxxz_xxxxx_0[j] + fr * tg_xxxxz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxx_0[j] - tg_xxxz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxz_xxxx_1[j];

                    tg_xxxxxz_xxxxy_0[j] = pb_x * tg_xxxxz_xxxxy_0[j] + fr * tg_xxxxz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxy_0[j] - tg_xxxz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxz_xxxy_1[j];

                    tg_xxxxxz_xxxxz_0[j] = pb_x * tg_xxxxz_xxxxz_0[j] + fr * tg_xxxxz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxxz_0[j] - tg_xxxz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxz_xxxz_1[j];

                    tg_xxxxxz_xxxyy_0[j] = pb_x * tg_xxxxz_xxxyy_0[j] + fr * tg_xxxxz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxyy_0[j] - tg_xxxz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxz_xxyy_1[j];

                    tg_xxxxxz_xxxyz_0[j] = pb_x * tg_xxxxz_xxxyz_0[j] + fr * tg_xxxxz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxyz_0[j] - tg_xxxz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxz_xxyz_1[j];

                    tg_xxxxxz_xxxzz_0[j] = pb_x * tg_xxxxz_xxxzz_0[j] + fr * tg_xxxxz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxxzz_0[j] - tg_xxxz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxz_xxzz_1[j];

                    tg_xxxxxz_xxyyy_0[j] = pb_x * tg_xxxxz_xxyyy_0[j] + fr * tg_xxxxz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxyyy_0[j] - tg_xxxz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xyyy_1[j];

                    tg_xxxxxz_xxyyz_0[j] = pb_x * tg_xxxxz_xxyyz_0[j] + fr * tg_xxxxz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxyyz_0[j] - tg_xxxz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xyyz_1[j];

                    tg_xxxxxz_xxyzz_0[j] = pb_x * tg_xxxxz_xxyzz_0[j] + fr * tg_xxxxz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxyzz_0[j] - tg_xxxz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xyzz_1[j];

                    tg_xxxxxz_xxzzz_0[j] = pb_x * tg_xxxxz_xxzzz_0[j] + fr * tg_xxxxz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xxzzz_0[j] - tg_xxxz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxz_xzzz_1[j];

                    tg_xxxxxz_xyyyy_0[j] = pb_x * tg_xxxxz_xyyyy_0[j] + fr * tg_xxxxz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyyyy_0[j] - tg_xxxz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yyyy_1[j];

                    tg_xxxxxz_xyyyz_0[j] = pb_x * tg_xxxxz_xyyyz_0[j] + fr * tg_xxxxz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyyyz_0[j] - tg_xxxz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yyyz_1[j];

                    tg_xxxxxz_xyyzz_0[j] = pb_x * tg_xxxxz_xyyzz_0[j] + fr * tg_xxxxz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyyzz_0[j] - tg_xxxz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yyzz_1[j];

                    tg_xxxxxz_xyzzz_0[j] = pb_x * tg_xxxxz_xyzzz_0[j] + fr * tg_xxxxz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xyzzz_0[j] - tg_xxxz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_yzzz_1[j];

                    tg_xxxxxz_xzzzz_0[j] = pb_x * tg_xxxxz_xzzzz_0[j] + fr * tg_xxxxz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_xzzzz_0[j] - tg_xxxz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxz_zzzz_1[j];

                    tg_xxxxxz_yyyyy_0[j] = pb_x * tg_xxxxz_yyyyy_0[j] + fr * tg_xxxxz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyyyy_0[j] - tg_xxxz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxz_yyyyz_0[j] = pb_x * tg_xxxxz_yyyyz_0[j] + fr * tg_xxxxz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyyyz_0[j] - tg_xxxz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxz_yyyzz_0[j] = pb_x * tg_xxxxz_yyyzz_0[j] + fr * tg_xxxxz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyyzz_0[j] - tg_xxxz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxz_yyzzz_0[j] = pb_x * tg_xxxxz_yyzzz_0[j] + fr * tg_xxxxz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yyzzz_0[j] - tg_xxxz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxz_yzzzz_0[j] = pb_x * tg_xxxxz_yzzzz_0[j] + fr * tg_xxxxz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_yzzzz_0[j] - tg_xxxz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxz_zzzzz_0[j] = pb_x * tg_xxxxz_zzzzz_0[j] + fr * tg_xxxxz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxz_zzzzz_0[j] - tg_xxxz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyy_xxxxx_0[j] = pb_x * tg_xxxyy_xxxxx_0[j] + fr * tg_xxxyy_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxx_0[j] - tg_xxyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyy_xxxx_1[j];

                    tg_xxxxyy_xxxxy_0[j] = pb_x * tg_xxxyy_xxxxy_0[j] + fr * tg_xxxyy_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxy_0[j] - tg_xxyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyy_xxxy_1[j];

                    tg_xxxxyy_xxxxz_0[j] = pb_x * tg_xxxyy_xxxxz_0[j] + fr * tg_xxxyy_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxxz_0[j] - tg_xxyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyy_xxxz_1[j];

                    tg_xxxxyy_xxxyy_0[j] = pb_x * tg_xxxyy_xxxyy_0[j] + fr * tg_xxxyy_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxyy_0[j] - tg_xxyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyy_xxyy_1[j];

                    tg_xxxxyy_xxxyz_0[j] = pb_x * tg_xxxyy_xxxyz_0[j] + fr * tg_xxxyy_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxyz_0[j] - tg_xxyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyy_xxyz_1[j];

                    tg_xxxxyy_xxxzz_0[j] = pb_x * tg_xxxyy_xxxzz_0[j] + fr * tg_xxxyy_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxxzz_0[j] - tg_xxyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyy_xxzz_1[j];

                    tg_xxxxyy_xxyyy_0[j] = pb_x * tg_xxxyy_xxyyy_0[j] + fr * tg_xxxyy_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxyyy_0[j] - tg_xxyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xyyy_1[j];

                    tg_xxxxyy_xxyyz_0[j] = pb_x * tg_xxxyy_xxyyz_0[j] + fr * tg_xxxyy_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxyyz_0[j] - tg_xxyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xyyz_1[j];

                    tg_xxxxyy_xxyzz_0[j] = pb_x * tg_xxxyy_xxyzz_0[j] + fr * tg_xxxyy_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxyzz_0[j] - tg_xxyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xyzz_1[j];

                    tg_xxxxyy_xxzzz_0[j] = pb_x * tg_xxxyy_xxzzz_0[j] + fr * tg_xxxyy_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xxzzz_0[j] - tg_xxyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyy_xzzz_1[j];

                    tg_xxxxyy_xyyyy_0[j] = pb_x * tg_xxxyy_xyyyy_0[j] + fr * tg_xxxyy_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyyyy_0[j] - tg_xxyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yyyy_1[j];

                    tg_xxxxyy_xyyyz_0[j] = pb_x * tg_xxxyy_xyyyz_0[j] + fr * tg_xxxyy_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyyyz_0[j] - tg_xxyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yyyz_1[j];

                    tg_xxxxyy_xyyzz_0[j] = pb_x * tg_xxxyy_xyyzz_0[j] + fr * tg_xxxyy_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyyzz_0[j] - tg_xxyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yyzz_1[j];

                    tg_xxxxyy_xyzzz_0[j] = pb_x * tg_xxxyy_xyzzz_0[j] + fr * tg_xxxyy_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xyzzz_0[j] - tg_xxyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_yzzz_1[j];

                    tg_xxxxyy_xzzzz_0[j] = pb_x * tg_xxxyy_xzzzz_0[j] + fr * tg_xxxyy_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_xzzzz_0[j] - tg_xxyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyy_zzzz_1[j];

                    tg_xxxxyy_yyyyy_0[j] = pb_x * tg_xxxyy_yyyyy_0[j] + fr * tg_xxxyy_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyyyy_0[j] - tg_xxyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyy_yyyyz_0[j] = pb_x * tg_xxxyy_yyyyz_0[j] + fr * tg_xxxyy_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyyyz_0[j] - tg_xxyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyy_yyyzz_0[j] = pb_x * tg_xxxyy_yyyzz_0[j] + fr * tg_xxxyy_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyyzz_0[j] - tg_xxyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyy_yyzzz_0[j] = pb_x * tg_xxxyy_yyzzz_0[j] + fr * tg_xxxyy_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yyzzz_0[j] - tg_xxyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyy_yzzzz_0[j] = pb_x * tg_xxxyy_yzzzz_0[j] + fr * tg_xxxyy_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_yzzzz_0[j] - tg_xxyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyy_zzzzz_0[j] = pb_x * tg_xxxyy_zzzzz_0[j] + fr * tg_xxxyy_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyy_zzzzz_0[j] - tg_xxyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyz_xxxxx_0[j] = pb_x * tg_xxxyz_xxxxx_0[j] + fr * tg_xxxyz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxx_0[j] - tg_xxyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyz_xxxx_1[j];

                    tg_xxxxyz_xxxxy_0[j] = pb_x * tg_xxxyz_xxxxy_0[j] + fr * tg_xxxyz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxy_0[j] - tg_xxyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyz_xxxy_1[j];

                    tg_xxxxyz_xxxxz_0[j] = pb_x * tg_xxxyz_xxxxz_0[j] + fr * tg_xxxyz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxxz_0[j] - tg_xxyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyz_xxxz_1[j];

                    tg_xxxxyz_xxxyy_0[j] = pb_x * tg_xxxyz_xxxyy_0[j] + fr * tg_xxxyz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxyy_0[j] - tg_xxyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyz_xxyy_1[j];

                    tg_xxxxyz_xxxyz_0[j] = pb_x * tg_xxxyz_xxxyz_0[j] + fr * tg_xxxyz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxyz_0[j] - tg_xxyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyz_xxyz_1[j];

                    tg_xxxxyz_xxxzz_0[j] = pb_x * tg_xxxyz_xxxzz_0[j] + fr * tg_xxxyz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxxzz_0[j] - tg_xxyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyz_xxzz_1[j];

                    tg_xxxxyz_xxyyy_0[j] = pb_x * tg_xxxyz_xxyyy_0[j] + fr * tg_xxxyz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxyyy_0[j] - tg_xxyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xyyy_1[j];

                    tg_xxxxyz_xxyyz_0[j] = pb_x * tg_xxxyz_xxyyz_0[j] + fr * tg_xxxyz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxyyz_0[j] - tg_xxyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xyyz_1[j];

                    tg_xxxxyz_xxyzz_0[j] = pb_x * tg_xxxyz_xxyzz_0[j] + fr * tg_xxxyz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxyzz_0[j] - tg_xxyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xyzz_1[j];

                    tg_xxxxyz_xxzzz_0[j] = pb_x * tg_xxxyz_xxzzz_0[j] + fr * tg_xxxyz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xxzzz_0[j] - tg_xxyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyz_xzzz_1[j];

                    tg_xxxxyz_xyyyy_0[j] = pb_x * tg_xxxyz_xyyyy_0[j] + fr * tg_xxxyz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyyyy_0[j] - tg_xxyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yyyy_1[j];

                    tg_xxxxyz_xyyyz_0[j] = pb_x * tg_xxxyz_xyyyz_0[j] + fr * tg_xxxyz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyyyz_0[j] - tg_xxyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yyyz_1[j];

                    tg_xxxxyz_xyyzz_0[j] = pb_x * tg_xxxyz_xyyzz_0[j] + fr * tg_xxxyz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyyzz_0[j] - tg_xxyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yyzz_1[j];

                    tg_xxxxyz_xyzzz_0[j] = pb_x * tg_xxxyz_xyzzz_0[j] + fr * tg_xxxyz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xyzzz_0[j] - tg_xxyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_yzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISH_98_196(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxzzz_xxxxy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 190); 

                auto tg_xxzzz_xxxxz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 191); 

                auto tg_xxzzz_xxxyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 192); 

                auto tg_xxzzz_xxxyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 193); 

                auto tg_xxzzz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 194); 

                auto tg_xxzzz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 195); 

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

                auto tg_xxzzz_xxxxy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 190); 

                auto tg_xxzzz_xxxxz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 191); 

                auto tg_xxzzz_xxxyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 192); 

                auto tg_xxzzz_xxxyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 193); 

                auto tg_xxzzz_xxxzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 194); 

                auto tg_xxzzz_xxyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 195); 

                auto tg_xxyz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 98); 

                auto tg_xxyz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 99); 

                auto tg_xxyz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 100); 

                auto tg_xxyz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 101); 

                auto tg_xxyz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 102); 

                auto tg_xxyz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 103); 

                auto tg_xxyz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 104); 

                auto tg_xxzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 105); 

                auto tg_xxzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 106); 

                auto tg_xxzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 107); 

                auto tg_xxzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 108); 

                auto tg_xxzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 109); 

                auto tg_xxzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 110); 

                auto tg_xxzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 111); 

                auto tg_xxzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 112); 

                auto tg_xxzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 113); 

                auto tg_xxzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 114); 

                auto tg_xxzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 115); 

                auto tg_xxzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 116); 

                auto tg_xxzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 117); 

                auto tg_xxzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 118); 

                auto tg_xxzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 119); 

                auto tg_xxzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 120); 

                auto tg_xxzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 121); 

                auto tg_xxzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 122); 

                auto tg_xxzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 123); 

                auto tg_xxzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 124); 

                auto tg_xxzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 125); 

                auto tg_xyyy_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 126); 

                auto tg_xyyy_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 127); 

                auto tg_xyyy_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 128); 

                auto tg_xyyy_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 129); 

                auto tg_xyyy_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 130); 

                auto tg_xyyy_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 131); 

                auto tg_xyyy_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 132); 

                auto tg_xyyy_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 133); 

                auto tg_xyyy_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 134); 

                auto tg_xyyy_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 135); 

                auto tg_xyyy_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 136); 

                auto tg_xyyy_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 137); 

                auto tg_xyyy_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 138); 

                auto tg_xyyy_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 139); 

                auto tg_xyyy_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 140); 

                auto tg_xyyy_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 141); 

                auto tg_xyyy_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 142); 

                auto tg_xyyy_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 143); 

                auto tg_xyyy_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 144); 

                auto tg_xyyy_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 145); 

                auto tg_xyyy_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 146); 

                auto tg_xyyz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 147); 

                auto tg_xyyz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 148); 

                auto tg_xyyz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 149); 

                auto tg_xyyz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 150); 

                auto tg_xyyz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 151); 

                auto tg_xyyz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 152); 

                auto tg_xyyz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 153); 

                auto tg_xyyz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 154); 

                auto tg_xyyz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 155); 

                auto tg_xyyz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 156); 

                auto tg_xyyz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 157); 

                auto tg_xyyz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 158); 

                auto tg_xyyz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 159); 

                auto tg_xyyz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 160); 

                auto tg_xyyz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 161); 

                auto tg_xyyz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 162); 

                auto tg_xyyz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 163); 

                auto tg_xyyz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 164); 

                auto tg_xyyz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 165); 

                auto tg_xyyz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 166); 

                auto tg_xyyz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 167); 

                auto tg_xyzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 168); 

                auto tg_xyzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 169); 

                auto tg_xyzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 170); 

                auto tg_xyzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 171); 

                auto tg_xyzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 172); 

                auto tg_xyzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 173); 

                auto tg_xyzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 174); 

                auto tg_xyzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 175); 

                auto tg_xyzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 176); 

                auto tg_xyzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 177); 

                auto tg_xyzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 178); 

                auto tg_xyzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 179); 

                auto tg_xyzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 180); 

                auto tg_xyzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 181); 

                auto tg_xyzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 182); 

                auto tg_xyzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 183); 

                auto tg_xyzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 184); 

                auto tg_xyzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 185); 

                auto tg_xyzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 186); 

                auto tg_xyzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 187); 

                auto tg_xyzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 188); 

                auto tg_xzzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 189); 

                auto tg_xzzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 190); 

                auto tg_xzzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 191); 

                auto tg_xzzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 192); 

                auto tg_xzzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 193); 

                auto tg_xzzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 194); 

                auto tg_xzzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 195); 

                auto tg_xxyz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 98); 

                auto tg_xxyz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 99); 

                auto tg_xxyz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 100); 

                auto tg_xxyz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 101); 

                auto tg_xxyz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 102); 

                auto tg_xxyz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 103); 

                auto tg_xxyz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 104); 

                auto tg_xxzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 105); 

                auto tg_xxzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 106); 

                auto tg_xxzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 107); 

                auto tg_xxzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 108); 

                auto tg_xxzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 109); 

                auto tg_xxzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 110); 

                auto tg_xxzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 111); 

                auto tg_xxzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 112); 

                auto tg_xxzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 113); 

                auto tg_xxzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 114); 

                auto tg_xxzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 115); 

                auto tg_xxzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 116); 

                auto tg_xxzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 117); 

                auto tg_xxzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 118); 

                auto tg_xxzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 119); 

                auto tg_xxzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 120); 

                auto tg_xxzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 121); 

                auto tg_xxzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 122); 

                auto tg_xxzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 123); 

                auto tg_xxzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 124); 

                auto tg_xxzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 125); 

                auto tg_xyyy_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 126); 

                auto tg_xyyy_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 127); 

                auto tg_xyyy_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 128); 

                auto tg_xyyy_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 129); 

                auto tg_xyyy_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 130); 

                auto tg_xyyy_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 131); 

                auto tg_xyyy_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 132); 

                auto tg_xyyy_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 133); 

                auto tg_xyyy_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 134); 

                auto tg_xyyy_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 135); 

                auto tg_xyyy_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 136); 

                auto tg_xyyy_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 137); 

                auto tg_xyyy_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 138); 

                auto tg_xyyy_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 139); 

                auto tg_xyyy_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 140); 

                auto tg_xyyy_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 141); 

                auto tg_xyyy_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 142); 

                auto tg_xyyy_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 143); 

                auto tg_xyyy_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 144); 

                auto tg_xyyy_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 145); 

                auto tg_xyyy_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 146); 

                auto tg_xyyz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 147); 

                auto tg_xyyz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 148); 

                auto tg_xyyz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 149); 

                auto tg_xyyz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 150); 

                auto tg_xyyz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 151); 

                auto tg_xyyz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 152); 

                auto tg_xyyz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 153); 

                auto tg_xyyz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 154); 

                auto tg_xyyz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 155); 

                auto tg_xyyz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 156); 

                auto tg_xyyz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 157); 

                auto tg_xyyz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 158); 

                auto tg_xyyz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 159); 

                auto tg_xyyz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 160); 

                auto tg_xyyz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 161); 

                auto tg_xyyz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 162); 

                auto tg_xyyz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 163); 

                auto tg_xyyz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 164); 

                auto tg_xyyz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 165); 

                auto tg_xyyz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 166); 

                auto tg_xyyz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 167); 

                auto tg_xyzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 168); 

                auto tg_xyzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 169); 

                auto tg_xyzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 170); 

                auto tg_xyzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 171); 

                auto tg_xyzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 172); 

                auto tg_xyzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 173); 

                auto tg_xyzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 174); 

                auto tg_xyzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 175); 

                auto tg_xyzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 176); 

                auto tg_xyzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 177); 

                auto tg_xyzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 178); 

                auto tg_xyzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 179); 

                auto tg_xyzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 180); 

                auto tg_xyzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 181); 

                auto tg_xyzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 182); 

                auto tg_xyzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 183); 

                auto tg_xyzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 184); 

                auto tg_xyzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 185); 

                auto tg_xyzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 186); 

                auto tg_xyzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 187); 

                auto tg_xyzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 188); 

                auto tg_xzzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 189); 

                auto tg_xzzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 190); 

                auto tg_xzzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 191); 

                auto tg_xzzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 192); 

                auto tg_xzzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 193); 

                auto tg_xzzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 194); 

                auto tg_xzzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 195); 

                auto tg_xxxyz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 74); 

                auto tg_xxxzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 75); 

                auto tg_xxxzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 76); 

                auto tg_xxxzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 77); 

                auto tg_xxxzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 78); 

                auto tg_xxxzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 79); 

                auto tg_xxxzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 80); 

                auto tg_xxxzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 81); 

                auto tg_xxxzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 82); 

                auto tg_xxxzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 83); 

                auto tg_xxxzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 84); 

                auto tg_xxxzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 85); 

                auto tg_xxxzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 86); 

                auto tg_xxxzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 87); 

                auto tg_xxxzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 88); 

                auto tg_xxxzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 89); 

                auto tg_xxyyy_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 90); 

                auto tg_xxyyy_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 91); 

                auto tg_xxyyy_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 92); 

                auto tg_xxyyy_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 93); 

                auto tg_xxyyy_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 94); 

                auto tg_xxyyy_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 95); 

                auto tg_xxyyy_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 96); 

                auto tg_xxyyy_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 97); 

                auto tg_xxyyy_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 98); 

                auto tg_xxyyy_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 99); 

                auto tg_xxyyy_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 100); 

                auto tg_xxyyy_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 101); 

                auto tg_xxyyy_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 102); 

                auto tg_xxyyy_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 103); 

                auto tg_xxyyy_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 104); 

                auto tg_xxyyz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 105); 

                auto tg_xxyyz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 106); 

                auto tg_xxyyz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 107); 

                auto tg_xxyyz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 108); 

                auto tg_xxyyz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 109); 

                auto tg_xxyyz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 110); 

                auto tg_xxyyz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 111); 

                auto tg_xxyyz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 112); 

                auto tg_xxyyz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 113); 

                auto tg_xxyyz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 114); 

                auto tg_xxyyz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 115); 

                auto tg_xxyyz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 116); 

                auto tg_xxyyz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 117); 

                auto tg_xxyyz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 118); 

                auto tg_xxyyz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 119); 

                auto tg_xxyzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 120); 

                auto tg_xxyzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 121); 

                auto tg_xxyzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 122); 

                auto tg_xxyzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 123); 

                auto tg_xxyzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 124); 

                auto tg_xxyzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 125); 

                auto tg_xxyzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 126); 

                auto tg_xxyzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 127); 

                auto tg_xxyzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 128); 

                auto tg_xxyzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 129); 

                auto tg_xxyzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 130); 

                auto tg_xxyzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 131); 

                auto tg_xxyzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 132); 

                auto tg_xxyzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 133); 

                auto tg_xxyzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 134); 

                auto tg_xxzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 135); 

                auto tg_xxzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 136); 

                auto tg_xxzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 137); 

                auto tg_xxzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 138); 

                auto tg_xxzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 139); 

                auto tg_xxzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 140); 

                auto tg_xxzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 141); 

                // set up pointers to integrals

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

                auto tg_xxxzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 190); 

                auto tg_xxxzzz_xxxxz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 191); 

                auto tg_xxxzzz_xxxyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 192); 

                auto tg_xxxzzz_xxxyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 193); 

                auto tg_xxxzzz_xxxzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 194); 

                auto tg_xxxzzz_xxyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 195); 

                // Batch of Integrals (98,196)

                #pragma omp simd aligned(fxn, fza, tg_xxxxyz_xzzzz_0, tg_xxxxyz_yyyyy_0, tg_xxxxyz_yyyyz_0, \
                                         tg_xxxxyz_yyyzz_0, tg_xxxxyz_yyzzz_0, tg_xxxxyz_yzzzz_0, tg_xxxxyz_zzzzz_0, \
                                         tg_xxxxzz_xxxxx_0, tg_xxxxzz_xxxxy_0, tg_xxxxzz_xxxxz_0, tg_xxxxzz_xxxyy_0, \
                                         tg_xxxxzz_xxxyz_0, tg_xxxxzz_xxxzz_0, tg_xxxxzz_xxyyy_0, tg_xxxxzz_xxyyz_0, \
                                         tg_xxxxzz_xxyzz_0, tg_xxxxzz_xxzzz_0, tg_xxxxzz_xyyyy_0, tg_xxxxzz_xyyyz_0, \
                                         tg_xxxxzz_xyyzz_0, tg_xxxxzz_xyzzz_0, tg_xxxxzz_xzzzz_0, tg_xxxxzz_yyyyy_0, \
                                         tg_xxxxzz_yyyyz_0, tg_xxxxzz_yyyzz_0, tg_xxxxzz_yyzzz_0, tg_xxxxzz_yzzzz_0, \
                                         tg_xxxxzz_zzzzz_0, tg_xxxyyy_xxxxx_0, tg_xxxyyy_xxxxy_0, tg_xxxyyy_xxxxz_0, \
                                         tg_xxxyyy_xxxyy_0, tg_xxxyyy_xxxyz_0, tg_xxxyyy_xxxzz_0, tg_xxxyyy_xxyyy_0, \
                                         tg_xxxyyy_xxyyz_0, tg_xxxyyy_xxyzz_0, tg_xxxyyy_xxzzz_0, tg_xxxyyy_xyyyy_0, \
                                         tg_xxxyyy_xyyyz_0, tg_xxxyyy_xyyzz_0, tg_xxxyyy_xyzzz_0, tg_xxxyyy_xzzzz_0, \
                                         tg_xxxyyy_yyyyy_0, tg_xxxyyy_yyyyz_0, tg_xxxyyy_yyyzz_0, tg_xxxyyy_yyzzz_0, \
                                         tg_xxxyyy_yzzzz_0, tg_xxxyyy_zzzzz_0, tg_xxxyyz_xxxxx_0, tg_xxxyyz_xxxxy_0, \
                                         tg_xxxyyz_xxxxz_0, tg_xxxyyz_xxxyy_0, tg_xxxyyz_xxxyz_0, tg_xxxyyz_xxxzz_0, \
                                         tg_xxxyyz_xxyyy_0, tg_xxxyyz_xxyyz_0, tg_xxxyyz_xxyzz_0, tg_xxxyyz_xxzzz_0, \
                                         tg_xxxyyz_xyyyy_0, tg_xxxyyz_xyyyz_0, tg_xxxyyz_xyyzz_0, tg_xxxyyz_xyzzz_0, \
                                         tg_xxxyyz_xzzzz_0, tg_xxxyyz_yyyyy_0, tg_xxxyyz_yyyyz_0, tg_xxxyyz_yyyzz_0, \
                                         tg_xxxyyz_yyzzz_0, tg_xxxyyz_yzzzz_0, tg_xxxyyz_zzzzz_0, tg_xxxyz_xzzzz_0, \
                                         tg_xxxyz_xzzzz_1, tg_xxxyz_yyyyy_0, tg_xxxyz_yyyyy_1, tg_xxxyz_yyyyz_0, \
                                         tg_xxxyz_yyyyz_1, tg_xxxyz_yyyzz_0, tg_xxxyz_yyyzz_1, tg_xxxyz_yyzzz_0, \
                                         tg_xxxyz_yyzzz_1, tg_xxxyz_yzzzz_0, tg_xxxyz_yzzzz_1, tg_xxxyz_zzzz_1, \
                                         tg_xxxyz_zzzzz_0, tg_xxxyz_zzzzz_1, tg_xxxyzz_xxxxx_0, tg_xxxyzz_xxxxy_0, \
                                         tg_xxxyzz_xxxxz_0, tg_xxxyzz_xxxyy_0, tg_xxxyzz_xxxyz_0, tg_xxxyzz_xxxzz_0, \
                                         tg_xxxyzz_xxyyy_0, tg_xxxyzz_xxyyz_0, tg_xxxyzz_xxyzz_0, tg_xxxyzz_xxzzz_0, \
                                         tg_xxxyzz_xyyyy_0, tg_xxxyzz_xyyyz_0, tg_xxxyzz_xyyzz_0, tg_xxxyzz_xyzzz_0, \
                                         tg_xxxyzz_xzzzz_0, tg_xxxyzz_yyyyy_0, tg_xxxyzz_yyyyz_0, tg_xxxyzz_yyyzz_0, \
                                         tg_xxxyzz_yyzzz_0, tg_xxxyzz_yzzzz_0, tg_xxxyzz_zzzzz_0, tg_xxxzz_xxxx_1, \
                                         tg_xxxzz_xxxxx_0, tg_xxxzz_xxxxx_1, tg_xxxzz_xxxxy_0, tg_xxxzz_xxxxy_1, \
                                         tg_xxxzz_xxxxz_0, tg_xxxzz_xxxxz_1, tg_xxxzz_xxxy_1, tg_xxxzz_xxxyy_0, \
                                         tg_xxxzz_xxxyy_1, tg_xxxzz_xxxyz_0, tg_xxxzz_xxxyz_1, tg_xxxzz_xxxz_1, \
                                         tg_xxxzz_xxxzz_0, tg_xxxzz_xxxzz_1, tg_xxxzz_xxyy_1, tg_xxxzz_xxyyy_0, \
                                         tg_xxxzz_xxyyy_1, tg_xxxzz_xxyyz_0, tg_xxxzz_xxyyz_1, tg_xxxzz_xxyz_1, \
                                         tg_xxxzz_xxyzz_0, tg_xxxzz_xxyzz_1, tg_xxxzz_xxzz_1, tg_xxxzz_xxzzz_0, \
                                         tg_xxxzz_xxzzz_1, tg_xxxzz_xyyy_1, tg_xxxzz_xyyyy_0, tg_xxxzz_xyyyy_1, \
                                         tg_xxxzz_xyyyz_0, tg_xxxzz_xyyyz_1, tg_xxxzz_xyyz_1, tg_xxxzz_xyyzz_0, \
                                         tg_xxxzz_xyyzz_1, tg_xxxzz_xyzz_1, tg_xxxzz_xyzzz_0, tg_xxxzz_xyzzz_1, \
                                         tg_xxxzz_xzzz_1, tg_xxxzz_xzzzz_0, tg_xxxzz_xzzzz_1, tg_xxxzz_yyyy_1, \
                                         tg_xxxzz_yyyyy_0, tg_xxxzz_yyyyy_1, tg_xxxzz_yyyyz_0, tg_xxxzz_yyyyz_1, \
                                         tg_xxxzz_yyyz_1, tg_xxxzz_yyyzz_0, tg_xxxzz_yyyzz_1, tg_xxxzz_yyzz_1, \
                                         tg_xxxzz_yyzzz_0, tg_xxxzz_yyzzz_1, tg_xxxzz_yzzz_1, tg_xxxzz_yzzzz_0, \
                                         tg_xxxzz_yzzzz_1, tg_xxxzz_zzzz_1, tg_xxxzz_zzzzz_0, tg_xxxzz_zzzzz_1, \
                                         tg_xxxzzz_xxxxx_0, tg_xxxzzz_xxxxy_0, tg_xxxzzz_xxxxz_0, tg_xxxzzz_xxxyy_0, \
                                         tg_xxxzzz_xxxyz_0, tg_xxxzzz_xxxzz_0, tg_xxxzzz_xxyyy_0, tg_xxyyy_xxxx_1, \
                                         tg_xxyyy_xxxxx_0, tg_xxyyy_xxxxx_1, tg_xxyyy_xxxxy_0, tg_xxyyy_xxxxy_1, \
                                         tg_xxyyy_xxxxz_0, tg_xxyyy_xxxxz_1, tg_xxyyy_xxxy_1, tg_xxyyy_xxxyy_0, \
                                         tg_xxyyy_xxxyy_1, tg_xxyyy_xxxyz_0, tg_xxyyy_xxxyz_1, tg_xxyyy_xxxz_1, \
                                         tg_xxyyy_xxxzz_0, tg_xxyyy_xxxzz_1, tg_xxyyy_xxyy_1, tg_xxyyy_xxyyy_0, \
                                         tg_xxyyy_xxyyy_1, tg_xxyyy_xxyyz_0, tg_xxyyy_xxyyz_1, tg_xxyyy_xxyz_1, \
                                         tg_xxyyy_xxyzz_0, tg_xxyyy_xxyzz_1, tg_xxyyy_xxzz_1, tg_xxyyy_xxzzz_0, \
                                         tg_xxyyy_xxzzz_1, tg_xxyyy_xyyy_1, tg_xxyyy_xyyyy_0, tg_xxyyy_xyyyy_1, \
                                         tg_xxyyy_xyyyz_0, tg_xxyyy_xyyyz_1, tg_xxyyy_xyyz_1, tg_xxyyy_xyyzz_0, \
                                         tg_xxyyy_xyyzz_1, tg_xxyyy_xyzz_1, tg_xxyyy_xyzzz_0, tg_xxyyy_xyzzz_1, \
                                         tg_xxyyy_xzzz_1, tg_xxyyy_xzzzz_0, tg_xxyyy_xzzzz_1, tg_xxyyy_yyyy_1, \
                                         tg_xxyyy_yyyyy_0, tg_xxyyy_yyyyy_1, tg_xxyyy_yyyyz_0, tg_xxyyy_yyyyz_1, \
                                         tg_xxyyy_yyyz_1, tg_xxyyy_yyyzz_0, tg_xxyyy_yyyzz_1, tg_xxyyy_yyzz_1, \
                                         tg_xxyyy_yyzzz_0, tg_xxyyy_yyzzz_1, tg_xxyyy_yzzz_1, tg_xxyyy_yzzzz_0, \
                                         tg_xxyyy_yzzzz_1, tg_xxyyy_zzzz_1, tg_xxyyy_zzzzz_0, tg_xxyyy_zzzzz_1, \
                                         tg_xxyyz_xxxx_1, tg_xxyyz_xxxxx_0, tg_xxyyz_xxxxx_1, tg_xxyyz_xxxxy_0, \
                                         tg_xxyyz_xxxxy_1, tg_xxyyz_xxxxz_0, tg_xxyyz_xxxxz_1, tg_xxyyz_xxxy_1, \
                                         tg_xxyyz_xxxyy_0, tg_xxyyz_xxxyy_1, tg_xxyyz_xxxyz_0, tg_xxyyz_xxxyz_1, \
                                         tg_xxyyz_xxxz_1, tg_xxyyz_xxxzz_0, tg_xxyyz_xxxzz_1, tg_xxyyz_xxyy_1, \
                                         tg_xxyyz_xxyyy_0, tg_xxyyz_xxyyy_1, tg_xxyyz_xxyyz_0, tg_xxyyz_xxyyz_1, \
                                         tg_xxyyz_xxyz_1, tg_xxyyz_xxyzz_0, tg_xxyyz_xxyzz_1, tg_xxyyz_xxzz_1, \
                                         tg_xxyyz_xxzzz_0, tg_xxyyz_xxzzz_1, tg_xxyyz_xyyy_1, tg_xxyyz_xyyyy_0, \
                                         tg_xxyyz_xyyyy_1, tg_xxyyz_xyyyz_0, tg_xxyyz_xyyyz_1, tg_xxyyz_xyyz_1, \
                                         tg_xxyyz_xyyzz_0, tg_xxyyz_xyyzz_1, tg_xxyyz_xyzz_1, tg_xxyyz_xyzzz_0, \
                                         tg_xxyyz_xyzzz_1, tg_xxyyz_xzzz_1, tg_xxyyz_xzzzz_0, tg_xxyyz_xzzzz_1, \
                                         tg_xxyyz_yyyy_1, tg_xxyyz_yyyyy_0, tg_xxyyz_yyyyy_1, tg_xxyyz_yyyyz_0, \
                                         tg_xxyyz_yyyyz_1, tg_xxyyz_yyyz_1, tg_xxyyz_yyyzz_0, tg_xxyyz_yyyzz_1, \
                                         tg_xxyyz_yyzz_1, tg_xxyyz_yyzzz_0, tg_xxyyz_yyzzz_1, tg_xxyyz_yzzz_1, \
                                         tg_xxyyz_yzzzz_0, tg_xxyyz_yzzzz_1, tg_xxyyz_zzzz_1, tg_xxyyz_zzzzz_0, \
                                         tg_xxyyz_zzzzz_1, tg_xxyz_xzzzz_0, tg_xxyz_xzzzz_1, tg_xxyz_yyyyy_0, tg_xxyz_yyyyy_1, \
                                         tg_xxyz_yyyyz_0, tg_xxyz_yyyyz_1, tg_xxyz_yyyzz_0, tg_xxyz_yyyzz_1, tg_xxyz_yyzzz_0, \
                                         tg_xxyz_yyzzz_1, tg_xxyz_yzzzz_0, tg_xxyz_yzzzz_1, tg_xxyz_zzzzz_0, tg_xxyz_zzzzz_1, \
                                         tg_xxyzz_xxxx_1, tg_xxyzz_xxxxx_0, tg_xxyzz_xxxxx_1, tg_xxyzz_xxxxy_0, \
                                         tg_xxyzz_xxxxy_1, tg_xxyzz_xxxxz_0, tg_xxyzz_xxxxz_1, tg_xxyzz_xxxy_1, \
                                         tg_xxyzz_xxxyy_0, tg_xxyzz_xxxyy_1, tg_xxyzz_xxxyz_0, tg_xxyzz_xxxyz_1, \
                                         tg_xxyzz_xxxz_1, tg_xxyzz_xxxzz_0, tg_xxyzz_xxxzz_1, tg_xxyzz_xxyy_1, \
                                         tg_xxyzz_xxyyy_0, tg_xxyzz_xxyyy_1, tg_xxyzz_xxyyz_0, tg_xxyzz_xxyyz_1, \
                                         tg_xxyzz_xxyz_1, tg_xxyzz_xxyzz_0, tg_xxyzz_xxyzz_1, tg_xxyzz_xxzz_1, \
                                         tg_xxyzz_xxzzz_0, tg_xxyzz_xxzzz_1, tg_xxyzz_xyyy_1, tg_xxyzz_xyyyy_0, \
                                         tg_xxyzz_xyyyy_1, tg_xxyzz_xyyyz_0, tg_xxyzz_xyyyz_1, tg_xxyzz_xyyz_1, \
                                         tg_xxyzz_xyyzz_0, tg_xxyzz_xyyzz_1, tg_xxyzz_xyzz_1, tg_xxyzz_xyzzz_0, \
                                         tg_xxyzz_xyzzz_1, tg_xxyzz_xzzz_1, tg_xxyzz_xzzzz_0, tg_xxyzz_xzzzz_1, \
                                         tg_xxyzz_yyyy_1, tg_xxyzz_yyyyy_0, tg_xxyzz_yyyyy_1, tg_xxyzz_yyyyz_0, \
                                         tg_xxyzz_yyyyz_1, tg_xxyzz_yyyz_1, tg_xxyzz_yyyzz_0, tg_xxyzz_yyyzz_1, \
                                         tg_xxyzz_yyzz_1, tg_xxyzz_yyzzz_0, tg_xxyzz_yyzzz_1, tg_xxyzz_yzzz_1, \
                                         tg_xxyzz_yzzzz_0, tg_xxyzz_yzzzz_1, tg_xxyzz_zzzz_1, tg_xxyzz_zzzzz_0, \
                                         tg_xxyzz_zzzzz_1, tg_xxzz_xxxxx_0, tg_xxzz_xxxxx_1, tg_xxzz_xxxxy_0, tg_xxzz_xxxxy_1, \
                                         tg_xxzz_xxxxz_0, tg_xxzz_xxxxz_1, tg_xxzz_xxxyy_0, tg_xxzz_xxxyy_1, tg_xxzz_xxxyz_0, \
                                         tg_xxzz_xxxyz_1, tg_xxzz_xxxzz_0, tg_xxzz_xxxzz_1, tg_xxzz_xxyyy_0, tg_xxzz_xxyyy_1, \
                                         tg_xxzz_xxyyz_0, tg_xxzz_xxyyz_1, tg_xxzz_xxyzz_0, tg_xxzz_xxyzz_1, tg_xxzz_xxzzz_0, \
                                         tg_xxzz_xxzzz_1, tg_xxzz_xyyyy_0, tg_xxzz_xyyyy_1, tg_xxzz_xyyyz_0, tg_xxzz_xyyyz_1, \
                                         tg_xxzz_xyyzz_0, tg_xxzz_xyyzz_1, tg_xxzz_xyzzz_0, tg_xxzz_xyzzz_1, tg_xxzz_xzzzz_0, \
                                         tg_xxzz_xzzzz_1, tg_xxzz_yyyyy_0, tg_xxzz_yyyyy_1, tg_xxzz_yyyyz_0, tg_xxzz_yyyyz_1, \
                                         tg_xxzz_yyyzz_0, tg_xxzz_yyyzz_1, tg_xxzz_yyzzz_0, tg_xxzz_yyzzz_1, tg_xxzz_yzzzz_0, \
                                         tg_xxzz_yzzzz_1, tg_xxzz_zzzzz_0, tg_xxzz_zzzzz_1, tg_xxzzz_xxxx_1, \
                                         tg_xxzzz_xxxxx_0, tg_xxzzz_xxxxx_1, tg_xxzzz_xxxxy_0, tg_xxzzz_xxxxy_1, \
                                         tg_xxzzz_xxxxz_0, tg_xxzzz_xxxxz_1, tg_xxzzz_xxxy_1, tg_xxzzz_xxxyy_0, \
                                         tg_xxzzz_xxxyy_1, tg_xxzzz_xxxyz_0, tg_xxzzz_xxxyz_1, tg_xxzzz_xxxz_1, \
                                         tg_xxzzz_xxxzz_0, tg_xxzzz_xxxzz_1, tg_xxzzz_xxyy_1, tg_xxzzz_xxyyy_0, \
                                         tg_xxzzz_xxyyy_1, tg_xxzzz_xxyz_1, tg_xxzzz_xxzz_1, tg_xxzzz_xyyy_1, tg_xyyy_xxxxx_0, \
                                         tg_xyyy_xxxxx_1, tg_xyyy_xxxxy_0, tg_xyyy_xxxxy_1, tg_xyyy_xxxxz_0, tg_xyyy_xxxxz_1, \
                                         tg_xyyy_xxxyy_0, tg_xyyy_xxxyy_1, tg_xyyy_xxxyz_0, tg_xyyy_xxxyz_1, tg_xyyy_xxxzz_0, \
                                         tg_xyyy_xxxzz_1, tg_xyyy_xxyyy_0, tg_xyyy_xxyyy_1, tg_xyyy_xxyyz_0, tg_xyyy_xxyyz_1, \
                                         tg_xyyy_xxyzz_0, tg_xyyy_xxyzz_1, tg_xyyy_xxzzz_0, tg_xyyy_xxzzz_1, tg_xyyy_xyyyy_0, \
                                         tg_xyyy_xyyyy_1, tg_xyyy_xyyyz_0, tg_xyyy_xyyyz_1, tg_xyyy_xyyzz_0, tg_xyyy_xyyzz_1, \
                                         tg_xyyy_xyzzz_0, tg_xyyy_xyzzz_1, tg_xyyy_xzzzz_0, tg_xyyy_xzzzz_1, tg_xyyy_yyyyy_0, \
                                         tg_xyyy_yyyyy_1, tg_xyyy_yyyyz_0, tg_xyyy_yyyyz_1, tg_xyyy_yyyzz_0, tg_xyyy_yyyzz_1, \
                                         tg_xyyy_yyzzz_0, tg_xyyy_yyzzz_1, tg_xyyy_yzzzz_0, tg_xyyy_yzzzz_1, tg_xyyy_zzzzz_0, \
                                         tg_xyyy_zzzzz_1, tg_xyyz_xxxxx_0, tg_xyyz_xxxxx_1, tg_xyyz_xxxxy_0, tg_xyyz_xxxxy_1, \
                                         tg_xyyz_xxxxz_0, tg_xyyz_xxxxz_1, tg_xyyz_xxxyy_0, tg_xyyz_xxxyy_1, tg_xyyz_xxxyz_0, \
                                         tg_xyyz_xxxyz_1, tg_xyyz_xxxzz_0, tg_xyyz_xxxzz_1, tg_xyyz_xxyyy_0, tg_xyyz_xxyyy_1, \
                                         tg_xyyz_xxyyz_0, tg_xyyz_xxyyz_1, tg_xyyz_xxyzz_0, tg_xyyz_xxyzz_1, tg_xyyz_xxzzz_0, \
                                         tg_xyyz_xxzzz_1, tg_xyyz_xyyyy_0, tg_xyyz_xyyyy_1, tg_xyyz_xyyyz_0, tg_xyyz_xyyyz_1, \
                                         tg_xyyz_xyyzz_0, tg_xyyz_xyyzz_1, tg_xyyz_xyzzz_0, tg_xyyz_xyzzz_1, tg_xyyz_xzzzz_0, \
                                         tg_xyyz_xzzzz_1, tg_xyyz_yyyyy_0, tg_xyyz_yyyyy_1, tg_xyyz_yyyyz_0, tg_xyyz_yyyyz_1, \
                                         tg_xyyz_yyyzz_0, tg_xyyz_yyyzz_1, tg_xyyz_yyzzz_0, tg_xyyz_yyzzz_1, tg_xyyz_yzzzz_0, \
                                         tg_xyyz_yzzzz_1, tg_xyyz_zzzzz_0, tg_xyyz_zzzzz_1, tg_xyzz_xxxxx_0, tg_xyzz_xxxxx_1, \
                                         tg_xyzz_xxxxy_0, tg_xyzz_xxxxy_1, tg_xyzz_xxxxz_0, tg_xyzz_xxxxz_1, tg_xyzz_xxxyy_0, \
                                         tg_xyzz_xxxyy_1, tg_xyzz_xxxyz_0, tg_xyzz_xxxyz_1, tg_xyzz_xxxzz_0, tg_xyzz_xxxzz_1, \
                                         tg_xyzz_xxyyy_0, tg_xyzz_xxyyy_1, tg_xyzz_xxyyz_0, tg_xyzz_xxyyz_1, tg_xyzz_xxyzz_0, \
                                         tg_xyzz_xxyzz_1, tg_xyzz_xxzzz_0, tg_xyzz_xxzzz_1, tg_xyzz_xyyyy_0, tg_xyzz_xyyyy_1, \
                                         tg_xyzz_xyyyz_0, tg_xyzz_xyyyz_1, tg_xyzz_xyyzz_0, tg_xyzz_xyyzz_1, tg_xyzz_xyzzz_0, \
                                         tg_xyzz_xyzzz_1, tg_xyzz_xzzzz_0, tg_xyzz_xzzzz_1, tg_xyzz_yyyyy_0, tg_xyzz_yyyyy_1, \
                                         tg_xyzz_yyyyz_0, tg_xyzz_yyyyz_1, tg_xyzz_yyyzz_0, tg_xyzz_yyyzz_1, tg_xyzz_yyzzz_0, \
                                         tg_xyzz_yyzzz_1, tg_xyzz_yzzzz_0, tg_xyzz_yzzzz_1, tg_xyzz_zzzzz_0, tg_xyzz_zzzzz_1, \
                                         tg_xzzz_xxxxx_0, tg_xzzz_xxxxx_1, tg_xzzz_xxxxy_0, tg_xzzz_xxxxy_1, tg_xzzz_xxxxz_0, \
                                         tg_xzzz_xxxxz_1, tg_xzzz_xxxyy_0, tg_xzzz_xxxyy_1, tg_xzzz_xxxyz_0, tg_xzzz_xxxyz_1, \
                                         tg_xzzz_xxxzz_0, tg_xzzz_xxxzz_1, tg_xzzz_xxyyy_0, tg_xzzz_xxyyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxyz_xzzzz_0[j] = pb_x * tg_xxxyz_xzzzz_0[j] + fr * tg_xxxyz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_xzzzz_0[j] - tg_xxyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyz_zzzz_1[j];

                    tg_xxxxyz_yyyyy_0[j] = pb_x * tg_xxxyz_yyyyy_0[j] + fr * tg_xxxyz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyyyy_0[j] - tg_xxyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyz_yyyyz_0[j] = pb_x * tg_xxxyz_yyyyz_0[j] + fr * tg_xxxyz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyyyz_0[j] - tg_xxyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyz_yyyzz_0[j] = pb_x * tg_xxxyz_yyyzz_0[j] + fr * tg_xxxyz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyyzz_0[j] - tg_xxyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyz_yyzzz_0[j] = pb_x * tg_xxxyz_yyzzz_0[j] + fr * tg_xxxyz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yyzzz_0[j] - tg_xxyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyz_yzzzz_0[j] = pb_x * tg_xxxyz_yzzzz_0[j] + fr * tg_xxxyz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_yzzzz_0[j] - tg_xxyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyz_zzzzz_0[j] = pb_x * tg_xxxyz_zzzzz_0[j] + fr * tg_xxxyz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyz_zzzzz_0[j] - tg_xxyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxzz_xxxxx_0[j] = pb_x * tg_xxxzz_xxxxx_0[j] + fr * tg_xxxzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxx_0[j] - tg_xxzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxzz_xxxx_1[j];

                    tg_xxxxzz_xxxxy_0[j] = pb_x * tg_xxxzz_xxxxy_0[j] + fr * tg_xxxzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxy_0[j] - tg_xxzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzz_xxxy_1[j];

                    tg_xxxxzz_xxxxz_0[j] = pb_x * tg_xxxzz_xxxxz_0[j] + fr * tg_xxxzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxxz_0[j] - tg_xxzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzz_xxxz_1[j];

                    tg_xxxxzz_xxxyy_0[j] = pb_x * tg_xxxzz_xxxyy_0[j] + fr * tg_xxxzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxyy_0[j] - tg_xxzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzz_xxyy_1[j];

                    tg_xxxxzz_xxxyz_0[j] = pb_x * tg_xxxzz_xxxyz_0[j] + fr * tg_xxxzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxyz_0[j] - tg_xxzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzz_xxyz_1[j];

                    tg_xxxxzz_xxxzz_0[j] = pb_x * tg_xxxzz_xxxzz_0[j] + fr * tg_xxxzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxxzz_0[j] - tg_xxzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzz_xxzz_1[j];

                    tg_xxxxzz_xxyyy_0[j] = pb_x * tg_xxxzz_xxyyy_0[j] + fr * tg_xxxzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxyyy_0[j] - tg_xxzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xyyy_1[j];

                    tg_xxxxzz_xxyyz_0[j] = pb_x * tg_xxxzz_xxyyz_0[j] + fr * tg_xxxzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxyyz_0[j] - tg_xxzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xyyz_1[j];

                    tg_xxxxzz_xxyzz_0[j] = pb_x * tg_xxxzz_xxyzz_0[j] + fr * tg_xxxzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxyzz_0[j] - tg_xxzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xyzz_1[j];

                    tg_xxxxzz_xxzzz_0[j] = pb_x * tg_xxxzz_xxzzz_0[j] + fr * tg_xxxzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xxzzz_0[j] - tg_xxzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzz_xzzz_1[j];

                    tg_xxxxzz_xyyyy_0[j] = pb_x * tg_xxxzz_xyyyy_0[j] + fr * tg_xxxzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyyyy_0[j] - tg_xxzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yyyy_1[j];

                    tg_xxxxzz_xyyyz_0[j] = pb_x * tg_xxxzz_xyyyz_0[j] + fr * tg_xxxzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyyyz_0[j] - tg_xxzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yyyz_1[j];

                    tg_xxxxzz_xyyzz_0[j] = pb_x * tg_xxxzz_xyyzz_0[j] + fr * tg_xxxzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyyzz_0[j] - tg_xxzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yyzz_1[j];

                    tg_xxxxzz_xyzzz_0[j] = pb_x * tg_xxxzz_xyzzz_0[j] + fr * tg_xxxzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xyzzz_0[j] - tg_xxzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_yzzz_1[j];

                    tg_xxxxzz_xzzzz_0[j] = pb_x * tg_xxxzz_xzzzz_0[j] + fr * tg_xxxzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_xzzzz_0[j] - tg_xxzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzz_zzzz_1[j];

                    tg_xxxxzz_yyyyy_0[j] = pb_x * tg_xxxzz_yyyyy_0[j] + fr * tg_xxxzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyyyy_0[j] - tg_xxzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxzz_yyyyz_0[j] = pb_x * tg_xxxzz_yyyyz_0[j] + fr * tg_xxxzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyyyz_0[j] - tg_xxzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxzz_yyyzz_0[j] = pb_x * tg_xxxzz_yyyzz_0[j] + fr * tg_xxxzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyyzz_0[j] - tg_xxzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxzz_yyzzz_0[j] = pb_x * tg_xxxzz_yyzzz_0[j] + fr * tg_xxxzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yyzzz_0[j] - tg_xxzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxzz_yzzzz_0[j] = pb_x * tg_xxxzz_yzzzz_0[j] + fr * tg_xxxzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_yzzzz_0[j] - tg_xxzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxzz_zzzzz_0[j] = pb_x * tg_xxxzz_zzzzz_0[j] + fr * tg_xxxzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzz_zzzzz_0[j] - tg_xxzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyy_xxxxx_0[j] = pb_x * tg_xxyyy_xxxxx_0[j] + fr * tg_xxyyy_xxxxx_1[j] + fl1_fx * (tg_xyyy_xxxxx_0[j] - tg_xyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyy_xxxx_1[j];

                    tg_xxxyyy_xxxxy_0[j] = pb_x * tg_xxyyy_xxxxy_0[j] + fr * tg_xxyyy_xxxxy_1[j] + fl1_fx * (tg_xyyy_xxxxy_0[j] - tg_xyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyy_xxxy_1[j];

                    tg_xxxyyy_xxxxz_0[j] = pb_x * tg_xxyyy_xxxxz_0[j] + fr * tg_xxyyy_xxxxz_1[j] + fl1_fx * (tg_xyyy_xxxxz_0[j] - tg_xyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyy_xxxz_1[j];

                    tg_xxxyyy_xxxyy_0[j] = pb_x * tg_xxyyy_xxxyy_0[j] + fr * tg_xxyyy_xxxyy_1[j] + fl1_fx * (tg_xyyy_xxxyy_0[j] - tg_xyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyy_xxyy_1[j];

                    tg_xxxyyy_xxxyz_0[j] = pb_x * tg_xxyyy_xxxyz_0[j] + fr * tg_xxyyy_xxxyz_1[j] + fl1_fx * (tg_xyyy_xxxyz_0[j] - tg_xyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyy_xxyz_1[j];

                    tg_xxxyyy_xxxzz_0[j] = pb_x * tg_xxyyy_xxxzz_0[j] + fr * tg_xxyyy_xxxzz_1[j] + fl1_fx * (tg_xyyy_xxxzz_0[j] - tg_xyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyy_xxzz_1[j];

                    tg_xxxyyy_xxyyy_0[j] = pb_x * tg_xxyyy_xxyyy_0[j] + fr * tg_xxyyy_xxyyy_1[j] + fl1_fx * (tg_xyyy_xxyyy_0[j] - tg_xyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xyyy_1[j];

                    tg_xxxyyy_xxyyz_0[j] = pb_x * tg_xxyyy_xxyyz_0[j] + fr * tg_xxyyy_xxyyz_1[j] + fl1_fx * (tg_xyyy_xxyyz_0[j] - tg_xyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xyyz_1[j];

                    tg_xxxyyy_xxyzz_0[j] = pb_x * tg_xxyyy_xxyzz_0[j] + fr * tg_xxyyy_xxyzz_1[j] + fl1_fx * (tg_xyyy_xxyzz_0[j] - tg_xyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xyzz_1[j];

                    tg_xxxyyy_xxzzz_0[j] = pb_x * tg_xxyyy_xxzzz_0[j] + fr * tg_xxyyy_xxzzz_1[j] + fl1_fx * (tg_xyyy_xxzzz_0[j] - tg_xyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyy_xzzz_1[j];

                    tg_xxxyyy_xyyyy_0[j] = pb_x * tg_xxyyy_xyyyy_0[j] + fr * tg_xxyyy_xyyyy_1[j] + fl1_fx * (tg_xyyy_xyyyy_0[j] - tg_xyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yyyy_1[j];

                    tg_xxxyyy_xyyyz_0[j] = pb_x * tg_xxyyy_xyyyz_0[j] + fr * tg_xxyyy_xyyyz_1[j] + fl1_fx * (tg_xyyy_xyyyz_0[j] - tg_xyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yyyz_1[j];

                    tg_xxxyyy_xyyzz_0[j] = pb_x * tg_xxyyy_xyyzz_0[j] + fr * tg_xxyyy_xyyzz_1[j] + fl1_fx * (tg_xyyy_xyyzz_0[j] - tg_xyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yyzz_1[j];

                    tg_xxxyyy_xyzzz_0[j] = pb_x * tg_xxyyy_xyzzz_0[j] + fr * tg_xxyyy_xyzzz_1[j] + fl1_fx * (tg_xyyy_xyzzz_0[j] - tg_xyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_yzzz_1[j];

                    tg_xxxyyy_xzzzz_0[j] = pb_x * tg_xxyyy_xzzzz_0[j] + fr * tg_xxyyy_xzzzz_1[j] + fl1_fx * (tg_xyyy_xzzzz_0[j] - tg_xyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyy_zzzz_1[j];

                    tg_xxxyyy_yyyyy_0[j] = pb_x * tg_xxyyy_yyyyy_0[j] + fr * tg_xxyyy_yyyyy_1[j] + fl1_fx * (tg_xyyy_yyyyy_0[j] - tg_xyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyy_yyyyz_0[j] = pb_x * tg_xxyyy_yyyyz_0[j] + fr * tg_xxyyy_yyyyz_1[j] + fl1_fx * (tg_xyyy_yyyyz_0[j] - tg_xyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyy_yyyzz_0[j] = pb_x * tg_xxyyy_yyyzz_0[j] + fr * tg_xxyyy_yyyzz_1[j] + fl1_fx * (tg_xyyy_yyyzz_0[j] - tg_xyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyy_yyzzz_0[j] = pb_x * tg_xxyyy_yyzzz_0[j] + fr * tg_xxyyy_yyzzz_1[j] + fl1_fx * (tg_xyyy_yyzzz_0[j] - tg_xyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyy_yzzzz_0[j] = pb_x * tg_xxyyy_yzzzz_0[j] + fr * tg_xxyyy_yzzzz_1[j] + fl1_fx * (tg_xyyy_yzzzz_0[j] - tg_xyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyy_zzzzz_0[j] = pb_x * tg_xxyyy_zzzzz_0[j] + fr * tg_xxyyy_zzzzz_1[j] + fl1_fx * (tg_xyyy_zzzzz_0[j] - tg_xyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyz_xxxxx_0[j] = pb_x * tg_xxyyz_xxxxx_0[j] + fr * tg_xxyyz_xxxxx_1[j] + fl1_fx * (tg_xyyz_xxxxx_0[j] - tg_xyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyz_xxxx_1[j];

                    tg_xxxyyz_xxxxy_0[j] = pb_x * tg_xxyyz_xxxxy_0[j] + fr * tg_xxyyz_xxxxy_1[j] + fl1_fx * (tg_xyyz_xxxxy_0[j] - tg_xyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyz_xxxy_1[j];

                    tg_xxxyyz_xxxxz_0[j] = pb_x * tg_xxyyz_xxxxz_0[j] + fr * tg_xxyyz_xxxxz_1[j] + fl1_fx * (tg_xyyz_xxxxz_0[j] - tg_xyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyz_xxxz_1[j];

                    tg_xxxyyz_xxxyy_0[j] = pb_x * tg_xxyyz_xxxyy_0[j] + fr * tg_xxyyz_xxxyy_1[j] + fl1_fx * (tg_xyyz_xxxyy_0[j] - tg_xyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyz_xxyy_1[j];

                    tg_xxxyyz_xxxyz_0[j] = pb_x * tg_xxyyz_xxxyz_0[j] + fr * tg_xxyyz_xxxyz_1[j] + fl1_fx * (tg_xyyz_xxxyz_0[j] - tg_xyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyz_xxyz_1[j];

                    tg_xxxyyz_xxxzz_0[j] = pb_x * tg_xxyyz_xxxzz_0[j] + fr * tg_xxyyz_xxxzz_1[j] + fl1_fx * (tg_xyyz_xxxzz_0[j] - tg_xyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyz_xxzz_1[j];

                    tg_xxxyyz_xxyyy_0[j] = pb_x * tg_xxyyz_xxyyy_0[j] + fr * tg_xxyyz_xxyyy_1[j] + fl1_fx * (tg_xyyz_xxyyy_0[j] - tg_xyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xyyy_1[j];

                    tg_xxxyyz_xxyyz_0[j] = pb_x * tg_xxyyz_xxyyz_0[j] + fr * tg_xxyyz_xxyyz_1[j] + fl1_fx * (tg_xyyz_xxyyz_0[j] - tg_xyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xyyz_1[j];

                    tg_xxxyyz_xxyzz_0[j] = pb_x * tg_xxyyz_xxyzz_0[j] + fr * tg_xxyyz_xxyzz_1[j] + fl1_fx * (tg_xyyz_xxyzz_0[j] - tg_xyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xyzz_1[j];

                    tg_xxxyyz_xxzzz_0[j] = pb_x * tg_xxyyz_xxzzz_0[j] + fr * tg_xxyyz_xxzzz_1[j] + fl1_fx * (tg_xyyz_xxzzz_0[j] - tg_xyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyz_xzzz_1[j];

                    tg_xxxyyz_xyyyy_0[j] = pb_x * tg_xxyyz_xyyyy_0[j] + fr * tg_xxyyz_xyyyy_1[j] + fl1_fx * (tg_xyyz_xyyyy_0[j] - tg_xyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yyyy_1[j];

                    tg_xxxyyz_xyyyz_0[j] = pb_x * tg_xxyyz_xyyyz_0[j] + fr * tg_xxyyz_xyyyz_1[j] + fl1_fx * (tg_xyyz_xyyyz_0[j] - tg_xyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yyyz_1[j];

                    tg_xxxyyz_xyyzz_0[j] = pb_x * tg_xxyyz_xyyzz_0[j] + fr * tg_xxyyz_xyyzz_1[j] + fl1_fx * (tg_xyyz_xyyzz_0[j] - tg_xyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yyzz_1[j];

                    tg_xxxyyz_xyzzz_0[j] = pb_x * tg_xxyyz_xyzzz_0[j] + fr * tg_xxyyz_xyzzz_1[j] + fl1_fx * (tg_xyyz_xyzzz_0[j] - tg_xyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_yzzz_1[j];

                    tg_xxxyyz_xzzzz_0[j] = pb_x * tg_xxyyz_xzzzz_0[j] + fr * tg_xxyyz_xzzzz_1[j] + fl1_fx * (tg_xyyz_xzzzz_0[j] - tg_xyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyz_zzzz_1[j];

                    tg_xxxyyz_yyyyy_0[j] = pb_x * tg_xxyyz_yyyyy_0[j] + fr * tg_xxyyz_yyyyy_1[j] + fl1_fx * (tg_xyyz_yyyyy_0[j] - tg_xyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyz_yyyyz_0[j] = pb_x * tg_xxyyz_yyyyz_0[j] + fr * tg_xxyyz_yyyyz_1[j] + fl1_fx * (tg_xyyz_yyyyz_0[j] - tg_xyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyz_yyyzz_0[j] = pb_x * tg_xxyyz_yyyzz_0[j] + fr * tg_xxyyz_yyyzz_1[j] + fl1_fx * (tg_xyyz_yyyzz_0[j] - tg_xyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyz_yyzzz_0[j] = pb_x * tg_xxyyz_yyzzz_0[j] + fr * tg_xxyyz_yyzzz_1[j] + fl1_fx * (tg_xyyz_yyzzz_0[j] - tg_xyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyz_yzzzz_0[j] = pb_x * tg_xxyyz_yzzzz_0[j] + fr * tg_xxyyz_yzzzz_1[j] + fl1_fx * (tg_xyyz_yzzzz_0[j] - tg_xyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyz_zzzzz_0[j] = pb_x * tg_xxyyz_zzzzz_0[j] + fr * tg_xxyyz_zzzzz_1[j] + fl1_fx * (tg_xyyz_zzzzz_0[j] - tg_xyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyzz_xxxxx_0[j] = pb_x * tg_xxyzz_xxxxx_0[j] + fr * tg_xxyzz_xxxxx_1[j] + fl1_fx * (tg_xyzz_xxxxx_0[j] - tg_xyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyzz_xxxx_1[j];

                    tg_xxxyzz_xxxxy_0[j] = pb_x * tg_xxyzz_xxxxy_0[j] + fr * tg_xxyzz_xxxxy_1[j] + fl1_fx * (tg_xyzz_xxxxy_0[j] - tg_xyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzz_xxxy_1[j];

                    tg_xxxyzz_xxxxz_0[j] = pb_x * tg_xxyzz_xxxxz_0[j] + fr * tg_xxyzz_xxxxz_1[j] + fl1_fx * (tg_xyzz_xxxxz_0[j] - tg_xyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzz_xxxz_1[j];

                    tg_xxxyzz_xxxyy_0[j] = pb_x * tg_xxyzz_xxxyy_0[j] + fr * tg_xxyzz_xxxyy_1[j] + fl1_fx * (tg_xyzz_xxxyy_0[j] - tg_xyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzz_xxyy_1[j];

                    tg_xxxyzz_xxxyz_0[j] = pb_x * tg_xxyzz_xxxyz_0[j] + fr * tg_xxyzz_xxxyz_1[j] + fl1_fx * (tg_xyzz_xxxyz_0[j] - tg_xyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzz_xxyz_1[j];

                    tg_xxxyzz_xxxzz_0[j] = pb_x * tg_xxyzz_xxxzz_0[j] + fr * tg_xxyzz_xxxzz_1[j] + fl1_fx * (tg_xyzz_xxxzz_0[j] - tg_xyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzz_xxzz_1[j];

                    tg_xxxyzz_xxyyy_0[j] = pb_x * tg_xxyzz_xxyyy_0[j] + fr * tg_xxyzz_xxyyy_1[j] + fl1_fx * (tg_xyzz_xxyyy_0[j] - tg_xyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xyyy_1[j];

                    tg_xxxyzz_xxyyz_0[j] = pb_x * tg_xxyzz_xxyyz_0[j] + fr * tg_xxyzz_xxyyz_1[j] + fl1_fx * (tg_xyzz_xxyyz_0[j] - tg_xyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xyyz_1[j];

                    tg_xxxyzz_xxyzz_0[j] = pb_x * tg_xxyzz_xxyzz_0[j] + fr * tg_xxyzz_xxyzz_1[j] + fl1_fx * (tg_xyzz_xxyzz_0[j] - tg_xyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xyzz_1[j];

                    tg_xxxyzz_xxzzz_0[j] = pb_x * tg_xxyzz_xxzzz_0[j] + fr * tg_xxyzz_xxzzz_1[j] + fl1_fx * (tg_xyzz_xxzzz_0[j] - tg_xyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzz_xzzz_1[j];

                    tg_xxxyzz_xyyyy_0[j] = pb_x * tg_xxyzz_xyyyy_0[j] + fr * tg_xxyzz_xyyyy_1[j] + fl1_fx * (tg_xyzz_xyyyy_0[j] - tg_xyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yyyy_1[j];

                    tg_xxxyzz_xyyyz_0[j] = pb_x * tg_xxyzz_xyyyz_0[j] + fr * tg_xxyzz_xyyyz_1[j] + fl1_fx * (tg_xyzz_xyyyz_0[j] - tg_xyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yyyz_1[j];

                    tg_xxxyzz_xyyzz_0[j] = pb_x * tg_xxyzz_xyyzz_0[j] + fr * tg_xxyzz_xyyzz_1[j] + fl1_fx * (tg_xyzz_xyyzz_0[j] - tg_xyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yyzz_1[j];

                    tg_xxxyzz_xyzzz_0[j] = pb_x * tg_xxyzz_xyzzz_0[j] + fr * tg_xxyzz_xyzzz_1[j] + fl1_fx * (tg_xyzz_xyzzz_0[j] - tg_xyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_yzzz_1[j];

                    tg_xxxyzz_xzzzz_0[j] = pb_x * tg_xxyzz_xzzzz_0[j] + fr * tg_xxyzz_xzzzz_1[j] + fl1_fx * (tg_xyzz_xzzzz_0[j] - tg_xyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzz_zzzz_1[j];

                    tg_xxxyzz_yyyyy_0[j] = pb_x * tg_xxyzz_yyyyy_0[j] + fr * tg_xxyzz_yyyyy_1[j] + fl1_fx * (tg_xyzz_yyyyy_0[j] - tg_xyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyzz_yyyyz_0[j] = pb_x * tg_xxyzz_yyyyz_0[j] + fr * tg_xxyzz_yyyyz_1[j] + fl1_fx * (tg_xyzz_yyyyz_0[j] - tg_xyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyzz_yyyzz_0[j] = pb_x * tg_xxyzz_yyyzz_0[j] + fr * tg_xxyzz_yyyzz_1[j] + fl1_fx * (tg_xyzz_yyyzz_0[j] - tg_xyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyzz_yyzzz_0[j] = pb_x * tg_xxyzz_yyzzz_0[j] + fr * tg_xxyzz_yyzzz_1[j] + fl1_fx * (tg_xyzz_yyzzz_0[j] - tg_xyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyzz_yzzzz_0[j] = pb_x * tg_xxyzz_yzzzz_0[j] + fr * tg_xxyzz_yzzzz_1[j] + fl1_fx * (tg_xyzz_yzzzz_0[j] - tg_xyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyzz_zzzzz_0[j] = pb_x * tg_xxyzz_zzzzz_0[j] + fr * tg_xxyzz_zzzzz_1[j] + fl1_fx * (tg_xyzz_zzzzz_0[j] - tg_xyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxzzz_xxxxx_0[j] = pb_x * tg_xxzzz_xxxxx_0[j] + fr * tg_xxzzz_xxxxx_1[j] + fl1_fx * (tg_xzzz_xxxxx_0[j] - tg_xzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzzz_xxxx_1[j];

                    tg_xxxzzz_xxxxy_0[j] = pb_x * tg_xxzzz_xxxxy_0[j] + fr * tg_xxzzz_xxxxy_1[j] + fl1_fx * (tg_xzzz_xxxxy_0[j] - tg_xzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzz_xxxy_1[j];

                    tg_xxxzzz_xxxxz_0[j] = pb_x * tg_xxzzz_xxxxz_0[j] + fr * tg_xxzzz_xxxxz_1[j] + fl1_fx * (tg_xzzz_xxxxz_0[j] - tg_xzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzz_xxxz_1[j];

                    tg_xxxzzz_xxxyy_0[j] = pb_x * tg_xxzzz_xxxyy_0[j] + fr * tg_xxzzz_xxxyy_1[j] + fl1_fx * (tg_xzzz_xxxyy_0[j] - tg_xzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzz_xxyy_1[j];

                    tg_xxxzzz_xxxyz_0[j] = pb_x * tg_xxzzz_xxxyz_0[j] + fr * tg_xxzzz_xxxyz_1[j] + fl1_fx * (tg_xzzz_xxxyz_0[j] - tg_xzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzz_xxyz_1[j];

                    tg_xxxzzz_xxxzz_0[j] = pb_x * tg_xxzzz_xxxzz_0[j] + fr * tg_xxzzz_xxxzz_1[j] + fl1_fx * (tg_xzzz_xxxzz_0[j] - tg_xzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzz_xxzz_1[j];

                    tg_xxxzzz_xxyyy_0[j] = pb_x * tg_xxzzz_xxyyy_0[j] + fr * tg_xxzzz_xxyyy_1[j] + fl1_fx * (tg_xzzz_xxyyy_0[j] - tg_xzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISH_196_294(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xyzzz_xyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 285); 

                auto tg_xyzzz_xyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 286); 

                auto tg_xyzzz_xzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 287); 

                auto tg_xyzzz_yyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 288); 

                auto tg_xyzzz_yyyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 289); 

                auto tg_xyzzz_yyyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 290); 

                auto tg_xyzzz_yyzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 291); 

                auto tg_xyzzz_yzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 292); 

                auto tg_xyzzz_zzzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 293); 

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

                auto tg_xyzzz_xyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 285); 

                auto tg_xyzzz_xyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 286); 

                auto tg_xyzzz_xzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 287); 

                auto tg_xyzzz_yyyyy_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 288); 

                auto tg_xyzzz_yyyyz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 289); 

                auto tg_xyzzz_yyyzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 290); 

                auto tg_xyzzz_yyzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 291); 

                auto tg_xyzzz_yzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 292); 

                auto tg_xyzzz_zzzzz_1 = primBuffer[pidx_g_5_5_m1].data(441 * idx + 293); 

                auto tg_xzzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 196); 

                auto tg_xzzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 197); 

                auto tg_xzzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 198); 

                auto tg_xzzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 199); 

                auto tg_xzzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 200); 

                auto tg_xzzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 201); 

                auto tg_xzzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 202); 

                auto tg_xzzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 203); 

                auto tg_xzzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 204); 

                auto tg_xzzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 205); 

                auto tg_xzzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 206); 

                auto tg_xzzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 207); 

                auto tg_xzzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 208); 

                auto tg_xzzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 209); 

                auto tg_yyyy_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 210); 

                auto tg_yyyy_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 211); 

                auto tg_yyyy_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 212); 

                auto tg_yyyy_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 213); 

                auto tg_yyyy_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 214); 

                auto tg_yyyy_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 215); 

                auto tg_yyyy_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 216); 

                auto tg_yyyy_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 217); 

                auto tg_yyyy_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 218); 

                auto tg_yyyy_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 219); 

                auto tg_yyyy_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 220); 

                auto tg_yyyy_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 221); 

                auto tg_yyyy_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 222); 

                auto tg_yyyy_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 223); 

                auto tg_yyyy_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 224); 

                auto tg_yyyy_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 225); 

                auto tg_yyyy_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 226); 

                auto tg_yyyy_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 227); 

                auto tg_yyyy_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 228); 

                auto tg_yyyy_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 229); 

                auto tg_yyyy_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 230); 

                auto tg_yyyz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 231); 

                auto tg_yyyz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 232); 

                auto tg_yyyz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 233); 

                auto tg_yyyz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 234); 

                auto tg_yyyz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 235); 

                auto tg_yyyz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 236); 

                auto tg_yyyz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 237); 

                auto tg_yyyz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 238); 

                auto tg_yyyz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 239); 

                auto tg_yyyz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 240); 

                auto tg_yyyz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 241); 

                auto tg_yyyz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 242); 

                auto tg_yyyz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 243); 

                auto tg_yyyz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 244); 

                auto tg_yyyz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 245); 

                auto tg_yyyz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 246); 

                auto tg_yyyz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 247); 

                auto tg_yyyz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 248); 

                auto tg_yyyz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 249); 

                auto tg_yyyz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 250); 

                auto tg_yyyz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 251); 

                auto tg_yyzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 252); 

                auto tg_yyzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 253); 

                auto tg_yyzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 254); 

                auto tg_yyzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 255); 

                auto tg_yyzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 256); 

                auto tg_yyzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 257); 

                auto tg_yyzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 258); 

                auto tg_yyzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 259); 

                auto tg_yyzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 260); 

                auto tg_yyzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 261); 

                auto tg_yyzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 262); 

                auto tg_yyzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 263); 

                auto tg_yyzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 264); 

                auto tg_yyzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 265); 

                auto tg_yyzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 266); 

                auto tg_yyzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 267); 

                auto tg_yyzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 268); 

                auto tg_yyzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 269); 

                auto tg_yyzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 270); 

                auto tg_yyzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 271); 

                auto tg_yyzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 272); 

                auto tg_yzzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 273); 

                auto tg_yzzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 274); 

                auto tg_yzzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 275); 

                auto tg_yzzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 276); 

                auto tg_yzzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 277); 

                auto tg_yzzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 278); 

                auto tg_yzzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 279); 

                auto tg_yzzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 280); 

                auto tg_yzzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 281); 

                auto tg_yzzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 282); 

                auto tg_yzzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 283); 

                auto tg_yzzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 284); 

                auto tg_yzzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 285); 

                auto tg_yzzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 286); 

                auto tg_yzzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 287); 

                auto tg_yzzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 288); 

                auto tg_yzzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 289); 

                auto tg_yzzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 290); 

                auto tg_yzzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 291); 

                auto tg_yzzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 292); 

                auto tg_yzzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 293); 

                auto tg_xzzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 196); 

                auto tg_xzzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 197); 

                auto tg_xzzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 198); 

                auto tg_xzzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 199); 

                auto tg_xzzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 200); 

                auto tg_xzzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 201); 

                auto tg_xzzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 202); 

                auto tg_xzzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 203); 

                auto tg_xzzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 204); 

                auto tg_xzzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 205); 

                auto tg_xzzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 206); 

                auto tg_xzzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 207); 

                auto tg_xzzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 208); 

                auto tg_xzzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 209); 

                auto tg_yyyy_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 210); 

                auto tg_yyyy_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 211); 

                auto tg_yyyy_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 212); 

                auto tg_yyyy_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 213); 

                auto tg_yyyy_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 214); 

                auto tg_yyyy_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 215); 

                auto tg_yyyy_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 216); 

                auto tg_yyyy_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 217); 

                auto tg_yyyy_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 218); 

                auto tg_yyyy_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 219); 

                auto tg_yyyy_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 220); 

                auto tg_yyyy_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 221); 

                auto tg_yyyy_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 222); 

                auto tg_yyyy_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 223); 

                auto tg_yyyy_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 224); 

                auto tg_yyyy_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 225); 

                auto tg_yyyy_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 226); 

                auto tg_yyyy_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 227); 

                auto tg_yyyy_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 228); 

                auto tg_yyyy_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 229); 

                auto tg_yyyy_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 230); 

                auto tg_yyyz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 231); 

                auto tg_yyyz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 232); 

                auto tg_yyyz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 233); 

                auto tg_yyyz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 234); 

                auto tg_yyyz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 235); 

                auto tg_yyyz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 236); 

                auto tg_yyyz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 237); 

                auto tg_yyyz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 238); 

                auto tg_yyyz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 239); 

                auto tg_yyyz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 240); 

                auto tg_yyyz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 241); 

                auto tg_yyyz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 242); 

                auto tg_yyyz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 243); 

                auto tg_yyyz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 244); 

                auto tg_yyyz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 245); 

                auto tg_yyyz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 246); 

                auto tg_yyyz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 247); 

                auto tg_yyyz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 248); 

                auto tg_yyyz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 249); 

                auto tg_yyyz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 250); 

                auto tg_yyyz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 251); 

                auto tg_yyzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 252); 

                auto tg_yyzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 253); 

                auto tg_yyzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 254); 

                auto tg_yyzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 255); 

                auto tg_yyzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 256); 

                auto tg_yyzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 257); 

                auto tg_yyzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 258); 

                auto tg_yyzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 259); 

                auto tg_yyzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 260); 

                auto tg_yyzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 261); 

                auto tg_yyzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 262); 

                auto tg_yyzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 263); 

                auto tg_yyzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 264); 

                auto tg_yyzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 265); 

                auto tg_yyzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 266); 

                auto tg_yyzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 267); 

                auto tg_yyzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 268); 

                auto tg_yyzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 269); 

                auto tg_yyzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 270); 

                auto tg_yyzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 271); 

                auto tg_yyzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 272); 

                auto tg_yzzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 273); 

                auto tg_yzzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 274); 

                auto tg_yzzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 275); 

                auto tg_yzzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 276); 

                auto tg_yzzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 277); 

                auto tg_yzzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 278); 

                auto tg_yzzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 279); 

                auto tg_yzzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 280); 

                auto tg_yzzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 281); 

                auto tg_yzzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 282); 

                auto tg_yzzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 283); 

                auto tg_yzzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 284); 

                auto tg_yzzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 285); 

                auto tg_yzzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 286); 

                auto tg_yzzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 287); 

                auto tg_yzzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 288); 

                auto tg_yzzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 289); 

                auto tg_yzzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 290); 

                auto tg_yzzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 291); 

                auto tg_yzzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 292); 

                auto tg_yzzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 293); 

                auto tg_xxzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 142); 

                auto tg_xxzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 143); 

                auto tg_xxzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 144); 

                auto tg_xxzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 145); 

                auto tg_xxzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 146); 

                auto tg_xxzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 147); 

                auto tg_xxzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 148); 

                auto tg_xxzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 149); 

                auto tg_xyyyy_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 150); 

                auto tg_xyyyy_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 151); 

                auto tg_xyyyy_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 152); 

                auto tg_xyyyy_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 153); 

                auto tg_xyyyy_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 154); 

                auto tg_xyyyy_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 155); 

                auto tg_xyyyy_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 156); 

                auto tg_xyyyy_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 157); 

                auto tg_xyyyy_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 158); 

                auto tg_xyyyy_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 159); 

                auto tg_xyyyy_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 160); 

                auto tg_xyyyy_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 161); 

                auto tg_xyyyy_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 162); 

                auto tg_xyyyy_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 163); 

                auto tg_xyyyy_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 164); 

                auto tg_xyyyz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 165); 

                auto tg_xyyyz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 166); 

                auto tg_xyyyz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 167); 

                auto tg_xyyyz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 168); 

                auto tg_xyyyz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 169); 

                auto tg_xyyyz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 170); 

                auto tg_xyyyz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 171); 

                auto tg_xyyyz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 172); 

                auto tg_xyyyz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 173); 

                auto tg_xyyyz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 174); 

                auto tg_xyyyz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 175); 

                auto tg_xyyyz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 176); 

                auto tg_xyyyz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 177); 

                auto tg_xyyyz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 178); 

                auto tg_xyyyz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 179); 

                auto tg_xyyzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 180); 

                auto tg_xyyzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 181); 

                auto tg_xyyzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 182); 

                auto tg_xyyzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 183); 

                auto tg_xyyzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 184); 

                auto tg_xyyzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 185); 

                auto tg_xyyzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 186); 

                auto tg_xyyzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 187); 

                auto tg_xyyzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 188); 

                auto tg_xyyzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 189); 

                auto tg_xyyzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 190); 

                auto tg_xyyzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 191); 

                auto tg_xyyzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 192); 

                auto tg_xyyzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 193); 

                auto tg_xyyzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 194); 

                auto tg_xyzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 195); 

                auto tg_xyzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 196); 

                auto tg_xyzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 197); 

                auto tg_xyzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 198); 

                auto tg_xyzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 199); 

                auto tg_xyzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 200); 

                auto tg_xyzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 201); 

                auto tg_xyzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 202); 

                auto tg_xyzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 203); 

                auto tg_xyzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 204); 

                auto tg_xyzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 205); 

                auto tg_xyzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 206); 

                auto tg_xyzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 207); 

                auto tg_xyzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 208); 

                auto tg_xyzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 209); 

                // set up pointers to integrals

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

                auto tg_xxyzzz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 285); 

                auto tg_xxyzzz_xyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 286); 

                auto tg_xxyzzz_xzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 287); 

                auto tg_xxyzzz_yyyyy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 288); 

                auto tg_xxyzzz_yyyyz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 289); 

                auto tg_xxyzzz_yyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 290); 

                auto tg_xxyzzz_yyzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 291); 

                auto tg_xxyzzz_yzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 292); 

                auto tg_xxyzzz_zzzzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 293); 

                // Batch of Integrals (196,294)

                #pragma omp simd aligned(fxn, fza, tg_xxxzzz_xxyyz_0, tg_xxxzzz_xxyzz_0, tg_xxxzzz_xxzzz_0, \
                                         tg_xxxzzz_xyyyy_0, tg_xxxzzz_xyyyz_0, tg_xxxzzz_xyyzz_0, tg_xxxzzz_xyzzz_0, \
                                         tg_xxxzzz_xzzzz_0, tg_xxxzzz_yyyyy_0, tg_xxxzzz_yyyyz_0, tg_xxxzzz_yyyzz_0, \
                                         tg_xxxzzz_yyzzz_0, tg_xxxzzz_yzzzz_0, tg_xxxzzz_zzzzz_0, tg_xxyyyy_xxxxx_0, \
                                         tg_xxyyyy_xxxxy_0, tg_xxyyyy_xxxxz_0, tg_xxyyyy_xxxyy_0, tg_xxyyyy_xxxyz_0, \
                                         tg_xxyyyy_xxxzz_0, tg_xxyyyy_xxyyy_0, tg_xxyyyy_xxyyz_0, tg_xxyyyy_xxyzz_0, \
                                         tg_xxyyyy_xxzzz_0, tg_xxyyyy_xyyyy_0, tg_xxyyyy_xyyyz_0, tg_xxyyyy_xyyzz_0, \
                                         tg_xxyyyy_xyzzz_0, tg_xxyyyy_xzzzz_0, tg_xxyyyy_yyyyy_0, tg_xxyyyy_yyyyz_0, \
                                         tg_xxyyyy_yyyzz_0, tg_xxyyyy_yyzzz_0, tg_xxyyyy_yzzzz_0, tg_xxyyyy_zzzzz_0, \
                                         tg_xxyyyz_xxxxx_0, tg_xxyyyz_xxxxy_0, tg_xxyyyz_xxxxz_0, tg_xxyyyz_xxxyy_0, \
                                         tg_xxyyyz_xxxyz_0, tg_xxyyyz_xxxzz_0, tg_xxyyyz_xxyyy_0, tg_xxyyyz_xxyyz_0, \
                                         tg_xxyyyz_xxyzz_0, tg_xxyyyz_xxzzz_0, tg_xxyyyz_xyyyy_0, tg_xxyyyz_xyyyz_0, \
                                         tg_xxyyyz_xyyzz_0, tg_xxyyyz_xyzzz_0, tg_xxyyyz_xzzzz_0, tg_xxyyyz_yyyyy_0, \
                                         tg_xxyyyz_yyyyz_0, tg_xxyyyz_yyyzz_0, tg_xxyyyz_yyzzz_0, tg_xxyyyz_yzzzz_0, \
                                         tg_xxyyyz_zzzzz_0, tg_xxyyzz_xxxxx_0, tg_xxyyzz_xxxxy_0, tg_xxyyzz_xxxxz_0, \
                                         tg_xxyyzz_xxxyy_0, tg_xxyyzz_xxxyz_0, tg_xxyyzz_xxxzz_0, tg_xxyyzz_xxyyy_0, \
                                         tg_xxyyzz_xxyyz_0, tg_xxyyzz_xxyzz_0, tg_xxyyzz_xxzzz_0, tg_xxyyzz_xyyyy_0, \
                                         tg_xxyyzz_xyyyz_0, tg_xxyyzz_xyyzz_0, tg_xxyyzz_xyzzz_0, tg_xxyyzz_xzzzz_0, \
                                         tg_xxyyzz_yyyyy_0, tg_xxyyzz_yyyyz_0, tg_xxyyzz_yyyzz_0, tg_xxyyzz_yyzzz_0, \
                                         tg_xxyyzz_yzzzz_0, tg_xxyyzz_zzzzz_0, tg_xxyzzz_xxxxx_0, tg_xxyzzz_xxxxy_0, \
                                         tg_xxyzzz_xxxxz_0, tg_xxyzzz_xxxyy_0, tg_xxyzzz_xxxyz_0, tg_xxyzzz_xxxzz_0, \
                                         tg_xxyzzz_xxyyy_0, tg_xxyzzz_xxyyz_0, tg_xxyzzz_xxyzz_0, tg_xxyzzz_xxzzz_0, \
                                         tg_xxyzzz_xyyyy_0, tg_xxyzzz_xyyyz_0, tg_xxyzzz_xyyzz_0, tg_xxyzzz_xyzzz_0, \
                                         tg_xxyzzz_xzzzz_0, tg_xxyzzz_yyyyy_0, tg_xxyzzz_yyyyz_0, tg_xxyzzz_yyyzz_0, \
                                         tg_xxyzzz_yyzzz_0, tg_xxyzzz_yzzzz_0, tg_xxyzzz_zzzzz_0, tg_xxzzz_xxyyz_0, \
                                         tg_xxzzz_xxyyz_1, tg_xxzzz_xxyzz_0, tg_xxzzz_xxyzz_1, tg_xxzzz_xxzzz_0, \
                                         tg_xxzzz_xxzzz_1, tg_xxzzz_xyyyy_0, tg_xxzzz_xyyyy_1, tg_xxzzz_xyyyz_0, \
                                         tg_xxzzz_xyyyz_1, tg_xxzzz_xyyz_1, tg_xxzzz_xyyzz_0, tg_xxzzz_xyyzz_1, \
                                         tg_xxzzz_xyzz_1, tg_xxzzz_xyzzz_0, tg_xxzzz_xyzzz_1, tg_xxzzz_xzzz_1, \
                                         tg_xxzzz_xzzzz_0, tg_xxzzz_xzzzz_1, tg_xxzzz_yyyy_1, tg_xxzzz_yyyyy_0, \
                                         tg_xxzzz_yyyyy_1, tg_xxzzz_yyyyz_0, tg_xxzzz_yyyyz_1, tg_xxzzz_yyyz_1, \
                                         tg_xxzzz_yyyzz_0, tg_xxzzz_yyyzz_1, tg_xxzzz_yyzz_1, tg_xxzzz_yyzzz_0, \
                                         tg_xxzzz_yyzzz_1, tg_xxzzz_yzzz_1, tg_xxzzz_yzzzz_0, tg_xxzzz_yzzzz_1, \
                                         tg_xxzzz_zzzz_1, tg_xxzzz_zzzzz_0, tg_xxzzz_zzzzz_1, tg_xyyyy_xxxx_1, \
                                         tg_xyyyy_xxxxx_0, tg_xyyyy_xxxxx_1, tg_xyyyy_xxxxy_0, tg_xyyyy_xxxxy_1, \
                                         tg_xyyyy_xxxxz_0, tg_xyyyy_xxxxz_1, tg_xyyyy_xxxy_1, tg_xyyyy_xxxyy_0, \
                                         tg_xyyyy_xxxyy_1, tg_xyyyy_xxxyz_0, tg_xyyyy_xxxyz_1, tg_xyyyy_xxxz_1, \
                                         tg_xyyyy_xxxzz_0, tg_xyyyy_xxxzz_1, tg_xyyyy_xxyy_1, tg_xyyyy_xxyyy_0, \
                                         tg_xyyyy_xxyyy_1, tg_xyyyy_xxyyz_0, tg_xyyyy_xxyyz_1, tg_xyyyy_xxyz_1, \
                                         tg_xyyyy_xxyzz_0, tg_xyyyy_xxyzz_1, tg_xyyyy_xxzz_1, tg_xyyyy_xxzzz_0, \
                                         tg_xyyyy_xxzzz_1, tg_xyyyy_xyyy_1, tg_xyyyy_xyyyy_0, tg_xyyyy_xyyyy_1, \
                                         tg_xyyyy_xyyyz_0, tg_xyyyy_xyyyz_1, tg_xyyyy_xyyz_1, tg_xyyyy_xyyzz_0, \
                                         tg_xyyyy_xyyzz_1, tg_xyyyy_xyzz_1, tg_xyyyy_xyzzz_0, tg_xyyyy_xyzzz_1, \
                                         tg_xyyyy_xzzz_1, tg_xyyyy_xzzzz_0, tg_xyyyy_xzzzz_1, tg_xyyyy_yyyy_1, \
                                         tg_xyyyy_yyyyy_0, tg_xyyyy_yyyyy_1, tg_xyyyy_yyyyz_0, tg_xyyyy_yyyyz_1, \
                                         tg_xyyyy_yyyz_1, tg_xyyyy_yyyzz_0, tg_xyyyy_yyyzz_1, tg_xyyyy_yyzz_1, \
                                         tg_xyyyy_yyzzz_0, tg_xyyyy_yyzzz_1, tg_xyyyy_yzzz_1, tg_xyyyy_yzzzz_0, \
                                         tg_xyyyy_yzzzz_1, tg_xyyyy_zzzz_1, tg_xyyyy_zzzzz_0, tg_xyyyy_zzzzz_1, \
                                         tg_xyyyz_xxxx_1, tg_xyyyz_xxxxx_0, tg_xyyyz_xxxxx_1, tg_xyyyz_xxxxy_0, \
                                         tg_xyyyz_xxxxy_1, tg_xyyyz_xxxxz_0, tg_xyyyz_xxxxz_1, tg_xyyyz_xxxy_1, \
                                         tg_xyyyz_xxxyy_0, tg_xyyyz_xxxyy_1, tg_xyyyz_xxxyz_0, tg_xyyyz_xxxyz_1, \
                                         tg_xyyyz_xxxz_1, tg_xyyyz_xxxzz_0, tg_xyyyz_xxxzz_1, tg_xyyyz_xxyy_1, \
                                         tg_xyyyz_xxyyy_0, tg_xyyyz_xxyyy_1, tg_xyyyz_xxyyz_0, tg_xyyyz_xxyyz_1, \
                                         tg_xyyyz_xxyz_1, tg_xyyyz_xxyzz_0, tg_xyyyz_xxyzz_1, tg_xyyyz_xxzz_1, \
                                         tg_xyyyz_xxzzz_0, tg_xyyyz_xxzzz_1, tg_xyyyz_xyyy_1, tg_xyyyz_xyyyy_0, \
                                         tg_xyyyz_xyyyy_1, tg_xyyyz_xyyyz_0, tg_xyyyz_xyyyz_1, tg_xyyyz_xyyz_1, \
                                         tg_xyyyz_xyyzz_0, tg_xyyyz_xyyzz_1, tg_xyyyz_xyzz_1, tg_xyyyz_xyzzz_0, \
                                         tg_xyyyz_xyzzz_1, tg_xyyyz_xzzz_1, tg_xyyyz_xzzzz_0, tg_xyyyz_xzzzz_1, \
                                         tg_xyyyz_yyyy_1, tg_xyyyz_yyyyy_0, tg_xyyyz_yyyyy_1, tg_xyyyz_yyyyz_0, \
                                         tg_xyyyz_yyyyz_1, tg_xyyyz_yyyz_1, tg_xyyyz_yyyzz_0, tg_xyyyz_yyyzz_1, \
                                         tg_xyyyz_yyzz_1, tg_xyyyz_yyzzz_0, tg_xyyyz_yyzzz_1, tg_xyyyz_yzzz_1, \
                                         tg_xyyyz_yzzzz_0, tg_xyyyz_yzzzz_1, tg_xyyyz_zzzz_1, tg_xyyyz_zzzzz_0, \
                                         tg_xyyyz_zzzzz_1, tg_xyyzz_xxxx_1, tg_xyyzz_xxxxx_0, tg_xyyzz_xxxxx_1, \
                                         tg_xyyzz_xxxxy_0, tg_xyyzz_xxxxy_1, tg_xyyzz_xxxxz_0, tg_xyyzz_xxxxz_1, \
                                         tg_xyyzz_xxxy_1, tg_xyyzz_xxxyy_0, tg_xyyzz_xxxyy_1, tg_xyyzz_xxxyz_0, \
                                         tg_xyyzz_xxxyz_1, tg_xyyzz_xxxz_1, tg_xyyzz_xxxzz_0, tg_xyyzz_xxxzz_1, \
                                         tg_xyyzz_xxyy_1, tg_xyyzz_xxyyy_0, tg_xyyzz_xxyyy_1, tg_xyyzz_xxyyz_0, \
                                         tg_xyyzz_xxyyz_1, tg_xyyzz_xxyz_1, tg_xyyzz_xxyzz_0, tg_xyyzz_xxyzz_1, \
                                         tg_xyyzz_xxzz_1, tg_xyyzz_xxzzz_0, tg_xyyzz_xxzzz_1, tg_xyyzz_xyyy_1, \
                                         tg_xyyzz_xyyyy_0, tg_xyyzz_xyyyy_1, tg_xyyzz_xyyyz_0, tg_xyyzz_xyyyz_1, \
                                         tg_xyyzz_xyyz_1, tg_xyyzz_xyyzz_0, tg_xyyzz_xyyzz_1, tg_xyyzz_xyzz_1, \
                                         tg_xyyzz_xyzzz_0, tg_xyyzz_xyzzz_1, tg_xyyzz_xzzz_1, tg_xyyzz_xzzzz_0, \
                                         tg_xyyzz_xzzzz_1, tg_xyyzz_yyyy_1, tg_xyyzz_yyyyy_0, tg_xyyzz_yyyyy_1, \
                                         tg_xyyzz_yyyyz_0, tg_xyyzz_yyyyz_1, tg_xyyzz_yyyz_1, tg_xyyzz_yyyzz_0, \
                                         tg_xyyzz_yyyzz_1, tg_xyyzz_yyzz_1, tg_xyyzz_yyzzz_0, tg_xyyzz_yyzzz_1, \
                                         tg_xyyzz_yzzz_1, tg_xyyzz_yzzzz_0, tg_xyyzz_yzzzz_1, tg_xyyzz_zzzz_1, \
                                         tg_xyyzz_zzzzz_0, tg_xyyzz_zzzzz_1, tg_xyzzz_xxxx_1, tg_xyzzz_xxxxx_0, \
                                         tg_xyzzz_xxxxx_1, tg_xyzzz_xxxxy_0, tg_xyzzz_xxxxy_1, tg_xyzzz_xxxxz_0, \
                                         tg_xyzzz_xxxxz_1, tg_xyzzz_xxxy_1, tg_xyzzz_xxxyy_0, tg_xyzzz_xxxyy_1, \
                                         tg_xyzzz_xxxyz_0, tg_xyzzz_xxxyz_1, tg_xyzzz_xxxz_1, tg_xyzzz_xxxzz_0, \
                                         tg_xyzzz_xxxzz_1, tg_xyzzz_xxyy_1, tg_xyzzz_xxyyy_0, tg_xyzzz_xxyyy_1, \
                                         tg_xyzzz_xxyyz_0, tg_xyzzz_xxyyz_1, tg_xyzzz_xxyz_1, tg_xyzzz_xxyzz_0, \
                                         tg_xyzzz_xxyzz_1, tg_xyzzz_xxzz_1, tg_xyzzz_xxzzz_0, tg_xyzzz_xxzzz_1, \
                                         tg_xyzzz_xyyy_1, tg_xyzzz_xyyyy_0, tg_xyzzz_xyyyy_1, tg_xyzzz_xyyyz_0, \
                                         tg_xyzzz_xyyyz_1, tg_xyzzz_xyyz_1, tg_xyzzz_xyyzz_0, tg_xyzzz_xyyzz_1, \
                                         tg_xyzzz_xyzz_1, tg_xyzzz_xyzzz_0, tg_xyzzz_xyzzz_1, tg_xyzzz_xzzz_1, \
                                         tg_xyzzz_xzzzz_0, tg_xyzzz_xzzzz_1, tg_xyzzz_yyyy_1, tg_xyzzz_yyyyy_0, \
                                         tg_xyzzz_yyyyy_1, tg_xyzzz_yyyyz_0, tg_xyzzz_yyyyz_1, tg_xyzzz_yyyz_1, \
                                         tg_xyzzz_yyyzz_0, tg_xyzzz_yyyzz_1, tg_xyzzz_yyzz_1, tg_xyzzz_yyzzz_0, \
                                         tg_xyzzz_yyzzz_1, tg_xyzzz_yzzz_1, tg_xyzzz_yzzzz_0, tg_xyzzz_yzzzz_1, \
                                         tg_xyzzz_zzzz_1, tg_xyzzz_zzzzz_0, tg_xyzzz_zzzzz_1, tg_xzzz_xxyyz_0, \
                                         tg_xzzz_xxyyz_1, tg_xzzz_xxyzz_0, tg_xzzz_xxyzz_1, tg_xzzz_xxzzz_0, tg_xzzz_xxzzz_1, \
                                         tg_xzzz_xyyyy_0, tg_xzzz_xyyyy_1, tg_xzzz_xyyyz_0, tg_xzzz_xyyyz_1, tg_xzzz_xyyzz_0, \
                                         tg_xzzz_xyyzz_1, tg_xzzz_xyzzz_0, tg_xzzz_xyzzz_1, tg_xzzz_xzzzz_0, tg_xzzz_xzzzz_1, \
                                         tg_xzzz_yyyyy_0, tg_xzzz_yyyyy_1, tg_xzzz_yyyyz_0, tg_xzzz_yyyyz_1, tg_xzzz_yyyzz_0, \
                                         tg_xzzz_yyyzz_1, tg_xzzz_yyzzz_0, tg_xzzz_yyzzz_1, tg_xzzz_yzzzz_0, tg_xzzz_yzzzz_1, \
                                         tg_xzzz_zzzzz_0, tg_xzzz_zzzzz_1, tg_yyyy_xxxxx_0, tg_yyyy_xxxxx_1, tg_yyyy_xxxxy_0, \
                                         tg_yyyy_xxxxy_1, tg_yyyy_xxxxz_0, tg_yyyy_xxxxz_1, tg_yyyy_xxxyy_0, tg_yyyy_xxxyy_1, \
                                         tg_yyyy_xxxyz_0, tg_yyyy_xxxyz_1, tg_yyyy_xxxzz_0, tg_yyyy_xxxzz_1, tg_yyyy_xxyyy_0, \
                                         tg_yyyy_xxyyy_1, tg_yyyy_xxyyz_0, tg_yyyy_xxyyz_1, tg_yyyy_xxyzz_0, tg_yyyy_xxyzz_1, \
                                         tg_yyyy_xxzzz_0, tg_yyyy_xxzzz_1, tg_yyyy_xyyyy_0, tg_yyyy_xyyyy_1, tg_yyyy_xyyyz_0, \
                                         tg_yyyy_xyyyz_1, tg_yyyy_xyyzz_0, tg_yyyy_xyyzz_1, tg_yyyy_xyzzz_0, tg_yyyy_xyzzz_1, \
                                         tg_yyyy_xzzzz_0, tg_yyyy_xzzzz_1, tg_yyyy_yyyyy_0, tg_yyyy_yyyyy_1, tg_yyyy_yyyyz_0, \
                                         tg_yyyy_yyyyz_1, tg_yyyy_yyyzz_0, tg_yyyy_yyyzz_1, tg_yyyy_yyzzz_0, tg_yyyy_yyzzz_1, \
                                         tg_yyyy_yzzzz_0, tg_yyyy_yzzzz_1, tg_yyyy_zzzzz_0, tg_yyyy_zzzzz_1, tg_yyyz_xxxxx_0, \
                                         tg_yyyz_xxxxx_1, tg_yyyz_xxxxy_0, tg_yyyz_xxxxy_1, tg_yyyz_xxxxz_0, tg_yyyz_xxxxz_1, \
                                         tg_yyyz_xxxyy_0, tg_yyyz_xxxyy_1, tg_yyyz_xxxyz_0, tg_yyyz_xxxyz_1, tg_yyyz_xxxzz_0, \
                                         tg_yyyz_xxxzz_1, tg_yyyz_xxyyy_0, tg_yyyz_xxyyy_1, tg_yyyz_xxyyz_0, tg_yyyz_xxyyz_1, \
                                         tg_yyyz_xxyzz_0, tg_yyyz_xxyzz_1, tg_yyyz_xxzzz_0, tg_yyyz_xxzzz_1, tg_yyyz_xyyyy_0, \
                                         tg_yyyz_xyyyy_1, tg_yyyz_xyyyz_0, tg_yyyz_xyyyz_1, tg_yyyz_xyyzz_0, tg_yyyz_xyyzz_1, \
                                         tg_yyyz_xyzzz_0, tg_yyyz_xyzzz_1, tg_yyyz_xzzzz_0, tg_yyyz_xzzzz_1, tg_yyyz_yyyyy_0, \
                                         tg_yyyz_yyyyy_1, tg_yyyz_yyyyz_0, tg_yyyz_yyyyz_1, tg_yyyz_yyyzz_0, tg_yyyz_yyyzz_1, \
                                         tg_yyyz_yyzzz_0, tg_yyyz_yyzzz_1, tg_yyyz_yzzzz_0, tg_yyyz_yzzzz_1, tg_yyyz_zzzzz_0, \
                                         tg_yyyz_zzzzz_1, tg_yyzz_xxxxx_0, tg_yyzz_xxxxx_1, tg_yyzz_xxxxy_0, tg_yyzz_xxxxy_1, \
                                         tg_yyzz_xxxxz_0, tg_yyzz_xxxxz_1, tg_yyzz_xxxyy_0, tg_yyzz_xxxyy_1, tg_yyzz_xxxyz_0, \
                                         tg_yyzz_xxxyz_1, tg_yyzz_xxxzz_0, tg_yyzz_xxxzz_1, tg_yyzz_xxyyy_0, tg_yyzz_xxyyy_1, \
                                         tg_yyzz_xxyyz_0, tg_yyzz_xxyyz_1, tg_yyzz_xxyzz_0, tg_yyzz_xxyzz_1, tg_yyzz_xxzzz_0, \
                                         tg_yyzz_xxzzz_1, tg_yyzz_xyyyy_0, tg_yyzz_xyyyy_1, tg_yyzz_xyyyz_0, tg_yyzz_xyyyz_1, \
                                         tg_yyzz_xyyzz_0, tg_yyzz_xyyzz_1, tg_yyzz_xyzzz_0, tg_yyzz_xyzzz_1, tg_yyzz_xzzzz_0, \
                                         tg_yyzz_xzzzz_1, tg_yyzz_yyyyy_0, tg_yyzz_yyyyy_1, tg_yyzz_yyyyz_0, tg_yyzz_yyyyz_1, \
                                         tg_yyzz_yyyzz_0, tg_yyzz_yyyzz_1, tg_yyzz_yyzzz_0, tg_yyzz_yyzzz_1, tg_yyzz_yzzzz_0, \
                                         tg_yyzz_yzzzz_1, tg_yyzz_zzzzz_0, tg_yyzz_zzzzz_1, tg_yzzz_xxxxx_0, tg_yzzz_xxxxx_1, \
                                         tg_yzzz_xxxxy_0, tg_yzzz_xxxxy_1, tg_yzzz_xxxxz_0, tg_yzzz_xxxxz_1, tg_yzzz_xxxyy_0, \
                                         tg_yzzz_xxxyy_1, tg_yzzz_xxxyz_0, tg_yzzz_xxxyz_1, tg_yzzz_xxxzz_0, tg_yzzz_xxxzz_1, \
                                         tg_yzzz_xxyyy_0, tg_yzzz_xxyyy_1, tg_yzzz_xxyyz_0, tg_yzzz_xxyyz_1, tg_yzzz_xxyzz_0, \
                                         tg_yzzz_xxyzz_1, tg_yzzz_xxzzz_0, tg_yzzz_xxzzz_1, tg_yzzz_xyyyy_0, tg_yzzz_xyyyy_1, \
                                         tg_yzzz_xyyyz_0, tg_yzzz_xyyyz_1, tg_yzzz_xyyzz_0, tg_yzzz_xyyzz_1, tg_yzzz_xyzzz_0, \
                                         tg_yzzz_xyzzz_1, tg_yzzz_xzzzz_0, tg_yzzz_xzzzz_1, tg_yzzz_yyyyy_0, tg_yzzz_yyyyy_1, \
                                         tg_yzzz_yyyyz_0, tg_yzzz_yyyyz_1, tg_yzzz_yyyzz_0, tg_yzzz_yyyzz_1, tg_yzzz_yyzzz_0, \
                                         tg_yzzz_yyzzz_1, tg_yzzz_yzzzz_0, tg_yzzz_yzzzz_1, tg_yzzz_zzzzz_0, tg_yzzz_zzzzz_1, \
                                         wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxzzz_xxyyz_0[j] = pb_x * tg_xxzzz_xxyyz_0[j] + fr * tg_xxzzz_xxyyz_1[j] + fl1_fx * (tg_xzzz_xxyyz_0[j] - tg_xzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xyyz_1[j];

                    tg_xxxzzz_xxyzz_0[j] = pb_x * tg_xxzzz_xxyzz_0[j] + fr * tg_xxzzz_xxyzz_1[j] + fl1_fx * (tg_xzzz_xxyzz_0[j] - tg_xzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xyzz_1[j];

                    tg_xxxzzz_xxzzz_0[j] = pb_x * tg_xxzzz_xxzzz_0[j] + fr * tg_xxzzz_xxzzz_1[j] + fl1_fx * (tg_xzzz_xxzzz_0[j] - tg_xzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzz_xzzz_1[j];

                    tg_xxxzzz_xyyyy_0[j] = pb_x * tg_xxzzz_xyyyy_0[j] + fr * tg_xxzzz_xyyyy_1[j] + fl1_fx * (tg_xzzz_xyyyy_0[j] - tg_xzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yyyy_1[j];

                    tg_xxxzzz_xyyyz_0[j] = pb_x * tg_xxzzz_xyyyz_0[j] + fr * tg_xxzzz_xyyyz_1[j] + fl1_fx * (tg_xzzz_xyyyz_0[j] - tg_xzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yyyz_1[j];

                    tg_xxxzzz_xyyzz_0[j] = pb_x * tg_xxzzz_xyyzz_0[j] + fr * tg_xxzzz_xyyzz_1[j] + fl1_fx * (tg_xzzz_xyyzz_0[j] - tg_xzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yyzz_1[j];

                    tg_xxxzzz_xyzzz_0[j] = pb_x * tg_xxzzz_xyzzz_0[j] + fr * tg_xxzzz_xyzzz_1[j] + fl1_fx * (tg_xzzz_xyzzz_0[j] - tg_xzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_yzzz_1[j];

                    tg_xxxzzz_xzzzz_0[j] = pb_x * tg_xxzzz_xzzzz_0[j] + fr * tg_xxzzz_xzzzz_1[j] + fl1_fx * (tg_xzzz_xzzzz_0[j] - tg_xzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzz_zzzz_1[j];

                    tg_xxxzzz_yyyyy_0[j] = pb_x * tg_xxzzz_yyyyy_0[j] + fr * tg_xxzzz_yyyyy_1[j] + fl1_fx * (tg_xzzz_yyyyy_0[j] - tg_xzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxzzz_yyyyz_0[j] = pb_x * tg_xxzzz_yyyyz_0[j] + fr * tg_xxzzz_yyyyz_1[j] + fl1_fx * (tg_xzzz_yyyyz_0[j] - tg_xzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxzzz_yyyzz_0[j] = pb_x * tg_xxzzz_yyyzz_0[j] + fr * tg_xxzzz_yyyzz_1[j] + fl1_fx * (tg_xzzz_yyyzz_0[j] - tg_xzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxzzz_yyzzz_0[j] = pb_x * tg_xxzzz_yyzzz_0[j] + fr * tg_xxzzz_yyzzz_1[j] + fl1_fx * (tg_xzzz_yyzzz_0[j] - tg_xzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxzzz_yzzzz_0[j] = pb_x * tg_xxzzz_yzzzz_0[j] + fr * tg_xxzzz_yzzzz_1[j] + fl1_fx * (tg_xzzz_yzzzz_0[j] - tg_xzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxzzz_zzzzz_0[j] = pb_x * tg_xxzzz_zzzzz_0[j] + fr * tg_xxzzz_zzzzz_1[j] + fl1_fx * (tg_xzzz_zzzzz_0[j] - tg_xzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyy_xxxxx_0[j] = pb_x * tg_xyyyy_xxxxx_0[j] + fr * tg_xyyyy_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxx_0[j] - tg_yyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyy_xxxx_1[j];

                    tg_xxyyyy_xxxxy_0[j] = pb_x * tg_xyyyy_xxxxy_0[j] + fr * tg_xyyyy_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxy_0[j] - tg_yyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyy_xxxy_1[j];

                    tg_xxyyyy_xxxxz_0[j] = pb_x * tg_xyyyy_xxxxz_0[j] + fr * tg_xyyyy_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxxz_0[j] - tg_yyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyy_xxxz_1[j];

                    tg_xxyyyy_xxxyy_0[j] = pb_x * tg_xyyyy_xxxyy_0[j] + fr * tg_xyyyy_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxyy_0[j] - tg_yyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyy_xxyy_1[j];

                    tg_xxyyyy_xxxyz_0[j] = pb_x * tg_xyyyy_xxxyz_0[j] + fr * tg_xyyyy_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxyz_0[j] - tg_yyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyy_xxyz_1[j];

                    tg_xxyyyy_xxxzz_0[j] = pb_x * tg_xyyyy_xxxzz_0[j] + fr * tg_xyyyy_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxxzz_0[j] - tg_yyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyy_xxzz_1[j];

                    tg_xxyyyy_xxyyy_0[j] = pb_x * tg_xyyyy_xxyyy_0[j] + fr * tg_xyyyy_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxyyy_0[j] - tg_yyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xyyy_1[j];

                    tg_xxyyyy_xxyyz_0[j] = pb_x * tg_xyyyy_xxyyz_0[j] + fr * tg_xyyyy_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxyyz_0[j] - tg_yyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xyyz_1[j];

                    tg_xxyyyy_xxyzz_0[j] = pb_x * tg_xyyyy_xxyzz_0[j] + fr * tg_xyyyy_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxyzz_0[j] - tg_yyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xyzz_1[j];

                    tg_xxyyyy_xxzzz_0[j] = pb_x * tg_xyyyy_xxzzz_0[j] + fr * tg_xyyyy_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xxzzz_0[j] - tg_yyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyy_xzzz_1[j];

                    tg_xxyyyy_xyyyy_0[j] = pb_x * tg_xyyyy_xyyyy_0[j] + fr * tg_xyyyy_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyyyy_0[j] - tg_yyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yyyy_1[j];

                    tg_xxyyyy_xyyyz_0[j] = pb_x * tg_xyyyy_xyyyz_0[j] + fr * tg_xyyyy_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyyyz_0[j] - tg_yyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yyyz_1[j];

                    tg_xxyyyy_xyyzz_0[j] = pb_x * tg_xyyyy_xyyzz_0[j] + fr * tg_xyyyy_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyyzz_0[j] - tg_yyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yyzz_1[j];

                    tg_xxyyyy_xyzzz_0[j] = pb_x * tg_xyyyy_xyzzz_0[j] + fr * tg_xyyyy_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xyzzz_0[j] - tg_yyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_yzzz_1[j];

                    tg_xxyyyy_xzzzz_0[j] = pb_x * tg_xyyyy_xzzzz_0[j] + fr * tg_xyyyy_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_xzzzz_0[j] - tg_yyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyy_zzzz_1[j];

                    tg_xxyyyy_yyyyy_0[j] = pb_x * tg_xyyyy_yyyyy_0[j] + fr * tg_xyyyy_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyyyy_0[j] - tg_yyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyy_yyyyz_0[j] = pb_x * tg_xyyyy_yyyyz_0[j] + fr * tg_xyyyy_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyyyz_0[j] - tg_yyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyy_yyyzz_0[j] = pb_x * tg_xyyyy_yyyzz_0[j] + fr * tg_xyyyy_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyyzz_0[j] - tg_yyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyy_yyzzz_0[j] = pb_x * tg_xyyyy_yyzzz_0[j] + fr * tg_xyyyy_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yyzzz_0[j] - tg_yyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyy_yzzzz_0[j] = pb_x * tg_xyyyy_yzzzz_0[j] + fr * tg_xyyyy_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_yzzzz_0[j] - tg_yyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyy_zzzzz_0[j] = pb_x * tg_xyyyy_zzzzz_0[j] + fr * tg_xyyyy_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyy_zzzzz_0[j] - tg_yyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyz_xxxxx_0[j] = pb_x * tg_xyyyz_xxxxx_0[j] + fr * tg_xyyyz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxx_0[j] - tg_yyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyz_xxxx_1[j];

                    tg_xxyyyz_xxxxy_0[j] = pb_x * tg_xyyyz_xxxxy_0[j] + fr * tg_xyyyz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxy_0[j] - tg_yyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyz_xxxy_1[j];

                    tg_xxyyyz_xxxxz_0[j] = pb_x * tg_xyyyz_xxxxz_0[j] + fr * tg_xyyyz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxxz_0[j] - tg_yyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyz_xxxz_1[j];

                    tg_xxyyyz_xxxyy_0[j] = pb_x * tg_xyyyz_xxxyy_0[j] + fr * tg_xyyyz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxyy_0[j] - tg_yyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyz_xxyy_1[j];

                    tg_xxyyyz_xxxyz_0[j] = pb_x * tg_xyyyz_xxxyz_0[j] + fr * tg_xyyyz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxyz_0[j] - tg_yyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyz_xxyz_1[j];

                    tg_xxyyyz_xxxzz_0[j] = pb_x * tg_xyyyz_xxxzz_0[j] + fr * tg_xyyyz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxxzz_0[j] - tg_yyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyz_xxzz_1[j];

                    tg_xxyyyz_xxyyy_0[j] = pb_x * tg_xyyyz_xxyyy_0[j] + fr * tg_xyyyz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxyyy_0[j] - tg_yyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xyyy_1[j];

                    tg_xxyyyz_xxyyz_0[j] = pb_x * tg_xyyyz_xxyyz_0[j] + fr * tg_xyyyz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxyyz_0[j] - tg_yyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xyyz_1[j];

                    tg_xxyyyz_xxyzz_0[j] = pb_x * tg_xyyyz_xxyzz_0[j] + fr * tg_xyyyz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxyzz_0[j] - tg_yyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xyzz_1[j];

                    tg_xxyyyz_xxzzz_0[j] = pb_x * tg_xyyyz_xxzzz_0[j] + fr * tg_xyyyz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xxzzz_0[j] - tg_yyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyz_xzzz_1[j];

                    tg_xxyyyz_xyyyy_0[j] = pb_x * tg_xyyyz_xyyyy_0[j] + fr * tg_xyyyz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyyyy_0[j] - tg_yyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yyyy_1[j];

                    tg_xxyyyz_xyyyz_0[j] = pb_x * tg_xyyyz_xyyyz_0[j] + fr * tg_xyyyz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyyyz_0[j] - tg_yyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yyyz_1[j];

                    tg_xxyyyz_xyyzz_0[j] = pb_x * tg_xyyyz_xyyzz_0[j] + fr * tg_xyyyz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyyzz_0[j] - tg_yyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yyzz_1[j];

                    tg_xxyyyz_xyzzz_0[j] = pb_x * tg_xyyyz_xyzzz_0[j] + fr * tg_xyyyz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xyzzz_0[j] - tg_yyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_yzzz_1[j];

                    tg_xxyyyz_xzzzz_0[j] = pb_x * tg_xyyyz_xzzzz_0[j] + fr * tg_xyyyz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_xzzzz_0[j] - tg_yyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyz_zzzz_1[j];

                    tg_xxyyyz_yyyyy_0[j] = pb_x * tg_xyyyz_yyyyy_0[j] + fr * tg_xyyyz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyyyy_0[j] - tg_yyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyz_yyyyz_0[j] = pb_x * tg_xyyyz_yyyyz_0[j] + fr * tg_xyyyz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyyyz_0[j] - tg_yyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyz_yyyzz_0[j] = pb_x * tg_xyyyz_yyyzz_0[j] + fr * tg_xyyyz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyyzz_0[j] - tg_yyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyz_yyzzz_0[j] = pb_x * tg_xyyyz_yyzzz_0[j] + fr * tg_xyyyz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yyzzz_0[j] - tg_yyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyz_yzzzz_0[j] = pb_x * tg_xyyyz_yzzzz_0[j] + fr * tg_xyyyz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_yzzzz_0[j] - tg_yyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyz_zzzzz_0[j] = pb_x * tg_xyyyz_zzzzz_0[j] + fr * tg_xyyyz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyz_zzzzz_0[j] - tg_yyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyzz_xxxxx_0[j] = pb_x * tg_xyyzz_xxxxx_0[j] + fr * tg_xyyzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxx_0[j] - tg_yyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyzz_xxxx_1[j];

                    tg_xxyyzz_xxxxy_0[j] = pb_x * tg_xyyzz_xxxxy_0[j] + fr * tg_xyyzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxy_0[j] - tg_yyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzz_xxxy_1[j];

                    tg_xxyyzz_xxxxz_0[j] = pb_x * tg_xyyzz_xxxxz_0[j] + fr * tg_xyyzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxxz_0[j] - tg_yyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzz_xxxz_1[j];

                    tg_xxyyzz_xxxyy_0[j] = pb_x * tg_xyyzz_xxxyy_0[j] + fr * tg_xyyzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxyy_0[j] - tg_yyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzz_xxyy_1[j];

                    tg_xxyyzz_xxxyz_0[j] = pb_x * tg_xyyzz_xxxyz_0[j] + fr * tg_xyyzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxyz_0[j] - tg_yyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzz_xxyz_1[j];

                    tg_xxyyzz_xxxzz_0[j] = pb_x * tg_xyyzz_xxxzz_0[j] + fr * tg_xyyzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxxzz_0[j] - tg_yyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzz_xxzz_1[j];

                    tg_xxyyzz_xxyyy_0[j] = pb_x * tg_xyyzz_xxyyy_0[j] + fr * tg_xyyzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxyyy_0[j] - tg_yyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xyyy_1[j];

                    tg_xxyyzz_xxyyz_0[j] = pb_x * tg_xyyzz_xxyyz_0[j] + fr * tg_xyyzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxyyz_0[j] - tg_yyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xyyz_1[j];

                    tg_xxyyzz_xxyzz_0[j] = pb_x * tg_xyyzz_xxyzz_0[j] + fr * tg_xyyzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxyzz_0[j] - tg_yyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xyzz_1[j];

                    tg_xxyyzz_xxzzz_0[j] = pb_x * tg_xyyzz_xxzzz_0[j] + fr * tg_xyyzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xxzzz_0[j] - tg_yyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzz_xzzz_1[j];

                    tg_xxyyzz_xyyyy_0[j] = pb_x * tg_xyyzz_xyyyy_0[j] + fr * tg_xyyzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyyyy_0[j] - tg_yyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yyyy_1[j];

                    tg_xxyyzz_xyyyz_0[j] = pb_x * tg_xyyzz_xyyyz_0[j] + fr * tg_xyyzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyyyz_0[j] - tg_yyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yyyz_1[j];

                    tg_xxyyzz_xyyzz_0[j] = pb_x * tg_xyyzz_xyyzz_0[j] + fr * tg_xyyzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyyzz_0[j] - tg_yyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yyzz_1[j];

                    tg_xxyyzz_xyzzz_0[j] = pb_x * tg_xyyzz_xyzzz_0[j] + fr * tg_xyyzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xyzzz_0[j] - tg_yyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_yzzz_1[j];

                    tg_xxyyzz_xzzzz_0[j] = pb_x * tg_xyyzz_xzzzz_0[j] + fr * tg_xyyzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_xzzzz_0[j] - tg_yyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzz_zzzz_1[j];

                    tg_xxyyzz_yyyyy_0[j] = pb_x * tg_xyyzz_yyyyy_0[j] + fr * tg_xyyzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyyyy_0[j] - tg_yyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyzz_yyyyz_0[j] = pb_x * tg_xyyzz_yyyyz_0[j] + fr * tg_xyyzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyyyz_0[j] - tg_yyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyzz_yyyzz_0[j] = pb_x * tg_xyyzz_yyyzz_0[j] + fr * tg_xyyzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyyzz_0[j] - tg_yyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyzz_yyzzz_0[j] = pb_x * tg_xyyzz_yyzzz_0[j] + fr * tg_xyyzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yyzzz_0[j] - tg_yyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyzz_yzzzz_0[j] = pb_x * tg_xyyzz_yzzzz_0[j] + fr * tg_xyyzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_yzzzz_0[j] - tg_yyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyzz_zzzzz_0[j] = pb_x * tg_xyyzz_zzzzz_0[j] + fr * tg_xyyzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzz_zzzzz_0[j] - tg_yyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyzzz_xxxxx_0[j] = pb_x * tg_xyzzz_xxxxx_0[j] + fr * tg_xyzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxx_0[j] - tg_yzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzzz_xxxx_1[j];

                    tg_xxyzzz_xxxxy_0[j] = pb_x * tg_xyzzz_xxxxy_0[j] + fr * tg_xyzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxy_0[j] - tg_yzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzz_xxxy_1[j];

                    tg_xxyzzz_xxxxz_0[j] = pb_x * tg_xyzzz_xxxxz_0[j] + fr * tg_xyzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxxz_0[j] - tg_yzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzz_xxxz_1[j];

                    tg_xxyzzz_xxxyy_0[j] = pb_x * tg_xyzzz_xxxyy_0[j] + fr * tg_xyzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxyy_0[j] - tg_yzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzz_xxyy_1[j];

                    tg_xxyzzz_xxxyz_0[j] = pb_x * tg_xyzzz_xxxyz_0[j] + fr * tg_xyzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxyz_0[j] - tg_yzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzz_xxyz_1[j];

                    tg_xxyzzz_xxxzz_0[j] = pb_x * tg_xyzzz_xxxzz_0[j] + fr * tg_xyzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxxzz_0[j] - tg_yzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzz_xxzz_1[j];

                    tg_xxyzzz_xxyyy_0[j] = pb_x * tg_xyzzz_xxyyy_0[j] + fr * tg_xyzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxyyy_0[j] - tg_yzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xyyy_1[j];

                    tg_xxyzzz_xxyyz_0[j] = pb_x * tg_xyzzz_xxyyz_0[j] + fr * tg_xyzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxyyz_0[j] - tg_yzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xyyz_1[j];

                    tg_xxyzzz_xxyzz_0[j] = pb_x * tg_xyzzz_xxyzz_0[j] + fr * tg_xyzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxyzz_0[j] - tg_yzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xyzz_1[j];

                    tg_xxyzzz_xxzzz_0[j] = pb_x * tg_xyzzz_xxzzz_0[j] + fr * tg_xyzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xxzzz_0[j] - tg_yzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzz_xzzz_1[j];

                    tg_xxyzzz_xyyyy_0[j] = pb_x * tg_xyzzz_xyyyy_0[j] + fr * tg_xyzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyyyy_0[j] - tg_yzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yyyy_1[j];

                    tg_xxyzzz_xyyyz_0[j] = pb_x * tg_xyzzz_xyyyz_0[j] + fr * tg_xyzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyyyz_0[j] - tg_yzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yyyz_1[j];

                    tg_xxyzzz_xyyzz_0[j] = pb_x * tg_xyzzz_xyyzz_0[j] + fr * tg_xyzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyyzz_0[j] - tg_yzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yyzz_1[j];

                    tg_xxyzzz_xyzzz_0[j] = pb_x * tg_xyzzz_xyzzz_0[j] + fr * tg_xyzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xyzzz_0[j] - tg_yzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_yzzz_1[j];

                    tg_xxyzzz_xzzzz_0[j] = pb_x * tg_xyzzz_xzzzz_0[j] + fr * tg_xyzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_xzzzz_0[j] - tg_yzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzz_zzzz_1[j];

                    tg_xxyzzz_yyyyy_0[j] = pb_x * tg_xyzzz_yyyyy_0[j] + fr * tg_xyzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyyyy_0[j] - tg_yzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyzzz_yyyyz_0[j] = pb_x * tg_xyzzz_yyyyz_0[j] + fr * tg_xyzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyyyz_0[j] - tg_yzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyzzz_yyyzz_0[j] = pb_x * tg_xyzzz_yyyzz_0[j] + fr * tg_xyzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyyzz_0[j] - tg_yzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyzzz_yyzzz_0[j] = pb_x * tg_xyzzz_yyzzz_0[j] + fr * tg_xyzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yyzzz_0[j] - tg_yzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyzzz_yzzzz_0[j] = pb_x * tg_xyzzz_yzzzz_0[j] + fr * tg_xyzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_yzzzz_0[j] - tg_yzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyzzz_zzzzz_0[j] = pb_x * tg_xyzzz_zzzzz_0[j] + fr * tg_xyzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzz_zzzzz_0[j] - tg_yzzz_zzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISH_294_392(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_zzzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 294); 

                auto tg_zzzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 295); 

                auto tg_zzzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 296); 

                auto tg_zzzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 297); 

                auto tg_zzzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 298); 

                auto tg_zzzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 299); 

                auto tg_zzzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 300); 

                auto tg_zzzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 301); 

                auto tg_zzzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 302); 

                auto tg_zzzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 303); 

                auto tg_zzzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 304); 

                auto tg_zzzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 305); 

                auto tg_zzzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 306); 

                auto tg_zzzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 307); 

                auto tg_zzzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 308); 

                auto tg_zzzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 309); 

                auto tg_zzzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 310); 

                auto tg_zzzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 311); 

                auto tg_zzzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 312); 

                auto tg_zzzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 313); 

                auto tg_zzzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 314); 

                auto tg_zzzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 294); 

                auto tg_zzzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 295); 

                auto tg_zzzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 296); 

                auto tg_zzzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 297); 

                auto tg_zzzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 298); 

                auto tg_zzzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 299); 

                auto tg_zzzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 300); 

                auto tg_zzzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 301); 

                auto tg_zzzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 302); 

                auto tg_zzzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 303); 

                auto tg_zzzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 304); 

                auto tg_zzzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 305); 

                auto tg_zzzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 306); 

                auto tg_zzzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 307); 

                auto tg_zzzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 308); 

                auto tg_zzzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 309); 

                auto tg_zzzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 310); 

                auto tg_zzzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 311); 

                auto tg_zzzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 312); 

                auto tg_zzzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 313); 

                auto tg_zzzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 314); 

                auto tg_xzzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 210); 

                auto tg_xzzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 211); 

                auto tg_xzzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 212); 

                auto tg_xzzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 213); 

                auto tg_xzzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 214); 

                auto tg_xzzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 215); 

                auto tg_xzzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 216); 

                auto tg_xzzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 217); 

                auto tg_xzzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 218); 

                auto tg_xzzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 219); 

                auto tg_xzzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 220); 

                auto tg_xzzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 221); 

                auto tg_xzzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 222); 

                auto tg_xzzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 223); 

                auto tg_xzzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 224); 

                auto tg_yyyyy_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 225); 

                auto tg_yyyyy_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 226); 

                auto tg_yyyyy_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 227); 

                auto tg_yyyyy_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 228); 

                auto tg_yyyyy_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 229); 

                auto tg_yyyyy_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 230); 

                auto tg_yyyyy_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 231); 

                auto tg_yyyyy_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 232); 

                auto tg_yyyyy_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 233); 

                auto tg_yyyyy_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 234); 

                auto tg_yyyyy_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 235); 

                auto tg_yyyyy_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 236); 

                auto tg_yyyyy_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 237); 

                auto tg_yyyyy_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 238); 

                auto tg_yyyyy_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 239); 

                auto tg_yyyyz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 240); 

                auto tg_yyyyz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 241); 

                auto tg_yyyyz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 242); 

                auto tg_yyyyz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 243); 

                auto tg_yyyyz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 244); 

                auto tg_yyyyz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 245); 

                auto tg_yyyyz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 246); 

                auto tg_yyyyz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 247); 

                auto tg_yyyyz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 248); 

                auto tg_yyyyz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 249); 

                auto tg_yyyyz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 250); 

                auto tg_yyyyz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 251); 

                auto tg_yyyyz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 252); 

                auto tg_yyyyz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 253); 

                auto tg_yyyyz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 254); 

                auto tg_yyyzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 255); 

                auto tg_yyyzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 256); 

                auto tg_yyyzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 257); 

                auto tg_yyyzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 258); 

                auto tg_yyyzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 259); 

                auto tg_yyyzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 260); 

                auto tg_yyyzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 261); 

                auto tg_yyyzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 262); 

                auto tg_yyyzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 263); 

                auto tg_yyyzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 264); 

                auto tg_yyyzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 265); 

                auto tg_yyyzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 266); 

                auto tg_yyyzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 267); 

                auto tg_yyyzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 268); 

                auto tg_yyyzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 269); 

                auto tg_yyzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 270); 

                auto tg_yyzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 271); 

                auto tg_yyzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 272); 

                auto tg_yyzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 273); 

                auto tg_yyzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 274); 

                auto tg_yyzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 275); 

                auto tg_yyzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 276); 

                auto tg_yyzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 277); 

                auto tg_yyzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 278); 

                auto tg_yyzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 279); 

                auto tg_yyzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 280); 

                auto tg_yyzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 281); 

                auto tg_yyzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 282); 

                auto tg_yyzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 283); 

                // set up pointers to integrals

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

                // Batch of Integrals (294,392)

                #pragma omp simd aligned(fxn, fza, tg_xxzzzz_xxxxx_0, tg_xxzzzz_xxxxy_0, tg_xxzzzz_xxxxz_0, \
                                         tg_xxzzzz_xxxyy_0, tg_xxzzzz_xxxyz_0, tg_xxzzzz_xxxzz_0, tg_xxzzzz_xxyyy_0, \
                                         tg_xxzzzz_xxyyz_0, tg_xxzzzz_xxyzz_0, tg_xxzzzz_xxzzz_0, tg_xxzzzz_xyyyy_0, \
                                         tg_xxzzzz_xyyyz_0, tg_xxzzzz_xyyzz_0, tg_xxzzzz_xyzzz_0, tg_xxzzzz_xzzzz_0, \
                                         tg_xxzzzz_yyyyy_0, tg_xxzzzz_yyyyz_0, tg_xxzzzz_yyyzz_0, tg_xxzzzz_yyzzz_0, \
                                         tg_xxzzzz_yzzzz_0, tg_xxzzzz_zzzzz_0, tg_xyyyyy_xxxxx_0, tg_xyyyyy_xxxxy_0, \
                                         tg_xyyyyy_xxxxz_0, tg_xyyyyy_xxxyy_0, tg_xyyyyy_xxxyz_0, tg_xyyyyy_xxxzz_0, \
                                         tg_xyyyyy_xxyyy_0, tg_xyyyyy_xxyyz_0, tg_xyyyyy_xxyzz_0, tg_xyyyyy_xxzzz_0, \
                                         tg_xyyyyy_xyyyy_0, tg_xyyyyy_xyyyz_0, tg_xyyyyy_xyyzz_0, tg_xyyyyy_xyzzz_0, \
                                         tg_xyyyyy_xzzzz_0, tg_xyyyyy_yyyyy_0, tg_xyyyyy_yyyyz_0, tg_xyyyyy_yyyzz_0, \
                                         tg_xyyyyy_yyzzz_0, tg_xyyyyy_yzzzz_0, tg_xyyyyy_zzzzz_0, tg_xyyyyz_xxxxx_0, \
                                         tg_xyyyyz_xxxxy_0, tg_xyyyyz_xxxxz_0, tg_xyyyyz_xxxyy_0, tg_xyyyyz_xxxyz_0, \
                                         tg_xyyyyz_xxxzz_0, tg_xyyyyz_xxyyy_0, tg_xyyyyz_xxyyz_0, tg_xyyyyz_xxyzz_0, \
                                         tg_xyyyyz_xxzzz_0, tg_xyyyyz_xyyyy_0, tg_xyyyyz_xyyyz_0, tg_xyyyyz_xyyzz_0, \
                                         tg_xyyyyz_xyzzz_0, tg_xyyyyz_xzzzz_0, tg_xyyyyz_yyyyy_0, tg_xyyyyz_yyyyz_0, \
                                         tg_xyyyyz_yyyzz_0, tg_xyyyyz_yyzzz_0, tg_xyyyyz_yzzzz_0, tg_xyyyyz_zzzzz_0, \
                                         tg_xyyyzz_xxxxx_0, tg_xyyyzz_xxxxy_0, tg_xyyyzz_xxxxz_0, tg_xyyyzz_xxxyy_0, \
                                         tg_xyyyzz_xxxyz_0, tg_xyyyzz_xxxzz_0, tg_xyyyzz_xxyyy_0, tg_xyyyzz_xxyyz_0, \
                                         tg_xyyyzz_xxyzz_0, tg_xyyyzz_xxzzz_0, tg_xyyyzz_xyyyy_0, tg_xyyyzz_xyyyz_0, \
                                         tg_xyyyzz_xyyzz_0, tg_xyyyzz_xyzzz_0, tg_xyyyzz_xzzzz_0, tg_xyyyzz_yyyyy_0, \
                                         tg_xyyyzz_yyyyz_0, tg_xyyyzz_yyyzz_0, tg_xyyyzz_yyzzz_0, tg_xyyyzz_yzzzz_0, \
                                         tg_xyyyzz_zzzzz_0, tg_xyyzzz_xxxxx_0, tg_xyyzzz_xxxxy_0, tg_xyyzzz_xxxxz_0, \
                                         tg_xyyzzz_xxxyy_0, tg_xyyzzz_xxxyz_0, tg_xyyzzz_xxxzz_0, tg_xyyzzz_xxyyy_0, \
                                         tg_xyyzzz_xxyyz_0, tg_xyyzzz_xxyzz_0, tg_xyyzzz_xxzzz_0, tg_xyyzzz_xyyyy_0, \
                                         tg_xyyzzz_xyyyz_0, tg_xyyzzz_xyyzz_0, tg_xyyzzz_xyzzz_0, tg_xzzzz_xxxx_1, \
                                         tg_xzzzz_xxxxx_0, tg_xzzzz_xxxxx_1, tg_xzzzz_xxxxy_0, tg_xzzzz_xxxxy_1, \
                                         tg_xzzzz_xxxxz_0, tg_xzzzz_xxxxz_1, tg_xzzzz_xxxy_1, tg_xzzzz_xxxyy_0, \
                                         tg_xzzzz_xxxyy_1, tg_xzzzz_xxxyz_0, tg_xzzzz_xxxyz_1, tg_xzzzz_xxxz_1, \
                                         tg_xzzzz_xxxzz_0, tg_xzzzz_xxxzz_1, tg_xzzzz_xxyy_1, tg_xzzzz_xxyyy_0, \
                                         tg_xzzzz_xxyyy_1, tg_xzzzz_xxyyz_0, tg_xzzzz_xxyyz_1, tg_xzzzz_xxyz_1, \
                                         tg_xzzzz_xxyzz_0, tg_xzzzz_xxyzz_1, tg_xzzzz_xxzz_1, tg_xzzzz_xxzzz_0, \
                                         tg_xzzzz_xxzzz_1, tg_xzzzz_xyyy_1, tg_xzzzz_xyyyy_0, tg_xzzzz_xyyyy_1, \
                                         tg_xzzzz_xyyyz_0, tg_xzzzz_xyyyz_1, tg_xzzzz_xyyz_1, tg_xzzzz_xyyzz_0, \
                                         tg_xzzzz_xyyzz_1, tg_xzzzz_xyzz_1, tg_xzzzz_xyzzz_0, tg_xzzzz_xyzzz_1, \
                                         tg_xzzzz_xzzz_1, tg_xzzzz_xzzzz_0, tg_xzzzz_xzzzz_1, tg_xzzzz_yyyy_1, \
                                         tg_xzzzz_yyyyy_0, tg_xzzzz_yyyyy_1, tg_xzzzz_yyyyz_0, tg_xzzzz_yyyyz_1, \
                                         tg_xzzzz_yyyz_1, tg_xzzzz_yyyzz_0, tg_xzzzz_yyyzz_1, tg_xzzzz_yyzz_1, \
                                         tg_xzzzz_yyzzz_0, tg_xzzzz_yyzzz_1, tg_xzzzz_yzzz_1, tg_xzzzz_yzzzz_0, \
                                         tg_xzzzz_yzzzz_1, tg_xzzzz_zzzz_1, tg_xzzzz_zzzzz_0, tg_xzzzz_zzzzz_1, \
                                         tg_yyyyy_xxxx_1, tg_yyyyy_xxxxx_0, tg_yyyyy_xxxxx_1, tg_yyyyy_xxxxy_0, \
                                         tg_yyyyy_xxxxy_1, tg_yyyyy_xxxxz_0, tg_yyyyy_xxxxz_1, tg_yyyyy_xxxy_1, \
                                         tg_yyyyy_xxxyy_0, tg_yyyyy_xxxyy_1, tg_yyyyy_xxxyz_0, tg_yyyyy_xxxyz_1, \
                                         tg_yyyyy_xxxz_1, tg_yyyyy_xxxzz_0, tg_yyyyy_xxxzz_1, tg_yyyyy_xxyy_1, \
                                         tg_yyyyy_xxyyy_0, tg_yyyyy_xxyyy_1, tg_yyyyy_xxyyz_0, tg_yyyyy_xxyyz_1, \
                                         tg_yyyyy_xxyz_1, tg_yyyyy_xxyzz_0, tg_yyyyy_xxyzz_1, tg_yyyyy_xxzz_1, \
                                         tg_yyyyy_xxzzz_0, tg_yyyyy_xxzzz_1, tg_yyyyy_xyyy_1, tg_yyyyy_xyyyy_0, \
                                         tg_yyyyy_xyyyy_1, tg_yyyyy_xyyyz_0, tg_yyyyy_xyyyz_1, tg_yyyyy_xyyz_1, \
                                         tg_yyyyy_xyyzz_0, tg_yyyyy_xyyzz_1, tg_yyyyy_xyzz_1, tg_yyyyy_xyzzz_0, \
                                         tg_yyyyy_xyzzz_1, tg_yyyyy_xzzz_1, tg_yyyyy_xzzzz_0, tg_yyyyy_xzzzz_1, \
                                         tg_yyyyy_yyyy_1, tg_yyyyy_yyyyy_0, tg_yyyyy_yyyyy_1, tg_yyyyy_yyyyz_0, \
                                         tg_yyyyy_yyyyz_1, tg_yyyyy_yyyz_1, tg_yyyyy_yyyzz_0, tg_yyyyy_yyyzz_1, \
                                         tg_yyyyy_yyzz_1, tg_yyyyy_yyzzz_0, tg_yyyyy_yyzzz_1, tg_yyyyy_yzzz_1, \
                                         tg_yyyyy_yzzzz_0, tg_yyyyy_yzzzz_1, tg_yyyyy_zzzz_1, tg_yyyyy_zzzzz_0, \
                                         tg_yyyyy_zzzzz_1, tg_yyyyz_xxxx_1, tg_yyyyz_xxxxx_0, tg_yyyyz_xxxxx_1, \
                                         tg_yyyyz_xxxxy_0, tg_yyyyz_xxxxy_1, tg_yyyyz_xxxxz_0, tg_yyyyz_xxxxz_1, \
                                         tg_yyyyz_xxxy_1, tg_yyyyz_xxxyy_0, tg_yyyyz_xxxyy_1, tg_yyyyz_xxxyz_0, \
                                         tg_yyyyz_xxxyz_1, tg_yyyyz_xxxz_1, tg_yyyyz_xxxzz_0, tg_yyyyz_xxxzz_1, \
                                         tg_yyyyz_xxyy_1, tg_yyyyz_xxyyy_0, tg_yyyyz_xxyyy_1, tg_yyyyz_xxyyz_0, \
                                         tg_yyyyz_xxyyz_1, tg_yyyyz_xxyz_1, tg_yyyyz_xxyzz_0, tg_yyyyz_xxyzz_1, \
                                         tg_yyyyz_xxzz_1, tg_yyyyz_xxzzz_0, tg_yyyyz_xxzzz_1, tg_yyyyz_xyyy_1, \
                                         tg_yyyyz_xyyyy_0, tg_yyyyz_xyyyy_1, tg_yyyyz_xyyyz_0, tg_yyyyz_xyyyz_1, \
                                         tg_yyyyz_xyyz_1, tg_yyyyz_xyyzz_0, tg_yyyyz_xyyzz_1, tg_yyyyz_xyzz_1, \
                                         tg_yyyyz_xyzzz_0, tg_yyyyz_xyzzz_1, tg_yyyyz_xzzz_1, tg_yyyyz_xzzzz_0, \
                                         tg_yyyyz_xzzzz_1, tg_yyyyz_yyyy_1, tg_yyyyz_yyyyy_0, tg_yyyyz_yyyyy_1, \
                                         tg_yyyyz_yyyyz_0, tg_yyyyz_yyyyz_1, tg_yyyyz_yyyz_1, tg_yyyyz_yyyzz_0, \
                                         tg_yyyyz_yyyzz_1, tg_yyyyz_yyzz_1, tg_yyyyz_yyzzz_0, tg_yyyyz_yyzzz_1, \
                                         tg_yyyyz_yzzz_1, tg_yyyyz_yzzzz_0, tg_yyyyz_yzzzz_1, tg_yyyyz_zzzz_1, \
                                         tg_yyyyz_zzzzz_0, tg_yyyyz_zzzzz_1, tg_yyyzz_xxxx_1, tg_yyyzz_xxxxx_0, \
                                         tg_yyyzz_xxxxx_1, tg_yyyzz_xxxxy_0, tg_yyyzz_xxxxy_1, tg_yyyzz_xxxxz_0, \
                                         tg_yyyzz_xxxxz_1, tg_yyyzz_xxxy_1, tg_yyyzz_xxxyy_0, tg_yyyzz_xxxyy_1, \
                                         tg_yyyzz_xxxyz_0, tg_yyyzz_xxxyz_1, tg_yyyzz_xxxz_1, tg_yyyzz_xxxzz_0, \
                                         tg_yyyzz_xxxzz_1, tg_yyyzz_xxyy_1, tg_yyyzz_xxyyy_0, tg_yyyzz_xxyyy_1, \
                                         tg_yyyzz_xxyyz_0, tg_yyyzz_xxyyz_1, tg_yyyzz_xxyz_1, tg_yyyzz_xxyzz_0, \
                                         tg_yyyzz_xxyzz_1, tg_yyyzz_xxzz_1, tg_yyyzz_xxzzz_0, tg_yyyzz_xxzzz_1, \
                                         tg_yyyzz_xyyy_1, tg_yyyzz_xyyyy_0, tg_yyyzz_xyyyy_1, tg_yyyzz_xyyyz_0, \
                                         tg_yyyzz_xyyyz_1, tg_yyyzz_xyyz_1, tg_yyyzz_xyyzz_0, tg_yyyzz_xyyzz_1, \
                                         tg_yyyzz_xyzz_1, tg_yyyzz_xyzzz_0, tg_yyyzz_xyzzz_1, tg_yyyzz_xzzz_1, \
                                         tg_yyyzz_xzzzz_0, tg_yyyzz_xzzzz_1, tg_yyyzz_yyyy_1, tg_yyyzz_yyyyy_0, \
                                         tg_yyyzz_yyyyy_1, tg_yyyzz_yyyyz_0, tg_yyyzz_yyyyz_1, tg_yyyzz_yyyz_1, \
                                         tg_yyyzz_yyyzz_0, tg_yyyzz_yyyzz_1, tg_yyyzz_yyzz_1, tg_yyyzz_yyzzz_0, \
                                         tg_yyyzz_yyzzz_1, tg_yyyzz_yzzz_1, tg_yyyzz_yzzzz_0, tg_yyyzz_yzzzz_1, \
                                         tg_yyyzz_zzzz_1, tg_yyyzz_zzzzz_0, tg_yyyzz_zzzzz_1, tg_yyzzz_xxxx_1, \
                                         tg_yyzzz_xxxxx_0, tg_yyzzz_xxxxx_1, tg_yyzzz_xxxxy_0, tg_yyzzz_xxxxy_1, \
                                         tg_yyzzz_xxxxz_0, tg_yyzzz_xxxxz_1, tg_yyzzz_xxxy_1, tg_yyzzz_xxxyy_0, \
                                         tg_yyzzz_xxxyy_1, tg_yyzzz_xxxyz_0, tg_yyzzz_xxxyz_1, tg_yyzzz_xxxz_1, \
                                         tg_yyzzz_xxxzz_0, tg_yyzzz_xxxzz_1, tg_yyzzz_xxyy_1, tg_yyzzz_xxyyy_0, \
                                         tg_yyzzz_xxyyy_1, tg_yyzzz_xxyyz_0, tg_yyzzz_xxyyz_1, tg_yyzzz_xxyz_1, \
                                         tg_yyzzz_xxyzz_0, tg_yyzzz_xxyzz_1, tg_yyzzz_xxzz_1, tg_yyzzz_xxzzz_0, \
                                         tg_yyzzz_xxzzz_1, tg_yyzzz_xyyy_1, tg_yyzzz_xyyyy_0, tg_yyzzz_xyyyy_1, \
                                         tg_yyzzz_xyyyz_0, tg_yyzzz_xyyyz_1, tg_yyzzz_xyyz_1, tg_yyzzz_xyyzz_0, \
                                         tg_yyzzz_xyyzz_1, tg_yyzzz_xyzz_1, tg_yyzzz_xyzzz_0, tg_yyzzz_xyzzz_1, \
                                         tg_yyzzz_xzzz_1, tg_yyzzz_yyyy_1, tg_yyzzz_yyyz_1, tg_yyzzz_yyzz_1, tg_yyzzz_yzzz_1, \
                                         tg_zzzz_xxxxx_0, tg_zzzz_xxxxx_1, tg_zzzz_xxxxy_0, tg_zzzz_xxxxy_1, tg_zzzz_xxxxz_0, \
                                         tg_zzzz_xxxxz_1, tg_zzzz_xxxyy_0, tg_zzzz_xxxyy_1, tg_zzzz_xxxyz_0, tg_zzzz_xxxyz_1, \
                                         tg_zzzz_xxxzz_0, tg_zzzz_xxxzz_1, tg_zzzz_xxyyy_0, tg_zzzz_xxyyy_1, tg_zzzz_xxyyz_0, \
                                         tg_zzzz_xxyyz_1, tg_zzzz_xxyzz_0, tg_zzzz_xxyzz_1, tg_zzzz_xxzzz_0, tg_zzzz_xxzzz_1, \
                                         tg_zzzz_xyyyy_0, tg_zzzz_xyyyy_1, tg_zzzz_xyyyz_0, tg_zzzz_xyyyz_1, tg_zzzz_xyyzz_0, \
                                         tg_zzzz_xyyzz_1, tg_zzzz_xyzzz_0, tg_zzzz_xyzzz_1, tg_zzzz_xzzzz_0, tg_zzzz_xzzzz_1, \
                                         tg_zzzz_yyyyy_0, tg_zzzz_yyyyy_1, tg_zzzz_yyyyz_0, tg_zzzz_yyyyz_1, tg_zzzz_yyyzz_0, \
                                         tg_zzzz_yyyzz_1, tg_zzzz_yyzzz_0, tg_zzzz_yyzzz_1, tg_zzzz_yzzzz_0, tg_zzzz_yzzzz_1, \
                                         tg_zzzz_zzzzz_0, tg_zzzz_zzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxzzzz_xxxxx_0[j] = pb_x * tg_xzzzz_xxxxx_0[j] + fr * tg_xzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxx_0[j] - tg_zzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzzz_xxxx_1[j];

                    tg_xxzzzz_xxxxy_0[j] = pb_x * tg_xzzzz_xxxxy_0[j] + fr * tg_xzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxy_0[j] - tg_zzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzz_xxxy_1[j];

                    tg_xxzzzz_xxxxz_0[j] = pb_x * tg_xzzzz_xxxxz_0[j] + fr * tg_xzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxz_0[j] - tg_zzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzz_xxxz_1[j];

                    tg_xxzzzz_xxxyy_0[j] = pb_x * tg_xzzzz_xxxyy_0[j] + fr * tg_xzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyy_0[j] - tg_zzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzz_xxyy_1[j];

                    tg_xxzzzz_xxxyz_0[j] = pb_x * tg_xzzzz_xxxyz_0[j] + fr * tg_xzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyz_0[j] - tg_zzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzz_xxyz_1[j];

                    tg_xxzzzz_xxxzz_0[j] = pb_x * tg_xzzzz_xxxzz_0[j] + fr * tg_xzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxzz_0[j] - tg_zzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzz_xxzz_1[j];

                    tg_xxzzzz_xxyyy_0[j] = pb_x * tg_xzzzz_xxyyy_0[j] + fr * tg_xzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyy_0[j] - tg_zzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xyyy_1[j];

                    tg_xxzzzz_xxyyz_0[j] = pb_x * tg_xzzzz_xxyyz_0[j] + fr * tg_xzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyz_0[j] - tg_zzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xyyz_1[j];

                    tg_xxzzzz_xxyzz_0[j] = pb_x * tg_xzzzz_xxyzz_0[j] + fr * tg_xzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyzz_0[j] - tg_zzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xyzz_1[j];

                    tg_xxzzzz_xxzzz_0[j] = pb_x * tg_xzzzz_xxzzz_0[j] + fr * tg_xzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxzzz_0[j] - tg_zzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzz_xzzz_1[j];

                    tg_xxzzzz_xyyyy_0[j] = pb_x * tg_xzzzz_xyyyy_0[j] + fr * tg_xzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyy_0[j] - tg_zzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yyyy_1[j];

                    tg_xxzzzz_xyyyz_0[j] = pb_x * tg_xzzzz_xyyyz_0[j] + fr * tg_xzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyz_0[j] - tg_zzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yyyz_1[j];

                    tg_xxzzzz_xyyzz_0[j] = pb_x * tg_xzzzz_xyyzz_0[j] + fr * tg_xzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyzz_0[j] - tg_zzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yyzz_1[j];

                    tg_xxzzzz_xyzzz_0[j] = pb_x * tg_xzzzz_xyzzz_0[j] + fr * tg_xzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyzzz_0[j] - tg_zzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_yzzz_1[j];

                    tg_xxzzzz_xzzzz_0[j] = pb_x * tg_xzzzz_xzzzz_0[j] + fr * tg_xzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xzzzz_0[j] - tg_zzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzz_zzzz_1[j];

                    tg_xxzzzz_yyyyy_0[j] = pb_x * tg_xzzzz_yyyyy_0[j] + fr * tg_xzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyy_0[j] - tg_zzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxzzzz_yyyyz_0[j] = pb_x * tg_xzzzz_yyyyz_0[j] + fr * tg_xzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyz_0[j] - tg_zzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxzzzz_yyyzz_0[j] = pb_x * tg_xzzzz_yyyzz_0[j] + fr * tg_xzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyzz_0[j] - tg_zzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxzzzz_yyzzz_0[j] = pb_x * tg_xzzzz_yyzzz_0[j] + fr * tg_xzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyzzz_0[j] - tg_zzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxzzzz_yzzzz_0[j] = pb_x * tg_xzzzz_yzzzz_0[j] + fr * tg_xzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yzzzz_0[j] - tg_zzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxzzzz_zzzzz_0[j] = pb_x * tg_xzzzz_zzzzz_0[j] + fr * tg_xzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_zzzzz_0[j] - tg_zzzz_zzzzz_1[j] * fl1_fza);

                    tg_xyyyyy_xxxxx_0[j] = pb_x * tg_yyyyy_xxxxx_0[j] + fr * tg_yyyyy_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyy_xxxx_1[j];

                    tg_xyyyyy_xxxxy_0[j] = pb_x * tg_yyyyy_xxxxy_0[j] + fr * tg_yyyyy_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyy_xxxy_1[j];

                    tg_xyyyyy_xxxxz_0[j] = pb_x * tg_yyyyy_xxxxz_0[j] + fr * tg_yyyyy_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyy_xxxz_1[j];

                    tg_xyyyyy_xxxyy_0[j] = pb_x * tg_yyyyy_xxxyy_0[j] + fr * tg_yyyyy_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxyy_1[j];

                    tg_xyyyyy_xxxyz_0[j] = pb_x * tg_yyyyy_xxxyz_0[j] + fr * tg_yyyyy_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxyz_1[j];

                    tg_xyyyyy_xxxzz_0[j] = pb_x * tg_yyyyy_xxxzz_0[j] + fr * tg_yyyyy_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyy_xxzz_1[j];

                    tg_xyyyyy_xxyyy_0[j] = pb_x * tg_yyyyy_xxyyy_0[j] + fr * tg_yyyyy_xxyyy_1[j] + fl1_fxn * tg_yyyyy_xyyy_1[j];

                    tg_xyyyyy_xxyyz_0[j] = pb_x * tg_yyyyy_xxyyz_0[j] + fr * tg_yyyyy_xxyyz_1[j] + fl1_fxn * tg_yyyyy_xyyz_1[j];

                    tg_xyyyyy_xxyzz_0[j] = pb_x * tg_yyyyy_xxyzz_0[j] + fr * tg_yyyyy_xxyzz_1[j] + fl1_fxn * tg_yyyyy_xyzz_1[j];

                    tg_xyyyyy_xxzzz_0[j] = pb_x * tg_yyyyy_xxzzz_0[j] + fr * tg_yyyyy_xxzzz_1[j] + fl1_fxn * tg_yyyyy_xzzz_1[j];

                    tg_xyyyyy_xyyyy_0[j] = pb_x * tg_yyyyy_xyyyy_0[j] + fr * tg_yyyyy_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyyy_1[j];

                    tg_xyyyyy_xyyyz_0[j] = pb_x * tg_yyyyy_xyyyz_0[j] + fr * tg_yyyyy_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyyz_1[j];

                    tg_xyyyyy_xyyzz_0[j] = pb_x * tg_yyyyy_xyyzz_0[j] + fr * tg_yyyyy_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yyzz_1[j];

                    tg_xyyyyy_xyzzz_0[j] = pb_x * tg_yyyyy_xyzzz_0[j] + fr * tg_yyyyy_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_yzzz_1[j];

                    tg_xyyyyy_xzzzz_0[j] = pb_x * tg_yyyyy_xzzzz_0[j] + fr * tg_yyyyy_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyy_zzzz_1[j];

                    tg_xyyyyy_yyyyy_0[j] = pb_x * tg_yyyyy_yyyyy_0[j] + fr * tg_yyyyy_yyyyy_1[j];

                    tg_xyyyyy_yyyyz_0[j] = pb_x * tg_yyyyy_yyyyz_0[j] + fr * tg_yyyyy_yyyyz_1[j];

                    tg_xyyyyy_yyyzz_0[j] = pb_x * tg_yyyyy_yyyzz_0[j] + fr * tg_yyyyy_yyyzz_1[j];

                    tg_xyyyyy_yyzzz_0[j] = pb_x * tg_yyyyy_yyzzz_0[j] + fr * tg_yyyyy_yyzzz_1[j];

                    tg_xyyyyy_yzzzz_0[j] = pb_x * tg_yyyyy_yzzzz_0[j] + fr * tg_yyyyy_yzzzz_1[j];

                    tg_xyyyyy_zzzzz_0[j] = pb_x * tg_yyyyy_zzzzz_0[j] + fr * tg_yyyyy_zzzzz_1[j];

                    tg_xyyyyz_xxxxx_0[j] = pb_x * tg_yyyyz_xxxxx_0[j] + fr * tg_yyyyz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyz_xxxx_1[j];

                    tg_xyyyyz_xxxxy_0[j] = pb_x * tg_yyyyz_xxxxy_0[j] + fr * tg_yyyyz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyz_xxxy_1[j];

                    tg_xyyyyz_xxxxz_0[j] = pb_x * tg_yyyyz_xxxxz_0[j] + fr * tg_yyyyz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyz_xxxz_1[j];

                    tg_xyyyyz_xxxyy_0[j] = pb_x * tg_yyyyz_xxxyy_0[j] + fr * tg_yyyyz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxyy_1[j];

                    tg_xyyyyz_xxxyz_0[j] = pb_x * tg_yyyyz_xxxyz_0[j] + fr * tg_yyyyz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxyz_1[j];

                    tg_xyyyyz_xxxzz_0[j] = pb_x * tg_yyyyz_xxxzz_0[j] + fr * tg_yyyyz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyz_xxzz_1[j];

                    tg_xyyyyz_xxyyy_0[j] = pb_x * tg_yyyyz_xxyyy_0[j] + fr * tg_yyyyz_xxyyy_1[j] + fl1_fxn * tg_yyyyz_xyyy_1[j];

                    tg_xyyyyz_xxyyz_0[j] = pb_x * tg_yyyyz_xxyyz_0[j] + fr * tg_yyyyz_xxyyz_1[j] + fl1_fxn * tg_yyyyz_xyyz_1[j];

                    tg_xyyyyz_xxyzz_0[j] = pb_x * tg_yyyyz_xxyzz_0[j] + fr * tg_yyyyz_xxyzz_1[j] + fl1_fxn * tg_yyyyz_xyzz_1[j];

                    tg_xyyyyz_xxzzz_0[j] = pb_x * tg_yyyyz_xxzzz_0[j] + fr * tg_yyyyz_xxzzz_1[j] + fl1_fxn * tg_yyyyz_xzzz_1[j];

                    tg_xyyyyz_xyyyy_0[j] = pb_x * tg_yyyyz_xyyyy_0[j] + fr * tg_yyyyz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyyy_1[j];

                    tg_xyyyyz_xyyyz_0[j] = pb_x * tg_yyyyz_xyyyz_0[j] + fr * tg_yyyyz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyyz_1[j];

                    tg_xyyyyz_xyyzz_0[j] = pb_x * tg_yyyyz_xyyzz_0[j] + fr * tg_yyyyz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yyzz_1[j];

                    tg_xyyyyz_xyzzz_0[j] = pb_x * tg_yyyyz_xyzzz_0[j] + fr * tg_yyyyz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_yzzz_1[j];

                    tg_xyyyyz_xzzzz_0[j] = pb_x * tg_yyyyz_xzzzz_0[j] + fr * tg_yyyyz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyz_zzzz_1[j];

                    tg_xyyyyz_yyyyy_0[j] = pb_x * tg_yyyyz_yyyyy_0[j] + fr * tg_yyyyz_yyyyy_1[j];

                    tg_xyyyyz_yyyyz_0[j] = pb_x * tg_yyyyz_yyyyz_0[j] + fr * tg_yyyyz_yyyyz_1[j];

                    tg_xyyyyz_yyyzz_0[j] = pb_x * tg_yyyyz_yyyzz_0[j] + fr * tg_yyyyz_yyyzz_1[j];

                    tg_xyyyyz_yyzzz_0[j] = pb_x * tg_yyyyz_yyzzz_0[j] + fr * tg_yyyyz_yyzzz_1[j];

                    tg_xyyyyz_yzzzz_0[j] = pb_x * tg_yyyyz_yzzzz_0[j] + fr * tg_yyyyz_yzzzz_1[j];

                    tg_xyyyyz_zzzzz_0[j] = pb_x * tg_yyyyz_zzzzz_0[j] + fr * tg_yyyyz_zzzzz_1[j];

                    tg_xyyyzz_xxxxx_0[j] = pb_x * tg_yyyzz_xxxxx_0[j] + fr * tg_yyyzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyzz_xxxx_1[j];

                    tg_xyyyzz_xxxxy_0[j] = pb_x * tg_yyyzz_xxxxy_0[j] + fr * tg_yyyzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyzz_xxxy_1[j];

                    tg_xyyyzz_xxxxz_0[j] = pb_x * tg_yyyzz_xxxxz_0[j] + fr * tg_yyyzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyzz_xxxz_1[j];

                    tg_xyyyzz_xxxyy_0[j] = pb_x * tg_yyyzz_xxxyy_0[j] + fr * tg_yyyzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxyy_1[j];

                    tg_xyyyzz_xxxyz_0[j] = pb_x * tg_yyyzz_xxxyz_0[j] + fr * tg_yyyzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxyz_1[j];

                    tg_xyyyzz_xxxzz_0[j] = pb_x * tg_yyyzz_xxxzz_0[j] + fr * tg_yyyzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyzz_xxzz_1[j];

                    tg_xyyyzz_xxyyy_0[j] = pb_x * tg_yyyzz_xxyyy_0[j] + fr * tg_yyyzz_xxyyy_1[j] + fl1_fxn * tg_yyyzz_xyyy_1[j];

                    tg_xyyyzz_xxyyz_0[j] = pb_x * tg_yyyzz_xxyyz_0[j] + fr * tg_yyyzz_xxyyz_1[j] + fl1_fxn * tg_yyyzz_xyyz_1[j];

                    tg_xyyyzz_xxyzz_0[j] = pb_x * tg_yyyzz_xxyzz_0[j] + fr * tg_yyyzz_xxyzz_1[j] + fl1_fxn * tg_yyyzz_xyzz_1[j];

                    tg_xyyyzz_xxzzz_0[j] = pb_x * tg_yyyzz_xxzzz_0[j] + fr * tg_yyyzz_xxzzz_1[j] + fl1_fxn * tg_yyyzz_xzzz_1[j];

                    tg_xyyyzz_xyyyy_0[j] = pb_x * tg_yyyzz_xyyyy_0[j] + fr * tg_yyyzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyyy_1[j];

                    tg_xyyyzz_xyyyz_0[j] = pb_x * tg_yyyzz_xyyyz_0[j] + fr * tg_yyyzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyyz_1[j];

                    tg_xyyyzz_xyyzz_0[j] = pb_x * tg_yyyzz_xyyzz_0[j] + fr * tg_yyyzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yyzz_1[j];

                    tg_xyyyzz_xyzzz_0[j] = pb_x * tg_yyyzz_xyzzz_0[j] + fr * tg_yyyzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_yzzz_1[j];

                    tg_xyyyzz_xzzzz_0[j] = pb_x * tg_yyyzz_xzzzz_0[j] + fr * tg_yyyzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzz_zzzz_1[j];

                    tg_xyyyzz_yyyyy_0[j] = pb_x * tg_yyyzz_yyyyy_0[j] + fr * tg_yyyzz_yyyyy_1[j];

                    tg_xyyyzz_yyyyz_0[j] = pb_x * tg_yyyzz_yyyyz_0[j] + fr * tg_yyyzz_yyyyz_1[j];

                    tg_xyyyzz_yyyzz_0[j] = pb_x * tg_yyyzz_yyyzz_0[j] + fr * tg_yyyzz_yyyzz_1[j];

                    tg_xyyyzz_yyzzz_0[j] = pb_x * tg_yyyzz_yyzzz_0[j] + fr * tg_yyyzz_yyzzz_1[j];

                    tg_xyyyzz_yzzzz_0[j] = pb_x * tg_yyyzz_yzzzz_0[j] + fr * tg_yyyzz_yzzzz_1[j];

                    tg_xyyyzz_zzzzz_0[j] = pb_x * tg_yyyzz_zzzzz_0[j] + fr * tg_yyyzz_zzzzz_1[j];

                    tg_xyyzzz_xxxxx_0[j] = pb_x * tg_yyzzz_xxxxx_0[j] + fr * tg_yyzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyzzz_xxxx_1[j];

                    tg_xyyzzz_xxxxy_0[j] = pb_x * tg_yyzzz_xxxxy_0[j] + fr * tg_yyzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyzzz_xxxy_1[j];

                    tg_xyyzzz_xxxxz_0[j] = pb_x * tg_yyzzz_xxxxz_0[j] + fr * tg_yyzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyzzz_xxxz_1[j];

                    tg_xyyzzz_xxxyy_0[j] = pb_x * tg_yyzzz_xxxyy_0[j] + fr * tg_yyzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxyy_1[j];

                    tg_xyyzzz_xxxyz_0[j] = pb_x * tg_yyzzz_xxxyz_0[j] + fr * tg_yyzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxyz_1[j];

                    tg_xyyzzz_xxxzz_0[j] = pb_x * tg_yyzzz_xxxzz_0[j] + fr * tg_yyzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyzzz_xxzz_1[j];

                    tg_xyyzzz_xxyyy_0[j] = pb_x * tg_yyzzz_xxyyy_0[j] + fr * tg_yyzzz_xxyyy_1[j] + fl1_fxn * tg_yyzzz_xyyy_1[j];

                    tg_xyyzzz_xxyyz_0[j] = pb_x * tg_yyzzz_xxyyz_0[j] + fr * tg_yyzzz_xxyyz_1[j] + fl1_fxn * tg_yyzzz_xyyz_1[j];

                    tg_xyyzzz_xxyzz_0[j] = pb_x * tg_yyzzz_xxyzz_0[j] + fr * tg_yyzzz_xxyzz_1[j] + fl1_fxn * tg_yyzzz_xyzz_1[j];

                    tg_xyyzzz_xxzzz_0[j] = pb_x * tg_yyzzz_xxzzz_0[j] + fr * tg_yyzzz_xxzzz_1[j] + fl1_fxn * tg_yyzzz_xzzz_1[j];

                    tg_xyyzzz_xyyyy_0[j] = pb_x * tg_yyzzz_xyyyy_0[j] + fr * tg_yyzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyyy_1[j];

                    tg_xyyzzz_xyyyz_0[j] = pb_x * tg_yyzzz_xyyyz_0[j] + fr * tg_yyzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyyz_1[j];

                    tg_xyyzzz_xyyzz_0[j] = pb_x * tg_yyzzz_xyyzz_0[j] + fr * tg_yyzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yyzz_1[j];

                    tg_xyyzzz_xyzzz_0[j] = pb_x * tg_yyzzz_xyzzz_0[j] + fr * tg_yyzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_yzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISH_392_490(      CMemBlock2D<double>* primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyyy_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 210); 

                auto tg_yyyy_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 211); 

                auto tg_yyyy_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 212); 

                auto tg_yyyy_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 213); 

                auto tg_yyyy_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 214); 

                auto tg_yyyy_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 215); 

                auto tg_yyyy_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 216); 

                auto tg_yyyy_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 217); 

                auto tg_yyyy_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 218); 

                auto tg_yyyy_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 219); 

                auto tg_yyyy_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 220); 

                auto tg_yyyy_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 221); 

                auto tg_yyyy_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 222); 

                auto tg_yyyy_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 223); 

                auto tg_yyyy_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 224); 

                auto tg_yyyy_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 225); 

                auto tg_yyyy_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 226); 

                auto tg_yyyy_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 227); 

                auto tg_yyyy_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 228); 

                auto tg_yyyy_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 229); 

                auto tg_yyyy_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 230); 

                auto tg_yyyz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 231); 

                auto tg_yyyz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 232); 

                auto tg_yyyz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 233); 

                auto tg_yyyz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 234); 

                auto tg_yyyz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 235); 

                auto tg_yyyz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 236); 

                auto tg_yyyz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 237); 

                auto tg_yyyz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 238); 

                auto tg_yyyz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 239); 

                auto tg_yyyz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 240); 

                auto tg_yyyz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 241); 

                auto tg_yyyz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 242); 

                auto tg_yyyz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 243); 

                auto tg_yyyz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 244); 

                auto tg_yyyz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 245); 

                auto tg_yyyz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 246); 

                auto tg_yyyz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 247); 

                auto tg_yyyz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 248); 

                auto tg_yyyz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 249); 

                auto tg_yyyz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 250); 

                auto tg_yyyz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 251); 

                auto tg_yyzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 252); 

                auto tg_yyzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 253); 

                auto tg_yyzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 254); 

                auto tg_yyzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 255); 

                auto tg_yyzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 256); 

                auto tg_yyzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 257); 

                auto tg_yyzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 258); 

                auto tg_yyyy_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 210); 

                auto tg_yyyy_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 211); 

                auto tg_yyyy_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 212); 

                auto tg_yyyy_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 213); 

                auto tg_yyyy_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 214); 

                auto tg_yyyy_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 215); 

                auto tg_yyyy_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 216); 

                auto tg_yyyy_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 217); 

                auto tg_yyyy_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 218); 

                auto tg_yyyy_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 219); 

                auto tg_yyyy_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 220); 

                auto tg_yyyy_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 221); 

                auto tg_yyyy_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 222); 

                auto tg_yyyy_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 223); 

                auto tg_yyyy_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 224); 

                auto tg_yyyy_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 225); 

                auto tg_yyyy_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 226); 

                auto tg_yyyy_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 227); 

                auto tg_yyyy_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 228); 

                auto tg_yyyy_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 229); 

                auto tg_yyyy_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 230); 

                auto tg_yyyz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 231); 

                auto tg_yyyz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 232); 

                auto tg_yyyz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 233); 

                auto tg_yyyz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 234); 

                auto tg_yyyz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 235); 

                auto tg_yyyz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 236); 

                auto tg_yyyz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 237); 

                auto tg_yyyz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 238); 

                auto tg_yyyz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 239); 

                auto tg_yyyz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 240); 

                auto tg_yyyz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 241); 

                auto tg_yyyz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 242); 

                auto tg_yyyz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 243); 

                auto tg_yyyz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 244); 

                auto tg_yyyz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 245); 

                auto tg_yyyz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 246); 

                auto tg_yyyz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 247); 

                auto tg_yyyz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 248); 

                auto tg_yyyz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 249); 

                auto tg_yyyz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 250); 

                auto tg_yyyz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 251); 

                auto tg_yyzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 252); 

                auto tg_yyzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 253); 

                auto tg_yyzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 254); 

                auto tg_yyzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 255); 

                auto tg_yyzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 256); 

                auto tg_yyzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 257); 

                auto tg_yyzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 258); 

                auto tg_yyyyy_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 225); 

                auto tg_yyyyy_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 226); 

                auto tg_yyyyy_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 227); 

                auto tg_yyyyy_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 228); 

                auto tg_yyyyy_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 229); 

                auto tg_yyyyy_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 230); 

                auto tg_yyyyy_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 231); 

                auto tg_yyyyy_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 232); 

                auto tg_yyyyy_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 233); 

                auto tg_yyyyy_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 234); 

                auto tg_yyyyy_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 235); 

                auto tg_yyyyy_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 236); 

                auto tg_yyyyy_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 237); 

                auto tg_yyyyy_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 238); 

                auto tg_yyyyy_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 239); 

                auto tg_yyyyz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 240); 

                auto tg_yyyyz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 241); 

                auto tg_yyyyz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 242); 

                auto tg_yyyyz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 243); 

                auto tg_yyyyz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 244); 

                auto tg_yyyyz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 245); 

                auto tg_yyyyz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 246); 

                auto tg_yyyyz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 247); 

                auto tg_yyyyz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 248); 

                auto tg_yyyyz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 249); 

                auto tg_yyyyz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 250); 

                auto tg_yyyyz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 251); 

                auto tg_yyyyz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 252); 

                auto tg_yyyyz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 253); 

                auto tg_yyyyz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 254); 

                auto tg_yyyzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 255); 

                auto tg_yyyzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 256); 

                auto tg_yyyzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 257); 

                auto tg_yyyzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 258); 

                auto tg_yyzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 284); 

                auto tg_yzzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 285); 

                auto tg_yzzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 286); 

                auto tg_yzzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 287); 

                auto tg_yzzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 288); 

                auto tg_yzzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 289); 

                auto tg_yzzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 290); 

                auto tg_yzzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 291); 

                auto tg_yzzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 292); 

                auto tg_yzzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 293); 

                auto tg_yzzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 294); 

                auto tg_yzzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 295); 

                auto tg_yzzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 296); 

                auto tg_yzzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 297); 

                auto tg_yzzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 298); 

                auto tg_yzzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 299); 

                auto tg_zzzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 300); 

                auto tg_zzzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 301); 

                auto tg_zzzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 302); 

                auto tg_zzzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 303); 

                auto tg_zzzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 304); 

                auto tg_zzzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 305); 

                auto tg_zzzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 306); 

                auto tg_zzzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 307); 

                auto tg_zzzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 308); 

                auto tg_zzzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 309); 

                auto tg_zzzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 310); 

                auto tg_zzzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 311); 

                auto tg_zzzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 312); 

                auto tg_zzzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 313); 

                auto tg_zzzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 314); 

                // set up pointers to integrals

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

                // Batch of Integrals (392,490)

                #pragma omp simd aligned(fxn, fza, tg_xyyzzz_xzzzz_0, tg_xyyzzz_yyyyy_0, tg_xyyzzz_yyyyz_0, \
                                         tg_xyyzzz_yyyzz_0, tg_xyyzzz_yyzzz_0, tg_xyyzzz_yzzzz_0, tg_xyyzzz_zzzzz_0, \
                                         tg_xyzzzz_xxxxx_0, tg_xyzzzz_xxxxy_0, tg_xyzzzz_xxxxz_0, tg_xyzzzz_xxxyy_0, \
                                         tg_xyzzzz_xxxyz_0, tg_xyzzzz_xxxzz_0, tg_xyzzzz_xxyyy_0, tg_xyzzzz_xxyyz_0, \
                                         tg_xyzzzz_xxyzz_0, tg_xyzzzz_xxzzz_0, tg_xyzzzz_xyyyy_0, tg_xyzzzz_xyyyz_0, \
                                         tg_xyzzzz_xyyzz_0, tg_xyzzzz_xyzzz_0, tg_xyzzzz_xzzzz_0, tg_xyzzzz_yyyyy_0, \
                                         tg_xyzzzz_yyyyz_0, tg_xyzzzz_yyyzz_0, tg_xyzzzz_yyzzz_0, tg_xyzzzz_yzzzz_0, \
                                         tg_xyzzzz_zzzzz_0, tg_xzzzzz_xxxxx_0, tg_xzzzzz_xxxxy_0, tg_xzzzzz_xxxxz_0, \
                                         tg_xzzzzz_xxxyy_0, tg_xzzzzz_xxxyz_0, tg_xzzzzz_xxxzz_0, tg_xzzzzz_xxyyy_0, \
                                         tg_xzzzzz_xxyyz_0, tg_xzzzzz_xxyzz_0, tg_xzzzzz_xxzzz_0, tg_xzzzzz_xyyyy_0, \
                                         tg_xzzzzz_xyyyz_0, tg_xzzzzz_xyyzz_0, tg_xzzzzz_xyzzz_0, tg_xzzzzz_xzzzz_0, \
                                         tg_xzzzzz_yyyyy_0, tg_xzzzzz_yyyyz_0, tg_xzzzzz_yyyzz_0, tg_xzzzzz_yyzzz_0, \
                                         tg_xzzzzz_yzzzz_0, tg_xzzzzz_zzzzz_0, tg_yyyy_xxxxx_0, tg_yyyy_xxxxx_1, \
                                         tg_yyyy_xxxxy_0, tg_yyyy_xxxxy_1, tg_yyyy_xxxxz_0, tg_yyyy_xxxxz_1, tg_yyyy_xxxyy_0, \
                                         tg_yyyy_xxxyy_1, tg_yyyy_xxxyz_0, tg_yyyy_xxxyz_1, tg_yyyy_xxxzz_0, tg_yyyy_xxxzz_1, \
                                         tg_yyyy_xxyyy_0, tg_yyyy_xxyyy_1, tg_yyyy_xxyyz_0, tg_yyyy_xxyyz_1, tg_yyyy_xxyzz_0, \
                                         tg_yyyy_xxyzz_1, tg_yyyy_xxzzz_0, tg_yyyy_xxzzz_1, tg_yyyy_xyyyy_0, tg_yyyy_xyyyy_1, \
                                         tg_yyyy_xyyyz_0, tg_yyyy_xyyyz_1, tg_yyyy_xyyzz_0, tg_yyyy_xyyzz_1, tg_yyyy_xyzzz_0, \
                                         tg_yyyy_xyzzz_1, tg_yyyy_xzzzz_0, tg_yyyy_xzzzz_1, tg_yyyy_yyyyy_0, tg_yyyy_yyyyy_1, \
                                         tg_yyyy_yyyyz_0, tg_yyyy_yyyyz_1, tg_yyyy_yyyzz_0, tg_yyyy_yyyzz_1, tg_yyyy_yyzzz_0, \
                                         tg_yyyy_yyzzz_1, tg_yyyy_yzzzz_0, tg_yyyy_yzzzz_1, tg_yyyy_zzzzz_0, tg_yyyy_zzzzz_1, \
                                         tg_yyyyy_xxxx_1, tg_yyyyy_xxxxx_0, tg_yyyyy_xxxxx_1, tg_yyyyy_xxxxy_0, \
                                         tg_yyyyy_xxxxy_1, tg_yyyyy_xxxxz_0, tg_yyyyy_xxxxz_1, tg_yyyyy_xxxy_1, \
                                         tg_yyyyy_xxxyy_0, tg_yyyyy_xxxyy_1, tg_yyyyy_xxxyz_0, tg_yyyyy_xxxyz_1, \
                                         tg_yyyyy_xxxz_1, tg_yyyyy_xxxzz_0, tg_yyyyy_xxxzz_1, tg_yyyyy_xxyy_1, \
                                         tg_yyyyy_xxyyy_0, tg_yyyyy_xxyyy_1, tg_yyyyy_xxyyz_0, tg_yyyyy_xxyyz_1, \
                                         tg_yyyyy_xxyz_1, tg_yyyyy_xxyzz_0, tg_yyyyy_xxyzz_1, tg_yyyyy_xxzz_1, \
                                         tg_yyyyy_xxzzz_0, tg_yyyyy_xxzzz_1, tg_yyyyy_xyyy_1, tg_yyyyy_xyyyy_0, \
                                         tg_yyyyy_xyyyy_1, tg_yyyyy_xyyyz_0, tg_yyyyy_xyyyz_1, tg_yyyyy_xyyz_1, \
                                         tg_yyyyy_xyyzz_0, tg_yyyyy_xyyzz_1, tg_yyyyy_xyzz_1, tg_yyyyy_xyzzz_0, \
                                         tg_yyyyy_xyzzz_1, tg_yyyyy_xzzz_1, tg_yyyyy_xzzzz_0, tg_yyyyy_xzzzz_1, \
                                         tg_yyyyy_yyyy_1, tg_yyyyy_yyyyy_0, tg_yyyyy_yyyyy_1, tg_yyyyy_yyyyz_0, \
                                         tg_yyyyy_yyyyz_1, tg_yyyyy_yyyz_1, tg_yyyyy_yyyzz_0, tg_yyyyy_yyyzz_1, \
                                         tg_yyyyy_yyzz_1, tg_yyyyy_yyzzz_0, tg_yyyyy_yyzzz_1, tg_yyyyy_yzzz_1, \
                                         tg_yyyyy_yzzzz_0, tg_yyyyy_yzzzz_1, tg_yyyyy_zzzz_1, tg_yyyyy_zzzzz_0, \
                                         tg_yyyyy_zzzzz_1, tg_yyyyyy_xxxxx_0, tg_yyyyyy_xxxxy_0, tg_yyyyyy_xxxxz_0, \
                                         tg_yyyyyy_xxxyy_0, tg_yyyyyy_xxxyz_0, tg_yyyyyy_xxxzz_0, tg_yyyyyy_xxyyy_0, \
                                         tg_yyyyyy_xxyyz_0, tg_yyyyyy_xxyzz_0, tg_yyyyyy_xxzzz_0, tg_yyyyyy_xyyyy_0, \
                                         tg_yyyyyy_xyyyz_0, tg_yyyyyy_xyyzz_0, tg_yyyyyy_xyzzz_0, tg_yyyyyy_xzzzz_0, \
                                         tg_yyyyyy_yyyyy_0, tg_yyyyyy_yyyyz_0, tg_yyyyyy_yyyzz_0, tg_yyyyyy_yyzzz_0, \
                                         tg_yyyyyy_yzzzz_0, tg_yyyyyy_zzzzz_0, tg_yyyyyz_xxxxx_0, tg_yyyyyz_xxxxy_0, \
                                         tg_yyyyyz_xxxxz_0, tg_yyyyyz_xxxyy_0, tg_yyyyyz_xxxyz_0, tg_yyyyyz_xxxzz_0, \
                                         tg_yyyyyz_xxyyy_0, tg_yyyyyz_xxyyz_0, tg_yyyyyz_xxyzz_0, tg_yyyyyz_xxzzz_0, \
                                         tg_yyyyyz_xyyyy_0, tg_yyyyyz_xyyyz_0, tg_yyyyyz_xyyzz_0, tg_yyyyyz_xyzzz_0, \
                                         tg_yyyyyz_xzzzz_0, tg_yyyyyz_yyyyy_0, tg_yyyyyz_yyyyz_0, tg_yyyyyz_yyyzz_0, \
                                         tg_yyyyyz_yyzzz_0, tg_yyyyyz_yzzzz_0, tg_yyyyyz_zzzzz_0, tg_yyyyz_xxxx_1, \
                                         tg_yyyyz_xxxxx_0, tg_yyyyz_xxxxx_1, tg_yyyyz_xxxxy_0, tg_yyyyz_xxxxy_1, \
                                         tg_yyyyz_xxxxz_0, tg_yyyyz_xxxxz_1, tg_yyyyz_xxxy_1, tg_yyyyz_xxxyy_0, \
                                         tg_yyyyz_xxxyy_1, tg_yyyyz_xxxyz_0, tg_yyyyz_xxxyz_1, tg_yyyyz_xxxz_1, \
                                         tg_yyyyz_xxxzz_0, tg_yyyyz_xxxzz_1, tg_yyyyz_xxyy_1, tg_yyyyz_xxyyy_0, \
                                         tg_yyyyz_xxyyy_1, tg_yyyyz_xxyyz_0, tg_yyyyz_xxyyz_1, tg_yyyyz_xxyz_1, \
                                         tg_yyyyz_xxyzz_0, tg_yyyyz_xxyzz_1, tg_yyyyz_xxzz_1, tg_yyyyz_xxzzz_0, \
                                         tg_yyyyz_xxzzz_1, tg_yyyyz_xyyy_1, tg_yyyyz_xyyyy_0, tg_yyyyz_xyyyy_1, \
                                         tg_yyyyz_xyyyz_0, tg_yyyyz_xyyyz_1, tg_yyyyz_xyyz_1, tg_yyyyz_xyyzz_0, \
                                         tg_yyyyz_xyyzz_1, tg_yyyyz_xyzz_1, tg_yyyyz_xyzzz_0, tg_yyyyz_xyzzz_1, \
                                         tg_yyyyz_xzzz_1, tg_yyyyz_xzzzz_0, tg_yyyyz_xzzzz_1, tg_yyyyz_yyyy_1, \
                                         tg_yyyyz_yyyyy_0, tg_yyyyz_yyyyy_1, tg_yyyyz_yyyyz_0, tg_yyyyz_yyyyz_1, \
                                         tg_yyyyz_yyyz_1, tg_yyyyz_yyyzz_0, tg_yyyyz_yyyzz_1, tg_yyyyz_yyzz_1, \
                                         tg_yyyyz_yyzzz_0, tg_yyyyz_yyzzz_1, tg_yyyyz_yzzz_1, tg_yyyyz_yzzzz_0, \
                                         tg_yyyyz_yzzzz_1, tg_yyyyz_zzzz_1, tg_yyyyz_zzzzz_0, tg_yyyyz_zzzzz_1, \
                                         tg_yyyyzz_xxxxx_0, tg_yyyyzz_xxxxy_0, tg_yyyyzz_xxxxz_0, tg_yyyyzz_xxxyy_0, \
                                         tg_yyyyzz_xxxyz_0, tg_yyyyzz_xxxzz_0, tg_yyyyzz_xxyyy_0, tg_yyyz_xxxxx_0, \
                                         tg_yyyz_xxxxx_1, tg_yyyz_xxxxy_0, tg_yyyz_xxxxy_1, tg_yyyz_xxxxz_0, tg_yyyz_xxxxz_1, \
                                         tg_yyyz_xxxyy_0, tg_yyyz_xxxyy_1, tg_yyyz_xxxyz_0, tg_yyyz_xxxyz_1, tg_yyyz_xxxzz_0, \
                                         tg_yyyz_xxxzz_1, tg_yyyz_xxyyy_0, tg_yyyz_xxyyy_1, tg_yyyz_xxyyz_0, tg_yyyz_xxyyz_1, \
                                         tg_yyyz_xxyzz_0, tg_yyyz_xxyzz_1, tg_yyyz_xxzzz_0, tg_yyyz_xxzzz_1, tg_yyyz_xyyyy_0, \
                                         tg_yyyz_xyyyy_1, tg_yyyz_xyyyz_0, tg_yyyz_xyyyz_1, tg_yyyz_xyyzz_0, tg_yyyz_xyyzz_1, \
                                         tg_yyyz_xyzzz_0, tg_yyyz_xyzzz_1, tg_yyyz_xzzzz_0, tg_yyyz_xzzzz_1, tg_yyyz_yyyyy_0, \
                                         tg_yyyz_yyyyy_1, tg_yyyz_yyyyz_0, tg_yyyz_yyyyz_1, tg_yyyz_yyyzz_0, tg_yyyz_yyyzz_1, \
                                         tg_yyyz_yyzzz_0, tg_yyyz_yyzzz_1, tg_yyyz_yzzzz_0, tg_yyyz_yzzzz_1, tg_yyyz_zzzzz_0, \
                                         tg_yyyz_zzzzz_1, tg_yyyzz_xxxx_1, tg_yyyzz_xxxxx_0, tg_yyyzz_xxxxx_1, \
                                         tg_yyyzz_xxxxy_0, tg_yyyzz_xxxxy_1, tg_yyyzz_xxxxz_0, tg_yyyzz_xxxxz_1, \
                                         tg_yyyzz_xxxy_1, tg_yyyzz_xxxyy_0, tg_yyyzz_xxxyy_1, tg_yyyzz_xxxyz_0, \
                                         tg_yyyzz_xxxyz_1, tg_yyyzz_xxxz_1, tg_yyyzz_xxxzz_0, tg_yyyzz_xxxzz_1, \
                                         tg_yyyzz_xxyy_1, tg_yyyzz_xxyyy_0, tg_yyyzz_xxyyy_1, tg_yyzz_xxxxx_0, \
                                         tg_yyzz_xxxxx_1, tg_yyzz_xxxxy_0, tg_yyzz_xxxxy_1, tg_yyzz_xxxxz_0, tg_yyzz_xxxxz_1, \
                                         tg_yyzz_xxxyy_0, tg_yyzz_xxxyy_1, tg_yyzz_xxxyz_0, tg_yyzz_xxxyz_1, tg_yyzz_xxxzz_0, \
                                         tg_yyzz_xxxzz_1, tg_yyzz_xxyyy_0, tg_yyzz_xxyyy_1, tg_yyzzz_xzzzz_0, \
                                         tg_yyzzz_xzzzz_1, tg_yyzzz_yyyyy_0, tg_yyzzz_yyyyy_1, tg_yyzzz_yyyyz_0, \
                                         tg_yyzzz_yyyyz_1, tg_yyzzz_yyyzz_0, tg_yyzzz_yyyzz_1, tg_yyzzz_yyzzz_0, \
                                         tg_yyzzz_yyzzz_1, tg_yyzzz_yzzzz_0, tg_yyzzz_yzzzz_1, tg_yyzzz_zzzz_1, \
                                         tg_yyzzz_zzzzz_0, tg_yyzzz_zzzzz_1, tg_yzzzz_xxxx_1, tg_yzzzz_xxxxx_0, \
                                         tg_yzzzz_xxxxx_1, tg_yzzzz_xxxxy_0, tg_yzzzz_xxxxy_1, tg_yzzzz_xxxxz_0, \
                                         tg_yzzzz_xxxxz_1, tg_yzzzz_xxxy_1, tg_yzzzz_xxxyy_0, tg_yzzzz_xxxyy_1, \
                                         tg_yzzzz_xxxyz_0, tg_yzzzz_xxxyz_1, tg_yzzzz_xxxz_1, tg_yzzzz_xxxzz_0, \
                                         tg_yzzzz_xxxzz_1, tg_yzzzz_xxyy_1, tg_yzzzz_xxyyy_0, tg_yzzzz_xxyyy_1, \
                                         tg_yzzzz_xxyyz_0, tg_yzzzz_xxyyz_1, tg_yzzzz_xxyz_1, tg_yzzzz_xxyzz_0, \
                                         tg_yzzzz_xxyzz_1, tg_yzzzz_xxzz_1, tg_yzzzz_xxzzz_0, tg_yzzzz_xxzzz_1, \
                                         tg_yzzzz_xyyy_1, tg_yzzzz_xyyyy_0, tg_yzzzz_xyyyy_1, tg_yzzzz_xyyyz_0, \
                                         tg_yzzzz_xyyyz_1, tg_yzzzz_xyyz_1, tg_yzzzz_xyyzz_0, tg_yzzzz_xyyzz_1, \
                                         tg_yzzzz_xyzz_1, tg_yzzzz_xyzzz_0, tg_yzzzz_xyzzz_1, tg_yzzzz_xzzz_1, \
                                         tg_yzzzz_xzzzz_0, tg_yzzzz_xzzzz_1, tg_yzzzz_yyyy_1, tg_yzzzz_yyyyy_0, \
                                         tg_yzzzz_yyyyy_1, tg_yzzzz_yyyyz_0, tg_yzzzz_yyyyz_1, tg_yzzzz_yyyz_1, \
                                         tg_yzzzz_yyyzz_0, tg_yzzzz_yyyzz_1, tg_yzzzz_yyzz_1, tg_yzzzz_yyzzz_0, \
                                         tg_yzzzz_yyzzz_1, tg_yzzzz_yzzz_1, tg_yzzzz_yzzzz_0, tg_yzzzz_yzzzz_1, \
                                         tg_yzzzz_zzzz_1, tg_yzzzz_zzzzz_0, tg_yzzzz_zzzzz_1, tg_zzzzz_xxxx_1, \
                                         tg_zzzzz_xxxxx_0, tg_zzzzz_xxxxx_1, tg_zzzzz_xxxxy_0, tg_zzzzz_xxxxy_1, \
                                         tg_zzzzz_xxxxz_0, tg_zzzzz_xxxxz_1, tg_zzzzz_xxxy_1, tg_zzzzz_xxxyy_0, \
                                         tg_zzzzz_xxxyy_1, tg_zzzzz_xxxyz_0, tg_zzzzz_xxxyz_1, tg_zzzzz_xxxz_1, \
                                         tg_zzzzz_xxxzz_0, tg_zzzzz_xxxzz_1, tg_zzzzz_xxyy_1, tg_zzzzz_xxyyy_0, \
                                         tg_zzzzz_xxyyy_1, tg_zzzzz_xxyyz_0, tg_zzzzz_xxyyz_1, tg_zzzzz_xxyz_1, \
                                         tg_zzzzz_xxyzz_0, tg_zzzzz_xxyzz_1, tg_zzzzz_xxzz_1, tg_zzzzz_xxzzz_0, \
                                         tg_zzzzz_xxzzz_1, tg_zzzzz_xyyy_1, tg_zzzzz_xyyyy_0, tg_zzzzz_xyyyy_1, \
                                         tg_zzzzz_xyyyz_0, tg_zzzzz_xyyyz_1, tg_zzzzz_xyyz_1, tg_zzzzz_xyyzz_0, \
                                         tg_zzzzz_xyyzz_1, tg_zzzzz_xyzz_1, tg_zzzzz_xyzzz_0, tg_zzzzz_xyzzz_1, \
                                         tg_zzzzz_xzzz_1, tg_zzzzz_xzzzz_0, tg_zzzzz_xzzzz_1, tg_zzzzz_yyyy_1, \
                                         tg_zzzzz_yyyyy_0, tg_zzzzz_yyyyy_1, tg_zzzzz_yyyyz_0, tg_zzzzz_yyyyz_1, \
                                         tg_zzzzz_yyyz_1, tg_zzzzz_yyyzz_0, tg_zzzzz_yyyzz_1, tg_zzzzz_yyzz_1, \
                                         tg_zzzzz_yyzzz_0, tg_zzzzz_yyzzz_1, tg_zzzzz_yzzz_1, tg_zzzzz_yzzzz_0, \
                                         tg_zzzzz_yzzzz_1, tg_zzzzz_zzzz_1, tg_zzzzz_zzzzz_0, tg_zzzzz_zzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyyzzz_xzzzz_0[j] = pb_x * tg_yyzzz_xzzzz_0[j] + fr * tg_yyzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzz_zzzz_1[j];

                    tg_xyyzzz_yyyyy_0[j] = pb_x * tg_yyzzz_yyyyy_0[j] + fr * tg_yyzzz_yyyyy_1[j];

                    tg_xyyzzz_yyyyz_0[j] = pb_x * tg_yyzzz_yyyyz_0[j] + fr * tg_yyzzz_yyyyz_1[j];

                    tg_xyyzzz_yyyzz_0[j] = pb_x * tg_yyzzz_yyyzz_0[j] + fr * tg_yyzzz_yyyzz_1[j];

                    tg_xyyzzz_yyzzz_0[j] = pb_x * tg_yyzzz_yyzzz_0[j] + fr * tg_yyzzz_yyzzz_1[j];

                    tg_xyyzzz_yzzzz_0[j] = pb_x * tg_yyzzz_yzzzz_0[j] + fr * tg_yyzzz_yzzzz_1[j];

                    tg_xyyzzz_zzzzz_0[j] = pb_x * tg_yyzzz_zzzzz_0[j] + fr * tg_yyzzz_zzzzz_1[j];

                    tg_xyzzzz_xxxxx_0[j] = pb_x * tg_yzzzz_xxxxx_0[j] + fr * tg_yzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yzzzz_xxxx_1[j];

                    tg_xyzzzz_xxxxy_0[j] = pb_x * tg_yzzzz_xxxxy_0[j] + fr * tg_yzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yzzzz_xxxy_1[j];

                    tg_xyzzzz_xxxxz_0[j] = pb_x * tg_yzzzz_xxxxz_0[j] + fr * tg_yzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yzzzz_xxxz_1[j];

                    tg_xyzzzz_xxxyy_0[j] = pb_x * tg_yzzzz_xxxyy_0[j] + fr * tg_yzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxyy_1[j];

                    tg_xyzzzz_xxxyz_0[j] = pb_x * tg_yzzzz_xxxyz_0[j] + fr * tg_yzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxyz_1[j];

                    tg_xyzzzz_xxxzz_0[j] = pb_x * tg_yzzzz_xxxzz_0[j] + fr * tg_yzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yzzzz_xxzz_1[j];

                    tg_xyzzzz_xxyyy_0[j] = pb_x * tg_yzzzz_xxyyy_0[j] + fr * tg_yzzzz_xxyyy_1[j] + fl1_fxn * tg_yzzzz_xyyy_1[j];

                    tg_xyzzzz_xxyyz_0[j] = pb_x * tg_yzzzz_xxyyz_0[j] + fr * tg_yzzzz_xxyyz_1[j] + fl1_fxn * tg_yzzzz_xyyz_1[j];

                    tg_xyzzzz_xxyzz_0[j] = pb_x * tg_yzzzz_xxyzz_0[j] + fr * tg_yzzzz_xxyzz_1[j] + fl1_fxn * tg_yzzzz_xyzz_1[j];

                    tg_xyzzzz_xxzzz_0[j] = pb_x * tg_yzzzz_xxzzz_0[j] + fr * tg_yzzzz_xxzzz_1[j] + fl1_fxn * tg_yzzzz_xzzz_1[j];

                    tg_xyzzzz_xyyyy_0[j] = pb_x * tg_yzzzz_xyyyy_0[j] + fr * tg_yzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyyy_1[j];

                    tg_xyzzzz_xyyyz_0[j] = pb_x * tg_yzzzz_xyyyz_0[j] + fr * tg_yzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyyz_1[j];

                    tg_xyzzzz_xyyzz_0[j] = pb_x * tg_yzzzz_xyyzz_0[j] + fr * tg_yzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yyzz_1[j];

                    tg_xyzzzz_xyzzz_0[j] = pb_x * tg_yzzzz_xyzzz_0[j] + fr * tg_yzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_yzzz_1[j];

                    tg_xyzzzz_xzzzz_0[j] = pb_x * tg_yzzzz_xzzzz_0[j] + fr * tg_yzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzz_zzzz_1[j];

                    tg_xyzzzz_yyyyy_0[j] = pb_x * tg_yzzzz_yyyyy_0[j] + fr * tg_yzzzz_yyyyy_1[j];

                    tg_xyzzzz_yyyyz_0[j] = pb_x * tg_yzzzz_yyyyz_0[j] + fr * tg_yzzzz_yyyyz_1[j];

                    tg_xyzzzz_yyyzz_0[j] = pb_x * tg_yzzzz_yyyzz_0[j] + fr * tg_yzzzz_yyyzz_1[j];

                    tg_xyzzzz_yyzzz_0[j] = pb_x * tg_yzzzz_yyzzz_0[j] + fr * tg_yzzzz_yyzzz_1[j];

                    tg_xyzzzz_yzzzz_0[j] = pb_x * tg_yzzzz_yzzzz_0[j] + fr * tg_yzzzz_yzzzz_1[j];

                    tg_xyzzzz_zzzzz_0[j] = pb_x * tg_yzzzz_zzzzz_0[j] + fr * tg_yzzzz_zzzzz_1[j];

                    tg_xzzzzz_xxxxx_0[j] = pb_x * tg_zzzzz_xxxxx_0[j] + fr * tg_zzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_zzzzz_xxxx_1[j];

                    tg_xzzzzz_xxxxy_0[j] = pb_x * tg_zzzzz_xxxxy_0[j] + fr * tg_zzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xxxy_1[j];

                    tg_xzzzzz_xxxxz_0[j] = pb_x * tg_zzzzz_xxxxz_0[j] + fr * tg_zzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xxxz_1[j];

                    tg_xzzzzz_xxxyy_0[j] = pb_x * tg_zzzzz_xxxyy_0[j] + fr * tg_zzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxyy_1[j];

                    tg_xzzzzz_xxxyz_0[j] = pb_x * tg_zzzzz_xxxyz_0[j] + fr * tg_zzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxyz_1[j];

                    tg_xzzzzz_xxxzz_0[j] = pb_x * tg_zzzzz_xxxzz_0[j] + fr * tg_zzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxzz_1[j];

                    tg_xzzzzz_xxyyy_0[j] = pb_x * tg_zzzzz_xxyyy_0[j] + fr * tg_zzzzz_xxyyy_1[j] + fl1_fxn * tg_zzzzz_xyyy_1[j];

                    tg_xzzzzz_xxyyz_0[j] = pb_x * tg_zzzzz_xxyyz_0[j] + fr * tg_zzzzz_xxyyz_1[j] + fl1_fxn * tg_zzzzz_xyyz_1[j];

                    tg_xzzzzz_xxyzz_0[j] = pb_x * tg_zzzzz_xxyzz_0[j] + fr * tg_zzzzz_xxyzz_1[j] + fl1_fxn * tg_zzzzz_xyzz_1[j];

                    tg_xzzzzz_xxzzz_0[j] = pb_x * tg_zzzzz_xxzzz_0[j] + fr * tg_zzzzz_xxzzz_1[j] + fl1_fxn * tg_zzzzz_xzzz_1[j];

                    tg_xzzzzz_xyyyy_0[j] = pb_x * tg_zzzzz_xyyyy_0[j] + fr * tg_zzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyyy_1[j];

                    tg_xzzzzz_xyyyz_0[j] = pb_x * tg_zzzzz_xyyyz_0[j] + fr * tg_zzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyyz_1[j];

                    tg_xzzzzz_xyyzz_0[j] = pb_x * tg_zzzzz_xyyzz_0[j] + fr * tg_zzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yyzz_1[j];

                    tg_xzzzzz_xyzzz_0[j] = pb_x * tg_zzzzz_xyzzz_0[j] + fr * tg_zzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_yzzz_1[j];

                    tg_xzzzzz_xzzzz_0[j] = pb_x * tg_zzzzz_xzzzz_0[j] + fr * tg_zzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_zzzz_1[j];

                    tg_xzzzzz_yyyyy_0[j] = pb_x * tg_zzzzz_yyyyy_0[j] + fr * tg_zzzzz_yyyyy_1[j];

                    tg_xzzzzz_yyyyz_0[j] = pb_x * tg_zzzzz_yyyyz_0[j] + fr * tg_zzzzz_yyyyz_1[j];

                    tg_xzzzzz_yyyzz_0[j] = pb_x * tg_zzzzz_yyyzz_0[j] + fr * tg_zzzzz_yyyzz_1[j];

                    tg_xzzzzz_yyzzz_0[j] = pb_x * tg_zzzzz_yyzzz_0[j] + fr * tg_zzzzz_yyzzz_1[j];

                    tg_xzzzzz_yzzzz_0[j] = pb_x * tg_zzzzz_yzzzz_0[j] + fr * tg_zzzzz_yzzzz_1[j];

                    tg_xzzzzz_zzzzz_0[j] = pb_x * tg_zzzzz_zzzzz_0[j] + fr * tg_zzzzz_zzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyyy_xxxxx_0[j] = pb_y * tg_yyyyy_xxxxx_0[j] + fr * tg_yyyyy_xxxxx_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxx_0[j] - tg_yyyy_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyy_xxxxy_0[j] = pb_y * tg_yyyyy_xxxxy_0[j] + fr * tg_yyyyy_xxxxy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxy_0[j] - tg_yyyy_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xxxx_1[j];

                    tg_yyyyyy_xxxxz_0[j] = pb_y * tg_yyyyy_xxxxz_0[j] + fr * tg_yyyyy_xxxxz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxxz_0[j] - tg_yyyy_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyy_xxxyy_0[j] = pb_y * tg_yyyyy_xxxyy_0[j] + fr * tg_yyyyy_xxxyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxyy_0[j] - tg_yyyy_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_xxxy_1[j];

                    tg_yyyyyy_xxxyz_0[j] = pb_y * tg_yyyyy_xxxyz_0[j] + fr * tg_yyyyy_xxxyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxyz_0[j] - tg_yyyy_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xxxz_1[j];

                    tg_yyyyyy_xxxzz_0[j] = pb_y * tg_yyyyy_xxxzz_0[j] + fr * tg_yyyyy_xxxzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxxzz_0[j] - tg_yyyy_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyy_xxyyy_0[j] = pb_y * tg_yyyyy_xxyyy_0[j] + fr * tg_yyyyy_xxyyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxyyy_0[j] - tg_yyyy_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyy_xxyy_1[j];

                    tg_yyyyyy_xxyyz_0[j] = pb_y * tg_yyyyy_xxyyz_0[j] + fr * tg_yyyyy_xxyyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxyyz_0[j] - tg_yyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_xxyz_1[j];

                    tg_yyyyyy_xxyzz_0[j] = pb_y * tg_yyyyy_xxyzz_0[j] + fr * tg_yyyyy_xxyzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxyzz_0[j] - tg_yyyy_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xxzz_1[j];

                    tg_yyyyyy_xxzzz_0[j] = pb_y * tg_yyyyy_xxzzz_0[j] + fr * tg_yyyyy_xxzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xxzzz_0[j] - tg_yyyy_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyy_xyyyy_0[j] = pb_y * tg_yyyyy_xyyyy_0[j] + fr * tg_yyyyy_xyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyyyy_0[j] - tg_yyyy_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyy_xyyy_1[j];

                    tg_yyyyyy_xyyyz_0[j] = pb_y * tg_yyyyy_xyyyz_0[j] + fr * tg_yyyyy_xyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyyyz_0[j] - tg_yyyy_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyy_xyyz_1[j];

                    tg_yyyyyy_xyyzz_0[j] = pb_y * tg_yyyyy_xyyzz_0[j] + fr * tg_yyyyy_xyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyyzz_0[j] - tg_yyyy_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_xyzz_1[j];

                    tg_yyyyyy_xyzzz_0[j] = pb_y * tg_yyyyy_xyzzz_0[j] + fr * tg_yyyyy_xyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xyzzz_0[j] - tg_yyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_xzzz_1[j];

                    tg_yyyyyy_xzzzz_0[j] = pb_y * tg_yyyyy_xzzzz_0[j] + fr * tg_yyyyy_xzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_xzzzz_0[j] - tg_yyyy_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyy_yyyyy_0[j] = pb_y * tg_yyyyy_yyyyy_0[j] + fr * tg_yyyyy_yyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyyyy_0[j] - tg_yyyy_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyy_yyyy_1[j];

                    tg_yyyyyy_yyyyz_0[j] = pb_y * tg_yyyyy_yyyyz_0[j] + fr * tg_yyyyy_yyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyyyz_0[j] - tg_yyyy_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyy_yyyz_1[j];

                    tg_yyyyyy_yyyzz_0[j] = pb_y * tg_yyyyy_yyyzz_0[j] + fr * tg_yyyyy_yyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyyzz_0[j] - tg_yyyy_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyy_yyzz_1[j];

                    tg_yyyyyy_yyzzz_0[j] = pb_y * tg_yyyyy_yyzzz_0[j] + fr * tg_yyyyy_yyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yyzzz_0[j] - tg_yyyy_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyy_yzzz_1[j];

                    tg_yyyyyy_yzzzz_0[j] = pb_y * tg_yyyyy_yzzzz_0[j] + fr * tg_yyyyy_yzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_yzzzz_0[j] - tg_yyyy_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyy_zzzz_1[j];

                    tg_yyyyyy_zzzzz_0[j] = pb_y * tg_yyyyy_zzzzz_0[j] + fr * tg_yyyyy_zzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyy_zzzzz_0[j] - tg_yyyy_zzzzz_1[j] * fl1_fza);

                    tg_yyyyyz_xxxxx_0[j] = pb_y * tg_yyyyz_xxxxx_0[j] + fr * tg_yyyyz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxx_0[j] - tg_yyyz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyz_xxxxy_0[j] = pb_y * tg_yyyyz_xxxxy_0[j] + fr * tg_yyyyz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxy_0[j] - tg_yyyz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xxxx_1[j];

                    tg_yyyyyz_xxxxz_0[j] = pb_y * tg_yyyyz_xxxxz_0[j] + fr * tg_yyyyz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxxz_0[j] - tg_yyyz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyz_xxxyy_0[j] = pb_y * tg_yyyyz_xxxyy_0[j] + fr * tg_yyyyz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxyy_0[j] - tg_yyyz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_xxxy_1[j];

                    tg_yyyyyz_xxxyz_0[j] = pb_y * tg_yyyyz_xxxyz_0[j] + fr * tg_yyyyz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxyz_0[j] - tg_yyyz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xxxz_1[j];

                    tg_yyyyyz_xxxzz_0[j] = pb_y * tg_yyyyz_xxxzz_0[j] + fr * tg_yyyyz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxxzz_0[j] - tg_yyyz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyz_xxyyy_0[j] = pb_y * tg_yyyyz_xxyyy_0[j] + fr * tg_yyyyz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxyyy_0[j] - tg_yyyz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyz_xxyy_1[j];

                    tg_yyyyyz_xxyyz_0[j] = pb_y * tg_yyyyz_xxyyz_0[j] + fr * tg_yyyyz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxyyz_0[j] - tg_yyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_xxyz_1[j];

                    tg_yyyyyz_xxyzz_0[j] = pb_y * tg_yyyyz_xxyzz_0[j] + fr * tg_yyyyz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxyzz_0[j] - tg_yyyz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xxzz_1[j];

                    tg_yyyyyz_xxzzz_0[j] = pb_y * tg_yyyyz_xxzzz_0[j] + fr * tg_yyyyz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xxzzz_0[j] - tg_yyyz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyz_xyyyy_0[j] = pb_y * tg_yyyyz_xyyyy_0[j] + fr * tg_yyyyz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyyyy_0[j] - tg_yyyz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyz_xyyy_1[j];

                    tg_yyyyyz_xyyyz_0[j] = pb_y * tg_yyyyz_xyyyz_0[j] + fr * tg_yyyyz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyyyz_0[j] - tg_yyyz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyz_xyyz_1[j];

                    tg_yyyyyz_xyyzz_0[j] = pb_y * tg_yyyyz_xyyzz_0[j] + fr * tg_yyyyz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyyzz_0[j] - tg_yyyz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_xyzz_1[j];

                    tg_yyyyyz_xyzzz_0[j] = pb_y * tg_yyyyz_xyzzz_0[j] + fr * tg_yyyyz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xyzzz_0[j] - tg_yyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_xzzz_1[j];

                    tg_yyyyyz_xzzzz_0[j] = pb_y * tg_yyyyz_xzzzz_0[j] + fr * tg_yyyyz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_xzzzz_0[j] - tg_yyyz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyz_yyyyy_0[j] = pb_y * tg_yyyyz_yyyyy_0[j] + fr * tg_yyyyz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyyyy_0[j] - tg_yyyz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyz_yyyy_1[j];

                    tg_yyyyyz_yyyyz_0[j] = pb_y * tg_yyyyz_yyyyz_0[j] + fr * tg_yyyyz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyyyz_0[j] - tg_yyyz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyz_yyyz_1[j];

                    tg_yyyyyz_yyyzz_0[j] = pb_y * tg_yyyyz_yyyzz_0[j] + fr * tg_yyyyz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyyzz_0[j] - tg_yyyz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyz_yyzz_1[j];

                    tg_yyyyyz_yyzzz_0[j] = pb_y * tg_yyyyz_yyzzz_0[j] + fr * tg_yyyyz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yyzzz_0[j] - tg_yyyz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyz_yzzz_1[j];

                    tg_yyyyyz_yzzzz_0[j] = pb_y * tg_yyyyz_yzzzz_0[j] + fr * tg_yyyyz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_yzzzz_0[j] - tg_yyyz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyz_zzzz_1[j];

                    tg_yyyyyz_zzzzz_0[j] = pb_y * tg_yyyyz_zzzzz_0[j] + fr * tg_yyyyz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyz_zzzzz_0[j] - tg_yyyz_zzzzz_1[j] * fl1_fza);

                    tg_yyyyzz_xxxxx_0[j] = pb_y * tg_yyyzz_xxxxx_0[j] + fr * tg_yyyzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxx_0[j] - tg_yyzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyzz_xxxxy_0[j] = pb_y * tg_yyyzz_xxxxy_0[j] + fr * tg_yyyzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxy_0[j] - tg_yyzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xxxx_1[j];

                    tg_yyyyzz_xxxxz_0[j] = pb_y * tg_yyyzz_xxxxz_0[j] + fr * tg_yyyzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxxz_0[j] - tg_yyzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyzz_xxxyy_0[j] = pb_y * tg_yyyzz_xxxyy_0[j] + fr * tg_yyyzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxyy_0[j] - tg_yyzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_xxxy_1[j];

                    tg_yyyyzz_xxxyz_0[j] = pb_y * tg_yyyzz_xxxyz_0[j] + fr * tg_yyyzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxyz_0[j] - tg_yyzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xxxz_1[j];

                    tg_yyyyzz_xxxzz_0[j] = pb_y * tg_yyyzz_xxxzz_0[j] + fr * tg_yyyzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxxzz_0[j] - tg_yyzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyzz_xxyyy_0[j] = pb_y * tg_yyyzz_xxyyy_0[j] + fr * tg_yyyzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxyyy_0[j] - tg_yyzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzz_xxyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISH_490_588(      CMemBlock2D<double>* primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {6, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 259); 

                auto tg_yyzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 260); 

                auto tg_yyzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 261); 

                auto tg_yyzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 262); 

                auto tg_yyzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 263); 

                auto tg_yyzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 264); 

                auto tg_yyzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 265); 

                auto tg_yyzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 266); 

                auto tg_yyzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 267); 

                auto tg_yyzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 268); 

                auto tg_yyzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 269); 

                auto tg_yyzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 270); 

                auto tg_yyzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 271); 

                auto tg_yyzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 272); 

                auto tg_yzzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 273); 

                auto tg_yzzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 274); 

                auto tg_yzzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 275); 

                auto tg_yzzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 276); 

                auto tg_yzzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 277); 

                auto tg_yzzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 278); 

                auto tg_yzzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 279); 

                auto tg_yzzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 280); 

                auto tg_yzzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 281); 

                auto tg_yzzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 282); 

                auto tg_yzzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 283); 

                auto tg_yzzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 284); 

                auto tg_yzzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 285); 

                auto tg_yzzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 286); 

                auto tg_yzzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 287); 

                auto tg_yzzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 288); 

                auto tg_yzzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 289); 

                auto tg_yzzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 290); 

                auto tg_yzzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 291); 

                auto tg_yzzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 292); 

                auto tg_yzzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 293); 

                auto tg_zzzz_xxxxx_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 294); 

                auto tg_zzzz_xxxxy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 295); 

                auto tg_zzzz_xxxxz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 296); 

                auto tg_zzzz_xxxyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 297); 

                auto tg_zzzz_xxxyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 298); 

                auto tg_zzzz_xxxzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 299); 

                auto tg_zzzz_xxyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 300); 

                auto tg_zzzz_xxyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 301); 

                auto tg_zzzz_xxyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 302); 

                auto tg_zzzz_xxzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 303); 

                auto tg_zzzz_xyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 304); 

                auto tg_zzzz_xyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 305); 

                auto tg_zzzz_xyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 306); 

                auto tg_zzzz_xyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 307); 

                auto tg_zzzz_xzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 308); 

                auto tg_zzzz_yyyyy_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 309); 

                auto tg_zzzz_yyyyz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 310); 

                auto tg_zzzz_yyyzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 311); 

                auto tg_zzzz_yyzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 312); 

                auto tg_zzzz_yzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 313); 

                auto tg_zzzz_zzzzz_0 = primBuffer[pidx_g_4_5_m0].data(315 * idx + 314); 

                auto tg_yyzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 259); 

                auto tg_yyzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 260); 

                auto tg_yyzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 261); 

                auto tg_yyzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 262); 

                auto tg_yyzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 263); 

                auto tg_yyzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 264); 

                auto tg_yyzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 265); 

                auto tg_yyzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 266); 

                auto tg_yyzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 267); 

                auto tg_yyzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 268); 

                auto tg_yyzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 269); 

                auto tg_yyzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 270); 

                auto tg_yyzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 271); 

                auto tg_yyzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 272); 

                auto tg_yzzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 273); 

                auto tg_yzzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 274); 

                auto tg_yzzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 275); 

                auto tg_yzzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 276); 

                auto tg_yzzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 277); 

                auto tg_yzzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 278); 

                auto tg_yzzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 279); 

                auto tg_yzzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 280); 

                auto tg_yzzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 281); 

                auto tg_yzzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 282); 

                auto tg_yzzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 283); 

                auto tg_yzzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 284); 

                auto tg_yzzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 285); 

                auto tg_yzzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 286); 

                auto tg_yzzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 287); 

                auto tg_yzzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 288); 

                auto tg_yzzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 289); 

                auto tg_yzzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 290); 

                auto tg_yzzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 291); 

                auto tg_yzzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 292); 

                auto tg_yzzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 293); 

                auto tg_zzzz_xxxxx_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 294); 

                auto tg_zzzz_xxxxy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 295); 

                auto tg_zzzz_xxxxz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 296); 

                auto tg_zzzz_xxxyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 297); 

                auto tg_zzzz_xxxyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 298); 

                auto tg_zzzz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 299); 

                auto tg_zzzz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 300); 

                auto tg_zzzz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 301); 

                auto tg_zzzz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 302); 

                auto tg_zzzz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 303); 

                auto tg_zzzz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 304); 

                auto tg_zzzz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 305); 

                auto tg_zzzz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 306); 

                auto tg_zzzz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 307); 

                auto tg_zzzz_xzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 308); 

                auto tg_zzzz_yyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 309); 

                auto tg_zzzz_yyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 310); 

                auto tg_zzzz_yyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 311); 

                auto tg_zzzz_yyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 312); 

                auto tg_zzzz_yzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 313); 

                auto tg_zzzz_zzzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 314); 

                auto tg_yyyzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 259); 

                auto tg_yyyzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 260); 

                auto tg_yyyzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 261); 

                auto tg_yyyzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 262); 

                auto tg_yyyzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 263); 

                auto tg_yyyzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 264); 

                auto tg_yyyzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 265); 

                auto tg_yyyzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 266); 

                auto tg_yyyzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 267); 

                auto tg_yyyzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 268); 

                auto tg_yyyzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 269); 

                auto tg_yyzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 270); 

                auto tg_yyzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 271); 

                auto tg_yyzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 272); 

                auto tg_yyzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 273); 

                auto tg_yyzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 274); 

                auto tg_yyzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 275); 

                auto tg_yyzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 276); 

                auto tg_yyzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 277); 

                auto tg_yyzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 278); 

                auto tg_yyzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 279); 

                auto tg_yyzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 280); 

                auto tg_yyzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 281); 

                auto tg_yyzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 282); 

                auto tg_yyzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 283); 

                auto tg_yyzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 284); 

                auto tg_yzzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 285); 

                auto tg_yzzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 286); 

                auto tg_yzzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 287); 

                auto tg_yzzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 288); 

                auto tg_yzzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 289); 

                auto tg_yzzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 290); 

                auto tg_yzzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 291); 

                auto tg_yzzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 292); 

                auto tg_yzzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 293); 

                auto tg_yzzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 294); 

                auto tg_yzzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 295); 

                auto tg_yzzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 296); 

                auto tg_yzzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 297); 

                auto tg_yzzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 298); 

                auto tg_yzzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 299); 

                auto tg_zzzzz_xxxx_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 300); 

                auto tg_zzzzz_xxxy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 301); 

                auto tg_zzzzz_xxxz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 302); 

                auto tg_zzzzz_xxyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 303); 

                auto tg_zzzzz_xxyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 304); 

                auto tg_zzzzz_xxzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 305); 

                auto tg_zzzzz_xyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 306); 

                auto tg_zzzzz_xyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 307); 

                auto tg_zzzzz_xyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 308); 

                auto tg_zzzzz_xzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 309); 

                auto tg_zzzzz_yyyy_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 310); 

                auto tg_zzzzz_yyyz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 311); 

                auto tg_zzzzz_yyzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 312); 

                auto tg_zzzzz_yzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 313); 

                auto tg_zzzzz_zzzz_1 = primBuffer[pidx_g_5_4_m1].data(315 * idx + 314); 

                // set up pointers to integrals

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

                // Batch of Integrals (490,588)

                #pragma omp simd aligned(fxn, fza, tg_yyyyzz_xxyyz_0, tg_yyyyzz_xxyzz_0, tg_yyyyzz_xxzzz_0, \
                                         tg_yyyyzz_xyyyy_0, tg_yyyyzz_xyyyz_0, tg_yyyyzz_xyyzz_0, tg_yyyyzz_xyzzz_0, \
                                         tg_yyyyzz_xzzzz_0, tg_yyyyzz_yyyyy_0, tg_yyyyzz_yyyyz_0, tg_yyyyzz_yyyzz_0, \
                                         tg_yyyyzz_yyzzz_0, tg_yyyyzz_yzzzz_0, tg_yyyyzz_zzzzz_0, tg_yyyzz_xxyyz_0, \
                                         tg_yyyzz_xxyyz_1, tg_yyyzz_xxyz_1, tg_yyyzz_xxyzz_0, tg_yyyzz_xxyzz_1, \
                                         tg_yyyzz_xxzz_1, tg_yyyzz_xxzzz_0, tg_yyyzz_xxzzz_1, tg_yyyzz_xyyy_1, \
                                         tg_yyyzz_xyyyy_0, tg_yyyzz_xyyyy_1, tg_yyyzz_xyyyz_0, tg_yyyzz_xyyyz_1, \
                                         tg_yyyzz_xyyz_1, tg_yyyzz_xyyzz_0, tg_yyyzz_xyyzz_1, tg_yyyzz_xyzz_1, \
                                         tg_yyyzz_xyzzz_0, tg_yyyzz_xyzzz_1, tg_yyyzz_xzzz_1, tg_yyyzz_xzzzz_0, \
                                         tg_yyyzz_xzzzz_1, tg_yyyzz_yyyy_1, tg_yyyzz_yyyyy_0, tg_yyyzz_yyyyy_1, \
                                         tg_yyyzz_yyyyz_0, tg_yyyzz_yyyyz_1, tg_yyyzz_yyyz_1, tg_yyyzz_yyyzz_0, \
                                         tg_yyyzz_yyyzz_1, tg_yyyzz_yyzz_1, tg_yyyzz_yyzzz_0, tg_yyyzz_yyzzz_1, \
                                         tg_yyyzz_yzzz_1, tg_yyyzz_yzzzz_0, tg_yyyzz_yzzzz_1, tg_yyyzz_zzzz_1, \
                                         tg_yyyzz_zzzzz_0, tg_yyyzz_zzzzz_1, tg_yyyzzz_xxxxx_0, tg_yyyzzz_xxxxy_0, \
                                         tg_yyyzzz_xxxxz_0, tg_yyyzzz_xxxyy_0, tg_yyyzzz_xxxyz_0, tg_yyyzzz_xxxzz_0, \
                                         tg_yyyzzz_xxyyy_0, tg_yyyzzz_xxyyz_0, tg_yyyzzz_xxyzz_0, tg_yyyzzz_xxzzz_0, \
                                         tg_yyyzzz_xyyyy_0, tg_yyyzzz_xyyyz_0, tg_yyyzzz_xyyzz_0, tg_yyyzzz_xyzzz_0, \
                                         tg_yyyzzz_xzzzz_0, tg_yyyzzz_yyyyy_0, tg_yyyzzz_yyyyz_0, tg_yyyzzz_yyyzz_0, \
                                         tg_yyyzzz_yyzzz_0, tg_yyyzzz_yzzzz_0, tg_yyyzzz_zzzzz_0, tg_yyzz_xxyyz_0, \
                                         tg_yyzz_xxyyz_1, tg_yyzz_xxyzz_0, tg_yyzz_xxyzz_1, tg_yyzz_xxzzz_0, tg_yyzz_xxzzz_1, \
                                         tg_yyzz_xyyyy_0, tg_yyzz_xyyyy_1, tg_yyzz_xyyyz_0, tg_yyzz_xyyyz_1, tg_yyzz_xyyzz_0, \
                                         tg_yyzz_xyyzz_1, tg_yyzz_xyzzz_0, tg_yyzz_xyzzz_1, tg_yyzz_xzzzz_0, tg_yyzz_xzzzz_1, \
                                         tg_yyzz_yyyyy_0, tg_yyzz_yyyyy_1, tg_yyzz_yyyyz_0, tg_yyzz_yyyyz_1, tg_yyzz_yyyzz_0, \
                                         tg_yyzz_yyyzz_1, tg_yyzz_yyzzz_0, tg_yyzz_yyzzz_1, tg_yyzz_yzzzz_0, tg_yyzz_yzzzz_1, \
                                         tg_yyzz_zzzzz_0, tg_yyzz_zzzzz_1, tg_yyzzz_xxxx_1, tg_yyzzz_xxxxx_0, \
                                         tg_yyzzz_xxxxx_1, tg_yyzzz_xxxxy_0, tg_yyzzz_xxxxy_1, tg_yyzzz_xxxxz_0, \
                                         tg_yyzzz_xxxxz_1, tg_yyzzz_xxxy_1, tg_yyzzz_xxxyy_0, tg_yyzzz_xxxyy_1, \
                                         tg_yyzzz_xxxyz_0, tg_yyzzz_xxxyz_1, tg_yyzzz_xxxz_1, tg_yyzzz_xxxzz_0, \
                                         tg_yyzzz_xxxzz_1, tg_yyzzz_xxyy_1, tg_yyzzz_xxyyy_0, tg_yyzzz_xxyyy_1, \
                                         tg_yyzzz_xxyyz_0, tg_yyzzz_xxyyz_1, tg_yyzzz_xxyz_1, tg_yyzzz_xxyzz_0, \
                                         tg_yyzzz_xxyzz_1, tg_yyzzz_xxzz_1, tg_yyzzz_xxzzz_0, tg_yyzzz_xxzzz_1, \
                                         tg_yyzzz_xyyy_1, tg_yyzzz_xyyyy_0, tg_yyzzz_xyyyy_1, tg_yyzzz_xyyyz_0, \
                                         tg_yyzzz_xyyyz_1, tg_yyzzz_xyyz_1, tg_yyzzz_xyyzz_0, tg_yyzzz_xyyzz_1, \
                                         tg_yyzzz_xyzz_1, tg_yyzzz_xyzzz_0, tg_yyzzz_xyzzz_1, tg_yyzzz_xzzz_1, \
                                         tg_yyzzz_xzzzz_0, tg_yyzzz_xzzzz_1, tg_yyzzz_yyyy_1, tg_yyzzz_yyyyy_0, \
                                         tg_yyzzz_yyyyy_1, tg_yyzzz_yyyyz_0, tg_yyzzz_yyyyz_1, tg_yyzzz_yyyz_1, \
                                         tg_yyzzz_yyyzz_0, tg_yyzzz_yyyzz_1, tg_yyzzz_yyzz_1, tg_yyzzz_yyzzz_0, \
                                         tg_yyzzz_yyzzz_1, tg_yyzzz_yzzz_1, tg_yyzzz_yzzzz_0, tg_yyzzz_yzzzz_1, \
                                         tg_yyzzz_zzzz_1, tg_yyzzz_zzzzz_0, tg_yyzzz_zzzzz_1, tg_yyzzzz_xxxxx_0, \
                                         tg_yyzzzz_xxxxy_0, tg_yyzzzz_xxxxz_0, tg_yyzzzz_xxxyy_0, tg_yyzzzz_xxxyz_0, \
                                         tg_yyzzzz_xxxzz_0, tg_yyzzzz_xxyyy_0, tg_yyzzzz_xxyyz_0, tg_yyzzzz_xxyzz_0, \
                                         tg_yyzzzz_xxzzz_0, tg_yyzzzz_xyyyy_0, tg_yyzzzz_xyyyz_0, tg_yyzzzz_xyyzz_0, \
                                         tg_yyzzzz_xyzzz_0, tg_yyzzzz_xzzzz_0, tg_yyzzzz_yyyyy_0, tg_yyzzzz_yyyyz_0, \
                                         tg_yyzzzz_yyyzz_0, tg_yyzzzz_yyzzz_0, tg_yyzzzz_yzzzz_0, tg_yyzzzz_zzzzz_0, \
                                         tg_yzzz_xxxxx_0, tg_yzzz_xxxxx_1, tg_yzzz_xxxxy_0, tg_yzzz_xxxxy_1, tg_yzzz_xxxxz_0, \
                                         tg_yzzz_xxxxz_1, tg_yzzz_xxxyy_0, tg_yzzz_xxxyy_1, tg_yzzz_xxxyz_0, tg_yzzz_xxxyz_1, \
                                         tg_yzzz_xxxzz_0, tg_yzzz_xxxzz_1, tg_yzzz_xxyyy_0, tg_yzzz_xxyyy_1, tg_yzzz_xxyyz_0, \
                                         tg_yzzz_xxyyz_1, tg_yzzz_xxyzz_0, tg_yzzz_xxyzz_1, tg_yzzz_xxzzz_0, tg_yzzz_xxzzz_1, \
                                         tg_yzzz_xyyyy_0, tg_yzzz_xyyyy_1, tg_yzzz_xyyyz_0, tg_yzzz_xyyyz_1, tg_yzzz_xyyzz_0, \
                                         tg_yzzz_xyyzz_1, tg_yzzz_xyzzz_0, tg_yzzz_xyzzz_1, tg_yzzz_xzzzz_0, tg_yzzz_xzzzz_1, \
                                         tg_yzzz_yyyyy_0, tg_yzzz_yyyyy_1, tg_yzzz_yyyyz_0, tg_yzzz_yyyyz_1, tg_yzzz_yyyzz_0, \
                                         tg_yzzz_yyyzz_1, tg_yzzz_yyzzz_0, tg_yzzz_yyzzz_1, tg_yzzz_yzzzz_0, tg_yzzz_yzzzz_1, \
                                         tg_yzzz_zzzzz_0, tg_yzzz_zzzzz_1, tg_yzzzz_xxxx_1, tg_yzzzz_xxxxx_0, \
                                         tg_yzzzz_xxxxx_1, tg_yzzzz_xxxxy_0, tg_yzzzz_xxxxy_1, tg_yzzzz_xxxxz_0, \
                                         tg_yzzzz_xxxxz_1, tg_yzzzz_xxxy_1, tg_yzzzz_xxxyy_0, tg_yzzzz_xxxyy_1, \
                                         tg_yzzzz_xxxyz_0, tg_yzzzz_xxxyz_1, tg_yzzzz_xxxz_1, tg_yzzzz_xxxzz_0, \
                                         tg_yzzzz_xxxzz_1, tg_yzzzz_xxyy_1, tg_yzzzz_xxyyy_0, tg_yzzzz_xxyyy_1, \
                                         tg_yzzzz_xxyyz_0, tg_yzzzz_xxyyz_1, tg_yzzzz_xxyz_1, tg_yzzzz_xxyzz_0, \
                                         tg_yzzzz_xxyzz_1, tg_yzzzz_xxzz_1, tg_yzzzz_xxzzz_0, tg_yzzzz_xxzzz_1, \
                                         tg_yzzzz_xyyy_1, tg_yzzzz_xyyyy_0, tg_yzzzz_xyyyy_1, tg_yzzzz_xyyyz_0, \
                                         tg_yzzzz_xyyyz_1, tg_yzzzz_xyyz_1, tg_yzzzz_xyyzz_0, tg_yzzzz_xyyzz_1, \
                                         tg_yzzzz_xyzz_1, tg_yzzzz_xyzzz_0, tg_yzzzz_xyzzz_1, tg_yzzzz_xzzz_1, \
                                         tg_yzzzz_xzzzz_0, tg_yzzzz_xzzzz_1, tg_yzzzz_yyyy_1, tg_yzzzz_yyyyy_0, \
                                         tg_yzzzz_yyyyy_1, tg_yzzzz_yyyyz_0, tg_yzzzz_yyyyz_1, tg_yzzzz_yyyz_1, \
                                         tg_yzzzz_yyyzz_0, tg_yzzzz_yyyzz_1, tg_yzzzz_yyzz_1, tg_yzzzz_yyzzz_0, \
                                         tg_yzzzz_yyzzz_1, tg_yzzzz_yzzz_1, tg_yzzzz_yzzzz_0, tg_yzzzz_yzzzz_1, \
                                         tg_yzzzz_zzzz_1, tg_yzzzz_zzzzz_0, tg_yzzzz_zzzzz_1, tg_yzzzzz_xxxxx_0, \
                                         tg_yzzzzz_xxxxy_0, tg_yzzzzz_xxxxz_0, tg_yzzzzz_xxxyy_0, tg_yzzzzz_xxxyz_0, \
                                         tg_yzzzzz_xxxzz_0, tg_yzzzzz_xxyyy_0, tg_yzzzzz_xxyyz_0, tg_yzzzzz_xxyzz_0, \
                                         tg_yzzzzz_xxzzz_0, tg_yzzzzz_xyyyy_0, tg_yzzzzz_xyyyz_0, tg_yzzzzz_xyyzz_0, \
                                         tg_yzzzzz_xyzzz_0, tg_yzzzzz_xzzzz_0, tg_yzzzzz_yyyyy_0, tg_yzzzzz_yyyyz_0, \
                                         tg_yzzzzz_yyyzz_0, tg_yzzzzz_yyzzz_0, tg_yzzzzz_yzzzz_0, tg_yzzzzz_zzzzz_0, \
                                         tg_zzzz_xxxxx_0, tg_zzzz_xxxxx_1, tg_zzzz_xxxxy_0, tg_zzzz_xxxxy_1, tg_zzzz_xxxxz_0, \
                                         tg_zzzz_xxxxz_1, tg_zzzz_xxxyy_0, tg_zzzz_xxxyy_1, tg_zzzz_xxxyz_0, tg_zzzz_xxxyz_1, \
                                         tg_zzzz_xxxzz_0, tg_zzzz_xxxzz_1, tg_zzzz_xxyyy_0, tg_zzzz_xxyyy_1, tg_zzzz_xxyyz_0, \
                                         tg_zzzz_xxyyz_1, tg_zzzz_xxyzz_0, tg_zzzz_xxyzz_1, tg_zzzz_xxzzz_0, tg_zzzz_xxzzz_1, \
                                         tg_zzzz_xyyyy_0, tg_zzzz_xyyyy_1, tg_zzzz_xyyyz_0, tg_zzzz_xyyyz_1, tg_zzzz_xyyzz_0, \
                                         tg_zzzz_xyyzz_1, tg_zzzz_xyzzz_0, tg_zzzz_xyzzz_1, tg_zzzz_xzzzz_0, tg_zzzz_xzzzz_1, \
                                         tg_zzzz_yyyyy_0, tg_zzzz_yyyyy_1, tg_zzzz_yyyyz_0, tg_zzzz_yyyyz_1, tg_zzzz_yyyzz_0, \
                                         tg_zzzz_yyyzz_1, tg_zzzz_yyzzz_0, tg_zzzz_yyzzz_1, tg_zzzz_yzzzz_0, tg_zzzz_yzzzz_1, \
                                         tg_zzzz_zzzzz_0, tg_zzzz_zzzzz_1, tg_zzzzz_xxxx_1, tg_zzzzz_xxxxx_0, \
                                         tg_zzzzz_xxxxx_1, tg_zzzzz_xxxxy_0, tg_zzzzz_xxxxy_1, tg_zzzzz_xxxxz_0, \
                                         tg_zzzzz_xxxxz_1, tg_zzzzz_xxxy_1, tg_zzzzz_xxxyy_0, tg_zzzzz_xxxyy_1, \
                                         tg_zzzzz_xxxyz_0, tg_zzzzz_xxxyz_1, tg_zzzzz_xxxz_1, tg_zzzzz_xxxzz_0, \
                                         tg_zzzzz_xxxzz_1, tg_zzzzz_xxyy_1, tg_zzzzz_xxyyy_0, tg_zzzzz_xxyyy_1, \
                                         tg_zzzzz_xxyyz_0, tg_zzzzz_xxyyz_1, tg_zzzzz_xxyz_1, tg_zzzzz_xxyzz_0, \
                                         tg_zzzzz_xxyzz_1, tg_zzzzz_xxzz_1, tg_zzzzz_xxzzz_0, tg_zzzzz_xxzzz_1, \
                                         tg_zzzzz_xyyy_1, tg_zzzzz_xyyyy_0, tg_zzzzz_xyyyy_1, tg_zzzzz_xyyyz_0, \
                                         tg_zzzzz_xyyyz_1, tg_zzzzz_xyyz_1, tg_zzzzz_xyyzz_0, tg_zzzzz_xyyzz_1, \
                                         tg_zzzzz_xyzz_1, tg_zzzzz_xyzzz_0, tg_zzzzz_xyzzz_1, tg_zzzzz_xzzz_1, \
                                         tg_zzzzz_xzzzz_0, tg_zzzzz_xzzzz_1, tg_zzzzz_yyyy_1, tg_zzzzz_yyyyy_0, \
                                         tg_zzzzz_yyyyy_1, tg_zzzzz_yyyyz_0, tg_zzzzz_yyyyz_1, tg_zzzzz_yyyz_1, \
                                         tg_zzzzz_yyyzz_0, tg_zzzzz_yyyzz_1, tg_zzzzz_yyzz_1, tg_zzzzz_yyzzz_0, \
                                         tg_zzzzz_yyzzz_1, tg_zzzzz_yzzz_1, tg_zzzzz_yzzzz_0, tg_zzzzz_yzzzz_1, \
                                         tg_zzzzz_zzzz_1, tg_zzzzz_zzzzz_0, tg_zzzzz_zzzzz_1, tg_zzzzzz_xxxxx_0, \
                                         tg_zzzzzz_xxxxy_0, tg_zzzzzz_xxxxz_0, tg_zzzzzz_xxxyy_0, tg_zzzzzz_xxxyz_0, \
                                         tg_zzzzzz_xxxzz_0, tg_zzzzzz_xxyyy_0, tg_zzzzzz_xxyyz_0, tg_zzzzzz_xxyzz_0, \
                                         tg_zzzzzz_xxzzz_0, tg_zzzzzz_xyyyy_0, tg_zzzzzz_xyyyz_0, tg_zzzzzz_xyyzz_0, \
                                         tg_zzzzzz_xyzzz_0, tg_zzzzzz_xzzzz_0, tg_zzzzzz_yyyyy_0, tg_zzzzzz_yyyyz_0, \
                                         tg_zzzzzz_yyyzz_0, tg_zzzzzz_yyzzz_0, tg_zzzzzz_yzzzz_0, tg_zzzzzz_zzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyzz_xxyyz_0[j] = pb_y * tg_yyyzz_xxyyz_0[j] + fr * tg_yyyzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxyyz_0[j] - tg_yyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_xxyz_1[j];

                    tg_yyyyzz_xxyzz_0[j] = pb_y * tg_yyyzz_xxyzz_0[j] + fr * tg_yyyzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxyzz_0[j] - tg_yyzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xxzz_1[j];

                    tg_yyyyzz_xxzzz_0[j] = pb_y * tg_yyyzz_xxzzz_0[j] + fr * tg_yyyzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xxzzz_0[j] - tg_yyzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyzz_xyyyy_0[j] = pb_y * tg_yyyzz_xyyyy_0[j] + fr * tg_yyyzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyyyy_0[j] - tg_yyzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzz_xyyy_1[j];

                    tg_yyyyzz_xyyyz_0[j] = pb_y * tg_yyyzz_xyyyz_0[j] + fr * tg_yyyzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyyyz_0[j] - tg_yyzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzz_xyyz_1[j];

                    tg_yyyyzz_xyyzz_0[j] = pb_y * tg_yyyzz_xyyzz_0[j] + fr * tg_yyyzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyyzz_0[j] - tg_yyzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_xyzz_1[j];

                    tg_yyyyzz_xyzzz_0[j] = pb_y * tg_yyyzz_xyzzz_0[j] + fr * tg_yyyzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xyzzz_0[j] - tg_yyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_xzzz_1[j];

                    tg_yyyyzz_xzzzz_0[j] = pb_y * tg_yyyzz_xzzzz_0[j] + fr * tg_yyyzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_xzzzz_0[j] - tg_yyzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyzz_yyyyy_0[j] = pb_y * tg_yyyzz_yyyyy_0[j] + fr * tg_yyyzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyyyy_0[j] - tg_yyzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyzz_yyyy_1[j];

                    tg_yyyyzz_yyyyz_0[j] = pb_y * tg_yyyzz_yyyyz_0[j] + fr * tg_yyyzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyyyz_0[j] - tg_yyzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzz_yyyz_1[j];

                    tg_yyyyzz_yyyzz_0[j] = pb_y * tg_yyyzz_yyyzz_0[j] + fr * tg_yyyzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyyzz_0[j] - tg_yyzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzz_yyzz_1[j];

                    tg_yyyyzz_yyzzz_0[j] = pb_y * tg_yyyzz_yyzzz_0[j] + fr * tg_yyyzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yyzzz_0[j] - tg_yyzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzz_yzzz_1[j];

                    tg_yyyyzz_yzzzz_0[j] = pb_y * tg_yyyzz_yzzzz_0[j] + fr * tg_yyyzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_yzzzz_0[j] - tg_yyzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzz_zzzz_1[j];

                    tg_yyyyzz_zzzzz_0[j] = pb_y * tg_yyyzz_zzzzz_0[j] + fr * tg_yyyzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzz_zzzzz_0[j] - tg_yyzz_zzzzz_1[j] * fl1_fza);

                    tg_yyyzzz_xxxxx_0[j] = pb_y * tg_yyzzz_xxxxx_0[j] + fr * tg_yyzzz_xxxxx_1[j] + fl1_fx * (tg_yzzz_xxxxx_0[j] - tg_yzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyzzz_xxxxy_0[j] = pb_y * tg_yyzzz_xxxxy_0[j] + fr * tg_yyzzz_xxxxy_1[j] + fl1_fx * (tg_yzzz_xxxxy_0[j] - tg_yzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xxxx_1[j];

                    tg_yyyzzz_xxxxz_0[j] = pb_y * tg_yyzzz_xxxxz_0[j] + fr * tg_yyzzz_xxxxz_1[j] + fl1_fx * (tg_yzzz_xxxxz_0[j] - tg_yzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyzzz_xxxyy_0[j] = pb_y * tg_yyzzz_xxxyy_0[j] + fr * tg_yyzzz_xxxyy_1[j] + fl1_fx * (tg_yzzz_xxxyy_0[j] - tg_yzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_xxxy_1[j];

                    tg_yyyzzz_xxxyz_0[j] = pb_y * tg_yyzzz_xxxyz_0[j] + fr * tg_yyzzz_xxxyz_1[j] + fl1_fx * (tg_yzzz_xxxyz_0[j] - tg_yzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xxxz_1[j];

                    tg_yyyzzz_xxxzz_0[j] = pb_y * tg_yyzzz_xxxzz_0[j] + fr * tg_yyzzz_xxxzz_1[j] + fl1_fx * (tg_yzzz_xxxzz_0[j] - tg_yzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyzzz_xxyyy_0[j] = pb_y * tg_yyzzz_xxyyy_0[j] + fr * tg_yyzzz_xxyyy_1[j] + fl1_fx * (tg_yzzz_xxyyy_0[j] - tg_yzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzz_xxyy_1[j];

                    tg_yyyzzz_xxyyz_0[j] = pb_y * tg_yyzzz_xxyyz_0[j] + fr * tg_yyzzz_xxyyz_1[j] + fl1_fx * (tg_yzzz_xxyyz_0[j] - tg_yzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_xxyz_1[j];

                    tg_yyyzzz_xxyzz_0[j] = pb_y * tg_yyzzz_xxyzz_0[j] + fr * tg_yyzzz_xxyzz_1[j] + fl1_fx * (tg_yzzz_xxyzz_0[j] - tg_yzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xxzz_1[j];

                    tg_yyyzzz_xxzzz_0[j] = pb_y * tg_yyzzz_xxzzz_0[j] + fr * tg_yyzzz_xxzzz_1[j] + fl1_fx * (tg_yzzz_xxzzz_0[j] - tg_yzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyzzz_xyyyy_0[j] = pb_y * tg_yyzzz_xyyyy_0[j] + fr * tg_yyzzz_xyyyy_1[j] + fl1_fx * (tg_yzzz_xyyyy_0[j] - tg_yzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzz_xyyy_1[j];

                    tg_yyyzzz_xyyyz_0[j] = pb_y * tg_yyzzz_xyyyz_0[j] + fr * tg_yyzzz_xyyyz_1[j] + fl1_fx * (tg_yzzz_xyyyz_0[j] - tg_yzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzz_xyyz_1[j];

                    tg_yyyzzz_xyyzz_0[j] = pb_y * tg_yyzzz_xyyzz_0[j] + fr * tg_yyzzz_xyyzz_1[j] + fl1_fx * (tg_yzzz_xyyzz_0[j] - tg_yzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_xyzz_1[j];

                    tg_yyyzzz_xyzzz_0[j] = pb_y * tg_yyzzz_xyzzz_0[j] + fr * tg_yyzzz_xyzzz_1[j] + fl1_fx * (tg_yzzz_xyzzz_0[j] - tg_yzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_xzzz_1[j];

                    tg_yyyzzz_xzzzz_0[j] = pb_y * tg_yyzzz_xzzzz_0[j] + fr * tg_yyzzz_xzzzz_1[j] + fl1_fx * (tg_yzzz_xzzzz_0[j] - tg_yzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyzzz_yyyyy_0[j] = pb_y * tg_yyzzz_yyyyy_0[j] + fr * tg_yyzzz_yyyyy_1[j] + fl1_fx * (tg_yzzz_yyyyy_0[j] - tg_yzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzzz_yyyy_1[j];

                    tg_yyyzzz_yyyyz_0[j] = pb_y * tg_yyzzz_yyyyz_0[j] + fr * tg_yyzzz_yyyyz_1[j] + fl1_fx * (tg_yzzz_yyyyz_0[j] - tg_yzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzz_yyyz_1[j];

                    tg_yyyzzz_yyyzz_0[j] = pb_y * tg_yyzzz_yyyzz_0[j] + fr * tg_yyzzz_yyyzz_1[j] + fl1_fx * (tg_yzzz_yyyzz_0[j] - tg_yzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzz_yyzz_1[j];

                    tg_yyyzzz_yyzzz_0[j] = pb_y * tg_yyzzz_yyzzz_0[j] + fr * tg_yyzzz_yyzzz_1[j] + fl1_fx * (tg_yzzz_yyzzz_0[j] - tg_yzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzz_yzzz_1[j];

                    tg_yyyzzz_yzzzz_0[j] = pb_y * tg_yyzzz_yzzzz_0[j] + fr * tg_yyzzz_yzzzz_1[j] + fl1_fx * (tg_yzzz_yzzzz_0[j] - tg_yzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzz_zzzz_1[j];

                    tg_yyyzzz_zzzzz_0[j] = pb_y * tg_yyzzz_zzzzz_0[j] + fr * tg_yyzzz_zzzzz_1[j] + fl1_fx * (tg_yzzz_zzzzz_0[j] - tg_yzzz_zzzzz_1[j] * fl1_fza);

                    tg_yyzzzz_xxxxx_0[j] = pb_y * tg_yzzzz_xxxxx_0[j] + fr * tg_yzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxx_0[j] - tg_zzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyzzzz_xxxxy_0[j] = pb_y * tg_yzzzz_xxxxy_0[j] + fr * tg_yzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxy_0[j] - tg_zzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xxxx_1[j];

                    tg_yyzzzz_xxxxz_0[j] = pb_y * tg_yzzzz_xxxxz_0[j] + fr * tg_yzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxxz_0[j] - tg_zzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyzzzz_xxxyy_0[j] = pb_y * tg_yzzzz_xxxyy_0[j] + fr * tg_yzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyy_0[j] - tg_zzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_xxxy_1[j];

                    tg_yyzzzz_xxxyz_0[j] = pb_y * tg_yzzzz_xxxyz_0[j] + fr * tg_yzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxyz_0[j] - tg_zzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xxxz_1[j];

                    tg_yyzzzz_xxxzz_0[j] = pb_y * tg_yzzzz_xxxzz_0[j] + fr * tg_yzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxxzz_0[j] - tg_zzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyzzzz_xxyyy_0[j] = pb_y * tg_yzzzz_xxyyy_0[j] + fr * tg_yzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyy_0[j] - tg_zzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzz_xxyy_1[j];

                    tg_yyzzzz_xxyyz_0[j] = pb_y * tg_yzzzz_xxyyz_0[j] + fr * tg_yzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyyz_0[j] - tg_zzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_xxyz_1[j];

                    tg_yyzzzz_xxyzz_0[j] = pb_y * tg_yzzzz_xxyzz_0[j] + fr * tg_yzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxyzz_0[j] - tg_zzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xxzz_1[j];

                    tg_yyzzzz_xxzzz_0[j] = pb_y * tg_yzzzz_xxzzz_0[j] + fr * tg_yzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xxzzz_0[j] - tg_zzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyzzzz_xyyyy_0[j] = pb_y * tg_yzzzz_xyyyy_0[j] + fr * tg_yzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyy_0[j] - tg_zzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzz_xyyy_1[j];

                    tg_yyzzzz_xyyyz_0[j] = pb_y * tg_yzzzz_xyyyz_0[j] + fr * tg_yzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyyz_0[j] - tg_zzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzz_xyyz_1[j];

                    tg_yyzzzz_xyyzz_0[j] = pb_y * tg_yzzzz_xyyzz_0[j] + fr * tg_yzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyyzz_0[j] - tg_zzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_xyzz_1[j];

                    tg_yyzzzz_xyzzz_0[j] = pb_y * tg_yzzzz_xyzzz_0[j] + fr * tg_yzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xyzzz_0[j] - tg_zzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_xzzz_1[j];

                    tg_yyzzzz_xzzzz_0[j] = pb_y * tg_yzzzz_xzzzz_0[j] + fr * tg_yzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_xzzzz_0[j] - tg_zzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyzzzz_yyyyy_0[j] = pb_y * tg_yzzzz_yyyyy_0[j] + fr * tg_yzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyy_0[j] - tg_zzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzzz_yyyy_1[j];

                    tg_yyzzzz_yyyyz_0[j] = pb_y * tg_yzzzz_yyyyz_0[j] + fr * tg_yzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyyz_0[j] - tg_zzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzz_yyyz_1[j];

                    tg_yyzzzz_yyyzz_0[j] = pb_y * tg_yzzzz_yyyzz_0[j] + fr * tg_yzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyyzz_0[j] - tg_zzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzz_yyzz_1[j];

                    tg_yyzzzz_yyzzz_0[j] = pb_y * tg_yzzzz_yyzzz_0[j] + fr * tg_yzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yyzzz_0[j] - tg_zzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzz_yzzz_1[j];

                    tg_yyzzzz_yzzzz_0[j] = pb_y * tg_yzzzz_yzzzz_0[j] + fr * tg_yzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_yzzzz_0[j] - tg_zzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzz_zzzz_1[j];

                    tg_yyzzzz_zzzzz_0[j] = pb_y * tg_yzzzz_zzzzz_0[j] + fr * tg_yzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzz_zzzzz_0[j] - tg_zzzz_zzzzz_1[j] * fl1_fza);

                    tg_yzzzzz_xxxxx_0[j] = pb_y * tg_zzzzz_xxxxx_0[j] + fr * tg_zzzzz_xxxxx_1[j];

                    tg_yzzzzz_xxxxy_0[j] = pb_y * tg_zzzzz_xxxxy_0[j] + fr * tg_zzzzz_xxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxxx_1[j];

                    tg_yzzzzz_xxxxz_0[j] = pb_y * tg_zzzzz_xxxxz_0[j] + fr * tg_zzzzz_xxxxz_1[j];

                    tg_yzzzzz_xxxyy_0[j] = pb_y * tg_zzzzz_xxxyy_0[j] + fr * tg_zzzzz_xxxyy_1[j] + fl1_fxn * tg_zzzzz_xxxy_1[j];

                    tg_yzzzzz_xxxyz_0[j] = pb_y * tg_zzzzz_xxxyz_0[j] + fr * tg_zzzzz_xxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxxz_1[j];

                    tg_yzzzzz_xxxzz_0[j] = pb_y * tg_zzzzz_xxxzz_0[j] + fr * tg_zzzzz_xxxzz_1[j];

                    tg_yzzzzz_xxyyy_0[j] = pb_y * tg_zzzzz_xxyyy_0[j] + fr * tg_zzzzz_xxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xxyy_1[j];

                    tg_yzzzzz_xxyyz_0[j] = pb_y * tg_zzzzz_xxyyz_0[j] + fr * tg_zzzzz_xxyyz_1[j] + fl1_fxn * tg_zzzzz_xxyz_1[j];

                    tg_yzzzzz_xxyzz_0[j] = pb_y * tg_zzzzz_xxyzz_0[j] + fr * tg_zzzzz_xxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xxzz_1[j];

                    tg_yzzzzz_xxzzz_0[j] = pb_y * tg_zzzzz_xxzzz_0[j] + fr * tg_zzzzz_xxzzz_1[j];

                    tg_yzzzzz_xyyyy_0[j] = pb_y * tg_zzzzz_xyyyy_0[j] + fr * tg_zzzzz_xyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzz_xyyy_1[j];

                    tg_yzzzzz_xyyyz_0[j] = pb_y * tg_zzzzz_xyyyz_0[j] + fr * tg_zzzzz_xyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_xyyz_1[j];

                    tg_yzzzzz_xyyzz_0[j] = pb_y * tg_zzzzz_xyyzz_0[j] + fr * tg_zzzzz_xyyzz_1[j] + fl1_fxn * tg_zzzzz_xyzz_1[j];

                    tg_yzzzzz_xyzzz_0[j] = pb_y * tg_zzzzz_xyzzz_0[j] + fr * tg_zzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_xzzz_1[j];

                    tg_yzzzzz_xzzzz_0[j] = pb_y * tg_zzzzz_xzzzz_0[j] + fr * tg_zzzzz_xzzzz_1[j];

                    tg_yzzzzz_yyyyy_0[j] = pb_y * tg_zzzzz_yyyyy_0[j] + fr * tg_zzzzz_yyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzzz_yyyy_1[j];

                    tg_yzzzzz_yyyyz_0[j] = pb_y * tg_zzzzz_yyyyz_0[j] + fr * tg_zzzzz_yyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzzz_yyyz_1[j];

                    tg_yzzzzz_yyyzz_0[j] = pb_y * tg_zzzzz_yyyzz_0[j] + fr * tg_zzzzz_yyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzzz_yyzz_1[j];

                    tg_yzzzzz_yyzzz_0[j] = pb_y * tg_zzzzz_yyzzz_0[j] + fr * tg_zzzzz_yyzzz_1[j] + fl1_fxn * tg_zzzzz_yzzz_1[j];

                    tg_yzzzzz_yzzzz_0[j] = pb_y * tg_zzzzz_yzzzz_0[j] + fr * tg_zzzzz_yzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzz_zzzz_1[j];

                    tg_yzzzzz_zzzzz_0[j] = pb_y * tg_zzzzz_zzzzz_0[j] + fr * tg_zzzzz_zzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzz_xxxxx_0[j] = pb_z * tg_zzzzz_xxxxx_0[j] + fr * tg_zzzzz_xxxxx_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxx_0[j] - tg_zzzz_xxxxx_1[j] * fl1_fza);

                    tg_zzzzzz_xxxxy_0[j] = pb_z * tg_zzzzz_xxxxy_0[j] + fr * tg_zzzzz_xxxxy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxy_0[j] - tg_zzzz_xxxxy_1[j] * fl1_fza);

                    tg_zzzzzz_xxxxz_0[j] = pb_z * tg_zzzzz_xxxxz_0[j] + fr * tg_zzzzz_xxxxz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxxz_0[j] - tg_zzzz_xxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xxxx_1[j];

                    tg_zzzzzz_xxxyy_0[j] = pb_z * tg_zzzzz_xxxyy_0[j] + fr * tg_zzzzz_xxxyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxyy_0[j] - tg_zzzz_xxxyy_1[j] * fl1_fza);

                    tg_zzzzzz_xxxyz_0[j] = pb_z * tg_zzzzz_xxxyz_0[j] + fr * tg_zzzzz_xxxyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxyz_0[j] - tg_zzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xxxy_1[j];

                    tg_zzzzzz_xxxzz_0[j] = pb_z * tg_zzzzz_xxxzz_0[j] + fr * tg_zzzzz_xxxzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxxzz_0[j] - tg_zzzz_xxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_xxxz_1[j];

                    tg_zzzzzz_xxyyy_0[j] = pb_z * tg_zzzzz_xxyyy_0[j] + fr * tg_zzzzz_xxyyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxyyy_0[j] - tg_zzzz_xxyyy_1[j] * fl1_fza);

                    tg_zzzzzz_xxyyz_0[j] = pb_z * tg_zzzzz_xxyyz_0[j] + fr * tg_zzzzz_xxyyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxyyz_0[j] - tg_zzzz_xxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xxyy_1[j];

                    tg_zzzzzz_xxyzz_0[j] = pb_z * tg_zzzzz_xxyzz_0[j] + fr * tg_zzzzz_xxyzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxyzz_0[j] - tg_zzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_xxyz_1[j];

                    tg_zzzzzz_xxzzz_0[j] = pb_z * tg_zzzzz_xxzzz_0[j] + fr * tg_zzzzz_xxzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xxzzz_0[j] - tg_zzzz_xxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzz_xxzz_1[j];

                    tg_zzzzzz_xyyyy_0[j] = pb_z * tg_zzzzz_xyyyy_0[j] + fr * tg_zzzzz_xyyyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyyyy_0[j] - tg_zzzz_xyyyy_1[j] * fl1_fza);

                    tg_zzzzzz_xyyyz_0[j] = pb_z * tg_zzzzz_xyyyz_0[j] + fr * tg_zzzzz_xyyyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyyyz_0[j] - tg_zzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_xyyy_1[j];

                    tg_zzzzzz_xyyzz_0[j] = pb_z * tg_zzzzz_xyyzz_0[j] + fr * tg_zzzzz_xyyzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyyzz_0[j] - tg_zzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_xyyz_1[j];

                    tg_zzzzzz_xyzzz_0[j] = pb_z * tg_zzzzz_xyzzz_0[j] + fr * tg_zzzzz_xyzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xyzzz_0[j] - tg_zzzz_xyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzz_xyzz_1[j];

                    tg_zzzzzz_xzzzz_0[j] = pb_z * tg_zzzzz_xzzzz_0[j] + fr * tg_zzzzz_xzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_xzzzz_0[j] - tg_zzzz_xzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzz_xzzz_1[j];

                    tg_zzzzzz_yyyyy_0[j] = pb_z * tg_zzzzz_yyyyy_0[j] + fr * tg_zzzzz_yyyyy_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyyyy_0[j] - tg_zzzz_yyyyy_1[j] * fl1_fza);

                    tg_zzzzzz_yyyyz_0[j] = pb_z * tg_zzzzz_yyyyz_0[j] + fr * tg_zzzzz_yyyyz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyyyz_0[j] - tg_zzzz_yyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzz_yyyy_1[j];

                    tg_zzzzzz_yyyzz_0[j] = pb_z * tg_zzzzz_yyyzz_0[j] + fr * tg_zzzzz_yyyzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyyzz_0[j] - tg_zzzz_yyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzz_yyyz_1[j];

                    tg_zzzzzz_yyzzz_0[j] = pb_z * tg_zzzzz_yyzzz_0[j] + fr * tg_zzzzz_yyzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yyzzz_0[j] - tg_zzzz_yyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzz_yyzz_1[j];

                    tg_zzzzzz_yzzzz_0[j] = pb_z * tg_zzzzz_yzzzz_0[j] + fr * tg_zzzzz_yzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_yzzzz_0[j] - tg_zzzz_yzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzz_yzzz_1[j];

                    tg_zzzzzz_zzzzz_0[j] = pb_z * tg_zzzzz_zzzzz_0[j] + fr * tg_zzzzz_zzzzz_1[j] + 2.5 * fl1_fx * (tg_zzzz_zzzzz_0[j] - tg_zzzz_zzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

