//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForGH.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSGSH(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSGSH_0_79(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSGSH_79_158(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSGSH_158_237(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSH_237_315(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSGSH_0_79(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,79)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxx_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx); 

                auto tg_xxx_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 1); 

                auto tg_xxx_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 2); 

                auto tg_xxx_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 3); 

                auto tg_xxx_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 4); 

                auto tg_xxx_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 5); 

                auto tg_xxx_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 6); 

                auto tg_xxx_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 7); 

                auto tg_xxx_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 8); 

                auto tg_xxx_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 9); 

                auto tg_xxx_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 10); 

                auto tg_xxx_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 11); 

                auto tg_xxx_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 12); 

                auto tg_xxx_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 13); 

                auto tg_xxx_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 14); 

                auto tg_xxx_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 15); 

                auto tg_xxx_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 16); 

                auto tg_xxx_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 17); 

                auto tg_xxx_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 18); 

                auto tg_xxx_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 19); 

                auto tg_xxx_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 20); 

                auto tg_xxy_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 21); 

                auto tg_xxy_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 22); 

                auto tg_xxy_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 23); 

                auto tg_xxy_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 24); 

                auto tg_xxy_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 25); 

                auto tg_xxy_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 26); 

                auto tg_xxy_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 27); 

                auto tg_xxy_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 28); 

                auto tg_xxy_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 29); 

                auto tg_xxy_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 30); 

                auto tg_xxy_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 31); 

                auto tg_xxy_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 32); 

                auto tg_xxy_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 33); 

                auto tg_xxy_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 34); 

                auto tg_xxy_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 35); 

                auto tg_xxy_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 36); 

                auto tg_xxy_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 37); 

                auto tg_xxy_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 38); 

                auto tg_xxy_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 39); 

                auto tg_xxy_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 40); 

                auto tg_xxy_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 41); 

                auto tg_xxz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 42); 

                auto tg_xxz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 43); 

                auto tg_xxz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 44); 

                auto tg_xxz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 45); 

                auto tg_xxz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 46); 

                auto tg_xxz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 47); 

                auto tg_xxz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 48); 

                auto tg_xxz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 49); 

                auto tg_xxz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 50); 

                auto tg_xxz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 51); 

                auto tg_xxz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 52); 

                auto tg_xxz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 53); 

                auto tg_xxz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 54); 

                auto tg_xxz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 55); 

                auto tg_xxz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 56); 

                auto tg_xxz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 57); 

                auto tg_xxz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 58); 

                auto tg_xxz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 59); 

                auto tg_xxz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 60); 

                auto tg_xxz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 61); 

                auto tg_xxz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 62); 

                auto tg_xyy_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 63); 

                auto tg_xyy_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 64); 

                auto tg_xyy_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 65); 

                auto tg_xyy_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 66); 

                auto tg_xyy_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 67); 

                auto tg_xyy_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 68); 

                auto tg_xyy_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 69); 

                auto tg_xyy_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 70); 

                auto tg_xyy_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 71); 

                auto tg_xyy_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 72); 

                auto tg_xyy_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 73); 

                auto tg_xyy_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 74); 

                auto tg_xyy_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 75); 

                auto tg_xyy_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 76); 

                auto tg_xyy_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 77); 

                auto tg_xyy_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 78); 

                auto tg_xxx_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx); 

                auto tg_xxx_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 1); 

                auto tg_xxx_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 2); 

                auto tg_xxx_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 3); 

                auto tg_xxx_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 4); 

                auto tg_xxx_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 5); 

                auto tg_xxx_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 6); 

                auto tg_xxx_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 7); 

                auto tg_xxx_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 8); 

                auto tg_xxx_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 9); 

                auto tg_xxx_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 10); 

                auto tg_xxx_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 11); 

                auto tg_xxx_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 12); 

                auto tg_xxx_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 13); 

                auto tg_xxx_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 14); 

                auto tg_xxx_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 15); 

                auto tg_xxx_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 16); 

                auto tg_xxx_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 17); 

                auto tg_xxx_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 18); 

                auto tg_xxx_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 19); 

                auto tg_xxx_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 20); 

                auto tg_xxy_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 21); 

                auto tg_xxy_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 22); 

                auto tg_xxy_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 23); 

                auto tg_xxy_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 24); 

                auto tg_xxy_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 25); 

                auto tg_xxy_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 26); 

                auto tg_xxy_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 27); 

                auto tg_xxy_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 28); 

                auto tg_xxy_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 29); 

                auto tg_xxy_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 30); 

                auto tg_xxy_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 31); 

                auto tg_xxy_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 32); 

                auto tg_xxy_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 33); 

                auto tg_xxy_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 34); 

                auto tg_xxy_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 35); 

                auto tg_xxy_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 36); 

                auto tg_xxy_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 37); 

                auto tg_xxy_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 38); 

                auto tg_xxy_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 39); 

                auto tg_xxy_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 40); 

                auto tg_xxy_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 41); 

                auto tg_xxz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 42); 

                auto tg_xxz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 43); 

                auto tg_xxz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 44); 

                auto tg_xxz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 45); 

                auto tg_xxz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 46); 

                auto tg_xxz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 47); 

                auto tg_xxz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 48); 

                auto tg_xxz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 49); 

                auto tg_xxz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 50); 

                auto tg_xxz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 51); 

                auto tg_xxz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 52); 

                auto tg_xxz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 53); 

                auto tg_xxz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 54); 

                auto tg_xxz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 55); 

                auto tg_xxz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 56); 

                auto tg_xxz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 57); 

                auto tg_xxz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 58); 

                auto tg_xxz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 59); 

                auto tg_xxz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 60); 

                auto tg_xxz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 61); 

                auto tg_xxz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 62); 

                auto tg_xyy_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 63); 

                auto tg_xyy_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 64); 

                auto tg_xyy_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 65); 

                auto tg_xyy_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 66); 

                auto tg_xyy_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 67); 

                auto tg_xyy_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 68); 

                auto tg_xyy_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 69); 

                auto tg_xyy_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 70); 

                auto tg_xyy_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 71); 

                auto tg_xyy_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 72); 

                auto tg_xyy_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 73); 

                auto tg_xyy_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 74); 

                auto tg_xyy_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 75); 

                auto tg_xyy_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 76); 

                auto tg_xyy_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 77); 

                auto tg_xyy_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 78); 

                auto tg_xx_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx); 

                auto tg_xx_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 1); 

                auto tg_xx_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 2); 

                auto tg_xx_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 3); 

                auto tg_xx_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 4); 

                auto tg_xx_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 5); 

                auto tg_xx_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 6); 

                auto tg_xx_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 7); 

                auto tg_xx_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 8); 

                auto tg_xx_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 9); 

                auto tg_xx_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 10); 

                auto tg_xx_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 11); 

                auto tg_xx_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 12); 

                auto tg_xx_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 13); 

                auto tg_xx_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 14); 

                auto tg_xx_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 15); 

                auto tg_xx_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 16); 

                auto tg_xx_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 17); 

                auto tg_xx_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 18); 

                auto tg_xx_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 19); 

                auto tg_xx_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 20); 

                auto tg_xy_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 21); 

                auto tg_xy_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 22); 

                auto tg_xy_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 23); 

                auto tg_xy_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 24); 

                auto tg_xy_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 25); 

                auto tg_xy_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 26); 

                auto tg_xy_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 27); 

                auto tg_xy_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 28); 

                auto tg_xy_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 29); 

                auto tg_xy_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 30); 

                auto tg_xy_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 31); 

                auto tg_xy_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 32); 

                auto tg_xy_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 33); 

                auto tg_xy_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 34); 

                auto tg_xy_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 35); 

                auto tg_xy_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 36); 

                auto tg_xy_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 37); 

                auto tg_xy_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 38); 

                auto tg_xy_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 39); 

                auto tg_xy_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 40); 

                auto tg_xy_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 41); 

                auto tg_xz_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 42); 

                auto tg_xz_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 43); 

                auto tg_xz_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 44); 

                auto tg_xz_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 45); 

                auto tg_xz_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 46); 

                auto tg_xz_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 47); 

                auto tg_xz_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 48); 

                auto tg_xz_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 49); 

                auto tg_xz_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 50); 

                auto tg_xz_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 51); 

                auto tg_xz_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 52); 

                auto tg_xz_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 53); 

                auto tg_xz_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 54); 

                auto tg_xz_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 55); 

                auto tg_xz_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 56); 

                auto tg_xz_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 57); 

                auto tg_xz_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 58); 

                auto tg_xz_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 59); 

                auto tg_xz_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 60); 

                auto tg_xz_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 61); 

                auto tg_xz_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 62); 

                auto tg_yy_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 63); 

                auto tg_yy_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 64); 

                auto tg_yy_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 65); 

                auto tg_yy_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 66); 

                auto tg_yy_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 67); 

                auto tg_yy_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 68); 

                auto tg_yy_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 69); 

                auto tg_yy_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 70); 

                auto tg_yy_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 71); 

                auto tg_yy_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 72); 

                auto tg_yy_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 73); 

                auto tg_yy_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 74); 

                auto tg_yy_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 75); 

                auto tg_yy_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 76); 

                auto tg_yy_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 77); 

                auto tg_yy_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 78); 

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

                auto tg_yy_xyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 73); 

                auto tg_yy_xyyyz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 74); 

                auto tg_yy_xyyzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 75); 

                auto tg_yy_xyzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 76); 

                auto tg_yy_xzzzz_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 77); 

                auto tg_yy_yyyyy_1 = primBuffer[pidx_g_2_5_m1].data(126 * idx + 78); 

                auto tg_xxx_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx); 

                auto tg_xxx_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 1); 

                auto tg_xxx_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 2); 

                auto tg_xxx_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 3); 

                auto tg_xxx_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 4); 

                auto tg_xxx_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 5); 

                auto tg_xxx_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 6); 

                auto tg_xxx_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 7); 

                auto tg_xxx_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 8); 

                auto tg_xxx_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 9); 

                auto tg_xxx_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 10); 

                auto tg_xxx_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 11); 

                auto tg_xxx_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 12); 

                auto tg_xxx_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 13); 

                auto tg_xxx_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 14); 

                auto tg_xxy_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 15); 

                auto tg_xxy_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 16); 

                auto tg_xxy_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 17); 

                auto tg_xxy_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 18); 

                auto tg_xxy_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 19); 

                auto tg_xxy_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 20); 

                auto tg_xxy_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 21); 

                auto tg_xxy_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 22); 

                auto tg_xxy_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 23); 

                auto tg_xxy_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 24); 

                auto tg_xxy_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 25); 

                auto tg_xxy_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 26); 

                auto tg_xxy_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 27); 

                auto tg_xxy_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 28); 

                auto tg_xxy_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 29); 

                auto tg_xxz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 30); 

                auto tg_xxz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 31); 

                auto tg_xxz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 32); 

                auto tg_xxz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 33); 

                auto tg_xxz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 34); 

                auto tg_xxz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 35); 

                auto tg_xxz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 36); 

                auto tg_xxz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 37); 

                auto tg_xxz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 38); 

                auto tg_xxz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 39); 

                auto tg_xxz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 40); 

                auto tg_xxz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 41); 

                auto tg_xxz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 42); 

                auto tg_xxz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 43); 

                auto tg_xxz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 44); 

                auto tg_xyy_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 45); 

                auto tg_xyy_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 46); 

                auto tg_xyy_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 47); 

                auto tg_xyy_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 48); 

                auto tg_xyy_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 49); 

                auto tg_xyy_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 50); 

                auto tg_xyy_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 51); 

                auto tg_xyy_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 52); 

                auto tg_xyy_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 53); 

                auto tg_xyy_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 54); 

                auto tg_xyy_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 55); 

                auto tg_xyy_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 56); 

                auto tg_xyy_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 57); 

                auto tg_xyy_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 58); 

                auto tg_xyy_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 59); 

                // set up pointers to integrals

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

                // Batch of Integrals (0,79)

                #pragma omp simd aligned(fxn, fza, tg_xx_xxxxx_0, tg_xx_xxxxx_1, tg_xx_xxxxy_0, tg_xx_xxxxy_1, \
                                         tg_xx_xxxxz_0, tg_xx_xxxxz_1, tg_xx_xxxyy_0, tg_xx_xxxyy_1, tg_xx_xxxyz_0, \
                                         tg_xx_xxxyz_1, tg_xx_xxxzz_0, tg_xx_xxxzz_1, tg_xx_xxyyy_0, tg_xx_xxyyy_1, \
                                         tg_xx_xxyyz_0, tg_xx_xxyyz_1, tg_xx_xxyzz_0, tg_xx_xxyzz_1, tg_xx_xxzzz_0, \
                                         tg_xx_xxzzz_1, tg_xx_xyyyy_0, tg_xx_xyyyy_1, tg_xx_xyyyz_0, tg_xx_xyyyz_1, \
                                         tg_xx_xyyzz_0, tg_xx_xyyzz_1, tg_xx_xyzzz_0, tg_xx_xyzzz_1, tg_xx_xzzzz_0, \
                                         tg_xx_xzzzz_1, tg_xx_yyyyy_0, tg_xx_yyyyy_1, tg_xx_yyyyz_0, tg_xx_yyyyz_1, \
                                         tg_xx_yyyzz_0, tg_xx_yyyzz_1, tg_xx_yyzzz_0, tg_xx_yyzzz_1, tg_xx_yzzzz_0, \
                                         tg_xx_yzzzz_1, tg_xx_zzzzz_0, tg_xx_zzzzz_1, tg_xxx_xxxx_1, tg_xxx_xxxxx_0, \
                                         tg_xxx_xxxxx_1, tg_xxx_xxxxy_0, tg_xxx_xxxxy_1, tg_xxx_xxxxz_0, tg_xxx_xxxxz_1, \
                                         tg_xxx_xxxy_1, tg_xxx_xxxyy_0, tg_xxx_xxxyy_1, tg_xxx_xxxyz_0, tg_xxx_xxxyz_1, \
                                         tg_xxx_xxxz_1, tg_xxx_xxxzz_0, tg_xxx_xxxzz_1, tg_xxx_xxyy_1, tg_xxx_xxyyy_0, \
                                         tg_xxx_xxyyy_1, tg_xxx_xxyyz_0, tg_xxx_xxyyz_1, tg_xxx_xxyz_1, tg_xxx_xxyzz_0, \
                                         tg_xxx_xxyzz_1, tg_xxx_xxzz_1, tg_xxx_xxzzz_0, tg_xxx_xxzzz_1, tg_xxx_xyyy_1, \
                                         tg_xxx_xyyyy_0, tg_xxx_xyyyy_1, tg_xxx_xyyyz_0, tg_xxx_xyyyz_1, tg_xxx_xyyz_1, \
                                         tg_xxx_xyyzz_0, tg_xxx_xyyzz_1, tg_xxx_xyzz_1, tg_xxx_xyzzz_0, tg_xxx_xyzzz_1, \
                                         tg_xxx_xzzz_1, tg_xxx_xzzzz_0, tg_xxx_xzzzz_1, tg_xxx_yyyy_1, tg_xxx_yyyyy_0, \
                                         tg_xxx_yyyyy_1, tg_xxx_yyyyz_0, tg_xxx_yyyyz_1, tg_xxx_yyyz_1, tg_xxx_yyyzz_0, \
                                         tg_xxx_yyyzz_1, tg_xxx_yyzz_1, tg_xxx_yyzzz_0, tg_xxx_yyzzz_1, tg_xxx_yzzz_1, \
                                         tg_xxx_yzzzz_0, tg_xxx_yzzzz_1, tg_xxx_zzzz_1, tg_xxx_zzzzz_0, tg_xxx_zzzzz_1, \
                                         tg_xxxx_xxxxx_0, tg_xxxx_xxxxy_0, tg_xxxx_xxxxz_0, tg_xxxx_xxxyy_0, tg_xxxx_xxxyz_0, \
                                         tg_xxxx_xxxzz_0, tg_xxxx_xxyyy_0, tg_xxxx_xxyyz_0, tg_xxxx_xxyzz_0, tg_xxxx_xxzzz_0, \
                                         tg_xxxx_xyyyy_0, tg_xxxx_xyyyz_0, tg_xxxx_xyyzz_0, tg_xxxx_xyzzz_0, tg_xxxx_xzzzz_0, \
                                         tg_xxxx_yyyyy_0, tg_xxxx_yyyyz_0, tg_xxxx_yyyzz_0, tg_xxxx_yyzzz_0, tg_xxxx_yzzzz_0, \
                                         tg_xxxx_zzzzz_0, tg_xxxy_xxxxx_0, tg_xxxy_xxxxy_0, tg_xxxy_xxxxz_0, tg_xxxy_xxxyy_0, \
                                         tg_xxxy_xxxyz_0, tg_xxxy_xxxzz_0, tg_xxxy_xxyyy_0, tg_xxxy_xxyyz_0, tg_xxxy_xxyzz_0, \
                                         tg_xxxy_xxzzz_0, tg_xxxy_xyyyy_0, tg_xxxy_xyyyz_0, tg_xxxy_xyyzz_0, tg_xxxy_xyzzz_0, \
                                         tg_xxxy_xzzzz_0, tg_xxxy_yyyyy_0, tg_xxxy_yyyyz_0, tg_xxxy_yyyzz_0, tg_xxxy_yyzzz_0, \
                                         tg_xxxy_yzzzz_0, tg_xxxy_zzzzz_0, tg_xxxz_xxxxx_0, tg_xxxz_xxxxy_0, tg_xxxz_xxxxz_0, \
                                         tg_xxxz_xxxyy_0, tg_xxxz_xxxyz_0, tg_xxxz_xxxzz_0, tg_xxxz_xxyyy_0, tg_xxxz_xxyyz_0, \
                                         tg_xxxz_xxyzz_0, tg_xxxz_xxzzz_0, tg_xxxz_xyyyy_0, tg_xxxz_xyyyz_0, tg_xxxz_xyyzz_0, \
                                         tg_xxxz_xyzzz_0, tg_xxxz_xzzzz_0, tg_xxxz_yyyyy_0, tg_xxxz_yyyyz_0, tg_xxxz_yyyzz_0, \
                                         tg_xxxz_yyzzz_0, tg_xxxz_yzzzz_0, tg_xxxz_zzzzz_0, tg_xxy_xxxx_1, tg_xxy_xxxxx_0, \
                                         tg_xxy_xxxxx_1, tg_xxy_xxxxy_0, tg_xxy_xxxxy_1, tg_xxy_xxxxz_0, tg_xxy_xxxxz_1, \
                                         tg_xxy_xxxy_1, tg_xxy_xxxyy_0, tg_xxy_xxxyy_1, tg_xxy_xxxyz_0, tg_xxy_xxxyz_1, \
                                         tg_xxy_xxxz_1, tg_xxy_xxxzz_0, tg_xxy_xxxzz_1, tg_xxy_xxyy_1, tg_xxy_xxyyy_0, \
                                         tg_xxy_xxyyy_1, tg_xxy_xxyyz_0, tg_xxy_xxyyz_1, tg_xxy_xxyz_1, tg_xxy_xxyzz_0, \
                                         tg_xxy_xxyzz_1, tg_xxy_xxzz_1, tg_xxy_xxzzz_0, tg_xxy_xxzzz_1, tg_xxy_xyyy_1, \
                                         tg_xxy_xyyyy_0, tg_xxy_xyyyy_1, tg_xxy_xyyyz_0, tg_xxy_xyyyz_1, tg_xxy_xyyz_1, \
                                         tg_xxy_xyyzz_0, tg_xxy_xyyzz_1, tg_xxy_xyzz_1, tg_xxy_xyzzz_0, tg_xxy_xyzzz_1, \
                                         tg_xxy_xzzz_1, tg_xxy_xzzzz_0, tg_xxy_xzzzz_1, tg_xxy_yyyy_1, tg_xxy_yyyyy_0, \
                                         tg_xxy_yyyyy_1, tg_xxy_yyyyz_0, tg_xxy_yyyyz_1, tg_xxy_yyyz_1, tg_xxy_yyyzz_0, \
                                         tg_xxy_yyyzz_1, tg_xxy_yyzz_1, tg_xxy_yyzzz_0, tg_xxy_yyzzz_1, tg_xxy_yzzz_1, \
                                         tg_xxy_yzzzz_0, tg_xxy_yzzzz_1, tg_xxy_zzzz_1, tg_xxy_zzzzz_0, tg_xxy_zzzzz_1, \
                                         tg_xxyy_xxxxx_0, tg_xxyy_xxxxy_0, tg_xxyy_xxxxz_0, tg_xxyy_xxxyy_0, tg_xxyy_xxxyz_0, \
                                         tg_xxyy_xxxzz_0, tg_xxyy_xxyyy_0, tg_xxyy_xxyyz_0, tg_xxyy_xxyzz_0, tg_xxyy_xxzzz_0, \
                                         tg_xxyy_xyyyy_0, tg_xxyy_xyyyz_0, tg_xxyy_xyyzz_0, tg_xxyy_xyzzz_0, tg_xxyy_xzzzz_0, \
                                         tg_xxyy_yyyyy_0, tg_xxz_xxxx_1, tg_xxz_xxxxx_0, tg_xxz_xxxxx_1, tg_xxz_xxxxy_0, \
                                         tg_xxz_xxxxy_1, tg_xxz_xxxxz_0, tg_xxz_xxxxz_1, tg_xxz_xxxy_1, tg_xxz_xxxyy_0, \
                                         tg_xxz_xxxyy_1, tg_xxz_xxxyz_0, tg_xxz_xxxyz_1, tg_xxz_xxxz_1, tg_xxz_xxxzz_0, \
                                         tg_xxz_xxxzz_1, tg_xxz_xxyy_1, tg_xxz_xxyyy_0, tg_xxz_xxyyy_1, tg_xxz_xxyyz_0, \
                                         tg_xxz_xxyyz_1, tg_xxz_xxyz_1, tg_xxz_xxyzz_0, tg_xxz_xxyzz_1, tg_xxz_xxzz_1, \
                                         tg_xxz_xxzzz_0, tg_xxz_xxzzz_1, tg_xxz_xyyy_1, tg_xxz_xyyyy_0, tg_xxz_xyyyy_1, \
                                         tg_xxz_xyyyz_0, tg_xxz_xyyyz_1, tg_xxz_xyyz_1, tg_xxz_xyyzz_0, tg_xxz_xyyzz_1, \
                                         tg_xxz_xyzz_1, tg_xxz_xyzzz_0, tg_xxz_xyzzz_1, tg_xxz_xzzz_1, tg_xxz_xzzzz_0, \
                                         tg_xxz_xzzzz_1, tg_xxz_yyyy_1, tg_xxz_yyyyy_0, tg_xxz_yyyyy_1, tg_xxz_yyyyz_0, \
                                         tg_xxz_yyyyz_1, tg_xxz_yyyz_1, tg_xxz_yyyzz_0, tg_xxz_yyyzz_1, tg_xxz_yyzz_1, \
                                         tg_xxz_yyzzz_0, tg_xxz_yyzzz_1, tg_xxz_yzzz_1, tg_xxz_yzzzz_0, tg_xxz_yzzzz_1, \
                                         tg_xxz_zzzz_1, tg_xxz_zzzzz_0, tg_xxz_zzzzz_1, tg_xy_xxxxx_0, tg_xy_xxxxx_1, \
                                         tg_xy_xxxxy_0, tg_xy_xxxxy_1, tg_xy_xxxxz_0, tg_xy_xxxxz_1, tg_xy_xxxyy_0, \
                                         tg_xy_xxxyy_1, tg_xy_xxxyz_0, tg_xy_xxxyz_1, tg_xy_xxxzz_0, tg_xy_xxxzz_1, \
                                         tg_xy_xxyyy_0, tg_xy_xxyyy_1, tg_xy_xxyyz_0, tg_xy_xxyyz_1, tg_xy_xxyzz_0, \
                                         tg_xy_xxyzz_1, tg_xy_xxzzz_0, tg_xy_xxzzz_1, tg_xy_xyyyy_0, tg_xy_xyyyy_1, \
                                         tg_xy_xyyyz_0, tg_xy_xyyyz_1, tg_xy_xyyzz_0, tg_xy_xyyzz_1, tg_xy_xyzzz_0, \
                                         tg_xy_xyzzz_1, tg_xy_xzzzz_0, tg_xy_xzzzz_1, tg_xy_yyyyy_0, tg_xy_yyyyy_1, \
                                         tg_xy_yyyyz_0, tg_xy_yyyyz_1, tg_xy_yyyzz_0, tg_xy_yyyzz_1, tg_xy_yyzzz_0, \
                                         tg_xy_yyzzz_1, tg_xy_yzzzz_0, tg_xy_yzzzz_1, tg_xy_zzzzz_0, tg_xy_zzzzz_1, \
                                         tg_xyy_xxxx_1, tg_xyy_xxxxx_0, tg_xyy_xxxxx_1, tg_xyy_xxxxy_0, tg_xyy_xxxxy_1, \
                                         tg_xyy_xxxxz_0, tg_xyy_xxxxz_1, tg_xyy_xxxy_1, tg_xyy_xxxyy_0, tg_xyy_xxxyy_1, \
                                         tg_xyy_xxxyz_0, tg_xyy_xxxyz_1, tg_xyy_xxxz_1, tg_xyy_xxxzz_0, tg_xyy_xxxzz_1, \
                                         tg_xyy_xxyy_1, tg_xyy_xxyyy_0, tg_xyy_xxyyy_1, tg_xyy_xxyyz_0, tg_xyy_xxyyz_1, \
                                         tg_xyy_xxyz_1, tg_xyy_xxyzz_0, tg_xyy_xxyzz_1, tg_xyy_xxzz_1, tg_xyy_xxzzz_0, \
                                         tg_xyy_xxzzz_1, tg_xyy_xyyy_1, tg_xyy_xyyyy_0, tg_xyy_xyyyy_1, tg_xyy_xyyyz_0, \
                                         tg_xyy_xyyyz_1, tg_xyy_xyyz_1, tg_xyy_xyyzz_0, tg_xyy_xyyzz_1, tg_xyy_xyzz_1, \
                                         tg_xyy_xyzzz_0, tg_xyy_xyzzz_1, tg_xyy_xzzz_1, tg_xyy_xzzzz_0, tg_xyy_xzzzz_1, \
                                         tg_xyy_yyyy_1, tg_xyy_yyyyy_0, tg_xyy_yyyyy_1, tg_xyy_yyyz_1, tg_xyy_yyzz_1, \
                                         tg_xyy_yzzz_1, tg_xyy_zzzz_1, tg_xz_xxxxx_0, tg_xz_xxxxx_1, tg_xz_xxxxy_0, \
                                         tg_xz_xxxxy_1, tg_xz_xxxxz_0, tg_xz_xxxxz_1, tg_xz_xxxyy_0, tg_xz_xxxyy_1, \
                                         tg_xz_xxxyz_0, tg_xz_xxxyz_1, tg_xz_xxxzz_0, tg_xz_xxxzz_1, tg_xz_xxyyy_0, \
                                         tg_xz_xxyyy_1, tg_xz_xxyyz_0, tg_xz_xxyyz_1, tg_xz_xxyzz_0, tg_xz_xxyzz_1, \
                                         tg_xz_xxzzz_0, tg_xz_xxzzz_1, tg_xz_xyyyy_0, tg_xz_xyyyy_1, tg_xz_xyyyz_0, \
                                         tg_xz_xyyyz_1, tg_xz_xyyzz_0, tg_xz_xyyzz_1, tg_xz_xyzzz_0, tg_xz_xyzzz_1, \
                                         tg_xz_xzzzz_0, tg_xz_xzzzz_1, tg_xz_yyyyy_0, tg_xz_yyyyy_1, tg_xz_yyyyz_0, \
                                         tg_xz_yyyyz_1, tg_xz_yyyzz_0, tg_xz_yyyzz_1, tg_xz_yyzzz_0, tg_xz_yyzzz_1, \
                                         tg_xz_yzzzz_0, tg_xz_yzzzz_1, tg_xz_zzzzz_0, tg_xz_zzzzz_1, tg_yy_xxxxx_0, \
                                         tg_yy_xxxxx_1, tg_yy_xxxxy_0, tg_yy_xxxxy_1, tg_yy_xxxxz_0, tg_yy_xxxxz_1, \
                                         tg_yy_xxxyy_0, tg_yy_xxxyy_1, tg_yy_xxxyz_0, tg_yy_xxxyz_1, tg_yy_xxxzz_0, \
                                         tg_yy_xxxzz_1, tg_yy_xxyyy_0, tg_yy_xxyyy_1, tg_yy_xxyyz_0, tg_yy_xxyyz_1, \
                                         tg_yy_xxyzz_0, tg_yy_xxyzz_1, tg_yy_xxzzz_0, tg_yy_xxzzz_1, tg_yy_xyyyy_0, \
                                         tg_yy_xyyyy_1, tg_yy_xyyyz_0, tg_yy_xyyyz_1, tg_yy_xyyzz_0, tg_yy_xyyzz_1, \
                                         tg_yy_xyzzz_0, tg_yy_xyzzz_1, tg_yy_xzzzz_0, tg_yy_xzzzz_1, tg_yy_yyyyy_0, \
                                         tg_yy_yyyyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxx_xxxxx_0[j] = pb_x * tg_xxx_xxxxx_0[j] + fr * tg_xxx_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxx_0[j] - tg_xx_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxx_xxxx_1[j];

                    tg_xxxx_xxxxy_0[j] = pb_x * tg_xxx_xxxxy_0[j] + fr * tg_xxx_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxy_0[j] - tg_xx_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxy_1[j];

                    tg_xxxx_xxxxz_0[j] = pb_x * tg_xxx_xxxxz_0[j] + fr * tg_xxx_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxz_0[j] - tg_xx_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxz_1[j];

                    tg_xxxx_xxxyy_0[j] = pb_x * tg_xxx_xxxyy_0[j] + fr * tg_xxx_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyy_0[j] - tg_xx_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyy_1[j];

                    tg_xxxx_xxxyz_0[j] = pb_x * tg_xxx_xxxyz_0[j] + fr * tg_xxx_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyz_0[j] - tg_xx_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyz_1[j];

                    tg_xxxx_xxxzz_0[j] = pb_x * tg_xxx_xxxzz_0[j] + fr * tg_xxx_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxzz_0[j] - tg_xx_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxzz_1[j];

                    tg_xxxx_xxyyy_0[j] = pb_x * tg_xxx_xxyyy_0[j] + fr * tg_xxx_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyy_0[j] - tg_xx_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyy_1[j];

                    tg_xxxx_xxyyz_0[j] = pb_x * tg_xxx_xxyyz_0[j] + fr * tg_xxx_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyz_0[j] - tg_xx_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyz_1[j];

                    tg_xxxx_xxyzz_0[j] = pb_x * tg_xxx_xxyzz_0[j] + fr * tg_xxx_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyzz_0[j] - tg_xx_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyzz_1[j];

                    tg_xxxx_xxzzz_0[j] = pb_x * tg_xxx_xxzzz_0[j] + fr * tg_xxx_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxzzz_0[j] - tg_xx_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xzzz_1[j];

                    tg_xxxx_xyyyy_0[j] = pb_x * tg_xxx_xyyyy_0[j] + fr * tg_xxx_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyy_0[j] - tg_xx_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyy_1[j];

                    tg_xxxx_xyyyz_0[j] = pb_x * tg_xxx_xyyyz_0[j] + fr * tg_xxx_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyz_0[j] - tg_xx_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyz_1[j];

                    tg_xxxx_xyyzz_0[j] = pb_x * tg_xxx_xyyzz_0[j] + fr * tg_xxx_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyzz_0[j] - tg_xx_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyzz_1[j];

                    tg_xxxx_xyzzz_0[j] = pb_x * tg_xxx_xyzzz_0[j] + fr * tg_xxx_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyzzz_0[j] - tg_xx_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yzzz_1[j];

                    tg_xxxx_xzzzz_0[j] = pb_x * tg_xxx_xzzzz_0[j] + fr * tg_xxx_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xzzzz_0[j] - tg_xx_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_zzzz_1[j];

                    tg_xxxx_yyyyy_0[j] = pb_x * tg_xxx_yyyyy_0[j] + fr * tg_xxx_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyy_0[j] - tg_xx_yyyyy_1[j] * fl1_fza);

                    tg_xxxx_yyyyz_0[j] = pb_x * tg_xxx_yyyyz_0[j] + fr * tg_xxx_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyz_0[j] - tg_xx_yyyyz_1[j] * fl1_fza);

                    tg_xxxx_yyyzz_0[j] = pb_x * tg_xxx_yyyzz_0[j] + fr * tg_xxx_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyzz_0[j] - tg_xx_yyyzz_1[j] * fl1_fza);

                    tg_xxxx_yyzzz_0[j] = pb_x * tg_xxx_yyzzz_0[j] + fr * tg_xxx_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyzzz_0[j] - tg_xx_yyzzz_1[j] * fl1_fza);

                    tg_xxxx_yzzzz_0[j] = pb_x * tg_xxx_yzzzz_0[j] + fr * tg_xxx_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yzzzz_0[j] - tg_xx_yzzzz_1[j] * fl1_fza);

                    tg_xxxx_zzzzz_0[j] = pb_x * tg_xxx_zzzzz_0[j] + fr * tg_xxx_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_zzzzz_0[j] - tg_xx_zzzzz_1[j] * fl1_fza);

                    tg_xxxy_xxxxx_0[j] = pb_x * tg_xxy_xxxxx_0[j] + fr * tg_xxy_xxxxx_1[j] + fl1_fx * (tg_xy_xxxxx_0[j] - tg_xy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxy_xxxx_1[j];

                    tg_xxxy_xxxxy_0[j] = pb_x * tg_xxy_xxxxy_0[j] + fr * tg_xxy_xxxxy_1[j] + fl1_fx * (tg_xy_xxxxy_0[j] - tg_xy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxy_1[j];

                    tg_xxxy_xxxxz_0[j] = pb_x * tg_xxy_xxxxz_0[j] + fr * tg_xxy_xxxxz_1[j] + fl1_fx * (tg_xy_xxxxz_0[j] - tg_xy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxz_1[j];

                    tg_xxxy_xxxyy_0[j] = pb_x * tg_xxy_xxxyy_0[j] + fr * tg_xxy_xxxyy_1[j] + fl1_fx * (tg_xy_xxxyy_0[j] - tg_xy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyy_1[j];

                    tg_xxxy_xxxyz_0[j] = pb_x * tg_xxy_xxxyz_0[j] + fr * tg_xxy_xxxyz_1[j] + fl1_fx * (tg_xy_xxxyz_0[j] - tg_xy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyz_1[j];

                    tg_xxxy_xxxzz_0[j] = pb_x * tg_xxy_xxxzz_0[j] + fr * tg_xxy_xxxzz_1[j] + fl1_fx * (tg_xy_xxxzz_0[j] - tg_xy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxzz_1[j];

                    tg_xxxy_xxyyy_0[j] = pb_x * tg_xxy_xxyyy_0[j] + fr * tg_xxy_xxyyy_1[j] + fl1_fx * (tg_xy_xxyyy_0[j] - tg_xy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyy_1[j];

                    tg_xxxy_xxyyz_0[j] = pb_x * tg_xxy_xxyyz_0[j] + fr * tg_xxy_xxyyz_1[j] + fl1_fx * (tg_xy_xxyyz_0[j] - tg_xy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyz_1[j];

                    tg_xxxy_xxyzz_0[j] = pb_x * tg_xxy_xxyzz_0[j] + fr * tg_xxy_xxyzz_1[j] + fl1_fx * (tg_xy_xxyzz_0[j] - tg_xy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyzz_1[j];

                    tg_xxxy_xxzzz_0[j] = pb_x * tg_xxy_xxzzz_0[j] + fr * tg_xxy_xxzzz_1[j] + fl1_fx * (tg_xy_xxzzz_0[j] - tg_xy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xzzz_1[j];

                    tg_xxxy_xyyyy_0[j] = pb_x * tg_xxy_xyyyy_0[j] + fr * tg_xxy_xyyyy_1[j] + fl1_fx * (tg_xy_xyyyy_0[j] - tg_xy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyy_1[j];

                    tg_xxxy_xyyyz_0[j] = pb_x * tg_xxy_xyyyz_0[j] + fr * tg_xxy_xyyyz_1[j] + fl1_fx * (tg_xy_xyyyz_0[j] - tg_xy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyz_1[j];

                    tg_xxxy_xyyzz_0[j] = pb_x * tg_xxy_xyyzz_0[j] + fr * tg_xxy_xyyzz_1[j] + fl1_fx * (tg_xy_xyyzz_0[j] - tg_xy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyzz_1[j];

                    tg_xxxy_xyzzz_0[j] = pb_x * tg_xxy_xyzzz_0[j] + fr * tg_xxy_xyzzz_1[j] + fl1_fx * (tg_xy_xyzzz_0[j] - tg_xy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yzzz_1[j];

                    tg_xxxy_xzzzz_0[j] = pb_x * tg_xxy_xzzzz_0[j] + fr * tg_xxy_xzzzz_1[j] + fl1_fx * (tg_xy_xzzzz_0[j] - tg_xy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_zzzz_1[j];

                    tg_xxxy_yyyyy_0[j] = pb_x * tg_xxy_yyyyy_0[j] + fr * tg_xxy_yyyyy_1[j] + fl1_fx * (tg_xy_yyyyy_0[j] - tg_xy_yyyyy_1[j] * fl1_fza);

                    tg_xxxy_yyyyz_0[j] = pb_x * tg_xxy_yyyyz_0[j] + fr * tg_xxy_yyyyz_1[j] + fl1_fx * (tg_xy_yyyyz_0[j] - tg_xy_yyyyz_1[j] * fl1_fza);

                    tg_xxxy_yyyzz_0[j] = pb_x * tg_xxy_yyyzz_0[j] + fr * tg_xxy_yyyzz_1[j] + fl1_fx * (tg_xy_yyyzz_0[j] - tg_xy_yyyzz_1[j] * fl1_fza);

                    tg_xxxy_yyzzz_0[j] = pb_x * tg_xxy_yyzzz_0[j] + fr * tg_xxy_yyzzz_1[j] + fl1_fx * (tg_xy_yyzzz_0[j] - tg_xy_yyzzz_1[j] * fl1_fza);

                    tg_xxxy_yzzzz_0[j] = pb_x * tg_xxy_yzzzz_0[j] + fr * tg_xxy_yzzzz_1[j] + fl1_fx * (tg_xy_yzzzz_0[j] - tg_xy_yzzzz_1[j] * fl1_fza);

                    tg_xxxy_zzzzz_0[j] = pb_x * tg_xxy_zzzzz_0[j] + fr * tg_xxy_zzzzz_1[j] + fl1_fx * (tg_xy_zzzzz_0[j] - tg_xy_zzzzz_1[j] * fl1_fza);

                    tg_xxxz_xxxxx_0[j] = pb_x * tg_xxz_xxxxx_0[j] + fr * tg_xxz_xxxxx_1[j] + fl1_fx * (tg_xz_xxxxx_0[j] - tg_xz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxz_xxxx_1[j];

                    tg_xxxz_xxxxy_0[j] = pb_x * tg_xxz_xxxxy_0[j] + fr * tg_xxz_xxxxy_1[j] + fl1_fx * (tg_xz_xxxxy_0[j] - tg_xz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxy_1[j];

                    tg_xxxz_xxxxz_0[j] = pb_x * tg_xxz_xxxxz_0[j] + fr * tg_xxz_xxxxz_1[j] + fl1_fx * (tg_xz_xxxxz_0[j] - tg_xz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxz_1[j];

                    tg_xxxz_xxxyy_0[j] = pb_x * tg_xxz_xxxyy_0[j] + fr * tg_xxz_xxxyy_1[j] + fl1_fx * (tg_xz_xxxyy_0[j] - tg_xz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyy_1[j];

                    tg_xxxz_xxxyz_0[j] = pb_x * tg_xxz_xxxyz_0[j] + fr * tg_xxz_xxxyz_1[j] + fl1_fx * (tg_xz_xxxyz_0[j] - tg_xz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyz_1[j];

                    tg_xxxz_xxxzz_0[j] = pb_x * tg_xxz_xxxzz_0[j] + fr * tg_xxz_xxxzz_1[j] + fl1_fx * (tg_xz_xxxzz_0[j] - tg_xz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxzz_1[j];

                    tg_xxxz_xxyyy_0[j] = pb_x * tg_xxz_xxyyy_0[j] + fr * tg_xxz_xxyyy_1[j] + fl1_fx * (tg_xz_xxyyy_0[j] - tg_xz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyy_1[j];

                    tg_xxxz_xxyyz_0[j] = pb_x * tg_xxz_xxyyz_0[j] + fr * tg_xxz_xxyyz_1[j] + fl1_fx * (tg_xz_xxyyz_0[j] - tg_xz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyz_1[j];

                    tg_xxxz_xxyzz_0[j] = pb_x * tg_xxz_xxyzz_0[j] + fr * tg_xxz_xxyzz_1[j] + fl1_fx * (tg_xz_xxyzz_0[j] - tg_xz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyzz_1[j];

                    tg_xxxz_xxzzz_0[j] = pb_x * tg_xxz_xxzzz_0[j] + fr * tg_xxz_xxzzz_1[j] + fl1_fx * (tg_xz_xxzzz_0[j] - tg_xz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xzzz_1[j];

                    tg_xxxz_xyyyy_0[j] = pb_x * tg_xxz_xyyyy_0[j] + fr * tg_xxz_xyyyy_1[j] + fl1_fx * (tg_xz_xyyyy_0[j] - tg_xz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyy_1[j];

                    tg_xxxz_xyyyz_0[j] = pb_x * tg_xxz_xyyyz_0[j] + fr * tg_xxz_xyyyz_1[j] + fl1_fx * (tg_xz_xyyyz_0[j] - tg_xz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyz_1[j];

                    tg_xxxz_xyyzz_0[j] = pb_x * tg_xxz_xyyzz_0[j] + fr * tg_xxz_xyyzz_1[j] + fl1_fx * (tg_xz_xyyzz_0[j] - tg_xz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyzz_1[j];

                    tg_xxxz_xyzzz_0[j] = pb_x * tg_xxz_xyzzz_0[j] + fr * tg_xxz_xyzzz_1[j] + fl1_fx * (tg_xz_xyzzz_0[j] - tg_xz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yzzz_1[j];

                    tg_xxxz_xzzzz_0[j] = pb_x * tg_xxz_xzzzz_0[j] + fr * tg_xxz_xzzzz_1[j] + fl1_fx * (tg_xz_xzzzz_0[j] - tg_xz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_zzzz_1[j];

                    tg_xxxz_yyyyy_0[j] = pb_x * tg_xxz_yyyyy_0[j] + fr * tg_xxz_yyyyy_1[j] + fl1_fx * (tg_xz_yyyyy_0[j] - tg_xz_yyyyy_1[j] * fl1_fza);

                    tg_xxxz_yyyyz_0[j] = pb_x * tg_xxz_yyyyz_0[j] + fr * tg_xxz_yyyyz_1[j] + fl1_fx * (tg_xz_yyyyz_0[j] - tg_xz_yyyyz_1[j] * fl1_fza);

                    tg_xxxz_yyyzz_0[j] = pb_x * tg_xxz_yyyzz_0[j] + fr * tg_xxz_yyyzz_1[j] + fl1_fx * (tg_xz_yyyzz_0[j] - tg_xz_yyyzz_1[j] * fl1_fza);

                    tg_xxxz_yyzzz_0[j] = pb_x * tg_xxz_yyzzz_0[j] + fr * tg_xxz_yyzzz_1[j] + fl1_fx * (tg_xz_yyzzz_0[j] - tg_xz_yyzzz_1[j] * fl1_fza);

                    tg_xxxz_yzzzz_0[j] = pb_x * tg_xxz_yzzzz_0[j] + fr * tg_xxz_yzzzz_1[j] + fl1_fx * (tg_xz_yzzzz_0[j] - tg_xz_yzzzz_1[j] * fl1_fza);

                    tg_xxxz_zzzzz_0[j] = pb_x * tg_xxz_zzzzz_0[j] + fr * tg_xxz_zzzzz_1[j] + fl1_fx * (tg_xz_zzzzz_0[j] - tg_xz_zzzzz_1[j] * fl1_fza);

                    tg_xxyy_xxxxx_0[j] = pb_x * tg_xyy_xxxxx_0[j] + fr * tg_xyy_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxx_0[j] - tg_yy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyy_xxxx_1[j];

                    tg_xxyy_xxxxy_0[j] = pb_x * tg_xyy_xxxxy_0[j] + fr * tg_xyy_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxy_0[j] - tg_yy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxy_1[j];

                    tg_xxyy_xxxxz_0[j] = pb_x * tg_xyy_xxxxz_0[j] + fr * tg_xyy_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxz_0[j] - tg_yy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxz_1[j];

                    tg_xxyy_xxxyy_0[j] = pb_x * tg_xyy_xxxyy_0[j] + fr * tg_xyy_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyy_0[j] - tg_yy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyy_1[j];

                    tg_xxyy_xxxyz_0[j] = pb_x * tg_xyy_xxxyz_0[j] + fr * tg_xyy_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyz_0[j] - tg_yy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyz_1[j];

                    tg_xxyy_xxxzz_0[j] = pb_x * tg_xyy_xxxzz_0[j] + fr * tg_xyy_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxzz_0[j] - tg_yy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxzz_1[j];

                    tg_xxyy_xxyyy_0[j] = pb_x * tg_xyy_xxyyy_0[j] + fr * tg_xyy_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyy_0[j] - tg_yy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyy_1[j];

                    tg_xxyy_xxyyz_0[j] = pb_x * tg_xyy_xxyyz_0[j] + fr * tg_xyy_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyz_0[j] - tg_yy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyz_1[j];

                    tg_xxyy_xxyzz_0[j] = pb_x * tg_xyy_xxyzz_0[j] + fr * tg_xyy_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyzz_0[j] - tg_yy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyzz_1[j];

                    tg_xxyy_xxzzz_0[j] = pb_x * tg_xyy_xxzzz_0[j] + fr * tg_xyy_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxzzz_0[j] - tg_yy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xzzz_1[j];

                    tg_xxyy_xyyyy_0[j] = pb_x * tg_xyy_xyyyy_0[j] + fr * tg_xyy_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyy_0[j] - tg_yy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyy_1[j];

                    tg_xxyy_xyyyz_0[j] = pb_x * tg_xyy_xyyyz_0[j] + fr * tg_xyy_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyz_0[j] - tg_yy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyz_1[j];

                    tg_xxyy_xyyzz_0[j] = pb_x * tg_xyy_xyyzz_0[j] + fr * tg_xyy_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyzz_0[j] - tg_yy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyzz_1[j];

                    tg_xxyy_xyzzz_0[j] = pb_x * tg_xyy_xyzzz_0[j] + fr * tg_xyy_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyzzz_0[j] - tg_yy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yzzz_1[j];

                    tg_xxyy_xzzzz_0[j] = pb_x * tg_xyy_xzzzz_0[j] + fr * tg_xyy_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xzzzz_0[j] - tg_yy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_zzzz_1[j];

                    tg_xxyy_yyyyy_0[j] = pb_x * tg_xyy_yyyyy_0[j] + fr * tg_xyy_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyy_0[j] - tg_yy_yyyyy_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSH_79_158(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (79,158)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xyy_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 79); 

                auto tg_xyy_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 80); 

                auto tg_xyy_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 81); 

                auto tg_xyy_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 82); 

                auto tg_xyy_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 83); 

                auto tg_xyz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 84); 

                auto tg_xyz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 85); 

                auto tg_xyz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 86); 

                auto tg_xyz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 87); 

                auto tg_xyz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 88); 

                auto tg_xyz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 89); 

                auto tg_xyz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 90); 

                auto tg_xyz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 91); 

                auto tg_xyz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 92); 

                auto tg_xyz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 93); 

                auto tg_xyz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 94); 

                auto tg_xyz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 95); 

                auto tg_xyz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 96); 

                auto tg_xyz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 97); 

                auto tg_xyz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 98); 

                auto tg_xyz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 99); 

                auto tg_xyz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 100); 

                auto tg_xyz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 101); 

                auto tg_xyz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 102); 

                auto tg_xyz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 103); 

                auto tg_xyz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 104); 

                auto tg_xzz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 105); 

                auto tg_xzz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 106); 

                auto tg_xzz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 107); 

                auto tg_xzz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 108); 

                auto tg_xzz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 109); 

                auto tg_xzz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 110); 

                auto tg_xzz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 111); 

                auto tg_xzz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 112); 

                auto tg_xzz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 113); 

                auto tg_xzz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 114); 

                auto tg_xzz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 115); 

                auto tg_xzz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 116); 

                auto tg_xzz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 117); 

                auto tg_xzz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 118); 

                auto tg_xzz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 119); 

                auto tg_xzz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 120); 

                auto tg_xzz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 121); 

                auto tg_xzz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 122); 

                auto tg_xzz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 123); 

                auto tg_xzz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 124); 

                auto tg_xzz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 125); 

                auto tg_yyy_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 126); 

                auto tg_yyy_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 127); 

                auto tg_yyy_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 128); 

                auto tg_yyy_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 129); 

                auto tg_yyy_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 130); 

                auto tg_yyy_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 131); 

                auto tg_yyy_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 132); 

                auto tg_yyy_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 133); 

                auto tg_yyy_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 134); 

                auto tg_yyy_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 135); 

                auto tg_yyy_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 136); 

                auto tg_yyy_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 137); 

                auto tg_yyy_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 138); 

                auto tg_yyy_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 139); 

                auto tg_yyy_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 140); 

                auto tg_yyy_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 141); 

                auto tg_yyy_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 142); 

                auto tg_yyy_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 143); 

                auto tg_yyy_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 144); 

                auto tg_yyy_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 145); 

                auto tg_yyy_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 146); 

                auto tg_yyz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 147); 

                auto tg_yyz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 148); 

                auto tg_yyz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 149); 

                auto tg_yyz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 150); 

                auto tg_yyz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 151); 

                auto tg_yyz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 152); 

                auto tg_yyz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 153); 

                auto tg_yyz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 154); 

                auto tg_yyz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 155); 

                auto tg_yyz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 156); 

                auto tg_yyz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 157); 

                auto tg_xyy_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 79); 

                auto tg_xyy_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 80); 

                auto tg_xyy_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 81); 

                auto tg_xyy_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 82); 

                auto tg_xyy_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 83); 

                auto tg_xyz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 84); 

                auto tg_xyz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 85); 

                auto tg_xyz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 86); 

                auto tg_xyz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 87); 

                auto tg_xyz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 88); 

                auto tg_xyz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 89); 

                auto tg_xyz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 90); 

                auto tg_xyz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 91); 

                auto tg_xyz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 92); 

                auto tg_xyz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 93); 

                auto tg_xyz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 94); 

                auto tg_xyz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 95); 

                auto tg_xyz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 96); 

                auto tg_xyz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 97); 

                auto tg_xyz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 98); 

                auto tg_xyz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 99); 

                auto tg_xyz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 100); 

                auto tg_xyz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 101); 

                auto tg_xyz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 102); 

                auto tg_xyz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 103); 

                auto tg_xyz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 104); 

                auto tg_xzz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 105); 

                auto tg_xzz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 106); 

                auto tg_xzz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 107); 

                auto tg_xzz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 108); 

                auto tg_xzz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 109); 

                auto tg_xzz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 110); 

                auto tg_xzz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 111); 

                auto tg_xzz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 112); 

                auto tg_xzz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 113); 

                auto tg_xzz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 114); 

                auto tg_xzz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 115); 

                auto tg_xzz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 116); 

                auto tg_xzz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 117); 

                auto tg_xzz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 118); 

                auto tg_xzz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 119); 

                auto tg_xzz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 120); 

                auto tg_xzz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 121); 

                auto tg_xzz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 122); 

                auto tg_xzz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 123); 

                auto tg_xzz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 124); 

                auto tg_xzz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 125); 

                auto tg_yyy_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 126); 

                auto tg_yyy_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 127); 

                auto tg_yyy_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 128); 

                auto tg_yyy_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 129); 

                auto tg_yyy_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 130); 

                auto tg_yyy_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 131); 

                auto tg_yyy_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 132); 

                auto tg_yyy_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 133); 

                auto tg_yyy_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 134); 

                auto tg_yyy_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 135); 

                auto tg_yyy_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 136); 

                auto tg_yyy_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 137); 

                auto tg_yyy_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 138); 

                auto tg_yyy_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 139); 

                auto tg_yyy_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 140); 

                auto tg_yyy_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 141); 

                auto tg_yyy_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 142); 

                auto tg_yyy_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 143); 

                auto tg_yyy_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 144); 

                auto tg_yyy_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 145); 

                auto tg_yyy_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 146); 

                auto tg_yyz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 147); 

                auto tg_yyz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 148); 

                auto tg_yyz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 149); 

                auto tg_yyz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 150); 

                auto tg_yyz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 151); 

                auto tg_yyz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 152); 

                auto tg_yyz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 153); 

                auto tg_yyz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 154); 

                auto tg_yyz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 155); 

                auto tg_yyz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 156); 

                auto tg_yyz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 157); 

                auto tg_yy_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 79); 

                auto tg_yy_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 80); 

                auto tg_yy_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 81); 

                auto tg_yy_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 82); 

                auto tg_yy_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 83); 

                auto tg_yz_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 84); 

                auto tg_yz_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 85); 

                auto tg_yz_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 86); 

                auto tg_yz_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 87); 

                auto tg_yz_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 88); 

                auto tg_yz_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 89); 

                auto tg_yz_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 90); 

                auto tg_yz_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 91); 

                auto tg_yz_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 92); 

                auto tg_yz_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 93); 

                auto tg_yz_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 94); 

                auto tg_yz_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 95); 

                auto tg_yz_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 96); 

                auto tg_yz_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 97); 

                auto tg_yz_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 98); 

                auto tg_yz_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 99); 

                auto tg_yz_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 100); 

                auto tg_yz_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 101); 

                auto tg_yz_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 102); 

                auto tg_yz_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 103); 

                auto tg_yz_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 104); 

                auto tg_zz_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 105); 

                auto tg_zz_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 106); 

                auto tg_zz_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 107); 

                auto tg_zz_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 108); 

                auto tg_zz_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 109); 

                auto tg_zz_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 110); 

                auto tg_zz_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 111); 

                auto tg_zz_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 112); 

                auto tg_zz_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 113); 

                auto tg_zz_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 114); 

                auto tg_zz_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 115); 

                auto tg_zz_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 116); 

                auto tg_zz_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 117); 

                auto tg_zz_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 118); 

                auto tg_zz_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 119); 

                auto tg_zz_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 120); 

                auto tg_zz_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 121); 

                auto tg_zz_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 122); 

                auto tg_zz_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 123); 

                auto tg_zz_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 124); 

                auto tg_zz_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 125); 

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

                auto tg_xyz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 60); 

                auto tg_xyz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 61); 

                auto tg_xyz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 62); 

                auto tg_xyz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 63); 

                auto tg_xyz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 64); 

                auto tg_xyz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 65); 

                auto tg_xyz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 66); 

                auto tg_xyz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 67); 

                auto tg_xyz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 68); 

                auto tg_xyz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 69); 

                auto tg_xyz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 70); 

                auto tg_xyz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 71); 

                auto tg_xyz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 72); 

                auto tg_xyz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 73); 

                auto tg_xyz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 74); 

                auto tg_xzz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 75); 

                auto tg_xzz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 76); 

                auto tg_xzz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 77); 

                auto tg_xzz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 78); 

                auto tg_xzz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 79); 

                auto tg_xzz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 80); 

                auto tg_xzz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 81); 

                auto tg_xzz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 82); 

                auto tg_xzz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 83); 

                auto tg_xzz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 84); 

                auto tg_xzz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 85); 

                auto tg_xzz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 86); 

                auto tg_xzz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 87); 

                auto tg_xzz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 88); 

                auto tg_xzz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 89); 

                auto tg_yyy_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 90); 

                auto tg_yyy_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 91); 

                auto tg_yyy_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 92); 

                auto tg_yyy_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 93); 

                auto tg_yyy_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 94); 

                auto tg_yyy_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 95); 

                auto tg_yyy_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 96); 

                auto tg_yyy_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 97); 

                auto tg_yyy_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 98); 

                auto tg_yyy_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 99); 

                auto tg_yyy_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 100); 

                auto tg_yyy_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 101); 

                auto tg_yyy_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 102); 

                auto tg_yyy_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 103); 

                auto tg_yyy_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 104); 

                auto tg_yyz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 105); 

                auto tg_yyz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 106); 

                auto tg_yyz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 107); 

                auto tg_yyz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 108); 

                auto tg_yyz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 109); 

                auto tg_yyz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 110); 

                auto tg_yyz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 111); 

                auto tg_yyz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 112); 

                auto tg_yyz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 113); 

                auto tg_yyz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 114); 

                auto tg_yyz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 115); 

                // set up pointers to integrals

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

                // Batch of Integrals (79,158)

                #pragma omp simd aligned(fxn, fza, tg_xxyy_yyyyz_0, tg_xxyy_yyyzz_0, tg_xxyy_yyzzz_0, \
                                         tg_xxyy_yzzzz_0, tg_xxyy_zzzzz_0, tg_xxyz_xxxxx_0, tg_xxyz_xxxxy_0, tg_xxyz_xxxxz_0, \
                                         tg_xxyz_xxxyy_0, tg_xxyz_xxxyz_0, tg_xxyz_xxxzz_0, tg_xxyz_xxyyy_0, tg_xxyz_xxyyz_0, \
                                         tg_xxyz_xxyzz_0, tg_xxyz_xxzzz_0, tg_xxyz_xyyyy_0, tg_xxyz_xyyyz_0, tg_xxyz_xyyzz_0, \
                                         tg_xxyz_xyzzz_0, tg_xxyz_xzzzz_0, tg_xxyz_yyyyy_0, tg_xxyz_yyyyz_0, tg_xxyz_yyyzz_0, \
                                         tg_xxyz_yyzzz_0, tg_xxyz_yzzzz_0, tg_xxyz_zzzzz_0, tg_xxzz_xxxxx_0, tg_xxzz_xxxxy_0, \
                                         tg_xxzz_xxxxz_0, tg_xxzz_xxxyy_0, tg_xxzz_xxxyz_0, tg_xxzz_xxxzz_0, tg_xxzz_xxyyy_0, \
                                         tg_xxzz_xxyyz_0, tg_xxzz_xxyzz_0, tg_xxzz_xxzzz_0, tg_xxzz_xyyyy_0, tg_xxzz_xyyyz_0, \
                                         tg_xxzz_xyyzz_0, tg_xxzz_xyzzz_0, tg_xxzz_xzzzz_0, tg_xxzz_yyyyy_0, tg_xxzz_yyyyz_0, \
                                         tg_xxzz_yyyzz_0, tg_xxzz_yyzzz_0, tg_xxzz_yzzzz_0, tg_xxzz_zzzzz_0, tg_xyy_yyyyz_0, \
                                         tg_xyy_yyyyz_1, tg_xyy_yyyzz_0, tg_xyy_yyyzz_1, tg_xyy_yyzzz_0, tg_xyy_yyzzz_1, \
                                         tg_xyy_yzzzz_0, tg_xyy_yzzzz_1, tg_xyy_zzzzz_0, tg_xyy_zzzzz_1, tg_xyyy_xxxxx_0, \
                                         tg_xyyy_xxxxy_0, tg_xyyy_xxxxz_0, tg_xyyy_xxxyy_0, tg_xyyy_xxxyz_0, tg_xyyy_xxxzz_0, \
                                         tg_xyyy_xxyyy_0, tg_xyyy_xxyyz_0, tg_xyyy_xxyzz_0, tg_xyyy_xxzzz_0, tg_xyyy_xyyyy_0, \
                                         tg_xyyy_xyyyz_0, tg_xyyy_xyyzz_0, tg_xyyy_xyzzz_0, tg_xyyy_xzzzz_0, tg_xyyy_yyyyy_0, \
                                         tg_xyyy_yyyyz_0, tg_xyyy_yyyzz_0, tg_xyyy_yyzzz_0, tg_xyyy_yzzzz_0, tg_xyyy_zzzzz_0, \
                                         tg_xyyz_xxxxx_0, tg_xyyz_xxxxy_0, tg_xyyz_xxxxz_0, tg_xyyz_xxxyy_0, tg_xyyz_xxxyz_0, \
                                         tg_xyyz_xxxzz_0, tg_xyyz_xxyyy_0, tg_xyyz_xxyyz_0, tg_xyyz_xxyzz_0, tg_xyyz_xxzzz_0, \
                                         tg_xyyz_xyyyy_0, tg_xyz_xxxx_1, tg_xyz_xxxxx_0, tg_xyz_xxxxx_1, tg_xyz_xxxxy_0, \
                                         tg_xyz_xxxxy_1, tg_xyz_xxxxz_0, tg_xyz_xxxxz_1, tg_xyz_xxxy_1, tg_xyz_xxxyy_0, \
                                         tg_xyz_xxxyy_1, tg_xyz_xxxyz_0, tg_xyz_xxxyz_1, tg_xyz_xxxz_1, tg_xyz_xxxzz_0, \
                                         tg_xyz_xxxzz_1, tg_xyz_xxyy_1, tg_xyz_xxyyy_0, tg_xyz_xxyyy_1, tg_xyz_xxyyz_0, \
                                         tg_xyz_xxyyz_1, tg_xyz_xxyz_1, tg_xyz_xxyzz_0, tg_xyz_xxyzz_1, tg_xyz_xxzz_1, \
                                         tg_xyz_xxzzz_0, tg_xyz_xxzzz_1, tg_xyz_xyyy_1, tg_xyz_xyyyy_0, tg_xyz_xyyyy_1, \
                                         tg_xyz_xyyyz_0, tg_xyz_xyyyz_1, tg_xyz_xyyz_1, tg_xyz_xyyzz_0, tg_xyz_xyyzz_1, \
                                         tg_xyz_xyzz_1, tg_xyz_xyzzz_0, tg_xyz_xyzzz_1, tg_xyz_xzzz_1, tg_xyz_xzzzz_0, \
                                         tg_xyz_xzzzz_1, tg_xyz_yyyy_1, tg_xyz_yyyyy_0, tg_xyz_yyyyy_1, tg_xyz_yyyyz_0, \
                                         tg_xyz_yyyyz_1, tg_xyz_yyyz_1, tg_xyz_yyyzz_0, tg_xyz_yyyzz_1, tg_xyz_yyzz_1, \
                                         tg_xyz_yyzzz_0, tg_xyz_yyzzz_1, tg_xyz_yzzz_1, tg_xyz_yzzzz_0, tg_xyz_yzzzz_1, \
                                         tg_xyz_zzzz_1, tg_xyz_zzzzz_0, tg_xyz_zzzzz_1, tg_xzz_xxxx_1, tg_xzz_xxxxx_0, \
                                         tg_xzz_xxxxx_1, tg_xzz_xxxxy_0, tg_xzz_xxxxy_1, tg_xzz_xxxxz_0, tg_xzz_xxxxz_1, \
                                         tg_xzz_xxxy_1, tg_xzz_xxxyy_0, tg_xzz_xxxyy_1, tg_xzz_xxxyz_0, tg_xzz_xxxyz_1, \
                                         tg_xzz_xxxz_1, tg_xzz_xxxzz_0, tg_xzz_xxxzz_1, tg_xzz_xxyy_1, tg_xzz_xxyyy_0, \
                                         tg_xzz_xxyyy_1, tg_xzz_xxyyz_0, tg_xzz_xxyyz_1, tg_xzz_xxyz_1, tg_xzz_xxyzz_0, \
                                         tg_xzz_xxyzz_1, tg_xzz_xxzz_1, tg_xzz_xxzzz_0, tg_xzz_xxzzz_1, tg_xzz_xyyy_1, \
                                         tg_xzz_xyyyy_0, tg_xzz_xyyyy_1, tg_xzz_xyyyz_0, tg_xzz_xyyyz_1, tg_xzz_xyyz_1, \
                                         tg_xzz_xyyzz_0, tg_xzz_xyyzz_1, tg_xzz_xyzz_1, tg_xzz_xyzzz_0, tg_xzz_xyzzz_1, \
                                         tg_xzz_xzzz_1, tg_xzz_xzzzz_0, tg_xzz_xzzzz_1, tg_xzz_yyyy_1, tg_xzz_yyyyy_0, \
                                         tg_xzz_yyyyy_1, tg_xzz_yyyyz_0, tg_xzz_yyyyz_1, tg_xzz_yyyz_1, tg_xzz_yyyzz_0, \
                                         tg_xzz_yyyzz_1, tg_xzz_yyzz_1, tg_xzz_yyzzz_0, tg_xzz_yyzzz_1, tg_xzz_yzzz_1, \
                                         tg_xzz_yzzzz_0, tg_xzz_yzzzz_1, tg_xzz_zzzz_1, tg_xzz_zzzzz_0, tg_xzz_zzzzz_1, \
                                         tg_yy_yyyyz_0, tg_yy_yyyyz_1, tg_yy_yyyzz_0, tg_yy_yyyzz_1, tg_yy_yyzzz_0, \
                                         tg_yy_yyzzz_1, tg_yy_yzzzz_0, tg_yy_yzzzz_1, tg_yy_zzzzz_0, tg_yy_zzzzz_1, \
                                         tg_yyy_xxxx_1, tg_yyy_xxxxx_0, tg_yyy_xxxxx_1, tg_yyy_xxxxy_0, tg_yyy_xxxxy_1, \
                                         tg_yyy_xxxxz_0, tg_yyy_xxxxz_1, tg_yyy_xxxy_1, tg_yyy_xxxyy_0, tg_yyy_xxxyy_1, \
                                         tg_yyy_xxxyz_0, tg_yyy_xxxyz_1, tg_yyy_xxxz_1, tg_yyy_xxxzz_0, tg_yyy_xxxzz_1, \
                                         tg_yyy_xxyy_1, tg_yyy_xxyyy_0, tg_yyy_xxyyy_1, tg_yyy_xxyyz_0, tg_yyy_xxyyz_1, \
                                         tg_yyy_xxyz_1, tg_yyy_xxyzz_0, tg_yyy_xxyzz_1, tg_yyy_xxzz_1, tg_yyy_xxzzz_0, \
                                         tg_yyy_xxzzz_1, tg_yyy_xyyy_1, tg_yyy_xyyyy_0, tg_yyy_xyyyy_1, tg_yyy_xyyyz_0, \
                                         tg_yyy_xyyyz_1, tg_yyy_xyyz_1, tg_yyy_xyyzz_0, tg_yyy_xyyzz_1, tg_yyy_xyzz_1, \
                                         tg_yyy_xyzzz_0, tg_yyy_xyzzz_1, tg_yyy_xzzz_1, tg_yyy_xzzzz_0, tg_yyy_xzzzz_1, \
                                         tg_yyy_yyyy_1, tg_yyy_yyyyy_0, tg_yyy_yyyyy_1, tg_yyy_yyyyz_0, tg_yyy_yyyyz_1, \
                                         tg_yyy_yyyz_1, tg_yyy_yyyzz_0, tg_yyy_yyyzz_1, tg_yyy_yyzz_1, tg_yyy_yyzzz_0, \
                                         tg_yyy_yyzzz_1, tg_yyy_yzzz_1, tg_yyy_yzzzz_0, tg_yyy_yzzzz_1, tg_yyy_zzzz_1, \
                                         tg_yyy_zzzzz_0, tg_yyy_zzzzz_1, tg_yyz_xxxx_1, tg_yyz_xxxxx_0, tg_yyz_xxxxx_1, \
                                         tg_yyz_xxxxy_0, tg_yyz_xxxxy_1, tg_yyz_xxxxz_0, tg_yyz_xxxxz_1, tg_yyz_xxxy_1, \
                                         tg_yyz_xxxyy_0, tg_yyz_xxxyy_1, tg_yyz_xxxyz_0, tg_yyz_xxxyz_1, tg_yyz_xxxz_1, \
                                         tg_yyz_xxxzz_0, tg_yyz_xxxzz_1, tg_yyz_xxyy_1, tg_yyz_xxyyy_0, tg_yyz_xxyyy_1, \
                                         tg_yyz_xxyyz_0, tg_yyz_xxyyz_1, tg_yyz_xxyz_1, tg_yyz_xxyzz_0, tg_yyz_xxyzz_1, \
                                         tg_yyz_xxzz_1, tg_yyz_xxzzz_0, tg_yyz_xxzzz_1, tg_yyz_xyyy_1, tg_yyz_xyyyy_0, \
                                         tg_yyz_xyyyy_1, tg_yyz_xyyz_1, tg_yyz_xyzz_1, tg_yyz_xzzz_1, tg_yyz_yyyy_1, \
                                         tg_yz_xxxxx_0, tg_yz_xxxxx_1, tg_yz_xxxxy_0, tg_yz_xxxxy_1, tg_yz_xxxxz_0, \
                                         tg_yz_xxxxz_1, tg_yz_xxxyy_0, tg_yz_xxxyy_1, tg_yz_xxxyz_0, tg_yz_xxxyz_1, \
                                         tg_yz_xxxzz_0, tg_yz_xxxzz_1, tg_yz_xxyyy_0, tg_yz_xxyyy_1, tg_yz_xxyyz_0, \
                                         tg_yz_xxyyz_1, tg_yz_xxyzz_0, tg_yz_xxyzz_1, tg_yz_xxzzz_0, tg_yz_xxzzz_1, \
                                         tg_yz_xyyyy_0, tg_yz_xyyyy_1, tg_yz_xyyyz_0, tg_yz_xyyyz_1, tg_yz_xyyzz_0, \
                                         tg_yz_xyyzz_1, tg_yz_xyzzz_0, tg_yz_xyzzz_1, tg_yz_xzzzz_0, tg_yz_xzzzz_1, \
                                         tg_yz_yyyyy_0, tg_yz_yyyyy_1, tg_yz_yyyyz_0, tg_yz_yyyyz_1, tg_yz_yyyzz_0, \
                                         tg_yz_yyyzz_1, tg_yz_yyzzz_0, tg_yz_yyzzz_1, tg_yz_yzzzz_0, tg_yz_yzzzz_1, \
                                         tg_yz_zzzzz_0, tg_yz_zzzzz_1, tg_zz_xxxxx_0, tg_zz_xxxxx_1, tg_zz_xxxxy_0, \
                                         tg_zz_xxxxy_1, tg_zz_xxxxz_0, tg_zz_xxxxz_1, tg_zz_xxxyy_0, tg_zz_xxxyy_1, \
                                         tg_zz_xxxyz_0, tg_zz_xxxyz_1, tg_zz_xxxzz_0, tg_zz_xxxzz_1, tg_zz_xxyyy_0, \
                                         tg_zz_xxyyy_1, tg_zz_xxyyz_0, tg_zz_xxyyz_1, tg_zz_xxyzz_0, tg_zz_xxyzz_1, \
                                         tg_zz_xxzzz_0, tg_zz_xxzzz_1, tg_zz_xyyyy_0, tg_zz_xyyyy_1, tg_zz_xyyyz_0, \
                                         tg_zz_xyyyz_1, tg_zz_xyyzz_0, tg_zz_xyyzz_1, tg_zz_xyzzz_0, tg_zz_xyzzz_1, \
                                         tg_zz_xzzzz_0, tg_zz_xzzzz_1, tg_zz_yyyyy_0, tg_zz_yyyyy_1, tg_zz_yyyyz_0, \
                                         tg_zz_yyyyz_1, tg_zz_yyyzz_0, tg_zz_yyyzz_1, tg_zz_yyzzz_0, tg_zz_yyzzz_1, \
                                         tg_zz_yzzzz_0, tg_zz_yzzzz_1, tg_zz_zzzzz_0, tg_zz_zzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyy_yyyyz_0[j] = pb_x * tg_xyy_yyyyz_0[j] + fr * tg_xyy_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyz_0[j] - tg_yy_yyyyz_1[j] * fl1_fza);

                    tg_xxyy_yyyzz_0[j] = pb_x * tg_xyy_yyyzz_0[j] + fr * tg_xyy_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyzz_0[j] - tg_yy_yyyzz_1[j] * fl1_fza);

                    tg_xxyy_yyzzz_0[j] = pb_x * tg_xyy_yyzzz_0[j] + fr * tg_xyy_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyzzz_0[j] - tg_yy_yyzzz_1[j] * fl1_fza);

                    tg_xxyy_yzzzz_0[j] = pb_x * tg_xyy_yzzzz_0[j] + fr * tg_xyy_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yzzzz_0[j] - tg_yy_yzzzz_1[j] * fl1_fza);

                    tg_xxyy_zzzzz_0[j] = pb_x * tg_xyy_zzzzz_0[j] + fr * tg_xyy_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_zzzzz_0[j] - tg_yy_zzzzz_1[j] * fl1_fza);

                    tg_xxyz_xxxxx_0[j] = pb_x * tg_xyz_xxxxx_0[j] + fr * tg_xyz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxx_0[j] - tg_yz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyz_xxxx_1[j];

                    tg_xxyz_xxxxy_0[j] = pb_x * tg_xyz_xxxxy_0[j] + fr * tg_xyz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxy_0[j] - tg_yz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxy_1[j];

                    tg_xxyz_xxxxz_0[j] = pb_x * tg_xyz_xxxxz_0[j] + fr * tg_xyz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxz_0[j] - tg_yz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxz_1[j];

                    tg_xxyz_xxxyy_0[j] = pb_x * tg_xyz_xxxyy_0[j] + fr * tg_xyz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyy_0[j] - tg_yz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyy_1[j];

                    tg_xxyz_xxxyz_0[j] = pb_x * tg_xyz_xxxyz_0[j] + fr * tg_xyz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyz_0[j] - tg_yz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyz_1[j];

                    tg_xxyz_xxxzz_0[j] = pb_x * tg_xyz_xxxzz_0[j] + fr * tg_xyz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxzz_0[j] - tg_yz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxzz_1[j];

                    tg_xxyz_xxyyy_0[j] = pb_x * tg_xyz_xxyyy_0[j] + fr * tg_xyz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyy_0[j] - tg_yz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyy_1[j];

                    tg_xxyz_xxyyz_0[j] = pb_x * tg_xyz_xxyyz_0[j] + fr * tg_xyz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyz_0[j] - tg_yz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyz_1[j];

                    tg_xxyz_xxyzz_0[j] = pb_x * tg_xyz_xxyzz_0[j] + fr * tg_xyz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyzz_0[j] - tg_yz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyzz_1[j];

                    tg_xxyz_xxzzz_0[j] = pb_x * tg_xyz_xxzzz_0[j] + fr * tg_xyz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxzzz_0[j] - tg_yz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xzzz_1[j];

                    tg_xxyz_xyyyy_0[j] = pb_x * tg_xyz_xyyyy_0[j] + fr * tg_xyz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyy_0[j] - tg_yz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyy_1[j];

                    tg_xxyz_xyyyz_0[j] = pb_x * tg_xyz_xyyyz_0[j] + fr * tg_xyz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyz_0[j] - tg_yz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyz_1[j];

                    tg_xxyz_xyyzz_0[j] = pb_x * tg_xyz_xyyzz_0[j] + fr * tg_xyz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyzz_0[j] - tg_yz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyzz_1[j];

                    tg_xxyz_xyzzz_0[j] = pb_x * tg_xyz_xyzzz_0[j] + fr * tg_xyz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyzzz_0[j] - tg_yz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yzzz_1[j];

                    tg_xxyz_xzzzz_0[j] = pb_x * tg_xyz_xzzzz_0[j] + fr * tg_xyz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xzzzz_0[j] - tg_yz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_zzzz_1[j];

                    tg_xxyz_yyyyy_0[j] = pb_x * tg_xyz_yyyyy_0[j] + fr * tg_xyz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyy_0[j] - tg_yz_yyyyy_1[j] * fl1_fza);

                    tg_xxyz_yyyyz_0[j] = pb_x * tg_xyz_yyyyz_0[j] + fr * tg_xyz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyz_0[j] - tg_yz_yyyyz_1[j] * fl1_fza);

                    tg_xxyz_yyyzz_0[j] = pb_x * tg_xyz_yyyzz_0[j] + fr * tg_xyz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyzz_0[j] - tg_yz_yyyzz_1[j] * fl1_fza);

                    tg_xxyz_yyzzz_0[j] = pb_x * tg_xyz_yyzzz_0[j] + fr * tg_xyz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyzzz_0[j] - tg_yz_yyzzz_1[j] * fl1_fza);

                    tg_xxyz_yzzzz_0[j] = pb_x * tg_xyz_yzzzz_0[j] + fr * tg_xyz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yzzzz_0[j] - tg_yz_yzzzz_1[j] * fl1_fza);

                    tg_xxyz_zzzzz_0[j] = pb_x * tg_xyz_zzzzz_0[j] + fr * tg_xyz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_zzzzz_0[j] - tg_yz_zzzzz_1[j] * fl1_fza);

                    tg_xxzz_xxxxx_0[j] = pb_x * tg_xzz_xxxxx_0[j] + fr * tg_xzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxx_0[j] - tg_zz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzz_xxxx_1[j];

                    tg_xxzz_xxxxy_0[j] = pb_x * tg_xzz_xxxxy_0[j] + fr * tg_xzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxy_0[j] - tg_zz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxy_1[j];

                    tg_xxzz_xxxxz_0[j] = pb_x * tg_xzz_xxxxz_0[j] + fr * tg_xzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxz_0[j] - tg_zz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxz_1[j];

                    tg_xxzz_xxxyy_0[j] = pb_x * tg_xzz_xxxyy_0[j] + fr * tg_xzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyy_0[j] - tg_zz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyy_1[j];

                    tg_xxzz_xxxyz_0[j] = pb_x * tg_xzz_xxxyz_0[j] + fr * tg_xzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyz_0[j] - tg_zz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyz_1[j];

                    tg_xxzz_xxxzz_0[j] = pb_x * tg_xzz_xxxzz_0[j] + fr * tg_xzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxzz_0[j] - tg_zz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxzz_1[j];

                    tg_xxzz_xxyyy_0[j] = pb_x * tg_xzz_xxyyy_0[j] + fr * tg_xzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyy_0[j] - tg_zz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyy_1[j];

                    tg_xxzz_xxyyz_0[j] = pb_x * tg_xzz_xxyyz_0[j] + fr * tg_xzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyz_0[j] - tg_zz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyz_1[j];

                    tg_xxzz_xxyzz_0[j] = pb_x * tg_xzz_xxyzz_0[j] + fr * tg_xzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyzz_0[j] - tg_zz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyzz_1[j];

                    tg_xxzz_xxzzz_0[j] = pb_x * tg_xzz_xxzzz_0[j] + fr * tg_xzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxzzz_0[j] - tg_zz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xzzz_1[j];

                    tg_xxzz_xyyyy_0[j] = pb_x * tg_xzz_xyyyy_0[j] + fr * tg_xzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyy_0[j] - tg_zz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyy_1[j];

                    tg_xxzz_xyyyz_0[j] = pb_x * tg_xzz_xyyyz_0[j] + fr * tg_xzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyz_0[j] - tg_zz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyz_1[j];

                    tg_xxzz_xyyzz_0[j] = pb_x * tg_xzz_xyyzz_0[j] + fr * tg_xzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyzz_0[j] - tg_zz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyzz_1[j];

                    tg_xxzz_xyzzz_0[j] = pb_x * tg_xzz_xyzzz_0[j] + fr * tg_xzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyzzz_0[j] - tg_zz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yzzz_1[j];

                    tg_xxzz_xzzzz_0[j] = pb_x * tg_xzz_xzzzz_0[j] + fr * tg_xzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzzzz_0[j] - tg_zz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_zzzz_1[j];

                    tg_xxzz_yyyyy_0[j] = pb_x * tg_xzz_yyyyy_0[j] + fr * tg_xzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyy_0[j] - tg_zz_yyyyy_1[j] * fl1_fza);

                    tg_xxzz_yyyyz_0[j] = pb_x * tg_xzz_yyyyz_0[j] + fr * tg_xzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyz_0[j] - tg_zz_yyyyz_1[j] * fl1_fza);

                    tg_xxzz_yyyzz_0[j] = pb_x * tg_xzz_yyyzz_0[j] + fr * tg_xzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyzz_0[j] - tg_zz_yyyzz_1[j] * fl1_fza);

                    tg_xxzz_yyzzz_0[j] = pb_x * tg_xzz_yyzzz_0[j] + fr * tg_xzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyzzz_0[j] - tg_zz_yyzzz_1[j] * fl1_fza);

                    tg_xxzz_yzzzz_0[j] = pb_x * tg_xzz_yzzzz_0[j] + fr * tg_xzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzzzz_0[j] - tg_zz_yzzzz_1[j] * fl1_fza);

                    tg_xxzz_zzzzz_0[j] = pb_x * tg_xzz_zzzzz_0[j] + fr * tg_xzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzzzz_0[j] - tg_zz_zzzzz_1[j] * fl1_fza);

                    tg_xyyy_xxxxx_0[j] = pb_x * tg_yyy_xxxxx_0[j] + fr * tg_yyy_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxx_1[j];

                    tg_xyyy_xxxxy_0[j] = pb_x * tg_yyy_xxxxy_0[j] + fr * tg_yyy_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxy_1[j];

                    tg_xyyy_xxxxz_0[j] = pb_x * tg_yyy_xxxxz_0[j] + fr * tg_yyy_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxz_1[j];

                    tg_xyyy_xxxyy_0[j] = pb_x * tg_yyy_xxxyy_0[j] + fr * tg_yyy_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyy_1[j];

                    tg_xyyy_xxxyz_0[j] = pb_x * tg_yyy_xxxyz_0[j] + fr * tg_yyy_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyz_1[j];

                    tg_xyyy_xxxzz_0[j] = pb_x * tg_yyy_xxxzz_0[j] + fr * tg_yyy_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxzz_1[j];

                    tg_xyyy_xxyyy_0[j] = pb_x * tg_yyy_xxyyy_0[j] + fr * tg_yyy_xxyyy_1[j] + fl1_fxn * tg_yyy_xyyy_1[j];

                    tg_xyyy_xxyyz_0[j] = pb_x * tg_yyy_xxyyz_0[j] + fr * tg_yyy_xxyyz_1[j] + fl1_fxn * tg_yyy_xyyz_1[j];

                    tg_xyyy_xxyzz_0[j] = pb_x * tg_yyy_xxyzz_0[j] + fr * tg_yyy_xxyzz_1[j] + fl1_fxn * tg_yyy_xyzz_1[j];

                    tg_xyyy_xxzzz_0[j] = pb_x * tg_yyy_xxzzz_0[j] + fr * tg_yyy_xxzzz_1[j] + fl1_fxn * tg_yyy_xzzz_1[j];

                    tg_xyyy_xyyyy_0[j] = pb_x * tg_yyy_xyyyy_0[j] + fr * tg_yyy_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyy_1[j];

                    tg_xyyy_xyyyz_0[j] = pb_x * tg_yyy_xyyyz_0[j] + fr * tg_yyy_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyz_1[j];

                    tg_xyyy_xyyzz_0[j] = pb_x * tg_yyy_xyyzz_0[j] + fr * tg_yyy_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyzz_1[j];

                    tg_xyyy_xyzzz_0[j] = pb_x * tg_yyy_xyzzz_0[j] + fr * tg_yyy_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yzzz_1[j];

                    tg_xyyy_xzzzz_0[j] = pb_x * tg_yyy_xzzzz_0[j] + fr * tg_yyy_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zzzz_1[j];

                    tg_xyyy_yyyyy_0[j] = pb_x * tg_yyy_yyyyy_0[j] + fr * tg_yyy_yyyyy_1[j];

                    tg_xyyy_yyyyz_0[j] = pb_x * tg_yyy_yyyyz_0[j] + fr * tg_yyy_yyyyz_1[j];

                    tg_xyyy_yyyzz_0[j] = pb_x * tg_yyy_yyyzz_0[j] + fr * tg_yyy_yyyzz_1[j];

                    tg_xyyy_yyzzz_0[j] = pb_x * tg_yyy_yyzzz_0[j] + fr * tg_yyy_yyzzz_1[j];

                    tg_xyyy_yzzzz_0[j] = pb_x * tg_yyy_yzzzz_0[j] + fr * tg_yyy_yzzzz_1[j];

                    tg_xyyy_zzzzz_0[j] = pb_x * tg_yyy_zzzzz_0[j] + fr * tg_yyy_zzzzz_1[j];

                    tg_xyyz_xxxxx_0[j] = pb_x * tg_yyz_xxxxx_0[j] + fr * tg_yyz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxx_1[j];

                    tg_xyyz_xxxxy_0[j] = pb_x * tg_yyz_xxxxy_0[j] + fr * tg_yyz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxy_1[j];

                    tg_xyyz_xxxxz_0[j] = pb_x * tg_yyz_xxxxz_0[j] + fr * tg_yyz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxz_1[j];

                    tg_xyyz_xxxyy_0[j] = pb_x * tg_yyz_xxxyy_0[j] + fr * tg_yyz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyy_1[j];

                    tg_xyyz_xxxyz_0[j] = pb_x * tg_yyz_xxxyz_0[j] + fr * tg_yyz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyz_1[j];

                    tg_xyyz_xxxzz_0[j] = pb_x * tg_yyz_xxxzz_0[j] + fr * tg_yyz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxzz_1[j];

                    tg_xyyz_xxyyy_0[j] = pb_x * tg_yyz_xxyyy_0[j] + fr * tg_yyz_xxyyy_1[j] + fl1_fxn * tg_yyz_xyyy_1[j];

                    tg_xyyz_xxyyz_0[j] = pb_x * tg_yyz_xxyyz_0[j] + fr * tg_yyz_xxyyz_1[j] + fl1_fxn * tg_yyz_xyyz_1[j];

                    tg_xyyz_xxyzz_0[j] = pb_x * tg_yyz_xxyzz_0[j] + fr * tg_yyz_xxyzz_1[j] + fl1_fxn * tg_yyz_xyzz_1[j];

                    tg_xyyz_xxzzz_0[j] = pb_x * tg_yyz_xxzzz_0[j] + fr * tg_yyz_xxzzz_1[j] + fl1_fxn * tg_yyz_xzzz_1[j];

                    tg_xyyz_xyyyy_0[j] = pb_x * tg_yyz_xyyyy_0[j] + fr * tg_yyz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSH_158_237(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (158,237)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_yyy_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 126); 

                auto tg_yyy_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 127); 

                auto tg_yyy_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 128); 

                auto tg_yyy_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 129); 

                auto tg_yyy_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 130); 

                auto tg_yyy_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 131); 

                auto tg_yyy_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 132); 

                auto tg_yyy_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 133); 

                auto tg_yyy_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 134); 

                auto tg_yyy_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 135); 

                auto tg_yyy_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 136); 

                auto tg_yyy_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 137); 

                auto tg_yyy_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 138); 

                auto tg_yyy_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 139); 

                auto tg_yyy_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 140); 

                auto tg_yyy_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 141); 

                auto tg_yyy_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 142); 

                auto tg_yyy_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 143); 

                auto tg_yyy_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 144); 

                auto tg_yyy_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 145); 

                auto tg_yyy_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 146); 

                auto tg_yyz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 147); 

                auto tg_yyz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 148); 

                auto tg_yyz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 149); 

                auto tg_yyz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 150); 

                auto tg_yyz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 151); 

                auto tg_yyz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 152); 

                auto tg_yyz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 158); 

                auto tg_yyz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 159); 

                auto tg_yyz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 160); 

                auto tg_yyz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 161); 

                auto tg_yyz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 162); 

                auto tg_yyz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 163); 

                auto tg_yyz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 164); 

                auto tg_yyz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 165); 

                auto tg_yyz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 166); 

                auto tg_yyz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 167); 

                auto tg_yzz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 168); 

                auto tg_yzz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 169); 

                auto tg_yzz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 170); 

                auto tg_yzz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 171); 

                auto tg_yzz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 172); 

                auto tg_yzz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 173); 

                auto tg_yzz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 174); 

                auto tg_yzz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 175); 

                auto tg_yzz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 176); 

                auto tg_yzz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 177); 

                auto tg_yzz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 178); 

                auto tg_yzz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 179); 

                auto tg_yzz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 180); 

                auto tg_yzz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 181); 

                auto tg_yzz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 182); 

                auto tg_yzz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 183); 

                auto tg_yzz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 184); 

                auto tg_yzz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 185); 

                auto tg_yzz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 186); 

                auto tg_yzz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 187); 

                auto tg_yzz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 188); 

                auto tg_zzz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 189); 

                auto tg_zzz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 190); 

                auto tg_zzz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 191); 

                auto tg_zzz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 192); 

                auto tg_zzz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 193); 

                auto tg_zzz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 194); 

                auto tg_zzz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 195); 

                auto tg_zzz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 196); 

                auto tg_zzz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 197); 

                auto tg_zzz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 198); 

                auto tg_zzz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 199); 

                auto tg_zzz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 200); 

                auto tg_zzz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 201); 

                auto tg_zzz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 202); 

                auto tg_zzz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 203); 

                auto tg_zzz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 204); 

                auto tg_zzz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 205); 

                auto tg_zzz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 206); 

                auto tg_zzz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 207); 

                auto tg_zzz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 208); 

                auto tg_zzz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 209); 

                auto tg_yyy_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 126); 

                auto tg_yyy_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 127); 

                auto tg_yyy_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 128); 

                auto tg_yyy_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 129); 

                auto tg_yyy_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 130); 

                auto tg_yyy_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 131); 

                auto tg_yyy_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 132); 

                auto tg_yyy_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 133); 

                auto tg_yyy_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 134); 

                auto tg_yyy_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 135); 

                auto tg_yyy_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 136); 

                auto tg_yyy_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 137); 

                auto tg_yyy_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 138); 

                auto tg_yyy_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 139); 

                auto tg_yyy_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 140); 

                auto tg_yyy_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 141); 

                auto tg_yyy_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 142); 

                auto tg_yyy_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 143); 

                auto tg_yyy_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 144); 

                auto tg_yyy_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 145); 

                auto tg_yyy_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 146); 

                auto tg_yyz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 147); 

                auto tg_yyz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 148); 

                auto tg_yyz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 149); 

                auto tg_yyz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 150); 

                auto tg_yyz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 151); 

                auto tg_yyz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 152); 

                auto tg_yyz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 158); 

                auto tg_yyz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 159); 

                auto tg_yyz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 160); 

                auto tg_yyz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 161); 

                auto tg_yyz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 162); 

                auto tg_yyz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 163); 

                auto tg_yyz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 164); 

                auto tg_yyz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 165); 

                auto tg_yyz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 166); 

                auto tg_yyz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 167); 

                auto tg_yzz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 168); 

                auto tg_yzz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 169); 

                auto tg_yzz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 170); 

                auto tg_yzz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 171); 

                auto tg_yzz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 172); 

                auto tg_yzz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 173); 

                auto tg_yzz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 174); 

                auto tg_yzz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 175); 

                auto tg_yzz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 176); 

                auto tg_yzz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 177); 

                auto tg_yzz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 178); 

                auto tg_yzz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 179); 

                auto tg_yzz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 180); 

                auto tg_yzz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 181); 

                auto tg_yzz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 182); 

                auto tg_yzz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 183); 

                auto tg_yzz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 184); 

                auto tg_yzz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 185); 

                auto tg_yzz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 186); 

                auto tg_yzz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 187); 

                auto tg_yzz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 188); 

                auto tg_zzz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 189); 

                auto tg_zzz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 190); 

                auto tg_zzz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 191); 

                auto tg_zzz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 192); 

                auto tg_zzz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 193); 

                auto tg_zzz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 194); 

                auto tg_zzz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 195); 

                auto tg_zzz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 196); 

                auto tg_zzz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 197); 

                auto tg_zzz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 198); 

                auto tg_zzz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 199); 

                auto tg_zzz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 200); 

                auto tg_zzz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 201); 

                auto tg_zzz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 202); 

                auto tg_zzz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 203); 

                auto tg_zzz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 204); 

                auto tg_zzz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 205); 

                auto tg_zzz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 206); 

                auto tg_zzz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 207); 

                auto tg_zzz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 208); 

                auto tg_zzz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 209); 

                auto tg_yy_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 63); 

                auto tg_yy_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 64); 

                auto tg_yy_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 65); 

                auto tg_yy_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 66); 

                auto tg_yy_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 67); 

                auto tg_yy_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 68); 

                auto tg_yy_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 69); 

                auto tg_yy_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 70); 

                auto tg_yy_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 71); 

                auto tg_yy_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 72); 

                auto tg_yy_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 73); 

                auto tg_yy_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 74); 

                auto tg_yy_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 75); 

                auto tg_yy_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 76); 

                auto tg_yy_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 77); 

                auto tg_yy_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 78); 

                auto tg_yy_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 79); 

                auto tg_yy_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 80); 

                auto tg_yy_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 81); 

                auto tg_yy_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 82); 

                auto tg_yy_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 83); 

                auto tg_yz_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 84); 

                auto tg_yz_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 85); 

                auto tg_yz_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 86); 

                auto tg_yz_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 87); 

                auto tg_yz_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 88); 

                auto tg_yz_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 89); 

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

                auto tg_yyy_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 90); 

                auto tg_yyy_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 91); 

                auto tg_yyy_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 92); 

                auto tg_yyy_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 93); 

                auto tg_yyy_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 94); 

                auto tg_yyy_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 95); 

                auto tg_yyy_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 96); 

                auto tg_yyy_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 97); 

                auto tg_yyy_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 98); 

                auto tg_yyy_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 99); 

                auto tg_yyy_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 100); 

                auto tg_yyy_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 101); 

                auto tg_yyy_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 102); 

                auto tg_yyy_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 103); 

                auto tg_yyy_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 104); 

                auto tg_yyz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 105); 

                auto tg_yyz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 106); 

                auto tg_yyz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 107); 

                auto tg_yyz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 116); 

                auto tg_yyz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 117); 

                auto tg_yyz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 118); 

                auto tg_yyz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 119); 

                auto tg_yzz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 120); 

                auto tg_yzz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 121); 

                auto tg_yzz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 122); 

                auto tg_yzz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 123); 

                auto tg_yzz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 124); 

                auto tg_yzz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 125); 

                auto tg_yzz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 126); 

                auto tg_yzz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 127); 

                auto tg_yzz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 128); 

                auto tg_yzz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 129); 

                auto tg_yzz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 130); 

                auto tg_yzz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 131); 

                auto tg_yzz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 132); 

                auto tg_yzz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 133); 

                auto tg_yzz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 134); 

                auto tg_zzz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 135); 

                auto tg_zzz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 136); 

                auto tg_zzz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 137); 

                auto tg_zzz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 138); 

                auto tg_zzz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 139); 

                auto tg_zzz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 140); 

                auto tg_zzz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 141); 

                auto tg_zzz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 142); 

                auto tg_zzz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 143); 

                auto tg_zzz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 144); 

                auto tg_zzz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 145); 

                auto tg_zzz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 146); 

                auto tg_zzz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 147); 

                auto tg_zzz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 148); 

                auto tg_zzz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 149); 

                // set up pointers to integrals

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

                // Batch of Integrals (158,237)

                #pragma omp simd aligned(fxn, fza, tg_xyyz_xyyyz_0, tg_xyyz_xyyzz_0, tg_xyyz_xyzzz_0, \
                                         tg_xyyz_xzzzz_0, tg_xyyz_yyyyy_0, tg_xyyz_yyyyz_0, tg_xyyz_yyyzz_0, tg_xyyz_yyzzz_0, \
                                         tg_xyyz_yzzzz_0, tg_xyyz_zzzzz_0, tg_xyzz_xxxxx_0, tg_xyzz_xxxxy_0, tg_xyzz_xxxxz_0, \
                                         tg_xyzz_xxxyy_0, tg_xyzz_xxxyz_0, tg_xyzz_xxxzz_0, tg_xyzz_xxyyy_0, tg_xyzz_xxyyz_0, \
                                         tg_xyzz_xxyzz_0, tg_xyzz_xxzzz_0, tg_xyzz_xyyyy_0, tg_xyzz_xyyyz_0, tg_xyzz_xyyzz_0, \
                                         tg_xyzz_xyzzz_0, tg_xyzz_xzzzz_0, tg_xyzz_yyyyy_0, tg_xyzz_yyyyz_0, tg_xyzz_yyyzz_0, \
                                         tg_xyzz_yyzzz_0, tg_xyzz_yzzzz_0, tg_xyzz_zzzzz_0, tg_xzzz_xxxxx_0, tg_xzzz_xxxxy_0, \
                                         tg_xzzz_xxxxz_0, tg_xzzz_xxxyy_0, tg_xzzz_xxxyz_0, tg_xzzz_xxxzz_0, tg_xzzz_xxyyy_0, \
                                         tg_xzzz_xxyyz_0, tg_xzzz_xxyzz_0, tg_xzzz_xxzzz_0, tg_xzzz_xyyyy_0, tg_xzzz_xyyyz_0, \
                                         tg_xzzz_xyyzz_0, tg_xzzz_xyzzz_0, tg_xzzz_xzzzz_0, tg_xzzz_yyyyy_0, tg_xzzz_yyyyz_0, \
                                         tg_xzzz_yyyzz_0, tg_xzzz_yyzzz_0, tg_xzzz_yzzzz_0, tg_xzzz_zzzzz_0, tg_yy_xxxxx_0, \
                                         tg_yy_xxxxx_1, tg_yy_xxxxy_0, tg_yy_xxxxy_1, tg_yy_xxxxz_0, tg_yy_xxxxz_1, \
                                         tg_yy_xxxyy_0, tg_yy_xxxyy_1, tg_yy_xxxyz_0, tg_yy_xxxyz_1, tg_yy_xxxzz_0, \
                                         tg_yy_xxxzz_1, tg_yy_xxyyy_0, tg_yy_xxyyy_1, tg_yy_xxyyz_0, tg_yy_xxyyz_1, \
                                         tg_yy_xxyzz_0, tg_yy_xxyzz_1, tg_yy_xxzzz_0, tg_yy_xxzzz_1, tg_yy_xyyyy_0, \
                                         tg_yy_xyyyy_1, tg_yy_xyyyz_0, tg_yy_xyyyz_1, tg_yy_xyyzz_0, tg_yy_xyyzz_1, \
                                         tg_yy_xyzzz_0, tg_yy_xyzzz_1, tg_yy_xzzzz_0, tg_yy_xzzzz_1, tg_yy_yyyyy_0, \
                                         tg_yy_yyyyy_1, tg_yy_yyyyz_0, tg_yy_yyyyz_1, tg_yy_yyyzz_0, tg_yy_yyyzz_1, \
                                         tg_yy_yyzzz_0, tg_yy_yyzzz_1, tg_yy_yzzzz_0, tg_yy_yzzzz_1, tg_yy_zzzzz_0, \
                                         tg_yy_zzzzz_1, tg_yyy_xxxx_1, tg_yyy_xxxxx_0, tg_yyy_xxxxx_1, tg_yyy_xxxxy_0, \
                                         tg_yyy_xxxxy_1, tg_yyy_xxxxz_0, tg_yyy_xxxxz_1, tg_yyy_xxxy_1, tg_yyy_xxxyy_0, \
                                         tg_yyy_xxxyy_1, tg_yyy_xxxyz_0, tg_yyy_xxxyz_1, tg_yyy_xxxz_1, tg_yyy_xxxzz_0, \
                                         tg_yyy_xxxzz_1, tg_yyy_xxyy_1, tg_yyy_xxyyy_0, tg_yyy_xxyyy_1, tg_yyy_xxyyz_0, \
                                         tg_yyy_xxyyz_1, tg_yyy_xxyz_1, tg_yyy_xxyzz_0, tg_yyy_xxyzz_1, tg_yyy_xxzz_1, \
                                         tg_yyy_xxzzz_0, tg_yyy_xxzzz_1, tg_yyy_xyyy_1, tg_yyy_xyyyy_0, tg_yyy_xyyyy_1, \
                                         tg_yyy_xyyyz_0, tg_yyy_xyyyz_1, tg_yyy_xyyz_1, tg_yyy_xyyzz_0, tg_yyy_xyyzz_1, \
                                         tg_yyy_xyzz_1, tg_yyy_xyzzz_0, tg_yyy_xyzzz_1, tg_yyy_xzzz_1, tg_yyy_xzzzz_0, \
                                         tg_yyy_xzzzz_1, tg_yyy_yyyy_1, tg_yyy_yyyyy_0, tg_yyy_yyyyy_1, tg_yyy_yyyyz_0, \
                                         tg_yyy_yyyyz_1, tg_yyy_yyyz_1, tg_yyy_yyyzz_0, tg_yyy_yyyzz_1, tg_yyy_yyzz_1, \
                                         tg_yyy_yyzzz_0, tg_yyy_yyzzz_1, tg_yyy_yzzz_1, tg_yyy_yzzzz_0, tg_yyy_yzzzz_1, \
                                         tg_yyy_zzzz_1, tg_yyy_zzzzz_0, tg_yyy_zzzzz_1, tg_yyyy_xxxxx_0, tg_yyyy_xxxxy_0, \
                                         tg_yyyy_xxxxz_0, tg_yyyy_xxxyy_0, tg_yyyy_xxxyz_0, tg_yyyy_xxxzz_0, tg_yyyy_xxyyy_0, \
                                         tg_yyyy_xxyyz_0, tg_yyyy_xxyzz_0, tg_yyyy_xxzzz_0, tg_yyyy_xyyyy_0, tg_yyyy_xyyyz_0, \
                                         tg_yyyy_xyyzz_0, tg_yyyy_xyzzz_0, tg_yyyy_xzzzz_0, tg_yyyy_yyyyy_0, tg_yyyy_yyyyz_0, \
                                         tg_yyyy_yyyzz_0, tg_yyyy_yyzzz_0, tg_yyyy_yzzzz_0, tg_yyyy_zzzzz_0, tg_yyyz_xxxxx_0, \
                                         tg_yyyz_xxxxy_0, tg_yyyz_xxxxz_0, tg_yyyz_xxxyy_0, tg_yyyz_xxxyz_0, tg_yyyz_xxxzz_0, \
                                         tg_yyz_xxxx_1, tg_yyz_xxxxx_0, tg_yyz_xxxxx_1, tg_yyz_xxxxy_0, tg_yyz_xxxxy_1, \
                                         tg_yyz_xxxxz_0, tg_yyz_xxxxz_1, tg_yyz_xxxy_1, tg_yyz_xxxyy_0, tg_yyz_xxxyy_1, \
                                         tg_yyz_xxxyz_0, tg_yyz_xxxyz_1, tg_yyz_xxxz_1, tg_yyz_xxxzz_0, tg_yyz_xxxzz_1, \
                                         tg_yyz_xyyyz_0, tg_yyz_xyyyz_1, tg_yyz_xyyzz_0, tg_yyz_xyyzz_1, tg_yyz_xyzzz_0, \
                                         tg_yyz_xyzzz_1, tg_yyz_xzzzz_0, tg_yyz_xzzzz_1, tg_yyz_yyyyy_0, tg_yyz_yyyyy_1, \
                                         tg_yyz_yyyyz_0, tg_yyz_yyyyz_1, tg_yyz_yyyz_1, tg_yyz_yyyzz_0, tg_yyz_yyyzz_1, \
                                         tg_yyz_yyzz_1, tg_yyz_yyzzz_0, tg_yyz_yyzzz_1, tg_yyz_yzzz_1, tg_yyz_yzzzz_0, \
                                         tg_yyz_yzzzz_1, tg_yyz_zzzz_1, tg_yyz_zzzzz_0, tg_yyz_zzzzz_1, tg_yz_xxxxx_0, \
                                         tg_yz_xxxxx_1, tg_yz_xxxxy_0, tg_yz_xxxxy_1, tg_yz_xxxxz_0, tg_yz_xxxxz_1, \
                                         tg_yz_xxxyy_0, tg_yz_xxxyy_1, tg_yz_xxxyz_0, tg_yz_xxxyz_1, tg_yz_xxxzz_0, \
                                         tg_yz_xxxzz_1, tg_yzz_xxxx_1, tg_yzz_xxxxx_0, tg_yzz_xxxxx_1, tg_yzz_xxxxy_0, \
                                         tg_yzz_xxxxy_1, tg_yzz_xxxxz_0, tg_yzz_xxxxz_1, tg_yzz_xxxy_1, tg_yzz_xxxyy_0, \
                                         tg_yzz_xxxyy_1, tg_yzz_xxxyz_0, tg_yzz_xxxyz_1, tg_yzz_xxxz_1, tg_yzz_xxxzz_0, \
                                         tg_yzz_xxxzz_1, tg_yzz_xxyy_1, tg_yzz_xxyyy_0, tg_yzz_xxyyy_1, tg_yzz_xxyyz_0, \
                                         tg_yzz_xxyyz_1, tg_yzz_xxyz_1, tg_yzz_xxyzz_0, tg_yzz_xxyzz_1, tg_yzz_xxzz_1, \
                                         tg_yzz_xxzzz_0, tg_yzz_xxzzz_1, tg_yzz_xyyy_1, tg_yzz_xyyyy_0, tg_yzz_xyyyy_1, \
                                         tg_yzz_xyyyz_0, tg_yzz_xyyyz_1, tg_yzz_xyyz_1, tg_yzz_xyyzz_0, tg_yzz_xyyzz_1, \
                                         tg_yzz_xyzz_1, tg_yzz_xyzzz_0, tg_yzz_xyzzz_1, tg_yzz_xzzz_1, tg_yzz_xzzzz_0, \
                                         tg_yzz_xzzzz_1, tg_yzz_yyyy_1, tg_yzz_yyyyy_0, tg_yzz_yyyyy_1, tg_yzz_yyyyz_0, \
                                         tg_yzz_yyyyz_1, tg_yzz_yyyz_1, tg_yzz_yyyzz_0, tg_yzz_yyyzz_1, tg_yzz_yyzz_1, \
                                         tg_yzz_yyzzz_0, tg_yzz_yyzzz_1, tg_yzz_yzzz_1, tg_yzz_yzzzz_0, tg_yzz_yzzzz_1, \
                                         tg_yzz_zzzz_1, tg_yzz_zzzzz_0, tg_yzz_zzzzz_1, tg_zzz_xxxx_1, tg_zzz_xxxxx_0, \
                                         tg_zzz_xxxxx_1, tg_zzz_xxxxy_0, tg_zzz_xxxxy_1, tg_zzz_xxxxz_0, tg_zzz_xxxxz_1, \
                                         tg_zzz_xxxy_1, tg_zzz_xxxyy_0, tg_zzz_xxxyy_1, tg_zzz_xxxyz_0, tg_zzz_xxxyz_1, \
                                         tg_zzz_xxxz_1, tg_zzz_xxxzz_0, tg_zzz_xxxzz_1, tg_zzz_xxyy_1, tg_zzz_xxyyy_0, \
                                         tg_zzz_xxyyy_1, tg_zzz_xxyyz_0, tg_zzz_xxyyz_1, tg_zzz_xxyz_1, tg_zzz_xxyzz_0, \
                                         tg_zzz_xxyzz_1, tg_zzz_xxzz_1, tg_zzz_xxzzz_0, tg_zzz_xxzzz_1, tg_zzz_xyyy_1, \
                                         tg_zzz_xyyyy_0, tg_zzz_xyyyy_1, tg_zzz_xyyyz_0, tg_zzz_xyyyz_1, tg_zzz_xyyz_1, \
                                         tg_zzz_xyyzz_0, tg_zzz_xyyzz_1, tg_zzz_xyzz_1, tg_zzz_xyzzz_0, tg_zzz_xyzzz_1, \
                                         tg_zzz_xzzz_1, tg_zzz_xzzzz_0, tg_zzz_xzzzz_1, tg_zzz_yyyy_1, tg_zzz_yyyyy_0, \
                                         tg_zzz_yyyyy_1, tg_zzz_yyyyz_0, tg_zzz_yyyyz_1, tg_zzz_yyyz_1, tg_zzz_yyyzz_0, \
                                         tg_zzz_yyyzz_1, tg_zzz_yyzz_1, tg_zzz_yyzzz_0, tg_zzz_yyzzz_1, tg_zzz_yzzz_1, \
                                         tg_zzz_yzzzz_0, tg_zzz_yzzzz_1, tg_zzz_zzzz_1, tg_zzz_zzzzz_0, tg_zzz_zzzzz_1, wp_x, \
                                         wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyyz_xyyyz_0[j] = pb_x * tg_yyz_xyyyz_0[j] + fr * tg_yyz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyz_1[j];

                    tg_xyyz_xyyzz_0[j] = pb_x * tg_yyz_xyyzz_0[j] + fr * tg_yyz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyzz_1[j];

                    tg_xyyz_xyzzz_0[j] = pb_x * tg_yyz_xyzzz_0[j] + fr * tg_yyz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yzzz_1[j];

                    tg_xyyz_xzzzz_0[j] = pb_x * tg_yyz_xzzzz_0[j] + fr * tg_yyz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zzzz_1[j];

                    tg_xyyz_yyyyy_0[j] = pb_x * tg_yyz_yyyyy_0[j] + fr * tg_yyz_yyyyy_1[j];

                    tg_xyyz_yyyyz_0[j] = pb_x * tg_yyz_yyyyz_0[j] + fr * tg_yyz_yyyyz_1[j];

                    tg_xyyz_yyyzz_0[j] = pb_x * tg_yyz_yyyzz_0[j] + fr * tg_yyz_yyyzz_1[j];

                    tg_xyyz_yyzzz_0[j] = pb_x * tg_yyz_yyzzz_0[j] + fr * tg_yyz_yyzzz_1[j];

                    tg_xyyz_yzzzz_0[j] = pb_x * tg_yyz_yzzzz_0[j] + fr * tg_yyz_yzzzz_1[j];

                    tg_xyyz_zzzzz_0[j] = pb_x * tg_yyz_zzzzz_0[j] + fr * tg_yyz_zzzzz_1[j];

                    tg_xyzz_xxxxx_0[j] = pb_x * tg_yzz_xxxxx_0[j] + fr * tg_yzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxx_1[j];

                    tg_xyzz_xxxxy_0[j] = pb_x * tg_yzz_xxxxy_0[j] + fr * tg_yzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxy_1[j];

                    tg_xyzz_xxxxz_0[j] = pb_x * tg_yzz_xxxxz_0[j] + fr * tg_yzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxz_1[j];

                    tg_xyzz_xxxyy_0[j] = pb_x * tg_yzz_xxxyy_0[j] + fr * tg_yzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyy_1[j];

                    tg_xyzz_xxxyz_0[j] = pb_x * tg_yzz_xxxyz_0[j] + fr * tg_yzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyz_1[j];

                    tg_xyzz_xxxzz_0[j] = pb_x * tg_yzz_xxxzz_0[j] + fr * tg_yzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxzz_1[j];

                    tg_xyzz_xxyyy_0[j] = pb_x * tg_yzz_xxyyy_0[j] + fr * tg_yzz_xxyyy_1[j] + fl1_fxn * tg_yzz_xyyy_1[j];

                    tg_xyzz_xxyyz_0[j] = pb_x * tg_yzz_xxyyz_0[j] + fr * tg_yzz_xxyyz_1[j] + fl1_fxn * tg_yzz_xyyz_1[j];

                    tg_xyzz_xxyzz_0[j] = pb_x * tg_yzz_xxyzz_0[j] + fr * tg_yzz_xxyzz_1[j] + fl1_fxn * tg_yzz_xyzz_1[j];

                    tg_xyzz_xxzzz_0[j] = pb_x * tg_yzz_xxzzz_0[j] + fr * tg_yzz_xxzzz_1[j] + fl1_fxn * tg_yzz_xzzz_1[j];

                    tg_xyzz_xyyyy_0[j] = pb_x * tg_yzz_xyyyy_0[j] + fr * tg_yzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyy_1[j];

                    tg_xyzz_xyyyz_0[j] = pb_x * tg_yzz_xyyyz_0[j] + fr * tg_yzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyz_1[j];

                    tg_xyzz_xyyzz_0[j] = pb_x * tg_yzz_xyyzz_0[j] + fr * tg_yzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyzz_1[j];

                    tg_xyzz_xyzzz_0[j] = pb_x * tg_yzz_xyzzz_0[j] + fr * tg_yzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yzzz_1[j];

                    tg_xyzz_xzzzz_0[j] = pb_x * tg_yzz_xzzzz_0[j] + fr * tg_yzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zzzz_1[j];

                    tg_xyzz_yyyyy_0[j] = pb_x * tg_yzz_yyyyy_0[j] + fr * tg_yzz_yyyyy_1[j];

                    tg_xyzz_yyyyz_0[j] = pb_x * tg_yzz_yyyyz_0[j] + fr * tg_yzz_yyyyz_1[j];

                    tg_xyzz_yyyzz_0[j] = pb_x * tg_yzz_yyyzz_0[j] + fr * tg_yzz_yyyzz_1[j];

                    tg_xyzz_yyzzz_0[j] = pb_x * tg_yzz_yyzzz_0[j] + fr * tg_yzz_yyzzz_1[j];

                    tg_xyzz_yzzzz_0[j] = pb_x * tg_yzz_yzzzz_0[j] + fr * tg_yzz_yzzzz_1[j];

                    tg_xyzz_zzzzz_0[j] = pb_x * tg_yzz_zzzzz_0[j] + fr * tg_yzz_zzzzz_1[j];

                    tg_xzzz_xxxxx_0[j] = pb_x * tg_zzz_xxxxx_0[j] + fr * tg_zzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxx_1[j];

                    tg_xzzz_xxxxy_0[j] = pb_x * tg_zzz_xxxxy_0[j] + fr * tg_zzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxy_1[j];

                    tg_xzzz_xxxxz_0[j] = pb_x * tg_zzz_xxxxz_0[j] + fr * tg_zzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxz_1[j];

                    tg_xzzz_xxxyy_0[j] = pb_x * tg_zzz_xxxyy_0[j] + fr * tg_zzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyy_1[j];

                    tg_xzzz_xxxyz_0[j] = pb_x * tg_zzz_xxxyz_0[j] + fr * tg_zzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyz_1[j];

                    tg_xzzz_xxxzz_0[j] = pb_x * tg_zzz_xxxzz_0[j] + fr * tg_zzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxzz_1[j];

                    tg_xzzz_xxyyy_0[j] = pb_x * tg_zzz_xxyyy_0[j] + fr * tg_zzz_xxyyy_1[j] + fl1_fxn * tg_zzz_xyyy_1[j];

                    tg_xzzz_xxyyz_0[j] = pb_x * tg_zzz_xxyyz_0[j] + fr * tg_zzz_xxyyz_1[j] + fl1_fxn * tg_zzz_xyyz_1[j];

                    tg_xzzz_xxyzz_0[j] = pb_x * tg_zzz_xxyzz_0[j] + fr * tg_zzz_xxyzz_1[j] + fl1_fxn * tg_zzz_xyzz_1[j];

                    tg_xzzz_xxzzz_0[j] = pb_x * tg_zzz_xxzzz_0[j] + fr * tg_zzz_xxzzz_1[j] + fl1_fxn * tg_zzz_xzzz_1[j];

                    tg_xzzz_xyyyy_0[j] = pb_x * tg_zzz_xyyyy_0[j] + fr * tg_zzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyy_1[j];

                    tg_xzzz_xyyyz_0[j] = pb_x * tg_zzz_xyyyz_0[j] + fr * tg_zzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyz_1[j];

                    tg_xzzz_xyyzz_0[j] = pb_x * tg_zzz_xyyzz_0[j] + fr * tg_zzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyzz_1[j];

                    tg_xzzz_xyzzz_0[j] = pb_x * tg_zzz_xyzzz_0[j] + fr * tg_zzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yzzz_1[j];

                    tg_xzzz_xzzzz_0[j] = pb_x * tg_zzz_xzzzz_0[j] + fr * tg_zzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzz_1[j];

                    tg_xzzz_yyyyy_0[j] = pb_x * tg_zzz_yyyyy_0[j] + fr * tg_zzz_yyyyy_1[j];

                    tg_xzzz_yyyyz_0[j] = pb_x * tg_zzz_yyyyz_0[j] + fr * tg_zzz_yyyyz_1[j];

                    tg_xzzz_yyyzz_0[j] = pb_x * tg_zzz_yyyzz_0[j] + fr * tg_zzz_yyyzz_1[j];

                    tg_xzzz_yyzzz_0[j] = pb_x * tg_zzz_yyzzz_0[j] + fr * tg_zzz_yyzzz_1[j];

                    tg_xzzz_yzzzz_0[j] = pb_x * tg_zzz_yzzzz_0[j] + fr * tg_zzz_yzzzz_1[j];

                    tg_xzzz_zzzzz_0[j] = pb_x * tg_zzz_zzzzz_0[j] + fr * tg_zzz_zzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyy_xxxxx_0[j] = pb_y * tg_yyy_xxxxx_0[j] + fr * tg_yyy_xxxxx_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxx_0[j] - tg_yy_xxxxx_1[j] * fl1_fza);

                    tg_yyyy_xxxxy_0[j] = pb_y * tg_yyy_xxxxy_0[j] + fr * tg_yyy_xxxxy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxy_0[j] - tg_yy_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxx_1[j];

                    tg_yyyy_xxxxz_0[j] = pb_y * tg_yyy_xxxxz_0[j] + fr * tg_yyy_xxxxz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxz_0[j] - tg_yy_xxxxz_1[j] * fl1_fza);

                    tg_yyyy_xxxyy_0[j] = pb_y * tg_yyy_xxxyy_0[j] + fr * tg_yyy_xxxyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyy_0[j] - tg_yy_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxxy_1[j];

                    tg_yyyy_xxxyz_0[j] = pb_y * tg_yyy_xxxyz_0[j] + fr * tg_yyy_xxxyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyz_0[j] - tg_yy_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxz_1[j];

                    tg_yyyy_xxxzz_0[j] = pb_y * tg_yyy_xxxzz_0[j] + fr * tg_yyy_xxxzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxzz_0[j] - tg_yy_xxxzz_1[j] * fl1_fza);

                    tg_yyyy_xxyyy_0[j] = pb_y * tg_yyy_xxyyy_0[j] + fr * tg_yyy_xxyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyy_0[j] - tg_yy_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xxyy_1[j];

                    tg_yyyy_xxyyz_0[j] = pb_y * tg_yyy_xxyyz_0[j] + fr * tg_yyy_xxyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyz_0[j] - tg_yy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxyz_1[j];

                    tg_yyyy_xxyzz_0[j] = pb_y * tg_yyy_xxyzz_0[j] + fr * tg_yyy_xxyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyzz_0[j] - tg_yy_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxzz_1[j];

                    tg_yyyy_xxzzz_0[j] = pb_y * tg_yyy_xxzzz_0[j] + fr * tg_yyy_xxzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxzzz_0[j] - tg_yy_xxzzz_1[j] * fl1_fza);

                    tg_yyyy_xyyyy_0[j] = pb_y * tg_yyy_xyyyy_0[j] + fr * tg_yyy_xyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyy_0[j] - tg_yy_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_xyyy_1[j];

                    tg_yyyy_xyyyz_0[j] = pb_y * tg_yyy_xyyyz_0[j] + fr * tg_yyy_xyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyz_0[j] - tg_yy_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xyyz_1[j];

                    tg_yyyy_xyyzz_0[j] = pb_y * tg_yyy_xyyzz_0[j] + fr * tg_yyy_xyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyzz_0[j] - tg_yy_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xyzz_1[j];

                    tg_yyyy_xyzzz_0[j] = pb_y * tg_yyy_xyzzz_0[j] + fr * tg_yyy_xyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyzzz_0[j] - tg_yy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xzzz_1[j];

                    tg_yyyy_xzzzz_0[j] = pb_y * tg_yyy_xzzzz_0[j] + fr * tg_yyy_xzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xzzzz_0[j] - tg_yy_xzzzz_1[j] * fl1_fza);

                    tg_yyyy_yyyyy_0[j] = pb_y * tg_yyy_yyyyy_0[j] + fr * tg_yyy_yyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyy_0[j] - tg_yy_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyy_yyyy_1[j];

                    tg_yyyy_yyyyz_0[j] = pb_y * tg_yyy_yyyyz_0[j] + fr * tg_yyy_yyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyz_0[j] - tg_yy_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_yyyz_1[j];

                    tg_yyyy_yyyzz_0[j] = pb_y * tg_yyy_yyyzz_0[j] + fr * tg_yyy_yyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyzz_0[j] - tg_yy_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_yyzz_1[j];

                    tg_yyyy_yyzzz_0[j] = pb_y * tg_yyy_yyzzz_0[j] + fr * tg_yyy_yyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyzzz_0[j] - tg_yy_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_yzzz_1[j];

                    tg_yyyy_yzzzz_0[j] = pb_y * tg_yyy_yzzzz_0[j] + fr * tg_yyy_yzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yzzzz_0[j] - tg_yy_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_zzzz_1[j];

                    tg_yyyy_zzzzz_0[j] = pb_y * tg_yyy_zzzzz_0[j] + fr * tg_yyy_zzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_zzzzz_0[j] - tg_yy_zzzzz_1[j] * fl1_fza);

                    tg_yyyz_xxxxx_0[j] = pb_y * tg_yyz_xxxxx_0[j] + fr * tg_yyz_xxxxx_1[j] + fl1_fx * (tg_yz_xxxxx_0[j] - tg_yz_xxxxx_1[j] * fl1_fza);

                    tg_yyyz_xxxxy_0[j] = pb_y * tg_yyz_xxxxy_0[j] + fr * tg_yyz_xxxxy_1[j] + fl1_fx * (tg_yz_xxxxy_0[j] - tg_yz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxx_1[j];

                    tg_yyyz_xxxxz_0[j] = pb_y * tg_yyz_xxxxz_0[j] + fr * tg_yyz_xxxxz_1[j] + fl1_fx * (tg_yz_xxxxz_0[j] - tg_yz_xxxxz_1[j] * fl1_fza);

                    tg_yyyz_xxxyy_0[j] = pb_y * tg_yyz_xxxyy_0[j] + fr * tg_yyz_xxxyy_1[j] + fl1_fx * (tg_yz_xxxyy_0[j] - tg_yz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxxy_1[j];

                    tg_yyyz_xxxyz_0[j] = pb_y * tg_yyz_xxxyz_0[j] + fr * tg_yyz_xxxyz_1[j] + fl1_fx * (tg_yz_xxxyz_0[j] - tg_yz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxz_1[j];

                    tg_yyyz_xxxzz_0[j] = pb_y * tg_yyz_xxxzz_0[j] + fr * tg_yyz_xxxzz_1[j] + fl1_fx * (tg_yz_xxxzz_0[j] - tg_yz_xxxzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSH_237_315(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (237,315)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 153); 

                auto tg_yyz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 154); 

                auto tg_yyz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 155); 

                auto tg_yyz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 156); 

                auto tg_yyz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 157); 

                auto tg_yyz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 158); 

                auto tg_yyz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 159); 

                auto tg_yyz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 160); 

                auto tg_yyz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 161); 

                auto tg_yyz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 162); 

                auto tg_yyz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 163); 

                auto tg_yyz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 164); 

                auto tg_yyz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 165); 

                auto tg_yyz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 166); 

                auto tg_yyz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 167); 

                auto tg_yzz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 168); 

                auto tg_yzz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 169); 

                auto tg_yzz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 170); 

                auto tg_yzz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 171); 

                auto tg_yzz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 172); 

                auto tg_yzz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 173); 

                auto tg_yzz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 174); 

                auto tg_yzz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 175); 

                auto tg_yzz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 176); 

                auto tg_yzz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 177); 

                auto tg_yzz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 178); 

                auto tg_yzz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 179); 

                auto tg_yzz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 180); 

                auto tg_yzz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 181); 

                auto tg_yzz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 182); 

                auto tg_yzz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 183); 

                auto tg_yzz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 184); 

                auto tg_yzz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 185); 

                auto tg_yzz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 186); 

                auto tg_yzz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 187); 

                auto tg_yzz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 188); 

                auto tg_zzz_xxxxx_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 189); 

                auto tg_zzz_xxxxy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 190); 

                auto tg_zzz_xxxxz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 191); 

                auto tg_zzz_xxxyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 192); 

                auto tg_zzz_xxxyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 193); 

                auto tg_zzz_xxxzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 194); 

                auto tg_zzz_xxyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 195); 

                auto tg_zzz_xxyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 196); 

                auto tg_zzz_xxyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 197); 

                auto tg_zzz_xxzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 198); 

                auto tg_zzz_xyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 199); 

                auto tg_zzz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 200); 

                auto tg_zzz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 201); 

                auto tg_zzz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 202); 

                auto tg_zzz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 203); 

                auto tg_zzz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 204); 

                auto tg_zzz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 205); 

                auto tg_zzz_yyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 206); 

                auto tg_zzz_yyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 207); 

                auto tg_zzz_yzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 208); 

                auto tg_zzz_zzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 209); 

                auto tg_yyz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 153); 

                auto tg_yyz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 154); 

                auto tg_yyz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 155); 

                auto tg_yyz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 156); 

                auto tg_yyz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 157); 

                auto tg_yyz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 158); 

                auto tg_yyz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 159); 

                auto tg_yyz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 160); 

                auto tg_yyz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 161); 

                auto tg_yyz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 162); 

                auto tg_yyz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 163); 

                auto tg_yyz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 164); 

                auto tg_yyz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 165); 

                auto tg_yyz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 166); 

                auto tg_yyz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 167); 

                auto tg_yzz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 168); 

                auto tg_yzz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 169); 

                auto tg_yzz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 170); 

                auto tg_yzz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 171); 

                auto tg_yzz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 172); 

                auto tg_yzz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 173); 

                auto tg_yzz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 174); 

                auto tg_yzz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 175); 

                auto tg_yzz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 176); 

                auto tg_yzz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 177); 

                auto tg_yzz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 178); 

                auto tg_yzz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 179); 

                auto tg_yzz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 180); 

                auto tg_yzz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 181); 

                auto tg_yzz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 182); 

                auto tg_yzz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 183); 

                auto tg_yzz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 184); 

                auto tg_yzz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 185); 

                auto tg_yzz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 186); 

                auto tg_yzz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 187); 

                auto tg_yzz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 188); 

                auto tg_zzz_xxxxx_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 189); 

                auto tg_zzz_xxxxy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 190); 

                auto tg_zzz_xxxxz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 191); 

                auto tg_zzz_xxxyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 192); 

                auto tg_zzz_xxxyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 193); 

                auto tg_zzz_xxxzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 194); 

                auto tg_zzz_xxyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 195); 

                auto tg_zzz_xxyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 196); 

                auto tg_zzz_xxyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 197); 

                auto tg_zzz_xxzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 198); 

                auto tg_zzz_xyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 199); 

                auto tg_zzz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 200); 

                auto tg_zzz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 201); 

                auto tg_zzz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 202); 

                auto tg_zzz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 203); 

                auto tg_zzz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 204); 

                auto tg_zzz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 205); 

                auto tg_zzz_yyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 206); 

                auto tg_zzz_yyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 207); 

                auto tg_zzz_yzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 208); 

                auto tg_zzz_zzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 209); 

                auto tg_yz_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 90); 

                auto tg_yz_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 91); 

                auto tg_yz_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 92); 

                auto tg_yz_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 93); 

                auto tg_yz_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 94); 

                auto tg_yz_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 95); 

                auto tg_yz_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 96); 

                auto tg_yz_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 97); 

                auto tg_yz_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 98); 

                auto tg_yz_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 99); 

                auto tg_yz_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 100); 

                auto tg_yz_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 101); 

                auto tg_yz_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 102); 

                auto tg_yz_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 103); 

                auto tg_yz_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 104); 

                auto tg_zz_xxxxx_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 105); 

                auto tg_zz_xxxxy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 106); 

                auto tg_zz_xxxxz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 107); 

                auto tg_zz_xxxyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 108); 

                auto tg_zz_xxxyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 109); 

                auto tg_zz_xxxzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 110); 

                auto tg_zz_xxyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 111); 

                auto tg_zz_xxyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 112); 

                auto tg_zz_xxyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 113); 

                auto tg_zz_xxzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 114); 

                auto tg_zz_xyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 115); 

                auto tg_zz_xyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 116); 

                auto tg_zz_xyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 117); 

                auto tg_zz_xyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 118); 

                auto tg_zz_xzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 119); 

                auto tg_zz_yyyyy_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 120); 

                auto tg_zz_yyyyz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 121); 

                auto tg_zz_yyyzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 122); 

                auto tg_zz_yyzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 123); 

                auto tg_zz_yzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 124); 

                auto tg_zz_zzzzz_0 = primBuffer[pidx_g_2_5_m0].data(126 * idx + 125); 

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

                auto tg_yyz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 108); 

                auto tg_yyz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 109); 

                auto tg_yyz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 110); 

                auto tg_yyz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 111); 

                auto tg_yyz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 112); 

                auto tg_yyz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 113); 

                auto tg_yyz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 114); 

                auto tg_yyz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 115); 

                auto tg_yyz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 116); 

                auto tg_yyz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 117); 

                auto tg_yyz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 118); 

                auto tg_yyz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 119); 

                auto tg_yzz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 120); 

                auto tg_yzz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 121); 

                auto tg_yzz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 122); 

                auto tg_yzz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 123); 

                auto tg_yzz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 124); 

                auto tg_yzz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 125); 

                auto tg_yzz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 126); 

                auto tg_yzz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 127); 

                auto tg_yzz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 128); 

                auto tg_yzz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 129); 

                auto tg_yzz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 130); 

                auto tg_yzz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 131); 

                auto tg_yzz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 132); 

                auto tg_yzz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 133); 

                auto tg_yzz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 134); 

                auto tg_zzz_xxxx_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 135); 

                auto tg_zzz_xxxy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 136); 

                auto tg_zzz_xxxz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 137); 

                auto tg_zzz_xxyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 138); 

                auto tg_zzz_xxyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 139); 

                auto tg_zzz_xxzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 140); 

                auto tg_zzz_xyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 141); 

                auto tg_zzz_xyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 142); 

                auto tg_zzz_xyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 143); 

                auto tg_zzz_xzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 144); 

                auto tg_zzz_yyyy_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 145); 

                auto tg_zzz_yyyz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 146); 

                auto tg_zzz_yyzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 147); 

                auto tg_zzz_yzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 148); 

                auto tg_zzz_zzzz_1 = primBuffer[pidx_g_3_4_m1].data(150 * idx + 149); 

                // set up pointers to integrals

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

                // Batch of Integrals (237,315)

                #pragma omp simd aligned(fxn, fza, tg_yyyz_xxyyy_0, tg_yyyz_xxyyz_0, tg_yyyz_xxyzz_0, \
                                         tg_yyyz_xxzzz_0, tg_yyyz_xyyyy_0, tg_yyyz_xyyyz_0, tg_yyyz_xyyzz_0, tg_yyyz_xyzzz_0, \
                                         tg_yyyz_xzzzz_0, tg_yyyz_yyyyy_0, tg_yyyz_yyyyz_0, tg_yyyz_yyyzz_0, tg_yyyz_yyzzz_0, \
                                         tg_yyyz_yzzzz_0, tg_yyyz_zzzzz_0, tg_yyz_xxyy_1, tg_yyz_xxyyy_0, tg_yyz_xxyyy_1, \
                                         tg_yyz_xxyyz_0, tg_yyz_xxyyz_1, tg_yyz_xxyz_1, tg_yyz_xxyzz_0, tg_yyz_xxyzz_1, \
                                         tg_yyz_xxzz_1, tg_yyz_xxzzz_0, tg_yyz_xxzzz_1, tg_yyz_xyyy_1, tg_yyz_xyyyy_0, \
                                         tg_yyz_xyyyy_1, tg_yyz_xyyyz_0, tg_yyz_xyyyz_1, tg_yyz_xyyz_1, tg_yyz_xyyzz_0, \
                                         tg_yyz_xyyzz_1, tg_yyz_xyzz_1, tg_yyz_xyzzz_0, tg_yyz_xyzzz_1, tg_yyz_xzzz_1, \
                                         tg_yyz_xzzzz_0, tg_yyz_xzzzz_1, tg_yyz_yyyy_1, tg_yyz_yyyyy_0, tg_yyz_yyyyy_1, \
                                         tg_yyz_yyyyz_0, tg_yyz_yyyyz_1, tg_yyz_yyyz_1, tg_yyz_yyyzz_0, tg_yyz_yyyzz_1, \
                                         tg_yyz_yyzz_1, tg_yyz_yyzzz_0, tg_yyz_yyzzz_1, tg_yyz_yzzz_1, tg_yyz_yzzzz_0, \
                                         tg_yyz_yzzzz_1, tg_yyz_zzzz_1, tg_yyz_zzzzz_0, tg_yyz_zzzzz_1, tg_yyzz_xxxxx_0, \
                                         tg_yyzz_xxxxy_0, tg_yyzz_xxxxz_0, tg_yyzz_xxxyy_0, tg_yyzz_xxxyz_0, tg_yyzz_xxxzz_0, \
                                         tg_yyzz_xxyyy_0, tg_yyzz_xxyyz_0, tg_yyzz_xxyzz_0, tg_yyzz_xxzzz_0, tg_yyzz_xyyyy_0, \
                                         tg_yyzz_xyyyz_0, tg_yyzz_xyyzz_0, tg_yyzz_xyzzz_0, tg_yyzz_xzzzz_0, tg_yyzz_yyyyy_0, \
                                         tg_yyzz_yyyyz_0, tg_yyzz_yyyzz_0, tg_yyzz_yyzzz_0, tg_yyzz_yzzzz_0, tg_yyzz_zzzzz_0, \
                                         tg_yz_xxyyy_0, tg_yz_xxyyy_1, tg_yz_xxyyz_0, tg_yz_xxyyz_1, tg_yz_xxyzz_0, \
                                         tg_yz_xxyzz_1, tg_yz_xxzzz_0, tg_yz_xxzzz_1, tg_yz_xyyyy_0, tg_yz_xyyyy_1, \
                                         tg_yz_xyyyz_0, tg_yz_xyyyz_1, tg_yz_xyyzz_0, tg_yz_xyyzz_1, tg_yz_xyzzz_0, \
                                         tg_yz_xyzzz_1, tg_yz_xzzzz_0, tg_yz_xzzzz_1, tg_yz_yyyyy_0, tg_yz_yyyyy_1, \
                                         tg_yz_yyyyz_0, tg_yz_yyyyz_1, tg_yz_yyyzz_0, tg_yz_yyyzz_1, tg_yz_yyzzz_0, \
                                         tg_yz_yyzzz_1, tg_yz_yzzzz_0, tg_yz_yzzzz_1, tg_yz_zzzzz_0, tg_yz_zzzzz_1, \
                                         tg_yzz_xxxx_1, tg_yzz_xxxxx_0, tg_yzz_xxxxx_1, tg_yzz_xxxxy_0, tg_yzz_xxxxy_1, \
                                         tg_yzz_xxxxz_0, tg_yzz_xxxxz_1, tg_yzz_xxxy_1, tg_yzz_xxxyy_0, tg_yzz_xxxyy_1, \
                                         tg_yzz_xxxyz_0, tg_yzz_xxxyz_1, tg_yzz_xxxz_1, tg_yzz_xxxzz_0, tg_yzz_xxxzz_1, \
                                         tg_yzz_xxyy_1, tg_yzz_xxyyy_0, tg_yzz_xxyyy_1, tg_yzz_xxyyz_0, tg_yzz_xxyyz_1, \
                                         tg_yzz_xxyz_1, tg_yzz_xxyzz_0, tg_yzz_xxyzz_1, tg_yzz_xxzz_1, tg_yzz_xxzzz_0, \
                                         tg_yzz_xxzzz_1, tg_yzz_xyyy_1, tg_yzz_xyyyy_0, tg_yzz_xyyyy_1, tg_yzz_xyyyz_0, \
                                         tg_yzz_xyyyz_1, tg_yzz_xyyz_1, tg_yzz_xyyzz_0, tg_yzz_xyyzz_1, tg_yzz_xyzz_1, \
                                         tg_yzz_xyzzz_0, tg_yzz_xyzzz_1, tg_yzz_xzzz_1, tg_yzz_xzzzz_0, tg_yzz_xzzzz_1, \
                                         tg_yzz_yyyy_1, tg_yzz_yyyyy_0, tg_yzz_yyyyy_1, tg_yzz_yyyyz_0, tg_yzz_yyyyz_1, \
                                         tg_yzz_yyyz_1, tg_yzz_yyyzz_0, tg_yzz_yyyzz_1, tg_yzz_yyzz_1, tg_yzz_yyzzz_0, \
                                         tg_yzz_yyzzz_1, tg_yzz_yzzz_1, tg_yzz_yzzzz_0, tg_yzz_yzzzz_1, tg_yzz_zzzz_1, \
                                         tg_yzz_zzzzz_0, tg_yzz_zzzzz_1, tg_yzzz_xxxxx_0, tg_yzzz_xxxxy_0, tg_yzzz_xxxxz_0, \
                                         tg_yzzz_xxxyy_0, tg_yzzz_xxxyz_0, tg_yzzz_xxxzz_0, tg_yzzz_xxyyy_0, tg_yzzz_xxyyz_0, \
                                         tg_yzzz_xxyzz_0, tg_yzzz_xxzzz_0, tg_yzzz_xyyyy_0, tg_yzzz_xyyyz_0, tg_yzzz_xyyzz_0, \
                                         tg_yzzz_xyzzz_0, tg_yzzz_xzzzz_0, tg_yzzz_yyyyy_0, tg_yzzz_yyyyz_0, tg_yzzz_yyyzz_0, \
                                         tg_yzzz_yyzzz_0, tg_yzzz_yzzzz_0, tg_yzzz_zzzzz_0, tg_zz_xxxxx_0, tg_zz_xxxxx_1, \
                                         tg_zz_xxxxy_0, tg_zz_xxxxy_1, tg_zz_xxxxz_0, tg_zz_xxxxz_1, tg_zz_xxxyy_0, \
                                         tg_zz_xxxyy_1, tg_zz_xxxyz_0, tg_zz_xxxyz_1, tg_zz_xxxzz_0, tg_zz_xxxzz_1, \
                                         tg_zz_xxyyy_0, tg_zz_xxyyy_1, tg_zz_xxyyz_0, tg_zz_xxyyz_1, tg_zz_xxyzz_0, \
                                         tg_zz_xxyzz_1, tg_zz_xxzzz_0, tg_zz_xxzzz_1, tg_zz_xyyyy_0, tg_zz_xyyyy_1, \
                                         tg_zz_xyyyz_0, tg_zz_xyyyz_1, tg_zz_xyyzz_0, tg_zz_xyyzz_1, tg_zz_xyzzz_0, \
                                         tg_zz_xyzzz_1, tg_zz_xzzzz_0, tg_zz_xzzzz_1, tg_zz_yyyyy_0, tg_zz_yyyyy_1, \
                                         tg_zz_yyyyz_0, tg_zz_yyyyz_1, tg_zz_yyyzz_0, tg_zz_yyyzz_1, tg_zz_yyzzz_0, \
                                         tg_zz_yyzzz_1, tg_zz_yzzzz_0, tg_zz_yzzzz_1, tg_zz_zzzzz_0, tg_zz_zzzzz_1, \
                                         tg_zzz_xxxx_1, tg_zzz_xxxxx_0, tg_zzz_xxxxx_1, tg_zzz_xxxxy_0, tg_zzz_xxxxy_1, \
                                         tg_zzz_xxxxz_0, tg_zzz_xxxxz_1, tg_zzz_xxxy_1, tg_zzz_xxxyy_0, tg_zzz_xxxyy_1, \
                                         tg_zzz_xxxyz_0, tg_zzz_xxxyz_1, tg_zzz_xxxz_1, tg_zzz_xxxzz_0, tg_zzz_xxxzz_1, \
                                         tg_zzz_xxyy_1, tg_zzz_xxyyy_0, tg_zzz_xxyyy_1, tg_zzz_xxyyz_0, tg_zzz_xxyyz_1, \
                                         tg_zzz_xxyz_1, tg_zzz_xxyzz_0, tg_zzz_xxyzz_1, tg_zzz_xxzz_1, tg_zzz_xxzzz_0, \
                                         tg_zzz_xxzzz_1, tg_zzz_xyyy_1, tg_zzz_xyyyy_0, tg_zzz_xyyyy_1, tg_zzz_xyyyz_0, \
                                         tg_zzz_xyyyz_1, tg_zzz_xyyz_1, tg_zzz_xyyzz_0, tg_zzz_xyyzz_1, tg_zzz_xyzz_1, \
                                         tg_zzz_xyzzz_0, tg_zzz_xyzzz_1, tg_zzz_xzzz_1, tg_zzz_xzzzz_0, tg_zzz_xzzzz_1, \
                                         tg_zzz_yyyy_1, tg_zzz_yyyyy_0, tg_zzz_yyyyy_1, tg_zzz_yyyyz_0, tg_zzz_yyyyz_1, \
                                         tg_zzz_yyyz_1, tg_zzz_yyyzz_0, tg_zzz_yyyzz_1, tg_zzz_yyzz_1, tg_zzz_yyzzz_0, \
                                         tg_zzz_yyzzz_1, tg_zzz_yzzz_1, tg_zzz_yzzzz_0, tg_zzz_yzzzz_1, tg_zzz_zzzz_1, \
                                         tg_zzz_zzzzz_0, tg_zzz_zzzzz_1, tg_zzzz_xxxxx_0, tg_zzzz_xxxxy_0, tg_zzzz_xxxxz_0, \
                                         tg_zzzz_xxxyy_0, tg_zzzz_xxxyz_0, tg_zzzz_xxxzz_0, tg_zzzz_xxyyy_0, tg_zzzz_xxyyz_0, \
                                         tg_zzzz_xxyzz_0, tg_zzzz_xxzzz_0, tg_zzzz_xyyyy_0, tg_zzzz_xyyyz_0, tg_zzzz_xyyzz_0, \
                                         tg_zzzz_xyzzz_0, tg_zzzz_xzzzz_0, tg_zzzz_yyyyy_0, tg_zzzz_yyyyz_0, tg_zzzz_yyyzz_0, \
                                         tg_zzzz_yyzzz_0, tg_zzzz_yzzzz_0, tg_zzzz_zzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyz_xxyyy_0[j] = pb_y * tg_yyz_xxyyy_0[j] + fr * tg_yyz_xxyyy_1[j] + fl1_fx * (tg_yz_xxyyy_0[j] - tg_yz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xxyy_1[j];

                    tg_yyyz_xxyyz_0[j] = pb_y * tg_yyz_xxyyz_0[j] + fr * tg_yyz_xxyyz_1[j] + fl1_fx * (tg_yz_xxyyz_0[j] - tg_yz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxyz_1[j];

                    tg_yyyz_xxyzz_0[j] = pb_y * tg_yyz_xxyzz_0[j] + fr * tg_yyz_xxyzz_1[j] + fl1_fx * (tg_yz_xxyzz_0[j] - tg_yz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxzz_1[j];

                    tg_yyyz_xxzzz_0[j] = pb_y * tg_yyz_xxzzz_0[j] + fr * tg_yyz_xxzzz_1[j] + fl1_fx * (tg_yz_xxzzz_0[j] - tg_yz_xxzzz_1[j] * fl1_fza);

                    tg_yyyz_xyyyy_0[j] = pb_y * tg_yyz_xyyyy_0[j] + fr * tg_yyz_xyyyy_1[j] + fl1_fx * (tg_yz_xyyyy_0[j] - tg_yz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_xyyy_1[j];

                    tg_yyyz_xyyyz_0[j] = pb_y * tg_yyz_xyyyz_0[j] + fr * tg_yyz_xyyyz_1[j] + fl1_fx * (tg_yz_xyyyz_0[j] - tg_yz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xyyz_1[j];

                    tg_yyyz_xyyzz_0[j] = pb_y * tg_yyz_xyyzz_0[j] + fr * tg_yyz_xyyzz_1[j] + fl1_fx * (tg_yz_xyyzz_0[j] - tg_yz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xyzz_1[j];

                    tg_yyyz_xyzzz_0[j] = pb_y * tg_yyz_xyzzz_0[j] + fr * tg_yyz_xyzzz_1[j] + fl1_fx * (tg_yz_xyzzz_0[j] - tg_yz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xzzz_1[j];

                    tg_yyyz_xzzzz_0[j] = pb_y * tg_yyz_xzzzz_0[j] + fr * tg_yyz_xzzzz_1[j] + fl1_fx * (tg_yz_xzzzz_0[j] - tg_yz_xzzzz_1[j] * fl1_fza);

                    tg_yyyz_yyyyy_0[j] = pb_y * tg_yyz_yyyyy_0[j] + fr * tg_yyz_yyyyy_1[j] + fl1_fx * (tg_yz_yyyyy_0[j] - tg_yz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyz_yyyy_1[j];

                    tg_yyyz_yyyyz_0[j] = pb_y * tg_yyz_yyyyz_0[j] + fr * tg_yyz_yyyyz_1[j] + fl1_fx * (tg_yz_yyyyz_0[j] - tg_yz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_yyyz_1[j];

                    tg_yyyz_yyyzz_0[j] = pb_y * tg_yyz_yyyzz_0[j] + fr * tg_yyz_yyyzz_1[j] + fl1_fx * (tg_yz_yyyzz_0[j] - tg_yz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_yyzz_1[j];

                    tg_yyyz_yyzzz_0[j] = pb_y * tg_yyz_yyzzz_0[j] + fr * tg_yyz_yyzzz_1[j] + fl1_fx * (tg_yz_yyzzz_0[j] - tg_yz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_yzzz_1[j];

                    tg_yyyz_yzzzz_0[j] = pb_y * tg_yyz_yzzzz_0[j] + fr * tg_yyz_yzzzz_1[j] + fl1_fx * (tg_yz_yzzzz_0[j] - tg_yz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_zzzz_1[j];

                    tg_yyyz_zzzzz_0[j] = pb_y * tg_yyz_zzzzz_0[j] + fr * tg_yyz_zzzzz_1[j] + fl1_fx * (tg_yz_zzzzz_0[j] - tg_yz_zzzzz_1[j] * fl1_fza);

                    tg_yyzz_xxxxx_0[j] = pb_y * tg_yzz_xxxxx_0[j] + fr * tg_yzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxx_0[j] - tg_zz_xxxxx_1[j] * fl1_fza);

                    tg_yyzz_xxxxy_0[j] = pb_y * tg_yzz_xxxxy_0[j] + fr * tg_yzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxy_0[j] - tg_zz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxx_1[j];

                    tg_yyzz_xxxxz_0[j] = pb_y * tg_yzz_xxxxz_0[j] + fr * tg_yzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxz_0[j] - tg_zz_xxxxz_1[j] * fl1_fza);

                    tg_yyzz_xxxyy_0[j] = pb_y * tg_yzz_xxxyy_0[j] + fr * tg_yzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyy_0[j] - tg_zz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxxy_1[j];

                    tg_yyzz_xxxyz_0[j] = pb_y * tg_yzz_xxxyz_0[j] + fr * tg_yzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyz_0[j] - tg_zz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxz_1[j];

                    tg_yyzz_xxxzz_0[j] = pb_y * tg_yzz_xxxzz_0[j] + fr * tg_yzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxzz_0[j] - tg_zz_xxxzz_1[j] * fl1_fza);

                    tg_yyzz_xxyyy_0[j] = pb_y * tg_yzz_xxyyy_0[j] + fr * tg_yzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyy_0[j] - tg_zz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xxyy_1[j];

                    tg_yyzz_xxyyz_0[j] = pb_y * tg_yzz_xxyyz_0[j] + fr * tg_yzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyz_0[j] - tg_zz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxyz_1[j];

                    tg_yyzz_xxyzz_0[j] = pb_y * tg_yzz_xxyzz_0[j] + fr * tg_yzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyzz_0[j] - tg_zz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxzz_1[j];

                    tg_yyzz_xxzzz_0[j] = pb_y * tg_yzz_xxzzz_0[j] + fr * tg_yzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxzzz_0[j] - tg_zz_xxzzz_1[j] * fl1_fza);

                    tg_yyzz_xyyyy_0[j] = pb_y * tg_yzz_xyyyy_0[j] + fr * tg_yzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyy_0[j] - tg_zz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_xyyy_1[j];

                    tg_yyzz_xyyyz_0[j] = pb_y * tg_yzz_xyyyz_0[j] + fr * tg_yzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyz_0[j] - tg_zz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xyyz_1[j];

                    tg_yyzz_xyyzz_0[j] = pb_y * tg_yzz_xyyzz_0[j] + fr * tg_yzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyzz_0[j] - tg_zz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xyzz_1[j];

                    tg_yyzz_xyzzz_0[j] = pb_y * tg_yzz_xyzzz_0[j] + fr * tg_yzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyzzz_0[j] - tg_zz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xzzz_1[j];

                    tg_yyzz_xzzzz_0[j] = pb_y * tg_yzz_xzzzz_0[j] + fr * tg_yzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzzzz_0[j] - tg_zz_xzzzz_1[j] * fl1_fza);

                    tg_yyzz_yyyyy_0[j] = pb_y * tg_yzz_yyyyy_0[j] + fr * tg_yzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyy_0[j] - tg_zz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzz_yyyy_1[j];

                    tg_yyzz_yyyyz_0[j] = pb_y * tg_yzz_yyyyz_0[j] + fr * tg_yzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyz_0[j] - tg_zz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_yyyz_1[j];

                    tg_yyzz_yyyzz_0[j] = pb_y * tg_yzz_yyyzz_0[j] + fr * tg_yzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyzz_0[j] - tg_zz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_yyzz_1[j];

                    tg_yyzz_yyzzz_0[j] = pb_y * tg_yzz_yyzzz_0[j] + fr * tg_yzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyzzz_0[j] - tg_zz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_yzzz_1[j];

                    tg_yyzz_yzzzz_0[j] = pb_y * tg_yzz_yzzzz_0[j] + fr * tg_yzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzzzz_0[j] - tg_zz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_zzzz_1[j];

                    tg_yyzz_zzzzz_0[j] = pb_y * tg_yzz_zzzzz_0[j] + fr * tg_yzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzzzz_0[j] - tg_zz_zzzzz_1[j] * fl1_fza);

                    tg_yzzz_xxxxx_0[j] = pb_y * tg_zzz_xxxxx_0[j] + fr * tg_zzz_xxxxx_1[j];

                    tg_yzzz_xxxxy_0[j] = pb_y * tg_zzz_xxxxy_0[j] + fr * tg_zzz_xxxxy_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxx_1[j];

                    tg_yzzz_xxxxz_0[j] = pb_y * tg_zzz_xxxxz_0[j] + fr * tg_zzz_xxxxz_1[j];

                    tg_yzzz_xxxyy_0[j] = pb_y * tg_zzz_xxxyy_0[j] + fr * tg_zzz_xxxyy_1[j] + fl1_fxn * tg_zzz_xxxy_1[j];

                    tg_yzzz_xxxyz_0[j] = pb_y * tg_zzz_xxxyz_0[j] + fr * tg_zzz_xxxyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxz_1[j];

                    tg_yzzz_xxxzz_0[j] = pb_y * tg_zzz_xxxzz_0[j] + fr * tg_zzz_xxxzz_1[j];

                    tg_yzzz_xxyyy_0[j] = pb_y * tg_zzz_xxyyy_0[j] + fr * tg_zzz_xxyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyy_1[j];

                    tg_yzzz_xxyyz_0[j] = pb_y * tg_zzz_xxyyz_0[j] + fr * tg_zzz_xxyyz_1[j] + fl1_fxn * tg_zzz_xxyz_1[j];

                    tg_yzzz_xxyzz_0[j] = pb_y * tg_zzz_xxyzz_0[j] + fr * tg_zzz_xxyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxzz_1[j];

                    tg_yzzz_xxzzz_0[j] = pb_y * tg_zzz_xxzzz_0[j] + fr * tg_zzz_xxzzz_1[j];

                    tg_yzzz_xyyyy_0[j] = pb_y * tg_zzz_xyyyy_0[j] + fr * tg_zzz_xyyyy_1[j] + 2.0 * fl1_fxn * tg_zzz_xyyy_1[j];

                    tg_yzzz_xyyyz_0[j] = pb_y * tg_zzz_xyyyz_0[j] + fr * tg_zzz_xyyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xyyz_1[j];

                    tg_yzzz_xyyzz_0[j] = pb_y * tg_zzz_xyyzz_0[j] + fr * tg_zzz_xyyzz_1[j] + fl1_fxn * tg_zzz_xyzz_1[j];

                    tg_yzzz_xyzzz_0[j] = pb_y * tg_zzz_xyzzz_0[j] + fr * tg_zzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xzzz_1[j];

                    tg_yzzz_xzzzz_0[j] = pb_y * tg_zzz_xzzzz_0[j] + fr * tg_zzz_xzzzz_1[j];

                    tg_yzzz_yyyyy_0[j] = pb_y * tg_zzz_yyyyy_0[j] + fr * tg_zzz_yyyyy_1[j] + 2.5 * fl1_fxn * tg_zzz_yyyy_1[j];

                    tg_yzzz_yyyyz_0[j] = pb_y * tg_zzz_yyyyz_0[j] + fr * tg_zzz_yyyyz_1[j] + 2.0 * fl1_fxn * tg_zzz_yyyz_1[j];

                    tg_yzzz_yyyzz_0[j] = pb_y * tg_zzz_yyyzz_0[j] + fr * tg_zzz_yyyzz_1[j] + 1.5 * fl1_fxn * tg_zzz_yyzz_1[j];

                    tg_yzzz_yyzzz_0[j] = pb_y * tg_zzz_yyzzz_0[j] + fr * tg_zzz_yyzzz_1[j] + fl1_fxn * tg_zzz_yzzz_1[j];

                    tg_yzzz_yzzzz_0[j] = pb_y * tg_zzz_yzzzz_0[j] + fr * tg_zzz_yzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzz_1[j];

                    tg_yzzz_zzzzz_0[j] = pb_y * tg_zzz_zzzzz_0[j] + fr * tg_zzz_zzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzz_xxxxx_0[j] = pb_z * tg_zzz_xxxxx_0[j] + fr * tg_zzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxx_0[j] - tg_zz_xxxxx_1[j] * fl1_fza);

                    tg_zzzz_xxxxy_0[j] = pb_z * tg_zzz_xxxxy_0[j] + fr * tg_zzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxy_0[j] - tg_zz_xxxxy_1[j] * fl1_fza);

                    tg_zzzz_xxxxz_0[j] = pb_z * tg_zzz_xxxxz_0[j] + fr * tg_zzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxz_0[j] - tg_zz_xxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxx_1[j];

                    tg_zzzz_xxxyy_0[j] = pb_z * tg_zzz_xxxyy_0[j] + fr * tg_zzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyy_0[j] - tg_zz_xxxyy_1[j] * fl1_fza);

                    tg_zzzz_xxxyz_0[j] = pb_z * tg_zzz_xxxyz_0[j] + fr * tg_zzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyz_0[j] - tg_zz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxy_1[j];

                    tg_zzzz_xxxzz_0[j] = pb_z * tg_zzz_xxxzz_0[j] + fr * tg_zzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxzz_0[j] - tg_zz_xxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxxz_1[j];

                    tg_zzzz_xxyyy_0[j] = pb_z * tg_zzz_xxyyy_0[j] + fr * tg_zzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyy_0[j] - tg_zz_xxyyy_1[j] * fl1_fza);

                    tg_zzzz_xxyyz_0[j] = pb_z * tg_zzz_xxyyz_0[j] + fr * tg_zzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyz_0[j] - tg_zz_xxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxyy_1[j];

                    tg_zzzz_xxyzz_0[j] = pb_z * tg_zzz_xxyzz_0[j] + fr * tg_zzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyzz_0[j] - tg_zz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxyz_1[j];

                    tg_zzzz_xxzzz_0[j] = pb_z * tg_zzz_xxzzz_0[j] + fr * tg_zzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxzzz_0[j] - tg_zz_xxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xxzz_1[j];

                    tg_zzzz_xyyyy_0[j] = pb_z * tg_zzz_xyyyy_0[j] + fr * tg_zzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyy_0[j] - tg_zz_xyyyy_1[j] * fl1_fza);

                    tg_zzzz_xyyyz_0[j] = pb_z * tg_zzz_xyyyz_0[j] + fr * tg_zzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyz_0[j] - tg_zz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xyyy_1[j];

                    tg_zzzz_xyyzz_0[j] = pb_z * tg_zzz_xyyzz_0[j] + fr * tg_zzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyzz_0[j] - tg_zz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xyyz_1[j];

                    tg_zzzz_xyzzz_0[j] = pb_z * tg_zzz_xyzzz_0[j] + fr * tg_zzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyzzz_0[j] - tg_zz_xyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xyzz_1[j];

                    tg_zzzz_xzzzz_0[j] = pb_z * tg_zzz_xzzzz_0[j] + fr * tg_zzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xzzzz_0[j] - tg_zz_xzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_xzzz_1[j];

                    tg_zzzz_yyyyy_0[j] = pb_z * tg_zzz_yyyyy_0[j] + fr * tg_zzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyy_0[j] - tg_zz_yyyyy_1[j] * fl1_fza);

                    tg_zzzz_yyyyz_0[j] = pb_z * tg_zzz_yyyyz_0[j] + fr * tg_zzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyz_0[j] - tg_zz_yyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_yyyy_1[j];

                    tg_zzzz_yyyzz_0[j] = pb_z * tg_zzz_yyyzz_0[j] + fr * tg_zzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyzz_0[j] - tg_zz_yyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_yyyz_1[j];

                    tg_zzzz_yyzzz_0[j] = pb_z * tg_zzz_yyzzz_0[j] + fr * tg_zzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyzzz_0[j] - tg_zz_yyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_yyzz_1[j];

                    tg_zzzz_yzzzz_0[j] = pb_z * tg_zzz_yzzzz_0[j] + fr * tg_zzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yzzzz_0[j] - tg_zz_yzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_yzzz_1[j];

                    tg_zzzz_zzzzz_0[j] = pb_z * tg_zzz_zzzzz_0[j] + fr * tg_zzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_zzzzz_0[j] - tg_zz_zzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

