//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForFH.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSFSH(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSFSH_0_70(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSFSH_70_140(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSFSH_140_210(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSFSH_0_70(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,70)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xx_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx); 

                auto tg_xx_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 1); 

                auto tg_xx_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 2); 

                auto tg_xx_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 3); 

                auto tg_xx_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 4); 

                auto tg_xx_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 5); 

                auto tg_xx_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 6); 

                auto tg_xx_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 7); 

                auto tg_xx_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 8); 

                auto tg_xx_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 9); 

                auto tg_xx_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 10); 

                auto tg_xx_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 11); 

                auto tg_xx_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 12); 

                auto tg_xx_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 13); 

                auto tg_xx_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 14); 

                auto tg_xx_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 15); 

                auto tg_xx_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 16); 

                auto tg_xx_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 17); 

                auto tg_xx_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 18); 

                auto tg_xx_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 19); 

                auto tg_xx_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 20); 

                auto tg_xy_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 21); 

                auto tg_xy_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 22); 

                auto tg_xy_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 23); 

                auto tg_xy_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 24); 

                auto tg_xy_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 25); 

                auto tg_xy_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 26); 

                auto tg_xy_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 27); 

                auto tg_xy_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 28); 

                auto tg_xy_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 29); 

                auto tg_xy_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 30); 

                auto tg_xy_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 31); 

                auto tg_xy_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 32); 

                auto tg_xy_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 33); 

                auto tg_xy_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 34); 

                auto tg_xy_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 35); 

                auto tg_xy_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 36); 

                auto tg_xy_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 37); 

                auto tg_xy_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 38); 

                auto tg_xy_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 39); 

                auto tg_xy_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 40); 

                auto tg_xy_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 41); 

                auto tg_xz_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 42); 

                auto tg_xz_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 43); 

                auto tg_xz_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 44); 

                auto tg_xz_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 45); 

                auto tg_xz_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 46); 

                auto tg_xz_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 47); 

                auto tg_xz_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 48); 

                auto tg_xz_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 49); 

                auto tg_xz_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 50); 

                auto tg_xz_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 51); 

                auto tg_xz_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 52); 

                auto tg_xz_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 53); 

                auto tg_xz_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 54); 

                auto tg_xz_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 55); 

                auto tg_xz_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 56); 

                auto tg_xz_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 57); 

                auto tg_xz_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 58); 

                auto tg_xz_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 59); 

                auto tg_xz_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 60); 

                auto tg_xz_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 61); 

                auto tg_xz_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 62); 

                auto tg_yy_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 63); 

                auto tg_yy_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 64); 

                auto tg_yy_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 65); 

                auto tg_yy_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 66); 

                auto tg_yy_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 67); 

                auto tg_yy_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 68); 

                auto tg_yy_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 69); 

                auto tg_xx_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx); 

                auto tg_xx_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 1); 

                auto tg_xx_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 2); 

                auto tg_xx_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 3); 

                auto tg_xx_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 4); 

                auto tg_xx_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 5); 

                auto tg_xx_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 6); 

                auto tg_xx_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 7); 

                auto tg_xx_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 8); 

                auto tg_xx_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 9); 

                auto tg_xx_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 10); 

                auto tg_xx_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 11); 

                auto tg_xx_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 12); 

                auto tg_xx_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 13); 

                auto tg_xx_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 14); 

                auto tg_xx_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 15); 

                auto tg_xx_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 16); 

                auto tg_xx_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 17); 

                auto tg_xx_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 18); 

                auto tg_xx_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 19); 

                auto tg_xx_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 20); 

                auto tg_xy_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 21); 

                auto tg_xy_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 22); 

                auto tg_xy_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 23); 

                auto tg_xy_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 24); 

                auto tg_xy_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 25); 

                auto tg_xy_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 26); 

                auto tg_xy_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 27); 

                auto tg_xy_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 28); 

                auto tg_xy_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 29); 

                auto tg_xy_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 30); 

                auto tg_xy_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 31); 

                auto tg_xy_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 32); 

                auto tg_xy_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 33); 

                auto tg_xy_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 34); 

                auto tg_xy_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 35); 

                auto tg_xy_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 36); 

                auto tg_xy_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 37); 

                auto tg_xy_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 38); 

                auto tg_xy_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 39); 

                auto tg_xy_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 40); 

                auto tg_xy_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 41); 

                auto tg_xz_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 42); 

                auto tg_xz_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 43); 

                auto tg_xz_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 44); 

                auto tg_xz_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 45); 

                auto tg_xz_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 46); 

                auto tg_xz_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 47); 

                auto tg_xz_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 48); 

                auto tg_xz_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 49); 

                auto tg_xz_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 50); 

                auto tg_xz_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 51); 

                auto tg_xz_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 52); 

                auto tg_xz_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 53); 

                auto tg_xz_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 54); 

                auto tg_xz_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 55); 

                auto tg_xz_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 56); 

                auto tg_xz_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 57); 

                auto tg_xz_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 58); 

                auto tg_xz_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 59); 

                auto tg_xz_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 60); 

                auto tg_xz_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 61); 

                auto tg_xz_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 62); 

                auto tg_yy_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 63); 

                auto tg_yy_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 64); 

                auto tg_yy_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 65); 

                auto tg_yy_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 66); 

                auto tg_yy_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 67); 

                auto tg_yy_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 68); 

                auto tg_yy_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 69); 

                auto tg_x_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx); 

                auto tg_x_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 1); 

                auto tg_x_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 2); 

                auto tg_x_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 3); 

                auto tg_x_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 4); 

                auto tg_x_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 5); 

                auto tg_x_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 6); 

                auto tg_x_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 7); 

                auto tg_x_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 8); 

                auto tg_x_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 9); 

                auto tg_x_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 10); 

                auto tg_x_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 11); 

                auto tg_x_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 12); 

                auto tg_x_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 13); 

                auto tg_x_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 14); 

                auto tg_x_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 15); 

                auto tg_x_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 16); 

                auto tg_x_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 17); 

                auto tg_x_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 18); 

                auto tg_x_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 19); 

                auto tg_x_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 20); 

                auto tg_y_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 21); 

                auto tg_y_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 22); 

                auto tg_y_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 23); 

                auto tg_y_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 24); 

                auto tg_y_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 25); 

                auto tg_y_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 26); 

                auto tg_y_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 27); 

                auto tg_y_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 28); 

                auto tg_y_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 29); 

                auto tg_y_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 30); 

                auto tg_y_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 31); 

                auto tg_y_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 32); 

                auto tg_y_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 33); 

                auto tg_y_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 34); 

                auto tg_y_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 35); 

                auto tg_y_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 36); 

                auto tg_y_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 37); 

                auto tg_y_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 38); 

                auto tg_y_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 39); 

                auto tg_y_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 40); 

                auto tg_y_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 41); 

                auto tg_z_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 42); 

                auto tg_z_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 43); 

                auto tg_z_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 44); 

                auto tg_z_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 45); 

                auto tg_z_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 46); 

                auto tg_z_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 47); 

                auto tg_z_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 48); 

                auto tg_z_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 49); 

                auto tg_z_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 50); 

                auto tg_z_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 51); 

                auto tg_z_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 52); 

                auto tg_z_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 53); 

                auto tg_z_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 54); 

                auto tg_z_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 55); 

                auto tg_z_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 56); 

                auto tg_z_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 57); 

                auto tg_z_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 58); 

                auto tg_z_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 59); 

                auto tg_z_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 60); 

                auto tg_z_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 61); 

                auto tg_z_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 62); 

                auto tg_x_xxxxx_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx); 

                auto tg_x_xxxxy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 1); 

                auto tg_x_xxxxz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 2); 

                auto tg_x_xxxyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 3); 

                auto tg_x_xxxyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 4); 

                auto tg_x_xxxzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 5); 

                auto tg_x_xxyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 6); 

                auto tg_x_xxyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 7); 

                auto tg_x_xxyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 8); 

                auto tg_x_xxzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 9); 

                auto tg_x_xyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 10); 

                auto tg_x_xyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 11); 

                auto tg_x_xyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 12); 

                auto tg_x_xyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 13); 

                auto tg_x_xzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 14); 

                auto tg_x_yyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 15); 

                auto tg_x_yyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 16); 

                auto tg_x_yyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 17); 

                auto tg_x_yyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 18); 

                auto tg_x_yzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 19); 

                auto tg_x_zzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 20); 

                auto tg_y_xxxxx_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 21); 

                auto tg_y_xxxxy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 22); 

                auto tg_y_xxxxz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 23); 

                auto tg_y_xxxyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 24); 

                auto tg_y_xxxyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 25); 

                auto tg_y_xxxzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 26); 

                auto tg_y_xxyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 27); 

                auto tg_y_xxyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 28); 

                auto tg_y_xxyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 29); 

                auto tg_y_xxzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 30); 

                auto tg_y_xyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 31); 

                auto tg_y_xyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 32); 

                auto tg_y_xyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 33); 

                auto tg_y_xyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 34); 

                auto tg_y_xzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 35); 

                auto tg_y_yyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 36); 

                auto tg_y_yyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 37); 

                auto tg_y_yyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 38); 

                auto tg_y_yyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 39); 

                auto tg_y_yzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 40); 

                auto tg_y_zzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 41); 

                auto tg_z_xxxxx_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 42); 

                auto tg_z_xxxxy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 43); 

                auto tg_z_xxxxz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 44); 

                auto tg_z_xxxyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 45); 

                auto tg_z_xxxyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 46); 

                auto tg_z_xxxzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 47); 

                auto tg_z_xxyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 48); 

                auto tg_z_xxyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 49); 

                auto tg_z_xxyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 50); 

                auto tg_z_xxzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 51); 

                auto tg_z_xyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 52); 

                auto tg_z_xyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 53); 

                auto tg_z_xyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 54); 

                auto tg_z_xyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 55); 

                auto tg_z_xzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 56); 

                auto tg_z_yyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 57); 

                auto tg_z_yyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 58); 

                auto tg_z_yyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 59); 

                auto tg_z_yyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 60); 

                auto tg_z_yzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 61); 

                auto tg_z_zzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 62); 

                auto tg_xx_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx); 

                auto tg_xx_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 1); 

                auto tg_xx_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 2); 

                auto tg_xx_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 3); 

                auto tg_xx_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 4); 

                auto tg_xx_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 5); 

                auto tg_xx_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 6); 

                auto tg_xx_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 7); 

                auto tg_xx_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 8); 

                auto tg_xx_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 9); 

                auto tg_xx_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 10); 

                auto tg_xx_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 11); 

                auto tg_xx_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 12); 

                auto tg_xx_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 13); 

                auto tg_xx_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 14); 

                auto tg_xy_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 15); 

                auto tg_xy_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 16); 

                auto tg_xy_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 17); 

                auto tg_xy_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 18); 

                auto tg_xy_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 19); 

                auto tg_xy_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 20); 

                auto tg_xy_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 21); 

                auto tg_xy_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 22); 

                auto tg_xy_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 23); 

                auto tg_xy_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 24); 

                auto tg_xy_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 25); 

                auto tg_xy_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 26); 

                auto tg_xy_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 27); 

                auto tg_xy_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 28); 

                auto tg_xy_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 29); 

                auto tg_xz_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 30); 

                auto tg_xz_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 31); 

                auto tg_xz_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 32); 

                auto tg_xz_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 33); 

                auto tg_xz_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 34); 

                auto tg_xz_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 35); 

                auto tg_xz_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 36); 

                auto tg_xz_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 37); 

                auto tg_xz_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 38); 

                auto tg_xz_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 39); 

                auto tg_xz_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 40); 

                auto tg_xz_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 41); 

                auto tg_xz_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 42); 

                auto tg_xz_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 43); 

                auto tg_xz_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 44); 

                auto tg_yy_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 45); 

                auto tg_yy_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 46); 

                auto tg_yy_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 47); 

                auto tg_yy_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 48); 

                auto tg_yy_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 49); 

                auto tg_yy_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 50); 

                auto tg_yy_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 51); 

                // set up pointers to integrals

                auto tg_xxx_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx); 

                auto tg_xxx_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 1); 

                auto tg_xxx_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 2); 

                auto tg_xxx_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 3); 

                auto tg_xxx_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 4); 

                auto tg_xxx_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 5); 

                auto tg_xxx_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 6); 

                auto tg_xxx_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 7); 

                auto tg_xxx_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 8); 

                auto tg_xxx_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 9); 

                auto tg_xxx_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 10); 

                auto tg_xxx_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 11); 

                auto tg_xxx_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 12); 

                auto tg_xxx_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 13); 

                auto tg_xxx_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 14); 

                auto tg_xxx_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 15); 

                auto tg_xxx_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 16); 

                auto tg_xxx_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 17); 

                auto tg_xxx_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 18); 

                auto tg_xxx_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 19); 

                auto tg_xxx_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 20); 

                auto tg_xxy_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 21); 

                auto tg_xxy_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 22); 

                auto tg_xxy_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 23); 

                auto tg_xxy_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 24); 

                auto tg_xxy_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 25); 

                auto tg_xxy_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 26); 

                auto tg_xxy_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 27); 

                auto tg_xxy_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 28); 

                auto tg_xxy_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 29); 

                auto tg_xxy_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 30); 

                auto tg_xxy_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 31); 

                auto tg_xxy_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 32); 

                auto tg_xxy_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 33); 

                auto tg_xxy_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 34); 

                auto tg_xxy_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 35); 

                auto tg_xxy_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 36); 

                auto tg_xxy_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 37); 

                auto tg_xxy_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 38); 

                auto tg_xxy_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 39); 

                auto tg_xxy_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 40); 

                auto tg_xxy_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 41); 

                auto tg_xxz_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 42); 

                auto tg_xxz_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 43); 

                auto tg_xxz_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 44); 

                auto tg_xxz_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 45); 

                auto tg_xxz_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 46); 

                auto tg_xxz_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 47); 

                auto tg_xxz_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 48); 

                auto tg_xxz_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 49); 

                auto tg_xxz_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 50); 

                auto tg_xxz_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 51); 

                auto tg_xxz_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 52); 

                auto tg_xxz_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 53); 

                auto tg_xxz_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 54); 

                auto tg_xxz_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 55); 

                auto tg_xxz_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 56); 

                auto tg_xxz_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 57); 

                auto tg_xxz_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 58); 

                auto tg_xxz_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 59); 

                auto tg_xxz_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 60); 

                auto tg_xxz_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 61); 

                auto tg_xxz_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 62); 

                auto tg_xyy_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 63); 

                auto tg_xyy_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 64); 

                auto tg_xyy_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 65); 

                auto tg_xyy_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 66); 

                auto tg_xyy_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 67); 

                auto tg_xyy_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 68); 

                auto tg_xyy_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 69); 

                // Batch of Integrals (0,70)

                #pragma omp simd aligned(fxn, fza, tg_x_xxxxx_0, tg_x_xxxxx_1, tg_x_xxxxy_0, tg_x_xxxxy_1, \
                                         tg_x_xxxxz_0, tg_x_xxxxz_1, tg_x_xxxyy_0, tg_x_xxxyy_1, tg_x_xxxyz_0, tg_x_xxxyz_1, \
                                         tg_x_xxxzz_0, tg_x_xxxzz_1, tg_x_xxyyy_0, tg_x_xxyyy_1, tg_x_xxyyz_0, tg_x_xxyyz_1, \
                                         tg_x_xxyzz_0, tg_x_xxyzz_1, tg_x_xxzzz_0, tg_x_xxzzz_1, tg_x_xyyyy_0, tg_x_xyyyy_1, \
                                         tg_x_xyyyz_0, tg_x_xyyyz_1, tg_x_xyyzz_0, tg_x_xyyzz_1, tg_x_xyzzz_0, tg_x_xyzzz_1, \
                                         tg_x_xzzzz_0, tg_x_xzzzz_1, tg_x_yyyyy_0, tg_x_yyyyy_1, tg_x_yyyyz_0, tg_x_yyyyz_1, \
                                         tg_x_yyyzz_0, tg_x_yyyzz_1, tg_x_yyzzz_0, tg_x_yyzzz_1, tg_x_yzzzz_0, tg_x_yzzzz_1, \
                                         tg_x_zzzzz_0, tg_x_zzzzz_1, tg_xx_xxxx_1, tg_xx_xxxxx_0, tg_xx_xxxxx_1, \
                                         tg_xx_xxxxy_0, tg_xx_xxxxy_1, tg_xx_xxxxz_0, tg_xx_xxxxz_1, tg_xx_xxxy_1, \
                                         tg_xx_xxxyy_0, tg_xx_xxxyy_1, tg_xx_xxxyz_0, tg_xx_xxxyz_1, tg_xx_xxxz_1, \
                                         tg_xx_xxxzz_0, tg_xx_xxxzz_1, tg_xx_xxyy_1, tg_xx_xxyyy_0, tg_xx_xxyyy_1, \
                                         tg_xx_xxyyz_0, tg_xx_xxyyz_1, tg_xx_xxyz_1, tg_xx_xxyzz_0, tg_xx_xxyzz_1, \
                                         tg_xx_xxzz_1, tg_xx_xxzzz_0, tg_xx_xxzzz_1, tg_xx_xyyy_1, tg_xx_xyyyy_0, \
                                         tg_xx_xyyyy_1, tg_xx_xyyyz_0, tg_xx_xyyyz_1, tg_xx_xyyz_1, tg_xx_xyyzz_0, \
                                         tg_xx_xyyzz_1, tg_xx_xyzz_1, tg_xx_xyzzz_0, tg_xx_xyzzz_1, tg_xx_xzzz_1, \
                                         tg_xx_xzzzz_0, tg_xx_xzzzz_1, tg_xx_yyyy_1, tg_xx_yyyyy_0, tg_xx_yyyyy_1, \
                                         tg_xx_yyyyz_0, tg_xx_yyyyz_1, tg_xx_yyyz_1, tg_xx_yyyzz_0, tg_xx_yyyzz_1, \
                                         tg_xx_yyzz_1, tg_xx_yyzzz_0, tg_xx_yyzzz_1, tg_xx_yzzz_1, tg_xx_yzzzz_0, \
                                         tg_xx_yzzzz_1, tg_xx_zzzz_1, tg_xx_zzzzz_0, tg_xx_zzzzz_1, tg_xxx_xxxxx_0, \
                                         tg_xxx_xxxxy_0, tg_xxx_xxxxz_0, tg_xxx_xxxyy_0, tg_xxx_xxxyz_0, tg_xxx_xxxzz_0, \
                                         tg_xxx_xxyyy_0, tg_xxx_xxyyz_0, tg_xxx_xxyzz_0, tg_xxx_xxzzz_0, tg_xxx_xyyyy_0, \
                                         tg_xxx_xyyyz_0, tg_xxx_xyyzz_0, tg_xxx_xyzzz_0, tg_xxx_xzzzz_0, tg_xxx_yyyyy_0, \
                                         tg_xxx_yyyyz_0, tg_xxx_yyyzz_0, tg_xxx_yyzzz_0, tg_xxx_yzzzz_0, tg_xxx_zzzzz_0, \
                                         tg_xxy_xxxxx_0, tg_xxy_xxxxy_0, tg_xxy_xxxxz_0, tg_xxy_xxxyy_0, tg_xxy_xxxyz_0, \
                                         tg_xxy_xxxzz_0, tg_xxy_xxyyy_0, tg_xxy_xxyyz_0, tg_xxy_xxyzz_0, tg_xxy_xxzzz_0, \
                                         tg_xxy_xyyyy_0, tg_xxy_xyyyz_0, tg_xxy_xyyzz_0, tg_xxy_xyzzz_0, tg_xxy_xzzzz_0, \
                                         tg_xxy_yyyyy_0, tg_xxy_yyyyz_0, tg_xxy_yyyzz_0, tg_xxy_yyzzz_0, tg_xxy_yzzzz_0, \
                                         tg_xxy_zzzzz_0, tg_xxz_xxxxx_0, tg_xxz_xxxxy_0, tg_xxz_xxxxz_0, tg_xxz_xxxyy_0, \
                                         tg_xxz_xxxyz_0, tg_xxz_xxxzz_0, tg_xxz_xxyyy_0, tg_xxz_xxyyz_0, tg_xxz_xxyzz_0, \
                                         tg_xxz_xxzzz_0, tg_xxz_xyyyy_0, tg_xxz_xyyyz_0, tg_xxz_xyyzz_0, tg_xxz_xyzzz_0, \
                                         tg_xxz_xzzzz_0, tg_xxz_yyyyy_0, tg_xxz_yyyyz_0, tg_xxz_yyyzz_0, tg_xxz_yyzzz_0, \
                                         tg_xxz_yzzzz_0, tg_xxz_zzzzz_0, tg_xy_xxxx_1, tg_xy_xxxxx_0, tg_xy_xxxxx_1, \
                                         tg_xy_xxxxy_0, tg_xy_xxxxy_1, tg_xy_xxxxz_0, tg_xy_xxxxz_1, tg_xy_xxxy_1, \
                                         tg_xy_xxxyy_0, tg_xy_xxxyy_1, tg_xy_xxxyz_0, tg_xy_xxxyz_1, tg_xy_xxxz_1, \
                                         tg_xy_xxxzz_0, tg_xy_xxxzz_1, tg_xy_xxyy_1, tg_xy_xxyyy_0, tg_xy_xxyyy_1, \
                                         tg_xy_xxyyz_0, tg_xy_xxyyz_1, tg_xy_xxyz_1, tg_xy_xxyzz_0, tg_xy_xxyzz_1, \
                                         tg_xy_xxzz_1, tg_xy_xxzzz_0, tg_xy_xxzzz_1, tg_xy_xyyy_1, tg_xy_xyyyy_0, \
                                         tg_xy_xyyyy_1, tg_xy_xyyyz_0, tg_xy_xyyyz_1, tg_xy_xyyz_1, tg_xy_xyyzz_0, \
                                         tg_xy_xyyzz_1, tg_xy_xyzz_1, tg_xy_xyzzz_0, tg_xy_xyzzz_1, tg_xy_xzzz_1, \
                                         tg_xy_xzzzz_0, tg_xy_xzzzz_1, tg_xy_yyyy_1, tg_xy_yyyyy_0, tg_xy_yyyyy_1, \
                                         tg_xy_yyyyz_0, tg_xy_yyyyz_1, tg_xy_yyyz_1, tg_xy_yyyzz_0, tg_xy_yyyzz_1, \
                                         tg_xy_yyzz_1, tg_xy_yyzzz_0, tg_xy_yyzzz_1, tg_xy_yzzz_1, tg_xy_yzzzz_0, \
                                         tg_xy_yzzzz_1, tg_xy_zzzz_1, tg_xy_zzzzz_0, tg_xy_zzzzz_1, tg_xyy_xxxxx_0, \
                                         tg_xyy_xxxxy_0, tg_xyy_xxxxz_0, tg_xyy_xxxyy_0, tg_xyy_xxxyz_0, tg_xyy_xxxzz_0, \
                                         tg_xyy_xxyyy_0, tg_xz_xxxx_1, tg_xz_xxxxx_0, tg_xz_xxxxx_1, tg_xz_xxxxy_0, \
                                         tg_xz_xxxxy_1, tg_xz_xxxxz_0, tg_xz_xxxxz_1, tg_xz_xxxy_1, tg_xz_xxxyy_0, \
                                         tg_xz_xxxyy_1, tg_xz_xxxyz_0, tg_xz_xxxyz_1, tg_xz_xxxz_1, tg_xz_xxxzz_0, \
                                         tg_xz_xxxzz_1, tg_xz_xxyy_1, tg_xz_xxyyy_0, tg_xz_xxyyy_1, tg_xz_xxyyz_0, \
                                         tg_xz_xxyyz_1, tg_xz_xxyz_1, tg_xz_xxyzz_0, tg_xz_xxyzz_1, tg_xz_xxzz_1, \
                                         tg_xz_xxzzz_0, tg_xz_xxzzz_1, tg_xz_xyyy_1, tg_xz_xyyyy_0, tg_xz_xyyyy_1, \
                                         tg_xz_xyyyz_0, tg_xz_xyyyz_1, tg_xz_xyyz_1, tg_xz_xyyzz_0, tg_xz_xyyzz_1, \
                                         tg_xz_xyzz_1, tg_xz_xyzzz_0, tg_xz_xyzzz_1, tg_xz_xzzz_1, tg_xz_xzzzz_0, \
                                         tg_xz_xzzzz_1, tg_xz_yyyy_1, tg_xz_yyyyy_0, tg_xz_yyyyy_1, tg_xz_yyyyz_0, \
                                         tg_xz_yyyyz_1, tg_xz_yyyz_1, tg_xz_yyyzz_0, tg_xz_yyyzz_1, tg_xz_yyzz_1, \
                                         tg_xz_yyzzz_0, tg_xz_yyzzz_1, tg_xz_yzzz_1, tg_xz_yzzzz_0, tg_xz_yzzzz_1, \
                                         tg_xz_zzzz_1, tg_xz_zzzzz_0, tg_xz_zzzzz_1, tg_y_xxxxx_0, tg_y_xxxxx_1, \
                                         tg_y_xxxxy_0, tg_y_xxxxy_1, tg_y_xxxxz_0, tg_y_xxxxz_1, tg_y_xxxyy_0, tg_y_xxxyy_1, \
                                         tg_y_xxxyz_0, tg_y_xxxyz_1, tg_y_xxxzz_0, tg_y_xxxzz_1, tg_y_xxyyy_0, tg_y_xxyyy_1, \
                                         tg_y_xxyyz_0, tg_y_xxyyz_1, tg_y_xxyzz_0, tg_y_xxyzz_1, tg_y_xxzzz_0, tg_y_xxzzz_1, \
                                         tg_y_xyyyy_0, tg_y_xyyyy_1, tg_y_xyyyz_0, tg_y_xyyyz_1, tg_y_xyyzz_0, tg_y_xyyzz_1, \
                                         tg_y_xyzzz_0, tg_y_xyzzz_1, tg_y_xzzzz_0, tg_y_xzzzz_1, tg_y_yyyyy_0, tg_y_yyyyy_1, \
                                         tg_y_yyyyz_0, tg_y_yyyyz_1, tg_y_yyyzz_0, tg_y_yyyzz_1, tg_y_yyzzz_0, tg_y_yyzzz_1, \
                                         tg_y_yzzzz_0, tg_y_yzzzz_1, tg_y_zzzzz_0, tg_y_zzzzz_1, tg_yy_xxxx_1, \
                                         tg_yy_xxxxx_0, tg_yy_xxxxx_1, tg_yy_xxxxy_0, tg_yy_xxxxy_1, tg_yy_xxxxz_0, \
                                         tg_yy_xxxxz_1, tg_yy_xxxy_1, tg_yy_xxxyy_0, tg_yy_xxxyy_1, tg_yy_xxxyz_0, \
                                         tg_yy_xxxyz_1, tg_yy_xxxz_1, tg_yy_xxxzz_0, tg_yy_xxxzz_1, tg_yy_xxyy_1, \
                                         tg_yy_xxyyy_0, tg_yy_xxyyy_1, tg_yy_xxyz_1, tg_yy_xxzz_1, tg_yy_xyyy_1, \
                                         tg_z_xxxxx_0, tg_z_xxxxx_1, tg_z_xxxxy_0, tg_z_xxxxy_1, tg_z_xxxxz_0, tg_z_xxxxz_1, \
                                         tg_z_xxxyy_0, tg_z_xxxyy_1, tg_z_xxxyz_0, tg_z_xxxyz_1, tg_z_xxxzz_0, tg_z_xxxzz_1, \
                                         tg_z_xxyyy_0, tg_z_xxyyy_1, tg_z_xxyyz_0, tg_z_xxyyz_1, tg_z_xxyzz_0, tg_z_xxyzz_1, \
                                         tg_z_xxzzz_0, tg_z_xxzzz_1, tg_z_xyyyy_0, tg_z_xyyyy_1, tg_z_xyyyz_0, tg_z_xyyyz_1, \
                                         tg_z_xyyzz_0, tg_z_xyyzz_1, tg_z_xyzzz_0, tg_z_xyzzz_1, tg_z_xzzzz_0, tg_z_xzzzz_1, \
                                         tg_z_yyyyy_0, tg_z_yyyyy_1, tg_z_yyyyz_0, tg_z_yyyyz_1, tg_z_yyyzz_0, tg_z_yyyzz_1, \
                                         tg_z_yyzzz_0, tg_z_yyzzz_1, tg_z_yzzzz_0, tg_z_yzzzz_1, tg_z_zzzzz_0, tg_z_zzzzz_1, \
                                         wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxx_xxxxx_0[j] = pb_x * tg_xx_xxxxx_0[j] + wp_x[j] * tg_xx_xxxxx_1[j] + fl1_fx * tg_x_xxxxx_0[j] - fl1_fx * fl1_fza * tg_x_xxxxx_1[j] + 2.5 * fl1_fxn * tg_xx_xxxx_1[j];

                    tg_xxx_xxxxy_0[j] = pb_x * tg_xx_xxxxy_0[j] + wp_x[j] * tg_xx_xxxxy_1[j] + fl1_fx * tg_x_xxxxy_0[j] - fl1_fx * fl1_fza * tg_x_xxxxy_1[j] + 2.0 * fl1_fxn * tg_xx_xxxy_1[j];

                    tg_xxx_xxxxz_0[j] = pb_x * tg_xx_xxxxz_0[j] + wp_x[j] * tg_xx_xxxxz_1[j] + fl1_fx * tg_x_xxxxz_0[j] - fl1_fx * fl1_fza * tg_x_xxxxz_1[j] + 2.0 * fl1_fxn * tg_xx_xxxz_1[j];

                    tg_xxx_xxxyy_0[j] = pb_x * tg_xx_xxxyy_0[j] + wp_x[j] * tg_xx_xxxyy_1[j] + fl1_fx * tg_x_xxxyy_0[j] - fl1_fx * fl1_fza * tg_x_xxxyy_1[j] + 1.5 * fl1_fxn * tg_xx_xxyy_1[j];

                    tg_xxx_xxxyz_0[j] = pb_x * tg_xx_xxxyz_0[j] + wp_x[j] * tg_xx_xxxyz_1[j] + fl1_fx * tg_x_xxxyz_0[j] - fl1_fx * fl1_fza * tg_x_xxxyz_1[j] + 1.5 * fl1_fxn * tg_xx_xxyz_1[j];

                    tg_xxx_xxxzz_0[j] = pb_x * tg_xx_xxxzz_0[j] + wp_x[j] * tg_xx_xxxzz_1[j] + fl1_fx * tg_x_xxxzz_0[j] - fl1_fx * fl1_fza * tg_x_xxxzz_1[j] + 1.5 * fl1_fxn * tg_xx_xxzz_1[j];

                    tg_xxx_xxyyy_0[j] = pb_x * tg_xx_xxyyy_0[j] + wp_x[j] * tg_xx_xxyyy_1[j] + fl1_fx * tg_x_xxyyy_0[j] - fl1_fx * fl1_fza * tg_x_xxyyy_1[j] + fl1_fxn * tg_xx_xyyy_1[j];

                    tg_xxx_xxyyz_0[j] = pb_x * tg_xx_xxyyz_0[j] + wp_x[j] * tg_xx_xxyyz_1[j] + fl1_fx * tg_x_xxyyz_0[j] - fl1_fx * fl1_fza * tg_x_xxyyz_1[j] + fl1_fxn * tg_xx_xyyz_1[j];

                    tg_xxx_xxyzz_0[j] = pb_x * tg_xx_xxyzz_0[j] + wp_x[j] * tg_xx_xxyzz_1[j] + fl1_fx * tg_x_xxyzz_0[j] - fl1_fx * fl1_fza * tg_x_xxyzz_1[j] + fl1_fxn * tg_xx_xyzz_1[j];

                    tg_xxx_xxzzz_0[j] = pb_x * tg_xx_xxzzz_0[j] + wp_x[j] * tg_xx_xxzzz_1[j] + fl1_fx * tg_x_xxzzz_0[j] - fl1_fx * fl1_fza * tg_x_xxzzz_1[j] + fl1_fxn * tg_xx_xzzz_1[j];

                    tg_xxx_xyyyy_0[j] = pb_x * tg_xx_xyyyy_0[j] + wp_x[j] * tg_xx_xyyyy_1[j] + fl1_fx * tg_x_xyyyy_0[j] - fl1_fx * fl1_fza * tg_x_xyyyy_1[j] + 0.5 * fl1_fxn * tg_xx_yyyy_1[j];

                    tg_xxx_xyyyz_0[j] = pb_x * tg_xx_xyyyz_0[j] + wp_x[j] * tg_xx_xyyyz_1[j] + fl1_fx * tg_x_xyyyz_0[j] - fl1_fx * fl1_fza * tg_x_xyyyz_1[j] + 0.5 * fl1_fxn * tg_xx_yyyz_1[j];

                    tg_xxx_xyyzz_0[j] = pb_x * tg_xx_xyyzz_0[j] + wp_x[j] * tg_xx_xyyzz_1[j] + fl1_fx * tg_x_xyyzz_0[j] - fl1_fx * fl1_fza * tg_x_xyyzz_1[j] + 0.5 * fl1_fxn * tg_xx_yyzz_1[j];

                    tg_xxx_xyzzz_0[j] = pb_x * tg_xx_xyzzz_0[j] + wp_x[j] * tg_xx_xyzzz_1[j] + fl1_fx * tg_x_xyzzz_0[j] - fl1_fx * fl1_fza * tg_x_xyzzz_1[j] + 0.5 * fl1_fxn * tg_xx_yzzz_1[j];

                    tg_xxx_xzzzz_0[j] = pb_x * tg_xx_xzzzz_0[j] + wp_x[j] * tg_xx_xzzzz_1[j] + fl1_fx * tg_x_xzzzz_0[j] - fl1_fx * fl1_fza * tg_x_xzzzz_1[j] + 0.5 * fl1_fxn * tg_xx_zzzz_1[j];

                    tg_xxx_yyyyy_0[j] = pb_x * tg_xx_yyyyy_0[j] + wp_x[j] * tg_xx_yyyyy_1[j] + fl1_fx * tg_x_yyyyy_0[j] - fl1_fx * fl1_fza * tg_x_yyyyy_1[j];

                    tg_xxx_yyyyz_0[j] = pb_x * tg_xx_yyyyz_0[j] + wp_x[j] * tg_xx_yyyyz_1[j] + fl1_fx * tg_x_yyyyz_0[j] - fl1_fx * fl1_fza * tg_x_yyyyz_1[j];

                    tg_xxx_yyyzz_0[j] = pb_x * tg_xx_yyyzz_0[j] + wp_x[j] * tg_xx_yyyzz_1[j] + fl1_fx * tg_x_yyyzz_0[j] - fl1_fx * fl1_fza * tg_x_yyyzz_1[j];

                    tg_xxx_yyzzz_0[j] = pb_x * tg_xx_yyzzz_0[j] + wp_x[j] * tg_xx_yyzzz_1[j] + fl1_fx * tg_x_yyzzz_0[j] - fl1_fx * fl1_fza * tg_x_yyzzz_1[j];

                    tg_xxx_yzzzz_0[j] = pb_x * tg_xx_yzzzz_0[j] + wp_x[j] * tg_xx_yzzzz_1[j] + fl1_fx * tg_x_yzzzz_0[j] - fl1_fx * fl1_fza * tg_x_yzzzz_1[j];

                    tg_xxx_zzzzz_0[j] = pb_x * tg_xx_zzzzz_0[j] + wp_x[j] * tg_xx_zzzzz_1[j] + fl1_fx * tg_x_zzzzz_0[j] - fl1_fx * fl1_fza * tg_x_zzzzz_1[j];

                    tg_xxy_xxxxx_0[j] = pb_x * tg_xy_xxxxx_0[j] + wp_x[j] * tg_xy_xxxxx_1[j] + 0.5 * fl1_fx * tg_y_xxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxx_1[j] + 2.5 * fl1_fxn * tg_xy_xxxx_1[j];

                    tg_xxy_xxxxy_0[j] = pb_x * tg_xy_xxxxy_0[j] + wp_x[j] * tg_xy_xxxxy_1[j] + 0.5 * fl1_fx * tg_y_xxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxy_1[j] + 2.0 * fl1_fxn * tg_xy_xxxy_1[j];

                    tg_xxy_xxxxz_0[j] = pb_x * tg_xy_xxxxz_0[j] + wp_x[j] * tg_xy_xxxxz_1[j] + 0.5 * fl1_fx * tg_y_xxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxxz_1[j] + 2.0 * fl1_fxn * tg_xy_xxxz_1[j];

                    tg_xxy_xxxyy_0[j] = pb_x * tg_xy_xxxyy_0[j] + wp_x[j] * tg_xy_xxxyy_1[j] + 0.5 * fl1_fx * tg_y_xxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyy_1[j] + 1.5 * fl1_fxn * tg_xy_xxyy_1[j];

                    tg_xxy_xxxyz_0[j] = pb_x * tg_xy_xxxyz_0[j] + wp_x[j] * tg_xy_xxxyz_1[j] + 0.5 * fl1_fx * tg_y_xxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxyz_1[j] + 1.5 * fl1_fxn * tg_xy_xxyz_1[j];

                    tg_xxy_xxxzz_0[j] = pb_x * tg_xy_xxxzz_0[j] + wp_x[j] * tg_xy_xxxzz_1[j] + 0.5 * fl1_fx * tg_y_xxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxxzz_1[j] + 1.5 * fl1_fxn * tg_xy_xxzz_1[j];

                    tg_xxy_xxyyy_0[j] = pb_x * tg_xy_xxyyy_0[j] + wp_x[j] * tg_xy_xxyyy_1[j] + 0.5 * fl1_fx * tg_y_xxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyy_1[j] + fl1_fxn * tg_xy_xyyy_1[j];

                    tg_xxy_xxyyz_0[j] = pb_x * tg_xy_xxyyz_0[j] + wp_x[j] * tg_xy_xxyyz_1[j] + 0.5 * fl1_fx * tg_y_xxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyyz_1[j] + fl1_fxn * tg_xy_xyyz_1[j];

                    tg_xxy_xxyzz_0[j] = pb_x * tg_xy_xxyzz_0[j] + wp_x[j] * tg_xy_xxyzz_1[j] + 0.5 * fl1_fx * tg_y_xxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxyzz_1[j] + fl1_fxn * tg_xy_xyzz_1[j];

                    tg_xxy_xxzzz_0[j] = pb_x * tg_xy_xxzzz_0[j] + wp_x[j] * tg_xy_xxzzz_1[j] + 0.5 * fl1_fx * tg_y_xxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xxzzz_1[j] + fl1_fxn * tg_xy_xzzz_1[j];

                    tg_xxy_xyyyy_0[j] = pb_x * tg_xy_xyyyy_0[j] + wp_x[j] * tg_xy_xyyyy_1[j] + 0.5 * fl1_fx * tg_y_xyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyy_1[j] + 0.5 * fl1_fxn * tg_xy_yyyy_1[j];

                    tg_xxy_xyyyz_0[j] = pb_x * tg_xy_xyyyz_0[j] + wp_x[j] * tg_xy_xyyyz_1[j] + 0.5 * fl1_fx * tg_y_xyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyyz_1[j] + 0.5 * fl1_fxn * tg_xy_yyyz_1[j];

                    tg_xxy_xyyzz_0[j] = pb_x * tg_xy_xyyzz_0[j] + wp_x[j] * tg_xy_xyyzz_1[j] + 0.5 * fl1_fx * tg_y_xyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyyzz_1[j] + 0.5 * fl1_fxn * tg_xy_yyzz_1[j];

                    tg_xxy_xyzzz_0[j] = pb_x * tg_xy_xyzzz_0[j] + wp_x[j] * tg_xy_xyzzz_1[j] + 0.5 * fl1_fx * tg_y_xyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xyzzz_1[j] + 0.5 * fl1_fxn * tg_xy_yzzz_1[j];

                    tg_xxy_xzzzz_0[j] = pb_x * tg_xy_xzzzz_0[j] + wp_x[j] * tg_xy_xzzzz_1[j] + 0.5 * fl1_fx * tg_y_xzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_xzzzz_1[j] + 0.5 * fl1_fxn * tg_xy_zzzz_1[j];

                    tg_xxy_yyyyy_0[j] = pb_x * tg_xy_yyyyy_0[j] + wp_x[j] * tg_xy_yyyyy_1[j] + 0.5 * fl1_fx * tg_y_yyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyy_1[j];

                    tg_xxy_yyyyz_0[j] = pb_x * tg_xy_yyyyz_0[j] + wp_x[j] * tg_xy_yyyyz_1[j] + 0.5 * fl1_fx * tg_y_yyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyyz_1[j];

                    tg_xxy_yyyzz_0[j] = pb_x * tg_xy_yyyzz_0[j] + wp_x[j] * tg_xy_yyyzz_1[j] + 0.5 * fl1_fx * tg_y_yyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyyzz_1[j];

                    tg_xxy_yyzzz_0[j] = pb_x * tg_xy_yyzzz_0[j] + wp_x[j] * tg_xy_yyzzz_1[j] + 0.5 * fl1_fx * tg_y_yyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yyzzz_1[j];

                    tg_xxy_yzzzz_0[j] = pb_x * tg_xy_yzzzz_0[j] + wp_x[j] * tg_xy_yzzzz_1[j] + 0.5 * fl1_fx * tg_y_yzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_yzzzz_1[j];

                    tg_xxy_zzzzz_0[j] = pb_x * tg_xy_zzzzz_0[j] + wp_x[j] * tg_xy_zzzzz_1[j] + 0.5 * fl1_fx * tg_y_zzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_zzzzz_1[j];

                    tg_xxz_xxxxx_0[j] = pb_x * tg_xz_xxxxx_0[j] + wp_x[j] * tg_xz_xxxxx_1[j] + 0.5 * fl1_fx * tg_z_xxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxx_1[j] + 2.5 * fl1_fxn * tg_xz_xxxx_1[j];

                    tg_xxz_xxxxy_0[j] = pb_x * tg_xz_xxxxy_0[j] + wp_x[j] * tg_xz_xxxxy_1[j] + 0.5 * fl1_fx * tg_z_xxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxy_1[j] + 2.0 * fl1_fxn * tg_xz_xxxy_1[j];

                    tg_xxz_xxxxz_0[j] = pb_x * tg_xz_xxxxz_0[j] + wp_x[j] * tg_xz_xxxxz_1[j] + 0.5 * fl1_fx * tg_z_xxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxz_1[j] + 2.0 * fl1_fxn * tg_xz_xxxz_1[j];

                    tg_xxz_xxxyy_0[j] = pb_x * tg_xz_xxxyy_0[j] + wp_x[j] * tg_xz_xxxyy_1[j] + 0.5 * fl1_fx * tg_z_xxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyy_1[j] + 1.5 * fl1_fxn * tg_xz_xxyy_1[j];

                    tg_xxz_xxxyz_0[j] = pb_x * tg_xz_xxxyz_0[j] + wp_x[j] * tg_xz_xxxyz_1[j] + 0.5 * fl1_fx * tg_z_xxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyz_1[j] + 1.5 * fl1_fxn * tg_xz_xxyz_1[j];

                    tg_xxz_xxxzz_0[j] = pb_x * tg_xz_xxxzz_0[j] + wp_x[j] * tg_xz_xxxzz_1[j] + 0.5 * fl1_fx * tg_z_xxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxzz_1[j] + 1.5 * fl1_fxn * tg_xz_xxzz_1[j];

                    tg_xxz_xxyyy_0[j] = pb_x * tg_xz_xxyyy_0[j] + wp_x[j] * tg_xz_xxyyy_1[j] + 0.5 * fl1_fx * tg_z_xxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyy_1[j] + fl1_fxn * tg_xz_xyyy_1[j];

                    tg_xxz_xxyyz_0[j] = pb_x * tg_xz_xxyyz_0[j] + wp_x[j] * tg_xz_xxyyz_1[j] + 0.5 * fl1_fx * tg_z_xxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyz_1[j] + fl1_fxn * tg_xz_xyyz_1[j];

                    tg_xxz_xxyzz_0[j] = pb_x * tg_xz_xxyzz_0[j] + wp_x[j] * tg_xz_xxyzz_1[j] + 0.5 * fl1_fx * tg_z_xxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyzz_1[j] + fl1_fxn * tg_xz_xyzz_1[j];

                    tg_xxz_xxzzz_0[j] = pb_x * tg_xz_xxzzz_0[j] + wp_x[j] * tg_xz_xxzzz_1[j] + 0.5 * fl1_fx * tg_z_xxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxzzz_1[j] + fl1_fxn * tg_xz_xzzz_1[j];

                    tg_xxz_xyyyy_0[j] = pb_x * tg_xz_xyyyy_0[j] + wp_x[j] * tg_xz_xyyyy_1[j] + 0.5 * fl1_fx * tg_z_xyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyy_1[j] + 0.5 * fl1_fxn * tg_xz_yyyy_1[j];

                    tg_xxz_xyyyz_0[j] = pb_x * tg_xz_xyyyz_0[j] + wp_x[j] * tg_xz_xyyyz_1[j] + 0.5 * fl1_fx * tg_z_xyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyz_1[j] + 0.5 * fl1_fxn * tg_xz_yyyz_1[j];

                    tg_xxz_xyyzz_0[j] = pb_x * tg_xz_xyyzz_0[j] + wp_x[j] * tg_xz_xyyzz_1[j] + 0.5 * fl1_fx * tg_z_xyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyzz_1[j] + 0.5 * fl1_fxn * tg_xz_yyzz_1[j];

                    tg_xxz_xyzzz_0[j] = pb_x * tg_xz_xyzzz_0[j] + wp_x[j] * tg_xz_xyzzz_1[j] + 0.5 * fl1_fx * tg_z_xyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyzzz_1[j] + 0.5 * fl1_fxn * tg_xz_yzzz_1[j];

                    tg_xxz_xzzzz_0[j] = pb_x * tg_xz_xzzzz_0[j] + wp_x[j] * tg_xz_xzzzz_1[j] + 0.5 * fl1_fx * tg_z_xzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xzzzz_1[j] + 0.5 * fl1_fxn * tg_xz_zzzz_1[j];

                    tg_xxz_yyyyy_0[j] = pb_x * tg_xz_yyyyy_0[j] + wp_x[j] * tg_xz_yyyyy_1[j] + 0.5 * fl1_fx * tg_z_yyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyy_1[j];

                    tg_xxz_yyyyz_0[j] = pb_x * tg_xz_yyyyz_0[j] + wp_x[j] * tg_xz_yyyyz_1[j] + 0.5 * fl1_fx * tg_z_yyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyz_1[j];

                    tg_xxz_yyyzz_0[j] = pb_x * tg_xz_yyyzz_0[j] + wp_x[j] * tg_xz_yyyzz_1[j] + 0.5 * fl1_fx * tg_z_yyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyzz_1[j];

                    tg_xxz_yyzzz_0[j] = pb_x * tg_xz_yyzzz_0[j] + wp_x[j] * tg_xz_yyzzz_1[j] + 0.5 * fl1_fx * tg_z_yyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyzzz_1[j];

                    tg_xxz_yzzzz_0[j] = pb_x * tg_xz_yzzzz_0[j] + wp_x[j] * tg_xz_yzzzz_1[j] + 0.5 * fl1_fx * tg_z_yzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yzzzz_1[j];

                    tg_xxz_zzzzz_0[j] = pb_x * tg_xz_zzzzz_0[j] + wp_x[j] * tg_xz_zzzzz_1[j] + 0.5 * fl1_fx * tg_z_zzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_zzzzz_1[j];

                    tg_xyy_xxxxx_0[j] = pb_x * tg_yy_xxxxx_0[j] + wp_x[j] * tg_yy_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yy_xxxx_1[j];

                    tg_xyy_xxxxy_0[j] = pb_x * tg_yy_xxxxy_0[j] + wp_x[j] * tg_yy_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yy_xxxy_1[j];

                    tg_xyy_xxxxz_0[j] = pb_x * tg_yy_xxxxz_0[j] + wp_x[j] * tg_yy_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yy_xxxz_1[j];

                    tg_xyy_xxxyy_0[j] = pb_x * tg_yy_xxxyy_0[j] + wp_x[j] * tg_yy_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yy_xxyy_1[j];

                    tg_xyy_xxxyz_0[j] = pb_x * tg_yy_xxxyz_0[j] + wp_x[j] * tg_yy_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yy_xxyz_1[j];

                    tg_xyy_xxxzz_0[j] = pb_x * tg_yy_xxxzz_0[j] + wp_x[j] * tg_yy_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yy_xxzz_1[j];

                    tg_xyy_xxyyy_0[j] = pb_x * tg_yy_xxyyy_0[j] + wp_x[j] * tg_yy_xxyyy_1[j] + fl1_fxn * tg_yy_xyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSH_70_140(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (70,140)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yy_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 63); 

                auto tg_yy_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 64); 

                auto tg_yy_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 65); 

                auto tg_yy_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 66); 

                auto tg_yy_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 67); 

                auto tg_yy_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 68); 

                auto tg_yy_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 69); 

                auto tg_yy_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 70); 

                auto tg_yy_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 71); 

                auto tg_yy_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 72); 

                auto tg_yy_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 73); 

                auto tg_yy_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 74); 

                auto tg_yy_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 75); 

                auto tg_yy_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 76); 

                auto tg_yy_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 77); 

                auto tg_yy_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 78); 

                auto tg_yy_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 79); 

                auto tg_yy_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 80); 

                auto tg_yy_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 81); 

                auto tg_yy_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 82); 

                auto tg_yy_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 83); 

                auto tg_yz_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 84); 

                auto tg_yz_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 85); 

                auto tg_yz_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 86); 

                auto tg_yz_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 87); 

                auto tg_yz_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 88); 

                auto tg_yz_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 89); 

                auto tg_yz_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 90); 

                auto tg_yz_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 91); 

                auto tg_yz_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 92); 

                auto tg_yz_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 93); 

                auto tg_yz_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 94); 

                auto tg_yz_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 95); 

                auto tg_yz_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 96); 

                auto tg_yz_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 97); 

                auto tg_yz_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 98); 

                auto tg_yz_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 99); 

                auto tg_yz_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 100); 

                auto tg_yz_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 101); 

                auto tg_yz_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 102); 

                auto tg_yz_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 103); 

                auto tg_yz_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 104); 

                auto tg_zz_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 105); 

                auto tg_zz_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 106); 

                auto tg_zz_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 107); 

                auto tg_zz_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 108); 

                auto tg_zz_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 109); 

                auto tg_zz_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 110); 

                auto tg_zz_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 111); 

                auto tg_zz_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 112); 

                auto tg_zz_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 113); 

                auto tg_zz_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 114); 

                auto tg_zz_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 115); 

                auto tg_zz_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 116); 

                auto tg_zz_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 117); 

                auto tg_zz_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 118); 

                auto tg_zz_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 119); 

                auto tg_zz_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 120); 

                auto tg_zz_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 121); 

                auto tg_zz_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 122); 

                auto tg_zz_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 123); 

                auto tg_zz_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 124); 

                auto tg_zz_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 125); 

                auto tg_yy_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 63); 

                auto tg_yy_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 64); 

                auto tg_yy_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 65); 

                auto tg_yy_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 66); 

                auto tg_yy_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 67); 

                auto tg_yy_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 68); 

                auto tg_yy_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 69); 

                auto tg_yy_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 70); 

                auto tg_yy_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 71); 

                auto tg_yy_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 72); 

                auto tg_yy_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 73); 

                auto tg_yy_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 74); 

                auto tg_yy_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 75); 

                auto tg_yy_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 76); 

                auto tg_yy_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 77); 

                auto tg_yy_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 78); 

                auto tg_yy_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 79); 

                auto tg_yy_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 80); 

                auto tg_yy_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 81); 

                auto tg_yy_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 82); 

                auto tg_yy_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 83); 

                auto tg_yz_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 84); 

                auto tg_yz_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 85); 

                auto tg_yz_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 86); 

                auto tg_yz_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 87); 

                auto tg_yz_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 88); 

                auto tg_yz_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 89); 

                auto tg_yz_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 90); 

                auto tg_yz_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 91); 

                auto tg_yz_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 92); 

                auto tg_yz_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 93); 

                auto tg_yz_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 94); 

                auto tg_yz_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 95); 

                auto tg_yz_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 96); 

                auto tg_yz_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 97); 

                auto tg_yz_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 98); 

                auto tg_yz_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 99); 

                auto tg_yz_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 100); 

                auto tg_yz_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 101); 

                auto tg_yz_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 102); 

                auto tg_yz_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 103); 

                auto tg_yz_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 104); 

                auto tg_zz_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 105); 

                auto tg_zz_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 106); 

                auto tg_zz_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 107); 

                auto tg_zz_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 108); 

                auto tg_zz_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 109); 

                auto tg_zz_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 110); 

                auto tg_zz_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 111); 

                auto tg_zz_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 112); 

                auto tg_zz_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 113); 

                auto tg_zz_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 114); 

                auto tg_zz_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 115); 

                auto tg_zz_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 116); 

                auto tg_zz_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 117); 

                auto tg_zz_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 118); 

                auto tg_zz_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 119); 

                auto tg_zz_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 120); 

                auto tg_zz_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 121); 

                auto tg_zz_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 122); 

                auto tg_zz_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 123); 

                auto tg_zz_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 124); 

                auto tg_zz_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 125); 

                auto tg_y_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 21); 

                auto tg_y_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 22); 

                auto tg_y_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 23); 

                auto tg_y_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 24); 

                auto tg_y_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 25); 

                auto tg_y_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 26); 

                auto tg_y_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 27); 

                auto tg_y_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 28); 

                auto tg_y_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 29); 

                auto tg_y_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 30); 

                auto tg_y_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 31); 

                auto tg_y_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 32); 

                auto tg_y_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 33); 

                auto tg_y_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 34); 

                auto tg_y_xxxxx_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 21); 

                auto tg_y_xxxxy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 22); 

                auto tg_y_xxxxz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 23); 

                auto tg_y_xxxyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 24); 

                auto tg_y_xxxyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 25); 

                auto tg_y_xxxzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 26); 

                auto tg_y_xxyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 27); 

                auto tg_y_xxyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 28); 

                auto tg_y_xxyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 29); 

                auto tg_y_xxzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 30); 

                auto tg_y_xyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 31); 

                auto tg_y_xyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 32); 

                auto tg_y_xyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 33); 

                auto tg_y_xyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 34); 

                auto tg_yy_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 45); 

                auto tg_yy_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 46); 

                auto tg_yy_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 47); 

                auto tg_yy_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 48); 

                auto tg_yy_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 49); 

                auto tg_yy_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 50); 

                auto tg_yy_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 51); 

                auto tg_yy_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 52); 

                auto tg_yy_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 53); 

                auto tg_yy_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 54); 

                auto tg_yy_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 55); 

                auto tg_yy_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 56); 

                auto tg_yy_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 57); 

                auto tg_yy_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 58); 

                auto tg_yy_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 59); 

                auto tg_yz_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 60); 

                auto tg_yz_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 61); 

                auto tg_yz_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 62); 

                auto tg_yz_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 63); 

                auto tg_yz_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 64); 

                auto tg_yz_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 65); 

                auto tg_yz_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 66); 

                auto tg_yz_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 67); 

                auto tg_yz_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 68); 

                auto tg_yz_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 69); 

                auto tg_yz_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 70); 

                auto tg_yz_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 71); 

                auto tg_yz_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 72); 

                auto tg_yz_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 73); 

                auto tg_yz_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 74); 

                auto tg_zz_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 75); 

                auto tg_zz_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 76); 

                auto tg_zz_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 77); 

                auto tg_zz_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 78); 

                auto tg_zz_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 79); 

                auto tg_zz_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 80); 

                auto tg_zz_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 81); 

                auto tg_zz_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 82); 

                auto tg_zz_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 83); 

                auto tg_zz_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 84); 

                auto tg_zz_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 85); 

                auto tg_zz_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 86); 

                auto tg_zz_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 87); 

                auto tg_zz_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 88); 

                auto tg_zz_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 89); 

                // set up pointers to integrals

                auto tg_xyy_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 70); 

                auto tg_xyy_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 71); 

                auto tg_xyy_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 72); 

                auto tg_xyy_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 73); 

                auto tg_xyy_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 74); 

                auto tg_xyy_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 75); 

                auto tg_xyy_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 76); 

                auto tg_xyy_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 77); 

                auto tg_xyy_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 78); 

                auto tg_xyy_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 79); 

                auto tg_xyy_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 80); 

                auto tg_xyy_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 81); 

                auto tg_xyy_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 82); 

                auto tg_xyy_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 83); 

                auto tg_xyz_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 84); 

                auto tg_xyz_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 85); 

                auto tg_xyz_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 86); 

                auto tg_xyz_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 87); 

                auto tg_xyz_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 88); 

                auto tg_xyz_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 89); 

                auto tg_xyz_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 90); 

                auto tg_xyz_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 91); 

                auto tg_xyz_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 92); 

                auto tg_xyz_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 93); 

                auto tg_xyz_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 94); 

                auto tg_xyz_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 95); 

                auto tg_xyz_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 96); 

                auto tg_xyz_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 97); 

                auto tg_xyz_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 98); 

                auto tg_xyz_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 99); 

                auto tg_xyz_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 100); 

                auto tg_xyz_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 101); 

                auto tg_xyz_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 102); 

                auto tg_xyz_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 103); 

                auto tg_xyz_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 104); 

                auto tg_xzz_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 105); 

                auto tg_xzz_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 106); 

                auto tg_xzz_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 107); 

                auto tg_xzz_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 108); 

                auto tg_xzz_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 109); 

                auto tg_xzz_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 110); 

                auto tg_xzz_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 111); 

                auto tg_xzz_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 112); 

                auto tg_xzz_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 113); 

                auto tg_xzz_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 114); 

                auto tg_xzz_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 115); 

                auto tg_xzz_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 116); 

                auto tg_xzz_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 117); 

                auto tg_xzz_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 118); 

                auto tg_xzz_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 119); 

                auto tg_xzz_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 120); 

                auto tg_xzz_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 121); 

                auto tg_xzz_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 122); 

                auto tg_xzz_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 123); 

                auto tg_xzz_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 124); 

                auto tg_xzz_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 125); 

                auto tg_yyy_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 126); 

                auto tg_yyy_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 127); 

                auto tg_yyy_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 128); 

                auto tg_yyy_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 129); 

                auto tg_yyy_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 130); 

                auto tg_yyy_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 131); 

                auto tg_yyy_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 132); 

                auto tg_yyy_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 133); 

                auto tg_yyy_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 134); 

                auto tg_yyy_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 135); 

                auto tg_yyy_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 136); 

                auto tg_yyy_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 137); 

                auto tg_yyy_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 138); 

                auto tg_yyy_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 139); 

                // Batch of Integrals (70,140)

                #pragma omp simd aligned(fxn, fza, tg_xyy_xxyyz_0, tg_xyy_xxyzz_0, tg_xyy_xxzzz_0, \
                                         tg_xyy_xyyyy_0, tg_xyy_xyyyz_0, tg_xyy_xyyzz_0, tg_xyy_xyzzz_0, tg_xyy_xzzzz_0, \
                                         tg_xyy_yyyyy_0, tg_xyy_yyyyz_0, tg_xyy_yyyzz_0, tg_xyy_yyzzz_0, tg_xyy_yzzzz_0, \
                                         tg_xyy_zzzzz_0, tg_xyz_xxxxx_0, tg_xyz_xxxxy_0, tg_xyz_xxxxz_0, tg_xyz_xxxyy_0, \
                                         tg_xyz_xxxyz_0, tg_xyz_xxxzz_0, tg_xyz_xxyyy_0, tg_xyz_xxyyz_0, tg_xyz_xxyzz_0, \
                                         tg_xyz_xxzzz_0, tg_xyz_xyyyy_0, tg_xyz_xyyyz_0, tg_xyz_xyyzz_0, tg_xyz_xyzzz_0, \
                                         tg_xyz_xzzzz_0, tg_xyz_yyyyy_0, tg_xyz_yyyyz_0, tg_xyz_yyyzz_0, tg_xyz_yyzzz_0, \
                                         tg_xyz_yzzzz_0, tg_xyz_zzzzz_0, tg_xzz_xxxxx_0, tg_xzz_xxxxy_0, tg_xzz_xxxxz_0, \
                                         tg_xzz_xxxyy_0, tg_xzz_xxxyz_0, tg_xzz_xxxzz_0, tg_xzz_xxyyy_0, tg_xzz_xxyyz_0, \
                                         tg_xzz_xxyzz_0, tg_xzz_xxzzz_0, tg_xzz_xyyyy_0, tg_xzz_xyyyz_0, tg_xzz_xyyzz_0, \
                                         tg_xzz_xyzzz_0, tg_xzz_xzzzz_0, tg_xzz_yyyyy_0, tg_xzz_yyyyz_0, tg_xzz_yyyzz_0, \
                                         tg_xzz_yyzzz_0, tg_xzz_yzzzz_0, tg_xzz_zzzzz_0, tg_y_xxxxx_0, tg_y_xxxxx_1, \
                                         tg_y_xxxxy_0, tg_y_xxxxy_1, tg_y_xxxxz_0, tg_y_xxxxz_1, tg_y_xxxyy_0, tg_y_xxxyy_1, \
                                         tg_y_xxxyz_0, tg_y_xxxyz_1, tg_y_xxxzz_0, tg_y_xxxzz_1, tg_y_xxyyy_0, tg_y_xxyyy_1, \
                                         tg_y_xxyyz_0, tg_y_xxyyz_1, tg_y_xxyzz_0, tg_y_xxyzz_1, tg_y_xxzzz_0, tg_y_xxzzz_1, \
                                         tg_y_xyyyy_0, tg_y_xyyyy_1, tg_y_xyyyz_0, tg_y_xyyyz_1, tg_y_xyyzz_0, tg_y_xyyzz_1, \
                                         tg_y_xyzzz_0, tg_y_xyzzz_1, tg_yy_xxxx_1, tg_yy_xxxxx_0, tg_yy_xxxxx_1, \
                                         tg_yy_xxxxy_0, tg_yy_xxxxy_1, tg_yy_xxxxz_0, tg_yy_xxxxz_1, tg_yy_xxxy_1, \
                                         tg_yy_xxxyy_0, tg_yy_xxxyy_1, tg_yy_xxxyz_0, tg_yy_xxxyz_1, tg_yy_xxxz_1, \
                                         tg_yy_xxxzz_0, tg_yy_xxxzz_1, tg_yy_xxyy_1, tg_yy_xxyyy_0, tg_yy_xxyyy_1, \
                                         tg_yy_xxyyz_0, tg_yy_xxyyz_1, tg_yy_xxyz_1, tg_yy_xxyzz_0, tg_yy_xxyzz_1, \
                                         tg_yy_xxzz_1, tg_yy_xxzzz_0, tg_yy_xxzzz_1, tg_yy_xyyy_1, tg_yy_xyyyy_0, \
                                         tg_yy_xyyyy_1, tg_yy_xyyyz_0, tg_yy_xyyyz_1, tg_yy_xyyz_1, tg_yy_xyyzz_0, \
                                         tg_yy_xyyzz_1, tg_yy_xyzz_1, tg_yy_xyzzz_0, tg_yy_xyzzz_1, tg_yy_xzzz_1, \
                                         tg_yy_xzzzz_0, tg_yy_xzzzz_1, tg_yy_yyyy_1, tg_yy_yyyyy_0, tg_yy_yyyyy_1, \
                                         tg_yy_yyyyz_0, tg_yy_yyyyz_1, tg_yy_yyyz_1, tg_yy_yyyzz_0, tg_yy_yyyzz_1, \
                                         tg_yy_yyzz_1, tg_yy_yyzzz_0, tg_yy_yyzzz_1, tg_yy_yzzz_1, tg_yy_yzzzz_0, \
                                         tg_yy_yzzzz_1, tg_yy_zzzz_1, tg_yy_zzzzz_0, tg_yy_zzzzz_1, tg_yyy_xxxxx_0, \
                                         tg_yyy_xxxxy_0, tg_yyy_xxxxz_0, tg_yyy_xxxyy_0, tg_yyy_xxxyz_0, tg_yyy_xxxzz_0, \
                                         tg_yyy_xxyyy_0, tg_yyy_xxyyz_0, tg_yyy_xxyzz_0, tg_yyy_xxzzz_0, tg_yyy_xyyyy_0, \
                                         tg_yyy_xyyyz_0, tg_yyy_xyyzz_0, tg_yyy_xyzzz_0, tg_yz_xxxx_1, tg_yz_xxxxx_0, \
                                         tg_yz_xxxxx_1, tg_yz_xxxxy_0, tg_yz_xxxxy_1, tg_yz_xxxxz_0, tg_yz_xxxxz_1, \
                                         tg_yz_xxxy_1, tg_yz_xxxyy_0, tg_yz_xxxyy_1, tg_yz_xxxyz_0, tg_yz_xxxyz_1, \
                                         tg_yz_xxxz_1, tg_yz_xxxzz_0, tg_yz_xxxzz_1, tg_yz_xxyy_1, tg_yz_xxyyy_0, \
                                         tg_yz_xxyyy_1, tg_yz_xxyyz_0, tg_yz_xxyyz_1, tg_yz_xxyz_1, tg_yz_xxyzz_0, \
                                         tg_yz_xxyzz_1, tg_yz_xxzz_1, tg_yz_xxzzz_0, tg_yz_xxzzz_1, tg_yz_xyyy_1, \
                                         tg_yz_xyyyy_0, tg_yz_xyyyy_1, tg_yz_xyyyz_0, tg_yz_xyyyz_1, tg_yz_xyyz_1, \
                                         tg_yz_xyyzz_0, tg_yz_xyyzz_1, tg_yz_xyzz_1, tg_yz_xyzzz_0, tg_yz_xyzzz_1, \
                                         tg_yz_xzzz_1, tg_yz_xzzzz_0, tg_yz_xzzzz_1, tg_yz_yyyy_1, tg_yz_yyyyy_0, \
                                         tg_yz_yyyyy_1, tg_yz_yyyyz_0, tg_yz_yyyyz_1, tg_yz_yyyz_1, tg_yz_yyyzz_0, \
                                         tg_yz_yyyzz_1, tg_yz_yyzz_1, tg_yz_yyzzz_0, tg_yz_yyzzz_1, tg_yz_yzzz_1, \
                                         tg_yz_yzzzz_0, tg_yz_yzzzz_1, tg_yz_zzzz_1, tg_yz_zzzzz_0, tg_yz_zzzzz_1, \
                                         tg_zz_xxxx_1, tg_zz_xxxxx_0, tg_zz_xxxxx_1, tg_zz_xxxxy_0, tg_zz_xxxxy_1, \
                                         tg_zz_xxxxz_0, tg_zz_xxxxz_1, tg_zz_xxxy_1, tg_zz_xxxyy_0, tg_zz_xxxyy_1, \
                                         tg_zz_xxxyz_0, tg_zz_xxxyz_1, tg_zz_xxxz_1, tg_zz_xxxzz_0, tg_zz_xxxzz_1, \
                                         tg_zz_xxyy_1, tg_zz_xxyyy_0, tg_zz_xxyyy_1, tg_zz_xxyyz_0, tg_zz_xxyyz_1, \
                                         tg_zz_xxyz_1, tg_zz_xxyzz_0, tg_zz_xxyzz_1, tg_zz_xxzz_1, tg_zz_xxzzz_0, \
                                         tg_zz_xxzzz_1, tg_zz_xyyy_1, tg_zz_xyyyy_0, tg_zz_xyyyy_1, tg_zz_xyyyz_0, \
                                         tg_zz_xyyyz_1, tg_zz_xyyz_1, tg_zz_xyyzz_0, tg_zz_xyyzz_1, tg_zz_xyzz_1, \
                                         tg_zz_xyzzz_0, tg_zz_xyzzz_1, tg_zz_xzzz_1, tg_zz_xzzzz_0, tg_zz_xzzzz_1, \
                                         tg_zz_yyyy_1, tg_zz_yyyyy_0, tg_zz_yyyyy_1, tg_zz_yyyyz_0, tg_zz_yyyyz_1, \
                                         tg_zz_yyyz_1, tg_zz_yyyzz_0, tg_zz_yyyzz_1, tg_zz_yyzz_1, tg_zz_yyzzz_0, \
                                         tg_zz_yyzzz_1, tg_zz_yzzz_1, tg_zz_yzzzz_0, tg_zz_yzzzz_1, tg_zz_zzzz_1, \
                                         tg_zz_zzzzz_0, tg_zz_zzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xyy_xxyyz_0[j] = pb_x * tg_yy_xxyyz_0[j] + wp_x[j] * tg_yy_xxyyz_1[j] + fl1_fxn * tg_yy_xyyz_1[j];

                    tg_xyy_xxyzz_0[j] = pb_x * tg_yy_xxyzz_0[j] + wp_x[j] * tg_yy_xxyzz_1[j] + fl1_fxn * tg_yy_xyzz_1[j];

                    tg_xyy_xxzzz_0[j] = pb_x * tg_yy_xxzzz_0[j] + wp_x[j] * tg_yy_xxzzz_1[j] + fl1_fxn * tg_yy_xzzz_1[j];

                    tg_xyy_xyyyy_0[j] = pb_x * tg_yy_xyyyy_0[j] + wp_x[j] * tg_yy_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yy_yyyy_1[j];

                    tg_xyy_xyyyz_0[j] = pb_x * tg_yy_xyyyz_0[j] + wp_x[j] * tg_yy_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yy_yyyz_1[j];

                    tg_xyy_xyyzz_0[j] = pb_x * tg_yy_xyyzz_0[j] + wp_x[j] * tg_yy_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yy_yyzz_1[j];

                    tg_xyy_xyzzz_0[j] = pb_x * tg_yy_xyzzz_0[j] + wp_x[j] * tg_yy_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yy_yzzz_1[j];

                    tg_xyy_xzzzz_0[j] = pb_x * tg_yy_xzzzz_0[j] + wp_x[j] * tg_yy_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzzz_1[j];

                    tg_xyy_yyyyy_0[j] = pb_x * tg_yy_yyyyy_0[j] + wp_x[j] * tg_yy_yyyyy_1[j];

                    tg_xyy_yyyyz_0[j] = pb_x * tg_yy_yyyyz_0[j] + wp_x[j] * tg_yy_yyyyz_1[j];

                    tg_xyy_yyyzz_0[j] = pb_x * tg_yy_yyyzz_0[j] + wp_x[j] * tg_yy_yyyzz_1[j];

                    tg_xyy_yyzzz_0[j] = pb_x * tg_yy_yyzzz_0[j] + wp_x[j] * tg_yy_yyzzz_1[j];

                    tg_xyy_yzzzz_0[j] = pb_x * tg_yy_yzzzz_0[j] + wp_x[j] * tg_yy_yzzzz_1[j];

                    tg_xyy_zzzzz_0[j] = pb_x * tg_yy_zzzzz_0[j] + wp_x[j] * tg_yy_zzzzz_1[j];

                    tg_xyz_xxxxx_0[j] = pb_x * tg_yz_xxxxx_0[j] + wp_x[j] * tg_yz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yz_xxxx_1[j];

                    tg_xyz_xxxxy_0[j] = pb_x * tg_yz_xxxxy_0[j] + wp_x[j] * tg_yz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yz_xxxy_1[j];

                    tg_xyz_xxxxz_0[j] = pb_x * tg_yz_xxxxz_0[j] + wp_x[j] * tg_yz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yz_xxxz_1[j];

                    tg_xyz_xxxyy_0[j] = pb_x * tg_yz_xxxyy_0[j] + wp_x[j] * tg_yz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yz_xxyy_1[j];

                    tg_xyz_xxxyz_0[j] = pb_x * tg_yz_xxxyz_0[j] + wp_x[j] * tg_yz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yz_xxyz_1[j];

                    tg_xyz_xxxzz_0[j] = pb_x * tg_yz_xxxzz_0[j] + wp_x[j] * tg_yz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yz_xxzz_1[j];

                    tg_xyz_xxyyy_0[j] = pb_x * tg_yz_xxyyy_0[j] + wp_x[j] * tg_yz_xxyyy_1[j] + fl1_fxn * tg_yz_xyyy_1[j];

                    tg_xyz_xxyyz_0[j] = pb_x * tg_yz_xxyyz_0[j] + wp_x[j] * tg_yz_xxyyz_1[j] + fl1_fxn * tg_yz_xyyz_1[j];

                    tg_xyz_xxyzz_0[j] = pb_x * tg_yz_xxyzz_0[j] + wp_x[j] * tg_yz_xxyzz_1[j] + fl1_fxn * tg_yz_xyzz_1[j];

                    tg_xyz_xxzzz_0[j] = pb_x * tg_yz_xxzzz_0[j] + wp_x[j] * tg_yz_xxzzz_1[j] + fl1_fxn * tg_yz_xzzz_1[j];

                    tg_xyz_xyyyy_0[j] = pb_x * tg_yz_xyyyy_0[j] + wp_x[j] * tg_yz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yz_yyyy_1[j];

                    tg_xyz_xyyyz_0[j] = pb_x * tg_yz_xyyyz_0[j] + wp_x[j] * tg_yz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yz_yyyz_1[j];

                    tg_xyz_xyyzz_0[j] = pb_x * tg_yz_xyyzz_0[j] + wp_x[j] * tg_yz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yz_yyzz_1[j];

                    tg_xyz_xyzzz_0[j] = pb_x * tg_yz_xyzzz_0[j] + wp_x[j] * tg_yz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yz_yzzz_1[j];

                    tg_xyz_xzzzz_0[j] = pb_x * tg_yz_xzzzz_0[j] + wp_x[j] * tg_yz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzzz_1[j];

                    tg_xyz_yyyyy_0[j] = pb_x * tg_yz_yyyyy_0[j] + wp_x[j] * tg_yz_yyyyy_1[j];

                    tg_xyz_yyyyz_0[j] = pb_x * tg_yz_yyyyz_0[j] + wp_x[j] * tg_yz_yyyyz_1[j];

                    tg_xyz_yyyzz_0[j] = pb_x * tg_yz_yyyzz_0[j] + wp_x[j] * tg_yz_yyyzz_1[j];

                    tg_xyz_yyzzz_0[j] = pb_x * tg_yz_yyzzz_0[j] + wp_x[j] * tg_yz_yyzzz_1[j];

                    tg_xyz_yzzzz_0[j] = pb_x * tg_yz_yzzzz_0[j] + wp_x[j] * tg_yz_yzzzz_1[j];

                    tg_xyz_zzzzz_0[j] = pb_x * tg_yz_zzzzz_0[j] + wp_x[j] * tg_yz_zzzzz_1[j];

                    tg_xzz_xxxxx_0[j] = pb_x * tg_zz_xxxxx_0[j] + wp_x[j] * tg_zz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_zz_xxxx_1[j];

                    tg_xzz_xxxxy_0[j] = pb_x * tg_zz_xxxxy_0[j] + wp_x[j] * tg_zz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_zz_xxxy_1[j];

                    tg_xzz_xxxxz_0[j] = pb_x * tg_zz_xxxxz_0[j] + wp_x[j] * tg_zz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_zz_xxxz_1[j];

                    tg_xzz_xxxyy_0[j] = pb_x * tg_zz_xxxyy_0[j] + wp_x[j] * tg_zz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxyy_1[j];

                    tg_xzz_xxxyz_0[j] = pb_x * tg_zz_xxxyz_0[j] + wp_x[j] * tg_zz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_zz_xxyz_1[j];

                    tg_xzz_xxxzz_0[j] = pb_x * tg_zz_xxxzz_0[j] + wp_x[j] * tg_zz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxzz_1[j];

                    tg_xzz_xxyyy_0[j] = pb_x * tg_zz_xxyyy_0[j] + wp_x[j] * tg_zz_xxyyy_1[j] + fl1_fxn * tg_zz_xyyy_1[j];

                    tg_xzz_xxyyz_0[j] = pb_x * tg_zz_xxyyz_0[j] + wp_x[j] * tg_zz_xxyyz_1[j] + fl1_fxn * tg_zz_xyyz_1[j];

                    tg_xzz_xxyzz_0[j] = pb_x * tg_zz_xxyzz_0[j] + wp_x[j] * tg_zz_xxyzz_1[j] + fl1_fxn * tg_zz_xyzz_1[j];

                    tg_xzz_xxzzz_0[j] = pb_x * tg_zz_xxzzz_0[j] + wp_x[j] * tg_zz_xxzzz_1[j] + fl1_fxn * tg_zz_xzzz_1[j];

                    tg_xzz_xyyyy_0[j] = pb_x * tg_zz_xyyyy_0[j] + wp_x[j] * tg_zz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_zz_yyyy_1[j];

                    tg_xzz_xyyyz_0[j] = pb_x * tg_zz_xyyyz_0[j] + wp_x[j] * tg_zz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyz_1[j];

                    tg_xzz_xyyzz_0[j] = pb_x * tg_zz_xyyzz_0[j] + wp_x[j] * tg_zz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_zz_yyzz_1[j];

                    tg_xzz_xyzzz_0[j] = pb_x * tg_zz_xyzzz_0[j] + wp_x[j] * tg_zz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_yzzz_1[j];

                    tg_xzz_xzzzz_0[j] = pb_x * tg_zz_xzzzz_0[j] + wp_x[j] * tg_zz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzz_1[j];

                    tg_xzz_yyyyy_0[j] = pb_x * tg_zz_yyyyy_0[j] + wp_x[j] * tg_zz_yyyyy_1[j];

                    tg_xzz_yyyyz_0[j] = pb_x * tg_zz_yyyyz_0[j] + wp_x[j] * tg_zz_yyyyz_1[j];

                    tg_xzz_yyyzz_0[j] = pb_x * tg_zz_yyyzz_0[j] + wp_x[j] * tg_zz_yyyzz_1[j];

                    tg_xzz_yyzzz_0[j] = pb_x * tg_zz_yyzzz_0[j] + wp_x[j] * tg_zz_yyzzz_1[j];

                    tg_xzz_yzzzz_0[j] = pb_x * tg_zz_yzzzz_0[j] + wp_x[j] * tg_zz_yzzzz_1[j];

                    tg_xzz_zzzzz_0[j] = pb_x * tg_zz_zzzzz_0[j] + wp_x[j] * tg_zz_zzzzz_1[j];

                    tg_yyy_xxxxx_0[j] = pb_y * tg_yy_xxxxx_0[j] + wp_y[j] * tg_yy_xxxxx_1[j] + fl1_fx * tg_y_xxxxx_0[j] - fl1_fx * fl1_fza * tg_y_xxxxx_1[j];

                    tg_yyy_xxxxy_0[j] = pb_y * tg_yy_xxxxy_0[j] + wp_y[j] * tg_yy_xxxxy_1[j] + fl1_fx * tg_y_xxxxy_0[j] - fl1_fx * fl1_fza * tg_y_xxxxy_1[j] + 0.5 * fl1_fxn * tg_yy_xxxx_1[j];

                    tg_yyy_xxxxz_0[j] = pb_y * tg_yy_xxxxz_0[j] + wp_y[j] * tg_yy_xxxxz_1[j] + fl1_fx * tg_y_xxxxz_0[j] - fl1_fx * fl1_fza * tg_y_xxxxz_1[j];

                    tg_yyy_xxxyy_0[j] = pb_y * tg_yy_xxxyy_0[j] + wp_y[j] * tg_yy_xxxyy_1[j] + fl1_fx * tg_y_xxxyy_0[j] - fl1_fx * fl1_fza * tg_y_xxxyy_1[j] + fl1_fxn * tg_yy_xxxy_1[j];

                    tg_yyy_xxxyz_0[j] = pb_y * tg_yy_xxxyz_0[j] + wp_y[j] * tg_yy_xxxyz_1[j] + fl1_fx * tg_y_xxxyz_0[j] - fl1_fx * fl1_fza * tg_y_xxxyz_1[j] + 0.5 * fl1_fxn * tg_yy_xxxz_1[j];

                    tg_yyy_xxxzz_0[j] = pb_y * tg_yy_xxxzz_0[j] + wp_y[j] * tg_yy_xxxzz_1[j] + fl1_fx * tg_y_xxxzz_0[j] - fl1_fx * fl1_fza * tg_y_xxxzz_1[j];

                    tg_yyy_xxyyy_0[j] = pb_y * tg_yy_xxyyy_0[j] + wp_y[j] * tg_yy_xxyyy_1[j] + fl1_fx * tg_y_xxyyy_0[j] - fl1_fx * fl1_fza * tg_y_xxyyy_1[j] + 1.5 * fl1_fxn * tg_yy_xxyy_1[j];

                    tg_yyy_xxyyz_0[j] = pb_y * tg_yy_xxyyz_0[j] + wp_y[j] * tg_yy_xxyyz_1[j] + fl1_fx * tg_y_xxyyz_0[j] - fl1_fx * fl1_fza * tg_y_xxyyz_1[j] + fl1_fxn * tg_yy_xxyz_1[j];

                    tg_yyy_xxyzz_0[j] = pb_y * tg_yy_xxyzz_0[j] + wp_y[j] * tg_yy_xxyzz_1[j] + fl1_fx * tg_y_xxyzz_0[j] - fl1_fx * fl1_fza * tg_y_xxyzz_1[j] + 0.5 * fl1_fxn * tg_yy_xxzz_1[j];

                    tg_yyy_xxzzz_0[j] = pb_y * tg_yy_xxzzz_0[j] + wp_y[j] * tg_yy_xxzzz_1[j] + fl1_fx * tg_y_xxzzz_0[j] - fl1_fx * fl1_fza * tg_y_xxzzz_1[j];

                    tg_yyy_xyyyy_0[j] = pb_y * tg_yy_xyyyy_0[j] + wp_y[j] * tg_yy_xyyyy_1[j] + fl1_fx * tg_y_xyyyy_0[j] - fl1_fx * fl1_fza * tg_y_xyyyy_1[j] + 2.0 * fl1_fxn * tg_yy_xyyy_1[j];

                    tg_yyy_xyyyz_0[j] = pb_y * tg_yy_xyyyz_0[j] + wp_y[j] * tg_yy_xyyyz_1[j] + fl1_fx * tg_y_xyyyz_0[j] - fl1_fx * fl1_fza * tg_y_xyyyz_1[j] + 1.5 * fl1_fxn * tg_yy_xyyz_1[j];

                    tg_yyy_xyyzz_0[j] = pb_y * tg_yy_xyyzz_0[j] + wp_y[j] * tg_yy_xyyzz_1[j] + fl1_fx * tg_y_xyyzz_0[j] - fl1_fx * fl1_fza * tg_y_xyyzz_1[j] + fl1_fxn * tg_yy_xyzz_1[j];

                    tg_yyy_xyzzz_0[j] = pb_y * tg_yy_xyzzz_0[j] + wp_y[j] * tg_yy_xyzzz_1[j] + fl1_fx * tg_y_xyzzz_0[j] - fl1_fx * fl1_fza * tg_y_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yy_xzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSH_140_210(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (140,210)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yy_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 77); 

                auto tg_yy_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 78); 

                auto tg_yy_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 79); 

                auto tg_yy_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 80); 

                auto tg_yy_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 81); 

                auto tg_yy_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 82); 

                auto tg_yy_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 83); 

                auto tg_yz_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 84); 

                auto tg_yz_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 85); 

                auto tg_yz_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 86); 

                auto tg_yz_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 87); 

                auto tg_yz_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 88); 

                auto tg_yz_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 89); 

                auto tg_yz_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 90); 

                auto tg_yz_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 91); 

                auto tg_yz_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 92); 

                auto tg_yz_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 93); 

                auto tg_yz_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 94); 

                auto tg_yz_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 95); 

                auto tg_yz_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 96); 

                auto tg_yz_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 97); 

                auto tg_yz_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 98); 

                auto tg_yz_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 99); 

                auto tg_yz_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 100); 

                auto tg_yz_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 101); 

                auto tg_yz_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 102); 

                auto tg_yz_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 103); 

                auto tg_yz_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 104); 

                auto tg_zz_xxxxx_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 105); 

                auto tg_zz_xxxxy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 106); 

                auto tg_zz_xxxxz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 107); 

                auto tg_zz_xxxyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 108); 

                auto tg_zz_xxxyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 109); 

                auto tg_zz_xxxzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 110); 

                auto tg_zz_xxyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 111); 

                auto tg_zz_xxyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 112); 

                auto tg_zz_xxyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 113); 

                auto tg_zz_xxzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 114); 

                auto tg_zz_xyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 115); 

                auto tg_zz_xyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 116); 

                auto tg_zz_xyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 117); 

                auto tg_zz_xyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 118); 

                auto tg_zz_xzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 119); 

                auto tg_zz_yyyyy_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 120); 

                auto tg_zz_yyyyz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 121); 

                auto tg_zz_yyyzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 122); 

                auto tg_zz_yyzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 123); 

                auto tg_zz_yzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 124); 

                auto tg_zz_zzzzz_0 = primBuffer.data(pidx_g_2_5_m0 + 126 * idx + 125); 

                auto tg_yy_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 77); 

                auto tg_yy_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 78); 

                auto tg_yy_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 79); 

                auto tg_yy_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 80); 

                auto tg_yy_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 81); 

                auto tg_yy_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 82); 

                auto tg_yy_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 83); 

                auto tg_yz_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 84); 

                auto tg_yz_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 85); 

                auto tg_yz_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 86); 

                auto tg_yz_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 87); 

                auto tg_yz_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 88); 

                auto tg_yz_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 89); 

                auto tg_yz_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 90); 

                auto tg_yz_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 91); 

                auto tg_yz_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 92); 

                auto tg_yz_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 93); 

                auto tg_yz_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 94); 

                auto tg_yz_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 95); 

                auto tg_yz_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 96); 

                auto tg_yz_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 97); 

                auto tg_yz_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 98); 

                auto tg_yz_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 99); 

                auto tg_yz_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 100); 

                auto tg_yz_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 101); 

                auto tg_yz_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 102); 

                auto tg_yz_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 103); 

                auto tg_yz_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 104); 

                auto tg_zz_xxxxx_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 105); 

                auto tg_zz_xxxxy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 106); 

                auto tg_zz_xxxxz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 107); 

                auto tg_zz_xxxyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 108); 

                auto tg_zz_xxxyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 109); 

                auto tg_zz_xxxzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 110); 

                auto tg_zz_xxyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 111); 

                auto tg_zz_xxyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 112); 

                auto tg_zz_xxyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 113); 

                auto tg_zz_xxzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 114); 

                auto tg_zz_xyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 115); 

                auto tg_zz_xyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 116); 

                auto tg_zz_xyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 117); 

                auto tg_zz_xyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 118); 

                auto tg_zz_xzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 119); 

                auto tg_zz_yyyyy_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 120); 

                auto tg_zz_yyyyz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 121); 

                auto tg_zz_yyyzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 122); 

                auto tg_zz_yyzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 123); 

                auto tg_zz_yzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 124); 

                auto tg_zz_zzzzz_1 = primBuffer.data(pidx_g_2_5_m1 + 126 * idx + 125); 

                auto tg_y_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 35); 

                auto tg_y_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 36); 

                auto tg_y_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 37); 

                auto tg_y_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 38); 

                auto tg_y_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 39); 

                auto tg_y_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 40); 

                auto tg_y_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 41); 

                auto tg_z_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 42); 

                auto tg_z_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 43); 

                auto tg_z_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 44); 

                auto tg_z_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 45); 

                auto tg_z_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 46); 

                auto tg_z_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 47); 

                auto tg_z_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 48); 

                auto tg_z_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 49); 

                auto tg_z_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 50); 

                auto tg_z_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 51); 

                auto tg_z_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 52); 

                auto tg_z_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 53); 

                auto tg_z_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 54); 

                auto tg_z_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 55); 

                auto tg_z_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 56); 

                auto tg_z_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 57); 

                auto tg_z_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 58); 

                auto tg_z_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 59); 

                auto tg_z_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 60); 

                auto tg_z_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 61); 

                auto tg_z_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 62); 

                auto tg_y_xzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 35); 

                auto tg_y_yyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 36); 

                auto tg_y_yyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 37); 

                auto tg_y_yyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 38); 

                auto tg_y_yyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 39); 

                auto tg_y_yzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 40); 

                auto tg_y_zzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 41); 

                auto tg_z_xxxxx_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 42); 

                auto tg_z_xxxxy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 43); 

                auto tg_z_xxxxz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 44); 

                auto tg_z_xxxyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 45); 

                auto tg_z_xxxyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 46); 

                auto tg_z_xxxzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 47); 

                auto tg_z_xxyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 48); 

                auto tg_z_xxyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 49); 

                auto tg_z_xxyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 50); 

                auto tg_z_xxzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 51); 

                auto tg_z_xyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 52); 

                auto tg_z_xyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 53); 

                auto tg_z_xyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 54); 

                auto tg_z_xyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 55); 

                auto tg_z_xzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 56); 

                auto tg_z_yyyyy_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 57); 

                auto tg_z_yyyyz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 58); 

                auto tg_z_yyyzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 59); 

                auto tg_z_yyzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 60); 

                auto tg_z_yzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 61); 

                auto tg_z_zzzzz_1 = primBuffer.data(pidx_g_1_5_m1 + 63 * idx + 62); 

                auto tg_yy_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 55); 

                auto tg_yy_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 56); 

                auto tg_yy_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 57); 

                auto tg_yy_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 58); 

                auto tg_yy_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 59); 

                auto tg_yz_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 60); 

                auto tg_yz_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 61); 

                auto tg_yz_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 62); 

                auto tg_yz_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 63); 

                auto tg_yz_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 64); 

                auto tg_yz_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 65); 

                auto tg_yz_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 66); 

                auto tg_yz_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 67); 

                auto tg_yz_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 68); 

                auto tg_yz_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 69); 

                auto tg_yz_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 70); 

                auto tg_yz_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 71); 

                auto tg_yz_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 72); 

                auto tg_yz_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 73); 

                auto tg_yz_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 74); 

                auto tg_zz_xxxx_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 75); 

                auto tg_zz_xxxy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 76); 

                auto tg_zz_xxxz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 77); 

                auto tg_zz_xxyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 78); 

                auto tg_zz_xxyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 79); 

                auto tg_zz_xxzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 80); 

                auto tg_zz_xyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 81); 

                auto tg_zz_xyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 82); 

                auto tg_zz_xyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 83); 

                auto tg_zz_xzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 84); 

                auto tg_zz_yyyy_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 85); 

                auto tg_zz_yyyz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 86); 

                auto tg_zz_yyzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 87); 

                auto tg_zz_yzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 88); 

                auto tg_zz_zzzz_1 = primBuffer.data(pidx_g_2_4_m1 + 90 * idx + 89); 

                // set up pointers to integrals

                auto tg_yyy_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 140); 

                auto tg_yyy_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 141); 

                auto tg_yyy_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 142); 

                auto tg_yyy_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 143); 

                auto tg_yyy_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 144); 

                auto tg_yyy_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 145); 

                auto tg_yyy_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 146); 

                auto tg_yyz_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 147); 

                auto tg_yyz_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 148); 

                auto tg_yyz_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 149); 

                auto tg_yyz_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 150); 

                auto tg_yyz_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 151); 

                auto tg_yyz_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 152); 

                auto tg_yyz_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 153); 

                auto tg_yyz_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 154); 

                auto tg_yyz_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 155); 

                auto tg_yyz_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 156); 

                auto tg_yyz_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 157); 

                auto tg_yyz_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 158); 

                auto tg_yyz_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 159); 

                auto tg_yyz_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 160); 

                auto tg_yyz_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 161); 

                auto tg_yyz_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 162); 

                auto tg_yyz_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 163); 

                auto tg_yyz_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 164); 

                auto tg_yyz_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 165); 

                auto tg_yyz_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 166); 

                auto tg_yyz_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 167); 

                auto tg_yzz_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 168); 

                auto tg_yzz_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 169); 

                auto tg_yzz_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 170); 

                auto tg_yzz_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 171); 

                auto tg_yzz_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 172); 

                auto tg_yzz_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 173); 

                auto tg_yzz_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 174); 

                auto tg_yzz_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 175); 

                auto tg_yzz_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 176); 

                auto tg_yzz_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 177); 

                auto tg_yzz_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 178); 

                auto tg_yzz_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 179); 

                auto tg_yzz_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 180); 

                auto tg_yzz_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 181); 

                auto tg_yzz_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 182); 

                auto tg_yzz_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 183); 

                auto tg_yzz_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 184); 

                auto tg_yzz_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 185); 

                auto tg_yzz_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 186); 

                auto tg_yzz_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 187); 

                auto tg_yzz_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 188); 

                auto tg_zzz_xxxxx_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 189); 

                auto tg_zzz_xxxxy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 190); 

                auto tg_zzz_xxxxz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 191); 

                auto tg_zzz_xxxyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 192); 

                auto tg_zzz_xxxyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 193); 

                auto tg_zzz_xxxzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 194); 

                auto tg_zzz_xxyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 195); 

                auto tg_zzz_xxyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 196); 

                auto tg_zzz_xxyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 197); 

                auto tg_zzz_xxzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 198); 

                auto tg_zzz_xyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 199); 

                auto tg_zzz_xyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 200); 

                auto tg_zzz_xyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 201); 

                auto tg_zzz_xyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 202); 

                auto tg_zzz_xzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 203); 

                auto tg_zzz_yyyyy_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 204); 

                auto tg_zzz_yyyyz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 205); 

                auto tg_zzz_yyyzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 206); 

                auto tg_zzz_yyzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 207); 

                auto tg_zzz_yzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 208); 

                auto tg_zzz_zzzzz_0 = primBuffer.data(pidx_g_3_5_m0 + 210 * idx + 209); 

                // Batch of Integrals (140,210)

                #pragma omp simd aligned(fxn, fza, tg_y_xzzzz_0, tg_y_xzzzz_1, tg_y_yyyyy_0, tg_y_yyyyy_1, \
                                         tg_y_yyyyz_0, tg_y_yyyyz_1, tg_y_yyyzz_0, tg_y_yyyzz_1, tg_y_yyzzz_0, tg_y_yyzzz_1, \
                                         tg_y_yzzzz_0, tg_y_yzzzz_1, tg_y_zzzzz_0, tg_y_zzzzz_1, tg_yy_xzzzz_0, \
                                         tg_yy_xzzzz_1, tg_yy_yyyy_1, tg_yy_yyyyy_0, tg_yy_yyyyy_1, tg_yy_yyyyz_0, \
                                         tg_yy_yyyyz_1, tg_yy_yyyz_1, tg_yy_yyyzz_0, tg_yy_yyyzz_1, tg_yy_yyzz_1, \
                                         tg_yy_yyzzz_0, tg_yy_yyzzz_1, tg_yy_yzzz_1, tg_yy_yzzzz_0, tg_yy_yzzzz_1, \
                                         tg_yy_zzzz_1, tg_yy_zzzzz_0, tg_yy_zzzzz_1, tg_yyy_xzzzz_0, tg_yyy_yyyyy_0, \
                                         tg_yyy_yyyyz_0, tg_yyy_yyyzz_0, tg_yyy_yyzzz_0, tg_yyy_yzzzz_0, tg_yyy_zzzzz_0, \
                                         tg_yyz_xxxxx_0, tg_yyz_xxxxy_0, tg_yyz_xxxxz_0, tg_yyz_xxxyy_0, tg_yyz_xxxyz_0, \
                                         tg_yyz_xxxzz_0, tg_yyz_xxyyy_0, tg_yyz_xxyyz_0, tg_yyz_xxyzz_0, tg_yyz_xxzzz_0, \
                                         tg_yyz_xyyyy_0, tg_yyz_xyyyz_0, tg_yyz_xyyzz_0, tg_yyz_xyzzz_0, tg_yyz_xzzzz_0, \
                                         tg_yyz_yyyyy_0, tg_yyz_yyyyz_0, tg_yyz_yyyzz_0, tg_yyz_yyzzz_0, tg_yyz_yzzzz_0, \
                                         tg_yyz_zzzzz_0, tg_yz_xxxx_1, tg_yz_xxxxx_0, tg_yz_xxxxx_1, tg_yz_xxxxy_0, \
                                         tg_yz_xxxxy_1, tg_yz_xxxxz_0, tg_yz_xxxxz_1, tg_yz_xxxy_1, tg_yz_xxxyy_0, \
                                         tg_yz_xxxyy_1, tg_yz_xxxyz_0, tg_yz_xxxyz_1, tg_yz_xxxz_1, tg_yz_xxxzz_0, \
                                         tg_yz_xxxzz_1, tg_yz_xxyy_1, tg_yz_xxyyy_0, tg_yz_xxyyy_1, tg_yz_xxyyz_0, \
                                         tg_yz_xxyyz_1, tg_yz_xxyz_1, tg_yz_xxyzz_0, tg_yz_xxyzz_1, tg_yz_xxzz_1, \
                                         tg_yz_xxzzz_0, tg_yz_xxzzz_1, tg_yz_xyyy_1, tg_yz_xyyyy_0, tg_yz_xyyyy_1, \
                                         tg_yz_xyyyz_0, tg_yz_xyyyz_1, tg_yz_xyyz_1, tg_yz_xyyzz_0, tg_yz_xyyzz_1, \
                                         tg_yz_xyzz_1, tg_yz_xyzzz_0, tg_yz_xyzzz_1, tg_yz_xzzz_1, tg_yz_xzzzz_0, \
                                         tg_yz_xzzzz_1, tg_yz_yyyy_1, tg_yz_yyyyy_0, tg_yz_yyyyy_1, tg_yz_yyyyz_0, \
                                         tg_yz_yyyyz_1, tg_yz_yyyz_1, tg_yz_yyyzz_0, tg_yz_yyyzz_1, tg_yz_yyzz_1, \
                                         tg_yz_yyzzz_0, tg_yz_yyzzz_1, tg_yz_yzzz_1, tg_yz_yzzzz_0, tg_yz_yzzzz_1, \
                                         tg_yz_zzzz_1, tg_yz_zzzzz_0, tg_yz_zzzzz_1, tg_yzz_xxxxx_0, tg_yzz_xxxxy_0, \
                                         tg_yzz_xxxxz_0, tg_yzz_xxxyy_0, tg_yzz_xxxyz_0, tg_yzz_xxxzz_0, tg_yzz_xxyyy_0, \
                                         tg_yzz_xxyyz_0, tg_yzz_xxyzz_0, tg_yzz_xxzzz_0, tg_yzz_xyyyy_0, tg_yzz_xyyyz_0, \
                                         tg_yzz_xyyzz_0, tg_yzz_xyzzz_0, tg_yzz_xzzzz_0, tg_yzz_yyyyy_0, tg_yzz_yyyyz_0, \
                                         tg_yzz_yyyzz_0, tg_yzz_yyzzz_0, tg_yzz_yzzzz_0, tg_yzz_zzzzz_0, tg_z_xxxxx_0, \
                                         tg_z_xxxxx_1, tg_z_xxxxy_0, tg_z_xxxxy_1, tg_z_xxxxz_0, tg_z_xxxxz_1, tg_z_xxxyy_0, \
                                         tg_z_xxxyy_1, tg_z_xxxyz_0, tg_z_xxxyz_1, tg_z_xxxzz_0, tg_z_xxxzz_1, tg_z_xxyyy_0, \
                                         tg_z_xxyyy_1, tg_z_xxyyz_0, tg_z_xxyyz_1, tg_z_xxyzz_0, tg_z_xxyzz_1, tg_z_xxzzz_0, \
                                         tg_z_xxzzz_1, tg_z_xyyyy_0, tg_z_xyyyy_1, tg_z_xyyyz_0, tg_z_xyyyz_1, tg_z_xyyzz_0, \
                                         tg_z_xyyzz_1, tg_z_xyzzz_0, tg_z_xyzzz_1, tg_z_xzzzz_0, tg_z_xzzzz_1, tg_z_yyyyy_0, \
                                         tg_z_yyyyy_1, tg_z_yyyyz_0, tg_z_yyyyz_1, tg_z_yyyzz_0, tg_z_yyyzz_1, tg_z_yyzzz_0, \
                                         tg_z_yyzzz_1, tg_z_yzzzz_0, tg_z_yzzzz_1, tg_z_zzzzz_0, tg_z_zzzzz_1, tg_zz_xxxx_1, \
                                         tg_zz_xxxxx_0, tg_zz_xxxxx_1, tg_zz_xxxxy_0, tg_zz_xxxxy_1, tg_zz_xxxxz_0, \
                                         tg_zz_xxxxz_1, tg_zz_xxxy_1, tg_zz_xxxyy_0, tg_zz_xxxyy_1, tg_zz_xxxyz_0, \
                                         tg_zz_xxxyz_1, tg_zz_xxxz_1, tg_zz_xxxzz_0, tg_zz_xxxzz_1, tg_zz_xxyy_1, \
                                         tg_zz_xxyyy_0, tg_zz_xxyyy_1, tg_zz_xxyyz_0, tg_zz_xxyyz_1, tg_zz_xxyz_1, \
                                         tg_zz_xxyzz_0, tg_zz_xxyzz_1, tg_zz_xxzz_1, tg_zz_xxzzz_0, tg_zz_xxzzz_1, \
                                         tg_zz_xyyy_1, tg_zz_xyyyy_0, tg_zz_xyyyy_1, tg_zz_xyyyz_0, tg_zz_xyyyz_1, \
                                         tg_zz_xyyz_1, tg_zz_xyyzz_0, tg_zz_xyyzz_1, tg_zz_xyzz_1, tg_zz_xyzzz_0, \
                                         tg_zz_xyzzz_1, tg_zz_xzzz_1, tg_zz_xzzzz_0, tg_zz_xzzzz_1, tg_zz_yyyy_1, \
                                         tg_zz_yyyyy_0, tg_zz_yyyyy_1, tg_zz_yyyyz_0, tg_zz_yyyyz_1, tg_zz_yyyz_1, \
                                         tg_zz_yyyzz_0, tg_zz_yyyzz_1, tg_zz_yyzz_1, tg_zz_yyzzz_0, tg_zz_yyzzz_1, \
                                         tg_zz_yzzz_1, tg_zz_yzzzz_0, tg_zz_yzzzz_1, tg_zz_zzzz_1, tg_zz_zzzzz_0, \
                                         tg_zz_zzzzz_1, tg_zzz_xxxxx_0, tg_zzz_xxxxy_0, tg_zzz_xxxxz_0, tg_zzz_xxxyy_0, \
                                         tg_zzz_xxxyz_0, tg_zzz_xxxzz_0, tg_zzz_xxyyy_0, tg_zzz_xxyyz_0, tg_zzz_xxyzz_0, \
                                         tg_zzz_xxzzz_0, tg_zzz_xyyyy_0, tg_zzz_xyyyz_0, tg_zzz_xyyzz_0, tg_zzz_xyzzz_0, \
                                         tg_zzz_xzzzz_0, tg_zzz_yyyyy_0, tg_zzz_yyyyz_0, tg_zzz_yyyzz_0, tg_zzz_yyzzz_0, \
                                         tg_zzz_yzzzz_0, tg_zzz_zzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyy_xzzzz_0[j] = pb_y * tg_yy_xzzzz_0[j] + wp_y[j] * tg_yy_xzzzz_1[j] + fl1_fx * tg_y_xzzzz_0[j] - fl1_fx * fl1_fza * tg_y_xzzzz_1[j];

                    tg_yyy_yyyyy_0[j] = pb_y * tg_yy_yyyyy_0[j] + wp_y[j] * tg_yy_yyyyy_1[j] + fl1_fx * tg_y_yyyyy_0[j] - fl1_fx * fl1_fza * tg_y_yyyyy_1[j] + 2.5 * fl1_fxn * tg_yy_yyyy_1[j];

                    tg_yyy_yyyyz_0[j] = pb_y * tg_yy_yyyyz_0[j] + wp_y[j] * tg_yy_yyyyz_1[j] + fl1_fx * tg_y_yyyyz_0[j] - fl1_fx * fl1_fza * tg_y_yyyyz_1[j] + 2.0 * fl1_fxn * tg_yy_yyyz_1[j];

                    tg_yyy_yyyzz_0[j] = pb_y * tg_yy_yyyzz_0[j] + wp_y[j] * tg_yy_yyyzz_1[j] + fl1_fx * tg_y_yyyzz_0[j] - fl1_fx * fl1_fza * tg_y_yyyzz_1[j] + 1.5 * fl1_fxn * tg_yy_yyzz_1[j];

                    tg_yyy_yyzzz_0[j] = pb_y * tg_yy_yyzzz_0[j] + wp_y[j] * tg_yy_yyzzz_1[j] + fl1_fx * tg_y_yyzzz_0[j] - fl1_fx * fl1_fza * tg_y_yyzzz_1[j] + fl1_fxn * tg_yy_yzzz_1[j];

                    tg_yyy_yzzzz_0[j] = pb_y * tg_yy_yzzzz_0[j] + wp_y[j] * tg_yy_yzzzz_1[j] + fl1_fx * tg_y_yzzzz_0[j] - fl1_fx * fl1_fza * tg_y_yzzzz_1[j] + 0.5 * fl1_fxn * tg_yy_zzzz_1[j];

                    tg_yyy_zzzzz_0[j] = pb_y * tg_yy_zzzzz_0[j] + wp_y[j] * tg_yy_zzzzz_1[j] + fl1_fx * tg_y_zzzzz_0[j] - fl1_fx * fl1_fza * tg_y_zzzzz_1[j];

                    tg_yyz_xxxxx_0[j] = pb_y * tg_yz_xxxxx_0[j] + wp_y[j] * tg_yz_xxxxx_1[j] + 0.5 * fl1_fx * tg_z_xxxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxx_1[j];

                    tg_yyz_xxxxy_0[j] = pb_y * tg_yz_xxxxy_0[j] + wp_y[j] * tg_yz_xxxxy_1[j] + 0.5 * fl1_fx * tg_z_xxxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxy_1[j] + 0.5 * fl1_fxn * tg_yz_xxxx_1[j];

                    tg_yyz_xxxxz_0[j] = pb_y * tg_yz_xxxxz_0[j] + wp_y[j] * tg_yz_xxxxz_1[j] + 0.5 * fl1_fx * tg_z_xxxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxxz_1[j];

                    tg_yyz_xxxyy_0[j] = pb_y * tg_yz_xxxyy_0[j] + wp_y[j] * tg_yz_xxxyy_1[j] + 0.5 * fl1_fx * tg_z_xxxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyy_1[j] + fl1_fxn * tg_yz_xxxy_1[j];

                    tg_yyz_xxxyz_0[j] = pb_y * tg_yz_xxxyz_0[j] + wp_y[j] * tg_yz_xxxyz_1[j] + 0.5 * fl1_fx * tg_z_xxxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxyz_1[j] + 0.5 * fl1_fxn * tg_yz_xxxz_1[j];

                    tg_yyz_xxxzz_0[j] = pb_y * tg_yz_xxxzz_0[j] + wp_y[j] * tg_yz_xxxzz_1[j] + 0.5 * fl1_fx * tg_z_xxxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxxzz_1[j];

                    tg_yyz_xxyyy_0[j] = pb_y * tg_yz_xxyyy_0[j] + wp_y[j] * tg_yz_xxyyy_1[j] + 0.5 * fl1_fx * tg_z_xxyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyy_1[j] + 1.5 * fl1_fxn * tg_yz_xxyy_1[j];

                    tg_yyz_xxyyz_0[j] = pb_y * tg_yz_xxyyz_0[j] + wp_y[j] * tg_yz_xxyyz_1[j] + 0.5 * fl1_fx * tg_z_xxyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyyz_1[j] + fl1_fxn * tg_yz_xxyz_1[j];

                    tg_yyz_xxyzz_0[j] = pb_y * tg_yz_xxyzz_0[j] + wp_y[j] * tg_yz_xxyzz_1[j] + 0.5 * fl1_fx * tg_z_xxyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxyzz_1[j] + 0.5 * fl1_fxn * tg_yz_xxzz_1[j];

                    tg_yyz_xxzzz_0[j] = pb_y * tg_yz_xxzzz_0[j] + wp_y[j] * tg_yz_xxzzz_1[j] + 0.5 * fl1_fx * tg_z_xxzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xxzzz_1[j];

                    tg_yyz_xyyyy_0[j] = pb_y * tg_yz_xyyyy_0[j] + wp_y[j] * tg_yz_xyyyy_1[j] + 0.5 * fl1_fx * tg_z_xyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyy_1[j] + 2.0 * fl1_fxn * tg_yz_xyyy_1[j];

                    tg_yyz_xyyyz_0[j] = pb_y * tg_yz_xyyyz_0[j] + wp_y[j] * tg_yz_xyyyz_1[j] + 0.5 * fl1_fx * tg_z_xyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyyz_1[j] + 1.5 * fl1_fxn * tg_yz_xyyz_1[j];

                    tg_yyz_xyyzz_0[j] = pb_y * tg_yz_xyyzz_0[j] + wp_y[j] * tg_yz_xyyzz_1[j] + 0.5 * fl1_fx * tg_z_xyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyyzz_1[j] + fl1_fxn * tg_yz_xyzz_1[j];

                    tg_yyz_xyzzz_0[j] = pb_y * tg_yz_xyzzz_0[j] + wp_y[j] * tg_yz_xyzzz_1[j] + 0.5 * fl1_fx * tg_z_xyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yz_xzzz_1[j];

                    tg_yyz_xzzzz_0[j] = pb_y * tg_yz_xzzzz_0[j] + wp_y[j] * tg_yz_xzzzz_1[j] + 0.5 * fl1_fx * tg_z_xzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_xzzzz_1[j];

                    tg_yyz_yyyyy_0[j] = pb_y * tg_yz_yyyyy_0[j] + wp_y[j] * tg_yz_yyyyy_1[j] + 0.5 * fl1_fx * tg_z_yyyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyy_1[j] + 2.5 * fl1_fxn * tg_yz_yyyy_1[j];

                    tg_yyz_yyyyz_0[j] = pb_y * tg_yz_yyyyz_0[j] + wp_y[j] * tg_yz_yyyyz_1[j] + 0.5 * fl1_fx * tg_z_yyyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyyz_1[j] + 2.0 * fl1_fxn * tg_yz_yyyz_1[j];

                    tg_yyz_yyyzz_0[j] = pb_y * tg_yz_yyyzz_0[j] + wp_y[j] * tg_yz_yyyzz_1[j] + 0.5 * fl1_fx * tg_z_yyyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyyzz_1[j] + 1.5 * fl1_fxn * tg_yz_yyzz_1[j];

                    tg_yyz_yyzzz_0[j] = pb_y * tg_yz_yyzzz_0[j] + wp_y[j] * tg_yz_yyzzz_1[j] + 0.5 * fl1_fx * tg_z_yyzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yyzzz_1[j] + fl1_fxn * tg_yz_yzzz_1[j];

                    tg_yyz_yzzzz_0[j] = pb_y * tg_yz_yzzzz_0[j] + wp_y[j] * tg_yz_yzzzz_1[j] + 0.5 * fl1_fx * tg_z_yzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_yzzzz_1[j] + 0.5 * fl1_fxn * tg_yz_zzzz_1[j];

                    tg_yyz_zzzzz_0[j] = pb_y * tg_yz_zzzzz_0[j] + wp_y[j] * tg_yz_zzzzz_1[j] + 0.5 * fl1_fx * tg_z_zzzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_zzzzz_1[j];

                    tg_yzz_xxxxx_0[j] = pb_y * tg_zz_xxxxx_0[j] + wp_y[j] * tg_zz_xxxxx_1[j];

                    tg_yzz_xxxxy_0[j] = pb_y * tg_zz_xxxxy_0[j] + wp_y[j] * tg_zz_xxxxy_1[j] + 0.5 * fl1_fxn * tg_zz_xxxx_1[j];

                    tg_yzz_xxxxz_0[j] = pb_y * tg_zz_xxxxz_0[j] + wp_y[j] * tg_zz_xxxxz_1[j];

                    tg_yzz_xxxyy_0[j] = pb_y * tg_zz_xxxyy_0[j] + wp_y[j] * tg_zz_xxxyy_1[j] + fl1_fxn * tg_zz_xxxy_1[j];

                    tg_yzz_xxxyz_0[j] = pb_y * tg_zz_xxxyz_0[j] + wp_y[j] * tg_zz_xxxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxz_1[j];

                    tg_yzz_xxxzz_0[j] = pb_y * tg_zz_xxxzz_0[j] + wp_y[j] * tg_zz_xxxzz_1[j];

                    tg_yzz_xxyyy_0[j] = pb_y * tg_zz_xxyyy_0[j] + wp_y[j] * tg_zz_xxyyy_1[j] + 1.5 * fl1_fxn * tg_zz_xxyy_1[j];

                    tg_yzz_xxyyz_0[j] = pb_y * tg_zz_xxyyz_0[j] + wp_y[j] * tg_zz_xxyyz_1[j] + fl1_fxn * tg_zz_xxyz_1[j];

                    tg_yzz_xxyzz_0[j] = pb_y * tg_zz_xxyzz_0[j] + wp_y[j] * tg_zz_xxyzz_1[j] + 0.5 * fl1_fxn * tg_zz_xxzz_1[j];

                    tg_yzz_xxzzz_0[j] = pb_y * tg_zz_xxzzz_0[j] + wp_y[j] * tg_zz_xxzzz_1[j];

                    tg_yzz_xyyyy_0[j] = pb_y * tg_zz_xyyyy_0[j] + wp_y[j] * tg_zz_xyyyy_1[j] + 2.0 * fl1_fxn * tg_zz_xyyy_1[j];

                    tg_yzz_xyyyz_0[j] = pb_y * tg_zz_xyyyz_0[j] + wp_y[j] * tg_zz_xyyyz_1[j] + 1.5 * fl1_fxn * tg_zz_xyyz_1[j];

                    tg_yzz_xyyzz_0[j] = pb_y * tg_zz_xyyzz_0[j] + wp_y[j] * tg_zz_xyyzz_1[j] + fl1_fxn * tg_zz_xyzz_1[j];

                    tg_yzz_xyzzz_0[j] = pb_y * tg_zz_xyzzz_0[j] + wp_y[j] * tg_zz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zz_xzzz_1[j];

                    tg_yzz_xzzzz_0[j] = pb_y * tg_zz_xzzzz_0[j] + wp_y[j] * tg_zz_xzzzz_1[j];

                    tg_yzz_yyyyy_0[j] = pb_y * tg_zz_yyyyy_0[j] + wp_y[j] * tg_zz_yyyyy_1[j] + 2.5 * fl1_fxn * tg_zz_yyyy_1[j];

                    tg_yzz_yyyyz_0[j] = pb_y * tg_zz_yyyyz_0[j] + wp_y[j] * tg_zz_yyyyz_1[j] + 2.0 * fl1_fxn * tg_zz_yyyz_1[j];

                    tg_yzz_yyyzz_0[j] = pb_y * tg_zz_yyyzz_0[j] + wp_y[j] * tg_zz_yyyzz_1[j] + 1.5 * fl1_fxn * tg_zz_yyzz_1[j];

                    tg_yzz_yyzzz_0[j] = pb_y * tg_zz_yyzzz_0[j] + wp_y[j] * tg_zz_yyzzz_1[j] + fl1_fxn * tg_zz_yzzz_1[j];

                    tg_yzz_yzzzz_0[j] = pb_y * tg_zz_yzzzz_0[j] + wp_y[j] * tg_zz_yzzzz_1[j] + 0.5 * fl1_fxn * tg_zz_zzzz_1[j];

                    tg_yzz_zzzzz_0[j] = pb_y * tg_zz_zzzzz_0[j] + wp_y[j] * tg_zz_zzzzz_1[j];

                    tg_zzz_xxxxx_0[j] = pb_z * tg_zz_xxxxx_0[j] + wp_z[j] * tg_zz_xxxxx_1[j] + fl1_fx * tg_z_xxxxx_0[j] - fl1_fx * fl1_fza * tg_z_xxxxx_1[j];

                    tg_zzz_xxxxy_0[j] = pb_z * tg_zz_xxxxy_0[j] + wp_z[j] * tg_zz_xxxxy_1[j] + fl1_fx * tg_z_xxxxy_0[j] - fl1_fx * fl1_fza * tg_z_xxxxy_1[j];

                    tg_zzz_xxxxz_0[j] = pb_z * tg_zz_xxxxz_0[j] + wp_z[j] * tg_zz_xxxxz_1[j] + fl1_fx * tg_z_xxxxz_0[j] - fl1_fx * fl1_fza * tg_z_xxxxz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxx_1[j];

                    tg_zzz_xxxyy_0[j] = pb_z * tg_zz_xxxyy_0[j] + wp_z[j] * tg_zz_xxxyy_1[j] + fl1_fx * tg_z_xxxyy_0[j] - fl1_fx * fl1_fza * tg_z_xxxyy_1[j];

                    tg_zzz_xxxyz_0[j] = pb_z * tg_zz_xxxyz_0[j] + wp_z[j] * tg_zz_xxxyz_1[j] + fl1_fx * tg_z_xxxyz_0[j] - fl1_fx * fl1_fza * tg_z_xxxyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxxy_1[j];

                    tg_zzz_xxxzz_0[j] = pb_z * tg_zz_xxxzz_0[j] + wp_z[j] * tg_zz_xxxzz_1[j] + fl1_fx * tg_z_xxxzz_0[j] - fl1_fx * fl1_fza * tg_z_xxxzz_1[j] + fl1_fxn * tg_zz_xxxz_1[j];

                    tg_zzz_xxyyy_0[j] = pb_z * tg_zz_xxyyy_0[j] + wp_z[j] * tg_zz_xxyyy_1[j] + fl1_fx * tg_z_xxyyy_0[j] - fl1_fx * fl1_fza * tg_z_xxyyy_1[j];

                    tg_zzz_xxyyz_0[j] = pb_z * tg_zz_xxyyz_0[j] + wp_z[j] * tg_zz_xxyyz_1[j] + fl1_fx * tg_z_xxyyz_0[j] - fl1_fx * fl1_fza * tg_z_xxyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xxyy_1[j];

                    tg_zzz_xxyzz_0[j] = pb_z * tg_zz_xxyzz_0[j] + wp_z[j] * tg_zz_xxyzz_1[j] + fl1_fx * tg_z_xxyzz_0[j] - fl1_fx * fl1_fza * tg_z_xxyzz_1[j] + fl1_fxn * tg_zz_xxyz_1[j];

                    tg_zzz_xxzzz_0[j] = pb_z * tg_zz_xxzzz_0[j] + wp_z[j] * tg_zz_xxzzz_1[j] + fl1_fx * tg_z_xxzzz_0[j] - fl1_fx * fl1_fza * tg_z_xxzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xxzz_1[j];

                    tg_zzz_xyyyy_0[j] = pb_z * tg_zz_xyyyy_0[j] + wp_z[j] * tg_zz_xyyyy_1[j] + fl1_fx * tg_z_xyyyy_0[j] - fl1_fx * fl1_fza * tg_z_xyyyy_1[j];

                    tg_zzz_xyyyz_0[j] = pb_z * tg_zz_xyyyz_0[j] + wp_z[j] * tg_zz_xyyyz_1[j] + fl1_fx * tg_z_xyyyz_0[j] - fl1_fx * fl1_fza * tg_z_xyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_xyyy_1[j];

                    tg_zzz_xyyzz_0[j] = pb_z * tg_zz_xyyzz_0[j] + wp_z[j] * tg_zz_xyyzz_1[j] + fl1_fx * tg_z_xyyzz_0[j] - fl1_fx * fl1_fza * tg_z_xyyzz_1[j] + fl1_fxn * tg_zz_xyyz_1[j];

                    tg_zzz_xyzzz_0[j] = pb_z * tg_zz_xyzzz_0[j] + wp_z[j] * tg_zz_xyzzz_1[j] + fl1_fx * tg_z_xyzzz_0[j] - fl1_fx * fl1_fza * tg_z_xyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_xyzz_1[j];

                    tg_zzz_xzzzz_0[j] = pb_z * tg_zz_xzzzz_0[j] + wp_z[j] * tg_zz_xzzzz_1[j] + fl1_fx * tg_z_xzzzz_0[j] - fl1_fx * fl1_fza * tg_z_xzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_xzzz_1[j];

                    tg_zzz_yyyyy_0[j] = pb_z * tg_zz_yyyyy_0[j] + wp_z[j] * tg_zz_yyyyy_1[j] + fl1_fx * tg_z_yyyyy_0[j] - fl1_fx * fl1_fza * tg_z_yyyyy_1[j];

                    tg_zzz_yyyyz_0[j] = pb_z * tg_zz_yyyyz_0[j] + wp_z[j] * tg_zz_yyyyz_1[j] + fl1_fx * tg_z_yyyyz_0[j] - fl1_fx * fl1_fza * tg_z_yyyyz_1[j] + 0.5 * fl1_fxn * tg_zz_yyyy_1[j];

                    tg_zzz_yyyzz_0[j] = pb_z * tg_zz_yyyzz_0[j] + wp_z[j] * tg_zz_yyyzz_1[j] + fl1_fx * tg_z_yyyzz_0[j] - fl1_fx * fl1_fza * tg_z_yyyzz_1[j] + fl1_fxn * tg_zz_yyyz_1[j];

                    tg_zzz_yyzzz_0[j] = pb_z * tg_zz_yyzzz_0[j] + wp_z[j] * tg_zz_yyzzz_1[j] + fl1_fx * tg_z_yyzzz_0[j] - fl1_fx * fl1_fza * tg_z_yyzzz_1[j] + 1.5 * fl1_fxn * tg_zz_yyzz_1[j];

                    tg_zzz_yzzzz_0[j] = pb_z * tg_zz_yzzzz_0[j] + wp_z[j] * tg_zz_yzzzz_1[j] + fl1_fx * tg_z_yzzzz_0[j] - fl1_fx * fl1_fza * tg_z_yzzzz_1[j] + 2.0 * fl1_fxn * tg_zz_yzzz_1[j];

                    tg_zzz_zzzzz_0[j] = pb_z * tg_zz_zzzzz_0[j] + wp_z[j] * tg_zz_zzzzz_1[j] + fl1_fx * tg_z_zzzzz_0[j] - fl1_fx * fl1_fza * tg_z_zzzzz_1[j] + 2.5 * fl1_fxn * tg_zz_zzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

