//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "ElectronRepulsionRecFuncForLH.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSLSH(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSLSH_0_95(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_95_190(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_190_285(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_285_380(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_380_475(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_475_569(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_569_663(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_663_757(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_757_851(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSH_851_945(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSLSH_0_95(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxxxxxx_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx); 

                auto tg_xxxxxxx_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 1); 

                auto tg_xxxxxxx_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 2); 

                auto tg_xxxxxxx_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 3); 

                auto tg_xxxxxxx_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 4); 

                auto tg_xxxxxxx_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 5); 

                auto tg_xxxxxxx_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 6); 

                auto tg_xxxxxxx_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 7); 

                auto tg_xxxxxxx_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 8); 

                auto tg_xxxxxxx_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 9); 

                auto tg_xxxxxxx_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 10); 

                auto tg_xxxxxxx_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 11); 

                auto tg_xxxxxxx_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 12); 

                auto tg_xxxxxxx_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 13); 

                auto tg_xxxxxxx_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 14); 

                auto tg_xxxxxxx_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 15); 

                auto tg_xxxxxxx_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 16); 

                auto tg_xxxxxxx_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 17); 

                auto tg_xxxxxxx_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 18); 

                auto tg_xxxxxxx_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 19); 

                auto tg_xxxxxxx_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 20); 

                auto tg_xxxxxxy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 21); 

                auto tg_xxxxxxy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 22); 

                auto tg_xxxxxxy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 23); 

                auto tg_xxxxxxy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 24); 

                auto tg_xxxxxxy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 25); 

                auto tg_xxxxxxy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 26); 

                auto tg_xxxxxxy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 27); 

                auto tg_xxxxxxy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 28); 

                auto tg_xxxxxxy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 29); 

                auto tg_xxxxxxy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 30); 

                auto tg_xxxxxxy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 31); 

                auto tg_xxxxxxy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 32); 

                auto tg_xxxxxxy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 33); 

                auto tg_xxxxxxy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 34); 

                auto tg_xxxxxxy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 35); 

                auto tg_xxxxxxy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 36); 

                auto tg_xxxxxxy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 37); 

                auto tg_xxxxxxy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 38); 

                auto tg_xxxxxxy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 39); 

                auto tg_xxxxxxy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 40); 

                auto tg_xxxxxxy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 41); 

                auto tg_xxxxxxz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 42); 

                auto tg_xxxxxxz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 43); 

                auto tg_xxxxxxz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 44); 

                auto tg_xxxxxxz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 45); 

                auto tg_xxxxxxz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 46); 

                auto tg_xxxxxxz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 47); 

                auto tg_xxxxxxz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 48); 

                auto tg_xxxxxxz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 49); 

                auto tg_xxxxxxz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 50); 

                auto tg_xxxxxxz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 51); 

                auto tg_xxxxxxz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 52); 

                auto tg_xxxxxxz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 53); 

                auto tg_xxxxxxz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 54); 

                auto tg_xxxxxxz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 55); 

                auto tg_xxxxxxz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 56); 

                auto tg_xxxxxxz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 57); 

                auto tg_xxxxxxz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 58); 

                auto tg_xxxxxxz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 59); 

                auto tg_xxxxxxz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 60); 

                auto tg_xxxxxxz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 61); 

                auto tg_xxxxxxz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 62); 

                auto tg_xxxxxyy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 63); 

                auto tg_xxxxxyy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 64); 

                auto tg_xxxxxyy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 65); 

                auto tg_xxxxxyy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 66); 

                auto tg_xxxxxyy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 67); 

                auto tg_xxxxxyy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 68); 

                auto tg_xxxxxyy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 69); 

                auto tg_xxxxxyy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 70); 

                auto tg_xxxxxyy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 71); 

                auto tg_xxxxxyy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 72); 

                auto tg_xxxxxyy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 73); 

                auto tg_xxxxxyy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 74); 

                auto tg_xxxxxyy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 75); 

                auto tg_xxxxxyy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 76); 

                auto tg_xxxxxyy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 77); 

                auto tg_xxxxxyy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 78); 

                auto tg_xxxxxyy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 79); 

                auto tg_xxxxxyy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 80); 

                auto tg_xxxxxyy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 81); 

                auto tg_xxxxxyy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 82); 

                auto tg_xxxxxyy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 83); 

                auto tg_xxxxxyz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 84); 

                auto tg_xxxxxyz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 85); 

                auto tg_xxxxxyz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 86); 

                auto tg_xxxxxyz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 87); 

                auto tg_xxxxxyz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 88); 

                auto tg_xxxxxyz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 89); 

                auto tg_xxxxxyz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 90); 

                auto tg_xxxxxyz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 91); 

                auto tg_xxxxxyz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 92); 

                auto tg_xxxxxyz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 93); 

                auto tg_xxxxxyz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 94); 

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

                auto tg_xxxxxxx_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx); 

                auto tg_xxxxxxx_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 1); 

                auto tg_xxxxxxx_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 2); 

                auto tg_xxxxxxx_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 3); 

                auto tg_xxxxxxx_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 4); 

                auto tg_xxxxxxx_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 5); 

                auto tg_xxxxxxx_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 6); 

                auto tg_xxxxxxx_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 7); 

                auto tg_xxxxxxx_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 8); 

                auto tg_xxxxxxx_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 9); 

                auto tg_xxxxxxx_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 10); 

                auto tg_xxxxxxx_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 11); 

                auto tg_xxxxxxx_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 12); 

                auto tg_xxxxxxx_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 13); 

                auto tg_xxxxxxx_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 14); 

                auto tg_xxxxxxy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 15); 

                auto tg_xxxxxxy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 16); 

                auto tg_xxxxxxy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 17); 

                auto tg_xxxxxxy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 18); 

                auto tg_xxxxxxy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 19); 

                auto tg_xxxxxxy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 20); 

                auto tg_xxxxxxy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 21); 

                auto tg_xxxxxxy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 22); 

                auto tg_xxxxxxy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 23); 

                auto tg_xxxxxxy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 24); 

                auto tg_xxxxxxy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 25); 

                auto tg_xxxxxxy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 26); 

                auto tg_xxxxxxy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 27); 

                auto tg_xxxxxxy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 28); 

                auto tg_xxxxxxy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 29); 

                auto tg_xxxxxxz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 30); 

                auto tg_xxxxxxz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 31); 

                auto tg_xxxxxxz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 32); 

                auto tg_xxxxxxz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 33); 

                auto tg_xxxxxxz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 34); 

                auto tg_xxxxxxz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 35); 

                auto tg_xxxxxxz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 36); 

                auto tg_xxxxxxz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 37); 

                auto tg_xxxxxxz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 38); 

                auto tg_xxxxxxz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 39); 

                auto tg_xxxxxxz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 40); 

                auto tg_xxxxxxz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 41); 

                auto tg_xxxxxxz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 42); 

                auto tg_xxxxxxz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 43); 

                auto tg_xxxxxxz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 44); 

                auto tg_xxxxxyy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 45); 

                auto tg_xxxxxyy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 46); 

                auto tg_xxxxxyy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 47); 

                auto tg_xxxxxyy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 48); 

                auto tg_xxxxxyy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 49); 

                auto tg_xxxxxyy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 50); 

                auto tg_xxxxxyy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 51); 

                auto tg_xxxxxyy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 52); 

                auto tg_xxxxxyy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 53); 

                auto tg_xxxxxyy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 54); 

                auto tg_xxxxxyy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 55); 

                auto tg_xxxxxyy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 56); 

                auto tg_xxxxxyy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 57); 

                auto tg_xxxxxyy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 58); 

                auto tg_xxxxxyy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 59); 

                auto tg_xxxxxyz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 60); 

                auto tg_xxxxxyz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 61); 

                auto tg_xxxxxyz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 62); 

                auto tg_xxxxxyz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 63); 

                auto tg_xxxxxyz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 64); 

                auto tg_xxxxxyz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 65); 

                auto tg_xxxxxyz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 66); 

                auto tg_xxxxxyz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 67); 

                auto tg_xxxxxyz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 68); 

                auto tg_xxxxxyz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 69); 

                auto tg_xxxxxyz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 70); 

                // set up pointers to integrals

                auto tg_xxxxxxxx_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx); 

                auto tg_xxxxxxxx_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 1); 

                auto tg_xxxxxxxx_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 2); 

                auto tg_xxxxxxxx_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 3); 

                auto tg_xxxxxxxx_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 4); 

                auto tg_xxxxxxxx_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 5); 

                auto tg_xxxxxxxx_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 6); 

                auto tg_xxxxxxxx_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 7); 

                auto tg_xxxxxxxx_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 8); 

                auto tg_xxxxxxxx_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 9); 

                auto tg_xxxxxxxx_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 10); 

                auto tg_xxxxxxxx_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 11); 

                auto tg_xxxxxxxx_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 12); 

                auto tg_xxxxxxxx_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 13); 

                auto tg_xxxxxxxx_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 14); 

                auto tg_xxxxxxxx_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 15); 

                auto tg_xxxxxxxx_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 16); 

                auto tg_xxxxxxxx_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 17); 

                auto tg_xxxxxxxx_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 18); 

                auto tg_xxxxxxxx_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 19); 

                auto tg_xxxxxxxx_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 20); 

                auto tg_xxxxxxxy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 21); 

                auto tg_xxxxxxxy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 22); 

                auto tg_xxxxxxxy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 23); 

                auto tg_xxxxxxxy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 24); 

                auto tg_xxxxxxxy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 25); 

                auto tg_xxxxxxxy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 26); 

                auto tg_xxxxxxxy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 27); 

                auto tg_xxxxxxxy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 28); 

                auto tg_xxxxxxxy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 29); 

                auto tg_xxxxxxxy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 30); 

                auto tg_xxxxxxxy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 31); 

                auto tg_xxxxxxxy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 32); 

                auto tg_xxxxxxxy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 33); 

                auto tg_xxxxxxxy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 34); 

                auto tg_xxxxxxxy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 35); 

                auto tg_xxxxxxxy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 36); 

                auto tg_xxxxxxxy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 37); 

                auto tg_xxxxxxxy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 38); 

                auto tg_xxxxxxxy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 39); 

                auto tg_xxxxxxxy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 40); 

                auto tg_xxxxxxxy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 41); 

                auto tg_xxxxxxxz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 42); 

                auto tg_xxxxxxxz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 43); 

                auto tg_xxxxxxxz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 44); 

                auto tg_xxxxxxxz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 45); 

                auto tg_xxxxxxxz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 46); 

                auto tg_xxxxxxxz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 47); 

                auto tg_xxxxxxxz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 48); 

                auto tg_xxxxxxxz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 49); 

                auto tg_xxxxxxxz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 50); 

                auto tg_xxxxxxxz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 51); 

                auto tg_xxxxxxxz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 52); 

                auto tg_xxxxxxxz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 53); 

                auto tg_xxxxxxxz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 54); 

                auto tg_xxxxxxxz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 55); 

                auto tg_xxxxxxxz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 56); 

                auto tg_xxxxxxxz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 57); 

                auto tg_xxxxxxxz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 58); 

                auto tg_xxxxxxxz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 59); 

                auto tg_xxxxxxxz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 60); 

                auto tg_xxxxxxxz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 61); 

                auto tg_xxxxxxxz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 62); 

                auto tg_xxxxxxyy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 63); 

                auto tg_xxxxxxyy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 64); 

                auto tg_xxxxxxyy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 65); 

                auto tg_xxxxxxyy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 66); 

                auto tg_xxxxxxyy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 67); 

                auto tg_xxxxxxyy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 68); 

                auto tg_xxxxxxyy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 69); 

                auto tg_xxxxxxyy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 70); 

                auto tg_xxxxxxyy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 71); 

                auto tg_xxxxxxyy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 72); 

                auto tg_xxxxxxyy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 73); 

                auto tg_xxxxxxyy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 74); 

                auto tg_xxxxxxyy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 75); 

                auto tg_xxxxxxyy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 76); 

                auto tg_xxxxxxyy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 77); 

                auto tg_xxxxxxyy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 78); 

                auto tg_xxxxxxyy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 79); 

                auto tg_xxxxxxyy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 80); 

                auto tg_xxxxxxyy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 81); 

                auto tg_xxxxxxyy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 82); 

                auto tg_xxxxxxyy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 83); 

                auto tg_xxxxxxyz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 84); 

                auto tg_xxxxxxyz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 85); 

                auto tg_xxxxxxyz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 86); 

                auto tg_xxxxxxyz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 87); 

                auto tg_xxxxxxyz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 88); 

                auto tg_xxxxxxyz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 89); 

                auto tg_xxxxxxyz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 90); 

                auto tg_xxxxxxyz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 91); 

                auto tg_xxxxxxyz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 92); 

                auto tg_xxxxxxyz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 93); 

                auto tg_xxxxxxyz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 94); 

                // Batch of Integrals (0,95)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxx_xxxxx_0, tg_xxxxxx_xxxxx_1, tg_xxxxxx_xxxxy_0, \
                                         tg_xxxxxx_xxxxy_1, tg_xxxxxx_xxxxz_0, tg_xxxxxx_xxxxz_1, tg_xxxxxx_xxxyy_0, \
                                         tg_xxxxxx_xxxyy_1, tg_xxxxxx_xxxyz_0, tg_xxxxxx_xxxyz_1, tg_xxxxxx_xxxzz_0, \
                                         tg_xxxxxx_xxxzz_1, tg_xxxxxx_xxyyy_0, tg_xxxxxx_xxyyy_1, tg_xxxxxx_xxyyz_0, \
                                         tg_xxxxxx_xxyyz_1, tg_xxxxxx_xxyzz_0, tg_xxxxxx_xxyzz_1, tg_xxxxxx_xxzzz_0, \
                                         tg_xxxxxx_xxzzz_1, tg_xxxxxx_xyyyy_0, tg_xxxxxx_xyyyy_1, tg_xxxxxx_xyyyz_0, \
                                         tg_xxxxxx_xyyyz_1, tg_xxxxxx_xyyzz_0, tg_xxxxxx_xyyzz_1, tg_xxxxxx_xyzzz_0, \
                                         tg_xxxxxx_xyzzz_1, tg_xxxxxx_xzzzz_0, tg_xxxxxx_xzzzz_1, tg_xxxxxx_yyyyy_0, \
                                         tg_xxxxxx_yyyyy_1, tg_xxxxxx_yyyyz_0, tg_xxxxxx_yyyyz_1, tg_xxxxxx_yyyzz_0, \
                                         tg_xxxxxx_yyyzz_1, tg_xxxxxx_yyzzz_0, tg_xxxxxx_yyzzz_1, tg_xxxxxx_yzzzz_0, \
                                         tg_xxxxxx_yzzzz_1, tg_xxxxxx_zzzzz_0, tg_xxxxxx_zzzzz_1, tg_xxxxxxx_xxxx_1, \
                                         tg_xxxxxxx_xxxxx_0, tg_xxxxxxx_xxxxx_1, tg_xxxxxxx_xxxxy_0, tg_xxxxxxx_xxxxy_1, \
                                         tg_xxxxxxx_xxxxz_0, tg_xxxxxxx_xxxxz_1, tg_xxxxxxx_xxxy_1, tg_xxxxxxx_xxxyy_0, \
                                         tg_xxxxxxx_xxxyy_1, tg_xxxxxxx_xxxyz_0, tg_xxxxxxx_xxxyz_1, tg_xxxxxxx_xxxz_1, \
                                         tg_xxxxxxx_xxxzz_0, tg_xxxxxxx_xxxzz_1, tg_xxxxxxx_xxyy_1, tg_xxxxxxx_xxyyy_0, \
                                         tg_xxxxxxx_xxyyy_1, tg_xxxxxxx_xxyyz_0, tg_xxxxxxx_xxyyz_1, tg_xxxxxxx_xxyz_1, \
                                         tg_xxxxxxx_xxyzz_0, tg_xxxxxxx_xxyzz_1, tg_xxxxxxx_xxzz_1, tg_xxxxxxx_xxzzz_0, \
                                         tg_xxxxxxx_xxzzz_1, tg_xxxxxxx_xyyy_1, tg_xxxxxxx_xyyyy_0, tg_xxxxxxx_xyyyy_1, \
                                         tg_xxxxxxx_xyyyz_0, tg_xxxxxxx_xyyyz_1, tg_xxxxxxx_xyyz_1, tg_xxxxxxx_xyyzz_0, \
                                         tg_xxxxxxx_xyyzz_1, tg_xxxxxxx_xyzz_1, tg_xxxxxxx_xyzzz_0, tg_xxxxxxx_xyzzz_1, \
                                         tg_xxxxxxx_xzzz_1, tg_xxxxxxx_xzzzz_0, tg_xxxxxxx_xzzzz_1, tg_xxxxxxx_yyyy_1, \
                                         tg_xxxxxxx_yyyyy_0, tg_xxxxxxx_yyyyy_1, tg_xxxxxxx_yyyyz_0, tg_xxxxxxx_yyyyz_1, \
                                         tg_xxxxxxx_yyyz_1, tg_xxxxxxx_yyyzz_0, tg_xxxxxxx_yyyzz_1, tg_xxxxxxx_yyzz_1, \
                                         tg_xxxxxxx_yyzzz_0, tg_xxxxxxx_yyzzz_1, tg_xxxxxxx_yzzz_1, tg_xxxxxxx_yzzzz_0, \
                                         tg_xxxxxxx_yzzzz_1, tg_xxxxxxx_zzzz_1, tg_xxxxxxx_zzzzz_0, tg_xxxxxxx_zzzzz_1, \
                                         tg_xxxxxxxx_xxxxx_0, tg_xxxxxxxx_xxxxy_0, tg_xxxxxxxx_xxxxz_0, tg_xxxxxxxx_xxxyy_0, \
                                         tg_xxxxxxxx_xxxyz_0, tg_xxxxxxxx_xxxzz_0, tg_xxxxxxxx_xxyyy_0, tg_xxxxxxxx_xxyyz_0, \
                                         tg_xxxxxxxx_xxyzz_0, tg_xxxxxxxx_xxzzz_0, tg_xxxxxxxx_xyyyy_0, tg_xxxxxxxx_xyyyz_0, \
                                         tg_xxxxxxxx_xyyzz_0, tg_xxxxxxxx_xyzzz_0, tg_xxxxxxxx_xzzzz_0, tg_xxxxxxxx_yyyyy_0, \
                                         tg_xxxxxxxx_yyyyz_0, tg_xxxxxxxx_yyyzz_0, tg_xxxxxxxx_yyzzz_0, tg_xxxxxxxx_yzzzz_0, \
                                         tg_xxxxxxxx_zzzzz_0, tg_xxxxxxxy_xxxxx_0, tg_xxxxxxxy_xxxxy_0, tg_xxxxxxxy_xxxxz_0, \
                                         tg_xxxxxxxy_xxxyy_0, tg_xxxxxxxy_xxxyz_0, tg_xxxxxxxy_xxxzz_0, tg_xxxxxxxy_xxyyy_0, \
                                         tg_xxxxxxxy_xxyyz_0, tg_xxxxxxxy_xxyzz_0, tg_xxxxxxxy_xxzzz_0, tg_xxxxxxxy_xyyyy_0, \
                                         tg_xxxxxxxy_xyyyz_0, tg_xxxxxxxy_xyyzz_0, tg_xxxxxxxy_xyzzz_0, tg_xxxxxxxy_xzzzz_0, \
                                         tg_xxxxxxxy_yyyyy_0, tg_xxxxxxxy_yyyyz_0, tg_xxxxxxxy_yyyzz_0, tg_xxxxxxxy_yyzzz_0, \
                                         tg_xxxxxxxy_yzzzz_0, tg_xxxxxxxy_zzzzz_0, tg_xxxxxxxz_xxxxx_0, tg_xxxxxxxz_xxxxy_0, \
                                         tg_xxxxxxxz_xxxxz_0, tg_xxxxxxxz_xxxyy_0, tg_xxxxxxxz_xxxyz_0, tg_xxxxxxxz_xxxzz_0, \
                                         tg_xxxxxxxz_xxyyy_0, tg_xxxxxxxz_xxyyz_0, tg_xxxxxxxz_xxyzz_0, tg_xxxxxxxz_xxzzz_0, \
                                         tg_xxxxxxxz_xyyyy_0, tg_xxxxxxxz_xyyyz_0, tg_xxxxxxxz_xyyzz_0, tg_xxxxxxxz_xyzzz_0, \
                                         tg_xxxxxxxz_xzzzz_0, tg_xxxxxxxz_yyyyy_0, tg_xxxxxxxz_yyyyz_0, tg_xxxxxxxz_yyyzz_0, \
                                         tg_xxxxxxxz_yyzzz_0, tg_xxxxxxxz_yzzzz_0, tg_xxxxxxxz_zzzzz_0, tg_xxxxxxy_xxxx_1, \
                                         tg_xxxxxxy_xxxxx_0, tg_xxxxxxy_xxxxx_1, tg_xxxxxxy_xxxxy_0, tg_xxxxxxy_xxxxy_1, \
                                         tg_xxxxxxy_xxxxz_0, tg_xxxxxxy_xxxxz_1, tg_xxxxxxy_xxxy_1, tg_xxxxxxy_xxxyy_0, \
                                         tg_xxxxxxy_xxxyy_1, tg_xxxxxxy_xxxyz_0, tg_xxxxxxy_xxxyz_1, tg_xxxxxxy_xxxz_1, \
                                         tg_xxxxxxy_xxxzz_0, tg_xxxxxxy_xxxzz_1, tg_xxxxxxy_xxyy_1, tg_xxxxxxy_xxyyy_0, \
                                         tg_xxxxxxy_xxyyy_1, tg_xxxxxxy_xxyyz_0, tg_xxxxxxy_xxyyz_1, tg_xxxxxxy_xxyz_1, \
                                         tg_xxxxxxy_xxyzz_0, tg_xxxxxxy_xxyzz_1, tg_xxxxxxy_xxzz_1, tg_xxxxxxy_xxzzz_0, \
                                         tg_xxxxxxy_xxzzz_1, tg_xxxxxxy_xyyy_1, tg_xxxxxxy_xyyyy_0, tg_xxxxxxy_xyyyy_1, \
                                         tg_xxxxxxy_xyyyz_0, tg_xxxxxxy_xyyyz_1, tg_xxxxxxy_xyyz_1, tg_xxxxxxy_xyyzz_0, \
                                         tg_xxxxxxy_xyyzz_1, tg_xxxxxxy_xyzz_1, tg_xxxxxxy_xyzzz_0, tg_xxxxxxy_xyzzz_1, \
                                         tg_xxxxxxy_xzzz_1, tg_xxxxxxy_xzzzz_0, tg_xxxxxxy_xzzzz_1, tg_xxxxxxy_yyyy_1, \
                                         tg_xxxxxxy_yyyyy_0, tg_xxxxxxy_yyyyy_1, tg_xxxxxxy_yyyyz_0, tg_xxxxxxy_yyyyz_1, \
                                         tg_xxxxxxy_yyyz_1, tg_xxxxxxy_yyyzz_0, tg_xxxxxxy_yyyzz_1, tg_xxxxxxy_yyzz_1, \
                                         tg_xxxxxxy_yyzzz_0, tg_xxxxxxy_yyzzz_1, tg_xxxxxxy_yzzz_1, tg_xxxxxxy_yzzzz_0, \
                                         tg_xxxxxxy_yzzzz_1, tg_xxxxxxy_zzzz_1, tg_xxxxxxy_zzzzz_0, tg_xxxxxxy_zzzzz_1, \
                                         tg_xxxxxxyy_xxxxx_0, tg_xxxxxxyy_xxxxy_0, tg_xxxxxxyy_xxxxz_0, tg_xxxxxxyy_xxxyy_0, \
                                         tg_xxxxxxyy_xxxyz_0, tg_xxxxxxyy_xxxzz_0, tg_xxxxxxyy_xxyyy_0, tg_xxxxxxyy_xxyyz_0, \
                                         tg_xxxxxxyy_xxyzz_0, tg_xxxxxxyy_xxzzz_0, tg_xxxxxxyy_xyyyy_0, tg_xxxxxxyy_xyyyz_0, \
                                         tg_xxxxxxyy_xyyzz_0, tg_xxxxxxyy_xyzzz_0, tg_xxxxxxyy_xzzzz_0, tg_xxxxxxyy_yyyyy_0, \
                                         tg_xxxxxxyy_yyyyz_0, tg_xxxxxxyy_yyyzz_0, tg_xxxxxxyy_yyzzz_0, tg_xxxxxxyy_yzzzz_0, \
                                         tg_xxxxxxyy_zzzzz_0, tg_xxxxxxyz_xxxxx_0, tg_xxxxxxyz_xxxxy_0, tg_xxxxxxyz_xxxxz_0, \
                                         tg_xxxxxxyz_xxxyy_0, tg_xxxxxxyz_xxxyz_0, tg_xxxxxxyz_xxxzz_0, tg_xxxxxxyz_xxyyy_0, \
                                         tg_xxxxxxyz_xxyyz_0, tg_xxxxxxyz_xxyzz_0, tg_xxxxxxyz_xxzzz_0, tg_xxxxxxyz_xyyyy_0, \
                                         tg_xxxxxxz_xxxx_1, tg_xxxxxxz_xxxxx_0, tg_xxxxxxz_xxxxx_1, tg_xxxxxxz_xxxxy_0, \
                                         tg_xxxxxxz_xxxxy_1, tg_xxxxxxz_xxxxz_0, tg_xxxxxxz_xxxxz_1, tg_xxxxxxz_xxxy_1, \
                                         tg_xxxxxxz_xxxyy_0, tg_xxxxxxz_xxxyy_1, tg_xxxxxxz_xxxyz_0, tg_xxxxxxz_xxxyz_1, \
                                         tg_xxxxxxz_xxxz_1, tg_xxxxxxz_xxxzz_0, tg_xxxxxxz_xxxzz_1, tg_xxxxxxz_xxyy_1, \
                                         tg_xxxxxxz_xxyyy_0, tg_xxxxxxz_xxyyy_1, tg_xxxxxxz_xxyyz_0, tg_xxxxxxz_xxyyz_1, \
                                         tg_xxxxxxz_xxyz_1, tg_xxxxxxz_xxyzz_0, tg_xxxxxxz_xxyzz_1, tg_xxxxxxz_xxzz_1, \
                                         tg_xxxxxxz_xxzzz_0, tg_xxxxxxz_xxzzz_1, tg_xxxxxxz_xyyy_1, tg_xxxxxxz_xyyyy_0, \
                                         tg_xxxxxxz_xyyyy_1, tg_xxxxxxz_xyyyz_0, tg_xxxxxxz_xyyyz_1, tg_xxxxxxz_xyyz_1, \
                                         tg_xxxxxxz_xyyzz_0, tg_xxxxxxz_xyyzz_1, tg_xxxxxxz_xyzz_1, tg_xxxxxxz_xyzzz_0, \
                                         tg_xxxxxxz_xyzzz_1, tg_xxxxxxz_xzzz_1, tg_xxxxxxz_xzzzz_0, tg_xxxxxxz_xzzzz_1, \
                                         tg_xxxxxxz_yyyy_1, tg_xxxxxxz_yyyyy_0, tg_xxxxxxz_yyyyy_1, tg_xxxxxxz_yyyyz_0, \
                                         tg_xxxxxxz_yyyyz_1, tg_xxxxxxz_yyyz_1, tg_xxxxxxz_yyyzz_0, tg_xxxxxxz_yyyzz_1, \
                                         tg_xxxxxxz_yyzz_1, tg_xxxxxxz_yyzzz_0, tg_xxxxxxz_yyzzz_1, tg_xxxxxxz_yzzz_1, \
                                         tg_xxxxxxz_yzzzz_0, tg_xxxxxxz_yzzzz_1, tg_xxxxxxz_zzzz_1, tg_xxxxxxz_zzzzz_0, \
                                         tg_xxxxxxz_zzzzz_1, tg_xxxxxy_xxxxx_0, tg_xxxxxy_xxxxx_1, tg_xxxxxy_xxxxy_0, \
                                         tg_xxxxxy_xxxxy_1, tg_xxxxxy_xxxxz_0, tg_xxxxxy_xxxxz_1, tg_xxxxxy_xxxyy_0, \
                                         tg_xxxxxy_xxxyy_1, tg_xxxxxy_xxxyz_0, tg_xxxxxy_xxxyz_1, tg_xxxxxy_xxxzz_0, \
                                         tg_xxxxxy_xxxzz_1, tg_xxxxxy_xxyyy_0, tg_xxxxxy_xxyyy_1, tg_xxxxxy_xxyyz_0, \
                                         tg_xxxxxy_xxyyz_1, tg_xxxxxy_xxyzz_0, tg_xxxxxy_xxyzz_1, tg_xxxxxy_xxzzz_0, \
                                         tg_xxxxxy_xxzzz_1, tg_xxxxxy_xyyyy_0, tg_xxxxxy_xyyyy_1, tg_xxxxxy_xyyyz_0, \
                                         tg_xxxxxy_xyyyz_1, tg_xxxxxy_xyyzz_0, tg_xxxxxy_xyyzz_1, tg_xxxxxy_xyzzz_0, \
                                         tg_xxxxxy_xyzzz_1, tg_xxxxxy_xzzzz_0, tg_xxxxxy_xzzzz_1, tg_xxxxxy_yyyyy_0, \
                                         tg_xxxxxy_yyyyy_1, tg_xxxxxy_yyyyz_0, tg_xxxxxy_yyyyz_1, tg_xxxxxy_yyyzz_0, \
                                         tg_xxxxxy_yyyzz_1, tg_xxxxxy_yyzzz_0, tg_xxxxxy_yyzzz_1, tg_xxxxxy_yzzzz_0, \
                                         tg_xxxxxy_yzzzz_1, tg_xxxxxy_zzzzz_0, tg_xxxxxy_zzzzz_1, tg_xxxxxyy_xxxx_1, \
                                         tg_xxxxxyy_xxxxx_0, tg_xxxxxyy_xxxxx_1, tg_xxxxxyy_xxxxy_0, tg_xxxxxyy_xxxxy_1, \
                                         tg_xxxxxyy_xxxxz_0, tg_xxxxxyy_xxxxz_1, tg_xxxxxyy_xxxy_1, tg_xxxxxyy_xxxyy_0, \
                                         tg_xxxxxyy_xxxyy_1, tg_xxxxxyy_xxxyz_0, tg_xxxxxyy_xxxyz_1, tg_xxxxxyy_xxxz_1, \
                                         tg_xxxxxyy_xxxzz_0, tg_xxxxxyy_xxxzz_1, tg_xxxxxyy_xxyy_1, tg_xxxxxyy_xxyyy_0, \
                                         tg_xxxxxyy_xxyyy_1, tg_xxxxxyy_xxyyz_0, tg_xxxxxyy_xxyyz_1, tg_xxxxxyy_xxyz_1, \
                                         tg_xxxxxyy_xxyzz_0, tg_xxxxxyy_xxyzz_1, tg_xxxxxyy_xxzz_1, tg_xxxxxyy_xxzzz_0, \
                                         tg_xxxxxyy_xxzzz_1, tg_xxxxxyy_xyyy_1, tg_xxxxxyy_xyyyy_0, tg_xxxxxyy_xyyyy_1, \
                                         tg_xxxxxyy_xyyyz_0, tg_xxxxxyy_xyyyz_1, tg_xxxxxyy_xyyz_1, tg_xxxxxyy_xyyzz_0, \
                                         tg_xxxxxyy_xyyzz_1, tg_xxxxxyy_xyzz_1, tg_xxxxxyy_xyzzz_0, tg_xxxxxyy_xyzzz_1, \
                                         tg_xxxxxyy_xzzz_1, tg_xxxxxyy_xzzzz_0, tg_xxxxxyy_xzzzz_1, tg_xxxxxyy_yyyy_1, \
                                         tg_xxxxxyy_yyyyy_0, tg_xxxxxyy_yyyyy_1, tg_xxxxxyy_yyyyz_0, tg_xxxxxyy_yyyyz_1, \
                                         tg_xxxxxyy_yyyz_1, tg_xxxxxyy_yyyzz_0, tg_xxxxxyy_yyyzz_1, tg_xxxxxyy_yyzz_1, \
                                         tg_xxxxxyy_yyzzz_0, tg_xxxxxyy_yyzzz_1, tg_xxxxxyy_yzzz_1, tg_xxxxxyy_yzzzz_0, \
                                         tg_xxxxxyy_yzzzz_1, tg_xxxxxyy_zzzz_1, tg_xxxxxyy_zzzzz_0, tg_xxxxxyy_zzzzz_1, \
                                         tg_xxxxxyz_xxxx_1, tg_xxxxxyz_xxxxx_0, tg_xxxxxyz_xxxxx_1, tg_xxxxxyz_xxxxy_0, \
                                         tg_xxxxxyz_xxxxy_1, tg_xxxxxyz_xxxxz_0, tg_xxxxxyz_xxxxz_1, tg_xxxxxyz_xxxy_1, \
                                         tg_xxxxxyz_xxxyy_0, tg_xxxxxyz_xxxyy_1, tg_xxxxxyz_xxxyz_0, tg_xxxxxyz_xxxyz_1, \
                                         tg_xxxxxyz_xxxz_1, tg_xxxxxyz_xxxzz_0, tg_xxxxxyz_xxxzz_1, tg_xxxxxyz_xxyy_1, \
                                         tg_xxxxxyz_xxyyy_0, tg_xxxxxyz_xxyyy_1, tg_xxxxxyz_xxyyz_0, tg_xxxxxyz_xxyyz_1, \
                                         tg_xxxxxyz_xxyz_1, tg_xxxxxyz_xxyzz_0, tg_xxxxxyz_xxyzz_1, tg_xxxxxyz_xxzz_1, \
                                         tg_xxxxxyz_xxzzz_0, tg_xxxxxyz_xxzzz_1, tg_xxxxxyz_xyyy_1, tg_xxxxxyz_xyyyy_0, \
                                         tg_xxxxxyz_xyyyy_1, tg_xxxxxyz_xyyz_1, tg_xxxxxyz_xyzz_1, tg_xxxxxyz_xzzz_1, \
                                         tg_xxxxxyz_yyyy_1, tg_xxxxxz_xxxxx_0, tg_xxxxxz_xxxxx_1, tg_xxxxxz_xxxxy_0, \
                                         tg_xxxxxz_xxxxy_1, tg_xxxxxz_xxxxz_0, tg_xxxxxz_xxxxz_1, tg_xxxxxz_xxxyy_0, \
                                         tg_xxxxxz_xxxyy_1, tg_xxxxxz_xxxyz_0, tg_xxxxxz_xxxyz_1, tg_xxxxxz_xxxzz_0, \
                                         tg_xxxxxz_xxxzz_1, tg_xxxxxz_xxyyy_0, tg_xxxxxz_xxyyy_1, tg_xxxxxz_xxyyz_0, \
                                         tg_xxxxxz_xxyyz_1, tg_xxxxxz_xxyzz_0, tg_xxxxxz_xxyzz_1, tg_xxxxxz_xxzzz_0, \
                                         tg_xxxxxz_xxzzz_1, tg_xxxxxz_xyyyy_0, tg_xxxxxz_xyyyy_1, tg_xxxxxz_xyyyz_0, \
                                         tg_xxxxxz_xyyyz_1, tg_xxxxxz_xyyzz_0, tg_xxxxxz_xyyzz_1, tg_xxxxxz_xyzzz_0, \
                                         tg_xxxxxz_xyzzz_1, tg_xxxxxz_xzzzz_0, tg_xxxxxz_xzzzz_1, tg_xxxxxz_yyyyy_0, \
                                         tg_xxxxxz_yyyyy_1, tg_xxxxxz_yyyyz_0, tg_xxxxxz_yyyyz_1, tg_xxxxxz_yyyzz_0, \
                                         tg_xxxxxz_yyyzz_1, tg_xxxxxz_yyzzz_0, tg_xxxxxz_yyzzz_1, tg_xxxxxz_yzzzz_0, \
                                         tg_xxxxxz_yzzzz_1, tg_xxxxxz_zzzzz_0, tg_xxxxxz_zzzzz_1, tg_xxxxyy_xxxxx_0, \
                                         tg_xxxxyy_xxxxx_1, tg_xxxxyy_xxxxy_0, tg_xxxxyy_xxxxy_1, tg_xxxxyy_xxxxz_0, \
                                         tg_xxxxyy_xxxxz_1, tg_xxxxyy_xxxyy_0, tg_xxxxyy_xxxyy_1, tg_xxxxyy_xxxyz_0, \
                                         tg_xxxxyy_xxxyz_1, tg_xxxxyy_xxxzz_0, tg_xxxxyy_xxxzz_1, tg_xxxxyy_xxyyy_0, \
                                         tg_xxxxyy_xxyyy_1, tg_xxxxyy_xxyyz_0, tg_xxxxyy_xxyyz_1, tg_xxxxyy_xxyzz_0, \
                                         tg_xxxxyy_xxyzz_1, tg_xxxxyy_xxzzz_0, tg_xxxxyy_xxzzz_1, tg_xxxxyy_xyyyy_0, \
                                         tg_xxxxyy_xyyyy_1, tg_xxxxyy_xyyyz_0, tg_xxxxyy_xyyyz_1, tg_xxxxyy_xyyzz_0, \
                                         tg_xxxxyy_xyyzz_1, tg_xxxxyy_xyzzz_0, tg_xxxxyy_xyzzz_1, tg_xxxxyy_xzzzz_0, \
                                         tg_xxxxyy_xzzzz_1, tg_xxxxyy_yyyyy_0, tg_xxxxyy_yyyyy_1, tg_xxxxyy_yyyyz_0, \
                                         tg_xxxxyy_yyyyz_1, tg_xxxxyy_yyyzz_0, tg_xxxxyy_yyyzz_1, tg_xxxxyy_yyzzz_0, \
                                         tg_xxxxyy_yyzzz_1, tg_xxxxyy_yzzzz_0, tg_xxxxyy_yzzzz_1, tg_xxxxyy_zzzzz_0, \
                                         tg_xxxxyy_zzzzz_1, tg_xxxxyz_xxxxx_0, tg_xxxxyz_xxxxx_1, tg_xxxxyz_xxxxy_0, \
                                         tg_xxxxyz_xxxxy_1, tg_xxxxyz_xxxxz_0, tg_xxxxyz_xxxxz_1, tg_xxxxyz_xxxyy_0, \
                                         tg_xxxxyz_xxxyy_1, tg_xxxxyz_xxxyz_0, tg_xxxxyz_xxxyz_1, tg_xxxxyz_xxxzz_0, \
                                         tg_xxxxyz_xxxzz_1, tg_xxxxyz_xxyyy_0, tg_xxxxyz_xxyyy_1, tg_xxxxyz_xxyyz_0, \
                                         tg_xxxxyz_xxyyz_1, tg_xxxxyz_xxyzz_0, tg_xxxxyz_xxyzz_1, tg_xxxxyz_xxzzz_0, \
                                         tg_xxxxyz_xxzzz_1, tg_xxxxyz_xyyyy_0, tg_xxxxyz_xyyyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxxxx_xxxxx_0[j] = pb_x * tg_xxxxxxx_xxxxx_0[j] + fr * tg_xxxxxxx_xxxxx_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxxxx_0[j] - tg_xxxxxx_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxxx_xxxx_1[j];

                    tg_xxxxxxxx_xxxxy_0[j] = pb_x * tg_xxxxxxx_xxxxy_0[j] + fr * tg_xxxxxxx_xxxxy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxxxy_0[j] - tg_xxxxxx_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxxx_xxxy_1[j];

                    tg_xxxxxxxx_xxxxz_0[j] = pb_x * tg_xxxxxxx_xxxxz_0[j] + fr * tg_xxxxxxx_xxxxz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxxxz_0[j] - tg_xxxxxx_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxxx_xxxz_1[j];

                    tg_xxxxxxxx_xxxyy_0[j] = pb_x * tg_xxxxxxx_xxxyy_0[j] + fr * tg_xxxxxxx_xxxyy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxxyy_0[j] - tg_xxxxxx_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxx_xxyy_1[j];

                    tg_xxxxxxxx_xxxyz_0[j] = pb_x * tg_xxxxxxx_xxxyz_0[j] + fr * tg_xxxxxxx_xxxyz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxxyz_0[j] - tg_xxxxxx_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxx_xxyz_1[j];

                    tg_xxxxxxxx_xxxzz_0[j] = pb_x * tg_xxxxxxx_xxxzz_0[j] + fr * tg_xxxxxxx_xxxzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxxzz_0[j] - tg_xxxxxx_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxx_xxzz_1[j];

                    tg_xxxxxxxx_xxyyy_0[j] = pb_x * tg_xxxxxxx_xxyyy_0[j] + fr * tg_xxxxxxx_xxyyy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxyyy_0[j] - tg_xxxxxx_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxx_xyyy_1[j];

                    tg_xxxxxxxx_xxyyz_0[j] = pb_x * tg_xxxxxxx_xxyyz_0[j] + fr * tg_xxxxxxx_xxyyz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxyyz_0[j] - tg_xxxxxx_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxx_xyyz_1[j];

                    tg_xxxxxxxx_xxyzz_0[j] = pb_x * tg_xxxxxxx_xxyzz_0[j] + fr * tg_xxxxxxx_xxyzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxyzz_0[j] - tg_xxxxxx_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxx_xyzz_1[j];

                    tg_xxxxxxxx_xxzzz_0[j] = pb_x * tg_xxxxxxx_xxzzz_0[j] + fr * tg_xxxxxxx_xxzzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxzzz_0[j] - tg_xxxxxx_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxx_xzzz_1[j];

                    tg_xxxxxxxx_xyyyy_0[j] = pb_x * tg_xxxxxxx_xyyyy_0[j] + fr * tg_xxxxxxx_xyyyy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xyyyy_0[j] - tg_xxxxxx_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_yyyy_1[j];

                    tg_xxxxxxxx_xyyyz_0[j] = pb_x * tg_xxxxxxx_xyyyz_0[j] + fr * tg_xxxxxxx_xyyyz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xyyyz_0[j] - tg_xxxxxx_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_yyyz_1[j];

                    tg_xxxxxxxx_xyyzz_0[j] = pb_x * tg_xxxxxxx_xyyzz_0[j] + fr * tg_xxxxxxx_xyyzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xyyzz_0[j] - tg_xxxxxx_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_yyzz_1[j];

                    tg_xxxxxxxx_xyzzz_0[j] = pb_x * tg_xxxxxxx_xyzzz_0[j] + fr * tg_xxxxxxx_xyzzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xyzzz_0[j] - tg_xxxxxx_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_yzzz_1[j];

                    tg_xxxxxxxx_xzzzz_0[j] = pb_x * tg_xxxxxxx_xzzzz_0[j] + fr * tg_xxxxxxx_xzzzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xzzzz_0[j] - tg_xxxxxx_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_zzzz_1[j];

                    tg_xxxxxxxx_yyyyy_0[j] = pb_x * tg_xxxxxxx_yyyyy_0[j] + fr * tg_xxxxxxx_yyyyy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yyyyy_0[j] - tg_xxxxxx_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxxx_yyyyz_0[j] = pb_x * tg_xxxxxxx_yyyyz_0[j] + fr * tg_xxxxxxx_yyyyz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yyyyz_0[j] - tg_xxxxxx_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxxx_yyyzz_0[j] = pb_x * tg_xxxxxxx_yyyzz_0[j] + fr * tg_xxxxxxx_yyyzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yyyzz_0[j] - tg_xxxxxx_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxxx_yyzzz_0[j] = pb_x * tg_xxxxxxx_yyzzz_0[j] + fr * tg_xxxxxxx_yyzzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yyzzz_0[j] - tg_xxxxxx_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxxx_yzzzz_0[j] = pb_x * tg_xxxxxxx_yzzzz_0[j] + fr * tg_xxxxxxx_yzzzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yzzzz_0[j] - tg_xxxxxx_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxxx_zzzzz_0[j] = pb_x * tg_xxxxxxx_zzzzz_0[j] + fr * tg_xxxxxxx_zzzzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_zzzzz_0[j] - tg_xxxxxx_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxxxy_xxxxx_0[j] = pb_x * tg_xxxxxxy_xxxxx_0[j] + fr * tg_xxxxxxy_xxxxx_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxxxx_0[j] - tg_xxxxxy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxxy_xxxx_1[j];

                    tg_xxxxxxxy_xxxxy_0[j] = pb_x * tg_xxxxxxy_xxxxy_0[j] + fr * tg_xxxxxxy_xxxxy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxxxy_0[j] - tg_xxxxxy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxxy_xxxy_1[j];

                    tg_xxxxxxxy_xxxxz_0[j] = pb_x * tg_xxxxxxy_xxxxz_0[j] + fr * tg_xxxxxxy_xxxxz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxxxz_0[j] - tg_xxxxxy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxxy_xxxz_1[j];

                    tg_xxxxxxxy_xxxyy_0[j] = pb_x * tg_xxxxxxy_xxxyy_0[j] + fr * tg_xxxxxxy_xxxyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxxyy_0[j] - tg_xxxxxy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxy_xxyy_1[j];

                    tg_xxxxxxxy_xxxyz_0[j] = pb_x * tg_xxxxxxy_xxxyz_0[j] + fr * tg_xxxxxxy_xxxyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxxyz_0[j] - tg_xxxxxy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxy_xxyz_1[j];

                    tg_xxxxxxxy_xxxzz_0[j] = pb_x * tg_xxxxxxy_xxxzz_0[j] + fr * tg_xxxxxxy_xxxzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxxzz_0[j] - tg_xxxxxy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxy_xxzz_1[j];

                    tg_xxxxxxxy_xxyyy_0[j] = pb_x * tg_xxxxxxy_xxyyy_0[j] + fr * tg_xxxxxxy_xxyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxyyy_0[j] - tg_xxxxxy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxy_xyyy_1[j];

                    tg_xxxxxxxy_xxyyz_0[j] = pb_x * tg_xxxxxxy_xxyyz_0[j] + fr * tg_xxxxxxy_xxyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxyyz_0[j] - tg_xxxxxy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxy_xyyz_1[j];

                    tg_xxxxxxxy_xxyzz_0[j] = pb_x * tg_xxxxxxy_xxyzz_0[j] + fr * tg_xxxxxxy_xxyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxyzz_0[j] - tg_xxxxxy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxy_xyzz_1[j];

                    tg_xxxxxxxy_xxzzz_0[j] = pb_x * tg_xxxxxxy_xxzzz_0[j] + fr * tg_xxxxxxy_xxzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxzzz_0[j] - tg_xxxxxy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxy_xzzz_1[j];

                    tg_xxxxxxxy_xyyyy_0[j] = pb_x * tg_xxxxxxy_xyyyy_0[j] + fr * tg_xxxxxxy_xyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xyyyy_0[j] - tg_xxxxxy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_yyyy_1[j];

                    tg_xxxxxxxy_xyyyz_0[j] = pb_x * tg_xxxxxxy_xyyyz_0[j] + fr * tg_xxxxxxy_xyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xyyyz_0[j] - tg_xxxxxy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_yyyz_1[j];

                    tg_xxxxxxxy_xyyzz_0[j] = pb_x * tg_xxxxxxy_xyyzz_0[j] + fr * tg_xxxxxxy_xyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xyyzz_0[j] - tg_xxxxxy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_yyzz_1[j];

                    tg_xxxxxxxy_xyzzz_0[j] = pb_x * tg_xxxxxxy_xyzzz_0[j] + fr * tg_xxxxxxy_xyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xyzzz_0[j] - tg_xxxxxy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_yzzz_1[j];

                    tg_xxxxxxxy_xzzzz_0[j] = pb_x * tg_xxxxxxy_xzzzz_0[j] + fr * tg_xxxxxxy_xzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xzzzz_0[j] - tg_xxxxxy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_zzzz_1[j];

                    tg_xxxxxxxy_yyyyy_0[j] = pb_x * tg_xxxxxxy_yyyyy_0[j] + fr * tg_xxxxxxy_yyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yyyyy_0[j] - tg_xxxxxy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxxy_yyyyz_0[j] = pb_x * tg_xxxxxxy_yyyyz_0[j] + fr * tg_xxxxxxy_yyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yyyyz_0[j] - tg_xxxxxy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxxy_yyyzz_0[j] = pb_x * tg_xxxxxxy_yyyzz_0[j] + fr * tg_xxxxxxy_yyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yyyzz_0[j] - tg_xxxxxy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxxy_yyzzz_0[j] = pb_x * tg_xxxxxxy_yyzzz_0[j] + fr * tg_xxxxxxy_yyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yyzzz_0[j] - tg_xxxxxy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxxy_yzzzz_0[j] = pb_x * tg_xxxxxxy_yzzzz_0[j] + fr * tg_xxxxxxy_yzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yzzzz_0[j] - tg_xxxxxy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxxy_zzzzz_0[j] = pb_x * tg_xxxxxxy_zzzzz_0[j] + fr * tg_xxxxxxy_zzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_zzzzz_0[j] - tg_xxxxxy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxxxz_xxxxx_0[j] = pb_x * tg_xxxxxxz_xxxxx_0[j] + fr * tg_xxxxxxz_xxxxx_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxxxx_0[j] - tg_xxxxxz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxxz_xxxx_1[j];

                    tg_xxxxxxxz_xxxxy_0[j] = pb_x * tg_xxxxxxz_xxxxy_0[j] + fr * tg_xxxxxxz_xxxxy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxxxy_0[j] - tg_xxxxxz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxxz_xxxy_1[j];

                    tg_xxxxxxxz_xxxxz_0[j] = pb_x * tg_xxxxxxz_xxxxz_0[j] + fr * tg_xxxxxxz_xxxxz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxxxz_0[j] - tg_xxxxxz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxxz_xxxz_1[j];

                    tg_xxxxxxxz_xxxyy_0[j] = pb_x * tg_xxxxxxz_xxxyy_0[j] + fr * tg_xxxxxxz_xxxyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxxyy_0[j] - tg_xxxxxz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxz_xxyy_1[j];

                    tg_xxxxxxxz_xxxyz_0[j] = pb_x * tg_xxxxxxz_xxxyz_0[j] + fr * tg_xxxxxxz_xxxyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxxyz_0[j] - tg_xxxxxz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxz_xxyz_1[j];

                    tg_xxxxxxxz_xxxzz_0[j] = pb_x * tg_xxxxxxz_xxxzz_0[j] + fr * tg_xxxxxxz_xxxzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxxzz_0[j] - tg_xxxxxz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxz_xxzz_1[j];

                    tg_xxxxxxxz_xxyyy_0[j] = pb_x * tg_xxxxxxz_xxyyy_0[j] + fr * tg_xxxxxxz_xxyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxyyy_0[j] - tg_xxxxxz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxz_xyyy_1[j];

                    tg_xxxxxxxz_xxyyz_0[j] = pb_x * tg_xxxxxxz_xxyyz_0[j] + fr * tg_xxxxxxz_xxyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxyyz_0[j] - tg_xxxxxz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxz_xyyz_1[j];

                    tg_xxxxxxxz_xxyzz_0[j] = pb_x * tg_xxxxxxz_xxyzz_0[j] + fr * tg_xxxxxxz_xxyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxyzz_0[j] - tg_xxxxxz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxz_xyzz_1[j];

                    tg_xxxxxxxz_xxzzz_0[j] = pb_x * tg_xxxxxxz_xxzzz_0[j] + fr * tg_xxxxxxz_xxzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxzzz_0[j] - tg_xxxxxz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxz_xzzz_1[j];

                    tg_xxxxxxxz_xyyyy_0[j] = pb_x * tg_xxxxxxz_xyyyy_0[j] + fr * tg_xxxxxxz_xyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xyyyy_0[j] - tg_xxxxxz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_yyyy_1[j];

                    tg_xxxxxxxz_xyyyz_0[j] = pb_x * tg_xxxxxxz_xyyyz_0[j] + fr * tg_xxxxxxz_xyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xyyyz_0[j] - tg_xxxxxz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_yyyz_1[j];

                    tg_xxxxxxxz_xyyzz_0[j] = pb_x * tg_xxxxxxz_xyyzz_0[j] + fr * tg_xxxxxxz_xyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xyyzz_0[j] - tg_xxxxxz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_yyzz_1[j];

                    tg_xxxxxxxz_xyzzz_0[j] = pb_x * tg_xxxxxxz_xyzzz_0[j] + fr * tg_xxxxxxz_xyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xyzzz_0[j] - tg_xxxxxz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_yzzz_1[j];

                    tg_xxxxxxxz_xzzzz_0[j] = pb_x * tg_xxxxxxz_xzzzz_0[j] + fr * tg_xxxxxxz_xzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xzzzz_0[j] - tg_xxxxxz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_zzzz_1[j];

                    tg_xxxxxxxz_yyyyy_0[j] = pb_x * tg_xxxxxxz_yyyyy_0[j] + fr * tg_xxxxxxz_yyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yyyyy_0[j] - tg_xxxxxz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxxz_yyyyz_0[j] = pb_x * tg_xxxxxxz_yyyyz_0[j] + fr * tg_xxxxxxz_yyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yyyyz_0[j] - tg_xxxxxz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxxz_yyyzz_0[j] = pb_x * tg_xxxxxxz_yyyzz_0[j] + fr * tg_xxxxxxz_yyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yyyzz_0[j] - tg_xxxxxz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxxz_yyzzz_0[j] = pb_x * tg_xxxxxxz_yyzzz_0[j] + fr * tg_xxxxxxz_yyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yyzzz_0[j] - tg_xxxxxz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxxz_yzzzz_0[j] = pb_x * tg_xxxxxxz_yzzzz_0[j] + fr * tg_xxxxxxz_yzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yzzzz_0[j] - tg_xxxxxz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxxz_zzzzz_0[j] = pb_x * tg_xxxxxxz_zzzzz_0[j] + fr * tg_xxxxxxz_zzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_zzzzz_0[j] - tg_xxxxxz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxxyy_xxxxx_0[j] = pb_x * tg_xxxxxyy_xxxxx_0[j] + fr * tg_xxxxxyy_xxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxxxx_0[j] - tg_xxxxyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxyy_xxxx_1[j];

                    tg_xxxxxxyy_xxxxy_0[j] = pb_x * tg_xxxxxyy_xxxxy_0[j] + fr * tg_xxxxxyy_xxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxxxy_0[j] - tg_xxxxyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxyy_xxxy_1[j];

                    tg_xxxxxxyy_xxxxz_0[j] = pb_x * tg_xxxxxyy_xxxxz_0[j] + fr * tg_xxxxxyy_xxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxxxz_0[j] - tg_xxxxyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxyy_xxxz_1[j];

                    tg_xxxxxxyy_xxxyy_0[j] = pb_x * tg_xxxxxyy_xxxyy_0[j] + fr * tg_xxxxxyy_xxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxxyy_0[j] - tg_xxxxyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyy_xxyy_1[j];

                    tg_xxxxxxyy_xxxyz_0[j] = pb_x * tg_xxxxxyy_xxxyz_0[j] + fr * tg_xxxxxyy_xxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxxyz_0[j] - tg_xxxxyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyy_xxyz_1[j];

                    tg_xxxxxxyy_xxxzz_0[j] = pb_x * tg_xxxxxyy_xxxzz_0[j] + fr * tg_xxxxxyy_xxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxxzz_0[j] - tg_xxxxyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyy_xxzz_1[j];

                    tg_xxxxxxyy_xxyyy_0[j] = pb_x * tg_xxxxxyy_xxyyy_0[j] + fr * tg_xxxxxyy_xxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxyyy_0[j] - tg_xxxxyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyy_xyyy_1[j];

                    tg_xxxxxxyy_xxyyz_0[j] = pb_x * tg_xxxxxyy_xxyyz_0[j] + fr * tg_xxxxxyy_xxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxyyz_0[j] - tg_xxxxyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyy_xyyz_1[j];

                    tg_xxxxxxyy_xxyzz_0[j] = pb_x * tg_xxxxxyy_xxyzz_0[j] + fr * tg_xxxxxyy_xxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxyzz_0[j] - tg_xxxxyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyy_xyzz_1[j];

                    tg_xxxxxxyy_xxzzz_0[j] = pb_x * tg_xxxxxyy_xxzzz_0[j] + fr * tg_xxxxxyy_xxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxzzz_0[j] - tg_xxxxyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyy_xzzz_1[j];

                    tg_xxxxxxyy_xyyyy_0[j] = pb_x * tg_xxxxxyy_xyyyy_0[j] + fr * tg_xxxxxyy_xyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xyyyy_0[j] - tg_xxxxyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_yyyy_1[j];

                    tg_xxxxxxyy_xyyyz_0[j] = pb_x * tg_xxxxxyy_xyyyz_0[j] + fr * tg_xxxxxyy_xyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xyyyz_0[j] - tg_xxxxyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_yyyz_1[j];

                    tg_xxxxxxyy_xyyzz_0[j] = pb_x * tg_xxxxxyy_xyyzz_0[j] + fr * tg_xxxxxyy_xyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xyyzz_0[j] - tg_xxxxyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_yyzz_1[j];

                    tg_xxxxxxyy_xyzzz_0[j] = pb_x * tg_xxxxxyy_xyzzz_0[j] + fr * tg_xxxxxyy_xyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xyzzz_0[j] - tg_xxxxyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_yzzz_1[j];

                    tg_xxxxxxyy_xzzzz_0[j] = pb_x * tg_xxxxxyy_xzzzz_0[j] + fr * tg_xxxxxyy_xzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xzzzz_0[j] - tg_xxxxyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_zzzz_1[j];

                    tg_xxxxxxyy_yyyyy_0[j] = pb_x * tg_xxxxxyy_yyyyy_0[j] + fr * tg_xxxxxyy_yyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yyyyy_0[j] - tg_xxxxyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxyy_yyyyz_0[j] = pb_x * tg_xxxxxyy_yyyyz_0[j] + fr * tg_xxxxxyy_yyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yyyyz_0[j] - tg_xxxxyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxyy_yyyzz_0[j] = pb_x * tg_xxxxxyy_yyyzz_0[j] + fr * tg_xxxxxyy_yyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yyyzz_0[j] - tg_xxxxyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxyy_yyzzz_0[j] = pb_x * tg_xxxxxyy_yyzzz_0[j] + fr * tg_xxxxxyy_yyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yyzzz_0[j] - tg_xxxxyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxyy_yzzzz_0[j] = pb_x * tg_xxxxxyy_yzzzz_0[j] + fr * tg_xxxxxyy_yzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yzzzz_0[j] - tg_xxxxyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxyy_zzzzz_0[j] = pb_x * tg_xxxxxyy_zzzzz_0[j] + fr * tg_xxxxxyy_zzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_zzzzz_0[j] - tg_xxxxyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxxyz_xxxxx_0[j] = pb_x * tg_xxxxxyz_xxxxx_0[j] + fr * tg_xxxxxyz_xxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxxxx_0[j] - tg_xxxxyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxyz_xxxx_1[j];

                    tg_xxxxxxyz_xxxxy_0[j] = pb_x * tg_xxxxxyz_xxxxy_0[j] + fr * tg_xxxxxyz_xxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxxxy_0[j] - tg_xxxxyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxyz_xxxy_1[j];

                    tg_xxxxxxyz_xxxxz_0[j] = pb_x * tg_xxxxxyz_xxxxz_0[j] + fr * tg_xxxxxyz_xxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxxxz_0[j] - tg_xxxxyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxyz_xxxz_1[j];

                    tg_xxxxxxyz_xxxyy_0[j] = pb_x * tg_xxxxxyz_xxxyy_0[j] + fr * tg_xxxxxyz_xxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxxyy_0[j] - tg_xxxxyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyz_xxyy_1[j];

                    tg_xxxxxxyz_xxxyz_0[j] = pb_x * tg_xxxxxyz_xxxyz_0[j] + fr * tg_xxxxxyz_xxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxxyz_0[j] - tg_xxxxyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyz_xxyz_1[j];

                    tg_xxxxxxyz_xxxzz_0[j] = pb_x * tg_xxxxxyz_xxxzz_0[j] + fr * tg_xxxxxyz_xxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxxzz_0[j] - tg_xxxxyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyz_xxzz_1[j];

                    tg_xxxxxxyz_xxyyy_0[j] = pb_x * tg_xxxxxyz_xxyyy_0[j] + fr * tg_xxxxxyz_xxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxyyy_0[j] - tg_xxxxyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyz_xyyy_1[j];

                    tg_xxxxxxyz_xxyyz_0[j] = pb_x * tg_xxxxxyz_xxyyz_0[j] + fr * tg_xxxxxyz_xxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxyyz_0[j] - tg_xxxxyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyz_xyyz_1[j];

                    tg_xxxxxxyz_xxyzz_0[j] = pb_x * tg_xxxxxyz_xxyzz_0[j] + fr * tg_xxxxxyz_xxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxyzz_0[j] - tg_xxxxyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyz_xyzz_1[j];

                    tg_xxxxxxyz_xxzzz_0[j] = pb_x * tg_xxxxxyz_xxzzz_0[j] + fr * tg_xxxxxyz_xxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxzzz_0[j] - tg_xxxxyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyz_xzzz_1[j];

                    tg_xxxxxxyz_xyyyy_0[j] = pb_x * tg_xxxxxyz_xyyyy_0[j] + fr * tg_xxxxxyz_xyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xyyyy_0[j] - tg_xxxxyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_yyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_95_190(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxxxxyz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 95); 

                auto tg_xxxxxyz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 96); 

                auto tg_xxxxxyz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 97); 

                auto tg_xxxxxyz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 98); 

                auto tg_xxxxxyz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 99); 

                auto tg_xxxxxyz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 100); 

                auto tg_xxxxxyz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 101); 

                auto tg_xxxxxyz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 102); 

                auto tg_xxxxxyz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 103); 

                auto tg_xxxxxyz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 104); 

                auto tg_xxxxxzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 105); 

                auto tg_xxxxxzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 106); 

                auto tg_xxxxxzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 107); 

                auto tg_xxxxxzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 108); 

                auto tg_xxxxxzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 109); 

                auto tg_xxxxxzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 110); 

                auto tg_xxxxxzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 111); 

                auto tg_xxxxxzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 112); 

                auto tg_xxxxxzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 113); 

                auto tg_xxxxxzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 114); 

                auto tg_xxxxxzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 115); 

                auto tg_xxxxxzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 116); 

                auto tg_xxxxxzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 117); 

                auto tg_xxxxxzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 118); 

                auto tg_xxxxxzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 119); 

                auto tg_xxxxxzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 120); 

                auto tg_xxxxxzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 121); 

                auto tg_xxxxxzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 122); 

                auto tg_xxxxxzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 123); 

                auto tg_xxxxxzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 124); 

                auto tg_xxxxxzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 125); 

                auto tg_xxxxyyy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 126); 

                auto tg_xxxxyyy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 127); 

                auto tg_xxxxyyy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 128); 

                auto tg_xxxxyyy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 129); 

                auto tg_xxxxyyy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 130); 

                auto tg_xxxxyyy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 131); 

                auto tg_xxxxyyy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 132); 

                auto tg_xxxxyyy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 133); 

                auto tg_xxxxyyy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 134); 

                auto tg_xxxxyyy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 135); 

                auto tg_xxxxyyy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 136); 

                auto tg_xxxxyyy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 137); 

                auto tg_xxxxyyy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 138); 

                auto tg_xxxxyyy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 139); 

                auto tg_xxxxyyy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 140); 

                auto tg_xxxxyyy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 141); 

                auto tg_xxxxyyy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 142); 

                auto tg_xxxxyyy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 143); 

                auto tg_xxxxyyy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 144); 

                auto tg_xxxxyyy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 145); 

                auto tg_xxxxyyy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 146); 

                auto tg_xxxxyyz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 147); 

                auto tg_xxxxyyz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 148); 

                auto tg_xxxxyyz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 149); 

                auto tg_xxxxyyz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 150); 

                auto tg_xxxxyyz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 151); 

                auto tg_xxxxyyz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 152); 

                auto tg_xxxxyyz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 153); 

                auto tg_xxxxyyz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 154); 

                auto tg_xxxxyyz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 155); 

                auto tg_xxxxyyz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 156); 

                auto tg_xxxxyyz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 157); 

                auto tg_xxxxyyz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 158); 

                auto tg_xxxxyyz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 159); 

                auto tg_xxxxyyz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 160); 

                auto tg_xxxxyyz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 161); 

                auto tg_xxxxyyz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 162); 

                auto tg_xxxxyyz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 163); 

                auto tg_xxxxyyz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 164); 

                auto tg_xxxxyyz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 165); 

                auto tg_xxxxyyz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 166); 

                auto tg_xxxxyyz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 167); 

                auto tg_xxxxyzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 168); 

                auto tg_xxxxyzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 169); 

                auto tg_xxxxyzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 170); 

                auto tg_xxxxyzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 171); 

                auto tg_xxxxyzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 172); 

                auto tg_xxxxyzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 173); 

                auto tg_xxxxyzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 174); 

                auto tg_xxxxyzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 175); 

                auto tg_xxxxyzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 176); 

                auto tg_xxxxyzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 177); 

                auto tg_xxxxyzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 178); 

                auto tg_xxxxyzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 179); 

                auto tg_xxxxyzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 180); 

                auto tg_xxxxyzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 181); 

                auto tg_xxxxyzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 182); 

                auto tg_xxxxyzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 183); 

                auto tg_xxxxyzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 184); 

                auto tg_xxxxyzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 185); 

                auto tg_xxxxyzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 186); 

                auto tg_xxxxyzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 187); 

                auto tg_xxxxyzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 188); 

                auto tg_xxxxzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 189); 

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

                auto tg_xxxxxyz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 71); 

                auto tg_xxxxxyz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 72); 

                auto tg_xxxxxyz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 73); 

                auto tg_xxxxxyz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 74); 

                auto tg_xxxxxzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 75); 

                auto tg_xxxxxzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 76); 

                auto tg_xxxxxzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 77); 

                auto tg_xxxxxzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 78); 

                auto tg_xxxxxzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 79); 

                auto tg_xxxxxzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 80); 

                auto tg_xxxxxzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 81); 

                auto tg_xxxxxzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 82); 

                auto tg_xxxxxzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 83); 

                auto tg_xxxxxzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 84); 

                auto tg_xxxxxzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 85); 

                auto tg_xxxxxzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 86); 

                auto tg_xxxxxzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 87); 

                auto tg_xxxxxzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 88); 

                auto tg_xxxxxzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 89); 

                auto tg_xxxxyyy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 90); 

                auto tg_xxxxyyy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 91); 

                auto tg_xxxxyyy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 92); 

                auto tg_xxxxyyy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 93); 

                auto tg_xxxxyyy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 94); 

                auto tg_xxxxyyy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 95); 

                auto tg_xxxxyyy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 96); 

                auto tg_xxxxyyy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 97); 

                auto tg_xxxxyyy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 98); 

                auto tg_xxxxyyy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 99); 

                auto tg_xxxxyyy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 100); 

                auto tg_xxxxyyy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 101); 

                auto tg_xxxxyyy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 102); 

                auto tg_xxxxyyy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 103); 

                auto tg_xxxxyyy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 104); 

                auto tg_xxxxyyz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 105); 

                auto tg_xxxxyyz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 106); 

                auto tg_xxxxyyz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 107); 

                auto tg_xxxxyyz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 108); 

                auto tg_xxxxyyz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 109); 

                auto tg_xxxxyyz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 110); 

                auto tg_xxxxyyz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 111); 

                auto tg_xxxxyyz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 112); 

                auto tg_xxxxyyz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 113); 

                auto tg_xxxxyyz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 114); 

                auto tg_xxxxyyz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 115); 

                auto tg_xxxxyyz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 116); 

                auto tg_xxxxyyz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 117); 

                auto tg_xxxxyyz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 118); 

                auto tg_xxxxyyz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 119); 

                auto tg_xxxxyzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 120); 

                auto tg_xxxxyzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 121); 

                auto tg_xxxxyzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 122); 

                auto tg_xxxxyzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 123); 

                auto tg_xxxxyzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 124); 

                auto tg_xxxxyzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 125); 

                auto tg_xxxxyzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 126); 

                auto tg_xxxxyzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 127); 

                auto tg_xxxxyzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 128); 

                auto tg_xxxxyzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 129); 

                auto tg_xxxxyzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 130); 

                auto tg_xxxxyzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 131); 

                auto tg_xxxxyzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 132); 

                auto tg_xxxxyzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 133); 

                auto tg_xxxxyzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 134); 

                auto tg_xxxxzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 135); 

                // set up pointers to integrals

                auto tg_xxxxxxyz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 95); 

                auto tg_xxxxxxyz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 96); 

                auto tg_xxxxxxyz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 97); 

                auto tg_xxxxxxyz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 98); 

                auto tg_xxxxxxyz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 99); 

                auto tg_xxxxxxyz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 100); 

                auto tg_xxxxxxyz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 101); 

                auto tg_xxxxxxyz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 102); 

                auto tg_xxxxxxyz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 103); 

                auto tg_xxxxxxyz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 104); 

                auto tg_xxxxxxzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 105); 

                auto tg_xxxxxxzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 106); 

                auto tg_xxxxxxzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 107); 

                auto tg_xxxxxxzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 108); 

                auto tg_xxxxxxzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 109); 

                auto tg_xxxxxxzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 110); 

                auto tg_xxxxxxzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 111); 

                auto tg_xxxxxxzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 112); 

                auto tg_xxxxxxzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 113); 

                auto tg_xxxxxxzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 114); 

                auto tg_xxxxxxzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 115); 

                auto tg_xxxxxxzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 116); 

                auto tg_xxxxxxzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 117); 

                auto tg_xxxxxxzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 118); 

                auto tg_xxxxxxzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 119); 

                auto tg_xxxxxxzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 120); 

                auto tg_xxxxxxzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 121); 

                auto tg_xxxxxxzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 122); 

                auto tg_xxxxxxzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 123); 

                auto tg_xxxxxxzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 124); 

                auto tg_xxxxxxzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 125); 

                auto tg_xxxxxyyy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 126); 

                auto tg_xxxxxyyy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 127); 

                auto tg_xxxxxyyy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 128); 

                auto tg_xxxxxyyy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 129); 

                auto tg_xxxxxyyy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 130); 

                auto tg_xxxxxyyy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 131); 

                auto tg_xxxxxyyy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 132); 

                auto tg_xxxxxyyy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 133); 

                auto tg_xxxxxyyy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 134); 

                auto tg_xxxxxyyy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 135); 

                auto tg_xxxxxyyy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 136); 

                auto tg_xxxxxyyy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 137); 

                auto tg_xxxxxyyy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 138); 

                auto tg_xxxxxyyy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 139); 

                auto tg_xxxxxyyy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 140); 

                auto tg_xxxxxyyy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 141); 

                auto tg_xxxxxyyy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 142); 

                auto tg_xxxxxyyy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 143); 

                auto tg_xxxxxyyy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 144); 

                auto tg_xxxxxyyy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 145); 

                auto tg_xxxxxyyy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 146); 

                auto tg_xxxxxyyz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 147); 

                auto tg_xxxxxyyz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 148); 

                auto tg_xxxxxyyz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 149); 

                auto tg_xxxxxyyz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 150); 

                auto tg_xxxxxyyz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 151); 

                auto tg_xxxxxyyz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 152); 

                auto tg_xxxxxyyz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 153); 

                auto tg_xxxxxyyz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 154); 

                auto tg_xxxxxyyz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 155); 

                auto tg_xxxxxyyz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 156); 

                auto tg_xxxxxyyz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 157); 

                auto tg_xxxxxyyz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 158); 

                auto tg_xxxxxyyz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 159); 

                auto tg_xxxxxyyz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 160); 

                auto tg_xxxxxyyz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 161); 

                auto tg_xxxxxyyz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 162); 

                auto tg_xxxxxyyz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 163); 

                auto tg_xxxxxyyz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 164); 

                auto tg_xxxxxyyz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 165); 

                auto tg_xxxxxyyz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 166); 

                auto tg_xxxxxyyz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 167); 

                auto tg_xxxxxyzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 168); 

                auto tg_xxxxxyzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 169); 

                auto tg_xxxxxyzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 170); 

                auto tg_xxxxxyzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 171); 

                auto tg_xxxxxyzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 172); 

                auto tg_xxxxxyzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 173); 

                auto tg_xxxxxyzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 174); 

                auto tg_xxxxxyzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 175); 

                auto tg_xxxxxyzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 176); 

                auto tg_xxxxxyzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 177); 

                auto tg_xxxxxyzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 178); 

                auto tg_xxxxxyzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 179); 

                auto tg_xxxxxyzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 180); 

                auto tg_xxxxxyzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 181); 

                auto tg_xxxxxyzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 182); 

                auto tg_xxxxxyzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 183); 

                auto tg_xxxxxyzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 184); 

                auto tg_xxxxxyzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 185); 

                auto tg_xxxxxyzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 186); 

                auto tg_xxxxxyzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 187); 

                auto tg_xxxxxyzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 188); 

                auto tg_xxxxxzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 189); 

                // Batch of Integrals (95,190)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxxyz_xyyyz_0, tg_xxxxxxyz_xyyzz_0, \
                                         tg_xxxxxxyz_xyzzz_0, tg_xxxxxxyz_xzzzz_0, tg_xxxxxxyz_yyyyy_0, tg_xxxxxxyz_yyyyz_0, \
                                         tg_xxxxxxyz_yyyzz_0, tg_xxxxxxyz_yyzzz_0, tg_xxxxxxyz_yzzzz_0, tg_xxxxxxyz_zzzzz_0, \
                                         tg_xxxxxxzz_xxxxx_0, tg_xxxxxxzz_xxxxy_0, tg_xxxxxxzz_xxxxz_0, tg_xxxxxxzz_xxxyy_0, \
                                         tg_xxxxxxzz_xxxyz_0, tg_xxxxxxzz_xxxzz_0, tg_xxxxxxzz_xxyyy_0, tg_xxxxxxzz_xxyyz_0, \
                                         tg_xxxxxxzz_xxyzz_0, tg_xxxxxxzz_xxzzz_0, tg_xxxxxxzz_xyyyy_0, tg_xxxxxxzz_xyyyz_0, \
                                         tg_xxxxxxzz_xyyzz_0, tg_xxxxxxzz_xyzzz_0, tg_xxxxxxzz_xzzzz_0, tg_xxxxxxzz_yyyyy_0, \
                                         tg_xxxxxxzz_yyyyz_0, tg_xxxxxxzz_yyyzz_0, tg_xxxxxxzz_yyzzz_0, tg_xxxxxxzz_yzzzz_0, \
                                         tg_xxxxxxzz_zzzzz_0, tg_xxxxxyyy_xxxxx_0, tg_xxxxxyyy_xxxxy_0, tg_xxxxxyyy_xxxxz_0, \
                                         tg_xxxxxyyy_xxxyy_0, tg_xxxxxyyy_xxxyz_0, tg_xxxxxyyy_xxxzz_0, tg_xxxxxyyy_xxyyy_0, \
                                         tg_xxxxxyyy_xxyyz_0, tg_xxxxxyyy_xxyzz_0, tg_xxxxxyyy_xxzzz_0, tg_xxxxxyyy_xyyyy_0, \
                                         tg_xxxxxyyy_xyyyz_0, tg_xxxxxyyy_xyyzz_0, tg_xxxxxyyy_xyzzz_0, tg_xxxxxyyy_xzzzz_0, \
                                         tg_xxxxxyyy_yyyyy_0, tg_xxxxxyyy_yyyyz_0, tg_xxxxxyyy_yyyzz_0, tg_xxxxxyyy_yyzzz_0, \
                                         tg_xxxxxyyy_yzzzz_0, tg_xxxxxyyy_zzzzz_0, tg_xxxxxyyz_xxxxx_0, tg_xxxxxyyz_xxxxy_0, \
                                         tg_xxxxxyyz_xxxxz_0, tg_xxxxxyyz_xxxyy_0, tg_xxxxxyyz_xxxyz_0, tg_xxxxxyyz_xxxzz_0, \
                                         tg_xxxxxyyz_xxyyy_0, tg_xxxxxyyz_xxyyz_0, tg_xxxxxyyz_xxyzz_0, tg_xxxxxyyz_xxzzz_0, \
                                         tg_xxxxxyyz_xyyyy_0, tg_xxxxxyyz_xyyyz_0, tg_xxxxxyyz_xyyzz_0, tg_xxxxxyyz_xyzzz_0, \
                                         tg_xxxxxyyz_xzzzz_0, tg_xxxxxyyz_yyyyy_0, tg_xxxxxyyz_yyyyz_0, tg_xxxxxyyz_yyyzz_0, \
                                         tg_xxxxxyyz_yyzzz_0, tg_xxxxxyyz_yzzzz_0, tg_xxxxxyyz_zzzzz_0, tg_xxxxxyz_xyyyz_0, \
                                         tg_xxxxxyz_xyyyz_1, tg_xxxxxyz_xyyzz_0, tg_xxxxxyz_xyyzz_1, tg_xxxxxyz_xyzzz_0, \
                                         tg_xxxxxyz_xyzzz_1, tg_xxxxxyz_xzzzz_0, tg_xxxxxyz_xzzzz_1, tg_xxxxxyz_yyyyy_0, \
                                         tg_xxxxxyz_yyyyy_1, tg_xxxxxyz_yyyyz_0, tg_xxxxxyz_yyyyz_1, tg_xxxxxyz_yyyz_1, \
                                         tg_xxxxxyz_yyyzz_0, tg_xxxxxyz_yyyzz_1, tg_xxxxxyz_yyzz_1, tg_xxxxxyz_yyzzz_0, \
                                         tg_xxxxxyz_yyzzz_1, tg_xxxxxyz_yzzz_1, tg_xxxxxyz_yzzzz_0, tg_xxxxxyz_yzzzz_1, \
                                         tg_xxxxxyz_zzzz_1, tg_xxxxxyz_zzzzz_0, tg_xxxxxyz_zzzzz_1, tg_xxxxxyzz_xxxxx_0, \
                                         tg_xxxxxyzz_xxxxy_0, tg_xxxxxyzz_xxxxz_0, tg_xxxxxyzz_xxxyy_0, tg_xxxxxyzz_xxxyz_0, \
                                         tg_xxxxxyzz_xxxzz_0, tg_xxxxxyzz_xxyyy_0, tg_xxxxxyzz_xxyyz_0, tg_xxxxxyzz_xxyzz_0, \
                                         tg_xxxxxyzz_xxzzz_0, tg_xxxxxyzz_xyyyy_0, tg_xxxxxyzz_xyyyz_0, tg_xxxxxyzz_xyyzz_0, \
                                         tg_xxxxxyzz_xyzzz_0, tg_xxxxxyzz_xzzzz_0, tg_xxxxxyzz_yyyyy_0, tg_xxxxxyzz_yyyyz_0, \
                                         tg_xxxxxyzz_yyyzz_0, tg_xxxxxyzz_yyzzz_0, tg_xxxxxyzz_yzzzz_0, tg_xxxxxyzz_zzzzz_0, \
                                         tg_xxxxxzz_xxxx_1, tg_xxxxxzz_xxxxx_0, tg_xxxxxzz_xxxxx_1, tg_xxxxxzz_xxxxy_0, \
                                         tg_xxxxxzz_xxxxy_1, tg_xxxxxzz_xxxxz_0, tg_xxxxxzz_xxxxz_1, tg_xxxxxzz_xxxy_1, \
                                         tg_xxxxxzz_xxxyy_0, tg_xxxxxzz_xxxyy_1, tg_xxxxxzz_xxxyz_0, tg_xxxxxzz_xxxyz_1, \
                                         tg_xxxxxzz_xxxz_1, tg_xxxxxzz_xxxzz_0, tg_xxxxxzz_xxxzz_1, tg_xxxxxzz_xxyy_1, \
                                         tg_xxxxxzz_xxyyy_0, tg_xxxxxzz_xxyyy_1, tg_xxxxxzz_xxyyz_0, tg_xxxxxzz_xxyyz_1, \
                                         tg_xxxxxzz_xxyz_1, tg_xxxxxzz_xxyzz_0, tg_xxxxxzz_xxyzz_1, tg_xxxxxzz_xxzz_1, \
                                         tg_xxxxxzz_xxzzz_0, tg_xxxxxzz_xxzzz_1, tg_xxxxxzz_xyyy_1, tg_xxxxxzz_xyyyy_0, \
                                         tg_xxxxxzz_xyyyy_1, tg_xxxxxzz_xyyyz_0, tg_xxxxxzz_xyyyz_1, tg_xxxxxzz_xyyz_1, \
                                         tg_xxxxxzz_xyyzz_0, tg_xxxxxzz_xyyzz_1, tg_xxxxxzz_xyzz_1, tg_xxxxxzz_xyzzz_0, \
                                         tg_xxxxxzz_xyzzz_1, tg_xxxxxzz_xzzz_1, tg_xxxxxzz_xzzzz_0, tg_xxxxxzz_xzzzz_1, \
                                         tg_xxxxxzz_yyyy_1, tg_xxxxxzz_yyyyy_0, tg_xxxxxzz_yyyyy_1, tg_xxxxxzz_yyyyz_0, \
                                         tg_xxxxxzz_yyyyz_1, tg_xxxxxzz_yyyz_1, tg_xxxxxzz_yyyzz_0, tg_xxxxxzz_yyyzz_1, \
                                         tg_xxxxxzz_yyzz_1, tg_xxxxxzz_yyzzz_0, tg_xxxxxzz_yyzzz_1, tg_xxxxxzz_yzzz_1, \
                                         tg_xxxxxzz_yzzzz_0, tg_xxxxxzz_yzzzz_1, tg_xxxxxzz_zzzz_1, tg_xxxxxzz_zzzzz_0, \
                                         tg_xxxxxzz_zzzzz_1, tg_xxxxxzzz_xxxxx_0, tg_xxxxyyy_xxxx_1, tg_xxxxyyy_xxxxx_0, \
                                         tg_xxxxyyy_xxxxx_1, tg_xxxxyyy_xxxxy_0, tg_xxxxyyy_xxxxy_1, tg_xxxxyyy_xxxxz_0, \
                                         tg_xxxxyyy_xxxxz_1, tg_xxxxyyy_xxxy_1, tg_xxxxyyy_xxxyy_0, tg_xxxxyyy_xxxyy_1, \
                                         tg_xxxxyyy_xxxyz_0, tg_xxxxyyy_xxxyz_1, tg_xxxxyyy_xxxz_1, tg_xxxxyyy_xxxzz_0, \
                                         tg_xxxxyyy_xxxzz_1, tg_xxxxyyy_xxyy_1, tg_xxxxyyy_xxyyy_0, tg_xxxxyyy_xxyyy_1, \
                                         tg_xxxxyyy_xxyyz_0, tg_xxxxyyy_xxyyz_1, tg_xxxxyyy_xxyz_1, tg_xxxxyyy_xxyzz_0, \
                                         tg_xxxxyyy_xxyzz_1, tg_xxxxyyy_xxzz_1, tg_xxxxyyy_xxzzz_0, tg_xxxxyyy_xxzzz_1, \
                                         tg_xxxxyyy_xyyy_1, tg_xxxxyyy_xyyyy_0, tg_xxxxyyy_xyyyy_1, tg_xxxxyyy_xyyyz_0, \
                                         tg_xxxxyyy_xyyyz_1, tg_xxxxyyy_xyyz_1, tg_xxxxyyy_xyyzz_0, tg_xxxxyyy_xyyzz_1, \
                                         tg_xxxxyyy_xyzz_1, tg_xxxxyyy_xyzzz_0, tg_xxxxyyy_xyzzz_1, tg_xxxxyyy_xzzz_1, \
                                         tg_xxxxyyy_xzzzz_0, tg_xxxxyyy_xzzzz_1, tg_xxxxyyy_yyyy_1, tg_xxxxyyy_yyyyy_0, \
                                         tg_xxxxyyy_yyyyy_1, tg_xxxxyyy_yyyyz_0, tg_xxxxyyy_yyyyz_1, tg_xxxxyyy_yyyz_1, \
                                         tg_xxxxyyy_yyyzz_0, tg_xxxxyyy_yyyzz_1, tg_xxxxyyy_yyzz_1, tg_xxxxyyy_yyzzz_0, \
                                         tg_xxxxyyy_yyzzz_1, tg_xxxxyyy_yzzz_1, tg_xxxxyyy_yzzzz_0, tg_xxxxyyy_yzzzz_1, \
                                         tg_xxxxyyy_zzzz_1, tg_xxxxyyy_zzzzz_0, tg_xxxxyyy_zzzzz_1, tg_xxxxyyz_xxxx_1, \
                                         tg_xxxxyyz_xxxxx_0, tg_xxxxyyz_xxxxx_1, tg_xxxxyyz_xxxxy_0, tg_xxxxyyz_xxxxy_1, \
                                         tg_xxxxyyz_xxxxz_0, tg_xxxxyyz_xxxxz_1, tg_xxxxyyz_xxxy_1, tg_xxxxyyz_xxxyy_0, \
                                         tg_xxxxyyz_xxxyy_1, tg_xxxxyyz_xxxyz_0, tg_xxxxyyz_xxxyz_1, tg_xxxxyyz_xxxz_1, \
                                         tg_xxxxyyz_xxxzz_0, tg_xxxxyyz_xxxzz_1, tg_xxxxyyz_xxyy_1, tg_xxxxyyz_xxyyy_0, \
                                         tg_xxxxyyz_xxyyy_1, tg_xxxxyyz_xxyyz_0, tg_xxxxyyz_xxyyz_1, tg_xxxxyyz_xxyz_1, \
                                         tg_xxxxyyz_xxyzz_0, tg_xxxxyyz_xxyzz_1, tg_xxxxyyz_xxzz_1, tg_xxxxyyz_xxzzz_0, \
                                         tg_xxxxyyz_xxzzz_1, tg_xxxxyyz_xyyy_1, tg_xxxxyyz_xyyyy_0, tg_xxxxyyz_xyyyy_1, \
                                         tg_xxxxyyz_xyyyz_0, tg_xxxxyyz_xyyyz_1, tg_xxxxyyz_xyyz_1, tg_xxxxyyz_xyyzz_0, \
                                         tg_xxxxyyz_xyyzz_1, tg_xxxxyyz_xyzz_1, tg_xxxxyyz_xyzzz_0, tg_xxxxyyz_xyzzz_1, \
                                         tg_xxxxyyz_xzzz_1, tg_xxxxyyz_xzzzz_0, tg_xxxxyyz_xzzzz_1, tg_xxxxyyz_yyyy_1, \
                                         tg_xxxxyyz_yyyyy_0, tg_xxxxyyz_yyyyy_1, tg_xxxxyyz_yyyyz_0, tg_xxxxyyz_yyyyz_1, \
                                         tg_xxxxyyz_yyyz_1, tg_xxxxyyz_yyyzz_0, tg_xxxxyyz_yyyzz_1, tg_xxxxyyz_yyzz_1, \
                                         tg_xxxxyyz_yyzzz_0, tg_xxxxyyz_yyzzz_1, tg_xxxxyyz_yzzz_1, tg_xxxxyyz_yzzzz_0, \
                                         tg_xxxxyyz_yzzzz_1, tg_xxxxyyz_zzzz_1, tg_xxxxyyz_zzzzz_0, tg_xxxxyyz_zzzzz_1, \
                                         tg_xxxxyz_xyyyz_0, tg_xxxxyz_xyyyz_1, tg_xxxxyz_xyyzz_0, tg_xxxxyz_xyyzz_1, \
                                         tg_xxxxyz_xyzzz_0, tg_xxxxyz_xyzzz_1, tg_xxxxyz_xzzzz_0, tg_xxxxyz_xzzzz_1, \
                                         tg_xxxxyz_yyyyy_0, tg_xxxxyz_yyyyy_1, tg_xxxxyz_yyyyz_0, tg_xxxxyz_yyyyz_1, \
                                         tg_xxxxyz_yyyzz_0, tg_xxxxyz_yyyzz_1, tg_xxxxyz_yyzzz_0, tg_xxxxyz_yyzzz_1, \
                                         tg_xxxxyz_yzzzz_0, tg_xxxxyz_yzzzz_1, tg_xxxxyz_zzzzz_0, tg_xxxxyz_zzzzz_1, \
                                         tg_xxxxyzz_xxxx_1, tg_xxxxyzz_xxxxx_0, tg_xxxxyzz_xxxxx_1, tg_xxxxyzz_xxxxy_0, \
                                         tg_xxxxyzz_xxxxy_1, tg_xxxxyzz_xxxxz_0, tg_xxxxyzz_xxxxz_1, tg_xxxxyzz_xxxy_1, \
                                         tg_xxxxyzz_xxxyy_0, tg_xxxxyzz_xxxyy_1, tg_xxxxyzz_xxxyz_0, tg_xxxxyzz_xxxyz_1, \
                                         tg_xxxxyzz_xxxz_1, tg_xxxxyzz_xxxzz_0, tg_xxxxyzz_xxxzz_1, tg_xxxxyzz_xxyy_1, \
                                         tg_xxxxyzz_xxyyy_0, tg_xxxxyzz_xxyyy_1, tg_xxxxyzz_xxyyz_0, tg_xxxxyzz_xxyyz_1, \
                                         tg_xxxxyzz_xxyz_1, tg_xxxxyzz_xxyzz_0, tg_xxxxyzz_xxyzz_1, tg_xxxxyzz_xxzz_1, \
                                         tg_xxxxyzz_xxzzz_0, tg_xxxxyzz_xxzzz_1, tg_xxxxyzz_xyyy_1, tg_xxxxyzz_xyyyy_0, \
                                         tg_xxxxyzz_xyyyy_1, tg_xxxxyzz_xyyyz_0, tg_xxxxyzz_xyyyz_1, tg_xxxxyzz_xyyz_1, \
                                         tg_xxxxyzz_xyyzz_0, tg_xxxxyzz_xyyzz_1, tg_xxxxyzz_xyzz_1, tg_xxxxyzz_xyzzz_0, \
                                         tg_xxxxyzz_xyzzz_1, tg_xxxxyzz_xzzz_1, tg_xxxxyzz_xzzzz_0, tg_xxxxyzz_xzzzz_1, \
                                         tg_xxxxyzz_yyyy_1, tg_xxxxyzz_yyyyy_0, tg_xxxxyzz_yyyyy_1, tg_xxxxyzz_yyyyz_0, \
                                         tg_xxxxyzz_yyyyz_1, tg_xxxxyzz_yyyz_1, tg_xxxxyzz_yyyzz_0, tg_xxxxyzz_yyyzz_1, \
                                         tg_xxxxyzz_yyzz_1, tg_xxxxyzz_yyzzz_0, tg_xxxxyzz_yyzzz_1, tg_xxxxyzz_yzzz_1, \
                                         tg_xxxxyzz_yzzzz_0, tg_xxxxyzz_yzzzz_1, tg_xxxxyzz_zzzz_1, tg_xxxxyzz_zzzzz_0, \
                                         tg_xxxxyzz_zzzzz_1, tg_xxxxzz_xxxxx_0, tg_xxxxzz_xxxxx_1, tg_xxxxzz_xxxxy_0, \
                                         tg_xxxxzz_xxxxy_1, tg_xxxxzz_xxxxz_0, tg_xxxxzz_xxxxz_1, tg_xxxxzz_xxxyy_0, \
                                         tg_xxxxzz_xxxyy_1, tg_xxxxzz_xxxyz_0, tg_xxxxzz_xxxyz_1, tg_xxxxzz_xxxzz_0, \
                                         tg_xxxxzz_xxxzz_1, tg_xxxxzz_xxyyy_0, tg_xxxxzz_xxyyy_1, tg_xxxxzz_xxyyz_0, \
                                         tg_xxxxzz_xxyyz_1, tg_xxxxzz_xxyzz_0, tg_xxxxzz_xxyzz_1, tg_xxxxzz_xxzzz_0, \
                                         tg_xxxxzz_xxzzz_1, tg_xxxxzz_xyyyy_0, tg_xxxxzz_xyyyy_1, tg_xxxxzz_xyyyz_0, \
                                         tg_xxxxzz_xyyyz_1, tg_xxxxzz_xyyzz_0, tg_xxxxzz_xyyzz_1, tg_xxxxzz_xyzzz_0, \
                                         tg_xxxxzz_xyzzz_1, tg_xxxxzz_xzzzz_0, tg_xxxxzz_xzzzz_1, tg_xxxxzz_yyyyy_0, \
                                         tg_xxxxzz_yyyyy_1, tg_xxxxzz_yyyyz_0, tg_xxxxzz_yyyyz_1, tg_xxxxzz_yyyzz_0, \
                                         tg_xxxxzz_yyyzz_1, tg_xxxxzz_yyzzz_0, tg_xxxxzz_yyzzz_1, tg_xxxxzz_yzzzz_0, \
                                         tg_xxxxzz_yzzzz_1, tg_xxxxzz_zzzzz_0, tg_xxxxzz_zzzzz_1, tg_xxxxzzz_xxxx_1, \
                                         tg_xxxxzzz_xxxxx_0, tg_xxxxzzz_xxxxx_1, tg_xxxyyy_xxxxx_0, tg_xxxyyy_xxxxx_1, \
                                         tg_xxxyyy_xxxxy_0, tg_xxxyyy_xxxxy_1, tg_xxxyyy_xxxxz_0, tg_xxxyyy_xxxxz_1, \
                                         tg_xxxyyy_xxxyy_0, tg_xxxyyy_xxxyy_1, tg_xxxyyy_xxxyz_0, tg_xxxyyy_xxxyz_1, \
                                         tg_xxxyyy_xxxzz_0, tg_xxxyyy_xxxzz_1, tg_xxxyyy_xxyyy_0, tg_xxxyyy_xxyyy_1, \
                                         tg_xxxyyy_xxyyz_0, tg_xxxyyy_xxyyz_1, tg_xxxyyy_xxyzz_0, tg_xxxyyy_xxyzz_1, \
                                         tg_xxxyyy_xxzzz_0, tg_xxxyyy_xxzzz_1, tg_xxxyyy_xyyyy_0, tg_xxxyyy_xyyyy_1, \
                                         tg_xxxyyy_xyyyz_0, tg_xxxyyy_xyyyz_1, tg_xxxyyy_xyyzz_0, tg_xxxyyy_xyyzz_1, \
                                         tg_xxxyyy_xyzzz_0, tg_xxxyyy_xyzzz_1, tg_xxxyyy_xzzzz_0, tg_xxxyyy_xzzzz_1, \
                                         tg_xxxyyy_yyyyy_0, tg_xxxyyy_yyyyy_1, tg_xxxyyy_yyyyz_0, tg_xxxyyy_yyyyz_1, \
                                         tg_xxxyyy_yyyzz_0, tg_xxxyyy_yyyzz_1, tg_xxxyyy_yyzzz_0, tg_xxxyyy_yyzzz_1, \
                                         tg_xxxyyy_yzzzz_0, tg_xxxyyy_yzzzz_1, tg_xxxyyy_zzzzz_0, tg_xxxyyy_zzzzz_1, \
                                         tg_xxxyyz_xxxxx_0, tg_xxxyyz_xxxxx_1, tg_xxxyyz_xxxxy_0, tg_xxxyyz_xxxxy_1, \
                                         tg_xxxyyz_xxxxz_0, tg_xxxyyz_xxxxz_1, tg_xxxyyz_xxxyy_0, tg_xxxyyz_xxxyy_1, \
                                         tg_xxxyyz_xxxyz_0, tg_xxxyyz_xxxyz_1, tg_xxxyyz_xxxzz_0, tg_xxxyyz_xxxzz_1, \
                                         tg_xxxyyz_xxyyy_0, tg_xxxyyz_xxyyy_1, tg_xxxyyz_xxyyz_0, tg_xxxyyz_xxyyz_1, \
                                         tg_xxxyyz_xxyzz_0, tg_xxxyyz_xxyzz_1, tg_xxxyyz_xxzzz_0, tg_xxxyyz_xxzzz_1, \
                                         tg_xxxyyz_xyyyy_0, tg_xxxyyz_xyyyy_1, tg_xxxyyz_xyyyz_0, tg_xxxyyz_xyyyz_1, \
                                         tg_xxxyyz_xyyzz_0, tg_xxxyyz_xyyzz_1, tg_xxxyyz_xyzzz_0, tg_xxxyyz_xyzzz_1, \
                                         tg_xxxyyz_xzzzz_0, tg_xxxyyz_xzzzz_1, tg_xxxyyz_yyyyy_0, tg_xxxyyz_yyyyy_1, \
                                         tg_xxxyyz_yyyyz_0, tg_xxxyyz_yyyyz_1, tg_xxxyyz_yyyzz_0, tg_xxxyyz_yyyzz_1, \
                                         tg_xxxyyz_yyzzz_0, tg_xxxyyz_yyzzz_1, tg_xxxyyz_yzzzz_0, tg_xxxyyz_yzzzz_1, \
                                         tg_xxxyyz_zzzzz_0, tg_xxxyyz_zzzzz_1, tg_xxxyzz_xxxxx_0, tg_xxxyzz_xxxxx_1, \
                                         tg_xxxyzz_xxxxy_0, tg_xxxyzz_xxxxy_1, tg_xxxyzz_xxxxz_0, tg_xxxyzz_xxxxz_1, \
                                         tg_xxxyzz_xxxyy_0, tg_xxxyzz_xxxyy_1, tg_xxxyzz_xxxyz_0, tg_xxxyzz_xxxyz_1, \
                                         tg_xxxyzz_xxxzz_0, tg_xxxyzz_xxxzz_1, tg_xxxyzz_xxyyy_0, tg_xxxyzz_xxyyy_1, \
                                         tg_xxxyzz_xxyyz_0, tg_xxxyzz_xxyyz_1, tg_xxxyzz_xxyzz_0, tg_xxxyzz_xxyzz_1, \
                                         tg_xxxyzz_xxzzz_0, tg_xxxyzz_xxzzz_1, tg_xxxyzz_xyyyy_0, tg_xxxyzz_xyyyy_1, \
                                         tg_xxxyzz_xyyyz_0, tg_xxxyzz_xyyyz_1, tg_xxxyzz_xyyzz_0, tg_xxxyzz_xyyzz_1, \
                                         tg_xxxyzz_xyzzz_0, tg_xxxyzz_xyzzz_1, tg_xxxyzz_xzzzz_0, tg_xxxyzz_xzzzz_1, \
                                         tg_xxxyzz_yyyyy_0, tg_xxxyzz_yyyyy_1, tg_xxxyzz_yyyyz_0, tg_xxxyzz_yyyyz_1, \
                                         tg_xxxyzz_yyyzz_0, tg_xxxyzz_yyyzz_1, tg_xxxyzz_yyzzz_0, tg_xxxyzz_yyzzz_1, \
                                         tg_xxxyzz_yzzzz_0, tg_xxxyzz_yzzzz_1, tg_xxxyzz_zzzzz_0, tg_xxxyzz_zzzzz_1, \
                                         tg_xxxzzz_xxxxx_0, tg_xxxzzz_xxxxx_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxxyz_xyyyz_0[j] = pb_x * tg_xxxxxyz_xyyyz_0[j] + fr * tg_xxxxxyz_xyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xyyyz_0[j] - tg_xxxxyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_yyyz_1[j];

                    tg_xxxxxxyz_xyyzz_0[j] = pb_x * tg_xxxxxyz_xyyzz_0[j] + fr * tg_xxxxxyz_xyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xyyzz_0[j] - tg_xxxxyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_yyzz_1[j];

                    tg_xxxxxxyz_xyzzz_0[j] = pb_x * tg_xxxxxyz_xyzzz_0[j] + fr * tg_xxxxxyz_xyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xyzzz_0[j] - tg_xxxxyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_yzzz_1[j];

                    tg_xxxxxxyz_xzzzz_0[j] = pb_x * tg_xxxxxyz_xzzzz_0[j] + fr * tg_xxxxxyz_xzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xzzzz_0[j] - tg_xxxxyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_zzzz_1[j];

                    tg_xxxxxxyz_yyyyy_0[j] = pb_x * tg_xxxxxyz_yyyyy_0[j] + fr * tg_xxxxxyz_yyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yyyyy_0[j] - tg_xxxxyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxyz_yyyyz_0[j] = pb_x * tg_xxxxxyz_yyyyz_0[j] + fr * tg_xxxxxyz_yyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yyyyz_0[j] - tg_xxxxyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxyz_yyyzz_0[j] = pb_x * tg_xxxxxyz_yyyzz_0[j] + fr * tg_xxxxxyz_yyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yyyzz_0[j] - tg_xxxxyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxyz_yyzzz_0[j] = pb_x * tg_xxxxxyz_yyzzz_0[j] + fr * tg_xxxxxyz_yyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yyzzz_0[j] - tg_xxxxyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxyz_yzzzz_0[j] = pb_x * tg_xxxxxyz_yzzzz_0[j] + fr * tg_xxxxxyz_yzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yzzzz_0[j] - tg_xxxxyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxyz_zzzzz_0[j] = pb_x * tg_xxxxxyz_zzzzz_0[j] + fr * tg_xxxxxyz_zzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_zzzzz_0[j] - tg_xxxxyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxxzz_xxxxx_0[j] = pb_x * tg_xxxxxzz_xxxxx_0[j] + fr * tg_xxxxxzz_xxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxxxx_0[j] - tg_xxxxzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxzz_xxxx_1[j];

                    tg_xxxxxxzz_xxxxy_0[j] = pb_x * tg_xxxxxzz_xxxxy_0[j] + fr * tg_xxxxxzz_xxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxxxy_0[j] - tg_xxxxzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxzz_xxxy_1[j];

                    tg_xxxxxxzz_xxxxz_0[j] = pb_x * tg_xxxxxzz_xxxxz_0[j] + fr * tg_xxxxxzz_xxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxxxz_0[j] - tg_xxxxzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxzz_xxxz_1[j];

                    tg_xxxxxxzz_xxxyy_0[j] = pb_x * tg_xxxxxzz_xxxyy_0[j] + fr * tg_xxxxxzz_xxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxxyy_0[j] - tg_xxxxzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxzz_xxyy_1[j];

                    tg_xxxxxxzz_xxxyz_0[j] = pb_x * tg_xxxxxzz_xxxyz_0[j] + fr * tg_xxxxxzz_xxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxxyz_0[j] - tg_xxxxzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxzz_xxyz_1[j];

                    tg_xxxxxxzz_xxxzz_0[j] = pb_x * tg_xxxxxzz_xxxzz_0[j] + fr * tg_xxxxxzz_xxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxxzz_0[j] - tg_xxxxzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxzz_xxzz_1[j];

                    tg_xxxxxxzz_xxyyy_0[j] = pb_x * tg_xxxxxzz_xxyyy_0[j] + fr * tg_xxxxxzz_xxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxyyy_0[j] - tg_xxxxzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxzz_xyyy_1[j];

                    tg_xxxxxxzz_xxyyz_0[j] = pb_x * tg_xxxxxzz_xxyyz_0[j] + fr * tg_xxxxxzz_xxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxyyz_0[j] - tg_xxxxzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxzz_xyyz_1[j];

                    tg_xxxxxxzz_xxyzz_0[j] = pb_x * tg_xxxxxzz_xxyzz_0[j] + fr * tg_xxxxxzz_xxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxyzz_0[j] - tg_xxxxzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxzz_xyzz_1[j];

                    tg_xxxxxxzz_xxzzz_0[j] = pb_x * tg_xxxxxzz_xxzzz_0[j] + fr * tg_xxxxxzz_xxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxzzz_0[j] - tg_xxxxzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxzz_xzzz_1[j];

                    tg_xxxxxxzz_xyyyy_0[j] = pb_x * tg_xxxxxzz_xyyyy_0[j] + fr * tg_xxxxxzz_xyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xyyyy_0[j] - tg_xxxxzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_yyyy_1[j];

                    tg_xxxxxxzz_xyyyz_0[j] = pb_x * tg_xxxxxzz_xyyyz_0[j] + fr * tg_xxxxxzz_xyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xyyyz_0[j] - tg_xxxxzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_yyyz_1[j];

                    tg_xxxxxxzz_xyyzz_0[j] = pb_x * tg_xxxxxzz_xyyzz_0[j] + fr * tg_xxxxxzz_xyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xyyzz_0[j] - tg_xxxxzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_yyzz_1[j];

                    tg_xxxxxxzz_xyzzz_0[j] = pb_x * tg_xxxxxzz_xyzzz_0[j] + fr * tg_xxxxxzz_xyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xyzzz_0[j] - tg_xxxxzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_yzzz_1[j];

                    tg_xxxxxxzz_xzzzz_0[j] = pb_x * tg_xxxxxzz_xzzzz_0[j] + fr * tg_xxxxxzz_xzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xzzzz_0[j] - tg_xxxxzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_zzzz_1[j];

                    tg_xxxxxxzz_yyyyy_0[j] = pb_x * tg_xxxxxzz_yyyyy_0[j] + fr * tg_xxxxxzz_yyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yyyyy_0[j] - tg_xxxxzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxxzz_yyyyz_0[j] = pb_x * tg_xxxxxzz_yyyyz_0[j] + fr * tg_xxxxxzz_yyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yyyyz_0[j] - tg_xxxxzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxxzz_yyyzz_0[j] = pb_x * tg_xxxxxzz_yyyzz_0[j] + fr * tg_xxxxxzz_yyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yyyzz_0[j] - tg_xxxxzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxxzz_yyzzz_0[j] = pb_x * tg_xxxxxzz_yyzzz_0[j] + fr * tg_xxxxxzz_yyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yyzzz_0[j] - tg_xxxxzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxxzz_yzzzz_0[j] = pb_x * tg_xxxxxzz_yzzzz_0[j] + fr * tg_xxxxxzz_yzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yzzzz_0[j] - tg_xxxxzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxxzz_zzzzz_0[j] = pb_x * tg_xxxxxzz_zzzzz_0[j] + fr * tg_xxxxxzz_zzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_zzzzz_0[j] - tg_xxxxzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxyyy_xxxxx_0[j] = pb_x * tg_xxxxyyy_xxxxx_0[j] + fr * tg_xxxxyyy_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxxxx_0[j] - tg_xxxyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyyy_xxxx_1[j];

                    tg_xxxxxyyy_xxxxy_0[j] = pb_x * tg_xxxxyyy_xxxxy_0[j] + fr * tg_xxxxyyy_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxxxy_0[j] - tg_xxxyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyyy_xxxy_1[j];

                    tg_xxxxxyyy_xxxxz_0[j] = pb_x * tg_xxxxyyy_xxxxz_0[j] + fr * tg_xxxxyyy_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxxxz_0[j] - tg_xxxyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyyy_xxxz_1[j];

                    tg_xxxxxyyy_xxxyy_0[j] = pb_x * tg_xxxxyyy_xxxyy_0[j] + fr * tg_xxxxyyy_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxxyy_0[j] - tg_xxxyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyy_xxyy_1[j];

                    tg_xxxxxyyy_xxxyz_0[j] = pb_x * tg_xxxxyyy_xxxyz_0[j] + fr * tg_xxxxyyy_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxxyz_0[j] - tg_xxxyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyy_xxyz_1[j];

                    tg_xxxxxyyy_xxxzz_0[j] = pb_x * tg_xxxxyyy_xxxzz_0[j] + fr * tg_xxxxyyy_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxxzz_0[j] - tg_xxxyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyy_xxzz_1[j];

                    tg_xxxxxyyy_xxyyy_0[j] = pb_x * tg_xxxxyyy_xxyyy_0[j] + fr * tg_xxxxyyy_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxyyy_0[j] - tg_xxxyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyy_xyyy_1[j];

                    tg_xxxxxyyy_xxyyz_0[j] = pb_x * tg_xxxxyyy_xxyyz_0[j] + fr * tg_xxxxyyy_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxyyz_0[j] - tg_xxxyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyy_xyyz_1[j];

                    tg_xxxxxyyy_xxyzz_0[j] = pb_x * tg_xxxxyyy_xxyzz_0[j] + fr * tg_xxxxyyy_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxyzz_0[j] - tg_xxxyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyy_xyzz_1[j];

                    tg_xxxxxyyy_xxzzz_0[j] = pb_x * tg_xxxxyyy_xxzzz_0[j] + fr * tg_xxxxyyy_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxzzz_0[j] - tg_xxxyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyy_xzzz_1[j];

                    tg_xxxxxyyy_xyyyy_0[j] = pb_x * tg_xxxxyyy_xyyyy_0[j] + fr * tg_xxxxyyy_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xyyyy_0[j] - tg_xxxyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_yyyy_1[j];

                    tg_xxxxxyyy_xyyyz_0[j] = pb_x * tg_xxxxyyy_xyyyz_0[j] + fr * tg_xxxxyyy_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xyyyz_0[j] - tg_xxxyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_yyyz_1[j];

                    tg_xxxxxyyy_xyyzz_0[j] = pb_x * tg_xxxxyyy_xyyzz_0[j] + fr * tg_xxxxyyy_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xyyzz_0[j] - tg_xxxyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_yyzz_1[j];

                    tg_xxxxxyyy_xyzzz_0[j] = pb_x * tg_xxxxyyy_xyzzz_0[j] + fr * tg_xxxxyyy_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xyzzz_0[j] - tg_xxxyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_yzzz_1[j];

                    tg_xxxxxyyy_xzzzz_0[j] = pb_x * tg_xxxxyyy_xzzzz_0[j] + fr * tg_xxxxyyy_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xzzzz_0[j] - tg_xxxyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_zzzz_1[j];

                    tg_xxxxxyyy_yyyyy_0[j] = pb_x * tg_xxxxyyy_yyyyy_0[j] + fr * tg_xxxxyyy_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yyyyy_0[j] - tg_xxxyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxyyy_yyyyz_0[j] = pb_x * tg_xxxxyyy_yyyyz_0[j] + fr * tg_xxxxyyy_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yyyyz_0[j] - tg_xxxyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxyyy_yyyzz_0[j] = pb_x * tg_xxxxyyy_yyyzz_0[j] + fr * tg_xxxxyyy_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yyyzz_0[j] - tg_xxxyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxyyy_yyzzz_0[j] = pb_x * tg_xxxxyyy_yyzzz_0[j] + fr * tg_xxxxyyy_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yyzzz_0[j] - tg_xxxyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxyyy_yzzzz_0[j] = pb_x * tg_xxxxyyy_yzzzz_0[j] + fr * tg_xxxxyyy_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yzzzz_0[j] - tg_xxxyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxyyy_zzzzz_0[j] = pb_x * tg_xxxxyyy_zzzzz_0[j] + fr * tg_xxxxyyy_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_zzzzz_0[j] - tg_xxxyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxyyz_xxxxx_0[j] = pb_x * tg_xxxxyyz_xxxxx_0[j] + fr * tg_xxxxyyz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxxxx_0[j] - tg_xxxyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyyz_xxxx_1[j];

                    tg_xxxxxyyz_xxxxy_0[j] = pb_x * tg_xxxxyyz_xxxxy_0[j] + fr * tg_xxxxyyz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxxxy_0[j] - tg_xxxyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyyz_xxxy_1[j];

                    tg_xxxxxyyz_xxxxz_0[j] = pb_x * tg_xxxxyyz_xxxxz_0[j] + fr * tg_xxxxyyz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxxxz_0[j] - tg_xxxyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyyz_xxxz_1[j];

                    tg_xxxxxyyz_xxxyy_0[j] = pb_x * tg_xxxxyyz_xxxyy_0[j] + fr * tg_xxxxyyz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxxyy_0[j] - tg_xxxyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyz_xxyy_1[j];

                    tg_xxxxxyyz_xxxyz_0[j] = pb_x * tg_xxxxyyz_xxxyz_0[j] + fr * tg_xxxxyyz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxxyz_0[j] - tg_xxxyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyz_xxyz_1[j];

                    tg_xxxxxyyz_xxxzz_0[j] = pb_x * tg_xxxxyyz_xxxzz_0[j] + fr * tg_xxxxyyz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxxzz_0[j] - tg_xxxyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyz_xxzz_1[j];

                    tg_xxxxxyyz_xxyyy_0[j] = pb_x * tg_xxxxyyz_xxyyy_0[j] + fr * tg_xxxxyyz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxyyy_0[j] - tg_xxxyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyz_xyyy_1[j];

                    tg_xxxxxyyz_xxyyz_0[j] = pb_x * tg_xxxxyyz_xxyyz_0[j] + fr * tg_xxxxyyz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxyyz_0[j] - tg_xxxyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyz_xyyz_1[j];

                    tg_xxxxxyyz_xxyzz_0[j] = pb_x * tg_xxxxyyz_xxyzz_0[j] + fr * tg_xxxxyyz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxyzz_0[j] - tg_xxxyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyz_xyzz_1[j];

                    tg_xxxxxyyz_xxzzz_0[j] = pb_x * tg_xxxxyyz_xxzzz_0[j] + fr * tg_xxxxyyz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxzzz_0[j] - tg_xxxyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyz_xzzz_1[j];

                    tg_xxxxxyyz_xyyyy_0[j] = pb_x * tg_xxxxyyz_xyyyy_0[j] + fr * tg_xxxxyyz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xyyyy_0[j] - tg_xxxyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_yyyy_1[j];

                    tg_xxxxxyyz_xyyyz_0[j] = pb_x * tg_xxxxyyz_xyyyz_0[j] + fr * tg_xxxxyyz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xyyyz_0[j] - tg_xxxyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_yyyz_1[j];

                    tg_xxxxxyyz_xyyzz_0[j] = pb_x * tg_xxxxyyz_xyyzz_0[j] + fr * tg_xxxxyyz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xyyzz_0[j] - tg_xxxyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_yyzz_1[j];

                    tg_xxxxxyyz_xyzzz_0[j] = pb_x * tg_xxxxyyz_xyzzz_0[j] + fr * tg_xxxxyyz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xyzzz_0[j] - tg_xxxyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_yzzz_1[j];

                    tg_xxxxxyyz_xzzzz_0[j] = pb_x * tg_xxxxyyz_xzzzz_0[j] + fr * tg_xxxxyyz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xzzzz_0[j] - tg_xxxyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_zzzz_1[j];

                    tg_xxxxxyyz_yyyyy_0[j] = pb_x * tg_xxxxyyz_yyyyy_0[j] + fr * tg_xxxxyyz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yyyyy_0[j] - tg_xxxyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxyyz_yyyyz_0[j] = pb_x * tg_xxxxyyz_yyyyz_0[j] + fr * tg_xxxxyyz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yyyyz_0[j] - tg_xxxyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxyyz_yyyzz_0[j] = pb_x * tg_xxxxyyz_yyyzz_0[j] + fr * tg_xxxxyyz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yyyzz_0[j] - tg_xxxyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxyyz_yyzzz_0[j] = pb_x * tg_xxxxyyz_yyzzz_0[j] + fr * tg_xxxxyyz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yyzzz_0[j] - tg_xxxyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxyyz_yzzzz_0[j] = pb_x * tg_xxxxyyz_yzzzz_0[j] + fr * tg_xxxxyyz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yzzzz_0[j] - tg_xxxyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxyyz_zzzzz_0[j] = pb_x * tg_xxxxyyz_zzzzz_0[j] + fr * tg_xxxxyyz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_zzzzz_0[j] - tg_xxxyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxyzz_xxxxx_0[j] = pb_x * tg_xxxxyzz_xxxxx_0[j] + fr * tg_xxxxyzz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxxxx_0[j] - tg_xxxyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyzz_xxxx_1[j];

                    tg_xxxxxyzz_xxxxy_0[j] = pb_x * tg_xxxxyzz_xxxxy_0[j] + fr * tg_xxxxyzz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxxxy_0[j] - tg_xxxyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyzz_xxxy_1[j];

                    tg_xxxxxyzz_xxxxz_0[j] = pb_x * tg_xxxxyzz_xxxxz_0[j] + fr * tg_xxxxyzz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxxxz_0[j] - tg_xxxyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyzz_xxxz_1[j];

                    tg_xxxxxyzz_xxxyy_0[j] = pb_x * tg_xxxxyzz_xxxyy_0[j] + fr * tg_xxxxyzz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxxyy_0[j] - tg_xxxyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyzz_xxyy_1[j];

                    tg_xxxxxyzz_xxxyz_0[j] = pb_x * tg_xxxxyzz_xxxyz_0[j] + fr * tg_xxxxyzz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxxyz_0[j] - tg_xxxyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyzz_xxyz_1[j];

                    tg_xxxxxyzz_xxxzz_0[j] = pb_x * tg_xxxxyzz_xxxzz_0[j] + fr * tg_xxxxyzz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxxzz_0[j] - tg_xxxyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyzz_xxzz_1[j];

                    tg_xxxxxyzz_xxyyy_0[j] = pb_x * tg_xxxxyzz_xxyyy_0[j] + fr * tg_xxxxyzz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxyyy_0[j] - tg_xxxyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyzz_xyyy_1[j];

                    tg_xxxxxyzz_xxyyz_0[j] = pb_x * tg_xxxxyzz_xxyyz_0[j] + fr * tg_xxxxyzz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxyyz_0[j] - tg_xxxyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyzz_xyyz_1[j];

                    tg_xxxxxyzz_xxyzz_0[j] = pb_x * tg_xxxxyzz_xxyzz_0[j] + fr * tg_xxxxyzz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxyzz_0[j] - tg_xxxyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyzz_xyzz_1[j];

                    tg_xxxxxyzz_xxzzz_0[j] = pb_x * tg_xxxxyzz_xxzzz_0[j] + fr * tg_xxxxyzz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxzzz_0[j] - tg_xxxyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyzz_xzzz_1[j];

                    tg_xxxxxyzz_xyyyy_0[j] = pb_x * tg_xxxxyzz_xyyyy_0[j] + fr * tg_xxxxyzz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xyyyy_0[j] - tg_xxxyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_yyyy_1[j];

                    tg_xxxxxyzz_xyyyz_0[j] = pb_x * tg_xxxxyzz_xyyyz_0[j] + fr * tg_xxxxyzz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xyyyz_0[j] - tg_xxxyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_yyyz_1[j];

                    tg_xxxxxyzz_xyyzz_0[j] = pb_x * tg_xxxxyzz_xyyzz_0[j] + fr * tg_xxxxyzz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xyyzz_0[j] - tg_xxxyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_yyzz_1[j];

                    tg_xxxxxyzz_xyzzz_0[j] = pb_x * tg_xxxxyzz_xyzzz_0[j] + fr * tg_xxxxyzz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xyzzz_0[j] - tg_xxxyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_yzzz_1[j];

                    tg_xxxxxyzz_xzzzz_0[j] = pb_x * tg_xxxxyzz_xzzzz_0[j] + fr * tg_xxxxyzz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xzzzz_0[j] - tg_xxxyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_zzzz_1[j];

                    tg_xxxxxyzz_yyyyy_0[j] = pb_x * tg_xxxxyzz_yyyyy_0[j] + fr * tg_xxxxyzz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yyyyy_0[j] - tg_xxxyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxyzz_yyyyz_0[j] = pb_x * tg_xxxxyzz_yyyyz_0[j] + fr * tg_xxxxyzz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yyyyz_0[j] - tg_xxxyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxyzz_yyyzz_0[j] = pb_x * tg_xxxxyzz_yyyzz_0[j] + fr * tg_xxxxyzz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yyyzz_0[j] - tg_xxxyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxyzz_yyzzz_0[j] = pb_x * tg_xxxxyzz_yyzzz_0[j] + fr * tg_xxxxyzz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yyzzz_0[j] - tg_xxxyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxyzz_yzzzz_0[j] = pb_x * tg_xxxxyzz_yzzzz_0[j] + fr * tg_xxxxyzz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yzzzz_0[j] - tg_xxxyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxyzz_zzzzz_0[j] = pb_x * tg_xxxxyzz_zzzzz_0[j] + fr * tg_xxxxyzz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_zzzzz_0[j] - tg_xxxyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxxzzz_xxxxx_0[j] = pb_x * tg_xxxxzzz_xxxxx_0[j] + fr * tg_xxxxzzz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxxxx_0[j] - tg_xxxzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxzzz_xxxx_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_190_285(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxxxzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 190); 

                auto tg_xxxxzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 191); 

                auto tg_xxxxzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 192); 

                auto tg_xxxxzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 193); 

                auto tg_xxxxzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 194); 

                auto tg_xxxxzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 195); 

                auto tg_xxxxzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 196); 

                auto tg_xxxxzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 197); 

                auto tg_xxxxzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 198); 

                auto tg_xxxxzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 199); 

                auto tg_xxxxzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 200); 

                auto tg_xxxxzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 201); 

                auto tg_xxxxzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 202); 

                auto tg_xxxxzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 203); 

                auto tg_xxxxzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 204); 

                auto tg_xxxxzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 205); 

                auto tg_xxxxzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 206); 

                auto tg_xxxxzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 207); 

                auto tg_xxxxzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 208); 

                auto tg_xxxxzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 209); 

                auto tg_xxxyyyy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 210); 

                auto tg_xxxyyyy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 211); 

                auto tg_xxxyyyy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 212); 

                auto tg_xxxyyyy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 213); 

                auto tg_xxxyyyy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 214); 

                auto tg_xxxyyyy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 215); 

                auto tg_xxxyyyy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 216); 

                auto tg_xxxyyyy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 217); 

                auto tg_xxxyyyy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 218); 

                auto tg_xxxyyyy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 219); 

                auto tg_xxxyyyy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 220); 

                auto tg_xxxyyyy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 221); 

                auto tg_xxxyyyy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 222); 

                auto tg_xxxyyyy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 223); 

                auto tg_xxxyyyy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 224); 

                auto tg_xxxyyyy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 225); 

                auto tg_xxxyyyy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 226); 

                auto tg_xxxyyyy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 227); 

                auto tg_xxxyyyy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 228); 

                auto tg_xxxyyyy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 229); 

                auto tg_xxxyyyy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 230); 

                auto tg_xxxyyyz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 231); 

                auto tg_xxxyyyz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 232); 

                auto tg_xxxyyyz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 233); 

                auto tg_xxxyyyz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 234); 

                auto tg_xxxyyyz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 235); 

                auto tg_xxxyyyz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 236); 

                auto tg_xxxyyyz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 237); 

                auto tg_xxxyyyz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 238); 

                auto tg_xxxyyyz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 239); 

                auto tg_xxxyyyz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 240); 

                auto tg_xxxyyyz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 241); 

                auto tg_xxxyyyz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 242); 

                auto tg_xxxyyyz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 243); 

                auto tg_xxxyyyz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 244); 

                auto tg_xxxyyyz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 245); 

                auto tg_xxxyyyz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 246); 

                auto tg_xxxyyyz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 247); 

                auto tg_xxxyyyz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 248); 

                auto tg_xxxyyyz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 249); 

                auto tg_xxxyyyz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 250); 

                auto tg_xxxyyyz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 251); 

                auto tg_xxxyyzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 252); 

                auto tg_xxxyyzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 253); 

                auto tg_xxxyyzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 254); 

                auto tg_xxxyyzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 255); 

                auto tg_xxxyyzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 256); 

                auto tg_xxxyyzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 257); 

                auto tg_xxxyyzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 258); 

                auto tg_xxxyyzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 259); 

                auto tg_xxxyyzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 260); 

                auto tg_xxxyyzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 261); 

                auto tg_xxxyyzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 262); 

                auto tg_xxxyyzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 263); 

                auto tg_xxxyyzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 264); 

                auto tg_xxxyyzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 265); 

                auto tg_xxxyyzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 266); 

                auto tg_xxxyyzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 267); 

                auto tg_xxxyyzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 268); 

                auto tg_xxxyyzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 269); 

                auto tg_xxxyyzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 270); 

                auto tg_xxxyyzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 271); 

                auto tg_xxxyyzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 272); 

                auto tg_xxxyzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 273); 

                auto tg_xxxyzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 274); 

                auto tg_xxxyzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 275); 

                auto tg_xxxyzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 276); 

                auto tg_xxxyzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 277); 

                auto tg_xxxyzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 278); 

                auto tg_xxxyzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 279); 

                auto tg_xxxyzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 280); 

                auto tg_xxxyzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 281); 

                auto tg_xxxyzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 282); 

                auto tg_xxxyzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 283); 

                auto tg_xxxyzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 284); 

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

                auto tg_xxxxzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 136); 

                auto tg_xxxxzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 137); 

                auto tg_xxxxzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 138); 

                auto tg_xxxxzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 139); 

                auto tg_xxxxzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 140); 

                auto tg_xxxxzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 141); 

                auto tg_xxxxzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 142); 

                auto tg_xxxxzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 143); 

                auto tg_xxxxzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 144); 

                auto tg_xxxxzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 145); 

                auto tg_xxxxzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 146); 

                auto tg_xxxxzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 147); 

                auto tg_xxxxzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 148); 

                auto tg_xxxxzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 149); 

                auto tg_xxxyyyy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 150); 

                auto tg_xxxyyyy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 151); 

                auto tg_xxxyyyy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 152); 

                auto tg_xxxyyyy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 153); 

                auto tg_xxxyyyy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 154); 

                auto tg_xxxyyyy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 155); 

                auto tg_xxxyyyy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 156); 

                auto tg_xxxyyyy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 157); 

                auto tg_xxxyyyy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 158); 

                auto tg_xxxyyyy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 159); 

                auto tg_xxxyyyy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 160); 

                auto tg_xxxyyyy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 161); 

                auto tg_xxxyyyy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 162); 

                auto tg_xxxyyyy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 163); 

                auto tg_xxxyyyy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 164); 

                auto tg_xxxyyyz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 165); 

                auto tg_xxxyyyz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 166); 

                auto tg_xxxyyyz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 167); 

                auto tg_xxxyyyz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 168); 

                auto tg_xxxyyyz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 169); 

                auto tg_xxxyyyz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 170); 

                auto tg_xxxyyyz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 171); 

                auto tg_xxxyyyz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 172); 

                auto tg_xxxyyyz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 173); 

                auto tg_xxxyyyz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 174); 

                auto tg_xxxyyyz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 175); 

                auto tg_xxxyyyz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 176); 

                auto tg_xxxyyyz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 177); 

                auto tg_xxxyyyz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 178); 

                auto tg_xxxyyyz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 179); 

                auto tg_xxxyyzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 180); 

                auto tg_xxxyyzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 181); 

                auto tg_xxxyyzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 182); 

                auto tg_xxxyyzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 183); 

                auto tg_xxxyyzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 184); 

                auto tg_xxxyyzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 185); 

                auto tg_xxxyyzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 186); 

                auto tg_xxxyyzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 187); 

                auto tg_xxxyyzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 188); 

                auto tg_xxxyyzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 189); 

                auto tg_xxxyyzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 190); 

                auto tg_xxxyyzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 191); 

                auto tg_xxxyyzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 192); 

                auto tg_xxxyyzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 193); 

                auto tg_xxxyyzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 194); 

                auto tg_xxxyzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 195); 

                auto tg_xxxyzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 196); 

                auto tg_xxxyzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 197); 

                auto tg_xxxyzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 198); 

                auto tg_xxxyzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 199); 

                auto tg_xxxyzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 200); 

                auto tg_xxxyzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 201); 

                auto tg_xxxyzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 202); 

                auto tg_xxxyzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 203); 

                auto tg_xxxyzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 204); 

                auto tg_xxxyzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 205); 

                auto tg_xxxyzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 206); 

                // set up pointers to integrals

                auto tg_xxxxxzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 190); 

                auto tg_xxxxxzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 191); 

                auto tg_xxxxxzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 192); 

                auto tg_xxxxxzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 193); 

                auto tg_xxxxxzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 194); 

                auto tg_xxxxxzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 195); 

                auto tg_xxxxxzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 196); 

                auto tg_xxxxxzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 197); 

                auto tg_xxxxxzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 198); 

                auto tg_xxxxxzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 199); 

                auto tg_xxxxxzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 200); 

                auto tg_xxxxxzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 201); 

                auto tg_xxxxxzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 202); 

                auto tg_xxxxxzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 203); 

                auto tg_xxxxxzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 204); 

                auto tg_xxxxxzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 205); 

                auto tg_xxxxxzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 206); 

                auto tg_xxxxxzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 207); 

                auto tg_xxxxxzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 208); 

                auto tg_xxxxxzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 209); 

                auto tg_xxxxyyyy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 210); 

                auto tg_xxxxyyyy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 211); 

                auto tg_xxxxyyyy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 212); 

                auto tg_xxxxyyyy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 213); 

                auto tg_xxxxyyyy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 214); 

                auto tg_xxxxyyyy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 215); 

                auto tg_xxxxyyyy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 216); 

                auto tg_xxxxyyyy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 217); 

                auto tg_xxxxyyyy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 218); 

                auto tg_xxxxyyyy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 219); 

                auto tg_xxxxyyyy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 220); 

                auto tg_xxxxyyyy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 221); 

                auto tg_xxxxyyyy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 222); 

                auto tg_xxxxyyyy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 223); 

                auto tg_xxxxyyyy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 224); 

                auto tg_xxxxyyyy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 225); 

                auto tg_xxxxyyyy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 226); 

                auto tg_xxxxyyyy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 227); 

                auto tg_xxxxyyyy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 228); 

                auto tg_xxxxyyyy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 229); 

                auto tg_xxxxyyyy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 230); 

                auto tg_xxxxyyyz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 231); 

                auto tg_xxxxyyyz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 232); 

                auto tg_xxxxyyyz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 233); 

                auto tg_xxxxyyyz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 234); 

                auto tg_xxxxyyyz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 235); 

                auto tg_xxxxyyyz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 236); 

                auto tg_xxxxyyyz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 237); 

                auto tg_xxxxyyyz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 238); 

                auto tg_xxxxyyyz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 239); 

                auto tg_xxxxyyyz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 240); 

                auto tg_xxxxyyyz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 241); 

                auto tg_xxxxyyyz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 242); 

                auto tg_xxxxyyyz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 243); 

                auto tg_xxxxyyyz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 244); 

                auto tg_xxxxyyyz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 245); 

                auto tg_xxxxyyyz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 246); 

                auto tg_xxxxyyyz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 247); 

                auto tg_xxxxyyyz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 248); 

                auto tg_xxxxyyyz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 249); 

                auto tg_xxxxyyyz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 250); 

                auto tg_xxxxyyyz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 251); 

                auto tg_xxxxyyzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 252); 

                auto tg_xxxxyyzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 253); 

                auto tg_xxxxyyzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 254); 

                auto tg_xxxxyyzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 255); 

                auto tg_xxxxyyzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 256); 

                auto tg_xxxxyyzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 257); 

                auto tg_xxxxyyzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 258); 

                auto tg_xxxxyyzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 259); 

                auto tg_xxxxyyzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 260); 

                auto tg_xxxxyyzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 261); 

                auto tg_xxxxyyzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 262); 

                auto tg_xxxxyyzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 263); 

                auto tg_xxxxyyzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 264); 

                auto tg_xxxxyyzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 265); 

                auto tg_xxxxyyzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 266); 

                auto tg_xxxxyyzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 267); 

                auto tg_xxxxyyzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 268); 

                auto tg_xxxxyyzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 269); 

                auto tg_xxxxyyzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 270); 

                auto tg_xxxxyyzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 271); 

                auto tg_xxxxyyzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 272); 

                auto tg_xxxxyzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 273); 

                auto tg_xxxxyzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 274); 

                auto tg_xxxxyzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 275); 

                auto tg_xxxxyzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 276); 

                auto tg_xxxxyzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 277); 

                auto tg_xxxxyzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 278); 

                auto tg_xxxxyzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 279); 

                auto tg_xxxxyzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 280); 

                auto tg_xxxxyzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 281); 

                auto tg_xxxxyzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 282); 

                auto tg_xxxxyzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 283); 

                auto tg_xxxxyzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 284); 

                // Batch of Integrals (190,285)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxzzz_xxxxy_0, tg_xxxxxzzz_xxxxz_0, \
                                         tg_xxxxxzzz_xxxyy_0, tg_xxxxxzzz_xxxyz_0, tg_xxxxxzzz_xxxzz_0, tg_xxxxxzzz_xxyyy_0, \
                                         tg_xxxxxzzz_xxyyz_0, tg_xxxxxzzz_xxyzz_0, tg_xxxxxzzz_xxzzz_0, tg_xxxxxzzz_xyyyy_0, \
                                         tg_xxxxxzzz_xyyyz_0, tg_xxxxxzzz_xyyzz_0, tg_xxxxxzzz_xyzzz_0, tg_xxxxxzzz_xzzzz_0, \
                                         tg_xxxxxzzz_yyyyy_0, tg_xxxxxzzz_yyyyz_0, tg_xxxxxzzz_yyyzz_0, tg_xxxxxzzz_yyzzz_0, \
                                         tg_xxxxxzzz_yzzzz_0, tg_xxxxxzzz_zzzzz_0, tg_xxxxyyyy_xxxxx_0, tg_xxxxyyyy_xxxxy_0, \
                                         tg_xxxxyyyy_xxxxz_0, tg_xxxxyyyy_xxxyy_0, tg_xxxxyyyy_xxxyz_0, tg_xxxxyyyy_xxxzz_0, \
                                         tg_xxxxyyyy_xxyyy_0, tg_xxxxyyyy_xxyyz_0, tg_xxxxyyyy_xxyzz_0, tg_xxxxyyyy_xxzzz_0, \
                                         tg_xxxxyyyy_xyyyy_0, tg_xxxxyyyy_xyyyz_0, tg_xxxxyyyy_xyyzz_0, tg_xxxxyyyy_xyzzz_0, \
                                         tg_xxxxyyyy_xzzzz_0, tg_xxxxyyyy_yyyyy_0, tg_xxxxyyyy_yyyyz_0, tg_xxxxyyyy_yyyzz_0, \
                                         tg_xxxxyyyy_yyzzz_0, tg_xxxxyyyy_yzzzz_0, tg_xxxxyyyy_zzzzz_0, tg_xxxxyyyz_xxxxx_0, \
                                         tg_xxxxyyyz_xxxxy_0, tg_xxxxyyyz_xxxxz_0, tg_xxxxyyyz_xxxyy_0, tg_xxxxyyyz_xxxyz_0, \
                                         tg_xxxxyyyz_xxxzz_0, tg_xxxxyyyz_xxyyy_0, tg_xxxxyyyz_xxyyz_0, tg_xxxxyyyz_xxyzz_0, \
                                         tg_xxxxyyyz_xxzzz_0, tg_xxxxyyyz_xyyyy_0, tg_xxxxyyyz_xyyyz_0, tg_xxxxyyyz_xyyzz_0, \
                                         tg_xxxxyyyz_xyzzz_0, tg_xxxxyyyz_xzzzz_0, tg_xxxxyyyz_yyyyy_0, tg_xxxxyyyz_yyyyz_0, \
                                         tg_xxxxyyyz_yyyzz_0, tg_xxxxyyyz_yyzzz_0, tg_xxxxyyyz_yzzzz_0, tg_xxxxyyyz_zzzzz_0, \
                                         tg_xxxxyyzz_xxxxx_0, tg_xxxxyyzz_xxxxy_0, tg_xxxxyyzz_xxxxz_0, tg_xxxxyyzz_xxxyy_0, \
                                         tg_xxxxyyzz_xxxyz_0, tg_xxxxyyzz_xxxzz_0, tg_xxxxyyzz_xxyyy_0, tg_xxxxyyzz_xxyyz_0, \
                                         tg_xxxxyyzz_xxyzz_0, tg_xxxxyyzz_xxzzz_0, tg_xxxxyyzz_xyyyy_0, tg_xxxxyyzz_xyyyz_0, \
                                         tg_xxxxyyzz_xyyzz_0, tg_xxxxyyzz_xyzzz_0, tg_xxxxyyzz_xzzzz_0, tg_xxxxyyzz_yyyyy_0, \
                                         tg_xxxxyyzz_yyyyz_0, tg_xxxxyyzz_yyyzz_0, tg_xxxxyyzz_yyzzz_0, tg_xxxxyyzz_yzzzz_0, \
                                         tg_xxxxyyzz_zzzzz_0, tg_xxxxyzzz_xxxxx_0, tg_xxxxyzzz_xxxxy_0, tg_xxxxyzzz_xxxxz_0, \
                                         tg_xxxxyzzz_xxxyy_0, tg_xxxxyzzz_xxxyz_0, tg_xxxxyzzz_xxxzz_0, tg_xxxxyzzz_xxyyy_0, \
                                         tg_xxxxyzzz_xxyyz_0, tg_xxxxyzzz_xxyzz_0, tg_xxxxyzzz_xxzzz_0, tg_xxxxyzzz_xyyyy_0, \
                                         tg_xxxxyzzz_xyyyz_0, tg_xxxxzzz_xxxxy_0, tg_xxxxzzz_xxxxy_1, tg_xxxxzzz_xxxxz_0, \
                                         tg_xxxxzzz_xxxxz_1, tg_xxxxzzz_xxxy_1, tg_xxxxzzz_xxxyy_0, tg_xxxxzzz_xxxyy_1, \
                                         tg_xxxxzzz_xxxyz_0, tg_xxxxzzz_xxxyz_1, tg_xxxxzzz_xxxz_1, tg_xxxxzzz_xxxzz_0, \
                                         tg_xxxxzzz_xxxzz_1, tg_xxxxzzz_xxyy_1, tg_xxxxzzz_xxyyy_0, tg_xxxxzzz_xxyyy_1, \
                                         tg_xxxxzzz_xxyyz_0, tg_xxxxzzz_xxyyz_1, tg_xxxxzzz_xxyz_1, tg_xxxxzzz_xxyzz_0, \
                                         tg_xxxxzzz_xxyzz_1, tg_xxxxzzz_xxzz_1, tg_xxxxzzz_xxzzz_0, tg_xxxxzzz_xxzzz_1, \
                                         tg_xxxxzzz_xyyy_1, tg_xxxxzzz_xyyyy_0, tg_xxxxzzz_xyyyy_1, tg_xxxxzzz_xyyyz_0, \
                                         tg_xxxxzzz_xyyyz_1, tg_xxxxzzz_xyyz_1, tg_xxxxzzz_xyyzz_0, tg_xxxxzzz_xyyzz_1, \
                                         tg_xxxxzzz_xyzz_1, tg_xxxxzzz_xyzzz_0, tg_xxxxzzz_xyzzz_1, tg_xxxxzzz_xzzz_1, \
                                         tg_xxxxzzz_xzzzz_0, tg_xxxxzzz_xzzzz_1, tg_xxxxzzz_yyyy_1, tg_xxxxzzz_yyyyy_0, \
                                         tg_xxxxzzz_yyyyy_1, tg_xxxxzzz_yyyyz_0, tg_xxxxzzz_yyyyz_1, tg_xxxxzzz_yyyz_1, \
                                         tg_xxxxzzz_yyyzz_0, tg_xxxxzzz_yyyzz_1, tg_xxxxzzz_yyzz_1, tg_xxxxzzz_yyzzz_0, \
                                         tg_xxxxzzz_yyzzz_1, tg_xxxxzzz_yzzz_1, tg_xxxxzzz_yzzzz_0, tg_xxxxzzz_yzzzz_1, \
                                         tg_xxxxzzz_zzzz_1, tg_xxxxzzz_zzzzz_0, tg_xxxxzzz_zzzzz_1, tg_xxxyyyy_xxxx_1, \
                                         tg_xxxyyyy_xxxxx_0, tg_xxxyyyy_xxxxx_1, tg_xxxyyyy_xxxxy_0, tg_xxxyyyy_xxxxy_1, \
                                         tg_xxxyyyy_xxxxz_0, tg_xxxyyyy_xxxxz_1, tg_xxxyyyy_xxxy_1, tg_xxxyyyy_xxxyy_0, \
                                         tg_xxxyyyy_xxxyy_1, tg_xxxyyyy_xxxyz_0, tg_xxxyyyy_xxxyz_1, tg_xxxyyyy_xxxz_1, \
                                         tg_xxxyyyy_xxxzz_0, tg_xxxyyyy_xxxzz_1, tg_xxxyyyy_xxyy_1, tg_xxxyyyy_xxyyy_0, \
                                         tg_xxxyyyy_xxyyy_1, tg_xxxyyyy_xxyyz_0, tg_xxxyyyy_xxyyz_1, tg_xxxyyyy_xxyz_1, \
                                         tg_xxxyyyy_xxyzz_0, tg_xxxyyyy_xxyzz_1, tg_xxxyyyy_xxzz_1, tg_xxxyyyy_xxzzz_0, \
                                         tg_xxxyyyy_xxzzz_1, tg_xxxyyyy_xyyy_1, tg_xxxyyyy_xyyyy_0, tg_xxxyyyy_xyyyy_1, \
                                         tg_xxxyyyy_xyyyz_0, tg_xxxyyyy_xyyyz_1, tg_xxxyyyy_xyyz_1, tg_xxxyyyy_xyyzz_0, \
                                         tg_xxxyyyy_xyyzz_1, tg_xxxyyyy_xyzz_1, tg_xxxyyyy_xyzzz_0, tg_xxxyyyy_xyzzz_1, \
                                         tg_xxxyyyy_xzzz_1, tg_xxxyyyy_xzzzz_0, tg_xxxyyyy_xzzzz_1, tg_xxxyyyy_yyyy_1, \
                                         tg_xxxyyyy_yyyyy_0, tg_xxxyyyy_yyyyy_1, tg_xxxyyyy_yyyyz_0, tg_xxxyyyy_yyyyz_1, \
                                         tg_xxxyyyy_yyyz_1, tg_xxxyyyy_yyyzz_0, tg_xxxyyyy_yyyzz_1, tg_xxxyyyy_yyzz_1, \
                                         tg_xxxyyyy_yyzzz_0, tg_xxxyyyy_yyzzz_1, tg_xxxyyyy_yzzz_1, tg_xxxyyyy_yzzzz_0, \
                                         tg_xxxyyyy_yzzzz_1, tg_xxxyyyy_zzzz_1, tg_xxxyyyy_zzzzz_0, tg_xxxyyyy_zzzzz_1, \
                                         tg_xxxyyyz_xxxx_1, tg_xxxyyyz_xxxxx_0, tg_xxxyyyz_xxxxx_1, tg_xxxyyyz_xxxxy_0, \
                                         tg_xxxyyyz_xxxxy_1, tg_xxxyyyz_xxxxz_0, tg_xxxyyyz_xxxxz_1, tg_xxxyyyz_xxxy_1, \
                                         tg_xxxyyyz_xxxyy_0, tg_xxxyyyz_xxxyy_1, tg_xxxyyyz_xxxyz_0, tg_xxxyyyz_xxxyz_1, \
                                         tg_xxxyyyz_xxxz_1, tg_xxxyyyz_xxxzz_0, tg_xxxyyyz_xxxzz_1, tg_xxxyyyz_xxyy_1, \
                                         tg_xxxyyyz_xxyyy_0, tg_xxxyyyz_xxyyy_1, tg_xxxyyyz_xxyyz_0, tg_xxxyyyz_xxyyz_1, \
                                         tg_xxxyyyz_xxyz_1, tg_xxxyyyz_xxyzz_0, tg_xxxyyyz_xxyzz_1, tg_xxxyyyz_xxzz_1, \
                                         tg_xxxyyyz_xxzzz_0, tg_xxxyyyz_xxzzz_1, tg_xxxyyyz_xyyy_1, tg_xxxyyyz_xyyyy_0, \
                                         tg_xxxyyyz_xyyyy_1, tg_xxxyyyz_xyyyz_0, tg_xxxyyyz_xyyyz_1, tg_xxxyyyz_xyyz_1, \
                                         tg_xxxyyyz_xyyzz_0, tg_xxxyyyz_xyyzz_1, tg_xxxyyyz_xyzz_1, tg_xxxyyyz_xyzzz_0, \
                                         tg_xxxyyyz_xyzzz_1, tg_xxxyyyz_xzzz_1, tg_xxxyyyz_xzzzz_0, tg_xxxyyyz_xzzzz_1, \
                                         tg_xxxyyyz_yyyy_1, tg_xxxyyyz_yyyyy_0, tg_xxxyyyz_yyyyy_1, tg_xxxyyyz_yyyyz_0, \
                                         tg_xxxyyyz_yyyyz_1, tg_xxxyyyz_yyyz_1, tg_xxxyyyz_yyyzz_0, tg_xxxyyyz_yyyzz_1, \
                                         tg_xxxyyyz_yyzz_1, tg_xxxyyyz_yyzzz_0, tg_xxxyyyz_yyzzz_1, tg_xxxyyyz_yzzz_1, \
                                         tg_xxxyyyz_yzzzz_0, tg_xxxyyyz_yzzzz_1, tg_xxxyyyz_zzzz_1, tg_xxxyyyz_zzzzz_0, \
                                         tg_xxxyyyz_zzzzz_1, tg_xxxyyzz_xxxx_1, tg_xxxyyzz_xxxxx_0, tg_xxxyyzz_xxxxx_1, \
                                         tg_xxxyyzz_xxxxy_0, tg_xxxyyzz_xxxxy_1, tg_xxxyyzz_xxxxz_0, tg_xxxyyzz_xxxxz_1, \
                                         tg_xxxyyzz_xxxy_1, tg_xxxyyzz_xxxyy_0, tg_xxxyyzz_xxxyy_1, tg_xxxyyzz_xxxyz_0, \
                                         tg_xxxyyzz_xxxyz_1, tg_xxxyyzz_xxxz_1, tg_xxxyyzz_xxxzz_0, tg_xxxyyzz_xxxzz_1, \
                                         tg_xxxyyzz_xxyy_1, tg_xxxyyzz_xxyyy_0, tg_xxxyyzz_xxyyy_1, tg_xxxyyzz_xxyyz_0, \
                                         tg_xxxyyzz_xxyyz_1, tg_xxxyyzz_xxyz_1, tg_xxxyyzz_xxyzz_0, tg_xxxyyzz_xxyzz_1, \
                                         tg_xxxyyzz_xxzz_1, tg_xxxyyzz_xxzzz_0, tg_xxxyyzz_xxzzz_1, tg_xxxyyzz_xyyy_1, \
                                         tg_xxxyyzz_xyyyy_0, tg_xxxyyzz_xyyyy_1, tg_xxxyyzz_xyyyz_0, tg_xxxyyzz_xyyyz_1, \
                                         tg_xxxyyzz_xyyz_1, tg_xxxyyzz_xyyzz_0, tg_xxxyyzz_xyyzz_1, tg_xxxyyzz_xyzz_1, \
                                         tg_xxxyyzz_xyzzz_0, tg_xxxyyzz_xyzzz_1, tg_xxxyyzz_xzzz_1, tg_xxxyyzz_xzzzz_0, \
                                         tg_xxxyyzz_xzzzz_1, tg_xxxyyzz_yyyy_1, tg_xxxyyzz_yyyyy_0, tg_xxxyyzz_yyyyy_1, \
                                         tg_xxxyyzz_yyyyz_0, tg_xxxyyzz_yyyyz_1, tg_xxxyyzz_yyyz_1, tg_xxxyyzz_yyyzz_0, \
                                         tg_xxxyyzz_yyyzz_1, tg_xxxyyzz_yyzz_1, tg_xxxyyzz_yyzzz_0, tg_xxxyyzz_yyzzz_1, \
                                         tg_xxxyyzz_yzzz_1, tg_xxxyyzz_yzzzz_0, tg_xxxyyzz_yzzzz_1, tg_xxxyyzz_zzzz_1, \
                                         tg_xxxyyzz_zzzzz_0, tg_xxxyyzz_zzzzz_1, tg_xxxyzzz_xxxx_1, tg_xxxyzzz_xxxxx_0, \
                                         tg_xxxyzzz_xxxxx_1, tg_xxxyzzz_xxxxy_0, tg_xxxyzzz_xxxxy_1, tg_xxxyzzz_xxxxz_0, \
                                         tg_xxxyzzz_xxxxz_1, tg_xxxyzzz_xxxy_1, tg_xxxyzzz_xxxyy_0, tg_xxxyzzz_xxxyy_1, \
                                         tg_xxxyzzz_xxxyz_0, tg_xxxyzzz_xxxyz_1, tg_xxxyzzz_xxxz_1, tg_xxxyzzz_xxxzz_0, \
                                         tg_xxxyzzz_xxxzz_1, tg_xxxyzzz_xxyy_1, tg_xxxyzzz_xxyyy_0, tg_xxxyzzz_xxyyy_1, \
                                         tg_xxxyzzz_xxyyz_0, tg_xxxyzzz_xxyyz_1, tg_xxxyzzz_xxyz_1, tg_xxxyzzz_xxyzz_0, \
                                         tg_xxxyzzz_xxyzz_1, tg_xxxyzzz_xxzz_1, tg_xxxyzzz_xxzzz_0, tg_xxxyzzz_xxzzz_1, \
                                         tg_xxxyzzz_xyyy_1, tg_xxxyzzz_xyyyy_0, tg_xxxyzzz_xyyyy_1, tg_xxxyzzz_xyyyz_0, \
                                         tg_xxxyzzz_xyyyz_1, tg_xxxyzzz_xyyz_1, tg_xxxyzzz_xyzz_1, tg_xxxyzzz_xzzz_1, \
                                         tg_xxxyzzz_yyyy_1, tg_xxxyzzz_yyyz_1, tg_xxxzzz_xxxxy_0, tg_xxxzzz_xxxxy_1, \
                                         tg_xxxzzz_xxxxz_0, tg_xxxzzz_xxxxz_1, tg_xxxzzz_xxxyy_0, tg_xxxzzz_xxxyy_1, \
                                         tg_xxxzzz_xxxyz_0, tg_xxxzzz_xxxyz_1, tg_xxxzzz_xxxzz_0, tg_xxxzzz_xxxzz_1, \
                                         tg_xxxzzz_xxyyy_0, tg_xxxzzz_xxyyy_1, tg_xxxzzz_xxyyz_0, tg_xxxzzz_xxyyz_1, \
                                         tg_xxxzzz_xxyzz_0, tg_xxxzzz_xxyzz_1, tg_xxxzzz_xxzzz_0, tg_xxxzzz_xxzzz_1, \
                                         tg_xxxzzz_xyyyy_0, tg_xxxzzz_xyyyy_1, tg_xxxzzz_xyyyz_0, tg_xxxzzz_xyyyz_1, \
                                         tg_xxxzzz_xyyzz_0, tg_xxxzzz_xyyzz_1, tg_xxxzzz_xyzzz_0, tg_xxxzzz_xyzzz_1, \
                                         tg_xxxzzz_xzzzz_0, tg_xxxzzz_xzzzz_1, tg_xxxzzz_yyyyy_0, tg_xxxzzz_yyyyy_1, \
                                         tg_xxxzzz_yyyyz_0, tg_xxxzzz_yyyyz_1, tg_xxxzzz_yyyzz_0, tg_xxxzzz_yyyzz_1, \
                                         tg_xxxzzz_yyzzz_0, tg_xxxzzz_yyzzz_1, tg_xxxzzz_yzzzz_0, tg_xxxzzz_yzzzz_1, \
                                         tg_xxxzzz_zzzzz_0, tg_xxxzzz_zzzzz_1, tg_xxyyyy_xxxxx_0, tg_xxyyyy_xxxxx_1, \
                                         tg_xxyyyy_xxxxy_0, tg_xxyyyy_xxxxy_1, tg_xxyyyy_xxxxz_0, tg_xxyyyy_xxxxz_1, \
                                         tg_xxyyyy_xxxyy_0, tg_xxyyyy_xxxyy_1, tg_xxyyyy_xxxyz_0, tg_xxyyyy_xxxyz_1, \
                                         tg_xxyyyy_xxxzz_0, tg_xxyyyy_xxxzz_1, tg_xxyyyy_xxyyy_0, tg_xxyyyy_xxyyy_1, \
                                         tg_xxyyyy_xxyyz_0, tg_xxyyyy_xxyyz_1, tg_xxyyyy_xxyzz_0, tg_xxyyyy_xxyzz_1, \
                                         tg_xxyyyy_xxzzz_0, tg_xxyyyy_xxzzz_1, tg_xxyyyy_xyyyy_0, tg_xxyyyy_xyyyy_1, \
                                         tg_xxyyyy_xyyyz_0, tg_xxyyyy_xyyyz_1, tg_xxyyyy_xyyzz_0, tg_xxyyyy_xyyzz_1, \
                                         tg_xxyyyy_xyzzz_0, tg_xxyyyy_xyzzz_1, tg_xxyyyy_xzzzz_0, tg_xxyyyy_xzzzz_1, \
                                         tg_xxyyyy_yyyyy_0, tg_xxyyyy_yyyyy_1, tg_xxyyyy_yyyyz_0, tg_xxyyyy_yyyyz_1, \
                                         tg_xxyyyy_yyyzz_0, tg_xxyyyy_yyyzz_1, tg_xxyyyy_yyzzz_0, tg_xxyyyy_yyzzz_1, \
                                         tg_xxyyyy_yzzzz_0, tg_xxyyyy_yzzzz_1, tg_xxyyyy_zzzzz_0, tg_xxyyyy_zzzzz_1, \
                                         tg_xxyyyz_xxxxx_0, tg_xxyyyz_xxxxx_1, tg_xxyyyz_xxxxy_0, tg_xxyyyz_xxxxy_1, \
                                         tg_xxyyyz_xxxxz_0, tg_xxyyyz_xxxxz_1, tg_xxyyyz_xxxyy_0, tg_xxyyyz_xxxyy_1, \
                                         tg_xxyyyz_xxxyz_0, tg_xxyyyz_xxxyz_1, tg_xxyyyz_xxxzz_0, tg_xxyyyz_xxxzz_1, \
                                         tg_xxyyyz_xxyyy_0, tg_xxyyyz_xxyyy_1, tg_xxyyyz_xxyyz_0, tg_xxyyyz_xxyyz_1, \
                                         tg_xxyyyz_xxyzz_0, tg_xxyyyz_xxyzz_1, tg_xxyyyz_xxzzz_0, tg_xxyyyz_xxzzz_1, \
                                         tg_xxyyyz_xyyyy_0, tg_xxyyyz_xyyyy_1, tg_xxyyyz_xyyyz_0, tg_xxyyyz_xyyyz_1, \
                                         tg_xxyyyz_xyyzz_0, tg_xxyyyz_xyyzz_1, tg_xxyyyz_xyzzz_0, tg_xxyyyz_xyzzz_1, \
                                         tg_xxyyyz_xzzzz_0, tg_xxyyyz_xzzzz_1, tg_xxyyyz_yyyyy_0, tg_xxyyyz_yyyyy_1, \
                                         tg_xxyyyz_yyyyz_0, tg_xxyyyz_yyyyz_1, tg_xxyyyz_yyyzz_0, tg_xxyyyz_yyyzz_1, \
                                         tg_xxyyyz_yyzzz_0, tg_xxyyyz_yyzzz_1, tg_xxyyyz_yzzzz_0, tg_xxyyyz_yzzzz_1, \
                                         tg_xxyyyz_zzzzz_0, tg_xxyyyz_zzzzz_1, tg_xxyyzz_xxxxx_0, tg_xxyyzz_xxxxx_1, \
                                         tg_xxyyzz_xxxxy_0, tg_xxyyzz_xxxxy_1, tg_xxyyzz_xxxxz_0, tg_xxyyzz_xxxxz_1, \
                                         tg_xxyyzz_xxxyy_0, tg_xxyyzz_xxxyy_1, tg_xxyyzz_xxxyz_0, tg_xxyyzz_xxxyz_1, \
                                         tg_xxyyzz_xxxzz_0, tg_xxyyzz_xxxzz_1, tg_xxyyzz_xxyyy_0, tg_xxyyzz_xxyyy_1, \
                                         tg_xxyyzz_xxyyz_0, tg_xxyyzz_xxyyz_1, tg_xxyyzz_xxyzz_0, tg_xxyyzz_xxyzz_1, \
                                         tg_xxyyzz_xxzzz_0, tg_xxyyzz_xxzzz_1, tg_xxyyzz_xyyyy_0, tg_xxyyzz_xyyyy_1, \
                                         tg_xxyyzz_xyyyz_0, tg_xxyyzz_xyyyz_1, tg_xxyyzz_xyyzz_0, tg_xxyyzz_xyyzz_1, \
                                         tg_xxyyzz_xyzzz_0, tg_xxyyzz_xyzzz_1, tg_xxyyzz_xzzzz_0, tg_xxyyzz_xzzzz_1, \
                                         tg_xxyyzz_yyyyy_0, tg_xxyyzz_yyyyy_1, tg_xxyyzz_yyyyz_0, tg_xxyyzz_yyyyz_1, \
                                         tg_xxyyzz_yyyzz_0, tg_xxyyzz_yyyzz_1, tg_xxyyzz_yyzzz_0, tg_xxyyzz_yyzzz_1, \
                                         tg_xxyyzz_yzzzz_0, tg_xxyyzz_yzzzz_1, tg_xxyyzz_zzzzz_0, tg_xxyyzz_zzzzz_1, \
                                         tg_xxyzzz_xxxxx_0, tg_xxyzzz_xxxxx_1, tg_xxyzzz_xxxxy_0, tg_xxyzzz_xxxxy_1, \
                                         tg_xxyzzz_xxxxz_0, tg_xxyzzz_xxxxz_1, tg_xxyzzz_xxxyy_0, tg_xxyzzz_xxxyy_1, \
                                         tg_xxyzzz_xxxyz_0, tg_xxyzzz_xxxyz_1, tg_xxyzzz_xxxzz_0, tg_xxyzzz_xxxzz_1, \
                                         tg_xxyzzz_xxyyy_0, tg_xxyzzz_xxyyy_1, tg_xxyzzz_xxyyz_0, tg_xxyzzz_xxyyz_1, \
                                         tg_xxyzzz_xxyzz_0, tg_xxyzzz_xxyzz_1, tg_xxyzzz_xxzzz_0, tg_xxyzzz_xxzzz_1, \
                                         tg_xxyzzz_xyyyy_0, tg_xxyzzz_xyyyy_1, tg_xxyzzz_xyyyz_0, tg_xxyzzz_xyyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxzzz_xxxxy_0[j] = pb_x * tg_xxxxzzz_xxxxy_0[j] + fr * tg_xxxxzzz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxxxy_0[j] - tg_xxxzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzzz_xxxy_1[j];

                    tg_xxxxxzzz_xxxxz_0[j] = pb_x * tg_xxxxzzz_xxxxz_0[j] + fr * tg_xxxxzzz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxxxz_0[j] - tg_xxxzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzzz_xxxz_1[j];

                    tg_xxxxxzzz_xxxyy_0[j] = pb_x * tg_xxxxzzz_xxxyy_0[j] + fr * tg_xxxxzzz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxxyy_0[j] - tg_xxxzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzzz_xxyy_1[j];

                    tg_xxxxxzzz_xxxyz_0[j] = pb_x * tg_xxxxzzz_xxxyz_0[j] + fr * tg_xxxxzzz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxxyz_0[j] - tg_xxxzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzzz_xxyz_1[j];

                    tg_xxxxxzzz_xxxzz_0[j] = pb_x * tg_xxxxzzz_xxxzz_0[j] + fr * tg_xxxxzzz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxxzz_0[j] - tg_xxxzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzzz_xxzz_1[j];

                    tg_xxxxxzzz_xxyyy_0[j] = pb_x * tg_xxxxzzz_xxyyy_0[j] + fr * tg_xxxxzzz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxyyy_0[j] - tg_xxxzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzzz_xyyy_1[j];

                    tg_xxxxxzzz_xxyyz_0[j] = pb_x * tg_xxxxzzz_xxyyz_0[j] + fr * tg_xxxxzzz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxyyz_0[j] - tg_xxxzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzzz_xyyz_1[j];

                    tg_xxxxxzzz_xxyzz_0[j] = pb_x * tg_xxxxzzz_xxyzz_0[j] + fr * tg_xxxxzzz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxyzz_0[j] - tg_xxxzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzzz_xyzz_1[j];

                    tg_xxxxxzzz_xxzzz_0[j] = pb_x * tg_xxxxzzz_xxzzz_0[j] + fr * tg_xxxxzzz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxzzz_0[j] - tg_xxxzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzzz_xzzz_1[j];

                    tg_xxxxxzzz_xyyyy_0[j] = pb_x * tg_xxxxzzz_xyyyy_0[j] + fr * tg_xxxxzzz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xyyyy_0[j] - tg_xxxzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_yyyy_1[j];

                    tg_xxxxxzzz_xyyyz_0[j] = pb_x * tg_xxxxzzz_xyyyz_0[j] + fr * tg_xxxxzzz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xyyyz_0[j] - tg_xxxzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_yyyz_1[j];

                    tg_xxxxxzzz_xyyzz_0[j] = pb_x * tg_xxxxzzz_xyyzz_0[j] + fr * tg_xxxxzzz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xyyzz_0[j] - tg_xxxzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_yyzz_1[j];

                    tg_xxxxxzzz_xyzzz_0[j] = pb_x * tg_xxxxzzz_xyzzz_0[j] + fr * tg_xxxxzzz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xyzzz_0[j] - tg_xxxzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_yzzz_1[j];

                    tg_xxxxxzzz_xzzzz_0[j] = pb_x * tg_xxxxzzz_xzzzz_0[j] + fr * tg_xxxxzzz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xzzzz_0[j] - tg_xxxzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_zzzz_1[j];

                    tg_xxxxxzzz_yyyyy_0[j] = pb_x * tg_xxxxzzz_yyyyy_0[j] + fr * tg_xxxxzzz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yyyyy_0[j] - tg_xxxzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxxzzz_yyyyz_0[j] = pb_x * tg_xxxxzzz_yyyyz_0[j] + fr * tg_xxxxzzz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yyyyz_0[j] - tg_xxxzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxxzzz_yyyzz_0[j] = pb_x * tg_xxxxzzz_yyyzz_0[j] + fr * tg_xxxxzzz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yyyzz_0[j] - tg_xxxzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxxzzz_yyzzz_0[j] = pb_x * tg_xxxxzzz_yyzzz_0[j] + fr * tg_xxxxzzz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yyzzz_0[j] - tg_xxxzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxxzzz_yzzzz_0[j] = pb_x * tg_xxxxzzz_yzzzz_0[j] + fr * tg_xxxxzzz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yzzzz_0[j] - tg_xxxzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxxzzz_zzzzz_0[j] = pb_x * tg_xxxxzzz_zzzzz_0[j] + fr * tg_xxxxzzz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_zzzzz_0[j] - tg_xxxzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyyyy_xxxxx_0[j] = pb_x * tg_xxxyyyy_xxxxx_0[j] + fr * tg_xxxyyyy_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxxxx_0[j] - tg_xxyyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyyy_xxxx_1[j];

                    tg_xxxxyyyy_xxxxy_0[j] = pb_x * tg_xxxyyyy_xxxxy_0[j] + fr * tg_xxxyyyy_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxxxy_0[j] - tg_xxyyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyyy_xxxy_1[j];

                    tg_xxxxyyyy_xxxxz_0[j] = pb_x * tg_xxxyyyy_xxxxz_0[j] + fr * tg_xxxyyyy_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxxxz_0[j] - tg_xxyyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyyy_xxxz_1[j];

                    tg_xxxxyyyy_xxxyy_0[j] = pb_x * tg_xxxyyyy_xxxyy_0[j] + fr * tg_xxxyyyy_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxxyy_0[j] - tg_xxyyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyy_xxyy_1[j];

                    tg_xxxxyyyy_xxxyz_0[j] = pb_x * tg_xxxyyyy_xxxyz_0[j] + fr * tg_xxxyyyy_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxxyz_0[j] - tg_xxyyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyy_xxyz_1[j];

                    tg_xxxxyyyy_xxxzz_0[j] = pb_x * tg_xxxyyyy_xxxzz_0[j] + fr * tg_xxxyyyy_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxxzz_0[j] - tg_xxyyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyy_xxzz_1[j];

                    tg_xxxxyyyy_xxyyy_0[j] = pb_x * tg_xxxyyyy_xxyyy_0[j] + fr * tg_xxxyyyy_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxyyy_0[j] - tg_xxyyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyy_xyyy_1[j];

                    tg_xxxxyyyy_xxyyz_0[j] = pb_x * tg_xxxyyyy_xxyyz_0[j] + fr * tg_xxxyyyy_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxyyz_0[j] - tg_xxyyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyy_xyyz_1[j];

                    tg_xxxxyyyy_xxyzz_0[j] = pb_x * tg_xxxyyyy_xxyzz_0[j] + fr * tg_xxxyyyy_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxyzz_0[j] - tg_xxyyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyy_xyzz_1[j];

                    tg_xxxxyyyy_xxzzz_0[j] = pb_x * tg_xxxyyyy_xxzzz_0[j] + fr * tg_xxxyyyy_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxzzz_0[j] - tg_xxyyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyy_xzzz_1[j];

                    tg_xxxxyyyy_xyyyy_0[j] = pb_x * tg_xxxyyyy_xyyyy_0[j] + fr * tg_xxxyyyy_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xyyyy_0[j] - tg_xxyyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_yyyy_1[j];

                    tg_xxxxyyyy_xyyyz_0[j] = pb_x * tg_xxxyyyy_xyyyz_0[j] + fr * tg_xxxyyyy_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xyyyz_0[j] - tg_xxyyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_yyyz_1[j];

                    tg_xxxxyyyy_xyyzz_0[j] = pb_x * tg_xxxyyyy_xyyzz_0[j] + fr * tg_xxxyyyy_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xyyzz_0[j] - tg_xxyyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_yyzz_1[j];

                    tg_xxxxyyyy_xyzzz_0[j] = pb_x * tg_xxxyyyy_xyzzz_0[j] + fr * tg_xxxyyyy_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xyzzz_0[j] - tg_xxyyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_yzzz_1[j];

                    tg_xxxxyyyy_xzzzz_0[j] = pb_x * tg_xxxyyyy_xzzzz_0[j] + fr * tg_xxxyyyy_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xzzzz_0[j] - tg_xxyyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_zzzz_1[j];

                    tg_xxxxyyyy_yyyyy_0[j] = pb_x * tg_xxxyyyy_yyyyy_0[j] + fr * tg_xxxyyyy_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yyyyy_0[j] - tg_xxyyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyyyy_yyyyz_0[j] = pb_x * tg_xxxyyyy_yyyyz_0[j] + fr * tg_xxxyyyy_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yyyyz_0[j] - tg_xxyyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyyyy_yyyzz_0[j] = pb_x * tg_xxxyyyy_yyyzz_0[j] + fr * tg_xxxyyyy_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yyyzz_0[j] - tg_xxyyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyyyy_yyzzz_0[j] = pb_x * tg_xxxyyyy_yyzzz_0[j] + fr * tg_xxxyyyy_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yyzzz_0[j] - tg_xxyyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyyyy_yzzzz_0[j] = pb_x * tg_xxxyyyy_yzzzz_0[j] + fr * tg_xxxyyyy_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yzzzz_0[j] - tg_xxyyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyyyy_zzzzz_0[j] = pb_x * tg_xxxyyyy_zzzzz_0[j] + fr * tg_xxxyyyy_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_zzzzz_0[j] - tg_xxyyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyyyz_xxxxx_0[j] = pb_x * tg_xxxyyyz_xxxxx_0[j] + fr * tg_xxxyyyz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxxxx_0[j] - tg_xxyyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyyz_xxxx_1[j];

                    tg_xxxxyyyz_xxxxy_0[j] = pb_x * tg_xxxyyyz_xxxxy_0[j] + fr * tg_xxxyyyz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxxxy_0[j] - tg_xxyyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyyz_xxxy_1[j];

                    tg_xxxxyyyz_xxxxz_0[j] = pb_x * tg_xxxyyyz_xxxxz_0[j] + fr * tg_xxxyyyz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxxxz_0[j] - tg_xxyyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyyz_xxxz_1[j];

                    tg_xxxxyyyz_xxxyy_0[j] = pb_x * tg_xxxyyyz_xxxyy_0[j] + fr * tg_xxxyyyz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxxyy_0[j] - tg_xxyyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyz_xxyy_1[j];

                    tg_xxxxyyyz_xxxyz_0[j] = pb_x * tg_xxxyyyz_xxxyz_0[j] + fr * tg_xxxyyyz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxxyz_0[j] - tg_xxyyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyz_xxyz_1[j];

                    tg_xxxxyyyz_xxxzz_0[j] = pb_x * tg_xxxyyyz_xxxzz_0[j] + fr * tg_xxxyyyz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxxzz_0[j] - tg_xxyyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyz_xxzz_1[j];

                    tg_xxxxyyyz_xxyyy_0[j] = pb_x * tg_xxxyyyz_xxyyy_0[j] + fr * tg_xxxyyyz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxyyy_0[j] - tg_xxyyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyz_xyyy_1[j];

                    tg_xxxxyyyz_xxyyz_0[j] = pb_x * tg_xxxyyyz_xxyyz_0[j] + fr * tg_xxxyyyz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxyyz_0[j] - tg_xxyyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyz_xyyz_1[j];

                    tg_xxxxyyyz_xxyzz_0[j] = pb_x * tg_xxxyyyz_xxyzz_0[j] + fr * tg_xxxyyyz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxyzz_0[j] - tg_xxyyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyz_xyzz_1[j];

                    tg_xxxxyyyz_xxzzz_0[j] = pb_x * tg_xxxyyyz_xxzzz_0[j] + fr * tg_xxxyyyz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxzzz_0[j] - tg_xxyyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyz_xzzz_1[j];

                    tg_xxxxyyyz_xyyyy_0[j] = pb_x * tg_xxxyyyz_xyyyy_0[j] + fr * tg_xxxyyyz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xyyyy_0[j] - tg_xxyyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_yyyy_1[j];

                    tg_xxxxyyyz_xyyyz_0[j] = pb_x * tg_xxxyyyz_xyyyz_0[j] + fr * tg_xxxyyyz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xyyyz_0[j] - tg_xxyyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_yyyz_1[j];

                    tg_xxxxyyyz_xyyzz_0[j] = pb_x * tg_xxxyyyz_xyyzz_0[j] + fr * tg_xxxyyyz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xyyzz_0[j] - tg_xxyyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_yyzz_1[j];

                    tg_xxxxyyyz_xyzzz_0[j] = pb_x * tg_xxxyyyz_xyzzz_0[j] + fr * tg_xxxyyyz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xyzzz_0[j] - tg_xxyyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_yzzz_1[j];

                    tg_xxxxyyyz_xzzzz_0[j] = pb_x * tg_xxxyyyz_xzzzz_0[j] + fr * tg_xxxyyyz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xzzzz_0[j] - tg_xxyyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_zzzz_1[j];

                    tg_xxxxyyyz_yyyyy_0[j] = pb_x * tg_xxxyyyz_yyyyy_0[j] + fr * tg_xxxyyyz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yyyyy_0[j] - tg_xxyyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyyyz_yyyyz_0[j] = pb_x * tg_xxxyyyz_yyyyz_0[j] + fr * tg_xxxyyyz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yyyyz_0[j] - tg_xxyyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyyyz_yyyzz_0[j] = pb_x * tg_xxxyyyz_yyyzz_0[j] + fr * tg_xxxyyyz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yyyzz_0[j] - tg_xxyyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyyyz_yyzzz_0[j] = pb_x * tg_xxxyyyz_yyzzz_0[j] + fr * tg_xxxyyyz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yyzzz_0[j] - tg_xxyyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyyyz_yzzzz_0[j] = pb_x * tg_xxxyyyz_yzzzz_0[j] + fr * tg_xxxyyyz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yzzzz_0[j] - tg_xxyyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyyyz_zzzzz_0[j] = pb_x * tg_xxxyyyz_zzzzz_0[j] + fr * tg_xxxyyyz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_zzzzz_0[j] - tg_xxyyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyyzz_xxxxx_0[j] = pb_x * tg_xxxyyzz_xxxxx_0[j] + fr * tg_xxxyyzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxxxx_0[j] - tg_xxyyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyzz_xxxx_1[j];

                    tg_xxxxyyzz_xxxxy_0[j] = pb_x * tg_xxxyyzz_xxxxy_0[j] + fr * tg_xxxyyzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxxxy_0[j] - tg_xxyyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyzz_xxxy_1[j];

                    tg_xxxxyyzz_xxxxz_0[j] = pb_x * tg_xxxyyzz_xxxxz_0[j] + fr * tg_xxxyyzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxxxz_0[j] - tg_xxyyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyzz_xxxz_1[j];

                    tg_xxxxyyzz_xxxyy_0[j] = pb_x * tg_xxxyyzz_xxxyy_0[j] + fr * tg_xxxyyzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxxyy_0[j] - tg_xxyyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyzz_xxyy_1[j];

                    tg_xxxxyyzz_xxxyz_0[j] = pb_x * tg_xxxyyzz_xxxyz_0[j] + fr * tg_xxxyyzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxxyz_0[j] - tg_xxyyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyzz_xxyz_1[j];

                    tg_xxxxyyzz_xxxzz_0[j] = pb_x * tg_xxxyyzz_xxxzz_0[j] + fr * tg_xxxyyzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxxzz_0[j] - tg_xxyyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyzz_xxzz_1[j];

                    tg_xxxxyyzz_xxyyy_0[j] = pb_x * tg_xxxyyzz_xxyyy_0[j] + fr * tg_xxxyyzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxyyy_0[j] - tg_xxyyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyzz_xyyy_1[j];

                    tg_xxxxyyzz_xxyyz_0[j] = pb_x * tg_xxxyyzz_xxyyz_0[j] + fr * tg_xxxyyzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxyyz_0[j] - tg_xxyyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyzz_xyyz_1[j];

                    tg_xxxxyyzz_xxyzz_0[j] = pb_x * tg_xxxyyzz_xxyzz_0[j] + fr * tg_xxxyyzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxyzz_0[j] - tg_xxyyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyzz_xyzz_1[j];

                    tg_xxxxyyzz_xxzzz_0[j] = pb_x * tg_xxxyyzz_xxzzz_0[j] + fr * tg_xxxyyzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxzzz_0[j] - tg_xxyyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyzz_xzzz_1[j];

                    tg_xxxxyyzz_xyyyy_0[j] = pb_x * tg_xxxyyzz_xyyyy_0[j] + fr * tg_xxxyyzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xyyyy_0[j] - tg_xxyyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_yyyy_1[j];

                    tg_xxxxyyzz_xyyyz_0[j] = pb_x * tg_xxxyyzz_xyyyz_0[j] + fr * tg_xxxyyzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xyyyz_0[j] - tg_xxyyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_yyyz_1[j];

                    tg_xxxxyyzz_xyyzz_0[j] = pb_x * tg_xxxyyzz_xyyzz_0[j] + fr * tg_xxxyyzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xyyzz_0[j] - tg_xxyyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_yyzz_1[j];

                    tg_xxxxyyzz_xyzzz_0[j] = pb_x * tg_xxxyyzz_xyzzz_0[j] + fr * tg_xxxyyzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xyzzz_0[j] - tg_xxyyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_yzzz_1[j];

                    tg_xxxxyyzz_xzzzz_0[j] = pb_x * tg_xxxyyzz_xzzzz_0[j] + fr * tg_xxxyyzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xzzzz_0[j] - tg_xxyyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_zzzz_1[j];

                    tg_xxxxyyzz_yyyyy_0[j] = pb_x * tg_xxxyyzz_yyyyy_0[j] + fr * tg_xxxyyzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yyyyy_0[j] - tg_xxyyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyyzz_yyyyz_0[j] = pb_x * tg_xxxyyzz_yyyyz_0[j] + fr * tg_xxxyyzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yyyyz_0[j] - tg_xxyyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyyzz_yyyzz_0[j] = pb_x * tg_xxxyyzz_yyyzz_0[j] + fr * tg_xxxyyzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yyyzz_0[j] - tg_xxyyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyyzz_yyzzz_0[j] = pb_x * tg_xxxyyzz_yyzzz_0[j] + fr * tg_xxxyyzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yyzzz_0[j] - tg_xxyyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyyzz_yzzzz_0[j] = pb_x * tg_xxxyyzz_yzzzz_0[j] + fr * tg_xxxyyzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yzzzz_0[j] - tg_xxyyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyyzz_zzzzz_0[j] = pb_x * tg_xxxyyzz_zzzzz_0[j] + fr * tg_xxxyyzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_zzzzz_0[j] - tg_xxyyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxyzzz_xxxxx_0[j] = pb_x * tg_xxxyzzz_xxxxx_0[j] + fr * tg_xxxyzzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxxxx_0[j] - tg_xxyzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyzzz_xxxx_1[j];

                    tg_xxxxyzzz_xxxxy_0[j] = pb_x * tg_xxxyzzz_xxxxy_0[j] + fr * tg_xxxyzzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxxxy_0[j] - tg_xxyzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzzz_xxxy_1[j];

                    tg_xxxxyzzz_xxxxz_0[j] = pb_x * tg_xxxyzzz_xxxxz_0[j] + fr * tg_xxxyzzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxxxz_0[j] - tg_xxyzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzzz_xxxz_1[j];

                    tg_xxxxyzzz_xxxyy_0[j] = pb_x * tg_xxxyzzz_xxxyy_0[j] + fr * tg_xxxyzzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxxyy_0[j] - tg_xxyzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzzz_xxyy_1[j];

                    tg_xxxxyzzz_xxxyz_0[j] = pb_x * tg_xxxyzzz_xxxyz_0[j] + fr * tg_xxxyzzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxxyz_0[j] - tg_xxyzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzzz_xxyz_1[j];

                    tg_xxxxyzzz_xxxzz_0[j] = pb_x * tg_xxxyzzz_xxxzz_0[j] + fr * tg_xxxyzzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxxzz_0[j] - tg_xxyzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzzz_xxzz_1[j];

                    tg_xxxxyzzz_xxyyy_0[j] = pb_x * tg_xxxyzzz_xxyyy_0[j] + fr * tg_xxxyzzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxyyy_0[j] - tg_xxyzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzzz_xyyy_1[j];

                    tg_xxxxyzzz_xxyyz_0[j] = pb_x * tg_xxxyzzz_xxyyz_0[j] + fr * tg_xxxyzzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxyyz_0[j] - tg_xxyzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzzz_xyyz_1[j];

                    tg_xxxxyzzz_xxyzz_0[j] = pb_x * tg_xxxyzzz_xxyzz_0[j] + fr * tg_xxxyzzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxyzz_0[j] - tg_xxyzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzzz_xyzz_1[j];

                    tg_xxxxyzzz_xxzzz_0[j] = pb_x * tg_xxxyzzz_xxzzz_0[j] + fr * tg_xxxyzzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxzzz_0[j] - tg_xxyzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzzz_xzzz_1[j];

                    tg_xxxxyzzz_xyyyy_0[j] = pb_x * tg_xxxyzzz_xyyyy_0[j] + fr * tg_xxxyzzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xyyyy_0[j] - tg_xxyzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_yyyy_1[j];

                    tg_xxxxyzzz_xyyyz_0[j] = pb_x * tg_xxxyzzz_xyyyz_0[j] + fr * tg_xxxyzzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xyyyz_0[j] - tg_xxyzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_285_380(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxxyzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 285); 

                auto tg_xxxyzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 286); 

                auto tg_xxxyzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 287); 

                auto tg_xxxyzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 288); 

                auto tg_xxxyzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 289); 

                auto tg_xxxyzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 290); 

                auto tg_xxxyzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 291); 

                auto tg_xxxyzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 292); 

                auto tg_xxxyzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 293); 

                auto tg_xxxzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 294); 

                auto tg_xxxzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 295); 

                auto tg_xxxzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 296); 

                auto tg_xxxzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 297); 

                auto tg_xxxzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 298); 

                auto tg_xxxzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 299); 

                auto tg_xxxzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 300); 

                auto tg_xxxzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 301); 

                auto tg_xxxzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 302); 

                auto tg_xxxzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 303); 

                auto tg_xxxzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 304); 

                auto tg_xxxzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 305); 

                auto tg_xxxzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 306); 

                auto tg_xxxzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 307); 

                auto tg_xxxzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 308); 

                auto tg_xxxzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 309); 

                auto tg_xxxzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 310); 

                auto tg_xxxzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 311); 

                auto tg_xxxzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 312); 

                auto tg_xxxzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 313); 

                auto tg_xxxzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 314); 

                auto tg_xxyyyyy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 315); 

                auto tg_xxyyyyy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 316); 

                auto tg_xxyyyyy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 317); 

                auto tg_xxyyyyy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 318); 

                auto tg_xxyyyyy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 319); 

                auto tg_xxyyyyy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 320); 

                auto tg_xxyyyyy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 321); 

                auto tg_xxyyyyy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 322); 

                auto tg_xxyyyyy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 323); 

                auto tg_xxyyyyy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 324); 

                auto tg_xxyyyyy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 325); 

                auto tg_xxyyyyy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 326); 

                auto tg_xxyyyyy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 327); 

                auto tg_xxyyyyy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 328); 

                auto tg_xxyyyyy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 329); 

                auto tg_xxyyyyy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 330); 

                auto tg_xxyyyyy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 331); 

                auto tg_xxyyyyy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 332); 

                auto tg_xxyyyyy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 333); 

                auto tg_xxyyyyy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 334); 

                auto tg_xxyyyyy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 335); 

                auto tg_xxyyyyz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 336); 

                auto tg_xxyyyyz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 337); 

                auto tg_xxyyyyz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 338); 

                auto tg_xxyyyyz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 339); 

                auto tg_xxyyyyz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 340); 

                auto tg_xxyyyyz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 341); 

                auto tg_xxyyyyz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 342); 

                auto tg_xxyyyyz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 343); 

                auto tg_xxyyyyz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 344); 

                auto tg_xxyyyyz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 345); 

                auto tg_xxyyyyz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 346); 

                auto tg_xxyyyyz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 347); 

                auto tg_xxyyyyz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 348); 

                auto tg_xxyyyyz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 349); 

                auto tg_xxyyyyz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 350); 

                auto tg_xxyyyyz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 351); 

                auto tg_xxyyyyz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 352); 

                auto tg_xxyyyyz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 353); 

                auto tg_xxyyyyz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 354); 

                auto tg_xxyyyyz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 355); 

                auto tg_xxyyyyz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 356); 

                auto tg_xxyyyzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 357); 

                auto tg_xxyyyzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 358); 

                auto tg_xxyyyzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 359); 

                auto tg_xxyyyzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 360); 

                auto tg_xxyyyzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 361); 

                auto tg_xxyyyzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 362); 

                auto tg_xxyyyzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 363); 

                auto tg_xxyyyzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 364); 

                auto tg_xxyyyzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 365); 

                auto tg_xxyyyzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 366); 

                auto tg_xxyyyzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 367); 

                auto tg_xxyyyzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 368); 

                auto tg_xxyyyzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 369); 

                auto tg_xxyyyzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 370); 

                auto tg_xxyyyzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 371); 

                auto tg_xxyyyzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 372); 

                auto tg_xxyyyzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 373); 

                auto tg_xxyyyzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 374); 

                auto tg_xxyyyzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 375); 

                auto tg_xxyyyzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 376); 

                auto tg_xxyyyzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 377); 

                auto tg_xxyyzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 378); 

                auto tg_xxyyzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 379); 

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

                auto tg_xxxyzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 207); 

                auto tg_xxxyzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 208); 

                auto tg_xxxyzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 209); 

                auto tg_xxxzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 210); 

                auto tg_xxxzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 211); 

                auto tg_xxxzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 212); 

                auto tg_xxxzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 213); 

                auto tg_xxxzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 214); 

                auto tg_xxxzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 215); 

                auto tg_xxxzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 216); 

                auto tg_xxxzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 217); 

                auto tg_xxxzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 218); 

                auto tg_xxxzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 219); 

                auto tg_xxxzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 220); 

                auto tg_xxxzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 221); 

                auto tg_xxxzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 222); 

                auto tg_xxxzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 223); 

                auto tg_xxxzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 224); 

                auto tg_xxyyyyy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 225); 

                auto tg_xxyyyyy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 226); 

                auto tg_xxyyyyy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 227); 

                auto tg_xxyyyyy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 228); 

                auto tg_xxyyyyy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 229); 

                auto tg_xxyyyyy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 230); 

                auto tg_xxyyyyy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 231); 

                auto tg_xxyyyyy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 232); 

                auto tg_xxyyyyy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 233); 

                auto tg_xxyyyyy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 234); 

                auto tg_xxyyyyy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 235); 

                auto tg_xxyyyyy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 236); 

                auto tg_xxyyyyy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 237); 

                auto tg_xxyyyyy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 238); 

                auto tg_xxyyyyy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 239); 

                auto tg_xxyyyyz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 240); 

                auto tg_xxyyyyz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 241); 

                auto tg_xxyyyyz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 242); 

                auto tg_xxyyyyz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 243); 

                auto tg_xxyyyyz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 244); 

                auto tg_xxyyyyz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 245); 

                auto tg_xxyyyyz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 246); 

                auto tg_xxyyyyz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 247); 

                auto tg_xxyyyyz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 248); 

                auto tg_xxyyyyz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 249); 

                auto tg_xxyyyyz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 250); 

                auto tg_xxyyyyz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 251); 

                auto tg_xxyyyyz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 252); 

                auto tg_xxyyyyz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 253); 

                auto tg_xxyyyyz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 254); 

                auto tg_xxyyyzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 255); 

                auto tg_xxyyyzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 256); 

                auto tg_xxyyyzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 257); 

                auto tg_xxyyyzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 258); 

                auto tg_xxyyyzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 259); 

                auto tg_xxyyyzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 260); 

                auto tg_xxyyyzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 261); 

                auto tg_xxyyyzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 262); 

                auto tg_xxyyyzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 263); 

                auto tg_xxyyyzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 264); 

                auto tg_xxyyyzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 265); 

                auto tg_xxyyyzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 266); 

                auto tg_xxyyyzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 267); 

                auto tg_xxyyyzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 268); 

                auto tg_xxyyyzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 269); 

                auto tg_xxyyzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 270); 

                auto tg_xxyyzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 271); 

                // set up pointers to integrals

                auto tg_xxxxyzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 285); 

                auto tg_xxxxyzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 286); 

                auto tg_xxxxyzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 287); 

                auto tg_xxxxyzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 288); 

                auto tg_xxxxyzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 289); 

                auto tg_xxxxyzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 290); 

                auto tg_xxxxyzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 291); 

                auto tg_xxxxyzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 292); 

                auto tg_xxxxyzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 293); 

                auto tg_xxxxzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 294); 

                auto tg_xxxxzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 295); 

                auto tg_xxxxzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 296); 

                auto tg_xxxxzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 297); 

                auto tg_xxxxzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 298); 

                auto tg_xxxxzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 299); 

                auto tg_xxxxzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 300); 

                auto tg_xxxxzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 301); 

                auto tg_xxxxzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 302); 

                auto tg_xxxxzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 303); 

                auto tg_xxxxzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 304); 

                auto tg_xxxxzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 305); 

                auto tg_xxxxzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 306); 

                auto tg_xxxxzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 307); 

                auto tg_xxxxzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 308); 

                auto tg_xxxxzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 309); 

                auto tg_xxxxzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 310); 

                auto tg_xxxxzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 311); 

                auto tg_xxxxzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 312); 

                auto tg_xxxxzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 313); 

                auto tg_xxxxzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 314); 

                auto tg_xxxyyyyy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 315); 

                auto tg_xxxyyyyy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 316); 

                auto tg_xxxyyyyy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 317); 

                auto tg_xxxyyyyy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 318); 

                auto tg_xxxyyyyy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 319); 

                auto tg_xxxyyyyy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 320); 

                auto tg_xxxyyyyy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 321); 

                auto tg_xxxyyyyy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 322); 

                auto tg_xxxyyyyy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 323); 

                auto tg_xxxyyyyy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 324); 

                auto tg_xxxyyyyy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 325); 

                auto tg_xxxyyyyy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 326); 

                auto tg_xxxyyyyy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 327); 

                auto tg_xxxyyyyy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 328); 

                auto tg_xxxyyyyy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 329); 

                auto tg_xxxyyyyy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 330); 

                auto tg_xxxyyyyy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 331); 

                auto tg_xxxyyyyy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 332); 

                auto tg_xxxyyyyy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 333); 

                auto tg_xxxyyyyy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 334); 

                auto tg_xxxyyyyy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 335); 

                auto tg_xxxyyyyz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 336); 

                auto tg_xxxyyyyz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 337); 

                auto tg_xxxyyyyz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 338); 

                auto tg_xxxyyyyz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 339); 

                auto tg_xxxyyyyz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 340); 

                auto tg_xxxyyyyz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 341); 

                auto tg_xxxyyyyz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 342); 

                auto tg_xxxyyyyz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 343); 

                auto tg_xxxyyyyz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 344); 

                auto tg_xxxyyyyz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 345); 

                auto tg_xxxyyyyz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 346); 

                auto tg_xxxyyyyz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 347); 

                auto tg_xxxyyyyz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 348); 

                auto tg_xxxyyyyz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 349); 

                auto tg_xxxyyyyz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 350); 

                auto tg_xxxyyyyz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 351); 

                auto tg_xxxyyyyz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 352); 

                auto tg_xxxyyyyz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 353); 

                auto tg_xxxyyyyz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 354); 

                auto tg_xxxyyyyz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 355); 

                auto tg_xxxyyyyz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 356); 

                auto tg_xxxyyyzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 357); 

                auto tg_xxxyyyzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 358); 

                auto tg_xxxyyyzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 359); 

                auto tg_xxxyyyzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 360); 

                auto tg_xxxyyyzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 361); 

                auto tg_xxxyyyzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 362); 

                auto tg_xxxyyyzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 363); 

                auto tg_xxxyyyzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 364); 

                auto tg_xxxyyyzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 365); 

                auto tg_xxxyyyzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 366); 

                auto tg_xxxyyyzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 367); 

                auto tg_xxxyyyzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 368); 

                auto tg_xxxyyyzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 369); 

                auto tg_xxxyyyzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 370); 

                auto tg_xxxyyyzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 371); 

                auto tg_xxxyyyzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 372); 

                auto tg_xxxyyyzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 373); 

                auto tg_xxxyyyzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 374); 

                auto tg_xxxyyyzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 375); 

                auto tg_xxxyyyzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 376); 

                auto tg_xxxyyyzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 377); 

                auto tg_xxxyyzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 378); 

                auto tg_xxxyyzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 379); 

                // Batch of Integrals (285,380)

                #pragma omp simd aligned(fxn, fza, tg_xxxxyzzz_xyyzz_0, tg_xxxxyzzz_xyzzz_0, \
                                         tg_xxxxyzzz_xzzzz_0, tg_xxxxyzzz_yyyyy_0, tg_xxxxyzzz_yyyyz_0, tg_xxxxyzzz_yyyzz_0, \
                                         tg_xxxxyzzz_yyzzz_0, tg_xxxxyzzz_yzzzz_0, tg_xxxxyzzz_zzzzz_0, tg_xxxxzzzz_xxxxx_0, \
                                         tg_xxxxzzzz_xxxxy_0, tg_xxxxzzzz_xxxxz_0, tg_xxxxzzzz_xxxyy_0, tg_xxxxzzzz_xxxyz_0, \
                                         tg_xxxxzzzz_xxxzz_0, tg_xxxxzzzz_xxyyy_0, tg_xxxxzzzz_xxyyz_0, tg_xxxxzzzz_xxyzz_0, \
                                         tg_xxxxzzzz_xxzzz_0, tg_xxxxzzzz_xyyyy_0, tg_xxxxzzzz_xyyyz_0, tg_xxxxzzzz_xyyzz_0, \
                                         tg_xxxxzzzz_xyzzz_0, tg_xxxxzzzz_xzzzz_0, tg_xxxxzzzz_yyyyy_0, tg_xxxxzzzz_yyyyz_0, \
                                         tg_xxxxzzzz_yyyzz_0, tg_xxxxzzzz_yyzzz_0, tg_xxxxzzzz_yzzzz_0, tg_xxxxzzzz_zzzzz_0, \
                                         tg_xxxyyyyy_xxxxx_0, tg_xxxyyyyy_xxxxy_0, tg_xxxyyyyy_xxxxz_0, tg_xxxyyyyy_xxxyy_0, \
                                         tg_xxxyyyyy_xxxyz_0, tg_xxxyyyyy_xxxzz_0, tg_xxxyyyyy_xxyyy_0, tg_xxxyyyyy_xxyyz_0, \
                                         tg_xxxyyyyy_xxyzz_0, tg_xxxyyyyy_xxzzz_0, tg_xxxyyyyy_xyyyy_0, tg_xxxyyyyy_xyyyz_0, \
                                         tg_xxxyyyyy_xyyzz_0, tg_xxxyyyyy_xyzzz_0, tg_xxxyyyyy_xzzzz_0, tg_xxxyyyyy_yyyyy_0, \
                                         tg_xxxyyyyy_yyyyz_0, tg_xxxyyyyy_yyyzz_0, tg_xxxyyyyy_yyzzz_0, tg_xxxyyyyy_yzzzz_0, \
                                         tg_xxxyyyyy_zzzzz_0, tg_xxxyyyyz_xxxxx_0, tg_xxxyyyyz_xxxxy_0, tg_xxxyyyyz_xxxxz_0, \
                                         tg_xxxyyyyz_xxxyy_0, tg_xxxyyyyz_xxxyz_0, tg_xxxyyyyz_xxxzz_0, tg_xxxyyyyz_xxyyy_0, \
                                         tg_xxxyyyyz_xxyyz_0, tg_xxxyyyyz_xxyzz_0, tg_xxxyyyyz_xxzzz_0, tg_xxxyyyyz_xyyyy_0, \
                                         tg_xxxyyyyz_xyyyz_0, tg_xxxyyyyz_xyyzz_0, tg_xxxyyyyz_xyzzz_0, tg_xxxyyyyz_xzzzz_0, \
                                         tg_xxxyyyyz_yyyyy_0, tg_xxxyyyyz_yyyyz_0, tg_xxxyyyyz_yyyzz_0, tg_xxxyyyyz_yyzzz_0, \
                                         tg_xxxyyyyz_yzzzz_0, tg_xxxyyyyz_zzzzz_0, tg_xxxyyyzz_xxxxx_0, tg_xxxyyyzz_xxxxy_0, \
                                         tg_xxxyyyzz_xxxxz_0, tg_xxxyyyzz_xxxyy_0, tg_xxxyyyzz_xxxyz_0, tg_xxxyyyzz_xxxzz_0, \
                                         tg_xxxyyyzz_xxyyy_0, tg_xxxyyyzz_xxyyz_0, tg_xxxyyyzz_xxyzz_0, tg_xxxyyyzz_xxzzz_0, \
                                         tg_xxxyyyzz_xyyyy_0, tg_xxxyyyzz_xyyyz_0, tg_xxxyyyzz_xyyzz_0, tg_xxxyyyzz_xyzzz_0, \
                                         tg_xxxyyyzz_xzzzz_0, tg_xxxyyyzz_yyyyy_0, tg_xxxyyyzz_yyyyz_0, tg_xxxyyyzz_yyyzz_0, \
                                         tg_xxxyyyzz_yyzzz_0, tg_xxxyyyzz_yzzzz_0, tg_xxxyyyzz_zzzzz_0, tg_xxxyyzzz_xxxxx_0, \
                                         tg_xxxyyzzz_xxxxy_0, tg_xxxyzzz_xyyzz_0, tg_xxxyzzz_xyyzz_1, tg_xxxyzzz_xyzzz_0, \
                                         tg_xxxyzzz_xyzzz_1, tg_xxxyzzz_xzzzz_0, tg_xxxyzzz_xzzzz_1, tg_xxxyzzz_yyyyy_0, \
                                         tg_xxxyzzz_yyyyy_1, tg_xxxyzzz_yyyyz_0, tg_xxxyzzz_yyyyz_1, tg_xxxyzzz_yyyzz_0, \
                                         tg_xxxyzzz_yyyzz_1, tg_xxxyzzz_yyzz_1, tg_xxxyzzz_yyzzz_0, tg_xxxyzzz_yyzzz_1, \
                                         tg_xxxyzzz_yzzz_1, tg_xxxyzzz_yzzzz_0, tg_xxxyzzz_yzzzz_1, tg_xxxyzzz_zzzz_1, \
                                         tg_xxxyzzz_zzzzz_0, tg_xxxyzzz_zzzzz_1, tg_xxxzzzz_xxxx_1, tg_xxxzzzz_xxxxx_0, \
                                         tg_xxxzzzz_xxxxx_1, tg_xxxzzzz_xxxxy_0, tg_xxxzzzz_xxxxy_1, tg_xxxzzzz_xxxxz_0, \
                                         tg_xxxzzzz_xxxxz_1, tg_xxxzzzz_xxxy_1, tg_xxxzzzz_xxxyy_0, tg_xxxzzzz_xxxyy_1, \
                                         tg_xxxzzzz_xxxyz_0, tg_xxxzzzz_xxxyz_1, tg_xxxzzzz_xxxz_1, tg_xxxzzzz_xxxzz_0, \
                                         tg_xxxzzzz_xxxzz_1, tg_xxxzzzz_xxyy_1, tg_xxxzzzz_xxyyy_0, tg_xxxzzzz_xxyyy_1, \
                                         tg_xxxzzzz_xxyyz_0, tg_xxxzzzz_xxyyz_1, tg_xxxzzzz_xxyz_1, tg_xxxzzzz_xxyzz_0, \
                                         tg_xxxzzzz_xxyzz_1, tg_xxxzzzz_xxzz_1, tg_xxxzzzz_xxzzz_0, tg_xxxzzzz_xxzzz_1, \
                                         tg_xxxzzzz_xyyy_1, tg_xxxzzzz_xyyyy_0, tg_xxxzzzz_xyyyy_1, tg_xxxzzzz_xyyyz_0, \
                                         tg_xxxzzzz_xyyyz_1, tg_xxxzzzz_xyyz_1, tg_xxxzzzz_xyyzz_0, tg_xxxzzzz_xyyzz_1, \
                                         tg_xxxzzzz_xyzz_1, tg_xxxzzzz_xyzzz_0, tg_xxxzzzz_xyzzz_1, tg_xxxzzzz_xzzz_1, \
                                         tg_xxxzzzz_xzzzz_0, tg_xxxzzzz_xzzzz_1, tg_xxxzzzz_yyyy_1, tg_xxxzzzz_yyyyy_0, \
                                         tg_xxxzzzz_yyyyy_1, tg_xxxzzzz_yyyyz_0, tg_xxxzzzz_yyyyz_1, tg_xxxzzzz_yyyz_1, \
                                         tg_xxxzzzz_yyyzz_0, tg_xxxzzzz_yyyzz_1, tg_xxxzzzz_yyzz_1, tg_xxxzzzz_yyzzz_0, \
                                         tg_xxxzzzz_yyzzz_1, tg_xxxzzzz_yzzz_1, tg_xxxzzzz_yzzzz_0, tg_xxxzzzz_yzzzz_1, \
                                         tg_xxxzzzz_zzzz_1, tg_xxxzzzz_zzzzz_0, tg_xxxzzzz_zzzzz_1, tg_xxyyyyy_xxxx_1, \
                                         tg_xxyyyyy_xxxxx_0, tg_xxyyyyy_xxxxx_1, tg_xxyyyyy_xxxxy_0, tg_xxyyyyy_xxxxy_1, \
                                         tg_xxyyyyy_xxxxz_0, tg_xxyyyyy_xxxxz_1, tg_xxyyyyy_xxxy_1, tg_xxyyyyy_xxxyy_0, \
                                         tg_xxyyyyy_xxxyy_1, tg_xxyyyyy_xxxyz_0, tg_xxyyyyy_xxxyz_1, tg_xxyyyyy_xxxz_1, \
                                         tg_xxyyyyy_xxxzz_0, tg_xxyyyyy_xxxzz_1, tg_xxyyyyy_xxyy_1, tg_xxyyyyy_xxyyy_0, \
                                         tg_xxyyyyy_xxyyy_1, tg_xxyyyyy_xxyyz_0, tg_xxyyyyy_xxyyz_1, tg_xxyyyyy_xxyz_1, \
                                         tg_xxyyyyy_xxyzz_0, tg_xxyyyyy_xxyzz_1, tg_xxyyyyy_xxzz_1, tg_xxyyyyy_xxzzz_0, \
                                         tg_xxyyyyy_xxzzz_1, tg_xxyyyyy_xyyy_1, tg_xxyyyyy_xyyyy_0, tg_xxyyyyy_xyyyy_1, \
                                         tg_xxyyyyy_xyyyz_0, tg_xxyyyyy_xyyyz_1, tg_xxyyyyy_xyyz_1, tg_xxyyyyy_xyyzz_0, \
                                         tg_xxyyyyy_xyyzz_1, tg_xxyyyyy_xyzz_1, tg_xxyyyyy_xyzzz_0, tg_xxyyyyy_xyzzz_1, \
                                         tg_xxyyyyy_xzzz_1, tg_xxyyyyy_xzzzz_0, tg_xxyyyyy_xzzzz_1, tg_xxyyyyy_yyyy_1, \
                                         tg_xxyyyyy_yyyyy_0, tg_xxyyyyy_yyyyy_1, tg_xxyyyyy_yyyyz_0, tg_xxyyyyy_yyyyz_1, \
                                         tg_xxyyyyy_yyyz_1, tg_xxyyyyy_yyyzz_0, tg_xxyyyyy_yyyzz_1, tg_xxyyyyy_yyzz_1, \
                                         tg_xxyyyyy_yyzzz_0, tg_xxyyyyy_yyzzz_1, tg_xxyyyyy_yzzz_1, tg_xxyyyyy_yzzzz_0, \
                                         tg_xxyyyyy_yzzzz_1, tg_xxyyyyy_zzzz_1, tg_xxyyyyy_zzzzz_0, tg_xxyyyyy_zzzzz_1, \
                                         tg_xxyyyyz_xxxx_1, tg_xxyyyyz_xxxxx_0, tg_xxyyyyz_xxxxx_1, tg_xxyyyyz_xxxxy_0, \
                                         tg_xxyyyyz_xxxxy_1, tg_xxyyyyz_xxxxz_0, tg_xxyyyyz_xxxxz_1, tg_xxyyyyz_xxxy_1, \
                                         tg_xxyyyyz_xxxyy_0, tg_xxyyyyz_xxxyy_1, tg_xxyyyyz_xxxyz_0, tg_xxyyyyz_xxxyz_1, \
                                         tg_xxyyyyz_xxxz_1, tg_xxyyyyz_xxxzz_0, tg_xxyyyyz_xxxzz_1, tg_xxyyyyz_xxyy_1, \
                                         tg_xxyyyyz_xxyyy_0, tg_xxyyyyz_xxyyy_1, tg_xxyyyyz_xxyyz_0, tg_xxyyyyz_xxyyz_1, \
                                         tg_xxyyyyz_xxyz_1, tg_xxyyyyz_xxyzz_0, tg_xxyyyyz_xxyzz_1, tg_xxyyyyz_xxzz_1, \
                                         tg_xxyyyyz_xxzzz_0, tg_xxyyyyz_xxzzz_1, tg_xxyyyyz_xyyy_1, tg_xxyyyyz_xyyyy_0, \
                                         tg_xxyyyyz_xyyyy_1, tg_xxyyyyz_xyyyz_0, tg_xxyyyyz_xyyyz_1, tg_xxyyyyz_xyyz_1, \
                                         tg_xxyyyyz_xyyzz_0, tg_xxyyyyz_xyyzz_1, tg_xxyyyyz_xyzz_1, tg_xxyyyyz_xyzzz_0, \
                                         tg_xxyyyyz_xyzzz_1, tg_xxyyyyz_xzzz_1, tg_xxyyyyz_xzzzz_0, tg_xxyyyyz_xzzzz_1, \
                                         tg_xxyyyyz_yyyy_1, tg_xxyyyyz_yyyyy_0, tg_xxyyyyz_yyyyy_1, tg_xxyyyyz_yyyyz_0, \
                                         tg_xxyyyyz_yyyyz_1, tg_xxyyyyz_yyyz_1, tg_xxyyyyz_yyyzz_0, tg_xxyyyyz_yyyzz_1, \
                                         tg_xxyyyyz_yyzz_1, tg_xxyyyyz_yyzzz_0, tg_xxyyyyz_yyzzz_1, tg_xxyyyyz_yzzz_1, \
                                         tg_xxyyyyz_yzzzz_0, tg_xxyyyyz_yzzzz_1, tg_xxyyyyz_zzzz_1, tg_xxyyyyz_zzzzz_0, \
                                         tg_xxyyyyz_zzzzz_1, tg_xxyyyzz_xxxx_1, tg_xxyyyzz_xxxxx_0, tg_xxyyyzz_xxxxx_1, \
                                         tg_xxyyyzz_xxxxy_0, tg_xxyyyzz_xxxxy_1, tg_xxyyyzz_xxxxz_0, tg_xxyyyzz_xxxxz_1, \
                                         tg_xxyyyzz_xxxy_1, tg_xxyyyzz_xxxyy_0, tg_xxyyyzz_xxxyy_1, tg_xxyyyzz_xxxyz_0, \
                                         tg_xxyyyzz_xxxyz_1, tg_xxyyyzz_xxxz_1, tg_xxyyyzz_xxxzz_0, tg_xxyyyzz_xxxzz_1, \
                                         tg_xxyyyzz_xxyy_1, tg_xxyyyzz_xxyyy_0, tg_xxyyyzz_xxyyy_1, tg_xxyyyzz_xxyyz_0, \
                                         tg_xxyyyzz_xxyyz_1, tg_xxyyyzz_xxyz_1, tg_xxyyyzz_xxyzz_0, tg_xxyyyzz_xxyzz_1, \
                                         tg_xxyyyzz_xxzz_1, tg_xxyyyzz_xxzzz_0, tg_xxyyyzz_xxzzz_1, tg_xxyyyzz_xyyy_1, \
                                         tg_xxyyyzz_xyyyy_0, tg_xxyyyzz_xyyyy_1, tg_xxyyyzz_xyyyz_0, tg_xxyyyzz_xyyyz_1, \
                                         tg_xxyyyzz_xyyz_1, tg_xxyyyzz_xyyzz_0, tg_xxyyyzz_xyyzz_1, tg_xxyyyzz_xyzz_1, \
                                         tg_xxyyyzz_xyzzz_0, tg_xxyyyzz_xyzzz_1, tg_xxyyyzz_xzzz_1, tg_xxyyyzz_xzzzz_0, \
                                         tg_xxyyyzz_xzzzz_1, tg_xxyyyzz_yyyy_1, tg_xxyyyzz_yyyyy_0, tg_xxyyyzz_yyyyy_1, \
                                         tg_xxyyyzz_yyyyz_0, tg_xxyyyzz_yyyyz_1, tg_xxyyyzz_yyyz_1, tg_xxyyyzz_yyyzz_0, \
                                         tg_xxyyyzz_yyyzz_1, tg_xxyyyzz_yyzz_1, tg_xxyyyzz_yyzzz_0, tg_xxyyyzz_yyzzz_1, \
                                         tg_xxyyyzz_yzzz_1, tg_xxyyyzz_yzzzz_0, tg_xxyyyzz_yzzzz_1, tg_xxyyyzz_zzzz_1, \
                                         tg_xxyyyzz_zzzzz_0, tg_xxyyyzz_zzzzz_1, tg_xxyyzzz_xxxx_1, tg_xxyyzzz_xxxxx_0, \
                                         tg_xxyyzzz_xxxxx_1, tg_xxyyzzz_xxxxy_0, tg_xxyyzzz_xxxxy_1, tg_xxyyzzz_xxxy_1, \
                                         tg_xxyzzz_xyyzz_0, tg_xxyzzz_xyyzz_1, tg_xxyzzz_xyzzz_0, tg_xxyzzz_xyzzz_1, \
                                         tg_xxyzzz_xzzzz_0, tg_xxyzzz_xzzzz_1, tg_xxyzzz_yyyyy_0, tg_xxyzzz_yyyyy_1, \
                                         tg_xxyzzz_yyyyz_0, tg_xxyzzz_yyyyz_1, tg_xxyzzz_yyyzz_0, tg_xxyzzz_yyyzz_1, \
                                         tg_xxyzzz_yyzzz_0, tg_xxyzzz_yyzzz_1, tg_xxyzzz_yzzzz_0, tg_xxyzzz_yzzzz_1, \
                                         tg_xxyzzz_zzzzz_0, tg_xxyzzz_zzzzz_1, tg_xxzzzz_xxxxx_0, tg_xxzzzz_xxxxx_1, \
                                         tg_xxzzzz_xxxxy_0, tg_xxzzzz_xxxxy_1, tg_xxzzzz_xxxxz_0, tg_xxzzzz_xxxxz_1, \
                                         tg_xxzzzz_xxxyy_0, tg_xxzzzz_xxxyy_1, tg_xxzzzz_xxxyz_0, tg_xxzzzz_xxxyz_1, \
                                         tg_xxzzzz_xxxzz_0, tg_xxzzzz_xxxzz_1, tg_xxzzzz_xxyyy_0, tg_xxzzzz_xxyyy_1, \
                                         tg_xxzzzz_xxyyz_0, tg_xxzzzz_xxyyz_1, tg_xxzzzz_xxyzz_0, tg_xxzzzz_xxyzz_1, \
                                         tg_xxzzzz_xxzzz_0, tg_xxzzzz_xxzzz_1, tg_xxzzzz_xyyyy_0, tg_xxzzzz_xyyyy_1, \
                                         tg_xxzzzz_xyyyz_0, tg_xxzzzz_xyyyz_1, tg_xxzzzz_xyyzz_0, tg_xxzzzz_xyyzz_1, \
                                         tg_xxzzzz_xyzzz_0, tg_xxzzzz_xyzzz_1, tg_xxzzzz_xzzzz_0, tg_xxzzzz_xzzzz_1, \
                                         tg_xxzzzz_yyyyy_0, tg_xxzzzz_yyyyy_1, tg_xxzzzz_yyyyz_0, tg_xxzzzz_yyyyz_1, \
                                         tg_xxzzzz_yyyzz_0, tg_xxzzzz_yyyzz_1, tg_xxzzzz_yyzzz_0, tg_xxzzzz_yyzzz_1, \
                                         tg_xxzzzz_yzzzz_0, tg_xxzzzz_yzzzz_1, tg_xxzzzz_zzzzz_0, tg_xxzzzz_zzzzz_1, \
                                         tg_xyyyyy_xxxxx_0, tg_xyyyyy_xxxxx_1, tg_xyyyyy_xxxxy_0, tg_xyyyyy_xxxxy_1, \
                                         tg_xyyyyy_xxxxz_0, tg_xyyyyy_xxxxz_1, tg_xyyyyy_xxxyy_0, tg_xyyyyy_xxxyy_1, \
                                         tg_xyyyyy_xxxyz_0, tg_xyyyyy_xxxyz_1, tg_xyyyyy_xxxzz_0, tg_xyyyyy_xxxzz_1, \
                                         tg_xyyyyy_xxyyy_0, tg_xyyyyy_xxyyy_1, tg_xyyyyy_xxyyz_0, tg_xyyyyy_xxyyz_1, \
                                         tg_xyyyyy_xxyzz_0, tg_xyyyyy_xxyzz_1, tg_xyyyyy_xxzzz_0, tg_xyyyyy_xxzzz_1, \
                                         tg_xyyyyy_xyyyy_0, tg_xyyyyy_xyyyy_1, tg_xyyyyy_xyyyz_0, tg_xyyyyy_xyyyz_1, \
                                         tg_xyyyyy_xyyzz_0, tg_xyyyyy_xyyzz_1, tg_xyyyyy_xyzzz_0, tg_xyyyyy_xyzzz_1, \
                                         tg_xyyyyy_xzzzz_0, tg_xyyyyy_xzzzz_1, tg_xyyyyy_yyyyy_0, tg_xyyyyy_yyyyy_1, \
                                         tg_xyyyyy_yyyyz_0, tg_xyyyyy_yyyyz_1, tg_xyyyyy_yyyzz_0, tg_xyyyyy_yyyzz_1, \
                                         tg_xyyyyy_yyzzz_0, tg_xyyyyy_yyzzz_1, tg_xyyyyy_yzzzz_0, tg_xyyyyy_yzzzz_1, \
                                         tg_xyyyyy_zzzzz_0, tg_xyyyyy_zzzzz_1, tg_xyyyyz_xxxxx_0, tg_xyyyyz_xxxxx_1, \
                                         tg_xyyyyz_xxxxy_0, tg_xyyyyz_xxxxy_1, tg_xyyyyz_xxxxz_0, tg_xyyyyz_xxxxz_1, \
                                         tg_xyyyyz_xxxyy_0, tg_xyyyyz_xxxyy_1, tg_xyyyyz_xxxyz_0, tg_xyyyyz_xxxyz_1, \
                                         tg_xyyyyz_xxxzz_0, tg_xyyyyz_xxxzz_1, tg_xyyyyz_xxyyy_0, tg_xyyyyz_xxyyy_1, \
                                         tg_xyyyyz_xxyyz_0, tg_xyyyyz_xxyyz_1, tg_xyyyyz_xxyzz_0, tg_xyyyyz_xxyzz_1, \
                                         tg_xyyyyz_xxzzz_0, tg_xyyyyz_xxzzz_1, tg_xyyyyz_xyyyy_0, tg_xyyyyz_xyyyy_1, \
                                         tg_xyyyyz_xyyyz_0, tg_xyyyyz_xyyyz_1, tg_xyyyyz_xyyzz_0, tg_xyyyyz_xyyzz_1, \
                                         tg_xyyyyz_xyzzz_0, tg_xyyyyz_xyzzz_1, tg_xyyyyz_xzzzz_0, tg_xyyyyz_xzzzz_1, \
                                         tg_xyyyyz_yyyyy_0, tg_xyyyyz_yyyyy_1, tg_xyyyyz_yyyyz_0, tg_xyyyyz_yyyyz_1, \
                                         tg_xyyyyz_yyyzz_0, tg_xyyyyz_yyyzz_1, tg_xyyyyz_yyzzz_0, tg_xyyyyz_yyzzz_1, \
                                         tg_xyyyyz_yzzzz_0, tg_xyyyyz_yzzzz_1, tg_xyyyyz_zzzzz_0, tg_xyyyyz_zzzzz_1, \
                                         tg_xyyyzz_xxxxx_0, tg_xyyyzz_xxxxx_1, tg_xyyyzz_xxxxy_0, tg_xyyyzz_xxxxy_1, \
                                         tg_xyyyzz_xxxxz_0, tg_xyyyzz_xxxxz_1, tg_xyyyzz_xxxyy_0, tg_xyyyzz_xxxyy_1, \
                                         tg_xyyyzz_xxxyz_0, tg_xyyyzz_xxxyz_1, tg_xyyyzz_xxxzz_0, tg_xyyyzz_xxxzz_1, \
                                         tg_xyyyzz_xxyyy_0, tg_xyyyzz_xxyyy_1, tg_xyyyzz_xxyyz_0, tg_xyyyzz_xxyyz_1, \
                                         tg_xyyyzz_xxyzz_0, tg_xyyyzz_xxyzz_1, tg_xyyyzz_xxzzz_0, tg_xyyyzz_xxzzz_1, \
                                         tg_xyyyzz_xyyyy_0, tg_xyyyzz_xyyyy_1, tg_xyyyzz_xyyyz_0, tg_xyyyzz_xyyyz_1, \
                                         tg_xyyyzz_xyyzz_0, tg_xyyyzz_xyyzz_1, tg_xyyyzz_xyzzz_0, tg_xyyyzz_xyzzz_1, \
                                         tg_xyyyzz_xzzzz_0, tg_xyyyzz_xzzzz_1, tg_xyyyzz_yyyyy_0, tg_xyyyzz_yyyyy_1, \
                                         tg_xyyyzz_yyyyz_0, tg_xyyyzz_yyyyz_1, tg_xyyyzz_yyyzz_0, tg_xyyyzz_yyyzz_1, \
                                         tg_xyyyzz_yyzzz_0, tg_xyyyzz_yyzzz_1, tg_xyyyzz_yzzzz_0, tg_xyyyzz_yzzzz_1, \
                                         tg_xyyyzz_zzzzz_0, tg_xyyyzz_zzzzz_1, tg_xyyzzz_xxxxx_0, tg_xyyzzz_xxxxx_1, \
                                         tg_xyyzzz_xxxxy_0, tg_xyyzzz_xxxxy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxyzzz_xyyzz_0[j] = pb_x * tg_xxxyzzz_xyyzz_0[j] + fr * tg_xxxyzzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xyyzz_0[j] - tg_xxyzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_yyzz_1[j];

                    tg_xxxxyzzz_xyzzz_0[j] = pb_x * tg_xxxyzzz_xyzzz_0[j] + fr * tg_xxxyzzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xyzzz_0[j] - tg_xxyzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_yzzz_1[j];

                    tg_xxxxyzzz_xzzzz_0[j] = pb_x * tg_xxxyzzz_xzzzz_0[j] + fr * tg_xxxyzzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xzzzz_0[j] - tg_xxyzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_zzzz_1[j];

                    tg_xxxxyzzz_yyyyy_0[j] = pb_x * tg_xxxyzzz_yyyyy_0[j] + fr * tg_xxxyzzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yyyyy_0[j] - tg_xxyzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxyzzz_yyyyz_0[j] = pb_x * tg_xxxyzzz_yyyyz_0[j] + fr * tg_xxxyzzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yyyyz_0[j] - tg_xxyzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxyzzz_yyyzz_0[j] = pb_x * tg_xxxyzzz_yyyzz_0[j] + fr * tg_xxxyzzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yyyzz_0[j] - tg_xxyzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxyzzz_yyzzz_0[j] = pb_x * tg_xxxyzzz_yyzzz_0[j] + fr * tg_xxxyzzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yyzzz_0[j] - tg_xxyzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxyzzz_yzzzz_0[j] = pb_x * tg_xxxyzzz_yzzzz_0[j] + fr * tg_xxxyzzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yzzzz_0[j] - tg_xxyzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxyzzz_zzzzz_0[j] = pb_x * tg_xxxyzzz_zzzzz_0[j] + fr * tg_xxxyzzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_zzzzz_0[j] - tg_xxyzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxxzzzz_xxxxx_0[j] = pb_x * tg_xxxzzzz_xxxxx_0[j] + fr * tg_xxxzzzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxxxx_0[j] - tg_xxzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxzzzz_xxxx_1[j];

                    tg_xxxxzzzz_xxxxy_0[j] = pb_x * tg_xxxzzzz_xxxxy_0[j] + fr * tg_xxxzzzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxxxy_0[j] - tg_xxzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzzz_xxxy_1[j];

                    tg_xxxxzzzz_xxxxz_0[j] = pb_x * tg_xxxzzzz_xxxxz_0[j] + fr * tg_xxxzzzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxxxz_0[j] - tg_xxzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzzz_xxxz_1[j];

                    tg_xxxxzzzz_xxxyy_0[j] = pb_x * tg_xxxzzzz_xxxyy_0[j] + fr * tg_xxxzzzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxxyy_0[j] - tg_xxzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzzz_xxyy_1[j];

                    tg_xxxxzzzz_xxxyz_0[j] = pb_x * tg_xxxzzzz_xxxyz_0[j] + fr * tg_xxxzzzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxxyz_0[j] - tg_xxzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzzz_xxyz_1[j];

                    tg_xxxxzzzz_xxxzz_0[j] = pb_x * tg_xxxzzzz_xxxzz_0[j] + fr * tg_xxxzzzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxxzz_0[j] - tg_xxzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzzz_xxzz_1[j];

                    tg_xxxxzzzz_xxyyy_0[j] = pb_x * tg_xxxzzzz_xxyyy_0[j] + fr * tg_xxxzzzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxyyy_0[j] - tg_xxzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzzz_xyyy_1[j];

                    tg_xxxxzzzz_xxyyz_0[j] = pb_x * tg_xxxzzzz_xxyyz_0[j] + fr * tg_xxxzzzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxyyz_0[j] - tg_xxzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzzz_xyyz_1[j];

                    tg_xxxxzzzz_xxyzz_0[j] = pb_x * tg_xxxzzzz_xxyzz_0[j] + fr * tg_xxxzzzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxyzz_0[j] - tg_xxzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzzz_xyzz_1[j];

                    tg_xxxxzzzz_xxzzz_0[j] = pb_x * tg_xxxzzzz_xxzzz_0[j] + fr * tg_xxxzzzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxzzz_0[j] - tg_xxzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzzz_xzzz_1[j];

                    tg_xxxxzzzz_xyyyy_0[j] = pb_x * tg_xxxzzzz_xyyyy_0[j] + fr * tg_xxxzzzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xyyyy_0[j] - tg_xxzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_yyyy_1[j];

                    tg_xxxxzzzz_xyyyz_0[j] = pb_x * tg_xxxzzzz_xyyyz_0[j] + fr * tg_xxxzzzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xyyyz_0[j] - tg_xxzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_yyyz_1[j];

                    tg_xxxxzzzz_xyyzz_0[j] = pb_x * tg_xxxzzzz_xyyzz_0[j] + fr * tg_xxxzzzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xyyzz_0[j] - tg_xxzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_yyzz_1[j];

                    tg_xxxxzzzz_xyzzz_0[j] = pb_x * tg_xxxzzzz_xyzzz_0[j] + fr * tg_xxxzzzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xyzzz_0[j] - tg_xxzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_yzzz_1[j];

                    tg_xxxxzzzz_xzzzz_0[j] = pb_x * tg_xxxzzzz_xzzzz_0[j] + fr * tg_xxxzzzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xzzzz_0[j] - tg_xxzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_zzzz_1[j];

                    tg_xxxxzzzz_yyyyy_0[j] = pb_x * tg_xxxzzzz_yyyyy_0[j] + fr * tg_xxxzzzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yyyyy_0[j] - tg_xxzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxzzzz_yyyyz_0[j] = pb_x * tg_xxxzzzz_yyyyz_0[j] + fr * tg_xxxzzzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yyyyz_0[j] - tg_xxzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxzzzz_yyyzz_0[j] = pb_x * tg_xxxzzzz_yyyzz_0[j] + fr * tg_xxxzzzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yyyzz_0[j] - tg_xxzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxzzzz_yyzzz_0[j] = pb_x * tg_xxxzzzz_yyzzz_0[j] + fr * tg_xxxzzzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yyzzz_0[j] - tg_xxzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxzzzz_yzzzz_0[j] = pb_x * tg_xxxzzzz_yzzzz_0[j] + fr * tg_xxxzzzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yzzzz_0[j] - tg_xxzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxzzzz_zzzzz_0[j] = pb_x * tg_xxxzzzz_zzzzz_0[j] + fr * tg_xxxzzzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_zzzzz_0[j] - tg_xxzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyyyy_xxxxx_0[j] = pb_x * tg_xxyyyyy_xxxxx_0[j] + fr * tg_xxyyyyy_xxxxx_1[j] + fl1_fx * (tg_xyyyyy_xxxxx_0[j] - tg_xyyyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyyy_xxxx_1[j];

                    tg_xxxyyyyy_xxxxy_0[j] = pb_x * tg_xxyyyyy_xxxxy_0[j] + fr * tg_xxyyyyy_xxxxy_1[j] + fl1_fx * (tg_xyyyyy_xxxxy_0[j] - tg_xyyyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyyy_xxxy_1[j];

                    tg_xxxyyyyy_xxxxz_0[j] = pb_x * tg_xxyyyyy_xxxxz_0[j] + fr * tg_xxyyyyy_xxxxz_1[j] + fl1_fx * (tg_xyyyyy_xxxxz_0[j] - tg_xyyyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyyy_xxxz_1[j];

                    tg_xxxyyyyy_xxxyy_0[j] = pb_x * tg_xxyyyyy_xxxyy_0[j] + fr * tg_xxyyyyy_xxxyy_1[j] + fl1_fx * (tg_xyyyyy_xxxyy_0[j] - tg_xyyyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyy_xxyy_1[j];

                    tg_xxxyyyyy_xxxyz_0[j] = pb_x * tg_xxyyyyy_xxxyz_0[j] + fr * tg_xxyyyyy_xxxyz_1[j] + fl1_fx * (tg_xyyyyy_xxxyz_0[j] - tg_xyyyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyy_xxyz_1[j];

                    tg_xxxyyyyy_xxxzz_0[j] = pb_x * tg_xxyyyyy_xxxzz_0[j] + fr * tg_xxyyyyy_xxxzz_1[j] + fl1_fx * (tg_xyyyyy_xxxzz_0[j] - tg_xyyyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyy_xxzz_1[j];

                    tg_xxxyyyyy_xxyyy_0[j] = pb_x * tg_xxyyyyy_xxyyy_0[j] + fr * tg_xxyyyyy_xxyyy_1[j] + fl1_fx * (tg_xyyyyy_xxyyy_0[j] - tg_xyyyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyy_xyyy_1[j];

                    tg_xxxyyyyy_xxyyz_0[j] = pb_x * tg_xxyyyyy_xxyyz_0[j] + fr * tg_xxyyyyy_xxyyz_1[j] + fl1_fx * (tg_xyyyyy_xxyyz_0[j] - tg_xyyyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyy_xyyz_1[j];

                    tg_xxxyyyyy_xxyzz_0[j] = pb_x * tg_xxyyyyy_xxyzz_0[j] + fr * tg_xxyyyyy_xxyzz_1[j] + fl1_fx * (tg_xyyyyy_xxyzz_0[j] - tg_xyyyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyy_xyzz_1[j];

                    tg_xxxyyyyy_xxzzz_0[j] = pb_x * tg_xxyyyyy_xxzzz_0[j] + fr * tg_xxyyyyy_xxzzz_1[j] + fl1_fx * (tg_xyyyyy_xxzzz_0[j] - tg_xyyyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyy_xzzz_1[j];

                    tg_xxxyyyyy_xyyyy_0[j] = pb_x * tg_xxyyyyy_xyyyy_0[j] + fr * tg_xxyyyyy_xyyyy_1[j] + fl1_fx * (tg_xyyyyy_xyyyy_0[j] - tg_xyyyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_yyyy_1[j];

                    tg_xxxyyyyy_xyyyz_0[j] = pb_x * tg_xxyyyyy_xyyyz_0[j] + fr * tg_xxyyyyy_xyyyz_1[j] + fl1_fx * (tg_xyyyyy_xyyyz_0[j] - tg_xyyyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_yyyz_1[j];

                    tg_xxxyyyyy_xyyzz_0[j] = pb_x * tg_xxyyyyy_xyyzz_0[j] + fr * tg_xxyyyyy_xyyzz_1[j] + fl1_fx * (tg_xyyyyy_xyyzz_0[j] - tg_xyyyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_yyzz_1[j];

                    tg_xxxyyyyy_xyzzz_0[j] = pb_x * tg_xxyyyyy_xyzzz_0[j] + fr * tg_xxyyyyy_xyzzz_1[j] + fl1_fx * (tg_xyyyyy_xyzzz_0[j] - tg_xyyyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_yzzz_1[j];

                    tg_xxxyyyyy_xzzzz_0[j] = pb_x * tg_xxyyyyy_xzzzz_0[j] + fr * tg_xxyyyyy_xzzzz_1[j] + fl1_fx * (tg_xyyyyy_xzzzz_0[j] - tg_xyyyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_zzzz_1[j];

                    tg_xxxyyyyy_yyyyy_0[j] = pb_x * tg_xxyyyyy_yyyyy_0[j] + fr * tg_xxyyyyy_yyyyy_1[j] + fl1_fx * (tg_xyyyyy_yyyyy_0[j] - tg_xyyyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyyyy_yyyyz_0[j] = pb_x * tg_xxyyyyy_yyyyz_0[j] + fr * tg_xxyyyyy_yyyyz_1[j] + fl1_fx * (tg_xyyyyy_yyyyz_0[j] - tg_xyyyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyyyy_yyyzz_0[j] = pb_x * tg_xxyyyyy_yyyzz_0[j] + fr * tg_xxyyyyy_yyyzz_1[j] + fl1_fx * (tg_xyyyyy_yyyzz_0[j] - tg_xyyyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyyyy_yyzzz_0[j] = pb_x * tg_xxyyyyy_yyzzz_0[j] + fr * tg_xxyyyyy_yyzzz_1[j] + fl1_fx * (tg_xyyyyy_yyzzz_0[j] - tg_xyyyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyyyy_yzzzz_0[j] = pb_x * tg_xxyyyyy_yzzzz_0[j] + fr * tg_xxyyyyy_yzzzz_1[j] + fl1_fx * (tg_xyyyyy_yzzzz_0[j] - tg_xyyyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyyyy_zzzzz_0[j] = pb_x * tg_xxyyyyy_zzzzz_0[j] + fr * tg_xxyyyyy_zzzzz_1[j] + fl1_fx * (tg_xyyyyy_zzzzz_0[j] - tg_xyyyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyyyz_xxxxx_0[j] = pb_x * tg_xxyyyyz_xxxxx_0[j] + fr * tg_xxyyyyz_xxxxx_1[j] + fl1_fx * (tg_xyyyyz_xxxxx_0[j] - tg_xyyyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyyz_xxxx_1[j];

                    tg_xxxyyyyz_xxxxy_0[j] = pb_x * tg_xxyyyyz_xxxxy_0[j] + fr * tg_xxyyyyz_xxxxy_1[j] + fl1_fx * (tg_xyyyyz_xxxxy_0[j] - tg_xyyyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyyz_xxxy_1[j];

                    tg_xxxyyyyz_xxxxz_0[j] = pb_x * tg_xxyyyyz_xxxxz_0[j] + fr * tg_xxyyyyz_xxxxz_1[j] + fl1_fx * (tg_xyyyyz_xxxxz_0[j] - tg_xyyyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyyz_xxxz_1[j];

                    tg_xxxyyyyz_xxxyy_0[j] = pb_x * tg_xxyyyyz_xxxyy_0[j] + fr * tg_xxyyyyz_xxxyy_1[j] + fl1_fx * (tg_xyyyyz_xxxyy_0[j] - tg_xyyyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyz_xxyy_1[j];

                    tg_xxxyyyyz_xxxyz_0[j] = pb_x * tg_xxyyyyz_xxxyz_0[j] + fr * tg_xxyyyyz_xxxyz_1[j] + fl1_fx * (tg_xyyyyz_xxxyz_0[j] - tg_xyyyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyz_xxyz_1[j];

                    tg_xxxyyyyz_xxxzz_0[j] = pb_x * tg_xxyyyyz_xxxzz_0[j] + fr * tg_xxyyyyz_xxxzz_1[j] + fl1_fx * (tg_xyyyyz_xxxzz_0[j] - tg_xyyyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyz_xxzz_1[j];

                    tg_xxxyyyyz_xxyyy_0[j] = pb_x * tg_xxyyyyz_xxyyy_0[j] + fr * tg_xxyyyyz_xxyyy_1[j] + fl1_fx * (tg_xyyyyz_xxyyy_0[j] - tg_xyyyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyz_xyyy_1[j];

                    tg_xxxyyyyz_xxyyz_0[j] = pb_x * tg_xxyyyyz_xxyyz_0[j] + fr * tg_xxyyyyz_xxyyz_1[j] + fl1_fx * (tg_xyyyyz_xxyyz_0[j] - tg_xyyyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyz_xyyz_1[j];

                    tg_xxxyyyyz_xxyzz_0[j] = pb_x * tg_xxyyyyz_xxyzz_0[j] + fr * tg_xxyyyyz_xxyzz_1[j] + fl1_fx * (tg_xyyyyz_xxyzz_0[j] - tg_xyyyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyz_xyzz_1[j];

                    tg_xxxyyyyz_xxzzz_0[j] = pb_x * tg_xxyyyyz_xxzzz_0[j] + fr * tg_xxyyyyz_xxzzz_1[j] + fl1_fx * (tg_xyyyyz_xxzzz_0[j] - tg_xyyyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyz_xzzz_1[j];

                    tg_xxxyyyyz_xyyyy_0[j] = pb_x * tg_xxyyyyz_xyyyy_0[j] + fr * tg_xxyyyyz_xyyyy_1[j] + fl1_fx * (tg_xyyyyz_xyyyy_0[j] - tg_xyyyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_yyyy_1[j];

                    tg_xxxyyyyz_xyyyz_0[j] = pb_x * tg_xxyyyyz_xyyyz_0[j] + fr * tg_xxyyyyz_xyyyz_1[j] + fl1_fx * (tg_xyyyyz_xyyyz_0[j] - tg_xyyyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_yyyz_1[j];

                    tg_xxxyyyyz_xyyzz_0[j] = pb_x * tg_xxyyyyz_xyyzz_0[j] + fr * tg_xxyyyyz_xyyzz_1[j] + fl1_fx * (tg_xyyyyz_xyyzz_0[j] - tg_xyyyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_yyzz_1[j];

                    tg_xxxyyyyz_xyzzz_0[j] = pb_x * tg_xxyyyyz_xyzzz_0[j] + fr * tg_xxyyyyz_xyzzz_1[j] + fl1_fx * (tg_xyyyyz_xyzzz_0[j] - tg_xyyyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_yzzz_1[j];

                    tg_xxxyyyyz_xzzzz_0[j] = pb_x * tg_xxyyyyz_xzzzz_0[j] + fr * tg_xxyyyyz_xzzzz_1[j] + fl1_fx * (tg_xyyyyz_xzzzz_0[j] - tg_xyyyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_zzzz_1[j];

                    tg_xxxyyyyz_yyyyy_0[j] = pb_x * tg_xxyyyyz_yyyyy_0[j] + fr * tg_xxyyyyz_yyyyy_1[j] + fl1_fx * (tg_xyyyyz_yyyyy_0[j] - tg_xyyyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyyyz_yyyyz_0[j] = pb_x * tg_xxyyyyz_yyyyz_0[j] + fr * tg_xxyyyyz_yyyyz_1[j] + fl1_fx * (tg_xyyyyz_yyyyz_0[j] - tg_xyyyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyyyz_yyyzz_0[j] = pb_x * tg_xxyyyyz_yyyzz_0[j] + fr * tg_xxyyyyz_yyyzz_1[j] + fl1_fx * (tg_xyyyyz_yyyzz_0[j] - tg_xyyyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyyyz_yyzzz_0[j] = pb_x * tg_xxyyyyz_yyzzz_0[j] + fr * tg_xxyyyyz_yyzzz_1[j] + fl1_fx * (tg_xyyyyz_yyzzz_0[j] - tg_xyyyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyyyz_yzzzz_0[j] = pb_x * tg_xxyyyyz_yzzzz_0[j] + fr * tg_xxyyyyz_yzzzz_1[j] + fl1_fx * (tg_xyyyyz_yzzzz_0[j] - tg_xyyyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyyyz_zzzzz_0[j] = pb_x * tg_xxyyyyz_zzzzz_0[j] + fr * tg_xxyyyyz_zzzzz_1[j] + fl1_fx * (tg_xyyyyz_zzzzz_0[j] - tg_xyyyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyyzz_xxxxx_0[j] = pb_x * tg_xxyyyzz_xxxxx_0[j] + fr * tg_xxyyyzz_xxxxx_1[j] + fl1_fx * (tg_xyyyzz_xxxxx_0[j] - tg_xyyyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyzz_xxxx_1[j];

                    tg_xxxyyyzz_xxxxy_0[j] = pb_x * tg_xxyyyzz_xxxxy_0[j] + fr * tg_xxyyyzz_xxxxy_1[j] + fl1_fx * (tg_xyyyzz_xxxxy_0[j] - tg_xyyyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyzz_xxxy_1[j];

                    tg_xxxyyyzz_xxxxz_0[j] = pb_x * tg_xxyyyzz_xxxxz_0[j] + fr * tg_xxyyyzz_xxxxz_1[j] + fl1_fx * (tg_xyyyzz_xxxxz_0[j] - tg_xyyyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyzz_xxxz_1[j];

                    tg_xxxyyyzz_xxxyy_0[j] = pb_x * tg_xxyyyzz_xxxyy_0[j] + fr * tg_xxyyyzz_xxxyy_1[j] + fl1_fx * (tg_xyyyzz_xxxyy_0[j] - tg_xyyyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyzz_xxyy_1[j];

                    tg_xxxyyyzz_xxxyz_0[j] = pb_x * tg_xxyyyzz_xxxyz_0[j] + fr * tg_xxyyyzz_xxxyz_1[j] + fl1_fx * (tg_xyyyzz_xxxyz_0[j] - tg_xyyyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyzz_xxyz_1[j];

                    tg_xxxyyyzz_xxxzz_0[j] = pb_x * tg_xxyyyzz_xxxzz_0[j] + fr * tg_xxyyyzz_xxxzz_1[j] + fl1_fx * (tg_xyyyzz_xxxzz_0[j] - tg_xyyyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyzz_xxzz_1[j];

                    tg_xxxyyyzz_xxyyy_0[j] = pb_x * tg_xxyyyzz_xxyyy_0[j] + fr * tg_xxyyyzz_xxyyy_1[j] + fl1_fx * (tg_xyyyzz_xxyyy_0[j] - tg_xyyyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyzz_xyyy_1[j];

                    tg_xxxyyyzz_xxyyz_0[j] = pb_x * tg_xxyyyzz_xxyyz_0[j] + fr * tg_xxyyyzz_xxyyz_1[j] + fl1_fx * (tg_xyyyzz_xxyyz_0[j] - tg_xyyyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyzz_xyyz_1[j];

                    tg_xxxyyyzz_xxyzz_0[j] = pb_x * tg_xxyyyzz_xxyzz_0[j] + fr * tg_xxyyyzz_xxyzz_1[j] + fl1_fx * (tg_xyyyzz_xxyzz_0[j] - tg_xyyyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyzz_xyzz_1[j];

                    tg_xxxyyyzz_xxzzz_0[j] = pb_x * tg_xxyyyzz_xxzzz_0[j] + fr * tg_xxyyyzz_xxzzz_1[j] + fl1_fx * (tg_xyyyzz_xxzzz_0[j] - tg_xyyyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyzz_xzzz_1[j];

                    tg_xxxyyyzz_xyyyy_0[j] = pb_x * tg_xxyyyzz_xyyyy_0[j] + fr * tg_xxyyyzz_xyyyy_1[j] + fl1_fx * (tg_xyyyzz_xyyyy_0[j] - tg_xyyyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_yyyy_1[j];

                    tg_xxxyyyzz_xyyyz_0[j] = pb_x * tg_xxyyyzz_xyyyz_0[j] + fr * tg_xxyyyzz_xyyyz_1[j] + fl1_fx * (tg_xyyyzz_xyyyz_0[j] - tg_xyyyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_yyyz_1[j];

                    tg_xxxyyyzz_xyyzz_0[j] = pb_x * tg_xxyyyzz_xyyzz_0[j] + fr * tg_xxyyyzz_xyyzz_1[j] + fl1_fx * (tg_xyyyzz_xyyzz_0[j] - tg_xyyyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_yyzz_1[j];

                    tg_xxxyyyzz_xyzzz_0[j] = pb_x * tg_xxyyyzz_xyzzz_0[j] + fr * tg_xxyyyzz_xyzzz_1[j] + fl1_fx * (tg_xyyyzz_xyzzz_0[j] - tg_xyyyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_yzzz_1[j];

                    tg_xxxyyyzz_xzzzz_0[j] = pb_x * tg_xxyyyzz_xzzzz_0[j] + fr * tg_xxyyyzz_xzzzz_1[j] + fl1_fx * (tg_xyyyzz_xzzzz_0[j] - tg_xyyyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_zzzz_1[j];

                    tg_xxxyyyzz_yyyyy_0[j] = pb_x * tg_xxyyyzz_yyyyy_0[j] + fr * tg_xxyyyzz_yyyyy_1[j] + fl1_fx * (tg_xyyyzz_yyyyy_0[j] - tg_xyyyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyyzz_yyyyz_0[j] = pb_x * tg_xxyyyzz_yyyyz_0[j] + fr * tg_xxyyyzz_yyyyz_1[j] + fl1_fx * (tg_xyyyzz_yyyyz_0[j] - tg_xyyyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyyzz_yyyzz_0[j] = pb_x * tg_xxyyyzz_yyyzz_0[j] + fr * tg_xxyyyzz_yyyzz_1[j] + fl1_fx * (tg_xyyyzz_yyyzz_0[j] - tg_xyyyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyyzz_yyzzz_0[j] = pb_x * tg_xxyyyzz_yyzzz_0[j] + fr * tg_xxyyyzz_yyzzz_1[j] + fl1_fx * (tg_xyyyzz_yyzzz_0[j] - tg_xyyyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyyzz_yzzzz_0[j] = pb_x * tg_xxyyyzz_yzzzz_0[j] + fr * tg_xxyyyzz_yzzzz_1[j] + fl1_fx * (tg_xyyyzz_yzzzz_0[j] - tg_xyyyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyyzz_zzzzz_0[j] = pb_x * tg_xxyyyzz_zzzzz_0[j] + fr * tg_xxyyyzz_zzzzz_1[j] + fl1_fx * (tg_xyyyzz_zzzzz_0[j] - tg_xyyyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyyzzz_xxxxx_0[j] = pb_x * tg_xxyyzzz_xxxxx_0[j] + fr * tg_xxyyzzz_xxxxx_1[j] + fl1_fx * (tg_xyyzzz_xxxxx_0[j] - tg_xyyzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyzzz_xxxx_1[j];

                    tg_xxxyyzzz_xxxxy_0[j] = pb_x * tg_xxyyzzz_xxxxy_0[j] + fr * tg_xxyyzzz_xxxxy_1[j] + fl1_fx * (tg_xyyzzz_xxxxy_0[j] - tg_xyyzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzzz_xxxy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_380_475(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xyyyyyz_xyyzz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 474); 

                auto tg_xxyyzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 380); 

                auto tg_xxyyzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 381); 

                auto tg_xxyyzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 382); 

                auto tg_xxyyzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 383); 

                auto tg_xxyyzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 384); 

                auto tg_xxyyzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 385); 

                auto tg_xxyyzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 386); 

                auto tg_xxyyzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 387); 

                auto tg_xxyyzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 388); 

                auto tg_xxyyzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 389); 

                auto tg_xxyyzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 390); 

                auto tg_xxyyzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 391); 

                auto tg_xxyyzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 392); 

                auto tg_xxyyzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 393); 

                auto tg_xxyyzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 394); 

                auto tg_xxyyzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 395); 

                auto tg_xxyyzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 396); 

                auto tg_xxyyzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 397); 

                auto tg_xxyyzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 398); 

                auto tg_xxyzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 399); 

                auto tg_xxyzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 400); 

                auto tg_xxyzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 401); 

                auto tg_xxyzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 402); 

                auto tg_xxyzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 403); 

                auto tg_xxyzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 404); 

                auto tg_xxyzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 405); 

                auto tg_xxyzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 406); 

                auto tg_xxyzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 407); 

                auto tg_xxyzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 408); 

                auto tg_xxyzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 409); 

                auto tg_xxyzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 410); 

                auto tg_xxyzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 411); 

                auto tg_xxyzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 412); 

                auto tg_xxyzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 413); 

                auto tg_xxyzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 414); 

                auto tg_xxyzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 415); 

                auto tg_xxyzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 416); 

                auto tg_xxyzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 417); 

                auto tg_xxyzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 418); 

                auto tg_xxyzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 419); 

                auto tg_xxzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 420); 

                auto tg_xxzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 421); 

                auto tg_xxzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 422); 

                auto tg_xxzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 423); 

                auto tg_xxzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 424); 

                auto tg_xxzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 425); 

                auto tg_xxzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 426); 

                auto tg_xxzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 427); 

                auto tg_xxzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 428); 

                auto tg_xxzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 429); 

                auto tg_xxzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 430); 

                auto tg_xxzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 431); 

                auto tg_xxzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 432); 

                auto tg_xxzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 433); 

                auto tg_xxzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 434); 

                auto tg_xxzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 435); 

                auto tg_xxzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 436); 

                auto tg_xxzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 437); 

                auto tg_xxzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 438); 

                auto tg_xxzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 439); 

                auto tg_xxzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 440); 

                auto tg_xyyyyyy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 441); 

                auto tg_xyyyyyy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 442); 

                auto tg_xyyyyyy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 443); 

                auto tg_xyyyyyy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 444); 

                auto tg_xyyyyyy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 445); 

                auto tg_xyyyyyy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 446); 

                auto tg_xyyyyyy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 447); 

                auto tg_xyyyyyy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 448); 

                auto tg_xyyyyyy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 449); 

                auto tg_xyyyyyy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 450); 

                auto tg_xyyyyyy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 451); 

                auto tg_xyyyyyy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 452); 

                auto tg_xyyyyyy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 453); 

                auto tg_xyyyyyy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 454); 

                auto tg_xyyyyyy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 455); 

                auto tg_xyyyyyy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 456); 

                auto tg_xyyyyyy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 457); 

                auto tg_xyyyyyy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 458); 

                auto tg_xyyyyyy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 459); 

                auto tg_xyyyyyy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 460); 

                auto tg_xyyyyyy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 461); 

                auto tg_xyyyyyz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 462); 

                auto tg_xyyyyyz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 463); 

                auto tg_xyyyyyz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 464); 

                auto tg_xyyyyyz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 465); 

                auto tg_xyyyyyz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 466); 

                auto tg_xyyyyyz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 467); 

                auto tg_xyyyyyz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 468); 

                auto tg_xyyyyyz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 469); 

                auto tg_xyyyyyz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 470); 

                auto tg_xyyyyyz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 471); 

                auto tg_xyyyyyz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 472); 

                auto tg_xyyyyyz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 473); 

                auto tg_xyyyyyz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 474); 

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

                auto tg_yyyyyz_xyyzz_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 474); 

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

                auto tg_yyyyyz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 474); 

                auto tg_xxyyzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 272); 

                auto tg_xxyyzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 273); 

                auto tg_xxyyzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 274); 

                auto tg_xxyyzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 275); 

                auto tg_xxyyzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 276); 

                auto tg_xxyyzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 277); 

                auto tg_xxyyzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 278); 

                auto tg_xxyyzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 279); 

                auto tg_xxyyzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 280); 

                auto tg_xxyyzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 281); 

                auto tg_xxyyzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 282); 

                auto tg_xxyyzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 283); 

                auto tg_xxyyzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 284); 

                auto tg_xxyzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 285); 

                auto tg_xxyzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 286); 

                auto tg_xxyzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 287); 

                auto tg_xxyzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 288); 

                auto tg_xxyzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 289); 

                auto tg_xxyzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 290); 

                auto tg_xxyzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 291); 

                auto tg_xxyzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 292); 

                auto tg_xxyzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 293); 

                auto tg_xxyzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 294); 

                auto tg_xxyzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 295); 

                auto tg_xxyzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 296); 

                auto tg_xxyzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 297); 

                auto tg_xxyzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 298); 

                auto tg_xxyzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 299); 

                auto tg_xxzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 300); 

                auto tg_xxzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 301); 

                auto tg_xxzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 302); 

                auto tg_xxzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 303); 

                auto tg_xxzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 304); 

                auto tg_xxzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 305); 

                auto tg_xxzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 306); 

                auto tg_xxzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 307); 

                auto tg_xxzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 308); 

                auto tg_xxzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 309); 

                auto tg_xxzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 310); 

                auto tg_xxzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 311); 

                auto tg_xxzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 312); 

                auto tg_xxzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 313); 

                auto tg_xxzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 314); 

                auto tg_xyyyyyy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 315); 

                auto tg_xyyyyyy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 316); 

                auto tg_xyyyyyy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 317); 

                auto tg_xyyyyyy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 318); 

                auto tg_xyyyyyy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 319); 

                auto tg_xyyyyyy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 320); 

                auto tg_xyyyyyy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 321); 

                auto tg_xyyyyyy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 322); 

                auto tg_xyyyyyy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 323); 

                auto tg_xyyyyyy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 324); 

                auto tg_xyyyyyy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 325); 

                auto tg_xyyyyyy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 326); 

                auto tg_xyyyyyy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 327); 

                auto tg_xyyyyyy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 328); 

                auto tg_xyyyyyy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 329); 

                auto tg_xyyyyyz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 330); 

                auto tg_xyyyyyz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 331); 

                auto tg_xyyyyyz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 332); 

                auto tg_xyyyyyz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 333); 

                auto tg_xyyyyyz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 334); 

                auto tg_xyyyyyz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 335); 

                auto tg_xyyyyyz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 336); 

                auto tg_xyyyyyz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 337); 

                auto tg_xyyyyyz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 338); 

                auto tg_xyyyyyz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 339); 

                auto tg_xyyyyyz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 340); 

                auto tg_xyyyyyz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 341); 

                auto tg_xyyyyyz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 342); 

                // set up pointers to integrals

                auto tg_xxxyyzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 380); 

                auto tg_xxxyyzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 381); 

                auto tg_xxxyyzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 382); 

                auto tg_xxxyyzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 383); 

                auto tg_xxxyyzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 384); 

                auto tg_xxxyyzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 385); 

                auto tg_xxxyyzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 386); 

                auto tg_xxxyyzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 387); 

                auto tg_xxxyyzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 388); 

                auto tg_xxxyyzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 389); 

                auto tg_xxxyyzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 390); 

                auto tg_xxxyyzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 391); 

                auto tg_xxxyyzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 392); 

                auto tg_xxxyyzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 393); 

                auto tg_xxxyyzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 394); 

                auto tg_xxxyyzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 395); 

                auto tg_xxxyyzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 396); 

                auto tg_xxxyyzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 397); 

                auto tg_xxxyyzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 398); 

                auto tg_xxxyzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 399); 

                auto tg_xxxyzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 400); 

                auto tg_xxxyzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 401); 

                auto tg_xxxyzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 402); 

                auto tg_xxxyzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 403); 

                auto tg_xxxyzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 404); 

                auto tg_xxxyzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 405); 

                auto tg_xxxyzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 406); 

                auto tg_xxxyzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 407); 

                auto tg_xxxyzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 408); 

                auto tg_xxxyzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 409); 

                auto tg_xxxyzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 410); 

                auto tg_xxxyzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 411); 

                auto tg_xxxyzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 412); 

                auto tg_xxxyzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 413); 

                auto tg_xxxyzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 414); 

                auto tg_xxxyzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 415); 

                auto tg_xxxyzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 416); 

                auto tg_xxxyzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 417); 

                auto tg_xxxyzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 418); 

                auto tg_xxxyzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 419); 

                auto tg_xxxzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 420); 

                auto tg_xxxzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 421); 

                auto tg_xxxzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 422); 

                auto tg_xxxzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 423); 

                auto tg_xxxzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 424); 

                auto tg_xxxzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 425); 

                auto tg_xxxzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 426); 

                auto tg_xxxzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 427); 

                auto tg_xxxzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 428); 

                auto tg_xxxzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 429); 

                auto tg_xxxzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 430); 

                auto tg_xxxzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 431); 

                auto tg_xxxzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 432); 

                auto tg_xxxzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 433); 

                auto tg_xxxzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 434); 

                auto tg_xxxzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 435); 

                auto tg_xxxzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 436); 

                auto tg_xxxzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 437); 

                auto tg_xxxzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 438); 

                auto tg_xxxzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 439); 

                auto tg_xxxzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 440); 

                auto tg_xxyyyyyy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 441); 

                auto tg_xxyyyyyy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 442); 

                auto tg_xxyyyyyy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 443); 

                auto tg_xxyyyyyy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 444); 

                auto tg_xxyyyyyy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 445); 

                auto tg_xxyyyyyy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 446); 

                auto tg_xxyyyyyy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 447); 

                auto tg_xxyyyyyy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 448); 

                auto tg_xxyyyyyy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 449); 

                auto tg_xxyyyyyy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 450); 

                auto tg_xxyyyyyy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 451); 

                auto tg_xxyyyyyy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 452); 

                auto tg_xxyyyyyy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 453); 

                auto tg_xxyyyyyy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 454); 

                auto tg_xxyyyyyy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 455); 

                auto tg_xxyyyyyy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 456); 

                auto tg_xxyyyyyy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 457); 

                auto tg_xxyyyyyy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 458); 

                auto tg_xxyyyyyy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 459); 

                auto tg_xxyyyyyy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 460); 

                auto tg_xxyyyyyy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 461); 

                auto tg_xxyyyyyz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 462); 

                auto tg_xxyyyyyz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 463); 

                auto tg_xxyyyyyz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 464); 

                auto tg_xxyyyyyz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 465); 

                auto tg_xxyyyyyz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 466); 

                auto tg_xxyyyyyz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 467); 

                auto tg_xxyyyyyz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 468); 

                auto tg_xxyyyyyz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 469); 

                auto tg_xxyyyyyz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 470); 

                auto tg_xxyyyyyz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 471); 

                auto tg_xxyyyyyz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 472); 

                auto tg_xxyyyyyz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 473); 

                auto tg_xxyyyyyz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 474); 

                // Batch of Integrals (380,475)

                #pragma omp simd aligned(fxn, fza, tg_xxxyyzzz_xxxxz_0, tg_xxxyyzzz_xxxyy_0, \
                                         tg_xxxyyzzz_xxxyz_0, tg_xxxyyzzz_xxxzz_0, tg_xxxyyzzz_xxyyy_0, tg_xxxyyzzz_xxyyz_0, \
                                         tg_xxxyyzzz_xxyzz_0, tg_xxxyyzzz_xxzzz_0, tg_xxxyyzzz_xyyyy_0, tg_xxxyyzzz_xyyyz_0, \
                                         tg_xxxyyzzz_xyyzz_0, tg_xxxyyzzz_xyzzz_0, tg_xxxyyzzz_xzzzz_0, tg_xxxyyzzz_yyyyy_0, \
                                         tg_xxxyyzzz_yyyyz_0, tg_xxxyyzzz_yyyzz_0, tg_xxxyyzzz_yyzzz_0, tg_xxxyyzzz_yzzzz_0, \
                                         tg_xxxyyzzz_zzzzz_0, tg_xxxyzzzz_xxxxx_0, tg_xxxyzzzz_xxxxy_0, tg_xxxyzzzz_xxxxz_0, \
                                         tg_xxxyzzzz_xxxyy_0, tg_xxxyzzzz_xxxyz_0, tg_xxxyzzzz_xxxzz_0, tg_xxxyzzzz_xxyyy_0, \
                                         tg_xxxyzzzz_xxyyz_0, tg_xxxyzzzz_xxyzz_0, tg_xxxyzzzz_xxzzz_0, tg_xxxyzzzz_xyyyy_0, \
                                         tg_xxxyzzzz_xyyyz_0, tg_xxxyzzzz_xyyzz_0, tg_xxxyzzzz_xyzzz_0, tg_xxxyzzzz_xzzzz_0, \
                                         tg_xxxyzzzz_yyyyy_0, tg_xxxyzzzz_yyyyz_0, tg_xxxyzzzz_yyyzz_0, tg_xxxyzzzz_yyzzz_0, \
                                         tg_xxxyzzzz_yzzzz_0, tg_xxxyzzzz_zzzzz_0, tg_xxxzzzzz_xxxxx_0, tg_xxxzzzzz_xxxxy_0, \
                                         tg_xxxzzzzz_xxxxz_0, tg_xxxzzzzz_xxxyy_0, tg_xxxzzzzz_xxxyz_0, tg_xxxzzzzz_xxxzz_0, \
                                         tg_xxxzzzzz_xxyyy_0, tg_xxxzzzzz_xxyyz_0, tg_xxxzzzzz_xxyzz_0, tg_xxxzzzzz_xxzzz_0, \
                                         tg_xxxzzzzz_xyyyy_0, tg_xxxzzzzz_xyyyz_0, tg_xxxzzzzz_xyyzz_0, tg_xxxzzzzz_xyzzz_0, \
                                         tg_xxxzzzzz_xzzzz_0, tg_xxxzzzzz_yyyyy_0, tg_xxxzzzzz_yyyyz_0, tg_xxxzzzzz_yyyzz_0, \
                                         tg_xxxzzzzz_yyzzz_0, tg_xxxzzzzz_yzzzz_0, tg_xxxzzzzz_zzzzz_0, tg_xxyyyyyy_xxxxx_0, \
                                         tg_xxyyyyyy_xxxxy_0, tg_xxyyyyyy_xxxxz_0, tg_xxyyyyyy_xxxyy_0, tg_xxyyyyyy_xxxyz_0, \
                                         tg_xxyyyyyy_xxxzz_0, tg_xxyyyyyy_xxyyy_0, tg_xxyyyyyy_xxyyz_0, tg_xxyyyyyy_xxyzz_0, \
                                         tg_xxyyyyyy_xxzzz_0, tg_xxyyyyyy_xyyyy_0, tg_xxyyyyyy_xyyyz_0, tg_xxyyyyyy_xyyzz_0, \
                                         tg_xxyyyyyy_xyzzz_0, tg_xxyyyyyy_xzzzz_0, tg_xxyyyyyy_yyyyy_0, tg_xxyyyyyy_yyyyz_0, \
                                         tg_xxyyyyyy_yyyzz_0, tg_xxyyyyyy_yyzzz_0, tg_xxyyyyyy_yzzzz_0, tg_xxyyyyyy_zzzzz_0, \
                                         tg_xxyyyyyz_xxxxx_0, tg_xxyyyyyz_xxxxy_0, tg_xxyyyyyz_xxxxz_0, tg_xxyyyyyz_xxxyy_0, \
                                         tg_xxyyyyyz_xxxyz_0, tg_xxyyyyyz_xxxzz_0, tg_xxyyyyyz_xxyyy_0, tg_xxyyyyyz_xxyyz_0, \
                                         tg_xxyyyyyz_xxyzz_0, tg_xxyyyyyz_xxzzz_0, tg_xxyyyyyz_xyyyy_0, tg_xxyyyyyz_xyyyz_0, \
                                         tg_xxyyyyyz_xyyzz_0, tg_xxyyzzz_xxxxz_0, tg_xxyyzzz_xxxxz_1, tg_xxyyzzz_xxxyy_0, \
                                         tg_xxyyzzz_xxxyy_1, tg_xxyyzzz_xxxyz_0, tg_xxyyzzz_xxxyz_1, tg_xxyyzzz_xxxz_1, \
                                         tg_xxyyzzz_xxxzz_0, tg_xxyyzzz_xxxzz_1, tg_xxyyzzz_xxyy_1, tg_xxyyzzz_xxyyy_0, \
                                         tg_xxyyzzz_xxyyy_1, tg_xxyyzzz_xxyyz_0, tg_xxyyzzz_xxyyz_1, tg_xxyyzzz_xxyz_1, \
                                         tg_xxyyzzz_xxyzz_0, tg_xxyyzzz_xxyzz_1, tg_xxyyzzz_xxzz_1, tg_xxyyzzz_xxzzz_0, \
                                         tg_xxyyzzz_xxzzz_1, tg_xxyyzzz_xyyy_1, tg_xxyyzzz_xyyyy_0, tg_xxyyzzz_xyyyy_1, \
                                         tg_xxyyzzz_xyyyz_0, tg_xxyyzzz_xyyyz_1, tg_xxyyzzz_xyyz_1, tg_xxyyzzz_xyyzz_0, \
                                         tg_xxyyzzz_xyyzz_1, tg_xxyyzzz_xyzz_1, tg_xxyyzzz_xyzzz_0, tg_xxyyzzz_xyzzz_1, \
                                         tg_xxyyzzz_xzzz_1, tg_xxyyzzz_xzzzz_0, tg_xxyyzzz_xzzzz_1, tg_xxyyzzz_yyyy_1, \
                                         tg_xxyyzzz_yyyyy_0, tg_xxyyzzz_yyyyy_1, tg_xxyyzzz_yyyyz_0, tg_xxyyzzz_yyyyz_1, \
                                         tg_xxyyzzz_yyyz_1, tg_xxyyzzz_yyyzz_0, tg_xxyyzzz_yyyzz_1, tg_xxyyzzz_yyzz_1, \
                                         tg_xxyyzzz_yyzzz_0, tg_xxyyzzz_yyzzz_1, tg_xxyyzzz_yzzz_1, tg_xxyyzzz_yzzzz_0, \
                                         tg_xxyyzzz_yzzzz_1, tg_xxyyzzz_zzzz_1, tg_xxyyzzz_zzzzz_0, tg_xxyyzzz_zzzzz_1, \
                                         tg_xxyzzzz_xxxx_1, tg_xxyzzzz_xxxxx_0, tg_xxyzzzz_xxxxx_1, tg_xxyzzzz_xxxxy_0, \
                                         tg_xxyzzzz_xxxxy_1, tg_xxyzzzz_xxxxz_0, tg_xxyzzzz_xxxxz_1, tg_xxyzzzz_xxxy_1, \
                                         tg_xxyzzzz_xxxyy_0, tg_xxyzzzz_xxxyy_1, tg_xxyzzzz_xxxyz_0, tg_xxyzzzz_xxxyz_1, \
                                         tg_xxyzzzz_xxxz_1, tg_xxyzzzz_xxxzz_0, tg_xxyzzzz_xxxzz_1, tg_xxyzzzz_xxyy_1, \
                                         tg_xxyzzzz_xxyyy_0, tg_xxyzzzz_xxyyy_1, tg_xxyzzzz_xxyyz_0, tg_xxyzzzz_xxyyz_1, \
                                         tg_xxyzzzz_xxyz_1, tg_xxyzzzz_xxyzz_0, tg_xxyzzzz_xxyzz_1, tg_xxyzzzz_xxzz_1, \
                                         tg_xxyzzzz_xxzzz_0, tg_xxyzzzz_xxzzz_1, tg_xxyzzzz_xyyy_1, tg_xxyzzzz_xyyyy_0, \
                                         tg_xxyzzzz_xyyyy_1, tg_xxyzzzz_xyyyz_0, tg_xxyzzzz_xyyyz_1, tg_xxyzzzz_xyyz_1, \
                                         tg_xxyzzzz_xyyzz_0, tg_xxyzzzz_xyyzz_1, tg_xxyzzzz_xyzz_1, tg_xxyzzzz_xyzzz_0, \
                                         tg_xxyzzzz_xyzzz_1, tg_xxyzzzz_xzzz_1, tg_xxyzzzz_xzzzz_0, tg_xxyzzzz_xzzzz_1, \
                                         tg_xxyzzzz_yyyy_1, tg_xxyzzzz_yyyyy_0, tg_xxyzzzz_yyyyy_1, tg_xxyzzzz_yyyyz_0, \
                                         tg_xxyzzzz_yyyyz_1, tg_xxyzzzz_yyyz_1, tg_xxyzzzz_yyyzz_0, tg_xxyzzzz_yyyzz_1, \
                                         tg_xxyzzzz_yyzz_1, tg_xxyzzzz_yyzzz_0, tg_xxyzzzz_yyzzz_1, tg_xxyzzzz_yzzz_1, \
                                         tg_xxyzzzz_yzzzz_0, tg_xxyzzzz_yzzzz_1, tg_xxyzzzz_zzzz_1, tg_xxyzzzz_zzzzz_0, \
                                         tg_xxyzzzz_zzzzz_1, tg_xxzzzzz_xxxx_1, tg_xxzzzzz_xxxxx_0, tg_xxzzzzz_xxxxx_1, \
                                         tg_xxzzzzz_xxxxy_0, tg_xxzzzzz_xxxxy_1, tg_xxzzzzz_xxxxz_0, tg_xxzzzzz_xxxxz_1, \
                                         tg_xxzzzzz_xxxy_1, tg_xxzzzzz_xxxyy_0, tg_xxzzzzz_xxxyy_1, tg_xxzzzzz_xxxyz_0, \
                                         tg_xxzzzzz_xxxyz_1, tg_xxzzzzz_xxxz_1, tg_xxzzzzz_xxxzz_0, tg_xxzzzzz_xxxzz_1, \
                                         tg_xxzzzzz_xxyy_1, tg_xxzzzzz_xxyyy_0, tg_xxzzzzz_xxyyy_1, tg_xxzzzzz_xxyyz_0, \
                                         tg_xxzzzzz_xxyyz_1, tg_xxzzzzz_xxyz_1, tg_xxzzzzz_xxyzz_0, tg_xxzzzzz_xxyzz_1, \
                                         tg_xxzzzzz_xxzz_1, tg_xxzzzzz_xxzzz_0, tg_xxzzzzz_xxzzz_1, tg_xxzzzzz_xyyy_1, \
                                         tg_xxzzzzz_xyyyy_0, tg_xxzzzzz_xyyyy_1, tg_xxzzzzz_xyyyz_0, tg_xxzzzzz_xyyyz_1, \
                                         tg_xxzzzzz_xyyz_1, tg_xxzzzzz_xyyzz_0, tg_xxzzzzz_xyyzz_1, tg_xxzzzzz_xyzz_1, \
                                         tg_xxzzzzz_xyzzz_0, tg_xxzzzzz_xyzzz_1, tg_xxzzzzz_xzzz_1, tg_xxzzzzz_xzzzz_0, \
                                         tg_xxzzzzz_xzzzz_1, tg_xxzzzzz_yyyy_1, tg_xxzzzzz_yyyyy_0, tg_xxzzzzz_yyyyy_1, \
                                         tg_xxzzzzz_yyyyz_0, tg_xxzzzzz_yyyyz_1, tg_xxzzzzz_yyyz_1, tg_xxzzzzz_yyyzz_0, \
                                         tg_xxzzzzz_yyyzz_1, tg_xxzzzzz_yyzz_1, tg_xxzzzzz_yyzzz_0, tg_xxzzzzz_yyzzz_1, \
                                         tg_xxzzzzz_yzzz_1, tg_xxzzzzz_yzzzz_0, tg_xxzzzzz_yzzzz_1, tg_xxzzzzz_zzzz_1, \
                                         tg_xxzzzzz_zzzzz_0, tg_xxzzzzz_zzzzz_1, tg_xyyyyyy_xxxx_1, tg_xyyyyyy_xxxxx_0, \
                                         tg_xyyyyyy_xxxxx_1, tg_xyyyyyy_xxxxy_0, tg_xyyyyyy_xxxxy_1, tg_xyyyyyy_xxxxz_0, \
                                         tg_xyyyyyy_xxxxz_1, tg_xyyyyyy_xxxy_1, tg_xyyyyyy_xxxyy_0, tg_xyyyyyy_xxxyy_1, \
                                         tg_xyyyyyy_xxxyz_0, tg_xyyyyyy_xxxyz_1, tg_xyyyyyy_xxxz_1, tg_xyyyyyy_xxxzz_0, \
                                         tg_xyyyyyy_xxxzz_1, tg_xyyyyyy_xxyy_1, tg_xyyyyyy_xxyyy_0, tg_xyyyyyy_xxyyy_1, \
                                         tg_xyyyyyy_xxyyz_0, tg_xyyyyyy_xxyyz_1, tg_xyyyyyy_xxyz_1, tg_xyyyyyy_xxyzz_0, \
                                         tg_xyyyyyy_xxyzz_1, tg_xyyyyyy_xxzz_1, tg_xyyyyyy_xxzzz_0, tg_xyyyyyy_xxzzz_1, \
                                         tg_xyyyyyy_xyyy_1, tg_xyyyyyy_xyyyy_0, tg_xyyyyyy_xyyyy_1, tg_xyyyyyy_xyyyz_0, \
                                         tg_xyyyyyy_xyyyz_1, tg_xyyyyyy_xyyz_1, tg_xyyyyyy_xyyzz_0, tg_xyyyyyy_xyyzz_1, \
                                         tg_xyyyyyy_xyzz_1, tg_xyyyyyy_xyzzz_0, tg_xyyyyyy_xyzzz_1, tg_xyyyyyy_xzzz_1, \
                                         tg_xyyyyyy_xzzzz_0, tg_xyyyyyy_xzzzz_1, tg_xyyyyyy_yyyy_1, tg_xyyyyyy_yyyyy_0, \
                                         tg_xyyyyyy_yyyyy_1, tg_xyyyyyy_yyyyz_0, tg_xyyyyyy_yyyyz_1, tg_xyyyyyy_yyyz_1, \
                                         tg_xyyyyyy_yyyzz_0, tg_xyyyyyy_yyyzz_1, tg_xyyyyyy_yyzz_1, tg_xyyyyyy_yyzzz_0, \
                                         tg_xyyyyyy_yyzzz_1, tg_xyyyyyy_yzzz_1, tg_xyyyyyy_yzzzz_0, tg_xyyyyyy_yzzzz_1, \
                                         tg_xyyyyyy_zzzz_1, tg_xyyyyyy_zzzzz_0, tg_xyyyyyy_zzzzz_1, tg_xyyyyyz_xxxx_1, \
                                         tg_xyyyyyz_xxxxx_0, tg_xyyyyyz_xxxxx_1, tg_xyyyyyz_xxxxy_0, tg_xyyyyyz_xxxxy_1, \
                                         tg_xyyyyyz_xxxxz_0, tg_xyyyyyz_xxxxz_1, tg_xyyyyyz_xxxy_1, tg_xyyyyyz_xxxyy_0, \
                                         tg_xyyyyyz_xxxyy_1, tg_xyyyyyz_xxxyz_0, tg_xyyyyyz_xxxyz_1, tg_xyyyyyz_xxxz_1, \
                                         tg_xyyyyyz_xxxzz_0, tg_xyyyyyz_xxxzz_1, tg_xyyyyyz_xxyy_1, tg_xyyyyyz_xxyyy_0, \
                                         tg_xyyyyyz_xxyyy_1, tg_xyyyyyz_xxyyz_0, tg_xyyyyyz_xxyyz_1, tg_xyyyyyz_xxyz_1, \
                                         tg_xyyyyyz_xxyzz_0, tg_xyyyyyz_xxyzz_1, tg_xyyyyyz_xxzz_1, tg_xyyyyyz_xxzzz_0, \
                                         tg_xyyyyyz_xxzzz_1, tg_xyyyyyz_xyyy_1, tg_xyyyyyz_xyyyy_0, tg_xyyyyyz_xyyyy_1, \
                                         tg_xyyyyyz_xyyyz_0, tg_xyyyyyz_xyyyz_1, tg_xyyyyyz_xyyz_1, tg_xyyyyyz_xyyzz_0, \
                                         tg_xyyyyyz_xyyzz_1, tg_xyyyyyz_xyzz_1, tg_xyyyyyz_xzzz_1, tg_xyyyyyz_yyyy_1, \
                                         tg_xyyyyyz_yyyz_1, tg_xyyyyyz_yyzz_1, tg_xyyzzz_xxxxz_0, tg_xyyzzz_xxxxz_1, \
                                         tg_xyyzzz_xxxyy_0, tg_xyyzzz_xxxyy_1, tg_xyyzzz_xxxyz_0, tg_xyyzzz_xxxyz_1, \
                                         tg_xyyzzz_xxxzz_0, tg_xyyzzz_xxxzz_1, tg_xyyzzz_xxyyy_0, tg_xyyzzz_xxyyy_1, \
                                         tg_xyyzzz_xxyyz_0, tg_xyyzzz_xxyyz_1, tg_xyyzzz_xxyzz_0, tg_xyyzzz_xxyzz_1, \
                                         tg_xyyzzz_xxzzz_0, tg_xyyzzz_xxzzz_1, tg_xyyzzz_xyyyy_0, tg_xyyzzz_xyyyy_1, \
                                         tg_xyyzzz_xyyyz_0, tg_xyyzzz_xyyyz_1, tg_xyyzzz_xyyzz_0, tg_xyyzzz_xyyzz_1, \
                                         tg_xyyzzz_xyzzz_0, tg_xyyzzz_xyzzz_1, tg_xyyzzz_xzzzz_0, tg_xyyzzz_xzzzz_1, \
                                         tg_xyyzzz_yyyyy_0, tg_xyyzzz_yyyyy_1, tg_xyyzzz_yyyyz_0, tg_xyyzzz_yyyyz_1, \
                                         tg_xyyzzz_yyyzz_0, tg_xyyzzz_yyyzz_1, tg_xyyzzz_yyzzz_0, tg_xyyzzz_yyzzz_1, \
                                         tg_xyyzzz_yzzzz_0, tg_xyyzzz_yzzzz_1, tg_xyyzzz_zzzzz_0, tg_xyyzzz_zzzzz_1, \
                                         tg_xyzzzz_xxxxx_0, tg_xyzzzz_xxxxx_1, tg_xyzzzz_xxxxy_0, tg_xyzzzz_xxxxy_1, \
                                         tg_xyzzzz_xxxxz_0, tg_xyzzzz_xxxxz_1, tg_xyzzzz_xxxyy_0, tg_xyzzzz_xxxyy_1, \
                                         tg_xyzzzz_xxxyz_0, tg_xyzzzz_xxxyz_1, tg_xyzzzz_xxxzz_0, tg_xyzzzz_xxxzz_1, \
                                         tg_xyzzzz_xxyyy_0, tg_xyzzzz_xxyyy_1, tg_xyzzzz_xxyyz_0, tg_xyzzzz_xxyyz_1, \
                                         tg_xyzzzz_xxyzz_0, tg_xyzzzz_xxyzz_1, tg_xyzzzz_xxzzz_0, tg_xyzzzz_xxzzz_1, \
                                         tg_xyzzzz_xyyyy_0, tg_xyzzzz_xyyyy_1, tg_xyzzzz_xyyyz_0, tg_xyzzzz_xyyyz_1, \
                                         tg_xyzzzz_xyyzz_0, tg_xyzzzz_xyyzz_1, tg_xyzzzz_xyzzz_0, tg_xyzzzz_xyzzz_1, \
                                         tg_xyzzzz_xzzzz_0, tg_xyzzzz_xzzzz_1, tg_xyzzzz_yyyyy_0, tg_xyzzzz_yyyyy_1, \
                                         tg_xyzzzz_yyyyz_0, tg_xyzzzz_yyyyz_1, tg_xyzzzz_yyyzz_0, tg_xyzzzz_yyyzz_1, \
                                         tg_xyzzzz_yyzzz_0, tg_xyzzzz_yyzzz_1, tg_xyzzzz_yzzzz_0, tg_xyzzzz_yzzzz_1, \
                                         tg_xyzzzz_zzzzz_0, tg_xyzzzz_zzzzz_1, tg_xzzzzz_xxxxx_0, tg_xzzzzz_xxxxx_1, \
                                         tg_xzzzzz_xxxxy_0, tg_xzzzzz_xxxxy_1, tg_xzzzzz_xxxxz_0, tg_xzzzzz_xxxxz_1, \
                                         tg_xzzzzz_xxxyy_0, tg_xzzzzz_xxxyy_1, tg_xzzzzz_xxxyz_0, tg_xzzzzz_xxxyz_1, \
                                         tg_xzzzzz_xxxzz_0, tg_xzzzzz_xxxzz_1, tg_xzzzzz_xxyyy_0, tg_xzzzzz_xxyyy_1, \
                                         tg_xzzzzz_xxyyz_0, tg_xzzzzz_xxyyz_1, tg_xzzzzz_xxyzz_0, tg_xzzzzz_xxyzz_1, \
                                         tg_xzzzzz_xxzzz_0, tg_xzzzzz_xxzzz_1, tg_xzzzzz_xyyyy_0, tg_xzzzzz_xyyyy_1, \
                                         tg_xzzzzz_xyyyz_0, tg_xzzzzz_xyyyz_1, tg_xzzzzz_xyyzz_0, tg_xzzzzz_xyyzz_1, \
                                         tg_xzzzzz_xyzzz_0, tg_xzzzzz_xyzzz_1, tg_xzzzzz_xzzzz_0, tg_xzzzzz_xzzzz_1, \
                                         tg_xzzzzz_yyyyy_0, tg_xzzzzz_yyyyy_1, tg_xzzzzz_yyyyz_0, tg_xzzzzz_yyyyz_1, \
                                         tg_xzzzzz_yyyzz_0, tg_xzzzzz_yyyzz_1, tg_xzzzzz_yyzzz_0, tg_xzzzzz_yyzzz_1, \
                                         tg_xzzzzz_yzzzz_0, tg_xzzzzz_yzzzz_1, tg_xzzzzz_zzzzz_0, tg_xzzzzz_zzzzz_1, \
                                         tg_yyyyyy_xxxxx_0, tg_yyyyyy_xxxxx_1, tg_yyyyyy_xxxxy_0, tg_yyyyyy_xxxxy_1, \
                                         tg_yyyyyy_xxxxz_0, tg_yyyyyy_xxxxz_1, tg_yyyyyy_xxxyy_0, tg_yyyyyy_xxxyy_1, \
                                         tg_yyyyyy_xxxyz_0, tg_yyyyyy_xxxyz_1, tg_yyyyyy_xxxzz_0, tg_yyyyyy_xxxzz_1, \
                                         tg_yyyyyy_xxyyy_0, tg_yyyyyy_xxyyy_1, tg_yyyyyy_xxyyz_0, tg_yyyyyy_xxyyz_1, \
                                         tg_yyyyyy_xxyzz_0, tg_yyyyyy_xxyzz_1, tg_yyyyyy_xxzzz_0, tg_yyyyyy_xxzzz_1, \
                                         tg_yyyyyy_xyyyy_0, tg_yyyyyy_xyyyy_1, tg_yyyyyy_xyyyz_0, tg_yyyyyy_xyyyz_1, \
                                         tg_yyyyyy_xyyzz_0, tg_yyyyyy_xyyzz_1, tg_yyyyyy_xyzzz_0, tg_yyyyyy_xyzzz_1, \
                                         tg_yyyyyy_xzzzz_0, tg_yyyyyy_xzzzz_1, tg_yyyyyy_yyyyy_0, tg_yyyyyy_yyyyy_1, \
                                         tg_yyyyyy_yyyyz_0, tg_yyyyyy_yyyyz_1, tg_yyyyyy_yyyzz_0, tg_yyyyyy_yyyzz_1, \
                                         tg_yyyyyy_yyzzz_0, tg_yyyyyy_yyzzz_1, tg_yyyyyy_yzzzz_0, tg_yyyyyy_yzzzz_1, \
                                         tg_yyyyyy_zzzzz_0, tg_yyyyyy_zzzzz_1, tg_yyyyyz_xxxxx_0, tg_yyyyyz_xxxxx_1, \
                                         tg_yyyyyz_xxxxy_0, tg_yyyyyz_xxxxy_1, tg_yyyyyz_xxxxz_0, tg_yyyyyz_xxxxz_1, \
                                         tg_yyyyyz_xxxyy_0, tg_yyyyyz_xxxyy_1, tg_yyyyyz_xxxyz_0, tg_yyyyyz_xxxyz_1, \
                                         tg_yyyyyz_xxxzz_0, tg_yyyyyz_xxxzz_1, tg_yyyyyz_xxyyy_0, tg_yyyyyz_xxyyy_1, \
                                         tg_yyyyyz_xxyyz_0, tg_yyyyyz_xxyyz_1, tg_yyyyyz_xxyzz_0, tg_yyyyyz_xxyzz_1, \
                                         tg_yyyyyz_xxzzz_0, tg_yyyyyz_xxzzz_1, tg_yyyyyz_xyyyy_0, tg_yyyyyz_xyyyy_1, \
                                         tg_yyyyyz_xyyyz_0, tg_yyyyyz_xyyyz_1, tg_yyyyyz_xyyzz_0, tg_yyyyyz_xyyzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyyzzz_xxxxz_0[j] = pb_x * tg_xxyyzzz_xxxxz_0[j] + fr * tg_xxyyzzz_xxxxz_1[j] + fl1_fx * (tg_xyyzzz_xxxxz_0[j] - tg_xyyzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzzz_xxxz_1[j];

                    tg_xxxyyzzz_xxxyy_0[j] = pb_x * tg_xxyyzzz_xxxyy_0[j] + fr * tg_xxyyzzz_xxxyy_1[j] + fl1_fx * (tg_xyyzzz_xxxyy_0[j] - tg_xyyzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzzz_xxyy_1[j];

                    tg_xxxyyzzz_xxxyz_0[j] = pb_x * tg_xxyyzzz_xxxyz_0[j] + fr * tg_xxyyzzz_xxxyz_1[j] + fl1_fx * (tg_xyyzzz_xxxyz_0[j] - tg_xyyzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzzz_xxyz_1[j];

                    tg_xxxyyzzz_xxxzz_0[j] = pb_x * tg_xxyyzzz_xxxzz_0[j] + fr * tg_xxyyzzz_xxxzz_1[j] + fl1_fx * (tg_xyyzzz_xxxzz_0[j] - tg_xyyzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzzz_xxzz_1[j];

                    tg_xxxyyzzz_xxyyy_0[j] = pb_x * tg_xxyyzzz_xxyyy_0[j] + fr * tg_xxyyzzz_xxyyy_1[j] + fl1_fx * (tg_xyyzzz_xxyyy_0[j] - tg_xyyzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzzz_xyyy_1[j];

                    tg_xxxyyzzz_xxyyz_0[j] = pb_x * tg_xxyyzzz_xxyyz_0[j] + fr * tg_xxyyzzz_xxyyz_1[j] + fl1_fx * (tg_xyyzzz_xxyyz_0[j] - tg_xyyzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzzz_xyyz_1[j];

                    tg_xxxyyzzz_xxyzz_0[j] = pb_x * tg_xxyyzzz_xxyzz_0[j] + fr * tg_xxyyzzz_xxyzz_1[j] + fl1_fx * (tg_xyyzzz_xxyzz_0[j] - tg_xyyzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzzz_xyzz_1[j];

                    tg_xxxyyzzz_xxzzz_0[j] = pb_x * tg_xxyyzzz_xxzzz_0[j] + fr * tg_xxyyzzz_xxzzz_1[j] + fl1_fx * (tg_xyyzzz_xxzzz_0[j] - tg_xyyzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzzz_xzzz_1[j];

                    tg_xxxyyzzz_xyyyy_0[j] = pb_x * tg_xxyyzzz_xyyyy_0[j] + fr * tg_xxyyzzz_xyyyy_1[j] + fl1_fx * (tg_xyyzzz_xyyyy_0[j] - tg_xyyzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_yyyy_1[j];

                    tg_xxxyyzzz_xyyyz_0[j] = pb_x * tg_xxyyzzz_xyyyz_0[j] + fr * tg_xxyyzzz_xyyyz_1[j] + fl1_fx * (tg_xyyzzz_xyyyz_0[j] - tg_xyyzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_yyyz_1[j];

                    tg_xxxyyzzz_xyyzz_0[j] = pb_x * tg_xxyyzzz_xyyzz_0[j] + fr * tg_xxyyzzz_xyyzz_1[j] + fl1_fx * (tg_xyyzzz_xyyzz_0[j] - tg_xyyzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_yyzz_1[j];

                    tg_xxxyyzzz_xyzzz_0[j] = pb_x * tg_xxyyzzz_xyzzz_0[j] + fr * tg_xxyyzzz_xyzzz_1[j] + fl1_fx * (tg_xyyzzz_xyzzz_0[j] - tg_xyyzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_yzzz_1[j];

                    tg_xxxyyzzz_xzzzz_0[j] = pb_x * tg_xxyyzzz_xzzzz_0[j] + fr * tg_xxyyzzz_xzzzz_1[j] + fl1_fx * (tg_xyyzzz_xzzzz_0[j] - tg_xyyzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_zzzz_1[j];

                    tg_xxxyyzzz_yyyyy_0[j] = pb_x * tg_xxyyzzz_yyyyy_0[j] + fr * tg_xxyyzzz_yyyyy_1[j] + fl1_fx * (tg_xyyzzz_yyyyy_0[j] - tg_xyyzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyyzzz_yyyyz_0[j] = pb_x * tg_xxyyzzz_yyyyz_0[j] + fr * tg_xxyyzzz_yyyyz_1[j] + fl1_fx * (tg_xyyzzz_yyyyz_0[j] - tg_xyyzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyyzzz_yyyzz_0[j] = pb_x * tg_xxyyzzz_yyyzz_0[j] + fr * tg_xxyyzzz_yyyzz_1[j] + fl1_fx * (tg_xyyzzz_yyyzz_0[j] - tg_xyyzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyyzzz_yyzzz_0[j] = pb_x * tg_xxyyzzz_yyzzz_0[j] + fr * tg_xxyyzzz_yyzzz_1[j] + fl1_fx * (tg_xyyzzz_yyzzz_0[j] - tg_xyyzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyyzzz_yzzzz_0[j] = pb_x * tg_xxyyzzz_yzzzz_0[j] + fr * tg_xxyyzzz_yzzzz_1[j] + fl1_fx * (tg_xyyzzz_yzzzz_0[j] - tg_xyyzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyyzzz_zzzzz_0[j] = pb_x * tg_xxyyzzz_zzzzz_0[j] + fr * tg_xxyyzzz_zzzzz_1[j] + fl1_fx * (tg_xyyzzz_zzzzz_0[j] - tg_xyyzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyzzzz_xxxxx_0[j] = pb_x * tg_xxyzzzz_xxxxx_0[j] + fr * tg_xxyzzzz_xxxxx_1[j] + fl1_fx * (tg_xyzzzz_xxxxx_0[j] - tg_xyzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyzzzz_xxxx_1[j];

                    tg_xxxyzzzz_xxxxy_0[j] = pb_x * tg_xxyzzzz_xxxxy_0[j] + fr * tg_xxyzzzz_xxxxy_1[j] + fl1_fx * (tg_xyzzzz_xxxxy_0[j] - tg_xyzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzzz_xxxy_1[j];

                    tg_xxxyzzzz_xxxxz_0[j] = pb_x * tg_xxyzzzz_xxxxz_0[j] + fr * tg_xxyzzzz_xxxxz_1[j] + fl1_fx * (tg_xyzzzz_xxxxz_0[j] - tg_xyzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzzz_xxxz_1[j];

                    tg_xxxyzzzz_xxxyy_0[j] = pb_x * tg_xxyzzzz_xxxyy_0[j] + fr * tg_xxyzzzz_xxxyy_1[j] + fl1_fx * (tg_xyzzzz_xxxyy_0[j] - tg_xyzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzzz_xxyy_1[j];

                    tg_xxxyzzzz_xxxyz_0[j] = pb_x * tg_xxyzzzz_xxxyz_0[j] + fr * tg_xxyzzzz_xxxyz_1[j] + fl1_fx * (tg_xyzzzz_xxxyz_0[j] - tg_xyzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzzz_xxyz_1[j];

                    tg_xxxyzzzz_xxxzz_0[j] = pb_x * tg_xxyzzzz_xxxzz_0[j] + fr * tg_xxyzzzz_xxxzz_1[j] + fl1_fx * (tg_xyzzzz_xxxzz_0[j] - tg_xyzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzzz_xxzz_1[j];

                    tg_xxxyzzzz_xxyyy_0[j] = pb_x * tg_xxyzzzz_xxyyy_0[j] + fr * tg_xxyzzzz_xxyyy_1[j] + fl1_fx * (tg_xyzzzz_xxyyy_0[j] - tg_xyzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzzz_xyyy_1[j];

                    tg_xxxyzzzz_xxyyz_0[j] = pb_x * tg_xxyzzzz_xxyyz_0[j] + fr * tg_xxyzzzz_xxyyz_1[j] + fl1_fx * (tg_xyzzzz_xxyyz_0[j] - tg_xyzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzzz_xyyz_1[j];

                    tg_xxxyzzzz_xxyzz_0[j] = pb_x * tg_xxyzzzz_xxyzz_0[j] + fr * tg_xxyzzzz_xxyzz_1[j] + fl1_fx * (tg_xyzzzz_xxyzz_0[j] - tg_xyzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzzz_xyzz_1[j];

                    tg_xxxyzzzz_xxzzz_0[j] = pb_x * tg_xxyzzzz_xxzzz_0[j] + fr * tg_xxyzzzz_xxzzz_1[j] + fl1_fx * (tg_xyzzzz_xxzzz_0[j] - tg_xyzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzzz_xzzz_1[j];

                    tg_xxxyzzzz_xyyyy_0[j] = pb_x * tg_xxyzzzz_xyyyy_0[j] + fr * tg_xxyzzzz_xyyyy_1[j] + fl1_fx * (tg_xyzzzz_xyyyy_0[j] - tg_xyzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_yyyy_1[j];

                    tg_xxxyzzzz_xyyyz_0[j] = pb_x * tg_xxyzzzz_xyyyz_0[j] + fr * tg_xxyzzzz_xyyyz_1[j] + fl1_fx * (tg_xyzzzz_xyyyz_0[j] - tg_xyzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_yyyz_1[j];

                    tg_xxxyzzzz_xyyzz_0[j] = pb_x * tg_xxyzzzz_xyyzz_0[j] + fr * tg_xxyzzzz_xyyzz_1[j] + fl1_fx * (tg_xyzzzz_xyyzz_0[j] - tg_xyzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_yyzz_1[j];

                    tg_xxxyzzzz_xyzzz_0[j] = pb_x * tg_xxyzzzz_xyzzz_0[j] + fr * tg_xxyzzzz_xyzzz_1[j] + fl1_fx * (tg_xyzzzz_xyzzz_0[j] - tg_xyzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_yzzz_1[j];

                    tg_xxxyzzzz_xzzzz_0[j] = pb_x * tg_xxyzzzz_xzzzz_0[j] + fr * tg_xxyzzzz_xzzzz_1[j] + fl1_fx * (tg_xyzzzz_xzzzz_0[j] - tg_xyzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_zzzz_1[j];

                    tg_xxxyzzzz_yyyyy_0[j] = pb_x * tg_xxyzzzz_yyyyy_0[j] + fr * tg_xxyzzzz_yyyyy_1[j] + fl1_fx * (tg_xyzzzz_yyyyy_0[j] - tg_xyzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyzzzz_yyyyz_0[j] = pb_x * tg_xxyzzzz_yyyyz_0[j] + fr * tg_xxyzzzz_yyyyz_1[j] + fl1_fx * (tg_xyzzzz_yyyyz_0[j] - tg_xyzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyzzzz_yyyzz_0[j] = pb_x * tg_xxyzzzz_yyyzz_0[j] + fr * tg_xxyzzzz_yyyzz_1[j] + fl1_fx * (tg_xyzzzz_yyyzz_0[j] - tg_xyzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyzzzz_yyzzz_0[j] = pb_x * tg_xxyzzzz_yyzzz_0[j] + fr * tg_xxyzzzz_yyzzz_1[j] + fl1_fx * (tg_xyzzzz_yyzzz_0[j] - tg_xyzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyzzzz_yzzzz_0[j] = pb_x * tg_xxyzzzz_yzzzz_0[j] + fr * tg_xxyzzzz_yzzzz_1[j] + fl1_fx * (tg_xyzzzz_yzzzz_0[j] - tg_xyzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyzzzz_zzzzz_0[j] = pb_x * tg_xxyzzzz_zzzzz_0[j] + fr * tg_xxyzzzz_zzzzz_1[j] + fl1_fx * (tg_xyzzzz_zzzzz_0[j] - tg_xyzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxxzzzzz_xxxxx_0[j] = pb_x * tg_xxzzzzz_xxxxx_0[j] + fr * tg_xxzzzzz_xxxxx_1[j] + fl1_fx * (tg_xzzzzz_xxxxx_0[j] - tg_xzzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzzzzz_xxxx_1[j];

                    tg_xxxzzzzz_xxxxy_0[j] = pb_x * tg_xxzzzzz_xxxxy_0[j] + fr * tg_xxzzzzz_xxxxy_1[j] + fl1_fx * (tg_xzzzzz_xxxxy_0[j] - tg_xzzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzzz_xxxy_1[j];

                    tg_xxxzzzzz_xxxxz_0[j] = pb_x * tg_xxzzzzz_xxxxz_0[j] + fr * tg_xxzzzzz_xxxxz_1[j] + fl1_fx * (tg_xzzzzz_xxxxz_0[j] - tg_xzzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzzz_xxxz_1[j];

                    tg_xxxzzzzz_xxxyy_0[j] = pb_x * tg_xxzzzzz_xxxyy_0[j] + fr * tg_xxzzzzz_xxxyy_1[j] + fl1_fx * (tg_xzzzzz_xxxyy_0[j] - tg_xzzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzzz_xxyy_1[j];

                    tg_xxxzzzzz_xxxyz_0[j] = pb_x * tg_xxzzzzz_xxxyz_0[j] + fr * tg_xxzzzzz_xxxyz_1[j] + fl1_fx * (tg_xzzzzz_xxxyz_0[j] - tg_xzzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzzz_xxyz_1[j];

                    tg_xxxzzzzz_xxxzz_0[j] = pb_x * tg_xxzzzzz_xxxzz_0[j] + fr * tg_xxzzzzz_xxxzz_1[j] + fl1_fx * (tg_xzzzzz_xxxzz_0[j] - tg_xzzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzzz_xxzz_1[j];

                    tg_xxxzzzzz_xxyyy_0[j] = pb_x * tg_xxzzzzz_xxyyy_0[j] + fr * tg_xxzzzzz_xxyyy_1[j] + fl1_fx * (tg_xzzzzz_xxyyy_0[j] - tg_xzzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzzz_xyyy_1[j];

                    tg_xxxzzzzz_xxyyz_0[j] = pb_x * tg_xxzzzzz_xxyyz_0[j] + fr * tg_xxzzzzz_xxyyz_1[j] + fl1_fx * (tg_xzzzzz_xxyyz_0[j] - tg_xzzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzzz_xyyz_1[j];

                    tg_xxxzzzzz_xxyzz_0[j] = pb_x * tg_xxzzzzz_xxyzz_0[j] + fr * tg_xxzzzzz_xxyzz_1[j] + fl1_fx * (tg_xzzzzz_xxyzz_0[j] - tg_xzzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzzz_xyzz_1[j];

                    tg_xxxzzzzz_xxzzz_0[j] = pb_x * tg_xxzzzzz_xxzzz_0[j] + fr * tg_xxzzzzz_xxzzz_1[j] + fl1_fx * (tg_xzzzzz_xxzzz_0[j] - tg_xzzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzzz_xzzz_1[j];

                    tg_xxxzzzzz_xyyyy_0[j] = pb_x * tg_xxzzzzz_xyyyy_0[j] + fr * tg_xxzzzzz_xyyyy_1[j] + fl1_fx * (tg_xzzzzz_xyyyy_0[j] - tg_xzzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_yyyy_1[j];

                    tg_xxxzzzzz_xyyyz_0[j] = pb_x * tg_xxzzzzz_xyyyz_0[j] + fr * tg_xxzzzzz_xyyyz_1[j] + fl1_fx * (tg_xzzzzz_xyyyz_0[j] - tg_xzzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_yyyz_1[j];

                    tg_xxxzzzzz_xyyzz_0[j] = pb_x * tg_xxzzzzz_xyyzz_0[j] + fr * tg_xxzzzzz_xyyzz_1[j] + fl1_fx * (tg_xzzzzz_xyyzz_0[j] - tg_xzzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_yyzz_1[j];

                    tg_xxxzzzzz_xyzzz_0[j] = pb_x * tg_xxzzzzz_xyzzz_0[j] + fr * tg_xxzzzzz_xyzzz_1[j] + fl1_fx * (tg_xzzzzz_xyzzz_0[j] - tg_xzzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_yzzz_1[j];

                    tg_xxxzzzzz_xzzzz_0[j] = pb_x * tg_xxzzzzz_xzzzz_0[j] + fr * tg_xxzzzzz_xzzzz_1[j] + fl1_fx * (tg_xzzzzz_xzzzz_0[j] - tg_xzzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_zzzz_1[j];

                    tg_xxxzzzzz_yyyyy_0[j] = pb_x * tg_xxzzzzz_yyyyy_0[j] + fr * tg_xxzzzzz_yyyyy_1[j] + fl1_fx * (tg_xzzzzz_yyyyy_0[j] - tg_xzzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxzzzzz_yyyyz_0[j] = pb_x * tg_xxzzzzz_yyyyz_0[j] + fr * tg_xxzzzzz_yyyyz_1[j] + fl1_fx * (tg_xzzzzz_yyyyz_0[j] - tg_xzzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxzzzzz_yyyzz_0[j] = pb_x * tg_xxzzzzz_yyyzz_0[j] + fr * tg_xxzzzzz_yyyzz_1[j] + fl1_fx * (tg_xzzzzz_yyyzz_0[j] - tg_xzzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxzzzzz_yyzzz_0[j] = pb_x * tg_xxzzzzz_yyzzz_0[j] + fr * tg_xxzzzzz_yyzzz_1[j] + fl1_fx * (tg_xzzzzz_yyzzz_0[j] - tg_xzzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxzzzzz_yzzzz_0[j] = pb_x * tg_xxzzzzz_yzzzz_0[j] + fr * tg_xxzzzzz_yzzzz_1[j] + fl1_fx * (tg_xzzzzz_yzzzz_0[j] - tg_xzzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxzzzzz_zzzzz_0[j] = pb_x * tg_xxzzzzz_zzzzz_0[j] + fr * tg_xxzzzzz_zzzzz_1[j] + fl1_fx * (tg_xzzzzz_zzzzz_0[j] - tg_xzzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyyyy_xxxxx_0[j] = pb_x * tg_xyyyyyy_xxxxx_0[j] + fr * tg_xyyyyyy_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxxxx_0[j] - tg_yyyyyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyyy_xxxx_1[j];

                    tg_xxyyyyyy_xxxxy_0[j] = pb_x * tg_xyyyyyy_xxxxy_0[j] + fr * tg_xyyyyyy_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxxxy_0[j] - tg_yyyyyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyyy_xxxy_1[j];

                    tg_xxyyyyyy_xxxxz_0[j] = pb_x * tg_xyyyyyy_xxxxz_0[j] + fr * tg_xyyyyyy_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxxxz_0[j] - tg_yyyyyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyyy_xxxz_1[j];

                    tg_xxyyyyyy_xxxyy_0[j] = pb_x * tg_xyyyyyy_xxxyy_0[j] + fr * tg_xyyyyyy_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxxyy_0[j] - tg_yyyyyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyy_xxyy_1[j];

                    tg_xxyyyyyy_xxxyz_0[j] = pb_x * tg_xyyyyyy_xxxyz_0[j] + fr * tg_xyyyyyy_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxxyz_0[j] - tg_yyyyyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyy_xxyz_1[j];

                    tg_xxyyyyyy_xxxzz_0[j] = pb_x * tg_xyyyyyy_xxxzz_0[j] + fr * tg_xyyyyyy_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxxzz_0[j] - tg_yyyyyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyy_xxzz_1[j];

                    tg_xxyyyyyy_xxyyy_0[j] = pb_x * tg_xyyyyyy_xxyyy_0[j] + fr * tg_xyyyyyy_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxyyy_0[j] - tg_yyyyyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyy_xyyy_1[j];

                    tg_xxyyyyyy_xxyyz_0[j] = pb_x * tg_xyyyyyy_xxyyz_0[j] + fr * tg_xyyyyyy_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxyyz_0[j] - tg_yyyyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyy_xyyz_1[j];

                    tg_xxyyyyyy_xxyzz_0[j] = pb_x * tg_xyyyyyy_xxyzz_0[j] + fr * tg_xyyyyyy_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxyzz_0[j] - tg_yyyyyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyy_xyzz_1[j];

                    tg_xxyyyyyy_xxzzz_0[j] = pb_x * tg_xyyyyyy_xxzzz_0[j] + fr * tg_xyyyyyy_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxzzz_0[j] - tg_yyyyyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyy_xzzz_1[j];

                    tg_xxyyyyyy_xyyyy_0[j] = pb_x * tg_xyyyyyy_xyyyy_0[j] + fr * tg_xyyyyyy_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xyyyy_0[j] - tg_yyyyyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_yyyy_1[j];

                    tg_xxyyyyyy_xyyyz_0[j] = pb_x * tg_xyyyyyy_xyyyz_0[j] + fr * tg_xyyyyyy_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xyyyz_0[j] - tg_yyyyyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_yyyz_1[j];

                    tg_xxyyyyyy_xyyzz_0[j] = pb_x * tg_xyyyyyy_xyyzz_0[j] + fr * tg_xyyyyyy_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xyyzz_0[j] - tg_yyyyyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_yyzz_1[j];

                    tg_xxyyyyyy_xyzzz_0[j] = pb_x * tg_xyyyyyy_xyzzz_0[j] + fr * tg_xyyyyyy_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xyzzz_0[j] - tg_yyyyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_yzzz_1[j];

                    tg_xxyyyyyy_xzzzz_0[j] = pb_x * tg_xyyyyyy_xzzzz_0[j] + fr * tg_xyyyyyy_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xzzzz_0[j] - tg_yyyyyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_zzzz_1[j];

                    tg_xxyyyyyy_yyyyy_0[j] = pb_x * tg_xyyyyyy_yyyyy_0[j] + fr * tg_xyyyyyy_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yyyyy_0[j] - tg_yyyyyy_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyyyy_yyyyz_0[j] = pb_x * tg_xyyyyyy_yyyyz_0[j] + fr * tg_xyyyyyy_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yyyyz_0[j] - tg_yyyyyy_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyyyy_yyyzz_0[j] = pb_x * tg_xyyyyyy_yyyzz_0[j] + fr * tg_xyyyyyy_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yyyzz_0[j] - tg_yyyyyy_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyyyy_yyzzz_0[j] = pb_x * tg_xyyyyyy_yyzzz_0[j] + fr * tg_xyyyyyy_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yyzzz_0[j] - tg_yyyyyy_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyyyy_yzzzz_0[j] = pb_x * tg_xyyyyyy_yzzzz_0[j] + fr * tg_xyyyyyy_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yzzzz_0[j] - tg_yyyyyy_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyyyy_zzzzz_0[j] = pb_x * tg_xyyyyyy_zzzzz_0[j] + fr * tg_xyyyyyy_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_zzzzz_0[j] - tg_yyyyyy_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyyyz_xxxxx_0[j] = pb_x * tg_xyyyyyz_xxxxx_0[j] + fr * tg_xyyyyyz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxxxx_0[j] - tg_yyyyyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyyz_xxxx_1[j];

                    tg_xxyyyyyz_xxxxy_0[j] = pb_x * tg_xyyyyyz_xxxxy_0[j] + fr * tg_xyyyyyz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxxxy_0[j] - tg_yyyyyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyyz_xxxy_1[j];

                    tg_xxyyyyyz_xxxxz_0[j] = pb_x * tg_xyyyyyz_xxxxz_0[j] + fr * tg_xyyyyyz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxxxz_0[j] - tg_yyyyyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyyz_xxxz_1[j];

                    tg_xxyyyyyz_xxxyy_0[j] = pb_x * tg_xyyyyyz_xxxyy_0[j] + fr * tg_xyyyyyz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxxyy_0[j] - tg_yyyyyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyz_xxyy_1[j];

                    tg_xxyyyyyz_xxxyz_0[j] = pb_x * tg_xyyyyyz_xxxyz_0[j] + fr * tg_xyyyyyz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxxyz_0[j] - tg_yyyyyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyz_xxyz_1[j];

                    tg_xxyyyyyz_xxxzz_0[j] = pb_x * tg_xyyyyyz_xxxzz_0[j] + fr * tg_xyyyyyz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxxzz_0[j] - tg_yyyyyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyz_xxzz_1[j];

                    tg_xxyyyyyz_xxyyy_0[j] = pb_x * tg_xyyyyyz_xxyyy_0[j] + fr * tg_xyyyyyz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxyyy_0[j] - tg_yyyyyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyz_xyyy_1[j];

                    tg_xxyyyyyz_xxyyz_0[j] = pb_x * tg_xyyyyyz_xxyyz_0[j] + fr * tg_xyyyyyz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxyyz_0[j] - tg_yyyyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyz_xyyz_1[j];

                    tg_xxyyyyyz_xxyzz_0[j] = pb_x * tg_xyyyyyz_xxyzz_0[j] + fr * tg_xyyyyyz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxyzz_0[j] - tg_yyyyyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyz_xyzz_1[j];

                    tg_xxyyyyyz_xxzzz_0[j] = pb_x * tg_xyyyyyz_xxzzz_0[j] + fr * tg_xyyyyyz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxzzz_0[j] - tg_yyyyyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyz_xzzz_1[j];

                    tg_xxyyyyyz_xyyyy_0[j] = pb_x * tg_xyyyyyz_xyyyy_0[j] + fr * tg_xyyyyyz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xyyyy_0[j] - tg_yyyyyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_yyyy_1[j];

                    tg_xxyyyyyz_xyyyz_0[j] = pb_x * tg_xyyyyyz_xyyyz_0[j] + fr * tg_xyyyyyz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xyyyz_0[j] - tg_yyyyyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_yyyz_1[j];

                    tg_xxyyyyyz_xyyzz_0[j] = pb_x * tg_xyyyyyz_xyyzz_0[j] + fr * tg_xyyyyyz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xyyzz_0[j] - tg_yyyyyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_yyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_475_569(      CMemBlock2D<double>* primBuffer,
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

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xzzzzzz_xxxxy_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 568); 

                auto tg_xyyyyyz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 475); 

                auto tg_xyyyyyz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 476); 

                auto tg_xyyyyyz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 477); 

                auto tg_xyyyyyz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 478); 

                auto tg_xyyyyyz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 479); 

                auto tg_xyyyyyz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 480); 

                auto tg_xyyyyyz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 481); 

                auto tg_xyyyyyz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 482); 

                auto tg_xyyyyzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 483); 

                auto tg_xyyyyzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 484); 

                auto tg_xyyyyzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 485); 

                auto tg_xyyyyzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 486); 

                auto tg_xyyyyzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 487); 

                auto tg_xyyyyzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 488); 

                auto tg_xyyyyzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 489); 

                auto tg_xyyyyzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 490); 

                auto tg_xyyyyzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 491); 

                auto tg_xyyyyzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 492); 

                auto tg_xyyyyzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 493); 

                auto tg_xyyyyzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 494); 

                auto tg_xyyyyzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 495); 

                auto tg_xyyyyzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 496); 

                auto tg_xyyyyzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 497); 

                auto tg_xyyyyzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 498); 

                auto tg_xyyyyzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 499); 

                auto tg_xyyyyzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 500); 

                auto tg_xyyyyzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 501); 

                auto tg_xyyyyzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 502); 

                auto tg_xyyyyzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 503); 

                auto tg_xyyyzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 504); 

                auto tg_xyyyzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 505); 

                auto tg_xyyyzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 506); 

                auto tg_xyyyzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 507); 

                auto tg_xyyyzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 508); 

                auto tg_xyyyzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 509); 

                auto tg_xyyyzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 510); 

                auto tg_xyyyzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 511); 

                auto tg_xyyyzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 512); 

                auto tg_xyyyzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 513); 

                auto tg_xyyyzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 514); 

                auto tg_xyyyzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 515); 

                auto tg_xyyyzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 516); 

                auto tg_xyyyzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 517); 

                auto tg_xyyyzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 518); 

                auto tg_xyyyzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 519); 

                auto tg_xyyyzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 520); 

                auto tg_xyyyzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 521); 

                auto tg_xyyyzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 522); 

                auto tg_xyyyzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 523); 

                auto tg_xyyyzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 524); 

                auto tg_xyyzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 525); 

                auto tg_xyyzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 526); 

                auto tg_xyyzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 527); 

                auto tg_xyyzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 528); 

                auto tg_xyyzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 529); 

                auto tg_xyyzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 530); 

                auto tg_xyyzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 531); 

                auto tg_xyyzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 532); 

                auto tg_xyyzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 533); 

                auto tg_xyyzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 534); 

                auto tg_xyyzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 535); 

                auto tg_xyyzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 536); 

                auto tg_xyyzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 537); 

                auto tg_xyyzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 538); 

                auto tg_xyyzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 539); 

                auto tg_xyyzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 540); 

                auto tg_xyyzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 541); 

                auto tg_xyyzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 542); 

                auto tg_xyyzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 543); 

                auto tg_xyyzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 544); 

                auto tg_xyyzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 545); 

                auto tg_xyzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 546); 

                auto tg_xyzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 547); 

                auto tg_xyzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 548); 

                auto tg_xyzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 549); 

                auto tg_xyzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 550); 

                auto tg_xyzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 551); 

                auto tg_xyzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 552); 

                auto tg_xyzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 553); 

                auto tg_xyzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 554); 

                auto tg_xyzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 555); 

                auto tg_xyzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 556); 

                auto tg_xyzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 557); 

                auto tg_xyzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 558); 

                auto tg_xyzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 559); 

                auto tg_xyzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 560); 

                auto tg_xyzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 561); 

                auto tg_xyzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 562); 

                auto tg_xyzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 563); 

                auto tg_xyzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 564); 

                auto tg_xyzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 565); 

                auto tg_xyzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 566); 

                auto tg_xzzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 567); 

                auto tg_xzzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 568); 

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

                auto tg_zzzzzz_xxxxy_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 568); 

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

                auto tg_zzzzzz_xxxxy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 568); 

                auto tg_xyyyyyz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 343); 

                auto tg_xyyyyyz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 344); 

                auto tg_xyyyyzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 345); 

                auto tg_xyyyyzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 346); 

                auto tg_xyyyyzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 347); 

                auto tg_xyyyyzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 348); 

                auto tg_xyyyyzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 349); 

                auto tg_xyyyyzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 350); 

                auto tg_xyyyyzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 351); 

                auto tg_xyyyyzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 352); 

                auto tg_xyyyyzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 353); 

                auto tg_xyyyyzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 354); 

                auto tg_xyyyyzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 355); 

                auto tg_xyyyyzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 356); 

                auto tg_xyyyyzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 357); 

                auto tg_xyyyyzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 358); 

                auto tg_xyyyyzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 359); 

                auto tg_xyyyzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 360); 

                auto tg_xyyyzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 361); 

                auto tg_xyyyzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 362); 

                auto tg_xyyyzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 363); 

                auto tg_xyyyzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 364); 

                auto tg_xyyyzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 365); 

                auto tg_xyyyzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 366); 

                auto tg_xyyyzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 367); 

                auto tg_xyyyzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 368); 

                auto tg_xyyyzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 369); 

                auto tg_xyyyzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 370); 

                auto tg_xyyyzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 371); 

                auto tg_xyyyzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 372); 

                auto tg_xyyyzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 373); 

                auto tg_xyyyzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 374); 

                auto tg_xyyzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 375); 

                auto tg_xyyzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 376); 

                auto tg_xyyzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 377); 

                auto tg_xyyzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 378); 

                auto tg_xyyzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 379); 

                auto tg_xyyzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 380); 

                auto tg_xyyzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 381); 

                auto tg_xyyzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 382); 

                auto tg_xyyzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 383); 

                auto tg_xyyzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 384); 

                auto tg_xyyzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 385); 

                auto tg_xyyzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 386); 

                auto tg_xyyzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 387); 

                auto tg_xyyzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 388); 

                auto tg_xyyzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 389); 

                auto tg_xyzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 390); 

                auto tg_xyzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 391); 

                auto tg_xyzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 392); 

                auto tg_xyzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 393); 

                auto tg_xyzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 394); 

                auto tg_xyzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 395); 

                auto tg_xyzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 396); 

                auto tg_xyzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 397); 

                auto tg_xyzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 398); 

                auto tg_xyzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 399); 

                auto tg_xyzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 400); 

                auto tg_xyzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 401); 

                auto tg_xyzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 402); 

                auto tg_xyzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 403); 

                auto tg_xyzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 404); 

                auto tg_xzzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 405); 

                auto tg_xzzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 406); 

                // set up pointers to integrals

                auto tg_xxyyyyyz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 475); 

                auto tg_xxyyyyyz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 476); 

                auto tg_xxyyyyyz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 477); 

                auto tg_xxyyyyyz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 478); 

                auto tg_xxyyyyyz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 479); 

                auto tg_xxyyyyyz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 480); 

                auto tg_xxyyyyyz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 481); 

                auto tg_xxyyyyyz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 482); 

                auto tg_xxyyyyzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 483); 

                auto tg_xxyyyyzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 484); 

                auto tg_xxyyyyzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 485); 

                auto tg_xxyyyyzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 486); 

                auto tg_xxyyyyzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 487); 

                auto tg_xxyyyyzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 488); 

                auto tg_xxyyyyzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 489); 

                auto tg_xxyyyyzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 490); 

                auto tg_xxyyyyzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 491); 

                auto tg_xxyyyyzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 492); 

                auto tg_xxyyyyzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 493); 

                auto tg_xxyyyyzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 494); 

                auto tg_xxyyyyzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 495); 

                auto tg_xxyyyyzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 496); 

                auto tg_xxyyyyzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 497); 

                auto tg_xxyyyyzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 498); 

                auto tg_xxyyyyzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 499); 

                auto tg_xxyyyyzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 500); 

                auto tg_xxyyyyzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 501); 

                auto tg_xxyyyyzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 502); 

                auto tg_xxyyyyzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 503); 

                auto tg_xxyyyzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 504); 

                auto tg_xxyyyzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 505); 

                auto tg_xxyyyzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 506); 

                auto tg_xxyyyzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 507); 

                auto tg_xxyyyzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 508); 

                auto tg_xxyyyzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 509); 

                auto tg_xxyyyzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 510); 

                auto tg_xxyyyzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 511); 

                auto tg_xxyyyzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 512); 

                auto tg_xxyyyzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 513); 

                auto tg_xxyyyzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 514); 

                auto tg_xxyyyzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 515); 

                auto tg_xxyyyzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 516); 

                auto tg_xxyyyzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 517); 

                auto tg_xxyyyzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 518); 

                auto tg_xxyyyzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 519); 

                auto tg_xxyyyzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 520); 

                auto tg_xxyyyzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 521); 

                auto tg_xxyyyzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 522); 

                auto tg_xxyyyzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 523); 

                auto tg_xxyyyzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 524); 

                auto tg_xxyyzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 525); 

                auto tg_xxyyzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 526); 

                auto tg_xxyyzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 527); 

                auto tg_xxyyzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 528); 

                auto tg_xxyyzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 529); 

                auto tg_xxyyzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 530); 

                auto tg_xxyyzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 531); 

                auto tg_xxyyzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 532); 

                auto tg_xxyyzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 533); 

                auto tg_xxyyzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 534); 

                auto tg_xxyyzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 535); 

                auto tg_xxyyzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 536); 

                auto tg_xxyyzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 537); 

                auto tg_xxyyzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 538); 

                auto tg_xxyyzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 539); 

                auto tg_xxyyzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 540); 

                auto tg_xxyyzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 541); 

                auto tg_xxyyzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 542); 

                auto tg_xxyyzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 543); 

                auto tg_xxyyzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 544); 

                auto tg_xxyyzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 545); 

                auto tg_xxyzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 546); 

                auto tg_xxyzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 547); 

                auto tg_xxyzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 548); 

                auto tg_xxyzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 549); 

                auto tg_xxyzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 550); 

                auto tg_xxyzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 551); 

                auto tg_xxyzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 552); 

                auto tg_xxyzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 553); 

                auto tg_xxyzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 554); 

                auto tg_xxyzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 555); 

                auto tg_xxyzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 556); 

                auto tg_xxyzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 557); 

                auto tg_xxyzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 558); 

                auto tg_xxyzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 559); 

                auto tg_xxyzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 560); 

                auto tg_xxyzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 561); 

                auto tg_xxyzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 562); 

                auto tg_xxyzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 563); 

                auto tg_xxyzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 564); 

                auto tg_xxyzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 565); 

                auto tg_xxyzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 566); 

                auto tg_xxzzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 567); 

                auto tg_xxzzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 568); 

                // Batch of Integrals (475,569)

                #pragma omp simd aligned(fxn, fza, tg_xxyyyyyz_xyzzz_0, tg_xxyyyyyz_xzzzz_0, \
                                         tg_xxyyyyyz_yyyyy_0, tg_xxyyyyyz_yyyyz_0, tg_xxyyyyyz_yyyzz_0, tg_xxyyyyyz_yyzzz_0, \
                                         tg_xxyyyyyz_yzzzz_0, tg_xxyyyyyz_zzzzz_0, tg_xxyyyyzz_xxxxx_0, tg_xxyyyyzz_xxxxy_0, \
                                         tg_xxyyyyzz_xxxxz_0, tg_xxyyyyzz_xxxyy_0, tg_xxyyyyzz_xxxyz_0, tg_xxyyyyzz_xxxzz_0, \
                                         tg_xxyyyyzz_xxyyy_0, tg_xxyyyyzz_xxyyz_0, tg_xxyyyyzz_xxyzz_0, tg_xxyyyyzz_xxzzz_0, \
                                         tg_xxyyyyzz_xyyyy_0, tg_xxyyyyzz_xyyyz_0, tg_xxyyyyzz_xyyzz_0, tg_xxyyyyzz_xyzzz_0, \
                                         tg_xxyyyyzz_xzzzz_0, tg_xxyyyyzz_yyyyy_0, tg_xxyyyyzz_yyyyz_0, tg_xxyyyyzz_yyyzz_0, \
                                         tg_xxyyyyzz_yyzzz_0, tg_xxyyyyzz_yzzzz_0, tg_xxyyyyzz_zzzzz_0, tg_xxyyyzzz_xxxxx_0, \
                                         tg_xxyyyzzz_xxxxy_0, tg_xxyyyzzz_xxxxz_0, tg_xxyyyzzz_xxxyy_0, tg_xxyyyzzz_xxxyz_0, \
                                         tg_xxyyyzzz_xxxzz_0, tg_xxyyyzzz_xxyyy_0, tg_xxyyyzzz_xxyyz_0, tg_xxyyyzzz_xxyzz_0, \
                                         tg_xxyyyzzz_xxzzz_0, tg_xxyyyzzz_xyyyy_0, tg_xxyyyzzz_xyyyz_0, tg_xxyyyzzz_xyyzz_0, \
                                         tg_xxyyyzzz_xyzzz_0, tg_xxyyyzzz_xzzzz_0, tg_xxyyyzzz_yyyyy_0, tg_xxyyyzzz_yyyyz_0, \
                                         tg_xxyyyzzz_yyyzz_0, tg_xxyyyzzz_yyzzz_0, tg_xxyyyzzz_yzzzz_0, tg_xxyyyzzz_zzzzz_0, \
                                         tg_xxyyzzzz_xxxxx_0, tg_xxyyzzzz_xxxxy_0, tg_xxyyzzzz_xxxxz_0, tg_xxyyzzzz_xxxyy_0, \
                                         tg_xxyyzzzz_xxxyz_0, tg_xxyyzzzz_xxxzz_0, tg_xxyyzzzz_xxyyy_0, tg_xxyyzzzz_xxyyz_0, \
                                         tg_xxyyzzzz_xxyzz_0, tg_xxyyzzzz_xxzzz_0, tg_xxyyzzzz_xyyyy_0, tg_xxyyzzzz_xyyyz_0, \
                                         tg_xxyyzzzz_xyyzz_0, tg_xxyyzzzz_xyzzz_0, tg_xxyyzzzz_xzzzz_0, tg_xxyyzzzz_yyyyy_0, \
                                         tg_xxyyzzzz_yyyyz_0, tg_xxyyzzzz_yyyzz_0, tg_xxyyzzzz_yyzzz_0, tg_xxyyzzzz_yzzzz_0, \
                                         tg_xxyyzzzz_zzzzz_0, tg_xxyzzzzz_xxxxx_0, tg_xxyzzzzz_xxxxy_0, tg_xxyzzzzz_xxxxz_0, \
                                         tg_xxyzzzzz_xxxyy_0, tg_xxyzzzzz_xxxyz_0, tg_xxyzzzzz_xxxzz_0, tg_xxyzzzzz_xxyyy_0, \
                                         tg_xxyzzzzz_xxyyz_0, tg_xxyzzzzz_xxyzz_0, tg_xxyzzzzz_xxzzz_0, tg_xxyzzzzz_xyyyy_0, \
                                         tg_xxyzzzzz_xyyyz_0, tg_xxyzzzzz_xyyzz_0, tg_xxyzzzzz_xyzzz_0, tg_xxyzzzzz_xzzzz_0, \
                                         tg_xxyzzzzz_yyyyy_0, tg_xxyzzzzz_yyyyz_0, tg_xxyzzzzz_yyyzz_0, tg_xxyzzzzz_yyzzz_0, \
                                         tg_xxyzzzzz_yzzzz_0, tg_xxyzzzzz_zzzzz_0, tg_xxzzzzzz_xxxxx_0, tg_xxzzzzzz_xxxxy_0, \
                                         tg_xyyyyyz_xyzzz_0, tg_xyyyyyz_xyzzz_1, tg_xyyyyyz_xzzzz_0, tg_xyyyyyz_xzzzz_1, \
                                         tg_xyyyyyz_yyyyy_0, tg_xyyyyyz_yyyyy_1, tg_xyyyyyz_yyyyz_0, tg_xyyyyyz_yyyyz_1, \
                                         tg_xyyyyyz_yyyzz_0, tg_xyyyyyz_yyyzz_1, tg_xyyyyyz_yyzzz_0, tg_xyyyyyz_yyzzz_1, \
                                         tg_xyyyyyz_yzzz_1, tg_xyyyyyz_yzzzz_0, tg_xyyyyyz_yzzzz_1, tg_xyyyyyz_zzzz_1, \
                                         tg_xyyyyyz_zzzzz_0, tg_xyyyyyz_zzzzz_1, tg_xyyyyzz_xxxx_1, tg_xyyyyzz_xxxxx_0, \
                                         tg_xyyyyzz_xxxxx_1, tg_xyyyyzz_xxxxy_0, tg_xyyyyzz_xxxxy_1, tg_xyyyyzz_xxxxz_0, \
                                         tg_xyyyyzz_xxxxz_1, tg_xyyyyzz_xxxy_1, tg_xyyyyzz_xxxyy_0, tg_xyyyyzz_xxxyy_1, \
                                         tg_xyyyyzz_xxxyz_0, tg_xyyyyzz_xxxyz_1, tg_xyyyyzz_xxxz_1, tg_xyyyyzz_xxxzz_0, \
                                         tg_xyyyyzz_xxxzz_1, tg_xyyyyzz_xxyy_1, tg_xyyyyzz_xxyyy_0, tg_xyyyyzz_xxyyy_1, \
                                         tg_xyyyyzz_xxyyz_0, tg_xyyyyzz_xxyyz_1, tg_xyyyyzz_xxyz_1, tg_xyyyyzz_xxyzz_0, \
                                         tg_xyyyyzz_xxyzz_1, tg_xyyyyzz_xxzz_1, tg_xyyyyzz_xxzzz_0, tg_xyyyyzz_xxzzz_1, \
                                         tg_xyyyyzz_xyyy_1, tg_xyyyyzz_xyyyy_0, tg_xyyyyzz_xyyyy_1, tg_xyyyyzz_xyyyz_0, \
                                         tg_xyyyyzz_xyyyz_1, tg_xyyyyzz_xyyz_1, tg_xyyyyzz_xyyzz_0, tg_xyyyyzz_xyyzz_1, \
                                         tg_xyyyyzz_xyzz_1, tg_xyyyyzz_xyzzz_0, tg_xyyyyzz_xyzzz_1, tg_xyyyyzz_xzzz_1, \
                                         tg_xyyyyzz_xzzzz_0, tg_xyyyyzz_xzzzz_1, tg_xyyyyzz_yyyy_1, tg_xyyyyzz_yyyyy_0, \
                                         tg_xyyyyzz_yyyyy_1, tg_xyyyyzz_yyyyz_0, tg_xyyyyzz_yyyyz_1, tg_xyyyyzz_yyyz_1, \
                                         tg_xyyyyzz_yyyzz_0, tg_xyyyyzz_yyyzz_1, tg_xyyyyzz_yyzz_1, tg_xyyyyzz_yyzzz_0, \
                                         tg_xyyyyzz_yyzzz_1, tg_xyyyyzz_yzzz_1, tg_xyyyyzz_yzzzz_0, tg_xyyyyzz_yzzzz_1, \
                                         tg_xyyyyzz_zzzz_1, tg_xyyyyzz_zzzzz_0, tg_xyyyyzz_zzzzz_1, tg_xyyyzzz_xxxx_1, \
                                         tg_xyyyzzz_xxxxx_0, tg_xyyyzzz_xxxxx_1, tg_xyyyzzz_xxxxy_0, tg_xyyyzzz_xxxxy_1, \
                                         tg_xyyyzzz_xxxxz_0, tg_xyyyzzz_xxxxz_1, tg_xyyyzzz_xxxy_1, tg_xyyyzzz_xxxyy_0, \
                                         tg_xyyyzzz_xxxyy_1, tg_xyyyzzz_xxxyz_0, tg_xyyyzzz_xxxyz_1, tg_xyyyzzz_xxxz_1, \
                                         tg_xyyyzzz_xxxzz_0, tg_xyyyzzz_xxxzz_1, tg_xyyyzzz_xxyy_1, tg_xyyyzzz_xxyyy_0, \
                                         tg_xyyyzzz_xxyyy_1, tg_xyyyzzz_xxyyz_0, tg_xyyyzzz_xxyyz_1, tg_xyyyzzz_xxyz_1, \
                                         tg_xyyyzzz_xxyzz_0, tg_xyyyzzz_xxyzz_1, tg_xyyyzzz_xxzz_1, tg_xyyyzzz_xxzzz_0, \
                                         tg_xyyyzzz_xxzzz_1, tg_xyyyzzz_xyyy_1, tg_xyyyzzz_xyyyy_0, tg_xyyyzzz_xyyyy_1, \
                                         tg_xyyyzzz_xyyyz_0, tg_xyyyzzz_xyyyz_1, tg_xyyyzzz_xyyz_1, tg_xyyyzzz_xyyzz_0, \
                                         tg_xyyyzzz_xyyzz_1, tg_xyyyzzz_xyzz_1, tg_xyyyzzz_xyzzz_0, tg_xyyyzzz_xyzzz_1, \
                                         tg_xyyyzzz_xzzz_1, tg_xyyyzzz_xzzzz_0, tg_xyyyzzz_xzzzz_1, tg_xyyyzzz_yyyy_1, \
                                         tg_xyyyzzz_yyyyy_0, tg_xyyyzzz_yyyyy_1, tg_xyyyzzz_yyyyz_0, tg_xyyyzzz_yyyyz_1, \
                                         tg_xyyyzzz_yyyz_1, tg_xyyyzzz_yyyzz_0, tg_xyyyzzz_yyyzz_1, tg_xyyyzzz_yyzz_1, \
                                         tg_xyyyzzz_yyzzz_0, tg_xyyyzzz_yyzzz_1, tg_xyyyzzz_yzzz_1, tg_xyyyzzz_yzzzz_0, \
                                         tg_xyyyzzz_yzzzz_1, tg_xyyyzzz_zzzz_1, tg_xyyyzzz_zzzzz_0, tg_xyyyzzz_zzzzz_1, \
                                         tg_xyyzzzz_xxxx_1, tg_xyyzzzz_xxxxx_0, tg_xyyzzzz_xxxxx_1, tg_xyyzzzz_xxxxy_0, \
                                         tg_xyyzzzz_xxxxy_1, tg_xyyzzzz_xxxxz_0, tg_xyyzzzz_xxxxz_1, tg_xyyzzzz_xxxy_1, \
                                         tg_xyyzzzz_xxxyy_0, tg_xyyzzzz_xxxyy_1, tg_xyyzzzz_xxxyz_0, tg_xyyzzzz_xxxyz_1, \
                                         tg_xyyzzzz_xxxz_1, tg_xyyzzzz_xxxzz_0, tg_xyyzzzz_xxxzz_1, tg_xyyzzzz_xxyy_1, \
                                         tg_xyyzzzz_xxyyy_0, tg_xyyzzzz_xxyyy_1, tg_xyyzzzz_xxyyz_0, tg_xyyzzzz_xxyyz_1, \
                                         tg_xyyzzzz_xxyz_1, tg_xyyzzzz_xxyzz_0, tg_xyyzzzz_xxyzz_1, tg_xyyzzzz_xxzz_1, \
                                         tg_xyyzzzz_xxzzz_0, tg_xyyzzzz_xxzzz_1, tg_xyyzzzz_xyyy_1, tg_xyyzzzz_xyyyy_0, \
                                         tg_xyyzzzz_xyyyy_1, tg_xyyzzzz_xyyyz_0, tg_xyyzzzz_xyyyz_1, tg_xyyzzzz_xyyz_1, \
                                         tg_xyyzzzz_xyyzz_0, tg_xyyzzzz_xyyzz_1, tg_xyyzzzz_xyzz_1, tg_xyyzzzz_xyzzz_0, \
                                         tg_xyyzzzz_xyzzz_1, tg_xyyzzzz_xzzz_1, tg_xyyzzzz_xzzzz_0, tg_xyyzzzz_xzzzz_1, \
                                         tg_xyyzzzz_yyyy_1, tg_xyyzzzz_yyyyy_0, tg_xyyzzzz_yyyyy_1, tg_xyyzzzz_yyyyz_0, \
                                         tg_xyyzzzz_yyyyz_1, tg_xyyzzzz_yyyz_1, tg_xyyzzzz_yyyzz_0, tg_xyyzzzz_yyyzz_1, \
                                         tg_xyyzzzz_yyzz_1, tg_xyyzzzz_yyzzz_0, tg_xyyzzzz_yyzzz_1, tg_xyyzzzz_yzzz_1, \
                                         tg_xyyzzzz_yzzzz_0, tg_xyyzzzz_yzzzz_1, tg_xyyzzzz_zzzz_1, tg_xyyzzzz_zzzzz_0, \
                                         tg_xyyzzzz_zzzzz_1, tg_xyzzzzz_xxxx_1, tg_xyzzzzz_xxxxx_0, tg_xyzzzzz_xxxxx_1, \
                                         tg_xyzzzzz_xxxxy_0, tg_xyzzzzz_xxxxy_1, tg_xyzzzzz_xxxxz_0, tg_xyzzzzz_xxxxz_1, \
                                         tg_xyzzzzz_xxxy_1, tg_xyzzzzz_xxxyy_0, tg_xyzzzzz_xxxyy_1, tg_xyzzzzz_xxxyz_0, \
                                         tg_xyzzzzz_xxxyz_1, tg_xyzzzzz_xxxz_1, tg_xyzzzzz_xxxzz_0, tg_xyzzzzz_xxxzz_1, \
                                         tg_xyzzzzz_xxyy_1, tg_xyzzzzz_xxyyy_0, tg_xyzzzzz_xxyyy_1, tg_xyzzzzz_xxyyz_0, \
                                         tg_xyzzzzz_xxyyz_1, tg_xyzzzzz_xxyz_1, tg_xyzzzzz_xxyzz_0, tg_xyzzzzz_xxyzz_1, \
                                         tg_xyzzzzz_xxzz_1, tg_xyzzzzz_xxzzz_0, tg_xyzzzzz_xxzzz_1, tg_xyzzzzz_xyyy_1, \
                                         tg_xyzzzzz_xyyyy_0, tg_xyzzzzz_xyyyy_1, tg_xyzzzzz_xyyyz_0, tg_xyzzzzz_xyyyz_1, \
                                         tg_xyzzzzz_xyyz_1, tg_xyzzzzz_xyyzz_0, tg_xyzzzzz_xyyzz_1, tg_xyzzzzz_xyzz_1, \
                                         tg_xyzzzzz_xyzzz_0, tg_xyzzzzz_xyzzz_1, tg_xyzzzzz_xzzz_1, tg_xyzzzzz_xzzzz_0, \
                                         tg_xyzzzzz_xzzzz_1, tg_xyzzzzz_yyyy_1, tg_xyzzzzz_yyyyy_0, tg_xyzzzzz_yyyyy_1, \
                                         tg_xyzzzzz_yyyyz_0, tg_xyzzzzz_yyyyz_1, tg_xyzzzzz_yyyz_1, tg_xyzzzzz_yyyzz_0, \
                                         tg_xyzzzzz_yyyzz_1, tg_xyzzzzz_yyzz_1, tg_xyzzzzz_yyzzz_0, tg_xyzzzzz_yyzzz_1, \
                                         tg_xyzzzzz_yzzz_1, tg_xyzzzzz_yzzzz_0, tg_xyzzzzz_yzzzz_1, tg_xyzzzzz_zzzz_1, \
                                         tg_xyzzzzz_zzzzz_0, tg_xyzzzzz_zzzzz_1, tg_xzzzzzz_xxxx_1, tg_xzzzzzz_xxxxx_0, \
                                         tg_xzzzzzz_xxxxx_1, tg_xzzzzzz_xxxxy_0, tg_xzzzzzz_xxxxy_1, tg_xzzzzzz_xxxy_1, \
                                         tg_yyyyyz_xyzzz_0, tg_yyyyyz_xyzzz_1, tg_yyyyyz_xzzzz_0, tg_yyyyyz_xzzzz_1, \
                                         tg_yyyyyz_yyyyy_0, tg_yyyyyz_yyyyy_1, tg_yyyyyz_yyyyz_0, tg_yyyyyz_yyyyz_1, \
                                         tg_yyyyyz_yyyzz_0, tg_yyyyyz_yyyzz_1, tg_yyyyyz_yyzzz_0, tg_yyyyyz_yyzzz_1, \
                                         tg_yyyyyz_yzzzz_0, tg_yyyyyz_yzzzz_1, tg_yyyyyz_zzzzz_0, tg_yyyyyz_zzzzz_1, \
                                         tg_yyyyzz_xxxxx_0, tg_yyyyzz_xxxxx_1, tg_yyyyzz_xxxxy_0, tg_yyyyzz_xxxxy_1, \
                                         tg_yyyyzz_xxxxz_0, tg_yyyyzz_xxxxz_1, tg_yyyyzz_xxxyy_0, tg_yyyyzz_xxxyy_1, \
                                         tg_yyyyzz_xxxyz_0, tg_yyyyzz_xxxyz_1, tg_yyyyzz_xxxzz_0, tg_yyyyzz_xxxzz_1, \
                                         tg_yyyyzz_xxyyy_0, tg_yyyyzz_xxyyy_1, tg_yyyyzz_xxyyz_0, tg_yyyyzz_xxyyz_1, \
                                         tg_yyyyzz_xxyzz_0, tg_yyyyzz_xxyzz_1, tg_yyyyzz_xxzzz_0, tg_yyyyzz_xxzzz_1, \
                                         tg_yyyyzz_xyyyy_0, tg_yyyyzz_xyyyy_1, tg_yyyyzz_xyyyz_0, tg_yyyyzz_xyyyz_1, \
                                         tg_yyyyzz_xyyzz_0, tg_yyyyzz_xyyzz_1, tg_yyyyzz_xyzzz_0, tg_yyyyzz_xyzzz_1, \
                                         tg_yyyyzz_xzzzz_0, tg_yyyyzz_xzzzz_1, tg_yyyyzz_yyyyy_0, tg_yyyyzz_yyyyy_1, \
                                         tg_yyyyzz_yyyyz_0, tg_yyyyzz_yyyyz_1, tg_yyyyzz_yyyzz_0, tg_yyyyzz_yyyzz_1, \
                                         tg_yyyyzz_yyzzz_0, tg_yyyyzz_yyzzz_1, tg_yyyyzz_yzzzz_0, tg_yyyyzz_yzzzz_1, \
                                         tg_yyyyzz_zzzzz_0, tg_yyyyzz_zzzzz_1, tg_yyyzzz_xxxxx_0, tg_yyyzzz_xxxxx_1, \
                                         tg_yyyzzz_xxxxy_0, tg_yyyzzz_xxxxy_1, tg_yyyzzz_xxxxz_0, tg_yyyzzz_xxxxz_1, \
                                         tg_yyyzzz_xxxyy_0, tg_yyyzzz_xxxyy_1, tg_yyyzzz_xxxyz_0, tg_yyyzzz_xxxyz_1, \
                                         tg_yyyzzz_xxxzz_0, tg_yyyzzz_xxxzz_1, tg_yyyzzz_xxyyy_0, tg_yyyzzz_xxyyy_1, \
                                         tg_yyyzzz_xxyyz_0, tg_yyyzzz_xxyyz_1, tg_yyyzzz_xxyzz_0, tg_yyyzzz_xxyzz_1, \
                                         tg_yyyzzz_xxzzz_0, tg_yyyzzz_xxzzz_1, tg_yyyzzz_xyyyy_0, tg_yyyzzz_xyyyy_1, \
                                         tg_yyyzzz_xyyyz_0, tg_yyyzzz_xyyyz_1, tg_yyyzzz_xyyzz_0, tg_yyyzzz_xyyzz_1, \
                                         tg_yyyzzz_xyzzz_0, tg_yyyzzz_xyzzz_1, tg_yyyzzz_xzzzz_0, tg_yyyzzz_xzzzz_1, \
                                         tg_yyyzzz_yyyyy_0, tg_yyyzzz_yyyyy_1, tg_yyyzzz_yyyyz_0, tg_yyyzzz_yyyyz_1, \
                                         tg_yyyzzz_yyyzz_0, tg_yyyzzz_yyyzz_1, tg_yyyzzz_yyzzz_0, tg_yyyzzz_yyzzz_1, \
                                         tg_yyyzzz_yzzzz_0, tg_yyyzzz_yzzzz_1, tg_yyyzzz_zzzzz_0, tg_yyyzzz_zzzzz_1, \
                                         tg_yyzzzz_xxxxx_0, tg_yyzzzz_xxxxx_1, tg_yyzzzz_xxxxy_0, tg_yyzzzz_xxxxy_1, \
                                         tg_yyzzzz_xxxxz_0, tg_yyzzzz_xxxxz_1, tg_yyzzzz_xxxyy_0, tg_yyzzzz_xxxyy_1, \
                                         tg_yyzzzz_xxxyz_0, tg_yyzzzz_xxxyz_1, tg_yyzzzz_xxxzz_0, tg_yyzzzz_xxxzz_1, \
                                         tg_yyzzzz_xxyyy_0, tg_yyzzzz_xxyyy_1, tg_yyzzzz_xxyyz_0, tg_yyzzzz_xxyyz_1, \
                                         tg_yyzzzz_xxyzz_0, tg_yyzzzz_xxyzz_1, tg_yyzzzz_xxzzz_0, tg_yyzzzz_xxzzz_1, \
                                         tg_yyzzzz_xyyyy_0, tg_yyzzzz_xyyyy_1, tg_yyzzzz_xyyyz_0, tg_yyzzzz_xyyyz_1, \
                                         tg_yyzzzz_xyyzz_0, tg_yyzzzz_xyyzz_1, tg_yyzzzz_xyzzz_0, tg_yyzzzz_xyzzz_1, \
                                         tg_yyzzzz_xzzzz_0, tg_yyzzzz_xzzzz_1, tg_yyzzzz_yyyyy_0, tg_yyzzzz_yyyyy_1, \
                                         tg_yyzzzz_yyyyz_0, tg_yyzzzz_yyyyz_1, tg_yyzzzz_yyyzz_0, tg_yyzzzz_yyyzz_1, \
                                         tg_yyzzzz_yyzzz_0, tg_yyzzzz_yyzzz_1, tg_yyzzzz_yzzzz_0, tg_yyzzzz_yzzzz_1, \
                                         tg_yyzzzz_zzzzz_0, tg_yyzzzz_zzzzz_1, tg_yzzzzz_xxxxx_0, tg_yzzzzz_xxxxx_1, \
                                         tg_yzzzzz_xxxxy_0, tg_yzzzzz_xxxxy_1, tg_yzzzzz_xxxxz_0, tg_yzzzzz_xxxxz_1, \
                                         tg_yzzzzz_xxxyy_0, tg_yzzzzz_xxxyy_1, tg_yzzzzz_xxxyz_0, tg_yzzzzz_xxxyz_1, \
                                         tg_yzzzzz_xxxzz_0, tg_yzzzzz_xxxzz_1, tg_yzzzzz_xxyyy_0, tg_yzzzzz_xxyyy_1, \
                                         tg_yzzzzz_xxyyz_0, tg_yzzzzz_xxyyz_1, tg_yzzzzz_xxyzz_0, tg_yzzzzz_xxyzz_1, \
                                         tg_yzzzzz_xxzzz_0, tg_yzzzzz_xxzzz_1, tg_yzzzzz_xyyyy_0, tg_yzzzzz_xyyyy_1, \
                                         tg_yzzzzz_xyyyz_0, tg_yzzzzz_xyyyz_1, tg_yzzzzz_xyyzz_0, tg_yzzzzz_xyyzz_1, \
                                         tg_yzzzzz_xyzzz_0, tg_yzzzzz_xyzzz_1, tg_yzzzzz_xzzzz_0, tg_yzzzzz_xzzzz_1, \
                                         tg_yzzzzz_yyyyy_0, tg_yzzzzz_yyyyy_1, tg_yzzzzz_yyyyz_0, tg_yzzzzz_yyyyz_1, \
                                         tg_yzzzzz_yyyzz_0, tg_yzzzzz_yyyzz_1, tg_yzzzzz_yyzzz_0, tg_yzzzzz_yyzzz_1, \
                                         tg_yzzzzz_yzzzz_0, tg_yzzzzz_yzzzz_1, tg_yzzzzz_zzzzz_0, tg_yzzzzz_zzzzz_1, \
                                         tg_zzzzzz_xxxxx_0, tg_zzzzzz_xxxxx_1, tg_zzzzzz_xxxxy_0, tg_zzzzzz_xxxxy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyyyyz_xyzzz_0[j] = pb_x * tg_xyyyyyz_xyzzz_0[j] + fr * tg_xyyyyyz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xyzzz_0[j] - tg_yyyyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_yzzz_1[j];

                    tg_xxyyyyyz_xzzzz_0[j] = pb_x * tg_xyyyyyz_xzzzz_0[j] + fr * tg_xyyyyyz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xzzzz_0[j] - tg_yyyyyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_zzzz_1[j];

                    tg_xxyyyyyz_yyyyy_0[j] = pb_x * tg_xyyyyyz_yyyyy_0[j] + fr * tg_xyyyyyz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yyyyy_0[j] - tg_yyyyyz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyyyz_yyyyz_0[j] = pb_x * tg_xyyyyyz_yyyyz_0[j] + fr * tg_xyyyyyz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yyyyz_0[j] - tg_yyyyyz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyyyz_yyyzz_0[j] = pb_x * tg_xyyyyyz_yyyzz_0[j] + fr * tg_xyyyyyz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yyyzz_0[j] - tg_yyyyyz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyyyz_yyzzz_0[j] = pb_x * tg_xyyyyyz_yyzzz_0[j] + fr * tg_xyyyyyz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yyzzz_0[j] - tg_yyyyyz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyyyz_yzzzz_0[j] = pb_x * tg_xyyyyyz_yzzzz_0[j] + fr * tg_xyyyyyz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yzzzz_0[j] - tg_yyyyyz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyyyz_zzzzz_0[j] = pb_x * tg_xyyyyyz_zzzzz_0[j] + fr * tg_xyyyyyz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_zzzzz_0[j] - tg_yyyyyz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyyzz_xxxxx_0[j] = pb_x * tg_xyyyyzz_xxxxx_0[j] + fr * tg_xyyyyzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxxxx_0[j] - tg_yyyyzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyzz_xxxx_1[j];

                    tg_xxyyyyzz_xxxxy_0[j] = pb_x * tg_xyyyyzz_xxxxy_0[j] + fr * tg_xyyyyzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxxxy_0[j] - tg_yyyyzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyzz_xxxy_1[j];

                    tg_xxyyyyzz_xxxxz_0[j] = pb_x * tg_xyyyyzz_xxxxz_0[j] + fr * tg_xyyyyzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxxxz_0[j] - tg_yyyyzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyzz_xxxz_1[j];

                    tg_xxyyyyzz_xxxyy_0[j] = pb_x * tg_xyyyyzz_xxxyy_0[j] + fr * tg_xyyyyzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxxyy_0[j] - tg_yyyyzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyzz_xxyy_1[j];

                    tg_xxyyyyzz_xxxyz_0[j] = pb_x * tg_xyyyyzz_xxxyz_0[j] + fr * tg_xyyyyzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxxyz_0[j] - tg_yyyyzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyzz_xxyz_1[j];

                    tg_xxyyyyzz_xxxzz_0[j] = pb_x * tg_xyyyyzz_xxxzz_0[j] + fr * tg_xyyyyzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxxzz_0[j] - tg_yyyyzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyzz_xxzz_1[j];

                    tg_xxyyyyzz_xxyyy_0[j] = pb_x * tg_xyyyyzz_xxyyy_0[j] + fr * tg_xyyyyzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxyyy_0[j] - tg_yyyyzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyzz_xyyy_1[j];

                    tg_xxyyyyzz_xxyyz_0[j] = pb_x * tg_xyyyyzz_xxyyz_0[j] + fr * tg_xyyyyzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxyyz_0[j] - tg_yyyyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyzz_xyyz_1[j];

                    tg_xxyyyyzz_xxyzz_0[j] = pb_x * tg_xyyyyzz_xxyzz_0[j] + fr * tg_xyyyyzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxyzz_0[j] - tg_yyyyzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyzz_xyzz_1[j];

                    tg_xxyyyyzz_xxzzz_0[j] = pb_x * tg_xyyyyzz_xxzzz_0[j] + fr * tg_xyyyyzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxzzz_0[j] - tg_yyyyzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyzz_xzzz_1[j];

                    tg_xxyyyyzz_xyyyy_0[j] = pb_x * tg_xyyyyzz_xyyyy_0[j] + fr * tg_xyyyyzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xyyyy_0[j] - tg_yyyyzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_yyyy_1[j];

                    tg_xxyyyyzz_xyyyz_0[j] = pb_x * tg_xyyyyzz_xyyyz_0[j] + fr * tg_xyyyyzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xyyyz_0[j] - tg_yyyyzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_yyyz_1[j];

                    tg_xxyyyyzz_xyyzz_0[j] = pb_x * tg_xyyyyzz_xyyzz_0[j] + fr * tg_xyyyyzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xyyzz_0[j] - tg_yyyyzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_yyzz_1[j];

                    tg_xxyyyyzz_xyzzz_0[j] = pb_x * tg_xyyyyzz_xyzzz_0[j] + fr * tg_xyyyyzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xyzzz_0[j] - tg_yyyyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_yzzz_1[j];

                    tg_xxyyyyzz_xzzzz_0[j] = pb_x * tg_xyyyyzz_xzzzz_0[j] + fr * tg_xyyyyzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xzzzz_0[j] - tg_yyyyzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_zzzz_1[j];

                    tg_xxyyyyzz_yyyyy_0[j] = pb_x * tg_xyyyyzz_yyyyy_0[j] + fr * tg_xyyyyzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yyyyy_0[j] - tg_yyyyzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyyzz_yyyyz_0[j] = pb_x * tg_xyyyyzz_yyyyz_0[j] + fr * tg_xyyyyzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yyyyz_0[j] - tg_yyyyzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyyzz_yyyzz_0[j] = pb_x * tg_xyyyyzz_yyyzz_0[j] + fr * tg_xyyyyzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yyyzz_0[j] - tg_yyyyzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyyzz_yyzzz_0[j] = pb_x * tg_xyyyyzz_yyzzz_0[j] + fr * tg_xyyyyzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yyzzz_0[j] - tg_yyyyzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyyzz_yzzzz_0[j] = pb_x * tg_xyyyyzz_yzzzz_0[j] + fr * tg_xyyyyzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yzzzz_0[j] - tg_yyyyzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyyzz_zzzzz_0[j] = pb_x * tg_xyyyyzz_zzzzz_0[j] + fr * tg_xyyyyzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_zzzzz_0[j] - tg_yyyyzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyyzzz_xxxxx_0[j] = pb_x * tg_xyyyzzz_xxxxx_0[j] + fr * tg_xyyyzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxxxx_0[j] - tg_yyyzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyzzz_xxxx_1[j];

                    tg_xxyyyzzz_xxxxy_0[j] = pb_x * tg_xyyyzzz_xxxxy_0[j] + fr * tg_xyyyzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxxxy_0[j] - tg_yyyzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzzz_xxxy_1[j];

                    tg_xxyyyzzz_xxxxz_0[j] = pb_x * tg_xyyyzzz_xxxxz_0[j] + fr * tg_xyyyzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxxxz_0[j] - tg_yyyzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzzz_xxxz_1[j];

                    tg_xxyyyzzz_xxxyy_0[j] = pb_x * tg_xyyyzzz_xxxyy_0[j] + fr * tg_xyyyzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxxyy_0[j] - tg_yyyzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzzz_xxyy_1[j];

                    tg_xxyyyzzz_xxxyz_0[j] = pb_x * tg_xyyyzzz_xxxyz_0[j] + fr * tg_xyyyzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxxyz_0[j] - tg_yyyzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzzz_xxyz_1[j];

                    tg_xxyyyzzz_xxxzz_0[j] = pb_x * tg_xyyyzzz_xxxzz_0[j] + fr * tg_xyyyzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxxzz_0[j] - tg_yyyzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzzz_xxzz_1[j];

                    tg_xxyyyzzz_xxyyy_0[j] = pb_x * tg_xyyyzzz_xxyyy_0[j] + fr * tg_xyyyzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxyyy_0[j] - tg_yyyzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzzz_xyyy_1[j];

                    tg_xxyyyzzz_xxyyz_0[j] = pb_x * tg_xyyyzzz_xxyyz_0[j] + fr * tg_xyyyzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxyyz_0[j] - tg_yyyzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzzz_xyyz_1[j];

                    tg_xxyyyzzz_xxyzz_0[j] = pb_x * tg_xyyyzzz_xxyzz_0[j] + fr * tg_xyyyzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxyzz_0[j] - tg_yyyzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzzz_xyzz_1[j];

                    tg_xxyyyzzz_xxzzz_0[j] = pb_x * tg_xyyyzzz_xxzzz_0[j] + fr * tg_xyyyzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxzzz_0[j] - tg_yyyzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzzz_xzzz_1[j];

                    tg_xxyyyzzz_xyyyy_0[j] = pb_x * tg_xyyyzzz_xyyyy_0[j] + fr * tg_xyyyzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xyyyy_0[j] - tg_yyyzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_yyyy_1[j];

                    tg_xxyyyzzz_xyyyz_0[j] = pb_x * tg_xyyyzzz_xyyyz_0[j] + fr * tg_xyyyzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xyyyz_0[j] - tg_yyyzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_yyyz_1[j];

                    tg_xxyyyzzz_xyyzz_0[j] = pb_x * tg_xyyyzzz_xyyzz_0[j] + fr * tg_xyyyzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xyyzz_0[j] - tg_yyyzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_yyzz_1[j];

                    tg_xxyyyzzz_xyzzz_0[j] = pb_x * tg_xyyyzzz_xyzzz_0[j] + fr * tg_xyyyzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xyzzz_0[j] - tg_yyyzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_yzzz_1[j];

                    tg_xxyyyzzz_xzzzz_0[j] = pb_x * tg_xyyyzzz_xzzzz_0[j] + fr * tg_xyyyzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xzzzz_0[j] - tg_yyyzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_zzzz_1[j];

                    tg_xxyyyzzz_yyyyy_0[j] = pb_x * tg_xyyyzzz_yyyyy_0[j] + fr * tg_xyyyzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yyyyy_0[j] - tg_yyyzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyyzzz_yyyyz_0[j] = pb_x * tg_xyyyzzz_yyyyz_0[j] + fr * tg_xyyyzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yyyyz_0[j] - tg_yyyzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyyzzz_yyyzz_0[j] = pb_x * tg_xyyyzzz_yyyzz_0[j] + fr * tg_xyyyzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yyyzz_0[j] - tg_yyyzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyyzzz_yyzzz_0[j] = pb_x * tg_xyyyzzz_yyzzz_0[j] + fr * tg_xyyyzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yyzzz_0[j] - tg_yyyzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyyzzz_yzzzz_0[j] = pb_x * tg_xyyyzzz_yzzzz_0[j] + fr * tg_xyyyzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yzzzz_0[j] - tg_yyyzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyyzzz_zzzzz_0[j] = pb_x * tg_xyyyzzz_zzzzz_0[j] + fr * tg_xyyyzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_zzzzz_0[j] - tg_yyyzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyzzzz_xxxxx_0[j] = pb_x * tg_xyyzzzz_xxxxx_0[j] + fr * tg_xyyzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxxxx_0[j] - tg_yyzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyzzzz_xxxx_1[j];

                    tg_xxyyzzzz_xxxxy_0[j] = pb_x * tg_xyyzzzz_xxxxy_0[j] + fr * tg_xyyzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxxxy_0[j] - tg_yyzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzzz_xxxy_1[j];

                    tg_xxyyzzzz_xxxxz_0[j] = pb_x * tg_xyyzzzz_xxxxz_0[j] + fr * tg_xyyzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxxxz_0[j] - tg_yyzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzzz_xxxz_1[j];

                    tg_xxyyzzzz_xxxyy_0[j] = pb_x * tg_xyyzzzz_xxxyy_0[j] + fr * tg_xyyzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxxyy_0[j] - tg_yyzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzzz_xxyy_1[j];

                    tg_xxyyzzzz_xxxyz_0[j] = pb_x * tg_xyyzzzz_xxxyz_0[j] + fr * tg_xyyzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxxyz_0[j] - tg_yyzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzzz_xxyz_1[j];

                    tg_xxyyzzzz_xxxzz_0[j] = pb_x * tg_xyyzzzz_xxxzz_0[j] + fr * tg_xyyzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxxzz_0[j] - tg_yyzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzzz_xxzz_1[j];

                    tg_xxyyzzzz_xxyyy_0[j] = pb_x * tg_xyyzzzz_xxyyy_0[j] + fr * tg_xyyzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxyyy_0[j] - tg_yyzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzzz_xyyy_1[j];

                    tg_xxyyzzzz_xxyyz_0[j] = pb_x * tg_xyyzzzz_xxyyz_0[j] + fr * tg_xyyzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxyyz_0[j] - tg_yyzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzzz_xyyz_1[j];

                    tg_xxyyzzzz_xxyzz_0[j] = pb_x * tg_xyyzzzz_xxyzz_0[j] + fr * tg_xyyzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxyzz_0[j] - tg_yyzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzzz_xyzz_1[j];

                    tg_xxyyzzzz_xxzzz_0[j] = pb_x * tg_xyyzzzz_xxzzz_0[j] + fr * tg_xyyzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxzzz_0[j] - tg_yyzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzzz_xzzz_1[j];

                    tg_xxyyzzzz_xyyyy_0[j] = pb_x * tg_xyyzzzz_xyyyy_0[j] + fr * tg_xyyzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xyyyy_0[j] - tg_yyzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_yyyy_1[j];

                    tg_xxyyzzzz_xyyyz_0[j] = pb_x * tg_xyyzzzz_xyyyz_0[j] + fr * tg_xyyzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xyyyz_0[j] - tg_yyzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_yyyz_1[j];

                    tg_xxyyzzzz_xyyzz_0[j] = pb_x * tg_xyyzzzz_xyyzz_0[j] + fr * tg_xyyzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xyyzz_0[j] - tg_yyzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_yyzz_1[j];

                    tg_xxyyzzzz_xyzzz_0[j] = pb_x * tg_xyyzzzz_xyzzz_0[j] + fr * tg_xyyzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xyzzz_0[j] - tg_yyzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_yzzz_1[j];

                    tg_xxyyzzzz_xzzzz_0[j] = pb_x * tg_xyyzzzz_xzzzz_0[j] + fr * tg_xyyzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xzzzz_0[j] - tg_yyzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_zzzz_1[j];

                    tg_xxyyzzzz_yyyyy_0[j] = pb_x * tg_xyyzzzz_yyyyy_0[j] + fr * tg_xyyzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yyyyy_0[j] - tg_yyzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyzzzz_yyyyz_0[j] = pb_x * tg_xyyzzzz_yyyyz_0[j] + fr * tg_xyyzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yyyyz_0[j] - tg_yyzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyzzzz_yyyzz_0[j] = pb_x * tg_xyyzzzz_yyyzz_0[j] + fr * tg_xyyzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yyyzz_0[j] - tg_yyzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyzzzz_yyzzz_0[j] = pb_x * tg_xyyzzzz_yyzzz_0[j] + fr * tg_xyyzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yyzzz_0[j] - tg_yyzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyzzzz_yzzzz_0[j] = pb_x * tg_xyyzzzz_yzzzz_0[j] + fr * tg_xyyzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yzzzz_0[j] - tg_yyzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyzzzz_zzzzz_0[j] = pb_x * tg_xyyzzzz_zzzzz_0[j] + fr * tg_xyyzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_zzzzz_0[j] - tg_yyzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyzzzzz_xxxxx_0[j] = pb_x * tg_xyzzzzz_xxxxx_0[j] + fr * tg_xyzzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxxxx_0[j] - tg_yzzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzzzzz_xxxx_1[j];

                    tg_xxyzzzzz_xxxxy_0[j] = pb_x * tg_xyzzzzz_xxxxy_0[j] + fr * tg_xyzzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxxxy_0[j] - tg_yzzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzzz_xxxy_1[j];

                    tg_xxyzzzzz_xxxxz_0[j] = pb_x * tg_xyzzzzz_xxxxz_0[j] + fr * tg_xyzzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxxxz_0[j] - tg_yzzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzzz_xxxz_1[j];

                    tg_xxyzzzzz_xxxyy_0[j] = pb_x * tg_xyzzzzz_xxxyy_0[j] + fr * tg_xyzzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxxyy_0[j] - tg_yzzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzzz_xxyy_1[j];

                    tg_xxyzzzzz_xxxyz_0[j] = pb_x * tg_xyzzzzz_xxxyz_0[j] + fr * tg_xyzzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxxyz_0[j] - tg_yzzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzzz_xxyz_1[j];

                    tg_xxyzzzzz_xxxzz_0[j] = pb_x * tg_xyzzzzz_xxxzz_0[j] + fr * tg_xyzzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxxzz_0[j] - tg_yzzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzzz_xxzz_1[j];

                    tg_xxyzzzzz_xxyyy_0[j] = pb_x * tg_xyzzzzz_xxyyy_0[j] + fr * tg_xyzzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxyyy_0[j] - tg_yzzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzzz_xyyy_1[j];

                    tg_xxyzzzzz_xxyyz_0[j] = pb_x * tg_xyzzzzz_xxyyz_0[j] + fr * tg_xyzzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxyyz_0[j] - tg_yzzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzzz_xyyz_1[j];

                    tg_xxyzzzzz_xxyzz_0[j] = pb_x * tg_xyzzzzz_xxyzz_0[j] + fr * tg_xyzzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxyzz_0[j] - tg_yzzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzzz_xyzz_1[j];

                    tg_xxyzzzzz_xxzzz_0[j] = pb_x * tg_xyzzzzz_xxzzz_0[j] + fr * tg_xyzzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxzzz_0[j] - tg_yzzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzzz_xzzz_1[j];

                    tg_xxyzzzzz_xyyyy_0[j] = pb_x * tg_xyzzzzz_xyyyy_0[j] + fr * tg_xyzzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xyyyy_0[j] - tg_yzzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_yyyy_1[j];

                    tg_xxyzzzzz_xyyyz_0[j] = pb_x * tg_xyzzzzz_xyyyz_0[j] + fr * tg_xyzzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xyyyz_0[j] - tg_yzzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_yyyz_1[j];

                    tg_xxyzzzzz_xyyzz_0[j] = pb_x * tg_xyzzzzz_xyyzz_0[j] + fr * tg_xyzzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xyyzz_0[j] - tg_yzzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_yyzz_1[j];

                    tg_xxyzzzzz_xyzzz_0[j] = pb_x * tg_xyzzzzz_xyzzz_0[j] + fr * tg_xyzzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xyzzz_0[j] - tg_yzzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_yzzz_1[j];

                    tg_xxyzzzzz_xzzzz_0[j] = pb_x * tg_xyzzzzz_xzzzz_0[j] + fr * tg_xyzzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xzzzz_0[j] - tg_yzzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_zzzz_1[j];

                    tg_xxyzzzzz_yyyyy_0[j] = pb_x * tg_xyzzzzz_yyyyy_0[j] + fr * tg_xyzzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yyyyy_0[j] - tg_yzzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyzzzzz_yyyyz_0[j] = pb_x * tg_xyzzzzz_yyyyz_0[j] + fr * tg_xyzzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yyyyz_0[j] - tg_yzzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyzzzzz_yyyzz_0[j] = pb_x * tg_xyzzzzz_yyyzz_0[j] + fr * tg_xyzzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yyyzz_0[j] - tg_yzzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyzzzzz_yyzzz_0[j] = pb_x * tg_xyzzzzz_yyzzz_0[j] + fr * tg_xyzzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yyzzz_0[j] - tg_yzzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyzzzzz_yzzzz_0[j] = pb_x * tg_xyzzzzz_yzzzz_0[j] + fr * tg_xyzzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yzzzz_0[j] - tg_yzzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyzzzzz_zzzzz_0[j] = pb_x * tg_xyzzzzz_zzzzz_0[j] + fr * tg_xyzzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_zzzzz_0[j] - tg_yzzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xxzzzzzz_xxxxx_0[j] = pb_x * tg_xzzzzzz_xxxxx_0[j] + fr * tg_xzzzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxxx_0[j] - tg_zzzzzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzzzzz_xxxx_1[j];

                    tg_xxzzzzzz_xxxxy_0[j] = pb_x * tg_xzzzzzz_xxxxy_0[j] + fr * tg_xzzzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxxy_0[j] - tg_zzzzzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzzz_xxxy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_569_663(      CMemBlock2D<double>* primBuffer,
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

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyyyzzz_xyyyz_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 662); 

                auto tg_xzzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 569); 

                auto tg_xzzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 570); 

                auto tg_xzzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 571); 

                auto tg_xzzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 572); 

                auto tg_xzzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 573); 

                auto tg_xzzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 574); 

                auto tg_xzzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 575); 

                auto tg_xzzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 576); 

                auto tg_xzzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 577); 

                auto tg_xzzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 578); 

                auto tg_xzzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 579); 

                auto tg_xzzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 580); 

                auto tg_xzzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 581); 

                auto tg_xzzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 582); 

                auto tg_xzzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 583); 

                auto tg_xzzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 584); 

                auto tg_xzzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 585); 

                auto tg_xzzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 586); 

                auto tg_xzzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 587); 

                auto tg_yyyyyyy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 588); 

                auto tg_yyyyyyy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 589); 

                auto tg_yyyyyyy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 590); 

                auto tg_yyyyyyy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 591); 

                auto tg_yyyyyyy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 592); 

                auto tg_yyyyyyy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 593); 

                auto tg_yyyyyyy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 594); 

                auto tg_yyyyyyy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 595); 

                auto tg_yyyyyyy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 596); 

                auto tg_yyyyyyy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 597); 

                auto tg_yyyyyyy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 598); 

                auto tg_yyyyyyy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 599); 

                auto tg_yyyyyyy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 600); 

                auto tg_yyyyyyy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 601); 

                auto tg_yyyyyyy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 602); 

                auto tg_yyyyyyy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 603); 

                auto tg_yyyyyyy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 604); 

                auto tg_yyyyyyy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 605); 

                auto tg_yyyyyyy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 606); 

                auto tg_yyyyyyy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 607); 

                auto tg_yyyyyyy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 608); 

                auto tg_yyyyyyz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 609); 

                auto tg_yyyyyyz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 610); 

                auto tg_yyyyyyz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 611); 

                auto tg_yyyyyyz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 612); 

                auto tg_yyyyyyz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 613); 

                auto tg_yyyyyyz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 614); 

                auto tg_yyyyyyz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 615); 

                auto tg_yyyyyyz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 616); 

                auto tg_yyyyyyz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 617); 

                auto tg_yyyyyyz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 618); 

                auto tg_yyyyyyz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 619); 

                auto tg_yyyyyyz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 620); 

                auto tg_yyyyyyz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 621); 

                auto tg_yyyyyyz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 622); 

                auto tg_yyyyyyz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 623); 

                auto tg_yyyyyyz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 624); 

                auto tg_yyyyyyz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 625); 

                auto tg_yyyyyyz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 626); 

                auto tg_yyyyyyz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 627); 

                auto tg_yyyyyyz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 628); 

                auto tg_yyyyyyz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 629); 

                auto tg_yyyyyzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 630); 

                auto tg_yyyyyzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 631); 

                auto tg_yyyyyzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 632); 

                auto tg_yyyyyzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 633); 

                auto tg_yyyyyzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 634); 

                auto tg_yyyyyzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 635); 

                auto tg_yyyyyzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 636); 

                auto tg_yyyyyzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 637); 

                auto tg_yyyyyzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 638); 

                auto tg_yyyyyzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 639); 

                auto tg_yyyyyzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 640); 

                auto tg_yyyyyzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 641); 

                auto tg_yyyyyzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 642); 

                auto tg_yyyyyzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 643); 

                auto tg_yyyyyzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 644); 

                auto tg_yyyyyzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 645); 

                auto tg_yyyyyzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 646); 

                auto tg_yyyyyzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 647); 

                auto tg_yyyyyzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 648); 

                auto tg_yyyyyzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 649); 

                auto tg_yyyyyzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 650); 

                auto tg_yyyyzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 651); 

                auto tg_yyyyzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 652); 

                auto tg_yyyyzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 653); 

                auto tg_yyyyzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 654); 

                auto tg_yyyyzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 655); 

                auto tg_yyyyzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 656); 

                auto tg_yyyyzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 657); 

                auto tg_yyyyzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 658); 

                auto tg_yyyyzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 659); 

                auto tg_yyyyzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 660); 

                auto tg_yyyyzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 661); 

                auto tg_yyyyzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 662); 

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

                auto tg_xzzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 407); 

                auto tg_xzzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 408); 

                auto tg_xzzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 409); 

                auto tg_xzzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 410); 

                auto tg_xzzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 411); 

                auto tg_xzzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 412); 

                auto tg_xzzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 413); 

                auto tg_xzzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 414); 

                auto tg_xzzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 415); 

                auto tg_xzzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 416); 

                auto tg_xzzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 417); 

                auto tg_xzzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 418); 

                auto tg_xzzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 419); 

                auto tg_yyyyyyy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 420); 

                auto tg_yyyyyyy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 421); 

                auto tg_yyyyyyy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 422); 

                auto tg_yyyyyyy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 423); 

                auto tg_yyyyyyy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 424); 

                auto tg_yyyyyyy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 425); 

                auto tg_yyyyyyy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 426); 

                auto tg_yyyyyyy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 427); 

                auto tg_yyyyyyy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 428); 

                auto tg_yyyyyyy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 429); 

                auto tg_yyyyyyy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 430); 

                auto tg_yyyyyyy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 431); 

                auto tg_yyyyyyy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 432); 

                auto tg_yyyyyyy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 433); 

                auto tg_yyyyyyy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 434); 

                auto tg_yyyyyyz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 435); 

                auto tg_yyyyyyz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 436); 

                auto tg_yyyyyyz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 437); 

                auto tg_yyyyyyz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 438); 

                auto tg_yyyyyyz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 439); 

                auto tg_yyyyyyz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 440); 

                auto tg_yyyyyyz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 441); 

                auto tg_yyyyyyz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 442); 

                auto tg_yyyyyyz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 443); 

                auto tg_yyyyyyz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 444); 

                auto tg_yyyyyyz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 445); 

                auto tg_yyyyyyz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 446); 

                auto tg_yyyyyyz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 447); 

                auto tg_yyyyyyz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 448); 

                auto tg_yyyyyyz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 449); 

                auto tg_yyyyyzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 450); 

                auto tg_yyyyyzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 451); 

                auto tg_yyyyyzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 452); 

                auto tg_yyyyyzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 453); 

                auto tg_yyyyyzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 454); 

                auto tg_yyyyyzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 455); 

                auto tg_yyyyyzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 456); 

                auto tg_yyyyyzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 457); 

                auto tg_yyyyyzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 458); 

                auto tg_yyyyyzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 459); 

                auto tg_yyyyyzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 460); 

                auto tg_yyyyyzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 461); 

                auto tg_yyyyyzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 462); 

                auto tg_yyyyyzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 463); 

                auto tg_yyyyyzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 464); 

                auto tg_yyyyzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 465); 

                auto tg_yyyyzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 466); 

                auto tg_yyyyzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 467); 

                auto tg_yyyyzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 468); 

                auto tg_yyyyzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 469); 

                auto tg_yyyyzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 470); 

                auto tg_yyyyzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 471); 

                auto tg_yyyyzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 472); 

                auto tg_yyyyzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 473); 

                auto tg_yyyyzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 474); 

                auto tg_yyyyzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 475); 

                auto tg_yyyyzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 476); 

                // set up pointers to integrals

                auto tg_xxzzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 569); 

                auto tg_xxzzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 570); 

                auto tg_xxzzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 571); 

                auto tg_xxzzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 572); 

                auto tg_xxzzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 573); 

                auto tg_xxzzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 574); 

                auto tg_xxzzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 575); 

                auto tg_xxzzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 576); 

                auto tg_xxzzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 577); 

                auto tg_xxzzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 578); 

                auto tg_xxzzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 579); 

                auto tg_xxzzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 580); 

                auto tg_xxzzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 581); 

                auto tg_xxzzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 582); 

                auto tg_xxzzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 583); 

                auto tg_xxzzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 584); 

                auto tg_xxzzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 585); 

                auto tg_xxzzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 586); 

                auto tg_xxzzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 587); 

                auto tg_xyyyyyyy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 588); 

                auto tg_xyyyyyyy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 589); 

                auto tg_xyyyyyyy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 590); 

                auto tg_xyyyyyyy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 591); 

                auto tg_xyyyyyyy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 592); 

                auto tg_xyyyyyyy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 593); 

                auto tg_xyyyyyyy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 594); 

                auto tg_xyyyyyyy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 595); 

                auto tg_xyyyyyyy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 596); 

                auto tg_xyyyyyyy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 597); 

                auto tg_xyyyyyyy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 598); 

                auto tg_xyyyyyyy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 599); 

                auto tg_xyyyyyyy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 600); 

                auto tg_xyyyyyyy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 601); 

                auto tg_xyyyyyyy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 602); 

                auto tg_xyyyyyyy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 603); 

                auto tg_xyyyyyyy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 604); 

                auto tg_xyyyyyyy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 605); 

                auto tg_xyyyyyyy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 606); 

                auto tg_xyyyyyyy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 607); 

                auto tg_xyyyyyyy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 608); 

                auto tg_xyyyyyyz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 609); 

                auto tg_xyyyyyyz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 610); 

                auto tg_xyyyyyyz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 611); 

                auto tg_xyyyyyyz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 612); 

                auto tg_xyyyyyyz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 613); 

                auto tg_xyyyyyyz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 614); 

                auto tg_xyyyyyyz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 615); 

                auto tg_xyyyyyyz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 616); 

                auto tg_xyyyyyyz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 617); 

                auto tg_xyyyyyyz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 618); 

                auto tg_xyyyyyyz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 619); 

                auto tg_xyyyyyyz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 620); 

                auto tg_xyyyyyyz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 621); 

                auto tg_xyyyyyyz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 622); 

                auto tg_xyyyyyyz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 623); 

                auto tg_xyyyyyyz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 624); 

                auto tg_xyyyyyyz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 625); 

                auto tg_xyyyyyyz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 626); 

                auto tg_xyyyyyyz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 627); 

                auto tg_xyyyyyyz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 628); 

                auto tg_xyyyyyyz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 629); 

                auto tg_xyyyyyzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 630); 

                auto tg_xyyyyyzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 631); 

                auto tg_xyyyyyzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 632); 

                auto tg_xyyyyyzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 633); 

                auto tg_xyyyyyzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 634); 

                auto tg_xyyyyyzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 635); 

                auto tg_xyyyyyzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 636); 

                auto tg_xyyyyyzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 637); 

                auto tg_xyyyyyzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 638); 

                auto tg_xyyyyyzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 639); 

                auto tg_xyyyyyzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 640); 

                auto tg_xyyyyyzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 641); 

                auto tg_xyyyyyzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 642); 

                auto tg_xyyyyyzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 643); 

                auto tg_xyyyyyzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 644); 

                auto tg_xyyyyyzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 645); 

                auto tg_xyyyyyzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 646); 

                auto tg_xyyyyyzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 647); 

                auto tg_xyyyyyzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 648); 

                auto tg_xyyyyyzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 649); 

                auto tg_xyyyyyzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 650); 

                auto tg_xyyyyzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 651); 

                auto tg_xyyyyzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 652); 

                auto tg_xyyyyzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 653); 

                auto tg_xyyyyzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 654); 

                auto tg_xyyyyzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 655); 

                auto tg_xyyyyzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 656); 

                auto tg_xyyyyzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 657); 

                auto tg_xyyyyzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 658); 

                auto tg_xyyyyzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 659); 

                auto tg_xyyyyzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 660); 

                auto tg_xyyyyzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 661); 

                auto tg_xyyyyzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 662); 

                // Batch of Integrals (569,663)

                #pragma omp simd aligned(fxn, fza, tg_xxzzzzzz_xxxxz_0, tg_xxzzzzzz_xxxyy_0, \
                                         tg_xxzzzzzz_xxxyz_0, tg_xxzzzzzz_xxxzz_0, tg_xxzzzzzz_xxyyy_0, tg_xxzzzzzz_xxyyz_0, \
                                         tg_xxzzzzzz_xxyzz_0, tg_xxzzzzzz_xxzzz_0, tg_xxzzzzzz_xyyyy_0, tg_xxzzzzzz_xyyyz_0, \
                                         tg_xxzzzzzz_xyyzz_0, tg_xxzzzzzz_xyzzz_0, tg_xxzzzzzz_xzzzz_0, tg_xxzzzzzz_yyyyy_0, \
                                         tg_xxzzzzzz_yyyyz_0, tg_xxzzzzzz_yyyzz_0, tg_xxzzzzzz_yyzzz_0, tg_xxzzzzzz_yzzzz_0, \
                                         tg_xxzzzzzz_zzzzz_0, tg_xyyyyyyy_xxxxx_0, tg_xyyyyyyy_xxxxy_0, tg_xyyyyyyy_xxxxz_0, \
                                         tg_xyyyyyyy_xxxyy_0, tg_xyyyyyyy_xxxyz_0, tg_xyyyyyyy_xxxzz_0, tg_xyyyyyyy_xxyyy_0, \
                                         tg_xyyyyyyy_xxyyz_0, tg_xyyyyyyy_xxyzz_0, tg_xyyyyyyy_xxzzz_0, tg_xyyyyyyy_xyyyy_0, \
                                         tg_xyyyyyyy_xyyyz_0, tg_xyyyyyyy_xyyzz_0, tg_xyyyyyyy_xyzzz_0, tg_xyyyyyyy_xzzzz_0, \
                                         tg_xyyyyyyy_yyyyy_0, tg_xyyyyyyy_yyyyz_0, tg_xyyyyyyy_yyyzz_0, tg_xyyyyyyy_yyzzz_0, \
                                         tg_xyyyyyyy_yzzzz_0, tg_xyyyyyyy_zzzzz_0, tg_xyyyyyyz_xxxxx_0, tg_xyyyyyyz_xxxxy_0, \
                                         tg_xyyyyyyz_xxxxz_0, tg_xyyyyyyz_xxxyy_0, tg_xyyyyyyz_xxxyz_0, tg_xyyyyyyz_xxxzz_0, \
                                         tg_xyyyyyyz_xxyyy_0, tg_xyyyyyyz_xxyyz_0, tg_xyyyyyyz_xxyzz_0, tg_xyyyyyyz_xxzzz_0, \
                                         tg_xyyyyyyz_xyyyy_0, tg_xyyyyyyz_xyyyz_0, tg_xyyyyyyz_xyyzz_0, tg_xyyyyyyz_xyzzz_0, \
                                         tg_xyyyyyyz_xzzzz_0, tg_xyyyyyyz_yyyyy_0, tg_xyyyyyyz_yyyyz_0, tg_xyyyyyyz_yyyzz_0, \
                                         tg_xyyyyyyz_yyzzz_0, tg_xyyyyyyz_yzzzz_0, tg_xyyyyyyz_zzzzz_0, tg_xyyyyyzz_xxxxx_0, \
                                         tg_xyyyyyzz_xxxxy_0, tg_xyyyyyzz_xxxxz_0, tg_xyyyyyzz_xxxyy_0, tg_xyyyyyzz_xxxyz_0, \
                                         tg_xyyyyyzz_xxxzz_0, tg_xyyyyyzz_xxyyy_0, tg_xyyyyyzz_xxyyz_0, tg_xyyyyyzz_xxyzz_0, \
                                         tg_xyyyyyzz_xxzzz_0, tg_xyyyyyzz_xyyyy_0, tg_xyyyyyzz_xyyyz_0, tg_xyyyyyzz_xyyzz_0, \
                                         tg_xyyyyyzz_xyzzz_0, tg_xyyyyyzz_xzzzz_0, tg_xyyyyyzz_yyyyy_0, tg_xyyyyyzz_yyyyz_0, \
                                         tg_xyyyyyzz_yyyzz_0, tg_xyyyyyzz_yyzzz_0, tg_xyyyyyzz_yzzzz_0, tg_xyyyyyzz_zzzzz_0, \
                                         tg_xyyyyzzz_xxxxx_0, tg_xyyyyzzz_xxxxy_0, tg_xyyyyzzz_xxxxz_0, tg_xyyyyzzz_xxxyy_0, \
                                         tg_xyyyyzzz_xxxyz_0, tg_xyyyyzzz_xxxzz_0, tg_xyyyyzzz_xxyyy_0, tg_xyyyyzzz_xxyyz_0, \
                                         tg_xyyyyzzz_xxyzz_0, tg_xyyyyzzz_xxzzz_0, tg_xyyyyzzz_xyyyy_0, tg_xyyyyzzz_xyyyz_0, \
                                         tg_xzzzzzz_xxxxz_0, tg_xzzzzzz_xxxxz_1, tg_xzzzzzz_xxxyy_0, tg_xzzzzzz_xxxyy_1, \
                                         tg_xzzzzzz_xxxyz_0, tg_xzzzzzz_xxxyz_1, tg_xzzzzzz_xxxz_1, tg_xzzzzzz_xxxzz_0, \
                                         tg_xzzzzzz_xxxzz_1, tg_xzzzzzz_xxyy_1, tg_xzzzzzz_xxyyy_0, tg_xzzzzzz_xxyyy_1, \
                                         tg_xzzzzzz_xxyyz_0, tg_xzzzzzz_xxyyz_1, tg_xzzzzzz_xxyz_1, tg_xzzzzzz_xxyzz_0, \
                                         tg_xzzzzzz_xxyzz_1, tg_xzzzzzz_xxzz_1, tg_xzzzzzz_xxzzz_0, tg_xzzzzzz_xxzzz_1, \
                                         tg_xzzzzzz_xyyy_1, tg_xzzzzzz_xyyyy_0, tg_xzzzzzz_xyyyy_1, tg_xzzzzzz_xyyyz_0, \
                                         tg_xzzzzzz_xyyyz_1, tg_xzzzzzz_xyyz_1, tg_xzzzzzz_xyyzz_0, tg_xzzzzzz_xyyzz_1, \
                                         tg_xzzzzzz_xyzz_1, tg_xzzzzzz_xyzzz_0, tg_xzzzzzz_xyzzz_1, tg_xzzzzzz_xzzz_1, \
                                         tg_xzzzzzz_xzzzz_0, tg_xzzzzzz_xzzzz_1, tg_xzzzzzz_yyyy_1, tg_xzzzzzz_yyyyy_0, \
                                         tg_xzzzzzz_yyyyy_1, tg_xzzzzzz_yyyyz_0, tg_xzzzzzz_yyyyz_1, tg_xzzzzzz_yyyz_1, \
                                         tg_xzzzzzz_yyyzz_0, tg_xzzzzzz_yyyzz_1, tg_xzzzzzz_yyzz_1, tg_xzzzzzz_yyzzz_0, \
                                         tg_xzzzzzz_yyzzz_1, tg_xzzzzzz_yzzz_1, tg_xzzzzzz_yzzzz_0, tg_xzzzzzz_yzzzz_1, \
                                         tg_xzzzzzz_zzzz_1, tg_xzzzzzz_zzzzz_0, tg_xzzzzzz_zzzzz_1, tg_yyyyyyy_xxxx_1, \
                                         tg_yyyyyyy_xxxxx_0, tg_yyyyyyy_xxxxx_1, tg_yyyyyyy_xxxxy_0, tg_yyyyyyy_xxxxy_1, \
                                         tg_yyyyyyy_xxxxz_0, tg_yyyyyyy_xxxxz_1, tg_yyyyyyy_xxxy_1, tg_yyyyyyy_xxxyy_0, \
                                         tg_yyyyyyy_xxxyy_1, tg_yyyyyyy_xxxyz_0, tg_yyyyyyy_xxxyz_1, tg_yyyyyyy_xxxz_1, \
                                         tg_yyyyyyy_xxxzz_0, tg_yyyyyyy_xxxzz_1, tg_yyyyyyy_xxyy_1, tg_yyyyyyy_xxyyy_0, \
                                         tg_yyyyyyy_xxyyy_1, tg_yyyyyyy_xxyyz_0, tg_yyyyyyy_xxyyz_1, tg_yyyyyyy_xxyz_1, \
                                         tg_yyyyyyy_xxyzz_0, tg_yyyyyyy_xxyzz_1, tg_yyyyyyy_xxzz_1, tg_yyyyyyy_xxzzz_0, \
                                         tg_yyyyyyy_xxzzz_1, tg_yyyyyyy_xyyy_1, tg_yyyyyyy_xyyyy_0, tg_yyyyyyy_xyyyy_1, \
                                         tg_yyyyyyy_xyyyz_0, tg_yyyyyyy_xyyyz_1, tg_yyyyyyy_xyyz_1, tg_yyyyyyy_xyyzz_0, \
                                         tg_yyyyyyy_xyyzz_1, tg_yyyyyyy_xyzz_1, tg_yyyyyyy_xyzzz_0, tg_yyyyyyy_xyzzz_1, \
                                         tg_yyyyyyy_xzzz_1, tg_yyyyyyy_xzzzz_0, tg_yyyyyyy_xzzzz_1, tg_yyyyyyy_yyyy_1, \
                                         tg_yyyyyyy_yyyyy_0, tg_yyyyyyy_yyyyy_1, tg_yyyyyyy_yyyyz_0, tg_yyyyyyy_yyyyz_1, \
                                         tg_yyyyyyy_yyyz_1, tg_yyyyyyy_yyyzz_0, tg_yyyyyyy_yyyzz_1, tg_yyyyyyy_yyzz_1, \
                                         tg_yyyyyyy_yyzzz_0, tg_yyyyyyy_yyzzz_1, tg_yyyyyyy_yzzz_1, tg_yyyyyyy_yzzzz_0, \
                                         tg_yyyyyyy_yzzzz_1, tg_yyyyyyy_zzzz_1, tg_yyyyyyy_zzzzz_0, tg_yyyyyyy_zzzzz_1, \
                                         tg_yyyyyyz_xxxx_1, tg_yyyyyyz_xxxxx_0, tg_yyyyyyz_xxxxx_1, tg_yyyyyyz_xxxxy_0, \
                                         tg_yyyyyyz_xxxxy_1, tg_yyyyyyz_xxxxz_0, tg_yyyyyyz_xxxxz_1, tg_yyyyyyz_xxxy_1, \
                                         tg_yyyyyyz_xxxyy_0, tg_yyyyyyz_xxxyy_1, tg_yyyyyyz_xxxyz_0, tg_yyyyyyz_xxxyz_1, \
                                         tg_yyyyyyz_xxxz_1, tg_yyyyyyz_xxxzz_0, tg_yyyyyyz_xxxzz_1, tg_yyyyyyz_xxyy_1, \
                                         tg_yyyyyyz_xxyyy_0, tg_yyyyyyz_xxyyy_1, tg_yyyyyyz_xxyyz_0, tg_yyyyyyz_xxyyz_1, \
                                         tg_yyyyyyz_xxyz_1, tg_yyyyyyz_xxyzz_0, tg_yyyyyyz_xxyzz_1, tg_yyyyyyz_xxzz_1, \
                                         tg_yyyyyyz_xxzzz_0, tg_yyyyyyz_xxzzz_1, tg_yyyyyyz_xyyy_1, tg_yyyyyyz_xyyyy_0, \
                                         tg_yyyyyyz_xyyyy_1, tg_yyyyyyz_xyyyz_0, tg_yyyyyyz_xyyyz_1, tg_yyyyyyz_xyyz_1, \
                                         tg_yyyyyyz_xyyzz_0, tg_yyyyyyz_xyyzz_1, tg_yyyyyyz_xyzz_1, tg_yyyyyyz_xyzzz_0, \
                                         tg_yyyyyyz_xyzzz_1, tg_yyyyyyz_xzzz_1, tg_yyyyyyz_xzzzz_0, tg_yyyyyyz_xzzzz_1, \
                                         tg_yyyyyyz_yyyy_1, tg_yyyyyyz_yyyyy_0, tg_yyyyyyz_yyyyy_1, tg_yyyyyyz_yyyyz_0, \
                                         tg_yyyyyyz_yyyyz_1, tg_yyyyyyz_yyyz_1, tg_yyyyyyz_yyyzz_0, tg_yyyyyyz_yyyzz_1, \
                                         tg_yyyyyyz_yyzz_1, tg_yyyyyyz_yyzzz_0, tg_yyyyyyz_yyzzz_1, tg_yyyyyyz_yzzz_1, \
                                         tg_yyyyyyz_yzzzz_0, tg_yyyyyyz_yzzzz_1, tg_yyyyyyz_zzzz_1, tg_yyyyyyz_zzzzz_0, \
                                         tg_yyyyyyz_zzzzz_1, tg_yyyyyzz_xxxx_1, tg_yyyyyzz_xxxxx_0, tg_yyyyyzz_xxxxx_1, \
                                         tg_yyyyyzz_xxxxy_0, tg_yyyyyzz_xxxxy_1, tg_yyyyyzz_xxxxz_0, tg_yyyyyzz_xxxxz_1, \
                                         tg_yyyyyzz_xxxy_1, tg_yyyyyzz_xxxyy_0, tg_yyyyyzz_xxxyy_1, tg_yyyyyzz_xxxyz_0, \
                                         tg_yyyyyzz_xxxyz_1, tg_yyyyyzz_xxxz_1, tg_yyyyyzz_xxxzz_0, tg_yyyyyzz_xxxzz_1, \
                                         tg_yyyyyzz_xxyy_1, tg_yyyyyzz_xxyyy_0, tg_yyyyyzz_xxyyy_1, tg_yyyyyzz_xxyyz_0, \
                                         tg_yyyyyzz_xxyyz_1, tg_yyyyyzz_xxyz_1, tg_yyyyyzz_xxyzz_0, tg_yyyyyzz_xxyzz_1, \
                                         tg_yyyyyzz_xxzz_1, tg_yyyyyzz_xxzzz_0, tg_yyyyyzz_xxzzz_1, tg_yyyyyzz_xyyy_1, \
                                         tg_yyyyyzz_xyyyy_0, tg_yyyyyzz_xyyyy_1, tg_yyyyyzz_xyyyz_0, tg_yyyyyzz_xyyyz_1, \
                                         tg_yyyyyzz_xyyz_1, tg_yyyyyzz_xyyzz_0, tg_yyyyyzz_xyyzz_1, tg_yyyyyzz_xyzz_1, \
                                         tg_yyyyyzz_xyzzz_0, tg_yyyyyzz_xyzzz_1, tg_yyyyyzz_xzzz_1, tg_yyyyyzz_xzzzz_0, \
                                         tg_yyyyyzz_xzzzz_1, tg_yyyyyzz_yyyy_1, tg_yyyyyzz_yyyyy_0, tg_yyyyyzz_yyyyy_1, \
                                         tg_yyyyyzz_yyyyz_0, tg_yyyyyzz_yyyyz_1, tg_yyyyyzz_yyyz_1, tg_yyyyyzz_yyyzz_0, \
                                         tg_yyyyyzz_yyyzz_1, tg_yyyyyzz_yyzz_1, tg_yyyyyzz_yyzzz_0, tg_yyyyyzz_yyzzz_1, \
                                         tg_yyyyyzz_yzzz_1, tg_yyyyyzz_yzzzz_0, tg_yyyyyzz_yzzzz_1, tg_yyyyyzz_zzzz_1, \
                                         tg_yyyyyzz_zzzzz_0, tg_yyyyyzz_zzzzz_1, tg_yyyyzzz_xxxx_1, tg_yyyyzzz_xxxxx_0, \
                                         tg_yyyyzzz_xxxxx_1, tg_yyyyzzz_xxxxy_0, tg_yyyyzzz_xxxxy_1, tg_yyyyzzz_xxxxz_0, \
                                         tg_yyyyzzz_xxxxz_1, tg_yyyyzzz_xxxy_1, tg_yyyyzzz_xxxyy_0, tg_yyyyzzz_xxxyy_1, \
                                         tg_yyyyzzz_xxxyz_0, tg_yyyyzzz_xxxyz_1, tg_yyyyzzz_xxxz_1, tg_yyyyzzz_xxxzz_0, \
                                         tg_yyyyzzz_xxxzz_1, tg_yyyyzzz_xxyy_1, tg_yyyyzzz_xxyyy_0, tg_yyyyzzz_xxyyy_1, \
                                         tg_yyyyzzz_xxyyz_0, tg_yyyyzzz_xxyyz_1, tg_yyyyzzz_xxyz_1, tg_yyyyzzz_xxyzz_0, \
                                         tg_yyyyzzz_xxyzz_1, tg_yyyyzzz_xxzz_1, tg_yyyyzzz_xxzzz_0, tg_yyyyzzz_xxzzz_1, \
                                         tg_yyyyzzz_xyyy_1, tg_yyyyzzz_xyyyy_0, tg_yyyyzzz_xyyyy_1, tg_yyyyzzz_xyyyz_0, \
                                         tg_yyyyzzz_xyyyz_1, tg_yyyyzzz_xyyz_1, tg_yyyyzzz_xyzz_1, tg_yyyyzzz_xzzz_1, \
                                         tg_yyyyzzz_yyyy_1, tg_yyyyzzz_yyyz_1, tg_zzzzzz_xxxxz_0, tg_zzzzzz_xxxxz_1, \
                                         tg_zzzzzz_xxxyy_0, tg_zzzzzz_xxxyy_1, tg_zzzzzz_xxxyz_0, tg_zzzzzz_xxxyz_1, \
                                         tg_zzzzzz_xxxzz_0, tg_zzzzzz_xxxzz_1, tg_zzzzzz_xxyyy_0, tg_zzzzzz_xxyyy_1, \
                                         tg_zzzzzz_xxyyz_0, tg_zzzzzz_xxyyz_1, tg_zzzzzz_xxyzz_0, tg_zzzzzz_xxyzz_1, \
                                         tg_zzzzzz_xxzzz_0, tg_zzzzzz_xxzzz_1, tg_zzzzzz_xyyyy_0, tg_zzzzzz_xyyyy_1, \
                                         tg_zzzzzz_xyyyz_0, tg_zzzzzz_xyyyz_1, tg_zzzzzz_xyyzz_0, tg_zzzzzz_xyyzz_1, \
                                         tg_zzzzzz_xyzzz_0, tg_zzzzzz_xyzzz_1, tg_zzzzzz_xzzzz_0, tg_zzzzzz_xzzzz_1, \
                                         tg_zzzzzz_yyyyy_0, tg_zzzzzz_yyyyy_1, tg_zzzzzz_yyyyz_0, tg_zzzzzz_yyyyz_1, \
                                         tg_zzzzzz_yyyzz_0, tg_zzzzzz_yyyzz_1, tg_zzzzzz_yyzzz_0, tg_zzzzzz_yyzzz_1, \
                                         tg_zzzzzz_yzzzz_0, tg_zzzzzz_yzzzz_1, tg_zzzzzz_zzzzz_0, tg_zzzzzz_zzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxzzzzzz_xxxxz_0[j] = pb_x * tg_xzzzzzz_xxxxz_0[j] + fr * tg_xzzzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxxz_0[j] - tg_zzzzzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzzz_xxxz_1[j];

                    tg_xxzzzzzz_xxxyy_0[j] = pb_x * tg_xzzzzzz_xxxyy_0[j] + fr * tg_xzzzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxyy_0[j] - tg_zzzzzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzzz_xxyy_1[j];

                    tg_xxzzzzzz_xxxyz_0[j] = pb_x * tg_xzzzzzz_xxxyz_0[j] + fr * tg_xzzzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxyz_0[j] - tg_zzzzzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzzz_xxyz_1[j];

                    tg_xxzzzzzz_xxxzz_0[j] = pb_x * tg_xzzzzzz_xxxzz_0[j] + fr * tg_xzzzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxzz_0[j] - tg_zzzzzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzzz_xxzz_1[j];

                    tg_xxzzzzzz_xxyyy_0[j] = pb_x * tg_xzzzzzz_xxyyy_0[j] + fr * tg_xzzzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxyyy_0[j] - tg_zzzzzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzzz_xyyy_1[j];

                    tg_xxzzzzzz_xxyyz_0[j] = pb_x * tg_xzzzzzz_xxyyz_0[j] + fr * tg_xzzzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxyyz_0[j] - tg_zzzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzzz_xyyz_1[j];

                    tg_xxzzzzzz_xxyzz_0[j] = pb_x * tg_xzzzzzz_xxyzz_0[j] + fr * tg_xzzzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxyzz_0[j] - tg_zzzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzzz_xyzz_1[j];

                    tg_xxzzzzzz_xxzzz_0[j] = pb_x * tg_xzzzzzz_xxzzz_0[j] + fr * tg_xzzzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxzzz_0[j] - tg_zzzzzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzzz_xzzz_1[j];

                    tg_xxzzzzzz_xyyyy_0[j] = pb_x * tg_xzzzzzz_xyyyy_0[j] + fr * tg_xzzzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyyyy_0[j] - tg_zzzzzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_yyyy_1[j];

                    tg_xxzzzzzz_xyyyz_0[j] = pb_x * tg_xzzzzzz_xyyyz_0[j] + fr * tg_xzzzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyyyz_0[j] - tg_zzzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_yyyz_1[j];

                    tg_xxzzzzzz_xyyzz_0[j] = pb_x * tg_xzzzzzz_xyyzz_0[j] + fr * tg_xzzzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyyzz_0[j] - tg_zzzzzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_yyzz_1[j];

                    tg_xxzzzzzz_xyzzz_0[j] = pb_x * tg_xzzzzzz_xyzzz_0[j] + fr * tg_xzzzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyzzz_0[j] - tg_zzzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_yzzz_1[j];

                    tg_xxzzzzzz_xzzzz_0[j] = pb_x * tg_xzzzzzz_xzzzz_0[j] + fr * tg_xzzzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xzzzz_0[j] - tg_zzzzzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_zzzz_1[j];

                    tg_xxzzzzzz_yyyyy_0[j] = pb_x * tg_xzzzzzz_yyyyy_0[j] + fr * tg_xzzzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyyyy_0[j] - tg_zzzzzz_yyyyy_1[j] * fl1_fza);

                    tg_xxzzzzzz_yyyyz_0[j] = pb_x * tg_xzzzzzz_yyyyz_0[j] + fr * tg_xzzzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyyyz_0[j] - tg_zzzzzz_yyyyz_1[j] * fl1_fza);

                    tg_xxzzzzzz_yyyzz_0[j] = pb_x * tg_xzzzzzz_yyyzz_0[j] + fr * tg_xzzzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyyzz_0[j] - tg_zzzzzz_yyyzz_1[j] * fl1_fza);

                    tg_xxzzzzzz_yyzzz_0[j] = pb_x * tg_xzzzzzz_yyzzz_0[j] + fr * tg_xzzzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyzzz_0[j] - tg_zzzzzz_yyzzz_1[j] * fl1_fza);

                    tg_xxzzzzzz_yzzzz_0[j] = pb_x * tg_xzzzzzz_yzzzz_0[j] + fr * tg_xzzzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yzzzz_0[j] - tg_zzzzzz_yzzzz_1[j] * fl1_fza);

                    tg_xxzzzzzz_zzzzz_0[j] = pb_x * tg_xzzzzzz_zzzzz_0[j] + fr * tg_xzzzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_zzzzz_0[j] - tg_zzzzzz_zzzzz_1[j] * fl1_fza);

                    tg_xyyyyyyy_xxxxx_0[j] = pb_x * tg_yyyyyyy_xxxxx_0[j] + fr * tg_yyyyyyy_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyyyy_xxxx_1[j];

                    tg_xyyyyyyy_xxxxy_0[j] = pb_x * tg_yyyyyyy_xxxxy_0[j] + fr * tg_yyyyyyy_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyyyy_xxxy_1[j];

                    tg_xyyyyyyy_xxxxz_0[j] = pb_x * tg_yyyyyyy_xxxxz_0[j] + fr * tg_yyyyyyy_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyyyy_xxxz_1[j];

                    tg_xyyyyyyy_xxxyy_0[j] = pb_x * tg_yyyyyyy_xxxyy_0[j] + fr * tg_yyyyyyy_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_xxyy_1[j];

                    tg_xyyyyyyy_xxxyz_0[j] = pb_x * tg_yyyyyyy_xxxyz_0[j] + fr * tg_yyyyyyy_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_xxyz_1[j];

                    tg_xyyyyyyy_xxxzz_0[j] = pb_x * tg_yyyyyyy_xxxzz_0[j] + fr * tg_yyyyyyy_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_xxzz_1[j];

                    tg_xyyyyyyy_xxyyy_0[j] = pb_x * tg_yyyyyyy_xxyyy_0[j] + fr * tg_yyyyyyy_xxyyy_1[j] + fl1_fxn * tg_yyyyyyy_xyyy_1[j];

                    tg_xyyyyyyy_xxyyz_0[j] = pb_x * tg_yyyyyyy_xxyyz_0[j] + fr * tg_yyyyyyy_xxyyz_1[j] + fl1_fxn * tg_yyyyyyy_xyyz_1[j];

                    tg_xyyyyyyy_xxyzz_0[j] = pb_x * tg_yyyyyyy_xxyzz_0[j] + fr * tg_yyyyyyy_xxyzz_1[j] + fl1_fxn * tg_yyyyyyy_xyzz_1[j];

                    tg_xyyyyyyy_xxzzz_0[j] = pb_x * tg_yyyyyyy_xxzzz_0[j] + fr * tg_yyyyyyy_xxzzz_1[j] + fl1_fxn * tg_yyyyyyy_xzzz_1[j];

                    tg_xyyyyyyy_xyyyy_0[j] = pb_x * tg_yyyyyyy_xyyyy_0[j] + fr * tg_yyyyyyy_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yyyy_1[j];

                    tg_xyyyyyyy_xyyyz_0[j] = pb_x * tg_yyyyyyy_xyyyz_0[j] + fr * tg_yyyyyyy_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yyyz_1[j];

                    tg_xyyyyyyy_xyyzz_0[j] = pb_x * tg_yyyyyyy_xyyzz_0[j] + fr * tg_yyyyyyy_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yyzz_1[j];

                    tg_xyyyyyyy_xyzzz_0[j] = pb_x * tg_yyyyyyy_xyzzz_0[j] + fr * tg_yyyyyyy_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yzzz_1[j];

                    tg_xyyyyyyy_xzzzz_0[j] = pb_x * tg_yyyyyyy_xzzzz_0[j] + fr * tg_yyyyyyy_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_zzzz_1[j];

                    tg_xyyyyyyy_yyyyy_0[j] = pb_x * tg_yyyyyyy_yyyyy_0[j] + fr * tg_yyyyyyy_yyyyy_1[j];

                    tg_xyyyyyyy_yyyyz_0[j] = pb_x * tg_yyyyyyy_yyyyz_0[j] + fr * tg_yyyyyyy_yyyyz_1[j];

                    tg_xyyyyyyy_yyyzz_0[j] = pb_x * tg_yyyyyyy_yyyzz_0[j] + fr * tg_yyyyyyy_yyyzz_1[j];

                    tg_xyyyyyyy_yyzzz_0[j] = pb_x * tg_yyyyyyy_yyzzz_0[j] + fr * tg_yyyyyyy_yyzzz_1[j];

                    tg_xyyyyyyy_yzzzz_0[j] = pb_x * tg_yyyyyyy_yzzzz_0[j] + fr * tg_yyyyyyy_yzzzz_1[j];

                    tg_xyyyyyyy_zzzzz_0[j] = pb_x * tg_yyyyyyy_zzzzz_0[j] + fr * tg_yyyyyyy_zzzzz_1[j];

                    tg_xyyyyyyz_xxxxx_0[j] = pb_x * tg_yyyyyyz_xxxxx_0[j] + fr * tg_yyyyyyz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyyyz_xxxx_1[j];

                    tg_xyyyyyyz_xxxxy_0[j] = pb_x * tg_yyyyyyz_xxxxy_0[j] + fr * tg_yyyyyyz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyyyz_xxxy_1[j];

                    tg_xyyyyyyz_xxxxz_0[j] = pb_x * tg_yyyyyyz_xxxxz_0[j] + fr * tg_yyyyyyz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyyyz_xxxz_1[j];

                    tg_xyyyyyyz_xxxyy_0[j] = pb_x * tg_yyyyyyz_xxxyy_0[j] + fr * tg_yyyyyyz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_xxyy_1[j];

                    tg_xyyyyyyz_xxxyz_0[j] = pb_x * tg_yyyyyyz_xxxyz_0[j] + fr * tg_yyyyyyz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_xxyz_1[j];

                    tg_xyyyyyyz_xxxzz_0[j] = pb_x * tg_yyyyyyz_xxxzz_0[j] + fr * tg_yyyyyyz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_xxzz_1[j];

                    tg_xyyyyyyz_xxyyy_0[j] = pb_x * tg_yyyyyyz_xxyyy_0[j] + fr * tg_yyyyyyz_xxyyy_1[j] + fl1_fxn * tg_yyyyyyz_xyyy_1[j];

                    tg_xyyyyyyz_xxyyz_0[j] = pb_x * tg_yyyyyyz_xxyyz_0[j] + fr * tg_yyyyyyz_xxyyz_1[j] + fl1_fxn * tg_yyyyyyz_xyyz_1[j];

                    tg_xyyyyyyz_xxyzz_0[j] = pb_x * tg_yyyyyyz_xxyzz_0[j] + fr * tg_yyyyyyz_xxyzz_1[j] + fl1_fxn * tg_yyyyyyz_xyzz_1[j];

                    tg_xyyyyyyz_xxzzz_0[j] = pb_x * tg_yyyyyyz_xxzzz_0[j] + fr * tg_yyyyyyz_xxzzz_1[j] + fl1_fxn * tg_yyyyyyz_xzzz_1[j];

                    tg_xyyyyyyz_xyyyy_0[j] = pb_x * tg_yyyyyyz_xyyyy_0[j] + fr * tg_yyyyyyz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yyyy_1[j];

                    tg_xyyyyyyz_xyyyz_0[j] = pb_x * tg_yyyyyyz_xyyyz_0[j] + fr * tg_yyyyyyz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yyyz_1[j];

                    tg_xyyyyyyz_xyyzz_0[j] = pb_x * tg_yyyyyyz_xyyzz_0[j] + fr * tg_yyyyyyz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yyzz_1[j];

                    tg_xyyyyyyz_xyzzz_0[j] = pb_x * tg_yyyyyyz_xyzzz_0[j] + fr * tg_yyyyyyz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yzzz_1[j];

                    tg_xyyyyyyz_xzzzz_0[j] = pb_x * tg_yyyyyyz_xzzzz_0[j] + fr * tg_yyyyyyz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_zzzz_1[j];

                    tg_xyyyyyyz_yyyyy_0[j] = pb_x * tg_yyyyyyz_yyyyy_0[j] + fr * tg_yyyyyyz_yyyyy_1[j];

                    tg_xyyyyyyz_yyyyz_0[j] = pb_x * tg_yyyyyyz_yyyyz_0[j] + fr * tg_yyyyyyz_yyyyz_1[j];

                    tg_xyyyyyyz_yyyzz_0[j] = pb_x * tg_yyyyyyz_yyyzz_0[j] + fr * tg_yyyyyyz_yyyzz_1[j];

                    tg_xyyyyyyz_yyzzz_0[j] = pb_x * tg_yyyyyyz_yyzzz_0[j] + fr * tg_yyyyyyz_yyzzz_1[j];

                    tg_xyyyyyyz_yzzzz_0[j] = pb_x * tg_yyyyyyz_yzzzz_0[j] + fr * tg_yyyyyyz_yzzzz_1[j];

                    tg_xyyyyyyz_zzzzz_0[j] = pb_x * tg_yyyyyyz_zzzzz_0[j] + fr * tg_yyyyyyz_zzzzz_1[j];

                    tg_xyyyyyzz_xxxxx_0[j] = pb_x * tg_yyyyyzz_xxxxx_0[j] + fr * tg_yyyyyzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyyzz_xxxx_1[j];

                    tg_xyyyyyzz_xxxxy_0[j] = pb_x * tg_yyyyyzz_xxxxy_0[j] + fr * tg_yyyyyzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyyzz_xxxy_1[j];

                    tg_xyyyyyzz_xxxxz_0[j] = pb_x * tg_yyyyyzz_xxxxz_0[j] + fr * tg_yyyyyzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyyzz_xxxz_1[j];

                    tg_xyyyyyzz_xxxyy_0[j] = pb_x * tg_yyyyyzz_xxxyy_0[j] + fr * tg_yyyyyzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_xxyy_1[j];

                    tg_xyyyyyzz_xxxyz_0[j] = pb_x * tg_yyyyyzz_xxxyz_0[j] + fr * tg_yyyyyzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_xxyz_1[j];

                    tg_xyyyyyzz_xxxzz_0[j] = pb_x * tg_yyyyyzz_xxxzz_0[j] + fr * tg_yyyyyzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_xxzz_1[j];

                    tg_xyyyyyzz_xxyyy_0[j] = pb_x * tg_yyyyyzz_xxyyy_0[j] + fr * tg_yyyyyzz_xxyyy_1[j] + fl1_fxn * tg_yyyyyzz_xyyy_1[j];

                    tg_xyyyyyzz_xxyyz_0[j] = pb_x * tg_yyyyyzz_xxyyz_0[j] + fr * tg_yyyyyzz_xxyyz_1[j] + fl1_fxn * tg_yyyyyzz_xyyz_1[j];

                    tg_xyyyyyzz_xxyzz_0[j] = pb_x * tg_yyyyyzz_xxyzz_0[j] + fr * tg_yyyyyzz_xxyzz_1[j] + fl1_fxn * tg_yyyyyzz_xyzz_1[j];

                    tg_xyyyyyzz_xxzzz_0[j] = pb_x * tg_yyyyyzz_xxzzz_0[j] + fr * tg_yyyyyzz_xxzzz_1[j] + fl1_fxn * tg_yyyyyzz_xzzz_1[j];

                    tg_xyyyyyzz_xyyyy_0[j] = pb_x * tg_yyyyyzz_xyyyy_0[j] + fr * tg_yyyyyzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yyyy_1[j];

                    tg_xyyyyyzz_xyyyz_0[j] = pb_x * tg_yyyyyzz_xyyyz_0[j] + fr * tg_yyyyyzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yyyz_1[j];

                    tg_xyyyyyzz_xyyzz_0[j] = pb_x * tg_yyyyyzz_xyyzz_0[j] + fr * tg_yyyyyzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yyzz_1[j];

                    tg_xyyyyyzz_xyzzz_0[j] = pb_x * tg_yyyyyzz_xyzzz_0[j] + fr * tg_yyyyyzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yzzz_1[j];

                    tg_xyyyyyzz_xzzzz_0[j] = pb_x * tg_yyyyyzz_xzzzz_0[j] + fr * tg_yyyyyzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_zzzz_1[j];

                    tg_xyyyyyzz_yyyyy_0[j] = pb_x * tg_yyyyyzz_yyyyy_0[j] + fr * tg_yyyyyzz_yyyyy_1[j];

                    tg_xyyyyyzz_yyyyz_0[j] = pb_x * tg_yyyyyzz_yyyyz_0[j] + fr * tg_yyyyyzz_yyyyz_1[j];

                    tg_xyyyyyzz_yyyzz_0[j] = pb_x * tg_yyyyyzz_yyyzz_0[j] + fr * tg_yyyyyzz_yyyzz_1[j];

                    tg_xyyyyyzz_yyzzz_0[j] = pb_x * tg_yyyyyzz_yyzzz_0[j] + fr * tg_yyyyyzz_yyzzz_1[j];

                    tg_xyyyyyzz_yzzzz_0[j] = pb_x * tg_yyyyyzz_yzzzz_0[j] + fr * tg_yyyyyzz_yzzzz_1[j];

                    tg_xyyyyyzz_zzzzz_0[j] = pb_x * tg_yyyyyzz_zzzzz_0[j] + fr * tg_yyyyyzz_zzzzz_1[j];

                    tg_xyyyyzzz_xxxxx_0[j] = pb_x * tg_yyyyzzz_xxxxx_0[j] + fr * tg_yyyyzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyyzzz_xxxx_1[j];

                    tg_xyyyyzzz_xxxxy_0[j] = pb_x * tg_yyyyzzz_xxxxy_0[j] + fr * tg_yyyyzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyyzzz_xxxy_1[j];

                    tg_xyyyyzzz_xxxxz_0[j] = pb_x * tg_yyyyzzz_xxxxz_0[j] + fr * tg_yyyyzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyyzzz_xxxz_1[j];

                    tg_xyyyyzzz_xxxyy_0[j] = pb_x * tg_yyyyzzz_xxxyy_0[j] + fr * tg_yyyyzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_xxyy_1[j];

                    tg_xyyyyzzz_xxxyz_0[j] = pb_x * tg_yyyyzzz_xxxyz_0[j] + fr * tg_yyyyzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_xxyz_1[j];

                    tg_xyyyyzzz_xxxzz_0[j] = pb_x * tg_yyyyzzz_xxxzz_0[j] + fr * tg_yyyyzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_xxzz_1[j];

                    tg_xyyyyzzz_xxyyy_0[j] = pb_x * tg_yyyyzzz_xxyyy_0[j] + fr * tg_yyyyzzz_xxyyy_1[j] + fl1_fxn * tg_yyyyzzz_xyyy_1[j];

                    tg_xyyyyzzz_xxyyz_0[j] = pb_x * tg_yyyyzzz_xxyyz_0[j] + fr * tg_yyyyzzz_xxyyz_1[j] + fl1_fxn * tg_yyyyzzz_xyyz_1[j];

                    tg_xyyyyzzz_xxyzz_0[j] = pb_x * tg_yyyyzzz_xxyzz_0[j] + fr * tg_yyyyzzz_xxyzz_1[j] + fl1_fxn * tg_yyyyzzz_xyzz_1[j];

                    tg_xyyyyzzz_xxzzz_0[j] = pb_x * tg_yyyyzzz_xxzzz_0[j] + fr * tg_yyyyzzz_xxzzz_1[j] + fl1_fxn * tg_yyyyzzz_xzzz_1[j];

                    tg_xyyyyzzz_xyyyy_0[j] = pb_x * tg_yyyyzzz_xyyyy_0[j] + fr * tg_yyyyzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yyyy_1[j];

                    tg_xyyyyzzz_xyyyz_0[j] = pb_x * tg_yyyyzzz_xyyyz_0[j] + fr * tg_yyyyzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_663_757(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyyyyyy_xxxxx_0 = primBuffer[pidx_g_7_5_m0].data(756 * idx + 588); 

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

                auto tg_yyyyyyy_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 588); 

                auto tg_yyyyzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 663); 

                auto tg_yyyyzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 664); 

                auto tg_yyyyzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 665); 

                auto tg_yyyyzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 666); 

                auto tg_yyyyzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 667); 

                auto tg_yyyyzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 668); 

                auto tg_yyyyzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 669); 

                auto tg_yyyyzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 670); 

                auto tg_yyyyzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 671); 

                auto tg_yyyzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 672); 

                auto tg_yyyzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 673); 

                auto tg_yyyzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 674); 

                auto tg_yyyzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 675); 

                auto tg_yyyzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 676); 

                auto tg_yyyzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 677); 

                auto tg_yyyzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 678); 

                auto tg_yyyzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 679); 

                auto tg_yyyzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 680); 

                auto tg_yyyzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 681); 

                auto tg_yyyzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 682); 

                auto tg_yyyzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 683); 

                auto tg_yyyzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 684); 

                auto tg_yyyzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 685); 

                auto tg_yyyzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 686); 

                auto tg_yyyzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 687); 

                auto tg_yyyzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 688); 

                auto tg_yyyzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 689); 

                auto tg_yyyzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 690); 

                auto tg_yyyzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 691); 

                auto tg_yyyzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 692); 

                auto tg_yyzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 693); 

                auto tg_yyzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 694); 

                auto tg_yyzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 695); 

                auto tg_yyzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 696); 

                auto tg_yyzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 697); 

                auto tg_yyzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 698); 

                auto tg_yyzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 699); 

                auto tg_yyzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 700); 

                auto tg_yyzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 701); 

                auto tg_yyzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 702); 

                auto tg_yyzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 703); 

                auto tg_yyzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 704); 

                auto tg_yyzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 705); 

                auto tg_yyzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 706); 

                auto tg_yyzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 707); 

                auto tg_yyzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 708); 

                auto tg_yyzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 709); 

                auto tg_yyzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 710); 

                auto tg_yyzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 711); 

                auto tg_yyzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 712); 

                auto tg_yyzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 713); 

                auto tg_yzzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 714); 

                auto tg_yzzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 715); 

                auto tg_yzzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 716); 

                auto tg_yzzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 717); 

                auto tg_yzzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 718); 

                auto tg_yzzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 719); 

                auto tg_yzzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 720); 

                auto tg_yzzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 721); 

                auto tg_yzzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 722); 

                auto tg_yzzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 723); 

                auto tg_yzzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 724); 

                auto tg_yzzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 725); 

                auto tg_yzzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 726); 

                auto tg_yzzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 727); 

                auto tg_yzzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 728); 

                auto tg_yzzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 729); 

                auto tg_yzzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 730); 

                auto tg_yzzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 731); 

                auto tg_yzzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 732); 

                auto tg_yzzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 733); 

                auto tg_yzzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 734); 

                auto tg_zzzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 735); 

                auto tg_zzzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 736); 

                auto tg_zzzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 737); 

                auto tg_zzzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 738); 

                auto tg_zzzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 739); 

                auto tg_zzzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 740); 

                auto tg_zzzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 741); 

                auto tg_zzzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 742); 

                auto tg_zzzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 743); 

                auto tg_zzzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 744); 

                auto tg_zzzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 745); 

                auto tg_zzzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 746); 

                auto tg_zzzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 747); 

                auto tg_zzzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 748); 

                auto tg_zzzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 749); 

                auto tg_zzzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 750); 

                auto tg_zzzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 751); 

                auto tg_zzzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 752); 

                auto tg_zzzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 753); 

                auto tg_zzzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 754); 

                auto tg_zzzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 755); 

                auto tg_yyyyyy_xxxxx_0 = primBuffer[pidx_g_6_5_m0].data(588 * idx + 441); 

                auto tg_yyyyyy_xxxxx_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 441); 

                auto tg_yyyyzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 477); 

                auto tg_yyyyzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 478); 

                auto tg_yyyyzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 479); 

                auto tg_yyyzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 480); 

                auto tg_yyyzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 481); 

                auto tg_yyyzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 482); 

                auto tg_yyyzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 483); 

                auto tg_yyyzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 484); 

                auto tg_yyyzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 485); 

                auto tg_yyyzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 486); 

                auto tg_yyyzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 487); 

                auto tg_yyyzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 488); 

                auto tg_yyyzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 489); 

                auto tg_yyyzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 490); 

                auto tg_yyyzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 491); 

                auto tg_yyyzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 492); 

                auto tg_yyyzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 493); 

                auto tg_yyyzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 494); 

                auto tg_yyzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 495); 

                auto tg_yyzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 496); 

                auto tg_yyzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 497); 

                auto tg_yyzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 498); 

                auto tg_yyzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 499); 

                auto tg_yyzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 500); 

                auto tg_yyzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 501); 

                auto tg_yyzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 502); 

                auto tg_yyzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 503); 

                auto tg_yyzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 504); 

                auto tg_yyzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 505); 

                auto tg_yyzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 506); 

                auto tg_yyzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 507); 

                auto tg_yyzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 508); 

                auto tg_yyzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 509); 

                auto tg_yzzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 510); 

                auto tg_yzzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 511); 

                auto tg_yzzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 512); 

                auto tg_yzzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 513); 

                auto tg_yzzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 514); 

                auto tg_yzzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 515); 

                auto tg_yzzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 516); 

                auto tg_yzzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 517); 

                auto tg_yzzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 518); 

                auto tg_yzzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 519); 

                auto tg_yzzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 520); 

                auto tg_yzzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 521); 

                auto tg_yzzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 522); 

                auto tg_yzzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 523); 

                auto tg_yzzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 524); 

                auto tg_zzzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 525); 

                auto tg_zzzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 526); 

                auto tg_zzzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 527); 

                auto tg_zzzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 528); 

                auto tg_zzzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 529); 

                auto tg_zzzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 530); 

                auto tg_zzzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 531); 

                auto tg_zzzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 532); 

                auto tg_zzzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 533); 

                auto tg_zzzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 534); 

                auto tg_zzzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 535); 

                auto tg_zzzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 536); 

                auto tg_zzzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 537); 

                auto tg_zzzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 538); 

                auto tg_zzzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 539); 

                // set up pointers to integrals

                auto tg_xyyyyzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 663); 

                auto tg_xyyyyzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 664); 

                auto tg_xyyyyzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 665); 

                auto tg_xyyyyzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 666); 

                auto tg_xyyyyzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 667); 

                auto tg_xyyyyzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 668); 

                auto tg_xyyyyzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 669); 

                auto tg_xyyyyzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 670); 

                auto tg_xyyyyzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 671); 

                auto tg_xyyyzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 672); 

                auto tg_xyyyzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 673); 

                auto tg_xyyyzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 674); 

                auto tg_xyyyzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 675); 

                auto tg_xyyyzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 676); 

                auto tg_xyyyzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 677); 

                auto tg_xyyyzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 678); 

                auto tg_xyyyzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 679); 

                auto tg_xyyyzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 680); 

                auto tg_xyyyzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 681); 

                auto tg_xyyyzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 682); 

                auto tg_xyyyzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 683); 

                auto tg_xyyyzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 684); 

                auto tg_xyyyzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 685); 

                auto tg_xyyyzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 686); 

                auto tg_xyyyzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 687); 

                auto tg_xyyyzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 688); 

                auto tg_xyyyzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 689); 

                auto tg_xyyyzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 690); 

                auto tg_xyyyzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 691); 

                auto tg_xyyyzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 692); 

                auto tg_xyyzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 693); 

                auto tg_xyyzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 694); 

                auto tg_xyyzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 695); 

                auto tg_xyyzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 696); 

                auto tg_xyyzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 697); 

                auto tg_xyyzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 698); 

                auto tg_xyyzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 699); 

                auto tg_xyyzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 700); 

                auto tg_xyyzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 701); 

                auto tg_xyyzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 702); 

                auto tg_xyyzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 703); 

                auto tg_xyyzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 704); 

                auto tg_xyyzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 705); 

                auto tg_xyyzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 706); 

                auto tg_xyyzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 707); 

                auto tg_xyyzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 708); 

                auto tg_xyyzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 709); 

                auto tg_xyyzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 710); 

                auto tg_xyyzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 711); 

                auto tg_xyyzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 712); 

                auto tg_xyyzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 713); 

                auto tg_xyzzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 714); 

                auto tg_xyzzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 715); 

                auto tg_xyzzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 716); 

                auto tg_xyzzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 717); 

                auto tg_xyzzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 718); 

                auto tg_xyzzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 719); 

                auto tg_xyzzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 720); 

                auto tg_xyzzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 721); 

                auto tg_xyzzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 722); 

                auto tg_xyzzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 723); 

                auto tg_xyzzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 724); 

                auto tg_xyzzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 725); 

                auto tg_xyzzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 726); 

                auto tg_xyzzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 727); 

                auto tg_xyzzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 728); 

                auto tg_xyzzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 729); 

                auto tg_xyzzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 730); 

                auto tg_xyzzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 731); 

                auto tg_xyzzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 732); 

                auto tg_xyzzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 733); 

                auto tg_xyzzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 734); 

                auto tg_xzzzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 735); 

                auto tg_xzzzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 736); 

                auto tg_xzzzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 737); 

                auto tg_xzzzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 738); 

                auto tg_xzzzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 739); 

                auto tg_xzzzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 740); 

                auto tg_xzzzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 741); 

                auto tg_xzzzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 742); 

                auto tg_xzzzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 743); 

                auto tg_xzzzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 744); 

                auto tg_xzzzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 745); 

                auto tg_xzzzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 746); 

                auto tg_xzzzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 747); 

                auto tg_xzzzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 748); 

                auto tg_xzzzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 749); 

                auto tg_xzzzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 750); 

                auto tg_xzzzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 751); 

                auto tg_xzzzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 752); 

                auto tg_xzzzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 753); 

                auto tg_xzzzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 754); 

                auto tg_xzzzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 755); 

                auto tg_yyyyyyyy_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 756); 

                // Batch of Integrals (663,757)

                #pragma omp simd aligned(fxn, fza, tg_xyyyyzzz_xyyzz_0, tg_xyyyyzzz_xyzzz_0, \
                                         tg_xyyyyzzz_xzzzz_0, tg_xyyyyzzz_yyyyy_0, tg_xyyyyzzz_yyyyz_0, tg_xyyyyzzz_yyyzz_0, \
                                         tg_xyyyyzzz_yyzzz_0, tg_xyyyyzzz_yzzzz_0, tg_xyyyyzzz_zzzzz_0, tg_xyyyzzzz_xxxxx_0, \
                                         tg_xyyyzzzz_xxxxy_0, tg_xyyyzzzz_xxxxz_0, tg_xyyyzzzz_xxxyy_0, tg_xyyyzzzz_xxxyz_0, \
                                         tg_xyyyzzzz_xxxzz_0, tg_xyyyzzzz_xxyyy_0, tg_xyyyzzzz_xxyyz_0, tg_xyyyzzzz_xxyzz_0, \
                                         tg_xyyyzzzz_xxzzz_0, tg_xyyyzzzz_xyyyy_0, tg_xyyyzzzz_xyyyz_0, tg_xyyyzzzz_xyyzz_0, \
                                         tg_xyyyzzzz_xyzzz_0, tg_xyyyzzzz_xzzzz_0, tg_xyyyzzzz_yyyyy_0, tg_xyyyzzzz_yyyyz_0, \
                                         tg_xyyyzzzz_yyyzz_0, tg_xyyyzzzz_yyzzz_0, tg_xyyyzzzz_yzzzz_0, tg_xyyyzzzz_zzzzz_0, \
                                         tg_xyyzzzzz_xxxxx_0, tg_xyyzzzzz_xxxxy_0, tg_xyyzzzzz_xxxxz_0, tg_xyyzzzzz_xxxyy_0, \
                                         tg_xyyzzzzz_xxxyz_0, tg_xyyzzzzz_xxxzz_0, tg_xyyzzzzz_xxyyy_0, tg_xyyzzzzz_xxyyz_0, \
                                         tg_xyyzzzzz_xxyzz_0, tg_xyyzzzzz_xxzzz_0, tg_xyyzzzzz_xyyyy_0, tg_xyyzzzzz_xyyyz_0, \
                                         tg_xyyzzzzz_xyyzz_0, tg_xyyzzzzz_xyzzz_0, tg_xyyzzzzz_xzzzz_0, tg_xyyzzzzz_yyyyy_0, \
                                         tg_xyyzzzzz_yyyyz_0, tg_xyyzzzzz_yyyzz_0, tg_xyyzzzzz_yyzzz_0, tg_xyyzzzzz_yzzzz_0, \
                                         tg_xyyzzzzz_zzzzz_0, tg_xyzzzzzz_xxxxx_0, tg_xyzzzzzz_xxxxy_0, tg_xyzzzzzz_xxxxz_0, \
                                         tg_xyzzzzzz_xxxyy_0, tg_xyzzzzzz_xxxyz_0, tg_xyzzzzzz_xxxzz_0, tg_xyzzzzzz_xxyyy_0, \
                                         tg_xyzzzzzz_xxyyz_0, tg_xyzzzzzz_xxyzz_0, tg_xyzzzzzz_xxzzz_0, tg_xyzzzzzz_xyyyy_0, \
                                         tg_xyzzzzzz_xyyyz_0, tg_xyzzzzzz_xyyzz_0, tg_xyzzzzzz_xyzzz_0, tg_xyzzzzzz_xzzzz_0, \
                                         tg_xyzzzzzz_yyyyy_0, tg_xyzzzzzz_yyyyz_0, tg_xyzzzzzz_yyyzz_0, tg_xyzzzzzz_yyzzz_0, \
                                         tg_xyzzzzzz_yzzzz_0, tg_xyzzzzzz_zzzzz_0, tg_xzzzzzzz_xxxxx_0, tg_xzzzzzzz_xxxxy_0, \
                                         tg_xzzzzzzz_xxxxz_0, tg_xzzzzzzz_xxxyy_0, tg_xzzzzzzz_xxxyz_0, tg_xzzzzzzz_xxxzz_0, \
                                         tg_xzzzzzzz_xxyyy_0, tg_xzzzzzzz_xxyyz_0, tg_xzzzzzzz_xxyzz_0, tg_xzzzzzzz_xxzzz_0, \
                                         tg_xzzzzzzz_xyyyy_0, tg_xzzzzzzz_xyyyz_0, tg_xzzzzzzz_xyyzz_0, tg_xzzzzzzz_xyzzz_0, \
                                         tg_xzzzzzzz_xzzzz_0, tg_xzzzzzzz_yyyyy_0, tg_xzzzzzzz_yyyyz_0, tg_xzzzzzzz_yyyzz_0, \
                                         tg_xzzzzzzz_yyzzz_0, tg_xzzzzzzz_yzzzz_0, tg_xzzzzzzz_zzzzz_0, tg_yyyyyy_xxxxx_0, \
                                         tg_yyyyyy_xxxxx_1, tg_yyyyyyy_xxxxx_0, tg_yyyyyyy_xxxxx_1, tg_yyyyyyyy_xxxxx_0, \
                                         tg_yyyyzzz_xyyzz_0, tg_yyyyzzz_xyyzz_1, tg_yyyyzzz_xyzzz_0, tg_yyyyzzz_xyzzz_1, \
                                         tg_yyyyzzz_xzzzz_0, tg_yyyyzzz_xzzzz_1, tg_yyyyzzz_yyyyy_0, tg_yyyyzzz_yyyyy_1, \
                                         tg_yyyyzzz_yyyyz_0, tg_yyyyzzz_yyyyz_1, tg_yyyyzzz_yyyzz_0, tg_yyyyzzz_yyyzz_1, \
                                         tg_yyyyzzz_yyzz_1, tg_yyyyzzz_yyzzz_0, tg_yyyyzzz_yyzzz_1, tg_yyyyzzz_yzzz_1, \
                                         tg_yyyyzzz_yzzzz_0, tg_yyyyzzz_yzzzz_1, tg_yyyyzzz_zzzz_1, tg_yyyyzzz_zzzzz_0, \
                                         tg_yyyyzzz_zzzzz_1, tg_yyyzzzz_xxxx_1, tg_yyyzzzz_xxxxx_0, tg_yyyzzzz_xxxxx_1, \
                                         tg_yyyzzzz_xxxxy_0, tg_yyyzzzz_xxxxy_1, tg_yyyzzzz_xxxxz_0, tg_yyyzzzz_xxxxz_1, \
                                         tg_yyyzzzz_xxxy_1, tg_yyyzzzz_xxxyy_0, tg_yyyzzzz_xxxyy_1, tg_yyyzzzz_xxxyz_0, \
                                         tg_yyyzzzz_xxxyz_1, tg_yyyzzzz_xxxz_1, tg_yyyzzzz_xxxzz_0, tg_yyyzzzz_xxxzz_1, \
                                         tg_yyyzzzz_xxyy_1, tg_yyyzzzz_xxyyy_0, tg_yyyzzzz_xxyyy_1, tg_yyyzzzz_xxyyz_0, \
                                         tg_yyyzzzz_xxyyz_1, tg_yyyzzzz_xxyz_1, tg_yyyzzzz_xxyzz_0, tg_yyyzzzz_xxyzz_1, \
                                         tg_yyyzzzz_xxzz_1, tg_yyyzzzz_xxzzz_0, tg_yyyzzzz_xxzzz_1, tg_yyyzzzz_xyyy_1, \
                                         tg_yyyzzzz_xyyyy_0, tg_yyyzzzz_xyyyy_1, tg_yyyzzzz_xyyyz_0, tg_yyyzzzz_xyyyz_1, \
                                         tg_yyyzzzz_xyyz_1, tg_yyyzzzz_xyyzz_0, tg_yyyzzzz_xyyzz_1, tg_yyyzzzz_xyzz_1, \
                                         tg_yyyzzzz_xyzzz_0, tg_yyyzzzz_xyzzz_1, tg_yyyzzzz_xzzz_1, tg_yyyzzzz_xzzzz_0, \
                                         tg_yyyzzzz_xzzzz_1, tg_yyyzzzz_yyyy_1, tg_yyyzzzz_yyyyy_0, tg_yyyzzzz_yyyyy_1, \
                                         tg_yyyzzzz_yyyyz_0, tg_yyyzzzz_yyyyz_1, tg_yyyzzzz_yyyz_1, tg_yyyzzzz_yyyzz_0, \
                                         tg_yyyzzzz_yyyzz_1, tg_yyyzzzz_yyzz_1, tg_yyyzzzz_yyzzz_0, tg_yyyzzzz_yyzzz_1, \
                                         tg_yyyzzzz_yzzz_1, tg_yyyzzzz_yzzzz_0, tg_yyyzzzz_yzzzz_1, tg_yyyzzzz_zzzz_1, \
                                         tg_yyyzzzz_zzzzz_0, tg_yyyzzzz_zzzzz_1, tg_yyzzzzz_xxxx_1, tg_yyzzzzz_xxxxx_0, \
                                         tg_yyzzzzz_xxxxx_1, tg_yyzzzzz_xxxxy_0, tg_yyzzzzz_xxxxy_1, tg_yyzzzzz_xxxxz_0, \
                                         tg_yyzzzzz_xxxxz_1, tg_yyzzzzz_xxxy_1, tg_yyzzzzz_xxxyy_0, tg_yyzzzzz_xxxyy_1, \
                                         tg_yyzzzzz_xxxyz_0, tg_yyzzzzz_xxxyz_1, tg_yyzzzzz_xxxz_1, tg_yyzzzzz_xxxzz_0, \
                                         tg_yyzzzzz_xxxzz_1, tg_yyzzzzz_xxyy_1, tg_yyzzzzz_xxyyy_0, tg_yyzzzzz_xxyyy_1, \
                                         tg_yyzzzzz_xxyyz_0, tg_yyzzzzz_xxyyz_1, tg_yyzzzzz_xxyz_1, tg_yyzzzzz_xxyzz_0, \
                                         tg_yyzzzzz_xxyzz_1, tg_yyzzzzz_xxzz_1, tg_yyzzzzz_xxzzz_0, tg_yyzzzzz_xxzzz_1, \
                                         tg_yyzzzzz_xyyy_1, tg_yyzzzzz_xyyyy_0, tg_yyzzzzz_xyyyy_1, tg_yyzzzzz_xyyyz_0, \
                                         tg_yyzzzzz_xyyyz_1, tg_yyzzzzz_xyyz_1, tg_yyzzzzz_xyyzz_0, tg_yyzzzzz_xyyzz_1, \
                                         tg_yyzzzzz_xyzz_1, tg_yyzzzzz_xyzzz_0, tg_yyzzzzz_xyzzz_1, tg_yyzzzzz_xzzz_1, \
                                         tg_yyzzzzz_xzzzz_0, tg_yyzzzzz_xzzzz_1, tg_yyzzzzz_yyyy_1, tg_yyzzzzz_yyyyy_0, \
                                         tg_yyzzzzz_yyyyy_1, tg_yyzzzzz_yyyyz_0, tg_yyzzzzz_yyyyz_1, tg_yyzzzzz_yyyz_1, \
                                         tg_yyzzzzz_yyyzz_0, tg_yyzzzzz_yyyzz_1, tg_yyzzzzz_yyzz_1, tg_yyzzzzz_yyzzz_0, \
                                         tg_yyzzzzz_yyzzz_1, tg_yyzzzzz_yzzz_1, tg_yyzzzzz_yzzzz_0, tg_yyzzzzz_yzzzz_1, \
                                         tg_yyzzzzz_zzzz_1, tg_yyzzzzz_zzzzz_0, tg_yyzzzzz_zzzzz_1, tg_yzzzzzz_xxxx_1, \
                                         tg_yzzzzzz_xxxxx_0, tg_yzzzzzz_xxxxx_1, tg_yzzzzzz_xxxxy_0, tg_yzzzzzz_xxxxy_1, \
                                         tg_yzzzzzz_xxxxz_0, tg_yzzzzzz_xxxxz_1, tg_yzzzzzz_xxxy_1, tg_yzzzzzz_xxxyy_0, \
                                         tg_yzzzzzz_xxxyy_1, tg_yzzzzzz_xxxyz_0, tg_yzzzzzz_xxxyz_1, tg_yzzzzzz_xxxz_1, \
                                         tg_yzzzzzz_xxxzz_0, tg_yzzzzzz_xxxzz_1, tg_yzzzzzz_xxyy_1, tg_yzzzzzz_xxyyy_0, \
                                         tg_yzzzzzz_xxyyy_1, tg_yzzzzzz_xxyyz_0, tg_yzzzzzz_xxyyz_1, tg_yzzzzzz_xxyz_1, \
                                         tg_yzzzzzz_xxyzz_0, tg_yzzzzzz_xxyzz_1, tg_yzzzzzz_xxzz_1, tg_yzzzzzz_xxzzz_0, \
                                         tg_yzzzzzz_xxzzz_1, tg_yzzzzzz_xyyy_1, tg_yzzzzzz_xyyyy_0, tg_yzzzzzz_xyyyy_1, \
                                         tg_yzzzzzz_xyyyz_0, tg_yzzzzzz_xyyyz_1, tg_yzzzzzz_xyyz_1, tg_yzzzzzz_xyyzz_0, \
                                         tg_yzzzzzz_xyyzz_1, tg_yzzzzzz_xyzz_1, tg_yzzzzzz_xyzzz_0, tg_yzzzzzz_xyzzz_1, \
                                         tg_yzzzzzz_xzzz_1, tg_yzzzzzz_xzzzz_0, tg_yzzzzzz_xzzzz_1, tg_yzzzzzz_yyyy_1, \
                                         tg_yzzzzzz_yyyyy_0, tg_yzzzzzz_yyyyy_1, tg_yzzzzzz_yyyyz_0, tg_yzzzzzz_yyyyz_1, \
                                         tg_yzzzzzz_yyyz_1, tg_yzzzzzz_yyyzz_0, tg_yzzzzzz_yyyzz_1, tg_yzzzzzz_yyzz_1, \
                                         tg_yzzzzzz_yyzzz_0, tg_yzzzzzz_yyzzz_1, tg_yzzzzzz_yzzz_1, tg_yzzzzzz_yzzzz_0, \
                                         tg_yzzzzzz_yzzzz_1, tg_yzzzzzz_zzzz_1, tg_yzzzzzz_zzzzz_0, tg_yzzzzzz_zzzzz_1, \
                                         tg_zzzzzzz_xxxx_1, tg_zzzzzzz_xxxxx_0, tg_zzzzzzz_xxxxx_1, tg_zzzzzzz_xxxxy_0, \
                                         tg_zzzzzzz_xxxxy_1, tg_zzzzzzz_xxxxz_0, tg_zzzzzzz_xxxxz_1, tg_zzzzzzz_xxxy_1, \
                                         tg_zzzzzzz_xxxyy_0, tg_zzzzzzz_xxxyy_1, tg_zzzzzzz_xxxyz_0, tg_zzzzzzz_xxxyz_1, \
                                         tg_zzzzzzz_xxxz_1, tg_zzzzzzz_xxxzz_0, tg_zzzzzzz_xxxzz_1, tg_zzzzzzz_xxyy_1, \
                                         tg_zzzzzzz_xxyyy_0, tg_zzzzzzz_xxyyy_1, tg_zzzzzzz_xxyyz_0, tg_zzzzzzz_xxyyz_1, \
                                         tg_zzzzzzz_xxyz_1, tg_zzzzzzz_xxyzz_0, tg_zzzzzzz_xxyzz_1, tg_zzzzzzz_xxzz_1, \
                                         tg_zzzzzzz_xxzzz_0, tg_zzzzzzz_xxzzz_1, tg_zzzzzzz_xyyy_1, tg_zzzzzzz_xyyyy_0, \
                                         tg_zzzzzzz_xyyyy_1, tg_zzzzzzz_xyyyz_0, tg_zzzzzzz_xyyyz_1, tg_zzzzzzz_xyyz_1, \
                                         tg_zzzzzzz_xyyzz_0, tg_zzzzzzz_xyyzz_1, tg_zzzzzzz_xyzz_1, tg_zzzzzzz_xyzzz_0, \
                                         tg_zzzzzzz_xyzzz_1, tg_zzzzzzz_xzzz_1, tg_zzzzzzz_xzzzz_0, tg_zzzzzzz_xzzzz_1, \
                                         tg_zzzzzzz_yyyy_1, tg_zzzzzzz_yyyyy_0, tg_zzzzzzz_yyyyy_1, tg_zzzzzzz_yyyyz_0, \
                                         tg_zzzzzzz_yyyyz_1, tg_zzzzzzz_yyyz_1, tg_zzzzzzz_yyyzz_0, tg_zzzzzzz_yyyzz_1, \
                                         tg_zzzzzzz_yyzz_1, tg_zzzzzzz_yyzzz_0, tg_zzzzzzz_yyzzz_1, tg_zzzzzzz_yzzz_1, \
                                         tg_zzzzzzz_yzzzz_0, tg_zzzzzzz_yzzzz_1, tg_zzzzzzz_zzzz_1, tg_zzzzzzz_zzzzz_0, \
                                         tg_zzzzzzz_zzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyyyyzzz_xyyzz_0[j] = pb_x * tg_yyyyzzz_xyyzz_0[j] + fr * tg_yyyyzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yyzz_1[j];

                    tg_xyyyyzzz_xyzzz_0[j] = pb_x * tg_yyyyzzz_xyzzz_0[j] + fr * tg_yyyyzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yzzz_1[j];

                    tg_xyyyyzzz_xzzzz_0[j] = pb_x * tg_yyyyzzz_xzzzz_0[j] + fr * tg_yyyyzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_zzzz_1[j];

                    tg_xyyyyzzz_yyyyy_0[j] = pb_x * tg_yyyyzzz_yyyyy_0[j] + fr * tg_yyyyzzz_yyyyy_1[j];

                    tg_xyyyyzzz_yyyyz_0[j] = pb_x * tg_yyyyzzz_yyyyz_0[j] + fr * tg_yyyyzzz_yyyyz_1[j];

                    tg_xyyyyzzz_yyyzz_0[j] = pb_x * tg_yyyyzzz_yyyzz_0[j] + fr * tg_yyyyzzz_yyyzz_1[j];

                    tg_xyyyyzzz_yyzzz_0[j] = pb_x * tg_yyyyzzz_yyzzz_0[j] + fr * tg_yyyyzzz_yyzzz_1[j];

                    tg_xyyyyzzz_yzzzz_0[j] = pb_x * tg_yyyyzzz_yzzzz_0[j] + fr * tg_yyyyzzz_yzzzz_1[j];

                    tg_xyyyyzzz_zzzzz_0[j] = pb_x * tg_yyyyzzz_zzzzz_0[j] + fr * tg_yyyyzzz_zzzzz_1[j];

                    tg_xyyyzzzz_xxxxx_0[j] = pb_x * tg_yyyzzzz_xxxxx_0[j] + fr * tg_yyyzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyzzzz_xxxx_1[j];

                    tg_xyyyzzzz_xxxxy_0[j] = pb_x * tg_yyyzzzz_xxxxy_0[j] + fr * tg_yyyzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyzzzz_xxxy_1[j];

                    tg_xyyyzzzz_xxxxz_0[j] = pb_x * tg_yyyzzzz_xxxxz_0[j] + fr * tg_yyyzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyzzzz_xxxz_1[j];

                    tg_xyyyzzzz_xxxyy_0[j] = pb_x * tg_yyyzzzz_xxxyy_0[j] + fr * tg_yyyzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_xxyy_1[j];

                    tg_xyyyzzzz_xxxyz_0[j] = pb_x * tg_yyyzzzz_xxxyz_0[j] + fr * tg_yyyzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_xxyz_1[j];

                    tg_xyyyzzzz_xxxzz_0[j] = pb_x * tg_yyyzzzz_xxxzz_0[j] + fr * tg_yyyzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_xxzz_1[j];

                    tg_xyyyzzzz_xxyyy_0[j] = pb_x * tg_yyyzzzz_xxyyy_0[j] + fr * tg_yyyzzzz_xxyyy_1[j] + fl1_fxn * tg_yyyzzzz_xyyy_1[j];

                    tg_xyyyzzzz_xxyyz_0[j] = pb_x * tg_yyyzzzz_xxyyz_0[j] + fr * tg_yyyzzzz_xxyyz_1[j] + fl1_fxn * tg_yyyzzzz_xyyz_1[j];

                    tg_xyyyzzzz_xxyzz_0[j] = pb_x * tg_yyyzzzz_xxyzz_0[j] + fr * tg_yyyzzzz_xxyzz_1[j] + fl1_fxn * tg_yyyzzzz_xyzz_1[j];

                    tg_xyyyzzzz_xxzzz_0[j] = pb_x * tg_yyyzzzz_xxzzz_0[j] + fr * tg_yyyzzzz_xxzzz_1[j] + fl1_fxn * tg_yyyzzzz_xzzz_1[j];

                    tg_xyyyzzzz_xyyyy_0[j] = pb_x * tg_yyyzzzz_xyyyy_0[j] + fr * tg_yyyzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yyyy_1[j];

                    tg_xyyyzzzz_xyyyz_0[j] = pb_x * tg_yyyzzzz_xyyyz_0[j] + fr * tg_yyyzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yyyz_1[j];

                    tg_xyyyzzzz_xyyzz_0[j] = pb_x * tg_yyyzzzz_xyyzz_0[j] + fr * tg_yyyzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yyzz_1[j];

                    tg_xyyyzzzz_xyzzz_0[j] = pb_x * tg_yyyzzzz_xyzzz_0[j] + fr * tg_yyyzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yzzz_1[j];

                    tg_xyyyzzzz_xzzzz_0[j] = pb_x * tg_yyyzzzz_xzzzz_0[j] + fr * tg_yyyzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_zzzz_1[j];

                    tg_xyyyzzzz_yyyyy_0[j] = pb_x * tg_yyyzzzz_yyyyy_0[j] + fr * tg_yyyzzzz_yyyyy_1[j];

                    tg_xyyyzzzz_yyyyz_0[j] = pb_x * tg_yyyzzzz_yyyyz_0[j] + fr * tg_yyyzzzz_yyyyz_1[j];

                    tg_xyyyzzzz_yyyzz_0[j] = pb_x * tg_yyyzzzz_yyyzz_0[j] + fr * tg_yyyzzzz_yyyzz_1[j];

                    tg_xyyyzzzz_yyzzz_0[j] = pb_x * tg_yyyzzzz_yyzzz_0[j] + fr * tg_yyyzzzz_yyzzz_1[j];

                    tg_xyyyzzzz_yzzzz_0[j] = pb_x * tg_yyyzzzz_yzzzz_0[j] + fr * tg_yyyzzzz_yzzzz_1[j];

                    tg_xyyyzzzz_zzzzz_0[j] = pb_x * tg_yyyzzzz_zzzzz_0[j] + fr * tg_yyyzzzz_zzzzz_1[j];

                    tg_xyyzzzzz_xxxxx_0[j] = pb_x * tg_yyzzzzz_xxxxx_0[j] + fr * tg_yyzzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyzzzzz_xxxx_1[j];

                    tg_xyyzzzzz_xxxxy_0[j] = pb_x * tg_yyzzzzz_xxxxy_0[j] + fr * tg_yyzzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyzzzzz_xxxy_1[j];

                    tg_xyyzzzzz_xxxxz_0[j] = pb_x * tg_yyzzzzz_xxxxz_0[j] + fr * tg_yyzzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyzzzzz_xxxz_1[j];

                    tg_xyyzzzzz_xxxyy_0[j] = pb_x * tg_yyzzzzz_xxxyy_0[j] + fr * tg_yyzzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_xxyy_1[j];

                    tg_xyyzzzzz_xxxyz_0[j] = pb_x * tg_yyzzzzz_xxxyz_0[j] + fr * tg_yyzzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_xxyz_1[j];

                    tg_xyyzzzzz_xxxzz_0[j] = pb_x * tg_yyzzzzz_xxxzz_0[j] + fr * tg_yyzzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_xxzz_1[j];

                    tg_xyyzzzzz_xxyyy_0[j] = pb_x * tg_yyzzzzz_xxyyy_0[j] + fr * tg_yyzzzzz_xxyyy_1[j] + fl1_fxn * tg_yyzzzzz_xyyy_1[j];

                    tg_xyyzzzzz_xxyyz_0[j] = pb_x * tg_yyzzzzz_xxyyz_0[j] + fr * tg_yyzzzzz_xxyyz_1[j] + fl1_fxn * tg_yyzzzzz_xyyz_1[j];

                    tg_xyyzzzzz_xxyzz_0[j] = pb_x * tg_yyzzzzz_xxyzz_0[j] + fr * tg_yyzzzzz_xxyzz_1[j] + fl1_fxn * tg_yyzzzzz_xyzz_1[j];

                    tg_xyyzzzzz_xxzzz_0[j] = pb_x * tg_yyzzzzz_xxzzz_0[j] + fr * tg_yyzzzzz_xxzzz_1[j] + fl1_fxn * tg_yyzzzzz_xzzz_1[j];

                    tg_xyyzzzzz_xyyyy_0[j] = pb_x * tg_yyzzzzz_xyyyy_0[j] + fr * tg_yyzzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yyyy_1[j];

                    tg_xyyzzzzz_xyyyz_0[j] = pb_x * tg_yyzzzzz_xyyyz_0[j] + fr * tg_yyzzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yyyz_1[j];

                    tg_xyyzzzzz_xyyzz_0[j] = pb_x * tg_yyzzzzz_xyyzz_0[j] + fr * tg_yyzzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yyzz_1[j];

                    tg_xyyzzzzz_xyzzz_0[j] = pb_x * tg_yyzzzzz_xyzzz_0[j] + fr * tg_yyzzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yzzz_1[j];

                    tg_xyyzzzzz_xzzzz_0[j] = pb_x * tg_yyzzzzz_xzzzz_0[j] + fr * tg_yyzzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_zzzz_1[j];

                    tg_xyyzzzzz_yyyyy_0[j] = pb_x * tg_yyzzzzz_yyyyy_0[j] + fr * tg_yyzzzzz_yyyyy_1[j];

                    tg_xyyzzzzz_yyyyz_0[j] = pb_x * tg_yyzzzzz_yyyyz_0[j] + fr * tg_yyzzzzz_yyyyz_1[j];

                    tg_xyyzzzzz_yyyzz_0[j] = pb_x * tg_yyzzzzz_yyyzz_0[j] + fr * tg_yyzzzzz_yyyzz_1[j];

                    tg_xyyzzzzz_yyzzz_0[j] = pb_x * tg_yyzzzzz_yyzzz_0[j] + fr * tg_yyzzzzz_yyzzz_1[j];

                    tg_xyyzzzzz_yzzzz_0[j] = pb_x * tg_yyzzzzz_yzzzz_0[j] + fr * tg_yyzzzzz_yzzzz_1[j];

                    tg_xyyzzzzz_zzzzz_0[j] = pb_x * tg_yyzzzzz_zzzzz_0[j] + fr * tg_yyzzzzz_zzzzz_1[j];

                    tg_xyzzzzzz_xxxxx_0[j] = pb_x * tg_yzzzzzz_xxxxx_0[j] + fr * tg_yzzzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yzzzzzz_xxxx_1[j];

                    tg_xyzzzzzz_xxxxy_0[j] = pb_x * tg_yzzzzzz_xxxxy_0[j] + fr * tg_yzzzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yzzzzzz_xxxy_1[j];

                    tg_xyzzzzzz_xxxxz_0[j] = pb_x * tg_yzzzzzz_xxxxz_0[j] + fr * tg_yzzzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yzzzzzz_xxxz_1[j];

                    tg_xyzzzzzz_xxxyy_0[j] = pb_x * tg_yzzzzzz_xxxyy_0[j] + fr * tg_yzzzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_xxyy_1[j];

                    tg_xyzzzzzz_xxxyz_0[j] = pb_x * tg_yzzzzzz_xxxyz_0[j] + fr * tg_yzzzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_xxyz_1[j];

                    tg_xyzzzzzz_xxxzz_0[j] = pb_x * tg_yzzzzzz_xxxzz_0[j] + fr * tg_yzzzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_xxzz_1[j];

                    tg_xyzzzzzz_xxyyy_0[j] = pb_x * tg_yzzzzzz_xxyyy_0[j] + fr * tg_yzzzzzz_xxyyy_1[j] + fl1_fxn * tg_yzzzzzz_xyyy_1[j];

                    tg_xyzzzzzz_xxyyz_0[j] = pb_x * tg_yzzzzzz_xxyyz_0[j] + fr * tg_yzzzzzz_xxyyz_1[j] + fl1_fxn * tg_yzzzzzz_xyyz_1[j];

                    tg_xyzzzzzz_xxyzz_0[j] = pb_x * tg_yzzzzzz_xxyzz_0[j] + fr * tg_yzzzzzz_xxyzz_1[j] + fl1_fxn * tg_yzzzzzz_xyzz_1[j];

                    tg_xyzzzzzz_xxzzz_0[j] = pb_x * tg_yzzzzzz_xxzzz_0[j] + fr * tg_yzzzzzz_xxzzz_1[j] + fl1_fxn * tg_yzzzzzz_xzzz_1[j];

                    tg_xyzzzzzz_xyyyy_0[j] = pb_x * tg_yzzzzzz_xyyyy_0[j] + fr * tg_yzzzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yyyy_1[j];

                    tg_xyzzzzzz_xyyyz_0[j] = pb_x * tg_yzzzzzz_xyyyz_0[j] + fr * tg_yzzzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yyyz_1[j];

                    tg_xyzzzzzz_xyyzz_0[j] = pb_x * tg_yzzzzzz_xyyzz_0[j] + fr * tg_yzzzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yyzz_1[j];

                    tg_xyzzzzzz_xyzzz_0[j] = pb_x * tg_yzzzzzz_xyzzz_0[j] + fr * tg_yzzzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yzzz_1[j];

                    tg_xyzzzzzz_xzzzz_0[j] = pb_x * tg_yzzzzzz_xzzzz_0[j] + fr * tg_yzzzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_zzzz_1[j];

                    tg_xyzzzzzz_yyyyy_0[j] = pb_x * tg_yzzzzzz_yyyyy_0[j] + fr * tg_yzzzzzz_yyyyy_1[j];

                    tg_xyzzzzzz_yyyyz_0[j] = pb_x * tg_yzzzzzz_yyyyz_0[j] + fr * tg_yzzzzzz_yyyyz_1[j];

                    tg_xyzzzzzz_yyyzz_0[j] = pb_x * tg_yzzzzzz_yyyzz_0[j] + fr * tg_yzzzzzz_yyyzz_1[j];

                    tg_xyzzzzzz_yyzzz_0[j] = pb_x * tg_yzzzzzz_yyzzz_0[j] + fr * tg_yzzzzzz_yyzzz_1[j];

                    tg_xyzzzzzz_yzzzz_0[j] = pb_x * tg_yzzzzzz_yzzzz_0[j] + fr * tg_yzzzzzz_yzzzz_1[j];

                    tg_xyzzzzzz_zzzzz_0[j] = pb_x * tg_yzzzzzz_zzzzz_0[j] + fr * tg_yzzzzzz_zzzzz_1[j];

                    tg_xzzzzzzz_xxxxx_0[j] = pb_x * tg_zzzzzzz_xxxxx_0[j] + fr * tg_zzzzzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_zzzzzzz_xxxx_1[j];

                    tg_xzzzzzzz_xxxxy_0[j] = pb_x * tg_zzzzzzz_xxxxy_0[j] + fr * tg_zzzzzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_zzzzzzz_xxxy_1[j];

                    tg_xzzzzzzz_xxxxz_0[j] = pb_x * tg_zzzzzzz_xxxxz_0[j] + fr * tg_zzzzzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_zzzzzzz_xxxz_1[j];

                    tg_xzzzzzzz_xxxyy_0[j] = pb_x * tg_zzzzzzz_xxxyy_0[j] + fr * tg_zzzzzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xxyy_1[j];

                    tg_xzzzzzzz_xxxyz_0[j] = pb_x * tg_zzzzzzz_xxxyz_0[j] + fr * tg_zzzzzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xxyz_1[j];

                    tg_xzzzzzzz_xxxzz_0[j] = pb_x * tg_zzzzzzz_xxxzz_0[j] + fr * tg_zzzzzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xxzz_1[j];

                    tg_xzzzzzzz_xxyyy_0[j] = pb_x * tg_zzzzzzz_xxyyy_0[j] + fr * tg_zzzzzzz_xxyyy_1[j] + fl1_fxn * tg_zzzzzzz_xyyy_1[j];

                    tg_xzzzzzzz_xxyyz_0[j] = pb_x * tg_zzzzzzz_xxyyz_0[j] + fr * tg_zzzzzzz_xxyyz_1[j] + fl1_fxn * tg_zzzzzzz_xyyz_1[j];

                    tg_xzzzzzzz_xxyzz_0[j] = pb_x * tg_zzzzzzz_xxyzz_0[j] + fr * tg_zzzzzzz_xxyzz_1[j] + fl1_fxn * tg_zzzzzzz_xyzz_1[j];

                    tg_xzzzzzzz_xxzzz_0[j] = pb_x * tg_zzzzzzz_xxzzz_0[j] + fr * tg_zzzzzzz_xxzzz_1[j] + fl1_fxn * tg_zzzzzzz_xzzz_1[j];

                    tg_xzzzzzzz_xyyyy_0[j] = pb_x * tg_zzzzzzz_xyyyy_0[j] + fr * tg_zzzzzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yyyy_1[j];

                    tg_xzzzzzzz_xyyyz_0[j] = pb_x * tg_zzzzzzz_xyyyz_0[j] + fr * tg_zzzzzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yyyz_1[j];

                    tg_xzzzzzzz_xyyzz_0[j] = pb_x * tg_zzzzzzz_xyyzz_0[j] + fr * tg_zzzzzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yyzz_1[j];

                    tg_xzzzzzzz_xyzzz_0[j] = pb_x * tg_zzzzzzz_xyzzz_0[j] + fr * tg_zzzzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yzzz_1[j];

                    tg_xzzzzzzz_xzzzz_0[j] = pb_x * tg_zzzzzzz_xzzzz_0[j] + fr * tg_zzzzzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_zzzz_1[j];

                    tg_xzzzzzzz_yyyyy_0[j] = pb_x * tg_zzzzzzz_yyyyy_0[j] + fr * tg_zzzzzzz_yyyyy_1[j];

                    tg_xzzzzzzz_yyyyz_0[j] = pb_x * tg_zzzzzzz_yyyyz_0[j] + fr * tg_zzzzzzz_yyyyz_1[j];

                    tg_xzzzzzzz_yyyzz_0[j] = pb_x * tg_zzzzzzz_yyyzz_0[j] + fr * tg_zzzzzzz_yyyzz_1[j];

                    tg_xzzzzzzz_yyzzz_0[j] = pb_x * tg_zzzzzzz_yyzzz_0[j] + fr * tg_zzzzzzz_yyzzz_1[j];

                    tg_xzzzzzzz_yzzzz_0[j] = pb_x * tg_zzzzzzz_yzzzz_0[j] + fr * tg_zzzzzzz_yzzzz_1[j];

                    tg_xzzzzzzz_zzzzz_0[j] = pb_x * tg_zzzzzzz_zzzzz_0[j] + fr * tg_zzzzzzz_zzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyyyyy_xxxxx_0[j] = pb_y * tg_yyyyyyy_xxxxx_0[j] + fr * tg_yyyyyyy_xxxxx_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxxxx_0[j] - tg_yyyyyy_xxxxx_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_757_851(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyyyyyy_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 589); 

                auto tg_yyyyyyy_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 590); 

                auto tg_yyyyyyy_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 591); 

                auto tg_yyyyyyy_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 592); 

                auto tg_yyyyyyy_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 593); 

                auto tg_yyyyyyy_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 594); 

                auto tg_yyyyyyy_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 595); 

                auto tg_yyyyyyy_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 596); 

                auto tg_yyyyyyy_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 597); 

                auto tg_yyyyyyy_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 598); 

                auto tg_yyyyyyy_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 599); 

                auto tg_yyyyyyy_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 600); 

                auto tg_yyyyyyy_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 601); 

                auto tg_yyyyyyy_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 602); 

                auto tg_yyyyyyy_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 603); 

                auto tg_yyyyyyy_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 604); 

                auto tg_yyyyyyy_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 605); 

                auto tg_yyyyyyy_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 606); 

                auto tg_yyyyyyy_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 607); 

                auto tg_yyyyyyy_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 608); 

                auto tg_yyyyyyz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 609); 

                auto tg_yyyyyyz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 610); 

                auto tg_yyyyyyz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 611); 

                auto tg_yyyyyyz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 612); 

                auto tg_yyyyyyz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 613); 

                auto tg_yyyyyyz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 614); 

                auto tg_yyyyyyz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 615); 

                auto tg_yyyyyyz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 616); 

                auto tg_yyyyyyz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 617); 

                auto tg_yyyyyyz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 618); 

                auto tg_yyyyyyz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 619); 

                auto tg_yyyyyyz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 620); 

                auto tg_yyyyyyz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 621); 

                auto tg_yyyyyyz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 622); 

                auto tg_yyyyyyz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 623); 

                auto tg_yyyyyyz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 624); 

                auto tg_yyyyyyz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 625); 

                auto tg_yyyyyyz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 626); 

                auto tg_yyyyyyz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 627); 

                auto tg_yyyyyyz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 628); 

                auto tg_yyyyyyz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 629); 

                auto tg_yyyyyzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 630); 

                auto tg_yyyyyzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 631); 

                auto tg_yyyyyzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 632); 

                auto tg_yyyyyzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 633); 

                auto tg_yyyyyzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 634); 

                auto tg_yyyyyzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 635); 

                auto tg_yyyyyzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 636); 

                auto tg_yyyyyzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 637); 

                auto tg_yyyyyzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 638); 

                auto tg_yyyyyzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 639); 

                auto tg_yyyyyzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 640); 

                auto tg_yyyyyzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 641); 

                auto tg_yyyyyzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 642); 

                auto tg_yyyyyzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 643); 

                auto tg_yyyyyzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 644); 

                auto tg_yyyyyzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 645); 

                auto tg_yyyyyzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 646); 

                auto tg_yyyyyzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 647); 

                auto tg_yyyyyzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 648); 

                auto tg_yyyyyzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 649); 

                auto tg_yyyyyzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 650); 

                auto tg_yyyyzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 651); 

                auto tg_yyyyzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 652); 

                auto tg_yyyyzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 653); 

                auto tg_yyyyzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 654); 

                auto tg_yyyyzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 655); 

                auto tg_yyyyzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 656); 

                auto tg_yyyyzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 657); 

                auto tg_yyyyzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 658); 

                auto tg_yyyyzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 659); 

                auto tg_yyyyzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 660); 

                auto tg_yyyyzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 661); 

                auto tg_yyyyzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 662); 

                auto tg_yyyyzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 663); 

                auto tg_yyyyzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 664); 

                auto tg_yyyyzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 665); 

                auto tg_yyyyzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 666); 

                auto tg_yyyyzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 667); 

                auto tg_yyyyzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 668); 

                auto tg_yyyyzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 669); 

                auto tg_yyyyzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 670); 

                auto tg_yyyyzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 671); 

                auto tg_yyyzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 672); 

                auto tg_yyyzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 673); 

                auto tg_yyyzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 674); 

                auto tg_yyyzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 675); 

                auto tg_yyyzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 676); 

                auto tg_yyyzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 677); 

                auto tg_yyyzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 678); 

                auto tg_yyyzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 679); 

                auto tg_yyyzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 680); 

                auto tg_yyyzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 681); 

                auto tg_yyyzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 682); 

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

                auto tg_yyyyyyy_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 420); 

                auto tg_yyyyyyy_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 421); 

                auto tg_yyyyyyy_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 422); 

                auto tg_yyyyyyy_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 423); 

                auto tg_yyyyyyy_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 424); 

                auto tg_yyyyyyy_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 425); 

                auto tg_yyyyyyy_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 426); 

                auto tg_yyyyyyy_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 427); 

                auto tg_yyyyyyy_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 428); 

                auto tg_yyyyyyy_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 429); 

                auto tg_yyyyyyy_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 430); 

                auto tg_yyyyyyy_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 431); 

                auto tg_yyyyyyy_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 432); 

                auto tg_yyyyyyy_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 433); 

                auto tg_yyyyyyy_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 434); 

                auto tg_yyyyyyz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 435); 

                auto tg_yyyyyyz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 436); 

                auto tg_yyyyyyz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 437); 

                auto tg_yyyyyyz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 438); 

                auto tg_yyyyyyz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 439); 

                auto tg_yyyyyyz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 440); 

                auto tg_yyyyyyz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 441); 

                auto tg_yyyyyyz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 442); 

                auto tg_yyyyyyz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 443); 

                auto tg_yyyyyyz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 444); 

                auto tg_yyyyyyz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 445); 

                auto tg_yyyyyyz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 446); 

                auto tg_yyyyyyz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 447); 

                auto tg_yyyyyyz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 448); 

                auto tg_yyyyyyz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 449); 

                auto tg_yyyyyzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 450); 

                auto tg_yyyyyzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 451); 

                auto tg_yyyyyzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 452); 

                auto tg_yyyyyzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 453); 

                auto tg_yyyyyzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 454); 

                auto tg_yyyyyzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 455); 

                auto tg_yyyyyzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 456); 

                auto tg_yyyyyzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 457); 

                auto tg_yyyyyzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 458); 

                auto tg_yyyyyzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 459); 

                auto tg_yyyyyzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 460); 

                auto tg_yyyyyzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 461); 

                auto tg_yyyyyzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 462); 

                auto tg_yyyyyzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 463); 

                auto tg_yyyyyzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 464); 

                auto tg_yyyyzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 465); 

                auto tg_yyyyzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 466); 

                auto tg_yyyyzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 467); 

                auto tg_yyyyzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 468); 

                auto tg_yyyyzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 469); 

                auto tg_yyyyzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 470); 

                auto tg_yyyyzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 471); 

                auto tg_yyyyzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 472); 

                auto tg_yyyyzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 473); 

                auto tg_yyyyzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 474); 

                auto tg_yyyyzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 475); 

                auto tg_yyyyzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 476); 

                auto tg_yyyyzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 477); 

                auto tg_yyyyzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 478); 

                auto tg_yyyyzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 479); 

                auto tg_yyyzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 480); 

                auto tg_yyyzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 481); 

                auto tg_yyyzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 482); 

                auto tg_yyyzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 483); 

                auto tg_yyyzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 484); 

                auto tg_yyyzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 485); 

                auto tg_yyyzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 486); 

                // set up pointers to integrals

                auto tg_yyyyyyyy_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 757); 

                auto tg_yyyyyyyy_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 758); 

                auto tg_yyyyyyyy_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 759); 

                auto tg_yyyyyyyy_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 760); 

                auto tg_yyyyyyyy_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 761); 

                auto tg_yyyyyyyy_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 762); 

                auto tg_yyyyyyyy_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 763); 

                auto tg_yyyyyyyy_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 764); 

                auto tg_yyyyyyyy_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 765); 

                auto tg_yyyyyyyy_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 766); 

                auto tg_yyyyyyyy_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 767); 

                auto tg_yyyyyyyy_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 768); 

                auto tg_yyyyyyyy_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 769); 

                auto tg_yyyyyyyy_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 770); 

                auto tg_yyyyyyyy_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 771); 

                auto tg_yyyyyyyy_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 772); 

                auto tg_yyyyyyyy_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 773); 

                auto tg_yyyyyyyy_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 774); 

                auto tg_yyyyyyyy_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 775); 

                auto tg_yyyyyyyy_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 776); 

                auto tg_yyyyyyyz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 777); 

                auto tg_yyyyyyyz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 778); 

                auto tg_yyyyyyyz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 779); 

                auto tg_yyyyyyyz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 780); 

                auto tg_yyyyyyyz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 781); 

                auto tg_yyyyyyyz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 782); 

                auto tg_yyyyyyyz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 783); 

                auto tg_yyyyyyyz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 784); 

                auto tg_yyyyyyyz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 785); 

                auto tg_yyyyyyyz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 786); 

                auto tg_yyyyyyyz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 787); 

                auto tg_yyyyyyyz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 788); 

                auto tg_yyyyyyyz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 789); 

                auto tg_yyyyyyyz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 790); 

                auto tg_yyyyyyyz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 791); 

                auto tg_yyyyyyyz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 792); 

                auto tg_yyyyyyyz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 793); 

                auto tg_yyyyyyyz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 794); 

                auto tg_yyyyyyyz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 795); 

                auto tg_yyyyyyyz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 796); 

                auto tg_yyyyyyyz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 797); 

                auto tg_yyyyyyzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 798); 

                auto tg_yyyyyyzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 799); 

                auto tg_yyyyyyzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 800); 

                auto tg_yyyyyyzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 801); 

                auto tg_yyyyyyzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 802); 

                auto tg_yyyyyyzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 803); 

                auto tg_yyyyyyzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 804); 

                auto tg_yyyyyyzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 805); 

                auto tg_yyyyyyzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 806); 

                auto tg_yyyyyyzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 807); 

                auto tg_yyyyyyzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 808); 

                auto tg_yyyyyyzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 809); 

                auto tg_yyyyyyzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 810); 

                auto tg_yyyyyyzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 811); 

                auto tg_yyyyyyzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 812); 

                auto tg_yyyyyyzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 813); 

                auto tg_yyyyyyzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 814); 

                auto tg_yyyyyyzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 815); 

                auto tg_yyyyyyzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 816); 

                auto tg_yyyyyyzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 817); 

                auto tg_yyyyyyzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 818); 

                auto tg_yyyyyzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 819); 

                auto tg_yyyyyzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 820); 

                auto tg_yyyyyzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 821); 

                auto tg_yyyyyzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 822); 

                auto tg_yyyyyzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 823); 

                auto tg_yyyyyzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 824); 

                auto tg_yyyyyzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 825); 

                auto tg_yyyyyzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 826); 

                auto tg_yyyyyzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 827); 

                auto tg_yyyyyzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 828); 

                auto tg_yyyyyzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 829); 

                auto tg_yyyyyzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 830); 

                auto tg_yyyyyzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 831); 

                auto tg_yyyyyzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 832); 

                auto tg_yyyyyzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 833); 

                auto tg_yyyyyzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 834); 

                auto tg_yyyyyzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 835); 

                auto tg_yyyyyzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 836); 

                auto tg_yyyyyzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 837); 

                auto tg_yyyyyzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 838); 

                auto tg_yyyyyzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 839); 

                auto tg_yyyyzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 840); 

                auto tg_yyyyzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 841); 

                auto tg_yyyyzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 842); 

                auto tg_yyyyzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 843); 

                auto tg_yyyyzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 844); 

                auto tg_yyyyzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 845); 

                auto tg_yyyyzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 846); 

                auto tg_yyyyzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 847); 

                auto tg_yyyyzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 848); 

                auto tg_yyyyzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 849); 

                auto tg_yyyyzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 850); 

                // Batch of Integrals (757,851)

                #pragma omp simd aligned(fxn, fza, tg_yyyyyy_xxxxy_0, tg_yyyyyy_xxxxy_1, tg_yyyyyy_xxxxz_0, \
                                         tg_yyyyyy_xxxxz_1, tg_yyyyyy_xxxyy_0, tg_yyyyyy_xxxyy_1, tg_yyyyyy_xxxyz_0, \
                                         tg_yyyyyy_xxxyz_1, tg_yyyyyy_xxxzz_0, tg_yyyyyy_xxxzz_1, tg_yyyyyy_xxyyy_0, \
                                         tg_yyyyyy_xxyyy_1, tg_yyyyyy_xxyyz_0, tg_yyyyyy_xxyyz_1, tg_yyyyyy_xxyzz_0, \
                                         tg_yyyyyy_xxyzz_1, tg_yyyyyy_xxzzz_0, tg_yyyyyy_xxzzz_1, tg_yyyyyy_xyyyy_0, \
                                         tg_yyyyyy_xyyyy_1, tg_yyyyyy_xyyyz_0, tg_yyyyyy_xyyyz_1, tg_yyyyyy_xyyzz_0, \
                                         tg_yyyyyy_xyyzz_1, tg_yyyyyy_xyzzz_0, tg_yyyyyy_xyzzz_1, tg_yyyyyy_xzzzz_0, \
                                         tg_yyyyyy_xzzzz_1, tg_yyyyyy_yyyyy_0, tg_yyyyyy_yyyyy_1, tg_yyyyyy_yyyyz_0, \
                                         tg_yyyyyy_yyyyz_1, tg_yyyyyy_yyyzz_0, tg_yyyyyy_yyyzz_1, tg_yyyyyy_yyzzz_0, \
                                         tg_yyyyyy_yyzzz_1, tg_yyyyyy_yzzzz_0, tg_yyyyyy_yzzzz_1, tg_yyyyyy_zzzzz_0, \
                                         tg_yyyyyy_zzzzz_1, tg_yyyyyyy_xxxx_1, tg_yyyyyyy_xxxxy_0, tg_yyyyyyy_xxxxy_1, \
                                         tg_yyyyyyy_xxxxz_0, tg_yyyyyyy_xxxxz_1, tg_yyyyyyy_xxxy_1, tg_yyyyyyy_xxxyy_0, \
                                         tg_yyyyyyy_xxxyy_1, tg_yyyyyyy_xxxyz_0, tg_yyyyyyy_xxxyz_1, tg_yyyyyyy_xxxz_1, \
                                         tg_yyyyyyy_xxxzz_0, tg_yyyyyyy_xxxzz_1, tg_yyyyyyy_xxyy_1, tg_yyyyyyy_xxyyy_0, \
                                         tg_yyyyyyy_xxyyy_1, tg_yyyyyyy_xxyyz_0, tg_yyyyyyy_xxyyz_1, tg_yyyyyyy_xxyz_1, \
                                         tg_yyyyyyy_xxyzz_0, tg_yyyyyyy_xxyzz_1, tg_yyyyyyy_xxzz_1, tg_yyyyyyy_xxzzz_0, \
                                         tg_yyyyyyy_xxzzz_1, tg_yyyyyyy_xyyy_1, tg_yyyyyyy_xyyyy_0, tg_yyyyyyy_xyyyy_1, \
                                         tg_yyyyyyy_xyyyz_0, tg_yyyyyyy_xyyyz_1, tg_yyyyyyy_xyyz_1, tg_yyyyyyy_xyyzz_0, \
                                         tg_yyyyyyy_xyyzz_1, tg_yyyyyyy_xyzz_1, tg_yyyyyyy_xyzzz_0, tg_yyyyyyy_xyzzz_1, \
                                         tg_yyyyyyy_xzzz_1, tg_yyyyyyy_xzzzz_0, tg_yyyyyyy_xzzzz_1, tg_yyyyyyy_yyyy_1, \
                                         tg_yyyyyyy_yyyyy_0, tg_yyyyyyy_yyyyy_1, tg_yyyyyyy_yyyyz_0, tg_yyyyyyy_yyyyz_1, \
                                         tg_yyyyyyy_yyyz_1, tg_yyyyyyy_yyyzz_0, tg_yyyyyyy_yyyzz_1, tg_yyyyyyy_yyzz_1, \
                                         tg_yyyyyyy_yyzzz_0, tg_yyyyyyy_yyzzz_1, tg_yyyyyyy_yzzz_1, tg_yyyyyyy_yzzzz_0, \
                                         tg_yyyyyyy_yzzzz_1, tg_yyyyyyy_zzzz_1, tg_yyyyyyy_zzzzz_0, tg_yyyyyyy_zzzzz_1, \
                                         tg_yyyyyyyy_xxxxy_0, tg_yyyyyyyy_xxxxz_0, tg_yyyyyyyy_xxxyy_0, tg_yyyyyyyy_xxxyz_0, \
                                         tg_yyyyyyyy_xxxzz_0, tg_yyyyyyyy_xxyyy_0, tg_yyyyyyyy_xxyyz_0, tg_yyyyyyyy_xxyzz_0, \
                                         tg_yyyyyyyy_xxzzz_0, tg_yyyyyyyy_xyyyy_0, tg_yyyyyyyy_xyyyz_0, tg_yyyyyyyy_xyyzz_0, \
                                         tg_yyyyyyyy_xyzzz_0, tg_yyyyyyyy_xzzzz_0, tg_yyyyyyyy_yyyyy_0, tg_yyyyyyyy_yyyyz_0, \
                                         tg_yyyyyyyy_yyyzz_0, tg_yyyyyyyy_yyzzz_0, tg_yyyyyyyy_yzzzz_0, tg_yyyyyyyy_zzzzz_0, \
                                         tg_yyyyyyyz_xxxxx_0, tg_yyyyyyyz_xxxxy_0, tg_yyyyyyyz_xxxxz_0, tg_yyyyyyyz_xxxyy_0, \
                                         tg_yyyyyyyz_xxxyz_0, tg_yyyyyyyz_xxxzz_0, tg_yyyyyyyz_xxyyy_0, tg_yyyyyyyz_xxyyz_0, \
                                         tg_yyyyyyyz_xxyzz_0, tg_yyyyyyyz_xxzzz_0, tg_yyyyyyyz_xyyyy_0, tg_yyyyyyyz_xyyyz_0, \
                                         tg_yyyyyyyz_xyyzz_0, tg_yyyyyyyz_xyzzz_0, tg_yyyyyyyz_xzzzz_0, tg_yyyyyyyz_yyyyy_0, \
                                         tg_yyyyyyyz_yyyyz_0, tg_yyyyyyyz_yyyzz_0, tg_yyyyyyyz_yyzzz_0, tg_yyyyyyyz_yzzzz_0, \
                                         tg_yyyyyyyz_zzzzz_0, tg_yyyyyyz_xxxx_1, tg_yyyyyyz_xxxxx_0, tg_yyyyyyz_xxxxx_1, \
                                         tg_yyyyyyz_xxxxy_0, tg_yyyyyyz_xxxxy_1, tg_yyyyyyz_xxxxz_0, tg_yyyyyyz_xxxxz_1, \
                                         tg_yyyyyyz_xxxy_1, tg_yyyyyyz_xxxyy_0, tg_yyyyyyz_xxxyy_1, tg_yyyyyyz_xxxyz_0, \
                                         tg_yyyyyyz_xxxyz_1, tg_yyyyyyz_xxxz_1, tg_yyyyyyz_xxxzz_0, tg_yyyyyyz_xxxzz_1, \
                                         tg_yyyyyyz_xxyy_1, tg_yyyyyyz_xxyyy_0, tg_yyyyyyz_xxyyy_1, tg_yyyyyyz_xxyyz_0, \
                                         tg_yyyyyyz_xxyyz_1, tg_yyyyyyz_xxyz_1, tg_yyyyyyz_xxyzz_0, tg_yyyyyyz_xxyzz_1, \
                                         tg_yyyyyyz_xxzz_1, tg_yyyyyyz_xxzzz_0, tg_yyyyyyz_xxzzz_1, tg_yyyyyyz_xyyy_1, \
                                         tg_yyyyyyz_xyyyy_0, tg_yyyyyyz_xyyyy_1, tg_yyyyyyz_xyyyz_0, tg_yyyyyyz_xyyyz_1, \
                                         tg_yyyyyyz_xyyz_1, tg_yyyyyyz_xyyzz_0, tg_yyyyyyz_xyyzz_1, tg_yyyyyyz_xyzz_1, \
                                         tg_yyyyyyz_xyzzz_0, tg_yyyyyyz_xyzzz_1, tg_yyyyyyz_xzzz_1, tg_yyyyyyz_xzzzz_0, \
                                         tg_yyyyyyz_xzzzz_1, tg_yyyyyyz_yyyy_1, tg_yyyyyyz_yyyyy_0, tg_yyyyyyz_yyyyy_1, \
                                         tg_yyyyyyz_yyyyz_0, tg_yyyyyyz_yyyyz_1, tg_yyyyyyz_yyyz_1, tg_yyyyyyz_yyyzz_0, \
                                         tg_yyyyyyz_yyyzz_1, tg_yyyyyyz_yyzz_1, tg_yyyyyyz_yyzzz_0, tg_yyyyyyz_yyzzz_1, \
                                         tg_yyyyyyz_yzzz_1, tg_yyyyyyz_yzzzz_0, tg_yyyyyyz_yzzzz_1, tg_yyyyyyz_zzzz_1, \
                                         tg_yyyyyyz_zzzzz_0, tg_yyyyyyz_zzzzz_1, tg_yyyyyyzz_xxxxx_0, tg_yyyyyyzz_xxxxy_0, \
                                         tg_yyyyyyzz_xxxxz_0, tg_yyyyyyzz_xxxyy_0, tg_yyyyyyzz_xxxyz_0, tg_yyyyyyzz_xxxzz_0, \
                                         tg_yyyyyyzz_xxyyy_0, tg_yyyyyyzz_xxyyz_0, tg_yyyyyyzz_xxyzz_0, tg_yyyyyyzz_xxzzz_0, \
                                         tg_yyyyyyzz_xyyyy_0, tg_yyyyyyzz_xyyyz_0, tg_yyyyyyzz_xyyzz_0, tg_yyyyyyzz_xyzzz_0, \
                                         tg_yyyyyyzz_xzzzz_0, tg_yyyyyyzz_yyyyy_0, tg_yyyyyyzz_yyyyz_0, tg_yyyyyyzz_yyyzz_0, \
                                         tg_yyyyyyzz_yyzzz_0, tg_yyyyyyzz_yzzzz_0, tg_yyyyyyzz_zzzzz_0, tg_yyyyyz_xxxxx_0, \
                                         tg_yyyyyz_xxxxx_1, tg_yyyyyz_xxxxy_0, tg_yyyyyz_xxxxy_1, tg_yyyyyz_xxxxz_0, \
                                         tg_yyyyyz_xxxxz_1, tg_yyyyyz_xxxyy_0, tg_yyyyyz_xxxyy_1, tg_yyyyyz_xxxyz_0, \
                                         tg_yyyyyz_xxxyz_1, tg_yyyyyz_xxxzz_0, tg_yyyyyz_xxxzz_1, tg_yyyyyz_xxyyy_0, \
                                         tg_yyyyyz_xxyyy_1, tg_yyyyyz_xxyyz_0, tg_yyyyyz_xxyyz_1, tg_yyyyyz_xxyzz_0, \
                                         tg_yyyyyz_xxyzz_1, tg_yyyyyz_xxzzz_0, tg_yyyyyz_xxzzz_1, tg_yyyyyz_xyyyy_0, \
                                         tg_yyyyyz_xyyyy_1, tg_yyyyyz_xyyyz_0, tg_yyyyyz_xyyyz_1, tg_yyyyyz_xyyzz_0, \
                                         tg_yyyyyz_xyyzz_1, tg_yyyyyz_xyzzz_0, tg_yyyyyz_xyzzz_1, tg_yyyyyz_xzzzz_0, \
                                         tg_yyyyyz_xzzzz_1, tg_yyyyyz_yyyyy_0, tg_yyyyyz_yyyyy_1, tg_yyyyyz_yyyyz_0, \
                                         tg_yyyyyz_yyyyz_1, tg_yyyyyz_yyyzz_0, tg_yyyyyz_yyyzz_1, tg_yyyyyz_yyzzz_0, \
                                         tg_yyyyyz_yyzzz_1, tg_yyyyyz_yzzzz_0, tg_yyyyyz_yzzzz_1, tg_yyyyyz_zzzzz_0, \
                                         tg_yyyyyz_zzzzz_1, tg_yyyyyzz_xxxx_1, tg_yyyyyzz_xxxxx_0, tg_yyyyyzz_xxxxx_1, \
                                         tg_yyyyyzz_xxxxy_0, tg_yyyyyzz_xxxxy_1, tg_yyyyyzz_xxxxz_0, tg_yyyyyzz_xxxxz_1, \
                                         tg_yyyyyzz_xxxy_1, tg_yyyyyzz_xxxyy_0, tg_yyyyyzz_xxxyy_1, tg_yyyyyzz_xxxyz_0, \
                                         tg_yyyyyzz_xxxyz_1, tg_yyyyyzz_xxxz_1, tg_yyyyyzz_xxxzz_0, tg_yyyyyzz_xxxzz_1, \
                                         tg_yyyyyzz_xxyy_1, tg_yyyyyzz_xxyyy_0, tg_yyyyyzz_xxyyy_1, tg_yyyyyzz_xxyyz_0, \
                                         tg_yyyyyzz_xxyyz_1, tg_yyyyyzz_xxyz_1, tg_yyyyyzz_xxyzz_0, tg_yyyyyzz_xxyzz_1, \
                                         tg_yyyyyzz_xxzz_1, tg_yyyyyzz_xxzzz_0, tg_yyyyyzz_xxzzz_1, tg_yyyyyzz_xyyy_1, \
                                         tg_yyyyyzz_xyyyy_0, tg_yyyyyzz_xyyyy_1, tg_yyyyyzz_xyyyz_0, tg_yyyyyzz_xyyyz_1, \
                                         tg_yyyyyzz_xyyz_1, tg_yyyyyzz_xyyzz_0, tg_yyyyyzz_xyyzz_1, tg_yyyyyzz_xyzz_1, \
                                         tg_yyyyyzz_xyzzz_0, tg_yyyyyzz_xyzzz_1, tg_yyyyyzz_xzzz_1, tg_yyyyyzz_xzzzz_0, \
                                         tg_yyyyyzz_xzzzz_1, tg_yyyyyzz_yyyy_1, tg_yyyyyzz_yyyyy_0, tg_yyyyyzz_yyyyy_1, \
                                         tg_yyyyyzz_yyyyz_0, tg_yyyyyzz_yyyyz_1, tg_yyyyyzz_yyyz_1, tg_yyyyyzz_yyyzz_0, \
                                         tg_yyyyyzz_yyyzz_1, tg_yyyyyzz_yyzz_1, tg_yyyyyzz_yyzzz_0, tg_yyyyyzz_yyzzz_1, \
                                         tg_yyyyyzz_yzzz_1, tg_yyyyyzz_yzzzz_0, tg_yyyyyzz_yzzzz_1, tg_yyyyyzz_zzzz_1, \
                                         tg_yyyyyzz_zzzzz_0, tg_yyyyyzz_zzzzz_1, tg_yyyyyzzz_xxxxx_0, tg_yyyyyzzz_xxxxy_0, \
                                         tg_yyyyyzzz_xxxxz_0, tg_yyyyyzzz_xxxyy_0, tg_yyyyyzzz_xxxyz_0, tg_yyyyyzzz_xxxzz_0, \
                                         tg_yyyyyzzz_xxyyy_0, tg_yyyyyzzz_xxyyz_0, tg_yyyyyzzz_xxyzz_0, tg_yyyyyzzz_xxzzz_0, \
                                         tg_yyyyyzzz_xyyyy_0, tg_yyyyyzzz_xyyyz_0, tg_yyyyyzzz_xyyzz_0, tg_yyyyyzzz_xyzzz_0, \
                                         tg_yyyyyzzz_xzzzz_0, tg_yyyyyzzz_yyyyy_0, tg_yyyyyzzz_yyyyz_0, tg_yyyyyzzz_yyyzz_0, \
                                         tg_yyyyyzzz_yyzzz_0, tg_yyyyyzzz_yzzzz_0, tg_yyyyyzzz_zzzzz_0, tg_yyyyzz_xxxxx_0, \
                                         tg_yyyyzz_xxxxx_1, tg_yyyyzz_xxxxy_0, tg_yyyyzz_xxxxy_1, tg_yyyyzz_xxxxz_0, \
                                         tg_yyyyzz_xxxxz_1, tg_yyyyzz_xxxyy_0, tg_yyyyzz_xxxyy_1, tg_yyyyzz_xxxyz_0, \
                                         tg_yyyyzz_xxxyz_1, tg_yyyyzz_xxxzz_0, tg_yyyyzz_xxxzz_1, tg_yyyyzz_xxyyy_0, \
                                         tg_yyyyzz_xxyyy_1, tg_yyyyzz_xxyyz_0, tg_yyyyzz_xxyyz_1, tg_yyyyzz_xxyzz_0, \
                                         tg_yyyyzz_xxyzz_1, tg_yyyyzz_xxzzz_0, tg_yyyyzz_xxzzz_1, tg_yyyyzz_xyyyy_0, \
                                         tg_yyyyzz_xyyyy_1, tg_yyyyzz_xyyyz_0, tg_yyyyzz_xyyyz_1, tg_yyyyzz_xyyzz_0, \
                                         tg_yyyyzz_xyyzz_1, tg_yyyyzz_xyzzz_0, tg_yyyyzz_xyzzz_1, tg_yyyyzz_xzzzz_0, \
                                         tg_yyyyzz_xzzzz_1, tg_yyyyzz_yyyyy_0, tg_yyyyzz_yyyyy_1, tg_yyyyzz_yyyyz_0, \
                                         tg_yyyyzz_yyyyz_1, tg_yyyyzz_yyyzz_0, tg_yyyyzz_yyyzz_1, tg_yyyyzz_yyzzz_0, \
                                         tg_yyyyzz_yyzzz_1, tg_yyyyzz_yzzzz_0, tg_yyyyzz_yzzzz_1, tg_yyyyzz_zzzzz_0, \
                                         tg_yyyyzz_zzzzz_1, tg_yyyyzzz_xxxx_1, tg_yyyyzzz_xxxxx_0, tg_yyyyzzz_xxxxx_1, \
                                         tg_yyyyzzz_xxxxy_0, tg_yyyyzzz_xxxxy_1, tg_yyyyzzz_xxxxz_0, tg_yyyyzzz_xxxxz_1, \
                                         tg_yyyyzzz_xxxy_1, tg_yyyyzzz_xxxyy_0, tg_yyyyzzz_xxxyy_1, tg_yyyyzzz_xxxyz_0, \
                                         tg_yyyyzzz_xxxyz_1, tg_yyyyzzz_xxxz_1, tg_yyyyzzz_xxxzz_0, tg_yyyyzzz_xxxzz_1, \
                                         tg_yyyyzzz_xxyy_1, tg_yyyyzzz_xxyyy_0, tg_yyyyzzz_xxyyy_1, tg_yyyyzzz_xxyyz_0, \
                                         tg_yyyyzzz_xxyyz_1, tg_yyyyzzz_xxyz_1, tg_yyyyzzz_xxyzz_0, tg_yyyyzzz_xxyzz_1, \
                                         tg_yyyyzzz_xxzz_1, tg_yyyyzzz_xxzzz_0, tg_yyyyzzz_xxzzz_1, tg_yyyyzzz_xyyy_1, \
                                         tg_yyyyzzz_xyyyy_0, tg_yyyyzzz_xyyyy_1, tg_yyyyzzz_xyyyz_0, tg_yyyyzzz_xyyyz_1, \
                                         tg_yyyyzzz_xyyz_1, tg_yyyyzzz_xyyzz_0, tg_yyyyzzz_xyyzz_1, tg_yyyyzzz_xyzz_1, \
                                         tg_yyyyzzz_xyzzz_0, tg_yyyyzzz_xyzzz_1, tg_yyyyzzz_xzzz_1, tg_yyyyzzz_xzzzz_0, \
                                         tg_yyyyzzz_xzzzz_1, tg_yyyyzzz_yyyy_1, tg_yyyyzzz_yyyyy_0, tg_yyyyzzz_yyyyy_1, \
                                         tg_yyyyzzz_yyyyz_0, tg_yyyyzzz_yyyyz_1, tg_yyyyzzz_yyyz_1, tg_yyyyzzz_yyyzz_0, \
                                         tg_yyyyzzz_yyyzz_1, tg_yyyyzzz_yyzz_1, tg_yyyyzzz_yyzzz_0, tg_yyyyzzz_yyzzz_1, \
                                         tg_yyyyzzz_yzzz_1, tg_yyyyzzz_yzzzz_0, tg_yyyyzzz_yzzzz_1, tg_yyyyzzz_zzzz_1, \
                                         tg_yyyyzzz_zzzzz_0, tg_yyyyzzz_zzzzz_1, tg_yyyyzzzz_xxxxx_0, tg_yyyyzzzz_xxxxy_0, \
                                         tg_yyyyzzzz_xxxxz_0, tg_yyyyzzzz_xxxyy_0, tg_yyyyzzzz_xxxyz_0, tg_yyyyzzzz_xxxzz_0, \
                                         tg_yyyyzzzz_xxyyy_0, tg_yyyyzzzz_xxyyz_0, tg_yyyyzzzz_xxyzz_0, tg_yyyyzzzz_xxzzz_0, \
                                         tg_yyyyzzzz_xyyyy_0, tg_yyyzzz_xxxxx_0, tg_yyyzzz_xxxxx_1, tg_yyyzzz_xxxxy_0, \
                                         tg_yyyzzz_xxxxy_1, tg_yyyzzz_xxxxz_0, tg_yyyzzz_xxxxz_1, tg_yyyzzz_xxxyy_0, \
                                         tg_yyyzzz_xxxyy_1, tg_yyyzzz_xxxyz_0, tg_yyyzzz_xxxyz_1, tg_yyyzzz_xxxzz_0, \
                                         tg_yyyzzz_xxxzz_1, tg_yyyzzz_xxyyy_0, tg_yyyzzz_xxyyy_1, tg_yyyzzz_xxyyz_0, \
                                         tg_yyyzzz_xxyyz_1, tg_yyyzzz_xxyzz_0, tg_yyyzzz_xxyzz_1, tg_yyyzzz_xxzzz_0, \
                                         tg_yyyzzz_xxzzz_1, tg_yyyzzz_xyyyy_0, tg_yyyzzz_xyyyy_1, tg_yyyzzz_xyyyz_0, \
                                         tg_yyyzzz_xyyyz_1, tg_yyyzzz_xyyzz_0, tg_yyyzzz_xyyzz_1, tg_yyyzzz_xyzzz_0, \
                                         tg_yyyzzz_xyzzz_1, tg_yyyzzz_xzzzz_0, tg_yyyzzz_xzzzz_1, tg_yyyzzz_yyyyy_0, \
                                         tg_yyyzzz_yyyyy_1, tg_yyyzzz_yyyyz_0, tg_yyyzzz_yyyyz_1, tg_yyyzzz_yyyzz_0, \
                                         tg_yyyzzz_yyyzz_1, tg_yyyzzz_yyzzz_0, tg_yyyzzz_yyzzz_1, tg_yyyzzz_yzzzz_0, \
                                         tg_yyyzzz_yzzzz_1, tg_yyyzzz_zzzzz_0, tg_yyyzzz_zzzzz_1, tg_yyyzzzz_xxxx_1, \
                                         tg_yyyzzzz_xxxxx_0, tg_yyyzzzz_xxxxx_1, tg_yyyzzzz_xxxxy_0, tg_yyyzzzz_xxxxy_1, \
                                         tg_yyyzzzz_xxxxz_0, tg_yyyzzzz_xxxxz_1, tg_yyyzzzz_xxxy_1, tg_yyyzzzz_xxxyy_0, \
                                         tg_yyyzzzz_xxxyy_1, tg_yyyzzzz_xxxyz_0, tg_yyyzzzz_xxxyz_1, tg_yyyzzzz_xxxz_1, \
                                         tg_yyyzzzz_xxxzz_0, tg_yyyzzzz_xxxzz_1, tg_yyyzzzz_xxyy_1, tg_yyyzzzz_xxyyy_0, \
                                         tg_yyyzzzz_xxyyy_1, tg_yyyzzzz_xxyyz_0, tg_yyyzzzz_xxyyz_1, tg_yyyzzzz_xxyz_1, \
                                         tg_yyyzzzz_xxyzz_0, tg_yyyzzzz_xxyzz_1, tg_yyyzzzz_xxzz_1, tg_yyyzzzz_xxzzz_0, \
                                         tg_yyyzzzz_xxzzz_1, tg_yyyzzzz_xyyy_1, tg_yyyzzzz_xyyyy_0, tg_yyyzzzz_xyyyy_1, \
                                         tg_yyzzzz_xxxxx_0, tg_yyzzzz_xxxxx_1, tg_yyzzzz_xxxxy_0, tg_yyzzzz_xxxxy_1, \
                                         tg_yyzzzz_xxxxz_0, tg_yyzzzz_xxxxz_1, tg_yyzzzz_xxxyy_0, tg_yyzzzz_xxxyy_1, \
                                         tg_yyzzzz_xxxyz_0, tg_yyzzzz_xxxyz_1, tg_yyzzzz_xxxzz_0, tg_yyzzzz_xxxzz_1, \
                                         tg_yyzzzz_xxyyy_0, tg_yyzzzz_xxyyy_1, tg_yyzzzz_xxyyz_0, tg_yyzzzz_xxyyz_1, \
                                         tg_yyzzzz_xxyzz_0, tg_yyzzzz_xxyzz_1, tg_yyzzzz_xxzzz_0, tg_yyzzzz_xxzzz_1, \
                                         tg_yyzzzz_xyyyy_0, tg_yyzzzz_xyyyy_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyyyyy_xxxxy_0[j] = pb_y * tg_yyyyyyy_xxxxy_0[j] + fr * tg_yyyyyyy_xxxxy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxxxy_0[j] - tg_yyyyyy_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_xxxx_1[j];

                    tg_yyyyyyyy_xxxxz_0[j] = pb_y * tg_yyyyyyy_xxxxz_0[j] + fr * tg_yyyyyyy_xxxxz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxxxz_0[j] - tg_yyyyyy_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyyyy_xxxyy_0[j] = pb_y * tg_yyyyyyy_xxxyy_0[j] + fr * tg_yyyyyyy_xxxyy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxxyy_0[j] - tg_yyyyyy_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyy_xxxy_1[j];

                    tg_yyyyyyyy_xxxyz_0[j] = pb_y * tg_yyyyyyy_xxxyz_0[j] + fr * tg_yyyyyyy_xxxyz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxxyz_0[j] - tg_yyyyyy_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_xxxz_1[j];

                    tg_yyyyyyyy_xxxzz_0[j] = pb_y * tg_yyyyyyy_xxxzz_0[j] + fr * tg_yyyyyyy_xxxzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxxzz_0[j] - tg_yyyyyy_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyyyy_xxyyy_0[j] = pb_y * tg_yyyyyyy_xxyyy_0[j] + fr * tg_yyyyyyy_xxyyy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxyyy_0[j] - tg_yyyyyy_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyy_xxyy_1[j];

                    tg_yyyyyyyy_xxyyz_0[j] = pb_y * tg_yyyyyyy_xxyyz_0[j] + fr * tg_yyyyyyy_xxyyz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxyyz_0[j] - tg_yyyyyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyy_xxyz_1[j];

                    tg_yyyyyyyy_xxyzz_0[j] = pb_y * tg_yyyyyyy_xxyzz_0[j] + fr * tg_yyyyyyy_xxyzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxyzz_0[j] - tg_yyyyyy_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_xxzz_1[j];

                    tg_yyyyyyyy_xxzzz_0[j] = pb_y * tg_yyyyyyy_xxzzz_0[j] + fr * tg_yyyyyyy_xxzzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxzzz_0[j] - tg_yyyyyy_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyyyy_xyyyy_0[j] = pb_y * tg_yyyyyyy_xyyyy_0[j] + fr * tg_yyyyyyy_xyyyy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xyyyy_0[j] - tg_yyyyyy_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyyy_xyyy_1[j];

                    tg_yyyyyyyy_xyyyz_0[j] = pb_y * tg_yyyyyyy_xyyyz_0[j] + fr * tg_yyyyyyy_xyyyz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xyyyz_0[j] - tg_yyyyyy_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyy_xyyz_1[j];

                    tg_yyyyyyyy_xyyzz_0[j] = pb_y * tg_yyyyyyy_xyyzz_0[j] + fr * tg_yyyyyyy_xyyzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xyyzz_0[j] - tg_yyyyyy_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyy_xyzz_1[j];

                    tg_yyyyyyyy_xyzzz_0[j] = pb_y * tg_yyyyyyy_xyzzz_0[j] + fr * tg_yyyyyyy_xyzzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xyzzz_0[j] - tg_yyyyyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_xzzz_1[j];

                    tg_yyyyyyyy_xzzzz_0[j] = pb_y * tg_yyyyyyy_xzzzz_0[j] + fr * tg_yyyyyyy_xzzzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xzzzz_0[j] - tg_yyyyyy_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyyyy_yyyyy_0[j] = pb_y * tg_yyyyyyy_yyyyy_0[j] + fr * tg_yyyyyyy_yyyyy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yyyyy_0[j] - tg_yyyyyy_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyyy_yyyy_1[j];

                    tg_yyyyyyyy_yyyyz_0[j] = pb_y * tg_yyyyyyy_yyyyz_0[j] + fr * tg_yyyyyyy_yyyyz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yyyyz_0[j] - tg_yyyyyy_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyyy_yyyz_1[j];

                    tg_yyyyyyyy_yyyzz_0[j] = pb_y * tg_yyyyyyy_yyyzz_0[j] + fr * tg_yyyyyyy_yyyzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yyyzz_0[j] - tg_yyyyyy_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyy_yyzz_1[j];

                    tg_yyyyyyyy_yyzzz_0[j] = pb_y * tg_yyyyyyy_yyzzz_0[j] + fr * tg_yyyyyyy_yyzzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yyzzz_0[j] - tg_yyyyyy_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyy_yzzz_1[j];

                    tg_yyyyyyyy_yzzzz_0[j] = pb_y * tg_yyyyyyy_yzzzz_0[j] + fr * tg_yyyyyyy_yzzzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yzzzz_0[j] - tg_yyyyyy_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_zzzz_1[j];

                    tg_yyyyyyyy_zzzzz_0[j] = pb_y * tg_yyyyyyy_zzzzz_0[j] + fr * tg_yyyyyyy_zzzzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_zzzzz_0[j] - tg_yyyyyy_zzzzz_1[j] * fl1_fza);

                    tg_yyyyyyyz_xxxxx_0[j] = pb_y * tg_yyyyyyz_xxxxx_0[j] + fr * tg_yyyyyyz_xxxxx_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxxxx_0[j] - tg_yyyyyz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyyyz_xxxxy_0[j] = pb_y * tg_yyyyyyz_xxxxy_0[j] + fr * tg_yyyyyyz_xxxxy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxxxy_0[j] - tg_yyyyyz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_xxxx_1[j];

                    tg_yyyyyyyz_xxxxz_0[j] = pb_y * tg_yyyyyyz_xxxxz_0[j] + fr * tg_yyyyyyz_xxxxz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxxxz_0[j] - tg_yyyyyz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyyyz_xxxyy_0[j] = pb_y * tg_yyyyyyz_xxxyy_0[j] + fr * tg_yyyyyyz_xxxyy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxxyy_0[j] - tg_yyyyyz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyz_xxxy_1[j];

                    tg_yyyyyyyz_xxxyz_0[j] = pb_y * tg_yyyyyyz_xxxyz_0[j] + fr * tg_yyyyyyz_xxxyz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxxyz_0[j] - tg_yyyyyz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_xxxz_1[j];

                    tg_yyyyyyyz_xxxzz_0[j] = pb_y * tg_yyyyyyz_xxxzz_0[j] + fr * tg_yyyyyyz_xxxzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxxzz_0[j] - tg_yyyyyz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyyyz_xxyyy_0[j] = pb_y * tg_yyyyyyz_xxyyy_0[j] + fr * tg_yyyyyyz_xxyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxyyy_0[j] - tg_yyyyyz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyz_xxyy_1[j];

                    tg_yyyyyyyz_xxyyz_0[j] = pb_y * tg_yyyyyyz_xxyyz_0[j] + fr * tg_yyyyyyz_xxyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxyyz_0[j] - tg_yyyyyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyz_xxyz_1[j];

                    tg_yyyyyyyz_xxyzz_0[j] = pb_y * tg_yyyyyyz_xxyzz_0[j] + fr * tg_yyyyyyz_xxyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxyzz_0[j] - tg_yyyyyz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_xxzz_1[j];

                    tg_yyyyyyyz_xxzzz_0[j] = pb_y * tg_yyyyyyz_xxzzz_0[j] + fr * tg_yyyyyyz_xxzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxzzz_0[j] - tg_yyyyyz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyyyz_xyyyy_0[j] = pb_y * tg_yyyyyyz_xyyyy_0[j] + fr * tg_yyyyyyz_xyyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xyyyy_0[j] - tg_yyyyyz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyyz_xyyy_1[j];

                    tg_yyyyyyyz_xyyyz_0[j] = pb_y * tg_yyyyyyz_xyyyz_0[j] + fr * tg_yyyyyyz_xyyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xyyyz_0[j] - tg_yyyyyz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyz_xyyz_1[j];

                    tg_yyyyyyyz_xyyzz_0[j] = pb_y * tg_yyyyyyz_xyyzz_0[j] + fr * tg_yyyyyyz_xyyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xyyzz_0[j] - tg_yyyyyz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyz_xyzz_1[j];

                    tg_yyyyyyyz_xyzzz_0[j] = pb_y * tg_yyyyyyz_xyzzz_0[j] + fr * tg_yyyyyyz_xyzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xyzzz_0[j] - tg_yyyyyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_xzzz_1[j];

                    tg_yyyyyyyz_xzzzz_0[j] = pb_y * tg_yyyyyyz_xzzzz_0[j] + fr * tg_yyyyyyz_xzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xzzzz_0[j] - tg_yyyyyz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyyyz_yyyyy_0[j] = pb_y * tg_yyyyyyz_yyyyy_0[j] + fr * tg_yyyyyyz_yyyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yyyyy_0[j] - tg_yyyyyz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyyz_yyyy_1[j];

                    tg_yyyyyyyz_yyyyz_0[j] = pb_y * tg_yyyyyyz_yyyyz_0[j] + fr * tg_yyyyyyz_yyyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yyyyz_0[j] - tg_yyyyyz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyyz_yyyz_1[j];

                    tg_yyyyyyyz_yyyzz_0[j] = pb_y * tg_yyyyyyz_yyyzz_0[j] + fr * tg_yyyyyyz_yyyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yyyzz_0[j] - tg_yyyyyz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyz_yyzz_1[j];

                    tg_yyyyyyyz_yyzzz_0[j] = pb_y * tg_yyyyyyz_yyzzz_0[j] + fr * tg_yyyyyyz_yyzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yyzzz_0[j] - tg_yyyyyz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyz_yzzz_1[j];

                    tg_yyyyyyyz_yzzzz_0[j] = pb_y * tg_yyyyyyz_yzzzz_0[j] + fr * tg_yyyyyyz_yzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yzzzz_0[j] - tg_yyyyyz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_zzzz_1[j];

                    tg_yyyyyyyz_zzzzz_0[j] = pb_y * tg_yyyyyyz_zzzzz_0[j] + fr * tg_yyyyyyz_zzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_zzzzz_0[j] - tg_yyyyyz_zzzzz_1[j] * fl1_fza);

                    tg_yyyyyyzz_xxxxx_0[j] = pb_y * tg_yyyyyzz_xxxxx_0[j] + fr * tg_yyyyyzz_xxxxx_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxxxx_0[j] - tg_yyyyzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyyzz_xxxxy_0[j] = pb_y * tg_yyyyyzz_xxxxy_0[j] + fr * tg_yyyyyzz_xxxxy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxxxy_0[j] - tg_yyyyzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_xxxx_1[j];

                    tg_yyyyyyzz_xxxxz_0[j] = pb_y * tg_yyyyyzz_xxxxz_0[j] + fr * tg_yyyyyzz_xxxxz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxxxz_0[j] - tg_yyyyzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyyzz_xxxyy_0[j] = pb_y * tg_yyyyyzz_xxxyy_0[j] + fr * tg_yyyyyzz_xxxyy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxxyy_0[j] - tg_yyyyzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyzz_xxxy_1[j];

                    tg_yyyyyyzz_xxxyz_0[j] = pb_y * tg_yyyyyzz_xxxyz_0[j] + fr * tg_yyyyyzz_xxxyz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxxyz_0[j] - tg_yyyyzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_xxxz_1[j];

                    tg_yyyyyyzz_xxxzz_0[j] = pb_y * tg_yyyyyzz_xxxzz_0[j] + fr * tg_yyyyyzz_xxxzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxxzz_0[j] - tg_yyyyzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyyzz_xxyyy_0[j] = pb_y * tg_yyyyyzz_xxyyy_0[j] + fr * tg_yyyyyzz_xxyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxyyy_0[j] - tg_yyyyzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyzz_xxyy_1[j];

                    tg_yyyyyyzz_xxyyz_0[j] = pb_y * tg_yyyyyzz_xxyyz_0[j] + fr * tg_yyyyyzz_xxyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxyyz_0[j] - tg_yyyyzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyzz_xxyz_1[j];

                    tg_yyyyyyzz_xxyzz_0[j] = pb_y * tg_yyyyyzz_xxyzz_0[j] + fr * tg_yyyyyzz_xxyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxyzz_0[j] - tg_yyyyzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_xxzz_1[j];

                    tg_yyyyyyzz_xxzzz_0[j] = pb_y * tg_yyyyyzz_xxzzz_0[j] + fr * tg_yyyyyzz_xxzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxzzz_0[j] - tg_yyyyzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyyzz_xyyyy_0[j] = pb_y * tg_yyyyyzz_xyyyy_0[j] + fr * tg_yyyyyzz_xyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xyyyy_0[j] - tg_yyyyzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyzz_xyyy_1[j];

                    tg_yyyyyyzz_xyyyz_0[j] = pb_y * tg_yyyyyzz_xyyyz_0[j] + fr * tg_yyyyyzz_xyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xyyyz_0[j] - tg_yyyyzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyzz_xyyz_1[j];

                    tg_yyyyyyzz_xyyzz_0[j] = pb_y * tg_yyyyyzz_xyyzz_0[j] + fr * tg_yyyyyzz_xyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xyyzz_0[j] - tg_yyyyzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyzz_xyzz_1[j];

                    tg_yyyyyyzz_xyzzz_0[j] = pb_y * tg_yyyyyzz_xyzzz_0[j] + fr * tg_yyyyyzz_xyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xyzzz_0[j] - tg_yyyyzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_xzzz_1[j];

                    tg_yyyyyyzz_xzzzz_0[j] = pb_y * tg_yyyyyzz_xzzzz_0[j] + fr * tg_yyyyyzz_xzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xzzzz_0[j] - tg_yyyyzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyyzz_yyyyy_0[j] = pb_y * tg_yyyyyzz_yyyyy_0[j] + fr * tg_yyyyyzz_yyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yyyyy_0[j] - tg_yyyyzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyzz_yyyy_1[j];

                    tg_yyyyyyzz_yyyyz_0[j] = pb_y * tg_yyyyyzz_yyyyz_0[j] + fr * tg_yyyyyzz_yyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yyyyz_0[j] - tg_yyyyzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyzz_yyyz_1[j];

                    tg_yyyyyyzz_yyyzz_0[j] = pb_y * tg_yyyyyzz_yyyzz_0[j] + fr * tg_yyyyyzz_yyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yyyzz_0[j] - tg_yyyyzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyzz_yyzz_1[j];

                    tg_yyyyyyzz_yyzzz_0[j] = pb_y * tg_yyyyyzz_yyzzz_0[j] + fr * tg_yyyyyzz_yyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yyzzz_0[j] - tg_yyyyzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyzz_yzzz_1[j];

                    tg_yyyyyyzz_yzzzz_0[j] = pb_y * tg_yyyyyzz_yzzzz_0[j] + fr * tg_yyyyyzz_yzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yzzzz_0[j] - tg_yyyyzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_zzzz_1[j];

                    tg_yyyyyyzz_zzzzz_0[j] = pb_y * tg_yyyyyzz_zzzzz_0[j] + fr * tg_yyyyyzz_zzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_zzzzz_0[j] - tg_yyyyzz_zzzzz_1[j] * fl1_fza);

                    tg_yyyyyzzz_xxxxx_0[j] = pb_y * tg_yyyyzzz_xxxxx_0[j] + fr * tg_yyyyzzz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxxxx_0[j] - tg_yyyzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyyzzz_xxxxy_0[j] = pb_y * tg_yyyyzzz_xxxxy_0[j] + fr * tg_yyyyzzz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxxxy_0[j] - tg_yyyzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_xxxx_1[j];

                    tg_yyyyyzzz_xxxxz_0[j] = pb_y * tg_yyyyzzz_xxxxz_0[j] + fr * tg_yyyyzzz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxxxz_0[j] - tg_yyyzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyyzzz_xxxyy_0[j] = pb_y * tg_yyyyzzz_xxxyy_0[j] + fr * tg_yyyyzzz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxxyy_0[j] - tg_yyyzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzzz_xxxy_1[j];

                    tg_yyyyyzzz_xxxyz_0[j] = pb_y * tg_yyyyzzz_xxxyz_0[j] + fr * tg_yyyyzzz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxxyz_0[j] - tg_yyyzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_xxxz_1[j];

                    tg_yyyyyzzz_xxxzz_0[j] = pb_y * tg_yyyyzzz_xxxzz_0[j] + fr * tg_yyyyzzz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxxzz_0[j] - tg_yyyzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyyzzz_xxyyy_0[j] = pb_y * tg_yyyyzzz_xxyyy_0[j] + fr * tg_yyyyzzz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxyyy_0[j] - tg_yyyzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzzz_xxyy_1[j];

                    tg_yyyyyzzz_xxyyz_0[j] = pb_y * tg_yyyyzzz_xxyyz_0[j] + fr * tg_yyyyzzz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxyyz_0[j] - tg_yyyzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzzz_xxyz_1[j];

                    tg_yyyyyzzz_xxyzz_0[j] = pb_y * tg_yyyyzzz_xxyzz_0[j] + fr * tg_yyyyzzz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxyzz_0[j] - tg_yyyzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_xxzz_1[j];

                    tg_yyyyyzzz_xxzzz_0[j] = pb_y * tg_yyyyzzz_xxzzz_0[j] + fr * tg_yyyyzzz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxzzz_0[j] - tg_yyyzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyyzzz_xyyyy_0[j] = pb_y * tg_yyyyzzz_xyyyy_0[j] + fr * tg_yyyyzzz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xyyyy_0[j] - tg_yyyzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzzz_xyyy_1[j];

                    tg_yyyyyzzz_xyyyz_0[j] = pb_y * tg_yyyyzzz_xyyyz_0[j] + fr * tg_yyyyzzz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xyyyz_0[j] - tg_yyyzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzzz_xyyz_1[j];

                    tg_yyyyyzzz_xyyzz_0[j] = pb_y * tg_yyyyzzz_xyyzz_0[j] + fr * tg_yyyyzzz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xyyzz_0[j] - tg_yyyzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzzz_xyzz_1[j];

                    tg_yyyyyzzz_xyzzz_0[j] = pb_y * tg_yyyyzzz_xyzzz_0[j] + fr * tg_yyyyzzz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xyzzz_0[j] - tg_yyyzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_xzzz_1[j];

                    tg_yyyyyzzz_xzzzz_0[j] = pb_y * tg_yyyyzzz_xzzzz_0[j] + fr * tg_yyyyzzz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xzzzz_0[j] - tg_yyyzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyyzzz_yyyyy_0[j] = pb_y * tg_yyyyzzz_yyyyy_0[j] + fr * tg_yyyyzzz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yyyyy_0[j] - tg_yyyzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyzzz_yyyy_1[j];

                    tg_yyyyyzzz_yyyyz_0[j] = pb_y * tg_yyyyzzz_yyyyz_0[j] + fr * tg_yyyyzzz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yyyyz_0[j] - tg_yyyzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzzz_yyyz_1[j];

                    tg_yyyyyzzz_yyyzz_0[j] = pb_y * tg_yyyyzzz_yyyzz_0[j] + fr * tg_yyyyzzz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yyyzz_0[j] - tg_yyyzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzzz_yyzz_1[j];

                    tg_yyyyyzzz_yyzzz_0[j] = pb_y * tg_yyyyzzz_yyzzz_0[j] + fr * tg_yyyyzzz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yyzzz_0[j] - tg_yyyzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzzz_yzzz_1[j];

                    tg_yyyyyzzz_yzzzz_0[j] = pb_y * tg_yyyyzzz_yzzzz_0[j] + fr * tg_yyyyzzz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yzzzz_0[j] - tg_yyyzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_zzzz_1[j];

                    tg_yyyyyzzz_zzzzz_0[j] = pb_y * tg_yyyyzzz_zzzzz_0[j] + fr * tg_yyyyzzz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_zzzzz_0[j] - tg_yyyzzz_zzzzz_1[j] * fl1_fza);

                    tg_yyyyzzzz_xxxxx_0[j] = pb_y * tg_yyyzzzz_xxxxx_0[j] + fr * tg_yyyzzzz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxxxx_0[j] - tg_yyzzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyzzzz_xxxxy_0[j] = pb_y * tg_yyyzzzz_xxxxy_0[j] + fr * tg_yyyzzzz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxxxy_0[j] - tg_yyzzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_xxxx_1[j];

                    tg_yyyyzzzz_xxxxz_0[j] = pb_y * tg_yyyzzzz_xxxxz_0[j] + fr * tg_yyyzzzz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxxxz_0[j] - tg_yyzzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyzzzz_xxxyy_0[j] = pb_y * tg_yyyzzzz_xxxyy_0[j] + fr * tg_yyyzzzz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxxyy_0[j] - tg_yyzzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzzz_xxxy_1[j];

                    tg_yyyyzzzz_xxxyz_0[j] = pb_y * tg_yyyzzzz_xxxyz_0[j] + fr * tg_yyyzzzz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxxyz_0[j] - tg_yyzzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_xxxz_1[j];

                    tg_yyyyzzzz_xxxzz_0[j] = pb_y * tg_yyyzzzz_xxxzz_0[j] + fr * tg_yyyzzzz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxxzz_0[j] - tg_yyzzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyzzzz_xxyyy_0[j] = pb_y * tg_yyyzzzz_xxyyy_0[j] + fr * tg_yyyzzzz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxyyy_0[j] - tg_yyzzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzzz_xxyy_1[j];

                    tg_yyyyzzzz_xxyyz_0[j] = pb_y * tg_yyyzzzz_xxyyz_0[j] + fr * tg_yyyzzzz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxyyz_0[j] - tg_yyzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzzz_xxyz_1[j];

                    tg_yyyyzzzz_xxyzz_0[j] = pb_y * tg_yyyzzzz_xxyzz_0[j] + fr * tg_yyyzzzz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxyzz_0[j] - tg_yyzzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_xxzz_1[j];

                    tg_yyyyzzzz_xxzzz_0[j] = pb_y * tg_yyyzzzz_xxzzz_0[j] + fr * tg_yyyzzzz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxzzz_0[j] - tg_yyzzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyzzzz_xyyyy_0[j] = pb_y * tg_yyyzzzz_xyyyy_0[j] + fr * tg_yyyzzzz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xyyyy_0[j] - tg_yyzzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzzz_xyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSH_851_945(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyyzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 683); 

                auto tg_yyyzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 684); 

                auto tg_yyyzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 685); 

                auto tg_yyyzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 686); 

                auto tg_yyyzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 687); 

                auto tg_yyyzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 688); 

                auto tg_yyyzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 689); 

                auto tg_yyyzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 690); 

                auto tg_yyyzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 691); 

                auto tg_yyyzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 692); 

                auto tg_yyzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 693); 

                auto tg_yyzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 694); 

                auto tg_yyzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 695); 

                auto tg_yyzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 696); 

                auto tg_yyzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 697); 

                auto tg_yyzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 698); 

                auto tg_yyzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 699); 

                auto tg_yyzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 700); 

                auto tg_yyzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 701); 

                auto tg_yyzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 702); 

                auto tg_yyzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 703); 

                auto tg_yyzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 704); 

                auto tg_yyzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 705); 

                auto tg_yyzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 706); 

                auto tg_yyzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 707); 

                auto tg_yyzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 708); 

                auto tg_yyzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 709); 

                auto tg_yyzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 710); 

                auto tg_yyzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 711); 

                auto tg_yyzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 712); 

                auto tg_yyzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 713); 

                auto tg_yzzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 714); 

                auto tg_yzzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 715); 

                auto tg_yzzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 716); 

                auto tg_yzzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 717); 

                auto tg_yzzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 718); 

                auto tg_yzzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 719); 

                auto tg_yzzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 720); 

                auto tg_yzzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 721); 

                auto tg_yzzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 722); 

                auto tg_yzzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 723); 

                auto tg_yzzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 724); 

                auto tg_yzzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 725); 

                auto tg_yzzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 726); 

                auto tg_yzzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 727); 

                auto tg_yzzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 728); 

                auto tg_yzzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 729); 

                auto tg_yzzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 730); 

                auto tg_yzzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 731); 

                auto tg_yzzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 732); 

                auto tg_yzzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 733); 

                auto tg_yzzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 734); 

                auto tg_zzzzzzz_xxxxx_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 735); 

                auto tg_zzzzzzz_xxxxy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 736); 

                auto tg_zzzzzzz_xxxxz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 737); 

                auto tg_zzzzzzz_xxxyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 738); 

                auto tg_zzzzzzz_xxxyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 739); 

                auto tg_zzzzzzz_xxxzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 740); 

                auto tg_zzzzzzz_xxyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 741); 

                auto tg_zzzzzzz_xxyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 742); 

                auto tg_zzzzzzz_xxyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 743); 

                auto tg_zzzzzzz_xxzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 744); 

                auto tg_zzzzzzz_xyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 745); 

                auto tg_zzzzzzz_xyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 746); 

                auto tg_zzzzzzz_xyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 747); 

                auto tg_zzzzzzz_xyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 748); 

                auto tg_zzzzzzz_xzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 749); 

                auto tg_zzzzzzz_yyyyy_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 750); 

                auto tg_zzzzzzz_yyyyz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 751); 

                auto tg_zzzzzzz_yyyzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 752); 

                auto tg_zzzzzzz_yyzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 753); 

                auto tg_zzzzzzz_yzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 754); 

                auto tg_zzzzzzz_zzzzz_1 = primBuffer[pidx_g_7_5_m1].data(756 * idx + 755); 

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

                auto tg_yyyzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 487); 

                auto tg_yyyzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 488); 

                auto tg_yyyzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 489); 

                auto tg_yyyzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 490); 

                auto tg_yyyzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 491); 

                auto tg_yyyzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 492); 

                auto tg_yyyzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 493); 

                auto tg_yyyzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 494); 

                auto tg_yyzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 495); 

                auto tg_yyzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 496); 

                auto tg_yyzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 497); 

                auto tg_yyzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 498); 

                auto tg_yyzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 499); 

                auto tg_yyzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 500); 

                auto tg_yyzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 501); 

                auto tg_yyzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 502); 

                auto tg_yyzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 503); 

                auto tg_yyzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 504); 

                auto tg_yyzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 505); 

                auto tg_yyzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 506); 

                auto tg_yyzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 507); 

                auto tg_yyzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 508); 

                auto tg_yyzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 509); 

                auto tg_yzzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 510); 

                auto tg_yzzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 511); 

                auto tg_yzzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 512); 

                auto tg_yzzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 513); 

                auto tg_yzzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 514); 

                auto tg_yzzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 515); 

                auto tg_yzzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 516); 

                auto tg_yzzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 517); 

                auto tg_yzzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 518); 

                auto tg_yzzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 519); 

                auto tg_yzzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 520); 

                auto tg_yzzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 521); 

                auto tg_yzzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 522); 

                auto tg_yzzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 523); 

                auto tg_yzzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 524); 

                auto tg_zzzzzzz_xxxx_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 525); 

                auto tg_zzzzzzz_xxxy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 526); 

                auto tg_zzzzzzz_xxxz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 527); 

                auto tg_zzzzzzz_xxyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 528); 

                auto tg_zzzzzzz_xxyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 529); 

                auto tg_zzzzzzz_xxzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 530); 

                auto tg_zzzzzzz_xyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 531); 

                auto tg_zzzzzzz_xyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 532); 

                auto tg_zzzzzzz_xyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 533); 

                auto tg_zzzzzzz_xzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 534); 

                auto tg_zzzzzzz_yyyy_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 535); 

                auto tg_zzzzzzz_yyyz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 536); 

                auto tg_zzzzzzz_yyzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 537); 

                auto tg_zzzzzzz_yzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 538); 

                auto tg_zzzzzzz_zzzz_1 = primBuffer[pidx_g_7_4_m1].data(540 * idx + 539); 

                // set up pointers to integrals

                auto tg_yyyyzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 851); 

                auto tg_yyyyzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 852); 

                auto tg_yyyyzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 853); 

                auto tg_yyyyzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 854); 

                auto tg_yyyyzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 855); 

                auto tg_yyyyzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 856); 

                auto tg_yyyyzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 857); 

                auto tg_yyyyzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 858); 

                auto tg_yyyyzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 859); 

                auto tg_yyyyzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 860); 

                auto tg_yyyzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 861); 

                auto tg_yyyzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 862); 

                auto tg_yyyzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 863); 

                auto tg_yyyzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 864); 

                auto tg_yyyzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 865); 

                auto tg_yyyzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 866); 

                auto tg_yyyzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 867); 

                auto tg_yyyzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 868); 

                auto tg_yyyzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 869); 

                auto tg_yyyzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 870); 

                auto tg_yyyzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 871); 

                auto tg_yyyzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 872); 

                auto tg_yyyzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 873); 

                auto tg_yyyzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 874); 

                auto tg_yyyzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 875); 

                auto tg_yyyzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 876); 

                auto tg_yyyzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 877); 

                auto tg_yyyzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 878); 

                auto tg_yyyzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 879); 

                auto tg_yyyzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 880); 

                auto tg_yyyzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 881); 

                auto tg_yyzzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 882); 

                auto tg_yyzzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 883); 

                auto tg_yyzzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 884); 

                auto tg_yyzzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 885); 

                auto tg_yyzzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 886); 

                auto tg_yyzzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 887); 

                auto tg_yyzzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 888); 

                auto tg_yyzzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 889); 

                auto tg_yyzzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 890); 

                auto tg_yyzzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 891); 

                auto tg_yyzzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 892); 

                auto tg_yyzzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 893); 

                auto tg_yyzzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 894); 

                auto tg_yyzzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 895); 

                auto tg_yyzzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 896); 

                auto tg_yyzzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 897); 

                auto tg_yyzzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 898); 

                auto tg_yyzzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 899); 

                auto tg_yyzzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 900); 

                auto tg_yyzzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 901); 

                auto tg_yyzzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 902); 

                auto tg_yzzzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 903); 

                auto tg_yzzzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 904); 

                auto tg_yzzzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 905); 

                auto tg_yzzzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 906); 

                auto tg_yzzzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 907); 

                auto tg_yzzzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 908); 

                auto tg_yzzzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 909); 

                auto tg_yzzzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 910); 

                auto tg_yzzzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 911); 

                auto tg_yzzzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 912); 

                auto tg_yzzzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 913); 

                auto tg_yzzzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 914); 

                auto tg_yzzzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 915); 

                auto tg_yzzzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 916); 

                auto tg_yzzzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 917); 

                auto tg_yzzzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 918); 

                auto tg_yzzzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 919); 

                auto tg_yzzzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 920); 

                auto tg_yzzzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 921); 

                auto tg_yzzzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 922); 

                auto tg_yzzzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 923); 

                auto tg_zzzzzzzz_xxxxx_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 924); 

                auto tg_zzzzzzzz_xxxxy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 925); 

                auto tg_zzzzzzzz_xxxxz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 926); 

                auto tg_zzzzzzzz_xxxyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 927); 

                auto tg_zzzzzzzz_xxxyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 928); 

                auto tg_zzzzzzzz_xxxzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 929); 

                auto tg_zzzzzzzz_xxyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 930); 

                auto tg_zzzzzzzz_xxyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 931); 

                auto tg_zzzzzzzz_xxyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 932); 

                auto tg_zzzzzzzz_xxzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 933); 

                auto tg_zzzzzzzz_xyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 934); 

                auto tg_zzzzzzzz_xyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 935); 

                auto tg_zzzzzzzz_xyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 936); 

                auto tg_zzzzzzzz_xyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 937); 

                auto tg_zzzzzzzz_xzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 938); 

                auto tg_zzzzzzzz_yyyyy_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 939); 

                auto tg_zzzzzzzz_yyyyz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 940); 

                auto tg_zzzzzzzz_yyyzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 941); 

                auto tg_zzzzzzzz_yyzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 942); 

                auto tg_zzzzzzzz_yzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 943); 

                auto tg_zzzzzzzz_zzzzz_0 = primBuffer[pidx_g_8_5_m0].data(945 * idx + 944); 

                // Batch of Integrals (851,945)

                #pragma omp simd aligned(fxn, fza, tg_yyyyzzzz_xyyyz_0, tg_yyyyzzzz_xyyzz_0, \
                                         tg_yyyyzzzz_xyzzz_0, tg_yyyyzzzz_xzzzz_0, tg_yyyyzzzz_yyyyy_0, tg_yyyyzzzz_yyyyz_0, \
                                         tg_yyyyzzzz_yyyzz_0, tg_yyyyzzzz_yyzzz_0, tg_yyyyzzzz_yzzzz_0, tg_yyyyzzzz_zzzzz_0, \
                                         tg_yyyzzzz_xyyyz_0, tg_yyyzzzz_xyyyz_1, tg_yyyzzzz_xyyz_1, tg_yyyzzzz_xyyzz_0, \
                                         tg_yyyzzzz_xyyzz_1, tg_yyyzzzz_xyzz_1, tg_yyyzzzz_xyzzz_0, tg_yyyzzzz_xyzzz_1, \
                                         tg_yyyzzzz_xzzz_1, tg_yyyzzzz_xzzzz_0, tg_yyyzzzz_xzzzz_1, tg_yyyzzzz_yyyy_1, \
                                         tg_yyyzzzz_yyyyy_0, tg_yyyzzzz_yyyyy_1, tg_yyyzzzz_yyyyz_0, tg_yyyzzzz_yyyyz_1, \
                                         tg_yyyzzzz_yyyz_1, tg_yyyzzzz_yyyzz_0, tg_yyyzzzz_yyyzz_1, tg_yyyzzzz_yyzz_1, \
                                         tg_yyyzzzz_yyzzz_0, tg_yyyzzzz_yyzzz_1, tg_yyyzzzz_yzzz_1, tg_yyyzzzz_yzzzz_0, \
                                         tg_yyyzzzz_yzzzz_1, tg_yyyzzzz_zzzz_1, tg_yyyzzzz_zzzzz_0, tg_yyyzzzz_zzzzz_1, \
                                         tg_yyyzzzzz_xxxxx_0, tg_yyyzzzzz_xxxxy_0, tg_yyyzzzzz_xxxxz_0, tg_yyyzzzzz_xxxyy_0, \
                                         tg_yyyzzzzz_xxxyz_0, tg_yyyzzzzz_xxxzz_0, tg_yyyzzzzz_xxyyy_0, tg_yyyzzzzz_xxyyz_0, \
                                         tg_yyyzzzzz_xxyzz_0, tg_yyyzzzzz_xxzzz_0, tg_yyyzzzzz_xyyyy_0, tg_yyyzzzzz_xyyyz_0, \
                                         tg_yyyzzzzz_xyyzz_0, tg_yyyzzzzz_xyzzz_0, tg_yyyzzzzz_xzzzz_0, tg_yyyzzzzz_yyyyy_0, \
                                         tg_yyyzzzzz_yyyyz_0, tg_yyyzzzzz_yyyzz_0, tg_yyyzzzzz_yyzzz_0, tg_yyyzzzzz_yzzzz_0, \
                                         tg_yyyzzzzz_zzzzz_0, tg_yyzzzz_xyyyz_0, tg_yyzzzz_xyyyz_1, tg_yyzzzz_xyyzz_0, \
                                         tg_yyzzzz_xyyzz_1, tg_yyzzzz_xyzzz_0, tg_yyzzzz_xyzzz_1, tg_yyzzzz_xzzzz_0, \
                                         tg_yyzzzz_xzzzz_1, tg_yyzzzz_yyyyy_0, tg_yyzzzz_yyyyy_1, tg_yyzzzz_yyyyz_0, \
                                         tg_yyzzzz_yyyyz_1, tg_yyzzzz_yyyzz_0, tg_yyzzzz_yyyzz_1, tg_yyzzzz_yyzzz_0, \
                                         tg_yyzzzz_yyzzz_1, tg_yyzzzz_yzzzz_0, tg_yyzzzz_yzzzz_1, tg_yyzzzz_zzzzz_0, \
                                         tg_yyzzzz_zzzzz_1, tg_yyzzzzz_xxxx_1, tg_yyzzzzz_xxxxx_0, tg_yyzzzzz_xxxxx_1, \
                                         tg_yyzzzzz_xxxxy_0, tg_yyzzzzz_xxxxy_1, tg_yyzzzzz_xxxxz_0, tg_yyzzzzz_xxxxz_1, \
                                         tg_yyzzzzz_xxxy_1, tg_yyzzzzz_xxxyy_0, tg_yyzzzzz_xxxyy_1, tg_yyzzzzz_xxxyz_0, \
                                         tg_yyzzzzz_xxxyz_1, tg_yyzzzzz_xxxz_1, tg_yyzzzzz_xxxzz_0, tg_yyzzzzz_xxxzz_1, \
                                         tg_yyzzzzz_xxyy_1, tg_yyzzzzz_xxyyy_0, tg_yyzzzzz_xxyyy_1, tg_yyzzzzz_xxyyz_0, \
                                         tg_yyzzzzz_xxyyz_1, tg_yyzzzzz_xxyz_1, tg_yyzzzzz_xxyzz_0, tg_yyzzzzz_xxyzz_1, \
                                         tg_yyzzzzz_xxzz_1, tg_yyzzzzz_xxzzz_0, tg_yyzzzzz_xxzzz_1, tg_yyzzzzz_xyyy_1, \
                                         tg_yyzzzzz_xyyyy_0, tg_yyzzzzz_xyyyy_1, tg_yyzzzzz_xyyyz_0, tg_yyzzzzz_xyyyz_1, \
                                         tg_yyzzzzz_xyyz_1, tg_yyzzzzz_xyyzz_0, tg_yyzzzzz_xyyzz_1, tg_yyzzzzz_xyzz_1, \
                                         tg_yyzzzzz_xyzzz_0, tg_yyzzzzz_xyzzz_1, tg_yyzzzzz_xzzz_1, tg_yyzzzzz_xzzzz_0, \
                                         tg_yyzzzzz_xzzzz_1, tg_yyzzzzz_yyyy_1, tg_yyzzzzz_yyyyy_0, tg_yyzzzzz_yyyyy_1, \
                                         tg_yyzzzzz_yyyyz_0, tg_yyzzzzz_yyyyz_1, tg_yyzzzzz_yyyz_1, tg_yyzzzzz_yyyzz_0, \
                                         tg_yyzzzzz_yyyzz_1, tg_yyzzzzz_yyzz_1, tg_yyzzzzz_yyzzz_0, tg_yyzzzzz_yyzzz_1, \
                                         tg_yyzzzzz_yzzz_1, tg_yyzzzzz_yzzzz_0, tg_yyzzzzz_yzzzz_1, tg_yyzzzzz_zzzz_1, \
                                         tg_yyzzzzz_zzzzz_0, tg_yyzzzzz_zzzzz_1, tg_yyzzzzzz_xxxxx_0, tg_yyzzzzzz_xxxxy_0, \
                                         tg_yyzzzzzz_xxxxz_0, tg_yyzzzzzz_xxxyy_0, tg_yyzzzzzz_xxxyz_0, tg_yyzzzzzz_xxxzz_0, \
                                         tg_yyzzzzzz_xxyyy_0, tg_yyzzzzzz_xxyyz_0, tg_yyzzzzzz_xxyzz_0, tg_yyzzzzzz_xxzzz_0, \
                                         tg_yyzzzzzz_xyyyy_0, tg_yyzzzzzz_xyyyz_0, tg_yyzzzzzz_xyyzz_0, tg_yyzzzzzz_xyzzz_0, \
                                         tg_yyzzzzzz_xzzzz_0, tg_yyzzzzzz_yyyyy_0, tg_yyzzzzzz_yyyyz_0, tg_yyzzzzzz_yyyzz_0, \
                                         tg_yyzzzzzz_yyzzz_0, tg_yyzzzzzz_yzzzz_0, tg_yyzzzzzz_zzzzz_0, tg_yzzzzz_xxxxx_0, \
                                         tg_yzzzzz_xxxxx_1, tg_yzzzzz_xxxxy_0, tg_yzzzzz_xxxxy_1, tg_yzzzzz_xxxxz_0, \
                                         tg_yzzzzz_xxxxz_1, tg_yzzzzz_xxxyy_0, tg_yzzzzz_xxxyy_1, tg_yzzzzz_xxxyz_0, \
                                         tg_yzzzzz_xxxyz_1, tg_yzzzzz_xxxzz_0, tg_yzzzzz_xxxzz_1, tg_yzzzzz_xxyyy_0, \
                                         tg_yzzzzz_xxyyy_1, tg_yzzzzz_xxyyz_0, tg_yzzzzz_xxyyz_1, tg_yzzzzz_xxyzz_0, \
                                         tg_yzzzzz_xxyzz_1, tg_yzzzzz_xxzzz_0, tg_yzzzzz_xxzzz_1, tg_yzzzzz_xyyyy_0, \
                                         tg_yzzzzz_xyyyy_1, tg_yzzzzz_xyyyz_0, tg_yzzzzz_xyyyz_1, tg_yzzzzz_xyyzz_0, \
                                         tg_yzzzzz_xyyzz_1, tg_yzzzzz_xyzzz_0, tg_yzzzzz_xyzzz_1, tg_yzzzzz_xzzzz_0, \
                                         tg_yzzzzz_xzzzz_1, tg_yzzzzz_yyyyy_0, tg_yzzzzz_yyyyy_1, tg_yzzzzz_yyyyz_0, \
                                         tg_yzzzzz_yyyyz_1, tg_yzzzzz_yyyzz_0, tg_yzzzzz_yyyzz_1, tg_yzzzzz_yyzzz_0, \
                                         tg_yzzzzz_yyzzz_1, tg_yzzzzz_yzzzz_0, tg_yzzzzz_yzzzz_1, tg_yzzzzz_zzzzz_0, \
                                         tg_yzzzzz_zzzzz_1, tg_yzzzzzz_xxxx_1, tg_yzzzzzz_xxxxx_0, tg_yzzzzzz_xxxxx_1, \
                                         tg_yzzzzzz_xxxxy_0, tg_yzzzzzz_xxxxy_1, tg_yzzzzzz_xxxxz_0, tg_yzzzzzz_xxxxz_1, \
                                         tg_yzzzzzz_xxxy_1, tg_yzzzzzz_xxxyy_0, tg_yzzzzzz_xxxyy_1, tg_yzzzzzz_xxxyz_0, \
                                         tg_yzzzzzz_xxxyz_1, tg_yzzzzzz_xxxz_1, tg_yzzzzzz_xxxzz_0, tg_yzzzzzz_xxxzz_1, \
                                         tg_yzzzzzz_xxyy_1, tg_yzzzzzz_xxyyy_0, tg_yzzzzzz_xxyyy_1, tg_yzzzzzz_xxyyz_0, \
                                         tg_yzzzzzz_xxyyz_1, tg_yzzzzzz_xxyz_1, tg_yzzzzzz_xxyzz_0, tg_yzzzzzz_xxyzz_1, \
                                         tg_yzzzzzz_xxzz_1, tg_yzzzzzz_xxzzz_0, tg_yzzzzzz_xxzzz_1, tg_yzzzzzz_xyyy_1, \
                                         tg_yzzzzzz_xyyyy_0, tg_yzzzzzz_xyyyy_1, tg_yzzzzzz_xyyyz_0, tg_yzzzzzz_xyyyz_1, \
                                         tg_yzzzzzz_xyyz_1, tg_yzzzzzz_xyyzz_0, tg_yzzzzzz_xyyzz_1, tg_yzzzzzz_xyzz_1, \
                                         tg_yzzzzzz_xyzzz_0, tg_yzzzzzz_xyzzz_1, tg_yzzzzzz_xzzz_1, tg_yzzzzzz_xzzzz_0, \
                                         tg_yzzzzzz_xzzzz_1, tg_yzzzzzz_yyyy_1, tg_yzzzzzz_yyyyy_0, tg_yzzzzzz_yyyyy_1, \
                                         tg_yzzzzzz_yyyyz_0, tg_yzzzzzz_yyyyz_1, tg_yzzzzzz_yyyz_1, tg_yzzzzzz_yyyzz_0, \
                                         tg_yzzzzzz_yyyzz_1, tg_yzzzzzz_yyzz_1, tg_yzzzzzz_yyzzz_0, tg_yzzzzzz_yyzzz_1, \
                                         tg_yzzzzzz_yzzz_1, tg_yzzzzzz_yzzzz_0, tg_yzzzzzz_yzzzz_1, tg_yzzzzzz_zzzz_1, \
                                         tg_yzzzzzz_zzzzz_0, tg_yzzzzzz_zzzzz_1, tg_yzzzzzzz_xxxxx_0, tg_yzzzzzzz_xxxxy_0, \
                                         tg_yzzzzzzz_xxxxz_0, tg_yzzzzzzz_xxxyy_0, tg_yzzzzzzz_xxxyz_0, tg_yzzzzzzz_xxxzz_0, \
                                         tg_yzzzzzzz_xxyyy_0, tg_yzzzzzzz_xxyyz_0, tg_yzzzzzzz_xxyzz_0, tg_yzzzzzzz_xxzzz_0, \
                                         tg_yzzzzzzz_xyyyy_0, tg_yzzzzzzz_xyyyz_0, tg_yzzzzzzz_xyyzz_0, tg_yzzzzzzz_xyzzz_0, \
                                         tg_yzzzzzzz_xzzzz_0, tg_yzzzzzzz_yyyyy_0, tg_yzzzzzzz_yyyyz_0, tg_yzzzzzzz_yyyzz_0, \
                                         tg_yzzzzzzz_yyzzz_0, tg_yzzzzzzz_yzzzz_0, tg_yzzzzzzz_zzzzz_0, tg_zzzzzz_xxxxx_0, \
                                         tg_zzzzzz_xxxxx_1, tg_zzzzzz_xxxxy_0, tg_zzzzzz_xxxxy_1, tg_zzzzzz_xxxxz_0, \
                                         tg_zzzzzz_xxxxz_1, tg_zzzzzz_xxxyy_0, tg_zzzzzz_xxxyy_1, tg_zzzzzz_xxxyz_0, \
                                         tg_zzzzzz_xxxyz_1, tg_zzzzzz_xxxzz_0, tg_zzzzzz_xxxzz_1, tg_zzzzzz_xxyyy_0, \
                                         tg_zzzzzz_xxyyy_1, tg_zzzzzz_xxyyz_0, tg_zzzzzz_xxyyz_1, tg_zzzzzz_xxyzz_0, \
                                         tg_zzzzzz_xxyzz_1, tg_zzzzzz_xxzzz_0, tg_zzzzzz_xxzzz_1, tg_zzzzzz_xyyyy_0, \
                                         tg_zzzzzz_xyyyy_1, tg_zzzzzz_xyyyz_0, tg_zzzzzz_xyyyz_1, tg_zzzzzz_xyyzz_0, \
                                         tg_zzzzzz_xyyzz_1, tg_zzzzzz_xyzzz_0, tg_zzzzzz_xyzzz_1, tg_zzzzzz_xzzzz_0, \
                                         tg_zzzzzz_xzzzz_1, tg_zzzzzz_yyyyy_0, tg_zzzzzz_yyyyy_1, tg_zzzzzz_yyyyz_0, \
                                         tg_zzzzzz_yyyyz_1, tg_zzzzzz_yyyzz_0, tg_zzzzzz_yyyzz_1, tg_zzzzzz_yyzzz_0, \
                                         tg_zzzzzz_yyzzz_1, tg_zzzzzz_yzzzz_0, tg_zzzzzz_yzzzz_1, tg_zzzzzz_zzzzz_0, \
                                         tg_zzzzzz_zzzzz_1, tg_zzzzzzz_xxxx_1, tg_zzzzzzz_xxxxx_0, tg_zzzzzzz_xxxxx_1, \
                                         tg_zzzzzzz_xxxxy_0, tg_zzzzzzz_xxxxy_1, tg_zzzzzzz_xxxxz_0, tg_zzzzzzz_xxxxz_1, \
                                         tg_zzzzzzz_xxxy_1, tg_zzzzzzz_xxxyy_0, tg_zzzzzzz_xxxyy_1, tg_zzzzzzz_xxxyz_0, \
                                         tg_zzzzzzz_xxxyz_1, tg_zzzzzzz_xxxz_1, tg_zzzzzzz_xxxzz_0, tg_zzzzzzz_xxxzz_1, \
                                         tg_zzzzzzz_xxyy_1, tg_zzzzzzz_xxyyy_0, tg_zzzzzzz_xxyyy_1, tg_zzzzzzz_xxyyz_0, \
                                         tg_zzzzzzz_xxyyz_1, tg_zzzzzzz_xxyz_1, tg_zzzzzzz_xxyzz_0, tg_zzzzzzz_xxyzz_1, \
                                         tg_zzzzzzz_xxzz_1, tg_zzzzzzz_xxzzz_0, tg_zzzzzzz_xxzzz_1, tg_zzzzzzz_xyyy_1, \
                                         tg_zzzzzzz_xyyyy_0, tg_zzzzzzz_xyyyy_1, tg_zzzzzzz_xyyyz_0, tg_zzzzzzz_xyyyz_1, \
                                         tg_zzzzzzz_xyyz_1, tg_zzzzzzz_xyyzz_0, tg_zzzzzzz_xyyzz_1, tg_zzzzzzz_xyzz_1, \
                                         tg_zzzzzzz_xyzzz_0, tg_zzzzzzz_xyzzz_1, tg_zzzzzzz_xzzz_1, tg_zzzzzzz_xzzzz_0, \
                                         tg_zzzzzzz_xzzzz_1, tg_zzzzzzz_yyyy_1, tg_zzzzzzz_yyyyy_0, tg_zzzzzzz_yyyyy_1, \
                                         tg_zzzzzzz_yyyyz_0, tg_zzzzzzz_yyyyz_1, tg_zzzzzzz_yyyz_1, tg_zzzzzzz_yyyzz_0, \
                                         tg_zzzzzzz_yyyzz_1, tg_zzzzzzz_yyzz_1, tg_zzzzzzz_yyzzz_0, tg_zzzzzzz_yyzzz_1, \
                                         tg_zzzzzzz_yzzz_1, tg_zzzzzzz_yzzzz_0, tg_zzzzzzz_yzzzz_1, tg_zzzzzzz_zzzz_1, \
                                         tg_zzzzzzz_zzzzz_0, tg_zzzzzzz_zzzzz_1, tg_zzzzzzzz_xxxxx_0, tg_zzzzzzzz_xxxxy_0, \
                                         tg_zzzzzzzz_xxxxz_0, tg_zzzzzzzz_xxxyy_0, tg_zzzzzzzz_xxxyz_0, tg_zzzzzzzz_xxxzz_0, \
                                         tg_zzzzzzzz_xxyyy_0, tg_zzzzzzzz_xxyyz_0, tg_zzzzzzzz_xxyzz_0, tg_zzzzzzzz_xxzzz_0, \
                                         tg_zzzzzzzz_xyyyy_0, tg_zzzzzzzz_xyyyz_0, tg_zzzzzzzz_xyyzz_0, tg_zzzzzzzz_xyzzz_0, \
                                         tg_zzzzzzzz_xzzzz_0, tg_zzzzzzzz_yyyyy_0, tg_zzzzzzzz_yyyyz_0, tg_zzzzzzzz_yyyzz_0, \
                                         tg_zzzzzzzz_yyzzz_0, tg_zzzzzzzz_yzzzz_0, tg_zzzzzzzz_zzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyzzzz_xyyyz_0[j] = pb_y * tg_yyyzzzz_xyyyz_0[j] + fr * tg_yyyzzzz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xyyyz_0[j] - tg_yyzzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzzz_xyyz_1[j];

                    tg_yyyyzzzz_xyyzz_0[j] = pb_y * tg_yyyzzzz_xyyzz_0[j] + fr * tg_yyyzzzz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xyyzz_0[j] - tg_yyzzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzzz_xyzz_1[j];

                    tg_yyyyzzzz_xyzzz_0[j] = pb_y * tg_yyyzzzz_xyzzz_0[j] + fr * tg_yyyzzzz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xyzzz_0[j] - tg_yyzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_xzzz_1[j];

                    tg_yyyyzzzz_xzzzz_0[j] = pb_y * tg_yyyzzzz_xzzzz_0[j] + fr * tg_yyyzzzz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xzzzz_0[j] - tg_yyzzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyzzzz_yyyyy_0[j] = pb_y * tg_yyyzzzz_yyyyy_0[j] + fr * tg_yyyzzzz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yyyyy_0[j] - tg_yyzzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyzzzz_yyyy_1[j];

                    tg_yyyyzzzz_yyyyz_0[j] = pb_y * tg_yyyzzzz_yyyyz_0[j] + fr * tg_yyyzzzz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yyyyz_0[j] - tg_yyzzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzzz_yyyz_1[j];

                    tg_yyyyzzzz_yyyzz_0[j] = pb_y * tg_yyyzzzz_yyyzz_0[j] + fr * tg_yyyzzzz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yyyzz_0[j] - tg_yyzzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzzz_yyzz_1[j];

                    tg_yyyyzzzz_yyzzz_0[j] = pb_y * tg_yyyzzzz_yyzzz_0[j] + fr * tg_yyyzzzz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yyzzz_0[j] - tg_yyzzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzzz_yzzz_1[j];

                    tg_yyyyzzzz_yzzzz_0[j] = pb_y * tg_yyyzzzz_yzzzz_0[j] + fr * tg_yyyzzzz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yzzzz_0[j] - tg_yyzzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_zzzz_1[j];

                    tg_yyyyzzzz_zzzzz_0[j] = pb_y * tg_yyyzzzz_zzzzz_0[j] + fr * tg_yyyzzzz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_zzzzz_0[j] - tg_yyzzzz_zzzzz_1[j] * fl1_fza);

                    tg_yyyzzzzz_xxxxx_0[j] = pb_y * tg_yyzzzzz_xxxxx_0[j] + fr * tg_yyzzzzz_xxxxx_1[j] + fl1_fx * (tg_yzzzzz_xxxxx_0[j] - tg_yzzzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyzzzzz_xxxxy_0[j] = pb_y * tg_yyzzzzz_xxxxy_0[j] + fr * tg_yyzzzzz_xxxxy_1[j] + fl1_fx * (tg_yzzzzz_xxxxy_0[j] - tg_yzzzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_xxxx_1[j];

                    tg_yyyzzzzz_xxxxz_0[j] = pb_y * tg_yyzzzzz_xxxxz_0[j] + fr * tg_yyzzzzz_xxxxz_1[j] + fl1_fx * (tg_yzzzzz_xxxxz_0[j] - tg_yzzzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyzzzzz_xxxyy_0[j] = pb_y * tg_yyzzzzz_xxxyy_0[j] + fr * tg_yyzzzzz_xxxyy_1[j] + fl1_fx * (tg_yzzzzz_xxxyy_0[j] - tg_yzzzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzzz_xxxy_1[j];

                    tg_yyyzzzzz_xxxyz_0[j] = pb_y * tg_yyzzzzz_xxxyz_0[j] + fr * tg_yyzzzzz_xxxyz_1[j] + fl1_fx * (tg_yzzzzz_xxxyz_0[j] - tg_yzzzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_xxxz_1[j];

                    tg_yyyzzzzz_xxxzz_0[j] = pb_y * tg_yyzzzzz_xxxzz_0[j] + fr * tg_yyzzzzz_xxxzz_1[j] + fl1_fx * (tg_yzzzzz_xxxzz_0[j] - tg_yzzzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyzzzzz_xxyyy_0[j] = pb_y * tg_yyzzzzz_xxyyy_0[j] + fr * tg_yyzzzzz_xxyyy_1[j] + fl1_fx * (tg_yzzzzz_xxyyy_0[j] - tg_yzzzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzzz_xxyy_1[j];

                    tg_yyyzzzzz_xxyyz_0[j] = pb_y * tg_yyzzzzz_xxyyz_0[j] + fr * tg_yyzzzzz_xxyyz_1[j] + fl1_fx * (tg_yzzzzz_xxyyz_0[j] - tg_yzzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzzz_xxyz_1[j];

                    tg_yyyzzzzz_xxyzz_0[j] = pb_y * tg_yyzzzzz_xxyzz_0[j] + fr * tg_yyzzzzz_xxyzz_1[j] + fl1_fx * (tg_yzzzzz_xxyzz_0[j] - tg_yzzzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_xxzz_1[j];

                    tg_yyyzzzzz_xxzzz_0[j] = pb_y * tg_yyzzzzz_xxzzz_0[j] + fr * tg_yyzzzzz_xxzzz_1[j] + fl1_fx * (tg_yzzzzz_xxzzz_0[j] - tg_yzzzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyzzzzz_xyyyy_0[j] = pb_y * tg_yyzzzzz_xyyyy_0[j] + fr * tg_yyzzzzz_xyyyy_1[j] + fl1_fx * (tg_yzzzzz_xyyyy_0[j] - tg_yzzzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzzz_xyyy_1[j];

                    tg_yyyzzzzz_xyyyz_0[j] = pb_y * tg_yyzzzzz_xyyyz_0[j] + fr * tg_yyzzzzz_xyyyz_1[j] + fl1_fx * (tg_yzzzzz_xyyyz_0[j] - tg_yzzzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzzz_xyyz_1[j];

                    tg_yyyzzzzz_xyyzz_0[j] = pb_y * tg_yyzzzzz_xyyzz_0[j] + fr * tg_yyzzzzz_xyyzz_1[j] + fl1_fx * (tg_yzzzzz_xyyzz_0[j] - tg_yzzzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzzz_xyzz_1[j];

                    tg_yyyzzzzz_xyzzz_0[j] = pb_y * tg_yyzzzzz_xyzzz_0[j] + fr * tg_yyzzzzz_xyzzz_1[j] + fl1_fx * (tg_yzzzzz_xyzzz_0[j] - tg_yzzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_xzzz_1[j];

                    tg_yyyzzzzz_xzzzz_0[j] = pb_y * tg_yyzzzzz_xzzzz_0[j] + fr * tg_yyzzzzz_xzzzz_1[j] + fl1_fx * (tg_yzzzzz_xzzzz_0[j] - tg_yzzzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyzzzzz_yyyyy_0[j] = pb_y * tg_yyzzzzz_yyyyy_0[j] + fr * tg_yyzzzzz_yyyyy_1[j] + fl1_fx * (tg_yzzzzz_yyyyy_0[j] - tg_yzzzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzzzzz_yyyy_1[j];

                    tg_yyyzzzzz_yyyyz_0[j] = pb_y * tg_yyzzzzz_yyyyz_0[j] + fr * tg_yyzzzzz_yyyyz_1[j] + fl1_fx * (tg_yzzzzz_yyyyz_0[j] - tg_yzzzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzzz_yyyz_1[j];

                    tg_yyyzzzzz_yyyzz_0[j] = pb_y * tg_yyzzzzz_yyyzz_0[j] + fr * tg_yyzzzzz_yyyzz_1[j] + fl1_fx * (tg_yzzzzz_yyyzz_0[j] - tg_yzzzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzzz_yyzz_1[j];

                    tg_yyyzzzzz_yyzzz_0[j] = pb_y * tg_yyzzzzz_yyzzz_0[j] + fr * tg_yyzzzzz_yyzzz_1[j] + fl1_fx * (tg_yzzzzz_yyzzz_0[j] - tg_yzzzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzzz_yzzz_1[j];

                    tg_yyyzzzzz_yzzzz_0[j] = pb_y * tg_yyzzzzz_yzzzz_0[j] + fr * tg_yyzzzzz_yzzzz_1[j] + fl1_fx * (tg_yzzzzz_yzzzz_0[j] - tg_yzzzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_zzzz_1[j];

                    tg_yyyzzzzz_zzzzz_0[j] = pb_y * tg_yyzzzzz_zzzzz_0[j] + fr * tg_yyzzzzz_zzzzz_1[j] + fl1_fx * (tg_yzzzzz_zzzzz_0[j] - tg_yzzzzz_zzzzz_1[j] * fl1_fza);

                    tg_yyzzzzzz_xxxxx_0[j] = pb_y * tg_yzzzzzz_xxxxx_0[j] + fr * tg_yzzzzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxxx_0[j] - tg_zzzzzz_xxxxx_1[j] * fl1_fza);

                    tg_yyzzzzzz_xxxxy_0[j] = pb_y * tg_yzzzzzz_xxxxy_0[j] + fr * tg_yzzzzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxxy_0[j] - tg_zzzzzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_xxxx_1[j];

                    tg_yyzzzzzz_xxxxz_0[j] = pb_y * tg_yzzzzzz_xxxxz_0[j] + fr * tg_yzzzzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxxz_0[j] - tg_zzzzzz_xxxxz_1[j] * fl1_fza);

                    tg_yyzzzzzz_xxxyy_0[j] = pb_y * tg_yzzzzzz_xxxyy_0[j] + fr * tg_yzzzzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxyy_0[j] - tg_zzzzzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzzz_xxxy_1[j];

                    tg_yyzzzzzz_xxxyz_0[j] = pb_y * tg_yzzzzzz_xxxyz_0[j] + fr * tg_yzzzzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxyz_0[j] - tg_zzzzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_xxxz_1[j];

                    tg_yyzzzzzz_xxxzz_0[j] = pb_y * tg_yzzzzzz_xxxzz_0[j] + fr * tg_yzzzzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxxzz_0[j] - tg_zzzzzz_xxxzz_1[j] * fl1_fza);

                    tg_yyzzzzzz_xxyyy_0[j] = pb_y * tg_yzzzzzz_xxyyy_0[j] + fr * tg_yzzzzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxyyy_0[j] - tg_zzzzzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzzz_xxyy_1[j];

                    tg_yyzzzzzz_xxyyz_0[j] = pb_y * tg_yzzzzzz_xxyyz_0[j] + fr * tg_yzzzzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxyyz_0[j] - tg_zzzzzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzzz_xxyz_1[j];

                    tg_yyzzzzzz_xxyzz_0[j] = pb_y * tg_yzzzzzz_xxyzz_0[j] + fr * tg_yzzzzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxyzz_0[j] - tg_zzzzzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_xxzz_1[j];

                    tg_yyzzzzzz_xxzzz_0[j] = pb_y * tg_yzzzzzz_xxzzz_0[j] + fr * tg_yzzzzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxzzz_0[j] - tg_zzzzzz_xxzzz_1[j] * fl1_fza);

                    tg_yyzzzzzz_xyyyy_0[j] = pb_y * tg_yzzzzzz_xyyyy_0[j] + fr * tg_yzzzzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyyyy_0[j] - tg_zzzzzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzzz_xyyy_1[j];

                    tg_yyzzzzzz_xyyyz_0[j] = pb_y * tg_yzzzzzz_xyyyz_0[j] + fr * tg_yzzzzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyyyz_0[j] - tg_zzzzzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzzz_xyyz_1[j];

                    tg_yyzzzzzz_xyyzz_0[j] = pb_y * tg_yzzzzzz_xyyzz_0[j] + fr * tg_yzzzzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyyzz_0[j] - tg_zzzzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzzz_xyzz_1[j];

                    tg_yyzzzzzz_xyzzz_0[j] = pb_y * tg_yzzzzzz_xyzzz_0[j] + fr * tg_yzzzzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyzzz_0[j] - tg_zzzzzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_xzzz_1[j];

                    tg_yyzzzzzz_xzzzz_0[j] = pb_y * tg_yzzzzzz_xzzzz_0[j] + fr * tg_yzzzzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xzzzz_0[j] - tg_zzzzzz_xzzzz_1[j] * fl1_fza);

                    tg_yyzzzzzz_yyyyy_0[j] = pb_y * tg_yzzzzzz_yyyyy_0[j] + fr * tg_yzzzzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyyyy_0[j] - tg_zzzzzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzzzzz_yyyy_1[j];

                    tg_yyzzzzzz_yyyyz_0[j] = pb_y * tg_yzzzzzz_yyyyz_0[j] + fr * tg_yzzzzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyyyz_0[j] - tg_zzzzzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzzz_yyyz_1[j];

                    tg_yyzzzzzz_yyyzz_0[j] = pb_y * tg_yzzzzzz_yyyzz_0[j] + fr * tg_yzzzzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyyzz_0[j] - tg_zzzzzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzzz_yyzz_1[j];

                    tg_yyzzzzzz_yyzzz_0[j] = pb_y * tg_yzzzzzz_yyzzz_0[j] + fr * tg_yzzzzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyzzz_0[j] - tg_zzzzzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzzz_yzzz_1[j];

                    tg_yyzzzzzz_yzzzz_0[j] = pb_y * tg_yzzzzzz_yzzzz_0[j] + fr * tg_yzzzzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yzzzz_0[j] - tg_zzzzzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_zzzz_1[j];

                    tg_yyzzzzzz_zzzzz_0[j] = pb_y * tg_yzzzzzz_zzzzz_0[j] + fr * tg_yzzzzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_zzzzz_0[j] - tg_zzzzzz_zzzzz_1[j] * fl1_fza);

                    tg_yzzzzzzz_xxxxx_0[j] = pb_y * tg_zzzzzzz_xxxxx_0[j] + fr * tg_zzzzzzz_xxxxx_1[j];

                    tg_yzzzzzzz_xxxxy_0[j] = pb_y * tg_zzzzzzz_xxxxy_0[j] + fr * tg_zzzzzzz_xxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xxxx_1[j];

                    tg_yzzzzzzz_xxxxz_0[j] = pb_y * tg_zzzzzzz_xxxxz_0[j] + fr * tg_zzzzzzz_xxxxz_1[j];

                    tg_yzzzzzzz_xxxyy_0[j] = pb_y * tg_zzzzzzz_xxxyy_0[j] + fr * tg_zzzzzzz_xxxyy_1[j] + fl1_fxn * tg_zzzzzzz_xxxy_1[j];

                    tg_yzzzzzzz_xxxyz_0[j] = pb_y * tg_zzzzzzz_xxxyz_0[j] + fr * tg_zzzzzzz_xxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xxxz_1[j];

                    tg_yzzzzzzz_xxxzz_0[j] = pb_y * tg_zzzzzzz_xxxzz_0[j] + fr * tg_zzzzzzz_xxxzz_1[j];

                    tg_yzzzzzzz_xxyyy_0[j] = pb_y * tg_zzzzzzz_xxyyy_0[j] + fr * tg_zzzzzzz_xxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xxyy_1[j];

                    tg_yzzzzzzz_xxyyz_0[j] = pb_y * tg_zzzzzzz_xxyyz_0[j] + fr * tg_zzzzzzz_xxyyz_1[j] + fl1_fxn * tg_zzzzzzz_xxyz_1[j];

                    tg_yzzzzzzz_xxyzz_0[j] = pb_y * tg_zzzzzzz_xxyzz_0[j] + fr * tg_zzzzzzz_xxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xxzz_1[j];

                    tg_yzzzzzzz_xxzzz_0[j] = pb_y * tg_zzzzzzz_xxzzz_0[j] + fr * tg_zzzzzzz_xxzzz_1[j];

                    tg_yzzzzzzz_xyyyy_0[j] = pb_y * tg_zzzzzzz_xyyyy_0[j] + fr * tg_zzzzzzz_xyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzzzz_xyyy_1[j];

                    tg_yzzzzzzz_xyyyz_0[j] = pb_y * tg_zzzzzzz_xyyyz_0[j] + fr * tg_zzzzzzz_xyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xyyz_1[j];

                    tg_yzzzzzzz_xyyzz_0[j] = pb_y * tg_zzzzzzz_xyyzz_0[j] + fr * tg_zzzzzzz_xyyzz_1[j] + fl1_fxn * tg_zzzzzzz_xyzz_1[j];

                    tg_yzzzzzzz_xyzzz_0[j] = pb_y * tg_zzzzzzz_xyzzz_0[j] + fr * tg_zzzzzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xzzz_1[j];

                    tg_yzzzzzzz_xzzzz_0[j] = pb_y * tg_zzzzzzz_xzzzz_0[j] + fr * tg_zzzzzzz_xzzzz_1[j];

                    tg_yzzzzzzz_yyyyy_0[j] = pb_y * tg_zzzzzzz_yyyyy_0[j] + fr * tg_zzzzzzz_yyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzzzzz_yyyy_1[j];

                    tg_yzzzzzzz_yyyyz_0[j] = pb_y * tg_zzzzzzz_yyyyz_0[j] + fr * tg_zzzzzzz_yyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzzzzz_yyyz_1[j];

                    tg_yzzzzzzz_yyyzz_0[j] = pb_y * tg_zzzzzzz_yyyzz_0[j] + fr * tg_zzzzzzz_yyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_yyzz_1[j];

                    tg_yzzzzzzz_yyzzz_0[j] = pb_y * tg_zzzzzzz_yyzzz_0[j] + fr * tg_zzzzzzz_yyzzz_1[j] + fl1_fxn * tg_zzzzzzz_yzzz_1[j];

                    tg_yzzzzzzz_yzzzz_0[j] = pb_y * tg_zzzzzzz_yzzzz_0[j] + fr * tg_zzzzzzz_yzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_zzzz_1[j];

                    tg_yzzzzzzz_zzzzz_0[j] = pb_y * tg_zzzzzzz_zzzzz_0[j] + fr * tg_zzzzzzz_zzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzzzz_xxxxx_0[j] = pb_z * tg_zzzzzzz_xxxxx_0[j] + fr * tg_zzzzzzz_xxxxx_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxxxx_0[j] - tg_zzzzzz_xxxxx_1[j] * fl1_fza);

                    tg_zzzzzzzz_xxxxy_0[j] = pb_z * tg_zzzzzzz_xxxxy_0[j] + fr * tg_zzzzzzz_xxxxy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxxxy_0[j] - tg_zzzzzz_xxxxy_1[j] * fl1_fza);

                    tg_zzzzzzzz_xxxxz_0[j] = pb_z * tg_zzzzzzz_xxxxz_0[j] + fr * tg_zzzzzzz_xxxxz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxxxz_0[j] - tg_zzzzzz_xxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_xxxx_1[j];

                    tg_zzzzzzzz_xxxyy_0[j] = pb_z * tg_zzzzzzz_xxxyy_0[j] + fr * tg_zzzzzzz_xxxyy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxxyy_0[j] - tg_zzzzzz_xxxyy_1[j] * fl1_fza);

                    tg_zzzzzzzz_xxxyz_0[j] = pb_z * tg_zzzzzzz_xxxyz_0[j] + fr * tg_zzzzzzz_xxxyz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxxyz_0[j] - tg_zzzzzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_xxxy_1[j];

                    tg_zzzzzzzz_xxxzz_0[j] = pb_z * tg_zzzzzzz_xxxzz_0[j] + fr * tg_zzzzzzz_xxxzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxxzz_0[j] - tg_zzzzzz_xxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzzz_xxxz_1[j];

                    tg_zzzzzzzz_xxyyy_0[j] = pb_z * tg_zzzzzzz_xxyyy_0[j] + fr * tg_zzzzzzz_xxyyy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxyyy_0[j] - tg_zzzzzz_xxyyy_1[j] * fl1_fza);

                    tg_zzzzzzzz_xxyyz_0[j] = pb_z * tg_zzzzzzz_xxyyz_0[j] + fr * tg_zzzzzzz_xxyyz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxyyz_0[j] - tg_zzzzzz_xxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_xxyy_1[j];

                    tg_zzzzzzzz_xxyzz_0[j] = pb_z * tg_zzzzzzz_xxyzz_0[j] + fr * tg_zzzzzzz_xxyzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxyzz_0[j] - tg_zzzzzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzzz_xxyz_1[j];

                    tg_zzzzzzzz_xxzzz_0[j] = pb_z * tg_zzzzzzz_xxzzz_0[j] + fr * tg_zzzzzzz_xxzzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxzzz_0[j] - tg_zzzzzz_xxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzzz_xxzz_1[j];

                    tg_zzzzzzzz_xyyyy_0[j] = pb_z * tg_zzzzzzz_xyyyy_0[j] + fr * tg_zzzzzzz_xyyyy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xyyyy_0[j] - tg_zzzzzz_xyyyy_1[j] * fl1_fza);

                    tg_zzzzzzzz_xyyyz_0[j] = pb_z * tg_zzzzzzz_xyyyz_0[j] + fr * tg_zzzzzzz_xyyyz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xyyyz_0[j] - tg_zzzzzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_xyyy_1[j];

                    tg_zzzzzzzz_xyyzz_0[j] = pb_z * tg_zzzzzzz_xyyzz_0[j] + fr * tg_zzzzzzz_xyyzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xyyzz_0[j] - tg_zzzzzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzzz_xyyz_1[j];

                    tg_zzzzzzzz_xyzzz_0[j] = pb_z * tg_zzzzzzz_xyzzz_0[j] + fr * tg_zzzzzzz_xyzzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xyzzz_0[j] - tg_zzzzzz_xyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzzz_xyzz_1[j];

                    tg_zzzzzzzz_xzzzz_0[j] = pb_z * tg_zzzzzzz_xzzzz_0[j] + fr * tg_zzzzzzz_xzzzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xzzzz_0[j] - tg_zzzzzz_xzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzzz_xzzz_1[j];

                    tg_zzzzzzzz_yyyyy_0[j] = pb_z * tg_zzzzzzz_yyyyy_0[j] + fr * tg_zzzzzzz_yyyyy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yyyyy_0[j] - tg_zzzzzz_yyyyy_1[j] * fl1_fza);

                    tg_zzzzzzzz_yyyyz_0[j] = pb_z * tg_zzzzzzz_yyyyz_0[j] + fr * tg_zzzzzzz_yyyyz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yyyyz_0[j] - tg_zzzzzz_yyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_yyyy_1[j];

                    tg_zzzzzzzz_yyyzz_0[j] = pb_z * tg_zzzzzzz_yyyzz_0[j] + fr * tg_zzzzzzz_yyyzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yyyzz_0[j] - tg_zzzzzz_yyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzzz_yyyz_1[j];

                    tg_zzzzzzzz_yyzzz_0[j] = pb_z * tg_zzzzzzz_yyzzz_0[j] + fr * tg_zzzzzzz_yyzzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yyzzz_0[j] - tg_zzzzzz_yyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzzz_yyzz_1[j];

                    tg_zzzzzzzz_yzzzz_0[j] = pb_z * tg_zzzzzzz_yzzzz_0[j] + fr * tg_zzzzzzz_yzzzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yzzzz_0[j] - tg_zzzzzz_yzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzzz_yzzz_1[j];

                    tg_zzzzzzzz_zzzzz_0[j] = pb_z * tg_zzzzzzz_zzzzz_0[j] + fr * tg_zzzzzzz_zzzzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_zzzzz_0[j] - tg_zzzzzz_zzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzzzzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

