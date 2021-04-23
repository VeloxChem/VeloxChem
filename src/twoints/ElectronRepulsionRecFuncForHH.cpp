//
//                           VELOXCHEM 1.0-RC2
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

#include "ElectronRepulsionRecFuncForHH.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSHSH(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSHSH_0_89(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSHSH_89_177(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSHSH_177_265(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSH_265_353(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSH_353_441(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSHSH_0_89(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,89)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxxx_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx); 

                auto tg_xxxx_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 1); 

                auto tg_xxxx_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 2); 

                auto tg_xxxx_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 3); 

                auto tg_xxxx_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 4); 

                auto tg_xxxx_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 5); 

                auto tg_xxxx_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 6); 

                auto tg_xxxx_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 7); 

                auto tg_xxxx_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 8); 

                auto tg_xxxx_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 9); 

                auto tg_xxxx_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 10); 

                auto tg_xxxx_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 11); 

                auto tg_xxxx_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 12); 

                auto tg_xxxx_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 13); 

                auto tg_xxxx_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 14); 

                auto tg_xxxy_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 15); 

                auto tg_xxxy_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 16); 

                auto tg_xxxy_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 17); 

                auto tg_xxxy_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 18); 

                auto tg_xxxy_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 19); 

                auto tg_xxxy_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 20); 

                auto tg_xxxy_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 21); 

                auto tg_xxxy_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 22); 

                auto tg_xxxy_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 23); 

                auto tg_xxxy_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 24); 

                auto tg_xxxy_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 25); 

                auto tg_xxxy_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 26); 

                auto tg_xxxy_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 27); 

                auto tg_xxxy_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 28); 

                auto tg_xxxy_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 29); 

                auto tg_xxxz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 30); 

                auto tg_xxxz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 31); 

                auto tg_xxxz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 32); 

                auto tg_xxxz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 33); 

                auto tg_xxxz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 34); 

                auto tg_xxxz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 35); 

                auto tg_xxxz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 36); 

                auto tg_xxxz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 37); 

                auto tg_xxxz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 38); 

                auto tg_xxxz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 39); 

                auto tg_xxxz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 40); 

                auto tg_xxxz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 41); 

                auto tg_xxxz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 42); 

                auto tg_xxxz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 43); 

                auto tg_xxxz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 44); 

                auto tg_xxyy_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 45); 

                auto tg_xxyy_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 46); 

                auto tg_xxyy_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 47); 

                auto tg_xxyy_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 48); 

                auto tg_xxyy_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 49); 

                auto tg_xxyy_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 50); 

                auto tg_xxyy_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 51); 

                auto tg_xxyy_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 52); 

                auto tg_xxyy_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 53); 

                auto tg_xxyy_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 54); 

                auto tg_xxyy_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 55); 

                auto tg_xxyy_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 56); 

                auto tg_xxyy_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 57); 

                auto tg_xxyy_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 58); 

                auto tg_xxyy_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 59); 

                auto tg_xxyz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 60); 

                auto tg_xxyz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 61); 

                auto tg_xxyz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 62); 

                auto tg_xxyz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 63); 

                auto tg_xxyz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 64); 

                // set up pointers to integrals

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

                // Batch of Integrals (0,89)

                #pragma omp simd aligned(fxn, fza, tg_xxx_xxxxx_0, tg_xxx_xxxxx_1, tg_xxx_xxxxy_0, \
                                         tg_xxx_xxxxy_1, tg_xxx_xxxxz_0, tg_xxx_xxxxz_1, tg_xxx_xxxyy_0, tg_xxx_xxxyy_1, \
                                         tg_xxx_xxxyz_0, tg_xxx_xxxyz_1, tg_xxx_xxxzz_0, tg_xxx_xxxzz_1, tg_xxx_xxyyy_0, \
                                         tg_xxx_xxyyy_1, tg_xxx_xxyyz_0, tg_xxx_xxyyz_1, tg_xxx_xxyzz_0, tg_xxx_xxyzz_1, \
                                         tg_xxx_xxzzz_0, tg_xxx_xxzzz_1, tg_xxx_xyyyy_0, tg_xxx_xyyyy_1, tg_xxx_xyyyz_0, \
                                         tg_xxx_xyyyz_1, tg_xxx_xyyzz_0, tg_xxx_xyyzz_1, tg_xxx_xyzzz_0, tg_xxx_xyzzz_1, \
                                         tg_xxx_xzzzz_0, tg_xxx_xzzzz_1, tg_xxx_yyyyy_0, tg_xxx_yyyyy_1, tg_xxx_yyyyz_0, \
                                         tg_xxx_yyyyz_1, tg_xxx_yyyzz_0, tg_xxx_yyyzz_1, tg_xxx_yyzzz_0, tg_xxx_yyzzz_1, \
                                         tg_xxx_yzzzz_0, tg_xxx_yzzzz_1, tg_xxx_zzzzz_0, tg_xxx_zzzzz_1, tg_xxxx_xxxx_1, \
                                         tg_xxxx_xxxxx_0, tg_xxxx_xxxxx_1, tg_xxxx_xxxxy_0, tg_xxxx_xxxxy_1, tg_xxxx_xxxxz_0, \
                                         tg_xxxx_xxxxz_1, tg_xxxx_xxxy_1, tg_xxxx_xxxyy_0, tg_xxxx_xxxyy_1, tg_xxxx_xxxyz_0, \
                                         tg_xxxx_xxxyz_1, tg_xxxx_xxxz_1, tg_xxxx_xxxzz_0, tg_xxxx_xxxzz_1, tg_xxxx_xxyy_1, \
                                         tg_xxxx_xxyyy_0, tg_xxxx_xxyyy_1, tg_xxxx_xxyyz_0, tg_xxxx_xxyyz_1, tg_xxxx_xxyz_1, \
                                         tg_xxxx_xxyzz_0, tg_xxxx_xxyzz_1, tg_xxxx_xxzz_1, tg_xxxx_xxzzz_0, tg_xxxx_xxzzz_1, \
                                         tg_xxxx_xyyy_1, tg_xxxx_xyyyy_0, tg_xxxx_xyyyy_1, tg_xxxx_xyyyz_0, tg_xxxx_xyyyz_1, \
                                         tg_xxxx_xyyz_1, tg_xxxx_xyyzz_0, tg_xxxx_xyyzz_1, tg_xxxx_xyzz_1, tg_xxxx_xyzzz_0, \
                                         tg_xxxx_xyzzz_1, tg_xxxx_xzzz_1, tg_xxxx_xzzzz_0, tg_xxxx_xzzzz_1, tg_xxxx_yyyy_1, \
                                         tg_xxxx_yyyyy_0, tg_xxxx_yyyyy_1, tg_xxxx_yyyyz_0, tg_xxxx_yyyyz_1, tg_xxxx_yyyz_1, \
                                         tg_xxxx_yyyzz_0, tg_xxxx_yyyzz_1, tg_xxxx_yyzz_1, tg_xxxx_yyzzz_0, tg_xxxx_yyzzz_1, \
                                         tg_xxxx_yzzz_1, tg_xxxx_yzzzz_0, tg_xxxx_yzzzz_1, tg_xxxx_zzzz_1, tg_xxxx_zzzzz_0, \
                                         tg_xxxx_zzzzz_1, tg_xxxxx_xxxxx_0, tg_xxxxx_xxxxy_0, tg_xxxxx_xxxxz_0, \
                                         tg_xxxxx_xxxyy_0, tg_xxxxx_xxxyz_0, tg_xxxxx_xxxzz_0, tg_xxxxx_xxyyy_0, \
                                         tg_xxxxx_xxyyz_0, tg_xxxxx_xxyzz_0, tg_xxxxx_xxzzz_0, tg_xxxxx_xyyyy_0, \
                                         tg_xxxxx_xyyyz_0, tg_xxxxx_xyyzz_0, tg_xxxxx_xyzzz_0, tg_xxxxx_xzzzz_0, \
                                         tg_xxxxx_yyyyy_0, tg_xxxxx_yyyyz_0, tg_xxxxx_yyyzz_0, tg_xxxxx_yyzzz_0, \
                                         tg_xxxxx_yzzzz_0, tg_xxxxx_zzzzz_0, tg_xxxxy_xxxxx_0, tg_xxxxy_xxxxy_0, \
                                         tg_xxxxy_xxxxz_0, tg_xxxxy_xxxyy_0, tg_xxxxy_xxxyz_0, tg_xxxxy_xxxzz_0, \
                                         tg_xxxxy_xxyyy_0, tg_xxxxy_xxyyz_0, tg_xxxxy_xxyzz_0, tg_xxxxy_xxzzz_0, \
                                         tg_xxxxy_xyyyy_0, tg_xxxxy_xyyyz_0, tg_xxxxy_xyyzz_0, tg_xxxxy_xyzzz_0, \
                                         tg_xxxxy_xzzzz_0, tg_xxxxy_yyyyy_0, tg_xxxxy_yyyyz_0, tg_xxxxy_yyyzz_0, \
                                         tg_xxxxy_yyzzz_0, tg_xxxxy_yzzzz_0, tg_xxxxy_zzzzz_0, tg_xxxxz_xxxxx_0, \
                                         tg_xxxxz_xxxxy_0, tg_xxxxz_xxxxz_0, tg_xxxxz_xxxyy_0, tg_xxxxz_xxxyz_0, \
                                         tg_xxxxz_xxxzz_0, tg_xxxxz_xxyyy_0, tg_xxxxz_xxyyz_0, tg_xxxxz_xxyzz_0, \
                                         tg_xxxxz_xxzzz_0, tg_xxxxz_xyyyy_0, tg_xxxxz_xyyyz_0, tg_xxxxz_xyyzz_0, \
                                         tg_xxxxz_xyzzz_0, tg_xxxxz_xzzzz_0, tg_xxxxz_yyyyy_0, tg_xxxxz_yyyyz_0, \
                                         tg_xxxxz_yyyzz_0, tg_xxxxz_yyzzz_0, tg_xxxxz_yzzzz_0, tg_xxxxz_zzzzz_0, \
                                         tg_xxxy_xxxx_1, tg_xxxy_xxxxx_0, tg_xxxy_xxxxx_1, tg_xxxy_xxxxy_0, tg_xxxy_xxxxy_1, \
                                         tg_xxxy_xxxxz_0, tg_xxxy_xxxxz_1, tg_xxxy_xxxy_1, tg_xxxy_xxxyy_0, tg_xxxy_xxxyy_1, \
                                         tg_xxxy_xxxyz_0, tg_xxxy_xxxyz_1, tg_xxxy_xxxz_1, tg_xxxy_xxxzz_0, tg_xxxy_xxxzz_1, \
                                         tg_xxxy_xxyy_1, tg_xxxy_xxyyy_0, tg_xxxy_xxyyy_1, tg_xxxy_xxyyz_0, tg_xxxy_xxyyz_1, \
                                         tg_xxxy_xxyz_1, tg_xxxy_xxyzz_0, tg_xxxy_xxyzz_1, tg_xxxy_xxzz_1, tg_xxxy_xxzzz_0, \
                                         tg_xxxy_xxzzz_1, tg_xxxy_xyyy_1, tg_xxxy_xyyyy_0, tg_xxxy_xyyyy_1, tg_xxxy_xyyyz_0, \
                                         tg_xxxy_xyyyz_1, tg_xxxy_xyyz_1, tg_xxxy_xyyzz_0, tg_xxxy_xyyzz_1, tg_xxxy_xyzz_1, \
                                         tg_xxxy_xyzzz_0, tg_xxxy_xyzzz_1, tg_xxxy_xzzz_1, tg_xxxy_xzzzz_0, tg_xxxy_xzzzz_1, \
                                         tg_xxxy_yyyy_1, tg_xxxy_yyyyy_0, tg_xxxy_yyyyy_1, tg_xxxy_yyyyz_0, tg_xxxy_yyyyz_1, \
                                         tg_xxxy_yyyz_1, tg_xxxy_yyyzz_0, tg_xxxy_yyyzz_1, tg_xxxy_yyzz_1, tg_xxxy_yyzzz_0, \
                                         tg_xxxy_yyzzz_1, tg_xxxy_yzzz_1, tg_xxxy_yzzzz_0, tg_xxxy_yzzzz_1, tg_xxxy_zzzz_1, \
                                         tg_xxxy_zzzzz_0, tg_xxxy_zzzzz_1, tg_xxxyy_xxxxx_0, tg_xxxyy_xxxxy_0, \
                                         tg_xxxyy_xxxxz_0, tg_xxxyy_xxxyy_0, tg_xxxyy_xxxyz_0, tg_xxxyy_xxxzz_0, \
                                         tg_xxxyy_xxyyy_0, tg_xxxyy_xxyyz_0, tg_xxxyy_xxyzz_0, tg_xxxyy_xxzzz_0, \
                                         tg_xxxyy_xyyyy_0, tg_xxxyy_xyyyz_0, tg_xxxyy_xyyzz_0, tg_xxxyy_xyzzz_0, \
                                         tg_xxxyy_xzzzz_0, tg_xxxyy_yyyyy_0, tg_xxxyy_yyyyz_0, tg_xxxyy_yyyzz_0, \
                                         tg_xxxyy_yyzzz_0, tg_xxxyy_yzzzz_0, tg_xxxyy_zzzzz_0, tg_xxxyz_xxxxx_0, \
                                         tg_xxxyz_xxxxy_0, tg_xxxyz_xxxxz_0, tg_xxxyz_xxxyy_0, tg_xxxyz_xxxyz_0, \
                                         tg_xxxz_xxxx_1, tg_xxxz_xxxxx_0, tg_xxxz_xxxxx_1, tg_xxxz_xxxxy_0, tg_xxxz_xxxxy_1, \
                                         tg_xxxz_xxxxz_0, tg_xxxz_xxxxz_1, tg_xxxz_xxxy_1, tg_xxxz_xxxyy_0, tg_xxxz_xxxyy_1, \
                                         tg_xxxz_xxxyz_0, tg_xxxz_xxxyz_1, tg_xxxz_xxxz_1, tg_xxxz_xxxzz_0, tg_xxxz_xxxzz_1, \
                                         tg_xxxz_xxyy_1, tg_xxxz_xxyyy_0, tg_xxxz_xxyyy_1, tg_xxxz_xxyyz_0, tg_xxxz_xxyyz_1, \
                                         tg_xxxz_xxyz_1, tg_xxxz_xxyzz_0, tg_xxxz_xxyzz_1, tg_xxxz_xxzz_1, tg_xxxz_xxzzz_0, \
                                         tg_xxxz_xxzzz_1, tg_xxxz_xyyy_1, tg_xxxz_xyyyy_0, tg_xxxz_xyyyy_1, tg_xxxz_xyyyz_0, \
                                         tg_xxxz_xyyyz_1, tg_xxxz_xyyz_1, tg_xxxz_xyyzz_0, tg_xxxz_xyyzz_1, tg_xxxz_xyzz_1, \
                                         tg_xxxz_xyzzz_0, tg_xxxz_xyzzz_1, tg_xxxz_xzzz_1, tg_xxxz_xzzzz_0, tg_xxxz_xzzzz_1, \
                                         tg_xxxz_yyyy_1, tg_xxxz_yyyyy_0, tg_xxxz_yyyyy_1, tg_xxxz_yyyyz_0, tg_xxxz_yyyyz_1, \
                                         tg_xxxz_yyyz_1, tg_xxxz_yyyzz_0, tg_xxxz_yyyzz_1, tg_xxxz_yyzz_1, tg_xxxz_yyzzz_0, \
                                         tg_xxxz_yyzzz_1, tg_xxxz_yzzz_1, tg_xxxz_yzzzz_0, tg_xxxz_yzzzz_1, tg_xxxz_zzzz_1, \
                                         tg_xxxz_zzzzz_0, tg_xxxz_zzzzz_1, tg_xxy_xxxxx_0, tg_xxy_xxxxx_1, tg_xxy_xxxxy_0, \
                                         tg_xxy_xxxxy_1, tg_xxy_xxxxz_0, tg_xxy_xxxxz_1, tg_xxy_xxxyy_0, tg_xxy_xxxyy_1, \
                                         tg_xxy_xxxyz_0, tg_xxy_xxxyz_1, tg_xxy_xxxzz_0, tg_xxy_xxxzz_1, tg_xxy_xxyyy_0, \
                                         tg_xxy_xxyyy_1, tg_xxy_xxyyz_0, tg_xxy_xxyyz_1, tg_xxy_xxyzz_0, tg_xxy_xxyzz_1, \
                                         tg_xxy_xxzzz_0, tg_xxy_xxzzz_1, tg_xxy_xyyyy_0, tg_xxy_xyyyy_1, tg_xxy_xyyyz_0, \
                                         tg_xxy_xyyyz_1, tg_xxy_xyyzz_0, tg_xxy_xyyzz_1, tg_xxy_xyzzz_0, tg_xxy_xyzzz_1, \
                                         tg_xxy_xzzzz_0, tg_xxy_xzzzz_1, tg_xxy_yyyyy_0, tg_xxy_yyyyy_1, tg_xxy_yyyyz_0, \
                                         tg_xxy_yyyyz_1, tg_xxy_yyyzz_0, tg_xxy_yyyzz_1, tg_xxy_yyzzz_0, tg_xxy_yyzzz_1, \
                                         tg_xxy_yzzzz_0, tg_xxy_yzzzz_1, tg_xxy_zzzzz_0, tg_xxy_zzzzz_1, tg_xxyy_xxxx_1, \
                                         tg_xxyy_xxxxx_0, tg_xxyy_xxxxx_1, tg_xxyy_xxxxy_0, tg_xxyy_xxxxy_1, tg_xxyy_xxxxz_0, \
                                         tg_xxyy_xxxxz_1, tg_xxyy_xxxy_1, tg_xxyy_xxxyy_0, tg_xxyy_xxxyy_1, tg_xxyy_xxxyz_0, \
                                         tg_xxyy_xxxyz_1, tg_xxyy_xxxz_1, tg_xxyy_xxxzz_0, tg_xxyy_xxxzz_1, tg_xxyy_xxyy_1, \
                                         tg_xxyy_xxyyy_0, tg_xxyy_xxyyy_1, tg_xxyy_xxyyz_0, tg_xxyy_xxyyz_1, tg_xxyy_xxyz_1, \
                                         tg_xxyy_xxyzz_0, tg_xxyy_xxyzz_1, tg_xxyy_xxzz_1, tg_xxyy_xxzzz_0, tg_xxyy_xxzzz_1, \
                                         tg_xxyy_xyyy_1, tg_xxyy_xyyyy_0, tg_xxyy_xyyyy_1, tg_xxyy_xyyyz_0, tg_xxyy_xyyyz_1, \
                                         tg_xxyy_xyyz_1, tg_xxyy_xyyzz_0, tg_xxyy_xyyzz_1, tg_xxyy_xyzz_1, tg_xxyy_xyzzz_0, \
                                         tg_xxyy_xyzzz_1, tg_xxyy_xzzz_1, tg_xxyy_xzzzz_0, tg_xxyy_xzzzz_1, tg_xxyy_yyyy_1, \
                                         tg_xxyy_yyyyy_0, tg_xxyy_yyyyy_1, tg_xxyy_yyyyz_0, tg_xxyy_yyyyz_1, tg_xxyy_yyyz_1, \
                                         tg_xxyy_yyyzz_0, tg_xxyy_yyyzz_1, tg_xxyy_yyzz_1, tg_xxyy_yyzzz_0, tg_xxyy_yyzzz_1, \
                                         tg_xxyy_yzzz_1, tg_xxyy_yzzzz_0, tg_xxyy_yzzzz_1, tg_xxyy_zzzz_1, tg_xxyy_zzzzz_0, \
                                         tg_xxyy_zzzzz_1, tg_xxyz_xxxx_1, tg_xxyz_xxxxx_0, tg_xxyz_xxxxx_1, tg_xxyz_xxxxy_0, \
                                         tg_xxyz_xxxxy_1, tg_xxyz_xxxxz_0, tg_xxyz_xxxxz_1, tg_xxyz_xxxy_1, tg_xxyz_xxxyy_0, \
                                         tg_xxyz_xxxyy_1, tg_xxyz_xxxyz_0, tg_xxyz_xxxyz_1, tg_xxyz_xxxz_1, tg_xxyz_xxyy_1, \
                                         tg_xxyz_xxyz_1, tg_xxz_xxxxx_0, tg_xxz_xxxxx_1, tg_xxz_xxxxy_0, tg_xxz_xxxxy_1, \
                                         tg_xxz_xxxxz_0, tg_xxz_xxxxz_1, tg_xxz_xxxyy_0, tg_xxz_xxxyy_1, tg_xxz_xxxyz_0, \
                                         tg_xxz_xxxyz_1, tg_xxz_xxxzz_0, tg_xxz_xxxzz_1, tg_xxz_xxyyy_0, tg_xxz_xxyyy_1, \
                                         tg_xxz_xxyyz_0, tg_xxz_xxyyz_1, tg_xxz_xxyzz_0, tg_xxz_xxyzz_1, tg_xxz_xxzzz_0, \
                                         tg_xxz_xxzzz_1, tg_xxz_xyyyy_0, tg_xxz_xyyyy_1, tg_xxz_xyyyz_0, tg_xxz_xyyyz_1, \
                                         tg_xxz_xyyzz_0, tg_xxz_xyyzz_1, tg_xxz_xyzzz_0, tg_xxz_xyzzz_1, tg_xxz_xzzzz_0, \
                                         tg_xxz_xzzzz_1, tg_xxz_yyyyy_0, tg_xxz_yyyyy_1, tg_xxz_yyyyz_0, tg_xxz_yyyyz_1, \
                                         tg_xxz_yyyzz_0, tg_xxz_yyyzz_1, tg_xxz_yyzzz_0, tg_xxz_yyzzz_1, tg_xxz_yzzzz_0, \
                                         tg_xxz_yzzzz_1, tg_xxz_zzzzz_0, tg_xxz_zzzzz_1, tg_xyy_xxxxx_0, tg_xyy_xxxxx_1, \
                                         tg_xyy_xxxxy_0, tg_xyy_xxxxy_1, tg_xyy_xxxxz_0, tg_xyy_xxxxz_1, tg_xyy_xxxyy_0, \
                                         tg_xyy_xxxyy_1, tg_xyy_xxxyz_0, tg_xyy_xxxyz_1, tg_xyy_xxxzz_0, tg_xyy_xxxzz_1, \
                                         tg_xyy_xxyyy_0, tg_xyy_xxyyy_1, tg_xyy_xxyyz_0, tg_xyy_xxyyz_1, tg_xyy_xxyzz_0, \
                                         tg_xyy_xxyzz_1, tg_xyy_xxzzz_0, tg_xyy_xxzzz_1, tg_xyy_xyyyy_0, tg_xyy_xyyyy_1, \
                                         tg_xyy_xyyyz_0, tg_xyy_xyyyz_1, tg_xyy_xyyzz_0, tg_xyy_xyyzz_1, tg_xyy_xyzzz_0, \
                                         tg_xyy_xyzzz_1, tg_xyy_xzzzz_0, tg_xyy_xzzzz_1, tg_xyy_yyyyy_0, tg_xyy_yyyyy_1, \
                                         tg_xyy_yyyyz_0, tg_xyy_yyyyz_1, tg_xyy_yyyzz_0, tg_xyy_yyyzz_1, tg_xyy_yyzzz_0, \
                                         tg_xyy_yyzzz_1, tg_xyy_yzzzz_0, tg_xyy_yzzzz_1, tg_xyy_zzzzz_0, tg_xyy_zzzzz_1, \
                                         tg_xyz_xxxxx_0, tg_xyz_xxxxx_1, tg_xyz_xxxxy_0, tg_xyz_xxxxy_1, tg_xyz_xxxxz_0, \
                                         tg_xyz_xxxxz_1, tg_xyz_xxxyy_0, tg_xyz_xxxyy_1, tg_xyz_xxxyz_0, tg_xyz_xxxyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxx_xxxxx_0[j] = pb_x * tg_xxxx_xxxxx_0[j] + fr * tg_xxxx_xxxxx_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxx_0[j] - tg_xxx_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxx_xxxx_1[j];

                    tg_xxxxx_xxxxy_0[j] = pb_x * tg_xxxx_xxxxy_0[j] + fr * tg_xxxx_xxxxy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxy_0[j] - tg_xxx_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxy_1[j];

                    tg_xxxxx_xxxxz_0[j] = pb_x * tg_xxxx_xxxxz_0[j] + fr * tg_xxxx_xxxxz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxz_0[j] - tg_xxx_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxz_1[j];

                    tg_xxxxx_xxxyy_0[j] = pb_x * tg_xxxx_xxxyy_0[j] + fr * tg_xxxx_xxxyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyy_0[j] - tg_xxx_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyy_1[j];

                    tg_xxxxx_xxxyz_0[j] = pb_x * tg_xxxx_xxxyz_0[j] + fr * tg_xxxx_xxxyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyz_0[j] - tg_xxx_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyz_1[j];

                    tg_xxxxx_xxxzz_0[j] = pb_x * tg_xxxx_xxxzz_0[j] + fr * tg_xxxx_xxxzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxzz_0[j] - tg_xxx_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxzz_1[j];

                    tg_xxxxx_xxyyy_0[j] = pb_x * tg_xxxx_xxyyy_0[j] + fr * tg_xxxx_xxyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyy_0[j] - tg_xxx_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyy_1[j];

                    tg_xxxxx_xxyyz_0[j] = pb_x * tg_xxxx_xxyyz_0[j] + fr * tg_xxxx_xxyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyz_0[j] - tg_xxx_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyz_1[j];

                    tg_xxxxx_xxyzz_0[j] = pb_x * tg_xxxx_xxyzz_0[j] + fr * tg_xxxx_xxyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyzz_0[j] - tg_xxx_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyzz_1[j];

                    tg_xxxxx_xxzzz_0[j] = pb_x * tg_xxxx_xxzzz_0[j] + fr * tg_xxxx_xxzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxzzz_0[j] - tg_xxx_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xzzz_1[j];

                    tg_xxxxx_xyyyy_0[j] = pb_x * tg_xxxx_xyyyy_0[j] + fr * tg_xxxx_xyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyy_0[j] - tg_xxx_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyy_1[j];

                    tg_xxxxx_xyyyz_0[j] = pb_x * tg_xxxx_xyyyz_0[j] + fr * tg_xxxx_xyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyz_0[j] - tg_xxx_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyz_1[j];

                    tg_xxxxx_xyyzz_0[j] = pb_x * tg_xxxx_xyyzz_0[j] + fr * tg_xxxx_xyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyzz_0[j] - tg_xxx_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyzz_1[j];

                    tg_xxxxx_xyzzz_0[j] = pb_x * tg_xxxx_xyzzz_0[j] + fr * tg_xxxx_xyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyzzz_0[j] - tg_xxx_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yzzz_1[j];

                    tg_xxxxx_xzzzz_0[j] = pb_x * tg_xxxx_xzzzz_0[j] + fr * tg_xxxx_xzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xzzzz_0[j] - tg_xxx_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_zzzz_1[j];

                    tg_xxxxx_yyyyy_0[j] = pb_x * tg_xxxx_yyyyy_0[j] + fr * tg_xxxx_yyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyy_0[j] - tg_xxx_yyyyy_1[j] * fl1_fza);

                    tg_xxxxx_yyyyz_0[j] = pb_x * tg_xxxx_yyyyz_0[j] + fr * tg_xxxx_yyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyz_0[j] - tg_xxx_yyyyz_1[j] * fl1_fza);

                    tg_xxxxx_yyyzz_0[j] = pb_x * tg_xxxx_yyyzz_0[j] + fr * tg_xxxx_yyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyzz_0[j] - tg_xxx_yyyzz_1[j] * fl1_fza);

                    tg_xxxxx_yyzzz_0[j] = pb_x * tg_xxxx_yyzzz_0[j] + fr * tg_xxxx_yyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyzzz_0[j] - tg_xxx_yyzzz_1[j] * fl1_fza);

                    tg_xxxxx_yzzzz_0[j] = pb_x * tg_xxxx_yzzzz_0[j] + fr * tg_xxxx_yzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yzzzz_0[j] - tg_xxx_yzzzz_1[j] * fl1_fza);

                    tg_xxxxx_zzzzz_0[j] = pb_x * tg_xxxx_zzzzz_0[j] + fr * tg_xxxx_zzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_zzzzz_0[j] - tg_xxx_zzzzz_1[j] * fl1_fza);

                    tg_xxxxy_xxxxx_0[j] = pb_x * tg_xxxy_xxxxx_0[j] + fr * tg_xxxy_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxx_0[j] - tg_xxy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxy_xxxx_1[j];

                    tg_xxxxy_xxxxy_0[j] = pb_x * tg_xxxy_xxxxy_0[j] + fr * tg_xxxy_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxy_0[j] - tg_xxy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxy_1[j];

                    tg_xxxxy_xxxxz_0[j] = pb_x * tg_xxxy_xxxxz_0[j] + fr * tg_xxxy_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxz_0[j] - tg_xxy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxz_1[j];

                    tg_xxxxy_xxxyy_0[j] = pb_x * tg_xxxy_xxxyy_0[j] + fr * tg_xxxy_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyy_0[j] - tg_xxy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyy_1[j];

                    tg_xxxxy_xxxyz_0[j] = pb_x * tg_xxxy_xxxyz_0[j] + fr * tg_xxxy_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyz_0[j] - tg_xxy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyz_1[j];

                    tg_xxxxy_xxxzz_0[j] = pb_x * tg_xxxy_xxxzz_0[j] + fr * tg_xxxy_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxzz_0[j] - tg_xxy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxzz_1[j];

                    tg_xxxxy_xxyyy_0[j] = pb_x * tg_xxxy_xxyyy_0[j] + fr * tg_xxxy_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyy_0[j] - tg_xxy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyy_1[j];

                    tg_xxxxy_xxyyz_0[j] = pb_x * tg_xxxy_xxyyz_0[j] + fr * tg_xxxy_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyz_0[j] - tg_xxy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyz_1[j];

                    tg_xxxxy_xxyzz_0[j] = pb_x * tg_xxxy_xxyzz_0[j] + fr * tg_xxxy_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyzz_0[j] - tg_xxy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyzz_1[j];

                    tg_xxxxy_xxzzz_0[j] = pb_x * tg_xxxy_xxzzz_0[j] + fr * tg_xxxy_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxzzz_0[j] - tg_xxy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xzzz_1[j];

                    tg_xxxxy_xyyyy_0[j] = pb_x * tg_xxxy_xyyyy_0[j] + fr * tg_xxxy_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyy_0[j] - tg_xxy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyy_1[j];

                    tg_xxxxy_xyyyz_0[j] = pb_x * tg_xxxy_xyyyz_0[j] + fr * tg_xxxy_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyz_0[j] - tg_xxy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyz_1[j];

                    tg_xxxxy_xyyzz_0[j] = pb_x * tg_xxxy_xyyzz_0[j] + fr * tg_xxxy_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyzz_0[j] - tg_xxy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyzz_1[j];

                    tg_xxxxy_xyzzz_0[j] = pb_x * tg_xxxy_xyzzz_0[j] + fr * tg_xxxy_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyzzz_0[j] - tg_xxy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yzzz_1[j];

                    tg_xxxxy_xzzzz_0[j] = pb_x * tg_xxxy_xzzzz_0[j] + fr * tg_xxxy_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xzzzz_0[j] - tg_xxy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_zzzz_1[j];

                    tg_xxxxy_yyyyy_0[j] = pb_x * tg_xxxy_yyyyy_0[j] + fr * tg_xxxy_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyy_0[j] - tg_xxy_yyyyy_1[j] * fl1_fza);

                    tg_xxxxy_yyyyz_0[j] = pb_x * tg_xxxy_yyyyz_0[j] + fr * tg_xxxy_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyz_0[j] - tg_xxy_yyyyz_1[j] * fl1_fza);

                    tg_xxxxy_yyyzz_0[j] = pb_x * tg_xxxy_yyyzz_0[j] + fr * tg_xxxy_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyzz_0[j] - tg_xxy_yyyzz_1[j] * fl1_fza);

                    tg_xxxxy_yyzzz_0[j] = pb_x * tg_xxxy_yyzzz_0[j] + fr * tg_xxxy_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyzzz_0[j] - tg_xxy_yyzzz_1[j] * fl1_fza);

                    tg_xxxxy_yzzzz_0[j] = pb_x * tg_xxxy_yzzzz_0[j] + fr * tg_xxxy_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yzzzz_0[j] - tg_xxy_yzzzz_1[j] * fl1_fza);

                    tg_xxxxy_zzzzz_0[j] = pb_x * tg_xxxy_zzzzz_0[j] + fr * tg_xxxy_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_zzzzz_0[j] - tg_xxy_zzzzz_1[j] * fl1_fza);

                    tg_xxxxz_xxxxx_0[j] = pb_x * tg_xxxz_xxxxx_0[j] + fr * tg_xxxz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxx_0[j] - tg_xxz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxz_xxxx_1[j];

                    tg_xxxxz_xxxxy_0[j] = pb_x * tg_xxxz_xxxxy_0[j] + fr * tg_xxxz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxy_0[j] - tg_xxz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxy_1[j];

                    tg_xxxxz_xxxxz_0[j] = pb_x * tg_xxxz_xxxxz_0[j] + fr * tg_xxxz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxz_0[j] - tg_xxz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxz_1[j];

                    tg_xxxxz_xxxyy_0[j] = pb_x * tg_xxxz_xxxyy_0[j] + fr * tg_xxxz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyy_0[j] - tg_xxz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyy_1[j];

                    tg_xxxxz_xxxyz_0[j] = pb_x * tg_xxxz_xxxyz_0[j] + fr * tg_xxxz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyz_0[j] - tg_xxz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyz_1[j];

                    tg_xxxxz_xxxzz_0[j] = pb_x * tg_xxxz_xxxzz_0[j] + fr * tg_xxxz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxzz_0[j] - tg_xxz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxzz_1[j];

                    tg_xxxxz_xxyyy_0[j] = pb_x * tg_xxxz_xxyyy_0[j] + fr * tg_xxxz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyy_0[j] - tg_xxz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyy_1[j];

                    tg_xxxxz_xxyyz_0[j] = pb_x * tg_xxxz_xxyyz_0[j] + fr * tg_xxxz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyz_0[j] - tg_xxz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyz_1[j];

                    tg_xxxxz_xxyzz_0[j] = pb_x * tg_xxxz_xxyzz_0[j] + fr * tg_xxxz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyzz_0[j] - tg_xxz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyzz_1[j];

                    tg_xxxxz_xxzzz_0[j] = pb_x * tg_xxxz_xxzzz_0[j] + fr * tg_xxxz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxzzz_0[j] - tg_xxz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xzzz_1[j];

                    tg_xxxxz_xyyyy_0[j] = pb_x * tg_xxxz_xyyyy_0[j] + fr * tg_xxxz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyy_0[j] - tg_xxz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyy_1[j];

                    tg_xxxxz_xyyyz_0[j] = pb_x * tg_xxxz_xyyyz_0[j] + fr * tg_xxxz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyz_0[j] - tg_xxz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyz_1[j];

                    tg_xxxxz_xyyzz_0[j] = pb_x * tg_xxxz_xyyzz_0[j] + fr * tg_xxxz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyzz_0[j] - tg_xxz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyzz_1[j];

                    tg_xxxxz_xyzzz_0[j] = pb_x * tg_xxxz_xyzzz_0[j] + fr * tg_xxxz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyzzz_0[j] - tg_xxz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yzzz_1[j];

                    tg_xxxxz_xzzzz_0[j] = pb_x * tg_xxxz_xzzzz_0[j] + fr * tg_xxxz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xzzzz_0[j] - tg_xxz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_zzzz_1[j];

                    tg_xxxxz_yyyyy_0[j] = pb_x * tg_xxxz_yyyyy_0[j] + fr * tg_xxxz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyy_0[j] - tg_xxz_yyyyy_1[j] * fl1_fza);

                    tg_xxxxz_yyyyz_0[j] = pb_x * tg_xxxz_yyyyz_0[j] + fr * tg_xxxz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyz_0[j] - tg_xxz_yyyyz_1[j] * fl1_fza);

                    tg_xxxxz_yyyzz_0[j] = pb_x * tg_xxxz_yyyzz_0[j] + fr * tg_xxxz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyzz_0[j] - tg_xxz_yyyzz_1[j] * fl1_fza);

                    tg_xxxxz_yyzzz_0[j] = pb_x * tg_xxxz_yyzzz_0[j] + fr * tg_xxxz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyzzz_0[j] - tg_xxz_yyzzz_1[j] * fl1_fza);

                    tg_xxxxz_yzzzz_0[j] = pb_x * tg_xxxz_yzzzz_0[j] + fr * tg_xxxz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yzzzz_0[j] - tg_xxz_yzzzz_1[j] * fl1_fza);

                    tg_xxxxz_zzzzz_0[j] = pb_x * tg_xxxz_zzzzz_0[j] + fr * tg_xxxz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_zzzzz_0[j] - tg_xxz_zzzzz_1[j] * fl1_fza);

                    tg_xxxyy_xxxxx_0[j] = pb_x * tg_xxyy_xxxxx_0[j] + fr * tg_xxyy_xxxxx_1[j] + fl1_fx * (tg_xyy_xxxxx_0[j] - tg_xyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyy_xxxx_1[j];

                    tg_xxxyy_xxxxy_0[j] = pb_x * tg_xxyy_xxxxy_0[j] + fr * tg_xxyy_xxxxy_1[j] + fl1_fx * (tg_xyy_xxxxy_0[j] - tg_xyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxy_1[j];

                    tg_xxxyy_xxxxz_0[j] = pb_x * tg_xxyy_xxxxz_0[j] + fr * tg_xxyy_xxxxz_1[j] + fl1_fx * (tg_xyy_xxxxz_0[j] - tg_xyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxz_1[j];

                    tg_xxxyy_xxxyy_0[j] = pb_x * tg_xxyy_xxxyy_0[j] + fr * tg_xxyy_xxxyy_1[j] + fl1_fx * (tg_xyy_xxxyy_0[j] - tg_xyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyy_1[j];

                    tg_xxxyy_xxxyz_0[j] = pb_x * tg_xxyy_xxxyz_0[j] + fr * tg_xxyy_xxxyz_1[j] + fl1_fx * (tg_xyy_xxxyz_0[j] - tg_xyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyz_1[j];

                    tg_xxxyy_xxxzz_0[j] = pb_x * tg_xxyy_xxxzz_0[j] + fr * tg_xxyy_xxxzz_1[j] + fl1_fx * (tg_xyy_xxxzz_0[j] - tg_xyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxzz_1[j];

                    tg_xxxyy_xxyyy_0[j] = pb_x * tg_xxyy_xxyyy_0[j] + fr * tg_xxyy_xxyyy_1[j] + fl1_fx * (tg_xyy_xxyyy_0[j] - tg_xyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyy_1[j];

                    tg_xxxyy_xxyyz_0[j] = pb_x * tg_xxyy_xxyyz_0[j] + fr * tg_xxyy_xxyyz_1[j] + fl1_fx * (tg_xyy_xxyyz_0[j] - tg_xyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyz_1[j];

                    tg_xxxyy_xxyzz_0[j] = pb_x * tg_xxyy_xxyzz_0[j] + fr * tg_xxyy_xxyzz_1[j] + fl1_fx * (tg_xyy_xxyzz_0[j] - tg_xyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyzz_1[j];

                    tg_xxxyy_xxzzz_0[j] = pb_x * tg_xxyy_xxzzz_0[j] + fr * tg_xxyy_xxzzz_1[j] + fl1_fx * (tg_xyy_xxzzz_0[j] - tg_xyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xzzz_1[j];

                    tg_xxxyy_xyyyy_0[j] = pb_x * tg_xxyy_xyyyy_0[j] + fr * tg_xxyy_xyyyy_1[j] + fl1_fx * (tg_xyy_xyyyy_0[j] - tg_xyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyy_1[j];

                    tg_xxxyy_xyyyz_0[j] = pb_x * tg_xxyy_xyyyz_0[j] + fr * tg_xxyy_xyyyz_1[j] + fl1_fx * (tg_xyy_xyyyz_0[j] - tg_xyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyz_1[j];

                    tg_xxxyy_xyyzz_0[j] = pb_x * tg_xxyy_xyyzz_0[j] + fr * tg_xxyy_xyyzz_1[j] + fl1_fx * (tg_xyy_xyyzz_0[j] - tg_xyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyzz_1[j];

                    tg_xxxyy_xyzzz_0[j] = pb_x * tg_xxyy_xyzzz_0[j] + fr * tg_xxyy_xyzzz_1[j] + fl1_fx * (tg_xyy_xyzzz_0[j] - tg_xyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yzzz_1[j];

                    tg_xxxyy_xzzzz_0[j] = pb_x * tg_xxyy_xzzzz_0[j] + fr * tg_xxyy_xzzzz_1[j] + fl1_fx * (tg_xyy_xzzzz_0[j] - tg_xyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_zzzz_1[j];

                    tg_xxxyy_yyyyy_0[j] = pb_x * tg_xxyy_yyyyy_0[j] + fr * tg_xxyy_yyyyy_1[j] + fl1_fx * (tg_xyy_yyyyy_0[j] - tg_xyy_yyyyy_1[j] * fl1_fza);

                    tg_xxxyy_yyyyz_0[j] = pb_x * tg_xxyy_yyyyz_0[j] + fr * tg_xxyy_yyyyz_1[j] + fl1_fx * (tg_xyy_yyyyz_0[j] - tg_xyy_yyyyz_1[j] * fl1_fza);

                    tg_xxxyy_yyyzz_0[j] = pb_x * tg_xxyy_yyyzz_0[j] + fr * tg_xxyy_yyyzz_1[j] + fl1_fx * (tg_xyy_yyyzz_0[j] - tg_xyy_yyyzz_1[j] * fl1_fza);

                    tg_xxxyy_yyzzz_0[j] = pb_x * tg_xxyy_yyzzz_0[j] + fr * tg_xxyy_yyzzz_1[j] + fl1_fx * (tg_xyy_yyzzz_0[j] - tg_xyy_yyzzz_1[j] * fl1_fza);

                    tg_xxxyy_yzzzz_0[j] = pb_x * tg_xxyy_yzzzz_0[j] + fr * tg_xxyy_yzzzz_1[j] + fl1_fx * (tg_xyy_yzzzz_0[j] - tg_xyy_yzzzz_1[j] * fl1_fza);

                    tg_xxxyy_zzzzz_0[j] = pb_x * tg_xxyy_zzzzz_0[j] + fr * tg_xxyy_zzzzz_1[j] + fl1_fx * (tg_xyy_zzzzz_0[j] - tg_xyy_zzzzz_1[j] * fl1_fza);

                    tg_xxxyz_xxxxx_0[j] = pb_x * tg_xxyz_xxxxx_0[j] + fr * tg_xxyz_xxxxx_1[j] + fl1_fx * (tg_xyz_xxxxx_0[j] - tg_xyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyz_xxxx_1[j];

                    tg_xxxyz_xxxxy_0[j] = pb_x * tg_xxyz_xxxxy_0[j] + fr * tg_xxyz_xxxxy_1[j] + fl1_fx * (tg_xyz_xxxxy_0[j] - tg_xyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxy_1[j];

                    tg_xxxyz_xxxxz_0[j] = pb_x * tg_xxyz_xxxxz_0[j] + fr * tg_xxyz_xxxxz_1[j] + fl1_fx * (tg_xyz_xxxxz_0[j] - tg_xyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxz_1[j];

                    tg_xxxyz_xxxyy_0[j] = pb_x * tg_xxyz_xxxyy_0[j] + fr * tg_xxyz_xxxyy_1[j] + fl1_fx * (tg_xyz_xxxyy_0[j] - tg_xyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyy_1[j];

                    tg_xxxyz_xxxyz_0[j] = pb_x * tg_xxyz_xxxyz_0[j] + fr * tg_xxyz_xxxyz_1[j] + fl1_fx * (tg_xyz_xxxyz_0[j] - tg_xyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSH_89_177(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (89,177)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xxyz_xxxzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 89); 

                auto tg_xxyz_xxyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 90); 

                auto tg_xxyz_xxyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 91); 

                auto tg_xxyz_xxyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 92); 

                auto tg_xxyz_xxzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 93); 

                auto tg_xxyz_xyyyy_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 94); 

                auto tg_xxyz_xyyyz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 95); 

                auto tg_xxyz_xyyzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 96); 

                auto tg_xxyz_xyzzz_1 = primBuffer[pidx_g_4_5_m1].data(315 * idx + 97); 

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

                auto tg_xxyz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 65); 

                auto tg_xxyz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 66); 

                auto tg_xxyz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 67); 

                auto tg_xxyz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 68); 

                auto tg_xxyz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 69); 

                auto tg_xxyz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 70); 

                auto tg_xxyz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 71); 

                auto tg_xxyz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 72); 

                auto tg_xxyz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 73); 

                auto tg_xxyz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 74); 

                auto tg_xxzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 75); 

                auto tg_xxzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 76); 

                auto tg_xxzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 77); 

                auto tg_xxzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 78); 

                auto tg_xxzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 79); 

                auto tg_xxzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 80); 

                auto tg_xxzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 81); 

                auto tg_xxzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 82); 

                auto tg_xxzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 83); 

                auto tg_xxzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 84); 

                auto tg_xxzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 85); 

                auto tg_xxzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 86); 

                auto tg_xxzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 87); 

                auto tg_xxzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 88); 

                auto tg_xxzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 89); 

                auto tg_xyyy_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 90); 

                auto tg_xyyy_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 91); 

                auto tg_xyyy_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 92); 

                auto tg_xyyy_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 93); 

                auto tg_xyyy_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 94); 

                auto tg_xyyy_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 95); 

                auto tg_xyyy_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 96); 

                auto tg_xyyy_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 97); 

                auto tg_xyyy_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 98); 

                auto tg_xyyy_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 99); 

                auto tg_xyyy_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 100); 

                auto tg_xyyy_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 101); 

                auto tg_xyyy_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 102); 

                auto tg_xyyy_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 103); 

                auto tg_xyyy_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 104); 

                auto tg_xyyz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 105); 

                auto tg_xyyz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 106); 

                auto tg_xyyz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 107); 

                auto tg_xyyz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 108); 

                auto tg_xyyz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 109); 

                auto tg_xyyz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 110); 

                auto tg_xyyz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 111); 

                auto tg_xyyz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 112); 

                auto tg_xyyz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 113); 

                auto tg_xyyz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 114); 

                auto tg_xyyz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 115); 

                auto tg_xyyz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 116); 

                auto tg_xyyz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 117); 

                auto tg_xyyz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 118); 

                auto tg_xyyz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 119); 

                auto tg_xyzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 120); 

                auto tg_xyzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 121); 

                auto tg_xyzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 122); 

                auto tg_xyzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 123); 

                auto tg_xyzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 124); 

                auto tg_xyzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 125); 

                auto tg_xyzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 126); 

                auto tg_xyzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 127); 

                auto tg_xyzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 128); 

                // set up pointers to integrals

                auto tg_xxxyz_xxxzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 89); 

                auto tg_xxxyz_xxyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 90); 

                auto tg_xxxyz_xxyyz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 91); 

                auto tg_xxxyz_xxyzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 92); 

                auto tg_xxxyz_xxzzz_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 93); 

                auto tg_xxxyz_xyyyy_0 = primBuffer[pidx_g_5_5_m0].data(441 * idx + 94); 

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

                // Batch of Integrals (89,177)

                #pragma omp simd aligned(fxn, fza, tg_xxxyz_xxxzz_0, tg_xxxyz_xxyyy_0, tg_xxxyz_xxyyz_0, \
                                         tg_xxxyz_xxyzz_0, tg_xxxyz_xxzzz_0, tg_xxxyz_xyyyy_0, tg_xxxyz_xyyyz_0, \
                                         tg_xxxyz_xyyzz_0, tg_xxxyz_xyzzz_0, tg_xxxyz_xzzzz_0, tg_xxxyz_yyyyy_0, \
                                         tg_xxxyz_yyyyz_0, tg_xxxyz_yyyzz_0, tg_xxxyz_yyzzz_0, tg_xxxyz_yzzzz_0, \
                                         tg_xxxyz_zzzzz_0, tg_xxxzz_xxxxx_0, tg_xxxzz_xxxxy_0, tg_xxxzz_xxxxz_0, \
                                         tg_xxxzz_xxxyy_0, tg_xxxzz_xxxyz_0, tg_xxxzz_xxxzz_0, tg_xxxzz_xxyyy_0, \
                                         tg_xxxzz_xxyyz_0, tg_xxxzz_xxyzz_0, tg_xxxzz_xxzzz_0, tg_xxxzz_xyyyy_0, \
                                         tg_xxxzz_xyyyz_0, tg_xxxzz_xyyzz_0, tg_xxxzz_xyzzz_0, tg_xxxzz_xzzzz_0, \
                                         tg_xxxzz_yyyyy_0, tg_xxxzz_yyyyz_0, tg_xxxzz_yyyzz_0, tg_xxxzz_yyzzz_0, \
                                         tg_xxxzz_yzzzz_0, tg_xxxzz_zzzzz_0, tg_xxyyy_xxxxx_0, tg_xxyyy_xxxxy_0, \
                                         tg_xxyyy_xxxxz_0, tg_xxyyy_xxxyy_0, tg_xxyyy_xxxyz_0, tg_xxyyy_xxxzz_0, \
                                         tg_xxyyy_xxyyy_0, tg_xxyyy_xxyyz_0, tg_xxyyy_xxyzz_0, tg_xxyyy_xxzzz_0, \
                                         tg_xxyyy_xyyyy_0, tg_xxyyy_xyyyz_0, tg_xxyyy_xyyzz_0, tg_xxyyy_xyzzz_0, \
                                         tg_xxyyy_xzzzz_0, tg_xxyyy_yyyyy_0, tg_xxyyy_yyyyz_0, tg_xxyyy_yyyzz_0, \
                                         tg_xxyyy_yyzzz_0, tg_xxyyy_yzzzz_0, tg_xxyyy_zzzzz_0, tg_xxyyz_xxxxx_0, \
                                         tg_xxyyz_xxxxy_0, tg_xxyyz_xxxxz_0, tg_xxyyz_xxxyy_0, tg_xxyyz_xxxyz_0, \
                                         tg_xxyyz_xxxzz_0, tg_xxyyz_xxyyy_0, tg_xxyyz_xxyyz_0, tg_xxyyz_xxyzz_0, \
                                         tg_xxyyz_xxzzz_0, tg_xxyyz_xyyyy_0, tg_xxyyz_xyyyz_0, tg_xxyyz_xyyzz_0, \
                                         tg_xxyyz_xyzzz_0, tg_xxyyz_xzzzz_0, tg_xxyyz_yyyyy_0, tg_xxyyz_yyyyz_0, \
                                         tg_xxyyz_yyyzz_0, tg_xxyyz_yyzzz_0, tg_xxyyz_yzzzz_0, tg_xxyyz_zzzzz_0, \
                                         tg_xxyz_xxxzz_0, tg_xxyz_xxxzz_1, tg_xxyz_xxyyy_0, tg_xxyz_xxyyy_1, tg_xxyz_xxyyz_0, \
                                         tg_xxyz_xxyyz_1, tg_xxyz_xxyzz_0, tg_xxyz_xxyzz_1, tg_xxyz_xxzz_1, tg_xxyz_xxzzz_0, \
                                         tg_xxyz_xxzzz_1, tg_xxyz_xyyy_1, tg_xxyz_xyyyy_0, tg_xxyz_xyyyy_1, tg_xxyz_xyyyz_0, \
                                         tg_xxyz_xyyyz_1, tg_xxyz_xyyz_1, tg_xxyz_xyyzz_0, tg_xxyz_xyyzz_1, tg_xxyz_xyzz_1, \
                                         tg_xxyz_xyzzz_0, tg_xxyz_xyzzz_1, tg_xxyz_xzzz_1, tg_xxyz_xzzzz_0, tg_xxyz_xzzzz_1, \
                                         tg_xxyz_yyyy_1, tg_xxyz_yyyyy_0, tg_xxyz_yyyyy_1, tg_xxyz_yyyyz_0, tg_xxyz_yyyyz_1, \
                                         tg_xxyz_yyyz_1, tg_xxyz_yyyzz_0, tg_xxyz_yyyzz_1, tg_xxyz_yyzz_1, tg_xxyz_yyzzz_0, \
                                         tg_xxyz_yyzzz_1, tg_xxyz_yzzz_1, tg_xxyz_yzzzz_0, tg_xxyz_yzzzz_1, tg_xxyz_zzzz_1, \
                                         tg_xxyz_zzzzz_0, tg_xxyz_zzzzz_1, tg_xxyzz_xxxxx_0, tg_xxyzz_xxxxy_0, \
                                         tg_xxyzz_xxxxz_0, tg_xxyzz_xxxyy_0, tg_xxyzz_xxxyz_0, tg_xxyzz_xxxzz_0, \
                                         tg_xxyzz_xxyyy_0, tg_xxyzz_xxyyz_0, tg_xxyzz_xxyzz_0, tg_xxzz_xxxx_1, \
                                         tg_xxzz_xxxxx_0, tg_xxzz_xxxxx_1, tg_xxzz_xxxxy_0, tg_xxzz_xxxxy_1, tg_xxzz_xxxxz_0, \
                                         tg_xxzz_xxxxz_1, tg_xxzz_xxxy_1, tg_xxzz_xxxyy_0, tg_xxzz_xxxyy_1, tg_xxzz_xxxyz_0, \
                                         tg_xxzz_xxxyz_1, tg_xxzz_xxxz_1, tg_xxzz_xxxzz_0, tg_xxzz_xxxzz_1, tg_xxzz_xxyy_1, \
                                         tg_xxzz_xxyyy_0, tg_xxzz_xxyyy_1, tg_xxzz_xxyyz_0, tg_xxzz_xxyyz_1, tg_xxzz_xxyz_1, \
                                         tg_xxzz_xxyzz_0, tg_xxzz_xxyzz_1, tg_xxzz_xxzz_1, tg_xxzz_xxzzz_0, tg_xxzz_xxzzz_1, \
                                         tg_xxzz_xyyy_1, tg_xxzz_xyyyy_0, tg_xxzz_xyyyy_1, tg_xxzz_xyyyz_0, tg_xxzz_xyyyz_1, \
                                         tg_xxzz_xyyz_1, tg_xxzz_xyyzz_0, tg_xxzz_xyyzz_1, tg_xxzz_xyzz_1, tg_xxzz_xyzzz_0, \
                                         tg_xxzz_xyzzz_1, tg_xxzz_xzzz_1, tg_xxzz_xzzzz_0, tg_xxzz_xzzzz_1, tg_xxzz_yyyy_1, \
                                         tg_xxzz_yyyyy_0, tg_xxzz_yyyyy_1, tg_xxzz_yyyyz_0, tg_xxzz_yyyyz_1, tg_xxzz_yyyz_1, \
                                         tg_xxzz_yyyzz_0, tg_xxzz_yyyzz_1, tg_xxzz_yyzz_1, tg_xxzz_yyzzz_0, tg_xxzz_yyzzz_1, \
                                         tg_xxzz_yzzz_1, tg_xxzz_yzzzz_0, tg_xxzz_yzzzz_1, tg_xxzz_zzzz_1, tg_xxzz_zzzzz_0, \
                                         tg_xxzz_zzzzz_1, tg_xyyy_xxxx_1, tg_xyyy_xxxxx_0, tg_xyyy_xxxxx_1, tg_xyyy_xxxxy_0, \
                                         tg_xyyy_xxxxy_1, tg_xyyy_xxxxz_0, tg_xyyy_xxxxz_1, tg_xyyy_xxxy_1, tg_xyyy_xxxyy_0, \
                                         tg_xyyy_xxxyy_1, tg_xyyy_xxxyz_0, tg_xyyy_xxxyz_1, tg_xyyy_xxxz_1, tg_xyyy_xxxzz_0, \
                                         tg_xyyy_xxxzz_1, tg_xyyy_xxyy_1, tg_xyyy_xxyyy_0, tg_xyyy_xxyyy_1, tg_xyyy_xxyyz_0, \
                                         tg_xyyy_xxyyz_1, tg_xyyy_xxyz_1, tg_xyyy_xxyzz_0, tg_xyyy_xxyzz_1, tg_xyyy_xxzz_1, \
                                         tg_xyyy_xxzzz_0, tg_xyyy_xxzzz_1, tg_xyyy_xyyy_1, tg_xyyy_xyyyy_0, tg_xyyy_xyyyy_1, \
                                         tg_xyyy_xyyyz_0, tg_xyyy_xyyyz_1, tg_xyyy_xyyz_1, tg_xyyy_xyyzz_0, tg_xyyy_xyyzz_1, \
                                         tg_xyyy_xyzz_1, tg_xyyy_xyzzz_0, tg_xyyy_xyzzz_1, tg_xyyy_xzzz_1, tg_xyyy_xzzzz_0, \
                                         tg_xyyy_xzzzz_1, tg_xyyy_yyyy_1, tg_xyyy_yyyyy_0, tg_xyyy_yyyyy_1, tg_xyyy_yyyyz_0, \
                                         tg_xyyy_yyyyz_1, tg_xyyy_yyyz_1, tg_xyyy_yyyzz_0, tg_xyyy_yyyzz_1, tg_xyyy_yyzz_1, \
                                         tg_xyyy_yyzzz_0, tg_xyyy_yyzzz_1, tg_xyyy_yzzz_1, tg_xyyy_yzzzz_0, tg_xyyy_yzzzz_1, \
                                         tg_xyyy_zzzz_1, tg_xyyy_zzzzz_0, tg_xyyy_zzzzz_1, tg_xyyz_xxxx_1, tg_xyyz_xxxxx_0, \
                                         tg_xyyz_xxxxx_1, tg_xyyz_xxxxy_0, tg_xyyz_xxxxy_1, tg_xyyz_xxxxz_0, tg_xyyz_xxxxz_1, \
                                         tg_xyyz_xxxy_1, tg_xyyz_xxxyy_0, tg_xyyz_xxxyy_1, tg_xyyz_xxxyz_0, tg_xyyz_xxxyz_1, \
                                         tg_xyyz_xxxz_1, tg_xyyz_xxxzz_0, tg_xyyz_xxxzz_1, tg_xyyz_xxyy_1, tg_xyyz_xxyyy_0, \
                                         tg_xyyz_xxyyy_1, tg_xyyz_xxyyz_0, tg_xyyz_xxyyz_1, tg_xyyz_xxyz_1, tg_xyyz_xxyzz_0, \
                                         tg_xyyz_xxyzz_1, tg_xyyz_xxzz_1, tg_xyyz_xxzzz_0, tg_xyyz_xxzzz_1, tg_xyyz_xyyy_1, \
                                         tg_xyyz_xyyyy_0, tg_xyyz_xyyyy_1, tg_xyyz_xyyyz_0, tg_xyyz_xyyyz_1, tg_xyyz_xyyz_1, \
                                         tg_xyyz_xyyzz_0, tg_xyyz_xyyzz_1, tg_xyyz_xyzz_1, tg_xyyz_xyzzz_0, tg_xyyz_xyzzz_1, \
                                         tg_xyyz_xzzz_1, tg_xyyz_xzzzz_0, tg_xyyz_xzzzz_1, tg_xyyz_yyyy_1, tg_xyyz_yyyyy_0, \
                                         tg_xyyz_yyyyy_1, tg_xyyz_yyyyz_0, tg_xyyz_yyyyz_1, tg_xyyz_yyyz_1, tg_xyyz_yyyzz_0, \
                                         tg_xyyz_yyyzz_1, tg_xyyz_yyzz_1, tg_xyyz_yyzzz_0, tg_xyyz_yyzzz_1, tg_xyyz_yzzz_1, \
                                         tg_xyyz_yzzzz_0, tg_xyyz_yzzzz_1, tg_xyyz_zzzz_1, tg_xyyz_zzzzz_0, tg_xyyz_zzzzz_1, \
                                         tg_xyz_xxxzz_0, tg_xyz_xxxzz_1, tg_xyz_xxyyy_0, tg_xyz_xxyyy_1, tg_xyz_xxyyz_0, \
                                         tg_xyz_xxyyz_1, tg_xyz_xxyzz_0, tg_xyz_xxyzz_1, tg_xyz_xxzzz_0, tg_xyz_xxzzz_1, \
                                         tg_xyz_xyyyy_0, tg_xyz_xyyyy_1, tg_xyz_xyyyz_0, tg_xyz_xyyyz_1, tg_xyz_xyyzz_0, \
                                         tg_xyz_xyyzz_1, tg_xyz_xyzzz_0, tg_xyz_xyzzz_1, tg_xyz_xzzzz_0, tg_xyz_xzzzz_1, \
                                         tg_xyz_yyyyy_0, tg_xyz_yyyyy_1, tg_xyz_yyyyz_0, tg_xyz_yyyyz_1, tg_xyz_yyyzz_0, \
                                         tg_xyz_yyyzz_1, tg_xyz_yyzzz_0, tg_xyz_yyzzz_1, tg_xyz_yzzzz_0, tg_xyz_yzzzz_1, \
                                         tg_xyz_zzzzz_0, tg_xyz_zzzzz_1, tg_xyzz_xxxx_1, tg_xyzz_xxxxx_0, tg_xyzz_xxxxx_1, \
                                         tg_xyzz_xxxxy_0, tg_xyzz_xxxxy_1, tg_xyzz_xxxxz_0, tg_xyzz_xxxxz_1, tg_xyzz_xxxy_1, \
                                         tg_xyzz_xxxyy_0, tg_xyzz_xxxyy_1, tg_xyzz_xxxyz_0, tg_xyzz_xxxyz_1, tg_xyzz_xxxz_1, \
                                         tg_xyzz_xxxzz_0, tg_xyzz_xxxzz_1, tg_xyzz_xxyy_1, tg_xyzz_xxyyy_0, tg_xyzz_xxyyy_1, \
                                         tg_xyzz_xxyyz_0, tg_xyzz_xxyyz_1, tg_xyzz_xxyz_1, tg_xyzz_xxyzz_0, tg_xyzz_xxyzz_1, \
                                         tg_xyzz_xxzz_1, tg_xyzz_xyyy_1, tg_xyzz_xyyz_1, tg_xyzz_xyzz_1, tg_xzz_xxxxx_0, \
                                         tg_xzz_xxxxx_1, tg_xzz_xxxxy_0, tg_xzz_xxxxy_1, tg_xzz_xxxxz_0, tg_xzz_xxxxz_1, \
                                         tg_xzz_xxxyy_0, tg_xzz_xxxyy_1, tg_xzz_xxxyz_0, tg_xzz_xxxyz_1, tg_xzz_xxxzz_0, \
                                         tg_xzz_xxxzz_1, tg_xzz_xxyyy_0, tg_xzz_xxyyy_1, tg_xzz_xxyyz_0, tg_xzz_xxyyz_1, \
                                         tg_xzz_xxyzz_0, tg_xzz_xxyzz_1, tg_xzz_xxzzz_0, tg_xzz_xxzzz_1, tg_xzz_xyyyy_0, \
                                         tg_xzz_xyyyy_1, tg_xzz_xyyyz_0, tg_xzz_xyyyz_1, tg_xzz_xyyzz_0, tg_xzz_xyyzz_1, \
                                         tg_xzz_xyzzz_0, tg_xzz_xyzzz_1, tg_xzz_xzzzz_0, tg_xzz_xzzzz_1, tg_xzz_yyyyy_0, \
                                         tg_xzz_yyyyy_1, tg_xzz_yyyyz_0, tg_xzz_yyyyz_1, tg_xzz_yyyzz_0, tg_xzz_yyyzz_1, \
                                         tg_xzz_yyzzz_0, tg_xzz_yyzzz_1, tg_xzz_yzzzz_0, tg_xzz_yzzzz_1, tg_xzz_zzzzz_0, \
                                         tg_xzz_zzzzz_1, tg_yyy_xxxxx_0, tg_yyy_xxxxx_1, tg_yyy_xxxxy_0, tg_yyy_xxxxy_1, \
                                         tg_yyy_xxxxz_0, tg_yyy_xxxxz_1, tg_yyy_xxxyy_0, tg_yyy_xxxyy_1, tg_yyy_xxxyz_0, \
                                         tg_yyy_xxxyz_1, tg_yyy_xxxzz_0, tg_yyy_xxxzz_1, tg_yyy_xxyyy_0, tg_yyy_xxyyy_1, \
                                         tg_yyy_xxyyz_0, tg_yyy_xxyyz_1, tg_yyy_xxyzz_0, tg_yyy_xxyzz_1, tg_yyy_xxzzz_0, \
                                         tg_yyy_xxzzz_1, tg_yyy_xyyyy_0, tg_yyy_xyyyy_1, tg_yyy_xyyyz_0, tg_yyy_xyyyz_1, \
                                         tg_yyy_xyyzz_0, tg_yyy_xyyzz_1, tg_yyy_xyzzz_0, tg_yyy_xyzzz_1, tg_yyy_xzzzz_0, \
                                         tg_yyy_xzzzz_1, tg_yyy_yyyyy_0, tg_yyy_yyyyy_1, tg_yyy_yyyyz_0, tg_yyy_yyyyz_1, \
                                         tg_yyy_yyyzz_0, tg_yyy_yyyzz_1, tg_yyy_yyzzz_0, tg_yyy_yyzzz_1, tg_yyy_yzzzz_0, \
                                         tg_yyy_yzzzz_1, tg_yyy_zzzzz_0, tg_yyy_zzzzz_1, tg_yyz_xxxxx_0, tg_yyz_xxxxx_1, \
                                         tg_yyz_xxxxy_0, tg_yyz_xxxxy_1, tg_yyz_xxxxz_0, tg_yyz_xxxxz_1, tg_yyz_xxxyy_0, \
                                         tg_yyz_xxxyy_1, tg_yyz_xxxyz_0, tg_yyz_xxxyz_1, tg_yyz_xxxzz_0, tg_yyz_xxxzz_1, \
                                         tg_yyz_xxyyy_0, tg_yyz_xxyyy_1, tg_yyz_xxyyz_0, tg_yyz_xxyyz_1, tg_yyz_xxyzz_0, \
                                         tg_yyz_xxyzz_1, tg_yyz_xxzzz_0, tg_yyz_xxzzz_1, tg_yyz_xyyyy_0, tg_yyz_xyyyy_1, \
                                         tg_yyz_xyyyz_0, tg_yyz_xyyyz_1, tg_yyz_xyyzz_0, tg_yyz_xyyzz_1, tg_yyz_xyzzz_0, \
                                         tg_yyz_xyzzz_1, tg_yyz_xzzzz_0, tg_yyz_xzzzz_1, tg_yyz_yyyyy_0, tg_yyz_yyyyy_1, \
                                         tg_yyz_yyyyz_0, tg_yyz_yyyyz_1, tg_yyz_yyyzz_0, tg_yyz_yyyzz_1, tg_yyz_yyzzz_0, \
                                         tg_yyz_yyzzz_1, tg_yyz_yzzzz_0, tg_yyz_yzzzz_1, tg_yyz_zzzzz_0, tg_yyz_zzzzz_1, \
                                         tg_yzz_xxxxx_0, tg_yzz_xxxxx_1, tg_yzz_xxxxy_0, tg_yzz_xxxxy_1, tg_yzz_xxxxz_0, \
                                         tg_yzz_xxxxz_1, tg_yzz_xxxyy_0, tg_yzz_xxxyy_1, tg_yzz_xxxyz_0, tg_yzz_xxxyz_1, \
                                         tg_yzz_xxxzz_0, tg_yzz_xxxzz_1, tg_yzz_xxyyy_0, tg_yzz_xxyyy_1, tg_yzz_xxyyz_0, \
                                         tg_yzz_xxyyz_1, tg_yzz_xxyzz_0, tg_yzz_xxyzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyz_xxxzz_0[j] = pb_x * tg_xxyz_xxxzz_0[j] + fr * tg_xxyz_xxxzz_1[j] + fl1_fx * (tg_xyz_xxxzz_0[j] - tg_xyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxzz_1[j];

                    tg_xxxyz_xxyyy_0[j] = pb_x * tg_xxyz_xxyyy_0[j] + fr * tg_xxyz_xxyyy_1[j] + fl1_fx * (tg_xyz_xxyyy_0[j] - tg_xyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyy_1[j];

                    tg_xxxyz_xxyyz_0[j] = pb_x * tg_xxyz_xxyyz_0[j] + fr * tg_xxyz_xxyyz_1[j] + fl1_fx * (tg_xyz_xxyyz_0[j] - tg_xyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyz_1[j];

                    tg_xxxyz_xxyzz_0[j] = pb_x * tg_xxyz_xxyzz_0[j] + fr * tg_xxyz_xxyzz_1[j] + fl1_fx * (tg_xyz_xxyzz_0[j] - tg_xyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyzz_1[j];

                    tg_xxxyz_xxzzz_0[j] = pb_x * tg_xxyz_xxzzz_0[j] + fr * tg_xxyz_xxzzz_1[j] + fl1_fx * (tg_xyz_xxzzz_0[j] - tg_xyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xzzz_1[j];

                    tg_xxxyz_xyyyy_0[j] = pb_x * tg_xxyz_xyyyy_0[j] + fr * tg_xxyz_xyyyy_1[j] + fl1_fx * (tg_xyz_xyyyy_0[j] - tg_xyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyy_1[j];

                    tg_xxxyz_xyyyz_0[j] = pb_x * tg_xxyz_xyyyz_0[j] + fr * tg_xxyz_xyyyz_1[j] + fl1_fx * (tg_xyz_xyyyz_0[j] - tg_xyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyz_1[j];

                    tg_xxxyz_xyyzz_0[j] = pb_x * tg_xxyz_xyyzz_0[j] + fr * tg_xxyz_xyyzz_1[j] + fl1_fx * (tg_xyz_xyyzz_0[j] - tg_xyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyzz_1[j];

                    tg_xxxyz_xyzzz_0[j] = pb_x * tg_xxyz_xyzzz_0[j] + fr * tg_xxyz_xyzzz_1[j] + fl1_fx * (tg_xyz_xyzzz_0[j] - tg_xyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yzzz_1[j];

                    tg_xxxyz_xzzzz_0[j] = pb_x * tg_xxyz_xzzzz_0[j] + fr * tg_xxyz_xzzzz_1[j] + fl1_fx * (tg_xyz_xzzzz_0[j] - tg_xyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_zzzz_1[j];

                    tg_xxxyz_yyyyy_0[j] = pb_x * tg_xxyz_yyyyy_0[j] + fr * tg_xxyz_yyyyy_1[j] + fl1_fx * (tg_xyz_yyyyy_0[j] - tg_xyz_yyyyy_1[j] * fl1_fza);

                    tg_xxxyz_yyyyz_0[j] = pb_x * tg_xxyz_yyyyz_0[j] + fr * tg_xxyz_yyyyz_1[j] + fl1_fx * (tg_xyz_yyyyz_0[j] - tg_xyz_yyyyz_1[j] * fl1_fza);

                    tg_xxxyz_yyyzz_0[j] = pb_x * tg_xxyz_yyyzz_0[j] + fr * tg_xxyz_yyyzz_1[j] + fl1_fx * (tg_xyz_yyyzz_0[j] - tg_xyz_yyyzz_1[j] * fl1_fza);

                    tg_xxxyz_yyzzz_0[j] = pb_x * tg_xxyz_yyzzz_0[j] + fr * tg_xxyz_yyzzz_1[j] + fl1_fx * (tg_xyz_yyzzz_0[j] - tg_xyz_yyzzz_1[j] * fl1_fza);

                    tg_xxxyz_yzzzz_0[j] = pb_x * tg_xxyz_yzzzz_0[j] + fr * tg_xxyz_yzzzz_1[j] + fl1_fx * (tg_xyz_yzzzz_0[j] - tg_xyz_yzzzz_1[j] * fl1_fza);

                    tg_xxxyz_zzzzz_0[j] = pb_x * tg_xxyz_zzzzz_0[j] + fr * tg_xxyz_zzzzz_1[j] + fl1_fx * (tg_xyz_zzzzz_0[j] - tg_xyz_zzzzz_1[j] * fl1_fza);

                    tg_xxxzz_xxxxx_0[j] = pb_x * tg_xxzz_xxxxx_0[j] + fr * tg_xxzz_xxxxx_1[j] + fl1_fx * (tg_xzz_xxxxx_0[j] - tg_xzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzz_xxxx_1[j];

                    tg_xxxzz_xxxxy_0[j] = pb_x * tg_xxzz_xxxxy_0[j] + fr * tg_xxzz_xxxxy_1[j] + fl1_fx * (tg_xzz_xxxxy_0[j] - tg_xzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxy_1[j];

                    tg_xxxzz_xxxxz_0[j] = pb_x * tg_xxzz_xxxxz_0[j] + fr * tg_xxzz_xxxxz_1[j] + fl1_fx * (tg_xzz_xxxxz_0[j] - tg_xzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxz_1[j];

                    tg_xxxzz_xxxyy_0[j] = pb_x * tg_xxzz_xxxyy_0[j] + fr * tg_xxzz_xxxyy_1[j] + fl1_fx * (tg_xzz_xxxyy_0[j] - tg_xzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyy_1[j];

                    tg_xxxzz_xxxyz_0[j] = pb_x * tg_xxzz_xxxyz_0[j] + fr * tg_xxzz_xxxyz_1[j] + fl1_fx * (tg_xzz_xxxyz_0[j] - tg_xzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyz_1[j];

                    tg_xxxzz_xxxzz_0[j] = pb_x * tg_xxzz_xxxzz_0[j] + fr * tg_xxzz_xxxzz_1[j] + fl1_fx * (tg_xzz_xxxzz_0[j] - tg_xzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxzz_1[j];

                    tg_xxxzz_xxyyy_0[j] = pb_x * tg_xxzz_xxyyy_0[j] + fr * tg_xxzz_xxyyy_1[j] + fl1_fx * (tg_xzz_xxyyy_0[j] - tg_xzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyy_1[j];

                    tg_xxxzz_xxyyz_0[j] = pb_x * tg_xxzz_xxyyz_0[j] + fr * tg_xxzz_xxyyz_1[j] + fl1_fx * (tg_xzz_xxyyz_0[j] - tg_xzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyz_1[j];

                    tg_xxxzz_xxyzz_0[j] = pb_x * tg_xxzz_xxyzz_0[j] + fr * tg_xxzz_xxyzz_1[j] + fl1_fx * (tg_xzz_xxyzz_0[j] - tg_xzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyzz_1[j];

                    tg_xxxzz_xxzzz_0[j] = pb_x * tg_xxzz_xxzzz_0[j] + fr * tg_xxzz_xxzzz_1[j] + fl1_fx * (tg_xzz_xxzzz_0[j] - tg_xzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xzzz_1[j];

                    tg_xxxzz_xyyyy_0[j] = pb_x * tg_xxzz_xyyyy_0[j] + fr * tg_xxzz_xyyyy_1[j] + fl1_fx * (tg_xzz_xyyyy_0[j] - tg_xzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyy_1[j];

                    tg_xxxzz_xyyyz_0[j] = pb_x * tg_xxzz_xyyyz_0[j] + fr * tg_xxzz_xyyyz_1[j] + fl1_fx * (tg_xzz_xyyyz_0[j] - tg_xzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyz_1[j];

                    tg_xxxzz_xyyzz_0[j] = pb_x * tg_xxzz_xyyzz_0[j] + fr * tg_xxzz_xyyzz_1[j] + fl1_fx * (tg_xzz_xyyzz_0[j] - tg_xzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyzz_1[j];

                    tg_xxxzz_xyzzz_0[j] = pb_x * tg_xxzz_xyzzz_0[j] + fr * tg_xxzz_xyzzz_1[j] + fl1_fx * (tg_xzz_xyzzz_0[j] - tg_xzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yzzz_1[j];

                    tg_xxxzz_xzzzz_0[j] = pb_x * tg_xxzz_xzzzz_0[j] + fr * tg_xxzz_xzzzz_1[j] + fl1_fx * (tg_xzz_xzzzz_0[j] - tg_xzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_zzzz_1[j];

                    tg_xxxzz_yyyyy_0[j] = pb_x * tg_xxzz_yyyyy_0[j] + fr * tg_xxzz_yyyyy_1[j] + fl1_fx * (tg_xzz_yyyyy_0[j] - tg_xzz_yyyyy_1[j] * fl1_fza);

                    tg_xxxzz_yyyyz_0[j] = pb_x * tg_xxzz_yyyyz_0[j] + fr * tg_xxzz_yyyyz_1[j] + fl1_fx * (tg_xzz_yyyyz_0[j] - tg_xzz_yyyyz_1[j] * fl1_fza);

                    tg_xxxzz_yyyzz_0[j] = pb_x * tg_xxzz_yyyzz_0[j] + fr * tg_xxzz_yyyzz_1[j] + fl1_fx * (tg_xzz_yyyzz_0[j] - tg_xzz_yyyzz_1[j] * fl1_fza);

                    tg_xxxzz_yyzzz_0[j] = pb_x * tg_xxzz_yyzzz_0[j] + fr * tg_xxzz_yyzzz_1[j] + fl1_fx * (tg_xzz_yyzzz_0[j] - tg_xzz_yyzzz_1[j] * fl1_fza);

                    tg_xxxzz_yzzzz_0[j] = pb_x * tg_xxzz_yzzzz_0[j] + fr * tg_xxzz_yzzzz_1[j] + fl1_fx * (tg_xzz_yzzzz_0[j] - tg_xzz_yzzzz_1[j] * fl1_fza);

                    tg_xxxzz_zzzzz_0[j] = pb_x * tg_xxzz_zzzzz_0[j] + fr * tg_xxzz_zzzzz_1[j] + fl1_fx * (tg_xzz_zzzzz_0[j] - tg_xzz_zzzzz_1[j] * fl1_fza);

                    tg_xxyyy_xxxxx_0[j] = pb_x * tg_xyyy_xxxxx_0[j] + fr * tg_xyyy_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxx_0[j] - tg_yyy_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyy_xxxx_1[j];

                    tg_xxyyy_xxxxy_0[j] = pb_x * tg_xyyy_xxxxy_0[j] + fr * tg_xyyy_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxy_0[j] - tg_yyy_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxy_1[j];

                    tg_xxyyy_xxxxz_0[j] = pb_x * tg_xyyy_xxxxz_0[j] + fr * tg_xyyy_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxz_0[j] - tg_yyy_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxz_1[j];

                    tg_xxyyy_xxxyy_0[j] = pb_x * tg_xyyy_xxxyy_0[j] + fr * tg_xyyy_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyy_0[j] - tg_yyy_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyy_1[j];

                    tg_xxyyy_xxxyz_0[j] = pb_x * tg_xyyy_xxxyz_0[j] + fr * tg_xyyy_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyz_0[j] - tg_yyy_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyz_1[j];

                    tg_xxyyy_xxxzz_0[j] = pb_x * tg_xyyy_xxxzz_0[j] + fr * tg_xyyy_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxzz_0[j] - tg_yyy_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxzz_1[j];

                    tg_xxyyy_xxyyy_0[j] = pb_x * tg_xyyy_xxyyy_0[j] + fr * tg_xyyy_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyy_0[j] - tg_yyy_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyy_1[j];

                    tg_xxyyy_xxyyz_0[j] = pb_x * tg_xyyy_xxyyz_0[j] + fr * tg_xyyy_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyz_0[j] - tg_yyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyz_1[j];

                    tg_xxyyy_xxyzz_0[j] = pb_x * tg_xyyy_xxyzz_0[j] + fr * tg_xyyy_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyzz_0[j] - tg_yyy_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyzz_1[j];

                    tg_xxyyy_xxzzz_0[j] = pb_x * tg_xyyy_xxzzz_0[j] + fr * tg_xyyy_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxzzz_0[j] - tg_yyy_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xzzz_1[j];

                    tg_xxyyy_xyyyy_0[j] = pb_x * tg_xyyy_xyyyy_0[j] + fr * tg_xyyy_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyy_0[j] - tg_yyy_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyy_1[j];

                    tg_xxyyy_xyyyz_0[j] = pb_x * tg_xyyy_xyyyz_0[j] + fr * tg_xyyy_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyz_0[j] - tg_yyy_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyz_1[j];

                    tg_xxyyy_xyyzz_0[j] = pb_x * tg_xyyy_xyyzz_0[j] + fr * tg_xyyy_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyzz_0[j] - tg_yyy_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyzz_1[j];

                    tg_xxyyy_xyzzz_0[j] = pb_x * tg_xyyy_xyzzz_0[j] + fr * tg_xyyy_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyzzz_0[j] - tg_yyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yzzz_1[j];

                    tg_xxyyy_xzzzz_0[j] = pb_x * tg_xyyy_xzzzz_0[j] + fr * tg_xyyy_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xzzzz_0[j] - tg_yyy_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_zzzz_1[j];

                    tg_xxyyy_yyyyy_0[j] = pb_x * tg_xyyy_yyyyy_0[j] + fr * tg_xyyy_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyy_0[j] - tg_yyy_yyyyy_1[j] * fl1_fza);

                    tg_xxyyy_yyyyz_0[j] = pb_x * tg_xyyy_yyyyz_0[j] + fr * tg_xyyy_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyz_0[j] - tg_yyy_yyyyz_1[j] * fl1_fza);

                    tg_xxyyy_yyyzz_0[j] = pb_x * tg_xyyy_yyyzz_0[j] + fr * tg_xyyy_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyzz_0[j] - tg_yyy_yyyzz_1[j] * fl1_fza);

                    tg_xxyyy_yyzzz_0[j] = pb_x * tg_xyyy_yyzzz_0[j] + fr * tg_xyyy_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyzzz_0[j] - tg_yyy_yyzzz_1[j] * fl1_fza);

                    tg_xxyyy_yzzzz_0[j] = pb_x * tg_xyyy_yzzzz_0[j] + fr * tg_xyyy_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yzzzz_0[j] - tg_yyy_yzzzz_1[j] * fl1_fza);

                    tg_xxyyy_zzzzz_0[j] = pb_x * tg_xyyy_zzzzz_0[j] + fr * tg_xyyy_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_zzzzz_0[j] - tg_yyy_zzzzz_1[j] * fl1_fza);

                    tg_xxyyz_xxxxx_0[j] = pb_x * tg_xyyz_xxxxx_0[j] + fr * tg_xyyz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxx_0[j] - tg_yyz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyz_xxxx_1[j];

                    tg_xxyyz_xxxxy_0[j] = pb_x * tg_xyyz_xxxxy_0[j] + fr * tg_xyyz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxy_0[j] - tg_yyz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxy_1[j];

                    tg_xxyyz_xxxxz_0[j] = pb_x * tg_xyyz_xxxxz_0[j] + fr * tg_xyyz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxz_0[j] - tg_yyz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxz_1[j];

                    tg_xxyyz_xxxyy_0[j] = pb_x * tg_xyyz_xxxyy_0[j] + fr * tg_xyyz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyy_0[j] - tg_yyz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyy_1[j];

                    tg_xxyyz_xxxyz_0[j] = pb_x * tg_xyyz_xxxyz_0[j] + fr * tg_xyyz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyz_0[j] - tg_yyz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyz_1[j];

                    tg_xxyyz_xxxzz_0[j] = pb_x * tg_xyyz_xxxzz_0[j] + fr * tg_xyyz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxzz_0[j] - tg_yyz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxzz_1[j];

                    tg_xxyyz_xxyyy_0[j] = pb_x * tg_xyyz_xxyyy_0[j] + fr * tg_xyyz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyy_0[j] - tg_yyz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyy_1[j];

                    tg_xxyyz_xxyyz_0[j] = pb_x * tg_xyyz_xxyyz_0[j] + fr * tg_xyyz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyz_0[j] - tg_yyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyz_1[j];

                    tg_xxyyz_xxyzz_0[j] = pb_x * tg_xyyz_xxyzz_0[j] + fr * tg_xyyz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyzz_0[j] - tg_yyz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyzz_1[j];

                    tg_xxyyz_xxzzz_0[j] = pb_x * tg_xyyz_xxzzz_0[j] + fr * tg_xyyz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxzzz_0[j] - tg_yyz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xzzz_1[j];

                    tg_xxyyz_xyyyy_0[j] = pb_x * tg_xyyz_xyyyy_0[j] + fr * tg_xyyz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyy_0[j] - tg_yyz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyy_1[j];

                    tg_xxyyz_xyyyz_0[j] = pb_x * tg_xyyz_xyyyz_0[j] + fr * tg_xyyz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyz_0[j] - tg_yyz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyz_1[j];

                    tg_xxyyz_xyyzz_0[j] = pb_x * tg_xyyz_xyyzz_0[j] + fr * tg_xyyz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyzz_0[j] - tg_yyz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyzz_1[j];

                    tg_xxyyz_xyzzz_0[j] = pb_x * tg_xyyz_xyzzz_0[j] + fr * tg_xyyz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyzzz_0[j] - tg_yyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yzzz_1[j];

                    tg_xxyyz_xzzzz_0[j] = pb_x * tg_xyyz_xzzzz_0[j] + fr * tg_xyyz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xzzzz_0[j] - tg_yyz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_zzzz_1[j];

                    tg_xxyyz_yyyyy_0[j] = pb_x * tg_xyyz_yyyyy_0[j] + fr * tg_xyyz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyy_0[j] - tg_yyz_yyyyy_1[j] * fl1_fza);

                    tg_xxyyz_yyyyz_0[j] = pb_x * tg_xyyz_yyyyz_0[j] + fr * tg_xyyz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyz_0[j] - tg_yyz_yyyyz_1[j] * fl1_fza);

                    tg_xxyyz_yyyzz_0[j] = pb_x * tg_xyyz_yyyzz_0[j] + fr * tg_xyyz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyzz_0[j] - tg_yyz_yyyzz_1[j] * fl1_fza);

                    tg_xxyyz_yyzzz_0[j] = pb_x * tg_xyyz_yyzzz_0[j] + fr * tg_xyyz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyzzz_0[j] - tg_yyz_yyzzz_1[j] * fl1_fza);

                    tg_xxyyz_yzzzz_0[j] = pb_x * tg_xyyz_yzzzz_0[j] + fr * tg_xyyz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yzzzz_0[j] - tg_yyz_yzzzz_1[j] * fl1_fza);

                    tg_xxyyz_zzzzz_0[j] = pb_x * tg_xyyz_zzzzz_0[j] + fr * tg_xyyz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_zzzzz_0[j] - tg_yyz_zzzzz_1[j] * fl1_fza);

                    tg_xxyzz_xxxxx_0[j] = pb_x * tg_xyzz_xxxxx_0[j] + fr * tg_xyzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxx_0[j] - tg_yzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzz_xxxx_1[j];

                    tg_xxyzz_xxxxy_0[j] = pb_x * tg_xyzz_xxxxy_0[j] + fr * tg_xyzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxy_0[j] - tg_yzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxy_1[j];

                    tg_xxyzz_xxxxz_0[j] = pb_x * tg_xyzz_xxxxz_0[j] + fr * tg_xyzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxz_0[j] - tg_yzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxz_1[j];

                    tg_xxyzz_xxxyy_0[j] = pb_x * tg_xyzz_xxxyy_0[j] + fr * tg_xyzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyy_0[j] - tg_yzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyy_1[j];

                    tg_xxyzz_xxxyz_0[j] = pb_x * tg_xyzz_xxxyz_0[j] + fr * tg_xyzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyz_0[j] - tg_yzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyz_1[j];

                    tg_xxyzz_xxxzz_0[j] = pb_x * tg_xyzz_xxxzz_0[j] + fr * tg_xyzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxzz_0[j] - tg_yzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxzz_1[j];

                    tg_xxyzz_xxyyy_0[j] = pb_x * tg_xyzz_xxyyy_0[j] + fr * tg_xyzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyy_0[j] - tg_yzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyy_1[j];

                    tg_xxyzz_xxyyz_0[j] = pb_x * tg_xyzz_xxyyz_0[j] + fr * tg_xyzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyz_0[j] - tg_yzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyz_1[j];

                    tg_xxyzz_xxyzz_0[j] = pb_x * tg_xyzz_xxyzz_0[j] + fr * tg_xyzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyzz_0[j] - tg_yzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSH_177_265(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (177,265)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_xyzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 129); 

                auto tg_xyzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 130); 

                auto tg_xyzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 131); 

                auto tg_xyzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 132); 

                auto tg_xyzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 133); 

                auto tg_xyzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 134); 

                auto tg_xzzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 135); 

                auto tg_xzzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 136); 

                auto tg_xzzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 137); 

                auto tg_xzzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 138); 

                auto tg_xzzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 139); 

                auto tg_xzzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 140); 

                auto tg_xzzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 141); 

                auto tg_xzzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 142); 

                auto tg_xzzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 143); 

                auto tg_xzzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 144); 

                auto tg_xzzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 145); 

                auto tg_xzzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 146); 

                auto tg_xzzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 147); 

                auto tg_xzzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 148); 

                auto tg_xzzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 149); 

                auto tg_yyyy_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 150); 

                auto tg_yyyy_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 151); 

                auto tg_yyyy_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 152); 

                auto tg_yyyy_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 153); 

                auto tg_yyyy_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 154); 

                auto tg_yyyy_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 155); 

                auto tg_yyyy_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 156); 

                auto tg_yyyy_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 157); 

                auto tg_yyyy_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 158); 

                auto tg_yyyy_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 159); 

                auto tg_yyyy_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 160); 

                auto tg_yyyy_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 161); 

                auto tg_yyyy_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 162); 

                auto tg_yyyy_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 163); 

                auto tg_yyyy_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 164); 

                auto tg_yyyz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 165); 

                auto tg_yyyz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 166); 

                auto tg_yyyz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 167); 

                auto tg_yyyz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 168); 

                auto tg_yyyz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 169); 

                auto tg_yyyz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 170); 

                auto tg_yyyz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 171); 

                auto tg_yyyz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 172); 

                auto tg_yyyz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 173); 

                auto tg_yyyz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 174); 

                auto tg_yyyz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 175); 

                auto tg_yyyz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 176); 

                auto tg_yyyz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 177); 

                auto tg_yyyz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 178); 

                auto tg_yyyz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 179); 

                auto tg_yyzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 180); 

                auto tg_yyzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 181); 

                auto tg_yyzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 182); 

                auto tg_yyzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 183); 

                auto tg_yyzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 184); 

                auto tg_yyzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 185); 

                auto tg_yyzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 186); 

                auto tg_yyzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 187); 

                auto tg_yyzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 188); 

                auto tg_yyzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 189); 

                auto tg_yyzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 190); 

                auto tg_yyzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 191); 

                auto tg_yyzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 192); 

                // set up pointers to integrals

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

                // Batch of Integrals (177,265)

                #pragma omp simd aligned(fxn, fza, tg_xxyzz_xxzzz_0, tg_xxyzz_xyyyy_0, tg_xxyzz_xyyyz_0, \
                                         tg_xxyzz_xyyzz_0, tg_xxyzz_xyzzz_0, tg_xxyzz_xzzzz_0, tg_xxyzz_yyyyy_0, \
                                         tg_xxyzz_yyyyz_0, tg_xxyzz_yyyzz_0, tg_xxyzz_yyzzz_0, tg_xxyzz_yzzzz_0, \
                                         tg_xxyzz_zzzzz_0, tg_xxzzz_xxxxx_0, tg_xxzzz_xxxxy_0, tg_xxzzz_xxxxz_0, \
                                         tg_xxzzz_xxxyy_0, tg_xxzzz_xxxyz_0, tg_xxzzz_xxxzz_0, tg_xxzzz_xxyyy_0, \
                                         tg_xxzzz_xxyyz_0, tg_xxzzz_xxyzz_0, tg_xxzzz_xxzzz_0, tg_xxzzz_xyyyy_0, \
                                         tg_xxzzz_xyyyz_0, tg_xxzzz_xyyzz_0, tg_xxzzz_xyzzz_0, tg_xxzzz_xzzzz_0, \
                                         tg_xxzzz_yyyyy_0, tg_xxzzz_yyyyz_0, tg_xxzzz_yyyzz_0, tg_xxzzz_yyzzz_0, \
                                         tg_xxzzz_yzzzz_0, tg_xxzzz_zzzzz_0, tg_xyyyy_xxxxx_0, tg_xyyyy_xxxxy_0, \
                                         tg_xyyyy_xxxxz_0, tg_xyyyy_xxxyy_0, tg_xyyyy_xxxyz_0, tg_xyyyy_xxxzz_0, \
                                         tg_xyyyy_xxyyy_0, tg_xyyyy_xxyyz_0, tg_xyyyy_xxyzz_0, tg_xyyyy_xxzzz_0, \
                                         tg_xyyyy_xyyyy_0, tg_xyyyy_xyyyz_0, tg_xyyyy_xyyzz_0, tg_xyyyy_xyzzz_0, \
                                         tg_xyyyy_xzzzz_0, tg_xyyyy_yyyyy_0, tg_xyyyy_yyyyz_0, tg_xyyyy_yyyzz_0, \
                                         tg_xyyyy_yyzzz_0, tg_xyyyy_yzzzz_0, tg_xyyyy_zzzzz_0, tg_xyyyz_xxxxx_0, \
                                         tg_xyyyz_xxxxy_0, tg_xyyyz_xxxxz_0, tg_xyyyz_xxxyy_0, tg_xyyyz_xxxyz_0, \
                                         tg_xyyyz_xxxzz_0, tg_xyyyz_xxyyy_0, tg_xyyyz_xxyyz_0, tg_xyyyz_xxyzz_0, \
                                         tg_xyyyz_xxzzz_0, tg_xyyyz_xyyyy_0, tg_xyyyz_xyyyz_0, tg_xyyyz_xyyzz_0, \
                                         tg_xyyyz_xyzzz_0, tg_xyyyz_xzzzz_0, tg_xyyyz_yyyyy_0, tg_xyyyz_yyyyz_0, \
                                         tg_xyyyz_yyyzz_0, tg_xyyyz_yyzzz_0, tg_xyyyz_yzzzz_0, tg_xyyyz_zzzzz_0, \
                                         tg_xyyzz_xxxxx_0, tg_xyyzz_xxxxy_0, tg_xyyzz_xxxxz_0, tg_xyyzz_xxxyy_0, \
                                         tg_xyyzz_xxxyz_0, tg_xyyzz_xxxzz_0, tg_xyyzz_xxyyy_0, tg_xyyzz_xxyyz_0, \
                                         tg_xyyzz_xxyzz_0, tg_xyyzz_xxzzz_0, tg_xyyzz_xyyyy_0, tg_xyyzz_xyyyz_0, \
                                         tg_xyyzz_xyyzz_0, tg_xyzz_xxzzz_0, tg_xyzz_xxzzz_1, tg_xyzz_xyyyy_0, tg_xyzz_xyyyy_1, \
                                         tg_xyzz_xyyyz_0, tg_xyzz_xyyyz_1, tg_xyzz_xyyzz_0, tg_xyzz_xyyzz_1, tg_xyzz_xyzzz_0, \
                                         tg_xyzz_xyzzz_1, tg_xyzz_xzzz_1, tg_xyzz_xzzzz_0, tg_xyzz_xzzzz_1, tg_xyzz_yyyy_1, \
                                         tg_xyzz_yyyyy_0, tg_xyzz_yyyyy_1, tg_xyzz_yyyyz_0, tg_xyzz_yyyyz_1, tg_xyzz_yyyz_1, \
                                         tg_xyzz_yyyzz_0, tg_xyzz_yyyzz_1, tg_xyzz_yyzz_1, tg_xyzz_yyzzz_0, tg_xyzz_yyzzz_1, \
                                         tg_xyzz_yzzz_1, tg_xyzz_yzzzz_0, tg_xyzz_yzzzz_1, tg_xyzz_zzzz_1, tg_xyzz_zzzzz_0, \
                                         tg_xyzz_zzzzz_1, tg_xzzz_xxxx_1, tg_xzzz_xxxxx_0, tg_xzzz_xxxxx_1, tg_xzzz_xxxxy_0, \
                                         tg_xzzz_xxxxy_1, tg_xzzz_xxxxz_0, tg_xzzz_xxxxz_1, tg_xzzz_xxxy_1, tg_xzzz_xxxyy_0, \
                                         tg_xzzz_xxxyy_1, tg_xzzz_xxxyz_0, tg_xzzz_xxxyz_1, tg_xzzz_xxxz_1, tg_xzzz_xxxzz_0, \
                                         tg_xzzz_xxxzz_1, tg_xzzz_xxyy_1, tg_xzzz_xxyyy_0, tg_xzzz_xxyyy_1, tg_xzzz_xxyyz_0, \
                                         tg_xzzz_xxyyz_1, tg_xzzz_xxyz_1, tg_xzzz_xxyzz_0, tg_xzzz_xxyzz_1, tg_xzzz_xxzz_1, \
                                         tg_xzzz_xxzzz_0, tg_xzzz_xxzzz_1, tg_xzzz_xyyy_1, tg_xzzz_xyyyy_0, tg_xzzz_xyyyy_1, \
                                         tg_xzzz_xyyyz_0, tg_xzzz_xyyyz_1, tg_xzzz_xyyz_1, tg_xzzz_xyyzz_0, tg_xzzz_xyyzz_1, \
                                         tg_xzzz_xyzz_1, tg_xzzz_xyzzz_0, tg_xzzz_xyzzz_1, tg_xzzz_xzzz_1, tg_xzzz_xzzzz_0, \
                                         tg_xzzz_xzzzz_1, tg_xzzz_yyyy_1, tg_xzzz_yyyyy_0, tg_xzzz_yyyyy_1, tg_xzzz_yyyyz_0, \
                                         tg_xzzz_yyyyz_1, tg_xzzz_yyyz_1, tg_xzzz_yyyzz_0, tg_xzzz_yyyzz_1, tg_xzzz_yyzz_1, \
                                         tg_xzzz_yyzzz_0, tg_xzzz_yyzzz_1, tg_xzzz_yzzz_1, tg_xzzz_yzzzz_0, tg_xzzz_yzzzz_1, \
                                         tg_xzzz_zzzz_1, tg_xzzz_zzzzz_0, tg_xzzz_zzzzz_1, tg_yyyy_xxxx_1, tg_yyyy_xxxxx_0, \
                                         tg_yyyy_xxxxx_1, tg_yyyy_xxxxy_0, tg_yyyy_xxxxy_1, tg_yyyy_xxxxz_0, tg_yyyy_xxxxz_1, \
                                         tg_yyyy_xxxy_1, tg_yyyy_xxxyy_0, tg_yyyy_xxxyy_1, tg_yyyy_xxxyz_0, tg_yyyy_xxxyz_1, \
                                         tg_yyyy_xxxz_1, tg_yyyy_xxxzz_0, tg_yyyy_xxxzz_1, tg_yyyy_xxyy_1, tg_yyyy_xxyyy_0, \
                                         tg_yyyy_xxyyy_1, tg_yyyy_xxyyz_0, tg_yyyy_xxyyz_1, tg_yyyy_xxyz_1, tg_yyyy_xxyzz_0, \
                                         tg_yyyy_xxyzz_1, tg_yyyy_xxzz_1, tg_yyyy_xxzzz_0, tg_yyyy_xxzzz_1, tg_yyyy_xyyy_1, \
                                         tg_yyyy_xyyyy_0, tg_yyyy_xyyyy_1, tg_yyyy_xyyyz_0, tg_yyyy_xyyyz_1, tg_yyyy_xyyz_1, \
                                         tg_yyyy_xyyzz_0, tg_yyyy_xyyzz_1, tg_yyyy_xyzz_1, tg_yyyy_xyzzz_0, tg_yyyy_xyzzz_1, \
                                         tg_yyyy_xzzz_1, tg_yyyy_xzzzz_0, tg_yyyy_xzzzz_1, tg_yyyy_yyyy_1, tg_yyyy_yyyyy_0, \
                                         tg_yyyy_yyyyy_1, tg_yyyy_yyyyz_0, tg_yyyy_yyyyz_1, tg_yyyy_yyyz_1, tg_yyyy_yyyzz_0, \
                                         tg_yyyy_yyyzz_1, tg_yyyy_yyzz_1, tg_yyyy_yyzzz_0, tg_yyyy_yyzzz_1, tg_yyyy_yzzz_1, \
                                         tg_yyyy_yzzzz_0, tg_yyyy_yzzzz_1, tg_yyyy_zzzz_1, tg_yyyy_zzzzz_0, tg_yyyy_zzzzz_1, \
                                         tg_yyyz_xxxx_1, tg_yyyz_xxxxx_0, tg_yyyz_xxxxx_1, tg_yyyz_xxxxy_0, tg_yyyz_xxxxy_1, \
                                         tg_yyyz_xxxxz_0, tg_yyyz_xxxxz_1, tg_yyyz_xxxy_1, tg_yyyz_xxxyy_0, tg_yyyz_xxxyy_1, \
                                         tg_yyyz_xxxyz_0, tg_yyyz_xxxyz_1, tg_yyyz_xxxz_1, tg_yyyz_xxxzz_0, tg_yyyz_xxxzz_1, \
                                         tg_yyyz_xxyy_1, tg_yyyz_xxyyy_0, tg_yyyz_xxyyy_1, tg_yyyz_xxyyz_0, tg_yyyz_xxyyz_1, \
                                         tg_yyyz_xxyz_1, tg_yyyz_xxyzz_0, tg_yyyz_xxyzz_1, tg_yyyz_xxzz_1, tg_yyyz_xxzzz_0, \
                                         tg_yyyz_xxzzz_1, tg_yyyz_xyyy_1, tg_yyyz_xyyyy_0, tg_yyyz_xyyyy_1, tg_yyyz_xyyyz_0, \
                                         tg_yyyz_xyyyz_1, tg_yyyz_xyyz_1, tg_yyyz_xyyzz_0, tg_yyyz_xyyzz_1, tg_yyyz_xyzz_1, \
                                         tg_yyyz_xyzzz_0, tg_yyyz_xyzzz_1, tg_yyyz_xzzz_1, tg_yyyz_xzzzz_0, tg_yyyz_xzzzz_1, \
                                         tg_yyyz_yyyy_1, tg_yyyz_yyyyy_0, tg_yyyz_yyyyy_1, tg_yyyz_yyyyz_0, tg_yyyz_yyyyz_1, \
                                         tg_yyyz_yyyz_1, tg_yyyz_yyyzz_0, tg_yyyz_yyyzz_1, tg_yyyz_yyzz_1, tg_yyyz_yyzzz_0, \
                                         tg_yyyz_yyzzz_1, tg_yyyz_yzzz_1, tg_yyyz_yzzzz_0, tg_yyyz_yzzzz_1, tg_yyyz_zzzz_1, \
                                         tg_yyyz_zzzzz_0, tg_yyyz_zzzzz_1, tg_yyzz_xxxx_1, tg_yyzz_xxxxx_0, tg_yyzz_xxxxx_1, \
                                         tg_yyzz_xxxxy_0, tg_yyzz_xxxxy_1, tg_yyzz_xxxxz_0, tg_yyzz_xxxxz_1, tg_yyzz_xxxy_1, \
                                         tg_yyzz_xxxyy_0, tg_yyzz_xxxyy_1, tg_yyzz_xxxyz_0, tg_yyzz_xxxyz_1, tg_yyzz_xxxz_1, \
                                         tg_yyzz_xxxzz_0, tg_yyzz_xxxzz_1, tg_yyzz_xxyy_1, tg_yyzz_xxyyy_0, tg_yyzz_xxyyy_1, \
                                         tg_yyzz_xxyyz_0, tg_yyzz_xxyyz_1, tg_yyzz_xxyz_1, tg_yyzz_xxyzz_0, tg_yyzz_xxyzz_1, \
                                         tg_yyzz_xxzz_1, tg_yyzz_xxzzz_0, tg_yyzz_xxzzz_1, tg_yyzz_xyyy_1, tg_yyzz_xyyyy_0, \
                                         tg_yyzz_xyyyy_1, tg_yyzz_xyyyz_0, tg_yyzz_xyyyz_1, tg_yyzz_xyyz_1, tg_yyzz_xyyzz_0, \
                                         tg_yyzz_xyyzz_1, tg_yyzz_xyzz_1, tg_yyzz_xzzz_1, tg_yyzz_yyyy_1, tg_yyzz_yyyz_1, \
                                         tg_yyzz_yyzz_1, tg_yzz_xxzzz_0, tg_yzz_xxzzz_1, tg_yzz_xyyyy_0, tg_yzz_xyyyy_1, \
                                         tg_yzz_xyyyz_0, tg_yzz_xyyyz_1, tg_yzz_xyyzz_0, tg_yzz_xyyzz_1, tg_yzz_xyzzz_0, \
                                         tg_yzz_xyzzz_1, tg_yzz_xzzzz_0, tg_yzz_xzzzz_1, tg_yzz_yyyyy_0, tg_yzz_yyyyy_1, \
                                         tg_yzz_yyyyz_0, tg_yzz_yyyyz_1, tg_yzz_yyyzz_0, tg_yzz_yyyzz_1, tg_yzz_yyzzz_0, \
                                         tg_yzz_yyzzz_1, tg_yzz_yzzzz_0, tg_yzz_yzzzz_1, tg_yzz_zzzzz_0, tg_yzz_zzzzz_1, \
                                         tg_zzz_xxxxx_0, tg_zzz_xxxxx_1, tg_zzz_xxxxy_0, tg_zzz_xxxxy_1, tg_zzz_xxxxz_0, \
                                         tg_zzz_xxxxz_1, tg_zzz_xxxyy_0, tg_zzz_xxxyy_1, tg_zzz_xxxyz_0, tg_zzz_xxxyz_1, \
                                         tg_zzz_xxxzz_0, tg_zzz_xxxzz_1, tg_zzz_xxyyy_0, tg_zzz_xxyyy_1, tg_zzz_xxyyz_0, \
                                         tg_zzz_xxyyz_1, tg_zzz_xxyzz_0, tg_zzz_xxyzz_1, tg_zzz_xxzzz_0, tg_zzz_xxzzz_1, \
                                         tg_zzz_xyyyy_0, tg_zzz_xyyyy_1, tg_zzz_xyyyz_0, tg_zzz_xyyyz_1, tg_zzz_xyyzz_0, \
                                         tg_zzz_xyyzz_1, tg_zzz_xyzzz_0, tg_zzz_xyzzz_1, tg_zzz_xzzzz_0, tg_zzz_xzzzz_1, \
                                         tg_zzz_yyyyy_0, tg_zzz_yyyyy_1, tg_zzz_yyyyz_0, tg_zzz_yyyyz_1, tg_zzz_yyyzz_0, \
                                         tg_zzz_yyyzz_1, tg_zzz_yyzzz_0, tg_zzz_yyzzz_1, tg_zzz_yzzzz_0, tg_zzz_yzzzz_1, \
                                         tg_zzz_zzzzz_0, tg_zzz_zzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyzz_xxzzz_0[j] = pb_x * tg_xyzz_xxzzz_0[j] + fr * tg_xyzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxzzz_0[j] - tg_yzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xzzz_1[j];

                    tg_xxyzz_xyyyy_0[j] = pb_x * tg_xyzz_xyyyy_0[j] + fr * tg_xyzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyy_0[j] - tg_yzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyy_1[j];

                    tg_xxyzz_xyyyz_0[j] = pb_x * tg_xyzz_xyyyz_0[j] + fr * tg_xyzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyz_0[j] - tg_yzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyz_1[j];

                    tg_xxyzz_xyyzz_0[j] = pb_x * tg_xyzz_xyyzz_0[j] + fr * tg_xyzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyzz_0[j] - tg_yzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyzz_1[j];

                    tg_xxyzz_xyzzz_0[j] = pb_x * tg_xyzz_xyzzz_0[j] + fr * tg_xyzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyzzz_0[j] - tg_yzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yzzz_1[j];

                    tg_xxyzz_xzzzz_0[j] = pb_x * tg_xyzz_xzzzz_0[j] + fr * tg_xyzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xzzzz_0[j] - tg_yzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_zzzz_1[j];

                    tg_xxyzz_yyyyy_0[j] = pb_x * tg_xyzz_yyyyy_0[j] + fr * tg_xyzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyy_0[j] - tg_yzz_yyyyy_1[j] * fl1_fza);

                    tg_xxyzz_yyyyz_0[j] = pb_x * tg_xyzz_yyyyz_0[j] + fr * tg_xyzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyz_0[j] - tg_yzz_yyyyz_1[j] * fl1_fza);

                    tg_xxyzz_yyyzz_0[j] = pb_x * tg_xyzz_yyyzz_0[j] + fr * tg_xyzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyzz_0[j] - tg_yzz_yyyzz_1[j] * fl1_fza);

                    tg_xxyzz_yyzzz_0[j] = pb_x * tg_xyzz_yyzzz_0[j] + fr * tg_xyzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyzzz_0[j] - tg_yzz_yyzzz_1[j] * fl1_fza);

                    tg_xxyzz_yzzzz_0[j] = pb_x * tg_xyzz_yzzzz_0[j] + fr * tg_xyzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yzzzz_0[j] - tg_yzz_yzzzz_1[j] * fl1_fza);

                    tg_xxyzz_zzzzz_0[j] = pb_x * tg_xyzz_zzzzz_0[j] + fr * tg_xyzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_zzzzz_0[j] - tg_yzz_zzzzz_1[j] * fl1_fza);

                    tg_xxzzz_xxxxx_0[j] = pb_x * tg_xzzz_xxxxx_0[j] + fr * tg_xzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxx_0[j] - tg_zzz_xxxxx_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzz_xxxx_1[j];

                    tg_xxzzz_xxxxy_0[j] = pb_x * tg_xzzz_xxxxy_0[j] + fr * tg_xzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxy_0[j] - tg_zzz_xxxxy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxy_1[j];

                    tg_xxzzz_xxxxz_0[j] = pb_x * tg_xzzz_xxxxz_0[j] + fr * tg_xzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxz_0[j] - tg_zzz_xxxxz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxz_1[j];

                    tg_xxzzz_xxxyy_0[j] = pb_x * tg_xzzz_xxxyy_0[j] + fr * tg_xzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyy_0[j] - tg_zzz_xxxyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyy_1[j];

                    tg_xxzzz_xxxyz_0[j] = pb_x * tg_xzzz_xxxyz_0[j] + fr * tg_xzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyz_0[j] - tg_zzz_xxxyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyz_1[j];

                    tg_xxzzz_xxxzz_0[j] = pb_x * tg_xzzz_xxxzz_0[j] + fr * tg_xzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxzz_0[j] - tg_zzz_xxxzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxzz_1[j];

                    tg_xxzzz_xxyyy_0[j] = pb_x * tg_xzzz_xxyyy_0[j] + fr * tg_xzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyy_0[j] - tg_zzz_xxyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyy_1[j];

                    tg_xxzzz_xxyyz_0[j] = pb_x * tg_xzzz_xxyyz_0[j] + fr * tg_xzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyz_0[j] - tg_zzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyz_1[j];

                    tg_xxzzz_xxyzz_0[j] = pb_x * tg_xzzz_xxyzz_0[j] + fr * tg_xzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyzz_0[j] - tg_zzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyzz_1[j];

                    tg_xxzzz_xxzzz_0[j] = pb_x * tg_xzzz_xxzzz_0[j] + fr * tg_xzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxzzz_0[j] - tg_zzz_xxzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xzzz_1[j];

                    tg_xxzzz_xyyyy_0[j] = pb_x * tg_xzzz_xyyyy_0[j] + fr * tg_xzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyy_0[j] - tg_zzz_xyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyy_1[j];

                    tg_xxzzz_xyyyz_0[j] = pb_x * tg_xzzz_xyyyz_0[j] + fr * tg_xzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyz_0[j] - tg_zzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyz_1[j];

                    tg_xxzzz_xyyzz_0[j] = pb_x * tg_xzzz_xyyzz_0[j] + fr * tg_xzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyzz_0[j] - tg_zzz_xyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyzz_1[j];

                    tg_xxzzz_xyzzz_0[j] = pb_x * tg_xzzz_xyzzz_0[j] + fr * tg_xzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyzzz_0[j] - tg_zzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yzzz_1[j];

                    tg_xxzzz_xzzzz_0[j] = pb_x * tg_xzzz_xzzzz_0[j] + fr * tg_xzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xzzzz_0[j] - tg_zzz_xzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_zzzz_1[j];

                    tg_xxzzz_yyyyy_0[j] = pb_x * tg_xzzz_yyyyy_0[j] + fr * tg_xzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyy_0[j] - tg_zzz_yyyyy_1[j] * fl1_fza);

                    tg_xxzzz_yyyyz_0[j] = pb_x * tg_xzzz_yyyyz_0[j] + fr * tg_xzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyz_0[j] - tg_zzz_yyyyz_1[j] * fl1_fza);

                    tg_xxzzz_yyyzz_0[j] = pb_x * tg_xzzz_yyyzz_0[j] + fr * tg_xzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyzz_0[j] - tg_zzz_yyyzz_1[j] * fl1_fza);

                    tg_xxzzz_yyzzz_0[j] = pb_x * tg_xzzz_yyzzz_0[j] + fr * tg_xzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyzzz_0[j] - tg_zzz_yyzzz_1[j] * fl1_fza);

                    tg_xxzzz_yzzzz_0[j] = pb_x * tg_xzzz_yzzzz_0[j] + fr * tg_xzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yzzzz_0[j] - tg_zzz_yzzzz_1[j] * fl1_fza);

                    tg_xxzzz_zzzzz_0[j] = pb_x * tg_xzzz_zzzzz_0[j] + fr * tg_xzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_zzzzz_0[j] - tg_zzz_zzzzz_1[j] * fl1_fza);

                    tg_xyyyy_xxxxx_0[j] = pb_x * tg_yyyy_xxxxx_0[j] + fr * tg_yyyy_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxx_1[j];

                    tg_xyyyy_xxxxy_0[j] = pb_x * tg_yyyy_xxxxy_0[j] + fr * tg_yyyy_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxy_1[j];

                    tg_xyyyy_xxxxz_0[j] = pb_x * tg_yyyy_xxxxz_0[j] + fr * tg_yyyy_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxz_1[j];

                    tg_xyyyy_xxxyy_0[j] = pb_x * tg_yyyy_xxxyy_0[j] + fr * tg_yyyy_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyy_1[j];

                    tg_xyyyy_xxxyz_0[j] = pb_x * tg_yyyy_xxxyz_0[j] + fr * tg_yyyy_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyz_1[j];

                    tg_xyyyy_xxxzz_0[j] = pb_x * tg_yyyy_xxxzz_0[j] + fr * tg_yyyy_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxzz_1[j];

                    tg_xyyyy_xxyyy_0[j] = pb_x * tg_yyyy_xxyyy_0[j] + fr * tg_yyyy_xxyyy_1[j] + fl1_fxn * tg_yyyy_xyyy_1[j];

                    tg_xyyyy_xxyyz_0[j] = pb_x * tg_yyyy_xxyyz_0[j] + fr * tg_yyyy_xxyyz_1[j] + fl1_fxn * tg_yyyy_xyyz_1[j];

                    tg_xyyyy_xxyzz_0[j] = pb_x * tg_yyyy_xxyzz_0[j] + fr * tg_yyyy_xxyzz_1[j] + fl1_fxn * tg_yyyy_xyzz_1[j];

                    tg_xyyyy_xxzzz_0[j] = pb_x * tg_yyyy_xxzzz_0[j] + fr * tg_yyyy_xxzzz_1[j] + fl1_fxn * tg_yyyy_xzzz_1[j];

                    tg_xyyyy_xyyyy_0[j] = pb_x * tg_yyyy_xyyyy_0[j] + fr * tg_yyyy_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyy_1[j];

                    tg_xyyyy_xyyyz_0[j] = pb_x * tg_yyyy_xyyyz_0[j] + fr * tg_yyyy_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyz_1[j];

                    tg_xyyyy_xyyzz_0[j] = pb_x * tg_yyyy_xyyzz_0[j] + fr * tg_yyyy_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyzz_1[j];

                    tg_xyyyy_xyzzz_0[j] = pb_x * tg_yyyy_xyzzz_0[j] + fr * tg_yyyy_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yzzz_1[j];

                    tg_xyyyy_xzzzz_0[j] = pb_x * tg_yyyy_xzzzz_0[j] + fr * tg_yyyy_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zzzz_1[j];

                    tg_xyyyy_yyyyy_0[j] = pb_x * tg_yyyy_yyyyy_0[j] + fr * tg_yyyy_yyyyy_1[j];

                    tg_xyyyy_yyyyz_0[j] = pb_x * tg_yyyy_yyyyz_0[j] + fr * tg_yyyy_yyyyz_1[j];

                    tg_xyyyy_yyyzz_0[j] = pb_x * tg_yyyy_yyyzz_0[j] + fr * tg_yyyy_yyyzz_1[j];

                    tg_xyyyy_yyzzz_0[j] = pb_x * tg_yyyy_yyzzz_0[j] + fr * tg_yyyy_yyzzz_1[j];

                    tg_xyyyy_yzzzz_0[j] = pb_x * tg_yyyy_yzzzz_0[j] + fr * tg_yyyy_yzzzz_1[j];

                    tg_xyyyy_zzzzz_0[j] = pb_x * tg_yyyy_zzzzz_0[j] + fr * tg_yyyy_zzzzz_1[j];

                    tg_xyyyz_xxxxx_0[j] = pb_x * tg_yyyz_xxxxx_0[j] + fr * tg_yyyz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxx_1[j];

                    tg_xyyyz_xxxxy_0[j] = pb_x * tg_yyyz_xxxxy_0[j] + fr * tg_yyyz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxy_1[j];

                    tg_xyyyz_xxxxz_0[j] = pb_x * tg_yyyz_xxxxz_0[j] + fr * tg_yyyz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxz_1[j];

                    tg_xyyyz_xxxyy_0[j] = pb_x * tg_yyyz_xxxyy_0[j] + fr * tg_yyyz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyy_1[j];

                    tg_xyyyz_xxxyz_0[j] = pb_x * tg_yyyz_xxxyz_0[j] + fr * tg_yyyz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyz_1[j];

                    tg_xyyyz_xxxzz_0[j] = pb_x * tg_yyyz_xxxzz_0[j] + fr * tg_yyyz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxzz_1[j];

                    tg_xyyyz_xxyyy_0[j] = pb_x * tg_yyyz_xxyyy_0[j] + fr * tg_yyyz_xxyyy_1[j] + fl1_fxn * tg_yyyz_xyyy_1[j];

                    tg_xyyyz_xxyyz_0[j] = pb_x * tg_yyyz_xxyyz_0[j] + fr * tg_yyyz_xxyyz_1[j] + fl1_fxn * tg_yyyz_xyyz_1[j];

                    tg_xyyyz_xxyzz_0[j] = pb_x * tg_yyyz_xxyzz_0[j] + fr * tg_yyyz_xxyzz_1[j] + fl1_fxn * tg_yyyz_xyzz_1[j];

                    tg_xyyyz_xxzzz_0[j] = pb_x * tg_yyyz_xxzzz_0[j] + fr * tg_yyyz_xxzzz_1[j] + fl1_fxn * tg_yyyz_xzzz_1[j];

                    tg_xyyyz_xyyyy_0[j] = pb_x * tg_yyyz_xyyyy_0[j] + fr * tg_yyyz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyy_1[j];

                    tg_xyyyz_xyyyz_0[j] = pb_x * tg_yyyz_xyyyz_0[j] + fr * tg_yyyz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyz_1[j];

                    tg_xyyyz_xyyzz_0[j] = pb_x * tg_yyyz_xyyzz_0[j] + fr * tg_yyyz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyzz_1[j];

                    tg_xyyyz_xyzzz_0[j] = pb_x * tg_yyyz_xyzzz_0[j] + fr * tg_yyyz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yzzz_1[j];

                    tg_xyyyz_xzzzz_0[j] = pb_x * tg_yyyz_xzzzz_0[j] + fr * tg_yyyz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zzzz_1[j];

                    tg_xyyyz_yyyyy_0[j] = pb_x * tg_yyyz_yyyyy_0[j] + fr * tg_yyyz_yyyyy_1[j];

                    tg_xyyyz_yyyyz_0[j] = pb_x * tg_yyyz_yyyyz_0[j] + fr * tg_yyyz_yyyyz_1[j];

                    tg_xyyyz_yyyzz_0[j] = pb_x * tg_yyyz_yyyzz_0[j] + fr * tg_yyyz_yyyzz_1[j];

                    tg_xyyyz_yyzzz_0[j] = pb_x * tg_yyyz_yyzzz_0[j] + fr * tg_yyyz_yyzzz_1[j];

                    tg_xyyyz_yzzzz_0[j] = pb_x * tg_yyyz_yzzzz_0[j] + fr * tg_yyyz_yzzzz_1[j];

                    tg_xyyyz_zzzzz_0[j] = pb_x * tg_yyyz_zzzzz_0[j] + fr * tg_yyyz_zzzzz_1[j];

                    tg_xyyzz_xxxxx_0[j] = pb_x * tg_yyzz_xxxxx_0[j] + fr * tg_yyzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxx_1[j];

                    tg_xyyzz_xxxxy_0[j] = pb_x * tg_yyzz_xxxxy_0[j] + fr * tg_yyzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxy_1[j];

                    tg_xyyzz_xxxxz_0[j] = pb_x * tg_yyzz_xxxxz_0[j] + fr * tg_yyzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxz_1[j];

                    tg_xyyzz_xxxyy_0[j] = pb_x * tg_yyzz_xxxyy_0[j] + fr * tg_yyzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyy_1[j];

                    tg_xyyzz_xxxyz_0[j] = pb_x * tg_yyzz_xxxyz_0[j] + fr * tg_yyzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyz_1[j];

                    tg_xyyzz_xxxzz_0[j] = pb_x * tg_yyzz_xxxzz_0[j] + fr * tg_yyzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxzz_1[j];

                    tg_xyyzz_xxyyy_0[j] = pb_x * tg_yyzz_xxyyy_0[j] + fr * tg_yyzz_xxyyy_1[j] + fl1_fxn * tg_yyzz_xyyy_1[j];

                    tg_xyyzz_xxyyz_0[j] = pb_x * tg_yyzz_xxyyz_0[j] + fr * tg_yyzz_xxyyz_1[j] + fl1_fxn * tg_yyzz_xyyz_1[j];

                    tg_xyyzz_xxyzz_0[j] = pb_x * tg_yyzz_xxyzz_0[j] + fr * tg_yyzz_xxyzz_1[j] + fl1_fxn * tg_yyzz_xyzz_1[j];

                    tg_xyyzz_xxzzz_0[j] = pb_x * tg_yyzz_xxzzz_0[j] + fr * tg_yyzz_xxzzz_1[j] + fl1_fxn * tg_yyzz_xzzz_1[j];

                    tg_xyyzz_xyyyy_0[j] = pb_x * tg_yyzz_xyyyy_0[j] + fr * tg_yyzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyy_1[j];

                    tg_xyyzz_xyyyz_0[j] = pb_x * tg_yyzz_xyyyz_0[j] + fr * tg_yyzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyz_1[j];

                    tg_xyyzz_xyyzz_0[j] = pb_x * tg_yyzz_xyyzz_0[j] + fr * tg_yyzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSH_265_353(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (265,353)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyz_xyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 158); 

                auto tg_yyz_xyyzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 159); 

                auto tg_yyz_xyzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 160); 

                auto tg_yyz_xzzzz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 161); 

                auto tg_yyz_yyyyy_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 162); 

                auto tg_yyz_yyyyz_0 = primBuffer[pidx_g_3_5_m0].data(210 * idx + 163); 

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

                auto tg_yyz_xyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 158); 

                auto tg_yyz_xyyzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 159); 

                auto tg_yyz_xyzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 160); 

                auto tg_yyz_xzzzz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 161); 

                auto tg_yyz_yyyyy_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 162); 

                auto tg_yyz_yyyyz_1 = primBuffer[pidx_g_3_5_m1].data(210 * idx + 163); 

                auto tg_yyyy_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 150); 

                auto tg_yyyy_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 151); 

                auto tg_yyyy_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 152); 

                auto tg_yyyy_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 153); 

                auto tg_yyyy_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 154); 

                auto tg_yyyy_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 155); 

                auto tg_yyyy_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 156); 

                auto tg_yyyy_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 157); 

                auto tg_yyyy_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 158); 

                auto tg_yyyy_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 159); 

                auto tg_yyyy_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 160); 

                auto tg_yyyy_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 161); 

                auto tg_yyyy_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 162); 

                auto tg_yyyy_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 163); 

                auto tg_yyyy_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 164); 

                auto tg_yyyz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 165); 

                auto tg_yyyz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 166); 

                auto tg_yyyz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 167); 

                auto tg_yyyz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 168); 

                auto tg_yyyz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 169); 

                auto tg_yyyz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 170); 

                auto tg_yyyz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 171); 

                auto tg_yyyz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 172); 

                auto tg_yyyz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 173); 

                auto tg_yyyz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 174); 

                auto tg_yyyz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 175); 

                auto tg_yyyz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 176); 

                auto tg_yyzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 193); 

                auto tg_yyzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 194); 

                auto tg_yzzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 195); 

                auto tg_yzzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 196); 

                auto tg_yzzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 197); 

                auto tg_yzzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 198); 

                auto tg_yzzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 199); 

                auto tg_yzzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 200); 

                auto tg_yzzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 201); 

                auto tg_yzzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 202); 

                auto tg_yzzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 203); 

                auto tg_yzzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 204); 

                auto tg_yzzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 205); 

                auto tg_yzzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 206); 

                auto tg_yzzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 207); 

                auto tg_yzzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 208); 

                auto tg_yzzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 209); 

                auto tg_zzzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 210); 

                auto tg_zzzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 211); 

                auto tg_zzzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 212); 

                auto tg_zzzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 213); 

                auto tg_zzzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 214); 

                auto tg_zzzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 215); 

                auto tg_zzzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 216); 

                auto tg_zzzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 217); 

                auto tg_zzzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 218); 

                auto tg_zzzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 219); 

                auto tg_zzzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 220); 

                auto tg_zzzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 221); 

                auto tg_zzzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 222); 

                auto tg_zzzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 223); 

                auto tg_zzzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 224); 

                // set up pointers to integrals

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

                // Batch of Integrals (265,353)

                #pragma omp simd aligned(fxn, fza, tg_xyyzz_xyzzz_0, tg_xyyzz_xzzzz_0, tg_xyyzz_yyyyy_0, \
                                         tg_xyyzz_yyyyz_0, tg_xyyzz_yyyzz_0, tg_xyyzz_yyzzz_0, tg_xyyzz_yzzzz_0, \
                                         tg_xyyzz_zzzzz_0, tg_xyzzz_xxxxx_0, tg_xyzzz_xxxxy_0, tg_xyzzz_xxxxz_0, \
                                         tg_xyzzz_xxxyy_0, tg_xyzzz_xxxyz_0, tg_xyzzz_xxxzz_0, tg_xyzzz_xxyyy_0, \
                                         tg_xyzzz_xxyyz_0, tg_xyzzz_xxyzz_0, tg_xyzzz_xxzzz_0, tg_xyzzz_xyyyy_0, \
                                         tg_xyzzz_xyyyz_0, tg_xyzzz_xyyzz_0, tg_xyzzz_xyzzz_0, tg_xyzzz_xzzzz_0, \
                                         tg_xyzzz_yyyyy_0, tg_xyzzz_yyyyz_0, tg_xyzzz_yyyzz_0, tg_xyzzz_yyzzz_0, \
                                         tg_xyzzz_yzzzz_0, tg_xyzzz_zzzzz_0, tg_xzzzz_xxxxx_0, tg_xzzzz_xxxxy_0, \
                                         tg_xzzzz_xxxxz_0, tg_xzzzz_xxxyy_0, tg_xzzzz_xxxyz_0, tg_xzzzz_xxxzz_0, \
                                         tg_xzzzz_xxyyy_0, tg_xzzzz_xxyyz_0, tg_xzzzz_xxyzz_0, tg_xzzzz_xxzzz_0, \
                                         tg_xzzzz_xyyyy_0, tg_xzzzz_xyyyz_0, tg_xzzzz_xyyzz_0, tg_xzzzz_xyzzz_0, \
                                         tg_xzzzz_xzzzz_0, tg_xzzzz_yyyyy_0, tg_xzzzz_yyyyz_0, tg_xzzzz_yyyzz_0, \
                                         tg_xzzzz_yyzzz_0, tg_xzzzz_yzzzz_0, tg_xzzzz_zzzzz_0, tg_yyy_xxxxx_0, tg_yyy_xxxxx_1, \
                                         tg_yyy_xxxxy_0, tg_yyy_xxxxy_1, tg_yyy_xxxxz_0, tg_yyy_xxxxz_1, tg_yyy_xxxyy_0, \
                                         tg_yyy_xxxyy_1, tg_yyy_xxxyz_0, tg_yyy_xxxyz_1, tg_yyy_xxxzz_0, tg_yyy_xxxzz_1, \
                                         tg_yyy_xxyyy_0, tg_yyy_xxyyy_1, tg_yyy_xxyyz_0, tg_yyy_xxyyz_1, tg_yyy_xxyzz_0, \
                                         tg_yyy_xxyzz_1, tg_yyy_xxzzz_0, tg_yyy_xxzzz_1, tg_yyy_xyyyy_0, tg_yyy_xyyyy_1, \
                                         tg_yyy_xyyyz_0, tg_yyy_xyyyz_1, tg_yyy_xyyzz_0, tg_yyy_xyyzz_1, tg_yyy_xyzzz_0, \
                                         tg_yyy_xyzzz_1, tg_yyy_xzzzz_0, tg_yyy_xzzzz_1, tg_yyy_yyyyy_0, tg_yyy_yyyyy_1, \
                                         tg_yyy_yyyyz_0, tg_yyy_yyyyz_1, tg_yyy_yyyzz_0, tg_yyy_yyyzz_1, tg_yyy_yyzzz_0, \
                                         tg_yyy_yyzzz_1, tg_yyy_yzzzz_0, tg_yyy_yzzzz_1, tg_yyy_zzzzz_0, tg_yyy_zzzzz_1, \
                                         tg_yyyy_xxxx_1, tg_yyyy_xxxxx_0, tg_yyyy_xxxxx_1, tg_yyyy_xxxxy_0, tg_yyyy_xxxxy_1, \
                                         tg_yyyy_xxxxz_0, tg_yyyy_xxxxz_1, tg_yyyy_xxxy_1, tg_yyyy_xxxyy_0, tg_yyyy_xxxyy_1, \
                                         tg_yyyy_xxxyz_0, tg_yyyy_xxxyz_1, tg_yyyy_xxxz_1, tg_yyyy_xxxzz_0, tg_yyyy_xxxzz_1, \
                                         tg_yyyy_xxyy_1, tg_yyyy_xxyyy_0, tg_yyyy_xxyyy_1, tg_yyyy_xxyyz_0, tg_yyyy_xxyyz_1, \
                                         tg_yyyy_xxyz_1, tg_yyyy_xxyzz_0, tg_yyyy_xxyzz_1, tg_yyyy_xxzz_1, tg_yyyy_xxzzz_0, \
                                         tg_yyyy_xxzzz_1, tg_yyyy_xyyy_1, tg_yyyy_xyyyy_0, tg_yyyy_xyyyy_1, tg_yyyy_xyyyz_0, \
                                         tg_yyyy_xyyyz_1, tg_yyyy_xyyz_1, tg_yyyy_xyyzz_0, tg_yyyy_xyyzz_1, tg_yyyy_xyzz_1, \
                                         tg_yyyy_xyzzz_0, tg_yyyy_xyzzz_1, tg_yyyy_xzzz_1, tg_yyyy_xzzzz_0, tg_yyyy_xzzzz_1, \
                                         tg_yyyy_yyyy_1, tg_yyyy_yyyyy_0, tg_yyyy_yyyyy_1, tg_yyyy_yyyyz_0, tg_yyyy_yyyyz_1, \
                                         tg_yyyy_yyyz_1, tg_yyyy_yyyzz_0, tg_yyyy_yyyzz_1, tg_yyyy_yyzz_1, tg_yyyy_yyzzz_0, \
                                         tg_yyyy_yyzzz_1, tg_yyyy_yzzz_1, tg_yyyy_yzzzz_0, tg_yyyy_yzzzz_1, tg_yyyy_zzzz_1, \
                                         tg_yyyy_zzzzz_0, tg_yyyy_zzzzz_1, tg_yyyyy_xxxxx_0, tg_yyyyy_xxxxy_0, \
                                         tg_yyyyy_xxxxz_0, tg_yyyyy_xxxyy_0, tg_yyyyy_xxxyz_0, tg_yyyyy_xxxzz_0, \
                                         tg_yyyyy_xxyyy_0, tg_yyyyy_xxyyz_0, tg_yyyyy_xxyzz_0, tg_yyyyy_xxzzz_0, \
                                         tg_yyyyy_xyyyy_0, tg_yyyyy_xyyyz_0, tg_yyyyy_xyyzz_0, tg_yyyyy_xyzzz_0, \
                                         tg_yyyyy_xzzzz_0, tg_yyyyy_yyyyy_0, tg_yyyyy_yyyyz_0, tg_yyyyy_yyyzz_0, \
                                         tg_yyyyy_yyzzz_0, tg_yyyyy_yzzzz_0, tg_yyyyy_zzzzz_0, tg_yyyyz_xxxxx_0, \
                                         tg_yyyyz_xxxxy_0, tg_yyyyz_xxxxz_0, tg_yyyyz_xxxyy_0, tg_yyyyz_xxxyz_0, \
                                         tg_yyyyz_xxxzz_0, tg_yyyyz_xxyyy_0, tg_yyyyz_xxyyz_0, tg_yyyyz_xxyzz_0, \
                                         tg_yyyyz_xxzzz_0, tg_yyyyz_xyyyy_0, tg_yyyyz_xyyyz_0, tg_yyyyz_xyyzz_0, \
                                         tg_yyyyz_xyzzz_0, tg_yyyyz_xzzzz_0, tg_yyyyz_yyyyy_0, tg_yyyyz_yyyyz_0, \
                                         tg_yyyz_xxxx_1, tg_yyyz_xxxxx_0, tg_yyyz_xxxxx_1, tg_yyyz_xxxxy_0, tg_yyyz_xxxxy_1, \
                                         tg_yyyz_xxxxz_0, tg_yyyz_xxxxz_1, tg_yyyz_xxxy_1, tg_yyyz_xxxyy_0, tg_yyyz_xxxyy_1, \
                                         tg_yyyz_xxxyz_0, tg_yyyz_xxxyz_1, tg_yyyz_xxxz_1, tg_yyyz_xxxzz_0, tg_yyyz_xxxzz_1, \
                                         tg_yyyz_xxyy_1, tg_yyyz_xxyyy_0, tg_yyyz_xxyyy_1, tg_yyyz_xxyyz_0, tg_yyyz_xxyyz_1, \
                                         tg_yyyz_xxyz_1, tg_yyyz_xxyzz_0, tg_yyyz_xxyzz_1, tg_yyyz_xxzz_1, tg_yyyz_xxzzz_0, \
                                         tg_yyyz_xxzzz_1, tg_yyyz_xyyy_1, tg_yyyz_xyyyy_0, tg_yyyz_xyyyy_1, tg_yyyz_xyyyz_0, \
                                         tg_yyyz_xyyyz_1, tg_yyyz_xyyz_1, tg_yyyz_xyyzz_0, tg_yyyz_xyyzz_1, tg_yyyz_xyzz_1, \
                                         tg_yyyz_xyzzz_0, tg_yyyz_xyzzz_1, tg_yyyz_xzzz_1, tg_yyyz_xzzzz_0, tg_yyyz_xzzzz_1, \
                                         tg_yyyz_yyyy_1, tg_yyyz_yyyyy_0, tg_yyyz_yyyyy_1, tg_yyyz_yyyyz_0, tg_yyyz_yyyyz_1, \
                                         tg_yyyz_yyyz_1, tg_yyz_xxxxx_0, tg_yyz_xxxxx_1, tg_yyz_xxxxy_0, tg_yyz_xxxxy_1, \
                                         tg_yyz_xxxxz_0, tg_yyz_xxxxz_1, tg_yyz_xxxyy_0, tg_yyz_xxxyy_1, tg_yyz_xxxyz_0, \
                                         tg_yyz_xxxyz_1, tg_yyz_xxxzz_0, tg_yyz_xxxzz_1, tg_yyz_xxyyy_0, tg_yyz_xxyyy_1, \
                                         tg_yyz_xxyyz_0, tg_yyz_xxyyz_1, tg_yyz_xxyzz_0, tg_yyz_xxyzz_1, tg_yyz_xxzzz_0, \
                                         tg_yyz_xxzzz_1, tg_yyz_xyyyy_0, tg_yyz_xyyyy_1, tg_yyz_xyyyz_0, tg_yyz_xyyyz_1, \
                                         tg_yyz_xyyzz_0, tg_yyz_xyyzz_1, tg_yyz_xyzzz_0, tg_yyz_xyzzz_1, tg_yyz_xzzzz_0, \
                                         tg_yyz_xzzzz_1, tg_yyz_yyyyy_0, tg_yyz_yyyyy_1, tg_yyz_yyyyz_0, tg_yyz_yyyyz_1, \
                                         tg_yyzz_xyzzz_0, tg_yyzz_xyzzz_1, tg_yyzz_xzzzz_0, tg_yyzz_xzzzz_1, tg_yyzz_yyyyy_0, \
                                         tg_yyzz_yyyyy_1, tg_yyzz_yyyyz_0, tg_yyzz_yyyyz_1, tg_yyzz_yyyzz_0, tg_yyzz_yyyzz_1, \
                                         tg_yyzz_yyzzz_0, tg_yyzz_yyzzz_1, tg_yyzz_yzzz_1, tg_yyzz_yzzzz_0, tg_yyzz_yzzzz_1, \
                                         tg_yyzz_zzzz_1, tg_yyzz_zzzzz_0, tg_yyzz_zzzzz_1, tg_yzzz_xxxx_1, tg_yzzz_xxxxx_0, \
                                         tg_yzzz_xxxxx_1, tg_yzzz_xxxxy_0, tg_yzzz_xxxxy_1, tg_yzzz_xxxxz_0, tg_yzzz_xxxxz_1, \
                                         tg_yzzz_xxxy_1, tg_yzzz_xxxyy_0, tg_yzzz_xxxyy_1, tg_yzzz_xxxyz_0, tg_yzzz_xxxyz_1, \
                                         tg_yzzz_xxxz_1, tg_yzzz_xxxzz_0, tg_yzzz_xxxzz_1, tg_yzzz_xxyy_1, tg_yzzz_xxyyy_0, \
                                         tg_yzzz_xxyyy_1, tg_yzzz_xxyyz_0, tg_yzzz_xxyyz_1, tg_yzzz_xxyz_1, tg_yzzz_xxyzz_0, \
                                         tg_yzzz_xxyzz_1, tg_yzzz_xxzz_1, tg_yzzz_xxzzz_0, tg_yzzz_xxzzz_1, tg_yzzz_xyyy_1, \
                                         tg_yzzz_xyyyy_0, tg_yzzz_xyyyy_1, tg_yzzz_xyyyz_0, tg_yzzz_xyyyz_1, tg_yzzz_xyyz_1, \
                                         tg_yzzz_xyyzz_0, tg_yzzz_xyyzz_1, tg_yzzz_xyzz_1, tg_yzzz_xyzzz_0, tg_yzzz_xyzzz_1, \
                                         tg_yzzz_xzzz_1, tg_yzzz_xzzzz_0, tg_yzzz_xzzzz_1, tg_yzzz_yyyy_1, tg_yzzz_yyyyy_0, \
                                         tg_yzzz_yyyyy_1, tg_yzzz_yyyyz_0, tg_yzzz_yyyyz_1, tg_yzzz_yyyz_1, tg_yzzz_yyyzz_0, \
                                         tg_yzzz_yyyzz_1, tg_yzzz_yyzz_1, tg_yzzz_yyzzz_0, tg_yzzz_yyzzz_1, tg_yzzz_yzzz_1, \
                                         tg_yzzz_yzzzz_0, tg_yzzz_yzzzz_1, tg_yzzz_zzzz_1, tg_yzzz_zzzzz_0, tg_yzzz_zzzzz_1, \
                                         tg_zzzz_xxxx_1, tg_zzzz_xxxxx_0, tg_zzzz_xxxxx_1, tg_zzzz_xxxxy_0, tg_zzzz_xxxxy_1, \
                                         tg_zzzz_xxxxz_0, tg_zzzz_xxxxz_1, tg_zzzz_xxxy_1, tg_zzzz_xxxyy_0, tg_zzzz_xxxyy_1, \
                                         tg_zzzz_xxxyz_0, tg_zzzz_xxxyz_1, tg_zzzz_xxxz_1, tg_zzzz_xxxzz_0, tg_zzzz_xxxzz_1, \
                                         tg_zzzz_xxyy_1, tg_zzzz_xxyyy_0, tg_zzzz_xxyyy_1, tg_zzzz_xxyyz_0, tg_zzzz_xxyyz_1, \
                                         tg_zzzz_xxyz_1, tg_zzzz_xxyzz_0, tg_zzzz_xxyzz_1, tg_zzzz_xxzz_1, tg_zzzz_xxzzz_0, \
                                         tg_zzzz_xxzzz_1, tg_zzzz_xyyy_1, tg_zzzz_xyyyy_0, tg_zzzz_xyyyy_1, tg_zzzz_xyyyz_0, \
                                         tg_zzzz_xyyyz_1, tg_zzzz_xyyz_1, tg_zzzz_xyyzz_0, tg_zzzz_xyyzz_1, tg_zzzz_xyzz_1, \
                                         tg_zzzz_xyzzz_0, tg_zzzz_xyzzz_1, tg_zzzz_xzzz_1, tg_zzzz_xzzzz_0, tg_zzzz_xzzzz_1, \
                                         tg_zzzz_yyyy_1, tg_zzzz_yyyyy_0, tg_zzzz_yyyyy_1, tg_zzzz_yyyyz_0, tg_zzzz_yyyyz_1, \
                                         tg_zzzz_yyyz_1, tg_zzzz_yyyzz_0, tg_zzzz_yyyzz_1, tg_zzzz_yyzz_1, tg_zzzz_yyzzz_0, \
                                         tg_zzzz_yyzzz_1, tg_zzzz_yzzz_1, tg_zzzz_yzzzz_0, tg_zzzz_yzzzz_1, tg_zzzz_zzzz_1, \
                                         tg_zzzz_zzzzz_0, tg_zzzz_zzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyyzz_xyzzz_0[j] = pb_x * tg_yyzz_xyzzz_0[j] + fr * tg_yyzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yzzz_1[j];

                    tg_xyyzz_xzzzz_0[j] = pb_x * tg_yyzz_xzzzz_0[j] + fr * tg_yyzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zzzz_1[j];

                    tg_xyyzz_yyyyy_0[j] = pb_x * tg_yyzz_yyyyy_0[j] + fr * tg_yyzz_yyyyy_1[j];

                    tg_xyyzz_yyyyz_0[j] = pb_x * tg_yyzz_yyyyz_0[j] + fr * tg_yyzz_yyyyz_1[j];

                    tg_xyyzz_yyyzz_0[j] = pb_x * tg_yyzz_yyyzz_0[j] + fr * tg_yyzz_yyyzz_1[j];

                    tg_xyyzz_yyzzz_0[j] = pb_x * tg_yyzz_yyzzz_0[j] + fr * tg_yyzz_yyzzz_1[j];

                    tg_xyyzz_yzzzz_0[j] = pb_x * tg_yyzz_yzzzz_0[j] + fr * tg_yyzz_yzzzz_1[j];

                    tg_xyyzz_zzzzz_0[j] = pb_x * tg_yyzz_zzzzz_0[j] + fr * tg_yyzz_zzzzz_1[j];

                    tg_xyzzz_xxxxx_0[j] = pb_x * tg_yzzz_xxxxx_0[j] + fr * tg_yzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxx_1[j];

                    tg_xyzzz_xxxxy_0[j] = pb_x * tg_yzzz_xxxxy_0[j] + fr * tg_yzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxy_1[j];

                    tg_xyzzz_xxxxz_0[j] = pb_x * tg_yzzz_xxxxz_0[j] + fr * tg_yzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxz_1[j];

                    tg_xyzzz_xxxyy_0[j] = pb_x * tg_yzzz_xxxyy_0[j] + fr * tg_yzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyy_1[j];

                    tg_xyzzz_xxxyz_0[j] = pb_x * tg_yzzz_xxxyz_0[j] + fr * tg_yzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyz_1[j];

                    tg_xyzzz_xxxzz_0[j] = pb_x * tg_yzzz_xxxzz_0[j] + fr * tg_yzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxzz_1[j];

                    tg_xyzzz_xxyyy_0[j] = pb_x * tg_yzzz_xxyyy_0[j] + fr * tg_yzzz_xxyyy_1[j] + fl1_fxn * tg_yzzz_xyyy_1[j];

                    tg_xyzzz_xxyyz_0[j] = pb_x * tg_yzzz_xxyyz_0[j] + fr * tg_yzzz_xxyyz_1[j] + fl1_fxn * tg_yzzz_xyyz_1[j];

                    tg_xyzzz_xxyzz_0[j] = pb_x * tg_yzzz_xxyzz_0[j] + fr * tg_yzzz_xxyzz_1[j] + fl1_fxn * tg_yzzz_xyzz_1[j];

                    tg_xyzzz_xxzzz_0[j] = pb_x * tg_yzzz_xxzzz_0[j] + fr * tg_yzzz_xxzzz_1[j] + fl1_fxn * tg_yzzz_xzzz_1[j];

                    tg_xyzzz_xyyyy_0[j] = pb_x * tg_yzzz_xyyyy_0[j] + fr * tg_yzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyy_1[j];

                    tg_xyzzz_xyyyz_0[j] = pb_x * tg_yzzz_xyyyz_0[j] + fr * tg_yzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyz_1[j];

                    tg_xyzzz_xyyzz_0[j] = pb_x * tg_yzzz_xyyzz_0[j] + fr * tg_yzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyzz_1[j];

                    tg_xyzzz_xyzzz_0[j] = pb_x * tg_yzzz_xyzzz_0[j] + fr * tg_yzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yzzz_1[j];

                    tg_xyzzz_xzzzz_0[j] = pb_x * tg_yzzz_xzzzz_0[j] + fr * tg_yzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zzzz_1[j];

                    tg_xyzzz_yyyyy_0[j] = pb_x * tg_yzzz_yyyyy_0[j] + fr * tg_yzzz_yyyyy_1[j];

                    tg_xyzzz_yyyyz_0[j] = pb_x * tg_yzzz_yyyyz_0[j] + fr * tg_yzzz_yyyyz_1[j];

                    tg_xyzzz_yyyzz_0[j] = pb_x * tg_yzzz_yyyzz_0[j] + fr * tg_yzzz_yyyzz_1[j];

                    tg_xyzzz_yyzzz_0[j] = pb_x * tg_yzzz_yyzzz_0[j] + fr * tg_yzzz_yyzzz_1[j];

                    tg_xyzzz_yzzzz_0[j] = pb_x * tg_yzzz_yzzzz_0[j] + fr * tg_yzzz_yzzzz_1[j];

                    tg_xyzzz_zzzzz_0[j] = pb_x * tg_yzzz_zzzzz_0[j] + fr * tg_yzzz_zzzzz_1[j];

                    tg_xzzzz_xxxxx_0[j] = pb_x * tg_zzzz_xxxxx_0[j] + fr * tg_zzzz_xxxxx_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxx_1[j];

                    tg_xzzzz_xxxxy_0[j] = pb_x * tg_zzzz_xxxxy_0[j] + fr * tg_zzzz_xxxxy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxy_1[j];

                    tg_xzzzz_xxxxz_0[j] = pb_x * tg_zzzz_xxxxz_0[j] + fr * tg_zzzz_xxxxz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxz_1[j];

                    tg_xzzzz_xxxyy_0[j] = pb_x * tg_zzzz_xxxyy_0[j] + fr * tg_zzzz_xxxyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyy_1[j];

                    tg_xzzzz_xxxyz_0[j] = pb_x * tg_zzzz_xxxyz_0[j] + fr * tg_zzzz_xxxyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyz_1[j];

                    tg_xzzzz_xxxzz_0[j] = pb_x * tg_zzzz_xxxzz_0[j] + fr * tg_zzzz_xxxzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxzz_1[j];

                    tg_xzzzz_xxyyy_0[j] = pb_x * tg_zzzz_xxyyy_0[j] + fr * tg_zzzz_xxyyy_1[j] + fl1_fxn * tg_zzzz_xyyy_1[j];

                    tg_xzzzz_xxyyz_0[j] = pb_x * tg_zzzz_xxyyz_0[j] + fr * tg_zzzz_xxyyz_1[j] + fl1_fxn * tg_zzzz_xyyz_1[j];

                    tg_xzzzz_xxyzz_0[j] = pb_x * tg_zzzz_xxyzz_0[j] + fr * tg_zzzz_xxyzz_1[j] + fl1_fxn * tg_zzzz_xyzz_1[j];

                    tg_xzzzz_xxzzz_0[j] = pb_x * tg_zzzz_xxzzz_0[j] + fr * tg_zzzz_xxzzz_1[j] + fl1_fxn * tg_zzzz_xzzz_1[j];

                    tg_xzzzz_xyyyy_0[j] = pb_x * tg_zzzz_xyyyy_0[j] + fr * tg_zzzz_xyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyy_1[j];

                    tg_xzzzz_xyyyz_0[j] = pb_x * tg_zzzz_xyyyz_0[j] + fr * tg_zzzz_xyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyz_1[j];

                    tg_xzzzz_xyyzz_0[j] = pb_x * tg_zzzz_xyyzz_0[j] + fr * tg_zzzz_xyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyzz_1[j];

                    tg_xzzzz_xyzzz_0[j] = pb_x * tg_zzzz_xyzzz_0[j] + fr * tg_zzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yzzz_1[j];

                    tg_xzzzz_xzzzz_0[j] = pb_x * tg_zzzz_xzzzz_0[j] + fr * tg_zzzz_xzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzz_1[j];

                    tg_xzzzz_yyyyy_0[j] = pb_x * tg_zzzz_yyyyy_0[j] + fr * tg_zzzz_yyyyy_1[j];

                    tg_xzzzz_yyyyz_0[j] = pb_x * tg_zzzz_yyyyz_0[j] + fr * tg_zzzz_yyyyz_1[j];

                    tg_xzzzz_yyyzz_0[j] = pb_x * tg_zzzz_yyyzz_0[j] + fr * tg_zzzz_yyyzz_1[j];

                    tg_xzzzz_yyzzz_0[j] = pb_x * tg_zzzz_yyzzz_0[j] + fr * tg_zzzz_yyzzz_1[j];

                    tg_xzzzz_yzzzz_0[j] = pb_x * tg_zzzz_yzzzz_0[j] + fr * tg_zzzz_yzzzz_1[j];

                    tg_xzzzz_zzzzz_0[j] = pb_x * tg_zzzz_zzzzz_0[j] + fr * tg_zzzz_zzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyy_xxxxx_0[j] = pb_y * tg_yyyy_xxxxx_0[j] + fr * tg_yyyy_xxxxx_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxx_0[j] - tg_yyy_xxxxx_1[j] * fl1_fza);

                    tg_yyyyy_xxxxy_0[j] = pb_y * tg_yyyy_xxxxy_0[j] + fr * tg_yyyy_xxxxy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxy_0[j] - tg_yyy_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxx_1[j];

                    tg_yyyyy_xxxxz_0[j] = pb_y * tg_yyyy_xxxxz_0[j] + fr * tg_yyyy_xxxxz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxz_0[j] - tg_yyy_xxxxz_1[j] * fl1_fza);

                    tg_yyyyy_xxxyy_0[j] = pb_y * tg_yyyy_xxxyy_0[j] + fr * tg_yyyy_xxxyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyy_0[j] - tg_yyy_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxxy_1[j];

                    tg_yyyyy_xxxyz_0[j] = pb_y * tg_yyyy_xxxyz_0[j] + fr * tg_yyyy_xxxyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyz_0[j] - tg_yyy_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxz_1[j];

                    tg_yyyyy_xxxzz_0[j] = pb_y * tg_yyyy_xxxzz_0[j] + fr * tg_yyyy_xxxzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxzz_0[j] - tg_yyy_xxxzz_1[j] * fl1_fza);

                    tg_yyyyy_xxyyy_0[j] = pb_y * tg_yyyy_xxyyy_0[j] + fr * tg_yyyy_xxyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyy_0[j] - tg_yyy_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xxyy_1[j];

                    tg_yyyyy_xxyyz_0[j] = pb_y * tg_yyyy_xxyyz_0[j] + fr * tg_yyyy_xxyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyz_0[j] - tg_yyy_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxyz_1[j];

                    tg_yyyyy_xxyzz_0[j] = pb_y * tg_yyyy_xxyzz_0[j] + fr * tg_yyyy_xxyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyzz_0[j] - tg_yyy_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxzz_1[j];

                    tg_yyyyy_xxzzz_0[j] = pb_y * tg_yyyy_xxzzz_0[j] + fr * tg_yyyy_xxzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxzzz_0[j] - tg_yyy_xxzzz_1[j] * fl1_fza);

                    tg_yyyyy_xyyyy_0[j] = pb_y * tg_yyyy_xyyyy_0[j] + fr * tg_yyyy_xyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyy_0[j] - tg_yyy_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_xyyy_1[j];

                    tg_yyyyy_xyyyz_0[j] = pb_y * tg_yyyy_xyyyz_0[j] + fr * tg_yyyy_xyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyz_0[j] - tg_yyy_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xyyz_1[j];

                    tg_yyyyy_xyyzz_0[j] = pb_y * tg_yyyy_xyyzz_0[j] + fr * tg_yyyy_xyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyzz_0[j] - tg_yyy_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xyzz_1[j];

                    tg_yyyyy_xyzzz_0[j] = pb_y * tg_yyyy_xyzzz_0[j] + fr * tg_yyyy_xyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyzzz_0[j] - tg_yyy_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xzzz_1[j];

                    tg_yyyyy_xzzzz_0[j] = pb_y * tg_yyyy_xzzzz_0[j] + fr * tg_yyyy_xzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xzzzz_0[j] - tg_yyy_xzzzz_1[j] * fl1_fza);

                    tg_yyyyy_yyyyy_0[j] = pb_y * tg_yyyy_yyyyy_0[j] + fr * tg_yyyy_yyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyy_0[j] - tg_yyy_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyy_yyyy_1[j];

                    tg_yyyyy_yyyyz_0[j] = pb_y * tg_yyyy_yyyyz_0[j] + fr * tg_yyyy_yyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyz_0[j] - tg_yyy_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_yyyz_1[j];

                    tg_yyyyy_yyyzz_0[j] = pb_y * tg_yyyy_yyyzz_0[j] + fr * tg_yyyy_yyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyzz_0[j] - tg_yyy_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_yyzz_1[j];

                    tg_yyyyy_yyzzz_0[j] = pb_y * tg_yyyy_yyzzz_0[j] + fr * tg_yyyy_yyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyzzz_0[j] - tg_yyy_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_yzzz_1[j];

                    tg_yyyyy_yzzzz_0[j] = pb_y * tg_yyyy_yzzzz_0[j] + fr * tg_yyyy_yzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yzzzz_0[j] - tg_yyy_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_zzzz_1[j];

                    tg_yyyyy_zzzzz_0[j] = pb_y * tg_yyyy_zzzzz_0[j] + fr * tg_yyyy_zzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_zzzzz_0[j] - tg_yyy_zzzzz_1[j] * fl1_fza);

                    tg_yyyyz_xxxxx_0[j] = pb_y * tg_yyyz_xxxxx_0[j] + fr * tg_yyyz_xxxxx_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxx_0[j] - tg_yyz_xxxxx_1[j] * fl1_fza);

                    tg_yyyyz_xxxxy_0[j] = pb_y * tg_yyyz_xxxxy_0[j] + fr * tg_yyyz_xxxxy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxy_0[j] - tg_yyz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxx_1[j];

                    tg_yyyyz_xxxxz_0[j] = pb_y * tg_yyyz_xxxxz_0[j] + fr * tg_yyyz_xxxxz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxz_0[j] - tg_yyz_xxxxz_1[j] * fl1_fza);

                    tg_yyyyz_xxxyy_0[j] = pb_y * tg_yyyz_xxxyy_0[j] + fr * tg_yyyz_xxxyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyy_0[j] - tg_yyz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxxy_1[j];

                    tg_yyyyz_xxxyz_0[j] = pb_y * tg_yyyz_xxxyz_0[j] + fr * tg_yyyz_xxxyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyz_0[j] - tg_yyz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxz_1[j];

                    tg_yyyyz_xxxzz_0[j] = pb_y * tg_yyyz_xxxzz_0[j] + fr * tg_yyyz_xxxzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxzz_0[j] - tg_yyz_xxxzz_1[j] * fl1_fza);

                    tg_yyyyz_xxyyy_0[j] = pb_y * tg_yyyz_xxyyy_0[j] + fr * tg_yyyz_xxyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyy_0[j] - tg_yyz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xxyy_1[j];

                    tg_yyyyz_xxyyz_0[j] = pb_y * tg_yyyz_xxyyz_0[j] + fr * tg_yyyz_xxyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyz_0[j] - tg_yyz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxyz_1[j];

                    tg_yyyyz_xxyzz_0[j] = pb_y * tg_yyyz_xxyzz_0[j] + fr * tg_yyyz_xxyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyzz_0[j] - tg_yyz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxzz_1[j];

                    tg_yyyyz_xxzzz_0[j] = pb_y * tg_yyyz_xxzzz_0[j] + fr * tg_yyyz_xxzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxzzz_0[j] - tg_yyz_xxzzz_1[j] * fl1_fza);

                    tg_yyyyz_xyyyy_0[j] = pb_y * tg_yyyz_xyyyy_0[j] + fr * tg_yyyz_xyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyy_0[j] - tg_yyz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_xyyy_1[j];

                    tg_yyyyz_xyyyz_0[j] = pb_y * tg_yyyz_xyyyz_0[j] + fr * tg_yyyz_xyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyz_0[j] - tg_yyz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xyyz_1[j];

                    tg_yyyyz_xyyzz_0[j] = pb_y * tg_yyyz_xyyzz_0[j] + fr * tg_yyyz_xyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyzz_0[j] - tg_yyz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xyzz_1[j];

                    tg_yyyyz_xyzzz_0[j] = pb_y * tg_yyyz_xyzzz_0[j] + fr * tg_yyyz_xyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyzzz_0[j] - tg_yyz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xzzz_1[j];

                    tg_yyyyz_xzzzz_0[j] = pb_y * tg_yyyz_xzzzz_0[j] + fr * tg_yyyz_xzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xzzzz_0[j] - tg_yyz_xzzzz_1[j] * fl1_fza);

                    tg_yyyyz_yyyyy_0[j] = pb_y * tg_yyyz_yyyyy_0[j] + fr * tg_yyyz_yyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyy_0[j] - tg_yyz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyz_yyyy_1[j];

                    tg_yyyyz_yyyyz_0[j] = pb_y * tg_yyyz_yyyyz_0[j] + fr * tg_yyyz_yyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyz_0[j] - tg_yyz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSH_353_441(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (353,441)

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
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tg_yyyz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 177); 

                auto tg_yyyz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 178); 

                auto tg_yyyz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 179); 

                auto tg_yyzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 180); 

                auto tg_yyzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 181); 

                auto tg_yyzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 182); 

                auto tg_yyzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 183); 

                auto tg_yyzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 184); 

                auto tg_yyzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 185); 

                auto tg_yyzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 186); 

                auto tg_yyzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 187); 

                auto tg_yyzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 188); 

                auto tg_yyzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 189); 

                auto tg_yyzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 190); 

                auto tg_yyzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 191); 

                auto tg_yyzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 192); 

                auto tg_yyzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 193); 

                auto tg_yyzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 194); 

                auto tg_yzzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 195); 

                auto tg_yzzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 196); 

                auto tg_yzzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 197); 

                auto tg_yzzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 198); 

                auto tg_yzzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 199); 

                auto tg_yzzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 200); 

                auto tg_yzzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 201); 

                auto tg_yzzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 202); 

                auto tg_yzzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 203); 

                auto tg_yzzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 204); 

                auto tg_yzzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 205); 

                auto tg_yzzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 206); 

                auto tg_yzzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 207); 

                auto tg_yzzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 208); 

                auto tg_yzzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 209); 

                auto tg_zzzz_xxxx_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 210); 

                auto tg_zzzz_xxxy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 211); 

                auto tg_zzzz_xxxz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 212); 

                auto tg_zzzz_xxyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 213); 

                auto tg_zzzz_xxyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 214); 

                auto tg_zzzz_xxzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 215); 

                auto tg_zzzz_xyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 216); 

                auto tg_zzzz_xyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 217); 

                auto tg_zzzz_xyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 218); 

                auto tg_zzzz_xzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 219); 

                auto tg_zzzz_yyyy_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 220); 

                auto tg_zzzz_yyyz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 221); 

                auto tg_zzzz_yyzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 222); 

                auto tg_zzzz_yzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 223); 

                auto tg_zzzz_zzzz_1 = primBuffer[pidx_g_4_4_m1].data(225 * idx + 224); 

                // set up pointers to integrals

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

                // Batch of Integrals (353,441)

                #pragma omp simd aligned(fxn, fza, tg_yyyyz_yyyzz_0, tg_yyyyz_yyzzz_0, tg_yyyyz_yzzzz_0, \
                                         tg_yyyyz_zzzzz_0, tg_yyyz_yyyzz_0, tg_yyyz_yyyzz_1, tg_yyyz_yyzz_1, tg_yyyz_yyzzz_0, \
                                         tg_yyyz_yyzzz_1, tg_yyyz_yzzz_1, tg_yyyz_yzzzz_0, tg_yyyz_yzzzz_1, tg_yyyz_zzzz_1, \
                                         tg_yyyz_zzzzz_0, tg_yyyz_zzzzz_1, tg_yyyzz_xxxxx_0, tg_yyyzz_xxxxy_0, \
                                         tg_yyyzz_xxxxz_0, tg_yyyzz_xxxyy_0, tg_yyyzz_xxxyz_0, tg_yyyzz_xxxzz_0, \
                                         tg_yyyzz_xxyyy_0, tg_yyyzz_xxyyz_0, tg_yyyzz_xxyzz_0, tg_yyyzz_xxzzz_0, \
                                         tg_yyyzz_xyyyy_0, tg_yyyzz_xyyyz_0, tg_yyyzz_xyyzz_0, tg_yyyzz_xyzzz_0, \
                                         tg_yyyzz_xzzzz_0, tg_yyyzz_yyyyy_0, tg_yyyzz_yyyyz_0, tg_yyyzz_yyyzz_0, \
                                         tg_yyyzz_yyzzz_0, tg_yyyzz_yzzzz_0, tg_yyyzz_zzzzz_0, tg_yyz_yyyzz_0, tg_yyz_yyyzz_1, \
                                         tg_yyz_yyzzz_0, tg_yyz_yyzzz_1, tg_yyz_yzzzz_0, tg_yyz_yzzzz_1, tg_yyz_zzzzz_0, \
                                         tg_yyz_zzzzz_1, tg_yyzz_xxxx_1, tg_yyzz_xxxxx_0, tg_yyzz_xxxxx_1, tg_yyzz_xxxxy_0, \
                                         tg_yyzz_xxxxy_1, tg_yyzz_xxxxz_0, tg_yyzz_xxxxz_1, tg_yyzz_xxxy_1, tg_yyzz_xxxyy_0, \
                                         tg_yyzz_xxxyy_1, tg_yyzz_xxxyz_0, tg_yyzz_xxxyz_1, tg_yyzz_xxxz_1, tg_yyzz_xxxzz_0, \
                                         tg_yyzz_xxxzz_1, tg_yyzz_xxyy_1, tg_yyzz_xxyyy_0, tg_yyzz_xxyyy_1, tg_yyzz_xxyyz_0, \
                                         tg_yyzz_xxyyz_1, tg_yyzz_xxyz_1, tg_yyzz_xxyzz_0, tg_yyzz_xxyzz_1, tg_yyzz_xxzz_1, \
                                         tg_yyzz_xxzzz_0, tg_yyzz_xxzzz_1, tg_yyzz_xyyy_1, tg_yyzz_xyyyy_0, tg_yyzz_xyyyy_1, \
                                         tg_yyzz_xyyyz_0, tg_yyzz_xyyyz_1, tg_yyzz_xyyz_1, tg_yyzz_xyyzz_0, tg_yyzz_xyyzz_1, \
                                         tg_yyzz_xyzz_1, tg_yyzz_xyzzz_0, tg_yyzz_xyzzz_1, tg_yyzz_xzzz_1, tg_yyzz_xzzzz_0, \
                                         tg_yyzz_xzzzz_1, tg_yyzz_yyyy_1, tg_yyzz_yyyyy_0, tg_yyzz_yyyyy_1, tg_yyzz_yyyyz_0, \
                                         tg_yyzz_yyyyz_1, tg_yyzz_yyyz_1, tg_yyzz_yyyzz_0, tg_yyzz_yyyzz_1, tg_yyzz_yyzz_1, \
                                         tg_yyzz_yyzzz_0, tg_yyzz_yyzzz_1, tg_yyzz_yzzz_1, tg_yyzz_yzzzz_0, tg_yyzz_yzzzz_1, \
                                         tg_yyzz_zzzz_1, tg_yyzz_zzzzz_0, tg_yyzz_zzzzz_1, tg_yyzzz_xxxxx_0, \
                                         tg_yyzzz_xxxxy_0, tg_yyzzz_xxxxz_0, tg_yyzzz_xxxyy_0, tg_yyzzz_xxxyz_0, \
                                         tg_yyzzz_xxxzz_0, tg_yyzzz_xxyyy_0, tg_yyzzz_xxyyz_0, tg_yyzzz_xxyzz_0, \
                                         tg_yyzzz_xxzzz_0, tg_yyzzz_xyyyy_0, tg_yyzzz_xyyyz_0, tg_yyzzz_xyyzz_0, \
                                         tg_yyzzz_xyzzz_0, tg_yyzzz_xzzzz_0, tg_yyzzz_yyyyy_0, tg_yyzzz_yyyyz_0, \
                                         tg_yyzzz_yyyzz_0, tg_yyzzz_yyzzz_0, tg_yyzzz_yzzzz_0, tg_yyzzz_zzzzz_0, \
                                         tg_yzz_xxxxx_0, tg_yzz_xxxxx_1, tg_yzz_xxxxy_0, tg_yzz_xxxxy_1, tg_yzz_xxxxz_0, \
                                         tg_yzz_xxxxz_1, tg_yzz_xxxyy_0, tg_yzz_xxxyy_1, tg_yzz_xxxyz_0, tg_yzz_xxxyz_1, \
                                         tg_yzz_xxxzz_0, tg_yzz_xxxzz_1, tg_yzz_xxyyy_0, tg_yzz_xxyyy_1, tg_yzz_xxyyz_0, \
                                         tg_yzz_xxyyz_1, tg_yzz_xxyzz_0, tg_yzz_xxyzz_1, tg_yzz_xxzzz_0, tg_yzz_xxzzz_1, \
                                         tg_yzz_xyyyy_0, tg_yzz_xyyyy_1, tg_yzz_xyyyz_0, tg_yzz_xyyyz_1, tg_yzz_xyyzz_0, \
                                         tg_yzz_xyyzz_1, tg_yzz_xyzzz_0, tg_yzz_xyzzz_1, tg_yzz_xzzzz_0, tg_yzz_xzzzz_1, \
                                         tg_yzz_yyyyy_0, tg_yzz_yyyyy_1, tg_yzz_yyyyz_0, tg_yzz_yyyyz_1, tg_yzz_yyyzz_0, \
                                         tg_yzz_yyyzz_1, tg_yzz_yyzzz_0, tg_yzz_yyzzz_1, tg_yzz_yzzzz_0, tg_yzz_yzzzz_1, \
                                         tg_yzz_zzzzz_0, tg_yzz_zzzzz_1, tg_yzzz_xxxx_1, tg_yzzz_xxxxx_0, tg_yzzz_xxxxx_1, \
                                         tg_yzzz_xxxxy_0, tg_yzzz_xxxxy_1, tg_yzzz_xxxxz_0, tg_yzzz_xxxxz_1, tg_yzzz_xxxy_1, \
                                         tg_yzzz_xxxyy_0, tg_yzzz_xxxyy_1, tg_yzzz_xxxyz_0, tg_yzzz_xxxyz_1, tg_yzzz_xxxz_1, \
                                         tg_yzzz_xxxzz_0, tg_yzzz_xxxzz_1, tg_yzzz_xxyy_1, tg_yzzz_xxyyy_0, tg_yzzz_xxyyy_1, \
                                         tg_yzzz_xxyyz_0, tg_yzzz_xxyyz_1, tg_yzzz_xxyz_1, tg_yzzz_xxyzz_0, tg_yzzz_xxyzz_1, \
                                         tg_yzzz_xxzz_1, tg_yzzz_xxzzz_0, tg_yzzz_xxzzz_1, tg_yzzz_xyyy_1, tg_yzzz_xyyyy_0, \
                                         tg_yzzz_xyyyy_1, tg_yzzz_xyyyz_0, tg_yzzz_xyyyz_1, tg_yzzz_xyyz_1, tg_yzzz_xyyzz_0, \
                                         tg_yzzz_xyyzz_1, tg_yzzz_xyzz_1, tg_yzzz_xyzzz_0, tg_yzzz_xyzzz_1, tg_yzzz_xzzz_1, \
                                         tg_yzzz_xzzzz_0, tg_yzzz_xzzzz_1, tg_yzzz_yyyy_1, tg_yzzz_yyyyy_0, tg_yzzz_yyyyy_1, \
                                         tg_yzzz_yyyyz_0, tg_yzzz_yyyyz_1, tg_yzzz_yyyz_1, tg_yzzz_yyyzz_0, tg_yzzz_yyyzz_1, \
                                         tg_yzzz_yyzz_1, tg_yzzz_yyzzz_0, tg_yzzz_yyzzz_1, tg_yzzz_yzzz_1, tg_yzzz_yzzzz_0, \
                                         tg_yzzz_yzzzz_1, tg_yzzz_zzzz_1, tg_yzzz_zzzzz_0, tg_yzzz_zzzzz_1, tg_yzzzz_xxxxx_0, \
                                         tg_yzzzz_xxxxy_0, tg_yzzzz_xxxxz_0, tg_yzzzz_xxxyy_0, tg_yzzzz_xxxyz_0, \
                                         tg_yzzzz_xxxzz_0, tg_yzzzz_xxyyy_0, tg_yzzzz_xxyyz_0, tg_yzzzz_xxyzz_0, \
                                         tg_yzzzz_xxzzz_0, tg_yzzzz_xyyyy_0, tg_yzzzz_xyyyz_0, tg_yzzzz_xyyzz_0, \
                                         tg_yzzzz_xyzzz_0, tg_yzzzz_xzzzz_0, tg_yzzzz_yyyyy_0, tg_yzzzz_yyyyz_0, \
                                         tg_yzzzz_yyyzz_0, tg_yzzzz_yyzzz_0, tg_yzzzz_yzzzz_0, tg_yzzzz_zzzzz_0, \
                                         tg_zzz_xxxxx_0, tg_zzz_xxxxx_1, tg_zzz_xxxxy_0, tg_zzz_xxxxy_1, tg_zzz_xxxxz_0, \
                                         tg_zzz_xxxxz_1, tg_zzz_xxxyy_0, tg_zzz_xxxyy_1, tg_zzz_xxxyz_0, tg_zzz_xxxyz_1, \
                                         tg_zzz_xxxzz_0, tg_zzz_xxxzz_1, tg_zzz_xxyyy_0, tg_zzz_xxyyy_1, tg_zzz_xxyyz_0, \
                                         tg_zzz_xxyyz_1, tg_zzz_xxyzz_0, tg_zzz_xxyzz_1, tg_zzz_xxzzz_0, tg_zzz_xxzzz_1, \
                                         tg_zzz_xyyyy_0, tg_zzz_xyyyy_1, tg_zzz_xyyyz_0, tg_zzz_xyyyz_1, tg_zzz_xyyzz_0, \
                                         tg_zzz_xyyzz_1, tg_zzz_xyzzz_0, tg_zzz_xyzzz_1, tg_zzz_xzzzz_0, tg_zzz_xzzzz_1, \
                                         tg_zzz_yyyyy_0, tg_zzz_yyyyy_1, tg_zzz_yyyyz_0, tg_zzz_yyyyz_1, tg_zzz_yyyzz_0, \
                                         tg_zzz_yyyzz_1, tg_zzz_yyzzz_0, tg_zzz_yyzzz_1, tg_zzz_yzzzz_0, tg_zzz_yzzzz_1, \
                                         tg_zzz_zzzzz_0, tg_zzz_zzzzz_1, tg_zzzz_xxxx_1, tg_zzzz_xxxxx_0, tg_zzzz_xxxxx_1, \
                                         tg_zzzz_xxxxy_0, tg_zzzz_xxxxy_1, tg_zzzz_xxxxz_0, tg_zzzz_xxxxz_1, tg_zzzz_xxxy_1, \
                                         tg_zzzz_xxxyy_0, tg_zzzz_xxxyy_1, tg_zzzz_xxxyz_0, tg_zzzz_xxxyz_1, tg_zzzz_xxxz_1, \
                                         tg_zzzz_xxxzz_0, tg_zzzz_xxxzz_1, tg_zzzz_xxyy_1, tg_zzzz_xxyyy_0, tg_zzzz_xxyyy_1, \
                                         tg_zzzz_xxyyz_0, tg_zzzz_xxyyz_1, tg_zzzz_xxyz_1, tg_zzzz_xxyzz_0, tg_zzzz_xxyzz_1, \
                                         tg_zzzz_xxzz_1, tg_zzzz_xxzzz_0, tg_zzzz_xxzzz_1, tg_zzzz_xyyy_1, tg_zzzz_xyyyy_0, \
                                         tg_zzzz_xyyyy_1, tg_zzzz_xyyyz_0, tg_zzzz_xyyyz_1, tg_zzzz_xyyz_1, tg_zzzz_xyyzz_0, \
                                         tg_zzzz_xyyzz_1, tg_zzzz_xyzz_1, tg_zzzz_xyzzz_0, tg_zzzz_xyzzz_1, tg_zzzz_xzzz_1, \
                                         tg_zzzz_xzzzz_0, tg_zzzz_xzzzz_1, tg_zzzz_yyyy_1, tg_zzzz_yyyyy_0, tg_zzzz_yyyyy_1, \
                                         tg_zzzz_yyyyz_0, tg_zzzz_yyyyz_1, tg_zzzz_yyyz_1, tg_zzzz_yyyzz_0, tg_zzzz_yyyzz_1, \
                                         tg_zzzz_yyzz_1, tg_zzzz_yyzzz_0, tg_zzzz_yyzzz_1, tg_zzzz_yzzz_1, tg_zzzz_yzzzz_0, \
                                         tg_zzzz_yzzzz_1, tg_zzzz_zzzz_1, tg_zzzz_zzzzz_0, tg_zzzz_zzzzz_1, tg_zzzzz_xxxxx_0, \
                                         tg_zzzzz_xxxxy_0, tg_zzzzz_xxxxz_0, tg_zzzzz_xxxyy_0, tg_zzzzz_xxxyz_0, \
                                         tg_zzzzz_xxxzz_0, tg_zzzzz_xxyyy_0, tg_zzzzz_xxyyz_0, tg_zzzzz_xxyzz_0, \
                                         tg_zzzzz_xxzzz_0, tg_zzzzz_xyyyy_0, tg_zzzzz_xyyyz_0, tg_zzzzz_xyyzz_0, \
                                         tg_zzzzz_xyzzz_0, tg_zzzzz_xzzzz_0, tg_zzzzz_yyyyy_0, tg_zzzzz_yyyyz_0, \
                                         tg_zzzzz_yyyzz_0, tg_zzzzz_yyzzz_0, tg_zzzzz_yzzzz_0, tg_zzzzz_zzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyz_yyyzz_0[j] = pb_y * tg_yyyz_yyyzz_0[j] + fr * tg_yyyz_yyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyzz_0[j] - tg_yyz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_yyzz_1[j];

                    tg_yyyyz_yyzzz_0[j] = pb_y * tg_yyyz_yyzzz_0[j] + fr * tg_yyyz_yyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyzzz_0[j] - tg_yyz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_yzzz_1[j];

                    tg_yyyyz_yzzzz_0[j] = pb_y * tg_yyyz_yzzzz_0[j] + fr * tg_yyyz_yzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yzzzz_0[j] - tg_yyz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_zzzz_1[j];

                    tg_yyyyz_zzzzz_0[j] = pb_y * tg_yyyz_zzzzz_0[j] + fr * tg_yyyz_zzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_zzzzz_0[j] - tg_yyz_zzzzz_1[j] * fl1_fza);

                    tg_yyyzz_xxxxx_0[j] = pb_y * tg_yyzz_xxxxx_0[j] + fr * tg_yyzz_xxxxx_1[j] + fl1_fx * (tg_yzz_xxxxx_0[j] - tg_yzz_xxxxx_1[j] * fl1_fza);

                    tg_yyyzz_xxxxy_0[j] = pb_y * tg_yyzz_xxxxy_0[j] + fr * tg_yyzz_xxxxy_1[j] + fl1_fx * (tg_yzz_xxxxy_0[j] - tg_yzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxx_1[j];

                    tg_yyyzz_xxxxz_0[j] = pb_y * tg_yyzz_xxxxz_0[j] + fr * tg_yyzz_xxxxz_1[j] + fl1_fx * (tg_yzz_xxxxz_0[j] - tg_yzz_xxxxz_1[j] * fl1_fza);

                    tg_yyyzz_xxxyy_0[j] = pb_y * tg_yyzz_xxxyy_0[j] + fr * tg_yyzz_xxxyy_1[j] + fl1_fx * (tg_yzz_xxxyy_0[j] - tg_yzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxxy_1[j];

                    tg_yyyzz_xxxyz_0[j] = pb_y * tg_yyzz_xxxyz_0[j] + fr * tg_yyzz_xxxyz_1[j] + fl1_fx * (tg_yzz_xxxyz_0[j] - tg_yzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxz_1[j];

                    tg_yyyzz_xxxzz_0[j] = pb_y * tg_yyzz_xxxzz_0[j] + fr * tg_yyzz_xxxzz_1[j] + fl1_fx * (tg_yzz_xxxzz_0[j] - tg_yzz_xxxzz_1[j] * fl1_fza);

                    tg_yyyzz_xxyyy_0[j] = pb_y * tg_yyzz_xxyyy_0[j] + fr * tg_yyzz_xxyyy_1[j] + fl1_fx * (tg_yzz_xxyyy_0[j] - tg_yzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xxyy_1[j];

                    tg_yyyzz_xxyyz_0[j] = pb_y * tg_yyzz_xxyyz_0[j] + fr * tg_yyzz_xxyyz_1[j] + fl1_fx * (tg_yzz_xxyyz_0[j] - tg_yzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxyz_1[j];

                    tg_yyyzz_xxyzz_0[j] = pb_y * tg_yyzz_xxyzz_0[j] + fr * tg_yyzz_xxyzz_1[j] + fl1_fx * (tg_yzz_xxyzz_0[j] - tg_yzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxzz_1[j];

                    tg_yyyzz_xxzzz_0[j] = pb_y * tg_yyzz_xxzzz_0[j] + fr * tg_yyzz_xxzzz_1[j] + fl1_fx * (tg_yzz_xxzzz_0[j] - tg_yzz_xxzzz_1[j] * fl1_fza);

                    tg_yyyzz_xyyyy_0[j] = pb_y * tg_yyzz_xyyyy_0[j] + fr * tg_yyzz_xyyyy_1[j] + fl1_fx * (tg_yzz_xyyyy_0[j] - tg_yzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_xyyy_1[j];

                    tg_yyyzz_xyyyz_0[j] = pb_y * tg_yyzz_xyyyz_0[j] + fr * tg_yyzz_xyyyz_1[j] + fl1_fx * (tg_yzz_xyyyz_0[j] - tg_yzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xyyz_1[j];

                    tg_yyyzz_xyyzz_0[j] = pb_y * tg_yyzz_xyyzz_0[j] + fr * tg_yyzz_xyyzz_1[j] + fl1_fx * (tg_yzz_xyyzz_0[j] - tg_yzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xyzz_1[j];

                    tg_yyyzz_xyzzz_0[j] = pb_y * tg_yyzz_xyzzz_0[j] + fr * tg_yyzz_xyzzz_1[j] + fl1_fx * (tg_yzz_xyzzz_0[j] - tg_yzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xzzz_1[j];

                    tg_yyyzz_xzzzz_0[j] = pb_y * tg_yyzz_xzzzz_0[j] + fr * tg_yyzz_xzzzz_1[j] + fl1_fx * (tg_yzz_xzzzz_0[j] - tg_yzz_xzzzz_1[j] * fl1_fza);

                    tg_yyyzz_yyyyy_0[j] = pb_y * tg_yyzz_yyyyy_0[j] + fr * tg_yyzz_yyyyy_1[j] + fl1_fx * (tg_yzz_yyyyy_0[j] - tg_yzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzz_yyyy_1[j];

                    tg_yyyzz_yyyyz_0[j] = pb_y * tg_yyzz_yyyyz_0[j] + fr * tg_yyzz_yyyyz_1[j] + fl1_fx * (tg_yzz_yyyyz_0[j] - tg_yzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_yyyz_1[j];

                    tg_yyyzz_yyyzz_0[j] = pb_y * tg_yyzz_yyyzz_0[j] + fr * tg_yyzz_yyyzz_1[j] + fl1_fx * (tg_yzz_yyyzz_0[j] - tg_yzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_yyzz_1[j];

                    tg_yyyzz_yyzzz_0[j] = pb_y * tg_yyzz_yyzzz_0[j] + fr * tg_yyzz_yyzzz_1[j] + fl1_fx * (tg_yzz_yyzzz_0[j] - tg_yzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_yzzz_1[j];

                    tg_yyyzz_yzzzz_0[j] = pb_y * tg_yyzz_yzzzz_0[j] + fr * tg_yyzz_yzzzz_1[j] + fl1_fx * (tg_yzz_yzzzz_0[j] - tg_yzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_zzzz_1[j];

                    tg_yyyzz_zzzzz_0[j] = pb_y * tg_yyzz_zzzzz_0[j] + fr * tg_yyzz_zzzzz_1[j] + fl1_fx * (tg_yzz_zzzzz_0[j] - tg_yzz_zzzzz_1[j] * fl1_fza);

                    tg_yyzzz_xxxxx_0[j] = pb_y * tg_yzzz_xxxxx_0[j] + fr * tg_yzzz_xxxxx_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxx_0[j] - tg_zzz_xxxxx_1[j] * fl1_fza);

                    tg_yyzzz_xxxxy_0[j] = pb_y * tg_yzzz_xxxxy_0[j] + fr * tg_yzzz_xxxxy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxy_0[j] - tg_zzz_xxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxx_1[j];

                    tg_yyzzz_xxxxz_0[j] = pb_y * tg_yzzz_xxxxz_0[j] + fr * tg_yzzz_xxxxz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxz_0[j] - tg_zzz_xxxxz_1[j] * fl1_fza);

                    tg_yyzzz_xxxyy_0[j] = pb_y * tg_yzzz_xxxyy_0[j] + fr * tg_yzzz_xxxyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyy_0[j] - tg_zzz_xxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxxy_1[j];

                    tg_yyzzz_xxxyz_0[j] = pb_y * tg_yzzz_xxxyz_0[j] + fr * tg_yzzz_xxxyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyz_0[j] - tg_zzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxz_1[j];

                    tg_yyzzz_xxxzz_0[j] = pb_y * tg_yzzz_xxxzz_0[j] + fr * tg_yzzz_xxxzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxzz_0[j] - tg_zzz_xxxzz_1[j] * fl1_fza);

                    tg_yyzzz_xxyyy_0[j] = pb_y * tg_yzzz_xxyyy_0[j] + fr * tg_yzzz_xxyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyy_0[j] - tg_zzz_xxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xxyy_1[j];

                    tg_yyzzz_xxyyz_0[j] = pb_y * tg_yzzz_xxyyz_0[j] + fr * tg_yzzz_xxyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyz_0[j] - tg_zzz_xxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxyz_1[j];

                    tg_yyzzz_xxyzz_0[j] = pb_y * tg_yzzz_xxyzz_0[j] + fr * tg_yzzz_xxyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyzz_0[j] - tg_zzz_xxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxzz_1[j];

                    tg_yyzzz_xxzzz_0[j] = pb_y * tg_yzzz_xxzzz_0[j] + fr * tg_yzzz_xxzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxzzz_0[j] - tg_zzz_xxzzz_1[j] * fl1_fza);

                    tg_yyzzz_xyyyy_0[j] = pb_y * tg_yzzz_xyyyy_0[j] + fr * tg_yzzz_xyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyy_0[j] - tg_zzz_xyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_xyyy_1[j];

                    tg_yyzzz_xyyyz_0[j] = pb_y * tg_yzzz_xyyyz_0[j] + fr * tg_yzzz_xyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyz_0[j] - tg_zzz_xyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xyyz_1[j];

                    tg_yyzzz_xyyzz_0[j] = pb_y * tg_yzzz_xyyzz_0[j] + fr * tg_yzzz_xyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyzz_0[j] - tg_zzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xyzz_1[j];

                    tg_yyzzz_xyzzz_0[j] = pb_y * tg_yzzz_xyzzz_0[j] + fr * tg_yzzz_xyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyzzz_0[j] - tg_zzz_xyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xzzz_1[j];

                    tg_yyzzz_xzzzz_0[j] = pb_y * tg_yzzz_xzzzz_0[j] + fr * tg_yzzz_xzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xzzzz_0[j] - tg_zzz_xzzzz_1[j] * fl1_fza);

                    tg_yyzzz_yyyyy_0[j] = pb_y * tg_yzzz_yyyyy_0[j] + fr * tg_yzzz_yyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyy_0[j] - tg_zzz_yyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzz_yyyy_1[j];

                    tg_yyzzz_yyyyz_0[j] = pb_y * tg_yzzz_yyyyz_0[j] + fr * tg_yzzz_yyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyz_0[j] - tg_zzz_yyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_yyyz_1[j];

                    tg_yyzzz_yyyzz_0[j] = pb_y * tg_yzzz_yyyzz_0[j] + fr * tg_yzzz_yyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyzz_0[j] - tg_zzz_yyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_yyzz_1[j];

                    tg_yyzzz_yyzzz_0[j] = pb_y * tg_yzzz_yyzzz_0[j] + fr * tg_yzzz_yyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyzzz_0[j] - tg_zzz_yyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_yzzz_1[j];

                    tg_yyzzz_yzzzz_0[j] = pb_y * tg_yzzz_yzzzz_0[j] + fr * tg_yzzz_yzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yzzzz_0[j] - tg_zzz_yzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_zzzz_1[j];

                    tg_yyzzz_zzzzz_0[j] = pb_y * tg_yzzz_zzzzz_0[j] + fr * tg_yzzz_zzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_zzzzz_0[j] - tg_zzz_zzzzz_1[j] * fl1_fza);

                    tg_yzzzz_xxxxx_0[j] = pb_y * tg_zzzz_xxxxx_0[j] + fr * tg_zzzz_xxxxx_1[j];

                    tg_yzzzz_xxxxy_0[j] = pb_y * tg_zzzz_xxxxy_0[j] + fr * tg_zzzz_xxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxx_1[j];

                    tg_yzzzz_xxxxz_0[j] = pb_y * tg_zzzz_xxxxz_0[j] + fr * tg_zzzz_xxxxz_1[j];

                    tg_yzzzz_xxxyy_0[j] = pb_y * tg_zzzz_xxxyy_0[j] + fr * tg_zzzz_xxxyy_1[j] + fl1_fxn * tg_zzzz_xxxy_1[j];

                    tg_yzzzz_xxxyz_0[j] = pb_y * tg_zzzz_xxxyz_0[j] + fr * tg_zzzz_xxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxz_1[j];

                    tg_yzzzz_xxxzz_0[j] = pb_y * tg_zzzz_xxxzz_0[j] + fr * tg_zzzz_xxxzz_1[j];

                    tg_yzzzz_xxyyy_0[j] = pb_y * tg_zzzz_xxyyy_0[j] + fr * tg_zzzz_xxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyy_1[j];

                    tg_yzzzz_xxyyz_0[j] = pb_y * tg_zzzz_xxyyz_0[j] + fr * tg_zzzz_xxyyz_1[j] + fl1_fxn * tg_zzzz_xxyz_1[j];

                    tg_yzzzz_xxyzz_0[j] = pb_y * tg_zzzz_xxyzz_0[j] + fr * tg_zzzz_xxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxzz_1[j];

                    tg_yzzzz_xxzzz_0[j] = pb_y * tg_zzzz_xxzzz_0[j] + fr * tg_zzzz_xxzzz_1[j];

                    tg_yzzzz_xyyyy_0[j] = pb_y * tg_zzzz_xyyyy_0[j] + fr * tg_zzzz_xyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xyyy_1[j];

                    tg_yzzzz_xyyyz_0[j] = pb_y * tg_zzzz_xyyyz_0[j] + fr * tg_zzzz_xyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xyyz_1[j];

                    tg_yzzzz_xyyzz_0[j] = pb_y * tg_zzzz_xyyzz_0[j] + fr * tg_zzzz_xyyzz_1[j] + fl1_fxn * tg_zzzz_xyzz_1[j];

                    tg_yzzzz_xyzzz_0[j] = pb_y * tg_zzzz_xyzzz_0[j] + fr * tg_zzzz_xyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xzzz_1[j];

                    tg_yzzzz_xzzzz_0[j] = pb_y * tg_zzzz_xzzzz_0[j] + fr * tg_zzzz_xzzzz_1[j];

                    tg_yzzzz_yyyyy_0[j] = pb_y * tg_zzzz_yyyyy_0[j] + fr * tg_zzzz_yyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzz_yyyy_1[j];

                    tg_yzzzz_yyyyz_0[j] = pb_y * tg_zzzz_yyyyz_0[j] + fr * tg_zzzz_yyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzz_yyyz_1[j];

                    tg_yzzzz_yyyzz_0[j] = pb_y * tg_zzzz_yyyzz_0[j] + fr * tg_zzzz_yyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_yyzz_1[j];

                    tg_yzzzz_yyzzz_0[j] = pb_y * tg_zzzz_yyzzz_0[j] + fr * tg_zzzz_yyzzz_1[j] + fl1_fxn * tg_zzzz_yzzz_1[j];

                    tg_yzzzz_yzzzz_0[j] = pb_y * tg_zzzz_yzzzz_0[j] + fr * tg_zzzz_yzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzz_1[j];

                    tg_yzzzz_zzzzz_0[j] = pb_y * tg_zzzz_zzzzz_0[j] + fr * tg_zzzz_zzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzz_xxxxx_0[j] = pb_z * tg_zzzz_xxxxx_0[j] + fr * tg_zzzz_xxxxx_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxx_0[j] - tg_zzz_xxxxx_1[j] * fl1_fza);

                    tg_zzzzz_xxxxy_0[j] = pb_z * tg_zzzz_xxxxy_0[j] + fr * tg_zzzz_xxxxy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxy_0[j] - tg_zzz_xxxxy_1[j] * fl1_fza);

                    tg_zzzzz_xxxxz_0[j] = pb_z * tg_zzzz_xxxxz_0[j] + fr * tg_zzzz_xxxxz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxz_0[j] - tg_zzz_xxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxx_1[j];

                    tg_zzzzz_xxxyy_0[j] = pb_z * tg_zzzz_xxxyy_0[j] + fr * tg_zzzz_xxxyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyy_0[j] - tg_zzz_xxxyy_1[j] * fl1_fza);

                    tg_zzzzz_xxxyz_0[j] = pb_z * tg_zzzz_xxxyz_0[j] + fr * tg_zzzz_xxxyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyz_0[j] - tg_zzz_xxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxy_1[j];

                    tg_zzzzz_xxxzz_0[j] = pb_z * tg_zzzz_xxxzz_0[j] + fr * tg_zzzz_xxxzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxzz_0[j] - tg_zzz_xxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxxz_1[j];

                    tg_zzzzz_xxyyy_0[j] = pb_z * tg_zzzz_xxyyy_0[j] + fr * tg_zzzz_xxyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyy_0[j] - tg_zzz_xxyyy_1[j] * fl1_fza);

                    tg_zzzzz_xxyyz_0[j] = pb_z * tg_zzzz_xxyyz_0[j] + fr * tg_zzzz_xxyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyz_0[j] - tg_zzz_xxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxyy_1[j];

                    tg_zzzzz_xxyzz_0[j] = pb_z * tg_zzzz_xxyzz_0[j] + fr * tg_zzzz_xxyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyzz_0[j] - tg_zzz_xxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxyz_1[j];

                    tg_zzzzz_xxzzz_0[j] = pb_z * tg_zzzz_xxzzz_0[j] + fr * tg_zzzz_xxzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxzzz_0[j] - tg_zzz_xxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xxzz_1[j];

                    tg_zzzzz_xyyyy_0[j] = pb_z * tg_zzzz_xyyyy_0[j] + fr * tg_zzzz_xyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyy_0[j] - tg_zzz_xyyyy_1[j] * fl1_fza);

                    tg_zzzzz_xyyyz_0[j] = pb_z * tg_zzzz_xyyyz_0[j] + fr * tg_zzzz_xyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyz_0[j] - tg_zzz_xyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xyyy_1[j];

                    tg_zzzzz_xyyzz_0[j] = pb_z * tg_zzzz_xyyzz_0[j] + fr * tg_zzzz_xyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyzz_0[j] - tg_zzz_xyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xyyz_1[j];

                    tg_zzzzz_xyzzz_0[j] = pb_z * tg_zzzz_xyzzz_0[j] + fr * tg_zzzz_xyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyzzz_0[j] - tg_zzz_xyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xyzz_1[j];

                    tg_zzzzz_xzzzz_0[j] = pb_z * tg_zzzz_xzzzz_0[j] + fr * tg_zzzz_xzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xzzzz_0[j] - tg_zzz_xzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_xzzz_1[j];

                    tg_zzzzz_yyyyy_0[j] = pb_z * tg_zzzz_yyyyy_0[j] + fr * tg_zzzz_yyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyy_0[j] - tg_zzz_yyyyy_1[j] * fl1_fza);

                    tg_zzzzz_yyyyz_0[j] = pb_z * tg_zzzz_yyyyz_0[j] + fr * tg_zzzz_yyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyz_0[j] - tg_zzz_yyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_yyyy_1[j];

                    tg_zzzzz_yyyzz_0[j] = pb_z * tg_zzzz_yyyzz_0[j] + fr * tg_zzzz_yyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyzz_0[j] - tg_zzz_yyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_yyyz_1[j];

                    tg_zzzzz_yyzzz_0[j] = pb_z * tg_zzzz_yyzzz_0[j] + fr * tg_zzzz_yyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyzzz_0[j] - tg_zzz_yyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_yyzz_1[j];

                    tg_zzzzz_yzzzz_0[j] = pb_z * tg_zzzz_yzzzz_0[j] + fr * tg_zzzz_yzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yzzzz_0[j] - tg_zzz_yzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_yzzz_1[j];

                    tg_zzzzz_zzzzz_0[j] = pb_z * tg_zzzz_zzzzz_0[j] + fr * tg_zzzz_zzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_zzzzz_0[j] - tg_zzz_zzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

