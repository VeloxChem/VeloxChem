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

#include "ElectronRepulsionRecFuncForKI.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSKSI(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSKSI_0_92(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_92_184(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_184_276(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_276_368(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_368_460(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_460_552(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_552_644(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_644_735(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_735_826(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_826_917(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSI_917_1008(primBuffer,
                                                          recursionMap,
                                                          osFactors,
                                                          wpDistances, 
                                                          braGtoPairsBlock,
                                                          ketGtoPairsBlock,
                                                          nKetPrimPairs,
                                                          iContrPair); 
    }

    void
    compElectronRepulsionForSKSI_0_92(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,92)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxxxxx_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx); 

                auto tg_xxxxxx_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 1); 

                auto tg_xxxxxx_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 2); 

                auto tg_xxxxxx_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 3); 

                auto tg_xxxxxx_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 4); 

                auto tg_xxxxxx_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 5); 

                auto tg_xxxxxx_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 6); 

                auto tg_xxxxxx_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 7); 

                auto tg_xxxxxx_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 8); 

                auto tg_xxxxxx_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 9); 

                auto tg_xxxxxx_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 10); 

                auto tg_xxxxxx_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 11); 

                auto tg_xxxxxx_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 12); 

                auto tg_xxxxxx_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 13); 

                auto tg_xxxxxx_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 14); 

                auto tg_xxxxxx_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 15); 

                auto tg_xxxxxx_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 16); 

                auto tg_xxxxxx_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 17); 

                auto tg_xxxxxx_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 18); 

                auto tg_xxxxxx_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 19); 

                auto tg_xxxxxx_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 20); 

                auto tg_xxxxxx_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 21); 

                auto tg_xxxxxx_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 22); 

                auto tg_xxxxxx_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 23); 

                auto tg_xxxxxx_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 24); 

                auto tg_xxxxxx_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 25); 

                auto tg_xxxxxx_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 26); 

                auto tg_xxxxxx_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 27); 

                auto tg_xxxxxy_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 28); 

                auto tg_xxxxxy_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 29); 

                auto tg_xxxxxy_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 30); 

                auto tg_xxxxxy_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 31); 

                auto tg_xxxxxy_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 32); 

                auto tg_xxxxxy_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 33); 

                auto tg_xxxxxy_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 34); 

                auto tg_xxxxxy_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 35); 

                auto tg_xxxxxy_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 36); 

                auto tg_xxxxxy_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 37); 

                auto tg_xxxxxy_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 38); 

                auto tg_xxxxxy_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 39); 

                auto tg_xxxxxy_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 40); 

                auto tg_xxxxxy_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 41); 

                auto tg_xxxxxy_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 42); 

                auto tg_xxxxxy_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 43); 

                auto tg_xxxxxy_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 44); 

                auto tg_xxxxxy_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 45); 

                auto tg_xxxxxy_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 46); 

                auto tg_xxxxxy_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 47); 

                auto tg_xxxxxy_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 48); 

                auto tg_xxxxxy_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 49); 

                auto tg_xxxxxy_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 50); 

                auto tg_xxxxxy_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 51); 

                auto tg_xxxxxy_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 52); 

                auto tg_xxxxxy_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 53); 

                auto tg_xxxxxy_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 54); 

                auto tg_xxxxxy_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 55); 

                auto tg_xxxxxz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 56); 

                auto tg_xxxxxz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 57); 

                auto tg_xxxxxz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 58); 

                auto tg_xxxxxz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 59); 

                auto tg_xxxxxz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 60); 

                auto tg_xxxxxz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 61); 

                auto tg_xxxxxz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 62); 

                auto tg_xxxxxz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 63); 

                auto tg_xxxxxz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 64); 

                auto tg_xxxxxz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 65); 

                auto tg_xxxxxz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 66); 

                auto tg_xxxxxz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 67); 

                auto tg_xxxxxz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 68); 

                auto tg_xxxxxz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 69); 

                auto tg_xxxxxz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 70); 

                auto tg_xxxxxz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 71); 

                auto tg_xxxxxz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 72); 

                auto tg_xxxxxz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 73); 

                auto tg_xxxxxz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 74); 

                auto tg_xxxxxz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 75); 

                auto tg_xxxxxz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 76); 

                auto tg_xxxxxz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 77); 

                auto tg_xxxxxz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 78); 

                auto tg_xxxxxz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 79); 

                auto tg_xxxxxz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 80); 

                auto tg_xxxxxz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 81); 

                auto tg_xxxxxz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 82); 

                auto tg_xxxxxz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 83); 

                auto tg_xxxxyy_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 84); 

                auto tg_xxxxyy_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 85); 

                auto tg_xxxxyy_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 86); 

                auto tg_xxxxyy_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 87); 

                auto tg_xxxxyy_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 88); 

                auto tg_xxxxyy_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 89); 

                auto tg_xxxxyy_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 90); 

                auto tg_xxxxyy_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 91); 

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

                // set up pointers to integrals

                auto tg_xxxxxxx_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx); 

                auto tg_xxxxxxx_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1); 

                auto tg_xxxxxxx_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 2); 

                auto tg_xxxxxxx_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 3); 

                auto tg_xxxxxxx_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 4); 

                auto tg_xxxxxxx_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 5); 

                auto tg_xxxxxxx_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 6); 

                auto tg_xxxxxxx_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 7); 

                auto tg_xxxxxxx_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 8); 

                auto tg_xxxxxxx_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 9); 

                auto tg_xxxxxxx_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 10); 

                auto tg_xxxxxxx_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 11); 

                auto tg_xxxxxxx_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 12); 

                auto tg_xxxxxxx_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 13); 

                auto tg_xxxxxxx_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 14); 

                auto tg_xxxxxxx_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 15); 

                auto tg_xxxxxxx_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 16); 

                auto tg_xxxxxxx_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 17); 

                auto tg_xxxxxxx_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 18); 

                auto tg_xxxxxxx_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 19); 

                auto tg_xxxxxxx_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 20); 

                auto tg_xxxxxxx_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 21); 

                auto tg_xxxxxxx_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 22); 

                auto tg_xxxxxxx_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 23); 

                auto tg_xxxxxxx_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 24); 

                auto tg_xxxxxxx_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 25); 

                auto tg_xxxxxxx_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 26); 

                auto tg_xxxxxxx_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 27); 

                auto tg_xxxxxxy_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 28); 

                auto tg_xxxxxxy_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 29); 

                auto tg_xxxxxxy_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 30); 

                auto tg_xxxxxxy_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 31); 

                auto tg_xxxxxxy_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 32); 

                auto tg_xxxxxxy_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 33); 

                auto tg_xxxxxxy_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 34); 

                auto tg_xxxxxxy_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 35); 

                auto tg_xxxxxxy_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 36); 

                auto tg_xxxxxxy_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 37); 

                auto tg_xxxxxxy_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 38); 

                auto tg_xxxxxxy_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 39); 

                auto tg_xxxxxxy_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 40); 

                auto tg_xxxxxxy_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 41); 

                auto tg_xxxxxxy_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 42); 

                auto tg_xxxxxxy_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 43); 

                auto tg_xxxxxxy_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 44); 

                auto tg_xxxxxxy_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 45); 

                auto tg_xxxxxxy_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 46); 

                auto tg_xxxxxxy_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 47); 

                auto tg_xxxxxxy_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 48); 

                auto tg_xxxxxxy_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 49); 

                auto tg_xxxxxxy_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 50); 

                auto tg_xxxxxxy_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 51); 

                auto tg_xxxxxxy_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 52); 

                auto tg_xxxxxxy_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 53); 

                auto tg_xxxxxxy_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 54); 

                auto tg_xxxxxxy_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 55); 

                auto tg_xxxxxxz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 56); 

                auto tg_xxxxxxz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 57); 

                auto tg_xxxxxxz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 58); 

                auto tg_xxxxxxz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 59); 

                auto tg_xxxxxxz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 60); 

                auto tg_xxxxxxz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 61); 

                auto tg_xxxxxxz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 62); 

                auto tg_xxxxxxz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 63); 

                auto tg_xxxxxxz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 64); 

                auto tg_xxxxxxz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 65); 

                auto tg_xxxxxxz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 66); 

                auto tg_xxxxxxz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 67); 

                auto tg_xxxxxxz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 68); 

                auto tg_xxxxxxz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 69); 

                auto tg_xxxxxxz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 70); 

                auto tg_xxxxxxz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 71); 

                auto tg_xxxxxxz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 72); 

                auto tg_xxxxxxz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 73); 

                auto tg_xxxxxxz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 74); 

                auto tg_xxxxxxz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 75); 

                auto tg_xxxxxxz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 76); 

                auto tg_xxxxxxz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 77); 

                auto tg_xxxxxxz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 78); 

                auto tg_xxxxxxz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 79); 

                auto tg_xxxxxxz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 80); 

                auto tg_xxxxxxz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 81); 

                auto tg_xxxxxxz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 82); 

                auto tg_xxxxxxz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 83); 

                auto tg_xxxxxyy_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 84); 

                auto tg_xxxxxyy_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 85); 

                auto tg_xxxxxyy_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 86); 

                auto tg_xxxxxyy_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 87); 

                auto tg_xxxxxyy_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 88); 

                auto tg_xxxxxyy_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 89); 

                auto tg_xxxxxyy_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 90); 

                auto tg_xxxxxyy_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 91); 

                // Batch of Integrals (0,92)

                #pragma omp simd aligned(fxn, fza, tg_xxxxx_xxxxxx_0, tg_xxxxx_xxxxxx_1, tg_xxxxx_xxxxxy_0, \
                                         tg_xxxxx_xxxxxy_1, tg_xxxxx_xxxxxz_0, tg_xxxxx_xxxxxz_1, tg_xxxxx_xxxxyy_0, \
                                         tg_xxxxx_xxxxyy_1, tg_xxxxx_xxxxyz_0, tg_xxxxx_xxxxyz_1, tg_xxxxx_xxxxzz_0, \
                                         tg_xxxxx_xxxxzz_1, tg_xxxxx_xxxyyy_0, tg_xxxxx_xxxyyy_1, tg_xxxxx_xxxyyz_0, \
                                         tg_xxxxx_xxxyyz_1, tg_xxxxx_xxxyzz_0, tg_xxxxx_xxxyzz_1, tg_xxxxx_xxxzzz_0, \
                                         tg_xxxxx_xxxzzz_1, tg_xxxxx_xxyyyy_0, tg_xxxxx_xxyyyy_1, tg_xxxxx_xxyyyz_0, \
                                         tg_xxxxx_xxyyyz_1, tg_xxxxx_xxyyzz_0, tg_xxxxx_xxyyzz_1, tg_xxxxx_xxyzzz_0, \
                                         tg_xxxxx_xxyzzz_1, tg_xxxxx_xxzzzz_0, tg_xxxxx_xxzzzz_1, tg_xxxxx_xyyyyy_0, \
                                         tg_xxxxx_xyyyyy_1, tg_xxxxx_xyyyyz_0, tg_xxxxx_xyyyyz_1, tg_xxxxx_xyyyzz_0, \
                                         tg_xxxxx_xyyyzz_1, tg_xxxxx_xyyzzz_0, tg_xxxxx_xyyzzz_1, tg_xxxxx_xyzzzz_0, \
                                         tg_xxxxx_xyzzzz_1, tg_xxxxx_xzzzzz_0, tg_xxxxx_xzzzzz_1, tg_xxxxx_yyyyyy_0, \
                                         tg_xxxxx_yyyyyy_1, tg_xxxxx_yyyyyz_0, tg_xxxxx_yyyyyz_1, tg_xxxxx_yyyyzz_0, \
                                         tg_xxxxx_yyyyzz_1, tg_xxxxx_yyyzzz_0, tg_xxxxx_yyyzzz_1, tg_xxxxx_yyzzzz_0, \
                                         tg_xxxxx_yyzzzz_1, tg_xxxxx_yzzzzz_0, tg_xxxxx_yzzzzz_1, tg_xxxxx_zzzzzz_0, \
                                         tg_xxxxx_zzzzzz_1, tg_xxxxxx_xxxxx_1, tg_xxxxxx_xxxxxx_0, tg_xxxxxx_xxxxxx_1, \
                                         tg_xxxxxx_xxxxxy_0, tg_xxxxxx_xxxxxy_1, tg_xxxxxx_xxxxxz_0, tg_xxxxxx_xxxxxz_1, \
                                         tg_xxxxxx_xxxxy_1, tg_xxxxxx_xxxxyy_0, tg_xxxxxx_xxxxyy_1, tg_xxxxxx_xxxxyz_0, \
                                         tg_xxxxxx_xxxxyz_1, tg_xxxxxx_xxxxz_1, tg_xxxxxx_xxxxzz_0, tg_xxxxxx_xxxxzz_1, \
                                         tg_xxxxxx_xxxyy_1, tg_xxxxxx_xxxyyy_0, tg_xxxxxx_xxxyyy_1, tg_xxxxxx_xxxyyz_0, \
                                         tg_xxxxxx_xxxyyz_1, tg_xxxxxx_xxxyz_1, tg_xxxxxx_xxxyzz_0, tg_xxxxxx_xxxyzz_1, \
                                         tg_xxxxxx_xxxzz_1, tg_xxxxxx_xxxzzz_0, tg_xxxxxx_xxxzzz_1, tg_xxxxxx_xxyyy_1, \
                                         tg_xxxxxx_xxyyyy_0, tg_xxxxxx_xxyyyy_1, tg_xxxxxx_xxyyyz_0, tg_xxxxxx_xxyyyz_1, \
                                         tg_xxxxxx_xxyyz_1, tg_xxxxxx_xxyyzz_0, tg_xxxxxx_xxyyzz_1, tg_xxxxxx_xxyzz_1, \
                                         tg_xxxxxx_xxyzzz_0, tg_xxxxxx_xxyzzz_1, tg_xxxxxx_xxzzz_1, tg_xxxxxx_xxzzzz_0, \
                                         tg_xxxxxx_xxzzzz_1, tg_xxxxxx_xyyyy_1, tg_xxxxxx_xyyyyy_0, tg_xxxxxx_xyyyyy_1, \
                                         tg_xxxxxx_xyyyyz_0, tg_xxxxxx_xyyyyz_1, tg_xxxxxx_xyyyz_1, tg_xxxxxx_xyyyzz_0, \
                                         tg_xxxxxx_xyyyzz_1, tg_xxxxxx_xyyzz_1, tg_xxxxxx_xyyzzz_0, tg_xxxxxx_xyyzzz_1, \
                                         tg_xxxxxx_xyzzz_1, tg_xxxxxx_xyzzzz_0, tg_xxxxxx_xyzzzz_1, tg_xxxxxx_xzzzz_1, \
                                         tg_xxxxxx_xzzzzz_0, tg_xxxxxx_xzzzzz_1, tg_xxxxxx_yyyyy_1, tg_xxxxxx_yyyyyy_0, \
                                         tg_xxxxxx_yyyyyy_1, tg_xxxxxx_yyyyyz_0, tg_xxxxxx_yyyyyz_1, tg_xxxxxx_yyyyz_1, \
                                         tg_xxxxxx_yyyyzz_0, tg_xxxxxx_yyyyzz_1, tg_xxxxxx_yyyzz_1, tg_xxxxxx_yyyzzz_0, \
                                         tg_xxxxxx_yyyzzz_1, tg_xxxxxx_yyzzz_1, tg_xxxxxx_yyzzzz_0, tg_xxxxxx_yyzzzz_1, \
                                         tg_xxxxxx_yzzzz_1, tg_xxxxxx_yzzzzz_0, tg_xxxxxx_yzzzzz_1, tg_xxxxxx_zzzzz_1, \
                                         tg_xxxxxx_zzzzzz_0, tg_xxxxxx_zzzzzz_1, tg_xxxxxxx_xxxxxx_0, tg_xxxxxxx_xxxxxy_0, \
                                         tg_xxxxxxx_xxxxxz_0, tg_xxxxxxx_xxxxyy_0, tg_xxxxxxx_xxxxyz_0, tg_xxxxxxx_xxxxzz_0, \
                                         tg_xxxxxxx_xxxyyy_0, tg_xxxxxxx_xxxyyz_0, tg_xxxxxxx_xxxyzz_0, tg_xxxxxxx_xxxzzz_0, \
                                         tg_xxxxxxx_xxyyyy_0, tg_xxxxxxx_xxyyyz_0, tg_xxxxxxx_xxyyzz_0, tg_xxxxxxx_xxyzzz_0, \
                                         tg_xxxxxxx_xxzzzz_0, tg_xxxxxxx_xyyyyy_0, tg_xxxxxxx_xyyyyz_0, tg_xxxxxxx_xyyyzz_0, \
                                         tg_xxxxxxx_xyyzzz_0, tg_xxxxxxx_xyzzzz_0, tg_xxxxxxx_xzzzzz_0, tg_xxxxxxx_yyyyyy_0, \
                                         tg_xxxxxxx_yyyyyz_0, tg_xxxxxxx_yyyyzz_0, tg_xxxxxxx_yyyzzz_0, tg_xxxxxxx_yyzzzz_0, \
                                         tg_xxxxxxx_yzzzzz_0, tg_xxxxxxx_zzzzzz_0, tg_xxxxxxy_xxxxxx_0, tg_xxxxxxy_xxxxxy_0, \
                                         tg_xxxxxxy_xxxxxz_0, tg_xxxxxxy_xxxxyy_0, tg_xxxxxxy_xxxxyz_0, tg_xxxxxxy_xxxxzz_0, \
                                         tg_xxxxxxy_xxxyyy_0, tg_xxxxxxy_xxxyyz_0, tg_xxxxxxy_xxxyzz_0, tg_xxxxxxy_xxxzzz_0, \
                                         tg_xxxxxxy_xxyyyy_0, tg_xxxxxxy_xxyyyz_0, tg_xxxxxxy_xxyyzz_0, tg_xxxxxxy_xxyzzz_0, \
                                         tg_xxxxxxy_xxzzzz_0, tg_xxxxxxy_xyyyyy_0, tg_xxxxxxy_xyyyyz_0, tg_xxxxxxy_xyyyzz_0, \
                                         tg_xxxxxxy_xyyzzz_0, tg_xxxxxxy_xyzzzz_0, tg_xxxxxxy_xzzzzz_0, tg_xxxxxxy_yyyyyy_0, \
                                         tg_xxxxxxy_yyyyyz_0, tg_xxxxxxy_yyyyzz_0, tg_xxxxxxy_yyyzzz_0, tg_xxxxxxy_yyzzzz_0, \
                                         tg_xxxxxxy_yzzzzz_0, tg_xxxxxxy_zzzzzz_0, tg_xxxxxxz_xxxxxx_0, tg_xxxxxxz_xxxxxy_0, \
                                         tg_xxxxxxz_xxxxxz_0, tg_xxxxxxz_xxxxyy_0, tg_xxxxxxz_xxxxyz_0, tg_xxxxxxz_xxxxzz_0, \
                                         tg_xxxxxxz_xxxyyy_0, tg_xxxxxxz_xxxyyz_0, tg_xxxxxxz_xxxyzz_0, tg_xxxxxxz_xxxzzz_0, \
                                         tg_xxxxxxz_xxyyyy_0, tg_xxxxxxz_xxyyyz_0, tg_xxxxxxz_xxyyzz_0, tg_xxxxxxz_xxyzzz_0, \
                                         tg_xxxxxxz_xxzzzz_0, tg_xxxxxxz_xyyyyy_0, tg_xxxxxxz_xyyyyz_0, tg_xxxxxxz_xyyyzz_0, \
                                         tg_xxxxxxz_xyyzzz_0, tg_xxxxxxz_xyzzzz_0, tg_xxxxxxz_xzzzzz_0, tg_xxxxxxz_yyyyyy_0, \
                                         tg_xxxxxxz_yyyyyz_0, tg_xxxxxxz_yyyyzz_0, tg_xxxxxxz_yyyzzz_0, tg_xxxxxxz_yyzzzz_0, \
                                         tg_xxxxxxz_yzzzzz_0, tg_xxxxxxz_zzzzzz_0, tg_xxxxxy_xxxxx_1, tg_xxxxxy_xxxxxx_0, \
                                         tg_xxxxxy_xxxxxx_1, tg_xxxxxy_xxxxxy_0, tg_xxxxxy_xxxxxy_1, tg_xxxxxy_xxxxxz_0, \
                                         tg_xxxxxy_xxxxxz_1, tg_xxxxxy_xxxxy_1, tg_xxxxxy_xxxxyy_0, tg_xxxxxy_xxxxyy_1, \
                                         tg_xxxxxy_xxxxyz_0, tg_xxxxxy_xxxxyz_1, tg_xxxxxy_xxxxz_1, tg_xxxxxy_xxxxzz_0, \
                                         tg_xxxxxy_xxxxzz_1, tg_xxxxxy_xxxyy_1, tg_xxxxxy_xxxyyy_0, tg_xxxxxy_xxxyyy_1, \
                                         tg_xxxxxy_xxxyyz_0, tg_xxxxxy_xxxyyz_1, tg_xxxxxy_xxxyz_1, tg_xxxxxy_xxxyzz_0, \
                                         tg_xxxxxy_xxxyzz_1, tg_xxxxxy_xxxzz_1, tg_xxxxxy_xxxzzz_0, tg_xxxxxy_xxxzzz_1, \
                                         tg_xxxxxy_xxyyy_1, tg_xxxxxy_xxyyyy_0, tg_xxxxxy_xxyyyy_1, tg_xxxxxy_xxyyyz_0, \
                                         tg_xxxxxy_xxyyyz_1, tg_xxxxxy_xxyyz_1, tg_xxxxxy_xxyyzz_0, tg_xxxxxy_xxyyzz_1, \
                                         tg_xxxxxy_xxyzz_1, tg_xxxxxy_xxyzzz_0, tg_xxxxxy_xxyzzz_1, tg_xxxxxy_xxzzz_1, \
                                         tg_xxxxxy_xxzzzz_0, tg_xxxxxy_xxzzzz_1, tg_xxxxxy_xyyyy_1, tg_xxxxxy_xyyyyy_0, \
                                         tg_xxxxxy_xyyyyy_1, tg_xxxxxy_xyyyyz_0, tg_xxxxxy_xyyyyz_1, tg_xxxxxy_xyyyz_1, \
                                         tg_xxxxxy_xyyyzz_0, tg_xxxxxy_xyyyzz_1, tg_xxxxxy_xyyzz_1, tg_xxxxxy_xyyzzz_0, \
                                         tg_xxxxxy_xyyzzz_1, tg_xxxxxy_xyzzz_1, tg_xxxxxy_xyzzzz_0, tg_xxxxxy_xyzzzz_1, \
                                         tg_xxxxxy_xzzzz_1, tg_xxxxxy_xzzzzz_0, tg_xxxxxy_xzzzzz_1, tg_xxxxxy_yyyyy_1, \
                                         tg_xxxxxy_yyyyyy_0, tg_xxxxxy_yyyyyy_1, tg_xxxxxy_yyyyyz_0, tg_xxxxxy_yyyyyz_1, \
                                         tg_xxxxxy_yyyyz_1, tg_xxxxxy_yyyyzz_0, tg_xxxxxy_yyyyzz_1, tg_xxxxxy_yyyzz_1, \
                                         tg_xxxxxy_yyyzzz_0, tg_xxxxxy_yyyzzz_1, tg_xxxxxy_yyzzz_1, tg_xxxxxy_yyzzzz_0, \
                                         tg_xxxxxy_yyzzzz_1, tg_xxxxxy_yzzzz_1, tg_xxxxxy_yzzzzz_0, tg_xxxxxy_yzzzzz_1, \
                                         tg_xxxxxy_zzzzz_1, tg_xxxxxy_zzzzzz_0, tg_xxxxxy_zzzzzz_1, tg_xxxxxyy_xxxxxx_0, \
                                         tg_xxxxxyy_xxxxxy_0, tg_xxxxxyy_xxxxxz_0, tg_xxxxxyy_xxxxyy_0, tg_xxxxxyy_xxxxyz_0, \
                                         tg_xxxxxyy_xxxxzz_0, tg_xxxxxyy_xxxyyy_0, tg_xxxxxyy_xxxyyz_0, tg_xxxxxz_xxxxx_1, \
                                         tg_xxxxxz_xxxxxx_0, tg_xxxxxz_xxxxxx_1, tg_xxxxxz_xxxxxy_0, tg_xxxxxz_xxxxxy_1, \
                                         tg_xxxxxz_xxxxxz_0, tg_xxxxxz_xxxxxz_1, tg_xxxxxz_xxxxy_1, tg_xxxxxz_xxxxyy_0, \
                                         tg_xxxxxz_xxxxyy_1, tg_xxxxxz_xxxxyz_0, tg_xxxxxz_xxxxyz_1, tg_xxxxxz_xxxxz_1, \
                                         tg_xxxxxz_xxxxzz_0, tg_xxxxxz_xxxxzz_1, tg_xxxxxz_xxxyy_1, tg_xxxxxz_xxxyyy_0, \
                                         tg_xxxxxz_xxxyyy_1, tg_xxxxxz_xxxyyz_0, tg_xxxxxz_xxxyyz_1, tg_xxxxxz_xxxyz_1, \
                                         tg_xxxxxz_xxxyzz_0, tg_xxxxxz_xxxyzz_1, tg_xxxxxz_xxxzz_1, tg_xxxxxz_xxxzzz_0, \
                                         tg_xxxxxz_xxxzzz_1, tg_xxxxxz_xxyyy_1, tg_xxxxxz_xxyyyy_0, tg_xxxxxz_xxyyyy_1, \
                                         tg_xxxxxz_xxyyyz_0, tg_xxxxxz_xxyyyz_1, tg_xxxxxz_xxyyz_1, tg_xxxxxz_xxyyzz_0, \
                                         tg_xxxxxz_xxyyzz_1, tg_xxxxxz_xxyzz_1, tg_xxxxxz_xxyzzz_0, tg_xxxxxz_xxyzzz_1, \
                                         tg_xxxxxz_xxzzz_1, tg_xxxxxz_xxzzzz_0, tg_xxxxxz_xxzzzz_1, tg_xxxxxz_xyyyy_1, \
                                         tg_xxxxxz_xyyyyy_0, tg_xxxxxz_xyyyyy_1, tg_xxxxxz_xyyyyz_0, tg_xxxxxz_xyyyyz_1, \
                                         tg_xxxxxz_xyyyz_1, tg_xxxxxz_xyyyzz_0, tg_xxxxxz_xyyyzz_1, tg_xxxxxz_xyyzz_1, \
                                         tg_xxxxxz_xyyzzz_0, tg_xxxxxz_xyyzzz_1, tg_xxxxxz_xyzzz_1, tg_xxxxxz_xyzzzz_0, \
                                         tg_xxxxxz_xyzzzz_1, tg_xxxxxz_xzzzz_1, tg_xxxxxz_xzzzzz_0, tg_xxxxxz_xzzzzz_1, \
                                         tg_xxxxxz_yyyyy_1, tg_xxxxxz_yyyyyy_0, tg_xxxxxz_yyyyyy_1, tg_xxxxxz_yyyyyz_0, \
                                         tg_xxxxxz_yyyyyz_1, tg_xxxxxz_yyyyz_1, tg_xxxxxz_yyyyzz_0, tg_xxxxxz_yyyyzz_1, \
                                         tg_xxxxxz_yyyzz_1, tg_xxxxxz_yyyzzz_0, tg_xxxxxz_yyyzzz_1, tg_xxxxxz_yyzzz_1, \
                                         tg_xxxxxz_yyzzzz_0, tg_xxxxxz_yyzzzz_1, tg_xxxxxz_yzzzz_1, tg_xxxxxz_yzzzzz_0, \
                                         tg_xxxxxz_yzzzzz_1, tg_xxxxxz_zzzzz_1, tg_xxxxxz_zzzzzz_0, tg_xxxxxz_zzzzzz_1, \
                                         tg_xxxxy_xxxxxx_0, tg_xxxxy_xxxxxx_1, tg_xxxxy_xxxxxy_0, tg_xxxxy_xxxxxy_1, \
                                         tg_xxxxy_xxxxxz_0, tg_xxxxy_xxxxxz_1, tg_xxxxy_xxxxyy_0, tg_xxxxy_xxxxyy_1, \
                                         tg_xxxxy_xxxxyz_0, tg_xxxxy_xxxxyz_1, tg_xxxxy_xxxxzz_0, tg_xxxxy_xxxxzz_1, \
                                         tg_xxxxy_xxxyyy_0, tg_xxxxy_xxxyyy_1, tg_xxxxy_xxxyyz_0, tg_xxxxy_xxxyyz_1, \
                                         tg_xxxxy_xxxyzz_0, tg_xxxxy_xxxyzz_1, tg_xxxxy_xxxzzz_0, tg_xxxxy_xxxzzz_1, \
                                         tg_xxxxy_xxyyyy_0, tg_xxxxy_xxyyyy_1, tg_xxxxy_xxyyyz_0, tg_xxxxy_xxyyyz_1, \
                                         tg_xxxxy_xxyyzz_0, tg_xxxxy_xxyyzz_1, tg_xxxxy_xxyzzz_0, tg_xxxxy_xxyzzz_1, \
                                         tg_xxxxy_xxzzzz_0, tg_xxxxy_xxzzzz_1, tg_xxxxy_xyyyyy_0, tg_xxxxy_xyyyyy_1, \
                                         tg_xxxxy_xyyyyz_0, tg_xxxxy_xyyyyz_1, tg_xxxxy_xyyyzz_0, tg_xxxxy_xyyyzz_1, \
                                         tg_xxxxy_xyyzzz_0, tg_xxxxy_xyyzzz_1, tg_xxxxy_xyzzzz_0, tg_xxxxy_xyzzzz_1, \
                                         tg_xxxxy_xzzzzz_0, tg_xxxxy_xzzzzz_1, tg_xxxxy_yyyyyy_0, tg_xxxxy_yyyyyy_1, \
                                         tg_xxxxy_yyyyyz_0, tg_xxxxy_yyyyyz_1, tg_xxxxy_yyyyzz_0, tg_xxxxy_yyyyzz_1, \
                                         tg_xxxxy_yyyzzz_0, tg_xxxxy_yyyzzz_1, tg_xxxxy_yyzzzz_0, tg_xxxxy_yyzzzz_1, \
                                         tg_xxxxy_yzzzzz_0, tg_xxxxy_yzzzzz_1, tg_xxxxy_zzzzzz_0, tg_xxxxy_zzzzzz_1, \
                                         tg_xxxxyy_xxxxx_1, tg_xxxxyy_xxxxxx_0, tg_xxxxyy_xxxxxx_1, tg_xxxxyy_xxxxxy_0, \
                                         tg_xxxxyy_xxxxxy_1, tg_xxxxyy_xxxxxz_0, tg_xxxxyy_xxxxxz_1, tg_xxxxyy_xxxxy_1, \
                                         tg_xxxxyy_xxxxyy_0, tg_xxxxyy_xxxxyy_1, tg_xxxxyy_xxxxyz_0, tg_xxxxyy_xxxxyz_1, \
                                         tg_xxxxyy_xxxxz_1, tg_xxxxyy_xxxxzz_0, tg_xxxxyy_xxxxzz_1, tg_xxxxyy_xxxyy_1, \
                                         tg_xxxxyy_xxxyyy_0, tg_xxxxyy_xxxyyy_1, tg_xxxxyy_xxxyyz_0, tg_xxxxyy_xxxyyz_1, \
                                         tg_xxxxyy_xxxyz_1, tg_xxxxyy_xxxzz_1, tg_xxxxyy_xxyyy_1, tg_xxxxyy_xxyyz_1, \
                                         tg_xxxxz_xxxxxx_0, tg_xxxxz_xxxxxx_1, tg_xxxxz_xxxxxy_0, tg_xxxxz_xxxxxy_1, \
                                         tg_xxxxz_xxxxxz_0, tg_xxxxz_xxxxxz_1, tg_xxxxz_xxxxyy_0, tg_xxxxz_xxxxyy_1, \
                                         tg_xxxxz_xxxxyz_0, tg_xxxxz_xxxxyz_1, tg_xxxxz_xxxxzz_0, tg_xxxxz_xxxxzz_1, \
                                         tg_xxxxz_xxxyyy_0, tg_xxxxz_xxxyyy_1, tg_xxxxz_xxxyyz_0, tg_xxxxz_xxxyyz_1, \
                                         tg_xxxxz_xxxyzz_0, tg_xxxxz_xxxyzz_1, tg_xxxxz_xxxzzz_0, tg_xxxxz_xxxzzz_1, \
                                         tg_xxxxz_xxyyyy_0, tg_xxxxz_xxyyyy_1, tg_xxxxz_xxyyyz_0, tg_xxxxz_xxyyyz_1, \
                                         tg_xxxxz_xxyyzz_0, tg_xxxxz_xxyyzz_1, tg_xxxxz_xxyzzz_0, tg_xxxxz_xxyzzz_1, \
                                         tg_xxxxz_xxzzzz_0, tg_xxxxz_xxzzzz_1, tg_xxxxz_xyyyyy_0, tg_xxxxz_xyyyyy_1, \
                                         tg_xxxxz_xyyyyz_0, tg_xxxxz_xyyyyz_1, tg_xxxxz_xyyyzz_0, tg_xxxxz_xyyyzz_1, \
                                         tg_xxxxz_xyyzzz_0, tg_xxxxz_xyyzzz_1, tg_xxxxz_xyzzzz_0, tg_xxxxz_xyzzzz_1, \
                                         tg_xxxxz_xzzzzz_0, tg_xxxxz_xzzzzz_1, tg_xxxxz_yyyyyy_0, tg_xxxxz_yyyyyy_1, \
                                         tg_xxxxz_yyyyyz_0, tg_xxxxz_yyyyyz_1, tg_xxxxz_yyyyzz_0, tg_xxxxz_yyyyzz_1, \
                                         tg_xxxxz_yyyzzz_0, tg_xxxxz_yyyzzz_1, tg_xxxxz_yyzzzz_0, tg_xxxxz_yyzzzz_1, \
                                         tg_xxxxz_yzzzzz_0, tg_xxxxz_yzzzzz_1, tg_xxxxz_zzzzzz_0, tg_xxxxz_zzzzzz_1, \
                                         tg_xxxyy_xxxxxx_0, tg_xxxyy_xxxxxx_1, tg_xxxyy_xxxxxy_0, tg_xxxyy_xxxxxy_1, \
                                         tg_xxxyy_xxxxxz_0, tg_xxxyy_xxxxxz_1, tg_xxxyy_xxxxyy_0, tg_xxxyy_xxxxyy_1, \
                                         tg_xxxyy_xxxxyz_0, tg_xxxyy_xxxxyz_1, tg_xxxyy_xxxxzz_0, tg_xxxyy_xxxxzz_1, \
                                         tg_xxxyy_xxxyyy_0, tg_xxxyy_xxxyyy_1, tg_xxxyy_xxxyyz_0, tg_xxxyy_xxxyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxxx_xxxxxx_0[j] = pb_x * tg_xxxxxx_xxxxxx_0[j] + fr * tg_xxxxxx_xxxxxx_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxxx_0[j] - tg_xxxxx_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxxx_xxxxx_1[j];

                    tg_xxxxxxx_xxxxxy_0[j] = pb_x * tg_xxxxxx_xxxxxy_0[j] + fr * tg_xxxxxx_xxxxxy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxxy_0[j] - tg_xxxxx_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxx_xxxxy_1[j];

                    tg_xxxxxxx_xxxxxz_0[j] = pb_x * tg_xxxxxx_xxxxxz_0[j] + fr * tg_xxxxxx_xxxxxz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxxz_0[j] - tg_xxxxx_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxx_xxxxz_1[j];

                    tg_xxxxxxx_xxxxyy_0[j] = pb_x * tg_xxxxxx_xxxxyy_0[j] + fr * tg_xxxxxx_xxxxyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxyy_0[j] - tg_xxxxx_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxx_xxxyy_1[j];

                    tg_xxxxxxx_xxxxyz_0[j] = pb_x * tg_xxxxxx_xxxxyz_0[j] + fr * tg_xxxxxx_xxxxyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxyz_0[j] - tg_xxxxx_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxx_xxxyz_1[j];

                    tg_xxxxxxx_xxxxzz_0[j] = pb_x * tg_xxxxxx_xxxxzz_0[j] + fr * tg_xxxxxx_xxxxzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxxzz_0[j] - tg_xxxxx_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxx_xxxzz_1[j];

                    tg_xxxxxxx_xxxyyy_0[j] = pb_x * tg_xxxxxx_xxxyyy_0[j] + fr * tg_xxxxxx_xxxyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxyyy_0[j] - tg_xxxxx_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxyyy_1[j];

                    tg_xxxxxxx_xxxyyz_0[j] = pb_x * tg_xxxxxx_xxxyyz_0[j] + fr * tg_xxxxxx_xxxyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxyyz_0[j] - tg_xxxxx_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxyyz_1[j];

                    tg_xxxxxxx_xxxyzz_0[j] = pb_x * tg_xxxxxx_xxxyzz_0[j] + fr * tg_xxxxxx_xxxyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxyzz_0[j] - tg_xxxxx_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxyzz_1[j];

                    tg_xxxxxxx_xxxzzz_0[j] = pb_x * tg_xxxxxx_xxxzzz_0[j] + fr * tg_xxxxxx_xxxzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxzzz_0[j] - tg_xxxxx_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxzzz_1[j];

                    tg_xxxxxxx_xxyyyy_0[j] = pb_x * tg_xxxxxx_xxyyyy_0[j] + fr * tg_xxxxxx_xxyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyyyy_0[j] - tg_xxxxx_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyyyy_1[j];

                    tg_xxxxxxx_xxyyyz_0[j] = pb_x * tg_xxxxxx_xxyyyz_0[j] + fr * tg_xxxxxx_xxyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyyyz_0[j] - tg_xxxxx_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyyyz_1[j];

                    tg_xxxxxxx_xxyyzz_0[j] = pb_x * tg_xxxxxx_xxyyzz_0[j] + fr * tg_xxxxxx_xxyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyyzz_0[j] - tg_xxxxx_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyyzz_1[j];

                    tg_xxxxxxx_xxyzzz_0[j] = pb_x * tg_xxxxxx_xxyzzz_0[j] + fr * tg_xxxxxx_xxyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyzzz_0[j] - tg_xxxxx_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyzzz_1[j];

                    tg_xxxxxxx_xxzzzz_0[j] = pb_x * tg_xxxxxx_xxzzzz_0[j] + fr * tg_xxxxxx_xxzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxzzzz_0[j] - tg_xxxxx_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xzzzz_1[j];

                    tg_xxxxxxx_xyyyyy_0[j] = pb_x * tg_xxxxxx_xyyyyy_0[j] + fr * tg_xxxxxx_xyyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyyyy_0[j] - tg_xxxxx_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyyyy_1[j];

                    tg_xxxxxxx_xyyyyz_0[j] = pb_x * tg_xxxxxx_xyyyyz_0[j] + fr * tg_xxxxxx_xyyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyyyz_0[j] - tg_xxxxx_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyyyz_1[j];

                    tg_xxxxxxx_xyyyzz_0[j] = pb_x * tg_xxxxxx_xyyyzz_0[j] + fr * tg_xxxxxx_xyyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyyzz_0[j] - tg_xxxxx_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyyzz_1[j];

                    tg_xxxxxxx_xyyzzz_0[j] = pb_x * tg_xxxxxx_xyyzzz_0[j] + fr * tg_xxxxxx_xyyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyzzz_0[j] - tg_xxxxx_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyzzz_1[j];

                    tg_xxxxxxx_xyzzzz_0[j] = pb_x * tg_xxxxxx_xyzzzz_0[j] + fr * tg_xxxxxx_xyzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyzzzz_0[j] - tg_xxxxx_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yzzzz_1[j];

                    tg_xxxxxxx_xzzzzz_0[j] = pb_x * tg_xxxxxx_xzzzzz_0[j] + fr * tg_xxxxxx_xzzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xzzzzz_0[j] - tg_xxxxx_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_zzzzz_1[j];

                    tg_xxxxxxx_yyyyyy_0[j] = pb_x * tg_xxxxxx_yyyyyy_0[j] + fr * tg_xxxxxx_yyyyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyyyy_0[j] - tg_xxxxx_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxxx_yyyyyz_0[j] = pb_x * tg_xxxxxx_yyyyyz_0[j] + fr * tg_xxxxxx_yyyyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyyyz_0[j] - tg_xxxxx_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxxx_yyyyzz_0[j] = pb_x * tg_xxxxxx_yyyyzz_0[j] + fr * tg_xxxxxx_yyyyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyyzz_0[j] - tg_xxxxx_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxxx_yyyzzz_0[j] = pb_x * tg_xxxxxx_yyyzzz_0[j] + fr * tg_xxxxxx_yyyzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyzzz_0[j] - tg_xxxxx_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxxx_yyzzzz_0[j] = pb_x * tg_xxxxxx_yyzzzz_0[j] + fr * tg_xxxxxx_yyzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyzzzz_0[j] - tg_xxxxx_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxxx_yzzzzz_0[j] = pb_x * tg_xxxxxx_yzzzzz_0[j] + fr * tg_xxxxxx_yzzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yzzzzz_0[j] - tg_xxxxx_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxxx_zzzzzz_0[j] = pb_x * tg_xxxxxx_zzzzzz_0[j] + fr * tg_xxxxxx_zzzzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_zzzzzz_0[j] - tg_xxxxx_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_xxxxxx_0[j] = pb_x * tg_xxxxxy_xxxxxx_0[j] + fr * tg_xxxxxy_xxxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxxx_0[j] - tg_xxxxy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxxy_xxxxx_1[j];

                    tg_xxxxxxy_xxxxxy_0[j] = pb_x * tg_xxxxxy_xxxxxy_0[j] + fr * tg_xxxxxy_xxxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxxy_0[j] - tg_xxxxy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxy_xxxxy_1[j];

                    tg_xxxxxxy_xxxxxz_0[j] = pb_x * tg_xxxxxy_xxxxxz_0[j] + fr * tg_xxxxxy_xxxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxxz_0[j] - tg_xxxxy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxy_xxxxz_1[j];

                    tg_xxxxxxy_xxxxyy_0[j] = pb_x * tg_xxxxxy_xxxxyy_0[j] + fr * tg_xxxxxy_xxxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxyy_0[j] - tg_xxxxy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxy_xxxyy_1[j];

                    tg_xxxxxxy_xxxxyz_0[j] = pb_x * tg_xxxxxy_xxxxyz_0[j] + fr * tg_xxxxxy_xxxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxyz_0[j] - tg_xxxxy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxy_xxxyz_1[j];

                    tg_xxxxxxy_xxxxzz_0[j] = pb_x * tg_xxxxxy_xxxxzz_0[j] + fr * tg_xxxxxy_xxxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxxzz_0[j] - tg_xxxxy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxy_xxxzz_1[j];

                    tg_xxxxxxy_xxxyyy_0[j] = pb_x * tg_xxxxxy_xxxyyy_0[j] + fr * tg_xxxxxy_xxxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxyyy_0[j] - tg_xxxxy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxyyy_1[j];

                    tg_xxxxxxy_xxxyyz_0[j] = pb_x * tg_xxxxxy_xxxyyz_0[j] + fr * tg_xxxxxy_xxxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxyyz_0[j] - tg_xxxxy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxyyz_1[j];

                    tg_xxxxxxy_xxxyzz_0[j] = pb_x * tg_xxxxxy_xxxyzz_0[j] + fr * tg_xxxxxy_xxxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxyzz_0[j] - tg_xxxxy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxyzz_1[j];

                    tg_xxxxxxy_xxxzzz_0[j] = pb_x * tg_xxxxxy_xxxzzz_0[j] + fr * tg_xxxxxy_xxxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxzzz_0[j] - tg_xxxxy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxzzz_1[j];

                    tg_xxxxxxy_xxyyyy_0[j] = pb_x * tg_xxxxxy_xxyyyy_0[j] + fr * tg_xxxxxy_xxyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyyyy_0[j] - tg_xxxxy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyyyy_1[j];

                    tg_xxxxxxy_xxyyyz_0[j] = pb_x * tg_xxxxxy_xxyyyz_0[j] + fr * tg_xxxxxy_xxyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyyyz_0[j] - tg_xxxxy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyyyz_1[j];

                    tg_xxxxxxy_xxyyzz_0[j] = pb_x * tg_xxxxxy_xxyyzz_0[j] + fr * tg_xxxxxy_xxyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyyzz_0[j] - tg_xxxxy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyyzz_1[j];

                    tg_xxxxxxy_xxyzzz_0[j] = pb_x * tg_xxxxxy_xxyzzz_0[j] + fr * tg_xxxxxy_xxyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyzzz_0[j] - tg_xxxxy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyzzz_1[j];

                    tg_xxxxxxy_xxzzzz_0[j] = pb_x * tg_xxxxxy_xxzzzz_0[j] + fr * tg_xxxxxy_xxzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxzzzz_0[j] - tg_xxxxy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xzzzz_1[j];

                    tg_xxxxxxy_xyyyyy_0[j] = pb_x * tg_xxxxxy_xyyyyy_0[j] + fr * tg_xxxxxy_xyyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyyyy_0[j] - tg_xxxxy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyyyy_1[j];

                    tg_xxxxxxy_xyyyyz_0[j] = pb_x * tg_xxxxxy_xyyyyz_0[j] + fr * tg_xxxxxy_xyyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyyyz_0[j] - tg_xxxxy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyyyz_1[j];

                    tg_xxxxxxy_xyyyzz_0[j] = pb_x * tg_xxxxxy_xyyyzz_0[j] + fr * tg_xxxxxy_xyyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyyzz_0[j] - tg_xxxxy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyyzz_1[j];

                    tg_xxxxxxy_xyyzzz_0[j] = pb_x * tg_xxxxxy_xyyzzz_0[j] + fr * tg_xxxxxy_xyyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyzzz_0[j] - tg_xxxxy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyzzz_1[j];

                    tg_xxxxxxy_xyzzzz_0[j] = pb_x * tg_xxxxxy_xyzzzz_0[j] + fr * tg_xxxxxy_xyzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyzzzz_0[j] - tg_xxxxy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yzzzz_1[j];

                    tg_xxxxxxy_xzzzzz_0[j] = pb_x * tg_xxxxxy_xzzzzz_0[j] + fr * tg_xxxxxy_xzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xzzzzz_0[j] - tg_xxxxy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_zzzzz_1[j];

                    tg_xxxxxxy_yyyyyy_0[j] = pb_x * tg_xxxxxy_yyyyyy_0[j] + fr * tg_xxxxxy_yyyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyyyy_0[j] - tg_xxxxy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxxy_yyyyyz_0[j] = pb_x * tg_xxxxxy_yyyyyz_0[j] + fr * tg_xxxxxy_yyyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyyyz_0[j] - tg_xxxxy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxxy_yyyyzz_0[j] = pb_x * tg_xxxxxy_yyyyzz_0[j] + fr * tg_xxxxxy_yyyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyyzz_0[j] - tg_xxxxy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxxy_yyyzzz_0[j] = pb_x * tg_xxxxxy_yyyzzz_0[j] + fr * tg_xxxxxy_yyyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyzzz_0[j] - tg_xxxxy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_yyzzzz_0[j] = pb_x * tg_xxxxxy_yyzzzz_0[j] + fr * tg_xxxxxy_yyzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyzzzz_0[j] - tg_xxxxy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_yzzzzz_0[j] = pb_x * tg_xxxxxy_yzzzzz_0[j] + fr * tg_xxxxxy_yzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yzzzzz_0[j] - tg_xxxxy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_zzzzzz_0[j] = pb_x * tg_xxxxxy_zzzzzz_0[j] + fr * tg_xxxxxy_zzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_zzzzzz_0[j] - tg_xxxxy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_xxxxxx_0[j] = pb_x * tg_xxxxxz_xxxxxx_0[j] + fr * tg_xxxxxz_xxxxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxxx_0[j] - tg_xxxxz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxxz_xxxxx_1[j];

                    tg_xxxxxxz_xxxxxy_0[j] = pb_x * tg_xxxxxz_xxxxxy_0[j] + fr * tg_xxxxxz_xxxxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxxy_0[j] - tg_xxxxz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxz_xxxxy_1[j];

                    tg_xxxxxxz_xxxxxz_0[j] = pb_x * tg_xxxxxz_xxxxxz_0[j] + fr * tg_xxxxxz_xxxxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxxz_0[j] - tg_xxxxz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxxz_xxxxz_1[j];

                    tg_xxxxxxz_xxxxyy_0[j] = pb_x * tg_xxxxxz_xxxxyy_0[j] + fr * tg_xxxxxz_xxxxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxyy_0[j] - tg_xxxxz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxz_xxxyy_1[j];

                    tg_xxxxxxz_xxxxyz_0[j] = pb_x * tg_xxxxxz_xxxxyz_0[j] + fr * tg_xxxxxz_xxxxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxyz_0[j] - tg_xxxxz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxz_xxxyz_1[j];

                    tg_xxxxxxz_xxxxzz_0[j] = pb_x * tg_xxxxxz_xxxxzz_0[j] + fr * tg_xxxxxz_xxxxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxxzz_0[j] - tg_xxxxz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxz_xxxzz_1[j];

                    tg_xxxxxxz_xxxyyy_0[j] = pb_x * tg_xxxxxz_xxxyyy_0[j] + fr * tg_xxxxxz_xxxyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxyyy_0[j] - tg_xxxxz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxyyy_1[j];

                    tg_xxxxxxz_xxxyyz_0[j] = pb_x * tg_xxxxxz_xxxyyz_0[j] + fr * tg_xxxxxz_xxxyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxyyz_0[j] - tg_xxxxz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxyyz_1[j];

                    tg_xxxxxxz_xxxyzz_0[j] = pb_x * tg_xxxxxz_xxxyzz_0[j] + fr * tg_xxxxxz_xxxyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxyzz_0[j] - tg_xxxxz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxyzz_1[j];

                    tg_xxxxxxz_xxxzzz_0[j] = pb_x * tg_xxxxxz_xxxzzz_0[j] + fr * tg_xxxxxz_xxxzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxzzz_0[j] - tg_xxxxz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxzzz_1[j];

                    tg_xxxxxxz_xxyyyy_0[j] = pb_x * tg_xxxxxz_xxyyyy_0[j] + fr * tg_xxxxxz_xxyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyyyy_0[j] - tg_xxxxz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyyyy_1[j];

                    tg_xxxxxxz_xxyyyz_0[j] = pb_x * tg_xxxxxz_xxyyyz_0[j] + fr * tg_xxxxxz_xxyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyyyz_0[j] - tg_xxxxz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyyyz_1[j];

                    tg_xxxxxxz_xxyyzz_0[j] = pb_x * tg_xxxxxz_xxyyzz_0[j] + fr * tg_xxxxxz_xxyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyyzz_0[j] - tg_xxxxz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyyzz_1[j];

                    tg_xxxxxxz_xxyzzz_0[j] = pb_x * tg_xxxxxz_xxyzzz_0[j] + fr * tg_xxxxxz_xxyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyzzz_0[j] - tg_xxxxz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyzzz_1[j];

                    tg_xxxxxxz_xxzzzz_0[j] = pb_x * tg_xxxxxz_xxzzzz_0[j] + fr * tg_xxxxxz_xxzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxzzzz_0[j] - tg_xxxxz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xzzzz_1[j];

                    tg_xxxxxxz_xyyyyy_0[j] = pb_x * tg_xxxxxz_xyyyyy_0[j] + fr * tg_xxxxxz_xyyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyyyy_0[j] - tg_xxxxz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyyyy_1[j];

                    tg_xxxxxxz_xyyyyz_0[j] = pb_x * tg_xxxxxz_xyyyyz_0[j] + fr * tg_xxxxxz_xyyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyyyz_0[j] - tg_xxxxz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyyyz_1[j];

                    tg_xxxxxxz_xyyyzz_0[j] = pb_x * tg_xxxxxz_xyyyzz_0[j] + fr * tg_xxxxxz_xyyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyyzz_0[j] - tg_xxxxz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyyzz_1[j];

                    tg_xxxxxxz_xyyzzz_0[j] = pb_x * tg_xxxxxz_xyyzzz_0[j] + fr * tg_xxxxxz_xyyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyzzz_0[j] - tg_xxxxz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyzzz_1[j];

                    tg_xxxxxxz_xyzzzz_0[j] = pb_x * tg_xxxxxz_xyzzzz_0[j] + fr * tg_xxxxxz_xyzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyzzzz_0[j] - tg_xxxxz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yzzzz_1[j];

                    tg_xxxxxxz_xzzzzz_0[j] = pb_x * tg_xxxxxz_xzzzzz_0[j] + fr * tg_xxxxxz_xzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xzzzzz_0[j] - tg_xxxxz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_zzzzz_1[j];

                    tg_xxxxxxz_yyyyyy_0[j] = pb_x * tg_xxxxxz_yyyyyy_0[j] + fr * tg_xxxxxz_yyyyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyyyy_0[j] - tg_xxxxz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxxz_yyyyyz_0[j] = pb_x * tg_xxxxxz_yyyyyz_0[j] + fr * tg_xxxxxz_yyyyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyyyz_0[j] - tg_xxxxz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxxz_yyyyzz_0[j] = pb_x * tg_xxxxxz_yyyyzz_0[j] + fr * tg_xxxxxz_yyyyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyyzz_0[j] - tg_xxxxz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxxz_yyyzzz_0[j] = pb_x * tg_xxxxxz_yyyzzz_0[j] + fr * tg_xxxxxz_yyyzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyzzz_0[j] - tg_xxxxz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_yyzzzz_0[j] = pb_x * tg_xxxxxz_yyzzzz_0[j] + fr * tg_xxxxxz_yyzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyzzzz_0[j] - tg_xxxxz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_yzzzzz_0[j] = pb_x * tg_xxxxxz_yzzzzz_0[j] + fr * tg_xxxxxz_yzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yzzzzz_0[j] - tg_xxxxz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_zzzzzz_0[j] = pb_x * tg_xxxxxz_zzzzzz_0[j] + fr * tg_xxxxxz_zzzzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_zzzzzz_0[j] - tg_xxxxz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_xxxxxx_0[j] = pb_x * tg_xxxxyy_xxxxxx_0[j] + fr * tg_xxxxyy_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxxx_0[j] - tg_xxxyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxyy_xxxxx_1[j];

                    tg_xxxxxyy_xxxxxy_0[j] = pb_x * tg_xxxxyy_xxxxxy_0[j] + fr * tg_xxxxyy_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxxy_0[j] - tg_xxxyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyy_xxxxy_1[j];

                    tg_xxxxxyy_xxxxxz_0[j] = pb_x * tg_xxxxyy_xxxxxz_0[j] + fr * tg_xxxxyy_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxxz_0[j] - tg_xxxyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyy_xxxxz_1[j];

                    tg_xxxxxyy_xxxxyy_0[j] = pb_x * tg_xxxxyy_xxxxyy_0[j] + fr * tg_xxxxyy_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxyy_0[j] - tg_xxxyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyy_xxxyy_1[j];

                    tg_xxxxxyy_xxxxyz_0[j] = pb_x * tg_xxxxyy_xxxxyz_0[j] + fr * tg_xxxxyy_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxyz_0[j] - tg_xxxyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyy_xxxyz_1[j];

                    tg_xxxxxyy_xxxxzz_0[j] = pb_x * tg_xxxxyy_xxxxzz_0[j] + fr * tg_xxxxyy_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxxzz_0[j] - tg_xxxyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyy_xxxzz_1[j];

                    tg_xxxxxyy_xxxyyy_0[j] = pb_x * tg_xxxxyy_xxxyyy_0[j] + fr * tg_xxxxyy_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxyyy_0[j] - tg_xxxyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxyyy_1[j];

                    tg_xxxxxyy_xxxyyz_0[j] = pb_x * tg_xxxxyy_xxxyyz_0[j] + fr * tg_xxxxyy_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxyyz_0[j] - tg_xxxyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_92_184(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (92,184)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto tg_xxxxyy_xxxyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 92); 

                auto tg_xxxxyy_xxxzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 93); 

                auto tg_xxxxyy_xxyyyy_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 94); 

                auto tg_xxxxyy_xxyyyz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 95); 

                auto tg_xxxxyy_xxyyzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 96); 

                auto tg_xxxxyy_xxyzzz_0 = primBuffer[pidx_g_6_6_m0].data(784 * idx + 97); 

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

                auto tg_xxxxyy_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 92); 

                auto tg_xxxxyy_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 93); 

                auto tg_xxxxyy_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 94); 

                auto tg_xxxxyy_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 95); 

                auto tg_xxxxyy_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 96); 

                auto tg_xxxxyy_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 97); 

                auto tg_xxxxyy_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 98); 

                auto tg_xxxxyy_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 99); 

                auto tg_xxxxyy_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 100); 

                auto tg_xxxxyy_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 101); 

                auto tg_xxxxyy_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 102); 

                auto tg_xxxxyy_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 103); 

                auto tg_xxxxyy_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 104); 

                auto tg_xxxxyy_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 105); 

                auto tg_xxxxyy_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 106); 

                auto tg_xxxxyy_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 107); 

                auto tg_xxxxyy_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 108); 

                auto tg_xxxxyy_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 109); 

                auto tg_xxxxyy_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 110); 

                auto tg_xxxxyy_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 111); 

                auto tg_xxxxyz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 112); 

                auto tg_xxxxyz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 113); 

                auto tg_xxxxyz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 114); 

                auto tg_xxxxyz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 115); 

                auto tg_xxxxyz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 116); 

                auto tg_xxxxyz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 117); 

                auto tg_xxxxyz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 118); 

                auto tg_xxxxyz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 119); 

                auto tg_xxxxyz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 120); 

                auto tg_xxxxyz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 121); 

                auto tg_xxxxyz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 122); 

                auto tg_xxxxyz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 123); 

                auto tg_xxxxyz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 124); 

                auto tg_xxxxyz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 125); 

                auto tg_xxxxyz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 126); 

                auto tg_xxxxyz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 127); 

                auto tg_xxxxyz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 128); 

                auto tg_xxxxyz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 129); 

                auto tg_xxxxyz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 130); 

                auto tg_xxxxyz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 131); 

                auto tg_xxxxyz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 132); 

                auto tg_xxxxyz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 133); 

                auto tg_xxxxyz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 134); 

                auto tg_xxxxyz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 135); 

                auto tg_xxxxyz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 136); 

                auto tg_xxxxyz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 137); 

                auto tg_xxxxyz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 138); 

                auto tg_xxxxyz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 139); 

                auto tg_xxxxzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 140); 

                auto tg_xxxxzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 141); 

                auto tg_xxxxzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 142); 

                auto tg_xxxxzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 143); 

                auto tg_xxxxzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 144); 

                auto tg_xxxxzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 145); 

                auto tg_xxxxzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 146); 

                auto tg_xxxxzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 147); 

                auto tg_xxxxzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 148); 

                auto tg_xxxxzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 149); 

                auto tg_xxxxzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 150); 

                auto tg_xxxxzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 151); 

                auto tg_xxxxzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 152); 

                auto tg_xxxxzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 153); 

                auto tg_xxxxzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 154); 

                auto tg_xxxxzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 155); 

                auto tg_xxxxzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 156); 

                auto tg_xxxxzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 157); 

                auto tg_xxxxzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 158); 

                auto tg_xxxxzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 159); 

                auto tg_xxxxzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 160); 

                auto tg_xxxxzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 161); 

                auto tg_xxxxzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 162); 

                auto tg_xxxxzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 163); 

                auto tg_xxxxzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 164); 

                auto tg_xxxxzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 165); 

                auto tg_xxxxzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 166); 

                auto tg_xxxxzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 167); 

                auto tg_xxxyyy_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 168); 

                auto tg_xxxyyy_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 169); 

                auto tg_xxxyyy_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 170); 

                auto tg_xxxyyy_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 171); 

                auto tg_xxxyyy_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 172); 

                auto tg_xxxyyy_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 173); 

                auto tg_xxxyyy_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 174); 

                auto tg_xxxyyy_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 175); 

                auto tg_xxxyyy_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 176); 

                auto tg_xxxyyy_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 177); 

                auto tg_xxxyyy_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 178); 

                auto tg_xxxyyy_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 179); 

                auto tg_xxxyyy_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 180); 

                auto tg_xxxyyy_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 181); 

                auto tg_xxxyyy_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 182); 

                auto tg_xxxyyy_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 183); 

                auto tg_xxxyy_xxxyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 92); 

                auto tg_xxxyy_xxxzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 93); 

                auto tg_xxxyy_xxyyyy_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 94); 

                auto tg_xxxyy_xxyyyz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 95); 

                auto tg_xxxyy_xxyyzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 96); 

                auto tg_xxxyy_xxyzzz_0 = primBuffer[pidx_g_5_6_m0].data(588 * idx + 97); 

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

                auto tg_xxxyy_xxxyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 92); 

                auto tg_xxxyy_xxxzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 93); 

                auto tg_xxxyy_xxyyyy_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 94); 

                auto tg_xxxyy_xxyyyz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 95); 

                auto tg_xxxyy_xxyyzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 96); 

                auto tg_xxxyy_xxyzzz_1 = primBuffer[pidx_g_5_6_m1].data(588 * idx + 97); 

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

                // set up pointers to integrals

                auto tg_xxxxxyy_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 92); 

                auto tg_xxxxxyy_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 93); 

                auto tg_xxxxxyy_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 94); 

                auto tg_xxxxxyy_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 95); 

                auto tg_xxxxxyy_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 96); 

                auto tg_xxxxxyy_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 97); 

                auto tg_xxxxxyy_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 98); 

                auto tg_xxxxxyy_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 99); 

                auto tg_xxxxxyy_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 100); 

                auto tg_xxxxxyy_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 101); 

                auto tg_xxxxxyy_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 102); 

                auto tg_xxxxxyy_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 103); 

                auto tg_xxxxxyy_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 104); 

                auto tg_xxxxxyy_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 105); 

                auto tg_xxxxxyy_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 106); 

                auto tg_xxxxxyy_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 107); 

                auto tg_xxxxxyy_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 108); 

                auto tg_xxxxxyy_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 109); 

                auto tg_xxxxxyy_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 110); 

                auto tg_xxxxxyy_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 111); 

                auto tg_xxxxxyz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 112); 

                auto tg_xxxxxyz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 113); 

                auto tg_xxxxxyz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 114); 

                auto tg_xxxxxyz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 115); 

                auto tg_xxxxxyz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 116); 

                auto tg_xxxxxyz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 117); 

                auto tg_xxxxxyz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 118); 

                auto tg_xxxxxyz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 119); 

                auto tg_xxxxxyz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 120); 

                auto tg_xxxxxyz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 121); 

                auto tg_xxxxxyz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 122); 

                auto tg_xxxxxyz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 123); 

                auto tg_xxxxxyz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 124); 

                auto tg_xxxxxyz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 125); 

                auto tg_xxxxxyz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 126); 

                auto tg_xxxxxyz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 127); 

                auto tg_xxxxxyz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 128); 

                auto tg_xxxxxyz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 129); 

                auto tg_xxxxxyz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 130); 

                auto tg_xxxxxyz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 131); 

                auto tg_xxxxxyz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 132); 

                auto tg_xxxxxyz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 133); 

                auto tg_xxxxxyz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 134); 

                auto tg_xxxxxyz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 135); 

                auto tg_xxxxxyz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 136); 

                auto tg_xxxxxyz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 137); 

                auto tg_xxxxxyz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 138); 

                auto tg_xxxxxyz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 139); 

                auto tg_xxxxxzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 140); 

                auto tg_xxxxxzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 141); 

                auto tg_xxxxxzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 142); 

                auto tg_xxxxxzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 143); 

                auto tg_xxxxxzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 144); 

                auto tg_xxxxxzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 145); 

                auto tg_xxxxxzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 146); 

                auto tg_xxxxxzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 147); 

                auto tg_xxxxxzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 148); 

                auto tg_xxxxxzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 149); 

                auto tg_xxxxxzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 150); 

                auto tg_xxxxxzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 151); 

                auto tg_xxxxxzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 152); 

                auto tg_xxxxxzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 153); 

                auto tg_xxxxxzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 154); 

                auto tg_xxxxxzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 155); 

                auto tg_xxxxxzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 156); 

                auto tg_xxxxxzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 157); 

                auto tg_xxxxxzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 158); 

                auto tg_xxxxxzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 159); 

                auto tg_xxxxxzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 160); 

                auto tg_xxxxxzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 161); 

                auto tg_xxxxxzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 162); 

                auto tg_xxxxxzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 163); 

                auto tg_xxxxxzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 164); 

                auto tg_xxxxxzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 165); 

                auto tg_xxxxxzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 166); 

                auto tg_xxxxxzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 167); 

                auto tg_xxxxyyy_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 168); 

                auto tg_xxxxyyy_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 169); 

                auto tg_xxxxyyy_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 170); 

                auto tg_xxxxyyy_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 171); 

                auto tg_xxxxyyy_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 172); 

                auto tg_xxxxyyy_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 173); 

                auto tg_xxxxyyy_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 174); 

                auto tg_xxxxyyy_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 175); 

                auto tg_xxxxyyy_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 176); 

                auto tg_xxxxyyy_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 177); 

                auto tg_xxxxyyy_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 178); 

                auto tg_xxxxyyy_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 179); 

                auto tg_xxxxyyy_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 180); 

                auto tg_xxxxyyy_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 181); 

                auto tg_xxxxyyy_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 182); 

                auto tg_xxxxyyy_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 183); 

                // Batch of Integrals (92,184)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxyy_xxxyzz_0, tg_xxxxxyy_xxxzzz_0, \
                                         tg_xxxxxyy_xxyyyy_0, tg_xxxxxyy_xxyyyz_0, tg_xxxxxyy_xxyyzz_0, tg_xxxxxyy_xxyzzz_0, \
                                         tg_xxxxxyy_xxzzzz_0, tg_xxxxxyy_xyyyyy_0, tg_xxxxxyy_xyyyyz_0, tg_xxxxxyy_xyyyzz_0, \
                                         tg_xxxxxyy_xyyzzz_0, tg_xxxxxyy_xyzzzz_0, tg_xxxxxyy_xzzzzz_0, tg_xxxxxyy_yyyyyy_0, \
                                         tg_xxxxxyy_yyyyyz_0, tg_xxxxxyy_yyyyzz_0, tg_xxxxxyy_yyyzzz_0, tg_xxxxxyy_yyzzzz_0, \
                                         tg_xxxxxyy_yzzzzz_0, tg_xxxxxyy_zzzzzz_0, tg_xxxxxyz_xxxxxx_0, tg_xxxxxyz_xxxxxy_0, \
                                         tg_xxxxxyz_xxxxxz_0, tg_xxxxxyz_xxxxyy_0, tg_xxxxxyz_xxxxyz_0, tg_xxxxxyz_xxxxzz_0, \
                                         tg_xxxxxyz_xxxyyy_0, tg_xxxxxyz_xxxyyz_0, tg_xxxxxyz_xxxyzz_0, tg_xxxxxyz_xxxzzz_0, \
                                         tg_xxxxxyz_xxyyyy_0, tg_xxxxxyz_xxyyyz_0, tg_xxxxxyz_xxyyzz_0, tg_xxxxxyz_xxyzzz_0, \
                                         tg_xxxxxyz_xxzzzz_0, tg_xxxxxyz_xyyyyy_0, tg_xxxxxyz_xyyyyz_0, tg_xxxxxyz_xyyyzz_0, \
                                         tg_xxxxxyz_xyyzzz_0, tg_xxxxxyz_xyzzzz_0, tg_xxxxxyz_xzzzzz_0, tg_xxxxxyz_yyyyyy_0, \
                                         tg_xxxxxyz_yyyyyz_0, tg_xxxxxyz_yyyyzz_0, tg_xxxxxyz_yyyzzz_0, tg_xxxxxyz_yyzzzz_0, \
                                         tg_xxxxxyz_yzzzzz_0, tg_xxxxxyz_zzzzzz_0, tg_xxxxxzz_xxxxxx_0, tg_xxxxxzz_xxxxxy_0, \
                                         tg_xxxxxzz_xxxxxz_0, tg_xxxxxzz_xxxxyy_0, tg_xxxxxzz_xxxxyz_0, tg_xxxxxzz_xxxxzz_0, \
                                         tg_xxxxxzz_xxxyyy_0, tg_xxxxxzz_xxxyyz_0, tg_xxxxxzz_xxxyzz_0, tg_xxxxxzz_xxxzzz_0, \
                                         tg_xxxxxzz_xxyyyy_0, tg_xxxxxzz_xxyyyz_0, tg_xxxxxzz_xxyyzz_0, tg_xxxxxzz_xxyzzz_0, \
                                         tg_xxxxxzz_xxzzzz_0, tg_xxxxxzz_xyyyyy_0, tg_xxxxxzz_xyyyyz_0, tg_xxxxxzz_xyyyzz_0, \
                                         tg_xxxxxzz_xyyzzz_0, tg_xxxxxzz_xyzzzz_0, tg_xxxxxzz_xzzzzz_0, tg_xxxxxzz_yyyyyy_0, \
                                         tg_xxxxxzz_yyyyyz_0, tg_xxxxxzz_yyyyzz_0, tg_xxxxxzz_yyyzzz_0, tg_xxxxxzz_yyzzzz_0, \
                                         tg_xxxxxzz_yzzzzz_0, tg_xxxxxzz_zzzzzz_0, tg_xxxxyy_xxxyzz_0, tg_xxxxyy_xxxyzz_1, \
                                         tg_xxxxyy_xxxzzz_0, tg_xxxxyy_xxxzzz_1, tg_xxxxyy_xxyyyy_0, tg_xxxxyy_xxyyyy_1, \
                                         tg_xxxxyy_xxyyyz_0, tg_xxxxyy_xxyyyz_1, tg_xxxxyy_xxyyzz_0, tg_xxxxyy_xxyyzz_1, \
                                         tg_xxxxyy_xxyzz_1, tg_xxxxyy_xxyzzz_0, tg_xxxxyy_xxyzzz_1, tg_xxxxyy_xxzzz_1, \
                                         tg_xxxxyy_xxzzzz_0, tg_xxxxyy_xxzzzz_1, tg_xxxxyy_xyyyy_1, tg_xxxxyy_xyyyyy_0, \
                                         tg_xxxxyy_xyyyyy_1, tg_xxxxyy_xyyyyz_0, tg_xxxxyy_xyyyyz_1, tg_xxxxyy_xyyyz_1, \
                                         tg_xxxxyy_xyyyzz_0, tg_xxxxyy_xyyyzz_1, tg_xxxxyy_xyyzz_1, tg_xxxxyy_xyyzzz_0, \
                                         tg_xxxxyy_xyyzzz_1, tg_xxxxyy_xyzzz_1, tg_xxxxyy_xyzzzz_0, tg_xxxxyy_xyzzzz_1, \
                                         tg_xxxxyy_xzzzz_1, tg_xxxxyy_xzzzzz_0, tg_xxxxyy_xzzzzz_1, tg_xxxxyy_yyyyy_1, \
                                         tg_xxxxyy_yyyyyy_0, tg_xxxxyy_yyyyyy_1, tg_xxxxyy_yyyyyz_0, tg_xxxxyy_yyyyyz_1, \
                                         tg_xxxxyy_yyyyz_1, tg_xxxxyy_yyyyzz_0, tg_xxxxyy_yyyyzz_1, tg_xxxxyy_yyyzz_1, \
                                         tg_xxxxyy_yyyzzz_0, tg_xxxxyy_yyyzzz_1, tg_xxxxyy_yyzzz_1, tg_xxxxyy_yyzzzz_0, \
                                         tg_xxxxyy_yyzzzz_1, tg_xxxxyy_yzzzz_1, tg_xxxxyy_yzzzzz_0, tg_xxxxyy_yzzzzz_1, \
                                         tg_xxxxyy_zzzzz_1, tg_xxxxyy_zzzzzz_0, tg_xxxxyy_zzzzzz_1, tg_xxxxyyy_xxxxxx_0, \
                                         tg_xxxxyyy_xxxxxy_0, tg_xxxxyyy_xxxxxz_0, tg_xxxxyyy_xxxxyy_0, tg_xxxxyyy_xxxxyz_0, \
                                         tg_xxxxyyy_xxxxzz_0, tg_xxxxyyy_xxxyyy_0, tg_xxxxyyy_xxxyyz_0, tg_xxxxyyy_xxxyzz_0, \
                                         tg_xxxxyyy_xxxzzz_0, tg_xxxxyyy_xxyyyy_0, tg_xxxxyyy_xxyyyz_0, tg_xxxxyyy_xxyyzz_0, \
                                         tg_xxxxyyy_xxyzzz_0, tg_xxxxyyy_xxzzzz_0, tg_xxxxyyy_xyyyyy_0, tg_xxxxyz_xxxxx_1, \
                                         tg_xxxxyz_xxxxxx_0, tg_xxxxyz_xxxxxx_1, tg_xxxxyz_xxxxxy_0, tg_xxxxyz_xxxxxy_1, \
                                         tg_xxxxyz_xxxxxz_0, tg_xxxxyz_xxxxxz_1, tg_xxxxyz_xxxxy_1, tg_xxxxyz_xxxxyy_0, \
                                         tg_xxxxyz_xxxxyy_1, tg_xxxxyz_xxxxyz_0, tg_xxxxyz_xxxxyz_1, tg_xxxxyz_xxxxz_1, \
                                         tg_xxxxyz_xxxxzz_0, tg_xxxxyz_xxxxzz_1, tg_xxxxyz_xxxyy_1, tg_xxxxyz_xxxyyy_0, \
                                         tg_xxxxyz_xxxyyy_1, tg_xxxxyz_xxxyyz_0, tg_xxxxyz_xxxyyz_1, tg_xxxxyz_xxxyz_1, \
                                         tg_xxxxyz_xxxyzz_0, tg_xxxxyz_xxxyzz_1, tg_xxxxyz_xxxzz_1, tg_xxxxyz_xxxzzz_0, \
                                         tg_xxxxyz_xxxzzz_1, tg_xxxxyz_xxyyy_1, tg_xxxxyz_xxyyyy_0, tg_xxxxyz_xxyyyy_1, \
                                         tg_xxxxyz_xxyyyz_0, tg_xxxxyz_xxyyyz_1, tg_xxxxyz_xxyyz_1, tg_xxxxyz_xxyyzz_0, \
                                         tg_xxxxyz_xxyyzz_1, tg_xxxxyz_xxyzz_1, tg_xxxxyz_xxyzzz_0, tg_xxxxyz_xxyzzz_1, \
                                         tg_xxxxyz_xxzzz_1, tg_xxxxyz_xxzzzz_0, tg_xxxxyz_xxzzzz_1, tg_xxxxyz_xyyyy_1, \
                                         tg_xxxxyz_xyyyyy_0, tg_xxxxyz_xyyyyy_1, tg_xxxxyz_xyyyyz_0, tg_xxxxyz_xyyyyz_1, \
                                         tg_xxxxyz_xyyyz_1, tg_xxxxyz_xyyyzz_0, tg_xxxxyz_xyyyzz_1, tg_xxxxyz_xyyzz_1, \
                                         tg_xxxxyz_xyyzzz_0, tg_xxxxyz_xyyzzz_1, tg_xxxxyz_xyzzz_1, tg_xxxxyz_xyzzzz_0, \
                                         tg_xxxxyz_xyzzzz_1, tg_xxxxyz_xzzzz_1, tg_xxxxyz_xzzzzz_0, tg_xxxxyz_xzzzzz_1, \
                                         tg_xxxxyz_yyyyy_1, tg_xxxxyz_yyyyyy_0, tg_xxxxyz_yyyyyy_1, tg_xxxxyz_yyyyyz_0, \
                                         tg_xxxxyz_yyyyyz_1, tg_xxxxyz_yyyyz_1, tg_xxxxyz_yyyyzz_0, tg_xxxxyz_yyyyzz_1, \
                                         tg_xxxxyz_yyyzz_1, tg_xxxxyz_yyyzzz_0, tg_xxxxyz_yyyzzz_1, tg_xxxxyz_yyzzz_1, \
                                         tg_xxxxyz_yyzzzz_0, tg_xxxxyz_yyzzzz_1, tg_xxxxyz_yzzzz_1, tg_xxxxyz_yzzzzz_0, \
                                         tg_xxxxyz_yzzzzz_1, tg_xxxxyz_zzzzz_1, tg_xxxxyz_zzzzzz_0, tg_xxxxyz_zzzzzz_1, \
                                         tg_xxxxzz_xxxxx_1, tg_xxxxzz_xxxxxx_0, tg_xxxxzz_xxxxxx_1, tg_xxxxzz_xxxxxy_0, \
                                         tg_xxxxzz_xxxxxy_1, tg_xxxxzz_xxxxxz_0, tg_xxxxzz_xxxxxz_1, tg_xxxxzz_xxxxy_1, \
                                         tg_xxxxzz_xxxxyy_0, tg_xxxxzz_xxxxyy_1, tg_xxxxzz_xxxxyz_0, tg_xxxxzz_xxxxyz_1, \
                                         tg_xxxxzz_xxxxz_1, tg_xxxxzz_xxxxzz_0, tg_xxxxzz_xxxxzz_1, tg_xxxxzz_xxxyy_1, \
                                         tg_xxxxzz_xxxyyy_0, tg_xxxxzz_xxxyyy_1, tg_xxxxzz_xxxyyz_0, tg_xxxxzz_xxxyyz_1, \
                                         tg_xxxxzz_xxxyz_1, tg_xxxxzz_xxxyzz_0, tg_xxxxzz_xxxyzz_1, tg_xxxxzz_xxxzz_1, \
                                         tg_xxxxzz_xxxzzz_0, tg_xxxxzz_xxxzzz_1, tg_xxxxzz_xxyyy_1, tg_xxxxzz_xxyyyy_0, \
                                         tg_xxxxzz_xxyyyy_1, tg_xxxxzz_xxyyyz_0, tg_xxxxzz_xxyyyz_1, tg_xxxxzz_xxyyz_1, \
                                         tg_xxxxzz_xxyyzz_0, tg_xxxxzz_xxyyzz_1, tg_xxxxzz_xxyzz_1, tg_xxxxzz_xxyzzz_0, \
                                         tg_xxxxzz_xxyzzz_1, tg_xxxxzz_xxzzz_1, tg_xxxxzz_xxzzzz_0, tg_xxxxzz_xxzzzz_1, \
                                         tg_xxxxzz_xyyyy_1, tg_xxxxzz_xyyyyy_0, tg_xxxxzz_xyyyyy_1, tg_xxxxzz_xyyyyz_0, \
                                         tg_xxxxzz_xyyyyz_1, tg_xxxxzz_xyyyz_1, tg_xxxxzz_xyyyzz_0, tg_xxxxzz_xyyyzz_1, \
                                         tg_xxxxzz_xyyzz_1, tg_xxxxzz_xyyzzz_0, tg_xxxxzz_xyyzzz_1, tg_xxxxzz_xyzzz_1, \
                                         tg_xxxxzz_xyzzzz_0, tg_xxxxzz_xyzzzz_1, tg_xxxxzz_xzzzz_1, tg_xxxxzz_xzzzzz_0, \
                                         tg_xxxxzz_xzzzzz_1, tg_xxxxzz_yyyyy_1, tg_xxxxzz_yyyyyy_0, tg_xxxxzz_yyyyyy_1, \
                                         tg_xxxxzz_yyyyyz_0, tg_xxxxzz_yyyyyz_1, tg_xxxxzz_yyyyz_1, tg_xxxxzz_yyyyzz_0, \
                                         tg_xxxxzz_yyyyzz_1, tg_xxxxzz_yyyzz_1, tg_xxxxzz_yyyzzz_0, tg_xxxxzz_yyyzzz_1, \
                                         tg_xxxxzz_yyzzz_1, tg_xxxxzz_yyzzzz_0, tg_xxxxzz_yyzzzz_1, tg_xxxxzz_yzzzz_1, \
                                         tg_xxxxzz_yzzzzz_0, tg_xxxxzz_yzzzzz_1, tg_xxxxzz_zzzzz_1, tg_xxxxzz_zzzzzz_0, \
                                         tg_xxxxzz_zzzzzz_1, tg_xxxyy_xxxyzz_0, tg_xxxyy_xxxyzz_1, tg_xxxyy_xxxzzz_0, \
                                         tg_xxxyy_xxxzzz_1, tg_xxxyy_xxyyyy_0, tg_xxxyy_xxyyyy_1, tg_xxxyy_xxyyyz_0, \
                                         tg_xxxyy_xxyyyz_1, tg_xxxyy_xxyyzz_0, tg_xxxyy_xxyyzz_1, tg_xxxyy_xxyzzz_0, \
                                         tg_xxxyy_xxyzzz_1, tg_xxxyy_xxzzzz_0, tg_xxxyy_xxzzzz_1, tg_xxxyy_xyyyyy_0, \
                                         tg_xxxyy_xyyyyy_1, tg_xxxyy_xyyyyz_0, tg_xxxyy_xyyyyz_1, tg_xxxyy_xyyyzz_0, \
                                         tg_xxxyy_xyyyzz_1, tg_xxxyy_xyyzzz_0, tg_xxxyy_xyyzzz_1, tg_xxxyy_xyzzzz_0, \
                                         tg_xxxyy_xyzzzz_1, tg_xxxyy_xzzzzz_0, tg_xxxyy_xzzzzz_1, tg_xxxyy_yyyyyy_0, \
                                         tg_xxxyy_yyyyyy_1, tg_xxxyy_yyyyyz_0, tg_xxxyy_yyyyyz_1, tg_xxxyy_yyyyzz_0, \
                                         tg_xxxyy_yyyyzz_1, tg_xxxyy_yyyzzz_0, tg_xxxyy_yyyzzz_1, tg_xxxyy_yyzzzz_0, \
                                         tg_xxxyy_yyzzzz_1, tg_xxxyy_yzzzzz_0, tg_xxxyy_yzzzzz_1, tg_xxxyy_zzzzzz_0, \
                                         tg_xxxyy_zzzzzz_1, tg_xxxyyy_xxxxx_1, tg_xxxyyy_xxxxxx_0, tg_xxxyyy_xxxxxx_1, \
                                         tg_xxxyyy_xxxxxy_0, tg_xxxyyy_xxxxxy_1, tg_xxxyyy_xxxxxz_0, tg_xxxyyy_xxxxxz_1, \
                                         tg_xxxyyy_xxxxy_1, tg_xxxyyy_xxxxyy_0, tg_xxxyyy_xxxxyy_1, tg_xxxyyy_xxxxyz_0, \
                                         tg_xxxyyy_xxxxyz_1, tg_xxxyyy_xxxxz_1, tg_xxxyyy_xxxxzz_0, tg_xxxyyy_xxxxzz_1, \
                                         tg_xxxyyy_xxxyy_1, tg_xxxyyy_xxxyyy_0, tg_xxxyyy_xxxyyy_1, tg_xxxyyy_xxxyyz_0, \
                                         tg_xxxyyy_xxxyyz_1, tg_xxxyyy_xxxyz_1, tg_xxxyyy_xxxyzz_0, tg_xxxyyy_xxxyzz_1, \
                                         tg_xxxyyy_xxxzz_1, tg_xxxyyy_xxxzzz_0, tg_xxxyyy_xxxzzz_1, tg_xxxyyy_xxyyy_1, \
                                         tg_xxxyyy_xxyyyy_0, tg_xxxyyy_xxyyyy_1, tg_xxxyyy_xxyyyz_0, tg_xxxyyy_xxyyyz_1, \
                                         tg_xxxyyy_xxyyz_1, tg_xxxyyy_xxyyzz_0, tg_xxxyyy_xxyyzz_1, tg_xxxyyy_xxyzz_1, \
                                         tg_xxxyyy_xxyzzz_0, tg_xxxyyy_xxyzzz_1, tg_xxxyyy_xxzzz_1, tg_xxxyyy_xxzzzz_0, \
                                         tg_xxxyyy_xxzzzz_1, tg_xxxyyy_xyyyy_1, tg_xxxyyy_xyyyyy_0, tg_xxxyyy_xyyyyy_1, \
                                         tg_xxxyyy_xyyyz_1, tg_xxxyyy_xyyzz_1, tg_xxxyyy_xyzzz_1, tg_xxxyyy_xzzzz_1, \
                                         tg_xxxyyy_yyyyy_1, tg_xxxyz_xxxxxx_0, tg_xxxyz_xxxxxx_1, tg_xxxyz_xxxxxy_0, \
                                         tg_xxxyz_xxxxxy_1, tg_xxxyz_xxxxxz_0, tg_xxxyz_xxxxxz_1, tg_xxxyz_xxxxyy_0, \
                                         tg_xxxyz_xxxxyy_1, tg_xxxyz_xxxxyz_0, tg_xxxyz_xxxxyz_1, tg_xxxyz_xxxxzz_0, \
                                         tg_xxxyz_xxxxzz_1, tg_xxxyz_xxxyyy_0, tg_xxxyz_xxxyyy_1, tg_xxxyz_xxxyyz_0, \
                                         tg_xxxyz_xxxyyz_1, tg_xxxyz_xxxyzz_0, tg_xxxyz_xxxyzz_1, tg_xxxyz_xxxzzz_0, \
                                         tg_xxxyz_xxxzzz_1, tg_xxxyz_xxyyyy_0, tg_xxxyz_xxyyyy_1, tg_xxxyz_xxyyyz_0, \
                                         tg_xxxyz_xxyyyz_1, tg_xxxyz_xxyyzz_0, tg_xxxyz_xxyyzz_1, tg_xxxyz_xxyzzz_0, \
                                         tg_xxxyz_xxyzzz_1, tg_xxxyz_xxzzzz_0, tg_xxxyz_xxzzzz_1, tg_xxxyz_xyyyyy_0, \
                                         tg_xxxyz_xyyyyy_1, tg_xxxyz_xyyyyz_0, tg_xxxyz_xyyyyz_1, tg_xxxyz_xyyyzz_0, \
                                         tg_xxxyz_xyyyzz_1, tg_xxxyz_xyyzzz_0, tg_xxxyz_xyyzzz_1, tg_xxxyz_xyzzzz_0, \
                                         tg_xxxyz_xyzzzz_1, tg_xxxyz_xzzzzz_0, tg_xxxyz_xzzzzz_1, tg_xxxyz_yyyyyy_0, \
                                         tg_xxxyz_yyyyyy_1, tg_xxxyz_yyyyyz_0, tg_xxxyz_yyyyyz_1, tg_xxxyz_yyyyzz_0, \
                                         tg_xxxyz_yyyyzz_1, tg_xxxyz_yyyzzz_0, tg_xxxyz_yyyzzz_1, tg_xxxyz_yyzzzz_0, \
                                         tg_xxxyz_yyzzzz_1, tg_xxxyz_yzzzzz_0, tg_xxxyz_yzzzzz_1, tg_xxxyz_zzzzzz_0, \
                                         tg_xxxyz_zzzzzz_1, tg_xxxzz_xxxxxx_0, tg_xxxzz_xxxxxx_1, tg_xxxzz_xxxxxy_0, \
                                         tg_xxxzz_xxxxxy_1, tg_xxxzz_xxxxxz_0, tg_xxxzz_xxxxxz_1, tg_xxxzz_xxxxyy_0, \
                                         tg_xxxzz_xxxxyy_1, tg_xxxzz_xxxxyz_0, tg_xxxzz_xxxxyz_1, tg_xxxzz_xxxxzz_0, \
                                         tg_xxxzz_xxxxzz_1, tg_xxxzz_xxxyyy_0, tg_xxxzz_xxxyyy_1, tg_xxxzz_xxxyyz_0, \
                                         tg_xxxzz_xxxyyz_1, tg_xxxzz_xxxyzz_0, tg_xxxzz_xxxyzz_1, tg_xxxzz_xxxzzz_0, \
                                         tg_xxxzz_xxxzzz_1, tg_xxxzz_xxyyyy_0, tg_xxxzz_xxyyyy_1, tg_xxxzz_xxyyyz_0, \
                                         tg_xxxzz_xxyyyz_1, tg_xxxzz_xxyyzz_0, tg_xxxzz_xxyyzz_1, tg_xxxzz_xxyzzz_0, \
                                         tg_xxxzz_xxyzzz_1, tg_xxxzz_xxzzzz_0, tg_xxxzz_xxzzzz_1, tg_xxxzz_xyyyyy_0, \
                                         tg_xxxzz_xyyyyy_1, tg_xxxzz_xyyyyz_0, tg_xxxzz_xyyyyz_1, tg_xxxzz_xyyyzz_0, \
                                         tg_xxxzz_xyyyzz_1, tg_xxxzz_xyyzzz_0, tg_xxxzz_xyyzzz_1, tg_xxxzz_xyzzzz_0, \
                                         tg_xxxzz_xyzzzz_1, tg_xxxzz_xzzzzz_0, tg_xxxzz_xzzzzz_1, tg_xxxzz_yyyyyy_0, \
                                         tg_xxxzz_yyyyyy_1, tg_xxxzz_yyyyyz_0, tg_xxxzz_yyyyyz_1, tg_xxxzz_yyyyzz_0, \
                                         tg_xxxzz_yyyyzz_1, tg_xxxzz_yyyzzz_0, tg_xxxzz_yyyzzz_1, tg_xxxzz_yyzzzz_0, \
                                         tg_xxxzz_yyzzzz_1, tg_xxxzz_yzzzzz_0, tg_xxxzz_yzzzzz_1, tg_xxxzz_zzzzzz_0, \
                                         tg_xxxzz_zzzzzz_1, tg_xxyyy_xxxxxx_0, tg_xxyyy_xxxxxx_1, tg_xxyyy_xxxxxy_0, \
                                         tg_xxyyy_xxxxxy_1, tg_xxyyy_xxxxxz_0, tg_xxyyy_xxxxxz_1, tg_xxyyy_xxxxyy_0, \
                                         tg_xxyyy_xxxxyy_1, tg_xxyyy_xxxxyz_0, tg_xxyyy_xxxxyz_1, tg_xxyyy_xxxxzz_0, \
                                         tg_xxyyy_xxxxzz_1, tg_xxyyy_xxxyyy_0, tg_xxyyy_xxxyyy_1, tg_xxyyy_xxxyyz_0, \
                                         tg_xxyyy_xxxyyz_1, tg_xxyyy_xxxyzz_0, tg_xxyyy_xxxyzz_1, tg_xxyyy_xxxzzz_0, \
                                         tg_xxyyy_xxxzzz_1, tg_xxyyy_xxyyyy_0, tg_xxyyy_xxyyyy_1, tg_xxyyy_xxyyyz_0, \
                                         tg_xxyyy_xxyyyz_1, tg_xxyyy_xxyyzz_0, tg_xxyyy_xxyyzz_1, tg_xxyyy_xxyzzz_0, \
                                         tg_xxyyy_xxyzzz_1, tg_xxyyy_xxzzzz_0, tg_xxyyy_xxzzzz_1, tg_xxyyy_xyyyyy_0, \
                                         tg_xxyyy_xyyyyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxyy_xxxyzz_0[j] = pb_x * tg_xxxxyy_xxxyzz_0[j] + fr * tg_xxxxyy_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxyzz_0[j] - tg_xxxyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxyzz_1[j];

                    tg_xxxxxyy_xxxzzz_0[j] = pb_x * tg_xxxxyy_xxxzzz_0[j] + fr * tg_xxxxyy_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxzzz_0[j] - tg_xxxyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxzzz_1[j];

                    tg_xxxxxyy_xxyyyy_0[j] = pb_x * tg_xxxxyy_xxyyyy_0[j] + fr * tg_xxxxyy_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyyyy_0[j] - tg_xxxyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyyyy_1[j];

                    tg_xxxxxyy_xxyyyz_0[j] = pb_x * tg_xxxxyy_xxyyyz_0[j] + fr * tg_xxxxyy_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyyyz_0[j] - tg_xxxyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyyyz_1[j];

                    tg_xxxxxyy_xxyyzz_0[j] = pb_x * tg_xxxxyy_xxyyzz_0[j] + fr * tg_xxxxyy_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyyzz_0[j] - tg_xxxyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyyzz_1[j];

                    tg_xxxxxyy_xxyzzz_0[j] = pb_x * tg_xxxxyy_xxyzzz_0[j] + fr * tg_xxxxyy_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyzzz_0[j] - tg_xxxyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyzzz_1[j];

                    tg_xxxxxyy_xxzzzz_0[j] = pb_x * tg_xxxxyy_xxzzzz_0[j] + fr * tg_xxxxyy_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxzzzz_0[j] - tg_xxxyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xzzzz_1[j];

                    tg_xxxxxyy_xyyyyy_0[j] = pb_x * tg_xxxxyy_xyyyyy_0[j] + fr * tg_xxxxyy_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyyyy_0[j] - tg_xxxyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyyyy_1[j];

                    tg_xxxxxyy_xyyyyz_0[j] = pb_x * tg_xxxxyy_xyyyyz_0[j] + fr * tg_xxxxyy_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyyyz_0[j] - tg_xxxyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyyyz_1[j];

                    tg_xxxxxyy_xyyyzz_0[j] = pb_x * tg_xxxxyy_xyyyzz_0[j] + fr * tg_xxxxyy_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyyzz_0[j] - tg_xxxyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyyzz_1[j];

                    tg_xxxxxyy_xyyzzz_0[j] = pb_x * tg_xxxxyy_xyyzzz_0[j] + fr * tg_xxxxyy_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyzzz_0[j] - tg_xxxyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyzzz_1[j];

                    tg_xxxxxyy_xyzzzz_0[j] = pb_x * tg_xxxxyy_xyzzzz_0[j] + fr * tg_xxxxyy_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyzzzz_0[j] - tg_xxxyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yzzzz_1[j];

                    tg_xxxxxyy_xzzzzz_0[j] = pb_x * tg_xxxxyy_xzzzzz_0[j] + fr * tg_xxxxyy_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xzzzzz_0[j] - tg_xxxyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_zzzzz_1[j];

                    tg_xxxxxyy_yyyyyy_0[j] = pb_x * tg_xxxxyy_yyyyyy_0[j] + fr * tg_xxxxyy_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyyyy_0[j] - tg_xxxyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxyy_yyyyyz_0[j] = pb_x * tg_xxxxyy_yyyyyz_0[j] + fr * tg_xxxxyy_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyyyz_0[j] - tg_xxxyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxyy_yyyyzz_0[j] = pb_x * tg_xxxxyy_yyyyzz_0[j] + fr * tg_xxxxyy_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyyzz_0[j] - tg_xxxyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxyy_yyyzzz_0[j] = pb_x * tg_xxxxyy_yyyzzz_0[j] + fr * tg_xxxxyy_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyzzz_0[j] - tg_xxxyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_yyzzzz_0[j] = pb_x * tg_xxxxyy_yyzzzz_0[j] + fr * tg_xxxxyy_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyzzzz_0[j] - tg_xxxyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_yzzzzz_0[j] = pb_x * tg_xxxxyy_yzzzzz_0[j] + fr * tg_xxxxyy_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yzzzzz_0[j] - tg_xxxyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_zzzzzz_0[j] = pb_x * tg_xxxxyy_zzzzzz_0[j] + fr * tg_xxxxyy_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_zzzzzz_0[j] - tg_xxxyy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_xxxxxx_0[j] = pb_x * tg_xxxxyz_xxxxxx_0[j] + fr * tg_xxxxyz_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxxx_0[j] - tg_xxxyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxyz_xxxxx_1[j];

                    tg_xxxxxyz_xxxxxy_0[j] = pb_x * tg_xxxxyz_xxxxxy_0[j] + fr * tg_xxxxyz_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxxy_0[j] - tg_xxxyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyz_xxxxy_1[j];

                    tg_xxxxxyz_xxxxxz_0[j] = pb_x * tg_xxxxyz_xxxxxz_0[j] + fr * tg_xxxxyz_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxxz_0[j] - tg_xxxyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxyz_xxxxz_1[j];

                    tg_xxxxxyz_xxxxyy_0[j] = pb_x * tg_xxxxyz_xxxxyy_0[j] + fr * tg_xxxxyz_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxyy_0[j] - tg_xxxyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyz_xxxyy_1[j];

                    tg_xxxxxyz_xxxxyz_0[j] = pb_x * tg_xxxxyz_xxxxyz_0[j] + fr * tg_xxxxyz_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxyz_0[j] - tg_xxxyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyz_xxxyz_1[j];

                    tg_xxxxxyz_xxxxzz_0[j] = pb_x * tg_xxxxyz_xxxxzz_0[j] + fr * tg_xxxxyz_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxxzz_0[j] - tg_xxxyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyz_xxxzz_1[j];

                    tg_xxxxxyz_xxxyyy_0[j] = pb_x * tg_xxxxyz_xxxyyy_0[j] + fr * tg_xxxxyz_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxyyy_0[j] - tg_xxxyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxyyy_1[j];

                    tg_xxxxxyz_xxxyyz_0[j] = pb_x * tg_xxxxyz_xxxyyz_0[j] + fr * tg_xxxxyz_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxyyz_0[j] - tg_xxxyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxyyz_1[j];

                    tg_xxxxxyz_xxxyzz_0[j] = pb_x * tg_xxxxyz_xxxyzz_0[j] + fr * tg_xxxxyz_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxyzz_0[j] - tg_xxxyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxyzz_1[j];

                    tg_xxxxxyz_xxxzzz_0[j] = pb_x * tg_xxxxyz_xxxzzz_0[j] + fr * tg_xxxxyz_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxzzz_0[j] - tg_xxxyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxzzz_1[j];

                    tg_xxxxxyz_xxyyyy_0[j] = pb_x * tg_xxxxyz_xxyyyy_0[j] + fr * tg_xxxxyz_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyyyy_0[j] - tg_xxxyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyyyy_1[j];

                    tg_xxxxxyz_xxyyyz_0[j] = pb_x * tg_xxxxyz_xxyyyz_0[j] + fr * tg_xxxxyz_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyyyz_0[j] - tg_xxxyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyyyz_1[j];

                    tg_xxxxxyz_xxyyzz_0[j] = pb_x * tg_xxxxyz_xxyyzz_0[j] + fr * tg_xxxxyz_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyyzz_0[j] - tg_xxxyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyyzz_1[j];

                    tg_xxxxxyz_xxyzzz_0[j] = pb_x * tg_xxxxyz_xxyzzz_0[j] + fr * tg_xxxxyz_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyzzz_0[j] - tg_xxxyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyzzz_1[j];

                    tg_xxxxxyz_xxzzzz_0[j] = pb_x * tg_xxxxyz_xxzzzz_0[j] + fr * tg_xxxxyz_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxzzzz_0[j] - tg_xxxyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xzzzz_1[j];

                    tg_xxxxxyz_xyyyyy_0[j] = pb_x * tg_xxxxyz_xyyyyy_0[j] + fr * tg_xxxxyz_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyyyy_0[j] - tg_xxxyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyyyy_1[j];

                    tg_xxxxxyz_xyyyyz_0[j] = pb_x * tg_xxxxyz_xyyyyz_0[j] + fr * tg_xxxxyz_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyyyz_0[j] - tg_xxxyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyyyz_1[j];

                    tg_xxxxxyz_xyyyzz_0[j] = pb_x * tg_xxxxyz_xyyyzz_0[j] + fr * tg_xxxxyz_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyyzz_0[j] - tg_xxxyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyyzz_1[j];

                    tg_xxxxxyz_xyyzzz_0[j] = pb_x * tg_xxxxyz_xyyzzz_0[j] + fr * tg_xxxxyz_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyzzz_0[j] - tg_xxxyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyzzz_1[j];

                    tg_xxxxxyz_xyzzzz_0[j] = pb_x * tg_xxxxyz_xyzzzz_0[j] + fr * tg_xxxxyz_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyzzzz_0[j] - tg_xxxyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yzzzz_1[j];

                    tg_xxxxxyz_xzzzzz_0[j] = pb_x * tg_xxxxyz_xzzzzz_0[j] + fr * tg_xxxxyz_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xzzzzz_0[j] - tg_xxxyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_zzzzz_1[j];

                    tg_xxxxxyz_yyyyyy_0[j] = pb_x * tg_xxxxyz_yyyyyy_0[j] + fr * tg_xxxxyz_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyyyy_0[j] - tg_xxxyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxyz_yyyyyz_0[j] = pb_x * tg_xxxxyz_yyyyyz_0[j] + fr * tg_xxxxyz_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyyyz_0[j] - tg_xxxyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxyz_yyyyzz_0[j] = pb_x * tg_xxxxyz_yyyyzz_0[j] + fr * tg_xxxxyz_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyyzz_0[j] - tg_xxxyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxyz_yyyzzz_0[j] = pb_x * tg_xxxxyz_yyyzzz_0[j] + fr * tg_xxxxyz_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyzzz_0[j] - tg_xxxyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_yyzzzz_0[j] = pb_x * tg_xxxxyz_yyzzzz_0[j] + fr * tg_xxxxyz_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyzzzz_0[j] - tg_xxxyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_yzzzzz_0[j] = pb_x * tg_xxxxyz_yzzzzz_0[j] + fr * tg_xxxxyz_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yzzzzz_0[j] - tg_xxxyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_zzzzzz_0[j] = pb_x * tg_xxxxyz_zzzzzz_0[j] + fr * tg_xxxxyz_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_zzzzzz_0[j] - tg_xxxyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_xxxxxx_0[j] = pb_x * tg_xxxxzz_xxxxxx_0[j] + fr * tg_xxxxzz_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxxx_0[j] - tg_xxxzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxxzz_xxxxx_1[j];

                    tg_xxxxxzz_xxxxxy_0[j] = pb_x * tg_xxxxzz_xxxxxy_0[j] + fr * tg_xxxxzz_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxxy_0[j] - tg_xxxzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxzz_xxxxy_1[j];

                    tg_xxxxxzz_xxxxxz_0[j] = pb_x * tg_xxxxzz_xxxxxz_0[j] + fr * tg_xxxxzz_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxxz_0[j] - tg_xxxzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxxzz_xxxxz_1[j];

                    tg_xxxxxzz_xxxxyy_0[j] = pb_x * tg_xxxxzz_xxxxyy_0[j] + fr * tg_xxxxzz_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxyy_0[j] - tg_xxxzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzz_xxxyy_1[j];

                    tg_xxxxxzz_xxxxyz_0[j] = pb_x * tg_xxxxzz_xxxxyz_0[j] + fr * tg_xxxxzz_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxyz_0[j] - tg_xxxzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzz_xxxyz_1[j];

                    tg_xxxxxzz_xxxxzz_0[j] = pb_x * tg_xxxxzz_xxxxzz_0[j] + fr * tg_xxxxzz_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxxzz_0[j] - tg_xxxzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzz_xxxzz_1[j];

                    tg_xxxxxzz_xxxyyy_0[j] = pb_x * tg_xxxxzz_xxxyyy_0[j] + fr * tg_xxxxzz_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxyyy_0[j] - tg_xxxzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxyyy_1[j];

                    tg_xxxxxzz_xxxyyz_0[j] = pb_x * tg_xxxxzz_xxxyyz_0[j] + fr * tg_xxxxzz_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxyyz_0[j] - tg_xxxzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxyyz_1[j];

                    tg_xxxxxzz_xxxyzz_0[j] = pb_x * tg_xxxxzz_xxxyzz_0[j] + fr * tg_xxxxzz_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxyzz_0[j] - tg_xxxzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxyzz_1[j];

                    tg_xxxxxzz_xxxzzz_0[j] = pb_x * tg_xxxxzz_xxxzzz_0[j] + fr * tg_xxxxzz_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxzzz_0[j] - tg_xxxzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxzzz_1[j];

                    tg_xxxxxzz_xxyyyy_0[j] = pb_x * tg_xxxxzz_xxyyyy_0[j] + fr * tg_xxxxzz_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyyyy_0[j] - tg_xxxzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyyyy_1[j];

                    tg_xxxxxzz_xxyyyz_0[j] = pb_x * tg_xxxxzz_xxyyyz_0[j] + fr * tg_xxxxzz_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyyyz_0[j] - tg_xxxzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyyyz_1[j];

                    tg_xxxxxzz_xxyyzz_0[j] = pb_x * tg_xxxxzz_xxyyzz_0[j] + fr * tg_xxxxzz_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyyzz_0[j] - tg_xxxzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyyzz_1[j];

                    tg_xxxxxzz_xxyzzz_0[j] = pb_x * tg_xxxxzz_xxyzzz_0[j] + fr * tg_xxxxzz_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyzzz_0[j] - tg_xxxzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyzzz_1[j];

                    tg_xxxxxzz_xxzzzz_0[j] = pb_x * tg_xxxxzz_xxzzzz_0[j] + fr * tg_xxxxzz_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxzzzz_0[j] - tg_xxxzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xzzzz_1[j];

                    tg_xxxxxzz_xyyyyy_0[j] = pb_x * tg_xxxxzz_xyyyyy_0[j] + fr * tg_xxxxzz_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyyyy_0[j] - tg_xxxzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyyyy_1[j];

                    tg_xxxxxzz_xyyyyz_0[j] = pb_x * tg_xxxxzz_xyyyyz_0[j] + fr * tg_xxxxzz_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyyyz_0[j] - tg_xxxzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyyyz_1[j];

                    tg_xxxxxzz_xyyyzz_0[j] = pb_x * tg_xxxxzz_xyyyzz_0[j] + fr * tg_xxxxzz_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyyzz_0[j] - tg_xxxzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyyzz_1[j];

                    tg_xxxxxzz_xyyzzz_0[j] = pb_x * tg_xxxxzz_xyyzzz_0[j] + fr * tg_xxxxzz_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyzzz_0[j] - tg_xxxzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyzzz_1[j];

                    tg_xxxxxzz_xyzzzz_0[j] = pb_x * tg_xxxxzz_xyzzzz_0[j] + fr * tg_xxxxzz_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyzzzz_0[j] - tg_xxxzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yzzzz_1[j];

                    tg_xxxxxzz_xzzzzz_0[j] = pb_x * tg_xxxxzz_xzzzzz_0[j] + fr * tg_xxxxzz_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xzzzzz_0[j] - tg_xxxzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_zzzzz_1[j];

                    tg_xxxxxzz_yyyyyy_0[j] = pb_x * tg_xxxxzz_yyyyyy_0[j] + fr * tg_xxxxzz_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyyyy_0[j] - tg_xxxzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxxzz_yyyyyz_0[j] = pb_x * tg_xxxxzz_yyyyyz_0[j] + fr * tg_xxxxzz_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyyyz_0[j] - tg_xxxzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxxzz_yyyyzz_0[j] = pb_x * tg_xxxxzz_yyyyzz_0[j] + fr * tg_xxxxzz_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyyzz_0[j] - tg_xxxzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxxzz_yyyzzz_0[j] = pb_x * tg_xxxxzz_yyyzzz_0[j] + fr * tg_xxxxzz_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyzzz_0[j] - tg_xxxzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_yyzzzz_0[j] = pb_x * tg_xxxxzz_yyzzzz_0[j] + fr * tg_xxxxzz_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyzzzz_0[j] - tg_xxxzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_yzzzzz_0[j] = pb_x * tg_xxxxzz_yzzzzz_0[j] + fr * tg_xxxxzz_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yzzzzz_0[j] - tg_xxxzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_zzzzzz_0[j] = pb_x * tg_xxxxzz_zzzzzz_0[j] + fr * tg_xxxxzz_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_zzzzzz_0[j] - tg_xxxzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_xxxxxx_0[j] = pb_x * tg_xxxyyy_xxxxxx_0[j] + fr * tg_xxxyyy_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxxx_0[j] - tg_xxyyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxyyy_xxxxx_1[j];

                    tg_xxxxyyy_xxxxxy_0[j] = pb_x * tg_xxxyyy_xxxxxy_0[j] + fr * tg_xxxyyy_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxxy_0[j] - tg_xxyyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyy_xxxxy_1[j];

                    tg_xxxxyyy_xxxxxz_0[j] = pb_x * tg_xxxyyy_xxxxxz_0[j] + fr * tg_xxxyyy_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxxz_0[j] - tg_xxyyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyy_xxxxz_1[j];

                    tg_xxxxyyy_xxxxyy_0[j] = pb_x * tg_xxxyyy_xxxxyy_0[j] + fr * tg_xxxyyy_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxyy_0[j] - tg_xxyyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyy_xxxyy_1[j];

                    tg_xxxxyyy_xxxxyz_0[j] = pb_x * tg_xxxyyy_xxxxyz_0[j] + fr * tg_xxxyyy_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxyz_0[j] - tg_xxyyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyy_xxxyz_1[j];

                    tg_xxxxyyy_xxxxzz_0[j] = pb_x * tg_xxxyyy_xxxxzz_0[j] + fr * tg_xxxyyy_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxxzz_0[j] - tg_xxyyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyy_xxxzz_1[j];

                    tg_xxxxyyy_xxxyyy_0[j] = pb_x * tg_xxxyyy_xxxyyy_0[j] + fr * tg_xxxyyy_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxyyy_0[j] - tg_xxyyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxyyy_1[j];

                    tg_xxxxyyy_xxxyyz_0[j] = pb_x * tg_xxxyyy_xxxyyz_0[j] + fr * tg_xxxyyy_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxyyz_0[j] - tg_xxyyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxyyz_1[j];

                    tg_xxxxyyy_xxxyzz_0[j] = pb_x * tg_xxxyyy_xxxyzz_0[j] + fr * tg_xxxyyy_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxyzz_0[j] - tg_xxyyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxyzz_1[j];

                    tg_xxxxyyy_xxxzzz_0[j] = pb_x * tg_xxxyyy_xxxzzz_0[j] + fr * tg_xxxyyy_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxzzz_0[j] - tg_xxyyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxzzz_1[j];

                    tg_xxxxyyy_xxyyyy_0[j] = pb_x * tg_xxxyyy_xxyyyy_0[j] + fr * tg_xxxyyy_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyyyy_0[j] - tg_xxyyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyyyy_1[j];

                    tg_xxxxyyy_xxyyyz_0[j] = pb_x * tg_xxxyyy_xxyyyz_0[j] + fr * tg_xxxyyy_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyyyz_0[j] - tg_xxyyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyyyz_1[j];

                    tg_xxxxyyy_xxyyzz_0[j] = pb_x * tg_xxxyyy_xxyyzz_0[j] + fr * tg_xxxyyy_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyyzz_0[j] - tg_xxyyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyyzz_1[j];

                    tg_xxxxyyy_xxyzzz_0[j] = pb_x * tg_xxxyyy_xxyzzz_0[j] + fr * tg_xxxyyy_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyzzz_0[j] - tg_xxyyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyzzz_1[j];

                    tg_xxxxyyy_xxzzzz_0[j] = pb_x * tg_xxxyyy_xxzzzz_0[j] + fr * tg_xxxyyy_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxzzzz_0[j] - tg_xxyyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xzzzz_1[j];

                    tg_xxxxyyy_xyyyyy_0[j] = pb_x * tg_xxxyyy_xyyyyy_0[j] + fr * tg_xxxyyy_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyyyy_0[j] - tg_xxyyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_184_276(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (184,276)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxxyyy_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 184); 

                auto tg_xxxyyy_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 185); 

                auto tg_xxxyyy_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 186); 

                auto tg_xxxyyy_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 187); 

                auto tg_xxxyyy_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 188); 

                auto tg_xxxyyy_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 189); 

                auto tg_xxxyyy_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 190); 

                auto tg_xxxyyy_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 191); 

                auto tg_xxxyyy_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 192); 

                auto tg_xxxyyy_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 193); 

                auto tg_xxxyyy_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 194); 

                auto tg_xxxyyy_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 195); 

                auto tg_xxxyyz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 196); 

                auto tg_xxxyyz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 197); 

                auto tg_xxxyyz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 198); 

                auto tg_xxxyyz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 199); 

                auto tg_xxxyyz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 200); 

                auto tg_xxxyyz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 201); 

                auto tg_xxxyyz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 202); 

                auto tg_xxxyyz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 203); 

                auto tg_xxxyyz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 204); 

                auto tg_xxxyyz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 205); 

                auto tg_xxxyyz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 206); 

                auto tg_xxxyyz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 207); 

                auto tg_xxxyyz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 208); 

                auto tg_xxxyyz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 209); 

                auto tg_xxxyyz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 210); 

                auto tg_xxxyyz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 211); 

                auto tg_xxxyyz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 212); 

                auto tg_xxxyyz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 213); 

                auto tg_xxxyyz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 214); 

                auto tg_xxxyyz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 215); 

                auto tg_xxxyyz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 216); 

                auto tg_xxxyyz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 217); 

                auto tg_xxxyyz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 218); 

                auto tg_xxxyyz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 219); 

                auto tg_xxxyyz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 220); 

                auto tg_xxxyyz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 221); 

                auto tg_xxxyyz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 222); 

                auto tg_xxxyyz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 223); 

                auto tg_xxxyzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 224); 

                auto tg_xxxyzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 225); 

                auto tg_xxxyzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 226); 

                auto tg_xxxyzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 227); 

                auto tg_xxxyzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 228); 

                auto tg_xxxyzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 229); 

                auto tg_xxxyzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 230); 

                auto tg_xxxyzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 231); 

                auto tg_xxxyzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 232); 

                auto tg_xxxyzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 233); 

                auto tg_xxxyzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 234); 

                auto tg_xxxyzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 235); 

                auto tg_xxxyzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 236); 

                auto tg_xxxyzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 237); 

                auto tg_xxxyzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 238); 

                auto tg_xxxyzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 239); 

                auto tg_xxxyzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 240); 

                auto tg_xxxyzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 241); 

                auto tg_xxxyzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 242); 

                auto tg_xxxyzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 243); 

                auto tg_xxxyzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 244); 

                auto tg_xxxyzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 245); 

                auto tg_xxxyzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 246); 

                auto tg_xxxyzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 247); 

                auto tg_xxxyzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 248); 

                auto tg_xxxyzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 249); 

                auto tg_xxxyzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 250); 

                auto tg_xxxyzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 251); 

                auto tg_xxxzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 252); 

                auto tg_xxxzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 253); 

                auto tg_xxxzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 254); 

                auto tg_xxxzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 255); 

                auto tg_xxxzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 256); 

                auto tg_xxxzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 257); 

                auto tg_xxxzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 258); 

                auto tg_xxxzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 259); 

                auto tg_xxxzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 260); 

                auto tg_xxxzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 261); 

                auto tg_xxxzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 262); 

                auto tg_xxxzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 263); 

                auto tg_xxxzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 264); 

                auto tg_xxxzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 265); 

                auto tg_xxxzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 266); 

                auto tg_xxxzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 267); 

                auto tg_xxxzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 268); 

                auto tg_xxxzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 269); 

                auto tg_xxxzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 270); 

                auto tg_xxxzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 271); 

                auto tg_xxxzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 272); 

                auto tg_xxxzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 273); 

                auto tg_xxxzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 274); 

                auto tg_xxxzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 275); 

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

                // set up pointers to integrals

                auto tg_xxxxyyy_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 184); 

                auto tg_xxxxyyy_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 185); 

                auto tg_xxxxyyy_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 186); 

                auto tg_xxxxyyy_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 187); 

                auto tg_xxxxyyy_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 188); 

                auto tg_xxxxyyy_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 189); 

                auto tg_xxxxyyy_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 190); 

                auto tg_xxxxyyy_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 191); 

                auto tg_xxxxyyy_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 192); 

                auto tg_xxxxyyy_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 193); 

                auto tg_xxxxyyy_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 194); 

                auto tg_xxxxyyy_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 195); 

                auto tg_xxxxyyz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 196); 

                auto tg_xxxxyyz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 197); 

                auto tg_xxxxyyz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 198); 

                auto tg_xxxxyyz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 199); 

                auto tg_xxxxyyz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 200); 

                auto tg_xxxxyyz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 201); 

                auto tg_xxxxyyz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 202); 

                auto tg_xxxxyyz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 203); 

                auto tg_xxxxyyz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 204); 

                auto tg_xxxxyyz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 205); 

                auto tg_xxxxyyz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 206); 

                auto tg_xxxxyyz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 207); 

                auto tg_xxxxyyz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 208); 

                auto tg_xxxxyyz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 209); 

                auto tg_xxxxyyz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 210); 

                auto tg_xxxxyyz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 211); 

                auto tg_xxxxyyz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 212); 

                auto tg_xxxxyyz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 213); 

                auto tg_xxxxyyz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 214); 

                auto tg_xxxxyyz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 215); 

                auto tg_xxxxyyz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 216); 

                auto tg_xxxxyyz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 217); 

                auto tg_xxxxyyz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 218); 

                auto tg_xxxxyyz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 219); 

                auto tg_xxxxyyz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 220); 

                auto tg_xxxxyyz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 221); 

                auto tg_xxxxyyz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 222); 

                auto tg_xxxxyyz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 223); 

                auto tg_xxxxyzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 224); 

                auto tg_xxxxyzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 225); 

                auto tg_xxxxyzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 226); 

                auto tg_xxxxyzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 227); 

                auto tg_xxxxyzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 228); 

                auto tg_xxxxyzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 229); 

                auto tg_xxxxyzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 230); 

                auto tg_xxxxyzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 231); 

                auto tg_xxxxyzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 232); 

                auto tg_xxxxyzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 233); 

                auto tg_xxxxyzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 234); 

                auto tg_xxxxyzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 235); 

                auto tg_xxxxyzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 236); 

                auto tg_xxxxyzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 237); 

                auto tg_xxxxyzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 238); 

                auto tg_xxxxyzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 239); 

                auto tg_xxxxyzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 240); 

                auto tg_xxxxyzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 241); 

                auto tg_xxxxyzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 242); 

                auto tg_xxxxyzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 243); 

                auto tg_xxxxyzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 244); 

                auto tg_xxxxyzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 245); 

                auto tg_xxxxyzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 246); 

                auto tg_xxxxyzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 247); 

                auto tg_xxxxyzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 248); 

                auto tg_xxxxyzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 249); 

                auto tg_xxxxyzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 250); 

                auto tg_xxxxyzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 251); 

                auto tg_xxxxzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 252); 

                auto tg_xxxxzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 253); 

                auto tg_xxxxzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 254); 

                auto tg_xxxxzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 255); 

                auto tg_xxxxzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 256); 

                auto tg_xxxxzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 257); 

                auto tg_xxxxzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 258); 

                auto tg_xxxxzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 259); 

                auto tg_xxxxzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 260); 

                auto tg_xxxxzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 261); 

                auto tg_xxxxzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 262); 

                auto tg_xxxxzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 263); 

                auto tg_xxxxzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 264); 

                auto tg_xxxxzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 265); 

                auto tg_xxxxzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 266); 

                auto tg_xxxxzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 267); 

                auto tg_xxxxzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 268); 

                auto tg_xxxxzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 269); 

                auto tg_xxxxzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 270); 

                auto tg_xxxxzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 271); 

                auto tg_xxxxzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 272); 

                auto tg_xxxxzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 273); 

                auto tg_xxxxzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 274); 

                auto tg_xxxxzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 275); 

                // Batch of Integrals (184,276)

                #pragma omp simd aligned(fxn, fza, tg_xxxxyyy_xyyyyz_0, tg_xxxxyyy_xyyyzz_0, \
                                         tg_xxxxyyy_xyyzzz_0, tg_xxxxyyy_xyzzzz_0, tg_xxxxyyy_xzzzzz_0, tg_xxxxyyy_yyyyyy_0, \
                                         tg_xxxxyyy_yyyyyz_0, tg_xxxxyyy_yyyyzz_0, tg_xxxxyyy_yyyzzz_0, tg_xxxxyyy_yyzzzz_0, \
                                         tg_xxxxyyy_yzzzzz_0, tg_xxxxyyy_zzzzzz_0, tg_xxxxyyz_xxxxxx_0, tg_xxxxyyz_xxxxxy_0, \
                                         tg_xxxxyyz_xxxxxz_0, tg_xxxxyyz_xxxxyy_0, tg_xxxxyyz_xxxxyz_0, tg_xxxxyyz_xxxxzz_0, \
                                         tg_xxxxyyz_xxxyyy_0, tg_xxxxyyz_xxxyyz_0, tg_xxxxyyz_xxxyzz_0, tg_xxxxyyz_xxxzzz_0, \
                                         tg_xxxxyyz_xxyyyy_0, tg_xxxxyyz_xxyyyz_0, tg_xxxxyyz_xxyyzz_0, tg_xxxxyyz_xxyzzz_0, \
                                         tg_xxxxyyz_xxzzzz_0, tg_xxxxyyz_xyyyyy_0, tg_xxxxyyz_xyyyyz_0, tg_xxxxyyz_xyyyzz_0, \
                                         tg_xxxxyyz_xyyzzz_0, tg_xxxxyyz_xyzzzz_0, tg_xxxxyyz_xzzzzz_0, tg_xxxxyyz_yyyyyy_0, \
                                         tg_xxxxyyz_yyyyyz_0, tg_xxxxyyz_yyyyzz_0, tg_xxxxyyz_yyyzzz_0, tg_xxxxyyz_yyzzzz_0, \
                                         tg_xxxxyyz_yzzzzz_0, tg_xxxxyyz_zzzzzz_0, tg_xxxxyzz_xxxxxx_0, tg_xxxxyzz_xxxxxy_0, \
                                         tg_xxxxyzz_xxxxxz_0, tg_xxxxyzz_xxxxyy_0, tg_xxxxyzz_xxxxyz_0, tg_xxxxyzz_xxxxzz_0, \
                                         tg_xxxxyzz_xxxyyy_0, tg_xxxxyzz_xxxyyz_0, tg_xxxxyzz_xxxyzz_0, tg_xxxxyzz_xxxzzz_0, \
                                         tg_xxxxyzz_xxyyyy_0, tg_xxxxyzz_xxyyyz_0, tg_xxxxyzz_xxyyzz_0, tg_xxxxyzz_xxyzzz_0, \
                                         tg_xxxxyzz_xxzzzz_0, tg_xxxxyzz_xyyyyy_0, tg_xxxxyzz_xyyyyz_0, tg_xxxxyzz_xyyyzz_0, \
                                         tg_xxxxyzz_xyyzzz_0, tg_xxxxyzz_xyzzzz_0, tg_xxxxyzz_xzzzzz_0, tg_xxxxyzz_yyyyyy_0, \
                                         tg_xxxxyzz_yyyyyz_0, tg_xxxxyzz_yyyyzz_0, tg_xxxxyzz_yyyzzz_0, tg_xxxxyzz_yyzzzz_0, \
                                         tg_xxxxyzz_yzzzzz_0, tg_xxxxyzz_zzzzzz_0, tg_xxxxzzz_xxxxxx_0, tg_xxxxzzz_xxxxxy_0, \
                                         tg_xxxxzzz_xxxxxz_0, tg_xxxxzzz_xxxxyy_0, tg_xxxxzzz_xxxxyz_0, tg_xxxxzzz_xxxxzz_0, \
                                         tg_xxxxzzz_xxxyyy_0, tg_xxxxzzz_xxxyyz_0, tg_xxxxzzz_xxxyzz_0, tg_xxxxzzz_xxxzzz_0, \
                                         tg_xxxxzzz_xxyyyy_0, tg_xxxxzzz_xxyyyz_0, tg_xxxxzzz_xxyyzz_0, tg_xxxxzzz_xxyzzz_0, \
                                         tg_xxxxzzz_xxzzzz_0, tg_xxxxzzz_xyyyyy_0, tg_xxxxzzz_xyyyyz_0, tg_xxxxzzz_xyyyzz_0, \
                                         tg_xxxxzzz_xyyzzz_0, tg_xxxxzzz_xyzzzz_0, tg_xxxxzzz_xzzzzz_0, tg_xxxxzzz_yyyyyy_0, \
                                         tg_xxxxzzz_yyyyyz_0, tg_xxxxzzz_yyyyzz_0, tg_xxxyyy_xyyyyz_0, tg_xxxyyy_xyyyyz_1, \
                                         tg_xxxyyy_xyyyzz_0, tg_xxxyyy_xyyyzz_1, tg_xxxyyy_xyyzzz_0, tg_xxxyyy_xyyzzz_1, \
                                         tg_xxxyyy_xyzzzz_0, tg_xxxyyy_xyzzzz_1, tg_xxxyyy_xzzzzz_0, tg_xxxyyy_xzzzzz_1, \
                                         tg_xxxyyy_yyyyyy_0, tg_xxxyyy_yyyyyy_1, tg_xxxyyy_yyyyyz_0, tg_xxxyyy_yyyyyz_1, \
                                         tg_xxxyyy_yyyyz_1, tg_xxxyyy_yyyyzz_0, tg_xxxyyy_yyyyzz_1, tg_xxxyyy_yyyzz_1, \
                                         tg_xxxyyy_yyyzzz_0, tg_xxxyyy_yyyzzz_1, tg_xxxyyy_yyzzz_1, tg_xxxyyy_yyzzzz_0, \
                                         tg_xxxyyy_yyzzzz_1, tg_xxxyyy_yzzzz_1, tg_xxxyyy_yzzzzz_0, tg_xxxyyy_yzzzzz_1, \
                                         tg_xxxyyy_zzzzz_1, tg_xxxyyy_zzzzzz_0, tg_xxxyyy_zzzzzz_1, tg_xxxyyz_xxxxx_1, \
                                         tg_xxxyyz_xxxxxx_0, tg_xxxyyz_xxxxxx_1, tg_xxxyyz_xxxxxy_0, tg_xxxyyz_xxxxxy_1, \
                                         tg_xxxyyz_xxxxxz_0, tg_xxxyyz_xxxxxz_1, tg_xxxyyz_xxxxy_1, tg_xxxyyz_xxxxyy_0, \
                                         tg_xxxyyz_xxxxyy_1, tg_xxxyyz_xxxxyz_0, tg_xxxyyz_xxxxyz_1, tg_xxxyyz_xxxxz_1, \
                                         tg_xxxyyz_xxxxzz_0, tg_xxxyyz_xxxxzz_1, tg_xxxyyz_xxxyy_1, tg_xxxyyz_xxxyyy_0, \
                                         tg_xxxyyz_xxxyyy_1, tg_xxxyyz_xxxyyz_0, tg_xxxyyz_xxxyyz_1, tg_xxxyyz_xxxyz_1, \
                                         tg_xxxyyz_xxxyzz_0, tg_xxxyyz_xxxyzz_1, tg_xxxyyz_xxxzz_1, tg_xxxyyz_xxxzzz_0, \
                                         tg_xxxyyz_xxxzzz_1, tg_xxxyyz_xxyyy_1, tg_xxxyyz_xxyyyy_0, tg_xxxyyz_xxyyyy_1, \
                                         tg_xxxyyz_xxyyyz_0, tg_xxxyyz_xxyyyz_1, tg_xxxyyz_xxyyz_1, tg_xxxyyz_xxyyzz_0, \
                                         tg_xxxyyz_xxyyzz_1, tg_xxxyyz_xxyzz_1, tg_xxxyyz_xxyzzz_0, tg_xxxyyz_xxyzzz_1, \
                                         tg_xxxyyz_xxzzz_1, tg_xxxyyz_xxzzzz_0, tg_xxxyyz_xxzzzz_1, tg_xxxyyz_xyyyy_1, \
                                         tg_xxxyyz_xyyyyy_0, tg_xxxyyz_xyyyyy_1, tg_xxxyyz_xyyyyz_0, tg_xxxyyz_xyyyyz_1, \
                                         tg_xxxyyz_xyyyz_1, tg_xxxyyz_xyyyzz_0, tg_xxxyyz_xyyyzz_1, tg_xxxyyz_xyyzz_1, \
                                         tg_xxxyyz_xyyzzz_0, tg_xxxyyz_xyyzzz_1, tg_xxxyyz_xyzzz_1, tg_xxxyyz_xyzzzz_0, \
                                         tg_xxxyyz_xyzzzz_1, tg_xxxyyz_xzzzz_1, tg_xxxyyz_xzzzzz_0, tg_xxxyyz_xzzzzz_1, \
                                         tg_xxxyyz_yyyyy_1, tg_xxxyyz_yyyyyy_0, tg_xxxyyz_yyyyyy_1, tg_xxxyyz_yyyyyz_0, \
                                         tg_xxxyyz_yyyyyz_1, tg_xxxyyz_yyyyz_1, tg_xxxyyz_yyyyzz_0, tg_xxxyyz_yyyyzz_1, \
                                         tg_xxxyyz_yyyzz_1, tg_xxxyyz_yyyzzz_0, tg_xxxyyz_yyyzzz_1, tg_xxxyyz_yyzzz_1, \
                                         tg_xxxyyz_yyzzzz_0, tg_xxxyyz_yyzzzz_1, tg_xxxyyz_yzzzz_1, tg_xxxyyz_yzzzzz_0, \
                                         tg_xxxyyz_yzzzzz_1, tg_xxxyyz_zzzzz_1, tg_xxxyyz_zzzzzz_0, tg_xxxyyz_zzzzzz_1, \
                                         tg_xxxyzz_xxxxx_1, tg_xxxyzz_xxxxxx_0, tg_xxxyzz_xxxxxx_1, tg_xxxyzz_xxxxxy_0, \
                                         tg_xxxyzz_xxxxxy_1, tg_xxxyzz_xxxxxz_0, tg_xxxyzz_xxxxxz_1, tg_xxxyzz_xxxxy_1, \
                                         tg_xxxyzz_xxxxyy_0, tg_xxxyzz_xxxxyy_1, tg_xxxyzz_xxxxyz_0, tg_xxxyzz_xxxxyz_1, \
                                         tg_xxxyzz_xxxxz_1, tg_xxxyzz_xxxxzz_0, tg_xxxyzz_xxxxzz_1, tg_xxxyzz_xxxyy_1, \
                                         tg_xxxyzz_xxxyyy_0, tg_xxxyzz_xxxyyy_1, tg_xxxyzz_xxxyyz_0, tg_xxxyzz_xxxyyz_1, \
                                         tg_xxxyzz_xxxyz_1, tg_xxxyzz_xxxyzz_0, tg_xxxyzz_xxxyzz_1, tg_xxxyzz_xxxzz_1, \
                                         tg_xxxyzz_xxxzzz_0, tg_xxxyzz_xxxzzz_1, tg_xxxyzz_xxyyy_1, tg_xxxyzz_xxyyyy_0, \
                                         tg_xxxyzz_xxyyyy_1, tg_xxxyzz_xxyyyz_0, tg_xxxyzz_xxyyyz_1, tg_xxxyzz_xxyyz_1, \
                                         tg_xxxyzz_xxyyzz_0, tg_xxxyzz_xxyyzz_1, tg_xxxyzz_xxyzz_1, tg_xxxyzz_xxyzzz_0, \
                                         tg_xxxyzz_xxyzzz_1, tg_xxxyzz_xxzzz_1, tg_xxxyzz_xxzzzz_0, tg_xxxyzz_xxzzzz_1, \
                                         tg_xxxyzz_xyyyy_1, tg_xxxyzz_xyyyyy_0, tg_xxxyzz_xyyyyy_1, tg_xxxyzz_xyyyyz_0, \
                                         tg_xxxyzz_xyyyyz_1, tg_xxxyzz_xyyyz_1, tg_xxxyzz_xyyyzz_0, tg_xxxyzz_xyyyzz_1, \
                                         tg_xxxyzz_xyyzz_1, tg_xxxyzz_xyyzzz_0, tg_xxxyzz_xyyzzz_1, tg_xxxyzz_xyzzz_1, \
                                         tg_xxxyzz_xyzzzz_0, tg_xxxyzz_xyzzzz_1, tg_xxxyzz_xzzzz_1, tg_xxxyzz_xzzzzz_0, \
                                         tg_xxxyzz_xzzzzz_1, tg_xxxyzz_yyyyy_1, tg_xxxyzz_yyyyyy_0, tg_xxxyzz_yyyyyy_1, \
                                         tg_xxxyzz_yyyyyz_0, tg_xxxyzz_yyyyyz_1, tg_xxxyzz_yyyyz_1, tg_xxxyzz_yyyyzz_0, \
                                         tg_xxxyzz_yyyyzz_1, tg_xxxyzz_yyyzz_1, tg_xxxyzz_yyyzzz_0, tg_xxxyzz_yyyzzz_1, \
                                         tg_xxxyzz_yyzzz_1, tg_xxxyzz_yyzzzz_0, tg_xxxyzz_yyzzzz_1, tg_xxxyzz_yzzzz_1, \
                                         tg_xxxyzz_yzzzzz_0, tg_xxxyzz_yzzzzz_1, tg_xxxyzz_zzzzz_1, tg_xxxyzz_zzzzzz_0, \
                                         tg_xxxyzz_zzzzzz_1, tg_xxxzzz_xxxxx_1, tg_xxxzzz_xxxxxx_0, tg_xxxzzz_xxxxxx_1, \
                                         tg_xxxzzz_xxxxxy_0, tg_xxxzzz_xxxxxy_1, tg_xxxzzz_xxxxxz_0, tg_xxxzzz_xxxxxz_1, \
                                         tg_xxxzzz_xxxxy_1, tg_xxxzzz_xxxxyy_0, tg_xxxzzz_xxxxyy_1, tg_xxxzzz_xxxxyz_0, \
                                         tg_xxxzzz_xxxxyz_1, tg_xxxzzz_xxxxz_1, tg_xxxzzz_xxxxzz_0, tg_xxxzzz_xxxxzz_1, \
                                         tg_xxxzzz_xxxyy_1, tg_xxxzzz_xxxyyy_0, tg_xxxzzz_xxxyyy_1, tg_xxxzzz_xxxyyz_0, \
                                         tg_xxxzzz_xxxyyz_1, tg_xxxzzz_xxxyz_1, tg_xxxzzz_xxxyzz_0, tg_xxxzzz_xxxyzz_1, \
                                         tg_xxxzzz_xxxzz_1, tg_xxxzzz_xxxzzz_0, tg_xxxzzz_xxxzzz_1, tg_xxxzzz_xxyyy_1, \
                                         tg_xxxzzz_xxyyyy_0, tg_xxxzzz_xxyyyy_1, tg_xxxzzz_xxyyyz_0, tg_xxxzzz_xxyyyz_1, \
                                         tg_xxxzzz_xxyyz_1, tg_xxxzzz_xxyyzz_0, tg_xxxzzz_xxyyzz_1, tg_xxxzzz_xxyzz_1, \
                                         tg_xxxzzz_xxyzzz_0, tg_xxxzzz_xxyzzz_1, tg_xxxzzz_xxzzz_1, tg_xxxzzz_xxzzzz_0, \
                                         tg_xxxzzz_xxzzzz_1, tg_xxxzzz_xyyyy_1, tg_xxxzzz_xyyyyy_0, tg_xxxzzz_xyyyyy_1, \
                                         tg_xxxzzz_xyyyyz_0, tg_xxxzzz_xyyyyz_1, tg_xxxzzz_xyyyz_1, tg_xxxzzz_xyyyzz_0, \
                                         tg_xxxzzz_xyyyzz_1, tg_xxxzzz_xyyzz_1, tg_xxxzzz_xyyzzz_0, tg_xxxzzz_xyyzzz_1, \
                                         tg_xxxzzz_xyzzz_1, tg_xxxzzz_xyzzzz_0, tg_xxxzzz_xyzzzz_1, tg_xxxzzz_xzzzz_1, \
                                         tg_xxxzzz_xzzzzz_0, tg_xxxzzz_xzzzzz_1, tg_xxxzzz_yyyyy_1, tg_xxxzzz_yyyyyy_0, \
                                         tg_xxxzzz_yyyyyy_1, tg_xxxzzz_yyyyyz_0, tg_xxxzzz_yyyyyz_1, tg_xxxzzz_yyyyz_1, \
                                         tg_xxxzzz_yyyyzz_0, tg_xxxzzz_yyyyzz_1, tg_xxxzzz_yyyzz_1, tg_xxxzzz_yyzzz_1, \
                                         tg_xxxzzz_yzzzz_1, tg_xxxzzz_zzzzz_1, tg_xxyyy_xyyyyz_0, tg_xxyyy_xyyyyz_1, \
                                         tg_xxyyy_xyyyzz_0, tg_xxyyy_xyyyzz_1, tg_xxyyy_xyyzzz_0, tg_xxyyy_xyyzzz_1, \
                                         tg_xxyyy_xyzzzz_0, tg_xxyyy_xyzzzz_1, tg_xxyyy_xzzzzz_0, tg_xxyyy_xzzzzz_1, \
                                         tg_xxyyy_yyyyyy_0, tg_xxyyy_yyyyyy_1, tg_xxyyy_yyyyyz_0, tg_xxyyy_yyyyyz_1, \
                                         tg_xxyyy_yyyyzz_0, tg_xxyyy_yyyyzz_1, tg_xxyyy_yyyzzz_0, tg_xxyyy_yyyzzz_1, \
                                         tg_xxyyy_yyzzzz_0, tg_xxyyy_yyzzzz_1, tg_xxyyy_yzzzzz_0, tg_xxyyy_yzzzzz_1, \
                                         tg_xxyyy_zzzzzz_0, tg_xxyyy_zzzzzz_1, tg_xxyyz_xxxxxx_0, tg_xxyyz_xxxxxx_1, \
                                         tg_xxyyz_xxxxxy_0, tg_xxyyz_xxxxxy_1, tg_xxyyz_xxxxxz_0, tg_xxyyz_xxxxxz_1, \
                                         tg_xxyyz_xxxxyy_0, tg_xxyyz_xxxxyy_1, tg_xxyyz_xxxxyz_0, tg_xxyyz_xxxxyz_1, \
                                         tg_xxyyz_xxxxzz_0, tg_xxyyz_xxxxzz_1, tg_xxyyz_xxxyyy_0, tg_xxyyz_xxxyyy_1, \
                                         tg_xxyyz_xxxyyz_0, tg_xxyyz_xxxyyz_1, tg_xxyyz_xxxyzz_0, tg_xxyyz_xxxyzz_1, \
                                         tg_xxyyz_xxxzzz_0, tg_xxyyz_xxxzzz_1, tg_xxyyz_xxyyyy_0, tg_xxyyz_xxyyyy_1, \
                                         tg_xxyyz_xxyyyz_0, tg_xxyyz_xxyyyz_1, tg_xxyyz_xxyyzz_0, tg_xxyyz_xxyyzz_1, \
                                         tg_xxyyz_xxyzzz_0, tg_xxyyz_xxyzzz_1, tg_xxyyz_xxzzzz_0, tg_xxyyz_xxzzzz_1, \
                                         tg_xxyyz_xyyyyy_0, tg_xxyyz_xyyyyy_1, tg_xxyyz_xyyyyz_0, tg_xxyyz_xyyyyz_1, \
                                         tg_xxyyz_xyyyzz_0, tg_xxyyz_xyyyzz_1, tg_xxyyz_xyyzzz_0, tg_xxyyz_xyyzzz_1, \
                                         tg_xxyyz_xyzzzz_0, tg_xxyyz_xyzzzz_1, tg_xxyyz_xzzzzz_0, tg_xxyyz_xzzzzz_1, \
                                         tg_xxyyz_yyyyyy_0, tg_xxyyz_yyyyyy_1, tg_xxyyz_yyyyyz_0, tg_xxyyz_yyyyyz_1, \
                                         tg_xxyyz_yyyyzz_0, tg_xxyyz_yyyyzz_1, tg_xxyyz_yyyzzz_0, tg_xxyyz_yyyzzz_1, \
                                         tg_xxyyz_yyzzzz_0, tg_xxyyz_yyzzzz_1, tg_xxyyz_yzzzzz_0, tg_xxyyz_yzzzzz_1, \
                                         tg_xxyyz_zzzzzz_0, tg_xxyyz_zzzzzz_1, tg_xxyzz_xxxxxx_0, tg_xxyzz_xxxxxx_1, \
                                         tg_xxyzz_xxxxxy_0, tg_xxyzz_xxxxxy_1, tg_xxyzz_xxxxxz_0, tg_xxyzz_xxxxxz_1, \
                                         tg_xxyzz_xxxxyy_0, tg_xxyzz_xxxxyy_1, tg_xxyzz_xxxxyz_0, tg_xxyzz_xxxxyz_1, \
                                         tg_xxyzz_xxxxzz_0, tg_xxyzz_xxxxzz_1, tg_xxyzz_xxxyyy_0, tg_xxyzz_xxxyyy_1, \
                                         tg_xxyzz_xxxyyz_0, tg_xxyzz_xxxyyz_1, tg_xxyzz_xxxyzz_0, tg_xxyzz_xxxyzz_1, \
                                         tg_xxyzz_xxxzzz_0, tg_xxyzz_xxxzzz_1, tg_xxyzz_xxyyyy_0, tg_xxyzz_xxyyyy_1, \
                                         tg_xxyzz_xxyyyz_0, tg_xxyzz_xxyyyz_1, tg_xxyzz_xxyyzz_0, tg_xxyzz_xxyyzz_1, \
                                         tg_xxyzz_xxyzzz_0, tg_xxyzz_xxyzzz_1, tg_xxyzz_xxzzzz_0, tg_xxyzz_xxzzzz_1, \
                                         tg_xxyzz_xyyyyy_0, tg_xxyzz_xyyyyy_1, tg_xxyzz_xyyyyz_0, tg_xxyzz_xyyyyz_1, \
                                         tg_xxyzz_xyyyzz_0, tg_xxyzz_xyyyzz_1, tg_xxyzz_xyyzzz_0, tg_xxyzz_xyyzzz_1, \
                                         tg_xxyzz_xyzzzz_0, tg_xxyzz_xyzzzz_1, tg_xxyzz_xzzzzz_0, tg_xxyzz_xzzzzz_1, \
                                         tg_xxyzz_yyyyyy_0, tg_xxyzz_yyyyyy_1, tg_xxyzz_yyyyyz_0, tg_xxyzz_yyyyyz_1, \
                                         tg_xxyzz_yyyyzz_0, tg_xxyzz_yyyyzz_1, tg_xxyzz_yyyzzz_0, tg_xxyzz_yyyzzz_1, \
                                         tg_xxyzz_yyzzzz_0, tg_xxyzz_yyzzzz_1, tg_xxyzz_yzzzzz_0, tg_xxyzz_yzzzzz_1, \
                                         tg_xxyzz_zzzzzz_0, tg_xxyzz_zzzzzz_1, tg_xxzzz_xxxxxx_0, tg_xxzzz_xxxxxx_1, \
                                         tg_xxzzz_xxxxxy_0, tg_xxzzz_xxxxxy_1, tg_xxzzz_xxxxxz_0, tg_xxzzz_xxxxxz_1, \
                                         tg_xxzzz_xxxxyy_0, tg_xxzzz_xxxxyy_1, tg_xxzzz_xxxxyz_0, tg_xxzzz_xxxxyz_1, \
                                         tg_xxzzz_xxxxzz_0, tg_xxzzz_xxxxzz_1, tg_xxzzz_xxxyyy_0, tg_xxzzz_xxxyyy_1, \
                                         tg_xxzzz_xxxyyz_0, tg_xxzzz_xxxyyz_1, tg_xxzzz_xxxyzz_0, tg_xxzzz_xxxyzz_1, \
                                         tg_xxzzz_xxxzzz_0, tg_xxzzz_xxxzzz_1, tg_xxzzz_xxyyyy_0, tg_xxzzz_xxyyyy_1, \
                                         tg_xxzzz_xxyyyz_0, tg_xxzzz_xxyyyz_1, tg_xxzzz_xxyyzz_0, tg_xxzzz_xxyyzz_1, \
                                         tg_xxzzz_xxyzzz_0, tg_xxzzz_xxyzzz_1, tg_xxzzz_xxzzzz_0, tg_xxzzz_xxzzzz_1, \
                                         tg_xxzzz_xyyyyy_0, tg_xxzzz_xyyyyy_1, tg_xxzzz_xyyyyz_0, tg_xxzzz_xyyyyz_1, \
                                         tg_xxzzz_xyyyzz_0, tg_xxzzz_xyyyzz_1, tg_xxzzz_xyyzzz_0, tg_xxzzz_xyyzzz_1, \
                                         tg_xxzzz_xyzzzz_0, tg_xxzzz_xyzzzz_1, tg_xxzzz_xzzzzz_0, tg_xxzzz_xzzzzz_1, \
                                         tg_xxzzz_yyyyyy_0, tg_xxzzz_yyyyyy_1, tg_xxzzz_yyyyyz_0, tg_xxzzz_yyyyyz_1, \
                                         tg_xxzzz_yyyyzz_0, tg_xxzzz_yyyyzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxyyy_xyyyyz_0[j] = pb_x * tg_xxxyyy_xyyyyz_0[j] + fr * tg_xxxyyy_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyyyz_0[j] - tg_xxyyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyyyz_1[j];

                    tg_xxxxyyy_xyyyzz_0[j] = pb_x * tg_xxxyyy_xyyyzz_0[j] + fr * tg_xxxyyy_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyyzz_0[j] - tg_xxyyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyyzz_1[j];

                    tg_xxxxyyy_xyyzzz_0[j] = pb_x * tg_xxxyyy_xyyzzz_0[j] + fr * tg_xxxyyy_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyzzz_0[j] - tg_xxyyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyzzz_1[j];

                    tg_xxxxyyy_xyzzzz_0[j] = pb_x * tg_xxxyyy_xyzzzz_0[j] + fr * tg_xxxyyy_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyzzzz_0[j] - tg_xxyyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yzzzz_1[j];

                    tg_xxxxyyy_xzzzzz_0[j] = pb_x * tg_xxxyyy_xzzzzz_0[j] + fr * tg_xxxyyy_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xzzzzz_0[j] - tg_xxyyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_zzzzz_1[j];

                    tg_xxxxyyy_yyyyyy_0[j] = pb_x * tg_xxxyyy_yyyyyy_0[j] + fr * tg_xxxyyy_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyyyy_0[j] - tg_xxyyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxyyy_yyyyyz_0[j] = pb_x * tg_xxxyyy_yyyyyz_0[j] + fr * tg_xxxyyy_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyyyz_0[j] - tg_xxyyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxyyy_yyyyzz_0[j] = pb_x * tg_xxxyyy_yyyyzz_0[j] + fr * tg_xxxyyy_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyyzz_0[j] - tg_xxyyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxyyy_yyyzzz_0[j] = pb_x * tg_xxxyyy_yyyzzz_0[j] + fr * tg_xxxyyy_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyzzz_0[j] - tg_xxyyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_yyzzzz_0[j] = pb_x * tg_xxxyyy_yyzzzz_0[j] + fr * tg_xxxyyy_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyzzzz_0[j] - tg_xxyyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_yzzzzz_0[j] = pb_x * tg_xxxyyy_yzzzzz_0[j] + fr * tg_xxxyyy_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yzzzzz_0[j] - tg_xxyyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_zzzzzz_0[j] = pb_x * tg_xxxyyy_zzzzzz_0[j] + fr * tg_xxxyyy_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_zzzzzz_0[j] - tg_xxyyy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_xxxxxx_0[j] = pb_x * tg_xxxyyz_xxxxxx_0[j] + fr * tg_xxxyyz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxxx_0[j] - tg_xxyyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxyyz_xxxxx_1[j];

                    tg_xxxxyyz_xxxxxy_0[j] = pb_x * tg_xxxyyz_xxxxxy_0[j] + fr * tg_xxxyyz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxxy_0[j] - tg_xxyyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyz_xxxxy_1[j];

                    tg_xxxxyyz_xxxxxz_0[j] = pb_x * tg_xxxyyz_xxxxxz_0[j] + fr * tg_xxxyyz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxxz_0[j] - tg_xxyyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyyz_xxxxz_1[j];

                    tg_xxxxyyz_xxxxyy_0[j] = pb_x * tg_xxxyyz_xxxxyy_0[j] + fr * tg_xxxyyz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxyy_0[j] - tg_xxyyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyz_xxxyy_1[j];

                    tg_xxxxyyz_xxxxyz_0[j] = pb_x * tg_xxxyyz_xxxxyz_0[j] + fr * tg_xxxyyz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxyz_0[j] - tg_xxyyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyz_xxxyz_1[j];

                    tg_xxxxyyz_xxxxzz_0[j] = pb_x * tg_xxxyyz_xxxxzz_0[j] + fr * tg_xxxyyz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxxzz_0[j] - tg_xxyyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyz_xxxzz_1[j];

                    tg_xxxxyyz_xxxyyy_0[j] = pb_x * tg_xxxyyz_xxxyyy_0[j] + fr * tg_xxxyyz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxyyy_0[j] - tg_xxyyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxyyy_1[j];

                    tg_xxxxyyz_xxxyyz_0[j] = pb_x * tg_xxxyyz_xxxyyz_0[j] + fr * tg_xxxyyz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxyyz_0[j] - tg_xxyyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxyyz_1[j];

                    tg_xxxxyyz_xxxyzz_0[j] = pb_x * tg_xxxyyz_xxxyzz_0[j] + fr * tg_xxxyyz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxyzz_0[j] - tg_xxyyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxyzz_1[j];

                    tg_xxxxyyz_xxxzzz_0[j] = pb_x * tg_xxxyyz_xxxzzz_0[j] + fr * tg_xxxyyz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxzzz_0[j] - tg_xxyyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxzzz_1[j];

                    tg_xxxxyyz_xxyyyy_0[j] = pb_x * tg_xxxyyz_xxyyyy_0[j] + fr * tg_xxxyyz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyyyy_0[j] - tg_xxyyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyyyy_1[j];

                    tg_xxxxyyz_xxyyyz_0[j] = pb_x * tg_xxxyyz_xxyyyz_0[j] + fr * tg_xxxyyz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyyyz_0[j] - tg_xxyyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyyyz_1[j];

                    tg_xxxxyyz_xxyyzz_0[j] = pb_x * tg_xxxyyz_xxyyzz_0[j] + fr * tg_xxxyyz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyyzz_0[j] - tg_xxyyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyyzz_1[j];

                    tg_xxxxyyz_xxyzzz_0[j] = pb_x * tg_xxxyyz_xxyzzz_0[j] + fr * tg_xxxyyz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyzzz_0[j] - tg_xxyyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyzzz_1[j];

                    tg_xxxxyyz_xxzzzz_0[j] = pb_x * tg_xxxyyz_xxzzzz_0[j] + fr * tg_xxxyyz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxzzzz_0[j] - tg_xxyyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xzzzz_1[j];

                    tg_xxxxyyz_xyyyyy_0[j] = pb_x * tg_xxxyyz_xyyyyy_0[j] + fr * tg_xxxyyz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyyyy_0[j] - tg_xxyyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyyyy_1[j];

                    tg_xxxxyyz_xyyyyz_0[j] = pb_x * tg_xxxyyz_xyyyyz_0[j] + fr * tg_xxxyyz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyyyz_0[j] - tg_xxyyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyyyz_1[j];

                    tg_xxxxyyz_xyyyzz_0[j] = pb_x * tg_xxxyyz_xyyyzz_0[j] + fr * tg_xxxyyz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyyzz_0[j] - tg_xxyyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyyzz_1[j];

                    tg_xxxxyyz_xyyzzz_0[j] = pb_x * tg_xxxyyz_xyyzzz_0[j] + fr * tg_xxxyyz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyzzz_0[j] - tg_xxyyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyzzz_1[j];

                    tg_xxxxyyz_xyzzzz_0[j] = pb_x * tg_xxxyyz_xyzzzz_0[j] + fr * tg_xxxyyz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyzzzz_0[j] - tg_xxyyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yzzzz_1[j];

                    tg_xxxxyyz_xzzzzz_0[j] = pb_x * tg_xxxyyz_xzzzzz_0[j] + fr * tg_xxxyyz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xzzzzz_0[j] - tg_xxyyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_zzzzz_1[j];

                    tg_xxxxyyz_yyyyyy_0[j] = pb_x * tg_xxxyyz_yyyyyy_0[j] + fr * tg_xxxyyz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyyyy_0[j] - tg_xxyyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxyyz_yyyyyz_0[j] = pb_x * tg_xxxyyz_yyyyyz_0[j] + fr * tg_xxxyyz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyyyz_0[j] - tg_xxyyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxyyz_yyyyzz_0[j] = pb_x * tg_xxxyyz_yyyyzz_0[j] + fr * tg_xxxyyz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyyzz_0[j] - tg_xxyyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxyyz_yyyzzz_0[j] = pb_x * tg_xxxyyz_yyyzzz_0[j] + fr * tg_xxxyyz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyzzz_0[j] - tg_xxyyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_yyzzzz_0[j] = pb_x * tg_xxxyyz_yyzzzz_0[j] + fr * tg_xxxyyz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyzzzz_0[j] - tg_xxyyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_yzzzzz_0[j] = pb_x * tg_xxxyyz_yzzzzz_0[j] + fr * tg_xxxyyz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yzzzzz_0[j] - tg_xxyyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_zzzzzz_0[j] = pb_x * tg_xxxyyz_zzzzzz_0[j] + fr * tg_xxxyyz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_zzzzzz_0[j] - tg_xxyyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_xxxxxx_0[j] = pb_x * tg_xxxyzz_xxxxxx_0[j] + fr * tg_xxxyzz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxxx_0[j] - tg_xxyzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxyzz_xxxxx_1[j];

                    tg_xxxxyzz_xxxxxy_0[j] = pb_x * tg_xxxyzz_xxxxxy_0[j] + fr * tg_xxxyzz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxxy_0[j] - tg_xxyzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyzz_xxxxy_1[j];

                    tg_xxxxyzz_xxxxxz_0[j] = pb_x * tg_xxxyzz_xxxxxz_0[j] + fr * tg_xxxyzz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxxz_0[j] - tg_xxyzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxyzz_xxxxz_1[j];

                    tg_xxxxyzz_xxxxyy_0[j] = pb_x * tg_xxxyzz_xxxxyy_0[j] + fr * tg_xxxyzz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxyy_0[j] - tg_xxyzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzz_xxxyy_1[j];

                    tg_xxxxyzz_xxxxyz_0[j] = pb_x * tg_xxxyzz_xxxxyz_0[j] + fr * tg_xxxyzz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxyz_0[j] - tg_xxyzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzz_xxxyz_1[j];

                    tg_xxxxyzz_xxxxzz_0[j] = pb_x * tg_xxxyzz_xxxxzz_0[j] + fr * tg_xxxyzz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxxzz_0[j] - tg_xxyzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzz_xxxzz_1[j];

                    tg_xxxxyzz_xxxyyy_0[j] = pb_x * tg_xxxyzz_xxxyyy_0[j] + fr * tg_xxxyzz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxyyy_0[j] - tg_xxyzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxyyy_1[j];

                    tg_xxxxyzz_xxxyyz_0[j] = pb_x * tg_xxxyzz_xxxyyz_0[j] + fr * tg_xxxyzz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxyyz_0[j] - tg_xxyzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxyyz_1[j];

                    tg_xxxxyzz_xxxyzz_0[j] = pb_x * tg_xxxyzz_xxxyzz_0[j] + fr * tg_xxxyzz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxyzz_0[j] - tg_xxyzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxyzz_1[j];

                    tg_xxxxyzz_xxxzzz_0[j] = pb_x * tg_xxxyzz_xxxzzz_0[j] + fr * tg_xxxyzz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxzzz_0[j] - tg_xxyzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxzzz_1[j];

                    tg_xxxxyzz_xxyyyy_0[j] = pb_x * tg_xxxyzz_xxyyyy_0[j] + fr * tg_xxxyzz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyyyy_0[j] - tg_xxyzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyyyy_1[j];

                    tg_xxxxyzz_xxyyyz_0[j] = pb_x * tg_xxxyzz_xxyyyz_0[j] + fr * tg_xxxyzz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyyyz_0[j] - tg_xxyzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyyyz_1[j];

                    tg_xxxxyzz_xxyyzz_0[j] = pb_x * tg_xxxyzz_xxyyzz_0[j] + fr * tg_xxxyzz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyyzz_0[j] - tg_xxyzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyyzz_1[j];

                    tg_xxxxyzz_xxyzzz_0[j] = pb_x * tg_xxxyzz_xxyzzz_0[j] + fr * tg_xxxyzz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyzzz_0[j] - tg_xxyzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyzzz_1[j];

                    tg_xxxxyzz_xxzzzz_0[j] = pb_x * tg_xxxyzz_xxzzzz_0[j] + fr * tg_xxxyzz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxzzzz_0[j] - tg_xxyzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xzzzz_1[j];

                    tg_xxxxyzz_xyyyyy_0[j] = pb_x * tg_xxxyzz_xyyyyy_0[j] + fr * tg_xxxyzz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyyyy_0[j] - tg_xxyzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyyyy_1[j];

                    tg_xxxxyzz_xyyyyz_0[j] = pb_x * tg_xxxyzz_xyyyyz_0[j] + fr * tg_xxxyzz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyyyz_0[j] - tg_xxyzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyyyz_1[j];

                    tg_xxxxyzz_xyyyzz_0[j] = pb_x * tg_xxxyzz_xyyyzz_0[j] + fr * tg_xxxyzz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyyzz_0[j] - tg_xxyzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyyzz_1[j];

                    tg_xxxxyzz_xyyzzz_0[j] = pb_x * tg_xxxyzz_xyyzzz_0[j] + fr * tg_xxxyzz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyzzz_0[j] - tg_xxyzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyzzz_1[j];

                    tg_xxxxyzz_xyzzzz_0[j] = pb_x * tg_xxxyzz_xyzzzz_0[j] + fr * tg_xxxyzz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyzzzz_0[j] - tg_xxyzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yzzzz_1[j];

                    tg_xxxxyzz_xzzzzz_0[j] = pb_x * tg_xxxyzz_xzzzzz_0[j] + fr * tg_xxxyzz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xzzzzz_0[j] - tg_xxyzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_zzzzz_1[j];

                    tg_xxxxyzz_yyyyyy_0[j] = pb_x * tg_xxxyzz_yyyyyy_0[j] + fr * tg_xxxyzz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyyyy_0[j] - tg_xxyzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxyzz_yyyyyz_0[j] = pb_x * tg_xxxyzz_yyyyyz_0[j] + fr * tg_xxxyzz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyyyz_0[j] - tg_xxyzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxyzz_yyyyzz_0[j] = pb_x * tg_xxxyzz_yyyyzz_0[j] + fr * tg_xxxyzz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyyzz_0[j] - tg_xxyzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxyzz_yyyzzz_0[j] = pb_x * tg_xxxyzz_yyyzzz_0[j] + fr * tg_xxxyzz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyzzz_0[j] - tg_xxyzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_yyzzzz_0[j] = pb_x * tg_xxxyzz_yyzzzz_0[j] + fr * tg_xxxyzz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyzzzz_0[j] - tg_xxyzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_yzzzzz_0[j] = pb_x * tg_xxxyzz_yzzzzz_0[j] + fr * tg_xxxyzz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yzzzzz_0[j] - tg_xxyzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_zzzzzz_0[j] = pb_x * tg_xxxyzz_zzzzzz_0[j] + fr * tg_xxxyzz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_zzzzzz_0[j] - tg_xxyzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_xxxxxx_0[j] = pb_x * tg_xxxzzz_xxxxxx_0[j] + fr * tg_xxxzzz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxxx_0[j] - tg_xxzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxzzz_xxxxx_1[j];

                    tg_xxxxzzz_xxxxxy_0[j] = pb_x * tg_xxxzzz_xxxxxy_0[j] + fr * tg_xxxzzz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxxy_0[j] - tg_xxzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxzzz_xxxxy_1[j];

                    tg_xxxxzzz_xxxxxz_0[j] = pb_x * tg_xxxzzz_xxxxxz_0[j] + fr * tg_xxxzzz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxxz_0[j] - tg_xxzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxzzz_xxxxz_1[j];

                    tg_xxxxzzz_xxxxyy_0[j] = pb_x * tg_xxxzzz_xxxxyy_0[j] + fr * tg_xxxzzz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxyy_0[j] - tg_xxzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzz_xxxyy_1[j];

                    tg_xxxxzzz_xxxxyz_0[j] = pb_x * tg_xxxzzz_xxxxyz_0[j] + fr * tg_xxxzzz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxyz_0[j] - tg_xxzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzz_xxxyz_1[j];

                    tg_xxxxzzz_xxxxzz_0[j] = pb_x * tg_xxxzzz_xxxxzz_0[j] + fr * tg_xxxzzz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxxzz_0[j] - tg_xxzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzz_xxxzz_1[j];

                    tg_xxxxzzz_xxxyyy_0[j] = pb_x * tg_xxxzzz_xxxyyy_0[j] + fr * tg_xxxzzz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxyyy_0[j] - tg_xxzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxyyy_1[j];

                    tg_xxxxzzz_xxxyyz_0[j] = pb_x * tg_xxxzzz_xxxyyz_0[j] + fr * tg_xxxzzz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxyyz_0[j] - tg_xxzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxyyz_1[j];

                    tg_xxxxzzz_xxxyzz_0[j] = pb_x * tg_xxxzzz_xxxyzz_0[j] + fr * tg_xxxzzz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxyzz_0[j] - tg_xxzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxyzz_1[j];

                    tg_xxxxzzz_xxxzzz_0[j] = pb_x * tg_xxxzzz_xxxzzz_0[j] + fr * tg_xxxzzz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxzzz_0[j] - tg_xxzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxzzz_1[j];

                    tg_xxxxzzz_xxyyyy_0[j] = pb_x * tg_xxxzzz_xxyyyy_0[j] + fr * tg_xxxzzz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyyyy_0[j] - tg_xxzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyyyy_1[j];

                    tg_xxxxzzz_xxyyyz_0[j] = pb_x * tg_xxxzzz_xxyyyz_0[j] + fr * tg_xxxzzz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyyyz_0[j] - tg_xxzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyyyz_1[j];

                    tg_xxxxzzz_xxyyzz_0[j] = pb_x * tg_xxxzzz_xxyyzz_0[j] + fr * tg_xxxzzz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyyzz_0[j] - tg_xxzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyyzz_1[j];

                    tg_xxxxzzz_xxyzzz_0[j] = pb_x * tg_xxxzzz_xxyzzz_0[j] + fr * tg_xxxzzz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyzzz_0[j] - tg_xxzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyzzz_1[j];

                    tg_xxxxzzz_xxzzzz_0[j] = pb_x * tg_xxxzzz_xxzzzz_0[j] + fr * tg_xxxzzz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxzzzz_0[j] - tg_xxzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xzzzz_1[j];

                    tg_xxxxzzz_xyyyyy_0[j] = pb_x * tg_xxxzzz_xyyyyy_0[j] + fr * tg_xxxzzz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyyyy_0[j] - tg_xxzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyyyy_1[j];

                    tg_xxxxzzz_xyyyyz_0[j] = pb_x * tg_xxxzzz_xyyyyz_0[j] + fr * tg_xxxzzz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyyyz_0[j] - tg_xxzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyyyz_1[j];

                    tg_xxxxzzz_xyyyzz_0[j] = pb_x * tg_xxxzzz_xyyyzz_0[j] + fr * tg_xxxzzz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyyzz_0[j] - tg_xxzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyyzz_1[j];

                    tg_xxxxzzz_xyyzzz_0[j] = pb_x * tg_xxxzzz_xyyzzz_0[j] + fr * tg_xxxzzz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyzzz_0[j] - tg_xxzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyzzz_1[j];

                    tg_xxxxzzz_xyzzzz_0[j] = pb_x * tg_xxxzzz_xyzzzz_0[j] + fr * tg_xxxzzz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyzzzz_0[j] - tg_xxzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yzzzz_1[j];

                    tg_xxxxzzz_xzzzzz_0[j] = pb_x * tg_xxxzzz_xzzzzz_0[j] + fr * tg_xxxzzz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xzzzzz_0[j] - tg_xxzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_zzzzz_1[j];

                    tg_xxxxzzz_yyyyyy_0[j] = pb_x * tg_xxxzzz_yyyyyy_0[j] + fr * tg_xxxzzz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyyyy_0[j] - tg_xxzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxzzz_yyyyyz_0[j] = pb_x * tg_xxxzzz_yyyyyz_0[j] + fr * tg_xxxzzz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyyyz_0[j] - tg_xxzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxzzz_yyyyzz_0[j] = pb_x * tg_xxxzzz_yyyyzz_0[j] + fr * tg_xxxzzz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyyzz_0[j] - tg_xxzzz_yyyyzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_276_368(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (276,368)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxxzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 276); 

                auto tg_xxxzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 277); 

                auto tg_xxxzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 278); 

                auto tg_xxxzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 279); 

                auto tg_xxyyyy_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 280); 

                auto tg_xxyyyy_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 281); 

                auto tg_xxyyyy_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 282); 

                auto tg_xxyyyy_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 283); 

                auto tg_xxyyyy_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 284); 

                auto tg_xxyyyy_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 285); 

                auto tg_xxyyyy_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 286); 

                auto tg_xxyyyy_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 287); 

                auto tg_xxyyyy_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 288); 

                auto tg_xxyyyy_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 289); 

                auto tg_xxyyyy_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 290); 

                auto tg_xxyyyy_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 291); 

                auto tg_xxyyyy_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 292); 

                auto tg_xxyyyy_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 293); 

                auto tg_xxyyyy_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 294); 

                auto tg_xxyyyy_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 295); 

                auto tg_xxyyyy_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 296); 

                auto tg_xxyyyy_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 297); 

                auto tg_xxyyyy_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 298); 

                auto tg_xxyyyy_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 299); 

                auto tg_xxyyyy_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 300); 

                auto tg_xxyyyy_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 301); 

                auto tg_xxyyyy_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 302); 

                auto tg_xxyyyy_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 303); 

                auto tg_xxyyyy_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 304); 

                auto tg_xxyyyy_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 305); 

                auto tg_xxyyyy_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 306); 

                auto tg_xxyyyy_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 307); 

                auto tg_xxyyyz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 308); 

                auto tg_xxyyyz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 309); 

                auto tg_xxyyyz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 310); 

                auto tg_xxyyyz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 311); 

                auto tg_xxyyyz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 312); 

                auto tg_xxyyyz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 313); 

                auto tg_xxyyyz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 314); 

                auto tg_xxyyyz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 315); 

                auto tg_xxyyyz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 316); 

                auto tg_xxyyyz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 317); 

                auto tg_xxyyyz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 318); 

                auto tg_xxyyyz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 319); 

                auto tg_xxyyyz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 320); 

                auto tg_xxyyyz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 321); 

                auto tg_xxyyyz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 322); 

                auto tg_xxyyyz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 323); 

                auto tg_xxyyyz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 324); 

                auto tg_xxyyyz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 325); 

                auto tg_xxyyyz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 326); 

                auto tg_xxyyyz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 327); 

                auto tg_xxyyyz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 328); 

                auto tg_xxyyyz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 329); 

                auto tg_xxyyyz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 330); 

                auto tg_xxyyyz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 331); 

                auto tg_xxyyyz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 332); 

                auto tg_xxyyyz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 333); 

                auto tg_xxyyyz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 334); 

                auto tg_xxyyyz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 335); 

                auto tg_xxyyzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 336); 

                auto tg_xxyyzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 337); 

                auto tg_xxyyzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 338); 

                auto tg_xxyyzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 339); 

                auto tg_xxyyzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 340); 

                auto tg_xxyyzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 341); 

                auto tg_xxyyzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 342); 

                auto tg_xxyyzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 343); 

                auto tg_xxyyzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 344); 

                auto tg_xxyyzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 345); 

                auto tg_xxyyzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 346); 

                auto tg_xxyyzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 347); 

                auto tg_xxyyzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 348); 

                auto tg_xxyyzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 349); 

                auto tg_xxyyzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 350); 

                auto tg_xxyyzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 351); 

                auto tg_xxyyzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 352); 

                auto tg_xxyyzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 353); 

                auto tg_xxyyzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 354); 

                auto tg_xxyyzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 355); 

                auto tg_xxyyzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 356); 

                auto tg_xxyyzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 357); 

                auto tg_xxyyzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 358); 

                auto tg_xxyyzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 359); 

                auto tg_xxyyzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 360); 

                auto tg_xxyyzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 361); 

                auto tg_xxyyzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 362); 

                auto tg_xxyyzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 363); 

                auto tg_xxyzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 364); 

                auto tg_xxyzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 365); 

                auto tg_xxyzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 366); 

                auto tg_xxyzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 367); 

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

                // set up pointers to integrals

                auto tg_xxxxzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 276); 

                auto tg_xxxxzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 277); 

                auto tg_xxxxzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 278); 

                auto tg_xxxxzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 279); 

                auto tg_xxxyyyy_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 280); 

                auto tg_xxxyyyy_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 281); 

                auto tg_xxxyyyy_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 282); 

                auto tg_xxxyyyy_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 283); 

                auto tg_xxxyyyy_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 284); 

                auto tg_xxxyyyy_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 285); 

                auto tg_xxxyyyy_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 286); 

                auto tg_xxxyyyy_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 287); 

                auto tg_xxxyyyy_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 288); 

                auto tg_xxxyyyy_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 289); 

                auto tg_xxxyyyy_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 290); 

                auto tg_xxxyyyy_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 291); 

                auto tg_xxxyyyy_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 292); 

                auto tg_xxxyyyy_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 293); 

                auto tg_xxxyyyy_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 294); 

                auto tg_xxxyyyy_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 295); 

                auto tg_xxxyyyy_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 296); 

                auto tg_xxxyyyy_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 297); 

                auto tg_xxxyyyy_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 298); 

                auto tg_xxxyyyy_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 299); 

                auto tg_xxxyyyy_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 300); 

                auto tg_xxxyyyy_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 301); 

                auto tg_xxxyyyy_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 302); 

                auto tg_xxxyyyy_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 303); 

                auto tg_xxxyyyy_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 304); 

                auto tg_xxxyyyy_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 305); 

                auto tg_xxxyyyy_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 306); 

                auto tg_xxxyyyy_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 307); 

                auto tg_xxxyyyz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 308); 

                auto tg_xxxyyyz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 309); 

                auto tg_xxxyyyz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 310); 

                auto tg_xxxyyyz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 311); 

                auto tg_xxxyyyz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 312); 

                auto tg_xxxyyyz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 313); 

                auto tg_xxxyyyz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 314); 

                auto tg_xxxyyyz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 315); 

                auto tg_xxxyyyz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 316); 

                auto tg_xxxyyyz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 317); 

                auto tg_xxxyyyz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 318); 

                auto tg_xxxyyyz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 319); 

                auto tg_xxxyyyz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 320); 

                auto tg_xxxyyyz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 321); 

                auto tg_xxxyyyz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 322); 

                auto tg_xxxyyyz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 323); 

                auto tg_xxxyyyz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 324); 

                auto tg_xxxyyyz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 325); 

                auto tg_xxxyyyz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 326); 

                auto tg_xxxyyyz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 327); 

                auto tg_xxxyyyz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 328); 

                auto tg_xxxyyyz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 329); 

                auto tg_xxxyyyz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 330); 

                auto tg_xxxyyyz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 331); 

                auto tg_xxxyyyz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 332); 

                auto tg_xxxyyyz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 333); 

                auto tg_xxxyyyz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 334); 

                auto tg_xxxyyyz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 335); 

                auto tg_xxxyyzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 336); 

                auto tg_xxxyyzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 337); 

                auto tg_xxxyyzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 338); 

                auto tg_xxxyyzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 339); 

                auto tg_xxxyyzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 340); 

                auto tg_xxxyyzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 341); 

                auto tg_xxxyyzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 342); 

                auto tg_xxxyyzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 343); 

                auto tg_xxxyyzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 344); 

                auto tg_xxxyyzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 345); 

                auto tg_xxxyyzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 346); 

                auto tg_xxxyyzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 347); 

                auto tg_xxxyyzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 348); 

                auto tg_xxxyyzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 349); 

                auto tg_xxxyyzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 350); 

                auto tg_xxxyyzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 351); 

                auto tg_xxxyyzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 352); 

                auto tg_xxxyyzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 353); 

                auto tg_xxxyyzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 354); 

                auto tg_xxxyyzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 355); 

                auto tg_xxxyyzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 356); 

                auto tg_xxxyyzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 357); 

                auto tg_xxxyyzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 358); 

                auto tg_xxxyyzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 359); 

                auto tg_xxxyyzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 360); 

                auto tg_xxxyyzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 361); 

                auto tg_xxxyyzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 362); 

                auto tg_xxxyyzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 363); 

                auto tg_xxxyzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 364); 

                auto tg_xxxyzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 365); 

                auto tg_xxxyzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 366); 

                auto tg_xxxyzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 367); 

                // Batch of Integrals (276,368)

                #pragma omp simd aligned(fxn, fza, tg_xxxxzzz_yyyzzz_0, tg_xxxxzzz_yyzzzz_0, \
                                         tg_xxxxzzz_yzzzzz_0, tg_xxxxzzz_zzzzzz_0, tg_xxxyyyy_xxxxxx_0, tg_xxxyyyy_xxxxxy_0, \
                                         tg_xxxyyyy_xxxxxz_0, tg_xxxyyyy_xxxxyy_0, tg_xxxyyyy_xxxxyz_0, tg_xxxyyyy_xxxxzz_0, \
                                         tg_xxxyyyy_xxxyyy_0, tg_xxxyyyy_xxxyyz_0, tg_xxxyyyy_xxxyzz_0, tg_xxxyyyy_xxxzzz_0, \
                                         tg_xxxyyyy_xxyyyy_0, tg_xxxyyyy_xxyyyz_0, tg_xxxyyyy_xxyyzz_0, tg_xxxyyyy_xxyzzz_0, \
                                         tg_xxxyyyy_xxzzzz_0, tg_xxxyyyy_xyyyyy_0, tg_xxxyyyy_xyyyyz_0, tg_xxxyyyy_xyyyzz_0, \
                                         tg_xxxyyyy_xyyzzz_0, tg_xxxyyyy_xyzzzz_0, tg_xxxyyyy_xzzzzz_0, tg_xxxyyyy_yyyyyy_0, \
                                         tg_xxxyyyy_yyyyyz_0, tg_xxxyyyy_yyyyzz_0, tg_xxxyyyy_yyyzzz_0, tg_xxxyyyy_yyzzzz_0, \
                                         tg_xxxyyyy_yzzzzz_0, tg_xxxyyyy_zzzzzz_0, tg_xxxyyyz_xxxxxx_0, tg_xxxyyyz_xxxxxy_0, \
                                         tg_xxxyyyz_xxxxxz_0, tg_xxxyyyz_xxxxyy_0, tg_xxxyyyz_xxxxyz_0, tg_xxxyyyz_xxxxzz_0, \
                                         tg_xxxyyyz_xxxyyy_0, tg_xxxyyyz_xxxyyz_0, tg_xxxyyyz_xxxyzz_0, tg_xxxyyyz_xxxzzz_0, \
                                         tg_xxxyyyz_xxyyyy_0, tg_xxxyyyz_xxyyyz_0, tg_xxxyyyz_xxyyzz_0, tg_xxxyyyz_xxyzzz_0, \
                                         tg_xxxyyyz_xxzzzz_0, tg_xxxyyyz_xyyyyy_0, tg_xxxyyyz_xyyyyz_0, tg_xxxyyyz_xyyyzz_0, \
                                         tg_xxxyyyz_xyyzzz_0, tg_xxxyyyz_xyzzzz_0, tg_xxxyyyz_xzzzzz_0, tg_xxxyyyz_yyyyyy_0, \
                                         tg_xxxyyyz_yyyyyz_0, tg_xxxyyyz_yyyyzz_0, tg_xxxyyyz_yyyzzz_0, tg_xxxyyyz_yyzzzz_0, \
                                         tg_xxxyyyz_yzzzzz_0, tg_xxxyyyz_zzzzzz_0, tg_xxxyyzz_xxxxxx_0, tg_xxxyyzz_xxxxxy_0, \
                                         tg_xxxyyzz_xxxxxz_0, tg_xxxyyzz_xxxxyy_0, tg_xxxyyzz_xxxxyz_0, tg_xxxyyzz_xxxxzz_0, \
                                         tg_xxxyyzz_xxxyyy_0, tg_xxxyyzz_xxxyyz_0, tg_xxxyyzz_xxxyzz_0, tg_xxxyyzz_xxxzzz_0, \
                                         tg_xxxyyzz_xxyyyy_0, tg_xxxyyzz_xxyyyz_0, tg_xxxyyzz_xxyyzz_0, tg_xxxyyzz_xxyzzz_0, \
                                         tg_xxxyyzz_xxzzzz_0, tg_xxxyyzz_xyyyyy_0, tg_xxxyyzz_xyyyyz_0, tg_xxxyyzz_xyyyzz_0, \
                                         tg_xxxyyzz_xyyzzz_0, tg_xxxyyzz_xyzzzz_0, tg_xxxyyzz_xzzzzz_0, tg_xxxyyzz_yyyyyy_0, \
                                         tg_xxxyyzz_yyyyyz_0, tg_xxxyyzz_yyyyzz_0, tg_xxxyyzz_yyyzzz_0, tg_xxxyyzz_yyzzzz_0, \
                                         tg_xxxyyzz_yzzzzz_0, tg_xxxyyzz_zzzzzz_0, tg_xxxyzzz_xxxxxx_0, tg_xxxyzzz_xxxxxy_0, \
                                         tg_xxxyzzz_xxxxxz_0, tg_xxxyzzz_xxxxyy_0, tg_xxxzzz_yyyzzz_0, tg_xxxzzz_yyyzzz_1, \
                                         tg_xxxzzz_yyzzzz_0, tg_xxxzzz_yyzzzz_1, tg_xxxzzz_yzzzzz_0, tg_xxxzzz_yzzzzz_1, \
                                         tg_xxxzzz_zzzzzz_0, tg_xxxzzz_zzzzzz_1, tg_xxyyyy_xxxxx_1, tg_xxyyyy_xxxxxx_0, \
                                         tg_xxyyyy_xxxxxx_1, tg_xxyyyy_xxxxxy_0, tg_xxyyyy_xxxxxy_1, tg_xxyyyy_xxxxxz_0, \
                                         tg_xxyyyy_xxxxxz_1, tg_xxyyyy_xxxxy_1, tg_xxyyyy_xxxxyy_0, tg_xxyyyy_xxxxyy_1, \
                                         tg_xxyyyy_xxxxyz_0, tg_xxyyyy_xxxxyz_1, tg_xxyyyy_xxxxz_1, tg_xxyyyy_xxxxzz_0, \
                                         tg_xxyyyy_xxxxzz_1, tg_xxyyyy_xxxyy_1, tg_xxyyyy_xxxyyy_0, tg_xxyyyy_xxxyyy_1, \
                                         tg_xxyyyy_xxxyyz_0, tg_xxyyyy_xxxyyz_1, tg_xxyyyy_xxxyz_1, tg_xxyyyy_xxxyzz_0, \
                                         tg_xxyyyy_xxxyzz_1, tg_xxyyyy_xxxzz_1, tg_xxyyyy_xxxzzz_0, tg_xxyyyy_xxxzzz_1, \
                                         tg_xxyyyy_xxyyy_1, tg_xxyyyy_xxyyyy_0, tg_xxyyyy_xxyyyy_1, tg_xxyyyy_xxyyyz_0, \
                                         tg_xxyyyy_xxyyyz_1, tg_xxyyyy_xxyyz_1, tg_xxyyyy_xxyyzz_0, tg_xxyyyy_xxyyzz_1, \
                                         tg_xxyyyy_xxyzz_1, tg_xxyyyy_xxyzzz_0, tg_xxyyyy_xxyzzz_1, tg_xxyyyy_xxzzz_1, \
                                         tg_xxyyyy_xxzzzz_0, tg_xxyyyy_xxzzzz_1, tg_xxyyyy_xyyyy_1, tg_xxyyyy_xyyyyy_0, \
                                         tg_xxyyyy_xyyyyy_1, tg_xxyyyy_xyyyyz_0, tg_xxyyyy_xyyyyz_1, tg_xxyyyy_xyyyz_1, \
                                         tg_xxyyyy_xyyyzz_0, tg_xxyyyy_xyyyzz_1, tg_xxyyyy_xyyzz_1, tg_xxyyyy_xyyzzz_0, \
                                         tg_xxyyyy_xyyzzz_1, tg_xxyyyy_xyzzz_1, tg_xxyyyy_xyzzzz_0, tg_xxyyyy_xyzzzz_1, \
                                         tg_xxyyyy_xzzzz_1, tg_xxyyyy_xzzzzz_0, tg_xxyyyy_xzzzzz_1, tg_xxyyyy_yyyyy_1, \
                                         tg_xxyyyy_yyyyyy_0, tg_xxyyyy_yyyyyy_1, tg_xxyyyy_yyyyyz_0, tg_xxyyyy_yyyyyz_1, \
                                         tg_xxyyyy_yyyyz_1, tg_xxyyyy_yyyyzz_0, tg_xxyyyy_yyyyzz_1, tg_xxyyyy_yyyzz_1, \
                                         tg_xxyyyy_yyyzzz_0, tg_xxyyyy_yyyzzz_1, tg_xxyyyy_yyzzz_1, tg_xxyyyy_yyzzzz_0, \
                                         tg_xxyyyy_yyzzzz_1, tg_xxyyyy_yzzzz_1, tg_xxyyyy_yzzzzz_0, tg_xxyyyy_yzzzzz_1, \
                                         tg_xxyyyy_zzzzz_1, tg_xxyyyy_zzzzzz_0, tg_xxyyyy_zzzzzz_1, tg_xxyyyz_xxxxx_1, \
                                         tg_xxyyyz_xxxxxx_0, tg_xxyyyz_xxxxxx_1, tg_xxyyyz_xxxxxy_0, tg_xxyyyz_xxxxxy_1, \
                                         tg_xxyyyz_xxxxxz_0, tg_xxyyyz_xxxxxz_1, tg_xxyyyz_xxxxy_1, tg_xxyyyz_xxxxyy_0, \
                                         tg_xxyyyz_xxxxyy_1, tg_xxyyyz_xxxxyz_0, tg_xxyyyz_xxxxyz_1, tg_xxyyyz_xxxxz_1, \
                                         tg_xxyyyz_xxxxzz_0, tg_xxyyyz_xxxxzz_1, tg_xxyyyz_xxxyy_1, tg_xxyyyz_xxxyyy_0, \
                                         tg_xxyyyz_xxxyyy_1, tg_xxyyyz_xxxyyz_0, tg_xxyyyz_xxxyyz_1, tg_xxyyyz_xxxyz_1, \
                                         tg_xxyyyz_xxxyzz_0, tg_xxyyyz_xxxyzz_1, tg_xxyyyz_xxxzz_1, tg_xxyyyz_xxxzzz_0, \
                                         tg_xxyyyz_xxxzzz_1, tg_xxyyyz_xxyyy_1, tg_xxyyyz_xxyyyy_0, tg_xxyyyz_xxyyyy_1, \
                                         tg_xxyyyz_xxyyyz_0, tg_xxyyyz_xxyyyz_1, tg_xxyyyz_xxyyz_1, tg_xxyyyz_xxyyzz_0, \
                                         tg_xxyyyz_xxyyzz_1, tg_xxyyyz_xxyzz_1, tg_xxyyyz_xxyzzz_0, tg_xxyyyz_xxyzzz_1, \
                                         tg_xxyyyz_xxzzz_1, tg_xxyyyz_xxzzzz_0, tg_xxyyyz_xxzzzz_1, tg_xxyyyz_xyyyy_1, \
                                         tg_xxyyyz_xyyyyy_0, tg_xxyyyz_xyyyyy_1, tg_xxyyyz_xyyyyz_0, tg_xxyyyz_xyyyyz_1, \
                                         tg_xxyyyz_xyyyz_1, tg_xxyyyz_xyyyzz_0, tg_xxyyyz_xyyyzz_1, tg_xxyyyz_xyyzz_1, \
                                         tg_xxyyyz_xyyzzz_0, tg_xxyyyz_xyyzzz_1, tg_xxyyyz_xyzzz_1, tg_xxyyyz_xyzzzz_0, \
                                         tg_xxyyyz_xyzzzz_1, tg_xxyyyz_xzzzz_1, tg_xxyyyz_xzzzzz_0, tg_xxyyyz_xzzzzz_1, \
                                         tg_xxyyyz_yyyyy_1, tg_xxyyyz_yyyyyy_0, tg_xxyyyz_yyyyyy_1, tg_xxyyyz_yyyyyz_0, \
                                         tg_xxyyyz_yyyyyz_1, tg_xxyyyz_yyyyz_1, tg_xxyyyz_yyyyzz_0, tg_xxyyyz_yyyyzz_1, \
                                         tg_xxyyyz_yyyzz_1, tg_xxyyyz_yyyzzz_0, tg_xxyyyz_yyyzzz_1, tg_xxyyyz_yyzzz_1, \
                                         tg_xxyyyz_yyzzzz_0, tg_xxyyyz_yyzzzz_1, tg_xxyyyz_yzzzz_1, tg_xxyyyz_yzzzzz_0, \
                                         tg_xxyyyz_yzzzzz_1, tg_xxyyyz_zzzzz_1, tg_xxyyyz_zzzzzz_0, tg_xxyyyz_zzzzzz_1, \
                                         tg_xxyyzz_xxxxx_1, tg_xxyyzz_xxxxxx_0, tg_xxyyzz_xxxxxx_1, tg_xxyyzz_xxxxxy_0, \
                                         tg_xxyyzz_xxxxxy_1, tg_xxyyzz_xxxxxz_0, tg_xxyyzz_xxxxxz_1, tg_xxyyzz_xxxxy_1, \
                                         tg_xxyyzz_xxxxyy_0, tg_xxyyzz_xxxxyy_1, tg_xxyyzz_xxxxyz_0, tg_xxyyzz_xxxxyz_1, \
                                         tg_xxyyzz_xxxxz_1, tg_xxyyzz_xxxxzz_0, tg_xxyyzz_xxxxzz_1, tg_xxyyzz_xxxyy_1, \
                                         tg_xxyyzz_xxxyyy_0, tg_xxyyzz_xxxyyy_1, tg_xxyyzz_xxxyyz_0, tg_xxyyzz_xxxyyz_1, \
                                         tg_xxyyzz_xxxyz_1, tg_xxyyzz_xxxyzz_0, tg_xxyyzz_xxxyzz_1, tg_xxyyzz_xxxzz_1, \
                                         tg_xxyyzz_xxxzzz_0, tg_xxyyzz_xxxzzz_1, tg_xxyyzz_xxyyy_1, tg_xxyyzz_xxyyyy_0, \
                                         tg_xxyyzz_xxyyyy_1, tg_xxyyzz_xxyyyz_0, tg_xxyyzz_xxyyyz_1, tg_xxyyzz_xxyyz_1, \
                                         tg_xxyyzz_xxyyzz_0, tg_xxyyzz_xxyyzz_1, tg_xxyyzz_xxyzz_1, tg_xxyyzz_xxyzzz_0, \
                                         tg_xxyyzz_xxyzzz_1, tg_xxyyzz_xxzzz_1, tg_xxyyzz_xxzzzz_0, tg_xxyyzz_xxzzzz_1, \
                                         tg_xxyyzz_xyyyy_1, tg_xxyyzz_xyyyyy_0, tg_xxyyzz_xyyyyy_1, tg_xxyyzz_xyyyyz_0, \
                                         tg_xxyyzz_xyyyyz_1, tg_xxyyzz_xyyyz_1, tg_xxyyzz_xyyyzz_0, tg_xxyyzz_xyyyzz_1, \
                                         tg_xxyyzz_xyyzz_1, tg_xxyyzz_xyyzzz_0, tg_xxyyzz_xyyzzz_1, tg_xxyyzz_xyzzz_1, \
                                         tg_xxyyzz_xyzzzz_0, tg_xxyyzz_xyzzzz_1, tg_xxyyzz_xzzzz_1, tg_xxyyzz_xzzzzz_0, \
                                         tg_xxyyzz_xzzzzz_1, tg_xxyyzz_yyyyy_1, tg_xxyyzz_yyyyyy_0, tg_xxyyzz_yyyyyy_1, \
                                         tg_xxyyzz_yyyyyz_0, tg_xxyyzz_yyyyyz_1, tg_xxyyzz_yyyyz_1, tg_xxyyzz_yyyyzz_0, \
                                         tg_xxyyzz_yyyyzz_1, tg_xxyyzz_yyyzz_1, tg_xxyyzz_yyyzzz_0, tg_xxyyzz_yyyzzz_1, \
                                         tg_xxyyzz_yyzzz_1, tg_xxyyzz_yyzzzz_0, tg_xxyyzz_yyzzzz_1, tg_xxyyzz_yzzzz_1, \
                                         tg_xxyyzz_yzzzzz_0, tg_xxyyzz_yzzzzz_1, tg_xxyyzz_zzzzz_1, tg_xxyyzz_zzzzzz_0, \
                                         tg_xxyyzz_zzzzzz_1, tg_xxyzzz_xxxxx_1, tg_xxyzzz_xxxxxx_0, tg_xxyzzz_xxxxxx_1, \
                                         tg_xxyzzz_xxxxxy_0, tg_xxyzzz_xxxxxy_1, tg_xxyzzz_xxxxxz_0, tg_xxyzzz_xxxxxz_1, \
                                         tg_xxyzzz_xxxxy_1, tg_xxyzzz_xxxxyy_0, tg_xxyzzz_xxxxyy_1, tg_xxyzzz_xxxxz_1, \
                                         tg_xxyzzz_xxxyy_1, tg_xxzzz_yyyzzz_0, tg_xxzzz_yyyzzz_1, tg_xxzzz_yyzzzz_0, \
                                         tg_xxzzz_yyzzzz_1, tg_xxzzz_yzzzzz_0, tg_xxzzz_yzzzzz_1, tg_xxzzz_zzzzzz_0, \
                                         tg_xxzzz_zzzzzz_1, tg_xyyyy_xxxxxx_0, tg_xyyyy_xxxxxx_1, tg_xyyyy_xxxxxy_0, \
                                         tg_xyyyy_xxxxxy_1, tg_xyyyy_xxxxxz_0, tg_xyyyy_xxxxxz_1, tg_xyyyy_xxxxyy_0, \
                                         tg_xyyyy_xxxxyy_1, tg_xyyyy_xxxxyz_0, tg_xyyyy_xxxxyz_1, tg_xyyyy_xxxxzz_0, \
                                         tg_xyyyy_xxxxzz_1, tg_xyyyy_xxxyyy_0, tg_xyyyy_xxxyyy_1, tg_xyyyy_xxxyyz_0, \
                                         tg_xyyyy_xxxyyz_1, tg_xyyyy_xxxyzz_0, tg_xyyyy_xxxyzz_1, tg_xyyyy_xxxzzz_0, \
                                         tg_xyyyy_xxxzzz_1, tg_xyyyy_xxyyyy_0, tg_xyyyy_xxyyyy_1, tg_xyyyy_xxyyyz_0, \
                                         tg_xyyyy_xxyyyz_1, tg_xyyyy_xxyyzz_0, tg_xyyyy_xxyyzz_1, tg_xyyyy_xxyzzz_0, \
                                         tg_xyyyy_xxyzzz_1, tg_xyyyy_xxzzzz_0, tg_xyyyy_xxzzzz_1, tg_xyyyy_xyyyyy_0, \
                                         tg_xyyyy_xyyyyy_1, tg_xyyyy_xyyyyz_0, tg_xyyyy_xyyyyz_1, tg_xyyyy_xyyyzz_0, \
                                         tg_xyyyy_xyyyzz_1, tg_xyyyy_xyyzzz_0, tg_xyyyy_xyyzzz_1, tg_xyyyy_xyzzzz_0, \
                                         tg_xyyyy_xyzzzz_1, tg_xyyyy_xzzzzz_0, tg_xyyyy_xzzzzz_1, tg_xyyyy_yyyyyy_0, \
                                         tg_xyyyy_yyyyyy_1, tg_xyyyy_yyyyyz_0, tg_xyyyy_yyyyyz_1, tg_xyyyy_yyyyzz_0, \
                                         tg_xyyyy_yyyyzz_1, tg_xyyyy_yyyzzz_0, tg_xyyyy_yyyzzz_1, tg_xyyyy_yyzzzz_0, \
                                         tg_xyyyy_yyzzzz_1, tg_xyyyy_yzzzzz_0, tg_xyyyy_yzzzzz_1, tg_xyyyy_zzzzzz_0, \
                                         tg_xyyyy_zzzzzz_1, tg_xyyyz_xxxxxx_0, tg_xyyyz_xxxxxx_1, tg_xyyyz_xxxxxy_0, \
                                         tg_xyyyz_xxxxxy_1, tg_xyyyz_xxxxxz_0, tg_xyyyz_xxxxxz_1, tg_xyyyz_xxxxyy_0, \
                                         tg_xyyyz_xxxxyy_1, tg_xyyyz_xxxxyz_0, tg_xyyyz_xxxxyz_1, tg_xyyyz_xxxxzz_0, \
                                         tg_xyyyz_xxxxzz_1, tg_xyyyz_xxxyyy_0, tg_xyyyz_xxxyyy_1, tg_xyyyz_xxxyyz_0, \
                                         tg_xyyyz_xxxyyz_1, tg_xyyyz_xxxyzz_0, tg_xyyyz_xxxyzz_1, tg_xyyyz_xxxzzz_0, \
                                         tg_xyyyz_xxxzzz_1, tg_xyyyz_xxyyyy_0, tg_xyyyz_xxyyyy_1, tg_xyyyz_xxyyyz_0, \
                                         tg_xyyyz_xxyyyz_1, tg_xyyyz_xxyyzz_0, tg_xyyyz_xxyyzz_1, tg_xyyyz_xxyzzz_0, \
                                         tg_xyyyz_xxyzzz_1, tg_xyyyz_xxzzzz_0, tg_xyyyz_xxzzzz_1, tg_xyyyz_xyyyyy_0, \
                                         tg_xyyyz_xyyyyy_1, tg_xyyyz_xyyyyz_0, tg_xyyyz_xyyyyz_1, tg_xyyyz_xyyyzz_0, \
                                         tg_xyyyz_xyyyzz_1, tg_xyyyz_xyyzzz_0, tg_xyyyz_xyyzzz_1, tg_xyyyz_xyzzzz_0, \
                                         tg_xyyyz_xyzzzz_1, tg_xyyyz_xzzzzz_0, tg_xyyyz_xzzzzz_1, tg_xyyyz_yyyyyy_0, \
                                         tg_xyyyz_yyyyyy_1, tg_xyyyz_yyyyyz_0, tg_xyyyz_yyyyyz_1, tg_xyyyz_yyyyzz_0, \
                                         tg_xyyyz_yyyyzz_1, tg_xyyyz_yyyzzz_0, tg_xyyyz_yyyzzz_1, tg_xyyyz_yyzzzz_0, \
                                         tg_xyyyz_yyzzzz_1, tg_xyyyz_yzzzzz_0, tg_xyyyz_yzzzzz_1, tg_xyyyz_zzzzzz_0, \
                                         tg_xyyyz_zzzzzz_1, tg_xyyzz_xxxxxx_0, tg_xyyzz_xxxxxx_1, tg_xyyzz_xxxxxy_0, \
                                         tg_xyyzz_xxxxxy_1, tg_xyyzz_xxxxxz_0, tg_xyyzz_xxxxxz_1, tg_xyyzz_xxxxyy_0, \
                                         tg_xyyzz_xxxxyy_1, tg_xyyzz_xxxxyz_0, tg_xyyzz_xxxxyz_1, tg_xyyzz_xxxxzz_0, \
                                         tg_xyyzz_xxxxzz_1, tg_xyyzz_xxxyyy_0, tg_xyyzz_xxxyyy_1, tg_xyyzz_xxxyyz_0, \
                                         tg_xyyzz_xxxyyz_1, tg_xyyzz_xxxyzz_0, tg_xyyzz_xxxyzz_1, tg_xyyzz_xxxzzz_0, \
                                         tg_xyyzz_xxxzzz_1, tg_xyyzz_xxyyyy_0, tg_xyyzz_xxyyyy_1, tg_xyyzz_xxyyyz_0, \
                                         tg_xyyzz_xxyyyz_1, tg_xyyzz_xxyyzz_0, tg_xyyzz_xxyyzz_1, tg_xyyzz_xxyzzz_0, \
                                         tg_xyyzz_xxyzzz_1, tg_xyyzz_xxzzzz_0, tg_xyyzz_xxzzzz_1, tg_xyyzz_xyyyyy_0, \
                                         tg_xyyzz_xyyyyy_1, tg_xyyzz_xyyyyz_0, tg_xyyzz_xyyyyz_1, tg_xyyzz_xyyyzz_0, \
                                         tg_xyyzz_xyyyzz_1, tg_xyyzz_xyyzzz_0, tg_xyyzz_xyyzzz_1, tg_xyyzz_xyzzzz_0, \
                                         tg_xyyzz_xyzzzz_1, tg_xyyzz_xzzzzz_0, tg_xyyzz_xzzzzz_1, tg_xyyzz_yyyyyy_0, \
                                         tg_xyyzz_yyyyyy_1, tg_xyyzz_yyyyyz_0, tg_xyyzz_yyyyyz_1, tg_xyyzz_yyyyzz_0, \
                                         tg_xyyzz_yyyyzz_1, tg_xyyzz_yyyzzz_0, tg_xyyzz_yyyzzz_1, tg_xyyzz_yyzzzz_0, \
                                         tg_xyyzz_yyzzzz_1, tg_xyyzz_yzzzzz_0, tg_xyyzz_yzzzzz_1, tg_xyyzz_zzzzzz_0, \
                                         tg_xyyzz_zzzzzz_1, tg_xyzzz_xxxxxx_0, tg_xyzzz_xxxxxx_1, tg_xyzzz_xxxxxy_0, \
                                         tg_xyzzz_xxxxxy_1, tg_xyzzz_xxxxxz_0, tg_xyzzz_xxxxxz_1, tg_xyzzz_xxxxyy_0, \
                                         tg_xyzzz_xxxxyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxzzz_yyyzzz_0[j] = pb_x * tg_xxxzzz_yyyzzz_0[j] + fr * tg_xxxzzz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyzzz_0[j] - tg_xxzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_yyzzzz_0[j] = pb_x * tg_xxxzzz_yyzzzz_0[j] + fr * tg_xxxzzz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyzzzz_0[j] - tg_xxzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_yzzzzz_0[j] = pb_x * tg_xxxzzz_yzzzzz_0[j] + fr * tg_xxxzzz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yzzzzz_0[j] - tg_xxzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_zzzzzz_0[j] = pb_x * tg_xxxzzz_zzzzzz_0[j] + fr * tg_xxxzzz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_zzzzzz_0[j] - tg_xxzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_xxxxxx_0[j] = pb_x * tg_xxyyyy_xxxxxx_0[j] + fr * tg_xxyyyy_xxxxxx_1[j] + fl1_fx * (tg_xyyyy_xxxxxx_0[j] - tg_xyyyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyyyy_xxxxx_1[j];

                    tg_xxxyyyy_xxxxxy_0[j] = pb_x * tg_xxyyyy_xxxxxy_0[j] + fr * tg_xxyyyy_xxxxxy_1[j] + fl1_fx * (tg_xyyyy_xxxxxy_0[j] - tg_xyyyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyy_xxxxy_1[j];

                    tg_xxxyyyy_xxxxxz_0[j] = pb_x * tg_xxyyyy_xxxxxz_0[j] + fr * tg_xxyyyy_xxxxxz_1[j] + fl1_fx * (tg_xyyyy_xxxxxz_0[j] - tg_xyyyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyy_xxxxz_1[j];

                    tg_xxxyyyy_xxxxyy_0[j] = pb_x * tg_xxyyyy_xxxxyy_0[j] + fr * tg_xxyyyy_xxxxyy_1[j] + fl1_fx * (tg_xyyyy_xxxxyy_0[j] - tg_xyyyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyy_xxxyy_1[j];

                    tg_xxxyyyy_xxxxyz_0[j] = pb_x * tg_xxyyyy_xxxxyz_0[j] + fr * tg_xxyyyy_xxxxyz_1[j] + fl1_fx * (tg_xyyyy_xxxxyz_0[j] - tg_xyyyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyy_xxxyz_1[j];

                    tg_xxxyyyy_xxxxzz_0[j] = pb_x * tg_xxyyyy_xxxxzz_0[j] + fr * tg_xxyyyy_xxxxzz_1[j] + fl1_fx * (tg_xyyyy_xxxxzz_0[j] - tg_xyyyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyy_xxxzz_1[j];

                    tg_xxxyyyy_xxxyyy_0[j] = pb_x * tg_xxyyyy_xxxyyy_0[j] + fr * tg_xxyyyy_xxxyyy_1[j] + fl1_fx * (tg_xyyyy_xxxyyy_0[j] - tg_xyyyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxyyy_1[j];

                    tg_xxxyyyy_xxxyyz_0[j] = pb_x * tg_xxyyyy_xxxyyz_0[j] + fr * tg_xxyyyy_xxxyyz_1[j] + fl1_fx * (tg_xyyyy_xxxyyz_0[j] - tg_xyyyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxyyz_1[j];

                    tg_xxxyyyy_xxxyzz_0[j] = pb_x * tg_xxyyyy_xxxyzz_0[j] + fr * tg_xxyyyy_xxxyzz_1[j] + fl1_fx * (tg_xyyyy_xxxyzz_0[j] - tg_xyyyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxyzz_1[j];

                    tg_xxxyyyy_xxxzzz_0[j] = pb_x * tg_xxyyyy_xxxzzz_0[j] + fr * tg_xxyyyy_xxxzzz_1[j] + fl1_fx * (tg_xyyyy_xxxzzz_0[j] - tg_xyyyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxzzz_1[j];

                    tg_xxxyyyy_xxyyyy_0[j] = pb_x * tg_xxyyyy_xxyyyy_0[j] + fr * tg_xxyyyy_xxyyyy_1[j] + fl1_fx * (tg_xyyyy_xxyyyy_0[j] - tg_xyyyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyyyy_1[j];

                    tg_xxxyyyy_xxyyyz_0[j] = pb_x * tg_xxyyyy_xxyyyz_0[j] + fr * tg_xxyyyy_xxyyyz_1[j] + fl1_fx * (tg_xyyyy_xxyyyz_0[j] - tg_xyyyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyyyz_1[j];

                    tg_xxxyyyy_xxyyzz_0[j] = pb_x * tg_xxyyyy_xxyyzz_0[j] + fr * tg_xxyyyy_xxyyzz_1[j] + fl1_fx * (tg_xyyyy_xxyyzz_0[j] - tg_xyyyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyyzz_1[j];

                    tg_xxxyyyy_xxyzzz_0[j] = pb_x * tg_xxyyyy_xxyzzz_0[j] + fr * tg_xxyyyy_xxyzzz_1[j] + fl1_fx * (tg_xyyyy_xxyzzz_0[j] - tg_xyyyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyzzz_1[j];

                    tg_xxxyyyy_xxzzzz_0[j] = pb_x * tg_xxyyyy_xxzzzz_0[j] + fr * tg_xxyyyy_xxzzzz_1[j] + fl1_fx * (tg_xyyyy_xxzzzz_0[j] - tg_xyyyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xzzzz_1[j];

                    tg_xxxyyyy_xyyyyy_0[j] = pb_x * tg_xxyyyy_xyyyyy_0[j] + fr * tg_xxyyyy_xyyyyy_1[j] + fl1_fx * (tg_xyyyy_xyyyyy_0[j] - tg_xyyyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyyyy_1[j];

                    tg_xxxyyyy_xyyyyz_0[j] = pb_x * tg_xxyyyy_xyyyyz_0[j] + fr * tg_xxyyyy_xyyyyz_1[j] + fl1_fx * (tg_xyyyy_xyyyyz_0[j] - tg_xyyyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyyyz_1[j];

                    tg_xxxyyyy_xyyyzz_0[j] = pb_x * tg_xxyyyy_xyyyzz_0[j] + fr * tg_xxyyyy_xyyyzz_1[j] + fl1_fx * (tg_xyyyy_xyyyzz_0[j] - tg_xyyyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyyzz_1[j];

                    tg_xxxyyyy_xyyzzz_0[j] = pb_x * tg_xxyyyy_xyyzzz_0[j] + fr * tg_xxyyyy_xyyzzz_1[j] + fl1_fx * (tg_xyyyy_xyyzzz_0[j] - tg_xyyyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyzzz_1[j];

                    tg_xxxyyyy_xyzzzz_0[j] = pb_x * tg_xxyyyy_xyzzzz_0[j] + fr * tg_xxyyyy_xyzzzz_1[j] + fl1_fx * (tg_xyyyy_xyzzzz_0[j] - tg_xyyyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yzzzz_1[j];

                    tg_xxxyyyy_xzzzzz_0[j] = pb_x * tg_xxyyyy_xzzzzz_0[j] + fr * tg_xxyyyy_xzzzzz_1[j] + fl1_fx * (tg_xyyyy_xzzzzz_0[j] - tg_xyyyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_zzzzz_1[j];

                    tg_xxxyyyy_yyyyyy_0[j] = pb_x * tg_xxyyyy_yyyyyy_0[j] + fr * tg_xxyyyy_yyyyyy_1[j] + fl1_fx * (tg_xyyyy_yyyyyy_0[j] - tg_xyyyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyyyy_yyyyyz_0[j] = pb_x * tg_xxyyyy_yyyyyz_0[j] + fr * tg_xxyyyy_yyyyyz_1[j] + fl1_fx * (tg_xyyyy_yyyyyz_0[j] - tg_xyyyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyyyy_yyyyzz_0[j] = pb_x * tg_xxyyyy_yyyyzz_0[j] + fr * tg_xxyyyy_yyyyzz_1[j] + fl1_fx * (tg_xyyyy_yyyyzz_0[j] - tg_xyyyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyyyy_yyyzzz_0[j] = pb_x * tg_xxyyyy_yyyzzz_0[j] + fr * tg_xxyyyy_yyyzzz_1[j] + fl1_fx * (tg_xyyyy_yyyzzz_0[j] - tg_xyyyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_yyzzzz_0[j] = pb_x * tg_xxyyyy_yyzzzz_0[j] + fr * tg_xxyyyy_yyzzzz_1[j] + fl1_fx * (tg_xyyyy_yyzzzz_0[j] - tg_xyyyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_yzzzzz_0[j] = pb_x * tg_xxyyyy_yzzzzz_0[j] + fr * tg_xxyyyy_yzzzzz_1[j] + fl1_fx * (tg_xyyyy_yzzzzz_0[j] - tg_xyyyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_zzzzzz_0[j] = pb_x * tg_xxyyyy_zzzzzz_0[j] + fr * tg_xxyyyy_zzzzzz_1[j] + fl1_fx * (tg_xyyyy_zzzzzz_0[j] - tg_xyyyy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_xxxxxx_0[j] = pb_x * tg_xxyyyz_xxxxxx_0[j] + fr * tg_xxyyyz_xxxxxx_1[j] + fl1_fx * (tg_xyyyz_xxxxxx_0[j] - tg_xyyyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyyyz_xxxxx_1[j];

                    tg_xxxyyyz_xxxxxy_0[j] = pb_x * tg_xxyyyz_xxxxxy_0[j] + fr * tg_xxyyyz_xxxxxy_1[j] + fl1_fx * (tg_xyyyz_xxxxxy_0[j] - tg_xyyyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyz_xxxxy_1[j];

                    tg_xxxyyyz_xxxxxz_0[j] = pb_x * tg_xxyyyz_xxxxxz_0[j] + fr * tg_xxyyyz_xxxxxz_1[j] + fl1_fx * (tg_xyyyz_xxxxxz_0[j] - tg_xyyyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyyz_xxxxz_1[j];

                    tg_xxxyyyz_xxxxyy_0[j] = pb_x * tg_xxyyyz_xxxxyy_0[j] + fr * tg_xxyyyz_xxxxyy_1[j] + fl1_fx * (tg_xyyyz_xxxxyy_0[j] - tg_xyyyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyz_xxxyy_1[j];

                    tg_xxxyyyz_xxxxyz_0[j] = pb_x * tg_xxyyyz_xxxxyz_0[j] + fr * tg_xxyyyz_xxxxyz_1[j] + fl1_fx * (tg_xyyyz_xxxxyz_0[j] - tg_xyyyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyz_xxxyz_1[j];

                    tg_xxxyyyz_xxxxzz_0[j] = pb_x * tg_xxyyyz_xxxxzz_0[j] + fr * tg_xxyyyz_xxxxzz_1[j] + fl1_fx * (tg_xyyyz_xxxxzz_0[j] - tg_xyyyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyz_xxxzz_1[j];

                    tg_xxxyyyz_xxxyyy_0[j] = pb_x * tg_xxyyyz_xxxyyy_0[j] + fr * tg_xxyyyz_xxxyyy_1[j] + fl1_fx * (tg_xyyyz_xxxyyy_0[j] - tg_xyyyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxyyy_1[j];

                    tg_xxxyyyz_xxxyyz_0[j] = pb_x * tg_xxyyyz_xxxyyz_0[j] + fr * tg_xxyyyz_xxxyyz_1[j] + fl1_fx * (tg_xyyyz_xxxyyz_0[j] - tg_xyyyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxyyz_1[j];

                    tg_xxxyyyz_xxxyzz_0[j] = pb_x * tg_xxyyyz_xxxyzz_0[j] + fr * tg_xxyyyz_xxxyzz_1[j] + fl1_fx * (tg_xyyyz_xxxyzz_0[j] - tg_xyyyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxyzz_1[j];

                    tg_xxxyyyz_xxxzzz_0[j] = pb_x * tg_xxyyyz_xxxzzz_0[j] + fr * tg_xxyyyz_xxxzzz_1[j] + fl1_fx * (tg_xyyyz_xxxzzz_0[j] - tg_xyyyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxzzz_1[j];

                    tg_xxxyyyz_xxyyyy_0[j] = pb_x * tg_xxyyyz_xxyyyy_0[j] + fr * tg_xxyyyz_xxyyyy_1[j] + fl1_fx * (tg_xyyyz_xxyyyy_0[j] - tg_xyyyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyyyy_1[j];

                    tg_xxxyyyz_xxyyyz_0[j] = pb_x * tg_xxyyyz_xxyyyz_0[j] + fr * tg_xxyyyz_xxyyyz_1[j] + fl1_fx * (tg_xyyyz_xxyyyz_0[j] - tg_xyyyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyyyz_1[j];

                    tg_xxxyyyz_xxyyzz_0[j] = pb_x * tg_xxyyyz_xxyyzz_0[j] + fr * tg_xxyyyz_xxyyzz_1[j] + fl1_fx * (tg_xyyyz_xxyyzz_0[j] - tg_xyyyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyyzz_1[j];

                    tg_xxxyyyz_xxyzzz_0[j] = pb_x * tg_xxyyyz_xxyzzz_0[j] + fr * tg_xxyyyz_xxyzzz_1[j] + fl1_fx * (tg_xyyyz_xxyzzz_0[j] - tg_xyyyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyzzz_1[j];

                    tg_xxxyyyz_xxzzzz_0[j] = pb_x * tg_xxyyyz_xxzzzz_0[j] + fr * tg_xxyyyz_xxzzzz_1[j] + fl1_fx * (tg_xyyyz_xxzzzz_0[j] - tg_xyyyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xzzzz_1[j];

                    tg_xxxyyyz_xyyyyy_0[j] = pb_x * tg_xxyyyz_xyyyyy_0[j] + fr * tg_xxyyyz_xyyyyy_1[j] + fl1_fx * (tg_xyyyz_xyyyyy_0[j] - tg_xyyyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyyyy_1[j];

                    tg_xxxyyyz_xyyyyz_0[j] = pb_x * tg_xxyyyz_xyyyyz_0[j] + fr * tg_xxyyyz_xyyyyz_1[j] + fl1_fx * (tg_xyyyz_xyyyyz_0[j] - tg_xyyyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyyyz_1[j];

                    tg_xxxyyyz_xyyyzz_0[j] = pb_x * tg_xxyyyz_xyyyzz_0[j] + fr * tg_xxyyyz_xyyyzz_1[j] + fl1_fx * (tg_xyyyz_xyyyzz_0[j] - tg_xyyyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyyzz_1[j];

                    tg_xxxyyyz_xyyzzz_0[j] = pb_x * tg_xxyyyz_xyyzzz_0[j] + fr * tg_xxyyyz_xyyzzz_1[j] + fl1_fx * (tg_xyyyz_xyyzzz_0[j] - tg_xyyyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyzzz_1[j];

                    tg_xxxyyyz_xyzzzz_0[j] = pb_x * tg_xxyyyz_xyzzzz_0[j] + fr * tg_xxyyyz_xyzzzz_1[j] + fl1_fx * (tg_xyyyz_xyzzzz_0[j] - tg_xyyyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yzzzz_1[j];

                    tg_xxxyyyz_xzzzzz_0[j] = pb_x * tg_xxyyyz_xzzzzz_0[j] + fr * tg_xxyyyz_xzzzzz_1[j] + fl1_fx * (tg_xyyyz_xzzzzz_0[j] - tg_xyyyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_zzzzz_1[j];

                    tg_xxxyyyz_yyyyyy_0[j] = pb_x * tg_xxyyyz_yyyyyy_0[j] + fr * tg_xxyyyz_yyyyyy_1[j] + fl1_fx * (tg_xyyyz_yyyyyy_0[j] - tg_xyyyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyyyz_yyyyyz_0[j] = pb_x * tg_xxyyyz_yyyyyz_0[j] + fr * tg_xxyyyz_yyyyyz_1[j] + fl1_fx * (tg_xyyyz_yyyyyz_0[j] - tg_xyyyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyyyz_yyyyzz_0[j] = pb_x * tg_xxyyyz_yyyyzz_0[j] + fr * tg_xxyyyz_yyyyzz_1[j] + fl1_fx * (tg_xyyyz_yyyyzz_0[j] - tg_xyyyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyyyz_yyyzzz_0[j] = pb_x * tg_xxyyyz_yyyzzz_0[j] + fr * tg_xxyyyz_yyyzzz_1[j] + fl1_fx * (tg_xyyyz_yyyzzz_0[j] - tg_xyyyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_yyzzzz_0[j] = pb_x * tg_xxyyyz_yyzzzz_0[j] + fr * tg_xxyyyz_yyzzzz_1[j] + fl1_fx * (tg_xyyyz_yyzzzz_0[j] - tg_xyyyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_yzzzzz_0[j] = pb_x * tg_xxyyyz_yzzzzz_0[j] + fr * tg_xxyyyz_yzzzzz_1[j] + fl1_fx * (tg_xyyyz_yzzzzz_0[j] - tg_xyyyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_zzzzzz_0[j] = pb_x * tg_xxyyyz_zzzzzz_0[j] + fr * tg_xxyyyz_zzzzzz_1[j] + fl1_fx * (tg_xyyyz_zzzzzz_0[j] - tg_xyyyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_xxxxxx_0[j] = pb_x * tg_xxyyzz_xxxxxx_0[j] + fr * tg_xxyyzz_xxxxxx_1[j] + fl1_fx * (tg_xyyzz_xxxxxx_0[j] - tg_xyyzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyyzz_xxxxx_1[j];

                    tg_xxxyyzz_xxxxxy_0[j] = pb_x * tg_xxyyzz_xxxxxy_0[j] + fr * tg_xxyyzz_xxxxxy_1[j] + fl1_fx * (tg_xyyzz_xxxxxy_0[j] - tg_xyyzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyzz_xxxxy_1[j];

                    tg_xxxyyzz_xxxxxz_0[j] = pb_x * tg_xxyyzz_xxxxxz_0[j] + fr * tg_xxyyzz_xxxxxz_1[j] + fl1_fx * (tg_xyyzz_xxxxxz_0[j] - tg_xyyzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyyzz_xxxxz_1[j];

                    tg_xxxyyzz_xxxxyy_0[j] = pb_x * tg_xxyyzz_xxxxyy_0[j] + fr * tg_xxyyzz_xxxxyy_1[j] + fl1_fx * (tg_xyyzz_xxxxyy_0[j] - tg_xyyzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzz_xxxyy_1[j];

                    tg_xxxyyzz_xxxxyz_0[j] = pb_x * tg_xxyyzz_xxxxyz_0[j] + fr * tg_xxyyzz_xxxxyz_1[j] + fl1_fx * (tg_xyyzz_xxxxyz_0[j] - tg_xyyzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzz_xxxyz_1[j];

                    tg_xxxyyzz_xxxxzz_0[j] = pb_x * tg_xxyyzz_xxxxzz_0[j] + fr * tg_xxyyzz_xxxxzz_1[j] + fl1_fx * (tg_xyyzz_xxxxzz_0[j] - tg_xyyzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzz_xxxzz_1[j];

                    tg_xxxyyzz_xxxyyy_0[j] = pb_x * tg_xxyyzz_xxxyyy_0[j] + fr * tg_xxyyzz_xxxyyy_1[j] + fl1_fx * (tg_xyyzz_xxxyyy_0[j] - tg_xyyzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxyyy_1[j];

                    tg_xxxyyzz_xxxyyz_0[j] = pb_x * tg_xxyyzz_xxxyyz_0[j] + fr * tg_xxyyzz_xxxyyz_1[j] + fl1_fx * (tg_xyyzz_xxxyyz_0[j] - tg_xyyzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxyyz_1[j];

                    tg_xxxyyzz_xxxyzz_0[j] = pb_x * tg_xxyyzz_xxxyzz_0[j] + fr * tg_xxyyzz_xxxyzz_1[j] + fl1_fx * (tg_xyyzz_xxxyzz_0[j] - tg_xyyzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxyzz_1[j];

                    tg_xxxyyzz_xxxzzz_0[j] = pb_x * tg_xxyyzz_xxxzzz_0[j] + fr * tg_xxyyzz_xxxzzz_1[j] + fl1_fx * (tg_xyyzz_xxxzzz_0[j] - tg_xyyzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxzzz_1[j];

                    tg_xxxyyzz_xxyyyy_0[j] = pb_x * tg_xxyyzz_xxyyyy_0[j] + fr * tg_xxyyzz_xxyyyy_1[j] + fl1_fx * (tg_xyyzz_xxyyyy_0[j] - tg_xyyzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyyyy_1[j];

                    tg_xxxyyzz_xxyyyz_0[j] = pb_x * tg_xxyyzz_xxyyyz_0[j] + fr * tg_xxyyzz_xxyyyz_1[j] + fl1_fx * (tg_xyyzz_xxyyyz_0[j] - tg_xyyzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyyyz_1[j];

                    tg_xxxyyzz_xxyyzz_0[j] = pb_x * tg_xxyyzz_xxyyzz_0[j] + fr * tg_xxyyzz_xxyyzz_1[j] + fl1_fx * (tg_xyyzz_xxyyzz_0[j] - tg_xyyzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyyzz_1[j];

                    tg_xxxyyzz_xxyzzz_0[j] = pb_x * tg_xxyyzz_xxyzzz_0[j] + fr * tg_xxyyzz_xxyzzz_1[j] + fl1_fx * (tg_xyyzz_xxyzzz_0[j] - tg_xyyzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyzzz_1[j];

                    tg_xxxyyzz_xxzzzz_0[j] = pb_x * tg_xxyyzz_xxzzzz_0[j] + fr * tg_xxyyzz_xxzzzz_1[j] + fl1_fx * (tg_xyyzz_xxzzzz_0[j] - tg_xyyzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xzzzz_1[j];

                    tg_xxxyyzz_xyyyyy_0[j] = pb_x * tg_xxyyzz_xyyyyy_0[j] + fr * tg_xxyyzz_xyyyyy_1[j] + fl1_fx * (tg_xyyzz_xyyyyy_0[j] - tg_xyyzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyyyy_1[j];

                    tg_xxxyyzz_xyyyyz_0[j] = pb_x * tg_xxyyzz_xyyyyz_0[j] + fr * tg_xxyyzz_xyyyyz_1[j] + fl1_fx * (tg_xyyzz_xyyyyz_0[j] - tg_xyyzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyyyz_1[j];

                    tg_xxxyyzz_xyyyzz_0[j] = pb_x * tg_xxyyzz_xyyyzz_0[j] + fr * tg_xxyyzz_xyyyzz_1[j] + fl1_fx * (tg_xyyzz_xyyyzz_0[j] - tg_xyyzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyyzz_1[j];

                    tg_xxxyyzz_xyyzzz_0[j] = pb_x * tg_xxyyzz_xyyzzz_0[j] + fr * tg_xxyyzz_xyyzzz_1[j] + fl1_fx * (tg_xyyzz_xyyzzz_0[j] - tg_xyyzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyzzz_1[j];

                    tg_xxxyyzz_xyzzzz_0[j] = pb_x * tg_xxyyzz_xyzzzz_0[j] + fr * tg_xxyyzz_xyzzzz_1[j] + fl1_fx * (tg_xyyzz_xyzzzz_0[j] - tg_xyyzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yzzzz_1[j];

                    tg_xxxyyzz_xzzzzz_0[j] = pb_x * tg_xxyyzz_xzzzzz_0[j] + fr * tg_xxyyzz_xzzzzz_1[j] + fl1_fx * (tg_xyyzz_xzzzzz_0[j] - tg_xyyzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_zzzzz_1[j];

                    tg_xxxyyzz_yyyyyy_0[j] = pb_x * tg_xxyyzz_yyyyyy_0[j] + fr * tg_xxyyzz_yyyyyy_1[j] + fl1_fx * (tg_xyyzz_yyyyyy_0[j] - tg_xyyzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyyzz_yyyyyz_0[j] = pb_x * tg_xxyyzz_yyyyyz_0[j] + fr * tg_xxyyzz_yyyyyz_1[j] + fl1_fx * (tg_xyyzz_yyyyyz_0[j] - tg_xyyzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyyzz_yyyyzz_0[j] = pb_x * tg_xxyyzz_yyyyzz_0[j] + fr * tg_xxyyzz_yyyyzz_1[j] + fl1_fx * (tg_xyyzz_yyyyzz_0[j] - tg_xyyzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyyzz_yyyzzz_0[j] = pb_x * tg_xxyyzz_yyyzzz_0[j] + fr * tg_xxyyzz_yyyzzz_1[j] + fl1_fx * (tg_xyyzz_yyyzzz_0[j] - tg_xyyzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_yyzzzz_0[j] = pb_x * tg_xxyyzz_yyzzzz_0[j] + fr * tg_xxyyzz_yyzzzz_1[j] + fl1_fx * (tg_xyyzz_yyzzzz_0[j] - tg_xyyzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_yzzzzz_0[j] = pb_x * tg_xxyyzz_yzzzzz_0[j] + fr * tg_xxyyzz_yzzzzz_1[j] + fl1_fx * (tg_xyyzz_yzzzzz_0[j] - tg_xyyzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_zzzzzz_0[j] = pb_x * tg_xxyyzz_zzzzzz_0[j] + fr * tg_xxyyzz_zzzzzz_1[j] + fl1_fx * (tg_xyyzz_zzzzzz_0[j] - tg_xyyzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_xxxxxx_0[j] = pb_x * tg_xxyzzz_xxxxxx_0[j] + fr * tg_xxyzzz_xxxxxx_1[j] + fl1_fx * (tg_xyzzz_xxxxxx_0[j] - tg_xyzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyzzz_xxxxx_1[j];

                    tg_xxxyzzz_xxxxxy_0[j] = pb_x * tg_xxyzzz_xxxxxy_0[j] + fr * tg_xxyzzz_xxxxxy_1[j] + fl1_fx * (tg_xyzzz_xxxxxy_0[j] - tg_xyzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyzzz_xxxxy_1[j];

                    tg_xxxyzzz_xxxxxz_0[j] = pb_x * tg_xxyzzz_xxxxxz_0[j] + fr * tg_xxyzzz_xxxxxz_1[j] + fl1_fx * (tg_xyzzz_xxxxxz_0[j] - tg_xyzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyzzz_xxxxz_1[j];

                    tg_xxxyzzz_xxxxyy_0[j] = pb_x * tg_xxyzzz_xxxxyy_0[j] + fr * tg_xxyzzz_xxxxyy_1[j] + fl1_fx * (tg_xyzzz_xxxxyy_0[j] - tg_xyzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzz_xxxyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_368_460(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (368,460)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxyzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 368); 

                auto tg_xxyzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 369); 

                auto tg_xxyzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 370); 

                auto tg_xxyzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 371); 

                auto tg_xxyzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 372); 

                auto tg_xxyzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 373); 

                auto tg_xxyzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 374); 

                auto tg_xxyzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 375); 

                auto tg_xxyzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 376); 

                auto tg_xxyzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 377); 

                auto tg_xxyzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 378); 

                auto tg_xxyzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 379); 

                auto tg_xxyzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 380); 

                auto tg_xxyzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 381); 

                auto tg_xxyzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 382); 

                auto tg_xxyzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 383); 

                auto tg_xxyzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 384); 

                auto tg_xxyzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 385); 

                auto tg_xxyzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 386); 

                auto tg_xxyzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 387); 

                auto tg_xxyzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 388); 

                auto tg_xxyzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 389); 

                auto tg_xxyzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 390); 

                auto tg_xxyzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 391); 

                auto tg_xxzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 392); 

                auto tg_xxzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 393); 

                auto tg_xxzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 394); 

                auto tg_xxzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 395); 

                auto tg_xxzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 396); 

                auto tg_xxzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 397); 

                auto tg_xxzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 398); 

                auto tg_xxzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 399); 

                auto tg_xxzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 400); 

                auto tg_xxzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 401); 

                auto tg_xxzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 402); 

                auto tg_xxzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 403); 

                auto tg_xxzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 404); 

                auto tg_xxzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 405); 

                auto tg_xxzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 406); 

                auto tg_xxzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 407); 

                auto tg_xxzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 408); 

                auto tg_xxzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 409); 

                auto tg_xxzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 410); 

                auto tg_xxzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 411); 

                auto tg_xxzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 412); 

                auto tg_xxzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 413); 

                auto tg_xxzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 414); 

                auto tg_xxzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 415); 

                auto tg_xxzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 416); 

                auto tg_xxzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 417); 

                auto tg_xxzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 418); 

                auto tg_xxzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 419); 

                auto tg_xyyyyy_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 420); 

                auto tg_xyyyyy_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 421); 

                auto tg_xyyyyy_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 422); 

                auto tg_xyyyyy_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 423); 

                auto tg_xyyyyy_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 424); 

                auto tg_xyyyyy_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 425); 

                auto tg_xyyyyy_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 426); 

                auto tg_xyyyyy_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 427); 

                auto tg_xyyyyy_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 428); 

                auto tg_xyyyyy_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 429); 

                auto tg_xyyyyy_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 430); 

                auto tg_xyyyyy_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 431); 

                auto tg_xyyyyy_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 432); 

                auto tg_xyyyyy_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 433); 

                auto tg_xyyyyy_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 434); 

                auto tg_xyyyyy_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 435); 

                auto tg_xyyyyy_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 436); 

                auto tg_xyyyyy_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 437); 

                auto tg_xyyyyy_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 438); 

                auto tg_xyyyyy_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 439); 

                auto tg_xyyyyy_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 440); 

                auto tg_xyyyyy_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 441); 

                auto tg_xyyyyy_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 442); 

                auto tg_xyyyyy_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 443); 

                auto tg_xyyyyy_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 444); 

                auto tg_xyyyyy_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 445); 

                auto tg_xyyyyy_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 446); 

                auto tg_xyyyyy_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 447); 

                auto tg_xyyyyz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 448); 

                auto tg_xyyyyz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 449); 

                auto tg_xyyyyz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 450); 

                auto tg_xyyyyz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 451); 

                auto tg_xyyyyz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 452); 

                auto tg_xyyyyz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 453); 

                auto tg_xyyyyz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 454); 

                auto tg_xyyyyz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 455); 

                auto tg_xyyyyz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 456); 

                auto tg_xyyyyz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 457); 

                auto tg_xyyyyz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 458); 

                auto tg_xyyyyz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 459); 

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

                auto tg_xxyzzz_xxxyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 277); 

                auto tg_xxyzzz_xxxzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 278); 

                auto tg_xxyzzz_xxyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 279); 

                auto tg_xxyzzz_xxyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 280); 

                auto tg_xxyzzz_xxyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 281); 

                auto tg_xxyzzz_xxzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 282); 

                auto tg_xxyzzz_xyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 283); 

                auto tg_xxyzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 284); 

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

                // set up pointers to integrals

                auto tg_xxxyzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 368); 

                auto tg_xxxyzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 369); 

                auto tg_xxxyzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 370); 

                auto tg_xxxyzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 371); 

                auto tg_xxxyzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 372); 

                auto tg_xxxyzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 373); 

                auto tg_xxxyzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 374); 

                auto tg_xxxyzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 375); 

                auto tg_xxxyzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 376); 

                auto tg_xxxyzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 377); 

                auto tg_xxxyzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 378); 

                auto tg_xxxyzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 379); 

                auto tg_xxxyzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 380); 

                auto tg_xxxyzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 381); 

                auto tg_xxxyzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 382); 

                auto tg_xxxyzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 383); 

                auto tg_xxxyzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 384); 

                auto tg_xxxyzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 385); 

                auto tg_xxxyzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 386); 

                auto tg_xxxyzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 387); 

                auto tg_xxxyzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 388); 

                auto tg_xxxyzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 389); 

                auto tg_xxxyzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 390); 

                auto tg_xxxyzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 391); 

                auto tg_xxxzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 392); 

                auto tg_xxxzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 393); 

                auto tg_xxxzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 394); 

                auto tg_xxxzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 395); 

                auto tg_xxxzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 396); 

                auto tg_xxxzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 397); 

                auto tg_xxxzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 398); 

                auto tg_xxxzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 399); 

                auto tg_xxxzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 400); 

                auto tg_xxxzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 401); 

                auto tg_xxxzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 402); 

                auto tg_xxxzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 403); 

                auto tg_xxxzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 404); 

                auto tg_xxxzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 405); 

                auto tg_xxxzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 406); 

                auto tg_xxxzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 407); 

                auto tg_xxxzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 408); 

                auto tg_xxxzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 409); 

                auto tg_xxxzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 410); 

                auto tg_xxxzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 411); 

                auto tg_xxxzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 412); 

                auto tg_xxxzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 413); 

                auto tg_xxxzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 414); 

                auto tg_xxxzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 415); 

                auto tg_xxxzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 416); 

                auto tg_xxxzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 417); 

                auto tg_xxxzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 418); 

                auto tg_xxxzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 419); 

                auto tg_xxyyyyy_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 420); 

                auto tg_xxyyyyy_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 421); 

                auto tg_xxyyyyy_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 422); 

                auto tg_xxyyyyy_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 423); 

                auto tg_xxyyyyy_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 424); 

                auto tg_xxyyyyy_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 425); 

                auto tg_xxyyyyy_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 426); 

                auto tg_xxyyyyy_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 427); 

                auto tg_xxyyyyy_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 428); 

                auto tg_xxyyyyy_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 429); 

                auto tg_xxyyyyy_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 430); 

                auto tg_xxyyyyy_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 431); 

                auto tg_xxyyyyy_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 432); 

                auto tg_xxyyyyy_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 433); 

                auto tg_xxyyyyy_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 434); 

                auto tg_xxyyyyy_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 435); 

                auto tg_xxyyyyy_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 436); 

                auto tg_xxyyyyy_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 437); 

                auto tg_xxyyyyy_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 438); 

                auto tg_xxyyyyy_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 439); 

                auto tg_xxyyyyy_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 440); 

                auto tg_xxyyyyy_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 441); 

                auto tg_xxyyyyy_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 442); 

                auto tg_xxyyyyy_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 443); 

                auto tg_xxyyyyy_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 444); 

                auto tg_xxyyyyy_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 445); 

                auto tg_xxyyyyy_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 446); 

                auto tg_xxyyyyy_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 447); 

                auto tg_xxyyyyz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 448); 

                auto tg_xxyyyyz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 449); 

                auto tg_xxyyyyz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 450); 

                auto tg_xxyyyyz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 451); 

                auto tg_xxyyyyz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 452); 

                auto tg_xxyyyyz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 453); 

                auto tg_xxyyyyz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 454); 

                auto tg_xxyyyyz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 455); 

                auto tg_xxyyyyz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 456); 

                auto tg_xxyyyyz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 457); 

                auto tg_xxyyyyz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 458); 

                auto tg_xxyyyyz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 459); 

                // Batch of Integrals (368,460)

                #pragma omp simd aligned(fxn, fza, tg_xxxyzzz_xxxxyz_0, tg_xxxyzzz_xxxxzz_0, \
                                         tg_xxxyzzz_xxxyyy_0, tg_xxxyzzz_xxxyyz_0, tg_xxxyzzz_xxxyzz_0, tg_xxxyzzz_xxxzzz_0, \
                                         tg_xxxyzzz_xxyyyy_0, tg_xxxyzzz_xxyyyz_0, tg_xxxyzzz_xxyyzz_0, tg_xxxyzzz_xxyzzz_0, \
                                         tg_xxxyzzz_xxzzzz_0, tg_xxxyzzz_xyyyyy_0, tg_xxxyzzz_xyyyyz_0, tg_xxxyzzz_xyyyzz_0, \
                                         tg_xxxyzzz_xyyzzz_0, tg_xxxyzzz_xyzzzz_0, tg_xxxyzzz_xzzzzz_0, tg_xxxyzzz_yyyyyy_0, \
                                         tg_xxxyzzz_yyyyyz_0, tg_xxxyzzz_yyyyzz_0, tg_xxxyzzz_yyyzzz_0, tg_xxxyzzz_yyzzzz_0, \
                                         tg_xxxyzzz_yzzzzz_0, tg_xxxyzzz_zzzzzz_0, tg_xxxzzzz_xxxxxx_0, tg_xxxzzzz_xxxxxy_0, \
                                         tg_xxxzzzz_xxxxxz_0, tg_xxxzzzz_xxxxyy_0, tg_xxxzzzz_xxxxyz_0, tg_xxxzzzz_xxxxzz_0, \
                                         tg_xxxzzzz_xxxyyy_0, tg_xxxzzzz_xxxyyz_0, tg_xxxzzzz_xxxyzz_0, tg_xxxzzzz_xxxzzz_0, \
                                         tg_xxxzzzz_xxyyyy_0, tg_xxxzzzz_xxyyyz_0, tg_xxxzzzz_xxyyzz_0, tg_xxxzzzz_xxyzzz_0, \
                                         tg_xxxzzzz_xxzzzz_0, tg_xxxzzzz_xyyyyy_0, tg_xxxzzzz_xyyyyz_0, tg_xxxzzzz_xyyyzz_0, \
                                         tg_xxxzzzz_xyyzzz_0, tg_xxxzzzz_xyzzzz_0, tg_xxxzzzz_xzzzzz_0, tg_xxxzzzz_yyyyyy_0, \
                                         tg_xxxzzzz_yyyyyz_0, tg_xxxzzzz_yyyyzz_0, tg_xxxzzzz_yyyzzz_0, tg_xxxzzzz_yyzzzz_0, \
                                         tg_xxxzzzz_yzzzzz_0, tg_xxxzzzz_zzzzzz_0, tg_xxyyyyy_xxxxxx_0, tg_xxyyyyy_xxxxxy_0, \
                                         tg_xxyyyyy_xxxxxz_0, tg_xxyyyyy_xxxxyy_0, tg_xxyyyyy_xxxxyz_0, tg_xxyyyyy_xxxxzz_0, \
                                         tg_xxyyyyy_xxxyyy_0, tg_xxyyyyy_xxxyyz_0, tg_xxyyyyy_xxxyzz_0, tg_xxyyyyy_xxxzzz_0, \
                                         tg_xxyyyyy_xxyyyy_0, tg_xxyyyyy_xxyyyz_0, tg_xxyyyyy_xxyyzz_0, tg_xxyyyyy_xxyzzz_0, \
                                         tg_xxyyyyy_xxzzzz_0, tg_xxyyyyy_xyyyyy_0, tg_xxyyyyy_xyyyyz_0, tg_xxyyyyy_xyyyzz_0, \
                                         tg_xxyyyyy_xyyzzz_0, tg_xxyyyyy_xyzzzz_0, tg_xxyyyyy_xzzzzz_0, tg_xxyyyyy_yyyyyy_0, \
                                         tg_xxyyyyy_yyyyyz_0, tg_xxyyyyy_yyyyzz_0, tg_xxyyyyy_yyyzzz_0, tg_xxyyyyy_yyzzzz_0, \
                                         tg_xxyyyyy_yzzzzz_0, tg_xxyyyyy_zzzzzz_0, tg_xxyyyyz_xxxxxx_0, tg_xxyyyyz_xxxxxy_0, \
                                         tg_xxyyyyz_xxxxxz_0, tg_xxyyyyz_xxxxyy_0, tg_xxyyyyz_xxxxyz_0, tg_xxyyyyz_xxxxzz_0, \
                                         tg_xxyyyyz_xxxyyy_0, tg_xxyyyyz_xxxyyz_0, tg_xxyyyyz_xxxyzz_0, tg_xxyyyyz_xxxzzz_0, \
                                         tg_xxyyyyz_xxyyyy_0, tg_xxyyyyz_xxyyyz_0, tg_xxyzzz_xxxxyz_0, tg_xxyzzz_xxxxyz_1, \
                                         tg_xxyzzz_xxxxzz_0, tg_xxyzzz_xxxxzz_1, tg_xxyzzz_xxxyyy_0, tg_xxyzzz_xxxyyy_1, \
                                         tg_xxyzzz_xxxyyz_0, tg_xxyzzz_xxxyyz_1, tg_xxyzzz_xxxyz_1, tg_xxyzzz_xxxyzz_0, \
                                         tg_xxyzzz_xxxyzz_1, tg_xxyzzz_xxxzz_1, tg_xxyzzz_xxxzzz_0, tg_xxyzzz_xxxzzz_1, \
                                         tg_xxyzzz_xxyyy_1, tg_xxyzzz_xxyyyy_0, tg_xxyzzz_xxyyyy_1, tg_xxyzzz_xxyyyz_0, \
                                         tg_xxyzzz_xxyyyz_1, tg_xxyzzz_xxyyz_1, tg_xxyzzz_xxyyzz_0, tg_xxyzzz_xxyyzz_1, \
                                         tg_xxyzzz_xxyzz_1, tg_xxyzzz_xxyzzz_0, tg_xxyzzz_xxyzzz_1, tg_xxyzzz_xxzzz_1, \
                                         tg_xxyzzz_xxzzzz_0, tg_xxyzzz_xxzzzz_1, tg_xxyzzz_xyyyy_1, tg_xxyzzz_xyyyyy_0, \
                                         tg_xxyzzz_xyyyyy_1, tg_xxyzzz_xyyyyz_0, tg_xxyzzz_xyyyyz_1, tg_xxyzzz_xyyyz_1, \
                                         tg_xxyzzz_xyyyzz_0, tg_xxyzzz_xyyyzz_1, tg_xxyzzz_xyyzz_1, tg_xxyzzz_xyyzzz_0, \
                                         tg_xxyzzz_xyyzzz_1, tg_xxyzzz_xyzzz_1, tg_xxyzzz_xyzzzz_0, tg_xxyzzz_xyzzzz_1, \
                                         tg_xxyzzz_xzzzz_1, tg_xxyzzz_xzzzzz_0, tg_xxyzzz_xzzzzz_1, tg_xxyzzz_yyyyy_1, \
                                         tg_xxyzzz_yyyyyy_0, tg_xxyzzz_yyyyyy_1, tg_xxyzzz_yyyyyz_0, tg_xxyzzz_yyyyyz_1, \
                                         tg_xxyzzz_yyyyz_1, tg_xxyzzz_yyyyzz_0, tg_xxyzzz_yyyyzz_1, tg_xxyzzz_yyyzz_1, \
                                         tg_xxyzzz_yyyzzz_0, tg_xxyzzz_yyyzzz_1, tg_xxyzzz_yyzzz_1, tg_xxyzzz_yyzzzz_0, \
                                         tg_xxyzzz_yyzzzz_1, tg_xxyzzz_yzzzz_1, tg_xxyzzz_yzzzzz_0, tg_xxyzzz_yzzzzz_1, \
                                         tg_xxyzzz_zzzzz_1, tg_xxyzzz_zzzzzz_0, tg_xxyzzz_zzzzzz_1, tg_xxzzzz_xxxxx_1, \
                                         tg_xxzzzz_xxxxxx_0, tg_xxzzzz_xxxxxx_1, tg_xxzzzz_xxxxxy_0, tg_xxzzzz_xxxxxy_1, \
                                         tg_xxzzzz_xxxxxz_0, tg_xxzzzz_xxxxxz_1, tg_xxzzzz_xxxxy_1, tg_xxzzzz_xxxxyy_0, \
                                         tg_xxzzzz_xxxxyy_1, tg_xxzzzz_xxxxyz_0, tg_xxzzzz_xxxxyz_1, tg_xxzzzz_xxxxz_1, \
                                         tg_xxzzzz_xxxxzz_0, tg_xxzzzz_xxxxzz_1, tg_xxzzzz_xxxyy_1, tg_xxzzzz_xxxyyy_0, \
                                         tg_xxzzzz_xxxyyy_1, tg_xxzzzz_xxxyyz_0, tg_xxzzzz_xxxyyz_1, tg_xxzzzz_xxxyz_1, \
                                         tg_xxzzzz_xxxyzz_0, tg_xxzzzz_xxxyzz_1, tg_xxzzzz_xxxzz_1, tg_xxzzzz_xxxzzz_0, \
                                         tg_xxzzzz_xxxzzz_1, tg_xxzzzz_xxyyy_1, tg_xxzzzz_xxyyyy_0, tg_xxzzzz_xxyyyy_1, \
                                         tg_xxzzzz_xxyyyz_0, tg_xxzzzz_xxyyyz_1, tg_xxzzzz_xxyyz_1, tg_xxzzzz_xxyyzz_0, \
                                         tg_xxzzzz_xxyyzz_1, tg_xxzzzz_xxyzz_1, tg_xxzzzz_xxyzzz_0, tg_xxzzzz_xxyzzz_1, \
                                         tg_xxzzzz_xxzzz_1, tg_xxzzzz_xxzzzz_0, tg_xxzzzz_xxzzzz_1, tg_xxzzzz_xyyyy_1, \
                                         tg_xxzzzz_xyyyyy_0, tg_xxzzzz_xyyyyy_1, tg_xxzzzz_xyyyyz_0, tg_xxzzzz_xyyyyz_1, \
                                         tg_xxzzzz_xyyyz_1, tg_xxzzzz_xyyyzz_0, tg_xxzzzz_xyyyzz_1, tg_xxzzzz_xyyzz_1, \
                                         tg_xxzzzz_xyyzzz_0, tg_xxzzzz_xyyzzz_1, tg_xxzzzz_xyzzz_1, tg_xxzzzz_xyzzzz_0, \
                                         tg_xxzzzz_xyzzzz_1, tg_xxzzzz_xzzzz_1, tg_xxzzzz_xzzzzz_0, tg_xxzzzz_xzzzzz_1, \
                                         tg_xxzzzz_yyyyy_1, tg_xxzzzz_yyyyyy_0, tg_xxzzzz_yyyyyy_1, tg_xxzzzz_yyyyyz_0, \
                                         tg_xxzzzz_yyyyyz_1, tg_xxzzzz_yyyyz_1, tg_xxzzzz_yyyyzz_0, tg_xxzzzz_yyyyzz_1, \
                                         tg_xxzzzz_yyyzz_1, tg_xxzzzz_yyyzzz_0, tg_xxzzzz_yyyzzz_1, tg_xxzzzz_yyzzz_1, \
                                         tg_xxzzzz_yyzzzz_0, tg_xxzzzz_yyzzzz_1, tg_xxzzzz_yzzzz_1, tg_xxzzzz_yzzzzz_0, \
                                         tg_xxzzzz_yzzzzz_1, tg_xxzzzz_zzzzz_1, tg_xxzzzz_zzzzzz_0, tg_xxzzzz_zzzzzz_1, \
                                         tg_xyyyyy_xxxxx_1, tg_xyyyyy_xxxxxx_0, tg_xyyyyy_xxxxxx_1, tg_xyyyyy_xxxxxy_0, \
                                         tg_xyyyyy_xxxxxy_1, tg_xyyyyy_xxxxxz_0, tg_xyyyyy_xxxxxz_1, tg_xyyyyy_xxxxy_1, \
                                         tg_xyyyyy_xxxxyy_0, tg_xyyyyy_xxxxyy_1, tg_xyyyyy_xxxxyz_0, tg_xyyyyy_xxxxyz_1, \
                                         tg_xyyyyy_xxxxz_1, tg_xyyyyy_xxxxzz_0, tg_xyyyyy_xxxxzz_1, tg_xyyyyy_xxxyy_1, \
                                         tg_xyyyyy_xxxyyy_0, tg_xyyyyy_xxxyyy_1, tg_xyyyyy_xxxyyz_0, tg_xyyyyy_xxxyyz_1, \
                                         tg_xyyyyy_xxxyz_1, tg_xyyyyy_xxxyzz_0, tg_xyyyyy_xxxyzz_1, tg_xyyyyy_xxxzz_1, \
                                         tg_xyyyyy_xxxzzz_0, tg_xyyyyy_xxxzzz_1, tg_xyyyyy_xxyyy_1, tg_xyyyyy_xxyyyy_0, \
                                         tg_xyyyyy_xxyyyy_1, tg_xyyyyy_xxyyyz_0, tg_xyyyyy_xxyyyz_1, tg_xyyyyy_xxyyz_1, \
                                         tg_xyyyyy_xxyyzz_0, tg_xyyyyy_xxyyzz_1, tg_xyyyyy_xxyzz_1, tg_xyyyyy_xxyzzz_0, \
                                         tg_xyyyyy_xxyzzz_1, tg_xyyyyy_xxzzz_1, tg_xyyyyy_xxzzzz_0, tg_xyyyyy_xxzzzz_1, \
                                         tg_xyyyyy_xyyyy_1, tg_xyyyyy_xyyyyy_0, tg_xyyyyy_xyyyyy_1, tg_xyyyyy_xyyyyz_0, \
                                         tg_xyyyyy_xyyyyz_1, tg_xyyyyy_xyyyz_1, tg_xyyyyy_xyyyzz_0, tg_xyyyyy_xyyyzz_1, \
                                         tg_xyyyyy_xyyzz_1, tg_xyyyyy_xyyzzz_0, tg_xyyyyy_xyyzzz_1, tg_xyyyyy_xyzzz_1, \
                                         tg_xyyyyy_xyzzzz_0, tg_xyyyyy_xyzzzz_1, tg_xyyyyy_xzzzz_1, tg_xyyyyy_xzzzzz_0, \
                                         tg_xyyyyy_xzzzzz_1, tg_xyyyyy_yyyyy_1, tg_xyyyyy_yyyyyy_0, tg_xyyyyy_yyyyyy_1, \
                                         tg_xyyyyy_yyyyyz_0, tg_xyyyyy_yyyyyz_1, tg_xyyyyy_yyyyz_1, tg_xyyyyy_yyyyzz_0, \
                                         tg_xyyyyy_yyyyzz_1, tg_xyyyyy_yyyzz_1, tg_xyyyyy_yyyzzz_0, tg_xyyyyy_yyyzzz_1, \
                                         tg_xyyyyy_yyzzz_1, tg_xyyyyy_yyzzzz_0, tg_xyyyyy_yyzzzz_1, tg_xyyyyy_yzzzz_1, \
                                         tg_xyyyyy_yzzzzz_0, tg_xyyyyy_yzzzzz_1, tg_xyyyyy_zzzzz_1, tg_xyyyyy_zzzzzz_0, \
                                         tg_xyyyyy_zzzzzz_1, tg_xyyyyz_xxxxx_1, tg_xyyyyz_xxxxxx_0, tg_xyyyyz_xxxxxx_1, \
                                         tg_xyyyyz_xxxxxy_0, tg_xyyyyz_xxxxxy_1, tg_xyyyyz_xxxxxz_0, tg_xyyyyz_xxxxxz_1, \
                                         tg_xyyyyz_xxxxy_1, tg_xyyyyz_xxxxyy_0, tg_xyyyyz_xxxxyy_1, tg_xyyyyz_xxxxyz_0, \
                                         tg_xyyyyz_xxxxyz_1, tg_xyyyyz_xxxxz_1, tg_xyyyyz_xxxxzz_0, tg_xyyyyz_xxxxzz_1, \
                                         tg_xyyyyz_xxxyy_1, tg_xyyyyz_xxxyyy_0, tg_xyyyyz_xxxyyy_1, tg_xyyyyz_xxxyyz_0, \
                                         tg_xyyyyz_xxxyyz_1, tg_xyyyyz_xxxyz_1, tg_xyyyyz_xxxyzz_0, tg_xyyyyz_xxxyzz_1, \
                                         tg_xyyyyz_xxxzz_1, tg_xyyyyz_xxxzzz_0, tg_xyyyyz_xxxzzz_1, tg_xyyyyz_xxyyy_1, \
                                         tg_xyyyyz_xxyyyy_0, tg_xyyyyz_xxyyyy_1, tg_xyyyyz_xxyyyz_0, tg_xyyyyz_xxyyyz_1, \
                                         tg_xyyyyz_xxyyz_1, tg_xyyyyz_xxyzz_1, tg_xyyyyz_xxzzz_1, tg_xyyyyz_xyyyy_1, \
                                         tg_xyyyyz_xyyyz_1, tg_xyzzz_xxxxyz_0, tg_xyzzz_xxxxyz_1, tg_xyzzz_xxxxzz_0, \
                                         tg_xyzzz_xxxxzz_1, tg_xyzzz_xxxyyy_0, tg_xyzzz_xxxyyy_1, tg_xyzzz_xxxyyz_0, \
                                         tg_xyzzz_xxxyyz_1, tg_xyzzz_xxxyzz_0, tg_xyzzz_xxxyzz_1, tg_xyzzz_xxxzzz_0, \
                                         tg_xyzzz_xxxzzz_1, tg_xyzzz_xxyyyy_0, tg_xyzzz_xxyyyy_1, tg_xyzzz_xxyyyz_0, \
                                         tg_xyzzz_xxyyyz_1, tg_xyzzz_xxyyzz_0, tg_xyzzz_xxyyzz_1, tg_xyzzz_xxyzzz_0, \
                                         tg_xyzzz_xxyzzz_1, tg_xyzzz_xxzzzz_0, tg_xyzzz_xxzzzz_1, tg_xyzzz_xyyyyy_0, \
                                         tg_xyzzz_xyyyyy_1, tg_xyzzz_xyyyyz_0, tg_xyzzz_xyyyyz_1, tg_xyzzz_xyyyzz_0, \
                                         tg_xyzzz_xyyyzz_1, tg_xyzzz_xyyzzz_0, tg_xyzzz_xyyzzz_1, tg_xyzzz_xyzzzz_0, \
                                         tg_xyzzz_xyzzzz_1, tg_xyzzz_xzzzzz_0, tg_xyzzz_xzzzzz_1, tg_xyzzz_yyyyyy_0, \
                                         tg_xyzzz_yyyyyy_1, tg_xyzzz_yyyyyz_0, tg_xyzzz_yyyyyz_1, tg_xyzzz_yyyyzz_0, \
                                         tg_xyzzz_yyyyzz_1, tg_xyzzz_yyyzzz_0, tg_xyzzz_yyyzzz_1, tg_xyzzz_yyzzzz_0, \
                                         tg_xyzzz_yyzzzz_1, tg_xyzzz_yzzzzz_0, tg_xyzzz_yzzzzz_1, tg_xyzzz_zzzzzz_0, \
                                         tg_xyzzz_zzzzzz_1, tg_xzzzz_xxxxxx_0, tg_xzzzz_xxxxxx_1, tg_xzzzz_xxxxxy_0, \
                                         tg_xzzzz_xxxxxy_1, tg_xzzzz_xxxxxz_0, tg_xzzzz_xxxxxz_1, tg_xzzzz_xxxxyy_0, \
                                         tg_xzzzz_xxxxyy_1, tg_xzzzz_xxxxyz_0, tg_xzzzz_xxxxyz_1, tg_xzzzz_xxxxzz_0, \
                                         tg_xzzzz_xxxxzz_1, tg_xzzzz_xxxyyy_0, tg_xzzzz_xxxyyy_1, tg_xzzzz_xxxyyz_0, \
                                         tg_xzzzz_xxxyyz_1, tg_xzzzz_xxxyzz_0, tg_xzzzz_xxxyzz_1, tg_xzzzz_xxxzzz_0, \
                                         tg_xzzzz_xxxzzz_1, tg_xzzzz_xxyyyy_0, tg_xzzzz_xxyyyy_1, tg_xzzzz_xxyyyz_0, \
                                         tg_xzzzz_xxyyyz_1, tg_xzzzz_xxyyzz_0, tg_xzzzz_xxyyzz_1, tg_xzzzz_xxyzzz_0, \
                                         tg_xzzzz_xxyzzz_1, tg_xzzzz_xxzzzz_0, tg_xzzzz_xxzzzz_1, tg_xzzzz_xyyyyy_0, \
                                         tg_xzzzz_xyyyyy_1, tg_xzzzz_xyyyyz_0, tg_xzzzz_xyyyyz_1, tg_xzzzz_xyyyzz_0, \
                                         tg_xzzzz_xyyyzz_1, tg_xzzzz_xyyzzz_0, tg_xzzzz_xyyzzz_1, tg_xzzzz_xyzzzz_0, \
                                         tg_xzzzz_xyzzzz_1, tg_xzzzz_xzzzzz_0, tg_xzzzz_xzzzzz_1, tg_xzzzz_yyyyyy_0, \
                                         tg_xzzzz_yyyyyy_1, tg_xzzzz_yyyyyz_0, tg_xzzzz_yyyyyz_1, tg_xzzzz_yyyyzz_0, \
                                         tg_xzzzz_yyyyzz_1, tg_xzzzz_yyyzzz_0, tg_xzzzz_yyyzzz_1, tg_xzzzz_yyzzzz_0, \
                                         tg_xzzzz_yyzzzz_1, tg_xzzzz_yzzzzz_0, tg_xzzzz_yzzzzz_1, tg_xzzzz_zzzzzz_0, \
                                         tg_xzzzz_zzzzzz_1, tg_yyyyy_xxxxxx_0, tg_yyyyy_xxxxxx_1, tg_yyyyy_xxxxxy_0, \
                                         tg_yyyyy_xxxxxy_1, tg_yyyyy_xxxxxz_0, tg_yyyyy_xxxxxz_1, tg_yyyyy_xxxxyy_0, \
                                         tg_yyyyy_xxxxyy_1, tg_yyyyy_xxxxyz_0, tg_yyyyy_xxxxyz_1, tg_yyyyy_xxxxzz_0, \
                                         tg_yyyyy_xxxxzz_1, tg_yyyyy_xxxyyy_0, tg_yyyyy_xxxyyy_1, tg_yyyyy_xxxyyz_0, \
                                         tg_yyyyy_xxxyyz_1, tg_yyyyy_xxxyzz_0, tg_yyyyy_xxxyzz_1, tg_yyyyy_xxxzzz_0, \
                                         tg_yyyyy_xxxzzz_1, tg_yyyyy_xxyyyy_0, tg_yyyyy_xxyyyy_1, tg_yyyyy_xxyyyz_0, \
                                         tg_yyyyy_xxyyyz_1, tg_yyyyy_xxyyzz_0, tg_yyyyy_xxyyzz_1, tg_yyyyy_xxyzzz_0, \
                                         tg_yyyyy_xxyzzz_1, tg_yyyyy_xxzzzz_0, tg_yyyyy_xxzzzz_1, tg_yyyyy_xyyyyy_0, \
                                         tg_yyyyy_xyyyyy_1, tg_yyyyy_xyyyyz_0, tg_yyyyy_xyyyyz_1, tg_yyyyy_xyyyzz_0, \
                                         tg_yyyyy_xyyyzz_1, tg_yyyyy_xyyzzz_0, tg_yyyyy_xyyzzz_1, tg_yyyyy_xyzzzz_0, \
                                         tg_yyyyy_xyzzzz_1, tg_yyyyy_xzzzzz_0, tg_yyyyy_xzzzzz_1, tg_yyyyy_yyyyyy_0, \
                                         tg_yyyyy_yyyyyy_1, tg_yyyyy_yyyyyz_0, tg_yyyyy_yyyyyz_1, tg_yyyyy_yyyyzz_0, \
                                         tg_yyyyy_yyyyzz_1, tg_yyyyy_yyyzzz_0, tg_yyyyy_yyyzzz_1, tg_yyyyy_yyzzzz_0, \
                                         tg_yyyyy_yyzzzz_1, tg_yyyyy_yzzzzz_0, tg_yyyyy_yzzzzz_1, tg_yyyyy_zzzzzz_0, \
                                         tg_yyyyy_zzzzzz_1, tg_yyyyz_xxxxxx_0, tg_yyyyz_xxxxxx_1, tg_yyyyz_xxxxxy_0, \
                                         tg_yyyyz_xxxxxy_1, tg_yyyyz_xxxxxz_0, tg_yyyyz_xxxxxz_1, tg_yyyyz_xxxxyy_0, \
                                         tg_yyyyz_xxxxyy_1, tg_yyyyz_xxxxyz_0, tg_yyyyz_xxxxyz_1, tg_yyyyz_xxxxzz_0, \
                                         tg_yyyyz_xxxxzz_1, tg_yyyyz_xxxyyy_0, tg_yyyyz_xxxyyy_1, tg_yyyyz_xxxyyz_0, \
                                         tg_yyyyz_xxxyyz_1, tg_yyyyz_xxxyzz_0, tg_yyyyz_xxxyzz_1, tg_yyyyz_xxxzzz_0, \
                                         tg_yyyyz_xxxzzz_1, tg_yyyyz_xxyyyy_0, tg_yyyyz_xxyyyy_1, tg_yyyyz_xxyyyz_0, \
                                         tg_yyyyz_xxyyyz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyzzz_xxxxyz_0[j] = pb_x * tg_xxyzzz_xxxxyz_0[j] + fr * tg_xxyzzz_xxxxyz_1[j] + fl1_fx * (tg_xyzzz_xxxxyz_0[j] - tg_xyzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzz_xxxyz_1[j];

                    tg_xxxyzzz_xxxxzz_0[j] = pb_x * tg_xxyzzz_xxxxzz_0[j] + fr * tg_xxyzzz_xxxxzz_1[j] + fl1_fx * (tg_xyzzz_xxxxzz_0[j] - tg_xyzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzz_xxxzz_1[j];

                    tg_xxxyzzz_xxxyyy_0[j] = pb_x * tg_xxyzzz_xxxyyy_0[j] + fr * tg_xxyzzz_xxxyyy_1[j] + fl1_fx * (tg_xyzzz_xxxyyy_0[j] - tg_xyzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxyyy_1[j];

                    tg_xxxyzzz_xxxyyz_0[j] = pb_x * tg_xxyzzz_xxxyyz_0[j] + fr * tg_xxyzzz_xxxyyz_1[j] + fl1_fx * (tg_xyzzz_xxxyyz_0[j] - tg_xyzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxyyz_1[j];

                    tg_xxxyzzz_xxxyzz_0[j] = pb_x * tg_xxyzzz_xxxyzz_0[j] + fr * tg_xxyzzz_xxxyzz_1[j] + fl1_fx * (tg_xyzzz_xxxyzz_0[j] - tg_xyzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxyzz_1[j];

                    tg_xxxyzzz_xxxzzz_0[j] = pb_x * tg_xxyzzz_xxxzzz_0[j] + fr * tg_xxyzzz_xxxzzz_1[j] + fl1_fx * (tg_xyzzz_xxxzzz_0[j] - tg_xyzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxzzz_1[j];

                    tg_xxxyzzz_xxyyyy_0[j] = pb_x * tg_xxyzzz_xxyyyy_0[j] + fr * tg_xxyzzz_xxyyyy_1[j] + fl1_fx * (tg_xyzzz_xxyyyy_0[j] - tg_xyzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyyyy_1[j];

                    tg_xxxyzzz_xxyyyz_0[j] = pb_x * tg_xxyzzz_xxyyyz_0[j] + fr * tg_xxyzzz_xxyyyz_1[j] + fl1_fx * (tg_xyzzz_xxyyyz_0[j] - tg_xyzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyyyz_1[j];

                    tg_xxxyzzz_xxyyzz_0[j] = pb_x * tg_xxyzzz_xxyyzz_0[j] + fr * tg_xxyzzz_xxyyzz_1[j] + fl1_fx * (tg_xyzzz_xxyyzz_0[j] - tg_xyzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyyzz_1[j];

                    tg_xxxyzzz_xxyzzz_0[j] = pb_x * tg_xxyzzz_xxyzzz_0[j] + fr * tg_xxyzzz_xxyzzz_1[j] + fl1_fx * (tg_xyzzz_xxyzzz_0[j] - tg_xyzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyzzz_1[j];

                    tg_xxxyzzz_xxzzzz_0[j] = pb_x * tg_xxyzzz_xxzzzz_0[j] + fr * tg_xxyzzz_xxzzzz_1[j] + fl1_fx * (tg_xyzzz_xxzzzz_0[j] - tg_xyzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xzzzz_1[j];

                    tg_xxxyzzz_xyyyyy_0[j] = pb_x * tg_xxyzzz_xyyyyy_0[j] + fr * tg_xxyzzz_xyyyyy_1[j] + fl1_fx * (tg_xyzzz_xyyyyy_0[j] - tg_xyzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyyyy_1[j];

                    tg_xxxyzzz_xyyyyz_0[j] = pb_x * tg_xxyzzz_xyyyyz_0[j] + fr * tg_xxyzzz_xyyyyz_1[j] + fl1_fx * (tg_xyzzz_xyyyyz_0[j] - tg_xyzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyyyz_1[j];

                    tg_xxxyzzz_xyyyzz_0[j] = pb_x * tg_xxyzzz_xyyyzz_0[j] + fr * tg_xxyzzz_xyyyzz_1[j] + fl1_fx * (tg_xyzzz_xyyyzz_0[j] - tg_xyzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyyzz_1[j];

                    tg_xxxyzzz_xyyzzz_0[j] = pb_x * tg_xxyzzz_xyyzzz_0[j] + fr * tg_xxyzzz_xyyzzz_1[j] + fl1_fx * (tg_xyzzz_xyyzzz_0[j] - tg_xyzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyzzz_1[j];

                    tg_xxxyzzz_xyzzzz_0[j] = pb_x * tg_xxyzzz_xyzzzz_0[j] + fr * tg_xxyzzz_xyzzzz_1[j] + fl1_fx * (tg_xyzzz_xyzzzz_0[j] - tg_xyzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yzzzz_1[j];

                    tg_xxxyzzz_xzzzzz_0[j] = pb_x * tg_xxyzzz_xzzzzz_0[j] + fr * tg_xxyzzz_xzzzzz_1[j] + fl1_fx * (tg_xyzzz_xzzzzz_0[j] - tg_xyzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_zzzzz_1[j];

                    tg_xxxyzzz_yyyyyy_0[j] = pb_x * tg_xxyzzz_yyyyyy_0[j] + fr * tg_xxyzzz_yyyyyy_1[j] + fl1_fx * (tg_xyzzz_yyyyyy_0[j] - tg_xyzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyzzz_yyyyyz_0[j] = pb_x * tg_xxyzzz_yyyyyz_0[j] + fr * tg_xxyzzz_yyyyyz_1[j] + fl1_fx * (tg_xyzzz_yyyyyz_0[j] - tg_xyzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyzzz_yyyyzz_0[j] = pb_x * tg_xxyzzz_yyyyzz_0[j] + fr * tg_xxyzzz_yyyyzz_1[j] + fl1_fx * (tg_xyzzz_yyyyzz_0[j] - tg_xyzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyzzz_yyyzzz_0[j] = pb_x * tg_xxyzzz_yyyzzz_0[j] + fr * tg_xxyzzz_yyyzzz_1[j] + fl1_fx * (tg_xyzzz_yyyzzz_0[j] - tg_xyzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_yyzzzz_0[j] = pb_x * tg_xxyzzz_yyzzzz_0[j] + fr * tg_xxyzzz_yyzzzz_1[j] + fl1_fx * (tg_xyzzz_yyzzzz_0[j] - tg_xyzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_yzzzzz_0[j] = pb_x * tg_xxyzzz_yzzzzz_0[j] + fr * tg_xxyzzz_yzzzzz_1[j] + fl1_fx * (tg_xyzzz_yzzzzz_0[j] - tg_xyzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_zzzzzz_0[j] = pb_x * tg_xxyzzz_zzzzzz_0[j] + fr * tg_xxyzzz_zzzzzz_1[j] + fl1_fx * (tg_xyzzz_zzzzzz_0[j] - tg_xyzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_xxxxxx_0[j] = pb_x * tg_xxzzzz_xxxxxx_0[j] + fr * tg_xxzzzz_xxxxxx_1[j] + fl1_fx * (tg_xzzzz_xxxxxx_0[j] - tg_xzzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxzzzz_xxxxx_1[j];

                    tg_xxxzzzz_xxxxxy_0[j] = pb_x * tg_xxzzzz_xxxxxy_0[j] + fr * tg_xxzzzz_xxxxxy_1[j] + fl1_fx * (tg_xzzzz_xxxxxy_0[j] - tg_xzzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzzzz_xxxxy_1[j];

                    tg_xxxzzzz_xxxxxz_0[j] = pb_x * tg_xxzzzz_xxxxxz_0[j] + fr * tg_xxzzzz_xxxxxz_1[j] + fl1_fx * (tg_xzzzz_xxxxxz_0[j] - tg_xzzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzzzz_xxxxz_1[j];

                    tg_xxxzzzz_xxxxyy_0[j] = pb_x * tg_xxzzzz_xxxxyy_0[j] + fr * tg_xxzzzz_xxxxyy_1[j] + fl1_fx * (tg_xzzzz_xxxxyy_0[j] - tg_xzzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzz_xxxyy_1[j];

                    tg_xxxzzzz_xxxxyz_0[j] = pb_x * tg_xxzzzz_xxxxyz_0[j] + fr * tg_xxzzzz_xxxxyz_1[j] + fl1_fx * (tg_xzzzz_xxxxyz_0[j] - tg_xzzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzz_xxxyz_1[j];

                    tg_xxxzzzz_xxxxzz_0[j] = pb_x * tg_xxzzzz_xxxxzz_0[j] + fr * tg_xxzzzz_xxxxzz_1[j] + fl1_fx * (tg_xzzzz_xxxxzz_0[j] - tg_xzzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzz_xxxzz_1[j];

                    tg_xxxzzzz_xxxyyy_0[j] = pb_x * tg_xxzzzz_xxxyyy_0[j] + fr * tg_xxzzzz_xxxyyy_1[j] + fl1_fx * (tg_xzzzz_xxxyyy_0[j] - tg_xzzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxyyy_1[j];

                    tg_xxxzzzz_xxxyyz_0[j] = pb_x * tg_xxzzzz_xxxyyz_0[j] + fr * tg_xxzzzz_xxxyyz_1[j] + fl1_fx * (tg_xzzzz_xxxyyz_0[j] - tg_xzzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxyyz_1[j];

                    tg_xxxzzzz_xxxyzz_0[j] = pb_x * tg_xxzzzz_xxxyzz_0[j] + fr * tg_xxzzzz_xxxyzz_1[j] + fl1_fx * (tg_xzzzz_xxxyzz_0[j] - tg_xzzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxyzz_1[j];

                    tg_xxxzzzz_xxxzzz_0[j] = pb_x * tg_xxzzzz_xxxzzz_0[j] + fr * tg_xxzzzz_xxxzzz_1[j] + fl1_fx * (tg_xzzzz_xxxzzz_0[j] - tg_xzzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxzzz_1[j];

                    tg_xxxzzzz_xxyyyy_0[j] = pb_x * tg_xxzzzz_xxyyyy_0[j] + fr * tg_xxzzzz_xxyyyy_1[j] + fl1_fx * (tg_xzzzz_xxyyyy_0[j] - tg_xzzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyyyy_1[j];

                    tg_xxxzzzz_xxyyyz_0[j] = pb_x * tg_xxzzzz_xxyyyz_0[j] + fr * tg_xxzzzz_xxyyyz_1[j] + fl1_fx * (tg_xzzzz_xxyyyz_0[j] - tg_xzzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyyyz_1[j];

                    tg_xxxzzzz_xxyyzz_0[j] = pb_x * tg_xxzzzz_xxyyzz_0[j] + fr * tg_xxzzzz_xxyyzz_1[j] + fl1_fx * (tg_xzzzz_xxyyzz_0[j] - tg_xzzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyyzz_1[j];

                    tg_xxxzzzz_xxyzzz_0[j] = pb_x * tg_xxzzzz_xxyzzz_0[j] + fr * tg_xxzzzz_xxyzzz_1[j] + fl1_fx * (tg_xzzzz_xxyzzz_0[j] - tg_xzzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyzzz_1[j];

                    tg_xxxzzzz_xxzzzz_0[j] = pb_x * tg_xxzzzz_xxzzzz_0[j] + fr * tg_xxzzzz_xxzzzz_1[j] + fl1_fx * (tg_xzzzz_xxzzzz_0[j] - tg_xzzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xzzzz_1[j];

                    tg_xxxzzzz_xyyyyy_0[j] = pb_x * tg_xxzzzz_xyyyyy_0[j] + fr * tg_xxzzzz_xyyyyy_1[j] + fl1_fx * (tg_xzzzz_xyyyyy_0[j] - tg_xzzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyyyy_1[j];

                    tg_xxxzzzz_xyyyyz_0[j] = pb_x * tg_xxzzzz_xyyyyz_0[j] + fr * tg_xxzzzz_xyyyyz_1[j] + fl1_fx * (tg_xzzzz_xyyyyz_0[j] - tg_xzzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyyyz_1[j];

                    tg_xxxzzzz_xyyyzz_0[j] = pb_x * tg_xxzzzz_xyyyzz_0[j] + fr * tg_xxzzzz_xyyyzz_1[j] + fl1_fx * (tg_xzzzz_xyyyzz_0[j] - tg_xzzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyyzz_1[j];

                    tg_xxxzzzz_xyyzzz_0[j] = pb_x * tg_xxzzzz_xyyzzz_0[j] + fr * tg_xxzzzz_xyyzzz_1[j] + fl1_fx * (tg_xzzzz_xyyzzz_0[j] - tg_xzzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyzzz_1[j];

                    tg_xxxzzzz_xyzzzz_0[j] = pb_x * tg_xxzzzz_xyzzzz_0[j] + fr * tg_xxzzzz_xyzzzz_1[j] + fl1_fx * (tg_xzzzz_xyzzzz_0[j] - tg_xzzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yzzzz_1[j];

                    tg_xxxzzzz_xzzzzz_0[j] = pb_x * tg_xxzzzz_xzzzzz_0[j] + fr * tg_xxzzzz_xzzzzz_1[j] + fl1_fx * (tg_xzzzz_xzzzzz_0[j] - tg_xzzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_zzzzz_1[j];

                    tg_xxxzzzz_yyyyyy_0[j] = pb_x * tg_xxzzzz_yyyyyy_0[j] + fr * tg_xxzzzz_yyyyyy_1[j] + fl1_fx * (tg_xzzzz_yyyyyy_0[j] - tg_xzzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxzzzz_yyyyyz_0[j] = pb_x * tg_xxzzzz_yyyyyz_0[j] + fr * tg_xxzzzz_yyyyyz_1[j] + fl1_fx * (tg_xzzzz_yyyyyz_0[j] - tg_xzzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxzzzz_yyyyzz_0[j] = pb_x * tg_xxzzzz_yyyyzz_0[j] + fr * tg_xxzzzz_yyyyzz_1[j] + fl1_fx * (tg_xzzzz_yyyyzz_0[j] - tg_xzzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxzzzz_yyyzzz_0[j] = pb_x * tg_xxzzzz_yyyzzz_0[j] + fr * tg_xxzzzz_yyyzzz_1[j] + fl1_fx * (tg_xzzzz_yyyzzz_0[j] - tg_xzzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_yyzzzz_0[j] = pb_x * tg_xxzzzz_yyzzzz_0[j] + fr * tg_xxzzzz_yyzzzz_1[j] + fl1_fx * (tg_xzzzz_yyzzzz_0[j] - tg_xzzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_yzzzzz_0[j] = pb_x * tg_xxzzzz_yzzzzz_0[j] + fr * tg_xxzzzz_yzzzzz_1[j] + fl1_fx * (tg_xzzzz_yzzzzz_0[j] - tg_xzzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_zzzzzz_0[j] = pb_x * tg_xxzzzz_zzzzzz_0[j] + fr * tg_xxzzzz_zzzzzz_1[j] + fl1_fx * (tg_xzzzz_zzzzzz_0[j] - tg_xzzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_xxxxxx_0[j] = pb_x * tg_xyyyyy_xxxxxx_0[j] + fr * tg_xyyyyy_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxxx_0[j] - tg_yyyyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyyyy_xxxxx_1[j];

                    tg_xxyyyyy_xxxxxy_0[j] = pb_x * tg_xyyyyy_xxxxxy_0[j] + fr * tg_xyyyyy_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxxy_0[j] - tg_yyyyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyy_xxxxy_1[j];

                    tg_xxyyyyy_xxxxxz_0[j] = pb_x * tg_xyyyyy_xxxxxz_0[j] + fr * tg_xyyyyy_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxxz_0[j] - tg_yyyyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyy_xxxxz_1[j];

                    tg_xxyyyyy_xxxxyy_0[j] = pb_x * tg_xyyyyy_xxxxyy_0[j] + fr * tg_xyyyyy_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxyy_0[j] - tg_yyyyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyy_xxxyy_1[j];

                    tg_xxyyyyy_xxxxyz_0[j] = pb_x * tg_xyyyyy_xxxxyz_0[j] + fr * tg_xyyyyy_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxyz_0[j] - tg_yyyyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyy_xxxyz_1[j];

                    tg_xxyyyyy_xxxxzz_0[j] = pb_x * tg_xyyyyy_xxxxzz_0[j] + fr * tg_xyyyyy_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxxzz_0[j] - tg_yyyyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyy_xxxzz_1[j];

                    tg_xxyyyyy_xxxyyy_0[j] = pb_x * tg_xyyyyy_xxxyyy_0[j] + fr * tg_xyyyyy_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxyyy_0[j] - tg_yyyyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxyyy_1[j];

                    tg_xxyyyyy_xxxyyz_0[j] = pb_x * tg_xyyyyy_xxxyyz_0[j] + fr * tg_xyyyyy_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxyyz_0[j] - tg_yyyyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxyyz_1[j];

                    tg_xxyyyyy_xxxyzz_0[j] = pb_x * tg_xyyyyy_xxxyzz_0[j] + fr * tg_xyyyyy_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxyzz_0[j] - tg_yyyyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxyzz_1[j];

                    tg_xxyyyyy_xxxzzz_0[j] = pb_x * tg_xyyyyy_xxxzzz_0[j] + fr * tg_xyyyyy_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxzzz_0[j] - tg_yyyyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxzzz_1[j];

                    tg_xxyyyyy_xxyyyy_0[j] = pb_x * tg_xyyyyy_xxyyyy_0[j] + fr * tg_xyyyyy_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyyyy_0[j] - tg_yyyyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyyyy_1[j];

                    tg_xxyyyyy_xxyyyz_0[j] = pb_x * tg_xyyyyy_xxyyyz_0[j] + fr * tg_xyyyyy_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyyyz_0[j] - tg_yyyyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyyyz_1[j];

                    tg_xxyyyyy_xxyyzz_0[j] = pb_x * tg_xyyyyy_xxyyzz_0[j] + fr * tg_xyyyyy_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyyzz_0[j] - tg_yyyyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyyzz_1[j];

                    tg_xxyyyyy_xxyzzz_0[j] = pb_x * tg_xyyyyy_xxyzzz_0[j] + fr * tg_xyyyyy_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyzzz_0[j] - tg_yyyyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyzzz_1[j];

                    tg_xxyyyyy_xxzzzz_0[j] = pb_x * tg_xyyyyy_xxzzzz_0[j] + fr * tg_xyyyyy_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxzzzz_0[j] - tg_yyyyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xzzzz_1[j];

                    tg_xxyyyyy_xyyyyy_0[j] = pb_x * tg_xyyyyy_xyyyyy_0[j] + fr * tg_xyyyyy_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyyyy_0[j] - tg_yyyyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyyyy_1[j];

                    tg_xxyyyyy_xyyyyz_0[j] = pb_x * tg_xyyyyy_xyyyyz_0[j] + fr * tg_xyyyyy_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyyyz_0[j] - tg_yyyyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyyyz_1[j];

                    tg_xxyyyyy_xyyyzz_0[j] = pb_x * tg_xyyyyy_xyyyzz_0[j] + fr * tg_xyyyyy_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyyzz_0[j] - tg_yyyyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyyzz_1[j];

                    tg_xxyyyyy_xyyzzz_0[j] = pb_x * tg_xyyyyy_xyyzzz_0[j] + fr * tg_xyyyyy_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyzzz_0[j] - tg_yyyyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyzzz_1[j];

                    tg_xxyyyyy_xyzzzz_0[j] = pb_x * tg_xyyyyy_xyzzzz_0[j] + fr * tg_xyyyyy_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyzzzz_0[j] - tg_yyyyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yzzzz_1[j];

                    tg_xxyyyyy_xzzzzz_0[j] = pb_x * tg_xyyyyy_xzzzzz_0[j] + fr * tg_xyyyyy_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xzzzzz_0[j] - tg_yyyyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_zzzzz_1[j];

                    tg_xxyyyyy_yyyyyy_0[j] = pb_x * tg_xyyyyy_yyyyyy_0[j] + fr * tg_xyyyyy_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyyyy_0[j] - tg_yyyyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyyyy_yyyyyz_0[j] = pb_x * tg_xyyyyy_yyyyyz_0[j] + fr * tg_xyyyyy_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyyyz_0[j] - tg_yyyyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyyyy_yyyyzz_0[j] = pb_x * tg_xyyyyy_yyyyzz_0[j] + fr * tg_xyyyyy_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyyzz_0[j] - tg_yyyyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyyyy_yyyzzz_0[j] = pb_x * tg_xyyyyy_yyyzzz_0[j] + fr * tg_xyyyyy_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyzzz_0[j] - tg_yyyyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_yyzzzz_0[j] = pb_x * tg_xyyyyy_yyzzzz_0[j] + fr * tg_xyyyyy_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyzzzz_0[j] - tg_yyyyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_yzzzzz_0[j] = pb_x * tg_xyyyyy_yzzzzz_0[j] + fr * tg_xyyyyy_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yzzzzz_0[j] - tg_yyyyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_zzzzzz_0[j] = pb_x * tg_xyyyyy_zzzzzz_0[j] + fr * tg_xyyyyy_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_zzzzzz_0[j] - tg_yyyyy_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_xxxxxx_0[j] = pb_x * tg_xyyyyz_xxxxxx_0[j] + fr * tg_xyyyyz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxxx_0[j] - tg_yyyyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyyyz_xxxxx_1[j];

                    tg_xxyyyyz_xxxxxy_0[j] = pb_x * tg_xyyyyz_xxxxxy_0[j] + fr * tg_xyyyyz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxxy_0[j] - tg_yyyyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyz_xxxxy_1[j];

                    tg_xxyyyyz_xxxxxz_0[j] = pb_x * tg_xyyyyz_xxxxxz_0[j] + fr * tg_xyyyyz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxxz_0[j] - tg_yyyyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyyz_xxxxz_1[j];

                    tg_xxyyyyz_xxxxyy_0[j] = pb_x * tg_xyyyyz_xxxxyy_0[j] + fr * tg_xyyyyz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxyy_0[j] - tg_yyyyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyz_xxxyy_1[j];

                    tg_xxyyyyz_xxxxyz_0[j] = pb_x * tg_xyyyyz_xxxxyz_0[j] + fr * tg_xyyyyz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxyz_0[j] - tg_yyyyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyz_xxxyz_1[j];

                    tg_xxyyyyz_xxxxzz_0[j] = pb_x * tg_xyyyyz_xxxxzz_0[j] + fr * tg_xyyyyz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxxzz_0[j] - tg_yyyyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyz_xxxzz_1[j];

                    tg_xxyyyyz_xxxyyy_0[j] = pb_x * tg_xyyyyz_xxxyyy_0[j] + fr * tg_xyyyyz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxyyy_0[j] - tg_yyyyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxyyy_1[j];

                    tg_xxyyyyz_xxxyyz_0[j] = pb_x * tg_xyyyyz_xxxyyz_0[j] + fr * tg_xyyyyz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxyyz_0[j] - tg_yyyyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxyyz_1[j];

                    tg_xxyyyyz_xxxyzz_0[j] = pb_x * tg_xyyyyz_xxxyzz_0[j] + fr * tg_xyyyyz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxyzz_0[j] - tg_yyyyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxyzz_1[j];

                    tg_xxyyyyz_xxxzzz_0[j] = pb_x * tg_xyyyyz_xxxzzz_0[j] + fr * tg_xyyyyz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxzzz_0[j] - tg_yyyyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxzzz_1[j];

                    tg_xxyyyyz_xxyyyy_0[j] = pb_x * tg_xyyyyz_xxyyyy_0[j] + fr * tg_xyyyyz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyyyy_0[j] - tg_yyyyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyyyy_1[j];

                    tg_xxyyyyz_xxyyyz_0[j] = pb_x * tg_xyyyyz_xxyyyz_0[j] + fr * tg_xyyyyz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyyyz_0[j] - tg_yyyyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_460_552(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (460,552)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xyyyyz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 460); 

                auto tg_xyyyyz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 461); 

                auto tg_xyyyyz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 462); 

                auto tg_xyyyyz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 463); 

                auto tg_xyyyyz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 464); 

                auto tg_xyyyyz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 465); 

                auto tg_xyyyyz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 466); 

                auto tg_xyyyyz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 467); 

                auto tg_xyyyyz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 468); 

                auto tg_xyyyyz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 469); 

                auto tg_xyyyyz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 470); 

                auto tg_xyyyyz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 471); 

                auto tg_xyyyyz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 472); 

                auto tg_xyyyyz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 473); 

                auto tg_xyyyyz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 474); 

                auto tg_xyyyyz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 475); 

                auto tg_xyyyzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 476); 

                auto tg_xyyyzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 477); 

                auto tg_xyyyzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 478); 

                auto tg_xyyyzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 479); 

                auto tg_xyyyzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 480); 

                auto tg_xyyyzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 481); 

                auto tg_xyyyzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 482); 

                auto tg_xyyyzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 483); 

                auto tg_xyyyzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 484); 

                auto tg_xyyyzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 485); 

                auto tg_xyyyzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 486); 

                auto tg_xyyyzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 487); 

                auto tg_xyyyzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 488); 

                auto tg_xyyyzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 489); 

                auto tg_xyyyzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 490); 

                auto tg_xyyyzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 491); 

                auto tg_xyyyzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 492); 

                auto tg_xyyyzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 493); 

                auto tg_xyyyzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 494); 

                auto tg_xyyyzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 495); 

                auto tg_xyyyzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 496); 

                auto tg_xyyyzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 497); 

                auto tg_xyyyzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 498); 

                auto tg_xyyyzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 499); 

                auto tg_xyyyzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 500); 

                auto tg_xyyyzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 501); 

                auto tg_xyyyzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 502); 

                auto tg_xyyyzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 503); 

                auto tg_xyyzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 504); 

                auto tg_xyyzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 505); 

                auto tg_xyyzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 506); 

                auto tg_xyyzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 507); 

                auto tg_xyyzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 508); 

                auto tg_xyyzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 509); 

                auto tg_xyyzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 510); 

                auto tg_xyyzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 511); 

                auto tg_xyyzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 512); 

                auto tg_xyyzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 513); 

                auto tg_xyyzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 514); 

                auto tg_xyyzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 515); 

                auto tg_xyyzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 516); 

                auto tg_xyyzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 517); 

                auto tg_xyyzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 518); 

                auto tg_xyyzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 519); 

                auto tg_xyyzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 520); 

                auto tg_xyyzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 521); 

                auto tg_xyyzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 522); 

                auto tg_xyyzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 523); 

                auto tg_xyyzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 524); 

                auto tg_xyyzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 525); 

                auto tg_xyyzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 526); 

                auto tg_xyyzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 527); 

                auto tg_xyyzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 528); 

                auto tg_xyyzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 529); 

                auto tg_xyyzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 530); 

                auto tg_xyyzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 531); 

                auto tg_xyzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 532); 

                auto tg_xyzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 533); 

                auto tg_xyzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 534); 

                auto tg_xyzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 535); 

                auto tg_xyzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 536); 

                auto tg_xyzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 537); 

                auto tg_xyzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 538); 

                auto tg_xyzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 539); 

                auto tg_xyzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 540); 

                auto tg_xyzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 541); 

                auto tg_xyzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 542); 

                auto tg_xyzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 543); 

                auto tg_xyzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 544); 

                auto tg_xyzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 545); 

                auto tg_xyzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 546); 

                auto tg_xyzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 547); 

                auto tg_xyzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 548); 

                auto tg_xyzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 549); 

                auto tg_xyzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 550); 

                auto tg_xyzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 551); 

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

                // set up pointers to integrals

                auto tg_xxyyyyz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 460); 

                auto tg_xxyyyyz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 461); 

                auto tg_xxyyyyz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 462); 

                auto tg_xxyyyyz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 463); 

                auto tg_xxyyyyz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 464); 

                auto tg_xxyyyyz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 465); 

                auto tg_xxyyyyz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 466); 

                auto tg_xxyyyyz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 467); 

                auto tg_xxyyyyz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 468); 

                auto tg_xxyyyyz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 469); 

                auto tg_xxyyyyz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 470); 

                auto tg_xxyyyyz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 471); 

                auto tg_xxyyyyz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 472); 

                auto tg_xxyyyyz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 473); 

                auto tg_xxyyyyz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 474); 

                auto tg_xxyyyyz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 475); 

                auto tg_xxyyyzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 476); 

                auto tg_xxyyyzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 477); 

                auto tg_xxyyyzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 478); 

                auto tg_xxyyyzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 479); 

                auto tg_xxyyyzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 480); 

                auto tg_xxyyyzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 481); 

                auto tg_xxyyyzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 482); 

                auto tg_xxyyyzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 483); 

                auto tg_xxyyyzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 484); 

                auto tg_xxyyyzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 485); 

                auto tg_xxyyyzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 486); 

                auto tg_xxyyyzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 487); 

                auto tg_xxyyyzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 488); 

                auto tg_xxyyyzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 489); 

                auto tg_xxyyyzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 490); 

                auto tg_xxyyyzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 491); 

                auto tg_xxyyyzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 492); 

                auto tg_xxyyyzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 493); 

                auto tg_xxyyyzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 494); 

                auto tg_xxyyyzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 495); 

                auto tg_xxyyyzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 496); 

                auto tg_xxyyyzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 497); 

                auto tg_xxyyyzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 498); 

                auto tg_xxyyyzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 499); 

                auto tg_xxyyyzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 500); 

                auto tg_xxyyyzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 501); 

                auto tg_xxyyyzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 502); 

                auto tg_xxyyyzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 503); 

                auto tg_xxyyzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 504); 

                auto tg_xxyyzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 505); 

                auto tg_xxyyzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 506); 

                auto tg_xxyyzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 507); 

                auto tg_xxyyzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 508); 

                auto tg_xxyyzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 509); 

                auto tg_xxyyzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 510); 

                auto tg_xxyyzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 511); 

                auto tg_xxyyzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 512); 

                auto tg_xxyyzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 513); 

                auto tg_xxyyzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 514); 

                auto tg_xxyyzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 515); 

                auto tg_xxyyzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 516); 

                auto tg_xxyyzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 517); 

                auto tg_xxyyzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 518); 

                auto tg_xxyyzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 519); 

                auto tg_xxyyzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 520); 

                auto tg_xxyyzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 521); 

                auto tg_xxyyzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 522); 

                auto tg_xxyyzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 523); 

                auto tg_xxyyzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 524); 

                auto tg_xxyyzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 525); 

                auto tg_xxyyzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 526); 

                auto tg_xxyyzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 527); 

                auto tg_xxyyzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 528); 

                auto tg_xxyyzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 529); 

                auto tg_xxyyzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 530); 

                auto tg_xxyyzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 531); 

                auto tg_xxyzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 532); 

                auto tg_xxyzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 533); 

                auto tg_xxyzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 534); 

                auto tg_xxyzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 535); 

                auto tg_xxyzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 536); 

                auto tg_xxyzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 537); 

                auto tg_xxyzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 538); 

                auto tg_xxyzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 539); 

                auto tg_xxyzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 540); 

                auto tg_xxyzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 541); 

                auto tg_xxyzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 542); 

                auto tg_xxyzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 543); 

                auto tg_xxyzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 544); 

                auto tg_xxyzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 545); 

                auto tg_xxyzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 546); 

                auto tg_xxyzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 547); 

                auto tg_xxyzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 548); 

                auto tg_xxyzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 549); 

                auto tg_xxyzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 550); 

                auto tg_xxyzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 551); 

                // Batch of Integrals (460,552)

                #pragma omp simd aligned(fxn, fza, tg_xxyyyyz_xxyyzz_0, tg_xxyyyyz_xxyzzz_0, \
                                         tg_xxyyyyz_xxzzzz_0, tg_xxyyyyz_xyyyyy_0, tg_xxyyyyz_xyyyyz_0, tg_xxyyyyz_xyyyzz_0, \
                                         tg_xxyyyyz_xyyzzz_0, tg_xxyyyyz_xyzzzz_0, tg_xxyyyyz_xzzzzz_0, tg_xxyyyyz_yyyyyy_0, \
                                         tg_xxyyyyz_yyyyyz_0, tg_xxyyyyz_yyyyzz_0, tg_xxyyyyz_yyyzzz_0, tg_xxyyyyz_yyzzzz_0, \
                                         tg_xxyyyyz_yzzzzz_0, tg_xxyyyyz_zzzzzz_0, tg_xxyyyzz_xxxxxx_0, tg_xxyyyzz_xxxxxy_0, \
                                         tg_xxyyyzz_xxxxxz_0, tg_xxyyyzz_xxxxyy_0, tg_xxyyyzz_xxxxyz_0, tg_xxyyyzz_xxxxzz_0, \
                                         tg_xxyyyzz_xxxyyy_0, tg_xxyyyzz_xxxyyz_0, tg_xxyyyzz_xxxyzz_0, tg_xxyyyzz_xxxzzz_0, \
                                         tg_xxyyyzz_xxyyyy_0, tg_xxyyyzz_xxyyyz_0, tg_xxyyyzz_xxyyzz_0, tg_xxyyyzz_xxyzzz_0, \
                                         tg_xxyyyzz_xxzzzz_0, tg_xxyyyzz_xyyyyy_0, tg_xxyyyzz_xyyyyz_0, tg_xxyyyzz_xyyyzz_0, \
                                         tg_xxyyyzz_xyyzzz_0, tg_xxyyyzz_xyzzzz_0, tg_xxyyyzz_xzzzzz_0, tg_xxyyyzz_yyyyyy_0, \
                                         tg_xxyyyzz_yyyyyz_0, tg_xxyyyzz_yyyyzz_0, tg_xxyyyzz_yyyzzz_0, tg_xxyyyzz_yyzzzz_0, \
                                         tg_xxyyyzz_yzzzzz_0, tg_xxyyyzz_zzzzzz_0, tg_xxyyzzz_xxxxxx_0, tg_xxyyzzz_xxxxxy_0, \
                                         tg_xxyyzzz_xxxxxz_0, tg_xxyyzzz_xxxxyy_0, tg_xxyyzzz_xxxxyz_0, tg_xxyyzzz_xxxxzz_0, \
                                         tg_xxyyzzz_xxxyyy_0, tg_xxyyzzz_xxxyyz_0, tg_xxyyzzz_xxxyzz_0, tg_xxyyzzz_xxxzzz_0, \
                                         tg_xxyyzzz_xxyyyy_0, tg_xxyyzzz_xxyyyz_0, tg_xxyyzzz_xxyyzz_0, tg_xxyyzzz_xxyzzz_0, \
                                         tg_xxyyzzz_xxzzzz_0, tg_xxyyzzz_xyyyyy_0, tg_xxyyzzz_xyyyyz_0, tg_xxyyzzz_xyyyzz_0, \
                                         tg_xxyyzzz_xyyzzz_0, tg_xxyyzzz_xyzzzz_0, tg_xxyyzzz_xzzzzz_0, tg_xxyyzzz_yyyyyy_0, \
                                         tg_xxyyzzz_yyyyyz_0, tg_xxyyzzz_yyyyzz_0, tg_xxyyzzz_yyyzzz_0, tg_xxyyzzz_yyzzzz_0, \
                                         tg_xxyyzzz_yzzzzz_0, tg_xxyyzzz_zzzzzz_0, tg_xxyzzzz_xxxxxx_0, tg_xxyzzzz_xxxxxy_0, \
                                         tg_xxyzzzz_xxxxxz_0, tg_xxyzzzz_xxxxyy_0, tg_xxyzzzz_xxxxyz_0, tg_xxyzzzz_xxxxzz_0, \
                                         tg_xxyzzzz_xxxyyy_0, tg_xxyzzzz_xxxyyz_0, tg_xxyzzzz_xxxyzz_0, tg_xxyzzzz_xxxzzz_0, \
                                         tg_xxyzzzz_xxyyyy_0, tg_xxyzzzz_xxyyyz_0, tg_xxyzzzz_xxyyzz_0, tg_xxyzzzz_xxyzzz_0, \
                                         tg_xxyzzzz_xxzzzz_0, tg_xxyzzzz_xyyyyy_0, tg_xxyzzzz_xyyyyz_0, tg_xxyzzzz_xyyyzz_0, \
                                         tg_xxyzzzz_xyyzzz_0, tg_xxyzzzz_xyzzzz_0, tg_xyyyyz_xxyyzz_0, tg_xyyyyz_xxyyzz_1, \
                                         tg_xyyyyz_xxyzzz_0, tg_xyyyyz_xxyzzz_1, tg_xyyyyz_xxzzzz_0, tg_xyyyyz_xxzzzz_1, \
                                         tg_xyyyyz_xyyyyy_0, tg_xyyyyz_xyyyyy_1, tg_xyyyyz_xyyyyz_0, tg_xyyyyz_xyyyyz_1, \
                                         tg_xyyyyz_xyyyzz_0, tg_xyyyyz_xyyyzz_1, tg_xyyyyz_xyyzz_1, tg_xyyyyz_xyyzzz_0, \
                                         tg_xyyyyz_xyyzzz_1, tg_xyyyyz_xyzzz_1, tg_xyyyyz_xyzzzz_0, tg_xyyyyz_xyzzzz_1, \
                                         tg_xyyyyz_xzzzz_1, tg_xyyyyz_xzzzzz_0, tg_xyyyyz_xzzzzz_1, tg_xyyyyz_yyyyy_1, \
                                         tg_xyyyyz_yyyyyy_0, tg_xyyyyz_yyyyyy_1, tg_xyyyyz_yyyyyz_0, tg_xyyyyz_yyyyyz_1, \
                                         tg_xyyyyz_yyyyz_1, tg_xyyyyz_yyyyzz_0, tg_xyyyyz_yyyyzz_1, tg_xyyyyz_yyyzz_1, \
                                         tg_xyyyyz_yyyzzz_0, tg_xyyyyz_yyyzzz_1, tg_xyyyyz_yyzzz_1, tg_xyyyyz_yyzzzz_0, \
                                         tg_xyyyyz_yyzzzz_1, tg_xyyyyz_yzzzz_1, tg_xyyyyz_yzzzzz_0, tg_xyyyyz_yzzzzz_1, \
                                         tg_xyyyyz_zzzzz_1, tg_xyyyyz_zzzzzz_0, tg_xyyyyz_zzzzzz_1, tg_xyyyzz_xxxxx_1, \
                                         tg_xyyyzz_xxxxxx_0, tg_xyyyzz_xxxxxx_1, tg_xyyyzz_xxxxxy_0, tg_xyyyzz_xxxxxy_1, \
                                         tg_xyyyzz_xxxxxz_0, tg_xyyyzz_xxxxxz_1, tg_xyyyzz_xxxxy_1, tg_xyyyzz_xxxxyy_0, \
                                         tg_xyyyzz_xxxxyy_1, tg_xyyyzz_xxxxyz_0, tg_xyyyzz_xxxxyz_1, tg_xyyyzz_xxxxz_1, \
                                         tg_xyyyzz_xxxxzz_0, tg_xyyyzz_xxxxzz_1, tg_xyyyzz_xxxyy_1, tg_xyyyzz_xxxyyy_0, \
                                         tg_xyyyzz_xxxyyy_1, tg_xyyyzz_xxxyyz_0, tg_xyyyzz_xxxyyz_1, tg_xyyyzz_xxxyz_1, \
                                         tg_xyyyzz_xxxyzz_0, tg_xyyyzz_xxxyzz_1, tg_xyyyzz_xxxzz_1, tg_xyyyzz_xxxzzz_0, \
                                         tg_xyyyzz_xxxzzz_1, tg_xyyyzz_xxyyy_1, tg_xyyyzz_xxyyyy_0, tg_xyyyzz_xxyyyy_1, \
                                         tg_xyyyzz_xxyyyz_0, tg_xyyyzz_xxyyyz_1, tg_xyyyzz_xxyyz_1, tg_xyyyzz_xxyyzz_0, \
                                         tg_xyyyzz_xxyyzz_1, tg_xyyyzz_xxyzz_1, tg_xyyyzz_xxyzzz_0, tg_xyyyzz_xxyzzz_1, \
                                         tg_xyyyzz_xxzzz_1, tg_xyyyzz_xxzzzz_0, tg_xyyyzz_xxzzzz_1, tg_xyyyzz_xyyyy_1, \
                                         tg_xyyyzz_xyyyyy_0, tg_xyyyzz_xyyyyy_1, tg_xyyyzz_xyyyyz_0, tg_xyyyzz_xyyyyz_1, \
                                         tg_xyyyzz_xyyyz_1, tg_xyyyzz_xyyyzz_0, tg_xyyyzz_xyyyzz_1, tg_xyyyzz_xyyzz_1, \
                                         tg_xyyyzz_xyyzzz_0, tg_xyyyzz_xyyzzz_1, tg_xyyyzz_xyzzz_1, tg_xyyyzz_xyzzzz_0, \
                                         tg_xyyyzz_xyzzzz_1, tg_xyyyzz_xzzzz_1, tg_xyyyzz_xzzzzz_0, tg_xyyyzz_xzzzzz_1, \
                                         tg_xyyyzz_yyyyy_1, tg_xyyyzz_yyyyyy_0, tg_xyyyzz_yyyyyy_1, tg_xyyyzz_yyyyyz_0, \
                                         tg_xyyyzz_yyyyyz_1, tg_xyyyzz_yyyyz_1, tg_xyyyzz_yyyyzz_0, tg_xyyyzz_yyyyzz_1, \
                                         tg_xyyyzz_yyyzz_1, tg_xyyyzz_yyyzzz_0, tg_xyyyzz_yyyzzz_1, tg_xyyyzz_yyzzz_1, \
                                         tg_xyyyzz_yyzzzz_0, tg_xyyyzz_yyzzzz_1, tg_xyyyzz_yzzzz_1, tg_xyyyzz_yzzzzz_0, \
                                         tg_xyyyzz_yzzzzz_1, tg_xyyyzz_zzzzz_1, tg_xyyyzz_zzzzzz_0, tg_xyyyzz_zzzzzz_1, \
                                         tg_xyyzzz_xxxxx_1, tg_xyyzzz_xxxxxx_0, tg_xyyzzz_xxxxxx_1, tg_xyyzzz_xxxxxy_0, \
                                         tg_xyyzzz_xxxxxy_1, tg_xyyzzz_xxxxxz_0, tg_xyyzzz_xxxxxz_1, tg_xyyzzz_xxxxy_1, \
                                         tg_xyyzzz_xxxxyy_0, tg_xyyzzz_xxxxyy_1, tg_xyyzzz_xxxxyz_0, tg_xyyzzz_xxxxyz_1, \
                                         tg_xyyzzz_xxxxz_1, tg_xyyzzz_xxxxzz_0, tg_xyyzzz_xxxxzz_1, tg_xyyzzz_xxxyy_1, \
                                         tg_xyyzzz_xxxyyy_0, tg_xyyzzz_xxxyyy_1, tg_xyyzzz_xxxyyz_0, tg_xyyzzz_xxxyyz_1, \
                                         tg_xyyzzz_xxxyz_1, tg_xyyzzz_xxxyzz_0, tg_xyyzzz_xxxyzz_1, tg_xyyzzz_xxxzz_1, \
                                         tg_xyyzzz_xxxzzz_0, tg_xyyzzz_xxxzzz_1, tg_xyyzzz_xxyyy_1, tg_xyyzzz_xxyyyy_0, \
                                         tg_xyyzzz_xxyyyy_1, tg_xyyzzz_xxyyyz_0, tg_xyyzzz_xxyyyz_1, tg_xyyzzz_xxyyz_1, \
                                         tg_xyyzzz_xxyyzz_0, tg_xyyzzz_xxyyzz_1, tg_xyyzzz_xxyzz_1, tg_xyyzzz_xxyzzz_0, \
                                         tg_xyyzzz_xxyzzz_1, tg_xyyzzz_xxzzz_1, tg_xyyzzz_xxzzzz_0, tg_xyyzzz_xxzzzz_1, \
                                         tg_xyyzzz_xyyyy_1, tg_xyyzzz_xyyyyy_0, tg_xyyzzz_xyyyyy_1, tg_xyyzzz_xyyyyz_0, \
                                         tg_xyyzzz_xyyyyz_1, tg_xyyzzz_xyyyz_1, tg_xyyzzz_xyyyzz_0, tg_xyyzzz_xyyyzz_1, \
                                         tg_xyyzzz_xyyzz_1, tg_xyyzzz_xyyzzz_0, tg_xyyzzz_xyyzzz_1, tg_xyyzzz_xyzzz_1, \
                                         tg_xyyzzz_xyzzzz_0, tg_xyyzzz_xyzzzz_1, tg_xyyzzz_xzzzz_1, tg_xyyzzz_xzzzzz_0, \
                                         tg_xyyzzz_xzzzzz_1, tg_xyyzzz_yyyyy_1, tg_xyyzzz_yyyyyy_0, tg_xyyzzz_yyyyyy_1, \
                                         tg_xyyzzz_yyyyyz_0, tg_xyyzzz_yyyyyz_1, tg_xyyzzz_yyyyz_1, tg_xyyzzz_yyyyzz_0, \
                                         tg_xyyzzz_yyyyzz_1, tg_xyyzzz_yyyzz_1, tg_xyyzzz_yyyzzz_0, tg_xyyzzz_yyyzzz_1, \
                                         tg_xyyzzz_yyzzz_1, tg_xyyzzz_yyzzzz_0, tg_xyyzzz_yyzzzz_1, tg_xyyzzz_yzzzz_1, \
                                         tg_xyyzzz_yzzzzz_0, tg_xyyzzz_yzzzzz_1, tg_xyyzzz_zzzzz_1, tg_xyyzzz_zzzzzz_0, \
                                         tg_xyyzzz_zzzzzz_1, tg_xyzzzz_xxxxx_1, tg_xyzzzz_xxxxxx_0, tg_xyzzzz_xxxxxx_1, \
                                         tg_xyzzzz_xxxxxy_0, tg_xyzzzz_xxxxxy_1, tg_xyzzzz_xxxxxz_0, tg_xyzzzz_xxxxxz_1, \
                                         tg_xyzzzz_xxxxy_1, tg_xyzzzz_xxxxyy_0, tg_xyzzzz_xxxxyy_1, tg_xyzzzz_xxxxyz_0, \
                                         tg_xyzzzz_xxxxyz_1, tg_xyzzzz_xxxxz_1, tg_xyzzzz_xxxxzz_0, tg_xyzzzz_xxxxzz_1, \
                                         tg_xyzzzz_xxxyy_1, tg_xyzzzz_xxxyyy_0, tg_xyzzzz_xxxyyy_1, tg_xyzzzz_xxxyyz_0, \
                                         tg_xyzzzz_xxxyyz_1, tg_xyzzzz_xxxyz_1, tg_xyzzzz_xxxyzz_0, tg_xyzzzz_xxxyzz_1, \
                                         tg_xyzzzz_xxxzz_1, tg_xyzzzz_xxxzzz_0, tg_xyzzzz_xxxzzz_1, tg_xyzzzz_xxyyy_1, \
                                         tg_xyzzzz_xxyyyy_0, tg_xyzzzz_xxyyyy_1, tg_xyzzzz_xxyyyz_0, tg_xyzzzz_xxyyyz_1, \
                                         tg_xyzzzz_xxyyz_1, tg_xyzzzz_xxyyzz_0, tg_xyzzzz_xxyyzz_1, tg_xyzzzz_xxyzz_1, \
                                         tg_xyzzzz_xxyzzz_0, tg_xyzzzz_xxyzzz_1, tg_xyzzzz_xxzzz_1, tg_xyzzzz_xxzzzz_0, \
                                         tg_xyzzzz_xxzzzz_1, tg_xyzzzz_xyyyy_1, tg_xyzzzz_xyyyyy_0, tg_xyzzzz_xyyyyy_1, \
                                         tg_xyzzzz_xyyyyz_0, tg_xyzzzz_xyyyyz_1, tg_xyzzzz_xyyyz_1, tg_xyzzzz_xyyyzz_0, \
                                         tg_xyzzzz_xyyyzz_1, tg_xyzzzz_xyyzz_1, tg_xyzzzz_xyyzzz_0, tg_xyzzzz_xyyzzz_1, \
                                         tg_xyzzzz_xyzzz_1, tg_xyzzzz_xyzzzz_0, tg_xyzzzz_xyzzzz_1, tg_xyzzzz_xzzzz_1, \
                                         tg_xyzzzz_yyyyy_1, tg_xyzzzz_yyyyz_1, tg_xyzzzz_yyyzz_1, tg_xyzzzz_yyzzz_1, \
                                         tg_xyzzzz_yzzzz_1, tg_yyyyz_xxyyzz_0, tg_yyyyz_xxyyzz_1, tg_yyyyz_xxyzzz_0, \
                                         tg_yyyyz_xxyzzz_1, tg_yyyyz_xxzzzz_0, tg_yyyyz_xxzzzz_1, tg_yyyyz_xyyyyy_0, \
                                         tg_yyyyz_xyyyyy_1, tg_yyyyz_xyyyyz_0, tg_yyyyz_xyyyyz_1, tg_yyyyz_xyyyzz_0, \
                                         tg_yyyyz_xyyyzz_1, tg_yyyyz_xyyzzz_0, tg_yyyyz_xyyzzz_1, tg_yyyyz_xyzzzz_0, \
                                         tg_yyyyz_xyzzzz_1, tg_yyyyz_xzzzzz_0, tg_yyyyz_xzzzzz_1, tg_yyyyz_yyyyyy_0, \
                                         tg_yyyyz_yyyyyy_1, tg_yyyyz_yyyyyz_0, tg_yyyyz_yyyyyz_1, tg_yyyyz_yyyyzz_0, \
                                         tg_yyyyz_yyyyzz_1, tg_yyyyz_yyyzzz_0, tg_yyyyz_yyyzzz_1, tg_yyyyz_yyzzzz_0, \
                                         tg_yyyyz_yyzzzz_1, tg_yyyyz_yzzzzz_0, tg_yyyyz_yzzzzz_1, tg_yyyyz_zzzzzz_0, \
                                         tg_yyyyz_zzzzzz_1, tg_yyyzz_xxxxxx_0, tg_yyyzz_xxxxxx_1, tg_yyyzz_xxxxxy_0, \
                                         tg_yyyzz_xxxxxy_1, tg_yyyzz_xxxxxz_0, tg_yyyzz_xxxxxz_1, tg_yyyzz_xxxxyy_0, \
                                         tg_yyyzz_xxxxyy_1, tg_yyyzz_xxxxyz_0, tg_yyyzz_xxxxyz_1, tg_yyyzz_xxxxzz_0, \
                                         tg_yyyzz_xxxxzz_1, tg_yyyzz_xxxyyy_0, tg_yyyzz_xxxyyy_1, tg_yyyzz_xxxyyz_0, \
                                         tg_yyyzz_xxxyyz_1, tg_yyyzz_xxxyzz_0, tg_yyyzz_xxxyzz_1, tg_yyyzz_xxxzzz_0, \
                                         tg_yyyzz_xxxzzz_1, tg_yyyzz_xxyyyy_0, tg_yyyzz_xxyyyy_1, tg_yyyzz_xxyyyz_0, \
                                         tg_yyyzz_xxyyyz_1, tg_yyyzz_xxyyzz_0, tg_yyyzz_xxyyzz_1, tg_yyyzz_xxyzzz_0, \
                                         tg_yyyzz_xxyzzz_1, tg_yyyzz_xxzzzz_0, tg_yyyzz_xxzzzz_1, tg_yyyzz_xyyyyy_0, \
                                         tg_yyyzz_xyyyyy_1, tg_yyyzz_xyyyyz_0, tg_yyyzz_xyyyyz_1, tg_yyyzz_xyyyzz_0, \
                                         tg_yyyzz_xyyyzz_1, tg_yyyzz_xyyzzz_0, tg_yyyzz_xyyzzz_1, tg_yyyzz_xyzzzz_0, \
                                         tg_yyyzz_xyzzzz_1, tg_yyyzz_xzzzzz_0, tg_yyyzz_xzzzzz_1, tg_yyyzz_yyyyyy_0, \
                                         tg_yyyzz_yyyyyy_1, tg_yyyzz_yyyyyz_0, tg_yyyzz_yyyyyz_1, tg_yyyzz_yyyyzz_0, \
                                         tg_yyyzz_yyyyzz_1, tg_yyyzz_yyyzzz_0, tg_yyyzz_yyyzzz_1, tg_yyyzz_yyzzzz_0, \
                                         tg_yyyzz_yyzzzz_1, tg_yyyzz_yzzzzz_0, tg_yyyzz_yzzzzz_1, tg_yyyzz_zzzzzz_0, \
                                         tg_yyyzz_zzzzzz_1, tg_yyzzz_xxxxxx_0, tg_yyzzz_xxxxxx_1, tg_yyzzz_xxxxxy_0, \
                                         tg_yyzzz_xxxxxy_1, tg_yyzzz_xxxxxz_0, tg_yyzzz_xxxxxz_1, tg_yyzzz_xxxxyy_0, \
                                         tg_yyzzz_xxxxyy_1, tg_yyzzz_xxxxyz_0, tg_yyzzz_xxxxyz_1, tg_yyzzz_xxxxzz_0, \
                                         tg_yyzzz_xxxxzz_1, tg_yyzzz_xxxyyy_0, tg_yyzzz_xxxyyy_1, tg_yyzzz_xxxyyz_0, \
                                         tg_yyzzz_xxxyyz_1, tg_yyzzz_xxxyzz_0, tg_yyzzz_xxxyzz_1, tg_yyzzz_xxxzzz_0, \
                                         tg_yyzzz_xxxzzz_1, tg_yyzzz_xxyyyy_0, tg_yyzzz_xxyyyy_1, tg_yyzzz_xxyyyz_0, \
                                         tg_yyzzz_xxyyyz_1, tg_yyzzz_xxyyzz_0, tg_yyzzz_xxyyzz_1, tg_yyzzz_xxyzzz_0, \
                                         tg_yyzzz_xxyzzz_1, tg_yyzzz_xxzzzz_0, tg_yyzzz_xxzzzz_1, tg_yyzzz_xyyyyy_0, \
                                         tg_yyzzz_xyyyyy_1, tg_yyzzz_xyyyyz_0, tg_yyzzz_xyyyyz_1, tg_yyzzz_xyyyzz_0, \
                                         tg_yyzzz_xyyyzz_1, tg_yyzzz_xyyzzz_0, tg_yyzzz_xyyzzz_1, tg_yyzzz_xyzzzz_0, \
                                         tg_yyzzz_xyzzzz_1, tg_yyzzz_xzzzzz_0, tg_yyzzz_xzzzzz_1, tg_yyzzz_yyyyyy_0, \
                                         tg_yyzzz_yyyyyy_1, tg_yyzzz_yyyyyz_0, tg_yyzzz_yyyyyz_1, tg_yyzzz_yyyyzz_0, \
                                         tg_yyzzz_yyyyzz_1, tg_yyzzz_yyyzzz_0, tg_yyzzz_yyyzzz_1, tg_yyzzz_yyzzzz_0, \
                                         tg_yyzzz_yyzzzz_1, tg_yyzzz_yzzzzz_0, tg_yyzzz_yzzzzz_1, tg_yyzzz_zzzzzz_0, \
                                         tg_yyzzz_zzzzzz_1, tg_yzzzz_xxxxxx_0, tg_yzzzz_xxxxxx_1, tg_yzzzz_xxxxxy_0, \
                                         tg_yzzzz_xxxxxy_1, tg_yzzzz_xxxxxz_0, tg_yzzzz_xxxxxz_1, tg_yzzzz_xxxxyy_0, \
                                         tg_yzzzz_xxxxyy_1, tg_yzzzz_xxxxyz_0, tg_yzzzz_xxxxyz_1, tg_yzzzz_xxxxzz_0, \
                                         tg_yzzzz_xxxxzz_1, tg_yzzzz_xxxyyy_0, tg_yzzzz_xxxyyy_1, tg_yzzzz_xxxyyz_0, \
                                         tg_yzzzz_xxxyyz_1, tg_yzzzz_xxxyzz_0, tg_yzzzz_xxxyzz_1, tg_yzzzz_xxxzzz_0, \
                                         tg_yzzzz_xxxzzz_1, tg_yzzzz_xxyyyy_0, tg_yzzzz_xxyyyy_1, tg_yzzzz_xxyyyz_0, \
                                         tg_yzzzz_xxyyyz_1, tg_yzzzz_xxyyzz_0, tg_yzzzz_xxyyzz_1, tg_yzzzz_xxyzzz_0, \
                                         tg_yzzzz_xxyzzz_1, tg_yzzzz_xxzzzz_0, tg_yzzzz_xxzzzz_1, tg_yzzzz_xyyyyy_0, \
                                         tg_yzzzz_xyyyyy_1, tg_yzzzz_xyyyyz_0, tg_yzzzz_xyyyyz_1, tg_yzzzz_xyyyzz_0, \
                                         tg_yzzzz_xyyyzz_1, tg_yzzzz_xyyzzz_0, tg_yzzzz_xyyzzz_1, tg_yzzzz_xyzzzz_0, \
                                         tg_yzzzz_xyzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyyyz_xxyyzz_0[j] = pb_x * tg_xyyyyz_xxyyzz_0[j] + fr * tg_xyyyyz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyyzz_0[j] - tg_yyyyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyyzz_1[j];

                    tg_xxyyyyz_xxyzzz_0[j] = pb_x * tg_xyyyyz_xxyzzz_0[j] + fr * tg_xyyyyz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyzzz_0[j] - tg_yyyyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyzzz_1[j];

                    tg_xxyyyyz_xxzzzz_0[j] = pb_x * tg_xyyyyz_xxzzzz_0[j] + fr * tg_xyyyyz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxzzzz_0[j] - tg_yyyyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xzzzz_1[j];

                    tg_xxyyyyz_xyyyyy_0[j] = pb_x * tg_xyyyyz_xyyyyy_0[j] + fr * tg_xyyyyz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyyyy_0[j] - tg_yyyyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyyyy_1[j];

                    tg_xxyyyyz_xyyyyz_0[j] = pb_x * tg_xyyyyz_xyyyyz_0[j] + fr * tg_xyyyyz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyyyz_0[j] - tg_yyyyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyyyz_1[j];

                    tg_xxyyyyz_xyyyzz_0[j] = pb_x * tg_xyyyyz_xyyyzz_0[j] + fr * tg_xyyyyz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyyzz_0[j] - tg_yyyyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyyzz_1[j];

                    tg_xxyyyyz_xyyzzz_0[j] = pb_x * tg_xyyyyz_xyyzzz_0[j] + fr * tg_xyyyyz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyzzz_0[j] - tg_yyyyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyzzz_1[j];

                    tg_xxyyyyz_xyzzzz_0[j] = pb_x * tg_xyyyyz_xyzzzz_0[j] + fr * tg_xyyyyz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyzzzz_0[j] - tg_yyyyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yzzzz_1[j];

                    tg_xxyyyyz_xzzzzz_0[j] = pb_x * tg_xyyyyz_xzzzzz_0[j] + fr * tg_xyyyyz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xzzzzz_0[j] - tg_yyyyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_zzzzz_1[j];

                    tg_xxyyyyz_yyyyyy_0[j] = pb_x * tg_xyyyyz_yyyyyy_0[j] + fr * tg_xyyyyz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyyyy_0[j] - tg_yyyyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyyyz_yyyyyz_0[j] = pb_x * tg_xyyyyz_yyyyyz_0[j] + fr * tg_xyyyyz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyyyz_0[j] - tg_yyyyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyyyz_yyyyzz_0[j] = pb_x * tg_xyyyyz_yyyyzz_0[j] + fr * tg_xyyyyz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyyzz_0[j] - tg_yyyyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyyyz_yyyzzz_0[j] = pb_x * tg_xyyyyz_yyyzzz_0[j] + fr * tg_xyyyyz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyzzz_0[j] - tg_yyyyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_yyzzzz_0[j] = pb_x * tg_xyyyyz_yyzzzz_0[j] + fr * tg_xyyyyz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyzzzz_0[j] - tg_yyyyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_yzzzzz_0[j] = pb_x * tg_xyyyyz_yzzzzz_0[j] + fr * tg_xyyyyz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yzzzzz_0[j] - tg_yyyyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_zzzzzz_0[j] = pb_x * tg_xyyyyz_zzzzzz_0[j] + fr * tg_xyyyyz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_zzzzzz_0[j] - tg_yyyyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_xxxxxx_0[j] = pb_x * tg_xyyyzz_xxxxxx_0[j] + fr * tg_xyyyzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxxx_0[j] - tg_yyyzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyyzz_xxxxx_1[j];

                    tg_xxyyyzz_xxxxxy_0[j] = pb_x * tg_xyyyzz_xxxxxy_0[j] + fr * tg_xyyyzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxxy_0[j] - tg_yyyzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyzz_xxxxy_1[j];

                    tg_xxyyyzz_xxxxxz_0[j] = pb_x * tg_xyyyzz_xxxxxz_0[j] + fr * tg_xyyyzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxxz_0[j] - tg_yyyzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyyzz_xxxxz_1[j];

                    tg_xxyyyzz_xxxxyy_0[j] = pb_x * tg_xyyyzz_xxxxyy_0[j] + fr * tg_xyyyzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxyy_0[j] - tg_yyyzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzz_xxxyy_1[j];

                    tg_xxyyyzz_xxxxyz_0[j] = pb_x * tg_xyyyzz_xxxxyz_0[j] + fr * tg_xyyyzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxyz_0[j] - tg_yyyzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzz_xxxyz_1[j];

                    tg_xxyyyzz_xxxxzz_0[j] = pb_x * tg_xyyyzz_xxxxzz_0[j] + fr * tg_xyyyzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxxzz_0[j] - tg_yyyzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzz_xxxzz_1[j];

                    tg_xxyyyzz_xxxyyy_0[j] = pb_x * tg_xyyyzz_xxxyyy_0[j] + fr * tg_xyyyzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxyyy_0[j] - tg_yyyzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxyyy_1[j];

                    tg_xxyyyzz_xxxyyz_0[j] = pb_x * tg_xyyyzz_xxxyyz_0[j] + fr * tg_xyyyzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxyyz_0[j] - tg_yyyzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxyyz_1[j];

                    tg_xxyyyzz_xxxyzz_0[j] = pb_x * tg_xyyyzz_xxxyzz_0[j] + fr * tg_xyyyzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxyzz_0[j] - tg_yyyzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxyzz_1[j];

                    tg_xxyyyzz_xxxzzz_0[j] = pb_x * tg_xyyyzz_xxxzzz_0[j] + fr * tg_xyyyzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxzzz_0[j] - tg_yyyzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxzzz_1[j];

                    tg_xxyyyzz_xxyyyy_0[j] = pb_x * tg_xyyyzz_xxyyyy_0[j] + fr * tg_xyyyzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyyyy_0[j] - tg_yyyzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyyyy_1[j];

                    tg_xxyyyzz_xxyyyz_0[j] = pb_x * tg_xyyyzz_xxyyyz_0[j] + fr * tg_xyyyzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyyyz_0[j] - tg_yyyzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyyyz_1[j];

                    tg_xxyyyzz_xxyyzz_0[j] = pb_x * tg_xyyyzz_xxyyzz_0[j] + fr * tg_xyyyzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyyzz_0[j] - tg_yyyzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyyzz_1[j];

                    tg_xxyyyzz_xxyzzz_0[j] = pb_x * tg_xyyyzz_xxyzzz_0[j] + fr * tg_xyyyzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyzzz_0[j] - tg_yyyzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyzzz_1[j];

                    tg_xxyyyzz_xxzzzz_0[j] = pb_x * tg_xyyyzz_xxzzzz_0[j] + fr * tg_xyyyzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxzzzz_0[j] - tg_yyyzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xzzzz_1[j];

                    tg_xxyyyzz_xyyyyy_0[j] = pb_x * tg_xyyyzz_xyyyyy_0[j] + fr * tg_xyyyzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyyyy_0[j] - tg_yyyzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyyyy_1[j];

                    tg_xxyyyzz_xyyyyz_0[j] = pb_x * tg_xyyyzz_xyyyyz_0[j] + fr * tg_xyyyzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyyyz_0[j] - tg_yyyzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyyyz_1[j];

                    tg_xxyyyzz_xyyyzz_0[j] = pb_x * tg_xyyyzz_xyyyzz_0[j] + fr * tg_xyyyzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyyzz_0[j] - tg_yyyzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyyzz_1[j];

                    tg_xxyyyzz_xyyzzz_0[j] = pb_x * tg_xyyyzz_xyyzzz_0[j] + fr * tg_xyyyzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyzzz_0[j] - tg_yyyzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyzzz_1[j];

                    tg_xxyyyzz_xyzzzz_0[j] = pb_x * tg_xyyyzz_xyzzzz_0[j] + fr * tg_xyyyzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyzzzz_0[j] - tg_yyyzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yzzzz_1[j];

                    tg_xxyyyzz_xzzzzz_0[j] = pb_x * tg_xyyyzz_xzzzzz_0[j] + fr * tg_xyyyzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xzzzzz_0[j] - tg_yyyzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_zzzzz_1[j];

                    tg_xxyyyzz_yyyyyy_0[j] = pb_x * tg_xyyyzz_yyyyyy_0[j] + fr * tg_xyyyzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyyyy_0[j] - tg_yyyzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyyzz_yyyyyz_0[j] = pb_x * tg_xyyyzz_yyyyyz_0[j] + fr * tg_xyyyzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyyyz_0[j] - tg_yyyzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyyzz_yyyyzz_0[j] = pb_x * tg_xyyyzz_yyyyzz_0[j] + fr * tg_xyyyzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyyzz_0[j] - tg_yyyzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyyzz_yyyzzz_0[j] = pb_x * tg_xyyyzz_yyyzzz_0[j] + fr * tg_xyyyzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyzzz_0[j] - tg_yyyzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_yyzzzz_0[j] = pb_x * tg_xyyyzz_yyzzzz_0[j] + fr * tg_xyyyzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyzzzz_0[j] - tg_yyyzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_yzzzzz_0[j] = pb_x * tg_xyyyzz_yzzzzz_0[j] + fr * tg_xyyyzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yzzzzz_0[j] - tg_yyyzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_zzzzzz_0[j] = pb_x * tg_xyyyzz_zzzzzz_0[j] + fr * tg_xyyyzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_zzzzzz_0[j] - tg_yyyzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_xxxxxx_0[j] = pb_x * tg_xyyzzz_xxxxxx_0[j] + fr * tg_xyyzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxxx_0[j] - tg_yyzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyzzz_xxxxx_1[j];

                    tg_xxyyzzz_xxxxxy_0[j] = pb_x * tg_xyyzzz_xxxxxy_0[j] + fr * tg_xyyzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxxy_0[j] - tg_yyzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyzzz_xxxxy_1[j];

                    tg_xxyyzzz_xxxxxz_0[j] = pb_x * tg_xyyzzz_xxxxxz_0[j] + fr * tg_xyyzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxxz_0[j] - tg_yyzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyzzz_xxxxz_1[j];

                    tg_xxyyzzz_xxxxyy_0[j] = pb_x * tg_xyyzzz_xxxxyy_0[j] + fr * tg_xyyzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxyy_0[j] - tg_yyzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzz_xxxyy_1[j];

                    tg_xxyyzzz_xxxxyz_0[j] = pb_x * tg_xyyzzz_xxxxyz_0[j] + fr * tg_xyyzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxyz_0[j] - tg_yyzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzz_xxxyz_1[j];

                    tg_xxyyzzz_xxxxzz_0[j] = pb_x * tg_xyyzzz_xxxxzz_0[j] + fr * tg_xyyzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxxzz_0[j] - tg_yyzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzz_xxxzz_1[j];

                    tg_xxyyzzz_xxxyyy_0[j] = pb_x * tg_xyyzzz_xxxyyy_0[j] + fr * tg_xyyzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxyyy_0[j] - tg_yyzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxyyy_1[j];

                    tg_xxyyzzz_xxxyyz_0[j] = pb_x * tg_xyyzzz_xxxyyz_0[j] + fr * tg_xyyzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxyyz_0[j] - tg_yyzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxyyz_1[j];

                    tg_xxyyzzz_xxxyzz_0[j] = pb_x * tg_xyyzzz_xxxyzz_0[j] + fr * tg_xyyzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxyzz_0[j] - tg_yyzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxyzz_1[j];

                    tg_xxyyzzz_xxxzzz_0[j] = pb_x * tg_xyyzzz_xxxzzz_0[j] + fr * tg_xyyzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxzzz_0[j] - tg_yyzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxzzz_1[j];

                    tg_xxyyzzz_xxyyyy_0[j] = pb_x * tg_xyyzzz_xxyyyy_0[j] + fr * tg_xyyzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyyyy_0[j] - tg_yyzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyyyy_1[j];

                    tg_xxyyzzz_xxyyyz_0[j] = pb_x * tg_xyyzzz_xxyyyz_0[j] + fr * tg_xyyzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyyyz_0[j] - tg_yyzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyyyz_1[j];

                    tg_xxyyzzz_xxyyzz_0[j] = pb_x * tg_xyyzzz_xxyyzz_0[j] + fr * tg_xyyzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyyzz_0[j] - tg_yyzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyyzz_1[j];

                    tg_xxyyzzz_xxyzzz_0[j] = pb_x * tg_xyyzzz_xxyzzz_0[j] + fr * tg_xyyzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyzzz_0[j] - tg_yyzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyzzz_1[j];

                    tg_xxyyzzz_xxzzzz_0[j] = pb_x * tg_xyyzzz_xxzzzz_0[j] + fr * tg_xyyzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxzzzz_0[j] - tg_yyzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xzzzz_1[j];

                    tg_xxyyzzz_xyyyyy_0[j] = pb_x * tg_xyyzzz_xyyyyy_0[j] + fr * tg_xyyzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyyyy_0[j] - tg_yyzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyyyy_1[j];

                    tg_xxyyzzz_xyyyyz_0[j] = pb_x * tg_xyyzzz_xyyyyz_0[j] + fr * tg_xyyzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyyyz_0[j] - tg_yyzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyyyz_1[j];

                    tg_xxyyzzz_xyyyzz_0[j] = pb_x * tg_xyyzzz_xyyyzz_0[j] + fr * tg_xyyzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyyzz_0[j] - tg_yyzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyyzz_1[j];

                    tg_xxyyzzz_xyyzzz_0[j] = pb_x * tg_xyyzzz_xyyzzz_0[j] + fr * tg_xyyzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyzzz_0[j] - tg_yyzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyzzz_1[j];

                    tg_xxyyzzz_xyzzzz_0[j] = pb_x * tg_xyyzzz_xyzzzz_0[j] + fr * tg_xyyzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyzzzz_0[j] - tg_yyzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yzzzz_1[j];

                    tg_xxyyzzz_xzzzzz_0[j] = pb_x * tg_xyyzzz_xzzzzz_0[j] + fr * tg_xyyzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xzzzzz_0[j] - tg_yyzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_zzzzz_1[j];

                    tg_xxyyzzz_yyyyyy_0[j] = pb_x * tg_xyyzzz_yyyyyy_0[j] + fr * tg_xyyzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyyyy_0[j] - tg_yyzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyzzz_yyyyyz_0[j] = pb_x * tg_xyyzzz_yyyyyz_0[j] + fr * tg_xyyzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyyyz_0[j] - tg_yyzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyzzz_yyyyzz_0[j] = pb_x * tg_xyyzzz_yyyyzz_0[j] + fr * tg_xyyzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyyzz_0[j] - tg_yyzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyzzz_yyyzzz_0[j] = pb_x * tg_xyyzzz_yyyzzz_0[j] + fr * tg_xyyzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyzzz_0[j] - tg_yyzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_yyzzzz_0[j] = pb_x * tg_xyyzzz_yyzzzz_0[j] + fr * tg_xyyzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyzzzz_0[j] - tg_yyzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_yzzzzz_0[j] = pb_x * tg_xyyzzz_yzzzzz_0[j] + fr * tg_xyyzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yzzzzz_0[j] - tg_yyzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_zzzzzz_0[j] = pb_x * tg_xyyzzz_zzzzzz_0[j] + fr * tg_xyyzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_zzzzzz_0[j] - tg_yyzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_xxxxxx_0[j] = pb_x * tg_xyzzzz_xxxxxx_0[j] + fr * tg_xyzzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxxx_0[j] - tg_yzzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyzzzz_xxxxx_1[j];

                    tg_xxyzzzz_xxxxxy_0[j] = pb_x * tg_xyzzzz_xxxxxy_0[j] + fr * tg_xyzzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxxy_0[j] - tg_yzzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzzzz_xxxxy_1[j];

                    tg_xxyzzzz_xxxxxz_0[j] = pb_x * tg_xyzzzz_xxxxxz_0[j] + fr * tg_xyzzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxxz_0[j] - tg_yzzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzzzz_xxxxz_1[j];

                    tg_xxyzzzz_xxxxyy_0[j] = pb_x * tg_xyzzzz_xxxxyy_0[j] + fr * tg_xyzzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxyy_0[j] - tg_yzzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzz_xxxyy_1[j];

                    tg_xxyzzzz_xxxxyz_0[j] = pb_x * tg_xyzzzz_xxxxyz_0[j] + fr * tg_xyzzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxyz_0[j] - tg_yzzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzz_xxxyz_1[j];

                    tg_xxyzzzz_xxxxzz_0[j] = pb_x * tg_xyzzzz_xxxxzz_0[j] + fr * tg_xyzzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxxzz_0[j] - tg_yzzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzz_xxxzz_1[j];

                    tg_xxyzzzz_xxxyyy_0[j] = pb_x * tg_xyzzzz_xxxyyy_0[j] + fr * tg_xyzzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxyyy_0[j] - tg_yzzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxyyy_1[j];

                    tg_xxyzzzz_xxxyyz_0[j] = pb_x * tg_xyzzzz_xxxyyz_0[j] + fr * tg_xyzzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxyyz_0[j] - tg_yzzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxyyz_1[j];

                    tg_xxyzzzz_xxxyzz_0[j] = pb_x * tg_xyzzzz_xxxyzz_0[j] + fr * tg_xyzzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxyzz_0[j] - tg_yzzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxyzz_1[j];

                    tg_xxyzzzz_xxxzzz_0[j] = pb_x * tg_xyzzzz_xxxzzz_0[j] + fr * tg_xyzzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxzzz_0[j] - tg_yzzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxzzz_1[j];

                    tg_xxyzzzz_xxyyyy_0[j] = pb_x * tg_xyzzzz_xxyyyy_0[j] + fr * tg_xyzzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyyyy_0[j] - tg_yzzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyyyy_1[j];

                    tg_xxyzzzz_xxyyyz_0[j] = pb_x * tg_xyzzzz_xxyyyz_0[j] + fr * tg_xyzzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyyyz_0[j] - tg_yzzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyyyz_1[j];

                    tg_xxyzzzz_xxyyzz_0[j] = pb_x * tg_xyzzzz_xxyyzz_0[j] + fr * tg_xyzzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyyzz_0[j] - tg_yzzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyyzz_1[j];

                    tg_xxyzzzz_xxyzzz_0[j] = pb_x * tg_xyzzzz_xxyzzz_0[j] + fr * tg_xyzzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyzzz_0[j] - tg_yzzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyzzz_1[j];

                    tg_xxyzzzz_xxzzzz_0[j] = pb_x * tg_xyzzzz_xxzzzz_0[j] + fr * tg_xyzzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxzzzz_0[j] - tg_yzzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xzzzz_1[j];

                    tg_xxyzzzz_xyyyyy_0[j] = pb_x * tg_xyzzzz_xyyyyy_0[j] + fr * tg_xyzzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyyyy_0[j] - tg_yzzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyyyy_1[j];

                    tg_xxyzzzz_xyyyyz_0[j] = pb_x * tg_xyzzzz_xyyyyz_0[j] + fr * tg_xyzzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyyyz_0[j] - tg_yzzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyyyz_1[j];

                    tg_xxyzzzz_xyyyzz_0[j] = pb_x * tg_xyzzzz_xyyyzz_0[j] + fr * tg_xyzzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyyzz_0[j] - tg_yzzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyyzz_1[j];

                    tg_xxyzzzz_xyyzzz_0[j] = pb_x * tg_xyzzzz_xyyzzz_0[j] + fr * tg_xyzzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyzzz_0[j] - tg_yzzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyzzz_1[j];

                    tg_xxyzzzz_xyzzzz_0[j] = pb_x * tg_xyzzzz_xyzzzz_0[j] + fr * tg_xyzzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyzzzz_0[j] - tg_yzzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_552_644(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (552,644)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xyzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 552); 

                auto tg_xyzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 553); 

                auto tg_xyzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 554); 

                auto tg_xyzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 555); 

                auto tg_xyzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 556); 

                auto tg_xyzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 557); 

                auto tg_xyzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 558); 

                auto tg_xyzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 559); 

                auto tg_xzzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 560); 

                auto tg_xzzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 561); 

                auto tg_xzzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 562); 

                auto tg_xzzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 563); 

                auto tg_xzzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 564); 

                auto tg_xzzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 565); 

                auto tg_xzzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 566); 

                auto tg_xzzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 567); 

                auto tg_xzzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 568); 

                auto tg_xzzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 569); 

                auto tg_xzzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 570); 

                auto tg_xzzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 571); 

                auto tg_xzzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 572); 

                auto tg_xzzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 573); 

                auto tg_xzzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 574); 

                auto tg_xzzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 575); 

                auto tg_xzzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 576); 

                auto tg_xzzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 577); 

                auto tg_xzzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 578); 

                auto tg_xzzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 579); 

                auto tg_xzzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 580); 

                auto tg_xzzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 581); 

                auto tg_xzzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 582); 

                auto tg_xzzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 583); 

                auto tg_xzzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 584); 

                auto tg_xzzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 585); 

                auto tg_xzzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 586); 

                auto tg_xzzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 587); 

                auto tg_yyyyyy_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 588); 

                auto tg_yyyyyy_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 589); 

                auto tg_yyyyyy_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 590); 

                auto tg_yyyyyy_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 591); 

                auto tg_yyyyyy_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 592); 

                auto tg_yyyyyy_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 593); 

                auto tg_yyyyyy_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 594); 

                auto tg_yyyyyy_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 595); 

                auto tg_yyyyyy_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 596); 

                auto tg_yyyyyy_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 597); 

                auto tg_yyyyyy_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 598); 

                auto tg_yyyyyy_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 599); 

                auto tg_yyyyyy_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 600); 

                auto tg_yyyyyy_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 601); 

                auto tg_yyyyyy_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 602); 

                auto tg_yyyyyy_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 603); 

                auto tg_yyyyyy_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 604); 

                auto tg_yyyyyy_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 605); 

                auto tg_yyyyyy_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 606); 

                auto tg_yyyyyy_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 607); 

                auto tg_yyyyyy_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 608); 

                auto tg_yyyyyy_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 609); 

                auto tg_yyyyyy_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 610); 

                auto tg_yyyyyy_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 611); 

                auto tg_yyyyyy_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 612); 

                auto tg_yyyyyy_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 613); 

                auto tg_yyyyyy_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 614); 

                auto tg_yyyyyy_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 615); 

                auto tg_yyyyyz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 616); 

                auto tg_yyyyyz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 617); 

                auto tg_yyyyyz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 618); 

                auto tg_yyyyyz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 619); 

                auto tg_yyyyyz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 620); 

                auto tg_yyyyyz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 621); 

                auto tg_yyyyyz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 622); 

                auto tg_yyyyyz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 623); 

                auto tg_yyyyyz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 624); 

                auto tg_yyyyyz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 625); 

                auto tg_yyyyyz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 626); 

                auto tg_yyyyyz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 627); 

                auto tg_yyyyyz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 628); 

                auto tg_yyyyyz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 629); 

                auto tg_yyyyyz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 630); 

                auto tg_yyyyyz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 631); 

                auto tg_yyyyyz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 632); 

                auto tg_yyyyyz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 633); 

                auto tg_yyyyyz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 634); 

                auto tg_yyyyyz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 635); 

                auto tg_yyyyyz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 636); 

                auto tg_yyyyyz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 637); 

                auto tg_yyyyyz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 638); 

                auto tg_yyyyyz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 639); 

                auto tg_yyyyyz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 640); 

                auto tg_yyyyyz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 641); 

                auto tg_yyyyyz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 642); 

                auto tg_yyyyyz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 643); 

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

                auto tg_yyyyyz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 475); 

                auto tg_yyyyyz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 476); 

                auto tg_yyyyyz_yyyyy_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 477); 

                auto tg_yyyyyz_yyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 478); 

                auto tg_yyyyyz_yyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 479); 

                auto tg_yyyyyz_yyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 480); 

                auto tg_yyyyyz_yzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 481); 

                auto tg_yyyyyz_zzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 482); 

                // set up pointers to integrals

                auto tg_xxyzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 552); 

                auto tg_xxyzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 553); 

                auto tg_xxyzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 554); 

                auto tg_xxyzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 555); 

                auto tg_xxyzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 556); 

                auto tg_xxyzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 557); 

                auto tg_xxyzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 558); 

                auto tg_xxyzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 559); 

                auto tg_xxzzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 560); 

                auto tg_xxzzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 561); 

                auto tg_xxzzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 562); 

                auto tg_xxzzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 563); 

                auto tg_xxzzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 564); 

                auto tg_xxzzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 565); 

                auto tg_xxzzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 566); 

                auto tg_xxzzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 567); 

                auto tg_xxzzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 568); 

                auto tg_xxzzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 569); 

                auto tg_xxzzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 570); 

                auto tg_xxzzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 571); 

                auto tg_xxzzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 572); 

                auto tg_xxzzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 573); 

                auto tg_xxzzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 574); 

                auto tg_xxzzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 575); 

                auto tg_xxzzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 576); 

                auto tg_xxzzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 577); 

                auto tg_xxzzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 578); 

                auto tg_xxzzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 579); 

                auto tg_xxzzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 580); 

                auto tg_xxzzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 581); 

                auto tg_xxzzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 582); 

                auto tg_xxzzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 583); 

                auto tg_xxzzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 584); 

                auto tg_xxzzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 585); 

                auto tg_xxzzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 586); 

                auto tg_xxzzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 587); 

                auto tg_xyyyyyy_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 588); 

                auto tg_xyyyyyy_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 589); 

                auto tg_xyyyyyy_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 590); 

                auto tg_xyyyyyy_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 591); 

                auto tg_xyyyyyy_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 592); 

                auto tg_xyyyyyy_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 593); 

                auto tg_xyyyyyy_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 594); 

                auto tg_xyyyyyy_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 595); 

                auto tg_xyyyyyy_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 596); 

                auto tg_xyyyyyy_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 597); 

                auto tg_xyyyyyy_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 598); 

                auto tg_xyyyyyy_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 599); 

                auto tg_xyyyyyy_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 600); 

                auto tg_xyyyyyy_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 601); 

                auto tg_xyyyyyy_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 602); 

                auto tg_xyyyyyy_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 603); 

                auto tg_xyyyyyy_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 604); 

                auto tg_xyyyyyy_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 605); 

                auto tg_xyyyyyy_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 606); 

                auto tg_xyyyyyy_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 607); 

                auto tg_xyyyyyy_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 608); 

                auto tg_xyyyyyy_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 609); 

                auto tg_xyyyyyy_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 610); 

                auto tg_xyyyyyy_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 611); 

                auto tg_xyyyyyy_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 612); 

                auto tg_xyyyyyy_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 613); 

                auto tg_xyyyyyy_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 614); 

                auto tg_xyyyyyy_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 615); 

                auto tg_xyyyyyz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 616); 

                auto tg_xyyyyyz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 617); 

                auto tg_xyyyyyz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 618); 

                auto tg_xyyyyyz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 619); 

                auto tg_xyyyyyz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 620); 

                auto tg_xyyyyyz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 621); 

                auto tg_xyyyyyz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 622); 

                auto tg_xyyyyyz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 623); 

                auto tg_xyyyyyz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 624); 

                auto tg_xyyyyyz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 625); 

                auto tg_xyyyyyz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 626); 

                auto tg_xyyyyyz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 627); 

                auto tg_xyyyyyz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 628); 

                auto tg_xyyyyyz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 629); 

                auto tg_xyyyyyz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 630); 

                auto tg_xyyyyyz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 631); 

                auto tg_xyyyyyz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 632); 

                auto tg_xyyyyyz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 633); 

                auto tg_xyyyyyz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 634); 

                auto tg_xyyyyyz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 635); 

                auto tg_xyyyyyz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 636); 

                auto tg_xyyyyyz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 637); 

                auto tg_xyyyyyz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 638); 

                auto tg_xyyyyyz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 639); 

                auto tg_xyyyyyz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 640); 

                auto tg_xyyyyyz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 641); 

                auto tg_xyyyyyz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 642); 

                auto tg_xyyyyyz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 643); 

                // Batch of Integrals (552,644)

                #pragma omp simd aligned(fxn, fza, tg_xxyzzzz_xzzzzz_0, tg_xxyzzzz_yyyyyy_0, \
                                         tg_xxyzzzz_yyyyyz_0, tg_xxyzzzz_yyyyzz_0, tg_xxyzzzz_yyyzzz_0, tg_xxyzzzz_yyzzzz_0, \
                                         tg_xxyzzzz_yzzzzz_0, tg_xxyzzzz_zzzzzz_0, tg_xxzzzzz_xxxxxx_0, tg_xxzzzzz_xxxxxy_0, \
                                         tg_xxzzzzz_xxxxxz_0, tg_xxzzzzz_xxxxyy_0, tg_xxzzzzz_xxxxyz_0, tg_xxzzzzz_xxxxzz_0, \
                                         tg_xxzzzzz_xxxyyy_0, tg_xxzzzzz_xxxyyz_0, tg_xxzzzzz_xxxyzz_0, tg_xxzzzzz_xxxzzz_0, \
                                         tg_xxzzzzz_xxyyyy_0, tg_xxzzzzz_xxyyyz_0, tg_xxzzzzz_xxyyzz_0, tg_xxzzzzz_xxyzzz_0, \
                                         tg_xxzzzzz_xxzzzz_0, tg_xxzzzzz_xyyyyy_0, tg_xxzzzzz_xyyyyz_0, tg_xxzzzzz_xyyyzz_0, \
                                         tg_xxzzzzz_xyyzzz_0, tg_xxzzzzz_xyzzzz_0, tg_xxzzzzz_xzzzzz_0, tg_xxzzzzz_yyyyyy_0, \
                                         tg_xxzzzzz_yyyyyz_0, tg_xxzzzzz_yyyyzz_0, tg_xxzzzzz_yyyzzz_0, tg_xxzzzzz_yyzzzz_0, \
                                         tg_xxzzzzz_yzzzzz_0, tg_xxzzzzz_zzzzzz_0, tg_xyyyyyy_xxxxxx_0, tg_xyyyyyy_xxxxxy_0, \
                                         tg_xyyyyyy_xxxxxz_0, tg_xyyyyyy_xxxxyy_0, tg_xyyyyyy_xxxxyz_0, tg_xyyyyyy_xxxxzz_0, \
                                         tg_xyyyyyy_xxxyyy_0, tg_xyyyyyy_xxxyyz_0, tg_xyyyyyy_xxxyzz_0, tg_xyyyyyy_xxxzzz_0, \
                                         tg_xyyyyyy_xxyyyy_0, tg_xyyyyyy_xxyyyz_0, tg_xyyyyyy_xxyyzz_0, tg_xyyyyyy_xxyzzz_0, \
                                         tg_xyyyyyy_xxzzzz_0, tg_xyyyyyy_xyyyyy_0, tg_xyyyyyy_xyyyyz_0, tg_xyyyyyy_xyyyzz_0, \
                                         tg_xyyyyyy_xyyzzz_0, tg_xyyyyyy_xyzzzz_0, tg_xyyyyyy_xzzzzz_0, tg_xyyyyyy_yyyyyy_0, \
                                         tg_xyyyyyy_yyyyyz_0, tg_xyyyyyy_yyyyzz_0, tg_xyyyyyy_yyyzzz_0, tg_xyyyyyy_yyzzzz_0, \
                                         tg_xyyyyyy_yzzzzz_0, tg_xyyyyyy_zzzzzz_0, tg_xyyyyyz_xxxxxx_0, tg_xyyyyyz_xxxxxy_0, \
                                         tg_xyyyyyz_xxxxxz_0, tg_xyyyyyz_xxxxyy_0, tg_xyyyyyz_xxxxyz_0, tg_xyyyyyz_xxxxzz_0, \
                                         tg_xyyyyyz_xxxyyy_0, tg_xyyyyyz_xxxyyz_0, tg_xyyyyyz_xxxyzz_0, tg_xyyyyyz_xxxzzz_0, \
                                         tg_xyyyyyz_xxyyyy_0, tg_xyyyyyz_xxyyyz_0, tg_xyyyyyz_xxyyzz_0, tg_xyyyyyz_xxyzzz_0, \
                                         tg_xyyyyyz_xxzzzz_0, tg_xyyyyyz_xyyyyy_0, tg_xyyyyyz_xyyyyz_0, tg_xyyyyyz_xyyyzz_0, \
                                         tg_xyyyyyz_xyyzzz_0, tg_xyyyyyz_xyzzzz_0, tg_xyyyyyz_xzzzzz_0, tg_xyyyyyz_yyyyyy_0, \
                                         tg_xyyyyyz_yyyyyz_0, tg_xyyyyyz_yyyyzz_0, tg_xyyyyyz_yyyzzz_0, tg_xyyyyyz_yyzzzz_0, \
                                         tg_xyyyyyz_yzzzzz_0, tg_xyyyyyz_zzzzzz_0, tg_xyzzzz_xzzzzz_0, tg_xyzzzz_xzzzzz_1, \
                                         tg_xyzzzz_yyyyyy_0, tg_xyzzzz_yyyyyy_1, tg_xyzzzz_yyyyyz_0, tg_xyzzzz_yyyyyz_1, \
                                         tg_xyzzzz_yyyyzz_0, tg_xyzzzz_yyyyzz_1, tg_xyzzzz_yyyzzz_0, tg_xyzzzz_yyyzzz_1, \
                                         tg_xyzzzz_yyzzzz_0, tg_xyzzzz_yyzzzz_1, tg_xyzzzz_yzzzzz_0, tg_xyzzzz_yzzzzz_1, \
                                         tg_xyzzzz_zzzzz_1, tg_xyzzzz_zzzzzz_0, tg_xyzzzz_zzzzzz_1, tg_xzzzzz_xxxxx_1, \
                                         tg_xzzzzz_xxxxxx_0, tg_xzzzzz_xxxxxx_1, tg_xzzzzz_xxxxxy_0, tg_xzzzzz_xxxxxy_1, \
                                         tg_xzzzzz_xxxxxz_0, tg_xzzzzz_xxxxxz_1, tg_xzzzzz_xxxxy_1, tg_xzzzzz_xxxxyy_0, \
                                         tg_xzzzzz_xxxxyy_1, tg_xzzzzz_xxxxyz_0, tg_xzzzzz_xxxxyz_1, tg_xzzzzz_xxxxz_1, \
                                         tg_xzzzzz_xxxxzz_0, tg_xzzzzz_xxxxzz_1, tg_xzzzzz_xxxyy_1, tg_xzzzzz_xxxyyy_0, \
                                         tg_xzzzzz_xxxyyy_1, tg_xzzzzz_xxxyyz_0, tg_xzzzzz_xxxyyz_1, tg_xzzzzz_xxxyz_1, \
                                         tg_xzzzzz_xxxyzz_0, tg_xzzzzz_xxxyzz_1, tg_xzzzzz_xxxzz_1, tg_xzzzzz_xxxzzz_0, \
                                         tg_xzzzzz_xxxzzz_1, tg_xzzzzz_xxyyy_1, tg_xzzzzz_xxyyyy_0, tg_xzzzzz_xxyyyy_1, \
                                         tg_xzzzzz_xxyyyz_0, tg_xzzzzz_xxyyyz_1, tg_xzzzzz_xxyyz_1, tg_xzzzzz_xxyyzz_0, \
                                         tg_xzzzzz_xxyyzz_1, tg_xzzzzz_xxyzz_1, tg_xzzzzz_xxyzzz_0, tg_xzzzzz_xxyzzz_1, \
                                         tg_xzzzzz_xxzzz_1, tg_xzzzzz_xxzzzz_0, tg_xzzzzz_xxzzzz_1, tg_xzzzzz_xyyyy_1, \
                                         tg_xzzzzz_xyyyyy_0, tg_xzzzzz_xyyyyy_1, tg_xzzzzz_xyyyyz_0, tg_xzzzzz_xyyyyz_1, \
                                         tg_xzzzzz_xyyyz_1, tg_xzzzzz_xyyyzz_0, tg_xzzzzz_xyyyzz_1, tg_xzzzzz_xyyzz_1, \
                                         tg_xzzzzz_xyyzzz_0, tg_xzzzzz_xyyzzz_1, tg_xzzzzz_xyzzz_1, tg_xzzzzz_xyzzzz_0, \
                                         tg_xzzzzz_xyzzzz_1, tg_xzzzzz_xzzzz_1, tg_xzzzzz_xzzzzz_0, tg_xzzzzz_xzzzzz_1, \
                                         tg_xzzzzz_yyyyy_1, tg_xzzzzz_yyyyyy_0, tg_xzzzzz_yyyyyy_1, tg_xzzzzz_yyyyyz_0, \
                                         tg_xzzzzz_yyyyyz_1, tg_xzzzzz_yyyyz_1, tg_xzzzzz_yyyyzz_0, tg_xzzzzz_yyyyzz_1, \
                                         tg_xzzzzz_yyyzz_1, tg_xzzzzz_yyyzzz_0, tg_xzzzzz_yyyzzz_1, tg_xzzzzz_yyzzz_1, \
                                         tg_xzzzzz_yyzzzz_0, tg_xzzzzz_yyzzzz_1, tg_xzzzzz_yzzzz_1, tg_xzzzzz_yzzzzz_0, \
                                         tg_xzzzzz_yzzzzz_1, tg_xzzzzz_zzzzz_1, tg_xzzzzz_zzzzzz_0, tg_xzzzzz_zzzzzz_1, \
                                         tg_yyyyyy_xxxxx_1, tg_yyyyyy_xxxxxx_0, tg_yyyyyy_xxxxxx_1, tg_yyyyyy_xxxxxy_0, \
                                         tg_yyyyyy_xxxxxy_1, tg_yyyyyy_xxxxxz_0, tg_yyyyyy_xxxxxz_1, tg_yyyyyy_xxxxy_1, \
                                         tg_yyyyyy_xxxxyy_0, tg_yyyyyy_xxxxyy_1, tg_yyyyyy_xxxxyz_0, tg_yyyyyy_xxxxyz_1, \
                                         tg_yyyyyy_xxxxz_1, tg_yyyyyy_xxxxzz_0, tg_yyyyyy_xxxxzz_1, tg_yyyyyy_xxxyy_1, \
                                         tg_yyyyyy_xxxyyy_0, tg_yyyyyy_xxxyyy_1, tg_yyyyyy_xxxyyz_0, tg_yyyyyy_xxxyyz_1, \
                                         tg_yyyyyy_xxxyz_1, tg_yyyyyy_xxxyzz_0, tg_yyyyyy_xxxyzz_1, tg_yyyyyy_xxxzz_1, \
                                         tg_yyyyyy_xxxzzz_0, tg_yyyyyy_xxxzzz_1, tg_yyyyyy_xxyyy_1, tg_yyyyyy_xxyyyy_0, \
                                         tg_yyyyyy_xxyyyy_1, tg_yyyyyy_xxyyyz_0, tg_yyyyyy_xxyyyz_1, tg_yyyyyy_xxyyz_1, \
                                         tg_yyyyyy_xxyyzz_0, tg_yyyyyy_xxyyzz_1, tg_yyyyyy_xxyzz_1, tg_yyyyyy_xxyzzz_0, \
                                         tg_yyyyyy_xxyzzz_1, tg_yyyyyy_xxzzz_1, tg_yyyyyy_xxzzzz_0, tg_yyyyyy_xxzzzz_1, \
                                         tg_yyyyyy_xyyyy_1, tg_yyyyyy_xyyyyy_0, tg_yyyyyy_xyyyyy_1, tg_yyyyyy_xyyyyz_0, \
                                         tg_yyyyyy_xyyyyz_1, tg_yyyyyy_xyyyz_1, tg_yyyyyy_xyyyzz_0, tg_yyyyyy_xyyyzz_1, \
                                         tg_yyyyyy_xyyzz_1, tg_yyyyyy_xyyzzz_0, tg_yyyyyy_xyyzzz_1, tg_yyyyyy_xyzzz_1, \
                                         tg_yyyyyy_xyzzzz_0, tg_yyyyyy_xyzzzz_1, tg_yyyyyy_xzzzz_1, tg_yyyyyy_xzzzzz_0, \
                                         tg_yyyyyy_xzzzzz_1, tg_yyyyyy_yyyyy_1, tg_yyyyyy_yyyyyy_0, tg_yyyyyy_yyyyyy_1, \
                                         tg_yyyyyy_yyyyyz_0, tg_yyyyyy_yyyyyz_1, tg_yyyyyy_yyyyz_1, tg_yyyyyy_yyyyzz_0, \
                                         tg_yyyyyy_yyyyzz_1, tg_yyyyyy_yyyzz_1, tg_yyyyyy_yyyzzz_0, tg_yyyyyy_yyyzzz_1, \
                                         tg_yyyyyy_yyzzz_1, tg_yyyyyy_yyzzzz_0, tg_yyyyyy_yyzzzz_1, tg_yyyyyy_yzzzz_1, \
                                         tg_yyyyyy_yzzzzz_0, tg_yyyyyy_yzzzzz_1, tg_yyyyyy_zzzzz_1, tg_yyyyyy_zzzzzz_0, \
                                         tg_yyyyyy_zzzzzz_1, tg_yyyyyz_xxxxx_1, tg_yyyyyz_xxxxxx_0, tg_yyyyyz_xxxxxx_1, \
                                         tg_yyyyyz_xxxxxy_0, tg_yyyyyz_xxxxxy_1, tg_yyyyyz_xxxxxz_0, tg_yyyyyz_xxxxxz_1, \
                                         tg_yyyyyz_xxxxy_1, tg_yyyyyz_xxxxyy_0, tg_yyyyyz_xxxxyy_1, tg_yyyyyz_xxxxyz_0, \
                                         tg_yyyyyz_xxxxyz_1, tg_yyyyyz_xxxxz_1, tg_yyyyyz_xxxxzz_0, tg_yyyyyz_xxxxzz_1, \
                                         tg_yyyyyz_xxxyy_1, tg_yyyyyz_xxxyyy_0, tg_yyyyyz_xxxyyy_1, tg_yyyyyz_xxxyyz_0, \
                                         tg_yyyyyz_xxxyyz_1, tg_yyyyyz_xxxyz_1, tg_yyyyyz_xxxyzz_0, tg_yyyyyz_xxxyzz_1, \
                                         tg_yyyyyz_xxxzz_1, tg_yyyyyz_xxxzzz_0, tg_yyyyyz_xxxzzz_1, tg_yyyyyz_xxyyy_1, \
                                         tg_yyyyyz_xxyyyy_0, tg_yyyyyz_xxyyyy_1, tg_yyyyyz_xxyyyz_0, tg_yyyyyz_xxyyyz_1, \
                                         tg_yyyyyz_xxyyz_1, tg_yyyyyz_xxyyzz_0, tg_yyyyyz_xxyyzz_1, tg_yyyyyz_xxyzz_1, \
                                         tg_yyyyyz_xxyzzz_0, tg_yyyyyz_xxyzzz_1, tg_yyyyyz_xxzzz_1, tg_yyyyyz_xxzzzz_0, \
                                         tg_yyyyyz_xxzzzz_1, tg_yyyyyz_xyyyy_1, tg_yyyyyz_xyyyyy_0, tg_yyyyyz_xyyyyy_1, \
                                         tg_yyyyyz_xyyyyz_0, tg_yyyyyz_xyyyyz_1, tg_yyyyyz_xyyyz_1, tg_yyyyyz_xyyyzz_0, \
                                         tg_yyyyyz_xyyyzz_1, tg_yyyyyz_xyyzz_1, tg_yyyyyz_xyyzzz_0, tg_yyyyyz_xyyzzz_1, \
                                         tg_yyyyyz_xyzzz_1, tg_yyyyyz_xyzzzz_0, tg_yyyyyz_xyzzzz_1, tg_yyyyyz_xzzzz_1, \
                                         tg_yyyyyz_xzzzzz_0, tg_yyyyyz_xzzzzz_1, tg_yyyyyz_yyyyy_1, tg_yyyyyz_yyyyyy_0, \
                                         tg_yyyyyz_yyyyyy_1, tg_yyyyyz_yyyyyz_0, tg_yyyyyz_yyyyyz_1, tg_yyyyyz_yyyyz_1, \
                                         tg_yyyyyz_yyyyzz_0, tg_yyyyyz_yyyyzz_1, tg_yyyyyz_yyyzz_1, tg_yyyyyz_yyyzzz_0, \
                                         tg_yyyyyz_yyyzzz_1, tg_yyyyyz_yyzzz_1, tg_yyyyyz_yyzzzz_0, tg_yyyyyz_yyzzzz_1, \
                                         tg_yyyyyz_yzzzz_1, tg_yyyyyz_yzzzzz_0, tg_yyyyyz_yzzzzz_1, tg_yyyyyz_zzzzz_1, \
                                         tg_yyyyyz_zzzzzz_0, tg_yyyyyz_zzzzzz_1, tg_yzzzz_xzzzzz_0, tg_yzzzz_xzzzzz_1, \
                                         tg_yzzzz_yyyyyy_0, tg_yzzzz_yyyyyy_1, tg_yzzzz_yyyyyz_0, tg_yzzzz_yyyyyz_1, \
                                         tg_yzzzz_yyyyzz_0, tg_yzzzz_yyyyzz_1, tg_yzzzz_yyyzzz_0, tg_yzzzz_yyyzzz_1, \
                                         tg_yzzzz_yyzzzz_0, tg_yzzzz_yyzzzz_1, tg_yzzzz_yzzzzz_0, tg_yzzzz_yzzzzz_1, \
                                         tg_yzzzz_zzzzzz_0, tg_yzzzz_zzzzzz_1, tg_zzzzz_xxxxxx_0, tg_zzzzz_xxxxxx_1, \
                                         tg_zzzzz_xxxxxy_0, tg_zzzzz_xxxxxy_1, tg_zzzzz_xxxxxz_0, tg_zzzzz_xxxxxz_1, \
                                         tg_zzzzz_xxxxyy_0, tg_zzzzz_xxxxyy_1, tg_zzzzz_xxxxyz_0, tg_zzzzz_xxxxyz_1, \
                                         tg_zzzzz_xxxxzz_0, tg_zzzzz_xxxxzz_1, tg_zzzzz_xxxyyy_0, tg_zzzzz_xxxyyy_1, \
                                         tg_zzzzz_xxxyyz_0, tg_zzzzz_xxxyyz_1, tg_zzzzz_xxxyzz_0, tg_zzzzz_xxxyzz_1, \
                                         tg_zzzzz_xxxzzz_0, tg_zzzzz_xxxzzz_1, tg_zzzzz_xxyyyy_0, tg_zzzzz_xxyyyy_1, \
                                         tg_zzzzz_xxyyyz_0, tg_zzzzz_xxyyyz_1, tg_zzzzz_xxyyzz_0, tg_zzzzz_xxyyzz_1, \
                                         tg_zzzzz_xxyzzz_0, tg_zzzzz_xxyzzz_1, tg_zzzzz_xxzzzz_0, tg_zzzzz_xxzzzz_1, \
                                         tg_zzzzz_xyyyyy_0, tg_zzzzz_xyyyyy_1, tg_zzzzz_xyyyyz_0, tg_zzzzz_xyyyyz_1, \
                                         tg_zzzzz_xyyyzz_0, tg_zzzzz_xyyyzz_1, tg_zzzzz_xyyzzz_0, tg_zzzzz_xyyzzz_1, \
                                         tg_zzzzz_xyzzzz_0, tg_zzzzz_xyzzzz_1, tg_zzzzz_xzzzzz_0, tg_zzzzz_xzzzzz_1, \
                                         tg_zzzzz_yyyyyy_0, tg_zzzzz_yyyyyy_1, tg_zzzzz_yyyyyz_0, tg_zzzzz_yyyyyz_1, \
                                         tg_zzzzz_yyyyzz_0, tg_zzzzz_yyyyzz_1, tg_zzzzz_yyyzzz_0, tg_zzzzz_yyyzzz_1, \
                                         tg_zzzzz_yyzzzz_0, tg_zzzzz_yyzzzz_1, tg_zzzzz_yzzzzz_0, tg_zzzzz_yzzzzz_1, \
                                         tg_zzzzz_zzzzzz_0, tg_zzzzz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyzzzz_xzzzzz_0[j] = pb_x * tg_xyzzzz_xzzzzz_0[j] + fr * tg_xyzzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xzzzzz_0[j] - tg_yzzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_zzzzz_1[j];

                    tg_xxyzzzz_yyyyyy_0[j] = pb_x * tg_xyzzzz_yyyyyy_0[j] + fr * tg_xyzzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyyyy_0[j] - tg_yzzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyzzzz_yyyyyz_0[j] = pb_x * tg_xyzzzz_yyyyyz_0[j] + fr * tg_xyzzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyyyz_0[j] - tg_yzzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyzzzz_yyyyzz_0[j] = pb_x * tg_xyzzzz_yyyyzz_0[j] + fr * tg_xyzzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyyzz_0[j] - tg_yzzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyzzzz_yyyzzz_0[j] = pb_x * tg_xyzzzz_yyyzzz_0[j] + fr * tg_xyzzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyzzz_0[j] - tg_yzzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_yyzzzz_0[j] = pb_x * tg_xyzzzz_yyzzzz_0[j] + fr * tg_xyzzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyzzzz_0[j] - tg_yzzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_yzzzzz_0[j] = pb_x * tg_xyzzzz_yzzzzz_0[j] + fr * tg_xyzzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yzzzzz_0[j] - tg_yzzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_zzzzzz_0[j] = pb_x * tg_xyzzzz_zzzzzz_0[j] + fr * tg_xyzzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_zzzzzz_0[j] - tg_yzzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_xxxxxx_0[j] = pb_x * tg_xzzzzz_xxxxxx_0[j] + fr * tg_xzzzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxxx_0[j] - tg_zzzzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzzzzz_xxxxx_1[j];

                    tg_xxzzzzz_xxxxxy_0[j] = pb_x * tg_xzzzzz_xxxxxy_0[j] + fr * tg_xzzzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxxy_0[j] - tg_zzzzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzzzz_xxxxy_1[j];

                    tg_xxzzzzz_xxxxxz_0[j] = pb_x * tg_xzzzzz_xxxxxz_0[j] + fr * tg_xzzzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxxz_0[j] - tg_zzzzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzzzz_xxxxz_1[j];

                    tg_xxzzzzz_xxxxyy_0[j] = pb_x * tg_xzzzzz_xxxxyy_0[j] + fr * tg_xzzzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxyy_0[j] - tg_zzzzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzz_xxxyy_1[j];

                    tg_xxzzzzz_xxxxyz_0[j] = pb_x * tg_xzzzzz_xxxxyz_0[j] + fr * tg_xzzzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxyz_0[j] - tg_zzzzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzz_xxxyz_1[j];

                    tg_xxzzzzz_xxxxzz_0[j] = pb_x * tg_xzzzzz_xxxxzz_0[j] + fr * tg_xzzzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxzz_0[j] - tg_zzzzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzz_xxxzz_1[j];

                    tg_xxzzzzz_xxxyyy_0[j] = pb_x * tg_xzzzzz_xxxyyy_0[j] + fr * tg_xzzzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyyy_0[j] - tg_zzzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxyyy_1[j];

                    tg_xxzzzzz_xxxyyz_0[j] = pb_x * tg_xzzzzz_xxxyyz_0[j] + fr * tg_xzzzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyyz_0[j] - tg_zzzzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxyyz_1[j];

                    tg_xxzzzzz_xxxyzz_0[j] = pb_x * tg_xzzzzz_xxxyzz_0[j] + fr * tg_xzzzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyzz_0[j] - tg_zzzzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxyzz_1[j];

                    tg_xxzzzzz_xxxzzz_0[j] = pb_x * tg_xzzzzz_xxxzzz_0[j] + fr * tg_xzzzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxzzz_0[j] - tg_zzzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxzzz_1[j];

                    tg_xxzzzzz_xxyyyy_0[j] = pb_x * tg_xzzzzz_xxyyyy_0[j] + fr * tg_xzzzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyyy_0[j] - tg_zzzzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyyyy_1[j];

                    tg_xxzzzzz_xxyyyz_0[j] = pb_x * tg_xzzzzz_xxyyyz_0[j] + fr * tg_xzzzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyyz_0[j] - tg_zzzzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyyyz_1[j];

                    tg_xxzzzzz_xxyyzz_0[j] = pb_x * tg_xzzzzz_xxyyzz_0[j] + fr * tg_xzzzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyzz_0[j] - tg_zzzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyyzz_1[j];

                    tg_xxzzzzz_xxyzzz_0[j] = pb_x * tg_xzzzzz_xxyzzz_0[j] + fr * tg_xzzzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyzzz_0[j] - tg_zzzzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyzzz_1[j];

                    tg_xxzzzzz_xxzzzz_0[j] = pb_x * tg_xzzzzz_xxzzzz_0[j] + fr * tg_xzzzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxzzzz_0[j] - tg_zzzzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xzzzz_1[j];

                    tg_xxzzzzz_xyyyyy_0[j] = pb_x * tg_xzzzzz_xyyyyy_0[j] + fr * tg_xzzzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyyy_0[j] - tg_zzzzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyyyy_1[j];

                    tg_xxzzzzz_xyyyyz_0[j] = pb_x * tg_xzzzzz_xyyyyz_0[j] + fr * tg_xzzzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyyz_0[j] - tg_zzzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyyyz_1[j];

                    tg_xxzzzzz_xyyyzz_0[j] = pb_x * tg_xzzzzz_xyyyzz_0[j] + fr * tg_xzzzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyzz_0[j] - tg_zzzzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyyzz_1[j];

                    tg_xxzzzzz_xyyzzz_0[j] = pb_x * tg_xzzzzz_xyyzzz_0[j] + fr * tg_xzzzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyzzz_0[j] - tg_zzzzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyzzz_1[j];

                    tg_xxzzzzz_xyzzzz_0[j] = pb_x * tg_xzzzzz_xyzzzz_0[j] + fr * tg_xzzzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyzzzz_0[j] - tg_zzzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yzzzz_1[j];

                    tg_xxzzzzz_xzzzzz_0[j] = pb_x * tg_xzzzzz_xzzzzz_0[j] + fr * tg_xzzzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzzzzz_0[j] - tg_zzzzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_zzzzz_1[j];

                    tg_xxzzzzz_yyyyyy_0[j] = pb_x * tg_xzzzzz_yyyyyy_0[j] + fr * tg_xzzzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyyy_0[j] - tg_zzzzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxzzzzz_yyyyyz_0[j] = pb_x * tg_xzzzzz_yyyyyz_0[j] + fr * tg_xzzzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyyz_0[j] - tg_zzzzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxzzzzz_yyyyzz_0[j] = pb_x * tg_xzzzzz_yyyyzz_0[j] + fr * tg_xzzzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyzz_0[j] - tg_zzzzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxzzzzz_yyyzzz_0[j] = pb_x * tg_xzzzzz_yyyzzz_0[j] + fr * tg_xzzzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyzzz_0[j] - tg_zzzzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_yyzzzz_0[j] = pb_x * tg_xzzzzz_yyzzzz_0[j] + fr * tg_xzzzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyzzzz_0[j] - tg_zzzzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_yzzzzz_0[j] = pb_x * tg_xzzzzz_yzzzzz_0[j] + fr * tg_xzzzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzzzzz_0[j] - tg_zzzzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_zzzzzz_0[j] = pb_x * tg_xzzzzz_zzzzzz_0[j] + fr * tg_xzzzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzzzzz_0[j] - tg_zzzzz_zzzzzz_1[j] * fl1_fza);

                    tg_xyyyyyy_xxxxxx_0[j] = pb_x * tg_yyyyyy_xxxxxx_0[j] + fr * tg_yyyyyy_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyyyy_xxxxx_1[j];

                    tg_xyyyyyy_xxxxxy_0[j] = pb_x * tg_yyyyyy_xxxxxy_0[j] + fr * tg_yyyyyy_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyyyy_xxxxy_1[j];

                    tg_xyyyyyy_xxxxxz_0[j] = pb_x * tg_yyyyyy_xxxxxz_0[j] + fr * tg_yyyyyy_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyyyy_xxxxz_1[j];

                    tg_xyyyyyy_xxxxyy_0[j] = pb_x * tg_yyyyyy_xxxxyy_0[j] + fr * tg_yyyyyy_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyyyy_xxxyy_1[j];

                    tg_xyyyyyy_xxxxyz_0[j] = pb_x * tg_yyyyyy_xxxxyz_0[j] + fr * tg_yyyyyy_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyyyy_xxxyz_1[j];

                    tg_xyyyyyy_xxxxzz_0[j] = pb_x * tg_yyyyyy_xxxxzz_0[j] + fr * tg_yyyyyy_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyyyy_xxxzz_1[j];

                    tg_xyyyyyy_xxxyyy_0[j] = pb_x * tg_yyyyyy_xxxyyy_0[j] + fr * tg_yyyyyy_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxyyy_1[j];

                    tg_xyyyyyy_xxxyyz_0[j] = pb_x * tg_yyyyyy_xxxyyz_0[j] + fr * tg_yyyyyy_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxyyz_1[j];

                    tg_xyyyyyy_xxxyzz_0[j] = pb_x * tg_yyyyyy_xxxyzz_0[j] + fr * tg_yyyyyy_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxyzz_1[j];

                    tg_xyyyyyy_xxxzzz_0[j] = pb_x * tg_yyyyyy_xxxzzz_0[j] + fr * tg_yyyyyy_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxzzz_1[j];

                    tg_xyyyyyy_xxyyyy_0[j] = pb_x * tg_yyyyyy_xxyyyy_0[j] + fr * tg_yyyyyy_xxyyyy_1[j] + fl1_fxn * tg_yyyyyy_xyyyy_1[j];

                    tg_xyyyyyy_xxyyyz_0[j] = pb_x * tg_yyyyyy_xxyyyz_0[j] + fr * tg_yyyyyy_xxyyyz_1[j] + fl1_fxn * tg_yyyyyy_xyyyz_1[j];

                    tg_xyyyyyy_xxyyzz_0[j] = pb_x * tg_yyyyyy_xxyyzz_0[j] + fr * tg_yyyyyy_xxyyzz_1[j] + fl1_fxn * tg_yyyyyy_xyyzz_1[j];

                    tg_xyyyyyy_xxyzzz_0[j] = pb_x * tg_yyyyyy_xxyzzz_0[j] + fr * tg_yyyyyy_xxyzzz_1[j] + fl1_fxn * tg_yyyyyy_xyzzz_1[j];

                    tg_xyyyyyy_xxzzzz_0[j] = pb_x * tg_yyyyyy_xxzzzz_0[j] + fr * tg_yyyyyy_xxzzzz_1[j] + fl1_fxn * tg_yyyyyy_xzzzz_1[j];

                    tg_xyyyyyy_xyyyyy_0[j] = pb_x * tg_yyyyyy_xyyyyy_0[j] + fr * tg_yyyyyy_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyyyy_1[j];

                    tg_xyyyyyy_xyyyyz_0[j] = pb_x * tg_yyyyyy_xyyyyz_0[j] + fr * tg_yyyyyy_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyyyz_1[j];

                    tg_xyyyyyy_xyyyzz_0[j] = pb_x * tg_yyyyyy_xyyyzz_0[j] + fr * tg_yyyyyy_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyyzz_1[j];

                    tg_xyyyyyy_xyyzzz_0[j] = pb_x * tg_yyyyyy_xyyzzz_0[j] + fr * tg_yyyyyy_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyzzz_1[j];

                    tg_xyyyyyy_xyzzzz_0[j] = pb_x * tg_yyyyyy_xyzzzz_0[j] + fr * tg_yyyyyy_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yzzzz_1[j];

                    tg_xyyyyyy_xzzzzz_0[j] = pb_x * tg_yyyyyy_xzzzzz_0[j] + fr * tg_yyyyyy_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_zzzzz_1[j];

                    tg_xyyyyyy_yyyyyy_0[j] = pb_x * tg_yyyyyy_yyyyyy_0[j] + fr * tg_yyyyyy_yyyyyy_1[j];

                    tg_xyyyyyy_yyyyyz_0[j] = pb_x * tg_yyyyyy_yyyyyz_0[j] + fr * tg_yyyyyy_yyyyyz_1[j];

                    tg_xyyyyyy_yyyyzz_0[j] = pb_x * tg_yyyyyy_yyyyzz_0[j] + fr * tg_yyyyyy_yyyyzz_1[j];

                    tg_xyyyyyy_yyyzzz_0[j] = pb_x * tg_yyyyyy_yyyzzz_0[j] + fr * tg_yyyyyy_yyyzzz_1[j];

                    tg_xyyyyyy_yyzzzz_0[j] = pb_x * tg_yyyyyy_yyzzzz_0[j] + fr * tg_yyyyyy_yyzzzz_1[j];

                    tg_xyyyyyy_yzzzzz_0[j] = pb_x * tg_yyyyyy_yzzzzz_0[j] + fr * tg_yyyyyy_yzzzzz_1[j];

                    tg_xyyyyyy_zzzzzz_0[j] = pb_x * tg_yyyyyy_zzzzzz_0[j] + fr * tg_yyyyyy_zzzzzz_1[j];

                    tg_xyyyyyz_xxxxxx_0[j] = pb_x * tg_yyyyyz_xxxxxx_0[j] + fr * tg_yyyyyz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyyyz_xxxxx_1[j];

                    tg_xyyyyyz_xxxxxy_0[j] = pb_x * tg_yyyyyz_xxxxxy_0[j] + fr * tg_yyyyyz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyyyz_xxxxy_1[j];

                    tg_xyyyyyz_xxxxxz_0[j] = pb_x * tg_yyyyyz_xxxxxz_0[j] + fr * tg_yyyyyz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyyyz_xxxxz_1[j];

                    tg_xyyyyyz_xxxxyy_0[j] = pb_x * tg_yyyyyz_xxxxyy_0[j] + fr * tg_yyyyyz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyyyz_xxxyy_1[j];

                    tg_xyyyyyz_xxxxyz_0[j] = pb_x * tg_yyyyyz_xxxxyz_0[j] + fr * tg_yyyyyz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyyyz_xxxyz_1[j];

                    tg_xyyyyyz_xxxxzz_0[j] = pb_x * tg_yyyyyz_xxxxzz_0[j] + fr * tg_yyyyyz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyyyz_xxxzz_1[j];

                    tg_xyyyyyz_xxxyyy_0[j] = pb_x * tg_yyyyyz_xxxyyy_0[j] + fr * tg_yyyyyz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxyyy_1[j];

                    tg_xyyyyyz_xxxyyz_0[j] = pb_x * tg_yyyyyz_xxxyyz_0[j] + fr * tg_yyyyyz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxyyz_1[j];

                    tg_xyyyyyz_xxxyzz_0[j] = pb_x * tg_yyyyyz_xxxyzz_0[j] + fr * tg_yyyyyz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxyzz_1[j];

                    tg_xyyyyyz_xxxzzz_0[j] = pb_x * tg_yyyyyz_xxxzzz_0[j] + fr * tg_yyyyyz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxzzz_1[j];

                    tg_xyyyyyz_xxyyyy_0[j] = pb_x * tg_yyyyyz_xxyyyy_0[j] + fr * tg_yyyyyz_xxyyyy_1[j] + fl1_fxn * tg_yyyyyz_xyyyy_1[j];

                    tg_xyyyyyz_xxyyyz_0[j] = pb_x * tg_yyyyyz_xxyyyz_0[j] + fr * tg_yyyyyz_xxyyyz_1[j] + fl1_fxn * tg_yyyyyz_xyyyz_1[j];

                    tg_xyyyyyz_xxyyzz_0[j] = pb_x * tg_yyyyyz_xxyyzz_0[j] + fr * tg_yyyyyz_xxyyzz_1[j] + fl1_fxn * tg_yyyyyz_xyyzz_1[j];

                    tg_xyyyyyz_xxyzzz_0[j] = pb_x * tg_yyyyyz_xxyzzz_0[j] + fr * tg_yyyyyz_xxyzzz_1[j] + fl1_fxn * tg_yyyyyz_xyzzz_1[j];

                    tg_xyyyyyz_xxzzzz_0[j] = pb_x * tg_yyyyyz_xxzzzz_0[j] + fr * tg_yyyyyz_xxzzzz_1[j] + fl1_fxn * tg_yyyyyz_xzzzz_1[j];

                    tg_xyyyyyz_xyyyyy_0[j] = pb_x * tg_yyyyyz_xyyyyy_0[j] + fr * tg_yyyyyz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyyyy_1[j];

                    tg_xyyyyyz_xyyyyz_0[j] = pb_x * tg_yyyyyz_xyyyyz_0[j] + fr * tg_yyyyyz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyyyz_1[j];

                    tg_xyyyyyz_xyyyzz_0[j] = pb_x * tg_yyyyyz_xyyyzz_0[j] + fr * tg_yyyyyz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyyzz_1[j];

                    tg_xyyyyyz_xyyzzz_0[j] = pb_x * tg_yyyyyz_xyyzzz_0[j] + fr * tg_yyyyyz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyzzz_1[j];

                    tg_xyyyyyz_xyzzzz_0[j] = pb_x * tg_yyyyyz_xyzzzz_0[j] + fr * tg_yyyyyz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yzzzz_1[j];

                    tg_xyyyyyz_xzzzzz_0[j] = pb_x * tg_yyyyyz_xzzzzz_0[j] + fr * tg_yyyyyz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_zzzzz_1[j];

                    tg_xyyyyyz_yyyyyy_0[j] = pb_x * tg_yyyyyz_yyyyyy_0[j] + fr * tg_yyyyyz_yyyyyy_1[j];

                    tg_xyyyyyz_yyyyyz_0[j] = pb_x * tg_yyyyyz_yyyyyz_0[j] + fr * tg_yyyyyz_yyyyyz_1[j];

                    tg_xyyyyyz_yyyyzz_0[j] = pb_x * tg_yyyyyz_yyyyzz_0[j] + fr * tg_yyyyyz_yyyyzz_1[j];

                    tg_xyyyyyz_yyyzzz_0[j] = pb_x * tg_yyyyyz_yyyzzz_0[j] + fr * tg_yyyyyz_yyyzzz_1[j];

                    tg_xyyyyyz_yyzzzz_0[j] = pb_x * tg_yyyyyz_yyzzzz_0[j] + fr * tg_yyyyyz_yyzzzz_1[j];

                    tg_xyyyyyz_yzzzzz_0[j] = pb_x * tg_yyyyyz_yzzzzz_0[j] + fr * tg_yyyyyz_yzzzzz_1[j];

                    tg_xyyyyyz_zzzzzz_0[j] = pb_x * tg_yyyyyz_zzzzzz_0[j] + fr * tg_yyyyyz_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_644_735(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (644,735)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {7, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyyyzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 644); 

                auto tg_yyyyzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 645); 

                auto tg_yyyyzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 646); 

                auto tg_yyyyzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 647); 

                auto tg_yyyyzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 648); 

                auto tg_yyyyzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 649); 

                auto tg_yyyyzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 650); 

                auto tg_yyyyzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 651); 

                auto tg_yyyyzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 652); 

                auto tg_yyyyzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 653); 

                auto tg_yyyyzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 654); 

                auto tg_yyyyzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 655); 

                auto tg_yyyyzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 656); 

                auto tg_yyyyzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 657); 

                auto tg_yyyyzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 658); 

                auto tg_yyyyzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 659); 

                auto tg_yyyyzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 660); 

                auto tg_yyyyzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 661); 

                auto tg_yyyyzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 662); 

                auto tg_yyyyzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 663); 

                auto tg_yyyyzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 664); 

                auto tg_yyyyzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 665); 

                auto tg_yyyyzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 666); 

                auto tg_yyyyzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 667); 

                auto tg_yyyyzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 668); 

                auto tg_yyyyzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 669); 

                auto tg_yyyyzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 670); 

                auto tg_yyyyzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 671); 

                auto tg_yyyzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 672); 

                auto tg_yyyzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 673); 

                auto tg_yyyzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 674); 

                auto tg_yyyzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 675); 

                auto tg_yyyzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 676); 

                auto tg_yyyzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 677); 

                auto tg_yyyzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 678); 

                auto tg_yyyzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 679); 

                auto tg_yyyzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 680); 

                auto tg_yyyzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 681); 

                auto tg_yyyzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 682); 

                auto tg_yyyzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 683); 

                auto tg_yyyzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 684); 

                auto tg_yyyzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 685); 

                auto tg_yyyzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 686); 

                auto tg_yyyzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 687); 

                auto tg_yyyzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 688); 

                auto tg_yyyzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 689); 

                auto tg_yyyzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 690); 

                auto tg_yyyzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 691); 

                auto tg_yyyzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 692); 

                auto tg_yyyzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 693); 

                auto tg_yyyzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 694); 

                auto tg_yyyzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 695); 

                auto tg_yyyzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 696); 

                auto tg_yyyzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 697); 

                auto tg_yyyzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 698); 

                auto tg_yyyzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 699); 

                auto tg_yyzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 700); 

                auto tg_yyzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 701); 

                auto tg_yyzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 702); 

                auto tg_yyzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 703); 

                auto tg_yyzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 704); 

                auto tg_yyzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 705); 

                auto tg_yyzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 706); 

                auto tg_yyzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 707); 

                auto tg_yyzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 708); 

                auto tg_yyzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 709); 

                auto tg_yyzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 710); 

                auto tg_yyzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 711); 

                auto tg_yyzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 712); 

                auto tg_yyzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 713); 

                auto tg_yyzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 714); 

                auto tg_yyzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 715); 

                auto tg_yyzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 716); 

                auto tg_yyzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 717); 

                auto tg_yyzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 718); 

                auto tg_yyzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 719); 

                auto tg_yyzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 720); 

                auto tg_yyzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 721); 

                auto tg_yyzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 722); 

                auto tg_yyzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 723); 

                auto tg_yyzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 724); 

                auto tg_yyzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 725); 

                auto tg_yyzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 726); 

                auto tg_yyzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 727); 

                auto tg_yzzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 728); 

                auto tg_yzzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 729); 

                auto tg_yzzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 730); 

                auto tg_yzzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 731); 

                auto tg_yzzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 732); 

                auto tg_yzzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 733); 

                auto tg_yzzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 734); 

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

                // set up pointers to integrals

                auto tg_xyyyyzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 644); 

                auto tg_xyyyyzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 645); 

                auto tg_xyyyyzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 646); 

                auto tg_xyyyyzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 647); 

                auto tg_xyyyyzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 648); 

                auto tg_xyyyyzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 649); 

                auto tg_xyyyyzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 650); 

                auto tg_xyyyyzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 651); 

                auto tg_xyyyyzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 652); 

                auto tg_xyyyyzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 653); 

                auto tg_xyyyyzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 654); 

                auto tg_xyyyyzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 655); 

                auto tg_xyyyyzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 656); 

                auto tg_xyyyyzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 657); 

                auto tg_xyyyyzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 658); 

                auto tg_xyyyyzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 659); 

                auto tg_xyyyyzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 660); 

                auto tg_xyyyyzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 661); 

                auto tg_xyyyyzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 662); 

                auto tg_xyyyyzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 663); 

                auto tg_xyyyyzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 664); 

                auto tg_xyyyyzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 665); 

                auto tg_xyyyyzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 666); 

                auto tg_xyyyyzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 667); 

                auto tg_xyyyyzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 668); 

                auto tg_xyyyyzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 669); 

                auto tg_xyyyyzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 670); 

                auto tg_xyyyyzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 671); 

                auto tg_xyyyzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 672); 

                auto tg_xyyyzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 673); 

                auto tg_xyyyzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 674); 

                auto tg_xyyyzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 675); 

                auto tg_xyyyzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 676); 

                auto tg_xyyyzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 677); 

                auto tg_xyyyzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 678); 

                auto tg_xyyyzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 679); 

                auto tg_xyyyzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 680); 

                auto tg_xyyyzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 681); 

                auto tg_xyyyzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 682); 

                auto tg_xyyyzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 683); 

                auto tg_xyyyzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 684); 

                auto tg_xyyyzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 685); 

                auto tg_xyyyzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 686); 

                auto tg_xyyyzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 687); 

                auto tg_xyyyzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 688); 

                auto tg_xyyyzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 689); 

                auto tg_xyyyzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 690); 

                auto tg_xyyyzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 691); 

                auto tg_xyyyzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 692); 

                auto tg_xyyyzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 693); 

                auto tg_xyyyzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 694); 

                auto tg_xyyyzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 695); 

                auto tg_xyyyzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 696); 

                auto tg_xyyyzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 697); 

                auto tg_xyyyzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 698); 

                auto tg_xyyyzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 699); 

                auto tg_xyyzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 700); 

                auto tg_xyyzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 701); 

                auto tg_xyyzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 702); 

                auto tg_xyyzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 703); 

                auto tg_xyyzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 704); 

                auto tg_xyyzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 705); 

                auto tg_xyyzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 706); 

                auto tg_xyyzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 707); 

                auto tg_xyyzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 708); 

                auto tg_xyyzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 709); 

                auto tg_xyyzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 710); 

                auto tg_xyyzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 711); 

                auto tg_xyyzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 712); 

                auto tg_xyyzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 713); 

                auto tg_xyyzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 714); 

                auto tg_xyyzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 715); 

                auto tg_xyyzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 716); 

                auto tg_xyyzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 717); 

                auto tg_xyyzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 718); 

                auto tg_xyyzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 719); 

                auto tg_xyyzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 720); 

                auto tg_xyyzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 721); 

                auto tg_xyyzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 722); 

                auto tg_xyyzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 723); 

                auto tg_xyyzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 724); 

                auto tg_xyyzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 725); 

                auto tg_xyyzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 726); 

                auto tg_xyyzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 727); 

                auto tg_xyzzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 728); 

                auto tg_xyzzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 729); 

                auto tg_xyzzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 730); 

                auto tg_xyzzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 731); 

                auto tg_xyzzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 732); 

                auto tg_xyzzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 733); 

                auto tg_xyzzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 734); 

                // Batch of Integrals (644,735)

                #pragma omp simd aligned(fxn, tg_xyyyyzz_xxxxxx_0, tg_xyyyyzz_xxxxxy_0, tg_xyyyyzz_xxxxxz_0, \
                                         tg_xyyyyzz_xxxxyy_0, tg_xyyyyzz_xxxxyz_0, tg_xyyyyzz_xxxxzz_0, tg_xyyyyzz_xxxyyy_0, \
                                         tg_xyyyyzz_xxxyyz_0, tg_xyyyyzz_xxxyzz_0, tg_xyyyyzz_xxxzzz_0, tg_xyyyyzz_xxyyyy_0, \
                                         tg_xyyyyzz_xxyyyz_0, tg_xyyyyzz_xxyyzz_0, tg_xyyyyzz_xxyzzz_0, tg_xyyyyzz_xxzzzz_0, \
                                         tg_xyyyyzz_xyyyyy_0, tg_xyyyyzz_xyyyyz_0, tg_xyyyyzz_xyyyzz_0, tg_xyyyyzz_xyyzzz_0, \
                                         tg_xyyyyzz_xyzzzz_0, tg_xyyyyzz_xzzzzz_0, tg_xyyyyzz_yyyyyy_0, tg_xyyyyzz_yyyyyz_0, \
                                         tg_xyyyyzz_yyyyzz_0, tg_xyyyyzz_yyyzzz_0, tg_xyyyyzz_yyzzzz_0, tg_xyyyyzz_yzzzzz_0, \
                                         tg_xyyyyzz_zzzzzz_0, tg_xyyyzzz_xxxxxx_0, tg_xyyyzzz_xxxxxy_0, tg_xyyyzzz_xxxxxz_0, \
                                         tg_xyyyzzz_xxxxyy_0, tg_xyyyzzz_xxxxyz_0, tg_xyyyzzz_xxxxzz_0, tg_xyyyzzz_xxxyyy_0, \
                                         tg_xyyyzzz_xxxyyz_0, tg_xyyyzzz_xxxyzz_0, tg_xyyyzzz_xxxzzz_0, tg_xyyyzzz_xxyyyy_0, \
                                         tg_xyyyzzz_xxyyyz_0, tg_xyyyzzz_xxyyzz_0, tg_xyyyzzz_xxyzzz_0, tg_xyyyzzz_xxzzzz_0, \
                                         tg_xyyyzzz_xyyyyy_0, tg_xyyyzzz_xyyyyz_0, tg_xyyyzzz_xyyyzz_0, tg_xyyyzzz_xyyzzz_0, \
                                         tg_xyyyzzz_xyzzzz_0, tg_xyyyzzz_xzzzzz_0, tg_xyyyzzz_yyyyyy_0, tg_xyyyzzz_yyyyyz_0, \
                                         tg_xyyyzzz_yyyyzz_0, tg_xyyyzzz_yyyzzz_0, tg_xyyyzzz_yyzzzz_0, tg_xyyyzzz_yzzzzz_0, \
                                         tg_xyyyzzz_zzzzzz_0, tg_xyyzzzz_xxxxxx_0, tg_xyyzzzz_xxxxxy_0, tg_xyyzzzz_xxxxxz_0, \
                                         tg_xyyzzzz_xxxxyy_0, tg_xyyzzzz_xxxxyz_0, tg_xyyzzzz_xxxxzz_0, tg_xyyzzzz_xxxyyy_0, \
                                         tg_xyyzzzz_xxxyyz_0, tg_xyyzzzz_xxxyzz_0, tg_xyyzzzz_xxxzzz_0, tg_xyyzzzz_xxyyyy_0, \
                                         tg_xyyzzzz_xxyyyz_0, tg_xyyzzzz_xxyyzz_0, tg_xyyzzzz_xxyzzz_0, tg_xyyzzzz_xxzzzz_0, \
                                         tg_xyyzzzz_xyyyyy_0, tg_xyyzzzz_xyyyyz_0, tg_xyyzzzz_xyyyzz_0, tg_xyyzzzz_xyyzzz_0, \
                                         tg_xyyzzzz_xyzzzz_0, tg_xyyzzzz_xzzzzz_0, tg_xyyzzzz_yyyyyy_0, tg_xyyzzzz_yyyyyz_0, \
                                         tg_xyyzzzz_yyyyzz_0, tg_xyyzzzz_yyyzzz_0, tg_xyyzzzz_yyzzzz_0, tg_xyyzzzz_yzzzzz_0, \
                                         tg_xyyzzzz_zzzzzz_0, tg_xyzzzzz_xxxxxx_0, tg_xyzzzzz_xxxxxy_0, tg_xyzzzzz_xxxxxz_0, \
                                         tg_xyzzzzz_xxxxyy_0, tg_xyzzzzz_xxxxyz_0, tg_xyzzzzz_xxxxzz_0, tg_xyzzzzz_xxxyyy_0, \
                                         tg_yyyyzz_xxxxx_1, tg_yyyyzz_xxxxxx_0, tg_yyyyzz_xxxxxx_1, tg_yyyyzz_xxxxxy_0, \
                                         tg_yyyyzz_xxxxxy_1, tg_yyyyzz_xxxxxz_0, tg_yyyyzz_xxxxxz_1, tg_yyyyzz_xxxxy_1, \
                                         tg_yyyyzz_xxxxyy_0, tg_yyyyzz_xxxxyy_1, tg_yyyyzz_xxxxyz_0, tg_yyyyzz_xxxxyz_1, \
                                         tg_yyyyzz_xxxxz_1, tg_yyyyzz_xxxxzz_0, tg_yyyyzz_xxxxzz_1, tg_yyyyzz_xxxyy_1, \
                                         tg_yyyyzz_xxxyyy_0, tg_yyyyzz_xxxyyy_1, tg_yyyyzz_xxxyyz_0, tg_yyyyzz_xxxyyz_1, \
                                         tg_yyyyzz_xxxyz_1, tg_yyyyzz_xxxyzz_0, tg_yyyyzz_xxxyzz_1, tg_yyyyzz_xxxzz_1, \
                                         tg_yyyyzz_xxxzzz_0, tg_yyyyzz_xxxzzz_1, tg_yyyyzz_xxyyy_1, tg_yyyyzz_xxyyyy_0, \
                                         tg_yyyyzz_xxyyyy_1, tg_yyyyzz_xxyyyz_0, tg_yyyyzz_xxyyyz_1, tg_yyyyzz_xxyyz_1, \
                                         tg_yyyyzz_xxyyzz_0, tg_yyyyzz_xxyyzz_1, tg_yyyyzz_xxyzz_1, tg_yyyyzz_xxyzzz_0, \
                                         tg_yyyyzz_xxyzzz_1, tg_yyyyzz_xxzzz_1, tg_yyyyzz_xxzzzz_0, tg_yyyyzz_xxzzzz_1, \
                                         tg_yyyyzz_xyyyy_1, tg_yyyyzz_xyyyyy_0, tg_yyyyzz_xyyyyy_1, tg_yyyyzz_xyyyyz_0, \
                                         tg_yyyyzz_xyyyyz_1, tg_yyyyzz_xyyyz_1, tg_yyyyzz_xyyyzz_0, tg_yyyyzz_xyyyzz_1, \
                                         tg_yyyyzz_xyyzz_1, tg_yyyyzz_xyyzzz_0, tg_yyyyzz_xyyzzz_1, tg_yyyyzz_xyzzz_1, \
                                         tg_yyyyzz_xyzzzz_0, tg_yyyyzz_xyzzzz_1, tg_yyyyzz_xzzzz_1, tg_yyyyzz_xzzzzz_0, \
                                         tg_yyyyzz_xzzzzz_1, tg_yyyyzz_yyyyy_1, tg_yyyyzz_yyyyyy_0, tg_yyyyzz_yyyyyy_1, \
                                         tg_yyyyzz_yyyyyz_0, tg_yyyyzz_yyyyyz_1, tg_yyyyzz_yyyyz_1, tg_yyyyzz_yyyyzz_0, \
                                         tg_yyyyzz_yyyyzz_1, tg_yyyyzz_yyyzz_1, tg_yyyyzz_yyyzzz_0, tg_yyyyzz_yyyzzz_1, \
                                         tg_yyyyzz_yyzzz_1, tg_yyyyzz_yyzzzz_0, tg_yyyyzz_yyzzzz_1, tg_yyyyzz_yzzzz_1, \
                                         tg_yyyyzz_yzzzzz_0, tg_yyyyzz_yzzzzz_1, tg_yyyyzz_zzzzz_1, tg_yyyyzz_zzzzzz_0, \
                                         tg_yyyyzz_zzzzzz_1, tg_yyyzzz_xxxxx_1, tg_yyyzzz_xxxxxx_0, tg_yyyzzz_xxxxxx_1, \
                                         tg_yyyzzz_xxxxxy_0, tg_yyyzzz_xxxxxy_1, tg_yyyzzz_xxxxxz_0, tg_yyyzzz_xxxxxz_1, \
                                         tg_yyyzzz_xxxxy_1, tg_yyyzzz_xxxxyy_0, tg_yyyzzz_xxxxyy_1, tg_yyyzzz_xxxxyz_0, \
                                         tg_yyyzzz_xxxxyz_1, tg_yyyzzz_xxxxz_1, tg_yyyzzz_xxxxzz_0, tg_yyyzzz_xxxxzz_1, \
                                         tg_yyyzzz_xxxyy_1, tg_yyyzzz_xxxyyy_0, tg_yyyzzz_xxxyyy_1, tg_yyyzzz_xxxyyz_0, \
                                         tg_yyyzzz_xxxyyz_1, tg_yyyzzz_xxxyz_1, tg_yyyzzz_xxxyzz_0, tg_yyyzzz_xxxyzz_1, \
                                         tg_yyyzzz_xxxzz_1, tg_yyyzzz_xxxzzz_0, tg_yyyzzz_xxxzzz_1, tg_yyyzzz_xxyyy_1, \
                                         tg_yyyzzz_xxyyyy_0, tg_yyyzzz_xxyyyy_1, tg_yyyzzz_xxyyyz_0, tg_yyyzzz_xxyyyz_1, \
                                         tg_yyyzzz_xxyyz_1, tg_yyyzzz_xxyyzz_0, tg_yyyzzz_xxyyzz_1, tg_yyyzzz_xxyzz_1, \
                                         tg_yyyzzz_xxyzzz_0, tg_yyyzzz_xxyzzz_1, tg_yyyzzz_xxzzz_1, tg_yyyzzz_xxzzzz_0, \
                                         tg_yyyzzz_xxzzzz_1, tg_yyyzzz_xyyyy_1, tg_yyyzzz_xyyyyy_0, tg_yyyzzz_xyyyyy_1, \
                                         tg_yyyzzz_xyyyyz_0, tg_yyyzzz_xyyyyz_1, tg_yyyzzz_xyyyz_1, tg_yyyzzz_xyyyzz_0, \
                                         tg_yyyzzz_xyyyzz_1, tg_yyyzzz_xyyzz_1, tg_yyyzzz_xyyzzz_0, tg_yyyzzz_xyyzzz_1, \
                                         tg_yyyzzz_xyzzz_1, tg_yyyzzz_xyzzzz_0, tg_yyyzzz_xyzzzz_1, tg_yyyzzz_xzzzz_1, \
                                         tg_yyyzzz_xzzzzz_0, tg_yyyzzz_xzzzzz_1, tg_yyyzzz_yyyyy_1, tg_yyyzzz_yyyyyy_0, \
                                         tg_yyyzzz_yyyyyy_1, tg_yyyzzz_yyyyyz_0, tg_yyyzzz_yyyyyz_1, tg_yyyzzz_yyyyz_1, \
                                         tg_yyyzzz_yyyyzz_0, tg_yyyzzz_yyyyzz_1, tg_yyyzzz_yyyzz_1, tg_yyyzzz_yyyzzz_0, \
                                         tg_yyyzzz_yyyzzz_1, tg_yyyzzz_yyzzz_1, tg_yyyzzz_yyzzzz_0, tg_yyyzzz_yyzzzz_1, \
                                         tg_yyyzzz_yzzzz_1, tg_yyyzzz_yzzzzz_0, tg_yyyzzz_yzzzzz_1, tg_yyyzzz_zzzzz_1, \
                                         tg_yyyzzz_zzzzzz_0, tg_yyyzzz_zzzzzz_1, tg_yyzzzz_xxxxx_1, tg_yyzzzz_xxxxxx_0, \
                                         tg_yyzzzz_xxxxxx_1, tg_yyzzzz_xxxxxy_0, tg_yyzzzz_xxxxxy_1, tg_yyzzzz_xxxxxz_0, \
                                         tg_yyzzzz_xxxxxz_1, tg_yyzzzz_xxxxy_1, tg_yyzzzz_xxxxyy_0, tg_yyzzzz_xxxxyy_1, \
                                         tg_yyzzzz_xxxxyz_0, tg_yyzzzz_xxxxyz_1, tg_yyzzzz_xxxxz_1, tg_yyzzzz_xxxxzz_0, \
                                         tg_yyzzzz_xxxxzz_1, tg_yyzzzz_xxxyy_1, tg_yyzzzz_xxxyyy_0, tg_yyzzzz_xxxyyy_1, \
                                         tg_yyzzzz_xxxyyz_0, tg_yyzzzz_xxxyyz_1, tg_yyzzzz_xxxyz_1, tg_yyzzzz_xxxyzz_0, \
                                         tg_yyzzzz_xxxyzz_1, tg_yyzzzz_xxxzz_1, tg_yyzzzz_xxxzzz_0, tg_yyzzzz_xxxzzz_1, \
                                         tg_yyzzzz_xxyyy_1, tg_yyzzzz_xxyyyy_0, tg_yyzzzz_xxyyyy_1, tg_yyzzzz_xxyyyz_0, \
                                         tg_yyzzzz_xxyyyz_1, tg_yyzzzz_xxyyz_1, tg_yyzzzz_xxyyzz_0, tg_yyzzzz_xxyyzz_1, \
                                         tg_yyzzzz_xxyzz_1, tg_yyzzzz_xxyzzz_0, tg_yyzzzz_xxyzzz_1, tg_yyzzzz_xxzzz_1, \
                                         tg_yyzzzz_xxzzzz_0, tg_yyzzzz_xxzzzz_1, tg_yyzzzz_xyyyy_1, tg_yyzzzz_xyyyyy_0, \
                                         tg_yyzzzz_xyyyyy_1, tg_yyzzzz_xyyyyz_0, tg_yyzzzz_xyyyyz_1, tg_yyzzzz_xyyyz_1, \
                                         tg_yyzzzz_xyyyzz_0, tg_yyzzzz_xyyyzz_1, tg_yyzzzz_xyyzz_1, tg_yyzzzz_xyyzzz_0, \
                                         tg_yyzzzz_xyyzzz_1, tg_yyzzzz_xyzzz_1, tg_yyzzzz_xyzzzz_0, tg_yyzzzz_xyzzzz_1, \
                                         tg_yyzzzz_xzzzz_1, tg_yyzzzz_xzzzzz_0, tg_yyzzzz_xzzzzz_1, tg_yyzzzz_yyyyy_1, \
                                         tg_yyzzzz_yyyyyy_0, tg_yyzzzz_yyyyyy_1, tg_yyzzzz_yyyyyz_0, tg_yyzzzz_yyyyyz_1, \
                                         tg_yyzzzz_yyyyz_1, tg_yyzzzz_yyyyzz_0, tg_yyzzzz_yyyyzz_1, tg_yyzzzz_yyyzz_1, \
                                         tg_yyzzzz_yyyzzz_0, tg_yyzzzz_yyyzzz_1, tg_yyzzzz_yyzzz_1, tg_yyzzzz_yyzzzz_0, \
                                         tg_yyzzzz_yyzzzz_1, tg_yyzzzz_yzzzz_1, tg_yyzzzz_yzzzzz_0, tg_yyzzzz_yzzzzz_1, \
                                         tg_yyzzzz_zzzzz_1, tg_yyzzzz_zzzzzz_0, tg_yyzzzz_zzzzzz_1, tg_yzzzzz_xxxxx_1, \
                                         tg_yzzzzz_xxxxxx_0, tg_yzzzzz_xxxxxx_1, tg_yzzzzz_xxxxxy_0, tg_yzzzzz_xxxxxy_1, \
                                         tg_yzzzzz_xxxxxz_0, tg_yzzzzz_xxxxxz_1, tg_yzzzzz_xxxxy_1, tg_yzzzzz_xxxxyy_0, \
                                         tg_yzzzzz_xxxxyy_1, tg_yzzzzz_xxxxyz_0, tg_yzzzzz_xxxxyz_1, tg_yzzzzz_xxxxz_1, \
                                         tg_yzzzzz_xxxxzz_0, tg_yzzzzz_xxxxzz_1, tg_yzzzzz_xxxyy_1, tg_yzzzzz_xxxyyy_0, \
                                         tg_yzzzzz_xxxyyy_1, tg_yzzzzz_xxxyz_1, tg_yzzzzz_xxxzz_1, tg_yzzzzz_xxyyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyyyzz_xxxxxx_0[j] = pb_x * tg_yyyyzz_xxxxxx_0[j] + fr * tg_yyyyzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyyzz_xxxxx_1[j];

                    tg_xyyyyzz_xxxxxy_0[j] = pb_x * tg_yyyyzz_xxxxxy_0[j] + fr * tg_yyyyzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyyzz_xxxxy_1[j];

                    tg_xyyyyzz_xxxxxz_0[j] = pb_x * tg_yyyyzz_xxxxxz_0[j] + fr * tg_yyyyzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyyzz_xxxxz_1[j];

                    tg_xyyyyzz_xxxxyy_0[j] = pb_x * tg_yyyyzz_xxxxyy_0[j] + fr * tg_yyyyzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyyzz_xxxyy_1[j];

                    tg_xyyyyzz_xxxxyz_0[j] = pb_x * tg_yyyyzz_xxxxyz_0[j] + fr * tg_yyyyzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyyzz_xxxyz_1[j];

                    tg_xyyyyzz_xxxxzz_0[j] = pb_x * tg_yyyyzz_xxxxzz_0[j] + fr * tg_yyyyzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyyzz_xxxzz_1[j];

                    tg_xyyyyzz_xxxyyy_0[j] = pb_x * tg_yyyyzz_xxxyyy_0[j] + fr * tg_yyyyzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxyyy_1[j];

                    tg_xyyyyzz_xxxyyz_0[j] = pb_x * tg_yyyyzz_xxxyyz_0[j] + fr * tg_yyyyzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxyyz_1[j];

                    tg_xyyyyzz_xxxyzz_0[j] = pb_x * tg_yyyyzz_xxxyzz_0[j] + fr * tg_yyyyzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxyzz_1[j];

                    tg_xyyyyzz_xxxzzz_0[j] = pb_x * tg_yyyyzz_xxxzzz_0[j] + fr * tg_yyyyzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxzzz_1[j];

                    tg_xyyyyzz_xxyyyy_0[j] = pb_x * tg_yyyyzz_xxyyyy_0[j] + fr * tg_yyyyzz_xxyyyy_1[j] + fl1_fxn * tg_yyyyzz_xyyyy_1[j];

                    tg_xyyyyzz_xxyyyz_0[j] = pb_x * tg_yyyyzz_xxyyyz_0[j] + fr * tg_yyyyzz_xxyyyz_1[j] + fl1_fxn * tg_yyyyzz_xyyyz_1[j];

                    tg_xyyyyzz_xxyyzz_0[j] = pb_x * tg_yyyyzz_xxyyzz_0[j] + fr * tg_yyyyzz_xxyyzz_1[j] + fl1_fxn * tg_yyyyzz_xyyzz_1[j];

                    tg_xyyyyzz_xxyzzz_0[j] = pb_x * tg_yyyyzz_xxyzzz_0[j] + fr * tg_yyyyzz_xxyzzz_1[j] + fl1_fxn * tg_yyyyzz_xyzzz_1[j];

                    tg_xyyyyzz_xxzzzz_0[j] = pb_x * tg_yyyyzz_xxzzzz_0[j] + fr * tg_yyyyzz_xxzzzz_1[j] + fl1_fxn * tg_yyyyzz_xzzzz_1[j];

                    tg_xyyyyzz_xyyyyy_0[j] = pb_x * tg_yyyyzz_xyyyyy_0[j] + fr * tg_yyyyzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyyyy_1[j];

                    tg_xyyyyzz_xyyyyz_0[j] = pb_x * tg_yyyyzz_xyyyyz_0[j] + fr * tg_yyyyzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyyyz_1[j];

                    tg_xyyyyzz_xyyyzz_0[j] = pb_x * tg_yyyyzz_xyyyzz_0[j] + fr * tg_yyyyzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyyzz_1[j];

                    tg_xyyyyzz_xyyzzz_0[j] = pb_x * tg_yyyyzz_xyyzzz_0[j] + fr * tg_yyyyzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyzzz_1[j];

                    tg_xyyyyzz_xyzzzz_0[j] = pb_x * tg_yyyyzz_xyzzzz_0[j] + fr * tg_yyyyzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yzzzz_1[j];

                    tg_xyyyyzz_xzzzzz_0[j] = pb_x * tg_yyyyzz_xzzzzz_0[j] + fr * tg_yyyyzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_zzzzz_1[j];

                    tg_xyyyyzz_yyyyyy_0[j] = pb_x * tg_yyyyzz_yyyyyy_0[j] + fr * tg_yyyyzz_yyyyyy_1[j];

                    tg_xyyyyzz_yyyyyz_0[j] = pb_x * tg_yyyyzz_yyyyyz_0[j] + fr * tg_yyyyzz_yyyyyz_1[j];

                    tg_xyyyyzz_yyyyzz_0[j] = pb_x * tg_yyyyzz_yyyyzz_0[j] + fr * tg_yyyyzz_yyyyzz_1[j];

                    tg_xyyyyzz_yyyzzz_0[j] = pb_x * tg_yyyyzz_yyyzzz_0[j] + fr * tg_yyyyzz_yyyzzz_1[j];

                    tg_xyyyyzz_yyzzzz_0[j] = pb_x * tg_yyyyzz_yyzzzz_0[j] + fr * tg_yyyyzz_yyzzzz_1[j];

                    tg_xyyyyzz_yzzzzz_0[j] = pb_x * tg_yyyyzz_yzzzzz_0[j] + fr * tg_yyyyzz_yzzzzz_1[j];

                    tg_xyyyyzz_zzzzzz_0[j] = pb_x * tg_yyyyzz_zzzzzz_0[j] + fr * tg_yyyyzz_zzzzzz_1[j];

                    tg_xyyyzzz_xxxxxx_0[j] = pb_x * tg_yyyzzz_xxxxxx_0[j] + fr * tg_yyyzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyzzz_xxxxx_1[j];

                    tg_xyyyzzz_xxxxxy_0[j] = pb_x * tg_yyyzzz_xxxxxy_0[j] + fr * tg_yyyzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyzzz_xxxxy_1[j];

                    tg_xyyyzzz_xxxxxz_0[j] = pb_x * tg_yyyzzz_xxxxxz_0[j] + fr * tg_yyyzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyzzz_xxxxz_1[j];

                    tg_xyyyzzz_xxxxyy_0[j] = pb_x * tg_yyyzzz_xxxxyy_0[j] + fr * tg_yyyzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyzzz_xxxyy_1[j];

                    tg_xyyyzzz_xxxxyz_0[j] = pb_x * tg_yyyzzz_xxxxyz_0[j] + fr * tg_yyyzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyzzz_xxxyz_1[j];

                    tg_xyyyzzz_xxxxzz_0[j] = pb_x * tg_yyyzzz_xxxxzz_0[j] + fr * tg_yyyzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyzzz_xxxzz_1[j];

                    tg_xyyyzzz_xxxyyy_0[j] = pb_x * tg_yyyzzz_xxxyyy_0[j] + fr * tg_yyyzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxyyy_1[j];

                    tg_xyyyzzz_xxxyyz_0[j] = pb_x * tg_yyyzzz_xxxyyz_0[j] + fr * tg_yyyzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxyyz_1[j];

                    tg_xyyyzzz_xxxyzz_0[j] = pb_x * tg_yyyzzz_xxxyzz_0[j] + fr * tg_yyyzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxyzz_1[j];

                    tg_xyyyzzz_xxxzzz_0[j] = pb_x * tg_yyyzzz_xxxzzz_0[j] + fr * tg_yyyzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxzzz_1[j];

                    tg_xyyyzzz_xxyyyy_0[j] = pb_x * tg_yyyzzz_xxyyyy_0[j] + fr * tg_yyyzzz_xxyyyy_1[j] + fl1_fxn * tg_yyyzzz_xyyyy_1[j];

                    tg_xyyyzzz_xxyyyz_0[j] = pb_x * tg_yyyzzz_xxyyyz_0[j] + fr * tg_yyyzzz_xxyyyz_1[j] + fl1_fxn * tg_yyyzzz_xyyyz_1[j];

                    tg_xyyyzzz_xxyyzz_0[j] = pb_x * tg_yyyzzz_xxyyzz_0[j] + fr * tg_yyyzzz_xxyyzz_1[j] + fl1_fxn * tg_yyyzzz_xyyzz_1[j];

                    tg_xyyyzzz_xxyzzz_0[j] = pb_x * tg_yyyzzz_xxyzzz_0[j] + fr * tg_yyyzzz_xxyzzz_1[j] + fl1_fxn * tg_yyyzzz_xyzzz_1[j];

                    tg_xyyyzzz_xxzzzz_0[j] = pb_x * tg_yyyzzz_xxzzzz_0[j] + fr * tg_yyyzzz_xxzzzz_1[j] + fl1_fxn * tg_yyyzzz_xzzzz_1[j];

                    tg_xyyyzzz_xyyyyy_0[j] = pb_x * tg_yyyzzz_xyyyyy_0[j] + fr * tg_yyyzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyyyy_1[j];

                    tg_xyyyzzz_xyyyyz_0[j] = pb_x * tg_yyyzzz_xyyyyz_0[j] + fr * tg_yyyzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyyyz_1[j];

                    tg_xyyyzzz_xyyyzz_0[j] = pb_x * tg_yyyzzz_xyyyzz_0[j] + fr * tg_yyyzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyyzz_1[j];

                    tg_xyyyzzz_xyyzzz_0[j] = pb_x * tg_yyyzzz_xyyzzz_0[j] + fr * tg_yyyzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyzzz_1[j];

                    tg_xyyyzzz_xyzzzz_0[j] = pb_x * tg_yyyzzz_xyzzzz_0[j] + fr * tg_yyyzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yzzzz_1[j];

                    tg_xyyyzzz_xzzzzz_0[j] = pb_x * tg_yyyzzz_xzzzzz_0[j] + fr * tg_yyyzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_zzzzz_1[j];

                    tg_xyyyzzz_yyyyyy_0[j] = pb_x * tg_yyyzzz_yyyyyy_0[j] + fr * tg_yyyzzz_yyyyyy_1[j];

                    tg_xyyyzzz_yyyyyz_0[j] = pb_x * tg_yyyzzz_yyyyyz_0[j] + fr * tg_yyyzzz_yyyyyz_1[j];

                    tg_xyyyzzz_yyyyzz_0[j] = pb_x * tg_yyyzzz_yyyyzz_0[j] + fr * tg_yyyzzz_yyyyzz_1[j];

                    tg_xyyyzzz_yyyzzz_0[j] = pb_x * tg_yyyzzz_yyyzzz_0[j] + fr * tg_yyyzzz_yyyzzz_1[j];

                    tg_xyyyzzz_yyzzzz_0[j] = pb_x * tg_yyyzzz_yyzzzz_0[j] + fr * tg_yyyzzz_yyzzzz_1[j];

                    tg_xyyyzzz_yzzzzz_0[j] = pb_x * tg_yyyzzz_yzzzzz_0[j] + fr * tg_yyyzzz_yzzzzz_1[j];

                    tg_xyyyzzz_zzzzzz_0[j] = pb_x * tg_yyyzzz_zzzzzz_0[j] + fr * tg_yyyzzz_zzzzzz_1[j];

                    tg_xyyzzzz_xxxxxx_0[j] = pb_x * tg_yyzzzz_xxxxxx_0[j] + fr * tg_yyzzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyzzzz_xxxxx_1[j];

                    tg_xyyzzzz_xxxxxy_0[j] = pb_x * tg_yyzzzz_xxxxxy_0[j] + fr * tg_yyzzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyzzzz_xxxxy_1[j];

                    tg_xyyzzzz_xxxxxz_0[j] = pb_x * tg_yyzzzz_xxxxxz_0[j] + fr * tg_yyzzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyzzzz_xxxxz_1[j];

                    tg_xyyzzzz_xxxxyy_0[j] = pb_x * tg_yyzzzz_xxxxyy_0[j] + fr * tg_yyzzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyzzzz_xxxyy_1[j];

                    tg_xyyzzzz_xxxxyz_0[j] = pb_x * tg_yyzzzz_xxxxyz_0[j] + fr * tg_yyzzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyzzzz_xxxyz_1[j];

                    tg_xyyzzzz_xxxxzz_0[j] = pb_x * tg_yyzzzz_xxxxzz_0[j] + fr * tg_yyzzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyzzzz_xxxzz_1[j];

                    tg_xyyzzzz_xxxyyy_0[j] = pb_x * tg_yyzzzz_xxxyyy_0[j] + fr * tg_yyzzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxyyy_1[j];

                    tg_xyyzzzz_xxxyyz_0[j] = pb_x * tg_yyzzzz_xxxyyz_0[j] + fr * tg_yyzzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxyyz_1[j];

                    tg_xyyzzzz_xxxyzz_0[j] = pb_x * tg_yyzzzz_xxxyzz_0[j] + fr * tg_yyzzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxyzz_1[j];

                    tg_xyyzzzz_xxxzzz_0[j] = pb_x * tg_yyzzzz_xxxzzz_0[j] + fr * tg_yyzzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxzzz_1[j];

                    tg_xyyzzzz_xxyyyy_0[j] = pb_x * tg_yyzzzz_xxyyyy_0[j] + fr * tg_yyzzzz_xxyyyy_1[j] + fl1_fxn * tg_yyzzzz_xyyyy_1[j];

                    tg_xyyzzzz_xxyyyz_0[j] = pb_x * tg_yyzzzz_xxyyyz_0[j] + fr * tg_yyzzzz_xxyyyz_1[j] + fl1_fxn * tg_yyzzzz_xyyyz_1[j];

                    tg_xyyzzzz_xxyyzz_0[j] = pb_x * tg_yyzzzz_xxyyzz_0[j] + fr * tg_yyzzzz_xxyyzz_1[j] + fl1_fxn * tg_yyzzzz_xyyzz_1[j];

                    tg_xyyzzzz_xxyzzz_0[j] = pb_x * tg_yyzzzz_xxyzzz_0[j] + fr * tg_yyzzzz_xxyzzz_1[j] + fl1_fxn * tg_yyzzzz_xyzzz_1[j];

                    tg_xyyzzzz_xxzzzz_0[j] = pb_x * tg_yyzzzz_xxzzzz_0[j] + fr * tg_yyzzzz_xxzzzz_1[j] + fl1_fxn * tg_yyzzzz_xzzzz_1[j];

                    tg_xyyzzzz_xyyyyy_0[j] = pb_x * tg_yyzzzz_xyyyyy_0[j] + fr * tg_yyzzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyyyy_1[j];

                    tg_xyyzzzz_xyyyyz_0[j] = pb_x * tg_yyzzzz_xyyyyz_0[j] + fr * tg_yyzzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyyyz_1[j];

                    tg_xyyzzzz_xyyyzz_0[j] = pb_x * tg_yyzzzz_xyyyzz_0[j] + fr * tg_yyzzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyyzz_1[j];

                    tg_xyyzzzz_xyyzzz_0[j] = pb_x * tg_yyzzzz_xyyzzz_0[j] + fr * tg_yyzzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyzzz_1[j];

                    tg_xyyzzzz_xyzzzz_0[j] = pb_x * tg_yyzzzz_xyzzzz_0[j] + fr * tg_yyzzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yzzzz_1[j];

                    tg_xyyzzzz_xzzzzz_0[j] = pb_x * tg_yyzzzz_xzzzzz_0[j] + fr * tg_yyzzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_zzzzz_1[j];

                    tg_xyyzzzz_yyyyyy_0[j] = pb_x * tg_yyzzzz_yyyyyy_0[j] + fr * tg_yyzzzz_yyyyyy_1[j];

                    tg_xyyzzzz_yyyyyz_0[j] = pb_x * tg_yyzzzz_yyyyyz_0[j] + fr * tg_yyzzzz_yyyyyz_1[j];

                    tg_xyyzzzz_yyyyzz_0[j] = pb_x * tg_yyzzzz_yyyyzz_0[j] + fr * tg_yyzzzz_yyyyzz_1[j];

                    tg_xyyzzzz_yyyzzz_0[j] = pb_x * tg_yyzzzz_yyyzzz_0[j] + fr * tg_yyzzzz_yyyzzz_1[j];

                    tg_xyyzzzz_yyzzzz_0[j] = pb_x * tg_yyzzzz_yyzzzz_0[j] + fr * tg_yyzzzz_yyzzzz_1[j];

                    tg_xyyzzzz_yzzzzz_0[j] = pb_x * tg_yyzzzz_yzzzzz_0[j] + fr * tg_yyzzzz_yzzzzz_1[j];

                    tg_xyyzzzz_zzzzzz_0[j] = pb_x * tg_yyzzzz_zzzzzz_0[j] + fr * tg_yyzzzz_zzzzzz_1[j];

                    tg_xyzzzzz_xxxxxx_0[j] = pb_x * tg_yzzzzz_xxxxxx_0[j] + fr * tg_yzzzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yzzzzz_xxxxx_1[j];

                    tg_xyzzzzz_xxxxxy_0[j] = pb_x * tg_yzzzzz_xxxxxy_0[j] + fr * tg_yzzzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yzzzzz_xxxxy_1[j];

                    tg_xyzzzzz_xxxxxz_0[j] = pb_x * tg_yzzzzz_xxxxxz_0[j] + fr * tg_yzzzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yzzzzz_xxxxz_1[j];

                    tg_xyzzzzz_xxxxyy_0[j] = pb_x * tg_yzzzzz_xxxxyy_0[j] + fr * tg_yzzzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yzzzzz_xxxyy_1[j];

                    tg_xyzzzzz_xxxxyz_0[j] = pb_x * tg_yzzzzz_xxxxyz_0[j] + fr * tg_yzzzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yzzzzz_xxxyz_1[j];

                    tg_xyzzzzz_xxxxzz_0[j] = pb_x * tg_yzzzzz_xxxxzz_0[j] + fr * tg_yzzzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yzzzzz_xxxzz_1[j];

                    tg_xyzzzzz_xxxyyy_0[j] = pb_x * tg_yzzzzz_xxxyyy_0[j] + fr * tg_yzzzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxyyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_735_826(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (735,826)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yyyyyy_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 588); 

                auto tg_yyyyyy_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 589); 

                auto tg_yyyyyy_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 590); 

                auto tg_yyyyyy_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 591); 

                auto tg_yyyyyy_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 592); 

                auto tg_yyyyyy_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 593); 

                auto tg_yyyyyy_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 594); 

                auto tg_yyyyyy_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 595); 

                auto tg_yyyyyy_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 596); 

                auto tg_yyyyyy_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 597); 

                auto tg_yyyyyy_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 598); 

                auto tg_yyyyyy_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 599); 

                auto tg_yyyyyy_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 600); 

                auto tg_yyyyyy_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 601); 

                auto tg_yyyyyy_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 602); 

                auto tg_yyyyyy_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 603); 

                auto tg_yyyyyy_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 604); 

                auto tg_yyyyyy_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 605); 

                auto tg_yyyyyy_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 606); 

                auto tg_yyyyyy_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 607); 

                auto tg_yyyyyy_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 608); 

                auto tg_yyyyyy_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 609); 

                auto tg_yyyyyy_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 610); 

                auto tg_yyyyyy_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 611); 

                auto tg_yyyyyy_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 612); 

                auto tg_yyyyyy_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 613); 

                auto tg_yyyyyy_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 614); 

                auto tg_yyyyyy_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 615); 

                auto tg_yyyyyz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 616); 

                auto tg_yyyyyz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 617); 

                auto tg_yyyyyz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 618); 

                auto tg_yyyyyz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 619); 

                auto tg_yyyyyz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 620); 

                auto tg_yyyyyz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 621); 

                auto tg_yyyyyz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 622); 

                auto tg_yyyyyz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 623); 

                auto tg_yyyyyz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 624); 

                auto tg_yyyyyz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 625); 

                auto tg_yyyyyz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 626); 

                auto tg_yyyyyz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 627); 

                auto tg_yyyyyz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 628); 

                auto tg_yyyyyz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 629); 

                auto tg_yzzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 735); 

                auto tg_yzzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 736); 

                auto tg_yzzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 737); 

                auto tg_yzzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 738); 

                auto tg_yzzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 739); 

                auto tg_yzzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 740); 

                auto tg_yzzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 741); 

                auto tg_yzzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 742); 

                auto tg_yzzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 743); 

                auto tg_yzzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 744); 

                auto tg_yzzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 745); 

                auto tg_yzzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 746); 

                auto tg_yzzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 747); 

                auto tg_yzzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 748); 

                auto tg_yzzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 749); 

                auto tg_yzzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 750); 

                auto tg_yzzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 751); 

                auto tg_yzzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 752); 

                auto tg_yzzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 753); 

                auto tg_yzzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 754); 

                auto tg_yzzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 755); 

                auto tg_zzzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 756); 

                auto tg_zzzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 757); 

                auto tg_zzzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 758); 

                auto tg_zzzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 759); 

                auto tg_zzzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 760); 

                auto tg_zzzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 761); 

                auto tg_zzzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 762); 

                auto tg_zzzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 763); 

                auto tg_zzzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 764); 

                auto tg_zzzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 765); 

                auto tg_zzzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 766); 

                auto tg_zzzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 767); 

                auto tg_zzzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 768); 

                auto tg_zzzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 769); 

                auto tg_zzzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 770); 

                auto tg_zzzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 771); 

                auto tg_zzzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 772); 

                auto tg_zzzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 773); 

                auto tg_zzzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 774); 

                auto tg_zzzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 775); 

                auto tg_zzzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 776); 

                auto tg_zzzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 777); 

                auto tg_zzzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 778); 

                auto tg_zzzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 779); 

                auto tg_zzzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 780); 

                auto tg_zzzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 781); 

                auto tg_zzzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 782); 

                auto tg_zzzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 783); 

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

                // set up pointers to integrals

                auto tg_xyzzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 735); 

                auto tg_xyzzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 736); 

                auto tg_xyzzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 737); 

                auto tg_xyzzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 738); 

                auto tg_xyzzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 739); 

                auto tg_xyzzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 740); 

                auto tg_xyzzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 741); 

                auto tg_xyzzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 742); 

                auto tg_xyzzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 743); 

                auto tg_xyzzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 744); 

                auto tg_xyzzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 745); 

                auto tg_xyzzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 746); 

                auto tg_xyzzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 747); 

                auto tg_xyzzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 748); 

                auto tg_xyzzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 749); 

                auto tg_xyzzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 750); 

                auto tg_xyzzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 751); 

                auto tg_xyzzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 752); 

                auto tg_xyzzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 753); 

                auto tg_xyzzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 754); 

                auto tg_xyzzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 755); 

                auto tg_xzzzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 756); 

                auto tg_xzzzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 757); 

                auto tg_xzzzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 758); 

                auto tg_xzzzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 759); 

                auto tg_xzzzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 760); 

                auto tg_xzzzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 761); 

                auto tg_xzzzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 762); 

                auto tg_xzzzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 763); 

                auto tg_xzzzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 764); 

                auto tg_xzzzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 765); 

                auto tg_xzzzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 766); 

                auto tg_xzzzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 767); 

                auto tg_xzzzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 768); 

                auto tg_xzzzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 769); 

                auto tg_xzzzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 770); 

                auto tg_xzzzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 771); 

                auto tg_xzzzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 772); 

                auto tg_xzzzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 773); 

                auto tg_xzzzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 774); 

                auto tg_xzzzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 775); 

                auto tg_xzzzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 776); 

                auto tg_xzzzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 777); 

                auto tg_xzzzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 778); 

                auto tg_xzzzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 779); 

                auto tg_xzzzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 780); 

                auto tg_xzzzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 781); 

                auto tg_xzzzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 782); 

                auto tg_xzzzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 783); 

                auto tg_yyyyyyy_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 784); 

                auto tg_yyyyyyy_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 785); 

                auto tg_yyyyyyy_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 786); 

                auto tg_yyyyyyy_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 787); 

                auto tg_yyyyyyy_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 788); 

                auto tg_yyyyyyy_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 789); 

                auto tg_yyyyyyy_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 790); 

                auto tg_yyyyyyy_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 791); 

                auto tg_yyyyyyy_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 792); 

                auto tg_yyyyyyy_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 793); 

                auto tg_yyyyyyy_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 794); 

                auto tg_yyyyyyy_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 795); 

                auto tg_yyyyyyy_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 796); 

                auto tg_yyyyyyy_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 797); 

                auto tg_yyyyyyy_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 798); 

                auto tg_yyyyyyy_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 799); 

                auto tg_yyyyyyy_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 800); 

                auto tg_yyyyyyy_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 801); 

                auto tg_yyyyyyy_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 802); 

                auto tg_yyyyyyy_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 803); 

                auto tg_yyyyyyy_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 804); 

                auto tg_yyyyyyy_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 805); 

                auto tg_yyyyyyy_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 806); 

                auto tg_yyyyyyy_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 807); 

                auto tg_yyyyyyy_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 808); 

                auto tg_yyyyyyy_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 809); 

                auto tg_yyyyyyy_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 810); 

                auto tg_yyyyyyy_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 811); 

                auto tg_yyyyyyz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 812); 

                auto tg_yyyyyyz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 813); 

                auto tg_yyyyyyz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 814); 

                auto tg_yyyyyyz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 815); 

                auto tg_yyyyyyz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 816); 

                auto tg_yyyyyyz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 817); 

                auto tg_yyyyyyz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 818); 

                auto tg_yyyyyyz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 819); 

                auto tg_yyyyyyz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 820); 

                auto tg_yyyyyyz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 821); 

                auto tg_yyyyyyz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 822); 

                auto tg_yyyyyyz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 823); 

                auto tg_yyyyyyz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 824); 

                auto tg_yyyyyyz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 825); 

                // Batch of Integrals (735,826)

                #pragma omp simd aligned(fxn, fza, tg_xyzzzzz_xxxyyz_0, tg_xyzzzzz_xxxyzz_0, \
                                         tg_xyzzzzz_xxxzzz_0, tg_xyzzzzz_xxyyyy_0, tg_xyzzzzz_xxyyyz_0, tg_xyzzzzz_xxyyzz_0, \
                                         tg_xyzzzzz_xxyzzz_0, tg_xyzzzzz_xxzzzz_0, tg_xyzzzzz_xyyyyy_0, tg_xyzzzzz_xyyyyz_0, \
                                         tg_xyzzzzz_xyyyzz_0, tg_xyzzzzz_xyyzzz_0, tg_xyzzzzz_xyzzzz_0, tg_xyzzzzz_xzzzzz_0, \
                                         tg_xyzzzzz_yyyyyy_0, tg_xyzzzzz_yyyyyz_0, tg_xyzzzzz_yyyyzz_0, tg_xyzzzzz_yyyzzz_0, \
                                         tg_xyzzzzz_yyzzzz_0, tg_xyzzzzz_yzzzzz_0, tg_xyzzzzz_zzzzzz_0, tg_xzzzzzz_xxxxxx_0, \
                                         tg_xzzzzzz_xxxxxy_0, tg_xzzzzzz_xxxxxz_0, tg_xzzzzzz_xxxxyy_0, tg_xzzzzzz_xxxxyz_0, \
                                         tg_xzzzzzz_xxxxzz_0, tg_xzzzzzz_xxxyyy_0, tg_xzzzzzz_xxxyyz_0, tg_xzzzzzz_xxxyzz_0, \
                                         tg_xzzzzzz_xxxzzz_0, tg_xzzzzzz_xxyyyy_0, tg_xzzzzzz_xxyyyz_0, tg_xzzzzzz_xxyyzz_0, \
                                         tg_xzzzzzz_xxyzzz_0, tg_xzzzzzz_xxzzzz_0, tg_xzzzzzz_xyyyyy_0, tg_xzzzzzz_xyyyyz_0, \
                                         tg_xzzzzzz_xyyyzz_0, tg_xzzzzzz_xyyzzz_0, tg_xzzzzzz_xyzzzz_0, tg_xzzzzzz_xzzzzz_0, \
                                         tg_xzzzzzz_yyyyyy_0, tg_xzzzzzz_yyyyyz_0, tg_xzzzzzz_yyyyzz_0, tg_xzzzzzz_yyyzzz_0, \
                                         tg_xzzzzzz_yyzzzz_0, tg_xzzzzzz_yzzzzz_0, tg_xzzzzzz_zzzzzz_0, tg_yyyyy_xxxxxx_0, \
                                         tg_yyyyy_xxxxxx_1, tg_yyyyy_xxxxxy_0, tg_yyyyy_xxxxxy_1, tg_yyyyy_xxxxxz_0, \
                                         tg_yyyyy_xxxxxz_1, tg_yyyyy_xxxxyy_0, tg_yyyyy_xxxxyy_1, tg_yyyyy_xxxxyz_0, \
                                         tg_yyyyy_xxxxyz_1, tg_yyyyy_xxxxzz_0, tg_yyyyy_xxxxzz_1, tg_yyyyy_xxxyyy_0, \
                                         tg_yyyyy_xxxyyy_1, tg_yyyyy_xxxyyz_0, tg_yyyyy_xxxyyz_1, tg_yyyyy_xxxyzz_0, \
                                         tg_yyyyy_xxxyzz_1, tg_yyyyy_xxxzzz_0, tg_yyyyy_xxxzzz_1, tg_yyyyy_xxyyyy_0, \
                                         tg_yyyyy_xxyyyy_1, tg_yyyyy_xxyyyz_0, tg_yyyyy_xxyyyz_1, tg_yyyyy_xxyyzz_0, \
                                         tg_yyyyy_xxyyzz_1, tg_yyyyy_xxyzzz_0, tg_yyyyy_xxyzzz_1, tg_yyyyy_xxzzzz_0, \
                                         tg_yyyyy_xxzzzz_1, tg_yyyyy_xyyyyy_0, tg_yyyyy_xyyyyy_1, tg_yyyyy_xyyyyz_0, \
                                         tg_yyyyy_xyyyyz_1, tg_yyyyy_xyyyzz_0, tg_yyyyy_xyyyzz_1, tg_yyyyy_xyyzzz_0, \
                                         tg_yyyyy_xyyzzz_1, tg_yyyyy_xyzzzz_0, tg_yyyyy_xyzzzz_1, tg_yyyyy_xzzzzz_0, \
                                         tg_yyyyy_xzzzzz_1, tg_yyyyy_yyyyyy_0, tg_yyyyy_yyyyyy_1, tg_yyyyy_yyyyyz_0, \
                                         tg_yyyyy_yyyyyz_1, tg_yyyyy_yyyyzz_0, tg_yyyyy_yyyyzz_1, tg_yyyyy_yyyzzz_0, \
                                         tg_yyyyy_yyyzzz_1, tg_yyyyy_yyzzzz_0, tg_yyyyy_yyzzzz_1, tg_yyyyy_yzzzzz_0, \
                                         tg_yyyyy_yzzzzz_1, tg_yyyyy_zzzzzz_0, tg_yyyyy_zzzzzz_1, tg_yyyyyy_xxxxx_1, \
                                         tg_yyyyyy_xxxxxx_0, tg_yyyyyy_xxxxxx_1, tg_yyyyyy_xxxxxy_0, tg_yyyyyy_xxxxxy_1, \
                                         tg_yyyyyy_xxxxxz_0, tg_yyyyyy_xxxxxz_1, tg_yyyyyy_xxxxy_1, tg_yyyyyy_xxxxyy_0, \
                                         tg_yyyyyy_xxxxyy_1, tg_yyyyyy_xxxxyz_0, tg_yyyyyy_xxxxyz_1, tg_yyyyyy_xxxxz_1, \
                                         tg_yyyyyy_xxxxzz_0, tg_yyyyyy_xxxxzz_1, tg_yyyyyy_xxxyy_1, tg_yyyyyy_xxxyyy_0, \
                                         tg_yyyyyy_xxxyyy_1, tg_yyyyyy_xxxyyz_0, tg_yyyyyy_xxxyyz_1, tg_yyyyyy_xxxyz_1, \
                                         tg_yyyyyy_xxxyzz_0, tg_yyyyyy_xxxyzz_1, tg_yyyyyy_xxxzz_1, tg_yyyyyy_xxxzzz_0, \
                                         tg_yyyyyy_xxxzzz_1, tg_yyyyyy_xxyyy_1, tg_yyyyyy_xxyyyy_0, tg_yyyyyy_xxyyyy_1, \
                                         tg_yyyyyy_xxyyyz_0, tg_yyyyyy_xxyyyz_1, tg_yyyyyy_xxyyz_1, tg_yyyyyy_xxyyzz_0, \
                                         tg_yyyyyy_xxyyzz_1, tg_yyyyyy_xxyzz_1, tg_yyyyyy_xxyzzz_0, tg_yyyyyy_xxyzzz_1, \
                                         tg_yyyyyy_xxzzz_1, tg_yyyyyy_xxzzzz_0, tg_yyyyyy_xxzzzz_1, tg_yyyyyy_xyyyy_1, \
                                         tg_yyyyyy_xyyyyy_0, tg_yyyyyy_xyyyyy_1, tg_yyyyyy_xyyyyz_0, tg_yyyyyy_xyyyyz_1, \
                                         tg_yyyyyy_xyyyz_1, tg_yyyyyy_xyyyzz_0, tg_yyyyyy_xyyyzz_1, tg_yyyyyy_xyyzz_1, \
                                         tg_yyyyyy_xyyzzz_0, tg_yyyyyy_xyyzzz_1, tg_yyyyyy_xyzzz_1, tg_yyyyyy_xyzzzz_0, \
                                         tg_yyyyyy_xyzzzz_1, tg_yyyyyy_xzzzz_1, tg_yyyyyy_xzzzzz_0, tg_yyyyyy_xzzzzz_1, \
                                         tg_yyyyyy_yyyyy_1, tg_yyyyyy_yyyyyy_0, tg_yyyyyy_yyyyyy_1, tg_yyyyyy_yyyyyz_0, \
                                         tg_yyyyyy_yyyyyz_1, tg_yyyyyy_yyyyz_1, tg_yyyyyy_yyyyzz_0, tg_yyyyyy_yyyyzz_1, \
                                         tg_yyyyyy_yyyzz_1, tg_yyyyyy_yyyzzz_0, tg_yyyyyy_yyyzzz_1, tg_yyyyyy_yyzzz_1, \
                                         tg_yyyyyy_yyzzzz_0, tg_yyyyyy_yyzzzz_1, tg_yyyyyy_yzzzz_1, tg_yyyyyy_yzzzzz_0, \
                                         tg_yyyyyy_yzzzzz_1, tg_yyyyyy_zzzzz_1, tg_yyyyyy_zzzzzz_0, tg_yyyyyy_zzzzzz_1, \
                                         tg_yyyyyyy_xxxxxx_0, tg_yyyyyyy_xxxxxy_0, tg_yyyyyyy_xxxxxz_0, tg_yyyyyyy_xxxxyy_0, \
                                         tg_yyyyyyy_xxxxyz_0, tg_yyyyyyy_xxxxzz_0, tg_yyyyyyy_xxxyyy_0, tg_yyyyyyy_xxxyyz_0, \
                                         tg_yyyyyyy_xxxyzz_0, tg_yyyyyyy_xxxzzz_0, tg_yyyyyyy_xxyyyy_0, tg_yyyyyyy_xxyyyz_0, \
                                         tg_yyyyyyy_xxyyzz_0, tg_yyyyyyy_xxyzzz_0, tg_yyyyyyy_xxzzzz_0, tg_yyyyyyy_xyyyyy_0, \
                                         tg_yyyyyyy_xyyyyz_0, tg_yyyyyyy_xyyyzz_0, tg_yyyyyyy_xyyzzz_0, tg_yyyyyyy_xyzzzz_0, \
                                         tg_yyyyyyy_xzzzzz_0, tg_yyyyyyy_yyyyyy_0, tg_yyyyyyy_yyyyyz_0, tg_yyyyyyy_yyyyzz_0, \
                                         tg_yyyyyyy_yyyzzz_0, tg_yyyyyyy_yyzzzz_0, tg_yyyyyyy_yzzzzz_0, tg_yyyyyyy_zzzzzz_0, \
                                         tg_yyyyyyz_xxxxxx_0, tg_yyyyyyz_xxxxxy_0, tg_yyyyyyz_xxxxxz_0, tg_yyyyyyz_xxxxyy_0, \
                                         tg_yyyyyyz_xxxxyz_0, tg_yyyyyyz_xxxxzz_0, tg_yyyyyyz_xxxyyy_0, tg_yyyyyyz_xxxyyz_0, \
                                         tg_yyyyyyz_xxxyzz_0, tg_yyyyyyz_xxxzzz_0, tg_yyyyyyz_xxyyyy_0, tg_yyyyyyz_xxyyyz_0, \
                                         tg_yyyyyyz_xxyyzz_0, tg_yyyyyyz_xxyzzz_0, tg_yyyyyz_xxxxx_1, tg_yyyyyz_xxxxxx_0, \
                                         tg_yyyyyz_xxxxxx_1, tg_yyyyyz_xxxxxy_0, tg_yyyyyz_xxxxxy_1, tg_yyyyyz_xxxxxz_0, \
                                         tg_yyyyyz_xxxxxz_1, tg_yyyyyz_xxxxy_1, tg_yyyyyz_xxxxyy_0, tg_yyyyyz_xxxxyy_1, \
                                         tg_yyyyyz_xxxxyz_0, tg_yyyyyz_xxxxyz_1, tg_yyyyyz_xxxxz_1, tg_yyyyyz_xxxxzz_0, \
                                         tg_yyyyyz_xxxxzz_1, tg_yyyyyz_xxxyy_1, tg_yyyyyz_xxxyyy_0, tg_yyyyyz_xxxyyy_1, \
                                         tg_yyyyyz_xxxyyz_0, tg_yyyyyz_xxxyyz_1, tg_yyyyyz_xxxyz_1, tg_yyyyyz_xxxyzz_0, \
                                         tg_yyyyyz_xxxyzz_1, tg_yyyyyz_xxxzz_1, tg_yyyyyz_xxxzzz_0, tg_yyyyyz_xxxzzz_1, \
                                         tg_yyyyyz_xxyyy_1, tg_yyyyyz_xxyyyy_0, tg_yyyyyz_xxyyyy_1, tg_yyyyyz_xxyyyz_0, \
                                         tg_yyyyyz_xxyyyz_1, tg_yyyyyz_xxyyz_1, tg_yyyyyz_xxyyzz_0, tg_yyyyyz_xxyyzz_1, \
                                         tg_yyyyyz_xxyzz_1, tg_yyyyyz_xxyzzz_0, tg_yyyyyz_xxyzzz_1, tg_yyyyyz_xxzzz_1, \
                                         tg_yyyyz_xxxxxx_0, tg_yyyyz_xxxxxx_1, tg_yyyyz_xxxxxy_0, tg_yyyyz_xxxxxy_1, \
                                         tg_yyyyz_xxxxxz_0, tg_yyyyz_xxxxxz_1, tg_yyyyz_xxxxyy_0, tg_yyyyz_xxxxyy_1, \
                                         tg_yyyyz_xxxxyz_0, tg_yyyyz_xxxxyz_1, tg_yyyyz_xxxxzz_0, tg_yyyyz_xxxxzz_1, \
                                         tg_yyyyz_xxxyyy_0, tg_yyyyz_xxxyyy_1, tg_yyyyz_xxxyyz_0, tg_yyyyz_xxxyyz_1, \
                                         tg_yyyyz_xxxyzz_0, tg_yyyyz_xxxyzz_1, tg_yyyyz_xxxzzz_0, tg_yyyyz_xxxzzz_1, \
                                         tg_yyyyz_xxyyyy_0, tg_yyyyz_xxyyyy_1, tg_yyyyz_xxyyyz_0, tg_yyyyz_xxyyyz_1, \
                                         tg_yyyyz_xxyyzz_0, tg_yyyyz_xxyyzz_1, tg_yyyyz_xxyzzz_0, tg_yyyyz_xxyzzz_1, \
                                         tg_yzzzzz_xxxyyz_0, tg_yzzzzz_xxxyyz_1, tg_yzzzzz_xxxyzz_0, tg_yzzzzz_xxxyzz_1, \
                                         tg_yzzzzz_xxxzzz_0, tg_yzzzzz_xxxzzz_1, tg_yzzzzz_xxyyyy_0, tg_yzzzzz_xxyyyy_1, \
                                         tg_yzzzzz_xxyyyz_0, tg_yzzzzz_xxyyyz_1, tg_yzzzzz_xxyyz_1, tg_yzzzzz_xxyyzz_0, \
                                         tg_yzzzzz_xxyyzz_1, tg_yzzzzz_xxyzz_1, tg_yzzzzz_xxyzzz_0, tg_yzzzzz_xxyzzz_1, \
                                         tg_yzzzzz_xxzzz_1, tg_yzzzzz_xxzzzz_0, tg_yzzzzz_xxzzzz_1, tg_yzzzzz_xyyyy_1, \
                                         tg_yzzzzz_xyyyyy_0, tg_yzzzzz_xyyyyy_1, tg_yzzzzz_xyyyyz_0, tg_yzzzzz_xyyyyz_1, \
                                         tg_yzzzzz_xyyyz_1, tg_yzzzzz_xyyyzz_0, tg_yzzzzz_xyyyzz_1, tg_yzzzzz_xyyzz_1, \
                                         tg_yzzzzz_xyyzzz_0, tg_yzzzzz_xyyzzz_1, tg_yzzzzz_xyzzz_1, tg_yzzzzz_xyzzzz_0, \
                                         tg_yzzzzz_xyzzzz_1, tg_yzzzzz_xzzzz_1, tg_yzzzzz_xzzzzz_0, tg_yzzzzz_xzzzzz_1, \
                                         tg_yzzzzz_yyyyy_1, tg_yzzzzz_yyyyyy_0, tg_yzzzzz_yyyyyy_1, tg_yzzzzz_yyyyyz_0, \
                                         tg_yzzzzz_yyyyyz_1, tg_yzzzzz_yyyyz_1, tg_yzzzzz_yyyyzz_0, tg_yzzzzz_yyyyzz_1, \
                                         tg_yzzzzz_yyyzz_1, tg_yzzzzz_yyyzzz_0, tg_yzzzzz_yyyzzz_1, tg_yzzzzz_yyzzz_1, \
                                         tg_yzzzzz_yyzzzz_0, tg_yzzzzz_yyzzzz_1, tg_yzzzzz_yzzzz_1, tg_yzzzzz_yzzzzz_0, \
                                         tg_yzzzzz_yzzzzz_1, tg_yzzzzz_zzzzz_1, tg_yzzzzz_zzzzzz_0, tg_yzzzzz_zzzzzz_1, \
                                         tg_zzzzzz_xxxxx_1, tg_zzzzzz_xxxxxx_0, tg_zzzzzz_xxxxxx_1, tg_zzzzzz_xxxxxy_0, \
                                         tg_zzzzzz_xxxxxy_1, tg_zzzzzz_xxxxxz_0, tg_zzzzzz_xxxxxz_1, tg_zzzzzz_xxxxy_1, \
                                         tg_zzzzzz_xxxxyy_0, tg_zzzzzz_xxxxyy_1, tg_zzzzzz_xxxxyz_0, tg_zzzzzz_xxxxyz_1, \
                                         tg_zzzzzz_xxxxz_1, tg_zzzzzz_xxxxzz_0, tg_zzzzzz_xxxxzz_1, tg_zzzzzz_xxxyy_1, \
                                         tg_zzzzzz_xxxyyy_0, tg_zzzzzz_xxxyyy_1, tg_zzzzzz_xxxyyz_0, tg_zzzzzz_xxxyyz_1, \
                                         tg_zzzzzz_xxxyz_1, tg_zzzzzz_xxxyzz_0, tg_zzzzzz_xxxyzz_1, tg_zzzzzz_xxxzz_1, \
                                         tg_zzzzzz_xxxzzz_0, tg_zzzzzz_xxxzzz_1, tg_zzzzzz_xxyyy_1, tg_zzzzzz_xxyyyy_0, \
                                         tg_zzzzzz_xxyyyy_1, tg_zzzzzz_xxyyyz_0, tg_zzzzzz_xxyyyz_1, tg_zzzzzz_xxyyz_1, \
                                         tg_zzzzzz_xxyyzz_0, tg_zzzzzz_xxyyzz_1, tg_zzzzzz_xxyzz_1, tg_zzzzzz_xxyzzz_0, \
                                         tg_zzzzzz_xxyzzz_1, tg_zzzzzz_xxzzz_1, tg_zzzzzz_xxzzzz_0, tg_zzzzzz_xxzzzz_1, \
                                         tg_zzzzzz_xyyyy_1, tg_zzzzzz_xyyyyy_0, tg_zzzzzz_xyyyyy_1, tg_zzzzzz_xyyyyz_0, \
                                         tg_zzzzzz_xyyyyz_1, tg_zzzzzz_xyyyz_1, tg_zzzzzz_xyyyzz_0, tg_zzzzzz_xyyyzz_1, \
                                         tg_zzzzzz_xyyzz_1, tg_zzzzzz_xyyzzz_0, tg_zzzzzz_xyyzzz_1, tg_zzzzzz_xyzzz_1, \
                                         tg_zzzzzz_xyzzzz_0, tg_zzzzzz_xyzzzz_1, tg_zzzzzz_xzzzz_1, tg_zzzzzz_xzzzzz_0, \
                                         tg_zzzzzz_xzzzzz_1, tg_zzzzzz_yyyyy_1, tg_zzzzzz_yyyyyy_0, tg_zzzzzz_yyyyyy_1, \
                                         tg_zzzzzz_yyyyyz_0, tg_zzzzzz_yyyyyz_1, tg_zzzzzz_yyyyz_1, tg_zzzzzz_yyyyzz_0, \
                                         tg_zzzzzz_yyyyzz_1, tg_zzzzzz_yyyzz_1, tg_zzzzzz_yyyzzz_0, tg_zzzzzz_yyyzzz_1, \
                                         tg_zzzzzz_yyzzz_1, tg_zzzzzz_yyzzzz_0, tg_zzzzzz_yyzzzz_1, tg_zzzzzz_yzzzz_1, \
                                         tg_zzzzzz_yzzzzz_0, tg_zzzzzz_yzzzzz_1, tg_zzzzzz_zzzzz_1, tg_zzzzzz_zzzzzz_0, \
                                         tg_zzzzzz_zzzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyzzzzz_xxxyyz_0[j] = pb_x * tg_yzzzzz_xxxyyz_0[j] + fr * tg_yzzzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxyyz_1[j];

                    tg_xyzzzzz_xxxyzz_0[j] = pb_x * tg_yzzzzz_xxxyzz_0[j] + fr * tg_yzzzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxyzz_1[j];

                    tg_xyzzzzz_xxxzzz_0[j] = pb_x * tg_yzzzzz_xxxzzz_0[j] + fr * tg_yzzzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxzzz_1[j];

                    tg_xyzzzzz_xxyyyy_0[j] = pb_x * tg_yzzzzz_xxyyyy_0[j] + fr * tg_yzzzzz_xxyyyy_1[j] + fl1_fxn * tg_yzzzzz_xyyyy_1[j];

                    tg_xyzzzzz_xxyyyz_0[j] = pb_x * tg_yzzzzz_xxyyyz_0[j] + fr * tg_yzzzzz_xxyyyz_1[j] + fl1_fxn * tg_yzzzzz_xyyyz_1[j];

                    tg_xyzzzzz_xxyyzz_0[j] = pb_x * tg_yzzzzz_xxyyzz_0[j] + fr * tg_yzzzzz_xxyyzz_1[j] + fl1_fxn * tg_yzzzzz_xyyzz_1[j];

                    tg_xyzzzzz_xxyzzz_0[j] = pb_x * tg_yzzzzz_xxyzzz_0[j] + fr * tg_yzzzzz_xxyzzz_1[j] + fl1_fxn * tg_yzzzzz_xyzzz_1[j];

                    tg_xyzzzzz_xxzzzz_0[j] = pb_x * tg_yzzzzz_xxzzzz_0[j] + fr * tg_yzzzzz_xxzzzz_1[j] + fl1_fxn * tg_yzzzzz_xzzzz_1[j];

                    tg_xyzzzzz_xyyyyy_0[j] = pb_x * tg_yzzzzz_xyyyyy_0[j] + fr * tg_yzzzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyyyy_1[j];

                    tg_xyzzzzz_xyyyyz_0[j] = pb_x * tg_yzzzzz_xyyyyz_0[j] + fr * tg_yzzzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyyyz_1[j];

                    tg_xyzzzzz_xyyyzz_0[j] = pb_x * tg_yzzzzz_xyyyzz_0[j] + fr * tg_yzzzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyyzz_1[j];

                    tg_xyzzzzz_xyyzzz_0[j] = pb_x * tg_yzzzzz_xyyzzz_0[j] + fr * tg_yzzzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyzzz_1[j];

                    tg_xyzzzzz_xyzzzz_0[j] = pb_x * tg_yzzzzz_xyzzzz_0[j] + fr * tg_yzzzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yzzzz_1[j];

                    tg_xyzzzzz_xzzzzz_0[j] = pb_x * tg_yzzzzz_xzzzzz_0[j] + fr * tg_yzzzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_zzzzz_1[j];

                    tg_xyzzzzz_yyyyyy_0[j] = pb_x * tg_yzzzzz_yyyyyy_0[j] + fr * tg_yzzzzz_yyyyyy_1[j];

                    tg_xyzzzzz_yyyyyz_0[j] = pb_x * tg_yzzzzz_yyyyyz_0[j] + fr * tg_yzzzzz_yyyyyz_1[j];

                    tg_xyzzzzz_yyyyzz_0[j] = pb_x * tg_yzzzzz_yyyyzz_0[j] + fr * tg_yzzzzz_yyyyzz_1[j];

                    tg_xyzzzzz_yyyzzz_0[j] = pb_x * tg_yzzzzz_yyyzzz_0[j] + fr * tg_yzzzzz_yyyzzz_1[j];

                    tg_xyzzzzz_yyzzzz_0[j] = pb_x * tg_yzzzzz_yyzzzz_0[j] + fr * tg_yzzzzz_yyzzzz_1[j];

                    tg_xyzzzzz_yzzzzz_0[j] = pb_x * tg_yzzzzz_yzzzzz_0[j] + fr * tg_yzzzzz_yzzzzz_1[j];

                    tg_xyzzzzz_zzzzzz_0[j] = pb_x * tg_yzzzzz_zzzzzz_0[j] + fr * tg_yzzzzz_zzzzzz_1[j];

                    tg_xzzzzzz_xxxxxx_0[j] = pb_x * tg_zzzzzz_xxxxxx_0[j] + fr * tg_zzzzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_zzzzzz_xxxxx_1[j];

                    tg_xzzzzzz_xxxxxy_0[j] = pb_x * tg_zzzzzz_xxxxxy_0[j] + fr * tg_zzzzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_zzzzzz_xxxxy_1[j];

                    tg_xzzzzzz_xxxxxz_0[j] = pb_x * tg_zzzzzz_xxxxxz_0[j] + fr * tg_zzzzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_zzzzzz_xxxxz_1[j];

                    tg_xzzzzzz_xxxxyy_0[j] = pb_x * tg_zzzzzz_xxxxyy_0[j] + fr * tg_zzzzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xxxyy_1[j];

                    tg_xzzzzzz_xxxxyz_0[j] = pb_x * tg_zzzzzz_xxxxyz_0[j] + fr * tg_zzzzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xxxyz_1[j];

                    tg_xzzzzzz_xxxxzz_0[j] = pb_x * tg_zzzzzz_xxxxzz_0[j] + fr * tg_zzzzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xxxzz_1[j];

                    tg_xzzzzzz_xxxyyy_0[j] = pb_x * tg_zzzzzz_xxxyyy_0[j] + fr * tg_zzzzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxyyy_1[j];

                    tg_xzzzzzz_xxxyyz_0[j] = pb_x * tg_zzzzzz_xxxyyz_0[j] + fr * tg_zzzzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxyyz_1[j];

                    tg_xzzzzzz_xxxyzz_0[j] = pb_x * tg_zzzzzz_xxxyzz_0[j] + fr * tg_zzzzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxyzz_1[j];

                    tg_xzzzzzz_xxxzzz_0[j] = pb_x * tg_zzzzzz_xxxzzz_0[j] + fr * tg_zzzzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxzzz_1[j];

                    tg_xzzzzzz_xxyyyy_0[j] = pb_x * tg_zzzzzz_xxyyyy_0[j] + fr * tg_zzzzzz_xxyyyy_1[j] + fl1_fxn * tg_zzzzzz_xyyyy_1[j];

                    tg_xzzzzzz_xxyyyz_0[j] = pb_x * tg_zzzzzz_xxyyyz_0[j] + fr * tg_zzzzzz_xxyyyz_1[j] + fl1_fxn * tg_zzzzzz_xyyyz_1[j];

                    tg_xzzzzzz_xxyyzz_0[j] = pb_x * tg_zzzzzz_xxyyzz_0[j] + fr * tg_zzzzzz_xxyyzz_1[j] + fl1_fxn * tg_zzzzzz_xyyzz_1[j];

                    tg_xzzzzzz_xxyzzz_0[j] = pb_x * tg_zzzzzz_xxyzzz_0[j] + fr * tg_zzzzzz_xxyzzz_1[j] + fl1_fxn * tg_zzzzzz_xyzzz_1[j];

                    tg_xzzzzzz_xxzzzz_0[j] = pb_x * tg_zzzzzz_xxzzzz_0[j] + fr * tg_zzzzzz_xxzzzz_1[j] + fl1_fxn * tg_zzzzzz_xzzzz_1[j];

                    tg_xzzzzzz_xyyyyy_0[j] = pb_x * tg_zzzzzz_xyyyyy_0[j] + fr * tg_zzzzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyyyy_1[j];

                    tg_xzzzzzz_xyyyyz_0[j] = pb_x * tg_zzzzzz_xyyyyz_0[j] + fr * tg_zzzzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyyyz_1[j];

                    tg_xzzzzzz_xyyyzz_0[j] = pb_x * tg_zzzzzz_xyyyzz_0[j] + fr * tg_zzzzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyyzz_1[j];

                    tg_xzzzzzz_xyyzzz_0[j] = pb_x * tg_zzzzzz_xyyzzz_0[j] + fr * tg_zzzzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyzzz_1[j];

                    tg_xzzzzzz_xyzzzz_0[j] = pb_x * tg_zzzzzz_xyzzzz_0[j] + fr * tg_zzzzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yzzzz_1[j];

                    tg_xzzzzzz_xzzzzz_0[j] = pb_x * tg_zzzzzz_xzzzzz_0[j] + fr * tg_zzzzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zzzzz_1[j];

                    tg_xzzzzzz_yyyyyy_0[j] = pb_x * tg_zzzzzz_yyyyyy_0[j] + fr * tg_zzzzzz_yyyyyy_1[j];

                    tg_xzzzzzz_yyyyyz_0[j] = pb_x * tg_zzzzzz_yyyyyz_0[j] + fr * tg_zzzzzz_yyyyyz_1[j];

                    tg_xzzzzzz_yyyyzz_0[j] = pb_x * tg_zzzzzz_yyyyzz_0[j] + fr * tg_zzzzzz_yyyyzz_1[j];

                    tg_xzzzzzz_yyyzzz_0[j] = pb_x * tg_zzzzzz_yyyzzz_0[j] + fr * tg_zzzzzz_yyyzzz_1[j];

                    tg_xzzzzzz_yyzzzz_0[j] = pb_x * tg_zzzzzz_yyzzzz_0[j] + fr * tg_zzzzzz_yyzzzz_1[j];

                    tg_xzzzzzz_yzzzzz_0[j] = pb_x * tg_zzzzzz_yzzzzz_0[j] + fr * tg_zzzzzz_yzzzzz_1[j];

                    tg_xzzzzzz_zzzzzz_0[j] = pb_x * tg_zzzzzz_zzzzzz_0[j] + fr * tg_zzzzzz_zzzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyyyy_xxxxxx_0[j] = pb_y * tg_yyyyyy_xxxxxx_0[j] + fr * tg_yyyyyy_xxxxxx_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxxx_0[j] - tg_yyyyy_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyyyy_xxxxxy_0[j] = pb_y * tg_yyyyyy_xxxxxy_0[j] + fr * tg_yyyyyy_xxxxxy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxxy_0[j] - tg_yyyyy_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxxxx_1[j];

                    tg_yyyyyyy_xxxxxz_0[j] = pb_y * tg_yyyyyy_xxxxxz_0[j] + fr * tg_yyyyyy_xxxxxz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxxz_0[j] - tg_yyyyy_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyyyy_xxxxyy_0[j] = pb_y * tg_yyyyyy_xxxxyy_0[j] + fr * tg_yyyyyy_xxxxyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxyy_0[j] - tg_yyyyy_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xxxxy_1[j];

                    tg_yyyyyyy_xxxxyz_0[j] = pb_y * tg_yyyyyy_xxxxyz_0[j] + fr * tg_yyyyyy_xxxxyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxyz_0[j] - tg_yyyyy_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxxxz_1[j];

                    tg_yyyyyyy_xxxxzz_0[j] = pb_y * tg_yyyyyy_xxxxzz_0[j] + fr * tg_yyyyyy_xxxxzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxxzz_0[j] - tg_yyyyy_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyyyy_xxxyyy_0[j] = pb_y * tg_yyyyyy_xxxyyy_0[j] + fr * tg_yyyyyy_xxxyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxyyy_0[j] - tg_yyyyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_xxxyy_1[j];

                    tg_yyyyyyy_xxxyyz_0[j] = pb_y * tg_yyyyyy_xxxyyz_0[j] + fr * tg_yyyyyy_xxxyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxyyz_0[j] - tg_yyyyy_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xxxyz_1[j];

                    tg_yyyyyyy_xxxyzz_0[j] = pb_y * tg_yyyyyy_xxxyzz_0[j] + fr * tg_yyyyyy_xxxyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxyzz_0[j] - tg_yyyyy_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxxzz_1[j];

                    tg_yyyyyyy_xxxzzz_0[j] = pb_y * tg_yyyyyy_xxxzzz_0[j] + fr * tg_yyyyyy_xxxzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxzzz_0[j] - tg_yyyyy_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyyyy_xxyyyy_0[j] = pb_y * tg_yyyyyy_xxyyyy_0[j] + fr * tg_yyyyyy_xxyyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyyyy_0[j] - tg_yyyyy_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyy_xxyyy_1[j];

                    tg_yyyyyyy_xxyyyz_0[j] = pb_y * tg_yyyyyy_xxyyyz_0[j] + fr * tg_yyyyyy_xxyyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyyyz_0[j] - tg_yyyyy_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_xxyyz_1[j];

                    tg_yyyyyyy_xxyyzz_0[j] = pb_y * tg_yyyyyy_xxyyzz_0[j] + fr * tg_yyyyyy_xxyyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyyzz_0[j] - tg_yyyyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xxyzz_1[j];

                    tg_yyyyyyy_xxyzzz_0[j] = pb_y * tg_yyyyyy_xxyzzz_0[j] + fr * tg_yyyyyy_xxyzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyzzz_0[j] - tg_yyyyy_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxzzz_1[j];

                    tg_yyyyyyy_xxzzzz_0[j] = pb_y * tg_yyyyyy_xxzzzz_0[j] + fr * tg_yyyyyy_xxzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxzzzz_0[j] - tg_yyyyy_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyyyy_xyyyyy_0[j] = pb_y * tg_yyyyyy_xyyyyy_0[j] + fr * tg_yyyyyy_xyyyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyyyy_0[j] - tg_yyyyy_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyy_xyyyy_1[j];

                    tg_yyyyyyy_xyyyyz_0[j] = pb_y * tg_yyyyyy_xyyyyz_0[j] + fr * tg_yyyyyy_xyyyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyyyz_0[j] - tg_yyyyy_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyy_xyyyz_1[j];

                    tg_yyyyyyy_xyyyzz_0[j] = pb_y * tg_yyyyyy_xyyyzz_0[j] + fr * tg_yyyyyy_xyyyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyyzz_0[j] - tg_yyyyy_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_xyyzz_1[j];

                    tg_yyyyyyy_xyyzzz_0[j] = pb_y * tg_yyyyyy_xyyzzz_0[j] + fr * tg_yyyyyy_xyyzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyzzz_0[j] - tg_yyyyy_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xyzzz_1[j];

                    tg_yyyyyyy_xyzzzz_0[j] = pb_y * tg_yyyyyy_xyzzzz_0[j] + fr * tg_yyyyyy_xyzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyzzzz_0[j] - tg_yyyyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xzzzz_1[j];

                    tg_yyyyyyy_xzzzzz_0[j] = pb_y * tg_yyyyyy_xzzzzz_0[j] + fr * tg_yyyyyy_xzzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xzzzzz_0[j] - tg_yyyyy_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyyyy_yyyyyy_0[j] = pb_y * tg_yyyyyy_yyyyyy_0[j] + fr * tg_yyyyyy_yyyyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyyyy_0[j] - tg_yyyyy_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyyyy_yyyyy_1[j];

                    tg_yyyyyyy_yyyyyz_0[j] = pb_y * tg_yyyyyy_yyyyyz_0[j] + fr * tg_yyyyyy_yyyyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyyyz_0[j] - tg_yyyyy_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyy_yyyyz_1[j];

                    tg_yyyyyyy_yyyyzz_0[j] = pb_y * tg_yyyyyy_yyyyzz_0[j] + fr * tg_yyyyyy_yyyyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyyzz_0[j] - tg_yyyyy_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyy_yyyzz_1[j];

                    tg_yyyyyyy_yyyzzz_0[j] = pb_y * tg_yyyyyy_yyyzzz_0[j] + fr * tg_yyyyyy_yyyzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyzzz_0[j] - tg_yyyyy_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_yyzzz_1[j];

                    tg_yyyyyyy_yyzzzz_0[j] = pb_y * tg_yyyyyy_yyzzzz_0[j] + fr * tg_yyyyyy_yyzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyzzzz_0[j] - tg_yyyyy_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_yzzzz_1[j];

                    tg_yyyyyyy_yzzzzz_0[j] = pb_y * tg_yyyyyy_yzzzzz_0[j] + fr * tg_yyyyyy_yzzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yzzzzz_0[j] - tg_yyyyy_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_zzzzz_1[j];

                    tg_yyyyyyy_zzzzzz_0[j] = pb_y * tg_yyyyyy_zzzzzz_0[j] + fr * tg_yyyyyy_zzzzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_zzzzzz_0[j] - tg_yyyyy_zzzzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxxxx_0[j] = pb_y * tg_yyyyyz_xxxxxx_0[j] + fr * tg_yyyyyz_xxxxxx_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxxx_0[j] - tg_yyyyz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxxxy_0[j] = pb_y * tg_yyyyyz_xxxxxy_0[j] + fr * tg_yyyyyz_xxxxxy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxxy_0[j] - tg_yyyyz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxxxx_1[j];

                    tg_yyyyyyz_xxxxxz_0[j] = pb_y * tg_yyyyyz_xxxxxz_0[j] + fr * tg_yyyyyz_xxxxxz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxxz_0[j] - tg_yyyyz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxxyy_0[j] = pb_y * tg_yyyyyz_xxxxyy_0[j] + fr * tg_yyyyyz_xxxxyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxyy_0[j] - tg_yyyyz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xxxxy_1[j];

                    tg_yyyyyyz_xxxxyz_0[j] = pb_y * tg_yyyyyz_xxxxyz_0[j] + fr * tg_yyyyyz_xxxxyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxyz_0[j] - tg_yyyyz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxxxz_1[j];

                    tg_yyyyyyz_xxxxzz_0[j] = pb_y * tg_yyyyyz_xxxxzz_0[j] + fr * tg_yyyyyz_xxxxzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxxzz_0[j] - tg_yyyyz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxyyy_0[j] = pb_y * tg_yyyyyz_xxxyyy_0[j] + fr * tg_yyyyyz_xxxyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxyyy_0[j] - tg_yyyyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_xxxyy_1[j];

                    tg_yyyyyyz_xxxyyz_0[j] = pb_y * tg_yyyyyz_xxxyyz_0[j] + fr * tg_yyyyyz_xxxyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxyyz_0[j] - tg_yyyyz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xxxyz_1[j];

                    tg_yyyyyyz_xxxyzz_0[j] = pb_y * tg_yyyyyz_xxxyzz_0[j] + fr * tg_yyyyyz_xxxyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxyzz_0[j] - tg_yyyyz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxxzz_1[j];

                    tg_yyyyyyz_xxxzzz_0[j] = pb_y * tg_yyyyyz_xxxzzz_0[j] + fr * tg_yyyyyz_xxxzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxzzz_0[j] - tg_yyyyz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxyyyy_0[j] = pb_y * tg_yyyyyz_xxyyyy_0[j] + fr * tg_yyyyyz_xxyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyyyy_0[j] - tg_yyyyz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyz_xxyyy_1[j];

                    tg_yyyyyyz_xxyyyz_0[j] = pb_y * tg_yyyyyz_xxyyyz_0[j] + fr * tg_yyyyyz_xxyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyyyz_0[j] - tg_yyyyz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_xxyyz_1[j];

                    tg_yyyyyyz_xxyyzz_0[j] = pb_y * tg_yyyyyz_xxyyzz_0[j] + fr * tg_yyyyyz_xxyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyyzz_0[j] - tg_yyyyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xxyzz_1[j];

                    tg_yyyyyyz_xxyzzz_0[j] = pb_y * tg_yyyyyz_xxyzzz_0[j] + fr * tg_yyyyyz_xxyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyzzz_0[j] - tg_yyyyz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_826_917(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (826,917)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {7, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

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

                auto tg_yyyyyz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 630); 

                auto tg_yyyyyz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 631); 

                auto tg_yyyyyz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 632); 

                auto tg_yyyyyz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 633); 

                auto tg_yyyyyz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 634); 

                auto tg_yyyyyz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 635); 

                auto tg_yyyyyz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 636); 

                auto tg_yyyyyz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 637); 

                auto tg_yyyyyz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 638); 

                auto tg_yyyyyz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 639); 

                auto tg_yyyyyz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 640); 

                auto tg_yyyyyz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 641); 

                auto tg_yyyyyz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 642); 

                auto tg_yyyyyz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 643); 

                auto tg_yyyyzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 644); 

                auto tg_yyyyzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 645); 

                auto tg_yyyyzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 646); 

                auto tg_yyyyzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 647); 

                auto tg_yyyyzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 648); 

                auto tg_yyyyzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 649); 

                auto tg_yyyyzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 650); 

                auto tg_yyyyzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 651); 

                auto tg_yyyyzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 652); 

                auto tg_yyyyzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 653); 

                auto tg_yyyyzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 654); 

                auto tg_yyyyzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 655); 

                auto tg_yyyyzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 656); 

                auto tg_yyyyzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 657); 

                auto tg_yyyyzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 658); 

                auto tg_yyyyzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 659); 

                auto tg_yyyyzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 660); 

                auto tg_yyyyzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 661); 

                auto tg_yyyyzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 662); 

                auto tg_yyyyzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 663); 

                auto tg_yyyyzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 664); 

                auto tg_yyyyzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 665); 

                auto tg_yyyyzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 666); 

                auto tg_yyyyzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 667); 

                auto tg_yyyyzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 668); 

                auto tg_yyyyzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 669); 

                auto tg_yyyyzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 670); 

                auto tg_yyyyzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 671); 

                auto tg_yyyzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 672); 

                auto tg_yyyzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 673); 

                auto tg_yyyzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 674); 

                auto tg_yyyzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 675); 

                auto tg_yyyzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 676); 

                auto tg_yyyzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 677); 

                auto tg_yyyzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 678); 

                auto tg_yyyzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 679); 

                auto tg_yyyzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 680); 

                auto tg_yyyzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 681); 

                auto tg_yyyzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 682); 

                auto tg_yyyzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 683); 

                auto tg_yyyzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 684); 

                auto tg_yyyzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 685); 

                auto tg_yyyzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 686); 

                auto tg_yyyzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 687); 

                auto tg_yyyzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 688); 

                auto tg_yyyzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 689); 

                auto tg_yyyzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 690); 

                auto tg_yyyzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 691); 

                auto tg_yyyzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 692); 

                auto tg_yyyzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 693); 

                auto tg_yyyzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 694); 

                auto tg_yyyzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 695); 

                auto tg_yyyzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 696); 

                auto tg_yyyzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 697); 

                auto tg_yyyzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 698); 

                auto tg_yyyzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 699); 

                auto tg_yyzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 700); 

                auto tg_yyzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 701); 

                auto tg_yyzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 702); 

                auto tg_yyzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 703); 

                auto tg_yyzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 704); 

                auto tg_yyzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 705); 

                auto tg_yyzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 706); 

                auto tg_yyzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 707); 

                auto tg_yyzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 708); 

                auto tg_yyzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 709); 

                auto tg_yyzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 710); 

                auto tg_yyzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 711); 

                auto tg_yyzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 712); 

                auto tg_yyzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 713); 

                auto tg_yyzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 714); 

                auto tg_yyzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 715); 

                auto tg_yyzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 716); 

                auto tg_yyzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 717); 

                auto tg_yyzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 718); 

                auto tg_yyzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 719); 

                auto tg_yyzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 720); 

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

                auto tg_yyzzzz_xyyyz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 536); 

                auto tg_yyzzzz_xyyzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 537); 

                auto tg_yyzzzz_xyzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 538); 

                auto tg_yyzzzz_xzzzz_1 = primBuffer[pidx_g_6_5_m1].data(588 * idx + 539); 

                // set up pointers to integrals

                auto tg_yyyyyyz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 826); 

                auto tg_yyyyyyz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 827); 

                auto tg_yyyyyyz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 828); 

                auto tg_yyyyyyz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 829); 

                auto tg_yyyyyyz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 830); 

                auto tg_yyyyyyz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 831); 

                auto tg_yyyyyyz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 832); 

                auto tg_yyyyyyz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 833); 

                auto tg_yyyyyyz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 834); 

                auto tg_yyyyyyz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 835); 

                auto tg_yyyyyyz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 836); 

                auto tg_yyyyyyz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 837); 

                auto tg_yyyyyyz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 838); 

                auto tg_yyyyyyz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 839); 

                auto tg_yyyyyzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 840); 

                auto tg_yyyyyzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 841); 

                auto tg_yyyyyzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 842); 

                auto tg_yyyyyzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 843); 

                auto tg_yyyyyzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 844); 

                auto tg_yyyyyzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 845); 

                auto tg_yyyyyzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 846); 

                auto tg_yyyyyzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 847); 

                auto tg_yyyyyzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 848); 

                auto tg_yyyyyzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 849); 

                auto tg_yyyyyzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 850); 

                auto tg_yyyyyzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 851); 

                auto tg_yyyyyzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 852); 

                auto tg_yyyyyzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 853); 

                auto tg_yyyyyzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 854); 

                auto tg_yyyyyzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 855); 

                auto tg_yyyyyzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 856); 

                auto tg_yyyyyzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 857); 

                auto tg_yyyyyzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 858); 

                auto tg_yyyyyzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 859); 

                auto tg_yyyyyzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 860); 

                auto tg_yyyyyzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 861); 

                auto tg_yyyyyzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 862); 

                auto tg_yyyyyzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 863); 

                auto tg_yyyyyzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 864); 

                auto tg_yyyyyzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 865); 

                auto tg_yyyyyzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 866); 

                auto tg_yyyyyzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 867); 

                auto tg_yyyyzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 868); 

                auto tg_yyyyzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 869); 

                auto tg_yyyyzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 870); 

                auto tg_yyyyzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 871); 

                auto tg_yyyyzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 872); 

                auto tg_yyyyzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 873); 

                auto tg_yyyyzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 874); 

                auto tg_yyyyzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 875); 

                auto tg_yyyyzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 876); 

                auto tg_yyyyzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 877); 

                auto tg_yyyyzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 878); 

                auto tg_yyyyzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 879); 

                auto tg_yyyyzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 880); 

                auto tg_yyyyzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 881); 

                auto tg_yyyyzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 882); 

                auto tg_yyyyzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 883); 

                auto tg_yyyyzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 884); 

                auto tg_yyyyzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 885); 

                auto tg_yyyyzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 886); 

                auto tg_yyyyzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 887); 

                auto tg_yyyyzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 888); 

                auto tg_yyyyzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 889); 

                auto tg_yyyyzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 890); 

                auto tg_yyyyzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 891); 

                auto tg_yyyyzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 892); 

                auto tg_yyyyzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 893); 

                auto tg_yyyyzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 894); 

                auto tg_yyyyzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 895); 

                auto tg_yyyzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 896); 

                auto tg_yyyzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 897); 

                auto tg_yyyzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 898); 

                auto tg_yyyzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 899); 

                auto tg_yyyzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 900); 

                auto tg_yyyzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 901); 

                auto tg_yyyzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 902); 

                auto tg_yyyzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 903); 

                auto tg_yyyzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 904); 

                auto tg_yyyzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 905); 

                auto tg_yyyzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 906); 

                auto tg_yyyzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 907); 

                auto tg_yyyzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 908); 

                auto tg_yyyzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 909); 

                auto tg_yyyzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 910); 

                auto tg_yyyzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 911); 

                auto tg_yyyzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 912); 

                auto tg_yyyzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 913); 

                auto tg_yyyzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 914); 

                auto tg_yyyzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 915); 

                auto tg_yyyzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 916); 

                // Batch of Integrals (826,917)

                #pragma omp simd aligned(fxn, fza, tg_yyyyyyz_xxzzzz_0, tg_yyyyyyz_xyyyyy_0, \
                                         tg_yyyyyyz_xyyyyz_0, tg_yyyyyyz_xyyyzz_0, tg_yyyyyyz_xyyzzz_0, tg_yyyyyyz_xyzzzz_0, \
                                         tg_yyyyyyz_xzzzzz_0, tg_yyyyyyz_yyyyyy_0, tg_yyyyyyz_yyyyyz_0, tg_yyyyyyz_yyyyzz_0, \
                                         tg_yyyyyyz_yyyzzz_0, tg_yyyyyyz_yyzzzz_0, tg_yyyyyyz_yzzzzz_0, tg_yyyyyyz_zzzzzz_0, \
                                         tg_yyyyyz_xxzzzz_0, tg_yyyyyz_xxzzzz_1, tg_yyyyyz_xyyyy_1, tg_yyyyyz_xyyyyy_0, \
                                         tg_yyyyyz_xyyyyy_1, tg_yyyyyz_xyyyyz_0, tg_yyyyyz_xyyyyz_1, tg_yyyyyz_xyyyz_1, \
                                         tg_yyyyyz_xyyyzz_0, tg_yyyyyz_xyyyzz_1, tg_yyyyyz_xyyzz_1, tg_yyyyyz_xyyzzz_0, \
                                         tg_yyyyyz_xyyzzz_1, tg_yyyyyz_xyzzz_1, tg_yyyyyz_xyzzzz_0, tg_yyyyyz_xyzzzz_1, \
                                         tg_yyyyyz_xzzzz_1, tg_yyyyyz_xzzzzz_0, tg_yyyyyz_xzzzzz_1, tg_yyyyyz_yyyyy_1, \
                                         tg_yyyyyz_yyyyyy_0, tg_yyyyyz_yyyyyy_1, tg_yyyyyz_yyyyyz_0, tg_yyyyyz_yyyyyz_1, \
                                         tg_yyyyyz_yyyyz_1, tg_yyyyyz_yyyyzz_0, tg_yyyyyz_yyyyzz_1, tg_yyyyyz_yyyzz_1, \
                                         tg_yyyyyz_yyyzzz_0, tg_yyyyyz_yyyzzz_1, tg_yyyyyz_yyzzz_1, tg_yyyyyz_yyzzzz_0, \
                                         tg_yyyyyz_yyzzzz_1, tg_yyyyyz_yzzzz_1, tg_yyyyyz_yzzzzz_0, tg_yyyyyz_yzzzzz_1, \
                                         tg_yyyyyz_zzzzz_1, tg_yyyyyz_zzzzzz_0, tg_yyyyyz_zzzzzz_1, tg_yyyyyzz_xxxxxx_0, \
                                         tg_yyyyyzz_xxxxxy_0, tg_yyyyyzz_xxxxxz_0, tg_yyyyyzz_xxxxyy_0, tg_yyyyyzz_xxxxyz_0, \
                                         tg_yyyyyzz_xxxxzz_0, tg_yyyyyzz_xxxyyy_0, tg_yyyyyzz_xxxyyz_0, tg_yyyyyzz_xxxyzz_0, \
                                         tg_yyyyyzz_xxxzzz_0, tg_yyyyyzz_xxyyyy_0, tg_yyyyyzz_xxyyyz_0, tg_yyyyyzz_xxyyzz_0, \
                                         tg_yyyyyzz_xxyzzz_0, tg_yyyyyzz_xxzzzz_0, tg_yyyyyzz_xyyyyy_0, tg_yyyyyzz_xyyyyz_0, \
                                         tg_yyyyyzz_xyyyzz_0, tg_yyyyyzz_xyyzzz_0, tg_yyyyyzz_xyzzzz_0, tg_yyyyyzz_xzzzzz_0, \
                                         tg_yyyyyzz_yyyyyy_0, tg_yyyyyzz_yyyyyz_0, tg_yyyyyzz_yyyyzz_0, tg_yyyyyzz_yyyzzz_0, \
                                         tg_yyyyyzz_yyzzzz_0, tg_yyyyyzz_yzzzzz_0, tg_yyyyyzz_zzzzzz_0, tg_yyyyz_xxzzzz_0, \
                                         tg_yyyyz_xxzzzz_1, tg_yyyyz_xyyyyy_0, tg_yyyyz_xyyyyy_1, tg_yyyyz_xyyyyz_0, \
                                         tg_yyyyz_xyyyyz_1, tg_yyyyz_xyyyzz_0, tg_yyyyz_xyyyzz_1, tg_yyyyz_xyyzzz_0, \
                                         tg_yyyyz_xyyzzz_1, tg_yyyyz_xyzzzz_0, tg_yyyyz_xyzzzz_1, tg_yyyyz_xzzzzz_0, \
                                         tg_yyyyz_xzzzzz_1, tg_yyyyz_yyyyyy_0, tg_yyyyz_yyyyyy_1, tg_yyyyz_yyyyyz_0, \
                                         tg_yyyyz_yyyyyz_1, tg_yyyyz_yyyyzz_0, tg_yyyyz_yyyyzz_1, tg_yyyyz_yyyzzz_0, \
                                         tg_yyyyz_yyyzzz_1, tg_yyyyz_yyzzzz_0, tg_yyyyz_yyzzzz_1, tg_yyyyz_yzzzzz_0, \
                                         tg_yyyyz_yzzzzz_1, tg_yyyyz_zzzzzz_0, tg_yyyyz_zzzzzz_1, tg_yyyyzz_xxxxx_1, \
                                         tg_yyyyzz_xxxxxx_0, tg_yyyyzz_xxxxxx_1, tg_yyyyzz_xxxxxy_0, tg_yyyyzz_xxxxxy_1, \
                                         tg_yyyyzz_xxxxxz_0, tg_yyyyzz_xxxxxz_1, tg_yyyyzz_xxxxy_1, tg_yyyyzz_xxxxyy_0, \
                                         tg_yyyyzz_xxxxyy_1, tg_yyyyzz_xxxxyz_0, tg_yyyyzz_xxxxyz_1, tg_yyyyzz_xxxxz_1, \
                                         tg_yyyyzz_xxxxzz_0, tg_yyyyzz_xxxxzz_1, tg_yyyyzz_xxxyy_1, tg_yyyyzz_xxxyyy_0, \
                                         tg_yyyyzz_xxxyyy_1, tg_yyyyzz_xxxyyz_0, tg_yyyyzz_xxxyyz_1, tg_yyyyzz_xxxyz_1, \
                                         tg_yyyyzz_xxxyzz_0, tg_yyyyzz_xxxyzz_1, tg_yyyyzz_xxxzz_1, tg_yyyyzz_xxxzzz_0, \
                                         tg_yyyyzz_xxxzzz_1, tg_yyyyzz_xxyyy_1, tg_yyyyzz_xxyyyy_0, tg_yyyyzz_xxyyyy_1, \
                                         tg_yyyyzz_xxyyyz_0, tg_yyyyzz_xxyyyz_1, tg_yyyyzz_xxyyz_1, tg_yyyyzz_xxyyzz_0, \
                                         tg_yyyyzz_xxyyzz_1, tg_yyyyzz_xxyzz_1, tg_yyyyzz_xxyzzz_0, tg_yyyyzz_xxyzzz_1, \
                                         tg_yyyyzz_xxzzz_1, tg_yyyyzz_xxzzzz_0, tg_yyyyzz_xxzzzz_1, tg_yyyyzz_xyyyy_1, \
                                         tg_yyyyzz_xyyyyy_0, tg_yyyyzz_xyyyyy_1, tg_yyyyzz_xyyyyz_0, tg_yyyyzz_xyyyyz_1, \
                                         tg_yyyyzz_xyyyz_1, tg_yyyyzz_xyyyzz_0, tg_yyyyzz_xyyyzz_1, tg_yyyyzz_xyyzz_1, \
                                         tg_yyyyzz_xyyzzz_0, tg_yyyyzz_xyyzzz_1, tg_yyyyzz_xyzzz_1, tg_yyyyzz_xyzzzz_0, \
                                         tg_yyyyzz_xyzzzz_1, tg_yyyyzz_xzzzz_1, tg_yyyyzz_xzzzzz_0, tg_yyyyzz_xzzzzz_1, \
                                         tg_yyyyzz_yyyyy_1, tg_yyyyzz_yyyyyy_0, tg_yyyyzz_yyyyyy_1, tg_yyyyzz_yyyyyz_0, \
                                         tg_yyyyzz_yyyyyz_1, tg_yyyyzz_yyyyz_1, tg_yyyyzz_yyyyzz_0, tg_yyyyzz_yyyyzz_1, \
                                         tg_yyyyzz_yyyzz_1, tg_yyyyzz_yyyzzz_0, tg_yyyyzz_yyyzzz_1, tg_yyyyzz_yyzzz_1, \
                                         tg_yyyyzz_yyzzzz_0, tg_yyyyzz_yyzzzz_1, tg_yyyyzz_yzzzz_1, tg_yyyyzz_yzzzzz_0, \
                                         tg_yyyyzz_yzzzzz_1, tg_yyyyzz_zzzzz_1, tg_yyyyzz_zzzzzz_0, tg_yyyyzz_zzzzzz_1, \
                                         tg_yyyyzzz_xxxxxx_0, tg_yyyyzzz_xxxxxy_0, tg_yyyyzzz_xxxxxz_0, tg_yyyyzzz_xxxxyy_0, \
                                         tg_yyyyzzz_xxxxyz_0, tg_yyyyzzz_xxxxzz_0, tg_yyyyzzz_xxxyyy_0, tg_yyyyzzz_xxxyyz_0, \
                                         tg_yyyyzzz_xxxyzz_0, tg_yyyyzzz_xxxzzz_0, tg_yyyyzzz_xxyyyy_0, tg_yyyyzzz_xxyyyz_0, \
                                         tg_yyyyzzz_xxyyzz_0, tg_yyyyzzz_xxyzzz_0, tg_yyyyzzz_xxzzzz_0, tg_yyyyzzz_xyyyyy_0, \
                                         tg_yyyyzzz_xyyyyz_0, tg_yyyyzzz_xyyyzz_0, tg_yyyyzzz_xyyzzz_0, tg_yyyyzzz_xyzzzz_0, \
                                         tg_yyyyzzz_xzzzzz_0, tg_yyyyzzz_yyyyyy_0, tg_yyyyzzz_yyyyyz_0, tg_yyyyzzz_yyyyzz_0, \
                                         tg_yyyyzzz_yyyzzz_0, tg_yyyyzzz_yyzzzz_0, tg_yyyyzzz_yzzzzz_0, tg_yyyyzzz_zzzzzz_0, \
                                         tg_yyyzz_xxxxxx_0, tg_yyyzz_xxxxxx_1, tg_yyyzz_xxxxxy_0, tg_yyyzz_xxxxxy_1, \
                                         tg_yyyzz_xxxxxz_0, tg_yyyzz_xxxxxz_1, tg_yyyzz_xxxxyy_0, tg_yyyzz_xxxxyy_1, \
                                         tg_yyyzz_xxxxyz_0, tg_yyyzz_xxxxyz_1, tg_yyyzz_xxxxzz_0, tg_yyyzz_xxxxzz_1, \
                                         tg_yyyzz_xxxyyy_0, tg_yyyzz_xxxyyy_1, tg_yyyzz_xxxyyz_0, tg_yyyzz_xxxyyz_1, \
                                         tg_yyyzz_xxxyzz_0, tg_yyyzz_xxxyzz_1, tg_yyyzz_xxxzzz_0, tg_yyyzz_xxxzzz_1, \
                                         tg_yyyzz_xxyyyy_0, tg_yyyzz_xxyyyy_1, tg_yyyzz_xxyyyz_0, tg_yyyzz_xxyyyz_1, \
                                         tg_yyyzz_xxyyzz_0, tg_yyyzz_xxyyzz_1, tg_yyyzz_xxyzzz_0, tg_yyyzz_xxyzzz_1, \
                                         tg_yyyzz_xxzzzz_0, tg_yyyzz_xxzzzz_1, tg_yyyzz_xyyyyy_0, tg_yyyzz_xyyyyy_1, \
                                         tg_yyyzz_xyyyyz_0, tg_yyyzz_xyyyyz_1, tg_yyyzz_xyyyzz_0, tg_yyyzz_xyyyzz_1, \
                                         tg_yyyzz_xyyzzz_0, tg_yyyzz_xyyzzz_1, tg_yyyzz_xyzzzz_0, tg_yyyzz_xyzzzz_1, \
                                         tg_yyyzz_xzzzzz_0, tg_yyyzz_xzzzzz_1, tg_yyyzz_yyyyyy_0, tg_yyyzz_yyyyyy_1, \
                                         tg_yyyzz_yyyyyz_0, tg_yyyzz_yyyyyz_1, tg_yyyzz_yyyyzz_0, tg_yyyzz_yyyyzz_1, \
                                         tg_yyyzz_yyyzzz_0, tg_yyyzz_yyyzzz_1, tg_yyyzz_yyzzzz_0, tg_yyyzz_yyzzzz_1, \
                                         tg_yyyzz_yzzzzz_0, tg_yyyzz_yzzzzz_1, tg_yyyzz_zzzzzz_0, tg_yyyzz_zzzzzz_1, \
                                         tg_yyyzzz_xxxxx_1, tg_yyyzzz_xxxxxx_0, tg_yyyzzz_xxxxxx_1, tg_yyyzzz_xxxxxy_0, \
                                         tg_yyyzzz_xxxxxy_1, tg_yyyzzz_xxxxxz_0, tg_yyyzzz_xxxxxz_1, tg_yyyzzz_xxxxy_1, \
                                         tg_yyyzzz_xxxxyy_0, tg_yyyzzz_xxxxyy_1, tg_yyyzzz_xxxxyz_0, tg_yyyzzz_xxxxyz_1, \
                                         tg_yyyzzz_xxxxz_1, tg_yyyzzz_xxxxzz_0, tg_yyyzzz_xxxxzz_1, tg_yyyzzz_xxxyy_1, \
                                         tg_yyyzzz_xxxyyy_0, tg_yyyzzz_xxxyyy_1, tg_yyyzzz_xxxyyz_0, tg_yyyzzz_xxxyyz_1, \
                                         tg_yyyzzz_xxxyz_1, tg_yyyzzz_xxxyzz_0, tg_yyyzzz_xxxyzz_1, tg_yyyzzz_xxxzz_1, \
                                         tg_yyyzzz_xxxzzz_0, tg_yyyzzz_xxxzzz_1, tg_yyyzzz_xxyyy_1, tg_yyyzzz_xxyyyy_0, \
                                         tg_yyyzzz_xxyyyy_1, tg_yyyzzz_xxyyyz_0, tg_yyyzzz_xxyyyz_1, tg_yyyzzz_xxyyz_1, \
                                         tg_yyyzzz_xxyyzz_0, tg_yyyzzz_xxyyzz_1, tg_yyyzzz_xxyzz_1, tg_yyyzzz_xxyzzz_0, \
                                         tg_yyyzzz_xxyzzz_1, tg_yyyzzz_xxzzz_1, tg_yyyzzz_xxzzzz_0, tg_yyyzzz_xxzzzz_1, \
                                         tg_yyyzzz_xyyyy_1, tg_yyyzzz_xyyyyy_0, tg_yyyzzz_xyyyyy_1, tg_yyyzzz_xyyyyz_0, \
                                         tg_yyyzzz_xyyyyz_1, tg_yyyzzz_xyyyz_1, tg_yyyzzz_xyyyzz_0, tg_yyyzzz_xyyyzz_1, \
                                         tg_yyyzzz_xyyzz_1, tg_yyyzzz_xyyzzz_0, tg_yyyzzz_xyyzzz_1, tg_yyyzzz_xyzzz_1, \
                                         tg_yyyzzz_xyzzzz_0, tg_yyyzzz_xyzzzz_1, tg_yyyzzz_xzzzz_1, tg_yyyzzz_xzzzzz_0, \
                                         tg_yyyzzz_xzzzzz_1, tg_yyyzzz_yyyyy_1, tg_yyyzzz_yyyyyy_0, tg_yyyzzz_yyyyyy_1, \
                                         tg_yyyzzz_yyyyyz_0, tg_yyyzzz_yyyyyz_1, tg_yyyzzz_yyyyz_1, tg_yyyzzz_yyyyzz_0, \
                                         tg_yyyzzz_yyyyzz_1, tg_yyyzzz_yyyzz_1, tg_yyyzzz_yyyzzz_0, tg_yyyzzz_yyyzzz_1, \
                                         tg_yyyzzz_yyzzz_1, tg_yyyzzz_yyzzzz_0, tg_yyyzzz_yyzzzz_1, tg_yyyzzz_yzzzz_1, \
                                         tg_yyyzzz_yzzzzz_0, tg_yyyzzz_yzzzzz_1, tg_yyyzzz_zzzzz_1, tg_yyyzzz_zzzzzz_0, \
                                         tg_yyyzzz_zzzzzz_1, tg_yyyzzzz_xxxxxx_0, tg_yyyzzzz_xxxxxy_0, tg_yyyzzzz_xxxxxz_0, \
                                         tg_yyyzzzz_xxxxyy_0, tg_yyyzzzz_xxxxyz_0, tg_yyyzzzz_xxxxzz_0, tg_yyyzzzz_xxxyyy_0, \
                                         tg_yyyzzzz_xxxyyz_0, tg_yyyzzzz_xxxyzz_0, tg_yyyzzzz_xxxzzz_0, tg_yyyzzzz_xxyyyy_0, \
                                         tg_yyyzzzz_xxyyyz_0, tg_yyyzzzz_xxyyzz_0, tg_yyyzzzz_xxyzzz_0, tg_yyyzzzz_xxzzzz_0, \
                                         tg_yyyzzzz_xyyyyy_0, tg_yyyzzzz_xyyyyz_0, tg_yyyzzzz_xyyyzz_0, tg_yyyzzzz_xyyzzz_0, \
                                         tg_yyyzzzz_xyzzzz_0, tg_yyyzzzz_xzzzzz_0, tg_yyzzz_xxxxxx_0, tg_yyzzz_xxxxxx_1, \
                                         tg_yyzzz_xxxxxy_0, tg_yyzzz_xxxxxy_1, tg_yyzzz_xxxxxz_0, tg_yyzzz_xxxxxz_1, \
                                         tg_yyzzz_xxxxyy_0, tg_yyzzz_xxxxyy_1, tg_yyzzz_xxxxyz_0, tg_yyzzz_xxxxyz_1, \
                                         tg_yyzzz_xxxxzz_0, tg_yyzzz_xxxxzz_1, tg_yyzzz_xxxyyy_0, tg_yyzzz_xxxyyy_1, \
                                         tg_yyzzz_xxxyyz_0, tg_yyzzz_xxxyyz_1, tg_yyzzz_xxxyzz_0, tg_yyzzz_xxxyzz_1, \
                                         tg_yyzzz_xxxzzz_0, tg_yyzzz_xxxzzz_1, tg_yyzzz_xxyyyy_0, tg_yyzzz_xxyyyy_1, \
                                         tg_yyzzz_xxyyyz_0, tg_yyzzz_xxyyyz_1, tg_yyzzz_xxyyzz_0, tg_yyzzz_xxyyzz_1, \
                                         tg_yyzzz_xxyzzz_0, tg_yyzzz_xxyzzz_1, tg_yyzzz_xxzzzz_0, tg_yyzzz_xxzzzz_1, \
                                         tg_yyzzz_xyyyyy_0, tg_yyzzz_xyyyyy_1, tg_yyzzz_xyyyyz_0, tg_yyzzz_xyyyyz_1, \
                                         tg_yyzzz_xyyyzz_0, tg_yyzzz_xyyyzz_1, tg_yyzzz_xyyzzz_0, tg_yyzzz_xyyzzz_1, \
                                         tg_yyzzz_xyzzzz_0, tg_yyzzz_xyzzzz_1, tg_yyzzz_xzzzzz_0, tg_yyzzz_xzzzzz_1, \
                                         tg_yyzzz_yyyyyy_0, tg_yyzzz_yyyyyy_1, tg_yyzzz_yyyyyz_0, tg_yyzzz_yyyyyz_1, \
                                         tg_yyzzz_yyyyzz_0, tg_yyzzz_yyyyzz_1, tg_yyzzz_yyyzzz_0, tg_yyzzz_yyyzzz_1, \
                                         tg_yyzzz_yyzzzz_0, tg_yyzzz_yyzzzz_1, tg_yyzzz_yzzzzz_0, tg_yyzzz_yzzzzz_1, \
                                         tg_yyzzz_zzzzzz_0, tg_yyzzz_zzzzzz_1, tg_yyzzzz_xxxxx_1, tg_yyzzzz_xxxxxx_0, \
                                         tg_yyzzzz_xxxxxx_1, tg_yyzzzz_xxxxxy_0, tg_yyzzzz_xxxxxy_1, tg_yyzzzz_xxxxxz_0, \
                                         tg_yyzzzz_xxxxxz_1, tg_yyzzzz_xxxxy_1, tg_yyzzzz_xxxxyy_0, tg_yyzzzz_xxxxyy_1, \
                                         tg_yyzzzz_xxxxyz_0, tg_yyzzzz_xxxxyz_1, tg_yyzzzz_xxxxz_1, tg_yyzzzz_xxxxzz_0, \
                                         tg_yyzzzz_xxxxzz_1, tg_yyzzzz_xxxyy_1, tg_yyzzzz_xxxyyy_0, tg_yyzzzz_xxxyyy_1, \
                                         tg_yyzzzz_xxxyyz_0, tg_yyzzzz_xxxyyz_1, tg_yyzzzz_xxxyz_1, tg_yyzzzz_xxxyzz_0, \
                                         tg_yyzzzz_xxxyzz_1, tg_yyzzzz_xxxzz_1, tg_yyzzzz_xxxzzz_0, tg_yyzzzz_xxxzzz_1, \
                                         tg_yyzzzz_xxyyy_1, tg_yyzzzz_xxyyyy_0, tg_yyzzzz_xxyyyy_1, tg_yyzzzz_xxyyyz_0, \
                                         tg_yyzzzz_xxyyyz_1, tg_yyzzzz_xxyyz_1, tg_yyzzzz_xxyyzz_0, tg_yyzzzz_xxyyzz_1, \
                                         tg_yyzzzz_xxyzz_1, tg_yyzzzz_xxyzzz_0, tg_yyzzzz_xxyzzz_1, tg_yyzzzz_xxzzz_1, \
                                         tg_yyzzzz_xxzzzz_0, tg_yyzzzz_xxzzzz_1, tg_yyzzzz_xyyyy_1, tg_yyzzzz_xyyyyy_0, \
                                         tg_yyzzzz_xyyyyy_1, tg_yyzzzz_xyyyyz_0, tg_yyzzzz_xyyyyz_1, tg_yyzzzz_xyyyz_1, \
                                         tg_yyzzzz_xyyyzz_0, tg_yyzzzz_xyyyzz_1, tg_yyzzzz_xyyzz_1, tg_yyzzzz_xyyzzz_0, \
                                         tg_yyzzzz_xyyzzz_1, tg_yyzzzz_xyzzz_1, tg_yyzzzz_xyzzzz_0, tg_yyzzzz_xyzzzz_1, \
                                         tg_yyzzzz_xzzzz_1, tg_yyzzzz_xzzzzz_0, tg_yyzzzz_xzzzzz_1, tg_yzzzz_xxxxxx_0, \
                                         tg_yzzzz_xxxxxx_1, tg_yzzzz_xxxxxy_0, tg_yzzzz_xxxxxy_1, tg_yzzzz_xxxxxz_0, \
                                         tg_yzzzz_xxxxxz_1, tg_yzzzz_xxxxyy_0, tg_yzzzz_xxxxyy_1, tg_yzzzz_xxxxyz_0, \
                                         tg_yzzzz_xxxxyz_1, tg_yzzzz_xxxxzz_0, tg_yzzzz_xxxxzz_1, tg_yzzzz_xxxyyy_0, \
                                         tg_yzzzz_xxxyyy_1, tg_yzzzz_xxxyyz_0, tg_yzzzz_xxxyyz_1, tg_yzzzz_xxxyzz_0, \
                                         tg_yzzzz_xxxyzz_1, tg_yzzzz_xxxzzz_0, tg_yzzzz_xxxzzz_1, tg_yzzzz_xxyyyy_0, \
                                         tg_yzzzz_xxyyyy_1, tg_yzzzz_xxyyyz_0, tg_yzzzz_xxyyyz_1, tg_yzzzz_xxyyzz_0, \
                                         tg_yzzzz_xxyyzz_1, tg_yzzzz_xxyzzz_0, tg_yzzzz_xxyzzz_1, tg_yzzzz_xxzzzz_0, \
                                         tg_yzzzz_xxzzzz_1, tg_yzzzz_xyyyyy_0, tg_yzzzz_xyyyyy_1, tg_yzzzz_xyyyyz_0, \
                                         tg_yzzzz_xyyyyz_1, tg_yzzzz_xyyyzz_0, tg_yzzzz_xyyyzz_1, tg_yzzzz_xyyzzz_0, \
                                         tg_yzzzz_xyyzzz_1, tg_yzzzz_xyzzzz_0, tg_yzzzz_xyzzzz_1, tg_yzzzz_xzzzzz_0, \
                                         tg_yzzzz_xzzzzz_1, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyyyz_xxzzzz_0[j] = pb_y * tg_yyyyyz_xxzzzz_0[j] + fr * tg_yyyyyz_xxzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxzzzz_0[j] - tg_yyyyz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xyyyyy_0[j] = pb_y * tg_yyyyyz_xyyyyy_0[j] + fr * tg_yyyyyz_xyyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyyyy_0[j] - tg_yyyyz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyz_xyyyy_1[j];

                    tg_yyyyyyz_xyyyyz_0[j] = pb_y * tg_yyyyyz_xyyyyz_0[j] + fr * tg_yyyyyz_xyyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyyyz_0[j] - tg_yyyyz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyz_xyyyz_1[j];

                    tg_yyyyyyz_xyyyzz_0[j] = pb_y * tg_yyyyyz_xyyyzz_0[j] + fr * tg_yyyyyz_xyyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyyzz_0[j] - tg_yyyyz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_xyyzz_1[j];

                    tg_yyyyyyz_xyyzzz_0[j] = pb_y * tg_yyyyyz_xyyzzz_0[j] + fr * tg_yyyyyz_xyyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyzzz_0[j] - tg_yyyyz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xyzzz_1[j];

                    tg_yyyyyyz_xyzzzz_0[j] = pb_y * tg_yyyyyz_xyzzzz_0[j] + fr * tg_yyyyyz_xyzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyzzzz_0[j] - tg_yyyyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xzzzz_1[j];

                    tg_yyyyyyz_xzzzzz_0[j] = pb_y * tg_yyyyyz_xzzzzz_0[j] + fr * tg_yyyyyz_xzzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xzzzzz_0[j] - tg_yyyyz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_yyyyyy_0[j] = pb_y * tg_yyyyyz_yyyyyy_0[j] + fr * tg_yyyyyz_yyyyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyyyy_0[j] - tg_yyyyz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyyyz_yyyyy_1[j];

                    tg_yyyyyyz_yyyyyz_0[j] = pb_y * tg_yyyyyz_yyyyyz_0[j] + fr * tg_yyyyyz_yyyyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyyyz_0[j] - tg_yyyyz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyyz_yyyyz_1[j];

                    tg_yyyyyyz_yyyyzz_0[j] = pb_y * tg_yyyyyz_yyyyzz_0[j] + fr * tg_yyyyyz_yyyyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyyzz_0[j] - tg_yyyyz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyz_yyyzz_1[j];

                    tg_yyyyyyz_yyyzzz_0[j] = pb_y * tg_yyyyyz_yyyzzz_0[j] + fr * tg_yyyyyz_yyyzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyzzz_0[j] - tg_yyyyz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_yyzzz_1[j];

                    tg_yyyyyyz_yyzzzz_0[j] = pb_y * tg_yyyyyz_yyzzzz_0[j] + fr * tg_yyyyyz_yyzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyzzzz_0[j] - tg_yyyyz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_yzzzz_1[j];

                    tg_yyyyyyz_yzzzzz_0[j] = pb_y * tg_yyyyyz_yzzzzz_0[j] + fr * tg_yyyyyz_yzzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yzzzzz_0[j] - tg_yyyyz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_zzzzz_1[j];

                    tg_yyyyyyz_zzzzzz_0[j] = pb_y * tg_yyyyyz_zzzzzz_0[j] + fr * tg_yyyyyz_zzzzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_zzzzzz_0[j] - tg_yyyyz_zzzzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxxxx_0[j] = pb_y * tg_yyyyzz_xxxxxx_0[j] + fr * tg_yyyyzz_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxxx_0[j] - tg_yyyzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxxxy_0[j] = pb_y * tg_yyyyzz_xxxxxy_0[j] + fr * tg_yyyyzz_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxxy_0[j] - tg_yyyzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxxxx_1[j];

                    tg_yyyyyzz_xxxxxz_0[j] = pb_y * tg_yyyyzz_xxxxxz_0[j] + fr * tg_yyyyzz_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxxz_0[j] - tg_yyyzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxxyy_0[j] = pb_y * tg_yyyyzz_xxxxyy_0[j] + fr * tg_yyyyzz_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxyy_0[j] - tg_yyyzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xxxxy_1[j];

                    tg_yyyyyzz_xxxxyz_0[j] = pb_y * tg_yyyyzz_xxxxyz_0[j] + fr * tg_yyyyzz_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxyz_0[j] - tg_yyyzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxxxz_1[j];

                    tg_yyyyyzz_xxxxzz_0[j] = pb_y * tg_yyyyzz_xxxxzz_0[j] + fr * tg_yyyyzz_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxxzz_0[j] - tg_yyyzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxyyy_0[j] = pb_y * tg_yyyyzz_xxxyyy_0[j] + fr * tg_yyyyzz_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxyyy_0[j] - tg_yyyzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_xxxyy_1[j];

                    tg_yyyyyzz_xxxyyz_0[j] = pb_y * tg_yyyyzz_xxxyyz_0[j] + fr * tg_yyyyzz_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxyyz_0[j] - tg_yyyzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xxxyz_1[j];

                    tg_yyyyyzz_xxxyzz_0[j] = pb_y * tg_yyyyzz_xxxyzz_0[j] + fr * tg_yyyyzz_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxyzz_0[j] - tg_yyyzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxxzz_1[j];

                    tg_yyyyyzz_xxxzzz_0[j] = pb_y * tg_yyyyzz_xxxzzz_0[j] + fr * tg_yyyyzz_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxzzz_0[j] - tg_yyyzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxyyyy_0[j] = pb_y * tg_yyyyzz_xxyyyy_0[j] + fr * tg_yyyyzz_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyyyy_0[j] - tg_yyyzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzz_xxyyy_1[j];

                    tg_yyyyyzz_xxyyyz_0[j] = pb_y * tg_yyyyzz_xxyyyz_0[j] + fr * tg_yyyyzz_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyyyz_0[j] - tg_yyyzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_xxyyz_1[j];

                    tg_yyyyyzz_xxyyzz_0[j] = pb_y * tg_yyyyzz_xxyyzz_0[j] + fr * tg_yyyyzz_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyyzz_0[j] - tg_yyyzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xxyzz_1[j];

                    tg_yyyyyzz_xxyzzz_0[j] = pb_y * tg_yyyyzz_xxyzzz_0[j] + fr * tg_yyyyzz_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyzzz_0[j] - tg_yyyzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxzzz_1[j];

                    tg_yyyyyzz_xxzzzz_0[j] = pb_y * tg_yyyyzz_xxzzzz_0[j] + fr * tg_yyyyzz_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxzzzz_0[j] - tg_yyyzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xyyyyy_0[j] = pb_y * tg_yyyyzz_xyyyyy_0[j] + fr * tg_yyyyzz_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyyyy_0[j] - tg_yyyzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyzz_xyyyy_1[j];

                    tg_yyyyyzz_xyyyyz_0[j] = pb_y * tg_yyyyzz_xyyyyz_0[j] + fr * tg_yyyyzz_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyyyz_0[j] - tg_yyyzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzz_xyyyz_1[j];

                    tg_yyyyyzz_xyyyzz_0[j] = pb_y * tg_yyyyzz_xyyyzz_0[j] + fr * tg_yyyyzz_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyyzz_0[j] - tg_yyyzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_xyyzz_1[j];

                    tg_yyyyyzz_xyyzzz_0[j] = pb_y * tg_yyyyzz_xyyzzz_0[j] + fr * tg_yyyyzz_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyzzz_0[j] - tg_yyyzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xyzzz_1[j];

                    tg_yyyyyzz_xyzzzz_0[j] = pb_y * tg_yyyyzz_xyzzzz_0[j] + fr * tg_yyyyzz_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyzzzz_0[j] - tg_yyyzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xzzzz_1[j];

                    tg_yyyyyzz_xzzzzz_0[j] = pb_y * tg_yyyyzz_xzzzzz_0[j] + fr * tg_yyyyzz_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xzzzzz_0[j] - tg_yyyzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_yyyyyy_0[j] = pb_y * tg_yyyyzz_yyyyyy_0[j] + fr * tg_yyyyzz_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyyyy_0[j] - tg_yyyzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyyzz_yyyyy_1[j];

                    tg_yyyyyzz_yyyyyz_0[j] = pb_y * tg_yyyyzz_yyyyyz_0[j] + fr * tg_yyyyzz_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyyyz_0[j] - tg_yyyzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyyzz_yyyyz_1[j];

                    tg_yyyyyzz_yyyyzz_0[j] = pb_y * tg_yyyyzz_yyyyzz_0[j] + fr * tg_yyyyzz_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyyzz_0[j] - tg_yyyzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzz_yyyzz_1[j];

                    tg_yyyyyzz_yyyzzz_0[j] = pb_y * tg_yyyyzz_yyyzzz_0[j] + fr * tg_yyyyzz_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyzzz_0[j] - tg_yyyzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_yyzzz_1[j];

                    tg_yyyyyzz_yyzzzz_0[j] = pb_y * tg_yyyyzz_yyzzzz_0[j] + fr * tg_yyyyzz_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyzzzz_0[j] - tg_yyyzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_yzzzz_1[j];

                    tg_yyyyyzz_yzzzzz_0[j] = pb_y * tg_yyyyzz_yzzzzz_0[j] + fr * tg_yyyyzz_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yzzzzz_0[j] - tg_yyyzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_zzzzz_1[j];

                    tg_yyyyyzz_zzzzzz_0[j] = pb_y * tg_yyyyzz_zzzzzz_0[j] + fr * tg_yyyyzz_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_zzzzzz_0[j] - tg_yyyzz_zzzzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxxxx_0[j] = pb_y * tg_yyyzzz_xxxxxx_0[j] + fr * tg_yyyzzz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxxx_0[j] - tg_yyzzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxxxy_0[j] = pb_y * tg_yyyzzz_xxxxxy_0[j] + fr * tg_yyyzzz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxxy_0[j] - tg_yyzzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxxxx_1[j];

                    tg_yyyyzzz_xxxxxz_0[j] = pb_y * tg_yyyzzz_xxxxxz_0[j] + fr * tg_yyyzzz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxxz_0[j] - tg_yyzzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxxyy_0[j] = pb_y * tg_yyyzzz_xxxxyy_0[j] + fr * tg_yyyzzz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxyy_0[j] - tg_yyzzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xxxxy_1[j];

                    tg_yyyyzzz_xxxxyz_0[j] = pb_y * tg_yyyzzz_xxxxyz_0[j] + fr * tg_yyyzzz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxyz_0[j] - tg_yyzzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxxxz_1[j];

                    tg_yyyyzzz_xxxxzz_0[j] = pb_y * tg_yyyzzz_xxxxzz_0[j] + fr * tg_yyyzzz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxxzz_0[j] - tg_yyzzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxyyy_0[j] = pb_y * tg_yyyzzz_xxxyyy_0[j] + fr * tg_yyyzzz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxyyy_0[j] - tg_yyzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_xxxyy_1[j];

                    tg_yyyyzzz_xxxyyz_0[j] = pb_y * tg_yyyzzz_xxxyyz_0[j] + fr * tg_yyyzzz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxyyz_0[j] - tg_yyzzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xxxyz_1[j];

                    tg_yyyyzzz_xxxyzz_0[j] = pb_y * tg_yyyzzz_xxxyzz_0[j] + fr * tg_yyyzzz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxyzz_0[j] - tg_yyzzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxxzz_1[j];

                    tg_yyyyzzz_xxxzzz_0[j] = pb_y * tg_yyyzzz_xxxzzz_0[j] + fr * tg_yyyzzz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxzzz_0[j] - tg_yyzzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxyyyy_0[j] = pb_y * tg_yyyzzz_xxyyyy_0[j] + fr * tg_yyyzzz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyyyy_0[j] - tg_yyzzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzz_xxyyy_1[j];

                    tg_yyyyzzz_xxyyyz_0[j] = pb_y * tg_yyyzzz_xxyyyz_0[j] + fr * tg_yyyzzz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyyyz_0[j] - tg_yyzzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_xxyyz_1[j];

                    tg_yyyyzzz_xxyyzz_0[j] = pb_y * tg_yyyzzz_xxyyzz_0[j] + fr * tg_yyyzzz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyyzz_0[j] - tg_yyzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xxyzz_1[j];

                    tg_yyyyzzz_xxyzzz_0[j] = pb_y * tg_yyyzzz_xxyzzz_0[j] + fr * tg_yyyzzz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyzzz_0[j] - tg_yyzzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxzzz_1[j];

                    tg_yyyyzzz_xxzzzz_0[j] = pb_y * tg_yyyzzz_xxzzzz_0[j] + fr * tg_yyyzzz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxzzzz_0[j] - tg_yyzzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xyyyyy_0[j] = pb_y * tg_yyyzzz_xyyyyy_0[j] + fr * tg_yyyzzz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyyyy_0[j] - tg_yyzzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyzzz_xyyyy_1[j];

                    tg_yyyyzzz_xyyyyz_0[j] = pb_y * tg_yyyzzz_xyyyyz_0[j] + fr * tg_yyyzzz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyyyz_0[j] - tg_yyzzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzz_xyyyz_1[j];

                    tg_yyyyzzz_xyyyzz_0[j] = pb_y * tg_yyyzzz_xyyyzz_0[j] + fr * tg_yyyzzz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyyzz_0[j] - tg_yyzzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_xyyzz_1[j];

                    tg_yyyyzzz_xyyzzz_0[j] = pb_y * tg_yyyzzz_xyyzzz_0[j] + fr * tg_yyyzzz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyzzz_0[j] - tg_yyzzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xyzzz_1[j];

                    tg_yyyyzzz_xyzzzz_0[j] = pb_y * tg_yyyzzz_xyzzzz_0[j] + fr * tg_yyyzzz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyzzzz_0[j] - tg_yyzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xzzzz_1[j];

                    tg_yyyyzzz_xzzzzz_0[j] = pb_y * tg_yyyzzz_xzzzzz_0[j] + fr * tg_yyyzzz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xzzzzz_0[j] - tg_yyzzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_yyyyyy_0[j] = pb_y * tg_yyyzzz_yyyyyy_0[j] + fr * tg_yyyzzz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyyyy_0[j] - tg_yyzzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyzzz_yyyyy_1[j];

                    tg_yyyyzzz_yyyyyz_0[j] = pb_y * tg_yyyzzz_yyyyyz_0[j] + fr * tg_yyyzzz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyyyz_0[j] - tg_yyzzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyzzz_yyyyz_1[j];

                    tg_yyyyzzz_yyyyzz_0[j] = pb_y * tg_yyyzzz_yyyyzz_0[j] + fr * tg_yyyzzz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyyzz_0[j] - tg_yyzzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzz_yyyzz_1[j];

                    tg_yyyyzzz_yyyzzz_0[j] = pb_y * tg_yyyzzz_yyyzzz_0[j] + fr * tg_yyyzzz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyzzz_0[j] - tg_yyzzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_yyzzz_1[j];

                    tg_yyyyzzz_yyzzzz_0[j] = pb_y * tg_yyyzzz_yyzzzz_0[j] + fr * tg_yyyzzz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyzzzz_0[j] - tg_yyzzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_yzzzz_1[j];

                    tg_yyyyzzz_yzzzzz_0[j] = pb_y * tg_yyyzzz_yzzzzz_0[j] + fr * tg_yyyzzz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yzzzzz_0[j] - tg_yyzzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_zzzzz_1[j];

                    tg_yyyyzzz_zzzzzz_0[j] = pb_y * tg_yyyzzz_zzzzzz_0[j] + fr * tg_yyyzzz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_zzzzzz_0[j] - tg_yyzzz_zzzzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxxxx_0[j] = pb_y * tg_yyzzzz_xxxxxx_0[j] + fr * tg_yyzzzz_xxxxxx_1[j] + fl1_fx * (tg_yzzzz_xxxxxx_0[j] - tg_yzzzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxxxy_0[j] = pb_y * tg_yyzzzz_xxxxxy_0[j] + fr * tg_yyzzzz_xxxxxy_1[j] + fl1_fx * (tg_yzzzz_xxxxxy_0[j] - tg_yzzzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxxxx_1[j];

                    tg_yyyzzzz_xxxxxz_0[j] = pb_y * tg_yyzzzz_xxxxxz_0[j] + fr * tg_yyzzzz_xxxxxz_1[j] + fl1_fx * (tg_yzzzz_xxxxxz_0[j] - tg_yzzzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxxyy_0[j] = pb_y * tg_yyzzzz_xxxxyy_0[j] + fr * tg_yyzzzz_xxxxyy_1[j] + fl1_fx * (tg_yzzzz_xxxxyy_0[j] - tg_yzzzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xxxxy_1[j];

                    tg_yyyzzzz_xxxxyz_0[j] = pb_y * tg_yyzzzz_xxxxyz_0[j] + fr * tg_yyzzzz_xxxxyz_1[j] + fl1_fx * (tg_yzzzz_xxxxyz_0[j] - tg_yzzzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxxxz_1[j];

                    tg_yyyzzzz_xxxxzz_0[j] = pb_y * tg_yyzzzz_xxxxzz_0[j] + fr * tg_yyzzzz_xxxxzz_1[j] + fl1_fx * (tg_yzzzz_xxxxzz_0[j] - tg_yzzzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxyyy_0[j] = pb_y * tg_yyzzzz_xxxyyy_0[j] + fr * tg_yyzzzz_xxxyyy_1[j] + fl1_fx * (tg_yzzzz_xxxyyy_0[j] - tg_yzzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_xxxyy_1[j];

                    tg_yyyzzzz_xxxyyz_0[j] = pb_y * tg_yyzzzz_xxxyyz_0[j] + fr * tg_yyzzzz_xxxyyz_1[j] + fl1_fx * (tg_yzzzz_xxxyyz_0[j] - tg_yzzzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xxxyz_1[j];

                    tg_yyyzzzz_xxxyzz_0[j] = pb_y * tg_yyzzzz_xxxyzz_0[j] + fr * tg_yyzzzz_xxxyzz_1[j] + fl1_fx * (tg_yzzzz_xxxyzz_0[j] - tg_yzzzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxxzz_1[j];

                    tg_yyyzzzz_xxxzzz_0[j] = pb_y * tg_yyzzzz_xxxzzz_0[j] + fr * tg_yyzzzz_xxxzzz_1[j] + fl1_fx * (tg_yzzzz_xxxzzz_0[j] - tg_yzzzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxyyyy_0[j] = pb_y * tg_yyzzzz_xxyyyy_0[j] + fr * tg_yyzzzz_xxyyyy_1[j] + fl1_fx * (tg_yzzzz_xxyyyy_0[j] - tg_yzzzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzz_xxyyy_1[j];

                    tg_yyyzzzz_xxyyyz_0[j] = pb_y * tg_yyzzzz_xxyyyz_0[j] + fr * tg_yyzzzz_xxyyyz_1[j] + fl1_fx * (tg_yzzzz_xxyyyz_0[j] - tg_yzzzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_xxyyz_1[j];

                    tg_yyyzzzz_xxyyzz_0[j] = pb_y * tg_yyzzzz_xxyyzz_0[j] + fr * tg_yyzzzz_xxyyzz_1[j] + fl1_fx * (tg_yzzzz_xxyyzz_0[j] - tg_yzzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xxyzz_1[j];

                    tg_yyyzzzz_xxyzzz_0[j] = pb_y * tg_yyzzzz_xxyzzz_0[j] + fr * tg_yyzzzz_xxyzzz_1[j] + fl1_fx * (tg_yzzzz_xxyzzz_0[j] - tg_yzzzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxzzz_1[j];

                    tg_yyyzzzz_xxzzzz_0[j] = pb_y * tg_yyzzzz_xxzzzz_0[j] + fr * tg_yyzzzz_xxzzzz_1[j] + fl1_fx * (tg_yzzzz_xxzzzz_0[j] - tg_yzzzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xyyyyy_0[j] = pb_y * tg_yyzzzz_xyyyyy_0[j] + fr * tg_yyzzzz_xyyyyy_1[j] + fl1_fx * (tg_yzzzz_xyyyyy_0[j] - tg_yzzzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzzzz_xyyyy_1[j];

                    tg_yyyzzzz_xyyyyz_0[j] = pb_y * tg_yyzzzz_xyyyyz_0[j] + fr * tg_yyzzzz_xyyyyz_1[j] + fl1_fx * (tg_yzzzz_xyyyyz_0[j] - tg_yzzzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzz_xyyyz_1[j];

                    tg_yyyzzzz_xyyyzz_0[j] = pb_y * tg_yyzzzz_xyyyzz_0[j] + fr * tg_yyzzzz_xyyyzz_1[j] + fl1_fx * (tg_yzzzz_xyyyzz_0[j] - tg_yzzzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_xyyzz_1[j];

                    tg_yyyzzzz_xyyzzz_0[j] = pb_y * tg_yyzzzz_xyyzzz_0[j] + fr * tg_yyzzzz_xyyzzz_1[j] + fl1_fx * (tg_yzzzz_xyyzzz_0[j] - tg_yzzzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xyzzz_1[j];

                    tg_yyyzzzz_xyzzzz_0[j] = pb_y * tg_yyzzzz_xyzzzz_0[j] + fr * tg_yyzzzz_xyzzzz_1[j] + fl1_fx * (tg_yzzzz_xyzzzz_0[j] - tg_yzzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xzzzz_1[j];

                    tg_yyyzzzz_xzzzzz_0[j] = pb_y * tg_yyzzzz_xzzzzz_0[j] + fr * tg_yyzzzz_xzzzzz_1[j] + fl1_fx * (tg_yzzzz_xzzzzz_0[j] - tg_yzzzz_xzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSI_917_1008(      CMemBlock2D<double>* primBuffer,
                                          const CRecursionMap&       recursionMap,
                                          const CMemBlock2D<double>& osFactors,
                                          const CMemBlock2D<double>& wpDistances,
                                          const CGtoPairsBlock&      braGtoPairsBlock,
                                          const CGtoPairsBlock&      ketGtoPairsBlock,
                                          const int32_t              nKetPrimPairs,
                                          const int32_t              iContrPair)
    {
        // Batch of Integrals (917,1008)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

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

                auto tg_yyzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 721); 

                auto tg_yyzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 722); 

                auto tg_yyzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 723); 

                auto tg_yyzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 724); 

                auto tg_yyzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 725); 

                auto tg_yyzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 726); 

                auto tg_yyzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 727); 

                auto tg_yzzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 728); 

                auto tg_yzzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 729); 

                auto tg_yzzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 730); 

                auto tg_yzzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 731); 

                auto tg_yzzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 732); 

                auto tg_yzzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 733); 

                auto tg_yzzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 734); 

                auto tg_yzzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 735); 

                auto tg_yzzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 736); 

                auto tg_yzzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 737); 

                auto tg_yzzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 738); 

                auto tg_yzzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 739); 

                auto tg_yzzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 740); 

                auto tg_yzzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 741); 

                auto tg_yzzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 742); 

                auto tg_yzzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 743); 

                auto tg_yzzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 744); 

                auto tg_yzzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 745); 

                auto tg_yzzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 746); 

                auto tg_yzzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 747); 

                auto tg_yzzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 748); 

                auto tg_yzzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 749); 

                auto tg_yzzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 750); 

                auto tg_yzzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 751); 

                auto tg_yzzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 752); 

                auto tg_yzzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 753); 

                auto tg_yzzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 754); 

                auto tg_yzzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 755); 

                auto tg_zzzzzz_xxxxxx_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 756); 

                auto tg_zzzzzz_xxxxxy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 757); 

                auto tg_zzzzzz_xxxxxz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 758); 

                auto tg_zzzzzz_xxxxyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 759); 

                auto tg_zzzzzz_xxxxyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 760); 

                auto tg_zzzzzz_xxxxzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 761); 

                auto tg_zzzzzz_xxxyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 762); 

                auto tg_zzzzzz_xxxyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 763); 

                auto tg_zzzzzz_xxxyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 764); 

                auto tg_zzzzzz_xxxzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 765); 

                auto tg_zzzzzz_xxyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 766); 

                auto tg_zzzzzz_xxyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 767); 

                auto tg_zzzzzz_xxyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 768); 

                auto tg_zzzzzz_xxyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 769); 

                auto tg_zzzzzz_xxzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 770); 

                auto tg_zzzzzz_xyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 771); 

                auto tg_zzzzzz_xyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 772); 

                auto tg_zzzzzz_xyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 773); 

                auto tg_zzzzzz_xyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 774); 

                auto tg_zzzzzz_xyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 775); 

                auto tg_zzzzzz_xzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 776); 

                auto tg_zzzzzz_yyyyyy_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 777); 

                auto tg_zzzzzz_yyyyyz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 778); 

                auto tg_zzzzzz_yyyyzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 779); 

                auto tg_zzzzzz_yyyzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 780); 

                auto tg_zzzzzz_yyzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 781); 

                auto tg_zzzzzz_yzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 782); 

                auto tg_zzzzzz_zzzzzz_1 = primBuffer[pidx_g_6_6_m1].data(784 * idx + 783); 

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

                // set up pointers to integrals

                auto tg_yyyzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 917); 

                auto tg_yyyzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 918); 

                auto tg_yyyzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 919); 

                auto tg_yyyzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 920); 

                auto tg_yyyzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 921); 

                auto tg_yyyzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 922); 

                auto tg_yyyzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 923); 

                auto tg_yyzzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 924); 

                auto tg_yyzzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 925); 

                auto tg_yyzzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 926); 

                auto tg_yyzzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 927); 

                auto tg_yyzzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 928); 

                auto tg_yyzzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 929); 

                auto tg_yyzzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 930); 

                auto tg_yyzzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 931); 

                auto tg_yyzzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 932); 

                auto tg_yyzzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 933); 

                auto tg_yyzzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 934); 

                auto tg_yyzzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 935); 

                auto tg_yyzzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 936); 

                auto tg_yyzzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 937); 

                auto tg_yyzzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 938); 

                auto tg_yyzzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 939); 

                auto tg_yyzzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 940); 

                auto tg_yyzzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 941); 

                auto tg_yyzzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 942); 

                auto tg_yyzzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 943); 

                auto tg_yyzzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 944); 

                auto tg_yyzzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 945); 

                auto tg_yyzzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 946); 

                auto tg_yyzzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 947); 

                auto tg_yyzzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 948); 

                auto tg_yyzzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 949); 

                auto tg_yyzzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 950); 

                auto tg_yyzzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 951); 

                auto tg_yzzzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 952); 

                auto tg_yzzzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 953); 

                auto tg_yzzzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 954); 

                auto tg_yzzzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 955); 

                auto tg_yzzzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 956); 

                auto tg_yzzzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 957); 

                auto tg_yzzzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 958); 

                auto tg_yzzzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 959); 

                auto tg_yzzzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 960); 

                auto tg_yzzzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 961); 

                auto tg_yzzzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 962); 

                auto tg_yzzzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 963); 

                auto tg_yzzzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 964); 

                auto tg_yzzzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 965); 

                auto tg_yzzzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 966); 

                auto tg_yzzzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 967); 

                auto tg_yzzzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 968); 

                auto tg_yzzzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 969); 

                auto tg_yzzzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 970); 

                auto tg_yzzzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 971); 

                auto tg_yzzzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 972); 

                auto tg_yzzzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 973); 

                auto tg_yzzzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 974); 

                auto tg_yzzzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 975); 

                auto tg_yzzzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 976); 

                auto tg_yzzzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 977); 

                auto tg_yzzzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 978); 

                auto tg_yzzzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 979); 

                auto tg_zzzzzzz_xxxxxx_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 980); 

                auto tg_zzzzzzz_xxxxxy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 981); 

                auto tg_zzzzzzz_xxxxxz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 982); 

                auto tg_zzzzzzz_xxxxyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 983); 

                auto tg_zzzzzzz_xxxxyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 984); 

                auto tg_zzzzzzz_xxxxzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 985); 

                auto tg_zzzzzzz_xxxyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 986); 

                auto tg_zzzzzzz_xxxyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 987); 

                auto tg_zzzzzzz_xxxyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 988); 

                auto tg_zzzzzzz_xxxzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 989); 

                auto tg_zzzzzzz_xxyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 990); 

                auto tg_zzzzzzz_xxyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 991); 

                auto tg_zzzzzzz_xxyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 992); 

                auto tg_zzzzzzz_xxyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 993); 

                auto tg_zzzzzzz_xxzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 994); 

                auto tg_zzzzzzz_xyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 995); 

                auto tg_zzzzzzz_xyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 996); 

                auto tg_zzzzzzz_xyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 997); 

                auto tg_zzzzzzz_xyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 998); 

                auto tg_zzzzzzz_xyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 999); 

                auto tg_zzzzzzz_xzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1000); 

                auto tg_zzzzzzz_yyyyyy_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1001); 

                auto tg_zzzzzzz_yyyyyz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1002); 

                auto tg_zzzzzzz_yyyyzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1003); 

                auto tg_zzzzzzz_yyyzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1004); 

                auto tg_zzzzzzz_yyzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1005); 

                auto tg_zzzzzzz_yzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1006); 

                auto tg_zzzzzzz_zzzzzz_0 = primBuffer[pidx_g_7_6_m0].data(1008 * idx + 1007); 

                // Batch of Integrals (917,1008)

                #pragma omp simd aligned(fxn, fza, tg_yyyzzzz_yyyyyy_0, tg_yyyzzzz_yyyyyz_0, \
                                         tg_yyyzzzz_yyyyzz_0, tg_yyyzzzz_yyyzzz_0, tg_yyyzzzz_yyzzzz_0, tg_yyyzzzz_yzzzzz_0, \
                                         tg_yyyzzzz_zzzzzz_0, tg_yyzzzz_yyyyy_1, tg_yyzzzz_yyyyyy_0, tg_yyzzzz_yyyyyy_1, \
                                         tg_yyzzzz_yyyyyz_0, tg_yyzzzz_yyyyyz_1, tg_yyzzzz_yyyyz_1, tg_yyzzzz_yyyyzz_0, \
                                         tg_yyzzzz_yyyyzz_1, tg_yyzzzz_yyyzz_1, tg_yyzzzz_yyyzzz_0, tg_yyzzzz_yyyzzz_1, \
                                         tg_yyzzzz_yyzzz_1, tg_yyzzzz_yyzzzz_0, tg_yyzzzz_yyzzzz_1, tg_yyzzzz_yzzzz_1, \
                                         tg_yyzzzz_yzzzzz_0, tg_yyzzzz_yzzzzz_1, tg_yyzzzz_zzzzz_1, tg_yyzzzz_zzzzzz_0, \
                                         tg_yyzzzz_zzzzzz_1, tg_yyzzzzz_xxxxxx_0, tg_yyzzzzz_xxxxxy_0, tg_yyzzzzz_xxxxxz_0, \
                                         tg_yyzzzzz_xxxxyy_0, tg_yyzzzzz_xxxxyz_0, tg_yyzzzzz_xxxxzz_0, tg_yyzzzzz_xxxyyy_0, \
                                         tg_yyzzzzz_xxxyyz_0, tg_yyzzzzz_xxxyzz_0, tg_yyzzzzz_xxxzzz_0, tg_yyzzzzz_xxyyyy_0, \
                                         tg_yyzzzzz_xxyyyz_0, tg_yyzzzzz_xxyyzz_0, tg_yyzzzzz_xxyzzz_0, tg_yyzzzzz_xxzzzz_0, \
                                         tg_yyzzzzz_xyyyyy_0, tg_yyzzzzz_xyyyyz_0, tg_yyzzzzz_xyyyzz_0, tg_yyzzzzz_xyyzzz_0, \
                                         tg_yyzzzzz_xyzzzz_0, tg_yyzzzzz_xzzzzz_0, tg_yyzzzzz_yyyyyy_0, tg_yyzzzzz_yyyyyz_0, \
                                         tg_yyzzzzz_yyyyzz_0, tg_yyzzzzz_yyyzzz_0, tg_yyzzzzz_yyzzzz_0, tg_yyzzzzz_yzzzzz_0, \
                                         tg_yyzzzzz_zzzzzz_0, tg_yzzzz_yyyyyy_0, tg_yzzzz_yyyyyy_1, tg_yzzzz_yyyyyz_0, \
                                         tg_yzzzz_yyyyyz_1, tg_yzzzz_yyyyzz_0, tg_yzzzz_yyyyzz_1, tg_yzzzz_yyyzzz_0, \
                                         tg_yzzzz_yyyzzz_1, tg_yzzzz_yyzzzz_0, tg_yzzzz_yyzzzz_1, tg_yzzzz_yzzzzz_0, \
                                         tg_yzzzz_yzzzzz_1, tg_yzzzz_zzzzzz_0, tg_yzzzz_zzzzzz_1, tg_yzzzzz_xxxxx_1, \
                                         tg_yzzzzz_xxxxxx_0, tg_yzzzzz_xxxxxx_1, tg_yzzzzz_xxxxxy_0, tg_yzzzzz_xxxxxy_1, \
                                         tg_yzzzzz_xxxxxz_0, tg_yzzzzz_xxxxxz_1, tg_yzzzzz_xxxxy_1, tg_yzzzzz_xxxxyy_0, \
                                         tg_yzzzzz_xxxxyy_1, tg_yzzzzz_xxxxyz_0, tg_yzzzzz_xxxxyz_1, tg_yzzzzz_xxxxz_1, \
                                         tg_yzzzzz_xxxxzz_0, tg_yzzzzz_xxxxzz_1, tg_yzzzzz_xxxyy_1, tg_yzzzzz_xxxyyy_0, \
                                         tg_yzzzzz_xxxyyy_1, tg_yzzzzz_xxxyyz_0, tg_yzzzzz_xxxyyz_1, tg_yzzzzz_xxxyz_1, \
                                         tg_yzzzzz_xxxyzz_0, tg_yzzzzz_xxxyzz_1, tg_yzzzzz_xxxzz_1, tg_yzzzzz_xxxzzz_0, \
                                         tg_yzzzzz_xxxzzz_1, tg_yzzzzz_xxyyy_1, tg_yzzzzz_xxyyyy_0, tg_yzzzzz_xxyyyy_1, \
                                         tg_yzzzzz_xxyyyz_0, tg_yzzzzz_xxyyyz_1, tg_yzzzzz_xxyyz_1, tg_yzzzzz_xxyyzz_0, \
                                         tg_yzzzzz_xxyyzz_1, tg_yzzzzz_xxyzz_1, tg_yzzzzz_xxyzzz_0, tg_yzzzzz_xxyzzz_1, \
                                         tg_yzzzzz_xxzzz_1, tg_yzzzzz_xxzzzz_0, tg_yzzzzz_xxzzzz_1, tg_yzzzzz_xyyyy_1, \
                                         tg_yzzzzz_xyyyyy_0, tg_yzzzzz_xyyyyy_1, tg_yzzzzz_xyyyyz_0, tg_yzzzzz_xyyyyz_1, \
                                         tg_yzzzzz_xyyyz_1, tg_yzzzzz_xyyyzz_0, tg_yzzzzz_xyyyzz_1, tg_yzzzzz_xyyzz_1, \
                                         tg_yzzzzz_xyyzzz_0, tg_yzzzzz_xyyzzz_1, tg_yzzzzz_xyzzz_1, tg_yzzzzz_xyzzzz_0, \
                                         tg_yzzzzz_xyzzzz_1, tg_yzzzzz_xzzzz_1, tg_yzzzzz_xzzzzz_0, tg_yzzzzz_xzzzzz_1, \
                                         tg_yzzzzz_yyyyy_1, tg_yzzzzz_yyyyyy_0, tg_yzzzzz_yyyyyy_1, tg_yzzzzz_yyyyyz_0, \
                                         tg_yzzzzz_yyyyyz_1, tg_yzzzzz_yyyyz_1, tg_yzzzzz_yyyyzz_0, tg_yzzzzz_yyyyzz_1, \
                                         tg_yzzzzz_yyyzz_1, tg_yzzzzz_yyyzzz_0, tg_yzzzzz_yyyzzz_1, tg_yzzzzz_yyzzz_1, \
                                         tg_yzzzzz_yyzzzz_0, tg_yzzzzz_yyzzzz_1, tg_yzzzzz_yzzzz_1, tg_yzzzzz_yzzzzz_0, \
                                         tg_yzzzzz_yzzzzz_1, tg_yzzzzz_zzzzz_1, tg_yzzzzz_zzzzzz_0, tg_yzzzzz_zzzzzz_1, \
                                         tg_yzzzzzz_xxxxxx_0, tg_yzzzzzz_xxxxxy_0, tg_yzzzzzz_xxxxxz_0, tg_yzzzzzz_xxxxyy_0, \
                                         tg_yzzzzzz_xxxxyz_0, tg_yzzzzzz_xxxxzz_0, tg_yzzzzzz_xxxyyy_0, tg_yzzzzzz_xxxyyz_0, \
                                         tg_yzzzzzz_xxxyzz_0, tg_yzzzzzz_xxxzzz_0, tg_yzzzzzz_xxyyyy_0, tg_yzzzzzz_xxyyyz_0, \
                                         tg_yzzzzzz_xxyyzz_0, tg_yzzzzzz_xxyzzz_0, tg_yzzzzzz_xxzzzz_0, tg_yzzzzzz_xyyyyy_0, \
                                         tg_yzzzzzz_xyyyyz_0, tg_yzzzzzz_xyyyzz_0, tg_yzzzzzz_xyyzzz_0, tg_yzzzzzz_xyzzzz_0, \
                                         tg_yzzzzzz_xzzzzz_0, tg_yzzzzzz_yyyyyy_0, tg_yzzzzzz_yyyyyz_0, tg_yzzzzzz_yyyyzz_0, \
                                         tg_yzzzzzz_yyyzzz_0, tg_yzzzzzz_yyzzzz_0, tg_yzzzzzz_yzzzzz_0, tg_yzzzzzz_zzzzzz_0, \
                                         tg_zzzzz_xxxxxx_0, tg_zzzzz_xxxxxx_1, tg_zzzzz_xxxxxy_0, tg_zzzzz_xxxxxy_1, \
                                         tg_zzzzz_xxxxxz_0, tg_zzzzz_xxxxxz_1, tg_zzzzz_xxxxyy_0, tg_zzzzz_xxxxyy_1, \
                                         tg_zzzzz_xxxxyz_0, tg_zzzzz_xxxxyz_1, tg_zzzzz_xxxxzz_0, tg_zzzzz_xxxxzz_1, \
                                         tg_zzzzz_xxxyyy_0, tg_zzzzz_xxxyyy_1, tg_zzzzz_xxxyyz_0, tg_zzzzz_xxxyyz_1, \
                                         tg_zzzzz_xxxyzz_0, tg_zzzzz_xxxyzz_1, tg_zzzzz_xxxzzz_0, tg_zzzzz_xxxzzz_1, \
                                         tg_zzzzz_xxyyyy_0, tg_zzzzz_xxyyyy_1, tg_zzzzz_xxyyyz_0, tg_zzzzz_xxyyyz_1, \
                                         tg_zzzzz_xxyyzz_0, tg_zzzzz_xxyyzz_1, tg_zzzzz_xxyzzz_0, tg_zzzzz_xxyzzz_1, \
                                         tg_zzzzz_xxzzzz_0, tg_zzzzz_xxzzzz_1, tg_zzzzz_xyyyyy_0, tg_zzzzz_xyyyyy_1, \
                                         tg_zzzzz_xyyyyz_0, tg_zzzzz_xyyyyz_1, tg_zzzzz_xyyyzz_0, tg_zzzzz_xyyyzz_1, \
                                         tg_zzzzz_xyyzzz_0, tg_zzzzz_xyyzzz_1, tg_zzzzz_xyzzzz_0, tg_zzzzz_xyzzzz_1, \
                                         tg_zzzzz_xzzzzz_0, tg_zzzzz_xzzzzz_1, tg_zzzzz_yyyyyy_0, tg_zzzzz_yyyyyy_1, \
                                         tg_zzzzz_yyyyyz_0, tg_zzzzz_yyyyyz_1, tg_zzzzz_yyyyzz_0, tg_zzzzz_yyyyzz_1, \
                                         tg_zzzzz_yyyzzz_0, tg_zzzzz_yyyzzz_1, tg_zzzzz_yyzzzz_0, tg_zzzzz_yyzzzz_1, \
                                         tg_zzzzz_yzzzzz_0, tg_zzzzz_yzzzzz_1, tg_zzzzz_zzzzzz_0, tg_zzzzz_zzzzzz_1, \
                                         tg_zzzzzz_xxxxx_1, tg_zzzzzz_xxxxxx_0, tg_zzzzzz_xxxxxx_1, tg_zzzzzz_xxxxxy_0, \
                                         tg_zzzzzz_xxxxxy_1, tg_zzzzzz_xxxxxz_0, tg_zzzzzz_xxxxxz_1, tg_zzzzzz_xxxxy_1, \
                                         tg_zzzzzz_xxxxyy_0, tg_zzzzzz_xxxxyy_1, tg_zzzzzz_xxxxyz_0, tg_zzzzzz_xxxxyz_1, \
                                         tg_zzzzzz_xxxxz_1, tg_zzzzzz_xxxxzz_0, tg_zzzzzz_xxxxzz_1, tg_zzzzzz_xxxyy_1, \
                                         tg_zzzzzz_xxxyyy_0, tg_zzzzzz_xxxyyy_1, tg_zzzzzz_xxxyyz_0, tg_zzzzzz_xxxyyz_1, \
                                         tg_zzzzzz_xxxyz_1, tg_zzzzzz_xxxyzz_0, tg_zzzzzz_xxxyzz_1, tg_zzzzzz_xxxzz_1, \
                                         tg_zzzzzz_xxxzzz_0, tg_zzzzzz_xxxzzz_1, tg_zzzzzz_xxyyy_1, tg_zzzzzz_xxyyyy_0, \
                                         tg_zzzzzz_xxyyyy_1, tg_zzzzzz_xxyyyz_0, tg_zzzzzz_xxyyyz_1, tg_zzzzzz_xxyyz_1, \
                                         tg_zzzzzz_xxyyzz_0, tg_zzzzzz_xxyyzz_1, tg_zzzzzz_xxyzz_1, tg_zzzzzz_xxyzzz_0, \
                                         tg_zzzzzz_xxyzzz_1, tg_zzzzzz_xxzzz_1, tg_zzzzzz_xxzzzz_0, tg_zzzzzz_xxzzzz_1, \
                                         tg_zzzzzz_xyyyy_1, tg_zzzzzz_xyyyyy_0, tg_zzzzzz_xyyyyy_1, tg_zzzzzz_xyyyyz_0, \
                                         tg_zzzzzz_xyyyyz_1, tg_zzzzzz_xyyyz_1, tg_zzzzzz_xyyyzz_0, tg_zzzzzz_xyyyzz_1, \
                                         tg_zzzzzz_xyyzz_1, tg_zzzzzz_xyyzzz_0, tg_zzzzzz_xyyzzz_1, tg_zzzzzz_xyzzz_1, \
                                         tg_zzzzzz_xyzzzz_0, tg_zzzzzz_xyzzzz_1, tg_zzzzzz_xzzzz_1, tg_zzzzzz_xzzzzz_0, \
                                         tg_zzzzzz_xzzzzz_1, tg_zzzzzz_yyyyy_1, tg_zzzzzz_yyyyyy_0, tg_zzzzzz_yyyyyy_1, \
                                         tg_zzzzzz_yyyyyz_0, tg_zzzzzz_yyyyyz_1, tg_zzzzzz_yyyyz_1, tg_zzzzzz_yyyyzz_0, \
                                         tg_zzzzzz_yyyyzz_1, tg_zzzzzz_yyyzz_1, tg_zzzzzz_yyyzzz_0, tg_zzzzzz_yyyzzz_1, \
                                         tg_zzzzzz_yyzzz_1, tg_zzzzzz_yyzzzz_0, tg_zzzzzz_yyzzzz_1, tg_zzzzzz_yzzzz_1, \
                                         tg_zzzzzz_yzzzzz_0, tg_zzzzzz_yzzzzz_1, tg_zzzzzz_zzzzz_1, tg_zzzzzz_zzzzzz_0, \
                                         tg_zzzzzz_zzzzzz_1, tg_zzzzzzz_xxxxxx_0, tg_zzzzzzz_xxxxxy_0, tg_zzzzzzz_xxxxxz_0, \
                                         tg_zzzzzzz_xxxxyy_0, tg_zzzzzzz_xxxxyz_0, tg_zzzzzzz_xxxxzz_0, tg_zzzzzzz_xxxyyy_0, \
                                         tg_zzzzzzz_xxxyyz_0, tg_zzzzzzz_xxxyzz_0, tg_zzzzzzz_xxxzzz_0, tg_zzzzzzz_xxyyyy_0, \
                                         tg_zzzzzzz_xxyyyz_0, tg_zzzzzzz_xxyyzz_0, tg_zzzzzzz_xxyzzz_0, tg_zzzzzzz_xxzzzz_0, \
                                         tg_zzzzzzz_xyyyyy_0, tg_zzzzzzz_xyyyyz_0, tg_zzzzzzz_xyyyzz_0, tg_zzzzzzz_xyyzzz_0, \
                                         tg_zzzzzzz_xyzzzz_0, tg_zzzzzzz_xzzzzz_0, tg_zzzzzzz_yyyyyy_0, tg_zzzzzzz_yyyyyz_0, \
                                         tg_zzzzzzz_yyyyzz_0, tg_zzzzzzz_yyyzzz_0, tg_zzzzzzz_yyzzzz_0, tg_zzzzzzz_yzzzzz_0, \
                                         tg_zzzzzzz_zzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyzzzz_yyyyyy_0[j] = pb_y * tg_yyzzzz_yyyyyy_0[j] + fr * tg_yyzzzz_yyyyyy_1[j] + fl1_fx * (tg_yzzzz_yyyyyy_0[j] - tg_yzzzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyzzzz_yyyyy_1[j];

                    tg_yyyzzzz_yyyyyz_0[j] = pb_y * tg_yyzzzz_yyyyyz_0[j] + fr * tg_yyzzzz_yyyyyz_1[j] + fl1_fx * (tg_yzzzz_yyyyyz_0[j] - tg_yzzzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzzzz_yyyyz_1[j];

                    tg_yyyzzzz_yyyyzz_0[j] = pb_y * tg_yyzzzz_yyyyzz_0[j] + fr * tg_yyzzzz_yyyyzz_1[j] + fl1_fx * (tg_yzzzz_yyyyzz_0[j] - tg_yzzzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzz_yyyzz_1[j];

                    tg_yyyzzzz_yyyzzz_0[j] = pb_y * tg_yyzzzz_yyyzzz_0[j] + fr * tg_yyzzzz_yyyzzz_1[j] + fl1_fx * (tg_yzzzz_yyyzzz_0[j] - tg_yzzzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_yyzzz_1[j];

                    tg_yyyzzzz_yyzzzz_0[j] = pb_y * tg_yyzzzz_yyzzzz_0[j] + fr * tg_yyzzzz_yyzzzz_1[j] + fl1_fx * (tg_yzzzz_yyzzzz_0[j] - tg_yzzzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_yzzzz_1[j];

                    tg_yyyzzzz_yzzzzz_0[j] = pb_y * tg_yyzzzz_yzzzzz_0[j] + fr * tg_yyzzzz_yzzzzz_1[j] + fl1_fx * (tg_yzzzz_yzzzzz_0[j] - tg_yzzzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_zzzzz_1[j];

                    tg_yyyzzzz_zzzzzz_0[j] = pb_y * tg_yyzzzz_zzzzzz_0[j] + fr * tg_yyzzzz_zzzzzz_1[j] + fl1_fx * (tg_yzzzz_zzzzzz_0[j] - tg_yzzzz_zzzzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxxxx_0[j] = pb_y * tg_yzzzzz_xxxxxx_0[j] + fr * tg_yzzzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxxx_0[j] - tg_zzzzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxxxy_0[j] = pb_y * tg_yzzzzz_xxxxxy_0[j] + fr * tg_yzzzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxxy_0[j] - tg_zzzzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxxxx_1[j];

                    tg_yyzzzzz_xxxxxz_0[j] = pb_y * tg_yzzzzz_xxxxxz_0[j] + fr * tg_yzzzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxxz_0[j] - tg_zzzzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxxyy_0[j] = pb_y * tg_yzzzzz_xxxxyy_0[j] + fr * tg_yzzzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxyy_0[j] - tg_zzzzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xxxxy_1[j];

                    tg_yyzzzzz_xxxxyz_0[j] = pb_y * tg_yzzzzz_xxxxyz_0[j] + fr * tg_yzzzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxyz_0[j] - tg_zzzzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxxxz_1[j];

                    tg_yyzzzzz_xxxxzz_0[j] = pb_y * tg_yzzzzz_xxxxzz_0[j] + fr * tg_yzzzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxxzz_0[j] - tg_zzzzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxyyy_0[j] = pb_y * tg_yzzzzz_xxxyyy_0[j] + fr * tg_yzzzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyyy_0[j] - tg_zzzzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_xxxyy_1[j];

                    tg_yyzzzzz_xxxyyz_0[j] = pb_y * tg_yzzzzz_xxxyyz_0[j] + fr * tg_yzzzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyyz_0[j] - tg_zzzzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xxxyz_1[j];

                    tg_yyzzzzz_xxxyzz_0[j] = pb_y * tg_yzzzzz_xxxyzz_0[j] + fr * tg_yzzzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxyzz_0[j] - tg_zzzzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxxzz_1[j];

                    tg_yyzzzzz_xxxzzz_0[j] = pb_y * tg_yzzzzz_xxxzzz_0[j] + fr * tg_yzzzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxzzz_0[j] - tg_zzzzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxyyyy_0[j] = pb_y * tg_yzzzzz_xxyyyy_0[j] + fr * tg_yzzzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyyy_0[j] - tg_zzzzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzz_xxyyy_1[j];

                    tg_yyzzzzz_xxyyyz_0[j] = pb_y * tg_yzzzzz_xxyyyz_0[j] + fr * tg_yzzzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyyz_0[j] - tg_zzzzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_xxyyz_1[j];

                    tg_yyzzzzz_xxyyzz_0[j] = pb_y * tg_yzzzzz_xxyyzz_0[j] + fr * tg_yzzzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyyzz_0[j] - tg_zzzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xxyzz_1[j];

                    tg_yyzzzzz_xxyzzz_0[j] = pb_y * tg_yzzzzz_xxyzzz_0[j] + fr * tg_yzzzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyzzz_0[j] - tg_zzzzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxzzz_1[j];

                    tg_yyzzzzz_xxzzzz_0[j] = pb_y * tg_yzzzzz_xxzzzz_0[j] + fr * tg_yzzzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxzzzz_0[j] - tg_zzzzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xyyyyy_0[j] = pb_y * tg_yzzzzz_xyyyyy_0[j] + fr * tg_yzzzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyyy_0[j] - tg_zzzzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzzzz_xyyyy_1[j];

                    tg_yyzzzzz_xyyyyz_0[j] = pb_y * tg_yzzzzz_xyyyyz_0[j] + fr * tg_yzzzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyyz_0[j] - tg_zzzzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzz_xyyyz_1[j];

                    tg_yyzzzzz_xyyyzz_0[j] = pb_y * tg_yzzzzz_xyyyzz_0[j] + fr * tg_yzzzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyyzz_0[j] - tg_zzzzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_xyyzz_1[j];

                    tg_yyzzzzz_xyyzzz_0[j] = pb_y * tg_yzzzzz_xyyzzz_0[j] + fr * tg_yzzzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyzzz_0[j] - tg_zzzzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xyzzz_1[j];

                    tg_yyzzzzz_xyzzzz_0[j] = pb_y * tg_yzzzzz_xyzzzz_0[j] + fr * tg_yzzzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyzzzz_0[j] - tg_zzzzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xzzzz_1[j];

                    tg_yyzzzzz_xzzzzz_0[j] = pb_y * tg_yzzzzz_xzzzzz_0[j] + fr * tg_yzzzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzzzzz_0[j] - tg_zzzzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_yyyyyy_0[j] = pb_y * tg_yzzzzz_yyyyyy_0[j] + fr * tg_yzzzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyyy_0[j] - tg_zzzzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzzzzz_yyyyy_1[j];

                    tg_yyzzzzz_yyyyyz_0[j] = pb_y * tg_yzzzzz_yyyyyz_0[j] + fr * tg_yzzzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyyz_0[j] - tg_zzzzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzzzz_yyyyz_1[j];

                    tg_yyzzzzz_yyyyzz_0[j] = pb_y * tg_yzzzzz_yyyyzz_0[j] + fr * tg_yzzzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyyzz_0[j] - tg_zzzzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzz_yyyzz_1[j];

                    tg_yyzzzzz_yyyzzz_0[j] = pb_y * tg_yzzzzz_yyyzzz_0[j] + fr * tg_yzzzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyzzz_0[j] - tg_zzzzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_yyzzz_1[j];

                    tg_yyzzzzz_yyzzzz_0[j] = pb_y * tg_yzzzzz_yyzzzz_0[j] + fr * tg_yzzzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyzzzz_0[j] - tg_zzzzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_yzzzz_1[j];

                    tg_yyzzzzz_yzzzzz_0[j] = pb_y * tg_yzzzzz_yzzzzz_0[j] + fr * tg_yzzzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzzzzz_0[j] - tg_zzzzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_zzzzz_1[j];

                    tg_yyzzzzz_zzzzzz_0[j] = pb_y * tg_yzzzzz_zzzzzz_0[j] + fr * tg_yzzzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzzzzz_0[j] - tg_zzzzz_zzzzzz_1[j] * fl1_fza);

                    tg_yzzzzzz_xxxxxx_0[j] = pb_y * tg_zzzzzz_xxxxxx_0[j] + fr * tg_zzzzzz_xxxxxx_1[j];

                    tg_yzzzzzz_xxxxxy_0[j] = pb_y * tg_zzzzzz_xxxxxy_0[j] + fr * tg_zzzzzz_xxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxxxx_1[j];

                    tg_yzzzzzz_xxxxxz_0[j] = pb_y * tg_zzzzzz_xxxxxz_0[j] + fr * tg_zzzzzz_xxxxxz_1[j];

                    tg_yzzzzzz_xxxxyy_0[j] = pb_y * tg_zzzzzz_xxxxyy_0[j] + fr * tg_zzzzzz_xxxxyy_1[j] + fl1_fxn * tg_zzzzzz_xxxxy_1[j];

                    tg_yzzzzzz_xxxxyz_0[j] = pb_y * tg_zzzzzz_xxxxyz_0[j] + fr * tg_zzzzzz_xxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxxxz_1[j];

                    tg_yzzzzzz_xxxxzz_0[j] = pb_y * tg_zzzzzz_xxxxzz_0[j] + fr * tg_zzzzzz_xxxxzz_1[j];

                    tg_yzzzzzz_xxxyyy_0[j] = pb_y * tg_zzzzzz_xxxyyy_0[j] + fr * tg_zzzzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxxyy_1[j];

                    tg_yzzzzzz_xxxyyz_0[j] = pb_y * tg_zzzzzz_xxxyyz_0[j] + fr * tg_zzzzzz_xxxyyz_1[j] + fl1_fxn * tg_zzzzzz_xxxyz_1[j];

                    tg_yzzzzzz_xxxyzz_0[j] = pb_y * tg_zzzzzz_xxxyzz_0[j] + fr * tg_zzzzzz_xxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxxzz_1[j];

                    tg_yzzzzzz_xxxzzz_0[j] = pb_y * tg_zzzzzz_xxxzzz_0[j] + fr * tg_zzzzzz_xxxzzz_1[j];

                    tg_yzzzzzz_xxyyyy_0[j] = pb_y * tg_zzzzzz_xxyyyy_0[j] + fr * tg_zzzzzz_xxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xxyyy_1[j];

                    tg_yzzzzzz_xxyyyz_0[j] = pb_y * tg_zzzzzz_xxyyyz_0[j] + fr * tg_zzzzzz_xxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxyyz_1[j];

                    tg_yzzzzzz_xxyyzz_0[j] = pb_y * tg_zzzzzz_xxyyzz_0[j] + fr * tg_zzzzzz_xxyyzz_1[j] + fl1_fxn * tg_zzzzzz_xxyzz_1[j];

                    tg_yzzzzzz_xxyzzz_0[j] = pb_y * tg_zzzzzz_xxyzzz_0[j] + fr * tg_zzzzzz_xxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxzzz_1[j];

                    tg_yzzzzzz_xxzzzz_0[j] = pb_y * tg_zzzzzz_xxzzzz_0[j] + fr * tg_zzzzzz_xxzzzz_1[j];

                    tg_yzzzzzz_xyyyyy_0[j] = pb_y * tg_zzzzzz_xyyyyy_0[j] + fr * tg_zzzzzz_xyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzzzz_xyyyy_1[j];

                    tg_yzzzzzz_xyyyyz_0[j] = pb_y * tg_zzzzzz_xyyyyz_0[j] + fr * tg_zzzzzz_xyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xyyyz_1[j];

                    tg_yzzzzzz_xyyyzz_0[j] = pb_y * tg_zzzzzz_xyyyzz_0[j] + fr * tg_zzzzzz_xyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xyyzz_1[j];

                    tg_yzzzzzz_xyyzzz_0[j] = pb_y * tg_zzzzzz_xyyzzz_0[j] + fr * tg_zzzzzz_xyyzzz_1[j] + fl1_fxn * tg_zzzzzz_xyzzz_1[j];

                    tg_yzzzzzz_xyzzzz_0[j] = pb_y * tg_zzzzzz_xyzzzz_0[j] + fr * tg_zzzzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xzzzz_1[j];

                    tg_yzzzzzz_xzzzzz_0[j] = pb_y * tg_zzzzzz_xzzzzz_0[j] + fr * tg_zzzzzz_xzzzzz_1[j];

                    tg_yzzzzzz_yyyyyy_0[j] = pb_y * tg_zzzzzz_yyyyyy_0[j] + fr * tg_zzzzzz_yyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzzzzz_yyyyy_1[j];

                    tg_yzzzzzz_yyyyyz_0[j] = pb_y * tg_zzzzzz_yyyyyz_0[j] + fr * tg_zzzzzz_yyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzzzzz_yyyyz_1[j];

                    tg_yzzzzzz_yyyyzz_0[j] = pb_y * tg_zzzzzz_yyyyzz_0[j] + fr * tg_zzzzzz_yyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_yyyzz_1[j];

                    tg_yzzzzzz_yyyzzz_0[j] = pb_y * tg_zzzzzz_yyyzzz_0[j] + fr * tg_zzzzzz_yyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_yyzzz_1[j];

                    tg_yzzzzzz_yyzzzz_0[j] = pb_y * tg_zzzzzz_yyzzzz_0[j] + fr * tg_zzzzzz_yyzzzz_1[j] + fl1_fxn * tg_zzzzzz_yzzzz_1[j];

                    tg_yzzzzzz_yzzzzz_0[j] = pb_y * tg_zzzzzz_yzzzzz_0[j] + fr * tg_zzzzzz_yzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zzzzz_1[j];

                    tg_yzzzzzz_zzzzzz_0[j] = pb_y * tg_zzzzzz_zzzzzz_0[j] + fr * tg_zzzzzz_zzzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzzz_xxxxxx_0[j] = pb_z * tg_zzzzzz_xxxxxx_0[j] + fr * tg_zzzzzz_xxxxxx_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxxx_0[j] - tg_zzzzz_xxxxxx_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxxxy_0[j] = pb_z * tg_zzzzzz_xxxxxy_0[j] + fr * tg_zzzzzz_xxxxxy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxxy_0[j] - tg_zzzzz_xxxxxy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxxxz_0[j] = pb_z * tg_zzzzzz_xxxxxz_0[j] + fr * tg_zzzzzz_xxxxxz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxxz_0[j] - tg_zzzzz_xxxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxxxx_1[j];

                    tg_zzzzzzz_xxxxyy_0[j] = pb_z * tg_zzzzzz_xxxxyy_0[j] + fr * tg_zzzzzz_xxxxyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxyy_0[j] - tg_zzzzz_xxxxyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxxyz_0[j] = pb_z * tg_zzzzzz_xxxxyz_0[j] + fr * tg_zzzzzz_xxxxyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxyz_0[j] - tg_zzzzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxxxy_1[j];

                    tg_zzzzzzz_xxxxzz_0[j] = pb_z * tg_zzzzzz_xxxxzz_0[j] + fr * tg_zzzzzz_xxxxzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxxzz_0[j] - tg_zzzzz_xxxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xxxxz_1[j];

                    tg_zzzzzzz_xxxyyy_0[j] = pb_z * tg_zzzzzz_xxxyyy_0[j] + fr * tg_zzzzzz_xxxyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxyyy_0[j] - tg_zzzzz_xxxyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxyyz_0[j] = pb_z * tg_zzzzzz_xxxyyz_0[j] + fr * tg_zzzzzz_xxxyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxyyz_0[j] - tg_zzzzz_xxxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxxyy_1[j];

                    tg_zzzzzzz_xxxyzz_0[j] = pb_z * tg_zzzzzz_xxxyzz_0[j] + fr * tg_zzzzzz_xxxyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxyzz_0[j] - tg_zzzzz_xxxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xxxyz_1[j];

                    tg_zzzzzzz_xxxzzz_0[j] = pb_z * tg_zzzzzz_xxxzzz_0[j] + fr * tg_zzzzzz_xxxzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxzzz_0[j] - tg_zzzzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_xxxzz_1[j];

                    tg_zzzzzzz_xxyyyy_0[j] = pb_z * tg_zzzzzz_xxyyyy_0[j] + fr * tg_zzzzzz_xxyyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyyyy_0[j] - tg_zzzzz_xxyyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxyyyz_0[j] = pb_z * tg_zzzzzz_xxyyyz_0[j] + fr * tg_zzzzzz_xxyyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyyyz_0[j] - tg_zzzzz_xxyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxyyy_1[j];

                    tg_zzzzzzz_xxyyzz_0[j] = pb_z * tg_zzzzzz_xxyyzz_0[j] + fr * tg_zzzzzz_xxyyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyyzz_0[j] - tg_zzzzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xxyyz_1[j];

                    tg_zzzzzzz_xxyzzz_0[j] = pb_z * tg_zzzzzz_xxyzzz_0[j] + fr * tg_zzzzzz_xxyzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyzzz_0[j] - tg_zzzzz_xxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_xxyzz_1[j];

                    tg_zzzzzzz_xxzzzz_0[j] = pb_z * tg_zzzzzz_xxzzzz_0[j] + fr * tg_zzzzzz_xxzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxzzzz_0[j] - tg_zzzzz_xxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzz_xxzzz_1[j];

                    tg_zzzzzzz_xyyyyy_0[j] = pb_z * tg_zzzzzz_xyyyyy_0[j] + fr * tg_zzzzzz_xyyyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyyyy_0[j] - tg_zzzzz_xyyyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xyyyyz_0[j] = pb_z * tg_zzzzzz_xyyyyz_0[j] + fr * tg_zzzzzz_xyyyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyyyz_0[j] - tg_zzzzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xyyyy_1[j];

                    tg_zzzzzzz_xyyyzz_0[j] = pb_z * tg_zzzzzz_xyyyzz_0[j] + fr * tg_zzzzzz_xyyyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyyzz_0[j] - tg_zzzzz_xyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xyyyz_1[j];

                    tg_zzzzzzz_xyyzzz_0[j] = pb_z * tg_zzzzzz_xyyzzz_0[j] + fr * tg_zzzzzz_xyyzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyzzz_0[j] - tg_zzzzz_xyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_xyyzz_1[j];

                    tg_zzzzzzz_xyzzzz_0[j] = pb_z * tg_zzzzzz_xyzzzz_0[j] + fr * tg_zzzzzz_xyzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyzzzz_0[j] - tg_zzzzz_xyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzz_xyzzz_1[j];

                    tg_zzzzzzz_xzzzzz_0[j] = pb_z * tg_zzzzzz_xzzzzz_0[j] + fr * tg_zzzzzz_xzzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xzzzzz_0[j] - tg_zzzzz_xzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzzzz_xzzzz_1[j];

                    tg_zzzzzzz_yyyyyy_0[j] = pb_z * tg_zzzzzz_yyyyyy_0[j] + fr * tg_zzzzzz_yyyyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyyyy_0[j] - tg_zzzzz_yyyyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_yyyyyz_0[j] = pb_z * tg_zzzzzz_yyyyyz_0[j] + fr * tg_zzzzzz_yyyyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyyyz_0[j] - tg_zzzzz_yyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_yyyyy_1[j];

                    tg_zzzzzzz_yyyyzz_0[j] = pb_z * tg_zzzzzz_yyyyzz_0[j] + fr * tg_zzzzzz_yyyyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyyzz_0[j] - tg_zzzzz_yyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_yyyyz_1[j];

                    tg_zzzzzzz_yyyzzz_0[j] = pb_z * tg_zzzzzz_yyyzzz_0[j] + fr * tg_zzzzzz_yyyzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyzzz_0[j] - tg_zzzzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_yyyzz_1[j];

                    tg_zzzzzzz_yyzzzz_0[j] = pb_z * tg_zzzzzz_yyzzzz_0[j] + fr * tg_zzzzzz_yyzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyzzzz_0[j] - tg_zzzzz_yyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzz_yyzzz_1[j];

                    tg_zzzzzzz_yzzzzz_0[j] = pb_z * tg_zzzzzz_yzzzzz_0[j] + fr * tg_zzzzzz_yzzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yzzzzz_0[j] - tg_zzzzz_yzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzzzz_yzzzz_1[j];

                    tg_zzzzzzz_zzzzzz_0[j] = pb_z * tg_zzzzzz_zzzzzz_0[j] + fr * tg_zzzzzz_zzzzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_zzzzzz_0[j] - tg_zzzzz_zzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzzzzz_zzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

