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

#include "ElectronRepulsionRecFuncForGI.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSGSI(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSGSI_0_84(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSGSI_84_168(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSGSI_168_252(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSI_252_336(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSGSI_336_420(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSGSI_0_84(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,84)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xxx_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx); 

                auto tg_xxx_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 1); 

                auto tg_xxx_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 2); 

                auto tg_xxx_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 3); 

                auto tg_xxx_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 4); 

                auto tg_xxx_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 5); 

                auto tg_xxx_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 6); 

                auto tg_xxx_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 7); 

                auto tg_xxx_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 8); 

                auto tg_xxx_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 9); 

                auto tg_xxx_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 10); 

                auto tg_xxx_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 11); 

                auto tg_xxx_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 12); 

                auto tg_xxx_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 13); 

                auto tg_xxx_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 14); 

                auto tg_xxx_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 15); 

                auto tg_xxx_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 16); 

                auto tg_xxx_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 17); 

                auto tg_xxx_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 18); 

                auto tg_xxx_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 19); 

                auto tg_xxx_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 20); 

                auto tg_xxx_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 21); 

                auto tg_xxx_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 22); 

                auto tg_xxx_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 23); 

                auto tg_xxx_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 24); 

                auto tg_xxx_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 25); 

                auto tg_xxx_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 26); 

                auto tg_xxx_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 27); 

                auto tg_xxy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 28); 

                auto tg_xxy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 29); 

                auto tg_xxy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 30); 

                auto tg_xxy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 31); 

                auto tg_xxy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 32); 

                auto tg_xxy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 33); 

                auto tg_xxy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 34); 

                auto tg_xxy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 35); 

                auto tg_xxy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 36); 

                auto tg_xxy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 37); 

                auto tg_xxy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 38); 

                auto tg_xxy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 39); 

                auto tg_xxy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 40); 

                auto tg_xxy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 41); 

                auto tg_xxy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 42); 

                auto tg_xxy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 43); 

                auto tg_xxy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 44); 

                auto tg_xxy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 45); 

                auto tg_xxy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 46); 

                auto tg_xxy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 47); 

                auto tg_xxy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 48); 

                auto tg_xxy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 49); 

                auto tg_xxy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 50); 

                auto tg_xxy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 51); 

                auto tg_xxy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 52); 

                auto tg_xxy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 53); 

                auto tg_xxy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 54); 

                auto tg_xxy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 55); 

                auto tg_xxz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 56); 

                auto tg_xxz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 57); 

                auto tg_xxz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 58); 

                auto tg_xxz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 59); 

                auto tg_xxz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 60); 

                auto tg_xxz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 61); 

                auto tg_xxz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 62); 

                auto tg_xxz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 63); 

                auto tg_xxz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 64); 

                auto tg_xxz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 65); 

                auto tg_xxz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 66); 

                auto tg_xxz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 67); 

                auto tg_xxz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 68); 

                auto tg_xxz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 69); 

                auto tg_xxz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 70); 

                auto tg_xxz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 71); 

                auto tg_xxz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 72); 

                auto tg_xxz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 73); 

                auto tg_xxz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 74); 

                auto tg_xxz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 75); 

                auto tg_xxz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 76); 

                auto tg_xxz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 77); 

                auto tg_xxz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 78); 

                auto tg_xxz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 79); 

                auto tg_xxz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 80); 

                auto tg_xxz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 81); 

                auto tg_xxz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 82); 

                auto tg_xxz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 83); 

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

                // set up pointers to integrals

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

                // Batch of Integrals (0,84)

                #pragma omp simd aligned(fxn, fza, tg_xx_xxxxxx_0, tg_xx_xxxxxx_1, tg_xx_xxxxxy_0, \
                                         tg_xx_xxxxxy_1, tg_xx_xxxxxz_0, tg_xx_xxxxxz_1, tg_xx_xxxxyy_0, tg_xx_xxxxyy_1, \
                                         tg_xx_xxxxyz_0, tg_xx_xxxxyz_1, tg_xx_xxxxzz_0, tg_xx_xxxxzz_1, tg_xx_xxxyyy_0, \
                                         tg_xx_xxxyyy_1, tg_xx_xxxyyz_0, tg_xx_xxxyyz_1, tg_xx_xxxyzz_0, tg_xx_xxxyzz_1, \
                                         tg_xx_xxxzzz_0, tg_xx_xxxzzz_1, tg_xx_xxyyyy_0, tg_xx_xxyyyy_1, tg_xx_xxyyyz_0, \
                                         tg_xx_xxyyyz_1, tg_xx_xxyyzz_0, tg_xx_xxyyzz_1, tg_xx_xxyzzz_0, tg_xx_xxyzzz_1, \
                                         tg_xx_xxzzzz_0, tg_xx_xxzzzz_1, tg_xx_xyyyyy_0, tg_xx_xyyyyy_1, tg_xx_xyyyyz_0, \
                                         tg_xx_xyyyyz_1, tg_xx_xyyyzz_0, tg_xx_xyyyzz_1, tg_xx_xyyzzz_0, tg_xx_xyyzzz_1, \
                                         tg_xx_xyzzzz_0, tg_xx_xyzzzz_1, tg_xx_xzzzzz_0, tg_xx_xzzzzz_1, tg_xx_yyyyyy_0, \
                                         tg_xx_yyyyyy_1, tg_xx_yyyyyz_0, tg_xx_yyyyyz_1, tg_xx_yyyyzz_0, tg_xx_yyyyzz_1, \
                                         tg_xx_yyyzzz_0, tg_xx_yyyzzz_1, tg_xx_yyzzzz_0, tg_xx_yyzzzz_1, tg_xx_yzzzzz_0, \
                                         tg_xx_yzzzzz_1, tg_xx_zzzzzz_0, tg_xx_zzzzzz_1, tg_xxx_xxxxx_1, tg_xxx_xxxxxx_0, \
                                         tg_xxx_xxxxxx_1, tg_xxx_xxxxxy_0, tg_xxx_xxxxxy_1, tg_xxx_xxxxxz_0, tg_xxx_xxxxxz_1, \
                                         tg_xxx_xxxxy_1, tg_xxx_xxxxyy_0, tg_xxx_xxxxyy_1, tg_xxx_xxxxyz_0, tg_xxx_xxxxyz_1, \
                                         tg_xxx_xxxxz_1, tg_xxx_xxxxzz_0, tg_xxx_xxxxzz_1, tg_xxx_xxxyy_1, tg_xxx_xxxyyy_0, \
                                         tg_xxx_xxxyyy_1, tg_xxx_xxxyyz_0, tg_xxx_xxxyyz_1, tg_xxx_xxxyz_1, tg_xxx_xxxyzz_0, \
                                         tg_xxx_xxxyzz_1, tg_xxx_xxxzz_1, tg_xxx_xxxzzz_0, tg_xxx_xxxzzz_1, tg_xxx_xxyyy_1, \
                                         tg_xxx_xxyyyy_0, tg_xxx_xxyyyy_1, tg_xxx_xxyyyz_0, tg_xxx_xxyyyz_1, tg_xxx_xxyyz_1, \
                                         tg_xxx_xxyyzz_0, tg_xxx_xxyyzz_1, tg_xxx_xxyzz_1, tg_xxx_xxyzzz_0, tg_xxx_xxyzzz_1, \
                                         tg_xxx_xxzzz_1, tg_xxx_xxzzzz_0, tg_xxx_xxzzzz_1, tg_xxx_xyyyy_1, tg_xxx_xyyyyy_0, \
                                         tg_xxx_xyyyyy_1, tg_xxx_xyyyyz_0, tg_xxx_xyyyyz_1, tg_xxx_xyyyz_1, tg_xxx_xyyyzz_0, \
                                         tg_xxx_xyyyzz_1, tg_xxx_xyyzz_1, tg_xxx_xyyzzz_0, tg_xxx_xyyzzz_1, tg_xxx_xyzzz_1, \
                                         tg_xxx_xyzzzz_0, tg_xxx_xyzzzz_1, tg_xxx_xzzzz_1, tg_xxx_xzzzzz_0, tg_xxx_xzzzzz_1, \
                                         tg_xxx_yyyyy_1, tg_xxx_yyyyyy_0, tg_xxx_yyyyyy_1, tg_xxx_yyyyyz_0, tg_xxx_yyyyyz_1, \
                                         tg_xxx_yyyyz_1, tg_xxx_yyyyzz_0, tg_xxx_yyyyzz_1, tg_xxx_yyyzz_1, tg_xxx_yyyzzz_0, \
                                         tg_xxx_yyyzzz_1, tg_xxx_yyzzz_1, tg_xxx_yyzzzz_0, tg_xxx_yyzzzz_1, tg_xxx_yzzzz_1, \
                                         tg_xxx_yzzzzz_0, tg_xxx_yzzzzz_1, tg_xxx_zzzzz_1, tg_xxx_zzzzzz_0, tg_xxx_zzzzzz_1, \
                                         tg_xxxx_xxxxxx_0, tg_xxxx_xxxxxy_0, tg_xxxx_xxxxxz_0, tg_xxxx_xxxxyy_0, \
                                         tg_xxxx_xxxxyz_0, tg_xxxx_xxxxzz_0, tg_xxxx_xxxyyy_0, tg_xxxx_xxxyyz_0, \
                                         tg_xxxx_xxxyzz_0, tg_xxxx_xxxzzz_0, tg_xxxx_xxyyyy_0, tg_xxxx_xxyyyz_0, \
                                         tg_xxxx_xxyyzz_0, tg_xxxx_xxyzzz_0, tg_xxxx_xxzzzz_0, tg_xxxx_xyyyyy_0, \
                                         tg_xxxx_xyyyyz_0, tg_xxxx_xyyyzz_0, tg_xxxx_xyyzzz_0, tg_xxxx_xyzzzz_0, \
                                         tg_xxxx_xzzzzz_0, tg_xxxx_yyyyyy_0, tg_xxxx_yyyyyz_0, tg_xxxx_yyyyzz_0, \
                                         tg_xxxx_yyyzzz_0, tg_xxxx_yyzzzz_0, tg_xxxx_yzzzzz_0, tg_xxxx_zzzzzz_0, \
                                         tg_xxxy_xxxxxx_0, tg_xxxy_xxxxxy_0, tg_xxxy_xxxxxz_0, tg_xxxy_xxxxyy_0, \
                                         tg_xxxy_xxxxyz_0, tg_xxxy_xxxxzz_0, tg_xxxy_xxxyyy_0, tg_xxxy_xxxyyz_0, \
                                         tg_xxxy_xxxyzz_0, tg_xxxy_xxxzzz_0, tg_xxxy_xxyyyy_0, tg_xxxy_xxyyyz_0, \
                                         tg_xxxy_xxyyzz_0, tg_xxxy_xxyzzz_0, tg_xxxy_xxzzzz_0, tg_xxxy_xyyyyy_0, \
                                         tg_xxxy_xyyyyz_0, tg_xxxy_xyyyzz_0, tg_xxxy_xyyzzz_0, tg_xxxy_xyzzzz_0, \
                                         tg_xxxy_xzzzzz_0, tg_xxxy_yyyyyy_0, tg_xxxy_yyyyyz_0, tg_xxxy_yyyyzz_0, \
                                         tg_xxxy_yyyzzz_0, tg_xxxy_yyzzzz_0, tg_xxxy_yzzzzz_0, tg_xxxy_zzzzzz_0, \
                                         tg_xxxz_xxxxxx_0, tg_xxxz_xxxxxy_0, tg_xxxz_xxxxxz_0, tg_xxxz_xxxxyy_0, \
                                         tg_xxxz_xxxxyz_0, tg_xxxz_xxxxzz_0, tg_xxxz_xxxyyy_0, tg_xxxz_xxxyyz_0, \
                                         tg_xxxz_xxxyzz_0, tg_xxxz_xxxzzz_0, tg_xxxz_xxyyyy_0, tg_xxxz_xxyyyz_0, \
                                         tg_xxxz_xxyyzz_0, tg_xxxz_xxyzzz_0, tg_xxxz_xxzzzz_0, tg_xxxz_xyyyyy_0, \
                                         tg_xxxz_xyyyyz_0, tg_xxxz_xyyyzz_0, tg_xxxz_xyyzzz_0, tg_xxxz_xyzzzz_0, \
                                         tg_xxxz_xzzzzz_0, tg_xxxz_yyyyyy_0, tg_xxxz_yyyyyz_0, tg_xxxz_yyyyzz_0, \
                                         tg_xxxz_yyyzzz_0, tg_xxxz_yyzzzz_0, tg_xxxz_yzzzzz_0, tg_xxxz_zzzzzz_0, \
                                         tg_xxy_xxxxx_1, tg_xxy_xxxxxx_0, tg_xxy_xxxxxx_1, tg_xxy_xxxxxy_0, tg_xxy_xxxxxy_1, \
                                         tg_xxy_xxxxxz_0, tg_xxy_xxxxxz_1, tg_xxy_xxxxy_1, tg_xxy_xxxxyy_0, tg_xxy_xxxxyy_1, \
                                         tg_xxy_xxxxyz_0, tg_xxy_xxxxyz_1, tg_xxy_xxxxz_1, tg_xxy_xxxxzz_0, tg_xxy_xxxxzz_1, \
                                         tg_xxy_xxxyy_1, tg_xxy_xxxyyy_0, tg_xxy_xxxyyy_1, tg_xxy_xxxyyz_0, tg_xxy_xxxyyz_1, \
                                         tg_xxy_xxxyz_1, tg_xxy_xxxyzz_0, tg_xxy_xxxyzz_1, tg_xxy_xxxzz_1, tg_xxy_xxxzzz_0, \
                                         tg_xxy_xxxzzz_1, tg_xxy_xxyyy_1, tg_xxy_xxyyyy_0, tg_xxy_xxyyyy_1, tg_xxy_xxyyyz_0, \
                                         tg_xxy_xxyyyz_1, tg_xxy_xxyyz_1, tg_xxy_xxyyzz_0, tg_xxy_xxyyzz_1, tg_xxy_xxyzz_1, \
                                         tg_xxy_xxyzzz_0, tg_xxy_xxyzzz_1, tg_xxy_xxzzz_1, tg_xxy_xxzzzz_0, tg_xxy_xxzzzz_1, \
                                         tg_xxy_xyyyy_1, tg_xxy_xyyyyy_0, tg_xxy_xyyyyy_1, tg_xxy_xyyyyz_0, tg_xxy_xyyyyz_1, \
                                         tg_xxy_xyyyz_1, tg_xxy_xyyyzz_0, tg_xxy_xyyyzz_1, tg_xxy_xyyzz_1, tg_xxy_xyyzzz_0, \
                                         tg_xxy_xyyzzz_1, tg_xxy_xyzzz_1, tg_xxy_xyzzzz_0, tg_xxy_xyzzzz_1, tg_xxy_xzzzz_1, \
                                         tg_xxy_xzzzzz_0, tg_xxy_xzzzzz_1, tg_xxy_yyyyy_1, tg_xxy_yyyyyy_0, tg_xxy_yyyyyy_1, \
                                         tg_xxy_yyyyyz_0, tg_xxy_yyyyyz_1, tg_xxy_yyyyz_1, tg_xxy_yyyyzz_0, tg_xxy_yyyyzz_1, \
                                         tg_xxy_yyyzz_1, tg_xxy_yyyzzz_0, tg_xxy_yyyzzz_1, tg_xxy_yyzzz_1, tg_xxy_yyzzzz_0, \
                                         tg_xxy_yyzzzz_1, tg_xxy_yzzzz_1, tg_xxy_yzzzzz_0, tg_xxy_yzzzzz_1, tg_xxy_zzzzz_1, \
                                         tg_xxy_zzzzzz_0, tg_xxy_zzzzzz_1, tg_xxz_xxxxx_1, tg_xxz_xxxxxx_0, tg_xxz_xxxxxx_1, \
                                         tg_xxz_xxxxxy_0, tg_xxz_xxxxxy_1, tg_xxz_xxxxxz_0, tg_xxz_xxxxxz_1, tg_xxz_xxxxy_1, \
                                         tg_xxz_xxxxyy_0, tg_xxz_xxxxyy_1, tg_xxz_xxxxyz_0, tg_xxz_xxxxyz_1, tg_xxz_xxxxz_1, \
                                         tg_xxz_xxxxzz_0, tg_xxz_xxxxzz_1, tg_xxz_xxxyy_1, tg_xxz_xxxyyy_0, tg_xxz_xxxyyy_1, \
                                         tg_xxz_xxxyyz_0, tg_xxz_xxxyyz_1, tg_xxz_xxxyz_1, tg_xxz_xxxyzz_0, tg_xxz_xxxyzz_1, \
                                         tg_xxz_xxxzz_1, tg_xxz_xxxzzz_0, tg_xxz_xxxzzz_1, tg_xxz_xxyyy_1, tg_xxz_xxyyyy_0, \
                                         tg_xxz_xxyyyy_1, tg_xxz_xxyyyz_0, tg_xxz_xxyyyz_1, tg_xxz_xxyyz_1, tg_xxz_xxyyzz_0, \
                                         tg_xxz_xxyyzz_1, tg_xxz_xxyzz_1, tg_xxz_xxyzzz_0, tg_xxz_xxyzzz_1, tg_xxz_xxzzz_1, \
                                         tg_xxz_xxzzzz_0, tg_xxz_xxzzzz_1, tg_xxz_xyyyy_1, tg_xxz_xyyyyy_0, tg_xxz_xyyyyy_1, \
                                         tg_xxz_xyyyyz_0, tg_xxz_xyyyyz_1, tg_xxz_xyyyz_1, tg_xxz_xyyyzz_0, tg_xxz_xyyyzz_1, \
                                         tg_xxz_xyyzz_1, tg_xxz_xyyzzz_0, tg_xxz_xyyzzz_1, tg_xxz_xyzzz_1, tg_xxz_xyzzzz_0, \
                                         tg_xxz_xyzzzz_1, tg_xxz_xzzzz_1, tg_xxz_xzzzzz_0, tg_xxz_xzzzzz_1, tg_xxz_yyyyy_1, \
                                         tg_xxz_yyyyyy_0, tg_xxz_yyyyyy_1, tg_xxz_yyyyyz_0, tg_xxz_yyyyyz_1, tg_xxz_yyyyz_1, \
                                         tg_xxz_yyyyzz_0, tg_xxz_yyyyzz_1, tg_xxz_yyyzz_1, tg_xxz_yyyzzz_0, tg_xxz_yyyzzz_1, \
                                         tg_xxz_yyzzz_1, tg_xxz_yyzzzz_0, tg_xxz_yyzzzz_1, tg_xxz_yzzzz_1, tg_xxz_yzzzzz_0, \
                                         tg_xxz_yzzzzz_1, tg_xxz_zzzzz_1, tg_xxz_zzzzzz_0, tg_xxz_zzzzzz_1, tg_xy_xxxxxx_0, \
                                         tg_xy_xxxxxx_1, tg_xy_xxxxxy_0, tg_xy_xxxxxy_1, tg_xy_xxxxxz_0, tg_xy_xxxxxz_1, \
                                         tg_xy_xxxxyy_0, tg_xy_xxxxyy_1, tg_xy_xxxxyz_0, tg_xy_xxxxyz_1, tg_xy_xxxxzz_0, \
                                         tg_xy_xxxxzz_1, tg_xy_xxxyyy_0, tg_xy_xxxyyy_1, tg_xy_xxxyyz_0, tg_xy_xxxyyz_1, \
                                         tg_xy_xxxyzz_0, tg_xy_xxxyzz_1, tg_xy_xxxzzz_0, tg_xy_xxxzzz_1, tg_xy_xxyyyy_0, \
                                         tg_xy_xxyyyy_1, tg_xy_xxyyyz_0, tg_xy_xxyyyz_1, tg_xy_xxyyzz_0, tg_xy_xxyyzz_1, \
                                         tg_xy_xxyzzz_0, tg_xy_xxyzzz_1, tg_xy_xxzzzz_0, tg_xy_xxzzzz_1, tg_xy_xyyyyy_0, \
                                         tg_xy_xyyyyy_1, tg_xy_xyyyyz_0, tg_xy_xyyyyz_1, tg_xy_xyyyzz_0, tg_xy_xyyyzz_1, \
                                         tg_xy_xyyzzz_0, tg_xy_xyyzzz_1, tg_xy_xyzzzz_0, tg_xy_xyzzzz_1, tg_xy_xzzzzz_0, \
                                         tg_xy_xzzzzz_1, tg_xy_yyyyyy_0, tg_xy_yyyyyy_1, tg_xy_yyyyyz_0, tg_xy_yyyyyz_1, \
                                         tg_xy_yyyyzz_0, tg_xy_yyyyzz_1, tg_xy_yyyzzz_0, tg_xy_yyyzzz_1, tg_xy_yyzzzz_0, \
                                         tg_xy_yyzzzz_1, tg_xy_yzzzzz_0, tg_xy_yzzzzz_1, tg_xy_zzzzzz_0, tg_xy_zzzzzz_1, \
                                         tg_xz_xxxxxx_0, tg_xz_xxxxxx_1, tg_xz_xxxxxy_0, tg_xz_xxxxxy_1, tg_xz_xxxxxz_0, \
                                         tg_xz_xxxxxz_1, tg_xz_xxxxyy_0, tg_xz_xxxxyy_1, tg_xz_xxxxyz_0, tg_xz_xxxxyz_1, \
                                         tg_xz_xxxxzz_0, tg_xz_xxxxzz_1, tg_xz_xxxyyy_0, tg_xz_xxxyyy_1, tg_xz_xxxyyz_0, \
                                         tg_xz_xxxyyz_1, tg_xz_xxxyzz_0, tg_xz_xxxyzz_1, tg_xz_xxxzzz_0, tg_xz_xxxzzz_1, \
                                         tg_xz_xxyyyy_0, tg_xz_xxyyyy_1, tg_xz_xxyyyz_0, tg_xz_xxyyyz_1, tg_xz_xxyyzz_0, \
                                         tg_xz_xxyyzz_1, tg_xz_xxyzzz_0, tg_xz_xxyzzz_1, tg_xz_xxzzzz_0, tg_xz_xxzzzz_1, \
                                         tg_xz_xyyyyy_0, tg_xz_xyyyyy_1, tg_xz_xyyyyz_0, tg_xz_xyyyyz_1, tg_xz_xyyyzz_0, \
                                         tg_xz_xyyyzz_1, tg_xz_xyyzzz_0, tg_xz_xyyzzz_1, tg_xz_xyzzzz_0, tg_xz_xyzzzz_1, \
                                         tg_xz_xzzzzz_0, tg_xz_xzzzzz_1, tg_xz_yyyyyy_0, tg_xz_yyyyyy_1, tg_xz_yyyyyz_0, \
                                         tg_xz_yyyyyz_1, tg_xz_yyyyzz_0, tg_xz_yyyyzz_1, tg_xz_yyyzzz_0, tg_xz_yyyzzz_1, \
                                         tg_xz_yyzzzz_0, tg_xz_yyzzzz_1, tg_xz_yzzzzz_0, tg_xz_yzzzzz_1, tg_xz_zzzzzz_0, \
                                         tg_xz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxx_xxxxxx_0[j] = pb_x * tg_xxx_xxxxxx_0[j] + fr * tg_xxx_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxx_0[j] - tg_xx_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxx_xxxxx_1[j];

                    tg_xxxx_xxxxxy_0[j] = pb_x * tg_xxx_xxxxxy_0[j] + fr * tg_xxx_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxy_0[j] - tg_xx_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxx_xxxxy_1[j];

                    tg_xxxx_xxxxxz_0[j] = pb_x * tg_xxx_xxxxxz_0[j] + fr * tg_xxx_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxxz_0[j] - tg_xx_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxx_xxxxz_1[j];

                    tg_xxxx_xxxxyy_0[j] = pb_x * tg_xxx_xxxxyy_0[j] + fr * tg_xxx_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxyy_0[j] - tg_xx_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxyy_1[j];

                    tg_xxxx_xxxxyz_0[j] = pb_x * tg_xxx_xxxxyz_0[j] + fr * tg_xxx_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxyz_0[j] - tg_xx_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxyz_1[j];

                    tg_xxxx_xxxxzz_0[j] = pb_x * tg_xxx_xxxxzz_0[j] + fr * tg_xxx_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxxzz_0[j] - tg_xx_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxx_xxxzz_1[j];

                    tg_xxxx_xxxyyy_0[j] = pb_x * tg_xxx_xxxyyy_0[j] + fr * tg_xxx_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyyy_0[j] - tg_xx_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyyy_1[j];

                    tg_xxxx_xxxyyz_0[j] = pb_x * tg_xxx_xxxyyz_0[j] + fr * tg_xxx_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyyz_0[j] - tg_xx_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyyz_1[j];

                    tg_xxxx_xxxyzz_0[j] = pb_x * tg_xxx_xxxyzz_0[j] + fr * tg_xxx_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxyzz_0[j] - tg_xx_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxyzz_1[j];

                    tg_xxxx_xxxzzz_0[j] = pb_x * tg_xxx_xxxzzz_0[j] + fr * tg_xxx_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxxzzz_0[j] - tg_xx_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xxzzz_1[j];

                    tg_xxxx_xxyyyy_0[j] = pb_x * tg_xxx_xxyyyy_0[j] + fr * tg_xxx_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyyy_0[j] - tg_xx_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyyy_1[j];

                    tg_xxxx_xxyyyz_0[j] = pb_x * tg_xxx_xxyyyz_0[j] + fr * tg_xxx_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyyz_0[j] - tg_xx_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyyz_1[j];

                    tg_xxxx_xxyyzz_0[j] = pb_x * tg_xxx_xxyyzz_0[j] + fr * tg_xxx_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyyzz_0[j] - tg_xx_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyyzz_1[j];

                    tg_xxxx_xxyzzz_0[j] = pb_x * tg_xxx_xxyzzz_0[j] + fr * tg_xxx_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxyzzz_0[j] - tg_xx_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xyzzz_1[j];

                    tg_xxxx_xxzzzz_0[j] = pb_x * tg_xxx_xxzzzz_0[j] + fr * tg_xxx_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xxzzzz_0[j] - tg_xx_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xzzzz_1[j];

                    tg_xxxx_xyyyyy_0[j] = pb_x * tg_xxx_xyyyyy_0[j] + fr * tg_xxx_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyyy_0[j] - tg_xx_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyyy_1[j];

                    tg_xxxx_xyyyyz_0[j] = pb_x * tg_xxx_xyyyyz_0[j] + fr * tg_xxx_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyyz_0[j] - tg_xx_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyyz_1[j];

                    tg_xxxx_xyyyzz_0[j] = pb_x * tg_xxx_xyyyzz_0[j] + fr * tg_xxx_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyyzz_0[j] - tg_xx_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyyzz_1[j];

                    tg_xxxx_xyyzzz_0[j] = pb_x * tg_xxx_xyyzzz_0[j] + fr * tg_xxx_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyyzzz_0[j] - tg_xx_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yyzzz_1[j];

                    tg_xxxx_xyzzzz_0[j] = pb_x * tg_xxx_xyzzzz_0[j] + fr * tg_xxx_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xyzzzz_0[j] - tg_xx_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yzzzz_1[j];

                    tg_xxxx_xzzzzz_0[j] = pb_x * tg_xxx_xzzzzz_0[j] + fr * tg_xxx_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_xzzzzz_0[j] - tg_xx_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_zzzzz_1[j];

                    tg_xxxx_yyyyyy_0[j] = pb_x * tg_xxx_yyyyyy_0[j] + fr * tg_xxx_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyyy_0[j] - tg_xx_yyyyyy_1[j] * fl1_fza);

                    tg_xxxx_yyyyyz_0[j] = pb_x * tg_xxx_yyyyyz_0[j] + fr * tg_xxx_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyyz_0[j] - tg_xx_yyyyyz_1[j] * fl1_fza);

                    tg_xxxx_yyyyzz_0[j] = pb_x * tg_xxx_yyyyzz_0[j] + fr * tg_xxx_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyyzz_0[j] - tg_xx_yyyyzz_1[j] * fl1_fza);

                    tg_xxxx_yyyzzz_0[j] = pb_x * tg_xxx_yyyzzz_0[j] + fr * tg_xxx_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyyzzz_0[j] - tg_xx_yyyzzz_1[j] * fl1_fza);

                    tg_xxxx_yyzzzz_0[j] = pb_x * tg_xxx_yyzzzz_0[j] + fr * tg_xxx_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yyzzzz_0[j] - tg_xx_yyzzzz_1[j] * fl1_fza);

                    tg_xxxx_yzzzzz_0[j] = pb_x * tg_xxx_yzzzzz_0[j] + fr * tg_xxx_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_yzzzzz_0[j] - tg_xx_yzzzzz_1[j] * fl1_fza);

                    tg_xxxx_zzzzzz_0[j] = pb_x * tg_xxx_zzzzzz_0[j] + fr * tg_xxx_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xx_zzzzzz_0[j] - tg_xx_zzzzzz_1[j] * fl1_fza);

                    tg_xxxy_xxxxxx_0[j] = pb_x * tg_xxy_xxxxxx_0[j] + fr * tg_xxy_xxxxxx_1[j] + fl1_fx * (tg_xy_xxxxxx_0[j] - tg_xy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxy_xxxxx_1[j];

                    tg_xxxy_xxxxxy_0[j] = pb_x * tg_xxy_xxxxxy_0[j] + fr * tg_xxy_xxxxxy_1[j] + fl1_fx * (tg_xy_xxxxxy_0[j] - tg_xy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxy_xxxxy_1[j];

                    tg_xxxy_xxxxxz_0[j] = pb_x * tg_xxy_xxxxxz_0[j] + fr * tg_xxy_xxxxxz_1[j] + fl1_fx * (tg_xy_xxxxxz_0[j] - tg_xy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxy_xxxxz_1[j];

                    tg_xxxy_xxxxyy_0[j] = pb_x * tg_xxy_xxxxyy_0[j] + fr * tg_xxy_xxxxyy_1[j] + fl1_fx * (tg_xy_xxxxyy_0[j] - tg_xy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxyy_1[j];

                    tg_xxxy_xxxxyz_0[j] = pb_x * tg_xxy_xxxxyz_0[j] + fr * tg_xxy_xxxxyz_1[j] + fl1_fx * (tg_xy_xxxxyz_0[j] - tg_xy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxyz_1[j];

                    tg_xxxy_xxxxzz_0[j] = pb_x * tg_xxy_xxxxzz_0[j] + fr * tg_xxy_xxxxzz_1[j] + fl1_fx * (tg_xy_xxxxzz_0[j] - tg_xy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxy_xxxzz_1[j];

                    tg_xxxy_xxxyyy_0[j] = pb_x * tg_xxy_xxxyyy_0[j] + fr * tg_xxy_xxxyyy_1[j] + fl1_fx * (tg_xy_xxxyyy_0[j] - tg_xy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyyy_1[j];

                    tg_xxxy_xxxyyz_0[j] = pb_x * tg_xxy_xxxyyz_0[j] + fr * tg_xxy_xxxyyz_1[j] + fl1_fx * (tg_xy_xxxyyz_0[j] - tg_xy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyyz_1[j];

                    tg_xxxy_xxxyzz_0[j] = pb_x * tg_xxy_xxxyzz_0[j] + fr * tg_xxy_xxxyzz_1[j] + fl1_fx * (tg_xy_xxxyzz_0[j] - tg_xy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxyzz_1[j];

                    tg_xxxy_xxxzzz_0[j] = pb_x * tg_xxy_xxxzzz_0[j] + fr * tg_xxy_xxxzzz_1[j] + fl1_fx * (tg_xy_xxxzzz_0[j] - tg_xy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xxzzz_1[j];

                    tg_xxxy_xxyyyy_0[j] = pb_x * tg_xxy_xxyyyy_0[j] + fr * tg_xxy_xxyyyy_1[j] + fl1_fx * (tg_xy_xxyyyy_0[j] - tg_xy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyyy_1[j];

                    tg_xxxy_xxyyyz_0[j] = pb_x * tg_xxy_xxyyyz_0[j] + fr * tg_xxy_xxyyyz_1[j] + fl1_fx * (tg_xy_xxyyyz_0[j] - tg_xy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyyz_1[j];

                    tg_xxxy_xxyyzz_0[j] = pb_x * tg_xxy_xxyyzz_0[j] + fr * tg_xxy_xxyyzz_1[j] + fl1_fx * (tg_xy_xxyyzz_0[j] - tg_xy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyyzz_1[j];

                    tg_xxxy_xxyzzz_0[j] = pb_x * tg_xxy_xxyzzz_0[j] + fr * tg_xxy_xxyzzz_1[j] + fl1_fx * (tg_xy_xxyzzz_0[j] - tg_xy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xyzzz_1[j];

                    tg_xxxy_xxzzzz_0[j] = pb_x * tg_xxy_xxzzzz_0[j] + fr * tg_xxy_xxzzzz_1[j] + fl1_fx * (tg_xy_xxzzzz_0[j] - tg_xy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xzzzz_1[j];

                    tg_xxxy_xyyyyy_0[j] = pb_x * tg_xxy_xyyyyy_0[j] + fr * tg_xxy_xyyyyy_1[j] + fl1_fx * (tg_xy_xyyyyy_0[j] - tg_xy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyyy_1[j];

                    tg_xxxy_xyyyyz_0[j] = pb_x * tg_xxy_xyyyyz_0[j] + fr * tg_xxy_xyyyyz_1[j] + fl1_fx * (tg_xy_xyyyyz_0[j] - tg_xy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyyz_1[j];

                    tg_xxxy_xyyyzz_0[j] = pb_x * tg_xxy_xyyyzz_0[j] + fr * tg_xxy_xyyyzz_1[j] + fl1_fx * (tg_xy_xyyyzz_0[j] - tg_xy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyyzz_1[j];

                    tg_xxxy_xyyzzz_0[j] = pb_x * tg_xxy_xyyzzz_0[j] + fr * tg_xxy_xyyzzz_1[j] + fl1_fx * (tg_xy_xyyzzz_0[j] - tg_xy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yyzzz_1[j];

                    tg_xxxy_xyzzzz_0[j] = pb_x * tg_xxy_xyzzzz_0[j] + fr * tg_xxy_xyzzzz_1[j] + fl1_fx * (tg_xy_xyzzzz_0[j] - tg_xy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yzzzz_1[j];

                    tg_xxxy_xzzzzz_0[j] = pb_x * tg_xxy_xzzzzz_0[j] + fr * tg_xxy_xzzzzz_1[j] + fl1_fx * (tg_xy_xzzzzz_0[j] - tg_xy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_zzzzz_1[j];

                    tg_xxxy_yyyyyy_0[j] = pb_x * tg_xxy_yyyyyy_0[j] + fr * tg_xxy_yyyyyy_1[j] + fl1_fx * (tg_xy_yyyyyy_0[j] - tg_xy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxy_yyyyyz_0[j] = pb_x * tg_xxy_yyyyyz_0[j] + fr * tg_xxy_yyyyyz_1[j] + fl1_fx * (tg_xy_yyyyyz_0[j] - tg_xy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxy_yyyyzz_0[j] = pb_x * tg_xxy_yyyyzz_0[j] + fr * tg_xxy_yyyyzz_1[j] + fl1_fx * (tg_xy_yyyyzz_0[j] - tg_xy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxy_yyyzzz_0[j] = pb_x * tg_xxy_yyyzzz_0[j] + fr * tg_xxy_yyyzzz_1[j] + fl1_fx * (tg_xy_yyyzzz_0[j] - tg_xy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxy_yyzzzz_0[j] = pb_x * tg_xxy_yyzzzz_0[j] + fr * tg_xxy_yyzzzz_1[j] + fl1_fx * (tg_xy_yyzzzz_0[j] - tg_xy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxy_yzzzzz_0[j] = pb_x * tg_xxy_yzzzzz_0[j] + fr * tg_xxy_yzzzzz_1[j] + fl1_fx * (tg_xy_yzzzzz_0[j] - tg_xy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxy_zzzzzz_0[j] = pb_x * tg_xxy_zzzzzz_0[j] + fr * tg_xxy_zzzzzz_1[j] + fl1_fx * (tg_xy_zzzzzz_0[j] - tg_xy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxz_xxxxxx_0[j] = pb_x * tg_xxz_xxxxxx_0[j] + fr * tg_xxz_xxxxxx_1[j] + fl1_fx * (tg_xz_xxxxxx_0[j] - tg_xz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxz_xxxxx_1[j];

                    tg_xxxz_xxxxxy_0[j] = pb_x * tg_xxz_xxxxxy_0[j] + fr * tg_xxz_xxxxxy_1[j] + fl1_fx * (tg_xz_xxxxxy_0[j] - tg_xz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxz_xxxxy_1[j];

                    tg_xxxz_xxxxxz_0[j] = pb_x * tg_xxz_xxxxxz_0[j] + fr * tg_xxz_xxxxxz_1[j] + fl1_fx * (tg_xz_xxxxxz_0[j] - tg_xz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxz_xxxxz_1[j];

                    tg_xxxz_xxxxyy_0[j] = pb_x * tg_xxz_xxxxyy_0[j] + fr * tg_xxz_xxxxyy_1[j] + fl1_fx * (tg_xz_xxxxyy_0[j] - tg_xz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxyy_1[j];

                    tg_xxxz_xxxxyz_0[j] = pb_x * tg_xxz_xxxxyz_0[j] + fr * tg_xxz_xxxxyz_1[j] + fl1_fx * (tg_xz_xxxxyz_0[j] - tg_xz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxyz_1[j];

                    tg_xxxz_xxxxzz_0[j] = pb_x * tg_xxz_xxxxzz_0[j] + fr * tg_xxz_xxxxzz_1[j] + fl1_fx * (tg_xz_xxxxzz_0[j] - tg_xz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxz_xxxzz_1[j];

                    tg_xxxz_xxxyyy_0[j] = pb_x * tg_xxz_xxxyyy_0[j] + fr * tg_xxz_xxxyyy_1[j] + fl1_fx * (tg_xz_xxxyyy_0[j] - tg_xz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyyy_1[j];

                    tg_xxxz_xxxyyz_0[j] = pb_x * tg_xxz_xxxyyz_0[j] + fr * tg_xxz_xxxyyz_1[j] + fl1_fx * (tg_xz_xxxyyz_0[j] - tg_xz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyyz_1[j];

                    tg_xxxz_xxxyzz_0[j] = pb_x * tg_xxz_xxxyzz_0[j] + fr * tg_xxz_xxxyzz_1[j] + fl1_fx * (tg_xz_xxxyzz_0[j] - tg_xz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxyzz_1[j];

                    tg_xxxz_xxxzzz_0[j] = pb_x * tg_xxz_xxxzzz_0[j] + fr * tg_xxz_xxxzzz_1[j] + fl1_fx * (tg_xz_xxxzzz_0[j] - tg_xz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xxzzz_1[j];

                    tg_xxxz_xxyyyy_0[j] = pb_x * tg_xxz_xxyyyy_0[j] + fr * tg_xxz_xxyyyy_1[j] + fl1_fx * (tg_xz_xxyyyy_0[j] - tg_xz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyyy_1[j];

                    tg_xxxz_xxyyyz_0[j] = pb_x * tg_xxz_xxyyyz_0[j] + fr * tg_xxz_xxyyyz_1[j] + fl1_fx * (tg_xz_xxyyyz_0[j] - tg_xz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyyz_1[j];

                    tg_xxxz_xxyyzz_0[j] = pb_x * tg_xxz_xxyyzz_0[j] + fr * tg_xxz_xxyyzz_1[j] + fl1_fx * (tg_xz_xxyyzz_0[j] - tg_xz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyyzz_1[j];

                    tg_xxxz_xxyzzz_0[j] = pb_x * tg_xxz_xxyzzz_0[j] + fr * tg_xxz_xxyzzz_1[j] + fl1_fx * (tg_xz_xxyzzz_0[j] - tg_xz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xyzzz_1[j];

                    tg_xxxz_xxzzzz_0[j] = pb_x * tg_xxz_xxzzzz_0[j] + fr * tg_xxz_xxzzzz_1[j] + fl1_fx * (tg_xz_xxzzzz_0[j] - tg_xz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xzzzz_1[j];

                    tg_xxxz_xyyyyy_0[j] = pb_x * tg_xxz_xyyyyy_0[j] + fr * tg_xxz_xyyyyy_1[j] + fl1_fx * (tg_xz_xyyyyy_0[j] - tg_xz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyyy_1[j];

                    tg_xxxz_xyyyyz_0[j] = pb_x * tg_xxz_xyyyyz_0[j] + fr * tg_xxz_xyyyyz_1[j] + fl1_fx * (tg_xz_xyyyyz_0[j] - tg_xz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyyz_1[j];

                    tg_xxxz_xyyyzz_0[j] = pb_x * tg_xxz_xyyyzz_0[j] + fr * tg_xxz_xyyyzz_1[j] + fl1_fx * (tg_xz_xyyyzz_0[j] - tg_xz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyyzz_1[j];

                    tg_xxxz_xyyzzz_0[j] = pb_x * tg_xxz_xyyzzz_0[j] + fr * tg_xxz_xyyzzz_1[j] + fl1_fx * (tg_xz_xyyzzz_0[j] - tg_xz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yyzzz_1[j];

                    tg_xxxz_xyzzzz_0[j] = pb_x * tg_xxz_xyzzzz_0[j] + fr * tg_xxz_xyzzzz_1[j] + fl1_fx * (tg_xz_xyzzzz_0[j] - tg_xz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yzzzz_1[j];

                    tg_xxxz_xzzzzz_0[j] = pb_x * tg_xxz_xzzzzz_0[j] + fr * tg_xxz_xzzzzz_1[j] + fl1_fx * (tg_xz_xzzzzz_0[j] - tg_xz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_zzzzz_1[j];

                    tg_xxxz_yyyyyy_0[j] = pb_x * tg_xxz_yyyyyy_0[j] + fr * tg_xxz_yyyyyy_1[j] + fl1_fx * (tg_xz_yyyyyy_0[j] - tg_xz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxz_yyyyyz_0[j] = pb_x * tg_xxz_yyyyyz_0[j] + fr * tg_xxz_yyyyyz_1[j] + fl1_fx * (tg_xz_yyyyyz_0[j] - tg_xz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxz_yyyyzz_0[j] = pb_x * tg_xxz_yyyyzz_0[j] + fr * tg_xxz_yyyyzz_1[j] + fl1_fx * (tg_xz_yyyyzz_0[j] - tg_xz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxz_yyyzzz_0[j] = pb_x * tg_xxz_yyyzzz_0[j] + fr * tg_xxz_yyyzzz_1[j] + fl1_fx * (tg_xz_yyyzzz_0[j] - tg_xz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxz_yyzzzz_0[j] = pb_x * tg_xxz_yyzzzz_0[j] + fr * tg_xxz_yyzzzz_1[j] + fl1_fx * (tg_xz_yyzzzz_0[j] - tg_xz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxz_yzzzzz_0[j] = pb_x * tg_xxz_yzzzzz_0[j] + fr * tg_xxz_yzzzzz_1[j] + fl1_fx * (tg_xz_yzzzzz_0[j] - tg_xz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxz_zzzzzz_0[j] = pb_x * tg_xxz_zzzzzz_0[j] + fr * tg_xxz_zzzzzz_1[j] + fl1_fx * (tg_xz_zzzzzz_0[j] - tg_xz_zzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSI_84_168(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (84,168)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tg_xyy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 84); 

                auto tg_xyy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 85); 

                auto tg_xyy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 86); 

                auto tg_xyy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 87); 

                auto tg_xyy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 88); 

                auto tg_xyy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 89); 

                auto tg_xyy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 90); 

                auto tg_xyy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 91); 

                auto tg_xyy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 92); 

                auto tg_xyy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 93); 

                auto tg_xyy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 94); 

                auto tg_xyy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 95); 

                auto tg_xyy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 96); 

                auto tg_xyy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 97); 

                auto tg_xyy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 98); 

                auto tg_xyy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 99); 

                auto tg_xyy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 100); 

                auto tg_xyy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 101); 

                auto tg_xyy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 102); 

                auto tg_xyy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 103); 

                auto tg_xyy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 104); 

                auto tg_xyy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 105); 

                auto tg_xyy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 106); 

                auto tg_xyy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 107); 

                auto tg_xyy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 108); 

                auto tg_xyy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 109); 

                auto tg_xyy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 110); 

                auto tg_xyy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 111); 

                auto tg_xyz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 112); 

                auto tg_xyz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 113); 

                auto tg_xyz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 114); 

                auto tg_xyz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 115); 

                auto tg_xyz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 116); 

                auto tg_xyz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 117); 

                auto tg_xyz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 118); 

                auto tg_xyz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 119); 

                auto tg_xyz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 120); 

                auto tg_xyz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 121); 

                auto tg_xyz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 122); 

                auto tg_xyz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 123); 

                auto tg_xyz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 124); 

                auto tg_xyz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 125); 

                auto tg_xyz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 126); 

                auto tg_xyz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 127); 

                auto tg_xyz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 128); 

                auto tg_xyz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 129); 

                auto tg_xyz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 130); 

                auto tg_xyz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 131); 

                auto tg_xyz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 132); 

                auto tg_xyz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 133); 

                auto tg_xyz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 134); 

                auto tg_xyz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 135); 

                auto tg_xyz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 136); 

                auto tg_xyz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 137); 

                auto tg_xyz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 138); 

                auto tg_xyz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 139); 

                auto tg_xzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 140); 

                auto tg_xzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 141); 

                auto tg_xzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 142); 

                auto tg_xzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 143); 

                auto tg_xzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 144); 

                auto tg_xzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 145); 

                auto tg_xzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 146); 

                auto tg_xzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 147); 

                auto tg_xzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 148); 

                auto tg_xzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 149); 

                auto tg_xzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 150); 

                auto tg_xzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 151); 

                auto tg_xzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 152); 

                auto tg_xzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 153); 

                auto tg_xzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 154); 

                auto tg_xzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 155); 

                auto tg_xzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 156); 

                auto tg_xzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 157); 

                auto tg_xzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 158); 

                auto tg_xzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 159); 

                auto tg_xzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 160); 

                auto tg_xzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 161); 

                auto tg_xzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 162); 

                auto tg_xzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 163); 

                auto tg_xzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 164); 

                auto tg_xzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 165); 

                auto tg_xzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 166); 

                auto tg_xzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 167); 

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

                // set up pointers to integrals

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

                // Batch of Integrals (84,168)

                #pragma omp simd aligned(fxn, fza, tg_xxyy_xxxxxx_0, tg_xxyy_xxxxxy_0, tg_xxyy_xxxxxz_0, \
                                         tg_xxyy_xxxxyy_0, tg_xxyy_xxxxyz_0, tg_xxyy_xxxxzz_0, tg_xxyy_xxxyyy_0, \
                                         tg_xxyy_xxxyyz_0, tg_xxyy_xxxyzz_0, tg_xxyy_xxxzzz_0, tg_xxyy_xxyyyy_0, \
                                         tg_xxyy_xxyyyz_0, tg_xxyy_xxyyzz_0, tg_xxyy_xxyzzz_0, tg_xxyy_xxzzzz_0, \
                                         tg_xxyy_xyyyyy_0, tg_xxyy_xyyyyz_0, tg_xxyy_xyyyzz_0, tg_xxyy_xyyzzz_0, \
                                         tg_xxyy_xyzzzz_0, tg_xxyy_xzzzzz_0, tg_xxyy_yyyyyy_0, tg_xxyy_yyyyyz_0, \
                                         tg_xxyy_yyyyzz_0, tg_xxyy_yyyzzz_0, tg_xxyy_yyzzzz_0, tg_xxyy_yzzzzz_0, \
                                         tg_xxyy_zzzzzz_0, tg_xxyz_xxxxxx_0, tg_xxyz_xxxxxy_0, tg_xxyz_xxxxxz_0, \
                                         tg_xxyz_xxxxyy_0, tg_xxyz_xxxxyz_0, tg_xxyz_xxxxzz_0, tg_xxyz_xxxyyy_0, \
                                         tg_xxyz_xxxyyz_0, tg_xxyz_xxxyzz_0, tg_xxyz_xxxzzz_0, tg_xxyz_xxyyyy_0, \
                                         tg_xxyz_xxyyyz_0, tg_xxyz_xxyyzz_0, tg_xxyz_xxyzzz_0, tg_xxyz_xxzzzz_0, \
                                         tg_xxyz_xyyyyy_0, tg_xxyz_xyyyyz_0, tg_xxyz_xyyyzz_0, tg_xxyz_xyyzzz_0, \
                                         tg_xxyz_xyzzzz_0, tg_xxyz_xzzzzz_0, tg_xxyz_yyyyyy_0, tg_xxyz_yyyyyz_0, \
                                         tg_xxyz_yyyyzz_0, tg_xxyz_yyyzzz_0, tg_xxyz_yyzzzz_0, tg_xxyz_yzzzzz_0, \
                                         tg_xxyz_zzzzzz_0, tg_xxzz_xxxxxx_0, tg_xxzz_xxxxxy_0, tg_xxzz_xxxxxz_0, \
                                         tg_xxzz_xxxxyy_0, tg_xxzz_xxxxyz_0, tg_xxzz_xxxxzz_0, tg_xxzz_xxxyyy_0, \
                                         tg_xxzz_xxxyyz_0, tg_xxzz_xxxyzz_0, tg_xxzz_xxxzzz_0, tg_xxzz_xxyyyy_0, \
                                         tg_xxzz_xxyyyz_0, tg_xxzz_xxyyzz_0, tg_xxzz_xxyzzz_0, tg_xxzz_xxzzzz_0, \
                                         tg_xxzz_xyyyyy_0, tg_xxzz_xyyyyz_0, tg_xxzz_xyyyzz_0, tg_xxzz_xyyzzz_0, \
                                         tg_xxzz_xyzzzz_0, tg_xxzz_xzzzzz_0, tg_xxzz_yyyyyy_0, tg_xxzz_yyyyyz_0, \
                                         tg_xxzz_yyyyzz_0, tg_xxzz_yyyzzz_0, tg_xxzz_yyzzzz_0, tg_xxzz_yzzzzz_0, \
                                         tg_xxzz_zzzzzz_0, tg_xyy_xxxxx_1, tg_xyy_xxxxxx_0, tg_xyy_xxxxxx_1, tg_xyy_xxxxxy_0, \
                                         tg_xyy_xxxxxy_1, tg_xyy_xxxxxz_0, tg_xyy_xxxxxz_1, tg_xyy_xxxxy_1, tg_xyy_xxxxyy_0, \
                                         tg_xyy_xxxxyy_1, tg_xyy_xxxxyz_0, tg_xyy_xxxxyz_1, tg_xyy_xxxxz_1, tg_xyy_xxxxzz_0, \
                                         tg_xyy_xxxxzz_1, tg_xyy_xxxyy_1, tg_xyy_xxxyyy_0, tg_xyy_xxxyyy_1, tg_xyy_xxxyyz_0, \
                                         tg_xyy_xxxyyz_1, tg_xyy_xxxyz_1, tg_xyy_xxxyzz_0, tg_xyy_xxxyzz_1, tg_xyy_xxxzz_1, \
                                         tg_xyy_xxxzzz_0, tg_xyy_xxxzzz_1, tg_xyy_xxyyy_1, tg_xyy_xxyyyy_0, tg_xyy_xxyyyy_1, \
                                         tg_xyy_xxyyyz_0, tg_xyy_xxyyyz_1, tg_xyy_xxyyz_1, tg_xyy_xxyyzz_0, tg_xyy_xxyyzz_1, \
                                         tg_xyy_xxyzz_1, tg_xyy_xxyzzz_0, tg_xyy_xxyzzz_1, tg_xyy_xxzzz_1, tg_xyy_xxzzzz_0, \
                                         tg_xyy_xxzzzz_1, tg_xyy_xyyyy_1, tg_xyy_xyyyyy_0, tg_xyy_xyyyyy_1, tg_xyy_xyyyyz_0, \
                                         tg_xyy_xyyyyz_1, tg_xyy_xyyyz_1, tg_xyy_xyyyzz_0, tg_xyy_xyyyzz_1, tg_xyy_xyyzz_1, \
                                         tg_xyy_xyyzzz_0, tg_xyy_xyyzzz_1, tg_xyy_xyzzz_1, tg_xyy_xyzzzz_0, tg_xyy_xyzzzz_1, \
                                         tg_xyy_xzzzz_1, tg_xyy_xzzzzz_0, tg_xyy_xzzzzz_1, tg_xyy_yyyyy_1, tg_xyy_yyyyyy_0, \
                                         tg_xyy_yyyyyy_1, tg_xyy_yyyyyz_0, tg_xyy_yyyyyz_1, tg_xyy_yyyyz_1, tg_xyy_yyyyzz_0, \
                                         tg_xyy_yyyyzz_1, tg_xyy_yyyzz_1, tg_xyy_yyyzzz_0, tg_xyy_yyyzzz_1, tg_xyy_yyzzz_1, \
                                         tg_xyy_yyzzzz_0, tg_xyy_yyzzzz_1, tg_xyy_yzzzz_1, tg_xyy_yzzzzz_0, tg_xyy_yzzzzz_1, \
                                         tg_xyy_zzzzz_1, tg_xyy_zzzzzz_0, tg_xyy_zzzzzz_1, tg_xyz_xxxxx_1, tg_xyz_xxxxxx_0, \
                                         tg_xyz_xxxxxx_1, tg_xyz_xxxxxy_0, tg_xyz_xxxxxy_1, tg_xyz_xxxxxz_0, tg_xyz_xxxxxz_1, \
                                         tg_xyz_xxxxy_1, tg_xyz_xxxxyy_0, tg_xyz_xxxxyy_1, tg_xyz_xxxxyz_0, tg_xyz_xxxxyz_1, \
                                         tg_xyz_xxxxz_1, tg_xyz_xxxxzz_0, tg_xyz_xxxxzz_1, tg_xyz_xxxyy_1, tg_xyz_xxxyyy_0, \
                                         tg_xyz_xxxyyy_1, tg_xyz_xxxyyz_0, tg_xyz_xxxyyz_1, tg_xyz_xxxyz_1, tg_xyz_xxxyzz_0, \
                                         tg_xyz_xxxyzz_1, tg_xyz_xxxzz_1, tg_xyz_xxxzzz_0, tg_xyz_xxxzzz_1, tg_xyz_xxyyy_1, \
                                         tg_xyz_xxyyyy_0, tg_xyz_xxyyyy_1, tg_xyz_xxyyyz_0, tg_xyz_xxyyyz_1, tg_xyz_xxyyz_1, \
                                         tg_xyz_xxyyzz_0, tg_xyz_xxyyzz_1, tg_xyz_xxyzz_1, tg_xyz_xxyzzz_0, tg_xyz_xxyzzz_1, \
                                         tg_xyz_xxzzz_1, tg_xyz_xxzzzz_0, tg_xyz_xxzzzz_1, tg_xyz_xyyyy_1, tg_xyz_xyyyyy_0, \
                                         tg_xyz_xyyyyy_1, tg_xyz_xyyyyz_0, tg_xyz_xyyyyz_1, tg_xyz_xyyyz_1, tg_xyz_xyyyzz_0, \
                                         tg_xyz_xyyyzz_1, tg_xyz_xyyzz_1, tg_xyz_xyyzzz_0, tg_xyz_xyyzzz_1, tg_xyz_xyzzz_1, \
                                         tg_xyz_xyzzzz_0, tg_xyz_xyzzzz_1, tg_xyz_xzzzz_1, tg_xyz_xzzzzz_0, tg_xyz_xzzzzz_1, \
                                         tg_xyz_yyyyy_1, tg_xyz_yyyyyy_0, tg_xyz_yyyyyy_1, tg_xyz_yyyyyz_0, tg_xyz_yyyyyz_1, \
                                         tg_xyz_yyyyz_1, tg_xyz_yyyyzz_0, tg_xyz_yyyyzz_1, tg_xyz_yyyzz_1, tg_xyz_yyyzzz_0, \
                                         tg_xyz_yyyzzz_1, tg_xyz_yyzzz_1, tg_xyz_yyzzzz_0, tg_xyz_yyzzzz_1, tg_xyz_yzzzz_1, \
                                         tg_xyz_yzzzzz_0, tg_xyz_yzzzzz_1, tg_xyz_zzzzz_1, tg_xyz_zzzzzz_0, tg_xyz_zzzzzz_1, \
                                         tg_xzz_xxxxx_1, tg_xzz_xxxxxx_0, tg_xzz_xxxxxx_1, tg_xzz_xxxxxy_0, tg_xzz_xxxxxy_1, \
                                         tg_xzz_xxxxxz_0, tg_xzz_xxxxxz_1, tg_xzz_xxxxy_1, tg_xzz_xxxxyy_0, tg_xzz_xxxxyy_1, \
                                         tg_xzz_xxxxyz_0, tg_xzz_xxxxyz_1, tg_xzz_xxxxz_1, tg_xzz_xxxxzz_0, tg_xzz_xxxxzz_1, \
                                         tg_xzz_xxxyy_1, tg_xzz_xxxyyy_0, tg_xzz_xxxyyy_1, tg_xzz_xxxyyz_0, tg_xzz_xxxyyz_1, \
                                         tg_xzz_xxxyz_1, tg_xzz_xxxyzz_0, tg_xzz_xxxyzz_1, tg_xzz_xxxzz_1, tg_xzz_xxxzzz_0, \
                                         tg_xzz_xxxzzz_1, tg_xzz_xxyyy_1, tg_xzz_xxyyyy_0, tg_xzz_xxyyyy_1, tg_xzz_xxyyyz_0, \
                                         tg_xzz_xxyyyz_1, tg_xzz_xxyyz_1, tg_xzz_xxyyzz_0, tg_xzz_xxyyzz_1, tg_xzz_xxyzz_1, \
                                         tg_xzz_xxyzzz_0, tg_xzz_xxyzzz_1, tg_xzz_xxzzz_1, tg_xzz_xxzzzz_0, tg_xzz_xxzzzz_1, \
                                         tg_xzz_xyyyy_1, tg_xzz_xyyyyy_0, tg_xzz_xyyyyy_1, tg_xzz_xyyyyz_0, tg_xzz_xyyyyz_1, \
                                         tg_xzz_xyyyz_1, tg_xzz_xyyyzz_0, tg_xzz_xyyyzz_1, tg_xzz_xyyzz_1, tg_xzz_xyyzzz_0, \
                                         tg_xzz_xyyzzz_1, tg_xzz_xyzzz_1, tg_xzz_xyzzzz_0, tg_xzz_xyzzzz_1, tg_xzz_xzzzz_1, \
                                         tg_xzz_xzzzzz_0, tg_xzz_xzzzzz_1, tg_xzz_yyyyy_1, tg_xzz_yyyyyy_0, tg_xzz_yyyyyy_1, \
                                         tg_xzz_yyyyyz_0, tg_xzz_yyyyyz_1, tg_xzz_yyyyz_1, tg_xzz_yyyyzz_0, tg_xzz_yyyyzz_1, \
                                         tg_xzz_yyyzz_1, tg_xzz_yyyzzz_0, tg_xzz_yyyzzz_1, tg_xzz_yyzzz_1, tg_xzz_yyzzzz_0, \
                                         tg_xzz_yyzzzz_1, tg_xzz_yzzzz_1, tg_xzz_yzzzzz_0, tg_xzz_yzzzzz_1, tg_xzz_zzzzz_1, \
                                         tg_xzz_zzzzzz_0, tg_xzz_zzzzzz_1, tg_yy_xxxxxx_0, tg_yy_xxxxxx_1, tg_yy_xxxxxy_0, \
                                         tg_yy_xxxxxy_1, tg_yy_xxxxxz_0, tg_yy_xxxxxz_1, tg_yy_xxxxyy_0, tg_yy_xxxxyy_1, \
                                         tg_yy_xxxxyz_0, tg_yy_xxxxyz_1, tg_yy_xxxxzz_0, tg_yy_xxxxzz_1, tg_yy_xxxyyy_0, \
                                         tg_yy_xxxyyy_1, tg_yy_xxxyyz_0, tg_yy_xxxyyz_1, tg_yy_xxxyzz_0, tg_yy_xxxyzz_1, \
                                         tg_yy_xxxzzz_0, tg_yy_xxxzzz_1, tg_yy_xxyyyy_0, tg_yy_xxyyyy_1, tg_yy_xxyyyz_0, \
                                         tg_yy_xxyyyz_1, tg_yy_xxyyzz_0, tg_yy_xxyyzz_1, tg_yy_xxyzzz_0, tg_yy_xxyzzz_1, \
                                         tg_yy_xxzzzz_0, tg_yy_xxzzzz_1, tg_yy_xyyyyy_0, tg_yy_xyyyyy_1, tg_yy_xyyyyz_0, \
                                         tg_yy_xyyyyz_1, tg_yy_xyyyzz_0, tg_yy_xyyyzz_1, tg_yy_xyyzzz_0, tg_yy_xyyzzz_1, \
                                         tg_yy_xyzzzz_0, tg_yy_xyzzzz_1, tg_yy_xzzzzz_0, tg_yy_xzzzzz_1, tg_yy_yyyyyy_0, \
                                         tg_yy_yyyyyy_1, tg_yy_yyyyyz_0, tg_yy_yyyyyz_1, tg_yy_yyyyzz_0, tg_yy_yyyyzz_1, \
                                         tg_yy_yyyzzz_0, tg_yy_yyyzzz_1, tg_yy_yyzzzz_0, tg_yy_yyzzzz_1, tg_yy_yzzzzz_0, \
                                         tg_yy_yzzzzz_1, tg_yy_zzzzzz_0, tg_yy_zzzzzz_1, tg_yz_xxxxxx_0, tg_yz_xxxxxx_1, \
                                         tg_yz_xxxxxy_0, tg_yz_xxxxxy_1, tg_yz_xxxxxz_0, tg_yz_xxxxxz_1, tg_yz_xxxxyy_0, \
                                         tg_yz_xxxxyy_1, tg_yz_xxxxyz_0, tg_yz_xxxxyz_1, tg_yz_xxxxzz_0, tg_yz_xxxxzz_1, \
                                         tg_yz_xxxyyy_0, tg_yz_xxxyyy_1, tg_yz_xxxyyz_0, tg_yz_xxxyyz_1, tg_yz_xxxyzz_0, \
                                         tg_yz_xxxyzz_1, tg_yz_xxxzzz_0, tg_yz_xxxzzz_1, tg_yz_xxyyyy_0, tg_yz_xxyyyy_1, \
                                         tg_yz_xxyyyz_0, tg_yz_xxyyyz_1, tg_yz_xxyyzz_0, tg_yz_xxyyzz_1, tg_yz_xxyzzz_0, \
                                         tg_yz_xxyzzz_1, tg_yz_xxzzzz_0, tg_yz_xxzzzz_1, tg_yz_xyyyyy_0, tg_yz_xyyyyy_1, \
                                         tg_yz_xyyyyz_0, tg_yz_xyyyyz_1, tg_yz_xyyyzz_0, tg_yz_xyyyzz_1, tg_yz_xyyzzz_0, \
                                         tg_yz_xyyzzz_1, tg_yz_xyzzzz_0, tg_yz_xyzzzz_1, tg_yz_xzzzzz_0, tg_yz_xzzzzz_1, \
                                         tg_yz_yyyyyy_0, tg_yz_yyyyyy_1, tg_yz_yyyyyz_0, tg_yz_yyyyyz_1, tg_yz_yyyyzz_0, \
                                         tg_yz_yyyyzz_1, tg_yz_yyyzzz_0, tg_yz_yyyzzz_1, tg_yz_yyzzzz_0, tg_yz_yyzzzz_1, \
                                         tg_yz_yzzzzz_0, tg_yz_yzzzzz_1, tg_yz_zzzzzz_0, tg_yz_zzzzzz_1, tg_zz_xxxxxx_0, \
                                         tg_zz_xxxxxx_1, tg_zz_xxxxxy_0, tg_zz_xxxxxy_1, tg_zz_xxxxxz_0, tg_zz_xxxxxz_1, \
                                         tg_zz_xxxxyy_0, tg_zz_xxxxyy_1, tg_zz_xxxxyz_0, tg_zz_xxxxyz_1, tg_zz_xxxxzz_0, \
                                         tg_zz_xxxxzz_1, tg_zz_xxxyyy_0, tg_zz_xxxyyy_1, tg_zz_xxxyyz_0, tg_zz_xxxyyz_1, \
                                         tg_zz_xxxyzz_0, tg_zz_xxxyzz_1, tg_zz_xxxzzz_0, tg_zz_xxxzzz_1, tg_zz_xxyyyy_0, \
                                         tg_zz_xxyyyy_1, tg_zz_xxyyyz_0, tg_zz_xxyyyz_1, tg_zz_xxyyzz_0, tg_zz_xxyyzz_1, \
                                         tg_zz_xxyzzz_0, tg_zz_xxyzzz_1, tg_zz_xxzzzz_0, tg_zz_xxzzzz_1, tg_zz_xyyyyy_0, \
                                         tg_zz_xyyyyy_1, tg_zz_xyyyyz_0, tg_zz_xyyyyz_1, tg_zz_xyyyzz_0, tg_zz_xyyyzz_1, \
                                         tg_zz_xyyzzz_0, tg_zz_xyyzzz_1, tg_zz_xyzzzz_0, tg_zz_xyzzzz_1, tg_zz_xzzzzz_0, \
                                         tg_zz_xzzzzz_1, tg_zz_yyyyyy_0, tg_zz_yyyyyy_1, tg_zz_yyyyyz_0, tg_zz_yyyyyz_1, \
                                         tg_zz_yyyyzz_0, tg_zz_yyyyzz_1, tg_zz_yyyzzz_0, tg_zz_yyyzzz_1, tg_zz_yyzzzz_0, \
                                         tg_zz_yyzzzz_1, tg_zz_yzzzzz_0, tg_zz_yzzzzz_1, tg_zz_zzzzzz_0, tg_zz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyy_xxxxxx_0[j] = pb_x * tg_xyy_xxxxxx_0[j] + fr * tg_xyy_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxx_0[j] - tg_yy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyy_xxxxx_1[j];

                    tg_xxyy_xxxxxy_0[j] = pb_x * tg_xyy_xxxxxy_0[j] + fr * tg_xyy_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxy_0[j] - tg_yy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyy_xxxxy_1[j];

                    tg_xxyy_xxxxxz_0[j] = pb_x * tg_xyy_xxxxxz_0[j] + fr * tg_xyy_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxxz_0[j] - tg_yy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyy_xxxxz_1[j];

                    tg_xxyy_xxxxyy_0[j] = pb_x * tg_xyy_xxxxyy_0[j] + fr * tg_xyy_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxyy_0[j] - tg_yy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxyy_1[j];

                    tg_xxyy_xxxxyz_0[j] = pb_x * tg_xyy_xxxxyz_0[j] + fr * tg_xyy_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxyz_0[j] - tg_yy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxyz_1[j];

                    tg_xxyy_xxxxzz_0[j] = pb_x * tg_xyy_xxxxzz_0[j] + fr * tg_xyy_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxxzz_0[j] - tg_yy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyy_xxxzz_1[j];

                    tg_xxyy_xxxyyy_0[j] = pb_x * tg_xyy_xxxyyy_0[j] + fr * tg_xyy_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyyy_0[j] - tg_yy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyyy_1[j];

                    tg_xxyy_xxxyyz_0[j] = pb_x * tg_xyy_xxxyyz_0[j] + fr * tg_xyy_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyyz_0[j] - tg_yy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyyz_1[j];

                    tg_xxyy_xxxyzz_0[j] = pb_x * tg_xyy_xxxyzz_0[j] + fr * tg_xyy_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxyzz_0[j] - tg_yy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxyzz_1[j];

                    tg_xxyy_xxxzzz_0[j] = pb_x * tg_xyy_xxxzzz_0[j] + fr * tg_xyy_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxxzzz_0[j] - tg_yy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xxzzz_1[j];

                    tg_xxyy_xxyyyy_0[j] = pb_x * tg_xyy_xxyyyy_0[j] + fr * tg_xyy_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyyy_0[j] - tg_yy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyyy_1[j];

                    tg_xxyy_xxyyyz_0[j] = pb_x * tg_xyy_xxyyyz_0[j] + fr * tg_xyy_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyyz_0[j] - tg_yy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyyz_1[j];

                    tg_xxyy_xxyyzz_0[j] = pb_x * tg_xyy_xxyyzz_0[j] + fr * tg_xyy_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyyzz_0[j] - tg_yy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyyzz_1[j];

                    tg_xxyy_xxyzzz_0[j] = pb_x * tg_xyy_xxyzzz_0[j] + fr * tg_xyy_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxyzzz_0[j] - tg_yy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xyzzz_1[j];

                    tg_xxyy_xxzzzz_0[j] = pb_x * tg_xyy_xxzzzz_0[j] + fr * tg_xyy_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xxzzzz_0[j] - tg_yy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xzzzz_1[j];

                    tg_xxyy_xyyyyy_0[j] = pb_x * tg_xyy_xyyyyy_0[j] + fr * tg_xyy_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyyy_0[j] - tg_yy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyyy_1[j];

                    tg_xxyy_xyyyyz_0[j] = pb_x * tg_xyy_xyyyyz_0[j] + fr * tg_xyy_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyyz_0[j] - tg_yy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyyz_1[j];

                    tg_xxyy_xyyyzz_0[j] = pb_x * tg_xyy_xyyyzz_0[j] + fr * tg_xyy_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyyzz_0[j] - tg_yy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyyzz_1[j];

                    tg_xxyy_xyyzzz_0[j] = pb_x * tg_xyy_xyyzzz_0[j] + fr * tg_xyy_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyyzzz_0[j] - tg_yy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yyzzz_1[j];

                    tg_xxyy_xyzzzz_0[j] = pb_x * tg_xyy_xyzzzz_0[j] + fr * tg_xyy_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xyzzzz_0[j] - tg_yy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yzzzz_1[j];

                    tg_xxyy_xzzzzz_0[j] = pb_x * tg_xyy_xzzzzz_0[j] + fr * tg_xyy_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_xzzzzz_0[j] - tg_yy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_zzzzz_1[j];

                    tg_xxyy_yyyyyy_0[j] = pb_x * tg_xyy_yyyyyy_0[j] + fr * tg_xyy_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyyy_0[j] - tg_yy_yyyyyy_1[j] * fl1_fza);

                    tg_xxyy_yyyyyz_0[j] = pb_x * tg_xyy_yyyyyz_0[j] + fr * tg_xyy_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyyz_0[j] - tg_yy_yyyyyz_1[j] * fl1_fza);

                    tg_xxyy_yyyyzz_0[j] = pb_x * tg_xyy_yyyyzz_0[j] + fr * tg_xyy_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyyzz_0[j] - tg_yy_yyyyzz_1[j] * fl1_fza);

                    tg_xxyy_yyyzzz_0[j] = pb_x * tg_xyy_yyyzzz_0[j] + fr * tg_xyy_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyyzzz_0[j] - tg_yy_yyyzzz_1[j] * fl1_fza);

                    tg_xxyy_yyzzzz_0[j] = pb_x * tg_xyy_yyzzzz_0[j] + fr * tg_xyy_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yyzzzz_0[j] - tg_yy_yyzzzz_1[j] * fl1_fza);

                    tg_xxyy_yzzzzz_0[j] = pb_x * tg_xyy_yzzzzz_0[j] + fr * tg_xyy_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_yzzzzz_0[j] - tg_yy_yzzzzz_1[j] * fl1_fza);

                    tg_xxyy_zzzzzz_0[j] = pb_x * tg_xyy_zzzzzz_0[j] + fr * tg_xyy_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yy_zzzzzz_0[j] - tg_yy_zzzzzz_1[j] * fl1_fza);

                    tg_xxyz_xxxxxx_0[j] = pb_x * tg_xyz_xxxxxx_0[j] + fr * tg_xyz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxx_0[j] - tg_yz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyz_xxxxx_1[j];

                    tg_xxyz_xxxxxy_0[j] = pb_x * tg_xyz_xxxxxy_0[j] + fr * tg_xyz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxy_0[j] - tg_yz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyz_xxxxy_1[j];

                    tg_xxyz_xxxxxz_0[j] = pb_x * tg_xyz_xxxxxz_0[j] + fr * tg_xyz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxxz_0[j] - tg_yz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyz_xxxxz_1[j];

                    tg_xxyz_xxxxyy_0[j] = pb_x * tg_xyz_xxxxyy_0[j] + fr * tg_xyz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxyy_0[j] - tg_yz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxyy_1[j];

                    tg_xxyz_xxxxyz_0[j] = pb_x * tg_xyz_xxxxyz_0[j] + fr * tg_xyz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxyz_0[j] - tg_yz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxyz_1[j];

                    tg_xxyz_xxxxzz_0[j] = pb_x * tg_xyz_xxxxzz_0[j] + fr * tg_xyz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxxzz_0[j] - tg_yz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyz_xxxzz_1[j];

                    tg_xxyz_xxxyyy_0[j] = pb_x * tg_xyz_xxxyyy_0[j] + fr * tg_xyz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyyy_0[j] - tg_yz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyyy_1[j];

                    tg_xxyz_xxxyyz_0[j] = pb_x * tg_xyz_xxxyyz_0[j] + fr * tg_xyz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyyz_0[j] - tg_yz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyyz_1[j];

                    tg_xxyz_xxxyzz_0[j] = pb_x * tg_xyz_xxxyzz_0[j] + fr * tg_xyz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxyzz_0[j] - tg_yz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxyzz_1[j];

                    tg_xxyz_xxxzzz_0[j] = pb_x * tg_xyz_xxxzzz_0[j] + fr * tg_xyz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxxzzz_0[j] - tg_yz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xxzzz_1[j];

                    tg_xxyz_xxyyyy_0[j] = pb_x * tg_xyz_xxyyyy_0[j] + fr * tg_xyz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyyy_0[j] - tg_yz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyyy_1[j];

                    tg_xxyz_xxyyyz_0[j] = pb_x * tg_xyz_xxyyyz_0[j] + fr * tg_xyz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyyz_0[j] - tg_yz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyyz_1[j];

                    tg_xxyz_xxyyzz_0[j] = pb_x * tg_xyz_xxyyzz_0[j] + fr * tg_xyz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyyzz_0[j] - tg_yz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyyzz_1[j];

                    tg_xxyz_xxyzzz_0[j] = pb_x * tg_xyz_xxyzzz_0[j] + fr * tg_xyz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxyzzz_0[j] - tg_yz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xyzzz_1[j];

                    tg_xxyz_xxzzzz_0[j] = pb_x * tg_xyz_xxzzzz_0[j] + fr * tg_xyz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xxzzzz_0[j] - tg_yz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xzzzz_1[j];

                    tg_xxyz_xyyyyy_0[j] = pb_x * tg_xyz_xyyyyy_0[j] + fr * tg_xyz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyyy_0[j] - tg_yz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyyy_1[j];

                    tg_xxyz_xyyyyz_0[j] = pb_x * tg_xyz_xyyyyz_0[j] + fr * tg_xyz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyyz_0[j] - tg_yz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyyz_1[j];

                    tg_xxyz_xyyyzz_0[j] = pb_x * tg_xyz_xyyyzz_0[j] + fr * tg_xyz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyyzz_0[j] - tg_yz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyyzz_1[j];

                    tg_xxyz_xyyzzz_0[j] = pb_x * tg_xyz_xyyzzz_0[j] + fr * tg_xyz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyyzzz_0[j] - tg_yz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yyzzz_1[j];

                    tg_xxyz_xyzzzz_0[j] = pb_x * tg_xyz_xyzzzz_0[j] + fr * tg_xyz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xyzzzz_0[j] - tg_yz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yzzzz_1[j];

                    tg_xxyz_xzzzzz_0[j] = pb_x * tg_xyz_xzzzzz_0[j] + fr * tg_xyz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_xzzzzz_0[j] - tg_yz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_zzzzz_1[j];

                    tg_xxyz_yyyyyy_0[j] = pb_x * tg_xyz_yyyyyy_0[j] + fr * tg_xyz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyyy_0[j] - tg_yz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyz_yyyyyz_0[j] = pb_x * tg_xyz_yyyyyz_0[j] + fr * tg_xyz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyyz_0[j] - tg_yz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyz_yyyyzz_0[j] = pb_x * tg_xyz_yyyyzz_0[j] + fr * tg_xyz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyyzz_0[j] - tg_yz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyz_yyyzzz_0[j] = pb_x * tg_xyz_yyyzzz_0[j] + fr * tg_xyz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyyzzz_0[j] - tg_yz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyz_yyzzzz_0[j] = pb_x * tg_xyz_yyzzzz_0[j] + fr * tg_xyz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yyzzzz_0[j] - tg_yz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyz_yzzzzz_0[j] = pb_x * tg_xyz_yzzzzz_0[j] + fr * tg_xyz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_yzzzzz_0[j] - tg_yz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyz_zzzzzz_0[j] = pb_x * tg_xyz_zzzzzz_0[j] + fr * tg_xyz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yz_zzzzzz_0[j] - tg_yz_zzzzzz_1[j] * fl1_fza);

                    tg_xxzz_xxxxxx_0[j] = pb_x * tg_xzz_xxxxxx_0[j] + fr * tg_xzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxx_0[j] - tg_zz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzz_xxxxx_1[j];

                    tg_xxzz_xxxxxy_0[j] = pb_x * tg_xzz_xxxxxy_0[j] + fr * tg_xzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxy_0[j] - tg_zz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzz_xxxxy_1[j];

                    tg_xxzz_xxxxxz_0[j] = pb_x * tg_xzz_xxxxxz_0[j] + fr * tg_xzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxz_0[j] - tg_zz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzz_xxxxz_1[j];

                    tg_xxzz_xxxxyy_0[j] = pb_x * tg_xzz_xxxxyy_0[j] + fr * tg_xzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyy_0[j] - tg_zz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxyy_1[j];

                    tg_xxzz_xxxxyz_0[j] = pb_x * tg_xzz_xxxxyz_0[j] + fr * tg_xzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyz_0[j] - tg_zz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxyz_1[j];

                    tg_xxzz_xxxxzz_0[j] = pb_x * tg_xzz_xxxxzz_0[j] + fr * tg_xzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxzz_0[j] - tg_zz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzz_xxxzz_1[j];

                    tg_xxzz_xxxyyy_0[j] = pb_x * tg_xzz_xxxyyy_0[j] + fr * tg_xzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyy_0[j] - tg_zz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyyy_1[j];

                    tg_xxzz_xxxyyz_0[j] = pb_x * tg_xzz_xxxyyz_0[j] + fr * tg_xzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyz_0[j] - tg_zz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyyz_1[j];

                    tg_xxzz_xxxyzz_0[j] = pb_x * tg_xzz_xxxyzz_0[j] + fr * tg_xzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyzz_0[j] - tg_zz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxyzz_1[j];

                    tg_xxzz_xxxzzz_0[j] = pb_x * tg_xzz_xxxzzz_0[j] + fr * tg_xzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxzzz_0[j] - tg_zz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xxzzz_1[j];

                    tg_xxzz_xxyyyy_0[j] = pb_x * tg_xzz_xxyyyy_0[j] + fr * tg_xzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyy_0[j] - tg_zz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyyy_1[j];

                    tg_xxzz_xxyyyz_0[j] = pb_x * tg_xzz_xxyyyz_0[j] + fr * tg_xzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyz_0[j] - tg_zz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyyz_1[j];

                    tg_xxzz_xxyyzz_0[j] = pb_x * tg_xzz_xxyyzz_0[j] + fr * tg_xzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyzz_0[j] - tg_zz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyyzz_1[j];

                    tg_xxzz_xxyzzz_0[j] = pb_x * tg_xzz_xxyzzz_0[j] + fr * tg_xzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyzzz_0[j] - tg_zz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xyzzz_1[j];

                    tg_xxzz_xxzzzz_0[j] = pb_x * tg_xzz_xxzzzz_0[j] + fr * tg_xzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxzzzz_0[j] - tg_zz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xzzzz_1[j];

                    tg_xxzz_xyyyyy_0[j] = pb_x * tg_xzz_xyyyyy_0[j] + fr * tg_xzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyy_0[j] - tg_zz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyyy_1[j];

                    tg_xxzz_xyyyyz_0[j] = pb_x * tg_xzz_xyyyyz_0[j] + fr * tg_xzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyz_0[j] - tg_zz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyyz_1[j];

                    tg_xxzz_xyyyzz_0[j] = pb_x * tg_xzz_xyyyzz_0[j] + fr * tg_xzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyzz_0[j] - tg_zz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyyzz_1[j];

                    tg_xxzz_xyyzzz_0[j] = pb_x * tg_xzz_xyyzzz_0[j] + fr * tg_xzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyzzz_0[j] - tg_zz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yyzzz_1[j];

                    tg_xxzz_xyzzzz_0[j] = pb_x * tg_xzz_xyzzzz_0[j] + fr * tg_xzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyzzzz_0[j] - tg_zz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yzzzz_1[j];

                    tg_xxzz_xzzzzz_0[j] = pb_x * tg_xzz_xzzzzz_0[j] + fr * tg_xzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzzzzz_0[j] - tg_zz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_zzzzz_1[j];

                    tg_xxzz_yyyyyy_0[j] = pb_x * tg_xzz_yyyyyy_0[j] + fr * tg_xzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyy_0[j] - tg_zz_yyyyyy_1[j] * fl1_fza);

                    tg_xxzz_yyyyyz_0[j] = pb_x * tg_xzz_yyyyyz_0[j] + fr * tg_xzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyz_0[j] - tg_zz_yyyyyz_1[j] * fl1_fza);

                    tg_xxzz_yyyyzz_0[j] = pb_x * tg_xzz_yyyyzz_0[j] + fr * tg_xzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyzz_0[j] - tg_zz_yyyyzz_1[j] * fl1_fza);

                    tg_xxzz_yyyzzz_0[j] = pb_x * tg_xzz_yyyzzz_0[j] + fr * tg_xzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyzzz_0[j] - tg_zz_yyyzzz_1[j] * fl1_fza);

                    tg_xxzz_yyzzzz_0[j] = pb_x * tg_xzz_yyzzzz_0[j] + fr * tg_xzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyzzzz_0[j] - tg_zz_yyzzzz_1[j] * fl1_fza);

                    tg_xxzz_yzzzzz_0[j] = pb_x * tg_xzz_yzzzzz_0[j] + fr * tg_xzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzzzzz_0[j] - tg_zz_yzzzzz_1[j] * fl1_fza);

                    tg_xxzz_zzzzzz_0[j] = pb_x * tg_xzz_zzzzzz_0[j] + fr * tg_xzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzzzzz_0[j] - tg_zz_zzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSI_168_252(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (168,252)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
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

                auto tg_yyy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 168); 

                auto tg_yyy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 169); 

                auto tg_yyy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 170); 

                auto tg_yyy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 171); 

                auto tg_yyy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 172); 

                auto tg_yyy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 173); 

                auto tg_yyy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 174); 

                auto tg_yyy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 175); 

                auto tg_yyy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 176); 

                auto tg_yyy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 177); 

                auto tg_yyy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 178); 

                auto tg_yyy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 179); 

                auto tg_yyy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 180); 

                auto tg_yyy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 181); 

                auto tg_yyy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 182); 

                auto tg_yyy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 183); 

                auto tg_yyy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 184); 

                auto tg_yyy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 185); 

                auto tg_yyy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 186); 

                auto tg_yyy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 187); 

                auto tg_yyy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 188); 

                auto tg_yyy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 189); 

                auto tg_yyy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 190); 

                auto tg_yyy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 191); 

                auto tg_yyy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 192); 

                auto tg_yyy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 193); 

                auto tg_yyy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 194); 

                auto tg_yyy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 195); 

                auto tg_yyz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 196); 

                auto tg_yyz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 197); 

                auto tg_yyz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 198); 

                auto tg_yyz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 199); 

                auto tg_yyz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 200); 

                auto tg_yyz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 201); 

                auto tg_yyz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 202); 

                auto tg_yyz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 203); 

                auto tg_yyz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 204); 

                auto tg_yyz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 205); 

                auto tg_yyz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 206); 

                auto tg_yyz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 207); 

                auto tg_yyz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 208); 

                auto tg_yyz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 209); 

                auto tg_yyz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 210); 

                auto tg_yyz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 211); 

                auto tg_yyz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 212); 

                auto tg_yyz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 213); 

                auto tg_yyz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 214); 

                auto tg_yyz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 215); 

                auto tg_yyz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 216); 

                auto tg_yyz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 217); 

                auto tg_yyz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 218); 

                auto tg_yyz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 219); 

                auto tg_yyz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 220); 

                auto tg_yyz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 221); 

                auto tg_yyz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 222); 

                auto tg_yyz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 223); 

                auto tg_yzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 224); 

                auto tg_yzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 225); 

                auto tg_yzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 226); 

                auto tg_yzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 227); 

                auto tg_yzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 228); 

                auto tg_yzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 229); 

                auto tg_yzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 230); 

                auto tg_yzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 231); 

                auto tg_yzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 232); 

                auto tg_yzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 233); 

                auto tg_yzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 234); 

                auto tg_yzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 235); 

                auto tg_yzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 236); 

                auto tg_yzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 237); 

                auto tg_yzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 238); 

                auto tg_yzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 239); 

                auto tg_yzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 240); 

                auto tg_yzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 241); 

                auto tg_yzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 242); 

                auto tg_yzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 243); 

                auto tg_yzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 244); 

                auto tg_yzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 245); 

                auto tg_yzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 246); 

                auto tg_yzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 247); 

                auto tg_yzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 248); 

                auto tg_yzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 249); 

                auto tg_yzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 250); 

                auto tg_yzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 251); 

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

                // set up pointers to integrals

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

                // Batch of Integrals (168,252)

                #pragma omp simd aligned(fxn, tg_xyyy_xxxxxx_0, tg_xyyy_xxxxxy_0, tg_xyyy_xxxxxz_0, \
                                         tg_xyyy_xxxxyy_0, tg_xyyy_xxxxyz_0, tg_xyyy_xxxxzz_0, tg_xyyy_xxxyyy_0, \
                                         tg_xyyy_xxxyyz_0, tg_xyyy_xxxyzz_0, tg_xyyy_xxxzzz_0, tg_xyyy_xxyyyy_0, \
                                         tg_xyyy_xxyyyz_0, tg_xyyy_xxyyzz_0, tg_xyyy_xxyzzz_0, tg_xyyy_xxzzzz_0, \
                                         tg_xyyy_xyyyyy_0, tg_xyyy_xyyyyz_0, tg_xyyy_xyyyzz_0, tg_xyyy_xyyzzz_0, \
                                         tg_xyyy_xyzzzz_0, tg_xyyy_xzzzzz_0, tg_xyyy_yyyyyy_0, tg_xyyy_yyyyyz_0, \
                                         tg_xyyy_yyyyzz_0, tg_xyyy_yyyzzz_0, tg_xyyy_yyzzzz_0, tg_xyyy_yzzzzz_0, \
                                         tg_xyyy_zzzzzz_0, tg_xyyz_xxxxxx_0, tg_xyyz_xxxxxy_0, tg_xyyz_xxxxxz_0, \
                                         tg_xyyz_xxxxyy_0, tg_xyyz_xxxxyz_0, tg_xyyz_xxxxzz_0, tg_xyyz_xxxyyy_0, \
                                         tg_xyyz_xxxyyz_0, tg_xyyz_xxxyzz_0, tg_xyyz_xxxzzz_0, tg_xyyz_xxyyyy_0, \
                                         tg_xyyz_xxyyyz_0, tg_xyyz_xxyyzz_0, tg_xyyz_xxyzzz_0, tg_xyyz_xxzzzz_0, \
                                         tg_xyyz_xyyyyy_0, tg_xyyz_xyyyyz_0, tg_xyyz_xyyyzz_0, tg_xyyz_xyyzzz_0, \
                                         tg_xyyz_xyzzzz_0, tg_xyyz_xzzzzz_0, tg_xyyz_yyyyyy_0, tg_xyyz_yyyyyz_0, \
                                         tg_xyyz_yyyyzz_0, tg_xyyz_yyyzzz_0, tg_xyyz_yyzzzz_0, tg_xyyz_yzzzzz_0, \
                                         tg_xyyz_zzzzzz_0, tg_xyzz_xxxxxx_0, tg_xyzz_xxxxxy_0, tg_xyzz_xxxxxz_0, \
                                         tg_xyzz_xxxxyy_0, tg_xyzz_xxxxyz_0, tg_xyzz_xxxxzz_0, tg_xyzz_xxxyyy_0, \
                                         tg_xyzz_xxxyyz_0, tg_xyzz_xxxyzz_0, tg_xyzz_xxxzzz_0, tg_xyzz_xxyyyy_0, \
                                         tg_xyzz_xxyyyz_0, tg_xyzz_xxyyzz_0, tg_xyzz_xxyzzz_0, tg_xyzz_xxzzzz_0, \
                                         tg_xyzz_xyyyyy_0, tg_xyzz_xyyyyz_0, tg_xyzz_xyyyzz_0, tg_xyzz_xyyzzz_0, \
                                         tg_xyzz_xyzzzz_0, tg_xyzz_xzzzzz_0, tg_xyzz_yyyyyy_0, tg_xyzz_yyyyyz_0, \
                                         tg_xyzz_yyyyzz_0, tg_xyzz_yyyzzz_0, tg_xyzz_yyzzzz_0, tg_xyzz_yzzzzz_0, \
                                         tg_xyzz_zzzzzz_0, tg_yyy_xxxxx_1, tg_yyy_xxxxxx_0, tg_yyy_xxxxxx_1, tg_yyy_xxxxxy_0, \
                                         tg_yyy_xxxxxy_1, tg_yyy_xxxxxz_0, tg_yyy_xxxxxz_1, tg_yyy_xxxxy_1, tg_yyy_xxxxyy_0, \
                                         tg_yyy_xxxxyy_1, tg_yyy_xxxxyz_0, tg_yyy_xxxxyz_1, tg_yyy_xxxxz_1, tg_yyy_xxxxzz_0, \
                                         tg_yyy_xxxxzz_1, tg_yyy_xxxyy_1, tg_yyy_xxxyyy_0, tg_yyy_xxxyyy_1, tg_yyy_xxxyyz_0, \
                                         tg_yyy_xxxyyz_1, tg_yyy_xxxyz_1, tg_yyy_xxxyzz_0, tg_yyy_xxxyzz_1, tg_yyy_xxxzz_1, \
                                         tg_yyy_xxxzzz_0, tg_yyy_xxxzzz_1, tg_yyy_xxyyy_1, tg_yyy_xxyyyy_0, tg_yyy_xxyyyy_1, \
                                         tg_yyy_xxyyyz_0, tg_yyy_xxyyyz_1, tg_yyy_xxyyz_1, tg_yyy_xxyyzz_0, tg_yyy_xxyyzz_1, \
                                         tg_yyy_xxyzz_1, tg_yyy_xxyzzz_0, tg_yyy_xxyzzz_1, tg_yyy_xxzzz_1, tg_yyy_xxzzzz_0, \
                                         tg_yyy_xxzzzz_1, tg_yyy_xyyyy_1, tg_yyy_xyyyyy_0, tg_yyy_xyyyyy_1, tg_yyy_xyyyyz_0, \
                                         tg_yyy_xyyyyz_1, tg_yyy_xyyyz_1, tg_yyy_xyyyzz_0, tg_yyy_xyyyzz_1, tg_yyy_xyyzz_1, \
                                         tg_yyy_xyyzzz_0, tg_yyy_xyyzzz_1, tg_yyy_xyzzz_1, tg_yyy_xyzzzz_0, tg_yyy_xyzzzz_1, \
                                         tg_yyy_xzzzz_1, tg_yyy_xzzzzz_0, tg_yyy_xzzzzz_1, tg_yyy_yyyyy_1, tg_yyy_yyyyyy_0, \
                                         tg_yyy_yyyyyy_1, tg_yyy_yyyyyz_0, tg_yyy_yyyyyz_1, tg_yyy_yyyyz_1, tg_yyy_yyyyzz_0, \
                                         tg_yyy_yyyyzz_1, tg_yyy_yyyzz_1, tg_yyy_yyyzzz_0, tg_yyy_yyyzzz_1, tg_yyy_yyzzz_1, \
                                         tg_yyy_yyzzzz_0, tg_yyy_yyzzzz_1, tg_yyy_yzzzz_1, tg_yyy_yzzzzz_0, tg_yyy_yzzzzz_1, \
                                         tg_yyy_zzzzz_1, tg_yyy_zzzzzz_0, tg_yyy_zzzzzz_1, tg_yyz_xxxxx_1, tg_yyz_xxxxxx_0, \
                                         tg_yyz_xxxxxx_1, tg_yyz_xxxxxy_0, tg_yyz_xxxxxy_1, tg_yyz_xxxxxz_0, tg_yyz_xxxxxz_1, \
                                         tg_yyz_xxxxy_1, tg_yyz_xxxxyy_0, tg_yyz_xxxxyy_1, tg_yyz_xxxxyz_0, tg_yyz_xxxxyz_1, \
                                         tg_yyz_xxxxz_1, tg_yyz_xxxxzz_0, tg_yyz_xxxxzz_1, tg_yyz_xxxyy_1, tg_yyz_xxxyyy_0, \
                                         tg_yyz_xxxyyy_1, tg_yyz_xxxyyz_0, tg_yyz_xxxyyz_1, tg_yyz_xxxyz_1, tg_yyz_xxxyzz_0, \
                                         tg_yyz_xxxyzz_1, tg_yyz_xxxzz_1, tg_yyz_xxxzzz_0, tg_yyz_xxxzzz_1, tg_yyz_xxyyy_1, \
                                         tg_yyz_xxyyyy_0, tg_yyz_xxyyyy_1, tg_yyz_xxyyyz_0, tg_yyz_xxyyyz_1, tg_yyz_xxyyz_1, \
                                         tg_yyz_xxyyzz_0, tg_yyz_xxyyzz_1, tg_yyz_xxyzz_1, tg_yyz_xxyzzz_0, tg_yyz_xxyzzz_1, \
                                         tg_yyz_xxzzz_1, tg_yyz_xxzzzz_0, tg_yyz_xxzzzz_1, tg_yyz_xyyyy_1, tg_yyz_xyyyyy_0, \
                                         tg_yyz_xyyyyy_1, tg_yyz_xyyyyz_0, tg_yyz_xyyyyz_1, tg_yyz_xyyyz_1, tg_yyz_xyyyzz_0, \
                                         tg_yyz_xyyyzz_1, tg_yyz_xyyzz_1, tg_yyz_xyyzzz_0, tg_yyz_xyyzzz_1, tg_yyz_xyzzz_1, \
                                         tg_yyz_xyzzzz_0, tg_yyz_xyzzzz_1, tg_yyz_xzzzz_1, tg_yyz_xzzzzz_0, tg_yyz_xzzzzz_1, \
                                         tg_yyz_yyyyy_1, tg_yyz_yyyyyy_0, tg_yyz_yyyyyy_1, tg_yyz_yyyyyz_0, tg_yyz_yyyyyz_1, \
                                         tg_yyz_yyyyz_1, tg_yyz_yyyyzz_0, tg_yyz_yyyyzz_1, tg_yyz_yyyzz_1, tg_yyz_yyyzzz_0, \
                                         tg_yyz_yyyzzz_1, tg_yyz_yyzzz_1, tg_yyz_yyzzzz_0, tg_yyz_yyzzzz_1, tg_yyz_yzzzz_1, \
                                         tg_yyz_yzzzzz_0, tg_yyz_yzzzzz_1, tg_yyz_zzzzz_1, tg_yyz_zzzzzz_0, tg_yyz_zzzzzz_1, \
                                         tg_yzz_xxxxx_1, tg_yzz_xxxxxx_0, tg_yzz_xxxxxx_1, tg_yzz_xxxxxy_0, tg_yzz_xxxxxy_1, \
                                         tg_yzz_xxxxxz_0, tg_yzz_xxxxxz_1, tg_yzz_xxxxy_1, tg_yzz_xxxxyy_0, tg_yzz_xxxxyy_1, \
                                         tg_yzz_xxxxyz_0, tg_yzz_xxxxyz_1, tg_yzz_xxxxz_1, tg_yzz_xxxxzz_0, tg_yzz_xxxxzz_1, \
                                         tg_yzz_xxxyy_1, tg_yzz_xxxyyy_0, tg_yzz_xxxyyy_1, tg_yzz_xxxyyz_0, tg_yzz_xxxyyz_1, \
                                         tg_yzz_xxxyz_1, tg_yzz_xxxyzz_0, tg_yzz_xxxyzz_1, tg_yzz_xxxzz_1, tg_yzz_xxxzzz_0, \
                                         tg_yzz_xxxzzz_1, tg_yzz_xxyyy_1, tg_yzz_xxyyyy_0, tg_yzz_xxyyyy_1, tg_yzz_xxyyyz_0, \
                                         tg_yzz_xxyyyz_1, tg_yzz_xxyyz_1, tg_yzz_xxyyzz_0, tg_yzz_xxyyzz_1, tg_yzz_xxyzz_1, \
                                         tg_yzz_xxyzzz_0, tg_yzz_xxyzzz_1, tg_yzz_xxzzz_1, tg_yzz_xxzzzz_0, tg_yzz_xxzzzz_1, \
                                         tg_yzz_xyyyy_1, tg_yzz_xyyyyy_0, tg_yzz_xyyyyy_1, tg_yzz_xyyyyz_0, tg_yzz_xyyyyz_1, \
                                         tg_yzz_xyyyz_1, tg_yzz_xyyyzz_0, tg_yzz_xyyyzz_1, tg_yzz_xyyzz_1, tg_yzz_xyyzzz_0, \
                                         tg_yzz_xyyzzz_1, tg_yzz_xyzzz_1, tg_yzz_xyzzzz_0, tg_yzz_xyzzzz_1, tg_yzz_xzzzz_1, \
                                         tg_yzz_xzzzzz_0, tg_yzz_xzzzzz_1, tg_yzz_yyyyy_1, tg_yzz_yyyyyy_0, tg_yzz_yyyyyy_1, \
                                         tg_yzz_yyyyyz_0, tg_yzz_yyyyyz_1, tg_yzz_yyyyz_1, tg_yzz_yyyyzz_0, tg_yzz_yyyyzz_1, \
                                         tg_yzz_yyyzz_1, tg_yzz_yyyzzz_0, tg_yzz_yyyzzz_1, tg_yzz_yyzzz_1, tg_yzz_yyzzzz_0, \
                                         tg_yzz_yyzzzz_1, tg_yzz_yzzzz_1, tg_yzz_yzzzzz_0, tg_yzz_yzzzzz_1, tg_yzz_zzzzz_1, \
                                         tg_yzz_zzzzzz_0, tg_yzz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyy_xxxxxx_0[j] = pb_x * tg_yyy_xxxxxx_0[j] + fr * tg_yyy_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyy_xxxxx_1[j];

                    tg_xyyy_xxxxxy_0[j] = pb_x * tg_yyy_xxxxxy_0[j] + fr * tg_yyy_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxy_1[j];

                    tg_xyyy_xxxxxz_0[j] = pb_x * tg_yyy_xxxxxz_0[j] + fr * tg_yyy_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyy_xxxxz_1[j];

                    tg_xyyy_xxxxyy_0[j] = pb_x * tg_yyy_xxxxyy_0[j] + fr * tg_yyy_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyy_1[j];

                    tg_xyyy_xxxxyz_0[j] = pb_x * tg_yyy_xxxxyz_0[j] + fr * tg_yyy_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxyz_1[j];

                    tg_xyyy_xxxxzz_0[j] = pb_x * tg_yyy_xxxxzz_0[j] + fr * tg_yyy_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyy_xxxzz_1[j];

                    tg_xyyy_xxxyyy_0[j] = pb_x * tg_yyy_xxxyyy_0[j] + fr * tg_yyy_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyy_1[j];

                    tg_xyyy_xxxyyz_0[j] = pb_x * tg_yyy_xxxyyz_0[j] + fr * tg_yyy_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyyz_1[j];

                    tg_xyyy_xxxyzz_0[j] = pb_x * tg_yyy_xxxyzz_0[j] + fr * tg_yyy_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxyzz_1[j];

                    tg_xyyy_xxxzzz_0[j] = pb_x * tg_yyy_xxxzzz_0[j] + fr * tg_yyy_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxzzz_1[j];

                    tg_xyyy_xxyyyy_0[j] = pb_x * tg_yyy_xxyyyy_0[j] + fr * tg_yyy_xxyyyy_1[j] + fl1_fxn * tg_yyy_xyyyy_1[j];

                    tg_xyyy_xxyyyz_0[j] = pb_x * tg_yyy_xxyyyz_0[j] + fr * tg_yyy_xxyyyz_1[j] + fl1_fxn * tg_yyy_xyyyz_1[j];

                    tg_xyyy_xxyyzz_0[j] = pb_x * tg_yyy_xxyyzz_0[j] + fr * tg_yyy_xxyyzz_1[j] + fl1_fxn * tg_yyy_xyyzz_1[j];

                    tg_xyyy_xxyzzz_0[j] = pb_x * tg_yyy_xxyzzz_0[j] + fr * tg_yyy_xxyzzz_1[j] + fl1_fxn * tg_yyy_xyzzz_1[j];

                    tg_xyyy_xxzzzz_0[j] = pb_x * tg_yyy_xxzzzz_0[j] + fr * tg_yyy_xxzzzz_1[j] + fl1_fxn * tg_yyy_xzzzz_1[j];

                    tg_xyyy_xyyyyy_0[j] = pb_x * tg_yyy_xyyyyy_0[j] + fr * tg_yyy_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyy_1[j];

                    tg_xyyy_xyyyyz_0[j] = pb_x * tg_yyy_xyyyyz_0[j] + fr * tg_yyy_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyyz_1[j];

                    tg_xyyy_xyyyzz_0[j] = pb_x * tg_yyy_xyyyzz_0[j] + fr * tg_yyy_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyyzz_1[j];

                    tg_xyyy_xyyzzz_0[j] = pb_x * tg_yyy_xyyzzz_0[j] + fr * tg_yyy_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyzzz_1[j];

                    tg_xyyy_xyzzzz_0[j] = pb_x * tg_yyy_xyzzzz_0[j] + fr * tg_yyy_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yzzzz_1[j];

                    tg_xyyy_xzzzzz_0[j] = pb_x * tg_yyy_xzzzzz_0[j] + fr * tg_yyy_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zzzzz_1[j];

                    tg_xyyy_yyyyyy_0[j] = pb_x * tg_yyy_yyyyyy_0[j] + fr * tg_yyy_yyyyyy_1[j];

                    tg_xyyy_yyyyyz_0[j] = pb_x * tg_yyy_yyyyyz_0[j] + fr * tg_yyy_yyyyyz_1[j];

                    tg_xyyy_yyyyzz_0[j] = pb_x * tg_yyy_yyyyzz_0[j] + fr * tg_yyy_yyyyzz_1[j];

                    tg_xyyy_yyyzzz_0[j] = pb_x * tg_yyy_yyyzzz_0[j] + fr * tg_yyy_yyyzzz_1[j];

                    tg_xyyy_yyzzzz_0[j] = pb_x * tg_yyy_yyzzzz_0[j] + fr * tg_yyy_yyzzzz_1[j];

                    tg_xyyy_yzzzzz_0[j] = pb_x * tg_yyy_yzzzzz_0[j] + fr * tg_yyy_yzzzzz_1[j];

                    tg_xyyy_zzzzzz_0[j] = pb_x * tg_yyy_zzzzzz_0[j] + fr * tg_yyy_zzzzzz_1[j];

                    tg_xyyz_xxxxxx_0[j] = pb_x * tg_yyz_xxxxxx_0[j] + fr * tg_yyz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyz_xxxxx_1[j];

                    tg_xyyz_xxxxxy_0[j] = pb_x * tg_yyz_xxxxxy_0[j] + fr * tg_yyz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxy_1[j];

                    tg_xyyz_xxxxxz_0[j] = pb_x * tg_yyz_xxxxxz_0[j] + fr * tg_yyz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyz_xxxxz_1[j];

                    tg_xyyz_xxxxyy_0[j] = pb_x * tg_yyz_xxxxyy_0[j] + fr * tg_yyz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyy_1[j];

                    tg_xyyz_xxxxyz_0[j] = pb_x * tg_yyz_xxxxyz_0[j] + fr * tg_yyz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxyz_1[j];

                    tg_xyyz_xxxxzz_0[j] = pb_x * tg_yyz_xxxxzz_0[j] + fr * tg_yyz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyz_xxxzz_1[j];

                    tg_xyyz_xxxyyy_0[j] = pb_x * tg_yyz_xxxyyy_0[j] + fr * tg_yyz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyy_1[j];

                    tg_xyyz_xxxyyz_0[j] = pb_x * tg_yyz_xxxyyz_0[j] + fr * tg_yyz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyyz_1[j];

                    tg_xyyz_xxxyzz_0[j] = pb_x * tg_yyz_xxxyzz_0[j] + fr * tg_yyz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxyzz_1[j];

                    tg_xyyz_xxxzzz_0[j] = pb_x * tg_yyz_xxxzzz_0[j] + fr * tg_yyz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxzzz_1[j];

                    tg_xyyz_xxyyyy_0[j] = pb_x * tg_yyz_xxyyyy_0[j] + fr * tg_yyz_xxyyyy_1[j] + fl1_fxn * tg_yyz_xyyyy_1[j];

                    tg_xyyz_xxyyyz_0[j] = pb_x * tg_yyz_xxyyyz_0[j] + fr * tg_yyz_xxyyyz_1[j] + fl1_fxn * tg_yyz_xyyyz_1[j];

                    tg_xyyz_xxyyzz_0[j] = pb_x * tg_yyz_xxyyzz_0[j] + fr * tg_yyz_xxyyzz_1[j] + fl1_fxn * tg_yyz_xyyzz_1[j];

                    tg_xyyz_xxyzzz_0[j] = pb_x * tg_yyz_xxyzzz_0[j] + fr * tg_yyz_xxyzzz_1[j] + fl1_fxn * tg_yyz_xyzzz_1[j];

                    tg_xyyz_xxzzzz_0[j] = pb_x * tg_yyz_xxzzzz_0[j] + fr * tg_yyz_xxzzzz_1[j] + fl1_fxn * tg_yyz_xzzzz_1[j];

                    tg_xyyz_xyyyyy_0[j] = pb_x * tg_yyz_xyyyyy_0[j] + fr * tg_yyz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyy_1[j];

                    tg_xyyz_xyyyyz_0[j] = pb_x * tg_yyz_xyyyyz_0[j] + fr * tg_yyz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyyz_1[j];

                    tg_xyyz_xyyyzz_0[j] = pb_x * tg_yyz_xyyyzz_0[j] + fr * tg_yyz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyyzz_1[j];

                    tg_xyyz_xyyzzz_0[j] = pb_x * tg_yyz_xyyzzz_0[j] + fr * tg_yyz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyzzz_1[j];

                    tg_xyyz_xyzzzz_0[j] = pb_x * tg_yyz_xyzzzz_0[j] + fr * tg_yyz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yzzzz_1[j];

                    tg_xyyz_xzzzzz_0[j] = pb_x * tg_yyz_xzzzzz_0[j] + fr * tg_yyz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zzzzz_1[j];

                    tg_xyyz_yyyyyy_0[j] = pb_x * tg_yyz_yyyyyy_0[j] + fr * tg_yyz_yyyyyy_1[j];

                    tg_xyyz_yyyyyz_0[j] = pb_x * tg_yyz_yyyyyz_0[j] + fr * tg_yyz_yyyyyz_1[j];

                    tg_xyyz_yyyyzz_0[j] = pb_x * tg_yyz_yyyyzz_0[j] + fr * tg_yyz_yyyyzz_1[j];

                    tg_xyyz_yyyzzz_0[j] = pb_x * tg_yyz_yyyzzz_0[j] + fr * tg_yyz_yyyzzz_1[j];

                    tg_xyyz_yyzzzz_0[j] = pb_x * tg_yyz_yyzzzz_0[j] + fr * tg_yyz_yyzzzz_1[j];

                    tg_xyyz_yzzzzz_0[j] = pb_x * tg_yyz_yzzzzz_0[j] + fr * tg_yyz_yzzzzz_1[j];

                    tg_xyyz_zzzzzz_0[j] = pb_x * tg_yyz_zzzzzz_0[j] + fr * tg_yyz_zzzzzz_1[j];

                    tg_xyzz_xxxxxx_0[j] = pb_x * tg_yzz_xxxxxx_0[j] + fr * tg_yzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yzz_xxxxx_1[j];

                    tg_xyzz_xxxxxy_0[j] = pb_x * tg_yzz_xxxxxy_0[j] + fr * tg_yzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxy_1[j];

                    tg_xyzz_xxxxxz_0[j] = pb_x * tg_yzz_xxxxxz_0[j] + fr * tg_yzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yzz_xxxxz_1[j];

                    tg_xyzz_xxxxyy_0[j] = pb_x * tg_yzz_xxxxyy_0[j] + fr * tg_yzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyy_1[j];

                    tg_xyzz_xxxxyz_0[j] = pb_x * tg_yzz_xxxxyz_0[j] + fr * tg_yzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxyz_1[j];

                    tg_xyzz_xxxxzz_0[j] = pb_x * tg_yzz_xxxxzz_0[j] + fr * tg_yzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yzz_xxxzz_1[j];

                    tg_xyzz_xxxyyy_0[j] = pb_x * tg_yzz_xxxyyy_0[j] + fr * tg_yzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyy_1[j];

                    tg_xyzz_xxxyyz_0[j] = pb_x * tg_yzz_xxxyyz_0[j] + fr * tg_yzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyyz_1[j];

                    tg_xyzz_xxxyzz_0[j] = pb_x * tg_yzz_xxxyzz_0[j] + fr * tg_yzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxyzz_1[j];

                    tg_xyzz_xxxzzz_0[j] = pb_x * tg_yzz_xxxzzz_0[j] + fr * tg_yzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxzzz_1[j];

                    tg_xyzz_xxyyyy_0[j] = pb_x * tg_yzz_xxyyyy_0[j] + fr * tg_yzz_xxyyyy_1[j] + fl1_fxn * tg_yzz_xyyyy_1[j];

                    tg_xyzz_xxyyyz_0[j] = pb_x * tg_yzz_xxyyyz_0[j] + fr * tg_yzz_xxyyyz_1[j] + fl1_fxn * tg_yzz_xyyyz_1[j];

                    tg_xyzz_xxyyzz_0[j] = pb_x * tg_yzz_xxyyzz_0[j] + fr * tg_yzz_xxyyzz_1[j] + fl1_fxn * tg_yzz_xyyzz_1[j];

                    tg_xyzz_xxyzzz_0[j] = pb_x * tg_yzz_xxyzzz_0[j] + fr * tg_yzz_xxyzzz_1[j] + fl1_fxn * tg_yzz_xyzzz_1[j];

                    tg_xyzz_xxzzzz_0[j] = pb_x * tg_yzz_xxzzzz_0[j] + fr * tg_yzz_xxzzzz_1[j] + fl1_fxn * tg_yzz_xzzzz_1[j];

                    tg_xyzz_xyyyyy_0[j] = pb_x * tg_yzz_xyyyyy_0[j] + fr * tg_yzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyy_1[j];

                    tg_xyzz_xyyyyz_0[j] = pb_x * tg_yzz_xyyyyz_0[j] + fr * tg_yzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyyz_1[j];

                    tg_xyzz_xyyyzz_0[j] = pb_x * tg_yzz_xyyyzz_0[j] + fr * tg_yzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyyzz_1[j];

                    tg_xyzz_xyyzzz_0[j] = pb_x * tg_yzz_xyyzzz_0[j] + fr * tg_yzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyzzz_1[j];

                    tg_xyzz_xyzzzz_0[j] = pb_x * tg_yzz_xyzzzz_0[j] + fr * tg_yzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yzzzz_1[j];

                    tg_xyzz_xzzzzz_0[j] = pb_x * tg_yzz_xzzzzz_0[j] + fr * tg_yzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zzzzz_1[j];

                    tg_xyzz_yyyyyy_0[j] = pb_x * tg_yzz_yyyyyy_0[j] + fr * tg_yzz_yyyyyy_1[j];

                    tg_xyzz_yyyyyz_0[j] = pb_x * tg_yzz_yyyyyz_0[j] + fr * tg_yzz_yyyyyz_1[j];

                    tg_xyzz_yyyyzz_0[j] = pb_x * tg_yzz_yyyyzz_0[j] + fr * tg_yzz_yyyyzz_1[j];

                    tg_xyzz_yyyzzz_0[j] = pb_x * tg_yzz_yyyzzz_0[j] + fr * tg_yzz_yyyzzz_1[j];

                    tg_xyzz_yyzzzz_0[j] = pb_x * tg_yzz_yyzzzz_0[j] + fr * tg_yzz_yyzzzz_1[j];

                    tg_xyzz_yzzzzz_0[j] = pb_x * tg_yzz_yzzzzz_0[j] + fr * tg_yzz_yzzzzz_1[j];

                    tg_xyzz_zzzzzz_0[j] = pb_x * tg_yzz_zzzzzz_0[j] + fr * tg_yzz_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSI_252_336(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (252,336)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
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

                auto tg_yyy_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 168); 

                auto tg_yyy_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 169); 

                auto tg_yyy_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 170); 

                auto tg_yyy_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 171); 

                auto tg_yyy_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 172); 

                auto tg_yyy_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 173); 

                auto tg_yyy_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 174); 

                auto tg_yyy_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 175); 

                auto tg_yyy_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 176); 

                auto tg_yyy_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 177); 

                auto tg_yyy_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 178); 

                auto tg_yyy_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 179); 

                auto tg_yyy_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 180); 

                auto tg_yyy_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 181); 

                auto tg_yyy_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 182); 

                auto tg_yyy_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 183); 

                auto tg_yyy_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 184); 

                auto tg_yyy_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 185); 

                auto tg_yyy_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 186); 

                auto tg_yyy_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 187); 

                auto tg_yyy_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 188); 

                auto tg_yyy_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 189); 

                auto tg_yyy_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 190); 

                auto tg_yyy_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 191); 

                auto tg_yyy_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 192); 

                auto tg_yyy_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 193); 

                auto tg_yyy_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 194); 

                auto tg_yyy_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 195); 

                auto tg_yyz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 196); 

                auto tg_yyz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 197); 

                auto tg_yyz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 198); 

                auto tg_yyz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 199); 

                auto tg_yyz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 200); 

                auto tg_yyz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 201); 

                auto tg_yyz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 202); 

                auto tg_yyz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 203); 

                auto tg_yyz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 204); 

                auto tg_yyz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 205); 

                auto tg_yyz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 206); 

                auto tg_yyz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 207); 

                auto tg_yyz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 208); 

                auto tg_yyz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 209); 

                auto tg_yyz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 210); 

                auto tg_yyz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 211); 

                auto tg_yyz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 212); 

                auto tg_yyz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 213); 

                auto tg_yyz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 214); 

                auto tg_yyz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 215); 

                auto tg_yyz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 216); 

                auto tg_yyz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 217); 

                auto tg_yyz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 218); 

                auto tg_yyz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 219); 

                auto tg_yyz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 220); 

                auto tg_yyz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 221); 

                auto tg_yyz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 222); 

                auto tg_yyz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 223); 

                auto tg_zzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 252); 

                auto tg_zzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 253); 

                auto tg_zzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 254); 

                auto tg_zzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 255); 

                auto tg_zzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 256); 

                auto tg_zzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 257); 

                auto tg_zzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 258); 

                auto tg_zzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 259); 

                auto tg_zzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 260); 

                auto tg_zzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 261); 

                auto tg_zzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 262); 

                auto tg_zzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 263); 

                auto tg_zzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 264); 

                auto tg_zzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 265); 

                auto tg_zzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 266); 

                auto tg_zzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 267); 

                auto tg_zzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 268); 

                auto tg_zzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 269); 

                auto tg_zzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 270); 

                auto tg_zzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 271); 

                auto tg_zzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 272); 

                auto tg_zzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 273); 

                auto tg_zzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 274); 

                auto tg_zzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 275); 

                auto tg_zzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 276); 

                auto tg_zzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 277); 

                auto tg_zzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 278); 

                auto tg_zzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 279); 

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

                // set up pointers to integrals

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

                // Batch of Integrals (252,336)

                #pragma omp simd aligned(fxn, fza, tg_xzzz_xxxxxx_0, tg_xzzz_xxxxxy_0, tg_xzzz_xxxxxz_0, \
                                         tg_xzzz_xxxxyy_0, tg_xzzz_xxxxyz_0, tg_xzzz_xxxxzz_0, tg_xzzz_xxxyyy_0, \
                                         tg_xzzz_xxxyyz_0, tg_xzzz_xxxyzz_0, tg_xzzz_xxxzzz_0, tg_xzzz_xxyyyy_0, \
                                         tg_xzzz_xxyyyz_0, tg_xzzz_xxyyzz_0, tg_xzzz_xxyzzz_0, tg_xzzz_xxzzzz_0, \
                                         tg_xzzz_xyyyyy_0, tg_xzzz_xyyyyz_0, tg_xzzz_xyyyzz_0, tg_xzzz_xyyzzz_0, \
                                         tg_xzzz_xyzzzz_0, tg_xzzz_xzzzzz_0, tg_xzzz_yyyyyy_0, tg_xzzz_yyyyyz_0, \
                                         tg_xzzz_yyyyzz_0, tg_xzzz_yyyzzz_0, tg_xzzz_yyzzzz_0, tg_xzzz_yzzzzz_0, \
                                         tg_xzzz_zzzzzz_0, tg_yy_xxxxxx_0, tg_yy_xxxxxx_1, tg_yy_xxxxxy_0, tg_yy_xxxxxy_1, \
                                         tg_yy_xxxxxz_0, tg_yy_xxxxxz_1, tg_yy_xxxxyy_0, tg_yy_xxxxyy_1, tg_yy_xxxxyz_0, \
                                         tg_yy_xxxxyz_1, tg_yy_xxxxzz_0, tg_yy_xxxxzz_1, tg_yy_xxxyyy_0, tg_yy_xxxyyy_1, \
                                         tg_yy_xxxyyz_0, tg_yy_xxxyyz_1, tg_yy_xxxyzz_0, tg_yy_xxxyzz_1, tg_yy_xxxzzz_0, \
                                         tg_yy_xxxzzz_1, tg_yy_xxyyyy_0, tg_yy_xxyyyy_1, tg_yy_xxyyyz_0, tg_yy_xxyyyz_1, \
                                         tg_yy_xxyyzz_0, tg_yy_xxyyzz_1, tg_yy_xxyzzz_0, tg_yy_xxyzzz_1, tg_yy_xxzzzz_0, \
                                         tg_yy_xxzzzz_1, tg_yy_xyyyyy_0, tg_yy_xyyyyy_1, tg_yy_xyyyyz_0, tg_yy_xyyyyz_1, \
                                         tg_yy_xyyyzz_0, tg_yy_xyyyzz_1, tg_yy_xyyzzz_0, tg_yy_xyyzzz_1, tg_yy_xyzzzz_0, \
                                         tg_yy_xyzzzz_1, tg_yy_xzzzzz_0, tg_yy_xzzzzz_1, tg_yy_yyyyyy_0, tg_yy_yyyyyy_1, \
                                         tg_yy_yyyyyz_0, tg_yy_yyyyyz_1, tg_yy_yyyyzz_0, tg_yy_yyyyzz_1, tg_yy_yyyzzz_0, \
                                         tg_yy_yyyzzz_1, tg_yy_yyzzzz_0, tg_yy_yyzzzz_1, tg_yy_yzzzzz_0, tg_yy_yzzzzz_1, \
                                         tg_yy_zzzzzz_0, tg_yy_zzzzzz_1, tg_yyy_xxxxx_1, tg_yyy_xxxxxx_0, tg_yyy_xxxxxx_1, \
                                         tg_yyy_xxxxxy_0, tg_yyy_xxxxxy_1, tg_yyy_xxxxxz_0, tg_yyy_xxxxxz_1, tg_yyy_xxxxy_1, \
                                         tg_yyy_xxxxyy_0, tg_yyy_xxxxyy_1, tg_yyy_xxxxyz_0, tg_yyy_xxxxyz_1, tg_yyy_xxxxz_1, \
                                         tg_yyy_xxxxzz_0, tg_yyy_xxxxzz_1, tg_yyy_xxxyy_1, tg_yyy_xxxyyy_0, tg_yyy_xxxyyy_1, \
                                         tg_yyy_xxxyyz_0, tg_yyy_xxxyyz_1, tg_yyy_xxxyz_1, tg_yyy_xxxyzz_0, tg_yyy_xxxyzz_1, \
                                         tg_yyy_xxxzz_1, tg_yyy_xxxzzz_0, tg_yyy_xxxzzz_1, tg_yyy_xxyyy_1, tg_yyy_xxyyyy_0, \
                                         tg_yyy_xxyyyy_1, tg_yyy_xxyyyz_0, tg_yyy_xxyyyz_1, tg_yyy_xxyyz_1, tg_yyy_xxyyzz_0, \
                                         tg_yyy_xxyyzz_1, tg_yyy_xxyzz_1, tg_yyy_xxyzzz_0, tg_yyy_xxyzzz_1, tg_yyy_xxzzz_1, \
                                         tg_yyy_xxzzzz_0, tg_yyy_xxzzzz_1, tg_yyy_xyyyy_1, tg_yyy_xyyyyy_0, tg_yyy_xyyyyy_1, \
                                         tg_yyy_xyyyyz_0, tg_yyy_xyyyyz_1, tg_yyy_xyyyz_1, tg_yyy_xyyyzz_0, tg_yyy_xyyyzz_1, \
                                         tg_yyy_xyyzz_1, tg_yyy_xyyzzz_0, tg_yyy_xyyzzz_1, tg_yyy_xyzzz_1, tg_yyy_xyzzzz_0, \
                                         tg_yyy_xyzzzz_1, tg_yyy_xzzzz_1, tg_yyy_xzzzzz_0, tg_yyy_xzzzzz_1, tg_yyy_yyyyy_1, \
                                         tg_yyy_yyyyyy_0, tg_yyy_yyyyyy_1, tg_yyy_yyyyyz_0, tg_yyy_yyyyyz_1, tg_yyy_yyyyz_1, \
                                         tg_yyy_yyyyzz_0, tg_yyy_yyyyzz_1, tg_yyy_yyyzz_1, tg_yyy_yyyzzz_0, tg_yyy_yyyzzz_1, \
                                         tg_yyy_yyzzz_1, tg_yyy_yyzzzz_0, tg_yyy_yyzzzz_1, tg_yyy_yzzzz_1, tg_yyy_yzzzzz_0, \
                                         tg_yyy_yzzzzz_1, tg_yyy_zzzzz_1, tg_yyy_zzzzzz_0, tg_yyy_zzzzzz_1, tg_yyyy_xxxxxx_0, \
                                         tg_yyyy_xxxxxy_0, tg_yyyy_xxxxxz_0, tg_yyyy_xxxxyy_0, tg_yyyy_xxxxyz_0, \
                                         tg_yyyy_xxxxzz_0, tg_yyyy_xxxyyy_0, tg_yyyy_xxxyyz_0, tg_yyyy_xxxyzz_0, \
                                         tg_yyyy_xxxzzz_0, tg_yyyy_xxyyyy_0, tg_yyyy_xxyyyz_0, tg_yyyy_xxyyzz_0, \
                                         tg_yyyy_xxyzzz_0, tg_yyyy_xxzzzz_0, tg_yyyy_xyyyyy_0, tg_yyyy_xyyyyz_0, \
                                         tg_yyyy_xyyyzz_0, tg_yyyy_xyyzzz_0, tg_yyyy_xyzzzz_0, tg_yyyy_xzzzzz_0, \
                                         tg_yyyy_yyyyyy_0, tg_yyyy_yyyyyz_0, tg_yyyy_yyyyzz_0, tg_yyyy_yyyzzz_0, \
                                         tg_yyyy_yyzzzz_0, tg_yyyy_yzzzzz_0, tg_yyyy_zzzzzz_0, tg_yyyz_xxxxxx_0, \
                                         tg_yyyz_xxxxxy_0, tg_yyyz_xxxxxz_0, tg_yyyz_xxxxyy_0, tg_yyyz_xxxxyz_0, \
                                         tg_yyyz_xxxxzz_0, tg_yyyz_xxxyyy_0, tg_yyyz_xxxyyz_0, tg_yyyz_xxxyzz_0, \
                                         tg_yyyz_xxxzzz_0, tg_yyyz_xxyyyy_0, tg_yyyz_xxyyyz_0, tg_yyyz_xxyyzz_0, \
                                         tg_yyyz_xxyzzz_0, tg_yyyz_xxzzzz_0, tg_yyyz_xyyyyy_0, tg_yyyz_xyyyyz_0, \
                                         tg_yyyz_xyyyzz_0, tg_yyyz_xyyzzz_0, tg_yyyz_xyzzzz_0, tg_yyyz_xzzzzz_0, \
                                         tg_yyyz_yyyyyy_0, tg_yyyz_yyyyyz_0, tg_yyyz_yyyyzz_0, tg_yyyz_yyyzzz_0, \
                                         tg_yyyz_yyzzzz_0, tg_yyyz_yzzzzz_0, tg_yyyz_zzzzzz_0, tg_yyz_xxxxx_1, \
                                         tg_yyz_xxxxxx_0, tg_yyz_xxxxxx_1, tg_yyz_xxxxxy_0, tg_yyz_xxxxxy_1, tg_yyz_xxxxxz_0, \
                                         tg_yyz_xxxxxz_1, tg_yyz_xxxxy_1, tg_yyz_xxxxyy_0, tg_yyz_xxxxyy_1, tg_yyz_xxxxyz_0, \
                                         tg_yyz_xxxxyz_1, tg_yyz_xxxxz_1, tg_yyz_xxxxzz_0, tg_yyz_xxxxzz_1, tg_yyz_xxxyy_1, \
                                         tg_yyz_xxxyyy_0, tg_yyz_xxxyyy_1, tg_yyz_xxxyyz_0, tg_yyz_xxxyyz_1, tg_yyz_xxxyz_1, \
                                         tg_yyz_xxxyzz_0, tg_yyz_xxxyzz_1, tg_yyz_xxxzz_1, tg_yyz_xxxzzz_0, tg_yyz_xxxzzz_1, \
                                         tg_yyz_xxyyy_1, tg_yyz_xxyyyy_0, tg_yyz_xxyyyy_1, tg_yyz_xxyyyz_0, tg_yyz_xxyyyz_1, \
                                         tg_yyz_xxyyz_1, tg_yyz_xxyyzz_0, tg_yyz_xxyyzz_1, tg_yyz_xxyzz_1, tg_yyz_xxyzzz_0, \
                                         tg_yyz_xxyzzz_1, tg_yyz_xxzzz_1, tg_yyz_xxzzzz_0, tg_yyz_xxzzzz_1, tg_yyz_xyyyy_1, \
                                         tg_yyz_xyyyyy_0, tg_yyz_xyyyyy_1, tg_yyz_xyyyyz_0, tg_yyz_xyyyyz_1, tg_yyz_xyyyz_1, \
                                         tg_yyz_xyyyzz_0, tg_yyz_xyyyzz_1, tg_yyz_xyyzz_1, tg_yyz_xyyzzz_0, tg_yyz_xyyzzz_1, \
                                         tg_yyz_xyzzz_1, tg_yyz_xyzzzz_0, tg_yyz_xyzzzz_1, tg_yyz_xzzzz_1, tg_yyz_xzzzzz_0, \
                                         tg_yyz_xzzzzz_1, tg_yyz_yyyyy_1, tg_yyz_yyyyyy_0, tg_yyz_yyyyyy_1, tg_yyz_yyyyyz_0, \
                                         tg_yyz_yyyyyz_1, tg_yyz_yyyyz_1, tg_yyz_yyyyzz_0, tg_yyz_yyyyzz_1, tg_yyz_yyyzz_1, \
                                         tg_yyz_yyyzzz_0, tg_yyz_yyyzzz_1, tg_yyz_yyzzz_1, tg_yyz_yyzzzz_0, tg_yyz_yyzzzz_1, \
                                         tg_yyz_yzzzz_1, tg_yyz_yzzzzz_0, tg_yyz_yzzzzz_1, tg_yyz_zzzzz_1, tg_yyz_zzzzzz_0, \
                                         tg_yyz_zzzzzz_1, tg_yz_xxxxxx_0, tg_yz_xxxxxx_1, tg_yz_xxxxxy_0, tg_yz_xxxxxy_1, \
                                         tg_yz_xxxxxz_0, tg_yz_xxxxxz_1, tg_yz_xxxxyy_0, tg_yz_xxxxyy_1, tg_yz_xxxxyz_0, \
                                         tg_yz_xxxxyz_1, tg_yz_xxxxzz_0, tg_yz_xxxxzz_1, tg_yz_xxxyyy_0, tg_yz_xxxyyy_1, \
                                         tg_yz_xxxyyz_0, tg_yz_xxxyyz_1, tg_yz_xxxyzz_0, tg_yz_xxxyzz_1, tg_yz_xxxzzz_0, \
                                         tg_yz_xxxzzz_1, tg_yz_xxyyyy_0, tg_yz_xxyyyy_1, tg_yz_xxyyyz_0, tg_yz_xxyyyz_1, \
                                         tg_yz_xxyyzz_0, tg_yz_xxyyzz_1, tg_yz_xxyzzz_0, tg_yz_xxyzzz_1, tg_yz_xxzzzz_0, \
                                         tg_yz_xxzzzz_1, tg_yz_xyyyyy_0, tg_yz_xyyyyy_1, tg_yz_xyyyyz_0, tg_yz_xyyyyz_1, \
                                         tg_yz_xyyyzz_0, tg_yz_xyyyzz_1, tg_yz_xyyzzz_0, tg_yz_xyyzzz_1, tg_yz_xyzzzz_0, \
                                         tg_yz_xyzzzz_1, tg_yz_xzzzzz_0, tg_yz_xzzzzz_1, tg_yz_yyyyyy_0, tg_yz_yyyyyy_1, \
                                         tg_yz_yyyyyz_0, tg_yz_yyyyyz_1, tg_yz_yyyyzz_0, tg_yz_yyyyzz_1, tg_yz_yyyzzz_0, \
                                         tg_yz_yyyzzz_1, tg_yz_yyzzzz_0, tg_yz_yyzzzz_1, tg_yz_yzzzzz_0, tg_yz_yzzzzz_1, \
                                         tg_yz_zzzzzz_0, tg_yz_zzzzzz_1, tg_zzz_xxxxx_1, tg_zzz_xxxxxx_0, tg_zzz_xxxxxx_1, \
                                         tg_zzz_xxxxxy_0, tg_zzz_xxxxxy_1, tg_zzz_xxxxxz_0, tg_zzz_xxxxxz_1, tg_zzz_xxxxy_1, \
                                         tg_zzz_xxxxyy_0, tg_zzz_xxxxyy_1, tg_zzz_xxxxyz_0, tg_zzz_xxxxyz_1, tg_zzz_xxxxz_1, \
                                         tg_zzz_xxxxzz_0, tg_zzz_xxxxzz_1, tg_zzz_xxxyy_1, tg_zzz_xxxyyy_0, tg_zzz_xxxyyy_1, \
                                         tg_zzz_xxxyyz_0, tg_zzz_xxxyyz_1, tg_zzz_xxxyz_1, tg_zzz_xxxyzz_0, tg_zzz_xxxyzz_1, \
                                         tg_zzz_xxxzz_1, tg_zzz_xxxzzz_0, tg_zzz_xxxzzz_1, tg_zzz_xxyyy_1, tg_zzz_xxyyyy_0, \
                                         tg_zzz_xxyyyy_1, tg_zzz_xxyyyz_0, tg_zzz_xxyyyz_1, tg_zzz_xxyyz_1, tg_zzz_xxyyzz_0, \
                                         tg_zzz_xxyyzz_1, tg_zzz_xxyzz_1, tg_zzz_xxyzzz_0, tg_zzz_xxyzzz_1, tg_zzz_xxzzz_1, \
                                         tg_zzz_xxzzzz_0, tg_zzz_xxzzzz_1, tg_zzz_xyyyy_1, tg_zzz_xyyyyy_0, tg_zzz_xyyyyy_1, \
                                         tg_zzz_xyyyyz_0, tg_zzz_xyyyyz_1, tg_zzz_xyyyz_1, tg_zzz_xyyyzz_0, tg_zzz_xyyyzz_1, \
                                         tg_zzz_xyyzz_1, tg_zzz_xyyzzz_0, tg_zzz_xyyzzz_1, tg_zzz_xyzzz_1, tg_zzz_xyzzzz_0, \
                                         tg_zzz_xyzzzz_1, tg_zzz_xzzzz_1, tg_zzz_xzzzzz_0, tg_zzz_xzzzzz_1, tg_zzz_yyyyy_1, \
                                         tg_zzz_yyyyyy_0, tg_zzz_yyyyyy_1, tg_zzz_yyyyyz_0, tg_zzz_yyyyyz_1, tg_zzz_yyyyz_1, \
                                         tg_zzz_yyyyzz_0, tg_zzz_yyyyzz_1, tg_zzz_yyyzz_1, tg_zzz_yyyzzz_0, tg_zzz_yyyzzz_1, \
                                         tg_zzz_yyzzz_1, tg_zzz_yyzzzz_0, tg_zzz_yyzzzz_1, tg_zzz_yzzzz_1, tg_zzz_yzzzzz_0, \
                                         tg_zzz_yzzzzz_1, tg_zzz_zzzzz_1, tg_zzz_zzzzzz_0, tg_zzz_zzzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xzzz_xxxxxx_0[j] = pb_x * tg_zzz_xxxxxx_0[j] + fr * tg_zzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_zzz_xxxxx_1[j];

                    tg_xzzz_xxxxxy_0[j] = pb_x * tg_zzz_xxxxxy_0[j] + fr * tg_zzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxy_1[j];

                    tg_xzzz_xxxxxz_0[j] = pb_x * tg_zzz_xxxxxz_0[j] + fr * tg_zzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_zzz_xxxxz_1[j];

                    tg_xzzz_xxxxyy_0[j] = pb_x * tg_zzz_xxxxyy_0[j] + fr * tg_zzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyy_1[j];

                    tg_xzzz_xxxxyz_0[j] = pb_x * tg_zzz_xxxxyz_0[j] + fr * tg_zzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxyz_1[j];

                    tg_xzzz_xxxxzz_0[j] = pb_x * tg_zzz_xxxxzz_0[j] + fr * tg_zzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_zzz_xxxzz_1[j];

                    tg_xzzz_xxxyyy_0[j] = pb_x * tg_zzz_xxxyyy_0[j] + fr * tg_zzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyy_1[j];

                    tg_xzzz_xxxyyz_0[j] = pb_x * tg_zzz_xxxyyz_0[j] + fr * tg_zzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyz_1[j];

                    tg_xzzz_xxxyzz_0[j] = pb_x * tg_zzz_xxxyzz_0[j] + fr * tg_zzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyzz_1[j];

                    tg_xzzz_xxxzzz_0[j] = pb_x * tg_zzz_xxxzzz_0[j] + fr * tg_zzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxzzz_1[j];

                    tg_xzzz_xxyyyy_0[j] = pb_x * tg_zzz_xxyyyy_0[j] + fr * tg_zzz_xxyyyy_1[j] + fl1_fxn * tg_zzz_xyyyy_1[j];

                    tg_xzzz_xxyyyz_0[j] = pb_x * tg_zzz_xxyyyz_0[j] + fr * tg_zzz_xxyyyz_1[j] + fl1_fxn * tg_zzz_xyyyz_1[j];

                    tg_xzzz_xxyyzz_0[j] = pb_x * tg_zzz_xxyyzz_0[j] + fr * tg_zzz_xxyyzz_1[j] + fl1_fxn * tg_zzz_xyyzz_1[j];

                    tg_xzzz_xxyzzz_0[j] = pb_x * tg_zzz_xxyzzz_0[j] + fr * tg_zzz_xxyzzz_1[j] + fl1_fxn * tg_zzz_xyzzz_1[j];

                    tg_xzzz_xxzzzz_0[j] = pb_x * tg_zzz_xxzzzz_0[j] + fr * tg_zzz_xxzzzz_1[j] + fl1_fxn * tg_zzz_xzzzz_1[j];

                    tg_xzzz_xyyyyy_0[j] = pb_x * tg_zzz_xyyyyy_0[j] + fr * tg_zzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyy_1[j];

                    tg_xzzz_xyyyyz_0[j] = pb_x * tg_zzz_xyyyyz_0[j] + fr * tg_zzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyyz_1[j];

                    tg_xzzz_xyyyzz_0[j] = pb_x * tg_zzz_xyyyzz_0[j] + fr * tg_zzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyyzz_1[j];

                    tg_xzzz_xyyzzz_0[j] = pb_x * tg_zzz_xyyzzz_0[j] + fr * tg_zzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyzzz_1[j];

                    tg_xzzz_xyzzzz_0[j] = pb_x * tg_zzz_xyzzzz_0[j] + fr * tg_zzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yzzzz_1[j];

                    tg_xzzz_xzzzzz_0[j] = pb_x * tg_zzz_xzzzzz_0[j] + fr * tg_zzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzzz_1[j];

                    tg_xzzz_yyyyyy_0[j] = pb_x * tg_zzz_yyyyyy_0[j] + fr * tg_zzz_yyyyyy_1[j];

                    tg_xzzz_yyyyyz_0[j] = pb_x * tg_zzz_yyyyyz_0[j] + fr * tg_zzz_yyyyyz_1[j];

                    tg_xzzz_yyyyzz_0[j] = pb_x * tg_zzz_yyyyzz_0[j] + fr * tg_zzz_yyyyzz_1[j];

                    tg_xzzz_yyyzzz_0[j] = pb_x * tg_zzz_yyyzzz_0[j] + fr * tg_zzz_yyyzzz_1[j];

                    tg_xzzz_yyzzzz_0[j] = pb_x * tg_zzz_yyzzzz_0[j] + fr * tg_zzz_yyzzzz_1[j];

                    tg_xzzz_yzzzzz_0[j] = pb_x * tg_zzz_yzzzzz_0[j] + fr * tg_zzz_yzzzzz_1[j];

                    tg_xzzz_zzzzzz_0[j] = pb_x * tg_zzz_zzzzzz_0[j] + fr * tg_zzz_zzzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyy_xxxxxx_0[j] = pb_y * tg_yyy_xxxxxx_0[j] + fr * tg_yyy_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxx_0[j] - tg_yy_xxxxxx_1[j] * fl1_fza);

                    tg_yyyy_xxxxxy_0[j] = pb_y * tg_yyy_xxxxxy_0[j] + fr * tg_yyy_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxy_0[j] - tg_yy_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxxx_1[j];

                    tg_yyyy_xxxxxz_0[j] = pb_y * tg_yyy_xxxxxz_0[j] + fr * tg_yyy_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxxz_0[j] - tg_yy_xxxxxz_1[j] * fl1_fza);

                    tg_yyyy_xxxxyy_0[j] = pb_y * tg_yyy_xxxxyy_0[j] + fr * tg_yyy_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxyy_0[j] - tg_yy_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxxxy_1[j];

                    tg_yyyy_xxxxyz_0[j] = pb_y * tg_yyy_xxxxyz_0[j] + fr * tg_yyy_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxyz_0[j] - tg_yy_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxxz_1[j];

                    tg_yyyy_xxxxzz_0[j] = pb_y * tg_yyy_xxxxzz_0[j] + fr * tg_yyy_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxxzz_0[j] - tg_yy_xxxxzz_1[j] * fl1_fza);

                    tg_yyyy_xxxyyy_0[j] = pb_y * tg_yyy_xxxyyy_0[j] + fr * tg_yyy_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyyy_0[j] - tg_yy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xxxyy_1[j];

                    tg_yyyy_xxxyyz_0[j] = pb_y * tg_yyy_xxxyyz_0[j] + fr * tg_yyy_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyyz_0[j] - tg_yy_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxxyz_1[j];

                    tg_yyyy_xxxyzz_0[j] = pb_y * tg_yyy_xxxyzz_0[j] + fr * tg_yyy_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxyzz_0[j] - tg_yy_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxxzz_1[j];

                    tg_yyyy_xxxzzz_0[j] = pb_y * tg_yyy_xxxzzz_0[j] + fr * tg_yyy_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxxzzz_0[j] - tg_yy_xxxzzz_1[j] * fl1_fza);

                    tg_yyyy_xxyyyy_0[j] = pb_y * tg_yyy_xxyyyy_0[j] + fr * tg_yyy_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyyy_0[j] - tg_yy_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_xxyyy_1[j];

                    tg_yyyy_xxyyyz_0[j] = pb_y * tg_yyy_xxyyyz_0[j] + fr * tg_yyy_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyyz_0[j] - tg_yy_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xxyyz_1[j];

                    tg_yyyy_xxyyzz_0[j] = pb_y * tg_yyy_xxyyzz_0[j] + fr * tg_yyy_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyyzz_0[j] - tg_yy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xxyzz_1[j];

                    tg_yyyy_xxyzzz_0[j] = pb_y * tg_yyy_xxyzzz_0[j] + fr * tg_yyy_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxyzzz_0[j] - tg_yy_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xxzzz_1[j];

                    tg_yyyy_xxzzzz_0[j] = pb_y * tg_yyy_xxzzzz_0[j] + fr * tg_yyy_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xxzzzz_0[j] - tg_yy_xxzzzz_1[j] * fl1_fza);

                    tg_yyyy_xyyyyy_0[j] = pb_y * tg_yyy_xyyyyy_0[j] + fr * tg_yyy_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyyy_0[j] - tg_yy_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyy_xyyyy_1[j];

                    tg_yyyy_xyyyyz_0[j] = pb_y * tg_yyy_xyyyyz_0[j] + fr * tg_yyy_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyyz_0[j] - tg_yy_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_xyyyz_1[j];

                    tg_yyyy_xyyyzz_0[j] = pb_y * tg_yyy_xyyyzz_0[j] + fr * tg_yyy_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyyzz_0[j] - tg_yy_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_xyyzz_1[j];

                    tg_yyyy_xyyzzz_0[j] = pb_y * tg_yyy_xyyzzz_0[j] + fr * tg_yyy_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyyzzz_0[j] - tg_yy_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xyzzz_1[j];

                    tg_yyyy_xyzzzz_0[j] = pb_y * tg_yyy_xyzzzz_0[j] + fr * tg_yyy_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xyzzzz_0[j] - tg_yy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xzzzz_1[j];

                    tg_yyyy_xzzzzz_0[j] = pb_y * tg_yyy_xzzzzz_0[j] + fr * tg_yyy_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_xzzzzz_0[j] - tg_yy_xzzzzz_1[j] * fl1_fza);

                    tg_yyyy_yyyyyy_0[j] = pb_y * tg_yyy_yyyyyy_0[j] + fr * tg_yyy_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyyy_0[j] - tg_yy_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyy_yyyyy_1[j];

                    tg_yyyy_yyyyyz_0[j] = pb_y * tg_yyy_yyyyyz_0[j] + fr * tg_yyy_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyyz_0[j] - tg_yy_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyy_yyyyz_1[j];

                    tg_yyyy_yyyyzz_0[j] = pb_y * tg_yyy_yyyyzz_0[j] + fr * tg_yyy_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyyzz_0[j] - tg_yy_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyy_yyyzz_1[j];

                    tg_yyyy_yyyzzz_0[j] = pb_y * tg_yyy_yyyzzz_0[j] + fr * tg_yyy_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyyzzz_0[j] - tg_yy_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_yyzzz_1[j];

                    tg_yyyy_yyzzzz_0[j] = pb_y * tg_yyy_yyzzzz_0[j] + fr * tg_yyy_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yyzzzz_0[j] - tg_yy_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_yzzzz_1[j];

                    tg_yyyy_yzzzzz_0[j] = pb_y * tg_yyy_yzzzzz_0[j] + fr * tg_yyy_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_yzzzzz_0[j] - tg_yy_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_zzzzz_1[j];

                    tg_yyyy_zzzzzz_0[j] = pb_y * tg_yyy_zzzzzz_0[j] + fr * tg_yyy_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_yy_zzzzzz_0[j] - tg_yy_zzzzzz_1[j] * fl1_fza);

                    tg_yyyz_xxxxxx_0[j] = pb_y * tg_yyz_xxxxxx_0[j] + fr * tg_yyz_xxxxxx_1[j] + fl1_fx * (tg_yz_xxxxxx_0[j] - tg_yz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyz_xxxxxy_0[j] = pb_y * tg_yyz_xxxxxy_0[j] + fr * tg_yyz_xxxxxy_1[j] + fl1_fx * (tg_yz_xxxxxy_0[j] - tg_yz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxxx_1[j];

                    tg_yyyz_xxxxxz_0[j] = pb_y * tg_yyz_xxxxxz_0[j] + fr * tg_yyz_xxxxxz_1[j] + fl1_fx * (tg_yz_xxxxxz_0[j] - tg_yz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyz_xxxxyy_0[j] = pb_y * tg_yyz_xxxxyy_0[j] + fr * tg_yyz_xxxxyy_1[j] + fl1_fx * (tg_yz_xxxxyy_0[j] - tg_yz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxxxy_1[j];

                    tg_yyyz_xxxxyz_0[j] = pb_y * tg_yyz_xxxxyz_0[j] + fr * tg_yyz_xxxxyz_1[j] + fl1_fx * (tg_yz_xxxxyz_0[j] - tg_yz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxxz_1[j];

                    tg_yyyz_xxxxzz_0[j] = pb_y * tg_yyz_xxxxzz_0[j] + fr * tg_yyz_xxxxzz_1[j] + fl1_fx * (tg_yz_xxxxzz_0[j] - tg_yz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyz_xxxyyy_0[j] = pb_y * tg_yyz_xxxyyy_0[j] + fr * tg_yyz_xxxyyy_1[j] + fl1_fx * (tg_yz_xxxyyy_0[j] - tg_yz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xxxyy_1[j];

                    tg_yyyz_xxxyyz_0[j] = pb_y * tg_yyz_xxxyyz_0[j] + fr * tg_yyz_xxxyyz_1[j] + fl1_fx * (tg_yz_xxxyyz_0[j] - tg_yz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxxyz_1[j];

                    tg_yyyz_xxxyzz_0[j] = pb_y * tg_yyz_xxxyzz_0[j] + fr * tg_yyz_xxxyzz_1[j] + fl1_fx * (tg_yz_xxxyzz_0[j] - tg_yz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxxzz_1[j];

                    tg_yyyz_xxxzzz_0[j] = pb_y * tg_yyz_xxxzzz_0[j] + fr * tg_yyz_xxxzzz_1[j] + fl1_fx * (tg_yz_xxxzzz_0[j] - tg_yz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyz_xxyyyy_0[j] = pb_y * tg_yyz_xxyyyy_0[j] + fr * tg_yyz_xxyyyy_1[j] + fl1_fx * (tg_yz_xxyyyy_0[j] - tg_yz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_xxyyy_1[j];

                    tg_yyyz_xxyyyz_0[j] = pb_y * tg_yyz_xxyyyz_0[j] + fr * tg_yyz_xxyyyz_1[j] + fl1_fx * (tg_yz_xxyyyz_0[j] - tg_yz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xxyyz_1[j];

                    tg_yyyz_xxyyzz_0[j] = pb_y * tg_yyz_xxyyzz_0[j] + fr * tg_yyz_xxyyzz_1[j] + fl1_fx * (tg_yz_xxyyzz_0[j] - tg_yz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xxyzz_1[j];

                    tg_yyyz_xxyzzz_0[j] = pb_y * tg_yyz_xxyzzz_0[j] + fr * tg_yyz_xxyzzz_1[j] + fl1_fx * (tg_yz_xxyzzz_0[j] - tg_yz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xxzzz_1[j];

                    tg_yyyz_xxzzzz_0[j] = pb_y * tg_yyz_xxzzzz_0[j] + fr * tg_yyz_xxzzzz_1[j] + fl1_fx * (tg_yz_xxzzzz_0[j] - tg_yz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyz_xyyyyy_0[j] = pb_y * tg_yyz_xyyyyy_0[j] + fr * tg_yyz_xyyyyy_1[j] + fl1_fx * (tg_yz_xyyyyy_0[j] - tg_yz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyz_xyyyy_1[j];

                    tg_yyyz_xyyyyz_0[j] = pb_y * tg_yyz_xyyyyz_0[j] + fr * tg_yyz_xyyyyz_1[j] + fl1_fx * (tg_yz_xyyyyz_0[j] - tg_yz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_xyyyz_1[j];

                    tg_yyyz_xyyyzz_0[j] = pb_y * tg_yyz_xyyyzz_0[j] + fr * tg_yyz_xyyyzz_1[j] + fl1_fx * (tg_yz_xyyyzz_0[j] - tg_yz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_xyyzz_1[j];

                    tg_yyyz_xyyzzz_0[j] = pb_y * tg_yyz_xyyzzz_0[j] + fr * tg_yyz_xyyzzz_1[j] + fl1_fx * (tg_yz_xyyzzz_0[j] - tg_yz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xyzzz_1[j];

                    tg_yyyz_xyzzzz_0[j] = pb_y * tg_yyz_xyzzzz_0[j] + fr * tg_yyz_xyzzzz_1[j] + fl1_fx * (tg_yz_xyzzzz_0[j] - tg_yz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xzzzz_1[j];

                    tg_yyyz_xzzzzz_0[j] = pb_y * tg_yyz_xzzzzz_0[j] + fr * tg_yyz_xzzzzz_1[j] + fl1_fx * (tg_yz_xzzzzz_0[j] - tg_yz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyz_yyyyyy_0[j] = pb_y * tg_yyz_yyyyyy_0[j] + fr * tg_yyz_yyyyyy_1[j] + fl1_fx * (tg_yz_yyyyyy_0[j] - tg_yz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyz_yyyyy_1[j];

                    tg_yyyz_yyyyyz_0[j] = pb_y * tg_yyz_yyyyyz_0[j] + fr * tg_yyz_yyyyyz_1[j] + fl1_fx * (tg_yz_yyyyyz_0[j] - tg_yz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyz_yyyyz_1[j];

                    tg_yyyz_yyyyzz_0[j] = pb_y * tg_yyz_yyyyzz_0[j] + fr * tg_yyz_yyyyzz_1[j] + fl1_fx * (tg_yz_yyyyzz_0[j] - tg_yz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyz_yyyzz_1[j];

                    tg_yyyz_yyyzzz_0[j] = pb_y * tg_yyz_yyyzzz_0[j] + fr * tg_yyz_yyyzzz_1[j] + fl1_fx * (tg_yz_yyyzzz_0[j] - tg_yz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_yyzzz_1[j];

                    tg_yyyz_yyzzzz_0[j] = pb_y * tg_yyz_yyzzzz_0[j] + fr * tg_yyz_yyzzzz_1[j] + fl1_fx * (tg_yz_yyzzzz_0[j] - tg_yz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_yzzzz_1[j];

                    tg_yyyz_yzzzzz_0[j] = pb_y * tg_yyz_yzzzzz_0[j] + fr * tg_yyz_yzzzzz_1[j] + fl1_fx * (tg_yz_yzzzzz_0[j] - tg_yz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_zzzzz_1[j];

                    tg_yyyz_zzzzzz_0[j] = pb_y * tg_yyz_zzzzzz_0[j] + fr * tg_yyz_zzzzzz_1[j] + fl1_fx * (tg_yz_zzzzzz_0[j] - tg_yz_zzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSI_336_420(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (336,420)

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
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

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

                auto tg_yzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 224); 

                auto tg_yzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 225); 

                auto tg_yzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 226); 

                auto tg_yzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 227); 

                auto tg_yzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 228); 

                auto tg_yzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 229); 

                auto tg_yzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 230); 

                auto tg_yzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 231); 

                auto tg_yzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 232); 

                auto tg_yzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 233); 

                auto tg_yzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 234); 

                auto tg_yzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 235); 

                auto tg_yzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 236); 

                auto tg_yzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 237); 

                auto tg_yzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 238); 

                auto tg_yzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 239); 

                auto tg_yzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 240); 

                auto tg_yzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 241); 

                auto tg_yzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 242); 

                auto tg_yzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 243); 

                auto tg_yzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 244); 

                auto tg_yzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 245); 

                auto tg_yzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 246); 

                auto tg_yzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 247); 

                auto tg_yzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 248); 

                auto tg_yzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 249); 

                auto tg_yzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 250); 

                auto tg_yzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 251); 

                auto tg_zzz_xxxxxx_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 252); 

                auto tg_zzz_xxxxxy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 253); 

                auto tg_zzz_xxxxxz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 254); 

                auto tg_zzz_xxxxyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 255); 

                auto tg_zzz_xxxxyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 256); 

                auto tg_zzz_xxxxzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 257); 

                auto tg_zzz_xxxyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 258); 

                auto tg_zzz_xxxyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 259); 

                auto tg_zzz_xxxyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 260); 

                auto tg_zzz_xxxzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 261); 

                auto tg_zzz_xxyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 262); 

                auto tg_zzz_xxyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 263); 

                auto tg_zzz_xxyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 264); 

                auto tg_zzz_xxyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 265); 

                auto tg_zzz_xxzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 266); 

                auto tg_zzz_xyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 267); 

                auto tg_zzz_xyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 268); 

                auto tg_zzz_xyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 269); 

                auto tg_zzz_xyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 270); 

                auto tg_zzz_xyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 271); 

                auto tg_zzz_xzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 272); 

                auto tg_zzz_yyyyyy_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 273); 

                auto tg_zzz_yyyyyz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 274); 

                auto tg_zzz_yyyyzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 275); 

                auto tg_zzz_yyyzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 276); 

                auto tg_zzz_yyzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 277); 

                auto tg_zzz_yzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 278); 

                auto tg_zzz_zzzzzz_1 = primBuffer[pidx_g_3_6_m1].data(280 * idx + 279); 

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

                // set up pointers to integrals

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

                // Batch of Integrals (336,420)

                #pragma omp simd aligned(fxn, fza, tg_yyzz_xxxxxx_0, tg_yyzz_xxxxxy_0, tg_yyzz_xxxxxz_0, \
                                         tg_yyzz_xxxxyy_0, tg_yyzz_xxxxyz_0, tg_yyzz_xxxxzz_0, tg_yyzz_xxxyyy_0, \
                                         tg_yyzz_xxxyyz_0, tg_yyzz_xxxyzz_0, tg_yyzz_xxxzzz_0, tg_yyzz_xxyyyy_0, \
                                         tg_yyzz_xxyyyz_0, tg_yyzz_xxyyzz_0, tg_yyzz_xxyzzz_0, tg_yyzz_xxzzzz_0, \
                                         tg_yyzz_xyyyyy_0, tg_yyzz_xyyyyz_0, tg_yyzz_xyyyzz_0, tg_yyzz_xyyzzz_0, \
                                         tg_yyzz_xyzzzz_0, tg_yyzz_xzzzzz_0, tg_yyzz_yyyyyy_0, tg_yyzz_yyyyyz_0, \
                                         tg_yyzz_yyyyzz_0, tg_yyzz_yyyzzz_0, tg_yyzz_yyzzzz_0, tg_yyzz_yzzzzz_0, \
                                         tg_yyzz_zzzzzz_0, tg_yzz_xxxxx_1, tg_yzz_xxxxxx_0, tg_yzz_xxxxxx_1, tg_yzz_xxxxxy_0, \
                                         tg_yzz_xxxxxy_1, tg_yzz_xxxxxz_0, tg_yzz_xxxxxz_1, tg_yzz_xxxxy_1, tg_yzz_xxxxyy_0, \
                                         tg_yzz_xxxxyy_1, tg_yzz_xxxxyz_0, tg_yzz_xxxxyz_1, tg_yzz_xxxxz_1, tg_yzz_xxxxzz_0, \
                                         tg_yzz_xxxxzz_1, tg_yzz_xxxyy_1, tg_yzz_xxxyyy_0, tg_yzz_xxxyyy_1, tg_yzz_xxxyyz_0, \
                                         tg_yzz_xxxyyz_1, tg_yzz_xxxyz_1, tg_yzz_xxxyzz_0, tg_yzz_xxxyzz_1, tg_yzz_xxxzz_1, \
                                         tg_yzz_xxxzzz_0, tg_yzz_xxxzzz_1, tg_yzz_xxyyy_1, tg_yzz_xxyyyy_0, tg_yzz_xxyyyy_1, \
                                         tg_yzz_xxyyyz_0, tg_yzz_xxyyyz_1, tg_yzz_xxyyz_1, tg_yzz_xxyyzz_0, tg_yzz_xxyyzz_1, \
                                         tg_yzz_xxyzz_1, tg_yzz_xxyzzz_0, tg_yzz_xxyzzz_1, tg_yzz_xxzzz_1, tg_yzz_xxzzzz_0, \
                                         tg_yzz_xxzzzz_1, tg_yzz_xyyyy_1, tg_yzz_xyyyyy_0, tg_yzz_xyyyyy_1, tg_yzz_xyyyyz_0, \
                                         tg_yzz_xyyyyz_1, tg_yzz_xyyyz_1, tg_yzz_xyyyzz_0, tg_yzz_xyyyzz_1, tg_yzz_xyyzz_1, \
                                         tg_yzz_xyyzzz_0, tg_yzz_xyyzzz_1, tg_yzz_xyzzz_1, tg_yzz_xyzzzz_0, tg_yzz_xyzzzz_1, \
                                         tg_yzz_xzzzz_1, tg_yzz_xzzzzz_0, tg_yzz_xzzzzz_1, tg_yzz_yyyyy_1, tg_yzz_yyyyyy_0, \
                                         tg_yzz_yyyyyy_1, tg_yzz_yyyyyz_0, tg_yzz_yyyyyz_1, tg_yzz_yyyyz_1, tg_yzz_yyyyzz_0, \
                                         tg_yzz_yyyyzz_1, tg_yzz_yyyzz_1, tg_yzz_yyyzzz_0, tg_yzz_yyyzzz_1, tg_yzz_yyzzz_1, \
                                         tg_yzz_yyzzzz_0, tg_yzz_yyzzzz_1, tg_yzz_yzzzz_1, tg_yzz_yzzzzz_0, tg_yzz_yzzzzz_1, \
                                         tg_yzz_zzzzz_1, tg_yzz_zzzzzz_0, tg_yzz_zzzzzz_1, tg_yzzz_xxxxxx_0, \
                                         tg_yzzz_xxxxxy_0, tg_yzzz_xxxxxz_0, tg_yzzz_xxxxyy_0, tg_yzzz_xxxxyz_0, \
                                         tg_yzzz_xxxxzz_0, tg_yzzz_xxxyyy_0, tg_yzzz_xxxyyz_0, tg_yzzz_xxxyzz_0, \
                                         tg_yzzz_xxxzzz_0, tg_yzzz_xxyyyy_0, tg_yzzz_xxyyyz_0, tg_yzzz_xxyyzz_0, \
                                         tg_yzzz_xxyzzz_0, tg_yzzz_xxzzzz_0, tg_yzzz_xyyyyy_0, tg_yzzz_xyyyyz_0, \
                                         tg_yzzz_xyyyzz_0, tg_yzzz_xyyzzz_0, tg_yzzz_xyzzzz_0, tg_yzzz_xzzzzz_0, \
                                         tg_yzzz_yyyyyy_0, tg_yzzz_yyyyyz_0, tg_yzzz_yyyyzz_0, tg_yzzz_yyyzzz_0, \
                                         tg_yzzz_yyzzzz_0, tg_yzzz_yzzzzz_0, tg_yzzz_zzzzzz_0, tg_zz_xxxxxx_0, tg_zz_xxxxxx_1, \
                                         tg_zz_xxxxxy_0, tg_zz_xxxxxy_1, tg_zz_xxxxxz_0, tg_zz_xxxxxz_1, tg_zz_xxxxyy_0, \
                                         tg_zz_xxxxyy_1, tg_zz_xxxxyz_0, tg_zz_xxxxyz_1, tg_zz_xxxxzz_0, tg_zz_xxxxzz_1, \
                                         tg_zz_xxxyyy_0, tg_zz_xxxyyy_1, tg_zz_xxxyyz_0, tg_zz_xxxyyz_1, tg_zz_xxxyzz_0, \
                                         tg_zz_xxxyzz_1, tg_zz_xxxzzz_0, tg_zz_xxxzzz_1, tg_zz_xxyyyy_0, tg_zz_xxyyyy_1, \
                                         tg_zz_xxyyyz_0, tg_zz_xxyyyz_1, tg_zz_xxyyzz_0, tg_zz_xxyyzz_1, tg_zz_xxyzzz_0, \
                                         tg_zz_xxyzzz_1, tg_zz_xxzzzz_0, tg_zz_xxzzzz_1, tg_zz_xyyyyy_0, tg_zz_xyyyyy_1, \
                                         tg_zz_xyyyyz_0, tg_zz_xyyyyz_1, tg_zz_xyyyzz_0, tg_zz_xyyyzz_1, tg_zz_xyyzzz_0, \
                                         tg_zz_xyyzzz_1, tg_zz_xyzzzz_0, tg_zz_xyzzzz_1, tg_zz_xzzzzz_0, tg_zz_xzzzzz_1, \
                                         tg_zz_yyyyyy_0, tg_zz_yyyyyy_1, tg_zz_yyyyyz_0, tg_zz_yyyyyz_1, tg_zz_yyyyzz_0, \
                                         tg_zz_yyyyzz_1, tg_zz_yyyzzz_0, tg_zz_yyyzzz_1, tg_zz_yyzzzz_0, tg_zz_yyzzzz_1, \
                                         tg_zz_yzzzzz_0, tg_zz_yzzzzz_1, tg_zz_zzzzzz_0, tg_zz_zzzzzz_1, tg_zzz_xxxxx_1, \
                                         tg_zzz_xxxxxx_0, tg_zzz_xxxxxx_1, tg_zzz_xxxxxy_0, tg_zzz_xxxxxy_1, tg_zzz_xxxxxz_0, \
                                         tg_zzz_xxxxxz_1, tg_zzz_xxxxy_1, tg_zzz_xxxxyy_0, tg_zzz_xxxxyy_1, tg_zzz_xxxxyz_0, \
                                         tg_zzz_xxxxyz_1, tg_zzz_xxxxz_1, tg_zzz_xxxxzz_0, tg_zzz_xxxxzz_1, tg_zzz_xxxyy_1, \
                                         tg_zzz_xxxyyy_0, tg_zzz_xxxyyy_1, tg_zzz_xxxyyz_0, tg_zzz_xxxyyz_1, tg_zzz_xxxyz_1, \
                                         tg_zzz_xxxyzz_0, tg_zzz_xxxyzz_1, tg_zzz_xxxzz_1, tg_zzz_xxxzzz_0, tg_zzz_xxxzzz_1, \
                                         tg_zzz_xxyyy_1, tg_zzz_xxyyyy_0, tg_zzz_xxyyyy_1, tg_zzz_xxyyyz_0, tg_zzz_xxyyyz_1, \
                                         tg_zzz_xxyyz_1, tg_zzz_xxyyzz_0, tg_zzz_xxyyzz_1, tg_zzz_xxyzz_1, tg_zzz_xxyzzz_0, \
                                         tg_zzz_xxyzzz_1, tg_zzz_xxzzz_1, tg_zzz_xxzzzz_0, tg_zzz_xxzzzz_1, tg_zzz_xyyyy_1, \
                                         tg_zzz_xyyyyy_0, tg_zzz_xyyyyy_1, tg_zzz_xyyyyz_0, tg_zzz_xyyyyz_1, tg_zzz_xyyyz_1, \
                                         tg_zzz_xyyyzz_0, tg_zzz_xyyyzz_1, tg_zzz_xyyzz_1, tg_zzz_xyyzzz_0, tg_zzz_xyyzzz_1, \
                                         tg_zzz_xyzzz_1, tg_zzz_xyzzzz_0, tg_zzz_xyzzzz_1, tg_zzz_xzzzz_1, tg_zzz_xzzzzz_0, \
                                         tg_zzz_xzzzzz_1, tg_zzz_yyyyy_1, tg_zzz_yyyyyy_0, tg_zzz_yyyyyy_1, tg_zzz_yyyyyz_0, \
                                         tg_zzz_yyyyyz_1, tg_zzz_yyyyz_1, tg_zzz_yyyyzz_0, tg_zzz_yyyyzz_1, tg_zzz_yyyzz_1, \
                                         tg_zzz_yyyzzz_0, tg_zzz_yyyzzz_1, tg_zzz_yyzzz_1, tg_zzz_yyzzzz_0, tg_zzz_yyzzzz_1, \
                                         tg_zzz_yzzzz_1, tg_zzz_yzzzzz_0, tg_zzz_yzzzzz_1, tg_zzz_zzzzz_1, tg_zzz_zzzzzz_0, \
                                         tg_zzz_zzzzzz_1, tg_zzzz_xxxxxx_0, tg_zzzz_xxxxxy_0, tg_zzzz_xxxxxz_0, \
                                         tg_zzzz_xxxxyy_0, tg_zzzz_xxxxyz_0, tg_zzzz_xxxxzz_0, tg_zzzz_xxxyyy_0, \
                                         tg_zzzz_xxxyyz_0, tg_zzzz_xxxyzz_0, tg_zzzz_xxxzzz_0, tg_zzzz_xxyyyy_0, \
                                         tg_zzzz_xxyyyz_0, tg_zzzz_xxyyzz_0, tg_zzzz_xxyzzz_0, tg_zzzz_xxzzzz_0, \
                                         tg_zzzz_xyyyyy_0, tg_zzzz_xyyyyz_0, tg_zzzz_xyyyzz_0, tg_zzzz_xyyzzz_0, \
                                         tg_zzzz_xyzzzz_0, tg_zzzz_xzzzzz_0, tg_zzzz_yyyyyy_0, tg_zzzz_yyyyyz_0, \
                                         tg_zzzz_yyyyzz_0, tg_zzzz_yyyzzz_0, tg_zzzz_yyzzzz_0, tg_zzzz_yzzzzz_0, \
                                         tg_zzzz_zzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyzz_xxxxxx_0[j] = pb_y * tg_yzz_xxxxxx_0[j] + fr * tg_yzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxx_0[j] - tg_zz_xxxxxx_1[j] * fl1_fza);

                    tg_yyzz_xxxxxy_0[j] = pb_y * tg_yzz_xxxxxy_0[j] + fr * tg_yzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxy_0[j] - tg_zz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxxx_1[j];

                    tg_yyzz_xxxxxz_0[j] = pb_y * tg_yzz_xxxxxz_0[j] + fr * tg_yzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxxz_0[j] - tg_zz_xxxxxz_1[j] * fl1_fza);

                    tg_yyzz_xxxxyy_0[j] = pb_y * tg_yzz_xxxxyy_0[j] + fr * tg_yzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyy_0[j] - tg_zz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxxxy_1[j];

                    tg_yyzz_xxxxyz_0[j] = pb_y * tg_yzz_xxxxyz_0[j] + fr * tg_yzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxyz_0[j] - tg_zz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxxz_1[j];

                    tg_yyzz_xxxxzz_0[j] = pb_y * tg_yzz_xxxxzz_0[j] + fr * tg_yzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxxzz_0[j] - tg_zz_xxxxzz_1[j] * fl1_fza);

                    tg_yyzz_xxxyyy_0[j] = pb_y * tg_yzz_xxxyyy_0[j] + fr * tg_yzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyy_0[j] - tg_zz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xxxyy_1[j];

                    tg_yyzz_xxxyyz_0[j] = pb_y * tg_yzz_xxxyyz_0[j] + fr * tg_yzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyyz_0[j] - tg_zz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxxyz_1[j];

                    tg_yyzz_xxxyzz_0[j] = pb_y * tg_yzz_xxxyzz_0[j] + fr * tg_yzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxyzz_0[j] - tg_zz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxxzz_1[j];

                    tg_yyzz_xxxzzz_0[j] = pb_y * tg_yzz_xxxzzz_0[j] + fr * tg_yzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxxzzz_0[j] - tg_zz_xxxzzz_1[j] * fl1_fza);

                    tg_yyzz_xxyyyy_0[j] = pb_y * tg_yzz_xxyyyy_0[j] + fr * tg_yzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyy_0[j] - tg_zz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_xxyyy_1[j];

                    tg_yyzz_xxyyyz_0[j] = pb_y * tg_yzz_xxyyyz_0[j] + fr * tg_yzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyyz_0[j] - tg_zz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xxyyz_1[j];

                    tg_yyzz_xxyyzz_0[j] = pb_y * tg_yzz_xxyyzz_0[j] + fr * tg_yzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyyzz_0[j] - tg_zz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xxyzz_1[j];

                    tg_yyzz_xxyzzz_0[j] = pb_y * tg_yzz_xxyzzz_0[j] + fr * tg_yzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxyzzz_0[j] - tg_zz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xxzzz_1[j];

                    tg_yyzz_xxzzzz_0[j] = pb_y * tg_yzz_xxzzzz_0[j] + fr * tg_yzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xxzzzz_0[j] - tg_zz_xxzzzz_1[j] * fl1_fza);

                    tg_yyzz_xyyyyy_0[j] = pb_y * tg_yzz_xyyyyy_0[j] + fr * tg_yzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyy_0[j] - tg_zz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzz_xyyyy_1[j];

                    tg_yyzz_xyyyyz_0[j] = pb_y * tg_yzz_xyyyyz_0[j] + fr * tg_yzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyyz_0[j] - tg_zz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_xyyyz_1[j];

                    tg_yyzz_xyyyzz_0[j] = pb_y * tg_yzz_xyyyzz_0[j] + fr * tg_yzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyyzz_0[j] - tg_zz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_xyyzz_1[j];

                    tg_yyzz_xyyzzz_0[j] = pb_y * tg_yzz_xyyzzz_0[j] + fr * tg_yzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyyzzz_0[j] - tg_zz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xyzzz_1[j];

                    tg_yyzz_xyzzzz_0[j] = pb_y * tg_yzz_xyzzzz_0[j] + fr * tg_yzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xyzzzz_0[j] - tg_zz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xzzzz_1[j];

                    tg_yyzz_xzzzzz_0[j] = pb_y * tg_yzz_xzzzzz_0[j] + fr * tg_yzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzzzzz_0[j] - tg_zz_xzzzzz_1[j] * fl1_fza);

                    tg_yyzz_yyyyyy_0[j] = pb_y * tg_yzz_yyyyyy_0[j] + fr * tg_yzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyy_0[j] - tg_zz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzz_yyyyy_1[j];

                    tg_yyzz_yyyyyz_0[j] = pb_y * tg_yzz_yyyyyz_0[j] + fr * tg_yzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyyz_0[j] - tg_zz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzz_yyyyz_1[j];

                    tg_yyzz_yyyyzz_0[j] = pb_y * tg_yzz_yyyyzz_0[j] + fr * tg_yzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyyzz_0[j] - tg_zz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzz_yyyzz_1[j];

                    tg_yyzz_yyyzzz_0[j] = pb_y * tg_yzz_yyyzzz_0[j] + fr * tg_yzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyyzzz_0[j] - tg_zz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_yyzzz_1[j];

                    tg_yyzz_yyzzzz_0[j] = pb_y * tg_yzz_yyzzzz_0[j] + fr * tg_yzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yyzzzz_0[j] - tg_zz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_yzzzz_1[j];

                    tg_yyzz_yzzzzz_0[j] = pb_y * tg_yzz_yzzzzz_0[j] + fr * tg_yzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzzzzz_0[j] - tg_zz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_zzzzz_1[j];

                    tg_yyzz_zzzzzz_0[j] = pb_y * tg_yzz_zzzzzz_0[j] + fr * tg_yzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzzzzz_0[j] - tg_zz_zzzzzz_1[j] * fl1_fza);

                    tg_yzzz_xxxxxx_0[j] = pb_y * tg_zzz_xxxxxx_0[j] + fr * tg_zzz_xxxxxx_1[j];

                    tg_yzzz_xxxxxy_0[j] = pb_y * tg_zzz_xxxxxy_0[j] + fr * tg_zzz_xxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxx_1[j];

                    tg_yzzz_xxxxxz_0[j] = pb_y * tg_zzz_xxxxxz_0[j] + fr * tg_zzz_xxxxxz_1[j];

                    tg_yzzz_xxxxyy_0[j] = pb_y * tg_zzz_xxxxyy_0[j] + fr * tg_zzz_xxxxyy_1[j] + fl1_fxn * tg_zzz_xxxxy_1[j];

                    tg_yzzz_xxxxyz_0[j] = pb_y * tg_zzz_xxxxyz_0[j] + fr * tg_zzz_xxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxxz_1[j];

                    tg_yzzz_xxxxzz_0[j] = pb_y * tg_zzz_xxxxzz_0[j] + fr * tg_zzz_xxxxzz_1[j];

                    tg_yzzz_xxxyyy_0[j] = pb_y * tg_zzz_xxxyyy_0[j] + fr * tg_zzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxxyy_1[j];

                    tg_yzzz_xxxyyz_0[j] = pb_y * tg_zzz_xxxyyz_0[j] + fr * tg_zzz_xxxyyz_1[j] + fl1_fxn * tg_zzz_xxxyz_1[j];

                    tg_yzzz_xxxyzz_0[j] = pb_y * tg_zzz_xxxyzz_0[j] + fr * tg_zzz_xxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxxzz_1[j];

                    tg_yzzz_xxxzzz_0[j] = pb_y * tg_zzz_xxxzzz_0[j] + fr * tg_zzz_xxxzzz_1[j];

                    tg_yzzz_xxyyyy_0[j] = pb_y * tg_zzz_xxyyyy_0[j] + fr * tg_zzz_xxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzz_xxyyy_1[j];

                    tg_yzzz_xxyyyz_0[j] = pb_y * tg_zzz_xxyyyz_0[j] + fr * tg_zzz_xxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxyyz_1[j];

                    tg_yzzz_xxyyzz_0[j] = pb_y * tg_zzz_xxyyzz_0[j] + fr * tg_zzz_xxyyzz_1[j] + fl1_fxn * tg_zzz_xxyzz_1[j];

                    tg_yzzz_xxyzzz_0[j] = pb_y * tg_zzz_xxyzzz_0[j] + fr * tg_zzz_xxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxzzz_1[j];

                    tg_yzzz_xxzzzz_0[j] = pb_y * tg_zzz_xxzzzz_0[j] + fr * tg_zzz_xxzzzz_1[j];

                    tg_yzzz_xyyyyy_0[j] = pb_y * tg_zzz_xyyyyy_0[j] + fr * tg_zzz_xyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzz_xyyyy_1[j];

                    tg_yzzz_xyyyyz_0[j] = pb_y * tg_zzz_xyyyyz_0[j] + fr * tg_zzz_xyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzz_xyyyz_1[j];

                    tg_yzzz_xyyyzz_0[j] = pb_y * tg_zzz_xyyyzz_0[j] + fr * tg_zzz_xyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xyyzz_1[j];

                    tg_yzzz_xyyzzz_0[j] = pb_y * tg_zzz_xyyzzz_0[j] + fr * tg_zzz_xyyzzz_1[j] + fl1_fxn * tg_zzz_xyzzz_1[j];

                    tg_yzzz_xyzzzz_0[j] = pb_y * tg_zzz_xyzzzz_0[j] + fr * tg_zzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xzzzz_1[j];

                    tg_yzzz_xzzzzz_0[j] = pb_y * tg_zzz_xzzzzz_0[j] + fr * tg_zzz_xzzzzz_1[j];

                    tg_yzzz_yyyyyy_0[j] = pb_y * tg_zzz_yyyyyy_0[j] + fr * tg_zzz_yyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzz_yyyyy_1[j];

                    tg_yzzz_yyyyyz_0[j] = pb_y * tg_zzz_yyyyyz_0[j] + fr * tg_zzz_yyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzz_yyyyz_1[j];

                    tg_yzzz_yyyyzz_0[j] = pb_y * tg_zzz_yyyyzz_0[j] + fr * tg_zzz_yyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzz_yyyzz_1[j];

                    tg_yzzz_yyyzzz_0[j] = pb_y * tg_zzz_yyyzzz_0[j] + fr * tg_zzz_yyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_yyzzz_1[j];

                    tg_yzzz_yyzzzz_0[j] = pb_y * tg_zzz_yyzzzz_0[j] + fr * tg_zzz_yyzzzz_1[j] + fl1_fxn * tg_zzz_yzzzz_1[j];

                    tg_yzzz_yzzzzz_0[j] = pb_y * tg_zzz_yzzzzz_0[j] + fr * tg_zzz_yzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzzzz_1[j];

                    tg_yzzz_zzzzzz_0[j] = pb_y * tg_zzz_zzzzzz_0[j] + fr * tg_zzz_zzzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzz_xxxxxx_0[j] = pb_z * tg_zzz_xxxxxx_0[j] + fr * tg_zzz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxx_0[j] - tg_zz_xxxxxx_1[j] * fl1_fza);

                    tg_zzzz_xxxxxy_0[j] = pb_z * tg_zzz_xxxxxy_0[j] + fr * tg_zzz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxy_0[j] - tg_zz_xxxxxy_1[j] * fl1_fza);

                    tg_zzzz_xxxxxz_0[j] = pb_z * tg_zzz_xxxxxz_0[j] + fr * tg_zzz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxxz_0[j] - tg_zz_xxxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxxx_1[j];

                    tg_zzzz_xxxxyy_0[j] = pb_z * tg_zzz_xxxxyy_0[j] + fr * tg_zzz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxyy_0[j] - tg_zz_xxxxyy_1[j] * fl1_fza);

                    tg_zzzz_xxxxyz_0[j] = pb_z * tg_zzz_xxxxyz_0[j] + fr * tg_zzz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxyz_0[j] - tg_zz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxxy_1[j];

                    tg_zzzz_xxxxzz_0[j] = pb_z * tg_zzz_xxxxzz_0[j] + fr * tg_zzz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxxzz_0[j] - tg_zz_xxxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxxxz_1[j];

                    tg_zzzz_xxxyyy_0[j] = pb_z * tg_zzz_xxxyyy_0[j] + fr * tg_zzz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyyy_0[j] - tg_zz_xxxyyy_1[j] * fl1_fza);

                    tg_zzzz_xxxyyz_0[j] = pb_z * tg_zzz_xxxyyz_0[j] + fr * tg_zzz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyyz_0[j] - tg_zz_xxxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxxyy_1[j];

                    tg_zzzz_xxxyzz_0[j] = pb_z * tg_zzz_xxxyzz_0[j] + fr * tg_zzz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxyzz_0[j] - tg_zz_xxxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxxyz_1[j];

                    tg_zzzz_xxxzzz_0[j] = pb_z * tg_zzz_xxxzzz_0[j] + fr * tg_zzz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxxzzz_0[j] - tg_zz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xxxzz_1[j];

                    tg_zzzz_xxyyyy_0[j] = pb_z * tg_zzz_xxyyyy_0[j] + fr * tg_zzz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyyy_0[j] - tg_zz_xxyyyy_1[j] * fl1_fza);

                    tg_zzzz_xxyyyz_0[j] = pb_z * tg_zzz_xxyyyz_0[j] + fr * tg_zzz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyyz_0[j] - tg_zz_xxyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xxyyy_1[j];

                    tg_zzzz_xxyyzz_0[j] = pb_z * tg_zzz_xxyyzz_0[j] + fr * tg_zzz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyyzz_0[j] - tg_zz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xxyyz_1[j];

                    tg_zzzz_xxyzzz_0[j] = pb_z * tg_zzz_xxyzzz_0[j] + fr * tg_zzz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxyzzz_0[j] - tg_zz_xxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xxyzz_1[j];

                    tg_zzzz_xxzzzz_0[j] = pb_z * tg_zzz_xxzzzz_0[j] + fr * tg_zzz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xxzzzz_0[j] - tg_zz_xxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_xxzzz_1[j];

                    tg_zzzz_xyyyyy_0[j] = pb_z * tg_zzz_xyyyyy_0[j] + fr * tg_zzz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyyy_0[j] - tg_zz_xyyyyy_1[j] * fl1_fza);

                    tg_zzzz_xyyyyz_0[j] = pb_z * tg_zzz_xyyyyz_0[j] + fr * tg_zzz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyyz_0[j] - tg_zz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xyyyy_1[j];

                    tg_zzzz_xyyyzz_0[j] = pb_z * tg_zzz_xyyyzz_0[j] + fr * tg_zzz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyyzz_0[j] - tg_zz_xyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xyyyz_1[j];

                    tg_zzzz_xyyzzz_0[j] = pb_z * tg_zzz_xyyzzz_0[j] + fr * tg_zzz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyyzzz_0[j] - tg_zz_xyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_xyyzz_1[j];

                    tg_zzzz_xyzzzz_0[j] = pb_z * tg_zzz_xyzzzz_0[j] + fr * tg_zzz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xyzzzz_0[j] - tg_zz_xyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_xyzzz_1[j];

                    tg_zzzz_xzzzzz_0[j] = pb_z * tg_zzz_xzzzzz_0[j] + fr * tg_zzz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_xzzzzz_0[j] - tg_zz_xzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzz_xzzzz_1[j];

                    tg_zzzz_yyyyyy_0[j] = pb_z * tg_zzz_yyyyyy_0[j] + fr * tg_zzz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyyy_0[j] - tg_zz_yyyyyy_1[j] * fl1_fza);

                    tg_zzzz_yyyyyz_0[j] = pb_z * tg_zzz_yyyyyz_0[j] + fr * tg_zzz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyyz_0[j] - tg_zz_yyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_yyyyy_1[j];

                    tg_zzzz_yyyyzz_0[j] = pb_z * tg_zzz_yyyyzz_0[j] + fr * tg_zzz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyyzz_0[j] - tg_zz_yyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_yyyyz_1[j];

                    tg_zzzz_yyyzzz_0[j] = pb_z * tg_zzz_yyyzzz_0[j] + fr * tg_zzz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyyzzz_0[j] - tg_zz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_yyyzz_1[j];

                    tg_zzzz_yyzzzz_0[j] = pb_z * tg_zzz_yyzzzz_0[j] + fr * tg_zzz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yyzzzz_0[j] - tg_zz_yyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzz_yyzzz_1[j];

                    tg_zzzz_yzzzzz_0[j] = pb_z * tg_zzz_yzzzzz_0[j] + fr * tg_zzz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_yzzzzz_0[j] - tg_zz_yzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzz_yzzzz_1[j];

                    tg_zzzz_zzzzzz_0[j] = pb_z * tg_zzz_zzzzzz_0[j] + fr * tg_zzz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_zz_zzzzzz_0[j] - tg_zz_zzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzz_zzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

