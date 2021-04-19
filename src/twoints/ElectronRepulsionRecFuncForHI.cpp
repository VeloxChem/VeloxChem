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

#include "ElectronRepulsionRecFuncForHI.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSHSI(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSHSI_0_98(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSHSI_98_196(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSHSI_196_294(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSI_294_392(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSI_392_490(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSI_490_588(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSHSI_0_98(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                // set up pointers to integrals

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

                // Batch of Integrals (0,98)

                #pragma omp simd aligned(fxn, fza, tg_xxx_xxxxxx_0, tg_xxx_xxxxxx_1, tg_xxx_xxxxxy_0, \
                                         tg_xxx_xxxxxy_1, tg_xxx_xxxxxz_0, tg_xxx_xxxxxz_1, tg_xxx_xxxxyy_0, tg_xxx_xxxxyy_1, \
                                         tg_xxx_xxxxyz_0, tg_xxx_xxxxyz_1, tg_xxx_xxxxzz_0, tg_xxx_xxxxzz_1, tg_xxx_xxxyyy_0, \
                                         tg_xxx_xxxyyy_1, tg_xxx_xxxyyz_0, tg_xxx_xxxyyz_1, tg_xxx_xxxyzz_0, tg_xxx_xxxyzz_1, \
                                         tg_xxx_xxxzzz_0, tg_xxx_xxxzzz_1, tg_xxx_xxyyyy_0, tg_xxx_xxyyyy_1, tg_xxx_xxyyyz_0, \
                                         tg_xxx_xxyyyz_1, tg_xxx_xxyyzz_0, tg_xxx_xxyyzz_1, tg_xxx_xxyzzz_0, tg_xxx_xxyzzz_1, \
                                         tg_xxx_xxzzzz_0, tg_xxx_xxzzzz_1, tg_xxx_xyyyyy_0, tg_xxx_xyyyyy_1, tg_xxx_xyyyyz_0, \
                                         tg_xxx_xyyyyz_1, tg_xxx_xyyyzz_0, tg_xxx_xyyyzz_1, tg_xxx_xyyzzz_0, tg_xxx_xyyzzz_1, \
                                         tg_xxx_xyzzzz_0, tg_xxx_xyzzzz_1, tg_xxx_xzzzzz_0, tg_xxx_xzzzzz_1, tg_xxx_yyyyyy_0, \
                                         tg_xxx_yyyyyy_1, tg_xxx_yyyyyz_0, tg_xxx_yyyyyz_1, tg_xxx_yyyyzz_0, tg_xxx_yyyyzz_1, \
                                         tg_xxx_yyyzzz_0, tg_xxx_yyyzzz_1, tg_xxx_yyzzzz_0, tg_xxx_yyzzzz_1, tg_xxx_yzzzzz_0, \
                                         tg_xxx_yzzzzz_1, tg_xxx_zzzzzz_0, tg_xxx_zzzzzz_1, tg_xxxx_xxxxx_1, \
                                         tg_xxxx_xxxxxx_0, tg_xxxx_xxxxxx_1, tg_xxxx_xxxxxy_0, tg_xxxx_xxxxxy_1, \
                                         tg_xxxx_xxxxxz_0, tg_xxxx_xxxxxz_1, tg_xxxx_xxxxy_1, tg_xxxx_xxxxyy_0, \
                                         tg_xxxx_xxxxyy_1, tg_xxxx_xxxxyz_0, tg_xxxx_xxxxyz_1, tg_xxxx_xxxxz_1, \
                                         tg_xxxx_xxxxzz_0, tg_xxxx_xxxxzz_1, tg_xxxx_xxxyy_1, tg_xxxx_xxxyyy_0, \
                                         tg_xxxx_xxxyyy_1, tg_xxxx_xxxyyz_0, tg_xxxx_xxxyyz_1, tg_xxxx_xxxyz_1, \
                                         tg_xxxx_xxxyzz_0, tg_xxxx_xxxyzz_1, tg_xxxx_xxxzz_1, tg_xxxx_xxxzzz_0, \
                                         tg_xxxx_xxxzzz_1, tg_xxxx_xxyyy_1, tg_xxxx_xxyyyy_0, tg_xxxx_xxyyyy_1, \
                                         tg_xxxx_xxyyyz_0, tg_xxxx_xxyyyz_1, tg_xxxx_xxyyz_1, tg_xxxx_xxyyzz_0, \
                                         tg_xxxx_xxyyzz_1, tg_xxxx_xxyzz_1, tg_xxxx_xxyzzz_0, tg_xxxx_xxyzzz_1, \
                                         tg_xxxx_xxzzz_1, tg_xxxx_xxzzzz_0, tg_xxxx_xxzzzz_1, tg_xxxx_xyyyy_1, \
                                         tg_xxxx_xyyyyy_0, tg_xxxx_xyyyyy_1, tg_xxxx_xyyyyz_0, tg_xxxx_xyyyyz_1, \
                                         tg_xxxx_xyyyz_1, tg_xxxx_xyyyzz_0, tg_xxxx_xyyyzz_1, tg_xxxx_xyyzz_1, \
                                         tg_xxxx_xyyzzz_0, tg_xxxx_xyyzzz_1, tg_xxxx_xyzzz_1, tg_xxxx_xyzzzz_0, \
                                         tg_xxxx_xyzzzz_1, tg_xxxx_xzzzz_1, tg_xxxx_xzzzzz_0, tg_xxxx_xzzzzz_1, \
                                         tg_xxxx_yyyyy_1, tg_xxxx_yyyyyy_0, tg_xxxx_yyyyyy_1, tg_xxxx_yyyyyz_0, \
                                         tg_xxxx_yyyyyz_1, tg_xxxx_yyyyz_1, tg_xxxx_yyyyzz_0, tg_xxxx_yyyyzz_1, \
                                         tg_xxxx_yyyzz_1, tg_xxxx_yyyzzz_0, tg_xxxx_yyyzzz_1, tg_xxxx_yyzzz_1, \
                                         tg_xxxx_yyzzzz_0, tg_xxxx_yyzzzz_1, tg_xxxx_yzzzz_1, tg_xxxx_yzzzzz_0, \
                                         tg_xxxx_yzzzzz_1, tg_xxxx_zzzzz_1, tg_xxxx_zzzzzz_0, tg_xxxx_zzzzzz_1, \
                                         tg_xxxxx_xxxxxx_0, tg_xxxxx_xxxxxy_0, tg_xxxxx_xxxxxz_0, tg_xxxxx_xxxxyy_0, \
                                         tg_xxxxx_xxxxyz_0, tg_xxxxx_xxxxzz_0, tg_xxxxx_xxxyyy_0, tg_xxxxx_xxxyyz_0, \
                                         tg_xxxxx_xxxyzz_0, tg_xxxxx_xxxzzz_0, tg_xxxxx_xxyyyy_0, tg_xxxxx_xxyyyz_0, \
                                         tg_xxxxx_xxyyzz_0, tg_xxxxx_xxyzzz_0, tg_xxxxx_xxzzzz_0, tg_xxxxx_xyyyyy_0, \
                                         tg_xxxxx_xyyyyz_0, tg_xxxxx_xyyyzz_0, tg_xxxxx_xyyzzz_0, tg_xxxxx_xyzzzz_0, \
                                         tg_xxxxx_xzzzzz_0, tg_xxxxx_yyyyyy_0, tg_xxxxx_yyyyyz_0, tg_xxxxx_yyyyzz_0, \
                                         tg_xxxxx_yyyzzz_0, tg_xxxxx_yyzzzz_0, tg_xxxxx_yzzzzz_0, tg_xxxxx_zzzzzz_0, \
                                         tg_xxxxy_xxxxxx_0, tg_xxxxy_xxxxxy_0, tg_xxxxy_xxxxxz_0, tg_xxxxy_xxxxyy_0, \
                                         tg_xxxxy_xxxxyz_0, tg_xxxxy_xxxxzz_0, tg_xxxxy_xxxyyy_0, tg_xxxxy_xxxyyz_0, \
                                         tg_xxxxy_xxxyzz_0, tg_xxxxy_xxxzzz_0, tg_xxxxy_xxyyyy_0, tg_xxxxy_xxyyyz_0, \
                                         tg_xxxxy_xxyyzz_0, tg_xxxxy_xxyzzz_0, tg_xxxxy_xxzzzz_0, tg_xxxxy_xyyyyy_0, \
                                         tg_xxxxy_xyyyyz_0, tg_xxxxy_xyyyzz_0, tg_xxxxy_xyyzzz_0, tg_xxxxy_xyzzzz_0, \
                                         tg_xxxxy_xzzzzz_0, tg_xxxxy_yyyyyy_0, tg_xxxxy_yyyyyz_0, tg_xxxxy_yyyyzz_0, \
                                         tg_xxxxy_yyyzzz_0, tg_xxxxy_yyzzzz_0, tg_xxxxy_yzzzzz_0, tg_xxxxy_zzzzzz_0, \
                                         tg_xxxxz_xxxxxx_0, tg_xxxxz_xxxxxy_0, tg_xxxxz_xxxxxz_0, tg_xxxxz_xxxxyy_0, \
                                         tg_xxxxz_xxxxyz_0, tg_xxxxz_xxxxzz_0, tg_xxxxz_xxxyyy_0, tg_xxxxz_xxxyyz_0, \
                                         tg_xxxxz_xxxyzz_0, tg_xxxxz_xxxzzz_0, tg_xxxxz_xxyyyy_0, tg_xxxxz_xxyyyz_0, \
                                         tg_xxxxz_xxyyzz_0, tg_xxxxz_xxyzzz_0, tg_xxxxz_xxzzzz_0, tg_xxxxz_xyyyyy_0, \
                                         tg_xxxxz_xyyyyz_0, tg_xxxxz_xyyyzz_0, tg_xxxxz_xyyzzz_0, tg_xxxxz_xyzzzz_0, \
                                         tg_xxxxz_xzzzzz_0, tg_xxxxz_yyyyyy_0, tg_xxxxz_yyyyyz_0, tg_xxxxz_yyyyzz_0, \
                                         tg_xxxxz_yyyzzz_0, tg_xxxxz_yyzzzz_0, tg_xxxxz_yzzzzz_0, tg_xxxxz_zzzzzz_0, \
                                         tg_xxxy_xxxxx_1, tg_xxxy_xxxxxx_0, tg_xxxy_xxxxxx_1, tg_xxxy_xxxxxy_0, \
                                         tg_xxxy_xxxxxy_1, tg_xxxy_xxxxxz_0, tg_xxxy_xxxxxz_1, tg_xxxy_xxxxy_1, \
                                         tg_xxxy_xxxxyy_0, tg_xxxy_xxxxyy_1, tg_xxxy_xxxxyz_0, tg_xxxy_xxxxyz_1, \
                                         tg_xxxy_xxxxz_1, tg_xxxy_xxxxzz_0, tg_xxxy_xxxxzz_1, tg_xxxy_xxxyy_1, \
                                         tg_xxxy_xxxyyy_0, tg_xxxy_xxxyyy_1, tg_xxxy_xxxyyz_0, tg_xxxy_xxxyyz_1, \
                                         tg_xxxy_xxxyz_1, tg_xxxy_xxxyzz_0, tg_xxxy_xxxyzz_1, tg_xxxy_xxxzz_1, \
                                         tg_xxxy_xxxzzz_0, tg_xxxy_xxxzzz_1, tg_xxxy_xxyyy_1, tg_xxxy_xxyyyy_0, \
                                         tg_xxxy_xxyyyy_1, tg_xxxy_xxyyyz_0, tg_xxxy_xxyyyz_1, tg_xxxy_xxyyz_1, \
                                         tg_xxxy_xxyyzz_0, tg_xxxy_xxyyzz_1, tg_xxxy_xxyzz_1, tg_xxxy_xxyzzz_0, \
                                         tg_xxxy_xxyzzz_1, tg_xxxy_xxzzz_1, tg_xxxy_xxzzzz_0, tg_xxxy_xxzzzz_1, \
                                         tg_xxxy_xyyyy_1, tg_xxxy_xyyyyy_0, tg_xxxy_xyyyyy_1, tg_xxxy_xyyyyz_0, \
                                         tg_xxxy_xyyyyz_1, tg_xxxy_xyyyz_1, tg_xxxy_xyyyzz_0, tg_xxxy_xyyyzz_1, \
                                         tg_xxxy_xyyzz_1, tg_xxxy_xyyzzz_0, tg_xxxy_xyyzzz_1, tg_xxxy_xyzzz_1, \
                                         tg_xxxy_xyzzzz_0, tg_xxxy_xyzzzz_1, tg_xxxy_xzzzz_1, tg_xxxy_xzzzzz_0, \
                                         tg_xxxy_xzzzzz_1, tg_xxxy_yyyyy_1, tg_xxxy_yyyyyy_0, tg_xxxy_yyyyyy_1, \
                                         tg_xxxy_yyyyyz_0, tg_xxxy_yyyyyz_1, tg_xxxy_yyyyz_1, tg_xxxy_yyyyzz_0, \
                                         tg_xxxy_yyyyzz_1, tg_xxxy_yyyzz_1, tg_xxxy_yyyzzz_0, tg_xxxy_yyyzzz_1, \
                                         tg_xxxy_yyzzz_1, tg_xxxy_yyzzzz_0, tg_xxxy_yyzzzz_1, tg_xxxy_yzzzz_1, \
                                         tg_xxxy_yzzzzz_0, tg_xxxy_yzzzzz_1, tg_xxxy_zzzzz_1, tg_xxxy_zzzzzz_0, \
                                         tg_xxxy_zzzzzz_1, tg_xxxyy_xxxxxx_0, tg_xxxyy_xxxxxy_0, tg_xxxyy_xxxxxz_0, \
                                         tg_xxxyy_xxxxyy_0, tg_xxxyy_xxxxyz_0, tg_xxxyy_xxxxzz_0, tg_xxxyy_xxxyyy_0, \
                                         tg_xxxyy_xxxyyz_0, tg_xxxyy_xxxyzz_0, tg_xxxyy_xxxzzz_0, tg_xxxyy_xxyyyy_0, \
                                         tg_xxxyy_xxyyyz_0, tg_xxxyy_xxyyzz_0, tg_xxxyy_xxyzzz_0, tg_xxxz_xxxxx_1, \
                                         tg_xxxz_xxxxxx_0, tg_xxxz_xxxxxx_1, tg_xxxz_xxxxxy_0, tg_xxxz_xxxxxy_1, \
                                         tg_xxxz_xxxxxz_0, tg_xxxz_xxxxxz_1, tg_xxxz_xxxxy_1, tg_xxxz_xxxxyy_0, \
                                         tg_xxxz_xxxxyy_1, tg_xxxz_xxxxyz_0, tg_xxxz_xxxxyz_1, tg_xxxz_xxxxz_1, \
                                         tg_xxxz_xxxxzz_0, tg_xxxz_xxxxzz_1, tg_xxxz_xxxyy_1, tg_xxxz_xxxyyy_0, \
                                         tg_xxxz_xxxyyy_1, tg_xxxz_xxxyyz_0, tg_xxxz_xxxyyz_1, tg_xxxz_xxxyz_1, \
                                         tg_xxxz_xxxyzz_0, tg_xxxz_xxxyzz_1, tg_xxxz_xxxzz_1, tg_xxxz_xxxzzz_0, \
                                         tg_xxxz_xxxzzz_1, tg_xxxz_xxyyy_1, tg_xxxz_xxyyyy_0, tg_xxxz_xxyyyy_1, \
                                         tg_xxxz_xxyyyz_0, tg_xxxz_xxyyyz_1, tg_xxxz_xxyyz_1, tg_xxxz_xxyyzz_0, \
                                         tg_xxxz_xxyyzz_1, tg_xxxz_xxyzz_1, tg_xxxz_xxyzzz_0, tg_xxxz_xxyzzz_1, \
                                         tg_xxxz_xxzzz_1, tg_xxxz_xxzzzz_0, tg_xxxz_xxzzzz_1, tg_xxxz_xyyyy_1, \
                                         tg_xxxz_xyyyyy_0, tg_xxxz_xyyyyy_1, tg_xxxz_xyyyyz_0, tg_xxxz_xyyyyz_1, \
                                         tg_xxxz_xyyyz_1, tg_xxxz_xyyyzz_0, tg_xxxz_xyyyzz_1, tg_xxxz_xyyzz_1, \
                                         tg_xxxz_xyyzzz_0, tg_xxxz_xyyzzz_1, tg_xxxz_xyzzz_1, tg_xxxz_xyzzzz_0, \
                                         tg_xxxz_xyzzzz_1, tg_xxxz_xzzzz_1, tg_xxxz_xzzzzz_0, tg_xxxz_xzzzzz_1, \
                                         tg_xxxz_yyyyy_1, tg_xxxz_yyyyyy_0, tg_xxxz_yyyyyy_1, tg_xxxz_yyyyyz_0, \
                                         tg_xxxz_yyyyyz_1, tg_xxxz_yyyyz_1, tg_xxxz_yyyyzz_0, tg_xxxz_yyyyzz_1, \
                                         tg_xxxz_yyyzz_1, tg_xxxz_yyyzzz_0, tg_xxxz_yyyzzz_1, tg_xxxz_yyzzz_1, \
                                         tg_xxxz_yyzzzz_0, tg_xxxz_yyzzzz_1, tg_xxxz_yzzzz_1, tg_xxxz_yzzzzz_0, \
                                         tg_xxxz_yzzzzz_1, tg_xxxz_zzzzz_1, tg_xxxz_zzzzzz_0, tg_xxxz_zzzzzz_1, \
                                         tg_xxy_xxxxxx_0, tg_xxy_xxxxxx_1, tg_xxy_xxxxxy_0, tg_xxy_xxxxxy_1, tg_xxy_xxxxxz_0, \
                                         tg_xxy_xxxxxz_1, tg_xxy_xxxxyy_0, tg_xxy_xxxxyy_1, tg_xxy_xxxxyz_0, tg_xxy_xxxxyz_1, \
                                         tg_xxy_xxxxzz_0, tg_xxy_xxxxzz_1, tg_xxy_xxxyyy_0, tg_xxy_xxxyyy_1, tg_xxy_xxxyyz_0, \
                                         tg_xxy_xxxyyz_1, tg_xxy_xxxyzz_0, tg_xxy_xxxyzz_1, tg_xxy_xxxzzz_0, tg_xxy_xxxzzz_1, \
                                         tg_xxy_xxyyyy_0, tg_xxy_xxyyyy_1, tg_xxy_xxyyyz_0, tg_xxy_xxyyyz_1, tg_xxy_xxyyzz_0, \
                                         tg_xxy_xxyyzz_1, tg_xxy_xxyzzz_0, tg_xxy_xxyzzz_1, tg_xxy_xxzzzz_0, tg_xxy_xxzzzz_1, \
                                         tg_xxy_xyyyyy_0, tg_xxy_xyyyyy_1, tg_xxy_xyyyyz_0, tg_xxy_xyyyyz_1, tg_xxy_xyyyzz_0, \
                                         tg_xxy_xyyyzz_1, tg_xxy_xyyzzz_0, tg_xxy_xyyzzz_1, tg_xxy_xyzzzz_0, tg_xxy_xyzzzz_1, \
                                         tg_xxy_xzzzzz_0, tg_xxy_xzzzzz_1, tg_xxy_yyyyyy_0, tg_xxy_yyyyyy_1, tg_xxy_yyyyyz_0, \
                                         tg_xxy_yyyyyz_1, tg_xxy_yyyyzz_0, tg_xxy_yyyyzz_1, tg_xxy_yyyzzz_0, tg_xxy_yyyzzz_1, \
                                         tg_xxy_yyzzzz_0, tg_xxy_yyzzzz_1, tg_xxy_yzzzzz_0, tg_xxy_yzzzzz_1, tg_xxy_zzzzzz_0, \
                                         tg_xxy_zzzzzz_1, tg_xxyy_xxxxx_1, tg_xxyy_xxxxxx_0, tg_xxyy_xxxxxx_1, \
                                         tg_xxyy_xxxxxy_0, tg_xxyy_xxxxxy_1, tg_xxyy_xxxxxz_0, tg_xxyy_xxxxxz_1, \
                                         tg_xxyy_xxxxy_1, tg_xxyy_xxxxyy_0, tg_xxyy_xxxxyy_1, tg_xxyy_xxxxyz_0, \
                                         tg_xxyy_xxxxyz_1, tg_xxyy_xxxxz_1, tg_xxyy_xxxxzz_0, tg_xxyy_xxxxzz_1, \
                                         tg_xxyy_xxxyy_1, tg_xxyy_xxxyyy_0, tg_xxyy_xxxyyy_1, tg_xxyy_xxxyyz_0, \
                                         tg_xxyy_xxxyyz_1, tg_xxyy_xxxyz_1, tg_xxyy_xxxyzz_0, tg_xxyy_xxxyzz_1, \
                                         tg_xxyy_xxxzz_1, tg_xxyy_xxxzzz_0, tg_xxyy_xxxzzz_1, tg_xxyy_xxyyy_1, \
                                         tg_xxyy_xxyyyy_0, tg_xxyy_xxyyyy_1, tg_xxyy_xxyyyz_0, tg_xxyy_xxyyyz_1, \
                                         tg_xxyy_xxyyz_1, tg_xxyy_xxyyzz_0, tg_xxyy_xxyyzz_1, tg_xxyy_xxyzz_1, \
                                         tg_xxyy_xxyzzz_0, tg_xxyy_xxyzzz_1, tg_xxyy_xxzzz_1, tg_xxyy_xyyyy_1, \
                                         tg_xxyy_xyyyz_1, tg_xxyy_xyyzz_1, tg_xxyy_xyzzz_1, tg_xxz_xxxxxx_0, tg_xxz_xxxxxx_1, \
                                         tg_xxz_xxxxxy_0, tg_xxz_xxxxxy_1, tg_xxz_xxxxxz_0, tg_xxz_xxxxxz_1, tg_xxz_xxxxyy_0, \
                                         tg_xxz_xxxxyy_1, tg_xxz_xxxxyz_0, tg_xxz_xxxxyz_1, tg_xxz_xxxxzz_0, tg_xxz_xxxxzz_1, \
                                         tg_xxz_xxxyyy_0, tg_xxz_xxxyyy_1, tg_xxz_xxxyyz_0, tg_xxz_xxxyyz_1, tg_xxz_xxxyzz_0, \
                                         tg_xxz_xxxyzz_1, tg_xxz_xxxzzz_0, tg_xxz_xxxzzz_1, tg_xxz_xxyyyy_0, tg_xxz_xxyyyy_1, \
                                         tg_xxz_xxyyyz_0, tg_xxz_xxyyyz_1, tg_xxz_xxyyzz_0, tg_xxz_xxyyzz_1, tg_xxz_xxyzzz_0, \
                                         tg_xxz_xxyzzz_1, tg_xxz_xxzzzz_0, tg_xxz_xxzzzz_1, tg_xxz_xyyyyy_0, tg_xxz_xyyyyy_1, \
                                         tg_xxz_xyyyyz_0, tg_xxz_xyyyyz_1, tg_xxz_xyyyzz_0, tg_xxz_xyyyzz_1, tg_xxz_xyyzzz_0, \
                                         tg_xxz_xyyzzz_1, tg_xxz_xyzzzz_0, tg_xxz_xyzzzz_1, tg_xxz_xzzzzz_0, tg_xxz_xzzzzz_1, \
                                         tg_xxz_yyyyyy_0, tg_xxz_yyyyyy_1, tg_xxz_yyyyyz_0, tg_xxz_yyyyyz_1, tg_xxz_yyyyzz_0, \
                                         tg_xxz_yyyyzz_1, tg_xxz_yyyzzz_0, tg_xxz_yyyzzz_1, tg_xxz_yyzzzz_0, tg_xxz_yyzzzz_1, \
                                         tg_xxz_yzzzzz_0, tg_xxz_yzzzzz_1, tg_xxz_zzzzzz_0, tg_xxz_zzzzzz_1, tg_xyy_xxxxxx_0, \
                                         tg_xyy_xxxxxx_1, tg_xyy_xxxxxy_0, tg_xyy_xxxxxy_1, tg_xyy_xxxxxz_0, tg_xyy_xxxxxz_1, \
                                         tg_xyy_xxxxyy_0, tg_xyy_xxxxyy_1, tg_xyy_xxxxyz_0, tg_xyy_xxxxyz_1, tg_xyy_xxxxzz_0, \
                                         tg_xyy_xxxxzz_1, tg_xyy_xxxyyy_0, tg_xyy_xxxyyy_1, tg_xyy_xxxyyz_0, tg_xyy_xxxyyz_1, \
                                         tg_xyy_xxxyzz_0, tg_xyy_xxxyzz_1, tg_xyy_xxxzzz_0, tg_xyy_xxxzzz_1, tg_xyy_xxyyyy_0, \
                                         tg_xyy_xxyyyy_1, tg_xyy_xxyyyz_0, tg_xyy_xxyyyz_1, tg_xyy_xxyyzz_0, tg_xyy_xxyyzz_1, \
                                         tg_xyy_xxyzzz_0, tg_xyy_xxyzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxx_xxxxxx_0[j] = pb_x * tg_xxxx_xxxxxx_0[j] + fr * tg_xxxx_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxx_0[j] - tg_xxx_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxx_xxxxx_1[j];

                    tg_xxxxx_xxxxxy_0[j] = pb_x * tg_xxxx_xxxxxy_0[j] + fr * tg_xxxx_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxy_0[j] - tg_xxx_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxx_xxxxy_1[j];

                    tg_xxxxx_xxxxxz_0[j] = pb_x * tg_xxxx_xxxxxz_0[j] + fr * tg_xxxx_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxxz_0[j] - tg_xxx_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxx_xxxxz_1[j];

                    tg_xxxxx_xxxxyy_0[j] = pb_x * tg_xxxx_xxxxyy_0[j] + fr * tg_xxxx_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxyy_0[j] - tg_xxx_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxyy_1[j];

                    tg_xxxxx_xxxxyz_0[j] = pb_x * tg_xxxx_xxxxyz_0[j] + fr * tg_xxxx_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxyz_0[j] - tg_xxx_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxyz_1[j];

                    tg_xxxxx_xxxxzz_0[j] = pb_x * tg_xxxx_xxxxzz_0[j] + fr * tg_xxxx_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxxzz_0[j] - tg_xxx_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxx_xxxzz_1[j];

                    tg_xxxxx_xxxyyy_0[j] = pb_x * tg_xxxx_xxxyyy_0[j] + fr * tg_xxxx_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyyy_0[j] - tg_xxx_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyyy_1[j];

                    tg_xxxxx_xxxyyz_0[j] = pb_x * tg_xxxx_xxxyyz_0[j] + fr * tg_xxxx_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyyz_0[j] - tg_xxx_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyyz_1[j];

                    tg_xxxxx_xxxyzz_0[j] = pb_x * tg_xxxx_xxxyzz_0[j] + fr * tg_xxxx_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxyzz_0[j] - tg_xxx_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxyzz_1[j];

                    tg_xxxxx_xxxzzz_0[j] = pb_x * tg_xxxx_xxxzzz_0[j] + fr * tg_xxxx_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxxzzz_0[j] - tg_xxx_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxx_xxzzz_1[j];

                    tg_xxxxx_xxyyyy_0[j] = pb_x * tg_xxxx_xxyyyy_0[j] + fr * tg_xxxx_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyyy_0[j] - tg_xxx_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyyy_1[j];

                    tg_xxxxx_xxyyyz_0[j] = pb_x * tg_xxxx_xxyyyz_0[j] + fr * tg_xxxx_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyyz_0[j] - tg_xxx_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyyz_1[j];

                    tg_xxxxx_xxyyzz_0[j] = pb_x * tg_xxxx_xxyyzz_0[j] + fr * tg_xxxx_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyyzz_0[j] - tg_xxx_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyyzz_1[j];

                    tg_xxxxx_xxyzzz_0[j] = pb_x * tg_xxxx_xxyzzz_0[j] + fr * tg_xxxx_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxyzzz_0[j] - tg_xxx_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xyzzz_1[j];

                    tg_xxxxx_xxzzzz_0[j] = pb_x * tg_xxxx_xxzzzz_0[j] + fr * tg_xxxx_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xxzzzz_0[j] - tg_xxx_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxx_xzzzz_1[j];

                    tg_xxxxx_xyyyyy_0[j] = pb_x * tg_xxxx_xyyyyy_0[j] + fr * tg_xxxx_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyyy_0[j] - tg_xxx_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyyy_1[j];

                    tg_xxxxx_xyyyyz_0[j] = pb_x * tg_xxxx_xyyyyz_0[j] + fr * tg_xxxx_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyyz_0[j] - tg_xxx_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyyz_1[j];

                    tg_xxxxx_xyyyzz_0[j] = pb_x * tg_xxxx_xyyyzz_0[j] + fr * tg_xxxx_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyyzz_0[j] - tg_xxx_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyyzz_1[j];

                    tg_xxxxx_xyyzzz_0[j] = pb_x * tg_xxxx_xyyzzz_0[j] + fr * tg_xxxx_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyyzzz_0[j] - tg_xxx_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yyzzz_1[j];

                    tg_xxxxx_xyzzzz_0[j] = pb_x * tg_xxxx_xyzzzz_0[j] + fr * tg_xxxx_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xyzzzz_0[j] - tg_xxx_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_yzzzz_1[j];

                    tg_xxxxx_xzzzzz_0[j] = pb_x * tg_xxxx_xzzzzz_0[j] + fr * tg_xxxx_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_xzzzzz_0[j] - tg_xxx_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxx_zzzzz_1[j];

                    tg_xxxxx_yyyyyy_0[j] = pb_x * tg_xxxx_yyyyyy_0[j] + fr * tg_xxxx_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyyy_0[j] - tg_xxx_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxx_yyyyyz_0[j] = pb_x * tg_xxxx_yyyyyz_0[j] + fr * tg_xxxx_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyyz_0[j] - tg_xxx_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxx_yyyyzz_0[j] = pb_x * tg_xxxx_yyyyzz_0[j] + fr * tg_xxxx_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyyzz_0[j] - tg_xxx_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxx_yyyzzz_0[j] = pb_x * tg_xxxx_yyyzzz_0[j] + fr * tg_xxxx_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyyzzz_0[j] - tg_xxx_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxx_yyzzzz_0[j] = pb_x * tg_xxxx_yyzzzz_0[j] + fr * tg_xxxx_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yyzzzz_0[j] - tg_xxx_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxx_yzzzzz_0[j] = pb_x * tg_xxxx_yzzzzz_0[j] + fr * tg_xxxx_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_yzzzzz_0[j] - tg_xxx_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxx_zzzzzz_0[j] = pb_x * tg_xxxx_zzzzzz_0[j] + fr * tg_xxxx_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_xxx_zzzzzz_0[j] - tg_xxx_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxy_xxxxxx_0[j] = pb_x * tg_xxxy_xxxxxx_0[j] + fr * tg_xxxy_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxx_0[j] - tg_xxy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxy_xxxxx_1[j];

                    tg_xxxxy_xxxxxy_0[j] = pb_x * tg_xxxy_xxxxxy_0[j] + fr * tg_xxxy_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxy_0[j] - tg_xxy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxy_xxxxy_1[j];

                    tg_xxxxy_xxxxxz_0[j] = pb_x * tg_xxxy_xxxxxz_0[j] + fr * tg_xxxy_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxxz_0[j] - tg_xxy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxy_xxxxz_1[j];

                    tg_xxxxy_xxxxyy_0[j] = pb_x * tg_xxxy_xxxxyy_0[j] + fr * tg_xxxy_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxyy_0[j] - tg_xxy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxyy_1[j];

                    tg_xxxxy_xxxxyz_0[j] = pb_x * tg_xxxy_xxxxyz_0[j] + fr * tg_xxxy_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxyz_0[j] - tg_xxy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxyz_1[j];

                    tg_xxxxy_xxxxzz_0[j] = pb_x * tg_xxxy_xxxxzz_0[j] + fr * tg_xxxy_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxxzz_0[j] - tg_xxy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxy_xxxzz_1[j];

                    tg_xxxxy_xxxyyy_0[j] = pb_x * tg_xxxy_xxxyyy_0[j] + fr * tg_xxxy_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyyy_0[j] - tg_xxy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyyy_1[j];

                    tg_xxxxy_xxxyyz_0[j] = pb_x * tg_xxxy_xxxyyz_0[j] + fr * tg_xxxy_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyyz_0[j] - tg_xxy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyyz_1[j];

                    tg_xxxxy_xxxyzz_0[j] = pb_x * tg_xxxy_xxxyzz_0[j] + fr * tg_xxxy_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxyzz_0[j] - tg_xxy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxyzz_1[j];

                    tg_xxxxy_xxxzzz_0[j] = pb_x * tg_xxxy_xxxzzz_0[j] + fr * tg_xxxy_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxxzzz_0[j] - tg_xxy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxy_xxzzz_1[j];

                    tg_xxxxy_xxyyyy_0[j] = pb_x * tg_xxxy_xxyyyy_0[j] + fr * tg_xxxy_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyyy_0[j] - tg_xxy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyyy_1[j];

                    tg_xxxxy_xxyyyz_0[j] = pb_x * tg_xxxy_xxyyyz_0[j] + fr * tg_xxxy_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyyz_0[j] - tg_xxy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyyz_1[j];

                    tg_xxxxy_xxyyzz_0[j] = pb_x * tg_xxxy_xxyyzz_0[j] + fr * tg_xxxy_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyyzz_0[j] - tg_xxy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyyzz_1[j];

                    tg_xxxxy_xxyzzz_0[j] = pb_x * tg_xxxy_xxyzzz_0[j] + fr * tg_xxxy_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxyzzz_0[j] - tg_xxy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xyzzz_1[j];

                    tg_xxxxy_xxzzzz_0[j] = pb_x * tg_xxxy_xxzzzz_0[j] + fr * tg_xxxy_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xxzzzz_0[j] - tg_xxy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxy_xzzzz_1[j];

                    tg_xxxxy_xyyyyy_0[j] = pb_x * tg_xxxy_xyyyyy_0[j] + fr * tg_xxxy_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyyy_0[j] - tg_xxy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyyy_1[j];

                    tg_xxxxy_xyyyyz_0[j] = pb_x * tg_xxxy_xyyyyz_0[j] + fr * tg_xxxy_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyyz_0[j] - tg_xxy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyyz_1[j];

                    tg_xxxxy_xyyyzz_0[j] = pb_x * tg_xxxy_xyyyzz_0[j] + fr * tg_xxxy_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyyzz_0[j] - tg_xxy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyyzz_1[j];

                    tg_xxxxy_xyyzzz_0[j] = pb_x * tg_xxxy_xyyzzz_0[j] + fr * tg_xxxy_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyyzzz_0[j] - tg_xxy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yyzzz_1[j];

                    tg_xxxxy_xyzzzz_0[j] = pb_x * tg_xxxy_xyzzzz_0[j] + fr * tg_xxxy_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xyzzzz_0[j] - tg_xxy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_yzzzz_1[j];

                    tg_xxxxy_xzzzzz_0[j] = pb_x * tg_xxxy_xzzzzz_0[j] + fr * tg_xxxy_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_xzzzzz_0[j] - tg_xxy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxy_zzzzz_1[j];

                    tg_xxxxy_yyyyyy_0[j] = pb_x * tg_xxxy_yyyyyy_0[j] + fr * tg_xxxy_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyyy_0[j] - tg_xxy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxy_yyyyyz_0[j] = pb_x * tg_xxxy_yyyyyz_0[j] + fr * tg_xxxy_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyyz_0[j] - tg_xxy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxy_yyyyzz_0[j] = pb_x * tg_xxxy_yyyyzz_0[j] + fr * tg_xxxy_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyyzz_0[j] - tg_xxy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxy_yyyzzz_0[j] = pb_x * tg_xxxy_yyyzzz_0[j] + fr * tg_xxxy_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyyzzz_0[j] - tg_xxy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxy_yyzzzz_0[j] = pb_x * tg_xxxy_yyzzzz_0[j] + fr * tg_xxxy_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yyzzzz_0[j] - tg_xxy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxy_yzzzzz_0[j] = pb_x * tg_xxxy_yzzzzz_0[j] + fr * tg_xxxy_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_yzzzzz_0[j] - tg_xxy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxy_zzzzzz_0[j] = pb_x * tg_xxxy_zzzzzz_0[j] + fr * tg_xxxy_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxy_zzzzzz_0[j] - tg_xxy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxxz_xxxxxx_0[j] = pb_x * tg_xxxz_xxxxxx_0[j] + fr * tg_xxxz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxx_0[j] - tg_xxz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxxz_xxxxx_1[j];

                    tg_xxxxz_xxxxxy_0[j] = pb_x * tg_xxxz_xxxxxy_0[j] + fr * tg_xxxz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxy_0[j] - tg_xxz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxz_xxxxy_1[j];

                    tg_xxxxz_xxxxxz_0[j] = pb_x * tg_xxxz_xxxxxz_0[j] + fr * tg_xxxz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxxz_0[j] - tg_xxz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxxz_xxxxz_1[j];

                    tg_xxxxz_xxxxyy_0[j] = pb_x * tg_xxxz_xxxxyy_0[j] + fr * tg_xxxz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxyy_0[j] - tg_xxz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxyy_1[j];

                    tg_xxxxz_xxxxyz_0[j] = pb_x * tg_xxxz_xxxxyz_0[j] + fr * tg_xxxz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxyz_0[j] - tg_xxz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxyz_1[j];

                    tg_xxxxz_xxxxzz_0[j] = pb_x * tg_xxxz_xxxxzz_0[j] + fr * tg_xxxz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxxzz_0[j] - tg_xxz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxz_xxxzz_1[j];

                    tg_xxxxz_xxxyyy_0[j] = pb_x * tg_xxxz_xxxyyy_0[j] + fr * tg_xxxz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyyy_0[j] - tg_xxz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyyy_1[j];

                    tg_xxxxz_xxxyyz_0[j] = pb_x * tg_xxxz_xxxyyz_0[j] + fr * tg_xxxz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyyz_0[j] - tg_xxz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyyz_1[j];

                    tg_xxxxz_xxxyzz_0[j] = pb_x * tg_xxxz_xxxyzz_0[j] + fr * tg_xxxz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxyzz_0[j] - tg_xxz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxyzz_1[j];

                    tg_xxxxz_xxxzzz_0[j] = pb_x * tg_xxxz_xxxzzz_0[j] + fr * tg_xxxz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxxzzz_0[j] - tg_xxz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxz_xxzzz_1[j];

                    tg_xxxxz_xxyyyy_0[j] = pb_x * tg_xxxz_xxyyyy_0[j] + fr * tg_xxxz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyyy_0[j] - tg_xxz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyyy_1[j];

                    tg_xxxxz_xxyyyz_0[j] = pb_x * tg_xxxz_xxyyyz_0[j] + fr * tg_xxxz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyyz_0[j] - tg_xxz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyyz_1[j];

                    tg_xxxxz_xxyyzz_0[j] = pb_x * tg_xxxz_xxyyzz_0[j] + fr * tg_xxxz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyyzz_0[j] - tg_xxz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyyzz_1[j];

                    tg_xxxxz_xxyzzz_0[j] = pb_x * tg_xxxz_xxyzzz_0[j] + fr * tg_xxxz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxyzzz_0[j] - tg_xxz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xyzzz_1[j];

                    tg_xxxxz_xxzzzz_0[j] = pb_x * tg_xxxz_xxzzzz_0[j] + fr * tg_xxxz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xxzzzz_0[j] - tg_xxz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxz_xzzzz_1[j];

                    tg_xxxxz_xyyyyy_0[j] = pb_x * tg_xxxz_xyyyyy_0[j] + fr * tg_xxxz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyyy_0[j] - tg_xxz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyyy_1[j];

                    tg_xxxxz_xyyyyz_0[j] = pb_x * tg_xxxz_xyyyyz_0[j] + fr * tg_xxxz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyyz_0[j] - tg_xxz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyyz_1[j];

                    tg_xxxxz_xyyyzz_0[j] = pb_x * tg_xxxz_xyyyzz_0[j] + fr * tg_xxxz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyyzz_0[j] - tg_xxz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyyzz_1[j];

                    tg_xxxxz_xyyzzz_0[j] = pb_x * tg_xxxz_xyyzzz_0[j] + fr * tg_xxxz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyyzzz_0[j] - tg_xxz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yyzzz_1[j];

                    tg_xxxxz_xyzzzz_0[j] = pb_x * tg_xxxz_xyzzzz_0[j] + fr * tg_xxxz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xyzzzz_0[j] - tg_xxz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_yzzzz_1[j];

                    tg_xxxxz_xzzzzz_0[j] = pb_x * tg_xxxz_xzzzzz_0[j] + fr * tg_xxxz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_xzzzzz_0[j] - tg_xxz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxz_zzzzz_1[j];

                    tg_xxxxz_yyyyyy_0[j] = pb_x * tg_xxxz_yyyyyy_0[j] + fr * tg_xxxz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyyy_0[j] - tg_xxz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxxz_yyyyyz_0[j] = pb_x * tg_xxxz_yyyyyz_0[j] + fr * tg_xxxz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyyz_0[j] - tg_xxz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxxz_yyyyzz_0[j] = pb_x * tg_xxxz_yyyyzz_0[j] + fr * tg_xxxz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyyzz_0[j] - tg_xxz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxxz_yyyzzz_0[j] = pb_x * tg_xxxz_yyyzzz_0[j] + fr * tg_xxxz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyyzzz_0[j] - tg_xxz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxxz_yyzzzz_0[j] = pb_x * tg_xxxz_yyzzzz_0[j] + fr * tg_xxxz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yyzzzz_0[j] - tg_xxz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxxz_yzzzzz_0[j] = pb_x * tg_xxxz_yzzzzz_0[j] + fr * tg_xxxz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_yzzzzz_0[j] - tg_xxz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxxz_zzzzzz_0[j] = pb_x * tg_xxxz_zzzzzz_0[j] + fr * tg_xxxz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_xxz_zzzzzz_0[j] - tg_xxz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyy_xxxxxx_0[j] = pb_x * tg_xxyy_xxxxxx_0[j] + fr * tg_xxyy_xxxxxx_1[j] + fl1_fx * (tg_xyy_xxxxxx_0[j] - tg_xyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyy_xxxxx_1[j];

                    tg_xxxyy_xxxxxy_0[j] = pb_x * tg_xxyy_xxxxxy_0[j] + fr * tg_xxyy_xxxxxy_1[j] + fl1_fx * (tg_xyy_xxxxxy_0[j] - tg_xyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyy_xxxxy_1[j];

                    tg_xxxyy_xxxxxz_0[j] = pb_x * tg_xxyy_xxxxxz_0[j] + fr * tg_xxyy_xxxxxz_1[j] + fl1_fx * (tg_xyy_xxxxxz_0[j] - tg_xyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyy_xxxxz_1[j];

                    tg_xxxyy_xxxxyy_0[j] = pb_x * tg_xxyy_xxxxyy_0[j] + fr * tg_xxyy_xxxxyy_1[j] + fl1_fx * (tg_xyy_xxxxyy_0[j] - tg_xyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxyy_1[j];

                    tg_xxxyy_xxxxyz_0[j] = pb_x * tg_xxyy_xxxxyz_0[j] + fr * tg_xxyy_xxxxyz_1[j] + fl1_fx * (tg_xyy_xxxxyz_0[j] - tg_xyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxyz_1[j];

                    tg_xxxyy_xxxxzz_0[j] = pb_x * tg_xxyy_xxxxzz_0[j] + fr * tg_xxyy_xxxxzz_1[j] + fl1_fx * (tg_xyy_xxxxzz_0[j] - tg_xyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyy_xxxzz_1[j];

                    tg_xxxyy_xxxyyy_0[j] = pb_x * tg_xxyy_xxxyyy_0[j] + fr * tg_xxyy_xxxyyy_1[j] + fl1_fx * (tg_xyy_xxxyyy_0[j] - tg_xyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyyy_1[j];

                    tg_xxxyy_xxxyyz_0[j] = pb_x * tg_xxyy_xxxyyz_0[j] + fr * tg_xxyy_xxxyyz_1[j] + fl1_fx * (tg_xyy_xxxyyz_0[j] - tg_xyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyyz_1[j];

                    tg_xxxyy_xxxyzz_0[j] = pb_x * tg_xxyy_xxxyzz_0[j] + fr * tg_xxyy_xxxyzz_1[j] + fl1_fx * (tg_xyy_xxxyzz_0[j] - tg_xyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxyzz_1[j];

                    tg_xxxyy_xxxzzz_0[j] = pb_x * tg_xxyy_xxxzzz_0[j] + fr * tg_xxyy_xxxzzz_1[j] + fl1_fx * (tg_xyy_xxxzzz_0[j] - tg_xyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyy_xxzzz_1[j];

                    tg_xxxyy_xxyyyy_0[j] = pb_x * tg_xxyy_xxyyyy_0[j] + fr * tg_xxyy_xxyyyy_1[j] + fl1_fx * (tg_xyy_xxyyyy_0[j] - tg_xyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyyy_1[j];

                    tg_xxxyy_xxyyyz_0[j] = pb_x * tg_xxyy_xxyyyz_0[j] + fr * tg_xxyy_xxyyyz_1[j] + fl1_fx * (tg_xyy_xxyyyz_0[j] - tg_xyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyyz_1[j];

                    tg_xxxyy_xxyyzz_0[j] = pb_x * tg_xxyy_xxyyzz_0[j] + fr * tg_xxyy_xxyyzz_1[j] + fl1_fx * (tg_xyy_xxyyzz_0[j] - tg_xyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyyzz_1[j];

                    tg_xxxyy_xxyzzz_0[j] = pb_x * tg_xxyy_xxyzzz_0[j] + fr * tg_xxyy_xxyzzz_1[j] + fl1_fx * (tg_xyy_xxyzzz_0[j] - tg_xyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xyzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSI_98_196(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                // set up pointers to integrals

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

                // Batch of Integrals (98,196)

                #pragma omp simd aligned(fxn, fza, tg_xxxyy_xxzzzz_0, tg_xxxyy_xyyyyy_0, tg_xxxyy_xyyyyz_0, \
                                         tg_xxxyy_xyyyzz_0, tg_xxxyy_xyyzzz_0, tg_xxxyy_xyzzzz_0, tg_xxxyy_xzzzzz_0, \
                                         tg_xxxyy_yyyyyy_0, tg_xxxyy_yyyyyz_0, tg_xxxyy_yyyyzz_0, tg_xxxyy_yyyzzz_0, \
                                         tg_xxxyy_yyzzzz_0, tg_xxxyy_yzzzzz_0, tg_xxxyy_zzzzzz_0, tg_xxxyz_xxxxxx_0, \
                                         tg_xxxyz_xxxxxy_0, tg_xxxyz_xxxxxz_0, tg_xxxyz_xxxxyy_0, tg_xxxyz_xxxxyz_0, \
                                         tg_xxxyz_xxxxzz_0, tg_xxxyz_xxxyyy_0, tg_xxxyz_xxxyyz_0, tg_xxxyz_xxxyzz_0, \
                                         tg_xxxyz_xxxzzz_0, tg_xxxyz_xxyyyy_0, tg_xxxyz_xxyyyz_0, tg_xxxyz_xxyyzz_0, \
                                         tg_xxxyz_xxyzzz_0, tg_xxxyz_xxzzzz_0, tg_xxxyz_xyyyyy_0, tg_xxxyz_xyyyyz_0, \
                                         tg_xxxyz_xyyyzz_0, tg_xxxyz_xyyzzz_0, tg_xxxyz_xyzzzz_0, tg_xxxyz_xzzzzz_0, \
                                         tg_xxxyz_yyyyyy_0, tg_xxxyz_yyyyyz_0, tg_xxxyz_yyyyzz_0, tg_xxxyz_yyyzzz_0, \
                                         tg_xxxyz_yyzzzz_0, tg_xxxyz_yzzzzz_0, tg_xxxyz_zzzzzz_0, tg_xxxzz_xxxxxx_0, \
                                         tg_xxxzz_xxxxxy_0, tg_xxxzz_xxxxxz_0, tg_xxxzz_xxxxyy_0, tg_xxxzz_xxxxyz_0, \
                                         tg_xxxzz_xxxxzz_0, tg_xxxzz_xxxyyy_0, tg_xxxzz_xxxyyz_0, tg_xxxzz_xxxyzz_0, \
                                         tg_xxxzz_xxxzzz_0, tg_xxxzz_xxyyyy_0, tg_xxxzz_xxyyyz_0, tg_xxxzz_xxyyzz_0, \
                                         tg_xxxzz_xxyzzz_0, tg_xxxzz_xxzzzz_0, tg_xxxzz_xyyyyy_0, tg_xxxzz_xyyyyz_0, \
                                         tg_xxxzz_xyyyzz_0, tg_xxxzz_xyyzzz_0, tg_xxxzz_xyzzzz_0, tg_xxxzz_xzzzzz_0, \
                                         tg_xxxzz_yyyyyy_0, tg_xxxzz_yyyyyz_0, tg_xxxzz_yyyyzz_0, tg_xxxzz_yyyzzz_0, \
                                         tg_xxxzz_yyzzzz_0, tg_xxxzz_yzzzzz_0, tg_xxxzz_zzzzzz_0, tg_xxyy_xxzzzz_0, \
                                         tg_xxyy_xxzzzz_1, tg_xxyy_xyyyyy_0, tg_xxyy_xyyyyy_1, tg_xxyy_xyyyyz_0, \
                                         tg_xxyy_xyyyyz_1, tg_xxyy_xyyyzz_0, tg_xxyy_xyyyzz_1, tg_xxyy_xyyzzz_0, \
                                         tg_xxyy_xyyzzz_1, tg_xxyy_xyzzzz_0, tg_xxyy_xyzzzz_1, tg_xxyy_xzzzz_1, \
                                         tg_xxyy_xzzzzz_0, tg_xxyy_xzzzzz_1, tg_xxyy_yyyyy_1, tg_xxyy_yyyyyy_0, \
                                         tg_xxyy_yyyyyy_1, tg_xxyy_yyyyyz_0, tg_xxyy_yyyyyz_1, tg_xxyy_yyyyz_1, \
                                         tg_xxyy_yyyyzz_0, tg_xxyy_yyyyzz_1, tg_xxyy_yyyzz_1, tg_xxyy_yyyzzz_0, \
                                         tg_xxyy_yyyzzz_1, tg_xxyy_yyzzz_1, tg_xxyy_yyzzzz_0, tg_xxyy_yyzzzz_1, \
                                         tg_xxyy_yzzzz_1, tg_xxyy_yzzzzz_0, tg_xxyy_yzzzzz_1, tg_xxyy_zzzzz_1, \
                                         tg_xxyy_zzzzzz_0, tg_xxyy_zzzzzz_1, tg_xxyyy_xxxxxx_0, tg_xxyyy_xxxxxy_0, \
                                         tg_xxyyy_xxxxxz_0, tg_xxyyy_xxxxyy_0, tg_xxyyy_xxxxyz_0, tg_xxyyy_xxxxzz_0, \
                                         tg_xxyyy_xxxyyy_0, tg_xxyyy_xxxyyz_0, tg_xxyyy_xxxyzz_0, tg_xxyyy_xxxzzz_0, \
                                         tg_xxyyy_xxyyyy_0, tg_xxyyy_xxyyyz_0, tg_xxyyy_xxyyzz_0, tg_xxyyy_xxyzzz_0, \
                                         tg_xxyyy_xxzzzz_0, tg_xxyyy_xyyyyy_0, tg_xxyyy_xyyyyz_0, tg_xxyyy_xyyyzz_0, \
                                         tg_xxyyy_xyyzzz_0, tg_xxyyy_xyzzzz_0, tg_xxyyy_xzzzzz_0, tg_xxyyy_yyyyyy_0, \
                                         tg_xxyyy_yyyyyz_0, tg_xxyyy_yyyyzz_0, tg_xxyyy_yyyzzz_0, tg_xxyyy_yyzzzz_0, \
                                         tg_xxyyy_yzzzzz_0, tg_xxyyy_zzzzzz_0, tg_xxyz_xxxxx_1, tg_xxyz_xxxxxx_0, \
                                         tg_xxyz_xxxxxx_1, tg_xxyz_xxxxxy_0, tg_xxyz_xxxxxy_1, tg_xxyz_xxxxxz_0, \
                                         tg_xxyz_xxxxxz_1, tg_xxyz_xxxxy_1, tg_xxyz_xxxxyy_0, tg_xxyz_xxxxyy_1, \
                                         tg_xxyz_xxxxyz_0, tg_xxyz_xxxxyz_1, tg_xxyz_xxxxz_1, tg_xxyz_xxxxzz_0, \
                                         tg_xxyz_xxxxzz_1, tg_xxyz_xxxyy_1, tg_xxyz_xxxyyy_0, tg_xxyz_xxxyyy_1, \
                                         tg_xxyz_xxxyyz_0, tg_xxyz_xxxyyz_1, tg_xxyz_xxxyz_1, tg_xxyz_xxxyzz_0, \
                                         tg_xxyz_xxxyzz_1, tg_xxyz_xxxzz_1, tg_xxyz_xxxzzz_0, tg_xxyz_xxxzzz_1, \
                                         tg_xxyz_xxyyy_1, tg_xxyz_xxyyyy_0, tg_xxyz_xxyyyy_1, tg_xxyz_xxyyyz_0, \
                                         tg_xxyz_xxyyyz_1, tg_xxyz_xxyyz_1, tg_xxyz_xxyyzz_0, tg_xxyz_xxyyzz_1, \
                                         tg_xxyz_xxyzz_1, tg_xxyz_xxyzzz_0, tg_xxyz_xxyzzz_1, tg_xxyz_xxzzz_1, \
                                         tg_xxyz_xxzzzz_0, tg_xxyz_xxzzzz_1, tg_xxyz_xyyyy_1, tg_xxyz_xyyyyy_0, \
                                         tg_xxyz_xyyyyy_1, tg_xxyz_xyyyyz_0, tg_xxyz_xyyyyz_1, tg_xxyz_xyyyz_1, \
                                         tg_xxyz_xyyyzz_0, tg_xxyz_xyyyzz_1, tg_xxyz_xyyzz_1, tg_xxyz_xyyzzz_0, \
                                         tg_xxyz_xyyzzz_1, tg_xxyz_xyzzz_1, tg_xxyz_xyzzzz_0, tg_xxyz_xyzzzz_1, \
                                         tg_xxyz_xzzzz_1, tg_xxyz_xzzzzz_0, tg_xxyz_xzzzzz_1, tg_xxyz_yyyyy_1, \
                                         tg_xxyz_yyyyyy_0, tg_xxyz_yyyyyy_1, tg_xxyz_yyyyyz_0, tg_xxyz_yyyyyz_1, \
                                         tg_xxyz_yyyyz_1, tg_xxyz_yyyyzz_0, tg_xxyz_yyyyzz_1, tg_xxyz_yyyzz_1, \
                                         tg_xxyz_yyyzzz_0, tg_xxyz_yyyzzz_1, tg_xxyz_yyzzz_1, tg_xxyz_yyzzzz_0, \
                                         tg_xxyz_yyzzzz_1, tg_xxyz_yzzzz_1, tg_xxyz_yzzzzz_0, tg_xxyz_yzzzzz_1, \
                                         tg_xxyz_zzzzz_1, tg_xxyz_zzzzzz_0, tg_xxyz_zzzzzz_1, tg_xxzz_xxxxx_1, \
                                         tg_xxzz_xxxxxx_0, tg_xxzz_xxxxxx_1, tg_xxzz_xxxxxy_0, tg_xxzz_xxxxxy_1, \
                                         tg_xxzz_xxxxxz_0, tg_xxzz_xxxxxz_1, tg_xxzz_xxxxy_1, tg_xxzz_xxxxyy_0, \
                                         tg_xxzz_xxxxyy_1, tg_xxzz_xxxxyz_0, tg_xxzz_xxxxyz_1, tg_xxzz_xxxxz_1, \
                                         tg_xxzz_xxxxzz_0, tg_xxzz_xxxxzz_1, tg_xxzz_xxxyy_1, tg_xxzz_xxxyyy_0, \
                                         tg_xxzz_xxxyyy_1, tg_xxzz_xxxyyz_0, tg_xxzz_xxxyyz_1, tg_xxzz_xxxyz_1, \
                                         tg_xxzz_xxxyzz_0, tg_xxzz_xxxyzz_1, tg_xxzz_xxxzz_1, tg_xxzz_xxxzzz_0, \
                                         tg_xxzz_xxxzzz_1, tg_xxzz_xxyyy_1, tg_xxzz_xxyyyy_0, tg_xxzz_xxyyyy_1, \
                                         tg_xxzz_xxyyyz_0, tg_xxzz_xxyyyz_1, tg_xxzz_xxyyz_1, tg_xxzz_xxyyzz_0, \
                                         tg_xxzz_xxyyzz_1, tg_xxzz_xxyzz_1, tg_xxzz_xxyzzz_0, tg_xxzz_xxyzzz_1, \
                                         tg_xxzz_xxzzz_1, tg_xxzz_xxzzzz_0, tg_xxzz_xxzzzz_1, tg_xxzz_xyyyy_1, \
                                         tg_xxzz_xyyyyy_0, tg_xxzz_xyyyyy_1, tg_xxzz_xyyyyz_0, tg_xxzz_xyyyyz_1, \
                                         tg_xxzz_xyyyz_1, tg_xxzz_xyyyzz_0, tg_xxzz_xyyyzz_1, tg_xxzz_xyyzz_1, \
                                         tg_xxzz_xyyzzz_0, tg_xxzz_xyyzzz_1, tg_xxzz_xyzzz_1, tg_xxzz_xyzzzz_0, \
                                         tg_xxzz_xyzzzz_1, tg_xxzz_xzzzz_1, tg_xxzz_xzzzzz_0, tg_xxzz_xzzzzz_1, \
                                         tg_xxzz_yyyyy_1, tg_xxzz_yyyyyy_0, tg_xxzz_yyyyyy_1, tg_xxzz_yyyyyz_0, \
                                         tg_xxzz_yyyyyz_1, tg_xxzz_yyyyz_1, tg_xxzz_yyyyzz_0, tg_xxzz_yyyyzz_1, \
                                         tg_xxzz_yyyzz_1, tg_xxzz_yyyzzz_0, tg_xxzz_yyyzzz_1, tg_xxzz_yyzzz_1, \
                                         tg_xxzz_yyzzzz_0, tg_xxzz_yyzzzz_1, tg_xxzz_yzzzz_1, tg_xxzz_yzzzzz_0, \
                                         tg_xxzz_yzzzzz_1, tg_xxzz_zzzzz_1, tg_xxzz_zzzzzz_0, tg_xxzz_zzzzzz_1, \
                                         tg_xyy_xxzzzz_0, tg_xyy_xxzzzz_1, tg_xyy_xyyyyy_0, tg_xyy_xyyyyy_1, tg_xyy_xyyyyz_0, \
                                         tg_xyy_xyyyyz_1, tg_xyy_xyyyzz_0, tg_xyy_xyyyzz_1, tg_xyy_xyyzzz_0, tg_xyy_xyyzzz_1, \
                                         tg_xyy_xyzzzz_0, tg_xyy_xyzzzz_1, tg_xyy_xzzzzz_0, tg_xyy_xzzzzz_1, tg_xyy_yyyyyy_0, \
                                         tg_xyy_yyyyyy_1, tg_xyy_yyyyyz_0, tg_xyy_yyyyyz_1, tg_xyy_yyyyzz_0, tg_xyy_yyyyzz_1, \
                                         tg_xyy_yyyzzz_0, tg_xyy_yyyzzz_1, tg_xyy_yyzzzz_0, tg_xyy_yyzzzz_1, tg_xyy_yzzzzz_0, \
                                         tg_xyy_yzzzzz_1, tg_xyy_zzzzzz_0, tg_xyy_zzzzzz_1, tg_xyyy_xxxxx_1, \
                                         tg_xyyy_xxxxxx_0, tg_xyyy_xxxxxx_1, tg_xyyy_xxxxxy_0, tg_xyyy_xxxxxy_1, \
                                         tg_xyyy_xxxxxz_0, tg_xyyy_xxxxxz_1, tg_xyyy_xxxxy_1, tg_xyyy_xxxxyy_0, \
                                         tg_xyyy_xxxxyy_1, tg_xyyy_xxxxyz_0, tg_xyyy_xxxxyz_1, tg_xyyy_xxxxz_1, \
                                         tg_xyyy_xxxxzz_0, tg_xyyy_xxxxzz_1, tg_xyyy_xxxyy_1, tg_xyyy_xxxyyy_0, \
                                         tg_xyyy_xxxyyy_1, tg_xyyy_xxxyyz_0, tg_xyyy_xxxyyz_1, tg_xyyy_xxxyz_1, \
                                         tg_xyyy_xxxyzz_0, tg_xyyy_xxxyzz_1, tg_xyyy_xxxzz_1, tg_xyyy_xxxzzz_0, \
                                         tg_xyyy_xxxzzz_1, tg_xyyy_xxyyy_1, tg_xyyy_xxyyyy_0, tg_xyyy_xxyyyy_1, \
                                         tg_xyyy_xxyyyz_0, tg_xyyy_xxyyyz_1, tg_xyyy_xxyyz_1, tg_xyyy_xxyyzz_0, \
                                         tg_xyyy_xxyyzz_1, tg_xyyy_xxyzz_1, tg_xyyy_xxyzzz_0, tg_xyyy_xxyzzz_1, \
                                         tg_xyyy_xxzzz_1, tg_xyyy_xxzzzz_0, tg_xyyy_xxzzzz_1, tg_xyyy_xyyyy_1, \
                                         tg_xyyy_xyyyyy_0, tg_xyyy_xyyyyy_1, tg_xyyy_xyyyyz_0, tg_xyyy_xyyyyz_1, \
                                         tg_xyyy_xyyyz_1, tg_xyyy_xyyyzz_0, tg_xyyy_xyyyzz_1, tg_xyyy_xyyzz_1, \
                                         tg_xyyy_xyyzzz_0, tg_xyyy_xyyzzz_1, tg_xyyy_xyzzz_1, tg_xyyy_xyzzzz_0, \
                                         tg_xyyy_xyzzzz_1, tg_xyyy_xzzzz_1, tg_xyyy_xzzzzz_0, tg_xyyy_xzzzzz_1, \
                                         tg_xyyy_yyyyy_1, tg_xyyy_yyyyyy_0, tg_xyyy_yyyyyy_1, tg_xyyy_yyyyyz_0, \
                                         tg_xyyy_yyyyyz_1, tg_xyyy_yyyyz_1, tg_xyyy_yyyyzz_0, tg_xyyy_yyyyzz_1, \
                                         tg_xyyy_yyyzz_1, tg_xyyy_yyyzzz_0, tg_xyyy_yyyzzz_1, tg_xyyy_yyzzz_1, \
                                         tg_xyyy_yyzzzz_0, tg_xyyy_yyzzzz_1, tg_xyyy_yzzzz_1, tg_xyyy_yzzzzz_0, \
                                         tg_xyyy_yzzzzz_1, tg_xyyy_zzzzz_1, tg_xyyy_zzzzzz_0, tg_xyyy_zzzzzz_1, \
                                         tg_xyz_xxxxxx_0, tg_xyz_xxxxxx_1, tg_xyz_xxxxxy_0, tg_xyz_xxxxxy_1, tg_xyz_xxxxxz_0, \
                                         tg_xyz_xxxxxz_1, tg_xyz_xxxxyy_0, tg_xyz_xxxxyy_1, tg_xyz_xxxxyz_0, tg_xyz_xxxxyz_1, \
                                         tg_xyz_xxxxzz_0, tg_xyz_xxxxzz_1, tg_xyz_xxxyyy_0, tg_xyz_xxxyyy_1, tg_xyz_xxxyyz_0, \
                                         tg_xyz_xxxyyz_1, tg_xyz_xxxyzz_0, tg_xyz_xxxyzz_1, tg_xyz_xxxzzz_0, tg_xyz_xxxzzz_1, \
                                         tg_xyz_xxyyyy_0, tg_xyz_xxyyyy_1, tg_xyz_xxyyyz_0, tg_xyz_xxyyyz_1, tg_xyz_xxyyzz_0, \
                                         tg_xyz_xxyyzz_1, tg_xyz_xxyzzz_0, tg_xyz_xxyzzz_1, tg_xyz_xxzzzz_0, tg_xyz_xxzzzz_1, \
                                         tg_xyz_xyyyyy_0, tg_xyz_xyyyyy_1, tg_xyz_xyyyyz_0, tg_xyz_xyyyyz_1, tg_xyz_xyyyzz_0, \
                                         tg_xyz_xyyyzz_1, tg_xyz_xyyzzz_0, tg_xyz_xyyzzz_1, tg_xyz_xyzzzz_0, tg_xyz_xyzzzz_1, \
                                         tg_xyz_xzzzzz_0, tg_xyz_xzzzzz_1, tg_xyz_yyyyyy_0, tg_xyz_yyyyyy_1, tg_xyz_yyyyyz_0, \
                                         tg_xyz_yyyyyz_1, tg_xyz_yyyyzz_0, tg_xyz_yyyyzz_1, tg_xyz_yyyzzz_0, tg_xyz_yyyzzz_1, \
                                         tg_xyz_yyzzzz_0, tg_xyz_yyzzzz_1, tg_xyz_yzzzzz_0, tg_xyz_yzzzzz_1, tg_xyz_zzzzzz_0, \
                                         tg_xyz_zzzzzz_1, tg_xzz_xxxxxx_0, tg_xzz_xxxxxx_1, tg_xzz_xxxxxy_0, tg_xzz_xxxxxy_1, \
                                         tg_xzz_xxxxxz_0, tg_xzz_xxxxxz_1, tg_xzz_xxxxyy_0, tg_xzz_xxxxyy_1, tg_xzz_xxxxyz_0, \
                                         tg_xzz_xxxxyz_1, tg_xzz_xxxxzz_0, tg_xzz_xxxxzz_1, tg_xzz_xxxyyy_0, tg_xzz_xxxyyy_1, \
                                         tg_xzz_xxxyyz_0, tg_xzz_xxxyyz_1, tg_xzz_xxxyzz_0, tg_xzz_xxxyzz_1, tg_xzz_xxxzzz_0, \
                                         tg_xzz_xxxzzz_1, tg_xzz_xxyyyy_0, tg_xzz_xxyyyy_1, tg_xzz_xxyyyz_0, tg_xzz_xxyyyz_1, \
                                         tg_xzz_xxyyzz_0, tg_xzz_xxyyzz_1, tg_xzz_xxyzzz_0, tg_xzz_xxyzzz_1, tg_xzz_xxzzzz_0, \
                                         tg_xzz_xxzzzz_1, tg_xzz_xyyyyy_0, tg_xzz_xyyyyy_1, tg_xzz_xyyyyz_0, tg_xzz_xyyyyz_1, \
                                         tg_xzz_xyyyzz_0, tg_xzz_xyyyzz_1, tg_xzz_xyyzzz_0, tg_xzz_xyyzzz_1, tg_xzz_xyzzzz_0, \
                                         tg_xzz_xyzzzz_1, tg_xzz_xzzzzz_0, tg_xzz_xzzzzz_1, tg_xzz_yyyyyy_0, tg_xzz_yyyyyy_1, \
                                         tg_xzz_yyyyyz_0, tg_xzz_yyyyyz_1, tg_xzz_yyyyzz_0, tg_xzz_yyyyzz_1, tg_xzz_yyyzzz_0, \
                                         tg_xzz_yyyzzz_1, tg_xzz_yyzzzz_0, tg_xzz_yyzzzz_1, tg_xzz_yzzzzz_0, tg_xzz_yzzzzz_1, \
                                         tg_xzz_zzzzzz_0, tg_xzz_zzzzzz_1, tg_yyy_xxxxxx_0, tg_yyy_xxxxxx_1, tg_yyy_xxxxxy_0, \
                                         tg_yyy_xxxxxy_1, tg_yyy_xxxxxz_0, tg_yyy_xxxxxz_1, tg_yyy_xxxxyy_0, tg_yyy_xxxxyy_1, \
                                         tg_yyy_xxxxyz_0, tg_yyy_xxxxyz_1, tg_yyy_xxxxzz_0, tg_yyy_xxxxzz_1, tg_yyy_xxxyyy_0, \
                                         tg_yyy_xxxyyy_1, tg_yyy_xxxyyz_0, tg_yyy_xxxyyz_1, tg_yyy_xxxyzz_0, tg_yyy_xxxyzz_1, \
                                         tg_yyy_xxxzzz_0, tg_yyy_xxxzzz_1, tg_yyy_xxyyyy_0, tg_yyy_xxyyyy_1, tg_yyy_xxyyyz_0, \
                                         tg_yyy_xxyyyz_1, tg_yyy_xxyyzz_0, tg_yyy_xxyyzz_1, tg_yyy_xxyzzz_0, tg_yyy_xxyzzz_1, \
                                         tg_yyy_xxzzzz_0, tg_yyy_xxzzzz_1, tg_yyy_xyyyyy_0, tg_yyy_xyyyyy_1, tg_yyy_xyyyyz_0, \
                                         tg_yyy_xyyyyz_1, tg_yyy_xyyyzz_0, tg_yyy_xyyyzz_1, tg_yyy_xyyzzz_0, tg_yyy_xyyzzz_1, \
                                         tg_yyy_xyzzzz_0, tg_yyy_xyzzzz_1, tg_yyy_xzzzzz_0, tg_yyy_xzzzzz_1, tg_yyy_yyyyyy_0, \
                                         tg_yyy_yyyyyy_1, tg_yyy_yyyyyz_0, tg_yyy_yyyyyz_1, tg_yyy_yyyyzz_0, tg_yyy_yyyyzz_1, \
                                         tg_yyy_yyyzzz_0, tg_yyy_yyyzzz_1, tg_yyy_yyzzzz_0, tg_yyy_yyzzzz_1, tg_yyy_yzzzzz_0, \
                                         tg_yyy_yzzzzz_1, tg_yyy_zzzzzz_0, tg_yyy_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyy_xxzzzz_0[j] = pb_x * tg_xxyy_xxzzzz_0[j] + fr * tg_xxyy_xxzzzz_1[j] + fl1_fx * (tg_xyy_xxzzzz_0[j] - tg_xyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyy_xzzzz_1[j];

                    tg_xxxyy_xyyyyy_0[j] = pb_x * tg_xxyy_xyyyyy_0[j] + fr * tg_xxyy_xyyyyy_1[j] + fl1_fx * (tg_xyy_xyyyyy_0[j] - tg_xyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyyy_1[j];

                    tg_xxxyy_xyyyyz_0[j] = pb_x * tg_xxyy_xyyyyz_0[j] + fr * tg_xxyy_xyyyyz_1[j] + fl1_fx * (tg_xyy_xyyyyz_0[j] - tg_xyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyyz_1[j];

                    tg_xxxyy_xyyyzz_0[j] = pb_x * tg_xxyy_xyyyzz_0[j] + fr * tg_xxyy_xyyyzz_1[j] + fl1_fx * (tg_xyy_xyyyzz_0[j] - tg_xyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyyzz_1[j];

                    tg_xxxyy_xyyzzz_0[j] = pb_x * tg_xxyy_xyyzzz_0[j] + fr * tg_xxyy_xyyzzz_1[j] + fl1_fx * (tg_xyy_xyyzzz_0[j] - tg_xyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yyzzz_1[j];

                    tg_xxxyy_xyzzzz_0[j] = pb_x * tg_xxyy_xyzzzz_0[j] + fr * tg_xxyy_xyzzzz_1[j] + fl1_fx * (tg_xyy_xyzzzz_0[j] - tg_xyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_yzzzz_1[j];

                    tg_xxxyy_xzzzzz_0[j] = pb_x * tg_xxyy_xzzzzz_0[j] + fr * tg_xxyy_xzzzzz_1[j] + fl1_fx * (tg_xyy_xzzzzz_0[j] - tg_xyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyy_zzzzz_1[j];

                    tg_xxxyy_yyyyyy_0[j] = pb_x * tg_xxyy_yyyyyy_0[j] + fr * tg_xxyy_yyyyyy_1[j] + fl1_fx * (tg_xyy_yyyyyy_0[j] - tg_xyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyy_yyyyyz_0[j] = pb_x * tg_xxyy_yyyyyz_0[j] + fr * tg_xxyy_yyyyyz_1[j] + fl1_fx * (tg_xyy_yyyyyz_0[j] - tg_xyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyy_yyyyzz_0[j] = pb_x * tg_xxyy_yyyyzz_0[j] + fr * tg_xxyy_yyyyzz_1[j] + fl1_fx * (tg_xyy_yyyyzz_0[j] - tg_xyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyy_yyyzzz_0[j] = pb_x * tg_xxyy_yyyzzz_0[j] + fr * tg_xxyy_yyyzzz_1[j] + fl1_fx * (tg_xyy_yyyzzz_0[j] - tg_xyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyy_yyzzzz_0[j] = pb_x * tg_xxyy_yyzzzz_0[j] + fr * tg_xxyy_yyzzzz_1[j] + fl1_fx * (tg_xyy_yyzzzz_0[j] - tg_xyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyy_yzzzzz_0[j] = pb_x * tg_xxyy_yzzzzz_0[j] + fr * tg_xxyy_yzzzzz_1[j] + fl1_fx * (tg_xyy_yzzzzz_0[j] - tg_xyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyy_zzzzzz_0[j] = pb_x * tg_xxyy_zzzzzz_0[j] + fr * tg_xxyy_zzzzzz_1[j] + fl1_fx * (tg_xyy_zzzzzz_0[j] - tg_xyy_zzzzzz_1[j] * fl1_fza);

                    tg_xxxyz_xxxxxx_0[j] = pb_x * tg_xxyz_xxxxxx_0[j] + fr * tg_xxyz_xxxxxx_1[j] + fl1_fx * (tg_xyz_xxxxxx_0[j] - tg_xyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxyz_xxxxx_1[j];

                    tg_xxxyz_xxxxxy_0[j] = pb_x * tg_xxyz_xxxxxy_0[j] + fr * tg_xxyz_xxxxxy_1[j] + fl1_fx * (tg_xyz_xxxxxy_0[j] - tg_xyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyz_xxxxy_1[j];

                    tg_xxxyz_xxxxxz_0[j] = pb_x * tg_xxyz_xxxxxz_0[j] + fr * tg_xxyz_xxxxxz_1[j] + fl1_fx * (tg_xyz_xxxxxz_0[j] - tg_xyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxyz_xxxxz_1[j];

                    tg_xxxyz_xxxxyy_0[j] = pb_x * tg_xxyz_xxxxyy_0[j] + fr * tg_xxyz_xxxxyy_1[j] + fl1_fx * (tg_xyz_xxxxyy_0[j] - tg_xyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxyy_1[j];

                    tg_xxxyz_xxxxyz_0[j] = pb_x * tg_xxyz_xxxxyz_0[j] + fr * tg_xxyz_xxxxyz_1[j] + fl1_fx * (tg_xyz_xxxxyz_0[j] - tg_xyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxyz_1[j];

                    tg_xxxyz_xxxxzz_0[j] = pb_x * tg_xxyz_xxxxzz_0[j] + fr * tg_xxyz_xxxxzz_1[j] + fl1_fx * (tg_xyz_xxxxzz_0[j] - tg_xyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyz_xxxzz_1[j];

                    tg_xxxyz_xxxyyy_0[j] = pb_x * tg_xxyz_xxxyyy_0[j] + fr * tg_xxyz_xxxyyy_1[j] + fl1_fx * (tg_xyz_xxxyyy_0[j] - tg_xyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyyy_1[j];

                    tg_xxxyz_xxxyyz_0[j] = pb_x * tg_xxyz_xxxyyz_0[j] + fr * tg_xxyz_xxxyyz_1[j] + fl1_fx * (tg_xyz_xxxyyz_0[j] - tg_xyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyyz_1[j];

                    tg_xxxyz_xxxyzz_0[j] = pb_x * tg_xxyz_xxxyzz_0[j] + fr * tg_xxyz_xxxyzz_1[j] + fl1_fx * (tg_xyz_xxxyzz_0[j] - tg_xyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxyzz_1[j];

                    tg_xxxyz_xxxzzz_0[j] = pb_x * tg_xxyz_xxxzzz_0[j] + fr * tg_xxyz_xxxzzz_1[j] + fl1_fx * (tg_xyz_xxxzzz_0[j] - tg_xyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyz_xxzzz_1[j];

                    tg_xxxyz_xxyyyy_0[j] = pb_x * tg_xxyz_xxyyyy_0[j] + fr * tg_xxyz_xxyyyy_1[j] + fl1_fx * (tg_xyz_xxyyyy_0[j] - tg_xyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyyy_1[j];

                    tg_xxxyz_xxyyyz_0[j] = pb_x * tg_xxyz_xxyyyz_0[j] + fr * tg_xxyz_xxyyyz_1[j] + fl1_fx * (tg_xyz_xxyyyz_0[j] - tg_xyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyyz_1[j];

                    tg_xxxyz_xxyyzz_0[j] = pb_x * tg_xxyz_xxyyzz_0[j] + fr * tg_xxyz_xxyyzz_1[j] + fl1_fx * (tg_xyz_xxyyzz_0[j] - tg_xyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyyzz_1[j];

                    tg_xxxyz_xxyzzz_0[j] = pb_x * tg_xxyz_xxyzzz_0[j] + fr * tg_xxyz_xxyzzz_1[j] + fl1_fx * (tg_xyz_xxyzzz_0[j] - tg_xyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xyzzz_1[j];

                    tg_xxxyz_xxzzzz_0[j] = pb_x * tg_xxyz_xxzzzz_0[j] + fr * tg_xxyz_xxzzzz_1[j] + fl1_fx * (tg_xyz_xxzzzz_0[j] - tg_xyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyz_xzzzz_1[j];

                    tg_xxxyz_xyyyyy_0[j] = pb_x * tg_xxyz_xyyyyy_0[j] + fr * tg_xxyz_xyyyyy_1[j] + fl1_fx * (tg_xyz_xyyyyy_0[j] - tg_xyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyyy_1[j];

                    tg_xxxyz_xyyyyz_0[j] = pb_x * tg_xxyz_xyyyyz_0[j] + fr * tg_xxyz_xyyyyz_1[j] + fl1_fx * (tg_xyz_xyyyyz_0[j] - tg_xyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyyz_1[j];

                    tg_xxxyz_xyyyzz_0[j] = pb_x * tg_xxyz_xyyyzz_0[j] + fr * tg_xxyz_xyyyzz_1[j] + fl1_fx * (tg_xyz_xyyyzz_0[j] - tg_xyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyyzz_1[j];

                    tg_xxxyz_xyyzzz_0[j] = pb_x * tg_xxyz_xyyzzz_0[j] + fr * tg_xxyz_xyyzzz_1[j] + fl1_fx * (tg_xyz_xyyzzz_0[j] - tg_xyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yyzzz_1[j];

                    tg_xxxyz_xyzzzz_0[j] = pb_x * tg_xxyz_xyzzzz_0[j] + fr * tg_xxyz_xyzzzz_1[j] + fl1_fx * (tg_xyz_xyzzzz_0[j] - tg_xyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_yzzzz_1[j];

                    tg_xxxyz_xzzzzz_0[j] = pb_x * tg_xxyz_xzzzzz_0[j] + fr * tg_xxyz_xzzzzz_1[j] + fl1_fx * (tg_xyz_xzzzzz_0[j] - tg_xyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyz_zzzzz_1[j];

                    tg_xxxyz_yyyyyy_0[j] = pb_x * tg_xxyz_yyyyyy_0[j] + fr * tg_xxyz_yyyyyy_1[j] + fl1_fx * (tg_xyz_yyyyyy_0[j] - tg_xyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxyz_yyyyyz_0[j] = pb_x * tg_xxyz_yyyyyz_0[j] + fr * tg_xxyz_yyyyyz_1[j] + fl1_fx * (tg_xyz_yyyyyz_0[j] - tg_xyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxyz_yyyyzz_0[j] = pb_x * tg_xxyz_yyyyzz_0[j] + fr * tg_xxyz_yyyyzz_1[j] + fl1_fx * (tg_xyz_yyyyzz_0[j] - tg_xyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxyz_yyyzzz_0[j] = pb_x * tg_xxyz_yyyzzz_0[j] + fr * tg_xxyz_yyyzzz_1[j] + fl1_fx * (tg_xyz_yyyzzz_0[j] - tg_xyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxyz_yyzzzz_0[j] = pb_x * tg_xxyz_yyzzzz_0[j] + fr * tg_xxyz_yyzzzz_1[j] + fl1_fx * (tg_xyz_yyzzzz_0[j] - tg_xyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxyz_yzzzzz_0[j] = pb_x * tg_xxyz_yzzzzz_0[j] + fr * tg_xxyz_yzzzzz_1[j] + fl1_fx * (tg_xyz_yzzzzz_0[j] - tg_xyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxyz_zzzzzz_0[j] = pb_x * tg_xxyz_zzzzzz_0[j] + fr * tg_xxyz_zzzzzz_1[j] + fl1_fx * (tg_xyz_zzzzzz_0[j] - tg_xyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxxzz_xxxxxx_0[j] = pb_x * tg_xxzz_xxxxxx_0[j] + fr * tg_xxzz_xxxxxx_1[j] + fl1_fx * (tg_xzz_xxxxxx_0[j] - tg_xzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xxzz_xxxxx_1[j];

                    tg_xxxzz_xxxxxy_0[j] = pb_x * tg_xxzz_xxxxxy_0[j] + fr * tg_xxzz_xxxxxy_1[j] + fl1_fx * (tg_xzz_xxxxxy_0[j] - tg_xzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzz_xxxxy_1[j];

                    tg_xxxzz_xxxxxz_0[j] = pb_x * tg_xxzz_xxxxxz_0[j] + fr * tg_xxzz_xxxxxz_1[j] + fl1_fx * (tg_xzz_xxxxxz_0[j] - tg_xzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xxzz_xxxxz_1[j];

                    tg_xxxzz_xxxxyy_0[j] = pb_x * tg_xxzz_xxxxyy_0[j] + fr * tg_xxzz_xxxxyy_1[j] + fl1_fx * (tg_xzz_xxxxyy_0[j] - tg_xzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxyy_1[j];

                    tg_xxxzz_xxxxyz_0[j] = pb_x * tg_xxzz_xxxxyz_0[j] + fr * tg_xxzz_xxxxyz_1[j] + fl1_fx * (tg_xzz_xxxxyz_0[j] - tg_xzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxyz_1[j];

                    tg_xxxzz_xxxxzz_0[j] = pb_x * tg_xxzz_xxxxzz_0[j] + fr * tg_xxzz_xxxxzz_1[j] + fl1_fx * (tg_xzz_xxxxzz_0[j] - tg_xzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzz_xxxzz_1[j];

                    tg_xxxzz_xxxyyy_0[j] = pb_x * tg_xxzz_xxxyyy_0[j] + fr * tg_xxzz_xxxyyy_1[j] + fl1_fx * (tg_xzz_xxxyyy_0[j] - tg_xzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyyy_1[j];

                    tg_xxxzz_xxxyyz_0[j] = pb_x * tg_xxzz_xxxyyz_0[j] + fr * tg_xxzz_xxxyyz_1[j] + fl1_fx * (tg_xzz_xxxyyz_0[j] - tg_xzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyyz_1[j];

                    tg_xxxzz_xxxyzz_0[j] = pb_x * tg_xxzz_xxxyzz_0[j] + fr * tg_xxzz_xxxyzz_1[j] + fl1_fx * (tg_xzz_xxxyzz_0[j] - tg_xzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxyzz_1[j];

                    tg_xxxzz_xxxzzz_0[j] = pb_x * tg_xxzz_xxxzzz_0[j] + fr * tg_xxzz_xxxzzz_1[j] + fl1_fx * (tg_xzz_xxxzzz_0[j] - tg_xzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzz_xxzzz_1[j];

                    tg_xxxzz_xxyyyy_0[j] = pb_x * tg_xxzz_xxyyyy_0[j] + fr * tg_xxzz_xxyyyy_1[j] + fl1_fx * (tg_xzz_xxyyyy_0[j] - tg_xzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyyy_1[j];

                    tg_xxxzz_xxyyyz_0[j] = pb_x * tg_xxzz_xxyyyz_0[j] + fr * tg_xxzz_xxyyyz_1[j] + fl1_fx * (tg_xzz_xxyyyz_0[j] - tg_xzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyyz_1[j];

                    tg_xxxzz_xxyyzz_0[j] = pb_x * tg_xxzz_xxyyzz_0[j] + fr * tg_xxzz_xxyyzz_1[j] + fl1_fx * (tg_xzz_xxyyzz_0[j] - tg_xzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyyzz_1[j];

                    tg_xxxzz_xxyzzz_0[j] = pb_x * tg_xxzz_xxyzzz_0[j] + fr * tg_xxzz_xxyzzz_1[j] + fl1_fx * (tg_xzz_xxyzzz_0[j] - tg_xzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xyzzz_1[j];

                    tg_xxxzz_xxzzzz_0[j] = pb_x * tg_xxzz_xxzzzz_0[j] + fr * tg_xxzz_xxzzzz_1[j] + fl1_fx * (tg_xzz_xxzzzz_0[j] - tg_xzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzz_xzzzz_1[j];

                    tg_xxxzz_xyyyyy_0[j] = pb_x * tg_xxzz_xyyyyy_0[j] + fr * tg_xxzz_xyyyyy_1[j] + fl1_fx * (tg_xzz_xyyyyy_0[j] - tg_xzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyyy_1[j];

                    tg_xxxzz_xyyyyz_0[j] = pb_x * tg_xxzz_xyyyyz_0[j] + fr * tg_xxzz_xyyyyz_1[j] + fl1_fx * (tg_xzz_xyyyyz_0[j] - tg_xzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyyz_1[j];

                    tg_xxxzz_xyyyzz_0[j] = pb_x * tg_xxzz_xyyyzz_0[j] + fr * tg_xxzz_xyyyzz_1[j] + fl1_fx * (tg_xzz_xyyyzz_0[j] - tg_xzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyyzz_1[j];

                    tg_xxxzz_xyyzzz_0[j] = pb_x * tg_xxzz_xyyzzz_0[j] + fr * tg_xxzz_xyyzzz_1[j] + fl1_fx * (tg_xzz_xyyzzz_0[j] - tg_xzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yyzzz_1[j];

                    tg_xxxzz_xyzzzz_0[j] = pb_x * tg_xxzz_xyzzzz_0[j] + fr * tg_xxzz_xyzzzz_1[j] + fl1_fx * (tg_xzz_xyzzzz_0[j] - tg_xzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_yzzzz_1[j];

                    tg_xxxzz_xzzzzz_0[j] = pb_x * tg_xxzz_xzzzzz_0[j] + fr * tg_xxzz_xzzzzz_1[j] + fl1_fx * (tg_xzz_xzzzzz_0[j] - tg_xzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzz_zzzzz_1[j];

                    tg_xxxzz_yyyyyy_0[j] = pb_x * tg_xxzz_yyyyyy_0[j] + fr * tg_xxzz_yyyyyy_1[j] + fl1_fx * (tg_xzz_yyyyyy_0[j] - tg_xzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxxzz_yyyyyz_0[j] = pb_x * tg_xxzz_yyyyyz_0[j] + fr * tg_xxzz_yyyyyz_1[j] + fl1_fx * (tg_xzz_yyyyyz_0[j] - tg_xzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxxzz_yyyyzz_0[j] = pb_x * tg_xxzz_yyyyzz_0[j] + fr * tg_xxzz_yyyyzz_1[j] + fl1_fx * (tg_xzz_yyyyzz_0[j] - tg_xzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxxzz_yyyzzz_0[j] = pb_x * tg_xxzz_yyyzzz_0[j] + fr * tg_xxzz_yyyzzz_1[j] + fl1_fx * (tg_xzz_yyyzzz_0[j] - tg_xzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxxzz_yyzzzz_0[j] = pb_x * tg_xxzz_yyzzzz_0[j] + fr * tg_xxzz_yyzzzz_1[j] + fl1_fx * (tg_xzz_yyzzzz_0[j] - tg_xzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxxzz_yzzzzz_0[j] = pb_x * tg_xxzz_yzzzzz_0[j] + fr * tg_xxzz_yzzzzz_1[j] + fl1_fx * (tg_xzz_yzzzzz_0[j] - tg_xzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxxzz_zzzzzz_0[j] = pb_x * tg_xxzz_zzzzzz_0[j] + fr * tg_xxzz_zzzzzz_1[j] + fl1_fx * (tg_xzz_zzzzzz_0[j] - tg_xzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyyy_xxxxxx_0[j] = pb_x * tg_xyyy_xxxxxx_0[j] + fr * tg_xyyy_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxx_0[j] - tg_yyy_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyy_xxxxx_1[j];

                    tg_xxyyy_xxxxxy_0[j] = pb_x * tg_xyyy_xxxxxy_0[j] + fr * tg_xyyy_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxy_0[j] - tg_yyy_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyy_xxxxy_1[j];

                    tg_xxyyy_xxxxxz_0[j] = pb_x * tg_xyyy_xxxxxz_0[j] + fr * tg_xyyy_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxxz_0[j] - tg_yyy_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyy_xxxxz_1[j];

                    tg_xxyyy_xxxxyy_0[j] = pb_x * tg_xyyy_xxxxyy_0[j] + fr * tg_xyyy_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxyy_0[j] - tg_yyy_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxyy_1[j];

                    tg_xxyyy_xxxxyz_0[j] = pb_x * tg_xyyy_xxxxyz_0[j] + fr * tg_xyyy_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxyz_0[j] - tg_yyy_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxyz_1[j];

                    tg_xxyyy_xxxxzz_0[j] = pb_x * tg_xyyy_xxxxzz_0[j] + fr * tg_xyyy_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxxzz_0[j] - tg_yyy_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyy_xxxzz_1[j];

                    tg_xxyyy_xxxyyy_0[j] = pb_x * tg_xyyy_xxxyyy_0[j] + fr * tg_xyyy_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyyy_0[j] - tg_yyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyyy_1[j];

                    tg_xxyyy_xxxyyz_0[j] = pb_x * tg_xyyy_xxxyyz_0[j] + fr * tg_xyyy_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyyz_0[j] - tg_yyy_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyyz_1[j];

                    tg_xxyyy_xxxyzz_0[j] = pb_x * tg_xyyy_xxxyzz_0[j] + fr * tg_xyyy_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxyzz_0[j] - tg_yyy_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxyzz_1[j];

                    tg_xxyyy_xxxzzz_0[j] = pb_x * tg_xyyy_xxxzzz_0[j] + fr * tg_xyyy_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxxzzz_0[j] - tg_yyy_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyy_xxzzz_1[j];

                    tg_xxyyy_xxyyyy_0[j] = pb_x * tg_xyyy_xxyyyy_0[j] + fr * tg_xyyy_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyyy_0[j] - tg_yyy_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyyy_1[j];

                    tg_xxyyy_xxyyyz_0[j] = pb_x * tg_xyyy_xxyyyz_0[j] + fr * tg_xyyy_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyyz_0[j] - tg_yyy_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyyz_1[j];

                    tg_xxyyy_xxyyzz_0[j] = pb_x * tg_xyyy_xxyyzz_0[j] + fr * tg_xyyy_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyyzz_0[j] - tg_yyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyyzz_1[j];

                    tg_xxyyy_xxyzzz_0[j] = pb_x * tg_xyyy_xxyzzz_0[j] + fr * tg_xyyy_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxyzzz_0[j] - tg_yyy_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xyzzz_1[j];

                    tg_xxyyy_xxzzzz_0[j] = pb_x * tg_xyyy_xxzzzz_0[j] + fr * tg_xyyy_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xxzzzz_0[j] - tg_yyy_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyy_xzzzz_1[j];

                    tg_xxyyy_xyyyyy_0[j] = pb_x * tg_xyyy_xyyyyy_0[j] + fr * tg_xyyy_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyyy_0[j] - tg_yyy_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyyy_1[j];

                    tg_xxyyy_xyyyyz_0[j] = pb_x * tg_xyyy_xyyyyz_0[j] + fr * tg_xyyy_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyyz_0[j] - tg_yyy_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyyz_1[j];

                    tg_xxyyy_xyyyzz_0[j] = pb_x * tg_xyyy_xyyyzz_0[j] + fr * tg_xyyy_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyyzz_0[j] - tg_yyy_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyyzz_1[j];

                    tg_xxyyy_xyyzzz_0[j] = pb_x * tg_xyyy_xyyzzz_0[j] + fr * tg_xyyy_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyyzzz_0[j] - tg_yyy_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yyzzz_1[j];

                    tg_xxyyy_xyzzzz_0[j] = pb_x * tg_xyyy_xyzzzz_0[j] + fr * tg_xyyy_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xyzzzz_0[j] - tg_yyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_yzzzz_1[j];

                    tg_xxyyy_xzzzzz_0[j] = pb_x * tg_xyyy_xzzzzz_0[j] + fr * tg_xyyy_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_xzzzzz_0[j] - tg_yyy_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyy_zzzzz_1[j];

                    tg_xxyyy_yyyyyy_0[j] = pb_x * tg_xyyy_yyyyyy_0[j] + fr * tg_xyyy_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyyy_0[j] - tg_yyy_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyy_yyyyyz_0[j] = pb_x * tg_xyyy_yyyyyz_0[j] + fr * tg_xyyy_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyyz_0[j] - tg_yyy_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyy_yyyyzz_0[j] = pb_x * tg_xyyy_yyyyzz_0[j] + fr * tg_xyyy_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyyzz_0[j] - tg_yyy_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyy_yyyzzz_0[j] = pb_x * tg_xyyy_yyyzzz_0[j] + fr * tg_xyyy_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyyzzz_0[j] - tg_yyy_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyy_yyzzzz_0[j] = pb_x * tg_xyyy_yyzzzz_0[j] + fr * tg_xyyy_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yyzzzz_0[j] - tg_yyy_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyy_yzzzzz_0[j] = pb_x * tg_xyyy_yzzzzz_0[j] + fr * tg_xyyy_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_yzzzzz_0[j] - tg_yyy_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyy_zzzzzz_0[j] = pb_x * tg_xyyy_zzzzzz_0[j] + fr * tg_xyyy_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyy_zzzzzz_0[j] - tg_yyy_zzzzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSI_196_294(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                // set up pointers to integrals

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

                // Batch of Integrals (196,294)

                #pragma omp simd aligned(fxn, fza, tg_xxyyz_xxxxxx_0, tg_xxyyz_xxxxxy_0, tg_xxyyz_xxxxxz_0, \
                                         tg_xxyyz_xxxxyy_0, tg_xxyyz_xxxxyz_0, tg_xxyyz_xxxxzz_0, tg_xxyyz_xxxyyy_0, \
                                         tg_xxyyz_xxxyyz_0, tg_xxyyz_xxxyzz_0, tg_xxyyz_xxxzzz_0, tg_xxyyz_xxyyyy_0, \
                                         tg_xxyyz_xxyyyz_0, tg_xxyyz_xxyyzz_0, tg_xxyyz_xxyzzz_0, tg_xxyyz_xxzzzz_0, \
                                         tg_xxyyz_xyyyyy_0, tg_xxyyz_xyyyyz_0, tg_xxyyz_xyyyzz_0, tg_xxyyz_xyyzzz_0, \
                                         tg_xxyyz_xyzzzz_0, tg_xxyyz_xzzzzz_0, tg_xxyyz_yyyyyy_0, tg_xxyyz_yyyyyz_0, \
                                         tg_xxyyz_yyyyzz_0, tg_xxyyz_yyyzzz_0, tg_xxyyz_yyzzzz_0, tg_xxyyz_yzzzzz_0, \
                                         tg_xxyyz_zzzzzz_0, tg_xxyzz_xxxxxx_0, tg_xxyzz_xxxxxy_0, tg_xxyzz_xxxxxz_0, \
                                         tg_xxyzz_xxxxyy_0, tg_xxyzz_xxxxyz_0, tg_xxyzz_xxxxzz_0, tg_xxyzz_xxxyyy_0, \
                                         tg_xxyzz_xxxyyz_0, tg_xxyzz_xxxyzz_0, tg_xxyzz_xxxzzz_0, tg_xxyzz_xxyyyy_0, \
                                         tg_xxyzz_xxyyyz_0, tg_xxyzz_xxyyzz_0, tg_xxyzz_xxyzzz_0, tg_xxyzz_xxzzzz_0, \
                                         tg_xxyzz_xyyyyy_0, tg_xxyzz_xyyyyz_0, tg_xxyzz_xyyyzz_0, tg_xxyzz_xyyzzz_0, \
                                         tg_xxyzz_xyzzzz_0, tg_xxyzz_xzzzzz_0, tg_xxyzz_yyyyyy_0, tg_xxyzz_yyyyyz_0, \
                                         tg_xxyzz_yyyyzz_0, tg_xxyzz_yyyzzz_0, tg_xxyzz_yyzzzz_0, tg_xxyzz_yzzzzz_0, \
                                         tg_xxyzz_zzzzzz_0, tg_xxzzz_xxxxxx_0, tg_xxzzz_xxxxxy_0, tg_xxzzz_xxxxxz_0, \
                                         tg_xxzzz_xxxxyy_0, tg_xxzzz_xxxxyz_0, tg_xxzzz_xxxxzz_0, tg_xxzzz_xxxyyy_0, \
                                         tg_xxzzz_xxxyyz_0, tg_xxzzz_xxxyzz_0, tg_xxzzz_xxxzzz_0, tg_xxzzz_xxyyyy_0, \
                                         tg_xxzzz_xxyyyz_0, tg_xxzzz_xxyyzz_0, tg_xxzzz_xxyzzz_0, tg_xxzzz_xxzzzz_0, \
                                         tg_xxzzz_xyyyyy_0, tg_xxzzz_xyyyyz_0, tg_xxzzz_xyyyzz_0, tg_xxzzz_xyyzzz_0, \
                                         tg_xxzzz_xyzzzz_0, tg_xxzzz_xzzzzz_0, tg_xxzzz_yyyyyy_0, tg_xxzzz_yyyyyz_0, \
                                         tg_xxzzz_yyyyzz_0, tg_xxzzz_yyyzzz_0, tg_xxzzz_yyzzzz_0, tg_xxzzz_yzzzzz_0, \
                                         tg_xxzzz_zzzzzz_0, tg_xyyyy_xxxxxx_0, tg_xyyyy_xxxxxy_0, tg_xyyyy_xxxxxz_0, \
                                         tg_xyyyy_xxxxyy_0, tg_xyyyy_xxxxyz_0, tg_xyyyy_xxxxzz_0, tg_xyyyy_xxxyyy_0, \
                                         tg_xyyyy_xxxyyz_0, tg_xyyyy_xxxyzz_0, tg_xyyyy_xxxzzz_0, tg_xyyyy_xxyyyy_0, \
                                         tg_xyyyy_xxyyyz_0, tg_xyyyy_xxyyzz_0, tg_xyyyy_xxyzzz_0, tg_xyyz_xxxxx_1, \
                                         tg_xyyz_xxxxxx_0, tg_xyyz_xxxxxx_1, tg_xyyz_xxxxxy_0, tg_xyyz_xxxxxy_1, \
                                         tg_xyyz_xxxxxz_0, tg_xyyz_xxxxxz_1, tg_xyyz_xxxxy_1, tg_xyyz_xxxxyy_0, \
                                         tg_xyyz_xxxxyy_1, tg_xyyz_xxxxyz_0, tg_xyyz_xxxxyz_1, tg_xyyz_xxxxz_1, \
                                         tg_xyyz_xxxxzz_0, tg_xyyz_xxxxzz_1, tg_xyyz_xxxyy_1, tg_xyyz_xxxyyy_0, \
                                         tg_xyyz_xxxyyy_1, tg_xyyz_xxxyyz_0, tg_xyyz_xxxyyz_1, tg_xyyz_xxxyz_1, \
                                         tg_xyyz_xxxyzz_0, tg_xyyz_xxxyzz_1, tg_xyyz_xxxzz_1, tg_xyyz_xxxzzz_0, \
                                         tg_xyyz_xxxzzz_1, tg_xyyz_xxyyy_1, tg_xyyz_xxyyyy_0, tg_xyyz_xxyyyy_1, \
                                         tg_xyyz_xxyyyz_0, tg_xyyz_xxyyyz_1, tg_xyyz_xxyyz_1, tg_xyyz_xxyyzz_0, \
                                         tg_xyyz_xxyyzz_1, tg_xyyz_xxyzz_1, tg_xyyz_xxyzzz_0, tg_xyyz_xxyzzz_1, \
                                         tg_xyyz_xxzzz_1, tg_xyyz_xxzzzz_0, tg_xyyz_xxzzzz_1, tg_xyyz_xyyyy_1, \
                                         tg_xyyz_xyyyyy_0, tg_xyyz_xyyyyy_1, tg_xyyz_xyyyyz_0, tg_xyyz_xyyyyz_1, \
                                         tg_xyyz_xyyyz_1, tg_xyyz_xyyyzz_0, tg_xyyz_xyyyzz_1, tg_xyyz_xyyzz_1, \
                                         tg_xyyz_xyyzzz_0, tg_xyyz_xyyzzz_1, tg_xyyz_xyzzz_1, tg_xyyz_xyzzzz_0, \
                                         tg_xyyz_xyzzzz_1, tg_xyyz_xzzzz_1, tg_xyyz_xzzzzz_0, tg_xyyz_xzzzzz_1, \
                                         tg_xyyz_yyyyy_1, tg_xyyz_yyyyyy_0, tg_xyyz_yyyyyy_1, tg_xyyz_yyyyyz_0, \
                                         tg_xyyz_yyyyyz_1, tg_xyyz_yyyyz_1, tg_xyyz_yyyyzz_0, tg_xyyz_yyyyzz_1, \
                                         tg_xyyz_yyyzz_1, tg_xyyz_yyyzzz_0, tg_xyyz_yyyzzz_1, tg_xyyz_yyzzz_1, \
                                         tg_xyyz_yyzzzz_0, tg_xyyz_yyzzzz_1, tg_xyyz_yzzzz_1, tg_xyyz_yzzzzz_0, \
                                         tg_xyyz_yzzzzz_1, tg_xyyz_zzzzz_1, tg_xyyz_zzzzzz_0, tg_xyyz_zzzzzz_1, \
                                         tg_xyzz_xxxxx_1, tg_xyzz_xxxxxx_0, tg_xyzz_xxxxxx_1, tg_xyzz_xxxxxy_0, \
                                         tg_xyzz_xxxxxy_1, tg_xyzz_xxxxxz_0, tg_xyzz_xxxxxz_1, tg_xyzz_xxxxy_1, \
                                         tg_xyzz_xxxxyy_0, tg_xyzz_xxxxyy_1, tg_xyzz_xxxxyz_0, tg_xyzz_xxxxyz_1, \
                                         tg_xyzz_xxxxz_1, tg_xyzz_xxxxzz_0, tg_xyzz_xxxxzz_1, tg_xyzz_xxxyy_1, \
                                         tg_xyzz_xxxyyy_0, tg_xyzz_xxxyyy_1, tg_xyzz_xxxyyz_0, tg_xyzz_xxxyyz_1, \
                                         tg_xyzz_xxxyz_1, tg_xyzz_xxxyzz_0, tg_xyzz_xxxyzz_1, tg_xyzz_xxxzz_1, \
                                         tg_xyzz_xxxzzz_0, tg_xyzz_xxxzzz_1, tg_xyzz_xxyyy_1, tg_xyzz_xxyyyy_0, \
                                         tg_xyzz_xxyyyy_1, tg_xyzz_xxyyyz_0, tg_xyzz_xxyyyz_1, tg_xyzz_xxyyz_1, \
                                         tg_xyzz_xxyyzz_0, tg_xyzz_xxyyzz_1, tg_xyzz_xxyzz_1, tg_xyzz_xxyzzz_0, \
                                         tg_xyzz_xxyzzz_1, tg_xyzz_xxzzz_1, tg_xyzz_xxzzzz_0, tg_xyzz_xxzzzz_1, \
                                         tg_xyzz_xyyyy_1, tg_xyzz_xyyyyy_0, tg_xyzz_xyyyyy_1, tg_xyzz_xyyyyz_0, \
                                         tg_xyzz_xyyyyz_1, tg_xyzz_xyyyz_1, tg_xyzz_xyyyzz_0, tg_xyzz_xyyyzz_1, \
                                         tg_xyzz_xyyzz_1, tg_xyzz_xyyzzz_0, tg_xyzz_xyyzzz_1, tg_xyzz_xyzzz_1, \
                                         tg_xyzz_xyzzzz_0, tg_xyzz_xyzzzz_1, tg_xyzz_xzzzz_1, tg_xyzz_xzzzzz_0, \
                                         tg_xyzz_xzzzzz_1, tg_xyzz_yyyyy_1, tg_xyzz_yyyyyy_0, tg_xyzz_yyyyyy_1, \
                                         tg_xyzz_yyyyyz_0, tg_xyzz_yyyyyz_1, tg_xyzz_yyyyz_1, tg_xyzz_yyyyzz_0, \
                                         tg_xyzz_yyyyzz_1, tg_xyzz_yyyzz_1, tg_xyzz_yyyzzz_0, tg_xyzz_yyyzzz_1, \
                                         tg_xyzz_yyzzz_1, tg_xyzz_yyzzzz_0, tg_xyzz_yyzzzz_1, tg_xyzz_yzzzz_1, \
                                         tg_xyzz_yzzzzz_0, tg_xyzz_yzzzzz_1, tg_xyzz_zzzzz_1, tg_xyzz_zzzzzz_0, \
                                         tg_xyzz_zzzzzz_1, tg_xzzz_xxxxx_1, tg_xzzz_xxxxxx_0, tg_xzzz_xxxxxx_1, \
                                         tg_xzzz_xxxxxy_0, tg_xzzz_xxxxxy_1, tg_xzzz_xxxxxz_0, tg_xzzz_xxxxxz_1, \
                                         tg_xzzz_xxxxy_1, tg_xzzz_xxxxyy_0, tg_xzzz_xxxxyy_1, tg_xzzz_xxxxyz_0, \
                                         tg_xzzz_xxxxyz_1, tg_xzzz_xxxxz_1, tg_xzzz_xxxxzz_0, tg_xzzz_xxxxzz_1, \
                                         tg_xzzz_xxxyy_1, tg_xzzz_xxxyyy_0, tg_xzzz_xxxyyy_1, tg_xzzz_xxxyyz_0, \
                                         tg_xzzz_xxxyyz_1, tg_xzzz_xxxyz_1, tg_xzzz_xxxyzz_0, tg_xzzz_xxxyzz_1, \
                                         tg_xzzz_xxxzz_1, tg_xzzz_xxxzzz_0, tg_xzzz_xxxzzz_1, tg_xzzz_xxyyy_1, \
                                         tg_xzzz_xxyyyy_0, tg_xzzz_xxyyyy_1, tg_xzzz_xxyyyz_0, tg_xzzz_xxyyyz_1, \
                                         tg_xzzz_xxyyz_1, tg_xzzz_xxyyzz_0, tg_xzzz_xxyyzz_1, tg_xzzz_xxyzz_1, \
                                         tg_xzzz_xxyzzz_0, tg_xzzz_xxyzzz_1, tg_xzzz_xxzzz_1, tg_xzzz_xxzzzz_0, \
                                         tg_xzzz_xxzzzz_1, tg_xzzz_xyyyy_1, tg_xzzz_xyyyyy_0, tg_xzzz_xyyyyy_1, \
                                         tg_xzzz_xyyyyz_0, tg_xzzz_xyyyyz_1, tg_xzzz_xyyyz_1, tg_xzzz_xyyyzz_0, \
                                         tg_xzzz_xyyyzz_1, tg_xzzz_xyyzz_1, tg_xzzz_xyyzzz_0, tg_xzzz_xyyzzz_1, \
                                         tg_xzzz_xyzzz_1, tg_xzzz_xyzzzz_0, tg_xzzz_xyzzzz_1, tg_xzzz_xzzzz_1, \
                                         tg_xzzz_xzzzzz_0, tg_xzzz_xzzzzz_1, tg_xzzz_yyyyy_1, tg_xzzz_yyyyyy_0, \
                                         tg_xzzz_yyyyyy_1, tg_xzzz_yyyyyz_0, tg_xzzz_yyyyyz_1, tg_xzzz_yyyyz_1, \
                                         tg_xzzz_yyyyzz_0, tg_xzzz_yyyyzz_1, tg_xzzz_yyyzz_1, tg_xzzz_yyyzzz_0, \
                                         tg_xzzz_yyyzzz_1, tg_xzzz_yyzzz_1, tg_xzzz_yyzzzz_0, tg_xzzz_yyzzzz_1, \
                                         tg_xzzz_yzzzz_1, tg_xzzz_yzzzzz_0, tg_xzzz_yzzzzz_1, tg_xzzz_zzzzz_1, \
                                         tg_xzzz_zzzzzz_0, tg_xzzz_zzzzzz_1, tg_yyyy_xxxxx_1, tg_yyyy_xxxxxx_0, \
                                         tg_yyyy_xxxxxx_1, tg_yyyy_xxxxxy_0, tg_yyyy_xxxxxy_1, tg_yyyy_xxxxxz_0, \
                                         tg_yyyy_xxxxxz_1, tg_yyyy_xxxxy_1, tg_yyyy_xxxxyy_0, tg_yyyy_xxxxyy_1, \
                                         tg_yyyy_xxxxyz_0, tg_yyyy_xxxxyz_1, tg_yyyy_xxxxz_1, tg_yyyy_xxxxzz_0, \
                                         tg_yyyy_xxxxzz_1, tg_yyyy_xxxyy_1, tg_yyyy_xxxyyy_0, tg_yyyy_xxxyyy_1, \
                                         tg_yyyy_xxxyyz_0, tg_yyyy_xxxyyz_1, tg_yyyy_xxxyz_1, tg_yyyy_xxxyzz_0, \
                                         tg_yyyy_xxxyzz_1, tg_yyyy_xxxzz_1, tg_yyyy_xxxzzz_0, tg_yyyy_xxxzzz_1, \
                                         tg_yyyy_xxyyy_1, tg_yyyy_xxyyyy_0, tg_yyyy_xxyyyy_1, tg_yyyy_xxyyyz_0, \
                                         tg_yyyy_xxyyyz_1, tg_yyyy_xxyyz_1, tg_yyyy_xxyyzz_0, tg_yyyy_xxyyzz_1, \
                                         tg_yyyy_xxyzz_1, tg_yyyy_xxyzzz_0, tg_yyyy_xxyzzz_1, tg_yyyy_xxzzz_1, \
                                         tg_yyyy_xyyyy_1, tg_yyyy_xyyyz_1, tg_yyyy_xyyzz_1, tg_yyyy_xyzzz_1, tg_yyz_xxxxxx_0, \
                                         tg_yyz_xxxxxx_1, tg_yyz_xxxxxy_0, tg_yyz_xxxxxy_1, tg_yyz_xxxxxz_0, tg_yyz_xxxxxz_1, \
                                         tg_yyz_xxxxyy_0, tg_yyz_xxxxyy_1, tg_yyz_xxxxyz_0, tg_yyz_xxxxyz_1, tg_yyz_xxxxzz_0, \
                                         tg_yyz_xxxxzz_1, tg_yyz_xxxyyy_0, tg_yyz_xxxyyy_1, tg_yyz_xxxyyz_0, tg_yyz_xxxyyz_1, \
                                         tg_yyz_xxxyzz_0, tg_yyz_xxxyzz_1, tg_yyz_xxxzzz_0, tg_yyz_xxxzzz_1, tg_yyz_xxyyyy_0, \
                                         tg_yyz_xxyyyy_1, tg_yyz_xxyyyz_0, tg_yyz_xxyyyz_1, tg_yyz_xxyyzz_0, tg_yyz_xxyyzz_1, \
                                         tg_yyz_xxyzzz_0, tg_yyz_xxyzzz_1, tg_yyz_xxzzzz_0, tg_yyz_xxzzzz_1, tg_yyz_xyyyyy_0, \
                                         tg_yyz_xyyyyy_1, tg_yyz_xyyyyz_0, tg_yyz_xyyyyz_1, tg_yyz_xyyyzz_0, tg_yyz_xyyyzz_1, \
                                         tg_yyz_xyyzzz_0, tg_yyz_xyyzzz_1, tg_yyz_xyzzzz_0, tg_yyz_xyzzzz_1, tg_yyz_xzzzzz_0, \
                                         tg_yyz_xzzzzz_1, tg_yyz_yyyyyy_0, tg_yyz_yyyyyy_1, tg_yyz_yyyyyz_0, tg_yyz_yyyyyz_1, \
                                         tg_yyz_yyyyzz_0, tg_yyz_yyyyzz_1, tg_yyz_yyyzzz_0, tg_yyz_yyyzzz_1, tg_yyz_yyzzzz_0, \
                                         tg_yyz_yyzzzz_1, tg_yyz_yzzzzz_0, tg_yyz_yzzzzz_1, tg_yyz_zzzzzz_0, tg_yyz_zzzzzz_1, \
                                         tg_yzz_xxxxxx_0, tg_yzz_xxxxxx_1, tg_yzz_xxxxxy_0, tg_yzz_xxxxxy_1, tg_yzz_xxxxxz_0, \
                                         tg_yzz_xxxxxz_1, tg_yzz_xxxxyy_0, tg_yzz_xxxxyy_1, tg_yzz_xxxxyz_0, tg_yzz_xxxxyz_1, \
                                         tg_yzz_xxxxzz_0, tg_yzz_xxxxzz_1, tg_yzz_xxxyyy_0, tg_yzz_xxxyyy_1, tg_yzz_xxxyyz_0, \
                                         tg_yzz_xxxyyz_1, tg_yzz_xxxyzz_0, tg_yzz_xxxyzz_1, tg_yzz_xxxzzz_0, tg_yzz_xxxzzz_1, \
                                         tg_yzz_xxyyyy_0, tg_yzz_xxyyyy_1, tg_yzz_xxyyyz_0, tg_yzz_xxyyyz_1, tg_yzz_xxyyzz_0, \
                                         tg_yzz_xxyyzz_1, tg_yzz_xxyzzz_0, tg_yzz_xxyzzz_1, tg_yzz_xxzzzz_0, tg_yzz_xxzzzz_1, \
                                         tg_yzz_xyyyyy_0, tg_yzz_xyyyyy_1, tg_yzz_xyyyyz_0, tg_yzz_xyyyyz_1, tg_yzz_xyyyzz_0, \
                                         tg_yzz_xyyyzz_1, tg_yzz_xyyzzz_0, tg_yzz_xyyzzz_1, tg_yzz_xyzzzz_0, tg_yzz_xyzzzz_1, \
                                         tg_yzz_xzzzzz_0, tg_yzz_xzzzzz_1, tg_yzz_yyyyyy_0, tg_yzz_yyyyyy_1, tg_yzz_yyyyyz_0, \
                                         tg_yzz_yyyyyz_1, tg_yzz_yyyyzz_0, tg_yzz_yyyyzz_1, tg_yzz_yyyzzz_0, tg_yzz_yyyzzz_1, \
                                         tg_yzz_yyzzzz_0, tg_yzz_yyzzzz_1, tg_yzz_yzzzzz_0, tg_yzz_yzzzzz_1, tg_yzz_zzzzzz_0, \
                                         tg_yzz_zzzzzz_1, tg_zzz_xxxxxx_0, tg_zzz_xxxxxx_1, tg_zzz_xxxxxy_0, tg_zzz_xxxxxy_1, \
                                         tg_zzz_xxxxxz_0, tg_zzz_xxxxxz_1, tg_zzz_xxxxyy_0, tg_zzz_xxxxyy_1, tg_zzz_xxxxyz_0, \
                                         tg_zzz_xxxxyz_1, tg_zzz_xxxxzz_0, tg_zzz_xxxxzz_1, tg_zzz_xxxyyy_0, tg_zzz_xxxyyy_1, \
                                         tg_zzz_xxxyyz_0, tg_zzz_xxxyyz_1, tg_zzz_xxxyzz_0, tg_zzz_xxxyzz_1, tg_zzz_xxxzzz_0, \
                                         tg_zzz_xxxzzz_1, tg_zzz_xxyyyy_0, tg_zzz_xxyyyy_1, tg_zzz_xxyyyz_0, tg_zzz_xxyyyz_1, \
                                         tg_zzz_xxyyzz_0, tg_zzz_xxyyzz_1, tg_zzz_xxyzzz_0, tg_zzz_xxyzzz_1, tg_zzz_xxzzzz_0, \
                                         tg_zzz_xxzzzz_1, tg_zzz_xyyyyy_0, tg_zzz_xyyyyy_1, tg_zzz_xyyyyz_0, tg_zzz_xyyyyz_1, \
                                         tg_zzz_xyyyzz_0, tg_zzz_xyyyzz_1, tg_zzz_xyyzzz_0, tg_zzz_xyyzzz_1, tg_zzz_xyzzzz_0, \
                                         tg_zzz_xyzzzz_1, tg_zzz_xzzzzz_0, tg_zzz_xzzzzz_1, tg_zzz_yyyyyy_0, tg_zzz_yyyyyy_1, \
                                         tg_zzz_yyyyyz_0, tg_zzz_yyyyyz_1, tg_zzz_yyyyzz_0, tg_zzz_yyyyzz_1, tg_zzz_yyyzzz_0, \
                                         tg_zzz_yyyzzz_1, tg_zzz_yyzzzz_0, tg_zzz_yyzzzz_1, tg_zzz_yzzzzz_0, tg_zzz_yzzzzz_1, \
                                         tg_zzz_zzzzzz_0, tg_zzz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyz_xxxxxx_0[j] = pb_x * tg_xyyz_xxxxxx_0[j] + fr * tg_xyyz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxx_0[j] - tg_yyz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyyz_xxxxx_1[j];

                    tg_xxyyz_xxxxxy_0[j] = pb_x * tg_xyyz_xxxxxy_0[j] + fr * tg_xyyz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxy_0[j] - tg_yyz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyz_xxxxy_1[j];

                    tg_xxyyz_xxxxxz_0[j] = pb_x * tg_xyyz_xxxxxz_0[j] + fr * tg_xyyz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxxz_0[j] - tg_yyz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyyz_xxxxz_1[j];

                    tg_xxyyz_xxxxyy_0[j] = pb_x * tg_xyyz_xxxxyy_0[j] + fr * tg_xyyz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxyy_0[j] - tg_yyz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxyy_1[j];

                    tg_xxyyz_xxxxyz_0[j] = pb_x * tg_xyyz_xxxxyz_0[j] + fr * tg_xyyz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxyz_0[j] - tg_yyz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxyz_1[j];

                    tg_xxyyz_xxxxzz_0[j] = pb_x * tg_xyyz_xxxxzz_0[j] + fr * tg_xyyz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxxzz_0[j] - tg_yyz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyz_xxxzz_1[j];

                    tg_xxyyz_xxxyyy_0[j] = pb_x * tg_xyyz_xxxyyy_0[j] + fr * tg_xyyz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyyy_0[j] - tg_yyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyyy_1[j];

                    tg_xxyyz_xxxyyz_0[j] = pb_x * tg_xyyz_xxxyyz_0[j] + fr * tg_xyyz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyyz_0[j] - tg_yyz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyyz_1[j];

                    tg_xxyyz_xxxyzz_0[j] = pb_x * tg_xyyz_xxxyzz_0[j] + fr * tg_xyyz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxyzz_0[j] - tg_yyz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxyzz_1[j];

                    tg_xxyyz_xxxzzz_0[j] = pb_x * tg_xyyz_xxxzzz_0[j] + fr * tg_xyyz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxxzzz_0[j] - tg_yyz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyz_xxzzz_1[j];

                    tg_xxyyz_xxyyyy_0[j] = pb_x * tg_xyyz_xxyyyy_0[j] + fr * tg_xyyz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyyy_0[j] - tg_yyz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyyy_1[j];

                    tg_xxyyz_xxyyyz_0[j] = pb_x * tg_xyyz_xxyyyz_0[j] + fr * tg_xyyz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyyz_0[j] - tg_yyz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyyz_1[j];

                    tg_xxyyz_xxyyzz_0[j] = pb_x * tg_xyyz_xxyyzz_0[j] + fr * tg_xyyz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyyzz_0[j] - tg_yyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyyzz_1[j];

                    tg_xxyyz_xxyzzz_0[j] = pb_x * tg_xyyz_xxyzzz_0[j] + fr * tg_xyyz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxyzzz_0[j] - tg_yyz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xyzzz_1[j];

                    tg_xxyyz_xxzzzz_0[j] = pb_x * tg_xyyz_xxzzzz_0[j] + fr * tg_xyyz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xxzzzz_0[j] - tg_yyz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyz_xzzzz_1[j];

                    tg_xxyyz_xyyyyy_0[j] = pb_x * tg_xyyz_xyyyyy_0[j] + fr * tg_xyyz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyyy_0[j] - tg_yyz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyyy_1[j];

                    tg_xxyyz_xyyyyz_0[j] = pb_x * tg_xyyz_xyyyyz_0[j] + fr * tg_xyyz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyyz_0[j] - tg_yyz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyyz_1[j];

                    tg_xxyyz_xyyyzz_0[j] = pb_x * tg_xyyz_xyyyzz_0[j] + fr * tg_xyyz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyyzz_0[j] - tg_yyz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyyzz_1[j];

                    tg_xxyyz_xyyzzz_0[j] = pb_x * tg_xyyz_xyyzzz_0[j] + fr * tg_xyyz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyyzzz_0[j] - tg_yyz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yyzzz_1[j];

                    tg_xxyyz_xyzzzz_0[j] = pb_x * tg_xyyz_xyzzzz_0[j] + fr * tg_xyyz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xyzzzz_0[j] - tg_yyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_yzzzz_1[j];

                    tg_xxyyz_xzzzzz_0[j] = pb_x * tg_xyyz_xzzzzz_0[j] + fr * tg_xyyz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_xzzzzz_0[j] - tg_yyz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyz_zzzzz_1[j];

                    tg_xxyyz_yyyyyy_0[j] = pb_x * tg_xyyz_yyyyyy_0[j] + fr * tg_xyyz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyyy_0[j] - tg_yyz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyyz_yyyyyz_0[j] = pb_x * tg_xyyz_yyyyyz_0[j] + fr * tg_xyyz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyyz_0[j] - tg_yyz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyyz_yyyyzz_0[j] = pb_x * tg_xyyz_yyyyzz_0[j] + fr * tg_xyyz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyyzz_0[j] - tg_yyz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyyz_yyyzzz_0[j] = pb_x * tg_xyyz_yyyzzz_0[j] + fr * tg_xyyz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyyzzz_0[j] - tg_yyz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyyz_yyzzzz_0[j] = pb_x * tg_xyyz_yyzzzz_0[j] + fr * tg_xyyz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yyzzzz_0[j] - tg_yyz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyyz_yzzzzz_0[j] = pb_x * tg_xyyz_yzzzzz_0[j] + fr * tg_xyyz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_yzzzzz_0[j] - tg_yyz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyyz_zzzzzz_0[j] = pb_x * tg_xyyz_zzzzzz_0[j] + fr * tg_xyyz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yyz_zzzzzz_0[j] - tg_yyz_zzzzzz_1[j] * fl1_fza);

                    tg_xxyzz_xxxxxx_0[j] = pb_x * tg_xyzz_xxxxxx_0[j] + fr * tg_xyzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxx_0[j] - tg_yzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xyzz_xxxxx_1[j];

                    tg_xxyzz_xxxxxy_0[j] = pb_x * tg_xyzz_xxxxxy_0[j] + fr * tg_xyzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxy_0[j] - tg_yzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzz_xxxxy_1[j];

                    tg_xxyzz_xxxxxz_0[j] = pb_x * tg_xyzz_xxxxxz_0[j] + fr * tg_xyzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxxz_0[j] - tg_yzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xyzz_xxxxz_1[j];

                    tg_xxyzz_xxxxyy_0[j] = pb_x * tg_xyzz_xxxxyy_0[j] + fr * tg_xyzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxyy_0[j] - tg_yzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxyy_1[j];

                    tg_xxyzz_xxxxyz_0[j] = pb_x * tg_xyzz_xxxxyz_0[j] + fr * tg_xyzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxyz_0[j] - tg_yzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxyz_1[j];

                    tg_xxyzz_xxxxzz_0[j] = pb_x * tg_xyzz_xxxxzz_0[j] + fr * tg_xyzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxxzz_0[j] - tg_yzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzz_xxxzz_1[j];

                    tg_xxyzz_xxxyyy_0[j] = pb_x * tg_xyzz_xxxyyy_0[j] + fr * tg_xyzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyyy_0[j] - tg_yzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyyy_1[j];

                    tg_xxyzz_xxxyyz_0[j] = pb_x * tg_xyzz_xxxyyz_0[j] + fr * tg_xyzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyyz_0[j] - tg_yzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyyz_1[j];

                    tg_xxyzz_xxxyzz_0[j] = pb_x * tg_xyzz_xxxyzz_0[j] + fr * tg_xyzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxyzz_0[j] - tg_yzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxyzz_1[j];

                    tg_xxyzz_xxxzzz_0[j] = pb_x * tg_xyzz_xxxzzz_0[j] + fr * tg_xyzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxxzzz_0[j] - tg_yzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzz_xxzzz_1[j];

                    tg_xxyzz_xxyyyy_0[j] = pb_x * tg_xyzz_xxyyyy_0[j] + fr * tg_xyzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyyy_0[j] - tg_yzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyyy_1[j];

                    tg_xxyzz_xxyyyz_0[j] = pb_x * tg_xyzz_xxyyyz_0[j] + fr * tg_xyzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyyz_0[j] - tg_yzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyyz_1[j];

                    tg_xxyzz_xxyyzz_0[j] = pb_x * tg_xyzz_xxyyzz_0[j] + fr * tg_xyzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyyzz_0[j] - tg_yzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyyzz_1[j];

                    tg_xxyzz_xxyzzz_0[j] = pb_x * tg_xyzz_xxyzzz_0[j] + fr * tg_xyzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxyzzz_0[j] - tg_yzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xyzzz_1[j];

                    tg_xxyzz_xxzzzz_0[j] = pb_x * tg_xyzz_xxzzzz_0[j] + fr * tg_xyzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xxzzzz_0[j] - tg_yzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzz_xzzzz_1[j];

                    tg_xxyzz_xyyyyy_0[j] = pb_x * tg_xyzz_xyyyyy_0[j] + fr * tg_xyzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyyy_0[j] - tg_yzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyyy_1[j];

                    tg_xxyzz_xyyyyz_0[j] = pb_x * tg_xyzz_xyyyyz_0[j] + fr * tg_xyzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyyz_0[j] - tg_yzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyyz_1[j];

                    tg_xxyzz_xyyyzz_0[j] = pb_x * tg_xyzz_xyyyzz_0[j] + fr * tg_xyzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyyzz_0[j] - tg_yzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyyzz_1[j];

                    tg_xxyzz_xyyzzz_0[j] = pb_x * tg_xyzz_xyyzzz_0[j] + fr * tg_xyzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyyzzz_0[j] - tg_yzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yyzzz_1[j];

                    tg_xxyzz_xyzzzz_0[j] = pb_x * tg_xyzz_xyzzzz_0[j] + fr * tg_xyzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xyzzzz_0[j] - tg_yzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_yzzzz_1[j];

                    tg_xxyzz_xzzzzz_0[j] = pb_x * tg_xyzz_xzzzzz_0[j] + fr * tg_xyzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_xzzzzz_0[j] - tg_yzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzz_zzzzz_1[j];

                    tg_xxyzz_yyyyyy_0[j] = pb_x * tg_xyzz_yyyyyy_0[j] + fr * tg_xyzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyyy_0[j] - tg_yzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxyzz_yyyyyz_0[j] = pb_x * tg_xyzz_yyyyyz_0[j] + fr * tg_xyzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyyz_0[j] - tg_yzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxyzz_yyyyzz_0[j] = pb_x * tg_xyzz_yyyyzz_0[j] + fr * tg_xyzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyyzz_0[j] - tg_yzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxyzz_yyyzzz_0[j] = pb_x * tg_xyzz_yyyzzz_0[j] + fr * tg_xyzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyyzzz_0[j] - tg_yzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxyzz_yyzzzz_0[j] = pb_x * tg_xyzz_yyzzzz_0[j] + fr * tg_xyzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yyzzzz_0[j] - tg_yzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxyzz_yzzzzz_0[j] = pb_x * tg_xyzz_yzzzzz_0[j] + fr * tg_xyzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_yzzzzz_0[j] - tg_yzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxyzz_zzzzzz_0[j] = pb_x * tg_xyzz_zzzzzz_0[j] + fr * tg_xyzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_yzz_zzzzzz_0[j] - tg_yzz_zzzzzz_1[j] * fl1_fza);

                    tg_xxzzz_xxxxxx_0[j] = pb_x * tg_xzzz_xxxxxx_0[j] + fr * tg_xzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxx_0[j] - tg_zzz_xxxxxx_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_xzzz_xxxxx_1[j];

                    tg_xxzzz_xxxxxy_0[j] = pb_x * tg_xzzz_xxxxxy_0[j] + fr * tg_xzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxy_0[j] - tg_zzz_xxxxxy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzz_xxxxy_1[j];

                    tg_xxzzz_xxxxxz_0[j] = pb_x * tg_xzzz_xxxxxz_0[j] + fr * tg_xzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxz_0[j] - tg_zzz_xxxxxz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_xzzz_xxxxz_1[j];

                    tg_xxzzz_xxxxyy_0[j] = pb_x * tg_xzzz_xxxxyy_0[j] + fr * tg_xzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyy_0[j] - tg_zzz_xxxxyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxyy_1[j];

                    tg_xxzzz_xxxxyz_0[j] = pb_x * tg_xzzz_xxxxyz_0[j] + fr * tg_xzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyz_0[j] - tg_zzz_xxxxyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxyz_1[j];

                    tg_xxzzz_xxxxzz_0[j] = pb_x * tg_xzzz_xxxxzz_0[j] + fr * tg_xzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxzz_0[j] - tg_zzz_xxxxzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzz_xxxzz_1[j];

                    tg_xxzzz_xxxyyy_0[j] = pb_x * tg_xzzz_xxxyyy_0[j] + fr * tg_xzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyy_0[j] - tg_zzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyyy_1[j];

                    tg_xxzzz_xxxyyz_0[j] = pb_x * tg_xzzz_xxxyyz_0[j] + fr * tg_xzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyz_0[j] - tg_zzz_xxxyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyyz_1[j];

                    tg_xxzzz_xxxyzz_0[j] = pb_x * tg_xzzz_xxxyzz_0[j] + fr * tg_xzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyzz_0[j] - tg_zzz_xxxyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxyzz_1[j];

                    tg_xxzzz_xxxzzz_0[j] = pb_x * tg_xzzz_xxxzzz_0[j] + fr * tg_xzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxzzz_0[j] - tg_zzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzz_xxzzz_1[j];

                    tg_xxzzz_xxyyyy_0[j] = pb_x * tg_xzzz_xxyyyy_0[j] + fr * tg_xzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyy_0[j] - tg_zzz_xxyyyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyyy_1[j];

                    tg_xxzzz_xxyyyz_0[j] = pb_x * tg_xzzz_xxyyyz_0[j] + fr * tg_xzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyz_0[j] - tg_zzz_xxyyyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyyz_1[j];

                    tg_xxzzz_xxyyzz_0[j] = pb_x * tg_xzzz_xxyyzz_0[j] + fr * tg_xzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyzz_0[j] - tg_zzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyyzz_1[j];

                    tg_xxzzz_xxyzzz_0[j] = pb_x * tg_xzzz_xxyzzz_0[j] + fr * tg_xzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyzzz_0[j] - tg_zzz_xxyzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xyzzz_1[j];

                    tg_xxzzz_xxzzzz_0[j] = pb_x * tg_xzzz_xxzzzz_0[j] + fr * tg_xzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxzzzz_0[j] - tg_zzz_xxzzzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzz_xzzzz_1[j];

                    tg_xxzzz_xyyyyy_0[j] = pb_x * tg_xzzz_xyyyyy_0[j] + fr * tg_xzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyy_0[j] - tg_zzz_xyyyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyyy_1[j];

                    tg_xxzzz_xyyyyz_0[j] = pb_x * tg_xzzz_xyyyyz_0[j] + fr * tg_xzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyz_0[j] - tg_zzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyyz_1[j];

                    tg_xxzzz_xyyyzz_0[j] = pb_x * tg_xzzz_xyyyzz_0[j] + fr * tg_xzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyzz_0[j] - tg_zzz_xyyyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyyzz_1[j];

                    tg_xxzzz_xyyzzz_0[j] = pb_x * tg_xzzz_xyyzzz_0[j] + fr * tg_xzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyzzz_0[j] - tg_zzz_xyyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yyzzz_1[j];

                    tg_xxzzz_xyzzzz_0[j] = pb_x * tg_xzzz_xyzzzz_0[j] + fr * tg_xzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyzzzz_0[j] - tg_zzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_yzzzz_1[j];

                    tg_xxzzz_xzzzzz_0[j] = pb_x * tg_xzzz_xzzzzz_0[j] + fr * tg_xzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xzzzzz_0[j] - tg_zzz_xzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzz_zzzzz_1[j];

                    tg_xxzzz_yyyyyy_0[j] = pb_x * tg_xzzz_yyyyyy_0[j] + fr * tg_xzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyy_0[j] - tg_zzz_yyyyyy_1[j] * fl1_fza);

                    tg_xxzzz_yyyyyz_0[j] = pb_x * tg_xzzz_yyyyyz_0[j] + fr * tg_xzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyz_0[j] - tg_zzz_yyyyyz_1[j] * fl1_fza);

                    tg_xxzzz_yyyyzz_0[j] = pb_x * tg_xzzz_yyyyzz_0[j] + fr * tg_xzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyzz_0[j] - tg_zzz_yyyyzz_1[j] * fl1_fza);

                    tg_xxzzz_yyyzzz_0[j] = pb_x * tg_xzzz_yyyzzz_0[j] + fr * tg_xzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyzzz_0[j] - tg_zzz_yyyzzz_1[j] * fl1_fza);

                    tg_xxzzz_yyzzzz_0[j] = pb_x * tg_xzzz_yyzzzz_0[j] + fr * tg_xzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyzzzz_0[j] - tg_zzz_yyzzzz_1[j] * fl1_fza);

                    tg_xxzzz_yzzzzz_0[j] = pb_x * tg_xzzz_yzzzzz_0[j] + fr * tg_xzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yzzzzz_0[j] - tg_zzz_yzzzzz_1[j] * fl1_fza);

                    tg_xxzzz_zzzzzz_0[j] = pb_x * tg_xzzz_zzzzzz_0[j] + fr * tg_xzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_zzzzzz_0[j] - tg_zzz_zzzzzz_1[j] * fl1_fza);

                    tg_xyyyy_xxxxxx_0[j] = pb_x * tg_yyyy_xxxxxx_0[j] + fr * tg_yyyy_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyy_xxxxx_1[j];

                    tg_xyyyy_xxxxxy_0[j] = pb_x * tg_yyyy_xxxxxy_0[j] + fr * tg_yyyy_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxy_1[j];

                    tg_xyyyy_xxxxxz_0[j] = pb_x * tg_yyyy_xxxxxz_0[j] + fr * tg_yyyy_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyy_xxxxz_1[j];

                    tg_xyyyy_xxxxyy_0[j] = pb_x * tg_yyyy_xxxxyy_0[j] + fr * tg_yyyy_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyy_1[j];

                    tg_xyyyy_xxxxyz_0[j] = pb_x * tg_yyyy_xxxxyz_0[j] + fr * tg_yyyy_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxyz_1[j];

                    tg_xyyyy_xxxxzz_0[j] = pb_x * tg_yyyy_xxxxzz_0[j] + fr * tg_yyyy_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxxzz_1[j];

                    tg_xyyyy_xxxyyy_0[j] = pb_x * tg_yyyy_xxxyyy_0[j] + fr * tg_yyyy_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyy_1[j];

                    tg_xyyyy_xxxyyz_0[j] = pb_x * tg_yyyy_xxxyyz_0[j] + fr * tg_yyyy_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyyz_1[j];

                    tg_xyyyy_xxxyzz_0[j] = pb_x * tg_yyyy_xxxyzz_0[j] + fr * tg_yyyy_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxyzz_1[j];

                    tg_xyyyy_xxxzzz_0[j] = pb_x * tg_yyyy_xxxzzz_0[j] + fr * tg_yyyy_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxzzz_1[j];

                    tg_xyyyy_xxyyyy_0[j] = pb_x * tg_yyyy_xxyyyy_0[j] + fr * tg_yyyy_xxyyyy_1[j] + fl1_fxn * tg_yyyy_xyyyy_1[j];

                    tg_xyyyy_xxyyyz_0[j] = pb_x * tg_yyyy_xxyyyz_0[j] + fr * tg_yyyy_xxyyyz_1[j] + fl1_fxn * tg_yyyy_xyyyz_1[j];

                    tg_xyyyy_xxyyzz_0[j] = pb_x * tg_yyyy_xxyyzz_0[j] + fr * tg_yyyy_xxyyzz_1[j] + fl1_fxn * tg_yyyy_xyyzz_1[j];

                    tg_xyyyy_xxyzzz_0[j] = pb_x * tg_yyyy_xxyzzz_0[j] + fr * tg_yyyy_xxyzzz_1[j] + fl1_fxn * tg_yyyy_xyzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSI_294_392(      CMemBlock2D<double>* primBuffer,
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

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
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

                // set up pointers to integrals

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

                // Batch of Integrals (294,392)

                #pragma omp simd aligned(fxn, tg_xyyyy_xxzzzz_0, tg_xyyyy_xyyyyy_0, tg_xyyyy_xyyyyz_0, \
                                         tg_xyyyy_xyyyzz_0, tg_xyyyy_xyyzzz_0, tg_xyyyy_xyzzzz_0, tg_xyyyy_xzzzzz_0, \
                                         tg_xyyyy_yyyyyy_0, tg_xyyyy_yyyyyz_0, tg_xyyyy_yyyyzz_0, tg_xyyyy_yyyzzz_0, \
                                         tg_xyyyy_yyzzzz_0, tg_xyyyy_yzzzzz_0, tg_xyyyy_zzzzzz_0, tg_xyyyz_xxxxxx_0, \
                                         tg_xyyyz_xxxxxy_0, tg_xyyyz_xxxxxz_0, tg_xyyyz_xxxxyy_0, tg_xyyyz_xxxxyz_0, \
                                         tg_xyyyz_xxxxzz_0, tg_xyyyz_xxxyyy_0, tg_xyyyz_xxxyyz_0, tg_xyyyz_xxxyzz_0, \
                                         tg_xyyyz_xxxzzz_0, tg_xyyyz_xxyyyy_0, tg_xyyyz_xxyyyz_0, tg_xyyyz_xxyyzz_0, \
                                         tg_xyyyz_xxyzzz_0, tg_xyyyz_xxzzzz_0, tg_xyyyz_xyyyyy_0, tg_xyyyz_xyyyyz_0, \
                                         tg_xyyyz_xyyyzz_0, tg_xyyyz_xyyzzz_0, tg_xyyyz_xyzzzz_0, tg_xyyyz_xzzzzz_0, \
                                         tg_xyyyz_yyyyyy_0, tg_xyyyz_yyyyyz_0, tg_xyyyz_yyyyzz_0, tg_xyyyz_yyyzzz_0, \
                                         tg_xyyyz_yyzzzz_0, tg_xyyyz_yzzzzz_0, tg_xyyyz_zzzzzz_0, tg_xyyzz_xxxxxx_0, \
                                         tg_xyyzz_xxxxxy_0, tg_xyyzz_xxxxxz_0, tg_xyyzz_xxxxyy_0, tg_xyyzz_xxxxyz_0, \
                                         tg_xyyzz_xxxxzz_0, tg_xyyzz_xxxyyy_0, tg_xyyzz_xxxyyz_0, tg_xyyzz_xxxyzz_0, \
                                         tg_xyyzz_xxxzzz_0, tg_xyyzz_xxyyyy_0, tg_xyyzz_xxyyyz_0, tg_xyyzz_xxyyzz_0, \
                                         tg_xyyzz_xxyzzz_0, tg_xyyzz_xxzzzz_0, tg_xyyzz_xyyyyy_0, tg_xyyzz_xyyyyz_0, \
                                         tg_xyyzz_xyyyzz_0, tg_xyyzz_xyyzzz_0, tg_xyyzz_xyzzzz_0, tg_xyyzz_xzzzzz_0, \
                                         tg_xyyzz_yyyyyy_0, tg_xyyzz_yyyyyz_0, tg_xyyzz_yyyyzz_0, tg_xyyzz_yyyzzz_0, \
                                         tg_xyyzz_yyzzzz_0, tg_xyyzz_yzzzzz_0, tg_xyyzz_zzzzzz_0, tg_xyzzz_xxxxxx_0, \
                                         tg_xyzzz_xxxxxy_0, tg_xyzzz_xxxxxz_0, tg_xyzzz_xxxxyy_0, tg_xyzzz_xxxxyz_0, \
                                         tg_xyzzz_xxxxzz_0, tg_xyzzz_xxxyyy_0, tg_xyzzz_xxxyyz_0, tg_xyzzz_xxxyzz_0, \
                                         tg_xyzzz_xxxzzz_0, tg_xyzzz_xxyyyy_0, tg_xyzzz_xxyyyz_0, tg_xyzzz_xxyyzz_0, \
                                         tg_xyzzz_xxyzzz_0, tg_xyzzz_xxzzzz_0, tg_xyzzz_xyyyyy_0, tg_xyzzz_xyyyyz_0, \
                                         tg_xyzzz_xyyyzz_0, tg_xyzzz_xyyzzz_0, tg_xyzzz_xyzzzz_0, tg_xyzzz_xzzzzz_0, \
                                         tg_xyzzz_yyyyyy_0, tg_xyzzz_yyyyyz_0, tg_xyzzz_yyyyzz_0, tg_xyzzz_yyyzzz_0, \
                                         tg_xyzzz_yyzzzz_0, tg_xyzzz_yzzzzz_0, tg_xyzzz_zzzzzz_0, tg_yyyy_xxzzzz_0, \
                                         tg_yyyy_xxzzzz_1, tg_yyyy_xyyyyy_0, tg_yyyy_xyyyyy_1, tg_yyyy_xyyyyz_0, \
                                         tg_yyyy_xyyyyz_1, tg_yyyy_xyyyzz_0, tg_yyyy_xyyyzz_1, tg_yyyy_xyyzzz_0, \
                                         tg_yyyy_xyyzzz_1, tg_yyyy_xyzzzz_0, tg_yyyy_xyzzzz_1, tg_yyyy_xzzzz_1, \
                                         tg_yyyy_xzzzzz_0, tg_yyyy_xzzzzz_1, tg_yyyy_yyyyy_1, tg_yyyy_yyyyyy_0, \
                                         tg_yyyy_yyyyyy_1, tg_yyyy_yyyyyz_0, tg_yyyy_yyyyyz_1, tg_yyyy_yyyyz_1, \
                                         tg_yyyy_yyyyzz_0, tg_yyyy_yyyyzz_1, tg_yyyy_yyyzz_1, tg_yyyy_yyyzzz_0, \
                                         tg_yyyy_yyyzzz_1, tg_yyyy_yyzzz_1, tg_yyyy_yyzzzz_0, tg_yyyy_yyzzzz_1, \
                                         tg_yyyy_yzzzz_1, tg_yyyy_yzzzzz_0, tg_yyyy_yzzzzz_1, tg_yyyy_zzzzz_1, \
                                         tg_yyyy_zzzzzz_0, tg_yyyy_zzzzzz_1, tg_yyyz_xxxxx_1, tg_yyyz_xxxxxx_0, \
                                         tg_yyyz_xxxxxx_1, tg_yyyz_xxxxxy_0, tg_yyyz_xxxxxy_1, tg_yyyz_xxxxxz_0, \
                                         tg_yyyz_xxxxxz_1, tg_yyyz_xxxxy_1, tg_yyyz_xxxxyy_0, tg_yyyz_xxxxyy_1, \
                                         tg_yyyz_xxxxyz_0, tg_yyyz_xxxxyz_1, tg_yyyz_xxxxz_1, tg_yyyz_xxxxzz_0, \
                                         tg_yyyz_xxxxzz_1, tg_yyyz_xxxyy_1, tg_yyyz_xxxyyy_0, tg_yyyz_xxxyyy_1, \
                                         tg_yyyz_xxxyyz_0, tg_yyyz_xxxyyz_1, tg_yyyz_xxxyz_1, tg_yyyz_xxxyzz_0, \
                                         tg_yyyz_xxxyzz_1, tg_yyyz_xxxzz_1, tg_yyyz_xxxzzz_0, tg_yyyz_xxxzzz_1, \
                                         tg_yyyz_xxyyy_1, tg_yyyz_xxyyyy_0, tg_yyyz_xxyyyy_1, tg_yyyz_xxyyyz_0, \
                                         tg_yyyz_xxyyyz_1, tg_yyyz_xxyyz_1, tg_yyyz_xxyyzz_0, tg_yyyz_xxyyzz_1, \
                                         tg_yyyz_xxyzz_1, tg_yyyz_xxyzzz_0, tg_yyyz_xxyzzz_1, tg_yyyz_xxzzz_1, \
                                         tg_yyyz_xxzzzz_0, tg_yyyz_xxzzzz_1, tg_yyyz_xyyyy_1, tg_yyyz_xyyyyy_0, \
                                         tg_yyyz_xyyyyy_1, tg_yyyz_xyyyyz_0, tg_yyyz_xyyyyz_1, tg_yyyz_xyyyz_1, \
                                         tg_yyyz_xyyyzz_0, tg_yyyz_xyyyzz_1, tg_yyyz_xyyzz_1, tg_yyyz_xyyzzz_0, \
                                         tg_yyyz_xyyzzz_1, tg_yyyz_xyzzz_1, tg_yyyz_xyzzzz_0, tg_yyyz_xyzzzz_1, \
                                         tg_yyyz_xzzzz_1, tg_yyyz_xzzzzz_0, tg_yyyz_xzzzzz_1, tg_yyyz_yyyyy_1, \
                                         tg_yyyz_yyyyyy_0, tg_yyyz_yyyyyy_1, tg_yyyz_yyyyyz_0, tg_yyyz_yyyyyz_1, \
                                         tg_yyyz_yyyyz_1, tg_yyyz_yyyyzz_0, tg_yyyz_yyyyzz_1, tg_yyyz_yyyzz_1, \
                                         tg_yyyz_yyyzzz_0, tg_yyyz_yyyzzz_1, tg_yyyz_yyzzz_1, tg_yyyz_yyzzzz_0, \
                                         tg_yyyz_yyzzzz_1, tg_yyyz_yzzzz_1, tg_yyyz_yzzzzz_0, tg_yyyz_yzzzzz_1, \
                                         tg_yyyz_zzzzz_1, tg_yyyz_zzzzzz_0, tg_yyyz_zzzzzz_1, tg_yyzz_xxxxx_1, \
                                         tg_yyzz_xxxxxx_0, tg_yyzz_xxxxxx_1, tg_yyzz_xxxxxy_0, tg_yyzz_xxxxxy_1, \
                                         tg_yyzz_xxxxxz_0, tg_yyzz_xxxxxz_1, tg_yyzz_xxxxy_1, tg_yyzz_xxxxyy_0, \
                                         tg_yyzz_xxxxyy_1, tg_yyzz_xxxxyz_0, tg_yyzz_xxxxyz_1, tg_yyzz_xxxxz_1, \
                                         tg_yyzz_xxxxzz_0, tg_yyzz_xxxxzz_1, tg_yyzz_xxxyy_1, tg_yyzz_xxxyyy_0, \
                                         tg_yyzz_xxxyyy_1, tg_yyzz_xxxyyz_0, tg_yyzz_xxxyyz_1, tg_yyzz_xxxyz_1, \
                                         tg_yyzz_xxxyzz_0, tg_yyzz_xxxyzz_1, tg_yyzz_xxxzz_1, tg_yyzz_xxxzzz_0, \
                                         tg_yyzz_xxxzzz_1, tg_yyzz_xxyyy_1, tg_yyzz_xxyyyy_0, tg_yyzz_xxyyyy_1, \
                                         tg_yyzz_xxyyyz_0, tg_yyzz_xxyyyz_1, tg_yyzz_xxyyz_1, tg_yyzz_xxyyzz_0, \
                                         tg_yyzz_xxyyzz_1, tg_yyzz_xxyzz_1, tg_yyzz_xxyzzz_0, tg_yyzz_xxyzzz_1, \
                                         tg_yyzz_xxzzz_1, tg_yyzz_xxzzzz_0, tg_yyzz_xxzzzz_1, tg_yyzz_xyyyy_1, \
                                         tg_yyzz_xyyyyy_0, tg_yyzz_xyyyyy_1, tg_yyzz_xyyyyz_0, tg_yyzz_xyyyyz_1, \
                                         tg_yyzz_xyyyz_1, tg_yyzz_xyyyzz_0, tg_yyzz_xyyyzz_1, tg_yyzz_xyyzz_1, \
                                         tg_yyzz_xyyzzz_0, tg_yyzz_xyyzzz_1, tg_yyzz_xyzzz_1, tg_yyzz_xyzzzz_0, \
                                         tg_yyzz_xyzzzz_1, tg_yyzz_xzzzz_1, tg_yyzz_xzzzzz_0, tg_yyzz_xzzzzz_1, \
                                         tg_yyzz_yyyyy_1, tg_yyzz_yyyyyy_0, tg_yyzz_yyyyyy_1, tg_yyzz_yyyyyz_0, \
                                         tg_yyzz_yyyyyz_1, tg_yyzz_yyyyz_1, tg_yyzz_yyyyzz_0, tg_yyzz_yyyyzz_1, \
                                         tg_yyzz_yyyzz_1, tg_yyzz_yyyzzz_0, tg_yyzz_yyyzzz_1, tg_yyzz_yyzzz_1, \
                                         tg_yyzz_yyzzzz_0, tg_yyzz_yyzzzz_1, tg_yyzz_yzzzz_1, tg_yyzz_yzzzzz_0, \
                                         tg_yyzz_yzzzzz_1, tg_yyzz_zzzzz_1, tg_yyzz_zzzzzz_0, tg_yyzz_zzzzzz_1, \
                                         tg_yzzz_xxxxx_1, tg_yzzz_xxxxxx_0, tg_yzzz_xxxxxx_1, tg_yzzz_xxxxxy_0, \
                                         tg_yzzz_xxxxxy_1, tg_yzzz_xxxxxz_0, tg_yzzz_xxxxxz_1, tg_yzzz_xxxxy_1, \
                                         tg_yzzz_xxxxyy_0, tg_yzzz_xxxxyy_1, tg_yzzz_xxxxyz_0, tg_yzzz_xxxxyz_1, \
                                         tg_yzzz_xxxxz_1, tg_yzzz_xxxxzz_0, tg_yzzz_xxxxzz_1, tg_yzzz_xxxyy_1, \
                                         tg_yzzz_xxxyyy_0, tg_yzzz_xxxyyy_1, tg_yzzz_xxxyyz_0, tg_yzzz_xxxyyz_1, \
                                         tg_yzzz_xxxyz_1, tg_yzzz_xxxyzz_0, tg_yzzz_xxxyzz_1, tg_yzzz_xxxzz_1, \
                                         tg_yzzz_xxxzzz_0, tg_yzzz_xxxzzz_1, tg_yzzz_xxyyy_1, tg_yzzz_xxyyyy_0, \
                                         tg_yzzz_xxyyyy_1, tg_yzzz_xxyyyz_0, tg_yzzz_xxyyyz_1, tg_yzzz_xxyyz_1, \
                                         tg_yzzz_xxyyzz_0, tg_yzzz_xxyyzz_1, tg_yzzz_xxyzz_1, tg_yzzz_xxyzzz_0, \
                                         tg_yzzz_xxyzzz_1, tg_yzzz_xxzzz_1, tg_yzzz_xxzzzz_0, tg_yzzz_xxzzzz_1, \
                                         tg_yzzz_xyyyy_1, tg_yzzz_xyyyyy_0, tg_yzzz_xyyyyy_1, tg_yzzz_xyyyyz_0, \
                                         tg_yzzz_xyyyyz_1, tg_yzzz_xyyyz_1, tg_yzzz_xyyyzz_0, tg_yzzz_xyyyzz_1, \
                                         tg_yzzz_xyyzz_1, tg_yzzz_xyyzzz_0, tg_yzzz_xyyzzz_1, tg_yzzz_xyzzz_1, \
                                         tg_yzzz_xyzzzz_0, tg_yzzz_xyzzzz_1, tg_yzzz_xzzzz_1, tg_yzzz_xzzzzz_0, \
                                         tg_yzzz_xzzzzz_1, tg_yzzz_yyyyy_1, tg_yzzz_yyyyyy_0, tg_yzzz_yyyyyy_1, \
                                         tg_yzzz_yyyyyz_0, tg_yzzz_yyyyyz_1, tg_yzzz_yyyyz_1, tg_yzzz_yyyyzz_0, \
                                         tg_yzzz_yyyyzz_1, tg_yzzz_yyyzz_1, tg_yzzz_yyyzzz_0, tg_yzzz_yyyzzz_1, \
                                         tg_yzzz_yyzzz_1, tg_yzzz_yyzzzz_0, tg_yzzz_yyzzzz_1, tg_yzzz_yzzzz_1, \
                                         tg_yzzz_yzzzzz_0, tg_yzzz_yzzzzz_1, tg_yzzz_zzzzz_1, tg_yzzz_zzzzzz_0, \
                                         tg_yzzz_zzzzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    double fr = wp_x[j]; 

                    tg_xyyyy_xxzzzz_0[j] = pb_x * tg_yyyy_xxzzzz_0[j] + fr * tg_yyyy_xxzzzz_1[j] + fl1_fxn * tg_yyyy_xzzzz_1[j];

                    tg_xyyyy_xyyyyy_0[j] = pb_x * tg_yyyy_xyyyyy_0[j] + fr * tg_yyyy_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyy_1[j];

                    tg_xyyyy_xyyyyz_0[j] = pb_x * tg_yyyy_xyyyyz_0[j] + fr * tg_yyyy_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyyz_1[j];

                    tg_xyyyy_xyyyzz_0[j] = pb_x * tg_yyyy_xyyyzz_0[j] + fr * tg_yyyy_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyyzz_1[j];

                    tg_xyyyy_xyyzzz_0[j] = pb_x * tg_yyyy_xyyzzz_0[j] + fr * tg_yyyy_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyzzz_1[j];

                    tg_xyyyy_xyzzzz_0[j] = pb_x * tg_yyyy_xyzzzz_0[j] + fr * tg_yyyy_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yzzzz_1[j];

                    tg_xyyyy_xzzzzz_0[j] = pb_x * tg_yyyy_xzzzzz_0[j] + fr * tg_yyyy_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zzzzz_1[j];

                    tg_xyyyy_yyyyyy_0[j] = pb_x * tg_yyyy_yyyyyy_0[j] + fr * tg_yyyy_yyyyyy_1[j];

                    tg_xyyyy_yyyyyz_0[j] = pb_x * tg_yyyy_yyyyyz_0[j] + fr * tg_yyyy_yyyyyz_1[j];

                    tg_xyyyy_yyyyzz_0[j] = pb_x * tg_yyyy_yyyyzz_0[j] + fr * tg_yyyy_yyyyzz_1[j];

                    tg_xyyyy_yyyzzz_0[j] = pb_x * tg_yyyy_yyyzzz_0[j] + fr * tg_yyyy_yyyzzz_1[j];

                    tg_xyyyy_yyzzzz_0[j] = pb_x * tg_yyyy_yyzzzz_0[j] + fr * tg_yyyy_yyzzzz_1[j];

                    tg_xyyyy_yzzzzz_0[j] = pb_x * tg_yyyy_yzzzzz_0[j] + fr * tg_yyyy_yzzzzz_1[j];

                    tg_xyyyy_zzzzzz_0[j] = pb_x * tg_yyyy_zzzzzz_0[j] + fr * tg_yyyy_zzzzzz_1[j];

                    tg_xyyyz_xxxxxx_0[j] = pb_x * tg_yyyz_xxxxxx_0[j] + fr * tg_yyyz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyyz_xxxxx_1[j];

                    tg_xyyyz_xxxxxy_0[j] = pb_x * tg_yyyz_xxxxxy_0[j] + fr * tg_yyyz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxy_1[j];

                    tg_xyyyz_xxxxxz_0[j] = pb_x * tg_yyyz_xxxxxz_0[j] + fr * tg_yyyz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyyz_xxxxz_1[j];

                    tg_xyyyz_xxxxyy_0[j] = pb_x * tg_yyyz_xxxxyy_0[j] + fr * tg_yyyz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyy_1[j];

                    tg_xyyyz_xxxxyz_0[j] = pb_x * tg_yyyz_xxxxyz_0[j] + fr * tg_yyyz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxyz_1[j];

                    tg_xyyyz_xxxxzz_0[j] = pb_x * tg_yyyz_xxxxzz_0[j] + fr * tg_yyyz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxxzz_1[j];

                    tg_xyyyz_xxxyyy_0[j] = pb_x * tg_yyyz_xxxyyy_0[j] + fr * tg_yyyz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyy_1[j];

                    tg_xyyyz_xxxyyz_0[j] = pb_x * tg_yyyz_xxxyyz_0[j] + fr * tg_yyyz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyyz_1[j];

                    tg_xyyyz_xxxyzz_0[j] = pb_x * tg_yyyz_xxxyzz_0[j] + fr * tg_yyyz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxyzz_1[j];

                    tg_xyyyz_xxxzzz_0[j] = pb_x * tg_yyyz_xxxzzz_0[j] + fr * tg_yyyz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxzzz_1[j];

                    tg_xyyyz_xxyyyy_0[j] = pb_x * tg_yyyz_xxyyyy_0[j] + fr * tg_yyyz_xxyyyy_1[j] + fl1_fxn * tg_yyyz_xyyyy_1[j];

                    tg_xyyyz_xxyyyz_0[j] = pb_x * tg_yyyz_xxyyyz_0[j] + fr * tg_yyyz_xxyyyz_1[j] + fl1_fxn * tg_yyyz_xyyyz_1[j];

                    tg_xyyyz_xxyyzz_0[j] = pb_x * tg_yyyz_xxyyzz_0[j] + fr * tg_yyyz_xxyyzz_1[j] + fl1_fxn * tg_yyyz_xyyzz_1[j];

                    tg_xyyyz_xxyzzz_0[j] = pb_x * tg_yyyz_xxyzzz_0[j] + fr * tg_yyyz_xxyzzz_1[j] + fl1_fxn * tg_yyyz_xyzzz_1[j];

                    tg_xyyyz_xxzzzz_0[j] = pb_x * tg_yyyz_xxzzzz_0[j] + fr * tg_yyyz_xxzzzz_1[j] + fl1_fxn * tg_yyyz_xzzzz_1[j];

                    tg_xyyyz_xyyyyy_0[j] = pb_x * tg_yyyz_xyyyyy_0[j] + fr * tg_yyyz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyy_1[j];

                    tg_xyyyz_xyyyyz_0[j] = pb_x * tg_yyyz_xyyyyz_0[j] + fr * tg_yyyz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyyz_1[j];

                    tg_xyyyz_xyyyzz_0[j] = pb_x * tg_yyyz_xyyyzz_0[j] + fr * tg_yyyz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyyzz_1[j];

                    tg_xyyyz_xyyzzz_0[j] = pb_x * tg_yyyz_xyyzzz_0[j] + fr * tg_yyyz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyzzz_1[j];

                    tg_xyyyz_xyzzzz_0[j] = pb_x * tg_yyyz_xyzzzz_0[j] + fr * tg_yyyz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yzzzz_1[j];

                    tg_xyyyz_xzzzzz_0[j] = pb_x * tg_yyyz_xzzzzz_0[j] + fr * tg_yyyz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zzzzz_1[j];

                    tg_xyyyz_yyyyyy_0[j] = pb_x * tg_yyyz_yyyyyy_0[j] + fr * tg_yyyz_yyyyyy_1[j];

                    tg_xyyyz_yyyyyz_0[j] = pb_x * tg_yyyz_yyyyyz_0[j] + fr * tg_yyyz_yyyyyz_1[j];

                    tg_xyyyz_yyyyzz_0[j] = pb_x * tg_yyyz_yyyyzz_0[j] + fr * tg_yyyz_yyyyzz_1[j];

                    tg_xyyyz_yyyzzz_0[j] = pb_x * tg_yyyz_yyyzzz_0[j] + fr * tg_yyyz_yyyzzz_1[j];

                    tg_xyyyz_yyzzzz_0[j] = pb_x * tg_yyyz_yyzzzz_0[j] + fr * tg_yyyz_yyzzzz_1[j];

                    tg_xyyyz_yzzzzz_0[j] = pb_x * tg_yyyz_yzzzzz_0[j] + fr * tg_yyyz_yzzzzz_1[j];

                    tg_xyyyz_zzzzzz_0[j] = pb_x * tg_yyyz_zzzzzz_0[j] + fr * tg_yyyz_zzzzzz_1[j];

                    tg_xyyzz_xxxxxx_0[j] = pb_x * tg_yyzz_xxxxxx_0[j] + fr * tg_yyzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yyzz_xxxxx_1[j];

                    tg_xyyzz_xxxxxy_0[j] = pb_x * tg_yyzz_xxxxxy_0[j] + fr * tg_yyzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxy_1[j];

                    tg_xyyzz_xxxxxz_0[j] = pb_x * tg_yyzz_xxxxxz_0[j] + fr * tg_yyzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yyzz_xxxxz_1[j];

                    tg_xyyzz_xxxxyy_0[j] = pb_x * tg_yyzz_xxxxyy_0[j] + fr * tg_yyzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyy_1[j];

                    tg_xyyzz_xxxxyz_0[j] = pb_x * tg_yyzz_xxxxyz_0[j] + fr * tg_yyzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxyz_1[j];

                    tg_xyyzz_xxxxzz_0[j] = pb_x * tg_yyzz_xxxxzz_0[j] + fr * tg_yyzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxxzz_1[j];

                    tg_xyyzz_xxxyyy_0[j] = pb_x * tg_yyzz_xxxyyy_0[j] + fr * tg_yyzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyy_1[j];

                    tg_xyyzz_xxxyyz_0[j] = pb_x * tg_yyzz_xxxyyz_0[j] + fr * tg_yyzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyyz_1[j];

                    tg_xyyzz_xxxyzz_0[j] = pb_x * tg_yyzz_xxxyzz_0[j] + fr * tg_yyzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxyzz_1[j];

                    tg_xyyzz_xxxzzz_0[j] = pb_x * tg_yyzz_xxxzzz_0[j] + fr * tg_yyzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxzzz_1[j];

                    tg_xyyzz_xxyyyy_0[j] = pb_x * tg_yyzz_xxyyyy_0[j] + fr * tg_yyzz_xxyyyy_1[j] + fl1_fxn * tg_yyzz_xyyyy_1[j];

                    tg_xyyzz_xxyyyz_0[j] = pb_x * tg_yyzz_xxyyyz_0[j] + fr * tg_yyzz_xxyyyz_1[j] + fl1_fxn * tg_yyzz_xyyyz_1[j];

                    tg_xyyzz_xxyyzz_0[j] = pb_x * tg_yyzz_xxyyzz_0[j] + fr * tg_yyzz_xxyyzz_1[j] + fl1_fxn * tg_yyzz_xyyzz_1[j];

                    tg_xyyzz_xxyzzz_0[j] = pb_x * tg_yyzz_xxyzzz_0[j] + fr * tg_yyzz_xxyzzz_1[j] + fl1_fxn * tg_yyzz_xyzzz_1[j];

                    tg_xyyzz_xxzzzz_0[j] = pb_x * tg_yyzz_xxzzzz_0[j] + fr * tg_yyzz_xxzzzz_1[j] + fl1_fxn * tg_yyzz_xzzzz_1[j];

                    tg_xyyzz_xyyyyy_0[j] = pb_x * tg_yyzz_xyyyyy_0[j] + fr * tg_yyzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyy_1[j];

                    tg_xyyzz_xyyyyz_0[j] = pb_x * tg_yyzz_xyyyyz_0[j] + fr * tg_yyzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyyz_1[j];

                    tg_xyyzz_xyyyzz_0[j] = pb_x * tg_yyzz_xyyyzz_0[j] + fr * tg_yyzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyyzz_1[j];

                    tg_xyyzz_xyyzzz_0[j] = pb_x * tg_yyzz_xyyzzz_0[j] + fr * tg_yyzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyzzz_1[j];

                    tg_xyyzz_xyzzzz_0[j] = pb_x * tg_yyzz_xyzzzz_0[j] + fr * tg_yyzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yzzzz_1[j];

                    tg_xyyzz_xzzzzz_0[j] = pb_x * tg_yyzz_xzzzzz_0[j] + fr * tg_yyzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zzzzz_1[j];

                    tg_xyyzz_yyyyyy_0[j] = pb_x * tg_yyzz_yyyyyy_0[j] + fr * tg_yyzz_yyyyyy_1[j];

                    tg_xyyzz_yyyyyz_0[j] = pb_x * tg_yyzz_yyyyyz_0[j] + fr * tg_yyzz_yyyyyz_1[j];

                    tg_xyyzz_yyyyzz_0[j] = pb_x * tg_yyzz_yyyyzz_0[j] + fr * tg_yyzz_yyyyzz_1[j];

                    tg_xyyzz_yyyzzz_0[j] = pb_x * tg_yyzz_yyyzzz_0[j] + fr * tg_yyzz_yyyzzz_1[j];

                    tg_xyyzz_yyzzzz_0[j] = pb_x * tg_yyzz_yyzzzz_0[j] + fr * tg_yyzz_yyzzzz_1[j];

                    tg_xyyzz_yzzzzz_0[j] = pb_x * tg_yyzz_yzzzzz_0[j] + fr * tg_yyzz_yzzzzz_1[j];

                    tg_xyyzz_zzzzzz_0[j] = pb_x * tg_yyzz_zzzzzz_0[j] + fr * tg_yyzz_zzzzzz_1[j];

                    tg_xyzzz_xxxxxx_0[j] = pb_x * tg_yzzz_xxxxxx_0[j] + fr * tg_yzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_yzzz_xxxxx_1[j];

                    tg_xyzzz_xxxxxy_0[j] = pb_x * tg_yzzz_xxxxxy_0[j] + fr * tg_yzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxy_1[j];

                    tg_xyzzz_xxxxxz_0[j] = pb_x * tg_yzzz_xxxxxz_0[j] + fr * tg_yzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_yzzz_xxxxz_1[j];

                    tg_xyzzz_xxxxyy_0[j] = pb_x * tg_yzzz_xxxxyy_0[j] + fr * tg_yzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyy_1[j];

                    tg_xyzzz_xxxxyz_0[j] = pb_x * tg_yzzz_xxxxyz_0[j] + fr * tg_yzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxyz_1[j];

                    tg_xyzzz_xxxxzz_0[j] = pb_x * tg_yzzz_xxxxzz_0[j] + fr * tg_yzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxxzz_1[j];

                    tg_xyzzz_xxxyyy_0[j] = pb_x * tg_yzzz_xxxyyy_0[j] + fr * tg_yzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyy_1[j];

                    tg_xyzzz_xxxyyz_0[j] = pb_x * tg_yzzz_xxxyyz_0[j] + fr * tg_yzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyyz_1[j];

                    tg_xyzzz_xxxyzz_0[j] = pb_x * tg_yzzz_xxxyzz_0[j] + fr * tg_yzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxyzz_1[j];

                    tg_xyzzz_xxxzzz_0[j] = pb_x * tg_yzzz_xxxzzz_0[j] + fr * tg_yzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxzzz_1[j];

                    tg_xyzzz_xxyyyy_0[j] = pb_x * tg_yzzz_xxyyyy_0[j] + fr * tg_yzzz_xxyyyy_1[j] + fl1_fxn * tg_yzzz_xyyyy_1[j];

                    tg_xyzzz_xxyyyz_0[j] = pb_x * tg_yzzz_xxyyyz_0[j] + fr * tg_yzzz_xxyyyz_1[j] + fl1_fxn * tg_yzzz_xyyyz_1[j];

                    tg_xyzzz_xxyyzz_0[j] = pb_x * tg_yzzz_xxyyzz_0[j] + fr * tg_yzzz_xxyyzz_1[j] + fl1_fxn * tg_yzzz_xyyzz_1[j];

                    tg_xyzzz_xxyzzz_0[j] = pb_x * tg_yzzz_xxyzzz_0[j] + fr * tg_yzzz_xxyzzz_1[j] + fl1_fxn * tg_yzzz_xyzzz_1[j];

                    tg_xyzzz_xxzzzz_0[j] = pb_x * tg_yzzz_xxzzzz_0[j] + fr * tg_yzzz_xxzzzz_1[j] + fl1_fxn * tg_yzzz_xzzzz_1[j];

                    tg_xyzzz_xyyyyy_0[j] = pb_x * tg_yzzz_xyyyyy_0[j] + fr * tg_yzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyy_1[j];

                    tg_xyzzz_xyyyyz_0[j] = pb_x * tg_yzzz_xyyyyz_0[j] + fr * tg_yzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyyz_1[j];

                    tg_xyzzz_xyyyzz_0[j] = pb_x * tg_yzzz_xyyyzz_0[j] + fr * tg_yzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyyzz_1[j];

                    tg_xyzzz_xyyzzz_0[j] = pb_x * tg_yzzz_xyyzzz_0[j] + fr * tg_yzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyzzz_1[j];

                    tg_xyzzz_xyzzzz_0[j] = pb_x * tg_yzzz_xyzzzz_0[j] + fr * tg_yzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yzzzz_1[j];

                    tg_xyzzz_xzzzzz_0[j] = pb_x * tg_yzzz_xzzzzz_0[j] + fr * tg_yzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zzzzz_1[j];

                    tg_xyzzz_yyyyyy_0[j] = pb_x * tg_yzzz_yyyyyy_0[j] + fr * tg_yzzz_yyyyyy_1[j];

                    tg_xyzzz_yyyyyz_0[j] = pb_x * tg_yzzz_yyyyyz_0[j] + fr * tg_yzzz_yyyyyz_1[j];

                    tg_xyzzz_yyyyzz_0[j] = pb_x * tg_yzzz_yyyyzz_0[j] + fr * tg_yzzz_yyyyzz_1[j];

                    tg_xyzzz_yyyzzz_0[j] = pb_x * tg_yzzz_yyyzzz_0[j] + fr * tg_yzzz_yyyzzz_1[j];

                    tg_xyzzz_yyzzzz_0[j] = pb_x * tg_yzzz_yyzzzz_0[j] + fr * tg_yzzz_yyzzzz_1[j];

                    tg_xyzzz_yzzzzz_0[j] = pb_x * tg_yzzz_yzzzzz_0[j] + fr * tg_yzzz_yzzzzz_1[j];

                    tg_xyzzz_zzzzzz_0[j] = pb_x * tg_yzzz_zzzzzz_0[j] + fr * tg_yzzz_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSI_392_490(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_4_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {5, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
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

                // set up pointers to integrals

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

                // Batch of Integrals (392,490)

                #pragma omp simd aligned(fxn, fza, tg_xzzzz_xxxxxx_0, tg_xzzzz_xxxxxy_0, tg_xzzzz_xxxxxz_0, \
                                         tg_xzzzz_xxxxyy_0, tg_xzzzz_xxxxyz_0, tg_xzzzz_xxxxzz_0, tg_xzzzz_xxxyyy_0, \
                                         tg_xzzzz_xxxyyz_0, tg_xzzzz_xxxyzz_0, tg_xzzzz_xxxzzz_0, tg_xzzzz_xxyyyy_0, \
                                         tg_xzzzz_xxyyyz_0, tg_xzzzz_xxyyzz_0, tg_xzzzz_xxyzzz_0, tg_xzzzz_xxzzzz_0, \
                                         tg_xzzzz_xyyyyy_0, tg_xzzzz_xyyyyz_0, tg_xzzzz_xyyyzz_0, tg_xzzzz_xyyzzz_0, \
                                         tg_xzzzz_xyzzzz_0, tg_xzzzz_xzzzzz_0, tg_xzzzz_yyyyyy_0, tg_xzzzz_yyyyyz_0, \
                                         tg_xzzzz_yyyyzz_0, tg_xzzzz_yyyzzz_0, tg_xzzzz_yyzzzz_0, tg_xzzzz_yzzzzz_0, \
                                         tg_xzzzz_zzzzzz_0, tg_yyy_xxxxxx_0, tg_yyy_xxxxxx_1, tg_yyy_xxxxxy_0, tg_yyy_xxxxxy_1, \
                                         tg_yyy_xxxxxz_0, tg_yyy_xxxxxz_1, tg_yyy_xxxxyy_0, tg_yyy_xxxxyy_1, tg_yyy_xxxxyz_0, \
                                         tg_yyy_xxxxyz_1, tg_yyy_xxxxzz_0, tg_yyy_xxxxzz_1, tg_yyy_xxxyyy_0, tg_yyy_xxxyyy_1, \
                                         tg_yyy_xxxyyz_0, tg_yyy_xxxyyz_1, tg_yyy_xxxyzz_0, tg_yyy_xxxyzz_1, tg_yyy_xxxzzz_0, \
                                         tg_yyy_xxxzzz_1, tg_yyy_xxyyyy_0, tg_yyy_xxyyyy_1, tg_yyy_xxyyyz_0, tg_yyy_xxyyyz_1, \
                                         tg_yyy_xxyyzz_0, tg_yyy_xxyyzz_1, tg_yyy_xxyzzz_0, tg_yyy_xxyzzz_1, tg_yyy_xxzzzz_0, \
                                         tg_yyy_xxzzzz_1, tg_yyy_xyyyyy_0, tg_yyy_xyyyyy_1, tg_yyy_xyyyyz_0, tg_yyy_xyyyyz_1, \
                                         tg_yyy_xyyyzz_0, tg_yyy_xyyyzz_1, tg_yyy_xyyzzz_0, tg_yyy_xyyzzz_1, tg_yyy_xyzzzz_0, \
                                         tg_yyy_xyzzzz_1, tg_yyy_xzzzzz_0, tg_yyy_xzzzzz_1, tg_yyy_yyyyyy_0, tg_yyy_yyyyyy_1, \
                                         tg_yyy_yyyyyz_0, tg_yyy_yyyyyz_1, tg_yyy_yyyyzz_0, tg_yyy_yyyyzz_1, tg_yyy_yyyzzz_0, \
                                         tg_yyy_yyyzzz_1, tg_yyy_yyzzzz_0, tg_yyy_yyzzzz_1, tg_yyy_yzzzzz_0, tg_yyy_yzzzzz_1, \
                                         tg_yyy_zzzzzz_0, tg_yyy_zzzzzz_1, tg_yyyy_xxxxx_1, tg_yyyy_xxxxxx_0, \
                                         tg_yyyy_xxxxxx_1, tg_yyyy_xxxxxy_0, tg_yyyy_xxxxxy_1, tg_yyyy_xxxxxz_0, \
                                         tg_yyyy_xxxxxz_1, tg_yyyy_xxxxy_1, tg_yyyy_xxxxyy_0, tg_yyyy_xxxxyy_1, \
                                         tg_yyyy_xxxxyz_0, tg_yyyy_xxxxyz_1, tg_yyyy_xxxxz_1, tg_yyyy_xxxxzz_0, \
                                         tg_yyyy_xxxxzz_1, tg_yyyy_xxxyy_1, tg_yyyy_xxxyyy_0, tg_yyyy_xxxyyy_1, \
                                         tg_yyyy_xxxyyz_0, tg_yyyy_xxxyyz_1, tg_yyyy_xxxyz_1, tg_yyyy_xxxyzz_0, \
                                         tg_yyyy_xxxyzz_1, tg_yyyy_xxxzz_1, tg_yyyy_xxxzzz_0, tg_yyyy_xxxzzz_1, \
                                         tg_yyyy_xxyyy_1, tg_yyyy_xxyyyy_0, tg_yyyy_xxyyyy_1, tg_yyyy_xxyyyz_0, \
                                         tg_yyyy_xxyyyz_1, tg_yyyy_xxyyz_1, tg_yyyy_xxyyzz_0, tg_yyyy_xxyyzz_1, \
                                         tg_yyyy_xxyzz_1, tg_yyyy_xxyzzz_0, tg_yyyy_xxyzzz_1, tg_yyyy_xxzzz_1, \
                                         tg_yyyy_xxzzzz_0, tg_yyyy_xxzzzz_1, tg_yyyy_xyyyy_1, tg_yyyy_xyyyyy_0, \
                                         tg_yyyy_xyyyyy_1, tg_yyyy_xyyyyz_0, tg_yyyy_xyyyyz_1, tg_yyyy_xyyyz_1, \
                                         tg_yyyy_xyyyzz_0, tg_yyyy_xyyyzz_1, tg_yyyy_xyyzz_1, tg_yyyy_xyyzzz_0, \
                                         tg_yyyy_xyyzzz_1, tg_yyyy_xyzzz_1, tg_yyyy_xyzzzz_0, tg_yyyy_xyzzzz_1, \
                                         tg_yyyy_xzzzz_1, tg_yyyy_xzzzzz_0, tg_yyyy_xzzzzz_1, tg_yyyy_yyyyy_1, \
                                         tg_yyyy_yyyyyy_0, tg_yyyy_yyyyyy_1, tg_yyyy_yyyyyz_0, tg_yyyy_yyyyyz_1, \
                                         tg_yyyy_yyyyz_1, tg_yyyy_yyyyzz_0, tg_yyyy_yyyyzz_1, tg_yyyy_yyyzz_1, \
                                         tg_yyyy_yyyzzz_0, tg_yyyy_yyyzzz_1, tg_yyyy_yyzzz_1, tg_yyyy_yyzzzz_0, \
                                         tg_yyyy_yyzzzz_1, tg_yyyy_yzzzz_1, tg_yyyy_yzzzzz_0, tg_yyyy_yzzzzz_1, \
                                         tg_yyyy_zzzzz_1, tg_yyyy_zzzzzz_0, tg_yyyy_zzzzzz_1, tg_yyyyy_xxxxxx_0, \
                                         tg_yyyyy_xxxxxy_0, tg_yyyyy_xxxxxz_0, tg_yyyyy_xxxxyy_0, tg_yyyyy_xxxxyz_0, \
                                         tg_yyyyy_xxxxzz_0, tg_yyyyy_xxxyyy_0, tg_yyyyy_xxxyyz_0, tg_yyyyy_xxxyzz_0, \
                                         tg_yyyyy_xxxzzz_0, tg_yyyyy_xxyyyy_0, tg_yyyyy_xxyyyz_0, tg_yyyyy_xxyyzz_0, \
                                         tg_yyyyy_xxyzzz_0, tg_yyyyy_xxzzzz_0, tg_yyyyy_xyyyyy_0, tg_yyyyy_xyyyyz_0, \
                                         tg_yyyyy_xyyyzz_0, tg_yyyyy_xyyzzz_0, tg_yyyyy_xyzzzz_0, tg_yyyyy_xzzzzz_0, \
                                         tg_yyyyy_yyyyyy_0, tg_yyyyy_yyyyyz_0, tg_yyyyy_yyyyzz_0, tg_yyyyy_yyyzzz_0, \
                                         tg_yyyyy_yyzzzz_0, tg_yyyyy_yzzzzz_0, tg_yyyyy_zzzzzz_0, tg_yyyyz_xxxxxx_0, \
                                         tg_yyyyz_xxxxxy_0, tg_yyyyz_xxxxxz_0, tg_yyyyz_xxxxyy_0, tg_yyyyz_xxxxyz_0, \
                                         tg_yyyyz_xxxxzz_0, tg_yyyyz_xxxyyy_0, tg_yyyyz_xxxyyz_0, tg_yyyyz_xxxyzz_0, \
                                         tg_yyyyz_xxxzzz_0, tg_yyyyz_xxyyyy_0, tg_yyyyz_xxyyyz_0, tg_yyyyz_xxyyzz_0, \
                                         tg_yyyyz_xxyzzz_0, tg_yyyyz_xxzzzz_0, tg_yyyyz_xyyyyy_0, tg_yyyyz_xyyyyz_0, \
                                         tg_yyyyz_xyyyzz_0, tg_yyyyz_xyyzzz_0, tg_yyyyz_xyzzzz_0, tg_yyyyz_xzzzzz_0, \
                                         tg_yyyyz_yyyyyy_0, tg_yyyyz_yyyyyz_0, tg_yyyyz_yyyyzz_0, tg_yyyyz_yyyzzz_0, \
                                         tg_yyyyz_yyzzzz_0, tg_yyyyz_yzzzzz_0, tg_yyyyz_zzzzzz_0, tg_yyyz_xxxxx_1, \
                                         tg_yyyz_xxxxxx_0, tg_yyyz_xxxxxx_1, tg_yyyz_xxxxxy_0, tg_yyyz_xxxxxy_1, \
                                         tg_yyyz_xxxxxz_0, tg_yyyz_xxxxxz_1, tg_yyyz_xxxxy_1, tg_yyyz_xxxxyy_0, \
                                         tg_yyyz_xxxxyy_1, tg_yyyz_xxxxyz_0, tg_yyyz_xxxxyz_1, tg_yyyz_xxxxz_1, \
                                         tg_yyyz_xxxxzz_0, tg_yyyz_xxxxzz_1, tg_yyyz_xxxyy_1, tg_yyyz_xxxyyy_0, \
                                         tg_yyyz_xxxyyy_1, tg_yyyz_xxxyyz_0, tg_yyyz_xxxyyz_1, tg_yyyz_xxxyz_1, \
                                         tg_yyyz_xxxyzz_0, tg_yyyz_xxxyzz_1, tg_yyyz_xxxzz_1, tg_yyyz_xxxzzz_0, \
                                         tg_yyyz_xxxzzz_1, tg_yyyz_xxyyy_1, tg_yyyz_xxyyyy_0, tg_yyyz_xxyyyy_1, \
                                         tg_yyyz_xxyyyz_0, tg_yyyz_xxyyyz_1, tg_yyyz_xxyyz_1, tg_yyyz_xxyyzz_0, \
                                         tg_yyyz_xxyyzz_1, tg_yyyz_xxyzz_1, tg_yyyz_xxyzzz_0, tg_yyyz_xxyzzz_1, \
                                         tg_yyyz_xxzzz_1, tg_yyyz_xxzzzz_0, tg_yyyz_xxzzzz_1, tg_yyyz_xyyyy_1, \
                                         tg_yyyz_xyyyyy_0, tg_yyyz_xyyyyy_1, tg_yyyz_xyyyyz_0, tg_yyyz_xyyyyz_1, \
                                         tg_yyyz_xyyyz_1, tg_yyyz_xyyyzz_0, tg_yyyz_xyyyzz_1, tg_yyyz_xyyzz_1, \
                                         tg_yyyz_xyyzzz_0, tg_yyyz_xyyzzz_1, tg_yyyz_xyzzz_1, tg_yyyz_xyzzzz_0, \
                                         tg_yyyz_xyzzzz_1, tg_yyyz_xzzzz_1, tg_yyyz_xzzzzz_0, tg_yyyz_xzzzzz_1, \
                                         tg_yyyz_yyyyy_1, tg_yyyz_yyyyyy_0, tg_yyyz_yyyyyy_1, tg_yyyz_yyyyyz_0, \
                                         tg_yyyz_yyyyyz_1, tg_yyyz_yyyyz_1, tg_yyyz_yyyyzz_0, tg_yyyz_yyyyzz_1, \
                                         tg_yyyz_yyyzz_1, tg_yyyz_yyyzzz_0, tg_yyyz_yyyzzz_1, tg_yyyz_yyzzz_1, \
                                         tg_yyyz_yyzzzz_0, tg_yyyz_yyzzzz_1, tg_yyyz_yzzzz_1, tg_yyyz_yzzzzz_0, \
                                         tg_yyyz_yzzzzz_1, tg_yyyz_zzzzz_1, tg_yyyz_zzzzzz_0, tg_yyyz_zzzzzz_1, \
                                         tg_yyyzz_xxxxxx_0, tg_yyyzz_xxxxxy_0, tg_yyyzz_xxxxxz_0, tg_yyyzz_xxxxyy_0, \
                                         tg_yyyzz_xxxxyz_0, tg_yyyzz_xxxxzz_0, tg_yyyzz_xxxyyy_0, tg_yyyzz_xxxyyz_0, \
                                         tg_yyyzz_xxxyzz_0, tg_yyyzz_xxxzzz_0, tg_yyyzz_xxyyyy_0, tg_yyyzz_xxyyyz_0, \
                                         tg_yyyzz_xxyyzz_0, tg_yyyzz_xxyzzz_0, tg_yyz_xxxxxx_0, tg_yyz_xxxxxx_1, \
                                         tg_yyz_xxxxxy_0, tg_yyz_xxxxxy_1, tg_yyz_xxxxxz_0, tg_yyz_xxxxxz_1, tg_yyz_xxxxyy_0, \
                                         tg_yyz_xxxxyy_1, tg_yyz_xxxxyz_0, tg_yyz_xxxxyz_1, tg_yyz_xxxxzz_0, tg_yyz_xxxxzz_1, \
                                         tg_yyz_xxxyyy_0, tg_yyz_xxxyyy_1, tg_yyz_xxxyyz_0, tg_yyz_xxxyyz_1, tg_yyz_xxxyzz_0, \
                                         tg_yyz_xxxyzz_1, tg_yyz_xxxzzz_0, tg_yyz_xxxzzz_1, tg_yyz_xxyyyy_0, tg_yyz_xxyyyy_1, \
                                         tg_yyz_xxyyyz_0, tg_yyz_xxyyyz_1, tg_yyz_xxyyzz_0, tg_yyz_xxyyzz_1, tg_yyz_xxyzzz_0, \
                                         tg_yyz_xxyzzz_1, tg_yyz_xxzzzz_0, tg_yyz_xxzzzz_1, tg_yyz_xyyyyy_0, tg_yyz_xyyyyy_1, \
                                         tg_yyz_xyyyyz_0, tg_yyz_xyyyyz_1, tg_yyz_xyyyzz_0, tg_yyz_xyyyzz_1, tg_yyz_xyyzzz_0, \
                                         tg_yyz_xyyzzz_1, tg_yyz_xyzzzz_0, tg_yyz_xyzzzz_1, tg_yyz_xzzzzz_0, tg_yyz_xzzzzz_1, \
                                         tg_yyz_yyyyyy_0, tg_yyz_yyyyyy_1, tg_yyz_yyyyyz_0, tg_yyz_yyyyyz_1, tg_yyz_yyyyzz_0, \
                                         tg_yyz_yyyyzz_1, tg_yyz_yyyzzz_0, tg_yyz_yyyzzz_1, tg_yyz_yyzzzz_0, tg_yyz_yyzzzz_1, \
                                         tg_yyz_yzzzzz_0, tg_yyz_yzzzzz_1, tg_yyz_zzzzzz_0, tg_yyz_zzzzzz_1, tg_yyzz_xxxxx_1, \
                                         tg_yyzz_xxxxxx_0, tg_yyzz_xxxxxx_1, tg_yyzz_xxxxxy_0, tg_yyzz_xxxxxy_1, \
                                         tg_yyzz_xxxxxz_0, tg_yyzz_xxxxxz_1, tg_yyzz_xxxxy_1, tg_yyzz_xxxxyy_0, \
                                         tg_yyzz_xxxxyy_1, tg_yyzz_xxxxyz_0, tg_yyzz_xxxxyz_1, tg_yyzz_xxxxz_1, \
                                         tg_yyzz_xxxxzz_0, tg_yyzz_xxxxzz_1, tg_yyzz_xxxyy_1, tg_yyzz_xxxyyy_0, \
                                         tg_yyzz_xxxyyy_1, tg_yyzz_xxxyyz_0, tg_yyzz_xxxyyz_1, tg_yyzz_xxxyz_1, \
                                         tg_yyzz_xxxyzz_0, tg_yyzz_xxxyzz_1, tg_yyzz_xxxzz_1, tg_yyzz_xxxzzz_0, \
                                         tg_yyzz_xxxzzz_1, tg_yyzz_xxyyy_1, tg_yyzz_xxyyyy_0, tg_yyzz_xxyyyy_1, \
                                         tg_yyzz_xxyyyz_0, tg_yyzz_xxyyyz_1, tg_yyzz_xxyyz_1, tg_yyzz_xxyyzz_0, \
                                         tg_yyzz_xxyyzz_1, tg_yyzz_xxyzz_1, tg_yyzz_xxyzzz_0, tg_yyzz_xxyzzz_1, \
                                         tg_yyzz_xxzzz_1, tg_yzz_xxxxxx_0, tg_yzz_xxxxxx_1, tg_yzz_xxxxxy_0, tg_yzz_xxxxxy_1, \
                                         tg_yzz_xxxxxz_0, tg_yzz_xxxxxz_1, tg_yzz_xxxxyy_0, tg_yzz_xxxxyy_1, tg_yzz_xxxxyz_0, \
                                         tg_yzz_xxxxyz_1, tg_yzz_xxxxzz_0, tg_yzz_xxxxzz_1, tg_yzz_xxxyyy_0, tg_yzz_xxxyyy_1, \
                                         tg_yzz_xxxyyz_0, tg_yzz_xxxyyz_1, tg_yzz_xxxyzz_0, tg_yzz_xxxyzz_1, tg_yzz_xxxzzz_0, \
                                         tg_yzz_xxxzzz_1, tg_yzz_xxyyyy_0, tg_yzz_xxyyyy_1, tg_yzz_xxyyyz_0, tg_yzz_xxyyyz_1, \
                                         tg_yzz_xxyyzz_0, tg_yzz_xxyyzz_1, tg_yzz_xxyzzz_0, tg_yzz_xxyzzz_1, tg_zzzz_xxxxx_1, \
                                         tg_zzzz_xxxxxx_0, tg_zzzz_xxxxxx_1, tg_zzzz_xxxxxy_0, tg_zzzz_xxxxxy_1, \
                                         tg_zzzz_xxxxxz_0, tg_zzzz_xxxxxz_1, tg_zzzz_xxxxy_1, tg_zzzz_xxxxyy_0, \
                                         tg_zzzz_xxxxyy_1, tg_zzzz_xxxxyz_0, tg_zzzz_xxxxyz_1, tg_zzzz_xxxxz_1, \
                                         tg_zzzz_xxxxzz_0, tg_zzzz_xxxxzz_1, tg_zzzz_xxxyy_1, tg_zzzz_xxxyyy_0, \
                                         tg_zzzz_xxxyyy_1, tg_zzzz_xxxyyz_0, tg_zzzz_xxxyyz_1, tg_zzzz_xxxyz_1, \
                                         tg_zzzz_xxxyzz_0, tg_zzzz_xxxyzz_1, tg_zzzz_xxxzz_1, tg_zzzz_xxxzzz_0, \
                                         tg_zzzz_xxxzzz_1, tg_zzzz_xxyyy_1, tg_zzzz_xxyyyy_0, tg_zzzz_xxyyyy_1, \
                                         tg_zzzz_xxyyyz_0, tg_zzzz_xxyyyz_1, tg_zzzz_xxyyz_1, tg_zzzz_xxyyzz_0, \
                                         tg_zzzz_xxyyzz_1, tg_zzzz_xxyzz_1, tg_zzzz_xxyzzz_0, tg_zzzz_xxyzzz_1, \
                                         tg_zzzz_xxzzz_1, tg_zzzz_xxzzzz_0, tg_zzzz_xxzzzz_1, tg_zzzz_xyyyy_1, \
                                         tg_zzzz_xyyyyy_0, tg_zzzz_xyyyyy_1, tg_zzzz_xyyyyz_0, tg_zzzz_xyyyyz_1, \
                                         tg_zzzz_xyyyz_1, tg_zzzz_xyyyzz_0, tg_zzzz_xyyyzz_1, tg_zzzz_xyyzz_1, \
                                         tg_zzzz_xyyzzz_0, tg_zzzz_xyyzzz_1, tg_zzzz_xyzzz_1, tg_zzzz_xyzzzz_0, \
                                         tg_zzzz_xyzzzz_1, tg_zzzz_xzzzz_1, tg_zzzz_xzzzzz_0, tg_zzzz_xzzzzz_1, \
                                         tg_zzzz_yyyyy_1, tg_zzzz_yyyyyy_0, tg_zzzz_yyyyyy_1, tg_zzzz_yyyyyz_0, \
                                         tg_zzzz_yyyyyz_1, tg_zzzz_yyyyz_1, tg_zzzz_yyyyzz_0, tg_zzzz_yyyyzz_1, \
                                         tg_zzzz_yyyzz_1, tg_zzzz_yyyzzz_0, tg_zzzz_yyyzzz_1, tg_zzzz_yyzzz_1, \
                                         tg_zzzz_yyzzzz_0, tg_zzzz_yyzzzz_1, tg_zzzz_yzzzz_1, tg_zzzz_yzzzzz_0, \
                                         tg_zzzz_yzzzzz_1, tg_zzzz_zzzzz_1, tg_zzzz_zzzzzz_0, tg_zzzz_zzzzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xzzzz_xxxxxx_0[j] = pb_x * tg_zzzz_xxxxxx_0[j] + fr * tg_zzzz_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_zzzz_xxxxx_1[j];

                    tg_xzzzz_xxxxxy_0[j] = pb_x * tg_zzzz_xxxxxy_0[j] + fr * tg_zzzz_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxy_1[j];

                    tg_xzzzz_xxxxxz_0[j] = pb_x * tg_zzzz_xxxxxz_0[j] + fr * tg_zzzz_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_zzzz_xxxxz_1[j];

                    tg_xzzzz_xxxxyy_0[j] = pb_x * tg_zzzz_xxxxyy_0[j] + fr * tg_zzzz_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyy_1[j];

                    tg_xzzzz_xxxxyz_0[j] = pb_x * tg_zzzz_xxxxyz_0[j] + fr * tg_zzzz_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxyz_1[j];

                    tg_xzzzz_xxxxzz_0[j] = pb_x * tg_zzzz_xxxxzz_0[j] + fr * tg_zzzz_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxxzz_1[j];

                    tg_xzzzz_xxxyyy_0[j] = pb_x * tg_zzzz_xxxyyy_0[j] + fr * tg_zzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyy_1[j];

                    tg_xzzzz_xxxyyz_0[j] = pb_x * tg_zzzz_xxxyyz_0[j] + fr * tg_zzzz_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyz_1[j];

                    tg_xzzzz_xxxyzz_0[j] = pb_x * tg_zzzz_xxxyzz_0[j] + fr * tg_zzzz_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyzz_1[j];

                    tg_xzzzz_xxxzzz_0[j] = pb_x * tg_zzzz_xxxzzz_0[j] + fr * tg_zzzz_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxzzz_1[j];

                    tg_xzzzz_xxyyyy_0[j] = pb_x * tg_zzzz_xxyyyy_0[j] + fr * tg_zzzz_xxyyyy_1[j] + fl1_fxn * tg_zzzz_xyyyy_1[j];

                    tg_xzzzz_xxyyyz_0[j] = pb_x * tg_zzzz_xxyyyz_0[j] + fr * tg_zzzz_xxyyyz_1[j] + fl1_fxn * tg_zzzz_xyyyz_1[j];

                    tg_xzzzz_xxyyzz_0[j] = pb_x * tg_zzzz_xxyyzz_0[j] + fr * tg_zzzz_xxyyzz_1[j] + fl1_fxn * tg_zzzz_xyyzz_1[j];

                    tg_xzzzz_xxyzzz_0[j] = pb_x * tg_zzzz_xxyzzz_0[j] + fr * tg_zzzz_xxyzzz_1[j] + fl1_fxn * tg_zzzz_xyzzz_1[j];

                    tg_xzzzz_xxzzzz_0[j] = pb_x * tg_zzzz_xxzzzz_0[j] + fr * tg_zzzz_xxzzzz_1[j] + fl1_fxn * tg_zzzz_xzzzz_1[j];

                    tg_xzzzz_xyyyyy_0[j] = pb_x * tg_zzzz_xyyyyy_0[j] + fr * tg_zzzz_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyy_1[j];

                    tg_xzzzz_xyyyyz_0[j] = pb_x * tg_zzzz_xyyyyz_0[j] + fr * tg_zzzz_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyyz_1[j];

                    tg_xzzzz_xyyyzz_0[j] = pb_x * tg_zzzz_xyyyzz_0[j] + fr * tg_zzzz_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyyzz_1[j];

                    tg_xzzzz_xyyzzz_0[j] = pb_x * tg_zzzz_xyyzzz_0[j] + fr * tg_zzzz_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyzzz_1[j];

                    tg_xzzzz_xyzzzz_0[j] = pb_x * tg_zzzz_xyzzzz_0[j] + fr * tg_zzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yzzzz_1[j];

                    tg_xzzzz_xzzzzz_0[j] = pb_x * tg_zzzz_xzzzzz_0[j] + fr * tg_zzzz_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzzz_1[j];

                    tg_xzzzz_yyyyyy_0[j] = pb_x * tg_zzzz_yyyyyy_0[j] + fr * tg_zzzz_yyyyyy_1[j];

                    tg_xzzzz_yyyyyz_0[j] = pb_x * tg_zzzz_yyyyyz_0[j] + fr * tg_zzzz_yyyyyz_1[j];

                    tg_xzzzz_yyyyzz_0[j] = pb_x * tg_zzzz_yyyyzz_0[j] + fr * tg_zzzz_yyyyzz_1[j];

                    tg_xzzzz_yyyzzz_0[j] = pb_x * tg_zzzz_yyyzzz_0[j] + fr * tg_zzzz_yyyzzz_1[j];

                    tg_xzzzz_yyzzzz_0[j] = pb_x * tg_zzzz_yyzzzz_0[j] + fr * tg_zzzz_yyzzzz_1[j];

                    tg_xzzzz_yzzzzz_0[j] = pb_x * tg_zzzz_yzzzzz_0[j] + fr * tg_zzzz_yzzzzz_1[j];

                    tg_xzzzz_zzzzzz_0[j] = pb_x * tg_zzzz_zzzzzz_0[j] + fr * tg_zzzz_zzzzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyy_xxxxxx_0[j] = pb_y * tg_yyyy_xxxxxx_0[j] + fr * tg_yyyy_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxx_0[j] - tg_yyy_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyy_xxxxxy_0[j] = pb_y * tg_yyyy_xxxxxy_0[j] + fr * tg_yyyy_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxy_0[j] - tg_yyy_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxxx_1[j];

                    tg_yyyyy_xxxxxz_0[j] = pb_y * tg_yyyy_xxxxxz_0[j] + fr * tg_yyyy_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxxz_0[j] - tg_yyy_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyy_xxxxyy_0[j] = pb_y * tg_yyyy_xxxxyy_0[j] + fr * tg_yyyy_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxyy_0[j] - tg_yyy_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxxxy_1[j];

                    tg_yyyyy_xxxxyz_0[j] = pb_y * tg_yyyy_xxxxyz_0[j] + fr * tg_yyyy_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxyz_0[j] - tg_yyy_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxxz_1[j];

                    tg_yyyyy_xxxxzz_0[j] = pb_y * tg_yyyy_xxxxzz_0[j] + fr * tg_yyyy_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxxzz_0[j] - tg_yyy_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyy_xxxyyy_0[j] = pb_y * tg_yyyy_xxxyyy_0[j] + fr * tg_yyyy_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyyy_0[j] - tg_yyy_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xxxyy_1[j];

                    tg_yyyyy_xxxyyz_0[j] = pb_y * tg_yyyy_xxxyyz_0[j] + fr * tg_yyyy_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyyz_0[j] - tg_yyy_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxxyz_1[j];

                    tg_yyyyy_xxxyzz_0[j] = pb_y * tg_yyyy_xxxyzz_0[j] + fr * tg_yyyy_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxyzz_0[j] - tg_yyy_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxxzz_1[j];

                    tg_yyyyy_xxxzzz_0[j] = pb_y * tg_yyyy_xxxzzz_0[j] + fr * tg_yyyy_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxxzzz_0[j] - tg_yyy_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyy_xxyyyy_0[j] = pb_y * tg_yyyy_xxyyyy_0[j] + fr * tg_yyyy_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyyy_0[j] - tg_yyy_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_xxyyy_1[j];

                    tg_yyyyy_xxyyyz_0[j] = pb_y * tg_yyyy_xxyyyz_0[j] + fr * tg_yyyy_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyyz_0[j] - tg_yyy_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xxyyz_1[j];

                    tg_yyyyy_xxyyzz_0[j] = pb_y * tg_yyyy_xxyyzz_0[j] + fr * tg_yyyy_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyyzz_0[j] - tg_yyy_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xxyzz_1[j];

                    tg_yyyyy_xxyzzz_0[j] = pb_y * tg_yyyy_xxyzzz_0[j] + fr * tg_yyyy_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxyzzz_0[j] - tg_yyy_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xxzzz_1[j];

                    tg_yyyyy_xxzzzz_0[j] = pb_y * tg_yyyy_xxzzzz_0[j] + fr * tg_yyyy_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xxzzzz_0[j] - tg_yyy_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyy_xyyyyy_0[j] = pb_y * tg_yyyy_xyyyyy_0[j] + fr * tg_yyyy_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyyy_0[j] - tg_yyy_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyy_xyyyy_1[j];

                    tg_yyyyy_xyyyyz_0[j] = pb_y * tg_yyyy_xyyyyz_0[j] + fr * tg_yyyy_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyyz_0[j] - tg_yyy_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_xyyyz_1[j];

                    tg_yyyyy_xyyyzz_0[j] = pb_y * tg_yyyy_xyyyzz_0[j] + fr * tg_yyyy_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyyzz_0[j] - tg_yyy_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_xyyzz_1[j];

                    tg_yyyyy_xyyzzz_0[j] = pb_y * tg_yyyy_xyyzzz_0[j] + fr * tg_yyyy_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyyzzz_0[j] - tg_yyy_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_xyzzz_1[j];

                    tg_yyyyy_xyzzzz_0[j] = pb_y * tg_yyyy_xyzzzz_0[j] + fr * tg_yyyy_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xyzzzz_0[j] - tg_yyy_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_xzzzz_1[j];

                    tg_yyyyy_xzzzzz_0[j] = pb_y * tg_yyyy_xzzzzz_0[j] + fr * tg_yyyy_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_xzzzzz_0[j] - tg_yyy_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyy_yyyyyy_0[j] = pb_y * tg_yyyy_yyyyyy_0[j] + fr * tg_yyyy_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyyy_0[j] - tg_yyy_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyy_yyyyy_1[j];

                    tg_yyyyy_yyyyyz_0[j] = pb_y * tg_yyyy_yyyyyz_0[j] + fr * tg_yyyy_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyyz_0[j] - tg_yyy_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyy_yyyyz_1[j];

                    tg_yyyyy_yyyyzz_0[j] = pb_y * tg_yyyy_yyyyzz_0[j] + fr * tg_yyyy_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyyzz_0[j] - tg_yyy_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyy_yyyzz_1[j];

                    tg_yyyyy_yyyzzz_0[j] = pb_y * tg_yyyy_yyyzzz_0[j] + fr * tg_yyyy_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyyzzz_0[j] - tg_yyy_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyy_yyzzz_1[j];

                    tg_yyyyy_yyzzzz_0[j] = pb_y * tg_yyyy_yyzzzz_0[j] + fr * tg_yyyy_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yyzzzz_0[j] - tg_yyy_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyy_yzzzz_1[j];

                    tg_yyyyy_yzzzzz_0[j] = pb_y * tg_yyyy_yzzzzz_0[j] + fr * tg_yyyy_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_yzzzzz_0[j] - tg_yyy_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyy_zzzzz_1[j];

                    tg_yyyyy_zzzzzz_0[j] = pb_y * tg_yyyy_zzzzzz_0[j] + fr * tg_yyyy_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_yyy_zzzzzz_0[j] - tg_yyy_zzzzzz_1[j] * fl1_fza);

                    tg_yyyyz_xxxxxx_0[j] = pb_y * tg_yyyz_xxxxxx_0[j] + fr * tg_yyyz_xxxxxx_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxx_0[j] - tg_yyz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyyz_xxxxxy_0[j] = pb_y * tg_yyyz_xxxxxy_0[j] + fr * tg_yyyz_xxxxxy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxy_0[j] - tg_yyz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxxx_1[j];

                    tg_yyyyz_xxxxxz_0[j] = pb_y * tg_yyyz_xxxxxz_0[j] + fr * tg_yyyz_xxxxxz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxxz_0[j] - tg_yyz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyyz_xxxxyy_0[j] = pb_y * tg_yyyz_xxxxyy_0[j] + fr * tg_yyyz_xxxxyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxyy_0[j] - tg_yyz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxxxy_1[j];

                    tg_yyyyz_xxxxyz_0[j] = pb_y * tg_yyyz_xxxxyz_0[j] + fr * tg_yyyz_xxxxyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxyz_0[j] - tg_yyz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxxz_1[j];

                    tg_yyyyz_xxxxzz_0[j] = pb_y * tg_yyyz_xxxxzz_0[j] + fr * tg_yyyz_xxxxzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxxzz_0[j] - tg_yyz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyyz_xxxyyy_0[j] = pb_y * tg_yyyz_xxxyyy_0[j] + fr * tg_yyyz_xxxyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyyy_0[j] - tg_yyz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xxxyy_1[j];

                    tg_yyyyz_xxxyyz_0[j] = pb_y * tg_yyyz_xxxyyz_0[j] + fr * tg_yyyz_xxxyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyyz_0[j] - tg_yyz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxxyz_1[j];

                    tg_yyyyz_xxxyzz_0[j] = pb_y * tg_yyyz_xxxyzz_0[j] + fr * tg_yyyz_xxxyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxyzz_0[j] - tg_yyz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxxzz_1[j];

                    tg_yyyyz_xxxzzz_0[j] = pb_y * tg_yyyz_xxxzzz_0[j] + fr * tg_yyyz_xxxzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxxzzz_0[j] - tg_yyz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyyz_xxyyyy_0[j] = pb_y * tg_yyyz_xxyyyy_0[j] + fr * tg_yyyz_xxyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyyy_0[j] - tg_yyz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_xxyyy_1[j];

                    tg_yyyyz_xxyyyz_0[j] = pb_y * tg_yyyz_xxyyyz_0[j] + fr * tg_yyyz_xxyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyyz_0[j] - tg_yyz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xxyyz_1[j];

                    tg_yyyyz_xxyyzz_0[j] = pb_y * tg_yyyz_xxyyzz_0[j] + fr * tg_yyyz_xxyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyyzz_0[j] - tg_yyz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xxyzz_1[j];

                    tg_yyyyz_xxyzzz_0[j] = pb_y * tg_yyyz_xxyzzz_0[j] + fr * tg_yyyz_xxyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxyzzz_0[j] - tg_yyz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xxzzz_1[j];

                    tg_yyyyz_xxzzzz_0[j] = pb_y * tg_yyyz_xxzzzz_0[j] + fr * tg_yyyz_xxzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xxzzzz_0[j] - tg_yyz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyyz_xyyyyy_0[j] = pb_y * tg_yyyz_xyyyyy_0[j] + fr * tg_yyyz_xyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyyy_0[j] - tg_yyz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyz_xyyyy_1[j];

                    tg_yyyyz_xyyyyz_0[j] = pb_y * tg_yyyz_xyyyyz_0[j] + fr * tg_yyyz_xyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyyz_0[j] - tg_yyz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_xyyyz_1[j];

                    tg_yyyyz_xyyyzz_0[j] = pb_y * tg_yyyz_xyyyzz_0[j] + fr * tg_yyyz_xyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyyzz_0[j] - tg_yyz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_xyyzz_1[j];

                    tg_yyyyz_xyyzzz_0[j] = pb_y * tg_yyyz_xyyzzz_0[j] + fr * tg_yyyz_xyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyyzzz_0[j] - tg_yyz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_xyzzz_1[j];

                    tg_yyyyz_xyzzzz_0[j] = pb_y * tg_yyyz_xyzzzz_0[j] + fr * tg_yyyz_xyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xyzzzz_0[j] - tg_yyz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_xzzzz_1[j];

                    tg_yyyyz_xzzzzz_0[j] = pb_y * tg_yyyz_xzzzzz_0[j] + fr * tg_yyyz_xzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_xzzzzz_0[j] - tg_yyz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyyz_yyyyyy_0[j] = pb_y * tg_yyyz_yyyyyy_0[j] + fr * tg_yyyz_yyyyyy_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyyy_0[j] - tg_yyz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyyz_yyyyy_1[j];

                    tg_yyyyz_yyyyyz_0[j] = pb_y * tg_yyyz_yyyyyz_0[j] + fr * tg_yyyz_yyyyyz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyyz_0[j] - tg_yyz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyyz_yyyyz_1[j];

                    tg_yyyyz_yyyyzz_0[j] = pb_y * tg_yyyz_yyyyzz_0[j] + fr * tg_yyyz_yyyyzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyyzz_0[j] - tg_yyz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyz_yyyzz_1[j];

                    tg_yyyyz_yyyzzz_0[j] = pb_y * tg_yyyz_yyyzzz_0[j] + fr * tg_yyyz_yyyzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyyzzz_0[j] - tg_yyz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyz_yyzzz_1[j];

                    tg_yyyyz_yyzzzz_0[j] = pb_y * tg_yyyz_yyzzzz_0[j] + fr * tg_yyyz_yyzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yyzzzz_0[j] - tg_yyz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyz_yzzzz_1[j];

                    tg_yyyyz_yzzzzz_0[j] = pb_y * tg_yyyz_yzzzzz_0[j] + fr * tg_yyyz_yzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_yzzzzz_0[j] - tg_yyz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyz_zzzzz_1[j];

                    tg_yyyyz_zzzzzz_0[j] = pb_y * tg_yyyz_zzzzzz_0[j] + fr * tg_yyyz_zzzzzz_1[j] + 1.5 * fl1_fx * (tg_yyz_zzzzzz_0[j] - tg_yyz_zzzzzz_1[j] * fl1_fza);

                    tg_yyyzz_xxxxxx_0[j] = pb_y * tg_yyzz_xxxxxx_0[j] + fr * tg_yyzz_xxxxxx_1[j] + fl1_fx * (tg_yzz_xxxxxx_0[j] - tg_yzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyyzz_xxxxxy_0[j] = pb_y * tg_yyzz_xxxxxy_0[j] + fr * tg_yyzz_xxxxxy_1[j] + fl1_fx * (tg_yzz_xxxxxy_0[j] - tg_yzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxxx_1[j];

                    tg_yyyzz_xxxxxz_0[j] = pb_y * tg_yyzz_xxxxxz_0[j] + fr * tg_yyzz_xxxxxz_1[j] + fl1_fx * (tg_yzz_xxxxxz_0[j] - tg_yzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyyzz_xxxxyy_0[j] = pb_y * tg_yyzz_xxxxyy_0[j] + fr * tg_yyzz_xxxxyy_1[j] + fl1_fx * (tg_yzz_xxxxyy_0[j] - tg_yzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxxxy_1[j];

                    tg_yyyzz_xxxxyz_0[j] = pb_y * tg_yyzz_xxxxyz_0[j] + fr * tg_yyzz_xxxxyz_1[j] + fl1_fx * (tg_yzz_xxxxyz_0[j] - tg_yzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxxz_1[j];

                    tg_yyyzz_xxxxzz_0[j] = pb_y * tg_yyzz_xxxxzz_0[j] + fr * tg_yyzz_xxxxzz_1[j] + fl1_fx * (tg_yzz_xxxxzz_0[j] - tg_yzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyyzz_xxxyyy_0[j] = pb_y * tg_yyzz_xxxyyy_0[j] + fr * tg_yyzz_xxxyyy_1[j] + fl1_fx * (tg_yzz_xxxyyy_0[j] - tg_yzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xxxyy_1[j];

                    tg_yyyzz_xxxyyz_0[j] = pb_y * tg_yyzz_xxxyyz_0[j] + fr * tg_yyzz_xxxyyz_1[j] + fl1_fx * (tg_yzz_xxxyyz_0[j] - tg_yzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxxyz_1[j];

                    tg_yyyzz_xxxyzz_0[j] = pb_y * tg_yyzz_xxxyzz_0[j] + fr * tg_yyzz_xxxyzz_1[j] + fl1_fx * (tg_yzz_xxxyzz_0[j] - tg_yzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxxzz_1[j];

                    tg_yyyzz_xxxzzz_0[j] = pb_y * tg_yyzz_xxxzzz_0[j] + fr * tg_yyzz_xxxzzz_1[j] + fl1_fx * (tg_yzz_xxxzzz_0[j] - tg_yzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyyzz_xxyyyy_0[j] = pb_y * tg_yyzz_xxyyyy_0[j] + fr * tg_yyzz_xxyyyy_1[j] + fl1_fx * (tg_yzz_xxyyyy_0[j] - tg_yzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_xxyyy_1[j];

                    tg_yyyzz_xxyyyz_0[j] = pb_y * tg_yyzz_xxyyyz_0[j] + fr * tg_yyzz_xxyyyz_1[j] + fl1_fx * (tg_yzz_xxyyyz_0[j] - tg_yzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xxyyz_1[j];

                    tg_yyyzz_xxyyzz_0[j] = pb_y * tg_yyzz_xxyyzz_0[j] + fr * tg_yyzz_xxyyzz_1[j] + fl1_fx * (tg_yzz_xxyyzz_0[j] - tg_yzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xxyzz_1[j];

                    tg_yyyzz_xxyzzz_0[j] = pb_y * tg_yyzz_xxyzzz_0[j] + fr * tg_yyzz_xxyzzz_1[j] + fl1_fx * (tg_yzz_xxyzzz_0[j] - tg_yzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xxzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSI_490_588(      CMemBlock2D<double>* primBuffer,
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
                                             {5, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_4_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {6, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

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

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

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

                // set up pointers to integrals

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

                // Batch of Integrals (490,588)

                #pragma omp simd aligned(fxn, fza, tg_yyyzz_xxzzzz_0, tg_yyyzz_xyyyyy_0, tg_yyyzz_xyyyyz_0, \
                                         tg_yyyzz_xyyyzz_0, tg_yyyzz_xyyzzz_0, tg_yyyzz_xyzzzz_0, tg_yyyzz_xzzzzz_0, \
                                         tg_yyyzz_yyyyyy_0, tg_yyyzz_yyyyyz_0, tg_yyyzz_yyyyzz_0, tg_yyyzz_yyyzzz_0, \
                                         tg_yyyzz_yyzzzz_0, tg_yyyzz_yzzzzz_0, tg_yyyzz_zzzzzz_0, tg_yyzz_xxzzzz_0, \
                                         tg_yyzz_xxzzzz_1, tg_yyzz_xyyyy_1, tg_yyzz_xyyyyy_0, tg_yyzz_xyyyyy_1, \
                                         tg_yyzz_xyyyyz_0, tg_yyzz_xyyyyz_1, tg_yyzz_xyyyz_1, tg_yyzz_xyyyzz_0, \
                                         tg_yyzz_xyyyzz_1, tg_yyzz_xyyzz_1, tg_yyzz_xyyzzz_0, tg_yyzz_xyyzzz_1, \
                                         tg_yyzz_xyzzz_1, tg_yyzz_xyzzzz_0, tg_yyzz_xyzzzz_1, tg_yyzz_xzzzz_1, \
                                         tg_yyzz_xzzzzz_0, tg_yyzz_xzzzzz_1, tg_yyzz_yyyyy_1, tg_yyzz_yyyyyy_0, \
                                         tg_yyzz_yyyyyy_1, tg_yyzz_yyyyyz_0, tg_yyzz_yyyyyz_1, tg_yyzz_yyyyz_1, \
                                         tg_yyzz_yyyyzz_0, tg_yyzz_yyyyzz_1, tg_yyzz_yyyzz_1, tg_yyzz_yyyzzz_0, \
                                         tg_yyzz_yyyzzz_1, tg_yyzz_yyzzz_1, tg_yyzz_yyzzzz_0, tg_yyzz_yyzzzz_1, \
                                         tg_yyzz_yzzzz_1, tg_yyzz_yzzzzz_0, tg_yyzz_yzzzzz_1, tg_yyzz_zzzzz_1, \
                                         tg_yyzz_zzzzzz_0, tg_yyzz_zzzzzz_1, tg_yyzzz_xxxxxx_0, tg_yyzzz_xxxxxy_0, \
                                         tg_yyzzz_xxxxxz_0, tg_yyzzz_xxxxyy_0, tg_yyzzz_xxxxyz_0, tg_yyzzz_xxxxzz_0, \
                                         tg_yyzzz_xxxyyy_0, tg_yyzzz_xxxyyz_0, tg_yyzzz_xxxyzz_0, tg_yyzzz_xxxzzz_0, \
                                         tg_yyzzz_xxyyyy_0, tg_yyzzz_xxyyyz_0, tg_yyzzz_xxyyzz_0, tg_yyzzz_xxyzzz_0, \
                                         tg_yyzzz_xxzzzz_0, tg_yyzzz_xyyyyy_0, tg_yyzzz_xyyyyz_0, tg_yyzzz_xyyyzz_0, \
                                         tg_yyzzz_xyyzzz_0, tg_yyzzz_xyzzzz_0, tg_yyzzz_xzzzzz_0, tg_yyzzz_yyyyyy_0, \
                                         tg_yyzzz_yyyyyz_0, tg_yyzzz_yyyyzz_0, tg_yyzzz_yyyzzz_0, tg_yyzzz_yyzzzz_0, \
                                         tg_yyzzz_yzzzzz_0, tg_yyzzz_zzzzzz_0, tg_yzz_xxzzzz_0, tg_yzz_xxzzzz_1, \
                                         tg_yzz_xyyyyy_0, tg_yzz_xyyyyy_1, tg_yzz_xyyyyz_0, tg_yzz_xyyyyz_1, tg_yzz_xyyyzz_0, \
                                         tg_yzz_xyyyzz_1, tg_yzz_xyyzzz_0, tg_yzz_xyyzzz_1, tg_yzz_xyzzzz_0, tg_yzz_xyzzzz_1, \
                                         tg_yzz_xzzzzz_0, tg_yzz_xzzzzz_1, tg_yzz_yyyyyy_0, tg_yzz_yyyyyy_1, tg_yzz_yyyyyz_0, \
                                         tg_yzz_yyyyyz_1, tg_yzz_yyyyzz_0, tg_yzz_yyyyzz_1, tg_yzz_yyyzzz_0, tg_yzz_yyyzzz_1, \
                                         tg_yzz_yyzzzz_0, tg_yzz_yyzzzz_1, tg_yzz_yzzzzz_0, tg_yzz_yzzzzz_1, tg_yzz_zzzzzz_0, \
                                         tg_yzz_zzzzzz_1, tg_yzzz_xxxxx_1, tg_yzzz_xxxxxx_0, tg_yzzz_xxxxxx_1, \
                                         tg_yzzz_xxxxxy_0, tg_yzzz_xxxxxy_1, tg_yzzz_xxxxxz_0, tg_yzzz_xxxxxz_1, \
                                         tg_yzzz_xxxxy_1, tg_yzzz_xxxxyy_0, tg_yzzz_xxxxyy_1, tg_yzzz_xxxxyz_0, \
                                         tg_yzzz_xxxxyz_1, tg_yzzz_xxxxz_1, tg_yzzz_xxxxzz_0, tg_yzzz_xxxxzz_1, \
                                         tg_yzzz_xxxyy_1, tg_yzzz_xxxyyy_0, tg_yzzz_xxxyyy_1, tg_yzzz_xxxyyz_0, \
                                         tg_yzzz_xxxyyz_1, tg_yzzz_xxxyz_1, tg_yzzz_xxxyzz_0, tg_yzzz_xxxyzz_1, \
                                         tg_yzzz_xxxzz_1, tg_yzzz_xxxzzz_0, tg_yzzz_xxxzzz_1, tg_yzzz_xxyyy_1, \
                                         tg_yzzz_xxyyyy_0, tg_yzzz_xxyyyy_1, tg_yzzz_xxyyyz_0, tg_yzzz_xxyyyz_1, \
                                         tg_yzzz_xxyyz_1, tg_yzzz_xxyyzz_0, tg_yzzz_xxyyzz_1, tg_yzzz_xxyzz_1, \
                                         tg_yzzz_xxyzzz_0, tg_yzzz_xxyzzz_1, tg_yzzz_xxzzz_1, tg_yzzz_xxzzzz_0, \
                                         tg_yzzz_xxzzzz_1, tg_yzzz_xyyyy_1, tg_yzzz_xyyyyy_0, tg_yzzz_xyyyyy_1, \
                                         tg_yzzz_xyyyyz_0, tg_yzzz_xyyyyz_1, tg_yzzz_xyyyz_1, tg_yzzz_xyyyzz_0, \
                                         tg_yzzz_xyyyzz_1, tg_yzzz_xyyzz_1, tg_yzzz_xyyzzz_0, tg_yzzz_xyyzzz_1, \
                                         tg_yzzz_xyzzz_1, tg_yzzz_xyzzzz_0, tg_yzzz_xyzzzz_1, tg_yzzz_xzzzz_1, \
                                         tg_yzzz_xzzzzz_0, tg_yzzz_xzzzzz_1, tg_yzzz_yyyyy_1, tg_yzzz_yyyyyy_0, \
                                         tg_yzzz_yyyyyy_1, tg_yzzz_yyyyyz_0, tg_yzzz_yyyyyz_1, tg_yzzz_yyyyz_1, \
                                         tg_yzzz_yyyyzz_0, tg_yzzz_yyyyzz_1, tg_yzzz_yyyzz_1, tg_yzzz_yyyzzz_0, \
                                         tg_yzzz_yyyzzz_1, tg_yzzz_yyzzz_1, tg_yzzz_yyzzzz_0, tg_yzzz_yyzzzz_1, \
                                         tg_yzzz_yzzzz_1, tg_yzzz_yzzzzz_0, tg_yzzz_yzzzzz_1, tg_yzzz_zzzzz_1, \
                                         tg_yzzz_zzzzzz_0, tg_yzzz_zzzzzz_1, tg_yzzzz_xxxxxx_0, tg_yzzzz_xxxxxy_0, \
                                         tg_yzzzz_xxxxxz_0, tg_yzzzz_xxxxyy_0, tg_yzzzz_xxxxyz_0, tg_yzzzz_xxxxzz_0, \
                                         tg_yzzzz_xxxyyy_0, tg_yzzzz_xxxyyz_0, tg_yzzzz_xxxyzz_0, tg_yzzzz_xxxzzz_0, \
                                         tg_yzzzz_xxyyyy_0, tg_yzzzz_xxyyyz_0, tg_yzzzz_xxyyzz_0, tg_yzzzz_xxyzzz_0, \
                                         tg_yzzzz_xxzzzz_0, tg_yzzzz_xyyyyy_0, tg_yzzzz_xyyyyz_0, tg_yzzzz_xyyyzz_0, \
                                         tg_yzzzz_xyyzzz_0, tg_yzzzz_xyzzzz_0, tg_yzzzz_xzzzzz_0, tg_yzzzz_yyyyyy_0, \
                                         tg_yzzzz_yyyyyz_0, tg_yzzzz_yyyyzz_0, tg_yzzzz_yyyzzz_0, tg_yzzzz_yyzzzz_0, \
                                         tg_yzzzz_yzzzzz_0, tg_yzzzz_zzzzzz_0, tg_zzz_xxxxxx_0, tg_zzz_xxxxxx_1, \
                                         tg_zzz_xxxxxy_0, tg_zzz_xxxxxy_1, tg_zzz_xxxxxz_0, tg_zzz_xxxxxz_1, tg_zzz_xxxxyy_0, \
                                         tg_zzz_xxxxyy_1, tg_zzz_xxxxyz_0, tg_zzz_xxxxyz_1, tg_zzz_xxxxzz_0, tg_zzz_xxxxzz_1, \
                                         tg_zzz_xxxyyy_0, tg_zzz_xxxyyy_1, tg_zzz_xxxyyz_0, tg_zzz_xxxyyz_1, tg_zzz_xxxyzz_0, \
                                         tg_zzz_xxxyzz_1, tg_zzz_xxxzzz_0, tg_zzz_xxxzzz_1, tg_zzz_xxyyyy_0, tg_zzz_xxyyyy_1, \
                                         tg_zzz_xxyyyz_0, tg_zzz_xxyyyz_1, tg_zzz_xxyyzz_0, tg_zzz_xxyyzz_1, tg_zzz_xxyzzz_0, \
                                         tg_zzz_xxyzzz_1, tg_zzz_xxzzzz_0, tg_zzz_xxzzzz_1, tg_zzz_xyyyyy_0, tg_zzz_xyyyyy_1, \
                                         tg_zzz_xyyyyz_0, tg_zzz_xyyyyz_1, tg_zzz_xyyyzz_0, tg_zzz_xyyyzz_1, tg_zzz_xyyzzz_0, \
                                         tg_zzz_xyyzzz_1, tg_zzz_xyzzzz_0, tg_zzz_xyzzzz_1, tg_zzz_xzzzzz_0, tg_zzz_xzzzzz_1, \
                                         tg_zzz_yyyyyy_0, tg_zzz_yyyyyy_1, tg_zzz_yyyyyz_0, tg_zzz_yyyyyz_1, tg_zzz_yyyyzz_0, \
                                         tg_zzz_yyyyzz_1, tg_zzz_yyyzzz_0, tg_zzz_yyyzzz_1, tg_zzz_yyzzzz_0, tg_zzz_yyzzzz_1, \
                                         tg_zzz_yzzzzz_0, tg_zzz_yzzzzz_1, tg_zzz_zzzzzz_0, tg_zzz_zzzzzz_1, tg_zzzz_xxxxx_1, \
                                         tg_zzzz_xxxxxx_0, tg_zzzz_xxxxxx_1, tg_zzzz_xxxxxy_0, tg_zzzz_xxxxxy_1, \
                                         tg_zzzz_xxxxxz_0, tg_zzzz_xxxxxz_1, tg_zzzz_xxxxy_1, tg_zzzz_xxxxyy_0, \
                                         tg_zzzz_xxxxyy_1, tg_zzzz_xxxxyz_0, tg_zzzz_xxxxyz_1, tg_zzzz_xxxxz_1, \
                                         tg_zzzz_xxxxzz_0, tg_zzzz_xxxxzz_1, tg_zzzz_xxxyy_1, tg_zzzz_xxxyyy_0, \
                                         tg_zzzz_xxxyyy_1, tg_zzzz_xxxyyz_0, tg_zzzz_xxxyyz_1, tg_zzzz_xxxyz_1, \
                                         tg_zzzz_xxxyzz_0, tg_zzzz_xxxyzz_1, tg_zzzz_xxxzz_1, tg_zzzz_xxxzzz_0, \
                                         tg_zzzz_xxxzzz_1, tg_zzzz_xxyyy_1, tg_zzzz_xxyyyy_0, tg_zzzz_xxyyyy_1, \
                                         tg_zzzz_xxyyyz_0, tg_zzzz_xxyyyz_1, tg_zzzz_xxyyz_1, tg_zzzz_xxyyzz_0, \
                                         tg_zzzz_xxyyzz_1, tg_zzzz_xxyzz_1, tg_zzzz_xxyzzz_0, tg_zzzz_xxyzzz_1, \
                                         tg_zzzz_xxzzz_1, tg_zzzz_xxzzzz_0, tg_zzzz_xxzzzz_1, tg_zzzz_xyyyy_1, \
                                         tg_zzzz_xyyyyy_0, tg_zzzz_xyyyyy_1, tg_zzzz_xyyyyz_0, tg_zzzz_xyyyyz_1, \
                                         tg_zzzz_xyyyz_1, tg_zzzz_xyyyzz_0, tg_zzzz_xyyyzz_1, tg_zzzz_xyyzz_1, \
                                         tg_zzzz_xyyzzz_0, tg_zzzz_xyyzzz_1, tg_zzzz_xyzzz_1, tg_zzzz_xyzzzz_0, \
                                         tg_zzzz_xyzzzz_1, tg_zzzz_xzzzz_1, tg_zzzz_xzzzzz_0, tg_zzzz_xzzzzz_1, \
                                         tg_zzzz_yyyyy_1, tg_zzzz_yyyyyy_0, tg_zzzz_yyyyyy_1, tg_zzzz_yyyyyz_0, \
                                         tg_zzzz_yyyyyz_1, tg_zzzz_yyyyz_1, tg_zzzz_yyyyzz_0, tg_zzzz_yyyyzz_1, \
                                         tg_zzzz_yyyzz_1, tg_zzzz_yyyzzz_0, tg_zzzz_yyyzzz_1, tg_zzzz_yyzzz_1, \
                                         tg_zzzz_yyzzzz_0, tg_zzzz_yyzzzz_1, tg_zzzz_yzzzz_1, tg_zzzz_yzzzzz_0, \
                                         tg_zzzz_yzzzzz_1, tg_zzzz_zzzzz_1, tg_zzzz_zzzzzz_0, tg_zzzz_zzzzzz_1, \
                                         tg_zzzzz_xxxxxx_0, tg_zzzzz_xxxxxy_0, tg_zzzzz_xxxxxz_0, tg_zzzzz_xxxxyy_0, \
                                         tg_zzzzz_xxxxyz_0, tg_zzzzz_xxxxzz_0, tg_zzzzz_xxxyyy_0, tg_zzzzz_xxxyyz_0, \
                                         tg_zzzzz_xxxyzz_0, tg_zzzzz_xxxzzz_0, tg_zzzzz_xxyyyy_0, tg_zzzzz_xxyyyz_0, \
                                         tg_zzzzz_xxyyzz_0, tg_zzzzz_xxyzzz_0, tg_zzzzz_xxzzzz_0, tg_zzzzz_xyyyyy_0, \
                                         tg_zzzzz_xyyyyz_0, tg_zzzzz_xyyyzz_0, tg_zzzzz_xyyzzz_0, tg_zzzzz_xyzzzz_0, \
                                         tg_zzzzz_xzzzzz_0, tg_zzzzz_yyyyyy_0, tg_zzzzz_yyyyyz_0, tg_zzzzz_yyyyzz_0, \
                                         tg_zzzzz_yyyzzz_0, tg_zzzzz_yyzzzz_0, tg_zzzzz_yzzzzz_0, tg_zzzzz_zzzzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyzz_xxzzzz_0[j] = pb_y * tg_yyzz_xxzzzz_0[j] + fr * tg_yyzz_xxzzzz_1[j] + fl1_fx * (tg_yzz_xxzzzz_0[j] - tg_yzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyyzz_xyyyyy_0[j] = pb_y * tg_yyzz_xyyyyy_0[j] + fr * tg_yyzz_xyyyyy_1[j] + fl1_fx * (tg_yzz_xyyyyy_0[j] - tg_yzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzz_xyyyy_1[j];

                    tg_yyyzz_xyyyyz_0[j] = pb_y * tg_yyzz_xyyyyz_0[j] + fr * tg_yyzz_xyyyyz_1[j] + fl1_fx * (tg_yzz_xyyyyz_0[j] - tg_yzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_xyyyz_1[j];

                    tg_yyyzz_xyyyzz_0[j] = pb_y * tg_yyzz_xyyyzz_0[j] + fr * tg_yyzz_xyyyzz_1[j] + fl1_fx * (tg_yzz_xyyyzz_0[j] - tg_yzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_xyyzz_1[j];

                    tg_yyyzz_xyyzzz_0[j] = pb_y * tg_yyzz_xyyzzz_0[j] + fr * tg_yyzz_xyyzzz_1[j] + fl1_fx * (tg_yzz_xyyzzz_0[j] - tg_yzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_xyzzz_1[j];

                    tg_yyyzz_xyzzzz_0[j] = pb_y * tg_yyzz_xyzzzz_0[j] + fr * tg_yyzz_xyzzzz_1[j] + fl1_fx * (tg_yzz_xyzzzz_0[j] - tg_yzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_xzzzz_1[j];

                    tg_yyyzz_xzzzzz_0[j] = pb_y * tg_yyzz_xzzzzz_0[j] + fr * tg_yyzz_xzzzzz_1[j] + fl1_fx * (tg_yzz_xzzzzz_0[j] - tg_yzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyyzz_yyyyyy_0[j] = pb_y * tg_yyzz_yyyyyy_0[j] + fr * tg_yyzz_yyyyyy_1[j] + fl1_fx * (tg_yzz_yyyyyy_0[j] - tg_yzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yyzz_yyyyy_1[j];

                    tg_yyyzz_yyyyyz_0[j] = pb_y * tg_yyzz_yyyyyz_0[j] + fr * tg_yyzz_yyyyyz_1[j] + fl1_fx * (tg_yzz_yyyyyz_0[j] - tg_yzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yyzz_yyyyz_1[j];

                    tg_yyyzz_yyyyzz_0[j] = pb_y * tg_yyzz_yyyyzz_0[j] + fr * tg_yyzz_yyyyzz_1[j] + fl1_fx * (tg_yzz_yyyyzz_0[j] - tg_yzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzz_yyyzz_1[j];

                    tg_yyyzz_yyyzzz_0[j] = pb_y * tg_yyzz_yyyzzz_0[j] + fr * tg_yyzz_yyyzzz_1[j] + fl1_fx * (tg_yzz_yyyzzz_0[j] - tg_yzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzz_yyzzz_1[j];

                    tg_yyyzz_yyzzzz_0[j] = pb_y * tg_yyzz_yyzzzz_0[j] + fr * tg_yyzz_yyzzzz_1[j] + fl1_fx * (tg_yzz_yyzzzz_0[j] - tg_yzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzz_yzzzz_1[j];

                    tg_yyyzz_yzzzzz_0[j] = pb_y * tg_yyzz_yzzzzz_0[j] + fr * tg_yyzz_yzzzzz_1[j] + fl1_fx * (tg_yzz_yzzzzz_0[j] - tg_yzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzz_zzzzz_1[j];

                    tg_yyyzz_zzzzzz_0[j] = pb_y * tg_yyzz_zzzzzz_0[j] + fr * tg_yyzz_zzzzzz_1[j] + fl1_fx * (tg_yzz_zzzzzz_0[j] - tg_yzz_zzzzzz_1[j] * fl1_fza);

                    tg_yyzzz_xxxxxx_0[j] = pb_y * tg_yzzz_xxxxxx_0[j] + fr * tg_yzzz_xxxxxx_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxx_0[j] - tg_zzz_xxxxxx_1[j] * fl1_fza);

                    tg_yyzzz_xxxxxy_0[j] = pb_y * tg_yzzz_xxxxxy_0[j] + fr * tg_yzzz_xxxxxy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxy_0[j] - tg_zzz_xxxxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxxx_1[j];

                    tg_yyzzz_xxxxxz_0[j] = pb_y * tg_yzzz_xxxxxz_0[j] + fr * tg_yzzz_xxxxxz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxxz_0[j] - tg_zzz_xxxxxz_1[j] * fl1_fza);

                    tg_yyzzz_xxxxyy_0[j] = pb_y * tg_yzzz_xxxxyy_0[j] + fr * tg_yzzz_xxxxyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyy_0[j] - tg_zzz_xxxxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxxxy_1[j];

                    tg_yyzzz_xxxxyz_0[j] = pb_y * tg_yzzz_xxxxyz_0[j] + fr * tg_yzzz_xxxxyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxyz_0[j] - tg_zzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxxz_1[j];

                    tg_yyzzz_xxxxzz_0[j] = pb_y * tg_yzzz_xxxxzz_0[j] + fr * tg_yzzz_xxxxzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxxzz_0[j] - tg_zzz_xxxxzz_1[j] * fl1_fza);

                    tg_yyzzz_xxxyyy_0[j] = pb_y * tg_yzzz_xxxyyy_0[j] + fr * tg_yzzz_xxxyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyy_0[j] - tg_zzz_xxxyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xxxyy_1[j];

                    tg_yyzzz_xxxyyz_0[j] = pb_y * tg_yzzz_xxxyyz_0[j] + fr * tg_yzzz_xxxyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyyz_0[j] - tg_zzz_xxxyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxxyz_1[j];

                    tg_yyzzz_xxxyzz_0[j] = pb_y * tg_yzzz_xxxyzz_0[j] + fr * tg_yzzz_xxxyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxyzz_0[j] - tg_zzz_xxxyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxxzz_1[j];

                    tg_yyzzz_xxxzzz_0[j] = pb_y * tg_yzzz_xxxzzz_0[j] + fr * tg_yzzz_xxxzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxxzzz_0[j] - tg_zzz_xxxzzz_1[j] * fl1_fza);

                    tg_yyzzz_xxyyyy_0[j] = pb_y * tg_yzzz_xxyyyy_0[j] + fr * tg_yzzz_xxyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyy_0[j] - tg_zzz_xxyyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_xxyyy_1[j];

                    tg_yyzzz_xxyyyz_0[j] = pb_y * tg_yzzz_xxyyyz_0[j] + fr * tg_yzzz_xxyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyyz_0[j] - tg_zzz_xxyyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xxyyz_1[j];

                    tg_yyzzz_xxyyzz_0[j] = pb_y * tg_yzzz_xxyyzz_0[j] + fr * tg_yzzz_xxyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyyzz_0[j] - tg_zzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xxyzz_1[j];

                    tg_yyzzz_xxyzzz_0[j] = pb_y * tg_yzzz_xxyzzz_0[j] + fr * tg_yzzz_xxyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxyzzz_0[j] - tg_zzz_xxyzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xxzzz_1[j];

                    tg_yyzzz_xxzzzz_0[j] = pb_y * tg_yzzz_xxzzzz_0[j] + fr * tg_yzzz_xxzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xxzzzz_0[j] - tg_zzz_xxzzzz_1[j] * fl1_fza);

                    tg_yyzzz_xyyyyy_0[j] = pb_y * tg_yzzz_xyyyyy_0[j] + fr * tg_yzzz_xyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyy_0[j] - tg_zzz_xyyyyy_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzz_xyyyy_1[j];

                    tg_yyzzz_xyyyyz_0[j] = pb_y * tg_yzzz_xyyyyz_0[j] + fr * tg_yzzz_xyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyyz_0[j] - tg_zzz_xyyyyz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_xyyyz_1[j];

                    tg_yyzzz_xyyyzz_0[j] = pb_y * tg_yzzz_xyyyzz_0[j] + fr * tg_yzzz_xyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyyzz_0[j] - tg_zzz_xyyyzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_xyyzz_1[j];

                    tg_yyzzz_xyyzzz_0[j] = pb_y * tg_yzzz_xyyzzz_0[j] + fr * tg_yzzz_xyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyyzzz_0[j] - tg_zzz_xyyzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_xyzzz_1[j];

                    tg_yyzzz_xyzzzz_0[j] = pb_y * tg_yzzz_xyzzzz_0[j] + fr * tg_yzzz_xyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xyzzzz_0[j] - tg_zzz_xyzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_xzzzz_1[j];

                    tg_yyzzz_xzzzzz_0[j] = pb_y * tg_yzzz_xzzzzz_0[j] + fr * tg_yzzz_xzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_xzzzzz_0[j] - tg_zzz_xzzzzz_1[j] * fl1_fza);

                    tg_yyzzz_yyyyyy_0[j] = pb_y * tg_yzzz_yyyyyy_0[j] + fr * tg_yzzz_yyyyyy_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyy_0[j] - tg_zzz_yyyyyy_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_yzzz_yyyyy_1[j];

                    tg_yyzzz_yyyyyz_0[j] = pb_y * tg_yzzz_yyyyyz_0[j] + fr * tg_yzzz_yyyyyz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyyz_0[j] - tg_zzz_yyyyyz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_yzzz_yyyyz_1[j];

                    tg_yyzzz_yyyyzz_0[j] = pb_y * tg_yzzz_yyyyzz_0[j] + fr * tg_yzzz_yyyyzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyyzz_0[j] - tg_zzz_yyyyzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzz_yyyzz_1[j];

                    tg_yyzzz_yyyzzz_0[j] = pb_y * tg_yzzz_yyyzzz_0[j] + fr * tg_yzzz_yyyzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyyzzz_0[j] - tg_zzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzz_yyzzz_1[j];

                    tg_yyzzz_yyzzzz_0[j] = pb_y * tg_yzzz_yyzzzz_0[j] + fr * tg_yzzz_yyzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yyzzzz_0[j] - tg_zzz_yyzzzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzz_yzzzz_1[j];

                    tg_yyzzz_yzzzzz_0[j] = pb_y * tg_yzzz_yzzzzz_0[j] + fr * tg_yzzz_yzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_yzzzzz_0[j] - tg_zzz_yzzzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzz_zzzzz_1[j];

                    tg_yyzzz_zzzzzz_0[j] = pb_y * tg_yzzz_zzzzzz_0[j] + fr * tg_yzzz_zzzzzz_1[j] + 0.5 * fl1_fx * (tg_zzz_zzzzzz_0[j] - tg_zzz_zzzzzz_1[j] * fl1_fza);

                    tg_yzzzz_xxxxxx_0[j] = pb_y * tg_zzzz_xxxxxx_0[j] + fr * tg_zzzz_xxxxxx_1[j];

                    tg_yzzzz_xxxxxy_0[j] = pb_y * tg_zzzz_xxxxxy_0[j] + fr * tg_zzzz_xxxxxy_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxx_1[j];

                    tg_yzzzz_xxxxxz_0[j] = pb_y * tg_zzzz_xxxxxz_0[j] + fr * tg_zzzz_xxxxxz_1[j];

                    tg_yzzzz_xxxxyy_0[j] = pb_y * tg_zzzz_xxxxyy_0[j] + fr * tg_zzzz_xxxxyy_1[j] + fl1_fxn * tg_zzzz_xxxxy_1[j];

                    tg_yzzzz_xxxxyz_0[j] = pb_y * tg_zzzz_xxxxyz_0[j] + fr * tg_zzzz_xxxxyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxxz_1[j];

                    tg_yzzzz_xxxxzz_0[j] = pb_y * tg_zzzz_xxxxzz_0[j] + fr * tg_zzzz_xxxxzz_1[j];

                    tg_yzzzz_xxxyyy_0[j] = pb_y * tg_zzzz_xxxyyy_0[j] + fr * tg_zzzz_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxxyy_1[j];

                    tg_yzzzz_xxxyyz_0[j] = pb_y * tg_zzzz_xxxyyz_0[j] + fr * tg_zzzz_xxxyyz_1[j] + fl1_fxn * tg_zzzz_xxxyz_1[j];

                    tg_yzzzz_xxxyzz_0[j] = pb_y * tg_zzzz_xxxyzz_0[j] + fr * tg_zzzz_xxxyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxxzz_1[j];

                    tg_yzzzz_xxxzzz_0[j] = pb_y * tg_zzzz_xxxzzz_0[j] + fr * tg_zzzz_xxxzzz_1[j];

                    tg_yzzzz_xxyyyy_0[j] = pb_y * tg_zzzz_xxyyyy_0[j] + fr * tg_zzzz_xxyyyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxyyy_1[j];

                    tg_yzzzz_xxyyyz_0[j] = pb_y * tg_zzzz_xxyyyz_0[j] + fr * tg_zzzz_xxyyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxyyz_1[j];

                    tg_yzzzz_xxyyzz_0[j] = pb_y * tg_zzzz_xxyyzz_0[j] + fr * tg_zzzz_xxyyzz_1[j] + fl1_fxn * tg_zzzz_xxyzz_1[j];

                    tg_yzzzz_xxyzzz_0[j] = pb_y * tg_zzzz_xxyzzz_0[j] + fr * tg_zzzz_xxyzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxzzz_1[j];

                    tg_yzzzz_xxzzzz_0[j] = pb_y * tg_zzzz_xxzzzz_0[j] + fr * tg_zzzz_xxzzzz_1[j];

                    tg_yzzzz_xyyyyy_0[j] = pb_y * tg_zzzz_xyyyyy_0[j] + fr * tg_zzzz_xyyyyy_1[j] + 2.5 * fl1_fxn * tg_zzzz_xyyyy_1[j];

                    tg_yzzzz_xyyyyz_0[j] = pb_y * tg_zzzz_xyyyyz_0[j] + fr * tg_zzzz_xyyyyz_1[j] + 2.0 * fl1_fxn * tg_zzzz_xyyyz_1[j];

                    tg_yzzzz_xyyyzz_0[j] = pb_y * tg_zzzz_xyyyzz_0[j] + fr * tg_zzzz_xyyyzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xyyzz_1[j];

                    tg_yzzzz_xyyzzz_0[j] = pb_y * tg_zzzz_xyyzzz_0[j] + fr * tg_zzzz_xyyzzz_1[j] + fl1_fxn * tg_zzzz_xyzzz_1[j];

                    tg_yzzzz_xyzzzz_0[j] = pb_y * tg_zzzz_xyzzzz_0[j] + fr * tg_zzzz_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xzzzz_1[j];

                    tg_yzzzz_xzzzzz_0[j] = pb_y * tg_zzzz_xzzzzz_0[j] + fr * tg_zzzz_xzzzzz_1[j];

                    tg_yzzzz_yyyyyy_0[j] = pb_y * tg_zzzz_yyyyyy_0[j] + fr * tg_zzzz_yyyyyy_1[j] + 3.0 * fl1_fxn * tg_zzzz_yyyyy_1[j];

                    tg_yzzzz_yyyyyz_0[j] = pb_y * tg_zzzz_yyyyyz_0[j] + fr * tg_zzzz_yyyyyz_1[j] + 2.5 * fl1_fxn * tg_zzzz_yyyyz_1[j];

                    tg_yzzzz_yyyyzz_0[j] = pb_y * tg_zzzz_yyyyzz_0[j] + fr * tg_zzzz_yyyyzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_yyyzz_1[j];

                    tg_yzzzz_yyyzzz_0[j] = pb_y * tg_zzzz_yyyzzz_0[j] + fr * tg_zzzz_yyyzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_yyzzz_1[j];

                    tg_yzzzz_yyzzzz_0[j] = pb_y * tg_zzzz_yyzzzz_0[j] + fr * tg_zzzz_yyzzzz_1[j] + fl1_fxn * tg_zzzz_yzzzz_1[j];

                    tg_yzzzz_yzzzzz_0[j] = pb_y * tg_zzzz_yzzzzz_0[j] + fr * tg_zzzz_yzzzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzzzz_1[j];

                    tg_yzzzz_zzzzzz_0[j] = pb_y * tg_zzzz_zzzzzz_0[j] + fr * tg_zzzz_zzzzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzz_xxxxxx_0[j] = pb_z * tg_zzzz_xxxxxx_0[j] + fr * tg_zzzz_xxxxxx_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxx_0[j] - tg_zzz_xxxxxx_1[j] * fl1_fza);

                    tg_zzzzz_xxxxxy_0[j] = pb_z * tg_zzzz_xxxxxy_0[j] + fr * tg_zzzz_xxxxxy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxy_0[j] - tg_zzz_xxxxxy_1[j] * fl1_fza);

                    tg_zzzzz_xxxxxz_0[j] = pb_z * tg_zzzz_xxxxxz_0[j] + fr * tg_zzzz_xxxxxz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxxz_0[j] - tg_zzz_xxxxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxxx_1[j];

                    tg_zzzzz_xxxxyy_0[j] = pb_z * tg_zzzz_xxxxyy_0[j] + fr * tg_zzzz_xxxxyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxyy_0[j] - tg_zzz_xxxxyy_1[j] * fl1_fza);

                    tg_zzzzz_xxxxyz_0[j] = pb_z * tg_zzzz_xxxxyz_0[j] + fr * tg_zzzz_xxxxyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxyz_0[j] - tg_zzz_xxxxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxxy_1[j];

                    tg_zzzzz_xxxxzz_0[j] = pb_z * tg_zzzz_xxxxzz_0[j] + fr * tg_zzzz_xxxxzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxxzz_0[j] - tg_zzz_xxxxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxxxz_1[j];

                    tg_zzzzz_xxxyyy_0[j] = pb_z * tg_zzzz_xxxyyy_0[j] + fr * tg_zzzz_xxxyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyyy_0[j] - tg_zzz_xxxyyy_1[j] * fl1_fza);

                    tg_zzzzz_xxxyyz_0[j] = pb_z * tg_zzzz_xxxyyz_0[j] + fr * tg_zzzz_xxxyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyyz_0[j] - tg_zzz_xxxyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxxyy_1[j];

                    tg_zzzzz_xxxyzz_0[j] = pb_z * tg_zzzz_xxxyzz_0[j] + fr * tg_zzzz_xxxyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxyzz_0[j] - tg_zzz_xxxyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxxyz_1[j];

                    tg_zzzzz_xxxzzz_0[j] = pb_z * tg_zzzz_xxxzzz_0[j] + fr * tg_zzzz_xxxzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxxzzz_0[j] - tg_zzz_xxxzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xxxzz_1[j];

                    tg_zzzzz_xxyyyy_0[j] = pb_z * tg_zzzz_xxyyyy_0[j] + fr * tg_zzzz_xxyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyyy_0[j] - tg_zzz_xxyyyy_1[j] * fl1_fza);

                    tg_zzzzz_xxyyyz_0[j] = pb_z * tg_zzzz_xxyyyz_0[j] + fr * tg_zzzz_xxyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyyz_0[j] - tg_zzz_xxyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xxyyy_1[j];

                    tg_zzzzz_xxyyzz_0[j] = pb_z * tg_zzzz_xxyyzz_0[j] + fr * tg_zzzz_xxyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyyzz_0[j] - tg_zzz_xxyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xxyyz_1[j];

                    tg_zzzzz_xxyzzz_0[j] = pb_z * tg_zzzz_xxyzzz_0[j] + fr * tg_zzzz_xxyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxyzzz_0[j] - tg_zzz_xxyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xxyzz_1[j];

                    tg_zzzzz_xxzzzz_0[j] = pb_z * tg_zzzz_xxzzzz_0[j] + fr * tg_zzzz_xxzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xxzzzz_0[j] - tg_zzz_xxzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_xxzzz_1[j];

                    tg_zzzzz_xyyyyy_0[j] = pb_z * tg_zzzz_xyyyyy_0[j] + fr * tg_zzzz_xyyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyyy_0[j] - tg_zzz_xyyyyy_1[j] * fl1_fza);

                    tg_zzzzz_xyyyyz_0[j] = pb_z * tg_zzzz_xyyyyz_0[j] + fr * tg_zzzz_xyyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyyz_0[j] - tg_zzz_xyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_xyyyy_1[j];

                    tg_zzzzz_xyyyzz_0[j] = pb_z * tg_zzzz_xyyyzz_0[j] + fr * tg_zzzz_xyyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyyzz_0[j] - tg_zzz_xyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_xyyyz_1[j];

                    tg_zzzzz_xyyzzz_0[j] = pb_z * tg_zzzz_xyyzzz_0[j] + fr * tg_zzzz_xyyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyyzzz_0[j] - tg_zzz_xyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_xyyzz_1[j];

                    tg_zzzzz_xyzzzz_0[j] = pb_z * tg_zzzz_xyzzzz_0[j] + fr * tg_zzzz_xyzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xyzzzz_0[j] - tg_zzz_xyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_xyzzz_1[j];

                    tg_zzzzz_xzzzzz_0[j] = pb_z * tg_zzzz_xzzzzz_0[j] + fr * tg_zzzz_xzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_xzzzzz_0[j] - tg_zzz_xzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzz_xzzzz_1[j];

                    tg_zzzzz_yyyyyy_0[j] = pb_z * tg_zzzz_yyyyyy_0[j] + fr * tg_zzzz_yyyyyy_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyyy_0[j] - tg_zzz_yyyyyy_1[j] * fl1_fza);

                    tg_zzzzz_yyyyyz_0[j] = pb_z * tg_zzzz_yyyyyz_0[j] + fr * tg_zzzz_yyyyyz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyyz_0[j] - tg_zzz_yyyyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzz_yyyyy_1[j];

                    tg_zzzzz_yyyyzz_0[j] = pb_z * tg_zzzz_yyyyzz_0[j] + fr * tg_zzzz_yyyyzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyyzz_0[j] - tg_zzz_yyyyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzz_yyyyz_1[j];

                    tg_zzzzz_yyyzzz_0[j] = pb_z * tg_zzzz_yyyzzz_0[j] + fr * tg_zzzz_yyyzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyyzzz_0[j] - tg_zzz_yyyzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzz_yyyzz_1[j];

                    tg_zzzzz_yyzzzz_0[j] = pb_z * tg_zzzz_yyzzzz_0[j] + fr * tg_zzzz_yyzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yyzzzz_0[j] - tg_zzz_yyzzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzz_yyzzz_1[j];

                    tg_zzzzz_yzzzzz_0[j] = pb_z * tg_zzzz_yzzzzz_0[j] + fr * tg_zzzz_yzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_yzzzzz_0[j] - tg_zzz_yzzzzz_1[j] * fl1_fza) + 2.5 * fl1_fxn * tg_zzzz_yzzzz_1[j];

                    tg_zzzzz_zzzzzz_0[j] = pb_z * tg_zzzz_zzzzzz_0[j] + fr * tg_zzzz_zzzzzz_1[j] + 2.0 * fl1_fx * (tg_zzz_zzzzzz_0[j] - tg_zzz_zzzzzz_1[j] * fl1_fza) + 3.0 * fl1_fxn * tg_zzzz_zzzzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

