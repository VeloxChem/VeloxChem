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

#include "ElectronRepulsionRecFuncForKG.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSKSG(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSKSG_0_90(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSKSG_90_180(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSKSG_180_270(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSG_270_360(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSG_360_450(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSG_450_540(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSKSG_0_90(      CMemBlock2D<double>* primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,90)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxxxx_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx); 

                auto tg_xxxxxx_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 1); 

                auto tg_xxxxxx_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 2); 

                auto tg_xxxxxx_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 3); 

                auto tg_xxxxxx_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 4); 

                auto tg_xxxxxx_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 5); 

                auto tg_xxxxxx_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 6); 

                auto tg_xxxxxx_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 7); 

                auto tg_xxxxxx_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 8); 

                auto tg_xxxxxx_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 9); 

                auto tg_xxxxxx_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 10); 

                auto tg_xxxxxx_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 11); 

                auto tg_xxxxxx_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 12); 

                auto tg_xxxxxx_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 13); 

                auto tg_xxxxxx_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 14); 

                auto tg_xxxxxy_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 15); 

                auto tg_xxxxxy_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 16); 

                auto tg_xxxxxy_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 17); 

                auto tg_xxxxxy_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 18); 

                auto tg_xxxxxy_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 19); 

                auto tg_xxxxxy_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 20); 

                auto tg_xxxxxy_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 21); 

                auto tg_xxxxxy_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 22); 

                auto tg_xxxxxy_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 23); 

                auto tg_xxxxxy_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 24); 

                auto tg_xxxxxy_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 25); 

                auto tg_xxxxxy_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 26); 

                auto tg_xxxxxy_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 27); 

                auto tg_xxxxxy_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 28); 

                auto tg_xxxxxy_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 29); 

                auto tg_xxxxxz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 30); 

                auto tg_xxxxxz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 31); 

                auto tg_xxxxxz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 32); 

                auto tg_xxxxxz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 33); 

                auto tg_xxxxxz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 34); 

                auto tg_xxxxxz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 35); 

                auto tg_xxxxxz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 36); 

                auto tg_xxxxxz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 37); 

                auto tg_xxxxxz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 38); 

                auto tg_xxxxxz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 39); 

                auto tg_xxxxxz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 40); 

                auto tg_xxxxxz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 41); 

                auto tg_xxxxxz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 42); 

                auto tg_xxxxxz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 43); 

                auto tg_xxxxxz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 44); 

                auto tg_xxxxyy_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 45); 

                auto tg_xxxxyy_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 46); 

                auto tg_xxxxyy_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 47); 

                auto tg_xxxxyy_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 48); 

                auto tg_xxxxyy_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 49); 

                auto tg_xxxxyy_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 50); 

                auto tg_xxxxyy_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 51); 

                auto tg_xxxxyy_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 52); 

                auto tg_xxxxyy_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 53); 

                auto tg_xxxxyy_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 54); 

                auto tg_xxxxyy_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 55); 

                auto tg_xxxxyy_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 56); 

                auto tg_xxxxyy_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 57); 

                auto tg_xxxxyy_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 58); 

                auto tg_xxxxyy_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 59); 

                auto tg_xxxxyz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 60); 

                auto tg_xxxxyz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 61); 

                auto tg_xxxxyz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 62); 

                auto tg_xxxxyz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 63); 

                auto tg_xxxxyz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 64); 

                auto tg_xxxxyz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 65); 

                auto tg_xxxxyz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 66); 

                auto tg_xxxxyz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 67); 

                auto tg_xxxxyz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 68); 

                auto tg_xxxxyz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 69); 

                auto tg_xxxxyz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 70); 

                auto tg_xxxxyz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 71); 

                auto tg_xxxxyz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 72); 

                auto tg_xxxxyz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 73); 

                auto tg_xxxxyz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 74); 

                auto tg_xxxxzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 75); 

                auto tg_xxxxzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 76); 

                auto tg_xxxxzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 77); 

                auto tg_xxxxzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 78); 

                auto tg_xxxxzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 79); 

                auto tg_xxxxzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 80); 

                auto tg_xxxxzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 81); 

                auto tg_xxxxzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 82); 

                auto tg_xxxxzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 83); 

                auto tg_xxxxzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 84); 

                auto tg_xxxxzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 85); 

                auto tg_xxxxzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 86); 

                auto tg_xxxxzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 87); 

                auto tg_xxxxzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 88); 

                auto tg_xxxxzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 89); 

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

                auto tg_xxxxx_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx); 

                auto tg_xxxxx_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 1); 

                auto tg_xxxxx_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 2); 

                auto tg_xxxxx_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 3); 

                auto tg_xxxxx_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 4); 

                auto tg_xxxxx_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 5); 

                auto tg_xxxxx_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 6); 

                auto tg_xxxxx_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 7); 

                auto tg_xxxxx_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 8); 

                auto tg_xxxxx_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 9); 

                auto tg_xxxxx_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 10); 

                auto tg_xxxxx_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 11); 

                auto tg_xxxxx_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 12); 

                auto tg_xxxxx_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 13); 

                auto tg_xxxxx_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 14); 

                auto tg_xxxxy_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 15); 

                auto tg_xxxxy_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 16); 

                auto tg_xxxxy_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 17); 

                auto tg_xxxxy_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 18); 

                auto tg_xxxxy_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 19); 

                auto tg_xxxxy_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 20); 

                auto tg_xxxxy_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 21); 

                auto tg_xxxxy_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 22); 

                auto tg_xxxxy_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 23); 

                auto tg_xxxxy_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 24); 

                auto tg_xxxxy_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 25); 

                auto tg_xxxxy_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 26); 

                auto tg_xxxxy_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 27); 

                auto tg_xxxxy_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 28); 

                auto tg_xxxxy_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 29); 

                auto tg_xxxxz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 30); 

                auto tg_xxxxz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 31); 

                auto tg_xxxxz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 32); 

                auto tg_xxxxz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 33); 

                auto tg_xxxxz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 34); 

                auto tg_xxxxz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 35); 

                auto tg_xxxxz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 36); 

                auto tg_xxxxz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 37); 

                auto tg_xxxxz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 38); 

                auto tg_xxxxz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 39); 

                auto tg_xxxxz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 40); 

                auto tg_xxxxz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 41); 

                auto tg_xxxxz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 42); 

                auto tg_xxxxz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 43); 

                auto tg_xxxxz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 44); 

                auto tg_xxxyy_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 45); 

                auto tg_xxxyy_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 46); 

                auto tg_xxxyy_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 47); 

                auto tg_xxxyy_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 48); 

                auto tg_xxxyy_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 49); 

                auto tg_xxxyy_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 50); 

                auto tg_xxxyy_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 51); 

                auto tg_xxxyy_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 52); 

                auto tg_xxxyy_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 53); 

                auto tg_xxxyy_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 54); 

                auto tg_xxxyy_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 55); 

                auto tg_xxxyy_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 56); 

                auto tg_xxxyy_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 57); 

                auto tg_xxxyy_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 58); 

                auto tg_xxxyy_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 59); 

                auto tg_xxxyz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 60); 

                auto tg_xxxyz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 61); 

                auto tg_xxxyz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 62); 

                auto tg_xxxyz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 63); 

                auto tg_xxxyz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 64); 

                auto tg_xxxyz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 65); 

                auto tg_xxxyz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 66); 

                auto tg_xxxyz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 67); 

                auto tg_xxxyz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 68); 

                auto tg_xxxyz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 69); 

                auto tg_xxxyz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 70); 

                auto tg_xxxyz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 71); 

                auto tg_xxxyz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 72); 

                auto tg_xxxyz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 73); 

                auto tg_xxxyz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 74); 

                auto tg_xxxzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 75); 

                auto tg_xxxzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 76); 

                auto tg_xxxzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 77); 

                auto tg_xxxzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 78); 

                auto tg_xxxzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 79); 

                auto tg_xxxzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 80); 

                auto tg_xxxzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 81); 

                auto tg_xxxzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 82); 

                auto tg_xxxzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 83); 

                auto tg_xxxzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 84); 

                auto tg_xxxzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 85); 

                auto tg_xxxzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 86); 

                auto tg_xxxzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 87); 

                auto tg_xxxzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 88); 

                auto tg_xxxzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 89); 

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

                auto tg_xxxxxx_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx); 

                auto tg_xxxxxx_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 1); 

                auto tg_xxxxxx_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 2); 

                auto tg_xxxxxx_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 3); 

                auto tg_xxxxxx_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 4); 

                auto tg_xxxxxx_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 5); 

                auto tg_xxxxxx_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 6); 

                auto tg_xxxxxx_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 7); 

                auto tg_xxxxxx_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 8); 

                auto tg_xxxxxx_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 9); 

                auto tg_xxxxxy_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 10); 

                auto tg_xxxxxy_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 11); 

                auto tg_xxxxxy_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 12); 

                auto tg_xxxxxy_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 13); 

                auto tg_xxxxxy_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 14); 

                auto tg_xxxxxy_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 15); 

                auto tg_xxxxxy_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 16); 

                auto tg_xxxxxy_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 17); 

                auto tg_xxxxxy_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 18); 

                auto tg_xxxxxy_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 19); 

                auto tg_xxxxxz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 20); 

                auto tg_xxxxxz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 21); 

                auto tg_xxxxxz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 22); 

                auto tg_xxxxxz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 23); 

                auto tg_xxxxxz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 24); 

                auto tg_xxxxxz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 25); 

                auto tg_xxxxxz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 26); 

                auto tg_xxxxxz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 27); 

                auto tg_xxxxxz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 28); 

                auto tg_xxxxxz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 29); 

                auto tg_xxxxyy_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 30); 

                auto tg_xxxxyy_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 31); 

                auto tg_xxxxyy_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 32); 

                auto tg_xxxxyy_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 33); 

                auto tg_xxxxyy_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 34); 

                auto tg_xxxxyy_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 35); 

                auto tg_xxxxyy_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 36); 

                auto tg_xxxxyy_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 37); 

                auto tg_xxxxyy_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 38); 

                auto tg_xxxxyy_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 39); 

                auto tg_xxxxyz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 40); 

                auto tg_xxxxyz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 41); 

                auto tg_xxxxyz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 42); 

                auto tg_xxxxyz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 43); 

                auto tg_xxxxyz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 44); 

                auto tg_xxxxyz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 45); 

                auto tg_xxxxyz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 46); 

                auto tg_xxxxyz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 47); 

                auto tg_xxxxyz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 48); 

                auto tg_xxxxyz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 49); 

                auto tg_xxxxzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 50); 

                auto tg_xxxxzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 51); 

                auto tg_xxxxzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 52); 

                auto tg_xxxxzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 53); 

                auto tg_xxxxzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 54); 

                auto tg_xxxxzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 55); 

                auto tg_xxxxzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 56); 

                auto tg_xxxxzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 57); 

                auto tg_xxxxzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 58); 

                auto tg_xxxxzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 59); 

                // set up pointers to integrals

                auto tg_xxxxxxx_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx); 

                auto tg_xxxxxxx_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 1); 

                auto tg_xxxxxxx_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 2); 

                auto tg_xxxxxxx_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 3); 

                auto tg_xxxxxxx_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 4); 

                auto tg_xxxxxxx_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 5); 

                auto tg_xxxxxxx_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 6); 

                auto tg_xxxxxxx_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 7); 

                auto tg_xxxxxxx_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 8); 

                auto tg_xxxxxxx_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 9); 

                auto tg_xxxxxxx_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 10); 

                auto tg_xxxxxxx_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 11); 

                auto tg_xxxxxxx_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 12); 

                auto tg_xxxxxxx_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 13); 

                auto tg_xxxxxxx_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 14); 

                auto tg_xxxxxxy_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 15); 

                auto tg_xxxxxxy_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 16); 

                auto tg_xxxxxxy_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 17); 

                auto tg_xxxxxxy_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 18); 

                auto tg_xxxxxxy_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 19); 

                auto tg_xxxxxxy_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 20); 

                auto tg_xxxxxxy_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 21); 

                auto tg_xxxxxxy_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 22); 

                auto tg_xxxxxxy_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 23); 

                auto tg_xxxxxxy_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 24); 

                auto tg_xxxxxxy_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 25); 

                auto tg_xxxxxxy_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 26); 

                auto tg_xxxxxxy_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 27); 

                auto tg_xxxxxxy_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 28); 

                auto tg_xxxxxxy_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 29); 

                auto tg_xxxxxxz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 30); 

                auto tg_xxxxxxz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 31); 

                auto tg_xxxxxxz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 32); 

                auto tg_xxxxxxz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 33); 

                auto tg_xxxxxxz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 34); 

                auto tg_xxxxxxz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 35); 

                auto tg_xxxxxxz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 36); 

                auto tg_xxxxxxz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 37); 

                auto tg_xxxxxxz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 38); 

                auto tg_xxxxxxz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 39); 

                auto tg_xxxxxxz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 40); 

                auto tg_xxxxxxz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 41); 

                auto tg_xxxxxxz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 42); 

                auto tg_xxxxxxz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 43); 

                auto tg_xxxxxxz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 44); 

                auto tg_xxxxxyy_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 45); 

                auto tg_xxxxxyy_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 46); 

                auto tg_xxxxxyy_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 47); 

                auto tg_xxxxxyy_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 48); 

                auto tg_xxxxxyy_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 49); 

                auto tg_xxxxxyy_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 50); 

                auto tg_xxxxxyy_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 51); 

                auto tg_xxxxxyy_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 52); 

                auto tg_xxxxxyy_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 53); 

                auto tg_xxxxxyy_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 54); 

                auto tg_xxxxxyy_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 55); 

                auto tg_xxxxxyy_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 56); 

                auto tg_xxxxxyy_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 57); 

                auto tg_xxxxxyy_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 58); 

                auto tg_xxxxxyy_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 59); 

                auto tg_xxxxxyz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 60); 

                auto tg_xxxxxyz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 61); 

                auto tg_xxxxxyz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 62); 

                auto tg_xxxxxyz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 63); 

                auto tg_xxxxxyz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 64); 

                auto tg_xxxxxyz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 65); 

                auto tg_xxxxxyz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 66); 

                auto tg_xxxxxyz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 67); 

                auto tg_xxxxxyz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 68); 

                auto tg_xxxxxyz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 69); 

                auto tg_xxxxxyz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 70); 

                auto tg_xxxxxyz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 71); 

                auto tg_xxxxxyz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 72); 

                auto tg_xxxxxyz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 73); 

                auto tg_xxxxxyz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 74); 

                auto tg_xxxxxzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 75); 

                auto tg_xxxxxzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 76); 

                auto tg_xxxxxzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 77); 

                auto tg_xxxxxzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 78); 

                auto tg_xxxxxzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 79); 

                auto tg_xxxxxzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 80); 

                auto tg_xxxxxzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 81); 

                auto tg_xxxxxzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 82); 

                auto tg_xxxxxzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 83); 

                auto tg_xxxxxzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 84); 

                auto tg_xxxxxzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 85); 

                auto tg_xxxxxzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 86); 

                auto tg_xxxxxzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 87); 

                auto tg_xxxxxzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 88); 

                auto tg_xxxxxzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 89); 

                // Batch of Integrals (0,90)

                #pragma omp simd aligned(fxn, fza, tg_xxxxx_xxxx_0, tg_xxxxx_xxxx_1, tg_xxxxx_xxxy_0, \
                                         tg_xxxxx_xxxy_1, tg_xxxxx_xxxz_0, tg_xxxxx_xxxz_1, tg_xxxxx_xxyy_0, tg_xxxxx_xxyy_1, \
                                         tg_xxxxx_xxyz_0, tg_xxxxx_xxyz_1, tg_xxxxx_xxzz_0, tg_xxxxx_xxzz_1, tg_xxxxx_xyyy_0, \
                                         tg_xxxxx_xyyy_1, tg_xxxxx_xyyz_0, tg_xxxxx_xyyz_1, tg_xxxxx_xyzz_0, tg_xxxxx_xyzz_1, \
                                         tg_xxxxx_xzzz_0, tg_xxxxx_xzzz_1, tg_xxxxx_yyyy_0, tg_xxxxx_yyyy_1, tg_xxxxx_yyyz_0, \
                                         tg_xxxxx_yyyz_1, tg_xxxxx_yyzz_0, tg_xxxxx_yyzz_1, tg_xxxxx_yzzz_0, tg_xxxxx_yzzz_1, \
                                         tg_xxxxx_zzzz_0, tg_xxxxx_zzzz_1, tg_xxxxxx_xxx_1, tg_xxxxxx_xxxx_0, \
                                         tg_xxxxxx_xxxx_1, tg_xxxxxx_xxxy_0, tg_xxxxxx_xxxy_1, tg_xxxxxx_xxxz_0, \
                                         tg_xxxxxx_xxxz_1, tg_xxxxxx_xxy_1, tg_xxxxxx_xxyy_0, tg_xxxxxx_xxyy_1, \
                                         tg_xxxxxx_xxyz_0, tg_xxxxxx_xxyz_1, tg_xxxxxx_xxz_1, tg_xxxxxx_xxzz_0, \
                                         tg_xxxxxx_xxzz_1, tg_xxxxxx_xyy_1, tg_xxxxxx_xyyy_0, tg_xxxxxx_xyyy_1, \
                                         tg_xxxxxx_xyyz_0, tg_xxxxxx_xyyz_1, tg_xxxxxx_xyz_1, tg_xxxxxx_xyzz_0, \
                                         tg_xxxxxx_xyzz_1, tg_xxxxxx_xzz_1, tg_xxxxxx_xzzz_0, tg_xxxxxx_xzzz_1, \
                                         tg_xxxxxx_yyy_1, tg_xxxxxx_yyyy_0, tg_xxxxxx_yyyy_1, tg_xxxxxx_yyyz_0, \
                                         tg_xxxxxx_yyyz_1, tg_xxxxxx_yyz_1, tg_xxxxxx_yyzz_0, tg_xxxxxx_yyzz_1, \
                                         tg_xxxxxx_yzz_1, tg_xxxxxx_yzzz_0, tg_xxxxxx_yzzz_1, tg_xxxxxx_zzz_1, \
                                         tg_xxxxxx_zzzz_0, tg_xxxxxx_zzzz_1, tg_xxxxxxx_xxxx_0, tg_xxxxxxx_xxxy_0, \
                                         tg_xxxxxxx_xxxz_0, tg_xxxxxxx_xxyy_0, tg_xxxxxxx_xxyz_0, tg_xxxxxxx_xxzz_0, \
                                         tg_xxxxxxx_xyyy_0, tg_xxxxxxx_xyyz_0, tg_xxxxxxx_xyzz_0, tg_xxxxxxx_xzzz_0, \
                                         tg_xxxxxxx_yyyy_0, tg_xxxxxxx_yyyz_0, tg_xxxxxxx_yyzz_0, tg_xxxxxxx_yzzz_0, \
                                         tg_xxxxxxx_zzzz_0, tg_xxxxxxy_xxxx_0, tg_xxxxxxy_xxxy_0, tg_xxxxxxy_xxxz_0, \
                                         tg_xxxxxxy_xxyy_0, tg_xxxxxxy_xxyz_0, tg_xxxxxxy_xxzz_0, tg_xxxxxxy_xyyy_0, \
                                         tg_xxxxxxy_xyyz_0, tg_xxxxxxy_xyzz_0, tg_xxxxxxy_xzzz_0, tg_xxxxxxy_yyyy_0, \
                                         tg_xxxxxxy_yyyz_0, tg_xxxxxxy_yyzz_0, tg_xxxxxxy_yzzz_0, tg_xxxxxxy_zzzz_0, \
                                         tg_xxxxxxz_xxxx_0, tg_xxxxxxz_xxxy_0, tg_xxxxxxz_xxxz_0, tg_xxxxxxz_xxyy_0, \
                                         tg_xxxxxxz_xxyz_0, tg_xxxxxxz_xxzz_0, tg_xxxxxxz_xyyy_0, tg_xxxxxxz_xyyz_0, \
                                         tg_xxxxxxz_xyzz_0, tg_xxxxxxz_xzzz_0, tg_xxxxxxz_yyyy_0, tg_xxxxxxz_yyyz_0, \
                                         tg_xxxxxxz_yyzz_0, tg_xxxxxxz_yzzz_0, tg_xxxxxxz_zzzz_0, tg_xxxxxy_xxx_1, \
                                         tg_xxxxxy_xxxx_0, tg_xxxxxy_xxxx_1, tg_xxxxxy_xxxy_0, tg_xxxxxy_xxxy_1, \
                                         tg_xxxxxy_xxxz_0, tg_xxxxxy_xxxz_1, tg_xxxxxy_xxy_1, tg_xxxxxy_xxyy_0, \
                                         tg_xxxxxy_xxyy_1, tg_xxxxxy_xxyz_0, tg_xxxxxy_xxyz_1, tg_xxxxxy_xxz_1, \
                                         tg_xxxxxy_xxzz_0, tg_xxxxxy_xxzz_1, tg_xxxxxy_xyy_1, tg_xxxxxy_xyyy_0, \
                                         tg_xxxxxy_xyyy_1, tg_xxxxxy_xyyz_0, tg_xxxxxy_xyyz_1, tg_xxxxxy_xyz_1, \
                                         tg_xxxxxy_xyzz_0, tg_xxxxxy_xyzz_1, tg_xxxxxy_xzz_1, tg_xxxxxy_xzzz_0, \
                                         tg_xxxxxy_xzzz_1, tg_xxxxxy_yyy_1, tg_xxxxxy_yyyy_0, tg_xxxxxy_yyyy_1, \
                                         tg_xxxxxy_yyyz_0, tg_xxxxxy_yyyz_1, tg_xxxxxy_yyz_1, tg_xxxxxy_yyzz_0, \
                                         tg_xxxxxy_yyzz_1, tg_xxxxxy_yzz_1, tg_xxxxxy_yzzz_0, tg_xxxxxy_yzzz_1, \
                                         tg_xxxxxy_zzz_1, tg_xxxxxy_zzzz_0, tg_xxxxxy_zzzz_1, tg_xxxxxyy_xxxx_0, \
                                         tg_xxxxxyy_xxxy_0, tg_xxxxxyy_xxxz_0, tg_xxxxxyy_xxyy_0, tg_xxxxxyy_xxyz_0, \
                                         tg_xxxxxyy_xxzz_0, tg_xxxxxyy_xyyy_0, tg_xxxxxyy_xyyz_0, tg_xxxxxyy_xyzz_0, \
                                         tg_xxxxxyy_xzzz_0, tg_xxxxxyy_yyyy_0, tg_xxxxxyy_yyyz_0, tg_xxxxxyy_yyzz_0, \
                                         tg_xxxxxyy_yzzz_0, tg_xxxxxyy_zzzz_0, tg_xxxxxyz_xxxx_0, tg_xxxxxyz_xxxy_0, \
                                         tg_xxxxxyz_xxxz_0, tg_xxxxxyz_xxyy_0, tg_xxxxxyz_xxyz_0, tg_xxxxxyz_xxzz_0, \
                                         tg_xxxxxyz_xyyy_0, tg_xxxxxyz_xyyz_0, tg_xxxxxyz_xyzz_0, tg_xxxxxyz_xzzz_0, \
                                         tg_xxxxxyz_yyyy_0, tg_xxxxxyz_yyyz_0, tg_xxxxxyz_yyzz_0, tg_xxxxxyz_yzzz_0, \
                                         tg_xxxxxyz_zzzz_0, tg_xxxxxz_xxx_1, tg_xxxxxz_xxxx_0, tg_xxxxxz_xxxx_1, \
                                         tg_xxxxxz_xxxy_0, tg_xxxxxz_xxxy_1, tg_xxxxxz_xxxz_0, tg_xxxxxz_xxxz_1, \
                                         tg_xxxxxz_xxy_1, tg_xxxxxz_xxyy_0, tg_xxxxxz_xxyy_1, tg_xxxxxz_xxyz_0, \
                                         tg_xxxxxz_xxyz_1, tg_xxxxxz_xxz_1, tg_xxxxxz_xxzz_0, tg_xxxxxz_xxzz_1, \
                                         tg_xxxxxz_xyy_1, tg_xxxxxz_xyyy_0, tg_xxxxxz_xyyy_1, tg_xxxxxz_xyyz_0, \
                                         tg_xxxxxz_xyyz_1, tg_xxxxxz_xyz_1, tg_xxxxxz_xyzz_0, tg_xxxxxz_xyzz_1, \
                                         tg_xxxxxz_xzz_1, tg_xxxxxz_xzzz_0, tg_xxxxxz_xzzz_1, tg_xxxxxz_yyy_1, \
                                         tg_xxxxxz_yyyy_0, tg_xxxxxz_yyyy_1, tg_xxxxxz_yyyz_0, tg_xxxxxz_yyyz_1, \
                                         tg_xxxxxz_yyz_1, tg_xxxxxz_yyzz_0, tg_xxxxxz_yyzz_1, tg_xxxxxz_yzz_1, \
                                         tg_xxxxxz_yzzz_0, tg_xxxxxz_yzzz_1, tg_xxxxxz_zzz_1, tg_xxxxxz_zzzz_0, \
                                         tg_xxxxxz_zzzz_1, tg_xxxxxzz_xxxx_0, tg_xxxxxzz_xxxy_0, tg_xxxxxzz_xxxz_0, \
                                         tg_xxxxxzz_xxyy_0, tg_xxxxxzz_xxyz_0, tg_xxxxxzz_xxzz_0, tg_xxxxxzz_xyyy_0, \
                                         tg_xxxxxzz_xyyz_0, tg_xxxxxzz_xyzz_0, tg_xxxxxzz_xzzz_0, tg_xxxxxzz_yyyy_0, \
                                         tg_xxxxxzz_yyyz_0, tg_xxxxxzz_yyzz_0, tg_xxxxxzz_yzzz_0, tg_xxxxxzz_zzzz_0, \
                                         tg_xxxxy_xxxx_0, tg_xxxxy_xxxx_1, tg_xxxxy_xxxy_0, tg_xxxxy_xxxy_1, tg_xxxxy_xxxz_0, \
                                         tg_xxxxy_xxxz_1, tg_xxxxy_xxyy_0, tg_xxxxy_xxyy_1, tg_xxxxy_xxyz_0, tg_xxxxy_xxyz_1, \
                                         tg_xxxxy_xxzz_0, tg_xxxxy_xxzz_1, tg_xxxxy_xyyy_0, tg_xxxxy_xyyy_1, tg_xxxxy_xyyz_0, \
                                         tg_xxxxy_xyyz_1, tg_xxxxy_xyzz_0, tg_xxxxy_xyzz_1, tg_xxxxy_xzzz_0, tg_xxxxy_xzzz_1, \
                                         tg_xxxxy_yyyy_0, tg_xxxxy_yyyy_1, tg_xxxxy_yyyz_0, tg_xxxxy_yyyz_1, tg_xxxxy_yyzz_0, \
                                         tg_xxxxy_yyzz_1, tg_xxxxy_yzzz_0, tg_xxxxy_yzzz_1, tg_xxxxy_zzzz_0, tg_xxxxy_zzzz_1, \
                                         tg_xxxxyy_xxx_1, tg_xxxxyy_xxxx_0, tg_xxxxyy_xxxx_1, tg_xxxxyy_xxxy_0, \
                                         tg_xxxxyy_xxxy_1, tg_xxxxyy_xxxz_0, tg_xxxxyy_xxxz_1, tg_xxxxyy_xxy_1, \
                                         tg_xxxxyy_xxyy_0, tg_xxxxyy_xxyy_1, tg_xxxxyy_xxyz_0, tg_xxxxyy_xxyz_1, \
                                         tg_xxxxyy_xxz_1, tg_xxxxyy_xxzz_0, tg_xxxxyy_xxzz_1, tg_xxxxyy_xyy_1, \
                                         tg_xxxxyy_xyyy_0, tg_xxxxyy_xyyy_1, tg_xxxxyy_xyyz_0, tg_xxxxyy_xyyz_1, \
                                         tg_xxxxyy_xyz_1, tg_xxxxyy_xyzz_0, tg_xxxxyy_xyzz_1, tg_xxxxyy_xzz_1, \
                                         tg_xxxxyy_xzzz_0, tg_xxxxyy_xzzz_1, tg_xxxxyy_yyy_1, tg_xxxxyy_yyyy_0, \
                                         tg_xxxxyy_yyyy_1, tg_xxxxyy_yyyz_0, tg_xxxxyy_yyyz_1, tg_xxxxyy_yyz_1, \
                                         tg_xxxxyy_yyzz_0, tg_xxxxyy_yyzz_1, tg_xxxxyy_yzz_1, tg_xxxxyy_yzzz_0, \
                                         tg_xxxxyy_yzzz_1, tg_xxxxyy_zzz_1, tg_xxxxyy_zzzz_0, tg_xxxxyy_zzzz_1, \
                                         tg_xxxxyz_xxx_1, tg_xxxxyz_xxxx_0, tg_xxxxyz_xxxx_1, tg_xxxxyz_xxxy_0, \
                                         tg_xxxxyz_xxxy_1, tg_xxxxyz_xxxz_0, tg_xxxxyz_xxxz_1, tg_xxxxyz_xxy_1, \
                                         tg_xxxxyz_xxyy_0, tg_xxxxyz_xxyy_1, tg_xxxxyz_xxyz_0, tg_xxxxyz_xxyz_1, \
                                         tg_xxxxyz_xxz_1, tg_xxxxyz_xxzz_0, tg_xxxxyz_xxzz_1, tg_xxxxyz_xyy_1, \
                                         tg_xxxxyz_xyyy_0, tg_xxxxyz_xyyy_1, tg_xxxxyz_xyyz_0, tg_xxxxyz_xyyz_1, \
                                         tg_xxxxyz_xyz_1, tg_xxxxyz_xyzz_0, tg_xxxxyz_xyzz_1, tg_xxxxyz_xzz_1, \
                                         tg_xxxxyz_xzzz_0, tg_xxxxyz_xzzz_1, tg_xxxxyz_yyy_1, tg_xxxxyz_yyyy_0, \
                                         tg_xxxxyz_yyyy_1, tg_xxxxyz_yyyz_0, tg_xxxxyz_yyyz_1, tg_xxxxyz_yyz_1, \
                                         tg_xxxxyz_yyzz_0, tg_xxxxyz_yyzz_1, tg_xxxxyz_yzz_1, tg_xxxxyz_yzzz_0, \
                                         tg_xxxxyz_yzzz_1, tg_xxxxyz_zzz_1, tg_xxxxyz_zzzz_0, tg_xxxxyz_zzzz_1, \
                                         tg_xxxxz_xxxx_0, tg_xxxxz_xxxx_1, tg_xxxxz_xxxy_0, tg_xxxxz_xxxy_1, tg_xxxxz_xxxz_0, \
                                         tg_xxxxz_xxxz_1, tg_xxxxz_xxyy_0, tg_xxxxz_xxyy_1, tg_xxxxz_xxyz_0, tg_xxxxz_xxyz_1, \
                                         tg_xxxxz_xxzz_0, tg_xxxxz_xxzz_1, tg_xxxxz_xyyy_0, tg_xxxxz_xyyy_1, tg_xxxxz_xyyz_0, \
                                         tg_xxxxz_xyyz_1, tg_xxxxz_xyzz_0, tg_xxxxz_xyzz_1, tg_xxxxz_xzzz_0, tg_xxxxz_xzzz_1, \
                                         tg_xxxxz_yyyy_0, tg_xxxxz_yyyy_1, tg_xxxxz_yyyz_0, tg_xxxxz_yyyz_1, tg_xxxxz_yyzz_0, \
                                         tg_xxxxz_yyzz_1, tg_xxxxz_yzzz_0, tg_xxxxz_yzzz_1, tg_xxxxz_zzzz_0, tg_xxxxz_zzzz_1, \
                                         tg_xxxxzz_xxx_1, tg_xxxxzz_xxxx_0, tg_xxxxzz_xxxx_1, tg_xxxxzz_xxxy_0, \
                                         tg_xxxxzz_xxxy_1, tg_xxxxzz_xxxz_0, tg_xxxxzz_xxxz_1, tg_xxxxzz_xxy_1, \
                                         tg_xxxxzz_xxyy_0, tg_xxxxzz_xxyy_1, tg_xxxxzz_xxyz_0, tg_xxxxzz_xxyz_1, \
                                         tg_xxxxzz_xxz_1, tg_xxxxzz_xxzz_0, tg_xxxxzz_xxzz_1, tg_xxxxzz_xyy_1, \
                                         tg_xxxxzz_xyyy_0, tg_xxxxzz_xyyy_1, tg_xxxxzz_xyyz_0, tg_xxxxzz_xyyz_1, \
                                         tg_xxxxzz_xyz_1, tg_xxxxzz_xyzz_0, tg_xxxxzz_xyzz_1, tg_xxxxzz_xzz_1, \
                                         tg_xxxxzz_xzzz_0, tg_xxxxzz_xzzz_1, tg_xxxxzz_yyy_1, tg_xxxxzz_yyyy_0, \
                                         tg_xxxxzz_yyyy_1, tg_xxxxzz_yyyz_0, tg_xxxxzz_yyyz_1, tg_xxxxzz_yyz_1, \
                                         tg_xxxxzz_yyzz_0, tg_xxxxzz_yyzz_1, tg_xxxxzz_yzz_1, tg_xxxxzz_yzzz_0, \
                                         tg_xxxxzz_yzzz_1, tg_xxxxzz_zzz_1, tg_xxxxzz_zzzz_0, tg_xxxxzz_zzzz_1, \
                                         tg_xxxyy_xxxx_0, tg_xxxyy_xxxx_1, tg_xxxyy_xxxy_0, tg_xxxyy_xxxy_1, tg_xxxyy_xxxz_0, \
                                         tg_xxxyy_xxxz_1, tg_xxxyy_xxyy_0, tg_xxxyy_xxyy_1, tg_xxxyy_xxyz_0, tg_xxxyy_xxyz_1, \
                                         tg_xxxyy_xxzz_0, tg_xxxyy_xxzz_1, tg_xxxyy_xyyy_0, tg_xxxyy_xyyy_1, tg_xxxyy_xyyz_0, \
                                         tg_xxxyy_xyyz_1, tg_xxxyy_xyzz_0, tg_xxxyy_xyzz_1, tg_xxxyy_xzzz_0, tg_xxxyy_xzzz_1, \
                                         tg_xxxyy_yyyy_0, tg_xxxyy_yyyy_1, tg_xxxyy_yyyz_0, tg_xxxyy_yyyz_1, tg_xxxyy_yyzz_0, \
                                         tg_xxxyy_yyzz_1, tg_xxxyy_yzzz_0, tg_xxxyy_yzzz_1, tg_xxxyy_zzzz_0, tg_xxxyy_zzzz_1, \
                                         tg_xxxyz_xxxx_0, tg_xxxyz_xxxx_1, tg_xxxyz_xxxy_0, tg_xxxyz_xxxy_1, tg_xxxyz_xxxz_0, \
                                         tg_xxxyz_xxxz_1, tg_xxxyz_xxyy_0, tg_xxxyz_xxyy_1, tg_xxxyz_xxyz_0, tg_xxxyz_xxyz_1, \
                                         tg_xxxyz_xxzz_0, tg_xxxyz_xxzz_1, tg_xxxyz_xyyy_0, tg_xxxyz_xyyy_1, tg_xxxyz_xyyz_0, \
                                         tg_xxxyz_xyyz_1, tg_xxxyz_xyzz_0, tg_xxxyz_xyzz_1, tg_xxxyz_xzzz_0, tg_xxxyz_xzzz_1, \
                                         tg_xxxyz_yyyy_0, tg_xxxyz_yyyy_1, tg_xxxyz_yyyz_0, tg_xxxyz_yyyz_1, tg_xxxyz_yyzz_0, \
                                         tg_xxxyz_yyzz_1, tg_xxxyz_yzzz_0, tg_xxxyz_yzzz_1, tg_xxxyz_zzzz_0, tg_xxxyz_zzzz_1, \
                                         tg_xxxzz_xxxx_0, tg_xxxzz_xxxx_1, tg_xxxzz_xxxy_0, tg_xxxzz_xxxy_1, tg_xxxzz_xxxz_0, \
                                         tg_xxxzz_xxxz_1, tg_xxxzz_xxyy_0, tg_xxxzz_xxyy_1, tg_xxxzz_xxyz_0, tg_xxxzz_xxyz_1, \
                                         tg_xxxzz_xxzz_0, tg_xxxzz_xxzz_1, tg_xxxzz_xyyy_0, tg_xxxzz_xyyy_1, tg_xxxzz_xyyz_0, \
                                         tg_xxxzz_xyyz_1, tg_xxxzz_xyzz_0, tg_xxxzz_xyzz_1, tg_xxxzz_xzzz_0, tg_xxxzz_xzzz_1, \
                                         tg_xxxzz_yyyy_0, tg_xxxzz_yyyy_1, tg_xxxzz_yyyz_0, tg_xxxzz_yyyz_1, tg_xxxzz_yyzz_0, \
                                         tg_xxxzz_yyzz_1, tg_xxxzz_yzzz_0, tg_xxxzz_yzzz_1, tg_xxxzz_zzzz_0, tg_xxxzz_zzzz_1, \
                                         wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxxx_xxxx_0[j] = pb_x * tg_xxxxxx_xxxx_0[j] + fr * tg_xxxxxx_xxxx_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxx_0[j] - tg_xxxxx_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxx_xxx_1[j];

                    tg_xxxxxxx_xxxy_0[j] = pb_x * tg_xxxxxx_xxxy_0[j] + fr * tg_xxxxxx_xxxy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxy_0[j] - tg_xxxxx_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxy_1[j];

                    tg_xxxxxxx_xxxz_0[j] = pb_x * tg_xxxxxx_xxxz_0[j] + fr * tg_xxxxxx_xxxz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxxz_0[j] - tg_xxxxx_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xxz_1[j];

                    tg_xxxxxxx_xxyy_0[j] = pb_x * tg_xxxxxx_xxyy_0[j] + fr * tg_xxxxxx_xxyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyy_0[j] - tg_xxxxx_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyy_1[j];

                    tg_xxxxxxx_xxyz_0[j] = pb_x * tg_xxxxxx_xxyz_0[j] + fr * tg_xxxxxx_xxyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxyz_0[j] - tg_xxxxx_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xyz_1[j];

                    tg_xxxxxxx_xxzz_0[j] = pb_x * tg_xxxxxx_xxzz_0[j] + fr * tg_xxxxxx_xxzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxzz_0[j] - tg_xxxxx_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xzz_1[j];

                    tg_xxxxxxx_xyyy_0[j] = pb_x * tg_xxxxxx_xyyy_0[j] + fr * tg_xxxxxx_xyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyy_0[j] - tg_xxxxx_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyy_1[j];

                    tg_xxxxxxx_xyyz_0[j] = pb_x * tg_xxxxxx_xyyz_0[j] + fr * tg_xxxxxx_xyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyyz_0[j] - tg_xxxxx_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yyz_1[j];

                    tg_xxxxxxx_xyzz_0[j] = pb_x * tg_xxxxxx_xyzz_0[j] + fr * tg_xxxxxx_xyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyzz_0[j] - tg_xxxxx_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yzz_1[j];

                    tg_xxxxxxx_xzzz_0[j] = pb_x * tg_xxxxxx_xzzz_0[j] + fr * tg_xxxxxx_xzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xzzz_0[j] - tg_xxxxx_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_zzz_1[j];

                    tg_xxxxxxx_yyyy_0[j] = pb_x * tg_xxxxxx_yyyy_0[j] + fr * tg_xxxxxx_yyyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyy_0[j] - tg_xxxxx_yyyy_1[j] * fl1_fza);

                    tg_xxxxxxx_yyyz_0[j] = pb_x * tg_xxxxxx_yyyz_0[j] + fr * tg_xxxxxx_yyyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyyz_0[j] - tg_xxxxx_yyyz_1[j] * fl1_fza);

                    tg_xxxxxxx_yyzz_0[j] = pb_x * tg_xxxxxx_yyzz_0[j] + fr * tg_xxxxxx_yyzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyzz_0[j] - tg_xxxxx_yyzz_1[j] * fl1_fza);

                    tg_xxxxxxx_yzzz_0[j] = pb_x * tg_xxxxxx_yzzz_0[j] + fr * tg_xxxxxx_yzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yzzz_0[j] - tg_xxxxx_yzzz_1[j] * fl1_fza);

                    tg_xxxxxxx_zzzz_0[j] = pb_x * tg_xxxxxx_zzzz_0[j] + fr * tg_xxxxxx_zzzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_zzzz_0[j] - tg_xxxxx_zzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_xxxx_0[j] = pb_x * tg_xxxxxy_xxxx_0[j] + fr * tg_xxxxxy_xxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxx_0[j] - tg_xxxxy_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxy_xxx_1[j];

                    tg_xxxxxxy_xxxy_0[j] = pb_x * tg_xxxxxy_xxxy_0[j] + fr * tg_xxxxxy_xxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxy_0[j] - tg_xxxxy_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxy_1[j];

                    tg_xxxxxxy_xxxz_0[j] = pb_x * tg_xxxxxy_xxxz_0[j] + fr * tg_xxxxxy_xxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxxz_0[j] - tg_xxxxy_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xxz_1[j];

                    tg_xxxxxxy_xxyy_0[j] = pb_x * tg_xxxxxy_xxyy_0[j] + fr * tg_xxxxxy_xxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyy_0[j] - tg_xxxxy_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyy_1[j];

                    tg_xxxxxxy_xxyz_0[j] = pb_x * tg_xxxxxy_xxyz_0[j] + fr * tg_xxxxxy_xxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxyz_0[j] - tg_xxxxy_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xyz_1[j];

                    tg_xxxxxxy_xxzz_0[j] = pb_x * tg_xxxxxy_xxzz_0[j] + fr * tg_xxxxxy_xxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxzz_0[j] - tg_xxxxy_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xzz_1[j];

                    tg_xxxxxxy_xyyy_0[j] = pb_x * tg_xxxxxy_xyyy_0[j] + fr * tg_xxxxxy_xyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyy_0[j] - tg_xxxxy_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyy_1[j];

                    tg_xxxxxxy_xyyz_0[j] = pb_x * tg_xxxxxy_xyyz_0[j] + fr * tg_xxxxxy_xyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyyz_0[j] - tg_xxxxy_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yyz_1[j];

                    tg_xxxxxxy_xyzz_0[j] = pb_x * tg_xxxxxy_xyzz_0[j] + fr * tg_xxxxxy_xyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyzz_0[j] - tg_xxxxy_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yzz_1[j];

                    tg_xxxxxxy_xzzz_0[j] = pb_x * tg_xxxxxy_xzzz_0[j] + fr * tg_xxxxxy_xzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xzzz_0[j] - tg_xxxxy_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_zzz_1[j];

                    tg_xxxxxxy_yyyy_0[j] = pb_x * tg_xxxxxy_yyyy_0[j] + fr * tg_xxxxxy_yyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyy_0[j] - tg_xxxxy_yyyy_1[j] * fl1_fza);

                    tg_xxxxxxy_yyyz_0[j] = pb_x * tg_xxxxxy_yyyz_0[j] + fr * tg_xxxxxy_yyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyyz_0[j] - tg_xxxxy_yyyz_1[j] * fl1_fza);

                    tg_xxxxxxy_yyzz_0[j] = pb_x * tg_xxxxxy_yyzz_0[j] + fr * tg_xxxxxy_yyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyzz_0[j] - tg_xxxxy_yyzz_1[j] * fl1_fza);

                    tg_xxxxxxy_yzzz_0[j] = pb_x * tg_xxxxxy_yzzz_0[j] + fr * tg_xxxxxy_yzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yzzz_0[j] - tg_xxxxy_yzzz_1[j] * fl1_fza);

                    tg_xxxxxxy_zzzz_0[j] = pb_x * tg_xxxxxy_zzzz_0[j] + fr * tg_xxxxxy_zzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_zzzz_0[j] - tg_xxxxy_zzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_xxxx_0[j] = pb_x * tg_xxxxxz_xxxx_0[j] + fr * tg_xxxxxz_xxxx_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxx_0[j] - tg_xxxxz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxxz_xxx_1[j];

                    tg_xxxxxxz_xxxy_0[j] = pb_x * tg_xxxxxz_xxxy_0[j] + fr * tg_xxxxxz_xxxy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxy_0[j] - tg_xxxxz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxy_1[j];

                    tg_xxxxxxz_xxxz_0[j] = pb_x * tg_xxxxxz_xxxz_0[j] + fr * tg_xxxxxz_xxxz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxxz_0[j] - tg_xxxxz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xxz_1[j];

                    tg_xxxxxxz_xxyy_0[j] = pb_x * tg_xxxxxz_xxyy_0[j] + fr * tg_xxxxxz_xxyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyy_0[j] - tg_xxxxz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyy_1[j];

                    tg_xxxxxxz_xxyz_0[j] = pb_x * tg_xxxxxz_xxyz_0[j] + fr * tg_xxxxxz_xxyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxyz_0[j] - tg_xxxxz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xyz_1[j];

                    tg_xxxxxxz_xxzz_0[j] = pb_x * tg_xxxxxz_xxzz_0[j] + fr * tg_xxxxxz_xxzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxzz_0[j] - tg_xxxxz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xzz_1[j];

                    tg_xxxxxxz_xyyy_0[j] = pb_x * tg_xxxxxz_xyyy_0[j] + fr * tg_xxxxxz_xyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyy_0[j] - tg_xxxxz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyy_1[j];

                    tg_xxxxxxz_xyyz_0[j] = pb_x * tg_xxxxxz_xyyz_0[j] + fr * tg_xxxxxz_xyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyyz_0[j] - tg_xxxxz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yyz_1[j];

                    tg_xxxxxxz_xyzz_0[j] = pb_x * tg_xxxxxz_xyzz_0[j] + fr * tg_xxxxxz_xyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyzz_0[j] - tg_xxxxz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yzz_1[j];

                    tg_xxxxxxz_xzzz_0[j] = pb_x * tg_xxxxxz_xzzz_0[j] + fr * tg_xxxxxz_xzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xzzz_0[j] - tg_xxxxz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_zzz_1[j];

                    tg_xxxxxxz_yyyy_0[j] = pb_x * tg_xxxxxz_yyyy_0[j] + fr * tg_xxxxxz_yyyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyy_0[j] - tg_xxxxz_yyyy_1[j] * fl1_fza);

                    tg_xxxxxxz_yyyz_0[j] = pb_x * tg_xxxxxz_yyyz_0[j] + fr * tg_xxxxxz_yyyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyyz_0[j] - tg_xxxxz_yyyz_1[j] * fl1_fza);

                    tg_xxxxxxz_yyzz_0[j] = pb_x * tg_xxxxxz_yyzz_0[j] + fr * tg_xxxxxz_yyzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyzz_0[j] - tg_xxxxz_yyzz_1[j] * fl1_fza);

                    tg_xxxxxxz_yzzz_0[j] = pb_x * tg_xxxxxz_yzzz_0[j] + fr * tg_xxxxxz_yzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yzzz_0[j] - tg_xxxxz_yzzz_1[j] * fl1_fza);

                    tg_xxxxxxz_zzzz_0[j] = pb_x * tg_xxxxxz_zzzz_0[j] + fr * tg_xxxxxz_zzzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_zzzz_0[j] - tg_xxxxz_zzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_xxxx_0[j] = pb_x * tg_xxxxyy_xxxx_0[j] + fr * tg_xxxxyy_xxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxx_0[j] - tg_xxxyy_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyy_xxx_1[j];

                    tg_xxxxxyy_xxxy_0[j] = pb_x * tg_xxxxyy_xxxy_0[j] + fr * tg_xxxxyy_xxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxy_0[j] - tg_xxxyy_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxy_1[j];

                    tg_xxxxxyy_xxxz_0[j] = pb_x * tg_xxxxyy_xxxz_0[j] + fr * tg_xxxxyy_xxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxxz_0[j] - tg_xxxyy_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xxz_1[j];

                    tg_xxxxxyy_xxyy_0[j] = pb_x * tg_xxxxyy_xxyy_0[j] + fr * tg_xxxxyy_xxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyy_0[j] - tg_xxxyy_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyy_1[j];

                    tg_xxxxxyy_xxyz_0[j] = pb_x * tg_xxxxyy_xxyz_0[j] + fr * tg_xxxxyy_xxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxyz_0[j] - tg_xxxyy_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xyz_1[j];

                    tg_xxxxxyy_xxzz_0[j] = pb_x * tg_xxxxyy_xxzz_0[j] + fr * tg_xxxxyy_xxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxzz_0[j] - tg_xxxyy_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xzz_1[j];

                    tg_xxxxxyy_xyyy_0[j] = pb_x * tg_xxxxyy_xyyy_0[j] + fr * tg_xxxxyy_xyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyy_0[j] - tg_xxxyy_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyy_1[j];

                    tg_xxxxxyy_xyyz_0[j] = pb_x * tg_xxxxyy_xyyz_0[j] + fr * tg_xxxxyy_xyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyyz_0[j] - tg_xxxyy_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yyz_1[j];

                    tg_xxxxxyy_xyzz_0[j] = pb_x * tg_xxxxyy_xyzz_0[j] + fr * tg_xxxxyy_xyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyzz_0[j] - tg_xxxyy_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yzz_1[j];

                    tg_xxxxxyy_xzzz_0[j] = pb_x * tg_xxxxyy_xzzz_0[j] + fr * tg_xxxxyy_xzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xzzz_0[j] - tg_xxxyy_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_zzz_1[j];

                    tg_xxxxxyy_yyyy_0[j] = pb_x * tg_xxxxyy_yyyy_0[j] + fr * tg_xxxxyy_yyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyy_0[j] - tg_xxxyy_yyyy_1[j] * fl1_fza);

                    tg_xxxxxyy_yyyz_0[j] = pb_x * tg_xxxxyy_yyyz_0[j] + fr * tg_xxxxyy_yyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyyz_0[j] - tg_xxxyy_yyyz_1[j] * fl1_fza);

                    tg_xxxxxyy_yyzz_0[j] = pb_x * tg_xxxxyy_yyzz_0[j] + fr * tg_xxxxyy_yyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyzz_0[j] - tg_xxxyy_yyzz_1[j] * fl1_fza);

                    tg_xxxxxyy_yzzz_0[j] = pb_x * tg_xxxxyy_yzzz_0[j] + fr * tg_xxxxyy_yzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yzzz_0[j] - tg_xxxyy_yzzz_1[j] * fl1_fza);

                    tg_xxxxxyy_zzzz_0[j] = pb_x * tg_xxxxyy_zzzz_0[j] + fr * tg_xxxxyy_zzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_zzzz_0[j] - tg_xxxyy_zzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_xxxx_0[j] = pb_x * tg_xxxxyz_xxxx_0[j] + fr * tg_xxxxyz_xxxx_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxx_0[j] - tg_xxxyz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxyz_xxx_1[j];

                    tg_xxxxxyz_xxxy_0[j] = pb_x * tg_xxxxyz_xxxy_0[j] + fr * tg_xxxxyz_xxxy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxy_0[j] - tg_xxxyz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxy_1[j];

                    tg_xxxxxyz_xxxz_0[j] = pb_x * tg_xxxxyz_xxxz_0[j] + fr * tg_xxxxyz_xxxz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxxz_0[j] - tg_xxxyz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xxz_1[j];

                    tg_xxxxxyz_xxyy_0[j] = pb_x * tg_xxxxyz_xxyy_0[j] + fr * tg_xxxxyz_xxyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyy_0[j] - tg_xxxyz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyy_1[j];

                    tg_xxxxxyz_xxyz_0[j] = pb_x * tg_xxxxyz_xxyz_0[j] + fr * tg_xxxxyz_xxyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxyz_0[j] - tg_xxxyz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xyz_1[j];

                    tg_xxxxxyz_xxzz_0[j] = pb_x * tg_xxxxyz_xxzz_0[j] + fr * tg_xxxxyz_xxzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxzz_0[j] - tg_xxxyz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xzz_1[j];

                    tg_xxxxxyz_xyyy_0[j] = pb_x * tg_xxxxyz_xyyy_0[j] + fr * tg_xxxxyz_xyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyy_0[j] - tg_xxxyz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyy_1[j];

                    tg_xxxxxyz_xyyz_0[j] = pb_x * tg_xxxxyz_xyyz_0[j] + fr * tg_xxxxyz_xyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyyz_0[j] - tg_xxxyz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yyz_1[j];

                    tg_xxxxxyz_xyzz_0[j] = pb_x * tg_xxxxyz_xyzz_0[j] + fr * tg_xxxxyz_xyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyzz_0[j] - tg_xxxyz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yzz_1[j];

                    tg_xxxxxyz_xzzz_0[j] = pb_x * tg_xxxxyz_xzzz_0[j] + fr * tg_xxxxyz_xzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xzzz_0[j] - tg_xxxyz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_zzz_1[j];

                    tg_xxxxxyz_yyyy_0[j] = pb_x * tg_xxxxyz_yyyy_0[j] + fr * tg_xxxxyz_yyyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyy_0[j] - tg_xxxyz_yyyy_1[j] * fl1_fza);

                    tg_xxxxxyz_yyyz_0[j] = pb_x * tg_xxxxyz_yyyz_0[j] + fr * tg_xxxxyz_yyyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyyz_0[j] - tg_xxxyz_yyyz_1[j] * fl1_fza);

                    tg_xxxxxyz_yyzz_0[j] = pb_x * tg_xxxxyz_yyzz_0[j] + fr * tg_xxxxyz_yyzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyzz_0[j] - tg_xxxyz_yyzz_1[j] * fl1_fza);

                    tg_xxxxxyz_yzzz_0[j] = pb_x * tg_xxxxyz_yzzz_0[j] + fr * tg_xxxxyz_yzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yzzz_0[j] - tg_xxxyz_yzzz_1[j] * fl1_fza);

                    tg_xxxxxyz_zzzz_0[j] = pb_x * tg_xxxxyz_zzzz_0[j] + fr * tg_xxxxyz_zzzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_zzzz_0[j] - tg_xxxyz_zzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_xxxx_0[j] = pb_x * tg_xxxxzz_xxxx_0[j] + fr * tg_xxxxzz_xxxx_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxx_0[j] - tg_xxxzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxxzz_xxx_1[j];

                    tg_xxxxxzz_xxxy_0[j] = pb_x * tg_xxxxzz_xxxy_0[j] + fr * tg_xxxxzz_xxxy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxy_0[j] - tg_xxxzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxy_1[j];

                    tg_xxxxxzz_xxxz_0[j] = pb_x * tg_xxxxzz_xxxz_0[j] + fr * tg_xxxxzz_xxxz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxxz_0[j] - tg_xxxzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xxz_1[j];

                    tg_xxxxxzz_xxyy_0[j] = pb_x * tg_xxxxzz_xxyy_0[j] + fr * tg_xxxxzz_xxyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyy_0[j] - tg_xxxzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyy_1[j];

                    tg_xxxxxzz_xxyz_0[j] = pb_x * tg_xxxxzz_xxyz_0[j] + fr * tg_xxxxzz_xxyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxyz_0[j] - tg_xxxzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xyz_1[j];

                    tg_xxxxxzz_xxzz_0[j] = pb_x * tg_xxxxzz_xxzz_0[j] + fr * tg_xxxxzz_xxzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxzz_0[j] - tg_xxxzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xzz_1[j];

                    tg_xxxxxzz_xyyy_0[j] = pb_x * tg_xxxxzz_xyyy_0[j] + fr * tg_xxxxzz_xyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyy_0[j] - tg_xxxzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyy_1[j];

                    tg_xxxxxzz_xyyz_0[j] = pb_x * tg_xxxxzz_xyyz_0[j] + fr * tg_xxxxzz_xyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyyz_0[j] - tg_xxxzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yyz_1[j];

                    tg_xxxxxzz_xyzz_0[j] = pb_x * tg_xxxxzz_xyzz_0[j] + fr * tg_xxxxzz_xyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyzz_0[j] - tg_xxxzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yzz_1[j];

                    tg_xxxxxzz_xzzz_0[j] = pb_x * tg_xxxxzz_xzzz_0[j] + fr * tg_xxxxzz_xzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xzzz_0[j] - tg_xxxzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_zzz_1[j];

                    tg_xxxxxzz_yyyy_0[j] = pb_x * tg_xxxxzz_yyyy_0[j] + fr * tg_xxxxzz_yyyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyy_0[j] - tg_xxxzz_yyyy_1[j] * fl1_fza);

                    tg_xxxxxzz_yyyz_0[j] = pb_x * tg_xxxxzz_yyyz_0[j] + fr * tg_xxxxzz_yyyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyyz_0[j] - tg_xxxzz_yyyz_1[j] * fl1_fza);

                    tg_xxxxxzz_yyzz_0[j] = pb_x * tg_xxxxzz_yyzz_0[j] + fr * tg_xxxxzz_yyzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyzz_0[j] - tg_xxxzz_yyzz_1[j] * fl1_fza);

                    tg_xxxxxzz_yzzz_0[j] = pb_x * tg_xxxxzz_yzzz_0[j] + fr * tg_xxxxzz_yzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yzzz_0[j] - tg_xxxzz_yzzz_1[j] * fl1_fza);

                    tg_xxxxxzz_zzzz_0[j] = pb_x * tg_xxxxzz_zzzz_0[j] + fr * tg_xxxxzz_zzzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_zzzz_0[j] - tg_xxxzz_zzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSG_90_180(      CMemBlock2D<double>* primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (90,180)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxyyy_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 90); 

                auto tg_xxxyyy_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 91); 

                auto tg_xxxyyy_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 92); 

                auto tg_xxxyyy_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 93); 

                auto tg_xxxyyy_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 94); 

                auto tg_xxxyyy_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 95); 

                auto tg_xxxyyy_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 96); 

                auto tg_xxxyyy_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 97); 

                auto tg_xxxyyy_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 98); 

                auto tg_xxxyyy_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 99); 

                auto tg_xxxyyy_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 100); 

                auto tg_xxxyyy_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 101); 

                auto tg_xxxyyy_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 102); 

                auto tg_xxxyyy_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 103); 

                auto tg_xxxyyy_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 104); 

                auto tg_xxxyyz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 105); 

                auto tg_xxxyyz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 106); 

                auto tg_xxxyyz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 107); 

                auto tg_xxxyyz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 108); 

                auto tg_xxxyyz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 109); 

                auto tg_xxxyyz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 110); 

                auto tg_xxxyyz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 111); 

                auto tg_xxxyyz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 112); 

                auto tg_xxxyyz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 113); 

                auto tg_xxxyyz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 114); 

                auto tg_xxxyyz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 115); 

                auto tg_xxxyyz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 116); 

                auto tg_xxxyyz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 117); 

                auto tg_xxxyyz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 118); 

                auto tg_xxxyyz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 119); 

                auto tg_xxxyzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 120); 

                auto tg_xxxyzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 121); 

                auto tg_xxxyzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 122); 

                auto tg_xxxyzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 123); 

                auto tg_xxxyzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 124); 

                auto tg_xxxyzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 125); 

                auto tg_xxxyzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 126); 

                auto tg_xxxyzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 127); 

                auto tg_xxxyzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 128); 

                auto tg_xxxyzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 129); 

                auto tg_xxxyzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 130); 

                auto tg_xxxyzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 131); 

                auto tg_xxxyzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 132); 

                auto tg_xxxyzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 133); 

                auto tg_xxxyzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 134); 

                auto tg_xxxzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 135); 

                auto tg_xxxzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 136); 

                auto tg_xxxzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 137); 

                auto tg_xxxzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 138); 

                auto tg_xxxzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 139); 

                auto tg_xxxzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 140); 

                auto tg_xxxzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 141); 

                auto tg_xxxzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 142); 

                auto tg_xxxzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 143); 

                auto tg_xxxzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 144); 

                auto tg_xxxzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 145); 

                auto tg_xxxzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 146); 

                auto tg_xxxzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 147); 

                auto tg_xxxzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 148); 

                auto tg_xxxzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 149); 

                auto tg_xxyyyy_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 150); 

                auto tg_xxyyyy_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 151); 

                auto tg_xxyyyy_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 152); 

                auto tg_xxyyyy_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 153); 

                auto tg_xxyyyy_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 154); 

                auto tg_xxyyyy_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 155); 

                auto tg_xxyyyy_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 156); 

                auto tg_xxyyyy_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 157); 

                auto tg_xxyyyy_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 158); 

                auto tg_xxyyyy_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 159); 

                auto tg_xxyyyy_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 160); 

                auto tg_xxyyyy_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 161); 

                auto tg_xxyyyy_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 162); 

                auto tg_xxyyyy_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 163); 

                auto tg_xxyyyy_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 164); 

                auto tg_xxyyyz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 165); 

                auto tg_xxyyyz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 166); 

                auto tg_xxyyyz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 167); 

                auto tg_xxyyyz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 168); 

                auto tg_xxyyyz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 169); 

                auto tg_xxyyyz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 170); 

                auto tg_xxyyyz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 171); 

                auto tg_xxyyyz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 172); 

                auto tg_xxyyyz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 173); 

                auto tg_xxyyyz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 174); 

                auto tg_xxyyyz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 175); 

                auto tg_xxyyyz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 176); 

                auto tg_xxyyyz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 177); 

                auto tg_xxyyyz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 178); 

                auto tg_xxyyyz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 179); 

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

                auto tg_xxyyy_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 90); 

                auto tg_xxyyy_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 91); 

                auto tg_xxyyy_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 92); 

                auto tg_xxyyy_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 93); 

                auto tg_xxyyy_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 94); 

                auto tg_xxyyy_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 95); 

                auto tg_xxyyy_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 96); 

                auto tg_xxyyy_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 97); 

                auto tg_xxyyy_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 98); 

                auto tg_xxyyy_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 99); 

                auto tg_xxyyy_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 100); 

                auto tg_xxyyy_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 101); 

                auto tg_xxyyy_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 102); 

                auto tg_xxyyy_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 103); 

                auto tg_xxyyy_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 104); 

                auto tg_xxyyz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 105); 

                auto tg_xxyyz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 106); 

                auto tg_xxyyz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 107); 

                auto tg_xxyyz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 108); 

                auto tg_xxyyz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 109); 

                auto tg_xxyyz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 110); 

                auto tg_xxyyz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 111); 

                auto tg_xxyyz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 112); 

                auto tg_xxyyz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 113); 

                auto tg_xxyyz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 114); 

                auto tg_xxyyz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 115); 

                auto tg_xxyyz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 116); 

                auto tg_xxyyz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 117); 

                auto tg_xxyyz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 118); 

                auto tg_xxyyz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 119); 

                auto tg_xxyzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 120); 

                auto tg_xxyzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 121); 

                auto tg_xxyzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 122); 

                auto tg_xxyzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 123); 

                auto tg_xxyzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 124); 

                auto tg_xxyzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 125); 

                auto tg_xxyzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 126); 

                auto tg_xxyzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 127); 

                auto tg_xxyzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 128); 

                auto tg_xxyzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 129); 

                auto tg_xxyzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 130); 

                auto tg_xxyzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 131); 

                auto tg_xxyzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 132); 

                auto tg_xxyzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 133); 

                auto tg_xxyzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 134); 

                auto tg_xxzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 135); 

                auto tg_xxzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 136); 

                auto tg_xxzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 137); 

                auto tg_xxzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 138); 

                auto tg_xxzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 139); 

                auto tg_xxzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 140); 

                auto tg_xxzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 141); 

                auto tg_xxzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 142); 

                auto tg_xxzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 143); 

                auto tg_xxzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 144); 

                auto tg_xxzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 145); 

                auto tg_xxzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 146); 

                auto tg_xxzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 147); 

                auto tg_xxzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 148); 

                auto tg_xxzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 149); 

                auto tg_xyyyy_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 150); 

                auto tg_xyyyy_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 151); 

                auto tg_xyyyy_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 152); 

                auto tg_xyyyy_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 153); 

                auto tg_xyyyy_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 154); 

                auto tg_xyyyy_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 155); 

                auto tg_xyyyy_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 156); 

                auto tg_xyyyy_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 157); 

                auto tg_xyyyy_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 158); 

                auto tg_xyyyy_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 159); 

                auto tg_xyyyy_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 160); 

                auto tg_xyyyy_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 161); 

                auto tg_xyyyy_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 162); 

                auto tg_xyyyy_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 163); 

                auto tg_xyyyy_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 164); 

                auto tg_xyyyz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 165); 

                auto tg_xyyyz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 166); 

                auto tg_xyyyz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 167); 

                auto tg_xyyyz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 168); 

                auto tg_xyyyz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 169); 

                auto tg_xyyyz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 170); 

                auto tg_xyyyz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 171); 

                auto tg_xyyyz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 172); 

                auto tg_xyyyz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 173); 

                auto tg_xyyyz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 174); 

                auto tg_xyyyz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 175); 

                auto tg_xyyyz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 176); 

                auto tg_xyyyz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 177); 

                auto tg_xyyyz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 178); 

                auto tg_xyyyz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 179); 

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

                auto tg_xxxyyy_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 60); 

                auto tg_xxxyyy_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 61); 

                auto tg_xxxyyy_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 62); 

                auto tg_xxxyyy_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 63); 

                auto tg_xxxyyy_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 64); 

                auto tg_xxxyyy_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 65); 

                auto tg_xxxyyy_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 66); 

                auto tg_xxxyyy_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 67); 

                auto tg_xxxyyy_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 68); 

                auto tg_xxxyyy_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 69); 

                auto tg_xxxyyz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 70); 

                auto tg_xxxyyz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 71); 

                auto tg_xxxyyz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 72); 

                auto tg_xxxyyz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 73); 

                auto tg_xxxyyz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 74); 

                auto tg_xxxyyz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 75); 

                auto tg_xxxyyz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 76); 

                auto tg_xxxyyz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 77); 

                auto tg_xxxyyz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 78); 

                auto tg_xxxyyz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 79); 

                auto tg_xxxyzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 80); 

                auto tg_xxxyzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 81); 

                auto tg_xxxyzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 82); 

                auto tg_xxxyzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 83); 

                auto tg_xxxyzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 84); 

                auto tg_xxxyzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 85); 

                auto tg_xxxyzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 86); 

                auto tg_xxxyzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 87); 

                auto tg_xxxyzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 88); 

                auto tg_xxxyzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 89); 

                auto tg_xxxzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 90); 

                auto tg_xxxzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 91); 

                auto tg_xxxzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 92); 

                auto tg_xxxzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 93); 

                auto tg_xxxzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 94); 

                auto tg_xxxzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 95); 

                auto tg_xxxzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 96); 

                auto tg_xxxzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 97); 

                auto tg_xxxzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 98); 

                auto tg_xxxzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 99); 

                auto tg_xxyyyy_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 100); 

                auto tg_xxyyyy_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 101); 

                auto tg_xxyyyy_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 102); 

                auto tg_xxyyyy_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 103); 

                auto tg_xxyyyy_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 104); 

                auto tg_xxyyyy_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 105); 

                auto tg_xxyyyy_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 106); 

                auto tg_xxyyyy_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 107); 

                auto tg_xxyyyy_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 108); 

                auto tg_xxyyyy_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 109); 

                auto tg_xxyyyz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 110); 

                auto tg_xxyyyz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 111); 

                auto tg_xxyyyz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 112); 

                auto tg_xxyyyz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 113); 

                auto tg_xxyyyz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 114); 

                auto tg_xxyyyz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 115); 

                auto tg_xxyyyz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 116); 

                auto tg_xxyyyz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 117); 

                auto tg_xxyyyz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 118); 

                auto tg_xxyyyz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 119); 

                // set up pointers to integrals

                auto tg_xxxxyyy_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 90); 

                auto tg_xxxxyyy_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 91); 

                auto tg_xxxxyyy_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 92); 

                auto tg_xxxxyyy_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 93); 

                auto tg_xxxxyyy_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 94); 

                auto tg_xxxxyyy_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 95); 

                auto tg_xxxxyyy_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 96); 

                auto tg_xxxxyyy_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 97); 

                auto tg_xxxxyyy_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 98); 

                auto tg_xxxxyyy_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 99); 

                auto tg_xxxxyyy_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 100); 

                auto tg_xxxxyyy_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 101); 

                auto tg_xxxxyyy_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 102); 

                auto tg_xxxxyyy_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 103); 

                auto tg_xxxxyyy_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 104); 

                auto tg_xxxxyyz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 105); 

                auto tg_xxxxyyz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 106); 

                auto tg_xxxxyyz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 107); 

                auto tg_xxxxyyz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 108); 

                auto tg_xxxxyyz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 109); 

                auto tg_xxxxyyz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 110); 

                auto tg_xxxxyyz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 111); 

                auto tg_xxxxyyz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 112); 

                auto tg_xxxxyyz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 113); 

                auto tg_xxxxyyz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 114); 

                auto tg_xxxxyyz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 115); 

                auto tg_xxxxyyz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 116); 

                auto tg_xxxxyyz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 117); 

                auto tg_xxxxyyz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 118); 

                auto tg_xxxxyyz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 119); 

                auto tg_xxxxyzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 120); 

                auto tg_xxxxyzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 121); 

                auto tg_xxxxyzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 122); 

                auto tg_xxxxyzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 123); 

                auto tg_xxxxyzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 124); 

                auto tg_xxxxyzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 125); 

                auto tg_xxxxyzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 126); 

                auto tg_xxxxyzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 127); 

                auto tg_xxxxyzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 128); 

                auto tg_xxxxyzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 129); 

                auto tg_xxxxyzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 130); 

                auto tg_xxxxyzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 131); 

                auto tg_xxxxyzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 132); 

                auto tg_xxxxyzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 133); 

                auto tg_xxxxyzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 134); 

                auto tg_xxxxzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 135); 

                auto tg_xxxxzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 136); 

                auto tg_xxxxzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 137); 

                auto tg_xxxxzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 138); 

                auto tg_xxxxzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 139); 

                auto tg_xxxxzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 140); 

                auto tg_xxxxzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 141); 

                auto tg_xxxxzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 142); 

                auto tg_xxxxzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 143); 

                auto tg_xxxxzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 144); 

                auto tg_xxxxzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 145); 

                auto tg_xxxxzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 146); 

                auto tg_xxxxzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 147); 

                auto tg_xxxxzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 148); 

                auto tg_xxxxzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 149); 

                auto tg_xxxyyyy_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 150); 

                auto tg_xxxyyyy_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 151); 

                auto tg_xxxyyyy_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 152); 

                auto tg_xxxyyyy_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 153); 

                auto tg_xxxyyyy_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 154); 

                auto tg_xxxyyyy_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 155); 

                auto tg_xxxyyyy_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 156); 

                auto tg_xxxyyyy_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 157); 

                auto tg_xxxyyyy_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 158); 

                auto tg_xxxyyyy_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 159); 

                auto tg_xxxyyyy_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 160); 

                auto tg_xxxyyyy_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 161); 

                auto tg_xxxyyyy_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 162); 

                auto tg_xxxyyyy_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 163); 

                auto tg_xxxyyyy_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 164); 

                auto tg_xxxyyyz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 165); 

                auto tg_xxxyyyz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 166); 

                auto tg_xxxyyyz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 167); 

                auto tg_xxxyyyz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 168); 

                auto tg_xxxyyyz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 169); 

                auto tg_xxxyyyz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 170); 

                auto tg_xxxyyyz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 171); 

                auto tg_xxxyyyz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 172); 

                auto tg_xxxyyyz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 173); 

                auto tg_xxxyyyz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 174); 

                auto tg_xxxyyyz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 175); 

                auto tg_xxxyyyz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 176); 

                auto tg_xxxyyyz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 177); 

                auto tg_xxxyyyz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 178); 

                auto tg_xxxyyyz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 179); 

                // Batch of Integrals (90,180)

                #pragma omp simd aligned(fxn, fza, tg_xxxxyyy_xxxx_0, tg_xxxxyyy_xxxy_0, tg_xxxxyyy_xxxz_0, \
                                         tg_xxxxyyy_xxyy_0, tg_xxxxyyy_xxyz_0, tg_xxxxyyy_xxzz_0, tg_xxxxyyy_xyyy_0, \
                                         tg_xxxxyyy_xyyz_0, tg_xxxxyyy_xyzz_0, tg_xxxxyyy_xzzz_0, tg_xxxxyyy_yyyy_0, \
                                         tg_xxxxyyy_yyyz_0, tg_xxxxyyy_yyzz_0, tg_xxxxyyy_yzzz_0, tg_xxxxyyy_zzzz_0, \
                                         tg_xxxxyyz_xxxx_0, tg_xxxxyyz_xxxy_0, tg_xxxxyyz_xxxz_0, tg_xxxxyyz_xxyy_0, \
                                         tg_xxxxyyz_xxyz_0, tg_xxxxyyz_xxzz_0, tg_xxxxyyz_xyyy_0, tg_xxxxyyz_xyyz_0, \
                                         tg_xxxxyyz_xyzz_0, tg_xxxxyyz_xzzz_0, tg_xxxxyyz_yyyy_0, tg_xxxxyyz_yyyz_0, \
                                         tg_xxxxyyz_yyzz_0, tg_xxxxyyz_yzzz_0, tg_xxxxyyz_zzzz_0, tg_xxxxyzz_xxxx_0, \
                                         tg_xxxxyzz_xxxy_0, tg_xxxxyzz_xxxz_0, tg_xxxxyzz_xxyy_0, tg_xxxxyzz_xxyz_0, \
                                         tg_xxxxyzz_xxzz_0, tg_xxxxyzz_xyyy_0, tg_xxxxyzz_xyyz_0, tg_xxxxyzz_xyzz_0, \
                                         tg_xxxxyzz_xzzz_0, tg_xxxxyzz_yyyy_0, tg_xxxxyzz_yyyz_0, tg_xxxxyzz_yyzz_0, \
                                         tg_xxxxyzz_yzzz_0, tg_xxxxyzz_zzzz_0, tg_xxxxzzz_xxxx_0, tg_xxxxzzz_xxxy_0, \
                                         tg_xxxxzzz_xxxz_0, tg_xxxxzzz_xxyy_0, tg_xxxxzzz_xxyz_0, tg_xxxxzzz_xxzz_0, \
                                         tg_xxxxzzz_xyyy_0, tg_xxxxzzz_xyyz_0, tg_xxxxzzz_xyzz_0, tg_xxxxzzz_xzzz_0, \
                                         tg_xxxxzzz_yyyy_0, tg_xxxxzzz_yyyz_0, tg_xxxxzzz_yyzz_0, tg_xxxxzzz_yzzz_0, \
                                         tg_xxxxzzz_zzzz_0, tg_xxxyyy_xxx_1, tg_xxxyyy_xxxx_0, tg_xxxyyy_xxxx_1, \
                                         tg_xxxyyy_xxxy_0, tg_xxxyyy_xxxy_1, tg_xxxyyy_xxxz_0, tg_xxxyyy_xxxz_1, \
                                         tg_xxxyyy_xxy_1, tg_xxxyyy_xxyy_0, tg_xxxyyy_xxyy_1, tg_xxxyyy_xxyz_0, \
                                         tg_xxxyyy_xxyz_1, tg_xxxyyy_xxz_1, tg_xxxyyy_xxzz_0, tg_xxxyyy_xxzz_1, \
                                         tg_xxxyyy_xyy_1, tg_xxxyyy_xyyy_0, tg_xxxyyy_xyyy_1, tg_xxxyyy_xyyz_0, \
                                         tg_xxxyyy_xyyz_1, tg_xxxyyy_xyz_1, tg_xxxyyy_xyzz_0, tg_xxxyyy_xyzz_1, \
                                         tg_xxxyyy_xzz_1, tg_xxxyyy_xzzz_0, tg_xxxyyy_xzzz_1, tg_xxxyyy_yyy_1, \
                                         tg_xxxyyy_yyyy_0, tg_xxxyyy_yyyy_1, tg_xxxyyy_yyyz_0, tg_xxxyyy_yyyz_1, \
                                         tg_xxxyyy_yyz_1, tg_xxxyyy_yyzz_0, tg_xxxyyy_yyzz_1, tg_xxxyyy_yzz_1, \
                                         tg_xxxyyy_yzzz_0, tg_xxxyyy_yzzz_1, tg_xxxyyy_zzz_1, tg_xxxyyy_zzzz_0, \
                                         tg_xxxyyy_zzzz_1, tg_xxxyyyy_xxxx_0, tg_xxxyyyy_xxxy_0, tg_xxxyyyy_xxxz_0, \
                                         tg_xxxyyyy_xxyy_0, tg_xxxyyyy_xxyz_0, tg_xxxyyyy_xxzz_0, tg_xxxyyyy_xyyy_0, \
                                         tg_xxxyyyy_xyyz_0, tg_xxxyyyy_xyzz_0, tg_xxxyyyy_xzzz_0, tg_xxxyyyy_yyyy_0, \
                                         tg_xxxyyyy_yyyz_0, tg_xxxyyyy_yyzz_0, tg_xxxyyyy_yzzz_0, tg_xxxyyyy_zzzz_0, \
                                         tg_xxxyyyz_xxxx_0, tg_xxxyyyz_xxxy_0, tg_xxxyyyz_xxxz_0, tg_xxxyyyz_xxyy_0, \
                                         tg_xxxyyyz_xxyz_0, tg_xxxyyyz_xxzz_0, tg_xxxyyyz_xyyy_0, tg_xxxyyyz_xyyz_0, \
                                         tg_xxxyyyz_xyzz_0, tg_xxxyyyz_xzzz_0, tg_xxxyyyz_yyyy_0, tg_xxxyyyz_yyyz_0, \
                                         tg_xxxyyyz_yyzz_0, tg_xxxyyyz_yzzz_0, tg_xxxyyyz_zzzz_0, tg_xxxyyz_xxx_1, \
                                         tg_xxxyyz_xxxx_0, tg_xxxyyz_xxxx_1, tg_xxxyyz_xxxy_0, tg_xxxyyz_xxxy_1, \
                                         tg_xxxyyz_xxxz_0, tg_xxxyyz_xxxz_1, tg_xxxyyz_xxy_1, tg_xxxyyz_xxyy_0, \
                                         tg_xxxyyz_xxyy_1, tg_xxxyyz_xxyz_0, tg_xxxyyz_xxyz_1, tg_xxxyyz_xxz_1, \
                                         tg_xxxyyz_xxzz_0, tg_xxxyyz_xxzz_1, tg_xxxyyz_xyy_1, tg_xxxyyz_xyyy_0, \
                                         tg_xxxyyz_xyyy_1, tg_xxxyyz_xyyz_0, tg_xxxyyz_xyyz_1, tg_xxxyyz_xyz_1, \
                                         tg_xxxyyz_xyzz_0, tg_xxxyyz_xyzz_1, tg_xxxyyz_xzz_1, tg_xxxyyz_xzzz_0, \
                                         tg_xxxyyz_xzzz_1, tg_xxxyyz_yyy_1, tg_xxxyyz_yyyy_0, tg_xxxyyz_yyyy_1, \
                                         tg_xxxyyz_yyyz_0, tg_xxxyyz_yyyz_1, tg_xxxyyz_yyz_1, tg_xxxyyz_yyzz_0, \
                                         tg_xxxyyz_yyzz_1, tg_xxxyyz_yzz_1, tg_xxxyyz_yzzz_0, tg_xxxyyz_yzzz_1, \
                                         tg_xxxyyz_zzz_1, tg_xxxyyz_zzzz_0, tg_xxxyyz_zzzz_1, tg_xxxyzz_xxx_1, \
                                         tg_xxxyzz_xxxx_0, tg_xxxyzz_xxxx_1, tg_xxxyzz_xxxy_0, tg_xxxyzz_xxxy_1, \
                                         tg_xxxyzz_xxxz_0, tg_xxxyzz_xxxz_1, tg_xxxyzz_xxy_1, tg_xxxyzz_xxyy_0, \
                                         tg_xxxyzz_xxyy_1, tg_xxxyzz_xxyz_0, tg_xxxyzz_xxyz_1, tg_xxxyzz_xxz_1, \
                                         tg_xxxyzz_xxzz_0, tg_xxxyzz_xxzz_1, tg_xxxyzz_xyy_1, tg_xxxyzz_xyyy_0, \
                                         tg_xxxyzz_xyyy_1, tg_xxxyzz_xyyz_0, tg_xxxyzz_xyyz_1, tg_xxxyzz_xyz_1, \
                                         tg_xxxyzz_xyzz_0, tg_xxxyzz_xyzz_1, tg_xxxyzz_xzz_1, tg_xxxyzz_xzzz_0, \
                                         tg_xxxyzz_xzzz_1, tg_xxxyzz_yyy_1, tg_xxxyzz_yyyy_0, tg_xxxyzz_yyyy_1, \
                                         tg_xxxyzz_yyyz_0, tg_xxxyzz_yyyz_1, tg_xxxyzz_yyz_1, tg_xxxyzz_yyzz_0, \
                                         tg_xxxyzz_yyzz_1, tg_xxxyzz_yzz_1, tg_xxxyzz_yzzz_0, tg_xxxyzz_yzzz_1, \
                                         tg_xxxyzz_zzz_1, tg_xxxyzz_zzzz_0, tg_xxxyzz_zzzz_1, tg_xxxzzz_xxx_1, \
                                         tg_xxxzzz_xxxx_0, tg_xxxzzz_xxxx_1, tg_xxxzzz_xxxy_0, tg_xxxzzz_xxxy_1, \
                                         tg_xxxzzz_xxxz_0, tg_xxxzzz_xxxz_1, tg_xxxzzz_xxy_1, tg_xxxzzz_xxyy_0, \
                                         tg_xxxzzz_xxyy_1, tg_xxxzzz_xxyz_0, tg_xxxzzz_xxyz_1, tg_xxxzzz_xxz_1, \
                                         tg_xxxzzz_xxzz_0, tg_xxxzzz_xxzz_1, tg_xxxzzz_xyy_1, tg_xxxzzz_xyyy_0, \
                                         tg_xxxzzz_xyyy_1, tg_xxxzzz_xyyz_0, tg_xxxzzz_xyyz_1, tg_xxxzzz_xyz_1, \
                                         tg_xxxzzz_xyzz_0, tg_xxxzzz_xyzz_1, tg_xxxzzz_xzz_1, tg_xxxzzz_xzzz_0, \
                                         tg_xxxzzz_xzzz_1, tg_xxxzzz_yyy_1, tg_xxxzzz_yyyy_0, tg_xxxzzz_yyyy_1, \
                                         tg_xxxzzz_yyyz_0, tg_xxxzzz_yyyz_1, tg_xxxzzz_yyz_1, tg_xxxzzz_yyzz_0, \
                                         tg_xxxzzz_yyzz_1, tg_xxxzzz_yzz_1, tg_xxxzzz_yzzz_0, tg_xxxzzz_yzzz_1, \
                                         tg_xxxzzz_zzz_1, tg_xxxzzz_zzzz_0, tg_xxxzzz_zzzz_1, tg_xxyyy_xxxx_0, \
                                         tg_xxyyy_xxxx_1, tg_xxyyy_xxxy_0, tg_xxyyy_xxxy_1, tg_xxyyy_xxxz_0, tg_xxyyy_xxxz_1, \
                                         tg_xxyyy_xxyy_0, tg_xxyyy_xxyy_1, tg_xxyyy_xxyz_0, tg_xxyyy_xxyz_1, tg_xxyyy_xxzz_0, \
                                         tg_xxyyy_xxzz_1, tg_xxyyy_xyyy_0, tg_xxyyy_xyyy_1, tg_xxyyy_xyyz_0, tg_xxyyy_xyyz_1, \
                                         tg_xxyyy_xyzz_0, tg_xxyyy_xyzz_1, tg_xxyyy_xzzz_0, tg_xxyyy_xzzz_1, tg_xxyyy_yyyy_0, \
                                         tg_xxyyy_yyyy_1, tg_xxyyy_yyyz_0, tg_xxyyy_yyyz_1, tg_xxyyy_yyzz_0, tg_xxyyy_yyzz_1, \
                                         tg_xxyyy_yzzz_0, tg_xxyyy_yzzz_1, tg_xxyyy_zzzz_0, tg_xxyyy_zzzz_1, tg_xxyyyy_xxx_1, \
                                         tg_xxyyyy_xxxx_0, tg_xxyyyy_xxxx_1, tg_xxyyyy_xxxy_0, tg_xxyyyy_xxxy_1, \
                                         tg_xxyyyy_xxxz_0, tg_xxyyyy_xxxz_1, tg_xxyyyy_xxy_1, tg_xxyyyy_xxyy_0, \
                                         tg_xxyyyy_xxyy_1, tg_xxyyyy_xxyz_0, tg_xxyyyy_xxyz_1, tg_xxyyyy_xxz_1, \
                                         tg_xxyyyy_xxzz_0, tg_xxyyyy_xxzz_1, tg_xxyyyy_xyy_1, tg_xxyyyy_xyyy_0, \
                                         tg_xxyyyy_xyyy_1, tg_xxyyyy_xyyz_0, tg_xxyyyy_xyyz_1, tg_xxyyyy_xyz_1, \
                                         tg_xxyyyy_xyzz_0, tg_xxyyyy_xyzz_1, tg_xxyyyy_xzz_1, tg_xxyyyy_xzzz_0, \
                                         tg_xxyyyy_xzzz_1, tg_xxyyyy_yyy_1, tg_xxyyyy_yyyy_0, tg_xxyyyy_yyyy_1, \
                                         tg_xxyyyy_yyyz_0, tg_xxyyyy_yyyz_1, tg_xxyyyy_yyz_1, tg_xxyyyy_yyzz_0, \
                                         tg_xxyyyy_yyzz_1, tg_xxyyyy_yzz_1, tg_xxyyyy_yzzz_0, tg_xxyyyy_yzzz_1, \
                                         tg_xxyyyy_zzz_1, tg_xxyyyy_zzzz_0, tg_xxyyyy_zzzz_1, tg_xxyyyz_xxx_1, \
                                         tg_xxyyyz_xxxx_0, tg_xxyyyz_xxxx_1, tg_xxyyyz_xxxy_0, tg_xxyyyz_xxxy_1, \
                                         tg_xxyyyz_xxxz_0, tg_xxyyyz_xxxz_1, tg_xxyyyz_xxy_1, tg_xxyyyz_xxyy_0, \
                                         tg_xxyyyz_xxyy_1, tg_xxyyyz_xxyz_0, tg_xxyyyz_xxyz_1, tg_xxyyyz_xxz_1, \
                                         tg_xxyyyz_xxzz_0, tg_xxyyyz_xxzz_1, tg_xxyyyz_xyy_1, tg_xxyyyz_xyyy_0, \
                                         tg_xxyyyz_xyyy_1, tg_xxyyyz_xyyz_0, tg_xxyyyz_xyyz_1, tg_xxyyyz_xyz_1, \
                                         tg_xxyyyz_xyzz_0, tg_xxyyyz_xyzz_1, tg_xxyyyz_xzz_1, tg_xxyyyz_xzzz_0, \
                                         tg_xxyyyz_xzzz_1, tg_xxyyyz_yyy_1, tg_xxyyyz_yyyy_0, tg_xxyyyz_yyyy_1, \
                                         tg_xxyyyz_yyyz_0, tg_xxyyyz_yyyz_1, tg_xxyyyz_yyz_1, tg_xxyyyz_yyzz_0, \
                                         tg_xxyyyz_yyzz_1, tg_xxyyyz_yzz_1, tg_xxyyyz_yzzz_0, tg_xxyyyz_yzzz_1, \
                                         tg_xxyyyz_zzz_1, tg_xxyyyz_zzzz_0, tg_xxyyyz_zzzz_1, tg_xxyyz_xxxx_0, \
                                         tg_xxyyz_xxxx_1, tg_xxyyz_xxxy_0, tg_xxyyz_xxxy_1, tg_xxyyz_xxxz_0, tg_xxyyz_xxxz_1, \
                                         tg_xxyyz_xxyy_0, tg_xxyyz_xxyy_1, tg_xxyyz_xxyz_0, tg_xxyyz_xxyz_1, tg_xxyyz_xxzz_0, \
                                         tg_xxyyz_xxzz_1, tg_xxyyz_xyyy_0, tg_xxyyz_xyyy_1, tg_xxyyz_xyyz_0, tg_xxyyz_xyyz_1, \
                                         tg_xxyyz_xyzz_0, tg_xxyyz_xyzz_1, tg_xxyyz_xzzz_0, tg_xxyyz_xzzz_1, tg_xxyyz_yyyy_0, \
                                         tg_xxyyz_yyyy_1, tg_xxyyz_yyyz_0, tg_xxyyz_yyyz_1, tg_xxyyz_yyzz_0, tg_xxyyz_yyzz_1, \
                                         tg_xxyyz_yzzz_0, tg_xxyyz_yzzz_1, tg_xxyyz_zzzz_0, tg_xxyyz_zzzz_1, tg_xxyzz_xxxx_0, \
                                         tg_xxyzz_xxxx_1, tg_xxyzz_xxxy_0, tg_xxyzz_xxxy_1, tg_xxyzz_xxxz_0, tg_xxyzz_xxxz_1, \
                                         tg_xxyzz_xxyy_0, tg_xxyzz_xxyy_1, tg_xxyzz_xxyz_0, tg_xxyzz_xxyz_1, tg_xxyzz_xxzz_0, \
                                         tg_xxyzz_xxzz_1, tg_xxyzz_xyyy_0, tg_xxyzz_xyyy_1, tg_xxyzz_xyyz_0, tg_xxyzz_xyyz_1, \
                                         tg_xxyzz_xyzz_0, tg_xxyzz_xyzz_1, tg_xxyzz_xzzz_0, tg_xxyzz_xzzz_1, tg_xxyzz_yyyy_0, \
                                         tg_xxyzz_yyyy_1, tg_xxyzz_yyyz_0, tg_xxyzz_yyyz_1, tg_xxyzz_yyzz_0, tg_xxyzz_yyzz_1, \
                                         tg_xxyzz_yzzz_0, tg_xxyzz_yzzz_1, tg_xxyzz_zzzz_0, tg_xxyzz_zzzz_1, tg_xxzzz_xxxx_0, \
                                         tg_xxzzz_xxxx_1, tg_xxzzz_xxxy_0, tg_xxzzz_xxxy_1, tg_xxzzz_xxxz_0, tg_xxzzz_xxxz_1, \
                                         tg_xxzzz_xxyy_0, tg_xxzzz_xxyy_1, tg_xxzzz_xxyz_0, tg_xxzzz_xxyz_1, tg_xxzzz_xxzz_0, \
                                         tg_xxzzz_xxzz_1, tg_xxzzz_xyyy_0, tg_xxzzz_xyyy_1, tg_xxzzz_xyyz_0, tg_xxzzz_xyyz_1, \
                                         tg_xxzzz_xyzz_0, tg_xxzzz_xyzz_1, tg_xxzzz_xzzz_0, tg_xxzzz_xzzz_1, tg_xxzzz_yyyy_0, \
                                         tg_xxzzz_yyyy_1, tg_xxzzz_yyyz_0, tg_xxzzz_yyyz_1, tg_xxzzz_yyzz_0, tg_xxzzz_yyzz_1, \
                                         tg_xxzzz_yzzz_0, tg_xxzzz_yzzz_1, tg_xxzzz_zzzz_0, tg_xxzzz_zzzz_1, tg_xyyyy_xxxx_0, \
                                         tg_xyyyy_xxxx_1, tg_xyyyy_xxxy_0, tg_xyyyy_xxxy_1, tg_xyyyy_xxxz_0, tg_xyyyy_xxxz_1, \
                                         tg_xyyyy_xxyy_0, tg_xyyyy_xxyy_1, tg_xyyyy_xxyz_0, tg_xyyyy_xxyz_1, tg_xyyyy_xxzz_0, \
                                         tg_xyyyy_xxzz_1, tg_xyyyy_xyyy_0, tg_xyyyy_xyyy_1, tg_xyyyy_xyyz_0, tg_xyyyy_xyyz_1, \
                                         tg_xyyyy_xyzz_0, tg_xyyyy_xyzz_1, tg_xyyyy_xzzz_0, tg_xyyyy_xzzz_1, tg_xyyyy_yyyy_0, \
                                         tg_xyyyy_yyyy_1, tg_xyyyy_yyyz_0, tg_xyyyy_yyyz_1, tg_xyyyy_yyzz_0, tg_xyyyy_yyzz_1, \
                                         tg_xyyyy_yzzz_0, tg_xyyyy_yzzz_1, tg_xyyyy_zzzz_0, tg_xyyyy_zzzz_1, tg_xyyyz_xxxx_0, \
                                         tg_xyyyz_xxxx_1, tg_xyyyz_xxxy_0, tg_xyyyz_xxxy_1, tg_xyyyz_xxxz_0, tg_xyyyz_xxxz_1, \
                                         tg_xyyyz_xxyy_0, tg_xyyyz_xxyy_1, tg_xyyyz_xxyz_0, tg_xyyyz_xxyz_1, tg_xyyyz_xxzz_0, \
                                         tg_xyyyz_xxzz_1, tg_xyyyz_xyyy_0, tg_xyyyz_xyyy_1, tg_xyyyz_xyyz_0, tg_xyyyz_xyyz_1, \
                                         tg_xyyyz_xyzz_0, tg_xyyyz_xyzz_1, tg_xyyyz_xzzz_0, tg_xyyyz_xzzz_1, tg_xyyyz_yyyy_0, \
                                         tg_xyyyz_yyyy_1, tg_xyyyz_yyyz_0, tg_xyyyz_yyyz_1, tg_xyyyz_yyzz_0, tg_xyyyz_yyzz_1, \
                                         tg_xyyyz_yzzz_0, tg_xyyyz_yzzz_1, tg_xyyyz_zzzz_0, tg_xyyyz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxyyy_xxxx_0[j] = pb_x * tg_xxxyyy_xxxx_0[j] + fr * tg_xxxyyy_xxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxx_0[j] - tg_xxyyy_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyy_xxx_1[j];

                    tg_xxxxyyy_xxxy_0[j] = pb_x * tg_xxxyyy_xxxy_0[j] + fr * tg_xxxyyy_xxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxy_0[j] - tg_xxyyy_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxy_1[j];

                    tg_xxxxyyy_xxxz_0[j] = pb_x * tg_xxxyyy_xxxz_0[j] + fr * tg_xxxyyy_xxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxxz_0[j] - tg_xxyyy_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xxz_1[j];

                    tg_xxxxyyy_xxyy_0[j] = pb_x * tg_xxxyyy_xxyy_0[j] + fr * tg_xxxyyy_xxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyy_0[j] - tg_xxyyy_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyy_1[j];

                    tg_xxxxyyy_xxyz_0[j] = pb_x * tg_xxxyyy_xxyz_0[j] + fr * tg_xxxyyy_xxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxyz_0[j] - tg_xxyyy_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xyz_1[j];

                    tg_xxxxyyy_xxzz_0[j] = pb_x * tg_xxxyyy_xxzz_0[j] + fr * tg_xxxyyy_xxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxzz_0[j] - tg_xxyyy_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xzz_1[j];

                    tg_xxxxyyy_xyyy_0[j] = pb_x * tg_xxxyyy_xyyy_0[j] + fr * tg_xxxyyy_xyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyy_0[j] - tg_xxyyy_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyy_1[j];

                    tg_xxxxyyy_xyyz_0[j] = pb_x * tg_xxxyyy_xyyz_0[j] + fr * tg_xxxyyy_xyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyyz_0[j] - tg_xxyyy_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yyz_1[j];

                    tg_xxxxyyy_xyzz_0[j] = pb_x * tg_xxxyyy_xyzz_0[j] + fr * tg_xxxyyy_xyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyzz_0[j] - tg_xxyyy_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yzz_1[j];

                    tg_xxxxyyy_xzzz_0[j] = pb_x * tg_xxxyyy_xzzz_0[j] + fr * tg_xxxyyy_xzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xzzz_0[j] - tg_xxyyy_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_zzz_1[j];

                    tg_xxxxyyy_yyyy_0[j] = pb_x * tg_xxxyyy_yyyy_0[j] + fr * tg_xxxyyy_yyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyy_0[j] - tg_xxyyy_yyyy_1[j] * fl1_fza);

                    tg_xxxxyyy_yyyz_0[j] = pb_x * tg_xxxyyy_yyyz_0[j] + fr * tg_xxxyyy_yyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyyz_0[j] - tg_xxyyy_yyyz_1[j] * fl1_fza);

                    tg_xxxxyyy_yyzz_0[j] = pb_x * tg_xxxyyy_yyzz_0[j] + fr * tg_xxxyyy_yyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyzz_0[j] - tg_xxyyy_yyzz_1[j] * fl1_fza);

                    tg_xxxxyyy_yzzz_0[j] = pb_x * tg_xxxyyy_yzzz_0[j] + fr * tg_xxxyyy_yzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yzzz_0[j] - tg_xxyyy_yzzz_1[j] * fl1_fza);

                    tg_xxxxyyy_zzzz_0[j] = pb_x * tg_xxxyyy_zzzz_0[j] + fr * tg_xxxyyy_zzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_zzzz_0[j] - tg_xxyyy_zzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_xxxx_0[j] = pb_x * tg_xxxyyz_xxxx_0[j] + fr * tg_xxxyyz_xxxx_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxx_0[j] - tg_xxyyz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyyz_xxx_1[j];

                    tg_xxxxyyz_xxxy_0[j] = pb_x * tg_xxxyyz_xxxy_0[j] + fr * tg_xxxyyz_xxxy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxy_0[j] - tg_xxyyz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxy_1[j];

                    tg_xxxxyyz_xxxz_0[j] = pb_x * tg_xxxyyz_xxxz_0[j] + fr * tg_xxxyyz_xxxz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxxz_0[j] - tg_xxyyz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xxz_1[j];

                    tg_xxxxyyz_xxyy_0[j] = pb_x * tg_xxxyyz_xxyy_0[j] + fr * tg_xxxyyz_xxyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyy_0[j] - tg_xxyyz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyy_1[j];

                    tg_xxxxyyz_xxyz_0[j] = pb_x * tg_xxxyyz_xxyz_0[j] + fr * tg_xxxyyz_xxyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxyz_0[j] - tg_xxyyz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xyz_1[j];

                    tg_xxxxyyz_xxzz_0[j] = pb_x * tg_xxxyyz_xxzz_0[j] + fr * tg_xxxyyz_xxzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxzz_0[j] - tg_xxyyz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xzz_1[j];

                    tg_xxxxyyz_xyyy_0[j] = pb_x * tg_xxxyyz_xyyy_0[j] + fr * tg_xxxyyz_xyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyy_0[j] - tg_xxyyz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyy_1[j];

                    tg_xxxxyyz_xyyz_0[j] = pb_x * tg_xxxyyz_xyyz_0[j] + fr * tg_xxxyyz_xyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyyz_0[j] - tg_xxyyz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yyz_1[j];

                    tg_xxxxyyz_xyzz_0[j] = pb_x * tg_xxxyyz_xyzz_0[j] + fr * tg_xxxyyz_xyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyzz_0[j] - tg_xxyyz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yzz_1[j];

                    tg_xxxxyyz_xzzz_0[j] = pb_x * tg_xxxyyz_xzzz_0[j] + fr * tg_xxxyyz_xzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xzzz_0[j] - tg_xxyyz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_zzz_1[j];

                    tg_xxxxyyz_yyyy_0[j] = pb_x * tg_xxxyyz_yyyy_0[j] + fr * tg_xxxyyz_yyyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyy_0[j] - tg_xxyyz_yyyy_1[j] * fl1_fza);

                    tg_xxxxyyz_yyyz_0[j] = pb_x * tg_xxxyyz_yyyz_0[j] + fr * tg_xxxyyz_yyyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyyz_0[j] - tg_xxyyz_yyyz_1[j] * fl1_fza);

                    tg_xxxxyyz_yyzz_0[j] = pb_x * tg_xxxyyz_yyzz_0[j] + fr * tg_xxxyyz_yyzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyzz_0[j] - tg_xxyyz_yyzz_1[j] * fl1_fza);

                    tg_xxxxyyz_yzzz_0[j] = pb_x * tg_xxxyyz_yzzz_0[j] + fr * tg_xxxyyz_yzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yzzz_0[j] - tg_xxyyz_yzzz_1[j] * fl1_fza);

                    tg_xxxxyyz_zzzz_0[j] = pb_x * tg_xxxyyz_zzzz_0[j] + fr * tg_xxxyyz_zzzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_zzzz_0[j] - tg_xxyyz_zzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_xxxx_0[j] = pb_x * tg_xxxyzz_xxxx_0[j] + fr * tg_xxxyzz_xxxx_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxx_0[j] - tg_xxyzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxyzz_xxx_1[j];

                    tg_xxxxyzz_xxxy_0[j] = pb_x * tg_xxxyzz_xxxy_0[j] + fr * tg_xxxyzz_xxxy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxy_0[j] - tg_xxyzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxy_1[j];

                    tg_xxxxyzz_xxxz_0[j] = pb_x * tg_xxxyzz_xxxz_0[j] + fr * tg_xxxyzz_xxxz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxxz_0[j] - tg_xxyzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xxz_1[j];

                    tg_xxxxyzz_xxyy_0[j] = pb_x * tg_xxxyzz_xxyy_0[j] + fr * tg_xxxyzz_xxyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyy_0[j] - tg_xxyzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyy_1[j];

                    tg_xxxxyzz_xxyz_0[j] = pb_x * tg_xxxyzz_xxyz_0[j] + fr * tg_xxxyzz_xxyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxyz_0[j] - tg_xxyzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xyz_1[j];

                    tg_xxxxyzz_xxzz_0[j] = pb_x * tg_xxxyzz_xxzz_0[j] + fr * tg_xxxyzz_xxzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxzz_0[j] - tg_xxyzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xzz_1[j];

                    tg_xxxxyzz_xyyy_0[j] = pb_x * tg_xxxyzz_xyyy_0[j] + fr * tg_xxxyzz_xyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyy_0[j] - tg_xxyzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyy_1[j];

                    tg_xxxxyzz_xyyz_0[j] = pb_x * tg_xxxyzz_xyyz_0[j] + fr * tg_xxxyzz_xyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyyz_0[j] - tg_xxyzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yyz_1[j];

                    tg_xxxxyzz_xyzz_0[j] = pb_x * tg_xxxyzz_xyzz_0[j] + fr * tg_xxxyzz_xyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyzz_0[j] - tg_xxyzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yzz_1[j];

                    tg_xxxxyzz_xzzz_0[j] = pb_x * tg_xxxyzz_xzzz_0[j] + fr * tg_xxxyzz_xzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xzzz_0[j] - tg_xxyzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_zzz_1[j];

                    tg_xxxxyzz_yyyy_0[j] = pb_x * tg_xxxyzz_yyyy_0[j] + fr * tg_xxxyzz_yyyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyy_0[j] - tg_xxyzz_yyyy_1[j] * fl1_fza);

                    tg_xxxxyzz_yyyz_0[j] = pb_x * tg_xxxyzz_yyyz_0[j] + fr * tg_xxxyzz_yyyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyyz_0[j] - tg_xxyzz_yyyz_1[j] * fl1_fza);

                    tg_xxxxyzz_yyzz_0[j] = pb_x * tg_xxxyzz_yyzz_0[j] + fr * tg_xxxyzz_yyzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyzz_0[j] - tg_xxyzz_yyzz_1[j] * fl1_fza);

                    tg_xxxxyzz_yzzz_0[j] = pb_x * tg_xxxyzz_yzzz_0[j] + fr * tg_xxxyzz_yzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yzzz_0[j] - tg_xxyzz_yzzz_1[j] * fl1_fza);

                    tg_xxxxyzz_zzzz_0[j] = pb_x * tg_xxxyzz_zzzz_0[j] + fr * tg_xxxyzz_zzzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_zzzz_0[j] - tg_xxyzz_zzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_xxxx_0[j] = pb_x * tg_xxxzzz_xxxx_0[j] + fr * tg_xxxzzz_xxxx_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxx_0[j] - tg_xxzzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxxzzz_xxx_1[j];

                    tg_xxxxzzz_xxxy_0[j] = pb_x * tg_xxxzzz_xxxy_0[j] + fr * tg_xxxzzz_xxxy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxy_0[j] - tg_xxzzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxy_1[j];

                    tg_xxxxzzz_xxxz_0[j] = pb_x * tg_xxxzzz_xxxz_0[j] + fr * tg_xxxzzz_xxxz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxxz_0[j] - tg_xxzzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xxz_1[j];

                    tg_xxxxzzz_xxyy_0[j] = pb_x * tg_xxxzzz_xxyy_0[j] + fr * tg_xxxzzz_xxyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyy_0[j] - tg_xxzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyy_1[j];

                    tg_xxxxzzz_xxyz_0[j] = pb_x * tg_xxxzzz_xxyz_0[j] + fr * tg_xxxzzz_xxyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxyz_0[j] - tg_xxzzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xyz_1[j];

                    tg_xxxxzzz_xxzz_0[j] = pb_x * tg_xxxzzz_xxzz_0[j] + fr * tg_xxxzzz_xxzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxzz_0[j] - tg_xxzzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xzz_1[j];

                    tg_xxxxzzz_xyyy_0[j] = pb_x * tg_xxxzzz_xyyy_0[j] + fr * tg_xxxzzz_xyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyy_0[j] - tg_xxzzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyy_1[j];

                    tg_xxxxzzz_xyyz_0[j] = pb_x * tg_xxxzzz_xyyz_0[j] + fr * tg_xxxzzz_xyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyyz_0[j] - tg_xxzzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yyz_1[j];

                    tg_xxxxzzz_xyzz_0[j] = pb_x * tg_xxxzzz_xyzz_0[j] + fr * tg_xxxzzz_xyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyzz_0[j] - tg_xxzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yzz_1[j];

                    tg_xxxxzzz_xzzz_0[j] = pb_x * tg_xxxzzz_xzzz_0[j] + fr * tg_xxxzzz_xzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xzzz_0[j] - tg_xxzzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_zzz_1[j];

                    tg_xxxxzzz_yyyy_0[j] = pb_x * tg_xxxzzz_yyyy_0[j] + fr * tg_xxxzzz_yyyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyy_0[j] - tg_xxzzz_yyyy_1[j] * fl1_fza);

                    tg_xxxxzzz_yyyz_0[j] = pb_x * tg_xxxzzz_yyyz_0[j] + fr * tg_xxxzzz_yyyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyyz_0[j] - tg_xxzzz_yyyz_1[j] * fl1_fza);

                    tg_xxxxzzz_yyzz_0[j] = pb_x * tg_xxxzzz_yyzz_0[j] + fr * tg_xxxzzz_yyzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyzz_0[j] - tg_xxzzz_yyzz_1[j] * fl1_fza);

                    tg_xxxxzzz_yzzz_0[j] = pb_x * tg_xxxzzz_yzzz_0[j] + fr * tg_xxxzzz_yzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yzzz_0[j] - tg_xxzzz_yzzz_1[j] * fl1_fza);

                    tg_xxxxzzz_zzzz_0[j] = pb_x * tg_xxxzzz_zzzz_0[j] + fr * tg_xxxzzz_zzzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_zzzz_0[j] - tg_xxzzz_zzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_xxxx_0[j] = pb_x * tg_xxyyyy_xxxx_0[j] + fr * tg_xxyyyy_xxxx_1[j] + fl1_fx * (tg_xyyyy_xxxx_0[j] - tg_xyyyy_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyy_xxx_1[j];

                    tg_xxxyyyy_xxxy_0[j] = pb_x * tg_xxyyyy_xxxy_0[j] + fr * tg_xxyyyy_xxxy_1[j] + fl1_fx * (tg_xyyyy_xxxy_0[j] - tg_xyyyy_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxy_1[j];

                    tg_xxxyyyy_xxxz_0[j] = pb_x * tg_xxyyyy_xxxz_0[j] + fr * tg_xxyyyy_xxxz_1[j] + fl1_fx * (tg_xyyyy_xxxz_0[j] - tg_xyyyy_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xxz_1[j];

                    tg_xxxyyyy_xxyy_0[j] = pb_x * tg_xxyyyy_xxyy_0[j] + fr * tg_xxyyyy_xxyy_1[j] + fl1_fx * (tg_xyyyy_xxyy_0[j] - tg_xyyyy_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyy_1[j];

                    tg_xxxyyyy_xxyz_0[j] = pb_x * tg_xxyyyy_xxyz_0[j] + fr * tg_xxyyyy_xxyz_1[j] + fl1_fx * (tg_xyyyy_xxyz_0[j] - tg_xyyyy_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xyz_1[j];

                    tg_xxxyyyy_xxzz_0[j] = pb_x * tg_xxyyyy_xxzz_0[j] + fr * tg_xxyyyy_xxzz_1[j] + fl1_fx * (tg_xyyyy_xxzz_0[j] - tg_xyyyy_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xzz_1[j];

                    tg_xxxyyyy_xyyy_0[j] = pb_x * tg_xxyyyy_xyyy_0[j] + fr * tg_xxyyyy_xyyy_1[j] + fl1_fx * (tg_xyyyy_xyyy_0[j] - tg_xyyyy_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyy_1[j];

                    tg_xxxyyyy_xyyz_0[j] = pb_x * tg_xxyyyy_xyyz_0[j] + fr * tg_xxyyyy_xyyz_1[j] + fl1_fx * (tg_xyyyy_xyyz_0[j] - tg_xyyyy_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yyz_1[j];

                    tg_xxxyyyy_xyzz_0[j] = pb_x * tg_xxyyyy_xyzz_0[j] + fr * tg_xxyyyy_xyzz_1[j] + fl1_fx * (tg_xyyyy_xyzz_0[j] - tg_xyyyy_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yzz_1[j];

                    tg_xxxyyyy_xzzz_0[j] = pb_x * tg_xxyyyy_xzzz_0[j] + fr * tg_xxyyyy_xzzz_1[j] + fl1_fx * (tg_xyyyy_xzzz_0[j] - tg_xyyyy_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_zzz_1[j];

                    tg_xxxyyyy_yyyy_0[j] = pb_x * tg_xxyyyy_yyyy_0[j] + fr * tg_xxyyyy_yyyy_1[j] + fl1_fx * (tg_xyyyy_yyyy_0[j] - tg_xyyyy_yyyy_1[j] * fl1_fza);

                    tg_xxxyyyy_yyyz_0[j] = pb_x * tg_xxyyyy_yyyz_0[j] + fr * tg_xxyyyy_yyyz_1[j] + fl1_fx * (tg_xyyyy_yyyz_0[j] - tg_xyyyy_yyyz_1[j] * fl1_fza);

                    tg_xxxyyyy_yyzz_0[j] = pb_x * tg_xxyyyy_yyzz_0[j] + fr * tg_xxyyyy_yyzz_1[j] + fl1_fx * (tg_xyyyy_yyzz_0[j] - tg_xyyyy_yyzz_1[j] * fl1_fza);

                    tg_xxxyyyy_yzzz_0[j] = pb_x * tg_xxyyyy_yzzz_0[j] + fr * tg_xxyyyy_yzzz_1[j] + fl1_fx * (tg_xyyyy_yzzz_0[j] - tg_xyyyy_yzzz_1[j] * fl1_fza);

                    tg_xxxyyyy_zzzz_0[j] = pb_x * tg_xxyyyy_zzzz_0[j] + fr * tg_xxyyyy_zzzz_1[j] + fl1_fx * (tg_xyyyy_zzzz_0[j] - tg_xyyyy_zzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_xxxx_0[j] = pb_x * tg_xxyyyz_xxxx_0[j] + fr * tg_xxyyyz_xxxx_1[j] + fl1_fx * (tg_xyyyz_xxxx_0[j] - tg_xyyyz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyyz_xxx_1[j];

                    tg_xxxyyyz_xxxy_0[j] = pb_x * tg_xxyyyz_xxxy_0[j] + fr * tg_xxyyyz_xxxy_1[j] + fl1_fx * (tg_xyyyz_xxxy_0[j] - tg_xyyyz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxy_1[j];

                    tg_xxxyyyz_xxxz_0[j] = pb_x * tg_xxyyyz_xxxz_0[j] + fr * tg_xxyyyz_xxxz_1[j] + fl1_fx * (tg_xyyyz_xxxz_0[j] - tg_xyyyz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xxz_1[j];

                    tg_xxxyyyz_xxyy_0[j] = pb_x * tg_xxyyyz_xxyy_0[j] + fr * tg_xxyyyz_xxyy_1[j] + fl1_fx * (tg_xyyyz_xxyy_0[j] - tg_xyyyz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyy_1[j];

                    tg_xxxyyyz_xxyz_0[j] = pb_x * tg_xxyyyz_xxyz_0[j] + fr * tg_xxyyyz_xxyz_1[j] + fl1_fx * (tg_xyyyz_xxyz_0[j] - tg_xyyyz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xyz_1[j];

                    tg_xxxyyyz_xxzz_0[j] = pb_x * tg_xxyyyz_xxzz_0[j] + fr * tg_xxyyyz_xxzz_1[j] + fl1_fx * (tg_xyyyz_xxzz_0[j] - tg_xyyyz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xzz_1[j];

                    tg_xxxyyyz_xyyy_0[j] = pb_x * tg_xxyyyz_xyyy_0[j] + fr * tg_xxyyyz_xyyy_1[j] + fl1_fx * (tg_xyyyz_xyyy_0[j] - tg_xyyyz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyy_1[j];

                    tg_xxxyyyz_xyyz_0[j] = pb_x * tg_xxyyyz_xyyz_0[j] + fr * tg_xxyyyz_xyyz_1[j] + fl1_fx * (tg_xyyyz_xyyz_0[j] - tg_xyyyz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yyz_1[j];

                    tg_xxxyyyz_xyzz_0[j] = pb_x * tg_xxyyyz_xyzz_0[j] + fr * tg_xxyyyz_xyzz_1[j] + fl1_fx * (tg_xyyyz_xyzz_0[j] - tg_xyyyz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yzz_1[j];

                    tg_xxxyyyz_xzzz_0[j] = pb_x * tg_xxyyyz_xzzz_0[j] + fr * tg_xxyyyz_xzzz_1[j] + fl1_fx * (tg_xyyyz_xzzz_0[j] - tg_xyyyz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_zzz_1[j];

                    tg_xxxyyyz_yyyy_0[j] = pb_x * tg_xxyyyz_yyyy_0[j] + fr * tg_xxyyyz_yyyy_1[j] + fl1_fx * (tg_xyyyz_yyyy_0[j] - tg_xyyyz_yyyy_1[j] * fl1_fza);

                    tg_xxxyyyz_yyyz_0[j] = pb_x * tg_xxyyyz_yyyz_0[j] + fr * tg_xxyyyz_yyyz_1[j] + fl1_fx * (tg_xyyyz_yyyz_0[j] - tg_xyyyz_yyyz_1[j] * fl1_fza);

                    tg_xxxyyyz_yyzz_0[j] = pb_x * tg_xxyyyz_yyzz_0[j] + fr * tg_xxyyyz_yyzz_1[j] + fl1_fx * (tg_xyyyz_yyzz_0[j] - tg_xyyyz_yyzz_1[j] * fl1_fza);

                    tg_xxxyyyz_yzzz_0[j] = pb_x * tg_xxyyyz_yzzz_0[j] + fr * tg_xxyyyz_yzzz_1[j] + fl1_fx * (tg_xyyyz_yzzz_0[j] - tg_xyyyz_yzzz_1[j] * fl1_fza);

                    tg_xxxyyyz_zzzz_0[j] = pb_x * tg_xxyyyz_zzzz_0[j] + fr * tg_xxyyyz_zzzz_1[j] + fl1_fx * (tg_xyyyz_zzzz_0[j] - tg_xyyyz_zzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSG_180_270(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (180,270)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxyyzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 180); 

                auto tg_xxyyzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 181); 

                auto tg_xxyyzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 182); 

                auto tg_xxyyzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 183); 

                auto tg_xxyyzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 184); 

                auto tg_xxyyzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 185); 

                auto tg_xxyyzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 186); 

                auto tg_xxyyzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 187); 

                auto tg_xxyyzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 188); 

                auto tg_xxyyzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 189); 

                auto tg_xxyyzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 190); 

                auto tg_xxyyzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 191); 

                auto tg_xxyyzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 192); 

                auto tg_xxyyzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 193); 

                auto tg_xxyyzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 194); 

                auto tg_xxyzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 195); 

                auto tg_xxyzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 196); 

                auto tg_xxyzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 197); 

                auto tg_xxyzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 198); 

                auto tg_xxyzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 199); 

                auto tg_xxyzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 200); 

                auto tg_xxyzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 201); 

                auto tg_xxyzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 202); 

                auto tg_xxyzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 203); 

                auto tg_xxyzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 204); 

                auto tg_xxyzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 205); 

                auto tg_xxyzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 206); 

                auto tg_xxyzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 207); 

                auto tg_xxyzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 208); 

                auto tg_xxyzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 209); 

                auto tg_xxzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 210); 

                auto tg_xxzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 211); 

                auto tg_xxzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 212); 

                auto tg_xxzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 213); 

                auto tg_xxzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 214); 

                auto tg_xxzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 215); 

                auto tg_xxzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 216); 

                auto tg_xxzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 217); 

                auto tg_xxzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 218); 

                auto tg_xxzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 219); 

                auto tg_xxzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 220); 

                auto tg_xxzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 221); 

                auto tg_xxzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 222); 

                auto tg_xxzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 223); 

                auto tg_xxzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 224); 

                auto tg_xyyyyy_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 225); 

                auto tg_xyyyyy_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 226); 

                auto tg_xyyyyy_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 227); 

                auto tg_xyyyyy_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 228); 

                auto tg_xyyyyy_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 229); 

                auto tg_xyyyyy_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 230); 

                auto tg_xyyyyy_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 231); 

                auto tg_xyyyyy_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 232); 

                auto tg_xyyyyy_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 233); 

                auto tg_xyyyyy_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 234); 

                auto tg_xyyyyy_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 235); 

                auto tg_xyyyyy_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 236); 

                auto tg_xyyyyy_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 237); 

                auto tg_xyyyyy_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 238); 

                auto tg_xyyyyy_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 239); 

                auto tg_xyyyyz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 240); 

                auto tg_xyyyyz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 241); 

                auto tg_xyyyyz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 242); 

                auto tg_xyyyyz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 243); 

                auto tg_xyyyyz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 244); 

                auto tg_xyyyyz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 245); 

                auto tg_xyyyyz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 246); 

                auto tg_xyyyyz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 247); 

                auto tg_xyyyyz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 248); 

                auto tg_xyyyyz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 249); 

                auto tg_xyyyyz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 250); 

                auto tg_xyyyyz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 251); 

                auto tg_xyyyyz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 252); 

                auto tg_xyyyyz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 253); 

                auto tg_xyyyyz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 254); 

                auto tg_xyyyzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 255); 

                auto tg_xyyyzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 256); 

                auto tg_xyyyzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 257); 

                auto tg_xyyyzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 258); 

                auto tg_xyyyzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 259); 

                auto tg_xyyyzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 260); 

                auto tg_xyyyzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 261); 

                auto tg_xyyyzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 262); 

                auto tg_xyyyzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 263); 

                auto tg_xyyyzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 264); 

                auto tg_xyyyzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 265); 

                auto tg_xyyyzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 266); 

                auto tg_xyyyzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 267); 

                auto tg_xyyyzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 268); 

                auto tg_xyyyzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 269); 

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

                auto tg_xyyzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 180); 

                auto tg_xyyzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 181); 

                auto tg_xyyzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 182); 

                auto tg_xyyzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 183); 

                auto tg_xyyzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 184); 

                auto tg_xyyzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 185); 

                auto tg_xyyzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 186); 

                auto tg_xyyzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 187); 

                auto tg_xyyzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 188); 

                auto tg_xyyzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 189); 

                auto tg_xyyzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 190); 

                auto tg_xyyzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 191); 

                auto tg_xyyzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 192); 

                auto tg_xyyzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 193); 

                auto tg_xyyzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 194); 

                auto tg_xyzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 195); 

                auto tg_xyzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 196); 

                auto tg_xyzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 197); 

                auto tg_xyzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 198); 

                auto tg_xyzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 199); 

                auto tg_xyzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 200); 

                auto tg_xyzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 201); 

                auto tg_xyzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 202); 

                auto tg_xyzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 203); 

                auto tg_xyzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 204); 

                auto tg_xyzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 205); 

                auto tg_xyzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 206); 

                auto tg_xyzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 207); 

                auto tg_xyzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 208); 

                auto tg_xyzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 209); 

                auto tg_xzzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 210); 

                auto tg_xzzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 211); 

                auto tg_xzzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 212); 

                auto tg_xzzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 213); 

                auto tg_xzzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 214); 

                auto tg_xzzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 215); 

                auto tg_xzzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 216); 

                auto tg_xzzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 217); 

                auto tg_xzzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 218); 

                auto tg_xzzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 219); 

                auto tg_xzzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 220); 

                auto tg_xzzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 221); 

                auto tg_xzzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 222); 

                auto tg_xzzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 223); 

                auto tg_xzzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 224); 

                auto tg_yyyyy_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 225); 

                auto tg_yyyyy_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 226); 

                auto tg_yyyyy_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 227); 

                auto tg_yyyyy_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 228); 

                auto tg_yyyyy_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 229); 

                auto tg_yyyyy_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 230); 

                auto tg_yyyyy_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 231); 

                auto tg_yyyyy_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 232); 

                auto tg_yyyyy_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 233); 

                auto tg_yyyyy_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 234); 

                auto tg_yyyyy_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 235); 

                auto tg_yyyyy_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 236); 

                auto tg_yyyyy_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 237); 

                auto tg_yyyyy_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 238); 

                auto tg_yyyyy_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 239); 

                auto tg_yyyyz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 240); 

                auto tg_yyyyz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 241); 

                auto tg_yyyyz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 242); 

                auto tg_yyyyz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 243); 

                auto tg_yyyyz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 244); 

                auto tg_yyyyz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 245); 

                auto tg_yyyyz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 246); 

                auto tg_yyyyz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 247); 

                auto tg_yyyyz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 248); 

                auto tg_yyyyz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 249); 

                auto tg_yyyyz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 250); 

                auto tg_yyyyz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 251); 

                auto tg_yyyyz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 252); 

                auto tg_yyyyz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 253); 

                auto tg_yyyyz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 254); 

                auto tg_yyyzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 255); 

                auto tg_yyyzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 256); 

                auto tg_yyyzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 257); 

                auto tg_yyyzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 258); 

                auto tg_yyyzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 259); 

                auto tg_yyyzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 260); 

                auto tg_yyyzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 261); 

                auto tg_yyyzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 262); 

                auto tg_yyyzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 263); 

                auto tg_yyyzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 264); 

                auto tg_yyyzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 265); 

                auto tg_yyyzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 266); 

                auto tg_yyyzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 267); 

                auto tg_yyyzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 268); 

                auto tg_yyyzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 269); 

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

                auto tg_xxyyzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 120); 

                auto tg_xxyyzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 121); 

                auto tg_xxyyzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 122); 

                auto tg_xxyyzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 123); 

                auto tg_xxyyzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 124); 

                auto tg_xxyyzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 125); 

                auto tg_xxyyzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 126); 

                auto tg_xxyyzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 127); 

                auto tg_xxyyzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 128); 

                auto tg_xxyyzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 129); 

                auto tg_xxyzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 130); 

                auto tg_xxyzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 131); 

                auto tg_xxyzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 132); 

                auto tg_xxyzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 133); 

                auto tg_xxyzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 134); 

                auto tg_xxyzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 135); 

                auto tg_xxyzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 136); 

                auto tg_xxyzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 137); 

                auto tg_xxyzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 138); 

                auto tg_xxyzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 139); 

                auto tg_xxzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 140); 

                auto tg_xxzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 141); 

                auto tg_xxzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 142); 

                auto tg_xxzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 143); 

                auto tg_xxzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 144); 

                auto tg_xxzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 145); 

                auto tg_xxzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 146); 

                auto tg_xxzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 147); 

                auto tg_xxzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 148); 

                auto tg_xxzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 149); 

                auto tg_xyyyyy_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 150); 

                auto tg_xyyyyy_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 151); 

                auto tg_xyyyyy_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 152); 

                auto tg_xyyyyy_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 153); 

                auto tg_xyyyyy_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 154); 

                auto tg_xyyyyy_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 155); 

                auto tg_xyyyyy_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 156); 

                auto tg_xyyyyy_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 157); 

                auto tg_xyyyyy_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 158); 

                auto tg_xyyyyy_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 159); 

                auto tg_xyyyyz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 160); 

                auto tg_xyyyyz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 161); 

                auto tg_xyyyyz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 162); 

                auto tg_xyyyyz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 163); 

                auto tg_xyyyyz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 164); 

                auto tg_xyyyyz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 165); 

                auto tg_xyyyyz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 166); 

                auto tg_xyyyyz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 167); 

                auto tg_xyyyyz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 168); 

                auto tg_xyyyyz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 169); 

                auto tg_xyyyzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 170); 

                auto tg_xyyyzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 171); 

                auto tg_xyyyzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 172); 

                auto tg_xyyyzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 173); 

                auto tg_xyyyzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 174); 

                auto tg_xyyyzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 175); 

                auto tg_xyyyzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 176); 

                auto tg_xyyyzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 177); 

                auto tg_xyyyzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 178); 

                auto tg_xyyyzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 179); 

                // set up pointers to integrals

                auto tg_xxxyyzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 180); 

                auto tg_xxxyyzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 181); 

                auto tg_xxxyyzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 182); 

                auto tg_xxxyyzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 183); 

                auto tg_xxxyyzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 184); 

                auto tg_xxxyyzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 185); 

                auto tg_xxxyyzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 186); 

                auto tg_xxxyyzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 187); 

                auto tg_xxxyyzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 188); 

                auto tg_xxxyyzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 189); 

                auto tg_xxxyyzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 190); 

                auto tg_xxxyyzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 191); 

                auto tg_xxxyyzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 192); 

                auto tg_xxxyyzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 193); 

                auto tg_xxxyyzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 194); 

                auto tg_xxxyzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 195); 

                auto tg_xxxyzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 196); 

                auto tg_xxxyzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 197); 

                auto tg_xxxyzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 198); 

                auto tg_xxxyzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 199); 

                auto tg_xxxyzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 200); 

                auto tg_xxxyzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 201); 

                auto tg_xxxyzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 202); 

                auto tg_xxxyzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 203); 

                auto tg_xxxyzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 204); 

                auto tg_xxxyzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 205); 

                auto tg_xxxyzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 206); 

                auto tg_xxxyzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 207); 

                auto tg_xxxyzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 208); 

                auto tg_xxxyzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 209); 

                auto tg_xxxzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 210); 

                auto tg_xxxzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 211); 

                auto tg_xxxzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 212); 

                auto tg_xxxzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 213); 

                auto tg_xxxzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 214); 

                auto tg_xxxzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 215); 

                auto tg_xxxzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 216); 

                auto tg_xxxzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 217); 

                auto tg_xxxzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 218); 

                auto tg_xxxzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 219); 

                auto tg_xxxzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 220); 

                auto tg_xxxzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 221); 

                auto tg_xxxzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 222); 

                auto tg_xxxzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 223); 

                auto tg_xxxzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 224); 

                auto tg_xxyyyyy_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 225); 

                auto tg_xxyyyyy_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 226); 

                auto tg_xxyyyyy_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 227); 

                auto tg_xxyyyyy_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 228); 

                auto tg_xxyyyyy_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 229); 

                auto tg_xxyyyyy_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 230); 

                auto tg_xxyyyyy_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 231); 

                auto tg_xxyyyyy_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 232); 

                auto tg_xxyyyyy_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 233); 

                auto tg_xxyyyyy_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 234); 

                auto tg_xxyyyyy_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 235); 

                auto tg_xxyyyyy_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 236); 

                auto tg_xxyyyyy_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 237); 

                auto tg_xxyyyyy_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 238); 

                auto tg_xxyyyyy_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 239); 

                auto tg_xxyyyyz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 240); 

                auto tg_xxyyyyz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 241); 

                auto tg_xxyyyyz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 242); 

                auto tg_xxyyyyz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 243); 

                auto tg_xxyyyyz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 244); 

                auto tg_xxyyyyz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 245); 

                auto tg_xxyyyyz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 246); 

                auto tg_xxyyyyz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 247); 

                auto tg_xxyyyyz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 248); 

                auto tg_xxyyyyz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 249); 

                auto tg_xxyyyyz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 250); 

                auto tg_xxyyyyz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 251); 

                auto tg_xxyyyyz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 252); 

                auto tg_xxyyyyz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 253); 

                auto tg_xxyyyyz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 254); 

                auto tg_xxyyyzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 255); 

                auto tg_xxyyyzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 256); 

                auto tg_xxyyyzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 257); 

                auto tg_xxyyyzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 258); 

                auto tg_xxyyyzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 259); 

                auto tg_xxyyyzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 260); 

                auto tg_xxyyyzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 261); 

                auto tg_xxyyyzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 262); 

                auto tg_xxyyyzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 263); 

                auto tg_xxyyyzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 264); 

                auto tg_xxyyyzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 265); 

                auto tg_xxyyyzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 266); 

                auto tg_xxyyyzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 267); 

                auto tg_xxyyyzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 268); 

                auto tg_xxyyyzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 269); 

                // Batch of Integrals (180,270)

                #pragma omp simd aligned(fxn, fza, tg_xxxyyzz_xxxx_0, tg_xxxyyzz_xxxy_0, tg_xxxyyzz_xxxz_0, \
                                         tg_xxxyyzz_xxyy_0, tg_xxxyyzz_xxyz_0, tg_xxxyyzz_xxzz_0, tg_xxxyyzz_xyyy_0, \
                                         tg_xxxyyzz_xyyz_0, tg_xxxyyzz_xyzz_0, tg_xxxyyzz_xzzz_0, tg_xxxyyzz_yyyy_0, \
                                         tg_xxxyyzz_yyyz_0, tg_xxxyyzz_yyzz_0, tg_xxxyyzz_yzzz_0, tg_xxxyyzz_zzzz_0, \
                                         tg_xxxyzzz_xxxx_0, tg_xxxyzzz_xxxy_0, tg_xxxyzzz_xxxz_0, tg_xxxyzzz_xxyy_0, \
                                         tg_xxxyzzz_xxyz_0, tg_xxxyzzz_xxzz_0, tg_xxxyzzz_xyyy_0, tg_xxxyzzz_xyyz_0, \
                                         tg_xxxyzzz_xyzz_0, tg_xxxyzzz_xzzz_0, tg_xxxyzzz_yyyy_0, tg_xxxyzzz_yyyz_0, \
                                         tg_xxxyzzz_yyzz_0, tg_xxxyzzz_yzzz_0, tg_xxxyzzz_zzzz_0, tg_xxxzzzz_xxxx_0, \
                                         tg_xxxzzzz_xxxy_0, tg_xxxzzzz_xxxz_0, tg_xxxzzzz_xxyy_0, tg_xxxzzzz_xxyz_0, \
                                         tg_xxxzzzz_xxzz_0, tg_xxxzzzz_xyyy_0, tg_xxxzzzz_xyyz_0, tg_xxxzzzz_xyzz_0, \
                                         tg_xxxzzzz_xzzz_0, tg_xxxzzzz_yyyy_0, tg_xxxzzzz_yyyz_0, tg_xxxzzzz_yyzz_0, \
                                         tg_xxxzzzz_yzzz_0, tg_xxxzzzz_zzzz_0, tg_xxyyyyy_xxxx_0, tg_xxyyyyy_xxxy_0, \
                                         tg_xxyyyyy_xxxz_0, tg_xxyyyyy_xxyy_0, tg_xxyyyyy_xxyz_0, tg_xxyyyyy_xxzz_0, \
                                         tg_xxyyyyy_xyyy_0, tg_xxyyyyy_xyyz_0, tg_xxyyyyy_xyzz_0, tg_xxyyyyy_xzzz_0, \
                                         tg_xxyyyyy_yyyy_0, tg_xxyyyyy_yyyz_0, tg_xxyyyyy_yyzz_0, tg_xxyyyyy_yzzz_0, \
                                         tg_xxyyyyy_zzzz_0, tg_xxyyyyz_xxxx_0, tg_xxyyyyz_xxxy_0, tg_xxyyyyz_xxxz_0, \
                                         tg_xxyyyyz_xxyy_0, tg_xxyyyyz_xxyz_0, tg_xxyyyyz_xxzz_0, tg_xxyyyyz_xyyy_0, \
                                         tg_xxyyyyz_xyyz_0, tg_xxyyyyz_xyzz_0, tg_xxyyyyz_xzzz_0, tg_xxyyyyz_yyyy_0, \
                                         tg_xxyyyyz_yyyz_0, tg_xxyyyyz_yyzz_0, tg_xxyyyyz_yzzz_0, tg_xxyyyyz_zzzz_0, \
                                         tg_xxyyyzz_xxxx_0, tg_xxyyyzz_xxxy_0, tg_xxyyyzz_xxxz_0, tg_xxyyyzz_xxyy_0, \
                                         tg_xxyyyzz_xxyz_0, tg_xxyyyzz_xxzz_0, tg_xxyyyzz_xyyy_0, tg_xxyyyzz_xyyz_0, \
                                         tg_xxyyyzz_xyzz_0, tg_xxyyyzz_xzzz_0, tg_xxyyyzz_yyyy_0, tg_xxyyyzz_yyyz_0, \
                                         tg_xxyyyzz_yyzz_0, tg_xxyyyzz_yzzz_0, tg_xxyyyzz_zzzz_0, tg_xxyyzz_xxx_1, \
                                         tg_xxyyzz_xxxx_0, tg_xxyyzz_xxxx_1, tg_xxyyzz_xxxy_0, tg_xxyyzz_xxxy_1, \
                                         tg_xxyyzz_xxxz_0, tg_xxyyzz_xxxz_1, tg_xxyyzz_xxy_1, tg_xxyyzz_xxyy_0, \
                                         tg_xxyyzz_xxyy_1, tg_xxyyzz_xxyz_0, tg_xxyyzz_xxyz_1, tg_xxyyzz_xxz_1, \
                                         tg_xxyyzz_xxzz_0, tg_xxyyzz_xxzz_1, tg_xxyyzz_xyy_1, tg_xxyyzz_xyyy_0, \
                                         tg_xxyyzz_xyyy_1, tg_xxyyzz_xyyz_0, tg_xxyyzz_xyyz_1, tg_xxyyzz_xyz_1, \
                                         tg_xxyyzz_xyzz_0, tg_xxyyzz_xyzz_1, tg_xxyyzz_xzz_1, tg_xxyyzz_xzzz_0, \
                                         tg_xxyyzz_xzzz_1, tg_xxyyzz_yyy_1, tg_xxyyzz_yyyy_0, tg_xxyyzz_yyyy_1, \
                                         tg_xxyyzz_yyyz_0, tg_xxyyzz_yyyz_1, tg_xxyyzz_yyz_1, tg_xxyyzz_yyzz_0, \
                                         tg_xxyyzz_yyzz_1, tg_xxyyzz_yzz_1, tg_xxyyzz_yzzz_0, tg_xxyyzz_yzzz_1, \
                                         tg_xxyyzz_zzz_1, tg_xxyyzz_zzzz_0, tg_xxyyzz_zzzz_1, tg_xxyzzz_xxx_1, \
                                         tg_xxyzzz_xxxx_0, tg_xxyzzz_xxxx_1, tg_xxyzzz_xxxy_0, tg_xxyzzz_xxxy_1, \
                                         tg_xxyzzz_xxxz_0, tg_xxyzzz_xxxz_1, tg_xxyzzz_xxy_1, tg_xxyzzz_xxyy_0, \
                                         tg_xxyzzz_xxyy_1, tg_xxyzzz_xxyz_0, tg_xxyzzz_xxyz_1, tg_xxyzzz_xxz_1, \
                                         tg_xxyzzz_xxzz_0, tg_xxyzzz_xxzz_1, tg_xxyzzz_xyy_1, tg_xxyzzz_xyyy_0, \
                                         tg_xxyzzz_xyyy_1, tg_xxyzzz_xyyz_0, tg_xxyzzz_xyyz_1, tg_xxyzzz_xyz_1, \
                                         tg_xxyzzz_xyzz_0, tg_xxyzzz_xyzz_1, tg_xxyzzz_xzz_1, tg_xxyzzz_xzzz_0, \
                                         tg_xxyzzz_xzzz_1, tg_xxyzzz_yyy_1, tg_xxyzzz_yyyy_0, tg_xxyzzz_yyyy_1, \
                                         tg_xxyzzz_yyyz_0, tg_xxyzzz_yyyz_1, tg_xxyzzz_yyz_1, tg_xxyzzz_yyzz_0, \
                                         tg_xxyzzz_yyzz_1, tg_xxyzzz_yzz_1, tg_xxyzzz_yzzz_0, tg_xxyzzz_yzzz_1, \
                                         tg_xxyzzz_zzz_1, tg_xxyzzz_zzzz_0, tg_xxyzzz_zzzz_1, tg_xxzzzz_xxx_1, \
                                         tg_xxzzzz_xxxx_0, tg_xxzzzz_xxxx_1, tg_xxzzzz_xxxy_0, tg_xxzzzz_xxxy_1, \
                                         tg_xxzzzz_xxxz_0, tg_xxzzzz_xxxz_1, tg_xxzzzz_xxy_1, tg_xxzzzz_xxyy_0, \
                                         tg_xxzzzz_xxyy_1, tg_xxzzzz_xxyz_0, tg_xxzzzz_xxyz_1, tg_xxzzzz_xxz_1, \
                                         tg_xxzzzz_xxzz_0, tg_xxzzzz_xxzz_1, tg_xxzzzz_xyy_1, tg_xxzzzz_xyyy_0, \
                                         tg_xxzzzz_xyyy_1, tg_xxzzzz_xyyz_0, tg_xxzzzz_xyyz_1, tg_xxzzzz_xyz_1, \
                                         tg_xxzzzz_xyzz_0, tg_xxzzzz_xyzz_1, tg_xxzzzz_xzz_1, tg_xxzzzz_xzzz_0, \
                                         tg_xxzzzz_xzzz_1, tg_xxzzzz_yyy_1, tg_xxzzzz_yyyy_0, tg_xxzzzz_yyyy_1, \
                                         tg_xxzzzz_yyyz_0, tg_xxzzzz_yyyz_1, tg_xxzzzz_yyz_1, tg_xxzzzz_yyzz_0, \
                                         tg_xxzzzz_yyzz_1, tg_xxzzzz_yzz_1, tg_xxzzzz_yzzz_0, tg_xxzzzz_yzzz_1, \
                                         tg_xxzzzz_zzz_1, tg_xxzzzz_zzzz_0, tg_xxzzzz_zzzz_1, tg_xyyyyy_xxx_1, \
                                         tg_xyyyyy_xxxx_0, tg_xyyyyy_xxxx_1, tg_xyyyyy_xxxy_0, tg_xyyyyy_xxxy_1, \
                                         tg_xyyyyy_xxxz_0, tg_xyyyyy_xxxz_1, tg_xyyyyy_xxy_1, tg_xyyyyy_xxyy_0, \
                                         tg_xyyyyy_xxyy_1, tg_xyyyyy_xxyz_0, tg_xyyyyy_xxyz_1, tg_xyyyyy_xxz_1, \
                                         tg_xyyyyy_xxzz_0, tg_xyyyyy_xxzz_1, tg_xyyyyy_xyy_1, tg_xyyyyy_xyyy_0, \
                                         tg_xyyyyy_xyyy_1, tg_xyyyyy_xyyz_0, tg_xyyyyy_xyyz_1, tg_xyyyyy_xyz_1, \
                                         tg_xyyyyy_xyzz_0, tg_xyyyyy_xyzz_1, tg_xyyyyy_xzz_1, tg_xyyyyy_xzzz_0, \
                                         tg_xyyyyy_xzzz_1, tg_xyyyyy_yyy_1, tg_xyyyyy_yyyy_0, tg_xyyyyy_yyyy_1, \
                                         tg_xyyyyy_yyyz_0, tg_xyyyyy_yyyz_1, tg_xyyyyy_yyz_1, tg_xyyyyy_yyzz_0, \
                                         tg_xyyyyy_yyzz_1, tg_xyyyyy_yzz_1, tg_xyyyyy_yzzz_0, tg_xyyyyy_yzzz_1, \
                                         tg_xyyyyy_zzz_1, tg_xyyyyy_zzzz_0, tg_xyyyyy_zzzz_1, tg_xyyyyz_xxx_1, \
                                         tg_xyyyyz_xxxx_0, tg_xyyyyz_xxxx_1, tg_xyyyyz_xxxy_0, tg_xyyyyz_xxxy_1, \
                                         tg_xyyyyz_xxxz_0, tg_xyyyyz_xxxz_1, tg_xyyyyz_xxy_1, tg_xyyyyz_xxyy_0, \
                                         tg_xyyyyz_xxyy_1, tg_xyyyyz_xxyz_0, tg_xyyyyz_xxyz_1, tg_xyyyyz_xxz_1, \
                                         tg_xyyyyz_xxzz_0, tg_xyyyyz_xxzz_1, tg_xyyyyz_xyy_1, tg_xyyyyz_xyyy_0, \
                                         tg_xyyyyz_xyyy_1, tg_xyyyyz_xyyz_0, tg_xyyyyz_xyyz_1, tg_xyyyyz_xyz_1, \
                                         tg_xyyyyz_xyzz_0, tg_xyyyyz_xyzz_1, tg_xyyyyz_xzz_1, tg_xyyyyz_xzzz_0, \
                                         tg_xyyyyz_xzzz_1, tg_xyyyyz_yyy_1, tg_xyyyyz_yyyy_0, tg_xyyyyz_yyyy_1, \
                                         tg_xyyyyz_yyyz_0, tg_xyyyyz_yyyz_1, tg_xyyyyz_yyz_1, tg_xyyyyz_yyzz_0, \
                                         tg_xyyyyz_yyzz_1, tg_xyyyyz_yzz_1, tg_xyyyyz_yzzz_0, tg_xyyyyz_yzzz_1, \
                                         tg_xyyyyz_zzz_1, tg_xyyyyz_zzzz_0, tg_xyyyyz_zzzz_1, tg_xyyyzz_xxx_1, \
                                         tg_xyyyzz_xxxx_0, tg_xyyyzz_xxxx_1, tg_xyyyzz_xxxy_0, tg_xyyyzz_xxxy_1, \
                                         tg_xyyyzz_xxxz_0, tg_xyyyzz_xxxz_1, tg_xyyyzz_xxy_1, tg_xyyyzz_xxyy_0, \
                                         tg_xyyyzz_xxyy_1, tg_xyyyzz_xxyz_0, tg_xyyyzz_xxyz_1, tg_xyyyzz_xxz_1, \
                                         tg_xyyyzz_xxzz_0, tg_xyyyzz_xxzz_1, tg_xyyyzz_xyy_1, tg_xyyyzz_xyyy_0, \
                                         tg_xyyyzz_xyyy_1, tg_xyyyzz_xyyz_0, tg_xyyyzz_xyyz_1, tg_xyyyzz_xyz_1, \
                                         tg_xyyyzz_xyzz_0, tg_xyyyzz_xyzz_1, tg_xyyyzz_xzz_1, tg_xyyyzz_xzzz_0, \
                                         tg_xyyyzz_xzzz_1, tg_xyyyzz_yyy_1, tg_xyyyzz_yyyy_0, tg_xyyyzz_yyyy_1, \
                                         tg_xyyyzz_yyyz_0, tg_xyyyzz_yyyz_1, tg_xyyyzz_yyz_1, tg_xyyyzz_yyzz_0, \
                                         tg_xyyyzz_yyzz_1, tg_xyyyzz_yzz_1, tg_xyyyzz_yzzz_0, tg_xyyyzz_yzzz_1, \
                                         tg_xyyyzz_zzz_1, tg_xyyyzz_zzzz_0, tg_xyyyzz_zzzz_1, tg_xyyzz_xxxx_0, \
                                         tg_xyyzz_xxxx_1, tg_xyyzz_xxxy_0, tg_xyyzz_xxxy_1, tg_xyyzz_xxxz_0, tg_xyyzz_xxxz_1, \
                                         tg_xyyzz_xxyy_0, tg_xyyzz_xxyy_1, tg_xyyzz_xxyz_0, tg_xyyzz_xxyz_1, tg_xyyzz_xxzz_0, \
                                         tg_xyyzz_xxzz_1, tg_xyyzz_xyyy_0, tg_xyyzz_xyyy_1, tg_xyyzz_xyyz_0, tg_xyyzz_xyyz_1, \
                                         tg_xyyzz_xyzz_0, tg_xyyzz_xyzz_1, tg_xyyzz_xzzz_0, tg_xyyzz_xzzz_1, tg_xyyzz_yyyy_0, \
                                         tg_xyyzz_yyyy_1, tg_xyyzz_yyyz_0, tg_xyyzz_yyyz_1, tg_xyyzz_yyzz_0, tg_xyyzz_yyzz_1, \
                                         tg_xyyzz_yzzz_0, tg_xyyzz_yzzz_1, tg_xyyzz_zzzz_0, tg_xyyzz_zzzz_1, tg_xyzzz_xxxx_0, \
                                         tg_xyzzz_xxxx_1, tg_xyzzz_xxxy_0, tg_xyzzz_xxxy_1, tg_xyzzz_xxxz_0, tg_xyzzz_xxxz_1, \
                                         tg_xyzzz_xxyy_0, tg_xyzzz_xxyy_1, tg_xyzzz_xxyz_0, tg_xyzzz_xxyz_1, tg_xyzzz_xxzz_0, \
                                         tg_xyzzz_xxzz_1, tg_xyzzz_xyyy_0, tg_xyzzz_xyyy_1, tg_xyzzz_xyyz_0, tg_xyzzz_xyyz_1, \
                                         tg_xyzzz_xyzz_0, tg_xyzzz_xyzz_1, tg_xyzzz_xzzz_0, tg_xyzzz_xzzz_1, tg_xyzzz_yyyy_0, \
                                         tg_xyzzz_yyyy_1, tg_xyzzz_yyyz_0, tg_xyzzz_yyyz_1, tg_xyzzz_yyzz_0, tg_xyzzz_yyzz_1, \
                                         tg_xyzzz_yzzz_0, tg_xyzzz_yzzz_1, tg_xyzzz_zzzz_0, tg_xyzzz_zzzz_1, tg_xzzzz_xxxx_0, \
                                         tg_xzzzz_xxxx_1, tg_xzzzz_xxxy_0, tg_xzzzz_xxxy_1, tg_xzzzz_xxxz_0, tg_xzzzz_xxxz_1, \
                                         tg_xzzzz_xxyy_0, tg_xzzzz_xxyy_1, tg_xzzzz_xxyz_0, tg_xzzzz_xxyz_1, tg_xzzzz_xxzz_0, \
                                         tg_xzzzz_xxzz_1, tg_xzzzz_xyyy_0, tg_xzzzz_xyyy_1, tg_xzzzz_xyyz_0, tg_xzzzz_xyyz_1, \
                                         tg_xzzzz_xyzz_0, tg_xzzzz_xyzz_1, tg_xzzzz_xzzz_0, tg_xzzzz_xzzz_1, tg_xzzzz_yyyy_0, \
                                         tg_xzzzz_yyyy_1, tg_xzzzz_yyyz_0, tg_xzzzz_yyyz_1, tg_xzzzz_yyzz_0, tg_xzzzz_yyzz_1, \
                                         tg_xzzzz_yzzz_0, tg_xzzzz_yzzz_1, tg_xzzzz_zzzz_0, tg_xzzzz_zzzz_1, tg_yyyyy_xxxx_0, \
                                         tg_yyyyy_xxxx_1, tg_yyyyy_xxxy_0, tg_yyyyy_xxxy_1, tg_yyyyy_xxxz_0, tg_yyyyy_xxxz_1, \
                                         tg_yyyyy_xxyy_0, tg_yyyyy_xxyy_1, tg_yyyyy_xxyz_0, tg_yyyyy_xxyz_1, tg_yyyyy_xxzz_0, \
                                         tg_yyyyy_xxzz_1, tg_yyyyy_xyyy_0, tg_yyyyy_xyyy_1, tg_yyyyy_xyyz_0, tg_yyyyy_xyyz_1, \
                                         tg_yyyyy_xyzz_0, tg_yyyyy_xyzz_1, tg_yyyyy_xzzz_0, tg_yyyyy_xzzz_1, tg_yyyyy_yyyy_0, \
                                         tg_yyyyy_yyyy_1, tg_yyyyy_yyyz_0, tg_yyyyy_yyyz_1, tg_yyyyy_yyzz_0, tg_yyyyy_yyzz_1, \
                                         tg_yyyyy_yzzz_0, tg_yyyyy_yzzz_1, tg_yyyyy_zzzz_0, tg_yyyyy_zzzz_1, tg_yyyyz_xxxx_0, \
                                         tg_yyyyz_xxxx_1, tg_yyyyz_xxxy_0, tg_yyyyz_xxxy_1, tg_yyyyz_xxxz_0, tg_yyyyz_xxxz_1, \
                                         tg_yyyyz_xxyy_0, tg_yyyyz_xxyy_1, tg_yyyyz_xxyz_0, tg_yyyyz_xxyz_1, tg_yyyyz_xxzz_0, \
                                         tg_yyyyz_xxzz_1, tg_yyyyz_xyyy_0, tg_yyyyz_xyyy_1, tg_yyyyz_xyyz_0, tg_yyyyz_xyyz_1, \
                                         tg_yyyyz_xyzz_0, tg_yyyyz_xyzz_1, tg_yyyyz_xzzz_0, tg_yyyyz_xzzz_1, tg_yyyyz_yyyy_0, \
                                         tg_yyyyz_yyyy_1, tg_yyyyz_yyyz_0, tg_yyyyz_yyyz_1, tg_yyyyz_yyzz_0, tg_yyyyz_yyzz_1, \
                                         tg_yyyyz_yzzz_0, tg_yyyyz_yzzz_1, tg_yyyyz_zzzz_0, tg_yyyyz_zzzz_1, tg_yyyzz_xxxx_0, \
                                         tg_yyyzz_xxxx_1, tg_yyyzz_xxxy_0, tg_yyyzz_xxxy_1, tg_yyyzz_xxxz_0, tg_yyyzz_xxxz_1, \
                                         tg_yyyzz_xxyy_0, tg_yyyzz_xxyy_1, tg_yyyzz_xxyz_0, tg_yyyzz_xxyz_1, tg_yyyzz_xxzz_0, \
                                         tg_yyyzz_xxzz_1, tg_yyyzz_xyyy_0, tg_yyyzz_xyyy_1, tg_yyyzz_xyyz_0, tg_yyyzz_xyyz_1, \
                                         tg_yyyzz_xyzz_0, tg_yyyzz_xyzz_1, tg_yyyzz_xzzz_0, tg_yyyzz_xzzz_1, tg_yyyzz_yyyy_0, \
                                         tg_yyyzz_yyyy_1, tg_yyyzz_yyyz_0, tg_yyyzz_yyyz_1, tg_yyyzz_yyzz_0, tg_yyyzz_yyzz_1, \
                                         tg_yyyzz_yzzz_0, tg_yyyzz_yzzz_1, tg_yyyzz_zzzz_0, tg_yyyzz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyyzz_xxxx_0[j] = pb_x * tg_xxyyzz_xxxx_0[j] + fr * tg_xxyyzz_xxxx_1[j] + fl1_fx * (tg_xyyzz_xxxx_0[j] - tg_xyyzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyyzz_xxx_1[j];

                    tg_xxxyyzz_xxxy_0[j] = pb_x * tg_xxyyzz_xxxy_0[j] + fr * tg_xxyyzz_xxxy_1[j] + fl1_fx * (tg_xyyzz_xxxy_0[j] - tg_xyyzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxy_1[j];

                    tg_xxxyyzz_xxxz_0[j] = pb_x * tg_xxyyzz_xxxz_0[j] + fr * tg_xxyyzz_xxxz_1[j] + fl1_fx * (tg_xyyzz_xxxz_0[j] - tg_xyyzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xxz_1[j];

                    tg_xxxyyzz_xxyy_0[j] = pb_x * tg_xxyyzz_xxyy_0[j] + fr * tg_xxyyzz_xxyy_1[j] + fl1_fx * (tg_xyyzz_xxyy_0[j] - tg_xyyzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyy_1[j];

                    tg_xxxyyzz_xxyz_0[j] = pb_x * tg_xxyyzz_xxyz_0[j] + fr * tg_xxyyzz_xxyz_1[j] + fl1_fx * (tg_xyyzz_xxyz_0[j] - tg_xyyzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xyz_1[j];

                    tg_xxxyyzz_xxzz_0[j] = pb_x * tg_xxyyzz_xxzz_0[j] + fr * tg_xxyyzz_xxzz_1[j] + fl1_fx * (tg_xyyzz_xxzz_0[j] - tg_xyyzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xzz_1[j];

                    tg_xxxyyzz_xyyy_0[j] = pb_x * tg_xxyyzz_xyyy_0[j] + fr * tg_xxyyzz_xyyy_1[j] + fl1_fx * (tg_xyyzz_xyyy_0[j] - tg_xyyzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyy_1[j];

                    tg_xxxyyzz_xyyz_0[j] = pb_x * tg_xxyyzz_xyyz_0[j] + fr * tg_xxyyzz_xyyz_1[j] + fl1_fx * (tg_xyyzz_xyyz_0[j] - tg_xyyzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yyz_1[j];

                    tg_xxxyyzz_xyzz_0[j] = pb_x * tg_xxyyzz_xyzz_0[j] + fr * tg_xxyyzz_xyzz_1[j] + fl1_fx * (tg_xyyzz_xyzz_0[j] - tg_xyyzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yzz_1[j];

                    tg_xxxyyzz_xzzz_0[j] = pb_x * tg_xxyyzz_xzzz_0[j] + fr * tg_xxyyzz_xzzz_1[j] + fl1_fx * (tg_xyyzz_xzzz_0[j] - tg_xyyzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_zzz_1[j];

                    tg_xxxyyzz_yyyy_0[j] = pb_x * tg_xxyyzz_yyyy_0[j] + fr * tg_xxyyzz_yyyy_1[j] + fl1_fx * (tg_xyyzz_yyyy_0[j] - tg_xyyzz_yyyy_1[j] * fl1_fza);

                    tg_xxxyyzz_yyyz_0[j] = pb_x * tg_xxyyzz_yyyz_0[j] + fr * tg_xxyyzz_yyyz_1[j] + fl1_fx * (tg_xyyzz_yyyz_0[j] - tg_xyyzz_yyyz_1[j] * fl1_fza);

                    tg_xxxyyzz_yyzz_0[j] = pb_x * tg_xxyyzz_yyzz_0[j] + fr * tg_xxyyzz_yyzz_1[j] + fl1_fx * (tg_xyyzz_yyzz_0[j] - tg_xyyzz_yyzz_1[j] * fl1_fza);

                    tg_xxxyyzz_yzzz_0[j] = pb_x * tg_xxyyzz_yzzz_0[j] + fr * tg_xxyyzz_yzzz_1[j] + fl1_fx * (tg_xyyzz_yzzz_0[j] - tg_xyyzz_yzzz_1[j] * fl1_fza);

                    tg_xxxyyzz_zzzz_0[j] = pb_x * tg_xxyyzz_zzzz_0[j] + fr * tg_xxyyzz_zzzz_1[j] + fl1_fx * (tg_xyyzz_zzzz_0[j] - tg_xyyzz_zzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_xxxx_0[j] = pb_x * tg_xxyzzz_xxxx_0[j] + fr * tg_xxyzzz_xxxx_1[j] + fl1_fx * (tg_xyzzz_xxxx_0[j] - tg_xyzzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxyzzz_xxx_1[j];

                    tg_xxxyzzz_xxxy_0[j] = pb_x * tg_xxyzzz_xxxy_0[j] + fr * tg_xxyzzz_xxxy_1[j] + fl1_fx * (tg_xyzzz_xxxy_0[j] - tg_xyzzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxy_1[j];

                    tg_xxxyzzz_xxxz_0[j] = pb_x * tg_xxyzzz_xxxz_0[j] + fr * tg_xxyzzz_xxxz_1[j] + fl1_fx * (tg_xyzzz_xxxz_0[j] - tg_xyzzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xxz_1[j];

                    tg_xxxyzzz_xxyy_0[j] = pb_x * tg_xxyzzz_xxyy_0[j] + fr * tg_xxyzzz_xxyy_1[j] + fl1_fx * (tg_xyzzz_xxyy_0[j] - tg_xyzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyy_1[j];

                    tg_xxxyzzz_xxyz_0[j] = pb_x * tg_xxyzzz_xxyz_0[j] + fr * tg_xxyzzz_xxyz_1[j] + fl1_fx * (tg_xyzzz_xxyz_0[j] - tg_xyzzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xyz_1[j];

                    tg_xxxyzzz_xxzz_0[j] = pb_x * tg_xxyzzz_xxzz_0[j] + fr * tg_xxyzzz_xxzz_1[j] + fl1_fx * (tg_xyzzz_xxzz_0[j] - tg_xyzzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xzz_1[j];

                    tg_xxxyzzz_xyyy_0[j] = pb_x * tg_xxyzzz_xyyy_0[j] + fr * tg_xxyzzz_xyyy_1[j] + fl1_fx * (tg_xyzzz_xyyy_0[j] - tg_xyzzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyy_1[j];

                    tg_xxxyzzz_xyyz_0[j] = pb_x * tg_xxyzzz_xyyz_0[j] + fr * tg_xxyzzz_xyyz_1[j] + fl1_fx * (tg_xyzzz_xyyz_0[j] - tg_xyzzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yyz_1[j];

                    tg_xxxyzzz_xyzz_0[j] = pb_x * tg_xxyzzz_xyzz_0[j] + fr * tg_xxyzzz_xyzz_1[j] + fl1_fx * (tg_xyzzz_xyzz_0[j] - tg_xyzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yzz_1[j];

                    tg_xxxyzzz_xzzz_0[j] = pb_x * tg_xxyzzz_xzzz_0[j] + fr * tg_xxyzzz_xzzz_1[j] + fl1_fx * (tg_xyzzz_xzzz_0[j] - tg_xyzzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_zzz_1[j];

                    tg_xxxyzzz_yyyy_0[j] = pb_x * tg_xxyzzz_yyyy_0[j] + fr * tg_xxyzzz_yyyy_1[j] + fl1_fx * (tg_xyzzz_yyyy_0[j] - tg_xyzzz_yyyy_1[j] * fl1_fza);

                    tg_xxxyzzz_yyyz_0[j] = pb_x * tg_xxyzzz_yyyz_0[j] + fr * tg_xxyzzz_yyyz_1[j] + fl1_fx * (tg_xyzzz_yyyz_0[j] - tg_xyzzz_yyyz_1[j] * fl1_fza);

                    tg_xxxyzzz_yyzz_0[j] = pb_x * tg_xxyzzz_yyzz_0[j] + fr * tg_xxyzzz_yyzz_1[j] + fl1_fx * (tg_xyzzz_yyzz_0[j] - tg_xyzzz_yyzz_1[j] * fl1_fza);

                    tg_xxxyzzz_yzzz_0[j] = pb_x * tg_xxyzzz_yzzz_0[j] + fr * tg_xxyzzz_yzzz_1[j] + fl1_fx * (tg_xyzzz_yzzz_0[j] - tg_xyzzz_yzzz_1[j] * fl1_fza);

                    tg_xxxyzzz_zzzz_0[j] = pb_x * tg_xxyzzz_zzzz_0[j] + fr * tg_xxyzzz_zzzz_1[j] + fl1_fx * (tg_xyzzz_zzzz_0[j] - tg_xyzzz_zzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_xxxx_0[j] = pb_x * tg_xxzzzz_xxxx_0[j] + fr * tg_xxzzzz_xxxx_1[j] + fl1_fx * (tg_xzzzz_xxxx_0[j] - tg_xzzzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xxzzzz_xxx_1[j];

                    tg_xxxzzzz_xxxy_0[j] = pb_x * tg_xxzzzz_xxxy_0[j] + fr * tg_xxzzzz_xxxy_1[j] + fl1_fx * (tg_xzzzz_xxxy_0[j] - tg_xzzzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxy_1[j];

                    tg_xxxzzzz_xxxz_0[j] = pb_x * tg_xxzzzz_xxxz_0[j] + fr * tg_xxzzzz_xxxz_1[j] + fl1_fx * (tg_xzzzz_xxxz_0[j] - tg_xzzzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xxz_1[j];

                    tg_xxxzzzz_xxyy_0[j] = pb_x * tg_xxzzzz_xxyy_0[j] + fr * tg_xxzzzz_xxyy_1[j] + fl1_fx * (tg_xzzzz_xxyy_0[j] - tg_xzzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyy_1[j];

                    tg_xxxzzzz_xxyz_0[j] = pb_x * tg_xxzzzz_xxyz_0[j] + fr * tg_xxzzzz_xxyz_1[j] + fl1_fx * (tg_xzzzz_xxyz_0[j] - tg_xzzzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xyz_1[j];

                    tg_xxxzzzz_xxzz_0[j] = pb_x * tg_xxzzzz_xxzz_0[j] + fr * tg_xxzzzz_xxzz_1[j] + fl1_fx * (tg_xzzzz_xxzz_0[j] - tg_xzzzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xzz_1[j];

                    tg_xxxzzzz_xyyy_0[j] = pb_x * tg_xxzzzz_xyyy_0[j] + fr * tg_xxzzzz_xyyy_1[j] + fl1_fx * (tg_xzzzz_xyyy_0[j] - tg_xzzzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyy_1[j];

                    tg_xxxzzzz_xyyz_0[j] = pb_x * tg_xxzzzz_xyyz_0[j] + fr * tg_xxzzzz_xyyz_1[j] + fl1_fx * (tg_xzzzz_xyyz_0[j] - tg_xzzzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yyz_1[j];

                    tg_xxxzzzz_xyzz_0[j] = pb_x * tg_xxzzzz_xyzz_0[j] + fr * tg_xxzzzz_xyzz_1[j] + fl1_fx * (tg_xzzzz_xyzz_0[j] - tg_xzzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yzz_1[j];

                    tg_xxxzzzz_xzzz_0[j] = pb_x * tg_xxzzzz_xzzz_0[j] + fr * tg_xxzzzz_xzzz_1[j] + fl1_fx * (tg_xzzzz_xzzz_0[j] - tg_xzzzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_zzz_1[j];

                    tg_xxxzzzz_yyyy_0[j] = pb_x * tg_xxzzzz_yyyy_0[j] + fr * tg_xxzzzz_yyyy_1[j] + fl1_fx * (tg_xzzzz_yyyy_0[j] - tg_xzzzz_yyyy_1[j] * fl1_fza);

                    tg_xxxzzzz_yyyz_0[j] = pb_x * tg_xxzzzz_yyyz_0[j] + fr * tg_xxzzzz_yyyz_1[j] + fl1_fx * (tg_xzzzz_yyyz_0[j] - tg_xzzzz_yyyz_1[j] * fl1_fza);

                    tg_xxxzzzz_yyzz_0[j] = pb_x * tg_xxzzzz_yyzz_0[j] + fr * tg_xxzzzz_yyzz_1[j] + fl1_fx * (tg_xzzzz_yyzz_0[j] - tg_xzzzz_yyzz_1[j] * fl1_fza);

                    tg_xxxzzzz_yzzz_0[j] = pb_x * tg_xxzzzz_yzzz_0[j] + fr * tg_xxzzzz_yzzz_1[j] + fl1_fx * (tg_xzzzz_yzzz_0[j] - tg_xzzzz_yzzz_1[j] * fl1_fza);

                    tg_xxxzzzz_zzzz_0[j] = pb_x * tg_xxzzzz_zzzz_0[j] + fr * tg_xxzzzz_zzzz_1[j] + fl1_fx * (tg_xzzzz_zzzz_0[j] - tg_xzzzz_zzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_xxxx_0[j] = pb_x * tg_xyyyyy_xxxx_0[j] + fr * tg_xyyyyy_xxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxx_0[j] - tg_yyyyy_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyy_xxx_1[j];

                    tg_xxyyyyy_xxxy_0[j] = pb_x * tg_xyyyyy_xxxy_0[j] + fr * tg_xyyyyy_xxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxy_0[j] - tg_yyyyy_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxy_1[j];

                    tg_xxyyyyy_xxxz_0[j] = pb_x * tg_xyyyyy_xxxz_0[j] + fr * tg_xyyyyy_xxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxxz_0[j] - tg_yyyyy_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xxz_1[j];

                    tg_xxyyyyy_xxyy_0[j] = pb_x * tg_xyyyyy_xxyy_0[j] + fr * tg_xyyyyy_xxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyy_0[j] - tg_yyyyy_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyy_1[j];

                    tg_xxyyyyy_xxyz_0[j] = pb_x * tg_xyyyyy_xxyz_0[j] + fr * tg_xyyyyy_xxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxyz_0[j] - tg_yyyyy_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xyz_1[j];

                    tg_xxyyyyy_xxzz_0[j] = pb_x * tg_xyyyyy_xxzz_0[j] + fr * tg_xyyyyy_xxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxzz_0[j] - tg_yyyyy_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xzz_1[j];

                    tg_xxyyyyy_xyyy_0[j] = pb_x * tg_xyyyyy_xyyy_0[j] + fr * tg_xyyyyy_xyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyy_0[j] - tg_yyyyy_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyy_1[j];

                    tg_xxyyyyy_xyyz_0[j] = pb_x * tg_xyyyyy_xyyz_0[j] + fr * tg_xyyyyy_xyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyyz_0[j] - tg_yyyyy_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yyz_1[j];

                    tg_xxyyyyy_xyzz_0[j] = pb_x * tg_xyyyyy_xyzz_0[j] + fr * tg_xyyyyy_xyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyzz_0[j] - tg_yyyyy_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yzz_1[j];

                    tg_xxyyyyy_xzzz_0[j] = pb_x * tg_xyyyyy_xzzz_0[j] + fr * tg_xyyyyy_xzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xzzz_0[j] - tg_yyyyy_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_zzz_1[j];

                    tg_xxyyyyy_yyyy_0[j] = pb_x * tg_xyyyyy_yyyy_0[j] + fr * tg_xyyyyy_yyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyy_0[j] - tg_yyyyy_yyyy_1[j] * fl1_fza);

                    tg_xxyyyyy_yyyz_0[j] = pb_x * tg_xyyyyy_yyyz_0[j] + fr * tg_xyyyyy_yyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyyz_0[j] - tg_yyyyy_yyyz_1[j] * fl1_fza);

                    tg_xxyyyyy_yyzz_0[j] = pb_x * tg_xyyyyy_yyzz_0[j] + fr * tg_xyyyyy_yyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyzz_0[j] - tg_yyyyy_yyzz_1[j] * fl1_fza);

                    tg_xxyyyyy_yzzz_0[j] = pb_x * tg_xyyyyy_yzzz_0[j] + fr * tg_xyyyyy_yzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yzzz_0[j] - tg_yyyyy_yzzz_1[j] * fl1_fza);

                    tg_xxyyyyy_zzzz_0[j] = pb_x * tg_xyyyyy_zzzz_0[j] + fr * tg_xyyyyy_zzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_zzzz_0[j] - tg_yyyyy_zzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_xxxx_0[j] = pb_x * tg_xyyyyz_xxxx_0[j] + fr * tg_xyyyyz_xxxx_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxx_0[j] - tg_yyyyz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyyz_xxx_1[j];

                    tg_xxyyyyz_xxxy_0[j] = pb_x * tg_xyyyyz_xxxy_0[j] + fr * tg_xyyyyz_xxxy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxy_0[j] - tg_yyyyz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxy_1[j];

                    tg_xxyyyyz_xxxz_0[j] = pb_x * tg_xyyyyz_xxxz_0[j] + fr * tg_xyyyyz_xxxz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxxz_0[j] - tg_yyyyz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xxz_1[j];

                    tg_xxyyyyz_xxyy_0[j] = pb_x * tg_xyyyyz_xxyy_0[j] + fr * tg_xyyyyz_xxyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyy_0[j] - tg_yyyyz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyy_1[j];

                    tg_xxyyyyz_xxyz_0[j] = pb_x * tg_xyyyyz_xxyz_0[j] + fr * tg_xyyyyz_xxyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxyz_0[j] - tg_yyyyz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xyz_1[j];

                    tg_xxyyyyz_xxzz_0[j] = pb_x * tg_xyyyyz_xxzz_0[j] + fr * tg_xyyyyz_xxzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxzz_0[j] - tg_yyyyz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xzz_1[j];

                    tg_xxyyyyz_xyyy_0[j] = pb_x * tg_xyyyyz_xyyy_0[j] + fr * tg_xyyyyz_xyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyy_0[j] - tg_yyyyz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyy_1[j];

                    tg_xxyyyyz_xyyz_0[j] = pb_x * tg_xyyyyz_xyyz_0[j] + fr * tg_xyyyyz_xyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyyz_0[j] - tg_yyyyz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yyz_1[j];

                    tg_xxyyyyz_xyzz_0[j] = pb_x * tg_xyyyyz_xyzz_0[j] + fr * tg_xyyyyz_xyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyzz_0[j] - tg_yyyyz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yzz_1[j];

                    tg_xxyyyyz_xzzz_0[j] = pb_x * tg_xyyyyz_xzzz_0[j] + fr * tg_xyyyyz_xzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xzzz_0[j] - tg_yyyyz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_zzz_1[j];

                    tg_xxyyyyz_yyyy_0[j] = pb_x * tg_xyyyyz_yyyy_0[j] + fr * tg_xyyyyz_yyyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyy_0[j] - tg_yyyyz_yyyy_1[j] * fl1_fza);

                    tg_xxyyyyz_yyyz_0[j] = pb_x * tg_xyyyyz_yyyz_0[j] + fr * tg_xyyyyz_yyyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyyz_0[j] - tg_yyyyz_yyyz_1[j] * fl1_fza);

                    tg_xxyyyyz_yyzz_0[j] = pb_x * tg_xyyyyz_yyzz_0[j] + fr * tg_xyyyyz_yyzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyzz_0[j] - tg_yyyyz_yyzz_1[j] * fl1_fza);

                    tg_xxyyyyz_yzzz_0[j] = pb_x * tg_xyyyyz_yzzz_0[j] + fr * tg_xyyyyz_yzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yzzz_0[j] - tg_yyyyz_yzzz_1[j] * fl1_fza);

                    tg_xxyyyyz_zzzz_0[j] = pb_x * tg_xyyyyz_zzzz_0[j] + fr * tg_xyyyyz_zzzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_zzzz_0[j] - tg_yyyyz_zzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_xxxx_0[j] = pb_x * tg_xyyyzz_xxxx_0[j] + fr * tg_xyyyzz_xxxx_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxx_0[j] - tg_yyyzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyyzz_xxx_1[j];

                    tg_xxyyyzz_xxxy_0[j] = pb_x * tg_xyyyzz_xxxy_0[j] + fr * tg_xyyyzz_xxxy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxy_0[j] - tg_yyyzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxy_1[j];

                    tg_xxyyyzz_xxxz_0[j] = pb_x * tg_xyyyzz_xxxz_0[j] + fr * tg_xyyyzz_xxxz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxxz_0[j] - tg_yyyzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xxz_1[j];

                    tg_xxyyyzz_xxyy_0[j] = pb_x * tg_xyyyzz_xxyy_0[j] + fr * tg_xyyyzz_xxyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyy_0[j] - tg_yyyzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyy_1[j];

                    tg_xxyyyzz_xxyz_0[j] = pb_x * tg_xyyyzz_xxyz_0[j] + fr * tg_xyyyzz_xxyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxyz_0[j] - tg_yyyzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xyz_1[j];

                    tg_xxyyyzz_xxzz_0[j] = pb_x * tg_xyyyzz_xxzz_0[j] + fr * tg_xyyyzz_xxzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxzz_0[j] - tg_yyyzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xzz_1[j];

                    tg_xxyyyzz_xyyy_0[j] = pb_x * tg_xyyyzz_xyyy_0[j] + fr * tg_xyyyzz_xyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyy_0[j] - tg_yyyzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyy_1[j];

                    tg_xxyyyzz_xyyz_0[j] = pb_x * tg_xyyyzz_xyyz_0[j] + fr * tg_xyyyzz_xyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyyz_0[j] - tg_yyyzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yyz_1[j];

                    tg_xxyyyzz_xyzz_0[j] = pb_x * tg_xyyyzz_xyzz_0[j] + fr * tg_xyyyzz_xyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyzz_0[j] - tg_yyyzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yzz_1[j];

                    tg_xxyyyzz_xzzz_0[j] = pb_x * tg_xyyyzz_xzzz_0[j] + fr * tg_xyyyzz_xzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xzzz_0[j] - tg_yyyzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_zzz_1[j];

                    tg_xxyyyzz_yyyy_0[j] = pb_x * tg_xyyyzz_yyyy_0[j] + fr * tg_xyyyzz_yyyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyy_0[j] - tg_yyyzz_yyyy_1[j] * fl1_fza);

                    tg_xxyyyzz_yyyz_0[j] = pb_x * tg_xyyyzz_yyyz_0[j] + fr * tg_xyyyzz_yyyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyyz_0[j] - tg_yyyzz_yyyz_1[j] * fl1_fza);

                    tg_xxyyyzz_yyzz_0[j] = pb_x * tg_xyyyzz_yyzz_0[j] + fr * tg_xyyyzz_yyzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyzz_0[j] - tg_yyyzz_yyzz_1[j] * fl1_fza);

                    tg_xxyyyzz_yzzz_0[j] = pb_x * tg_xyyyzz_yzzz_0[j] + fr * tg_xyyyzz_yzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yzzz_0[j] - tg_yyyzz_yzzz_1[j] * fl1_fza);

                    tg_xxyyyzz_zzzz_0[j] = pb_x * tg_xyyyzz_zzzz_0[j] + fr * tg_xyyyzz_zzzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_zzzz_0[j] - tg_yyyzz_zzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSG_270_360(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (270,360)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xyyzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 270); 

                auto tg_xyyzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 271); 

                auto tg_xyyzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 272); 

                auto tg_xyyzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 273); 

                auto tg_xyyzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 274); 

                auto tg_xyyzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 275); 

                auto tg_xyyzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 276); 

                auto tg_xyyzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 277); 

                auto tg_xyyzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 278); 

                auto tg_xyyzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 279); 

                auto tg_xyyzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 280); 

                auto tg_xyyzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 281); 

                auto tg_xyyzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 282); 

                auto tg_xyyzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 283); 

                auto tg_xyyzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 284); 

                auto tg_xyzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 285); 

                auto tg_xyzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 286); 

                auto tg_xyzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 287); 

                auto tg_xyzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 288); 

                auto tg_xyzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 289); 

                auto tg_xyzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 290); 

                auto tg_xyzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 291); 

                auto tg_xyzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 292); 

                auto tg_xyzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 293); 

                auto tg_xyzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 294); 

                auto tg_xyzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 295); 

                auto tg_xyzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 296); 

                auto tg_xyzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 297); 

                auto tg_xyzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 298); 

                auto tg_xyzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 299); 

                auto tg_xzzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 300); 

                auto tg_xzzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 301); 

                auto tg_xzzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 302); 

                auto tg_xzzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 303); 

                auto tg_xzzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 304); 

                auto tg_xzzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 305); 

                auto tg_xzzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 306); 

                auto tg_xzzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 307); 

                auto tg_xzzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 308); 

                auto tg_xzzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 309); 

                auto tg_xzzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 310); 

                auto tg_xzzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 311); 

                auto tg_xzzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 312); 

                auto tg_xzzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 313); 

                auto tg_xzzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 314); 

                auto tg_yyyyyy_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 315); 

                auto tg_yyyyyy_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 316); 

                auto tg_yyyyyy_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 317); 

                auto tg_yyyyyy_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 318); 

                auto tg_yyyyyy_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 319); 

                auto tg_yyyyyy_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 320); 

                auto tg_yyyyyy_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 321); 

                auto tg_yyyyyy_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 322); 

                auto tg_yyyyyy_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 323); 

                auto tg_yyyyyy_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 324); 

                auto tg_yyyyyy_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 325); 

                auto tg_yyyyyy_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 326); 

                auto tg_yyyyyy_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 327); 

                auto tg_yyyyyy_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 328); 

                auto tg_yyyyyy_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 329); 

                auto tg_yyyyyz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 330); 

                auto tg_yyyyyz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 331); 

                auto tg_yyyyyz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 332); 

                auto tg_yyyyyz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 333); 

                auto tg_yyyyyz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 334); 

                auto tg_yyyyyz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 335); 

                auto tg_yyyyyz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 336); 

                auto tg_yyyyyz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 337); 

                auto tg_yyyyyz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 338); 

                auto tg_yyyyyz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 339); 

                auto tg_yyyyyz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 340); 

                auto tg_yyyyyz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 341); 

                auto tg_yyyyyz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 342); 

                auto tg_yyyyyz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 343); 

                auto tg_yyyyyz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 344); 

                auto tg_yyyyzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 345); 

                auto tg_yyyyzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 346); 

                auto tg_yyyyzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 347); 

                auto tg_yyyyzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 348); 

                auto tg_yyyyzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 349); 

                auto tg_yyyyzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 350); 

                auto tg_yyyyzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 351); 

                auto tg_yyyyzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 352); 

                auto tg_yyyyzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 353); 

                auto tg_yyyyzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 354); 

                auto tg_yyyyzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 355); 

                auto tg_yyyyzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 356); 

                auto tg_yyyyzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 357); 

                auto tg_yyyyzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 358); 

                auto tg_yyyyzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 359); 

                auto tg_xyyzzz_xxxx_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 270); 

                auto tg_xyyzzz_xxxy_1 = primBuffer[pidx_g_6_4_m1].data(420 * idx + 271); 

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

                auto tg_yyzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 270); 

                auto tg_yyzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 271); 

                auto tg_yyzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 272); 

                auto tg_yyzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 273); 

                auto tg_yyzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 274); 

                auto tg_yyzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 275); 

                auto tg_yyzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 276); 

                auto tg_yyzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 277); 

                auto tg_yyzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 278); 

                auto tg_yyzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 279); 

                auto tg_yyzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 280); 

                auto tg_yyzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 281); 

                auto tg_yyzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 282); 

                auto tg_yyzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 283); 

                auto tg_yyzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 284); 

                auto tg_yzzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 285); 

                auto tg_yzzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 286); 

                auto tg_yzzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 287); 

                auto tg_yzzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 288); 

                auto tg_yzzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 289); 

                auto tg_yzzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 290); 

                auto tg_yzzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 291); 

                auto tg_yzzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 292); 

                auto tg_yzzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 293); 

                auto tg_yzzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 294); 

                auto tg_yzzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 295); 

                auto tg_yzzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 296); 

                auto tg_yzzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 297); 

                auto tg_yzzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 298); 

                auto tg_yzzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 299); 

                auto tg_zzzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 300); 

                auto tg_zzzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 301); 

                auto tg_zzzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 302); 

                auto tg_zzzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 303); 

                auto tg_zzzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 304); 

                auto tg_zzzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 305); 

                auto tg_zzzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 306); 

                auto tg_zzzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 307); 

                auto tg_zzzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 308); 

                auto tg_zzzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 309); 

                auto tg_zzzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 310); 

                auto tg_zzzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 311); 

                auto tg_zzzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 312); 

                auto tg_zzzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 313); 

                auto tg_zzzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 314); 

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

                auto tg_xyyzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 180); 

                auto tg_xyyzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 181); 

                auto tg_xyyzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 182); 

                auto tg_xyyzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 183); 

                auto tg_xyyzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 184); 

                auto tg_xyyzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 185); 

                auto tg_xyyzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 186); 

                auto tg_xyyzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 187); 

                auto tg_xyyzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 188); 

                auto tg_xyyzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 189); 

                auto tg_xyzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 190); 

                auto tg_xyzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 191); 

                auto tg_xyzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 192); 

                auto tg_xyzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 193); 

                auto tg_xyzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 194); 

                auto tg_xyzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 195); 

                auto tg_xyzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 196); 

                auto tg_xyzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 197); 

                auto tg_xyzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 198); 

                auto tg_xyzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 199); 

                auto tg_xzzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 200); 

                auto tg_xzzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 201); 

                auto tg_xzzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 202); 

                auto tg_xzzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 203); 

                auto tg_xzzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 204); 

                auto tg_xzzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 205); 

                auto tg_xzzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 206); 

                auto tg_xzzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 207); 

                auto tg_xzzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 208); 

                auto tg_xzzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 209); 

                auto tg_yyyyyy_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 210); 

                auto tg_yyyyyy_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 211); 

                auto tg_yyyyyy_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 212); 

                auto tg_yyyyyy_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 213); 

                auto tg_yyyyyy_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 214); 

                auto tg_yyyyyy_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 215); 

                auto tg_yyyyyy_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 216); 

                auto tg_yyyyyy_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 217); 

                auto tg_yyyyyy_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 218); 

                auto tg_yyyyyy_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 219); 

                auto tg_yyyyyz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 220); 

                auto tg_yyyyyz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 221); 

                auto tg_yyyyyz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 222); 

                auto tg_yyyyyz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 223); 

                auto tg_yyyyyz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 224); 

                auto tg_yyyyyz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 225); 

                auto tg_yyyyyz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 226); 

                auto tg_yyyyyz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 227); 

                auto tg_yyyyyz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 228); 

                auto tg_yyyyyz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 229); 

                auto tg_yyyyzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 230); 

                auto tg_yyyyzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 231); 

                auto tg_yyyyzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 232); 

                auto tg_yyyyzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 233); 

                auto tg_yyyyzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 234); 

                auto tg_yyyyzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 235); 

                auto tg_yyyyzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 236); 

                auto tg_yyyyzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 237); 

                auto tg_yyyyzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 238); 

                auto tg_yyyyzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 239); 

                // set up pointers to integrals

                auto tg_xxyyzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 270); 

                auto tg_xxyyzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 271); 

                auto tg_xxyyzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 272); 

                auto tg_xxyyzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 273); 

                auto tg_xxyyzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 274); 

                auto tg_xxyyzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 275); 

                auto tg_xxyyzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 276); 

                auto tg_xxyyzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 277); 

                auto tg_xxyyzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 278); 

                auto tg_xxyyzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 279); 

                auto tg_xxyyzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 280); 

                auto tg_xxyyzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 281); 

                auto tg_xxyyzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 282); 

                auto tg_xxyyzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 283); 

                auto tg_xxyyzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 284); 

                auto tg_xxyzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 285); 

                auto tg_xxyzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 286); 

                auto tg_xxyzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 287); 

                auto tg_xxyzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 288); 

                auto tg_xxyzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 289); 

                auto tg_xxyzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 290); 

                auto tg_xxyzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 291); 

                auto tg_xxyzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 292); 

                auto tg_xxyzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 293); 

                auto tg_xxyzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 294); 

                auto tg_xxyzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 295); 

                auto tg_xxyzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 296); 

                auto tg_xxyzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 297); 

                auto tg_xxyzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 298); 

                auto tg_xxyzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 299); 

                auto tg_xxzzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 300); 

                auto tg_xxzzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 301); 

                auto tg_xxzzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 302); 

                auto tg_xxzzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 303); 

                auto tg_xxzzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 304); 

                auto tg_xxzzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 305); 

                auto tg_xxzzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 306); 

                auto tg_xxzzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 307); 

                auto tg_xxzzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 308); 

                auto tg_xxzzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 309); 

                auto tg_xxzzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 310); 

                auto tg_xxzzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 311); 

                auto tg_xxzzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 312); 

                auto tg_xxzzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 313); 

                auto tg_xxzzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 314); 

                auto tg_xyyyyyy_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 315); 

                auto tg_xyyyyyy_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 316); 

                auto tg_xyyyyyy_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 317); 

                auto tg_xyyyyyy_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 318); 

                auto tg_xyyyyyy_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 319); 

                auto tg_xyyyyyy_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 320); 

                auto tg_xyyyyyy_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 321); 

                auto tg_xyyyyyy_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 322); 

                auto tg_xyyyyyy_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 323); 

                auto tg_xyyyyyy_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 324); 

                auto tg_xyyyyyy_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 325); 

                auto tg_xyyyyyy_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 326); 

                auto tg_xyyyyyy_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 327); 

                auto tg_xyyyyyy_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 328); 

                auto tg_xyyyyyy_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 329); 

                auto tg_xyyyyyz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 330); 

                auto tg_xyyyyyz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 331); 

                auto tg_xyyyyyz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 332); 

                auto tg_xyyyyyz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 333); 

                auto tg_xyyyyyz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 334); 

                auto tg_xyyyyyz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 335); 

                auto tg_xyyyyyz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 336); 

                auto tg_xyyyyyz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 337); 

                auto tg_xyyyyyz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 338); 

                auto tg_xyyyyyz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 339); 

                auto tg_xyyyyyz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 340); 

                auto tg_xyyyyyz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 341); 

                auto tg_xyyyyyz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 342); 

                auto tg_xyyyyyz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 343); 

                auto tg_xyyyyyz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 344); 

                auto tg_xyyyyzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 345); 

                auto tg_xyyyyzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 346); 

                auto tg_xyyyyzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 347); 

                auto tg_xyyyyzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 348); 

                auto tg_xyyyyzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 349); 

                auto tg_xyyyyzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 350); 

                auto tg_xyyyyzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 351); 

                auto tg_xyyyyzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 352); 

                auto tg_xyyyyzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 353); 

                auto tg_xyyyyzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 354); 

                auto tg_xyyyyzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 355); 

                auto tg_xyyyyzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 356); 

                auto tg_xyyyyzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 357); 

                auto tg_xyyyyzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 358); 

                auto tg_xyyyyzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 359); 

                // Batch of Integrals (270,360)

                #pragma omp simd aligned(fxn, fza, tg_xxyyzzz_xxxx_0, tg_xxyyzzz_xxxy_0, tg_xxyyzzz_xxxz_0, \
                                         tg_xxyyzzz_xxyy_0, tg_xxyyzzz_xxyz_0, tg_xxyyzzz_xxzz_0, tg_xxyyzzz_xyyy_0, \
                                         tg_xxyyzzz_xyyz_0, tg_xxyyzzz_xyzz_0, tg_xxyyzzz_xzzz_0, tg_xxyyzzz_yyyy_0, \
                                         tg_xxyyzzz_yyyz_0, tg_xxyyzzz_yyzz_0, tg_xxyyzzz_yzzz_0, tg_xxyyzzz_zzzz_0, \
                                         tg_xxyzzzz_xxxx_0, tg_xxyzzzz_xxxy_0, tg_xxyzzzz_xxxz_0, tg_xxyzzzz_xxyy_0, \
                                         tg_xxyzzzz_xxyz_0, tg_xxyzzzz_xxzz_0, tg_xxyzzzz_xyyy_0, tg_xxyzzzz_xyyz_0, \
                                         tg_xxyzzzz_xyzz_0, tg_xxyzzzz_xzzz_0, tg_xxyzzzz_yyyy_0, tg_xxyzzzz_yyyz_0, \
                                         tg_xxyzzzz_yyzz_0, tg_xxyzzzz_yzzz_0, tg_xxyzzzz_zzzz_0, tg_xxzzzzz_xxxx_0, \
                                         tg_xxzzzzz_xxxy_0, tg_xxzzzzz_xxxz_0, tg_xxzzzzz_xxyy_0, tg_xxzzzzz_xxyz_0, \
                                         tg_xxzzzzz_xxzz_0, tg_xxzzzzz_xyyy_0, tg_xxzzzzz_xyyz_0, tg_xxzzzzz_xyzz_0, \
                                         tg_xxzzzzz_xzzz_0, tg_xxzzzzz_yyyy_0, tg_xxzzzzz_yyyz_0, tg_xxzzzzz_yyzz_0, \
                                         tg_xxzzzzz_yzzz_0, tg_xxzzzzz_zzzz_0, tg_xyyyyyy_xxxx_0, tg_xyyyyyy_xxxy_0, \
                                         tg_xyyyyyy_xxxz_0, tg_xyyyyyy_xxyy_0, tg_xyyyyyy_xxyz_0, tg_xyyyyyy_xxzz_0, \
                                         tg_xyyyyyy_xyyy_0, tg_xyyyyyy_xyyz_0, tg_xyyyyyy_xyzz_0, tg_xyyyyyy_xzzz_0, \
                                         tg_xyyyyyy_yyyy_0, tg_xyyyyyy_yyyz_0, tg_xyyyyyy_yyzz_0, tg_xyyyyyy_yzzz_0, \
                                         tg_xyyyyyy_zzzz_0, tg_xyyyyyz_xxxx_0, tg_xyyyyyz_xxxy_0, tg_xyyyyyz_xxxz_0, \
                                         tg_xyyyyyz_xxyy_0, tg_xyyyyyz_xxyz_0, tg_xyyyyyz_xxzz_0, tg_xyyyyyz_xyyy_0, \
                                         tg_xyyyyyz_xyyz_0, tg_xyyyyyz_xyzz_0, tg_xyyyyyz_xzzz_0, tg_xyyyyyz_yyyy_0, \
                                         tg_xyyyyyz_yyyz_0, tg_xyyyyyz_yyzz_0, tg_xyyyyyz_yzzz_0, tg_xyyyyyz_zzzz_0, \
                                         tg_xyyyyzz_xxxx_0, tg_xyyyyzz_xxxy_0, tg_xyyyyzz_xxxz_0, tg_xyyyyzz_xxyy_0, \
                                         tg_xyyyyzz_xxyz_0, tg_xyyyyzz_xxzz_0, tg_xyyyyzz_xyyy_0, tg_xyyyyzz_xyyz_0, \
                                         tg_xyyyyzz_xyzz_0, tg_xyyyyzz_xzzz_0, tg_xyyyyzz_yyyy_0, tg_xyyyyzz_yyyz_0, \
                                         tg_xyyyyzz_yyzz_0, tg_xyyyyzz_yzzz_0, tg_xyyyyzz_zzzz_0, tg_xyyzzz_xxx_1, \
                                         tg_xyyzzz_xxxx_0, tg_xyyzzz_xxxx_1, tg_xyyzzz_xxxy_0, tg_xyyzzz_xxxy_1, \
                                         tg_xyyzzz_xxxz_0, tg_xyyzzz_xxxz_1, tg_xyyzzz_xxy_1, tg_xyyzzz_xxyy_0, \
                                         tg_xyyzzz_xxyy_1, tg_xyyzzz_xxyz_0, tg_xyyzzz_xxyz_1, tg_xyyzzz_xxz_1, \
                                         tg_xyyzzz_xxzz_0, tg_xyyzzz_xxzz_1, tg_xyyzzz_xyy_1, tg_xyyzzz_xyyy_0, \
                                         tg_xyyzzz_xyyy_1, tg_xyyzzz_xyyz_0, tg_xyyzzz_xyyz_1, tg_xyyzzz_xyz_1, \
                                         tg_xyyzzz_xyzz_0, tg_xyyzzz_xyzz_1, tg_xyyzzz_xzz_1, tg_xyyzzz_xzzz_0, \
                                         tg_xyyzzz_xzzz_1, tg_xyyzzz_yyy_1, tg_xyyzzz_yyyy_0, tg_xyyzzz_yyyy_1, \
                                         tg_xyyzzz_yyyz_0, tg_xyyzzz_yyyz_1, tg_xyyzzz_yyz_1, tg_xyyzzz_yyzz_0, \
                                         tg_xyyzzz_yyzz_1, tg_xyyzzz_yzz_1, tg_xyyzzz_yzzz_0, tg_xyyzzz_yzzz_1, \
                                         tg_xyyzzz_zzz_1, tg_xyyzzz_zzzz_0, tg_xyyzzz_zzzz_1, tg_xyzzzz_xxx_1, \
                                         tg_xyzzzz_xxxx_0, tg_xyzzzz_xxxx_1, tg_xyzzzz_xxxy_0, tg_xyzzzz_xxxy_1, \
                                         tg_xyzzzz_xxxz_0, tg_xyzzzz_xxxz_1, tg_xyzzzz_xxy_1, tg_xyzzzz_xxyy_0, \
                                         tg_xyzzzz_xxyy_1, tg_xyzzzz_xxyz_0, tg_xyzzzz_xxyz_1, tg_xyzzzz_xxz_1, \
                                         tg_xyzzzz_xxzz_0, tg_xyzzzz_xxzz_1, tg_xyzzzz_xyy_1, tg_xyzzzz_xyyy_0, \
                                         tg_xyzzzz_xyyy_1, tg_xyzzzz_xyyz_0, tg_xyzzzz_xyyz_1, tg_xyzzzz_xyz_1, \
                                         tg_xyzzzz_xyzz_0, tg_xyzzzz_xyzz_1, tg_xyzzzz_xzz_1, tg_xyzzzz_xzzz_0, \
                                         tg_xyzzzz_xzzz_1, tg_xyzzzz_yyy_1, tg_xyzzzz_yyyy_0, tg_xyzzzz_yyyy_1, \
                                         tg_xyzzzz_yyyz_0, tg_xyzzzz_yyyz_1, tg_xyzzzz_yyz_1, tg_xyzzzz_yyzz_0, \
                                         tg_xyzzzz_yyzz_1, tg_xyzzzz_yzz_1, tg_xyzzzz_yzzz_0, tg_xyzzzz_yzzz_1, \
                                         tg_xyzzzz_zzz_1, tg_xyzzzz_zzzz_0, tg_xyzzzz_zzzz_1, tg_xzzzzz_xxx_1, \
                                         tg_xzzzzz_xxxx_0, tg_xzzzzz_xxxx_1, tg_xzzzzz_xxxy_0, tg_xzzzzz_xxxy_1, \
                                         tg_xzzzzz_xxxz_0, tg_xzzzzz_xxxz_1, tg_xzzzzz_xxy_1, tg_xzzzzz_xxyy_0, \
                                         tg_xzzzzz_xxyy_1, tg_xzzzzz_xxyz_0, tg_xzzzzz_xxyz_1, tg_xzzzzz_xxz_1, \
                                         tg_xzzzzz_xxzz_0, tg_xzzzzz_xxzz_1, tg_xzzzzz_xyy_1, tg_xzzzzz_xyyy_0, \
                                         tg_xzzzzz_xyyy_1, tg_xzzzzz_xyyz_0, tg_xzzzzz_xyyz_1, tg_xzzzzz_xyz_1, \
                                         tg_xzzzzz_xyzz_0, tg_xzzzzz_xyzz_1, tg_xzzzzz_xzz_1, tg_xzzzzz_xzzz_0, \
                                         tg_xzzzzz_xzzz_1, tg_xzzzzz_yyy_1, tg_xzzzzz_yyyy_0, tg_xzzzzz_yyyy_1, \
                                         tg_xzzzzz_yyyz_0, tg_xzzzzz_yyyz_1, tg_xzzzzz_yyz_1, tg_xzzzzz_yyzz_0, \
                                         tg_xzzzzz_yyzz_1, tg_xzzzzz_yzz_1, tg_xzzzzz_yzzz_0, tg_xzzzzz_yzzz_1, \
                                         tg_xzzzzz_zzz_1, tg_xzzzzz_zzzz_0, tg_xzzzzz_zzzz_1, tg_yyyyyy_xxx_1, \
                                         tg_yyyyyy_xxxx_0, tg_yyyyyy_xxxx_1, tg_yyyyyy_xxxy_0, tg_yyyyyy_xxxy_1, \
                                         tg_yyyyyy_xxxz_0, tg_yyyyyy_xxxz_1, tg_yyyyyy_xxy_1, tg_yyyyyy_xxyy_0, \
                                         tg_yyyyyy_xxyy_1, tg_yyyyyy_xxyz_0, tg_yyyyyy_xxyz_1, tg_yyyyyy_xxz_1, \
                                         tg_yyyyyy_xxzz_0, tg_yyyyyy_xxzz_1, tg_yyyyyy_xyy_1, tg_yyyyyy_xyyy_0, \
                                         tg_yyyyyy_xyyy_1, tg_yyyyyy_xyyz_0, tg_yyyyyy_xyyz_1, tg_yyyyyy_xyz_1, \
                                         tg_yyyyyy_xyzz_0, tg_yyyyyy_xyzz_1, tg_yyyyyy_xzz_1, tg_yyyyyy_xzzz_0, \
                                         tg_yyyyyy_xzzz_1, tg_yyyyyy_yyy_1, tg_yyyyyy_yyyy_0, tg_yyyyyy_yyyy_1, \
                                         tg_yyyyyy_yyyz_0, tg_yyyyyy_yyyz_1, tg_yyyyyy_yyz_1, tg_yyyyyy_yyzz_0, \
                                         tg_yyyyyy_yyzz_1, tg_yyyyyy_yzz_1, tg_yyyyyy_yzzz_0, tg_yyyyyy_yzzz_1, \
                                         tg_yyyyyy_zzz_1, tg_yyyyyy_zzzz_0, tg_yyyyyy_zzzz_1, tg_yyyyyz_xxx_1, \
                                         tg_yyyyyz_xxxx_0, tg_yyyyyz_xxxx_1, tg_yyyyyz_xxxy_0, tg_yyyyyz_xxxy_1, \
                                         tg_yyyyyz_xxxz_0, tg_yyyyyz_xxxz_1, tg_yyyyyz_xxy_1, tg_yyyyyz_xxyy_0, \
                                         tg_yyyyyz_xxyy_1, tg_yyyyyz_xxyz_0, tg_yyyyyz_xxyz_1, tg_yyyyyz_xxz_1, \
                                         tg_yyyyyz_xxzz_0, tg_yyyyyz_xxzz_1, tg_yyyyyz_xyy_1, tg_yyyyyz_xyyy_0, \
                                         tg_yyyyyz_xyyy_1, tg_yyyyyz_xyyz_0, tg_yyyyyz_xyyz_1, tg_yyyyyz_xyz_1, \
                                         tg_yyyyyz_xyzz_0, tg_yyyyyz_xyzz_1, tg_yyyyyz_xzz_1, tg_yyyyyz_xzzz_0, \
                                         tg_yyyyyz_xzzz_1, tg_yyyyyz_yyy_1, tg_yyyyyz_yyyy_0, tg_yyyyyz_yyyy_1, \
                                         tg_yyyyyz_yyyz_0, tg_yyyyyz_yyyz_1, tg_yyyyyz_yyz_1, tg_yyyyyz_yyzz_0, \
                                         tg_yyyyyz_yyzz_1, tg_yyyyyz_yzz_1, tg_yyyyyz_yzzz_0, tg_yyyyyz_yzzz_1, \
                                         tg_yyyyyz_zzz_1, tg_yyyyyz_zzzz_0, tg_yyyyyz_zzzz_1, tg_yyyyzz_xxx_1, \
                                         tg_yyyyzz_xxxx_0, tg_yyyyzz_xxxx_1, tg_yyyyzz_xxxy_0, tg_yyyyzz_xxxy_1, \
                                         tg_yyyyzz_xxxz_0, tg_yyyyzz_xxxz_1, tg_yyyyzz_xxy_1, tg_yyyyzz_xxyy_0, \
                                         tg_yyyyzz_xxyy_1, tg_yyyyzz_xxyz_0, tg_yyyyzz_xxyz_1, tg_yyyyzz_xxz_1, \
                                         tg_yyyyzz_xxzz_0, tg_yyyyzz_xxzz_1, tg_yyyyzz_xyy_1, tg_yyyyzz_xyyy_0, \
                                         tg_yyyyzz_xyyy_1, tg_yyyyzz_xyyz_0, tg_yyyyzz_xyyz_1, tg_yyyyzz_xyz_1, \
                                         tg_yyyyzz_xyzz_0, tg_yyyyzz_xyzz_1, tg_yyyyzz_xzz_1, tg_yyyyzz_xzzz_0, \
                                         tg_yyyyzz_xzzz_1, tg_yyyyzz_yyy_1, tg_yyyyzz_yyyy_0, tg_yyyyzz_yyyy_1, \
                                         tg_yyyyzz_yyyz_0, tg_yyyyzz_yyyz_1, tg_yyyyzz_yyz_1, tg_yyyyzz_yyzz_0, \
                                         tg_yyyyzz_yyzz_1, tg_yyyyzz_yzz_1, tg_yyyyzz_yzzz_0, tg_yyyyzz_yzzz_1, \
                                         tg_yyyyzz_zzz_1, tg_yyyyzz_zzzz_0, tg_yyyyzz_zzzz_1, tg_yyzzz_xxxx_0, \
                                         tg_yyzzz_xxxx_1, tg_yyzzz_xxxy_0, tg_yyzzz_xxxy_1, tg_yyzzz_xxxz_0, tg_yyzzz_xxxz_1, \
                                         tg_yyzzz_xxyy_0, tg_yyzzz_xxyy_1, tg_yyzzz_xxyz_0, tg_yyzzz_xxyz_1, tg_yyzzz_xxzz_0, \
                                         tg_yyzzz_xxzz_1, tg_yyzzz_xyyy_0, tg_yyzzz_xyyy_1, tg_yyzzz_xyyz_0, tg_yyzzz_xyyz_1, \
                                         tg_yyzzz_xyzz_0, tg_yyzzz_xyzz_1, tg_yyzzz_xzzz_0, tg_yyzzz_xzzz_1, tg_yyzzz_yyyy_0, \
                                         tg_yyzzz_yyyy_1, tg_yyzzz_yyyz_0, tg_yyzzz_yyyz_1, tg_yyzzz_yyzz_0, tg_yyzzz_yyzz_1, \
                                         tg_yyzzz_yzzz_0, tg_yyzzz_yzzz_1, tg_yyzzz_zzzz_0, tg_yyzzz_zzzz_1, tg_yzzzz_xxxx_0, \
                                         tg_yzzzz_xxxx_1, tg_yzzzz_xxxy_0, tg_yzzzz_xxxy_1, tg_yzzzz_xxxz_0, tg_yzzzz_xxxz_1, \
                                         tg_yzzzz_xxyy_0, tg_yzzzz_xxyy_1, tg_yzzzz_xxyz_0, tg_yzzzz_xxyz_1, tg_yzzzz_xxzz_0, \
                                         tg_yzzzz_xxzz_1, tg_yzzzz_xyyy_0, tg_yzzzz_xyyy_1, tg_yzzzz_xyyz_0, tg_yzzzz_xyyz_1, \
                                         tg_yzzzz_xyzz_0, tg_yzzzz_xyzz_1, tg_yzzzz_xzzz_0, tg_yzzzz_xzzz_1, tg_yzzzz_yyyy_0, \
                                         tg_yzzzz_yyyy_1, tg_yzzzz_yyyz_0, tg_yzzzz_yyyz_1, tg_yzzzz_yyzz_0, tg_yzzzz_yyzz_1, \
                                         tg_yzzzz_yzzz_0, tg_yzzzz_yzzz_1, tg_yzzzz_zzzz_0, tg_yzzzz_zzzz_1, tg_zzzzz_xxxx_0, \
                                         tg_zzzzz_xxxx_1, tg_zzzzz_xxxy_0, tg_zzzzz_xxxy_1, tg_zzzzz_xxxz_0, tg_zzzzz_xxxz_1, \
                                         tg_zzzzz_xxyy_0, tg_zzzzz_xxyy_1, tg_zzzzz_xxyz_0, tg_zzzzz_xxyz_1, tg_zzzzz_xxzz_0, \
                                         tg_zzzzz_xxzz_1, tg_zzzzz_xyyy_0, tg_zzzzz_xyyy_1, tg_zzzzz_xyyz_0, tg_zzzzz_xyyz_1, \
                                         tg_zzzzz_xyzz_0, tg_zzzzz_xyzz_1, tg_zzzzz_xzzz_0, tg_zzzzz_xzzz_1, tg_zzzzz_yyyy_0, \
                                         tg_zzzzz_yyyy_1, tg_zzzzz_yyyz_0, tg_zzzzz_yyyz_1, tg_zzzzz_yyzz_0, tg_zzzzz_yyzz_1, \
                                         tg_zzzzz_yzzz_0, tg_zzzzz_yzzz_1, tg_zzzzz_zzzz_0, tg_zzzzz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyzzz_xxxx_0[j] = pb_x * tg_xyyzzz_xxxx_0[j] + fr * tg_xyyzzz_xxxx_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxx_0[j] - tg_yyzzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyyzzz_xxx_1[j];

                    tg_xxyyzzz_xxxy_0[j] = pb_x * tg_xyyzzz_xxxy_0[j] + fr * tg_xyyzzz_xxxy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxy_0[j] - tg_yyzzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxy_1[j];

                    tg_xxyyzzz_xxxz_0[j] = pb_x * tg_xyyzzz_xxxz_0[j] + fr * tg_xyyzzz_xxxz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxxz_0[j] - tg_yyzzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xxz_1[j];

                    tg_xxyyzzz_xxyy_0[j] = pb_x * tg_xyyzzz_xxyy_0[j] + fr * tg_xyyzzz_xxyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyy_0[j] - tg_yyzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyy_1[j];

                    tg_xxyyzzz_xxyz_0[j] = pb_x * tg_xyyzzz_xxyz_0[j] + fr * tg_xyyzzz_xxyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxyz_0[j] - tg_yyzzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xyz_1[j];

                    tg_xxyyzzz_xxzz_0[j] = pb_x * tg_xyyzzz_xxzz_0[j] + fr * tg_xyyzzz_xxzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxzz_0[j] - tg_yyzzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xzz_1[j];

                    tg_xxyyzzz_xyyy_0[j] = pb_x * tg_xyyzzz_xyyy_0[j] + fr * tg_xyyzzz_xyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyy_0[j] - tg_yyzzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyy_1[j];

                    tg_xxyyzzz_xyyz_0[j] = pb_x * tg_xyyzzz_xyyz_0[j] + fr * tg_xyyzzz_xyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyyz_0[j] - tg_yyzzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yyz_1[j];

                    tg_xxyyzzz_xyzz_0[j] = pb_x * tg_xyyzzz_xyzz_0[j] + fr * tg_xyyzzz_xyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyzz_0[j] - tg_yyzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yzz_1[j];

                    tg_xxyyzzz_xzzz_0[j] = pb_x * tg_xyyzzz_xzzz_0[j] + fr * tg_xyyzzz_xzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xzzz_0[j] - tg_yyzzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_zzz_1[j];

                    tg_xxyyzzz_yyyy_0[j] = pb_x * tg_xyyzzz_yyyy_0[j] + fr * tg_xyyzzz_yyyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyy_0[j] - tg_yyzzz_yyyy_1[j] * fl1_fza);

                    tg_xxyyzzz_yyyz_0[j] = pb_x * tg_xyyzzz_yyyz_0[j] + fr * tg_xyyzzz_yyyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyyz_0[j] - tg_yyzzz_yyyz_1[j] * fl1_fza);

                    tg_xxyyzzz_yyzz_0[j] = pb_x * tg_xyyzzz_yyzz_0[j] + fr * tg_xyyzzz_yyzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyzz_0[j] - tg_yyzzz_yyzz_1[j] * fl1_fza);

                    tg_xxyyzzz_yzzz_0[j] = pb_x * tg_xyyzzz_yzzz_0[j] + fr * tg_xyyzzz_yzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yzzz_0[j] - tg_yyzzz_yzzz_1[j] * fl1_fza);

                    tg_xxyyzzz_zzzz_0[j] = pb_x * tg_xyyzzz_zzzz_0[j] + fr * tg_xyyzzz_zzzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_zzzz_0[j] - tg_yyzzz_zzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_xxxx_0[j] = pb_x * tg_xyzzzz_xxxx_0[j] + fr * tg_xyzzzz_xxxx_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxx_0[j] - tg_yzzzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xyzzzz_xxx_1[j];

                    tg_xxyzzzz_xxxy_0[j] = pb_x * tg_xyzzzz_xxxy_0[j] + fr * tg_xyzzzz_xxxy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxy_0[j] - tg_yzzzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxy_1[j];

                    tg_xxyzzzz_xxxz_0[j] = pb_x * tg_xyzzzz_xxxz_0[j] + fr * tg_xyzzzz_xxxz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxxz_0[j] - tg_yzzzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xxz_1[j];

                    tg_xxyzzzz_xxyy_0[j] = pb_x * tg_xyzzzz_xxyy_0[j] + fr * tg_xyzzzz_xxyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyy_0[j] - tg_yzzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyy_1[j];

                    tg_xxyzzzz_xxyz_0[j] = pb_x * tg_xyzzzz_xxyz_0[j] + fr * tg_xyzzzz_xxyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxyz_0[j] - tg_yzzzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xyz_1[j];

                    tg_xxyzzzz_xxzz_0[j] = pb_x * tg_xyzzzz_xxzz_0[j] + fr * tg_xyzzzz_xxzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxzz_0[j] - tg_yzzzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xzz_1[j];

                    tg_xxyzzzz_xyyy_0[j] = pb_x * tg_xyzzzz_xyyy_0[j] + fr * tg_xyzzzz_xyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyy_0[j] - tg_yzzzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyy_1[j];

                    tg_xxyzzzz_xyyz_0[j] = pb_x * tg_xyzzzz_xyyz_0[j] + fr * tg_xyzzzz_xyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyyz_0[j] - tg_yzzzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yyz_1[j];

                    tg_xxyzzzz_xyzz_0[j] = pb_x * tg_xyzzzz_xyzz_0[j] + fr * tg_xyzzzz_xyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyzz_0[j] - tg_yzzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yzz_1[j];

                    tg_xxyzzzz_xzzz_0[j] = pb_x * tg_xyzzzz_xzzz_0[j] + fr * tg_xyzzzz_xzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xzzz_0[j] - tg_yzzzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_zzz_1[j];

                    tg_xxyzzzz_yyyy_0[j] = pb_x * tg_xyzzzz_yyyy_0[j] + fr * tg_xyzzzz_yyyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyy_0[j] - tg_yzzzz_yyyy_1[j] * fl1_fza);

                    tg_xxyzzzz_yyyz_0[j] = pb_x * tg_xyzzzz_yyyz_0[j] + fr * tg_xyzzzz_yyyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyyz_0[j] - tg_yzzzz_yyyz_1[j] * fl1_fza);

                    tg_xxyzzzz_yyzz_0[j] = pb_x * tg_xyzzzz_yyzz_0[j] + fr * tg_xyzzzz_yyzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyzz_0[j] - tg_yzzzz_yyzz_1[j] * fl1_fza);

                    tg_xxyzzzz_yzzz_0[j] = pb_x * tg_xyzzzz_yzzz_0[j] + fr * tg_xyzzzz_yzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yzzz_0[j] - tg_yzzzz_yzzz_1[j] * fl1_fza);

                    tg_xxyzzzz_zzzz_0[j] = pb_x * tg_xyzzzz_zzzz_0[j] + fr * tg_xyzzzz_zzzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_zzzz_0[j] - tg_yzzzz_zzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_xxxx_0[j] = pb_x * tg_xzzzzz_xxxx_0[j] + fr * tg_xzzzzz_xxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxx_0[j] - tg_zzzzz_xxxx_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_xzzzzz_xxx_1[j];

                    tg_xxzzzzz_xxxy_0[j] = pb_x * tg_xzzzzz_xxxy_0[j] + fr * tg_xzzzzz_xxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxy_0[j] - tg_zzzzz_xxxy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxy_1[j];

                    tg_xxzzzzz_xxxz_0[j] = pb_x * tg_xzzzzz_xxxz_0[j] + fr * tg_xzzzzz_xxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxz_0[j] - tg_zzzzz_xxxz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xxz_1[j];

                    tg_xxzzzzz_xxyy_0[j] = pb_x * tg_xzzzzz_xxyy_0[j] + fr * tg_xzzzzz_xxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyy_0[j] - tg_zzzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyy_1[j];

                    tg_xxzzzzz_xxyz_0[j] = pb_x * tg_xzzzzz_xxyz_0[j] + fr * tg_xzzzzz_xxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyz_0[j] - tg_zzzzz_xxyz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xyz_1[j];

                    tg_xxzzzzz_xxzz_0[j] = pb_x * tg_xzzzzz_xxzz_0[j] + fr * tg_xzzzzz_xxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxzz_0[j] - tg_zzzzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xzz_1[j];

                    tg_xxzzzzz_xyyy_0[j] = pb_x * tg_xzzzzz_xyyy_0[j] + fr * tg_xzzzzz_xyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyy_0[j] - tg_zzzzz_xyyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyy_1[j];

                    tg_xxzzzzz_xyyz_0[j] = pb_x * tg_xzzzzz_xyyz_0[j] + fr * tg_xzzzzz_xyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyz_0[j] - tg_zzzzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yyz_1[j];

                    tg_xxzzzzz_xyzz_0[j] = pb_x * tg_xzzzzz_xyzz_0[j] + fr * tg_xzzzzz_xyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyzz_0[j] - tg_zzzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yzz_1[j];

                    tg_xxzzzzz_xzzz_0[j] = pb_x * tg_xzzzzz_xzzz_0[j] + fr * tg_xzzzzz_xzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzzz_0[j] - tg_zzzzz_xzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_zzz_1[j];

                    tg_xxzzzzz_yyyy_0[j] = pb_x * tg_xzzzzz_yyyy_0[j] + fr * tg_xzzzzz_yyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyy_0[j] - tg_zzzzz_yyyy_1[j] * fl1_fza);

                    tg_xxzzzzz_yyyz_0[j] = pb_x * tg_xzzzzz_yyyz_0[j] + fr * tg_xzzzzz_yyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyz_0[j] - tg_zzzzz_yyyz_1[j] * fl1_fza);

                    tg_xxzzzzz_yyzz_0[j] = pb_x * tg_xzzzzz_yyzz_0[j] + fr * tg_xzzzzz_yyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyzz_0[j] - tg_zzzzz_yyzz_1[j] * fl1_fza);

                    tg_xxzzzzz_yzzz_0[j] = pb_x * tg_xzzzzz_yzzz_0[j] + fr * tg_xzzzzz_yzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzzz_0[j] - tg_zzzzz_yzzz_1[j] * fl1_fza);

                    tg_xxzzzzz_zzzz_0[j] = pb_x * tg_xzzzzz_zzzz_0[j] + fr * tg_xzzzzz_zzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzzz_0[j] - tg_zzzzz_zzzz_1[j] * fl1_fza);

                    tg_xyyyyyy_xxxx_0[j] = pb_x * tg_yyyyyy_xxxx_0[j] + fr * tg_yyyyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyyy_xxx_1[j];

                    tg_xyyyyyy_xxxy_0[j] = pb_x * tg_yyyyyy_xxxy_0[j] + fr * tg_yyyyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxy_1[j];

                    tg_xyyyyyy_xxxz_0[j] = pb_x * tg_yyyyyy_xxxz_0[j] + fr * tg_yyyyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xxz_1[j];

                    tg_xyyyyyy_xxyy_0[j] = pb_x * tg_yyyyyy_xxyy_0[j] + fr * tg_yyyyyy_xxyy_1[j] + fl1_fxn * tg_yyyyyy_xyy_1[j];

                    tg_xyyyyyy_xxyz_0[j] = pb_x * tg_yyyyyy_xxyz_0[j] + fr * tg_yyyyyy_xxyz_1[j] + fl1_fxn * tg_yyyyyy_xyz_1[j];

                    tg_xyyyyyy_xxzz_0[j] = pb_x * tg_yyyyyy_xxzz_0[j] + fr * tg_yyyyyy_xxzz_1[j] + fl1_fxn * tg_yyyyyy_xzz_1[j];

                    tg_xyyyyyy_xyyy_0[j] = pb_x * tg_yyyyyy_xyyy_0[j] + fr * tg_yyyyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyy_1[j];

                    tg_xyyyyyy_xyyz_0[j] = pb_x * tg_yyyyyy_xyyz_0[j] + fr * tg_yyyyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yyz_1[j];

                    tg_xyyyyyy_xyzz_0[j] = pb_x * tg_yyyyyy_xyzz_0[j] + fr * tg_yyyyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yzz_1[j];

                    tg_xyyyyyy_xzzz_0[j] = pb_x * tg_yyyyyy_xzzz_0[j] + fr * tg_yyyyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_zzz_1[j];

                    tg_xyyyyyy_yyyy_0[j] = pb_x * tg_yyyyyy_yyyy_0[j] + fr * tg_yyyyyy_yyyy_1[j];

                    tg_xyyyyyy_yyyz_0[j] = pb_x * tg_yyyyyy_yyyz_0[j] + fr * tg_yyyyyy_yyyz_1[j];

                    tg_xyyyyyy_yyzz_0[j] = pb_x * tg_yyyyyy_yyzz_0[j] + fr * tg_yyyyyy_yyzz_1[j];

                    tg_xyyyyyy_yzzz_0[j] = pb_x * tg_yyyyyy_yzzz_0[j] + fr * tg_yyyyyy_yzzz_1[j];

                    tg_xyyyyyy_zzzz_0[j] = pb_x * tg_yyyyyy_zzzz_0[j] + fr * tg_yyyyyy_zzzz_1[j];

                    tg_xyyyyyz_xxxx_0[j] = pb_x * tg_yyyyyz_xxxx_0[j] + fr * tg_yyyyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyyz_xxx_1[j];

                    tg_xyyyyyz_xxxy_0[j] = pb_x * tg_yyyyyz_xxxy_0[j] + fr * tg_yyyyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxy_1[j];

                    tg_xyyyyyz_xxxz_0[j] = pb_x * tg_yyyyyz_xxxz_0[j] + fr * tg_yyyyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xxz_1[j];

                    tg_xyyyyyz_xxyy_0[j] = pb_x * tg_yyyyyz_xxyy_0[j] + fr * tg_yyyyyz_xxyy_1[j] + fl1_fxn * tg_yyyyyz_xyy_1[j];

                    tg_xyyyyyz_xxyz_0[j] = pb_x * tg_yyyyyz_xxyz_0[j] + fr * tg_yyyyyz_xxyz_1[j] + fl1_fxn * tg_yyyyyz_xyz_1[j];

                    tg_xyyyyyz_xxzz_0[j] = pb_x * tg_yyyyyz_xxzz_0[j] + fr * tg_yyyyyz_xxzz_1[j] + fl1_fxn * tg_yyyyyz_xzz_1[j];

                    tg_xyyyyyz_xyyy_0[j] = pb_x * tg_yyyyyz_xyyy_0[j] + fr * tg_yyyyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyy_1[j];

                    tg_xyyyyyz_xyyz_0[j] = pb_x * tg_yyyyyz_xyyz_0[j] + fr * tg_yyyyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yyz_1[j];

                    tg_xyyyyyz_xyzz_0[j] = pb_x * tg_yyyyyz_xyzz_0[j] + fr * tg_yyyyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yzz_1[j];

                    tg_xyyyyyz_xzzz_0[j] = pb_x * tg_yyyyyz_xzzz_0[j] + fr * tg_yyyyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_zzz_1[j];

                    tg_xyyyyyz_yyyy_0[j] = pb_x * tg_yyyyyz_yyyy_0[j] + fr * tg_yyyyyz_yyyy_1[j];

                    tg_xyyyyyz_yyyz_0[j] = pb_x * tg_yyyyyz_yyyz_0[j] + fr * tg_yyyyyz_yyyz_1[j];

                    tg_xyyyyyz_yyzz_0[j] = pb_x * tg_yyyyyz_yyzz_0[j] + fr * tg_yyyyyz_yyzz_1[j];

                    tg_xyyyyyz_yzzz_0[j] = pb_x * tg_yyyyyz_yzzz_0[j] + fr * tg_yyyyyz_yzzz_1[j];

                    tg_xyyyyyz_zzzz_0[j] = pb_x * tg_yyyyyz_zzzz_0[j] + fr * tg_yyyyyz_zzzz_1[j];

                    tg_xyyyyzz_xxxx_0[j] = pb_x * tg_yyyyzz_xxxx_0[j] + fr * tg_yyyyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyyzz_xxx_1[j];

                    tg_xyyyyzz_xxxy_0[j] = pb_x * tg_yyyyzz_xxxy_0[j] + fr * tg_yyyyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxy_1[j];

                    tg_xyyyyzz_xxxz_0[j] = pb_x * tg_yyyyzz_xxxz_0[j] + fr * tg_yyyyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xxz_1[j];

                    tg_xyyyyzz_xxyy_0[j] = pb_x * tg_yyyyzz_xxyy_0[j] + fr * tg_yyyyzz_xxyy_1[j] + fl1_fxn * tg_yyyyzz_xyy_1[j];

                    tg_xyyyyzz_xxyz_0[j] = pb_x * tg_yyyyzz_xxyz_0[j] + fr * tg_yyyyzz_xxyz_1[j] + fl1_fxn * tg_yyyyzz_xyz_1[j];

                    tg_xyyyyzz_xxzz_0[j] = pb_x * tg_yyyyzz_xxzz_0[j] + fr * tg_yyyyzz_xxzz_1[j] + fl1_fxn * tg_yyyyzz_xzz_1[j];

                    tg_xyyyyzz_xyyy_0[j] = pb_x * tg_yyyyzz_xyyy_0[j] + fr * tg_yyyyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyy_1[j];

                    tg_xyyyyzz_xyyz_0[j] = pb_x * tg_yyyyzz_xyyz_0[j] + fr * tg_yyyyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yyz_1[j];

                    tg_xyyyyzz_xyzz_0[j] = pb_x * tg_yyyyzz_xyzz_0[j] + fr * tg_yyyyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yzz_1[j];

                    tg_xyyyyzz_xzzz_0[j] = pb_x * tg_yyyyzz_xzzz_0[j] + fr * tg_yyyyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_zzz_1[j];

                    tg_xyyyyzz_yyyy_0[j] = pb_x * tg_yyyyzz_yyyy_0[j] + fr * tg_yyyyzz_yyyy_1[j];

                    tg_xyyyyzz_yyyz_0[j] = pb_x * tg_yyyyzz_yyyz_0[j] + fr * tg_yyyyzz_yyyz_1[j];

                    tg_xyyyyzz_yyzz_0[j] = pb_x * tg_yyyyzz_yyzz_0[j] + fr * tg_yyyyzz_yyzz_1[j];

                    tg_xyyyyzz_yzzz_0[j] = pb_x * tg_yyyyzz_yzzz_0[j] + fr * tg_yyyyzz_yzzz_1[j];

                    tg_xyyyyzz_zzzz_0[j] = pb_x * tg_yyyyzz_zzzz_0[j] + fr * tg_yyyyzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSG_360_450(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (360,450)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_yyyyyy_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 315); 

                auto tg_yyyyyy_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 316); 

                auto tg_yyyyyy_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 317); 

                auto tg_yyyyyy_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 318); 

                auto tg_yyyyyy_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 319); 

                auto tg_yyyyyy_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 320); 

                auto tg_yyyyyy_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 321); 

                auto tg_yyyyyy_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 322); 

                auto tg_yyyyyy_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 323); 

                auto tg_yyyyyy_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 324); 

                auto tg_yyyyyy_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 325); 

                auto tg_yyyyyy_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 326); 

                auto tg_yyyyyy_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 327); 

                auto tg_yyyyyy_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 328); 

                auto tg_yyyyyy_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 329); 

                auto tg_yyyyyz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 330); 

                auto tg_yyyyyz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 331); 

                auto tg_yyyyyz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 332); 

                auto tg_yyyyyz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 333); 

                auto tg_yyyyyz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 334); 

                auto tg_yyyyyz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 335); 

                auto tg_yyyyyz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 336); 

                auto tg_yyyyyz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 337); 

                auto tg_yyyyyz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 338); 

                auto tg_yyyyyz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 339); 

                auto tg_yyyyyz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 340); 

                auto tg_yyyyyz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 341); 

                auto tg_yyyyyz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 342); 

                auto tg_yyyyyz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 343); 

                auto tg_yyyyyz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 344); 

                auto tg_yyyzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 360); 

                auto tg_yyyzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 361); 

                auto tg_yyyzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 362); 

                auto tg_yyyzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 363); 

                auto tg_yyyzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 364); 

                auto tg_yyyzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 365); 

                auto tg_yyyzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 366); 

                auto tg_yyyzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 367); 

                auto tg_yyyzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 368); 

                auto tg_yyyzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 369); 

                auto tg_yyyzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 370); 

                auto tg_yyyzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 371); 

                auto tg_yyyzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 372); 

                auto tg_yyyzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 373); 

                auto tg_yyyzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 374); 

                auto tg_yyzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 375); 

                auto tg_yyzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 376); 

                auto tg_yyzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 377); 

                auto tg_yyzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 378); 

                auto tg_yyzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 379); 

                auto tg_yyzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 380); 

                auto tg_yyzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 381); 

                auto tg_yyzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 382); 

                auto tg_yyzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 383); 

                auto tg_yyzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 384); 

                auto tg_yyzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 385); 

                auto tg_yyzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 386); 

                auto tg_yyzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 387); 

                auto tg_yyzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 388); 

                auto tg_yyzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 389); 

                auto tg_yzzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 390); 

                auto tg_yzzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 391); 

                auto tg_yzzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 392); 

                auto tg_yzzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 393); 

                auto tg_yzzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 394); 

                auto tg_yzzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 395); 

                auto tg_yzzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 396); 

                auto tg_yzzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 397); 

                auto tg_yzzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 398); 

                auto tg_yzzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 399); 

                auto tg_yzzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 400); 

                auto tg_yzzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 401); 

                auto tg_yzzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 402); 

                auto tg_yzzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 403); 

                auto tg_yzzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 404); 

                auto tg_zzzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 405); 

                auto tg_zzzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 406); 

                auto tg_zzzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 407); 

                auto tg_zzzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 408); 

                auto tg_zzzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 409); 

                auto tg_zzzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 410); 

                auto tg_zzzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 411); 

                auto tg_zzzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 412); 

                auto tg_zzzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 413); 

                auto tg_zzzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 414); 

                auto tg_zzzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 415); 

                auto tg_zzzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 416); 

                auto tg_zzzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 417); 

                auto tg_zzzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 418); 

                auto tg_zzzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 419); 

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

                auto tg_yyyyy_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 225); 

                auto tg_yyyyy_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 226); 

                auto tg_yyyyy_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 227); 

                auto tg_yyyyy_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 228); 

                auto tg_yyyyy_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 229); 

                auto tg_yyyyy_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 230); 

                auto tg_yyyyy_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 231); 

                auto tg_yyyyy_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 232); 

                auto tg_yyyyy_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 233); 

                auto tg_yyyyy_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 234); 

                auto tg_yyyyy_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 235); 

                auto tg_yyyyy_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 236); 

                auto tg_yyyyy_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 237); 

                auto tg_yyyyy_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 238); 

                auto tg_yyyyy_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 239); 

                auto tg_yyyyz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 240); 

                auto tg_yyyyz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 241); 

                auto tg_yyyyz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 242); 

                auto tg_yyyyz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 243); 

                auto tg_yyyyz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 244); 

                auto tg_yyyyz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 245); 

                auto tg_yyyyz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 246); 

                auto tg_yyyyz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 247); 

                auto tg_yyyyz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 248); 

                auto tg_yyyyz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 249); 

                auto tg_yyyyz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 250); 

                auto tg_yyyyz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 251); 

                auto tg_yyyyz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 252); 

                auto tg_yyyyz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 253); 

                auto tg_yyyyz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 254); 

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

                auto tg_yyyyyy_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 210); 

                auto tg_yyyyyy_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 211); 

                auto tg_yyyyyy_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 212); 

                auto tg_yyyyyy_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 213); 

                auto tg_yyyyyy_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 214); 

                auto tg_yyyyyy_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 215); 

                auto tg_yyyyyy_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 216); 

                auto tg_yyyyyy_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 217); 

                auto tg_yyyyyy_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 218); 

                auto tg_yyyyyy_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 219); 

                auto tg_yyyyyz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 220); 

                auto tg_yyyyyz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 221); 

                auto tg_yyyyyz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 222); 

                auto tg_yyyyyz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 223); 

                auto tg_yyyyyz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 224); 

                auto tg_yyyyyz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 225); 

                auto tg_yyyyyz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 226); 

                auto tg_yyyyyz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 227); 

                auto tg_yyyyyz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 228); 

                auto tg_yyyyyz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 229); 

                auto tg_yyyzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 240); 

                auto tg_yyyzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 241); 

                auto tg_yyyzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 242); 

                auto tg_yyyzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 243); 

                auto tg_yyyzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 244); 

                auto tg_yyyzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 245); 

                auto tg_yyyzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 246); 

                auto tg_yyyzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 247); 

                auto tg_yyyzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 248); 

                auto tg_yyyzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 249); 

                auto tg_yyzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 250); 

                auto tg_yyzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 251); 

                auto tg_yyzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 252); 

                auto tg_yyzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 253); 

                auto tg_yyzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 254); 

                auto tg_yyzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 255); 

                auto tg_yyzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 256); 

                auto tg_yyzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 257); 

                auto tg_yyzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 258); 

                auto tg_yyzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 259); 

                auto tg_yzzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 260); 

                auto tg_yzzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 261); 

                auto tg_yzzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 262); 

                auto tg_yzzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 263); 

                auto tg_yzzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 264); 

                auto tg_yzzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 265); 

                auto tg_yzzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 266); 

                auto tg_yzzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 267); 

                auto tg_yzzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 268); 

                auto tg_yzzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 269); 

                auto tg_zzzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 270); 

                auto tg_zzzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 271); 

                auto tg_zzzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 272); 

                auto tg_zzzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 273); 

                auto tg_zzzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 274); 

                auto tg_zzzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 275); 

                auto tg_zzzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 276); 

                auto tg_zzzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 277); 

                auto tg_zzzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 278); 

                auto tg_zzzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 279); 

                // set up pointers to integrals

                auto tg_xyyyzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 360); 

                auto tg_xyyyzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 361); 

                auto tg_xyyyzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 362); 

                auto tg_xyyyzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 363); 

                auto tg_xyyyzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 364); 

                auto tg_xyyyzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 365); 

                auto tg_xyyyzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 366); 

                auto tg_xyyyzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 367); 

                auto tg_xyyyzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 368); 

                auto tg_xyyyzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 369); 

                auto tg_xyyyzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 370); 

                auto tg_xyyyzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 371); 

                auto tg_xyyyzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 372); 

                auto tg_xyyyzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 373); 

                auto tg_xyyyzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 374); 

                auto tg_xyyzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 375); 

                auto tg_xyyzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 376); 

                auto tg_xyyzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 377); 

                auto tg_xyyzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 378); 

                auto tg_xyyzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 379); 

                auto tg_xyyzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 380); 

                auto tg_xyyzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 381); 

                auto tg_xyyzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 382); 

                auto tg_xyyzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 383); 

                auto tg_xyyzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 384); 

                auto tg_xyyzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 385); 

                auto tg_xyyzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 386); 

                auto tg_xyyzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 387); 

                auto tg_xyyzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 388); 

                auto tg_xyyzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 389); 

                auto tg_xyzzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 390); 

                auto tg_xyzzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 391); 

                auto tg_xyzzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 392); 

                auto tg_xyzzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 393); 

                auto tg_xyzzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 394); 

                auto tg_xyzzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 395); 

                auto tg_xyzzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 396); 

                auto tg_xyzzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 397); 

                auto tg_xyzzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 398); 

                auto tg_xyzzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 399); 

                auto tg_xyzzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 400); 

                auto tg_xyzzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 401); 

                auto tg_xyzzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 402); 

                auto tg_xyzzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 403); 

                auto tg_xyzzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 404); 

                auto tg_xzzzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 405); 

                auto tg_xzzzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 406); 

                auto tg_xzzzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 407); 

                auto tg_xzzzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 408); 

                auto tg_xzzzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 409); 

                auto tg_xzzzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 410); 

                auto tg_xzzzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 411); 

                auto tg_xzzzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 412); 

                auto tg_xzzzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 413); 

                auto tg_xzzzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 414); 

                auto tg_xzzzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 415); 

                auto tg_xzzzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 416); 

                auto tg_xzzzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 417); 

                auto tg_xzzzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 418); 

                auto tg_xzzzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 419); 

                auto tg_yyyyyyy_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 420); 

                auto tg_yyyyyyy_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 421); 

                auto tg_yyyyyyy_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 422); 

                auto tg_yyyyyyy_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 423); 

                auto tg_yyyyyyy_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 424); 

                auto tg_yyyyyyy_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 425); 

                auto tg_yyyyyyy_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 426); 

                auto tg_yyyyyyy_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 427); 

                auto tg_yyyyyyy_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 428); 

                auto tg_yyyyyyy_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 429); 

                auto tg_yyyyyyy_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 430); 

                auto tg_yyyyyyy_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 431); 

                auto tg_yyyyyyy_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 432); 

                auto tg_yyyyyyy_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 433); 

                auto tg_yyyyyyy_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 434); 

                auto tg_yyyyyyz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 435); 

                auto tg_yyyyyyz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 436); 

                auto tg_yyyyyyz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 437); 

                auto tg_yyyyyyz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 438); 

                auto tg_yyyyyyz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 439); 

                auto tg_yyyyyyz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 440); 

                auto tg_yyyyyyz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 441); 

                auto tg_yyyyyyz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 442); 

                auto tg_yyyyyyz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 443); 

                auto tg_yyyyyyz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 444); 

                auto tg_yyyyyyz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 445); 

                auto tg_yyyyyyz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 446); 

                auto tg_yyyyyyz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 447); 

                auto tg_yyyyyyz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 448); 

                auto tg_yyyyyyz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 449); 

                // Batch of Integrals (360,450)

                #pragma omp simd aligned(fxn, fza, tg_xyyyzzz_xxxx_0, tg_xyyyzzz_xxxy_0, tg_xyyyzzz_xxxz_0, \
                                         tg_xyyyzzz_xxyy_0, tg_xyyyzzz_xxyz_0, tg_xyyyzzz_xxzz_0, tg_xyyyzzz_xyyy_0, \
                                         tg_xyyyzzz_xyyz_0, tg_xyyyzzz_xyzz_0, tg_xyyyzzz_xzzz_0, tg_xyyyzzz_yyyy_0, \
                                         tg_xyyyzzz_yyyz_0, tg_xyyyzzz_yyzz_0, tg_xyyyzzz_yzzz_0, tg_xyyyzzz_zzzz_0, \
                                         tg_xyyzzzz_xxxx_0, tg_xyyzzzz_xxxy_0, tg_xyyzzzz_xxxz_0, tg_xyyzzzz_xxyy_0, \
                                         tg_xyyzzzz_xxyz_0, tg_xyyzzzz_xxzz_0, tg_xyyzzzz_xyyy_0, tg_xyyzzzz_xyyz_0, \
                                         tg_xyyzzzz_xyzz_0, tg_xyyzzzz_xzzz_0, tg_xyyzzzz_yyyy_0, tg_xyyzzzz_yyyz_0, \
                                         tg_xyyzzzz_yyzz_0, tg_xyyzzzz_yzzz_0, tg_xyyzzzz_zzzz_0, tg_xyzzzzz_xxxx_0, \
                                         tg_xyzzzzz_xxxy_0, tg_xyzzzzz_xxxz_0, tg_xyzzzzz_xxyy_0, tg_xyzzzzz_xxyz_0, \
                                         tg_xyzzzzz_xxzz_0, tg_xyzzzzz_xyyy_0, tg_xyzzzzz_xyyz_0, tg_xyzzzzz_xyzz_0, \
                                         tg_xyzzzzz_xzzz_0, tg_xyzzzzz_yyyy_0, tg_xyzzzzz_yyyz_0, tg_xyzzzzz_yyzz_0, \
                                         tg_xyzzzzz_yzzz_0, tg_xyzzzzz_zzzz_0, tg_xzzzzzz_xxxx_0, tg_xzzzzzz_xxxy_0, \
                                         tg_xzzzzzz_xxxz_0, tg_xzzzzzz_xxyy_0, tg_xzzzzzz_xxyz_0, tg_xzzzzzz_xxzz_0, \
                                         tg_xzzzzzz_xyyy_0, tg_xzzzzzz_xyyz_0, tg_xzzzzzz_xyzz_0, tg_xzzzzzz_xzzz_0, \
                                         tg_xzzzzzz_yyyy_0, tg_xzzzzzz_yyyz_0, tg_xzzzzzz_yyzz_0, tg_xzzzzzz_yzzz_0, \
                                         tg_xzzzzzz_zzzz_0, tg_yyyyy_xxxx_0, tg_yyyyy_xxxx_1, tg_yyyyy_xxxy_0, tg_yyyyy_xxxy_1, \
                                         tg_yyyyy_xxxz_0, tg_yyyyy_xxxz_1, tg_yyyyy_xxyy_0, tg_yyyyy_xxyy_1, tg_yyyyy_xxyz_0, \
                                         tg_yyyyy_xxyz_1, tg_yyyyy_xxzz_0, tg_yyyyy_xxzz_1, tg_yyyyy_xyyy_0, tg_yyyyy_xyyy_1, \
                                         tg_yyyyy_xyyz_0, tg_yyyyy_xyyz_1, tg_yyyyy_xyzz_0, tg_yyyyy_xyzz_1, tg_yyyyy_xzzz_0, \
                                         tg_yyyyy_xzzz_1, tg_yyyyy_yyyy_0, tg_yyyyy_yyyy_1, tg_yyyyy_yyyz_0, tg_yyyyy_yyyz_1, \
                                         tg_yyyyy_yyzz_0, tg_yyyyy_yyzz_1, tg_yyyyy_yzzz_0, tg_yyyyy_yzzz_1, tg_yyyyy_zzzz_0, \
                                         tg_yyyyy_zzzz_1, tg_yyyyyy_xxx_1, tg_yyyyyy_xxxx_0, tg_yyyyyy_xxxx_1, \
                                         tg_yyyyyy_xxxy_0, tg_yyyyyy_xxxy_1, tg_yyyyyy_xxxz_0, tg_yyyyyy_xxxz_1, \
                                         tg_yyyyyy_xxy_1, tg_yyyyyy_xxyy_0, tg_yyyyyy_xxyy_1, tg_yyyyyy_xxyz_0, \
                                         tg_yyyyyy_xxyz_1, tg_yyyyyy_xxz_1, tg_yyyyyy_xxzz_0, tg_yyyyyy_xxzz_1, \
                                         tg_yyyyyy_xyy_1, tg_yyyyyy_xyyy_0, tg_yyyyyy_xyyy_1, tg_yyyyyy_xyyz_0, \
                                         tg_yyyyyy_xyyz_1, tg_yyyyyy_xyz_1, tg_yyyyyy_xyzz_0, tg_yyyyyy_xyzz_1, \
                                         tg_yyyyyy_xzz_1, tg_yyyyyy_xzzz_0, tg_yyyyyy_xzzz_1, tg_yyyyyy_yyy_1, \
                                         tg_yyyyyy_yyyy_0, tg_yyyyyy_yyyy_1, tg_yyyyyy_yyyz_0, tg_yyyyyy_yyyz_1, \
                                         tg_yyyyyy_yyz_1, tg_yyyyyy_yyzz_0, tg_yyyyyy_yyzz_1, tg_yyyyyy_yzz_1, \
                                         tg_yyyyyy_yzzz_0, tg_yyyyyy_yzzz_1, tg_yyyyyy_zzz_1, tg_yyyyyy_zzzz_0, \
                                         tg_yyyyyy_zzzz_1, tg_yyyyyyy_xxxx_0, tg_yyyyyyy_xxxy_0, tg_yyyyyyy_xxxz_0, \
                                         tg_yyyyyyy_xxyy_0, tg_yyyyyyy_xxyz_0, tg_yyyyyyy_xxzz_0, tg_yyyyyyy_xyyy_0, \
                                         tg_yyyyyyy_xyyz_0, tg_yyyyyyy_xyzz_0, tg_yyyyyyy_xzzz_0, tg_yyyyyyy_yyyy_0, \
                                         tg_yyyyyyy_yyyz_0, tg_yyyyyyy_yyzz_0, tg_yyyyyyy_yzzz_0, tg_yyyyyyy_zzzz_0, \
                                         tg_yyyyyyz_xxxx_0, tg_yyyyyyz_xxxy_0, tg_yyyyyyz_xxxz_0, tg_yyyyyyz_xxyy_0, \
                                         tg_yyyyyyz_xxyz_0, tg_yyyyyyz_xxzz_0, tg_yyyyyyz_xyyy_0, tg_yyyyyyz_xyyz_0, \
                                         tg_yyyyyyz_xyzz_0, tg_yyyyyyz_xzzz_0, tg_yyyyyyz_yyyy_0, tg_yyyyyyz_yyyz_0, \
                                         tg_yyyyyyz_yyzz_0, tg_yyyyyyz_yzzz_0, tg_yyyyyyz_zzzz_0, tg_yyyyyz_xxx_1, \
                                         tg_yyyyyz_xxxx_0, tg_yyyyyz_xxxx_1, tg_yyyyyz_xxxy_0, tg_yyyyyz_xxxy_1, \
                                         tg_yyyyyz_xxxz_0, tg_yyyyyz_xxxz_1, tg_yyyyyz_xxy_1, tg_yyyyyz_xxyy_0, \
                                         tg_yyyyyz_xxyy_1, tg_yyyyyz_xxyz_0, tg_yyyyyz_xxyz_1, tg_yyyyyz_xxz_1, \
                                         tg_yyyyyz_xxzz_0, tg_yyyyyz_xxzz_1, tg_yyyyyz_xyy_1, tg_yyyyyz_xyyy_0, \
                                         tg_yyyyyz_xyyy_1, tg_yyyyyz_xyyz_0, tg_yyyyyz_xyyz_1, tg_yyyyyz_xyz_1, \
                                         tg_yyyyyz_xyzz_0, tg_yyyyyz_xyzz_1, tg_yyyyyz_xzz_1, tg_yyyyyz_xzzz_0, \
                                         tg_yyyyyz_xzzz_1, tg_yyyyyz_yyy_1, tg_yyyyyz_yyyy_0, tg_yyyyyz_yyyy_1, \
                                         tg_yyyyyz_yyyz_0, tg_yyyyyz_yyyz_1, tg_yyyyyz_yyz_1, tg_yyyyyz_yyzz_0, \
                                         tg_yyyyyz_yyzz_1, tg_yyyyyz_yzz_1, tg_yyyyyz_yzzz_0, tg_yyyyyz_yzzz_1, \
                                         tg_yyyyyz_zzz_1, tg_yyyyyz_zzzz_0, tg_yyyyyz_zzzz_1, tg_yyyyz_xxxx_0, \
                                         tg_yyyyz_xxxx_1, tg_yyyyz_xxxy_0, tg_yyyyz_xxxy_1, tg_yyyyz_xxxz_0, tg_yyyyz_xxxz_1, \
                                         tg_yyyyz_xxyy_0, tg_yyyyz_xxyy_1, tg_yyyyz_xxyz_0, tg_yyyyz_xxyz_1, tg_yyyyz_xxzz_0, \
                                         tg_yyyyz_xxzz_1, tg_yyyyz_xyyy_0, tg_yyyyz_xyyy_1, tg_yyyyz_xyyz_0, tg_yyyyz_xyyz_1, \
                                         tg_yyyyz_xyzz_0, tg_yyyyz_xyzz_1, tg_yyyyz_xzzz_0, tg_yyyyz_xzzz_1, tg_yyyyz_yyyy_0, \
                                         tg_yyyyz_yyyy_1, tg_yyyyz_yyyz_0, tg_yyyyz_yyyz_1, tg_yyyyz_yyzz_0, tg_yyyyz_yyzz_1, \
                                         tg_yyyyz_yzzz_0, tg_yyyyz_yzzz_1, tg_yyyyz_zzzz_0, tg_yyyyz_zzzz_1, tg_yyyzzz_xxx_1, \
                                         tg_yyyzzz_xxxx_0, tg_yyyzzz_xxxx_1, tg_yyyzzz_xxxy_0, tg_yyyzzz_xxxy_1, \
                                         tg_yyyzzz_xxxz_0, tg_yyyzzz_xxxz_1, tg_yyyzzz_xxy_1, tg_yyyzzz_xxyy_0, \
                                         tg_yyyzzz_xxyy_1, tg_yyyzzz_xxyz_0, tg_yyyzzz_xxyz_1, tg_yyyzzz_xxz_1, \
                                         tg_yyyzzz_xxzz_0, tg_yyyzzz_xxzz_1, tg_yyyzzz_xyy_1, tg_yyyzzz_xyyy_0, \
                                         tg_yyyzzz_xyyy_1, tg_yyyzzz_xyyz_0, tg_yyyzzz_xyyz_1, tg_yyyzzz_xyz_1, \
                                         tg_yyyzzz_xyzz_0, tg_yyyzzz_xyzz_1, tg_yyyzzz_xzz_1, tg_yyyzzz_xzzz_0, \
                                         tg_yyyzzz_xzzz_1, tg_yyyzzz_yyy_1, tg_yyyzzz_yyyy_0, tg_yyyzzz_yyyy_1, \
                                         tg_yyyzzz_yyyz_0, tg_yyyzzz_yyyz_1, tg_yyyzzz_yyz_1, tg_yyyzzz_yyzz_0, \
                                         tg_yyyzzz_yyzz_1, tg_yyyzzz_yzz_1, tg_yyyzzz_yzzz_0, tg_yyyzzz_yzzz_1, \
                                         tg_yyyzzz_zzz_1, tg_yyyzzz_zzzz_0, tg_yyyzzz_zzzz_1, tg_yyzzzz_xxx_1, \
                                         tg_yyzzzz_xxxx_0, tg_yyzzzz_xxxx_1, tg_yyzzzz_xxxy_0, tg_yyzzzz_xxxy_1, \
                                         tg_yyzzzz_xxxz_0, tg_yyzzzz_xxxz_1, tg_yyzzzz_xxy_1, tg_yyzzzz_xxyy_0, \
                                         tg_yyzzzz_xxyy_1, tg_yyzzzz_xxyz_0, tg_yyzzzz_xxyz_1, tg_yyzzzz_xxz_1, \
                                         tg_yyzzzz_xxzz_0, tg_yyzzzz_xxzz_1, tg_yyzzzz_xyy_1, tg_yyzzzz_xyyy_0, \
                                         tg_yyzzzz_xyyy_1, tg_yyzzzz_xyyz_0, tg_yyzzzz_xyyz_1, tg_yyzzzz_xyz_1, \
                                         tg_yyzzzz_xyzz_0, tg_yyzzzz_xyzz_1, tg_yyzzzz_xzz_1, tg_yyzzzz_xzzz_0, \
                                         tg_yyzzzz_xzzz_1, tg_yyzzzz_yyy_1, tg_yyzzzz_yyyy_0, tg_yyzzzz_yyyy_1, \
                                         tg_yyzzzz_yyyz_0, tg_yyzzzz_yyyz_1, tg_yyzzzz_yyz_1, tg_yyzzzz_yyzz_0, \
                                         tg_yyzzzz_yyzz_1, tg_yyzzzz_yzz_1, tg_yyzzzz_yzzz_0, tg_yyzzzz_yzzz_1, \
                                         tg_yyzzzz_zzz_1, tg_yyzzzz_zzzz_0, tg_yyzzzz_zzzz_1, tg_yzzzzz_xxx_1, \
                                         tg_yzzzzz_xxxx_0, tg_yzzzzz_xxxx_1, tg_yzzzzz_xxxy_0, tg_yzzzzz_xxxy_1, \
                                         tg_yzzzzz_xxxz_0, tg_yzzzzz_xxxz_1, tg_yzzzzz_xxy_1, tg_yzzzzz_xxyy_0, \
                                         tg_yzzzzz_xxyy_1, tg_yzzzzz_xxyz_0, tg_yzzzzz_xxyz_1, tg_yzzzzz_xxz_1, \
                                         tg_yzzzzz_xxzz_0, tg_yzzzzz_xxzz_1, tg_yzzzzz_xyy_1, tg_yzzzzz_xyyy_0, \
                                         tg_yzzzzz_xyyy_1, tg_yzzzzz_xyyz_0, tg_yzzzzz_xyyz_1, tg_yzzzzz_xyz_1, \
                                         tg_yzzzzz_xyzz_0, tg_yzzzzz_xyzz_1, tg_yzzzzz_xzz_1, tg_yzzzzz_xzzz_0, \
                                         tg_yzzzzz_xzzz_1, tg_yzzzzz_yyy_1, tg_yzzzzz_yyyy_0, tg_yzzzzz_yyyy_1, \
                                         tg_yzzzzz_yyyz_0, tg_yzzzzz_yyyz_1, tg_yzzzzz_yyz_1, tg_yzzzzz_yyzz_0, \
                                         tg_yzzzzz_yyzz_1, tg_yzzzzz_yzz_1, tg_yzzzzz_yzzz_0, tg_yzzzzz_yzzz_1, \
                                         tg_yzzzzz_zzz_1, tg_yzzzzz_zzzz_0, tg_yzzzzz_zzzz_1, tg_zzzzzz_xxx_1, \
                                         tg_zzzzzz_xxxx_0, tg_zzzzzz_xxxx_1, tg_zzzzzz_xxxy_0, tg_zzzzzz_xxxy_1, \
                                         tg_zzzzzz_xxxz_0, tg_zzzzzz_xxxz_1, tg_zzzzzz_xxy_1, tg_zzzzzz_xxyy_0, \
                                         tg_zzzzzz_xxyy_1, tg_zzzzzz_xxyz_0, tg_zzzzzz_xxyz_1, tg_zzzzzz_xxz_1, \
                                         tg_zzzzzz_xxzz_0, tg_zzzzzz_xxzz_1, tg_zzzzzz_xyy_1, tg_zzzzzz_xyyy_0, \
                                         tg_zzzzzz_xyyy_1, tg_zzzzzz_xyyz_0, tg_zzzzzz_xyyz_1, tg_zzzzzz_xyz_1, \
                                         tg_zzzzzz_xyzz_0, tg_zzzzzz_xyzz_1, tg_zzzzzz_xzz_1, tg_zzzzzz_xzzz_0, \
                                         tg_zzzzzz_xzzz_1, tg_zzzzzz_yyy_1, tg_zzzzzz_yyyy_0, tg_zzzzzz_yyyy_1, \
                                         tg_zzzzzz_yyyz_0, tg_zzzzzz_yyyz_1, tg_zzzzzz_yyz_1, tg_zzzzzz_yyzz_0, \
                                         tg_zzzzzz_yyzz_1, tg_zzzzzz_yzz_1, tg_zzzzzz_yzzz_0, tg_zzzzzz_yzzz_1, \
                                         tg_zzzzzz_zzz_1, tg_zzzzzz_zzzz_0, tg_zzzzzz_zzzz_1, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyyyzzz_xxxx_0[j] = pb_x * tg_yyyzzz_xxxx_0[j] + fr * tg_yyyzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyzzz_xxx_1[j];

                    tg_xyyyzzz_xxxy_0[j] = pb_x * tg_yyyzzz_xxxy_0[j] + fr * tg_yyyzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxy_1[j];

                    tg_xyyyzzz_xxxz_0[j] = pb_x * tg_yyyzzz_xxxz_0[j] + fr * tg_yyyzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xxz_1[j];

                    tg_xyyyzzz_xxyy_0[j] = pb_x * tg_yyyzzz_xxyy_0[j] + fr * tg_yyyzzz_xxyy_1[j] + fl1_fxn * tg_yyyzzz_xyy_1[j];

                    tg_xyyyzzz_xxyz_0[j] = pb_x * tg_yyyzzz_xxyz_0[j] + fr * tg_yyyzzz_xxyz_1[j] + fl1_fxn * tg_yyyzzz_xyz_1[j];

                    tg_xyyyzzz_xxzz_0[j] = pb_x * tg_yyyzzz_xxzz_0[j] + fr * tg_yyyzzz_xxzz_1[j] + fl1_fxn * tg_yyyzzz_xzz_1[j];

                    tg_xyyyzzz_xyyy_0[j] = pb_x * tg_yyyzzz_xyyy_0[j] + fr * tg_yyyzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyy_1[j];

                    tg_xyyyzzz_xyyz_0[j] = pb_x * tg_yyyzzz_xyyz_0[j] + fr * tg_yyyzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yyz_1[j];

                    tg_xyyyzzz_xyzz_0[j] = pb_x * tg_yyyzzz_xyzz_0[j] + fr * tg_yyyzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yzz_1[j];

                    tg_xyyyzzz_xzzz_0[j] = pb_x * tg_yyyzzz_xzzz_0[j] + fr * tg_yyyzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_zzz_1[j];

                    tg_xyyyzzz_yyyy_0[j] = pb_x * tg_yyyzzz_yyyy_0[j] + fr * tg_yyyzzz_yyyy_1[j];

                    tg_xyyyzzz_yyyz_0[j] = pb_x * tg_yyyzzz_yyyz_0[j] + fr * tg_yyyzzz_yyyz_1[j];

                    tg_xyyyzzz_yyzz_0[j] = pb_x * tg_yyyzzz_yyzz_0[j] + fr * tg_yyyzzz_yyzz_1[j];

                    tg_xyyyzzz_yzzz_0[j] = pb_x * tg_yyyzzz_yzzz_0[j] + fr * tg_yyyzzz_yzzz_1[j];

                    tg_xyyyzzz_zzzz_0[j] = pb_x * tg_yyyzzz_zzzz_0[j] + fr * tg_yyyzzz_zzzz_1[j];

                    tg_xyyzzzz_xxxx_0[j] = pb_x * tg_yyzzzz_xxxx_0[j] + fr * tg_yyzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyzzzz_xxx_1[j];

                    tg_xyyzzzz_xxxy_0[j] = pb_x * tg_yyzzzz_xxxy_0[j] + fr * tg_yyzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxy_1[j];

                    tg_xyyzzzz_xxxz_0[j] = pb_x * tg_yyzzzz_xxxz_0[j] + fr * tg_yyzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xxz_1[j];

                    tg_xyyzzzz_xxyy_0[j] = pb_x * tg_yyzzzz_xxyy_0[j] + fr * tg_yyzzzz_xxyy_1[j] + fl1_fxn * tg_yyzzzz_xyy_1[j];

                    tg_xyyzzzz_xxyz_0[j] = pb_x * tg_yyzzzz_xxyz_0[j] + fr * tg_yyzzzz_xxyz_1[j] + fl1_fxn * tg_yyzzzz_xyz_1[j];

                    tg_xyyzzzz_xxzz_0[j] = pb_x * tg_yyzzzz_xxzz_0[j] + fr * tg_yyzzzz_xxzz_1[j] + fl1_fxn * tg_yyzzzz_xzz_1[j];

                    tg_xyyzzzz_xyyy_0[j] = pb_x * tg_yyzzzz_xyyy_0[j] + fr * tg_yyzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyy_1[j];

                    tg_xyyzzzz_xyyz_0[j] = pb_x * tg_yyzzzz_xyyz_0[j] + fr * tg_yyzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yyz_1[j];

                    tg_xyyzzzz_xyzz_0[j] = pb_x * tg_yyzzzz_xyzz_0[j] + fr * tg_yyzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yzz_1[j];

                    tg_xyyzzzz_xzzz_0[j] = pb_x * tg_yyzzzz_xzzz_0[j] + fr * tg_yyzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_zzz_1[j];

                    tg_xyyzzzz_yyyy_0[j] = pb_x * tg_yyzzzz_yyyy_0[j] + fr * tg_yyzzzz_yyyy_1[j];

                    tg_xyyzzzz_yyyz_0[j] = pb_x * tg_yyzzzz_yyyz_0[j] + fr * tg_yyzzzz_yyyz_1[j];

                    tg_xyyzzzz_yyzz_0[j] = pb_x * tg_yyzzzz_yyzz_0[j] + fr * tg_yyzzzz_yyzz_1[j];

                    tg_xyyzzzz_yzzz_0[j] = pb_x * tg_yyzzzz_yzzz_0[j] + fr * tg_yyzzzz_yzzz_1[j];

                    tg_xyyzzzz_zzzz_0[j] = pb_x * tg_yyzzzz_zzzz_0[j] + fr * tg_yyzzzz_zzzz_1[j];

                    tg_xyzzzzz_xxxx_0[j] = pb_x * tg_yzzzzz_xxxx_0[j] + fr * tg_yzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yzzzzz_xxx_1[j];

                    tg_xyzzzzz_xxxy_0[j] = pb_x * tg_yzzzzz_xxxy_0[j] + fr * tg_yzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxy_1[j];

                    tg_xyzzzzz_xxxz_0[j] = pb_x * tg_yzzzzz_xxxz_0[j] + fr * tg_yzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xxz_1[j];

                    tg_xyzzzzz_xxyy_0[j] = pb_x * tg_yzzzzz_xxyy_0[j] + fr * tg_yzzzzz_xxyy_1[j] + fl1_fxn * tg_yzzzzz_xyy_1[j];

                    tg_xyzzzzz_xxyz_0[j] = pb_x * tg_yzzzzz_xxyz_0[j] + fr * tg_yzzzzz_xxyz_1[j] + fl1_fxn * tg_yzzzzz_xyz_1[j];

                    tg_xyzzzzz_xxzz_0[j] = pb_x * tg_yzzzzz_xxzz_0[j] + fr * tg_yzzzzz_xxzz_1[j] + fl1_fxn * tg_yzzzzz_xzz_1[j];

                    tg_xyzzzzz_xyyy_0[j] = pb_x * tg_yzzzzz_xyyy_0[j] + fr * tg_yzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyy_1[j];

                    tg_xyzzzzz_xyyz_0[j] = pb_x * tg_yzzzzz_xyyz_0[j] + fr * tg_yzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yyz_1[j];

                    tg_xyzzzzz_xyzz_0[j] = pb_x * tg_yzzzzz_xyzz_0[j] + fr * tg_yzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yzz_1[j];

                    tg_xyzzzzz_xzzz_0[j] = pb_x * tg_yzzzzz_xzzz_0[j] + fr * tg_yzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_zzz_1[j];

                    tg_xyzzzzz_yyyy_0[j] = pb_x * tg_yzzzzz_yyyy_0[j] + fr * tg_yzzzzz_yyyy_1[j];

                    tg_xyzzzzz_yyyz_0[j] = pb_x * tg_yzzzzz_yyyz_0[j] + fr * tg_yzzzzz_yyyz_1[j];

                    tg_xyzzzzz_yyzz_0[j] = pb_x * tg_yzzzzz_yyzz_0[j] + fr * tg_yzzzzz_yyzz_1[j];

                    tg_xyzzzzz_yzzz_0[j] = pb_x * tg_yzzzzz_yzzz_0[j] + fr * tg_yzzzzz_yzzz_1[j];

                    tg_xyzzzzz_zzzz_0[j] = pb_x * tg_yzzzzz_zzzz_0[j] + fr * tg_yzzzzz_zzzz_1[j];

                    tg_xzzzzzz_xxxx_0[j] = pb_x * tg_zzzzzz_xxxx_0[j] + fr * tg_zzzzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_xxx_1[j];

                    tg_xzzzzzz_xxxy_0[j] = pb_x * tg_zzzzzz_xxxy_0[j] + fr * tg_zzzzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxy_1[j];

                    tg_xzzzzzz_xxxz_0[j] = pb_x * tg_zzzzzz_xxxz_0[j] + fr * tg_zzzzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xxz_1[j];

                    tg_xzzzzzz_xxyy_0[j] = pb_x * tg_zzzzzz_xxyy_0[j] + fr * tg_zzzzzz_xxyy_1[j] + fl1_fxn * tg_zzzzzz_xyy_1[j];

                    tg_xzzzzzz_xxyz_0[j] = pb_x * tg_zzzzzz_xxyz_0[j] + fr * tg_zzzzzz_xxyz_1[j] + fl1_fxn * tg_zzzzzz_xyz_1[j];

                    tg_xzzzzzz_xxzz_0[j] = pb_x * tg_zzzzzz_xxzz_0[j] + fr * tg_zzzzzz_xxzz_1[j] + fl1_fxn * tg_zzzzzz_xzz_1[j];

                    tg_xzzzzzz_xyyy_0[j] = pb_x * tg_zzzzzz_xyyy_0[j] + fr * tg_zzzzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyy_1[j];

                    tg_xzzzzzz_xyyz_0[j] = pb_x * tg_zzzzzz_xyyz_0[j] + fr * tg_zzzzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yyz_1[j];

                    tg_xzzzzzz_xyzz_0[j] = pb_x * tg_zzzzzz_xyzz_0[j] + fr * tg_zzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yzz_1[j];

                    tg_xzzzzzz_xzzz_0[j] = pb_x * tg_zzzzzz_xzzz_0[j] + fr * tg_zzzzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zzz_1[j];

                    tg_xzzzzzz_yyyy_0[j] = pb_x * tg_zzzzzz_yyyy_0[j] + fr * tg_zzzzzz_yyyy_1[j];

                    tg_xzzzzzz_yyyz_0[j] = pb_x * tg_zzzzzz_yyyz_0[j] + fr * tg_zzzzzz_yyyz_1[j];

                    tg_xzzzzzz_yyzz_0[j] = pb_x * tg_zzzzzz_yyzz_0[j] + fr * tg_zzzzzz_yyzz_1[j];

                    tg_xzzzzzz_yzzz_0[j] = pb_x * tg_zzzzzz_yzzz_0[j] + fr * tg_zzzzzz_yzzz_1[j];

                    tg_xzzzzzz_zzzz_0[j] = pb_x * tg_zzzzzz_zzzz_0[j] + fr * tg_zzzzzz_zzzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyyyy_xxxx_0[j] = pb_y * tg_yyyyyy_xxxx_0[j] + fr * tg_yyyyyy_xxxx_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxx_0[j] - tg_yyyyy_xxxx_1[j] * fl1_fza);

                    tg_yyyyyyy_xxxy_0[j] = pb_y * tg_yyyyyy_xxxy_0[j] + fr * tg_yyyyyy_xxxy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxy_0[j] - tg_yyyyy_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxx_1[j];

                    tg_yyyyyyy_xxxz_0[j] = pb_y * tg_yyyyyy_xxxz_0[j] + fr * tg_yyyyyy_xxxz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxxz_0[j] - tg_yyyyy_xxxz_1[j] * fl1_fza);

                    tg_yyyyyyy_xxyy_0[j] = pb_y * tg_yyyyyy_xxyy_0[j] + fr * tg_yyyyyy_xxyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyy_0[j] - tg_yyyyy_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xxy_1[j];

                    tg_yyyyyyy_xxyz_0[j] = pb_y * tg_yyyyyy_xxyz_0[j] + fr * tg_yyyyyy_xxyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxyz_0[j] - tg_yyyyy_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xxz_1[j];

                    tg_yyyyyyy_xxzz_0[j] = pb_y * tg_yyyyyy_xxzz_0[j] + fr * tg_yyyyyy_xxzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxzz_0[j] - tg_yyyyy_xxzz_1[j] * fl1_fza);

                    tg_yyyyyyy_xyyy_0[j] = pb_y * tg_yyyyyy_xyyy_0[j] + fr * tg_yyyyyy_xyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyy_0[j] - tg_yyyyy_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_xyy_1[j];

                    tg_yyyyyyy_xyyz_0[j] = pb_y * tg_yyyyyy_xyyz_0[j] + fr * tg_yyyyyy_xyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyyz_0[j] - tg_yyyyy_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xyz_1[j];

                    tg_yyyyyyy_xyzz_0[j] = pb_y * tg_yyyyyy_xyzz_0[j] + fr * tg_yyyyyy_xyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyzz_0[j] - tg_yyyyy_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xzz_1[j];

                    tg_yyyyyyy_xzzz_0[j] = pb_y * tg_yyyyyy_xzzz_0[j] + fr * tg_yyyyyy_xzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xzzz_0[j] - tg_yyyyy_xzzz_1[j] * fl1_fza);

                    tg_yyyyyyy_yyyy_0[j] = pb_y * tg_yyyyyy_yyyy_0[j] + fr * tg_yyyyyy_yyyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyy_0[j] - tg_yyyyy_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyy_yyy_1[j];

                    tg_yyyyyyy_yyyz_0[j] = pb_y * tg_yyyyyy_yyyz_0[j] + fr * tg_yyyyyy_yyyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyyz_0[j] - tg_yyyyy_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_yyz_1[j];

                    tg_yyyyyyy_yyzz_0[j] = pb_y * tg_yyyyyy_yyzz_0[j] + fr * tg_yyyyyy_yyzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyzz_0[j] - tg_yyyyy_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_yzz_1[j];

                    tg_yyyyyyy_yzzz_0[j] = pb_y * tg_yyyyyy_yzzz_0[j] + fr * tg_yyyyyy_yzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yzzz_0[j] - tg_yyyyy_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_zzz_1[j];

                    tg_yyyyyyy_zzzz_0[j] = pb_y * tg_yyyyyy_zzzz_0[j] + fr * tg_yyyyyy_zzzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_zzzz_0[j] - tg_yyyyy_zzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxx_0[j] = pb_y * tg_yyyyyz_xxxx_0[j] + fr * tg_yyyyyz_xxxx_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxx_0[j] - tg_yyyyz_xxxx_1[j] * fl1_fza);

                    tg_yyyyyyz_xxxy_0[j] = pb_y * tg_yyyyyz_xxxy_0[j] + fr * tg_yyyyyz_xxxy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxy_0[j] - tg_yyyyz_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxx_1[j];

                    tg_yyyyyyz_xxxz_0[j] = pb_y * tg_yyyyyz_xxxz_0[j] + fr * tg_yyyyyz_xxxz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxxz_0[j] - tg_yyyyz_xxxz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxyy_0[j] = pb_y * tg_yyyyyz_xxyy_0[j] + fr * tg_yyyyyz_xxyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyy_0[j] - tg_yyyyz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xxy_1[j];

                    tg_yyyyyyz_xxyz_0[j] = pb_y * tg_yyyyyz_xxyz_0[j] + fr * tg_yyyyyz_xxyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxyz_0[j] - tg_yyyyz_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xxz_1[j];

                    tg_yyyyyyz_xxzz_0[j] = pb_y * tg_yyyyyz_xxzz_0[j] + fr * tg_yyyyyz_xxzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxzz_0[j] - tg_yyyyz_xxzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xyyy_0[j] = pb_y * tg_yyyyyz_xyyy_0[j] + fr * tg_yyyyyz_xyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyy_0[j] - tg_yyyyz_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_xyy_1[j];

                    tg_yyyyyyz_xyyz_0[j] = pb_y * tg_yyyyyz_xyyz_0[j] + fr * tg_yyyyyz_xyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyyz_0[j] - tg_yyyyz_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xyz_1[j];

                    tg_yyyyyyz_xyzz_0[j] = pb_y * tg_yyyyyz_xyzz_0[j] + fr * tg_yyyyyz_xyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyzz_0[j] - tg_yyyyz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xzz_1[j];

                    tg_yyyyyyz_xzzz_0[j] = pb_y * tg_yyyyyz_xzzz_0[j] + fr * tg_yyyyyz_xzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xzzz_0[j] - tg_yyyyz_xzzz_1[j] * fl1_fza);

                    tg_yyyyyyz_yyyy_0[j] = pb_y * tg_yyyyyz_yyyy_0[j] + fr * tg_yyyyyz_yyyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyy_0[j] - tg_yyyyz_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyyz_yyy_1[j];

                    tg_yyyyyyz_yyyz_0[j] = pb_y * tg_yyyyyz_yyyz_0[j] + fr * tg_yyyyyz_yyyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyyz_0[j] - tg_yyyyz_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_yyz_1[j];

                    tg_yyyyyyz_yyzz_0[j] = pb_y * tg_yyyyyz_yyzz_0[j] + fr * tg_yyyyyz_yyzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyzz_0[j] - tg_yyyyz_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_yzz_1[j];

                    tg_yyyyyyz_yzzz_0[j] = pb_y * tg_yyyyyz_yzzz_0[j] + fr * tg_yyyyyz_yzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yzzz_0[j] - tg_yyyyz_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_zzz_1[j];

                    tg_yyyyyyz_zzzz_0[j] = pb_y * tg_yyyyyz_zzzz_0[j] + fr * tg_yyyyyz_zzzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_zzzz_0[j] - tg_yyyyz_zzzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSG_450_540(      CMemBlock2D<double>* primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (450,540)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yyyyzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 345); 

                auto tg_yyyyzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 346); 

                auto tg_yyyyzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 347); 

                auto tg_yyyyzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 348); 

                auto tg_yyyyzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 349); 

                auto tg_yyyyzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 350); 

                auto tg_yyyyzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 351); 

                auto tg_yyyyzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 352); 

                auto tg_yyyyzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 353); 

                auto tg_yyyyzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 354); 

                auto tg_yyyyzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 355); 

                auto tg_yyyyzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 356); 

                auto tg_yyyyzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 357); 

                auto tg_yyyyzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 358); 

                auto tg_yyyyzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 359); 

                auto tg_yyyzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 360); 

                auto tg_yyyzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 361); 

                auto tg_yyyzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 362); 

                auto tg_yyyzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 363); 

                auto tg_yyyzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 364); 

                auto tg_yyyzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 365); 

                auto tg_yyyzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 366); 

                auto tg_yyyzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 367); 

                auto tg_yyyzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 368); 

                auto tg_yyyzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 369); 

                auto tg_yyyzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 370); 

                auto tg_yyyzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 371); 

                auto tg_yyyzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 372); 

                auto tg_yyyzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 373); 

                auto tg_yyyzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 374); 

                auto tg_yyzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 375); 

                auto tg_yyzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 376); 

                auto tg_yyzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 377); 

                auto tg_yyzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 378); 

                auto tg_yyzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 379); 

                auto tg_yyzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 380); 

                auto tg_yyzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 381); 

                auto tg_yyzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 382); 

                auto tg_yyzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 383); 

                auto tg_yyzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 384); 

                auto tg_yyzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 385); 

                auto tg_yyzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 386); 

                auto tg_yyzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 387); 

                auto tg_yyzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 388); 

                auto tg_yyzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 389); 

                auto tg_yzzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 390); 

                auto tg_yzzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 391); 

                auto tg_yzzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 392); 

                auto tg_yzzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 393); 

                auto tg_yzzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 394); 

                auto tg_yzzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 395); 

                auto tg_yzzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 396); 

                auto tg_yzzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 397); 

                auto tg_yzzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 398); 

                auto tg_yzzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 399); 

                auto tg_yzzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 400); 

                auto tg_yzzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 401); 

                auto tg_yzzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 402); 

                auto tg_yzzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 403); 

                auto tg_yzzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 404); 

                auto tg_zzzzzz_xxxx_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 405); 

                auto tg_zzzzzz_xxxy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 406); 

                auto tg_zzzzzz_xxxz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 407); 

                auto tg_zzzzzz_xxyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 408); 

                auto tg_zzzzzz_xxyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 409); 

                auto tg_zzzzzz_xxzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 410); 

                auto tg_zzzzzz_xyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 411); 

                auto tg_zzzzzz_xyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 412); 

                auto tg_zzzzzz_xyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 413); 

                auto tg_zzzzzz_xzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 414); 

                auto tg_zzzzzz_yyyy_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 415); 

                auto tg_zzzzzz_yyyz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 416); 

                auto tg_zzzzzz_yyzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 417); 

                auto tg_zzzzzz_yzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 418); 

                auto tg_zzzzzz_zzzz_0 = primBuffer[pidx_g_6_4_m0].data(420 * idx + 419); 

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

                auto tg_yyyzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 255); 

                auto tg_yyyzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 256); 

                auto tg_yyyzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 257); 

                auto tg_yyyzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 258); 

                auto tg_yyyzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 259); 

                auto tg_yyyzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 260); 

                auto tg_yyyzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 261); 

                auto tg_yyyzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 262); 

                auto tg_yyyzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 263); 

                auto tg_yyyzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 264); 

                auto tg_yyyzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 265); 

                auto tg_yyyzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 266); 

                auto tg_yyyzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 267); 

                auto tg_yyyzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 268); 

                auto tg_yyyzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 269); 

                auto tg_yyzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 270); 

                auto tg_yyzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 271); 

                auto tg_yyzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 272); 

                auto tg_yyzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 273); 

                auto tg_yyzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 274); 

                auto tg_yyzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 275); 

                auto tg_yyzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 276); 

                auto tg_yyzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 277); 

                auto tg_yyzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 278); 

                auto tg_yyzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 279); 

                auto tg_yyzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 280); 

                auto tg_yyzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 281); 

                auto tg_yyzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 282); 

                auto tg_yyzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 283); 

                auto tg_yyzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 284); 

                auto tg_yzzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 285); 

                auto tg_yzzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 286); 

                auto tg_yzzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 287); 

                auto tg_yzzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 288); 

                auto tg_yzzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 289); 

                auto tg_yzzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 290); 

                auto tg_yzzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 291); 

                auto tg_yzzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 292); 

                auto tg_yzzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 293); 

                auto tg_yzzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 294); 

                auto tg_yzzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 295); 

                auto tg_yzzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 296); 

                auto tg_yzzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 297); 

                auto tg_yzzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 298); 

                auto tg_yzzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 299); 

                auto tg_zzzzz_xxxx_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 300); 

                auto tg_zzzzz_xxxy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 301); 

                auto tg_zzzzz_xxxz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 302); 

                auto tg_zzzzz_xxyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 303); 

                auto tg_zzzzz_xxyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 304); 

                auto tg_zzzzz_xxzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 305); 

                auto tg_zzzzz_xyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 306); 

                auto tg_zzzzz_xyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 307); 

                auto tg_zzzzz_xyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 308); 

                auto tg_zzzzz_xzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 309); 

                auto tg_zzzzz_yyyy_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 310); 

                auto tg_zzzzz_yyyz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 311); 

                auto tg_zzzzz_yyzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 312); 

                auto tg_zzzzz_yzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 313); 

                auto tg_zzzzz_zzzz_0 = primBuffer[pidx_g_5_4_m0].data(315 * idx + 314); 

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

                auto tg_yyyyzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 230); 

                auto tg_yyyyzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 231); 

                auto tg_yyyyzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 232); 

                auto tg_yyyyzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 233); 

                auto tg_yyyyzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 234); 

                auto tg_yyyyzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 235); 

                auto tg_yyyyzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 236); 

                auto tg_yyyyzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 237); 

                auto tg_yyyyzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 238); 

                auto tg_yyyyzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 239); 

                auto tg_yyyzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 240); 

                auto tg_yyyzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 241); 

                auto tg_yyyzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 242); 

                auto tg_yyyzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 243); 

                auto tg_yyyzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 244); 

                auto tg_yyyzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 245); 

                auto tg_yyyzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 246); 

                auto tg_yyyzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 247); 

                auto tg_yyyzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 248); 

                auto tg_yyyzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 249); 

                auto tg_yyzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 250); 

                auto tg_yyzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 251); 

                auto tg_yyzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 252); 

                auto tg_yyzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 253); 

                auto tg_yyzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 254); 

                auto tg_yyzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 255); 

                auto tg_yyzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 256); 

                auto tg_yyzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 257); 

                auto tg_yyzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 258); 

                auto tg_yyzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 259); 

                auto tg_yzzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 260); 

                auto tg_yzzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 261); 

                auto tg_yzzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 262); 

                auto tg_yzzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 263); 

                auto tg_yzzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 264); 

                auto tg_yzzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 265); 

                auto tg_yzzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 266); 

                auto tg_yzzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 267); 

                auto tg_yzzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 268); 

                auto tg_yzzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 269); 

                auto tg_zzzzzz_xxx_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 270); 

                auto tg_zzzzzz_xxy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 271); 

                auto tg_zzzzzz_xxz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 272); 

                auto tg_zzzzzz_xyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 273); 

                auto tg_zzzzzz_xyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 274); 

                auto tg_zzzzzz_xzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 275); 

                auto tg_zzzzzz_yyy_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 276); 

                auto tg_zzzzzz_yyz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 277); 

                auto tg_zzzzzz_yzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 278); 

                auto tg_zzzzzz_zzz_1 = primBuffer[pidx_g_6_3_m1].data(280 * idx + 279); 

                // set up pointers to integrals

                auto tg_yyyyyzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 450); 

                auto tg_yyyyyzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 451); 

                auto tg_yyyyyzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 452); 

                auto tg_yyyyyzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 453); 

                auto tg_yyyyyzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 454); 

                auto tg_yyyyyzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 455); 

                auto tg_yyyyyzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 456); 

                auto tg_yyyyyzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 457); 

                auto tg_yyyyyzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 458); 

                auto tg_yyyyyzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 459); 

                auto tg_yyyyyzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 460); 

                auto tg_yyyyyzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 461); 

                auto tg_yyyyyzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 462); 

                auto tg_yyyyyzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 463); 

                auto tg_yyyyyzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 464); 

                auto tg_yyyyzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 465); 

                auto tg_yyyyzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 466); 

                auto tg_yyyyzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 467); 

                auto tg_yyyyzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 468); 

                auto tg_yyyyzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 469); 

                auto tg_yyyyzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 470); 

                auto tg_yyyyzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 471); 

                auto tg_yyyyzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 472); 

                auto tg_yyyyzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 473); 

                auto tg_yyyyzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 474); 

                auto tg_yyyyzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 475); 

                auto tg_yyyyzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 476); 

                auto tg_yyyyzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 477); 

                auto tg_yyyyzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 478); 

                auto tg_yyyyzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 479); 

                auto tg_yyyzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 480); 

                auto tg_yyyzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 481); 

                auto tg_yyyzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 482); 

                auto tg_yyyzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 483); 

                auto tg_yyyzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 484); 

                auto tg_yyyzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 485); 

                auto tg_yyyzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 486); 

                auto tg_yyyzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 487); 

                auto tg_yyyzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 488); 

                auto tg_yyyzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 489); 

                auto tg_yyyzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 490); 

                auto tg_yyyzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 491); 

                auto tg_yyyzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 492); 

                auto tg_yyyzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 493); 

                auto tg_yyyzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 494); 

                auto tg_yyzzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 495); 

                auto tg_yyzzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 496); 

                auto tg_yyzzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 497); 

                auto tg_yyzzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 498); 

                auto tg_yyzzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 499); 

                auto tg_yyzzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 500); 

                auto tg_yyzzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 501); 

                auto tg_yyzzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 502); 

                auto tg_yyzzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 503); 

                auto tg_yyzzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 504); 

                auto tg_yyzzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 505); 

                auto tg_yyzzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 506); 

                auto tg_yyzzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 507); 

                auto tg_yyzzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 508); 

                auto tg_yyzzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 509); 

                auto tg_yzzzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 510); 

                auto tg_yzzzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 511); 

                auto tg_yzzzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 512); 

                auto tg_yzzzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 513); 

                auto tg_yzzzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 514); 

                auto tg_yzzzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 515); 

                auto tg_yzzzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 516); 

                auto tg_yzzzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 517); 

                auto tg_yzzzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 518); 

                auto tg_yzzzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 519); 

                auto tg_yzzzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 520); 

                auto tg_yzzzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 521); 

                auto tg_yzzzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 522); 

                auto tg_yzzzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 523); 

                auto tg_yzzzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 524); 

                auto tg_zzzzzzz_xxxx_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 525); 

                auto tg_zzzzzzz_xxxy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 526); 

                auto tg_zzzzzzz_xxxz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 527); 

                auto tg_zzzzzzz_xxyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 528); 

                auto tg_zzzzzzz_xxyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 529); 

                auto tg_zzzzzzz_xxzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 530); 

                auto tg_zzzzzzz_xyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 531); 

                auto tg_zzzzzzz_xyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 532); 

                auto tg_zzzzzzz_xyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 533); 

                auto tg_zzzzzzz_xzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 534); 

                auto tg_zzzzzzz_yyyy_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 535); 

                auto tg_zzzzzzz_yyyz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 536); 

                auto tg_zzzzzzz_yyzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 537); 

                auto tg_zzzzzzz_yzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 538); 

                auto tg_zzzzzzz_zzzz_0 = primBuffer[pidx_g_7_4_m0].data(540 * idx + 539); 

                // Batch of Integrals (450,540)

                #pragma omp simd aligned(fxn, fza, tg_yyyyyzz_xxxx_0, tg_yyyyyzz_xxxy_0, tg_yyyyyzz_xxxz_0, \
                                         tg_yyyyyzz_xxyy_0, tg_yyyyyzz_xxyz_0, tg_yyyyyzz_xxzz_0, tg_yyyyyzz_xyyy_0, \
                                         tg_yyyyyzz_xyyz_0, tg_yyyyyzz_xyzz_0, tg_yyyyyzz_xzzz_0, tg_yyyyyzz_yyyy_0, \
                                         tg_yyyyyzz_yyyz_0, tg_yyyyyzz_yyzz_0, tg_yyyyyzz_yzzz_0, tg_yyyyyzz_zzzz_0, \
                                         tg_yyyyzz_xxx_1, tg_yyyyzz_xxxx_0, tg_yyyyzz_xxxx_1, tg_yyyyzz_xxxy_0, \
                                         tg_yyyyzz_xxxy_1, tg_yyyyzz_xxxz_0, tg_yyyyzz_xxxz_1, tg_yyyyzz_xxy_1, \
                                         tg_yyyyzz_xxyy_0, tg_yyyyzz_xxyy_1, tg_yyyyzz_xxyz_0, tg_yyyyzz_xxyz_1, \
                                         tg_yyyyzz_xxz_1, tg_yyyyzz_xxzz_0, tg_yyyyzz_xxzz_1, tg_yyyyzz_xyy_1, \
                                         tg_yyyyzz_xyyy_0, tg_yyyyzz_xyyy_1, tg_yyyyzz_xyyz_0, tg_yyyyzz_xyyz_1, \
                                         tg_yyyyzz_xyz_1, tg_yyyyzz_xyzz_0, tg_yyyyzz_xyzz_1, tg_yyyyzz_xzz_1, \
                                         tg_yyyyzz_xzzz_0, tg_yyyyzz_xzzz_1, tg_yyyyzz_yyy_1, tg_yyyyzz_yyyy_0, \
                                         tg_yyyyzz_yyyy_1, tg_yyyyzz_yyyz_0, tg_yyyyzz_yyyz_1, tg_yyyyzz_yyz_1, \
                                         tg_yyyyzz_yyzz_0, tg_yyyyzz_yyzz_1, tg_yyyyzz_yzz_1, tg_yyyyzz_yzzz_0, \
                                         tg_yyyyzz_yzzz_1, tg_yyyyzz_zzz_1, tg_yyyyzz_zzzz_0, tg_yyyyzz_zzzz_1, \
                                         tg_yyyyzzz_xxxx_0, tg_yyyyzzz_xxxy_0, tg_yyyyzzz_xxxz_0, tg_yyyyzzz_xxyy_0, \
                                         tg_yyyyzzz_xxyz_0, tg_yyyyzzz_xxzz_0, tg_yyyyzzz_xyyy_0, tg_yyyyzzz_xyyz_0, \
                                         tg_yyyyzzz_xyzz_0, tg_yyyyzzz_xzzz_0, tg_yyyyzzz_yyyy_0, tg_yyyyzzz_yyyz_0, \
                                         tg_yyyyzzz_yyzz_0, tg_yyyyzzz_yzzz_0, tg_yyyyzzz_zzzz_0, tg_yyyzz_xxxx_0, \
                                         tg_yyyzz_xxxx_1, tg_yyyzz_xxxy_0, tg_yyyzz_xxxy_1, tg_yyyzz_xxxz_0, tg_yyyzz_xxxz_1, \
                                         tg_yyyzz_xxyy_0, tg_yyyzz_xxyy_1, tg_yyyzz_xxyz_0, tg_yyyzz_xxyz_1, tg_yyyzz_xxzz_0, \
                                         tg_yyyzz_xxzz_1, tg_yyyzz_xyyy_0, tg_yyyzz_xyyy_1, tg_yyyzz_xyyz_0, tg_yyyzz_xyyz_1, \
                                         tg_yyyzz_xyzz_0, tg_yyyzz_xyzz_1, tg_yyyzz_xzzz_0, tg_yyyzz_xzzz_1, tg_yyyzz_yyyy_0, \
                                         tg_yyyzz_yyyy_1, tg_yyyzz_yyyz_0, tg_yyyzz_yyyz_1, tg_yyyzz_yyzz_0, tg_yyyzz_yyzz_1, \
                                         tg_yyyzz_yzzz_0, tg_yyyzz_yzzz_1, tg_yyyzz_zzzz_0, tg_yyyzz_zzzz_1, tg_yyyzzz_xxx_1, \
                                         tg_yyyzzz_xxxx_0, tg_yyyzzz_xxxx_1, tg_yyyzzz_xxxy_0, tg_yyyzzz_xxxy_1, \
                                         tg_yyyzzz_xxxz_0, tg_yyyzzz_xxxz_1, tg_yyyzzz_xxy_1, tg_yyyzzz_xxyy_0, \
                                         tg_yyyzzz_xxyy_1, tg_yyyzzz_xxyz_0, tg_yyyzzz_xxyz_1, tg_yyyzzz_xxz_1, \
                                         tg_yyyzzz_xxzz_0, tg_yyyzzz_xxzz_1, tg_yyyzzz_xyy_1, tg_yyyzzz_xyyy_0, \
                                         tg_yyyzzz_xyyy_1, tg_yyyzzz_xyyz_0, tg_yyyzzz_xyyz_1, tg_yyyzzz_xyz_1, \
                                         tg_yyyzzz_xyzz_0, tg_yyyzzz_xyzz_1, tg_yyyzzz_xzz_1, tg_yyyzzz_xzzz_0, \
                                         tg_yyyzzz_xzzz_1, tg_yyyzzz_yyy_1, tg_yyyzzz_yyyy_0, tg_yyyzzz_yyyy_1, \
                                         tg_yyyzzz_yyyz_0, tg_yyyzzz_yyyz_1, tg_yyyzzz_yyz_1, tg_yyyzzz_yyzz_0, \
                                         tg_yyyzzz_yyzz_1, tg_yyyzzz_yzz_1, tg_yyyzzz_yzzz_0, tg_yyyzzz_yzzz_1, \
                                         tg_yyyzzz_zzz_1, tg_yyyzzz_zzzz_0, tg_yyyzzz_zzzz_1, tg_yyyzzzz_xxxx_0, \
                                         tg_yyyzzzz_xxxy_0, tg_yyyzzzz_xxxz_0, tg_yyyzzzz_xxyy_0, tg_yyyzzzz_xxyz_0, \
                                         tg_yyyzzzz_xxzz_0, tg_yyyzzzz_xyyy_0, tg_yyyzzzz_xyyz_0, tg_yyyzzzz_xyzz_0, \
                                         tg_yyyzzzz_xzzz_0, tg_yyyzzzz_yyyy_0, tg_yyyzzzz_yyyz_0, tg_yyyzzzz_yyzz_0, \
                                         tg_yyyzzzz_yzzz_0, tg_yyyzzzz_zzzz_0, tg_yyzzz_xxxx_0, tg_yyzzz_xxxx_1, \
                                         tg_yyzzz_xxxy_0, tg_yyzzz_xxxy_1, tg_yyzzz_xxxz_0, tg_yyzzz_xxxz_1, tg_yyzzz_xxyy_0, \
                                         tg_yyzzz_xxyy_1, tg_yyzzz_xxyz_0, tg_yyzzz_xxyz_1, tg_yyzzz_xxzz_0, tg_yyzzz_xxzz_1, \
                                         tg_yyzzz_xyyy_0, tg_yyzzz_xyyy_1, tg_yyzzz_xyyz_0, tg_yyzzz_xyyz_1, tg_yyzzz_xyzz_0, \
                                         tg_yyzzz_xyzz_1, tg_yyzzz_xzzz_0, tg_yyzzz_xzzz_1, tg_yyzzz_yyyy_0, tg_yyzzz_yyyy_1, \
                                         tg_yyzzz_yyyz_0, tg_yyzzz_yyyz_1, tg_yyzzz_yyzz_0, tg_yyzzz_yyzz_1, tg_yyzzz_yzzz_0, \
                                         tg_yyzzz_yzzz_1, tg_yyzzz_zzzz_0, tg_yyzzz_zzzz_1, tg_yyzzzz_xxx_1, \
                                         tg_yyzzzz_xxxx_0, tg_yyzzzz_xxxx_1, tg_yyzzzz_xxxy_0, tg_yyzzzz_xxxy_1, \
                                         tg_yyzzzz_xxxz_0, tg_yyzzzz_xxxz_1, tg_yyzzzz_xxy_1, tg_yyzzzz_xxyy_0, \
                                         tg_yyzzzz_xxyy_1, tg_yyzzzz_xxyz_0, tg_yyzzzz_xxyz_1, tg_yyzzzz_xxz_1, \
                                         tg_yyzzzz_xxzz_0, tg_yyzzzz_xxzz_1, tg_yyzzzz_xyy_1, tg_yyzzzz_xyyy_0, \
                                         tg_yyzzzz_xyyy_1, tg_yyzzzz_xyyz_0, tg_yyzzzz_xyyz_1, tg_yyzzzz_xyz_1, \
                                         tg_yyzzzz_xyzz_0, tg_yyzzzz_xyzz_1, tg_yyzzzz_xzz_1, tg_yyzzzz_xzzz_0, \
                                         tg_yyzzzz_xzzz_1, tg_yyzzzz_yyy_1, tg_yyzzzz_yyyy_0, tg_yyzzzz_yyyy_1, \
                                         tg_yyzzzz_yyyz_0, tg_yyzzzz_yyyz_1, tg_yyzzzz_yyz_1, tg_yyzzzz_yyzz_0, \
                                         tg_yyzzzz_yyzz_1, tg_yyzzzz_yzz_1, tg_yyzzzz_yzzz_0, tg_yyzzzz_yzzz_1, \
                                         tg_yyzzzz_zzz_1, tg_yyzzzz_zzzz_0, tg_yyzzzz_zzzz_1, tg_yyzzzzz_xxxx_0, \
                                         tg_yyzzzzz_xxxy_0, tg_yyzzzzz_xxxz_0, tg_yyzzzzz_xxyy_0, tg_yyzzzzz_xxyz_0, \
                                         tg_yyzzzzz_xxzz_0, tg_yyzzzzz_xyyy_0, tg_yyzzzzz_xyyz_0, tg_yyzzzzz_xyzz_0, \
                                         tg_yyzzzzz_xzzz_0, tg_yyzzzzz_yyyy_0, tg_yyzzzzz_yyyz_0, tg_yyzzzzz_yyzz_0, \
                                         tg_yyzzzzz_yzzz_0, tg_yyzzzzz_zzzz_0, tg_yzzzz_xxxx_0, tg_yzzzz_xxxx_1, \
                                         tg_yzzzz_xxxy_0, tg_yzzzz_xxxy_1, tg_yzzzz_xxxz_0, tg_yzzzz_xxxz_1, tg_yzzzz_xxyy_0, \
                                         tg_yzzzz_xxyy_1, tg_yzzzz_xxyz_0, tg_yzzzz_xxyz_1, tg_yzzzz_xxzz_0, tg_yzzzz_xxzz_1, \
                                         tg_yzzzz_xyyy_0, tg_yzzzz_xyyy_1, tg_yzzzz_xyyz_0, tg_yzzzz_xyyz_1, tg_yzzzz_xyzz_0, \
                                         tg_yzzzz_xyzz_1, tg_yzzzz_xzzz_0, tg_yzzzz_xzzz_1, tg_yzzzz_yyyy_0, tg_yzzzz_yyyy_1, \
                                         tg_yzzzz_yyyz_0, tg_yzzzz_yyyz_1, tg_yzzzz_yyzz_0, tg_yzzzz_yyzz_1, tg_yzzzz_yzzz_0, \
                                         tg_yzzzz_yzzz_1, tg_yzzzz_zzzz_0, tg_yzzzz_zzzz_1, tg_yzzzzz_xxx_1, \
                                         tg_yzzzzz_xxxx_0, tg_yzzzzz_xxxx_1, tg_yzzzzz_xxxy_0, tg_yzzzzz_xxxy_1, \
                                         tg_yzzzzz_xxxz_0, tg_yzzzzz_xxxz_1, tg_yzzzzz_xxy_1, tg_yzzzzz_xxyy_0, \
                                         tg_yzzzzz_xxyy_1, tg_yzzzzz_xxyz_0, tg_yzzzzz_xxyz_1, tg_yzzzzz_xxz_1, \
                                         tg_yzzzzz_xxzz_0, tg_yzzzzz_xxzz_1, tg_yzzzzz_xyy_1, tg_yzzzzz_xyyy_0, \
                                         tg_yzzzzz_xyyy_1, tg_yzzzzz_xyyz_0, tg_yzzzzz_xyyz_1, tg_yzzzzz_xyz_1, \
                                         tg_yzzzzz_xyzz_0, tg_yzzzzz_xyzz_1, tg_yzzzzz_xzz_1, tg_yzzzzz_xzzz_0, \
                                         tg_yzzzzz_xzzz_1, tg_yzzzzz_yyy_1, tg_yzzzzz_yyyy_0, tg_yzzzzz_yyyy_1, \
                                         tg_yzzzzz_yyyz_0, tg_yzzzzz_yyyz_1, tg_yzzzzz_yyz_1, tg_yzzzzz_yyzz_0, \
                                         tg_yzzzzz_yyzz_1, tg_yzzzzz_yzz_1, tg_yzzzzz_yzzz_0, tg_yzzzzz_yzzz_1, \
                                         tg_yzzzzz_zzz_1, tg_yzzzzz_zzzz_0, tg_yzzzzz_zzzz_1, tg_yzzzzzz_xxxx_0, \
                                         tg_yzzzzzz_xxxy_0, tg_yzzzzzz_xxxz_0, tg_yzzzzzz_xxyy_0, tg_yzzzzzz_xxyz_0, \
                                         tg_yzzzzzz_xxzz_0, tg_yzzzzzz_xyyy_0, tg_yzzzzzz_xyyz_0, tg_yzzzzzz_xyzz_0, \
                                         tg_yzzzzzz_xzzz_0, tg_yzzzzzz_yyyy_0, tg_yzzzzzz_yyyz_0, tg_yzzzzzz_yyzz_0, \
                                         tg_yzzzzzz_yzzz_0, tg_yzzzzzz_zzzz_0, tg_zzzzz_xxxx_0, tg_zzzzz_xxxx_1, \
                                         tg_zzzzz_xxxy_0, tg_zzzzz_xxxy_1, tg_zzzzz_xxxz_0, tg_zzzzz_xxxz_1, tg_zzzzz_xxyy_0, \
                                         tg_zzzzz_xxyy_1, tg_zzzzz_xxyz_0, tg_zzzzz_xxyz_1, tg_zzzzz_xxzz_0, tg_zzzzz_xxzz_1, \
                                         tg_zzzzz_xyyy_0, tg_zzzzz_xyyy_1, tg_zzzzz_xyyz_0, tg_zzzzz_xyyz_1, tg_zzzzz_xyzz_0, \
                                         tg_zzzzz_xyzz_1, tg_zzzzz_xzzz_0, tg_zzzzz_xzzz_1, tg_zzzzz_yyyy_0, tg_zzzzz_yyyy_1, \
                                         tg_zzzzz_yyyz_0, tg_zzzzz_yyyz_1, tg_zzzzz_yyzz_0, tg_zzzzz_yyzz_1, tg_zzzzz_yzzz_0, \
                                         tg_zzzzz_yzzz_1, tg_zzzzz_zzzz_0, tg_zzzzz_zzzz_1, tg_zzzzzz_xxx_1, \
                                         tg_zzzzzz_xxxx_0, tg_zzzzzz_xxxx_1, tg_zzzzzz_xxxy_0, tg_zzzzzz_xxxy_1, \
                                         tg_zzzzzz_xxxz_0, tg_zzzzzz_xxxz_1, tg_zzzzzz_xxy_1, tg_zzzzzz_xxyy_0, \
                                         tg_zzzzzz_xxyy_1, tg_zzzzzz_xxyz_0, tg_zzzzzz_xxyz_1, tg_zzzzzz_xxz_1, \
                                         tg_zzzzzz_xxzz_0, tg_zzzzzz_xxzz_1, tg_zzzzzz_xyy_1, tg_zzzzzz_xyyy_0, \
                                         tg_zzzzzz_xyyy_1, tg_zzzzzz_xyyz_0, tg_zzzzzz_xyyz_1, tg_zzzzzz_xyz_1, \
                                         tg_zzzzzz_xyzz_0, tg_zzzzzz_xyzz_1, tg_zzzzzz_xzz_1, tg_zzzzzz_xzzz_0, \
                                         tg_zzzzzz_xzzz_1, tg_zzzzzz_yyy_1, tg_zzzzzz_yyyy_0, tg_zzzzzz_yyyy_1, \
                                         tg_zzzzzz_yyyz_0, tg_zzzzzz_yyyz_1, tg_zzzzzz_yyz_1, tg_zzzzzz_yyzz_0, \
                                         tg_zzzzzz_yyzz_1, tg_zzzzzz_yzz_1, tg_zzzzzz_yzzz_0, tg_zzzzzz_yzzz_1, \
                                         tg_zzzzzz_zzz_1, tg_zzzzzz_zzzz_0, tg_zzzzzz_zzzz_1, tg_zzzzzzz_xxxx_0, \
                                         tg_zzzzzzz_xxxy_0, tg_zzzzzzz_xxxz_0, tg_zzzzzzz_xxyy_0, tg_zzzzzzz_xxyz_0, \
                                         tg_zzzzzzz_xxzz_0, tg_zzzzzzz_xyyy_0, tg_zzzzzzz_xyyz_0, tg_zzzzzzz_xyzz_0, \
                                         tg_zzzzzzz_xzzz_0, tg_zzzzzzz_yyyy_0, tg_zzzzzzz_yyyz_0, tg_zzzzzzz_yyzz_0, \
                                         tg_zzzzzzz_yzzz_0, tg_zzzzzzz_zzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyyzz_xxxx_0[j] = pb_y * tg_yyyyzz_xxxx_0[j] + fr * tg_yyyyzz_xxxx_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxx_0[j] - tg_yyyzz_xxxx_1[j] * fl1_fza);

                    tg_yyyyyzz_xxxy_0[j] = pb_y * tg_yyyyzz_xxxy_0[j] + fr * tg_yyyyzz_xxxy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxy_0[j] - tg_yyyzz_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxx_1[j];

                    tg_yyyyyzz_xxxz_0[j] = pb_y * tg_yyyyzz_xxxz_0[j] + fr * tg_yyyyzz_xxxz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxxz_0[j] - tg_yyyzz_xxxz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxyy_0[j] = pb_y * tg_yyyyzz_xxyy_0[j] + fr * tg_yyyyzz_xxyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyy_0[j] - tg_yyyzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xxy_1[j];

                    tg_yyyyyzz_xxyz_0[j] = pb_y * tg_yyyyzz_xxyz_0[j] + fr * tg_yyyyzz_xxyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxyz_0[j] - tg_yyyzz_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xxz_1[j];

                    tg_yyyyyzz_xxzz_0[j] = pb_y * tg_yyyyzz_xxzz_0[j] + fr * tg_yyyyzz_xxzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxzz_0[j] - tg_yyyzz_xxzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xyyy_0[j] = pb_y * tg_yyyyzz_xyyy_0[j] + fr * tg_yyyyzz_xyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyy_0[j] - tg_yyyzz_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_xyy_1[j];

                    tg_yyyyyzz_xyyz_0[j] = pb_y * tg_yyyyzz_xyyz_0[j] + fr * tg_yyyyzz_xyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyyz_0[j] - tg_yyyzz_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xyz_1[j];

                    tg_yyyyyzz_xyzz_0[j] = pb_y * tg_yyyyzz_xyzz_0[j] + fr * tg_yyyyzz_xyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyzz_0[j] - tg_yyyzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xzz_1[j];

                    tg_yyyyyzz_xzzz_0[j] = pb_y * tg_yyyyzz_xzzz_0[j] + fr * tg_yyyyzz_xzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xzzz_0[j] - tg_yyyzz_xzzz_1[j] * fl1_fza);

                    tg_yyyyyzz_yyyy_0[j] = pb_y * tg_yyyyzz_yyyy_0[j] + fr * tg_yyyyzz_yyyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyy_0[j] - tg_yyyzz_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyyzz_yyy_1[j];

                    tg_yyyyyzz_yyyz_0[j] = pb_y * tg_yyyyzz_yyyz_0[j] + fr * tg_yyyyzz_yyyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyyz_0[j] - tg_yyyzz_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_yyz_1[j];

                    tg_yyyyyzz_yyzz_0[j] = pb_y * tg_yyyyzz_yyzz_0[j] + fr * tg_yyyyzz_yyzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyzz_0[j] - tg_yyyzz_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_yzz_1[j];

                    tg_yyyyyzz_yzzz_0[j] = pb_y * tg_yyyyzz_yzzz_0[j] + fr * tg_yyyyzz_yzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yzzz_0[j] - tg_yyyzz_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_zzz_1[j];

                    tg_yyyyyzz_zzzz_0[j] = pb_y * tg_yyyyzz_zzzz_0[j] + fr * tg_yyyyzz_zzzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_zzzz_0[j] - tg_yyyzz_zzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxx_0[j] = pb_y * tg_yyyzzz_xxxx_0[j] + fr * tg_yyyzzz_xxxx_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxx_0[j] - tg_yyzzz_xxxx_1[j] * fl1_fza);

                    tg_yyyyzzz_xxxy_0[j] = pb_y * tg_yyyzzz_xxxy_0[j] + fr * tg_yyyzzz_xxxy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxy_0[j] - tg_yyzzz_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxx_1[j];

                    tg_yyyyzzz_xxxz_0[j] = pb_y * tg_yyyzzz_xxxz_0[j] + fr * tg_yyyzzz_xxxz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxxz_0[j] - tg_yyzzz_xxxz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxyy_0[j] = pb_y * tg_yyyzzz_xxyy_0[j] + fr * tg_yyyzzz_xxyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyy_0[j] - tg_yyzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xxy_1[j];

                    tg_yyyyzzz_xxyz_0[j] = pb_y * tg_yyyzzz_xxyz_0[j] + fr * tg_yyyzzz_xxyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxyz_0[j] - tg_yyzzz_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xxz_1[j];

                    tg_yyyyzzz_xxzz_0[j] = pb_y * tg_yyyzzz_xxzz_0[j] + fr * tg_yyyzzz_xxzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxzz_0[j] - tg_yyzzz_xxzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xyyy_0[j] = pb_y * tg_yyyzzz_xyyy_0[j] + fr * tg_yyyzzz_xyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyy_0[j] - tg_yyzzz_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_xyy_1[j];

                    tg_yyyyzzz_xyyz_0[j] = pb_y * tg_yyyzzz_xyyz_0[j] + fr * tg_yyyzzz_xyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyyz_0[j] - tg_yyzzz_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xyz_1[j];

                    tg_yyyyzzz_xyzz_0[j] = pb_y * tg_yyyzzz_xyzz_0[j] + fr * tg_yyyzzz_xyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyzz_0[j] - tg_yyzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xzz_1[j];

                    tg_yyyyzzz_xzzz_0[j] = pb_y * tg_yyyzzz_xzzz_0[j] + fr * tg_yyyzzz_xzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xzzz_0[j] - tg_yyzzz_xzzz_1[j] * fl1_fza);

                    tg_yyyyzzz_yyyy_0[j] = pb_y * tg_yyyzzz_yyyy_0[j] + fr * tg_yyyzzz_yyyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyy_0[j] - tg_yyzzz_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyyzzz_yyy_1[j];

                    tg_yyyyzzz_yyyz_0[j] = pb_y * tg_yyyzzz_yyyz_0[j] + fr * tg_yyyzzz_yyyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyyz_0[j] - tg_yyzzz_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_yyz_1[j];

                    tg_yyyyzzz_yyzz_0[j] = pb_y * tg_yyyzzz_yyzz_0[j] + fr * tg_yyyzzz_yyzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyzz_0[j] - tg_yyzzz_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_yzz_1[j];

                    tg_yyyyzzz_yzzz_0[j] = pb_y * tg_yyyzzz_yzzz_0[j] + fr * tg_yyyzzz_yzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yzzz_0[j] - tg_yyzzz_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_zzz_1[j];

                    tg_yyyyzzz_zzzz_0[j] = pb_y * tg_yyyzzz_zzzz_0[j] + fr * tg_yyyzzz_zzzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_zzzz_0[j] - tg_yyzzz_zzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxx_0[j] = pb_y * tg_yyzzzz_xxxx_0[j] + fr * tg_yyzzzz_xxxx_1[j] + fl1_fx * (tg_yzzzz_xxxx_0[j] - tg_yzzzz_xxxx_1[j] * fl1_fza);

                    tg_yyyzzzz_xxxy_0[j] = pb_y * tg_yyzzzz_xxxy_0[j] + fr * tg_yyzzzz_xxxy_1[j] + fl1_fx * (tg_yzzzz_xxxy_0[j] - tg_yzzzz_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxx_1[j];

                    tg_yyyzzzz_xxxz_0[j] = pb_y * tg_yyzzzz_xxxz_0[j] + fr * tg_yyzzzz_xxxz_1[j] + fl1_fx * (tg_yzzzz_xxxz_0[j] - tg_yzzzz_xxxz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxyy_0[j] = pb_y * tg_yyzzzz_xxyy_0[j] + fr * tg_yyzzzz_xxyy_1[j] + fl1_fx * (tg_yzzzz_xxyy_0[j] - tg_yzzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xxy_1[j];

                    tg_yyyzzzz_xxyz_0[j] = pb_y * tg_yyzzzz_xxyz_0[j] + fr * tg_yyzzzz_xxyz_1[j] + fl1_fx * (tg_yzzzz_xxyz_0[j] - tg_yzzzz_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xxz_1[j];

                    tg_yyyzzzz_xxzz_0[j] = pb_y * tg_yyzzzz_xxzz_0[j] + fr * tg_yyzzzz_xxzz_1[j] + fl1_fx * (tg_yzzzz_xxzz_0[j] - tg_yzzzz_xxzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xyyy_0[j] = pb_y * tg_yyzzzz_xyyy_0[j] + fr * tg_yyzzzz_xyyy_1[j] + fl1_fx * (tg_yzzzz_xyyy_0[j] - tg_yzzzz_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_xyy_1[j];

                    tg_yyyzzzz_xyyz_0[j] = pb_y * tg_yyzzzz_xyyz_0[j] + fr * tg_yyzzzz_xyyz_1[j] + fl1_fx * (tg_yzzzz_xyyz_0[j] - tg_yzzzz_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xyz_1[j];

                    tg_yyyzzzz_xyzz_0[j] = pb_y * tg_yyzzzz_xyzz_0[j] + fr * tg_yyzzzz_xyzz_1[j] + fl1_fx * (tg_yzzzz_xyzz_0[j] - tg_yzzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xzz_1[j];

                    tg_yyyzzzz_xzzz_0[j] = pb_y * tg_yyzzzz_xzzz_0[j] + fr * tg_yyzzzz_xzzz_1[j] + fl1_fx * (tg_yzzzz_xzzz_0[j] - tg_yzzzz_xzzz_1[j] * fl1_fza);

                    tg_yyyzzzz_yyyy_0[j] = pb_y * tg_yyzzzz_yyyy_0[j] + fr * tg_yyzzzz_yyyy_1[j] + fl1_fx * (tg_yzzzz_yyyy_0[j] - tg_yzzzz_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yyzzzz_yyy_1[j];

                    tg_yyyzzzz_yyyz_0[j] = pb_y * tg_yyzzzz_yyyz_0[j] + fr * tg_yyzzzz_yyyz_1[j] + fl1_fx * (tg_yzzzz_yyyz_0[j] - tg_yzzzz_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_yyz_1[j];

                    tg_yyyzzzz_yyzz_0[j] = pb_y * tg_yyzzzz_yyzz_0[j] + fr * tg_yyzzzz_yyzz_1[j] + fl1_fx * (tg_yzzzz_yyzz_0[j] - tg_yzzzz_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_yzz_1[j];

                    tg_yyyzzzz_yzzz_0[j] = pb_y * tg_yyzzzz_yzzz_0[j] + fr * tg_yyzzzz_yzzz_1[j] + fl1_fx * (tg_yzzzz_yzzz_0[j] - tg_yzzzz_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_zzz_1[j];

                    tg_yyyzzzz_zzzz_0[j] = pb_y * tg_yyzzzz_zzzz_0[j] + fr * tg_yyzzzz_zzzz_1[j] + fl1_fx * (tg_yzzzz_zzzz_0[j] - tg_yzzzz_zzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxx_0[j] = pb_y * tg_yzzzzz_xxxx_0[j] + fr * tg_yzzzzz_xxxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxx_0[j] - tg_zzzzz_xxxx_1[j] * fl1_fza);

                    tg_yyzzzzz_xxxy_0[j] = pb_y * tg_yzzzzz_xxxy_0[j] + fr * tg_yzzzzz_xxxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxy_0[j] - tg_zzzzz_xxxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxx_1[j];

                    tg_yyzzzzz_xxxz_0[j] = pb_y * tg_yzzzzz_xxxz_0[j] + fr * tg_yzzzzz_xxxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxxz_0[j] - tg_zzzzz_xxxz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxyy_0[j] = pb_y * tg_yzzzzz_xxyy_0[j] + fr * tg_yzzzzz_xxyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyy_0[j] - tg_zzzzz_xxyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xxy_1[j];

                    tg_yyzzzzz_xxyz_0[j] = pb_y * tg_yzzzzz_xxyz_0[j] + fr * tg_yzzzzz_xxyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxyz_0[j] - tg_zzzzz_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xxz_1[j];

                    tg_yyzzzzz_xxzz_0[j] = pb_y * tg_yzzzzz_xxzz_0[j] + fr * tg_yzzzzz_xxzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxzz_0[j] - tg_zzzzz_xxzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xyyy_0[j] = pb_y * tg_yzzzzz_xyyy_0[j] + fr * tg_yzzzzz_xyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyy_0[j] - tg_zzzzz_xyyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_xyy_1[j];

                    tg_yyzzzzz_xyyz_0[j] = pb_y * tg_yzzzzz_xyyz_0[j] + fr * tg_yzzzzz_xyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyyz_0[j] - tg_zzzzz_xyyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xyz_1[j];

                    tg_yyzzzzz_xyzz_0[j] = pb_y * tg_yzzzzz_xyzz_0[j] + fr * tg_yzzzzz_xyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyzz_0[j] - tg_zzzzz_xyzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xzz_1[j];

                    tg_yyzzzzz_xzzz_0[j] = pb_y * tg_yzzzzz_xzzz_0[j] + fr * tg_yzzzzz_xzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzzz_0[j] - tg_zzzzz_xzzz_1[j] * fl1_fza);

                    tg_yyzzzzz_yyyy_0[j] = pb_y * tg_yzzzzz_yyyy_0[j] + fr * tg_yzzzzz_yyyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyy_0[j] - tg_zzzzz_yyyy_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_yzzzzz_yyy_1[j];

                    tg_yyzzzzz_yyyz_0[j] = pb_y * tg_yzzzzz_yyyz_0[j] + fr * tg_yzzzzz_yyyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyyz_0[j] - tg_zzzzz_yyyz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_yyz_1[j];

                    tg_yyzzzzz_yyzz_0[j] = pb_y * tg_yzzzzz_yyzz_0[j] + fr * tg_yzzzzz_yyzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyzz_0[j] - tg_zzzzz_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_yzz_1[j];

                    tg_yyzzzzz_yzzz_0[j] = pb_y * tg_yzzzzz_yzzz_0[j] + fr * tg_yzzzzz_yzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzzz_0[j] - tg_zzzzz_yzzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_zzz_1[j];

                    tg_yyzzzzz_zzzz_0[j] = pb_y * tg_yzzzzz_zzzz_0[j] + fr * tg_yzzzzz_zzzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzzz_0[j] - tg_zzzzz_zzzz_1[j] * fl1_fza);

                    tg_yzzzzzz_xxxx_0[j] = pb_y * tg_zzzzzz_xxxx_0[j] + fr * tg_zzzzzz_xxxx_1[j];

                    tg_yzzzzzz_xxxy_0[j] = pb_y * tg_zzzzzz_xxxy_0[j] + fr * tg_zzzzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxx_1[j];

                    tg_yzzzzzz_xxxz_0[j] = pb_y * tg_zzzzzz_xxxz_0[j] + fr * tg_zzzzzz_xxxz_1[j];

                    tg_yzzzzzz_xxyy_0[j] = pb_y * tg_zzzzzz_xxyy_0[j] + fr * tg_zzzzzz_xxyy_1[j] + fl1_fxn * tg_zzzzzz_xxy_1[j];

                    tg_yzzzzzz_xxyz_0[j] = pb_y * tg_zzzzzz_xxyz_0[j] + fr * tg_zzzzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xxz_1[j];

                    tg_yzzzzzz_xxzz_0[j] = pb_y * tg_zzzzzz_xxzz_0[j] + fr * tg_zzzzzz_xxzz_1[j];

                    tg_yzzzzzz_xyyy_0[j] = pb_y * tg_zzzzzz_xyyy_0[j] + fr * tg_zzzzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xyy_1[j];

                    tg_yzzzzzz_xyyz_0[j] = pb_y * tg_zzzzzz_xyyz_0[j] + fr * tg_zzzzzz_xyyz_1[j] + fl1_fxn * tg_zzzzzz_xyz_1[j];

                    tg_yzzzzzz_xyzz_0[j] = pb_y * tg_zzzzzz_xyzz_0[j] + fr * tg_zzzzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xzz_1[j];

                    tg_yzzzzzz_xzzz_0[j] = pb_y * tg_zzzzzz_xzzz_0[j] + fr * tg_zzzzzz_xzzz_1[j];

                    tg_yzzzzzz_yyyy_0[j] = pb_y * tg_zzzzzz_yyyy_0[j] + fr * tg_zzzzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_zzzzzz_yyy_1[j];

                    tg_yzzzzzz_yyyz_0[j] = pb_y * tg_zzzzzz_yyyz_0[j] + fr * tg_zzzzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_yyz_1[j];

                    tg_yzzzzzz_yyzz_0[j] = pb_y * tg_zzzzzz_yyzz_0[j] + fr * tg_zzzzzz_yyzz_1[j] + fl1_fxn * tg_zzzzzz_yzz_1[j];

                    tg_yzzzzzz_yzzz_0[j] = pb_y * tg_zzzzzz_yzzz_0[j] + fr * tg_zzzzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zzz_1[j];

                    tg_yzzzzzz_zzzz_0[j] = pb_y * tg_zzzzzz_zzzz_0[j] + fr * tg_zzzzzz_zzzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzzz_xxxx_0[j] = pb_z * tg_zzzzzz_xxxx_0[j] + fr * tg_zzzzzz_xxxx_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxx_0[j] - tg_zzzzz_xxxx_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxy_0[j] = pb_z * tg_zzzzzz_xxxy_0[j] + fr * tg_zzzzzz_xxxy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxy_0[j] - tg_zzzzz_xxxy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxxz_0[j] = pb_z * tg_zzzzzz_xxxz_0[j] + fr * tg_zzzzzz_xxxz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxxz_0[j] - tg_zzzzz_xxxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxx_1[j];

                    tg_zzzzzzz_xxyy_0[j] = pb_z * tg_zzzzzz_xxyy_0[j] + fr * tg_zzzzzz_xxyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyy_0[j] - tg_zzzzz_xxyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxyz_0[j] = pb_z * tg_zzzzzz_xxyz_0[j] + fr * tg_zzzzzz_xxyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxyz_0[j] - tg_zzzzz_xxyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xxy_1[j];

                    tg_zzzzzzz_xxzz_0[j] = pb_z * tg_zzzzzz_xxzz_0[j] + fr * tg_zzzzzz_xxzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxzz_0[j] - tg_zzzzz_xxzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xxz_1[j];

                    tg_zzzzzzz_xyyy_0[j] = pb_z * tg_zzzzzz_xyyy_0[j] + fr * tg_zzzzzz_xyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyy_0[j] - tg_zzzzz_xyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xyyz_0[j] = pb_z * tg_zzzzzz_xyyz_0[j] + fr * tg_zzzzzz_xyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyyz_0[j] - tg_zzzzz_xyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xyy_1[j];

                    tg_zzzzzzz_xyzz_0[j] = pb_z * tg_zzzzzz_xyzz_0[j] + fr * tg_zzzzzz_xyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyzz_0[j] - tg_zzzzz_xyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xyz_1[j];

                    tg_zzzzzzz_xzzz_0[j] = pb_z * tg_zzzzzz_xzzz_0[j] + fr * tg_zzzzzz_xzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xzzz_0[j] - tg_zzzzz_xzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_xzz_1[j];

                    tg_zzzzzzz_yyyy_0[j] = pb_z * tg_zzzzzz_yyyy_0[j] + fr * tg_zzzzzz_yyyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyy_0[j] - tg_zzzzz_yyyy_1[j] * fl1_fza);

                    tg_zzzzzzz_yyyz_0[j] = pb_z * tg_zzzzzz_yyyz_0[j] + fr * tg_zzzzzz_yyyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyyz_0[j] - tg_zzzzz_yyyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_yyy_1[j];

                    tg_zzzzzzz_yyzz_0[j] = pb_z * tg_zzzzzz_yyzz_0[j] + fr * tg_zzzzzz_yyzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyzz_0[j] - tg_zzzzz_yyzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_yyz_1[j];

                    tg_zzzzzzz_yzzz_0[j] = pb_z * tg_zzzzzz_yzzz_0[j] + fr * tg_zzzzzz_yzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yzzz_0[j] - tg_zzzzz_yzzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_yzz_1[j];

                    tg_zzzzzzz_zzzz_0[j] = pb_z * tg_zzzzzz_zzzz_0[j] + fr * tg_zzzzzz_zzzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_zzzz_0[j] - tg_zzzzz_zzzz_1[j] * fl1_fza) + 2.0 * fl1_fxn * tg_zzzzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

