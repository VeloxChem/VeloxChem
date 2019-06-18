//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForHG.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSHSG(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSHSG_0_79(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSHSG_79_158(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSHSG_158_237(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSHSG_237_315(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSHSG_0_79(      CMemBlock2D<double>& primBuffer,
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
                                             {5, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxxx_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx); 

                auto tg_xxxx_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 1); 

                auto tg_xxxx_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 2); 

                auto tg_xxxx_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 3); 

                auto tg_xxxx_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 4); 

                auto tg_xxxx_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 5); 

                auto tg_xxxx_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 6); 

                auto tg_xxxx_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 7); 

                auto tg_xxxx_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 8); 

                auto tg_xxxx_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 9); 

                auto tg_xxxx_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 10); 

                auto tg_xxxx_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 11); 

                auto tg_xxxx_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 12); 

                auto tg_xxxx_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 13); 

                auto tg_xxxx_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 14); 

                auto tg_xxxy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 15); 

                auto tg_xxxy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 16); 

                auto tg_xxxy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 17); 

                auto tg_xxxy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 18); 

                auto tg_xxxy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 19); 

                auto tg_xxxy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 20); 

                auto tg_xxxy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 21); 

                auto tg_xxxy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 22); 

                auto tg_xxxy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 23); 

                auto tg_xxxy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 24); 

                auto tg_xxxy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 25); 

                auto tg_xxxy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 26); 

                auto tg_xxxy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 27); 

                auto tg_xxxy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 28); 

                auto tg_xxxy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 29); 

                auto tg_xxxz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 30); 

                auto tg_xxxz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 31); 

                auto tg_xxxz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 32); 

                auto tg_xxxz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 33); 

                auto tg_xxxz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 34); 

                auto tg_xxxz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 35); 

                auto tg_xxxz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 36); 

                auto tg_xxxz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 37); 

                auto tg_xxxz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 38); 

                auto tg_xxxz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 39); 

                auto tg_xxxz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 40); 

                auto tg_xxxz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 41); 

                auto tg_xxxz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 42); 

                auto tg_xxxz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 43); 

                auto tg_xxxz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 44); 

                auto tg_xxyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 45); 

                auto tg_xxyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 46); 

                auto tg_xxyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 47); 

                auto tg_xxyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 48); 

                auto tg_xxyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 49); 

                auto tg_xxyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 50); 

                auto tg_xxyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 51); 

                auto tg_xxyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 52); 

                auto tg_xxyy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 53); 

                auto tg_xxyy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 54); 

                auto tg_xxyy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 55); 

                auto tg_xxyy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 56); 

                auto tg_xxyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 57); 

                auto tg_xxyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 58); 

                auto tg_xxyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 59); 

                auto tg_xxyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 60); 

                auto tg_xxyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 61); 

                auto tg_xxyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 62); 

                auto tg_xxyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 63); 

                auto tg_xxyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 64); 

                auto tg_xxyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 65); 

                auto tg_xxyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 66); 

                auto tg_xxyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 67); 

                auto tg_xxyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 68); 

                auto tg_xxyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 69); 

                auto tg_xxyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 70); 

                auto tg_xxyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 71); 

                auto tg_xxyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 72); 

                auto tg_xxyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 73); 

                auto tg_xxyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 74); 

                auto tg_xxzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 75); 

                auto tg_xxzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 76); 

                auto tg_xxzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 77); 

                auto tg_xxzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 78); 

                auto tg_xxxx_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx); 

                auto tg_xxxx_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 1); 

                auto tg_xxxx_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 2); 

                auto tg_xxxx_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 3); 

                auto tg_xxxx_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 4); 

                auto tg_xxxx_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 5); 

                auto tg_xxxx_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 6); 

                auto tg_xxxx_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 7); 

                auto tg_xxxx_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 8); 

                auto tg_xxxx_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 9); 

                auto tg_xxxx_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 10); 

                auto tg_xxxx_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 11); 

                auto tg_xxxx_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 12); 

                auto tg_xxxx_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 13); 

                auto tg_xxxx_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 14); 

                auto tg_xxxy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 15); 

                auto tg_xxxy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 16); 

                auto tg_xxxy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 17); 

                auto tg_xxxy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 18); 

                auto tg_xxxy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 19); 

                auto tg_xxxy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 20); 

                auto tg_xxxy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 21); 

                auto tg_xxxy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 22); 

                auto tg_xxxy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 23); 

                auto tg_xxxy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 24); 

                auto tg_xxxy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 25); 

                auto tg_xxxy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 26); 

                auto tg_xxxy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 27); 

                auto tg_xxxy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 28); 

                auto tg_xxxy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 29); 

                auto tg_xxxz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 30); 

                auto tg_xxxz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 31); 

                auto tg_xxxz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 32); 

                auto tg_xxxz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 33); 

                auto tg_xxxz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 34); 

                auto tg_xxxz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 35); 

                auto tg_xxxz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 36); 

                auto tg_xxxz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 37); 

                auto tg_xxxz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 38); 

                auto tg_xxxz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 39); 

                auto tg_xxxz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 40); 

                auto tg_xxxz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 41); 

                auto tg_xxxz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 42); 

                auto tg_xxxz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 43); 

                auto tg_xxxz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 44); 

                auto tg_xxyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 45); 

                auto tg_xxyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 46); 

                auto tg_xxyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 47); 

                auto tg_xxyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 48); 

                auto tg_xxyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 49); 

                auto tg_xxyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 50); 

                auto tg_xxyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 51); 

                auto tg_xxyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 52); 

                auto tg_xxyy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 53); 

                auto tg_xxyy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 54); 

                auto tg_xxyy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 55); 

                auto tg_xxyy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 56); 

                auto tg_xxyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 57); 

                auto tg_xxyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 58); 

                auto tg_xxyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 59); 

                auto tg_xxyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 60); 

                auto tg_xxyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 61); 

                auto tg_xxyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 62); 

                auto tg_xxyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 63); 

                auto tg_xxyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 64); 

                auto tg_xxyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 65); 

                auto tg_xxyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 66); 

                auto tg_xxyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 67); 

                auto tg_xxyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 68); 

                auto tg_xxyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 69); 

                auto tg_xxyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 70); 

                auto tg_xxyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 71); 

                auto tg_xxyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 72); 

                auto tg_xxyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 73); 

                auto tg_xxyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 74); 

                auto tg_xxzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 75); 

                auto tg_xxzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 76); 

                auto tg_xxzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 77); 

                auto tg_xxzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 78); 

                auto tg_xxx_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx); 

                auto tg_xxx_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 1); 

                auto tg_xxx_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 2); 

                auto tg_xxx_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 3); 

                auto tg_xxx_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 4); 

                auto tg_xxx_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 5); 

                auto tg_xxx_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 6); 

                auto tg_xxx_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 7); 

                auto tg_xxx_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 8); 

                auto tg_xxx_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 9); 

                auto tg_xxx_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 10); 

                auto tg_xxx_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 11); 

                auto tg_xxx_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 12); 

                auto tg_xxx_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 13); 

                auto tg_xxx_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 14); 

                auto tg_xxy_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 15); 

                auto tg_xxy_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 16); 

                auto tg_xxy_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 17); 

                auto tg_xxy_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 18); 

                auto tg_xxy_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 19); 

                auto tg_xxy_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 20); 

                auto tg_xxy_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 21); 

                auto tg_xxy_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 22); 

                auto tg_xxy_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 23); 

                auto tg_xxy_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 24); 

                auto tg_xxy_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 25); 

                auto tg_xxy_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 26); 

                auto tg_xxy_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 27); 

                auto tg_xxy_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 28); 

                auto tg_xxy_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 29); 

                auto tg_xxz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 30); 

                auto tg_xxz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 31); 

                auto tg_xxz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 32); 

                auto tg_xxz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 33); 

                auto tg_xxz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 34); 

                auto tg_xxz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 35); 

                auto tg_xxz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 36); 

                auto tg_xxz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 37); 

                auto tg_xxz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 38); 

                auto tg_xxz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 39); 

                auto tg_xxz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 40); 

                auto tg_xxz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 41); 

                auto tg_xxz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 42); 

                auto tg_xxz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 43); 

                auto tg_xxz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 44); 

                auto tg_xyy_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 45); 

                auto tg_xyy_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 46); 

                auto tg_xyy_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 47); 

                auto tg_xyy_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 48); 

                auto tg_xyy_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 49); 

                auto tg_xyy_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 50); 

                auto tg_xyy_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 51); 

                auto tg_xyy_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 52); 

                auto tg_xyy_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 53); 

                auto tg_xyy_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 54); 

                auto tg_xyy_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 55); 

                auto tg_xyy_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 56); 

                auto tg_xyy_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 57); 

                auto tg_xyy_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 58); 

                auto tg_xyy_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 59); 

                auto tg_xyz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 60); 

                auto tg_xyz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 61); 

                auto tg_xyz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 62); 

                auto tg_xyz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 63); 

                auto tg_xyz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 64); 

                auto tg_xyz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 65); 

                auto tg_xyz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 66); 

                auto tg_xyz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 67); 

                auto tg_xyz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 68); 

                auto tg_xyz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 69); 

                auto tg_xyz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 70); 

                auto tg_xyz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 71); 

                auto tg_xyz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 72); 

                auto tg_xyz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 73); 

                auto tg_xyz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 74); 

                auto tg_xzz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 75); 

                auto tg_xzz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 76); 

                auto tg_xzz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 77); 

                auto tg_xzz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 78); 

                auto tg_xxx_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx); 

                auto tg_xxx_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 1); 

                auto tg_xxx_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 2); 

                auto tg_xxx_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 3); 

                auto tg_xxx_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 4); 

                auto tg_xxx_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 5); 

                auto tg_xxx_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 6); 

                auto tg_xxx_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 7); 

                auto tg_xxx_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 8); 

                auto tg_xxx_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 9); 

                auto tg_xxx_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 10); 

                auto tg_xxx_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 11); 

                auto tg_xxx_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 12); 

                auto tg_xxx_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 13); 

                auto tg_xxx_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 14); 

                auto tg_xxy_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 15); 

                auto tg_xxy_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 16); 

                auto tg_xxy_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 17); 

                auto tg_xxy_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 18); 

                auto tg_xxy_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 19); 

                auto tg_xxy_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 20); 

                auto tg_xxy_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 21); 

                auto tg_xxy_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 22); 

                auto tg_xxy_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 23); 

                auto tg_xxy_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 24); 

                auto tg_xxy_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 25); 

                auto tg_xxy_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 26); 

                auto tg_xxy_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 27); 

                auto tg_xxy_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 28); 

                auto tg_xxy_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 29); 

                auto tg_xxz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 30); 

                auto tg_xxz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 31); 

                auto tg_xxz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 32); 

                auto tg_xxz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 33); 

                auto tg_xxz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 34); 

                auto tg_xxz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 35); 

                auto tg_xxz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 36); 

                auto tg_xxz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 37); 

                auto tg_xxz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 38); 

                auto tg_xxz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 39); 

                auto tg_xxz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 40); 

                auto tg_xxz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 41); 

                auto tg_xxz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 42); 

                auto tg_xxz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 43); 

                auto tg_xxz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 44); 

                auto tg_xyy_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 45); 

                auto tg_xyy_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 46); 

                auto tg_xyy_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 47); 

                auto tg_xyy_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 48); 

                auto tg_xyy_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 49); 

                auto tg_xyy_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 50); 

                auto tg_xyy_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 51); 

                auto tg_xyy_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 52); 

                auto tg_xyy_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 53); 

                auto tg_xyy_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 54); 

                auto tg_xyy_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 55); 

                auto tg_xyy_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 56); 

                auto tg_xyy_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 57); 

                auto tg_xyy_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 58); 

                auto tg_xyy_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 59); 

                auto tg_xyz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 60); 

                auto tg_xyz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 61); 

                auto tg_xyz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 62); 

                auto tg_xyz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 63); 

                auto tg_xyz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 64); 

                auto tg_xyz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 65); 

                auto tg_xyz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 66); 

                auto tg_xyz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 67); 

                auto tg_xyz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 68); 

                auto tg_xyz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 69); 

                auto tg_xyz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 70); 

                auto tg_xyz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 71); 

                auto tg_xyz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 72); 

                auto tg_xyz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 73); 

                auto tg_xyz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 74); 

                auto tg_xzz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 75); 

                auto tg_xzz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 76); 

                auto tg_xzz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 77); 

                auto tg_xzz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 78); 

                auto tg_xxxx_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx); 

                auto tg_xxxx_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 1); 

                auto tg_xxxx_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 2); 

                auto tg_xxxx_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 3); 

                auto tg_xxxx_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 4); 

                auto tg_xxxx_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 5); 

                auto tg_xxxx_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 6); 

                auto tg_xxxx_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 7); 

                auto tg_xxxx_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 8); 

                auto tg_xxxx_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 9); 

                auto tg_xxxy_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 10); 

                auto tg_xxxy_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 11); 

                auto tg_xxxy_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 12); 

                auto tg_xxxy_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 13); 

                auto tg_xxxy_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 14); 

                auto tg_xxxy_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 15); 

                auto tg_xxxy_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 16); 

                auto tg_xxxy_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 17); 

                auto tg_xxxy_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 18); 

                auto tg_xxxy_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 19); 

                auto tg_xxxz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 20); 

                auto tg_xxxz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 21); 

                auto tg_xxxz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 22); 

                auto tg_xxxz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 23); 

                auto tg_xxxz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 24); 

                auto tg_xxxz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 25); 

                auto tg_xxxz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 26); 

                auto tg_xxxz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 27); 

                auto tg_xxxz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 28); 

                auto tg_xxxz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 29); 

                auto tg_xxyy_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 30); 

                auto tg_xxyy_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 31); 

                auto tg_xxyy_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 32); 

                auto tg_xxyy_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 33); 

                auto tg_xxyy_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 34); 

                auto tg_xxyy_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 35); 

                auto tg_xxyy_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 36); 

                auto tg_xxyy_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 37); 

                auto tg_xxyy_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 38); 

                auto tg_xxyy_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 39); 

                auto tg_xxyz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 40); 

                auto tg_xxyz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 41); 

                auto tg_xxyz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 42); 

                auto tg_xxyz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 43); 

                auto tg_xxyz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 44); 

                auto tg_xxyz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 45); 

                auto tg_xxyz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 46); 

                auto tg_xxyz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 47); 

                auto tg_xxyz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 48); 

                auto tg_xxyz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 49); 

                auto tg_xxzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 50); 

                auto tg_xxzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 51); 

                auto tg_xxzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 52); 

                auto tg_xxzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 53); 

                // set up pointers to integrals

                auto tg_xxxxx_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx); 

                auto tg_xxxxx_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 1); 

                auto tg_xxxxx_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 2); 

                auto tg_xxxxx_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 3); 

                auto tg_xxxxx_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 4); 

                auto tg_xxxxx_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 5); 

                auto tg_xxxxx_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 6); 

                auto tg_xxxxx_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 7); 

                auto tg_xxxxx_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 8); 

                auto tg_xxxxx_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 9); 

                auto tg_xxxxx_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 10); 

                auto tg_xxxxx_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 11); 

                auto tg_xxxxx_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 12); 

                auto tg_xxxxx_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 13); 

                auto tg_xxxxx_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 14); 

                auto tg_xxxxy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 15); 

                auto tg_xxxxy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 16); 

                auto tg_xxxxy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 17); 

                auto tg_xxxxy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 18); 

                auto tg_xxxxy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 19); 

                auto tg_xxxxy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 20); 

                auto tg_xxxxy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 21); 

                auto tg_xxxxy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 22); 

                auto tg_xxxxy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 23); 

                auto tg_xxxxy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 24); 

                auto tg_xxxxy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 25); 

                auto tg_xxxxy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 26); 

                auto tg_xxxxy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 27); 

                auto tg_xxxxy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 28); 

                auto tg_xxxxy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 29); 

                auto tg_xxxxz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 30); 

                auto tg_xxxxz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 31); 

                auto tg_xxxxz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 32); 

                auto tg_xxxxz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 33); 

                auto tg_xxxxz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 34); 

                auto tg_xxxxz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 35); 

                auto tg_xxxxz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 36); 

                auto tg_xxxxz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 37); 

                auto tg_xxxxz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 38); 

                auto tg_xxxxz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 39); 

                auto tg_xxxxz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 40); 

                auto tg_xxxxz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 41); 

                auto tg_xxxxz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 42); 

                auto tg_xxxxz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 43); 

                auto tg_xxxxz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 44); 

                auto tg_xxxyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 45); 

                auto tg_xxxyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 46); 

                auto tg_xxxyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 47); 

                auto tg_xxxyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 48); 

                auto tg_xxxyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 49); 

                auto tg_xxxyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 50); 

                auto tg_xxxyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 51); 

                auto tg_xxxyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 52); 

                auto tg_xxxyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 53); 

                auto tg_xxxyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 54); 

                auto tg_xxxyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 55); 

                auto tg_xxxyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 56); 

                auto tg_xxxyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 57); 

                auto tg_xxxyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 58); 

                auto tg_xxxyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 59); 

                auto tg_xxxyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 60); 

                auto tg_xxxyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 61); 

                auto tg_xxxyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 62); 

                auto tg_xxxyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 63); 

                auto tg_xxxyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 64); 

                auto tg_xxxyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 65); 

                auto tg_xxxyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 66); 

                auto tg_xxxyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 67); 

                auto tg_xxxyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 68); 

                auto tg_xxxyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 69); 

                auto tg_xxxyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 70); 

                auto tg_xxxyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 71); 

                auto tg_xxxyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 72); 

                auto tg_xxxyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 73); 

                auto tg_xxxyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 74); 

                auto tg_xxxzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 75); 

                auto tg_xxxzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 76); 

                auto tg_xxxzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 77); 

                auto tg_xxxzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 78); 

                // Batch of Integrals (0,79)

                #pragma omp simd aligned(fxn, fza, tg_xxx_xxxx_0, tg_xxx_xxxx_1, tg_xxx_xxxy_0, tg_xxx_xxxy_1, \
                                         tg_xxx_xxxz_0, tg_xxx_xxxz_1, tg_xxx_xxyy_0, tg_xxx_xxyy_1, tg_xxx_xxyz_0, \
                                         tg_xxx_xxyz_1, tg_xxx_xxzz_0, tg_xxx_xxzz_1, tg_xxx_xyyy_0, tg_xxx_xyyy_1, \
                                         tg_xxx_xyyz_0, tg_xxx_xyyz_1, tg_xxx_xyzz_0, tg_xxx_xyzz_1, tg_xxx_xzzz_0, \
                                         tg_xxx_xzzz_1, tg_xxx_yyyy_0, tg_xxx_yyyy_1, tg_xxx_yyyz_0, tg_xxx_yyyz_1, \
                                         tg_xxx_yyzz_0, tg_xxx_yyzz_1, tg_xxx_yzzz_0, tg_xxx_yzzz_1, tg_xxx_zzzz_0, \
                                         tg_xxx_zzzz_1, tg_xxxx_xxx_1, tg_xxxx_xxxx_0, tg_xxxx_xxxx_1, tg_xxxx_xxxy_0, \
                                         tg_xxxx_xxxy_1, tg_xxxx_xxxz_0, tg_xxxx_xxxz_1, tg_xxxx_xxy_1, tg_xxxx_xxyy_0, \
                                         tg_xxxx_xxyy_1, tg_xxxx_xxyz_0, tg_xxxx_xxyz_1, tg_xxxx_xxz_1, tg_xxxx_xxzz_0, \
                                         tg_xxxx_xxzz_1, tg_xxxx_xyy_1, tg_xxxx_xyyy_0, tg_xxxx_xyyy_1, tg_xxxx_xyyz_0, \
                                         tg_xxxx_xyyz_1, tg_xxxx_xyz_1, tg_xxxx_xyzz_0, tg_xxxx_xyzz_1, tg_xxxx_xzz_1, \
                                         tg_xxxx_xzzz_0, tg_xxxx_xzzz_1, tg_xxxx_yyy_1, tg_xxxx_yyyy_0, tg_xxxx_yyyy_1, \
                                         tg_xxxx_yyyz_0, tg_xxxx_yyyz_1, tg_xxxx_yyz_1, tg_xxxx_yyzz_0, tg_xxxx_yyzz_1, \
                                         tg_xxxx_yzz_1, tg_xxxx_yzzz_0, tg_xxxx_yzzz_1, tg_xxxx_zzz_1, tg_xxxx_zzzz_0, \
                                         tg_xxxx_zzzz_1, tg_xxxxx_xxxx_0, tg_xxxxx_xxxy_0, tg_xxxxx_xxxz_0, tg_xxxxx_xxyy_0, \
                                         tg_xxxxx_xxyz_0, tg_xxxxx_xxzz_0, tg_xxxxx_xyyy_0, tg_xxxxx_xyyz_0, tg_xxxxx_xyzz_0, \
                                         tg_xxxxx_xzzz_0, tg_xxxxx_yyyy_0, tg_xxxxx_yyyz_0, tg_xxxxx_yyzz_0, tg_xxxxx_yzzz_0, \
                                         tg_xxxxx_zzzz_0, tg_xxxxy_xxxx_0, tg_xxxxy_xxxy_0, tg_xxxxy_xxxz_0, tg_xxxxy_xxyy_0, \
                                         tg_xxxxy_xxyz_0, tg_xxxxy_xxzz_0, tg_xxxxy_xyyy_0, tg_xxxxy_xyyz_0, tg_xxxxy_xyzz_0, \
                                         tg_xxxxy_xzzz_0, tg_xxxxy_yyyy_0, tg_xxxxy_yyyz_0, tg_xxxxy_yyzz_0, tg_xxxxy_yzzz_0, \
                                         tg_xxxxy_zzzz_0, tg_xxxxz_xxxx_0, tg_xxxxz_xxxy_0, tg_xxxxz_xxxz_0, tg_xxxxz_xxyy_0, \
                                         tg_xxxxz_xxyz_0, tg_xxxxz_xxzz_0, tg_xxxxz_xyyy_0, tg_xxxxz_xyyz_0, tg_xxxxz_xyzz_0, \
                                         tg_xxxxz_xzzz_0, tg_xxxxz_yyyy_0, tg_xxxxz_yyyz_0, tg_xxxxz_yyzz_0, tg_xxxxz_yzzz_0, \
                                         tg_xxxxz_zzzz_0, tg_xxxy_xxx_1, tg_xxxy_xxxx_0, tg_xxxy_xxxx_1, tg_xxxy_xxxy_0, \
                                         tg_xxxy_xxxy_1, tg_xxxy_xxxz_0, tg_xxxy_xxxz_1, tg_xxxy_xxy_1, tg_xxxy_xxyy_0, \
                                         tg_xxxy_xxyy_1, tg_xxxy_xxyz_0, tg_xxxy_xxyz_1, tg_xxxy_xxz_1, tg_xxxy_xxzz_0, \
                                         tg_xxxy_xxzz_1, tg_xxxy_xyy_1, tg_xxxy_xyyy_0, tg_xxxy_xyyy_1, tg_xxxy_xyyz_0, \
                                         tg_xxxy_xyyz_1, tg_xxxy_xyz_1, tg_xxxy_xyzz_0, tg_xxxy_xyzz_1, tg_xxxy_xzz_1, \
                                         tg_xxxy_xzzz_0, tg_xxxy_xzzz_1, tg_xxxy_yyy_1, tg_xxxy_yyyy_0, tg_xxxy_yyyy_1, \
                                         tg_xxxy_yyyz_0, tg_xxxy_yyyz_1, tg_xxxy_yyz_1, tg_xxxy_yyzz_0, tg_xxxy_yyzz_1, \
                                         tg_xxxy_yzz_1, tg_xxxy_yzzz_0, tg_xxxy_yzzz_1, tg_xxxy_zzz_1, tg_xxxy_zzzz_0, \
                                         tg_xxxy_zzzz_1, tg_xxxyy_xxxx_0, tg_xxxyy_xxxy_0, tg_xxxyy_xxxz_0, tg_xxxyy_xxyy_0, \
                                         tg_xxxyy_xxyz_0, tg_xxxyy_xxzz_0, tg_xxxyy_xyyy_0, tg_xxxyy_xyyz_0, tg_xxxyy_xyzz_0, \
                                         tg_xxxyy_xzzz_0, tg_xxxyy_yyyy_0, tg_xxxyy_yyyz_0, tg_xxxyy_yyzz_0, tg_xxxyy_yzzz_0, \
                                         tg_xxxyy_zzzz_0, tg_xxxyz_xxxx_0, tg_xxxyz_xxxy_0, tg_xxxyz_xxxz_0, tg_xxxyz_xxyy_0, \
                                         tg_xxxyz_xxyz_0, tg_xxxyz_xxzz_0, tg_xxxyz_xyyy_0, tg_xxxyz_xyyz_0, tg_xxxyz_xyzz_0, \
                                         tg_xxxyz_xzzz_0, tg_xxxyz_yyyy_0, tg_xxxyz_yyyz_0, tg_xxxyz_yyzz_0, tg_xxxyz_yzzz_0, \
                                         tg_xxxyz_zzzz_0, tg_xxxz_xxx_1, tg_xxxz_xxxx_0, tg_xxxz_xxxx_1, tg_xxxz_xxxy_0, \
                                         tg_xxxz_xxxy_1, tg_xxxz_xxxz_0, tg_xxxz_xxxz_1, tg_xxxz_xxy_1, tg_xxxz_xxyy_0, \
                                         tg_xxxz_xxyy_1, tg_xxxz_xxyz_0, tg_xxxz_xxyz_1, tg_xxxz_xxz_1, tg_xxxz_xxzz_0, \
                                         tg_xxxz_xxzz_1, tg_xxxz_xyy_1, tg_xxxz_xyyy_0, tg_xxxz_xyyy_1, tg_xxxz_xyyz_0, \
                                         tg_xxxz_xyyz_1, tg_xxxz_xyz_1, tg_xxxz_xyzz_0, tg_xxxz_xyzz_1, tg_xxxz_xzz_1, \
                                         tg_xxxz_xzzz_0, tg_xxxz_xzzz_1, tg_xxxz_yyy_1, tg_xxxz_yyyy_0, tg_xxxz_yyyy_1, \
                                         tg_xxxz_yyyz_0, tg_xxxz_yyyz_1, tg_xxxz_yyz_1, tg_xxxz_yyzz_0, tg_xxxz_yyzz_1, \
                                         tg_xxxz_yzz_1, tg_xxxz_yzzz_0, tg_xxxz_yzzz_1, tg_xxxz_zzz_1, tg_xxxz_zzzz_0, \
                                         tg_xxxz_zzzz_1, tg_xxxzz_xxxx_0, tg_xxxzz_xxxy_0, tg_xxxzz_xxxz_0, tg_xxxzz_xxyy_0, \
                                         tg_xxy_xxxx_0, tg_xxy_xxxx_1, tg_xxy_xxxy_0, tg_xxy_xxxy_1, tg_xxy_xxxz_0, \
                                         tg_xxy_xxxz_1, tg_xxy_xxyy_0, tg_xxy_xxyy_1, tg_xxy_xxyz_0, tg_xxy_xxyz_1, \
                                         tg_xxy_xxzz_0, tg_xxy_xxzz_1, tg_xxy_xyyy_0, tg_xxy_xyyy_1, tg_xxy_xyyz_0, \
                                         tg_xxy_xyyz_1, tg_xxy_xyzz_0, tg_xxy_xyzz_1, tg_xxy_xzzz_0, tg_xxy_xzzz_1, \
                                         tg_xxy_yyyy_0, tg_xxy_yyyy_1, tg_xxy_yyyz_0, tg_xxy_yyyz_1, tg_xxy_yyzz_0, \
                                         tg_xxy_yyzz_1, tg_xxy_yzzz_0, tg_xxy_yzzz_1, tg_xxy_zzzz_0, tg_xxy_zzzz_1, \
                                         tg_xxyy_xxx_1, tg_xxyy_xxxx_0, tg_xxyy_xxxx_1, tg_xxyy_xxxy_0, tg_xxyy_xxxy_1, \
                                         tg_xxyy_xxxz_0, tg_xxyy_xxxz_1, tg_xxyy_xxy_1, tg_xxyy_xxyy_0, tg_xxyy_xxyy_1, \
                                         tg_xxyy_xxyz_0, tg_xxyy_xxyz_1, tg_xxyy_xxz_1, tg_xxyy_xxzz_0, tg_xxyy_xxzz_1, \
                                         tg_xxyy_xyy_1, tg_xxyy_xyyy_0, tg_xxyy_xyyy_1, tg_xxyy_xyyz_0, tg_xxyy_xyyz_1, \
                                         tg_xxyy_xyz_1, tg_xxyy_xyzz_0, tg_xxyy_xyzz_1, tg_xxyy_xzz_1, tg_xxyy_xzzz_0, \
                                         tg_xxyy_xzzz_1, tg_xxyy_yyy_1, tg_xxyy_yyyy_0, tg_xxyy_yyyy_1, tg_xxyy_yyyz_0, \
                                         tg_xxyy_yyyz_1, tg_xxyy_yyz_1, tg_xxyy_yyzz_0, tg_xxyy_yyzz_1, tg_xxyy_yzz_1, \
                                         tg_xxyy_yzzz_0, tg_xxyy_yzzz_1, tg_xxyy_zzz_1, tg_xxyy_zzzz_0, tg_xxyy_zzzz_1, \
                                         tg_xxyz_xxx_1, tg_xxyz_xxxx_0, tg_xxyz_xxxx_1, tg_xxyz_xxxy_0, tg_xxyz_xxxy_1, \
                                         tg_xxyz_xxxz_0, tg_xxyz_xxxz_1, tg_xxyz_xxy_1, tg_xxyz_xxyy_0, tg_xxyz_xxyy_1, \
                                         tg_xxyz_xxyz_0, tg_xxyz_xxyz_1, tg_xxyz_xxz_1, tg_xxyz_xxzz_0, tg_xxyz_xxzz_1, \
                                         tg_xxyz_xyy_1, tg_xxyz_xyyy_0, tg_xxyz_xyyy_1, tg_xxyz_xyyz_0, tg_xxyz_xyyz_1, \
                                         tg_xxyz_xyz_1, tg_xxyz_xyzz_0, tg_xxyz_xyzz_1, tg_xxyz_xzz_1, tg_xxyz_xzzz_0, \
                                         tg_xxyz_xzzz_1, tg_xxyz_yyy_1, tg_xxyz_yyyy_0, tg_xxyz_yyyy_1, tg_xxyz_yyyz_0, \
                                         tg_xxyz_yyyz_1, tg_xxyz_yyz_1, tg_xxyz_yyzz_0, tg_xxyz_yyzz_1, tg_xxyz_yzz_1, \
                                         tg_xxyz_yzzz_0, tg_xxyz_yzzz_1, tg_xxyz_zzz_1, tg_xxyz_zzzz_0, tg_xxyz_zzzz_1, \
                                         tg_xxz_xxxx_0, tg_xxz_xxxx_1, tg_xxz_xxxy_0, tg_xxz_xxxy_1, tg_xxz_xxxz_0, \
                                         tg_xxz_xxxz_1, tg_xxz_xxyy_0, tg_xxz_xxyy_1, tg_xxz_xxyz_0, tg_xxz_xxyz_1, \
                                         tg_xxz_xxzz_0, tg_xxz_xxzz_1, tg_xxz_xyyy_0, tg_xxz_xyyy_1, tg_xxz_xyyz_0, \
                                         tg_xxz_xyyz_1, tg_xxz_xyzz_0, tg_xxz_xyzz_1, tg_xxz_xzzz_0, tg_xxz_xzzz_1, \
                                         tg_xxz_yyyy_0, tg_xxz_yyyy_1, tg_xxz_yyyz_0, tg_xxz_yyyz_1, tg_xxz_yyzz_0, \
                                         tg_xxz_yyzz_1, tg_xxz_yzzz_0, tg_xxz_yzzz_1, tg_xxz_zzzz_0, tg_xxz_zzzz_1, \
                                         tg_xxzz_xxx_1, tg_xxzz_xxxx_0, tg_xxzz_xxxx_1, tg_xxzz_xxxy_0, tg_xxzz_xxxy_1, \
                                         tg_xxzz_xxxz_0, tg_xxzz_xxxz_1, tg_xxzz_xxy_1, tg_xxzz_xxyy_0, tg_xxzz_xxyy_1, \
                                         tg_xxzz_xxz_1, tg_xxzz_xyy_1, tg_xyy_xxxx_0, tg_xyy_xxxx_1, tg_xyy_xxxy_0, \
                                         tg_xyy_xxxy_1, tg_xyy_xxxz_0, tg_xyy_xxxz_1, tg_xyy_xxyy_0, tg_xyy_xxyy_1, \
                                         tg_xyy_xxyz_0, tg_xyy_xxyz_1, tg_xyy_xxzz_0, tg_xyy_xxzz_1, tg_xyy_xyyy_0, \
                                         tg_xyy_xyyy_1, tg_xyy_xyyz_0, tg_xyy_xyyz_1, tg_xyy_xyzz_0, tg_xyy_xyzz_1, \
                                         tg_xyy_xzzz_0, tg_xyy_xzzz_1, tg_xyy_yyyy_0, tg_xyy_yyyy_1, tg_xyy_yyyz_0, \
                                         tg_xyy_yyyz_1, tg_xyy_yyzz_0, tg_xyy_yyzz_1, tg_xyy_yzzz_0, tg_xyy_yzzz_1, \
                                         tg_xyy_zzzz_0, tg_xyy_zzzz_1, tg_xyz_xxxx_0, tg_xyz_xxxx_1, tg_xyz_xxxy_0, \
                                         tg_xyz_xxxy_1, tg_xyz_xxxz_0, tg_xyz_xxxz_1, tg_xyz_xxyy_0, tg_xyz_xxyy_1, \
                                         tg_xyz_xxyz_0, tg_xyz_xxyz_1, tg_xyz_xxzz_0, tg_xyz_xxzz_1, tg_xyz_xyyy_0, \
                                         tg_xyz_xyyy_1, tg_xyz_xyyz_0, tg_xyz_xyyz_1, tg_xyz_xyzz_0, tg_xyz_xyzz_1, \
                                         tg_xyz_xzzz_0, tg_xyz_xzzz_1, tg_xyz_yyyy_0, tg_xyz_yyyy_1, tg_xyz_yyyz_0, \
                                         tg_xyz_yyyz_1, tg_xyz_yyzz_0, tg_xyz_yyzz_1, tg_xyz_yzzz_0, tg_xyz_yzzz_1, \
                                         tg_xyz_zzzz_0, tg_xyz_zzzz_1, tg_xzz_xxxx_0, tg_xzz_xxxx_1, tg_xzz_xxxy_0, \
                                         tg_xzz_xxxy_1, tg_xzz_xxxz_0, tg_xzz_xxxz_1, tg_xzz_xxyy_0, tg_xzz_xxyy_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxx_xxxx_0[j] = pb_x * tg_xxxx_xxxx_0[j] + wp_x[j] * tg_xxxx_xxxx_1[j] + 2.0 * fl1_fx * tg_xxx_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxx_xxx_1[j];

                    tg_xxxxx_xxxy_0[j] = pb_x * tg_xxxx_xxxy_0[j] + wp_x[j] * tg_xxxx_xxxy_1[j] + 2.0 * fl1_fx * tg_xxx_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxx_xxy_1[j];

                    tg_xxxxx_xxxz_0[j] = pb_x * tg_xxxx_xxxz_0[j] + wp_x[j] * tg_xxxx_xxxz_1[j] + 2.0 * fl1_fx * tg_xxx_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxx_xxz_1[j];

                    tg_xxxxx_xxyy_0[j] = pb_x * tg_xxxx_xxyy_0[j] + wp_x[j] * tg_xxxx_xxyy_1[j] + 2.0 * fl1_fx * tg_xxx_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxyy_1[j] + fl1_fxn * tg_xxxx_xyy_1[j];

                    tg_xxxxx_xxyz_0[j] = pb_x * tg_xxxx_xxyz_0[j] + wp_x[j] * tg_xxxx_xxyz_1[j] + 2.0 * fl1_fx * tg_xxx_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxyz_1[j] + fl1_fxn * tg_xxxx_xyz_1[j];

                    tg_xxxxx_xxzz_0[j] = pb_x * tg_xxxx_xxzz_0[j] + wp_x[j] * tg_xxxx_xxzz_1[j] + 2.0 * fl1_fx * tg_xxx_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxzz_1[j] + fl1_fxn * tg_xxxx_xzz_1[j];

                    tg_xxxxx_xyyy_0[j] = pb_x * tg_xxxx_xyyy_0[j] + wp_x[j] * tg_xxxx_xyyy_1[j] + 2.0 * fl1_fx * tg_xxx_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxx_yyy_1[j];

                    tg_xxxxx_xyyz_0[j] = pb_x * tg_xxxx_xyyz_0[j] + wp_x[j] * tg_xxxx_xyyz_1[j] + 2.0 * fl1_fx * tg_xxx_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yyz_1[j];

                    tg_xxxxx_xyzz_0[j] = pb_x * tg_xxxx_xyzz_0[j] + wp_x[j] * tg_xxxx_xyzz_1[j] + 2.0 * fl1_fx * tg_xxx_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yzz_1[j];

                    tg_xxxxx_xzzz_0[j] = pb_x * tg_xxxx_xzzz_0[j] + wp_x[j] * tg_xxxx_xzzz_1[j] + 2.0 * fl1_fx * tg_xxx_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_zzz_1[j];

                    tg_xxxxx_yyyy_0[j] = pb_x * tg_xxxx_yyyy_0[j] + wp_x[j] * tg_xxxx_yyyy_1[j] + 2.0 * fl1_fx * tg_xxx_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyyy_1[j];

                    tg_xxxxx_yyyz_0[j] = pb_x * tg_xxxx_yyyz_0[j] + wp_x[j] * tg_xxxx_yyyz_1[j] + 2.0 * fl1_fx * tg_xxx_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyyz_1[j];

                    tg_xxxxx_yyzz_0[j] = pb_x * tg_xxxx_yyzz_0[j] + wp_x[j] * tg_xxxx_yyzz_1[j] + 2.0 * fl1_fx * tg_xxx_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyzz_1[j];

                    tg_xxxxx_yzzz_0[j] = pb_x * tg_xxxx_yzzz_0[j] + wp_x[j] * tg_xxxx_yzzz_1[j] + 2.0 * fl1_fx * tg_xxx_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yzzz_1[j];

                    tg_xxxxx_zzzz_0[j] = pb_x * tg_xxxx_zzzz_0[j] + wp_x[j] * tg_xxxx_zzzz_1[j] + 2.0 * fl1_fx * tg_xxx_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_zzzz_1[j];

                    tg_xxxxy_xxxx_0[j] = pb_x * tg_xxxy_xxxx_0[j] + wp_x[j] * tg_xxxy_xxxx_1[j] + 1.5 * fl1_fx * tg_xxy_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxy_xxx_1[j];

                    tg_xxxxy_xxxy_0[j] = pb_x * tg_xxxy_xxxy_0[j] + wp_x[j] * tg_xxxy_xxxy_1[j] + 1.5 * fl1_fx * tg_xxy_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxy_xxy_1[j];

                    tg_xxxxy_xxxz_0[j] = pb_x * tg_xxxy_xxxz_0[j] + wp_x[j] * tg_xxxy_xxxz_1[j] + 1.5 * fl1_fx * tg_xxy_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxy_xxz_1[j];

                    tg_xxxxy_xxyy_0[j] = pb_x * tg_xxxy_xxyy_0[j] + wp_x[j] * tg_xxxy_xxyy_1[j] + 1.5 * fl1_fx * tg_xxy_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxyy_1[j] + fl1_fxn * tg_xxxy_xyy_1[j];

                    tg_xxxxy_xxyz_0[j] = pb_x * tg_xxxy_xxyz_0[j] + wp_x[j] * tg_xxxy_xxyz_1[j] + 1.5 * fl1_fx * tg_xxy_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxyz_1[j] + fl1_fxn * tg_xxxy_xyz_1[j];

                    tg_xxxxy_xxzz_0[j] = pb_x * tg_xxxy_xxzz_0[j] + wp_x[j] * tg_xxxy_xxzz_1[j] + 1.5 * fl1_fx * tg_xxy_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxzz_1[j] + fl1_fxn * tg_xxxy_xzz_1[j];

                    tg_xxxxy_xyyy_0[j] = pb_x * tg_xxxy_xyyy_0[j] + wp_x[j] * tg_xxxy_xyyy_1[j] + 1.5 * fl1_fx * tg_xxy_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxy_yyy_1[j];

                    tg_xxxxy_xyyz_0[j] = pb_x * tg_xxxy_xyyz_0[j] + wp_x[j] * tg_xxxy_xyyz_1[j] + 1.5 * fl1_fx * tg_xxy_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yyz_1[j];

                    tg_xxxxy_xyzz_0[j] = pb_x * tg_xxxy_xyzz_0[j] + wp_x[j] * tg_xxxy_xyzz_1[j] + 1.5 * fl1_fx * tg_xxy_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yzz_1[j];

                    tg_xxxxy_xzzz_0[j] = pb_x * tg_xxxy_xzzz_0[j] + wp_x[j] * tg_xxxy_xzzz_1[j] + 1.5 * fl1_fx * tg_xxy_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_zzz_1[j];

                    tg_xxxxy_yyyy_0[j] = pb_x * tg_xxxy_yyyy_0[j] + wp_x[j] * tg_xxxy_yyyy_1[j] + 1.5 * fl1_fx * tg_xxy_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyyy_1[j];

                    tg_xxxxy_yyyz_0[j] = pb_x * tg_xxxy_yyyz_0[j] + wp_x[j] * tg_xxxy_yyyz_1[j] + 1.5 * fl1_fx * tg_xxy_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyyz_1[j];

                    tg_xxxxy_yyzz_0[j] = pb_x * tg_xxxy_yyzz_0[j] + wp_x[j] * tg_xxxy_yyzz_1[j] + 1.5 * fl1_fx * tg_xxy_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyzz_1[j];

                    tg_xxxxy_yzzz_0[j] = pb_x * tg_xxxy_yzzz_0[j] + wp_x[j] * tg_xxxy_yzzz_1[j] + 1.5 * fl1_fx * tg_xxy_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yzzz_1[j];

                    tg_xxxxy_zzzz_0[j] = pb_x * tg_xxxy_zzzz_0[j] + wp_x[j] * tg_xxxy_zzzz_1[j] + 1.5 * fl1_fx * tg_xxy_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_zzzz_1[j];

                    tg_xxxxz_xxxx_0[j] = pb_x * tg_xxxz_xxxx_0[j] + wp_x[j] * tg_xxxz_xxxx_1[j] + 1.5 * fl1_fx * tg_xxz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxxz_xxx_1[j];

                    tg_xxxxz_xxxy_0[j] = pb_x * tg_xxxz_xxxy_0[j] + wp_x[j] * tg_xxxz_xxxy_1[j] + 1.5 * fl1_fx * tg_xxz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxxz_xxy_1[j];

                    tg_xxxxz_xxxz_0[j] = pb_x * tg_xxxz_xxxz_0[j] + wp_x[j] * tg_xxxz_xxxz_1[j] + 1.5 * fl1_fx * tg_xxz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxxz_xxz_1[j];

                    tg_xxxxz_xxyy_0[j] = pb_x * tg_xxxz_xxyy_0[j] + wp_x[j] * tg_xxxz_xxyy_1[j] + 1.5 * fl1_fx * tg_xxz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxyy_1[j] + fl1_fxn * tg_xxxz_xyy_1[j];

                    tg_xxxxz_xxyz_0[j] = pb_x * tg_xxxz_xxyz_0[j] + wp_x[j] * tg_xxxz_xxyz_1[j] + 1.5 * fl1_fx * tg_xxz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxyz_1[j] + fl1_fxn * tg_xxxz_xyz_1[j];

                    tg_xxxxz_xxzz_0[j] = pb_x * tg_xxxz_xxzz_0[j] + wp_x[j] * tg_xxxz_xxzz_1[j] + 1.5 * fl1_fx * tg_xxz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxzz_1[j] + fl1_fxn * tg_xxxz_xzz_1[j];

                    tg_xxxxz_xyyy_0[j] = pb_x * tg_xxxz_xyyy_0[j] + wp_x[j] * tg_xxxz_xyyy_1[j] + 1.5 * fl1_fx * tg_xxz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxxz_yyy_1[j];

                    tg_xxxxz_xyyz_0[j] = pb_x * tg_xxxz_xyyz_0[j] + wp_x[j] * tg_xxxz_xyyz_1[j] + 1.5 * fl1_fx * tg_xxz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yyz_1[j];

                    tg_xxxxz_xyzz_0[j] = pb_x * tg_xxxz_xyzz_0[j] + wp_x[j] * tg_xxxz_xyzz_1[j] + 1.5 * fl1_fx * tg_xxz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yzz_1[j];

                    tg_xxxxz_xzzz_0[j] = pb_x * tg_xxxz_xzzz_0[j] + wp_x[j] * tg_xxxz_xzzz_1[j] + 1.5 * fl1_fx * tg_xxz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_zzz_1[j];

                    tg_xxxxz_yyyy_0[j] = pb_x * tg_xxxz_yyyy_0[j] + wp_x[j] * tg_xxxz_yyyy_1[j] + 1.5 * fl1_fx * tg_xxz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyyy_1[j];

                    tg_xxxxz_yyyz_0[j] = pb_x * tg_xxxz_yyyz_0[j] + wp_x[j] * tg_xxxz_yyyz_1[j] + 1.5 * fl1_fx * tg_xxz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyyz_1[j];

                    tg_xxxxz_yyzz_0[j] = pb_x * tg_xxxz_yyzz_0[j] + wp_x[j] * tg_xxxz_yyzz_1[j] + 1.5 * fl1_fx * tg_xxz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyzz_1[j];

                    tg_xxxxz_yzzz_0[j] = pb_x * tg_xxxz_yzzz_0[j] + wp_x[j] * tg_xxxz_yzzz_1[j] + 1.5 * fl1_fx * tg_xxz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yzzz_1[j];

                    tg_xxxxz_zzzz_0[j] = pb_x * tg_xxxz_zzzz_0[j] + wp_x[j] * tg_xxxz_zzzz_1[j] + 1.5 * fl1_fx * tg_xxz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_zzzz_1[j];

                    tg_xxxyy_xxxx_0[j] = pb_x * tg_xxyy_xxxx_0[j] + wp_x[j] * tg_xxyy_xxxx_1[j] + fl1_fx * tg_xyy_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyy_xxx_1[j];

                    tg_xxxyy_xxxy_0[j] = pb_x * tg_xxyy_xxxy_0[j] + wp_x[j] * tg_xxyy_xxxy_1[j] + fl1_fx * tg_xyy_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyy_xxy_1[j];

                    tg_xxxyy_xxxz_0[j] = pb_x * tg_xxyy_xxxz_0[j] + wp_x[j] * tg_xxyy_xxxz_1[j] + fl1_fx * tg_xyy_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyy_xxz_1[j];

                    tg_xxxyy_xxyy_0[j] = pb_x * tg_xxyy_xxyy_0[j] + wp_x[j] * tg_xxyy_xxyy_1[j] + fl1_fx * tg_xyy_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxyy_1[j] + fl1_fxn * tg_xxyy_xyy_1[j];

                    tg_xxxyy_xxyz_0[j] = pb_x * tg_xxyy_xxyz_0[j] + wp_x[j] * tg_xxyy_xxyz_1[j] + fl1_fx * tg_xyy_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxyz_1[j] + fl1_fxn * tg_xxyy_xyz_1[j];

                    tg_xxxyy_xxzz_0[j] = pb_x * tg_xxyy_xxzz_0[j] + wp_x[j] * tg_xxyy_xxzz_1[j] + fl1_fx * tg_xyy_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxzz_1[j] + fl1_fxn * tg_xxyy_xzz_1[j];

                    tg_xxxyy_xyyy_0[j] = pb_x * tg_xxyy_xyyy_0[j] + wp_x[j] * tg_xxyy_xyyy_1[j] + fl1_fx * tg_xyy_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyy_yyy_1[j];

                    tg_xxxyy_xyyz_0[j] = pb_x * tg_xxyy_xyyz_0[j] + wp_x[j] * tg_xxyy_xyyz_1[j] + fl1_fx * tg_xyy_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yyz_1[j];

                    tg_xxxyy_xyzz_0[j] = pb_x * tg_xxyy_xyzz_0[j] + wp_x[j] * tg_xxyy_xyzz_1[j] + fl1_fx * tg_xyy_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yzz_1[j];

                    tg_xxxyy_xzzz_0[j] = pb_x * tg_xxyy_xzzz_0[j] + wp_x[j] * tg_xxyy_xzzz_1[j] + fl1_fx * tg_xyy_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_zzz_1[j];

                    tg_xxxyy_yyyy_0[j] = pb_x * tg_xxyy_yyyy_0[j] + wp_x[j] * tg_xxyy_yyyy_1[j] + fl1_fx * tg_xyy_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyy_yyyy_1[j];

                    tg_xxxyy_yyyz_0[j] = pb_x * tg_xxyy_yyyz_0[j] + wp_x[j] * tg_xxyy_yyyz_1[j] + fl1_fx * tg_xyy_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyyz_1[j];

                    tg_xxxyy_yyzz_0[j] = pb_x * tg_xxyy_yyzz_0[j] + wp_x[j] * tg_xxyy_yyzz_1[j] + fl1_fx * tg_xyy_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyzz_1[j];

                    tg_xxxyy_yzzz_0[j] = pb_x * tg_xxyy_yzzz_0[j] + wp_x[j] * tg_xxyy_yzzz_1[j] + fl1_fx * tg_xyy_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yzzz_1[j];

                    tg_xxxyy_zzzz_0[j] = pb_x * tg_xxyy_zzzz_0[j] + wp_x[j] * tg_xxyy_zzzz_1[j] + fl1_fx * tg_xyy_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyy_zzzz_1[j];

                    tg_xxxyz_xxxx_0[j] = pb_x * tg_xxyz_xxxx_0[j] + wp_x[j] * tg_xxyz_xxxx_1[j] + fl1_fx * tg_xyz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxyz_xxx_1[j];

                    tg_xxxyz_xxxy_0[j] = pb_x * tg_xxyz_xxxy_0[j] + wp_x[j] * tg_xxyz_xxxy_1[j] + fl1_fx * tg_xyz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxyz_xxy_1[j];

                    tg_xxxyz_xxxz_0[j] = pb_x * tg_xxyz_xxxz_0[j] + wp_x[j] * tg_xxyz_xxxz_1[j] + fl1_fx * tg_xyz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxyz_xxz_1[j];

                    tg_xxxyz_xxyy_0[j] = pb_x * tg_xxyz_xxyy_0[j] + wp_x[j] * tg_xxyz_xxyy_1[j] + fl1_fx * tg_xyz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxyy_1[j] + fl1_fxn * tg_xxyz_xyy_1[j];

                    tg_xxxyz_xxyz_0[j] = pb_x * tg_xxyz_xxyz_0[j] + wp_x[j] * tg_xxyz_xxyz_1[j] + fl1_fx * tg_xyz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxyz_1[j] + fl1_fxn * tg_xxyz_xyz_1[j];

                    tg_xxxyz_xxzz_0[j] = pb_x * tg_xxyz_xxzz_0[j] + wp_x[j] * tg_xxyz_xxzz_1[j] + fl1_fx * tg_xyz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxzz_1[j] + fl1_fxn * tg_xxyz_xzz_1[j];

                    tg_xxxyz_xyyy_0[j] = pb_x * tg_xxyz_xyyy_0[j] + wp_x[j] * tg_xxyz_xyyy_1[j] + fl1_fx * tg_xyz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxyz_yyy_1[j];

                    tg_xxxyz_xyyz_0[j] = pb_x * tg_xxyz_xyyz_0[j] + wp_x[j] * tg_xxyz_xyyz_1[j] + fl1_fx * tg_xyz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yyz_1[j];

                    tg_xxxyz_xyzz_0[j] = pb_x * tg_xxyz_xyzz_0[j] + wp_x[j] * tg_xxyz_xyzz_1[j] + fl1_fx * tg_xyz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yzz_1[j];

                    tg_xxxyz_xzzz_0[j] = pb_x * tg_xxyz_xzzz_0[j] + wp_x[j] * tg_xxyz_xzzz_1[j] + fl1_fx * tg_xyz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_zzz_1[j];

                    tg_xxxyz_yyyy_0[j] = pb_x * tg_xxyz_yyyy_0[j] + wp_x[j] * tg_xxyz_yyyy_1[j] + fl1_fx * tg_xyz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xyz_yyyy_1[j];

                    tg_xxxyz_yyyz_0[j] = pb_x * tg_xxyz_yyyz_0[j] + wp_x[j] * tg_xxyz_yyyz_1[j] + fl1_fx * tg_xyz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyyz_1[j];

                    tg_xxxyz_yyzz_0[j] = pb_x * tg_xxyz_yyzz_0[j] + wp_x[j] * tg_xxyz_yyzz_1[j] + fl1_fx * tg_xyz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyzz_1[j];

                    tg_xxxyz_yzzz_0[j] = pb_x * tg_xxyz_yzzz_0[j] + wp_x[j] * tg_xxyz_yzzz_1[j] + fl1_fx * tg_xyz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yzzz_1[j];

                    tg_xxxyz_zzzz_0[j] = pb_x * tg_xxyz_zzzz_0[j] + wp_x[j] * tg_xxyz_zzzz_1[j] + fl1_fx * tg_xyz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xyz_zzzz_1[j];

                    tg_xxxzz_xxxx_0[j] = pb_x * tg_xxzz_xxxx_0[j] + wp_x[j] * tg_xxzz_xxxx_1[j] + fl1_fx * tg_xzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxzz_xxx_1[j];

                    tg_xxxzz_xxxy_0[j] = pb_x * tg_xxzz_xxxy_0[j] + wp_x[j] * tg_xxzz_xxxy_1[j] + fl1_fx * tg_xzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxzz_xxy_1[j];

                    tg_xxxzz_xxxz_0[j] = pb_x * tg_xxzz_xxxz_0[j] + wp_x[j] * tg_xxzz_xxxz_1[j] + fl1_fx * tg_xzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxzz_xxz_1[j];

                    tg_xxxzz_xxyy_0[j] = pb_x * tg_xxzz_xxyy_0[j] + wp_x[j] * tg_xxzz_xxyy_1[j] + fl1_fx * tg_xzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxyy_1[j] + fl1_fxn * tg_xxzz_xyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSG_79_158(      CMemBlock2D<double>& primBuffer,
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
                                             {5, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xxzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 79); 

                auto tg_xxzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 80); 

                auto tg_xxzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 81); 

                auto tg_xxzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 82); 

                auto tg_xxzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 83); 

                auto tg_xxzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 84); 

                auto tg_xxzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 85); 

                auto tg_xxzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 86); 

                auto tg_xxzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 87); 

                auto tg_xxzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 88); 

                auto tg_xxzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 89); 

                auto tg_xyyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 90); 

                auto tg_xyyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 91); 

                auto tg_xyyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 92); 

                auto tg_xyyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 93); 

                auto tg_xyyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 94); 

                auto tg_xyyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 95); 

                auto tg_xyyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 96); 

                auto tg_xyyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 97); 

                auto tg_xyyy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 98); 

                auto tg_xyyy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 99); 

                auto tg_xyyy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 100); 

                auto tg_xyyy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 101); 

                auto tg_xyyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 102); 

                auto tg_xyyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 103); 

                auto tg_xyyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 104); 

                auto tg_xyyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 105); 

                auto tg_xyyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 106); 

                auto tg_xyyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 107); 

                auto tg_xyyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 108); 

                auto tg_xyyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 109); 

                auto tg_xyyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 110); 

                auto tg_xyyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 111); 

                auto tg_xyyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 112); 

                auto tg_xyyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 113); 

                auto tg_xyyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 114); 

                auto tg_xyyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 115); 

                auto tg_xyyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 116); 

                auto tg_xyyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 117); 

                auto tg_xyyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 118); 

                auto tg_xyyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 119); 

                auto tg_xyzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 120); 

                auto tg_xyzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 121); 

                auto tg_xyzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 122); 

                auto tg_xyzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 123); 

                auto tg_xyzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 124); 

                auto tg_xyzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 125); 

                auto tg_xyzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 126); 

                auto tg_xyzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 127); 

                auto tg_xyzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 128); 

                auto tg_xyzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 129); 

                auto tg_xyzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 130); 

                auto tg_xyzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 131); 

                auto tg_xyzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 132); 

                auto tg_xyzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 133); 

                auto tg_xyzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 134); 

                auto tg_xzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 135); 

                auto tg_xzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 136); 

                auto tg_xzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 137); 

                auto tg_xzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 138); 

                auto tg_xzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 139); 

                auto tg_xzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 140); 

                auto tg_xzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 141); 

                auto tg_xzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 142); 

                auto tg_xzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 143); 

                auto tg_xzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 144); 

                auto tg_xzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 145); 

                auto tg_xzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 146); 

                auto tg_xzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 147); 

                auto tg_xzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 148); 

                auto tg_xzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 149); 

                auto tg_yyyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 157); 

                auto tg_xxzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 79); 

                auto tg_xxzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 80); 

                auto tg_xxzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 81); 

                auto tg_xxzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 82); 

                auto tg_xxzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 83); 

                auto tg_xxzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 84); 

                auto tg_xxzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 85); 

                auto tg_xxzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 86); 

                auto tg_xxzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 87); 

                auto tg_xxzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 88); 

                auto tg_xxzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 89); 

                auto tg_xyyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 90); 

                auto tg_xyyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 91); 

                auto tg_xyyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 92); 

                auto tg_xyyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 93); 

                auto tg_xyyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 94); 

                auto tg_xyyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 95); 

                auto tg_xyyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 96); 

                auto tg_xyyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 97); 

                auto tg_xyyy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 98); 

                auto tg_xyyy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 99); 

                auto tg_xyyy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 100); 

                auto tg_xyyy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 101); 

                auto tg_xyyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 102); 

                auto tg_xyyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 103); 

                auto tg_xyyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 104); 

                auto tg_xyyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 105); 

                auto tg_xyyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 106); 

                auto tg_xyyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 107); 

                auto tg_xyyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 108); 

                auto tg_xyyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 109); 

                auto tg_xyyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 110); 

                auto tg_xyyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 111); 

                auto tg_xyyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 112); 

                auto tg_xyyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 113); 

                auto tg_xyyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 114); 

                auto tg_xyyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 115); 

                auto tg_xyyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 116); 

                auto tg_xyyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 117); 

                auto tg_xyyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 118); 

                auto tg_xyyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 119); 

                auto tg_xyzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 120); 

                auto tg_xyzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 121); 

                auto tg_xyzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 122); 

                auto tg_xyzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 123); 

                auto tg_xyzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 124); 

                auto tg_xyzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 125); 

                auto tg_xyzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 126); 

                auto tg_xyzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 127); 

                auto tg_xyzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 128); 

                auto tg_xyzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 129); 

                auto tg_xyzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 130); 

                auto tg_xyzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 131); 

                auto tg_xyzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 132); 

                auto tg_xyzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 133); 

                auto tg_xyzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 134); 

                auto tg_xzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 135); 

                auto tg_xzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 136); 

                auto tg_xzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 137); 

                auto tg_xzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 138); 

                auto tg_xzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 139); 

                auto tg_xzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 140); 

                auto tg_xzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 141); 

                auto tg_xzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 142); 

                auto tg_xzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 143); 

                auto tg_xzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 144); 

                auto tg_xzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 145); 

                auto tg_xzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 146); 

                auto tg_xzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 147); 

                auto tg_xzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 148); 

                auto tg_xzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 149); 

                auto tg_yyyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 157); 

                auto tg_xzz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 79); 

                auto tg_xzz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 80); 

                auto tg_xzz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 81); 

                auto tg_xzz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 82); 

                auto tg_xzz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 83); 

                auto tg_xzz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 84); 

                auto tg_xzz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 85); 

                auto tg_xzz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 86); 

                auto tg_xzz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 87); 

                auto tg_xzz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 88); 

                auto tg_xzz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 89); 

                auto tg_yyy_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 90); 

                auto tg_yyy_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 91); 

                auto tg_yyy_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 92); 

                auto tg_yyy_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 93); 

                auto tg_yyy_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 94); 

                auto tg_yyy_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 95); 

                auto tg_yyy_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 96); 

                auto tg_yyy_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 97); 

                auto tg_yyy_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 98); 

                auto tg_yyy_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 99); 

                auto tg_yyy_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 100); 

                auto tg_yyy_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 101); 

                auto tg_yyy_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 102); 

                auto tg_yyy_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 103); 

                auto tg_yyy_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 104); 

                auto tg_yyz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 105); 

                auto tg_yyz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 106); 

                auto tg_yyz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 107); 

                auto tg_yyz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 108); 

                auto tg_yyz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 109); 

                auto tg_yyz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 110); 

                auto tg_yyz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 111); 

                auto tg_yyz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 112); 

                auto tg_yyz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 113); 

                auto tg_yyz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 114); 

                auto tg_yyz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 115); 

                auto tg_yyz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 116); 

                auto tg_yyz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 117); 

                auto tg_yyz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 118); 

                auto tg_yyz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 119); 

                auto tg_yzz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 120); 

                auto tg_yzz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 121); 

                auto tg_yzz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 122); 

                auto tg_yzz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 123); 

                auto tg_yzz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 124); 

                auto tg_yzz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 125); 

                auto tg_yzz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 126); 

                auto tg_yzz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 127); 

                auto tg_yzz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 128); 

                auto tg_yzz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 129); 

                auto tg_yzz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 130); 

                auto tg_yzz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 131); 

                auto tg_yzz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 132); 

                auto tg_yzz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 133); 

                auto tg_yzz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 134); 

                auto tg_zzz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 135); 

                auto tg_zzz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 136); 

                auto tg_zzz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 137); 

                auto tg_zzz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 138); 

                auto tg_zzz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 139); 

                auto tg_zzz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 140); 

                auto tg_zzz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 141); 

                auto tg_zzz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 142); 

                auto tg_zzz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 143); 

                auto tg_zzz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 144); 

                auto tg_zzz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 145); 

                auto tg_zzz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 146); 

                auto tg_zzz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 147); 

                auto tg_zzz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 148); 

                auto tg_zzz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 149); 

                auto tg_xzz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 79); 

                auto tg_xzz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 80); 

                auto tg_xzz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 81); 

                auto tg_xzz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 82); 

                auto tg_xzz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 83); 

                auto tg_xzz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 84); 

                auto tg_xzz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 85); 

                auto tg_xzz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 86); 

                auto tg_xzz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 87); 

                auto tg_xzz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 88); 

                auto tg_xzz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 89); 

                auto tg_yyy_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 90); 

                auto tg_yyy_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 91); 

                auto tg_yyy_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 92); 

                auto tg_yyy_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 93); 

                auto tg_yyy_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 94); 

                auto tg_yyy_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 95); 

                auto tg_yyy_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 96); 

                auto tg_yyy_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 97); 

                auto tg_yyy_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 98); 

                auto tg_yyy_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 99); 

                auto tg_yyy_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 100); 

                auto tg_yyy_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 101); 

                auto tg_yyy_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 102); 

                auto tg_yyy_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 103); 

                auto tg_yyy_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 104); 

                auto tg_yyz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 105); 

                auto tg_yyz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 106); 

                auto tg_yyz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 107); 

                auto tg_yyz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 108); 

                auto tg_yyz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 109); 

                auto tg_yyz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 110); 

                auto tg_yyz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 111); 

                auto tg_yyz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 112); 

                auto tg_yyz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 113); 

                auto tg_yyz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 114); 

                auto tg_yyz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 115); 

                auto tg_yyz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 116); 

                auto tg_yyz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 117); 

                auto tg_yyz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 118); 

                auto tg_yyz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 119); 

                auto tg_yzz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 120); 

                auto tg_yzz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 121); 

                auto tg_yzz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 122); 

                auto tg_yzz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 123); 

                auto tg_yzz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 124); 

                auto tg_yzz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 125); 

                auto tg_yzz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 126); 

                auto tg_yzz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 127); 

                auto tg_yzz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 128); 

                auto tg_yzz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 129); 

                auto tg_yzz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 130); 

                auto tg_yzz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 131); 

                auto tg_yzz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 132); 

                auto tg_yzz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 133); 

                auto tg_yzz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 134); 

                auto tg_zzz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 135); 

                auto tg_zzz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 136); 

                auto tg_zzz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 137); 

                auto tg_zzz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 138); 

                auto tg_zzz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 139); 

                auto tg_zzz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 140); 

                auto tg_zzz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 141); 

                auto tg_zzz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 142); 

                auto tg_zzz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 143); 

                auto tg_zzz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 144); 

                auto tg_zzz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 145); 

                auto tg_zzz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 146); 

                auto tg_zzz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 147); 

                auto tg_zzz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 148); 

                auto tg_zzz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 149); 

                auto tg_xxzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 54); 

                auto tg_xxzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 55); 

                auto tg_xxzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 56); 

                auto tg_xxzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 57); 

                auto tg_xxzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 58); 

                auto tg_xxzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 59); 

                auto tg_xyyy_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 60); 

                auto tg_xyyy_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 61); 

                auto tg_xyyy_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 62); 

                auto tg_xyyy_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 63); 

                auto tg_xyyy_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 64); 

                auto tg_xyyy_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 65); 

                auto tg_xyyy_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 66); 

                auto tg_xyyy_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 67); 

                auto tg_xyyy_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 68); 

                auto tg_xyyy_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 69); 

                auto tg_xyyz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 70); 

                auto tg_xyyz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 71); 

                auto tg_xyyz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 72); 

                auto tg_xyyz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 73); 

                auto tg_xyyz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 74); 

                auto tg_xyyz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 75); 

                auto tg_xyyz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 76); 

                auto tg_xyyz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 77); 

                auto tg_xyyz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 78); 

                auto tg_xyyz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 79); 

                auto tg_xyzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 80); 

                auto tg_xyzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 81); 

                auto tg_xyzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 82); 

                auto tg_xyzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 83); 

                auto tg_xyzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 84); 

                auto tg_xyzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 85); 

                auto tg_xyzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 86); 

                auto tg_xyzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 87); 

                auto tg_xyzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 88); 

                auto tg_xyzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 89); 

                auto tg_xzzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 90); 

                auto tg_xzzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 91); 

                auto tg_xzzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 92); 

                auto tg_xzzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 93); 

                auto tg_xzzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 94); 

                auto tg_xzzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 95); 

                auto tg_xzzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 96); 

                auto tg_xzzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 97); 

                auto tg_xzzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 98); 

                auto tg_xzzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 99); 

                auto tg_yyyy_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 100); 

                auto tg_yyyy_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 101); 

                auto tg_yyyy_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 102); 

                auto tg_yyyy_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 103); 

                auto tg_yyyy_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 104); 

                auto tg_yyyy_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 105); 

                auto tg_yyyy_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 106); 

                auto tg_yyyy_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 107); 

                // set up pointers to integrals

                auto tg_xxxzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 79); 

                auto tg_xxxzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 80); 

                auto tg_xxxzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 81); 

                auto tg_xxxzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 82); 

                auto tg_xxxzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 83); 

                auto tg_xxxzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 84); 

                auto tg_xxxzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 85); 

                auto tg_xxxzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 86); 

                auto tg_xxxzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 87); 

                auto tg_xxxzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 88); 

                auto tg_xxxzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 89); 

                auto tg_xxyyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 90); 

                auto tg_xxyyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 91); 

                auto tg_xxyyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 92); 

                auto tg_xxyyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 93); 

                auto tg_xxyyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 94); 

                auto tg_xxyyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 95); 

                auto tg_xxyyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 96); 

                auto tg_xxyyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 97); 

                auto tg_xxyyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 98); 

                auto tg_xxyyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 99); 

                auto tg_xxyyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 100); 

                auto tg_xxyyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 101); 

                auto tg_xxyyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 102); 

                auto tg_xxyyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 103); 

                auto tg_xxyyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 104); 

                auto tg_xxyyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 105); 

                auto tg_xxyyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 106); 

                auto tg_xxyyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 107); 

                auto tg_xxyyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 108); 

                auto tg_xxyyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 109); 

                auto tg_xxyyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 110); 

                auto tg_xxyyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 111); 

                auto tg_xxyyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 112); 

                auto tg_xxyyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 113); 

                auto tg_xxyyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 114); 

                auto tg_xxyyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 115); 

                auto tg_xxyyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 116); 

                auto tg_xxyyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 117); 

                auto tg_xxyyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 118); 

                auto tg_xxyyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 119); 

                auto tg_xxyzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 120); 

                auto tg_xxyzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 121); 

                auto tg_xxyzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 122); 

                auto tg_xxyzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 123); 

                auto tg_xxyzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 124); 

                auto tg_xxyzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 125); 

                auto tg_xxyzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 126); 

                auto tg_xxyzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 127); 

                auto tg_xxyzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 128); 

                auto tg_xxyzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 129); 

                auto tg_xxyzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 130); 

                auto tg_xxyzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 131); 

                auto tg_xxyzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 132); 

                auto tg_xxyzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 133); 

                auto tg_xxyzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 134); 

                auto tg_xxzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 135); 

                auto tg_xxzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 136); 

                auto tg_xxzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 137); 

                auto tg_xxzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 138); 

                auto tg_xxzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 139); 

                auto tg_xxzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 140); 

                auto tg_xxzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 141); 

                auto tg_xxzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 142); 

                auto tg_xxzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 143); 

                auto tg_xxzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 144); 

                auto tg_xxzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 145); 

                auto tg_xxzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 146); 

                auto tg_xxzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 147); 

                auto tg_xxzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 148); 

                auto tg_xxzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 149); 

                auto tg_xyyyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 150); 

                auto tg_xyyyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 151); 

                auto tg_xyyyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 152); 

                auto tg_xyyyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 153); 

                auto tg_xyyyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 154); 

                auto tg_xyyyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 155); 

                auto tg_xyyyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 156); 

                auto tg_xyyyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 157); 

                // Batch of Integrals (79,158)

                #pragma omp simd aligned(fxn, fza, tg_xxxzz_xxyz_0, tg_xxxzz_xxzz_0, tg_xxxzz_xyyy_0, \
                                         tg_xxxzz_xyyz_0, tg_xxxzz_xyzz_0, tg_xxxzz_xzzz_0, tg_xxxzz_yyyy_0, tg_xxxzz_yyyz_0, \
                                         tg_xxxzz_yyzz_0, tg_xxxzz_yzzz_0, tg_xxxzz_zzzz_0, tg_xxyyy_xxxx_0, tg_xxyyy_xxxy_0, \
                                         tg_xxyyy_xxxz_0, tg_xxyyy_xxyy_0, tg_xxyyy_xxyz_0, tg_xxyyy_xxzz_0, tg_xxyyy_xyyy_0, \
                                         tg_xxyyy_xyyz_0, tg_xxyyy_xyzz_0, tg_xxyyy_xzzz_0, tg_xxyyy_yyyy_0, tg_xxyyy_yyyz_0, \
                                         tg_xxyyy_yyzz_0, tg_xxyyy_yzzz_0, tg_xxyyy_zzzz_0, tg_xxyyz_xxxx_0, tg_xxyyz_xxxy_0, \
                                         tg_xxyyz_xxxz_0, tg_xxyyz_xxyy_0, tg_xxyyz_xxyz_0, tg_xxyyz_xxzz_0, tg_xxyyz_xyyy_0, \
                                         tg_xxyyz_xyyz_0, tg_xxyyz_xyzz_0, tg_xxyyz_xzzz_0, tg_xxyyz_yyyy_0, tg_xxyyz_yyyz_0, \
                                         tg_xxyyz_yyzz_0, tg_xxyyz_yzzz_0, tg_xxyyz_zzzz_0, tg_xxyzz_xxxx_0, tg_xxyzz_xxxy_0, \
                                         tg_xxyzz_xxxz_0, tg_xxyzz_xxyy_0, tg_xxyzz_xxyz_0, tg_xxyzz_xxzz_0, tg_xxyzz_xyyy_0, \
                                         tg_xxyzz_xyyz_0, tg_xxyzz_xyzz_0, tg_xxyzz_xzzz_0, tg_xxyzz_yyyy_0, tg_xxyzz_yyyz_0, \
                                         tg_xxyzz_yyzz_0, tg_xxyzz_yzzz_0, tg_xxyzz_zzzz_0, tg_xxzz_xxyz_0, tg_xxzz_xxyz_1, \
                                         tg_xxzz_xxzz_0, tg_xxzz_xxzz_1, tg_xxzz_xyyy_0, tg_xxzz_xyyy_1, tg_xxzz_xyyz_0, \
                                         tg_xxzz_xyyz_1, tg_xxzz_xyz_1, tg_xxzz_xyzz_0, tg_xxzz_xyzz_1, tg_xxzz_xzz_1, \
                                         tg_xxzz_xzzz_0, tg_xxzz_xzzz_1, tg_xxzz_yyy_1, tg_xxzz_yyyy_0, tg_xxzz_yyyy_1, \
                                         tg_xxzz_yyyz_0, tg_xxzz_yyyz_1, tg_xxzz_yyz_1, tg_xxzz_yyzz_0, tg_xxzz_yyzz_1, \
                                         tg_xxzz_yzz_1, tg_xxzz_yzzz_0, tg_xxzz_yzzz_1, tg_xxzz_zzz_1, tg_xxzz_zzzz_0, \
                                         tg_xxzz_zzzz_1, tg_xxzzz_xxxx_0, tg_xxzzz_xxxy_0, tg_xxzzz_xxxz_0, tg_xxzzz_xxyy_0, \
                                         tg_xxzzz_xxyz_0, tg_xxzzz_xxzz_0, tg_xxzzz_xyyy_0, tg_xxzzz_xyyz_0, tg_xxzzz_xyzz_0, \
                                         tg_xxzzz_xzzz_0, tg_xxzzz_yyyy_0, tg_xxzzz_yyyz_0, tg_xxzzz_yyzz_0, tg_xxzzz_yzzz_0, \
                                         tg_xxzzz_zzzz_0, tg_xyyy_xxx_1, tg_xyyy_xxxx_0, tg_xyyy_xxxx_1, tg_xyyy_xxxy_0, \
                                         tg_xyyy_xxxy_1, tg_xyyy_xxxz_0, tg_xyyy_xxxz_1, tg_xyyy_xxy_1, tg_xyyy_xxyy_0, \
                                         tg_xyyy_xxyy_1, tg_xyyy_xxyz_0, tg_xyyy_xxyz_1, tg_xyyy_xxz_1, tg_xyyy_xxzz_0, \
                                         tg_xyyy_xxzz_1, tg_xyyy_xyy_1, tg_xyyy_xyyy_0, tg_xyyy_xyyy_1, tg_xyyy_xyyz_0, \
                                         tg_xyyy_xyyz_1, tg_xyyy_xyz_1, tg_xyyy_xyzz_0, tg_xyyy_xyzz_1, tg_xyyy_xzz_1, \
                                         tg_xyyy_xzzz_0, tg_xyyy_xzzz_1, tg_xyyy_yyy_1, tg_xyyy_yyyy_0, tg_xyyy_yyyy_1, \
                                         tg_xyyy_yyyz_0, tg_xyyy_yyyz_1, tg_xyyy_yyz_1, tg_xyyy_yyzz_0, tg_xyyy_yyzz_1, \
                                         tg_xyyy_yzz_1, tg_xyyy_yzzz_0, tg_xyyy_yzzz_1, tg_xyyy_zzz_1, tg_xyyy_zzzz_0, \
                                         tg_xyyy_zzzz_1, tg_xyyyy_xxxx_0, tg_xyyyy_xxxy_0, tg_xyyyy_xxxz_0, tg_xyyyy_xxyy_0, \
                                         tg_xyyyy_xxyz_0, tg_xyyyy_xxzz_0, tg_xyyyy_xyyy_0, tg_xyyyy_xyyz_0, tg_xyyz_xxx_1, \
                                         tg_xyyz_xxxx_0, tg_xyyz_xxxx_1, tg_xyyz_xxxy_0, tg_xyyz_xxxy_1, tg_xyyz_xxxz_0, \
                                         tg_xyyz_xxxz_1, tg_xyyz_xxy_1, tg_xyyz_xxyy_0, tg_xyyz_xxyy_1, tg_xyyz_xxyz_0, \
                                         tg_xyyz_xxyz_1, tg_xyyz_xxz_1, tg_xyyz_xxzz_0, tg_xyyz_xxzz_1, tg_xyyz_xyy_1, \
                                         tg_xyyz_xyyy_0, tg_xyyz_xyyy_1, tg_xyyz_xyyz_0, tg_xyyz_xyyz_1, tg_xyyz_xyz_1, \
                                         tg_xyyz_xyzz_0, tg_xyyz_xyzz_1, tg_xyyz_xzz_1, tg_xyyz_xzzz_0, tg_xyyz_xzzz_1, \
                                         tg_xyyz_yyy_1, tg_xyyz_yyyy_0, tg_xyyz_yyyy_1, tg_xyyz_yyyz_0, tg_xyyz_yyyz_1, \
                                         tg_xyyz_yyz_1, tg_xyyz_yyzz_0, tg_xyyz_yyzz_1, tg_xyyz_yzz_1, tg_xyyz_yzzz_0, \
                                         tg_xyyz_yzzz_1, tg_xyyz_zzz_1, tg_xyyz_zzzz_0, tg_xyyz_zzzz_1, tg_xyzz_xxx_1, \
                                         tg_xyzz_xxxx_0, tg_xyzz_xxxx_1, tg_xyzz_xxxy_0, tg_xyzz_xxxy_1, tg_xyzz_xxxz_0, \
                                         tg_xyzz_xxxz_1, tg_xyzz_xxy_1, tg_xyzz_xxyy_0, tg_xyzz_xxyy_1, tg_xyzz_xxyz_0, \
                                         tg_xyzz_xxyz_1, tg_xyzz_xxz_1, tg_xyzz_xxzz_0, tg_xyzz_xxzz_1, tg_xyzz_xyy_1, \
                                         tg_xyzz_xyyy_0, tg_xyzz_xyyy_1, tg_xyzz_xyyz_0, tg_xyzz_xyyz_1, tg_xyzz_xyz_1, \
                                         tg_xyzz_xyzz_0, tg_xyzz_xyzz_1, tg_xyzz_xzz_1, tg_xyzz_xzzz_0, tg_xyzz_xzzz_1, \
                                         tg_xyzz_yyy_1, tg_xyzz_yyyy_0, tg_xyzz_yyyy_1, tg_xyzz_yyyz_0, tg_xyzz_yyyz_1, \
                                         tg_xyzz_yyz_1, tg_xyzz_yyzz_0, tg_xyzz_yyzz_1, tg_xyzz_yzz_1, tg_xyzz_yzzz_0, \
                                         tg_xyzz_yzzz_1, tg_xyzz_zzz_1, tg_xyzz_zzzz_0, tg_xyzz_zzzz_1, tg_xzz_xxyz_0, \
                                         tg_xzz_xxyz_1, tg_xzz_xxzz_0, tg_xzz_xxzz_1, tg_xzz_xyyy_0, tg_xzz_xyyy_1, \
                                         tg_xzz_xyyz_0, tg_xzz_xyyz_1, tg_xzz_xyzz_0, tg_xzz_xyzz_1, tg_xzz_xzzz_0, \
                                         tg_xzz_xzzz_1, tg_xzz_yyyy_0, tg_xzz_yyyy_1, tg_xzz_yyyz_0, tg_xzz_yyyz_1, \
                                         tg_xzz_yyzz_0, tg_xzz_yyzz_1, tg_xzz_yzzz_0, tg_xzz_yzzz_1, tg_xzz_zzzz_0, \
                                         tg_xzz_zzzz_1, tg_xzzz_xxx_1, tg_xzzz_xxxx_0, tg_xzzz_xxxx_1, tg_xzzz_xxxy_0, \
                                         tg_xzzz_xxxy_1, tg_xzzz_xxxz_0, tg_xzzz_xxxz_1, tg_xzzz_xxy_1, tg_xzzz_xxyy_0, \
                                         tg_xzzz_xxyy_1, tg_xzzz_xxyz_0, tg_xzzz_xxyz_1, tg_xzzz_xxz_1, tg_xzzz_xxzz_0, \
                                         tg_xzzz_xxzz_1, tg_xzzz_xyy_1, tg_xzzz_xyyy_0, tg_xzzz_xyyy_1, tg_xzzz_xyyz_0, \
                                         tg_xzzz_xyyz_1, tg_xzzz_xyz_1, tg_xzzz_xyzz_0, tg_xzzz_xyzz_1, tg_xzzz_xzz_1, \
                                         tg_xzzz_xzzz_0, tg_xzzz_xzzz_1, tg_xzzz_yyy_1, tg_xzzz_yyyy_0, tg_xzzz_yyyy_1, \
                                         tg_xzzz_yyyz_0, tg_xzzz_yyyz_1, tg_xzzz_yyz_1, tg_xzzz_yyzz_0, tg_xzzz_yyzz_1, \
                                         tg_xzzz_yzz_1, tg_xzzz_yzzz_0, tg_xzzz_yzzz_1, tg_xzzz_zzz_1, tg_xzzz_zzzz_0, \
                                         tg_xzzz_zzzz_1, tg_yyy_xxxx_0, tg_yyy_xxxx_1, tg_yyy_xxxy_0, tg_yyy_xxxy_1, \
                                         tg_yyy_xxxz_0, tg_yyy_xxxz_1, tg_yyy_xxyy_0, tg_yyy_xxyy_1, tg_yyy_xxyz_0, \
                                         tg_yyy_xxyz_1, tg_yyy_xxzz_0, tg_yyy_xxzz_1, tg_yyy_xyyy_0, tg_yyy_xyyy_1, \
                                         tg_yyy_xyyz_0, tg_yyy_xyyz_1, tg_yyy_xyzz_0, tg_yyy_xyzz_1, tg_yyy_xzzz_0, \
                                         tg_yyy_xzzz_1, tg_yyy_yyyy_0, tg_yyy_yyyy_1, tg_yyy_yyyz_0, tg_yyy_yyyz_1, \
                                         tg_yyy_yyzz_0, tg_yyy_yyzz_1, tg_yyy_yzzz_0, tg_yyy_yzzz_1, tg_yyy_zzzz_0, \
                                         tg_yyy_zzzz_1, tg_yyyy_xxx_1, tg_yyyy_xxxx_0, tg_yyyy_xxxx_1, tg_yyyy_xxxy_0, \
                                         tg_yyyy_xxxy_1, tg_yyyy_xxxz_0, tg_yyyy_xxxz_1, tg_yyyy_xxy_1, tg_yyyy_xxyy_0, \
                                         tg_yyyy_xxyy_1, tg_yyyy_xxyz_0, tg_yyyy_xxyz_1, tg_yyyy_xxz_1, tg_yyyy_xxzz_0, \
                                         tg_yyyy_xxzz_1, tg_yyyy_xyy_1, tg_yyyy_xyyy_0, tg_yyyy_xyyy_1, tg_yyyy_xyyz_0, \
                                         tg_yyyy_xyyz_1, tg_yyyy_xyz_1, tg_yyyy_xzz_1, tg_yyyy_yyy_1, tg_yyyy_yyz_1, \
                                         tg_yyz_xxxx_0, tg_yyz_xxxx_1, tg_yyz_xxxy_0, tg_yyz_xxxy_1, tg_yyz_xxxz_0, \
                                         tg_yyz_xxxz_1, tg_yyz_xxyy_0, tg_yyz_xxyy_1, tg_yyz_xxyz_0, tg_yyz_xxyz_1, \
                                         tg_yyz_xxzz_0, tg_yyz_xxzz_1, tg_yyz_xyyy_0, tg_yyz_xyyy_1, tg_yyz_xyyz_0, \
                                         tg_yyz_xyyz_1, tg_yyz_xyzz_0, tg_yyz_xyzz_1, tg_yyz_xzzz_0, tg_yyz_xzzz_1, \
                                         tg_yyz_yyyy_0, tg_yyz_yyyy_1, tg_yyz_yyyz_0, tg_yyz_yyyz_1, tg_yyz_yyzz_0, \
                                         tg_yyz_yyzz_1, tg_yyz_yzzz_0, tg_yyz_yzzz_1, tg_yyz_zzzz_0, tg_yyz_zzzz_1, \
                                         tg_yzz_xxxx_0, tg_yzz_xxxx_1, tg_yzz_xxxy_0, tg_yzz_xxxy_1, tg_yzz_xxxz_0, \
                                         tg_yzz_xxxz_1, tg_yzz_xxyy_0, tg_yzz_xxyy_1, tg_yzz_xxyz_0, tg_yzz_xxyz_1, \
                                         tg_yzz_xxzz_0, tg_yzz_xxzz_1, tg_yzz_xyyy_0, tg_yzz_xyyy_1, tg_yzz_xyyz_0, \
                                         tg_yzz_xyyz_1, tg_yzz_xyzz_0, tg_yzz_xyzz_1, tg_yzz_xzzz_0, tg_yzz_xzzz_1, \
                                         tg_yzz_yyyy_0, tg_yzz_yyyy_1, tg_yzz_yyyz_0, tg_yzz_yyyz_1, tg_yzz_yyzz_0, \
                                         tg_yzz_yyzz_1, tg_yzz_yzzz_0, tg_yzz_yzzz_1, tg_yzz_zzzz_0, tg_yzz_zzzz_1, \
                                         tg_zzz_xxxx_0, tg_zzz_xxxx_1, tg_zzz_xxxy_0, tg_zzz_xxxy_1, tg_zzz_xxxz_0, \
                                         tg_zzz_xxxz_1, tg_zzz_xxyy_0, tg_zzz_xxyy_1, tg_zzz_xxyz_0, tg_zzz_xxyz_1, \
                                         tg_zzz_xxzz_0, tg_zzz_xxzz_1, tg_zzz_xyyy_0, tg_zzz_xyyy_1, tg_zzz_xyyz_0, \
                                         tg_zzz_xyyz_1, tg_zzz_xyzz_0, tg_zzz_xyzz_1, tg_zzz_xzzz_0, tg_zzz_xzzz_1, \
                                         tg_zzz_yyyy_0, tg_zzz_yyyy_1, tg_zzz_yyyz_0, tg_zzz_yyyz_1, tg_zzz_yyzz_0, \
                                         tg_zzz_yyzz_1, tg_zzz_yzzz_0, tg_zzz_yzzz_1, tg_zzz_zzzz_0, tg_zzz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxzz_xxyz_0[j] = pb_x * tg_xxzz_xxyz_0[j] + wp_x[j] * tg_xxzz_xxyz_1[j] + fl1_fx * tg_xzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxyz_1[j] + fl1_fxn * tg_xxzz_xyz_1[j];

                    tg_xxxzz_xxzz_0[j] = pb_x * tg_xxzz_xxzz_0[j] + wp_x[j] * tg_xxzz_xxzz_1[j] + fl1_fx * tg_xzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxzz_1[j] + fl1_fxn * tg_xxzz_xzz_1[j];

                    tg_xxxzz_xyyy_0[j] = pb_x * tg_xxzz_xyyy_0[j] + wp_x[j] * tg_xxzz_xyyy_1[j] + fl1_fx * tg_xzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxzz_yyy_1[j];

                    tg_xxxzz_xyyz_0[j] = pb_x * tg_xxzz_xyyz_0[j] + wp_x[j] * tg_xxzz_xyyz_1[j] + fl1_fx * tg_xzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yyz_1[j];

                    tg_xxxzz_xyzz_0[j] = pb_x * tg_xxzz_xyzz_0[j] + wp_x[j] * tg_xxzz_xyzz_1[j] + fl1_fx * tg_xzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yzz_1[j];

                    tg_xxxzz_xzzz_0[j] = pb_x * tg_xxzz_xzzz_0[j] + wp_x[j] * tg_xxzz_xzzz_1[j] + fl1_fx * tg_xzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_zzz_1[j];

                    tg_xxxzz_yyyy_0[j] = pb_x * tg_xxzz_yyyy_0[j] + wp_x[j] * tg_xxzz_yyyy_1[j] + fl1_fx * tg_xzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xzz_yyyy_1[j];

                    tg_xxxzz_yyyz_0[j] = pb_x * tg_xxzz_yyyz_0[j] + wp_x[j] * tg_xxzz_yyyz_1[j] + fl1_fx * tg_xzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyyz_1[j];

                    tg_xxxzz_yyzz_0[j] = pb_x * tg_xxzz_yyzz_0[j] + wp_x[j] * tg_xxzz_yyzz_1[j] + fl1_fx * tg_xzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyzz_1[j];

                    tg_xxxzz_yzzz_0[j] = pb_x * tg_xxzz_yzzz_0[j] + wp_x[j] * tg_xxzz_yzzz_1[j] + fl1_fx * tg_xzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yzzz_1[j];

                    tg_xxxzz_zzzz_0[j] = pb_x * tg_xxzz_zzzz_0[j] + wp_x[j] * tg_xxzz_zzzz_1[j] + fl1_fx * tg_xzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xzz_zzzz_1[j];

                    tg_xxyyy_xxxx_0[j] = pb_x * tg_xyyy_xxxx_0[j] + wp_x[j] * tg_xyyy_xxxx_1[j] + 0.5 * fl1_fx * tg_yyy_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyy_xxx_1[j];

                    tg_xxyyy_xxxy_0[j] = pb_x * tg_xyyy_xxxy_0[j] + wp_x[j] * tg_xyyy_xxxy_1[j] + 0.5 * fl1_fx * tg_yyy_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyy_xxy_1[j];

                    tg_xxyyy_xxxz_0[j] = pb_x * tg_xyyy_xxxz_0[j] + wp_x[j] * tg_xyyy_xxxz_1[j] + 0.5 * fl1_fx * tg_yyy_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyy_xxz_1[j];

                    tg_xxyyy_xxyy_0[j] = pb_x * tg_xyyy_xxyy_0[j] + wp_x[j] * tg_xyyy_xxyy_1[j] + 0.5 * fl1_fx * tg_yyy_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxyy_1[j] + fl1_fxn * tg_xyyy_xyy_1[j];

                    tg_xxyyy_xxyz_0[j] = pb_x * tg_xyyy_xxyz_0[j] + wp_x[j] * tg_xyyy_xxyz_1[j] + 0.5 * fl1_fx * tg_yyy_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxyz_1[j] + fl1_fxn * tg_xyyy_xyz_1[j];

                    tg_xxyyy_xxzz_0[j] = pb_x * tg_xyyy_xxzz_0[j] + wp_x[j] * tg_xyyy_xxzz_1[j] + 0.5 * fl1_fx * tg_yyy_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxzz_1[j] + fl1_fxn * tg_xyyy_xzz_1[j];

                    tg_xxyyy_xyyy_0[j] = pb_x * tg_xyyy_xyyy_0[j] + wp_x[j] * tg_xyyy_xyyy_1[j] + 0.5 * fl1_fx * tg_yyy_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyy_yyy_1[j];

                    tg_xxyyy_xyyz_0[j] = pb_x * tg_xyyy_xyyz_0[j] + wp_x[j] * tg_xyyy_xyyz_1[j] + 0.5 * fl1_fx * tg_yyy_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yyz_1[j];

                    tg_xxyyy_xyzz_0[j] = pb_x * tg_xyyy_xyzz_0[j] + wp_x[j] * tg_xyyy_xyzz_1[j] + 0.5 * fl1_fx * tg_yyy_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yzz_1[j];

                    tg_xxyyy_xzzz_0[j] = pb_x * tg_xyyy_xzzz_0[j] + wp_x[j] * tg_xyyy_xzzz_1[j] + 0.5 * fl1_fx * tg_yyy_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_zzz_1[j];

                    tg_xxyyy_yyyy_0[j] = pb_x * tg_xyyy_yyyy_0[j] + wp_x[j] * tg_xyyy_yyyy_1[j] + 0.5 * fl1_fx * tg_yyy_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyyy_1[j];

                    tg_xxyyy_yyyz_0[j] = pb_x * tg_xyyy_yyyz_0[j] + wp_x[j] * tg_xyyy_yyyz_1[j] + 0.5 * fl1_fx * tg_yyy_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyyz_1[j];

                    tg_xxyyy_yyzz_0[j] = pb_x * tg_xyyy_yyzz_0[j] + wp_x[j] * tg_xyyy_yyzz_1[j] + 0.5 * fl1_fx * tg_yyy_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyzz_1[j];

                    tg_xxyyy_yzzz_0[j] = pb_x * tg_xyyy_yzzz_0[j] + wp_x[j] * tg_xyyy_yzzz_1[j] + 0.5 * fl1_fx * tg_yyy_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yzzz_1[j];

                    tg_xxyyy_zzzz_0[j] = pb_x * tg_xyyy_zzzz_0[j] + wp_x[j] * tg_xyyy_zzzz_1[j] + 0.5 * fl1_fx * tg_yyy_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_zzzz_1[j];

                    tg_xxyyz_xxxx_0[j] = pb_x * tg_xyyz_xxxx_0[j] + wp_x[j] * tg_xyyz_xxxx_1[j] + 0.5 * fl1_fx * tg_yyz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyyz_xxx_1[j];

                    tg_xxyyz_xxxy_0[j] = pb_x * tg_xyyz_xxxy_0[j] + wp_x[j] * tg_xyyz_xxxy_1[j] + 0.5 * fl1_fx * tg_yyz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyyz_xxy_1[j];

                    tg_xxyyz_xxxz_0[j] = pb_x * tg_xyyz_xxxz_0[j] + wp_x[j] * tg_xyyz_xxxz_1[j] + 0.5 * fl1_fx * tg_yyz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyyz_xxz_1[j];

                    tg_xxyyz_xxyy_0[j] = pb_x * tg_xyyz_xxyy_0[j] + wp_x[j] * tg_xyyz_xxyy_1[j] + 0.5 * fl1_fx * tg_yyz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxyy_1[j] + fl1_fxn * tg_xyyz_xyy_1[j];

                    tg_xxyyz_xxyz_0[j] = pb_x * tg_xyyz_xxyz_0[j] + wp_x[j] * tg_xyyz_xxyz_1[j] + 0.5 * fl1_fx * tg_yyz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxyz_1[j] + fl1_fxn * tg_xyyz_xyz_1[j];

                    tg_xxyyz_xxzz_0[j] = pb_x * tg_xyyz_xxzz_0[j] + wp_x[j] * tg_xyyz_xxzz_1[j] + 0.5 * fl1_fx * tg_yyz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxzz_1[j] + fl1_fxn * tg_xyyz_xzz_1[j];

                    tg_xxyyz_xyyy_0[j] = pb_x * tg_xyyz_xyyy_0[j] + wp_x[j] * tg_xyyz_xyyy_1[j] + 0.5 * fl1_fx * tg_yyz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyyz_yyy_1[j];

                    tg_xxyyz_xyyz_0[j] = pb_x * tg_xyyz_xyyz_0[j] + wp_x[j] * tg_xyyz_xyyz_1[j] + 0.5 * fl1_fx * tg_yyz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yyz_1[j];

                    tg_xxyyz_xyzz_0[j] = pb_x * tg_xyyz_xyzz_0[j] + wp_x[j] * tg_xyyz_xyzz_1[j] + 0.5 * fl1_fx * tg_yyz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yzz_1[j];

                    tg_xxyyz_xzzz_0[j] = pb_x * tg_xyyz_xzzz_0[j] + wp_x[j] * tg_xyyz_xzzz_1[j] + 0.5 * fl1_fx * tg_yyz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_zzz_1[j];

                    tg_xxyyz_yyyy_0[j] = pb_x * tg_xyyz_yyyy_0[j] + wp_x[j] * tg_xyyz_yyyy_1[j] + 0.5 * fl1_fx * tg_yyz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyyy_1[j];

                    tg_xxyyz_yyyz_0[j] = pb_x * tg_xyyz_yyyz_0[j] + wp_x[j] * tg_xyyz_yyyz_1[j] + 0.5 * fl1_fx * tg_yyz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyyz_1[j];

                    tg_xxyyz_yyzz_0[j] = pb_x * tg_xyyz_yyzz_0[j] + wp_x[j] * tg_xyyz_yyzz_1[j] + 0.5 * fl1_fx * tg_yyz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyzz_1[j];

                    tg_xxyyz_yzzz_0[j] = pb_x * tg_xyyz_yzzz_0[j] + wp_x[j] * tg_xyyz_yzzz_1[j] + 0.5 * fl1_fx * tg_yyz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yzzz_1[j];

                    tg_xxyyz_zzzz_0[j] = pb_x * tg_xyyz_zzzz_0[j] + wp_x[j] * tg_xyyz_zzzz_1[j] + 0.5 * fl1_fx * tg_yyz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_zzzz_1[j];

                    tg_xxyzz_xxxx_0[j] = pb_x * tg_xyzz_xxxx_0[j] + wp_x[j] * tg_xyzz_xxxx_1[j] + 0.5 * fl1_fx * tg_yzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyzz_xxx_1[j];

                    tg_xxyzz_xxxy_0[j] = pb_x * tg_xyzz_xxxy_0[j] + wp_x[j] * tg_xyzz_xxxy_1[j] + 0.5 * fl1_fx * tg_yzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyzz_xxy_1[j];

                    tg_xxyzz_xxxz_0[j] = pb_x * tg_xyzz_xxxz_0[j] + wp_x[j] * tg_xyzz_xxxz_1[j] + 0.5 * fl1_fx * tg_yzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyzz_xxz_1[j];

                    tg_xxyzz_xxyy_0[j] = pb_x * tg_xyzz_xxyy_0[j] + wp_x[j] * tg_xyzz_xxyy_1[j] + 0.5 * fl1_fx * tg_yzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxyy_1[j] + fl1_fxn * tg_xyzz_xyy_1[j];

                    tg_xxyzz_xxyz_0[j] = pb_x * tg_xyzz_xxyz_0[j] + wp_x[j] * tg_xyzz_xxyz_1[j] + 0.5 * fl1_fx * tg_yzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxyz_1[j] + fl1_fxn * tg_xyzz_xyz_1[j];

                    tg_xxyzz_xxzz_0[j] = pb_x * tg_xyzz_xxzz_0[j] + wp_x[j] * tg_xyzz_xxzz_1[j] + 0.5 * fl1_fx * tg_yzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxzz_1[j] + fl1_fxn * tg_xyzz_xzz_1[j];

                    tg_xxyzz_xyyy_0[j] = pb_x * tg_xyzz_xyyy_0[j] + wp_x[j] * tg_xyzz_xyyy_1[j] + 0.5 * fl1_fx * tg_yzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyzz_yyy_1[j];

                    tg_xxyzz_xyyz_0[j] = pb_x * tg_xyzz_xyyz_0[j] + wp_x[j] * tg_xyzz_xyyz_1[j] + 0.5 * fl1_fx * tg_yzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yyz_1[j];

                    tg_xxyzz_xyzz_0[j] = pb_x * tg_xyzz_xyzz_0[j] + wp_x[j] * tg_xyzz_xyzz_1[j] + 0.5 * fl1_fx * tg_yzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yzz_1[j];

                    tg_xxyzz_xzzz_0[j] = pb_x * tg_xyzz_xzzz_0[j] + wp_x[j] * tg_xyzz_xzzz_1[j] + 0.5 * fl1_fx * tg_yzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_zzz_1[j];

                    tg_xxyzz_yyyy_0[j] = pb_x * tg_xyzz_yyyy_0[j] + wp_x[j] * tg_xyzz_yyyy_1[j] + 0.5 * fl1_fx * tg_yzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyyy_1[j];

                    tg_xxyzz_yyyz_0[j] = pb_x * tg_xyzz_yyyz_0[j] + wp_x[j] * tg_xyzz_yyyz_1[j] + 0.5 * fl1_fx * tg_yzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyyz_1[j];

                    tg_xxyzz_yyzz_0[j] = pb_x * tg_xyzz_yyzz_0[j] + wp_x[j] * tg_xyzz_yyzz_1[j] + 0.5 * fl1_fx * tg_yzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyzz_1[j];

                    tg_xxyzz_yzzz_0[j] = pb_x * tg_xyzz_yzzz_0[j] + wp_x[j] * tg_xyzz_yzzz_1[j] + 0.5 * fl1_fx * tg_yzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yzzz_1[j];

                    tg_xxyzz_zzzz_0[j] = pb_x * tg_xyzz_zzzz_0[j] + wp_x[j] * tg_xyzz_zzzz_1[j] + 0.5 * fl1_fx * tg_yzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_zzzz_1[j];

                    tg_xxzzz_xxxx_0[j] = pb_x * tg_xzzz_xxxx_0[j] + wp_x[j] * tg_xzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xzzz_xxx_1[j];

                    tg_xxzzz_xxxy_0[j] = pb_x * tg_xzzz_xxxy_0[j] + wp_x[j] * tg_xzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xzzz_xxy_1[j];

                    tg_xxzzz_xxxz_0[j] = pb_x * tg_xzzz_xxxz_0[j] + wp_x[j] * tg_xzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xzzz_xxz_1[j];

                    tg_xxzzz_xxyy_0[j] = pb_x * tg_xzzz_xxyy_0[j] + wp_x[j] * tg_xzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyy_1[j] + fl1_fxn * tg_xzzz_xyy_1[j];

                    tg_xxzzz_xxyz_0[j] = pb_x * tg_xzzz_xxyz_0[j] + wp_x[j] * tg_xzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyz_1[j] + fl1_fxn * tg_xzzz_xyz_1[j];

                    tg_xxzzz_xxzz_0[j] = pb_x * tg_xzzz_xxzz_0[j] + wp_x[j] * tg_xzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxzz_1[j] + fl1_fxn * tg_xzzz_xzz_1[j];

                    tg_xxzzz_xyyy_0[j] = pb_x * tg_xzzz_xyyy_0[j] + wp_x[j] * tg_xzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xzzz_yyy_1[j];

                    tg_xxzzz_xyyz_0[j] = pb_x * tg_xzzz_xyyz_0[j] + wp_x[j] * tg_xzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yyz_1[j];

                    tg_xxzzz_xyzz_0[j] = pb_x * tg_xzzz_xyzz_0[j] + wp_x[j] * tg_xzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yzz_1[j];

                    tg_xxzzz_xzzz_0[j] = pb_x * tg_xzzz_xzzz_0[j] + wp_x[j] * tg_xzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_zzz_1[j];

                    tg_xxzzz_yyyy_0[j] = pb_x * tg_xzzz_yyyy_0[j] + wp_x[j] * tg_xzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyy_1[j];

                    tg_xxzzz_yyyz_0[j] = pb_x * tg_xzzz_yyyz_0[j] + wp_x[j] * tg_xzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyz_1[j];

                    tg_xxzzz_yyzz_0[j] = pb_x * tg_xzzz_yyzz_0[j] + wp_x[j] * tg_xzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyzz_1[j];

                    tg_xxzzz_yzzz_0[j] = pb_x * tg_xzzz_yzzz_0[j] + wp_x[j] * tg_xzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yzzz_1[j];

                    tg_xxzzz_zzzz_0[j] = pb_x * tg_xzzz_zzzz_0[j] + wp_x[j] * tg_xzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_zzzz_1[j];

                    tg_xyyyy_xxxx_0[j] = pb_x * tg_yyyy_xxxx_0[j] + wp_x[j] * tg_yyyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyy_xxx_1[j];

                    tg_xyyyy_xxxy_0[j] = pb_x * tg_yyyy_xxxy_0[j] + wp_x[j] * tg_yyyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxy_1[j];

                    tg_xyyyy_xxxz_0[j] = pb_x * tg_yyyy_xxxz_0[j] + wp_x[j] * tg_yyyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyy_xxz_1[j];

                    tg_xyyyy_xxyy_0[j] = pb_x * tg_yyyy_xxyy_0[j] + wp_x[j] * tg_yyyy_xxyy_1[j] + fl1_fxn * tg_yyyy_xyy_1[j];

                    tg_xyyyy_xxyz_0[j] = pb_x * tg_yyyy_xxyz_0[j] + wp_x[j] * tg_yyyy_xxyz_1[j] + fl1_fxn * tg_yyyy_xyz_1[j];

                    tg_xyyyy_xxzz_0[j] = pb_x * tg_yyyy_xxzz_0[j] + wp_x[j] * tg_yyyy_xxzz_1[j] + fl1_fxn * tg_yyyy_xzz_1[j];

                    tg_xyyyy_xyyy_0[j] = pb_x * tg_yyyy_xyyy_0[j] + wp_x[j] * tg_yyyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyy_1[j];

                    tg_xyyyy_xyyz_0[j] = pb_x * tg_yyyy_xyyz_0[j] + wp_x[j] * tg_yyyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSG_158_237(      CMemBlock2D<double>& primBuffer,
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
                                             {5, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
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

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_yyyy_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 157); 

                auto tg_yyyy_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 158); 

                auto tg_yyyy_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 159); 

                auto tg_yyyy_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 160); 

                auto tg_yyyy_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 161); 

                auto tg_yyyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 162); 

                auto tg_yyyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 187); 

                auto tg_yyzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 208); 

                auto tg_yzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 224); 

                auto tg_yyyy_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 150); 

                auto tg_yyyy_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 151); 

                auto tg_yyyy_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 152); 

                auto tg_yyyy_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 153); 

                auto tg_yyyy_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 154); 

                auto tg_yyyy_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 155); 

                auto tg_yyyy_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 156); 

                auto tg_yyyy_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 157); 

                auto tg_yyyy_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 158); 

                auto tg_yyyy_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 159); 

                auto tg_yyyy_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 160); 

                auto tg_yyyy_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 161); 

                auto tg_yyyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 162); 

                auto tg_yyyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 187); 

                auto tg_yyzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 208); 

                auto tg_yzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 224); 

                auto tg_yyy_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 90); 

                auto tg_yyy_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 91); 

                auto tg_yyy_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 92); 

                auto tg_yyy_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 93); 

                auto tg_yyy_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 94); 

                auto tg_yyy_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 95); 

                auto tg_yyy_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 96); 

                auto tg_yyy_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 97); 

                auto tg_yyy_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 98); 

                auto tg_yyy_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 99); 

                auto tg_yyy_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 100); 

                auto tg_yyy_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 101); 

                auto tg_yyy_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 90); 

                auto tg_yyy_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 91); 

                auto tg_yyy_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 92); 

                auto tg_yyy_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 93); 

                auto tg_yyy_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 94); 

                auto tg_yyy_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 95); 

                auto tg_yyy_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 96); 

                auto tg_yyy_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 97); 

                auto tg_yyy_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 98); 

                auto tg_yyy_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 99); 

                auto tg_yyy_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 100); 

                auto tg_yyy_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 101); 

                auto tg_yyyy_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 100); 

                auto tg_yyyy_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 101); 

                auto tg_yyyy_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 102); 

                auto tg_yyyy_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 103); 

                auto tg_yyyy_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 104); 

                auto tg_yyyy_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 105); 

                auto tg_yyyy_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 106); 

                auto tg_yyyy_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 107); 

                auto tg_yyyy_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 108); 

                auto tg_yyyy_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 109); 

                auto tg_yyyz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 110); 

                auto tg_yyyz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 111); 

                auto tg_yyyz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 112); 

                auto tg_yyyz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 113); 

                auto tg_yyyz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 114); 

                auto tg_yyyz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 115); 

                auto tg_yyyz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 116); 

                auto tg_yyyz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 117); 

                auto tg_yyyz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 118); 

                auto tg_yyyz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 119); 

                auto tg_yyzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 120); 

                auto tg_yyzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 121); 

                auto tg_yyzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 122); 

                auto tg_yyzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 123); 

                auto tg_yyzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 124); 

                auto tg_yyzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 125); 

                auto tg_yyzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 126); 

                auto tg_yyzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 127); 

                auto tg_yyzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 128); 

                auto tg_yyzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 129); 

                auto tg_yzzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 130); 

                auto tg_yzzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 131); 

                auto tg_yzzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 132); 

                auto tg_yzzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 133); 

                auto tg_yzzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 134); 

                auto tg_yzzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 135); 

                auto tg_yzzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 136); 

                auto tg_yzzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 137); 

                auto tg_yzzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 138); 

                auto tg_yzzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 139); 

                auto tg_zzzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 140); 

                auto tg_zzzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 141); 

                auto tg_zzzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 142); 

                auto tg_zzzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 143); 

                auto tg_zzzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 144); 

                auto tg_zzzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 145); 

                auto tg_zzzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 146); 

                auto tg_zzzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 147); 

                auto tg_zzzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 148); 

                auto tg_zzzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 149); 

                // set up pointers to integrals

                auto tg_xyyyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 158); 

                auto tg_xyyyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 159); 

                auto tg_xyyyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 160); 

                auto tg_xyyyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 161); 

                auto tg_xyyyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 162); 

                auto tg_xyyyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 163); 

                auto tg_xyyyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 164); 

                auto tg_xyyyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 165); 

                auto tg_xyyyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 166); 

                auto tg_xyyyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 167); 

                auto tg_xyyyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 168); 

                auto tg_xyyyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 169); 

                auto tg_xyyyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 170); 

                auto tg_xyyyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 171); 

                auto tg_xyyyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 172); 

                auto tg_xyyyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 173); 

                auto tg_xyyyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 174); 

                auto tg_xyyyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 175); 

                auto tg_xyyyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 176); 

                auto tg_xyyyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 177); 

                auto tg_xyyyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 178); 

                auto tg_xyyyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 179); 

                auto tg_xyyzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 180); 

                auto tg_xyyzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 181); 

                auto tg_xyyzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 182); 

                auto tg_xyyzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 183); 

                auto tg_xyyzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 184); 

                auto tg_xyyzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 185); 

                auto tg_xyyzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 186); 

                auto tg_xyyzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 187); 

                auto tg_xyyzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 188); 

                auto tg_xyyzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 189); 

                auto tg_xyyzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 190); 

                auto tg_xyyzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 191); 

                auto tg_xyyzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 192); 

                auto tg_xyyzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 193); 

                auto tg_xyyzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 194); 

                auto tg_xyzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 195); 

                auto tg_xyzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 196); 

                auto tg_xyzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 197); 

                auto tg_xyzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 198); 

                auto tg_xyzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 199); 

                auto tg_xyzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 200); 

                auto tg_xyzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 201); 

                auto tg_xyzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 202); 

                auto tg_xyzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 203); 

                auto tg_xyzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 204); 

                auto tg_xyzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 205); 

                auto tg_xyzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 206); 

                auto tg_xyzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 207); 

                auto tg_xyzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 208); 

                auto tg_xyzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 209); 

                auto tg_xzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 210); 

                auto tg_xzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 211); 

                auto tg_xzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 212); 

                auto tg_xzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 213); 

                auto tg_xzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 214); 

                auto tg_xzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 215); 

                auto tg_xzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 216); 

                auto tg_xzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 217); 

                auto tg_xzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 218); 

                auto tg_xzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 219); 

                auto tg_xzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 220); 

                auto tg_xzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 221); 

                auto tg_xzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 222); 

                auto tg_xzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 223); 

                auto tg_xzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 224); 

                auto tg_yyyyy_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 225); 

                auto tg_yyyyy_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 226); 

                auto tg_yyyyy_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 227); 

                auto tg_yyyyy_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 228); 

                auto tg_yyyyy_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 229); 

                auto tg_yyyyy_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 230); 

                auto tg_yyyyy_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 231); 

                auto tg_yyyyy_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 232); 

                auto tg_yyyyy_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 233); 

                auto tg_yyyyy_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 234); 

                auto tg_yyyyy_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 235); 

                auto tg_yyyyy_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 236); 

                // Batch of Integrals (158,237)

                #pragma omp simd aligned(fxn, fza, tg_xyyyy_xyzz_0, tg_xyyyy_xzzz_0, tg_xyyyy_yyyy_0, \
                                         tg_xyyyy_yyyz_0, tg_xyyyy_yyzz_0, tg_xyyyy_yzzz_0, tg_xyyyy_zzzz_0, tg_xyyyz_xxxx_0, \
                                         tg_xyyyz_xxxy_0, tg_xyyyz_xxxz_0, tg_xyyyz_xxyy_0, tg_xyyyz_xxyz_0, tg_xyyyz_xxzz_0, \
                                         tg_xyyyz_xyyy_0, tg_xyyyz_xyyz_0, tg_xyyyz_xyzz_0, tg_xyyyz_xzzz_0, tg_xyyyz_yyyy_0, \
                                         tg_xyyyz_yyyz_0, tg_xyyyz_yyzz_0, tg_xyyyz_yzzz_0, tg_xyyyz_zzzz_0, tg_xyyzz_xxxx_0, \
                                         tg_xyyzz_xxxy_0, tg_xyyzz_xxxz_0, tg_xyyzz_xxyy_0, tg_xyyzz_xxyz_0, tg_xyyzz_xxzz_0, \
                                         tg_xyyzz_xyyy_0, tg_xyyzz_xyyz_0, tg_xyyzz_xyzz_0, tg_xyyzz_xzzz_0, tg_xyyzz_yyyy_0, \
                                         tg_xyyzz_yyyz_0, tg_xyyzz_yyzz_0, tg_xyyzz_yzzz_0, tg_xyyzz_zzzz_0, tg_xyzzz_xxxx_0, \
                                         tg_xyzzz_xxxy_0, tg_xyzzz_xxxz_0, tg_xyzzz_xxyy_0, tg_xyzzz_xxyz_0, tg_xyzzz_xxzz_0, \
                                         tg_xyzzz_xyyy_0, tg_xyzzz_xyyz_0, tg_xyzzz_xyzz_0, tg_xyzzz_xzzz_0, tg_xyzzz_yyyy_0, \
                                         tg_xyzzz_yyyz_0, tg_xyzzz_yyzz_0, tg_xyzzz_yzzz_0, tg_xyzzz_zzzz_0, tg_xzzzz_xxxx_0, \
                                         tg_xzzzz_xxxy_0, tg_xzzzz_xxxz_0, tg_xzzzz_xxyy_0, tg_xzzzz_xxyz_0, tg_xzzzz_xxzz_0, \
                                         tg_xzzzz_xyyy_0, tg_xzzzz_xyyz_0, tg_xzzzz_xyzz_0, tg_xzzzz_xzzz_0, tg_xzzzz_yyyy_0, \
                                         tg_xzzzz_yyyz_0, tg_xzzzz_yyzz_0, tg_xzzzz_yzzz_0, tg_xzzzz_zzzz_0, tg_yyy_xxxx_0, \
                                         tg_yyy_xxxx_1, tg_yyy_xxxy_0, tg_yyy_xxxy_1, tg_yyy_xxxz_0, tg_yyy_xxxz_1, \
                                         tg_yyy_xxyy_0, tg_yyy_xxyy_1, tg_yyy_xxyz_0, tg_yyy_xxyz_1, tg_yyy_xxzz_0, \
                                         tg_yyy_xxzz_1, tg_yyy_xyyy_0, tg_yyy_xyyy_1, tg_yyy_xyyz_0, tg_yyy_xyyz_1, \
                                         tg_yyy_xyzz_0, tg_yyy_xyzz_1, tg_yyy_xzzz_0, tg_yyy_xzzz_1, tg_yyy_yyyy_0, \
                                         tg_yyy_yyyy_1, tg_yyy_yyyz_0, tg_yyy_yyyz_1, tg_yyyy_xxx_1, tg_yyyy_xxxx_0, \
                                         tg_yyyy_xxxx_1, tg_yyyy_xxxy_0, tg_yyyy_xxxy_1, tg_yyyy_xxxz_0, tg_yyyy_xxxz_1, \
                                         tg_yyyy_xxy_1, tg_yyyy_xxyy_0, tg_yyyy_xxyy_1, tg_yyyy_xxyz_0, tg_yyyy_xxyz_1, \
                                         tg_yyyy_xxz_1, tg_yyyy_xxzz_0, tg_yyyy_xxzz_1, tg_yyyy_xyy_1, tg_yyyy_xyyy_0, \
                                         tg_yyyy_xyyy_1, tg_yyyy_xyyz_0, tg_yyyy_xyyz_1, tg_yyyy_xyz_1, tg_yyyy_xyzz_0, \
                                         tg_yyyy_xyzz_1, tg_yyyy_xzz_1, tg_yyyy_xzzz_0, tg_yyyy_xzzz_1, tg_yyyy_yyy_1, \
                                         tg_yyyy_yyyy_0, tg_yyyy_yyyy_1, tg_yyyy_yyyz_0, tg_yyyy_yyyz_1, tg_yyyy_yyz_1, \
                                         tg_yyyy_yyzz_0, tg_yyyy_yyzz_1, tg_yyyy_yzz_1, tg_yyyy_yzzz_0, tg_yyyy_yzzz_1, \
                                         tg_yyyy_zzz_1, tg_yyyy_zzzz_0, tg_yyyy_zzzz_1, tg_yyyyy_xxxx_0, tg_yyyyy_xxxy_0, \
                                         tg_yyyyy_xxxz_0, tg_yyyyy_xxyy_0, tg_yyyyy_xxyz_0, tg_yyyyy_xxzz_0, tg_yyyyy_xyyy_0, \
                                         tg_yyyyy_xyyz_0, tg_yyyyy_xyzz_0, tg_yyyyy_xzzz_0, tg_yyyyy_yyyy_0, tg_yyyyy_yyyz_0, \
                                         tg_yyyz_xxx_1, tg_yyyz_xxxx_0, tg_yyyz_xxxx_1, tg_yyyz_xxxy_0, tg_yyyz_xxxy_1, \
                                         tg_yyyz_xxxz_0, tg_yyyz_xxxz_1, tg_yyyz_xxy_1, tg_yyyz_xxyy_0, tg_yyyz_xxyy_1, \
                                         tg_yyyz_xxyz_0, tg_yyyz_xxyz_1, tg_yyyz_xxz_1, tg_yyyz_xxzz_0, tg_yyyz_xxzz_1, \
                                         tg_yyyz_xyy_1, tg_yyyz_xyyy_0, tg_yyyz_xyyy_1, tg_yyyz_xyyz_0, tg_yyyz_xyyz_1, \
                                         tg_yyyz_xyz_1, tg_yyyz_xyzz_0, tg_yyyz_xyzz_1, tg_yyyz_xzz_1, tg_yyyz_xzzz_0, \
                                         tg_yyyz_xzzz_1, tg_yyyz_yyy_1, tg_yyyz_yyyy_0, tg_yyyz_yyyy_1, tg_yyyz_yyyz_0, \
                                         tg_yyyz_yyyz_1, tg_yyyz_yyz_1, tg_yyyz_yyzz_0, tg_yyyz_yyzz_1, tg_yyyz_yzz_1, \
                                         tg_yyyz_yzzz_0, tg_yyyz_yzzz_1, tg_yyyz_zzz_1, tg_yyyz_zzzz_0, tg_yyyz_zzzz_1, \
                                         tg_yyzz_xxx_1, tg_yyzz_xxxx_0, tg_yyzz_xxxx_1, tg_yyzz_xxxy_0, tg_yyzz_xxxy_1, \
                                         tg_yyzz_xxxz_0, tg_yyzz_xxxz_1, tg_yyzz_xxy_1, tg_yyzz_xxyy_0, tg_yyzz_xxyy_1, \
                                         tg_yyzz_xxyz_0, tg_yyzz_xxyz_1, tg_yyzz_xxz_1, tg_yyzz_xxzz_0, tg_yyzz_xxzz_1, \
                                         tg_yyzz_xyy_1, tg_yyzz_xyyy_0, tg_yyzz_xyyy_1, tg_yyzz_xyyz_0, tg_yyzz_xyyz_1, \
                                         tg_yyzz_xyz_1, tg_yyzz_xyzz_0, tg_yyzz_xyzz_1, tg_yyzz_xzz_1, tg_yyzz_xzzz_0, \
                                         tg_yyzz_xzzz_1, tg_yyzz_yyy_1, tg_yyzz_yyyy_0, tg_yyzz_yyyy_1, tg_yyzz_yyyz_0, \
                                         tg_yyzz_yyyz_1, tg_yyzz_yyz_1, tg_yyzz_yyzz_0, tg_yyzz_yyzz_1, tg_yyzz_yzz_1, \
                                         tg_yyzz_yzzz_0, tg_yyzz_yzzz_1, tg_yyzz_zzz_1, tg_yyzz_zzzz_0, tg_yyzz_zzzz_1, \
                                         tg_yzzz_xxx_1, tg_yzzz_xxxx_0, tg_yzzz_xxxx_1, tg_yzzz_xxxy_0, tg_yzzz_xxxy_1, \
                                         tg_yzzz_xxxz_0, tg_yzzz_xxxz_1, tg_yzzz_xxy_1, tg_yzzz_xxyy_0, tg_yzzz_xxyy_1, \
                                         tg_yzzz_xxyz_0, tg_yzzz_xxyz_1, tg_yzzz_xxz_1, tg_yzzz_xxzz_0, tg_yzzz_xxzz_1, \
                                         tg_yzzz_xyy_1, tg_yzzz_xyyy_0, tg_yzzz_xyyy_1, tg_yzzz_xyyz_0, tg_yzzz_xyyz_1, \
                                         tg_yzzz_xyz_1, tg_yzzz_xyzz_0, tg_yzzz_xyzz_1, tg_yzzz_xzz_1, tg_yzzz_xzzz_0, \
                                         tg_yzzz_xzzz_1, tg_yzzz_yyy_1, tg_yzzz_yyyy_0, tg_yzzz_yyyy_1, tg_yzzz_yyyz_0, \
                                         tg_yzzz_yyyz_1, tg_yzzz_yyz_1, tg_yzzz_yyzz_0, tg_yzzz_yyzz_1, tg_yzzz_yzz_1, \
                                         tg_yzzz_yzzz_0, tg_yzzz_yzzz_1, tg_yzzz_zzz_1, tg_yzzz_zzzz_0, tg_yzzz_zzzz_1, \
                                         tg_zzzz_xxx_1, tg_zzzz_xxxx_0, tg_zzzz_xxxx_1, tg_zzzz_xxxy_0, tg_zzzz_xxxy_1, \
                                         tg_zzzz_xxxz_0, tg_zzzz_xxxz_1, tg_zzzz_xxy_1, tg_zzzz_xxyy_0, tg_zzzz_xxyy_1, \
                                         tg_zzzz_xxyz_0, tg_zzzz_xxyz_1, tg_zzzz_xxz_1, tg_zzzz_xxzz_0, tg_zzzz_xxzz_1, \
                                         tg_zzzz_xyy_1, tg_zzzz_xyyy_0, tg_zzzz_xyyy_1, tg_zzzz_xyyz_0, tg_zzzz_xyyz_1, \
                                         tg_zzzz_xyz_1, tg_zzzz_xyzz_0, tg_zzzz_xyzz_1, tg_zzzz_xzz_1, tg_zzzz_xzzz_0, \
                                         tg_zzzz_xzzz_1, tg_zzzz_yyy_1, tg_zzzz_yyyy_0, tg_zzzz_yyyy_1, tg_zzzz_yyyz_0, \
                                         tg_zzzz_yyyz_1, tg_zzzz_yyz_1, tg_zzzz_yyzz_0, tg_zzzz_yyzz_1, tg_zzzz_yzz_1, \
                                         tg_zzzz_yzzz_0, tg_zzzz_yzzz_1, tg_zzzz_zzz_1, tg_zzzz_zzzz_0, tg_zzzz_zzzz_1, wp_x, \
                                         wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xyyyy_xyzz_0[j] = pb_x * tg_yyyy_xyzz_0[j] + wp_x[j] * tg_yyyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yzz_1[j];

                    tg_xyyyy_xzzz_0[j] = pb_x * tg_yyyy_xzzz_0[j] + wp_x[j] * tg_yyyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zzz_1[j];

                    tg_xyyyy_yyyy_0[j] = pb_x * tg_yyyy_yyyy_0[j] + wp_x[j] * tg_yyyy_yyyy_1[j];

                    tg_xyyyy_yyyz_0[j] = pb_x * tg_yyyy_yyyz_0[j] + wp_x[j] * tg_yyyy_yyyz_1[j];

                    tg_xyyyy_yyzz_0[j] = pb_x * tg_yyyy_yyzz_0[j] + wp_x[j] * tg_yyyy_yyzz_1[j];

                    tg_xyyyy_yzzz_0[j] = pb_x * tg_yyyy_yzzz_0[j] + wp_x[j] * tg_yyyy_yzzz_1[j];

                    tg_xyyyy_zzzz_0[j] = pb_x * tg_yyyy_zzzz_0[j] + wp_x[j] * tg_yyyy_zzzz_1[j];

                    tg_xyyyz_xxxx_0[j] = pb_x * tg_yyyz_xxxx_0[j] + wp_x[j] * tg_yyyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyyz_xxx_1[j];

                    tg_xyyyz_xxxy_0[j] = pb_x * tg_yyyz_xxxy_0[j] + wp_x[j] * tg_yyyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxy_1[j];

                    tg_xyyyz_xxxz_0[j] = pb_x * tg_yyyz_xxxz_0[j] + wp_x[j] * tg_yyyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyyz_xxz_1[j];

                    tg_xyyyz_xxyy_0[j] = pb_x * tg_yyyz_xxyy_0[j] + wp_x[j] * tg_yyyz_xxyy_1[j] + fl1_fxn * tg_yyyz_xyy_1[j];

                    tg_xyyyz_xxyz_0[j] = pb_x * tg_yyyz_xxyz_0[j] + wp_x[j] * tg_yyyz_xxyz_1[j] + fl1_fxn * tg_yyyz_xyz_1[j];

                    tg_xyyyz_xxzz_0[j] = pb_x * tg_yyyz_xxzz_0[j] + wp_x[j] * tg_yyyz_xxzz_1[j] + fl1_fxn * tg_yyyz_xzz_1[j];

                    tg_xyyyz_xyyy_0[j] = pb_x * tg_yyyz_xyyy_0[j] + wp_x[j] * tg_yyyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyy_1[j];

                    tg_xyyyz_xyyz_0[j] = pb_x * tg_yyyz_xyyz_0[j] + wp_x[j] * tg_yyyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yyz_1[j];

                    tg_xyyyz_xyzz_0[j] = pb_x * tg_yyyz_xyzz_0[j] + wp_x[j] * tg_yyyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yzz_1[j];

                    tg_xyyyz_xzzz_0[j] = pb_x * tg_yyyz_xzzz_0[j] + wp_x[j] * tg_yyyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zzz_1[j];

                    tg_xyyyz_yyyy_0[j] = pb_x * tg_yyyz_yyyy_0[j] + wp_x[j] * tg_yyyz_yyyy_1[j];

                    tg_xyyyz_yyyz_0[j] = pb_x * tg_yyyz_yyyz_0[j] + wp_x[j] * tg_yyyz_yyyz_1[j];

                    tg_xyyyz_yyzz_0[j] = pb_x * tg_yyyz_yyzz_0[j] + wp_x[j] * tg_yyyz_yyzz_1[j];

                    tg_xyyyz_yzzz_0[j] = pb_x * tg_yyyz_yzzz_0[j] + wp_x[j] * tg_yyyz_yzzz_1[j];

                    tg_xyyyz_zzzz_0[j] = pb_x * tg_yyyz_zzzz_0[j] + wp_x[j] * tg_yyyz_zzzz_1[j];

                    tg_xyyzz_xxxx_0[j] = pb_x * tg_yyzz_xxxx_0[j] + wp_x[j] * tg_yyzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyzz_xxx_1[j];

                    tg_xyyzz_xxxy_0[j] = pb_x * tg_yyzz_xxxy_0[j] + wp_x[j] * tg_yyzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxy_1[j];

                    tg_xyyzz_xxxz_0[j] = pb_x * tg_yyzz_xxxz_0[j] + wp_x[j] * tg_yyzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyzz_xxz_1[j];

                    tg_xyyzz_xxyy_0[j] = pb_x * tg_yyzz_xxyy_0[j] + wp_x[j] * tg_yyzz_xxyy_1[j] + fl1_fxn * tg_yyzz_xyy_1[j];

                    tg_xyyzz_xxyz_0[j] = pb_x * tg_yyzz_xxyz_0[j] + wp_x[j] * tg_yyzz_xxyz_1[j] + fl1_fxn * tg_yyzz_xyz_1[j];

                    tg_xyyzz_xxzz_0[j] = pb_x * tg_yyzz_xxzz_0[j] + wp_x[j] * tg_yyzz_xxzz_1[j] + fl1_fxn * tg_yyzz_xzz_1[j];

                    tg_xyyzz_xyyy_0[j] = pb_x * tg_yyzz_xyyy_0[j] + wp_x[j] * tg_yyzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyy_1[j];

                    tg_xyyzz_xyyz_0[j] = pb_x * tg_yyzz_xyyz_0[j] + wp_x[j] * tg_yyzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yyz_1[j];

                    tg_xyyzz_xyzz_0[j] = pb_x * tg_yyzz_xyzz_0[j] + wp_x[j] * tg_yyzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yzz_1[j];

                    tg_xyyzz_xzzz_0[j] = pb_x * tg_yyzz_xzzz_0[j] + wp_x[j] * tg_yyzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zzz_1[j];

                    tg_xyyzz_yyyy_0[j] = pb_x * tg_yyzz_yyyy_0[j] + wp_x[j] * tg_yyzz_yyyy_1[j];

                    tg_xyyzz_yyyz_0[j] = pb_x * tg_yyzz_yyyz_0[j] + wp_x[j] * tg_yyzz_yyyz_1[j];

                    tg_xyyzz_yyzz_0[j] = pb_x * tg_yyzz_yyzz_0[j] + wp_x[j] * tg_yyzz_yyzz_1[j];

                    tg_xyyzz_yzzz_0[j] = pb_x * tg_yyzz_yzzz_0[j] + wp_x[j] * tg_yyzz_yzzz_1[j];

                    tg_xyyzz_zzzz_0[j] = pb_x * tg_yyzz_zzzz_0[j] + wp_x[j] * tg_yyzz_zzzz_1[j];

                    tg_xyzzz_xxxx_0[j] = pb_x * tg_yzzz_xxxx_0[j] + wp_x[j] * tg_yzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yzzz_xxx_1[j];

                    tg_xyzzz_xxxy_0[j] = pb_x * tg_yzzz_xxxy_0[j] + wp_x[j] * tg_yzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxy_1[j];

                    tg_xyzzz_xxxz_0[j] = pb_x * tg_yzzz_xxxz_0[j] + wp_x[j] * tg_yzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yzzz_xxz_1[j];

                    tg_xyzzz_xxyy_0[j] = pb_x * tg_yzzz_xxyy_0[j] + wp_x[j] * tg_yzzz_xxyy_1[j] + fl1_fxn * tg_yzzz_xyy_1[j];

                    tg_xyzzz_xxyz_0[j] = pb_x * tg_yzzz_xxyz_0[j] + wp_x[j] * tg_yzzz_xxyz_1[j] + fl1_fxn * tg_yzzz_xyz_1[j];

                    tg_xyzzz_xxzz_0[j] = pb_x * tg_yzzz_xxzz_0[j] + wp_x[j] * tg_yzzz_xxzz_1[j] + fl1_fxn * tg_yzzz_xzz_1[j];

                    tg_xyzzz_xyyy_0[j] = pb_x * tg_yzzz_xyyy_0[j] + wp_x[j] * tg_yzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyy_1[j];

                    tg_xyzzz_xyyz_0[j] = pb_x * tg_yzzz_xyyz_0[j] + wp_x[j] * tg_yzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yyz_1[j];

                    tg_xyzzz_xyzz_0[j] = pb_x * tg_yzzz_xyzz_0[j] + wp_x[j] * tg_yzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yzz_1[j];

                    tg_xyzzz_xzzz_0[j] = pb_x * tg_yzzz_xzzz_0[j] + wp_x[j] * tg_yzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zzz_1[j];

                    tg_xyzzz_yyyy_0[j] = pb_x * tg_yzzz_yyyy_0[j] + wp_x[j] * tg_yzzz_yyyy_1[j];

                    tg_xyzzz_yyyz_0[j] = pb_x * tg_yzzz_yyyz_0[j] + wp_x[j] * tg_yzzz_yyyz_1[j];

                    tg_xyzzz_yyzz_0[j] = pb_x * tg_yzzz_yyzz_0[j] + wp_x[j] * tg_yzzz_yyzz_1[j];

                    tg_xyzzz_yzzz_0[j] = pb_x * tg_yzzz_yzzz_0[j] + wp_x[j] * tg_yzzz_yzzz_1[j];

                    tg_xyzzz_zzzz_0[j] = pb_x * tg_yzzz_zzzz_0[j] + wp_x[j] * tg_yzzz_zzzz_1[j];

                    tg_xzzzz_xxxx_0[j] = pb_x * tg_zzzz_xxxx_0[j] + wp_x[j] * tg_zzzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_zzzz_xxx_1[j];

                    tg_xzzzz_xxxy_0[j] = pb_x * tg_zzzz_xxxy_0[j] + wp_x[j] * tg_zzzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxy_1[j];

                    tg_xzzzz_xxxz_0[j] = pb_x * tg_zzzz_xxxz_0[j] + wp_x[j] * tg_zzzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xxz_1[j];

                    tg_xzzzz_xxyy_0[j] = pb_x * tg_zzzz_xxyy_0[j] + wp_x[j] * tg_zzzz_xxyy_1[j] + fl1_fxn * tg_zzzz_xyy_1[j];

                    tg_xzzzz_xxyz_0[j] = pb_x * tg_zzzz_xxyz_0[j] + wp_x[j] * tg_zzzz_xxyz_1[j] + fl1_fxn * tg_zzzz_xyz_1[j];

                    tg_xzzzz_xxzz_0[j] = pb_x * tg_zzzz_xxzz_0[j] + wp_x[j] * tg_zzzz_xxzz_1[j] + fl1_fxn * tg_zzzz_xzz_1[j];

                    tg_xzzzz_xyyy_0[j] = pb_x * tg_zzzz_xyyy_0[j] + wp_x[j] * tg_zzzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyy_1[j];

                    tg_xzzzz_xyyz_0[j] = pb_x * tg_zzzz_xyyz_0[j] + wp_x[j] * tg_zzzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyz_1[j];

                    tg_xzzzz_xyzz_0[j] = pb_x * tg_zzzz_xyzz_0[j] + wp_x[j] * tg_zzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yzz_1[j];

                    tg_xzzzz_xzzz_0[j] = pb_x * tg_zzzz_xzzz_0[j] + wp_x[j] * tg_zzzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzz_1[j];

                    tg_xzzzz_yyyy_0[j] = pb_x * tg_zzzz_yyyy_0[j] + wp_x[j] * tg_zzzz_yyyy_1[j];

                    tg_xzzzz_yyyz_0[j] = pb_x * tg_zzzz_yyyz_0[j] + wp_x[j] * tg_zzzz_yyyz_1[j];

                    tg_xzzzz_yyzz_0[j] = pb_x * tg_zzzz_yyzz_0[j] + wp_x[j] * tg_zzzz_yyzz_1[j];

                    tg_xzzzz_yzzz_0[j] = pb_x * tg_zzzz_yzzz_0[j] + wp_x[j] * tg_zzzz_yzzz_1[j];

                    tg_xzzzz_zzzz_0[j] = pb_x * tg_zzzz_zzzz_0[j] + wp_x[j] * tg_zzzz_zzzz_1[j];

                    tg_yyyyy_xxxx_0[j] = pb_y * tg_yyyy_xxxx_0[j] + wp_y[j] * tg_yyyy_xxxx_1[j] + 2.0 * fl1_fx * tg_yyy_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxx_1[j];

                    tg_yyyyy_xxxy_0[j] = pb_y * tg_yyyy_xxxy_0[j] + wp_y[j] * tg_yyyy_xxxy_1[j] + 2.0 * fl1_fx * tg_yyy_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyy_xxx_1[j];

                    tg_yyyyy_xxxz_0[j] = pb_y * tg_yyyy_xxxz_0[j] + wp_y[j] * tg_yyyy_xxxz_1[j] + 2.0 * fl1_fx * tg_yyy_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxxz_1[j];

                    tg_yyyyy_xxyy_0[j] = pb_y * tg_yyyy_xxyy_0[j] + wp_y[j] * tg_yyyy_xxyy_1[j] + 2.0 * fl1_fx * tg_yyy_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxyy_1[j] + fl1_fxn * tg_yyyy_xxy_1[j];

                    tg_yyyyy_xxyz_0[j] = pb_y * tg_yyyy_xxyz_0[j] + wp_y[j] * tg_yyyy_xxyz_1[j] + 2.0 * fl1_fx * tg_yyy_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xxz_1[j];

                    tg_yyyyy_xxzz_0[j] = pb_y * tg_yyyy_xxzz_0[j] + wp_y[j] * tg_yyyy_xxzz_1[j] + 2.0 * fl1_fx * tg_yyy_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxzz_1[j];

                    tg_yyyyy_xyyy_0[j] = pb_y * tg_yyyy_xyyy_0[j] + wp_y[j] * tg_yyyy_xyyy_1[j] + 2.0 * fl1_fx * tg_yyy_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyy_xyy_1[j];

                    tg_yyyyy_xyyz_0[j] = pb_y * tg_yyyy_xyyz_0[j] + wp_y[j] * tg_yyyy_xyyz_1[j] + 2.0 * fl1_fx * tg_yyy_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyyz_1[j] + fl1_fxn * tg_yyyy_xyz_1[j];

                    tg_yyyyy_xyzz_0[j] = pb_y * tg_yyyy_xyzz_0[j] + wp_y[j] * tg_yyyy_xyzz_1[j] + 2.0 * fl1_fx * tg_yyy_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xzz_1[j];

                    tg_yyyyy_xzzz_0[j] = pb_y * tg_yyyy_xzzz_0[j] + wp_y[j] * tg_yyyy_xzzz_1[j] + 2.0 * fl1_fx * tg_yyy_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xzzz_1[j];

                    tg_yyyyy_yyyy_0[j] = pb_y * tg_yyyy_yyyy_0[j] + wp_y[j] * tg_yyyy_yyyy_1[j] + 2.0 * fl1_fx * tg_yyy_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyy_yyy_1[j];

                    tg_yyyyy_yyyz_0[j] = pb_y * tg_yyyy_yyyz_0[j] + wp_y[j] * tg_yyyy_yyyz_1[j] + 2.0 * fl1_fx * tg_yyy_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyy_yyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSG_237_315(      CMemBlock2D<double>& primBuffer,
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
                                             {5, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yyyy_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 162); 

                auto tg_yyyy_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 187); 

                auto tg_yyzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 208); 

                auto tg_yzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 224); 

                auto tg_yyyy_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 162); 

                auto tg_yyyy_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 163); 

                auto tg_yyyy_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 164); 

                auto tg_yyyz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 165); 

                auto tg_yyyz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 166); 

                auto tg_yyyz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 167); 

                auto tg_yyyz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 168); 

                auto tg_yyyz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 169); 

                auto tg_yyyz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 170); 

                auto tg_yyyz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 171); 

                auto tg_yyyz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 172); 

                auto tg_yyyz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 173); 

                auto tg_yyyz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 174); 

                auto tg_yyyz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 175); 

                auto tg_yyyz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 176); 

                auto tg_yyyz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 177); 

                auto tg_yyyz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 178); 

                auto tg_yyyz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 179); 

                auto tg_yyzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 180); 

                auto tg_yyzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 181); 

                auto tg_yyzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 182); 

                auto tg_yyzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 183); 

                auto tg_yyzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 184); 

                auto tg_yyzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 185); 

                auto tg_yyzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 186); 

                auto tg_yyzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 187); 

                auto tg_yyzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 188); 

                auto tg_yyzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 189); 

                auto tg_yyzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 190); 

                auto tg_yyzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 191); 

                auto tg_yyzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 192); 

                auto tg_yyzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 193); 

                auto tg_yyzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 194); 

                auto tg_yzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 195); 

                auto tg_yzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 196); 

                auto tg_yzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 197); 

                auto tg_yzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 198); 

                auto tg_yzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 199); 

                auto tg_yzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 200); 

                auto tg_yzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 201); 

                auto tg_yzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 202); 

                auto tg_yzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 203); 

                auto tg_yzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 204); 

                auto tg_yzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 205); 

                auto tg_yzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 206); 

                auto tg_yzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 207); 

                auto tg_yzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 208); 

                auto tg_yzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 209); 

                auto tg_zzzz_xxxx_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 210); 

                auto tg_zzzz_xxxy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 211); 

                auto tg_zzzz_xxxz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 212); 

                auto tg_zzzz_xxyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 213); 

                auto tg_zzzz_xxyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 214); 

                auto tg_zzzz_xxzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 215); 

                auto tg_zzzz_xyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 216); 

                auto tg_zzzz_xyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 217); 

                auto tg_zzzz_xyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 218); 

                auto tg_zzzz_xzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 219); 

                auto tg_zzzz_yyyy_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 220); 

                auto tg_zzzz_yyyz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 221); 

                auto tg_zzzz_yyzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 222); 

                auto tg_zzzz_yzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 223); 

                auto tg_zzzz_zzzz_1 = primBuffer.data(pidx_g_4_4_m1 + 225 * idx + 224); 

                auto tg_yyy_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 102); 

                auto tg_yyy_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 103); 

                auto tg_yyy_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 104); 

                auto tg_yyz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 105); 

                auto tg_yyz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 106); 

                auto tg_yyz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 107); 

                auto tg_yyz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 108); 

                auto tg_yyz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 109); 

                auto tg_yyz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 110); 

                auto tg_yyz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 111); 

                auto tg_yyz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 112); 

                auto tg_yyz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 113); 

                auto tg_yyz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 114); 

                auto tg_yyz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 115); 

                auto tg_yyz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 116); 

                auto tg_yyz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 117); 

                auto tg_yyz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 118); 

                auto tg_yyz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 119); 

                auto tg_yzz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 120); 

                auto tg_yzz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 121); 

                auto tg_yzz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 122); 

                auto tg_yzz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 123); 

                auto tg_yzz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 124); 

                auto tg_yzz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 125); 

                auto tg_yzz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 126); 

                auto tg_yzz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 127); 

                auto tg_yzz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 128); 

                auto tg_yzz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 129); 

                auto tg_yzz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 130); 

                auto tg_yzz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 131); 

                auto tg_yzz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 132); 

                auto tg_yzz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 133); 

                auto tg_yzz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 134); 

                auto tg_zzz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 135); 

                auto tg_zzz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 136); 

                auto tg_zzz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 137); 

                auto tg_zzz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 138); 

                auto tg_zzz_xxyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 139); 

                auto tg_zzz_xxzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 140); 

                auto tg_zzz_xyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 141); 

                auto tg_zzz_xyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 142); 

                auto tg_zzz_xyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 143); 

                auto tg_zzz_xzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 144); 

                auto tg_zzz_yyyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 145); 

                auto tg_zzz_yyyz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 146); 

                auto tg_zzz_yyzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 147); 

                auto tg_zzz_yzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 148); 

                auto tg_zzz_zzzz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 149); 

                auto tg_yyy_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 102); 

                auto tg_yyy_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 103); 

                auto tg_yyy_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 104); 

                auto tg_yyz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 105); 

                auto tg_yyz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 106); 

                auto tg_yyz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 107); 

                auto tg_yyz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 108); 

                auto tg_yyz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 109); 

                auto tg_yyz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 110); 

                auto tg_yyz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 111); 

                auto tg_yyz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 112); 

                auto tg_yyz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 113); 

                auto tg_yyz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 114); 

                auto tg_yyz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 115); 

                auto tg_yyz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 116); 

                auto tg_yyz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 117); 

                auto tg_yyz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 118); 

                auto tg_yyz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 119); 

                auto tg_yzz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 120); 

                auto tg_yzz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 121); 

                auto tg_yzz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 122); 

                auto tg_yzz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 123); 

                auto tg_yzz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 124); 

                auto tg_yzz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 125); 

                auto tg_yzz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 126); 

                auto tg_yzz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 127); 

                auto tg_yzz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 128); 

                auto tg_yzz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 129); 

                auto tg_yzz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 130); 

                auto tg_yzz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 131); 

                auto tg_yzz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 132); 

                auto tg_yzz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 133); 

                auto tg_yzz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 134); 

                auto tg_zzz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 135); 

                auto tg_zzz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 136); 

                auto tg_zzz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 137); 

                auto tg_zzz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 138); 

                auto tg_zzz_xxyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 139); 

                auto tg_zzz_xxzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 140); 

                auto tg_zzz_xyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 141); 

                auto tg_zzz_xyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 142); 

                auto tg_zzz_xyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 143); 

                auto tg_zzz_xzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 144); 

                auto tg_zzz_yyyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 145); 

                auto tg_zzz_yyyz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 146); 

                auto tg_zzz_yyzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 147); 

                auto tg_zzz_yzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 148); 

                auto tg_zzz_zzzz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 149); 

                auto tg_yyyy_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 108); 

                auto tg_yyyy_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 109); 

                auto tg_yyyz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 110); 

                auto tg_yyyz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 111); 

                auto tg_yyyz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 112); 

                auto tg_yyyz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 113); 

                auto tg_yyyz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 114); 

                auto tg_yyyz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 115); 

                auto tg_yyyz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 116); 

                auto tg_yyyz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 117); 

                auto tg_yyyz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 118); 

                auto tg_yyyz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 119); 

                auto tg_yyzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 120); 

                auto tg_yyzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 121); 

                auto tg_yyzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 122); 

                auto tg_yyzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 123); 

                auto tg_yyzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 124); 

                auto tg_yyzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 125); 

                auto tg_yyzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 126); 

                auto tg_yyzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 127); 

                auto tg_yyzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 128); 

                auto tg_yyzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 129); 

                auto tg_yzzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 130); 

                auto tg_yzzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 131); 

                auto tg_yzzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 132); 

                auto tg_yzzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 133); 

                auto tg_yzzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 134); 

                auto tg_yzzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 135); 

                auto tg_yzzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 136); 

                auto tg_yzzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 137); 

                auto tg_yzzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 138); 

                auto tg_yzzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 139); 

                auto tg_zzzz_xxx_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 140); 

                auto tg_zzzz_xxy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 141); 

                auto tg_zzzz_xxz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 142); 

                auto tg_zzzz_xyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 143); 

                auto tg_zzzz_xyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 144); 

                auto tg_zzzz_xzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 145); 

                auto tg_zzzz_yyy_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 146); 

                auto tg_zzzz_yyz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 147); 

                auto tg_zzzz_yzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 148); 

                auto tg_zzzz_zzz_1 = primBuffer.data(pidx_g_4_3_m1 + 150 * idx + 149); 

                // set up pointers to integrals

                auto tg_yyyyy_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 237); 

                auto tg_yyyyy_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 238); 

                auto tg_yyyyy_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 239); 

                auto tg_yyyyz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 240); 

                auto tg_yyyyz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 241); 

                auto tg_yyyyz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 242); 

                auto tg_yyyyz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 243); 

                auto tg_yyyyz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 244); 

                auto tg_yyyyz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 245); 

                auto tg_yyyyz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 246); 

                auto tg_yyyyz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 247); 

                auto tg_yyyyz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 248); 

                auto tg_yyyyz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 249); 

                auto tg_yyyyz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 250); 

                auto tg_yyyyz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 251); 

                auto tg_yyyyz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 252); 

                auto tg_yyyyz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 253); 

                auto tg_yyyyz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 254); 

                auto tg_yyyzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 255); 

                auto tg_yyyzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 256); 

                auto tg_yyyzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 257); 

                auto tg_yyyzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 258); 

                auto tg_yyyzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 259); 

                auto tg_yyyzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 260); 

                auto tg_yyyzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 261); 

                auto tg_yyyzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 262); 

                auto tg_yyyzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 263); 

                auto tg_yyyzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 264); 

                auto tg_yyyzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 265); 

                auto tg_yyyzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 266); 

                auto tg_yyyzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 267); 

                auto tg_yyyzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 268); 

                auto tg_yyyzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 269); 

                auto tg_yyzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 270); 

                auto tg_yyzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 271); 

                auto tg_yyzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 272); 

                auto tg_yyzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 273); 

                auto tg_yyzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 274); 

                auto tg_yyzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 275); 

                auto tg_yyzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 276); 

                auto tg_yyzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 277); 

                auto tg_yyzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 278); 

                auto tg_yyzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 279); 

                auto tg_yyzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 280); 

                auto tg_yyzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 281); 

                auto tg_yyzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 282); 

                auto tg_yyzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 283); 

                auto tg_yyzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 284); 

                auto tg_yzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 285); 

                auto tg_yzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 286); 

                auto tg_yzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 287); 

                auto tg_yzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 288); 

                auto tg_yzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 289); 

                auto tg_yzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 290); 

                auto tg_yzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 291); 

                auto tg_yzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 292); 

                auto tg_yzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 293); 

                auto tg_yzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 294); 

                auto tg_yzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 295); 

                auto tg_yzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 296); 

                auto tg_yzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 297); 

                auto tg_yzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 298); 

                auto tg_yzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 299); 

                auto tg_zzzzz_xxxx_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 300); 

                auto tg_zzzzz_xxxy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 301); 

                auto tg_zzzzz_xxxz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 302); 

                auto tg_zzzzz_xxyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 303); 

                auto tg_zzzzz_xxyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 304); 

                auto tg_zzzzz_xxzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 305); 

                auto tg_zzzzz_xyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 306); 

                auto tg_zzzzz_xyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 307); 

                auto tg_zzzzz_xyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 308); 

                auto tg_zzzzz_xzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 309); 

                auto tg_zzzzz_yyyy_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 310); 

                auto tg_zzzzz_yyyz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 311); 

                auto tg_zzzzz_yyzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 312); 

                auto tg_zzzzz_yzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 313); 

                auto tg_zzzzz_zzzz_0 = primBuffer.data(pidx_g_5_4_m0 + 315 * idx + 314); 

                // Batch of Integrals (237,315)

                #pragma omp simd aligned(fxn, fza, tg_yyy_yyzz_0, tg_yyy_yyzz_1, tg_yyy_yzzz_0, tg_yyy_yzzz_1, \
                                         tg_yyy_zzzz_0, tg_yyy_zzzz_1, tg_yyyy_yyzz_0, tg_yyyy_yyzz_1, tg_yyyy_yzz_1, \
                                         tg_yyyy_yzzz_0, tg_yyyy_yzzz_1, tg_yyyy_zzz_1, tg_yyyy_zzzz_0, tg_yyyy_zzzz_1, \
                                         tg_yyyyy_yyzz_0, tg_yyyyy_yzzz_0, tg_yyyyy_zzzz_0, tg_yyyyz_xxxx_0, tg_yyyyz_xxxy_0, \
                                         tg_yyyyz_xxxz_0, tg_yyyyz_xxyy_0, tg_yyyyz_xxyz_0, tg_yyyyz_xxzz_0, tg_yyyyz_xyyy_0, \
                                         tg_yyyyz_xyyz_0, tg_yyyyz_xyzz_0, tg_yyyyz_xzzz_0, tg_yyyyz_yyyy_0, tg_yyyyz_yyyz_0, \
                                         tg_yyyyz_yyzz_0, tg_yyyyz_yzzz_0, tg_yyyyz_zzzz_0, tg_yyyz_xxx_1, tg_yyyz_xxxx_0, \
                                         tg_yyyz_xxxx_1, tg_yyyz_xxxy_0, tg_yyyz_xxxy_1, tg_yyyz_xxxz_0, tg_yyyz_xxxz_1, \
                                         tg_yyyz_xxy_1, tg_yyyz_xxyy_0, tg_yyyz_xxyy_1, tg_yyyz_xxyz_0, tg_yyyz_xxyz_1, \
                                         tg_yyyz_xxz_1, tg_yyyz_xxzz_0, tg_yyyz_xxzz_1, tg_yyyz_xyy_1, tg_yyyz_xyyy_0, \
                                         tg_yyyz_xyyy_1, tg_yyyz_xyyz_0, tg_yyyz_xyyz_1, tg_yyyz_xyz_1, tg_yyyz_xyzz_0, \
                                         tg_yyyz_xyzz_1, tg_yyyz_xzz_1, tg_yyyz_xzzz_0, tg_yyyz_xzzz_1, tg_yyyz_yyy_1, \
                                         tg_yyyz_yyyy_0, tg_yyyz_yyyy_1, tg_yyyz_yyyz_0, tg_yyyz_yyyz_1, tg_yyyz_yyz_1, \
                                         tg_yyyz_yyzz_0, tg_yyyz_yyzz_1, tg_yyyz_yzz_1, tg_yyyz_yzzz_0, tg_yyyz_yzzz_1, \
                                         tg_yyyz_zzz_1, tg_yyyz_zzzz_0, tg_yyyz_zzzz_1, tg_yyyzz_xxxx_0, tg_yyyzz_xxxy_0, \
                                         tg_yyyzz_xxxz_0, tg_yyyzz_xxyy_0, tg_yyyzz_xxyz_0, tg_yyyzz_xxzz_0, tg_yyyzz_xyyy_0, \
                                         tg_yyyzz_xyyz_0, tg_yyyzz_xyzz_0, tg_yyyzz_xzzz_0, tg_yyyzz_yyyy_0, tg_yyyzz_yyyz_0, \
                                         tg_yyyzz_yyzz_0, tg_yyyzz_yzzz_0, tg_yyyzz_zzzz_0, tg_yyz_xxxx_0, tg_yyz_xxxx_1, \
                                         tg_yyz_xxxy_0, tg_yyz_xxxy_1, tg_yyz_xxxz_0, tg_yyz_xxxz_1, tg_yyz_xxyy_0, \
                                         tg_yyz_xxyy_1, tg_yyz_xxyz_0, tg_yyz_xxyz_1, tg_yyz_xxzz_0, tg_yyz_xxzz_1, \
                                         tg_yyz_xyyy_0, tg_yyz_xyyy_1, tg_yyz_xyyz_0, tg_yyz_xyyz_1, tg_yyz_xyzz_0, \
                                         tg_yyz_xyzz_1, tg_yyz_xzzz_0, tg_yyz_xzzz_1, tg_yyz_yyyy_0, tg_yyz_yyyy_1, \
                                         tg_yyz_yyyz_0, tg_yyz_yyyz_1, tg_yyz_yyzz_0, tg_yyz_yyzz_1, tg_yyz_yzzz_0, \
                                         tg_yyz_yzzz_1, tg_yyz_zzzz_0, tg_yyz_zzzz_1, tg_yyzz_xxx_1, tg_yyzz_xxxx_0, \
                                         tg_yyzz_xxxx_1, tg_yyzz_xxxy_0, tg_yyzz_xxxy_1, tg_yyzz_xxxz_0, tg_yyzz_xxxz_1, \
                                         tg_yyzz_xxy_1, tg_yyzz_xxyy_0, tg_yyzz_xxyy_1, tg_yyzz_xxyz_0, tg_yyzz_xxyz_1, \
                                         tg_yyzz_xxz_1, tg_yyzz_xxzz_0, tg_yyzz_xxzz_1, tg_yyzz_xyy_1, tg_yyzz_xyyy_0, \
                                         tg_yyzz_xyyy_1, tg_yyzz_xyyz_0, tg_yyzz_xyyz_1, tg_yyzz_xyz_1, tg_yyzz_xyzz_0, \
                                         tg_yyzz_xyzz_1, tg_yyzz_xzz_1, tg_yyzz_xzzz_0, tg_yyzz_xzzz_1, tg_yyzz_yyy_1, \
                                         tg_yyzz_yyyy_0, tg_yyzz_yyyy_1, tg_yyzz_yyyz_0, tg_yyzz_yyyz_1, tg_yyzz_yyz_1, \
                                         tg_yyzz_yyzz_0, tg_yyzz_yyzz_1, tg_yyzz_yzz_1, tg_yyzz_yzzz_0, tg_yyzz_yzzz_1, \
                                         tg_yyzz_zzz_1, tg_yyzz_zzzz_0, tg_yyzz_zzzz_1, tg_yyzzz_xxxx_0, tg_yyzzz_xxxy_0, \
                                         tg_yyzzz_xxxz_0, tg_yyzzz_xxyy_0, tg_yyzzz_xxyz_0, tg_yyzzz_xxzz_0, tg_yyzzz_xyyy_0, \
                                         tg_yyzzz_xyyz_0, tg_yyzzz_xyzz_0, tg_yyzzz_xzzz_0, tg_yyzzz_yyyy_0, tg_yyzzz_yyyz_0, \
                                         tg_yyzzz_yyzz_0, tg_yyzzz_yzzz_0, tg_yyzzz_zzzz_0, tg_yzz_xxxx_0, tg_yzz_xxxx_1, \
                                         tg_yzz_xxxy_0, tg_yzz_xxxy_1, tg_yzz_xxxz_0, tg_yzz_xxxz_1, tg_yzz_xxyy_0, \
                                         tg_yzz_xxyy_1, tg_yzz_xxyz_0, tg_yzz_xxyz_1, tg_yzz_xxzz_0, tg_yzz_xxzz_1, \
                                         tg_yzz_xyyy_0, tg_yzz_xyyy_1, tg_yzz_xyyz_0, tg_yzz_xyyz_1, tg_yzz_xyzz_0, \
                                         tg_yzz_xyzz_1, tg_yzz_xzzz_0, tg_yzz_xzzz_1, tg_yzz_yyyy_0, tg_yzz_yyyy_1, \
                                         tg_yzz_yyyz_0, tg_yzz_yyyz_1, tg_yzz_yyzz_0, tg_yzz_yyzz_1, tg_yzz_yzzz_0, \
                                         tg_yzz_yzzz_1, tg_yzz_zzzz_0, tg_yzz_zzzz_1, tg_yzzz_xxx_1, tg_yzzz_xxxx_0, \
                                         tg_yzzz_xxxx_1, tg_yzzz_xxxy_0, tg_yzzz_xxxy_1, tg_yzzz_xxxz_0, tg_yzzz_xxxz_1, \
                                         tg_yzzz_xxy_1, tg_yzzz_xxyy_0, tg_yzzz_xxyy_1, tg_yzzz_xxyz_0, tg_yzzz_xxyz_1, \
                                         tg_yzzz_xxz_1, tg_yzzz_xxzz_0, tg_yzzz_xxzz_1, tg_yzzz_xyy_1, tg_yzzz_xyyy_0, \
                                         tg_yzzz_xyyy_1, tg_yzzz_xyyz_0, tg_yzzz_xyyz_1, tg_yzzz_xyz_1, tg_yzzz_xyzz_0, \
                                         tg_yzzz_xyzz_1, tg_yzzz_xzz_1, tg_yzzz_xzzz_0, tg_yzzz_xzzz_1, tg_yzzz_yyy_1, \
                                         tg_yzzz_yyyy_0, tg_yzzz_yyyy_1, tg_yzzz_yyyz_0, tg_yzzz_yyyz_1, tg_yzzz_yyz_1, \
                                         tg_yzzz_yyzz_0, tg_yzzz_yyzz_1, tg_yzzz_yzz_1, tg_yzzz_yzzz_0, tg_yzzz_yzzz_1, \
                                         tg_yzzz_zzz_1, tg_yzzz_zzzz_0, tg_yzzz_zzzz_1, tg_yzzzz_xxxx_0, tg_yzzzz_xxxy_0, \
                                         tg_yzzzz_xxxz_0, tg_yzzzz_xxyy_0, tg_yzzzz_xxyz_0, tg_yzzzz_xxzz_0, tg_yzzzz_xyyy_0, \
                                         tg_yzzzz_xyyz_0, tg_yzzzz_xyzz_0, tg_yzzzz_xzzz_0, tg_yzzzz_yyyy_0, tg_yzzzz_yyyz_0, \
                                         tg_yzzzz_yyzz_0, tg_yzzzz_yzzz_0, tg_yzzzz_zzzz_0, tg_zzz_xxxx_0, tg_zzz_xxxx_1, \
                                         tg_zzz_xxxy_0, tg_zzz_xxxy_1, tg_zzz_xxxz_0, tg_zzz_xxxz_1, tg_zzz_xxyy_0, \
                                         tg_zzz_xxyy_1, tg_zzz_xxyz_0, tg_zzz_xxyz_1, tg_zzz_xxzz_0, tg_zzz_xxzz_1, \
                                         tg_zzz_xyyy_0, tg_zzz_xyyy_1, tg_zzz_xyyz_0, tg_zzz_xyyz_1, tg_zzz_xyzz_0, \
                                         tg_zzz_xyzz_1, tg_zzz_xzzz_0, tg_zzz_xzzz_1, tg_zzz_yyyy_0, tg_zzz_yyyy_1, \
                                         tg_zzz_yyyz_0, tg_zzz_yyyz_1, tg_zzz_yyzz_0, tg_zzz_yyzz_1, tg_zzz_yzzz_0, \
                                         tg_zzz_yzzz_1, tg_zzz_zzzz_0, tg_zzz_zzzz_1, tg_zzzz_xxx_1, tg_zzzz_xxxx_0, \
                                         tg_zzzz_xxxx_1, tg_zzzz_xxxy_0, tg_zzzz_xxxy_1, tg_zzzz_xxxz_0, tg_zzzz_xxxz_1, \
                                         tg_zzzz_xxy_1, tg_zzzz_xxyy_0, tg_zzzz_xxyy_1, tg_zzzz_xxyz_0, tg_zzzz_xxyz_1, \
                                         tg_zzzz_xxz_1, tg_zzzz_xxzz_0, tg_zzzz_xxzz_1, tg_zzzz_xyy_1, tg_zzzz_xyyy_0, \
                                         tg_zzzz_xyyy_1, tg_zzzz_xyyz_0, tg_zzzz_xyyz_1, tg_zzzz_xyz_1, tg_zzzz_xyzz_0, \
                                         tg_zzzz_xyzz_1, tg_zzzz_xzz_1, tg_zzzz_xzzz_0, tg_zzzz_xzzz_1, tg_zzzz_yyy_1, \
                                         tg_zzzz_yyyy_0, tg_zzzz_yyyy_1, tg_zzzz_yyyz_0, tg_zzzz_yyyz_1, tg_zzzz_yyz_1, \
                                         tg_zzzz_yyzz_0, tg_zzzz_yyzz_1, tg_zzzz_yzz_1, tg_zzzz_yzzz_0, tg_zzzz_yzzz_1, \
                                         tg_zzzz_zzz_1, tg_zzzz_zzzz_0, tg_zzzz_zzzz_1, tg_zzzzz_xxxx_0, tg_zzzzz_xxxy_0, \
                                         tg_zzzzz_xxxz_0, tg_zzzzz_xxyy_0, tg_zzzzz_xxyz_0, tg_zzzzz_xxzz_0, tg_zzzzz_xyyy_0, \
                                         tg_zzzzz_xyyz_0, tg_zzzzz_xyzz_0, tg_zzzzz_xzzz_0, tg_zzzzz_yyyy_0, tg_zzzzz_yyyz_0, \
                                         tg_zzzzz_yyzz_0, tg_zzzzz_yzzz_0, tg_zzzzz_zzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyyy_yyzz_0[j] = pb_y * tg_yyyy_yyzz_0[j] + wp_y[j] * tg_yyyy_yyzz_1[j] + 2.0 * fl1_fx * tg_yyy_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyzz_1[j] + fl1_fxn * tg_yyyy_yzz_1[j];

                    tg_yyyyy_yzzz_0[j] = pb_y * tg_yyyy_yzzz_0[j] + wp_y[j] * tg_yyyy_yzzz_1[j] + 2.0 * fl1_fx * tg_yyy_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zzz_1[j];

                    tg_yyyyy_zzzz_0[j] = pb_y * tg_yyyy_zzzz_0[j] + wp_y[j] * tg_yyyy_zzzz_1[j] + 2.0 * fl1_fx * tg_yyy_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_zzzz_1[j];

                    tg_yyyyz_xxxx_0[j] = pb_y * tg_yyyz_xxxx_0[j] + wp_y[j] * tg_yyyz_xxxx_1[j] + 1.5 * fl1_fx * tg_yyz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxx_1[j];

                    tg_yyyyz_xxxy_0[j] = pb_y * tg_yyyz_xxxy_0[j] + wp_y[j] * tg_yyyz_xxxy_1[j] + 1.5 * fl1_fx * tg_yyz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyyz_xxx_1[j];

                    tg_yyyyz_xxxz_0[j] = pb_y * tg_yyyz_xxxz_0[j] + wp_y[j] * tg_yyyz_xxxz_1[j] + 1.5 * fl1_fx * tg_yyz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxxz_1[j];

                    tg_yyyyz_xxyy_0[j] = pb_y * tg_yyyz_xxyy_0[j] + wp_y[j] * tg_yyyz_xxyy_1[j] + 1.5 * fl1_fx * tg_yyz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxyy_1[j] + fl1_fxn * tg_yyyz_xxy_1[j];

                    tg_yyyyz_xxyz_0[j] = pb_y * tg_yyyz_xxyz_0[j] + wp_y[j] * tg_yyyz_xxyz_1[j] + 1.5 * fl1_fx * tg_yyz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xxz_1[j];

                    tg_yyyyz_xxzz_0[j] = pb_y * tg_yyyz_xxzz_0[j] + wp_y[j] * tg_yyyz_xxzz_1[j] + 1.5 * fl1_fx * tg_yyz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxzz_1[j];

                    tg_yyyyz_xyyy_0[j] = pb_y * tg_yyyz_xyyy_0[j] + wp_y[j] * tg_yyyz_xyyy_1[j] + 1.5 * fl1_fx * tg_yyz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyyz_xyy_1[j];

                    tg_yyyyz_xyyz_0[j] = pb_y * tg_yyyz_xyyz_0[j] + wp_y[j] * tg_yyyz_xyyz_1[j] + 1.5 * fl1_fx * tg_yyz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyyz_1[j] + fl1_fxn * tg_yyyz_xyz_1[j];

                    tg_yyyyz_xyzz_0[j] = pb_y * tg_yyyz_xyzz_0[j] + wp_y[j] * tg_yyyz_xyzz_1[j] + 1.5 * fl1_fx * tg_yyz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xzz_1[j];

                    tg_yyyyz_xzzz_0[j] = pb_y * tg_yyyz_xzzz_0[j] + wp_y[j] * tg_yyyz_xzzz_1[j] + 1.5 * fl1_fx * tg_yyz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xzzz_1[j];

                    tg_yyyyz_yyyy_0[j] = pb_y * tg_yyyz_yyyy_0[j] + wp_y[j] * tg_yyyz_yyyy_1[j] + 1.5 * fl1_fx * tg_yyz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyyz_yyy_1[j];

                    tg_yyyyz_yyyz_0[j] = pb_y * tg_yyyz_yyyz_0[j] + wp_y[j] * tg_yyyz_yyyz_1[j] + 1.5 * fl1_fx * tg_yyz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyyz_yyz_1[j];

                    tg_yyyyz_yyzz_0[j] = pb_y * tg_yyyz_yyzz_0[j] + wp_y[j] * tg_yyyz_yyzz_1[j] + 1.5 * fl1_fx * tg_yyz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyzz_1[j] + fl1_fxn * tg_yyyz_yzz_1[j];

                    tg_yyyyz_yzzz_0[j] = pb_y * tg_yyyz_yzzz_0[j] + wp_y[j] * tg_yyyz_yzzz_1[j] + 1.5 * fl1_fx * tg_yyz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zzz_1[j];

                    tg_yyyyz_zzzz_0[j] = pb_y * tg_yyyz_zzzz_0[j] + wp_y[j] * tg_yyyz_zzzz_1[j] + 1.5 * fl1_fx * tg_yyz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_zzzz_1[j];

                    tg_yyyzz_xxxx_0[j] = pb_y * tg_yyzz_xxxx_0[j] + wp_y[j] * tg_yyzz_xxxx_1[j] + fl1_fx * tg_yzz_xxxx_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxx_1[j];

                    tg_yyyzz_xxxy_0[j] = pb_y * tg_yyzz_xxxy_0[j] + wp_y[j] * tg_yyzz_xxxy_1[j] + fl1_fx * tg_yzz_xxxy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyzz_xxx_1[j];

                    tg_yyyzz_xxxz_0[j] = pb_y * tg_yyzz_xxxz_0[j] + wp_y[j] * tg_yyzz_xxxz_1[j] + fl1_fx * tg_yzz_xxxz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxxz_1[j];

                    tg_yyyzz_xxyy_0[j] = pb_y * tg_yyzz_xxyy_0[j] + wp_y[j] * tg_yyzz_xxyy_1[j] + fl1_fx * tg_yzz_xxyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxyy_1[j] + fl1_fxn * tg_yyzz_xxy_1[j];

                    tg_yyyzz_xxyz_0[j] = pb_y * tg_yyzz_xxyz_0[j] + wp_y[j] * tg_yyzz_xxyz_1[j] + fl1_fx * tg_yzz_xxyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xxz_1[j];

                    tg_yyyzz_xxzz_0[j] = pb_y * tg_yyzz_xxzz_0[j] + wp_y[j] * tg_yyzz_xxzz_1[j] + fl1_fx * tg_yzz_xxzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxzz_1[j];

                    tg_yyyzz_xyyy_0[j] = pb_y * tg_yyzz_xyyy_0[j] + wp_y[j] * tg_yyzz_xyyy_1[j] + fl1_fx * tg_yzz_xyyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyzz_xyy_1[j];

                    tg_yyyzz_xyyz_0[j] = pb_y * tg_yyzz_xyyz_0[j] + wp_y[j] * tg_yyzz_xyyz_1[j] + fl1_fx * tg_yzz_xyyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyyz_1[j] + fl1_fxn * tg_yyzz_xyz_1[j];

                    tg_yyyzz_xyzz_0[j] = pb_y * tg_yyzz_xyzz_0[j] + wp_y[j] * tg_yyzz_xyzz_1[j] + fl1_fx * tg_yzz_xyzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xzz_1[j];

                    tg_yyyzz_xzzz_0[j] = pb_y * tg_yyzz_xzzz_0[j] + wp_y[j] * tg_yyzz_xzzz_1[j] + fl1_fx * tg_yzz_xzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xzzz_1[j];

                    tg_yyyzz_yyyy_0[j] = pb_y * tg_yyzz_yyyy_0[j] + wp_y[j] * tg_yyzz_yyyy_1[j] + fl1_fx * tg_yzz_yyyy_0[j] - fl1_fx * fl1_fza * tg_yzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyzz_yyy_1[j];

                    tg_yyyzz_yyyz_0[j] = pb_y * tg_yyzz_yyyz_0[j] + wp_y[j] * tg_yyzz_yyyz_1[j] + fl1_fx * tg_yzz_yyyz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyzz_yyz_1[j];

                    tg_yyyzz_yyzz_0[j] = pb_y * tg_yyzz_yyzz_0[j] + wp_y[j] * tg_yyzz_yyzz_1[j] + fl1_fx * tg_yzz_yyzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyzz_1[j] + fl1_fxn * tg_yyzz_yzz_1[j];

                    tg_yyyzz_yzzz_0[j] = pb_y * tg_yyzz_yzzz_0[j] + wp_y[j] * tg_yyzz_yzzz_1[j] + fl1_fx * tg_yzz_yzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zzz_1[j];

                    tg_yyyzz_zzzz_0[j] = pb_y * tg_yyzz_zzzz_0[j] + wp_y[j] * tg_yyzz_zzzz_1[j] + fl1_fx * tg_yzz_zzzz_0[j] - fl1_fx * fl1_fza * tg_yzz_zzzz_1[j];

                    tg_yyzzz_xxxx_0[j] = pb_y * tg_yzzz_xxxx_0[j] + wp_y[j] * tg_yzzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zzz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxx_1[j];

                    tg_yyzzz_xxxy_0[j] = pb_y * tg_yzzz_xxxy_0[j] + wp_y[j] * tg_yzzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zzz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yzzz_xxx_1[j];

                    tg_yyzzz_xxxz_0[j] = pb_y * tg_yzzz_xxxz_0[j] + wp_y[j] * tg_yzzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zzz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxxz_1[j];

                    tg_yyzzz_xxyy_0[j] = pb_y * tg_yzzz_xxyy_0[j] + wp_y[j] * tg_yzzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zzz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyy_1[j] + fl1_fxn * tg_yzzz_xxy_1[j];

                    tg_yyzzz_xxyz_0[j] = pb_y * tg_yzzz_xxyz_0[j] + wp_y[j] * tg_yzzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zzz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xxz_1[j];

                    tg_yyzzz_xxzz_0[j] = pb_y * tg_yzzz_xxzz_0[j] + wp_y[j] * tg_yzzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zzz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxzz_1[j];

                    tg_yyzzz_xyyy_0[j] = pb_y * tg_yzzz_xyyy_0[j] + wp_y[j] * tg_yzzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zzz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yzzz_xyy_1[j];

                    tg_yyzzz_xyyz_0[j] = pb_y * tg_yzzz_xyyz_0[j] + wp_y[j] * tg_yzzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zzz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyyz_1[j] + fl1_fxn * tg_yzzz_xyz_1[j];

                    tg_yyzzz_xyzz_0[j] = pb_y * tg_yzzz_xyzz_0[j] + wp_y[j] * tg_yzzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zzz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xzz_1[j];

                    tg_yyzzz_xzzz_0[j] = pb_y * tg_yzzz_xzzz_0[j] + wp_y[j] * tg_yzzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zzz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xzzz_1[j];

                    tg_yyzzz_yyyy_0[j] = pb_y * tg_yzzz_yyyy_0[j] + wp_y[j] * tg_yzzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zzz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yzzz_yyy_1[j];

                    tg_yyzzz_yyyz_0[j] = pb_y * tg_yzzz_yyyz_0[j] + wp_y[j] * tg_yzzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zzz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yzzz_yyz_1[j];

                    tg_yyzzz_yyzz_0[j] = pb_y * tg_yzzz_yyzz_0[j] + wp_y[j] * tg_yzzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zzz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyzz_1[j] + fl1_fxn * tg_yzzz_yzz_1[j];

                    tg_yyzzz_yzzz_0[j] = pb_y * tg_yzzz_yzzz_0[j] + wp_y[j] * tg_yzzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zzz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zzz_1[j];

                    tg_yyzzz_zzzz_0[j] = pb_y * tg_yzzz_zzzz_0[j] + wp_y[j] * tg_yzzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zzz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_zzzz_1[j];

                    tg_yzzzz_xxxx_0[j] = pb_y * tg_zzzz_xxxx_0[j] + wp_y[j] * tg_zzzz_xxxx_1[j];

                    tg_yzzzz_xxxy_0[j] = pb_y * tg_zzzz_xxxy_0[j] + wp_y[j] * tg_zzzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxx_1[j];

                    tg_yzzzz_xxxz_0[j] = pb_y * tg_zzzz_xxxz_0[j] + wp_y[j] * tg_zzzz_xxxz_1[j];

                    tg_yzzzz_xxyy_0[j] = pb_y * tg_zzzz_xxyy_0[j] + wp_y[j] * tg_zzzz_xxyy_1[j] + fl1_fxn * tg_zzzz_xxy_1[j];

                    tg_yzzzz_xxyz_0[j] = pb_y * tg_zzzz_xxyz_0[j] + wp_y[j] * tg_zzzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxz_1[j];

                    tg_yzzzz_xxzz_0[j] = pb_y * tg_zzzz_xxzz_0[j] + wp_y[j] * tg_zzzz_xxzz_1[j];

                    tg_yzzzz_xyyy_0[j] = pb_y * tg_zzzz_xyyy_0[j] + wp_y[j] * tg_zzzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_xyy_1[j];

                    tg_yzzzz_xyyz_0[j] = pb_y * tg_zzzz_xyyz_0[j] + wp_y[j] * tg_zzzz_xyyz_1[j] + fl1_fxn * tg_zzzz_xyz_1[j];

                    tg_yzzzz_xyzz_0[j] = pb_y * tg_zzzz_xyzz_0[j] + wp_y[j] * tg_zzzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xzz_1[j];

                    tg_yzzzz_xzzz_0[j] = pb_y * tg_zzzz_xzzz_0[j] + wp_y[j] * tg_zzzz_xzzz_1[j];

                    tg_yzzzz_yyyy_0[j] = pb_y * tg_zzzz_yyyy_0[j] + wp_y[j] * tg_zzzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_zzzz_yyy_1[j];

                    tg_yzzzz_yyyz_0[j] = pb_y * tg_zzzz_yyyz_0[j] + wp_y[j] * tg_zzzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_zzzz_yyz_1[j];

                    tg_yzzzz_yyzz_0[j] = pb_y * tg_zzzz_yyzz_0[j] + wp_y[j] * tg_zzzz_yyzz_1[j] + fl1_fxn * tg_zzzz_yzz_1[j];

                    tg_yzzzz_yzzz_0[j] = pb_y * tg_zzzz_yzzz_0[j] + wp_y[j] * tg_zzzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zzz_1[j];

                    tg_yzzzz_zzzz_0[j] = pb_y * tg_zzzz_zzzz_0[j] + wp_y[j] * tg_zzzz_zzzz_1[j];

                    tg_zzzzz_xxxx_0[j] = pb_z * tg_zzzz_xxxx_0[j] + wp_z[j] * tg_zzzz_xxxx_1[j] + 2.0 * fl1_fx * tg_zzz_xxxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxx_1[j];

                    tg_zzzzz_xxxy_0[j] = pb_z * tg_zzzz_xxxy_0[j] + wp_z[j] * tg_zzzz_xxxy_1[j] + 2.0 * fl1_fx * tg_zzz_xxxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxy_1[j];

                    tg_zzzzz_xxxz_0[j] = pb_z * tg_zzzz_xxxz_0[j] + wp_z[j] * tg_zzzz_xxxz_1[j] + 2.0 * fl1_fx * tg_zzz_xxxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxxz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxx_1[j];

                    tg_zzzzz_xxyy_0[j] = pb_z * tg_zzzz_xxyy_0[j] + wp_z[j] * tg_zzzz_xxyy_1[j] + 2.0 * fl1_fx * tg_zzz_xxyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxyy_1[j];

                    tg_zzzzz_xxyz_0[j] = pb_z * tg_zzzz_xxyz_0[j] + wp_z[j] * tg_zzzz_xxyz_1[j] + 2.0 * fl1_fx * tg_zzz_xxyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xxy_1[j];

                    tg_zzzzz_xxzz_0[j] = pb_z * tg_zzzz_xxzz_0[j] + wp_z[j] * tg_zzzz_xxzz_1[j] + 2.0 * fl1_fx * tg_zzz_xxzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxzz_1[j] + fl1_fxn * tg_zzzz_xxz_1[j];

                    tg_zzzzz_xyyy_0[j] = pb_z * tg_zzzz_xyyy_0[j] + wp_z[j] * tg_zzzz_xyyy_1[j] + 2.0 * fl1_fx * tg_zzz_xyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyyy_1[j];

                    tg_zzzzz_xyyz_0[j] = pb_z * tg_zzzz_xyyz_0[j] + wp_z[j] * tg_zzzz_xyyz_1[j] + 2.0 * fl1_fx * tg_zzz_xyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xyy_1[j];

                    tg_zzzzz_xyzz_0[j] = pb_z * tg_zzzz_xyzz_0[j] + wp_z[j] * tg_zzzz_xyzz_1[j] + 2.0 * fl1_fx * tg_zzz_xyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyzz_1[j] + fl1_fxn * tg_zzzz_xyz_1[j];

                    tg_zzzzz_xzzz_0[j] = pb_z * tg_zzzz_xzzz_0[j] + wp_z[j] * tg_zzzz_xzzz_1[j] + 2.0 * fl1_fx * tg_zzz_xzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_xzz_1[j];

                    tg_zzzzz_yyyy_0[j] = pb_z * tg_zzzz_yyyy_0[j] + wp_z[j] * tg_zzzz_yyyy_1[j] + 2.0 * fl1_fx * tg_zzz_yyyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyyy_1[j];

                    tg_zzzzz_yyyz_0[j] = pb_z * tg_zzzz_yyyz_0[j] + wp_z[j] * tg_zzzz_yyyz_1[j] + 2.0 * fl1_fx * tg_zzz_yyyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yyy_1[j];

                    tg_zzzzz_yyzz_0[j] = pb_z * tg_zzzz_yyzz_0[j] + wp_z[j] * tg_zzzz_yyzz_1[j] + 2.0 * fl1_fx * tg_zzz_yyzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyzz_1[j] + fl1_fxn * tg_zzzz_yyz_1[j];

                    tg_zzzzz_yzzz_0[j] = pb_z * tg_zzzz_yzzz_0[j] + wp_z[j] * tg_zzzz_yzzz_1[j] + 2.0 * fl1_fx * tg_zzz_yzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yzzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_yzz_1[j];

                    tg_zzzzz_zzzz_0[j] = pb_z * tg_zzzz_zzzz_0[j] + wp_z[j] * tg_zzzz_zzzz_1[j] + 2.0 * fl1_fx * tg_zzz_zzzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_zzzz_1[j] + 2.0 * fl1_fxn * tg_zzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

