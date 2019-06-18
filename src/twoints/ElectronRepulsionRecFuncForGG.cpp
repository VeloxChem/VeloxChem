//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForGG.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSGSG(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSGSG_0_75(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSGSG_75_150(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSGSG_150_225(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSGSG_0_75(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,75)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xx_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx); 

                auto tg_xx_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 1); 

                auto tg_xx_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 2); 

                auto tg_xx_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 3); 

                auto tg_xx_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 4); 

                auto tg_xx_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 5); 

                auto tg_xx_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 6); 

                auto tg_xx_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 7); 

                auto tg_xx_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 8); 

                auto tg_xx_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 9); 

                auto tg_xx_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 10); 

                auto tg_xx_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 11); 

                auto tg_xx_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 12); 

                auto tg_xx_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 13); 

                auto tg_xx_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 14); 

                auto tg_xy_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 15); 

                auto tg_xy_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 16); 

                auto tg_xy_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 17); 

                auto tg_xy_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 18); 

                auto tg_xy_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 19); 

                auto tg_xy_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 20); 

                auto tg_xy_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 21); 

                auto tg_xy_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 22); 

                auto tg_xy_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 23); 

                auto tg_xy_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 24); 

                auto tg_xy_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 25); 

                auto tg_xy_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 26); 

                auto tg_xy_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 27); 

                auto tg_xy_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 28); 

                auto tg_xy_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 29); 

                auto tg_xz_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 30); 

                auto tg_xz_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 31); 

                auto tg_xz_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 32); 

                auto tg_xz_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 33); 

                auto tg_xz_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 34); 

                auto tg_xz_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 35); 

                auto tg_xz_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 36); 

                auto tg_xz_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 37); 

                auto tg_xz_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 38); 

                auto tg_xz_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 39); 

                auto tg_xz_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 40); 

                auto tg_xz_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 41); 

                auto tg_xz_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 42); 

                auto tg_xz_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 43); 

                auto tg_xz_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 44); 

                auto tg_yy_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 45); 

                auto tg_yy_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 46); 

                auto tg_yy_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 47); 

                auto tg_yy_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 48); 

                auto tg_yy_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 49); 

                auto tg_yy_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 50); 

                auto tg_yy_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 51); 

                auto tg_yy_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 52); 

                auto tg_yy_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 53); 

                auto tg_yy_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 54); 

                auto tg_yy_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 55); 

                auto tg_yy_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 56); 

                auto tg_yy_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 57); 

                auto tg_yy_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 58); 

                auto tg_yy_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 59); 

                auto tg_yz_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 60); 

                auto tg_yz_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 61); 

                auto tg_yz_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 62); 

                auto tg_yz_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 63); 

                auto tg_yz_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 64); 

                auto tg_yz_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 65); 

                auto tg_yz_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 66); 

                auto tg_yz_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 67); 

                auto tg_yz_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 68); 

                auto tg_yz_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 69); 

                auto tg_yz_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 70); 

                auto tg_yz_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 71); 

                auto tg_yz_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 72); 

                auto tg_yz_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 73); 

                auto tg_yz_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 74); 

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

                auto tg_xxx_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx); 

                auto tg_xxx_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 1); 

                auto tg_xxx_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 2); 

                auto tg_xxx_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 3); 

                auto tg_xxx_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 4); 

                auto tg_xxx_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 5); 

                auto tg_xxx_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 6); 

                auto tg_xxx_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 7); 

                auto tg_xxx_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 8); 

                auto tg_xxx_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 9); 

                auto tg_xxy_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 10); 

                auto tg_xxy_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 11); 

                auto tg_xxy_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 12); 

                auto tg_xxy_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 13); 

                auto tg_xxy_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 14); 

                auto tg_xxy_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 15); 

                auto tg_xxy_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 16); 

                auto tg_xxy_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 17); 

                auto tg_xxy_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 18); 

                auto tg_xxy_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 19); 

                auto tg_xxz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 20); 

                auto tg_xxz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 21); 

                auto tg_xxz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 22); 

                auto tg_xxz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 23); 

                auto tg_xxz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 24); 

                auto tg_xxz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 25); 

                auto tg_xxz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 26); 

                auto tg_xxz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 27); 

                auto tg_xxz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 28); 

                auto tg_xxz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 29); 

                auto tg_xyy_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 30); 

                auto tg_xyy_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 31); 

                auto tg_xyy_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 32); 

                auto tg_xyy_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 33); 

                auto tg_xyy_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 34); 

                auto tg_xyy_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 35); 

                auto tg_xyy_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 36); 

                auto tg_xyy_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 37); 

                auto tg_xyy_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 38); 

                auto tg_xyy_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 39); 

                auto tg_xyz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 40); 

                auto tg_xyz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 41); 

                auto tg_xyz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 42); 

                auto tg_xyz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 43); 

                auto tg_xyz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 44); 

                auto tg_xyz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 45); 

                auto tg_xyz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 46); 

                auto tg_xyz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 47); 

                auto tg_xyz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 48); 

                auto tg_xyz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 49); 

                // set up pointers to integrals

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

                // Batch of Integrals (0,75)

                #pragma omp simd aligned(fxn, fza, tg_xx_xxxx_0, tg_xx_xxxx_1, tg_xx_xxxy_0, tg_xx_xxxy_1, \
                                         tg_xx_xxxz_0, tg_xx_xxxz_1, tg_xx_xxyy_0, tg_xx_xxyy_1, tg_xx_xxyz_0, tg_xx_xxyz_1, \
                                         tg_xx_xxzz_0, tg_xx_xxzz_1, tg_xx_xyyy_0, tg_xx_xyyy_1, tg_xx_xyyz_0, tg_xx_xyyz_1, \
                                         tg_xx_xyzz_0, tg_xx_xyzz_1, tg_xx_xzzz_0, tg_xx_xzzz_1, tg_xx_yyyy_0, tg_xx_yyyy_1, \
                                         tg_xx_yyyz_0, tg_xx_yyyz_1, tg_xx_yyzz_0, tg_xx_yyzz_1, tg_xx_yzzz_0, tg_xx_yzzz_1, \
                                         tg_xx_zzzz_0, tg_xx_zzzz_1, tg_xxx_xxx_1, tg_xxx_xxxx_0, tg_xxx_xxxx_1, \
                                         tg_xxx_xxxy_0, tg_xxx_xxxy_1, tg_xxx_xxxz_0, tg_xxx_xxxz_1, tg_xxx_xxy_1, \
                                         tg_xxx_xxyy_0, tg_xxx_xxyy_1, tg_xxx_xxyz_0, tg_xxx_xxyz_1, tg_xxx_xxz_1, \
                                         tg_xxx_xxzz_0, tg_xxx_xxzz_1, tg_xxx_xyy_1, tg_xxx_xyyy_0, tg_xxx_xyyy_1, \
                                         tg_xxx_xyyz_0, tg_xxx_xyyz_1, tg_xxx_xyz_1, tg_xxx_xyzz_0, tg_xxx_xyzz_1, \
                                         tg_xxx_xzz_1, tg_xxx_xzzz_0, tg_xxx_xzzz_1, tg_xxx_yyy_1, tg_xxx_yyyy_0, \
                                         tg_xxx_yyyy_1, tg_xxx_yyyz_0, tg_xxx_yyyz_1, tg_xxx_yyz_1, tg_xxx_yyzz_0, \
                                         tg_xxx_yyzz_1, tg_xxx_yzz_1, tg_xxx_yzzz_0, tg_xxx_yzzz_1, tg_xxx_zzz_1, \
                                         tg_xxx_zzzz_0, tg_xxx_zzzz_1, tg_xxxx_xxxx_0, tg_xxxx_xxxy_0, tg_xxxx_xxxz_0, \
                                         tg_xxxx_xxyy_0, tg_xxxx_xxyz_0, tg_xxxx_xxzz_0, tg_xxxx_xyyy_0, tg_xxxx_xyyz_0, \
                                         tg_xxxx_xyzz_0, tg_xxxx_xzzz_0, tg_xxxx_yyyy_0, tg_xxxx_yyyz_0, tg_xxxx_yyzz_0, \
                                         tg_xxxx_yzzz_0, tg_xxxx_zzzz_0, tg_xxxy_xxxx_0, tg_xxxy_xxxy_0, tg_xxxy_xxxz_0, \
                                         tg_xxxy_xxyy_0, tg_xxxy_xxyz_0, tg_xxxy_xxzz_0, tg_xxxy_xyyy_0, tg_xxxy_xyyz_0, \
                                         tg_xxxy_xyzz_0, tg_xxxy_xzzz_0, tg_xxxy_yyyy_0, tg_xxxy_yyyz_0, tg_xxxy_yyzz_0, \
                                         tg_xxxy_yzzz_0, tg_xxxy_zzzz_0, tg_xxxz_xxxx_0, tg_xxxz_xxxy_0, tg_xxxz_xxxz_0, \
                                         tg_xxxz_xxyy_0, tg_xxxz_xxyz_0, tg_xxxz_xxzz_0, tg_xxxz_xyyy_0, tg_xxxz_xyyz_0, \
                                         tg_xxxz_xyzz_0, tg_xxxz_xzzz_0, tg_xxxz_yyyy_0, tg_xxxz_yyyz_0, tg_xxxz_yyzz_0, \
                                         tg_xxxz_yzzz_0, tg_xxxz_zzzz_0, tg_xxy_xxx_1, tg_xxy_xxxx_0, tg_xxy_xxxx_1, \
                                         tg_xxy_xxxy_0, tg_xxy_xxxy_1, tg_xxy_xxxz_0, tg_xxy_xxxz_1, tg_xxy_xxy_1, \
                                         tg_xxy_xxyy_0, tg_xxy_xxyy_1, tg_xxy_xxyz_0, tg_xxy_xxyz_1, tg_xxy_xxz_1, \
                                         tg_xxy_xxzz_0, tg_xxy_xxzz_1, tg_xxy_xyy_1, tg_xxy_xyyy_0, tg_xxy_xyyy_1, \
                                         tg_xxy_xyyz_0, tg_xxy_xyyz_1, tg_xxy_xyz_1, tg_xxy_xyzz_0, tg_xxy_xyzz_1, \
                                         tg_xxy_xzz_1, tg_xxy_xzzz_0, tg_xxy_xzzz_1, tg_xxy_yyy_1, tg_xxy_yyyy_0, \
                                         tg_xxy_yyyy_1, tg_xxy_yyyz_0, tg_xxy_yyyz_1, tg_xxy_yyz_1, tg_xxy_yyzz_0, \
                                         tg_xxy_yyzz_1, tg_xxy_yzz_1, tg_xxy_yzzz_0, tg_xxy_yzzz_1, tg_xxy_zzz_1, \
                                         tg_xxy_zzzz_0, tg_xxy_zzzz_1, tg_xxyy_xxxx_0, tg_xxyy_xxxy_0, tg_xxyy_xxxz_0, \
                                         tg_xxyy_xxyy_0, tg_xxyy_xxyz_0, tg_xxyy_xxzz_0, tg_xxyy_xyyy_0, tg_xxyy_xyyz_0, \
                                         tg_xxyy_xyzz_0, tg_xxyy_xzzz_0, tg_xxyy_yyyy_0, tg_xxyy_yyyz_0, tg_xxyy_yyzz_0, \
                                         tg_xxyy_yzzz_0, tg_xxyy_zzzz_0, tg_xxyz_xxxx_0, tg_xxyz_xxxy_0, tg_xxyz_xxxz_0, \
                                         tg_xxyz_xxyy_0, tg_xxyz_xxyz_0, tg_xxyz_xxzz_0, tg_xxyz_xyyy_0, tg_xxyz_xyyz_0, \
                                         tg_xxyz_xyzz_0, tg_xxyz_xzzz_0, tg_xxyz_yyyy_0, tg_xxyz_yyyz_0, tg_xxyz_yyzz_0, \
                                         tg_xxyz_yzzz_0, tg_xxyz_zzzz_0, tg_xxz_xxx_1, tg_xxz_xxxx_0, tg_xxz_xxxx_1, \
                                         tg_xxz_xxxy_0, tg_xxz_xxxy_1, tg_xxz_xxxz_0, tg_xxz_xxxz_1, tg_xxz_xxy_1, \
                                         tg_xxz_xxyy_0, tg_xxz_xxyy_1, tg_xxz_xxyz_0, tg_xxz_xxyz_1, tg_xxz_xxz_1, \
                                         tg_xxz_xxzz_0, tg_xxz_xxzz_1, tg_xxz_xyy_1, tg_xxz_xyyy_0, tg_xxz_xyyy_1, \
                                         tg_xxz_xyyz_0, tg_xxz_xyyz_1, tg_xxz_xyz_1, tg_xxz_xyzz_0, tg_xxz_xyzz_1, \
                                         tg_xxz_xzz_1, tg_xxz_xzzz_0, tg_xxz_xzzz_1, tg_xxz_yyy_1, tg_xxz_yyyy_0, \
                                         tg_xxz_yyyy_1, tg_xxz_yyyz_0, tg_xxz_yyyz_1, tg_xxz_yyz_1, tg_xxz_yyzz_0, \
                                         tg_xxz_yyzz_1, tg_xxz_yzz_1, tg_xxz_yzzz_0, tg_xxz_yzzz_1, tg_xxz_zzz_1, \
                                         tg_xxz_zzzz_0, tg_xxz_zzzz_1, tg_xy_xxxx_0, tg_xy_xxxx_1, tg_xy_xxxy_0, \
                                         tg_xy_xxxy_1, tg_xy_xxxz_0, tg_xy_xxxz_1, tg_xy_xxyy_0, tg_xy_xxyy_1, tg_xy_xxyz_0, \
                                         tg_xy_xxyz_1, tg_xy_xxzz_0, tg_xy_xxzz_1, tg_xy_xyyy_0, tg_xy_xyyy_1, tg_xy_xyyz_0, \
                                         tg_xy_xyyz_1, tg_xy_xyzz_0, tg_xy_xyzz_1, tg_xy_xzzz_0, tg_xy_xzzz_1, tg_xy_yyyy_0, \
                                         tg_xy_yyyy_1, tg_xy_yyyz_0, tg_xy_yyyz_1, tg_xy_yyzz_0, tg_xy_yyzz_1, tg_xy_yzzz_0, \
                                         tg_xy_yzzz_1, tg_xy_zzzz_0, tg_xy_zzzz_1, tg_xyy_xxx_1, tg_xyy_xxxx_0, \
                                         tg_xyy_xxxx_1, tg_xyy_xxxy_0, tg_xyy_xxxy_1, tg_xyy_xxxz_0, tg_xyy_xxxz_1, \
                                         tg_xyy_xxy_1, tg_xyy_xxyy_0, tg_xyy_xxyy_1, tg_xyy_xxyz_0, tg_xyy_xxyz_1, \
                                         tg_xyy_xxz_1, tg_xyy_xxzz_0, tg_xyy_xxzz_1, tg_xyy_xyy_1, tg_xyy_xyyy_0, \
                                         tg_xyy_xyyy_1, tg_xyy_xyyz_0, tg_xyy_xyyz_1, tg_xyy_xyz_1, tg_xyy_xyzz_0, \
                                         tg_xyy_xyzz_1, tg_xyy_xzz_1, tg_xyy_xzzz_0, tg_xyy_xzzz_1, tg_xyy_yyy_1, \
                                         tg_xyy_yyyy_0, tg_xyy_yyyy_1, tg_xyy_yyyz_0, tg_xyy_yyyz_1, tg_xyy_yyz_1, \
                                         tg_xyy_yyzz_0, tg_xyy_yyzz_1, tg_xyy_yzz_1, tg_xyy_yzzz_0, tg_xyy_yzzz_1, \
                                         tg_xyy_zzz_1, tg_xyy_zzzz_0, tg_xyy_zzzz_1, tg_xyz_xxx_1, tg_xyz_xxxx_0, \
                                         tg_xyz_xxxx_1, tg_xyz_xxxy_0, tg_xyz_xxxy_1, tg_xyz_xxxz_0, tg_xyz_xxxz_1, \
                                         tg_xyz_xxy_1, tg_xyz_xxyy_0, tg_xyz_xxyy_1, tg_xyz_xxyz_0, tg_xyz_xxyz_1, \
                                         tg_xyz_xxz_1, tg_xyz_xxzz_0, tg_xyz_xxzz_1, tg_xyz_xyy_1, tg_xyz_xyyy_0, \
                                         tg_xyz_xyyy_1, tg_xyz_xyyz_0, tg_xyz_xyyz_1, tg_xyz_xyz_1, tg_xyz_xyzz_0, \
                                         tg_xyz_xyzz_1, tg_xyz_xzz_1, tg_xyz_xzzz_0, tg_xyz_xzzz_1, tg_xyz_yyy_1, \
                                         tg_xyz_yyyy_0, tg_xyz_yyyy_1, tg_xyz_yyyz_0, tg_xyz_yyyz_1, tg_xyz_yyz_1, \
                                         tg_xyz_yyzz_0, tg_xyz_yyzz_1, tg_xyz_yzz_1, tg_xyz_yzzz_0, tg_xyz_yzzz_1, \
                                         tg_xyz_zzz_1, tg_xyz_zzzz_0, tg_xyz_zzzz_1, tg_xz_xxxx_0, tg_xz_xxxx_1, \
                                         tg_xz_xxxy_0, tg_xz_xxxy_1, tg_xz_xxxz_0, tg_xz_xxxz_1, tg_xz_xxyy_0, tg_xz_xxyy_1, \
                                         tg_xz_xxyz_0, tg_xz_xxyz_1, tg_xz_xxzz_0, tg_xz_xxzz_1, tg_xz_xyyy_0, tg_xz_xyyy_1, \
                                         tg_xz_xyyz_0, tg_xz_xyyz_1, tg_xz_xyzz_0, tg_xz_xyzz_1, tg_xz_xzzz_0, tg_xz_xzzz_1, \
                                         tg_xz_yyyy_0, tg_xz_yyyy_1, tg_xz_yyyz_0, tg_xz_yyyz_1, tg_xz_yyzz_0, tg_xz_yyzz_1, \
                                         tg_xz_yzzz_0, tg_xz_yzzz_1, tg_xz_zzzz_0, tg_xz_zzzz_1, tg_yy_xxxx_0, tg_yy_xxxx_1, \
                                         tg_yy_xxxy_0, tg_yy_xxxy_1, tg_yy_xxxz_0, tg_yy_xxxz_1, tg_yy_xxyy_0, tg_yy_xxyy_1, \
                                         tg_yy_xxyz_0, tg_yy_xxyz_1, tg_yy_xxzz_0, tg_yy_xxzz_1, tg_yy_xyyy_0, tg_yy_xyyy_1, \
                                         tg_yy_xyyz_0, tg_yy_xyyz_1, tg_yy_xyzz_0, tg_yy_xyzz_1, tg_yy_xzzz_0, tg_yy_xzzz_1, \
                                         tg_yy_yyyy_0, tg_yy_yyyy_1, tg_yy_yyyz_0, tg_yy_yyyz_1, tg_yy_yyzz_0, tg_yy_yyzz_1, \
                                         tg_yy_yzzz_0, tg_yy_yzzz_1, tg_yy_zzzz_0, tg_yy_zzzz_1, tg_yz_xxxx_0, tg_yz_xxxx_1, \
                                         tg_yz_xxxy_0, tg_yz_xxxy_1, tg_yz_xxxz_0, tg_yz_xxxz_1, tg_yz_xxyy_0, tg_yz_xxyy_1, \
                                         tg_yz_xxyz_0, tg_yz_xxyz_1, tg_yz_xxzz_0, tg_yz_xxzz_1, tg_yz_xyyy_0, tg_yz_xyyy_1, \
                                         tg_yz_xyyz_0, tg_yz_xyyz_1, tg_yz_xyzz_0, tg_yz_xyzz_1, tg_yz_xzzz_0, tg_yz_xzzz_1, \
                                         tg_yz_yyyy_0, tg_yz_yyyy_1, tg_yz_yyyz_0, tg_yz_yyyz_1, tg_yz_yyzz_0, tg_yz_yyzz_1, \
                                         tg_yz_yzzz_0, tg_yz_yzzz_1, tg_yz_zzzz_0, tg_yz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxx_xxxx_0[j] = pb_x * tg_xxx_xxxx_0[j] + wp_x[j] * tg_xxx_xxxx_1[j] + 1.5 * fl1_fx * tg_xx_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxx_xxx_1[j];

                    tg_xxxx_xxxy_0[j] = pb_x * tg_xxx_xxxy_0[j] + wp_x[j] * tg_xxx_xxxy_1[j] + 1.5 * fl1_fx * tg_xx_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxx_xxy_1[j];

                    tg_xxxx_xxxz_0[j] = pb_x * tg_xxx_xxxz_0[j] + wp_x[j] * tg_xxx_xxxz_1[j] + 1.5 * fl1_fx * tg_xx_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxx_xxz_1[j];

                    tg_xxxx_xxyy_0[j] = pb_x * tg_xxx_xxyy_0[j] + wp_x[j] * tg_xxx_xxyy_1[j] + 1.5 * fl1_fx * tg_xx_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyy_1[j] + fl1_fxn * tg_xxx_xyy_1[j];

                    tg_xxxx_xxyz_0[j] = pb_x * tg_xxx_xxyz_0[j] + wp_x[j] * tg_xxx_xxyz_1[j] + 1.5 * fl1_fx * tg_xx_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxyz_1[j] + fl1_fxn * tg_xxx_xyz_1[j];

                    tg_xxxx_xxzz_0[j] = pb_x * tg_xxx_xxzz_0[j] + wp_x[j] * tg_xxx_xxzz_1[j] + 1.5 * fl1_fx * tg_xx_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xxzz_1[j] + fl1_fxn * tg_xxx_xzz_1[j];

                    tg_xxxx_xyyy_0[j] = pb_x * tg_xxx_xyyy_0[j] + wp_x[j] * tg_xxx_xyyy_1[j] + 1.5 * fl1_fx * tg_xx_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxx_yyy_1[j];

                    tg_xxxx_xyyz_0[j] = pb_x * tg_xxx_xyyz_0[j] + wp_x[j] * tg_xxx_xyyz_1[j] + 1.5 * fl1_fx * tg_xx_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxx_yyz_1[j];

                    tg_xxxx_xyzz_0[j] = pb_x * tg_xxx_xyzz_0[j] + wp_x[j] * tg_xxx_xyzz_1[j] + 1.5 * fl1_fx * tg_xx_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxx_yzz_1[j];

                    tg_xxxx_xzzz_0[j] = pb_x * tg_xxx_xzzz_0[j] + wp_x[j] * tg_xxx_xzzz_1[j] + 1.5 * fl1_fx * tg_xx_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxx_zzz_1[j];

                    tg_xxxx_yyyy_0[j] = pb_x * tg_xxx_yyyy_0[j] + wp_x[j] * tg_xxx_yyyy_1[j] + 1.5 * fl1_fx * tg_xx_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyy_1[j];

                    tg_xxxx_yyyz_0[j] = pb_x * tg_xxx_yyyz_0[j] + wp_x[j] * tg_xxx_yyyz_1[j] + 1.5 * fl1_fx * tg_xx_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyyz_1[j];

                    tg_xxxx_yyzz_0[j] = pb_x * tg_xxx_yyzz_0[j] + wp_x[j] * tg_xxx_yyzz_1[j] + 1.5 * fl1_fx * tg_xx_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yyzz_1[j];

                    tg_xxxx_yzzz_0[j] = pb_x * tg_xxx_yzzz_0[j] + wp_x[j] * tg_xxx_yzzz_1[j] + 1.5 * fl1_fx * tg_xx_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_yzzz_1[j];

                    tg_xxxx_zzzz_0[j] = pb_x * tg_xxx_zzzz_0[j] + wp_x[j] * tg_xxx_zzzz_1[j] + 1.5 * fl1_fx * tg_xx_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_zzzz_1[j];

                    tg_xxxy_xxxx_0[j] = pb_x * tg_xxy_xxxx_0[j] + wp_x[j] * tg_xxy_xxxx_1[j] + fl1_fx * tg_xy_xxxx_0[j] - fl1_fx * fl1_fza * tg_xy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxy_xxx_1[j];

                    tg_xxxy_xxxy_0[j] = pb_x * tg_xxy_xxxy_0[j] + wp_x[j] * tg_xxy_xxxy_1[j] + fl1_fx * tg_xy_xxxy_0[j] - fl1_fx * fl1_fza * tg_xy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxy_xxy_1[j];

                    tg_xxxy_xxxz_0[j] = pb_x * tg_xxy_xxxz_0[j] + wp_x[j] * tg_xxy_xxxz_1[j] + fl1_fx * tg_xy_xxxz_0[j] - fl1_fx * fl1_fza * tg_xy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxy_xxz_1[j];

                    tg_xxxy_xxyy_0[j] = pb_x * tg_xxy_xxyy_0[j] + wp_x[j] * tg_xxy_xxyy_1[j] + fl1_fx * tg_xy_xxyy_0[j] - fl1_fx * fl1_fza * tg_xy_xxyy_1[j] + fl1_fxn * tg_xxy_xyy_1[j];

                    tg_xxxy_xxyz_0[j] = pb_x * tg_xxy_xxyz_0[j] + wp_x[j] * tg_xxy_xxyz_1[j] + fl1_fx * tg_xy_xxyz_0[j] - fl1_fx * fl1_fza * tg_xy_xxyz_1[j] + fl1_fxn * tg_xxy_xyz_1[j];

                    tg_xxxy_xxzz_0[j] = pb_x * tg_xxy_xxzz_0[j] + wp_x[j] * tg_xxy_xxzz_1[j] + fl1_fx * tg_xy_xxzz_0[j] - fl1_fx * fl1_fza * tg_xy_xxzz_1[j] + fl1_fxn * tg_xxy_xzz_1[j];

                    tg_xxxy_xyyy_0[j] = pb_x * tg_xxy_xyyy_0[j] + wp_x[j] * tg_xxy_xyyy_1[j] + fl1_fx * tg_xy_xyyy_0[j] - fl1_fx * fl1_fza * tg_xy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxy_yyy_1[j];

                    tg_xxxy_xyyz_0[j] = pb_x * tg_xxy_xyyz_0[j] + wp_x[j] * tg_xxy_xyyz_1[j] + fl1_fx * tg_xy_xyyz_0[j] - fl1_fx * fl1_fza * tg_xy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxy_yyz_1[j];

                    tg_xxxy_xyzz_0[j] = pb_x * tg_xxy_xyzz_0[j] + wp_x[j] * tg_xxy_xyzz_1[j] + fl1_fx * tg_xy_xyzz_0[j] - fl1_fx * fl1_fza * tg_xy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxy_yzz_1[j];

                    tg_xxxy_xzzz_0[j] = pb_x * tg_xxy_xzzz_0[j] + wp_x[j] * tg_xxy_xzzz_1[j] + fl1_fx * tg_xy_xzzz_0[j] - fl1_fx * fl1_fza * tg_xy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxy_zzz_1[j];

                    tg_xxxy_yyyy_0[j] = pb_x * tg_xxy_yyyy_0[j] + wp_x[j] * tg_xxy_yyyy_1[j] + fl1_fx * tg_xy_yyyy_0[j] - fl1_fx * fl1_fza * tg_xy_yyyy_1[j];

                    tg_xxxy_yyyz_0[j] = pb_x * tg_xxy_yyyz_0[j] + wp_x[j] * tg_xxy_yyyz_1[j] + fl1_fx * tg_xy_yyyz_0[j] - fl1_fx * fl1_fza * tg_xy_yyyz_1[j];

                    tg_xxxy_yyzz_0[j] = pb_x * tg_xxy_yyzz_0[j] + wp_x[j] * tg_xxy_yyzz_1[j] + fl1_fx * tg_xy_yyzz_0[j] - fl1_fx * fl1_fza * tg_xy_yyzz_1[j];

                    tg_xxxy_yzzz_0[j] = pb_x * tg_xxy_yzzz_0[j] + wp_x[j] * tg_xxy_yzzz_1[j] + fl1_fx * tg_xy_yzzz_0[j] - fl1_fx * fl1_fza * tg_xy_yzzz_1[j];

                    tg_xxxy_zzzz_0[j] = pb_x * tg_xxy_zzzz_0[j] + wp_x[j] * tg_xxy_zzzz_1[j] + fl1_fx * tg_xy_zzzz_0[j] - fl1_fx * fl1_fza * tg_xy_zzzz_1[j];

                    tg_xxxz_xxxx_0[j] = pb_x * tg_xxz_xxxx_0[j] + wp_x[j] * tg_xxz_xxxx_1[j] + fl1_fx * tg_xz_xxxx_0[j] - fl1_fx * fl1_fza * tg_xz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xxz_xxx_1[j];

                    tg_xxxz_xxxy_0[j] = pb_x * tg_xxz_xxxy_0[j] + wp_x[j] * tg_xxz_xxxy_1[j] + fl1_fx * tg_xz_xxxy_0[j] - fl1_fx * fl1_fza * tg_xz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xxz_xxy_1[j];

                    tg_xxxz_xxxz_0[j] = pb_x * tg_xxz_xxxz_0[j] + wp_x[j] * tg_xxz_xxxz_1[j] + fl1_fx * tg_xz_xxxz_0[j] - fl1_fx * fl1_fza * tg_xz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xxz_xxz_1[j];

                    tg_xxxz_xxyy_0[j] = pb_x * tg_xxz_xxyy_0[j] + wp_x[j] * tg_xxz_xxyy_1[j] + fl1_fx * tg_xz_xxyy_0[j] - fl1_fx * fl1_fza * tg_xz_xxyy_1[j] + fl1_fxn * tg_xxz_xyy_1[j];

                    tg_xxxz_xxyz_0[j] = pb_x * tg_xxz_xxyz_0[j] + wp_x[j] * tg_xxz_xxyz_1[j] + fl1_fx * tg_xz_xxyz_0[j] - fl1_fx * fl1_fza * tg_xz_xxyz_1[j] + fl1_fxn * tg_xxz_xyz_1[j];

                    tg_xxxz_xxzz_0[j] = pb_x * tg_xxz_xxzz_0[j] + wp_x[j] * tg_xxz_xxzz_1[j] + fl1_fx * tg_xz_xxzz_0[j] - fl1_fx * fl1_fza * tg_xz_xxzz_1[j] + fl1_fxn * tg_xxz_xzz_1[j];

                    tg_xxxz_xyyy_0[j] = pb_x * tg_xxz_xyyy_0[j] + wp_x[j] * tg_xxz_xyyy_1[j] + fl1_fx * tg_xz_xyyy_0[j] - fl1_fx * fl1_fza * tg_xz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xxz_yyy_1[j];

                    tg_xxxz_xyyz_0[j] = pb_x * tg_xxz_xyyz_0[j] + wp_x[j] * tg_xxz_xyyz_1[j] + fl1_fx * tg_xz_xyyz_0[j] - fl1_fx * fl1_fza * tg_xz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xxz_yyz_1[j];

                    tg_xxxz_xyzz_0[j] = pb_x * tg_xxz_xyzz_0[j] + wp_x[j] * tg_xxz_xyzz_1[j] + fl1_fx * tg_xz_xyzz_0[j] - fl1_fx * fl1_fza * tg_xz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xxz_yzz_1[j];

                    tg_xxxz_xzzz_0[j] = pb_x * tg_xxz_xzzz_0[j] + wp_x[j] * tg_xxz_xzzz_1[j] + fl1_fx * tg_xz_xzzz_0[j] - fl1_fx * fl1_fza * tg_xz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xxz_zzz_1[j];

                    tg_xxxz_yyyy_0[j] = pb_x * tg_xxz_yyyy_0[j] + wp_x[j] * tg_xxz_yyyy_1[j] + fl1_fx * tg_xz_yyyy_0[j] - fl1_fx * fl1_fza * tg_xz_yyyy_1[j];

                    tg_xxxz_yyyz_0[j] = pb_x * tg_xxz_yyyz_0[j] + wp_x[j] * tg_xxz_yyyz_1[j] + fl1_fx * tg_xz_yyyz_0[j] - fl1_fx * fl1_fza * tg_xz_yyyz_1[j];

                    tg_xxxz_yyzz_0[j] = pb_x * tg_xxz_yyzz_0[j] + wp_x[j] * tg_xxz_yyzz_1[j] + fl1_fx * tg_xz_yyzz_0[j] - fl1_fx * fl1_fza * tg_xz_yyzz_1[j];

                    tg_xxxz_yzzz_0[j] = pb_x * tg_xxz_yzzz_0[j] + wp_x[j] * tg_xxz_yzzz_1[j] + fl1_fx * tg_xz_yzzz_0[j] - fl1_fx * fl1_fza * tg_xz_yzzz_1[j];

                    tg_xxxz_zzzz_0[j] = pb_x * tg_xxz_zzzz_0[j] + wp_x[j] * tg_xxz_zzzz_1[j] + fl1_fx * tg_xz_zzzz_0[j] - fl1_fx * fl1_fza * tg_xz_zzzz_1[j];

                    tg_xxyy_xxxx_0[j] = pb_x * tg_xyy_xxxx_0[j] + wp_x[j] * tg_xyy_xxxx_1[j] + 0.5 * fl1_fx * tg_yy_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyy_xxx_1[j];

                    tg_xxyy_xxxy_0[j] = pb_x * tg_xyy_xxxy_0[j] + wp_x[j] * tg_xyy_xxxy_1[j] + 0.5 * fl1_fx * tg_yy_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyy_xxy_1[j];

                    tg_xxyy_xxxz_0[j] = pb_x * tg_xyy_xxxz_0[j] + wp_x[j] * tg_xyy_xxxz_1[j] + 0.5 * fl1_fx * tg_yy_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyy_xxz_1[j];

                    tg_xxyy_xxyy_0[j] = pb_x * tg_xyy_xxyy_0[j] + wp_x[j] * tg_xyy_xxyy_1[j] + 0.5 * fl1_fx * tg_yy_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyy_1[j] + fl1_fxn * tg_xyy_xyy_1[j];

                    tg_xxyy_xxyz_0[j] = pb_x * tg_xyy_xxyz_0[j] + wp_x[j] * tg_xyy_xxyz_1[j] + 0.5 * fl1_fx * tg_yy_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxyz_1[j] + fl1_fxn * tg_xyy_xyz_1[j];

                    tg_xxyy_xxzz_0[j] = pb_x * tg_xyy_xxzz_0[j] + wp_x[j] * tg_xyy_xxzz_1[j] + 0.5 * fl1_fx * tg_yy_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xxzz_1[j] + fl1_fxn * tg_xyy_xzz_1[j];

                    tg_xxyy_xyyy_0[j] = pb_x * tg_xyy_xyyy_0[j] + wp_x[j] * tg_xyy_xyyy_1[j] + 0.5 * fl1_fx * tg_yy_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyy_yyy_1[j];

                    tg_xxyy_xyyz_0[j] = pb_x * tg_xyy_xyyz_0[j] + wp_x[j] * tg_xyy_xyyz_1[j] + 0.5 * fl1_fx * tg_yy_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyy_yyz_1[j];

                    tg_xxyy_xyzz_0[j] = pb_x * tg_xyy_xyzz_0[j] + wp_x[j] * tg_xyy_xyzz_1[j] + 0.5 * fl1_fx * tg_yy_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyy_yzz_1[j];

                    tg_xxyy_xzzz_0[j] = pb_x * tg_xyy_xzzz_0[j] + wp_x[j] * tg_xyy_xzzz_1[j] + 0.5 * fl1_fx * tg_yy_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyy_zzz_1[j];

                    tg_xxyy_yyyy_0[j] = pb_x * tg_xyy_yyyy_0[j] + wp_x[j] * tg_xyy_yyyy_1[j] + 0.5 * fl1_fx * tg_yy_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyy_1[j];

                    tg_xxyy_yyyz_0[j] = pb_x * tg_xyy_yyyz_0[j] + wp_x[j] * tg_xyy_yyyz_1[j] + 0.5 * fl1_fx * tg_yy_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyyz_1[j];

                    tg_xxyy_yyzz_0[j] = pb_x * tg_xyy_yyzz_0[j] + wp_x[j] * tg_xyy_yyzz_1[j] + 0.5 * fl1_fx * tg_yy_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yyzz_1[j];

                    tg_xxyy_yzzz_0[j] = pb_x * tg_xyy_yzzz_0[j] + wp_x[j] * tg_xyy_yzzz_1[j] + 0.5 * fl1_fx * tg_yy_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_yzzz_1[j];

                    tg_xxyy_zzzz_0[j] = pb_x * tg_xyy_zzzz_0[j] + wp_x[j] * tg_xyy_zzzz_1[j] + 0.5 * fl1_fx * tg_yy_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_zzzz_1[j];

                    tg_xxyz_xxxx_0[j] = pb_x * tg_xyz_xxxx_0[j] + wp_x[j] * tg_xyz_xxxx_1[j] + 0.5 * fl1_fx * tg_yz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xyz_xxx_1[j];

                    tg_xxyz_xxxy_0[j] = pb_x * tg_xyz_xxxy_0[j] + wp_x[j] * tg_xyz_xxxy_1[j] + 0.5 * fl1_fx * tg_yz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xyz_xxy_1[j];

                    tg_xxyz_xxxz_0[j] = pb_x * tg_xyz_xxxz_0[j] + wp_x[j] * tg_xyz_xxxz_1[j] + 0.5 * fl1_fx * tg_yz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xyz_xxz_1[j];

                    tg_xxyz_xxyy_0[j] = pb_x * tg_xyz_xxyy_0[j] + wp_x[j] * tg_xyz_xxyy_1[j] + 0.5 * fl1_fx * tg_yz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyy_1[j] + fl1_fxn * tg_xyz_xyy_1[j];

                    tg_xxyz_xxyz_0[j] = pb_x * tg_xyz_xxyz_0[j] + wp_x[j] * tg_xyz_xxyz_1[j] + 0.5 * fl1_fx * tg_yz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxyz_1[j] + fl1_fxn * tg_xyz_xyz_1[j];

                    tg_xxyz_xxzz_0[j] = pb_x * tg_xyz_xxzz_0[j] + wp_x[j] * tg_xyz_xxzz_1[j] + 0.5 * fl1_fx * tg_yz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xxzz_1[j] + fl1_fxn * tg_xyz_xzz_1[j];

                    tg_xxyz_xyyy_0[j] = pb_x * tg_xyz_xyyy_0[j] + wp_x[j] * tg_xyz_xyyy_1[j] + 0.5 * fl1_fx * tg_yz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xyz_yyy_1[j];

                    tg_xxyz_xyyz_0[j] = pb_x * tg_xyz_xyyz_0[j] + wp_x[j] * tg_xyz_xyyz_1[j] + 0.5 * fl1_fx * tg_yz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xyz_yyz_1[j];

                    tg_xxyz_xyzz_0[j] = pb_x * tg_xyz_xyzz_0[j] + wp_x[j] * tg_xyz_xyzz_1[j] + 0.5 * fl1_fx * tg_yz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xyz_yzz_1[j];

                    tg_xxyz_xzzz_0[j] = pb_x * tg_xyz_xzzz_0[j] + wp_x[j] * tg_xyz_xzzz_1[j] + 0.5 * fl1_fx * tg_yz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xyz_zzz_1[j];

                    tg_xxyz_yyyy_0[j] = pb_x * tg_xyz_yyyy_0[j] + wp_x[j] * tg_xyz_yyyy_1[j] + 0.5 * fl1_fx * tg_yz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyy_1[j];

                    tg_xxyz_yyyz_0[j] = pb_x * tg_xyz_yyyz_0[j] + wp_x[j] * tg_xyz_yyyz_1[j] + 0.5 * fl1_fx * tg_yz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyyz_1[j];

                    tg_xxyz_yyzz_0[j] = pb_x * tg_xyz_yyzz_0[j] + wp_x[j] * tg_xyz_yyzz_1[j] + 0.5 * fl1_fx * tg_yz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yyzz_1[j];

                    tg_xxyz_yzzz_0[j] = pb_x * tg_xyz_yzzz_0[j] + wp_x[j] * tg_xyz_yzzz_1[j] + 0.5 * fl1_fx * tg_yz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_yzzz_1[j];

                    tg_xxyz_zzzz_0[j] = pb_x * tg_xyz_zzzz_0[j] + wp_x[j] * tg_xyz_zzzz_1[j] + 0.5 * fl1_fx * tg_yz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSG_75_150(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (75,150)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_xzz_xxxx_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 75); 

                auto tg_xzz_xxxy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 76); 

                auto tg_xzz_xxxz_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 77); 

                auto tg_xzz_xxyy_0 = primBuffer.data(pidx_g_3_4_m0 + 150 * idx + 78); 

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

                auto tg_xzz_xxxx_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 75); 

                auto tg_xzz_xxxy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 76); 

                auto tg_xzz_xxxz_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 77); 

                auto tg_xzz_xxyy_1 = primBuffer.data(pidx_g_3_4_m1 + 150 * idx + 78); 

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

                auto tg_zz_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 75); 

                auto tg_zz_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 76); 

                auto tg_zz_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 77); 

                auto tg_zz_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 78); 

                auto tg_zz_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 79); 

                auto tg_zz_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 80); 

                auto tg_zz_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 81); 

                auto tg_zz_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 82); 

                auto tg_zz_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 83); 

                auto tg_zz_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 84); 

                auto tg_zz_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 85); 

                auto tg_zz_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 86); 

                auto tg_zz_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 87); 

                auto tg_zz_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 88); 

                auto tg_zz_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 89); 

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

                auto tg_xzz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 50); 

                auto tg_xzz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 51); 

                auto tg_xzz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 52); 

                auto tg_xzz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 53); 

                auto tg_xzz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 54); 

                auto tg_xzz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 55); 

                auto tg_xzz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 56); 

                auto tg_xzz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 57); 

                auto tg_xzz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 58); 

                auto tg_xzz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 59); 

                auto tg_yyy_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 60); 

                auto tg_yyy_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 61); 

                auto tg_yyy_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 62); 

                auto tg_yyy_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 63); 

                auto tg_yyy_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 64); 

                auto tg_yyy_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 65); 

                auto tg_yyy_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 66); 

                auto tg_yyy_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 67); 

                auto tg_yyy_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 68); 

                auto tg_yyy_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 69); 

                auto tg_yyz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 70); 

                auto tg_yyz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 71); 

                auto tg_yyz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 72); 

                auto tg_yyz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 73); 

                auto tg_yyz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 74); 

                auto tg_yyz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 75); 

                auto tg_yyz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 76); 

                auto tg_yyz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 77); 

                auto tg_yyz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 78); 

                auto tg_yyz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 79); 

                auto tg_yzz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 80); 

                auto tg_yzz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 81); 

                auto tg_yzz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 82); 

                auto tg_yzz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 83); 

                auto tg_yzz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 84); 

                auto tg_yzz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 85); 

                auto tg_yzz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 86); 

                auto tg_yzz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 87); 

                auto tg_yzz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 88); 

                auto tg_yzz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 89); 

                auto tg_zzz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 90); 

                auto tg_zzz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 91); 

                auto tg_zzz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 92); 

                auto tg_zzz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 93); 

                auto tg_zzz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 94); 

                auto tg_zzz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 95); 

                auto tg_zzz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 96); 

                auto tg_zzz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 97); 

                auto tg_zzz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 98); 

                auto tg_zzz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 99); 

                // set up pointers to integrals

                auto tg_xxzz_xxxx_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 75); 

                auto tg_xxzz_xxxy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 76); 

                auto tg_xxzz_xxxz_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 77); 

                auto tg_xxzz_xxyy_0 = primBuffer.data(pidx_g_4_4_m0 + 225 * idx + 78); 

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

                // Batch of Integrals (75,150)

                #pragma omp simd aligned(fxn, fza, tg_xxzz_xxxx_0, tg_xxzz_xxxy_0, tg_xxzz_xxxz_0, \
                                         tg_xxzz_xxyy_0, tg_xxzz_xxyz_0, tg_xxzz_xxzz_0, tg_xxzz_xyyy_0, tg_xxzz_xyyz_0, \
                                         tg_xxzz_xyzz_0, tg_xxzz_xzzz_0, tg_xxzz_yyyy_0, tg_xxzz_yyyz_0, tg_xxzz_yyzz_0, \
                                         tg_xxzz_yzzz_0, tg_xxzz_zzzz_0, tg_xyyy_xxxx_0, tg_xyyy_xxxy_0, tg_xyyy_xxxz_0, \
                                         tg_xyyy_xxyy_0, tg_xyyy_xxyz_0, tg_xyyy_xxzz_0, tg_xyyy_xyyy_0, tg_xyyy_xyyz_0, \
                                         tg_xyyy_xyzz_0, tg_xyyy_xzzz_0, tg_xyyy_yyyy_0, tg_xyyy_yyyz_0, tg_xyyy_yyzz_0, \
                                         tg_xyyy_yzzz_0, tg_xyyy_zzzz_0, tg_xyyz_xxxx_0, tg_xyyz_xxxy_0, tg_xyyz_xxxz_0, \
                                         tg_xyyz_xxyy_0, tg_xyyz_xxyz_0, tg_xyyz_xxzz_0, tg_xyyz_xyyy_0, tg_xyyz_xyyz_0, \
                                         tg_xyyz_xyzz_0, tg_xyyz_xzzz_0, tg_xyyz_yyyy_0, tg_xyyz_yyyz_0, tg_xyyz_yyzz_0, \
                                         tg_xyyz_yzzz_0, tg_xyyz_zzzz_0, tg_xyzz_xxxx_0, tg_xyzz_xxxy_0, tg_xyzz_xxxz_0, \
                                         tg_xyzz_xxyy_0, tg_xyzz_xxyz_0, tg_xyzz_xxzz_0, tg_xyzz_xyyy_0, tg_xyzz_xyyz_0, \
                                         tg_xyzz_xyzz_0, tg_xyzz_xzzz_0, tg_xyzz_yyyy_0, tg_xyzz_yyyz_0, tg_xyzz_yyzz_0, \
                                         tg_xyzz_yzzz_0, tg_xyzz_zzzz_0, tg_xzz_xxx_1, tg_xzz_xxxx_0, tg_xzz_xxxx_1, \
                                         tg_xzz_xxxy_0, tg_xzz_xxxy_1, tg_xzz_xxxz_0, tg_xzz_xxxz_1, tg_xzz_xxy_1, \
                                         tg_xzz_xxyy_0, tg_xzz_xxyy_1, tg_xzz_xxyz_0, tg_xzz_xxyz_1, tg_xzz_xxz_1, \
                                         tg_xzz_xxzz_0, tg_xzz_xxzz_1, tg_xzz_xyy_1, tg_xzz_xyyy_0, tg_xzz_xyyy_1, \
                                         tg_xzz_xyyz_0, tg_xzz_xyyz_1, tg_xzz_xyz_1, tg_xzz_xyzz_0, tg_xzz_xyzz_1, \
                                         tg_xzz_xzz_1, tg_xzz_xzzz_0, tg_xzz_xzzz_1, tg_xzz_yyy_1, tg_xzz_yyyy_0, \
                                         tg_xzz_yyyy_1, tg_xzz_yyyz_0, tg_xzz_yyyz_1, tg_xzz_yyz_1, tg_xzz_yyzz_0, \
                                         tg_xzz_yyzz_1, tg_xzz_yzz_1, tg_xzz_yzzz_0, tg_xzz_yzzz_1, tg_xzz_zzz_1, \
                                         tg_xzz_zzzz_0, tg_xzz_zzzz_1, tg_xzzz_xxxx_0, tg_xzzz_xxxy_0, tg_xzzz_xxxz_0, \
                                         tg_xzzz_xxyy_0, tg_xzzz_xxyz_0, tg_xzzz_xxzz_0, tg_xzzz_xyyy_0, tg_xzzz_xyyz_0, \
                                         tg_xzzz_xyzz_0, tg_xzzz_xzzz_0, tg_xzzz_yyyy_0, tg_xzzz_yyyz_0, tg_xzzz_yyzz_0, \
                                         tg_xzzz_yzzz_0, tg_xzzz_zzzz_0, tg_yyy_xxx_1, tg_yyy_xxxx_0, tg_yyy_xxxx_1, \
                                         tg_yyy_xxxy_0, tg_yyy_xxxy_1, tg_yyy_xxxz_0, tg_yyy_xxxz_1, tg_yyy_xxy_1, \
                                         tg_yyy_xxyy_0, tg_yyy_xxyy_1, tg_yyy_xxyz_0, tg_yyy_xxyz_1, tg_yyy_xxz_1, \
                                         tg_yyy_xxzz_0, tg_yyy_xxzz_1, tg_yyy_xyy_1, tg_yyy_xyyy_0, tg_yyy_xyyy_1, \
                                         tg_yyy_xyyz_0, tg_yyy_xyyz_1, tg_yyy_xyz_1, tg_yyy_xyzz_0, tg_yyy_xyzz_1, \
                                         tg_yyy_xzz_1, tg_yyy_xzzz_0, tg_yyy_xzzz_1, tg_yyy_yyy_1, tg_yyy_yyyy_0, \
                                         tg_yyy_yyyy_1, tg_yyy_yyyz_0, tg_yyy_yyyz_1, tg_yyy_yyz_1, tg_yyy_yyzz_0, \
                                         tg_yyy_yyzz_1, tg_yyy_yzz_1, tg_yyy_yzzz_0, tg_yyy_yzzz_1, tg_yyy_zzz_1, \
                                         tg_yyy_zzzz_0, tg_yyy_zzzz_1, tg_yyz_xxx_1, tg_yyz_xxxx_0, tg_yyz_xxxx_1, \
                                         tg_yyz_xxxy_0, tg_yyz_xxxy_1, tg_yyz_xxxz_0, tg_yyz_xxxz_1, tg_yyz_xxy_1, \
                                         tg_yyz_xxyy_0, tg_yyz_xxyy_1, tg_yyz_xxyz_0, tg_yyz_xxyz_1, tg_yyz_xxz_1, \
                                         tg_yyz_xxzz_0, tg_yyz_xxzz_1, tg_yyz_xyy_1, tg_yyz_xyyy_0, tg_yyz_xyyy_1, \
                                         tg_yyz_xyyz_0, tg_yyz_xyyz_1, tg_yyz_xyz_1, tg_yyz_xyzz_0, tg_yyz_xyzz_1, \
                                         tg_yyz_xzz_1, tg_yyz_xzzz_0, tg_yyz_xzzz_1, tg_yyz_yyy_1, tg_yyz_yyyy_0, \
                                         tg_yyz_yyyy_1, tg_yyz_yyyz_0, tg_yyz_yyyz_1, tg_yyz_yyz_1, tg_yyz_yyzz_0, \
                                         tg_yyz_yyzz_1, tg_yyz_yzz_1, tg_yyz_yzzz_0, tg_yyz_yzzz_1, tg_yyz_zzz_1, \
                                         tg_yyz_zzzz_0, tg_yyz_zzzz_1, tg_yzz_xxx_1, tg_yzz_xxxx_0, tg_yzz_xxxx_1, \
                                         tg_yzz_xxxy_0, tg_yzz_xxxy_1, tg_yzz_xxxz_0, tg_yzz_xxxz_1, tg_yzz_xxy_1, \
                                         tg_yzz_xxyy_0, tg_yzz_xxyy_1, tg_yzz_xxyz_0, tg_yzz_xxyz_1, tg_yzz_xxz_1, \
                                         tg_yzz_xxzz_0, tg_yzz_xxzz_1, tg_yzz_xyy_1, tg_yzz_xyyy_0, tg_yzz_xyyy_1, \
                                         tg_yzz_xyyz_0, tg_yzz_xyyz_1, tg_yzz_xyz_1, tg_yzz_xyzz_0, tg_yzz_xyzz_1, \
                                         tg_yzz_xzz_1, tg_yzz_xzzz_0, tg_yzz_xzzz_1, tg_yzz_yyy_1, tg_yzz_yyyy_0, \
                                         tg_yzz_yyyy_1, tg_yzz_yyyz_0, tg_yzz_yyyz_1, tg_yzz_yyz_1, tg_yzz_yyzz_0, \
                                         tg_yzz_yyzz_1, tg_yzz_yzz_1, tg_yzz_yzzz_0, tg_yzz_yzzz_1, tg_yzz_zzz_1, \
                                         tg_yzz_zzzz_0, tg_yzz_zzzz_1, tg_zz_xxxx_0, tg_zz_xxxx_1, tg_zz_xxxy_0, \
                                         tg_zz_xxxy_1, tg_zz_xxxz_0, tg_zz_xxxz_1, tg_zz_xxyy_0, tg_zz_xxyy_1, tg_zz_xxyz_0, \
                                         tg_zz_xxyz_1, tg_zz_xxzz_0, tg_zz_xxzz_1, tg_zz_xyyy_0, tg_zz_xyyy_1, tg_zz_xyyz_0, \
                                         tg_zz_xyyz_1, tg_zz_xyzz_0, tg_zz_xyzz_1, tg_zz_xzzz_0, tg_zz_xzzz_1, tg_zz_yyyy_0, \
                                         tg_zz_yyyy_1, tg_zz_yyyz_0, tg_zz_yyyz_1, tg_zz_yyzz_0, tg_zz_yyzz_1, tg_zz_yzzz_0, \
                                         tg_zz_yzzz_1, tg_zz_zzzz_0, tg_zz_zzzz_1, tg_zzz_xxx_1, tg_zzz_xxxx_0, \
                                         tg_zzz_xxxx_1, tg_zzz_xxxy_0, tg_zzz_xxxy_1, tg_zzz_xxxz_0, tg_zzz_xxxz_1, \
                                         tg_zzz_xxy_1, tg_zzz_xxyy_0, tg_zzz_xxyy_1, tg_zzz_xxyz_0, tg_zzz_xxyz_1, \
                                         tg_zzz_xxz_1, tg_zzz_xxzz_0, tg_zzz_xxzz_1, tg_zzz_xyy_1, tg_zzz_xyyy_0, \
                                         tg_zzz_xyyy_1, tg_zzz_xyyz_0, tg_zzz_xyyz_1, tg_zzz_xyz_1, tg_zzz_xyzz_0, \
                                         tg_zzz_xyzz_1, tg_zzz_xzz_1, tg_zzz_xzzz_0, tg_zzz_xzzz_1, tg_zzz_yyy_1, \
                                         tg_zzz_yyyy_0, tg_zzz_yyyy_1, tg_zzz_yyyz_0, tg_zzz_yyyz_1, tg_zzz_yyz_1, \
                                         tg_zzz_yyzz_0, tg_zzz_yyzz_1, tg_zzz_yzz_1, tg_zzz_yzzz_0, tg_zzz_yzzz_1, \
                                         tg_zzz_zzz_1, tg_zzz_zzzz_0, tg_zzz_zzzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxzz_xxxx_0[j] = pb_x * tg_xzz_xxxx_0[j] + wp_x[j] * tg_xzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxx_1[j] + 2.0 * fl1_fxn * tg_xzz_xxx_1[j];

                    tg_xxzz_xxxy_0[j] = pb_x * tg_xzz_xxxy_0[j] + wp_x[j] * tg_xzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxy_1[j] + 1.5 * fl1_fxn * tg_xzz_xxy_1[j];

                    tg_xxzz_xxxz_0[j] = pb_x * tg_xzz_xxxz_0[j] + wp_x[j] * tg_xzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxz_1[j] + 1.5 * fl1_fxn * tg_xzz_xxz_1[j];

                    tg_xxzz_xxyy_0[j] = pb_x * tg_xzz_xxyy_0[j] + wp_x[j] * tg_xzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyy_1[j] + fl1_fxn * tg_xzz_xyy_1[j];

                    tg_xxzz_xxyz_0[j] = pb_x * tg_xzz_xxyz_0[j] + wp_x[j] * tg_xzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyz_1[j] + fl1_fxn * tg_xzz_xyz_1[j];

                    tg_xxzz_xxzz_0[j] = pb_x * tg_xzz_xxzz_0[j] + wp_x[j] * tg_xzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxzz_1[j] + fl1_fxn * tg_xzz_xzz_1[j];

                    tg_xxzz_xyyy_0[j] = pb_x * tg_xzz_xyyy_0[j] + wp_x[j] * tg_xzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyy_1[j] + 0.5 * fl1_fxn * tg_xzz_yyy_1[j];

                    tg_xxzz_xyyz_0[j] = pb_x * tg_xzz_xyyz_0[j] + wp_x[j] * tg_xzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyz_1[j] + 0.5 * fl1_fxn * tg_xzz_yyz_1[j];

                    tg_xxzz_xyzz_0[j] = pb_x * tg_xzz_xyzz_0[j] + wp_x[j] * tg_xzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyzz_1[j] + 0.5 * fl1_fxn * tg_xzz_yzz_1[j];

                    tg_xxzz_xzzz_0[j] = pb_x * tg_xzz_xzzz_0[j] + wp_x[j] * tg_xzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xzzz_1[j] + 0.5 * fl1_fxn * tg_xzz_zzz_1[j];

                    tg_xxzz_yyyy_0[j] = pb_x * tg_xzz_yyyy_0[j] + wp_x[j] * tg_xzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyy_1[j];

                    tg_xxzz_yyyz_0[j] = pb_x * tg_xzz_yyyz_0[j] + wp_x[j] * tg_xzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyz_1[j];

                    tg_xxzz_yyzz_0[j] = pb_x * tg_xzz_yyzz_0[j] + wp_x[j] * tg_xzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyzz_1[j];

                    tg_xxzz_yzzz_0[j] = pb_x * tg_xzz_yzzz_0[j] + wp_x[j] * tg_xzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yzzz_1[j];

                    tg_xxzz_zzzz_0[j] = pb_x * tg_xzz_zzzz_0[j] + wp_x[j] * tg_xzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_zzzz_1[j];

                    tg_xyyy_xxxx_0[j] = pb_x * tg_yyy_xxxx_0[j] + wp_x[j] * tg_yyy_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyy_xxx_1[j];

                    tg_xyyy_xxxy_0[j] = pb_x * tg_yyy_xxxy_0[j] + wp_x[j] * tg_yyy_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyy_xxy_1[j];

                    tg_xyyy_xxxz_0[j] = pb_x * tg_yyy_xxxz_0[j] + wp_x[j] * tg_yyy_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyy_xxz_1[j];

                    tg_xyyy_xxyy_0[j] = pb_x * tg_yyy_xxyy_0[j] + wp_x[j] * tg_yyy_xxyy_1[j] + fl1_fxn * tg_yyy_xyy_1[j];

                    tg_xyyy_xxyz_0[j] = pb_x * tg_yyy_xxyz_0[j] + wp_x[j] * tg_yyy_xxyz_1[j] + fl1_fxn * tg_yyy_xyz_1[j];

                    tg_xyyy_xxzz_0[j] = pb_x * tg_yyy_xxzz_0[j] + wp_x[j] * tg_yyy_xxzz_1[j] + fl1_fxn * tg_yyy_xzz_1[j];

                    tg_xyyy_xyyy_0[j] = pb_x * tg_yyy_xyyy_0[j] + wp_x[j] * tg_yyy_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyy_yyy_1[j];

                    tg_xyyy_xyyz_0[j] = pb_x * tg_yyy_xyyz_0[j] + wp_x[j] * tg_yyy_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyy_yyz_1[j];

                    tg_xyyy_xyzz_0[j] = pb_x * tg_yyy_xyzz_0[j] + wp_x[j] * tg_yyy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyy_yzz_1[j];

                    tg_xyyy_xzzz_0[j] = pb_x * tg_yyy_xzzz_0[j] + wp_x[j] * tg_yyy_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zzz_1[j];

                    tg_xyyy_yyyy_0[j] = pb_x * tg_yyy_yyyy_0[j] + wp_x[j] * tg_yyy_yyyy_1[j];

                    tg_xyyy_yyyz_0[j] = pb_x * tg_yyy_yyyz_0[j] + wp_x[j] * tg_yyy_yyyz_1[j];

                    tg_xyyy_yyzz_0[j] = pb_x * tg_yyy_yyzz_0[j] + wp_x[j] * tg_yyy_yyzz_1[j];

                    tg_xyyy_yzzz_0[j] = pb_x * tg_yyy_yzzz_0[j] + wp_x[j] * tg_yyy_yzzz_1[j];

                    tg_xyyy_zzzz_0[j] = pb_x * tg_yyy_zzzz_0[j] + wp_x[j] * tg_yyy_zzzz_1[j];

                    tg_xyyz_xxxx_0[j] = pb_x * tg_yyz_xxxx_0[j] + wp_x[j] * tg_yyz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yyz_xxx_1[j];

                    tg_xyyz_xxxy_0[j] = pb_x * tg_yyz_xxxy_0[j] + wp_x[j] * tg_yyz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yyz_xxy_1[j];

                    tg_xyyz_xxxz_0[j] = pb_x * tg_yyz_xxxz_0[j] + wp_x[j] * tg_yyz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yyz_xxz_1[j];

                    tg_xyyz_xxyy_0[j] = pb_x * tg_yyz_xxyy_0[j] + wp_x[j] * tg_yyz_xxyy_1[j] + fl1_fxn * tg_yyz_xyy_1[j];

                    tg_xyyz_xxyz_0[j] = pb_x * tg_yyz_xxyz_0[j] + wp_x[j] * tg_yyz_xxyz_1[j] + fl1_fxn * tg_yyz_xyz_1[j];

                    tg_xyyz_xxzz_0[j] = pb_x * tg_yyz_xxzz_0[j] + wp_x[j] * tg_yyz_xxzz_1[j] + fl1_fxn * tg_yyz_xzz_1[j];

                    tg_xyyz_xyyy_0[j] = pb_x * tg_yyz_xyyy_0[j] + wp_x[j] * tg_yyz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yyz_yyy_1[j];

                    tg_xyyz_xyyz_0[j] = pb_x * tg_yyz_xyyz_0[j] + wp_x[j] * tg_yyz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yyz_yyz_1[j];

                    tg_xyyz_xyzz_0[j] = pb_x * tg_yyz_xyzz_0[j] + wp_x[j] * tg_yyz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyz_yzz_1[j];

                    tg_xyyz_xzzz_0[j] = pb_x * tg_yyz_xzzz_0[j] + wp_x[j] * tg_yyz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zzz_1[j];

                    tg_xyyz_yyyy_0[j] = pb_x * tg_yyz_yyyy_0[j] + wp_x[j] * tg_yyz_yyyy_1[j];

                    tg_xyyz_yyyz_0[j] = pb_x * tg_yyz_yyyz_0[j] + wp_x[j] * tg_yyz_yyyz_1[j];

                    tg_xyyz_yyzz_0[j] = pb_x * tg_yyz_yyzz_0[j] + wp_x[j] * tg_yyz_yyzz_1[j];

                    tg_xyyz_yzzz_0[j] = pb_x * tg_yyz_yzzz_0[j] + wp_x[j] * tg_yyz_yzzz_1[j];

                    tg_xyyz_zzzz_0[j] = pb_x * tg_yyz_zzzz_0[j] + wp_x[j] * tg_yyz_zzzz_1[j];

                    tg_xyzz_xxxx_0[j] = pb_x * tg_yzz_xxxx_0[j] + wp_x[j] * tg_yzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_yzz_xxx_1[j];

                    tg_xyzz_xxxy_0[j] = pb_x * tg_yzz_xxxy_0[j] + wp_x[j] * tg_yzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_yzz_xxy_1[j];

                    tg_xyzz_xxxz_0[j] = pb_x * tg_yzz_xxxz_0[j] + wp_x[j] * tg_yzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_yzz_xxz_1[j];

                    tg_xyzz_xxyy_0[j] = pb_x * tg_yzz_xxyy_0[j] + wp_x[j] * tg_yzz_xxyy_1[j] + fl1_fxn * tg_yzz_xyy_1[j];

                    tg_xyzz_xxyz_0[j] = pb_x * tg_yzz_xxyz_0[j] + wp_x[j] * tg_yzz_xxyz_1[j] + fl1_fxn * tg_yzz_xyz_1[j];

                    tg_xyzz_xxzz_0[j] = pb_x * tg_yzz_xxzz_0[j] + wp_x[j] * tg_yzz_xxzz_1[j] + fl1_fxn * tg_yzz_xzz_1[j];

                    tg_xyzz_xyyy_0[j] = pb_x * tg_yzz_xyyy_0[j] + wp_x[j] * tg_yzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_yzz_yyy_1[j];

                    tg_xyzz_xyyz_0[j] = pb_x * tg_yzz_xyyz_0[j] + wp_x[j] * tg_yzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_yzz_yyz_1[j];

                    tg_xyzz_xyzz_0[j] = pb_x * tg_yzz_xyzz_0[j] + wp_x[j] * tg_yzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzz_yzz_1[j];

                    tg_xyzz_xzzz_0[j] = pb_x * tg_yzz_xzzz_0[j] + wp_x[j] * tg_yzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zzz_1[j];

                    tg_xyzz_yyyy_0[j] = pb_x * tg_yzz_yyyy_0[j] + wp_x[j] * tg_yzz_yyyy_1[j];

                    tg_xyzz_yyyz_0[j] = pb_x * tg_yzz_yyyz_0[j] + wp_x[j] * tg_yzz_yyyz_1[j];

                    tg_xyzz_yyzz_0[j] = pb_x * tg_yzz_yyzz_0[j] + wp_x[j] * tg_yzz_yyzz_1[j];

                    tg_xyzz_yzzz_0[j] = pb_x * tg_yzz_yzzz_0[j] + wp_x[j] * tg_yzz_yzzz_1[j];

                    tg_xyzz_zzzz_0[j] = pb_x * tg_yzz_zzzz_0[j] + wp_x[j] * tg_yzz_zzzz_1[j];

                    tg_xzzz_xxxx_0[j] = pb_x * tg_zzz_xxxx_0[j] + wp_x[j] * tg_zzz_xxxx_1[j] + 2.0 * fl1_fxn * tg_zzz_xxx_1[j];

                    tg_xzzz_xxxy_0[j] = pb_x * tg_zzz_xxxy_0[j] + wp_x[j] * tg_zzz_xxxy_1[j] + 1.5 * fl1_fxn * tg_zzz_xxy_1[j];

                    tg_xzzz_xxxz_0[j] = pb_x * tg_zzz_xxxz_0[j] + wp_x[j] * tg_zzz_xxxz_1[j] + 1.5 * fl1_fxn * tg_zzz_xxz_1[j];

                    tg_xzzz_xxyy_0[j] = pb_x * tg_zzz_xxyy_0[j] + wp_x[j] * tg_zzz_xxyy_1[j] + fl1_fxn * tg_zzz_xyy_1[j];

                    tg_xzzz_xxyz_0[j] = pb_x * tg_zzz_xxyz_0[j] + wp_x[j] * tg_zzz_xxyz_1[j] + fl1_fxn * tg_zzz_xyz_1[j];

                    tg_xzzz_xxzz_0[j] = pb_x * tg_zzz_xxzz_0[j] + wp_x[j] * tg_zzz_xxzz_1[j] + fl1_fxn * tg_zzz_xzz_1[j];

                    tg_xzzz_xyyy_0[j] = pb_x * tg_zzz_xyyy_0[j] + wp_x[j] * tg_zzz_xyyy_1[j] + 0.5 * fl1_fxn * tg_zzz_yyy_1[j];

                    tg_xzzz_xyyz_0[j] = pb_x * tg_zzz_xyyz_0[j] + wp_x[j] * tg_zzz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyz_1[j];

                    tg_xzzz_xyzz_0[j] = pb_x * tg_zzz_xyzz_0[j] + wp_x[j] * tg_zzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_yzz_1[j];

                    tg_xzzz_xzzz_0[j] = pb_x * tg_zzz_xzzz_0[j] + wp_x[j] * tg_zzz_xzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzz_1[j];

                    tg_xzzz_yyyy_0[j] = pb_x * tg_zzz_yyyy_0[j] + wp_x[j] * tg_zzz_yyyy_1[j];

                    tg_xzzz_yyyz_0[j] = pb_x * tg_zzz_yyyz_0[j] + wp_x[j] * tg_zzz_yyyz_1[j];

                    tg_xzzz_yyzz_0[j] = pb_x * tg_zzz_yyzz_0[j] + wp_x[j] * tg_zzz_yyzz_1[j];

                    tg_xzzz_yzzz_0[j] = pb_x * tg_zzz_yzzz_0[j] + wp_x[j] * tg_zzz_yzzz_1[j];

                    tg_xzzz_zzzz_0[j] = pb_x * tg_zzz_zzzz_0[j] + wp_x[j] * tg_zzz_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSG_150_225(      CMemBlock2D<double>& primBuffer,
                                         const CRecursionMap&       recursionMap,
                                         const CMemBlock2D<double>& osFactors,
                                         const CMemBlock2D<double>& wpDistances,
                                         const CGtoPairsBlock&      braGtoPairsBlock,
                                         const CGtoPairsBlock&      ketGtoPairsBlock,
                                         const int32_t              nKetPrimPairs,
                                         const int32_t              iContrPair)
    {
        // Batch of Integrals (150,225)

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
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yy_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 45); 

                auto tg_yy_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 46); 

                auto tg_yy_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 47); 

                auto tg_yy_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 48); 

                auto tg_yy_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 49); 

                auto tg_yy_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 50); 

                auto tg_yy_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 51); 

                auto tg_yy_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 52); 

                auto tg_yy_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 53); 

                auto tg_yy_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 54); 

                auto tg_yy_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 55); 

                auto tg_yy_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 56); 

                auto tg_yy_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 57); 

                auto tg_yy_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 58); 

                auto tg_yy_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 59); 

                auto tg_yz_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 60); 

                auto tg_yz_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 61); 

                auto tg_yz_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 62); 

                auto tg_yz_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 63); 

                auto tg_yz_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 64); 

                auto tg_yz_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 65); 

                auto tg_yz_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 66); 

                auto tg_yz_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 67); 

                auto tg_yz_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 68); 

                auto tg_yz_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 69); 

                auto tg_yz_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 70); 

                auto tg_yz_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 71); 

                auto tg_yz_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 72); 

                auto tg_yz_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 73); 

                auto tg_yz_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 74); 

                auto tg_zz_xxxx_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 75); 

                auto tg_zz_xxxy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 76); 

                auto tg_zz_xxxz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 77); 

                auto tg_zz_xxyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 78); 

                auto tg_zz_xxyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 79); 

                auto tg_zz_xxzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 80); 

                auto tg_zz_xyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 81); 

                auto tg_zz_xyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 82); 

                auto tg_zz_xyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 83); 

                auto tg_zz_xzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 84); 

                auto tg_zz_yyyy_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 85); 

                auto tg_zz_yyyz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 86); 

                auto tg_zz_yyzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 87); 

                auto tg_zz_yzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 88); 

                auto tg_zz_zzzz_0 = primBuffer.data(pidx_g_2_4_m0 + 90 * idx + 89); 

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

                auto tg_yyy_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 60); 

                auto tg_yyy_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 61); 

                auto tg_yyy_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 62); 

                auto tg_yyy_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 63); 

                auto tg_yyy_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 64); 

                auto tg_yyy_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 65); 

                auto tg_yyy_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 66); 

                auto tg_yyy_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 67); 

                auto tg_yyy_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 68); 

                auto tg_yyy_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 69); 

                auto tg_yyz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 70); 

                auto tg_yyz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 71); 

                auto tg_yyz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 72); 

                auto tg_yyz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 73); 

                auto tg_yyz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 74); 

                auto tg_yyz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 75); 

                auto tg_yyz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 76); 

                auto tg_yyz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 77); 

                auto tg_yyz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 78); 

                auto tg_yyz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 79); 

                auto tg_yzz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 80); 

                auto tg_yzz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 81); 

                auto tg_yzz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 82); 

                auto tg_yzz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 83); 

                auto tg_yzz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 84); 

                auto tg_yzz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 85); 

                auto tg_yzz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 86); 

                auto tg_yzz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 87); 

                auto tg_yzz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 88); 

                auto tg_yzz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 89); 

                auto tg_zzz_xxx_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 90); 

                auto tg_zzz_xxy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 91); 

                auto tg_zzz_xxz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 92); 

                auto tg_zzz_xyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 93); 

                auto tg_zzz_xyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 94); 

                auto tg_zzz_xzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 95); 

                auto tg_zzz_yyy_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 96); 

                auto tg_zzz_yyz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 97); 

                auto tg_zzz_yzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 98); 

                auto tg_zzz_zzz_1 = primBuffer.data(pidx_g_3_3_m1 + 100 * idx + 99); 

                // set up pointers to integrals

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

                // Batch of Integrals (150,225)

                #pragma omp simd aligned(fxn, fza, tg_yy_xxxx_0, tg_yy_xxxx_1, tg_yy_xxxy_0, tg_yy_xxxy_1, \
                                         tg_yy_xxxz_0, tg_yy_xxxz_1, tg_yy_xxyy_0, tg_yy_xxyy_1, tg_yy_xxyz_0, tg_yy_xxyz_1, \
                                         tg_yy_xxzz_0, tg_yy_xxzz_1, tg_yy_xyyy_0, tg_yy_xyyy_1, tg_yy_xyyz_0, tg_yy_xyyz_1, \
                                         tg_yy_xyzz_0, tg_yy_xyzz_1, tg_yy_xzzz_0, tg_yy_xzzz_1, tg_yy_yyyy_0, tg_yy_yyyy_1, \
                                         tg_yy_yyyz_0, tg_yy_yyyz_1, tg_yy_yyzz_0, tg_yy_yyzz_1, tg_yy_yzzz_0, tg_yy_yzzz_1, \
                                         tg_yy_zzzz_0, tg_yy_zzzz_1, tg_yyy_xxx_1, tg_yyy_xxxx_0, tg_yyy_xxxx_1, \
                                         tg_yyy_xxxy_0, tg_yyy_xxxy_1, tg_yyy_xxxz_0, tg_yyy_xxxz_1, tg_yyy_xxy_1, \
                                         tg_yyy_xxyy_0, tg_yyy_xxyy_1, tg_yyy_xxyz_0, tg_yyy_xxyz_1, tg_yyy_xxz_1, \
                                         tg_yyy_xxzz_0, tg_yyy_xxzz_1, tg_yyy_xyy_1, tg_yyy_xyyy_0, tg_yyy_xyyy_1, \
                                         tg_yyy_xyyz_0, tg_yyy_xyyz_1, tg_yyy_xyz_1, tg_yyy_xyzz_0, tg_yyy_xyzz_1, \
                                         tg_yyy_xzz_1, tg_yyy_xzzz_0, tg_yyy_xzzz_1, tg_yyy_yyy_1, tg_yyy_yyyy_0, \
                                         tg_yyy_yyyy_1, tg_yyy_yyyz_0, tg_yyy_yyyz_1, tg_yyy_yyz_1, tg_yyy_yyzz_0, \
                                         tg_yyy_yyzz_1, tg_yyy_yzz_1, tg_yyy_yzzz_0, tg_yyy_yzzz_1, tg_yyy_zzz_1, \
                                         tg_yyy_zzzz_0, tg_yyy_zzzz_1, tg_yyyy_xxxx_0, tg_yyyy_xxxy_0, tg_yyyy_xxxz_0, \
                                         tg_yyyy_xxyy_0, tg_yyyy_xxyz_0, tg_yyyy_xxzz_0, tg_yyyy_xyyy_0, tg_yyyy_xyyz_0, \
                                         tg_yyyy_xyzz_0, tg_yyyy_xzzz_0, tg_yyyy_yyyy_0, tg_yyyy_yyyz_0, tg_yyyy_yyzz_0, \
                                         tg_yyyy_yzzz_0, tg_yyyy_zzzz_0, tg_yyyz_xxxx_0, tg_yyyz_xxxy_0, tg_yyyz_xxxz_0, \
                                         tg_yyyz_xxyy_0, tg_yyyz_xxyz_0, tg_yyyz_xxzz_0, tg_yyyz_xyyy_0, tg_yyyz_xyyz_0, \
                                         tg_yyyz_xyzz_0, tg_yyyz_xzzz_0, tg_yyyz_yyyy_0, tg_yyyz_yyyz_0, tg_yyyz_yyzz_0, \
                                         tg_yyyz_yzzz_0, tg_yyyz_zzzz_0, tg_yyz_xxx_1, tg_yyz_xxxx_0, tg_yyz_xxxx_1, \
                                         tg_yyz_xxxy_0, tg_yyz_xxxy_1, tg_yyz_xxxz_0, tg_yyz_xxxz_1, tg_yyz_xxy_1, \
                                         tg_yyz_xxyy_0, tg_yyz_xxyy_1, tg_yyz_xxyz_0, tg_yyz_xxyz_1, tg_yyz_xxz_1, \
                                         tg_yyz_xxzz_0, tg_yyz_xxzz_1, tg_yyz_xyy_1, tg_yyz_xyyy_0, tg_yyz_xyyy_1, \
                                         tg_yyz_xyyz_0, tg_yyz_xyyz_1, tg_yyz_xyz_1, tg_yyz_xyzz_0, tg_yyz_xyzz_1, \
                                         tg_yyz_xzz_1, tg_yyz_xzzz_0, tg_yyz_xzzz_1, tg_yyz_yyy_1, tg_yyz_yyyy_0, \
                                         tg_yyz_yyyy_1, tg_yyz_yyyz_0, tg_yyz_yyyz_1, tg_yyz_yyz_1, tg_yyz_yyzz_0, \
                                         tg_yyz_yyzz_1, tg_yyz_yzz_1, tg_yyz_yzzz_0, tg_yyz_yzzz_1, tg_yyz_zzz_1, \
                                         tg_yyz_zzzz_0, tg_yyz_zzzz_1, tg_yyzz_xxxx_0, tg_yyzz_xxxy_0, tg_yyzz_xxxz_0, \
                                         tg_yyzz_xxyy_0, tg_yyzz_xxyz_0, tg_yyzz_xxzz_0, tg_yyzz_xyyy_0, tg_yyzz_xyyz_0, \
                                         tg_yyzz_xyzz_0, tg_yyzz_xzzz_0, tg_yyzz_yyyy_0, tg_yyzz_yyyz_0, tg_yyzz_yyzz_0, \
                                         tg_yyzz_yzzz_0, tg_yyzz_zzzz_0, tg_yz_xxxx_0, tg_yz_xxxx_1, tg_yz_xxxy_0, \
                                         tg_yz_xxxy_1, tg_yz_xxxz_0, tg_yz_xxxz_1, tg_yz_xxyy_0, tg_yz_xxyy_1, tg_yz_xxyz_0, \
                                         tg_yz_xxyz_1, tg_yz_xxzz_0, tg_yz_xxzz_1, tg_yz_xyyy_0, tg_yz_xyyy_1, tg_yz_xyyz_0, \
                                         tg_yz_xyyz_1, tg_yz_xyzz_0, tg_yz_xyzz_1, tg_yz_xzzz_0, tg_yz_xzzz_1, tg_yz_yyyy_0, \
                                         tg_yz_yyyy_1, tg_yz_yyyz_0, tg_yz_yyyz_1, tg_yz_yyzz_0, tg_yz_yyzz_1, tg_yz_yzzz_0, \
                                         tg_yz_yzzz_1, tg_yz_zzzz_0, tg_yz_zzzz_1, tg_yzz_xxx_1, tg_yzz_xxxx_0, \
                                         tg_yzz_xxxx_1, tg_yzz_xxxy_0, tg_yzz_xxxy_1, tg_yzz_xxxz_0, tg_yzz_xxxz_1, \
                                         tg_yzz_xxy_1, tg_yzz_xxyy_0, tg_yzz_xxyy_1, tg_yzz_xxyz_0, tg_yzz_xxyz_1, \
                                         tg_yzz_xxz_1, tg_yzz_xxzz_0, tg_yzz_xxzz_1, tg_yzz_xyy_1, tg_yzz_xyyy_0, \
                                         tg_yzz_xyyy_1, tg_yzz_xyyz_0, tg_yzz_xyyz_1, tg_yzz_xyz_1, tg_yzz_xyzz_0, \
                                         tg_yzz_xyzz_1, tg_yzz_xzz_1, tg_yzz_xzzz_0, tg_yzz_xzzz_1, tg_yzz_yyy_1, \
                                         tg_yzz_yyyy_0, tg_yzz_yyyy_1, tg_yzz_yyyz_0, tg_yzz_yyyz_1, tg_yzz_yyz_1, \
                                         tg_yzz_yyzz_0, tg_yzz_yyzz_1, tg_yzz_yzz_1, tg_yzz_yzzz_0, tg_yzz_yzzz_1, \
                                         tg_yzz_zzz_1, tg_yzz_zzzz_0, tg_yzz_zzzz_1, tg_yzzz_xxxx_0, tg_yzzz_xxxy_0, \
                                         tg_yzzz_xxxz_0, tg_yzzz_xxyy_0, tg_yzzz_xxyz_0, tg_yzzz_xxzz_0, tg_yzzz_xyyy_0, \
                                         tg_yzzz_xyyz_0, tg_yzzz_xyzz_0, tg_yzzz_xzzz_0, tg_yzzz_yyyy_0, tg_yzzz_yyyz_0, \
                                         tg_yzzz_yyzz_0, tg_yzzz_yzzz_0, tg_yzzz_zzzz_0, tg_zz_xxxx_0, tg_zz_xxxx_1, \
                                         tg_zz_xxxy_0, tg_zz_xxxy_1, tg_zz_xxxz_0, tg_zz_xxxz_1, tg_zz_xxyy_0, tg_zz_xxyy_1, \
                                         tg_zz_xxyz_0, tg_zz_xxyz_1, tg_zz_xxzz_0, tg_zz_xxzz_1, tg_zz_xyyy_0, tg_zz_xyyy_1, \
                                         tg_zz_xyyz_0, tg_zz_xyyz_1, tg_zz_xyzz_0, tg_zz_xyzz_1, tg_zz_xzzz_0, tg_zz_xzzz_1, \
                                         tg_zz_yyyy_0, tg_zz_yyyy_1, tg_zz_yyyz_0, tg_zz_yyyz_1, tg_zz_yyzz_0, tg_zz_yyzz_1, \
                                         tg_zz_yzzz_0, tg_zz_yzzz_1, tg_zz_zzzz_0, tg_zz_zzzz_1, tg_zzz_xxx_1, \
                                         tg_zzz_xxxx_0, tg_zzz_xxxx_1, tg_zzz_xxxy_0, tg_zzz_xxxy_1, tg_zzz_xxxz_0, \
                                         tg_zzz_xxxz_1, tg_zzz_xxy_1, tg_zzz_xxyy_0, tg_zzz_xxyy_1, tg_zzz_xxyz_0, \
                                         tg_zzz_xxyz_1, tg_zzz_xxz_1, tg_zzz_xxzz_0, tg_zzz_xxzz_1, tg_zzz_xyy_1, \
                                         tg_zzz_xyyy_0, tg_zzz_xyyy_1, tg_zzz_xyyz_0, tg_zzz_xyyz_1, tg_zzz_xyz_1, \
                                         tg_zzz_xyzz_0, tg_zzz_xyzz_1, tg_zzz_xzz_1, tg_zzz_xzzz_0, tg_zzz_xzzz_1, \
                                         tg_zzz_yyy_1, tg_zzz_yyyy_0, tg_zzz_yyyy_1, tg_zzz_yyyz_0, tg_zzz_yyyz_1, \
                                         tg_zzz_yyz_1, tg_zzz_yyzz_0, tg_zzz_yyzz_1, tg_zzz_yzz_1, tg_zzz_yzzz_0, \
                                         tg_zzz_yzzz_1, tg_zzz_zzz_1, tg_zzz_zzzz_0, tg_zzz_zzzz_1, tg_zzzz_xxxx_0, \
                                         tg_zzzz_xxxy_0, tg_zzzz_xxxz_0, tg_zzzz_xxyy_0, tg_zzzz_xxyz_0, tg_zzzz_xxzz_0, \
                                         tg_zzzz_xyyy_0, tg_zzzz_xyyz_0, tg_zzzz_xyzz_0, tg_zzzz_xzzz_0, tg_zzzz_yyyy_0, \
                                         tg_zzzz_yyyz_0, tg_zzzz_yyzz_0, tg_zzzz_yzzz_0, tg_zzzz_zzzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_yyyy_xxxx_0[j] = pb_y * tg_yyy_xxxx_0[j] + wp_y[j] * tg_yyy_xxxx_1[j] + 1.5 * fl1_fx * tg_yy_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxx_1[j];

                    tg_yyyy_xxxy_0[j] = pb_y * tg_yyy_xxxy_0[j] + wp_y[j] * tg_yyy_xxxy_1[j] + 1.5 * fl1_fx * tg_yy_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyy_xxx_1[j];

                    tg_yyyy_xxxz_0[j] = pb_y * tg_yyy_xxxz_0[j] + wp_y[j] * tg_yyy_xxxz_1[j] + 1.5 * fl1_fx * tg_yy_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxxz_1[j];

                    tg_yyyy_xxyy_0[j] = pb_y * tg_yyy_xxyy_0[j] + wp_y[j] * tg_yyy_xxyy_1[j] + 1.5 * fl1_fx * tg_yy_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyy_1[j] + fl1_fxn * tg_yyy_xxy_1[j];

                    tg_yyyy_xxyz_0[j] = pb_y * tg_yyy_xxyz_0[j] + wp_y[j] * tg_yyy_xxyz_1[j] + 1.5 * fl1_fx * tg_yy_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyy_xxz_1[j];

                    tg_yyyy_xxzz_0[j] = pb_y * tg_yyy_xxzz_0[j] + wp_y[j] * tg_yyy_xxzz_1[j] + 1.5 * fl1_fx * tg_yy_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xxzz_1[j];

                    tg_yyyy_xyyy_0[j] = pb_y * tg_yyy_xyyy_0[j] + wp_y[j] * tg_yyy_xyyy_1[j] + 1.5 * fl1_fx * tg_yy_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyy_xyy_1[j];

                    tg_yyyy_xyyz_0[j] = pb_y * tg_yyy_xyyz_0[j] + wp_y[j] * tg_yyy_xyyz_1[j] + 1.5 * fl1_fx * tg_yy_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyyz_1[j] + fl1_fxn * tg_yyy_xyz_1[j];

                    tg_yyyy_xyzz_0[j] = pb_y * tg_yyy_xyzz_0[j] + wp_y[j] * tg_yyy_xyzz_1[j] + 1.5 * fl1_fx * tg_yy_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyy_xzz_1[j];

                    tg_yyyy_xzzz_0[j] = pb_y * tg_yyy_xzzz_0[j] + wp_y[j] * tg_yyy_xzzz_1[j] + 1.5 * fl1_fx * tg_yy_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_xzzz_1[j];

                    tg_yyyy_yyyy_0[j] = pb_y * tg_yyy_yyyy_0[j] + wp_y[j] * tg_yyy_yyyy_1[j] + 1.5 * fl1_fx * tg_yy_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyy_yyy_1[j];

                    tg_yyyy_yyyz_0[j] = pb_y * tg_yyy_yyyz_0[j] + wp_y[j] * tg_yyy_yyyz_1[j] + 1.5 * fl1_fx * tg_yy_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyy_yyz_1[j];

                    tg_yyyy_yyzz_0[j] = pb_y * tg_yyy_yyzz_0[j] + wp_y[j] * tg_yyy_yyzz_1[j] + 1.5 * fl1_fx * tg_yy_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yyzz_1[j] + fl1_fxn * tg_yyy_yzz_1[j];

                    tg_yyyy_yzzz_0[j] = pb_y * tg_yyy_yzzz_0[j] + wp_y[j] * tg_yyy_yzzz_1[j] + 1.5 * fl1_fx * tg_yy_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zzz_1[j];

                    tg_yyyy_zzzz_0[j] = pb_y * tg_yyy_zzzz_0[j] + wp_y[j] * tg_yyy_zzzz_1[j] + 1.5 * fl1_fx * tg_yy_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_zzzz_1[j];

                    tg_yyyz_xxxx_0[j] = pb_y * tg_yyz_xxxx_0[j] + wp_y[j] * tg_yyz_xxxx_1[j] + fl1_fx * tg_yz_xxxx_0[j] - fl1_fx * fl1_fza * tg_yz_xxxx_1[j];

                    tg_yyyz_xxxy_0[j] = pb_y * tg_yyz_xxxy_0[j] + wp_y[j] * tg_yyz_xxxy_1[j] + fl1_fx * tg_yz_xxxy_0[j] - fl1_fx * fl1_fza * tg_yz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yyz_xxx_1[j];

                    tg_yyyz_xxxz_0[j] = pb_y * tg_yyz_xxxz_0[j] + wp_y[j] * tg_yyz_xxxz_1[j] + fl1_fx * tg_yz_xxxz_0[j] - fl1_fx * fl1_fza * tg_yz_xxxz_1[j];

                    tg_yyyz_xxyy_0[j] = pb_y * tg_yyz_xxyy_0[j] + wp_y[j] * tg_yyz_xxyy_1[j] + fl1_fx * tg_yz_xxyy_0[j] - fl1_fx * fl1_fza * tg_yz_xxyy_1[j] + fl1_fxn * tg_yyz_xxy_1[j];

                    tg_yyyz_xxyz_0[j] = pb_y * tg_yyz_xxyz_0[j] + wp_y[j] * tg_yyz_xxyz_1[j] + fl1_fx * tg_yz_xxyz_0[j] - fl1_fx * fl1_fza * tg_yz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yyz_xxz_1[j];

                    tg_yyyz_xxzz_0[j] = pb_y * tg_yyz_xxzz_0[j] + wp_y[j] * tg_yyz_xxzz_1[j] + fl1_fx * tg_yz_xxzz_0[j] - fl1_fx * fl1_fza * tg_yz_xxzz_1[j];

                    tg_yyyz_xyyy_0[j] = pb_y * tg_yyz_xyyy_0[j] + wp_y[j] * tg_yyz_xyyy_1[j] + fl1_fx * tg_yz_xyyy_0[j] - fl1_fx * fl1_fza * tg_yz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yyz_xyy_1[j];

                    tg_yyyz_xyyz_0[j] = pb_y * tg_yyz_xyyz_0[j] + wp_y[j] * tg_yyz_xyyz_1[j] + fl1_fx * tg_yz_xyyz_0[j] - fl1_fx * fl1_fza * tg_yz_xyyz_1[j] + fl1_fxn * tg_yyz_xyz_1[j];

                    tg_yyyz_xyzz_0[j] = pb_y * tg_yyz_xyzz_0[j] + wp_y[j] * tg_yyz_xyzz_1[j] + fl1_fx * tg_yz_xyzz_0[j] - fl1_fx * fl1_fza * tg_yz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yyz_xzz_1[j];

                    tg_yyyz_xzzz_0[j] = pb_y * tg_yyz_xzzz_0[j] + wp_y[j] * tg_yyz_xzzz_1[j] + fl1_fx * tg_yz_xzzz_0[j] - fl1_fx * fl1_fza * tg_yz_xzzz_1[j];

                    tg_yyyz_yyyy_0[j] = pb_y * tg_yyz_yyyy_0[j] + wp_y[j] * tg_yyz_yyyy_1[j] + fl1_fx * tg_yz_yyyy_0[j] - fl1_fx * fl1_fza * tg_yz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yyz_yyy_1[j];

                    tg_yyyz_yyyz_0[j] = pb_y * tg_yyz_yyyz_0[j] + wp_y[j] * tg_yyz_yyyz_1[j] + fl1_fx * tg_yz_yyyz_0[j] - fl1_fx * fl1_fza * tg_yz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yyz_yyz_1[j];

                    tg_yyyz_yyzz_0[j] = pb_y * tg_yyz_yyzz_0[j] + wp_y[j] * tg_yyz_yyzz_1[j] + fl1_fx * tg_yz_yyzz_0[j] - fl1_fx * fl1_fza * tg_yz_yyzz_1[j] + fl1_fxn * tg_yyz_yzz_1[j];

                    tg_yyyz_yzzz_0[j] = pb_y * tg_yyz_yzzz_0[j] + wp_y[j] * tg_yyz_yzzz_1[j] + fl1_fx * tg_yz_yzzz_0[j] - fl1_fx * fl1_fza * tg_yz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zzz_1[j];

                    tg_yyyz_zzzz_0[j] = pb_y * tg_yyz_zzzz_0[j] + wp_y[j] * tg_yyz_zzzz_1[j] + fl1_fx * tg_yz_zzzz_0[j] - fl1_fx * fl1_fza * tg_yz_zzzz_1[j];

                    tg_yyzz_xxxx_0[j] = pb_y * tg_yzz_xxxx_0[j] + wp_y[j] * tg_yzz_xxxx_1[j] + 0.5 * fl1_fx * tg_zz_xxxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxx_1[j];

                    tg_yyzz_xxxy_0[j] = pb_y * tg_yzz_xxxy_0[j] + wp_y[j] * tg_yzz_xxxy_1[j] + 0.5 * fl1_fx * tg_zz_xxxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxy_1[j] + 0.5 * fl1_fxn * tg_yzz_xxx_1[j];

                    tg_yyzz_xxxz_0[j] = pb_y * tg_yzz_xxxz_0[j] + wp_y[j] * tg_yzz_xxxz_1[j] + 0.5 * fl1_fx * tg_zz_xxxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxxz_1[j];

                    tg_yyzz_xxyy_0[j] = pb_y * tg_yzz_xxyy_0[j] + wp_y[j] * tg_yzz_xxyy_1[j] + 0.5 * fl1_fx * tg_zz_xxyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyy_1[j] + fl1_fxn * tg_yzz_xxy_1[j];

                    tg_yyzz_xxyz_0[j] = pb_y * tg_yzz_xxyz_0[j] + wp_y[j] * tg_yzz_xxyz_1[j] + 0.5 * fl1_fx * tg_zz_xxyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxyz_1[j] + 0.5 * fl1_fxn * tg_yzz_xxz_1[j];

                    tg_yyzz_xxzz_0[j] = pb_y * tg_yzz_xxzz_0[j] + wp_y[j] * tg_yzz_xxzz_1[j] + 0.5 * fl1_fx * tg_zz_xxzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xxzz_1[j];

                    tg_yyzz_xyyy_0[j] = pb_y * tg_yzz_xyyy_0[j] + wp_y[j] * tg_yzz_xyyy_1[j] + 0.5 * fl1_fx * tg_zz_xyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyy_1[j] + 1.5 * fl1_fxn * tg_yzz_xyy_1[j];

                    tg_yyzz_xyyz_0[j] = pb_y * tg_yzz_xyyz_0[j] + wp_y[j] * tg_yzz_xyyz_1[j] + 0.5 * fl1_fx * tg_zz_xyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyyz_1[j] + fl1_fxn * tg_yzz_xyz_1[j];

                    tg_yyzz_xyzz_0[j] = pb_y * tg_yzz_xyzz_0[j] + wp_y[j] * tg_yzz_xyzz_1[j] + 0.5 * fl1_fx * tg_zz_xyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xyzz_1[j] + 0.5 * fl1_fxn * tg_yzz_xzz_1[j];

                    tg_yyzz_xzzz_0[j] = pb_y * tg_yzz_xzzz_0[j] + wp_y[j] * tg_yzz_xzzz_1[j] + 0.5 * fl1_fx * tg_zz_xzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_xzzz_1[j];

                    tg_yyzz_yyyy_0[j] = pb_y * tg_yzz_yyyy_0[j] + wp_y[j] * tg_yzz_yyyy_1[j] + 0.5 * fl1_fx * tg_zz_yyyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyy_1[j] + 2.0 * fl1_fxn * tg_yzz_yyy_1[j];

                    tg_yyzz_yyyz_0[j] = pb_y * tg_yzz_yyyz_0[j] + wp_y[j] * tg_yzz_yyyz_1[j] + 0.5 * fl1_fx * tg_zz_yyyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyyz_1[j] + 1.5 * fl1_fxn * tg_yzz_yyz_1[j];

                    tg_yyzz_yyzz_0[j] = pb_y * tg_yzz_yyzz_0[j] + wp_y[j] * tg_yzz_yyzz_1[j] + 0.5 * fl1_fx * tg_zz_yyzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yyzz_1[j] + fl1_fxn * tg_yzz_yzz_1[j];

                    tg_yyzz_yzzz_0[j] = pb_y * tg_yzz_yzzz_0[j] + wp_y[j] * tg_yzz_yzzz_1[j] + 0.5 * fl1_fx * tg_zz_yzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_yzzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zzz_1[j];

                    tg_yyzz_zzzz_0[j] = pb_y * tg_yzz_zzzz_0[j] + wp_y[j] * tg_yzz_zzzz_1[j] + 0.5 * fl1_fx * tg_zz_zzzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_zzzz_1[j];

                    tg_yzzz_xxxx_0[j] = pb_y * tg_zzz_xxxx_0[j] + wp_y[j] * tg_zzz_xxxx_1[j];

                    tg_yzzz_xxxy_0[j] = pb_y * tg_zzz_xxxy_0[j] + wp_y[j] * tg_zzz_xxxy_1[j] + 0.5 * fl1_fxn * tg_zzz_xxx_1[j];

                    tg_yzzz_xxxz_0[j] = pb_y * tg_zzz_xxxz_0[j] + wp_y[j] * tg_zzz_xxxz_1[j];

                    tg_yzzz_xxyy_0[j] = pb_y * tg_zzz_xxyy_0[j] + wp_y[j] * tg_zzz_xxyy_1[j] + fl1_fxn * tg_zzz_xxy_1[j];

                    tg_yzzz_xxyz_0[j] = pb_y * tg_zzz_xxyz_0[j] + wp_y[j] * tg_zzz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxz_1[j];

                    tg_yzzz_xxzz_0[j] = pb_y * tg_zzz_xxzz_0[j] + wp_y[j] * tg_zzz_xxzz_1[j];

                    tg_yzzz_xyyy_0[j] = pb_y * tg_zzz_xyyy_0[j] + wp_y[j] * tg_zzz_xyyy_1[j] + 1.5 * fl1_fxn * tg_zzz_xyy_1[j];

                    tg_yzzz_xyyz_0[j] = pb_y * tg_zzz_xyyz_0[j] + wp_y[j] * tg_zzz_xyyz_1[j] + fl1_fxn * tg_zzz_xyz_1[j];

                    tg_yzzz_xyzz_0[j] = pb_y * tg_zzz_xyzz_0[j] + wp_y[j] * tg_zzz_xyzz_1[j] + 0.5 * fl1_fxn * tg_zzz_xzz_1[j];

                    tg_yzzz_xzzz_0[j] = pb_y * tg_zzz_xzzz_0[j] + wp_y[j] * tg_zzz_xzzz_1[j];

                    tg_yzzz_yyyy_0[j] = pb_y * tg_zzz_yyyy_0[j] + wp_y[j] * tg_zzz_yyyy_1[j] + 2.0 * fl1_fxn * tg_zzz_yyy_1[j];

                    tg_yzzz_yyyz_0[j] = pb_y * tg_zzz_yyyz_0[j] + wp_y[j] * tg_zzz_yyyz_1[j] + 1.5 * fl1_fxn * tg_zzz_yyz_1[j];

                    tg_yzzz_yyzz_0[j] = pb_y * tg_zzz_yyzz_0[j] + wp_y[j] * tg_zzz_yyzz_1[j] + fl1_fxn * tg_zzz_yzz_1[j];

                    tg_yzzz_yzzz_0[j] = pb_y * tg_zzz_yzzz_0[j] + wp_y[j] * tg_zzz_yzzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zzz_1[j];

                    tg_yzzz_zzzz_0[j] = pb_y * tg_zzz_zzzz_0[j] + wp_y[j] * tg_zzz_zzzz_1[j];

                    tg_zzzz_xxxx_0[j] = pb_z * tg_zzz_xxxx_0[j] + wp_z[j] * tg_zzz_xxxx_1[j] + 1.5 * fl1_fx * tg_zz_xxxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxx_1[j];

                    tg_zzzz_xxxy_0[j] = pb_z * tg_zzz_xxxy_0[j] + wp_z[j] * tg_zzz_xxxy_1[j] + 1.5 * fl1_fx * tg_zz_xxxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxy_1[j];

                    tg_zzzz_xxxz_0[j] = pb_z * tg_zzz_xxxz_0[j] + wp_z[j] * tg_zzz_xxxz_1[j] + 1.5 * fl1_fx * tg_zz_xxxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxxz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxx_1[j];

                    tg_zzzz_xxyy_0[j] = pb_z * tg_zzz_xxyy_0[j] + wp_z[j] * tg_zzz_xxyy_1[j] + 1.5 * fl1_fx * tg_zz_xxyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyy_1[j];

                    tg_zzzz_xxyz_0[j] = pb_z * tg_zzz_xxyz_0[j] + wp_z[j] * tg_zzz_xxyz_1[j] + 1.5 * fl1_fx * tg_zz_xxyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xxy_1[j];

                    tg_zzzz_xxzz_0[j] = pb_z * tg_zzz_xxzz_0[j] + wp_z[j] * tg_zzz_xxzz_1[j] + 1.5 * fl1_fx * tg_zz_xxzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xxzz_1[j] + fl1_fxn * tg_zzz_xxz_1[j];

                    tg_zzzz_xyyy_0[j] = pb_z * tg_zzz_xyyy_0[j] + wp_z[j] * tg_zzz_xyyy_1[j] + 1.5 * fl1_fx * tg_zz_xyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyy_1[j];

                    tg_zzzz_xyyz_0[j] = pb_z * tg_zzz_xyyz_0[j] + wp_z[j] * tg_zzz_xyyz_1[j] + 1.5 * fl1_fx * tg_zz_xyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xyy_1[j];

                    tg_zzzz_xyzz_0[j] = pb_z * tg_zzz_xyzz_0[j] + wp_z[j] * tg_zzz_xyzz_1[j] + 1.5 * fl1_fx * tg_zz_xyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xyzz_1[j] + fl1_fxn * tg_zzz_xyz_1[j];

                    tg_zzzz_xzzz_0[j] = pb_z * tg_zzz_xzzz_0[j] + wp_z[j] * tg_zzz_xzzz_1[j] + 1.5 * fl1_fx * tg_zz_xzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_xzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_xzz_1[j];

                    tg_zzzz_yyyy_0[j] = pb_z * tg_zzz_yyyy_0[j] + wp_z[j] * tg_zzz_yyyy_1[j] + 1.5 * fl1_fx * tg_zz_yyyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyy_1[j];

                    tg_zzzz_yyyz_0[j] = pb_z * tg_zzz_yyyz_0[j] + wp_z[j] * tg_zzz_yyyz_1[j] + 1.5 * fl1_fx * tg_zz_yyyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yyy_1[j];

                    tg_zzzz_yyzz_0[j] = pb_z * tg_zzz_yyzz_0[j] + wp_z[j] * tg_zzz_yyzz_1[j] + 1.5 * fl1_fx * tg_zz_yyzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yyzz_1[j] + fl1_fxn * tg_zzz_yyz_1[j];

                    tg_zzzz_yzzz_0[j] = pb_z * tg_zzz_yzzz_0[j] + wp_z[j] * tg_zzz_yzzz_1[j] + 1.5 * fl1_fx * tg_zz_yzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_yzzz_1[j] + 1.5 * fl1_fxn * tg_zzz_yzz_1[j];

                    tg_zzzz_zzzz_0[j] = pb_z * tg_zzz_zzzz_0[j] + wp_z[j] * tg_zzz_zzzz_1[j] + 1.5 * fl1_fx * tg_zz_zzzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_zzzz_1[j] + 2.0 * fl1_fxn * tg_zzz_zzz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

