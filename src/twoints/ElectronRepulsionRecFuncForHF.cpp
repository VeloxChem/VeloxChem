//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForHF.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSHSF(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSHSF_0_70(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSHSF_70_140(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSHSF_140_210(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSHSF_0_70(      CMemBlock2D<double>& primBuffer,
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
                                             {5, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xxxx_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx); 

                auto tg_xxxx_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 1); 

                auto tg_xxxx_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 2); 

                auto tg_xxxx_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 3); 

                auto tg_xxxx_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 4); 

                auto tg_xxxx_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 5); 

                auto tg_xxxx_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 6); 

                auto tg_xxxx_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 7); 

                auto tg_xxxx_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 8); 

                auto tg_xxxx_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 9); 

                auto tg_xxxy_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 10); 

                auto tg_xxxy_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 11); 

                auto tg_xxxy_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 12); 

                auto tg_xxxy_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 13); 

                auto tg_xxxy_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 14); 

                auto tg_xxxy_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 15); 

                auto tg_xxxy_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 16); 

                auto tg_xxxy_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 17); 

                auto tg_xxxy_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 18); 

                auto tg_xxxy_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 19); 

                auto tg_xxxz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 20); 

                auto tg_xxxz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 21); 

                auto tg_xxxz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 22); 

                auto tg_xxxz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 23); 

                auto tg_xxxz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 24); 

                auto tg_xxxz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 25); 

                auto tg_xxxz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 26); 

                auto tg_xxxz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 27); 

                auto tg_xxxz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 28); 

                auto tg_xxxz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 29); 

                auto tg_xxyy_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 30); 

                auto tg_xxyy_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 31); 

                auto tg_xxyy_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 32); 

                auto tg_xxyy_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 33); 

                auto tg_xxyy_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 34); 

                auto tg_xxyy_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 35); 

                auto tg_xxyy_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 36); 

                auto tg_xxyy_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 37); 

                auto tg_xxyy_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 38); 

                auto tg_xxyy_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 39); 

                auto tg_xxyz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 40); 

                auto tg_xxyz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 41); 

                auto tg_xxyz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 42); 

                auto tg_xxyz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 43); 

                auto tg_xxyz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 44); 

                auto tg_xxyz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 45); 

                auto tg_xxyz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 46); 

                auto tg_xxyz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 47); 

                auto tg_xxyz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 48); 

                auto tg_xxyz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 49); 

                auto tg_xxzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 50); 

                auto tg_xxzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 51); 

                auto tg_xxzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 52); 

                auto tg_xxzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 53); 

                auto tg_xxzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 54); 

                auto tg_xxzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 55); 

                auto tg_xxzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 56); 

                auto tg_xxzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 57); 

                auto tg_xxzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 58); 

                auto tg_xxzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 59); 

                auto tg_xyyy_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 60); 

                auto tg_xyyy_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 61); 

                auto tg_xyyy_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 62); 

                auto tg_xyyy_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 63); 

                auto tg_xyyy_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 64); 

                auto tg_xyyy_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 65); 

                auto tg_xyyy_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 66); 

                auto tg_xyyy_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 67); 

                auto tg_xyyy_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 68); 

                auto tg_xyyy_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 69); 

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

                auto tg_xxx_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx); 

                auto tg_xxx_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 1); 

                auto tg_xxx_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 2); 

                auto tg_xxx_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 3); 

                auto tg_xxx_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 4); 

                auto tg_xxx_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 5); 

                auto tg_xxx_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 6); 

                auto tg_xxx_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 7); 

                auto tg_xxx_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 8); 

                auto tg_xxx_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 9); 

                auto tg_xxy_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 10); 

                auto tg_xxy_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 11); 

                auto tg_xxy_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 12); 

                auto tg_xxy_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 13); 

                auto tg_xxy_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 14); 

                auto tg_xxy_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 15); 

                auto tg_xxy_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 16); 

                auto tg_xxy_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 17); 

                auto tg_xxy_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 18); 

                auto tg_xxy_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 19); 

                auto tg_xxz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 20); 

                auto tg_xxz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 21); 

                auto tg_xxz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 22); 

                auto tg_xxz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 23); 

                auto tg_xxz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 24); 

                auto tg_xxz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 25); 

                auto tg_xxz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 26); 

                auto tg_xxz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 27); 

                auto tg_xxz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 28); 

                auto tg_xxz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 29); 

                auto tg_xyy_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 30); 

                auto tg_xyy_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 31); 

                auto tg_xyy_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 32); 

                auto tg_xyy_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 33); 

                auto tg_xyy_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 34); 

                auto tg_xyy_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 35); 

                auto tg_xyy_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 36); 

                auto tg_xyy_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 37); 

                auto tg_xyy_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 38); 

                auto tg_xyy_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 39); 

                auto tg_xyz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 40); 

                auto tg_xyz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 41); 

                auto tg_xyz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 42); 

                auto tg_xyz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 43); 

                auto tg_xyz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 44); 

                auto tg_xyz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 45); 

                auto tg_xyz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 46); 

                auto tg_xyz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 47); 

                auto tg_xyz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 48); 

                auto tg_xyz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 49); 

                auto tg_xzz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 50); 

                auto tg_xzz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 51); 

                auto tg_xzz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 52); 

                auto tg_xzz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 53); 

                auto tg_xzz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 54); 

                auto tg_xzz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 55); 

                auto tg_xzz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 56); 

                auto tg_xzz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 57); 

                auto tg_xzz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 58); 

                auto tg_xzz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 59); 

                auto tg_yyy_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 60); 

                auto tg_yyy_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 61); 

                auto tg_yyy_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 62); 

                auto tg_yyy_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 63); 

                auto tg_yyy_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 64); 

                auto tg_yyy_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 65); 

                auto tg_yyy_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 66); 

                auto tg_yyy_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 67); 

                auto tg_yyy_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 68); 

                auto tg_yyy_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 69); 

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

                auto tg_xxxx_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx); 

                auto tg_xxxx_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 1); 

                auto tg_xxxx_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 2); 

                auto tg_xxxx_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 3); 

                auto tg_xxxx_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 4); 

                auto tg_xxxx_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 5); 

                auto tg_xxxy_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 6); 

                auto tg_xxxy_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 7); 

                auto tg_xxxy_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 8); 

                auto tg_xxxy_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 9); 

                auto tg_xxxy_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 10); 

                auto tg_xxxy_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 11); 

                auto tg_xxxz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 12); 

                auto tg_xxxz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 13); 

                auto tg_xxxz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 14); 

                auto tg_xxxz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 15); 

                auto tg_xxxz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 16); 

                auto tg_xxxz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 17); 

                auto tg_xxyy_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 18); 

                auto tg_xxyy_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 19); 

                auto tg_xxyy_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 20); 

                auto tg_xxyy_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 21); 

                auto tg_xxyy_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 22); 

                auto tg_xxyy_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 23); 

                auto tg_xxyz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 24); 

                auto tg_xxyz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 25); 

                auto tg_xxyz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 26); 

                auto tg_xxyz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 27); 

                auto tg_xxyz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 28); 

                auto tg_xxyz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 29); 

                auto tg_xxzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 30); 

                auto tg_xxzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 31); 

                auto tg_xxzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 32); 

                auto tg_xxzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 33); 

                auto tg_xxzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 34); 

                auto tg_xxzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 35); 

                auto tg_xyyy_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 36); 

                auto tg_xyyy_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 37); 

                auto tg_xyyy_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 38); 

                auto tg_xyyy_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 39); 

                auto tg_xyyy_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 40); 

                auto tg_xyyy_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 41); 

                // set up pointers to integrals

                auto tg_xxxxx_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx); 

                auto tg_xxxxx_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 1); 

                auto tg_xxxxx_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 2); 

                auto tg_xxxxx_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 3); 

                auto tg_xxxxx_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 4); 

                auto tg_xxxxx_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 5); 

                auto tg_xxxxx_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 6); 

                auto tg_xxxxx_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 7); 

                auto tg_xxxxx_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 8); 

                auto tg_xxxxx_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 9); 

                auto tg_xxxxy_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 10); 

                auto tg_xxxxy_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 11); 

                auto tg_xxxxy_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 12); 

                auto tg_xxxxy_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 13); 

                auto tg_xxxxy_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 14); 

                auto tg_xxxxy_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 15); 

                auto tg_xxxxy_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 16); 

                auto tg_xxxxy_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 17); 

                auto tg_xxxxy_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 18); 

                auto tg_xxxxy_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 19); 

                auto tg_xxxxz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 20); 

                auto tg_xxxxz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 21); 

                auto tg_xxxxz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 22); 

                auto tg_xxxxz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 23); 

                auto tg_xxxxz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 24); 

                auto tg_xxxxz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 25); 

                auto tg_xxxxz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 26); 

                auto tg_xxxxz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 27); 

                auto tg_xxxxz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 28); 

                auto tg_xxxxz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 29); 

                auto tg_xxxyy_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 30); 

                auto tg_xxxyy_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 31); 

                auto tg_xxxyy_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 32); 

                auto tg_xxxyy_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 33); 

                auto tg_xxxyy_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 34); 

                auto tg_xxxyy_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 35); 

                auto tg_xxxyy_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 36); 

                auto tg_xxxyy_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 37); 

                auto tg_xxxyy_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 38); 

                auto tg_xxxyy_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 39); 

                auto tg_xxxyz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 40); 

                auto tg_xxxyz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 41); 

                auto tg_xxxyz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 42); 

                auto tg_xxxyz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 43); 

                auto tg_xxxyz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 44); 

                auto tg_xxxyz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 45); 

                auto tg_xxxyz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 46); 

                auto tg_xxxyz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 47); 

                auto tg_xxxyz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 48); 

                auto tg_xxxyz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 49); 

                auto tg_xxxzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 50); 

                auto tg_xxxzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 51); 

                auto tg_xxxzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 52); 

                auto tg_xxxzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 53); 

                auto tg_xxxzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 54); 

                auto tg_xxxzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 55); 

                auto tg_xxxzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 56); 

                auto tg_xxxzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 57); 

                auto tg_xxxzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 58); 

                auto tg_xxxzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 59); 

                auto tg_xxyyy_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 60); 

                auto tg_xxyyy_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 61); 

                auto tg_xxyyy_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 62); 

                auto tg_xxyyy_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 63); 

                auto tg_xxyyy_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 64); 

                auto tg_xxyyy_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 65); 

                auto tg_xxyyy_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 66); 

                auto tg_xxyyy_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 67); 

                auto tg_xxyyy_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 68); 

                auto tg_xxyyy_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 69); 

                // Batch of Integrals (0,70)

                #pragma omp simd aligned(fxn, fza, tg_xxx_xxx_0, tg_xxx_xxx_1, tg_xxx_xxy_0, tg_xxx_xxy_1, \
                                         tg_xxx_xxz_0, tg_xxx_xxz_1, tg_xxx_xyy_0, tg_xxx_xyy_1, tg_xxx_xyz_0, tg_xxx_xyz_1, \
                                         tg_xxx_xzz_0, tg_xxx_xzz_1, tg_xxx_yyy_0, tg_xxx_yyy_1, tg_xxx_yyz_0, tg_xxx_yyz_1, \
                                         tg_xxx_yzz_0, tg_xxx_yzz_1, tg_xxx_zzz_0, tg_xxx_zzz_1, tg_xxxx_xx_1, \
                                         tg_xxxx_xxx_0, tg_xxxx_xxx_1, tg_xxxx_xxy_0, tg_xxxx_xxy_1, tg_xxxx_xxz_0, \
                                         tg_xxxx_xxz_1, tg_xxxx_xy_1, tg_xxxx_xyy_0, tg_xxxx_xyy_1, tg_xxxx_xyz_0, \
                                         tg_xxxx_xyz_1, tg_xxxx_xz_1, tg_xxxx_xzz_0, tg_xxxx_xzz_1, tg_xxxx_yy_1, \
                                         tg_xxxx_yyy_0, tg_xxxx_yyy_1, tg_xxxx_yyz_0, tg_xxxx_yyz_1, tg_xxxx_yz_1, \
                                         tg_xxxx_yzz_0, tg_xxxx_yzz_1, tg_xxxx_zz_1, tg_xxxx_zzz_0, tg_xxxx_zzz_1, \
                                         tg_xxxxx_xxx_0, tg_xxxxx_xxy_0, tg_xxxxx_xxz_0, tg_xxxxx_xyy_0, tg_xxxxx_xyz_0, \
                                         tg_xxxxx_xzz_0, tg_xxxxx_yyy_0, tg_xxxxx_yyz_0, tg_xxxxx_yzz_0, tg_xxxxx_zzz_0, \
                                         tg_xxxxy_xxx_0, tg_xxxxy_xxy_0, tg_xxxxy_xxz_0, tg_xxxxy_xyy_0, tg_xxxxy_xyz_0, \
                                         tg_xxxxy_xzz_0, tg_xxxxy_yyy_0, tg_xxxxy_yyz_0, tg_xxxxy_yzz_0, tg_xxxxy_zzz_0, \
                                         tg_xxxxz_xxx_0, tg_xxxxz_xxy_0, tg_xxxxz_xxz_0, tg_xxxxz_xyy_0, tg_xxxxz_xyz_0, \
                                         tg_xxxxz_xzz_0, tg_xxxxz_yyy_0, tg_xxxxz_yyz_0, tg_xxxxz_yzz_0, tg_xxxxz_zzz_0, \
                                         tg_xxxy_xx_1, tg_xxxy_xxx_0, tg_xxxy_xxx_1, tg_xxxy_xxy_0, tg_xxxy_xxy_1, \
                                         tg_xxxy_xxz_0, tg_xxxy_xxz_1, tg_xxxy_xy_1, tg_xxxy_xyy_0, tg_xxxy_xyy_1, \
                                         tg_xxxy_xyz_0, tg_xxxy_xyz_1, tg_xxxy_xz_1, tg_xxxy_xzz_0, tg_xxxy_xzz_1, \
                                         tg_xxxy_yy_1, tg_xxxy_yyy_0, tg_xxxy_yyy_1, tg_xxxy_yyz_0, tg_xxxy_yyz_1, \
                                         tg_xxxy_yz_1, tg_xxxy_yzz_0, tg_xxxy_yzz_1, tg_xxxy_zz_1, tg_xxxy_zzz_0, \
                                         tg_xxxy_zzz_1, tg_xxxyy_xxx_0, tg_xxxyy_xxy_0, tg_xxxyy_xxz_0, tg_xxxyy_xyy_0, \
                                         tg_xxxyy_xyz_0, tg_xxxyy_xzz_0, tg_xxxyy_yyy_0, tg_xxxyy_yyz_0, tg_xxxyy_yzz_0, \
                                         tg_xxxyy_zzz_0, tg_xxxyz_xxx_0, tg_xxxyz_xxy_0, tg_xxxyz_xxz_0, tg_xxxyz_xyy_0, \
                                         tg_xxxyz_xyz_0, tg_xxxyz_xzz_0, tg_xxxyz_yyy_0, tg_xxxyz_yyz_0, tg_xxxyz_yzz_0, \
                                         tg_xxxyz_zzz_0, tg_xxxz_xx_1, tg_xxxz_xxx_0, tg_xxxz_xxx_1, tg_xxxz_xxy_0, \
                                         tg_xxxz_xxy_1, tg_xxxz_xxz_0, tg_xxxz_xxz_1, tg_xxxz_xy_1, tg_xxxz_xyy_0, \
                                         tg_xxxz_xyy_1, tg_xxxz_xyz_0, tg_xxxz_xyz_1, tg_xxxz_xz_1, tg_xxxz_xzz_0, \
                                         tg_xxxz_xzz_1, tg_xxxz_yy_1, tg_xxxz_yyy_0, tg_xxxz_yyy_1, tg_xxxz_yyz_0, \
                                         tg_xxxz_yyz_1, tg_xxxz_yz_1, tg_xxxz_yzz_0, tg_xxxz_yzz_1, tg_xxxz_zz_1, \
                                         tg_xxxz_zzz_0, tg_xxxz_zzz_1, tg_xxxzz_xxx_0, tg_xxxzz_xxy_0, tg_xxxzz_xxz_0, \
                                         tg_xxxzz_xyy_0, tg_xxxzz_xyz_0, tg_xxxzz_xzz_0, tg_xxxzz_yyy_0, tg_xxxzz_yyz_0, \
                                         tg_xxxzz_yzz_0, tg_xxxzz_zzz_0, tg_xxy_xxx_0, tg_xxy_xxx_1, tg_xxy_xxy_0, \
                                         tg_xxy_xxy_1, tg_xxy_xxz_0, tg_xxy_xxz_1, tg_xxy_xyy_0, tg_xxy_xyy_1, tg_xxy_xyz_0, \
                                         tg_xxy_xyz_1, tg_xxy_xzz_0, tg_xxy_xzz_1, tg_xxy_yyy_0, tg_xxy_yyy_1, tg_xxy_yyz_0, \
                                         tg_xxy_yyz_1, tg_xxy_yzz_0, tg_xxy_yzz_1, tg_xxy_zzz_0, tg_xxy_zzz_1, tg_xxyy_xx_1, \
                                         tg_xxyy_xxx_0, tg_xxyy_xxx_1, tg_xxyy_xxy_0, tg_xxyy_xxy_1, tg_xxyy_xxz_0, \
                                         tg_xxyy_xxz_1, tg_xxyy_xy_1, tg_xxyy_xyy_0, tg_xxyy_xyy_1, tg_xxyy_xyz_0, \
                                         tg_xxyy_xyz_1, tg_xxyy_xz_1, tg_xxyy_xzz_0, tg_xxyy_xzz_1, tg_xxyy_yy_1, \
                                         tg_xxyy_yyy_0, tg_xxyy_yyy_1, tg_xxyy_yyz_0, tg_xxyy_yyz_1, tg_xxyy_yz_1, \
                                         tg_xxyy_yzz_0, tg_xxyy_yzz_1, tg_xxyy_zz_1, tg_xxyy_zzz_0, tg_xxyy_zzz_1, \
                                         tg_xxyyy_xxx_0, tg_xxyyy_xxy_0, tg_xxyyy_xxz_0, tg_xxyyy_xyy_0, tg_xxyyy_xyz_0, \
                                         tg_xxyyy_xzz_0, tg_xxyyy_yyy_0, tg_xxyyy_yyz_0, tg_xxyyy_yzz_0, tg_xxyyy_zzz_0, \
                                         tg_xxyz_xx_1, tg_xxyz_xxx_0, tg_xxyz_xxx_1, tg_xxyz_xxy_0, tg_xxyz_xxy_1, \
                                         tg_xxyz_xxz_0, tg_xxyz_xxz_1, tg_xxyz_xy_1, tg_xxyz_xyy_0, tg_xxyz_xyy_1, \
                                         tg_xxyz_xyz_0, tg_xxyz_xyz_1, tg_xxyz_xz_1, tg_xxyz_xzz_0, tg_xxyz_xzz_1, \
                                         tg_xxyz_yy_1, tg_xxyz_yyy_0, tg_xxyz_yyy_1, tg_xxyz_yyz_0, tg_xxyz_yyz_1, \
                                         tg_xxyz_yz_1, tg_xxyz_yzz_0, tg_xxyz_yzz_1, tg_xxyz_zz_1, tg_xxyz_zzz_0, \
                                         tg_xxyz_zzz_1, tg_xxz_xxx_0, tg_xxz_xxx_1, tg_xxz_xxy_0, tg_xxz_xxy_1, tg_xxz_xxz_0, \
                                         tg_xxz_xxz_1, tg_xxz_xyy_0, tg_xxz_xyy_1, tg_xxz_xyz_0, tg_xxz_xyz_1, tg_xxz_xzz_0, \
                                         tg_xxz_xzz_1, tg_xxz_yyy_0, tg_xxz_yyy_1, tg_xxz_yyz_0, tg_xxz_yyz_1, tg_xxz_yzz_0, \
                                         tg_xxz_yzz_1, tg_xxz_zzz_0, tg_xxz_zzz_1, tg_xxzz_xx_1, tg_xxzz_xxx_0, \
                                         tg_xxzz_xxx_1, tg_xxzz_xxy_0, tg_xxzz_xxy_1, tg_xxzz_xxz_0, tg_xxzz_xxz_1, \
                                         tg_xxzz_xy_1, tg_xxzz_xyy_0, tg_xxzz_xyy_1, tg_xxzz_xyz_0, tg_xxzz_xyz_1, \
                                         tg_xxzz_xz_1, tg_xxzz_xzz_0, tg_xxzz_xzz_1, tg_xxzz_yy_1, tg_xxzz_yyy_0, \
                                         tg_xxzz_yyy_1, tg_xxzz_yyz_0, tg_xxzz_yyz_1, tg_xxzz_yz_1, tg_xxzz_yzz_0, \
                                         tg_xxzz_yzz_1, tg_xxzz_zz_1, tg_xxzz_zzz_0, tg_xxzz_zzz_1, tg_xyy_xxx_0, \
                                         tg_xyy_xxx_1, tg_xyy_xxy_0, tg_xyy_xxy_1, tg_xyy_xxz_0, tg_xyy_xxz_1, tg_xyy_xyy_0, \
                                         tg_xyy_xyy_1, tg_xyy_xyz_0, tg_xyy_xyz_1, tg_xyy_xzz_0, tg_xyy_xzz_1, tg_xyy_yyy_0, \
                                         tg_xyy_yyy_1, tg_xyy_yyz_0, tg_xyy_yyz_1, tg_xyy_yzz_0, tg_xyy_yzz_1, tg_xyy_zzz_0, \
                                         tg_xyy_zzz_1, tg_xyyy_xx_1, tg_xyyy_xxx_0, tg_xyyy_xxx_1, tg_xyyy_xxy_0, \
                                         tg_xyyy_xxy_1, tg_xyyy_xxz_0, tg_xyyy_xxz_1, tg_xyyy_xy_1, tg_xyyy_xyy_0, \
                                         tg_xyyy_xyy_1, tg_xyyy_xyz_0, tg_xyyy_xyz_1, tg_xyyy_xz_1, tg_xyyy_xzz_0, \
                                         tg_xyyy_xzz_1, tg_xyyy_yy_1, tg_xyyy_yyy_0, tg_xyyy_yyy_1, tg_xyyy_yyz_0, \
                                         tg_xyyy_yyz_1, tg_xyyy_yz_1, tg_xyyy_yzz_0, tg_xyyy_yzz_1, tg_xyyy_zz_1, \
                                         tg_xyyy_zzz_0, tg_xyyy_zzz_1, tg_xyz_xxx_0, tg_xyz_xxx_1, tg_xyz_xxy_0, \
                                         tg_xyz_xxy_1, tg_xyz_xxz_0, tg_xyz_xxz_1, tg_xyz_xyy_0, tg_xyz_xyy_1, tg_xyz_xyz_0, \
                                         tg_xyz_xyz_1, tg_xyz_xzz_0, tg_xyz_xzz_1, tg_xyz_yyy_0, tg_xyz_yyy_1, tg_xyz_yyz_0, \
                                         tg_xyz_yyz_1, tg_xyz_yzz_0, tg_xyz_yzz_1, tg_xyz_zzz_0, tg_xyz_zzz_1, tg_xzz_xxx_0, \
                                         tg_xzz_xxx_1, tg_xzz_xxy_0, tg_xzz_xxy_1, tg_xzz_xxz_0, tg_xzz_xxz_1, tg_xzz_xyy_0, \
                                         tg_xzz_xyy_1, tg_xzz_xyz_0, tg_xzz_xyz_1, tg_xzz_xzz_0, tg_xzz_xzz_1, tg_xzz_yyy_0, \
                                         tg_xzz_yyy_1, tg_xzz_yyz_0, tg_xzz_yyz_1, tg_xzz_yzz_0, tg_xzz_yzz_1, tg_xzz_zzz_0, \
                                         tg_xzz_zzz_1, tg_yyy_xxx_0, tg_yyy_xxx_1, tg_yyy_xxy_0, tg_yyy_xxy_1, tg_yyy_xxz_0, \
                                         tg_yyy_xxz_1, tg_yyy_xyy_0, tg_yyy_xyy_1, tg_yyy_xyz_0, tg_yyy_xyz_1, tg_yyy_xzz_0, \
                                         tg_yyy_xzz_1, tg_yyy_yyy_0, tg_yyy_yyy_1, tg_yyy_yyz_0, tg_yyy_yyz_1, tg_yyy_yzz_0, \
                                         tg_yyy_yzz_1, tg_yyy_zzz_0, tg_yyy_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxx_xxx_0[j] = pb_x * tg_xxxx_xxx_0[j] + wp_x[j] * tg_xxxx_xxx_1[j] + 2.0 * fl1_fx * tg_xxx_xxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxx_1[j] + 1.5 * fl1_fxn * tg_xxxx_xx_1[j];

                    tg_xxxxx_xxy_0[j] = pb_x * tg_xxxx_xxy_0[j] + wp_x[j] * tg_xxxx_xxy_1[j] + 2.0 * fl1_fx * tg_xxx_xxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxy_1[j] + fl1_fxn * tg_xxxx_xy_1[j];

                    tg_xxxxx_xxz_0[j] = pb_x * tg_xxxx_xxz_0[j] + wp_x[j] * tg_xxxx_xxz_1[j] + 2.0 * fl1_fx * tg_xxx_xxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xxz_1[j] + fl1_fxn * tg_xxxx_xz_1[j];

                    tg_xxxxx_xyy_0[j] = pb_x * tg_xxxx_xyy_0[j] + wp_x[j] * tg_xxxx_xyy_1[j] + 2.0 * fl1_fx * tg_xxx_xyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyy_1[j] + 0.5 * fl1_fxn * tg_xxxx_yy_1[j];

                    tg_xxxxx_xyz_0[j] = pb_x * tg_xxxx_xyz_0[j] + wp_x[j] * tg_xxxx_xyz_1[j] + 2.0 * fl1_fx * tg_xxx_xyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xyz_1[j] + 0.5 * fl1_fxn * tg_xxxx_yz_1[j];

                    tg_xxxxx_xzz_0[j] = pb_x * tg_xxxx_xzz_0[j] + wp_x[j] * tg_xxxx_xzz_1[j] + 2.0 * fl1_fx * tg_xxx_xzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_xzz_1[j] + 0.5 * fl1_fxn * tg_xxxx_zz_1[j];

                    tg_xxxxx_yyy_0[j] = pb_x * tg_xxxx_yyy_0[j] + wp_x[j] * tg_xxxx_yyy_1[j] + 2.0 * fl1_fx * tg_xxx_yyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyy_1[j];

                    tg_xxxxx_yyz_0[j] = pb_x * tg_xxxx_yyz_0[j] + wp_x[j] * tg_xxxx_yyz_1[j] + 2.0 * fl1_fx * tg_xxx_yyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yyz_1[j];

                    tg_xxxxx_yzz_0[j] = pb_x * tg_xxxx_yzz_0[j] + wp_x[j] * tg_xxxx_yzz_1[j] + 2.0 * fl1_fx * tg_xxx_yzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_yzz_1[j];

                    tg_xxxxx_zzz_0[j] = pb_x * tg_xxxx_zzz_0[j] + wp_x[j] * tg_xxxx_zzz_1[j] + 2.0 * fl1_fx * tg_xxx_zzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_zzz_1[j];

                    tg_xxxxy_xxx_0[j] = pb_x * tg_xxxy_xxx_0[j] + wp_x[j] * tg_xxxy_xxx_1[j] + 1.5 * fl1_fx * tg_xxy_xxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxx_1[j] + 1.5 * fl1_fxn * tg_xxxy_xx_1[j];

                    tg_xxxxy_xxy_0[j] = pb_x * tg_xxxy_xxy_0[j] + wp_x[j] * tg_xxxy_xxy_1[j] + 1.5 * fl1_fx * tg_xxy_xxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxy_1[j] + fl1_fxn * tg_xxxy_xy_1[j];

                    tg_xxxxy_xxz_0[j] = pb_x * tg_xxxy_xxz_0[j] + wp_x[j] * tg_xxxy_xxz_1[j] + 1.5 * fl1_fx * tg_xxy_xxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xxz_1[j] + fl1_fxn * tg_xxxy_xz_1[j];

                    tg_xxxxy_xyy_0[j] = pb_x * tg_xxxy_xyy_0[j] + wp_x[j] * tg_xxxy_xyy_1[j] + 1.5 * fl1_fx * tg_xxy_xyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyy_1[j] + 0.5 * fl1_fxn * tg_xxxy_yy_1[j];

                    tg_xxxxy_xyz_0[j] = pb_x * tg_xxxy_xyz_0[j] + wp_x[j] * tg_xxxy_xyz_1[j] + 1.5 * fl1_fx * tg_xxy_xyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xyz_1[j] + 0.5 * fl1_fxn * tg_xxxy_yz_1[j];

                    tg_xxxxy_xzz_0[j] = pb_x * tg_xxxy_xzz_0[j] + wp_x[j] * tg_xxxy_xzz_1[j] + 1.5 * fl1_fx * tg_xxy_xzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_xzz_1[j] + 0.5 * fl1_fxn * tg_xxxy_zz_1[j];

                    tg_xxxxy_yyy_0[j] = pb_x * tg_xxxy_yyy_0[j] + wp_x[j] * tg_xxxy_yyy_1[j] + 1.5 * fl1_fx * tg_xxy_yyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyy_1[j];

                    tg_xxxxy_yyz_0[j] = pb_x * tg_xxxy_yyz_0[j] + wp_x[j] * tg_xxxy_yyz_1[j] + 1.5 * fl1_fx * tg_xxy_yyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yyz_1[j];

                    tg_xxxxy_yzz_0[j] = pb_x * tg_xxxy_yzz_0[j] + wp_x[j] * tg_xxxy_yzz_1[j] + 1.5 * fl1_fx * tg_xxy_yzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_yzz_1[j];

                    tg_xxxxy_zzz_0[j] = pb_x * tg_xxxy_zzz_0[j] + wp_x[j] * tg_xxxy_zzz_1[j] + 1.5 * fl1_fx * tg_xxy_zzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_zzz_1[j];

                    tg_xxxxz_xxx_0[j] = pb_x * tg_xxxz_xxx_0[j] + wp_x[j] * tg_xxxz_xxx_1[j] + 1.5 * fl1_fx * tg_xxz_xxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxx_1[j] + 1.5 * fl1_fxn * tg_xxxz_xx_1[j];

                    tg_xxxxz_xxy_0[j] = pb_x * tg_xxxz_xxy_0[j] + wp_x[j] * tg_xxxz_xxy_1[j] + 1.5 * fl1_fx * tg_xxz_xxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxy_1[j] + fl1_fxn * tg_xxxz_xy_1[j];

                    tg_xxxxz_xxz_0[j] = pb_x * tg_xxxz_xxz_0[j] + wp_x[j] * tg_xxxz_xxz_1[j] + 1.5 * fl1_fx * tg_xxz_xxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xxz_1[j] + fl1_fxn * tg_xxxz_xz_1[j];

                    tg_xxxxz_xyy_0[j] = pb_x * tg_xxxz_xyy_0[j] + wp_x[j] * tg_xxxz_xyy_1[j] + 1.5 * fl1_fx * tg_xxz_xyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyy_1[j] + 0.5 * fl1_fxn * tg_xxxz_yy_1[j];

                    tg_xxxxz_xyz_0[j] = pb_x * tg_xxxz_xyz_0[j] + wp_x[j] * tg_xxxz_xyz_1[j] + 1.5 * fl1_fx * tg_xxz_xyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xyz_1[j] + 0.5 * fl1_fxn * tg_xxxz_yz_1[j];

                    tg_xxxxz_xzz_0[j] = pb_x * tg_xxxz_xzz_0[j] + wp_x[j] * tg_xxxz_xzz_1[j] + 1.5 * fl1_fx * tg_xxz_xzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_xzz_1[j] + 0.5 * fl1_fxn * tg_xxxz_zz_1[j];

                    tg_xxxxz_yyy_0[j] = pb_x * tg_xxxz_yyy_0[j] + wp_x[j] * tg_xxxz_yyy_1[j] + 1.5 * fl1_fx * tg_xxz_yyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyy_1[j];

                    tg_xxxxz_yyz_0[j] = pb_x * tg_xxxz_yyz_0[j] + wp_x[j] * tg_xxxz_yyz_1[j] + 1.5 * fl1_fx * tg_xxz_yyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yyz_1[j];

                    tg_xxxxz_yzz_0[j] = pb_x * tg_xxxz_yzz_0[j] + wp_x[j] * tg_xxxz_yzz_1[j] + 1.5 * fl1_fx * tg_xxz_yzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_yzz_1[j];

                    tg_xxxxz_zzz_0[j] = pb_x * tg_xxxz_zzz_0[j] + wp_x[j] * tg_xxxz_zzz_1[j] + 1.5 * fl1_fx * tg_xxz_zzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_zzz_1[j];

                    tg_xxxyy_xxx_0[j] = pb_x * tg_xxyy_xxx_0[j] + wp_x[j] * tg_xxyy_xxx_1[j] + fl1_fx * tg_xyy_xxx_0[j] - fl1_fx * fl1_fza * tg_xyy_xxx_1[j] + 1.5 * fl1_fxn * tg_xxyy_xx_1[j];

                    tg_xxxyy_xxy_0[j] = pb_x * tg_xxyy_xxy_0[j] + wp_x[j] * tg_xxyy_xxy_1[j] + fl1_fx * tg_xyy_xxy_0[j] - fl1_fx * fl1_fza * tg_xyy_xxy_1[j] + fl1_fxn * tg_xxyy_xy_1[j];

                    tg_xxxyy_xxz_0[j] = pb_x * tg_xxyy_xxz_0[j] + wp_x[j] * tg_xxyy_xxz_1[j] + fl1_fx * tg_xyy_xxz_0[j] - fl1_fx * fl1_fza * tg_xyy_xxz_1[j] + fl1_fxn * tg_xxyy_xz_1[j];

                    tg_xxxyy_xyy_0[j] = pb_x * tg_xxyy_xyy_0[j] + wp_x[j] * tg_xxyy_xyy_1[j] + fl1_fx * tg_xyy_xyy_0[j] - fl1_fx * fl1_fza * tg_xyy_xyy_1[j] + 0.5 * fl1_fxn * tg_xxyy_yy_1[j];

                    tg_xxxyy_xyz_0[j] = pb_x * tg_xxyy_xyz_0[j] + wp_x[j] * tg_xxyy_xyz_1[j] + fl1_fx * tg_xyy_xyz_0[j] - fl1_fx * fl1_fza * tg_xyy_xyz_1[j] + 0.5 * fl1_fxn * tg_xxyy_yz_1[j];

                    tg_xxxyy_xzz_0[j] = pb_x * tg_xxyy_xzz_0[j] + wp_x[j] * tg_xxyy_xzz_1[j] + fl1_fx * tg_xyy_xzz_0[j] - fl1_fx * fl1_fza * tg_xyy_xzz_1[j] + 0.5 * fl1_fxn * tg_xxyy_zz_1[j];

                    tg_xxxyy_yyy_0[j] = pb_x * tg_xxyy_yyy_0[j] + wp_x[j] * tg_xxyy_yyy_1[j] + fl1_fx * tg_xyy_yyy_0[j] - fl1_fx * fl1_fza * tg_xyy_yyy_1[j];

                    tg_xxxyy_yyz_0[j] = pb_x * tg_xxyy_yyz_0[j] + wp_x[j] * tg_xxyy_yyz_1[j] + fl1_fx * tg_xyy_yyz_0[j] - fl1_fx * fl1_fza * tg_xyy_yyz_1[j];

                    tg_xxxyy_yzz_0[j] = pb_x * tg_xxyy_yzz_0[j] + wp_x[j] * tg_xxyy_yzz_1[j] + fl1_fx * tg_xyy_yzz_0[j] - fl1_fx * fl1_fza * tg_xyy_yzz_1[j];

                    tg_xxxyy_zzz_0[j] = pb_x * tg_xxyy_zzz_0[j] + wp_x[j] * tg_xxyy_zzz_1[j] + fl1_fx * tg_xyy_zzz_0[j] - fl1_fx * fl1_fza * tg_xyy_zzz_1[j];

                    tg_xxxyz_xxx_0[j] = pb_x * tg_xxyz_xxx_0[j] + wp_x[j] * tg_xxyz_xxx_1[j] + fl1_fx * tg_xyz_xxx_0[j] - fl1_fx * fl1_fza * tg_xyz_xxx_1[j] + 1.5 * fl1_fxn * tg_xxyz_xx_1[j];

                    tg_xxxyz_xxy_0[j] = pb_x * tg_xxyz_xxy_0[j] + wp_x[j] * tg_xxyz_xxy_1[j] + fl1_fx * tg_xyz_xxy_0[j] - fl1_fx * fl1_fza * tg_xyz_xxy_1[j] + fl1_fxn * tg_xxyz_xy_1[j];

                    tg_xxxyz_xxz_0[j] = pb_x * tg_xxyz_xxz_0[j] + wp_x[j] * tg_xxyz_xxz_1[j] + fl1_fx * tg_xyz_xxz_0[j] - fl1_fx * fl1_fza * tg_xyz_xxz_1[j] + fl1_fxn * tg_xxyz_xz_1[j];

                    tg_xxxyz_xyy_0[j] = pb_x * tg_xxyz_xyy_0[j] + wp_x[j] * tg_xxyz_xyy_1[j] + fl1_fx * tg_xyz_xyy_0[j] - fl1_fx * fl1_fza * tg_xyz_xyy_1[j] + 0.5 * fl1_fxn * tg_xxyz_yy_1[j];

                    tg_xxxyz_xyz_0[j] = pb_x * tg_xxyz_xyz_0[j] + wp_x[j] * tg_xxyz_xyz_1[j] + fl1_fx * tg_xyz_xyz_0[j] - fl1_fx * fl1_fza * tg_xyz_xyz_1[j] + 0.5 * fl1_fxn * tg_xxyz_yz_1[j];

                    tg_xxxyz_xzz_0[j] = pb_x * tg_xxyz_xzz_0[j] + wp_x[j] * tg_xxyz_xzz_1[j] + fl1_fx * tg_xyz_xzz_0[j] - fl1_fx * fl1_fza * tg_xyz_xzz_1[j] + 0.5 * fl1_fxn * tg_xxyz_zz_1[j];

                    tg_xxxyz_yyy_0[j] = pb_x * tg_xxyz_yyy_0[j] + wp_x[j] * tg_xxyz_yyy_1[j] + fl1_fx * tg_xyz_yyy_0[j] - fl1_fx * fl1_fza * tg_xyz_yyy_1[j];

                    tg_xxxyz_yyz_0[j] = pb_x * tg_xxyz_yyz_0[j] + wp_x[j] * tg_xxyz_yyz_1[j] + fl1_fx * tg_xyz_yyz_0[j] - fl1_fx * fl1_fza * tg_xyz_yyz_1[j];

                    tg_xxxyz_yzz_0[j] = pb_x * tg_xxyz_yzz_0[j] + wp_x[j] * tg_xxyz_yzz_1[j] + fl1_fx * tg_xyz_yzz_0[j] - fl1_fx * fl1_fza * tg_xyz_yzz_1[j];

                    tg_xxxyz_zzz_0[j] = pb_x * tg_xxyz_zzz_0[j] + wp_x[j] * tg_xxyz_zzz_1[j] + fl1_fx * tg_xyz_zzz_0[j] - fl1_fx * fl1_fza * tg_xyz_zzz_1[j];

                    tg_xxxzz_xxx_0[j] = pb_x * tg_xxzz_xxx_0[j] + wp_x[j] * tg_xxzz_xxx_1[j] + fl1_fx * tg_xzz_xxx_0[j] - fl1_fx * fl1_fza * tg_xzz_xxx_1[j] + 1.5 * fl1_fxn * tg_xxzz_xx_1[j];

                    tg_xxxzz_xxy_0[j] = pb_x * tg_xxzz_xxy_0[j] + wp_x[j] * tg_xxzz_xxy_1[j] + fl1_fx * tg_xzz_xxy_0[j] - fl1_fx * fl1_fza * tg_xzz_xxy_1[j] + fl1_fxn * tg_xxzz_xy_1[j];

                    tg_xxxzz_xxz_0[j] = pb_x * tg_xxzz_xxz_0[j] + wp_x[j] * tg_xxzz_xxz_1[j] + fl1_fx * tg_xzz_xxz_0[j] - fl1_fx * fl1_fza * tg_xzz_xxz_1[j] + fl1_fxn * tg_xxzz_xz_1[j];

                    tg_xxxzz_xyy_0[j] = pb_x * tg_xxzz_xyy_0[j] + wp_x[j] * tg_xxzz_xyy_1[j] + fl1_fx * tg_xzz_xyy_0[j] - fl1_fx * fl1_fza * tg_xzz_xyy_1[j] + 0.5 * fl1_fxn * tg_xxzz_yy_1[j];

                    tg_xxxzz_xyz_0[j] = pb_x * tg_xxzz_xyz_0[j] + wp_x[j] * tg_xxzz_xyz_1[j] + fl1_fx * tg_xzz_xyz_0[j] - fl1_fx * fl1_fza * tg_xzz_xyz_1[j] + 0.5 * fl1_fxn * tg_xxzz_yz_1[j];

                    tg_xxxzz_xzz_0[j] = pb_x * tg_xxzz_xzz_0[j] + wp_x[j] * tg_xxzz_xzz_1[j] + fl1_fx * tg_xzz_xzz_0[j] - fl1_fx * fl1_fza * tg_xzz_xzz_1[j] + 0.5 * fl1_fxn * tg_xxzz_zz_1[j];

                    tg_xxxzz_yyy_0[j] = pb_x * tg_xxzz_yyy_0[j] + wp_x[j] * tg_xxzz_yyy_1[j] + fl1_fx * tg_xzz_yyy_0[j] - fl1_fx * fl1_fza * tg_xzz_yyy_1[j];

                    tg_xxxzz_yyz_0[j] = pb_x * tg_xxzz_yyz_0[j] + wp_x[j] * tg_xxzz_yyz_1[j] + fl1_fx * tg_xzz_yyz_0[j] - fl1_fx * fl1_fza * tg_xzz_yyz_1[j];

                    tg_xxxzz_yzz_0[j] = pb_x * tg_xxzz_yzz_0[j] + wp_x[j] * tg_xxzz_yzz_1[j] + fl1_fx * tg_xzz_yzz_0[j] - fl1_fx * fl1_fza * tg_xzz_yzz_1[j];

                    tg_xxxzz_zzz_0[j] = pb_x * tg_xxzz_zzz_0[j] + wp_x[j] * tg_xxzz_zzz_1[j] + fl1_fx * tg_xzz_zzz_0[j] - fl1_fx * fl1_fza * tg_xzz_zzz_1[j];

                    tg_xxyyy_xxx_0[j] = pb_x * tg_xyyy_xxx_0[j] + wp_x[j] * tg_xyyy_xxx_1[j] + 0.5 * fl1_fx * tg_yyy_xxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxx_1[j] + 1.5 * fl1_fxn * tg_xyyy_xx_1[j];

                    tg_xxyyy_xxy_0[j] = pb_x * tg_xyyy_xxy_0[j] + wp_x[j] * tg_xyyy_xxy_1[j] + 0.5 * fl1_fx * tg_yyy_xxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxy_1[j] + fl1_fxn * tg_xyyy_xy_1[j];

                    tg_xxyyy_xxz_0[j] = pb_x * tg_xyyy_xxz_0[j] + wp_x[j] * tg_xyyy_xxz_1[j] + 0.5 * fl1_fx * tg_yyy_xxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xxz_1[j] + fl1_fxn * tg_xyyy_xz_1[j];

                    tg_xxyyy_xyy_0[j] = pb_x * tg_xyyy_xyy_0[j] + wp_x[j] * tg_xyyy_xyy_1[j] + 0.5 * fl1_fx * tg_yyy_xyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyy_1[j] + 0.5 * fl1_fxn * tg_xyyy_yy_1[j];

                    tg_xxyyy_xyz_0[j] = pb_x * tg_xyyy_xyz_0[j] + wp_x[j] * tg_xyyy_xyz_1[j] + 0.5 * fl1_fx * tg_yyy_xyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xyz_1[j] + 0.5 * fl1_fxn * tg_xyyy_yz_1[j];

                    tg_xxyyy_xzz_0[j] = pb_x * tg_xyyy_xzz_0[j] + wp_x[j] * tg_xyyy_xzz_1[j] + 0.5 * fl1_fx * tg_yyy_xzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_xzz_1[j] + 0.5 * fl1_fxn * tg_xyyy_zz_1[j];

                    tg_xxyyy_yyy_0[j] = pb_x * tg_xyyy_yyy_0[j] + wp_x[j] * tg_xyyy_yyy_1[j] + 0.5 * fl1_fx * tg_yyy_yyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyy_1[j];

                    tg_xxyyy_yyz_0[j] = pb_x * tg_xyyy_yyz_0[j] + wp_x[j] * tg_xyyy_yyz_1[j] + 0.5 * fl1_fx * tg_yyy_yyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yyz_1[j];

                    tg_xxyyy_yzz_0[j] = pb_x * tg_xyyy_yzz_0[j] + wp_x[j] * tg_xyyy_yzz_1[j] + 0.5 * fl1_fx * tg_yyy_yzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_yzz_1[j];

                    tg_xxyyy_zzz_0[j] = pb_x * tg_xyyy_zzz_0[j] + wp_x[j] * tg_xyyy_zzz_1[j] + 0.5 * fl1_fx * tg_yyy_zzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSF_70_140(      CMemBlock2D<double>& primBuffer,
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

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xyyz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 70); 

                auto tg_xyyz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 71); 

                auto tg_xyyz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 72); 

                auto tg_xyyz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 73); 

                auto tg_xyyz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 74); 

                auto tg_xyyz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 75); 

                auto tg_xyyz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 76); 

                auto tg_xyyz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 77); 

                auto tg_xyyz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 78); 

                auto tg_xyyz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 79); 

                auto tg_xyzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 80); 

                auto tg_xyzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 81); 

                auto tg_xyzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 82); 

                auto tg_xyzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 83); 

                auto tg_xyzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 84); 

                auto tg_xyzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 85); 

                auto tg_xyzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 86); 

                auto tg_xyzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 87); 

                auto tg_xyzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 88); 

                auto tg_xyzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 89); 

                auto tg_xzzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 90); 

                auto tg_xzzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 91); 

                auto tg_xzzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 92); 

                auto tg_xzzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 93); 

                auto tg_xzzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 94); 

                auto tg_xzzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 95); 

                auto tg_xzzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 96); 

                auto tg_xzzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 97); 

                auto tg_xzzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 98); 

                auto tg_xzzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 99); 

                auto tg_yyyy_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 100); 

                auto tg_yyyy_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 101); 

                auto tg_yyyy_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 102); 

                auto tg_yyyy_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 103); 

                auto tg_yyyy_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 104); 

                auto tg_yyyy_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 105); 

                auto tg_yyyy_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 106); 

                auto tg_yyyy_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 107); 

                auto tg_yyyy_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 108); 

                auto tg_yyyy_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 109); 

                auto tg_yyyz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 110); 

                auto tg_yyyz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 111); 

                auto tg_yyyz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 112); 

                auto tg_yyyz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 113); 

                auto tg_yyyz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 114); 

                auto tg_yyyz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 115); 

                auto tg_yyyz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 116); 

                auto tg_yyyz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 117); 

                auto tg_yyyz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 118); 

                auto tg_yyyz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 119); 

                auto tg_yyzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 120); 

                auto tg_yyzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 121); 

                auto tg_yyzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 122); 

                auto tg_yyzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 123); 

                auto tg_yyzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 124); 

                auto tg_yyzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 125); 

                auto tg_yyzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 126); 

                auto tg_yyzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 127); 

                auto tg_yyzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 128); 

                auto tg_yyzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 129); 

                auto tg_yzzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 130); 

                auto tg_yzzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 131); 

                auto tg_yzzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 132); 

                auto tg_yzzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 133); 

                auto tg_yzzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 134); 

                auto tg_yzzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 135); 

                auto tg_yzzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 136); 

                auto tg_yzzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 137); 

                auto tg_yzzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 138); 

                auto tg_yzzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 139); 

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

                auto tg_yyz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 70); 

                auto tg_yyz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 71); 

                auto tg_yyz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 72); 

                auto tg_yyz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 73); 

                auto tg_yyz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 74); 

                auto tg_yyz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 75); 

                auto tg_yyz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 76); 

                auto tg_yyz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 77); 

                auto tg_yyz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 78); 

                auto tg_yyz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 79); 

                auto tg_yzz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 80); 

                auto tg_yzz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 81); 

                auto tg_yzz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 82); 

                auto tg_yzz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 83); 

                auto tg_yzz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 84); 

                auto tg_yzz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 85); 

                auto tg_yzz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 86); 

                auto tg_yzz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 87); 

                auto tg_yzz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 88); 

                auto tg_yzz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 89); 

                auto tg_zzz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 90); 

                auto tg_zzz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 91); 

                auto tg_zzz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 92); 

                auto tg_zzz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 93); 

                auto tg_zzz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 94); 

                auto tg_zzz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 95); 

                auto tg_zzz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 96); 

                auto tg_zzz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 97); 

                auto tg_zzz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 98); 

                auto tg_zzz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 99); 

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

                auto tg_xyyz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 42); 

                auto tg_xyyz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 43); 

                auto tg_xyyz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 44); 

                auto tg_xyyz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 45); 

                auto tg_xyyz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 46); 

                auto tg_xyyz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 47); 

                auto tg_xyzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 48); 

                auto tg_xyzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 49); 

                auto tg_xyzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 50); 

                auto tg_xyzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 51); 

                auto tg_xyzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 52); 

                auto tg_xyzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 53); 

                auto tg_xzzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 54); 

                auto tg_xzzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 55); 

                auto tg_xzzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 56); 

                auto tg_xzzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 57); 

                auto tg_xzzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 58); 

                auto tg_xzzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 59); 

                auto tg_yyyy_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 60); 

                auto tg_yyyy_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 61); 

                auto tg_yyyy_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 62); 

                auto tg_yyyy_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 63); 

                auto tg_yyyy_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 64); 

                auto tg_yyyy_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 65); 

                auto tg_yyyz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 66); 

                auto tg_yyyz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 67); 

                auto tg_yyyz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 68); 

                auto tg_yyyz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 69); 

                auto tg_yyyz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 70); 

                auto tg_yyyz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 71); 

                auto tg_yyzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 72); 

                auto tg_yyzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 73); 

                auto tg_yyzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 74); 

                auto tg_yyzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 75); 

                auto tg_yyzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 76); 

                auto tg_yyzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 77); 

                auto tg_yzzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 78); 

                auto tg_yzzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 79); 

                auto tg_yzzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 80); 

                auto tg_yzzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 81); 

                auto tg_yzzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 82); 

                auto tg_yzzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 83); 

                // set up pointers to integrals

                auto tg_xxyyz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 70); 

                auto tg_xxyyz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 71); 

                auto tg_xxyyz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 72); 

                auto tg_xxyyz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 73); 

                auto tg_xxyyz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 74); 

                auto tg_xxyyz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 75); 

                auto tg_xxyyz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 76); 

                auto tg_xxyyz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 77); 

                auto tg_xxyyz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 78); 

                auto tg_xxyyz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 79); 

                auto tg_xxyzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 80); 

                auto tg_xxyzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 81); 

                auto tg_xxyzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 82); 

                auto tg_xxyzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 83); 

                auto tg_xxyzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 84); 

                auto tg_xxyzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 85); 

                auto tg_xxyzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 86); 

                auto tg_xxyzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 87); 

                auto tg_xxyzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 88); 

                auto tg_xxyzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 89); 

                auto tg_xxzzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 90); 

                auto tg_xxzzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 91); 

                auto tg_xxzzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 92); 

                auto tg_xxzzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 93); 

                auto tg_xxzzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 94); 

                auto tg_xxzzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 95); 

                auto tg_xxzzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 96); 

                auto tg_xxzzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 97); 

                auto tg_xxzzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 98); 

                auto tg_xxzzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 99); 

                auto tg_xyyyy_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 100); 

                auto tg_xyyyy_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 101); 

                auto tg_xyyyy_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 102); 

                auto tg_xyyyy_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 103); 

                auto tg_xyyyy_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 104); 

                auto tg_xyyyy_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 105); 

                auto tg_xyyyy_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 106); 

                auto tg_xyyyy_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 107); 

                auto tg_xyyyy_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 108); 

                auto tg_xyyyy_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 109); 

                auto tg_xyyyz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 110); 

                auto tg_xyyyz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 111); 

                auto tg_xyyyz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 112); 

                auto tg_xyyyz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 113); 

                auto tg_xyyyz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 114); 

                auto tg_xyyyz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 115); 

                auto tg_xyyyz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 116); 

                auto tg_xyyyz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 117); 

                auto tg_xyyyz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 118); 

                auto tg_xyyyz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 119); 

                auto tg_xyyzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 120); 

                auto tg_xyyzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 121); 

                auto tg_xyyzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 122); 

                auto tg_xyyzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 123); 

                auto tg_xyyzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 124); 

                auto tg_xyyzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 125); 

                auto tg_xyyzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 126); 

                auto tg_xyyzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 127); 

                auto tg_xyyzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 128); 

                auto tg_xyyzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 129); 

                auto tg_xyzzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 130); 

                auto tg_xyzzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 131); 

                auto tg_xyzzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 132); 

                auto tg_xyzzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 133); 

                auto tg_xyzzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 134); 

                auto tg_xyzzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 135); 

                auto tg_xyzzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 136); 

                auto tg_xyzzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 137); 

                auto tg_xyzzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 138); 

                auto tg_xyzzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 139); 

                // Batch of Integrals (70,140)

                #pragma omp simd aligned(fxn, fza, tg_xxyyz_xxx_0, tg_xxyyz_xxy_0, tg_xxyyz_xxz_0, \
                                         tg_xxyyz_xyy_0, tg_xxyyz_xyz_0, tg_xxyyz_xzz_0, tg_xxyyz_yyy_0, tg_xxyyz_yyz_0, \
                                         tg_xxyyz_yzz_0, tg_xxyyz_zzz_0, tg_xxyzz_xxx_0, tg_xxyzz_xxy_0, tg_xxyzz_xxz_0, \
                                         tg_xxyzz_xyy_0, tg_xxyzz_xyz_0, tg_xxyzz_xzz_0, tg_xxyzz_yyy_0, tg_xxyzz_yyz_0, \
                                         tg_xxyzz_yzz_0, tg_xxyzz_zzz_0, tg_xxzzz_xxx_0, tg_xxzzz_xxy_0, tg_xxzzz_xxz_0, \
                                         tg_xxzzz_xyy_0, tg_xxzzz_xyz_0, tg_xxzzz_xzz_0, tg_xxzzz_yyy_0, tg_xxzzz_yyz_0, \
                                         tg_xxzzz_yzz_0, tg_xxzzz_zzz_0, tg_xyyyy_xxx_0, tg_xyyyy_xxy_0, tg_xyyyy_xxz_0, \
                                         tg_xyyyy_xyy_0, tg_xyyyy_xyz_0, tg_xyyyy_xzz_0, tg_xyyyy_yyy_0, tg_xyyyy_yyz_0, \
                                         tg_xyyyy_yzz_0, tg_xyyyy_zzz_0, tg_xyyyz_xxx_0, tg_xyyyz_xxy_0, tg_xyyyz_xxz_0, \
                                         tg_xyyyz_xyy_0, tg_xyyyz_xyz_0, tg_xyyyz_xzz_0, tg_xyyyz_yyy_0, tg_xyyyz_yyz_0, \
                                         tg_xyyyz_yzz_0, tg_xyyyz_zzz_0, tg_xyyz_xx_1, tg_xyyz_xxx_0, tg_xyyz_xxx_1, \
                                         tg_xyyz_xxy_0, tg_xyyz_xxy_1, tg_xyyz_xxz_0, tg_xyyz_xxz_1, tg_xyyz_xy_1, \
                                         tg_xyyz_xyy_0, tg_xyyz_xyy_1, tg_xyyz_xyz_0, tg_xyyz_xyz_1, tg_xyyz_xz_1, \
                                         tg_xyyz_xzz_0, tg_xyyz_xzz_1, tg_xyyz_yy_1, tg_xyyz_yyy_0, tg_xyyz_yyy_1, \
                                         tg_xyyz_yyz_0, tg_xyyz_yyz_1, tg_xyyz_yz_1, tg_xyyz_yzz_0, tg_xyyz_yzz_1, \
                                         tg_xyyz_zz_1, tg_xyyz_zzz_0, tg_xyyz_zzz_1, tg_xyyzz_xxx_0, tg_xyyzz_xxy_0, \
                                         tg_xyyzz_xxz_0, tg_xyyzz_xyy_0, tg_xyyzz_xyz_0, tg_xyyzz_xzz_0, tg_xyyzz_yyy_0, \
                                         tg_xyyzz_yyz_0, tg_xyyzz_yzz_0, tg_xyyzz_zzz_0, tg_xyzz_xx_1, tg_xyzz_xxx_0, \
                                         tg_xyzz_xxx_1, tg_xyzz_xxy_0, tg_xyzz_xxy_1, tg_xyzz_xxz_0, tg_xyzz_xxz_1, \
                                         tg_xyzz_xy_1, tg_xyzz_xyy_0, tg_xyzz_xyy_1, tg_xyzz_xyz_0, tg_xyzz_xyz_1, \
                                         tg_xyzz_xz_1, tg_xyzz_xzz_0, tg_xyzz_xzz_1, tg_xyzz_yy_1, tg_xyzz_yyy_0, \
                                         tg_xyzz_yyy_1, tg_xyzz_yyz_0, tg_xyzz_yyz_1, tg_xyzz_yz_1, tg_xyzz_yzz_0, \
                                         tg_xyzz_yzz_1, tg_xyzz_zz_1, tg_xyzz_zzz_0, tg_xyzz_zzz_1, tg_xyzzz_xxx_0, \
                                         tg_xyzzz_xxy_0, tg_xyzzz_xxz_0, tg_xyzzz_xyy_0, tg_xyzzz_xyz_0, tg_xyzzz_xzz_0, \
                                         tg_xyzzz_yyy_0, tg_xyzzz_yyz_0, tg_xyzzz_yzz_0, tg_xyzzz_zzz_0, tg_xzzz_xx_1, \
                                         tg_xzzz_xxx_0, tg_xzzz_xxx_1, tg_xzzz_xxy_0, tg_xzzz_xxy_1, tg_xzzz_xxz_0, \
                                         tg_xzzz_xxz_1, tg_xzzz_xy_1, tg_xzzz_xyy_0, tg_xzzz_xyy_1, tg_xzzz_xyz_0, \
                                         tg_xzzz_xyz_1, tg_xzzz_xz_1, tg_xzzz_xzz_0, tg_xzzz_xzz_1, tg_xzzz_yy_1, \
                                         tg_xzzz_yyy_0, tg_xzzz_yyy_1, tg_xzzz_yyz_0, tg_xzzz_yyz_1, tg_xzzz_yz_1, \
                                         tg_xzzz_yzz_0, tg_xzzz_yzz_1, tg_xzzz_zz_1, tg_xzzz_zzz_0, tg_xzzz_zzz_1, \
                                         tg_yyyy_xx_1, tg_yyyy_xxx_0, tg_yyyy_xxx_1, tg_yyyy_xxy_0, tg_yyyy_xxy_1, \
                                         tg_yyyy_xxz_0, tg_yyyy_xxz_1, tg_yyyy_xy_1, tg_yyyy_xyy_0, tg_yyyy_xyy_1, \
                                         tg_yyyy_xyz_0, tg_yyyy_xyz_1, tg_yyyy_xz_1, tg_yyyy_xzz_0, tg_yyyy_xzz_1, \
                                         tg_yyyy_yy_1, tg_yyyy_yyy_0, tg_yyyy_yyy_1, tg_yyyy_yyz_0, tg_yyyy_yyz_1, \
                                         tg_yyyy_yz_1, tg_yyyy_yzz_0, tg_yyyy_yzz_1, tg_yyyy_zz_1, tg_yyyy_zzz_0, \
                                         tg_yyyy_zzz_1, tg_yyyz_xx_1, tg_yyyz_xxx_0, tg_yyyz_xxx_1, tg_yyyz_xxy_0, \
                                         tg_yyyz_xxy_1, tg_yyyz_xxz_0, tg_yyyz_xxz_1, tg_yyyz_xy_1, tg_yyyz_xyy_0, \
                                         tg_yyyz_xyy_1, tg_yyyz_xyz_0, tg_yyyz_xyz_1, tg_yyyz_xz_1, tg_yyyz_xzz_0, \
                                         tg_yyyz_xzz_1, tg_yyyz_yy_1, tg_yyyz_yyy_0, tg_yyyz_yyy_1, tg_yyyz_yyz_0, \
                                         tg_yyyz_yyz_1, tg_yyyz_yz_1, tg_yyyz_yzz_0, tg_yyyz_yzz_1, tg_yyyz_zz_1, \
                                         tg_yyyz_zzz_0, tg_yyyz_zzz_1, tg_yyz_xxx_0, tg_yyz_xxx_1, tg_yyz_xxy_0, \
                                         tg_yyz_xxy_1, tg_yyz_xxz_0, tg_yyz_xxz_1, tg_yyz_xyy_0, tg_yyz_xyy_1, tg_yyz_xyz_0, \
                                         tg_yyz_xyz_1, tg_yyz_xzz_0, tg_yyz_xzz_1, tg_yyz_yyy_0, tg_yyz_yyy_1, tg_yyz_yyz_0, \
                                         tg_yyz_yyz_1, tg_yyz_yzz_0, tg_yyz_yzz_1, tg_yyz_zzz_0, tg_yyz_zzz_1, tg_yyzz_xx_1, \
                                         tg_yyzz_xxx_0, tg_yyzz_xxx_1, tg_yyzz_xxy_0, tg_yyzz_xxy_1, tg_yyzz_xxz_0, \
                                         tg_yyzz_xxz_1, tg_yyzz_xy_1, tg_yyzz_xyy_0, tg_yyzz_xyy_1, tg_yyzz_xyz_0, \
                                         tg_yyzz_xyz_1, tg_yyzz_xz_1, tg_yyzz_xzz_0, tg_yyzz_xzz_1, tg_yyzz_yy_1, \
                                         tg_yyzz_yyy_0, tg_yyzz_yyy_1, tg_yyzz_yyz_0, tg_yyzz_yyz_1, tg_yyzz_yz_1, \
                                         tg_yyzz_yzz_0, tg_yyzz_yzz_1, tg_yyzz_zz_1, tg_yyzz_zzz_0, tg_yyzz_zzz_1, \
                                         tg_yzz_xxx_0, tg_yzz_xxx_1, tg_yzz_xxy_0, tg_yzz_xxy_1, tg_yzz_xxz_0, tg_yzz_xxz_1, \
                                         tg_yzz_xyy_0, tg_yzz_xyy_1, tg_yzz_xyz_0, tg_yzz_xyz_1, tg_yzz_xzz_0, tg_yzz_xzz_1, \
                                         tg_yzz_yyy_0, tg_yzz_yyy_1, tg_yzz_yyz_0, tg_yzz_yyz_1, tg_yzz_yzz_0, tg_yzz_yzz_1, \
                                         tg_yzz_zzz_0, tg_yzz_zzz_1, tg_yzzz_xx_1, tg_yzzz_xxx_0, tg_yzzz_xxx_1, \
                                         tg_yzzz_xxy_0, tg_yzzz_xxy_1, tg_yzzz_xxz_0, tg_yzzz_xxz_1, tg_yzzz_xy_1, \
                                         tg_yzzz_xyy_0, tg_yzzz_xyy_1, tg_yzzz_xyz_0, tg_yzzz_xyz_1, tg_yzzz_xz_1, \
                                         tg_yzzz_xzz_0, tg_yzzz_xzz_1, tg_yzzz_yy_1, tg_yzzz_yyy_0, tg_yzzz_yyy_1, \
                                         tg_yzzz_yyz_0, tg_yzzz_yyz_1, tg_yzzz_yz_1, tg_yzzz_yzz_0, tg_yzzz_yzz_1, \
                                         tg_yzzz_zz_1, tg_yzzz_zzz_0, tg_yzzz_zzz_1, tg_zzz_xxx_0, tg_zzz_xxx_1, \
                                         tg_zzz_xxy_0, tg_zzz_xxy_1, tg_zzz_xxz_0, tg_zzz_xxz_1, tg_zzz_xyy_0, tg_zzz_xyy_1, \
                                         tg_zzz_xyz_0, tg_zzz_xyz_1, tg_zzz_xzz_0, tg_zzz_xzz_1, tg_zzz_yyy_0, tg_zzz_yyy_1, \
                                         tg_zzz_yyz_0, tg_zzz_yyz_1, tg_zzz_yzz_0, tg_zzz_yzz_1, tg_zzz_zzz_0, tg_zzz_zzz_1, \
                                         wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyyz_xxx_0[j] = pb_x * tg_xyyz_xxx_0[j] + wp_x[j] * tg_xyyz_xxx_1[j] + 0.5 * fl1_fx * tg_yyz_xxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxx_1[j] + 1.5 * fl1_fxn * tg_xyyz_xx_1[j];

                    tg_xxyyz_xxy_0[j] = pb_x * tg_xyyz_xxy_0[j] + wp_x[j] * tg_xyyz_xxy_1[j] + 0.5 * fl1_fx * tg_yyz_xxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxy_1[j] + fl1_fxn * tg_xyyz_xy_1[j];

                    tg_xxyyz_xxz_0[j] = pb_x * tg_xyyz_xxz_0[j] + wp_x[j] * tg_xyyz_xxz_1[j] + 0.5 * fl1_fx * tg_yyz_xxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xxz_1[j] + fl1_fxn * tg_xyyz_xz_1[j];

                    tg_xxyyz_xyy_0[j] = pb_x * tg_xyyz_xyy_0[j] + wp_x[j] * tg_xyyz_xyy_1[j] + 0.5 * fl1_fx * tg_yyz_xyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyy_1[j] + 0.5 * fl1_fxn * tg_xyyz_yy_1[j];

                    tg_xxyyz_xyz_0[j] = pb_x * tg_xyyz_xyz_0[j] + wp_x[j] * tg_xyyz_xyz_1[j] + 0.5 * fl1_fx * tg_yyz_xyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xyz_1[j] + 0.5 * fl1_fxn * tg_xyyz_yz_1[j];

                    tg_xxyyz_xzz_0[j] = pb_x * tg_xyyz_xzz_0[j] + wp_x[j] * tg_xyyz_xzz_1[j] + 0.5 * fl1_fx * tg_yyz_xzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_xzz_1[j] + 0.5 * fl1_fxn * tg_xyyz_zz_1[j];

                    tg_xxyyz_yyy_0[j] = pb_x * tg_xyyz_yyy_0[j] + wp_x[j] * tg_xyyz_yyy_1[j] + 0.5 * fl1_fx * tg_yyz_yyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyy_1[j];

                    tg_xxyyz_yyz_0[j] = pb_x * tg_xyyz_yyz_0[j] + wp_x[j] * tg_xyyz_yyz_1[j] + 0.5 * fl1_fx * tg_yyz_yyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yyz_1[j];

                    tg_xxyyz_yzz_0[j] = pb_x * tg_xyyz_yzz_0[j] + wp_x[j] * tg_xyyz_yzz_1[j] + 0.5 * fl1_fx * tg_yyz_yzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_yzz_1[j];

                    tg_xxyyz_zzz_0[j] = pb_x * tg_xyyz_zzz_0[j] + wp_x[j] * tg_xyyz_zzz_1[j] + 0.5 * fl1_fx * tg_yyz_zzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_zzz_1[j];

                    tg_xxyzz_xxx_0[j] = pb_x * tg_xyzz_xxx_0[j] + wp_x[j] * tg_xyzz_xxx_1[j] + 0.5 * fl1_fx * tg_yzz_xxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxx_1[j] + 1.5 * fl1_fxn * tg_xyzz_xx_1[j];

                    tg_xxyzz_xxy_0[j] = pb_x * tg_xyzz_xxy_0[j] + wp_x[j] * tg_xyzz_xxy_1[j] + 0.5 * fl1_fx * tg_yzz_xxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxy_1[j] + fl1_fxn * tg_xyzz_xy_1[j];

                    tg_xxyzz_xxz_0[j] = pb_x * tg_xyzz_xxz_0[j] + wp_x[j] * tg_xyzz_xxz_1[j] + 0.5 * fl1_fx * tg_yzz_xxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xxz_1[j] + fl1_fxn * tg_xyzz_xz_1[j];

                    tg_xxyzz_xyy_0[j] = pb_x * tg_xyzz_xyy_0[j] + wp_x[j] * tg_xyzz_xyy_1[j] + 0.5 * fl1_fx * tg_yzz_xyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyy_1[j] + 0.5 * fl1_fxn * tg_xyzz_yy_1[j];

                    tg_xxyzz_xyz_0[j] = pb_x * tg_xyzz_xyz_0[j] + wp_x[j] * tg_xyzz_xyz_1[j] + 0.5 * fl1_fx * tg_yzz_xyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xyz_1[j] + 0.5 * fl1_fxn * tg_xyzz_yz_1[j];

                    tg_xxyzz_xzz_0[j] = pb_x * tg_xyzz_xzz_0[j] + wp_x[j] * tg_xyzz_xzz_1[j] + 0.5 * fl1_fx * tg_yzz_xzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_xzz_1[j] + 0.5 * fl1_fxn * tg_xyzz_zz_1[j];

                    tg_xxyzz_yyy_0[j] = pb_x * tg_xyzz_yyy_0[j] + wp_x[j] * tg_xyzz_yyy_1[j] + 0.5 * fl1_fx * tg_yzz_yyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyy_1[j];

                    tg_xxyzz_yyz_0[j] = pb_x * tg_xyzz_yyz_0[j] + wp_x[j] * tg_xyzz_yyz_1[j] + 0.5 * fl1_fx * tg_yzz_yyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yyz_1[j];

                    tg_xxyzz_yzz_0[j] = pb_x * tg_xyzz_yzz_0[j] + wp_x[j] * tg_xyzz_yzz_1[j] + 0.5 * fl1_fx * tg_yzz_yzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_yzz_1[j];

                    tg_xxyzz_zzz_0[j] = pb_x * tg_xyzz_zzz_0[j] + wp_x[j] * tg_xyzz_zzz_1[j] + 0.5 * fl1_fx * tg_yzz_zzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_zzz_1[j];

                    tg_xxzzz_xxx_0[j] = pb_x * tg_xzzz_xxx_0[j] + wp_x[j] * tg_xzzz_xxx_1[j] + 0.5 * fl1_fx * tg_zzz_xxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxx_1[j] + 1.5 * fl1_fxn * tg_xzzz_xx_1[j];

                    tg_xxzzz_xxy_0[j] = pb_x * tg_xzzz_xxy_0[j] + wp_x[j] * tg_xzzz_xxy_1[j] + 0.5 * fl1_fx * tg_zzz_xxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxy_1[j] + fl1_fxn * tg_xzzz_xy_1[j];

                    tg_xxzzz_xxz_0[j] = pb_x * tg_xzzz_xxz_0[j] + wp_x[j] * tg_xzzz_xxz_1[j] + 0.5 * fl1_fx * tg_zzz_xxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxz_1[j] + fl1_fxn * tg_xzzz_xz_1[j];

                    tg_xxzzz_xyy_0[j] = pb_x * tg_xzzz_xyy_0[j] + wp_x[j] * tg_xzzz_xyy_1[j] + 0.5 * fl1_fx * tg_zzz_xyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyy_1[j] + 0.5 * fl1_fxn * tg_xzzz_yy_1[j];

                    tg_xxzzz_xyz_0[j] = pb_x * tg_xzzz_xyz_0[j] + wp_x[j] * tg_xzzz_xyz_1[j] + 0.5 * fl1_fx * tg_zzz_xyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyz_1[j] + 0.5 * fl1_fxn * tg_xzzz_yz_1[j];

                    tg_xxzzz_xzz_0[j] = pb_x * tg_xzzz_xzz_0[j] + wp_x[j] * tg_xzzz_xzz_1[j] + 0.5 * fl1_fx * tg_zzz_xzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xzz_1[j] + 0.5 * fl1_fxn * tg_xzzz_zz_1[j];

                    tg_xxzzz_yyy_0[j] = pb_x * tg_xzzz_yyy_0[j] + wp_x[j] * tg_xzzz_yyy_1[j] + 0.5 * fl1_fx * tg_zzz_yyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyy_1[j];

                    tg_xxzzz_yyz_0[j] = pb_x * tg_xzzz_yyz_0[j] + wp_x[j] * tg_xzzz_yyz_1[j] + 0.5 * fl1_fx * tg_zzz_yyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyz_1[j];

                    tg_xxzzz_yzz_0[j] = pb_x * tg_xzzz_yzz_0[j] + wp_x[j] * tg_xzzz_yzz_1[j] + 0.5 * fl1_fx * tg_zzz_yzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yzz_1[j];

                    tg_xxzzz_zzz_0[j] = pb_x * tg_xzzz_zzz_0[j] + wp_x[j] * tg_xzzz_zzz_1[j] + 0.5 * fl1_fx * tg_zzz_zzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_zzz_1[j];

                    tg_xyyyy_xxx_0[j] = pb_x * tg_yyyy_xxx_0[j] + wp_x[j] * tg_yyyy_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyy_xx_1[j];

                    tg_xyyyy_xxy_0[j] = pb_x * tg_yyyy_xxy_0[j] + wp_x[j] * tg_yyyy_xxy_1[j] + fl1_fxn * tg_yyyy_xy_1[j];

                    tg_xyyyy_xxz_0[j] = pb_x * tg_yyyy_xxz_0[j] + wp_x[j] * tg_yyyy_xxz_1[j] + fl1_fxn * tg_yyyy_xz_1[j];

                    tg_xyyyy_xyy_0[j] = pb_x * tg_yyyy_xyy_0[j] + wp_x[j] * tg_yyyy_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyy_yy_1[j];

                    tg_xyyyy_xyz_0[j] = pb_x * tg_yyyy_xyz_0[j] + wp_x[j] * tg_yyyy_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_yz_1[j];

                    tg_xyyyy_xzz_0[j] = pb_x * tg_yyyy_xzz_0[j] + wp_x[j] * tg_yyyy_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zz_1[j];

                    tg_xyyyy_yyy_0[j] = pb_x * tg_yyyy_yyy_0[j] + wp_x[j] * tg_yyyy_yyy_1[j];

                    tg_xyyyy_yyz_0[j] = pb_x * tg_yyyy_yyz_0[j] + wp_x[j] * tg_yyyy_yyz_1[j];

                    tg_xyyyy_yzz_0[j] = pb_x * tg_yyyy_yzz_0[j] + wp_x[j] * tg_yyyy_yzz_1[j];

                    tg_xyyyy_zzz_0[j] = pb_x * tg_yyyy_zzz_0[j] + wp_x[j] * tg_yyyy_zzz_1[j];

                    tg_xyyyz_xxx_0[j] = pb_x * tg_yyyz_xxx_0[j] + wp_x[j] * tg_yyyz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyz_xx_1[j];

                    tg_xyyyz_xxy_0[j] = pb_x * tg_yyyz_xxy_0[j] + wp_x[j] * tg_yyyz_xxy_1[j] + fl1_fxn * tg_yyyz_xy_1[j];

                    tg_xyyyz_xxz_0[j] = pb_x * tg_yyyz_xxz_0[j] + wp_x[j] * tg_yyyz_xxz_1[j] + fl1_fxn * tg_yyyz_xz_1[j];

                    tg_xyyyz_xyy_0[j] = pb_x * tg_yyyz_xyy_0[j] + wp_x[j] * tg_yyyz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyz_yy_1[j];

                    tg_xyyyz_xyz_0[j] = pb_x * tg_yyyz_xyz_0[j] + wp_x[j] * tg_yyyz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_yz_1[j];

                    tg_xyyyz_xzz_0[j] = pb_x * tg_yyyz_xzz_0[j] + wp_x[j] * tg_yyyz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zz_1[j];

                    tg_xyyyz_yyy_0[j] = pb_x * tg_yyyz_yyy_0[j] + wp_x[j] * tg_yyyz_yyy_1[j];

                    tg_xyyyz_yyz_0[j] = pb_x * tg_yyyz_yyz_0[j] + wp_x[j] * tg_yyyz_yyz_1[j];

                    tg_xyyyz_yzz_0[j] = pb_x * tg_yyyz_yzz_0[j] + wp_x[j] * tg_yyyz_yzz_1[j];

                    tg_xyyyz_zzz_0[j] = pb_x * tg_yyyz_zzz_0[j] + wp_x[j] * tg_yyyz_zzz_1[j];

                    tg_xyyzz_xxx_0[j] = pb_x * tg_yyzz_xxx_0[j] + wp_x[j] * tg_yyzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyzz_xx_1[j];

                    tg_xyyzz_xxy_0[j] = pb_x * tg_yyzz_xxy_0[j] + wp_x[j] * tg_yyzz_xxy_1[j] + fl1_fxn * tg_yyzz_xy_1[j];

                    tg_xyyzz_xxz_0[j] = pb_x * tg_yyzz_xxz_0[j] + wp_x[j] * tg_yyzz_xxz_1[j] + fl1_fxn * tg_yyzz_xz_1[j];

                    tg_xyyzz_xyy_0[j] = pb_x * tg_yyzz_xyy_0[j] + wp_x[j] * tg_yyzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyzz_yy_1[j];

                    tg_xyyzz_xyz_0[j] = pb_x * tg_yyzz_xyz_0[j] + wp_x[j] * tg_yyzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_yz_1[j];

                    tg_xyyzz_xzz_0[j] = pb_x * tg_yyzz_xzz_0[j] + wp_x[j] * tg_yyzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zz_1[j];

                    tg_xyyzz_yyy_0[j] = pb_x * tg_yyzz_yyy_0[j] + wp_x[j] * tg_yyzz_yyy_1[j];

                    tg_xyyzz_yyz_0[j] = pb_x * tg_yyzz_yyz_0[j] + wp_x[j] * tg_yyzz_yyz_1[j];

                    tg_xyyzz_yzz_0[j] = pb_x * tg_yyzz_yzz_0[j] + wp_x[j] * tg_yyzz_yzz_1[j];

                    tg_xyyzz_zzz_0[j] = pb_x * tg_yyzz_zzz_0[j] + wp_x[j] * tg_yyzz_zzz_1[j];

                    tg_xyzzz_xxx_0[j] = pb_x * tg_yzzz_xxx_0[j] + wp_x[j] * tg_yzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yzzz_xx_1[j];

                    tg_xyzzz_xxy_0[j] = pb_x * tg_yzzz_xxy_0[j] + wp_x[j] * tg_yzzz_xxy_1[j] + fl1_fxn * tg_yzzz_xy_1[j];

                    tg_xyzzz_xxz_0[j] = pb_x * tg_yzzz_xxz_0[j] + wp_x[j] * tg_yzzz_xxz_1[j] + fl1_fxn * tg_yzzz_xz_1[j];

                    tg_xyzzz_xyy_0[j] = pb_x * tg_yzzz_xyy_0[j] + wp_x[j] * tg_yzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yzzz_yy_1[j];

                    tg_xyzzz_xyz_0[j] = pb_x * tg_yzzz_xyz_0[j] + wp_x[j] * tg_yzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_yz_1[j];

                    tg_xyzzz_xzz_0[j] = pb_x * tg_yzzz_xzz_0[j] + wp_x[j] * tg_yzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zz_1[j];

                    tg_xyzzz_yyy_0[j] = pb_x * tg_yzzz_yyy_0[j] + wp_x[j] * tg_yzzz_yyy_1[j];

                    tg_xyzzz_yyz_0[j] = pb_x * tg_yzzz_yyz_0[j] + wp_x[j] * tg_yzzz_yyz_1[j];

                    tg_xyzzz_yzz_0[j] = pb_x * tg_yzzz_yzz_0[j] + wp_x[j] * tg_yzzz_yzz_1[j];

                    tg_xyzzz_zzz_0[j] = pb_x * tg_yzzz_zzz_0[j] + wp_x[j] * tg_yzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSF_140_210(      CMemBlock2D<double>& primBuffer,
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

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {5, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

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

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_yyyy_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 100); 

                auto tg_yyyy_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 101); 

                auto tg_yyyy_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 102); 

                auto tg_yyyy_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 103); 

                auto tg_yyyy_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 104); 

                auto tg_yyyy_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 105); 

                auto tg_yyyy_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 106); 

                auto tg_yyyy_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 107); 

                auto tg_yyyy_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 108); 

                auto tg_yyyy_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 109); 

                auto tg_yyyz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 110); 

                auto tg_yyyz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 111); 

                auto tg_yyyz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 112); 

                auto tg_yyyz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 113); 

                auto tg_yyyz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 114); 

                auto tg_yyyz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 115); 

                auto tg_yyyz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 116); 

                auto tg_yyyz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 117); 

                auto tg_yyyz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 118); 

                auto tg_yyyz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 119); 

                auto tg_yyzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 120); 

                auto tg_yyzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 121); 

                auto tg_yyzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 122); 

                auto tg_yyzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 123); 

                auto tg_yyzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 124); 

                auto tg_yyzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 125); 

                auto tg_yyzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 126); 

                auto tg_yyzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 127); 

                auto tg_yyzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 128); 

                auto tg_yyzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 129); 

                auto tg_yzzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 130); 

                auto tg_yzzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 131); 

                auto tg_yzzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 132); 

                auto tg_yzzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 133); 

                auto tg_yzzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 134); 

                auto tg_yzzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 135); 

                auto tg_yzzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 136); 

                auto tg_yzzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 137); 

                auto tg_yzzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 138); 

                auto tg_yzzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 139); 

                auto tg_zzzz_xxx_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 140); 

                auto tg_zzzz_xxy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 141); 

                auto tg_zzzz_xxz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 142); 

                auto tg_zzzz_xyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 143); 

                auto tg_zzzz_xyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 144); 

                auto tg_zzzz_xzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 145); 

                auto tg_zzzz_yyy_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 146); 

                auto tg_zzzz_yyz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 147); 

                auto tg_zzzz_yzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 148); 

                auto tg_zzzz_zzz_0 = primBuffer.data(pidx_g_4_3_m0 + 150 * idx + 149); 

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

                auto tg_yyy_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 60); 

                auto tg_yyy_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 61); 

                auto tg_yyy_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 62); 

                auto tg_yyy_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 63); 

                auto tg_yyy_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 64); 

                auto tg_yyy_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 65); 

                auto tg_yyy_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 66); 

                auto tg_yyy_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 67); 

                auto tg_yyy_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 68); 

                auto tg_yyy_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 69); 

                auto tg_yyz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 70); 

                auto tg_yyz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 71); 

                auto tg_yyz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 72); 

                auto tg_yyz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 73); 

                auto tg_yyz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 74); 

                auto tg_yyz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 75); 

                auto tg_yyz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 76); 

                auto tg_yyz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 77); 

                auto tg_yyz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 78); 

                auto tg_yyz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 79); 

                auto tg_yzz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 80); 

                auto tg_yzz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 81); 

                auto tg_yzz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 82); 

                auto tg_yzz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 83); 

                auto tg_yzz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 84); 

                auto tg_yzz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 85); 

                auto tg_yzz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 86); 

                auto tg_yzz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 87); 

                auto tg_yzz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 88); 

                auto tg_yzz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 89); 

                auto tg_zzz_xxx_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 90); 

                auto tg_zzz_xxy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 91); 

                auto tg_zzz_xxz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 92); 

                auto tg_zzz_xyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 93); 

                auto tg_zzz_xyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 94); 

                auto tg_zzz_xzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 95); 

                auto tg_zzz_yyy_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 96); 

                auto tg_zzz_yyz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 97); 

                auto tg_zzz_yzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 98); 

                auto tg_zzz_zzz_0 = primBuffer.data(pidx_g_3_3_m0 + 100 * idx + 99); 

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

                auto tg_yyyy_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 60); 

                auto tg_yyyy_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 61); 

                auto tg_yyyy_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 62); 

                auto tg_yyyy_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 63); 

                auto tg_yyyy_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 64); 

                auto tg_yyyy_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 65); 

                auto tg_yyyz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 66); 

                auto tg_yyyz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 67); 

                auto tg_yyyz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 68); 

                auto tg_yyyz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 69); 

                auto tg_yyyz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 70); 

                auto tg_yyyz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 71); 

                auto tg_yyzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 72); 

                auto tg_yyzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 73); 

                auto tg_yyzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 74); 

                auto tg_yyzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 75); 

                auto tg_yyzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 76); 

                auto tg_yyzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 77); 

                auto tg_yzzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 78); 

                auto tg_yzzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 79); 

                auto tg_yzzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 80); 

                auto tg_yzzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 81); 

                auto tg_yzzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 82); 

                auto tg_yzzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 83); 

                auto tg_zzzz_xx_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 84); 

                auto tg_zzzz_xy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 85); 

                auto tg_zzzz_xz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 86); 

                auto tg_zzzz_yy_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 87); 

                auto tg_zzzz_yz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 88); 

                auto tg_zzzz_zz_1 = primBuffer.data(pidx_g_4_2_m1 + 90 * idx + 89); 

                // set up pointers to integrals

                auto tg_xzzzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 140); 

                auto tg_xzzzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 141); 

                auto tg_xzzzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 142); 

                auto tg_xzzzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 143); 

                auto tg_xzzzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 144); 

                auto tg_xzzzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 145); 

                auto tg_xzzzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 146); 

                auto tg_xzzzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 147); 

                auto tg_xzzzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 148); 

                auto tg_xzzzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 149); 

                auto tg_yyyyy_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 150); 

                auto tg_yyyyy_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 151); 

                auto tg_yyyyy_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 152); 

                auto tg_yyyyy_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 153); 

                auto tg_yyyyy_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 154); 

                auto tg_yyyyy_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 155); 

                auto tg_yyyyy_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 156); 

                auto tg_yyyyy_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 157); 

                auto tg_yyyyy_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 158); 

                auto tg_yyyyy_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 159); 

                auto tg_yyyyz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 160); 

                auto tg_yyyyz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 161); 

                auto tg_yyyyz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 162); 

                auto tg_yyyyz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 163); 

                auto tg_yyyyz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 164); 

                auto tg_yyyyz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 165); 

                auto tg_yyyyz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 166); 

                auto tg_yyyyz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 167); 

                auto tg_yyyyz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 168); 

                auto tg_yyyyz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 169); 

                auto tg_yyyzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 170); 

                auto tg_yyyzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 171); 

                auto tg_yyyzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 172); 

                auto tg_yyyzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 173); 

                auto tg_yyyzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 174); 

                auto tg_yyyzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 175); 

                auto tg_yyyzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 176); 

                auto tg_yyyzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 177); 

                auto tg_yyyzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 178); 

                auto tg_yyyzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 179); 

                auto tg_yyzzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 180); 

                auto tg_yyzzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 181); 

                auto tg_yyzzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 182); 

                auto tg_yyzzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 183); 

                auto tg_yyzzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 184); 

                auto tg_yyzzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 185); 

                auto tg_yyzzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 186); 

                auto tg_yyzzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 187); 

                auto tg_yyzzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 188); 

                auto tg_yyzzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 189); 

                auto tg_yzzzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 190); 

                auto tg_yzzzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 191); 

                auto tg_yzzzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 192); 

                auto tg_yzzzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 193); 

                auto tg_yzzzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 194); 

                auto tg_yzzzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 195); 

                auto tg_yzzzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 196); 

                auto tg_yzzzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 197); 

                auto tg_yzzzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 198); 

                auto tg_yzzzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 199); 

                auto tg_zzzzz_xxx_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 200); 

                auto tg_zzzzz_xxy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 201); 

                auto tg_zzzzz_xxz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 202); 

                auto tg_zzzzz_xyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 203); 

                auto tg_zzzzz_xyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 204); 

                auto tg_zzzzz_xzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 205); 

                auto tg_zzzzz_yyy_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 206); 

                auto tg_zzzzz_yyz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 207); 

                auto tg_zzzzz_yzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 208); 

                auto tg_zzzzz_zzz_0 = primBuffer.data(pidx_g_5_3_m0 + 210 * idx + 209); 

                // Batch of Integrals (140,210)

                #pragma omp simd aligned(fxn, fza, tg_xzzzz_xxx_0, tg_xzzzz_xxy_0, tg_xzzzz_xxz_0, \
                                         tg_xzzzz_xyy_0, tg_xzzzz_xyz_0, tg_xzzzz_xzz_0, tg_xzzzz_yyy_0, tg_xzzzz_yyz_0, \
                                         tg_xzzzz_yzz_0, tg_xzzzz_zzz_0, tg_yyy_xxx_0, tg_yyy_xxx_1, tg_yyy_xxy_0, \
                                         tg_yyy_xxy_1, tg_yyy_xxz_0, tg_yyy_xxz_1, tg_yyy_xyy_0, tg_yyy_xyy_1, tg_yyy_xyz_0, \
                                         tg_yyy_xyz_1, tg_yyy_xzz_0, tg_yyy_xzz_1, tg_yyy_yyy_0, tg_yyy_yyy_1, tg_yyy_yyz_0, \
                                         tg_yyy_yyz_1, tg_yyy_yzz_0, tg_yyy_yzz_1, tg_yyy_zzz_0, tg_yyy_zzz_1, tg_yyyy_xx_1, \
                                         tg_yyyy_xxx_0, tg_yyyy_xxx_1, tg_yyyy_xxy_0, tg_yyyy_xxy_1, tg_yyyy_xxz_0, \
                                         tg_yyyy_xxz_1, tg_yyyy_xy_1, tg_yyyy_xyy_0, tg_yyyy_xyy_1, tg_yyyy_xyz_0, \
                                         tg_yyyy_xyz_1, tg_yyyy_xz_1, tg_yyyy_xzz_0, tg_yyyy_xzz_1, tg_yyyy_yy_1, \
                                         tg_yyyy_yyy_0, tg_yyyy_yyy_1, tg_yyyy_yyz_0, tg_yyyy_yyz_1, tg_yyyy_yz_1, \
                                         tg_yyyy_yzz_0, tg_yyyy_yzz_1, tg_yyyy_zz_1, tg_yyyy_zzz_0, tg_yyyy_zzz_1, \
                                         tg_yyyyy_xxx_0, tg_yyyyy_xxy_0, tg_yyyyy_xxz_0, tg_yyyyy_xyy_0, tg_yyyyy_xyz_0, \
                                         tg_yyyyy_xzz_0, tg_yyyyy_yyy_0, tg_yyyyy_yyz_0, tg_yyyyy_yzz_0, tg_yyyyy_zzz_0, \
                                         tg_yyyyz_xxx_0, tg_yyyyz_xxy_0, tg_yyyyz_xxz_0, tg_yyyyz_xyy_0, tg_yyyyz_xyz_0, \
                                         tg_yyyyz_xzz_0, tg_yyyyz_yyy_0, tg_yyyyz_yyz_0, tg_yyyyz_yzz_0, tg_yyyyz_zzz_0, \
                                         tg_yyyz_xx_1, tg_yyyz_xxx_0, tg_yyyz_xxx_1, tg_yyyz_xxy_0, tg_yyyz_xxy_1, \
                                         tg_yyyz_xxz_0, tg_yyyz_xxz_1, tg_yyyz_xy_1, tg_yyyz_xyy_0, tg_yyyz_xyy_1, \
                                         tg_yyyz_xyz_0, tg_yyyz_xyz_1, tg_yyyz_xz_1, tg_yyyz_xzz_0, tg_yyyz_xzz_1, \
                                         tg_yyyz_yy_1, tg_yyyz_yyy_0, tg_yyyz_yyy_1, tg_yyyz_yyz_0, tg_yyyz_yyz_1, \
                                         tg_yyyz_yz_1, tg_yyyz_yzz_0, tg_yyyz_yzz_1, tg_yyyz_zz_1, tg_yyyz_zzz_0, \
                                         tg_yyyz_zzz_1, tg_yyyzz_xxx_0, tg_yyyzz_xxy_0, tg_yyyzz_xxz_0, tg_yyyzz_xyy_0, \
                                         tg_yyyzz_xyz_0, tg_yyyzz_xzz_0, tg_yyyzz_yyy_0, tg_yyyzz_yyz_0, tg_yyyzz_yzz_0, \
                                         tg_yyyzz_zzz_0, tg_yyz_xxx_0, tg_yyz_xxx_1, tg_yyz_xxy_0, tg_yyz_xxy_1, tg_yyz_xxz_0, \
                                         tg_yyz_xxz_1, tg_yyz_xyy_0, tg_yyz_xyy_1, tg_yyz_xyz_0, tg_yyz_xyz_1, tg_yyz_xzz_0, \
                                         tg_yyz_xzz_1, tg_yyz_yyy_0, tg_yyz_yyy_1, tg_yyz_yyz_0, tg_yyz_yyz_1, tg_yyz_yzz_0, \
                                         tg_yyz_yzz_1, tg_yyz_zzz_0, tg_yyz_zzz_1, tg_yyzz_xx_1, tg_yyzz_xxx_0, \
                                         tg_yyzz_xxx_1, tg_yyzz_xxy_0, tg_yyzz_xxy_1, tg_yyzz_xxz_0, tg_yyzz_xxz_1, \
                                         tg_yyzz_xy_1, tg_yyzz_xyy_0, tg_yyzz_xyy_1, tg_yyzz_xyz_0, tg_yyzz_xyz_1, \
                                         tg_yyzz_xz_1, tg_yyzz_xzz_0, tg_yyzz_xzz_1, tg_yyzz_yy_1, tg_yyzz_yyy_0, \
                                         tg_yyzz_yyy_1, tg_yyzz_yyz_0, tg_yyzz_yyz_1, tg_yyzz_yz_1, tg_yyzz_yzz_0, \
                                         tg_yyzz_yzz_1, tg_yyzz_zz_1, tg_yyzz_zzz_0, tg_yyzz_zzz_1, tg_yyzzz_xxx_0, \
                                         tg_yyzzz_xxy_0, tg_yyzzz_xxz_0, tg_yyzzz_xyy_0, tg_yyzzz_xyz_0, tg_yyzzz_xzz_0, \
                                         tg_yyzzz_yyy_0, tg_yyzzz_yyz_0, tg_yyzzz_yzz_0, tg_yyzzz_zzz_0, tg_yzz_xxx_0, \
                                         tg_yzz_xxx_1, tg_yzz_xxy_0, tg_yzz_xxy_1, tg_yzz_xxz_0, tg_yzz_xxz_1, tg_yzz_xyy_0, \
                                         tg_yzz_xyy_1, tg_yzz_xyz_0, tg_yzz_xyz_1, tg_yzz_xzz_0, tg_yzz_xzz_1, tg_yzz_yyy_0, \
                                         tg_yzz_yyy_1, tg_yzz_yyz_0, tg_yzz_yyz_1, tg_yzz_yzz_0, tg_yzz_yzz_1, tg_yzz_zzz_0, \
                                         tg_yzz_zzz_1, tg_yzzz_xx_1, tg_yzzz_xxx_0, tg_yzzz_xxx_1, tg_yzzz_xxy_0, \
                                         tg_yzzz_xxy_1, tg_yzzz_xxz_0, tg_yzzz_xxz_1, tg_yzzz_xy_1, tg_yzzz_xyy_0, \
                                         tg_yzzz_xyy_1, tg_yzzz_xyz_0, tg_yzzz_xyz_1, tg_yzzz_xz_1, tg_yzzz_xzz_0, \
                                         tg_yzzz_xzz_1, tg_yzzz_yy_1, tg_yzzz_yyy_0, tg_yzzz_yyy_1, tg_yzzz_yyz_0, \
                                         tg_yzzz_yyz_1, tg_yzzz_yz_1, tg_yzzz_yzz_0, tg_yzzz_yzz_1, tg_yzzz_zz_1, \
                                         tg_yzzz_zzz_0, tg_yzzz_zzz_1, tg_yzzzz_xxx_0, tg_yzzzz_xxy_0, tg_yzzzz_xxz_0, \
                                         tg_yzzzz_xyy_0, tg_yzzzz_xyz_0, tg_yzzzz_xzz_0, tg_yzzzz_yyy_0, tg_yzzzz_yyz_0, \
                                         tg_yzzzz_yzz_0, tg_yzzzz_zzz_0, tg_zzz_xxx_0, tg_zzz_xxx_1, tg_zzz_xxy_0, \
                                         tg_zzz_xxy_1, tg_zzz_xxz_0, tg_zzz_xxz_1, tg_zzz_xyy_0, tg_zzz_xyy_1, tg_zzz_xyz_0, \
                                         tg_zzz_xyz_1, tg_zzz_xzz_0, tg_zzz_xzz_1, tg_zzz_yyy_0, tg_zzz_yyy_1, tg_zzz_yyz_0, \
                                         tg_zzz_yyz_1, tg_zzz_yzz_0, tg_zzz_yzz_1, tg_zzz_zzz_0, tg_zzz_zzz_1, tg_zzzz_xx_1, \
                                         tg_zzzz_xxx_0, tg_zzzz_xxx_1, tg_zzzz_xxy_0, tg_zzzz_xxy_1, tg_zzzz_xxz_0, \
                                         tg_zzzz_xxz_1, tg_zzzz_xy_1, tg_zzzz_xyy_0, tg_zzzz_xyy_1, tg_zzzz_xyz_0, \
                                         tg_zzzz_xyz_1, tg_zzzz_xz_1, tg_zzzz_xzz_0, tg_zzzz_xzz_1, tg_zzzz_yy_1, \
                                         tg_zzzz_yyy_0, tg_zzzz_yyy_1, tg_zzzz_yyz_0, tg_zzzz_yyz_1, tg_zzzz_yz_1, \
                                         tg_zzzz_yzz_0, tg_zzzz_yzz_1, tg_zzzz_zz_1, tg_zzzz_zzz_0, tg_zzzz_zzz_1, \
                                         tg_zzzzz_xxx_0, tg_zzzzz_xxy_0, tg_zzzzz_xxz_0, tg_zzzzz_xyy_0, tg_zzzzz_xyz_0, \
                                         tg_zzzzz_xzz_0, tg_zzzzz_yyy_0, tg_zzzzz_yyz_0, tg_zzzzz_yzz_0, tg_zzzzz_zzz_0, wp_x, \
                                         wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xzzzz_xxx_0[j] = pb_x * tg_zzzz_xxx_0[j] + wp_x[j] * tg_zzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_zzzz_xx_1[j];

                    tg_xzzzz_xxy_0[j] = pb_x * tg_zzzz_xxy_0[j] + wp_x[j] * tg_zzzz_xxy_1[j] + fl1_fxn * tg_zzzz_xy_1[j];

                    tg_xzzzz_xxz_0[j] = pb_x * tg_zzzz_xxz_0[j] + wp_x[j] * tg_zzzz_xxz_1[j] + fl1_fxn * tg_zzzz_xz_1[j];

                    tg_xzzzz_xyy_0[j] = pb_x * tg_zzzz_xyy_0[j] + wp_x[j] * tg_zzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_zzzz_yy_1[j];

                    tg_xzzzz_xyz_0[j] = pb_x * tg_zzzz_xyz_0[j] + wp_x[j] * tg_zzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yz_1[j];

                    tg_xzzzz_xzz_0[j] = pb_x * tg_zzzz_xzz_0[j] + wp_x[j] * tg_zzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zz_1[j];

                    tg_xzzzz_yyy_0[j] = pb_x * tg_zzzz_yyy_0[j] + wp_x[j] * tg_zzzz_yyy_1[j];

                    tg_xzzzz_yyz_0[j] = pb_x * tg_zzzz_yyz_0[j] + wp_x[j] * tg_zzzz_yyz_1[j];

                    tg_xzzzz_yzz_0[j] = pb_x * tg_zzzz_yzz_0[j] + wp_x[j] * tg_zzzz_yzz_1[j];

                    tg_xzzzz_zzz_0[j] = pb_x * tg_zzzz_zzz_0[j] + wp_x[j] * tg_zzzz_zzz_1[j];

                    tg_yyyyy_xxx_0[j] = pb_y * tg_yyyy_xxx_0[j] + wp_y[j] * tg_yyyy_xxx_1[j] + 2.0 * fl1_fx * tg_yyy_xxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxx_1[j];

                    tg_yyyyy_xxy_0[j] = pb_y * tg_yyyy_xxy_0[j] + wp_y[j] * tg_yyyy_xxy_1[j] + 2.0 * fl1_fx * tg_yyy_xxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxy_1[j] + 0.5 * fl1_fxn * tg_yyyy_xx_1[j];

                    tg_yyyyy_xxz_0[j] = pb_y * tg_yyyy_xxz_0[j] + wp_y[j] * tg_yyyy_xxz_1[j] + 2.0 * fl1_fx * tg_yyy_xxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xxz_1[j];

                    tg_yyyyy_xyy_0[j] = pb_y * tg_yyyy_xyy_0[j] + wp_y[j] * tg_yyyy_xyy_1[j] + 2.0 * fl1_fx * tg_yyy_xyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyy_1[j] + fl1_fxn * tg_yyyy_xy_1[j];

                    tg_yyyyy_xyz_0[j] = pb_y * tg_yyyy_xyz_0[j] + wp_y[j] * tg_yyyy_xyz_1[j] + 2.0 * fl1_fx * tg_yyy_xyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyy_xz_1[j];

                    tg_yyyyy_xzz_0[j] = pb_y * tg_yyyy_xzz_0[j] + wp_y[j] * tg_yyyy_xzz_1[j] + 2.0 * fl1_fx * tg_yyy_xzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_xzz_1[j];

                    tg_yyyyy_yyy_0[j] = pb_y * tg_yyyy_yyy_0[j] + wp_y[j] * tg_yyyy_yyy_1[j] + 2.0 * fl1_fx * tg_yyy_yyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyy_1[j] + 1.5 * fl1_fxn * tg_yyyy_yy_1[j];

                    tg_yyyyy_yyz_0[j] = pb_y * tg_yyyy_yyz_0[j] + wp_y[j] * tg_yyyy_yyz_1[j] + 2.0 * fl1_fx * tg_yyy_yyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yyz_1[j] + fl1_fxn * tg_yyyy_yz_1[j];

                    tg_yyyyy_yzz_0[j] = pb_y * tg_yyyy_yzz_0[j] + wp_y[j] * tg_yyyy_yzz_1[j] + 2.0 * fl1_fx * tg_yyy_yzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_yzz_1[j] + 0.5 * fl1_fxn * tg_yyyy_zz_1[j];

                    tg_yyyyy_zzz_0[j] = pb_y * tg_yyyy_zzz_0[j] + wp_y[j] * tg_yyyy_zzz_1[j] + 2.0 * fl1_fx * tg_yyy_zzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_zzz_1[j];

                    tg_yyyyz_xxx_0[j] = pb_y * tg_yyyz_xxx_0[j] + wp_y[j] * tg_yyyz_xxx_1[j] + 1.5 * fl1_fx * tg_yyz_xxx_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxx_1[j];

                    tg_yyyyz_xxy_0[j] = pb_y * tg_yyyz_xxy_0[j] + wp_y[j] * tg_yyyz_xxy_1[j] + 1.5 * fl1_fx * tg_yyz_xxy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxy_1[j] + 0.5 * fl1_fxn * tg_yyyz_xx_1[j];

                    tg_yyyyz_xxz_0[j] = pb_y * tg_yyyz_xxz_0[j] + wp_y[j] * tg_yyyz_xxz_1[j] + 1.5 * fl1_fx * tg_yyz_xxz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xxz_1[j];

                    tg_yyyyz_xyy_0[j] = pb_y * tg_yyyz_xyy_0[j] + wp_y[j] * tg_yyyz_xyy_1[j] + 1.5 * fl1_fx * tg_yyz_xyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyy_1[j] + fl1_fxn * tg_yyyz_xy_1[j];

                    tg_yyyyz_xyz_0[j] = pb_y * tg_yyyz_xyz_0[j] + wp_y[j] * tg_yyyz_xyz_1[j] + 1.5 * fl1_fx * tg_yyz_xyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyz_xz_1[j];

                    tg_yyyyz_xzz_0[j] = pb_y * tg_yyyz_xzz_0[j] + wp_y[j] * tg_yyyz_xzz_1[j] + 1.5 * fl1_fx * tg_yyz_xzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_xzz_1[j];

                    tg_yyyyz_yyy_0[j] = pb_y * tg_yyyz_yyy_0[j] + wp_y[j] * tg_yyyz_yyy_1[j] + 1.5 * fl1_fx * tg_yyz_yyy_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyy_1[j] + 1.5 * fl1_fxn * tg_yyyz_yy_1[j];

                    tg_yyyyz_yyz_0[j] = pb_y * tg_yyyz_yyz_0[j] + wp_y[j] * tg_yyyz_yyz_1[j] + 1.5 * fl1_fx * tg_yyz_yyz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yyz_1[j] + fl1_fxn * tg_yyyz_yz_1[j];

                    tg_yyyyz_yzz_0[j] = pb_y * tg_yyyz_yzz_0[j] + wp_y[j] * tg_yyyz_yzz_1[j] + 1.5 * fl1_fx * tg_yyz_yzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_yzz_1[j] + 0.5 * fl1_fxn * tg_yyyz_zz_1[j];

                    tg_yyyyz_zzz_0[j] = pb_y * tg_yyyz_zzz_0[j] + wp_y[j] * tg_yyyz_zzz_1[j] + 1.5 * fl1_fx * tg_yyz_zzz_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_zzz_1[j];

                    tg_yyyzz_xxx_0[j] = pb_y * tg_yyzz_xxx_0[j] + wp_y[j] * tg_yyzz_xxx_1[j] + fl1_fx * tg_yzz_xxx_0[j] - fl1_fx * fl1_fza * tg_yzz_xxx_1[j];

                    tg_yyyzz_xxy_0[j] = pb_y * tg_yyzz_xxy_0[j] + wp_y[j] * tg_yyzz_xxy_1[j] + fl1_fx * tg_yzz_xxy_0[j] - fl1_fx * fl1_fza * tg_yzz_xxy_1[j] + 0.5 * fl1_fxn * tg_yyzz_xx_1[j];

                    tg_yyyzz_xxz_0[j] = pb_y * tg_yyzz_xxz_0[j] + wp_y[j] * tg_yyzz_xxz_1[j] + fl1_fx * tg_yzz_xxz_0[j] - fl1_fx * fl1_fza * tg_yzz_xxz_1[j];

                    tg_yyyzz_xyy_0[j] = pb_y * tg_yyzz_xyy_0[j] + wp_y[j] * tg_yyzz_xyy_1[j] + fl1_fx * tg_yzz_xyy_0[j] - fl1_fx * fl1_fza * tg_yzz_xyy_1[j] + fl1_fxn * tg_yyzz_xy_1[j];

                    tg_yyyzz_xyz_0[j] = pb_y * tg_yyzz_xyz_0[j] + wp_y[j] * tg_yyzz_xyz_1[j] + fl1_fx * tg_yzz_xyz_0[j] - fl1_fx * fl1_fza * tg_yzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyzz_xz_1[j];

                    tg_yyyzz_xzz_0[j] = pb_y * tg_yyzz_xzz_0[j] + wp_y[j] * tg_yyzz_xzz_1[j] + fl1_fx * tg_yzz_xzz_0[j] - fl1_fx * fl1_fza * tg_yzz_xzz_1[j];

                    tg_yyyzz_yyy_0[j] = pb_y * tg_yyzz_yyy_0[j] + wp_y[j] * tg_yyzz_yyy_1[j] + fl1_fx * tg_yzz_yyy_0[j] - fl1_fx * fl1_fza * tg_yzz_yyy_1[j] + 1.5 * fl1_fxn * tg_yyzz_yy_1[j];

                    tg_yyyzz_yyz_0[j] = pb_y * tg_yyzz_yyz_0[j] + wp_y[j] * tg_yyzz_yyz_1[j] + fl1_fx * tg_yzz_yyz_0[j] - fl1_fx * fl1_fza * tg_yzz_yyz_1[j] + fl1_fxn * tg_yyzz_yz_1[j];

                    tg_yyyzz_yzz_0[j] = pb_y * tg_yyzz_yzz_0[j] + wp_y[j] * tg_yyzz_yzz_1[j] + fl1_fx * tg_yzz_yzz_0[j] - fl1_fx * fl1_fza * tg_yzz_yzz_1[j] + 0.5 * fl1_fxn * tg_yyzz_zz_1[j];

                    tg_yyyzz_zzz_0[j] = pb_y * tg_yyzz_zzz_0[j] + wp_y[j] * tg_yyzz_zzz_1[j] + fl1_fx * tg_yzz_zzz_0[j] - fl1_fx * fl1_fza * tg_yzz_zzz_1[j];

                    tg_yyzzz_xxx_0[j] = pb_y * tg_yzzz_xxx_0[j] + wp_y[j] * tg_yzzz_xxx_1[j] + 0.5 * fl1_fx * tg_zzz_xxx_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxx_1[j];

                    tg_yyzzz_xxy_0[j] = pb_y * tg_yzzz_xxy_0[j] + wp_y[j] * tg_yzzz_xxy_1[j] + 0.5 * fl1_fx * tg_zzz_xxy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxy_1[j] + 0.5 * fl1_fxn * tg_yzzz_xx_1[j];

                    tg_yyzzz_xxz_0[j] = pb_y * tg_yzzz_xxz_0[j] + wp_y[j] * tg_yzzz_xxz_1[j] + 0.5 * fl1_fx * tg_zzz_xxz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xxz_1[j];

                    tg_yyzzz_xyy_0[j] = pb_y * tg_yzzz_xyy_0[j] + wp_y[j] * tg_yzzz_xyy_1[j] + 0.5 * fl1_fx * tg_zzz_xyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyy_1[j] + fl1_fxn * tg_yzzz_xy_1[j];

                    tg_yyzzz_xyz_0[j] = pb_y * tg_yzzz_xyz_0[j] + wp_y[j] * tg_yzzz_xyz_1[j] + 0.5 * fl1_fx * tg_zzz_xyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yzzz_xz_1[j];

                    tg_yyzzz_xzz_0[j] = pb_y * tg_yzzz_xzz_0[j] + wp_y[j] * tg_yzzz_xzz_1[j] + 0.5 * fl1_fx * tg_zzz_xzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_xzz_1[j];

                    tg_yyzzz_yyy_0[j] = pb_y * tg_yzzz_yyy_0[j] + wp_y[j] * tg_yzzz_yyy_1[j] + 0.5 * fl1_fx * tg_zzz_yyy_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyy_1[j] + 1.5 * fl1_fxn * tg_yzzz_yy_1[j];

                    tg_yyzzz_yyz_0[j] = pb_y * tg_yzzz_yyz_0[j] + wp_y[j] * tg_yzzz_yyz_1[j] + 0.5 * fl1_fx * tg_zzz_yyz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yyz_1[j] + fl1_fxn * tg_yzzz_yz_1[j];

                    tg_yyzzz_yzz_0[j] = pb_y * tg_yzzz_yzz_0[j] + wp_y[j] * tg_yzzz_yzz_1[j] + 0.5 * fl1_fx * tg_zzz_yzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_yzz_1[j] + 0.5 * fl1_fxn * tg_yzzz_zz_1[j];

                    tg_yyzzz_zzz_0[j] = pb_y * tg_yzzz_zzz_0[j] + wp_y[j] * tg_yzzz_zzz_1[j] + 0.5 * fl1_fx * tg_zzz_zzz_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_zzz_1[j];

                    tg_yzzzz_xxx_0[j] = pb_y * tg_zzzz_xxx_0[j] + wp_y[j] * tg_zzzz_xxx_1[j];

                    tg_yzzzz_xxy_0[j] = pb_y * tg_zzzz_xxy_0[j] + wp_y[j] * tg_zzzz_xxy_1[j] + 0.5 * fl1_fxn * tg_zzzz_xx_1[j];

                    tg_yzzzz_xxz_0[j] = pb_y * tg_zzzz_xxz_0[j] + wp_y[j] * tg_zzzz_xxz_1[j];

                    tg_yzzzz_xyy_0[j] = pb_y * tg_zzzz_xyy_0[j] + wp_y[j] * tg_zzzz_xyy_1[j] + fl1_fxn * tg_zzzz_xy_1[j];

                    tg_yzzzz_xyz_0[j] = pb_y * tg_zzzz_xyz_0[j] + wp_y[j] * tg_zzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xz_1[j];

                    tg_yzzzz_xzz_0[j] = pb_y * tg_zzzz_xzz_0[j] + wp_y[j] * tg_zzzz_xzz_1[j];

                    tg_yzzzz_yyy_0[j] = pb_y * tg_zzzz_yyy_0[j] + wp_y[j] * tg_zzzz_yyy_1[j] + 1.5 * fl1_fxn * tg_zzzz_yy_1[j];

                    tg_yzzzz_yyz_0[j] = pb_y * tg_zzzz_yyz_0[j] + wp_y[j] * tg_zzzz_yyz_1[j] + fl1_fxn * tg_zzzz_yz_1[j];

                    tg_yzzzz_yzz_0[j] = pb_y * tg_zzzz_yzz_0[j] + wp_y[j] * tg_zzzz_yzz_1[j] + 0.5 * fl1_fxn * tg_zzzz_zz_1[j];

                    tg_yzzzz_zzz_0[j] = pb_y * tg_zzzz_zzz_0[j] + wp_y[j] * tg_zzzz_zzz_1[j];

                    tg_zzzzz_xxx_0[j] = pb_z * tg_zzzz_xxx_0[j] + wp_z[j] * tg_zzzz_xxx_1[j] + 2.0 * fl1_fx * tg_zzz_xxx_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxx_1[j];

                    tg_zzzzz_xxy_0[j] = pb_z * tg_zzzz_xxy_0[j] + wp_z[j] * tg_zzzz_xxy_1[j] + 2.0 * fl1_fx * tg_zzz_xxy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxy_1[j];

                    tg_zzzzz_xxz_0[j] = pb_z * tg_zzzz_xxz_0[j] + wp_z[j] * tg_zzzz_xxz_1[j] + 2.0 * fl1_fx * tg_zzz_xxz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xxz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xx_1[j];

                    tg_zzzzz_xyy_0[j] = pb_z * tg_zzzz_xyy_0[j] + wp_z[j] * tg_zzzz_xyy_1[j] + 2.0 * fl1_fx * tg_zzz_xyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyy_1[j];

                    tg_zzzzz_xyz_0[j] = pb_z * tg_zzzz_xyz_0[j] + wp_z[j] * tg_zzzz_xyz_1[j] + 2.0 * fl1_fx * tg_zzz_xyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_xy_1[j];

                    tg_zzzzz_xzz_0[j] = pb_z * tg_zzzz_xzz_0[j] + wp_z[j] * tg_zzzz_xzz_1[j] + 2.0 * fl1_fx * tg_zzz_xzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_xzz_1[j] + fl1_fxn * tg_zzzz_xz_1[j];

                    tg_zzzzz_yyy_0[j] = pb_z * tg_zzzz_yyy_0[j] + wp_z[j] * tg_zzzz_yyy_1[j] + 2.0 * fl1_fx * tg_zzz_yyy_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyy_1[j];

                    tg_zzzzz_yyz_0[j] = pb_z * tg_zzzz_yyz_0[j] + wp_z[j] * tg_zzzz_yyz_1[j] + 2.0 * fl1_fx * tg_zzz_yyz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yyz_1[j] + 0.5 * fl1_fxn * tg_zzzz_yy_1[j];

                    tg_zzzzz_yzz_0[j] = pb_z * tg_zzzz_yzz_0[j] + wp_z[j] * tg_zzzz_yzz_1[j] + 2.0 * fl1_fx * tg_zzz_yzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_yzz_1[j] + fl1_fxn * tg_zzzz_yz_1[j];

                    tg_zzzzz_zzz_0[j] = pb_z * tg_zzzz_zzz_0[j] + wp_z[j] * tg_zzzz_zzz_1[j] + 2.0 * fl1_fx * tg_zzz_zzz_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_zzz_1[j] + 1.5 * fl1_fxn * tg_zzzz_zz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

