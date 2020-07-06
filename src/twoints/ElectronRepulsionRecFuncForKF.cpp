//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForKF.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSKSF(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSKSF_0_90(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSKSF_90_180(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSKSF_180_270(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSKSF_270_360(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSKSF_0_90(      CMemBlock2D<double>* primBuffer,
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
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xxxxxx_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx); 

                auto tg_xxxxxx_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 1); 

                auto tg_xxxxxx_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 2); 

                auto tg_xxxxxx_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 3); 

                auto tg_xxxxxx_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 4); 

                auto tg_xxxxxx_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 5); 

                auto tg_xxxxxx_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 6); 

                auto tg_xxxxxx_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 7); 

                auto tg_xxxxxx_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 8); 

                auto tg_xxxxxx_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 9); 

                auto tg_xxxxxy_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 10); 

                auto tg_xxxxxy_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 11); 

                auto tg_xxxxxy_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 12); 

                auto tg_xxxxxy_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 13); 

                auto tg_xxxxxy_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 14); 

                auto tg_xxxxxy_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 15); 

                auto tg_xxxxxy_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 16); 

                auto tg_xxxxxy_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 17); 

                auto tg_xxxxxy_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 18); 

                auto tg_xxxxxy_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 19); 

                auto tg_xxxxxz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 20); 

                auto tg_xxxxxz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 21); 

                auto tg_xxxxxz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 22); 

                auto tg_xxxxxz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 23); 

                auto tg_xxxxxz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 24); 

                auto tg_xxxxxz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 25); 

                auto tg_xxxxxz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 26); 

                auto tg_xxxxxz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 27); 

                auto tg_xxxxxz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 28); 

                auto tg_xxxxxz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 29); 

                auto tg_xxxxyy_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 30); 

                auto tg_xxxxyy_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 31); 

                auto tg_xxxxyy_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 32); 

                auto tg_xxxxyy_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 33); 

                auto tg_xxxxyy_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 34); 

                auto tg_xxxxyy_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 35); 

                auto tg_xxxxyy_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 36); 

                auto tg_xxxxyy_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 37); 

                auto tg_xxxxyy_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 38); 

                auto tg_xxxxyy_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 39); 

                auto tg_xxxxyz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 40); 

                auto tg_xxxxyz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 41); 

                auto tg_xxxxyz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 42); 

                auto tg_xxxxyz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 43); 

                auto tg_xxxxyz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 44); 

                auto tg_xxxxyz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 45); 

                auto tg_xxxxyz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 46); 

                auto tg_xxxxyz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 47); 

                auto tg_xxxxyz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 48); 

                auto tg_xxxxyz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 49); 

                auto tg_xxxxzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 50); 

                auto tg_xxxxzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 51); 

                auto tg_xxxxzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 52); 

                auto tg_xxxxzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 53); 

                auto tg_xxxxzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 54); 

                auto tg_xxxxzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 55); 

                auto tg_xxxxzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 56); 

                auto tg_xxxxzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 57); 

                auto tg_xxxxzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 58); 

                auto tg_xxxxzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 59); 

                auto tg_xxxyyy_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 60); 

                auto tg_xxxyyy_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 61); 

                auto tg_xxxyyy_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 62); 

                auto tg_xxxyyy_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 63); 

                auto tg_xxxyyy_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 64); 

                auto tg_xxxyyy_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 65); 

                auto tg_xxxyyy_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 66); 

                auto tg_xxxyyy_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 67); 

                auto tg_xxxyyy_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 68); 

                auto tg_xxxyyy_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 69); 

                auto tg_xxxyyz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 70); 

                auto tg_xxxyyz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 71); 

                auto tg_xxxyyz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 72); 

                auto tg_xxxyyz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 73); 

                auto tg_xxxyyz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 74); 

                auto tg_xxxyyz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 75); 

                auto tg_xxxyyz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 76); 

                auto tg_xxxyyz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 77); 

                auto tg_xxxyyz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 78); 

                auto tg_xxxyyz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 79); 

                auto tg_xxxyzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 80); 

                auto tg_xxxyzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 81); 

                auto tg_xxxyzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 82); 

                auto tg_xxxyzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 83); 

                auto tg_xxxyzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 84); 

                auto tg_xxxyzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 85); 

                auto tg_xxxyzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 86); 

                auto tg_xxxyzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 87); 

                auto tg_xxxyzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 88); 

                auto tg_xxxyzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 89); 

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

                auto tg_xxxxx_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx); 

                auto tg_xxxxx_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 1); 

                auto tg_xxxxx_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 2); 

                auto tg_xxxxx_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 3); 

                auto tg_xxxxx_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 4); 

                auto tg_xxxxx_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 5); 

                auto tg_xxxxx_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 6); 

                auto tg_xxxxx_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 7); 

                auto tg_xxxxx_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 8); 

                auto tg_xxxxx_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 9); 

                auto tg_xxxxy_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 10); 

                auto tg_xxxxy_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 11); 

                auto tg_xxxxy_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 12); 

                auto tg_xxxxy_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 13); 

                auto tg_xxxxy_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 14); 

                auto tg_xxxxy_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 15); 

                auto tg_xxxxy_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 16); 

                auto tg_xxxxy_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 17); 

                auto tg_xxxxy_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 18); 

                auto tg_xxxxy_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 19); 

                auto tg_xxxxz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 20); 

                auto tg_xxxxz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 21); 

                auto tg_xxxxz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 22); 

                auto tg_xxxxz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 23); 

                auto tg_xxxxz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 24); 

                auto tg_xxxxz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 25); 

                auto tg_xxxxz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 26); 

                auto tg_xxxxz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 27); 

                auto tg_xxxxz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 28); 

                auto tg_xxxxz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 29); 

                auto tg_xxxyy_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 30); 

                auto tg_xxxyy_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 31); 

                auto tg_xxxyy_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 32); 

                auto tg_xxxyy_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 33); 

                auto tg_xxxyy_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 34); 

                auto tg_xxxyy_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 35); 

                auto tg_xxxyy_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 36); 

                auto tg_xxxyy_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 37); 

                auto tg_xxxyy_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 38); 

                auto tg_xxxyy_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 39); 

                auto tg_xxxyz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 40); 

                auto tg_xxxyz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 41); 

                auto tg_xxxyz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 42); 

                auto tg_xxxyz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 43); 

                auto tg_xxxyz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 44); 

                auto tg_xxxyz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 45); 

                auto tg_xxxyz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 46); 

                auto tg_xxxyz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 47); 

                auto tg_xxxyz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 48); 

                auto tg_xxxyz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 49); 

                auto tg_xxxzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 50); 

                auto tg_xxxzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 51); 

                auto tg_xxxzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 52); 

                auto tg_xxxzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 53); 

                auto tg_xxxzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 54); 

                auto tg_xxxzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 55); 

                auto tg_xxxzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 56); 

                auto tg_xxxzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 57); 

                auto tg_xxxzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 58); 

                auto tg_xxxzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 59); 

                auto tg_xxyyy_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 60); 

                auto tg_xxyyy_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 61); 

                auto tg_xxyyy_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 62); 

                auto tg_xxyyy_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 63); 

                auto tg_xxyyy_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 64); 

                auto tg_xxyyy_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 65); 

                auto tg_xxyyy_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 66); 

                auto tg_xxyyy_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 67); 

                auto tg_xxyyy_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 68); 

                auto tg_xxyyy_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 69); 

                auto tg_xxyyz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 70); 

                auto tg_xxyyz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 71); 

                auto tg_xxyyz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 72); 

                auto tg_xxyyz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 73); 

                auto tg_xxyyz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 74); 

                auto tg_xxyyz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 75); 

                auto tg_xxyyz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 76); 

                auto tg_xxyyz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 77); 

                auto tg_xxyyz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 78); 

                auto tg_xxyyz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 79); 

                auto tg_xxyzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 80); 

                auto tg_xxyzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 81); 

                auto tg_xxyzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 82); 

                auto tg_xxyzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 83); 

                auto tg_xxyzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 84); 

                auto tg_xxyzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 85); 

                auto tg_xxyzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 86); 

                auto tg_xxyzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 87); 

                auto tg_xxyzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 88); 

                auto tg_xxyzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 89); 

                auto tg_xxxxx_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx); 

                auto tg_xxxxx_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 1); 

                auto tg_xxxxx_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 2); 

                auto tg_xxxxx_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 3); 

                auto tg_xxxxx_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 4); 

                auto tg_xxxxx_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 5); 

                auto tg_xxxxx_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 6); 

                auto tg_xxxxx_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 7); 

                auto tg_xxxxx_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 8); 

                auto tg_xxxxx_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 9); 

                auto tg_xxxxy_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 10); 

                auto tg_xxxxy_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 11); 

                auto tg_xxxxy_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 12); 

                auto tg_xxxxy_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 13); 

                auto tg_xxxxy_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 14); 

                auto tg_xxxxy_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 15); 

                auto tg_xxxxy_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 16); 

                auto tg_xxxxy_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 17); 

                auto tg_xxxxy_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 18); 

                auto tg_xxxxy_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 19); 

                auto tg_xxxxz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 20); 

                auto tg_xxxxz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 21); 

                auto tg_xxxxz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 22); 

                auto tg_xxxxz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 23); 

                auto tg_xxxxz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 24); 

                auto tg_xxxxz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 25); 

                auto tg_xxxxz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 26); 

                auto tg_xxxxz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 27); 

                auto tg_xxxxz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 28); 

                auto tg_xxxxz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 29); 

                auto tg_xxxyy_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 30); 

                auto tg_xxxyy_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 31); 

                auto tg_xxxyy_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 32); 

                auto tg_xxxyy_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 33); 

                auto tg_xxxyy_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 34); 

                auto tg_xxxyy_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 35); 

                auto tg_xxxyy_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 36); 

                auto tg_xxxyy_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 37); 

                auto tg_xxxyy_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 38); 

                auto tg_xxxyy_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 39); 

                auto tg_xxxyz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 40); 

                auto tg_xxxyz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 41); 

                auto tg_xxxyz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 42); 

                auto tg_xxxyz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 43); 

                auto tg_xxxyz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 44); 

                auto tg_xxxyz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 45); 

                auto tg_xxxyz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 46); 

                auto tg_xxxyz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 47); 

                auto tg_xxxyz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 48); 

                auto tg_xxxyz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 49); 

                auto tg_xxxzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 50); 

                auto tg_xxxzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 51); 

                auto tg_xxxzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 52); 

                auto tg_xxxzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 53); 

                auto tg_xxxzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 54); 

                auto tg_xxxzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 55); 

                auto tg_xxxzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 56); 

                auto tg_xxxzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 57); 

                auto tg_xxxzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 58); 

                auto tg_xxxzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 59); 

                auto tg_xxyyy_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 60); 

                auto tg_xxyyy_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 61); 

                auto tg_xxyyy_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 62); 

                auto tg_xxyyy_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 63); 

                auto tg_xxyyy_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 64); 

                auto tg_xxyyy_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 65); 

                auto tg_xxyyy_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 66); 

                auto tg_xxyyy_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 67); 

                auto tg_xxyyy_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 68); 

                auto tg_xxyyy_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 69); 

                auto tg_xxyyz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 70); 

                auto tg_xxyyz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 71); 

                auto tg_xxyyz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 72); 

                auto tg_xxyyz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 73); 

                auto tg_xxyyz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 74); 

                auto tg_xxyyz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 75); 

                auto tg_xxyyz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 76); 

                auto tg_xxyyz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 77); 

                auto tg_xxyyz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 78); 

                auto tg_xxyyz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 79); 

                auto tg_xxyzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 80); 

                auto tg_xxyzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 81); 

                auto tg_xxyzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 82); 

                auto tg_xxyzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 83); 

                auto tg_xxyzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 84); 

                auto tg_xxyzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 85); 

                auto tg_xxyzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 86); 

                auto tg_xxyzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 87); 

                auto tg_xxyzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 88); 

                auto tg_xxyzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 89); 

                auto tg_xxxxxx_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx); 

                auto tg_xxxxxx_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 1); 

                auto tg_xxxxxx_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 2); 

                auto tg_xxxxxx_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 3); 

                auto tg_xxxxxx_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 4); 

                auto tg_xxxxxx_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 5); 

                auto tg_xxxxxy_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 6); 

                auto tg_xxxxxy_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 7); 

                auto tg_xxxxxy_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 8); 

                auto tg_xxxxxy_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 9); 

                auto tg_xxxxxy_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 10); 

                auto tg_xxxxxy_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 11); 

                auto tg_xxxxxz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 12); 

                auto tg_xxxxxz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 13); 

                auto tg_xxxxxz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 14); 

                auto tg_xxxxxz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 15); 

                auto tg_xxxxxz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 16); 

                auto tg_xxxxxz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 17); 

                auto tg_xxxxyy_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 18); 

                auto tg_xxxxyy_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 19); 

                auto tg_xxxxyy_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 20); 

                auto tg_xxxxyy_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 21); 

                auto tg_xxxxyy_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 22); 

                auto tg_xxxxyy_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 23); 

                auto tg_xxxxyz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 24); 

                auto tg_xxxxyz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 25); 

                auto tg_xxxxyz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 26); 

                auto tg_xxxxyz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 27); 

                auto tg_xxxxyz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 28); 

                auto tg_xxxxyz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 29); 

                auto tg_xxxxzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 30); 

                auto tg_xxxxzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 31); 

                auto tg_xxxxzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 32); 

                auto tg_xxxxzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 33); 

                auto tg_xxxxzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 34); 

                auto tg_xxxxzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 35); 

                auto tg_xxxyyy_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 36); 

                auto tg_xxxyyy_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 37); 

                auto tg_xxxyyy_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 38); 

                auto tg_xxxyyy_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 39); 

                auto tg_xxxyyy_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 40); 

                auto tg_xxxyyy_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 41); 

                auto tg_xxxyyz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 42); 

                auto tg_xxxyyz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 43); 

                auto tg_xxxyyz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 44); 

                auto tg_xxxyyz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 45); 

                auto tg_xxxyyz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 46); 

                auto tg_xxxyyz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 47); 

                auto tg_xxxyzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 48); 

                auto tg_xxxyzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 49); 

                auto tg_xxxyzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 50); 

                auto tg_xxxyzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 51); 

                auto tg_xxxyzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 52); 

                auto tg_xxxyzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 53); 

                // set up pointers to integrals

                auto tg_xxxxxxx_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx); 

                auto tg_xxxxxxx_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 1); 

                auto tg_xxxxxxx_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 2); 

                auto tg_xxxxxxx_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 3); 

                auto tg_xxxxxxx_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 4); 

                auto tg_xxxxxxx_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 5); 

                auto tg_xxxxxxx_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 6); 

                auto tg_xxxxxxx_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 7); 

                auto tg_xxxxxxx_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 8); 

                auto tg_xxxxxxx_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 9); 

                auto tg_xxxxxxy_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 10); 

                auto tg_xxxxxxy_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 11); 

                auto tg_xxxxxxy_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 12); 

                auto tg_xxxxxxy_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 13); 

                auto tg_xxxxxxy_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 14); 

                auto tg_xxxxxxy_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 15); 

                auto tg_xxxxxxy_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 16); 

                auto tg_xxxxxxy_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 17); 

                auto tg_xxxxxxy_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 18); 

                auto tg_xxxxxxy_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 19); 

                auto tg_xxxxxxz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 20); 

                auto tg_xxxxxxz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 21); 

                auto tg_xxxxxxz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 22); 

                auto tg_xxxxxxz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 23); 

                auto tg_xxxxxxz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 24); 

                auto tg_xxxxxxz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 25); 

                auto tg_xxxxxxz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 26); 

                auto tg_xxxxxxz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 27); 

                auto tg_xxxxxxz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 28); 

                auto tg_xxxxxxz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 29); 

                auto tg_xxxxxyy_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 30); 

                auto tg_xxxxxyy_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 31); 

                auto tg_xxxxxyy_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 32); 

                auto tg_xxxxxyy_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 33); 

                auto tg_xxxxxyy_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 34); 

                auto tg_xxxxxyy_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 35); 

                auto tg_xxxxxyy_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 36); 

                auto tg_xxxxxyy_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 37); 

                auto tg_xxxxxyy_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 38); 

                auto tg_xxxxxyy_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 39); 

                auto tg_xxxxxyz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 40); 

                auto tg_xxxxxyz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 41); 

                auto tg_xxxxxyz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 42); 

                auto tg_xxxxxyz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 43); 

                auto tg_xxxxxyz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 44); 

                auto tg_xxxxxyz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 45); 

                auto tg_xxxxxyz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 46); 

                auto tg_xxxxxyz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 47); 

                auto tg_xxxxxyz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 48); 

                auto tg_xxxxxyz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 49); 

                auto tg_xxxxxzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 50); 

                auto tg_xxxxxzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 51); 

                auto tg_xxxxxzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 52); 

                auto tg_xxxxxzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 53); 

                auto tg_xxxxxzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 54); 

                auto tg_xxxxxzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 55); 

                auto tg_xxxxxzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 56); 

                auto tg_xxxxxzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 57); 

                auto tg_xxxxxzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 58); 

                auto tg_xxxxxzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 59); 

                auto tg_xxxxyyy_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 60); 

                auto tg_xxxxyyy_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 61); 

                auto tg_xxxxyyy_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 62); 

                auto tg_xxxxyyy_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 63); 

                auto tg_xxxxyyy_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 64); 

                auto tg_xxxxyyy_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 65); 

                auto tg_xxxxyyy_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 66); 

                auto tg_xxxxyyy_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 67); 

                auto tg_xxxxyyy_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 68); 

                auto tg_xxxxyyy_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 69); 

                auto tg_xxxxyyz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 70); 

                auto tg_xxxxyyz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 71); 

                auto tg_xxxxyyz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 72); 

                auto tg_xxxxyyz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 73); 

                auto tg_xxxxyyz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 74); 

                auto tg_xxxxyyz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 75); 

                auto tg_xxxxyyz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 76); 

                auto tg_xxxxyyz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 77); 

                auto tg_xxxxyyz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 78); 

                auto tg_xxxxyyz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 79); 

                auto tg_xxxxyzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 80); 

                auto tg_xxxxyzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 81); 

                auto tg_xxxxyzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 82); 

                auto tg_xxxxyzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 83); 

                auto tg_xxxxyzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 84); 

                auto tg_xxxxyzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 85); 

                auto tg_xxxxyzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 86); 

                auto tg_xxxxyzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 87); 

                auto tg_xxxxyzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 88); 

                auto tg_xxxxyzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 89); 

                // Batch of Integrals (0,90)

                #pragma omp simd aligned(fxn, fza, tg_xxxxx_xxx_0, tg_xxxxx_xxx_1, tg_xxxxx_xxy_0, \
                                         tg_xxxxx_xxy_1, tg_xxxxx_xxz_0, tg_xxxxx_xxz_1, tg_xxxxx_xyy_0, tg_xxxxx_xyy_1, \
                                         tg_xxxxx_xyz_0, tg_xxxxx_xyz_1, tg_xxxxx_xzz_0, tg_xxxxx_xzz_1, tg_xxxxx_yyy_0, \
                                         tg_xxxxx_yyy_1, tg_xxxxx_yyz_0, tg_xxxxx_yyz_1, tg_xxxxx_yzz_0, tg_xxxxx_yzz_1, \
                                         tg_xxxxx_zzz_0, tg_xxxxx_zzz_1, tg_xxxxxx_xx_1, tg_xxxxxx_xxx_0, tg_xxxxxx_xxx_1, \
                                         tg_xxxxxx_xxy_0, tg_xxxxxx_xxy_1, tg_xxxxxx_xxz_0, tg_xxxxxx_xxz_1, tg_xxxxxx_xy_1, \
                                         tg_xxxxxx_xyy_0, tg_xxxxxx_xyy_1, tg_xxxxxx_xyz_0, tg_xxxxxx_xyz_1, tg_xxxxxx_xz_1, \
                                         tg_xxxxxx_xzz_0, tg_xxxxxx_xzz_1, tg_xxxxxx_yy_1, tg_xxxxxx_yyy_0, tg_xxxxxx_yyy_1, \
                                         tg_xxxxxx_yyz_0, tg_xxxxxx_yyz_1, tg_xxxxxx_yz_1, tg_xxxxxx_yzz_0, tg_xxxxxx_yzz_1, \
                                         tg_xxxxxx_zz_1, tg_xxxxxx_zzz_0, tg_xxxxxx_zzz_1, tg_xxxxxxx_xxx_0, \
                                         tg_xxxxxxx_xxy_0, tg_xxxxxxx_xxz_0, tg_xxxxxxx_xyy_0, tg_xxxxxxx_xyz_0, \
                                         tg_xxxxxxx_xzz_0, tg_xxxxxxx_yyy_0, tg_xxxxxxx_yyz_0, tg_xxxxxxx_yzz_0, \
                                         tg_xxxxxxx_zzz_0, tg_xxxxxxy_xxx_0, tg_xxxxxxy_xxy_0, tg_xxxxxxy_xxz_0, \
                                         tg_xxxxxxy_xyy_0, tg_xxxxxxy_xyz_0, tg_xxxxxxy_xzz_0, tg_xxxxxxy_yyy_0, \
                                         tg_xxxxxxy_yyz_0, tg_xxxxxxy_yzz_0, tg_xxxxxxy_zzz_0, tg_xxxxxxz_xxx_0, \
                                         tg_xxxxxxz_xxy_0, tg_xxxxxxz_xxz_0, tg_xxxxxxz_xyy_0, tg_xxxxxxz_xyz_0, \
                                         tg_xxxxxxz_xzz_0, tg_xxxxxxz_yyy_0, tg_xxxxxxz_yyz_0, tg_xxxxxxz_yzz_0, \
                                         tg_xxxxxxz_zzz_0, tg_xxxxxy_xx_1, tg_xxxxxy_xxx_0, tg_xxxxxy_xxx_1, tg_xxxxxy_xxy_0, \
                                         tg_xxxxxy_xxy_1, tg_xxxxxy_xxz_0, tg_xxxxxy_xxz_1, tg_xxxxxy_xy_1, tg_xxxxxy_xyy_0, \
                                         tg_xxxxxy_xyy_1, tg_xxxxxy_xyz_0, tg_xxxxxy_xyz_1, tg_xxxxxy_xz_1, tg_xxxxxy_xzz_0, \
                                         tg_xxxxxy_xzz_1, tg_xxxxxy_yy_1, tg_xxxxxy_yyy_0, tg_xxxxxy_yyy_1, tg_xxxxxy_yyz_0, \
                                         tg_xxxxxy_yyz_1, tg_xxxxxy_yz_1, tg_xxxxxy_yzz_0, tg_xxxxxy_yzz_1, tg_xxxxxy_zz_1, \
                                         tg_xxxxxy_zzz_0, tg_xxxxxy_zzz_1, tg_xxxxxyy_xxx_0, tg_xxxxxyy_xxy_0, \
                                         tg_xxxxxyy_xxz_0, tg_xxxxxyy_xyy_0, tg_xxxxxyy_xyz_0, tg_xxxxxyy_xzz_0, \
                                         tg_xxxxxyy_yyy_0, tg_xxxxxyy_yyz_0, tg_xxxxxyy_yzz_0, tg_xxxxxyy_zzz_0, \
                                         tg_xxxxxyz_xxx_0, tg_xxxxxyz_xxy_0, tg_xxxxxyz_xxz_0, tg_xxxxxyz_xyy_0, \
                                         tg_xxxxxyz_xyz_0, tg_xxxxxyz_xzz_0, tg_xxxxxyz_yyy_0, tg_xxxxxyz_yyz_0, \
                                         tg_xxxxxyz_yzz_0, tg_xxxxxyz_zzz_0, tg_xxxxxz_xx_1, tg_xxxxxz_xxx_0, tg_xxxxxz_xxx_1, \
                                         tg_xxxxxz_xxy_0, tg_xxxxxz_xxy_1, tg_xxxxxz_xxz_0, tg_xxxxxz_xxz_1, tg_xxxxxz_xy_1, \
                                         tg_xxxxxz_xyy_0, tg_xxxxxz_xyy_1, tg_xxxxxz_xyz_0, tg_xxxxxz_xyz_1, tg_xxxxxz_xz_1, \
                                         tg_xxxxxz_xzz_0, tg_xxxxxz_xzz_1, tg_xxxxxz_yy_1, tg_xxxxxz_yyy_0, tg_xxxxxz_yyy_1, \
                                         tg_xxxxxz_yyz_0, tg_xxxxxz_yyz_1, tg_xxxxxz_yz_1, tg_xxxxxz_yzz_0, tg_xxxxxz_yzz_1, \
                                         tg_xxxxxz_zz_1, tg_xxxxxz_zzz_0, tg_xxxxxz_zzz_1, tg_xxxxxzz_xxx_0, \
                                         tg_xxxxxzz_xxy_0, tg_xxxxxzz_xxz_0, tg_xxxxxzz_xyy_0, tg_xxxxxzz_xyz_0, \
                                         tg_xxxxxzz_xzz_0, tg_xxxxxzz_yyy_0, tg_xxxxxzz_yyz_0, tg_xxxxxzz_yzz_0, \
                                         tg_xxxxxzz_zzz_0, tg_xxxxy_xxx_0, tg_xxxxy_xxx_1, tg_xxxxy_xxy_0, tg_xxxxy_xxy_1, \
                                         tg_xxxxy_xxz_0, tg_xxxxy_xxz_1, tg_xxxxy_xyy_0, tg_xxxxy_xyy_1, tg_xxxxy_xyz_0, \
                                         tg_xxxxy_xyz_1, tg_xxxxy_xzz_0, tg_xxxxy_xzz_1, tg_xxxxy_yyy_0, tg_xxxxy_yyy_1, \
                                         tg_xxxxy_yyz_0, tg_xxxxy_yyz_1, tg_xxxxy_yzz_0, tg_xxxxy_yzz_1, tg_xxxxy_zzz_0, \
                                         tg_xxxxy_zzz_1, tg_xxxxyy_xx_1, tg_xxxxyy_xxx_0, tg_xxxxyy_xxx_1, tg_xxxxyy_xxy_0, \
                                         tg_xxxxyy_xxy_1, tg_xxxxyy_xxz_0, tg_xxxxyy_xxz_1, tg_xxxxyy_xy_1, tg_xxxxyy_xyy_0, \
                                         tg_xxxxyy_xyy_1, tg_xxxxyy_xyz_0, tg_xxxxyy_xyz_1, tg_xxxxyy_xz_1, tg_xxxxyy_xzz_0, \
                                         tg_xxxxyy_xzz_1, tg_xxxxyy_yy_1, tg_xxxxyy_yyy_0, tg_xxxxyy_yyy_1, tg_xxxxyy_yyz_0, \
                                         tg_xxxxyy_yyz_1, tg_xxxxyy_yz_1, tg_xxxxyy_yzz_0, tg_xxxxyy_yzz_1, tg_xxxxyy_zz_1, \
                                         tg_xxxxyy_zzz_0, tg_xxxxyy_zzz_1, tg_xxxxyyy_xxx_0, tg_xxxxyyy_xxy_0, \
                                         tg_xxxxyyy_xxz_0, tg_xxxxyyy_xyy_0, tg_xxxxyyy_xyz_0, tg_xxxxyyy_xzz_0, \
                                         tg_xxxxyyy_yyy_0, tg_xxxxyyy_yyz_0, tg_xxxxyyy_yzz_0, tg_xxxxyyy_zzz_0, \
                                         tg_xxxxyyz_xxx_0, tg_xxxxyyz_xxy_0, tg_xxxxyyz_xxz_0, tg_xxxxyyz_xyy_0, \
                                         tg_xxxxyyz_xyz_0, tg_xxxxyyz_xzz_0, tg_xxxxyyz_yyy_0, tg_xxxxyyz_yyz_0, \
                                         tg_xxxxyyz_yzz_0, tg_xxxxyyz_zzz_0, tg_xxxxyz_xx_1, tg_xxxxyz_xxx_0, tg_xxxxyz_xxx_1, \
                                         tg_xxxxyz_xxy_0, tg_xxxxyz_xxy_1, tg_xxxxyz_xxz_0, tg_xxxxyz_xxz_1, tg_xxxxyz_xy_1, \
                                         tg_xxxxyz_xyy_0, tg_xxxxyz_xyy_1, tg_xxxxyz_xyz_0, tg_xxxxyz_xyz_1, tg_xxxxyz_xz_1, \
                                         tg_xxxxyz_xzz_0, tg_xxxxyz_xzz_1, tg_xxxxyz_yy_1, tg_xxxxyz_yyy_0, tg_xxxxyz_yyy_1, \
                                         tg_xxxxyz_yyz_0, tg_xxxxyz_yyz_1, tg_xxxxyz_yz_1, tg_xxxxyz_yzz_0, tg_xxxxyz_yzz_1, \
                                         tg_xxxxyz_zz_1, tg_xxxxyz_zzz_0, tg_xxxxyz_zzz_1, tg_xxxxyzz_xxx_0, \
                                         tg_xxxxyzz_xxy_0, tg_xxxxyzz_xxz_0, tg_xxxxyzz_xyy_0, tg_xxxxyzz_xyz_0, \
                                         tg_xxxxyzz_xzz_0, tg_xxxxyzz_yyy_0, tg_xxxxyzz_yyz_0, tg_xxxxyzz_yzz_0, \
                                         tg_xxxxyzz_zzz_0, tg_xxxxz_xxx_0, tg_xxxxz_xxx_1, tg_xxxxz_xxy_0, tg_xxxxz_xxy_1, \
                                         tg_xxxxz_xxz_0, tg_xxxxz_xxz_1, tg_xxxxz_xyy_0, tg_xxxxz_xyy_1, tg_xxxxz_xyz_0, \
                                         tg_xxxxz_xyz_1, tg_xxxxz_xzz_0, tg_xxxxz_xzz_1, tg_xxxxz_yyy_0, tg_xxxxz_yyy_1, \
                                         tg_xxxxz_yyz_0, tg_xxxxz_yyz_1, tg_xxxxz_yzz_0, tg_xxxxz_yzz_1, tg_xxxxz_zzz_0, \
                                         tg_xxxxz_zzz_1, tg_xxxxzz_xx_1, tg_xxxxzz_xxx_0, tg_xxxxzz_xxx_1, tg_xxxxzz_xxy_0, \
                                         tg_xxxxzz_xxy_1, tg_xxxxzz_xxz_0, tg_xxxxzz_xxz_1, tg_xxxxzz_xy_1, tg_xxxxzz_xyy_0, \
                                         tg_xxxxzz_xyy_1, tg_xxxxzz_xyz_0, tg_xxxxzz_xyz_1, tg_xxxxzz_xz_1, tg_xxxxzz_xzz_0, \
                                         tg_xxxxzz_xzz_1, tg_xxxxzz_yy_1, tg_xxxxzz_yyy_0, tg_xxxxzz_yyy_1, tg_xxxxzz_yyz_0, \
                                         tg_xxxxzz_yyz_1, tg_xxxxzz_yz_1, tg_xxxxzz_yzz_0, tg_xxxxzz_yzz_1, tg_xxxxzz_zz_1, \
                                         tg_xxxxzz_zzz_0, tg_xxxxzz_zzz_1, tg_xxxyy_xxx_0, tg_xxxyy_xxx_1, tg_xxxyy_xxy_0, \
                                         tg_xxxyy_xxy_1, tg_xxxyy_xxz_0, tg_xxxyy_xxz_1, tg_xxxyy_xyy_0, tg_xxxyy_xyy_1, \
                                         tg_xxxyy_xyz_0, tg_xxxyy_xyz_1, tg_xxxyy_xzz_0, tg_xxxyy_xzz_1, tg_xxxyy_yyy_0, \
                                         tg_xxxyy_yyy_1, tg_xxxyy_yyz_0, tg_xxxyy_yyz_1, tg_xxxyy_yzz_0, tg_xxxyy_yzz_1, \
                                         tg_xxxyy_zzz_0, tg_xxxyy_zzz_1, tg_xxxyyy_xx_1, tg_xxxyyy_xxx_0, tg_xxxyyy_xxx_1, \
                                         tg_xxxyyy_xxy_0, tg_xxxyyy_xxy_1, tg_xxxyyy_xxz_0, tg_xxxyyy_xxz_1, tg_xxxyyy_xy_1, \
                                         tg_xxxyyy_xyy_0, tg_xxxyyy_xyy_1, tg_xxxyyy_xyz_0, tg_xxxyyy_xyz_1, tg_xxxyyy_xz_1, \
                                         tg_xxxyyy_xzz_0, tg_xxxyyy_xzz_1, tg_xxxyyy_yy_1, tg_xxxyyy_yyy_0, tg_xxxyyy_yyy_1, \
                                         tg_xxxyyy_yyz_0, tg_xxxyyy_yyz_1, tg_xxxyyy_yz_1, tg_xxxyyy_yzz_0, tg_xxxyyy_yzz_1, \
                                         tg_xxxyyy_zz_1, tg_xxxyyy_zzz_0, tg_xxxyyy_zzz_1, tg_xxxyyz_xx_1, tg_xxxyyz_xxx_0, \
                                         tg_xxxyyz_xxx_1, tg_xxxyyz_xxy_0, tg_xxxyyz_xxy_1, tg_xxxyyz_xxz_0, tg_xxxyyz_xxz_1, \
                                         tg_xxxyyz_xy_1, tg_xxxyyz_xyy_0, tg_xxxyyz_xyy_1, tg_xxxyyz_xyz_0, tg_xxxyyz_xyz_1, \
                                         tg_xxxyyz_xz_1, tg_xxxyyz_xzz_0, tg_xxxyyz_xzz_1, tg_xxxyyz_yy_1, tg_xxxyyz_yyy_0, \
                                         tg_xxxyyz_yyy_1, tg_xxxyyz_yyz_0, tg_xxxyyz_yyz_1, tg_xxxyyz_yz_1, tg_xxxyyz_yzz_0, \
                                         tg_xxxyyz_yzz_1, tg_xxxyyz_zz_1, tg_xxxyyz_zzz_0, tg_xxxyyz_zzz_1, tg_xxxyz_xxx_0, \
                                         tg_xxxyz_xxx_1, tg_xxxyz_xxy_0, tg_xxxyz_xxy_1, tg_xxxyz_xxz_0, tg_xxxyz_xxz_1, \
                                         tg_xxxyz_xyy_0, tg_xxxyz_xyy_1, tg_xxxyz_xyz_0, tg_xxxyz_xyz_1, tg_xxxyz_xzz_0, \
                                         tg_xxxyz_xzz_1, tg_xxxyz_yyy_0, tg_xxxyz_yyy_1, tg_xxxyz_yyz_0, tg_xxxyz_yyz_1, \
                                         tg_xxxyz_yzz_0, tg_xxxyz_yzz_1, tg_xxxyz_zzz_0, tg_xxxyz_zzz_1, tg_xxxyzz_xx_1, \
                                         tg_xxxyzz_xxx_0, tg_xxxyzz_xxx_1, tg_xxxyzz_xxy_0, tg_xxxyzz_xxy_1, tg_xxxyzz_xxz_0, \
                                         tg_xxxyzz_xxz_1, tg_xxxyzz_xy_1, tg_xxxyzz_xyy_0, tg_xxxyzz_xyy_1, tg_xxxyzz_xyz_0, \
                                         tg_xxxyzz_xyz_1, tg_xxxyzz_xz_1, tg_xxxyzz_xzz_0, tg_xxxyzz_xzz_1, tg_xxxyzz_yy_1, \
                                         tg_xxxyzz_yyy_0, tg_xxxyzz_yyy_1, tg_xxxyzz_yyz_0, tg_xxxyzz_yyz_1, tg_xxxyzz_yz_1, \
                                         tg_xxxyzz_yzz_0, tg_xxxyzz_yzz_1, tg_xxxyzz_zz_1, tg_xxxyzz_zzz_0, tg_xxxyzz_zzz_1, \
                                         tg_xxxzz_xxx_0, tg_xxxzz_xxx_1, tg_xxxzz_xxy_0, tg_xxxzz_xxy_1, tg_xxxzz_xxz_0, \
                                         tg_xxxzz_xxz_1, tg_xxxzz_xyy_0, tg_xxxzz_xyy_1, tg_xxxzz_xyz_0, tg_xxxzz_xyz_1, \
                                         tg_xxxzz_xzz_0, tg_xxxzz_xzz_1, tg_xxxzz_yyy_0, tg_xxxzz_yyy_1, tg_xxxzz_yyz_0, \
                                         tg_xxxzz_yyz_1, tg_xxxzz_yzz_0, tg_xxxzz_yzz_1, tg_xxxzz_zzz_0, tg_xxxzz_zzz_1, \
                                         tg_xxyyy_xxx_0, tg_xxyyy_xxx_1, tg_xxyyy_xxy_0, tg_xxyyy_xxy_1, tg_xxyyy_xxz_0, \
                                         tg_xxyyy_xxz_1, tg_xxyyy_xyy_0, tg_xxyyy_xyy_1, tg_xxyyy_xyz_0, tg_xxyyy_xyz_1, \
                                         tg_xxyyy_xzz_0, tg_xxyyy_xzz_1, tg_xxyyy_yyy_0, tg_xxyyy_yyy_1, tg_xxyyy_yyz_0, \
                                         tg_xxyyy_yyz_1, tg_xxyyy_yzz_0, tg_xxyyy_yzz_1, tg_xxyyy_zzz_0, tg_xxyyy_zzz_1, \
                                         tg_xxyyz_xxx_0, tg_xxyyz_xxx_1, tg_xxyyz_xxy_0, tg_xxyyz_xxy_1, tg_xxyyz_xxz_0, \
                                         tg_xxyyz_xxz_1, tg_xxyyz_xyy_0, tg_xxyyz_xyy_1, tg_xxyyz_xyz_0, tg_xxyyz_xyz_1, \
                                         tg_xxyyz_xzz_0, tg_xxyyz_xzz_1, tg_xxyyz_yyy_0, tg_xxyyz_yyy_1, tg_xxyyz_yyz_0, \
                                         tg_xxyyz_yyz_1, tg_xxyyz_yzz_0, tg_xxyyz_yzz_1, tg_xxyyz_zzz_0, tg_xxyyz_zzz_1, \
                                         tg_xxyzz_xxx_0, tg_xxyzz_xxx_1, tg_xxyzz_xxy_0, tg_xxyzz_xxy_1, tg_xxyzz_xxz_0, \
                                         tg_xxyzz_xxz_1, tg_xxyzz_xyy_0, tg_xxyzz_xyy_1, tg_xxyzz_xyz_0, tg_xxyzz_xyz_1, \
                                         tg_xxyzz_xzz_0, tg_xxyzz_xzz_1, tg_xxyzz_yyy_0, tg_xxyzz_yyy_1, tg_xxyzz_yyz_0, \
                                         tg_xxyzz_yyz_1, tg_xxyzz_yzz_0, tg_xxyzz_yzz_1, tg_xxyzz_zzz_0, tg_xxyzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxxx_xxx_0[j] = pb_x * tg_xxxxxx_xxx_0[j] + fr * tg_xxxxxx_xxx_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxx_0[j] - tg_xxxxx_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxx_xx_1[j];

                    tg_xxxxxxx_xxy_0[j] = pb_x * tg_xxxxxx_xxy_0[j] + fr * tg_xxxxxx_xxy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxy_0[j] - tg_xxxxx_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xy_1[j];

                    tg_xxxxxxx_xxz_0[j] = pb_x * tg_xxxxxx_xxz_0[j] + fr * tg_xxxxxx_xxz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xxz_0[j] - tg_xxxxx_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxx_xz_1[j];

                    tg_xxxxxxx_xyy_0[j] = pb_x * tg_xxxxxx_xyy_0[j] + fr * tg_xxxxxx_xyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyy_0[j] - tg_xxxxx_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yy_1[j];

                    tg_xxxxxxx_xyz_0[j] = pb_x * tg_xxxxxx_xyz_0[j] + fr * tg_xxxxxx_xyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xyz_0[j] - tg_xxxxx_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_yz_1[j];

                    tg_xxxxxxx_xzz_0[j] = pb_x * tg_xxxxxx_xzz_0[j] + fr * tg_xxxxxx_xzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_xzz_0[j] - tg_xxxxx_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxx_zz_1[j];

                    tg_xxxxxxx_yyy_0[j] = pb_x * tg_xxxxxx_yyy_0[j] + fr * tg_xxxxxx_yyy_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyy_0[j] - tg_xxxxx_yyy_1[j] * fl1_fza);

                    tg_xxxxxxx_yyz_0[j] = pb_x * tg_xxxxxx_yyz_0[j] + fr * tg_xxxxxx_yyz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yyz_0[j] - tg_xxxxx_yyz_1[j] * fl1_fza);

                    tg_xxxxxxx_yzz_0[j] = pb_x * tg_xxxxxx_yzz_0[j] + fr * tg_xxxxxx_yzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_yzz_0[j] - tg_xxxxx_yzz_1[j] * fl1_fza);

                    tg_xxxxxxx_zzz_0[j] = pb_x * tg_xxxxxx_zzz_0[j] + fr * tg_xxxxxx_zzz_1[j] + 3.0 * fl1_fx * (tg_xxxxx_zzz_0[j] - tg_xxxxx_zzz_1[j] * fl1_fza);

                    tg_xxxxxxy_xxx_0[j] = pb_x * tg_xxxxxy_xxx_0[j] + fr * tg_xxxxxy_xxx_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxx_0[j] - tg_xxxxy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxy_xx_1[j];

                    tg_xxxxxxy_xxy_0[j] = pb_x * tg_xxxxxy_xxy_0[j] + fr * tg_xxxxxy_xxy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxy_0[j] - tg_xxxxy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xy_1[j];

                    tg_xxxxxxy_xxz_0[j] = pb_x * tg_xxxxxy_xxz_0[j] + fr * tg_xxxxxy_xxz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xxz_0[j] - tg_xxxxy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxy_xz_1[j];

                    tg_xxxxxxy_xyy_0[j] = pb_x * tg_xxxxxy_xyy_0[j] + fr * tg_xxxxxy_xyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyy_0[j] - tg_xxxxy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yy_1[j];

                    tg_xxxxxxy_xyz_0[j] = pb_x * tg_xxxxxy_xyz_0[j] + fr * tg_xxxxxy_xyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xyz_0[j] - tg_xxxxy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_yz_1[j];

                    tg_xxxxxxy_xzz_0[j] = pb_x * tg_xxxxxy_xzz_0[j] + fr * tg_xxxxxy_xzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_xzz_0[j] - tg_xxxxy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxy_zz_1[j];

                    tg_xxxxxxy_yyy_0[j] = pb_x * tg_xxxxxy_yyy_0[j] + fr * tg_xxxxxy_yyy_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyy_0[j] - tg_xxxxy_yyy_1[j] * fl1_fza);

                    tg_xxxxxxy_yyz_0[j] = pb_x * tg_xxxxxy_yyz_0[j] + fr * tg_xxxxxy_yyz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yyz_0[j] - tg_xxxxy_yyz_1[j] * fl1_fza);

                    tg_xxxxxxy_yzz_0[j] = pb_x * tg_xxxxxy_yzz_0[j] + fr * tg_xxxxxy_yzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_yzz_0[j] - tg_xxxxy_yzz_1[j] * fl1_fza);

                    tg_xxxxxxy_zzz_0[j] = pb_x * tg_xxxxxy_zzz_0[j] + fr * tg_xxxxxy_zzz_1[j] + 2.5 * fl1_fx * (tg_xxxxy_zzz_0[j] - tg_xxxxy_zzz_1[j] * fl1_fza);

                    tg_xxxxxxz_xxx_0[j] = pb_x * tg_xxxxxz_xxx_0[j] + fr * tg_xxxxxz_xxx_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxx_0[j] - tg_xxxxz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxz_xx_1[j];

                    tg_xxxxxxz_xxy_0[j] = pb_x * tg_xxxxxz_xxy_0[j] + fr * tg_xxxxxz_xxy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxy_0[j] - tg_xxxxz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xy_1[j];

                    tg_xxxxxxz_xxz_0[j] = pb_x * tg_xxxxxz_xxz_0[j] + fr * tg_xxxxxz_xxz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xxz_0[j] - tg_xxxxz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxz_xz_1[j];

                    tg_xxxxxxz_xyy_0[j] = pb_x * tg_xxxxxz_xyy_0[j] + fr * tg_xxxxxz_xyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyy_0[j] - tg_xxxxz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yy_1[j];

                    tg_xxxxxxz_xyz_0[j] = pb_x * tg_xxxxxz_xyz_0[j] + fr * tg_xxxxxz_xyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xyz_0[j] - tg_xxxxz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_yz_1[j];

                    tg_xxxxxxz_xzz_0[j] = pb_x * tg_xxxxxz_xzz_0[j] + fr * tg_xxxxxz_xzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_xzz_0[j] - tg_xxxxz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxz_zz_1[j];

                    tg_xxxxxxz_yyy_0[j] = pb_x * tg_xxxxxz_yyy_0[j] + fr * tg_xxxxxz_yyy_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyy_0[j] - tg_xxxxz_yyy_1[j] * fl1_fza);

                    tg_xxxxxxz_yyz_0[j] = pb_x * tg_xxxxxz_yyz_0[j] + fr * tg_xxxxxz_yyz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yyz_0[j] - tg_xxxxz_yyz_1[j] * fl1_fza);

                    tg_xxxxxxz_yzz_0[j] = pb_x * tg_xxxxxz_yzz_0[j] + fr * tg_xxxxxz_yzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_yzz_0[j] - tg_xxxxz_yzz_1[j] * fl1_fza);

                    tg_xxxxxxz_zzz_0[j] = pb_x * tg_xxxxxz_zzz_0[j] + fr * tg_xxxxxz_zzz_1[j] + 2.5 * fl1_fx * (tg_xxxxz_zzz_0[j] - tg_xxxxz_zzz_1[j] * fl1_fza);

                    tg_xxxxxyy_xxx_0[j] = pb_x * tg_xxxxyy_xxx_0[j] + fr * tg_xxxxyy_xxx_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxx_0[j] - tg_xxxyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyy_xx_1[j];

                    tg_xxxxxyy_xxy_0[j] = pb_x * tg_xxxxyy_xxy_0[j] + fr * tg_xxxxyy_xxy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxy_0[j] - tg_xxxyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xy_1[j];

                    tg_xxxxxyy_xxz_0[j] = pb_x * tg_xxxxyy_xxz_0[j] + fr * tg_xxxxyy_xxz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xxz_0[j] - tg_xxxyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyy_xz_1[j];

                    tg_xxxxxyy_xyy_0[j] = pb_x * tg_xxxxyy_xyy_0[j] + fr * tg_xxxxyy_xyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyy_0[j] - tg_xxxyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yy_1[j];

                    tg_xxxxxyy_xyz_0[j] = pb_x * tg_xxxxyy_xyz_0[j] + fr * tg_xxxxyy_xyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xyz_0[j] - tg_xxxyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_yz_1[j];

                    tg_xxxxxyy_xzz_0[j] = pb_x * tg_xxxxyy_xzz_0[j] + fr * tg_xxxxyy_xzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_xzz_0[j] - tg_xxxyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyy_zz_1[j];

                    tg_xxxxxyy_yyy_0[j] = pb_x * tg_xxxxyy_yyy_0[j] + fr * tg_xxxxyy_yyy_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyy_0[j] - tg_xxxyy_yyy_1[j] * fl1_fza);

                    tg_xxxxxyy_yyz_0[j] = pb_x * tg_xxxxyy_yyz_0[j] + fr * tg_xxxxyy_yyz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yyz_0[j] - tg_xxxyy_yyz_1[j] * fl1_fza);

                    tg_xxxxxyy_yzz_0[j] = pb_x * tg_xxxxyy_yzz_0[j] + fr * tg_xxxxyy_yzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_yzz_0[j] - tg_xxxyy_yzz_1[j] * fl1_fza);

                    tg_xxxxxyy_zzz_0[j] = pb_x * tg_xxxxyy_zzz_0[j] + fr * tg_xxxxyy_zzz_1[j] + 2.0 * fl1_fx * (tg_xxxyy_zzz_0[j] - tg_xxxyy_zzz_1[j] * fl1_fza);

                    tg_xxxxxyz_xxx_0[j] = pb_x * tg_xxxxyz_xxx_0[j] + fr * tg_xxxxyz_xxx_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxx_0[j] - tg_xxxyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyz_xx_1[j];

                    tg_xxxxxyz_xxy_0[j] = pb_x * tg_xxxxyz_xxy_0[j] + fr * tg_xxxxyz_xxy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxy_0[j] - tg_xxxyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xy_1[j];

                    tg_xxxxxyz_xxz_0[j] = pb_x * tg_xxxxyz_xxz_0[j] + fr * tg_xxxxyz_xxz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xxz_0[j] - tg_xxxyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyz_xz_1[j];

                    tg_xxxxxyz_xyy_0[j] = pb_x * tg_xxxxyz_xyy_0[j] + fr * tg_xxxxyz_xyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyy_0[j] - tg_xxxyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yy_1[j];

                    tg_xxxxxyz_xyz_0[j] = pb_x * tg_xxxxyz_xyz_0[j] + fr * tg_xxxxyz_xyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xyz_0[j] - tg_xxxyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_yz_1[j];

                    tg_xxxxxyz_xzz_0[j] = pb_x * tg_xxxxyz_xzz_0[j] + fr * tg_xxxxyz_xzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_xzz_0[j] - tg_xxxyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyz_zz_1[j];

                    tg_xxxxxyz_yyy_0[j] = pb_x * tg_xxxxyz_yyy_0[j] + fr * tg_xxxxyz_yyy_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyy_0[j] - tg_xxxyz_yyy_1[j] * fl1_fza);

                    tg_xxxxxyz_yyz_0[j] = pb_x * tg_xxxxyz_yyz_0[j] + fr * tg_xxxxyz_yyz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yyz_0[j] - tg_xxxyz_yyz_1[j] * fl1_fza);

                    tg_xxxxxyz_yzz_0[j] = pb_x * tg_xxxxyz_yzz_0[j] + fr * tg_xxxxyz_yzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_yzz_0[j] - tg_xxxyz_yzz_1[j] * fl1_fza);

                    tg_xxxxxyz_zzz_0[j] = pb_x * tg_xxxxyz_zzz_0[j] + fr * tg_xxxxyz_zzz_1[j] + 2.0 * fl1_fx * (tg_xxxyz_zzz_0[j] - tg_xxxyz_zzz_1[j] * fl1_fza);

                    tg_xxxxxzz_xxx_0[j] = pb_x * tg_xxxxzz_xxx_0[j] + fr * tg_xxxxzz_xxx_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxx_0[j] - tg_xxxzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzz_xx_1[j];

                    tg_xxxxxzz_xxy_0[j] = pb_x * tg_xxxxzz_xxy_0[j] + fr * tg_xxxxzz_xxy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxy_0[j] - tg_xxxzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xy_1[j];

                    tg_xxxxxzz_xxz_0[j] = pb_x * tg_xxxxzz_xxz_0[j] + fr * tg_xxxxzz_xxz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xxz_0[j] - tg_xxxzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzz_xz_1[j];

                    tg_xxxxxzz_xyy_0[j] = pb_x * tg_xxxxzz_xyy_0[j] + fr * tg_xxxxzz_xyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyy_0[j] - tg_xxxzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yy_1[j];

                    tg_xxxxxzz_xyz_0[j] = pb_x * tg_xxxxzz_xyz_0[j] + fr * tg_xxxxzz_xyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xyz_0[j] - tg_xxxzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_yz_1[j];

                    tg_xxxxxzz_xzz_0[j] = pb_x * tg_xxxxzz_xzz_0[j] + fr * tg_xxxxzz_xzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_xzz_0[j] - tg_xxxzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzz_zz_1[j];

                    tg_xxxxxzz_yyy_0[j] = pb_x * tg_xxxxzz_yyy_0[j] + fr * tg_xxxxzz_yyy_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyy_0[j] - tg_xxxzz_yyy_1[j] * fl1_fza);

                    tg_xxxxxzz_yyz_0[j] = pb_x * tg_xxxxzz_yyz_0[j] + fr * tg_xxxxzz_yyz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yyz_0[j] - tg_xxxzz_yyz_1[j] * fl1_fza);

                    tg_xxxxxzz_yzz_0[j] = pb_x * tg_xxxxzz_yzz_0[j] + fr * tg_xxxxzz_yzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_yzz_0[j] - tg_xxxzz_yzz_1[j] * fl1_fza);

                    tg_xxxxxzz_zzz_0[j] = pb_x * tg_xxxxzz_zzz_0[j] + fr * tg_xxxxzz_zzz_1[j] + 2.0 * fl1_fx * (tg_xxxzz_zzz_0[j] - tg_xxxzz_zzz_1[j] * fl1_fza);

                    tg_xxxxyyy_xxx_0[j] = pb_x * tg_xxxyyy_xxx_0[j] + fr * tg_xxxyyy_xxx_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxx_0[j] - tg_xxyyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyy_xx_1[j];

                    tg_xxxxyyy_xxy_0[j] = pb_x * tg_xxxyyy_xxy_0[j] + fr * tg_xxxyyy_xxy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxy_0[j] - tg_xxyyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xy_1[j];

                    tg_xxxxyyy_xxz_0[j] = pb_x * tg_xxxyyy_xxz_0[j] + fr * tg_xxxyyy_xxz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xxz_0[j] - tg_xxyyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyy_xz_1[j];

                    tg_xxxxyyy_xyy_0[j] = pb_x * tg_xxxyyy_xyy_0[j] + fr * tg_xxxyyy_xyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyy_0[j] - tg_xxyyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yy_1[j];

                    tg_xxxxyyy_xyz_0[j] = pb_x * tg_xxxyyy_xyz_0[j] + fr * tg_xxxyyy_xyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xyz_0[j] - tg_xxyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_yz_1[j];

                    tg_xxxxyyy_xzz_0[j] = pb_x * tg_xxxyyy_xzz_0[j] + fr * tg_xxxyyy_xzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_xzz_0[j] - tg_xxyyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyy_zz_1[j];

                    tg_xxxxyyy_yyy_0[j] = pb_x * tg_xxxyyy_yyy_0[j] + fr * tg_xxxyyy_yyy_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyy_0[j] - tg_xxyyy_yyy_1[j] * fl1_fza);

                    tg_xxxxyyy_yyz_0[j] = pb_x * tg_xxxyyy_yyz_0[j] + fr * tg_xxxyyy_yyz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yyz_0[j] - tg_xxyyy_yyz_1[j] * fl1_fza);

                    tg_xxxxyyy_yzz_0[j] = pb_x * tg_xxxyyy_yzz_0[j] + fr * tg_xxxyyy_yzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_yzz_0[j] - tg_xxyyy_yzz_1[j] * fl1_fza);

                    tg_xxxxyyy_zzz_0[j] = pb_x * tg_xxxyyy_zzz_0[j] + fr * tg_xxxyyy_zzz_1[j] + 1.5 * fl1_fx * (tg_xxyyy_zzz_0[j] - tg_xxyyy_zzz_1[j] * fl1_fza);

                    tg_xxxxyyz_xxx_0[j] = pb_x * tg_xxxyyz_xxx_0[j] + fr * tg_xxxyyz_xxx_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxx_0[j] - tg_xxyyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyz_xx_1[j];

                    tg_xxxxyyz_xxy_0[j] = pb_x * tg_xxxyyz_xxy_0[j] + fr * tg_xxxyyz_xxy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxy_0[j] - tg_xxyyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xy_1[j];

                    tg_xxxxyyz_xxz_0[j] = pb_x * tg_xxxyyz_xxz_0[j] + fr * tg_xxxyyz_xxz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xxz_0[j] - tg_xxyyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyz_xz_1[j];

                    tg_xxxxyyz_xyy_0[j] = pb_x * tg_xxxyyz_xyy_0[j] + fr * tg_xxxyyz_xyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyy_0[j] - tg_xxyyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yy_1[j];

                    tg_xxxxyyz_xyz_0[j] = pb_x * tg_xxxyyz_xyz_0[j] + fr * tg_xxxyyz_xyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xyz_0[j] - tg_xxyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_yz_1[j];

                    tg_xxxxyyz_xzz_0[j] = pb_x * tg_xxxyyz_xzz_0[j] + fr * tg_xxxyyz_xzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_xzz_0[j] - tg_xxyyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyz_zz_1[j];

                    tg_xxxxyyz_yyy_0[j] = pb_x * tg_xxxyyz_yyy_0[j] + fr * tg_xxxyyz_yyy_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyy_0[j] - tg_xxyyz_yyy_1[j] * fl1_fza);

                    tg_xxxxyyz_yyz_0[j] = pb_x * tg_xxxyyz_yyz_0[j] + fr * tg_xxxyyz_yyz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yyz_0[j] - tg_xxyyz_yyz_1[j] * fl1_fza);

                    tg_xxxxyyz_yzz_0[j] = pb_x * tg_xxxyyz_yzz_0[j] + fr * tg_xxxyyz_yzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_yzz_0[j] - tg_xxyyz_yzz_1[j] * fl1_fza);

                    tg_xxxxyyz_zzz_0[j] = pb_x * tg_xxxyyz_zzz_0[j] + fr * tg_xxxyyz_zzz_1[j] + 1.5 * fl1_fx * (tg_xxyyz_zzz_0[j] - tg_xxyyz_zzz_1[j] * fl1_fza);

                    tg_xxxxyzz_xxx_0[j] = pb_x * tg_xxxyzz_xxx_0[j] + fr * tg_xxxyzz_xxx_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxx_0[j] - tg_xxyzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzz_xx_1[j];

                    tg_xxxxyzz_xxy_0[j] = pb_x * tg_xxxyzz_xxy_0[j] + fr * tg_xxxyzz_xxy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxy_0[j] - tg_xxyzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xy_1[j];

                    tg_xxxxyzz_xxz_0[j] = pb_x * tg_xxxyzz_xxz_0[j] + fr * tg_xxxyzz_xxz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xxz_0[j] - tg_xxyzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzz_xz_1[j];

                    tg_xxxxyzz_xyy_0[j] = pb_x * tg_xxxyzz_xyy_0[j] + fr * tg_xxxyzz_xyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyy_0[j] - tg_xxyzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yy_1[j];

                    tg_xxxxyzz_xyz_0[j] = pb_x * tg_xxxyzz_xyz_0[j] + fr * tg_xxxyzz_xyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xyz_0[j] - tg_xxyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_yz_1[j];

                    tg_xxxxyzz_xzz_0[j] = pb_x * tg_xxxyzz_xzz_0[j] + fr * tg_xxxyzz_xzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_xzz_0[j] - tg_xxyzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzz_zz_1[j];

                    tg_xxxxyzz_yyy_0[j] = pb_x * tg_xxxyzz_yyy_0[j] + fr * tg_xxxyzz_yyy_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyy_0[j] - tg_xxyzz_yyy_1[j] * fl1_fza);

                    tg_xxxxyzz_yyz_0[j] = pb_x * tg_xxxyzz_yyz_0[j] + fr * tg_xxxyzz_yyz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yyz_0[j] - tg_xxyzz_yyz_1[j] * fl1_fza);

                    tg_xxxxyzz_yzz_0[j] = pb_x * tg_xxxyzz_yzz_0[j] + fr * tg_xxxyzz_yzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_yzz_0[j] - tg_xxyzz_yzz_1[j] * fl1_fza);

                    tg_xxxxyzz_zzz_0[j] = pb_x * tg_xxxyzz_zzz_0[j] + fr * tg_xxxyzz_zzz_1[j] + 1.5 * fl1_fx * (tg_xxyzz_zzz_0[j] - tg_xxyzz_zzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSF_90_180(      CMemBlock2D<double>* primBuffer,
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
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xxxzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 90); 

                auto tg_xxxzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 91); 

                auto tg_xxxzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 92); 

                auto tg_xxxzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 93); 

                auto tg_xxxzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 94); 

                auto tg_xxxzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 95); 

                auto tg_xxxzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 96); 

                auto tg_xxxzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 97); 

                auto tg_xxxzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 98); 

                auto tg_xxxzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 99); 

                auto tg_xxyyyy_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 100); 

                auto tg_xxyyyy_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 101); 

                auto tg_xxyyyy_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 102); 

                auto tg_xxyyyy_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 103); 

                auto tg_xxyyyy_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 104); 

                auto tg_xxyyyy_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 105); 

                auto tg_xxyyyy_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 106); 

                auto tg_xxyyyy_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 107); 

                auto tg_xxyyyy_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 108); 

                auto tg_xxyyyy_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 109); 

                auto tg_xxyyyz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 110); 

                auto tg_xxyyyz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 111); 

                auto tg_xxyyyz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 112); 

                auto tg_xxyyyz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 113); 

                auto tg_xxyyyz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 114); 

                auto tg_xxyyyz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 115); 

                auto tg_xxyyyz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 116); 

                auto tg_xxyyyz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 117); 

                auto tg_xxyyyz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 118); 

                auto tg_xxyyyz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 119); 

                auto tg_xxyyzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 120); 

                auto tg_xxyyzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 121); 

                auto tg_xxyyzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 122); 

                auto tg_xxyyzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 123); 

                auto tg_xxyyzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 124); 

                auto tg_xxyyzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 125); 

                auto tg_xxyyzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 126); 

                auto tg_xxyyzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 127); 

                auto tg_xxyyzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 128); 

                auto tg_xxyyzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 129); 

                auto tg_xxyzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 130); 

                auto tg_xxyzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 131); 

                auto tg_xxyzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 132); 

                auto tg_xxyzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 133); 

                auto tg_xxyzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 134); 

                auto tg_xxyzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 135); 

                auto tg_xxyzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 136); 

                auto tg_xxyzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 137); 

                auto tg_xxyzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 138); 

                auto tg_xxyzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 139); 

                auto tg_xxzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 140); 

                auto tg_xxzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 141); 

                auto tg_xxzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 142); 

                auto tg_xxzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 143); 

                auto tg_xxzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 144); 

                auto tg_xxzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 145); 

                auto tg_xxzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 146); 

                auto tg_xxzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 147); 

                auto tg_xxzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 148); 

                auto tg_xxzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 149); 

                auto tg_xyyyyy_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 150); 

                auto tg_xyyyyy_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 151); 

                auto tg_xyyyyy_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 152); 

                auto tg_xyyyyy_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 153); 

                auto tg_xyyyyy_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 154); 

                auto tg_xyyyyy_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 155); 

                auto tg_xyyyyy_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 156); 

                auto tg_xyyyyy_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 157); 

                auto tg_xyyyyy_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 158); 

                auto tg_xyyyyy_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 159); 

                auto tg_xyyyyz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 160); 

                auto tg_xyyyyz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 161); 

                auto tg_xyyyyz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 162); 

                auto tg_xyyyyz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 163); 

                auto tg_xyyyyz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 164); 

                auto tg_xyyyyz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 165); 

                auto tg_xyyyyz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 166); 

                auto tg_xyyyyz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 167); 

                auto tg_xyyyyz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 168); 

                auto tg_xyyyyz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 169); 

                auto tg_xyyyzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 170); 

                auto tg_xyyyzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 171); 

                auto tg_xyyyzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 172); 

                auto tg_xyyyzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 173); 

                auto tg_xyyyzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 174); 

                auto tg_xyyyzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 175); 

                auto tg_xyyyzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 176); 

                auto tg_xyyyzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 177); 

                auto tg_xyyyzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 178); 

                auto tg_xyyyzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 179); 

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

                auto tg_xxzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 90); 

                auto tg_xxzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 91); 

                auto tg_xxzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 92); 

                auto tg_xxzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 93); 

                auto tg_xxzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 94); 

                auto tg_xxzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 95); 

                auto tg_xxzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 96); 

                auto tg_xxzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 97); 

                auto tg_xxzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 98); 

                auto tg_xxzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 99); 

                auto tg_xyyyy_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 100); 

                auto tg_xyyyy_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 101); 

                auto tg_xyyyy_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 102); 

                auto tg_xyyyy_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 103); 

                auto tg_xyyyy_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 104); 

                auto tg_xyyyy_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 105); 

                auto tg_xyyyy_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 106); 

                auto tg_xyyyy_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 107); 

                auto tg_xyyyy_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 108); 

                auto tg_xyyyy_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 109); 

                auto tg_xyyyz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 110); 

                auto tg_xyyyz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 111); 

                auto tg_xyyyz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 112); 

                auto tg_xyyyz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 113); 

                auto tg_xyyyz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 114); 

                auto tg_xyyyz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 115); 

                auto tg_xyyyz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 116); 

                auto tg_xyyyz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 117); 

                auto tg_xyyyz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 118); 

                auto tg_xyyyz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 119); 

                auto tg_xyyzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 120); 

                auto tg_xyyzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 121); 

                auto tg_xyyzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 122); 

                auto tg_xyyzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 123); 

                auto tg_xyyzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 124); 

                auto tg_xyyzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 125); 

                auto tg_xyyzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 126); 

                auto tg_xyyzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 127); 

                auto tg_xyyzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 128); 

                auto tg_xyyzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 129); 

                auto tg_xyzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 130); 

                auto tg_xyzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 131); 

                auto tg_xyzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 132); 

                auto tg_xyzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 133); 

                auto tg_xyzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 134); 

                auto tg_xyzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 135); 

                auto tg_xyzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 136); 

                auto tg_xyzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 137); 

                auto tg_xyzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 138); 

                auto tg_xyzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 139); 

                auto tg_xzzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 140); 

                auto tg_xzzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 141); 

                auto tg_xzzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 142); 

                auto tg_xzzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 143); 

                auto tg_xzzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 144); 

                auto tg_xzzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 145); 

                auto tg_xzzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 146); 

                auto tg_xzzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 147); 

                auto tg_xzzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 148); 

                auto tg_xzzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 149); 

                auto tg_yyyyy_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 150); 

                auto tg_yyyyy_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 151); 

                auto tg_yyyyy_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 152); 

                auto tg_yyyyy_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 153); 

                auto tg_yyyyy_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 154); 

                auto tg_yyyyy_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 155); 

                auto tg_yyyyy_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 156); 

                auto tg_yyyyy_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 157); 

                auto tg_yyyyy_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 158); 

                auto tg_yyyyy_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 159); 

                auto tg_yyyyz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 160); 

                auto tg_yyyyz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 161); 

                auto tg_yyyyz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 162); 

                auto tg_yyyyz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 163); 

                auto tg_yyyyz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 164); 

                auto tg_yyyyz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 165); 

                auto tg_yyyyz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 166); 

                auto tg_yyyyz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 167); 

                auto tg_yyyyz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 168); 

                auto tg_yyyyz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 169); 

                auto tg_yyyzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 170); 

                auto tg_yyyzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 171); 

                auto tg_yyyzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 172); 

                auto tg_yyyzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 173); 

                auto tg_yyyzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 174); 

                auto tg_yyyzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 175); 

                auto tg_yyyzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 176); 

                auto tg_yyyzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 177); 

                auto tg_yyyzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 178); 

                auto tg_yyyzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 179); 

                auto tg_xxzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 90); 

                auto tg_xxzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 91); 

                auto tg_xxzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 92); 

                auto tg_xxzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 93); 

                auto tg_xxzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 94); 

                auto tg_xxzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 95); 

                auto tg_xxzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 96); 

                auto tg_xxzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 97); 

                auto tg_xxzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 98); 

                auto tg_xxzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 99); 

                auto tg_xyyyy_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 100); 

                auto tg_xyyyy_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 101); 

                auto tg_xyyyy_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 102); 

                auto tg_xyyyy_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 103); 

                auto tg_xyyyy_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 104); 

                auto tg_xyyyy_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 105); 

                auto tg_xyyyy_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 106); 

                auto tg_xyyyy_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 107); 

                auto tg_xyyyy_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 108); 

                auto tg_xyyyy_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 109); 

                auto tg_xyyyz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 110); 

                auto tg_xyyyz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 111); 

                auto tg_xyyyz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 112); 

                auto tg_xyyyz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 113); 

                auto tg_xyyyz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 114); 

                auto tg_xyyyz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 115); 

                auto tg_xyyyz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 116); 

                auto tg_xyyyz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 117); 

                auto tg_xyyyz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 118); 

                auto tg_xyyyz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 119); 

                auto tg_xyyzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 120); 

                auto tg_xyyzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 121); 

                auto tg_xyyzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 122); 

                auto tg_xyyzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 123); 

                auto tg_xyyzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 124); 

                auto tg_xyyzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 125); 

                auto tg_xyyzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 126); 

                auto tg_xyyzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 127); 

                auto tg_xyyzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 128); 

                auto tg_xyyzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 129); 

                auto tg_xyzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 130); 

                auto tg_xyzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 131); 

                auto tg_xyzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 132); 

                auto tg_xyzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 133); 

                auto tg_xyzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 134); 

                auto tg_xyzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 135); 

                auto tg_xyzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 136); 

                auto tg_xyzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 137); 

                auto tg_xyzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 138); 

                auto tg_xyzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 139); 

                auto tg_xzzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 140); 

                auto tg_xzzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 141); 

                auto tg_xzzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 142); 

                auto tg_xzzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 143); 

                auto tg_xzzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 144); 

                auto tg_xzzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 145); 

                auto tg_xzzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 146); 

                auto tg_xzzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 147); 

                auto tg_xzzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 148); 

                auto tg_xzzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 149); 

                auto tg_yyyyy_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 150); 

                auto tg_yyyyy_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 151); 

                auto tg_yyyyy_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 152); 

                auto tg_yyyyy_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 153); 

                auto tg_yyyyy_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 154); 

                auto tg_yyyyy_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 155); 

                auto tg_yyyyy_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 156); 

                auto tg_yyyyy_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 157); 

                auto tg_yyyyy_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 158); 

                auto tg_yyyyy_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 159); 

                auto tg_yyyyz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 160); 

                auto tg_yyyyz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 161); 

                auto tg_yyyyz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 162); 

                auto tg_yyyyz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 163); 

                auto tg_yyyyz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 164); 

                auto tg_yyyyz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 165); 

                auto tg_yyyyz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 166); 

                auto tg_yyyyz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 167); 

                auto tg_yyyyz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 168); 

                auto tg_yyyyz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 169); 

                auto tg_yyyzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 170); 

                auto tg_yyyzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 171); 

                auto tg_yyyzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 172); 

                auto tg_yyyzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 173); 

                auto tg_yyyzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 174); 

                auto tg_yyyzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 175); 

                auto tg_yyyzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 176); 

                auto tg_yyyzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 177); 

                auto tg_yyyzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 178); 

                auto tg_yyyzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 179); 

                auto tg_xxxzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 54); 

                auto tg_xxxzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 55); 

                auto tg_xxxzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 56); 

                auto tg_xxxzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 57); 

                auto tg_xxxzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 58); 

                auto tg_xxxzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 59); 

                auto tg_xxyyyy_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 60); 

                auto tg_xxyyyy_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 61); 

                auto tg_xxyyyy_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 62); 

                auto tg_xxyyyy_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 63); 

                auto tg_xxyyyy_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 64); 

                auto tg_xxyyyy_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 65); 

                auto tg_xxyyyz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 66); 

                auto tg_xxyyyz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 67); 

                auto tg_xxyyyz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 68); 

                auto tg_xxyyyz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 69); 

                auto tg_xxyyyz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 70); 

                auto tg_xxyyyz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 71); 

                auto tg_xxyyzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 72); 

                auto tg_xxyyzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 73); 

                auto tg_xxyyzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 74); 

                auto tg_xxyyzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 75); 

                auto tg_xxyyzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 76); 

                auto tg_xxyyzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 77); 

                auto tg_xxyzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 78); 

                auto tg_xxyzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 79); 

                auto tg_xxyzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 80); 

                auto tg_xxyzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 81); 

                auto tg_xxyzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 82); 

                auto tg_xxyzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 83); 

                auto tg_xxzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 84); 

                auto tg_xxzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 85); 

                auto tg_xxzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 86); 

                auto tg_xxzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 87); 

                auto tg_xxzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 88); 

                auto tg_xxzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 89); 

                auto tg_xyyyyy_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 90); 

                auto tg_xyyyyy_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 91); 

                auto tg_xyyyyy_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 92); 

                auto tg_xyyyyy_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 93); 

                auto tg_xyyyyy_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 94); 

                auto tg_xyyyyy_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 95); 

                auto tg_xyyyyz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 96); 

                auto tg_xyyyyz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 97); 

                auto tg_xyyyyz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 98); 

                auto tg_xyyyyz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 99); 

                auto tg_xyyyyz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 100); 

                auto tg_xyyyyz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 101); 

                auto tg_xyyyzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 102); 

                auto tg_xyyyzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 103); 

                auto tg_xyyyzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 104); 

                auto tg_xyyyzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 105); 

                auto tg_xyyyzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 106); 

                auto tg_xyyyzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 107); 

                // set up pointers to integrals

                auto tg_xxxxzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 90); 

                auto tg_xxxxzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 91); 

                auto tg_xxxxzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 92); 

                auto tg_xxxxzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 93); 

                auto tg_xxxxzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 94); 

                auto tg_xxxxzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 95); 

                auto tg_xxxxzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 96); 

                auto tg_xxxxzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 97); 

                auto tg_xxxxzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 98); 

                auto tg_xxxxzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 99); 

                auto tg_xxxyyyy_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 100); 

                auto tg_xxxyyyy_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 101); 

                auto tg_xxxyyyy_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 102); 

                auto tg_xxxyyyy_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 103); 

                auto tg_xxxyyyy_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 104); 

                auto tg_xxxyyyy_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 105); 

                auto tg_xxxyyyy_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 106); 

                auto tg_xxxyyyy_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 107); 

                auto tg_xxxyyyy_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 108); 

                auto tg_xxxyyyy_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 109); 

                auto tg_xxxyyyz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 110); 

                auto tg_xxxyyyz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 111); 

                auto tg_xxxyyyz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 112); 

                auto tg_xxxyyyz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 113); 

                auto tg_xxxyyyz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 114); 

                auto tg_xxxyyyz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 115); 

                auto tg_xxxyyyz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 116); 

                auto tg_xxxyyyz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 117); 

                auto tg_xxxyyyz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 118); 

                auto tg_xxxyyyz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 119); 

                auto tg_xxxyyzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 120); 

                auto tg_xxxyyzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 121); 

                auto tg_xxxyyzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 122); 

                auto tg_xxxyyzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 123); 

                auto tg_xxxyyzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 124); 

                auto tg_xxxyyzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 125); 

                auto tg_xxxyyzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 126); 

                auto tg_xxxyyzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 127); 

                auto tg_xxxyyzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 128); 

                auto tg_xxxyyzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 129); 

                auto tg_xxxyzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 130); 

                auto tg_xxxyzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 131); 

                auto tg_xxxyzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 132); 

                auto tg_xxxyzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 133); 

                auto tg_xxxyzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 134); 

                auto tg_xxxyzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 135); 

                auto tg_xxxyzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 136); 

                auto tg_xxxyzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 137); 

                auto tg_xxxyzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 138); 

                auto tg_xxxyzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 139); 

                auto tg_xxxzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 140); 

                auto tg_xxxzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 141); 

                auto tg_xxxzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 142); 

                auto tg_xxxzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 143); 

                auto tg_xxxzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 144); 

                auto tg_xxxzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 145); 

                auto tg_xxxzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 146); 

                auto tg_xxxzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 147); 

                auto tg_xxxzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 148); 

                auto tg_xxxzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 149); 

                auto tg_xxyyyyy_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 150); 

                auto tg_xxyyyyy_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 151); 

                auto tg_xxyyyyy_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 152); 

                auto tg_xxyyyyy_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 153); 

                auto tg_xxyyyyy_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 154); 

                auto tg_xxyyyyy_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 155); 

                auto tg_xxyyyyy_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 156); 

                auto tg_xxyyyyy_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 157); 

                auto tg_xxyyyyy_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 158); 

                auto tg_xxyyyyy_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 159); 

                auto tg_xxyyyyz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 160); 

                auto tg_xxyyyyz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 161); 

                auto tg_xxyyyyz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 162); 

                auto tg_xxyyyyz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 163); 

                auto tg_xxyyyyz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 164); 

                auto tg_xxyyyyz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 165); 

                auto tg_xxyyyyz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 166); 

                auto tg_xxyyyyz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 167); 

                auto tg_xxyyyyz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 168); 

                auto tg_xxyyyyz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 169); 

                auto tg_xxyyyzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 170); 

                auto tg_xxyyyzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 171); 

                auto tg_xxyyyzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 172); 

                auto tg_xxyyyzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 173); 

                auto tg_xxyyyzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 174); 

                auto tg_xxyyyzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 175); 

                auto tg_xxyyyzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 176); 

                auto tg_xxyyyzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 177); 

                auto tg_xxyyyzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 178); 

                auto tg_xxyyyzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 179); 

                // Batch of Integrals (90,180)

                #pragma omp simd aligned(fxn, fza, tg_xxxxzzz_xxx_0, tg_xxxxzzz_xxy_0, tg_xxxxzzz_xxz_0, \
                                         tg_xxxxzzz_xyy_0, tg_xxxxzzz_xyz_0, tg_xxxxzzz_xzz_0, tg_xxxxzzz_yyy_0, \
                                         tg_xxxxzzz_yyz_0, tg_xxxxzzz_yzz_0, tg_xxxxzzz_zzz_0, tg_xxxyyyy_xxx_0, \
                                         tg_xxxyyyy_xxy_0, tg_xxxyyyy_xxz_0, tg_xxxyyyy_xyy_0, tg_xxxyyyy_xyz_0, \
                                         tg_xxxyyyy_xzz_0, tg_xxxyyyy_yyy_0, tg_xxxyyyy_yyz_0, tg_xxxyyyy_yzz_0, \
                                         tg_xxxyyyy_zzz_0, tg_xxxyyyz_xxx_0, tg_xxxyyyz_xxy_0, tg_xxxyyyz_xxz_0, \
                                         tg_xxxyyyz_xyy_0, tg_xxxyyyz_xyz_0, tg_xxxyyyz_xzz_0, tg_xxxyyyz_yyy_0, \
                                         tg_xxxyyyz_yyz_0, tg_xxxyyyz_yzz_0, tg_xxxyyyz_zzz_0, tg_xxxyyzz_xxx_0, \
                                         tg_xxxyyzz_xxy_0, tg_xxxyyzz_xxz_0, tg_xxxyyzz_xyy_0, tg_xxxyyzz_xyz_0, \
                                         tg_xxxyyzz_xzz_0, tg_xxxyyzz_yyy_0, tg_xxxyyzz_yyz_0, tg_xxxyyzz_yzz_0, \
                                         tg_xxxyyzz_zzz_0, tg_xxxyzzz_xxx_0, tg_xxxyzzz_xxy_0, tg_xxxyzzz_xxz_0, \
                                         tg_xxxyzzz_xyy_0, tg_xxxyzzz_xyz_0, tg_xxxyzzz_xzz_0, tg_xxxyzzz_yyy_0, \
                                         tg_xxxyzzz_yyz_0, tg_xxxyzzz_yzz_0, tg_xxxyzzz_zzz_0, tg_xxxzzz_xx_1, \
                                         tg_xxxzzz_xxx_0, tg_xxxzzz_xxx_1, tg_xxxzzz_xxy_0, tg_xxxzzz_xxy_1, tg_xxxzzz_xxz_0, \
                                         tg_xxxzzz_xxz_1, tg_xxxzzz_xy_1, tg_xxxzzz_xyy_0, tg_xxxzzz_xyy_1, tg_xxxzzz_xyz_0, \
                                         tg_xxxzzz_xyz_1, tg_xxxzzz_xz_1, tg_xxxzzz_xzz_0, tg_xxxzzz_xzz_1, tg_xxxzzz_yy_1, \
                                         tg_xxxzzz_yyy_0, tg_xxxzzz_yyy_1, tg_xxxzzz_yyz_0, tg_xxxzzz_yyz_1, tg_xxxzzz_yz_1, \
                                         tg_xxxzzz_yzz_0, tg_xxxzzz_yzz_1, tg_xxxzzz_zz_1, tg_xxxzzz_zzz_0, tg_xxxzzz_zzz_1, \
                                         tg_xxxzzzz_xxx_0, tg_xxxzzzz_xxy_0, tg_xxxzzzz_xxz_0, tg_xxxzzzz_xyy_0, \
                                         tg_xxxzzzz_xyz_0, tg_xxxzzzz_xzz_0, tg_xxxzzzz_yyy_0, tg_xxxzzzz_yyz_0, \
                                         tg_xxxzzzz_yzz_0, tg_xxxzzzz_zzz_0, tg_xxyyyy_xx_1, tg_xxyyyy_xxx_0, tg_xxyyyy_xxx_1, \
                                         tg_xxyyyy_xxy_0, tg_xxyyyy_xxy_1, tg_xxyyyy_xxz_0, tg_xxyyyy_xxz_1, tg_xxyyyy_xy_1, \
                                         tg_xxyyyy_xyy_0, tg_xxyyyy_xyy_1, tg_xxyyyy_xyz_0, tg_xxyyyy_xyz_1, tg_xxyyyy_xz_1, \
                                         tg_xxyyyy_xzz_0, tg_xxyyyy_xzz_1, tg_xxyyyy_yy_1, tg_xxyyyy_yyy_0, tg_xxyyyy_yyy_1, \
                                         tg_xxyyyy_yyz_0, tg_xxyyyy_yyz_1, tg_xxyyyy_yz_1, tg_xxyyyy_yzz_0, tg_xxyyyy_yzz_1, \
                                         tg_xxyyyy_zz_1, tg_xxyyyy_zzz_0, tg_xxyyyy_zzz_1, tg_xxyyyyy_xxx_0, \
                                         tg_xxyyyyy_xxy_0, tg_xxyyyyy_xxz_0, tg_xxyyyyy_xyy_0, tg_xxyyyyy_xyz_0, \
                                         tg_xxyyyyy_xzz_0, tg_xxyyyyy_yyy_0, tg_xxyyyyy_yyz_0, tg_xxyyyyy_yzz_0, \
                                         tg_xxyyyyy_zzz_0, tg_xxyyyyz_xxx_0, tg_xxyyyyz_xxy_0, tg_xxyyyyz_xxz_0, \
                                         tg_xxyyyyz_xyy_0, tg_xxyyyyz_xyz_0, tg_xxyyyyz_xzz_0, tg_xxyyyyz_yyy_0, \
                                         tg_xxyyyyz_yyz_0, tg_xxyyyyz_yzz_0, tg_xxyyyyz_zzz_0, tg_xxyyyz_xx_1, \
                                         tg_xxyyyz_xxx_0, tg_xxyyyz_xxx_1, tg_xxyyyz_xxy_0, tg_xxyyyz_xxy_1, tg_xxyyyz_xxz_0, \
                                         tg_xxyyyz_xxz_1, tg_xxyyyz_xy_1, tg_xxyyyz_xyy_0, tg_xxyyyz_xyy_1, tg_xxyyyz_xyz_0, \
                                         tg_xxyyyz_xyz_1, tg_xxyyyz_xz_1, tg_xxyyyz_xzz_0, tg_xxyyyz_xzz_1, tg_xxyyyz_yy_1, \
                                         tg_xxyyyz_yyy_0, tg_xxyyyz_yyy_1, tg_xxyyyz_yyz_0, tg_xxyyyz_yyz_1, tg_xxyyyz_yz_1, \
                                         tg_xxyyyz_yzz_0, tg_xxyyyz_yzz_1, tg_xxyyyz_zz_1, tg_xxyyyz_zzz_0, tg_xxyyyz_zzz_1, \
                                         tg_xxyyyzz_xxx_0, tg_xxyyyzz_xxy_0, tg_xxyyyzz_xxz_0, tg_xxyyyzz_xyy_0, \
                                         tg_xxyyyzz_xyz_0, tg_xxyyyzz_xzz_0, tg_xxyyyzz_yyy_0, tg_xxyyyzz_yyz_0, \
                                         tg_xxyyyzz_yzz_0, tg_xxyyyzz_zzz_0, tg_xxyyzz_xx_1, tg_xxyyzz_xxx_0, tg_xxyyzz_xxx_1, \
                                         tg_xxyyzz_xxy_0, tg_xxyyzz_xxy_1, tg_xxyyzz_xxz_0, tg_xxyyzz_xxz_1, tg_xxyyzz_xy_1, \
                                         tg_xxyyzz_xyy_0, tg_xxyyzz_xyy_1, tg_xxyyzz_xyz_0, tg_xxyyzz_xyz_1, tg_xxyyzz_xz_1, \
                                         tg_xxyyzz_xzz_0, tg_xxyyzz_xzz_1, tg_xxyyzz_yy_1, tg_xxyyzz_yyy_0, tg_xxyyzz_yyy_1, \
                                         tg_xxyyzz_yyz_0, tg_xxyyzz_yyz_1, tg_xxyyzz_yz_1, tg_xxyyzz_yzz_0, tg_xxyyzz_yzz_1, \
                                         tg_xxyyzz_zz_1, tg_xxyyzz_zzz_0, tg_xxyyzz_zzz_1, tg_xxyzzz_xx_1, tg_xxyzzz_xxx_0, \
                                         tg_xxyzzz_xxx_1, tg_xxyzzz_xxy_0, tg_xxyzzz_xxy_1, tg_xxyzzz_xxz_0, tg_xxyzzz_xxz_1, \
                                         tg_xxyzzz_xy_1, tg_xxyzzz_xyy_0, tg_xxyzzz_xyy_1, tg_xxyzzz_xyz_0, tg_xxyzzz_xyz_1, \
                                         tg_xxyzzz_xz_1, tg_xxyzzz_xzz_0, tg_xxyzzz_xzz_1, tg_xxyzzz_yy_1, tg_xxyzzz_yyy_0, \
                                         tg_xxyzzz_yyy_1, tg_xxyzzz_yyz_0, tg_xxyzzz_yyz_1, tg_xxyzzz_yz_1, tg_xxyzzz_yzz_0, \
                                         tg_xxyzzz_yzz_1, tg_xxyzzz_zz_1, tg_xxyzzz_zzz_0, tg_xxyzzz_zzz_1, tg_xxzzz_xxx_0, \
                                         tg_xxzzz_xxx_1, tg_xxzzz_xxy_0, tg_xxzzz_xxy_1, tg_xxzzz_xxz_0, tg_xxzzz_xxz_1, \
                                         tg_xxzzz_xyy_0, tg_xxzzz_xyy_1, tg_xxzzz_xyz_0, tg_xxzzz_xyz_1, tg_xxzzz_xzz_0, \
                                         tg_xxzzz_xzz_1, tg_xxzzz_yyy_0, tg_xxzzz_yyy_1, tg_xxzzz_yyz_0, tg_xxzzz_yyz_1, \
                                         tg_xxzzz_yzz_0, tg_xxzzz_yzz_1, tg_xxzzz_zzz_0, tg_xxzzz_zzz_1, tg_xxzzzz_xx_1, \
                                         tg_xxzzzz_xxx_0, tg_xxzzzz_xxx_1, tg_xxzzzz_xxy_0, tg_xxzzzz_xxy_1, tg_xxzzzz_xxz_0, \
                                         tg_xxzzzz_xxz_1, tg_xxzzzz_xy_1, tg_xxzzzz_xyy_0, tg_xxzzzz_xyy_1, tg_xxzzzz_xyz_0, \
                                         tg_xxzzzz_xyz_1, tg_xxzzzz_xz_1, tg_xxzzzz_xzz_0, tg_xxzzzz_xzz_1, tg_xxzzzz_yy_1, \
                                         tg_xxzzzz_yyy_0, tg_xxzzzz_yyy_1, tg_xxzzzz_yyz_0, tg_xxzzzz_yyz_1, tg_xxzzzz_yz_1, \
                                         tg_xxzzzz_yzz_0, tg_xxzzzz_yzz_1, tg_xxzzzz_zz_1, tg_xxzzzz_zzz_0, tg_xxzzzz_zzz_1, \
                                         tg_xyyyy_xxx_0, tg_xyyyy_xxx_1, tg_xyyyy_xxy_0, tg_xyyyy_xxy_1, tg_xyyyy_xxz_0, \
                                         tg_xyyyy_xxz_1, tg_xyyyy_xyy_0, tg_xyyyy_xyy_1, tg_xyyyy_xyz_0, tg_xyyyy_xyz_1, \
                                         tg_xyyyy_xzz_0, tg_xyyyy_xzz_1, tg_xyyyy_yyy_0, tg_xyyyy_yyy_1, tg_xyyyy_yyz_0, \
                                         tg_xyyyy_yyz_1, tg_xyyyy_yzz_0, tg_xyyyy_yzz_1, tg_xyyyy_zzz_0, tg_xyyyy_zzz_1, \
                                         tg_xyyyyy_xx_1, tg_xyyyyy_xxx_0, tg_xyyyyy_xxx_1, tg_xyyyyy_xxy_0, tg_xyyyyy_xxy_1, \
                                         tg_xyyyyy_xxz_0, tg_xyyyyy_xxz_1, tg_xyyyyy_xy_1, tg_xyyyyy_xyy_0, tg_xyyyyy_xyy_1, \
                                         tg_xyyyyy_xyz_0, tg_xyyyyy_xyz_1, tg_xyyyyy_xz_1, tg_xyyyyy_xzz_0, tg_xyyyyy_xzz_1, \
                                         tg_xyyyyy_yy_1, tg_xyyyyy_yyy_0, tg_xyyyyy_yyy_1, tg_xyyyyy_yyz_0, tg_xyyyyy_yyz_1, \
                                         tg_xyyyyy_yz_1, tg_xyyyyy_yzz_0, tg_xyyyyy_yzz_1, tg_xyyyyy_zz_1, tg_xyyyyy_zzz_0, \
                                         tg_xyyyyy_zzz_1, tg_xyyyyz_xx_1, tg_xyyyyz_xxx_0, tg_xyyyyz_xxx_1, tg_xyyyyz_xxy_0, \
                                         tg_xyyyyz_xxy_1, tg_xyyyyz_xxz_0, tg_xyyyyz_xxz_1, tg_xyyyyz_xy_1, tg_xyyyyz_xyy_0, \
                                         tg_xyyyyz_xyy_1, tg_xyyyyz_xyz_0, tg_xyyyyz_xyz_1, tg_xyyyyz_xz_1, tg_xyyyyz_xzz_0, \
                                         tg_xyyyyz_xzz_1, tg_xyyyyz_yy_1, tg_xyyyyz_yyy_0, tg_xyyyyz_yyy_1, tg_xyyyyz_yyz_0, \
                                         tg_xyyyyz_yyz_1, tg_xyyyyz_yz_1, tg_xyyyyz_yzz_0, tg_xyyyyz_yzz_1, tg_xyyyyz_zz_1, \
                                         tg_xyyyyz_zzz_0, tg_xyyyyz_zzz_1, tg_xyyyz_xxx_0, tg_xyyyz_xxx_1, tg_xyyyz_xxy_0, \
                                         tg_xyyyz_xxy_1, tg_xyyyz_xxz_0, tg_xyyyz_xxz_1, tg_xyyyz_xyy_0, tg_xyyyz_xyy_1, \
                                         tg_xyyyz_xyz_0, tg_xyyyz_xyz_1, tg_xyyyz_xzz_0, tg_xyyyz_xzz_1, tg_xyyyz_yyy_0, \
                                         tg_xyyyz_yyy_1, tg_xyyyz_yyz_0, tg_xyyyz_yyz_1, tg_xyyyz_yzz_0, tg_xyyyz_yzz_1, \
                                         tg_xyyyz_zzz_0, tg_xyyyz_zzz_1, tg_xyyyzz_xx_1, tg_xyyyzz_xxx_0, tg_xyyyzz_xxx_1, \
                                         tg_xyyyzz_xxy_0, tg_xyyyzz_xxy_1, tg_xyyyzz_xxz_0, tg_xyyyzz_xxz_1, tg_xyyyzz_xy_1, \
                                         tg_xyyyzz_xyy_0, tg_xyyyzz_xyy_1, tg_xyyyzz_xyz_0, tg_xyyyzz_xyz_1, tg_xyyyzz_xz_1, \
                                         tg_xyyyzz_xzz_0, tg_xyyyzz_xzz_1, tg_xyyyzz_yy_1, tg_xyyyzz_yyy_0, tg_xyyyzz_yyy_1, \
                                         tg_xyyyzz_yyz_0, tg_xyyyzz_yyz_1, tg_xyyyzz_yz_1, tg_xyyyzz_yzz_0, tg_xyyyzz_yzz_1, \
                                         tg_xyyyzz_zz_1, tg_xyyyzz_zzz_0, tg_xyyyzz_zzz_1, tg_xyyzz_xxx_0, tg_xyyzz_xxx_1, \
                                         tg_xyyzz_xxy_0, tg_xyyzz_xxy_1, tg_xyyzz_xxz_0, tg_xyyzz_xxz_1, tg_xyyzz_xyy_0, \
                                         tg_xyyzz_xyy_1, tg_xyyzz_xyz_0, tg_xyyzz_xyz_1, tg_xyyzz_xzz_0, tg_xyyzz_xzz_1, \
                                         tg_xyyzz_yyy_0, tg_xyyzz_yyy_1, tg_xyyzz_yyz_0, tg_xyyzz_yyz_1, tg_xyyzz_yzz_0, \
                                         tg_xyyzz_yzz_1, tg_xyyzz_zzz_0, tg_xyyzz_zzz_1, tg_xyzzz_xxx_0, tg_xyzzz_xxx_1, \
                                         tg_xyzzz_xxy_0, tg_xyzzz_xxy_1, tg_xyzzz_xxz_0, tg_xyzzz_xxz_1, tg_xyzzz_xyy_0, \
                                         tg_xyzzz_xyy_1, tg_xyzzz_xyz_0, tg_xyzzz_xyz_1, tg_xyzzz_xzz_0, tg_xyzzz_xzz_1, \
                                         tg_xyzzz_yyy_0, tg_xyzzz_yyy_1, tg_xyzzz_yyz_0, tg_xyzzz_yyz_1, tg_xyzzz_yzz_0, \
                                         tg_xyzzz_yzz_1, tg_xyzzz_zzz_0, tg_xyzzz_zzz_1, tg_xzzzz_xxx_0, tg_xzzzz_xxx_1, \
                                         tg_xzzzz_xxy_0, tg_xzzzz_xxy_1, tg_xzzzz_xxz_0, tg_xzzzz_xxz_1, tg_xzzzz_xyy_0, \
                                         tg_xzzzz_xyy_1, tg_xzzzz_xyz_0, tg_xzzzz_xyz_1, tg_xzzzz_xzz_0, tg_xzzzz_xzz_1, \
                                         tg_xzzzz_yyy_0, tg_xzzzz_yyy_1, tg_xzzzz_yyz_0, tg_xzzzz_yyz_1, tg_xzzzz_yzz_0, \
                                         tg_xzzzz_yzz_1, tg_xzzzz_zzz_0, tg_xzzzz_zzz_1, tg_yyyyy_xxx_0, tg_yyyyy_xxx_1, \
                                         tg_yyyyy_xxy_0, tg_yyyyy_xxy_1, tg_yyyyy_xxz_0, tg_yyyyy_xxz_1, tg_yyyyy_xyy_0, \
                                         tg_yyyyy_xyy_1, tg_yyyyy_xyz_0, tg_yyyyy_xyz_1, tg_yyyyy_xzz_0, tg_yyyyy_xzz_1, \
                                         tg_yyyyy_yyy_0, tg_yyyyy_yyy_1, tg_yyyyy_yyz_0, tg_yyyyy_yyz_1, tg_yyyyy_yzz_0, \
                                         tg_yyyyy_yzz_1, tg_yyyyy_zzz_0, tg_yyyyy_zzz_1, tg_yyyyz_xxx_0, tg_yyyyz_xxx_1, \
                                         tg_yyyyz_xxy_0, tg_yyyyz_xxy_1, tg_yyyyz_xxz_0, tg_yyyyz_xxz_1, tg_yyyyz_xyy_0, \
                                         tg_yyyyz_xyy_1, tg_yyyyz_xyz_0, tg_yyyyz_xyz_1, tg_yyyyz_xzz_0, tg_yyyyz_xzz_1, \
                                         tg_yyyyz_yyy_0, tg_yyyyz_yyy_1, tg_yyyyz_yyz_0, tg_yyyyz_yyz_1, tg_yyyyz_yzz_0, \
                                         tg_yyyyz_yzz_1, tg_yyyyz_zzz_0, tg_yyyyz_zzz_1, tg_yyyzz_xxx_0, tg_yyyzz_xxx_1, \
                                         tg_yyyzz_xxy_0, tg_yyyzz_xxy_1, tg_yyyzz_xxz_0, tg_yyyzz_xxz_1, tg_yyyzz_xyy_0, \
                                         tg_yyyzz_xyy_1, tg_yyyzz_xyz_0, tg_yyyzz_xyz_1, tg_yyyzz_xzz_0, tg_yyyzz_xzz_1, \
                                         tg_yyyzz_yyy_0, tg_yyyzz_yyy_1, tg_yyyzz_yyz_0, tg_yyyzz_yyz_1, tg_yyyzz_yzz_0, \
                                         tg_yyyzz_yzz_1, tg_yyyzz_zzz_0, tg_yyyzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxzzz_xxx_0[j] = pb_x * tg_xxxzzz_xxx_0[j] + fr * tg_xxxzzz_xxx_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxx_0[j] - tg_xxzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzz_xx_1[j];

                    tg_xxxxzzz_xxy_0[j] = pb_x * tg_xxxzzz_xxy_0[j] + fr * tg_xxxzzz_xxy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxy_0[j] - tg_xxzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xy_1[j];

                    tg_xxxxzzz_xxz_0[j] = pb_x * tg_xxxzzz_xxz_0[j] + fr * tg_xxxzzz_xxz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xxz_0[j] - tg_xxzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzz_xz_1[j];

                    tg_xxxxzzz_xyy_0[j] = pb_x * tg_xxxzzz_xyy_0[j] + fr * tg_xxxzzz_xyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyy_0[j] - tg_xxzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yy_1[j];

                    tg_xxxxzzz_xyz_0[j] = pb_x * tg_xxxzzz_xyz_0[j] + fr * tg_xxxzzz_xyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xyz_0[j] - tg_xxzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_yz_1[j];

                    tg_xxxxzzz_xzz_0[j] = pb_x * tg_xxxzzz_xzz_0[j] + fr * tg_xxxzzz_xzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_xzz_0[j] - tg_xxzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzz_zz_1[j];

                    tg_xxxxzzz_yyy_0[j] = pb_x * tg_xxxzzz_yyy_0[j] + fr * tg_xxxzzz_yyy_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyy_0[j] - tg_xxzzz_yyy_1[j] * fl1_fza);

                    tg_xxxxzzz_yyz_0[j] = pb_x * tg_xxxzzz_yyz_0[j] + fr * tg_xxxzzz_yyz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yyz_0[j] - tg_xxzzz_yyz_1[j] * fl1_fza);

                    tg_xxxxzzz_yzz_0[j] = pb_x * tg_xxxzzz_yzz_0[j] + fr * tg_xxxzzz_yzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_yzz_0[j] - tg_xxzzz_yzz_1[j] * fl1_fza);

                    tg_xxxxzzz_zzz_0[j] = pb_x * tg_xxxzzz_zzz_0[j] + fr * tg_xxxzzz_zzz_1[j] + 1.5 * fl1_fx * (tg_xxzzz_zzz_0[j] - tg_xxzzz_zzz_1[j] * fl1_fza);

                    tg_xxxyyyy_xxx_0[j] = pb_x * tg_xxyyyy_xxx_0[j] + fr * tg_xxyyyy_xxx_1[j] + fl1_fx * (tg_xyyyy_xxx_0[j] - tg_xyyyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyy_xx_1[j];

                    tg_xxxyyyy_xxy_0[j] = pb_x * tg_xxyyyy_xxy_0[j] + fr * tg_xxyyyy_xxy_1[j] + fl1_fx * (tg_xyyyy_xxy_0[j] - tg_xyyyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xy_1[j];

                    tg_xxxyyyy_xxz_0[j] = pb_x * tg_xxyyyy_xxz_0[j] + fr * tg_xxyyyy_xxz_1[j] + fl1_fx * (tg_xyyyy_xxz_0[j] - tg_xyyyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyy_xz_1[j];

                    tg_xxxyyyy_xyy_0[j] = pb_x * tg_xxyyyy_xyy_0[j] + fr * tg_xxyyyy_xyy_1[j] + fl1_fx * (tg_xyyyy_xyy_0[j] - tg_xyyyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yy_1[j];

                    tg_xxxyyyy_xyz_0[j] = pb_x * tg_xxyyyy_xyz_0[j] + fr * tg_xxyyyy_xyz_1[j] + fl1_fx * (tg_xyyyy_xyz_0[j] - tg_xyyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_yz_1[j];

                    tg_xxxyyyy_xzz_0[j] = pb_x * tg_xxyyyy_xzz_0[j] + fr * tg_xxyyyy_xzz_1[j] + fl1_fx * (tg_xyyyy_xzz_0[j] - tg_xyyyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyy_zz_1[j];

                    tg_xxxyyyy_yyy_0[j] = pb_x * tg_xxyyyy_yyy_0[j] + fr * tg_xxyyyy_yyy_1[j] + fl1_fx * (tg_xyyyy_yyy_0[j] - tg_xyyyy_yyy_1[j] * fl1_fza);

                    tg_xxxyyyy_yyz_0[j] = pb_x * tg_xxyyyy_yyz_0[j] + fr * tg_xxyyyy_yyz_1[j] + fl1_fx * (tg_xyyyy_yyz_0[j] - tg_xyyyy_yyz_1[j] * fl1_fza);

                    tg_xxxyyyy_yzz_0[j] = pb_x * tg_xxyyyy_yzz_0[j] + fr * tg_xxyyyy_yzz_1[j] + fl1_fx * (tg_xyyyy_yzz_0[j] - tg_xyyyy_yzz_1[j] * fl1_fza);

                    tg_xxxyyyy_zzz_0[j] = pb_x * tg_xxyyyy_zzz_0[j] + fr * tg_xxyyyy_zzz_1[j] + fl1_fx * (tg_xyyyy_zzz_0[j] - tg_xyyyy_zzz_1[j] * fl1_fza);

                    tg_xxxyyyz_xxx_0[j] = pb_x * tg_xxyyyz_xxx_0[j] + fr * tg_xxyyyz_xxx_1[j] + fl1_fx * (tg_xyyyz_xxx_0[j] - tg_xyyyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyz_xx_1[j];

                    tg_xxxyyyz_xxy_0[j] = pb_x * tg_xxyyyz_xxy_0[j] + fr * tg_xxyyyz_xxy_1[j] + fl1_fx * (tg_xyyyz_xxy_0[j] - tg_xyyyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xy_1[j];

                    tg_xxxyyyz_xxz_0[j] = pb_x * tg_xxyyyz_xxz_0[j] + fr * tg_xxyyyz_xxz_1[j] + fl1_fx * (tg_xyyyz_xxz_0[j] - tg_xyyyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyz_xz_1[j];

                    tg_xxxyyyz_xyy_0[j] = pb_x * tg_xxyyyz_xyy_0[j] + fr * tg_xxyyyz_xyy_1[j] + fl1_fx * (tg_xyyyz_xyy_0[j] - tg_xyyyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yy_1[j];

                    tg_xxxyyyz_xyz_0[j] = pb_x * tg_xxyyyz_xyz_0[j] + fr * tg_xxyyyz_xyz_1[j] + fl1_fx * (tg_xyyyz_xyz_0[j] - tg_xyyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_yz_1[j];

                    tg_xxxyyyz_xzz_0[j] = pb_x * tg_xxyyyz_xzz_0[j] + fr * tg_xxyyyz_xzz_1[j] + fl1_fx * (tg_xyyyz_xzz_0[j] - tg_xyyyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyz_zz_1[j];

                    tg_xxxyyyz_yyy_0[j] = pb_x * tg_xxyyyz_yyy_0[j] + fr * tg_xxyyyz_yyy_1[j] + fl1_fx * (tg_xyyyz_yyy_0[j] - tg_xyyyz_yyy_1[j] * fl1_fza);

                    tg_xxxyyyz_yyz_0[j] = pb_x * tg_xxyyyz_yyz_0[j] + fr * tg_xxyyyz_yyz_1[j] + fl1_fx * (tg_xyyyz_yyz_0[j] - tg_xyyyz_yyz_1[j] * fl1_fza);

                    tg_xxxyyyz_yzz_0[j] = pb_x * tg_xxyyyz_yzz_0[j] + fr * tg_xxyyyz_yzz_1[j] + fl1_fx * (tg_xyyyz_yzz_0[j] - tg_xyyyz_yzz_1[j] * fl1_fza);

                    tg_xxxyyyz_zzz_0[j] = pb_x * tg_xxyyyz_zzz_0[j] + fr * tg_xxyyyz_zzz_1[j] + fl1_fx * (tg_xyyyz_zzz_0[j] - tg_xyyyz_zzz_1[j] * fl1_fza);

                    tg_xxxyyzz_xxx_0[j] = pb_x * tg_xxyyzz_xxx_0[j] + fr * tg_xxyyzz_xxx_1[j] + fl1_fx * (tg_xyyzz_xxx_0[j] - tg_xyyzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzz_xx_1[j];

                    tg_xxxyyzz_xxy_0[j] = pb_x * tg_xxyyzz_xxy_0[j] + fr * tg_xxyyzz_xxy_1[j] + fl1_fx * (tg_xyyzz_xxy_0[j] - tg_xyyzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xy_1[j];

                    tg_xxxyyzz_xxz_0[j] = pb_x * tg_xxyyzz_xxz_0[j] + fr * tg_xxyyzz_xxz_1[j] + fl1_fx * (tg_xyyzz_xxz_0[j] - tg_xyyzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzz_xz_1[j];

                    tg_xxxyyzz_xyy_0[j] = pb_x * tg_xxyyzz_xyy_0[j] + fr * tg_xxyyzz_xyy_1[j] + fl1_fx * (tg_xyyzz_xyy_0[j] - tg_xyyzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yy_1[j];

                    tg_xxxyyzz_xyz_0[j] = pb_x * tg_xxyyzz_xyz_0[j] + fr * tg_xxyyzz_xyz_1[j] + fl1_fx * (tg_xyyzz_xyz_0[j] - tg_xyyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_yz_1[j];

                    tg_xxxyyzz_xzz_0[j] = pb_x * tg_xxyyzz_xzz_0[j] + fr * tg_xxyyzz_xzz_1[j] + fl1_fx * (tg_xyyzz_xzz_0[j] - tg_xyyzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzz_zz_1[j];

                    tg_xxxyyzz_yyy_0[j] = pb_x * tg_xxyyzz_yyy_0[j] + fr * tg_xxyyzz_yyy_1[j] + fl1_fx * (tg_xyyzz_yyy_0[j] - tg_xyyzz_yyy_1[j] * fl1_fza);

                    tg_xxxyyzz_yyz_0[j] = pb_x * tg_xxyyzz_yyz_0[j] + fr * tg_xxyyzz_yyz_1[j] + fl1_fx * (tg_xyyzz_yyz_0[j] - tg_xyyzz_yyz_1[j] * fl1_fza);

                    tg_xxxyyzz_yzz_0[j] = pb_x * tg_xxyyzz_yzz_0[j] + fr * tg_xxyyzz_yzz_1[j] + fl1_fx * (tg_xyyzz_yzz_0[j] - tg_xyyzz_yzz_1[j] * fl1_fza);

                    tg_xxxyyzz_zzz_0[j] = pb_x * tg_xxyyzz_zzz_0[j] + fr * tg_xxyyzz_zzz_1[j] + fl1_fx * (tg_xyyzz_zzz_0[j] - tg_xyyzz_zzz_1[j] * fl1_fza);

                    tg_xxxyzzz_xxx_0[j] = pb_x * tg_xxyzzz_xxx_0[j] + fr * tg_xxyzzz_xxx_1[j] + fl1_fx * (tg_xyzzz_xxx_0[j] - tg_xyzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzz_xx_1[j];

                    tg_xxxyzzz_xxy_0[j] = pb_x * tg_xxyzzz_xxy_0[j] + fr * tg_xxyzzz_xxy_1[j] + fl1_fx * (tg_xyzzz_xxy_0[j] - tg_xyzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xy_1[j];

                    tg_xxxyzzz_xxz_0[j] = pb_x * tg_xxyzzz_xxz_0[j] + fr * tg_xxyzzz_xxz_1[j] + fl1_fx * (tg_xyzzz_xxz_0[j] - tg_xyzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzz_xz_1[j];

                    tg_xxxyzzz_xyy_0[j] = pb_x * tg_xxyzzz_xyy_0[j] + fr * tg_xxyzzz_xyy_1[j] + fl1_fx * (tg_xyzzz_xyy_0[j] - tg_xyzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yy_1[j];

                    tg_xxxyzzz_xyz_0[j] = pb_x * tg_xxyzzz_xyz_0[j] + fr * tg_xxyzzz_xyz_1[j] + fl1_fx * (tg_xyzzz_xyz_0[j] - tg_xyzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_yz_1[j];

                    tg_xxxyzzz_xzz_0[j] = pb_x * tg_xxyzzz_xzz_0[j] + fr * tg_xxyzzz_xzz_1[j] + fl1_fx * (tg_xyzzz_xzz_0[j] - tg_xyzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzz_zz_1[j];

                    tg_xxxyzzz_yyy_0[j] = pb_x * tg_xxyzzz_yyy_0[j] + fr * tg_xxyzzz_yyy_1[j] + fl1_fx * (tg_xyzzz_yyy_0[j] - tg_xyzzz_yyy_1[j] * fl1_fza);

                    tg_xxxyzzz_yyz_0[j] = pb_x * tg_xxyzzz_yyz_0[j] + fr * tg_xxyzzz_yyz_1[j] + fl1_fx * (tg_xyzzz_yyz_0[j] - tg_xyzzz_yyz_1[j] * fl1_fza);

                    tg_xxxyzzz_yzz_0[j] = pb_x * tg_xxyzzz_yzz_0[j] + fr * tg_xxyzzz_yzz_1[j] + fl1_fx * (tg_xyzzz_yzz_0[j] - tg_xyzzz_yzz_1[j] * fl1_fza);

                    tg_xxxyzzz_zzz_0[j] = pb_x * tg_xxyzzz_zzz_0[j] + fr * tg_xxyzzz_zzz_1[j] + fl1_fx * (tg_xyzzz_zzz_0[j] - tg_xyzzz_zzz_1[j] * fl1_fza);

                    tg_xxxzzzz_xxx_0[j] = pb_x * tg_xxzzzz_xxx_0[j] + fr * tg_xxzzzz_xxx_1[j] + fl1_fx * (tg_xzzzz_xxx_0[j] - tg_xzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzz_xx_1[j];

                    tg_xxxzzzz_xxy_0[j] = pb_x * tg_xxzzzz_xxy_0[j] + fr * tg_xxzzzz_xxy_1[j] + fl1_fx * (tg_xzzzz_xxy_0[j] - tg_xzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xy_1[j];

                    tg_xxxzzzz_xxz_0[j] = pb_x * tg_xxzzzz_xxz_0[j] + fr * tg_xxzzzz_xxz_1[j] + fl1_fx * (tg_xzzzz_xxz_0[j] - tg_xzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzz_xz_1[j];

                    tg_xxxzzzz_xyy_0[j] = pb_x * tg_xxzzzz_xyy_0[j] + fr * tg_xxzzzz_xyy_1[j] + fl1_fx * (tg_xzzzz_xyy_0[j] - tg_xzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yy_1[j];

                    tg_xxxzzzz_xyz_0[j] = pb_x * tg_xxzzzz_xyz_0[j] + fr * tg_xxzzzz_xyz_1[j] + fl1_fx * (tg_xzzzz_xyz_0[j] - tg_xzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_yz_1[j];

                    tg_xxxzzzz_xzz_0[j] = pb_x * tg_xxzzzz_xzz_0[j] + fr * tg_xxzzzz_xzz_1[j] + fl1_fx * (tg_xzzzz_xzz_0[j] - tg_xzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzz_zz_1[j];

                    tg_xxxzzzz_yyy_0[j] = pb_x * tg_xxzzzz_yyy_0[j] + fr * tg_xxzzzz_yyy_1[j] + fl1_fx * (tg_xzzzz_yyy_0[j] - tg_xzzzz_yyy_1[j] * fl1_fza);

                    tg_xxxzzzz_yyz_0[j] = pb_x * tg_xxzzzz_yyz_0[j] + fr * tg_xxzzzz_yyz_1[j] + fl1_fx * (tg_xzzzz_yyz_0[j] - tg_xzzzz_yyz_1[j] * fl1_fza);

                    tg_xxxzzzz_yzz_0[j] = pb_x * tg_xxzzzz_yzz_0[j] + fr * tg_xxzzzz_yzz_1[j] + fl1_fx * (tg_xzzzz_yzz_0[j] - tg_xzzzz_yzz_1[j] * fl1_fza);

                    tg_xxxzzzz_zzz_0[j] = pb_x * tg_xxzzzz_zzz_0[j] + fr * tg_xxzzzz_zzz_1[j] + fl1_fx * (tg_xzzzz_zzz_0[j] - tg_xzzzz_zzz_1[j] * fl1_fza);

                    tg_xxyyyyy_xxx_0[j] = pb_x * tg_xyyyyy_xxx_0[j] + fr * tg_xyyyyy_xxx_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxx_0[j] - tg_yyyyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyy_xx_1[j];

                    tg_xxyyyyy_xxy_0[j] = pb_x * tg_xyyyyy_xxy_0[j] + fr * tg_xyyyyy_xxy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxy_0[j] - tg_yyyyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xy_1[j];

                    tg_xxyyyyy_xxz_0[j] = pb_x * tg_xyyyyy_xxz_0[j] + fr * tg_xyyyyy_xxz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xxz_0[j] - tg_yyyyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyy_xz_1[j];

                    tg_xxyyyyy_xyy_0[j] = pb_x * tg_xyyyyy_xyy_0[j] + fr * tg_xyyyyy_xyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyy_0[j] - tg_yyyyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yy_1[j];

                    tg_xxyyyyy_xyz_0[j] = pb_x * tg_xyyyyy_xyz_0[j] + fr * tg_xyyyyy_xyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xyz_0[j] - tg_yyyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_yz_1[j];

                    tg_xxyyyyy_xzz_0[j] = pb_x * tg_xyyyyy_xzz_0[j] + fr * tg_xyyyyy_xzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_xzz_0[j] - tg_yyyyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyy_zz_1[j];

                    tg_xxyyyyy_yyy_0[j] = pb_x * tg_xyyyyy_yyy_0[j] + fr * tg_xyyyyy_yyy_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyy_0[j] - tg_yyyyy_yyy_1[j] * fl1_fza);

                    tg_xxyyyyy_yyz_0[j] = pb_x * tg_xyyyyy_yyz_0[j] + fr * tg_xyyyyy_yyz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yyz_0[j] - tg_yyyyy_yyz_1[j] * fl1_fza);

                    tg_xxyyyyy_yzz_0[j] = pb_x * tg_xyyyyy_yzz_0[j] + fr * tg_xyyyyy_yzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_yzz_0[j] - tg_yyyyy_yzz_1[j] * fl1_fza);

                    tg_xxyyyyy_zzz_0[j] = pb_x * tg_xyyyyy_zzz_0[j] + fr * tg_xyyyyy_zzz_1[j] + 0.5 * fl1_fx * (tg_yyyyy_zzz_0[j] - tg_yyyyy_zzz_1[j] * fl1_fza);

                    tg_xxyyyyz_xxx_0[j] = pb_x * tg_xyyyyz_xxx_0[j] + fr * tg_xyyyyz_xxx_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxx_0[j] - tg_yyyyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyz_xx_1[j];

                    tg_xxyyyyz_xxy_0[j] = pb_x * tg_xyyyyz_xxy_0[j] + fr * tg_xyyyyz_xxy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxy_0[j] - tg_yyyyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xy_1[j];

                    tg_xxyyyyz_xxz_0[j] = pb_x * tg_xyyyyz_xxz_0[j] + fr * tg_xyyyyz_xxz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xxz_0[j] - tg_yyyyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyz_xz_1[j];

                    tg_xxyyyyz_xyy_0[j] = pb_x * tg_xyyyyz_xyy_0[j] + fr * tg_xyyyyz_xyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyy_0[j] - tg_yyyyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yy_1[j];

                    tg_xxyyyyz_xyz_0[j] = pb_x * tg_xyyyyz_xyz_0[j] + fr * tg_xyyyyz_xyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xyz_0[j] - tg_yyyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_yz_1[j];

                    tg_xxyyyyz_xzz_0[j] = pb_x * tg_xyyyyz_xzz_0[j] + fr * tg_xyyyyz_xzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_xzz_0[j] - tg_yyyyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyz_zz_1[j];

                    tg_xxyyyyz_yyy_0[j] = pb_x * tg_xyyyyz_yyy_0[j] + fr * tg_xyyyyz_yyy_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyy_0[j] - tg_yyyyz_yyy_1[j] * fl1_fza);

                    tg_xxyyyyz_yyz_0[j] = pb_x * tg_xyyyyz_yyz_0[j] + fr * tg_xyyyyz_yyz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yyz_0[j] - tg_yyyyz_yyz_1[j] * fl1_fza);

                    tg_xxyyyyz_yzz_0[j] = pb_x * tg_xyyyyz_yzz_0[j] + fr * tg_xyyyyz_yzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_yzz_0[j] - tg_yyyyz_yzz_1[j] * fl1_fza);

                    tg_xxyyyyz_zzz_0[j] = pb_x * tg_xyyyyz_zzz_0[j] + fr * tg_xyyyyz_zzz_1[j] + 0.5 * fl1_fx * (tg_yyyyz_zzz_0[j] - tg_yyyyz_zzz_1[j] * fl1_fza);

                    tg_xxyyyzz_xxx_0[j] = pb_x * tg_xyyyzz_xxx_0[j] + fr * tg_xyyyzz_xxx_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxx_0[j] - tg_yyyzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzz_xx_1[j];

                    tg_xxyyyzz_xxy_0[j] = pb_x * tg_xyyyzz_xxy_0[j] + fr * tg_xyyyzz_xxy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxy_0[j] - tg_yyyzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xy_1[j];

                    tg_xxyyyzz_xxz_0[j] = pb_x * tg_xyyyzz_xxz_0[j] + fr * tg_xyyyzz_xxz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xxz_0[j] - tg_yyyzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzz_xz_1[j];

                    tg_xxyyyzz_xyy_0[j] = pb_x * tg_xyyyzz_xyy_0[j] + fr * tg_xyyyzz_xyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyy_0[j] - tg_yyyzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yy_1[j];

                    tg_xxyyyzz_xyz_0[j] = pb_x * tg_xyyyzz_xyz_0[j] + fr * tg_xyyyzz_xyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xyz_0[j] - tg_yyyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_yz_1[j];

                    tg_xxyyyzz_xzz_0[j] = pb_x * tg_xyyyzz_xzz_0[j] + fr * tg_xyyyzz_xzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_xzz_0[j] - tg_yyyzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzz_zz_1[j];

                    tg_xxyyyzz_yyy_0[j] = pb_x * tg_xyyyzz_yyy_0[j] + fr * tg_xyyyzz_yyy_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyy_0[j] - tg_yyyzz_yyy_1[j] * fl1_fza);

                    tg_xxyyyzz_yyz_0[j] = pb_x * tg_xyyyzz_yyz_0[j] + fr * tg_xyyyzz_yyz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yyz_0[j] - tg_yyyzz_yyz_1[j] * fl1_fza);

                    tg_xxyyyzz_yzz_0[j] = pb_x * tg_xyyyzz_yzz_0[j] + fr * tg_xyyyzz_yzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_yzz_0[j] - tg_yyyzz_yzz_1[j] * fl1_fza);

                    tg_xxyyyzz_zzz_0[j] = pb_x * tg_xyyyzz_zzz_0[j] + fr * tg_xyyyzz_zzz_1[j] + 0.5 * fl1_fx * (tg_yyyzz_zzz_0[j] - tg_yyyzz_zzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSF_180_270(      CMemBlock2D<double>* primBuffer,
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
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xyyzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 180); 

                auto tg_xyyzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 181); 

                auto tg_xyyzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 182); 

                auto tg_xyyzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 183); 

                auto tg_xyyzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 184); 

                auto tg_xyyzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 185); 

                auto tg_xyyzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 186); 

                auto tg_xyyzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 187); 

                auto tg_xyyzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 188); 

                auto tg_xyyzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 189); 

                auto tg_xyzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 190); 

                auto tg_xyzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 191); 

                auto tg_xyzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 192); 

                auto tg_xyzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 193); 

                auto tg_xyzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 194); 

                auto tg_xyzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 195); 

                auto tg_xyzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 196); 

                auto tg_xyzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 197); 

                auto tg_xyzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 198); 

                auto tg_xyzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 199); 

                auto tg_xzzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 200); 

                auto tg_xzzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 201); 

                auto tg_xzzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 202); 

                auto tg_xzzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 203); 

                auto tg_xzzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 204); 

                auto tg_xzzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 205); 

                auto tg_xzzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 206); 

                auto tg_xzzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 207); 

                auto tg_xzzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 208); 

                auto tg_xzzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 209); 

                auto tg_yyyyyy_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 210); 

                auto tg_yyyyyy_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 211); 

                auto tg_yyyyyy_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 212); 

                auto tg_yyyyyy_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 213); 

                auto tg_yyyyyy_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 214); 

                auto tg_yyyyyy_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 215); 

                auto tg_yyyyyy_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 216); 

                auto tg_yyyyyy_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 217); 

                auto tg_yyyyyy_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 218); 

                auto tg_yyyyyy_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 219); 

                auto tg_yyyyyz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 220); 

                auto tg_yyyyyz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 221); 

                auto tg_yyyyyz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 222); 

                auto tg_yyyyyz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 223); 

                auto tg_yyyyyz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 224); 

                auto tg_yyyyyz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 225); 

                auto tg_yyyyyz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 226); 

                auto tg_yyyyyz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 227); 

                auto tg_yyyyyz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 228); 

                auto tg_yyyyyz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 229); 

                auto tg_yyyyzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 230); 

                auto tg_yyyyzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 231); 

                auto tg_yyyyzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 232); 

                auto tg_yyyyzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 233); 

                auto tg_yyyyzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 234); 

                auto tg_yyyyzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 235); 

                auto tg_yyyyzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 236); 

                auto tg_yyyyzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 237); 

                auto tg_yyyyzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 238); 

                auto tg_yyyyzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 239); 

                auto tg_yyyzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 240); 

                auto tg_yyyzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 241); 

                auto tg_yyyzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 242); 

                auto tg_yyyzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 243); 

                auto tg_yyyzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 244); 

                auto tg_yyyzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 245); 

                auto tg_yyyzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 246); 

                auto tg_yyyzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 247); 

                auto tg_yyyzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 248); 

                auto tg_yyyzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 249); 

                auto tg_yyzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 250); 

                auto tg_yyzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 251); 

                auto tg_yyzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 252); 

                auto tg_yyzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 253); 

                auto tg_yyzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 254); 

                auto tg_yyzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 255); 

                auto tg_yyzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 256); 

                auto tg_yyzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 257); 

                auto tg_yyzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 258); 

                auto tg_yyzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 259); 

                auto tg_yzzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 260); 

                auto tg_yzzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 261); 

                auto tg_yzzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 262); 

                auto tg_yzzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 263); 

                auto tg_yzzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 264); 

                auto tg_yzzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 265); 

                auto tg_yzzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 266); 

                auto tg_yzzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 267); 

                auto tg_yzzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 268); 

                auto tg_yzzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 269); 

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

                auto tg_yyzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 180); 

                auto tg_yyzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 181); 

                auto tg_yyzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 182); 

                auto tg_yyzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 183); 

                auto tg_yyzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 184); 

                auto tg_yyzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 185); 

                auto tg_yyzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 186); 

                auto tg_yyzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 187); 

                auto tg_yyzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 188); 

                auto tg_yyzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 189); 

                auto tg_yzzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 190); 

                auto tg_yzzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 191); 

                auto tg_yzzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 192); 

                auto tg_yzzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 193); 

                auto tg_yzzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 194); 

                auto tg_yzzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 195); 

                auto tg_yzzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 196); 

                auto tg_yzzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 197); 

                auto tg_yzzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 198); 

                auto tg_yzzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 199); 

                auto tg_zzzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 200); 

                auto tg_zzzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 201); 

                auto tg_zzzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 202); 

                auto tg_zzzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 203); 

                auto tg_zzzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 204); 

                auto tg_zzzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 205); 

                auto tg_zzzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 206); 

                auto tg_zzzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 207); 

                auto tg_zzzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 208); 

                auto tg_zzzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 209); 

                auto tg_yyzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 180); 

                auto tg_yyzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 181); 

                auto tg_yyzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 182); 

                auto tg_yyzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 183); 

                auto tg_yyzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 184); 

                auto tg_yyzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 185); 

                auto tg_yyzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 186); 

                auto tg_yyzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 187); 

                auto tg_yyzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 188); 

                auto tg_yyzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 189); 

                auto tg_yzzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 190); 

                auto tg_yzzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 191); 

                auto tg_yzzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 192); 

                auto tg_yzzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 193); 

                auto tg_yzzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 194); 

                auto tg_yzzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 195); 

                auto tg_yzzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 196); 

                auto tg_yzzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 197); 

                auto tg_yzzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 198); 

                auto tg_yzzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 199); 

                auto tg_zzzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 200); 

                auto tg_zzzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 201); 

                auto tg_zzzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 202); 

                auto tg_zzzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 203); 

                auto tg_zzzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 204); 

                auto tg_zzzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 205); 

                auto tg_zzzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 206); 

                auto tg_zzzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 207); 

                auto tg_zzzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 208); 

                auto tg_zzzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 209); 

                auto tg_xyyzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 108); 

                auto tg_xyyzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 109); 

                auto tg_xyyzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 110); 

                auto tg_xyyzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 111); 

                auto tg_xyyzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 112); 

                auto tg_xyyzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 113); 

                auto tg_xyzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 114); 

                auto tg_xyzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 115); 

                auto tg_xyzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 116); 

                auto tg_xyzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 117); 

                auto tg_xyzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 118); 

                auto tg_xyzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 119); 

                auto tg_xzzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 120); 

                auto tg_xzzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 121); 

                auto tg_xzzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 122); 

                auto tg_xzzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 123); 

                auto tg_xzzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 124); 

                auto tg_xzzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 125); 

                auto tg_yyyyyy_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 126); 

                auto tg_yyyyyy_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 127); 

                auto tg_yyyyyy_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 128); 

                auto tg_yyyyyy_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 129); 

                auto tg_yyyyyy_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 130); 

                auto tg_yyyyyy_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 131); 

                auto tg_yyyyyz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 132); 

                auto tg_yyyyyz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 133); 

                auto tg_yyyyyz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 134); 

                auto tg_yyyyyz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 135); 

                auto tg_yyyyyz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 136); 

                auto tg_yyyyyz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 137); 

                auto tg_yyyyzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 138); 

                auto tg_yyyyzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 139); 

                auto tg_yyyyzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 140); 

                auto tg_yyyyzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 141); 

                auto tg_yyyyzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 142); 

                auto tg_yyyyzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 143); 

                auto tg_yyyzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 144); 

                auto tg_yyyzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 145); 

                auto tg_yyyzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 146); 

                auto tg_yyyzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 147); 

                auto tg_yyyzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 148); 

                auto tg_yyyzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 149); 

                auto tg_yyzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 150); 

                auto tg_yyzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 151); 

                auto tg_yyzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 152); 

                auto tg_yyzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 153); 

                auto tg_yyzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 154); 

                auto tg_yyzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 155); 

                auto tg_yzzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 156); 

                auto tg_yzzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 157); 

                auto tg_yzzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 158); 

                auto tg_yzzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 159); 

                auto tg_yzzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 160); 

                auto tg_yzzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 161); 

                // set up pointers to integrals

                auto tg_xxyyzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 180); 

                auto tg_xxyyzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 181); 

                auto tg_xxyyzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 182); 

                auto tg_xxyyzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 183); 

                auto tg_xxyyzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 184); 

                auto tg_xxyyzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 185); 

                auto tg_xxyyzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 186); 

                auto tg_xxyyzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 187); 

                auto tg_xxyyzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 188); 

                auto tg_xxyyzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 189); 

                auto tg_xxyzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 190); 

                auto tg_xxyzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 191); 

                auto tg_xxyzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 192); 

                auto tg_xxyzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 193); 

                auto tg_xxyzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 194); 

                auto tg_xxyzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 195); 

                auto tg_xxyzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 196); 

                auto tg_xxyzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 197); 

                auto tg_xxyzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 198); 

                auto tg_xxyzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 199); 

                auto tg_xxzzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 200); 

                auto tg_xxzzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 201); 

                auto tg_xxzzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 202); 

                auto tg_xxzzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 203); 

                auto tg_xxzzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 204); 

                auto tg_xxzzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 205); 

                auto tg_xxzzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 206); 

                auto tg_xxzzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 207); 

                auto tg_xxzzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 208); 

                auto tg_xxzzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 209); 

                auto tg_xyyyyyy_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 210); 

                auto tg_xyyyyyy_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 211); 

                auto tg_xyyyyyy_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 212); 

                auto tg_xyyyyyy_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 213); 

                auto tg_xyyyyyy_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 214); 

                auto tg_xyyyyyy_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 215); 

                auto tg_xyyyyyy_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 216); 

                auto tg_xyyyyyy_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 217); 

                auto tg_xyyyyyy_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 218); 

                auto tg_xyyyyyy_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 219); 

                auto tg_xyyyyyz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 220); 

                auto tg_xyyyyyz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 221); 

                auto tg_xyyyyyz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 222); 

                auto tg_xyyyyyz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 223); 

                auto tg_xyyyyyz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 224); 

                auto tg_xyyyyyz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 225); 

                auto tg_xyyyyyz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 226); 

                auto tg_xyyyyyz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 227); 

                auto tg_xyyyyyz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 228); 

                auto tg_xyyyyyz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 229); 

                auto tg_xyyyyzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 230); 

                auto tg_xyyyyzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 231); 

                auto tg_xyyyyzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 232); 

                auto tg_xyyyyzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 233); 

                auto tg_xyyyyzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 234); 

                auto tg_xyyyyzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 235); 

                auto tg_xyyyyzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 236); 

                auto tg_xyyyyzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 237); 

                auto tg_xyyyyzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 238); 

                auto tg_xyyyyzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 239); 

                auto tg_xyyyzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 240); 

                auto tg_xyyyzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 241); 

                auto tg_xyyyzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 242); 

                auto tg_xyyyzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 243); 

                auto tg_xyyyzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 244); 

                auto tg_xyyyzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 245); 

                auto tg_xyyyzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 246); 

                auto tg_xyyyzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 247); 

                auto tg_xyyyzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 248); 

                auto tg_xyyyzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 249); 

                auto tg_xyyzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 250); 

                auto tg_xyyzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 251); 

                auto tg_xyyzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 252); 

                auto tg_xyyzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 253); 

                auto tg_xyyzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 254); 

                auto tg_xyyzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 255); 

                auto tg_xyyzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 256); 

                auto tg_xyyzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 257); 

                auto tg_xyyzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 258); 

                auto tg_xyyzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 259); 

                auto tg_xyzzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 260); 

                auto tg_xyzzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 261); 

                auto tg_xyzzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 262); 

                auto tg_xyzzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 263); 

                auto tg_xyzzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 264); 

                auto tg_xyzzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 265); 

                auto tg_xyzzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 266); 

                auto tg_xyzzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 267); 

                auto tg_xyzzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 268); 

                auto tg_xyzzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 269); 

                // Batch of Integrals (180,270)

                #pragma omp simd aligned(fxn, fza, tg_xxyyzzz_xxx_0, tg_xxyyzzz_xxy_0, tg_xxyyzzz_xxz_0, \
                                         tg_xxyyzzz_xyy_0, tg_xxyyzzz_xyz_0, tg_xxyyzzz_xzz_0, tg_xxyyzzz_yyy_0, \
                                         tg_xxyyzzz_yyz_0, tg_xxyyzzz_yzz_0, tg_xxyyzzz_zzz_0, tg_xxyzzzz_xxx_0, \
                                         tg_xxyzzzz_xxy_0, tg_xxyzzzz_xxz_0, tg_xxyzzzz_xyy_0, tg_xxyzzzz_xyz_0, \
                                         tg_xxyzzzz_xzz_0, tg_xxyzzzz_yyy_0, tg_xxyzzzz_yyz_0, tg_xxyzzzz_yzz_0, \
                                         tg_xxyzzzz_zzz_0, tg_xxzzzzz_xxx_0, tg_xxzzzzz_xxy_0, tg_xxzzzzz_xxz_0, \
                                         tg_xxzzzzz_xyy_0, tg_xxzzzzz_xyz_0, tg_xxzzzzz_xzz_0, tg_xxzzzzz_yyy_0, \
                                         tg_xxzzzzz_yyz_0, tg_xxzzzzz_yzz_0, tg_xxzzzzz_zzz_0, tg_xyyyyyy_xxx_0, \
                                         tg_xyyyyyy_xxy_0, tg_xyyyyyy_xxz_0, tg_xyyyyyy_xyy_0, tg_xyyyyyy_xyz_0, \
                                         tg_xyyyyyy_xzz_0, tg_xyyyyyy_yyy_0, tg_xyyyyyy_yyz_0, tg_xyyyyyy_yzz_0, \
                                         tg_xyyyyyy_zzz_0, tg_xyyyyyz_xxx_0, tg_xyyyyyz_xxy_0, tg_xyyyyyz_xxz_0, \
                                         tg_xyyyyyz_xyy_0, tg_xyyyyyz_xyz_0, tg_xyyyyyz_xzz_0, tg_xyyyyyz_yyy_0, \
                                         tg_xyyyyyz_yyz_0, tg_xyyyyyz_yzz_0, tg_xyyyyyz_zzz_0, tg_xyyyyzz_xxx_0, \
                                         tg_xyyyyzz_xxy_0, tg_xyyyyzz_xxz_0, tg_xyyyyzz_xyy_0, tg_xyyyyzz_xyz_0, \
                                         tg_xyyyyzz_xzz_0, tg_xyyyyzz_yyy_0, tg_xyyyyzz_yyz_0, tg_xyyyyzz_yzz_0, \
                                         tg_xyyyyzz_zzz_0, tg_xyyyzzz_xxx_0, tg_xyyyzzz_xxy_0, tg_xyyyzzz_xxz_0, \
                                         tg_xyyyzzz_xyy_0, tg_xyyyzzz_xyz_0, tg_xyyyzzz_xzz_0, tg_xyyyzzz_yyy_0, \
                                         tg_xyyyzzz_yyz_0, tg_xyyyzzz_yzz_0, tg_xyyyzzz_zzz_0, tg_xyyzzz_xx_1, \
                                         tg_xyyzzz_xxx_0, tg_xyyzzz_xxx_1, tg_xyyzzz_xxy_0, tg_xyyzzz_xxy_1, tg_xyyzzz_xxz_0, \
                                         tg_xyyzzz_xxz_1, tg_xyyzzz_xy_1, tg_xyyzzz_xyy_0, tg_xyyzzz_xyy_1, tg_xyyzzz_xyz_0, \
                                         tg_xyyzzz_xyz_1, tg_xyyzzz_xz_1, tg_xyyzzz_xzz_0, tg_xyyzzz_xzz_1, tg_xyyzzz_yy_1, \
                                         tg_xyyzzz_yyy_0, tg_xyyzzz_yyy_1, tg_xyyzzz_yyz_0, tg_xyyzzz_yyz_1, tg_xyyzzz_yz_1, \
                                         tg_xyyzzz_yzz_0, tg_xyyzzz_yzz_1, tg_xyyzzz_zz_1, tg_xyyzzz_zzz_0, tg_xyyzzz_zzz_1, \
                                         tg_xyyzzzz_xxx_0, tg_xyyzzzz_xxy_0, tg_xyyzzzz_xxz_0, tg_xyyzzzz_xyy_0, \
                                         tg_xyyzzzz_xyz_0, tg_xyyzzzz_xzz_0, tg_xyyzzzz_yyy_0, tg_xyyzzzz_yyz_0, \
                                         tg_xyyzzzz_yzz_0, tg_xyyzzzz_zzz_0, tg_xyzzzz_xx_1, tg_xyzzzz_xxx_0, tg_xyzzzz_xxx_1, \
                                         tg_xyzzzz_xxy_0, tg_xyzzzz_xxy_1, tg_xyzzzz_xxz_0, tg_xyzzzz_xxz_1, tg_xyzzzz_xy_1, \
                                         tg_xyzzzz_xyy_0, tg_xyzzzz_xyy_1, tg_xyzzzz_xyz_0, tg_xyzzzz_xyz_1, tg_xyzzzz_xz_1, \
                                         tg_xyzzzz_xzz_0, tg_xyzzzz_xzz_1, tg_xyzzzz_yy_1, tg_xyzzzz_yyy_0, tg_xyzzzz_yyy_1, \
                                         tg_xyzzzz_yyz_0, tg_xyzzzz_yyz_1, tg_xyzzzz_yz_1, tg_xyzzzz_yzz_0, tg_xyzzzz_yzz_1, \
                                         tg_xyzzzz_zz_1, tg_xyzzzz_zzz_0, tg_xyzzzz_zzz_1, tg_xyzzzzz_xxx_0, \
                                         tg_xyzzzzz_xxy_0, tg_xyzzzzz_xxz_0, tg_xyzzzzz_xyy_0, tg_xyzzzzz_xyz_0, \
                                         tg_xyzzzzz_xzz_0, tg_xyzzzzz_yyy_0, tg_xyzzzzz_yyz_0, tg_xyzzzzz_yzz_0, \
                                         tg_xyzzzzz_zzz_0, tg_xzzzzz_xx_1, tg_xzzzzz_xxx_0, tg_xzzzzz_xxx_1, tg_xzzzzz_xxy_0, \
                                         tg_xzzzzz_xxy_1, tg_xzzzzz_xxz_0, tg_xzzzzz_xxz_1, tg_xzzzzz_xy_1, tg_xzzzzz_xyy_0, \
                                         tg_xzzzzz_xyy_1, tg_xzzzzz_xyz_0, tg_xzzzzz_xyz_1, tg_xzzzzz_xz_1, tg_xzzzzz_xzz_0, \
                                         tg_xzzzzz_xzz_1, tg_xzzzzz_yy_1, tg_xzzzzz_yyy_0, tg_xzzzzz_yyy_1, tg_xzzzzz_yyz_0, \
                                         tg_xzzzzz_yyz_1, tg_xzzzzz_yz_1, tg_xzzzzz_yzz_0, tg_xzzzzz_yzz_1, tg_xzzzzz_zz_1, \
                                         tg_xzzzzz_zzz_0, tg_xzzzzz_zzz_1, tg_yyyyyy_xx_1, tg_yyyyyy_xxx_0, tg_yyyyyy_xxx_1, \
                                         tg_yyyyyy_xxy_0, tg_yyyyyy_xxy_1, tg_yyyyyy_xxz_0, tg_yyyyyy_xxz_1, tg_yyyyyy_xy_1, \
                                         tg_yyyyyy_xyy_0, tg_yyyyyy_xyy_1, tg_yyyyyy_xyz_0, tg_yyyyyy_xyz_1, tg_yyyyyy_xz_1, \
                                         tg_yyyyyy_xzz_0, tg_yyyyyy_xzz_1, tg_yyyyyy_yy_1, tg_yyyyyy_yyy_0, tg_yyyyyy_yyy_1, \
                                         tg_yyyyyy_yyz_0, tg_yyyyyy_yyz_1, tg_yyyyyy_yz_1, tg_yyyyyy_yzz_0, tg_yyyyyy_yzz_1, \
                                         tg_yyyyyy_zz_1, tg_yyyyyy_zzz_0, tg_yyyyyy_zzz_1, tg_yyyyyz_xx_1, tg_yyyyyz_xxx_0, \
                                         tg_yyyyyz_xxx_1, tg_yyyyyz_xxy_0, tg_yyyyyz_xxy_1, tg_yyyyyz_xxz_0, tg_yyyyyz_xxz_1, \
                                         tg_yyyyyz_xy_1, tg_yyyyyz_xyy_0, tg_yyyyyz_xyy_1, tg_yyyyyz_xyz_0, tg_yyyyyz_xyz_1, \
                                         tg_yyyyyz_xz_1, tg_yyyyyz_xzz_0, tg_yyyyyz_xzz_1, tg_yyyyyz_yy_1, tg_yyyyyz_yyy_0, \
                                         tg_yyyyyz_yyy_1, tg_yyyyyz_yyz_0, tg_yyyyyz_yyz_1, tg_yyyyyz_yz_1, tg_yyyyyz_yzz_0, \
                                         tg_yyyyyz_yzz_1, tg_yyyyyz_zz_1, tg_yyyyyz_zzz_0, tg_yyyyyz_zzz_1, tg_yyyyzz_xx_1, \
                                         tg_yyyyzz_xxx_0, tg_yyyyzz_xxx_1, tg_yyyyzz_xxy_0, tg_yyyyzz_xxy_1, tg_yyyyzz_xxz_0, \
                                         tg_yyyyzz_xxz_1, tg_yyyyzz_xy_1, tg_yyyyzz_xyy_0, tg_yyyyzz_xyy_1, tg_yyyyzz_xyz_0, \
                                         tg_yyyyzz_xyz_1, tg_yyyyzz_xz_1, tg_yyyyzz_xzz_0, tg_yyyyzz_xzz_1, tg_yyyyzz_yy_1, \
                                         tg_yyyyzz_yyy_0, tg_yyyyzz_yyy_1, tg_yyyyzz_yyz_0, tg_yyyyzz_yyz_1, tg_yyyyzz_yz_1, \
                                         tg_yyyyzz_yzz_0, tg_yyyyzz_yzz_1, tg_yyyyzz_zz_1, tg_yyyyzz_zzz_0, tg_yyyyzz_zzz_1, \
                                         tg_yyyzzz_xx_1, tg_yyyzzz_xxx_0, tg_yyyzzz_xxx_1, tg_yyyzzz_xxy_0, tg_yyyzzz_xxy_1, \
                                         tg_yyyzzz_xxz_0, tg_yyyzzz_xxz_1, tg_yyyzzz_xy_1, tg_yyyzzz_xyy_0, tg_yyyzzz_xyy_1, \
                                         tg_yyyzzz_xyz_0, tg_yyyzzz_xyz_1, tg_yyyzzz_xz_1, tg_yyyzzz_xzz_0, tg_yyyzzz_xzz_1, \
                                         tg_yyyzzz_yy_1, tg_yyyzzz_yyy_0, tg_yyyzzz_yyy_1, tg_yyyzzz_yyz_0, tg_yyyzzz_yyz_1, \
                                         tg_yyyzzz_yz_1, tg_yyyzzz_yzz_0, tg_yyyzzz_yzz_1, tg_yyyzzz_zz_1, tg_yyyzzz_zzz_0, \
                                         tg_yyyzzz_zzz_1, tg_yyzzz_xxx_0, tg_yyzzz_xxx_1, tg_yyzzz_xxy_0, tg_yyzzz_xxy_1, \
                                         tg_yyzzz_xxz_0, tg_yyzzz_xxz_1, tg_yyzzz_xyy_0, tg_yyzzz_xyy_1, tg_yyzzz_xyz_0, \
                                         tg_yyzzz_xyz_1, tg_yyzzz_xzz_0, tg_yyzzz_xzz_1, tg_yyzzz_yyy_0, tg_yyzzz_yyy_1, \
                                         tg_yyzzz_yyz_0, tg_yyzzz_yyz_1, tg_yyzzz_yzz_0, tg_yyzzz_yzz_1, tg_yyzzz_zzz_0, \
                                         tg_yyzzz_zzz_1, tg_yyzzzz_xx_1, tg_yyzzzz_xxx_0, tg_yyzzzz_xxx_1, tg_yyzzzz_xxy_0, \
                                         tg_yyzzzz_xxy_1, tg_yyzzzz_xxz_0, tg_yyzzzz_xxz_1, tg_yyzzzz_xy_1, tg_yyzzzz_xyy_0, \
                                         tg_yyzzzz_xyy_1, tg_yyzzzz_xyz_0, tg_yyzzzz_xyz_1, tg_yyzzzz_xz_1, tg_yyzzzz_xzz_0, \
                                         tg_yyzzzz_xzz_1, tg_yyzzzz_yy_1, tg_yyzzzz_yyy_0, tg_yyzzzz_yyy_1, tg_yyzzzz_yyz_0, \
                                         tg_yyzzzz_yyz_1, tg_yyzzzz_yz_1, tg_yyzzzz_yzz_0, tg_yyzzzz_yzz_1, tg_yyzzzz_zz_1, \
                                         tg_yyzzzz_zzz_0, tg_yyzzzz_zzz_1, tg_yzzzz_xxx_0, tg_yzzzz_xxx_1, tg_yzzzz_xxy_0, \
                                         tg_yzzzz_xxy_1, tg_yzzzz_xxz_0, tg_yzzzz_xxz_1, tg_yzzzz_xyy_0, tg_yzzzz_xyy_1, \
                                         tg_yzzzz_xyz_0, tg_yzzzz_xyz_1, tg_yzzzz_xzz_0, tg_yzzzz_xzz_1, tg_yzzzz_yyy_0, \
                                         tg_yzzzz_yyy_1, tg_yzzzz_yyz_0, tg_yzzzz_yyz_1, tg_yzzzz_yzz_0, tg_yzzzz_yzz_1, \
                                         tg_yzzzz_zzz_0, tg_yzzzz_zzz_1, tg_yzzzzz_xx_1, tg_yzzzzz_xxx_0, tg_yzzzzz_xxx_1, \
                                         tg_yzzzzz_xxy_0, tg_yzzzzz_xxy_1, tg_yzzzzz_xxz_0, tg_yzzzzz_xxz_1, tg_yzzzzz_xy_1, \
                                         tg_yzzzzz_xyy_0, tg_yzzzzz_xyy_1, tg_yzzzzz_xyz_0, tg_yzzzzz_xyz_1, tg_yzzzzz_xz_1, \
                                         tg_yzzzzz_xzz_0, tg_yzzzzz_xzz_1, tg_yzzzzz_yy_1, tg_yzzzzz_yyy_0, tg_yzzzzz_yyy_1, \
                                         tg_yzzzzz_yyz_0, tg_yzzzzz_yyz_1, tg_yzzzzz_yz_1, tg_yzzzzz_yzz_0, tg_yzzzzz_yzz_1, \
                                         tg_yzzzzz_zz_1, tg_yzzzzz_zzz_0, tg_yzzzzz_zzz_1, tg_zzzzz_xxx_0, tg_zzzzz_xxx_1, \
                                         tg_zzzzz_xxy_0, tg_zzzzz_xxy_1, tg_zzzzz_xxz_0, tg_zzzzz_xxz_1, tg_zzzzz_xyy_0, \
                                         tg_zzzzz_xyy_1, tg_zzzzz_xyz_0, tg_zzzzz_xyz_1, tg_zzzzz_xzz_0, tg_zzzzz_xzz_1, \
                                         tg_zzzzz_yyy_0, tg_zzzzz_yyy_1, tg_zzzzz_yyz_0, tg_zzzzz_yyz_1, tg_zzzzz_yzz_0, \
                                         tg_zzzzz_yzz_1, tg_zzzzz_zzz_0, tg_zzzzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxyyzzz_xxx_0[j] = pb_x * tg_xyyzzz_xxx_0[j] + fr * tg_xyyzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxx_0[j] - tg_yyzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzz_xx_1[j];

                    tg_xxyyzzz_xxy_0[j] = pb_x * tg_xyyzzz_xxy_0[j] + fr * tg_xyyzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxy_0[j] - tg_yyzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xy_1[j];

                    tg_xxyyzzz_xxz_0[j] = pb_x * tg_xyyzzz_xxz_0[j] + fr * tg_xyyzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xxz_0[j] - tg_yyzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzz_xz_1[j];

                    tg_xxyyzzz_xyy_0[j] = pb_x * tg_xyyzzz_xyy_0[j] + fr * tg_xyyzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyy_0[j] - tg_yyzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yy_1[j];

                    tg_xxyyzzz_xyz_0[j] = pb_x * tg_xyyzzz_xyz_0[j] + fr * tg_xyyzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xyz_0[j] - tg_yyzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_yz_1[j];

                    tg_xxyyzzz_xzz_0[j] = pb_x * tg_xyyzzz_xzz_0[j] + fr * tg_xyyzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_xzz_0[j] - tg_yyzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzz_zz_1[j];

                    tg_xxyyzzz_yyy_0[j] = pb_x * tg_xyyzzz_yyy_0[j] + fr * tg_xyyzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyy_0[j] - tg_yyzzz_yyy_1[j] * fl1_fza);

                    tg_xxyyzzz_yyz_0[j] = pb_x * tg_xyyzzz_yyz_0[j] + fr * tg_xyyzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yyz_0[j] - tg_yyzzz_yyz_1[j] * fl1_fza);

                    tg_xxyyzzz_yzz_0[j] = pb_x * tg_xyyzzz_yzz_0[j] + fr * tg_xyyzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_yzz_0[j] - tg_yyzzz_yzz_1[j] * fl1_fza);

                    tg_xxyyzzz_zzz_0[j] = pb_x * tg_xyyzzz_zzz_0[j] + fr * tg_xyyzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_yyzzz_zzz_0[j] - tg_yyzzz_zzz_1[j] * fl1_fza);

                    tg_xxyzzzz_xxx_0[j] = pb_x * tg_xyzzzz_xxx_0[j] + fr * tg_xyzzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxx_0[j] - tg_yzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzz_xx_1[j];

                    tg_xxyzzzz_xxy_0[j] = pb_x * tg_xyzzzz_xxy_0[j] + fr * tg_xyzzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxy_0[j] - tg_yzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xy_1[j];

                    tg_xxyzzzz_xxz_0[j] = pb_x * tg_xyzzzz_xxz_0[j] + fr * tg_xyzzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xxz_0[j] - tg_yzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzz_xz_1[j];

                    tg_xxyzzzz_xyy_0[j] = pb_x * tg_xyzzzz_xyy_0[j] + fr * tg_xyzzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyy_0[j] - tg_yzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yy_1[j];

                    tg_xxyzzzz_xyz_0[j] = pb_x * tg_xyzzzz_xyz_0[j] + fr * tg_xyzzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xyz_0[j] - tg_yzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_yz_1[j];

                    tg_xxyzzzz_xzz_0[j] = pb_x * tg_xyzzzz_xzz_0[j] + fr * tg_xyzzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_xzz_0[j] - tg_yzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzz_zz_1[j];

                    tg_xxyzzzz_yyy_0[j] = pb_x * tg_xyzzzz_yyy_0[j] + fr * tg_xyzzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyy_0[j] - tg_yzzzz_yyy_1[j] * fl1_fza);

                    tg_xxyzzzz_yyz_0[j] = pb_x * tg_xyzzzz_yyz_0[j] + fr * tg_xyzzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yyz_0[j] - tg_yzzzz_yyz_1[j] * fl1_fza);

                    tg_xxyzzzz_yzz_0[j] = pb_x * tg_xyzzzz_yzz_0[j] + fr * tg_xyzzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_yzz_0[j] - tg_yzzzz_yzz_1[j] * fl1_fza);

                    tg_xxyzzzz_zzz_0[j] = pb_x * tg_xyzzzz_zzz_0[j] + fr * tg_xyzzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_yzzzz_zzz_0[j] - tg_yzzzz_zzz_1[j] * fl1_fza);

                    tg_xxzzzzz_xxx_0[j] = pb_x * tg_xzzzzz_xxx_0[j] + fr * tg_xzzzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxx_0[j] - tg_zzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzz_xx_1[j];

                    tg_xxzzzzz_xxy_0[j] = pb_x * tg_xzzzzz_xxy_0[j] + fr * tg_xzzzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxy_0[j] - tg_zzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xy_1[j];

                    tg_xxzzzzz_xxz_0[j] = pb_x * tg_xzzzzz_xxz_0[j] + fr * tg_xzzzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxz_0[j] - tg_zzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzz_xz_1[j];

                    tg_xxzzzzz_xyy_0[j] = pb_x * tg_xzzzzz_xyy_0[j] + fr * tg_xzzzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyy_0[j] - tg_zzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yy_1[j];

                    tg_xxzzzzz_xyz_0[j] = pb_x * tg_xzzzzz_xyz_0[j] + fr * tg_xzzzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyz_0[j] - tg_zzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_yz_1[j];

                    tg_xxzzzzz_xzz_0[j] = pb_x * tg_xzzzzz_xzz_0[j] + fr * tg_xzzzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzz_0[j] - tg_zzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzz_zz_1[j];

                    tg_xxzzzzz_yyy_0[j] = pb_x * tg_xzzzzz_yyy_0[j] + fr * tg_xzzzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyy_0[j] - tg_zzzzz_yyy_1[j] * fl1_fza);

                    tg_xxzzzzz_yyz_0[j] = pb_x * tg_xzzzzz_yyz_0[j] + fr * tg_xzzzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyz_0[j] - tg_zzzzz_yyz_1[j] * fl1_fza);

                    tg_xxzzzzz_yzz_0[j] = pb_x * tg_xzzzzz_yzz_0[j] + fr * tg_xzzzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzz_0[j] - tg_zzzzz_yzz_1[j] * fl1_fza);

                    tg_xxzzzzz_zzz_0[j] = pb_x * tg_xzzzzz_zzz_0[j] + fr * tg_xzzzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzz_0[j] - tg_zzzzz_zzz_1[j] * fl1_fza);

                    tg_xyyyyyy_xxx_0[j] = pb_x * tg_yyyyyy_xxx_0[j] + fr * tg_yyyyyy_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyyyy_xx_1[j];

                    tg_xyyyyyy_xxy_0[j] = pb_x * tg_yyyyyy_xxy_0[j] + fr * tg_yyyyyy_xxy_1[j] + fl1_fxn * tg_yyyyyy_xy_1[j];

                    tg_xyyyyyy_xxz_0[j] = pb_x * tg_yyyyyy_xxz_0[j] + fr * tg_yyyyyy_xxz_1[j] + fl1_fxn * tg_yyyyyy_xz_1[j];

                    tg_xyyyyyy_xyy_0[j] = pb_x * tg_yyyyyy_xyy_0[j] + fr * tg_yyyyyy_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yy_1[j];

                    tg_xyyyyyy_xyz_0[j] = pb_x * tg_yyyyyy_xyz_0[j] + fr * tg_yyyyyy_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_yz_1[j];

                    tg_xyyyyyy_xzz_0[j] = pb_x * tg_yyyyyy_xzz_0[j] + fr * tg_yyyyyy_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_zz_1[j];

                    tg_xyyyyyy_yyy_0[j] = pb_x * tg_yyyyyy_yyy_0[j] + fr * tg_yyyyyy_yyy_1[j];

                    tg_xyyyyyy_yyz_0[j] = pb_x * tg_yyyyyy_yyz_0[j] + fr * tg_yyyyyy_yyz_1[j];

                    tg_xyyyyyy_yzz_0[j] = pb_x * tg_yyyyyy_yzz_0[j] + fr * tg_yyyyyy_yzz_1[j];

                    tg_xyyyyyy_zzz_0[j] = pb_x * tg_yyyyyy_zzz_0[j] + fr * tg_yyyyyy_zzz_1[j];

                    tg_xyyyyyz_xxx_0[j] = pb_x * tg_yyyyyz_xxx_0[j] + fr * tg_yyyyyz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyyyz_xx_1[j];

                    tg_xyyyyyz_xxy_0[j] = pb_x * tg_yyyyyz_xxy_0[j] + fr * tg_yyyyyz_xxy_1[j] + fl1_fxn * tg_yyyyyz_xy_1[j];

                    tg_xyyyyyz_xxz_0[j] = pb_x * tg_yyyyyz_xxz_0[j] + fr * tg_yyyyyz_xxz_1[j] + fl1_fxn * tg_yyyyyz_xz_1[j];

                    tg_xyyyyyz_xyy_0[j] = pb_x * tg_yyyyyz_xyy_0[j] + fr * tg_yyyyyz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yy_1[j];

                    tg_xyyyyyz_xyz_0[j] = pb_x * tg_yyyyyz_xyz_0[j] + fr * tg_yyyyyz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_yz_1[j];

                    tg_xyyyyyz_xzz_0[j] = pb_x * tg_yyyyyz_xzz_0[j] + fr * tg_yyyyyz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_zz_1[j];

                    tg_xyyyyyz_yyy_0[j] = pb_x * tg_yyyyyz_yyy_0[j] + fr * tg_yyyyyz_yyy_1[j];

                    tg_xyyyyyz_yyz_0[j] = pb_x * tg_yyyyyz_yyz_0[j] + fr * tg_yyyyyz_yyz_1[j];

                    tg_xyyyyyz_yzz_0[j] = pb_x * tg_yyyyyz_yzz_0[j] + fr * tg_yyyyyz_yzz_1[j];

                    tg_xyyyyyz_zzz_0[j] = pb_x * tg_yyyyyz_zzz_0[j] + fr * tg_yyyyyz_zzz_1[j];

                    tg_xyyyyzz_xxx_0[j] = pb_x * tg_yyyyzz_xxx_0[j] + fr * tg_yyyyzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyyzz_xx_1[j];

                    tg_xyyyyzz_xxy_0[j] = pb_x * tg_yyyyzz_xxy_0[j] + fr * tg_yyyyzz_xxy_1[j] + fl1_fxn * tg_yyyyzz_xy_1[j];

                    tg_xyyyyzz_xxz_0[j] = pb_x * tg_yyyyzz_xxz_0[j] + fr * tg_yyyyzz_xxz_1[j] + fl1_fxn * tg_yyyyzz_xz_1[j];

                    tg_xyyyyzz_xyy_0[j] = pb_x * tg_yyyyzz_xyy_0[j] + fr * tg_yyyyzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yy_1[j];

                    tg_xyyyyzz_xyz_0[j] = pb_x * tg_yyyyzz_xyz_0[j] + fr * tg_yyyyzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_yz_1[j];

                    tg_xyyyyzz_xzz_0[j] = pb_x * tg_yyyyzz_xzz_0[j] + fr * tg_yyyyzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_zz_1[j];

                    tg_xyyyyzz_yyy_0[j] = pb_x * tg_yyyyzz_yyy_0[j] + fr * tg_yyyyzz_yyy_1[j];

                    tg_xyyyyzz_yyz_0[j] = pb_x * tg_yyyyzz_yyz_0[j] + fr * tg_yyyyzz_yyz_1[j];

                    tg_xyyyyzz_yzz_0[j] = pb_x * tg_yyyyzz_yzz_0[j] + fr * tg_yyyyzz_yzz_1[j];

                    tg_xyyyyzz_zzz_0[j] = pb_x * tg_yyyyzz_zzz_0[j] + fr * tg_yyyyzz_zzz_1[j];

                    tg_xyyyzzz_xxx_0[j] = pb_x * tg_yyyzzz_xxx_0[j] + fr * tg_yyyzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyzzz_xx_1[j];

                    tg_xyyyzzz_xxy_0[j] = pb_x * tg_yyyzzz_xxy_0[j] + fr * tg_yyyzzz_xxy_1[j] + fl1_fxn * tg_yyyzzz_xy_1[j];

                    tg_xyyyzzz_xxz_0[j] = pb_x * tg_yyyzzz_xxz_0[j] + fr * tg_yyyzzz_xxz_1[j] + fl1_fxn * tg_yyyzzz_xz_1[j];

                    tg_xyyyzzz_xyy_0[j] = pb_x * tg_yyyzzz_xyy_0[j] + fr * tg_yyyzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yy_1[j];

                    tg_xyyyzzz_xyz_0[j] = pb_x * tg_yyyzzz_xyz_0[j] + fr * tg_yyyzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_yz_1[j];

                    tg_xyyyzzz_xzz_0[j] = pb_x * tg_yyyzzz_xzz_0[j] + fr * tg_yyyzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_zz_1[j];

                    tg_xyyyzzz_yyy_0[j] = pb_x * tg_yyyzzz_yyy_0[j] + fr * tg_yyyzzz_yyy_1[j];

                    tg_xyyyzzz_yyz_0[j] = pb_x * tg_yyyzzz_yyz_0[j] + fr * tg_yyyzzz_yyz_1[j];

                    tg_xyyyzzz_yzz_0[j] = pb_x * tg_yyyzzz_yzz_0[j] + fr * tg_yyyzzz_yzz_1[j];

                    tg_xyyyzzz_zzz_0[j] = pb_x * tg_yyyzzz_zzz_0[j] + fr * tg_yyyzzz_zzz_1[j];

                    tg_xyyzzzz_xxx_0[j] = pb_x * tg_yyzzzz_xxx_0[j] + fr * tg_yyzzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyzzzz_xx_1[j];

                    tg_xyyzzzz_xxy_0[j] = pb_x * tg_yyzzzz_xxy_0[j] + fr * tg_yyzzzz_xxy_1[j] + fl1_fxn * tg_yyzzzz_xy_1[j];

                    tg_xyyzzzz_xxz_0[j] = pb_x * tg_yyzzzz_xxz_0[j] + fr * tg_yyzzzz_xxz_1[j] + fl1_fxn * tg_yyzzzz_xz_1[j];

                    tg_xyyzzzz_xyy_0[j] = pb_x * tg_yyzzzz_xyy_0[j] + fr * tg_yyzzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yy_1[j];

                    tg_xyyzzzz_xyz_0[j] = pb_x * tg_yyzzzz_xyz_0[j] + fr * tg_yyzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_yz_1[j];

                    tg_xyyzzzz_xzz_0[j] = pb_x * tg_yyzzzz_xzz_0[j] + fr * tg_yyzzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_zz_1[j];

                    tg_xyyzzzz_yyy_0[j] = pb_x * tg_yyzzzz_yyy_0[j] + fr * tg_yyzzzz_yyy_1[j];

                    tg_xyyzzzz_yyz_0[j] = pb_x * tg_yyzzzz_yyz_0[j] + fr * tg_yyzzzz_yyz_1[j];

                    tg_xyyzzzz_yzz_0[j] = pb_x * tg_yyzzzz_yzz_0[j] + fr * tg_yyzzzz_yzz_1[j];

                    tg_xyyzzzz_zzz_0[j] = pb_x * tg_yyzzzz_zzz_0[j] + fr * tg_yyzzzz_zzz_1[j];

                    tg_xyzzzzz_xxx_0[j] = pb_x * tg_yzzzzz_xxx_0[j] + fr * tg_yzzzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yzzzzz_xx_1[j];

                    tg_xyzzzzz_xxy_0[j] = pb_x * tg_yzzzzz_xxy_0[j] + fr * tg_yzzzzz_xxy_1[j] + fl1_fxn * tg_yzzzzz_xy_1[j];

                    tg_xyzzzzz_xxz_0[j] = pb_x * tg_yzzzzz_xxz_0[j] + fr * tg_yzzzzz_xxz_1[j] + fl1_fxn * tg_yzzzzz_xz_1[j];

                    tg_xyzzzzz_xyy_0[j] = pb_x * tg_yzzzzz_xyy_0[j] + fr * tg_yzzzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yy_1[j];

                    tg_xyzzzzz_xyz_0[j] = pb_x * tg_yzzzzz_xyz_0[j] + fr * tg_yzzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_yz_1[j];

                    tg_xyzzzzz_xzz_0[j] = pb_x * tg_yzzzzz_xzz_0[j] + fr * tg_yzzzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_zz_1[j];

                    tg_xyzzzzz_yyy_0[j] = pb_x * tg_yzzzzz_yyy_0[j] + fr * tg_yzzzzz_yyy_1[j];

                    tg_xyzzzzz_yyz_0[j] = pb_x * tg_yzzzzz_yyz_0[j] + fr * tg_yzzzzz_yyz_1[j];

                    tg_xyzzzzz_yzz_0[j] = pb_x * tg_yzzzzz_yzz_0[j] + fr * tg_yzzzzz_yzz_1[j];

                    tg_xyzzzzz_zzz_0[j] = pb_x * tg_yzzzzz_zzz_0[j] + fr * tg_yzzzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSF_270_360(      CMemBlock2D<double>* primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {7, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {2, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_5_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_5_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {5, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yyyyyy_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 210); 

                auto tg_yyyyyy_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 211); 

                auto tg_yyyyyy_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 212); 

                auto tg_yyyyyy_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 213); 

                auto tg_yyyyyy_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 214); 

                auto tg_yyyyyy_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 215); 

                auto tg_yyyyyy_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 216); 

                auto tg_yyyyyy_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 217); 

                auto tg_yyyyyy_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 218); 

                auto tg_yyyyyy_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 219); 

                auto tg_yyyyyz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 220); 

                auto tg_yyyyyz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 221); 

                auto tg_yyyyyz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 222); 

                auto tg_yyyyyz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 223); 

                auto tg_yyyyyz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 224); 

                auto tg_yyyyyz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 225); 

                auto tg_yyyyyz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 226); 

                auto tg_yyyyyz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 227); 

                auto tg_yyyyyz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 228); 

                auto tg_yyyyyz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 229); 

                auto tg_yyyyzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 230); 

                auto tg_yyyyzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 231); 

                auto tg_yyyyzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 232); 

                auto tg_yyyyzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 233); 

                auto tg_yyyyzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 234); 

                auto tg_yyyyzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 235); 

                auto tg_yyyyzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 236); 

                auto tg_yyyyzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 237); 

                auto tg_yyyyzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 238); 

                auto tg_yyyyzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 239); 

                auto tg_yyyzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 240); 

                auto tg_yyyzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 241); 

                auto tg_yyyzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 242); 

                auto tg_yyyzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 243); 

                auto tg_yyyzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 244); 

                auto tg_yyyzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 245); 

                auto tg_yyyzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 246); 

                auto tg_yyyzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 247); 

                auto tg_yyyzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 248); 

                auto tg_yyyzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 249); 

                auto tg_yyzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 250); 

                auto tg_yyzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 251); 

                auto tg_yyzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 252); 

                auto tg_yyzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 253); 

                auto tg_yyzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 254); 

                auto tg_yyzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 255); 

                auto tg_yyzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 256); 

                auto tg_yyzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 257); 

                auto tg_yyzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 258); 

                auto tg_yyzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 259); 

                auto tg_yzzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 260); 

                auto tg_yzzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 261); 

                auto tg_yzzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 262); 

                auto tg_yzzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 263); 

                auto tg_yzzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 264); 

                auto tg_yzzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 265); 

                auto tg_yzzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 266); 

                auto tg_yzzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 267); 

                auto tg_yzzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 268); 

                auto tg_yzzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 269); 

                auto tg_zzzzzz_xxx_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 270); 

                auto tg_zzzzzz_xxy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 271); 

                auto tg_zzzzzz_xxz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 272); 

                auto tg_zzzzzz_xyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 273); 

                auto tg_zzzzzz_xyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 274); 

                auto tg_zzzzzz_xzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 275); 

                auto tg_zzzzzz_yyy_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 276); 

                auto tg_zzzzzz_yyz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 277); 

                auto tg_zzzzzz_yzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 278); 

                auto tg_zzzzzz_zzz_0 = primBuffer[pidx_g_6_3_m0].data(280 * idx + 279); 

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

                auto tg_yyyyy_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 150); 

                auto tg_yyyyy_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 151); 

                auto tg_yyyyy_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 152); 

                auto tg_yyyyy_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 153); 

                auto tg_yyyyy_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 154); 

                auto tg_yyyyy_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 155); 

                auto tg_yyyyy_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 156); 

                auto tg_yyyyy_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 157); 

                auto tg_yyyyy_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 158); 

                auto tg_yyyyy_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 159); 

                auto tg_yyyyz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 160); 

                auto tg_yyyyz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 161); 

                auto tg_yyyyz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 162); 

                auto tg_yyyyz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 163); 

                auto tg_yyyyz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 164); 

                auto tg_yyyyz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 165); 

                auto tg_yyyyz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 166); 

                auto tg_yyyyz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 167); 

                auto tg_yyyyz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 168); 

                auto tg_yyyyz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 169); 

                auto tg_yyyzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 170); 

                auto tg_yyyzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 171); 

                auto tg_yyyzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 172); 

                auto tg_yyyzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 173); 

                auto tg_yyyzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 174); 

                auto tg_yyyzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 175); 

                auto tg_yyyzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 176); 

                auto tg_yyyzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 177); 

                auto tg_yyyzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 178); 

                auto tg_yyyzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 179); 

                auto tg_yyzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 180); 

                auto tg_yyzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 181); 

                auto tg_yyzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 182); 

                auto tg_yyzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 183); 

                auto tg_yyzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 184); 

                auto tg_yyzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 185); 

                auto tg_yyzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 186); 

                auto tg_yyzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 187); 

                auto tg_yyzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 188); 

                auto tg_yyzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 189); 

                auto tg_yzzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 190); 

                auto tg_yzzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 191); 

                auto tg_yzzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 192); 

                auto tg_yzzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 193); 

                auto tg_yzzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 194); 

                auto tg_yzzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 195); 

                auto tg_yzzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 196); 

                auto tg_yzzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 197); 

                auto tg_yzzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 198); 

                auto tg_yzzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 199); 

                auto tg_zzzzz_xxx_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 200); 

                auto tg_zzzzz_xxy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 201); 

                auto tg_zzzzz_xxz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 202); 

                auto tg_zzzzz_xyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 203); 

                auto tg_zzzzz_xyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 204); 

                auto tg_zzzzz_xzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 205); 

                auto tg_zzzzz_yyy_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 206); 

                auto tg_zzzzz_yyz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 207); 

                auto tg_zzzzz_yzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 208); 

                auto tg_zzzzz_zzz_0 = primBuffer[pidx_g_5_3_m0].data(210 * idx + 209); 

                auto tg_yyyyy_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 150); 

                auto tg_yyyyy_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 151); 

                auto tg_yyyyy_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 152); 

                auto tg_yyyyy_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 153); 

                auto tg_yyyyy_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 154); 

                auto tg_yyyyy_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 155); 

                auto tg_yyyyy_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 156); 

                auto tg_yyyyy_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 157); 

                auto tg_yyyyy_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 158); 

                auto tg_yyyyy_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 159); 

                auto tg_yyyyz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 160); 

                auto tg_yyyyz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 161); 

                auto tg_yyyyz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 162); 

                auto tg_yyyyz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 163); 

                auto tg_yyyyz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 164); 

                auto tg_yyyyz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 165); 

                auto tg_yyyyz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 166); 

                auto tg_yyyyz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 167); 

                auto tg_yyyyz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 168); 

                auto tg_yyyyz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 169); 

                auto tg_yyyzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 170); 

                auto tg_yyyzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 171); 

                auto tg_yyyzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 172); 

                auto tg_yyyzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 173); 

                auto tg_yyyzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 174); 

                auto tg_yyyzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 175); 

                auto tg_yyyzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 176); 

                auto tg_yyyzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 177); 

                auto tg_yyyzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 178); 

                auto tg_yyyzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 179); 

                auto tg_yyzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 180); 

                auto tg_yyzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 181); 

                auto tg_yyzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 182); 

                auto tg_yyzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 183); 

                auto tg_yyzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 184); 

                auto tg_yyzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 185); 

                auto tg_yyzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 186); 

                auto tg_yyzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 187); 

                auto tg_yyzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 188); 

                auto tg_yyzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 189); 

                auto tg_yzzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 190); 

                auto tg_yzzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 191); 

                auto tg_yzzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 192); 

                auto tg_yzzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 193); 

                auto tg_yzzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 194); 

                auto tg_yzzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 195); 

                auto tg_yzzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 196); 

                auto tg_yzzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 197); 

                auto tg_yzzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 198); 

                auto tg_yzzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 199); 

                auto tg_zzzzz_xxx_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 200); 

                auto tg_zzzzz_xxy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 201); 

                auto tg_zzzzz_xxz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 202); 

                auto tg_zzzzz_xyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 203); 

                auto tg_zzzzz_xyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 204); 

                auto tg_zzzzz_xzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 205); 

                auto tg_zzzzz_yyy_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 206); 

                auto tg_zzzzz_yyz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 207); 

                auto tg_zzzzz_yzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 208); 

                auto tg_zzzzz_zzz_1 = primBuffer[pidx_g_5_3_m1].data(210 * idx + 209); 

                auto tg_yyyyyy_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 126); 

                auto tg_yyyyyy_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 127); 

                auto tg_yyyyyy_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 128); 

                auto tg_yyyyyy_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 129); 

                auto tg_yyyyyy_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 130); 

                auto tg_yyyyyy_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 131); 

                auto tg_yyyyyz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 132); 

                auto tg_yyyyyz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 133); 

                auto tg_yyyyyz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 134); 

                auto tg_yyyyyz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 135); 

                auto tg_yyyyyz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 136); 

                auto tg_yyyyyz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 137); 

                auto tg_yyyyzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 138); 

                auto tg_yyyyzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 139); 

                auto tg_yyyyzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 140); 

                auto tg_yyyyzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 141); 

                auto tg_yyyyzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 142); 

                auto tg_yyyyzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 143); 

                auto tg_yyyzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 144); 

                auto tg_yyyzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 145); 

                auto tg_yyyzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 146); 

                auto tg_yyyzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 147); 

                auto tg_yyyzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 148); 

                auto tg_yyyzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 149); 

                auto tg_yyzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 150); 

                auto tg_yyzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 151); 

                auto tg_yyzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 152); 

                auto tg_yyzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 153); 

                auto tg_yyzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 154); 

                auto tg_yyzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 155); 

                auto tg_yzzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 156); 

                auto tg_yzzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 157); 

                auto tg_yzzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 158); 

                auto tg_yzzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 159); 

                auto tg_yzzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 160); 

                auto tg_yzzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 161); 

                auto tg_zzzzzz_xx_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 162); 

                auto tg_zzzzzz_xy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 163); 

                auto tg_zzzzzz_xz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 164); 

                auto tg_zzzzzz_yy_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 165); 

                auto tg_zzzzzz_yz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 166); 

                auto tg_zzzzzz_zz_1 = primBuffer[pidx_g_6_2_m1].data(168 * idx + 167); 

                // set up pointers to integrals

                auto tg_xzzzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 270); 

                auto tg_xzzzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 271); 

                auto tg_xzzzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 272); 

                auto tg_xzzzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 273); 

                auto tg_xzzzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 274); 

                auto tg_xzzzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 275); 

                auto tg_xzzzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 276); 

                auto tg_xzzzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 277); 

                auto tg_xzzzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 278); 

                auto tg_xzzzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 279); 

                auto tg_yyyyyyy_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 280); 

                auto tg_yyyyyyy_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 281); 

                auto tg_yyyyyyy_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 282); 

                auto tg_yyyyyyy_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 283); 

                auto tg_yyyyyyy_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 284); 

                auto tg_yyyyyyy_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 285); 

                auto tg_yyyyyyy_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 286); 

                auto tg_yyyyyyy_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 287); 

                auto tg_yyyyyyy_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 288); 

                auto tg_yyyyyyy_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 289); 

                auto tg_yyyyyyz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 290); 

                auto tg_yyyyyyz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 291); 

                auto tg_yyyyyyz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 292); 

                auto tg_yyyyyyz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 293); 

                auto tg_yyyyyyz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 294); 

                auto tg_yyyyyyz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 295); 

                auto tg_yyyyyyz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 296); 

                auto tg_yyyyyyz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 297); 

                auto tg_yyyyyyz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 298); 

                auto tg_yyyyyyz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 299); 

                auto tg_yyyyyzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 300); 

                auto tg_yyyyyzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 301); 

                auto tg_yyyyyzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 302); 

                auto tg_yyyyyzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 303); 

                auto tg_yyyyyzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 304); 

                auto tg_yyyyyzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 305); 

                auto tg_yyyyyzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 306); 

                auto tg_yyyyyzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 307); 

                auto tg_yyyyyzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 308); 

                auto tg_yyyyyzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 309); 

                auto tg_yyyyzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 310); 

                auto tg_yyyyzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 311); 

                auto tg_yyyyzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 312); 

                auto tg_yyyyzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 313); 

                auto tg_yyyyzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 314); 

                auto tg_yyyyzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 315); 

                auto tg_yyyyzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 316); 

                auto tg_yyyyzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 317); 

                auto tg_yyyyzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 318); 

                auto tg_yyyyzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 319); 

                auto tg_yyyzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 320); 

                auto tg_yyyzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 321); 

                auto tg_yyyzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 322); 

                auto tg_yyyzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 323); 

                auto tg_yyyzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 324); 

                auto tg_yyyzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 325); 

                auto tg_yyyzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 326); 

                auto tg_yyyzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 327); 

                auto tg_yyyzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 328); 

                auto tg_yyyzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 329); 

                auto tg_yyzzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 330); 

                auto tg_yyzzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 331); 

                auto tg_yyzzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 332); 

                auto tg_yyzzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 333); 

                auto tg_yyzzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 334); 

                auto tg_yyzzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 335); 

                auto tg_yyzzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 336); 

                auto tg_yyzzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 337); 

                auto tg_yyzzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 338); 

                auto tg_yyzzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 339); 

                auto tg_yzzzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 340); 

                auto tg_yzzzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 341); 

                auto tg_yzzzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 342); 

                auto tg_yzzzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 343); 

                auto tg_yzzzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 344); 

                auto tg_yzzzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 345); 

                auto tg_yzzzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 346); 

                auto tg_yzzzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 347); 

                auto tg_yzzzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 348); 

                auto tg_yzzzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 349); 

                auto tg_zzzzzzz_xxx_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 350); 

                auto tg_zzzzzzz_xxy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 351); 

                auto tg_zzzzzzz_xxz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 352); 

                auto tg_zzzzzzz_xyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 353); 

                auto tg_zzzzzzz_xyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 354); 

                auto tg_zzzzzzz_xzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 355); 

                auto tg_zzzzzzz_yyy_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 356); 

                auto tg_zzzzzzz_yyz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 357); 

                auto tg_zzzzzzz_yzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 358); 

                auto tg_zzzzzzz_zzz_0 = primBuffer[pidx_g_7_3_m0].data(360 * idx + 359); 

                // Batch of Integrals (270,360)

                #pragma omp simd aligned(fxn, fza, tg_xzzzzzz_xxx_0, tg_xzzzzzz_xxy_0, tg_xzzzzzz_xxz_0, \
                                         tg_xzzzzzz_xyy_0, tg_xzzzzzz_xyz_0, tg_xzzzzzz_xzz_0, tg_xzzzzzz_yyy_0, \
                                         tg_xzzzzzz_yyz_0, tg_xzzzzzz_yzz_0, tg_xzzzzzz_zzz_0, tg_yyyyy_xxx_0, tg_yyyyy_xxx_1, \
                                         tg_yyyyy_xxy_0, tg_yyyyy_xxy_1, tg_yyyyy_xxz_0, tg_yyyyy_xxz_1, tg_yyyyy_xyy_0, \
                                         tg_yyyyy_xyy_1, tg_yyyyy_xyz_0, tg_yyyyy_xyz_1, tg_yyyyy_xzz_0, tg_yyyyy_xzz_1, \
                                         tg_yyyyy_yyy_0, tg_yyyyy_yyy_1, tg_yyyyy_yyz_0, tg_yyyyy_yyz_1, tg_yyyyy_yzz_0, \
                                         tg_yyyyy_yzz_1, tg_yyyyy_zzz_0, tg_yyyyy_zzz_1, tg_yyyyyy_xx_1, tg_yyyyyy_xxx_0, \
                                         tg_yyyyyy_xxx_1, tg_yyyyyy_xxy_0, tg_yyyyyy_xxy_1, tg_yyyyyy_xxz_0, tg_yyyyyy_xxz_1, \
                                         tg_yyyyyy_xy_1, tg_yyyyyy_xyy_0, tg_yyyyyy_xyy_1, tg_yyyyyy_xyz_0, tg_yyyyyy_xyz_1, \
                                         tg_yyyyyy_xz_1, tg_yyyyyy_xzz_0, tg_yyyyyy_xzz_1, tg_yyyyyy_yy_1, tg_yyyyyy_yyy_0, \
                                         tg_yyyyyy_yyy_1, tg_yyyyyy_yyz_0, tg_yyyyyy_yyz_1, tg_yyyyyy_yz_1, tg_yyyyyy_yzz_0, \
                                         tg_yyyyyy_yzz_1, tg_yyyyyy_zz_1, tg_yyyyyy_zzz_0, tg_yyyyyy_zzz_1, tg_yyyyyyy_xxx_0, \
                                         tg_yyyyyyy_xxy_0, tg_yyyyyyy_xxz_0, tg_yyyyyyy_xyy_0, tg_yyyyyyy_xyz_0, \
                                         tg_yyyyyyy_xzz_0, tg_yyyyyyy_yyy_0, tg_yyyyyyy_yyz_0, tg_yyyyyyy_yzz_0, \
                                         tg_yyyyyyy_zzz_0, tg_yyyyyyz_xxx_0, tg_yyyyyyz_xxy_0, tg_yyyyyyz_xxz_0, \
                                         tg_yyyyyyz_xyy_0, tg_yyyyyyz_xyz_0, tg_yyyyyyz_xzz_0, tg_yyyyyyz_yyy_0, \
                                         tg_yyyyyyz_yyz_0, tg_yyyyyyz_yzz_0, tg_yyyyyyz_zzz_0, tg_yyyyyz_xx_1, \
                                         tg_yyyyyz_xxx_0, tg_yyyyyz_xxx_1, tg_yyyyyz_xxy_0, tg_yyyyyz_xxy_1, tg_yyyyyz_xxz_0, \
                                         tg_yyyyyz_xxz_1, tg_yyyyyz_xy_1, tg_yyyyyz_xyy_0, tg_yyyyyz_xyy_1, tg_yyyyyz_xyz_0, \
                                         tg_yyyyyz_xyz_1, tg_yyyyyz_xz_1, tg_yyyyyz_xzz_0, tg_yyyyyz_xzz_1, tg_yyyyyz_yy_1, \
                                         tg_yyyyyz_yyy_0, tg_yyyyyz_yyy_1, tg_yyyyyz_yyz_0, tg_yyyyyz_yyz_1, tg_yyyyyz_yz_1, \
                                         tg_yyyyyz_yzz_0, tg_yyyyyz_yzz_1, tg_yyyyyz_zz_1, tg_yyyyyz_zzz_0, tg_yyyyyz_zzz_1, \
                                         tg_yyyyyzz_xxx_0, tg_yyyyyzz_xxy_0, tg_yyyyyzz_xxz_0, tg_yyyyyzz_xyy_0, \
                                         tg_yyyyyzz_xyz_0, tg_yyyyyzz_xzz_0, tg_yyyyyzz_yyy_0, tg_yyyyyzz_yyz_0, \
                                         tg_yyyyyzz_yzz_0, tg_yyyyyzz_zzz_0, tg_yyyyz_xxx_0, tg_yyyyz_xxx_1, tg_yyyyz_xxy_0, \
                                         tg_yyyyz_xxy_1, tg_yyyyz_xxz_0, tg_yyyyz_xxz_1, tg_yyyyz_xyy_0, tg_yyyyz_xyy_1, \
                                         tg_yyyyz_xyz_0, tg_yyyyz_xyz_1, tg_yyyyz_xzz_0, tg_yyyyz_xzz_1, tg_yyyyz_yyy_0, \
                                         tg_yyyyz_yyy_1, tg_yyyyz_yyz_0, tg_yyyyz_yyz_1, tg_yyyyz_yzz_0, tg_yyyyz_yzz_1, \
                                         tg_yyyyz_zzz_0, tg_yyyyz_zzz_1, tg_yyyyzz_xx_1, tg_yyyyzz_xxx_0, tg_yyyyzz_xxx_1, \
                                         tg_yyyyzz_xxy_0, tg_yyyyzz_xxy_1, tg_yyyyzz_xxz_0, tg_yyyyzz_xxz_1, tg_yyyyzz_xy_1, \
                                         tg_yyyyzz_xyy_0, tg_yyyyzz_xyy_1, tg_yyyyzz_xyz_0, tg_yyyyzz_xyz_1, tg_yyyyzz_xz_1, \
                                         tg_yyyyzz_xzz_0, tg_yyyyzz_xzz_1, tg_yyyyzz_yy_1, tg_yyyyzz_yyy_0, tg_yyyyzz_yyy_1, \
                                         tg_yyyyzz_yyz_0, tg_yyyyzz_yyz_1, tg_yyyyzz_yz_1, tg_yyyyzz_yzz_0, tg_yyyyzz_yzz_1, \
                                         tg_yyyyzz_zz_1, tg_yyyyzz_zzz_0, tg_yyyyzz_zzz_1, tg_yyyyzzz_xxx_0, \
                                         tg_yyyyzzz_xxy_0, tg_yyyyzzz_xxz_0, tg_yyyyzzz_xyy_0, tg_yyyyzzz_xyz_0, \
                                         tg_yyyyzzz_xzz_0, tg_yyyyzzz_yyy_0, tg_yyyyzzz_yyz_0, tg_yyyyzzz_yzz_0, \
                                         tg_yyyyzzz_zzz_0, tg_yyyzz_xxx_0, tg_yyyzz_xxx_1, tg_yyyzz_xxy_0, tg_yyyzz_xxy_1, \
                                         tg_yyyzz_xxz_0, tg_yyyzz_xxz_1, tg_yyyzz_xyy_0, tg_yyyzz_xyy_1, tg_yyyzz_xyz_0, \
                                         tg_yyyzz_xyz_1, tg_yyyzz_xzz_0, tg_yyyzz_xzz_1, tg_yyyzz_yyy_0, tg_yyyzz_yyy_1, \
                                         tg_yyyzz_yyz_0, tg_yyyzz_yyz_1, tg_yyyzz_yzz_0, tg_yyyzz_yzz_1, tg_yyyzz_zzz_0, \
                                         tg_yyyzz_zzz_1, tg_yyyzzz_xx_1, tg_yyyzzz_xxx_0, tg_yyyzzz_xxx_1, tg_yyyzzz_xxy_0, \
                                         tg_yyyzzz_xxy_1, tg_yyyzzz_xxz_0, tg_yyyzzz_xxz_1, tg_yyyzzz_xy_1, tg_yyyzzz_xyy_0, \
                                         tg_yyyzzz_xyy_1, tg_yyyzzz_xyz_0, tg_yyyzzz_xyz_1, tg_yyyzzz_xz_1, tg_yyyzzz_xzz_0, \
                                         tg_yyyzzz_xzz_1, tg_yyyzzz_yy_1, tg_yyyzzz_yyy_0, tg_yyyzzz_yyy_1, tg_yyyzzz_yyz_0, \
                                         tg_yyyzzz_yyz_1, tg_yyyzzz_yz_1, tg_yyyzzz_yzz_0, tg_yyyzzz_yzz_1, tg_yyyzzz_zz_1, \
                                         tg_yyyzzz_zzz_0, tg_yyyzzz_zzz_1, tg_yyyzzzz_xxx_0, tg_yyyzzzz_xxy_0, \
                                         tg_yyyzzzz_xxz_0, tg_yyyzzzz_xyy_0, tg_yyyzzzz_xyz_0, tg_yyyzzzz_xzz_0, \
                                         tg_yyyzzzz_yyy_0, tg_yyyzzzz_yyz_0, tg_yyyzzzz_yzz_0, tg_yyyzzzz_zzz_0, \
                                         tg_yyzzz_xxx_0, tg_yyzzz_xxx_1, tg_yyzzz_xxy_0, tg_yyzzz_xxy_1, tg_yyzzz_xxz_0, \
                                         tg_yyzzz_xxz_1, tg_yyzzz_xyy_0, tg_yyzzz_xyy_1, tg_yyzzz_xyz_0, tg_yyzzz_xyz_1, \
                                         tg_yyzzz_xzz_0, tg_yyzzz_xzz_1, tg_yyzzz_yyy_0, tg_yyzzz_yyy_1, tg_yyzzz_yyz_0, \
                                         tg_yyzzz_yyz_1, tg_yyzzz_yzz_0, tg_yyzzz_yzz_1, tg_yyzzz_zzz_0, tg_yyzzz_zzz_1, \
                                         tg_yyzzzz_xx_1, tg_yyzzzz_xxx_0, tg_yyzzzz_xxx_1, tg_yyzzzz_xxy_0, tg_yyzzzz_xxy_1, \
                                         tg_yyzzzz_xxz_0, tg_yyzzzz_xxz_1, tg_yyzzzz_xy_1, tg_yyzzzz_xyy_0, tg_yyzzzz_xyy_1, \
                                         tg_yyzzzz_xyz_0, tg_yyzzzz_xyz_1, tg_yyzzzz_xz_1, tg_yyzzzz_xzz_0, tg_yyzzzz_xzz_1, \
                                         tg_yyzzzz_yy_1, tg_yyzzzz_yyy_0, tg_yyzzzz_yyy_1, tg_yyzzzz_yyz_0, tg_yyzzzz_yyz_1, \
                                         tg_yyzzzz_yz_1, tg_yyzzzz_yzz_0, tg_yyzzzz_yzz_1, tg_yyzzzz_zz_1, tg_yyzzzz_zzz_0, \
                                         tg_yyzzzz_zzz_1, tg_yyzzzzz_xxx_0, tg_yyzzzzz_xxy_0, tg_yyzzzzz_xxz_0, \
                                         tg_yyzzzzz_xyy_0, tg_yyzzzzz_xyz_0, tg_yyzzzzz_xzz_0, tg_yyzzzzz_yyy_0, \
                                         tg_yyzzzzz_yyz_0, tg_yyzzzzz_yzz_0, tg_yyzzzzz_zzz_0, tg_yzzzz_xxx_0, tg_yzzzz_xxx_1, \
                                         tg_yzzzz_xxy_0, tg_yzzzz_xxy_1, tg_yzzzz_xxz_0, tg_yzzzz_xxz_1, tg_yzzzz_xyy_0, \
                                         tg_yzzzz_xyy_1, tg_yzzzz_xyz_0, tg_yzzzz_xyz_1, tg_yzzzz_xzz_0, tg_yzzzz_xzz_1, \
                                         tg_yzzzz_yyy_0, tg_yzzzz_yyy_1, tg_yzzzz_yyz_0, tg_yzzzz_yyz_1, tg_yzzzz_yzz_0, \
                                         tg_yzzzz_yzz_1, tg_yzzzz_zzz_0, tg_yzzzz_zzz_1, tg_yzzzzz_xx_1, tg_yzzzzz_xxx_0, \
                                         tg_yzzzzz_xxx_1, tg_yzzzzz_xxy_0, tg_yzzzzz_xxy_1, tg_yzzzzz_xxz_0, tg_yzzzzz_xxz_1, \
                                         tg_yzzzzz_xy_1, tg_yzzzzz_xyy_0, tg_yzzzzz_xyy_1, tg_yzzzzz_xyz_0, tg_yzzzzz_xyz_1, \
                                         tg_yzzzzz_xz_1, tg_yzzzzz_xzz_0, tg_yzzzzz_xzz_1, tg_yzzzzz_yy_1, tg_yzzzzz_yyy_0, \
                                         tg_yzzzzz_yyy_1, tg_yzzzzz_yyz_0, tg_yzzzzz_yyz_1, tg_yzzzzz_yz_1, tg_yzzzzz_yzz_0, \
                                         tg_yzzzzz_yzz_1, tg_yzzzzz_zz_1, tg_yzzzzz_zzz_0, tg_yzzzzz_zzz_1, tg_yzzzzzz_xxx_0, \
                                         tg_yzzzzzz_xxy_0, tg_yzzzzzz_xxz_0, tg_yzzzzzz_xyy_0, tg_yzzzzzz_xyz_0, \
                                         tg_yzzzzzz_xzz_0, tg_yzzzzzz_yyy_0, tg_yzzzzzz_yyz_0, tg_yzzzzzz_yzz_0, \
                                         tg_yzzzzzz_zzz_0, tg_zzzzz_xxx_0, tg_zzzzz_xxx_1, tg_zzzzz_xxy_0, tg_zzzzz_xxy_1, \
                                         tg_zzzzz_xxz_0, tg_zzzzz_xxz_1, tg_zzzzz_xyy_0, tg_zzzzz_xyy_1, tg_zzzzz_xyz_0, \
                                         tg_zzzzz_xyz_1, tg_zzzzz_xzz_0, tg_zzzzz_xzz_1, tg_zzzzz_yyy_0, tg_zzzzz_yyy_1, \
                                         tg_zzzzz_yyz_0, tg_zzzzz_yyz_1, tg_zzzzz_yzz_0, tg_zzzzz_yzz_1, tg_zzzzz_zzz_0, \
                                         tg_zzzzz_zzz_1, tg_zzzzzz_xx_1, tg_zzzzzz_xxx_0, tg_zzzzzz_xxx_1, tg_zzzzzz_xxy_0, \
                                         tg_zzzzzz_xxy_1, tg_zzzzzz_xxz_0, tg_zzzzzz_xxz_1, tg_zzzzzz_xy_1, tg_zzzzzz_xyy_0, \
                                         tg_zzzzzz_xyy_1, tg_zzzzzz_xyz_0, tg_zzzzzz_xyz_1, tg_zzzzzz_xz_1, tg_zzzzzz_xzz_0, \
                                         tg_zzzzzz_xzz_1, tg_zzzzzz_yy_1, tg_zzzzzz_yyy_0, tg_zzzzzz_yyy_1, tg_zzzzzz_yyz_0, \
                                         tg_zzzzzz_yyz_1, tg_zzzzzz_yz_1, tg_zzzzzz_yzz_0, tg_zzzzzz_yzz_1, tg_zzzzzz_zz_1, \
                                         tg_zzzzzz_zzz_0, tg_zzzzzz_zzz_1, tg_zzzzzzz_xxx_0, tg_zzzzzzz_xxy_0, \
                                         tg_zzzzzzz_xxz_0, tg_zzzzzzz_xyy_0, tg_zzzzzzz_xyz_0, tg_zzzzzzz_xzz_0, \
                                         tg_zzzzzzz_yyy_0, tg_zzzzzzz_yyz_0, tg_zzzzzzz_yzz_0, tg_zzzzzzz_zzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xzzzzzz_xxx_0[j] = pb_x * tg_zzzzzz_xxx_0[j] + fr * tg_zzzzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_xx_1[j];

                    tg_xzzzzzz_xxy_0[j] = pb_x * tg_zzzzzz_xxy_0[j] + fr * tg_zzzzzz_xxy_1[j] + fl1_fxn * tg_zzzzzz_xy_1[j];

                    tg_xzzzzzz_xxz_0[j] = pb_x * tg_zzzzzz_xxz_0[j] + fr * tg_zzzzzz_xxz_1[j] + fl1_fxn * tg_zzzzzz_xz_1[j];

                    tg_xzzzzzz_xyy_0[j] = pb_x * tg_zzzzzz_xyy_0[j] + fr * tg_zzzzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yy_1[j];

                    tg_xzzzzzz_xyz_0[j] = pb_x * tg_zzzzzz_xyz_0[j] + fr * tg_zzzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_yz_1[j];

                    tg_xzzzzzz_xzz_0[j] = pb_x * tg_zzzzzz_xzz_0[j] + fr * tg_zzzzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zz_1[j];

                    tg_xzzzzzz_yyy_0[j] = pb_x * tg_zzzzzz_yyy_0[j] + fr * tg_zzzzzz_yyy_1[j];

                    tg_xzzzzzz_yyz_0[j] = pb_x * tg_zzzzzz_yyz_0[j] + fr * tg_zzzzzz_yyz_1[j];

                    tg_xzzzzzz_yzz_0[j] = pb_x * tg_zzzzzz_yzz_0[j] + fr * tg_zzzzzz_yzz_1[j];

                    tg_xzzzzzz_zzz_0[j] = pb_x * tg_zzzzzz_zzz_0[j] + fr * tg_zzzzzz_zzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyyyyy_xxx_0[j] = pb_y * tg_yyyyyy_xxx_0[j] + fr * tg_yyyyyy_xxx_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxx_0[j] - tg_yyyyy_xxx_1[j] * fl1_fza);

                    tg_yyyyyyy_xxy_0[j] = pb_y * tg_yyyyyy_xxy_0[j] + fr * tg_yyyyyy_xxy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxy_0[j] - tg_yyyyy_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xx_1[j];

                    tg_yyyyyyy_xxz_0[j] = pb_y * tg_yyyyyy_xxz_0[j] + fr * tg_yyyyyy_xxz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xxz_0[j] - tg_yyyyy_xxz_1[j] * fl1_fza);

                    tg_yyyyyyy_xyy_0[j] = pb_y * tg_yyyyyy_xyy_0[j] + fr * tg_yyyyyy_xyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyy_0[j] - tg_yyyyy_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_xy_1[j];

                    tg_yyyyyyy_xyz_0[j] = pb_y * tg_yyyyyy_xyz_0[j] + fr * tg_yyyyyy_xyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xyz_0[j] - tg_yyyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_xz_1[j];

                    tg_yyyyyyy_xzz_0[j] = pb_y * tg_yyyyyy_xzz_0[j] + fr * tg_yyyyyy_xzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_xzz_0[j] - tg_yyyyy_xzz_1[j] * fl1_fza);

                    tg_yyyyyyy_yyy_0[j] = pb_y * tg_yyyyyy_yyy_0[j] + fr * tg_yyyyyy_yyy_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyy_0[j] - tg_yyyyy_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyy_yy_1[j];

                    tg_yyyyyyy_yyz_0[j] = pb_y * tg_yyyyyy_yyz_0[j] + fr * tg_yyyyyy_yyz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yyz_0[j] - tg_yyyyy_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyy_yz_1[j];

                    tg_yyyyyyy_yzz_0[j] = pb_y * tg_yyyyyy_yzz_0[j] + fr * tg_yyyyyy_yzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_yzz_0[j] - tg_yyyyy_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyy_zz_1[j];

                    tg_yyyyyyy_zzz_0[j] = pb_y * tg_yyyyyy_zzz_0[j] + fr * tg_yyyyyy_zzz_1[j] + 3.0 * fl1_fx * (tg_yyyyy_zzz_0[j] - tg_yyyyy_zzz_1[j] * fl1_fza);

                    tg_yyyyyyz_xxx_0[j] = pb_y * tg_yyyyyz_xxx_0[j] + fr * tg_yyyyyz_xxx_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxx_0[j] - tg_yyyyz_xxx_1[j] * fl1_fza);

                    tg_yyyyyyz_xxy_0[j] = pb_y * tg_yyyyyz_xxy_0[j] + fr * tg_yyyyyz_xxy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxy_0[j] - tg_yyyyz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xx_1[j];

                    tg_yyyyyyz_xxz_0[j] = pb_y * tg_yyyyyz_xxz_0[j] + fr * tg_yyyyyz_xxz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xxz_0[j] - tg_yyyyz_xxz_1[j] * fl1_fza);

                    tg_yyyyyyz_xyy_0[j] = pb_y * tg_yyyyyz_xyy_0[j] + fr * tg_yyyyyz_xyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyy_0[j] - tg_yyyyz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_xy_1[j];

                    tg_yyyyyyz_xyz_0[j] = pb_y * tg_yyyyyz_xyz_0[j] + fr * tg_yyyyyz_xyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xyz_0[j] - tg_yyyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_xz_1[j];

                    tg_yyyyyyz_xzz_0[j] = pb_y * tg_yyyyyz_xzz_0[j] + fr * tg_yyyyyz_xzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_xzz_0[j] - tg_yyyyz_xzz_1[j] * fl1_fza);

                    tg_yyyyyyz_yyy_0[j] = pb_y * tg_yyyyyz_yyy_0[j] + fr * tg_yyyyyz_yyy_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyy_0[j] - tg_yyyyz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyz_yy_1[j];

                    tg_yyyyyyz_yyz_0[j] = pb_y * tg_yyyyyz_yyz_0[j] + fr * tg_yyyyyz_yyz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yyz_0[j] - tg_yyyyz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyz_yz_1[j];

                    tg_yyyyyyz_yzz_0[j] = pb_y * tg_yyyyyz_yzz_0[j] + fr * tg_yyyyyz_yzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_yzz_0[j] - tg_yyyyz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyz_zz_1[j];

                    tg_yyyyyyz_zzz_0[j] = pb_y * tg_yyyyyz_zzz_0[j] + fr * tg_yyyyyz_zzz_1[j] + 2.5 * fl1_fx * (tg_yyyyz_zzz_0[j] - tg_yyyyz_zzz_1[j] * fl1_fza);

                    tg_yyyyyzz_xxx_0[j] = pb_y * tg_yyyyzz_xxx_0[j] + fr * tg_yyyyzz_xxx_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxx_0[j] - tg_yyyzz_xxx_1[j] * fl1_fza);

                    tg_yyyyyzz_xxy_0[j] = pb_y * tg_yyyyzz_xxy_0[j] + fr * tg_yyyyzz_xxy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxy_0[j] - tg_yyyzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xx_1[j];

                    tg_yyyyyzz_xxz_0[j] = pb_y * tg_yyyyzz_xxz_0[j] + fr * tg_yyyyzz_xxz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xxz_0[j] - tg_yyyzz_xxz_1[j] * fl1_fza);

                    tg_yyyyyzz_xyy_0[j] = pb_y * tg_yyyyzz_xyy_0[j] + fr * tg_yyyyzz_xyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyy_0[j] - tg_yyyzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_xy_1[j];

                    tg_yyyyyzz_xyz_0[j] = pb_y * tg_yyyyzz_xyz_0[j] + fr * tg_yyyyzz_xyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xyz_0[j] - tg_yyyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_xz_1[j];

                    tg_yyyyyzz_xzz_0[j] = pb_y * tg_yyyyzz_xzz_0[j] + fr * tg_yyyyzz_xzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_xzz_0[j] - tg_yyyzz_xzz_1[j] * fl1_fza);

                    tg_yyyyyzz_yyy_0[j] = pb_y * tg_yyyyzz_yyy_0[j] + fr * tg_yyyyzz_yyy_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyy_0[j] - tg_yyyzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzz_yy_1[j];

                    tg_yyyyyzz_yyz_0[j] = pb_y * tg_yyyyzz_yyz_0[j] + fr * tg_yyyyzz_yyz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yyz_0[j] - tg_yyyzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzz_yz_1[j];

                    tg_yyyyyzz_yzz_0[j] = pb_y * tg_yyyyzz_yzz_0[j] + fr * tg_yyyyzz_yzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_yzz_0[j] - tg_yyyzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzz_zz_1[j];

                    tg_yyyyyzz_zzz_0[j] = pb_y * tg_yyyyzz_zzz_0[j] + fr * tg_yyyyzz_zzz_1[j] + 2.0 * fl1_fx * (tg_yyyzz_zzz_0[j] - tg_yyyzz_zzz_1[j] * fl1_fza);

                    tg_yyyyzzz_xxx_0[j] = pb_y * tg_yyyzzz_xxx_0[j] + fr * tg_yyyzzz_xxx_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxx_0[j] - tg_yyzzz_xxx_1[j] * fl1_fza);

                    tg_yyyyzzz_xxy_0[j] = pb_y * tg_yyyzzz_xxy_0[j] + fr * tg_yyyzzz_xxy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxy_0[j] - tg_yyzzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xx_1[j];

                    tg_yyyyzzz_xxz_0[j] = pb_y * tg_yyyzzz_xxz_0[j] + fr * tg_yyyzzz_xxz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xxz_0[j] - tg_yyzzz_xxz_1[j] * fl1_fza);

                    tg_yyyyzzz_xyy_0[j] = pb_y * tg_yyyzzz_xyy_0[j] + fr * tg_yyyzzz_xyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyy_0[j] - tg_yyzzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_xy_1[j];

                    tg_yyyyzzz_xyz_0[j] = pb_y * tg_yyyzzz_xyz_0[j] + fr * tg_yyyzzz_xyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xyz_0[j] - tg_yyzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_xz_1[j];

                    tg_yyyyzzz_xzz_0[j] = pb_y * tg_yyyzzz_xzz_0[j] + fr * tg_yyyzzz_xzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_xzz_0[j] - tg_yyzzz_xzz_1[j] * fl1_fza);

                    tg_yyyyzzz_yyy_0[j] = pb_y * tg_yyyzzz_yyy_0[j] + fr * tg_yyyzzz_yyy_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyy_0[j] - tg_yyzzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzz_yy_1[j];

                    tg_yyyyzzz_yyz_0[j] = pb_y * tg_yyyzzz_yyz_0[j] + fr * tg_yyyzzz_yyz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yyz_0[j] - tg_yyzzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzz_yz_1[j];

                    tg_yyyyzzz_yzz_0[j] = pb_y * tg_yyyzzz_yzz_0[j] + fr * tg_yyyzzz_yzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_yzz_0[j] - tg_yyzzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzz_zz_1[j];

                    tg_yyyyzzz_zzz_0[j] = pb_y * tg_yyyzzz_zzz_0[j] + fr * tg_yyyzzz_zzz_1[j] + 1.5 * fl1_fx * (tg_yyzzz_zzz_0[j] - tg_yyzzz_zzz_1[j] * fl1_fza);

                    tg_yyyzzzz_xxx_0[j] = pb_y * tg_yyzzzz_xxx_0[j] + fr * tg_yyzzzz_xxx_1[j] + fl1_fx * (tg_yzzzz_xxx_0[j] - tg_yzzzz_xxx_1[j] * fl1_fza);

                    tg_yyyzzzz_xxy_0[j] = pb_y * tg_yyzzzz_xxy_0[j] + fr * tg_yyzzzz_xxy_1[j] + fl1_fx * (tg_yzzzz_xxy_0[j] - tg_yzzzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xx_1[j];

                    tg_yyyzzzz_xxz_0[j] = pb_y * tg_yyzzzz_xxz_0[j] + fr * tg_yyzzzz_xxz_1[j] + fl1_fx * (tg_yzzzz_xxz_0[j] - tg_yzzzz_xxz_1[j] * fl1_fza);

                    tg_yyyzzzz_xyy_0[j] = pb_y * tg_yyzzzz_xyy_0[j] + fr * tg_yyzzzz_xyy_1[j] + fl1_fx * (tg_yzzzz_xyy_0[j] - tg_yzzzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_xy_1[j];

                    tg_yyyzzzz_xyz_0[j] = pb_y * tg_yyzzzz_xyz_0[j] + fr * tg_yyzzzz_xyz_1[j] + fl1_fx * (tg_yzzzz_xyz_0[j] - tg_yzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_xz_1[j];

                    tg_yyyzzzz_xzz_0[j] = pb_y * tg_yyzzzz_xzz_0[j] + fr * tg_yyzzzz_xzz_1[j] + fl1_fx * (tg_yzzzz_xzz_0[j] - tg_yzzzz_xzz_1[j] * fl1_fza);

                    tg_yyyzzzz_yyy_0[j] = pb_y * tg_yyzzzz_yyy_0[j] + fr * tg_yyzzzz_yyy_1[j] + fl1_fx * (tg_yzzzz_yyy_0[j] - tg_yzzzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzz_yy_1[j];

                    tg_yyyzzzz_yyz_0[j] = pb_y * tg_yyzzzz_yyz_0[j] + fr * tg_yyzzzz_yyz_1[j] + fl1_fx * (tg_yzzzz_yyz_0[j] - tg_yzzzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzz_yz_1[j];

                    tg_yyyzzzz_yzz_0[j] = pb_y * tg_yyzzzz_yzz_0[j] + fr * tg_yyzzzz_yzz_1[j] + fl1_fx * (tg_yzzzz_yzz_0[j] - tg_yzzzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzz_zz_1[j];

                    tg_yyyzzzz_zzz_0[j] = pb_y * tg_yyzzzz_zzz_0[j] + fr * tg_yyzzzz_zzz_1[j] + fl1_fx * (tg_yzzzz_zzz_0[j] - tg_yzzzz_zzz_1[j] * fl1_fza);

                    tg_yyzzzzz_xxx_0[j] = pb_y * tg_yzzzzz_xxx_0[j] + fr * tg_yzzzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxx_0[j] - tg_zzzzz_xxx_1[j] * fl1_fza);

                    tg_yyzzzzz_xxy_0[j] = pb_y * tg_yzzzzz_xxy_0[j] + fr * tg_yzzzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxy_0[j] - tg_zzzzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xx_1[j];

                    tg_yyzzzzz_xxz_0[j] = pb_y * tg_yzzzzz_xxz_0[j] + fr * tg_yzzzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xxz_0[j] - tg_zzzzz_xxz_1[j] * fl1_fza);

                    tg_yyzzzzz_xyy_0[j] = pb_y * tg_yzzzzz_xyy_0[j] + fr * tg_yzzzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyy_0[j] - tg_zzzzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_xy_1[j];

                    tg_yyzzzzz_xyz_0[j] = pb_y * tg_yzzzzz_xyz_0[j] + fr * tg_yzzzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xyz_0[j] - tg_zzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_xz_1[j];

                    tg_yyzzzzz_xzz_0[j] = pb_y * tg_yzzzzz_xzz_0[j] + fr * tg_yzzzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_xzz_0[j] - tg_zzzzz_xzz_1[j] * fl1_fza);

                    tg_yyzzzzz_yyy_0[j] = pb_y * tg_yzzzzz_yyy_0[j] + fr * tg_yzzzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyy_0[j] - tg_zzzzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzz_yy_1[j];

                    tg_yyzzzzz_yyz_0[j] = pb_y * tg_yzzzzz_yyz_0[j] + fr * tg_yzzzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yyz_0[j] - tg_zzzzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzz_yz_1[j];

                    tg_yyzzzzz_yzz_0[j] = pb_y * tg_yzzzzz_yzz_0[j] + fr * tg_yzzzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_yzz_0[j] - tg_zzzzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzz_zz_1[j];

                    tg_yyzzzzz_zzz_0[j] = pb_y * tg_yzzzzz_zzz_0[j] + fr * tg_yzzzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_zzzzz_zzz_0[j] - tg_zzzzz_zzz_1[j] * fl1_fza);

                    tg_yzzzzzz_xxx_0[j] = pb_y * tg_zzzzzz_xxx_0[j] + fr * tg_zzzzzz_xxx_1[j];

                    tg_yzzzzzz_xxy_0[j] = pb_y * tg_zzzzzz_xxy_0[j] + fr * tg_zzzzzz_xxy_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xx_1[j];

                    tg_yzzzzzz_xxz_0[j] = pb_y * tg_zzzzzz_xxz_0[j] + fr * tg_zzzzzz_xxz_1[j];

                    tg_yzzzzzz_xyy_0[j] = pb_y * tg_zzzzzz_xyy_0[j] + fr * tg_zzzzzz_xyy_1[j] + fl1_fxn * tg_zzzzzz_xy_1[j];

                    tg_yzzzzzz_xyz_0[j] = pb_y * tg_zzzzzz_xyz_0[j] + fr * tg_zzzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_xz_1[j];

                    tg_yzzzzzz_xzz_0[j] = pb_y * tg_zzzzzz_xzz_0[j] + fr * tg_zzzzzz_xzz_1[j];

                    tg_yzzzzzz_yyy_0[j] = pb_y * tg_zzzzzz_yyy_0[j] + fr * tg_zzzzzz_yyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzz_yy_1[j];

                    tg_yzzzzzz_yyz_0[j] = pb_y * tg_zzzzzz_yyz_0[j] + fr * tg_zzzzzz_yyz_1[j] + fl1_fxn * tg_zzzzzz_yz_1[j];

                    tg_yzzzzzz_yzz_0[j] = pb_y * tg_zzzzzz_yzz_0[j] + fr * tg_zzzzzz_yzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_zz_1[j];

                    tg_yzzzzzz_zzz_0[j] = pb_y * tg_zzzzzz_zzz_0[j] + fr * tg_zzzzzz_zzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzzz_xxx_0[j] = pb_z * tg_zzzzzz_xxx_0[j] + fr * tg_zzzzzz_xxx_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxx_0[j] - tg_zzzzz_xxx_1[j] * fl1_fza);

                    tg_zzzzzzz_xxy_0[j] = pb_z * tg_zzzzzz_xxy_0[j] + fr * tg_zzzzzz_xxy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxy_0[j] - tg_zzzzz_xxy_1[j] * fl1_fza);

                    tg_zzzzzzz_xxz_0[j] = pb_z * tg_zzzzzz_xxz_0[j] + fr * tg_zzzzzz_xxz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xxz_0[j] - tg_zzzzz_xxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xx_1[j];

                    tg_zzzzzzz_xyy_0[j] = pb_z * tg_zzzzzz_xyy_0[j] + fr * tg_zzzzzz_xyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyy_0[j] - tg_zzzzz_xyy_1[j] * fl1_fza);

                    tg_zzzzzzz_xyz_0[j] = pb_z * tg_zzzzzz_xyz_0[j] + fr * tg_zzzzzz_xyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xyz_0[j] - tg_zzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_xy_1[j];

                    tg_zzzzzzz_xzz_0[j] = pb_z * tg_zzzzzz_xzz_0[j] + fr * tg_zzzzzz_xzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_xzz_0[j] - tg_zzzzz_xzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_xz_1[j];

                    tg_zzzzzzz_yyy_0[j] = pb_z * tg_zzzzzz_yyy_0[j] + fr * tg_zzzzzz_yyy_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyy_0[j] - tg_zzzzz_yyy_1[j] * fl1_fza);

                    tg_zzzzzzz_yyz_0[j] = pb_z * tg_zzzzzz_yyz_0[j] + fr * tg_zzzzzz_yyz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yyz_0[j] - tg_zzzzz_yyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzz_yy_1[j];

                    tg_zzzzzzz_yzz_0[j] = pb_z * tg_zzzzzz_yzz_0[j] + fr * tg_zzzzzz_yzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_yzz_0[j] - tg_zzzzz_yzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzz_yz_1[j];

                    tg_zzzzzzz_zzz_0[j] = pb_z * tg_zzzzzz_zzz_0[j] + fr * tg_zzzzzz_zzz_1[j] + 3.0 * fl1_fx * (tg_zzzzz_zzz_0[j] - tg_zzzzz_zzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzz_zz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

