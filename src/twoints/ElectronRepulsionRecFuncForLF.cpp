//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForLF.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSLSF(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSLSF_0_90(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSLSF_90_180(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 

        erirecfunc::compElectronRepulsionForSLSF_180_270(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSF_270_360(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 

        erirecfunc::compElectronRepulsionForSLSF_360_450(primBuffer,
                                                         recursionMap,
                                                         osFactors,
                                                         wpDistances, 
                                                         braGtoPairsBlock,
                                                         ketGtoPairsBlock,
                                                         nKetPrimPairs,
                                                         iContrPair); 
    }

    void
    compElectronRepulsionForSLSF_0_90(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xxxxxxx_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx); 

                auto tg_xxxxxxx_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 1); 

                auto tg_xxxxxxx_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 2); 

                auto tg_xxxxxxx_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 3); 

                auto tg_xxxxxxx_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 4); 

                auto tg_xxxxxxx_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 5); 

                auto tg_xxxxxxx_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 6); 

                auto tg_xxxxxxx_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 7); 

                auto tg_xxxxxxx_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 8); 

                auto tg_xxxxxxx_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 9); 

                auto tg_xxxxxxy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 10); 

                auto tg_xxxxxxy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 11); 

                auto tg_xxxxxxy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 12); 

                auto tg_xxxxxxy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 13); 

                auto tg_xxxxxxy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 14); 

                auto tg_xxxxxxy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 15); 

                auto tg_xxxxxxy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 16); 

                auto tg_xxxxxxy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 17); 

                auto tg_xxxxxxy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 18); 

                auto tg_xxxxxxy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 19); 

                auto tg_xxxxxxz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 20); 

                auto tg_xxxxxxz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 21); 

                auto tg_xxxxxxz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 22); 

                auto tg_xxxxxxz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 23); 

                auto tg_xxxxxxz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 24); 

                auto tg_xxxxxxz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 25); 

                auto tg_xxxxxxz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 26); 

                auto tg_xxxxxxz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 27); 

                auto tg_xxxxxxz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 28); 

                auto tg_xxxxxxz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 29); 

                auto tg_xxxxxyy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 30); 

                auto tg_xxxxxyy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 31); 

                auto tg_xxxxxyy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 32); 

                auto tg_xxxxxyy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 33); 

                auto tg_xxxxxyy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 34); 

                auto tg_xxxxxyy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 35); 

                auto tg_xxxxxyy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 36); 

                auto tg_xxxxxyy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 37); 

                auto tg_xxxxxyy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 38); 

                auto tg_xxxxxyy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 39); 

                auto tg_xxxxxyz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 40); 

                auto tg_xxxxxyz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 41); 

                auto tg_xxxxxyz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 42); 

                auto tg_xxxxxyz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 43); 

                auto tg_xxxxxyz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 44); 

                auto tg_xxxxxyz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 45); 

                auto tg_xxxxxyz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 46); 

                auto tg_xxxxxyz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 47); 

                auto tg_xxxxxyz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 48); 

                auto tg_xxxxxyz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 49); 

                auto tg_xxxxxzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 50); 

                auto tg_xxxxxzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 51); 

                auto tg_xxxxxzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 52); 

                auto tg_xxxxxzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 53); 

                auto tg_xxxxxzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 54); 

                auto tg_xxxxxzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 55); 

                auto tg_xxxxxzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 56); 

                auto tg_xxxxxzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 57); 

                auto tg_xxxxxzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 58); 

                auto tg_xxxxxzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 59); 

                auto tg_xxxxyyy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 60); 

                auto tg_xxxxyyy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 61); 

                auto tg_xxxxyyy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 62); 

                auto tg_xxxxyyy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 63); 

                auto tg_xxxxyyy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 64); 

                auto tg_xxxxyyy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 65); 

                auto tg_xxxxyyy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 66); 

                auto tg_xxxxyyy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 67); 

                auto tg_xxxxyyy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 68); 

                auto tg_xxxxyyy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 69); 

                auto tg_xxxxyyz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 70); 

                auto tg_xxxxyyz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 71); 

                auto tg_xxxxyyz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 72); 

                auto tg_xxxxyyz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 73); 

                auto tg_xxxxyyz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 74); 

                auto tg_xxxxyyz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 75); 

                auto tg_xxxxyyz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 76); 

                auto tg_xxxxyyz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 77); 

                auto tg_xxxxyyz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 78); 

                auto tg_xxxxyyz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 79); 

                auto tg_xxxxyzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 80); 

                auto tg_xxxxyzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 81); 

                auto tg_xxxxyzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 82); 

                auto tg_xxxxyzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 83); 

                auto tg_xxxxyzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 84); 

                auto tg_xxxxyzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 85); 

                auto tg_xxxxyzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 86); 

                auto tg_xxxxyzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 87); 

                auto tg_xxxxyzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 88); 

                auto tg_xxxxyzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 89); 

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

                auto tg_xxxxxxx_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx); 

                auto tg_xxxxxxx_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 1); 

                auto tg_xxxxxxx_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 2); 

                auto tg_xxxxxxx_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 3); 

                auto tg_xxxxxxx_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 4); 

                auto tg_xxxxxxx_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 5); 

                auto tg_xxxxxxy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 6); 

                auto tg_xxxxxxy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 7); 

                auto tg_xxxxxxy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 8); 

                auto tg_xxxxxxy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 9); 

                auto tg_xxxxxxy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 10); 

                auto tg_xxxxxxy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 11); 

                auto tg_xxxxxxz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 12); 

                auto tg_xxxxxxz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 13); 

                auto tg_xxxxxxz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 14); 

                auto tg_xxxxxxz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 15); 

                auto tg_xxxxxxz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 16); 

                auto tg_xxxxxxz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 17); 

                auto tg_xxxxxyy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 18); 

                auto tg_xxxxxyy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 19); 

                auto tg_xxxxxyy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 20); 

                auto tg_xxxxxyy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 21); 

                auto tg_xxxxxyy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 22); 

                auto tg_xxxxxyy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 23); 

                auto tg_xxxxxyz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 24); 

                auto tg_xxxxxyz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 25); 

                auto tg_xxxxxyz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 26); 

                auto tg_xxxxxyz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 27); 

                auto tg_xxxxxyz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 28); 

                auto tg_xxxxxyz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 29); 

                auto tg_xxxxxzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 30); 

                auto tg_xxxxxzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 31); 

                auto tg_xxxxxzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 32); 

                auto tg_xxxxxzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 33); 

                auto tg_xxxxxzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 34); 

                auto tg_xxxxxzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 35); 

                auto tg_xxxxyyy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 36); 

                auto tg_xxxxyyy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 37); 

                auto tg_xxxxyyy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 38); 

                auto tg_xxxxyyy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 39); 

                auto tg_xxxxyyy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 40); 

                auto tg_xxxxyyy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 41); 

                auto tg_xxxxyyz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 42); 

                auto tg_xxxxyyz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 43); 

                auto tg_xxxxyyz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 44); 

                auto tg_xxxxyyz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 45); 

                auto tg_xxxxyyz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 46); 

                auto tg_xxxxyyz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 47); 

                auto tg_xxxxyzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 48); 

                auto tg_xxxxyzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 49); 

                auto tg_xxxxyzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 50); 

                auto tg_xxxxyzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 51); 

                auto tg_xxxxyzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 52); 

                auto tg_xxxxyzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 53); 

                // set up pointers to integrals

                auto tg_xxxxxxxx_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx); 

                auto tg_xxxxxxxx_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 1); 

                auto tg_xxxxxxxx_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 2); 

                auto tg_xxxxxxxx_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 3); 

                auto tg_xxxxxxxx_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 4); 

                auto tg_xxxxxxxx_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 5); 

                auto tg_xxxxxxxx_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 6); 

                auto tg_xxxxxxxx_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 7); 

                auto tg_xxxxxxxx_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 8); 

                auto tg_xxxxxxxx_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 9); 

                auto tg_xxxxxxxy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 10); 

                auto tg_xxxxxxxy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 11); 

                auto tg_xxxxxxxy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 12); 

                auto tg_xxxxxxxy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 13); 

                auto tg_xxxxxxxy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 14); 

                auto tg_xxxxxxxy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 15); 

                auto tg_xxxxxxxy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 16); 

                auto tg_xxxxxxxy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 17); 

                auto tg_xxxxxxxy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 18); 

                auto tg_xxxxxxxy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 19); 

                auto tg_xxxxxxxz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 20); 

                auto tg_xxxxxxxz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 21); 

                auto tg_xxxxxxxz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 22); 

                auto tg_xxxxxxxz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 23); 

                auto tg_xxxxxxxz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 24); 

                auto tg_xxxxxxxz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 25); 

                auto tg_xxxxxxxz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 26); 

                auto tg_xxxxxxxz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 27); 

                auto tg_xxxxxxxz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 28); 

                auto tg_xxxxxxxz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 29); 

                auto tg_xxxxxxyy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 30); 

                auto tg_xxxxxxyy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 31); 

                auto tg_xxxxxxyy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 32); 

                auto tg_xxxxxxyy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 33); 

                auto tg_xxxxxxyy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 34); 

                auto tg_xxxxxxyy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 35); 

                auto tg_xxxxxxyy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 36); 

                auto tg_xxxxxxyy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 37); 

                auto tg_xxxxxxyy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 38); 

                auto tg_xxxxxxyy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 39); 

                auto tg_xxxxxxyz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 40); 

                auto tg_xxxxxxyz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 41); 

                auto tg_xxxxxxyz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 42); 

                auto tg_xxxxxxyz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 43); 

                auto tg_xxxxxxyz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 44); 

                auto tg_xxxxxxyz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 45); 

                auto tg_xxxxxxyz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 46); 

                auto tg_xxxxxxyz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 47); 

                auto tg_xxxxxxyz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 48); 

                auto tg_xxxxxxyz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 49); 

                auto tg_xxxxxxzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 50); 

                auto tg_xxxxxxzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 51); 

                auto tg_xxxxxxzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 52); 

                auto tg_xxxxxxzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 53); 

                auto tg_xxxxxxzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 54); 

                auto tg_xxxxxxzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 55); 

                auto tg_xxxxxxzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 56); 

                auto tg_xxxxxxzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 57); 

                auto tg_xxxxxxzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 58); 

                auto tg_xxxxxxzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 59); 

                auto tg_xxxxxyyy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 60); 

                auto tg_xxxxxyyy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 61); 

                auto tg_xxxxxyyy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 62); 

                auto tg_xxxxxyyy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 63); 

                auto tg_xxxxxyyy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 64); 

                auto tg_xxxxxyyy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 65); 

                auto tg_xxxxxyyy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 66); 

                auto tg_xxxxxyyy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 67); 

                auto tg_xxxxxyyy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 68); 

                auto tg_xxxxxyyy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 69); 

                auto tg_xxxxxyyz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 70); 

                auto tg_xxxxxyyz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 71); 

                auto tg_xxxxxyyz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 72); 

                auto tg_xxxxxyyz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 73); 

                auto tg_xxxxxyyz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 74); 

                auto tg_xxxxxyyz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 75); 

                auto tg_xxxxxyyz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 76); 

                auto tg_xxxxxyyz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 77); 

                auto tg_xxxxxyyz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 78); 

                auto tg_xxxxxyyz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 79); 

                auto tg_xxxxxyzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 80); 

                auto tg_xxxxxyzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 81); 

                auto tg_xxxxxyzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 82); 

                auto tg_xxxxxyzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 83); 

                auto tg_xxxxxyzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 84); 

                auto tg_xxxxxyzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 85); 

                auto tg_xxxxxyzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 86); 

                auto tg_xxxxxyzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 87); 

                auto tg_xxxxxyzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 88); 

                auto tg_xxxxxyzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 89); 

                // Batch of Integrals (0,90)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxx_xxx_0, tg_xxxxxx_xxx_1, tg_xxxxxx_xxy_0, \
                                         tg_xxxxxx_xxy_1, tg_xxxxxx_xxz_0, tg_xxxxxx_xxz_1, tg_xxxxxx_xyy_0, tg_xxxxxx_xyy_1, \
                                         tg_xxxxxx_xyz_0, tg_xxxxxx_xyz_1, tg_xxxxxx_xzz_0, tg_xxxxxx_xzz_1, tg_xxxxxx_yyy_0, \
                                         tg_xxxxxx_yyy_1, tg_xxxxxx_yyz_0, tg_xxxxxx_yyz_1, tg_xxxxxx_yzz_0, tg_xxxxxx_yzz_1, \
                                         tg_xxxxxx_zzz_0, tg_xxxxxx_zzz_1, tg_xxxxxxx_xx_1, tg_xxxxxxx_xxx_0, \
                                         tg_xxxxxxx_xxx_1, tg_xxxxxxx_xxy_0, tg_xxxxxxx_xxy_1, tg_xxxxxxx_xxz_0, \
                                         tg_xxxxxxx_xxz_1, tg_xxxxxxx_xy_1, tg_xxxxxxx_xyy_0, tg_xxxxxxx_xyy_1, \
                                         tg_xxxxxxx_xyz_0, tg_xxxxxxx_xyz_1, tg_xxxxxxx_xz_1, tg_xxxxxxx_xzz_0, \
                                         tg_xxxxxxx_xzz_1, tg_xxxxxxx_yy_1, tg_xxxxxxx_yyy_0, tg_xxxxxxx_yyy_1, \
                                         tg_xxxxxxx_yyz_0, tg_xxxxxxx_yyz_1, tg_xxxxxxx_yz_1, tg_xxxxxxx_yzz_0, \
                                         tg_xxxxxxx_yzz_1, tg_xxxxxxx_zz_1, tg_xxxxxxx_zzz_0, tg_xxxxxxx_zzz_1, \
                                         tg_xxxxxxxx_xxx_0, tg_xxxxxxxx_xxy_0, tg_xxxxxxxx_xxz_0, tg_xxxxxxxx_xyy_0, \
                                         tg_xxxxxxxx_xyz_0, tg_xxxxxxxx_xzz_0, tg_xxxxxxxx_yyy_0, tg_xxxxxxxx_yyz_0, \
                                         tg_xxxxxxxx_yzz_0, tg_xxxxxxxx_zzz_0, tg_xxxxxxxy_xxx_0, tg_xxxxxxxy_xxy_0, \
                                         tg_xxxxxxxy_xxz_0, tg_xxxxxxxy_xyy_0, tg_xxxxxxxy_xyz_0, tg_xxxxxxxy_xzz_0, \
                                         tg_xxxxxxxy_yyy_0, tg_xxxxxxxy_yyz_0, tg_xxxxxxxy_yzz_0, tg_xxxxxxxy_zzz_0, \
                                         tg_xxxxxxxz_xxx_0, tg_xxxxxxxz_xxy_0, tg_xxxxxxxz_xxz_0, tg_xxxxxxxz_xyy_0, \
                                         tg_xxxxxxxz_xyz_0, tg_xxxxxxxz_xzz_0, tg_xxxxxxxz_yyy_0, tg_xxxxxxxz_yyz_0, \
                                         tg_xxxxxxxz_yzz_0, tg_xxxxxxxz_zzz_0, tg_xxxxxxy_xx_1, tg_xxxxxxy_xxx_0, \
                                         tg_xxxxxxy_xxx_1, tg_xxxxxxy_xxy_0, tg_xxxxxxy_xxy_1, tg_xxxxxxy_xxz_0, \
                                         tg_xxxxxxy_xxz_1, tg_xxxxxxy_xy_1, tg_xxxxxxy_xyy_0, tg_xxxxxxy_xyy_1, \
                                         tg_xxxxxxy_xyz_0, tg_xxxxxxy_xyz_1, tg_xxxxxxy_xz_1, tg_xxxxxxy_xzz_0, \
                                         tg_xxxxxxy_xzz_1, tg_xxxxxxy_yy_1, tg_xxxxxxy_yyy_0, tg_xxxxxxy_yyy_1, \
                                         tg_xxxxxxy_yyz_0, tg_xxxxxxy_yyz_1, tg_xxxxxxy_yz_1, tg_xxxxxxy_yzz_0, \
                                         tg_xxxxxxy_yzz_1, tg_xxxxxxy_zz_1, tg_xxxxxxy_zzz_0, tg_xxxxxxy_zzz_1, \
                                         tg_xxxxxxyy_xxx_0, tg_xxxxxxyy_xxy_0, tg_xxxxxxyy_xxz_0, tg_xxxxxxyy_xyy_0, \
                                         tg_xxxxxxyy_xyz_0, tg_xxxxxxyy_xzz_0, tg_xxxxxxyy_yyy_0, tg_xxxxxxyy_yyz_0, \
                                         tg_xxxxxxyy_yzz_0, tg_xxxxxxyy_zzz_0, tg_xxxxxxyz_xxx_0, tg_xxxxxxyz_xxy_0, \
                                         tg_xxxxxxyz_xxz_0, tg_xxxxxxyz_xyy_0, tg_xxxxxxyz_xyz_0, tg_xxxxxxyz_xzz_0, \
                                         tg_xxxxxxyz_yyy_0, tg_xxxxxxyz_yyz_0, tg_xxxxxxyz_yzz_0, tg_xxxxxxyz_zzz_0, \
                                         tg_xxxxxxz_xx_1, tg_xxxxxxz_xxx_0, tg_xxxxxxz_xxx_1, tg_xxxxxxz_xxy_0, \
                                         tg_xxxxxxz_xxy_1, tg_xxxxxxz_xxz_0, tg_xxxxxxz_xxz_1, tg_xxxxxxz_xy_1, \
                                         tg_xxxxxxz_xyy_0, tg_xxxxxxz_xyy_1, tg_xxxxxxz_xyz_0, tg_xxxxxxz_xyz_1, \
                                         tg_xxxxxxz_xz_1, tg_xxxxxxz_xzz_0, tg_xxxxxxz_xzz_1, tg_xxxxxxz_yy_1, \
                                         tg_xxxxxxz_yyy_0, tg_xxxxxxz_yyy_1, tg_xxxxxxz_yyz_0, tg_xxxxxxz_yyz_1, \
                                         tg_xxxxxxz_yz_1, tg_xxxxxxz_yzz_0, tg_xxxxxxz_yzz_1, tg_xxxxxxz_zz_1, \
                                         tg_xxxxxxz_zzz_0, tg_xxxxxxz_zzz_1, tg_xxxxxxzz_xxx_0, tg_xxxxxxzz_xxy_0, \
                                         tg_xxxxxxzz_xxz_0, tg_xxxxxxzz_xyy_0, tg_xxxxxxzz_xyz_0, tg_xxxxxxzz_xzz_0, \
                                         tg_xxxxxxzz_yyy_0, tg_xxxxxxzz_yyz_0, tg_xxxxxxzz_yzz_0, tg_xxxxxxzz_zzz_0, \
                                         tg_xxxxxy_xxx_0, tg_xxxxxy_xxx_1, tg_xxxxxy_xxy_0, tg_xxxxxy_xxy_1, tg_xxxxxy_xxz_0, \
                                         tg_xxxxxy_xxz_1, tg_xxxxxy_xyy_0, tg_xxxxxy_xyy_1, tg_xxxxxy_xyz_0, tg_xxxxxy_xyz_1, \
                                         tg_xxxxxy_xzz_0, tg_xxxxxy_xzz_1, tg_xxxxxy_yyy_0, tg_xxxxxy_yyy_1, tg_xxxxxy_yyz_0, \
                                         tg_xxxxxy_yyz_1, tg_xxxxxy_yzz_0, tg_xxxxxy_yzz_1, tg_xxxxxy_zzz_0, tg_xxxxxy_zzz_1, \
                                         tg_xxxxxyy_xx_1, tg_xxxxxyy_xxx_0, tg_xxxxxyy_xxx_1, tg_xxxxxyy_xxy_0, \
                                         tg_xxxxxyy_xxy_1, tg_xxxxxyy_xxz_0, tg_xxxxxyy_xxz_1, tg_xxxxxyy_xy_1, \
                                         tg_xxxxxyy_xyy_0, tg_xxxxxyy_xyy_1, tg_xxxxxyy_xyz_0, tg_xxxxxyy_xyz_1, \
                                         tg_xxxxxyy_xz_1, tg_xxxxxyy_xzz_0, tg_xxxxxyy_xzz_1, tg_xxxxxyy_yy_1, \
                                         tg_xxxxxyy_yyy_0, tg_xxxxxyy_yyy_1, tg_xxxxxyy_yyz_0, tg_xxxxxyy_yyz_1, \
                                         tg_xxxxxyy_yz_1, tg_xxxxxyy_yzz_0, tg_xxxxxyy_yzz_1, tg_xxxxxyy_zz_1, \
                                         tg_xxxxxyy_zzz_0, tg_xxxxxyy_zzz_1, tg_xxxxxyyy_xxx_0, tg_xxxxxyyy_xxy_0, \
                                         tg_xxxxxyyy_xxz_0, tg_xxxxxyyy_xyy_0, tg_xxxxxyyy_xyz_0, tg_xxxxxyyy_xzz_0, \
                                         tg_xxxxxyyy_yyy_0, tg_xxxxxyyy_yyz_0, tg_xxxxxyyy_yzz_0, tg_xxxxxyyy_zzz_0, \
                                         tg_xxxxxyyz_xxx_0, tg_xxxxxyyz_xxy_0, tg_xxxxxyyz_xxz_0, tg_xxxxxyyz_xyy_0, \
                                         tg_xxxxxyyz_xyz_0, tg_xxxxxyyz_xzz_0, tg_xxxxxyyz_yyy_0, tg_xxxxxyyz_yyz_0, \
                                         tg_xxxxxyyz_yzz_0, tg_xxxxxyyz_zzz_0, tg_xxxxxyz_xx_1, tg_xxxxxyz_xxx_0, \
                                         tg_xxxxxyz_xxx_1, tg_xxxxxyz_xxy_0, tg_xxxxxyz_xxy_1, tg_xxxxxyz_xxz_0, \
                                         tg_xxxxxyz_xxz_1, tg_xxxxxyz_xy_1, tg_xxxxxyz_xyy_0, tg_xxxxxyz_xyy_1, \
                                         tg_xxxxxyz_xyz_0, tg_xxxxxyz_xyz_1, tg_xxxxxyz_xz_1, tg_xxxxxyz_xzz_0, \
                                         tg_xxxxxyz_xzz_1, tg_xxxxxyz_yy_1, tg_xxxxxyz_yyy_0, tg_xxxxxyz_yyy_1, \
                                         tg_xxxxxyz_yyz_0, tg_xxxxxyz_yyz_1, tg_xxxxxyz_yz_1, tg_xxxxxyz_yzz_0, \
                                         tg_xxxxxyz_yzz_1, tg_xxxxxyz_zz_1, tg_xxxxxyz_zzz_0, tg_xxxxxyz_zzz_1, \
                                         tg_xxxxxyzz_xxx_0, tg_xxxxxyzz_xxy_0, tg_xxxxxyzz_xxz_0, tg_xxxxxyzz_xyy_0, \
                                         tg_xxxxxyzz_xyz_0, tg_xxxxxyzz_xzz_0, tg_xxxxxyzz_yyy_0, tg_xxxxxyzz_yyz_0, \
                                         tg_xxxxxyzz_yzz_0, tg_xxxxxyzz_zzz_0, tg_xxxxxz_xxx_0, tg_xxxxxz_xxx_1, \
                                         tg_xxxxxz_xxy_0, tg_xxxxxz_xxy_1, tg_xxxxxz_xxz_0, tg_xxxxxz_xxz_1, tg_xxxxxz_xyy_0, \
                                         tg_xxxxxz_xyy_1, tg_xxxxxz_xyz_0, tg_xxxxxz_xyz_1, tg_xxxxxz_xzz_0, tg_xxxxxz_xzz_1, \
                                         tg_xxxxxz_yyy_0, tg_xxxxxz_yyy_1, tg_xxxxxz_yyz_0, tg_xxxxxz_yyz_1, tg_xxxxxz_yzz_0, \
                                         tg_xxxxxz_yzz_1, tg_xxxxxz_zzz_0, tg_xxxxxz_zzz_1, tg_xxxxxzz_xx_1, \
                                         tg_xxxxxzz_xxx_0, tg_xxxxxzz_xxx_1, tg_xxxxxzz_xxy_0, tg_xxxxxzz_xxy_1, \
                                         tg_xxxxxzz_xxz_0, tg_xxxxxzz_xxz_1, tg_xxxxxzz_xy_1, tg_xxxxxzz_xyy_0, \
                                         tg_xxxxxzz_xyy_1, tg_xxxxxzz_xyz_0, tg_xxxxxzz_xyz_1, tg_xxxxxzz_xz_1, \
                                         tg_xxxxxzz_xzz_0, tg_xxxxxzz_xzz_1, tg_xxxxxzz_yy_1, tg_xxxxxzz_yyy_0, \
                                         tg_xxxxxzz_yyy_1, tg_xxxxxzz_yyz_0, tg_xxxxxzz_yyz_1, tg_xxxxxzz_yz_1, \
                                         tg_xxxxxzz_yzz_0, tg_xxxxxzz_yzz_1, tg_xxxxxzz_zz_1, tg_xxxxxzz_zzz_0, \
                                         tg_xxxxxzz_zzz_1, tg_xxxxyy_xxx_0, tg_xxxxyy_xxx_1, tg_xxxxyy_xxy_0, tg_xxxxyy_xxy_1, \
                                         tg_xxxxyy_xxz_0, tg_xxxxyy_xxz_1, tg_xxxxyy_xyy_0, tg_xxxxyy_xyy_1, tg_xxxxyy_xyz_0, \
                                         tg_xxxxyy_xyz_1, tg_xxxxyy_xzz_0, tg_xxxxyy_xzz_1, tg_xxxxyy_yyy_0, tg_xxxxyy_yyy_1, \
                                         tg_xxxxyy_yyz_0, tg_xxxxyy_yyz_1, tg_xxxxyy_yzz_0, tg_xxxxyy_yzz_1, tg_xxxxyy_zzz_0, \
                                         tg_xxxxyy_zzz_1, tg_xxxxyyy_xx_1, tg_xxxxyyy_xxx_0, tg_xxxxyyy_xxx_1, \
                                         tg_xxxxyyy_xxy_0, tg_xxxxyyy_xxy_1, tg_xxxxyyy_xxz_0, tg_xxxxyyy_xxz_1, \
                                         tg_xxxxyyy_xy_1, tg_xxxxyyy_xyy_0, tg_xxxxyyy_xyy_1, tg_xxxxyyy_xyz_0, \
                                         tg_xxxxyyy_xyz_1, tg_xxxxyyy_xz_1, tg_xxxxyyy_xzz_0, tg_xxxxyyy_xzz_1, \
                                         tg_xxxxyyy_yy_1, tg_xxxxyyy_yyy_0, tg_xxxxyyy_yyy_1, tg_xxxxyyy_yyz_0, \
                                         tg_xxxxyyy_yyz_1, tg_xxxxyyy_yz_1, tg_xxxxyyy_yzz_0, tg_xxxxyyy_yzz_1, \
                                         tg_xxxxyyy_zz_1, tg_xxxxyyy_zzz_0, tg_xxxxyyy_zzz_1, tg_xxxxyyz_xx_1, \
                                         tg_xxxxyyz_xxx_0, tg_xxxxyyz_xxx_1, tg_xxxxyyz_xxy_0, tg_xxxxyyz_xxy_1, \
                                         tg_xxxxyyz_xxz_0, tg_xxxxyyz_xxz_1, tg_xxxxyyz_xy_1, tg_xxxxyyz_xyy_0, \
                                         tg_xxxxyyz_xyy_1, tg_xxxxyyz_xyz_0, tg_xxxxyyz_xyz_1, tg_xxxxyyz_xz_1, \
                                         tg_xxxxyyz_xzz_0, tg_xxxxyyz_xzz_1, tg_xxxxyyz_yy_1, tg_xxxxyyz_yyy_0, \
                                         tg_xxxxyyz_yyy_1, tg_xxxxyyz_yyz_0, tg_xxxxyyz_yyz_1, tg_xxxxyyz_yz_1, \
                                         tg_xxxxyyz_yzz_0, tg_xxxxyyz_yzz_1, tg_xxxxyyz_zz_1, tg_xxxxyyz_zzz_0, \
                                         tg_xxxxyyz_zzz_1, tg_xxxxyz_xxx_0, tg_xxxxyz_xxx_1, tg_xxxxyz_xxy_0, tg_xxxxyz_xxy_1, \
                                         tg_xxxxyz_xxz_0, tg_xxxxyz_xxz_1, tg_xxxxyz_xyy_0, tg_xxxxyz_xyy_1, tg_xxxxyz_xyz_0, \
                                         tg_xxxxyz_xyz_1, tg_xxxxyz_xzz_0, tg_xxxxyz_xzz_1, tg_xxxxyz_yyy_0, tg_xxxxyz_yyy_1, \
                                         tg_xxxxyz_yyz_0, tg_xxxxyz_yyz_1, tg_xxxxyz_yzz_0, tg_xxxxyz_yzz_1, tg_xxxxyz_zzz_0, \
                                         tg_xxxxyz_zzz_1, tg_xxxxyzz_xx_1, tg_xxxxyzz_xxx_0, tg_xxxxyzz_xxx_1, \
                                         tg_xxxxyzz_xxy_0, tg_xxxxyzz_xxy_1, tg_xxxxyzz_xxz_0, tg_xxxxyzz_xxz_1, \
                                         tg_xxxxyzz_xy_1, tg_xxxxyzz_xyy_0, tg_xxxxyzz_xyy_1, tg_xxxxyzz_xyz_0, \
                                         tg_xxxxyzz_xyz_1, tg_xxxxyzz_xz_1, tg_xxxxyzz_xzz_0, tg_xxxxyzz_xzz_1, \
                                         tg_xxxxyzz_yy_1, tg_xxxxyzz_yyy_0, tg_xxxxyzz_yyy_1, tg_xxxxyzz_yyz_0, \
                                         tg_xxxxyzz_yyz_1, tg_xxxxyzz_yz_1, tg_xxxxyzz_yzz_0, tg_xxxxyzz_yzz_1, \
                                         tg_xxxxyzz_zz_1, tg_xxxxyzz_zzz_0, tg_xxxxyzz_zzz_1, tg_xxxxzz_xxx_0, \
                                         tg_xxxxzz_xxx_1, tg_xxxxzz_xxy_0, tg_xxxxzz_xxy_1, tg_xxxxzz_xxz_0, tg_xxxxzz_xxz_1, \
                                         tg_xxxxzz_xyy_0, tg_xxxxzz_xyy_1, tg_xxxxzz_xyz_0, tg_xxxxzz_xyz_1, tg_xxxxzz_xzz_0, \
                                         tg_xxxxzz_xzz_1, tg_xxxxzz_yyy_0, tg_xxxxzz_yyy_1, tg_xxxxzz_yyz_0, tg_xxxxzz_yyz_1, \
                                         tg_xxxxzz_yzz_0, tg_xxxxzz_yzz_1, tg_xxxxzz_zzz_0, tg_xxxxzz_zzz_1, tg_xxxyyy_xxx_0, \
                                         tg_xxxyyy_xxx_1, tg_xxxyyy_xxy_0, tg_xxxyyy_xxy_1, tg_xxxyyy_xxz_0, tg_xxxyyy_xxz_1, \
                                         tg_xxxyyy_xyy_0, tg_xxxyyy_xyy_1, tg_xxxyyy_xyz_0, tg_xxxyyy_xyz_1, tg_xxxyyy_xzz_0, \
                                         tg_xxxyyy_xzz_1, tg_xxxyyy_yyy_0, tg_xxxyyy_yyy_1, tg_xxxyyy_yyz_0, tg_xxxyyy_yyz_1, \
                                         tg_xxxyyy_yzz_0, tg_xxxyyy_yzz_1, tg_xxxyyy_zzz_0, tg_xxxyyy_zzz_1, tg_xxxyyz_xxx_0, \
                                         tg_xxxyyz_xxx_1, tg_xxxyyz_xxy_0, tg_xxxyyz_xxy_1, tg_xxxyyz_xxz_0, tg_xxxyyz_xxz_1, \
                                         tg_xxxyyz_xyy_0, tg_xxxyyz_xyy_1, tg_xxxyyz_xyz_0, tg_xxxyyz_xyz_1, tg_xxxyyz_xzz_0, \
                                         tg_xxxyyz_xzz_1, tg_xxxyyz_yyy_0, tg_xxxyyz_yyy_1, tg_xxxyyz_yyz_0, tg_xxxyyz_yyz_1, \
                                         tg_xxxyyz_yzz_0, tg_xxxyyz_yzz_1, tg_xxxyyz_zzz_0, tg_xxxyyz_zzz_1, tg_xxxyzz_xxx_0, \
                                         tg_xxxyzz_xxx_1, tg_xxxyzz_xxy_0, tg_xxxyzz_xxy_1, tg_xxxyzz_xxz_0, tg_xxxyzz_xxz_1, \
                                         tg_xxxyzz_xyy_0, tg_xxxyzz_xyy_1, tg_xxxyzz_xyz_0, tg_xxxyzz_xyz_1, tg_xxxyzz_xzz_0, \
                                         tg_xxxyzz_xzz_1, tg_xxxyzz_yyy_0, tg_xxxyzz_yyy_1, tg_xxxyzz_yyz_0, tg_xxxyzz_yyz_1, \
                                         tg_xxxyzz_yzz_0, tg_xxxyzz_yzz_1, tg_xxxyzz_zzz_0, tg_xxxyzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxxxx_xxx_0[j] = pb_x * tg_xxxxxxx_xxx_0[j] + fr * tg_xxxxxxx_xxx_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxx_0[j] - tg_xxxxxx_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxx_xx_1[j];

                    tg_xxxxxxxx_xxy_0[j] = pb_x * tg_xxxxxxx_xxy_0[j] + fr * tg_xxxxxxx_xxy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxy_0[j] - tg_xxxxxx_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxx_xy_1[j];

                    tg_xxxxxxxx_xxz_0[j] = pb_x * tg_xxxxxxx_xxz_0[j] + fr * tg_xxxxxxx_xxz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xxz_0[j] - tg_xxxxxx_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxx_xz_1[j];

                    tg_xxxxxxxx_xyy_0[j] = pb_x * tg_xxxxxxx_xyy_0[j] + fr * tg_xxxxxxx_xyy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xyy_0[j] - tg_xxxxxx_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_yy_1[j];

                    tg_xxxxxxxx_xyz_0[j] = pb_x * tg_xxxxxxx_xyz_0[j] + fr * tg_xxxxxxx_xyz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xyz_0[j] - tg_xxxxxx_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_yz_1[j];

                    tg_xxxxxxxx_xzz_0[j] = pb_x * tg_xxxxxxx_xzz_0[j] + fr * tg_xxxxxxx_xzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_xzz_0[j] - tg_xxxxxx_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxx_zz_1[j];

                    tg_xxxxxxxx_yyy_0[j] = pb_x * tg_xxxxxxx_yyy_0[j] + fr * tg_xxxxxxx_yyy_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yyy_0[j] - tg_xxxxxx_yyy_1[j] * fl1_fza);

                    tg_xxxxxxxx_yyz_0[j] = pb_x * tg_xxxxxxx_yyz_0[j] + fr * tg_xxxxxxx_yyz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yyz_0[j] - tg_xxxxxx_yyz_1[j] * fl1_fza);

                    tg_xxxxxxxx_yzz_0[j] = pb_x * tg_xxxxxxx_yzz_0[j] + fr * tg_xxxxxxx_yzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_yzz_0[j] - tg_xxxxxx_yzz_1[j] * fl1_fza);

                    tg_xxxxxxxx_zzz_0[j] = pb_x * tg_xxxxxxx_zzz_0[j] + fr * tg_xxxxxxx_zzz_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_zzz_0[j] - tg_xxxxxx_zzz_1[j] * fl1_fza);

                    tg_xxxxxxxy_xxx_0[j] = pb_x * tg_xxxxxxy_xxx_0[j] + fr * tg_xxxxxxy_xxx_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxx_0[j] - tg_xxxxxy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxy_xx_1[j];

                    tg_xxxxxxxy_xxy_0[j] = pb_x * tg_xxxxxxy_xxy_0[j] + fr * tg_xxxxxxy_xxy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxy_0[j] - tg_xxxxxy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxy_xy_1[j];

                    tg_xxxxxxxy_xxz_0[j] = pb_x * tg_xxxxxxy_xxz_0[j] + fr * tg_xxxxxxy_xxz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xxz_0[j] - tg_xxxxxy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxy_xz_1[j];

                    tg_xxxxxxxy_xyy_0[j] = pb_x * tg_xxxxxxy_xyy_0[j] + fr * tg_xxxxxxy_xyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xyy_0[j] - tg_xxxxxy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_yy_1[j];

                    tg_xxxxxxxy_xyz_0[j] = pb_x * tg_xxxxxxy_xyz_0[j] + fr * tg_xxxxxxy_xyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xyz_0[j] - tg_xxxxxy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_yz_1[j];

                    tg_xxxxxxxy_xzz_0[j] = pb_x * tg_xxxxxxy_xzz_0[j] + fr * tg_xxxxxxy_xzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_xzz_0[j] - tg_xxxxxy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxy_zz_1[j];

                    tg_xxxxxxxy_yyy_0[j] = pb_x * tg_xxxxxxy_yyy_0[j] + fr * tg_xxxxxxy_yyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yyy_0[j] - tg_xxxxxy_yyy_1[j] * fl1_fza);

                    tg_xxxxxxxy_yyz_0[j] = pb_x * tg_xxxxxxy_yyz_0[j] + fr * tg_xxxxxxy_yyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yyz_0[j] - tg_xxxxxy_yyz_1[j] * fl1_fza);

                    tg_xxxxxxxy_yzz_0[j] = pb_x * tg_xxxxxxy_yzz_0[j] + fr * tg_xxxxxxy_yzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_yzz_0[j] - tg_xxxxxy_yzz_1[j] * fl1_fza);

                    tg_xxxxxxxy_zzz_0[j] = pb_x * tg_xxxxxxy_zzz_0[j] + fr * tg_xxxxxxy_zzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_zzz_0[j] - tg_xxxxxy_zzz_1[j] * fl1_fza);

                    tg_xxxxxxxz_xxx_0[j] = pb_x * tg_xxxxxxz_xxx_0[j] + fr * tg_xxxxxxz_xxx_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxx_0[j] - tg_xxxxxz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxxz_xx_1[j];

                    tg_xxxxxxxz_xxy_0[j] = pb_x * tg_xxxxxxz_xxy_0[j] + fr * tg_xxxxxxz_xxy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxy_0[j] - tg_xxxxxz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxz_xy_1[j];

                    tg_xxxxxxxz_xxz_0[j] = pb_x * tg_xxxxxxz_xxz_0[j] + fr * tg_xxxxxxz_xxz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xxz_0[j] - tg_xxxxxz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxxz_xz_1[j];

                    tg_xxxxxxxz_xyy_0[j] = pb_x * tg_xxxxxxz_xyy_0[j] + fr * tg_xxxxxxz_xyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xyy_0[j] - tg_xxxxxz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_yy_1[j];

                    tg_xxxxxxxz_xyz_0[j] = pb_x * tg_xxxxxxz_xyz_0[j] + fr * tg_xxxxxxz_xyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xyz_0[j] - tg_xxxxxz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_yz_1[j];

                    tg_xxxxxxxz_xzz_0[j] = pb_x * tg_xxxxxxz_xzz_0[j] + fr * tg_xxxxxxz_xzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_xzz_0[j] - tg_xxxxxz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxxz_zz_1[j];

                    tg_xxxxxxxz_yyy_0[j] = pb_x * tg_xxxxxxz_yyy_0[j] + fr * tg_xxxxxxz_yyy_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yyy_0[j] - tg_xxxxxz_yyy_1[j] * fl1_fza);

                    tg_xxxxxxxz_yyz_0[j] = pb_x * tg_xxxxxxz_yyz_0[j] + fr * tg_xxxxxxz_yyz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yyz_0[j] - tg_xxxxxz_yyz_1[j] * fl1_fza);

                    tg_xxxxxxxz_yzz_0[j] = pb_x * tg_xxxxxxz_yzz_0[j] + fr * tg_xxxxxxz_yzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_yzz_0[j] - tg_xxxxxz_yzz_1[j] * fl1_fza);

                    tg_xxxxxxxz_zzz_0[j] = pb_x * tg_xxxxxxz_zzz_0[j] + fr * tg_xxxxxxz_zzz_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_zzz_0[j] - tg_xxxxxz_zzz_1[j] * fl1_fza);

                    tg_xxxxxxyy_xxx_0[j] = pb_x * tg_xxxxxyy_xxx_0[j] + fr * tg_xxxxxyy_xxx_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxx_0[j] - tg_xxxxyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyy_xx_1[j];

                    tg_xxxxxxyy_xxy_0[j] = pb_x * tg_xxxxxyy_xxy_0[j] + fr * tg_xxxxxyy_xxy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxy_0[j] - tg_xxxxyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyy_xy_1[j];

                    tg_xxxxxxyy_xxz_0[j] = pb_x * tg_xxxxxyy_xxz_0[j] + fr * tg_xxxxxyy_xxz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xxz_0[j] - tg_xxxxyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyy_xz_1[j];

                    tg_xxxxxxyy_xyy_0[j] = pb_x * tg_xxxxxyy_xyy_0[j] + fr * tg_xxxxxyy_xyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xyy_0[j] - tg_xxxxyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_yy_1[j];

                    tg_xxxxxxyy_xyz_0[j] = pb_x * tg_xxxxxyy_xyz_0[j] + fr * tg_xxxxxyy_xyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xyz_0[j] - tg_xxxxyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_yz_1[j];

                    tg_xxxxxxyy_xzz_0[j] = pb_x * tg_xxxxxyy_xzz_0[j] + fr * tg_xxxxxyy_xzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_xzz_0[j] - tg_xxxxyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyy_zz_1[j];

                    tg_xxxxxxyy_yyy_0[j] = pb_x * tg_xxxxxyy_yyy_0[j] + fr * tg_xxxxxyy_yyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yyy_0[j] - tg_xxxxyy_yyy_1[j] * fl1_fza);

                    tg_xxxxxxyy_yyz_0[j] = pb_x * tg_xxxxxyy_yyz_0[j] + fr * tg_xxxxxyy_yyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yyz_0[j] - tg_xxxxyy_yyz_1[j] * fl1_fza);

                    tg_xxxxxxyy_yzz_0[j] = pb_x * tg_xxxxxyy_yzz_0[j] + fr * tg_xxxxxyy_yzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_yzz_0[j] - tg_xxxxyy_yzz_1[j] * fl1_fza);

                    tg_xxxxxxyy_zzz_0[j] = pb_x * tg_xxxxxyy_zzz_0[j] + fr * tg_xxxxxyy_zzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_zzz_0[j] - tg_xxxxyy_zzz_1[j] * fl1_fza);

                    tg_xxxxxxyz_xxx_0[j] = pb_x * tg_xxxxxyz_xxx_0[j] + fr * tg_xxxxxyz_xxx_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxx_0[j] - tg_xxxxyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxyz_xx_1[j];

                    tg_xxxxxxyz_xxy_0[j] = pb_x * tg_xxxxxyz_xxy_0[j] + fr * tg_xxxxxyz_xxy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxy_0[j] - tg_xxxxyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyz_xy_1[j];

                    tg_xxxxxxyz_xxz_0[j] = pb_x * tg_xxxxxyz_xxz_0[j] + fr * tg_xxxxxyz_xxz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xxz_0[j] - tg_xxxxyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxyz_xz_1[j];

                    tg_xxxxxxyz_xyy_0[j] = pb_x * tg_xxxxxyz_xyy_0[j] + fr * tg_xxxxxyz_xyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xyy_0[j] - tg_xxxxyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_yy_1[j];

                    tg_xxxxxxyz_xyz_0[j] = pb_x * tg_xxxxxyz_xyz_0[j] + fr * tg_xxxxxyz_xyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xyz_0[j] - tg_xxxxyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_yz_1[j];

                    tg_xxxxxxyz_xzz_0[j] = pb_x * tg_xxxxxyz_xzz_0[j] + fr * tg_xxxxxyz_xzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_xzz_0[j] - tg_xxxxyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxyz_zz_1[j];

                    tg_xxxxxxyz_yyy_0[j] = pb_x * tg_xxxxxyz_yyy_0[j] + fr * tg_xxxxxyz_yyy_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yyy_0[j] - tg_xxxxyz_yyy_1[j] * fl1_fza);

                    tg_xxxxxxyz_yyz_0[j] = pb_x * tg_xxxxxyz_yyz_0[j] + fr * tg_xxxxxyz_yyz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yyz_0[j] - tg_xxxxyz_yyz_1[j] * fl1_fza);

                    tg_xxxxxxyz_yzz_0[j] = pb_x * tg_xxxxxyz_yzz_0[j] + fr * tg_xxxxxyz_yzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_yzz_0[j] - tg_xxxxyz_yzz_1[j] * fl1_fza);

                    tg_xxxxxxyz_zzz_0[j] = pb_x * tg_xxxxxyz_zzz_0[j] + fr * tg_xxxxxyz_zzz_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_zzz_0[j] - tg_xxxxyz_zzz_1[j] * fl1_fza);

                    tg_xxxxxxzz_xxx_0[j] = pb_x * tg_xxxxxzz_xxx_0[j] + fr * tg_xxxxxzz_xxx_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxx_0[j] - tg_xxxxzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxxzz_xx_1[j];

                    tg_xxxxxxzz_xxy_0[j] = pb_x * tg_xxxxxzz_xxy_0[j] + fr * tg_xxxxxzz_xxy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxy_0[j] - tg_xxxxzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxzz_xy_1[j];

                    tg_xxxxxxzz_xxz_0[j] = pb_x * tg_xxxxxzz_xxz_0[j] + fr * tg_xxxxxzz_xxz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xxz_0[j] - tg_xxxxzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxxzz_xz_1[j];

                    tg_xxxxxxzz_xyy_0[j] = pb_x * tg_xxxxxzz_xyy_0[j] + fr * tg_xxxxxzz_xyy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xyy_0[j] - tg_xxxxzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_yy_1[j];

                    tg_xxxxxxzz_xyz_0[j] = pb_x * tg_xxxxxzz_xyz_0[j] + fr * tg_xxxxxzz_xyz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xyz_0[j] - tg_xxxxzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_yz_1[j];

                    tg_xxxxxxzz_xzz_0[j] = pb_x * tg_xxxxxzz_xzz_0[j] + fr * tg_xxxxxzz_xzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_xzz_0[j] - tg_xxxxzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxxzz_zz_1[j];

                    tg_xxxxxxzz_yyy_0[j] = pb_x * tg_xxxxxzz_yyy_0[j] + fr * tg_xxxxxzz_yyy_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yyy_0[j] - tg_xxxxzz_yyy_1[j] * fl1_fza);

                    tg_xxxxxxzz_yyz_0[j] = pb_x * tg_xxxxxzz_yyz_0[j] + fr * tg_xxxxxzz_yyz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yyz_0[j] - tg_xxxxzz_yyz_1[j] * fl1_fza);

                    tg_xxxxxxzz_yzz_0[j] = pb_x * tg_xxxxxzz_yzz_0[j] + fr * tg_xxxxxzz_yzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_yzz_0[j] - tg_xxxxzz_yzz_1[j] * fl1_fza);

                    tg_xxxxxxzz_zzz_0[j] = pb_x * tg_xxxxxzz_zzz_0[j] + fr * tg_xxxxxzz_zzz_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_zzz_0[j] - tg_xxxxzz_zzz_1[j] * fl1_fza);

                    tg_xxxxxyyy_xxx_0[j] = pb_x * tg_xxxxyyy_xxx_0[j] + fr * tg_xxxxyyy_xxx_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxx_0[j] - tg_xxxyyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyy_xx_1[j];

                    tg_xxxxxyyy_xxy_0[j] = pb_x * tg_xxxxyyy_xxy_0[j] + fr * tg_xxxxyyy_xxy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxy_0[j] - tg_xxxyyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyy_xy_1[j];

                    tg_xxxxxyyy_xxz_0[j] = pb_x * tg_xxxxyyy_xxz_0[j] + fr * tg_xxxxyyy_xxz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xxz_0[j] - tg_xxxyyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyy_xz_1[j];

                    tg_xxxxxyyy_xyy_0[j] = pb_x * tg_xxxxyyy_xyy_0[j] + fr * tg_xxxxyyy_xyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xyy_0[j] - tg_xxxyyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_yy_1[j];

                    tg_xxxxxyyy_xyz_0[j] = pb_x * tg_xxxxyyy_xyz_0[j] + fr * tg_xxxxyyy_xyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xyz_0[j] - tg_xxxyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_yz_1[j];

                    tg_xxxxxyyy_xzz_0[j] = pb_x * tg_xxxxyyy_xzz_0[j] + fr * tg_xxxxyyy_xzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_xzz_0[j] - tg_xxxyyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyy_zz_1[j];

                    tg_xxxxxyyy_yyy_0[j] = pb_x * tg_xxxxyyy_yyy_0[j] + fr * tg_xxxxyyy_yyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yyy_0[j] - tg_xxxyyy_yyy_1[j] * fl1_fza);

                    tg_xxxxxyyy_yyz_0[j] = pb_x * tg_xxxxyyy_yyz_0[j] + fr * tg_xxxxyyy_yyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yyz_0[j] - tg_xxxyyy_yyz_1[j] * fl1_fza);

                    tg_xxxxxyyy_yzz_0[j] = pb_x * tg_xxxxyyy_yzz_0[j] + fr * tg_xxxxyyy_yzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_yzz_0[j] - tg_xxxyyy_yzz_1[j] * fl1_fza);

                    tg_xxxxxyyy_zzz_0[j] = pb_x * tg_xxxxyyy_zzz_0[j] + fr * tg_xxxxyyy_zzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_zzz_0[j] - tg_xxxyyy_zzz_1[j] * fl1_fza);

                    tg_xxxxxyyz_xxx_0[j] = pb_x * tg_xxxxyyz_xxx_0[j] + fr * tg_xxxxyyz_xxx_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxx_0[j] - tg_xxxyyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyyz_xx_1[j];

                    tg_xxxxxyyz_xxy_0[j] = pb_x * tg_xxxxyyz_xxy_0[j] + fr * tg_xxxxyyz_xxy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxy_0[j] - tg_xxxyyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyz_xy_1[j];

                    tg_xxxxxyyz_xxz_0[j] = pb_x * tg_xxxxyyz_xxz_0[j] + fr * tg_xxxxyyz_xxz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xxz_0[j] - tg_xxxyyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyyz_xz_1[j];

                    tg_xxxxxyyz_xyy_0[j] = pb_x * tg_xxxxyyz_xyy_0[j] + fr * tg_xxxxyyz_xyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xyy_0[j] - tg_xxxyyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_yy_1[j];

                    tg_xxxxxyyz_xyz_0[j] = pb_x * tg_xxxxyyz_xyz_0[j] + fr * tg_xxxxyyz_xyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xyz_0[j] - tg_xxxyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_yz_1[j];

                    tg_xxxxxyyz_xzz_0[j] = pb_x * tg_xxxxyyz_xzz_0[j] + fr * tg_xxxxyyz_xzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_xzz_0[j] - tg_xxxyyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyyz_zz_1[j];

                    tg_xxxxxyyz_yyy_0[j] = pb_x * tg_xxxxyyz_yyy_0[j] + fr * tg_xxxxyyz_yyy_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yyy_0[j] - tg_xxxyyz_yyy_1[j] * fl1_fza);

                    tg_xxxxxyyz_yyz_0[j] = pb_x * tg_xxxxyyz_yyz_0[j] + fr * tg_xxxxyyz_yyz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yyz_0[j] - tg_xxxyyz_yyz_1[j] * fl1_fza);

                    tg_xxxxxyyz_yzz_0[j] = pb_x * tg_xxxxyyz_yzz_0[j] + fr * tg_xxxxyyz_yzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_yzz_0[j] - tg_xxxyyz_yzz_1[j] * fl1_fza);

                    tg_xxxxxyyz_zzz_0[j] = pb_x * tg_xxxxyyz_zzz_0[j] + fr * tg_xxxxyyz_zzz_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_zzz_0[j] - tg_xxxyyz_zzz_1[j] * fl1_fza);

                    tg_xxxxxyzz_xxx_0[j] = pb_x * tg_xxxxyzz_xxx_0[j] + fr * tg_xxxxyzz_xxx_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxx_0[j] - tg_xxxyzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxyzz_xx_1[j];

                    tg_xxxxxyzz_xxy_0[j] = pb_x * tg_xxxxyzz_xxy_0[j] + fr * tg_xxxxyzz_xxy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxy_0[j] - tg_xxxyzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyzz_xy_1[j];

                    tg_xxxxxyzz_xxz_0[j] = pb_x * tg_xxxxyzz_xxz_0[j] + fr * tg_xxxxyzz_xxz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xxz_0[j] - tg_xxxyzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxyzz_xz_1[j];

                    tg_xxxxxyzz_xyy_0[j] = pb_x * tg_xxxxyzz_xyy_0[j] + fr * tg_xxxxyzz_xyy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xyy_0[j] - tg_xxxyzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_yy_1[j];

                    tg_xxxxxyzz_xyz_0[j] = pb_x * tg_xxxxyzz_xyz_0[j] + fr * tg_xxxxyzz_xyz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xyz_0[j] - tg_xxxyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_yz_1[j];

                    tg_xxxxxyzz_xzz_0[j] = pb_x * tg_xxxxyzz_xzz_0[j] + fr * tg_xxxxyzz_xzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_xzz_0[j] - tg_xxxyzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxyzz_zz_1[j];

                    tg_xxxxxyzz_yyy_0[j] = pb_x * tg_xxxxyzz_yyy_0[j] + fr * tg_xxxxyzz_yyy_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yyy_0[j] - tg_xxxyzz_yyy_1[j] * fl1_fza);

                    tg_xxxxxyzz_yyz_0[j] = pb_x * tg_xxxxyzz_yyz_0[j] + fr * tg_xxxxyzz_yyz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yyz_0[j] - tg_xxxyzz_yyz_1[j] * fl1_fza);

                    tg_xxxxxyzz_yzz_0[j] = pb_x * tg_xxxxyzz_yzz_0[j] + fr * tg_xxxxyzz_yzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_yzz_0[j] - tg_xxxyzz_yzz_1[j] * fl1_fza);

                    tg_xxxxxyzz_zzz_0[j] = pb_x * tg_xxxxyzz_zzz_0[j] + fr * tg_xxxxyzz_zzz_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_zzz_0[j] - tg_xxxyzz_zzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSF_90_180(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xxxxzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 90); 

                auto tg_xxxxzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 91); 

                auto tg_xxxxzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 92); 

                auto tg_xxxxzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 93); 

                auto tg_xxxxzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 94); 

                auto tg_xxxxzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 95); 

                auto tg_xxxxzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 96); 

                auto tg_xxxxzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 97); 

                auto tg_xxxxzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 98); 

                auto tg_xxxxzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 99); 

                auto tg_xxxyyyy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 100); 

                auto tg_xxxyyyy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 101); 

                auto tg_xxxyyyy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 102); 

                auto tg_xxxyyyy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 103); 

                auto tg_xxxyyyy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 104); 

                auto tg_xxxyyyy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 105); 

                auto tg_xxxyyyy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 106); 

                auto tg_xxxyyyy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 107); 

                auto tg_xxxyyyy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 108); 

                auto tg_xxxyyyy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 109); 

                auto tg_xxxyyyz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 110); 

                auto tg_xxxyyyz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 111); 

                auto tg_xxxyyyz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 112); 

                auto tg_xxxyyyz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 113); 

                auto tg_xxxyyyz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 114); 

                auto tg_xxxyyyz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 115); 

                auto tg_xxxyyyz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 116); 

                auto tg_xxxyyyz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 117); 

                auto tg_xxxyyyz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 118); 

                auto tg_xxxyyyz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 119); 

                auto tg_xxxyyzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 120); 

                auto tg_xxxyyzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 121); 

                auto tg_xxxyyzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 122); 

                auto tg_xxxyyzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 123); 

                auto tg_xxxyyzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 124); 

                auto tg_xxxyyzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 125); 

                auto tg_xxxyyzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 126); 

                auto tg_xxxyyzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 127); 

                auto tg_xxxyyzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 128); 

                auto tg_xxxyyzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 129); 

                auto tg_xxxyzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 130); 

                auto tg_xxxyzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 131); 

                auto tg_xxxyzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 132); 

                auto tg_xxxyzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 133); 

                auto tg_xxxyzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 134); 

                auto tg_xxxyzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 135); 

                auto tg_xxxyzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 136); 

                auto tg_xxxyzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 137); 

                auto tg_xxxyzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 138); 

                auto tg_xxxyzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 139); 

                auto tg_xxxzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 140); 

                auto tg_xxxzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 141); 

                auto tg_xxxzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 142); 

                auto tg_xxxzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 143); 

                auto tg_xxxzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 144); 

                auto tg_xxxzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 145); 

                auto tg_xxxzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 146); 

                auto tg_xxxzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 147); 

                auto tg_xxxzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 148); 

                auto tg_xxxzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 149); 

                auto tg_xxyyyyy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 150); 

                auto tg_xxyyyyy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 151); 

                auto tg_xxyyyyy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 152); 

                auto tg_xxyyyyy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 153); 

                auto tg_xxyyyyy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 154); 

                auto tg_xxyyyyy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 155); 

                auto tg_xxyyyyy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 156); 

                auto tg_xxyyyyy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 157); 

                auto tg_xxyyyyy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 158); 

                auto tg_xxyyyyy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 159); 

                auto tg_xxyyyyz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 160); 

                auto tg_xxyyyyz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 161); 

                auto tg_xxyyyyz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 162); 

                auto tg_xxyyyyz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 163); 

                auto tg_xxyyyyz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 164); 

                auto tg_xxyyyyz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 165); 

                auto tg_xxyyyyz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 166); 

                auto tg_xxyyyyz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 167); 

                auto tg_xxyyyyz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 168); 

                auto tg_xxyyyyz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 169); 

                auto tg_xxyyyzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 170); 

                auto tg_xxyyyzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 171); 

                auto tg_xxyyyzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 172); 

                auto tg_xxyyyzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 173); 

                auto tg_xxyyyzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 174); 

                auto tg_xxyyyzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 175); 

                auto tg_xxyyyzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 176); 

                auto tg_xxyyyzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 177); 

                auto tg_xxyyyzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 178); 

                auto tg_xxyyyzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 179); 

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

                auto tg_xxxxzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 54); 

                auto tg_xxxxzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 55); 

                auto tg_xxxxzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 56); 

                auto tg_xxxxzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 57); 

                auto tg_xxxxzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 58); 

                auto tg_xxxxzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 59); 

                auto tg_xxxyyyy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 60); 

                auto tg_xxxyyyy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 61); 

                auto tg_xxxyyyy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 62); 

                auto tg_xxxyyyy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 63); 

                auto tg_xxxyyyy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 64); 

                auto tg_xxxyyyy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 65); 

                auto tg_xxxyyyz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 66); 

                auto tg_xxxyyyz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 67); 

                auto tg_xxxyyyz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 68); 

                auto tg_xxxyyyz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 69); 

                auto tg_xxxyyyz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 70); 

                auto tg_xxxyyyz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 71); 

                auto tg_xxxyyzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 72); 

                auto tg_xxxyyzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 73); 

                auto tg_xxxyyzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 74); 

                auto tg_xxxyyzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 75); 

                auto tg_xxxyyzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 76); 

                auto tg_xxxyyzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 77); 

                auto tg_xxxyzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 78); 

                auto tg_xxxyzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 79); 

                auto tg_xxxyzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 80); 

                auto tg_xxxyzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 81); 

                auto tg_xxxyzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 82); 

                auto tg_xxxyzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 83); 

                auto tg_xxxzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 84); 

                auto tg_xxxzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 85); 

                auto tg_xxxzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 86); 

                auto tg_xxxzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 87); 

                auto tg_xxxzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 88); 

                auto tg_xxxzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 89); 

                auto tg_xxyyyyy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 90); 

                auto tg_xxyyyyy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 91); 

                auto tg_xxyyyyy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 92); 

                auto tg_xxyyyyy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 93); 

                auto tg_xxyyyyy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 94); 

                auto tg_xxyyyyy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 95); 

                auto tg_xxyyyyz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 96); 

                auto tg_xxyyyyz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 97); 

                auto tg_xxyyyyz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 98); 

                auto tg_xxyyyyz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 99); 

                auto tg_xxyyyyz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 100); 

                auto tg_xxyyyyz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 101); 

                auto tg_xxyyyzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 102); 

                auto tg_xxyyyzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 103); 

                auto tg_xxyyyzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 104); 

                auto tg_xxyyyzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 105); 

                auto tg_xxyyyzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 106); 

                auto tg_xxyyyzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 107); 

                // set up pointers to integrals

                auto tg_xxxxxzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 90); 

                auto tg_xxxxxzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 91); 

                auto tg_xxxxxzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 92); 

                auto tg_xxxxxzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 93); 

                auto tg_xxxxxzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 94); 

                auto tg_xxxxxzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 95); 

                auto tg_xxxxxzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 96); 

                auto tg_xxxxxzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 97); 

                auto tg_xxxxxzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 98); 

                auto tg_xxxxxzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 99); 

                auto tg_xxxxyyyy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 100); 

                auto tg_xxxxyyyy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 101); 

                auto tg_xxxxyyyy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 102); 

                auto tg_xxxxyyyy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 103); 

                auto tg_xxxxyyyy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 104); 

                auto tg_xxxxyyyy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 105); 

                auto tg_xxxxyyyy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 106); 

                auto tg_xxxxyyyy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 107); 

                auto tg_xxxxyyyy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 108); 

                auto tg_xxxxyyyy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 109); 

                auto tg_xxxxyyyz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 110); 

                auto tg_xxxxyyyz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 111); 

                auto tg_xxxxyyyz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 112); 

                auto tg_xxxxyyyz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 113); 

                auto tg_xxxxyyyz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 114); 

                auto tg_xxxxyyyz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 115); 

                auto tg_xxxxyyyz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 116); 

                auto tg_xxxxyyyz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 117); 

                auto tg_xxxxyyyz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 118); 

                auto tg_xxxxyyyz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 119); 

                auto tg_xxxxyyzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 120); 

                auto tg_xxxxyyzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 121); 

                auto tg_xxxxyyzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 122); 

                auto tg_xxxxyyzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 123); 

                auto tg_xxxxyyzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 124); 

                auto tg_xxxxyyzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 125); 

                auto tg_xxxxyyzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 126); 

                auto tg_xxxxyyzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 127); 

                auto tg_xxxxyyzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 128); 

                auto tg_xxxxyyzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 129); 

                auto tg_xxxxyzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 130); 

                auto tg_xxxxyzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 131); 

                auto tg_xxxxyzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 132); 

                auto tg_xxxxyzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 133); 

                auto tg_xxxxyzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 134); 

                auto tg_xxxxyzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 135); 

                auto tg_xxxxyzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 136); 

                auto tg_xxxxyzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 137); 

                auto tg_xxxxyzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 138); 

                auto tg_xxxxyzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 139); 

                auto tg_xxxxzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 140); 

                auto tg_xxxxzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 141); 

                auto tg_xxxxzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 142); 

                auto tg_xxxxzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 143); 

                auto tg_xxxxzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 144); 

                auto tg_xxxxzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 145); 

                auto tg_xxxxzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 146); 

                auto tg_xxxxzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 147); 

                auto tg_xxxxzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 148); 

                auto tg_xxxxzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 149); 

                auto tg_xxxyyyyy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 150); 

                auto tg_xxxyyyyy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 151); 

                auto tg_xxxyyyyy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 152); 

                auto tg_xxxyyyyy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 153); 

                auto tg_xxxyyyyy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 154); 

                auto tg_xxxyyyyy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 155); 

                auto tg_xxxyyyyy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 156); 

                auto tg_xxxyyyyy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 157); 

                auto tg_xxxyyyyy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 158); 

                auto tg_xxxyyyyy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 159); 

                auto tg_xxxyyyyz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 160); 

                auto tg_xxxyyyyz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 161); 

                auto tg_xxxyyyyz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 162); 

                auto tg_xxxyyyyz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 163); 

                auto tg_xxxyyyyz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 164); 

                auto tg_xxxyyyyz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 165); 

                auto tg_xxxyyyyz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 166); 

                auto tg_xxxyyyyz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 167); 

                auto tg_xxxyyyyz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 168); 

                auto tg_xxxyyyyz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 169); 

                auto tg_xxxyyyzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 170); 

                auto tg_xxxyyyzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 171); 

                auto tg_xxxyyyzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 172); 

                auto tg_xxxyyyzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 173); 

                auto tg_xxxyyyzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 174); 

                auto tg_xxxyyyzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 175); 

                auto tg_xxxyyyzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 176); 

                auto tg_xxxyyyzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 177); 

                auto tg_xxxyyyzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 178); 

                auto tg_xxxyyyzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 179); 

                // Batch of Integrals (90,180)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxzzz_xxx_0, tg_xxxxxzzz_xxy_0, tg_xxxxxzzz_xxz_0, \
                                         tg_xxxxxzzz_xyy_0, tg_xxxxxzzz_xyz_0, tg_xxxxxzzz_xzz_0, tg_xxxxxzzz_yyy_0, \
                                         tg_xxxxxzzz_yyz_0, tg_xxxxxzzz_yzz_0, tg_xxxxxzzz_zzz_0, tg_xxxxyyyy_xxx_0, \
                                         tg_xxxxyyyy_xxy_0, tg_xxxxyyyy_xxz_0, tg_xxxxyyyy_xyy_0, tg_xxxxyyyy_xyz_0, \
                                         tg_xxxxyyyy_xzz_0, tg_xxxxyyyy_yyy_0, tg_xxxxyyyy_yyz_0, tg_xxxxyyyy_yzz_0, \
                                         tg_xxxxyyyy_zzz_0, tg_xxxxyyyz_xxx_0, tg_xxxxyyyz_xxy_0, tg_xxxxyyyz_xxz_0, \
                                         tg_xxxxyyyz_xyy_0, tg_xxxxyyyz_xyz_0, tg_xxxxyyyz_xzz_0, tg_xxxxyyyz_yyy_0, \
                                         tg_xxxxyyyz_yyz_0, tg_xxxxyyyz_yzz_0, tg_xxxxyyyz_zzz_0, tg_xxxxyyzz_xxx_0, \
                                         tg_xxxxyyzz_xxy_0, tg_xxxxyyzz_xxz_0, tg_xxxxyyzz_xyy_0, tg_xxxxyyzz_xyz_0, \
                                         tg_xxxxyyzz_xzz_0, tg_xxxxyyzz_yyy_0, tg_xxxxyyzz_yyz_0, tg_xxxxyyzz_yzz_0, \
                                         tg_xxxxyyzz_zzz_0, tg_xxxxyzzz_xxx_0, tg_xxxxyzzz_xxy_0, tg_xxxxyzzz_xxz_0, \
                                         tg_xxxxyzzz_xyy_0, tg_xxxxyzzz_xyz_0, tg_xxxxyzzz_xzz_0, tg_xxxxyzzz_yyy_0, \
                                         tg_xxxxyzzz_yyz_0, tg_xxxxyzzz_yzz_0, tg_xxxxyzzz_zzz_0, tg_xxxxzzz_xx_1, \
                                         tg_xxxxzzz_xxx_0, tg_xxxxzzz_xxx_1, tg_xxxxzzz_xxy_0, tg_xxxxzzz_xxy_1, \
                                         tg_xxxxzzz_xxz_0, tg_xxxxzzz_xxz_1, tg_xxxxzzz_xy_1, tg_xxxxzzz_xyy_0, \
                                         tg_xxxxzzz_xyy_1, tg_xxxxzzz_xyz_0, tg_xxxxzzz_xyz_1, tg_xxxxzzz_xz_1, \
                                         tg_xxxxzzz_xzz_0, tg_xxxxzzz_xzz_1, tg_xxxxzzz_yy_1, tg_xxxxzzz_yyy_0, \
                                         tg_xxxxzzz_yyy_1, tg_xxxxzzz_yyz_0, tg_xxxxzzz_yyz_1, tg_xxxxzzz_yz_1, \
                                         tg_xxxxzzz_yzz_0, tg_xxxxzzz_yzz_1, tg_xxxxzzz_zz_1, tg_xxxxzzz_zzz_0, \
                                         tg_xxxxzzz_zzz_1, tg_xxxxzzzz_xxx_0, tg_xxxxzzzz_xxy_0, tg_xxxxzzzz_xxz_0, \
                                         tg_xxxxzzzz_xyy_0, tg_xxxxzzzz_xyz_0, tg_xxxxzzzz_xzz_0, tg_xxxxzzzz_yyy_0, \
                                         tg_xxxxzzzz_yyz_0, tg_xxxxzzzz_yzz_0, tg_xxxxzzzz_zzz_0, tg_xxxyyyy_xx_1, \
                                         tg_xxxyyyy_xxx_0, tg_xxxyyyy_xxx_1, tg_xxxyyyy_xxy_0, tg_xxxyyyy_xxy_1, \
                                         tg_xxxyyyy_xxz_0, tg_xxxyyyy_xxz_1, tg_xxxyyyy_xy_1, tg_xxxyyyy_xyy_0, \
                                         tg_xxxyyyy_xyy_1, tg_xxxyyyy_xyz_0, tg_xxxyyyy_xyz_1, tg_xxxyyyy_xz_1, \
                                         tg_xxxyyyy_xzz_0, tg_xxxyyyy_xzz_1, tg_xxxyyyy_yy_1, tg_xxxyyyy_yyy_0, \
                                         tg_xxxyyyy_yyy_1, tg_xxxyyyy_yyz_0, tg_xxxyyyy_yyz_1, tg_xxxyyyy_yz_1, \
                                         tg_xxxyyyy_yzz_0, tg_xxxyyyy_yzz_1, tg_xxxyyyy_zz_1, tg_xxxyyyy_zzz_0, \
                                         tg_xxxyyyy_zzz_1, tg_xxxyyyyy_xxx_0, tg_xxxyyyyy_xxy_0, tg_xxxyyyyy_xxz_0, \
                                         tg_xxxyyyyy_xyy_0, tg_xxxyyyyy_xyz_0, tg_xxxyyyyy_xzz_0, tg_xxxyyyyy_yyy_0, \
                                         tg_xxxyyyyy_yyz_0, tg_xxxyyyyy_yzz_0, tg_xxxyyyyy_zzz_0, tg_xxxyyyyz_xxx_0, \
                                         tg_xxxyyyyz_xxy_0, tg_xxxyyyyz_xxz_0, tg_xxxyyyyz_xyy_0, tg_xxxyyyyz_xyz_0, \
                                         tg_xxxyyyyz_xzz_0, tg_xxxyyyyz_yyy_0, tg_xxxyyyyz_yyz_0, tg_xxxyyyyz_yzz_0, \
                                         tg_xxxyyyyz_zzz_0, tg_xxxyyyz_xx_1, tg_xxxyyyz_xxx_0, tg_xxxyyyz_xxx_1, \
                                         tg_xxxyyyz_xxy_0, tg_xxxyyyz_xxy_1, tg_xxxyyyz_xxz_0, tg_xxxyyyz_xxz_1, \
                                         tg_xxxyyyz_xy_1, tg_xxxyyyz_xyy_0, tg_xxxyyyz_xyy_1, tg_xxxyyyz_xyz_0, \
                                         tg_xxxyyyz_xyz_1, tg_xxxyyyz_xz_1, tg_xxxyyyz_xzz_0, tg_xxxyyyz_xzz_1, \
                                         tg_xxxyyyz_yy_1, tg_xxxyyyz_yyy_0, tg_xxxyyyz_yyy_1, tg_xxxyyyz_yyz_0, \
                                         tg_xxxyyyz_yyz_1, tg_xxxyyyz_yz_1, tg_xxxyyyz_yzz_0, tg_xxxyyyz_yzz_1, \
                                         tg_xxxyyyz_zz_1, tg_xxxyyyz_zzz_0, tg_xxxyyyz_zzz_1, tg_xxxyyyzz_xxx_0, \
                                         tg_xxxyyyzz_xxy_0, tg_xxxyyyzz_xxz_0, tg_xxxyyyzz_xyy_0, tg_xxxyyyzz_xyz_0, \
                                         tg_xxxyyyzz_xzz_0, tg_xxxyyyzz_yyy_0, tg_xxxyyyzz_yyz_0, tg_xxxyyyzz_yzz_0, \
                                         tg_xxxyyyzz_zzz_0, tg_xxxyyzz_xx_1, tg_xxxyyzz_xxx_0, tg_xxxyyzz_xxx_1, \
                                         tg_xxxyyzz_xxy_0, tg_xxxyyzz_xxy_1, tg_xxxyyzz_xxz_0, tg_xxxyyzz_xxz_1, \
                                         tg_xxxyyzz_xy_1, tg_xxxyyzz_xyy_0, tg_xxxyyzz_xyy_1, tg_xxxyyzz_xyz_0, \
                                         tg_xxxyyzz_xyz_1, tg_xxxyyzz_xz_1, tg_xxxyyzz_xzz_0, tg_xxxyyzz_xzz_1, \
                                         tg_xxxyyzz_yy_1, tg_xxxyyzz_yyy_0, tg_xxxyyzz_yyy_1, tg_xxxyyzz_yyz_0, \
                                         tg_xxxyyzz_yyz_1, tg_xxxyyzz_yz_1, tg_xxxyyzz_yzz_0, tg_xxxyyzz_yzz_1, \
                                         tg_xxxyyzz_zz_1, tg_xxxyyzz_zzz_0, tg_xxxyyzz_zzz_1, tg_xxxyzzz_xx_1, \
                                         tg_xxxyzzz_xxx_0, tg_xxxyzzz_xxx_1, tg_xxxyzzz_xxy_0, tg_xxxyzzz_xxy_1, \
                                         tg_xxxyzzz_xxz_0, tg_xxxyzzz_xxz_1, tg_xxxyzzz_xy_1, tg_xxxyzzz_xyy_0, \
                                         tg_xxxyzzz_xyy_1, tg_xxxyzzz_xyz_0, tg_xxxyzzz_xyz_1, tg_xxxyzzz_xz_1, \
                                         tg_xxxyzzz_xzz_0, tg_xxxyzzz_xzz_1, tg_xxxyzzz_yy_1, tg_xxxyzzz_yyy_0, \
                                         tg_xxxyzzz_yyy_1, tg_xxxyzzz_yyz_0, tg_xxxyzzz_yyz_1, tg_xxxyzzz_yz_1, \
                                         tg_xxxyzzz_yzz_0, tg_xxxyzzz_yzz_1, tg_xxxyzzz_zz_1, tg_xxxyzzz_zzz_0, \
                                         tg_xxxyzzz_zzz_1, tg_xxxzzz_xxx_0, tg_xxxzzz_xxx_1, tg_xxxzzz_xxy_0, tg_xxxzzz_xxy_1, \
                                         tg_xxxzzz_xxz_0, tg_xxxzzz_xxz_1, tg_xxxzzz_xyy_0, tg_xxxzzz_xyy_1, tg_xxxzzz_xyz_0, \
                                         tg_xxxzzz_xyz_1, tg_xxxzzz_xzz_0, tg_xxxzzz_xzz_1, tg_xxxzzz_yyy_0, tg_xxxzzz_yyy_1, \
                                         tg_xxxzzz_yyz_0, tg_xxxzzz_yyz_1, tg_xxxzzz_yzz_0, tg_xxxzzz_yzz_1, tg_xxxzzz_zzz_0, \
                                         tg_xxxzzz_zzz_1, tg_xxxzzzz_xx_1, tg_xxxzzzz_xxx_0, tg_xxxzzzz_xxx_1, \
                                         tg_xxxzzzz_xxy_0, tg_xxxzzzz_xxy_1, tg_xxxzzzz_xxz_0, tg_xxxzzzz_xxz_1, \
                                         tg_xxxzzzz_xy_1, tg_xxxzzzz_xyy_0, tg_xxxzzzz_xyy_1, tg_xxxzzzz_xyz_0, \
                                         tg_xxxzzzz_xyz_1, tg_xxxzzzz_xz_1, tg_xxxzzzz_xzz_0, tg_xxxzzzz_xzz_1, \
                                         tg_xxxzzzz_yy_1, tg_xxxzzzz_yyy_0, tg_xxxzzzz_yyy_1, tg_xxxzzzz_yyz_0, \
                                         tg_xxxzzzz_yyz_1, tg_xxxzzzz_yz_1, tg_xxxzzzz_yzz_0, tg_xxxzzzz_yzz_1, \
                                         tg_xxxzzzz_zz_1, tg_xxxzzzz_zzz_0, tg_xxxzzzz_zzz_1, tg_xxyyyy_xxx_0, \
                                         tg_xxyyyy_xxx_1, tg_xxyyyy_xxy_0, tg_xxyyyy_xxy_1, tg_xxyyyy_xxz_0, tg_xxyyyy_xxz_1, \
                                         tg_xxyyyy_xyy_0, tg_xxyyyy_xyy_1, tg_xxyyyy_xyz_0, tg_xxyyyy_xyz_1, tg_xxyyyy_xzz_0, \
                                         tg_xxyyyy_xzz_1, tg_xxyyyy_yyy_0, tg_xxyyyy_yyy_1, tg_xxyyyy_yyz_0, tg_xxyyyy_yyz_1, \
                                         tg_xxyyyy_yzz_0, tg_xxyyyy_yzz_1, tg_xxyyyy_zzz_0, tg_xxyyyy_zzz_1, tg_xxyyyyy_xx_1, \
                                         tg_xxyyyyy_xxx_0, tg_xxyyyyy_xxx_1, tg_xxyyyyy_xxy_0, tg_xxyyyyy_xxy_1, \
                                         tg_xxyyyyy_xxz_0, tg_xxyyyyy_xxz_1, tg_xxyyyyy_xy_1, tg_xxyyyyy_xyy_0, \
                                         tg_xxyyyyy_xyy_1, tg_xxyyyyy_xyz_0, tg_xxyyyyy_xyz_1, tg_xxyyyyy_xz_1, \
                                         tg_xxyyyyy_xzz_0, tg_xxyyyyy_xzz_1, tg_xxyyyyy_yy_1, tg_xxyyyyy_yyy_0, \
                                         tg_xxyyyyy_yyy_1, tg_xxyyyyy_yyz_0, tg_xxyyyyy_yyz_1, tg_xxyyyyy_yz_1, \
                                         tg_xxyyyyy_yzz_0, tg_xxyyyyy_yzz_1, tg_xxyyyyy_zz_1, tg_xxyyyyy_zzz_0, \
                                         tg_xxyyyyy_zzz_1, tg_xxyyyyz_xx_1, tg_xxyyyyz_xxx_0, tg_xxyyyyz_xxx_1, \
                                         tg_xxyyyyz_xxy_0, tg_xxyyyyz_xxy_1, tg_xxyyyyz_xxz_0, tg_xxyyyyz_xxz_1, \
                                         tg_xxyyyyz_xy_1, tg_xxyyyyz_xyy_0, tg_xxyyyyz_xyy_1, tg_xxyyyyz_xyz_0, \
                                         tg_xxyyyyz_xyz_1, tg_xxyyyyz_xz_1, tg_xxyyyyz_xzz_0, tg_xxyyyyz_xzz_1, \
                                         tg_xxyyyyz_yy_1, tg_xxyyyyz_yyy_0, tg_xxyyyyz_yyy_1, tg_xxyyyyz_yyz_0, \
                                         tg_xxyyyyz_yyz_1, tg_xxyyyyz_yz_1, tg_xxyyyyz_yzz_0, tg_xxyyyyz_yzz_1, \
                                         tg_xxyyyyz_zz_1, tg_xxyyyyz_zzz_0, tg_xxyyyyz_zzz_1, tg_xxyyyz_xxx_0, \
                                         tg_xxyyyz_xxx_1, tg_xxyyyz_xxy_0, tg_xxyyyz_xxy_1, tg_xxyyyz_xxz_0, tg_xxyyyz_xxz_1, \
                                         tg_xxyyyz_xyy_0, tg_xxyyyz_xyy_1, tg_xxyyyz_xyz_0, tg_xxyyyz_xyz_1, tg_xxyyyz_xzz_0, \
                                         tg_xxyyyz_xzz_1, tg_xxyyyz_yyy_0, tg_xxyyyz_yyy_1, tg_xxyyyz_yyz_0, tg_xxyyyz_yyz_1, \
                                         tg_xxyyyz_yzz_0, tg_xxyyyz_yzz_1, tg_xxyyyz_zzz_0, tg_xxyyyz_zzz_1, tg_xxyyyzz_xx_1, \
                                         tg_xxyyyzz_xxx_0, tg_xxyyyzz_xxx_1, tg_xxyyyzz_xxy_0, tg_xxyyyzz_xxy_1, \
                                         tg_xxyyyzz_xxz_0, tg_xxyyyzz_xxz_1, tg_xxyyyzz_xy_1, tg_xxyyyzz_xyy_0, \
                                         tg_xxyyyzz_xyy_1, tg_xxyyyzz_xyz_0, tg_xxyyyzz_xyz_1, tg_xxyyyzz_xz_1, \
                                         tg_xxyyyzz_xzz_0, tg_xxyyyzz_xzz_1, tg_xxyyyzz_yy_1, tg_xxyyyzz_yyy_0, \
                                         tg_xxyyyzz_yyy_1, tg_xxyyyzz_yyz_0, tg_xxyyyzz_yyz_1, tg_xxyyyzz_yz_1, \
                                         tg_xxyyyzz_yzz_0, tg_xxyyyzz_yzz_1, tg_xxyyyzz_zz_1, tg_xxyyyzz_zzz_0, \
                                         tg_xxyyyzz_zzz_1, tg_xxyyzz_xxx_0, tg_xxyyzz_xxx_1, tg_xxyyzz_xxy_0, tg_xxyyzz_xxy_1, \
                                         tg_xxyyzz_xxz_0, tg_xxyyzz_xxz_1, tg_xxyyzz_xyy_0, tg_xxyyzz_xyy_1, tg_xxyyzz_xyz_0, \
                                         tg_xxyyzz_xyz_1, tg_xxyyzz_xzz_0, tg_xxyyzz_xzz_1, tg_xxyyzz_yyy_0, tg_xxyyzz_yyy_1, \
                                         tg_xxyyzz_yyz_0, tg_xxyyzz_yyz_1, tg_xxyyzz_yzz_0, tg_xxyyzz_yzz_1, tg_xxyyzz_zzz_0, \
                                         tg_xxyyzz_zzz_1, tg_xxyzzz_xxx_0, tg_xxyzzz_xxx_1, tg_xxyzzz_xxy_0, tg_xxyzzz_xxy_1, \
                                         tg_xxyzzz_xxz_0, tg_xxyzzz_xxz_1, tg_xxyzzz_xyy_0, tg_xxyzzz_xyy_1, tg_xxyzzz_xyz_0, \
                                         tg_xxyzzz_xyz_1, tg_xxyzzz_xzz_0, tg_xxyzzz_xzz_1, tg_xxyzzz_yyy_0, tg_xxyzzz_yyy_1, \
                                         tg_xxyzzz_yyz_0, tg_xxyzzz_yyz_1, tg_xxyzzz_yzz_0, tg_xxyzzz_yzz_1, tg_xxyzzz_zzz_0, \
                                         tg_xxyzzz_zzz_1, tg_xxzzzz_xxx_0, tg_xxzzzz_xxx_1, tg_xxzzzz_xxy_0, tg_xxzzzz_xxy_1, \
                                         tg_xxzzzz_xxz_0, tg_xxzzzz_xxz_1, tg_xxzzzz_xyy_0, tg_xxzzzz_xyy_1, tg_xxzzzz_xyz_0, \
                                         tg_xxzzzz_xyz_1, tg_xxzzzz_xzz_0, tg_xxzzzz_xzz_1, tg_xxzzzz_yyy_0, tg_xxzzzz_yyy_1, \
                                         tg_xxzzzz_yyz_0, tg_xxzzzz_yyz_1, tg_xxzzzz_yzz_0, tg_xxzzzz_yzz_1, tg_xxzzzz_zzz_0, \
                                         tg_xxzzzz_zzz_1, tg_xyyyyy_xxx_0, tg_xyyyyy_xxx_1, tg_xyyyyy_xxy_0, tg_xyyyyy_xxy_1, \
                                         tg_xyyyyy_xxz_0, tg_xyyyyy_xxz_1, tg_xyyyyy_xyy_0, tg_xyyyyy_xyy_1, tg_xyyyyy_xyz_0, \
                                         tg_xyyyyy_xyz_1, tg_xyyyyy_xzz_0, tg_xyyyyy_xzz_1, tg_xyyyyy_yyy_0, tg_xyyyyy_yyy_1, \
                                         tg_xyyyyy_yyz_0, tg_xyyyyy_yyz_1, tg_xyyyyy_yzz_0, tg_xyyyyy_yzz_1, tg_xyyyyy_zzz_0, \
                                         tg_xyyyyy_zzz_1, tg_xyyyyz_xxx_0, tg_xyyyyz_xxx_1, tg_xyyyyz_xxy_0, tg_xyyyyz_xxy_1, \
                                         tg_xyyyyz_xxz_0, tg_xyyyyz_xxz_1, tg_xyyyyz_xyy_0, tg_xyyyyz_xyy_1, tg_xyyyyz_xyz_0, \
                                         tg_xyyyyz_xyz_1, tg_xyyyyz_xzz_0, tg_xyyyyz_xzz_1, tg_xyyyyz_yyy_0, tg_xyyyyz_yyy_1, \
                                         tg_xyyyyz_yyz_0, tg_xyyyyz_yyz_1, tg_xyyyyz_yzz_0, tg_xyyyyz_yzz_1, tg_xyyyyz_zzz_0, \
                                         tg_xyyyyz_zzz_1, tg_xyyyzz_xxx_0, tg_xyyyzz_xxx_1, tg_xyyyzz_xxy_0, tg_xyyyzz_xxy_1, \
                                         tg_xyyyzz_xxz_0, tg_xyyyzz_xxz_1, tg_xyyyzz_xyy_0, tg_xyyyzz_xyy_1, tg_xyyyzz_xyz_0, \
                                         tg_xyyyzz_xyz_1, tg_xyyyzz_xzz_0, tg_xyyyzz_xzz_1, tg_xyyyzz_yyy_0, tg_xyyyzz_yyy_1, \
                                         tg_xyyyzz_yyz_0, tg_xyyyzz_yyz_1, tg_xyyyzz_yzz_0, tg_xyyyzz_yzz_1, tg_xyyyzz_zzz_0, \
                                         tg_xyyyzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxxxzzz_xxx_0[j] = pb_x * tg_xxxxzzz_xxx_0[j] + fr * tg_xxxxzzz_xxx_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxx_0[j] - tg_xxxzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxxzzz_xx_1[j];

                    tg_xxxxxzzz_xxy_0[j] = pb_x * tg_xxxxzzz_xxy_0[j] + fr * tg_xxxxzzz_xxy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxy_0[j] - tg_xxxzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzzz_xy_1[j];

                    tg_xxxxxzzz_xxz_0[j] = pb_x * tg_xxxxzzz_xxz_0[j] + fr * tg_xxxxzzz_xxz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xxz_0[j] - tg_xxxzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxxzzz_xz_1[j];

                    tg_xxxxxzzz_xyy_0[j] = pb_x * tg_xxxxzzz_xyy_0[j] + fr * tg_xxxxzzz_xyy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xyy_0[j] - tg_xxxzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_yy_1[j];

                    tg_xxxxxzzz_xyz_0[j] = pb_x * tg_xxxxzzz_xyz_0[j] + fr * tg_xxxxzzz_xyz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xyz_0[j] - tg_xxxzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_yz_1[j];

                    tg_xxxxxzzz_xzz_0[j] = pb_x * tg_xxxxzzz_xzz_0[j] + fr * tg_xxxxzzz_xzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_xzz_0[j] - tg_xxxzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxxzzz_zz_1[j];

                    tg_xxxxxzzz_yyy_0[j] = pb_x * tg_xxxxzzz_yyy_0[j] + fr * tg_xxxxzzz_yyy_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yyy_0[j] - tg_xxxzzz_yyy_1[j] * fl1_fza);

                    tg_xxxxxzzz_yyz_0[j] = pb_x * tg_xxxxzzz_yyz_0[j] + fr * tg_xxxxzzz_yyz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yyz_0[j] - tg_xxxzzz_yyz_1[j] * fl1_fza);

                    tg_xxxxxzzz_yzz_0[j] = pb_x * tg_xxxxzzz_yzz_0[j] + fr * tg_xxxxzzz_yzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_yzz_0[j] - tg_xxxzzz_yzz_1[j] * fl1_fza);

                    tg_xxxxxzzz_zzz_0[j] = pb_x * tg_xxxxzzz_zzz_0[j] + fr * tg_xxxxzzz_zzz_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_zzz_0[j] - tg_xxxzzz_zzz_1[j] * fl1_fza);

                    tg_xxxxyyyy_xxx_0[j] = pb_x * tg_xxxyyyy_xxx_0[j] + fr * tg_xxxyyyy_xxx_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxx_0[j] - tg_xxyyyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyy_xx_1[j];

                    tg_xxxxyyyy_xxy_0[j] = pb_x * tg_xxxyyyy_xxy_0[j] + fr * tg_xxxyyyy_xxy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxy_0[j] - tg_xxyyyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyy_xy_1[j];

                    tg_xxxxyyyy_xxz_0[j] = pb_x * tg_xxxyyyy_xxz_0[j] + fr * tg_xxxyyyy_xxz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xxz_0[j] - tg_xxyyyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyy_xz_1[j];

                    tg_xxxxyyyy_xyy_0[j] = pb_x * tg_xxxyyyy_xyy_0[j] + fr * tg_xxxyyyy_xyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xyy_0[j] - tg_xxyyyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_yy_1[j];

                    tg_xxxxyyyy_xyz_0[j] = pb_x * tg_xxxyyyy_xyz_0[j] + fr * tg_xxxyyyy_xyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xyz_0[j] - tg_xxyyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_yz_1[j];

                    tg_xxxxyyyy_xzz_0[j] = pb_x * tg_xxxyyyy_xzz_0[j] + fr * tg_xxxyyyy_xzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_xzz_0[j] - tg_xxyyyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyy_zz_1[j];

                    tg_xxxxyyyy_yyy_0[j] = pb_x * tg_xxxyyyy_yyy_0[j] + fr * tg_xxxyyyy_yyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yyy_0[j] - tg_xxyyyy_yyy_1[j] * fl1_fza);

                    tg_xxxxyyyy_yyz_0[j] = pb_x * tg_xxxyyyy_yyz_0[j] + fr * tg_xxxyyyy_yyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yyz_0[j] - tg_xxyyyy_yyz_1[j] * fl1_fza);

                    tg_xxxxyyyy_yzz_0[j] = pb_x * tg_xxxyyyy_yzz_0[j] + fr * tg_xxxyyyy_yzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_yzz_0[j] - tg_xxyyyy_yzz_1[j] * fl1_fza);

                    tg_xxxxyyyy_zzz_0[j] = pb_x * tg_xxxyyyy_zzz_0[j] + fr * tg_xxxyyyy_zzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_zzz_0[j] - tg_xxyyyy_zzz_1[j] * fl1_fza);

                    tg_xxxxyyyz_xxx_0[j] = pb_x * tg_xxxyyyz_xxx_0[j] + fr * tg_xxxyyyz_xxx_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxx_0[j] - tg_xxyyyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyyz_xx_1[j];

                    tg_xxxxyyyz_xxy_0[j] = pb_x * tg_xxxyyyz_xxy_0[j] + fr * tg_xxxyyyz_xxy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxy_0[j] - tg_xxyyyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyz_xy_1[j];

                    tg_xxxxyyyz_xxz_0[j] = pb_x * tg_xxxyyyz_xxz_0[j] + fr * tg_xxxyyyz_xxz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xxz_0[j] - tg_xxyyyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyyz_xz_1[j];

                    tg_xxxxyyyz_xyy_0[j] = pb_x * tg_xxxyyyz_xyy_0[j] + fr * tg_xxxyyyz_xyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xyy_0[j] - tg_xxyyyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_yy_1[j];

                    tg_xxxxyyyz_xyz_0[j] = pb_x * tg_xxxyyyz_xyz_0[j] + fr * tg_xxxyyyz_xyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xyz_0[j] - tg_xxyyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_yz_1[j];

                    tg_xxxxyyyz_xzz_0[j] = pb_x * tg_xxxyyyz_xzz_0[j] + fr * tg_xxxyyyz_xzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_xzz_0[j] - tg_xxyyyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyyz_zz_1[j];

                    tg_xxxxyyyz_yyy_0[j] = pb_x * tg_xxxyyyz_yyy_0[j] + fr * tg_xxxyyyz_yyy_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yyy_0[j] - tg_xxyyyz_yyy_1[j] * fl1_fza);

                    tg_xxxxyyyz_yyz_0[j] = pb_x * tg_xxxyyyz_yyz_0[j] + fr * tg_xxxyyyz_yyz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yyz_0[j] - tg_xxyyyz_yyz_1[j] * fl1_fza);

                    tg_xxxxyyyz_yzz_0[j] = pb_x * tg_xxxyyyz_yzz_0[j] + fr * tg_xxxyyyz_yzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_yzz_0[j] - tg_xxyyyz_yzz_1[j] * fl1_fza);

                    tg_xxxxyyyz_zzz_0[j] = pb_x * tg_xxxyyyz_zzz_0[j] + fr * tg_xxxyyyz_zzz_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_zzz_0[j] - tg_xxyyyz_zzz_1[j] * fl1_fza);

                    tg_xxxxyyzz_xxx_0[j] = pb_x * tg_xxxyyzz_xxx_0[j] + fr * tg_xxxyyzz_xxx_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxx_0[j] - tg_xxyyzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyyzz_xx_1[j];

                    tg_xxxxyyzz_xxy_0[j] = pb_x * tg_xxxyyzz_xxy_0[j] + fr * tg_xxxyyzz_xxy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxy_0[j] - tg_xxyyzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyzz_xy_1[j];

                    tg_xxxxyyzz_xxz_0[j] = pb_x * tg_xxxyyzz_xxz_0[j] + fr * tg_xxxyyzz_xxz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xxz_0[j] - tg_xxyyzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyyzz_xz_1[j];

                    tg_xxxxyyzz_xyy_0[j] = pb_x * tg_xxxyyzz_xyy_0[j] + fr * tg_xxxyyzz_xyy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xyy_0[j] - tg_xxyyzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_yy_1[j];

                    tg_xxxxyyzz_xyz_0[j] = pb_x * tg_xxxyyzz_xyz_0[j] + fr * tg_xxxyyzz_xyz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xyz_0[j] - tg_xxyyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_yz_1[j];

                    tg_xxxxyyzz_xzz_0[j] = pb_x * tg_xxxyyzz_xzz_0[j] + fr * tg_xxxyyzz_xzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_xzz_0[j] - tg_xxyyzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyyzz_zz_1[j];

                    tg_xxxxyyzz_yyy_0[j] = pb_x * tg_xxxyyzz_yyy_0[j] + fr * tg_xxxyyzz_yyy_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yyy_0[j] - tg_xxyyzz_yyy_1[j] * fl1_fza);

                    tg_xxxxyyzz_yyz_0[j] = pb_x * tg_xxxyyzz_yyz_0[j] + fr * tg_xxxyyzz_yyz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yyz_0[j] - tg_xxyyzz_yyz_1[j] * fl1_fza);

                    tg_xxxxyyzz_yzz_0[j] = pb_x * tg_xxxyyzz_yzz_0[j] + fr * tg_xxxyyzz_yzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_yzz_0[j] - tg_xxyyzz_yzz_1[j] * fl1_fza);

                    tg_xxxxyyzz_zzz_0[j] = pb_x * tg_xxxyyzz_zzz_0[j] + fr * tg_xxxyyzz_zzz_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_zzz_0[j] - tg_xxyyzz_zzz_1[j] * fl1_fza);

                    tg_xxxxyzzz_xxx_0[j] = pb_x * tg_xxxyzzz_xxx_0[j] + fr * tg_xxxyzzz_xxx_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxx_0[j] - tg_xxyzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxyzzz_xx_1[j];

                    tg_xxxxyzzz_xxy_0[j] = pb_x * tg_xxxyzzz_xxy_0[j] + fr * tg_xxxyzzz_xxy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxy_0[j] - tg_xxyzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzzz_xy_1[j];

                    tg_xxxxyzzz_xxz_0[j] = pb_x * tg_xxxyzzz_xxz_0[j] + fr * tg_xxxyzzz_xxz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xxz_0[j] - tg_xxyzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxyzzz_xz_1[j];

                    tg_xxxxyzzz_xyy_0[j] = pb_x * tg_xxxyzzz_xyy_0[j] + fr * tg_xxxyzzz_xyy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xyy_0[j] - tg_xxyzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_yy_1[j];

                    tg_xxxxyzzz_xyz_0[j] = pb_x * tg_xxxyzzz_xyz_0[j] + fr * tg_xxxyzzz_xyz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xyz_0[j] - tg_xxyzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_yz_1[j];

                    tg_xxxxyzzz_xzz_0[j] = pb_x * tg_xxxyzzz_xzz_0[j] + fr * tg_xxxyzzz_xzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_xzz_0[j] - tg_xxyzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxyzzz_zz_1[j];

                    tg_xxxxyzzz_yyy_0[j] = pb_x * tg_xxxyzzz_yyy_0[j] + fr * tg_xxxyzzz_yyy_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yyy_0[j] - tg_xxyzzz_yyy_1[j] * fl1_fza);

                    tg_xxxxyzzz_yyz_0[j] = pb_x * tg_xxxyzzz_yyz_0[j] + fr * tg_xxxyzzz_yyz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yyz_0[j] - tg_xxyzzz_yyz_1[j] * fl1_fza);

                    tg_xxxxyzzz_yzz_0[j] = pb_x * tg_xxxyzzz_yzz_0[j] + fr * tg_xxxyzzz_yzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_yzz_0[j] - tg_xxyzzz_yzz_1[j] * fl1_fza);

                    tg_xxxxyzzz_zzz_0[j] = pb_x * tg_xxxyzzz_zzz_0[j] + fr * tg_xxxyzzz_zzz_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_zzz_0[j] - tg_xxyzzz_zzz_1[j] * fl1_fza);

                    tg_xxxxzzzz_xxx_0[j] = pb_x * tg_xxxzzzz_xxx_0[j] + fr * tg_xxxzzzz_xxx_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxx_0[j] - tg_xxzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxxzzzz_xx_1[j];

                    tg_xxxxzzzz_xxy_0[j] = pb_x * tg_xxxzzzz_xxy_0[j] + fr * tg_xxxzzzz_xxy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxy_0[j] - tg_xxzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzzz_xy_1[j];

                    tg_xxxxzzzz_xxz_0[j] = pb_x * tg_xxxzzzz_xxz_0[j] + fr * tg_xxxzzzz_xxz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xxz_0[j] - tg_xxzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxxzzzz_xz_1[j];

                    tg_xxxxzzzz_xyy_0[j] = pb_x * tg_xxxzzzz_xyy_0[j] + fr * tg_xxxzzzz_xyy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xyy_0[j] - tg_xxzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_yy_1[j];

                    tg_xxxxzzzz_xyz_0[j] = pb_x * tg_xxxzzzz_xyz_0[j] + fr * tg_xxxzzzz_xyz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xyz_0[j] - tg_xxzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_yz_1[j];

                    tg_xxxxzzzz_xzz_0[j] = pb_x * tg_xxxzzzz_xzz_0[j] + fr * tg_xxxzzzz_xzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_xzz_0[j] - tg_xxzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxxzzzz_zz_1[j];

                    tg_xxxxzzzz_yyy_0[j] = pb_x * tg_xxxzzzz_yyy_0[j] + fr * tg_xxxzzzz_yyy_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yyy_0[j] - tg_xxzzzz_yyy_1[j] * fl1_fza);

                    tg_xxxxzzzz_yyz_0[j] = pb_x * tg_xxxzzzz_yyz_0[j] + fr * tg_xxxzzzz_yyz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yyz_0[j] - tg_xxzzzz_yyz_1[j] * fl1_fza);

                    tg_xxxxzzzz_yzz_0[j] = pb_x * tg_xxxzzzz_yzz_0[j] + fr * tg_xxxzzzz_yzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_yzz_0[j] - tg_xxzzzz_yzz_1[j] * fl1_fza);

                    tg_xxxxzzzz_zzz_0[j] = pb_x * tg_xxxzzzz_zzz_0[j] + fr * tg_xxxzzzz_zzz_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_zzz_0[j] - tg_xxzzzz_zzz_1[j] * fl1_fza);

                    tg_xxxyyyyy_xxx_0[j] = pb_x * tg_xxyyyyy_xxx_0[j] + fr * tg_xxyyyyy_xxx_1[j] + fl1_fx * (tg_xyyyyy_xxx_0[j] - tg_xyyyyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyy_xx_1[j];

                    tg_xxxyyyyy_xxy_0[j] = pb_x * tg_xxyyyyy_xxy_0[j] + fr * tg_xxyyyyy_xxy_1[j] + fl1_fx * (tg_xyyyyy_xxy_0[j] - tg_xyyyyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyy_xy_1[j];

                    tg_xxxyyyyy_xxz_0[j] = pb_x * tg_xxyyyyy_xxz_0[j] + fr * tg_xxyyyyy_xxz_1[j] + fl1_fx * (tg_xyyyyy_xxz_0[j] - tg_xyyyyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyy_xz_1[j];

                    tg_xxxyyyyy_xyy_0[j] = pb_x * tg_xxyyyyy_xyy_0[j] + fr * tg_xxyyyyy_xyy_1[j] + fl1_fx * (tg_xyyyyy_xyy_0[j] - tg_xyyyyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_yy_1[j];

                    tg_xxxyyyyy_xyz_0[j] = pb_x * tg_xxyyyyy_xyz_0[j] + fr * tg_xxyyyyy_xyz_1[j] + fl1_fx * (tg_xyyyyy_xyz_0[j] - tg_xyyyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_yz_1[j];

                    tg_xxxyyyyy_xzz_0[j] = pb_x * tg_xxyyyyy_xzz_0[j] + fr * tg_xxyyyyy_xzz_1[j] + fl1_fx * (tg_xyyyyy_xzz_0[j] - tg_xyyyyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyy_zz_1[j];

                    tg_xxxyyyyy_yyy_0[j] = pb_x * tg_xxyyyyy_yyy_0[j] + fr * tg_xxyyyyy_yyy_1[j] + fl1_fx * (tg_xyyyyy_yyy_0[j] - tg_xyyyyy_yyy_1[j] * fl1_fza);

                    tg_xxxyyyyy_yyz_0[j] = pb_x * tg_xxyyyyy_yyz_0[j] + fr * tg_xxyyyyy_yyz_1[j] + fl1_fx * (tg_xyyyyy_yyz_0[j] - tg_xyyyyy_yyz_1[j] * fl1_fza);

                    tg_xxxyyyyy_yzz_0[j] = pb_x * tg_xxyyyyy_yzz_0[j] + fr * tg_xxyyyyy_yzz_1[j] + fl1_fx * (tg_xyyyyy_yzz_0[j] - tg_xyyyyy_yzz_1[j] * fl1_fza);

                    tg_xxxyyyyy_zzz_0[j] = pb_x * tg_xxyyyyy_zzz_0[j] + fr * tg_xxyyyyy_zzz_1[j] + fl1_fx * (tg_xyyyyy_zzz_0[j] - tg_xyyyyy_zzz_1[j] * fl1_fza);

                    tg_xxxyyyyz_xxx_0[j] = pb_x * tg_xxyyyyz_xxx_0[j] + fr * tg_xxyyyyz_xxx_1[j] + fl1_fx * (tg_xyyyyz_xxx_0[j] - tg_xyyyyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyyz_xx_1[j];

                    tg_xxxyyyyz_xxy_0[j] = pb_x * tg_xxyyyyz_xxy_0[j] + fr * tg_xxyyyyz_xxy_1[j] + fl1_fx * (tg_xyyyyz_xxy_0[j] - tg_xyyyyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyz_xy_1[j];

                    tg_xxxyyyyz_xxz_0[j] = pb_x * tg_xxyyyyz_xxz_0[j] + fr * tg_xxyyyyz_xxz_1[j] + fl1_fx * (tg_xyyyyz_xxz_0[j] - tg_xyyyyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyyz_xz_1[j];

                    tg_xxxyyyyz_xyy_0[j] = pb_x * tg_xxyyyyz_xyy_0[j] + fr * tg_xxyyyyz_xyy_1[j] + fl1_fx * (tg_xyyyyz_xyy_0[j] - tg_xyyyyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_yy_1[j];

                    tg_xxxyyyyz_xyz_0[j] = pb_x * tg_xxyyyyz_xyz_0[j] + fr * tg_xxyyyyz_xyz_1[j] + fl1_fx * (tg_xyyyyz_xyz_0[j] - tg_xyyyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_yz_1[j];

                    tg_xxxyyyyz_xzz_0[j] = pb_x * tg_xxyyyyz_xzz_0[j] + fr * tg_xxyyyyz_xzz_1[j] + fl1_fx * (tg_xyyyyz_xzz_0[j] - tg_xyyyyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyyz_zz_1[j];

                    tg_xxxyyyyz_yyy_0[j] = pb_x * tg_xxyyyyz_yyy_0[j] + fr * tg_xxyyyyz_yyy_1[j] + fl1_fx * (tg_xyyyyz_yyy_0[j] - tg_xyyyyz_yyy_1[j] * fl1_fza);

                    tg_xxxyyyyz_yyz_0[j] = pb_x * tg_xxyyyyz_yyz_0[j] + fr * tg_xxyyyyz_yyz_1[j] + fl1_fx * (tg_xyyyyz_yyz_0[j] - tg_xyyyyz_yyz_1[j] * fl1_fza);

                    tg_xxxyyyyz_yzz_0[j] = pb_x * tg_xxyyyyz_yzz_0[j] + fr * tg_xxyyyyz_yzz_1[j] + fl1_fx * (tg_xyyyyz_yzz_0[j] - tg_xyyyyz_yzz_1[j] * fl1_fza);

                    tg_xxxyyyyz_zzz_0[j] = pb_x * tg_xxyyyyz_zzz_0[j] + fr * tg_xxyyyyz_zzz_1[j] + fl1_fx * (tg_xyyyyz_zzz_0[j] - tg_xyyyyz_zzz_1[j] * fl1_fza);

                    tg_xxxyyyzz_xxx_0[j] = pb_x * tg_xxyyyzz_xxx_0[j] + fr * tg_xxyyyzz_xxx_1[j] + fl1_fx * (tg_xyyyzz_xxx_0[j] - tg_xyyyzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyyzz_xx_1[j];

                    tg_xxxyyyzz_xxy_0[j] = pb_x * tg_xxyyyzz_xxy_0[j] + fr * tg_xxyyyzz_xxy_1[j] + fl1_fx * (tg_xyyyzz_xxy_0[j] - tg_xyyyzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyzz_xy_1[j];

                    tg_xxxyyyzz_xxz_0[j] = pb_x * tg_xxyyyzz_xxz_0[j] + fr * tg_xxyyyzz_xxz_1[j] + fl1_fx * (tg_xyyyzz_xxz_0[j] - tg_xyyyzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyyzz_xz_1[j];

                    tg_xxxyyyzz_xyy_0[j] = pb_x * tg_xxyyyzz_xyy_0[j] + fr * tg_xxyyyzz_xyy_1[j] + fl1_fx * (tg_xyyyzz_xyy_0[j] - tg_xyyyzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_yy_1[j];

                    tg_xxxyyyzz_xyz_0[j] = pb_x * tg_xxyyyzz_xyz_0[j] + fr * tg_xxyyyzz_xyz_1[j] + fl1_fx * (tg_xyyyzz_xyz_0[j] - tg_xyyyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_yz_1[j];

                    tg_xxxyyyzz_xzz_0[j] = pb_x * tg_xxyyyzz_xzz_0[j] + fr * tg_xxyyyzz_xzz_1[j] + fl1_fx * (tg_xyyyzz_xzz_0[j] - tg_xyyyzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyyzz_zz_1[j];

                    tg_xxxyyyzz_yyy_0[j] = pb_x * tg_xxyyyzz_yyy_0[j] + fr * tg_xxyyyzz_yyy_1[j] + fl1_fx * (tg_xyyyzz_yyy_0[j] - tg_xyyyzz_yyy_1[j] * fl1_fza);

                    tg_xxxyyyzz_yyz_0[j] = pb_x * tg_xxyyyzz_yyz_0[j] + fr * tg_xxyyyzz_yyz_1[j] + fl1_fx * (tg_xyyyzz_yyz_0[j] - tg_xyyyzz_yyz_1[j] * fl1_fza);

                    tg_xxxyyyzz_yzz_0[j] = pb_x * tg_xxyyyzz_yzz_0[j] + fr * tg_xxyyyzz_yzz_1[j] + fl1_fx * (tg_xyyyzz_yzz_0[j] - tg_xyyyzz_yzz_1[j] * fl1_fza);

                    tg_xxxyyyzz_zzz_0[j] = pb_x * tg_xxyyyzz_zzz_0[j] + fr * tg_xxyyyzz_zzz_1[j] + fl1_fx * (tg_xyyyzz_zzz_0[j] - tg_xyyyzz_zzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSF_180_270(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xxyyzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 180); 

                auto tg_xxyyzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 181); 

                auto tg_xxyyzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 182); 

                auto tg_xxyyzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 183); 

                auto tg_xxyyzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 184); 

                auto tg_xxyyzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 185); 

                auto tg_xxyyzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 186); 

                auto tg_xxyyzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 187); 

                auto tg_xxyyzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 188); 

                auto tg_xxyyzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 189); 

                auto tg_xxyzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 190); 

                auto tg_xxyzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 191); 

                auto tg_xxyzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 192); 

                auto tg_xxyzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 193); 

                auto tg_xxyzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 194); 

                auto tg_xxyzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 195); 

                auto tg_xxyzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 196); 

                auto tg_xxyzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 197); 

                auto tg_xxyzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 198); 

                auto tg_xxyzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 199); 

                auto tg_xxzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 200); 

                auto tg_xxzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 201); 

                auto tg_xxzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 202); 

                auto tg_xxzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 203); 

                auto tg_xxzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 204); 

                auto tg_xxzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 205); 

                auto tg_xxzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 206); 

                auto tg_xxzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 207); 

                auto tg_xxzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 208); 

                auto tg_xxzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 209); 

                auto tg_xyyyyyy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 210); 

                auto tg_xyyyyyy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 211); 

                auto tg_xyyyyyy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 212); 

                auto tg_xyyyyyy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 213); 

                auto tg_xyyyyyy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 214); 

                auto tg_xyyyyyy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 215); 

                auto tg_xyyyyyy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 216); 

                auto tg_xyyyyyy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 217); 

                auto tg_xyyyyyy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 218); 

                auto tg_xyyyyyy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 219); 

                auto tg_xyyyyyz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 220); 

                auto tg_xyyyyyz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 221); 

                auto tg_xyyyyyz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 222); 

                auto tg_xyyyyyz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 223); 

                auto tg_xyyyyyz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 224); 

                auto tg_xyyyyyz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 225); 

                auto tg_xyyyyyz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 226); 

                auto tg_xyyyyyz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 227); 

                auto tg_xyyyyyz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 228); 

                auto tg_xyyyyyz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 229); 

                auto tg_xyyyyzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 230); 

                auto tg_xyyyyzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 231); 

                auto tg_xyyyyzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 232); 

                auto tg_xyyyyzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 233); 

                auto tg_xyyyyzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 234); 

                auto tg_xyyyyzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 235); 

                auto tg_xyyyyzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 236); 

                auto tg_xyyyyzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 237); 

                auto tg_xyyyyzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 238); 

                auto tg_xyyyyzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 239); 

                auto tg_xyyyzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 240); 

                auto tg_xyyyzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 241); 

                auto tg_xyyyzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 242); 

                auto tg_xyyyzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 243); 

                auto tg_xyyyzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 244); 

                auto tg_xyyyzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 245); 

                auto tg_xyyyzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 246); 

                auto tg_xyyyzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 247); 

                auto tg_xyyyzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 248); 

                auto tg_xyyyzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 249); 

                auto tg_xyyzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 250); 

                auto tg_xyyzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 251); 

                auto tg_xyyzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 252); 

                auto tg_xyyzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 253); 

                auto tg_xyyzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 254); 

                auto tg_xyyzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 255); 

                auto tg_xyyzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 256); 

                auto tg_xyyzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 257); 

                auto tg_xyyzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 258); 

                auto tg_xyyzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 259); 

                auto tg_xyzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 260); 

                auto tg_xyzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 261); 

                auto tg_xyzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 262); 

                auto tg_xyzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 263); 

                auto tg_xyzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 264); 

                auto tg_xyzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 265); 

                auto tg_xyzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 266); 

                auto tg_xyzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 267); 

                auto tg_xyzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 268); 

                auto tg_xyzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 269); 

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

                auto tg_xxyyzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 108); 

                auto tg_xxyyzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 109); 

                auto tg_xxyyzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 110); 

                auto tg_xxyyzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 111); 

                auto tg_xxyyzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 112); 

                auto tg_xxyyzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 113); 

                auto tg_xxyzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 114); 

                auto tg_xxyzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 115); 

                auto tg_xxyzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 116); 

                auto tg_xxyzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 117); 

                auto tg_xxyzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 118); 

                auto tg_xxyzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 119); 

                auto tg_xxzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 120); 

                auto tg_xxzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 121); 

                auto tg_xxzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 122); 

                auto tg_xxzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 123); 

                auto tg_xxzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 124); 

                auto tg_xxzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 125); 

                auto tg_xyyyyyy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 126); 

                auto tg_xyyyyyy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 127); 

                auto tg_xyyyyyy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 128); 

                auto tg_xyyyyyy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 129); 

                auto tg_xyyyyyy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 130); 

                auto tg_xyyyyyy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 131); 

                auto tg_xyyyyyz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 132); 

                auto tg_xyyyyyz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 133); 

                auto tg_xyyyyyz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 134); 

                auto tg_xyyyyyz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 135); 

                auto tg_xyyyyyz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 136); 

                auto tg_xyyyyyz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 137); 

                auto tg_xyyyyzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 138); 

                auto tg_xyyyyzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 139); 

                auto tg_xyyyyzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 140); 

                auto tg_xyyyyzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 141); 

                auto tg_xyyyyzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 142); 

                auto tg_xyyyyzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 143); 

                auto tg_xyyyzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 144); 

                auto tg_xyyyzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 145); 

                auto tg_xyyyzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 146); 

                auto tg_xyyyzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 147); 

                auto tg_xyyyzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 148); 

                auto tg_xyyyzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 149); 

                auto tg_xyyzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 150); 

                auto tg_xyyzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 151); 

                auto tg_xyyzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 152); 

                auto tg_xyyzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 153); 

                auto tg_xyyzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 154); 

                auto tg_xyyzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 155); 

                auto tg_xyzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 156); 

                auto tg_xyzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 157); 

                auto tg_xyzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 158); 

                auto tg_xyzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 159); 

                auto tg_xyzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 160); 

                auto tg_xyzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 161); 

                // set up pointers to integrals

                auto tg_xxxyyzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 180); 

                auto tg_xxxyyzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 181); 

                auto tg_xxxyyzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 182); 

                auto tg_xxxyyzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 183); 

                auto tg_xxxyyzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 184); 

                auto tg_xxxyyzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 185); 

                auto tg_xxxyyzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 186); 

                auto tg_xxxyyzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 187); 

                auto tg_xxxyyzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 188); 

                auto tg_xxxyyzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 189); 

                auto tg_xxxyzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 190); 

                auto tg_xxxyzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 191); 

                auto tg_xxxyzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 192); 

                auto tg_xxxyzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 193); 

                auto tg_xxxyzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 194); 

                auto tg_xxxyzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 195); 

                auto tg_xxxyzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 196); 

                auto tg_xxxyzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 197); 

                auto tg_xxxyzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 198); 

                auto tg_xxxyzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 199); 

                auto tg_xxxzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 200); 

                auto tg_xxxzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 201); 

                auto tg_xxxzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 202); 

                auto tg_xxxzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 203); 

                auto tg_xxxzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 204); 

                auto tg_xxxzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 205); 

                auto tg_xxxzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 206); 

                auto tg_xxxzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 207); 

                auto tg_xxxzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 208); 

                auto tg_xxxzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 209); 

                auto tg_xxyyyyyy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 210); 

                auto tg_xxyyyyyy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 211); 

                auto tg_xxyyyyyy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 212); 

                auto tg_xxyyyyyy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 213); 

                auto tg_xxyyyyyy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 214); 

                auto tg_xxyyyyyy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 215); 

                auto tg_xxyyyyyy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 216); 

                auto tg_xxyyyyyy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 217); 

                auto tg_xxyyyyyy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 218); 

                auto tg_xxyyyyyy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 219); 

                auto tg_xxyyyyyz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 220); 

                auto tg_xxyyyyyz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 221); 

                auto tg_xxyyyyyz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 222); 

                auto tg_xxyyyyyz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 223); 

                auto tg_xxyyyyyz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 224); 

                auto tg_xxyyyyyz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 225); 

                auto tg_xxyyyyyz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 226); 

                auto tg_xxyyyyyz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 227); 

                auto tg_xxyyyyyz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 228); 

                auto tg_xxyyyyyz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 229); 

                auto tg_xxyyyyzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 230); 

                auto tg_xxyyyyzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 231); 

                auto tg_xxyyyyzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 232); 

                auto tg_xxyyyyzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 233); 

                auto tg_xxyyyyzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 234); 

                auto tg_xxyyyyzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 235); 

                auto tg_xxyyyyzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 236); 

                auto tg_xxyyyyzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 237); 

                auto tg_xxyyyyzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 238); 

                auto tg_xxyyyyzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 239); 

                auto tg_xxyyyzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 240); 

                auto tg_xxyyyzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 241); 

                auto tg_xxyyyzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 242); 

                auto tg_xxyyyzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 243); 

                auto tg_xxyyyzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 244); 

                auto tg_xxyyyzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 245); 

                auto tg_xxyyyzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 246); 

                auto tg_xxyyyzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 247); 

                auto tg_xxyyyzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 248); 

                auto tg_xxyyyzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 249); 

                auto tg_xxyyzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 250); 

                auto tg_xxyyzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 251); 

                auto tg_xxyyzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 252); 

                auto tg_xxyyzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 253); 

                auto tg_xxyyzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 254); 

                auto tg_xxyyzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 255); 

                auto tg_xxyyzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 256); 

                auto tg_xxyyzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 257); 

                auto tg_xxyyzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 258); 

                auto tg_xxyyzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 259); 

                auto tg_xxyzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 260); 

                auto tg_xxyzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 261); 

                auto tg_xxyzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 262); 

                auto tg_xxyzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 263); 

                auto tg_xxyzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 264); 

                auto tg_xxyzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 265); 

                auto tg_xxyzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 266); 

                auto tg_xxyzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 267); 

                auto tg_xxyzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 268); 

                auto tg_xxyzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 269); 

                // Batch of Integrals (180,270)

                #pragma omp simd aligned(fxn, fza, tg_xxxyyzzz_xxx_0, tg_xxxyyzzz_xxy_0, tg_xxxyyzzz_xxz_0, \
                                         tg_xxxyyzzz_xyy_0, tg_xxxyyzzz_xyz_0, tg_xxxyyzzz_xzz_0, tg_xxxyyzzz_yyy_0, \
                                         tg_xxxyyzzz_yyz_0, tg_xxxyyzzz_yzz_0, tg_xxxyyzzz_zzz_0, tg_xxxyzzzz_xxx_0, \
                                         tg_xxxyzzzz_xxy_0, tg_xxxyzzzz_xxz_0, tg_xxxyzzzz_xyy_0, tg_xxxyzzzz_xyz_0, \
                                         tg_xxxyzzzz_xzz_0, tg_xxxyzzzz_yyy_0, tg_xxxyzzzz_yyz_0, tg_xxxyzzzz_yzz_0, \
                                         tg_xxxyzzzz_zzz_0, tg_xxxzzzzz_xxx_0, tg_xxxzzzzz_xxy_0, tg_xxxzzzzz_xxz_0, \
                                         tg_xxxzzzzz_xyy_0, tg_xxxzzzzz_xyz_0, tg_xxxzzzzz_xzz_0, tg_xxxzzzzz_yyy_0, \
                                         tg_xxxzzzzz_yyz_0, tg_xxxzzzzz_yzz_0, tg_xxxzzzzz_zzz_0, tg_xxyyyyyy_xxx_0, \
                                         tg_xxyyyyyy_xxy_0, tg_xxyyyyyy_xxz_0, tg_xxyyyyyy_xyy_0, tg_xxyyyyyy_xyz_0, \
                                         tg_xxyyyyyy_xzz_0, tg_xxyyyyyy_yyy_0, tg_xxyyyyyy_yyz_0, tg_xxyyyyyy_yzz_0, \
                                         tg_xxyyyyyy_zzz_0, tg_xxyyyyyz_xxx_0, tg_xxyyyyyz_xxy_0, tg_xxyyyyyz_xxz_0, \
                                         tg_xxyyyyyz_xyy_0, tg_xxyyyyyz_xyz_0, tg_xxyyyyyz_xzz_0, tg_xxyyyyyz_yyy_0, \
                                         tg_xxyyyyyz_yyz_0, tg_xxyyyyyz_yzz_0, tg_xxyyyyyz_zzz_0, tg_xxyyyyzz_xxx_0, \
                                         tg_xxyyyyzz_xxy_0, tg_xxyyyyzz_xxz_0, tg_xxyyyyzz_xyy_0, tg_xxyyyyzz_xyz_0, \
                                         tg_xxyyyyzz_xzz_0, tg_xxyyyyzz_yyy_0, tg_xxyyyyzz_yyz_0, tg_xxyyyyzz_yzz_0, \
                                         tg_xxyyyyzz_zzz_0, tg_xxyyyzzz_xxx_0, tg_xxyyyzzz_xxy_0, tg_xxyyyzzz_xxz_0, \
                                         tg_xxyyyzzz_xyy_0, tg_xxyyyzzz_xyz_0, tg_xxyyyzzz_xzz_0, tg_xxyyyzzz_yyy_0, \
                                         tg_xxyyyzzz_yyz_0, tg_xxyyyzzz_yzz_0, tg_xxyyyzzz_zzz_0, tg_xxyyzzz_xx_1, \
                                         tg_xxyyzzz_xxx_0, tg_xxyyzzz_xxx_1, tg_xxyyzzz_xxy_0, tg_xxyyzzz_xxy_1, \
                                         tg_xxyyzzz_xxz_0, tg_xxyyzzz_xxz_1, tg_xxyyzzz_xy_1, tg_xxyyzzz_xyy_0, \
                                         tg_xxyyzzz_xyy_1, tg_xxyyzzz_xyz_0, tg_xxyyzzz_xyz_1, tg_xxyyzzz_xz_1, \
                                         tg_xxyyzzz_xzz_0, tg_xxyyzzz_xzz_1, tg_xxyyzzz_yy_1, tg_xxyyzzz_yyy_0, \
                                         tg_xxyyzzz_yyy_1, tg_xxyyzzz_yyz_0, tg_xxyyzzz_yyz_1, tg_xxyyzzz_yz_1, \
                                         tg_xxyyzzz_yzz_0, tg_xxyyzzz_yzz_1, tg_xxyyzzz_zz_1, tg_xxyyzzz_zzz_0, \
                                         tg_xxyyzzz_zzz_1, tg_xxyyzzzz_xxx_0, tg_xxyyzzzz_xxy_0, tg_xxyyzzzz_xxz_0, \
                                         tg_xxyyzzzz_xyy_0, tg_xxyyzzzz_xyz_0, tg_xxyyzzzz_xzz_0, tg_xxyyzzzz_yyy_0, \
                                         tg_xxyyzzzz_yyz_0, tg_xxyyzzzz_yzz_0, tg_xxyyzzzz_zzz_0, tg_xxyzzzz_xx_1, \
                                         tg_xxyzzzz_xxx_0, tg_xxyzzzz_xxx_1, tg_xxyzzzz_xxy_0, tg_xxyzzzz_xxy_1, \
                                         tg_xxyzzzz_xxz_0, tg_xxyzzzz_xxz_1, tg_xxyzzzz_xy_1, tg_xxyzzzz_xyy_0, \
                                         tg_xxyzzzz_xyy_1, tg_xxyzzzz_xyz_0, tg_xxyzzzz_xyz_1, tg_xxyzzzz_xz_1, \
                                         tg_xxyzzzz_xzz_0, tg_xxyzzzz_xzz_1, tg_xxyzzzz_yy_1, tg_xxyzzzz_yyy_0, \
                                         tg_xxyzzzz_yyy_1, tg_xxyzzzz_yyz_0, tg_xxyzzzz_yyz_1, tg_xxyzzzz_yz_1, \
                                         tg_xxyzzzz_yzz_0, tg_xxyzzzz_yzz_1, tg_xxyzzzz_zz_1, tg_xxyzzzz_zzz_0, \
                                         tg_xxyzzzz_zzz_1, tg_xxyzzzzz_xxx_0, tg_xxyzzzzz_xxy_0, tg_xxyzzzzz_xxz_0, \
                                         tg_xxyzzzzz_xyy_0, tg_xxyzzzzz_xyz_0, tg_xxyzzzzz_xzz_0, tg_xxyzzzzz_yyy_0, \
                                         tg_xxyzzzzz_yyz_0, tg_xxyzzzzz_yzz_0, tg_xxyzzzzz_zzz_0, tg_xxzzzzz_xx_1, \
                                         tg_xxzzzzz_xxx_0, tg_xxzzzzz_xxx_1, tg_xxzzzzz_xxy_0, tg_xxzzzzz_xxy_1, \
                                         tg_xxzzzzz_xxz_0, tg_xxzzzzz_xxz_1, tg_xxzzzzz_xy_1, tg_xxzzzzz_xyy_0, \
                                         tg_xxzzzzz_xyy_1, tg_xxzzzzz_xyz_0, tg_xxzzzzz_xyz_1, tg_xxzzzzz_xz_1, \
                                         tg_xxzzzzz_xzz_0, tg_xxzzzzz_xzz_1, tg_xxzzzzz_yy_1, tg_xxzzzzz_yyy_0, \
                                         tg_xxzzzzz_yyy_1, tg_xxzzzzz_yyz_0, tg_xxzzzzz_yyz_1, tg_xxzzzzz_yz_1, \
                                         tg_xxzzzzz_yzz_0, tg_xxzzzzz_yzz_1, tg_xxzzzzz_zz_1, tg_xxzzzzz_zzz_0, \
                                         tg_xxzzzzz_zzz_1, tg_xyyyyyy_xx_1, tg_xyyyyyy_xxx_0, tg_xyyyyyy_xxx_1, \
                                         tg_xyyyyyy_xxy_0, tg_xyyyyyy_xxy_1, tg_xyyyyyy_xxz_0, tg_xyyyyyy_xxz_1, \
                                         tg_xyyyyyy_xy_1, tg_xyyyyyy_xyy_0, tg_xyyyyyy_xyy_1, tg_xyyyyyy_xyz_0, \
                                         tg_xyyyyyy_xyz_1, tg_xyyyyyy_xz_1, tg_xyyyyyy_xzz_0, tg_xyyyyyy_xzz_1, \
                                         tg_xyyyyyy_yy_1, tg_xyyyyyy_yyy_0, tg_xyyyyyy_yyy_1, tg_xyyyyyy_yyz_0, \
                                         tg_xyyyyyy_yyz_1, tg_xyyyyyy_yz_1, tg_xyyyyyy_yzz_0, tg_xyyyyyy_yzz_1, \
                                         tg_xyyyyyy_zz_1, tg_xyyyyyy_zzz_0, tg_xyyyyyy_zzz_1, tg_xyyyyyz_xx_1, \
                                         tg_xyyyyyz_xxx_0, tg_xyyyyyz_xxx_1, tg_xyyyyyz_xxy_0, tg_xyyyyyz_xxy_1, \
                                         tg_xyyyyyz_xxz_0, tg_xyyyyyz_xxz_1, tg_xyyyyyz_xy_1, tg_xyyyyyz_xyy_0, \
                                         tg_xyyyyyz_xyy_1, tg_xyyyyyz_xyz_0, tg_xyyyyyz_xyz_1, tg_xyyyyyz_xz_1, \
                                         tg_xyyyyyz_xzz_0, tg_xyyyyyz_xzz_1, tg_xyyyyyz_yy_1, tg_xyyyyyz_yyy_0, \
                                         tg_xyyyyyz_yyy_1, tg_xyyyyyz_yyz_0, tg_xyyyyyz_yyz_1, tg_xyyyyyz_yz_1, \
                                         tg_xyyyyyz_yzz_0, tg_xyyyyyz_yzz_1, tg_xyyyyyz_zz_1, tg_xyyyyyz_zzz_0, \
                                         tg_xyyyyyz_zzz_1, tg_xyyyyzz_xx_1, tg_xyyyyzz_xxx_0, tg_xyyyyzz_xxx_1, \
                                         tg_xyyyyzz_xxy_0, tg_xyyyyzz_xxy_1, tg_xyyyyzz_xxz_0, tg_xyyyyzz_xxz_1, \
                                         tg_xyyyyzz_xy_1, tg_xyyyyzz_xyy_0, tg_xyyyyzz_xyy_1, tg_xyyyyzz_xyz_0, \
                                         tg_xyyyyzz_xyz_1, tg_xyyyyzz_xz_1, tg_xyyyyzz_xzz_0, tg_xyyyyzz_xzz_1, \
                                         tg_xyyyyzz_yy_1, tg_xyyyyzz_yyy_0, tg_xyyyyzz_yyy_1, tg_xyyyyzz_yyz_0, \
                                         tg_xyyyyzz_yyz_1, tg_xyyyyzz_yz_1, tg_xyyyyzz_yzz_0, tg_xyyyyzz_yzz_1, \
                                         tg_xyyyyzz_zz_1, tg_xyyyyzz_zzz_0, tg_xyyyyzz_zzz_1, tg_xyyyzzz_xx_1, \
                                         tg_xyyyzzz_xxx_0, tg_xyyyzzz_xxx_1, tg_xyyyzzz_xxy_0, tg_xyyyzzz_xxy_1, \
                                         tg_xyyyzzz_xxz_0, tg_xyyyzzz_xxz_1, tg_xyyyzzz_xy_1, tg_xyyyzzz_xyy_0, \
                                         tg_xyyyzzz_xyy_1, tg_xyyyzzz_xyz_0, tg_xyyyzzz_xyz_1, tg_xyyyzzz_xz_1, \
                                         tg_xyyyzzz_xzz_0, tg_xyyyzzz_xzz_1, tg_xyyyzzz_yy_1, tg_xyyyzzz_yyy_0, \
                                         tg_xyyyzzz_yyy_1, tg_xyyyzzz_yyz_0, tg_xyyyzzz_yyz_1, tg_xyyyzzz_yz_1, \
                                         tg_xyyyzzz_yzz_0, tg_xyyyzzz_yzz_1, tg_xyyyzzz_zz_1, tg_xyyyzzz_zzz_0, \
                                         tg_xyyyzzz_zzz_1, tg_xyyzzz_xxx_0, tg_xyyzzz_xxx_1, tg_xyyzzz_xxy_0, tg_xyyzzz_xxy_1, \
                                         tg_xyyzzz_xxz_0, tg_xyyzzz_xxz_1, tg_xyyzzz_xyy_0, tg_xyyzzz_xyy_1, tg_xyyzzz_xyz_0, \
                                         tg_xyyzzz_xyz_1, tg_xyyzzz_xzz_0, tg_xyyzzz_xzz_1, tg_xyyzzz_yyy_0, tg_xyyzzz_yyy_1, \
                                         tg_xyyzzz_yyz_0, tg_xyyzzz_yyz_1, tg_xyyzzz_yzz_0, tg_xyyzzz_yzz_1, tg_xyyzzz_zzz_0, \
                                         tg_xyyzzz_zzz_1, tg_xyyzzzz_xx_1, tg_xyyzzzz_xxx_0, tg_xyyzzzz_xxx_1, \
                                         tg_xyyzzzz_xxy_0, tg_xyyzzzz_xxy_1, tg_xyyzzzz_xxz_0, tg_xyyzzzz_xxz_1, \
                                         tg_xyyzzzz_xy_1, tg_xyyzzzz_xyy_0, tg_xyyzzzz_xyy_1, tg_xyyzzzz_xyz_0, \
                                         tg_xyyzzzz_xyz_1, tg_xyyzzzz_xz_1, tg_xyyzzzz_xzz_0, tg_xyyzzzz_xzz_1, \
                                         tg_xyyzzzz_yy_1, tg_xyyzzzz_yyy_0, tg_xyyzzzz_yyy_1, tg_xyyzzzz_yyz_0, \
                                         tg_xyyzzzz_yyz_1, tg_xyyzzzz_yz_1, tg_xyyzzzz_yzz_0, tg_xyyzzzz_yzz_1, \
                                         tg_xyyzzzz_zz_1, tg_xyyzzzz_zzz_0, tg_xyyzzzz_zzz_1, tg_xyzzzz_xxx_0, \
                                         tg_xyzzzz_xxx_1, tg_xyzzzz_xxy_0, tg_xyzzzz_xxy_1, tg_xyzzzz_xxz_0, tg_xyzzzz_xxz_1, \
                                         tg_xyzzzz_xyy_0, tg_xyzzzz_xyy_1, tg_xyzzzz_xyz_0, tg_xyzzzz_xyz_1, tg_xyzzzz_xzz_0, \
                                         tg_xyzzzz_xzz_1, tg_xyzzzz_yyy_0, tg_xyzzzz_yyy_1, tg_xyzzzz_yyz_0, tg_xyzzzz_yyz_1, \
                                         tg_xyzzzz_yzz_0, tg_xyzzzz_yzz_1, tg_xyzzzz_zzz_0, tg_xyzzzz_zzz_1, tg_xyzzzzz_xx_1, \
                                         tg_xyzzzzz_xxx_0, tg_xyzzzzz_xxx_1, tg_xyzzzzz_xxy_0, tg_xyzzzzz_xxy_1, \
                                         tg_xyzzzzz_xxz_0, tg_xyzzzzz_xxz_1, tg_xyzzzzz_xy_1, tg_xyzzzzz_xyy_0, \
                                         tg_xyzzzzz_xyy_1, tg_xyzzzzz_xyz_0, tg_xyzzzzz_xyz_1, tg_xyzzzzz_xz_1, \
                                         tg_xyzzzzz_xzz_0, tg_xyzzzzz_xzz_1, tg_xyzzzzz_yy_1, tg_xyzzzzz_yyy_0, \
                                         tg_xyzzzzz_yyy_1, tg_xyzzzzz_yyz_0, tg_xyzzzzz_yyz_1, tg_xyzzzzz_yz_1, \
                                         tg_xyzzzzz_yzz_0, tg_xyzzzzz_yzz_1, tg_xyzzzzz_zz_1, tg_xyzzzzz_zzz_0, \
                                         tg_xyzzzzz_zzz_1, tg_xzzzzz_xxx_0, tg_xzzzzz_xxx_1, tg_xzzzzz_xxy_0, tg_xzzzzz_xxy_1, \
                                         tg_xzzzzz_xxz_0, tg_xzzzzz_xxz_1, tg_xzzzzz_xyy_0, tg_xzzzzz_xyy_1, tg_xzzzzz_xyz_0, \
                                         tg_xzzzzz_xyz_1, tg_xzzzzz_xzz_0, tg_xzzzzz_xzz_1, tg_xzzzzz_yyy_0, tg_xzzzzz_yyy_1, \
                                         tg_xzzzzz_yyz_0, tg_xzzzzz_yyz_1, tg_xzzzzz_yzz_0, tg_xzzzzz_yzz_1, tg_xzzzzz_zzz_0, \
                                         tg_xzzzzz_zzz_1, tg_yyyyyy_xxx_0, tg_yyyyyy_xxx_1, tg_yyyyyy_xxy_0, tg_yyyyyy_xxy_1, \
                                         tg_yyyyyy_xxz_0, tg_yyyyyy_xxz_1, tg_yyyyyy_xyy_0, tg_yyyyyy_xyy_1, tg_yyyyyy_xyz_0, \
                                         tg_yyyyyy_xyz_1, tg_yyyyyy_xzz_0, tg_yyyyyy_xzz_1, tg_yyyyyy_yyy_0, tg_yyyyyy_yyy_1, \
                                         tg_yyyyyy_yyz_0, tg_yyyyyy_yyz_1, tg_yyyyyy_yzz_0, tg_yyyyyy_yzz_1, tg_yyyyyy_zzz_0, \
                                         tg_yyyyyy_zzz_1, tg_yyyyyz_xxx_0, tg_yyyyyz_xxx_1, tg_yyyyyz_xxy_0, tg_yyyyyz_xxy_1, \
                                         tg_yyyyyz_xxz_0, tg_yyyyyz_xxz_1, tg_yyyyyz_xyy_0, tg_yyyyyz_xyy_1, tg_yyyyyz_xyz_0, \
                                         tg_yyyyyz_xyz_1, tg_yyyyyz_xzz_0, tg_yyyyyz_xzz_1, tg_yyyyyz_yyy_0, tg_yyyyyz_yyy_1, \
                                         tg_yyyyyz_yyz_0, tg_yyyyyz_yyz_1, tg_yyyyyz_yzz_0, tg_yyyyyz_yzz_1, tg_yyyyyz_zzz_0, \
                                         tg_yyyyyz_zzz_1, tg_yyyyzz_xxx_0, tg_yyyyzz_xxx_1, tg_yyyyzz_xxy_0, tg_yyyyzz_xxy_1, \
                                         tg_yyyyzz_xxz_0, tg_yyyyzz_xxz_1, tg_yyyyzz_xyy_0, tg_yyyyzz_xyy_1, tg_yyyyzz_xyz_0, \
                                         tg_yyyyzz_xyz_1, tg_yyyyzz_xzz_0, tg_yyyyzz_xzz_1, tg_yyyyzz_yyy_0, tg_yyyyzz_yyy_1, \
                                         tg_yyyyzz_yyz_0, tg_yyyyzz_yyz_1, tg_yyyyzz_yzz_0, tg_yyyyzz_yzz_1, tg_yyyyzz_zzz_0, \
                                         tg_yyyyzz_zzz_1, tg_yyyzzz_xxx_0, tg_yyyzzz_xxx_1, tg_yyyzzz_xxy_0, tg_yyyzzz_xxy_1, \
                                         tg_yyyzzz_xxz_0, tg_yyyzzz_xxz_1, tg_yyyzzz_xyy_0, tg_yyyzzz_xyy_1, tg_yyyzzz_xyz_0, \
                                         tg_yyyzzz_xyz_1, tg_yyyzzz_xzz_0, tg_yyyzzz_xzz_1, tg_yyyzzz_yyy_0, tg_yyyzzz_yyy_1, \
                                         tg_yyyzzz_yyz_0, tg_yyyzzz_yyz_1, tg_yyyzzz_yzz_0, tg_yyyzzz_yzz_1, tg_yyyzzz_zzz_0, \
                                         tg_yyyzzz_zzz_1, tg_yyzzzz_xxx_0, tg_yyzzzz_xxx_1, tg_yyzzzz_xxy_0, tg_yyzzzz_xxy_1, \
                                         tg_yyzzzz_xxz_0, tg_yyzzzz_xxz_1, tg_yyzzzz_xyy_0, tg_yyzzzz_xyy_1, tg_yyzzzz_xyz_0, \
                                         tg_yyzzzz_xyz_1, tg_yyzzzz_xzz_0, tg_yyzzzz_xzz_1, tg_yyzzzz_yyy_0, tg_yyzzzz_yyy_1, \
                                         tg_yyzzzz_yyz_0, tg_yyzzzz_yyz_1, tg_yyzzzz_yzz_0, tg_yyzzzz_yzz_1, tg_yyzzzz_zzz_0, \
                                         tg_yyzzzz_zzz_1, tg_yzzzzz_xxx_0, tg_yzzzzz_xxx_1, tg_yzzzzz_xxy_0, tg_yzzzzz_xxy_1, \
                                         tg_yzzzzz_xxz_0, tg_yzzzzz_xxz_1, tg_yzzzzz_xyy_0, tg_yzzzzz_xyy_1, tg_yzzzzz_xyz_0, \
                                         tg_yzzzzz_xyz_1, tg_yzzzzz_xzz_0, tg_yzzzzz_xzz_1, tg_yzzzzz_yyy_0, tg_yzzzzz_yyy_1, \
                                         tg_yzzzzz_yyz_0, tg_yzzzzz_yyz_1, tg_yzzzzz_yzz_0, tg_yzzzzz_yzz_1, tg_yzzzzz_zzz_0, \
                                         tg_yzzzzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxyyzzz_xxx_0[j] = pb_x * tg_xxyyzzz_xxx_0[j] + fr * tg_xxyyzzz_xxx_1[j] + fl1_fx * (tg_xyyzzz_xxx_0[j] - tg_xyyzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyyzzz_xx_1[j];

                    tg_xxxyyzzz_xxy_0[j] = pb_x * tg_xxyyzzz_xxy_0[j] + fr * tg_xxyyzzz_xxy_1[j] + fl1_fx * (tg_xyyzzz_xxy_0[j] - tg_xyyzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzzz_xy_1[j];

                    tg_xxxyyzzz_xxz_0[j] = pb_x * tg_xxyyzzz_xxz_0[j] + fr * tg_xxyyzzz_xxz_1[j] + fl1_fx * (tg_xyyzzz_xxz_0[j] - tg_xyyzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyyzzz_xz_1[j];

                    tg_xxxyyzzz_xyy_0[j] = pb_x * tg_xxyyzzz_xyy_0[j] + fr * tg_xxyyzzz_xyy_1[j] + fl1_fx * (tg_xyyzzz_xyy_0[j] - tg_xyyzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_yy_1[j];

                    tg_xxxyyzzz_xyz_0[j] = pb_x * tg_xxyyzzz_xyz_0[j] + fr * tg_xxyyzzz_xyz_1[j] + fl1_fx * (tg_xyyzzz_xyz_0[j] - tg_xyyzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_yz_1[j];

                    tg_xxxyyzzz_xzz_0[j] = pb_x * tg_xxyyzzz_xzz_0[j] + fr * tg_xxyyzzz_xzz_1[j] + fl1_fx * (tg_xyyzzz_xzz_0[j] - tg_xyyzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyyzzz_zz_1[j];

                    tg_xxxyyzzz_yyy_0[j] = pb_x * tg_xxyyzzz_yyy_0[j] + fr * tg_xxyyzzz_yyy_1[j] + fl1_fx * (tg_xyyzzz_yyy_0[j] - tg_xyyzzz_yyy_1[j] * fl1_fza);

                    tg_xxxyyzzz_yyz_0[j] = pb_x * tg_xxyyzzz_yyz_0[j] + fr * tg_xxyyzzz_yyz_1[j] + fl1_fx * (tg_xyyzzz_yyz_0[j] - tg_xyyzzz_yyz_1[j] * fl1_fza);

                    tg_xxxyyzzz_yzz_0[j] = pb_x * tg_xxyyzzz_yzz_0[j] + fr * tg_xxyyzzz_yzz_1[j] + fl1_fx * (tg_xyyzzz_yzz_0[j] - tg_xyyzzz_yzz_1[j] * fl1_fza);

                    tg_xxxyyzzz_zzz_0[j] = pb_x * tg_xxyyzzz_zzz_0[j] + fr * tg_xxyyzzz_zzz_1[j] + fl1_fx * (tg_xyyzzz_zzz_0[j] - tg_xyyzzz_zzz_1[j] * fl1_fza);

                    tg_xxxyzzzz_xxx_0[j] = pb_x * tg_xxyzzzz_xxx_0[j] + fr * tg_xxyzzzz_xxx_1[j] + fl1_fx * (tg_xyzzzz_xxx_0[j] - tg_xyzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxyzzzz_xx_1[j];

                    tg_xxxyzzzz_xxy_0[j] = pb_x * tg_xxyzzzz_xxy_0[j] + fr * tg_xxyzzzz_xxy_1[j] + fl1_fx * (tg_xyzzzz_xxy_0[j] - tg_xyzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzzz_xy_1[j];

                    tg_xxxyzzzz_xxz_0[j] = pb_x * tg_xxyzzzz_xxz_0[j] + fr * tg_xxyzzzz_xxz_1[j] + fl1_fx * (tg_xyzzzz_xxz_0[j] - tg_xyzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxyzzzz_xz_1[j];

                    tg_xxxyzzzz_xyy_0[j] = pb_x * tg_xxyzzzz_xyy_0[j] + fr * tg_xxyzzzz_xyy_1[j] + fl1_fx * (tg_xyzzzz_xyy_0[j] - tg_xyzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_yy_1[j];

                    tg_xxxyzzzz_xyz_0[j] = pb_x * tg_xxyzzzz_xyz_0[j] + fr * tg_xxyzzzz_xyz_1[j] + fl1_fx * (tg_xyzzzz_xyz_0[j] - tg_xyzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_yz_1[j];

                    tg_xxxyzzzz_xzz_0[j] = pb_x * tg_xxyzzzz_xzz_0[j] + fr * tg_xxyzzzz_xzz_1[j] + fl1_fx * (tg_xyzzzz_xzz_0[j] - tg_xyzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxyzzzz_zz_1[j];

                    tg_xxxyzzzz_yyy_0[j] = pb_x * tg_xxyzzzz_yyy_0[j] + fr * tg_xxyzzzz_yyy_1[j] + fl1_fx * (tg_xyzzzz_yyy_0[j] - tg_xyzzzz_yyy_1[j] * fl1_fza);

                    tg_xxxyzzzz_yyz_0[j] = pb_x * tg_xxyzzzz_yyz_0[j] + fr * tg_xxyzzzz_yyz_1[j] + fl1_fx * (tg_xyzzzz_yyz_0[j] - tg_xyzzzz_yyz_1[j] * fl1_fza);

                    tg_xxxyzzzz_yzz_0[j] = pb_x * tg_xxyzzzz_yzz_0[j] + fr * tg_xxyzzzz_yzz_1[j] + fl1_fx * (tg_xyzzzz_yzz_0[j] - tg_xyzzzz_yzz_1[j] * fl1_fza);

                    tg_xxxyzzzz_zzz_0[j] = pb_x * tg_xxyzzzz_zzz_0[j] + fr * tg_xxyzzzz_zzz_1[j] + fl1_fx * (tg_xyzzzz_zzz_0[j] - tg_xyzzzz_zzz_1[j] * fl1_fza);

                    tg_xxxzzzzz_xxx_0[j] = pb_x * tg_xxzzzzz_xxx_0[j] + fr * tg_xxzzzzz_xxx_1[j] + fl1_fx * (tg_xzzzzz_xxx_0[j] - tg_xzzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxzzzzz_xx_1[j];

                    tg_xxxzzzzz_xxy_0[j] = pb_x * tg_xxzzzzz_xxy_0[j] + fr * tg_xxzzzzz_xxy_1[j] + fl1_fx * (tg_xzzzzz_xxy_0[j] - tg_xzzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzzz_xy_1[j];

                    tg_xxxzzzzz_xxz_0[j] = pb_x * tg_xxzzzzz_xxz_0[j] + fr * tg_xxzzzzz_xxz_1[j] + fl1_fx * (tg_xzzzzz_xxz_0[j] - tg_xzzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxzzzzz_xz_1[j];

                    tg_xxxzzzzz_xyy_0[j] = pb_x * tg_xxzzzzz_xyy_0[j] + fr * tg_xxzzzzz_xyy_1[j] + fl1_fx * (tg_xzzzzz_xyy_0[j] - tg_xzzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_yy_1[j];

                    tg_xxxzzzzz_xyz_0[j] = pb_x * tg_xxzzzzz_xyz_0[j] + fr * tg_xxzzzzz_xyz_1[j] + fl1_fx * (tg_xzzzzz_xyz_0[j] - tg_xzzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_yz_1[j];

                    tg_xxxzzzzz_xzz_0[j] = pb_x * tg_xxzzzzz_xzz_0[j] + fr * tg_xxzzzzz_xzz_1[j] + fl1_fx * (tg_xzzzzz_xzz_0[j] - tg_xzzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxzzzzz_zz_1[j];

                    tg_xxxzzzzz_yyy_0[j] = pb_x * tg_xxzzzzz_yyy_0[j] + fr * tg_xxzzzzz_yyy_1[j] + fl1_fx * (tg_xzzzzz_yyy_0[j] - tg_xzzzzz_yyy_1[j] * fl1_fza);

                    tg_xxxzzzzz_yyz_0[j] = pb_x * tg_xxzzzzz_yyz_0[j] + fr * tg_xxzzzzz_yyz_1[j] + fl1_fx * (tg_xzzzzz_yyz_0[j] - tg_xzzzzz_yyz_1[j] * fl1_fza);

                    tg_xxxzzzzz_yzz_0[j] = pb_x * tg_xxzzzzz_yzz_0[j] + fr * tg_xxzzzzz_yzz_1[j] + fl1_fx * (tg_xzzzzz_yzz_0[j] - tg_xzzzzz_yzz_1[j] * fl1_fza);

                    tg_xxxzzzzz_zzz_0[j] = pb_x * tg_xxzzzzz_zzz_0[j] + fr * tg_xxzzzzz_zzz_1[j] + fl1_fx * (tg_xzzzzz_zzz_0[j] - tg_xzzzzz_zzz_1[j] * fl1_fza);

                    tg_xxyyyyyy_xxx_0[j] = pb_x * tg_xyyyyyy_xxx_0[j] + fr * tg_xyyyyyy_xxx_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxx_0[j] - tg_yyyyyy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyy_xx_1[j];

                    tg_xxyyyyyy_xxy_0[j] = pb_x * tg_xyyyyyy_xxy_0[j] + fr * tg_xyyyyyy_xxy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxy_0[j] - tg_yyyyyy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyy_xy_1[j];

                    tg_xxyyyyyy_xxz_0[j] = pb_x * tg_xyyyyyy_xxz_0[j] + fr * tg_xyyyyyy_xxz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xxz_0[j] - tg_yyyyyy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyy_xz_1[j];

                    tg_xxyyyyyy_xyy_0[j] = pb_x * tg_xyyyyyy_xyy_0[j] + fr * tg_xyyyyyy_xyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xyy_0[j] - tg_yyyyyy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_yy_1[j];

                    tg_xxyyyyyy_xyz_0[j] = pb_x * tg_xyyyyyy_xyz_0[j] + fr * tg_xyyyyyy_xyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xyz_0[j] - tg_yyyyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_yz_1[j];

                    tg_xxyyyyyy_xzz_0[j] = pb_x * tg_xyyyyyy_xzz_0[j] + fr * tg_xyyyyyy_xzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_xzz_0[j] - tg_yyyyyy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyy_zz_1[j];

                    tg_xxyyyyyy_yyy_0[j] = pb_x * tg_xyyyyyy_yyy_0[j] + fr * tg_xyyyyyy_yyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yyy_0[j] - tg_yyyyyy_yyy_1[j] * fl1_fza);

                    tg_xxyyyyyy_yyz_0[j] = pb_x * tg_xyyyyyy_yyz_0[j] + fr * tg_xyyyyyy_yyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yyz_0[j] - tg_yyyyyy_yyz_1[j] * fl1_fza);

                    tg_xxyyyyyy_yzz_0[j] = pb_x * tg_xyyyyyy_yzz_0[j] + fr * tg_xyyyyyy_yzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_yzz_0[j] - tg_yyyyyy_yzz_1[j] * fl1_fza);

                    tg_xxyyyyyy_zzz_0[j] = pb_x * tg_xyyyyyy_zzz_0[j] + fr * tg_xyyyyyy_zzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_zzz_0[j] - tg_yyyyyy_zzz_1[j] * fl1_fza);

                    tg_xxyyyyyz_xxx_0[j] = pb_x * tg_xyyyyyz_xxx_0[j] + fr * tg_xyyyyyz_xxx_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxx_0[j] - tg_yyyyyz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyyz_xx_1[j];

                    tg_xxyyyyyz_xxy_0[j] = pb_x * tg_xyyyyyz_xxy_0[j] + fr * tg_xyyyyyz_xxy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxy_0[j] - tg_yyyyyz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyz_xy_1[j];

                    tg_xxyyyyyz_xxz_0[j] = pb_x * tg_xyyyyyz_xxz_0[j] + fr * tg_xyyyyyz_xxz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xxz_0[j] - tg_yyyyyz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyyz_xz_1[j];

                    tg_xxyyyyyz_xyy_0[j] = pb_x * tg_xyyyyyz_xyy_0[j] + fr * tg_xyyyyyz_xyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xyy_0[j] - tg_yyyyyz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_yy_1[j];

                    tg_xxyyyyyz_xyz_0[j] = pb_x * tg_xyyyyyz_xyz_0[j] + fr * tg_xyyyyyz_xyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xyz_0[j] - tg_yyyyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_yz_1[j];

                    tg_xxyyyyyz_xzz_0[j] = pb_x * tg_xyyyyyz_xzz_0[j] + fr * tg_xyyyyyz_xzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_xzz_0[j] - tg_yyyyyz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyyz_zz_1[j];

                    tg_xxyyyyyz_yyy_0[j] = pb_x * tg_xyyyyyz_yyy_0[j] + fr * tg_xyyyyyz_yyy_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yyy_0[j] - tg_yyyyyz_yyy_1[j] * fl1_fza);

                    tg_xxyyyyyz_yyz_0[j] = pb_x * tg_xyyyyyz_yyz_0[j] + fr * tg_xyyyyyz_yyz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yyz_0[j] - tg_yyyyyz_yyz_1[j] * fl1_fza);

                    tg_xxyyyyyz_yzz_0[j] = pb_x * tg_xyyyyyz_yzz_0[j] + fr * tg_xyyyyyz_yzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_yzz_0[j] - tg_yyyyyz_yzz_1[j] * fl1_fza);

                    tg_xxyyyyyz_zzz_0[j] = pb_x * tg_xyyyyyz_zzz_0[j] + fr * tg_xyyyyyz_zzz_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_zzz_0[j] - tg_yyyyyz_zzz_1[j] * fl1_fza);

                    tg_xxyyyyzz_xxx_0[j] = pb_x * tg_xyyyyzz_xxx_0[j] + fr * tg_xyyyyzz_xxx_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxx_0[j] - tg_yyyyzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyyzz_xx_1[j];

                    tg_xxyyyyzz_xxy_0[j] = pb_x * tg_xyyyyzz_xxy_0[j] + fr * tg_xyyyyzz_xxy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxy_0[j] - tg_yyyyzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyzz_xy_1[j];

                    tg_xxyyyyzz_xxz_0[j] = pb_x * tg_xyyyyzz_xxz_0[j] + fr * tg_xyyyyzz_xxz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xxz_0[j] - tg_yyyyzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyyzz_xz_1[j];

                    tg_xxyyyyzz_xyy_0[j] = pb_x * tg_xyyyyzz_xyy_0[j] + fr * tg_xyyyyzz_xyy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xyy_0[j] - tg_yyyyzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_yy_1[j];

                    tg_xxyyyyzz_xyz_0[j] = pb_x * tg_xyyyyzz_xyz_0[j] + fr * tg_xyyyyzz_xyz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xyz_0[j] - tg_yyyyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_yz_1[j];

                    tg_xxyyyyzz_xzz_0[j] = pb_x * tg_xyyyyzz_xzz_0[j] + fr * tg_xyyyyzz_xzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_xzz_0[j] - tg_yyyyzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyyzz_zz_1[j];

                    tg_xxyyyyzz_yyy_0[j] = pb_x * tg_xyyyyzz_yyy_0[j] + fr * tg_xyyyyzz_yyy_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yyy_0[j] - tg_yyyyzz_yyy_1[j] * fl1_fza);

                    tg_xxyyyyzz_yyz_0[j] = pb_x * tg_xyyyyzz_yyz_0[j] + fr * tg_xyyyyzz_yyz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yyz_0[j] - tg_yyyyzz_yyz_1[j] * fl1_fza);

                    tg_xxyyyyzz_yzz_0[j] = pb_x * tg_xyyyyzz_yzz_0[j] + fr * tg_xyyyyzz_yzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_yzz_0[j] - tg_yyyyzz_yzz_1[j] * fl1_fza);

                    tg_xxyyyyzz_zzz_0[j] = pb_x * tg_xyyyyzz_zzz_0[j] + fr * tg_xyyyyzz_zzz_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_zzz_0[j] - tg_yyyyzz_zzz_1[j] * fl1_fza);

                    tg_xxyyyzzz_xxx_0[j] = pb_x * tg_xyyyzzz_xxx_0[j] + fr * tg_xyyyzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxx_0[j] - tg_yyyzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyyzzz_xx_1[j];

                    tg_xxyyyzzz_xxy_0[j] = pb_x * tg_xyyyzzz_xxy_0[j] + fr * tg_xyyyzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxy_0[j] - tg_yyyzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzzz_xy_1[j];

                    tg_xxyyyzzz_xxz_0[j] = pb_x * tg_xyyyzzz_xxz_0[j] + fr * tg_xyyyzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xxz_0[j] - tg_yyyzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyyzzz_xz_1[j];

                    tg_xxyyyzzz_xyy_0[j] = pb_x * tg_xyyyzzz_xyy_0[j] + fr * tg_xyyyzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xyy_0[j] - tg_yyyzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_yy_1[j];

                    tg_xxyyyzzz_xyz_0[j] = pb_x * tg_xyyyzzz_xyz_0[j] + fr * tg_xyyyzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xyz_0[j] - tg_yyyzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_yz_1[j];

                    tg_xxyyyzzz_xzz_0[j] = pb_x * tg_xyyyzzz_xzz_0[j] + fr * tg_xyyyzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_xzz_0[j] - tg_yyyzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyyzzz_zz_1[j];

                    tg_xxyyyzzz_yyy_0[j] = pb_x * tg_xyyyzzz_yyy_0[j] + fr * tg_xyyyzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yyy_0[j] - tg_yyyzzz_yyy_1[j] * fl1_fza);

                    tg_xxyyyzzz_yyz_0[j] = pb_x * tg_xyyyzzz_yyz_0[j] + fr * tg_xyyyzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yyz_0[j] - tg_yyyzzz_yyz_1[j] * fl1_fza);

                    tg_xxyyyzzz_yzz_0[j] = pb_x * tg_xyyyzzz_yzz_0[j] + fr * tg_xyyyzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_yzz_0[j] - tg_yyyzzz_yzz_1[j] * fl1_fza);

                    tg_xxyyyzzz_zzz_0[j] = pb_x * tg_xyyyzzz_zzz_0[j] + fr * tg_xyyyzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_zzz_0[j] - tg_yyyzzz_zzz_1[j] * fl1_fza);

                    tg_xxyyzzzz_xxx_0[j] = pb_x * tg_xyyzzzz_xxx_0[j] + fr * tg_xyyzzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxx_0[j] - tg_yyzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyyzzzz_xx_1[j];

                    tg_xxyyzzzz_xxy_0[j] = pb_x * tg_xyyzzzz_xxy_0[j] + fr * tg_xyyzzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxy_0[j] - tg_yyzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzzz_xy_1[j];

                    tg_xxyyzzzz_xxz_0[j] = pb_x * tg_xyyzzzz_xxz_0[j] + fr * tg_xyyzzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xxz_0[j] - tg_yyzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyyzzzz_xz_1[j];

                    tg_xxyyzzzz_xyy_0[j] = pb_x * tg_xyyzzzz_xyy_0[j] + fr * tg_xyyzzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xyy_0[j] - tg_yyzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_yy_1[j];

                    tg_xxyyzzzz_xyz_0[j] = pb_x * tg_xyyzzzz_xyz_0[j] + fr * tg_xyyzzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xyz_0[j] - tg_yyzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_yz_1[j];

                    tg_xxyyzzzz_xzz_0[j] = pb_x * tg_xyyzzzz_xzz_0[j] + fr * tg_xyyzzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_xzz_0[j] - tg_yyzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyyzzzz_zz_1[j];

                    tg_xxyyzzzz_yyy_0[j] = pb_x * tg_xyyzzzz_yyy_0[j] + fr * tg_xyyzzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yyy_0[j] - tg_yyzzzz_yyy_1[j] * fl1_fza);

                    tg_xxyyzzzz_yyz_0[j] = pb_x * tg_xyyzzzz_yyz_0[j] + fr * tg_xyyzzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yyz_0[j] - tg_yyzzzz_yyz_1[j] * fl1_fza);

                    tg_xxyyzzzz_yzz_0[j] = pb_x * tg_xyyzzzz_yzz_0[j] + fr * tg_xyyzzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_yzz_0[j] - tg_yyzzzz_yzz_1[j] * fl1_fza);

                    tg_xxyyzzzz_zzz_0[j] = pb_x * tg_xyyzzzz_zzz_0[j] + fr * tg_xyyzzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_zzz_0[j] - tg_yyzzzz_zzz_1[j] * fl1_fza);

                    tg_xxyzzzzz_xxx_0[j] = pb_x * tg_xyzzzzz_xxx_0[j] + fr * tg_xyzzzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxx_0[j] - tg_yzzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyzzzzz_xx_1[j];

                    tg_xxyzzzzz_xxy_0[j] = pb_x * tg_xyzzzzz_xxy_0[j] + fr * tg_xyzzzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxy_0[j] - tg_yzzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzzz_xy_1[j];

                    tg_xxyzzzzz_xxz_0[j] = pb_x * tg_xyzzzzz_xxz_0[j] + fr * tg_xyzzzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xxz_0[j] - tg_yzzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyzzzzz_xz_1[j];

                    tg_xxyzzzzz_xyy_0[j] = pb_x * tg_xyzzzzz_xyy_0[j] + fr * tg_xyzzzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xyy_0[j] - tg_yzzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_yy_1[j];

                    tg_xxyzzzzz_xyz_0[j] = pb_x * tg_xyzzzzz_xyz_0[j] + fr * tg_xyzzzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xyz_0[j] - tg_yzzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_yz_1[j];

                    tg_xxyzzzzz_xzz_0[j] = pb_x * tg_xyzzzzz_xzz_0[j] + fr * tg_xyzzzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_xzz_0[j] - tg_yzzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyzzzzz_zz_1[j];

                    tg_xxyzzzzz_yyy_0[j] = pb_x * tg_xyzzzzz_yyy_0[j] + fr * tg_xyzzzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yyy_0[j] - tg_yzzzzz_yyy_1[j] * fl1_fza);

                    tg_xxyzzzzz_yyz_0[j] = pb_x * tg_xyzzzzz_yyz_0[j] + fr * tg_xyzzzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yyz_0[j] - tg_yzzzzz_yyz_1[j] * fl1_fza);

                    tg_xxyzzzzz_yzz_0[j] = pb_x * tg_xyzzzzz_yzz_0[j] + fr * tg_xyzzzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_yzz_0[j] - tg_yzzzzz_yzz_1[j] * fl1_fza);

                    tg_xxyzzzzz_zzz_0[j] = pb_x * tg_xyzzzzz_zzz_0[j] + fr * tg_xyzzzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_zzz_0[j] - tg_yzzzzz_zzz_1[j] * fl1_fza);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSF_270_360(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xzzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 270); 

                auto tg_xzzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 271); 

                auto tg_xzzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 272); 

                auto tg_xzzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 273); 

                auto tg_xzzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 274); 

                auto tg_xzzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 275); 

                auto tg_xzzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 276); 

                auto tg_xzzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 277); 

                auto tg_xzzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 278); 

                auto tg_xzzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 279); 

                auto tg_yyyyyyy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 280); 

                auto tg_yyyyyyy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 281); 

                auto tg_yyyyyyy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 282); 

                auto tg_yyyyyyy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 283); 

                auto tg_yyyyyyy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 284); 

                auto tg_yyyyyyy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 285); 

                auto tg_yyyyyyy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 286); 

                auto tg_yyyyyyy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 287); 

                auto tg_yyyyyyy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 288); 

                auto tg_yyyyyyy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 289); 

                auto tg_yyyyyyz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 290); 

                auto tg_yyyyyyz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 291); 

                auto tg_yyyyyyz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 292); 

                auto tg_yyyyyyz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 293); 

                auto tg_yyyyyyz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 294); 

                auto tg_yyyyyyz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 295); 

                auto tg_yyyyyyz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 296); 

                auto tg_yyyyyyz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 297); 

                auto tg_yyyyyyz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 298); 

                auto tg_yyyyyyz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 299); 

                auto tg_yyyyyzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 300); 

                auto tg_yyyyyzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 301); 

                auto tg_yyyyyzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 302); 

                auto tg_yyyyyzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 303); 

                auto tg_yyyyyzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 304); 

                auto tg_yyyyyzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 305); 

                auto tg_yyyyyzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 306); 

                auto tg_yyyyyzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 307); 

                auto tg_yyyyyzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 308); 

                auto tg_yyyyyzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 309); 

                auto tg_yyyyzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 310); 

                auto tg_yyyyzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 311); 

                auto tg_yyyyzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 312); 

                auto tg_yyyyzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 313); 

                auto tg_yyyyzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 314); 

                auto tg_yyyyzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 315); 

                auto tg_yyyyzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 316); 

                auto tg_yyyyzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 317); 

                auto tg_yyyyzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 318); 

                auto tg_yyyyzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 319); 

                auto tg_yyyzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 320); 

                auto tg_yyyzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 321); 

                auto tg_yyyzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 322); 

                auto tg_yyyzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 323); 

                auto tg_yyyzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 324); 

                auto tg_yyyzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 325); 

                auto tg_yyyzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 326); 

                auto tg_yyyzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 327); 

                auto tg_yyyzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 328); 

                auto tg_yyyzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 329); 

                auto tg_yyzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 330); 

                auto tg_yyzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 331); 

                auto tg_yyzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 332); 

                auto tg_yyzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 333); 

                auto tg_yyzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 334); 

                auto tg_yyzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 335); 

                auto tg_yyzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 336); 

                auto tg_yyzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 337); 

                auto tg_yyzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 338); 

                auto tg_yyzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 339); 

                auto tg_yzzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 340); 

                auto tg_yzzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 341); 

                auto tg_yzzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 342); 

                auto tg_yzzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 343); 

                auto tg_yzzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 344); 

                auto tg_yzzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 345); 

                auto tg_yzzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 346); 

                auto tg_yzzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 347); 

                auto tg_yzzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 348); 

                auto tg_yzzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 349); 

                auto tg_zzzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 350); 

                auto tg_zzzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 351); 

                auto tg_zzzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 352); 

                auto tg_zzzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 353); 

                auto tg_zzzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 354); 

                auto tg_zzzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 355); 

                auto tg_zzzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 356); 

                auto tg_zzzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 357); 

                auto tg_zzzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 358); 

                auto tg_zzzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 359); 

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

                auto tg_xzzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 162); 

                auto tg_xzzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 163); 

                auto tg_xzzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 164); 

                auto tg_xzzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 165); 

                auto tg_xzzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 166); 

                auto tg_xzzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 167); 

                auto tg_yyyyyyy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 168); 

                auto tg_yyyyyyy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 169); 

                auto tg_yyyyyyy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 170); 

                auto tg_yyyyyyy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 171); 

                auto tg_yyyyyyy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 172); 

                auto tg_yyyyyyy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 173); 

                auto tg_yyyyyyz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 174); 

                auto tg_yyyyyyz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 175); 

                auto tg_yyyyyyz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 176); 

                auto tg_yyyyyyz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 177); 

                auto tg_yyyyyyz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 178); 

                auto tg_yyyyyyz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 179); 

                auto tg_yyyyyzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 180); 

                auto tg_yyyyyzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 181); 

                auto tg_yyyyyzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 182); 

                auto tg_yyyyyzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 183); 

                auto tg_yyyyyzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 184); 

                auto tg_yyyyyzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 185); 

                auto tg_yyyyzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 186); 

                auto tg_yyyyzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 187); 

                auto tg_yyyyzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 188); 

                auto tg_yyyyzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 189); 

                auto tg_yyyyzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 190); 

                auto tg_yyyyzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 191); 

                auto tg_yyyzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 192); 

                auto tg_yyyzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 193); 

                auto tg_yyyzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 194); 

                auto tg_yyyzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 195); 

                auto tg_yyyzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 196); 

                auto tg_yyyzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 197); 

                auto tg_yyzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 198); 

                auto tg_yyzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 199); 

                auto tg_yyzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 200); 

                auto tg_yyzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 201); 

                auto tg_yyzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 202); 

                auto tg_yyzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 203); 

                auto tg_yzzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 204); 

                auto tg_yzzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 205); 

                auto tg_yzzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 206); 

                auto tg_yzzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 207); 

                auto tg_yzzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 208); 

                auto tg_yzzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 209); 

                auto tg_zzzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 210); 

                auto tg_zzzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 211); 

                auto tg_zzzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 212); 

                auto tg_zzzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 213); 

                auto tg_zzzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 214); 

                auto tg_zzzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 215); 

                // set up pointers to integrals

                auto tg_xxzzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 270); 

                auto tg_xxzzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 271); 

                auto tg_xxzzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 272); 

                auto tg_xxzzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 273); 

                auto tg_xxzzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 274); 

                auto tg_xxzzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 275); 

                auto tg_xxzzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 276); 

                auto tg_xxzzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 277); 

                auto tg_xxzzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 278); 

                auto tg_xxzzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 279); 

                auto tg_xyyyyyyy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 280); 

                auto tg_xyyyyyyy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 281); 

                auto tg_xyyyyyyy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 282); 

                auto tg_xyyyyyyy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 283); 

                auto tg_xyyyyyyy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 284); 

                auto tg_xyyyyyyy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 285); 

                auto tg_xyyyyyyy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 286); 

                auto tg_xyyyyyyy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 287); 

                auto tg_xyyyyyyy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 288); 

                auto tg_xyyyyyyy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 289); 

                auto tg_xyyyyyyz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 290); 

                auto tg_xyyyyyyz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 291); 

                auto tg_xyyyyyyz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 292); 

                auto tg_xyyyyyyz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 293); 

                auto tg_xyyyyyyz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 294); 

                auto tg_xyyyyyyz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 295); 

                auto tg_xyyyyyyz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 296); 

                auto tg_xyyyyyyz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 297); 

                auto tg_xyyyyyyz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 298); 

                auto tg_xyyyyyyz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 299); 

                auto tg_xyyyyyzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 300); 

                auto tg_xyyyyyzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 301); 

                auto tg_xyyyyyzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 302); 

                auto tg_xyyyyyzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 303); 

                auto tg_xyyyyyzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 304); 

                auto tg_xyyyyyzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 305); 

                auto tg_xyyyyyzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 306); 

                auto tg_xyyyyyzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 307); 

                auto tg_xyyyyyzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 308); 

                auto tg_xyyyyyzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 309); 

                auto tg_xyyyyzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 310); 

                auto tg_xyyyyzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 311); 

                auto tg_xyyyyzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 312); 

                auto tg_xyyyyzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 313); 

                auto tg_xyyyyzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 314); 

                auto tg_xyyyyzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 315); 

                auto tg_xyyyyzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 316); 

                auto tg_xyyyyzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 317); 

                auto tg_xyyyyzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 318); 

                auto tg_xyyyyzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 319); 

                auto tg_xyyyzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 320); 

                auto tg_xyyyzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 321); 

                auto tg_xyyyzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 322); 

                auto tg_xyyyzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 323); 

                auto tg_xyyyzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 324); 

                auto tg_xyyyzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 325); 

                auto tg_xyyyzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 326); 

                auto tg_xyyyzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 327); 

                auto tg_xyyyzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 328); 

                auto tg_xyyyzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 329); 

                auto tg_xyyzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 330); 

                auto tg_xyyzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 331); 

                auto tg_xyyzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 332); 

                auto tg_xyyzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 333); 

                auto tg_xyyzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 334); 

                auto tg_xyyzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 335); 

                auto tg_xyyzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 336); 

                auto tg_xyyzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 337); 

                auto tg_xyyzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 338); 

                auto tg_xyyzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 339); 

                auto tg_xyzzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 340); 

                auto tg_xyzzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 341); 

                auto tg_xyzzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 342); 

                auto tg_xyzzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 343); 

                auto tg_xyzzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 344); 

                auto tg_xyzzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 345); 

                auto tg_xyzzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 346); 

                auto tg_xyzzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 347); 

                auto tg_xyzzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 348); 

                auto tg_xyzzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 349); 

                auto tg_xzzzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 350); 

                auto tg_xzzzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 351); 

                auto tg_xzzzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 352); 

                auto tg_xzzzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 353); 

                auto tg_xzzzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 354); 

                auto tg_xzzzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 355); 

                auto tg_xzzzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 356); 

                auto tg_xzzzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 357); 

                auto tg_xzzzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 358); 

                auto tg_xzzzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 359); 

                // Batch of Integrals (270,360)

                #pragma omp simd aligned(fxn, fza, tg_xxzzzzzz_xxx_0, tg_xxzzzzzz_xxy_0, tg_xxzzzzzz_xxz_0, \
                                         tg_xxzzzzzz_xyy_0, tg_xxzzzzzz_xyz_0, tg_xxzzzzzz_xzz_0, tg_xxzzzzzz_yyy_0, \
                                         tg_xxzzzzzz_yyz_0, tg_xxzzzzzz_yzz_0, tg_xxzzzzzz_zzz_0, tg_xyyyyyyy_xxx_0, \
                                         tg_xyyyyyyy_xxy_0, tg_xyyyyyyy_xxz_0, tg_xyyyyyyy_xyy_0, tg_xyyyyyyy_xyz_0, \
                                         tg_xyyyyyyy_xzz_0, tg_xyyyyyyy_yyy_0, tg_xyyyyyyy_yyz_0, tg_xyyyyyyy_yzz_0, \
                                         tg_xyyyyyyy_zzz_0, tg_xyyyyyyz_xxx_0, tg_xyyyyyyz_xxy_0, tg_xyyyyyyz_xxz_0, \
                                         tg_xyyyyyyz_xyy_0, tg_xyyyyyyz_xyz_0, tg_xyyyyyyz_xzz_0, tg_xyyyyyyz_yyy_0, \
                                         tg_xyyyyyyz_yyz_0, tg_xyyyyyyz_yzz_0, tg_xyyyyyyz_zzz_0, tg_xyyyyyzz_xxx_0, \
                                         tg_xyyyyyzz_xxy_0, tg_xyyyyyzz_xxz_0, tg_xyyyyyzz_xyy_0, tg_xyyyyyzz_xyz_0, \
                                         tg_xyyyyyzz_xzz_0, tg_xyyyyyzz_yyy_0, tg_xyyyyyzz_yyz_0, tg_xyyyyyzz_yzz_0, \
                                         tg_xyyyyyzz_zzz_0, tg_xyyyyzzz_xxx_0, tg_xyyyyzzz_xxy_0, tg_xyyyyzzz_xxz_0, \
                                         tg_xyyyyzzz_xyy_0, tg_xyyyyzzz_xyz_0, tg_xyyyyzzz_xzz_0, tg_xyyyyzzz_yyy_0, \
                                         tg_xyyyyzzz_yyz_0, tg_xyyyyzzz_yzz_0, tg_xyyyyzzz_zzz_0, tg_xyyyzzzz_xxx_0, \
                                         tg_xyyyzzzz_xxy_0, tg_xyyyzzzz_xxz_0, tg_xyyyzzzz_xyy_0, tg_xyyyzzzz_xyz_0, \
                                         tg_xyyyzzzz_xzz_0, tg_xyyyzzzz_yyy_0, tg_xyyyzzzz_yyz_0, tg_xyyyzzzz_yzz_0, \
                                         tg_xyyyzzzz_zzz_0, tg_xyyzzzzz_xxx_0, tg_xyyzzzzz_xxy_0, tg_xyyzzzzz_xxz_0, \
                                         tg_xyyzzzzz_xyy_0, tg_xyyzzzzz_xyz_0, tg_xyyzzzzz_xzz_0, tg_xyyzzzzz_yyy_0, \
                                         tg_xyyzzzzz_yyz_0, tg_xyyzzzzz_yzz_0, tg_xyyzzzzz_zzz_0, tg_xyzzzzzz_xxx_0, \
                                         tg_xyzzzzzz_xxy_0, tg_xyzzzzzz_xxz_0, tg_xyzzzzzz_xyy_0, tg_xyzzzzzz_xyz_0, \
                                         tg_xyzzzzzz_xzz_0, tg_xyzzzzzz_yyy_0, tg_xyzzzzzz_yyz_0, tg_xyzzzzzz_yzz_0, \
                                         tg_xyzzzzzz_zzz_0, tg_xzzzzzz_xx_1, tg_xzzzzzz_xxx_0, tg_xzzzzzz_xxx_1, \
                                         tg_xzzzzzz_xxy_0, tg_xzzzzzz_xxy_1, tg_xzzzzzz_xxz_0, tg_xzzzzzz_xxz_1, \
                                         tg_xzzzzzz_xy_1, tg_xzzzzzz_xyy_0, tg_xzzzzzz_xyy_1, tg_xzzzzzz_xyz_0, \
                                         tg_xzzzzzz_xyz_1, tg_xzzzzzz_xz_1, tg_xzzzzzz_xzz_0, tg_xzzzzzz_xzz_1, \
                                         tg_xzzzzzz_yy_1, tg_xzzzzzz_yyy_0, tg_xzzzzzz_yyy_1, tg_xzzzzzz_yyz_0, \
                                         tg_xzzzzzz_yyz_1, tg_xzzzzzz_yz_1, tg_xzzzzzz_yzz_0, tg_xzzzzzz_yzz_1, \
                                         tg_xzzzzzz_zz_1, tg_xzzzzzz_zzz_0, tg_xzzzzzz_zzz_1, tg_xzzzzzzz_xxx_0, \
                                         tg_xzzzzzzz_xxy_0, tg_xzzzzzzz_xxz_0, tg_xzzzzzzz_xyy_0, tg_xzzzzzzz_xyz_0, \
                                         tg_xzzzzzzz_xzz_0, tg_xzzzzzzz_yyy_0, tg_xzzzzzzz_yyz_0, tg_xzzzzzzz_yzz_0, \
                                         tg_xzzzzzzz_zzz_0, tg_yyyyyyy_xx_1, tg_yyyyyyy_xxx_0, tg_yyyyyyy_xxx_1, \
                                         tg_yyyyyyy_xxy_0, tg_yyyyyyy_xxy_1, tg_yyyyyyy_xxz_0, tg_yyyyyyy_xxz_1, \
                                         tg_yyyyyyy_xy_1, tg_yyyyyyy_xyy_0, tg_yyyyyyy_xyy_1, tg_yyyyyyy_xyz_0, \
                                         tg_yyyyyyy_xyz_1, tg_yyyyyyy_xz_1, tg_yyyyyyy_xzz_0, tg_yyyyyyy_xzz_1, \
                                         tg_yyyyyyy_yy_1, tg_yyyyyyy_yyy_0, tg_yyyyyyy_yyy_1, tg_yyyyyyy_yyz_0, \
                                         tg_yyyyyyy_yyz_1, tg_yyyyyyy_yz_1, tg_yyyyyyy_yzz_0, tg_yyyyyyy_yzz_1, \
                                         tg_yyyyyyy_zz_1, tg_yyyyyyy_zzz_0, tg_yyyyyyy_zzz_1, tg_yyyyyyz_xx_1, \
                                         tg_yyyyyyz_xxx_0, tg_yyyyyyz_xxx_1, tg_yyyyyyz_xxy_0, tg_yyyyyyz_xxy_1, \
                                         tg_yyyyyyz_xxz_0, tg_yyyyyyz_xxz_1, tg_yyyyyyz_xy_1, tg_yyyyyyz_xyy_0, \
                                         tg_yyyyyyz_xyy_1, tg_yyyyyyz_xyz_0, tg_yyyyyyz_xyz_1, tg_yyyyyyz_xz_1, \
                                         tg_yyyyyyz_xzz_0, tg_yyyyyyz_xzz_1, tg_yyyyyyz_yy_1, tg_yyyyyyz_yyy_0, \
                                         tg_yyyyyyz_yyy_1, tg_yyyyyyz_yyz_0, tg_yyyyyyz_yyz_1, tg_yyyyyyz_yz_1, \
                                         tg_yyyyyyz_yzz_0, tg_yyyyyyz_yzz_1, tg_yyyyyyz_zz_1, tg_yyyyyyz_zzz_0, \
                                         tg_yyyyyyz_zzz_1, tg_yyyyyzz_xx_1, tg_yyyyyzz_xxx_0, tg_yyyyyzz_xxx_1, \
                                         tg_yyyyyzz_xxy_0, tg_yyyyyzz_xxy_1, tg_yyyyyzz_xxz_0, tg_yyyyyzz_xxz_1, \
                                         tg_yyyyyzz_xy_1, tg_yyyyyzz_xyy_0, tg_yyyyyzz_xyy_1, tg_yyyyyzz_xyz_0, \
                                         tg_yyyyyzz_xyz_1, tg_yyyyyzz_xz_1, tg_yyyyyzz_xzz_0, tg_yyyyyzz_xzz_1, \
                                         tg_yyyyyzz_yy_1, tg_yyyyyzz_yyy_0, tg_yyyyyzz_yyy_1, tg_yyyyyzz_yyz_0, \
                                         tg_yyyyyzz_yyz_1, tg_yyyyyzz_yz_1, tg_yyyyyzz_yzz_0, tg_yyyyyzz_yzz_1, \
                                         tg_yyyyyzz_zz_1, tg_yyyyyzz_zzz_0, tg_yyyyyzz_zzz_1, tg_yyyyzzz_xx_1, \
                                         tg_yyyyzzz_xxx_0, tg_yyyyzzz_xxx_1, tg_yyyyzzz_xxy_0, tg_yyyyzzz_xxy_1, \
                                         tg_yyyyzzz_xxz_0, tg_yyyyzzz_xxz_1, tg_yyyyzzz_xy_1, tg_yyyyzzz_xyy_0, \
                                         tg_yyyyzzz_xyy_1, tg_yyyyzzz_xyz_0, tg_yyyyzzz_xyz_1, tg_yyyyzzz_xz_1, \
                                         tg_yyyyzzz_xzz_0, tg_yyyyzzz_xzz_1, tg_yyyyzzz_yy_1, tg_yyyyzzz_yyy_0, \
                                         tg_yyyyzzz_yyy_1, tg_yyyyzzz_yyz_0, tg_yyyyzzz_yyz_1, tg_yyyyzzz_yz_1, \
                                         tg_yyyyzzz_yzz_0, tg_yyyyzzz_yzz_1, tg_yyyyzzz_zz_1, tg_yyyyzzz_zzz_0, \
                                         tg_yyyyzzz_zzz_1, tg_yyyzzzz_xx_1, tg_yyyzzzz_xxx_0, tg_yyyzzzz_xxx_1, \
                                         tg_yyyzzzz_xxy_0, tg_yyyzzzz_xxy_1, tg_yyyzzzz_xxz_0, tg_yyyzzzz_xxz_1, \
                                         tg_yyyzzzz_xy_1, tg_yyyzzzz_xyy_0, tg_yyyzzzz_xyy_1, tg_yyyzzzz_xyz_0, \
                                         tg_yyyzzzz_xyz_1, tg_yyyzzzz_xz_1, tg_yyyzzzz_xzz_0, tg_yyyzzzz_xzz_1, \
                                         tg_yyyzzzz_yy_1, tg_yyyzzzz_yyy_0, tg_yyyzzzz_yyy_1, tg_yyyzzzz_yyz_0, \
                                         tg_yyyzzzz_yyz_1, tg_yyyzzzz_yz_1, tg_yyyzzzz_yzz_0, tg_yyyzzzz_yzz_1, \
                                         tg_yyyzzzz_zz_1, tg_yyyzzzz_zzz_0, tg_yyyzzzz_zzz_1, tg_yyzzzzz_xx_1, \
                                         tg_yyzzzzz_xxx_0, tg_yyzzzzz_xxx_1, tg_yyzzzzz_xxy_0, tg_yyzzzzz_xxy_1, \
                                         tg_yyzzzzz_xxz_0, tg_yyzzzzz_xxz_1, tg_yyzzzzz_xy_1, tg_yyzzzzz_xyy_0, \
                                         tg_yyzzzzz_xyy_1, tg_yyzzzzz_xyz_0, tg_yyzzzzz_xyz_1, tg_yyzzzzz_xz_1, \
                                         tg_yyzzzzz_xzz_0, tg_yyzzzzz_xzz_1, tg_yyzzzzz_yy_1, tg_yyzzzzz_yyy_0, \
                                         tg_yyzzzzz_yyy_1, tg_yyzzzzz_yyz_0, tg_yyzzzzz_yyz_1, tg_yyzzzzz_yz_1, \
                                         tg_yyzzzzz_yzz_0, tg_yyzzzzz_yzz_1, tg_yyzzzzz_zz_1, tg_yyzzzzz_zzz_0, \
                                         tg_yyzzzzz_zzz_1, tg_yzzzzzz_xx_1, tg_yzzzzzz_xxx_0, tg_yzzzzzz_xxx_1, \
                                         tg_yzzzzzz_xxy_0, tg_yzzzzzz_xxy_1, tg_yzzzzzz_xxz_0, tg_yzzzzzz_xxz_1, \
                                         tg_yzzzzzz_xy_1, tg_yzzzzzz_xyy_0, tg_yzzzzzz_xyy_1, tg_yzzzzzz_xyz_0, \
                                         tg_yzzzzzz_xyz_1, tg_yzzzzzz_xz_1, tg_yzzzzzz_xzz_0, tg_yzzzzzz_xzz_1, \
                                         tg_yzzzzzz_yy_1, tg_yzzzzzz_yyy_0, tg_yzzzzzz_yyy_1, tg_yzzzzzz_yyz_0, \
                                         tg_yzzzzzz_yyz_1, tg_yzzzzzz_yz_1, tg_yzzzzzz_yzz_0, tg_yzzzzzz_yzz_1, \
                                         tg_yzzzzzz_zz_1, tg_yzzzzzz_zzz_0, tg_yzzzzzz_zzz_1, tg_zzzzzz_xxx_0, \
                                         tg_zzzzzz_xxx_1, tg_zzzzzz_xxy_0, tg_zzzzzz_xxy_1, tg_zzzzzz_xxz_0, tg_zzzzzz_xxz_1, \
                                         tg_zzzzzz_xyy_0, tg_zzzzzz_xyy_1, tg_zzzzzz_xyz_0, tg_zzzzzz_xyz_1, tg_zzzzzz_xzz_0, \
                                         tg_zzzzzz_xzz_1, tg_zzzzzz_yyy_0, tg_zzzzzz_yyy_1, tg_zzzzzz_yyz_0, tg_zzzzzz_yyz_1, \
                                         tg_zzzzzz_yzz_0, tg_zzzzzz_yzz_1, tg_zzzzzz_zzz_0, tg_zzzzzz_zzz_1, tg_zzzzzzz_xx_1, \
                                         tg_zzzzzzz_xxx_0, tg_zzzzzzz_xxx_1, tg_zzzzzzz_xxy_0, tg_zzzzzzz_xxy_1, \
                                         tg_zzzzzzz_xxz_0, tg_zzzzzzz_xxz_1, tg_zzzzzzz_xy_1, tg_zzzzzzz_xyy_0, \
                                         tg_zzzzzzz_xyy_1, tg_zzzzzzz_xyz_0, tg_zzzzzzz_xyz_1, tg_zzzzzzz_xz_1, \
                                         tg_zzzzzzz_xzz_0, tg_zzzzzzz_xzz_1, tg_zzzzzzz_yy_1, tg_zzzzzzz_yyy_0, \
                                         tg_zzzzzzz_yyy_1, tg_zzzzzzz_yyz_0, tg_zzzzzzz_yyz_1, tg_zzzzzzz_yz_1, \
                                         tg_zzzzzzz_yzz_0, tg_zzzzzzz_yzz_1, tg_zzzzzzz_zz_1, tg_zzzzzzz_zzz_0, \
                                         tg_zzzzzzz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxzzzzzz_xxx_0[j] = pb_x * tg_xzzzzzz_xxx_0[j] + fr * tg_xzzzzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxx_0[j] - tg_zzzzzz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzzzzzz_xx_1[j];

                    tg_xxzzzzzz_xxy_0[j] = pb_x * tg_xzzzzzz_xxy_0[j] + fr * tg_xzzzzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxy_0[j] - tg_zzzzzz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzzz_xy_1[j];

                    tg_xxzzzzzz_xxz_0[j] = pb_x * tg_xzzzzzz_xxz_0[j] + fr * tg_xzzzzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxz_0[j] - tg_zzzzzz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xzzzzzz_xz_1[j];

                    tg_xxzzzzzz_xyy_0[j] = pb_x * tg_xzzzzzz_xyy_0[j] + fr * tg_xzzzzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyy_0[j] - tg_zzzzzz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_yy_1[j];

                    tg_xxzzzzzz_xyz_0[j] = pb_x * tg_xzzzzzz_xyz_0[j] + fr * tg_xzzzzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyz_0[j] - tg_zzzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_yz_1[j];

                    tg_xxzzzzzz_xzz_0[j] = pb_x * tg_xzzzzzz_xzz_0[j] + fr * tg_xzzzzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xzz_0[j] - tg_zzzzzz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzzzzzz_zz_1[j];

                    tg_xxzzzzzz_yyy_0[j] = pb_x * tg_xzzzzzz_yyy_0[j] + fr * tg_xzzzzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyy_0[j] - tg_zzzzzz_yyy_1[j] * fl1_fza);

                    tg_xxzzzzzz_yyz_0[j] = pb_x * tg_xzzzzzz_yyz_0[j] + fr * tg_xzzzzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyz_0[j] - tg_zzzzzz_yyz_1[j] * fl1_fza);

                    tg_xxzzzzzz_yzz_0[j] = pb_x * tg_xzzzzzz_yzz_0[j] + fr * tg_xzzzzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yzz_0[j] - tg_zzzzzz_yzz_1[j] * fl1_fza);

                    tg_xxzzzzzz_zzz_0[j] = pb_x * tg_xzzzzzz_zzz_0[j] + fr * tg_xzzzzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_zzz_0[j] - tg_zzzzzz_zzz_1[j] * fl1_fza);

                    tg_xyyyyyyy_xxx_0[j] = pb_x * tg_yyyyyyy_xxx_0[j] + fr * tg_yyyyyyy_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyyyyy_xx_1[j];

                    tg_xyyyyyyy_xxy_0[j] = pb_x * tg_yyyyyyy_xxy_0[j] + fr * tg_yyyyyyy_xxy_1[j] + fl1_fxn * tg_yyyyyyy_xy_1[j];

                    tg_xyyyyyyy_xxz_0[j] = pb_x * tg_yyyyyyy_xxz_0[j] + fr * tg_yyyyyyy_xxz_1[j] + fl1_fxn * tg_yyyyyyy_xz_1[j];

                    tg_xyyyyyyy_xyy_0[j] = pb_x * tg_yyyyyyy_xyy_0[j] + fr * tg_yyyyyyy_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yy_1[j];

                    tg_xyyyyyyy_xyz_0[j] = pb_x * tg_yyyyyyy_xyz_0[j] + fr * tg_yyyyyyy_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_yz_1[j];

                    tg_xyyyyyyy_xzz_0[j] = pb_x * tg_yyyyyyy_xzz_0[j] + fr * tg_yyyyyyy_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_zz_1[j];

                    tg_xyyyyyyy_yyy_0[j] = pb_x * tg_yyyyyyy_yyy_0[j] + fr * tg_yyyyyyy_yyy_1[j];

                    tg_xyyyyyyy_yyz_0[j] = pb_x * tg_yyyyyyy_yyz_0[j] + fr * tg_yyyyyyy_yyz_1[j];

                    tg_xyyyyyyy_yzz_0[j] = pb_x * tg_yyyyyyy_yzz_0[j] + fr * tg_yyyyyyy_yzz_1[j];

                    tg_xyyyyyyy_zzz_0[j] = pb_x * tg_yyyyyyy_zzz_0[j] + fr * tg_yyyyyyy_zzz_1[j];

                    tg_xyyyyyyz_xxx_0[j] = pb_x * tg_yyyyyyz_xxx_0[j] + fr * tg_yyyyyyz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyyyyz_xx_1[j];

                    tg_xyyyyyyz_xxy_0[j] = pb_x * tg_yyyyyyz_xxy_0[j] + fr * tg_yyyyyyz_xxy_1[j] + fl1_fxn * tg_yyyyyyz_xy_1[j];

                    tg_xyyyyyyz_xxz_0[j] = pb_x * tg_yyyyyyz_xxz_0[j] + fr * tg_yyyyyyz_xxz_1[j] + fl1_fxn * tg_yyyyyyz_xz_1[j];

                    tg_xyyyyyyz_xyy_0[j] = pb_x * tg_yyyyyyz_xyy_0[j] + fr * tg_yyyyyyz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yy_1[j];

                    tg_xyyyyyyz_xyz_0[j] = pb_x * tg_yyyyyyz_xyz_0[j] + fr * tg_yyyyyyz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_yz_1[j];

                    tg_xyyyyyyz_xzz_0[j] = pb_x * tg_yyyyyyz_xzz_0[j] + fr * tg_yyyyyyz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_zz_1[j];

                    tg_xyyyyyyz_yyy_0[j] = pb_x * tg_yyyyyyz_yyy_0[j] + fr * tg_yyyyyyz_yyy_1[j];

                    tg_xyyyyyyz_yyz_0[j] = pb_x * tg_yyyyyyz_yyz_0[j] + fr * tg_yyyyyyz_yyz_1[j];

                    tg_xyyyyyyz_yzz_0[j] = pb_x * tg_yyyyyyz_yzz_0[j] + fr * tg_yyyyyyz_yzz_1[j];

                    tg_xyyyyyyz_zzz_0[j] = pb_x * tg_yyyyyyz_zzz_0[j] + fr * tg_yyyyyyz_zzz_1[j];

                    tg_xyyyyyzz_xxx_0[j] = pb_x * tg_yyyyyzz_xxx_0[j] + fr * tg_yyyyyzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyyyzz_xx_1[j];

                    tg_xyyyyyzz_xxy_0[j] = pb_x * tg_yyyyyzz_xxy_0[j] + fr * tg_yyyyyzz_xxy_1[j] + fl1_fxn * tg_yyyyyzz_xy_1[j];

                    tg_xyyyyyzz_xxz_0[j] = pb_x * tg_yyyyyzz_xxz_0[j] + fr * tg_yyyyyzz_xxz_1[j] + fl1_fxn * tg_yyyyyzz_xz_1[j];

                    tg_xyyyyyzz_xyy_0[j] = pb_x * tg_yyyyyzz_xyy_0[j] + fr * tg_yyyyyzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yy_1[j];

                    tg_xyyyyyzz_xyz_0[j] = pb_x * tg_yyyyyzz_xyz_0[j] + fr * tg_yyyyyzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_yz_1[j];

                    tg_xyyyyyzz_xzz_0[j] = pb_x * tg_yyyyyzz_xzz_0[j] + fr * tg_yyyyyzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_zz_1[j];

                    tg_xyyyyyzz_yyy_0[j] = pb_x * tg_yyyyyzz_yyy_0[j] + fr * tg_yyyyyzz_yyy_1[j];

                    tg_xyyyyyzz_yyz_0[j] = pb_x * tg_yyyyyzz_yyz_0[j] + fr * tg_yyyyyzz_yyz_1[j];

                    tg_xyyyyyzz_yzz_0[j] = pb_x * tg_yyyyyzz_yzz_0[j] + fr * tg_yyyyyzz_yzz_1[j];

                    tg_xyyyyyzz_zzz_0[j] = pb_x * tg_yyyyyzz_zzz_0[j] + fr * tg_yyyyyzz_zzz_1[j];

                    tg_xyyyyzzz_xxx_0[j] = pb_x * tg_yyyyzzz_xxx_0[j] + fr * tg_yyyyzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyyzzz_xx_1[j];

                    tg_xyyyyzzz_xxy_0[j] = pb_x * tg_yyyyzzz_xxy_0[j] + fr * tg_yyyyzzz_xxy_1[j] + fl1_fxn * tg_yyyyzzz_xy_1[j];

                    tg_xyyyyzzz_xxz_0[j] = pb_x * tg_yyyyzzz_xxz_0[j] + fr * tg_yyyyzzz_xxz_1[j] + fl1_fxn * tg_yyyyzzz_xz_1[j];

                    tg_xyyyyzzz_xyy_0[j] = pb_x * tg_yyyyzzz_xyy_0[j] + fr * tg_yyyyzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yy_1[j];

                    tg_xyyyyzzz_xyz_0[j] = pb_x * tg_yyyyzzz_xyz_0[j] + fr * tg_yyyyzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_yz_1[j];

                    tg_xyyyyzzz_xzz_0[j] = pb_x * tg_yyyyzzz_xzz_0[j] + fr * tg_yyyyzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_zz_1[j];

                    tg_xyyyyzzz_yyy_0[j] = pb_x * tg_yyyyzzz_yyy_0[j] + fr * tg_yyyyzzz_yyy_1[j];

                    tg_xyyyyzzz_yyz_0[j] = pb_x * tg_yyyyzzz_yyz_0[j] + fr * tg_yyyyzzz_yyz_1[j];

                    tg_xyyyyzzz_yzz_0[j] = pb_x * tg_yyyyzzz_yzz_0[j] + fr * tg_yyyyzzz_yzz_1[j];

                    tg_xyyyyzzz_zzz_0[j] = pb_x * tg_yyyyzzz_zzz_0[j] + fr * tg_yyyyzzz_zzz_1[j];

                    tg_xyyyzzzz_xxx_0[j] = pb_x * tg_yyyzzzz_xxx_0[j] + fr * tg_yyyzzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyyzzzz_xx_1[j];

                    tg_xyyyzzzz_xxy_0[j] = pb_x * tg_yyyzzzz_xxy_0[j] + fr * tg_yyyzzzz_xxy_1[j] + fl1_fxn * tg_yyyzzzz_xy_1[j];

                    tg_xyyyzzzz_xxz_0[j] = pb_x * tg_yyyzzzz_xxz_0[j] + fr * tg_yyyzzzz_xxz_1[j] + fl1_fxn * tg_yyyzzzz_xz_1[j];

                    tg_xyyyzzzz_xyy_0[j] = pb_x * tg_yyyzzzz_xyy_0[j] + fr * tg_yyyzzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yy_1[j];

                    tg_xyyyzzzz_xyz_0[j] = pb_x * tg_yyyzzzz_xyz_0[j] + fr * tg_yyyzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_yz_1[j];

                    tg_xyyyzzzz_xzz_0[j] = pb_x * tg_yyyzzzz_xzz_0[j] + fr * tg_yyyzzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_zz_1[j];

                    tg_xyyyzzzz_yyy_0[j] = pb_x * tg_yyyzzzz_yyy_0[j] + fr * tg_yyyzzzz_yyy_1[j];

                    tg_xyyyzzzz_yyz_0[j] = pb_x * tg_yyyzzzz_yyz_0[j] + fr * tg_yyyzzzz_yyz_1[j];

                    tg_xyyyzzzz_yzz_0[j] = pb_x * tg_yyyzzzz_yzz_0[j] + fr * tg_yyyzzzz_yzz_1[j];

                    tg_xyyyzzzz_zzz_0[j] = pb_x * tg_yyyzzzz_zzz_0[j] + fr * tg_yyyzzzz_zzz_1[j];

                    tg_xyyzzzzz_xxx_0[j] = pb_x * tg_yyzzzzz_xxx_0[j] + fr * tg_yyzzzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyzzzzz_xx_1[j];

                    tg_xyyzzzzz_xxy_0[j] = pb_x * tg_yyzzzzz_xxy_0[j] + fr * tg_yyzzzzz_xxy_1[j] + fl1_fxn * tg_yyzzzzz_xy_1[j];

                    tg_xyyzzzzz_xxz_0[j] = pb_x * tg_yyzzzzz_xxz_0[j] + fr * tg_yyzzzzz_xxz_1[j] + fl1_fxn * tg_yyzzzzz_xz_1[j];

                    tg_xyyzzzzz_xyy_0[j] = pb_x * tg_yyzzzzz_xyy_0[j] + fr * tg_yyzzzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yy_1[j];

                    tg_xyyzzzzz_xyz_0[j] = pb_x * tg_yyzzzzz_xyz_0[j] + fr * tg_yyzzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_yz_1[j];

                    tg_xyyzzzzz_xzz_0[j] = pb_x * tg_yyzzzzz_xzz_0[j] + fr * tg_yyzzzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_zz_1[j];

                    tg_xyyzzzzz_yyy_0[j] = pb_x * tg_yyzzzzz_yyy_0[j] + fr * tg_yyzzzzz_yyy_1[j];

                    tg_xyyzzzzz_yyz_0[j] = pb_x * tg_yyzzzzz_yyz_0[j] + fr * tg_yyzzzzz_yyz_1[j];

                    tg_xyyzzzzz_yzz_0[j] = pb_x * tg_yyzzzzz_yzz_0[j] + fr * tg_yyzzzzz_yzz_1[j];

                    tg_xyyzzzzz_zzz_0[j] = pb_x * tg_yyzzzzz_zzz_0[j] + fr * tg_yyzzzzz_zzz_1[j];

                    tg_xyzzzzzz_xxx_0[j] = pb_x * tg_yzzzzzz_xxx_0[j] + fr * tg_yzzzzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yzzzzzz_xx_1[j];

                    tg_xyzzzzzz_xxy_0[j] = pb_x * tg_yzzzzzz_xxy_0[j] + fr * tg_yzzzzzz_xxy_1[j] + fl1_fxn * tg_yzzzzzz_xy_1[j];

                    tg_xyzzzzzz_xxz_0[j] = pb_x * tg_yzzzzzz_xxz_0[j] + fr * tg_yzzzzzz_xxz_1[j] + fl1_fxn * tg_yzzzzzz_xz_1[j];

                    tg_xyzzzzzz_xyy_0[j] = pb_x * tg_yzzzzzz_xyy_0[j] + fr * tg_yzzzzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yy_1[j];

                    tg_xyzzzzzz_xyz_0[j] = pb_x * tg_yzzzzzz_xyz_0[j] + fr * tg_yzzzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_yz_1[j];

                    tg_xyzzzzzz_xzz_0[j] = pb_x * tg_yzzzzzz_xzz_0[j] + fr * tg_yzzzzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_zz_1[j];

                    tg_xyzzzzzz_yyy_0[j] = pb_x * tg_yzzzzzz_yyy_0[j] + fr * tg_yzzzzzz_yyy_1[j];

                    tg_xyzzzzzz_yyz_0[j] = pb_x * tg_yzzzzzz_yyz_0[j] + fr * tg_yzzzzzz_yyz_1[j];

                    tg_xyzzzzzz_yzz_0[j] = pb_x * tg_yzzzzzz_yzz_0[j] + fr * tg_yzzzzzz_yzz_1[j];

                    tg_xyzzzzzz_zzz_0[j] = pb_x * tg_yzzzzzz_zzz_0[j] + fr * tg_yzzzzzz_zzz_1[j];

                    tg_xzzzzzzz_xxx_0[j] = pb_x * tg_zzzzzzz_xxx_0[j] + fr * tg_zzzzzzz_xxx_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_xx_1[j];

                    tg_xzzzzzzz_xxy_0[j] = pb_x * tg_zzzzzzz_xxy_0[j] + fr * tg_zzzzzzz_xxy_1[j] + fl1_fxn * tg_zzzzzzz_xy_1[j];

                    tg_xzzzzzzz_xxz_0[j] = pb_x * tg_zzzzzzz_xxz_0[j] + fr * tg_zzzzzzz_xxz_1[j] + fl1_fxn * tg_zzzzzzz_xz_1[j];

                    tg_xzzzzzzz_xyy_0[j] = pb_x * tg_zzzzzzz_xyy_0[j] + fr * tg_zzzzzzz_xyy_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yy_1[j];

                    tg_xzzzzzzz_xyz_0[j] = pb_x * tg_zzzzzzz_xyz_0[j] + fr * tg_zzzzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_yz_1[j];

                    tg_xzzzzzzz_xzz_0[j] = pb_x * tg_zzzzzzz_xzz_0[j] + fr * tg_zzzzzzz_xzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_zz_1[j];

                    tg_xzzzzzzz_yyy_0[j] = pb_x * tg_zzzzzzz_yyy_0[j] + fr * tg_zzzzzzz_yyy_1[j];

                    tg_xzzzzzzz_yyz_0[j] = pb_x * tg_zzzzzzz_yyz_0[j] + fr * tg_zzzzzzz_yyz_1[j];

                    tg_xzzzzzzz_yzz_0[j] = pb_x * tg_zzzzzzz_yzz_0[j] + fr * tg_zzzzzzz_yzz_1[j];

                    tg_xzzzzzzz_zzz_0[j] = pb_x * tg_zzzzzzz_zzz_0[j] + fr * tg_zzzzzzz_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSF_360_450(      CMemBlock2D<double>* primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {8, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {8, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_7_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_6_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_6_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {6, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_7_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {7, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_yyyyyyy_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 280); 

                auto tg_yyyyyyy_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 281); 

                auto tg_yyyyyyy_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 282); 

                auto tg_yyyyyyy_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 283); 

                auto tg_yyyyyyy_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 284); 

                auto tg_yyyyyyy_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 285); 

                auto tg_yyyyyyy_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 286); 

                auto tg_yyyyyyy_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 287); 

                auto tg_yyyyyyy_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 288); 

                auto tg_yyyyyyy_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 289); 

                auto tg_yyyyyyz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 290); 

                auto tg_yyyyyyz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 291); 

                auto tg_yyyyyyz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 292); 

                auto tg_yyyyyyz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 293); 

                auto tg_yyyyyyz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 294); 

                auto tg_yyyyyyz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 295); 

                auto tg_yyyyyyz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 296); 

                auto tg_yyyyyyz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 297); 

                auto tg_yyyyyyz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 298); 

                auto tg_yyyyyyz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 299); 

                auto tg_yyyyyzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 300); 

                auto tg_yyyyyzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 301); 

                auto tg_yyyyyzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 302); 

                auto tg_yyyyyzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 303); 

                auto tg_yyyyyzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 304); 

                auto tg_yyyyyzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 305); 

                auto tg_yyyyyzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 306); 

                auto tg_yyyyyzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 307); 

                auto tg_yyyyyzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 308); 

                auto tg_yyyyyzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 309); 

                auto tg_yyyyzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 310); 

                auto tg_yyyyzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 311); 

                auto tg_yyyyzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 312); 

                auto tg_yyyyzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 313); 

                auto tg_yyyyzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 314); 

                auto tg_yyyyzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 315); 

                auto tg_yyyyzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 316); 

                auto tg_yyyyzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 317); 

                auto tg_yyyyzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 318); 

                auto tg_yyyyzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 319); 

                auto tg_yyyzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 320); 

                auto tg_yyyzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 321); 

                auto tg_yyyzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 322); 

                auto tg_yyyzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 323); 

                auto tg_yyyzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 324); 

                auto tg_yyyzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 325); 

                auto tg_yyyzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 326); 

                auto tg_yyyzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 327); 

                auto tg_yyyzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 328); 

                auto tg_yyyzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 329); 

                auto tg_yyzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 330); 

                auto tg_yyzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 331); 

                auto tg_yyzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 332); 

                auto tg_yyzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 333); 

                auto tg_yyzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 334); 

                auto tg_yyzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 335); 

                auto tg_yyzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 336); 

                auto tg_yyzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 337); 

                auto tg_yyzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 338); 

                auto tg_yyzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 339); 

                auto tg_yzzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 340); 

                auto tg_yzzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 341); 

                auto tg_yzzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 342); 

                auto tg_yzzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 343); 

                auto tg_yzzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 344); 

                auto tg_yzzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 345); 

                auto tg_yzzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 346); 

                auto tg_yzzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 347); 

                auto tg_yzzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 348); 

                auto tg_yzzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 349); 

                auto tg_zzzzzzz_xxx_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 350); 

                auto tg_zzzzzzz_xxy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 351); 

                auto tg_zzzzzzz_xxz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 352); 

                auto tg_zzzzzzz_xyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 353); 

                auto tg_zzzzzzz_xyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 354); 

                auto tg_zzzzzzz_xzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 355); 

                auto tg_zzzzzzz_yyy_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 356); 

                auto tg_zzzzzzz_yyz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 357); 

                auto tg_zzzzzzz_yzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 358); 

                auto tg_zzzzzzz_zzz_1 = primBuffer[pidx_g_7_3_m1].data(360 * idx + 359); 

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

                auto tg_yyyyyyy_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 168); 

                auto tg_yyyyyyy_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 169); 

                auto tg_yyyyyyy_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 170); 

                auto tg_yyyyyyy_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 171); 

                auto tg_yyyyyyy_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 172); 

                auto tg_yyyyyyy_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 173); 

                auto tg_yyyyyyz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 174); 

                auto tg_yyyyyyz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 175); 

                auto tg_yyyyyyz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 176); 

                auto tg_yyyyyyz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 177); 

                auto tg_yyyyyyz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 178); 

                auto tg_yyyyyyz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 179); 

                auto tg_yyyyyzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 180); 

                auto tg_yyyyyzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 181); 

                auto tg_yyyyyzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 182); 

                auto tg_yyyyyzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 183); 

                auto tg_yyyyyzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 184); 

                auto tg_yyyyyzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 185); 

                auto tg_yyyyzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 186); 

                auto tg_yyyyzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 187); 

                auto tg_yyyyzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 188); 

                auto tg_yyyyzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 189); 

                auto tg_yyyyzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 190); 

                auto tg_yyyyzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 191); 

                auto tg_yyyzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 192); 

                auto tg_yyyzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 193); 

                auto tg_yyyzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 194); 

                auto tg_yyyzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 195); 

                auto tg_yyyzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 196); 

                auto tg_yyyzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 197); 

                auto tg_yyzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 198); 

                auto tg_yyzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 199); 

                auto tg_yyzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 200); 

                auto tg_yyzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 201); 

                auto tg_yyzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 202); 

                auto tg_yyzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 203); 

                auto tg_yzzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 204); 

                auto tg_yzzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 205); 

                auto tg_yzzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 206); 

                auto tg_yzzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 207); 

                auto tg_yzzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 208); 

                auto tg_yzzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 209); 

                auto tg_zzzzzzz_xx_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 210); 

                auto tg_zzzzzzz_xy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 211); 

                auto tg_zzzzzzz_xz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 212); 

                auto tg_zzzzzzz_yy_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 213); 

                auto tg_zzzzzzz_yz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 214); 

                auto tg_zzzzzzz_zz_1 = primBuffer[pidx_g_7_2_m1].data(216 * idx + 215); 

                // set up pointers to integrals

                auto tg_yyyyyyyy_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 360); 

                auto tg_yyyyyyyy_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 361); 

                auto tg_yyyyyyyy_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 362); 

                auto tg_yyyyyyyy_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 363); 

                auto tg_yyyyyyyy_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 364); 

                auto tg_yyyyyyyy_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 365); 

                auto tg_yyyyyyyy_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 366); 

                auto tg_yyyyyyyy_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 367); 

                auto tg_yyyyyyyy_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 368); 

                auto tg_yyyyyyyy_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 369); 

                auto tg_yyyyyyyz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 370); 

                auto tg_yyyyyyyz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 371); 

                auto tg_yyyyyyyz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 372); 

                auto tg_yyyyyyyz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 373); 

                auto tg_yyyyyyyz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 374); 

                auto tg_yyyyyyyz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 375); 

                auto tg_yyyyyyyz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 376); 

                auto tg_yyyyyyyz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 377); 

                auto tg_yyyyyyyz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 378); 

                auto tg_yyyyyyyz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 379); 

                auto tg_yyyyyyzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 380); 

                auto tg_yyyyyyzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 381); 

                auto tg_yyyyyyzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 382); 

                auto tg_yyyyyyzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 383); 

                auto tg_yyyyyyzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 384); 

                auto tg_yyyyyyzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 385); 

                auto tg_yyyyyyzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 386); 

                auto tg_yyyyyyzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 387); 

                auto tg_yyyyyyzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 388); 

                auto tg_yyyyyyzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 389); 

                auto tg_yyyyyzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 390); 

                auto tg_yyyyyzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 391); 

                auto tg_yyyyyzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 392); 

                auto tg_yyyyyzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 393); 

                auto tg_yyyyyzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 394); 

                auto tg_yyyyyzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 395); 

                auto tg_yyyyyzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 396); 

                auto tg_yyyyyzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 397); 

                auto tg_yyyyyzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 398); 

                auto tg_yyyyyzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 399); 

                auto tg_yyyyzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 400); 

                auto tg_yyyyzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 401); 

                auto tg_yyyyzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 402); 

                auto tg_yyyyzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 403); 

                auto tg_yyyyzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 404); 

                auto tg_yyyyzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 405); 

                auto tg_yyyyzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 406); 

                auto tg_yyyyzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 407); 

                auto tg_yyyyzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 408); 

                auto tg_yyyyzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 409); 

                auto tg_yyyzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 410); 

                auto tg_yyyzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 411); 

                auto tg_yyyzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 412); 

                auto tg_yyyzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 413); 

                auto tg_yyyzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 414); 

                auto tg_yyyzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 415); 

                auto tg_yyyzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 416); 

                auto tg_yyyzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 417); 

                auto tg_yyyzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 418); 

                auto tg_yyyzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 419); 

                auto tg_yyzzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 420); 

                auto tg_yyzzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 421); 

                auto tg_yyzzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 422); 

                auto tg_yyzzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 423); 

                auto tg_yyzzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 424); 

                auto tg_yyzzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 425); 

                auto tg_yyzzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 426); 

                auto tg_yyzzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 427); 

                auto tg_yyzzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 428); 

                auto tg_yyzzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 429); 

                auto tg_yzzzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 430); 

                auto tg_yzzzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 431); 

                auto tg_yzzzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 432); 

                auto tg_yzzzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 433); 

                auto tg_yzzzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 434); 

                auto tg_yzzzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 435); 

                auto tg_yzzzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 436); 

                auto tg_yzzzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 437); 

                auto tg_yzzzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 438); 

                auto tg_yzzzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 439); 

                auto tg_zzzzzzzz_xxx_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 440); 

                auto tg_zzzzzzzz_xxy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 441); 

                auto tg_zzzzzzzz_xxz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 442); 

                auto tg_zzzzzzzz_xyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 443); 

                auto tg_zzzzzzzz_xyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 444); 

                auto tg_zzzzzzzz_xzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 445); 

                auto tg_zzzzzzzz_yyy_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 446); 

                auto tg_zzzzzzzz_yyz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 447); 

                auto tg_zzzzzzzz_yzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 448); 

                auto tg_zzzzzzzz_zzz_0 = primBuffer[pidx_g_8_3_m0].data(450 * idx + 449); 

                // Batch of Integrals (360,450)

                #pragma omp simd aligned(fxn, fza, tg_yyyyyy_xxx_0, tg_yyyyyy_xxx_1, tg_yyyyyy_xxy_0, \
                                         tg_yyyyyy_xxy_1, tg_yyyyyy_xxz_0, tg_yyyyyy_xxz_1, tg_yyyyyy_xyy_0, tg_yyyyyy_xyy_1, \
                                         tg_yyyyyy_xyz_0, tg_yyyyyy_xyz_1, tg_yyyyyy_xzz_0, tg_yyyyyy_xzz_1, tg_yyyyyy_yyy_0, \
                                         tg_yyyyyy_yyy_1, tg_yyyyyy_yyz_0, tg_yyyyyy_yyz_1, tg_yyyyyy_yzz_0, tg_yyyyyy_yzz_1, \
                                         tg_yyyyyy_zzz_0, tg_yyyyyy_zzz_1, tg_yyyyyyy_xx_1, tg_yyyyyyy_xxx_0, \
                                         tg_yyyyyyy_xxx_1, tg_yyyyyyy_xxy_0, tg_yyyyyyy_xxy_1, tg_yyyyyyy_xxz_0, \
                                         tg_yyyyyyy_xxz_1, tg_yyyyyyy_xy_1, tg_yyyyyyy_xyy_0, tg_yyyyyyy_xyy_1, \
                                         tg_yyyyyyy_xyz_0, tg_yyyyyyy_xyz_1, tg_yyyyyyy_xz_1, tg_yyyyyyy_xzz_0, \
                                         tg_yyyyyyy_xzz_1, tg_yyyyyyy_yy_1, tg_yyyyyyy_yyy_0, tg_yyyyyyy_yyy_1, \
                                         tg_yyyyyyy_yyz_0, tg_yyyyyyy_yyz_1, tg_yyyyyyy_yz_1, tg_yyyyyyy_yzz_0, \
                                         tg_yyyyyyy_yzz_1, tg_yyyyyyy_zz_1, tg_yyyyyyy_zzz_0, tg_yyyyyyy_zzz_1, \
                                         tg_yyyyyyyy_xxx_0, tg_yyyyyyyy_xxy_0, tg_yyyyyyyy_xxz_0, tg_yyyyyyyy_xyy_0, \
                                         tg_yyyyyyyy_xyz_0, tg_yyyyyyyy_xzz_0, tg_yyyyyyyy_yyy_0, tg_yyyyyyyy_yyz_0, \
                                         tg_yyyyyyyy_yzz_0, tg_yyyyyyyy_zzz_0, tg_yyyyyyyz_xxx_0, tg_yyyyyyyz_xxy_0, \
                                         tg_yyyyyyyz_xxz_0, tg_yyyyyyyz_xyy_0, tg_yyyyyyyz_xyz_0, tg_yyyyyyyz_xzz_0, \
                                         tg_yyyyyyyz_yyy_0, tg_yyyyyyyz_yyz_0, tg_yyyyyyyz_yzz_0, tg_yyyyyyyz_zzz_0, \
                                         tg_yyyyyyz_xx_1, tg_yyyyyyz_xxx_0, tg_yyyyyyz_xxx_1, tg_yyyyyyz_xxy_0, \
                                         tg_yyyyyyz_xxy_1, tg_yyyyyyz_xxz_0, tg_yyyyyyz_xxz_1, tg_yyyyyyz_xy_1, \
                                         tg_yyyyyyz_xyy_0, tg_yyyyyyz_xyy_1, tg_yyyyyyz_xyz_0, tg_yyyyyyz_xyz_1, \
                                         tg_yyyyyyz_xz_1, tg_yyyyyyz_xzz_0, tg_yyyyyyz_xzz_1, tg_yyyyyyz_yy_1, \
                                         tg_yyyyyyz_yyy_0, tg_yyyyyyz_yyy_1, tg_yyyyyyz_yyz_0, tg_yyyyyyz_yyz_1, \
                                         tg_yyyyyyz_yz_1, tg_yyyyyyz_yzz_0, tg_yyyyyyz_yzz_1, tg_yyyyyyz_zz_1, \
                                         tg_yyyyyyz_zzz_0, tg_yyyyyyz_zzz_1, tg_yyyyyyzz_xxx_0, tg_yyyyyyzz_xxy_0, \
                                         tg_yyyyyyzz_xxz_0, tg_yyyyyyzz_xyy_0, tg_yyyyyyzz_xyz_0, tg_yyyyyyzz_xzz_0, \
                                         tg_yyyyyyzz_yyy_0, tg_yyyyyyzz_yyz_0, tg_yyyyyyzz_yzz_0, tg_yyyyyyzz_zzz_0, \
                                         tg_yyyyyz_xxx_0, tg_yyyyyz_xxx_1, tg_yyyyyz_xxy_0, tg_yyyyyz_xxy_1, tg_yyyyyz_xxz_0, \
                                         tg_yyyyyz_xxz_1, tg_yyyyyz_xyy_0, tg_yyyyyz_xyy_1, tg_yyyyyz_xyz_0, tg_yyyyyz_xyz_1, \
                                         tg_yyyyyz_xzz_0, tg_yyyyyz_xzz_1, tg_yyyyyz_yyy_0, tg_yyyyyz_yyy_1, tg_yyyyyz_yyz_0, \
                                         tg_yyyyyz_yyz_1, tg_yyyyyz_yzz_0, tg_yyyyyz_yzz_1, tg_yyyyyz_zzz_0, tg_yyyyyz_zzz_1, \
                                         tg_yyyyyzz_xx_1, tg_yyyyyzz_xxx_0, tg_yyyyyzz_xxx_1, tg_yyyyyzz_xxy_0, \
                                         tg_yyyyyzz_xxy_1, tg_yyyyyzz_xxz_0, tg_yyyyyzz_xxz_1, tg_yyyyyzz_xy_1, \
                                         tg_yyyyyzz_xyy_0, tg_yyyyyzz_xyy_1, tg_yyyyyzz_xyz_0, tg_yyyyyzz_xyz_1, \
                                         tg_yyyyyzz_xz_1, tg_yyyyyzz_xzz_0, tg_yyyyyzz_xzz_1, tg_yyyyyzz_yy_1, \
                                         tg_yyyyyzz_yyy_0, tg_yyyyyzz_yyy_1, tg_yyyyyzz_yyz_0, tg_yyyyyzz_yyz_1, \
                                         tg_yyyyyzz_yz_1, tg_yyyyyzz_yzz_0, tg_yyyyyzz_yzz_1, tg_yyyyyzz_zz_1, \
                                         tg_yyyyyzz_zzz_0, tg_yyyyyzz_zzz_1, tg_yyyyyzzz_xxx_0, tg_yyyyyzzz_xxy_0, \
                                         tg_yyyyyzzz_xxz_0, tg_yyyyyzzz_xyy_0, tg_yyyyyzzz_xyz_0, tg_yyyyyzzz_xzz_0, \
                                         tg_yyyyyzzz_yyy_0, tg_yyyyyzzz_yyz_0, tg_yyyyyzzz_yzz_0, tg_yyyyyzzz_zzz_0, \
                                         tg_yyyyzz_xxx_0, tg_yyyyzz_xxx_1, tg_yyyyzz_xxy_0, tg_yyyyzz_xxy_1, tg_yyyyzz_xxz_0, \
                                         tg_yyyyzz_xxz_1, tg_yyyyzz_xyy_0, tg_yyyyzz_xyy_1, tg_yyyyzz_xyz_0, tg_yyyyzz_xyz_1, \
                                         tg_yyyyzz_xzz_0, tg_yyyyzz_xzz_1, tg_yyyyzz_yyy_0, tg_yyyyzz_yyy_1, tg_yyyyzz_yyz_0, \
                                         tg_yyyyzz_yyz_1, tg_yyyyzz_yzz_0, tg_yyyyzz_yzz_1, tg_yyyyzz_zzz_0, tg_yyyyzz_zzz_1, \
                                         tg_yyyyzzz_xx_1, tg_yyyyzzz_xxx_0, tg_yyyyzzz_xxx_1, tg_yyyyzzz_xxy_0, \
                                         tg_yyyyzzz_xxy_1, tg_yyyyzzz_xxz_0, tg_yyyyzzz_xxz_1, tg_yyyyzzz_xy_1, \
                                         tg_yyyyzzz_xyy_0, tg_yyyyzzz_xyy_1, tg_yyyyzzz_xyz_0, tg_yyyyzzz_xyz_1, \
                                         tg_yyyyzzz_xz_1, tg_yyyyzzz_xzz_0, tg_yyyyzzz_xzz_1, tg_yyyyzzz_yy_1, \
                                         tg_yyyyzzz_yyy_0, tg_yyyyzzz_yyy_1, tg_yyyyzzz_yyz_0, tg_yyyyzzz_yyz_1, \
                                         tg_yyyyzzz_yz_1, tg_yyyyzzz_yzz_0, tg_yyyyzzz_yzz_1, tg_yyyyzzz_zz_1, \
                                         tg_yyyyzzz_zzz_0, tg_yyyyzzz_zzz_1, tg_yyyyzzzz_xxx_0, tg_yyyyzzzz_xxy_0, \
                                         tg_yyyyzzzz_xxz_0, tg_yyyyzzzz_xyy_0, tg_yyyyzzzz_xyz_0, tg_yyyyzzzz_xzz_0, \
                                         tg_yyyyzzzz_yyy_0, tg_yyyyzzzz_yyz_0, tg_yyyyzzzz_yzz_0, tg_yyyyzzzz_zzz_0, \
                                         tg_yyyzzz_xxx_0, tg_yyyzzz_xxx_1, tg_yyyzzz_xxy_0, tg_yyyzzz_xxy_1, tg_yyyzzz_xxz_0, \
                                         tg_yyyzzz_xxz_1, tg_yyyzzz_xyy_0, tg_yyyzzz_xyy_1, tg_yyyzzz_xyz_0, tg_yyyzzz_xyz_1, \
                                         tg_yyyzzz_xzz_0, tg_yyyzzz_xzz_1, tg_yyyzzz_yyy_0, tg_yyyzzz_yyy_1, tg_yyyzzz_yyz_0, \
                                         tg_yyyzzz_yyz_1, tg_yyyzzz_yzz_0, tg_yyyzzz_yzz_1, tg_yyyzzz_zzz_0, tg_yyyzzz_zzz_1, \
                                         tg_yyyzzzz_xx_1, tg_yyyzzzz_xxx_0, tg_yyyzzzz_xxx_1, tg_yyyzzzz_xxy_0, \
                                         tg_yyyzzzz_xxy_1, tg_yyyzzzz_xxz_0, tg_yyyzzzz_xxz_1, tg_yyyzzzz_xy_1, \
                                         tg_yyyzzzz_xyy_0, tg_yyyzzzz_xyy_1, tg_yyyzzzz_xyz_0, tg_yyyzzzz_xyz_1, \
                                         tg_yyyzzzz_xz_1, tg_yyyzzzz_xzz_0, tg_yyyzzzz_xzz_1, tg_yyyzzzz_yy_1, \
                                         tg_yyyzzzz_yyy_0, tg_yyyzzzz_yyy_1, tg_yyyzzzz_yyz_0, tg_yyyzzzz_yyz_1, \
                                         tg_yyyzzzz_yz_1, tg_yyyzzzz_yzz_0, tg_yyyzzzz_yzz_1, tg_yyyzzzz_zz_1, \
                                         tg_yyyzzzz_zzz_0, tg_yyyzzzz_zzz_1, tg_yyyzzzzz_xxx_0, tg_yyyzzzzz_xxy_0, \
                                         tg_yyyzzzzz_xxz_0, tg_yyyzzzzz_xyy_0, tg_yyyzzzzz_xyz_0, tg_yyyzzzzz_xzz_0, \
                                         tg_yyyzzzzz_yyy_0, tg_yyyzzzzz_yyz_0, tg_yyyzzzzz_yzz_0, tg_yyyzzzzz_zzz_0, \
                                         tg_yyzzzz_xxx_0, tg_yyzzzz_xxx_1, tg_yyzzzz_xxy_0, tg_yyzzzz_xxy_1, tg_yyzzzz_xxz_0, \
                                         tg_yyzzzz_xxz_1, tg_yyzzzz_xyy_0, tg_yyzzzz_xyy_1, tg_yyzzzz_xyz_0, tg_yyzzzz_xyz_1, \
                                         tg_yyzzzz_xzz_0, tg_yyzzzz_xzz_1, tg_yyzzzz_yyy_0, tg_yyzzzz_yyy_1, tg_yyzzzz_yyz_0, \
                                         tg_yyzzzz_yyz_1, tg_yyzzzz_yzz_0, tg_yyzzzz_yzz_1, tg_yyzzzz_zzz_0, tg_yyzzzz_zzz_1, \
                                         tg_yyzzzzz_xx_1, tg_yyzzzzz_xxx_0, tg_yyzzzzz_xxx_1, tg_yyzzzzz_xxy_0, \
                                         tg_yyzzzzz_xxy_1, tg_yyzzzzz_xxz_0, tg_yyzzzzz_xxz_1, tg_yyzzzzz_xy_1, \
                                         tg_yyzzzzz_xyy_0, tg_yyzzzzz_xyy_1, tg_yyzzzzz_xyz_0, tg_yyzzzzz_xyz_1, \
                                         tg_yyzzzzz_xz_1, tg_yyzzzzz_xzz_0, tg_yyzzzzz_xzz_1, tg_yyzzzzz_yy_1, \
                                         tg_yyzzzzz_yyy_0, tg_yyzzzzz_yyy_1, tg_yyzzzzz_yyz_0, tg_yyzzzzz_yyz_1, \
                                         tg_yyzzzzz_yz_1, tg_yyzzzzz_yzz_0, tg_yyzzzzz_yzz_1, tg_yyzzzzz_zz_1, \
                                         tg_yyzzzzz_zzz_0, tg_yyzzzzz_zzz_1, tg_yyzzzzzz_xxx_0, tg_yyzzzzzz_xxy_0, \
                                         tg_yyzzzzzz_xxz_0, tg_yyzzzzzz_xyy_0, tg_yyzzzzzz_xyz_0, tg_yyzzzzzz_xzz_0, \
                                         tg_yyzzzzzz_yyy_0, tg_yyzzzzzz_yyz_0, tg_yyzzzzzz_yzz_0, tg_yyzzzzzz_zzz_0, \
                                         tg_yzzzzz_xxx_0, tg_yzzzzz_xxx_1, tg_yzzzzz_xxy_0, tg_yzzzzz_xxy_1, tg_yzzzzz_xxz_0, \
                                         tg_yzzzzz_xxz_1, tg_yzzzzz_xyy_0, tg_yzzzzz_xyy_1, tg_yzzzzz_xyz_0, tg_yzzzzz_xyz_1, \
                                         tg_yzzzzz_xzz_0, tg_yzzzzz_xzz_1, tg_yzzzzz_yyy_0, tg_yzzzzz_yyy_1, tg_yzzzzz_yyz_0, \
                                         tg_yzzzzz_yyz_1, tg_yzzzzz_yzz_0, tg_yzzzzz_yzz_1, tg_yzzzzz_zzz_0, tg_yzzzzz_zzz_1, \
                                         tg_yzzzzzz_xx_1, tg_yzzzzzz_xxx_0, tg_yzzzzzz_xxx_1, tg_yzzzzzz_xxy_0, \
                                         tg_yzzzzzz_xxy_1, tg_yzzzzzz_xxz_0, tg_yzzzzzz_xxz_1, tg_yzzzzzz_xy_1, \
                                         tg_yzzzzzz_xyy_0, tg_yzzzzzz_xyy_1, tg_yzzzzzz_xyz_0, tg_yzzzzzz_xyz_1, \
                                         tg_yzzzzzz_xz_1, tg_yzzzzzz_xzz_0, tg_yzzzzzz_xzz_1, tg_yzzzzzz_yy_1, \
                                         tg_yzzzzzz_yyy_0, tg_yzzzzzz_yyy_1, tg_yzzzzzz_yyz_0, tg_yzzzzzz_yyz_1, \
                                         tg_yzzzzzz_yz_1, tg_yzzzzzz_yzz_0, tg_yzzzzzz_yzz_1, tg_yzzzzzz_zz_1, \
                                         tg_yzzzzzz_zzz_0, tg_yzzzzzz_zzz_1, tg_yzzzzzzz_xxx_0, tg_yzzzzzzz_xxy_0, \
                                         tg_yzzzzzzz_xxz_0, tg_yzzzzzzz_xyy_0, tg_yzzzzzzz_xyz_0, tg_yzzzzzzz_xzz_0, \
                                         tg_yzzzzzzz_yyy_0, tg_yzzzzzzz_yyz_0, tg_yzzzzzzz_yzz_0, tg_yzzzzzzz_zzz_0, \
                                         tg_zzzzzz_xxx_0, tg_zzzzzz_xxx_1, tg_zzzzzz_xxy_0, tg_zzzzzz_xxy_1, tg_zzzzzz_xxz_0, \
                                         tg_zzzzzz_xxz_1, tg_zzzzzz_xyy_0, tg_zzzzzz_xyy_1, tg_zzzzzz_xyz_0, tg_zzzzzz_xyz_1, \
                                         tg_zzzzzz_xzz_0, tg_zzzzzz_xzz_1, tg_zzzzzz_yyy_0, tg_zzzzzz_yyy_1, tg_zzzzzz_yyz_0, \
                                         tg_zzzzzz_yyz_1, tg_zzzzzz_yzz_0, tg_zzzzzz_yzz_1, tg_zzzzzz_zzz_0, tg_zzzzzz_zzz_1, \
                                         tg_zzzzzzz_xx_1, tg_zzzzzzz_xxx_0, tg_zzzzzzz_xxx_1, tg_zzzzzzz_xxy_0, \
                                         tg_zzzzzzz_xxy_1, tg_zzzzzzz_xxz_0, tg_zzzzzzz_xxz_1, tg_zzzzzzz_xy_1, \
                                         tg_zzzzzzz_xyy_0, tg_zzzzzzz_xyy_1, tg_zzzzzzz_xyz_0, tg_zzzzzzz_xyz_1, \
                                         tg_zzzzzzz_xz_1, tg_zzzzzzz_xzz_0, tg_zzzzzzz_xzz_1, tg_zzzzzzz_yy_1, \
                                         tg_zzzzzzz_yyy_0, tg_zzzzzzz_yyy_1, tg_zzzzzzz_yyz_0, tg_zzzzzzz_yyz_1, \
                                         tg_zzzzzzz_yz_1, tg_zzzzzzz_yzz_0, tg_zzzzzzz_yzz_1, tg_zzzzzzz_zz_1, \
                                         tg_zzzzzzz_zzz_0, tg_zzzzzzz_zzz_1, tg_zzzzzzzz_xxx_0, tg_zzzzzzzz_xxy_0, \
                                         tg_zzzzzzzz_xxz_0, tg_zzzzzzzz_xyy_0, tg_zzzzzzzz_xyz_0, tg_zzzzzzzz_xzz_0, \
                                         tg_zzzzzzzz_yyy_0, tg_zzzzzzzz_yyz_0, tg_zzzzzzzz_yzz_0, tg_zzzzzzzz_zzz_0, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_y[j]; 

                    tg_yyyyyyyy_xxx_0[j] = pb_y * tg_yyyyyyy_xxx_0[j] + fr * tg_yyyyyyy_xxx_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxx_0[j] - tg_yyyyyy_xxx_1[j] * fl1_fza);

                    tg_yyyyyyyy_xxy_0[j] = pb_y * tg_yyyyyyy_xxy_0[j] + fr * tg_yyyyyyy_xxy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxy_0[j] - tg_yyyyyy_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_xx_1[j];

                    tg_yyyyyyyy_xxz_0[j] = pb_y * tg_yyyyyyy_xxz_0[j] + fr * tg_yyyyyyy_xxz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xxz_0[j] - tg_yyyyyy_xxz_1[j] * fl1_fza);

                    tg_yyyyyyyy_xyy_0[j] = pb_y * tg_yyyyyyy_xyy_0[j] + fr * tg_yyyyyyy_xyy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xyy_0[j] - tg_yyyyyy_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyy_xy_1[j];

                    tg_yyyyyyyy_xyz_0[j] = pb_y * tg_yyyyyyy_xyz_0[j] + fr * tg_yyyyyyy_xyz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xyz_0[j] - tg_yyyyyy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_xz_1[j];

                    tg_yyyyyyyy_xzz_0[j] = pb_y * tg_yyyyyyy_xzz_0[j] + fr * tg_yyyyyyy_xzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_xzz_0[j] - tg_yyyyyy_xzz_1[j] * fl1_fza);

                    tg_yyyyyyyy_yyy_0[j] = pb_y * tg_yyyyyyy_yyy_0[j] + fr * tg_yyyyyyy_yyy_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yyy_0[j] - tg_yyyyyy_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyy_yy_1[j];

                    tg_yyyyyyyy_yyz_0[j] = pb_y * tg_yyyyyyy_yyz_0[j] + fr * tg_yyyyyyy_yyz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yyz_0[j] - tg_yyyyyy_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyy_yz_1[j];

                    tg_yyyyyyyy_yzz_0[j] = pb_y * tg_yyyyyyy_yzz_0[j] + fr * tg_yyyyyyy_yzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_yzz_0[j] - tg_yyyyyy_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyy_zz_1[j];

                    tg_yyyyyyyy_zzz_0[j] = pb_y * tg_yyyyyyy_zzz_0[j] + fr * tg_yyyyyyy_zzz_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_zzz_0[j] - tg_yyyyyy_zzz_1[j] * fl1_fza);

                    tg_yyyyyyyz_xxx_0[j] = pb_y * tg_yyyyyyz_xxx_0[j] + fr * tg_yyyyyyz_xxx_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxx_0[j] - tg_yyyyyz_xxx_1[j] * fl1_fza);

                    tg_yyyyyyyz_xxy_0[j] = pb_y * tg_yyyyyyz_xxy_0[j] + fr * tg_yyyyyyz_xxy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxy_0[j] - tg_yyyyyz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_xx_1[j];

                    tg_yyyyyyyz_xxz_0[j] = pb_y * tg_yyyyyyz_xxz_0[j] + fr * tg_yyyyyyz_xxz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xxz_0[j] - tg_yyyyyz_xxz_1[j] * fl1_fza);

                    tg_yyyyyyyz_xyy_0[j] = pb_y * tg_yyyyyyz_xyy_0[j] + fr * tg_yyyyyyz_xyy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xyy_0[j] - tg_yyyyyz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyz_xy_1[j];

                    tg_yyyyyyyz_xyz_0[j] = pb_y * tg_yyyyyyz_xyz_0[j] + fr * tg_yyyyyyz_xyz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xyz_0[j] - tg_yyyyyz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_xz_1[j];

                    tg_yyyyyyyz_xzz_0[j] = pb_y * tg_yyyyyyz_xzz_0[j] + fr * tg_yyyyyyz_xzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_xzz_0[j] - tg_yyyyyz_xzz_1[j] * fl1_fza);

                    tg_yyyyyyyz_yyy_0[j] = pb_y * tg_yyyyyyz_yyy_0[j] + fr * tg_yyyyyyz_yyy_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yyy_0[j] - tg_yyyyyz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyyz_yy_1[j];

                    tg_yyyyyyyz_yyz_0[j] = pb_y * tg_yyyyyyz_yyz_0[j] + fr * tg_yyyyyyz_yyz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yyz_0[j] - tg_yyyyyz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyyz_yz_1[j];

                    tg_yyyyyyyz_yzz_0[j] = pb_y * tg_yyyyyyz_yzz_0[j] + fr * tg_yyyyyyz_yzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_yzz_0[j] - tg_yyyyyz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyyz_zz_1[j];

                    tg_yyyyyyyz_zzz_0[j] = pb_y * tg_yyyyyyz_zzz_0[j] + fr * tg_yyyyyyz_zzz_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_zzz_0[j] - tg_yyyyyz_zzz_1[j] * fl1_fza);

                    tg_yyyyyyzz_xxx_0[j] = pb_y * tg_yyyyyzz_xxx_0[j] + fr * tg_yyyyyzz_xxx_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxx_0[j] - tg_yyyyzz_xxx_1[j] * fl1_fza);

                    tg_yyyyyyzz_xxy_0[j] = pb_y * tg_yyyyyzz_xxy_0[j] + fr * tg_yyyyyzz_xxy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxy_0[j] - tg_yyyyzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_xx_1[j];

                    tg_yyyyyyzz_xxz_0[j] = pb_y * tg_yyyyyzz_xxz_0[j] + fr * tg_yyyyyzz_xxz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xxz_0[j] - tg_yyyyzz_xxz_1[j] * fl1_fza);

                    tg_yyyyyyzz_xyy_0[j] = pb_y * tg_yyyyyzz_xyy_0[j] + fr * tg_yyyyyzz_xyy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xyy_0[j] - tg_yyyyzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyzz_xy_1[j];

                    tg_yyyyyyzz_xyz_0[j] = pb_y * tg_yyyyyzz_xyz_0[j] + fr * tg_yyyyyzz_xyz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xyz_0[j] - tg_yyyyzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_xz_1[j];

                    tg_yyyyyyzz_xzz_0[j] = pb_y * tg_yyyyyzz_xzz_0[j] + fr * tg_yyyyyzz_xzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_xzz_0[j] - tg_yyyyzz_xzz_1[j] * fl1_fza);

                    tg_yyyyyyzz_yyy_0[j] = pb_y * tg_yyyyyzz_yyy_0[j] + fr * tg_yyyyyzz_yyy_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yyy_0[j] - tg_yyyyzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyyzz_yy_1[j];

                    tg_yyyyyyzz_yyz_0[j] = pb_y * tg_yyyyyzz_yyz_0[j] + fr * tg_yyyyyzz_yyz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yyz_0[j] - tg_yyyyzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyyzz_yz_1[j];

                    tg_yyyyyyzz_yzz_0[j] = pb_y * tg_yyyyyzz_yzz_0[j] + fr * tg_yyyyyzz_yzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_yzz_0[j] - tg_yyyyzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyyzz_zz_1[j];

                    tg_yyyyyyzz_zzz_0[j] = pb_y * tg_yyyyyzz_zzz_0[j] + fr * tg_yyyyyzz_zzz_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_zzz_0[j] - tg_yyyyzz_zzz_1[j] * fl1_fza);

                    tg_yyyyyzzz_xxx_0[j] = pb_y * tg_yyyyzzz_xxx_0[j] + fr * tg_yyyyzzz_xxx_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxx_0[j] - tg_yyyzzz_xxx_1[j] * fl1_fza);

                    tg_yyyyyzzz_xxy_0[j] = pb_y * tg_yyyyzzz_xxy_0[j] + fr * tg_yyyyzzz_xxy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxy_0[j] - tg_yyyzzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_xx_1[j];

                    tg_yyyyyzzz_xxz_0[j] = pb_y * tg_yyyyzzz_xxz_0[j] + fr * tg_yyyyzzz_xxz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xxz_0[j] - tg_yyyzzz_xxz_1[j] * fl1_fza);

                    tg_yyyyyzzz_xyy_0[j] = pb_y * tg_yyyyzzz_xyy_0[j] + fr * tg_yyyyzzz_xyy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xyy_0[j] - tg_yyyzzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzzz_xy_1[j];

                    tg_yyyyyzzz_xyz_0[j] = pb_y * tg_yyyyzzz_xyz_0[j] + fr * tg_yyyyzzz_xyz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xyz_0[j] - tg_yyyzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_xz_1[j];

                    tg_yyyyyzzz_xzz_0[j] = pb_y * tg_yyyyzzz_xzz_0[j] + fr * tg_yyyyzzz_xzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_xzz_0[j] - tg_yyyzzz_xzz_1[j] * fl1_fza);

                    tg_yyyyyzzz_yyy_0[j] = pb_y * tg_yyyyzzz_yyy_0[j] + fr * tg_yyyyzzz_yyy_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yyy_0[j] - tg_yyyzzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyyzzz_yy_1[j];

                    tg_yyyyyzzz_yyz_0[j] = pb_y * tg_yyyyzzz_yyz_0[j] + fr * tg_yyyyzzz_yyz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yyz_0[j] - tg_yyyzzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyyzzz_yz_1[j];

                    tg_yyyyyzzz_yzz_0[j] = pb_y * tg_yyyyzzz_yzz_0[j] + fr * tg_yyyyzzz_yzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_yzz_0[j] - tg_yyyzzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyyzzz_zz_1[j];

                    tg_yyyyyzzz_zzz_0[j] = pb_y * tg_yyyyzzz_zzz_0[j] + fr * tg_yyyyzzz_zzz_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_zzz_0[j] - tg_yyyzzz_zzz_1[j] * fl1_fza);

                    tg_yyyyzzzz_xxx_0[j] = pb_y * tg_yyyzzzz_xxx_0[j] + fr * tg_yyyzzzz_xxx_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxx_0[j] - tg_yyzzzz_xxx_1[j] * fl1_fza);

                    tg_yyyyzzzz_xxy_0[j] = pb_y * tg_yyyzzzz_xxy_0[j] + fr * tg_yyyzzzz_xxy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxy_0[j] - tg_yyzzzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_xx_1[j];

                    tg_yyyyzzzz_xxz_0[j] = pb_y * tg_yyyzzzz_xxz_0[j] + fr * tg_yyyzzzz_xxz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xxz_0[j] - tg_yyzzzz_xxz_1[j] * fl1_fza);

                    tg_yyyyzzzz_xyy_0[j] = pb_y * tg_yyyzzzz_xyy_0[j] + fr * tg_yyyzzzz_xyy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xyy_0[j] - tg_yyzzzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzzz_xy_1[j];

                    tg_yyyyzzzz_xyz_0[j] = pb_y * tg_yyyzzzz_xyz_0[j] + fr * tg_yyyzzzz_xyz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xyz_0[j] - tg_yyzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_xz_1[j];

                    tg_yyyyzzzz_xzz_0[j] = pb_y * tg_yyyzzzz_xzz_0[j] + fr * tg_yyyzzzz_xzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_xzz_0[j] - tg_yyzzzz_xzz_1[j] * fl1_fza);

                    tg_yyyyzzzz_yyy_0[j] = pb_y * tg_yyyzzzz_yyy_0[j] + fr * tg_yyyzzzz_yyy_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yyy_0[j] - tg_yyzzzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyyzzzz_yy_1[j];

                    tg_yyyyzzzz_yyz_0[j] = pb_y * tg_yyyzzzz_yyz_0[j] + fr * tg_yyyzzzz_yyz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yyz_0[j] - tg_yyzzzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyyzzzz_yz_1[j];

                    tg_yyyyzzzz_yzz_0[j] = pb_y * tg_yyyzzzz_yzz_0[j] + fr * tg_yyyzzzz_yzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_yzz_0[j] - tg_yyzzzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyyzzzz_zz_1[j];

                    tg_yyyyzzzz_zzz_0[j] = pb_y * tg_yyyzzzz_zzz_0[j] + fr * tg_yyyzzzz_zzz_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_zzz_0[j] - tg_yyzzzz_zzz_1[j] * fl1_fza);

                    tg_yyyzzzzz_xxx_0[j] = pb_y * tg_yyzzzzz_xxx_0[j] + fr * tg_yyzzzzz_xxx_1[j] + fl1_fx * (tg_yzzzzz_xxx_0[j] - tg_yzzzzz_xxx_1[j] * fl1_fza);

                    tg_yyyzzzzz_xxy_0[j] = pb_y * tg_yyzzzzz_xxy_0[j] + fr * tg_yyzzzzz_xxy_1[j] + fl1_fx * (tg_yzzzzz_xxy_0[j] - tg_yzzzzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_xx_1[j];

                    tg_yyyzzzzz_xxz_0[j] = pb_y * tg_yyzzzzz_xxz_0[j] + fr * tg_yyzzzzz_xxz_1[j] + fl1_fx * (tg_yzzzzz_xxz_0[j] - tg_yzzzzz_xxz_1[j] * fl1_fza);

                    tg_yyyzzzzz_xyy_0[j] = pb_y * tg_yyzzzzz_xyy_0[j] + fr * tg_yyzzzzz_xyy_1[j] + fl1_fx * (tg_yzzzzz_xyy_0[j] - tg_yzzzzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzzz_xy_1[j];

                    tg_yyyzzzzz_xyz_0[j] = pb_y * tg_yyzzzzz_xyz_0[j] + fr * tg_yyzzzzz_xyz_1[j] + fl1_fx * (tg_yzzzzz_xyz_0[j] - tg_yzzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_xz_1[j];

                    tg_yyyzzzzz_xzz_0[j] = pb_y * tg_yyzzzzz_xzz_0[j] + fr * tg_yyzzzzz_xzz_1[j] + fl1_fx * (tg_yzzzzz_xzz_0[j] - tg_yzzzzz_xzz_1[j] * fl1_fza);

                    tg_yyyzzzzz_yyy_0[j] = pb_y * tg_yyzzzzz_yyy_0[j] + fr * tg_yyzzzzz_yyy_1[j] + fl1_fx * (tg_yzzzzz_yyy_0[j] - tg_yzzzzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyzzzzz_yy_1[j];

                    tg_yyyzzzzz_yyz_0[j] = pb_y * tg_yyzzzzz_yyz_0[j] + fr * tg_yyzzzzz_yyz_1[j] + fl1_fx * (tg_yzzzzz_yyz_0[j] - tg_yzzzzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyzzzzz_yz_1[j];

                    tg_yyyzzzzz_yzz_0[j] = pb_y * tg_yyzzzzz_yzz_0[j] + fr * tg_yyzzzzz_yzz_1[j] + fl1_fx * (tg_yzzzzz_yzz_0[j] - tg_yzzzzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyzzzzz_zz_1[j];

                    tg_yyyzzzzz_zzz_0[j] = pb_y * tg_yyzzzzz_zzz_0[j] + fr * tg_yyzzzzz_zzz_1[j] + fl1_fx * (tg_yzzzzz_zzz_0[j] - tg_yzzzzz_zzz_1[j] * fl1_fza);

                    tg_yyzzzzzz_xxx_0[j] = pb_y * tg_yzzzzzz_xxx_0[j] + fr * tg_yzzzzzz_xxx_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxx_0[j] - tg_zzzzzz_xxx_1[j] * fl1_fza);

                    tg_yyzzzzzz_xxy_0[j] = pb_y * tg_yzzzzzz_xxy_0[j] + fr * tg_yzzzzzz_xxy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxy_0[j] - tg_zzzzzz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_xx_1[j];

                    tg_yyzzzzzz_xxz_0[j] = pb_y * tg_yzzzzzz_xxz_0[j] + fr * tg_yzzzzzz_xxz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xxz_0[j] - tg_zzzzzz_xxz_1[j] * fl1_fza);

                    tg_yyzzzzzz_xyy_0[j] = pb_y * tg_yzzzzzz_xyy_0[j] + fr * tg_yzzzzzz_xyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyy_0[j] - tg_zzzzzz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzzz_xy_1[j];

                    tg_yyzzzzzz_xyz_0[j] = pb_y * tg_yzzzzzz_xyz_0[j] + fr * tg_yzzzzzz_xyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xyz_0[j] - tg_zzzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_xz_1[j];

                    tg_yyzzzzzz_xzz_0[j] = pb_y * tg_yzzzzzz_xzz_0[j] + fr * tg_yzzzzzz_xzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_xzz_0[j] - tg_zzzzzz_xzz_1[j] * fl1_fza);

                    tg_yyzzzzzz_yyy_0[j] = pb_y * tg_yzzzzzz_yyy_0[j] + fr * tg_yzzzzzz_yyy_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyy_0[j] - tg_zzzzzz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzzzzzz_yy_1[j];

                    tg_yyzzzzzz_yyz_0[j] = pb_y * tg_yzzzzzz_yyz_0[j] + fr * tg_yzzzzzz_yyz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yyz_0[j] - tg_zzzzzz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yzzzzzz_yz_1[j];

                    tg_yyzzzzzz_yzz_0[j] = pb_y * tg_yzzzzzz_yzz_0[j] + fr * tg_yzzzzzz_yzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_yzz_0[j] - tg_zzzzzz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzzzzzz_zz_1[j];

                    tg_yyzzzzzz_zzz_0[j] = pb_y * tg_yzzzzzz_zzz_0[j] + fr * tg_yzzzzzz_zzz_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_zzz_0[j] - tg_zzzzzz_zzz_1[j] * fl1_fza);

                    tg_yzzzzzzz_xxx_0[j] = pb_y * tg_zzzzzzz_xxx_0[j] + fr * tg_zzzzzzz_xxx_1[j];

                    tg_yzzzzzzz_xxy_0[j] = pb_y * tg_zzzzzzz_xxy_0[j] + fr * tg_zzzzzzz_xxy_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xx_1[j];

                    tg_yzzzzzzz_xxz_0[j] = pb_y * tg_zzzzzzz_xxz_0[j] + fr * tg_zzzzzzz_xxz_1[j];

                    tg_yzzzzzzz_xyy_0[j] = pb_y * tg_zzzzzzz_xyy_0[j] + fr * tg_zzzzzzz_xyy_1[j] + fl1_fxn * tg_zzzzzzz_xy_1[j];

                    tg_yzzzzzzz_xyz_0[j] = pb_y * tg_zzzzzzz_xyz_0[j] + fr * tg_zzzzzzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_xz_1[j];

                    tg_yzzzzzzz_xzz_0[j] = pb_y * tg_zzzzzzz_xzz_0[j] + fr * tg_zzzzzzz_xzz_1[j];

                    tg_yzzzzzzz_yyy_0[j] = pb_y * tg_zzzzzzz_yyy_0[j] + fr * tg_zzzzzzz_yyy_1[j] + 1.5 * fl1_fxn * tg_zzzzzzz_yy_1[j];

                    tg_yzzzzzzz_yyz_0[j] = pb_y * tg_zzzzzzz_yyz_0[j] + fr * tg_zzzzzzz_yyz_1[j] + fl1_fxn * tg_zzzzzzz_yz_1[j];

                    tg_yzzzzzzz_yzz_0[j] = pb_y * tg_zzzzzzz_yzz_0[j] + fr * tg_zzzzzzz_yzz_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_zz_1[j];

                    tg_yzzzzzzz_zzz_0[j] = pb_y * tg_zzzzzzz_zzz_0[j] + fr * tg_zzzzzzz_zzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzzzzzz_xxx_0[j] = pb_z * tg_zzzzzzz_xxx_0[j] + fr * tg_zzzzzzz_xxx_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxx_0[j] - tg_zzzzzz_xxx_1[j] * fl1_fza);

                    tg_zzzzzzzz_xxy_0[j] = pb_z * tg_zzzzzzz_xxy_0[j] + fr * tg_zzzzzzz_xxy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxy_0[j] - tg_zzzzzz_xxy_1[j] * fl1_fza);

                    tg_zzzzzzzz_xxz_0[j] = pb_z * tg_zzzzzzz_xxz_0[j] + fr * tg_zzzzzzz_xxz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xxz_0[j] - tg_zzzzzz_xxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_xx_1[j];

                    tg_zzzzzzzz_xyy_0[j] = pb_z * tg_zzzzzzz_xyy_0[j] + fr * tg_zzzzzzz_xyy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xyy_0[j] - tg_zzzzzz_xyy_1[j] * fl1_fza);

                    tg_zzzzzzzz_xyz_0[j] = pb_z * tg_zzzzzzz_xyz_0[j] + fr * tg_zzzzzzz_xyz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xyz_0[j] - tg_zzzzzz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_xy_1[j];

                    tg_zzzzzzzz_xzz_0[j] = pb_z * tg_zzzzzzz_xzz_0[j] + fr * tg_zzzzzzz_xzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_xzz_0[j] - tg_zzzzzz_xzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzzz_xz_1[j];

                    tg_zzzzzzzz_yyy_0[j] = pb_z * tg_zzzzzzz_yyy_0[j] + fr * tg_zzzzzzz_yyy_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yyy_0[j] - tg_zzzzzz_yyy_1[j] * fl1_fza);

                    tg_zzzzzzzz_yyz_0[j] = pb_z * tg_zzzzzzz_yyz_0[j] + fr * tg_zzzzzzz_yyz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yyz_0[j] - tg_zzzzzz_yyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzzzzzz_yy_1[j];

                    tg_zzzzzzzz_yzz_0[j] = pb_z * tg_zzzzzzz_yzz_0[j] + fr * tg_zzzzzzz_yzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_yzz_0[j] - tg_zzzzzz_yzz_1[j] * fl1_fza) + fl1_fxn * tg_zzzzzzz_yz_1[j];

                    tg_zzzzzzzz_zzz_0[j] = pb_z * tg_zzzzzzz_zzz_0[j] + fr * tg_zzzzzzz_zzz_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_zzz_0[j] - tg_zzzzzz_zzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzzzzzz_zz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

