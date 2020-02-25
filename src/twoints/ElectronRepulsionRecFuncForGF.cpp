//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionRecFuncForGF.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSGSF(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSGSF_0_75(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSGSF_75_150(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForSGSF_0_75(      CMemBlock2D<double>* primBuffer,
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
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tg_xxx_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx); 

                auto tg_xxx_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 1); 

                auto tg_xxx_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 2); 

                auto tg_xxx_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 3); 

                auto tg_xxx_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 4); 

                auto tg_xxx_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 5); 

                auto tg_xxx_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 6); 

                auto tg_xxx_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 7); 

                auto tg_xxx_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 8); 

                auto tg_xxx_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 9); 

                auto tg_xxy_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 10); 

                auto tg_xxy_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 11); 

                auto tg_xxy_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 12); 

                auto tg_xxy_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 13); 

                auto tg_xxy_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 14); 

                auto tg_xxy_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 15); 

                auto tg_xxy_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 16); 

                auto tg_xxy_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 17); 

                auto tg_xxy_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 18); 

                auto tg_xxy_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 19); 

                auto tg_xxz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 20); 

                auto tg_xxz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 21); 

                auto tg_xxz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 22); 

                auto tg_xxz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 23); 

                auto tg_xxz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 24); 

                auto tg_xxz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 25); 

                auto tg_xxz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 26); 

                auto tg_xxz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 27); 

                auto tg_xxz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 28); 

                auto tg_xxz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 29); 

                auto tg_xyy_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 30); 

                auto tg_xyy_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 31); 

                auto tg_xyy_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 32); 

                auto tg_xyy_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 33); 

                auto tg_xyy_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 34); 

                auto tg_xyy_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 35); 

                auto tg_xyy_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 36); 

                auto tg_xyy_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 37); 

                auto tg_xyy_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 38); 

                auto tg_xyy_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 39); 

                auto tg_xyz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 40); 

                auto tg_xyz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 41); 

                auto tg_xyz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 42); 

                auto tg_xyz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 43); 

                auto tg_xyz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 44); 

                auto tg_xyz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 45); 

                auto tg_xyz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 46); 

                auto tg_xyz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 47); 

                auto tg_xyz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 48); 

                auto tg_xyz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 49); 

                auto tg_xzz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 50); 

                auto tg_xzz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 51); 

                auto tg_xzz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 52); 

                auto tg_xzz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 53); 

                auto tg_xzz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 54); 

                auto tg_xzz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 55); 

                auto tg_xzz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 56); 

                auto tg_xzz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 57); 

                auto tg_xzz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 58); 

                auto tg_xzz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 59); 

                auto tg_yyy_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 60); 

                auto tg_yyy_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 61); 

                auto tg_yyy_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 62); 

                auto tg_yyy_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 63); 

                auto tg_yyy_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 64); 

                auto tg_yyy_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 65); 

                auto tg_yyy_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 66); 

                auto tg_yyy_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 67); 

                auto tg_yyy_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 68); 

                auto tg_yyy_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 69); 

                auto tg_yyz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 70); 

                auto tg_yyz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 71); 

                auto tg_yyz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 72); 

                auto tg_yyz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 73); 

                auto tg_yyz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 74); 

                auto tg_xxx_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx); 

                auto tg_xxx_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 1); 

                auto tg_xxx_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 2); 

                auto tg_xxx_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 3); 

                auto tg_xxx_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 4); 

                auto tg_xxx_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 5); 

                auto tg_xxx_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 6); 

                auto tg_xxx_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 7); 

                auto tg_xxx_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 8); 

                auto tg_xxx_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 9); 

                auto tg_xxy_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 10); 

                auto tg_xxy_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 11); 

                auto tg_xxy_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 12); 

                auto tg_xxy_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 13); 

                auto tg_xxy_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 14); 

                auto tg_xxy_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 15); 

                auto tg_xxy_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 16); 

                auto tg_xxy_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 17); 

                auto tg_xxy_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 18); 

                auto tg_xxy_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 19); 

                auto tg_xxz_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 20); 

                auto tg_xxz_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 21); 

                auto tg_xxz_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 22); 

                auto tg_xxz_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 23); 

                auto tg_xxz_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 24); 

                auto tg_xxz_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 25); 

                auto tg_xxz_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 26); 

                auto tg_xxz_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 27); 

                auto tg_xxz_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 28); 

                auto tg_xxz_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 29); 

                auto tg_xyy_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 30); 

                auto tg_xyy_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 31); 

                auto tg_xyy_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 32); 

                auto tg_xyy_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 33); 

                auto tg_xyy_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 34); 

                auto tg_xyy_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 35); 

                auto tg_xyy_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 36); 

                auto tg_xyy_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 37); 

                auto tg_xyy_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 38); 

                auto tg_xyy_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 39); 

                auto tg_xyz_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 40); 

                auto tg_xyz_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 41); 

                auto tg_xyz_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 42); 

                auto tg_xyz_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 43); 

                auto tg_xyz_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 44); 

                auto tg_xyz_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 45); 

                auto tg_xyz_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 46); 

                auto tg_xyz_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 47); 

                auto tg_xyz_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 48); 

                auto tg_xyz_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 49); 

                auto tg_xzz_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 50); 

                auto tg_xzz_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 51); 

                auto tg_xzz_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 52); 

                auto tg_xzz_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 53); 

                auto tg_xzz_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 54); 

                auto tg_xzz_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 55); 

                auto tg_xzz_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 56); 

                auto tg_xzz_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 57); 

                auto tg_xzz_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 58); 

                auto tg_xzz_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 59); 

                auto tg_yyy_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 60); 

                auto tg_yyy_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 61); 

                auto tg_yyy_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 62); 

                auto tg_yyy_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 63); 

                auto tg_yyy_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 64); 

                auto tg_yyy_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 65); 

                auto tg_yyy_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 66); 

                auto tg_yyy_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 67); 

                auto tg_yyy_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 68); 

                auto tg_yyy_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 69); 

                auto tg_yyz_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 70); 

                auto tg_yyz_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 71); 

                auto tg_yyz_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 72); 

                auto tg_yyz_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 73); 

                auto tg_yyz_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 74); 

                auto tg_xx_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx); 

                auto tg_xx_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 1); 

                auto tg_xx_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 2); 

                auto tg_xx_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 3); 

                auto tg_xx_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 4); 

                auto tg_xx_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 5); 

                auto tg_xx_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 6); 

                auto tg_xx_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 7); 

                auto tg_xx_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 8); 

                auto tg_xx_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 9); 

                auto tg_xy_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 10); 

                auto tg_xy_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 11); 

                auto tg_xy_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 12); 

                auto tg_xy_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 13); 

                auto tg_xy_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 14); 

                auto tg_xy_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 15); 

                auto tg_xy_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 16); 

                auto tg_xy_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 17); 

                auto tg_xy_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 18); 

                auto tg_xy_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 19); 

                auto tg_xz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 20); 

                auto tg_xz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 21); 

                auto tg_xz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 22); 

                auto tg_xz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 23); 

                auto tg_xz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 24); 

                auto tg_xz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 25); 

                auto tg_xz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 26); 

                auto tg_xz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 27); 

                auto tg_xz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 28); 

                auto tg_xz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 29); 

                auto tg_yy_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 30); 

                auto tg_yy_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 31); 

                auto tg_yy_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 32); 

                auto tg_yy_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 33); 

                auto tg_yy_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 34); 

                auto tg_yy_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 35); 

                auto tg_yy_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 36); 

                auto tg_yy_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 37); 

                auto tg_yy_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 38); 

                auto tg_yy_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 39); 

                auto tg_yz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 40); 

                auto tg_yz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 41); 

                auto tg_yz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 42); 

                auto tg_yz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 43); 

                auto tg_yz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 44); 

                auto tg_yz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 45); 

                auto tg_yz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 46); 

                auto tg_yz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 47); 

                auto tg_yz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 48); 

                auto tg_yz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 49); 

                auto tg_zz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 50); 

                auto tg_zz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 51); 

                auto tg_zz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 52); 

                auto tg_zz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 53); 

                auto tg_zz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 54); 

                auto tg_zz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 55); 

                auto tg_zz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 56); 

                auto tg_zz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 57); 

                auto tg_zz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 58); 

                auto tg_zz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 59); 

                auto tg_xx_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx); 

                auto tg_xx_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 1); 

                auto tg_xx_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 2); 

                auto tg_xx_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 3); 

                auto tg_xx_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 4); 

                auto tg_xx_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 5); 

                auto tg_xx_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 6); 

                auto tg_xx_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 7); 

                auto tg_xx_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 8); 

                auto tg_xx_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 9); 

                auto tg_xy_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 10); 

                auto tg_xy_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 11); 

                auto tg_xy_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 12); 

                auto tg_xy_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 13); 

                auto tg_xy_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 14); 

                auto tg_xy_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 15); 

                auto tg_xy_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 16); 

                auto tg_xy_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 17); 

                auto tg_xy_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 18); 

                auto tg_xy_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 19); 

                auto tg_xz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 20); 

                auto tg_xz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 21); 

                auto tg_xz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 22); 

                auto tg_xz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 23); 

                auto tg_xz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 24); 

                auto tg_xz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 25); 

                auto tg_xz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 26); 

                auto tg_xz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 27); 

                auto tg_xz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 28); 

                auto tg_xz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 29); 

                auto tg_yy_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 30); 

                auto tg_yy_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 31); 

                auto tg_yy_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 32); 

                auto tg_yy_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 33); 

                auto tg_yy_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 34); 

                auto tg_yy_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 35); 

                auto tg_yy_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 36); 

                auto tg_yy_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 37); 

                auto tg_yy_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 38); 

                auto tg_yy_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 39); 

                auto tg_yz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 40); 

                auto tg_yz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 41); 

                auto tg_yz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 42); 

                auto tg_yz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 43); 

                auto tg_yz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 44); 

                auto tg_yz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 45); 

                auto tg_yz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 46); 

                auto tg_yz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 47); 

                auto tg_yz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 48); 

                auto tg_yz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 49); 

                auto tg_zz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 50); 

                auto tg_zz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 51); 

                auto tg_zz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 52); 

                auto tg_zz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 53); 

                auto tg_zz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 54); 

                auto tg_zz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 55); 

                auto tg_zz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 56); 

                auto tg_zz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 57); 

                auto tg_zz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 58); 

                auto tg_zz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 59); 

                auto tg_xxx_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx); 

                auto tg_xxx_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 1); 

                auto tg_xxx_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 2); 

                auto tg_xxx_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 3); 

                auto tg_xxx_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 4); 

                auto tg_xxx_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 5); 

                auto tg_xxy_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 6); 

                auto tg_xxy_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 7); 

                auto tg_xxy_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 8); 

                auto tg_xxy_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 9); 

                auto tg_xxy_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 10); 

                auto tg_xxy_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 11); 

                auto tg_xxz_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 12); 

                auto tg_xxz_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 13); 

                auto tg_xxz_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 14); 

                auto tg_xxz_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 15); 

                auto tg_xxz_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 16); 

                auto tg_xxz_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 17); 

                auto tg_xyy_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 18); 

                auto tg_xyy_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 19); 

                auto tg_xyy_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 20); 

                auto tg_xyy_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 21); 

                auto tg_xyy_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 22); 

                auto tg_xyy_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 23); 

                auto tg_xyz_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 24); 

                auto tg_xyz_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 25); 

                auto tg_xyz_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 26); 

                auto tg_xyz_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 27); 

                auto tg_xyz_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 28); 

                auto tg_xyz_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 29); 

                auto tg_xzz_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 30); 

                auto tg_xzz_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 31); 

                auto tg_xzz_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 32); 

                auto tg_xzz_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 33); 

                auto tg_xzz_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 34); 

                auto tg_xzz_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 35); 

                auto tg_yyy_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 36); 

                auto tg_yyy_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 37); 

                auto tg_yyy_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 38); 

                auto tg_yyy_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 39); 

                auto tg_yyy_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 40); 

                auto tg_yyy_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 41); 

                auto tg_yyz_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 42); 

                auto tg_yyz_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 43); 

                auto tg_yyz_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 44); 

                auto tg_yyz_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 45); 

                auto tg_yyz_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 46); 

                // set up pointers to integrals

                auto tg_xxxx_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx); 

                auto tg_xxxx_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 1); 

                auto tg_xxxx_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 2); 

                auto tg_xxxx_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 3); 

                auto tg_xxxx_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 4); 

                auto tg_xxxx_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 5); 

                auto tg_xxxx_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 6); 

                auto tg_xxxx_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 7); 

                auto tg_xxxx_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 8); 

                auto tg_xxxx_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 9); 

                auto tg_xxxy_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 10); 

                auto tg_xxxy_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 11); 

                auto tg_xxxy_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 12); 

                auto tg_xxxy_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 13); 

                auto tg_xxxy_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 14); 

                auto tg_xxxy_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 15); 

                auto tg_xxxy_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 16); 

                auto tg_xxxy_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 17); 

                auto tg_xxxy_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 18); 

                auto tg_xxxy_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 19); 

                auto tg_xxxz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 20); 

                auto tg_xxxz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 21); 

                auto tg_xxxz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 22); 

                auto tg_xxxz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 23); 

                auto tg_xxxz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 24); 

                auto tg_xxxz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 25); 

                auto tg_xxxz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 26); 

                auto tg_xxxz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 27); 

                auto tg_xxxz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 28); 

                auto tg_xxxz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 29); 

                auto tg_xxyy_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 30); 

                auto tg_xxyy_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 31); 

                auto tg_xxyy_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 32); 

                auto tg_xxyy_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 33); 

                auto tg_xxyy_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 34); 

                auto tg_xxyy_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 35); 

                auto tg_xxyy_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 36); 

                auto tg_xxyy_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 37); 

                auto tg_xxyy_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 38); 

                auto tg_xxyy_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 39); 

                auto tg_xxyz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 40); 

                auto tg_xxyz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 41); 

                auto tg_xxyz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 42); 

                auto tg_xxyz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 43); 

                auto tg_xxyz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 44); 

                auto tg_xxyz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 45); 

                auto tg_xxyz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 46); 

                auto tg_xxyz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 47); 

                auto tg_xxyz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 48); 

                auto tg_xxyz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 49); 

                auto tg_xxzz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 50); 

                auto tg_xxzz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 51); 

                auto tg_xxzz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 52); 

                auto tg_xxzz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 53); 

                auto tg_xxzz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 54); 

                auto tg_xxzz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 55); 

                auto tg_xxzz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 56); 

                auto tg_xxzz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 57); 

                auto tg_xxzz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 58); 

                auto tg_xxzz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 59); 

                auto tg_xyyy_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 60); 

                auto tg_xyyy_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 61); 

                auto tg_xyyy_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 62); 

                auto tg_xyyy_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 63); 

                auto tg_xyyy_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 64); 

                auto tg_xyyy_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 65); 

                auto tg_xyyy_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 66); 

                auto tg_xyyy_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 67); 

                auto tg_xyyy_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 68); 

                auto tg_xyyy_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 69); 

                auto tg_xyyz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 70); 

                auto tg_xyyz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 71); 

                auto tg_xyyz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 72); 

                auto tg_xyyz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 73); 

                auto tg_xyyz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 74); 

                // Batch of Integrals (0,75)

                #pragma omp simd aligned(fxn, fza, tg_xx_xxx_0, tg_xx_xxx_1, tg_xx_xxy_0, tg_xx_xxy_1, \
                                         tg_xx_xxz_0, tg_xx_xxz_1, tg_xx_xyy_0, tg_xx_xyy_1, tg_xx_xyz_0, tg_xx_xyz_1, \
                                         tg_xx_xzz_0, tg_xx_xzz_1, tg_xx_yyy_0, tg_xx_yyy_1, tg_xx_yyz_0, tg_xx_yyz_1, \
                                         tg_xx_yzz_0, tg_xx_yzz_1, tg_xx_zzz_0, tg_xx_zzz_1, tg_xxx_xx_1, tg_xxx_xxx_0, \
                                         tg_xxx_xxx_1, tg_xxx_xxy_0, tg_xxx_xxy_1, tg_xxx_xxz_0, tg_xxx_xxz_1, tg_xxx_xy_1, \
                                         tg_xxx_xyy_0, tg_xxx_xyy_1, tg_xxx_xyz_0, tg_xxx_xyz_1, tg_xxx_xz_1, tg_xxx_xzz_0, \
                                         tg_xxx_xzz_1, tg_xxx_yy_1, tg_xxx_yyy_0, tg_xxx_yyy_1, tg_xxx_yyz_0, tg_xxx_yyz_1, \
                                         tg_xxx_yz_1, tg_xxx_yzz_0, tg_xxx_yzz_1, tg_xxx_zz_1, tg_xxx_zzz_0, tg_xxx_zzz_1, \
                                         tg_xxxx_xxx_0, tg_xxxx_xxy_0, tg_xxxx_xxz_0, tg_xxxx_xyy_0, tg_xxxx_xyz_0, \
                                         tg_xxxx_xzz_0, tg_xxxx_yyy_0, tg_xxxx_yyz_0, tg_xxxx_yzz_0, tg_xxxx_zzz_0, \
                                         tg_xxxy_xxx_0, tg_xxxy_xxy_0, tg_xxxy_xxz_0, tg_xxxy_xyy_0, tg_xxxy_xyz_0, \
                                         tg_xxxy_xzz_0, tg_xxxy_yyy_0, tg_xxxy_yyz_0, tg_xxxy_yzz_0, tg_xxxy_zzz_0, \
                                         tg_xxxz_xxx_0, tg_xxxz_xxy_0, tg_xxxz_xxz_0, tg_xxxz_xyy_0, tg_xxxz_xyz_0, \
                                         tg_xxxz_xzz_0, tg_xxxz_yyy_0, tg_xxxz_yyz_0, tg_xxxz_yzz_0, tg_xxxz_zzz_0, \
                                         tg_xxy_xx_1, tg_xxy_xxx_0, tg_xxy_xxx_1, tg_xxy_xxy_0, tg_xxy_xxy_1, tg_xxy_xxz_0, \
                                         tg_xxy_xxz_1, tg_xxy_xy_1, tg_xxy_xyy_0, tg_xxy_xyy_1, tg_xxy_xyz_0, tg_xxy_xyz_1, \
                                         tg_xxy_xz_1, tg_xxy_xzz_0, tg_xxy_xzz_1, tg_xxy_yy_1, tg_xxy_yyy_0, tg_xxy_yyy_1, \
                                         tg_xxy_yyz_0, tg_xxy_yyz_1, tg_xxy_yz_1, tg_xxy_yzz_0, tg_xxy_yzz_1, tg_xxy_zz_1, \
                                         tg_xxy_zzz_0, tg_xxy_zzz_1, tg_xxyy_xxx_0, tg_xxyy_xxy_0, tg_xxyy_xxz_0, \
                                         tg_xxyy_xyy_0, tg_xxyy_xyz_0, tg_xxyy_xzz_0, tg_xxyy_yyy_0, tg_xxyy_yyz_0, \
                                         tg_xxyy_yzz_0, tg_xxyy_zzz_0, tg_xxyz_xxx_0, tg_xxyz_xxy_0, tg_xxyz_xxz_0, \
                                         tg_xxyz_xyy_0, tg_xxyz_xyz_0, tg_xxyz_xzz_0, tg_xxyz_yyy_0, tg_xxyz_yyz_0, \
                                         tg_xxyz_yzz_0, tg_xxyz_zzz_0, tg_xxz_xx_1, tg_xxz_xxx_0, tg_xxz_xxx_1, tg_xxz_xxy_0, \
                                         tg_xxz_xxy_1, tg_xxz_xxz_0, tg_xxz_xxz_1, tg_xxz_xy_1, tg_xxz_xyy_0, tg_xxz_xyy_1, \
                                         tg_xxz_xyz_0, tg_xxz_xyz_1, tg_xxz_xz_1, tg_xxz_xzz_0, tg_xxz_xzz_1, tg_xxz_yy_1, \
                                         tg_xxz_yyy_0, tg_xxz_yyy_1, tg_xxz_yyz_0, tg_xxz_yyz_1, tg_xxz_yz_1, tg_xxz_yzz_0, \
                                         tg_xxz_yzz_1, tg_xxz_zz_1, tg_xxz_zzz_0, tg_xxz_zzz_1, tg_xxzz_xxx_0, \
                                         tg_xxzz_xxy_0, tg_xxzz_xxz_0, tg_xxzz_xyy_0, tg_xxzz_xyz_0, tg_xxzz_xzz_0, \
                                         tg_xxzz_yyy_0, tg_xxzz_yyz_0, tg_xxzz_yzz_0, tg_xxzz_zzz_0, tg_xy_xxx_0, \
                                         tg_xy_xxx_1, tg_xy_xxy_0, tg_xy_xxy_1, tg_xy_xxz_0, tg_xy_xxz_1, tg_xy_xyy_0, \
                                         tg_xy_xyy_1, tg_xy_xyz_0, tg_xy_xyz_1, tg_xy_xzz_0, tg_xy_xzz_1, tg_xy_yyy_0, \
                                         tg_xy_yyy_1, tg_xy_yyz_0, tg_xy_yyz_1, tg_xy_yzz_0, tg_xy_yzz_1, tg_xy_zzz_0, \
                                         tg_xy_zzz_1, tg_xyy_xx_1, tg_xyy_xxx_0, tg_xyy_xxx_1, tg_xyy_xxy_0, tg_xyy_xxy_1, \
                                         tg_xyy_xxz_0, tg_xyy_xxz_1, tg_xyy_xy_1, tg_xyy_xyy_0, tg_xyy_xyy_1, tg_xyy_xyz_0, \
                                         tg_xyy_xyz_1, tg_xyy_xz_1, tg_xyy_xzz_0, tg_xyy_xzz_1, tg_xyy_yy_1, tg_xyy_yyy_0, \
                                         tg_xyy_yyy_1, tg_xyy_yyz_0, tg_xyy_yyz_1, tg_xyy_yz_1, tg_xyy_yzz_0, tg_xyy_yzz_1, \
                                         tg_xyy_zz_1, tg_xyy_zzz_0, tg_xyy_zzz_1, tg_xyyy_xxx_0, tg_xyyy_xxy_0, \
                                         tg_xyyy_xxz_0, tg_xyyy_xyy_0, tg_xyyy_xyz_0, tg_xyyy_xzz_0, tg_xyyy_yyy_0, \
                                         tg_xyyy_yyz_0, tg_xyyy_yzz_0, tg_xyyy_zzz_0, tg_xyyz_xxx_0, tg_xyyz_xxy_0, \
                                         tg_xyyz_xxz_0, tg_xyyz_xyy_0, tg_xyyz_xyz_0, tg_xyz_xx_1, tg_xyz_xxx_0, \
                                         tg_xyz_xxx_1, tg_xyz_xxy_0, tg_xyz_xxy_1, tg_xyz_xxz_0, tg_xyz_xxz_1, tg_xyz_xy_1, \
                                         tg_xyz_xyy_0, tg_xyz_xyy_1, tg_xyz_xyz_0, tg_xyz_xyz_1, tg_xyz_xz_1, tg_xyz_xzz_0, \
                                         tg_xyz_xzz_1, tg_xyz_yy_1, tg_xyz_yyy_0, tg_xyz_yyy_1, tg_xyz_yyz_0, tg_xyz_yyz_1, \
                                         tg_xyz_yz_1, tg_xyz_yzz_0, tg_xyz_yzz_1, tg_xyz_zz_1, tg_xyz_zzz_0, tg_xyz_zzz_1, \
                                         tg_xz_xxx_0, tg_xz_xxx_1, tg_xz_xxy_0, tg_xz_xxy_1, tg_xz_xxz_0, tg_xz_xxz_1, \
                                         tg_xz_xyy_0, tg_xz_xyy_1, tg_xz_xyz_0, tg_xz_xyz_1, tg_xz_xzz_0, tg_xz_xzz_1, \
                                         tg_xz_yyy_0, tg_xz_yyy_1, tg_xz_yyz_0, tg_xz_yyz_1, tg_xz_yzz_0, tg_xz_yzz_1, \
                                         tg_xz_zzz_0, tg_xz_zzz_1, tg_xzz_xx_1, tg_xzz_xxx_0, tg_xzz_xxx_1, tg_xzz_xxy_0, \
                                         tg_xzz_xxy_1, tg_xzz_xxz_0, tg_xzz_xxz_1, tg_xzz_xy_1, tg_xzz_xyy_0, tg_xzz_xyy_1, \
                                         tg_xzz_xyz_0, tg_xzz_xyz_1, tg_xzz_xz_1, tg_xzz_xzz_0, tg_xzz_xzz_1, tg_xzz_yy_1, \
                                         tg_xzz_yyy_0, tg_xzz_yyy_1, tg_xzz_yyz_0, tg_xzz_yyz_1, tg_xzz_yz_1, tg_xzz_yzz_0, \
                                         tg_xzz_yzz_1, tg_xzz_zz_1, tg_xzz_zzz_0, tg_xzz_zzz_1, tg_yy_xxx_0, tg_yy_xxx_1, \
                                         tg_yy_xxy_0, tg_yy_xxy_1, tg_yy_xxz_0, tg_yy_xxz_1, tg_yy_xyy_0, tg_yy_xyy_1, \
                                         tg_yy_xyz_0, tg_yy_xyz_1, tg_yy_xzz_0, tg_yy_xzz_1, tg_yy_yyy_0, tg_yy_yyy_1, \
                                         tg_yy_yyz_0, tg_yy_yyz_1, tg_yy_yzz_0, tg_yy_yzz_1, tg_yy_zzz_0, tg_yy_zzz_1, \
                                         tg_yyy_xx_1, tg_yyy_xxx_0, tg_yyy_xxx_1, tg_yyy_xxy_0, tg_yyy_xxy_1, tg_yyy_xxz_0, \
                                         tg_yyy_xxz_1, tg_yyy_xy_1, tg_yyy_xyy_0, tg_yyy_xyy_1, tg_yyy_xyz_0, tg_yyy_xyz_1, \
                                         tg_yyy_xz_1, tg_yyy_xzz_0, tg_yyy_xzz_1, tg_yyy_yy_1, tg_yyy_yyy_0, tg_yyy_yyy_1, \
                                         tg_yyy_yyz_0, tg_yyy_yyz_1, tg_yyy_yz_1, tg_yyy_yzz_0, tg_yyy_yzz_1, tg_yyy_zz_1, \
                                         tg_yyy_zzz_0, tg_yyy_zzz_1, tg_yyz_xx_1, tg_yyz_xxx_0, tg_yyz_xxx_1, tg_yyz_xxy_0, \
                                         tg_yyz_xxy_1, tg_yyz_xxz_0, tg_yyz_xxz_1, tg_yyz_xy_1, tg_yyz_xyy_0, tg_yyz_xyy_1, \
                                         tg_yyz_xyz_0, tg_yyz_xyz_1, tg_yyz_xz_1, tg_yyz_yy_1, tg_yyz_yz_1, tg_yz_xxx_0, \
                                         tg_yz_xxx_1, tg_yz_xxy_0, tg_yz_xxy_1, tg_yz_xxz_0, tg_yz_xxz_1, tg_yz_xyy_0, \
                                         tg_yz_xyy_1, tg_yz_xyz_0, tg_yz_xyz_1, tg_yz_xzz_0, tg_yz_xzz_1, tg_yz_yyy_0, \
                                         tg_yz_yyy_1, tg_yz_yyz_0, tg_yz_yyz_1, tg_yz_yzz_0, tg_yz_yzz_1, tg_yz_zzz_0, \
                                         tg_yz_zzz_1, tg_zz_xxx_0, tg_zz_xxx_1, tg_zz_xxy_0, tg_zz_xxy_1, tg_zz_xxz_0, \
                                         tg_zz_xxz_1, tg_zz_xyy_0, tg_zz_xyy_1, tg_zz_xyz_0, tg_zz_xyz_1, tg_zz_xzz_0, \
                                         tg_zz_xzz_1, tg_zz_yyy_0, tg_zz_yyy_1, tg_zz_yyz_0, tg_zz_yyz_1, tg_zz_yzz_0, \
                                         tg_zz_yzz_1, tg_zz_zzz_0, tg_zz_zzz_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xxxx_xxx_0[j] = pb_x * tg_xxx_xxx_0[j] + fr * tg_xxx_xxx_1[j] + 1.5 * fl1_fx * (tg_xx_xxx_0[j] - tg_xx_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxx_xx_1[j];

                    tg_xxxx_xxy_0[j] = pb_x * tg_xxx_xxy_0[j] + fr * tg_xxx_xxy_1[j] + 1.5 * fl1_fx * (tg_xx_xxy_0[j] - tg_xx_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xy_1[j];

                    tg_xxxx_xxz_0[j] = pb_x * tg_xxx_xxz_0[j] + fr * tg_xxx_xxz_1[j] + 1.5 * fl1_fx * (tg_xx_xxz_0[j] - tg_xx_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxx_xz_1[j];

                    tg_xxxx_xyy_0[j] = pb_x * tg_xxx_xyy_0[j] + fr * tg_xxx_xyy_1[j] + 1.5 * fl1_fx * (tg_xx_xyy_0[j] - tg_xx_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yy_1[j];

                    tg_xxxx_xyz_0[j] = pb_x * tg_xxx_xyz_0[j] + fr * tg_xxx_xyz_1[j] + 1.5 * fl1_fx * (tg_xx_xyz_0[j] - tg_xx_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_yz_1[j];

                    tg_xxxx_xzz_0[j] = pb_x * tg_xxx_xzz_0[j] + fr * tg_xxx_xzz_1[j] + 1.5 * fl1_fx * (tg_xx_xzz_0[j] - tg_xx_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxx_zz_1[j];

                    tg_xxxx_yyy_0[j] = pb_x * tg_xxx_yyy_0[j] + fr * tg_xxx_yyy_1[j] + 1.5 * fl1_fx * (tg_xx_yyy_0[j] - tg_xx_yyy_1[j] * fl1_fza);

                    tg_xxxx_yyz_0[j] = pb_x * tg_xxx_yyz_0[j] + fr * tg_xxx_yyz_1[j] + 1.5 * fl1_fx * (tg_xx_yyz_0[j] - tg_xx_yyz_1[j] * fl1_fza);

                    tg_xxxx_yzz_0[j] = pb_x * tg_xxx_yzz_0[j] + fr * tg_xxx_yzz_1[j] + 1.5 * fl1_fx * (tg_xx_yzz_0[j] - tg_xx_yzz_1[j] * fl1_fza);

                    tg_xxxx_zzz_0[j] = pb_x * tg_xxx_zzz_0[j] + fr * tg_xxx_zzz_1[j] + 1.5 * fl1_fx * (tg_xx_zzz_0[j] - tg_xx_zzz_1[j] * fl1_fza);

                    tg_xxxy_xxx_0[j] = pb_x * tg_xxy_xxx_0[j] + fr * tg_xxy_xxx_1[j] + fl1_fx * (tg_xy_xxx_0[j] - tg_xy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxy_xx_1[j];

                    tg_xxxy_xxy_0[j] = pb_x * tg_xxy_xxy_0[j] + fr * tg_xxy_xxy_1[j] + fl1_fx * (tg_xy_xxy_0[j] - tg_xy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xy_1[j];

                    tg_xxxy_xxz_0[j] = pb_x * tg_xxy_xxz_0[j] + fr * tg_xxy_xxz_1[j] + fl1_fx * (tg_xy_xxz_0[j] - tg_xy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxy_xz_1[j];

                    tg_xxxy_xyy_0[j] = pb_x * tg_xxy_xyy_0[j] + fr * tg_xxy_xyy_1[j] + fl1_fx * (tg_xy_xyy_0[j] - tg_xy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yy_1[j];

                    tg_xxxy_xyz_0[j] = pb_x * tg_xxy_xyz_0[j] + fr * tg_xxy_xyz_1[j] + fl1_fx * (tg_xy_xyz_0[j] - tg_xy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_yz_1[j];

                    tg_xxxy_xzz_0[j] = pb_x * tg_xxy_xzz_0[j] + fr * tg_xxy_xzz_1[j] + fl1_fx * (tg_xy_xzz_0[j] - tg_xy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxy_zz_1[j];

                    tg_xxxy_yyy_0[j] = pb_x * tg_xxy_yyy_0[j] + fr * tg_xxy_yyy_1[j] + fl1_fx * (tg_xy_yyy_0[j] - tg_xy_yyy_1[j] * fl1_fza);

                    tg_xxxy_yyz_0[j] = pb_x * tg_xxy_yyz_0[j] + fr * tg_xxy_yyz_1[j] + fl1_fx * (tg_xy_yyz_0[j] - tg_xy_yyz_1[j] * fl1_fza);

                    tg_xxxy_yzz_0[j] = pb_x * tg_xxy_yzz_0[j] + fr * tg_xxy_yzz_1[j] + fl1_fx * (tg_xy_yzz_0[j] - tg_xy_yzz_1[j] * fl1_fza);

                    tg_xxxy_zzz_0[j] = pb_x * tg_xxy_zzz_0[j] + fr * tg_xxy_zzz_1[j] + fl1_fx * (tg_xy_zzz_0[j] - tg_xy_zzz_1[j] * fl1_fza);

                    tg_xxxz_xxx_0[j] = pb_x * tg_xxz_xxx_0[j] + fr * tg_xxz_xxx_1[j] + fl1_fx * (tg_xz_xxx_0[j] - tg_xz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xxz_xx_1[j];

                    tg_xxxz_xxy_0[j] = pb_x * tg_xxz_xxy_0[j] + fr * tg_xxz_xxy_1[j] + fl1_fx * (tg_xz_xxy_0[j] - tg_xz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xy_1[j];

                    tg_xxxz_xxz_0[j] = pb_x * tg_xxz_xxz_0[j] + fr * tg_xxz_xxz_1[j] + fl1_fx * (tg_xz_xxz_0[j] - tg_xz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xxz_xz_1[j];

                    tg_xxxz_xyy_0[j] = pb_x * tg_xxz_xyy_0[j] + fr * tg_xxz_xyy_1[j] + fl1_fx * (tg_xz_xyy_0[j] - tg_xz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yy_1[j];

                    tg_xxxz_xyz_0[j] = pb_x * tg_xxz_xyz_0[j] + fr * tg_xxz_xyz_1[j] + fl1_fx * (tg_xz_xyz_0[j] - tg_xz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_yz_1[j];

                    tg_xxxz_xzz_0[j] = pb_x * tg_xxz_xzz_0[j] + fr * tg_xxz_xzz_1[j] + fl1_fx * (tg_xz_xzz_0[j] - tg_xz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xxz_zz_1[j];

                    tg_xxxz_yyy_0[j] = pb_x * tg_xxz_yyy_0[j] + fr * tg_xxz_yyy_1[j] + fl1_fx * (tg_xz_yyy_0[j] - tg_xz_yyy_1[j] * fl1_fza);

                    tg_xxxz_yyz_0[j] = pb_x * tg_xxz_yyz_0[j] + fr * tg_xxz_yyz_1[j] + fl1_fx * (tg_xz_yyz_0[j] - tg_xz_yyz_1[j] * fl1_fza);

                    tg_xxxz_yzz_0[j] = pb_x * tg_xxz_yzz_0[j] + fr * tg_xxz_yzz_1[j] + fl1_fx * (tg_xz_yzz_0[j] - tg_xz_yzz_1[j] * fl1_fza);

                    tg_xxxz_zzz_0[j] = pb_x * tg_xxz_zzz_0[j] + fr * tg_xxz_zzz_1[j] + fl1_fx * (tg_xz_zzz_0[j] - tg_xz_zzz_1[j] * fl1_fza);

                    tg_xxyy_xxx_0[j] = pb_x * tg_xyy_xxx_0[j] + fr * tg_xyy_xxx_1[j] + 0.5 * fl1_fx * (tg_yy_xxx_0[j] - tg_yy_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyy_xx_1[j];

                    tg_xxyy_xxy_0[j] = pb_x * tg_xyy_xxy_0[j] + fr * tg_xyy_xxy_1[j] + 0.5 * fl1_fx * (tg_yy_xxy_0[j] - tg_yy_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xy_1[j];

                    tg_xxyy_xxz_0[j] = pb_x * tg_xyy_xxz_0[j] + fr * tg_xyy_xxz_1[j] + 0.5 * fl1_fx * (tg_yy_xxz_0[j] - tg_yy_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyy_xz_1[j];

                    tg_xxyy_xyy_0[j] = pb_x * tg_xyy_xyy_0[j] + fr * tg_xyy_xyy_1[j] + 0.5 * fl1_fx * (tg_yy_xyy_0[j] - tg_yy_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yy_1[j];

                    tg_xxyy_xyz_0[j] = pb_x * tg_xyy_xyz_0[j] + fr * tg_xyy_xyz_1[j] + 0.5 * fl1_fx * (tg_yy_xyz_0[j] - tg_yy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_yz_1[j];

                    tg_xxyy_xzz_0[j] = pb_x * tg_xyy_xzz_0[j] + fr * tg_xyy_xzz_1[j] + 0.5 * fl1_fx * (tg_yy_xzz_0[j] - tg_yy_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyy_zz_1[j];

                    tg_xxyy_yyy_0[j] = pb_x * tg_xyy_yyy_0[j] + fr * tg_xyy_yyy_1[j] + 0.5 * fl1_fx * (tg_yy_yyy_0[j] - tg_yy_yyy_1[j] * fl1_fza);

                    tg_xxyy_yyz_0[j] = pb_x * tg_xyy_yyz_0[j] + fr * tg_xyy_yyz_1[j] + 0.5 * fl1_fx * (tg_yy_yyz_0[j] - tg_yy_yyz_1[j] * fl1_fza);

                    tg_xxyy_yzz_0[j] = pb_x * tg_xyy_yzz_0[j] + fr * tg_xyy_yzz_1[j] + 0.5 * fl1_fx * (tg_yy_yzz_0[j] - tg_yy_yzz_1[j] * fl1_fza);

                    tg_xxyy_zzz_0[j] = pb_x * tg_xyy_zzz_0[j] + fr * tg_xyy_zzz_1[j] + 0.5 * fl1_fx * (tg_yy_zzz_0[j] - tg_yy_zzz_1[j] * fl1_fza);

                    tg_xxyz_xxx_0[j] = pb_x * tg_xyz_xxx_0[j] + fr * tg_xyz_xxx_1[j] + 0.5 * fl1_fx * (tg_yz_xxx_0[j] - tg_yz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xyz_xx_1[j];

                    tg_xxyz_xxy_0[j] = pb_x * tg_xyz_xxy_0[j] + fr * tg_xyz_xxy_1[j] + 0.5 * fl1_fx * (tg_yz_xxy_0[j] - tg_yz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xy_1[j];

                    tg_xxyz_xxz_0[j] = pb_x * tg_xyz_xxz_0[j] + fr * tg_xyz_xxz_1[j] + 0.5 * fl1_fx * (tg_yz_xxz_0[j] - tg_yz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xyz_xz_1[j];

                    tg_xxyz_xyy_0[j] = pb_x * tg_xyz_xyy_0[j] + fr * tg_xyz_xyy_1[j] + 0.5 * fl1_fx * (tg_yz_xyy_0[j] - tg_yz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yy_1[j];

                    tg_xxyz_xyz_0[j] = pb_x * tg_xyz_xyz_0[j] + fr * tg_xyz_xyz_1[j] + 0.5 * fl1_fx * (tg_yz_xyz_0[j] - tg_yz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_yz_1[j];

                    tg_xxyz_xzz_0[j] = pb_x * tg_xyz_xzz_0[j] + fr * tg_xyz_xzz_1[j] + 0.5 * fl1_fx * (tg_yz_xzz_0[j] - tg_yz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xyz_zz_1[j];

                    tg_xxyz_yyy_0[j] = pb_x * tg_xyz_yyy_0[j] + fr * tg_xyz_yyy_1[j] + 0.5 * fl1_fx * (tg_yz_yyy_0[j] - tg_yz_yyy_1[j] * fl1_fza);

                    tg_xxyz_yyz_0[j] = pb_x * tg_xyz_yyz_0[j] + fr * tg_xyz_yyz_1[j] + 0.5 * fl1_fx * (tg_yz_yyz_0[j] - tg_yz_yyz_1[j] * fl1_fza);

                    tg_xxyz_yzz_0[j] = pb_x * tg_xyz_yzz_0[j] + fr * tg_xyz_yzz_1[j] + 0.5 * fl1_fx * (tg_yz_yzz_0[j] - tg_yz_yzz_1[j] * fl1_fza);

                    tg_xxyz_zzz_0[j] = pb_x * tg_xyz_zzz_0[j] + fr * tg_xyz_zzz_1[j] + 0.5 * fl1_fx * (tg_yz_zzz_0[j] - tg_yz_zzz_1[j] * fl1_fza);

                    tg_xxzz_xxx_0[j] = pb_x * tg_xzz_xxx_0[j] + fr * tg_xzz_xxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxx_0[j] - tg_zz_xxx_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_xzz_xx_1[j];

                    tg_xxzz_xxy_0[j] = pb_x * tg_xzz_xxy_0[j] + fr * tg_xzz_xxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxy_0[j] - tg_zz_xxy_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xy_1[j];

                    tg_xxzz_xxz_0[j] = pb_x * tg_xzz_xxz_0[j] + fr * tg_xzz_xxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxz_0[j] - tg_zz_xxz_1[j] * fl1_fza) + fl1_fxn * tg_xzz_xz_1[j];

                    tg_xxzz_xyy_0[j] = pb_x * tg_xzz_xyy_0[j] + fr * tg_xzz_xyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyy_0[j] - tg_zz_xyy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yy_1[j];

                    tg_xxzz_xyz_0[j] = pb_x * tg_xzz_xyz_0[j] + fr * tg_xzz_xyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyz_0[j] - tg_zz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_yz_1[j];

                    tg_xxzz_xzz_0[j] = pb_x * tg_xzz_xzz_0[j] + fr * tg_xzz_xzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzz_0[j] - tg_zz_xzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_xzz_zz_1[j];

                    tg_xxzz_yyy_0[j] = pb_x * tg_xzz_yyy_0[j] + fr * tg_xzz_yyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyy_0[j] - tg_zz_yyy_1[j] * fl1_fza);

                    tg_xxzz_yyz_0[j] = pb_x * tg_xzz_yyz_0[j] + fr * tg_xzz_yyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyz_0[j] - tg_zz_yyz_1[j] * fl1_fza);

                    tg_xxzz_yzz_0[j] = pb_x * tg_xzz_yzz_0[j] + fr * tg_xzz_yzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzz_0[j] - tg_zz_yzz_1[j] * fl1_fza);

                    tg_xxzz_zzz_0[j] = pb_x * tg_xzz_zzz_0[j] + fr * tg_xzz_zzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzz_0[j] - tg_zz_zzz_1[j] * fl1_fza);

                    tg_xyyy_xxx_0[j] = pb_x * tg_yyy_xxx_0[j] + fr * tg_yyy_xxx_1[j] + 1.5 * fl1_fxn * tg_yyy_xx_1[j];

                    tg_xyyy_xxy_0[j] = pb_x * tg_yyy_xxy_0[j] + fr * tg_yyy_xxy_1[j] + fl1_fxn * tg_yyy_xy_1[j];

                    tg_xyyy_xxz_0[j] = pb_x * tg_yyy_xxz_0[j] + fr * tg_yyy_xxz_1[j] + fl1_fxn * tg_yyy_xz_1[j];

                    tg_xyyy_xyy_0[j] = pb_x * tg_yyy_xyy_0[j] + fr * tg_yyy_xyy_1[j] + 0.5 * fl1_fxn * tg_yyy_yy_1[j];

                    tg_xyyy_xyz_0[j] = pb_x * tg_yyy_xyz_0[j] + fr * tg_yyy_xyz_1[j] + 0.5 * fl1_fxn * tg_yyy_yz_1[j];

                    tg_xyyy_xzz_0[j] = pb_x * tg_yyy_xzz_0[j] + fr * tg_yyy_xzz_1[j] + 0.5 * fl1_fxn * tg_yyy_zz_1[j];

                    tg_xyyy_yyy_0[j] = pb_x * tg_yyy_yyy_0[j] + fr * tg_yyy_yyy_1[j];

                    tg_xyyy_yyz_0[j] = pb_x * tg_yyy_yyz_0[j] + fr * tg_yyy_yyz_1[j];

                    tg_xyyy_yzz_0[j] = pb_x * tg_yyy_yzz_0[j] + fr * tg_yyy_yzz_1[j];

                    tg_xyyy_zzz_0[j] = pb_x * tg_yyy_zzz_0[j] + fr * tg_yyy_zzz_1[j];

                    tg_xyyz_xxx_0[j] = pb_x * tg_yyz_xxx_0[j] + fr * tg_yyz_xxx_1[j] + 1.5 * fl1_fxn * tg_yyz_xx_1[j];

                    tg_xyyz_xxy_0[j] = pb_x * tg_yyz_xxy_0[j] + fr * tg_yyz_xxy_1[j] + fl1_fxn * tg_yyz_xy_1[j];

                    tg_xyyz_xxz_0[j] = pb_x * tg_yyz_xxz_0[j] + fr * tg_yyz_xxz_1[j] + fl1_fxn * tg_yyz_xz_1[j];

                    tg_xyyz_xyy_0[j] = pb_x * tg_yyz_xyy_0[j] + fr * tg_yyz_xyy_1[j] + 0.5 * fl1_fxn * tg_yyz_yy_1[j];

                    tg_xyyz_xyz_0[j] = pb_x * tg_yyz_xyz_0[j] + fr * tg_yyz_xyz_1[j] + 0.5 * fl1_fxn * tg_yyz_yz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSF_75_150(      CMemBlock2D<double>* primBuffer,
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

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto b_fx = braGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {4, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_3_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_3_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                                   1, 1, iord + 1));

            auto pidx_g_2_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                                   1, 1, iord));

            auto pidx_g_2_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                                   {2, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tg_yyy_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 60); 

                auto tg_yyy_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 61); 

                auto tg_yyy_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 62); 

                auto tg_yyy_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 63); 

                auto tg_yyy_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 64); 

                auto tg_yyy_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 65); 

                auto tg_yyy_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 66); 

                auto tg_yyy_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 67); 

                auto tg_yyy_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 68); 

                auto tg_yyy_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 69); 

                auto tg_yyz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 70); 

                auto tg_yyz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 71); 

                auto tg_yyz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 72); 

                auto tg_yyz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 73); 

                auto tg_yyz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 74); 

                auto tg_yyz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 75); 

                auto tg_yyz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 76); 

                auto tg_yyz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 77); 

                auto tg_yyz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 78); 

                auto tg_yyz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 79); 

                auto tg_yzz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 80); 

                auto tg_yzz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 81); 

                auto tg_yzz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 82); 

                auto tg_yzz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 83); 

                auto tg_yzz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 84); 

                auto tg_yzz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 85); 

                auto tg_yzz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 86); 

                auto tg_yzz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 87); 

                auto tg_yzz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 88); 

                auto tg_yzz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 89); 

                auto tg_zzz_xxx_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 90); 

                auto tg_zzz_xxy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 91); 

                auto tg_zzz_xxz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 92); 

                auto tg_zzz_xyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 93); 

                auto tg_zzz_xyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 94); 

                auto tg_zzz_xzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 95); 

                auto tg_zzz_yyy_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 96); 

                auto tg_zzz_yyz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 97); 

                auto tg_zzz_yzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 98); 

                auto tg_zzz_zzz_0 = primBuffer[pidx_g_3_3_m0].data(100 * idx + 99); 

                auto tg_yyy_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 60); 

                auto tg_yyy_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 61); 

                auto tg_yyy_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 62); 

                auto tg_yyy_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 63); 

                auto tg_yyy_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 64); 

                auto tg_yyy_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 65); 

                auto tg_yyy_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 66); 

                auto tg_yyy_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 67); 

                auto tg_yyy_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 68); 

                auto tg_yyy_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 69); 

                auto tg_yyz_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 70); 

                auto tg_yyz_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 71); 

                auto tg_yyz_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 72); 

                auto tg_yyz_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 73); 

                auto tg_yyz_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 74); 

                auto tg_yyz_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 75); 

                auto tg_yyz_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 76); 

                auto tg_yyz_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 77); 

                auto tg_yyz_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 78); 

                auto tg_yyz_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 79); 

                auto tg_yzz_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 80); 

                auto tg_yzz_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 81); 

                auto tg_yzz_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 82); 

                auto tg_yzz_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 83); 

                auto tg_yzz_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 84); 

                auto tg_yzz_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 85); 

                auto tg_yzz_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 86); 

                auto tg_yzz_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 87); 

                auto tg_yzz_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 88); 

                auto tg_yzz_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 89); 

                auto tg_zzz_xxx_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 90); 

                auto tg_zzz_xxy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 91); 

                auto tg_zzz_xxz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 92); 

                auto tg_zzz_xyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 93); 

                auto tg_zzz_xyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 94); 

                auto tg_zzz_xzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 95); 

                auto tg_zzz_yyy_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 96); 

                auto tg_zzz_yyz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 97); 

                auto tg_zzz_yzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 98); 

                auto tg_zzz_zzz_1 = primBuffer[pidx_g_3_3_m1].data(100 * idx + 99); 

                auto tg_yy_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 30); 

                auto tg_yy_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 31); 

                auto tg_yy_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 32); 

                auto tg_yy_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 33); 

                auto tg_yy_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 34); 

                auto tg_yy_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 35); 

                auto tg_yy_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 36); 

                auto tg_yy_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 37); 

                auto tg_yy_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 38); 

                auto tg_yy_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 39); 

                auto tg_yz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 40); 

                auto tg_yz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 41); 

                auto tg_yz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 42); 

                auto tg_yz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 43); 

                auto tg_yz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 44); 

                auto tg_yz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 45); 

                auto tg_yz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 46); 

                auto tg_yz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 47); 

                auto tg_yz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 48); 

                auto tg_yz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 49); 

                auto tg_zz_xxx_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 50); 

                auto tg_zz_xxy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 51); 

                auto tg_zz_xxz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 52); 

                auto tg_zz_xyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 53); 

                auto tg_zz_xyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 54); 

                auto tg_zz_xzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 55); 

                auto tg_zz_yyy_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 56); 

                auto tg_zz_yyz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 57); 

                auto tg_zz_yzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 58); 

                auto tg_zz_zzz_0 = primBuffer[pidx_g_2_3_m0].data(60 * idx + 59); 

                auto tg_yy_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 30); 

                auto tg_yy_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 31); 

                auto tg_yy_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 32); 

                auto tg_yy_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 33); 

                auto tg_yy_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 34); 

                auto tg_yy_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 35); 

                auto tg_yy_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 36); 

                auto tg_yy_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 37); 

                auto tg_yy_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 38); 

                auto tg_yy_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 39); 

                auto tg_yz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 40); 

                auto tg_yz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 41); 

                auto tg_yz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 42); 

                auto tg_yz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 43); 

                auto tg_yz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 44); 

                auto tg_yz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 45); 

                auto tg_yz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 46); 

                auto tg_yz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 47); 

                auto tg_yz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 48); 

                auto tg_yz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 49); 

                auto tg_zz_xxx_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 50); 

                auto tg_zz_xxy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 51); 

                auto tg_zz_xxz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 52); 

                auto tg_zz_xyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 53); 

                auto tg_zz_xyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 54); 

                auto tg_zz_xzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 55); 

                auto tg_zz_yyy_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 56); 

                auto tg_zz_yyz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 57); 

                auto tg_zz_yzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 58); 

                auto tg_zz_zzz_1 = primBuffer[pidx_g_2_3_m1].data(60 * idx + 59); 

                auto tg_yyy_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 36); 

                auto tg_yyy_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 37); 

                auto tg_yyy_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 38); 

                auto tg_yyy_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 39); 

                auto tg_yyy_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 40); 

                auto tg_yyy_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 41); 

                auto tg_yyz_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 42); 

                auto tg_yyz_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 43); 

                auto tg_yyz_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 44); 

                auto tg_yyz_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 45); 

                auto tg_yyz_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 46); 

                auto tg_yyz_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 47); 

                auto tg_yzz_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 48); 

                auto tg_yzz_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 49); 

                auto tg_yzz_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 50); 

                auto tg_yzz_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 51); 

                auto tg_yzz_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 52); 

                auto tg_yzz_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 53); 

                auto tg_zzz_xx_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 54); 

                auto tg_zzz_xy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 55); 

                auto tg_zzz_xz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 56); 

                auto tg_zzz_yy_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 57); 

                auto tg_zzz_yz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 58); 

                auto tg_zzz_zz_1 = primBuffer[pidx_g_3_2_m1].data(60 * idx + 59); 

                // set up pointers to integrals

                auto tg_xyyz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 75); 

                auto tg_xyyz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 76); 

                auto tg_xyyz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 77); 

                auto tg_xyyz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 78); 

                auto tg_xyyz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 79); 

                auto tg_xyzz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 80); 

                auto tg_xyzz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 81); 

                auto tg_xyzz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 82); 

                auto tg_xyzz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 83); 

                auto tg_xyzz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 84); 

                auto tg_xyzz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 85); 

                auto tg_xyzz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 86); 

                auto tg_xyzz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 87); 

                auto tg_xyzz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 88); 

                auto tg_xyzz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 89); 

                auto tg_xzzz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 90); 

                auto tg_xzzz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 91); 

                auto tg_xzzz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 92); 

                auto tg_xzzz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 93); 

                auto tg_xzzz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 94); 

                auto tg_xzzz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 95); 

                auto tg_xzzz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 96); 

                auto tg_xzzz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 97); 

                auto tg_xzzz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 98); 

                auto tg_xzzz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 99); 

                auto tg_yyyy_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 100); 

                auto tg_yyyy_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 101); 

                auto tg_yyyy_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 102); 

                auto tg_yyyy_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 103); 

                auto tg_yyyy_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 104); 

                auto tg_yyyy_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 105); 

                auto tg_yyyy_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 106); 

                auto tg_yyyy_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 107); 

                auto tg_yyyy_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 108); 

                auto tg_yyyy_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 109); 

                auto tg_yyyz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 110); 

                auto tg_yyyz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 111); 

                auto tg_yyyz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 112); 

                auto tg_yyyz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 113); 

                auto tg_yyyz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 114); 

                auto tg_yyyz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 115); 

                auto tg_yyyz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 116); 

                auto tg_yyyz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 117); 

                auto tg_yyyz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 118); 

                auto tg_yyyz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 119); 

                auto tg_yyzz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 120); 

                auto tg_yyzz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 121); 

                auto tg_yyzz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 122); 

                auto tg_yyzz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 123); 

                auto tg_yyzz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 124); 

                auto tg_yyzz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 125); 

                auto tg_yyzz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 126); 

                auto tg_yyzz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 127); 

                auto tg_yyzz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 128); 

                auto tg_yyzz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 129); 

                auto tg_yzzz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 130); 

                auto tg_yzzz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 131); 

                auto tg_yzzz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 132); 

                auto tg_yzzz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 133); 

                auto tg_yzzz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 134); 

                auto tg_yzzz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 135); 

                auto tg_yzzz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 136); 

                auto tg_yzzz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 137); 

                auto tg_yzzz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 138); 

                auto tg_yzzz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 139); 

                auto tg_zzzz_xxx_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 140); 

                auto tg_zzzz_xxy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 141); 

                auto tg_zzzz_xxz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 142); 

                auto tg_zzzz_xyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 143); 

                auto tg_zzzz_xyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 144); 

                auto tg_zzzz_xzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 145); 

                auto tg_zzzz_yyy_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 146); 

                auto tg_zzzz_yyz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 147); 

                auto tg_zzzz_yzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 148); 

                auto tg_zzzz_zzz_0 = primBuffer[pidx_g_4_3_m0].data(150 * idx + 149); 

                // Batch of Integrals (75,150)

                #pragma omp simd aligned(fxn, fza, tg_xyyz_xzz_0, tg_xyyz_yyy_0, tg_xyyz_yyz_0, tg_xyyz_yzz_0, \
                                         tg_xyyz_zzz_0, tg_xyzz_xxx_0, tg_xyzz_xxy_0, tg_xyzz_xxz_0, tg_xyzz_xyy_0, \
                                         tg_xyzz_xyz_0, tg_xyzz_xzz_0, tg_xyzz_yyy_0, tg_xyzz_yyz_0, tg_xyzz_yzz_0, \
                                         tg_xyzz_zzz_0, tg_xzzz_xxx_0, tg_xzzz_xxy_0, tg_xzzz_xxz_0, tg_xzzz_xyy_0, \
                                         tg_xzzz_xyz_0, tg_xzzz_xzz_0, tg_xzzz_yyy_0, tg_xzzz_yyz_0, tg_xzzz_yzz_0, \
                                         tg_xzzz_zzz_0, tg_yy_xxx_0, tg_yy_xxx_1, tg_yy_xxy_0, tg_yy_xxy_1, tg_yy_xxz_0, \
                                         tg_yy_xxz_1, tg_yy_xyy_0, tg_yy_xyy_1, tg_yy_xyz_0, tg_yy_xyz_1, tg_yy_xzz_0, \
                                         tg_yy_xzz_1, tg_yy_yyy_0, tg_yy_yyy_1, tg_yy_yyz_0, tg_yy_yyz_1, tg_yy_yzz_0, \
                                         tg_yy_yzz_1, tg_yy_zzz_0, tg_yy_zzz_1, tg_yyy_xx_1, tg_yyy_xxx_0, tg_yyy_xxx_1, \
                                         tg_yyy_xxy_0, tg_yyy_xxy_1, tg_yyy_xxz_0, tg_yyy_xxz_1, tg_yyy_xy_1, tg_yyy_xyy_0, \
                                         tg_yyy_xyy_1, tg_yyy_xyz_0, tg_yyy_xyz_1, tg_yyy_xz_1, tg_yyy_xzz_0, tg_yyy_xzz_1, \
                                         tg_yyy_yy_1, tg_yyy_yyy_0, tg_yyy_yyy_1, tg_yyy_yyz_0, tg_yyy_yyz_1, tg_yyy_yz_1, \
                                         tg_yyy_yzz_0, tg_yyy_yzz_1, tg_yyy_zz_1, tg_yyy_zzz_0, tg_yyy_zzz_1, tg_yyyy_xxx_0, \
                                         tg_yyyy_xxy_0, tg_yyyy_xxz_0, tg_yyyy_xyy_0, tg_yyyy_xyz_0, tg_yyyy_xzz_0, \
                                         tg_yyyy_yyy_0, tg_yyyy_yyz_0, tg_yyyy_yzz_0, tg_yyyy_zzz_0, tg_yyyz_xxx_0, \
                                         tg_yyyz_xxy_0, tg_yyyz_xxz_0, tg_yyyz_xyy_0, tg_yyyz_xyz_0, tg_yyyz_xzz_0, \
                                         tg_yyyz_yyy_0, tg_yyyz_yyz_0, tg_yyyz_yzz_0, tg_yyyz_zzz_0, tg_yyz_xx_1, \
                                         tg_yyz_xxx_0, tg_yyz_xxx_1, tg_yyz_xxy_0, tg_yyz_xxy_1, tg_yyz_xxz_0, tg_yyz_xxz_1, \
                                         tg_yyz_xy_1, tg_yyz_xyy_0, tg_yyz_xyy_1, tg_yyz_xyz_0, tg_yyz_xyz_1, tg_yyz_xz_1, \
                                         tg_yyz_xzz_0, tg_yyz_xzz_1, tg_yyz_yy_1, tg_yyz_yyy_0, tg_yyz_yyy_1, tg_yyz_yyz_0, \
                                         tg_yyz_yyz_1, tg_yyz_yz_1, tg_yyz_yzz_0, tg_yyz_yzz_1, tg_yyz_zz_1, tg_yyz_zzz_0, \
                                         tg_yyz_zzz_1, tg_yyzz_xxx_0, tg_yyzz_xxy_0, tg_yyzz_xxz_0, tg_yyzz_xyy_0, \
                                         tg_yyzz_xyz_0, tg_yyzz_xzz_0, tg_yyzz_yyy_0, tg_yyzz_yyz_0, tg_yyzz_yzz_0, \
                                         tg_yyzz_zzz_0, tg_yz_xxx_0, tg_yz_xxx_1, tg_yz_xxy_0, tg_yz_xxy_1, tg_yz_xxz_0, \
                                         tg_yz_xxz_1, tg_yz_xyy_0, tg_yz_xyy_1, tg_yz_xyz_0, tg_yz_xyz_1, tg_yz_xzz_0, \
                                         tg_yz_xzz_1, tg_yz_yyy_0, tg_yz_yyy_1, tg_yz_yyz_0, tg_yz_yyz_1, tg_yz_yzz_0, \
                                         tg_yz_yzz_1, tg_yz_zzz_0, tg_yz_zzz_1, tg_yzz_xx_1, tg_yzz_xxx_0, tg_yzz_xxx_1, \
                                         tg_yzz_xxy_0, tg_yzz_xxy_1, tg_yzz_xxz_0, tg_yzz_xxz_1, tg_yzz_xy_1, tg_yzz_xyy_0, \
                                         tg_yzz_xyy_1, tg_yzz_xyz_0, tg_yzz_xyz_1, tg_yzz_xz_1, tg_yzz_xzz_0, tg_yzz_xzz_1, \
                                         tg_yzz_yy_1, tg_yzz_yyy_0, tg_yzz_yyy_1, tg_yzz_yyz_0, tg_yzz_yyz_1, tg_yzz_yz_1, \
                                         tg_yzz_yzz_0, tg_yzz_yzz_1, tg_yzz_zz_1, tg_yzz_zzz_0, tg_yzz_zzz_1, tg_yzzz_xxx_0, \
                                         tg_yzzz_xxy_0, tg_yzzz_xxz_0, tg_yzzz_xyy_0, tg_yzzz_xyz_0, tg_yzzz_xzz_0, \
                                         tg_yzzz_yyy_0, tg_yzzz_yyz_0, tg_yzzz_yzz_0, tg_yzzz_zzz_0, tg_zz_xxx_0, \
                                         tg_zz_xxx_1, tg_zz_xxy_0, tg_zz_xxy_1, tg_zz_xxz_0, tg_zz_xxz_1, tg_zz_xyy_0, \
                                         tg_zz_xyy_1, tg_zz_xyz_0, tg_zz_xyz_1, tg_zz_xzz_0, tg_zz_xzz_1, tg_zz_yyy_0, \
                                         tg_zz_yyy_1, tg_zz_yyz_0, tg_zz_yyz_1, tg_zz_yzz_0, tg_zz_yzz_1, tg_zz_zzz_0, \
                                         tg_zz_zzz_1, tg_zzz_xx_1, tg_zzz_xxx_0, tg_zzz_xxx_1, tg_zzz_xxy_0, tg_zzz_xxy_1, \
                                         tg_zzz_xxz_0, tg_zzz_xxz_1, tg_zzz_xy_1, tg_zzz_xyy_0, tg_zzz_xyy_1, tg_zzz_xyz_0, \
                                         tg_zzz_xyz_1, tg_zzz_xz_1, tg_zzz_xzz_0, tg_zzz_xzz_1, tg_zzz_yy_1, tg_zzz_yyy_0, \
                                         tg_zzz_yyy_1, tg_zzz_yyz_0, tg_zzz_yyz_1, tg_zzz_yz_1, tg_zzz_yzz_0, tg_zzz_yzz_1, \
                                         tg_zzz_zz_1, tg_zzz_zzz_0, tg_zzz_zzz_1, tg_zzzz_xxx_0, tg_zzzz_xxy_0, \
                                         tg_zzzz_xxz_0, tg_zzzz_xyy_0, tg_zzzz_xyz_0, tg_zzzz_xzz_0, tg_zzzz_yyy_0, \
                                         tg_zzzz_yyz_0, tg_zzzz_yzz_0, tg_zzzz_zzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    double fr = wp_x[j]; 

                    tg_xyyz_xzz_0[j] = pb_x * tg_yyz_xzz_0[j] + fr * tg_yyz_xzz_1[j] + 0.5 * fl1_fxn * tg_yyz_zz_1[j];

                    tg_xyyz_yyy_0[j] = pb_x * tg_yyz_yyy_0[j] + fr * tg_yyz_yyy_1[j];

                    tg_xyyz_yyz_0[j] = pb_x * tg_yyz_yyz_0[j] + fr * tg_yyz_yyz_1[j];

                    tg_xyyz_yzz_0[j] = pb_x * tg_yyz_yzz_0[j] + fr * tg_yyz_yzz_1[j];

                    tg_xyyz_zzz_0[j] = pb_x * tg_yyz_zzz_0[j] + fr * tg_yyz_zzz_1[j];

                    tg_xyzz_xxx_0[j] = pb_x * tg_yzz_xxx_0[j] + fr * tg_yzz_xxx_1[j] + 1.5 * fl1_fxn * tg_yzz_xx_1[j];

                    tg_xyzz_xxy_0[j] = pb_x * tg_yzz_xxy_0[j] + fr * tg_yzz_xxy_1[j] + fl1_fxn * tg_yzz_xy_1[j];

                    tg_xyzz_xxz_0[j] = pb_x * tg_yzz_xxz_0[j] + fr * tg_yzz_xxz_1[j] + fl1_fxn * tg_yzz_xz_1[j];

                    tg_xyzz_xyy_0[j] = pb_x * tg_yzz_xyy_0[j] + fr * tg_yzz_xyy_1[j] + 0.5 * fl1_fxn * tg_yzz_yy_1[j];

                    tg_xyzz_xyz_0[j] = pb_x * tg_yzz_xyz_0[j] + fr * tg_yzz_xyz_1[j] + 0.5 * fl1_fxn * tg_yzz_yz_1[j];

                    tg_xyzz_xzz_0[j] = pb_x * tg_yzz_xzz_0[j] + fr * tg_yzz_xzz_1[j] + 0.5 * fl1_fxn * tg_yzz_zz_1[j];

                    tg_xyzz_yyy_0[j] = pb_x * tg_yzz_yyy_0[j] + fr * tg_yzz_yyy_1[j];

                    tg_xyzz_yyz_0[j] = pb_x * tg_yzz_yyz_0[j] + fr * tg_yzz_yyz_1[j];

                    tg_xyzz_yzz_0[j] = pb_x * tg_yzz_yzz_0[j] + fr * tg_yzz_yzz_1[j];

                    tg_xyzz_zzz_0[j] = pb_x * tg_yzz_zzz_0[j] + fr * tg_yzz_zzz_1[j];

                    tg_xzzz_xxx_0[j] = pb_x * tg_zzz_xxx_0[j] + fr * tg_zzz_xxx_1[j] + 1.5 * fl1_fxn * tg_zzz_xx_1[j];

                    tg_xzzz_xxy_0[j] = pb_x * tg_zzz_xxy_0[j] + fr * tg_zzz_xxy_1[j] + fl1_fxn * tg_zzz_xy_1[j];

                    tg_xzzz_xxz_0[j] = pb_x * tg_zzz_xxz_0[j] + fr * tg_zzz_xxz_1[j] + fl1_fxn * tg_zzz_xz_1[j];

                    tg_xzzz_xyy_0[j] = pb_x * tg_zzz_xyy_0[j] + fr * tg_zzz_xyy_1[j] + 0.5 * fl1_fxn * tg_zzz_yy_1[j];

                    tg_xzzz_xyz_0[j] = pb_x * tg_zzz_xyz_0[j] + fr * tg_zzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzz_yz_1[j];

                    tg_xzzz_xzz_0[j] = pb_x * tg_zzz_xzz_0[j] + fr * tg_zzz_xzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zz_1[j];

                    tg_xzzz_yyy_0[j] = pb_x * tg_zzz_yyy_0[j] + fr * tg_zzz_yyy_1[j];

                    tg_xzzz_yyz_0[j] = pb_x * tg_zzz_yyz_0[j] + fr * tg_zzz_yyz_1[j];

                    tg_xzzz_yzz_0[j] = pb_x * tg_zzz_yzz_0[j] + fr * tg_zzz_yzz_1[j];

                    tg_xzzz_zzz_0[j] = pb_x * tg_zzz_zzz_0[j] + fr * tg_zzz_zzz_1[j];

                    fr = wp_y[j]; 

                    tg_yyyy_xxx_0[j] = pb_y * tg_yyy_xxx_0[j] + fr * tg_yyy_xxx_1[j] + 1.5 * fl1_fx * (tg_yy_xxx_0[j] - tg_yy_xxx_1[j] * fl1_fza);

                    tg_yyyy_xxy_0[j] = pb_y * tg_yyy_xxy_0[j] + fr * tg_yyy_xxy_1[j] + 1.5 * fl1_fx * (tg_yy_xxy_0[j] - tg_yy_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xx_1[j];

                    tg_yyyy_xxz_0[j] = pb_y * tg_yyy_xxz_0[j] + fr * tg_yyy_xxz_1[j] + 1.5 * fl1_fx * (tg_yy_xxz_0[j] - tg_yy_xxz_1[j] * fl1_fza);

                    tg_yyyy_xyy_0[j] = pb_y * tg_yyy_xyy_0[j] + fr * tg_yyy_xyy_1[j] + 1.5 * fl1_fx * (tg_yy_xyy_0[j] - tg_yy_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyy_xy_1[j];

                    tg_yyyy_xyz_0[j] = pb_y * tg_yyy_xyz_0[j] + fr * tg_yyy_xyz_1[j] + 1.5 * fl1_fx * (tg_yy_xyz_0[j] - tg_yy_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_xz_1[j];

                    tg_yyyy_xzz_0[j] = pb_y * tg_yyy_xzz_0[j] + fr * tg_yyy_xzz_1[j] + 1.5 * fl1_fx * (tg_yy_xzz_0[j] - tg_yy_xzz_1[j] * fl1_fza);

                    tg_yyyy_yyy_0[j] = pb_y * tg_yyy_yyy_0[j] + fr * tg_yyy_yyy_1[j] + 1.5 * fl1_fx * (tg_yy_yyy_0[j] - tg_yy_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyy_yy_1[j];

                    tg_yyyy_yyz_0[j] = pb_y * tg_yyy_yyz_0[j] + fr * tg_yyy_yyz_1[j] + 1.5 * fl1_fx * (tg_yy_yyz_0[j] - tg_yy_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyy_yz_1[j];

                    tg_yyyy_yzz_0[j] = pb_y * tg_yyy_yzz_0[j] + fr * tg_yyy_yzz_1[j] + 1.5 * fl1_fx * (tg_yy_yzz_0[j] - tg_yy_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyy_zz_1[j];

                    tg_yyyy_zzz_0[j] = pb_y * tg_yyy_zzz_0[j] + fr * tg_yyy_zzz_1[j] + 1.5 * fl1_fx * (tg_yy_zzz_0[j] - tg_yy_zzz_1[j] * fl1_fza);

                    tg_yyyz_xxx_0[j] = pb_y * tg_yyz_xxx_0[j] + fr * tg_yyz_xxx_1[j] + fl1_fx * (tg_yz_xxx_0[j] - tg_yz_xxx_1[j] * fl1_fza);

                    tg_yyyz_xxy_0[j] = pb_y * tg_yyz_xxy_0[j] + fr * tg_yyz_xxy_1[j] + fl1_fx * (tg_yz_xxy_0[j] - tg_yz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xx_1[j];

                    tg_yyyz_xxz_0[j] = pb_y * tg_yyz_xxz_0[j] + fr * tg_yyz_xxz_1[j] + fl1_fx * (tg_yz_xxz_0[j] - tg_yz_xxz_1[j] * fl1_fza);

                    tg_yyyz_xyy_0[j] = pb_y * tg_yyz_xyy_0[j] + fr * tg_yyz_xyy_1[j] + fl1_fx * (tg_yz_xyy_0[j] - tg_yz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yyz_xy_1[j];

                    tg_yyyz_xyz_0[j] = pb_y * tg_yyz_xyz_0[j] + fr * tg_yyz_xyz_1[j] + fl1_fx * (tg_yz_xyz_0[j] - tg_yz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_xz_1[j];

                    tg_yyyz_xzz_0[j] = pb_y * tg_yyz_xzz_0[j] + fr * tg_yyz_xzz_1[j] + fl1_fx * (tg_yz_xzz_0[j] - tg_yz_xzz_1[j] * fl1_fza);

                    tg_yyyz_yyy_0[j] = pb_y * tg_yyz_yyy_0[j] + fr * tg_yyz_yyy_1[j] + fl1_fx * (tg_yz_yyy_0[j] - tg_yz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yyz_yy_1[j];

                    tg_yyyz_yyz_0[j] = pb_y * tg_yyz_yyz_0[j] + fr * tg_yyz_yyz_1[j] + fl1_fx * (tg_yz_yyz_0[j] - tg_yz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yyz_yz_1[j];

                    tg_yyyz_yzz_0[j] = pb_y * tg_yyz_yzz_0[j] + fr * tg_yyz_yzz_1[j] + fl1_fx * (tg_yz_yzz_0[j] - tg_yz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yyz_zz_1[j];

                    tg_yyyz_zzz_0[j] = pb_y * tg_yyz_zzz_0[j] + fr * tg_yyz_zzz_1[j] + fl1_fx * (tg_yz_zzz_0[j] - tg_yz_zzz_1[j] * fl1_fza);

                    tg_yyzz_xxx_0[j] = pb_y * tg_yzz_xxx_0[j] + fr * tg_yzz_xxx_1[j] + 0.5 * fl1_fx * (tg_zz_xxx_0[j] - tg_zz_xxx_1[j] * fl1_fza);

                    tg_yyzz_xxy_0[j] = pb_y * tg_yzz_xxy_0[j] + fr * tg_yzz_xxy_1[j] + 0.5 * fl1_fx * (tg_zz_xxy_0[j] - tg_zz_xxy_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xx_1[j];

                    tg_yyzz_xxz_0[j] = pb_y * tg_yzz_xxz_0[j] + fr * tg_yzz_xxz_1[j] + 0.5 * fl1_fx * (tg_zz_xxz_0[j] - tg_zz_xxz_1[j] * fl1_fza);

                    tg_yyzz_xyy_0[j] = pb_y * tg_yzz_xyy_0[j] + fr * tg_yzz_xyy_1[j] + 0.5 * fl1_fx * (tg_zz_xyy_0[j] - tg_zz_xyy_1[j] * fl1_fza) + fl1_fxn * tg_yzz_xy_1[j];

                    tg_yyzz_xyz_0[j] = pb_y * tg_yzz_xyz_0[j] + fr * tg_yzz_xyz_1[j] + 0.5 * fl1_fx * (tg_zz_xyz_0[j] - tg_zz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_xz_1[j];

                    tg_yyzz_xzz_0[j] = pb_y * tg_yzz_xzz_0[j] + fr * tg_yzz_xzz_1[j] + 0.5 * fl1_fx * (tg_zz_xzz_0[j] - tg_zz_xzz_1[j] * fl1_fza);

                    tg_yyzz_yyy_0[j] = pb_y * tg_yzz_yyy_0[j] + fr * tg_yzz_yyy_1[j] + 0.5 * fl1_fx * (tg_zz_yyy_0[j] - tg_zz_yyy_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_yzz_yy_1[j];

                    tg_yyzz_yyz_0[j] = pb_y * tg_yzz_yyz_0[j] + fr * tg_yzz_yyz_1[j] + 0.5 * fl1_fx * (tg_zz_yyz_0[j] - tg_zz_yyz_1[j] * fl1_fza) + fl1_fxn * tg_yzz_yz_1[j];

                    tg_yyzz_yzz_0[j] = pb_y * tg_yzz_yzz_0[j] + fr * tg_yzz_yzz_1[j] + 0.5 * fl1_fx * (tg_zz_yzz_0[j] - tg_zz_yzz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_yzz_zz_1[j];

                    tg_yyzz_zzz_0[j] = pb_y * tg_yzz_zzz_0[j] + fr * tg_yzz_zzz_1[j] + 0.5 * fl1_fx * (tg_zz_zzz_0[j] - tg_zz_zzz_1[j] * fl1_fza);

                    tg_yzzz_xxx_0[j] = pb_y * tg_zzz_xxx_0[j] + fr * tg_zzz_xxx_1[j];

                    tg_yzzz_xxy_0[j] = pb_y * tg_zzz_xxy_0[j] + fr * tg_zzz_xxy_1[j] + 0.5 * fl1_fxn * tg_zzz_xx_1[j];

                    tg_yzzz_xxz_0[j] = pb_y * tg_zzz_xxz_0[j] + fr * tg_zzz_xxz_1[j];

                    tg_yzzz_xyy_0[j] = pb_y * tg_zzz_xyy_0[j] + fr * tg_zzz_xyy_1[j] + fl1_fxn * tg_zzz_xy_1[j];

                    tg_yzzz_xyz_0[j] = pb_y * tg_zzz_xyz_0[j] + fr * tg_zzz_xyz_1[j] + 0.5 * fl1_fxn * tg_zzz_xz_1[j];

                    tg_yzzz_xzz_0[j] = pb_y * tg_zzz_xzz_0[j] + fr * tg_zzz_xzz_1[j];

                    tg_yzzz_yyy_0[j] = pb_y * tg_zzz_yyy_0[j] + fr * tg_zzz_yyy_1[j] + 1.5 * fl1_fxn * tg_zzz_yy_1[j];

                    tg_yzzz_yyz_0[j] = pb_y * tg_zzz_yyz_0[j] + fr * tg_zzz_yyz_1[j] + fl1_fxn * tg_zzz_yz_1[j];

                    tg_yzzz_yzz_0[j] = pb_y * tg_zzz_yzz_0[j] + fr * tg_zzz_yzz_1[j] + 0.5 * fl1_fxn * tg_zzz_zz_1[j];

                    tg_yzzz_zzz_0[j] = pb_y * tg_zzz_zzz_0[j] + fr * tg_zzz_zzz_1[j];

                    fr = wp_z[j]; 

                    tg_zzzz_xxx_0[j] = pb_z * tg_zzz_xxx_0[j] + fr * tg_zzz_xxx_1[j] + 1.5 * fl1_fx * (tg_zz_xxx_0[j] - tg_zz_xxx_1[j] * fl1_fza);

                    tg_zzzz_xxy_0[j] = pb_z * tg_zzz_xxy_0[j] + fr * tg_zzz_xxy_1[j] + 1.5 * fl1_fx * (tg_zz_xxy_0[j] - tg_zz_xxy_1[j] * fl1_fza);

                    tg_zzzz_xxz_0[j] = pb_z * tg_zzz_xxz_0[j] + fr * tg_zzz_xxz_1[j] + 1.5 * fl1_fx * (tg_zz_xxz_0[j] - tg_zz_xxz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xx_1[j];

                    tg_zzzz_xyy_0[j] = pb_z * tg_zzz_xyy_0[j] + fr * tg_zzz_xyy_1[j] + 1.5 * fl1_fx * (tg_zz_xyy_0[j] - tg_zz_xyy_1[j] * fl1_fza);

                    tg_zzzz_xyz_0[j] = pb_z * tg_zzz_xyz_0[j] + fr * tg_zzz_xyz_1[j] + 1.5 * fl1_fx * (tg_zz_xyz_0[j] - tg_zz_xyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_xy_1[j];

                    tg_zzzz_xzz_0[j] = pb_z * tg_zzz_xzz_0[j] + fr * tg_zzz_xzz_1[j] + 1.5 * fl1_fx * (tg_zz_xzz_0[j] - tg_zz_xzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_xz_1[j];

                    tg_zzzz_yyy_0[j] = pb_z * tg_zzz_yyy_0[j] + fr * tg_zzz_yyy_1[j] + 1.5 * fl1_fx * (tg_zz_yyy_0[j] - tg_zz_yyy_1[j] * fl1_fza);

                    tg_zzzz_yyz_0[j] = pb_z * tg_zzz_yyz_0[j] + fr * tg_zzz_yyz_1[j] + 1.5 * fl1_fx * (tg_zz_yyz_0[j] - tg_zz_yyz_1[j] * fl1_fza) + 0.5 * fl1_fxn * tg_zzz_yy_1[j];

                    tg_zzzz_yzz_0[j] = pb_z * tg_zzz_yzz_0[j] + fr * tg_zzz_yzz_1[j] + 1.5 * fl1_fx * (tg_zz_yzz_0[j] - tg_zz_yzz_1[j] * fl1_fza) + fl1_fxn * tg_zzz_yz_1[j];

                    tg_zzzz_zzz_0[j] = pb_z * tg_zzz_zzz_0[j] + fr * tg_zzz_zzz_1[j] + 1.5 * fl1_fx * (tg_zz_zzz_0[j] - tg_zz_zzz_1[j] * fl1_fza) + 1.5 * fl1_fxn * tg_zzz_zz_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

