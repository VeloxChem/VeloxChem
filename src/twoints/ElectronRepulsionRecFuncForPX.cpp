//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForPX.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSPSP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {0, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

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

                auto tg_0_x_0 = primBuffer.data(pidx_g_0_1_m0 + 3 * idx); 

                auto tg_0_y_0 = primBuffer.data(pidx_g_0_1_m0 + 3 * idx + 1); 

                auto tg_0_z_0 = primBuffer.data(pidx_g_0_1_m0 + 3 * idx + 2); 

                auto tg_0_x_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx); 

                auto tg_0_y_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx + 1); 

                auto tg_0_z_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx + 2); 

                auto tg_0_0_1 = primBuffer.data(pidx_g_0_0_m1 + idx); 

                // set up pointers to integrals

                auto tg_x_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx); 

                auto tg_x_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 1); 

                auto tg_x_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 2); 

                auto tg_y_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 3); 

                auto tg_y_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 4); 

                auto tg_y_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 5); 

                auto tg_z_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 6); 

                auto tg_z_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 7); 

                auto tg_z_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 8); 

                #pragma omp simd aligned(fxn, tg_0_0_1, tg_0_x_0, tg_0_x_1, tg_0_y_0, tg_0_y_1, tg_0_z_0, tg_0_z_1, \
                                         tg_x_x_0, tg_x_y_0, tg_x_z_0, tg_y_x_0, tg_y_y_0, tg_y_z_0, tg_z_x_0, tg_z_y_0, \
                                         tg_z_z_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_x_0[j] = pb_x * tg_0_x_0[j] + wp_x[j] * tg_0_x_1[j] + 0.5 * fl1_fxn * tg_0_0_1[j];

                    tg_x_y_0[j] = pb_x * tg_0_y_0[j] + wp_x[j] * tg_0_y_1[j];

                    tg_x_z_0[j] = pb_x * tg_0_z_0[j] + wp_x[j] * tg_0_z_1[j];

                    tg_y_x_0[j] = pb_y * tg_0_x_0[j] + wp_y[j] * tg_0_x_1[j];

                    tg_y_y_0[j] = pb_y * tg_0_y_0[j] + wp_y[j] * tg_0_y_1[j] + 0.5 * fl1_fxn * tg_0_0_1[j];

                    tg_y_z_0[j] = pb_y * tg_0_z_0[j] + wp_y[j] * tg_0_z_1[j];

                    tg_z_x_0[j] = pb_z * tg_0_x_0[j] + wp_z[j] * tg_0_x_1[j];

                    tg_z_y_0[j] = pb_z * tg_0_y_0[j] + wp_z[j] * tg_0_y_1[j];

                    tg_z_z_0[j] = pb_z * tg_0_z_0[j] + wp_z[j] * tg_0_z_1[j] + 0.5 * fl1_fxn * tg_0_0_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSD(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

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

                auto tg_0_xx_0 = primBuffer.data(pidx_g_0_2_m0 + 6 * idx); 

                auto tg_0_xy_0 = primBuffer.data(pidx_g_0_2_m0 + 6 * idx + 1); 

                auto tg_0_xz_0 = primBuffer.data(pidx_g_0_2_m0 + 6 * idx + 2); 

                auto tg_0_yy_0 = primBuffer.data(pidx_g_0_2_m0 + 6 * idx + 3); 

                auto tg_0_yz_0 = primBuffer.data(pidx_g_0_2_m0 + 6 * idx + 4); 

                auto tg_0_zz_0 = primBuffer.data(pidx_g_0_2_m0 + 6 * idx + 5); 

                auto tg_0_xx_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx); 

                auto tg_0_xy_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 1); 

                auto tg_0_xz_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 2); 

                auto tg_0_yy_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 3); 

                auto tg_0_yz_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 4); 

                auto tg_0_zz_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 5); 

                auto tg_0_x_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx); 

                auto tg_0_y_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx + 1); 

                auto tg_0_z_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx + 2); 

                // set up pointers to integrals

                auto tg_x_xx_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx); 

                auto tg_x_xy_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 1); 

                auto tg_x_xz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 2); 

                auto tg_x_yy_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 3); 

                auto tg_x_yz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 4); 

                auto tg_x_zz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 5); 

                auto tg_y_xx_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 6); 

                auto tg_y_xy_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 7); 

                auto tg_y_xz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 8); 

                auto tg_y_yy_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 9); 

                auto tg_y_yz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 10); 

                auto tg_y_zz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 11); 

                auto tg_z_xx_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 12); 

                auto tg_z_xy_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 13); 

                auto tg_z_xz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 14); 

                auto tg_z_yy_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 15); 

                auto tg_z_yz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 16); 

                auto tg_z_zz_0 = primBuffer.data(pidx_g_1_2_m0 + 18 * idx + 17); 

                #pragma omp simd aligned(fxn, tg_0_x_1, tg_0_xx_0, tg_0_xx_1, tg_0_xy_0, tg_0_xy_1, tg_0_xz_0, \
                                         tg_0_xz_1, tg_0_y_1, tg_0_yy_0, tg_0_yy_1, tg_0_yz_0, tg_0_yz_1, tg_0_z_1, \
                                         tg_0_zz_0, tg_0_zz_1, tg_x_xx_0, tg_x_xy_0, tg_x_xz_0, tg_x_yy_0, tg_x_yz_0, \
                                         tg_x_zz_0, tg_y_xx_0, tg_y_xy_0, tg_y_xz_0, tg_y_yy_0, tg_y_yz_0, tg_y_zz_0, \
                                         tg_z_xx_0, tg_z_xy_0, tg_z_xz_0, tg_z_yy_0, tg_z_yz_0, tg_z_zz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_xx_0[j] = pb_x * tg_0_xx_0[j] + wp_x[j] * tg_0_xx_1[j] + fl1_fxn * tg_0_x_1[j];

                    tg_x_xy_0[j] = pb_x * tg_0_xy_0[j] + wp_x[j] * tg_0_xy_1[j] + 0.5 * fl1_fxn * tg_0_y_1[j];

                    tg_x_xz_0[j] = pb_x * tg_0_xz_0[j] + wp_x[j] * tg_0_xz_1[j] + 0.5 * fl1_fxn * tg_0_z_1[j];

                    tg_x_yy_0[j] = pb_x * tg_0_yy_0[j] + wp_x[j] * tg_0_yy_1[j];

                    tg_x_yz_0[j] = pb_x * tg_0_yz_0[j] + wp_x[j] * tg_0_yz_1[j];

                    tg_x_zz_0[j] = pb_x * tg_0_zz_0[j] + wp_x[j] * tg_0_zz_1[j];

                    tg_y_xx_0[j] = pb_y * tg_0_xx_0[j] + wp_y[j] * tg_0_xx_1[j];

                    tg_y_xy_0[j] = pb_y * tg_0_xy_0[j] + wp_y[j] * tg_0_xy_1[j] + 0.5 * fl1_fxn * tg_0_x_1[j];

                    tg_y_xz_0[j] = pb_y * tg_0_xz_0[j] + wp_y[j] * tg_0_xz_1[j];

                    tg_y_yy_0[j] = pb_y * tg_0_yy_0[j] + wp_y[j] * tg_0_yy_1[j] + fl1_fxn * tg_0_y_1[j];

                    tg_y_yz_0[j] = pb_y * tg_0_yz_0[j] + wp_y[j] * tg_0_yz_1[j] + 0.5 * fl1_fxn * tg_0_z_1[j];

                    tg_y_zz_0[j] = pb_y * tg_0_zz_0[j] + wp_y[j] * tg_0_zz_1[j];

                    tg_z_xx_0[j] = pb_z * tg_0_xx_0[j] + wp_z[j] * tg_0_xx_1[j];

                    tg_z_xy_0[j] = pb_z * tg_0_xy_0[j] + wp_z[j] * tg_0_xy_1[j];

                    tg_z_xz_0[j] = pb_z * tg_0_xz_0[j] + wp_z[j] * tg_0_xz_1[j] + 0.5 * fl1_fxn * tg_0_x_1[j];

                    tg_z_yy_0[j] = pb_z * tg_0_yy_0[j] + wp_z[j] * tg_0_yy_1[j];

                    tg_z_yz_0[j] = pb_z * tg_0_yz_0[j] + wp_z[j] * tg_0_yz_1[j] + 0.5 * fl1_fxn * tg_0_y_1[j];

                    tg_z_zz_0[j] = pb_z * tg_0_zz_0[j] + wp_z[j] * tg_0_zz_1[j] + fl1_fxn * tg_0_z_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSDSP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
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
                                             {2, -1, -1, -1},
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_2_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_x_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx); 

                auto tg_x_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 1); 

                auto tg_x_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 2); 

                auto tg_y_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 3); 

                auto tg_y_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 4); 

                auto tg_y_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 5); 

                auto tg_z_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 6); 

                auto tg_z_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 7); 

                auto tg_z_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 8); 

                auto tg_x_x_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx); 

                auto tg_x_y_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 1); 

                auto tg_x_z_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 2); 

                auto tg_y_x_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 3); 

                auto tg_y_y_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 4); 

                auto tg_y_z_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 5); 

                auto tg_z_x_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 6); 

                auto tg_z_y_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 7); 

                auto tg_z_z_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 8); 

                auto tg_0_x_0 = primBuffer.data(pidx_g_0_1_m0 + 3 * idx); 

                auto tg_0_y_0 = primBuffer.data(pidx_g_0_1_m0 + 3 * idx + 1); 

                auto tg_0_z_0 = primBuffer.data(pidx_g_0_1_m0 + 3 * idx + 2); 

                auto tg_0_x_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx); 

                auto tg_0_y_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx + 1); 

                auto tg_0_z_1 = primBuffer.data(pidx_g_0_1_m1 + 3 * idx + 2); 

                auto tg_x_0_1 = primBuffer.data(pidx_g_1_0_m1 + 3 * idx); 

                auto tg_y_0_1 = primBuffer.data(pidx_g_1_0_m1 + 3 * idx + 1); 

                auto tg_z_0_1 = primBuffer.data(pidx_g_1_0_m1 + 3 * idx + 2); 

                // set up pointers to integrals

                auto tg_xx_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx); 

                auto tg_xx_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 1); 

                auto tg_xx_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 2); 

                auto tg_xy_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 3); 

                auto tg_xy_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 4); 

                auto tg_xy_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 5); 

                auto tg_xz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 6); 

                auto tg_xz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 7); 

                auto tg_xz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 8); 

                auto tg_yy_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 9); 

                auto tg_yy_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 10); 

                auto tg_yy_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 11); 

                auto tg_yz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 12); 

                auto tg_yz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 13); 

                auto tg_yz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 14); 

                auto tg_zz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 15); 

                auto tg_zz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 16); 

                auto tg_zz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 17); 

                #pragma omp simd aligned(fxn, fza, tg_0_x_0, tg_0_x_1, tg_0_y_0, tg_0_y_1, tg_0_z_0, tg_0_z_1, \
                                         tg_x_0_1, tg_x_x_0, tg_x_x_1, tg_x_y_0, tg_x_y_1, tg_x_z_0, tg_x_z_1, tg_xx_x_0, \
                                         tg_xx_y_0, tg_xx_z_0, tg_xy_x_0, tg_xy_y_0, tg_xy_z_0, tg_xz_x_0, tg_xz_y_0, \
                                         tg_xz_z_0, tg_y_0_1, tg_y_x_0, tg_y_x_1, tg_y_y_0, tg_y_y_1, tg_y_z_0, tg_y_z_1, \
                                         tg_yy_x_0, tg_yy_y_0, tg_yy_z_0, tg_yz_x_0, tg_yz_y_0, tg_yz_z_0, tg_z_0_1, \
                                         tg_z_x_0, tg_z_x_1, tg_z_y_0, tg_z_y_1, tg_z_z_0, tg_z_z_1, tg_zz_x_0, tg_zz_y_0, \
                                         tg_zz_z_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xx_x_0[j] = pb_x * tg_x_x_0[j] + wp_x[j] * tg_x_x_1[j] + 0.5 * fl1_fx * tg_0_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_x_1[j] + 0.5 * fl1_fxn * tg_x_0_1[j];

                    tg_xx_y_0[j] = pb_x * tg_x_y_0[j] + wp_x[j] * tg_x_y_1[j] + 0.5 * fl1_fx * tg_0_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_y_1[j];

                    tg_xx_z_0[j] = pb_x * tg_x_z_0[j] + wp_x[j] * tg_x_z_1[j] + 0.5 * fl1_fx * tg_0_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_z_1[j];

                    tg_xy_x_0[j] = pb_x * tg_y_x_0[j] + wp_x[j] * tg_y_x_1[j] + 0.5 * fl1_fxn * tg_y_0_1[j];

                    tg_xy_y_0[j] = pb_x * tg_y_y_0[j] + wp_x[j] * tg_y_y_1[j];

                    tg_xy_z_0[j] = pb_x * tg_y_z_0[j] + wp_x[j] * tg_y_z_1[j];

                    tg_xz_x_0[j] = pb_x * tg_z_x_0[j] + wp_x[j] * tg_z_x_1[j] + 0.5 * fl1_fxn * tg_z_0_1[j];

                    tg_xz_y_0[j] = pb_x * tg_z_y_0[j] + wp_x[j] * tg_z_y_1[j];

                    tg_xz_z_0[j] = pb_x * tg_z_z_0[j] + wp_x[j] * tg_z_z_1[j];

                    tg_yy_x_0[j] = pb_y * tg_y_x_0[j] + wp_y[j] * tg_y_x_1[j] + 0.5 * fl1_fx * tg_0_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_x_1[j];

                    tg_yy_y_0[j] = pb_y * tg_y_y_0[j] + wp_y[j] * tg_y_y_1[j] + 0.5 * fl1_fx * tg_0_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_y_1[j] + 0.5 * fl1_fxn * tg_y_0_1[j];

                    tg_yy_z_0[j] = pb_y * tg_y_z_0[j] + wp_y[j] * tg_y_z_1[j] + 0.5 * fl1_fx * tg_0_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_z_1[j];

                    tg_yz_x_0[j] = pb_y * tg_z_x_0[j] + wp_y[j] * tg_z_x_1[j];

                    tg_yz_y_0[j] = pb_y * tg_z_y_0[j] + wp_y[j] * tg_z_y_1[j] + 0.5 * fl1_fxn * tg_z_0_1[j];

                    tg_yz_z_0[j] = pb_y * tg_z_z_0[j] + wp_y[j] * tg_z_z_1[j];

                    tg_zz_x_0[j] = pb_z * tg_z_x_0[j] + wp_z[j] * tg_z_x_1[j] + 0.5 * fl1_fx * tg_0_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_x_1[j];

                    tg_zz_y_0[j] = pb_z * tg_z_y_0[j] + wp_z[j] * tg_z_y_1[j] + 0.5 * fl1_fx * tg_0_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_y_1[j];

                    tg_zz_z_0[j] = pb_z * tg_z_z_0[j] + wp_z[j] * tg_z_z_1[j] + 0.5 * fl1_fx * tg_0_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_0_z_1[j] + 0.5 * fl1_fxn * tg_z_0_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSF(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

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

                auto tg_0_xxx_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx); 

                auto tg_0_xxy_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 1); 

                auto tg_0_xxz_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 2); 

                auto tg_0_xyy_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 3); 

                auto tg_0_xyz_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 4); 

                auto tg_0_xzz_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 5); 

                auto tg_0_yyy_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 6); 

                auto tg_0_yyz_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 7); 

                auto tg_0_yzz_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 8); 

                auto tg_0_zzz_0 = primBuffer.data(pidx_g_0_3_m0 + 10 * idx + 9); 

                auto tg_0_xxx_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx); 

                auto tg_0_xxy_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 1); 

                auto tg_0_xxz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 2); 

                auto tg_0_xyy_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 3); 

                auto tg_0_xyz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 4); 

                auto tg_0_xzz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 5); 

                auto tg_0_yyy_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 6); 

                auto tg_0_yyz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 7); 

                auto tg_0_yzz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 8); 

                auto tg_0_zzz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 9); 

                auto tg_0_xx_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx); 

                auto tg_0_xy_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 1); 

                auto tg_0_xz_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 2); 

                auto tg_0_yy_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 3); 

                auto tg_0_yz_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 4); 

                auto tg_0_zz_1 = primBuffer.data(pidx_g_0_2_m1 + 6 * idx + 5); 

                // set up pointers to integrals

                auto tg_x_xxx_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx); 

                auto tg_x_xxy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 1); 

                auto tg_x_xxz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 2); 

                auto tg_x_xyy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 3); 

                auto tg_x_xyz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 4); 

                auto tg_x_xzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 5); 

                auto tg_x_yyy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 6); 

                auto tg_x_yyz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 7); 

                auto tg_x_yzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 8); 

                auto tg_x_zzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 9); 

                auto tg_y_xxx_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 10); 

                auto tg_y_xxy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 11); 

                auto tg_y_xxz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 12); 

                auto tg_y_xyy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 13); 

                auto tg_y_xyz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 14); 

                auto tg_y_xzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 15); 

                auto tg_y_yyy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 16); 

                auto tg_y_yyz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 17); 

                auto tg_y_yzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 18); 

                auto tg_y_zzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 19); 

                auto tg_z_xxx_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 20); 

                auto tg_z_xxy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 21); 

                auto tg_z_xxz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 22); 

                auto tg_z_xyy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 23); 

                auto tg_z_xyz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 24); 

                auto tg_z_xzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 25); 

                auto tg_z_yyy_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 26); 

                auto tg_z_yyz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 27); 

                auto tg_z_yzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 28); 

                auto tg_z_zzz_0 = primBuffer.data(pidx_g_1_3_m0 + 30 * idx + 29); 

                #pragma omp simd aligned(fxn, tg_0_xx_1, tg_0_xxx_0, tg_0_xxx_1, tg_0_xxy_0, tg_0_xxy_1, \
                                         tg_0_xxz_0, tg_0_xxz_1, tg_0_xy_1, tg_0_xyy_0, tg_0_xyy_1, tg_0_xyz_0, tg_0_xyz_1, \
                                         tg_0_xz_1, tg_0_xzz_0, tg_0_xzz_1, tg_0_yy_1, tg_0_yyy_0, tg_0_yyy_1, tg_0_yyz_0, \
                                         tg_0_yyz_1, tg_0_yz_1, tg_0_yzz_0, tg_0_yzz_1, tg_0_zz_1, tg_0_zzz_0, tg_0_zzz_1, \
                                         tg_x_xxx_0, tg_x_xxy_0, tg_x_xxz_0, tg_x_xyy_0, tg_x_xyz_0, tg_x_xzz_0, tg_x_yyy_0, \
                                         tg_x_yyz_0, tg_x_yzz_0, tg_x_zzz_0, tg_y_xxx_0, tg_y_xxy_0, tg_y_xxz_0, tg_y_xyy_0, \
                                         tg_y_xyz_0, tg_y_xzz_0, tg_y_yyy_0, tg_y_yyz_0, tg_y_yzz_0, tg_y_zzz_0, tg_z_xxx_0, \
                                         tg_z_xxy_0, tg_z_xxz_0, tg_z_xyy_0, tg_z_xyz_0, tg_z_xzz_0, tg_z_yyy_0, tg_z_yyz_0, \
                                         tg_z_yzz_0, tg_z_zzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_xxx_0[j] = pb_x * tg_0_xxx_0[j] + wp_x[j] * tg_0_xxx_1[j] + 1.5 * fl1_fxn * tg_0_xx_1[j];

                    tg_x_xxy_0[j] = pb_x * tg_0_xxy_0[j] + wp_x[j] * tg_0_xxy_1[j] + fl1_fxn * tg_0_xy_1[j];

                    tg_x_xxz_0[j] = pb_x * tg_0_xxz_0[j] + wp_x[j] * tg_0_xxz_1[j] + fl1_fxn * tg_0_xz_1[j];

                    tg_x_xyy_0[j] = pb_x * tg_0_xyy_0[j] + wp_x[j] * tg_0_xyy_1[j] + 0.5 * fl1_fxn * tg_0_yy_1[j];

                    tg_x_xyz_0[j] = pb_x * tg_0_xyz_0[j] + wp_x[j] * tg_0_xyz_1[j] + 0.5 * fl1_fxn * tg_0_yz_1[j];

                    tg_x_xzz_0[j] = pb_x * tg_0_xzz_0[j] + wp_x[j] * tg_0_xzz_1[j] + 0.5 * fl1_fxn * tg_0_zz_1[j];

                    tg_x_yyy_0[j] = pb_x * tg_0_yyy_0[j] + wp_x[j] * tg_0_yyy_1[j];

                    tg_x_yyz_0[j] = pb_x * tg_0_yyz_0[j] + wp_x[j] * tg_0_yyz_1[j];

                    tg_x_yzz_0[j] = pb_x * tg_0_yzz_0[j] + wp_x[j] * tg_0_yzz_1[j];

                    tg_x_zzz_0[j] = pb_x * tg_0_zzz_0[j] + wp_x[j] * tg_0_zzz_1[j];

                    tg_y_xxx_0[j] = pb_y * tg_0_xxx_0[j] + wp_y[j] * tg_0_xxx_1[j];

                    tg_y_xxy_0[j] = pb_y * tg_0_xxy_0[j] + wp_y[j] * tg_0_xxy_1[j] + 0.5 * fl1_fxn * tg_0_xx_1[j];

                    tg_y_xxz_0[j] = pb_y * tg_0_xxz_0[j] + wp_y[j] * tg_0_xxz_1[j];

                    tg_y_xyy_0[j] = pb_y * tg_0_xyy_0[j] + wp_y[j] * tg_0_xyy_1[j] + fl1_fxn * tg_0_xy_1[j];

                    tg_y_xyz_0[j] = pb_y * tg_0_xyz_0[j] + wp_y[j] * tg_0_xyz_1[j] + 0.5 * fl1_fxn * tg_0_xz_1[j];

                    tg_y_xzz_0[j] = pb_y * tg_0_xzz_0[j] + wp_y[j] * tg_0_xzz_1[j];

                    tg_y_yyy_0[j] = pb_y * tg_0_yyy_0[j] + wp_y[j] * tg_0_yyy_1[j] + 1.5 * fl1_fxn * tg_0_yy_1[j];

                    tg_y_yyz_0[j] = pb_y * tg_0_yyz_0[j] + wp_y[j] * tg_0_yyz_1[j] + fl1_fxn * tg_0_yz_1[j];

                    tg_y_yzz_0[j] = pb_y * tg_0_yzz_0[j] + wp_y[j] * tg_0_yzz_1[j] + 0.5 * fl1_fxn * tg_0_zz_1[j];

                    tg_y_zzz_0[j] = pb_y * tg_0_zzz_0[j] + wp_y[j] * tg_0_zzz_1[j];

                    tg_z_xxx_0[j] = pb_z * tg_0_xxx_0[j] + wp_z[j] * tg_0_xxx_1[j];

                    tg_z_xxy_0[j] = pb_z * tg_0_xxy_0[j] + wp_z[j] * tg_0_xxy_1[j];

                    tg_z_xxz_0[j] = pb_z * tg_0_xxz_0[j] + wp_z[j] * tg_0_xxz_1[j] + 0.5 * fl1_fxn * tg_0_xx_1[j];

                    tg_z_xyy_0[j] = pb_z * tg_0_xyy_0[j] + wp_z[j] * tg_0_xyy_1[j];

                    tg_z_xyz_0[j] = pb_z * tg_0_xyz_0[j] + wp_z[j] * tg_0_xyz_1[j] + 0.5 * fl1_fxn * tg_0_xy_1[j];

                    tg_z_xzz_0[j] = pb_z * tg_0_xzz_0[j] + wp_z[j] * tg_0_xzz_1[j] + fl1_fxn * tg_0_xz_1[j];

                    tg_z_yyy_0[j] = pb_z * tg_0_yyy_0[j] + wp_z[j] * tg_0_yyy_1[j];

                    tg_z_yyz_0[j] = pb_z * tg_0_yyz_0[j] + wp_z[j] * tg_0_yyz_1[j] + 0.5 * fl1_fxn * tg_0_yy_1[j];

                    tg_z_yzz_0[j] = pb_z * tg_0_yzz_0[j] + wp_z[j] * tg_0_yzz_1[j] + fl1_fxn * tg_0_yz_1[j];

                    tg_z_zzz_0[j] = pb_z * tg_0_zzz_0[j] + wp_z[j] * tg_0_zzz_1[j] + 1.5 * fl1_fxn * tg_0_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
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
                                             {3, -1, -1, -1},
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_1_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xx_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx); 

                auto tg_xx_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 1); 

                auto tg_xx_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 2); 

                auto tg_xy_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 3); 

                auto tg_xy_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 4); 

                auto tg_xy_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 5); 

                auto tg_xz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 6); 

                auto tg_xz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 7); 

                auto tg_xz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 8); 

                auto tg_yy_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 9); 

                auto tg_yy_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 10); 

                auto tg_yy_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 11); 

                auto tg_yz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 12); 

                auto tg_yz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 13); 

                auto tg_yz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 14); 

                auto tg_zz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 15); 

                auto tg_zz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 16); 

                auto tg_zz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 17); 

                auto tg_xx_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx); 

                auto tg_xx_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 1); 

                auto tg_xx_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 2); 

                auto tg_xy_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 3); 

                auto tg_xy_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 4); 

                auto tg_xy_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 5); 

                auto tg_xz_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 6); 

                auto tg_xz_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 7); 

                auto tg_xz_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 8); 

                auto tg_yy_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 9); 

                auto tg_yy_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 10); 

                auto tg_yy_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 11); 

                auto tg_yz_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 12); 

                auto tg_yz_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 13); 

                auto tg_yz_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 14); 

                auto tg_zz_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 15); 

                auto tg_zz_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 16); 

                auto tg_zz_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 17); 

                auto tg_x_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx); 

                auto tg_x_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 1); 

                auto tg_x_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 2); 

                auto tg_y_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 3); 

                auto tg_y_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 4); 

                auto tg_y_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 5); 

                auto tg_z_x_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 6); 

                auto tg_z_y_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 7); 

                auto tg_z_z_0 = primBuffer.data(pidx_g_1_1_m0 + 9 * idx + 8); 

                auto tg_x_x_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx); 

                auto tg_x_y_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 1); 

                auto tg_x_z_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 2); 

                auto tg_y_x_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 3); 

                auto tg_y_y_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 4); 

                auto tg_y_z_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 5); 

                auto tg_z_x_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 6); 

                auto tg_z_y_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 7); 

                auto tg_z_z_1 = primBuffer.data(pidx_g_1_1_m1 + 9 * idx + 8); 

                auto tg_xx_0_1 = primBuffer.data(pidx_g_2_0_m1 + 6 * idx); 

                auto tg_xy_0_1 = primBuffer.data(pidx_g_2_0_m1 + 6 * idx + 1); 

                auto tg_xz_0_1 = primBuffer.data(pidx_g_2_0_m1 + 6 * idx + 2); 

                auto tg_yy_0_1 = primBuffer.data(pidx_g_2_0_m1 + 6 * idx + 3); 

                auto tg_yz_0_1 = primBuffer.data(pidx_g_2_0_m1 + 6 * idx + 4); 

                auto tg_zz_0_1 = primBuffer.data(pidx_g_2_0_m1 + 6 * idx + 5); 

                // set up pointers to integrals

                auto tg_xxx_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx); 

                auto tg_xxx_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 1); 

                auto tg_xxx_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 2); 

                auto tg_xxy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 3); 

                auto tg_xxy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 4); 

                auto tg_xxy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 5); 

                auto tg_xxz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 6); 

                auto tg_xxz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 7); 

                auto tg_xxz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 8); 

                auto tg_xyy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 9); 

                auto tg_xyy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 10); 

                auto tg_xyy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 11); 

                auto tg_xyz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 12); 

                auto tg_xyz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 13); 

                auto tg_xyz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 14); 

                auto tg_xzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 15); 

                auto tg_xzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 16); 

                auto tg_xzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 17); 

                auto tg_yyy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 18); 

                auto tg_yyy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 19); 

                auto tg_yyy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 20); 

                auto tg_yyz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 21); 

                auto tg_yyz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 22); 

                auto tg_yyz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 23); 

                auto tg_yzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 24); 

                auto tg_yzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 25); 

                auto tg_yzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 26); 

                auto tg_zzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 27); 

                auto tg_zzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 28); 

                auto tg_zzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 29); 

                #pragma omp simd aligned(fxn, fza, tg_x_x_0, tg_x_x_1, tg_x_y_0, tg_x_y_1, tg_x_z_0, tg_x_z_1, \
                                         tg_xx_0_1, tg_xx_x_0, tg_xx_x_1, tg_xx_y_0, tg_xx_y_1, tg_xx_z_0, tg_xx_z_1, \
                                         tg_xxx_x_0, tg_xxx_y_0, tg_xxx_z_0, tg_xxy_x_0, tg_xxy_y_0, tg_xxy_z_0, tg_xxz_x_0, \
                                         tg_xxz_y_0, tg_xxz_z_0, tg_xy_0_1, tg_xy_x_0, tg_xy_x_1, tg_xy_y_0, tg_xy_y_1, \
                                         tg_xy_z_0, tg_xy_z_1, tg_xyy_x_0, tg_xyy_y_0, tg_xyy_z_0, tg_xyz_x_0, tg_xyz_y_0, \
                                         tg_xyz_z_0, tg_xz_0_1, tg_xz_x_0, tg_xz_x_1, tg_xz_y_0, tg_xz_y_1, tg_xz_z_0, \
                                         tg_xz_z_1, tg_xzz_x_0, tg_xzz_y_0, tg_xzz_z_0, tg_y_x_0, tg_y_x_1, tg_y_y_0, \
                                         tg_y_y_1, tg_y_z_0, tg_y_z_1, tg_yy_0_1, tg_yy_x_0, tg_yy_x_1, tg_yy_y_0, \
                                         tg_yy_y_1, tg_yy_z_0, tg_yy_z_1, tg_yyy_x_0, tg_yyy_y_0, tg_yyy_z_0, tg_yyz_x_0, \
                                         tg_yyz_y_0, tg_yyz_z_0, tg_yz_0_1, tg_yz_x_0, tg_yz_x_1, tg_yz_y_0, tg_yz_y_1, \
                                         tg_yz_z_0, tg_yz_z_1, tg_yzz_x_0, tg_yzz_y_0, tg_yzz_z_0, tg_z_x_0, tg_z_x_1, \
                                         tg_z_y_0, tg_z_y_1, tg_z_z_0, tg_z_z_1, tg_zz_0_1, tg_zz_x_0, tg_zz_x_1, tg_zz_y_0, \
                                         tg_zz_y_1, tg_zz_z_0, tg_zz_z_1, tg_zzz_x_0, tg_zzz_y_0, tg_zzz_z_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxx_x_0[j] = pb_x * tg_xx_x_0[j] + wp_x[j] * tg_xx_x_1[j] + fl1_fx * tg_x_x_0[j] - fl1_fx * fl1_fza * tg_x_x_1[j] + 0.5 * fl1_fxn * tg_xx_0_1[j];

                    tg_xxx_y_0[j] = pb_x * tg_xx_y_0[j] + wp_x[j] * tg_xx_y_1[j] + fl1_fx * tg_x_y_0[j] - fl1_fx * fl1_fza * tg_x_y_1[j];

                    tg_xxx_z_0[j] = pb_x * tg_xx_z_0[j] + wp_x[j] * tg_xx_z_1[j] + fl1_fx * tg_x_z_0[j] - fl1_fx * fl1_fza * tg_x_z_1[j];

                    tg_xxy_x_0[j] = pb_x * tg_xy_x_0[j] + wp_x[j] * tg_xy_x_1[j] + 0.5 * fl1_fx * tg_y_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_x_1[j] + 0.5 * fl1_fxn * tg_xy_0_1[j];

                    tg_xxy_y_0[j] = pb_x * tg_xy_y_0[j] + wp_x[j] * tg_xy_y_1[j] + 0.5 * fl1_fx * tg_y_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_y_1[j];

                    tg_xxy_z_0[j] = pb_x * tg_xy_z_0[j] + wp_x[j] * tg_xy_z_1[j] + 0.5 * fl1_fx * tg_y_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_y_z_1[j];

                    tg_xxz_x_0[j] = pb_x * tg_xz_x_0[j] + wp_x[j] * tg_xz_x_1[j] + 0.5 * fl1_fx * tg_z_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_x_1[j] + 0.5 * fl1_fxn * tg_xz_0_1[j];

                    tg_xxz_y_0[j] = pb_x * tg_xz_y_0[j] + wp_x[j] * tg_xz_y_1[j] + 0.5 * fl1_fx * tg_z_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_y_1[j];

                    tg_xxz_z_0[j] = pb_x * tg_xz_z_0[j] + wp_x[j] * tg_xz_z_1[j] + 0.5 * fl1_fx * tg_z_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_z_1[j];

                    tg_xyy_x_0[j] = pb_x * tg_yy_x_0[j] + wp_x[j] * tg_yy_x_1[j] + 0.5 * fl1_fxn * tg_yy_0_1[j];

                    tg_xyy_y_0[j] = pb_x * tg_yy_y_0[j] + wp_x[j] * tg_yy_y_1[j];

                    tg_xyy_z_0[j] = pb_x * tg_yy_z_0[j] + wp_x[j] * tg_yy_z_1[j];

                    tg_xyz_x_0[j] = pb_x * tg_yz_x_0[j] + wp_x[j] * tg_yz_x_1[j] + 0.5 * fl1_fxn * tg_yz_0_1[j];

                    tg_xyz_y_0[j] = pb_x * tg_yz_y_0[j] + wp_x[j] * tg_yz_y_1[j];

                    tg_xyz_z_0[j] = pb_x * tg_yz_z_0[j] + wp_x[j] * tg_yz_z_1[j];

                    tg_xzz_x_0[j] = pb_x * tg_zz_x_0[j] + wp_x[j] * tg_zz_x_1[j] + 0.5 * fl1_fxn * tg_zz_0_1[j];

                    tg_xzz_y_0[j] = pb_x * tg_zz_y_0[j] + wp_x[j] * tg_zz_y_1[j];

                    tg_xzz_z_0[j] = pb_x * tg_zz_z_0[j] + wp_x[j] * tg_zz_z_1[j];

                    tg_yyy_x_0[j] = pb_y * tg_yy_x_0[j] + wp_y[j] * tg_yy_x_1[j] + fl1_fx * tg_y_x_0[j] - fl1_fx * fl1_fza * tg_y_x_1[j];

                    tg_yyy_y_0[j] = pb_y * tg_yy_y_0[j] + wp_y[j] * tg_yy_y_1[j] + fl1_fx * tg_y_y_0[j] - fl1_fx * fl1_fza * tg_y_y_1[j] + 0.5 * fl1_fxn * tg_yy_0_1[j];

                    tg_yyy_z_0[j] = pb_y * tg_yy_z_0[j] + wp_y[j] * tg_yy_z_1[j] + fl1_fx * tg_y_z_0[j] - fl1_fx * fl1_fza * tg_y_z_1[j];

                    tg_yyz_x_0[j] = pb_y * tg_yz_x_0[j] + wp_y[j] * tg_yz_x_1[j] + 0.5 * fl1_fx * tg_z_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_x_1[j];

                    tg_yyz_y_0[j] = pb_y * tg_yz_y_0[j] + wp_y[j] * tg_yz_y_1[j] + 0.5 * fl1_fx * tg_z_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_y_1[j] + 0.5 * fl1_fxn * tg_yz_0_1[j];

                    tg_yyz_z_0[j] = pb_y * tg_yz_z_0[j] + wp_y[j] * tg_yz_z_1[j] + 0.5 * fl1_fx * tg_z_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_z_z_1[j];

                    tg_yzz_x_0[j] = pb_y * tg_zz_x_0[j] + wp_y[j] * tg_zz_x_1[j];

                    tg_yzz_y_0[j] = pb_y * tg_zz_y_0[j] + wp_y[j] * tg_zz_y_1[j] + 0.5 * fl1_fxn * tg_zz_0_1[j];

                    tg_yzz_z_0[j] = pb_y * tg_zz_z_0[j] + wp_y[j] * tg_zz_z_1[j];

                    tg_zzz_x_0[j] = pb_z * tg_zz_x_0[j] + wp_z[j] * tg_zz_x_1[j] + fl1_fx * tg_z_x_0[j] - fl1_fx * fl1_fza * tg_z_x_1[j];

                    tg_zzz_y_0[j] = pb_z * tg_zz_y_0[j] + wp_z[j] * tg_zz_y_1[j] + fl1_fx * tg_z_y_0[j] - fl1_fx * fl1_fza * tg_z_y_1[j];

                    tg_zzz_z_0[j] = pb_z * tg_zz_z_0[j] + wp_z[j] * tg_zz_z_1[j] + fl1_fx * tg_z_z_0[j] - fl1_fx * fl1_fza * tg_z_z_1[j] + 0.5 * fl1_fxn * tg_zz_0_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSG(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

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

                auto tg_0_xxxx_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx); 

                auto tg_0_xxxy_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 1); 

                auto tg_0_xxxz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 2); 

                auto tg_0_xxyy_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 3); 

                auto tg_0_xxyz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 4); 

                auto tg_0_xxzz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 5); 

                auto tg_0_xyyy_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 6); 

                auto tg_0_xyyz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 7); 

                auto tg_0_xyzz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 8); 

                auto tg_0_xzzz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 9); 

                auto tg_0_yyyy_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 10); 

                auto tg_0_yyyz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 11); 

                auto tg_0_yyzz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 12); 

                auto tg_0_yzzz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 13); 

                auto tg_0_zzzz_0 = primBuffer.data(pidx_g_0_4_m0 + 15 * idx + 14); 

                auto tg_0_xxxx_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx); 

                auto tg_0_xxxy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 1); 

                auto tg_0_xxxz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 2); 

                auto tg_0_xxyy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 3); 

                auto tg_0_xxyz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 4); 

                auto tg_0_xxzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 5); 

                auto tg_0_xyyy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 6); 

                auto tg_0_xyyz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 7); 

                auto tg_0_xyzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 8); 

                auto tg_0_xzzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 9); 

                auto tg_0_yyyy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 10); 

                auto tg_0_yyyz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 11); 

                auto tg_0_yyzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 12); 

                auto tg_0_yzzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 13); 

                auto tg_0_zzzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 14); 

                auto tg_0_xxx_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx); 

                auto tg_0_xxy_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 1); 

                auto tg_0_xxz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 2); 

                auto tg_0_xyy_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 3); 

                auto tg_0_xyz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 4); 

                auto tg_0_xzz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 5); 

                auto tg_0_yyy_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 6); 

                auto tg_0_yyz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 7); 

                auto tg_0_yzz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 8); 

                auto tg_0_zzz_1 = primBuffer.data(pidx_g_0_3_m1 + 10 * idx + 9); 

                // set up pointers to integrals

                auto tg_x_xxxx_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx); 

                auto tg_x_xxxy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 1); 

                auto tg_x_xxxz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 2); 

                auto tg_x_xxyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 3); 

                auto tg_x_xxyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 4); 

                auto tg_x_xxzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 5); 

                auto tg_x_xyyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 6); 

                auto tg_x_xyyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 7); 

                auto tg_x_xyzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 8); 

                auto tg_x_xzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 9); 

                auto tg_x_yyyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 10); 

                auto tg_x_yyyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 11); 

                auto tg_x_yyzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 12); 

                auto tg_x_yzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 13); 

                auto tg_x_zzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 14); 

                auto tg_y_xxxx_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 15); 

                auto tg_y_xxxy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 16); 

                auto tg_y_xxxz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 17); 

                auto tg_y_xxyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 18); 

                auto tg_y_xxyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 19); 

                auto tg_y_xxzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 20); 

                auto tg_y_xyyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 21); 

                auto tg_y_xyyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 22); 

                auto tg_y_xyzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 23); 

                auto tg_y_xzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 24); 

                auto tg_y_yyyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 25); 

                auto tg_y_yyyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 26); 

                auto tg_y_yyzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 27); 

                auto tg_y_yzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 28); 

                auto tg_y_zzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 29); 

                auto tg_z_xxxx_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 30); 

                auto tg_z_xxxy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 31); 

                auto tg_z_xxxz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 32); 

                auto tg_z_xxyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 33); 

                auto tg_z_xxyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 34); 

                auto tg_z_xxzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 35); 

                auto tg_z_xyyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 36); 

                auto tg_z_xyyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 37); 

                auto tg_z_xyzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 38); 

                auto tg_z_xzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 39); 

                auto tg_z_yyyy_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 40); 

                auto tg_z_yyyz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 41); 

                auto tg_z_yyzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 42); 

                auto tg_z_yzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 43); 

                auto tg_z_zzzz_0 = primBuffer.data(pidx_g_1_4_m0 + 45 * idx + 44); 

                #pragma omp simd aligned(fxn, tg_0_xxx_1, tg_0_xxxx_0, tg_0_xxxx_1, tg_0_xxxy_0, tg_0_xxxy_1, \
                                         tg_0_xxxz_0, tg_0_xxxz_1, tg_0_xxy_1, tg_0_xxyy_0, tg_0_xxyy_1, tg_0_xxyz_0, \
                                         tg_0_xxyz_1, tg_0_xxz_1, tg_0_xxzz_0, tg_0_xxzz_1, tg_0_xyy_1, tg_0_xyyy_0, \
                                         tg_0_xyyy_1, tg_0_xyyz_0, tg_0_xyyz_1, tg_0_xyz_1, tg_0_xyzz_0, tg_0_xyzz_1, \
                                         tg_0_xzz_1, tg_0_xzzz_0, tg_0_xzzz_1, tg_0_yyy_1, tg_0_yyyy_0, tg_0_yyyy_1, \
                                         tg_0_yyyz_0, tg_0_yyyz_1, tg_0_yyz_1, tg_0_yyzz_0, tg_0_yyzz_1, tg_0_yzz_1, \
                                         tg_0_yzzz_0, tg_0_yzzz_1, tg_0_zzz_1, tg_0_zzzz_0, tg_0_zzzz_1, tg_x_xxxx_0, \
                                         tg_x_xxxy_0, tg_x_xxxz_0, tg_x_xxyy_0, tg_x_xxyz_0, tg_x_xxzz_0, tg_x_xyyy_0, \
                                         tg_x_xyyz_0, tg_x_xyzz_0, tg_x_xzzz_0, tg_x_yyyy_0, tg_x_yyyz_0, tg_x_yyzz_0, \
                                         tg_x_yzzz_0, tg_x_zzzz_0, tg_y_xxxx_0, tg_y_xxxy_0, tg_y_xxxz_0, tg_y_xxyy_0, \
                                         tg_y_xxyz_0, tg_y_xxzz_0, tg_y_xyyy_0, tg_y_xyyz_0, tg_y_xyzz_0, tg_y_xzzz_0, \
                                         tg_y_yyyy_0, tg_y_yyyz_0, tg_y_yyzz_0, tg_y_yzzz_0, tg_y_zzzz_0, tg_z_xxxx_0, \
                                         tg_z_xxxy_0, tg_z_xxxz_0, tg_z_xxyy_0, tg_z_xxyz_0, tg_z_xxzz_0, tg_z_xyyy_0, \
                                         tg_z_xyyz_0, tg_z_xyzz_0, tg_z_xzzz_0, tg_z_yyyy_0, tg_z_yyyz_0, tg_z_yyzz_0, \
                                         tg_z_yzzz_0, tg_z_zzzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_xxxx_0[j] = pb_x * tg_0_xxxx_0[j] + wp_x[j] * tg_0_xxxx_1[j] + 2.0 * fl1_fxn * tg_0_xxx_1[j];

                    tg_x_xxxy_0[j] = pb_x * tg_0_xxxy_0[j] + wp_x[j] * tg_0_xxxy_1[j] + 1.5 * fl1_fxn * tg_0_xxy_1[j];

                    tg_x_xxxz_0[j] = pb_x * tg_0_xxxz_0[j] + wp_x[j] * tg_0_xxxz_1[j] + 1.5 * fl1_fxn * tg_0_xxz_1[j];

                    tg_x_xxyy_0[j] = pb_x * tg_0_xxyy_0[j] + wp_x[j] * tg_0_xxyy_1[j] + fl1_fxn * tg_0_xyy_1[j];

                    tg_x_xxyz_0[j] = pb_x * tg_0_xxyz_0[j] + wp_x[j] * tg_0_xxyz_1[j] + fl1_fxn * tg_0_xyz_1[j];

                    tg_x_xxzz_0[j] = pb_x * tg_0_xxzz_0[j] + wp_x[j] * tg_0_xxzz_1[j] + fl1_fxn * tg_0_xzz_1[j];

                    tg_x_xyyy_0[j] = pb_x * tg_0_xyyy_0[j] + wp_x[j] * tg_0_xyyy_1[j] + 0.5 * fl1_fxn * tg_0_yyy_1[j];

                    tg_x_xyyz_0[j] = pb_x * tg_0_xyyz_0[j] + wp_x[j] * tg_0_xyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyz_1[j];

                    tg_x_xyzz_0[j] = pb_x * tg_0_xyzz_0[j] + wp_x[j] * tg_0_xyzz_1[j] + 0.5 * fl1_fxn * tg_0_yzz_1[j];

                    tg_x_xzzz_0[j] = pb_x * tg_0_xzzz_0[j] + wp_x[j] * tg_0_xzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzz_1[j];

                    tg_x_yyyy_0[j] = pb_x * tg_0_yyyy_0[j] + wp_x[j] * tg_0_yyyy_1[j];

                    tg_x_yyyz_0[j] = pb_x * tg_0_yyyz_0[j] + wp_x[j] * tg_0_yyyz_1[j];

                    tg_x_yyzz_0[j] = pb_x * tg_0_yyzz_0[j] + wp_x[j] * tg_0_yyzz_1[j];

                    tg_x_yzzz_0[j] = pb_x * tg_0_yzzz_0[j] + wp_x[j] * tg_0_yzzz_1[j];

                    tg_x_zzzz_0[j] = pb_x * tg_0_zzzz_0[j] + wp_x[j] * tg_0_zzzz_1[j];

                    tg_y_xxxx_0[j] = pb_y * tg_0_xxxx_0[j] + wp_y[j] * tg_0_xxxx_1[j];

                    tg_y_xxxy_0[j] = pb_y * tg_0_xxxy_0[j] + wp_y[j] * tg_0_xxxy_1[j] + 0.5 * fl1_fxn * tg_0_xxx_1[j];

                    tg_y_xxxz_0[j] = pb_y * tg_0_xxxz_0[j] + wp_y[j] * tg_0_xxxz_1[j];

                    tg_y_xxyy_0[j] = pb_y * tg_0_xxyy_0[j] + wp_y[j] * tg_0_xxyy_1[j] + fl1_fxn * tg_0_xxy_1[j];

                    tg_y_xxyz_0[j] = pb_y * tg_0_xxyz_0[j] + wp_y[j] * tg_0_xxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxz_1[j];

                    tg_y_xxzz_0[j] = pb_y * tg_0_xxzz_0[j] + wp_y[j] * tg_0_xxzz_1[j];

                    tg_y_xyyy_0[j] = pb_y * tg_0_xyyy_0[j] + wp_y[j] * tg_0_xyyy_1[j] + 1.5 * fl1_fxn * tg_0_xyy_1[j];

                    tg_y_xyyz_0[j] = pb_y * tg_0_xyyz_0[j] + wp_y[j] * tg_0_xyyz_1[j] + fl1_fxn * tg_0_xyz_1[j];

                    tg_y_xyzz_0[j] = pb_y * tg_0_xyzz_0[j] + wp_y[j] * tg_0_xyzz_1[j] + 0.5 * fl1_fxn * tg_0_xzz_1[j];

                    tg_y_xzzz_0[j] = pb_y * tg_0_xzzz_0[j] + wp_y[j] * tg_0_xzzz_1[j];

                    tg_y_yyyy_0[j] = pb_y * tg_0_yyyy_0[j] + wp_y[j] * tg_0_yyyy_1[j] + 2.0 * fl1_fxn * tg_0_yyy_1[j];

                    tg_y_yyyz_0[j] = pb_y * tg_0_yyyz_0[j] + wp_y[j] * tg_0_yyyz_1[j] + 1.5 * fl1_fxn * tg_0_yyz_1[j];

                    tg_y_yyzz_0[j] = pb_y * tg_0_yyzz_0[j] + wp_y[j] * tg_0_yyzz_1[j] + fl1_fxn * tg_0_yzz_1[j];

                    tg_y_yzzz_0[j] = pb_y * tg_0_yzzz_0[j] + wp_y[j] * tg_0_yzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzz_1[j];

                    tg_y_zzzz_0[j] = pb_y * tg_0_zzzz_0[j] + wp_y[j] * tg_0_zzzz_1[j];

                    tg_z_xxxx_0[j] = pb_z * tg_0_xxxx_0[j] + wp_z[j] * tg_0_xxxx_1[j];

                    tg_z_xxxy_0[j] = pb_z * tg_0_xxxy_0[j] + wp_z[j] * tg_0_xxxy_1[j];

                    tg_z_xxxz_0[j] = pb_z * tg_0_xxxz_0[j] + wp_z[j] * tg_0_xxxz_1[j] + 0.5 * fl1_fxn * tg_0_xxx_1[j];

                    tg_z_xxyy_0[j] = pb_z * tg_0_xxyy_0[j] + wp_z[j] * tg_0_xxyy_1[j];

                    tg_z_xxyz_0[j] = pb_z * tg_0_xxyz_0[j] + wp_z[j] * tg_0_xxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxy_1[j];

                    tg_z_xxzz_0[j] = pb_z * tg_0_xxzz_0[j] + wp_z[j] * tg_0_xxzz_1[j] + fl1_fxn * tg_0_xxz_1[j];

                    tg_z_xyyy_0[j] = pb_z * tg_0_xyyy_0[j] + wp_z[j] * tg_0_xyyy_1[j];

                    tg_z_xyyz_0[j] = pb_z * tg_0_xyyz_0[j] + wp_z[j] * tg_0_xyyz_1[j] + 0.5 * fl1_fxn * tg_0_xyy_1[j];

                    tg_z_xyzz_0[j] = pb_z * tg_0_xyzz_0[j] + wp_z[j] * tg_0_xyzz_1[j] + fl1_fxn * tg_0_xyz_1[j];

                    tg_z_xzzz_0[j] = pb_z * tg_0_xzzz_0[j] + wp_z[j] * tg_0_xzzz_1[j] + 1.5 * fl1_fxn * tg_0_xzz_1[j];

                    tg_z_yyyy_0[j] = pb_z * tg_0_yyyy_0[j] + wp_z[j] * tg_0_yyyy_1[j];

                    tg_z_yyyz_0[j] = pb_z * tg_0_yyyz_0[j] + wp_z[j] * tg_0_yyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyy_1[j];

                    tg_z_yyzz_0[j] = pb_z * tg_0_yyzz_0[j] + wp_z[j] * tg_0_yyzz_1[j] + fl1_fxn * tg_0_yyz_1[j];

                    tg_z_yzzz_0[j] = pb_z * tg_0_yzzz_0[j] + wp_z[j] * tg_0_yzzz_1[j] + 1.5 * fl1_fxn * tg_0_yzz_1[j];

                    tg_z_zzzz_0[j] = pb_z * tg_0_zzzz_0[j] + wp_z[j] * tg_0_zzzz_1[j] + 2.0 * fl1_fxn * tg_0_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
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
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_2_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xxx_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx); 

                auto tg_xxx_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 1); 

                auto tg_xxx_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 2); 

                auto tg_xxy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 3); 

                auto tg_xxy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 4); 

                auto tg_xxy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 5); 

                auto tg_xxz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 6); 

                auto tg_xxz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 7); 

                auto tg_xxz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 8); 

                auto tg_xyy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 9); 

                auto tg_xyy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 10); 

                auto tg_xyy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 11); 

                auto tg_xyz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 12); 

                auto tg_xyz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 13); 

                auto tg_xyz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 14); 

                auto tg_xzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 15); 

                auto tg_xzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 16); 

                auto tg_xzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 17); 

                auto tg_yyy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 18); 

                auto tg_yyy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 19); 

                auto tg_yyy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 20); 

                auto tg_yyz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 21); 

                auto tg_yyz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 22); 

                auto tg_yyz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 23); 

                auto tg_yzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 24); 

                auto tg_yzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 25); 

                auto tg_yzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 26); 

                auto tg_zzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 27); 

                auto tg_zzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 28); 

                auto tg_zzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 29); 

                auto tg_xxx_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx); 

                auto tg_xxx_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 1); 

                auto tg_xxx_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 2); 

                auto tg_xxy_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 3); 

                auto tg_xxy_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 4); 

                auto tg_xxy_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 5); 

                auto tg_xxz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 6); 

                auto tg_xxz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 7); 

                auto tg_xxz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 8); 

                auto tg_xyy_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 9); 

                auto tg_xyy_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 10); 

                auto tg_xyy_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 11); 

                auto tg_xyz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 12); 

                auto tg_xyz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 13); 

                auto tg_xyz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 14); 

                auto tg_xzz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 15); 

                auto tg_xzz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 16); 

                auto tg_xzz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 17); 

                auto tg_yyy_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 18); 

                auto tg_yyy_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 19); 

                auto tg_yyy_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 20); 

                auto tg_yyz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 21); 

                auto tg_yyz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 22); 

                auto tg_yyz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 23); 

                auto tg_yzz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 24); 

                auto tg_yzz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 25); 

                auto tg_yzz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 26); 

                auto tg_zzz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 27); 

                auto tg_zzz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 28); 

                auto tg_zzz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 29); 

                auto tg_xx_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx); 

                auto tg_xx_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 1); 

                auto tg_xx_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 2); 

                auto tg_xy_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 3); 

                auto tg_xy_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 4); 

                auto tg_xy_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 5); 

                auto tg_xz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 6); 

                auto tg_xz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 7); 

                auto tg_xz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 8); 

                auto tg_yy_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 9); 

                auto tg_yy_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 10); 

                auto tg_yy_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 11); 

                auto tg_yz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 12); 

                auto tg_yz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 13); 

                auto tg_yz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 14); 

                auto tg_zz_x_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 15); 

                auto tg_zz_y_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 16); 

                auto tg_zz_z_0 = primBuffer.data(pidx_g_2_1_m0 + 18 * idx + 17); 

                auto tg_xx_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx); 

                auto tg_xx_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 1); 

                auto tg_xx_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 2); 

                auto tg_xy_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 3); 

                auto tg_xy_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 4); 

                auto tg_xy_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 5); 

                auto tg_xz_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 6); 

                auto tg_xz_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 7); 

                auto tg_xz_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 8); 

                auto tg_yy_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 9); 

                auto tg_yy_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 10); 

                auto tg_yy_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 11); 

                auto tg_yz_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 12); 

                auto tg_yz_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 13); 

                auto tg_yz_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 14); 

                auto tg_zz_x_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 15); 

                auto tg_zz_y_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 16); 

                auto tg_zz_z_1 = primBuffer.data(pidx_g_2_1_m1 + 18 * idx + 17); 

                auto tg_xxx_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx); 

                auto tg_xxy_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 1); 

                auto tg_xxz_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 2); 

                auto tg_xyy_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 3); 

                auto tg_xyz_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 4); 

                auto tg_xzz_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 5); 

                auto tg_yyy_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 6); 

                auto tg_yyz_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 7); 

                auto tg_yzz_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 8); 

                auto tg_zzz_0_1 = primBuffer.data(pidx_g_3_0_m1 + 10 * idx + 9); 

                // set up pointers to integrals

                auto tg_xxxx_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx); 

                auto tg_xxxx_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 1); 

                auto tg_xxxx_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 2); 

                auto tg_xxxy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 3); 

                auto tg_xxxy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 4); 

                auto tg_xxxy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 5); 

                auto tg_xxxz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 6); 

                auto tg_xxxz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 7); 

                auto tg_xxxz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 8); 

                auto tg_xxyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 9); 

                auto tg_xxyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 10); 

                auto tg_xxyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 11); 

                auto tg_xxyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 12); 

                auto tg_xxyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 13); 

                auto tg_xxyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 14); 

                auto tg_xxzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 15); 

                auto tg_xxzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 16); 

                auto tg_xxzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 17); 

                auto tg_xyyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 18); 

                auto tg_xyyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 19); 

                auto tg_xyyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 20); 

                auto tg_xyyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 21); 

                auto tg_xyyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 22); 

                auto tg_xyyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 23); 

                auto tg_xyzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 24); 

                auto tg_xyzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 25); 

                auto tg_xyzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 26); 

                auto tg_xzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 27); 

                auto tg_xzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 28); 

                auto tg_xzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 29); 

                auto tg_yyyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 30); 

                auto tg_yyyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 31); 

                auto tg_yyyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 32); 

                auto tg_yyyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 33); 

                auto tg_yyyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 34); 

                auto tg_yyyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 35); 

                auto tg_yyzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 36); 

                auto tg_yyzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 37); 

                auto tg_yyzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 38); 

                auto tg_yzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 39); 

                auto tg_yzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 40); 

                auto tg_yzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 41); 

                auto tg_zzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 42); 

                auto tg_zzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 43); 

                auto tg_zzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 44); 

                #pragma omp simd aligned(fxn, fza, tg_xx_x_0, tg_xx_x_1, tg_xx_y_0, tg_xx_y_1, tg_xx_z_0, tg_xx_z_1, \
                                         tg_xxx_0_1, tg_xxx_x_0, tg_xxx_x_1, tg_xxx_y_0, tg_xxx_y_1, tg_xxx_z_0, tg_xxx_z_1, \
                                         tg_xxxx_x_0, tg_xxxx_y_0, tg_xxxx_z_0, tg_xxxy_x_0, tg_xxxy_y_0, tg_xxxy_z_0, \
                                         tg_xxxz_x_0, tg_xxxz_y_0, tg_xxxz_z_0, tg_xxy_0_1, tg_xxy_x_0, tg_xxy_x_1, \
                                         tg_xxy_y_0, tg_xxy_y_1, tg_xxy_z_0, tg_xxy_z_1, tg_xxyy_x_0, tg_xxyy_y_0, \
                                         tg_xxyy_z_0, tg_xxyz_x_0, tg_xxyz_y_0, tg_xxyz_z_0, tg_xxz_0_1, tg_xxz_x_0, \
                                         tg_xxz_x_1, tg_xxz_y_0, tg_xxz_y_1, tg_xxz_z_0, tg_xxz_z_1, tg_xxzz_x_0, \
                                         tg_xxzz_y_0, tg_xxzz_z_0, tg_xy_x_0, tg_xy_x_1, tg_xy_y_0, tg_xy_y_1, tg_xy_z_0, \
                                         tg_xy_z_1, tg_xyy_0_1, tg_xyy_x_0, tg_xyy_x_1, tg_xyy_y_0, tg_xyy_y_1, tg_xyy_z_0, \
                                         tg_xyy_z_1, tg_xyyy_x_0, tg_xyyy_y_0, tg_xyyy_z_0, tg_xyyz_x_0, tg_xyyz_y_0, \
                                         tg_xyyz_z_0, tg_xyz_0_1, tg_xyz_x_0, tg_xyz_x_1, tg_xyz_y_0, tg_xyz_y_1, tg_xyz_z_0, \
                                         tg_xyz_z_1, tg_xyzz_x_0, tg_xyzz_y_0, tg_xyzz_z_0, tg_xz_x_0, tg_xz_x_1, tg_xz_y_0, \
                                         tg_xz_y_1, tg_xz_z_0, tg_xz_z_1, tg_xzz_0_1, tg_xzz_x_0, tg_xzz_x_1, tg_xzz_y_0, \
                                         tg_xzz_y_1, tg_xzz_z_0, tg_xzz_z_1, tg_xzzz_x_0, tg_xzzz_y_0, tg_xzzz_z_0, \
                                         tg_yy_x_0, tg_yy_x_1, tg_yy_y_0, tg_yy_y_1, tg_yy_z_0, tg_yy_z_1, tg_yyy_0_1, \
                                         tg_yyy_x_0, tg_yyy_x_1, tg_yyy_y_0, tg_yyy_y_1, tg_yyy_z_0, tg_yyy_z_1, \
                                         tg_yyyy_x_0, tg_yyyy_y_0, tg_yyyy_z_0, tg_yyyz_x_0, tg_yyyz_y_0, tg_yyyz_z_0, \
                                         tg_yyz_0_1, tg_yyz_x_0, tg_yyz_x_1, tg_yyz_y_0, tg_yyz_y_1, tg_yyz_z_0, tg_yyz_z_1, \
                                         tg_yyzz_x_0, tg_yyzz_y_0, tg_yyzz_z_0, tg_yz_x_0, tg_yz_x_1, tg_yz_y_0, tg_yz_y_1, \
                                         tg_yz_z_0, tg_yz_z_1, tg_yzz_0_1, tg_yzz_x_0, tg_yzz_x_1, tg_yzz_y_0, tg_yzz_y_1, \
                                         tg_yzz_z_0, tg_yzz_z_1, tg_yzzz_x_0, tg_yzzz_y_0, tg_yzzz_z_0, tg_zz_x_0, \
                                         tg_zz_x_1, tg_zz_y_0, tg_zz_y_1, tg_zz_z_0, tg_zz_z_1, tg_zzz_0_1, tg_zzz_x_0, \
                                         tg_zzz_x_1, tg_zzz_y_0, tg_zzz_y_1, tg_zzz_z_0, tg_zzz_z_1, tg_zzzz_x_0, \
                                         tg_zzzz_y_0, tg_zzzz_z_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxx_x_0[j] = pb_x * tg_xxx_x_0[j] + wp_x[j] * tg_xxx_x_1[j] + 1.5 * fl1_fx * tg_xx_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_x_1[j] + 0.5 * fl1_fxn * tg_xxx_0_1[j];

                    tg_xxxx_y_0[j] = pb_x * tg_xxx_y_0[j] + wp_x[j] * tg_xxx_y_1[j] + 1.5 * fl1_fx * tg_xx_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_y_1[j];

                    tg_xxxx_z_0[j] = pb_x * tg_xxx_z_0[j] + wp_x[j] * tg_xxx_z_1[j] + 1.5 * fl1_fx * tg_xx_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xx_z_1[j];

                    tg_xxxy_x_0[j] = pb_x * tg_xxy_x_0[j] + wp_x[j] * tg_xxy_x_1[j] + fl1_fx * tg_xy_x_0[j] - fl1_fx * fl1_fza * tg_xy_x_1[j] + 0.5 * fl1_fxn * tg_xxy_0_1[j];

                    tg_xxxy_y_0[j] = pb_x * tg_xxy_y_0[j] + wp_x[j] * tg_xxy_y_1[j] + fl1_fx * tg_xy_y_0[j] - fl1_fx * fl1_fza * tg_xy_y_1[j];

                    tg_xxxy_z_0[j] = pb_x * tg_xxy_z_0[j] + wp_x[j] * tg_xxy_z_1[j] + fl1_fx * tg_xy_z_0[j] - fl1_fx * fl1_fza * tg_xy_z_1[j];

                    tg_xxxz_x_0[j] = pb_x * tg_xxz_x_0[j] + wp_x[j] * tg_xxz_x_1[j] + fl1_fx * tg_xz_x_0[j] - fl1_fx * fl1_fza * tg_xz_x_1[j] + 0.5 * fl1_fxn * tg_xxz_0_1[j];

                    tg_xxxz_y_0[j] = pb_x * tg_xxz_y_0[j] + wp_x[j] * tg_xxz_y_1[j] + fl1_fx * tg_xz_y_0[j] - fl1_fx * fl1_fza * tg_xz_y_1[j];

                    tg_xxxz_z_0[j] = pb_x * tg_xxz_z_0[j] + wp_x[j] * tg_xxz_z_1[j] + fl1_fx * tg_xz_z_0[j] - fl1_fx * fl1_fza * tg_xz_z_1[j];

                    tg_xxyy_x_0[j] = pb_x * tg_xyy_x_0[j] + wp_x[j] * tg_xyy_x_1[j] + 0.5 * fl1_fx * tg_yy_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_x_1[j] + 0.5 * fl1_fxn * tg_xyy_0_1[j];

                    tg_xxyy_y_0[j] = pb_x * tg_xyy_y_0[j] + wp_x[j] * tg_xyy_y_1[j] + 0.5 * fl1_fx * tg_yy_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_y_1[j];

                    tg_xxyy_z_0[j] = pb_x * tg_xyy_z_0[j] + wp_x[j] * tg_xyy_z_1[j] + 0.5 * fl1_fx * tg_yy_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yy_z_1[j];

                    tg_xxyz_x_0[j] = pb_x * tg_xyz_x_0[j] + wp_x[j] * tg_xyz_x_1[j] + 0.5 * fl1_fx * tg_yz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_x_1[j] + 0.5 * fl1_fxn * tg_xyz_0_1[j];

                    tg_xxyz_y_0[j] = pb_x * tg_xyz_y_0[j] + wp_x[j] * tg_xyz_y_1[j] + 0.5 * fl1_fx * tg_yz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_y_1[j];

                    tg_xxyz_z_0[j] = pb_x * tg_xyz_z_0[j] + wp_x[j] * tg_xyz_z_1[j] + 0.5 * fl1_fx * tg_yz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yz_z_1[j];

                    tg_xxzz_x_0[j] = pb_x * tg_xzz_x_0[j] + wp_x[j] * tg_xzz_x_1[j] + 0.5 * fl1_fx * tg_zz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_x_1[j] + 0.5 * fl1_fxn * tg_xzz_0_1[j];

                    tg_xxzz_y_0[j] = pb_x * tg_xzz_y_0[j] + wp_x[j] * tg_xzz_y_1[j] + 0.5 * fl1_fx * tg_zz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_y_1[j];

                    tg_xxzz_z_0[j] = pb_x * tg_xzz_z_0[j] + wp_x[j] * tg_xzz_z_1[j] + 0.5 * fl1_fx * tg_zz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_z_1[j];

                    tg_xyyy_x_0[j] = pb_x * tg_yyy_x_0[j] + wp_x[j] * tg_yyy_x_1[j] + 0.5 * fl1_fxn * tg_yyy_0_1[j];

                    tg_xyyy_y_0[j] = pb_x * tg_yyy_y_0[j] + wp_x[j] * tg_yyy_y_1[j];

                    tg_xyyy_z_0[j] = pb_x * tg_yyy_z_0[j] + wp_x[j] * tg_yyy_z_1[j];

                    tg_xyyz_x_0[j] = pb_x * tg_yyz_x_0[j] + wp_x[j] * tg_yyz_x_1[j] + 0.5 * fl1_fxn * tg_yyz_0_1[j];

                    tg_xyyz_y_0[j] = pb_x * tg_yyz_y_0[j] + wp_x[j] * tg_yyz_y_1[j];

                    tg_xyyz_z_0[j] = pb_x * tg_yyz_z_0[j] + wp_x[j] * tg_yyz_z_1[j];

                    tg_xyzz_x_0[j] = pb_x * tg_yzz_x_0[j] + wp_x[j] * tg_yzz_x_1[j] + 0.5 * fl1_fxn * tg_yzz_0_1[j];

                    tg_xyzz_y_0[j] = pb_x * tg_yzz_y_0[j] + wp_x[j] * tg_yzz_y_1[j];

                    tg_xyzz_z_0[j] = pb_x * tg_yzz_z_0[j] + wp_x[j] * tg_yzz_z_1[j];

                    tg_xzzz_x_0[j] = pb_x * tg_zzz_x_0[j] + wp_x[j] * tg_zzz_x_1[j] + 0.5 * fl1_fxn * tg_zzz_0_1[j];

                    tg_xzzz_y_0[j] = pb_x * tg_zzz_y_0[j] + wp_x[j] * tg_zzz_y_1[j];

                    tg_xzzz_z_0[j] = pb_x * tg_zzz_z_0[j] + wp_x[j] * tg_zzz_z_1[j];

                    tg_yyyy_x_0[j] = pb_y * tg_yyy_x_0[j] + wp_y[j] * tg_yyy_x_1[j] + 1.5 * fl1_fx * tg_yy_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_x_1[j];

                    tg_yyyy_y_0[j] = pb_y * tg_yyy_y_0[j] + wp_y[j] * tg_yyy_y_1[j] + 1.5 * fl1_fx * tg_yy_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_y_1[j] + 0.5 * fl1_fxn * tg_yyy_0_1[j];

                    tg_yyyy_z_0[j] = pb_y * tg_yyy_z_0[j] + wp_y[j] * tg_yyy_z_1[j] + 1.5 * fl1_fx * tg_yy_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yy_z_1[j];

                    tg_yyyz_x_0[j] = pb_y * tg_yyz_x_0[j] + wp_y[j] * tg_yyz_x_1[j] + fl1_fx * tg_yz_x_0[j] - fl1_fx * fl1_fza * tg_yz_x_1[j];

                    tg_yyyz_y_0[j] = pb_y * tg_yyz_y_0[j] + wp_y[j] * tg_yyz_y_1[j] + fl1_fx * tg_yz_y_0[j] - fl1_fx * fl1_fza * tg_yz_y_1[j] + 0.5 * fl1_fxn * tg_yyz_0_1[j];

                    tg_yyyz_z_0[j] = pb_y * tg_yyz_z_0[j] + wp_y[j] * tg_yyz_z_1[j] + fl1_fx * tg_yz_z_0[j] - fl1_fx * fl1_fza * tg_yz_z_1[j];

                    tg_yyzz_x_0[j] = pb_y * tg_yzz_x_0[j] + wp_y[j] * tg_yzz_x_1[j] + 0.5 * fl1_fx * tg_zz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_x_1[j];

                    tg_yyzz_y_0[j] = pb_y * tg_yzz_y_0[j] + wp_y[j] * tg_yzz_y_1[j] + 0.5 * fl1_fx * tg_zz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_y_1[j] + 0.5 * fl1_fxn * tg_yzz_0_1[j];

                    tg_yyzz_z_0[j] = pb_y * tg_yzz_z_0[j] + wp_y[j] * tg_yzz_z_1[j] + 0.5 * fl1_fx * tg_zz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zz_z_1[j];

                    tg_yzzz_x_0[j] = pb_y * tg_zzz_x_0[j] + wp_y[j] * tg_zzz_x_1[j];

                    tg_yzzz_y_0[j] = pb_y * tg_zzz_y_0[j] + wp_y[j] * tg_zzz_y_1[j] + 0.5 * fl1_fxn * tg_zzz_0_1[j];

                    tg_yzzz_z_0[j] = pb_y * tg_zzz_z_0[j] + wp_y[j] * tg_zzz_z_1[j];

                    tg_zzzz_x_0[j] = pb_z * tg_zzz_x_0[j] + wp_z[j] * tg_zzz_x_1[j] + 1.5 * fl1_fx * tg_zz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_x_1[j];

                    tg_zzzz_y_0[j] = pb_z * tg_zzz_y_0[j] + wp_z[j] * tg_zzz_y_1[j] + 1.5 * fl1_fx * tg_zz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_y_1[j];

                    tg_zzzz_z_0[j] = pb_z * tg_zzz_z_0[j] + wp_z[j] * tg_zzz_z_1[j] + 1.5 * fl1_fx * tg_zz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_zz_z_1[j] + 0.5 * fl1_fxn * tg_zzz_0_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSH(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_5_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

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

                auto tg_0_xxxxx_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx); 

                auto tg_0_xxxxy_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 1); 

                auto tg_0_xxxxz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 2); 

                auto tg_0_xxxyy_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 3); 

                auto tg_0_xxxyz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 4); 

                auto tg_0_xxxzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 5); 

                auto tg_0_xxyyy_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 6); 

                auto tg_0_xxyyz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 7); 

                auto tg_0_xxyzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 8); 

                auto tg_0_xxzzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 9); 

                auto tg_0_xyyyy_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 10); 

                auto tg_0_xyyyz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 11); 

                auto tg_0_xyyzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 12); 

                auto tg_0_xyzzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 13); 

                auto tg_0_xzzzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 14); 

                auto tg_0_yyyyy_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 15); 

                auto tg_0_yyyyz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 16); 

                auto tg_0_yyyzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 17); 

                auto tg_0_yyzzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 18); 

                auto tg_0_yzzzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 19); 

                auto tg_0_zzzzz_0 = primBuffer.data(pidx_g_0_5_m0 + 21 * idx + 20); 

                auto tg_0_xxxxx_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx); 

                auto tg_0_xxxxy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 1); 

                auto tg_0_xxxxz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 2); 

                auto tg_0_xxxyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 3); 

                auto tg_0_xxxyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 4); 

                auto tg_0_xxxzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 5); 

                auto tg_0_xxyyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 6); 

                auto tg_0_xxyyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 7); 

                auto tg_0_xxyzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 8); 

                auto tg_0_xxzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 9); 

                auto tg_0_xyyyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 10); 

                auto tg_0_xyyyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 11); 

                auto tg_0_xyyzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 12); 

                auto tg_0_xyzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 13); 

                auto tg_0_xzzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 14); 

                auto tg_0_yyyyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 15); 

                auto tg_0_yyyyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 16); 

                auto tg_0_yyyzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 17); 

                auto tg_0_yyzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 18); 

                auto tg_0_yzzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 19); 

                auto tg_0_zzzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 20); 

                auto tg_0_xxxx_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx); 

                auto tg_0_xxxy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 1); 

                auto tg_0_xxxz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 2); 

                auto tg_0_xxyy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 3); 

                auto tg_0_xxyz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 4); 

                auto tg_0_xxzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 5); 

                auto tg_0_xyyy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 6); 

                auto tg_0_xyyz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 7); 

                auto tg_0_xyzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 8); 

                auto tg_0_xzzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 9); 

                auto tg_0_yyyy_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 10); 

                auto tg_0_yyyz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 11); 

                auto tg_0_yyzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 12); 

                auto tg_0_yzzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 13); 

                auto tg_0_zzzz_1 = primBuffer.data(pidx_g_0_4_m1 + 15 * idx + 14); 

                // set up pointers to integrals

                auto tg_x_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx); 

                auto tg_x_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 1); 

                auto tg_x_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 2); 

                auto tg_x_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 3); 

                auto tg_x_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 4); 

                auto tg_x_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 5); 

                auto tg_x_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 6); 

                auto tg_x_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 7); 

                auto tg_x_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 8); 

                auto tg_x_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 9); 

                auto tg_x_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 10); 

                auto tg_x_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 11); 

                auto tg_x_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 12); 

                auto tg_x_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 13); 

                auto tg_x_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 14); 

                auto tg_x_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 15); 

                auto tg_x_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 16); 

                auto tg_x_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 17); 

                auto tg_x_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 18); 

                auto tg_x_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 19); 

                auto tg_x_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 20); 

                auto tg_y_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 21); 

                auto tg_y_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 22); 

                auto tg_y_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 23); 

                auto tg_y_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 24); 

                auto tg_y_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 25); 

                auto tg_y_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 26); 

                auto tg_y_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 27); 

                auto tg_y_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 28); 

                auto tg_y_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 29); 

                auto tg_y_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 30); 

                auto tg_y_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 31); 

                auto tg_y_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 32); 

                auto tg_y_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 33); 

                auto tg_y_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 34); 

                auto tg_y_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 35); 

                auto tg_y_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 36); 

                auto tg_y_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 37); 

                auto tg_y_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 38); 

                auto tg_y_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 39); 

                auto tg_y_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 40); 

                auto tg_y_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 41); 

                auto tg_z_xxxxx_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 42); 

                auto tg_z_xxxxy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 43); 

                auto tg_z_xxxxz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 44); 

                auto tg_z_xxxyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 45); 

                auto tg_z_xxxyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 46); 

                auto tg_z_xxxzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 47); 

                auto tg_z_xxyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 48); 

                auto tg_z_xxyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 49); 

                auto tg_z_xxyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 50); 

                auto tg_z_xxzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 51); 

                auto tg_z_xyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 52); 

                auto tg_z_xyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 53); 

                auto tg_z_xyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 54); 

                auto tg_z_xyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 55); 

                auto tg_z_xzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 56); 

                auto tg_z_yyyyy_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 57); 

                auto tg_z_yyyyz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 58); 

                auto tg_z_yyyzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 59); 

                auto tg_z_yyzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 60); 

                auto tg_z_yzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 61); 

                auto tg_z_zzzzz_0 = primBuffer.data(pidx_g_1_5_m0 + 63 * idx + 62); 

                #pragma omp simd aligned(fxn, tg_0_xxxx_1, tg_0_xxxxx_0, tg_0_xxxxx_1, tg_0_xxxxy_0, \
                                         tg_0_xxxxy_1, tg_0_xxxxz_0, tg_0_xxxxz_1, tg_0_xxxy_1, tg_0_xxxyy_0, tg_0_xxxyy_1, \
                                         tg_0_xxxyz_0, tg_0_xxxyz_1, tg_0_xxxz_1, tg_0_xxxzz_0, tg_0_xxxzz_1, tg_0_xxyy_1, \
                                         tg_0_xxyyy_0, tg_0_xxyyy_1, tg_0_xxyyz_0, tg_0_xxyyz_1, tg_0_xxyz_1, tg_0_xxyzz_0, \
                                         tg_0_xxyzz_1, tg_0_xxzz_1, tg_0_xxzzz_0, tg_0_xxzzz_1, tg_0_xyyy_1, tg_0_xyyyy_0, \
                                         tg_0_xyyyy_1, tg_0_xyyyz_0, tg_0_xyyyz_1, tg_0_xyyz_1, tg_0_xyyzz_0, tg_0_xyyzz_1, \
                                         tg_0_xyzz_1, tg_0_xyzzz_0, tg_0_xyzzz_1, tg_0_xzzz_1, tg_0_xzzzz_0, tg_0_xzzzz_1, \
                                         tg_0_yyyy_1, tg_0_yyyyy_0, tg_0_yyyyy_1, tg_0_yyyyz_0, tg_0_yyyyz_1, tg_0_yyyz_1, \
                                         tg_0_yyyzz_0, tg_0_yyyzz_1, tg_0_yyzz_1, tg_0_yyzzz_0, tg_0_yyzzz_1, tg_0_yzzz_1, \
                                         tg_0_yzzzz_0, tg_0_yzzzz_1, tg_0_zzzz_1, tg_0_zzzzz_0, tg_0_zzzzz_1, tg_x_xxxxx_0, \
                                         tg_x_xxxxy_0, tg_x_xxxxz_0, tg_x_xxxyy_0, tg_x_xxxyz_0, tg_x_xxxzz_0, tg_x_xxyyy_0, \
                                         tg_x_xxyyz_0, tg_x_xxyzz_0, tg_x_xxzzz_0, tg_x_xyyyy_0, tg_x_xyyyz_0, tg_x_xyyzz_0, \
                                         tg_x_xyzzz_0, tg_x_xzzzz_0, tg_x_yyyyy_0, tg_x_yyyyz_0, tg_x_yyyzz_0, tg_x_yyzzz_0, \
                                         tg_x_yzzzz_0, tg_x_zzzzz_0, tg_y_xxxxx_0, tg_y_xxxxy_0, tg_y_xxxxz_0, tg_y_xxxyy_0, \
                                         tg_y_xxxyz_0, tg_y_xxxzz_0, tg_y_xxyyy_0, tg_y_xxyyz_0, tg_y_xxyzz_0, tg_y_xxzzz_0, \
                                         tg_y_xyyyy_0, tg_y_xyyyz_0, tg_y_xyyzz_0, tg_y_xyzzz_0, tg_y_xzzzz_0, tg_y_yyyyy_0, \
                                         tg_y_yyyyz_0, tg_y_yyyzz_0, tg_y_yyzzz_0, tg_y_yzzzz_0, tg_y_zzzzz_0, tg_z_xxxxx_0, \
                                         tg_z_xxxxy_0, tg_z_xxxxz_0, tg_z_xxxyy_0, tg_z_xxxyz_0, tg_z_xxxzz_0, tg_z_xxyyy_0, \
                                         tg_z_xxyyz_0, tg_z_xxyzz_0, tg_z_xxzzz_0, tg_z_xyyyy_0, tg_z_xyyyz_0, tg_z_xyyzz_0, \
                                         tg_z_xyzzz_0, tg_z_xzzzz_0, tg_z_yyyyy_0, tg_z_yyyyz_0, tg_z_yyyzz_0, tg_z_yyzzz_0, \
                                         tg_z_yzzzz_0, tg_z_zzzzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_xxxxx_0[j] = pb_x * tg_0_xxxxx_0[j] + wp_x[j] * tg_0_xxxxx_1[j] + 2.5 * fl1_fxn * tg_0_xxxx_1[j];

                    tg_x_xxxxy_0[j] = pb_x * tg_0_xxxxy_0[j] + wp_x[j] * tg_0_xxxxy_1[j] + 2.0 * fl1_fxn * tg_0_xxxy_1[j];

                    tg_x_xxxxz_0[j] = pb_x * tg_0_xxxxz_0[j] + wp_x[j] * tg_0_xxxxz_1[j] + 2.0 * fl1_fxn * tg_0_xxxz_1[j];

                    tg_x_xxxyy_0[j] = pb_x * tg_0_xxxyy_0[j] + wp_x[j] * tg_0_xxxyy_1[j] + 1.5 * fl1_fxn * tg_0_xxyy_1[j];

                    tg_x_xxxyz_0[j] = pb_x * tg_0_xxxyz_0[j] + wp_x[j] * tg_0_xxxyz_1[j] + 1.5 * fl1_fxn * tg_0_xxyz_1[j];

                    tg_x_xxxzz_0[j] = pb_x * tg_0_xxxzz_0[j] + wp_x[j] * tg_0_xxxzz_1[j] + 1.5 * fl1_fxn * tg_0_xxzz_1[j];

                    tg_x_xxyyy_0[j] = pb_x * tg_0_xxyyy_0[j] + wp_x[j] * tg_0_xxyyy_1[j] + fl1_fxn * tg_0_xyyy_1[j];

                    tg_x_xxyyz_0[j] = pb_x * tg_0_xxyyz_0[j] + wp_x[j] * tg_0_xxyyz_1[j] + fl1_fxn * tg_0_xyyz_1[j];

                    tg_x_xxyzz_0[j] = pb_x * tg_0_xxyzz_0[j] + wp_x[j] * tg_0_xxyzz_1[j] + fl1_fxn * tg_0_xyzz_1[j];

                    tg_x_xxzzz_0[j] = pb_x * tg_0_xxzzz_0[j] + wp_x[j] * tg_0_xxzzz_1[j] + fl1_fxn * tg_0_xzzz_1[j];

                    tg_x_xyyyy_0[j] = pb_x * tg_0_xyyyy_0[j] + wp_x[j] * tg_0_xyyyy_1[j] + 0.5 * fl1_fxn * tg_0_yyyy_1[j];

                    tg_x_xyyyz_0[j] = pb_x * tg_0_xyyyz_0[j] + wp_x[j] * tg_0_xyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyz_1[j];

                    tg_x_xyyzz_0[j] = pb_x * tg_0_xyyzz_0[j] + wp_x[j] * tg_0_xyyzz_1[j] + 0.5 * fl1_fxn * tg_0_yyzz_1[j];

                    tg_x_xyzzz_0[j] = pb_x * tg_0_xyzzz_0[j] + wp_x[j] * tg_0_xyzzz_1[j] + 0.5 * fl1_fxn * tg_0_yzzz_1[j];

                    tg_x_xzzzz_0[j] = pb_x * tg_0_xzzzz_0[j] + wp_x[j] * tg_0_xzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzz_1[j];

                    tg_x_yyyyy_0[j] = pb_x * tg_0_yyyyy_0[j] + wp_x[j] * tg_0_yyyyy_1[j];

                    tg_x_yyyyz_0[j] = pb_x * tg_0_yyyyz_0[j] + wp_x[j] * tg_0_yyyyz_1[j];

                    tg_x_yyyzz_0[j] = pb_x * tg_0_yyyzz_0[j] + wp_x[j] * tg_0_yyyzz_1[j];

                    tg_x_yyzzz_0[j] = pb_x * tg_0_yyzzz_0[j] + wp_x[j] * tg_0_yyzzz_1[j];

                    tg_x_yzzzz_0[j] = pb_x * tg_0_yzzzz_0[j] + wp_x[j] * tg_0_yzzzz_1[j];

                    tg_x_zzzzz_0[j] = pb_x * tg_0_zzzzz_0[j] + wp_x[j] * tg_0_zzzzz_1[j];

                    tg_y_xxxxx_0[j] = pb_y * tg_0_xxxxx_0[j] + wp_y[j] * tg_0_xxxxx_1[j];

                    tg_y_xxxxy_0[j] = pb_y * tg_0_xxxxy_0[j] + wp_y[j] * tg_0_xxxxy_1[j] + 0.5 * fl1_fxn * tg_0_xxxx_1[j];

                    tg_y_xxxxz_0[j] = pb_y * tg_0_xxxxz_0[j] + wp_y[j] * tg_0_xxxxz_1[j];

                    tg_y_xxxyy_0[j] = pb_y * tg_0_xxxyy_0[j] + wp_y[j] * tg_0_xxxyy_1[j] + fl1_fxn * tg_0_xxxy_1[j];

                    tg_y_xxxyz_0[j] = pb_y * tg_0_xxxyz_0[j] + wp_y[j] * tg_0_xxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxz_1[j];

                    tg_y_xxxzz_0[j] = pb_y * tg_0_xxxzz_0[j] + wp_y[j] * tg_0_xxxzz_1[j];

                    tg_y_xxyyy_0[j] = pb_y * tg_0_xxyyy_0[j] + wp_y[j] * tg_0_xxyyy_1[j] + 1.5 * fl1_fxn * tg_0_xxyy_1[j];

                    tg_y_xxyyz_0[j] = pb_y * tg_0_xxyyz_0[j] + wp_y[j] * tg_0_xxyyz_1[j] + fl1_fxn * tg_0_xxyz_1[j];

                    tg_y_xxyzz_0[j] = pb_y * tg_0_xxyzz_0[j] + wp_y[j] * tg_0_xxyzz_1[j] + 0.5 * fl1_fxn * tg_0_xxzz_1[j];

                    tg_y_xxzzz_0[j] = pb_y * tg_0_xxzzz_0[j] + wp_y[j] * tg_0_xxzzz_1[j];

                    tg_y_xyyyy_0[j] = pb_y * tg_0_xyyyy_0[j] + wp_y[j] * tg_0_xyyyy_1[j] + 2.0 * fl1_fxn * tg_0_xyyy_1[j];

                    tg_y_xyyyz_0[j] = pb_y * tg_0_xyyyz_0[j] + wp_y[j] * tg_0_xyyyz_1[j] + 1.5 * fl1_fxn * tg_0_xyyz_1[j];

                    tg_y_xyyzz_0[j] = pb_y * tg_0_xyyzz_0[j] + wp_y[j] * tg_0_xyyzz_1[j] + fl1_fxn * tg_0_xyzz_1[j];

                    tg_y_xyzzz_0[j] = pb_y * tg_0_xyzzz_0[j] + wp_y[j] * tg_0_xyzzz_1[j] + 0.5 * fl1_fxn * tg_0_xzzz_1[j];

                    tg_y_xzzzz_0[j] = pb_y * tg_0_xzzzz_0[j] + wp_y[j] * tg_0_xzzzz_1[j];

                    tg_y_yyyyy_0[j] = pb_y * tg_0_yyyyy_0[j] + wp_y[j] * tg_0_yyyyy_1[j] + 2.5 * fl1_fxn * tg_0_yyyy_1[j];

                    tg_y_yyyyz_0[j] = pb_y * tg_0_yyyyz_0[j] + wp_y[j] * tg_0_yyyyz_1[j] + 2.0 * fl1_fxn * tg_0_yyyz_1[j];

                    tg_y_yyyzz_0[j] = pb_y * tg_0_yyyzz_0[j] + wp_y[j] * tg_0_yyyzz_1[j] + 1.5 * fl1_fxn * tg_0_yyzz_1[j];

                    tg_y_yyzzz_0[j] = pb_y * tg_0_yyzzz_0[j] + wp_y[j] * tg_0_yyzzz_1[j] + fl1_fxn * tg_0_yzzz_1[j];

                    tg_y_yzzzz_0[j] = pb_y * tg_0_yzzzz_0[j] + wp_y[j] * tg_0_yzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzz_1[j];

                    tg_y_zzzzz_0[j] = pb_y * tg_0_zzzzz_0[j] + wp_y[j] * tg_0_zzzzz_1[j];

                    tg_z_xxxxx_0[j] = pb_z * tg_0_xxxxx_0[j] + wp_z[j] * tg_0_xxxxx_1[j];

                    tg_z_xxxxy_0[j] = pb_z * tg_0_xxxxy_0[j] + wp_z[j] * tg_0_xxxxy_1[j];

                    tg_z_xxxxz_0[j] = pb_z * tg_0_xxxxz_0[j] + wp_z[j] * tg_0_xxxxz_1[j] + 0.5 * fl1_fxn * tg_0_xxxx_1[j];

                    tg_z_xxxyy_0[j] = pb_z * tg_0_xxxyy_0[j] + wp_z[j] * tg_0_xxxyy_1[j];

                    tg_z_xxxyz_0[j] = pb_z * tg_0_xxxyz_0[j] + wp_z[j] * tg_0_xxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxy_1[j];

                    tg_z_xxxzz_0[j] = pb_z * tg_0_xxxzz_0[j] + wp_z[j] * tg_0_xxxzz_1[j] + fl1_fxn * tg_0_xxxz_1[j];

                    tg_z_xxyyy_0[j] = pb_z * tg_0_xxyyy_0[j] + wp_z[j] * tg_0_xxyyy_1[j];

                    tg_z_xxyyz_0[j] = pb_z * tg_0_xxyyz_0[j] + wp_z[j] * tg_0_xxyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxyy_1[j];

                    tg_z_xxyzz_0[j] = pb_z * tg_0_xxyzz_0[j] + wp_z[j] * tg_0_xxyzz_1[j] + fl1_fxn * tg_0_xxyz_1[j];

                    tg_z_xxzzz_0[j] = pb_z * tg_0_xxzzz_0[j] + wp_z[j] * tg_0_xxzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxzz_1[j];

                    tg_z_xyyyy_0[j] = pb_z * tg_0_xyyyy_0[j] + wp_z[j] * tg_0_xyyyy_1[j];

                    tg_z_xyyyz_0[j] = pb_z * tg_0_xyyyz_0[j] + wp_z[j] * tg_0_xyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xyyy_1[j];

                    tg_z_xyyzz_0[j] = pb_z * tg_0_xyyzz_0[j] + wp_z[j] * tg_0_xyyzz_1[j] + fl1_fxn * tg_0_xyyz_1[j];

                    tg_z_xyzzz_0[j] = pb_z * tg_0_xyzzz_0[j] + wp_z[j] * tg_0_xyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xyzz_1[j];

                    tg_z_xzzzz_0[j] = pb_z * tg_0_xzzzz_0[j] + wp_z[j] * tg_0_xzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xzzz_1[j];

                    tg_z_yyyyy_0[j] = pb_z * tg_0_yyyyy_0[j] + wp_z[j] * tg_0_yyyyy_1[j];

                    tg_z_yyyyz_0[j] = pb_z * tg_0_yyyyz_0[j] + wp_z[j] * tg_0_yyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyy_1[j];

                    tg_z_yyyzz_0[j] = pb_z * tg_0_yyyzz_0[j] + wp_z[j] * tg_0_yyyzz_1[j] + fl1_fxn * tg_0_yyyz_1[j];

                    tg_z_yyzzz_0[j] = pb_z * tg_0_yyzzz_0[j] + wp_z[j] * tg_0_yyzzz_1[j] + 1.5 * fl1_fxn * tg_0_yyzz_1[j];

                    tg_z_yzzzz_0[j] = pb_z * tg_0_yzzzz_0[j] + wp_z[j] * tg_0_yzzzz_1[j] + 2.0 * fl1_fxn * tg_0_yzzz_1[j];

                    tg_z_zzzzz_0[j] = pb_z * tg_0_zzzzz_0[j] + wp_z[j] * tg_0_zzzzz_1[j] + 2.5 * fl1_fxn * tg_0_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
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
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xxxx_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx); 

                auto tg_xxxx_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 1); 

                auto tg_xxxx_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 2); 

                auto tg_xxxy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 3); 

                auto tg_xxxy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 4); 

                auto tg_xxxy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 5); 

                auto tg_xxxz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 6); 

                auto tg_xxxz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 7); 

                auto tg_xxxz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 8); 

                auto tg_xxyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 9); 

                auto tg_xxyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 10); 

                auto tg_xxyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 11); 

                auto tg_xxyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 12); 

                auto tg_xxyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 13); 

                auto tg_xxyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 14); 

                auto tg_xxzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 15); 

                auto tg_xxzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 16); 

                auto tg_xxzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 17); 

                auto tg_xyyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 18); 

                auto tg_xyyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 19); 

                auto tg_xyyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 20); 

                auto tg_xyyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 21); 

                auto tg_xyyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 22); 

                auto tg_xyyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 23); 

                auto tg_xyzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 24); 

                auto tg_xyzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 25); 

                auto tg_xyzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 26); 

                auto tg_xzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 27); 

                auto tg_xzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 28); 

                auto tg_xzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 29); 

                auto tg_yyyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 30); 

                auto tg_yyyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 31); 

                auto tg_yyyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 32); 

                auto tg_yyyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 33); 

                auto tg_yyyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 34); 

                auto tg_yyyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 35); 

                auto tg_yyzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 36); 

                auto tg_yyzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 37); 

                auto tg_yyzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 38); 

                auto tg_yzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 39); 

                auto tg_yzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 40); 

                auto tg_yzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 41); 

                auto tg_zzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 42); 

                auto tg_zzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 43); 

                auto tg_zzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 44); 

                auto tg_xxxx_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx); 

                auto tg_xxxx_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 1); 

                auto tg_xxxx_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 2); 

                auto tg_xxxy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 3); 

                auto tg_xxxy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 4); 

                auto tg_xxxy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 5); 

                auto tg_xxxz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 6); 

                auto tg_xxxz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 7); 

                auto tg_xxxz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 8); 

                auto tg_xxyy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 9); 

                auto tg_xxyy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 10); 

                auto tg_xxyy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 11); 

                auto tg_xxyz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 12); 

                auto tg_xxyz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 13); 

                auto tg_xxyz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 14); 

                auto tg_xxzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 15); 

                auto tg_xxzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 16); 

                auto tg_xxzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 17); 

                auto tg_xyyy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 18); 

                auto tg_xyyy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 19); 

                auto tg_xyyy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 20); 

                auto tg_xyyz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 21); 

                auto tg_xyyz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 22); 

                auto tg_xyyz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 23); 

                auto tg_xyzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 24); 

                auto tg_xyzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 25); 

                auto tg_xyzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 26); 

                auto tg_xzzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 27); 

                auto tg_xzzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 28); 

                auto tg_xzzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 29); 

                auto tg_yyyy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 30); 

                auto tg_yyyy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 31); 

                auto tg_yyyy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 32); 

                auto tg_yyyz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 33); 

                auto tg_yyyz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 34); 

                auto tg_yyyz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 35); 

                auto tg_yyzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 36); 

                auto tg_yyzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 37); 

                auto tg_yyzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 38); 

                auto tg_yzzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 39); 

                auto tg_yzzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 40); 

                auto tg_yzzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 41); 

                auto tg_zzzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 42); 

                auto tg_zzzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 43); 

                auto tg_zzzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 44); 

                auto tg_xxx_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx); 

                auto tg_xxx_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 1); 

                auto tg_xxx_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 2); 

                auto tg_xxy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 3); 

                auto tg_xxy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 4); 

                auto tg_xxy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 5); 

                auto tg_xxz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 6); 

                auto tg_xxz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 7); 

                auto tg_xxz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 8); 

                auto tg_xyy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 9); 

                auto tg_xyy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 10); 

                auto tg_xyy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 11); 

                auto tg_xyz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 12); 

                auto tg_xyz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 13); 

                auto tg_xyz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 14); 

                auto tg_xzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 15); 

                auto tg_xzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 16); 

                auto tg_xzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 17); 

                auto tg_yyy_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 18); 

                auto tg_yyy_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 19); 

                auto tg_yyy_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 20); 

                auto tg_yyz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 21); 

                auto tg_yyz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 22); 

                auto tg_yyz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 23); 

                auto tg_yzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 24); 

                auto tg_yzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 25); 

                auto tg_yzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 26); 

                auto tg_zzz_x_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 27); 

                auto tg_zzz_y_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 28); 

                auto tg_zzz_z_0 = primBuffer.data(pidx_g_3_1_m0 + 30 * idx + 29); 

                auto tg_xxx_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx); 

                auto tg_xxx_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 1); 

                auto tg_xxx_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 2); 

                auto tg_xxy_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 3); 

                auto tg_xxy_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 4); 

                auto tg_xxy_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 5); 

                auto tg_xxz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 6); 

                auto tg_xxz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 7); 

                auto tg_xxz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 8); 

                auto tg_xyy_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 9); 

                auto tg_xyy_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 10); 

                auto tg_xyy_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 11); 

                auto tg_xyz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 12); 

                auto tg_xyz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 13); 

                auto tg_xyz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 14); 

                auto tg_xzz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 15); 

                auto tg_xzz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 16); 

                auto tg_xzz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 17); 

                auto tg_yyy_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 18); 

                auto tg_yyy_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 19); 

                auto tg_yyy_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 20); 

                auto tg_yyz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 21); 

                auto tg_yyz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 22); 

                auto tg_yyz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 23); 

                auto tg_yzz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 24); 

                auto tg_yzz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 25); 

                auto tg_yzz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 26); 

                auto tg_zzz_x_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 27); 

                auto tg_zzz_y_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 28); 

                auto tg_zzz_z_1 = primBuffer.data(pidx_g_3_1_m1 + 30 * idx + 29); 

                auto tg_xxxx_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx); 

                auto tg_xxxy_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 1); 

                auto tg_xxxz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 2); 

                auto tg_xxyy_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 3); 

                auto tg_xxyz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 4); 

                auto tg_xxzz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 5); 

                auto tg_xyyy_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 6); 

                auto tg_xyyz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 7); 

                auto tg_xyzz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 8); 

                auto tg_xzzz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 9); 

                auto tg_yyyy_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 10); 

                auto tg_yyyz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 11); 

                auto tg_yyzz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 12); 

                auto tg_yzzz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 13); 

                auto tg_zzzz_0_1 = primBuffer.data(pidx_g_4_0_m1 + 15 * idx + 14); 

                // set up pointers to integrals

                auto tg_xxxxx_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx); 

                auto tg_xxxxx_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 1); 

                auto tg_xxxxx_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 2); 

                auto tg_xxxxy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 3); 

                auto tg_xxxxy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 4); 

                auto tg_xxxxy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 5); 

                auto tg_xxxxz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 6); 

                auto tg_xxxxz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 7); 

                auto tg_xxxxz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 8); 

                auto tg_xxxyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 9); 

                auto tg_xxxyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 10); 

                auto tg_xxxyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 11); 

                auto tg_xxxyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 12); 

                auto tg_xxxyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 13); 

                auto tg_xxxyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 14); 

                auto tg_xxxzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 15); 

                auto tg_xxxzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 16); 

                auto tg_xxxzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 17); 

                auto tg_xxyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 18); 

                auto tg_xxyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 19); 

                auto tg_xxyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 20); 

                auto tg_xxyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 21); 

                auto tg_xxyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 22); 

                auto tg_xxyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 23); 

                auto tg_xxyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 24); 

                auto tg_xxyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 25); 

                auto tg_xxyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 26); 

                auto tg_xxzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 27); 

                auto tg_xxzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 28); 

                auto tg_xxzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 29); 

                auto tg_xyyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 30); 

                auto tg_xyyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 31); 

                auto tg_xyyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 32); 

                auto tg_xyyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 33); 

                auto tg_xyyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 34); 

                auto tg_xyyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 35); 

                auto tg_xyyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 36); 

                auto tg_xyyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 37); 

                auto tg_xyyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 38); 

                auto tg_xyzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 39); 

                auto tg_xyzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 40); 

                auto tg_xyzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 41); 

                auto tg_xzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 42); 

                auto tg_xzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 43); 

                auto tg_xzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 44); 

                auto tg_yyyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 45); 

                auto tg_yyyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 46); 

                auto tg_yyyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 47); 

                auto tg_yyyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 48); 

                auto tg_yyyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 49); 

                auto tg_yyyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 50); 

                auto tg_yyyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 51); 

                auto tg_yyyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 52); 

                auto tg_yyyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 53); 

                auto tg_yyzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 54); 

                auto tg_yyzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 55); 

                auto tg_yyzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 56); 

                auto tg_yzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 57); 

                auto tg_yzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 58); 

                auto tg_yzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 59); 

                auto tg_zzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 60); 

                auto tg_zzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 61); 

                auto tg_zzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 62); 

                #pragma omp simd aligned(fxn, fza, tg_xxx_x_0, tg_xxx_x_1, tg_xxx_y_0, tg_xxx_y_1, tg_xxx_z_0, \
                                         tg_xxx_z_1, tg_xxxx_0_1, tg_xxxx_x_0, tg_xxxx_x_1, tg_xxxx_y_0, tg_xxxx_y_1, \
                                         tg_xxxx_z_0, tg_xxxx_z_1, tg_xxxxx_x_0, tg_xxxxx_y_0, tg_xxxxx_z_0, tg_xxxxy_x_0, \
                                         tg_xxxxy_y_0, tg_xxxxy_z_0, tg_xxxxz_x_0, tg_xxxxz_y_0, tg_xxxxz_z_0, tg_xxxy_0_1, \
                                         tg_xxxy_x_0, tg_xxxy_x_1, tg_xxxy_y_0, tg_xxxy_y_1, tg_xxxy_z_0, tg_xxxy_z_1, \
                                         tg_xxxyy_x_0, tg_xxxyy_y_0, tg_xxxyy_z_0, tg_xxxyz_x_0, tg_xxxyz_y_0, tg_xxxyz_z_0, \
                                         tg_xxxz_0_1, tg_xxxz_x_0, tg_xxxz_x_1, tg_xxxz_y_0, tg_xxxz_y_1, tg_xxxz_z_0, \
                                         tg_xxxz_z_1, tg_xxxzz_x_0, tg_xxxzz_y_0, tg_xxxzz_z_0, tg_xxy_x_0, tg_xxy_x_1, \
                                         tg_xxy_y_0, tg_xxy_y_1, tg_xxy_z_0, tg_xxy_z_1, tg_xxyy_0_1, tg_xxyy_x_0, \
                                         tg_xxyy_x_1, tg_xxyy_y_0, tg_xxyy_y_1, tg_xxyy_z_0, tg_xxyy_z_1, tg_xxyyy_x_0, \
                                         tg_xxyyy_y_0, tg_xxyyy_z_0, tg_xxyyz_x_0, tg_xxyyz_y_0, tg_xxyyz_z_0, tg_xxyz_0_1, \
                                         tg_xxyz_x_0, tg_xxyz_x_1, tg_xxyz_y_0, tg_xxyz_y_1, tg_xxyz_z_0, tg_xxyz_z_1, \
                                         tg_xxyzz_x_0, tg_xxyzz_y_0, tg_xxyzz_z_0, tg_xxz_x_0, tg_xxz_x_1, tg_xxz_y_0, \
                                         tg_xxz_y_1, tg_xxz_z_0, tg_xxz_z_1, tg_xxzz_0_1, tg_xxzz_x_0, tg_xxzz_x_1, \
                                         tg_xxzz_y_0, tg_xxzz_y_1, tg_xxzz_z_0, tg_xxzz_z_1, tg_xxzzz_x_0, tg_xxzzz_y_0, \
                                         tg_xxzzz_z_0, tg_xyy_x_0, tg_xyy_x_1, tg_xyy_y_0, tg_xyy_y_1, tg_xyy_z_0, tg_xyy_z_1, \
                                         tg_xyyy_0_1, tg_xyyy_x_0, tg_xyyy_x_1, tg_xyyy_y_0, tg_xyyy_y_1, tg_xyyy_z_0, \
                                         tg_xyyy_z_1, tg_xyyyy_x_0, tg_xyyyy_y_0, tg_xyyyy_z_0, tg_xyyyz_x_0, tg_xyyyz_y_0, \
                                         tg_xyyyz_z_0, tg_xyyz_0_1, tg_xyyz_x_0, tg_xyyz_x_1, tg_xyyz_y_0, tg_xyyz_y_1, \
                                         tg_xyyz_z_0, tg_xyyz_z_1, tg_xyyzz_x_0, tg_xyyzz_y_0, tg_xyyzz_z_0, tg_xyz_x_0, \
                                         tg_xyz_x_1, tg_xyz_y_0, tg_xyz_y_1, tg_xyz_z_0, tg_xyz_z_1, tg_xyzz_0_1, \
                                         tg_xyzz_x_0, tg_xyzz_x_1, tg_xyzz_y_0, tg_xyzz_y_1, tg_xyzz_z_0, tg_xyzz_z_1, \
                                         tg_xyzzz_x_0, tg_xyzzz_y_0, tg_xyzzz_z_0, tg_xzz_x_0, tg_xzz_x_1, tg_xzz_y_0, \
                                         tg_xzz_y_1, tg_xzz_z_0, tg_xzz_z_1, tg_xzzz_0_1, tg_xzzz_x_0, tg_xzzz_x_1, \
                                         tg_xzzz_y_0, tg_xzzz_y_1, tg_xzzz_z_0, tg_xzzz_z_1, tg_xzzzz_x_0, tg_xzzzz_y_0, \
                                         tg_xzzzz_z_0, tg_yyy_x_0, tg_yyy_x_1, tg_yyy_y_0, tg_yyy_y_1, tg_yyy_z_0, tg_yyy_z_1, \
                                         tg_yyyy_0_1, tg_yyyy_x_0, tg_yyyy_x_1, tg_yyyy_y_0, tg_yyyy_y_1, tg_yyyy_z_0, \
                                         tg_yyyy_z_1, tg_yyyyy_x_0, tg_yyyyy_y_0, tg_yyyyy_z_0, tg_yyyyz_x_0, tg_yyyyz_y_0, \
                                         tg_yyyyz_z_0, tg_yyyz_0_1, tg_yyyz_x_0, tg_yyyz_x_1, tg_yyyz_y_0, tg_yyyz_y_1, \
                                         tg_yyyz_z_0, tg_yyyz_z_1, tg_yyyzz_x_0, tg_yyyzz_y_0, tg_yyyzz_z_0, tg_yyz_x_0, \
                                         tg_yyz_x_1, tg_yyz_y_0, tg_yyz_y_1, tg_yyz_z_0, tg_yyz_z_1, tg_yyzz_0_1, \
                                         tg_yyzz_x_0, tg_yyzz_x_1, tg_yyzz_y_0, tg_yyzz_y_1, tg_yyzz_z_0, tg_yyzz_z_1, \
                                         tg_yyzzz_x_0, tg_yyzzz_y_0, tg_yyzzz_z_0, tg_yzz_x_0, tg_yzz_x_1, tg_yzz_y_0, \
                                         tg_yzz_y_1, tg_yzz_z_0, tg_yzz_z_1, tg_yzzz_0_1, tg_yzzz_x_0, tg_yzzz_x_1, \
                                         tg_yzzz_y_0, tg_yzzz_y_1, tg_yzzz_z_0, tg_yzzz_z_1, tg_yzzzz_x_0, tg_yzzzz_y_0, \
                                         tg_yzzzz_z_0, tg_zzz_x_0, tg_zzz_x_1, tg_zzz_y_0, tg_zzz_y_1, tg_zzz_z_0, tg_zzz_z_1, \
                                         tg_zzzz_0_1, tg_zzzz_x_0, tg_zzzz_x_1, tg_zzzz_y_0, tg_zzzz_y_1, tg_zzzz_z_0, \
                                         tg_zzzz_z_1, tg_zzzzz_x_0, tg_zzzzz_y_0, tg_zzzzz_z_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxx_x_0[j] = pb_x * tg_xxxx_x_0[j] + wp_x[j] * tg_xxxx_x_1[j] + 2.0 * fl1_fx * tg_xxx_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_x_1[j] + 0.5 * fl1_fxn * tg_xxxx_0_1[j];

                    tg_xxxxx_y_0[j] = pb_x * tg_xxxx_y_0[j] + wp_x[j] * tg_xxxx_y_1[j] + 2.0 * fl1_fx * tg_xxx_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_y_1[j];

                    tg_xxxxx_z_0[j] = pb_x * tg_xxxx_z_0[j] + wp_x[j] * tg_xxxx_z_1[j] + 2.0 * fl1_fx * tg_xxx_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxx_z_1[j];

                    tg_xxxxy_x_0[j] = pb_x * tg_xxxy_x_0[j] + wp_x[j] * tg_xxxy_x_1[j] + 1.5 * fl1_fx * tg_xxy_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_x_1[j] + 0.5 * fl1_fxn * tg_xxxy_0_1[j];

                    tg_xxxxy_y_0[j] = pb_x * tg_xxxy_y_0[j] + wp_x[j] * tg_xxxy_y_1[j] + 1.5 * fl1_fx * tg_xxy_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_y_1[j];

                    tg_xxxxy_z_0[j] = pb_x * tg_xxxy_z_0[j] + wp_x[j] * tg_xxxy_z_1[j] + 1.5 * fl1_fx * tg_xxy_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxy_z_1[j];

                    tg_xxxxz_x_0[j] = pb_x * tg_xxxz_x_0[j] + wp_x[j] * tg_xxxz_x_1[j] + 1.5 * fl1_fx * tg_xxz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_x_1[j] + 0.5 * fl1_fxn * tg_xxxz_0_1[j];

                    tg_xxxxz_y_0[j] = pb_x * tg_xxxz_y_0[j] + wp_x[j] * tg_xxxz_y_1[j] + 1.5 * fl1_fx * tg_xxz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_y_1[j];

                    tg_xxxxz_z_0[j] = pb_x * tg_xxxz_z_0[j] + wp_x[j] * tg_xxxz_z_1[j] + 1.5 * fl1_fx * tg_xxz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxz_z_1[j];

                    tg_xxxyy_x_0[j] = pb_x * tg_xxyy_x_0[j] + wp_x[j] * tg_xxyy_x_1[j] + fl1_fx * tg_xyy_x_0[j] - fl1_fx * fl1_fza * tg_xyy_x_1[j] + 0.5 * fl1_fxn * tg_xxyy_0_1[j];

                    tg_xxxyy_y_0[j] = pb_x * tg_xxyy_y_0[j] + wp_x[j] * tg_xxyy_y_1[j] + fl1_fx * tg_xyy_y_0[j] - fl1_fx * fl1_fza * tg_xyy_y_1[j];

                    tg_xxxyy_z_0[j] = pb_x * tg_xxyy_z_0[j] + wp_x[j] * tg_xxyy_z_1[j] + fl1_fx * tg_xyy_z_0[j] - fl1_fx * fl1_fza * tg_xyy_z_1[j];

                    tg_xxxyz_x_0[j] = pb_x * tg_xxyz_x_0[j] + wp_x[j] * tg_xxyz_x_1[j] + fl1_fx * tg_xyz_x_0[j] - fl1_fx * fl1_fza * tg_xyz_x_1[j] + 0.5 * fl1_fxn * tg_xxyz_0_1[j];

                    tg_xxxyz_y_0[j] = pb_x * tg_xxyz_y_0[j] + wp_x[j] * tg_xxyz_y_1[j] + fl1_fx * tg_xyz_y_0[j] - fl1_fx * fl1_fza * tg_xyz_y_1[j];

                    tg_xxxyz_z_0[j] = pb_x * tg_xxyz_z_0[j] + wp_x[j] * tg_xxyz_z_1[j] + fl1_fx * tg_xyz_z_0[j] - fl1_fx * fl1_fza * tg_xyz_z_1[j];

                    tg_xxxzz_x_0[j] = pb_x * tg_xxzz_x_0[j] + wp_x[j] * tg_xxzz_x_1[j] + fl1_fx * tg_xzz_x_0[j] - fl1_fx * fl1_fza * tg_xzz_x_1[j] + 0.5 * fl1_fxn * tg_xxzz_0_1[j];

                    tg_xxxzz_y_0[j] = pb_x * tg_xxzz_y_0[j] + wp_x[j] * tg_xxzz_y_1[j] + fl1_fx * tg_xzz_y_0[j] - fl1_fx * fl1_fza * tg_xzz_y_1[j];

                    tg_xxxzz_z_0[j] = pb_x * tg_xxzz_z_0[j] + wp_x[j] * tg_xxzz_z_1[j] + fl1_fx * tg_xzz_z_0[j] - fl1_fx * fl1_fza * tg_xzz_z_1[j];

                    tg_xxyyy_x_0[j] = pb_x * tg_xyyy_x_0[j] + wp_x[j] * tg_xyyy_x_1[j] + 0.5 * fl1_fx * tg_yyy_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_x_1[j] + 0.5 * fl1_fxn * tg_xyyy_0_1[j];

                    tg_xxyyy_y_0[j] = pb_x * tg_xyyy_y_0[j] + wp_x[j] * tg_xyyy_y_1[j] + 0.5 * fl1_fx * tg_yyy_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_y_1[j];

                    tg_xxyyy_z_0[j] = pb_x * tg_xyyy_z_0[j] + wp_x[j] * tg_xyyy_z_1[j] + 0.5 * fl1_fx * tg_yyy_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyy_z_1[j];

                    tg_xxyyz_x_0[j] = pb_x * tg_xyyz_x_0[j] + wp_x[j] * tg_xyyz_x_1[j] + 0.5 * fl1_fx * tg_yyz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_x_1[j] + 0.5 * fl1_fxn * tg_xyyz_0_1[j];

                    tg_xxyyz_y_0[j] = pb_x * tg_xyyz_y_0[j] + wp_x[j] * tg_xyyz_y_1[j] + 0.5 * fl1_fx * tg_yyz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_y_1[j];

                    tg_xxyyz_z_0[j] = pb_x * tg_xyyz_z_0[j] + wp_x[j] * tg_xyyz_z_1[j] + 0.5 * fl1_fx * tg_yyz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyz_z_1[j];

                    tg_xxyzz_x_0[j] = pb_x * tg_xyzz_x_0[j] + wp_x[j] * tg_xyzz_x_1[j] + 0.5 * fl1_fx * tg_yzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_x_1[j] + 0.5 * fl1_fxn * tg_xyzz_0_1[j];

                    tg_xxyzz_y_0[j] = pb_x * tg_xyzz_y_0[j] + wp_x[j] * tg_xyzz_y_1[j] + 0.5 * fl1_fx * tg_yzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_y_1[j];

                    tg_xxyzz_z_0[j] = pb_x * tg_xyzz_z_0[j] + wp_x[j] * tg_xyzz_z_1[j] + 0.5 * fl1_fx * tg_yzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzz_z_1[j];

                    tg_xxzzz_x_0[j] = pb_x * tg_xzzz_x_0[j] + wp_x[j] * tg_xzzz_x_1[j] + 0.5 * fl1_fx * tg_zzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_x_1[j] + 0.5 * fl1_fxn * tg_xzzz_0_1[j];

                    tg_xxzzz_y_0[j] = pb_x * tg_xzzz_y_0[j] + wp_x[j] * tg_xzzz_y_1[j] + 0.5 * fl1_fx * tg_zzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_y_1[j];

                    tg_xxzzz_z_0[j] = pb_x * tg_xzzz_z_0[j] + wp_x[j] * tg_xzzz_z_1[j] + 0.5 * fl1_fx * tg_zzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_z_1[j];

                    tg_xyyyy_x_0[j] = pb_x * tg_yyyy_x_0[j] + wp_x[j] * tg_yyyy_x_1[j] + 0.5 * fl1_fxn * tg_yyyy_0_1[j];

                    tg_xyyyy_y_0[j] = pb_x * tg_yyyy_y_0[j] + wp_x[j] * tg_yyyy_y_1[j];

                    tg_xyyyy_z_0[j] = pb_x * tg_yyyy_z_0[j] + wp_x[j] * tg_yyyy_z_1[j];

                    tg_xyyyz_x_0[j] = pb_x * tg_yyyz_x_0[j] + wp_x[j] * tg_yyyz_x_1[j] + 0.5 * fl1_fxn * tg_yyyz_0_1[j];

                    tg_xyyyz_y_0[j] = pb_x * tg_yyyz_y_0[j] + wp_x[j] * tg_yyyz_y_1[j];

                    tg_xyyyz_z_0[j] = pb_x * tg_yyyz_z_0[j] + wp_x[j] * tg_yyyz_z_1[j];

                    tg_xyyzz_x_0[j] = pb_x * tg_yyzz_x_0[j] + wp_x[j] * tg_yyzz_x_1[j] + 0.5 * fl1_fxn * tg_yyzz_0_1[j];

                    tg_xyyzz_y_0[j] = pb_x * tg_yyzz_y_0[j] + wp_x[j] * tg_yyzz_y_1[j];

                    tg_xyyzz_z_0[j] = pb_x * tg_yyzz_z_0[j] + wp_x[j] * tg_yyzz_z_1[j];

                    tg_xyzzz_x_0[j] = pb_x * tg_yzzz_x_0[j] + wp_x[j] * tg_yzzz_x_1[j] + 0.5 * fl1_fxn * tg_yzzz_0_1[j];

                    tg_xyzzz_y_0[j] = pb_x * tg_yzzz_y_0[j] + wp_x[j] * tg_yzzz_y_1[j];

                    tg_xyzzz_z_0[j] = pb_x * tg_yzzz_z_0[j] + wp_x[j] * tg_yzzz_z_1[j];

                    tg_xzzzz_x_0[j] = pb_x * tg_zzzz_x_0[j] + wp_x[j] * tg_zzzz_x_1[j] + 0.5 * fl1_fxn * tg_zzzz_0_1[j];

                    tg_xzzzz_y_0[j] = pb_x * tg_zzzz_y_0[j] + wp_x[j] * tg_zzzz_y_1[j];

                    tg_xzzzz_z_0[j] = pb_x * tg_zzzz_z_0[j] + wp_x[j] * tg_zzzz_z_1[j];

                    tg_yyyyy_x_0[j] = pb_y * tg_yyyy_x_0[j] + wp_y[j] * tg_yyyy_x_1[j] + 2.0 * fl1_fx * tg_yyy_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_x_1[j];

                    tg_yyyyy_y_0[j] = pb_y * tg_yyyy_y_0[j] + wp_y[j] * tg_yyyy_y_1[j] + 2.0 * fl1_fx * tg_yyy_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_y_1[j] + 0.5 * fl1_fxn * tg_yyyy_0_1[j];

                    tg_yyyyy_z_0[j] = pb_y * tg_yyyy_z_0[j] + wp_y[j] * tg_yyyy_z_1[j] + 2.0 * fl1_fx * tg_yyy_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyy_z_1[j];

                    tg_yyyyz_x_0[j] = pb_y * tg_yyyz_x_0[j] + wp_y[j] * tg_yyyz_x_1[j] + 1.5 * fl1_fx * tg_yyz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_x_1[j];

                    tg_yyyyz_y_0[j] = pb_y * tg_yyyz_y_0[j] + wp_y[j] * tg_yyyz_y_1[j] + 1.5 * fl1_fx * tg_yyz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_y_1[j] + 0.5 * fl1_fxn * tg_yyyz_0_1[j];

                    tg_yyyyz_z_0[j] = pb_y * tg_yyyz_z_0[j] + wp_y[j] * tg_yyyz_z_1[j] + 1.5 * fl1_fx * tg_yyz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyz_z_1[j];

                    tg_yyyzz_x_0[j] = pb_y * tg_yyzz_x_0[j] + wp_y[j] * tg_yyzz_x_1[j] + fl1_fx * tg_yzz_x_0[j] - fl1_fx * fl1_fza * tg_yzz_x_1[j];

                    tg_yyyzz_y_0[j] = pb_y * tg_yyzz_y_0[j] + wp_y[j] * tg_yyzz_y_1[j] + fl1_fx * tg_yzz_y_0[j] - fl1_fx * fl1_fza * tg_yzz_y_1[j] + 0.5 * fl1_fxn * tg_yyzz_0_1[j];

                    tg_yyyzz_z_0[j] = pb_y * tg_yyzz_z_0[j] + wp_y[j] * tg_yyzz_z_1[j] + fl1_fx * tg_yzz_z_0[j] - fl1_fx * fl1_fza * tg_yzz_z_1[j];

                    tg_yyzzz_x_0[j] = pb_y * tg_yzzz_x_0[j] + wp_y[j] * tg_yzzz_x_1[j] + 0.5 * fl1_fx * tg_zzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_x_1[j];

                    tg_yyzzz_y_0[j] = pb_y * tg_yzzz_y_0[j] + wp_y[j] * tg_yzzz_y_1[j] + 0.5 * fl1_fx * tg_zzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_y_1[j] + 0.5 * fl1_fxn * tg_yzzz_0_1[j];

                    tg_yyzzz_z_0[j] = pb_y * tg_yzzz_z_0[j] + wp_y[j] * tg_yzzz_z_1[j] + 0.5 * fl1_fx * tg_zzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzz_z_1[j];

                    tg_yzzzz_x_0[j] = pb_y * tg_zzzz_x_0[j] + wp_y[j] * tg_zzzz_x_1[j];

                    tg_yzzzz_y_0[j] = pb_y * tg_zzzz_y_0[j] + wp_y[j] * tg_zzzz_y_1[j] + 0.5 * fl1_fxn * tg_zzzz_0_1[j];

                    tg_yzzzz_z_0[j] = pb_y * tg_zzzz_z_0[j] + wp_y[j] * tg_zzzz_z_1[j];

                    tg_zzzzz_x_0[j] = pb_z * tg_zzzz_x_0[j] + wp_z[j] * tg_zzzz_x_1[j] + 2.0 * fl1_fx * tg_zzz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_x_1[j];

                    tg_zzzzz_y_0[j] = pb_z * tg_zzzz_y_0[j] + wp_z[j] * tg_zzzz_y_1[j] + 2.0 * fl1_fx * tg_zzz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_y_1[j];

                    tg_zzzzz_z_0[j] = pb_z * tg_zzzz_z_0[j] + wp_z[j] * tg_zzzz_z_1[j] + 2.0 * fl1_fx * tg_zzz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_zzz_z_1[j] + 0.5 * fl1_fxn * tg_zzzz_0_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSI(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {6, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_6_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {6, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {6, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_5_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {5, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

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

                auto tg_0_xxxxxx_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx); 

                auto tg_0_xxxxxy_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 1); 

                auto tg_0_xxxxxz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 2); 

                auto tg_0_xxxxyy_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 3); 

                auto tg_0_xxxxyz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 4); 

                auto tg_0_xxxxzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 5); 

                auto tg_0_xxxyyy_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 6); 

                auto tg_0_xxxyyz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 7); 

                auto tg_0_xxxyzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 8); 

                auto tg_0_xxxzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 9); 

                auto tg_0_xxyyyy_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 10); 

                auto tg_0_xxyyyz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 11); 

                auto tg_0_xxyyzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 12); 

                auto tg_0_xxyzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 13); 

                auto tg_0_xxzzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 14); 

                auto tg_0_xyyyyy_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 15); 

                auto tg_0_xyyyyz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 16); 

                auto tg_0_xyyyzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 17); 

                auto tg_0_xyyzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 18); 

                auto tg_0_xyzzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 19); 

                auto tg_0_xzzzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 20); 

                auto tg_0_yyyyyy_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 21); 

                auto tg_0_yyyyyz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 22); 

                auto tg_0_yyyyzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 23); 

                auto tg_0_yyyzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 24); 

                auto tg_0_yyzzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 25); 

                auto tg_0_yzzzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 26); 

                auto tg_0_zzzzzz_0 = primBuffer.data(pidx_g_0_6_m0 + 28 * idx + 27); 

                auto tg_0_xxxxxx_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx); 

                auto tg_0_xxxxxy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 1); 

                auto tg_0_xxxxxz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 2); 

                auto tg_0_xxxxyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 3); 

                auto tg_0_xxxxyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 4); 

                auto tg_0_xxxxzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 5); 

                auto tg_0_xxxyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 6); 

                auto tg_0_xxxyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 7); 

                auto tg_0_xxxyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 8); 

                auto tg_0_xxxzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 9); 

                auto tg_0_xxyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 10); 

                auto tg_0_xxyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 11); 

                auto tg_0_xxyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 12); 

                auto tg_0_xxyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 13); 

                auto tg_0_xxzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 14); 

                auto tg_0_xyyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 15); 

                auto tg_0_xyyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 16); 

                auto tg_0_xyyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 17); 

                auto tg_0_xyyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 18); 

                auto tg_0_xyzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 19); 

                auto tg_0_xzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 20); 

                auto tg_0_yyyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 21); 

                auto tg_0_yyyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 22); 

                auto tg_0_yyyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 23); 

                auto tg_0_yyyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 24); 

                auto tg_0_yyzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 25); 

                auto tg_0_yzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 26); 

                auto tg_0_zzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 27); 

                auto tg_0_xxxxx_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx); 

                auto tg_0_xxxxy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 1); 

                auto tg_0_xxxxz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 2); 

                auto tg_0_xxxyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 3); 

                auto tg_0_xxxyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 4); 

                auto tg_0_xxxzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 5); 

                auto tg_0_xxyyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 6); 

                auto tg_0_xxyyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 7); 

                auto tg_0_xxyzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 8); 

                auto tg_0_xxzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 9); 

                auto tg_0_xyyyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 10); 

                auto tg_0_xyyyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 11); 

                auto tg_0_xyyzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 12); 

                auto tg_0_xyzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 13); 

                auto tg_0_xzzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 14); 

                auto tg_0_yyyyy_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 15); 

                auto tg_0_yyyyz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 16); 

                auto tg_0_yyyzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 17); 

                auto tg_0_yyzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 18); 

                auto tg_0_yzzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 19); 

                auto tg_0_zzzzz_1 = primBuffer.data(pidx_g_0_5_m1 + 21 * idx + 20); 

                // set up pointers to integrals

                auto tg_x_xxxxxx_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx); 

                auto tg_x_xxxxxy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 1); 

                auto tg_x_xxxxxz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 2); 

                auto tg_x_xxxxyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 3); 

                auto tg_x_xxxxyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 4); 

                auto tg_x_xxxxzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 5); 

                auto tg_x_xxxyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 6); 

                auto tg_x_xxxyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 7); 

                auto tg_x_xxxyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 8); 

                auto tg_x_xxxzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 9); 

                auto tg_x_xxyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 10); 

                auto tg_x_xxyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 11); 

                auto tg_x_xxyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 12); 

                auto tg_x_xxyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 13); 

                auto tg_x_xxzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 14); 

                auto tg_x_xyyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 15); 

                auto tg_x_xyyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 16); 

                auto tg_x_xyyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 17); 

                auto tg_x_xyyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 18); 

                auto tg_x_xyzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 19); 

                auto tg_x_xzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 20); 

                auto tg_x_yyyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 21); 

                auto tg_x_yyyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 22); 

                auto tg_x_yyyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 23); 

                auto tg_x_yyyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 24); 

                auto tg_x_yyzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 25); 

                auto tg_x_yzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 26); 

                auto tg_x_zzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 27); 

                auto tg_y_xxxxxx_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 28); 

                auto tg_y_xxxxxy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 29); 

                auto tg_y_xxxxxz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 30); 

                auto tg_y_xxxxyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 31); 

                auto tg_y_xxxxyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 32); 

                auto tg_y_xxxxzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 33); 

                auto tg_y_xxxyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 34); 

                auto tg_y_xxxyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 35); 

                auto tg_y_xxxyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 36); 

                auto tg_y_xxxzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 37); 

                auto tg_y_xxyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 38); 

                auto tg_y_xxyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 39); 

                auto tg_y_xxyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 40); 

                auto tg_y_xxyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 41); 

                auto tg_y_xxzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 42); 

                auto tg_y_xyyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 43); 

                auto tg_y_xyyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 44); 

                auto tg_y_xyyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 45); 

                auto tg_y_xyyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 46); 

                auto tg_y_xyzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 47); 

                auto tg_y_xzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 48); 

                auto tg_y_yyyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 49); 

                auto tg_y_yyyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 50); 

                auto tg_y_yyyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 51); 

                auto tg_y_yyyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 52); 

                auto tg_y_yyzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 53); 

                auto tg_y_yzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 54); 

                auto tg_y_zzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 55); 

                auto tg_z_xxxxxx_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 56); 

                auto tg_z_xxxxxy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 57); 

                auto tg_z_xxxxxz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 58); 

                auto tg_z_xxxxyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 59); 

                auto tg_z_xxxxyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 60); 

                auto tg_z_xxxxzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 61); 

                auto tg_z_xxxyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 62); 

                auto tg_z_xxxyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 63); 

                auto tg_z_xxxyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 64); 

                auto tg_z_xxxzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 65); 

                auto tg_z_xxyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 66); 

                auto tg_z_xxyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 67); 

                auto tg_z_xxyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 68); 

                auto tg_z_xxyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 69); 

                auto tg_z_xxzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 70); 

                auto tg_z_xyyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 71); 

                auto tg_z_xyyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 72); 

                auto tg_z_xyyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 73); 

                auto tg_z_xyyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 74); 

                auto tg_z_xyzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 75); 

                auto tg_z_xzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 76); 

                auto tg_z_yyyyyy_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 77); 

                auto tg_z_yyyyyz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 78); 

                auto tg_z_yyyyzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 79); 

                auto tg_z_yyyzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 80); 

                auto tg_z_yyzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 81); 

                auto tg_z_yzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 82); 

                auto tg_z_zzzzzz_0 = primBuffer.data(pidx_g_1_6_m0 + 84 * idx + 83); 

                #pragma omp simd aligned(fxn, tg_0_xxxxx_1, tg_0_xxxxxx_0, tg_0_xxxxxx_1, tg_0_xxxxxy_0, \
                                         tg_0_xxxxxy_1, tg_0_xxxxxz_0, tg_0_xxxxxz_1, tg_0_xxxxy_1, tg_0_xxxxyy_0, \
                                         tg_0_xxxxyy_1, tg_0_xxxxyz_0, tg_0_xxxxyz_1, tg_0_xxxxz_1, tg_0_xxxxzz_0, \
                                         tg_0_xxxxzz_1, tg_0_xxxyy_1, tg_0_xxxyyy_0, tg_0_xxxyyy_1, tg_0_xxxyyz_0, \
                                         tg_0_xxxyyz_1, tg_0_xxxyz_1, tg_0_xxxyzz_0, tg_0_xxxyzz_1, tg_0_xxxzz_1, \
                                         tg_0_xxxzzz_0, tg_0_xxxzzz_1, tg_0_xxyyy_1, tg_0_xxyyyy_0, tg_0_xxyyyy_1, \
                                         tg_0_xxyyyz_0, tg_0_xxyyyz_1, tg_0_xxyyz_1, tg_0_xxyyzz_0, tg_0_xxyyzz_1, \
                                         tg_0_xxyzz_1, tg_0_xxyzzz_0, tg_0_xxyzzz_1, tg_0_xxzzz_1, tg_0_xxzzzz_0, \
                                         tg_0_xxzzzz_1, tg_0_xyyyy_1, tg_0_xyyyyy_0, tg_0_xyyyyy_1, tg_0_xyyyyz_0, \
                                         tg_0_xyyyyz_1, tg_0_xyyyz_1, tg_0_xyyyzz_0, tg_0_xyyyzz_1, tg_0_xyyzz_1, \
                                         tg_0_xyyzzz_0, tg_0_xyyzzz_1, tg_0_xyzzz_1, tg_0_xyzzzz_0, tg_0_xyzzzz_1, \
                                         tg_0_xzzzz_1, tg_0_xzzzzz_0, tg_0_xzzzzz_1, tg_0_yyyyy_1, tg_0_yyyyyy_0, \
                                         tg_0_yyyyyy_1, tg_0_yyyyyz_0, tg_0_yyyyyz_1, tg_0_yyyyz_1, tg_0_yyyyzz_0, \
                                         tg_0_yyyyzz_1, tg_0_yyyzz_1, tg_0_yyyzzz_0, tg_0_yyyzzz_1, tg_0_yyzzz_1, \
                                         tg_0_yyzzzz_0, tg_0_yyzzzz_1, tg_0_yzzzz_1, tg_0_yzzzzz_0, tg_0_yzzzzz_1, \
                                         tg_0_zzzzz_1, tg_0_zzzzzz_0, tg_0_zzzzzz_1, tg_x_xxxxxx_0, tg_x_xxxxxy_0, \
                                         tg_x_xxxxxz_0, tg_x_xxxxyy_0, tg_x_xxxxyz_0, tg_x_xxxxzz_0, tg_x_xxxyyy_0, \
                                         tg_x_xxxyyz_0, tg_x_xxxyzz_0, tg_x_xxxzzz_0, tg_x_xxyyyy_0, tg_x_xxyyyz_0, \
                                         tg_x_xxyyzz_0, tg_x_xxyzzz_0, tg_x_xxzzzz_0, tg_x_xyyyyy_0, tg_x_xyyyyz_0, \
                                         tg_x_xyyyzz_0, tg_x_xyyzzz_0, tg_x_xyzzzz_0, tg_x_xzzzzz_0, tg_x_yyyyyy_0, \
                                         tg_x_yyyyyz_0, tg_x_yyyyzz_0, tg_x_yyyzzz_0, tg_x_yyzzzz_0, tg_x_yzzzzz_0, \
                                         tg_x_zzzzzz_0, tg_y_xxxxxx_0, tg_y_xxxxxy_0, tg_y_xxxxxz_0, tg_y_xxxxyy_0, \
                                         tg_y_xxxxyz_0, tg_y_xxxxzz_0, tg_y_xxxyyy_0, tg_y_xxxyyz_0, tg_y_xxxyzz_0, \
                                         tg_y_xxxzzz_0, tg_y_xxyyyy_0, tg_y_xxyyyz_0, tg_y_xxyyzz_0, tg_y_xxyzzz_0, \
                                         tg_y_xxzzzz_0, tg_y_xyyyyy_0, tg_y_xyyyyz_0, tg_y_xyyyzz_0, tg_y_xyyzzz_0, \
                                         tg_y_xyzzzz_0, tg_y_xzzzzz_0, tg_y_yyyyyy_0, tg_y_yyyyyz_0, tg_y_yyyyzz_0, \
                                         tg_y_yyyzzz_0, tg_y_yyzzzz_0, tg_y_yzzzzz_0, tg_y_zzzzzz_0, tg_z_xxxxxx_0, \
                                         tg_z_xxxxxy_0, tg_z_xxxxxz_0, tg_z_xxxxyy_0, tg_z_xxxxyz_0, tg_z_xxxxzz_0, \
                                         tg_z_xxxyyy_0, tg_z_xxxyyz_0, tg_z_xxxyzz_0, tg_z_xxxzzz_0, tg_z_xxyyyy_0, \
                                         tg_z_xxyyyz_0, tg_z_xxyyzz_0, tg_z_xxyzzz_0, tg_z_xxzzzz_0, tg_z_xyyyyy_0, \
                                         tg_z_xyyyyz_0, tg_z_xyyyzz_0, tg_z_xyyzzz_0, tg_z_xyzzzz_0, tg_z_xzzzzz_0, \
                                         tg_z_yyyyyy_0, tg_z_yyyyyz_0, tg_z_yyyyzz_0, tg_z_yyyzzz_0, tg_z_yyzzzz_0, \
                                         tg_z_yzzzzz_0, tg_z_zzzzzz_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_xxxxxx_0[j] = pb_x * tg_0_xxxxxx_0[j] + wp_x[j] * tg_0_xxxxxx_1[j] + 3.0 * fl1_fxn * tg_0_xxxxx_1[j];

                    tg_x_xxxxxy_0[j] = pb_x * tg_0_xxxxxy_0[j] + wp_x[j] * tg_0_xxxxxy_1[j] + 2.5 * fl1_fxn * tg_0_xxxxy_1[j];

                    tg_x_xxxxxz_0[j] = pb_x * tg_0_xxxxxz_0[j] + wp_x[j] * tg_0_xxxxxz_1[j] + 2.5 * fl1_fxn * tg_0_xxxxz_1[j];

                    tg_x_xxxxyy_0[j] = pb_x * tg_0_xxxxyy_0[j] + wp_x[j] * tg_0_xxxxyy_1[j] + 2.0 * fl1_fxn * tg_0_xxxyy_1[j];

                    tg_x_xxxxyz_0[j] = pb_x * tg_0_xxxxyz_0[j] + wp_x[j] * tg_0_xxxxyz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyz_1[j];

                    tg_x_xxxxzz_0[j] = pb_x * tg_0_xxxxzz_0[j] + wp_x[j] * tg_0_xxxxzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxzz_1[j];

                    tg_x_xxxyyy_0[j] = pb_x * tg_0_xxxyyy_0[j] + wp_x[j] * tg_0_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_0_xxyyy_1[j];

                    tg_x_xxxyyz_0[j] = pb_x * tg_0_xxxyyz_0[j] + wp_x[j] * tg_0_xxxyyz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyz_1[j];

                    tg_x_xxxyzz_0[j] = pb_x * tg_0_xxxyzz_0[j] + wp_x[j] * tg_0_xxxyzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyzz_1[j];

                    tg_x_xxxzzz_0[j] = pb_x * tg_0_xxxzzz_0[j] + wp_x[j] * tg_0_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxzzz_1[j];

                    tg_x_xxyyyy_0[j] = pb_x * tg_0_xxyyyy_0[j] + wp_x[j] * tg_0_xxyyyy_1[j] + fl1_fxn * tg_0_xyyyy_1[j];

                    tg_x_xxyyyz_0[j] = pb_x * tg_0_xxyyyz_0[j] + wp_x[j] * tg_0_xxyyyz_1[j] + fl1_fxn * tg_0_xyyyz_1[j];

                    tg_x_xxyyzz_0[j] = pb_x * tg_0_xxyyzz_0[j] + wp_x[j] * tg_0_xxyyzz_1[j] + fl1_fxn * tg_0_xyyzz_1[j];

                    tg_x_xxyzzz_0[j] = pb_x * tg_0_xxyzzz_0[j] + wp_x[j] * tg_0_xxyzzz_1[j] + fl1_fxn * tg_0_xyzzz_1[j];

                    tg_x_xxzzzz_0[j] = pb_x * tg_0_xxzzzz_0[j] + wp_x[j] * tg_0_xxzzzz_1[j] + fl1_fxn * tg_0_xzzzz_1[j];

                    tg_x_xyyyyy_0[j] = pb_x * tg_0_xyyyyy_0[j] + wp_x[j] * tg_0_xyyyyy_1[j] + 0.5 * fl1_fxn * tg_0_yyyyy_1[j];

                    tg_x_xyyyyz_0[j] = pb_x * tg_0_xyyyyz_0[j] + wp_x[j] * tg_0_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyz_1[j];

                    tg_x_xyyyzz_0[j] = pb_x * tg_0_xyyyzz_0[j] + wp_x[j] * tg_0_xyyyzz_1[j] + 0.5 * fl1_fxn * tg_0_yyyzz_1[j];

                    tg_x_xyyzzz_0[j] = pb_x * tg_0_xyyzzz_0[j] + wp_x[j] * tg_0_xyyzzz_1[j] + 0.5 * fl1_fxn * tg_0_yyzzz_1[j];

                    tg_x_xyzzzz_0[j] = pb_x * tg_0_xyzzzz_0[j] + wp_x[j] * tg_0_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_0_yzzzz_1[j];

                    tg_x_xzzzzz_0[j] = pb_x * tg_0_xzzzzz_0[j] + wp_x[j] * tg_0_xzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzzz_1[j];

                    tg_x_yyyyyy_0[j] = pb_x * tg_0_yyyyyy_0[j] + wp_x[j] * tg_0_yyyyyy_1[j];

                    tg_x_yyyyyz_0[j] = pb_x * tg_0_yyyyyz_0[j] + wp_x[j] * tg_0_yyyyyz_1[j];

                    tg_x_yyyyzz_0[j] = pb_x * tg_0_yyyyzz_0[j] + wp_x[j] * tg_0_yyyyzz_1[j];

                    tg_x_yyyzzz_0[j] = pb_x * tg_0_yyyzzz_0[j] + wp_x[j] * tg_0_yyyzzz_1[j];

                    tg_x_yyzzzz_0[j] = pb_x * tg_0_yyzzzz_0[j] + wp_x[j] * tg_0_yyzzzz_1[j];

                    tg_x_yzzzzz_0[j] = pb_x * tg_0_yzzzzz_0[j] + wp_x[j] * tg_0_yzzzzz_1[j];

                    tg_x_zzzzzz_0[j] = pb_x * tg_0_zzzzzz_0[j] + wp_x[j] * tg_0_zzzzzz_1[j];

                    tg_y_xxxxxx_0[j] = pb_y * tg_0_xxxxxx_0[j] + wp_y[j] * tg_0_xxxxxx_1[j];

                    tg_y_xxxxxy_0[j] = pb_y * tg_0_xxxxxy_0[j] + wp_y[j] * tg_0_xxxxxy_1[j] + 0.5 * fl1_fxn * tg_0_xxxxx_1[j];

                    tg_y_xxxxxz_0[j] = pb_y * tg_0_xxxxxz_0[j] + wp_y[j] * tg_0_xxxxxz_1[j];

                    tg_y_xxxxyy_0[j] = pb_y * tg_0_xxxxyy_0[j] + wp_y[j] * tg_0_xxxxyy_1[j] + fl1_fxn * tg_0_xxxxy_1[j];

                    tg_y_xxxxyz_0[j] = pb_y * tg_0_xxxxyz_0[j] + wp_y[j] * tg_0_xxxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxz_1[j];

                    tg_y_xxxxzz_0[j] = pb_y * tg_0_xxxxzz_0[j] + wp_y[j] * tg_0_xxxxzz_1[j];

                    tg_y_xxxyyy_0[j] = pb_y * tg_0_xxxyyy_0[j] + wp_y[j] * tg_0_xxxyyy_1[j] + 1.5 * fl1_fxn * tg_0_xxxyy_1[j];

                    tg_y_xxxyyz_0[j] = pb_y * tg_0_xxxyyz_0[j] + wp_y[j] * tg_0_xxxyyz_1[j] + fl1_fxn * tg_0_xxxyz_1[j];

                    tg_y_xxxyzz_0[j] = pb_y * tg_0_xxxyzz_0[j] + wp_y[j] * tg_0_xxxyzz_1[j] + 0.5 * fl1_fxn * tg_0_xxxzz_1[j];

                    tg_y_xxxzzz_0[j] = pb_y * tg_0_xxxzzz_0[j] + wp_y[j] * tg_0_xxxzzz_1[j];

                    tg_y_xxyyyy_0[j] = pb_y * tg_0_xxyyyy_0[j] + wp_y[j] * tg_0_xxyyyy_1[j] + 2.0 * fl1_fxn * tg_0_xxyyy_1[j];

                    tg_y_xxyyyz_0[j] = pb_y * tg_0_xxyyyz_0[j] + wp_y[j] * tg_0_xxyyyz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyz_1[j];

                    tg_y_xxyyzz_0[j] = pb_y * tg_0_xxyyzz_0[j] + wp_y[j] * tg_0_xxyyzz_1[j] + fl1_fxn * tg_0_xxyzz_1[j];

                    tg_y_xxyzzz_0[j] = pb_y * tg_0_xxyzzz_0[j] + wp_y[j] * tg_0_xxyzzz_1[j] + 0.5 * fl1_fxn * tg_0_xxzzz_1[j];

                    tg_y_xxzzzz_0[j] = pb_y * tg_0_xxzzzz_0[j] + wp_y[j] * tg_0_xxzzzz_1[j];

                    tg_y_xyyyyy_0[j] = pb_y * tg_0_xyyyyy_0[j] + wp_y[j] * tg_0_xyyyyy_1[j] + 2.5 * fl1_fxn * tg_0_xyyyy_1[j];

                    tg_y_xyyyyz_0[j] = pb_y * tg_0_xyyyyz_0[j] + wp_y[j] * tg_0_xyyyyz_1[j] + 2.0 * fl1_fxn * tg_0_xyyyz_1[j];

                    tg_y_xyyyzz_0[j] = pb_y * tg_0_xyyyzz_0[j] + wp_y[j] * tg_0_xyyyzz_1[j] + 1.5 * fl1_fxn * tg_0_xyyzz_1[j];

                    tg_y_xyyzzz_0[j] = pb_y * tg_0_xyyzzz_0[j] + wp_y[j] * tg_0_xyyzzz_1[j] + fl1_fxn * tg_0_xyzzz_1[j];

                    tg_y_xyzzzz_0[j] = pb_y * tg_0_xyzzzz_0[j] + wp_y[j] * tg_0_xyzzzz_1[j] + 0.5 * fl1_fxn * tg_0_xzzzz_1[j];

                    tg_y_xzzzzz_0[j] = pb_y * tg_0_xzzzzz_0[j] + wp_y[j] * tg_0_xzzzzz_1[j];

                    tg_y_yyyyyy_0[j] = pb_y * tg_0_yyyyyy_0[j] + wp_y[j] * tg_0_yyyyyy_1[j] + 3.0 * fl1_fxn * tg_0_yyyyy_1[j];

                    tg_y_yyyyyz_0[j] = pb_y * tg_0_yyyyyz_0[j] + wp_y[j] * tg_0_yyyyyz_1[j] + 2.5 * fl1_fxn * tg_0_yyyyz_1[j];

                    tg_y_yyyyzz_0[j] = pb_y * tg_0_yyyyzz_0[j] + wp_y[j] * tg_0_yyyyzz_1[j] + 2.0 * fl1_fxn * tg_0_yyyzz_1[j];

                    tg_y_yyyzzz_0[j] = pb_y * tg_0_yyyzzz_0[j] + wp_y[j] * tg_0_yyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_yyzzz_1[j];

                    tg_y_yyzzzz_0[j] = pb_y * tg_0_yyzzzz_0[j] + wp_y[j] * tg_0_yyzzzz_1[j] + fl1_fxn * tg_0_yzzzz_1[j];

                    tg_y_yzzzzz_0[j] = pb_y * tg_0_yzzzzz_0[j] + wp_y[j] * tg_0_yzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzzz_1[j];

                    tg_y_zzzzzz_0[j] = pb_y * tg_0_zzzzzz_0[j] + wp_y[j] * tg_0_zzzzzz_1[j];

                    tg_z_xxxxxx_0[j] = pb_z * tg_0_xxxxxx_0[j] + wp_z[j] * tg_0_xxxxxx_1[j];

                    tg_z_xxxxxy_0[j] = pb_z * tg_0_xxxxxy_0[j] + wp_z[j] * tg_0_xxxxxy_1[j];

                    tg_z_xxxxxz_0[j] = pb_z * tg_0_xxxxxz_0[j] + wp_z[j] * tg_0_xxxxxz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxx_1[j];

                    tg_z_xxxxyy_0[j] = pb_z * tg_0_xxxxyy_0[j] + wp_z[j] * tg_0_xxxxyy_1[j];

                    tg_z_xxxxyz_0[j] = pb_z * tg_0_xxxxyz_0[j] + wp_z[j] * tg_0_xxxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxy_1[j];

                    tg_z_xxxxzz_0[j] = pb_z * tg_0_xxxxzz_0[j] + wp_z[j] * tg_0_xxxxzz_1[j] + fl1_fxn * tg_0_xxxxz_1[j];

                    tg_z_xxxyyy_0[j] = pb_z * tg_0_xxxyyy_0[j] + wp_z[j] * tg_0_xxxyyy_1[j];

                    tg_z_xxxyyz_0[j] = pb_z * tg_0_xxxyyz_0[j] + wp_z[j] * tg_0_xxxyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxyy_1[j];

                    tg_z_xxxyzz_0[j] = pb_z * tg_0_xxxyzz_0[j] + wp_z[j] * tg_0_xxxyzz_1[j] + fl1_fxn * tg_0_xxxyz_1[j];

                    tg_z_xxxzzz_0[j] = pb_z * tg_0_xxxzzz_0[j] + wp_z[j] * tg_0_xxxzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxxzz_1[j];

                    tg_z_xxyyyy_0[j] = pb_z * tg_0_xxyyyy_0[j] + wp_z[j] * tg_0_xxyyyy_1[j];

                    tg_z_xxyyyz_0[j] = pb_z * tg_0_xxyyyz_0[j] + wp_z[j] * tg_0_xxyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxyyy_1[j];

                    tg_z_xxyyzz_0[j] = pb_z * tg_0_xxyyzz_0[j] + wp_z[j] * tg_0_xxyyzz_1[j] + fl1_fxn * tg_0_xxyyz_1[j];

                    tg_z_xxyzzz_0[j] = pb_z * tg_0_xxyzzz_0[j] + wp_z[j] * tg_0_xxyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyzz_1[j];

                    tg_z_xxzzzz_0[j] = pb_z * tg_0_xxzzzz_0[j] + wp_z[j] * tg_0_xxzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxzzz_1[j];

                    tg_z_xyyyyy_0[j] = pb_z * tg_0_xyyyyy_0[j] + wp_z[j] * tg_0_xyyyyy_1[j];

                    tg_z_xyyyyz_0[j] = pb_z * tg_0_xyyyyz_0[j] + wp_z[j] * tg_0_xyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xyyyy_1[j];

                    tg_z_xyyyzz_0[j] = pb_z * tg_0_xyyyzz_0[j] + wp_z[j] * tg_0_xyyyzz_1[j] + fl1_fxn * tg_0_xyyyz_1[j];

                    tg_z_xyyzzz_0[j] = pb_z * tg_0_xyyzzz_0[j] + wp_z[j] * tg_0_xyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xyyzz_1[j];

                    tg_z_xyzzzz_0[j] = pb_z * tg_0_xyzzzz_0[j] + wp_z[j] * tg_0_xyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xyzzz_1[j];

                    tg_z_xzzzzz_0[j] = pb_z * tg_0_xzzzzz_0[j] + wp_z[j] * tg_0_xzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_xzzzz_1[j];

                    tg_z_yyyyyy_0[j] = pb_z * tg_0_yyyyyy_0[j] + wp_z[j] * tg_0_yyyyyy_1[j];

                    tg_z_yyyyyz_0[j] = pb_z * tg_0_yyyyyz_0[j] + wp_z[j] * tg_0_yyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyy_1[j];

                    tg_z_yyyyzz_0[j] = pb_z * tg_0_yyyyzz_0[j] + wp_z[j] * tg_0_yyyyzz_1[j] + fl1_fxn * tg_0_yyyyz_1[j];

                    tg_z_yyyzzz_0[j] = pb_z * tg_0_yyyzzz_0[j] + wp_z[j] * tg_0_yyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_yyyzz_1[j];

                    tg_z_yyzzzz_0[j] = pb_z * tg_0_yyzzzz_0[j] + wp_z[j] * tg_0_yyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_yyzzz_1[j];

                    tg_z_yzzzzz_0[j] = pb_z * tg_0_yzzzzz_0[j] + wp_z[j] * tg_0_yzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_yzzzz_1[j];

                    tg_z_zzzzzz_0[j] = pb_z * tg_0_zzzzzz_0[j] + wp_z[j] * tg_0_zzzzzz_1[j] + 3.0 * fl1_fxn * tg_0_zzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
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
                                             {6, -1, -1, -1},
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_4_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {4, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xxxxx_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx); 

                auto tg_xxxxx_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 1); 

                auto tg_xxxxx_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 2); 

                auto tg_xxxxy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 3); 

                auto tg_xxxxy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 4); 

                auto tg_xxxxy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 5); 

                auto tg_xxxxz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 6); 

                auto tg_xxxxz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 7); 

                auto tg_xxxxz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 8); 

                auto tg_xxxyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 9); 

                auto tg_xxxyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 10); 

                auto tg_xxxyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 11); 

                auto tg_xxxyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 12); 

                auto tg_xxxyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 13); 

                auto tg_xxxyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 14); 

                auto tg_xxxzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 15); 

                auto tg_xxxzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 16); 

                auto tg_xxxzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 17); 

                auto tg_xxyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 18); 

                auto tg_xxyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 19); 

                auto tg_xxyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 20); 

                auto tg_xxyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 21); 

                auto tg_xxyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 22); 

                auto tg_xxyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 23); 

                auto tg_xxyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 24); 

                auto tg_xxyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 25); 

                auto tg_xxyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 26); 

                auto tg_xxzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 27); 

                auto tg_xxzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 28); 

                auto tg_xxzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 29); 

                auto tg_xyyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 30); 

                auto tg_xyyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 31); 

                auto tg_xyyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 32); 

                auto tg_xyyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 33); 

                auto tg_xyyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 34); 

                auto tg_xyyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 35); 

                auto tg_xyyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 36); 

                auto tg_xyyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 37); 

                auto tg_xyyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 38); 

                auto tg_xyzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 39); 

                auto tg_xyzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 40); 

                auto tg_xyzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 41); 

                auto tg_xzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 42); 

                auto tg_xzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 43); 

                auto tg_xzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 44); 

                auto tg_yyyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 45); 

                auto tg_yyyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 46); 

                auto tg_yyyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 47); 

                auto tg_yyyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 48); 

                auto tg_yyyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 49); 

                auto tg_yyyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 50); 

                auto tg_yyyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 51); 

                auto tg_yyyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 52); 

                auto tg_yyyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 53); 

                auto tg_yyzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 54); 

                auto tg_yyzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 55); 

                auto tg_yyzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 56); 

                auto tg_yzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 57); 

                auto tg_yzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 58); 

                auto tg_yzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 59); 

                auto tg_zzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 60); 

                auto tg_zzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 61); 

                auto tg_zzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 62); 

                auto tg_xxxxx_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx); 

                auto tg_xxxxx_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 1); 

                auto tg_xxxxx_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 2); 

                auto tg_xxxxy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 3); 

                auto tg_xxxxy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 4); 

                auto tg_xxxxy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 5); 

                auto tg_xxxxz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 6); 

                auto tg_xxxxz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 7); 

                auto tg_xxxxz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 8); 

                auto tg_xxxyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 9); 

                auto tg_xxxyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 10); 

                auto tg_xxxyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 11); 

                auto tg_xxxyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 12); 

                auto tg_xxxyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 13); 

                auto tg_xxxyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 14); 

                auto tg_xxxzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 15); 

                auto tg_xxxzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 16); 

                auto tg_xxxzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 17); 

                auto tg_xxyyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 18); 

                auto tg_xxyyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 19); 

                auto tg_xxyyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 20); 

                auto tg_xxyyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 21); 

                auto tg_xxyyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 22); 

                auto tg_xxyyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 23); 

                auto tg_xxyzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 24); 

                auto tg_xxyzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 25); 

                auto tg_xxyzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 26); 

                auto tg_xxzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 27); 

                auto tg_xxzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 28); 

                auto tg_xxzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 29); 

                auto tg_xyyyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 30); 

                auto tg_xyyyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 31); 

                auto tg_xyyyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 32); 

                auto tg_xyyyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 33); 

                auto tg_xyyyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 34); 

                auto tg_xyyyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 35); 

                auto tg_xyyzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 36); 

                auto tg_xyyzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 37); 

                auto tg_xyyzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 38); 

                auto tg_xyzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 39); 

                auto tg_xyzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 40); 

                auto tg_xyzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 41); 

                auto tg_xzzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 42); 

                auto tg_xzzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 43); 

                auto tg_xzzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 44); 

                auto tg_yyyyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 45); 

                auto tg_yyyyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 46); 

                auto tg_yyyyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 47); 

                auto tg_yyyyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 48); 

                auto tg_yyyyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 49); 

                auto tg_yyyyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 50); 

                auto tg_yyyzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 51); 

                auto tg_yyyzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 52); 

                auto tg_yyyzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 53); 

                auto tg_yyzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 54); 

                auto tg_yyzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 55); 

                auto tg_yyzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 56); 

                auto tg_yzzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 57); 

                auto tg_yzzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 58); 

                auto tg_yzzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 59); 

                auto tg_zzzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 60); 

                auto tg_zzzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 61); 

                auto tg_zzzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 62); 

                auto tg_xxxx_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx); 

                auto tg_xxxx_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 1); 

                auto tg_xxxx_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 2); 

                auto tg_xxxy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 3); 

                auto tg_xxxy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 4); 

                auto tg_xxxy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 5); 

                auto tg_xxxz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 6); 

                auto tg_xxxz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 7); 

                auto tg_xxxz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 8); 

                auto tg_xxyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 9); 

                auto tg_xxyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 10); 

                auto tg_xxyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 11); 

                auto tg_xxyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 12); 

                auto tg_xxyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 13); 

                auto tg_xxyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 14); 

                auto tg_xxzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 15); 

                auto tg_xxzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 16); 

                auto tg_xxzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 17); 

                auto tg_xyyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 18); 

                auto tg_xyyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 19); 

                auto tg_xyyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 20); 

                auto tg_xyyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 21); 

                auto tg_xyyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 22); 

                auto tg_xyyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 23); 

                auto tg_xyzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 24); 

                auto tg_xyzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 25); 

                auto tg_xyzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 26); 

                auto tg_xzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 27); 

                auto tg_xzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 28); 

                auto tg_xzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 29); 

                auto tg_yyyy_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 30); 

                auto tg_yyyy_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 31); 

                auto tg_yyyy_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 32); 

                auto tg_yyyz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 33); 

                auto tg_yyyz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 34); 

                auto tg_yyyz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 35); 

                auto tg_yyzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 36); 

                auto tg_yyzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 37); 

                auto tg_yyzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 38); 

                auto tg_yzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 39); 

                auto tg_yzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 40); 

                auto tg_yzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 41); 

                auto tg_zzzz_x_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 42); 

                auto tg_zzzz_y_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 43); 

                auto tg_zzzz_z_0 = primBuffer.data(pidx_g_4_1_m0 + 45 * idx + 44); 

                auto tg_xxxx_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx); 

                auto tg_xxxx_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 1); 

                auto tg_xxxx_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 2); 

                auto tg_xxxy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 3); 

                auto tg_xxxy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 4); 

                auto tg_xxxy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 5); 

                auto tg_xxxz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 6); 

                auto tg_xxxz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 7); 

                auto tg_xxxz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 8); 

                auto tg_xxyy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 9); 

                auto tg_xxyy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 10); 

                auto tg_xxyy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 11); 

                auto tg_xxyz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 12); 

                auto tg_xxyz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 13); 

                auto tg_xxyz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 14); 

                auto tg_xxzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 15); 

                auto tg_xxzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 16); 

                auto tg_xxzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 17); 

                auto tg_xyyy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 18); 

                auto tg_xyyy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 19); 

                auto tg_xyyy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 20); 

                auto tg_xyyz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 21); 

                auto tg_xyyz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 22); 

                auto tg_xyyz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 23); 

                auto tg_xyzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 24); 

                auto tg_xyzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 25); 

                auto tg_xyzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 26); 

                auto tg_xzzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 27); 

                auto tg_xzzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 28); 

                auto tg_xzzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 29); 

                auto tg_yyyy_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 30); 

                auto tg_yyyy_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 31); 

                auto tg_yyyy_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 32); 

                auto tg_yyyz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 33); 

                auto tg_yyyz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 34); 

                auto tg_yyyz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 35); 

                auto tg_yyzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 36); 

                auto tg_yyzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 37); 

                auto tg_yyzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 38); 

                auto tg_yzzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 39); 

                auto tg_yzzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 40); 

                auto tg_yzzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 41); 

                auto tg_zzzz_x_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 42); 

                auto tg_zzzz_y_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 43); 

                auto tg_zzzz_z_1 = primBuffer.data(pidx_g_4_1_m1 + 45 * idx + 44); 

                auto tg_xxxxx_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx); 

                auto tg_xxxxy_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 1); 

                auto tg_xxxxz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 2); 

                auto tg_xxxyy_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 3); 

                auto tg_xxxyz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 4); 

                auto tg_xxxzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 5); 

                auto tg_xxyyy_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 6); 

                auto tg_xxyyz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 7); 

                auto tg_xxyzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 8); 

                auto tg_xxzzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 9); 

                auto tg_xyyyy_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 10); 

                auto tg_xyyyz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 11); 

                auto tg_xyyzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 12); 

                auto tg_xyzzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 13); 

                auto tg_xzzzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 14); 

                auto tg_yyyyy_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 15); 

                auto tg_yyyyz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 16); 

                auto tg_yyyzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 17); 

                auto tg_yyzzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 18); 

                auto tg_yzzzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 19); 

                auto tg_zzzzz_0_1 = primBuffer.data(pidx_g_5_0_m1 + 21 * idx + 20); 

                // set up pointers to integrals

                auto tg_xxxxxx_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx); 

                auto tg_xxxxxx_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 1); 

                auto tg_xxxxxx_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 2); 

                auto tg_xxxxxy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 3); 

                auto tg_xxxxxy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 4); 

                auto tg_xxxxxy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 5); 

                auto tg_xxxxxz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 6); 

                auto tg_xxxxxz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 7); 

                auto tg_xxxxxz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 8); 

                auto tg_xxxxyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 9); 

                auto tg_xxxxyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 10); 

                auto tg_xxxxyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 11); 

                auto tg_xxxxyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 12); 

                auto tg_xxxxyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 13); 

                auto tg_xxxxyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 14); 

                auto tg_xxxxzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 15); 

                auto tg_xxxxzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 16); 

                auto tg_xxxxzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 17); 

                auto tg_xxxyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 18); 

                auto tg_xxxyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 19); 

                auto tg_xxxyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 20); 

                auto tg_xxxyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 21); 

                auto tg_xxxyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 22); 

                auto tg_xxxyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 23); 

                auto tg_xxxyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 24); 

                auto tg_xxxyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 25); 

                auto tg_xxxyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 26); 

                auto tg_xxxzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 27); 

                auto tg_xxxzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 28); 

                auto tg_xxxzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 29); 

                auto tg_xxyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 30); 

                auto tg_xxyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 31); 

                auto tg_xxyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 32); 

                auto tg_xxyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 33); 

                auto tg_xxyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 34); 

                auto tg_xxyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 35); 

                auto tg_xxyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 36); 

                auto tg_xxyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 37); 

                auto tg_xxyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 38); 

                auto tg_xxyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 39); 

                auto tg_xxyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 40); 

                auto tg_xxyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 41); 

                auto tg_xxzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 42); 

                auto tg_xxzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 43); 

                auto tg_xxzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 44); 

                auto tg_xyyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 45); 

                auto tg_xyyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 46); 

                auto tg_xyyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 47); 

                auto tg_xyyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 48); 

                auto tg_xyyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 49); 

                auto tg_xyyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 50); 

                auto tg_xyyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 51); 

                auto tg_xyyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 52); 

                auto tg_xyyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 53); 

                auto tg_xyyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 54); 

                auto tg_xyyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 55); 

                auto tg_xyyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 56); 

                auto tg_xyzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 57); 

                auto tg_xyzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 58); 

                auto tg_xyzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 59); 

                auto tg_xzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 60); 

                auto tg_xzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 61); 

                auto tg_xzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 62); 

                auto tg_yyyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 63); 

                auto tg_yyyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 64); 

                auto tg_yyyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 65); 

                auto tg_yyyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 66); 

                auto tg_yyyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 67); 

                auto tg_yyyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 68); 

                auto tg_yyyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 69); 

                auto tg_yyyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 70); 

                auto tg_yyyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 71); 

                auto tg_yyyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 72); 

                auto tg_yyyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 73); 

                auto tg_yyyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 74); 

                auto tg_yyzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 75); 

                auto tg_yyzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 76); 

                auto tg_yyzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 77); 

                auto tg_yzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 78); 

                auto tg_yzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 79); 

                auto tg_yzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 80); 

                auto tg_zzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 81); 

                auto tg_zzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 82); 

                auto tg_zzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 83); 

                #pragma omp simd aligned(fxn, fza, tg_xxxx_x_0, tg_xxxx_x_1, tg_xxxx_y_0, tg_xxxx_y_1, \
                                         tg_xxxx_z_0, tg_xxxx_z_1, tg_xxxxx_0_1, tg_xxxxx_x_0, tg_xxxxx_x_1, tg_xxxxx_y_0, \
                                         tg_xxxxx_y_1, tg_xxxxx_z_0, tg_xxxxx_z_1, tg_xxxxxx_x_0, tg_xxxxxx_y_0, \
                                         tg_xxxxxx_z_0, tg_xxxxxy_x_0, tg_xxxxxy_y_0, tg_xxxxxy_z_0, tg_xxxxxz_x_0, \
                                         tg_xxxxxz_y_0, tg_xxxxxz_z_0, tg_xxxxy_0_1, tg_xxxxy_x_0, tg_xxxxy_x_1, \
                                         tg_xxxxy_y_0, tg_xxxxy_y_1, tg_xxxxy_z_0, tg_xxxxy_z_1, tg_xxxxyy_x_0, \
                                         tg_xxxxyy_y_0, tg_xxxxyy_z_0, tg_xxxxyz_x_0, tg_xxxxyz_y_0, tg_xxxxyz_z_0, \
                                         tg_xxxxz_0_1, tg_xxxxz_x_0, tg_xxxxz_x_1, tg_xxxxz_y_0, tg_xxxxz_y_1, tg_xxxxz_z_0, \
                                         tg_xxxxz_z_1, tg_xxxxzz_x_0, tg_xxxxzz_y_0, tg_xxxxzz_z_0, tg_xxxy_x_0, \
                                         tg_xxxy_x_1, tg_xxxy_y_0, tg_xxxy_y_1, tg_xxxy_z_0, tg_xxxy_z_1, tg_xxxyy_0_1, \
                                         tg_xxxyy_x_0, tg_xxxyy_x_1, tg_xxxyy_y_0, tg_xxxyy_y_1, tg_xxxyy_z_0, tg_xxxyy_z_1, \
                                         tg_xxxyyy_x_0, tg_xxxyyy_y_0, tg_xxxyyy_z_0, tg_xxxyyz_x_0, tg_xxxyyz_y_0, \
                                         tg_xxxyyz_z_0, tg_xxxyz_0_1, tg_xxxyz_x_0, tg_xxxyz_x_1, tg_xxxyz_y_0, tg_xxxyz_y_1, \
                                         tg_xxxyz_z_0, tg_xxxyz_z_1, tg_xxxyzz_x_0, tg_xxxyzz_y_0, tg_xxxyzz_z_0, \
                                         tg_xxxz_x_0, tg_xxxz_x_1, tg_xxxz_y_0, tg_xxxz_y_1, tg_xxxz_z_0, tg_xxxz_z_1, \
                                         tg_xxxzz_0_1, tg_xxxzz_x_0, tg_xxxzz_x_1, tg_xxxzz_y_0, tg_xxxzz_y_1, tg_xxxzz_z_0, \
                                         tg_xxxzz_z_1, tg_xxxzzz_x_0, tg_xxxzzz_y_0, tg_xxxzzz_z_0, tg_xxyy_x_0, \
                                         tg_xxyy_x_1, tg_xxyy_y_0, tg_xxyy_y_1, tg_xxyy_z_0, tg_xxyy_z_1, tg_xxyyy_0_1, \
                                         tg_xxyyy_x_0, tg_xxyyy_x_1, tg_xxyyy_y_0, tg_xxyyy_y_1, tg_xxyyy_z_0, tg_xxyyy_z_1, \
                                         tg_xxyyyy_x_0, tg_xxyyyy_y_0, tg_xxyyyy_z_0, tg_xxyyyz_x_0, tg_xxyyyz_y_0, \
                                         tg_xxyyyz_z_0, tg_xxyyz_0_1, tg_xxyyz_x_0, tg_xxyyz_x_1, tg_xxyyz_y_0, tg_xxyyz_y_1, \
                                         tg_xxyyz_z_0, tg_xxyyz_z_1, tg_xxyyzz_x_0, tg_xxyyzz_y_0, tg_xxyyzz_z_0, \
                                         tg_xxyz_x_0, tg_xxyz_x_1, tg_xxyz_y_0, tg_xxyz_y_1, tg_xxyz_z_0, tg_xxyz_z_1, \
                                         tg_xxyzz_0_1, tg_xxyzz_x_0, tg_xxyzz_x_1, tg_xxyzz_y_0, tg_xxyzz_y_1, tg_xxyzz_z_0, \
                                         tg_xxyzz_z_1, tg_xxyzzz_x_0, tg_xxyzzz_y_0, tg_xxyzzz_z_0, tg_xxzz_x_0, \
                                         tg_xxzz_x_1, tg_xxzz_y_0, tg_xxzz_y_1, tg_xxzz_z_0, tg_xxzz_z_1, tg_xxzzz_0_1, \
                                         tg_xxzzz_x_0, tg_xxzzz_x_1, tg_xxzzz_y_0, tg_xxzzz_y_1, tg_xxzzz_z_0, tg_xxzzz_z_1, \
                                         tg_xxzzzz_x_0, tg_xxzzzz_y_0, tg_xxzzzz_z_0, tg_xyyy_x_0, tg_xyyy_x_1, tg_xyyy_y_0, \
                                         tg_xyyy_y_1, tg_xyyy_z_0, tg_xyyy_z_1, tg_xyyyy_0_1, tg_xyyyy_x_0, tg_xyyyy_x_1, \
                                         tg_xyyyy_y_0, tg_xyyyy_y_1, tg_xyyyy_z_0, tg_xyyyy_z_1, tg_xyyyyy_x_0, \
                                         tg_xyyyyy_y_0, tg_xyyyyy_z_0, tg_xyyyyz_x_0, tg_xyyyyz_y_0, tg_xyyyyz_z_0, \
                                         tg_xyyyz_0_1, tg_xyyyz_x_0, tg_xyyyz_x_1, tg_xyyyz_y_0, tg_xyyyz_y_1, tg_xyyyz_z_0, \
                                         tg_xyyyz_z_1, tg_xyyyzz_x_0, tg_xyyyzz_y_0, tg_xyyyzz_z_0, tg_xyyz_x_0, \
                                         tg_xyyz_x_1, tg_xyyz_y_0, tg_xyyz_y_1, tg_xyyz_z_0, tg_xyyz_z_1, tg_xyyzz_0_1, \
                                         tg_xyyzz_x_0, tg_xyyzz_x_1, tg_xyyzz_y_0, tg_xyyzz_y_1, tg_xyyzz_z_0, tg_xyyzz_z_1, \
                                         tg_xyyzzz_x_0, tg_xyyzzz_y_0, tg_xyyzzz_z_0, tg_xyzz_x_0, tg_xyzz_x_1, tg_xyzz_y_0, \
                                         tg_xyzz_y_1, tg_xyzz_z_0, tg_xyzz_z_1, tg_xyzzz_0_1, tg_xyzzz_x_0, tg_xyzzz_x_1, \
                                         tg_xyzzz_y_0, tg_xyzzz_y_1, tg_xyzzz_z_0, tg_xyzzz_z_1, tg_xyzzzz_x_0, \
                                         tg_xyzzzz_y_0, tg_xyzzzz_z_0, tg_xzzz_x_0, tg_xzzz_x_1, tg_xzzz_y_0, tg_xzzz_y_1, \
                                         tg_xzzz_z_0, tg_xzzz_z_1, tg_xzzzz_0_1, tg_xzzzz_x_0, tg_xzzzz_x_1, tg_xzzzz_y_0, \
                                         tg_xzzzz_y_1, tg_xzzzz_z_0, tg_xzzzz_z_1, tg_xzzzzz_x_0, tg_xzzzzz_y_0, \
                                         tg_xzzzzz_z_0, tg_yyyy_x_0, tg_yyyy_x_1, tg_yyyy_y_0, tg_yyyy_y_1, tg_yyyy_z_0, \
                                         tg_yyyy_z_1, tg_yyyyy_0_1, tg_yyyyy_x_0, tg_yyyyy_x_1, tg_yyyyy_y_0, tg_yyyyy_y_1, \
                                         tg_yyyyy_z_0, tg_yyyyy_z_1, tg_yyyyyy_x_0, tg_yyyyyy_y_0, tg_yyyyyy_z_0, \
                                         tg_yyyyyz_x_0, tg_yyyyyz_y_0, tg_yyyyyz_z_0, tg_yyyyz_0_1, tg_yyyyz_x_0, \
                                         tg_yyyyz_x_1, tg_yyyyz_y_0, tg_yyyyz_y_1, tg_yyyyz_z_0, tg_yyyyz_z_1, \
                                         tg_yyyyzz_x_0, tg_yyyyzz_y_0, tg_yyyyzz_z_0, tg_yyyz_x_0, tg_yyyz_x_1, tg_yyyz_y_0, \
                                         tg_yyyz_y_1, tg_yyyz_z_0, tg_yyyz_z_1, tg_yyyzz_0_1, tg_yyyzz_x_0, tg_yyyzz_x_1, \
                                         tg_yyyzz_y_0, tg_yyyzz_y_1, tg_yyyzz_z_0, tg_yyyzz_z_1, tg_yyyzzz_x_0, \
                                         tg_yyyzzz_y_0, tg_yyyzzz_z_0, tg_yyzz_x_0, tg_yyzz_x_1, tg_yyzz_y_0, tg_yyzz_y_1, \
                                         tg_yyzz_z_0, tg_yyzz_z_1, tg_yyzzz_0_1, tg_yyzzz_x_0, tg_yyzzz_x_1, tg_yyzzz_y_0, \
                                         tg_yyzzz_y_1, tg_yyzzz_z_0, tg_yyzzz_z_1, tg_yyzzzz_x_0, tg_yyzzzz_y_0, \
                                         tg_yyzzzz_z_0, tg_yzzz_x_0, tg_yzzz_x_1, tg_yzzz_y_0, tg_yzzz_y_1, tg_yzzz_z_0, \
                                         tg_yzzz_z_1, tg_yzzzz_0_1, tg_yzzzz_x_0, tg_yzzzz_x_1, tg_yzzzz_y_0, tg_yzzzz_y_1, \
                                         tg_yzzzz_z_0, tg_yzzzz_z_1, tg_yzzzzz_x_0, tg_yzzzzz_y_0, tg_yzzzzz_z_0, \
                                         tg_zzzz_x_0, tg_zzzz_x_1, tg_zzzz_y_0, tg_zzzz_y_1, tg_zzzz_z_0, tg_zzzz_z_1, \
                                         tg_zzzzz_0_1, tg_zzzzz_x_0, tg_zzzzz_x_1, tg_zzzzz_y_0, tg_zzzzz_y_1, tg_zzzzz_z_0, \
                                         tg_zzzzz_z_1, tg_zzzzzz_x_0, tg_zzzzzz_y_0, tg_zzzzzz_z_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxx_x_0[j] = pb_x * tg_xxxxx_x_0[j] + wp_x[j] * tg_xxxxx_x_1[j] + 2.5 * fl1_fx * tg_xxxx_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_x_1[j] + 0.5 * fl1_fxn * tg_xxxxx_0_1[j];

                    tg_xxxxxx_y_0[j] = pb_x * tg_xxxxx_y_0[j] + wp_x[j] * tg_xxxxx_y_1[j] + 2.5 * fl1_fx * tg_xxxx_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_y_1[j];

                    tg_xxxxxx_z_0[j] = pb_x * tg_xxxxx_z_0[j] + wp_x[j] * tg_xxxxx_z_1[j] + 2.5 * fl1_fx * tg_xxxx_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxx_z_1[j];

                    tg_xxxxxy_x_0[j] = pb_x * tg_xxxxy_x_0[j] + wp_x[j] * tg_xxxxy_x_1[j] + 2.0 * fl1_fx * tg_xxxy_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_x_1[j] + 0.5 * fl1_fxn * tg_xxxxy_0_1[j];

                    tg_xxxxxy_y_0[j] = pb_x * tg_xxxxy_y_0[j] + wp_x[j] * tg_xxxxy_y_1[j] + 2.0 * fl1_fx * tg_xxxy_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_y_1[j];

                    tg_xxxxxy_z_0[j] = pb_x * tg_xxxxy_z_0[j] + wp_x[j] * tg_xxxxy_z_1[j] + 2.0 * fl1_fx * tg_xxxy_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxy_z_1[j];

                    tg_xxxxxz_x_0[j] = pb_x * tg_xxxxz_x_0[j] + wp_x[j] * tg_xxxxz_x_1[j] + 2.0 * fl1_fx * tg_xxxz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxz_0_1[j];

                    tg_xxxxxz_y_0[j] = pb_x * tg_xxxxz_y_0[j] + wp_x[j] * tg_xxxxz_y_1[j] + 2.0 * fl1_fx * tg_xxxz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_y_1[j];

                    tg_xxxxxz_z_0[j] = pb_x * tg_xxxxz_z_0[j] + wp_x[j] * tg_xxxxz_z_1[j] + 2.0 * fl1_fx * tg_xxxz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxz_z_1[j];

                    tg_xxxxyy_x_0[j] = pb_x * tg_xxxyy_x_0[j] + wp_x[j] * tg_xxxyy_x_1[j] + 1.5 * fl1_fx * tg_xxyy_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_x_1[j] + 0.5 * fl1_fxn * tg_xxxyy_0_1[j];

                    tg_xxxxyy_y_0[j] = pb_x * tg_xxxyy_y_0[j] + wp_x[j] * tg_xxxyy_y_1[j] + 1.5 * fl1_fx * tg_xxyy_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_y_1[j];

                    tg_xxxxyy_z_0[j] = pb_x * tg_xxxyy_z_0[j] + wp_x[j] * tg_xxxyy_z_1[j] + 1.5 * fl1_fx * tg_xxyy_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyy_z_1[j];

                    tg_xxxxyz_x_0[j] = pb_x * tg_xxxyz_x_0[j] + wp_x[j] * tg_xxxyz_x_1[j] + 1.5 * fl1_fx * tg_xxyz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_x_1[j] + 0.5 * fl1_fxn * tg_xxxyz_0_1[j];

                    tg_xxxxyz_y_0[j] = pb_x * tg_xxxyz_y_0[j] + wp_x[j] * tg_xxxyz_y_1[j] + 1.5 * fl1_fx * tg_xxyz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_y_1[j];

                    tg_xxxxyz_z_0[j] = pb_x * tg_xxxyz_z_0[j] + wp_x[j] * tg_xxxyz_z_1[j] + 1.5 * fl1_fx * tg_xxyz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyz_z_1[j];

                    tg_xxxxzz_x_0[j] = pb_x * tg_xxxzz_x_0[j] + wp_x[j] * tg_xxxzz_x_1[j] + 1.5 * fl1_fx * tg_xxzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxzz_0_1[j];

                    tg_xxxxzz_y_0[j] = pb_x * tg_xxxzz_y_0[j] + wp_x[j] * tg_xxxzz_y_1[j] + 1.5 * fl1_fx * tg_xxzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_y_1[j];

                    tg_xxxxzz_z_0[j] = pb_x * tg_xxxzz_z_0[j] + wp_x[j] * tg_xxxzz_z_1[j] + 1.5 * fl1_fx * tg_xxzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzz_z_1[j];

                    tg_xxxyyy_x_0[j] = pb_x * tg_xxyyy_x_0[j] + wp_x[j] * tg_xxyyy_x_1[j] + fl1_fx * tg_xyyy_x_0[j] - fl1_fx * fl1_fza * tg_xyyy_x_1[j] + 0.5 * fl1_fxn * tg_xxyyy_0_1[j];

                    tg_xxxyyy_y_0[j] = pb_x * tg_xxyyy_y_0[j] + wp_x[j] * tg_xxyyy_y_1[j] + fl1_fx * tg_xyyy_y_0[j] - fl1_fx * fl1_fza * tg_xyyy_y_1[j];

                    tg_xxxyyy_z_0[j] = pb_x * tg_xxyyy_z_0[j] + wp_x[j] * tg_xxyyy_z_1[j] + fl1_fx * tg_xyyy_z_0[j] - fl1_fx * fl1_fza * tg_xyyy_z_1[j];

                    tg_xxxyyz_x_0[j] = pb_x * tg_xxyyz_x_0[j] + wp_x[j] * tg_xxyyz_x_1[j] + fl1_fx * tg_xyyz_x_0[j] - fl1_fx * fl1_fza * tg_xyyz_x_1[j] + 0.5 * fl1_fxn * tg_xxyyz_0_1[j];

                    tg_xxxyyz_y_0[j] = pb_x * tg_xxyyz_y_0[j] + wp_x[j] * tg_xxyyz_y_1[j] + fl1_fx * tg_xyyz_y_0[j] - fl1_fx * fl1_fza * tg_xyyz_y_1[j];

                    tg_xxxyyz_z_0[j] = pb_x * tg_xxyyz_z_0[j] + wp_x[j] * tg_xxyyz_z_1[j] + fl1_fx * tg_xyyz_z_0[j] - fl1_fx * fl1_fza * tg_xyyz_z_1[j];

                    tg_xxxyzz_x_0[j] = pb_x * tg_xxyzz_x_0[j] + wp_x[j] * tg_xxyzz_x_1[j] + fl1_fx * tg_xyzz_x_0[j] - fl1_fx * fl1_fza * tg_xyzz_x_1[j] + 0.5 * fl1_fxn * tg_xxyzz_0_1[j];

                    tg_xxxyzz_y_0[j] = pb_x * tg_xxyzz_y_0[j] + wp_x[j] * tg_xxyzz_y_1[j] + fl1_fx * tg_xyzz_y_0[j] - fl1_fx * fl1_fza * tg_xyzz_y_1[j];

                    tg_xxxyzz_z_0[j] = pb_x * tg_xxyzz_z_0[j] + wp_x[j] * tg_xxyzz_z_1[j] + fl1_fx * tg_xyzz_z_0[j] - fl1_fx * fl1_fza * tg_xyzz_z_1[j];

                    tg_xxxzzz_x_0[j] = pb_x * tg_xxzzz_x_0[j] + wp_x[j] * tg_xxzzz_x_1[j] + fl1_fx * tg_xzzz_x_0[j] - fl1_fx * fl1_fza * tg_xzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxzzz_0_1[j];

                    tg_xxxzzz_y_0[j] = pb_x * tg_xxzzz_y_0[j] + wp_x[j] * tg_xxzzz_y_1[j] + fl1_fx * tg_xzzz_y_0[j] - fl1_fx * fl1_fza * tg_xzzz_y_1[j];

                    tg_xxxzzz_z_0[j] = pb_x * tg_xxzzz_z_0[j] + wp_x[j] * tg_xxzzz_z_1[j] + fl1_fx * tg_xzzz_z_0[j] - fl1_fx * fl1_fza * tg_xzzz_z_1[j];

                    tg_xxyyyy_x_0[j] = pb_x * tg_xyyyy_x_0[j] + wp_x[j] * tg_xyyyy_x_1[j] + 0.5 * fl1_fx * tg_yyyy_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_x_1[j] + 0.5 * fl1_fxn * tg_xyyyy_0_1[j];

                    tg_xxyyyy_y_0[j] = pb_x * tg_xyyyy_y_0[j] + wp_x[j] * tg_xyyyy_y_1[j] + 0.5 * fl1_fx * tg_yyyy_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_y_1[j];

                    tg_xxyyyy_z_0[j] = pb_x * tg_xyyyy_z_0[j] + wp_x[j] * tg_xyyyy_z_1[j] + 0.5 * fl1_fx * tg_yyyy_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyy_z_1[j];

                    tg_xxyyyz_x_0[j] = pb_x * tg_xyyyz_x_0[j] + wp_x[j] * tg_xyyyz_x_1[j] + 0.5 * fl1_fx * tg_yyyz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_x_1[j] + 0.5 * fl1_fxn * tg_xyyyz_0_1[j];

                    tg_xxyyyz_y_0[j] = pb_x * tg_xyyyz_y_0[j] + wp_x[j] * tg_xyyyz_y_1[j] + 0.5 * fl1_fx * tg_yyyz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_y_1[j];

                    tg_xxyyyz_z_0[j] = pb_x * tg_xyyyz_z_0[j] + wp_x[j] * tg_xyyyz_z_1[j] + 0.5 * fl1_fx * tg_yyyz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyz_z_1[j];

                    tg_xxyyzz_x_0[j] = pb_x * tg_xyyzz_x_0[j] + wp_x[j] * tg_xyyzz_x_1[j] + 0.5 * fl1_fx * tg_yyzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_x_1[j] + 0.5 * fl1_fxn * tg_xyyzz_0_1[j];

                    tg_xxyyzz_y_0[j] = pb_x * tg_xyyzz_y_0[j] + wp_x[j] * tg_xyyzz_y_1[j] + 0.5 * fl1_fx * tg_yyzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_y_1[j];

                    tg_xxyyzz_z_0[j] = pb_x * tg_xyyzz_z_0[j] + wp_x[j] * tg_xyyzz_z_1[j] + 0.5 * fl1_fx * tg_yyzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzz_z_1[j];

                    tg_xxyzzz_x_0[j] = pb_x * tg_xyzzz_x_0[j] + wp_x[j] * tg_xyzzz_x_1[j] + 0.5 * fl1_fx * tg_yzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_x_1[j] + 0.5 * fl1_fxn * tg_xyzzz_0_1[j];

                    tg_xxyzzz_y_0[j] = pb_x * tg_xyzzz_y_0[j] + wp_x[j] * tg_xyzzz_y_1[j] + 0.5 * fl1_fx * tg_yzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_y_1[j];

                    tg_xxyzzz_z_0[j] = pb_x * tg_xyzzz_z_0[j] + wp_x[j] * tg_xyzzz_z_1[j] + 0.5 * fl1_fx * tg_yzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzz_z_1[j];

                    tg_xxzzzz_x_0[j] = pb_x * tg_xzzzz_x_0[j] + wp_x[j] * tg_xzzzz_x_1[j] + 0.5 * fl1_fx * tg_zzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_x_1[j] + 0.5 * fl1_fxn * tg_xzzzz_0_1[j];

                    tg_xxzzzz_y_0[j] = pb_x * tg_xzzzz_y_0[j] + wp_x[j] * tg_xzzzz_y_1[j] + 0.5 * fl1_fx * tg_zzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_y_1[j];

                    tg_xxzzzz_z_0[j] = pb_x * tg_xzzzz_z_0[j] + wp_x[j] * tg_xzzzz_z_1[j] + 0.5 * fl1_fx * tg_zzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_z_1[j];

                    tg_xyyyyy_x_0[j] = pb_x * tg_yyyyy_x_0[j] + wp_x[j] * tg_yyyyy_x_1[j] + 0.5 * fl1_fxn * tg_yyyyy_0_1[j];

                    tg_xyyyyy_y_0[j] = pb_x * tg_yyyyy_y_0[j] + wp_x[j] * tg_yyyyy_y_1[j];

                    tg_xyyyyy_z_0[j] = pb_x * tg_yyyyy_z_0[j] + wp_x[j] * tg_yyyyy_z_1[j];

                    tg_xyyyyz_x_0[j] = pb_x * tg_yyyyz_x_0[j] + wp_x[j] * tg_yyyyz_x_1[j] + 0.5 * fl1_fxn * tg_yyyyz_0_1[j];

                    tg_xyyyyz_y_0[j] = pb_x * tg_yyyyz_y_0[j] + wp_x[j] * tg_yyyyz_y_1[j];

                    tg_xyyyyz_z_0[j] = pb_x * tg_yyyyz_z_0[j] + wp_x[j] * tg_yyyyz_z_1[j];

                    tg_xyyyzz_x_0[j] = pb_x * tg_yyyzz_x_0[j] + wp_x[j] * tg_yyyzz_x_1[j] + 0.5 * fl1_fxn * tg_yyyzz_0_1[j];

                    tg_xyyyzz_y_0[j] = pb_x * tg_yyyzz_y_0[j] + wp_x[j] * tg_yyyzz_y_1[j];

                    tg_xyyyzz_z_0[j] = pb_x * tg_yyyzz_z_0[j] + wp_x[j] * tg_yyyzz_z_1[j];

                    tg_xyyzzz_x_0[j] = pb_x * tg_yyzzz_x_0[j] + wp_x[j] * tg_yyzzz_x_1[j] + 0.5 * fl1_fxn * tg_yyzzz_0_1[j];

                    tg_xyyzzz_y_0[j] = pb_x * tg_yyzzz_y_0[j] + wp_x[j] * tg_yyzzz_y_1[j];

                    tg_xyyzzz_z_0[j] = pb_x * tg_yyzzz_z_0[j] + wp_x[j] * tg_yyzzz_z_1[j];

                    tg_xyzzzz_x_0[j] = pb_x * tg_yzzzz_x_0[j] + wp_x[j] * tg_yzzzz_x_1[j] + 0.5 * fl1_fxn * tg_yzzzz_0_1[j];

                    tg_xyzzzz_y_0[j] = pb_x * tg_yzzzz_y_0[j] + wp_x[j] * tg_yzzzz_y_1[j];

                    tg_xyzzzz_z_0[j] = pb_x * tg_yzzzz_z_0[j] + wp_x[j] * tg_yzzzz_z_1[j];

                    tg_xzzzzz_x_0[j] = pb_x * tg_zzzzz_x_0[j] + wp_x[j] * tg_zzzzz_x_1[j] + 0.5 * fl1_fxn * tg_zzzzz_0_1[j];

                    tg_xzzzzz_y_0[j] = pb_x * tg_zzzzz_y_0[j] + wp_x[j] * tg_zzzzz_y_1[j];

                    tg_xzzzzz_z_0[j] = pb_x * tg_zzzzz_z_0[j] + wp_x[j] * tg_zzzzz_z_1[j];

                    tg_yyyyyy_x_0[j] = pb_y * tg_yyyyy_x_0[j] + wp_y[j] * tg_yyyyy_x_1[j] + 2.5 * fl1_fx * tg_yyyy_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_x_1[j];

                    tg_yyyyyy_y_0[j] = pb_y * tg_yyyyy_y_0[j] + wp_y[j] * tg_yyyyy_y_1[j] + 2.5 * fl1_fx * tg_yyyy_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_y_1[j] + 0.5 * fl1_fxn * tg_yyyyy_0_1[j];

                    tg_yyyyyy_z_0[j] = pb_y * tg_yyyyy_z_0[j] + wp_y[j] * tg_yyyyy_z_1[j] + 2.5 * fl1_fx * tg_yyyy_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyy_z_1[j];

                    tg_yyyyyz_x_0[j] = pb_y * tg_yyyyz_x_0[j] + wp_y[j] * tg_yyyyz_x_1[j] + 2.0 * fl1_fx * tg_yyyz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_x_1[j];

                    tg_yyyyyz_y_0[j] = pb_y * tg_yyyyz_y_0[j] + wp_y[j] * tg_yyyyz_y_1[j] + 2.0 * fl1_fx * tg_yyyz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_y_1[j] + 0.5 * fl1_fxn * tg_yyyyz_0_1[j];

                    tg_yyyyyz_z_0[j] = pb_y * tg_yyyyz_z_0[j] + wp_y[j] * tg_yyyyz_z_1[j] + 2.0 * fl1_fx * tg_yyyz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyz_z_1[j];

                    tg_yyyyzz_x_0[j] = pb_y * tg_yyyzz_x_0[j] + wp_y[j] * tg_yyyzz_x_1[j] + 1.5 * fl1_fx * tg_yyzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_x_1[j];

                    tg_yyyyzz_y_0[j] = pb_y * tg_yyyzz_y_0[j] + wp_y[j] * tg_yyyzz_y_1[j] + 1.5 * fl1_fx * tg_yyzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_y_1[j] + 0.5 * fl1_fxn * tg_yyyzz_0_1[j];

                    tg_yyyyzz_z_0[j] = pb_y * tg_yyyzz_z_0[j] + wp_y[j] * tg_yyyzz_z_1[j] + 1.5 * fl1_fx * tg_yyzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzz_z_1[j];

                    tg_yyyzzz_x_0[j] = pb_y * tg_yyzzz_x_0[j] + wp_y[j] * tg_yyzzz_x_1[j] + fl1_fx * tg_yzzz_x_0[j] - fl1_fx * fl1_fza * tg_yzzz_x_1[j];

                    tg_yyyzzz_y_0[j] = pb_y * tg_yyzzz_y_0[j] + wp_y[j] * tg_yyzzz_y_1[j] + fl1_fx * tg_yzzz_y_0[j] - fl1_fx * fl1_fza * tg_yzzz_y_1[j] + 0.5 * fl1_fxn * tg_yyzzz_0_1[j];

                    tg_yyyzzz_z_0[j] = pb_y * tg_yyzzz_z_0[j] + wp_y[j] * tg_yyzzz_z_1[j] + fl1_fx * tg_yzzz_z_0[j] - fl1_fx * fl1_fza * tg_yzzz_z_1[j];

                    tg_yyzzzz_x_0[j] = pb_y * tg_yzzzz_x_0[j] + wp_y[j] * tg_yzzzz_x_1[j] + 0.5 * fl1_fx * tg_zzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_x_1[j];

                    tg_yyzzzz_y_0[j] = pb_y * tg_yzzzz_y_0[j] + wp_y[j] * tg_yzzzz_y_1[j] + 0.5 * fl1_fx * tg_zzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_y_1[j] + 0.5 * fl1_fxn * tg_yzzzz_0_1[j];

                    tg_yyzzzz_z_0[j] = pb_y * tg_yzzzz_z_0[j] + wp_y[j] * tg_yzzzz_z_1[j] + 0.5 * fl1_fx * tg_zzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzz_z_1[j];

                    tg_yzzzzz_x_0[j] = pb_y * tg_zzzzz_x_0[j] + wp_y[j] * tg_zzzzz_x_1[j];

                    tg_yzzzzz_y_0[j] = pb_y * tg_zzzzz_y_0[j] + wp_y[j] * tg_zzzzz_y_1[j] + 0.5 * fl1_fxn * tg_zzzzz_0_1[j];

                    tg_yzzzzz_z_0[j] = pb_y * tg_zzzzz_z_0[j] + wp_y[j] * tg_zzzzz_z_1[j];

                    tg_zzzzzz_x_0[j] = pb_z * tg_zzzzz_x_0[j] + wp_z[j] * tg_zzzzz_x_1[j] + 2.5 * fl1_fx * tg_zzzz_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_x_1[j];

                    tg_zzzzzz_y_0[j] = pb_z * tg_zzzzz_y_0[j] + wp_z[j] * tg_zzzzz_y_1[j] + 2.5 * fl1_fx * tg_zzzz_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_y_1[j];

                    tg_zzzzzz_z_0[j] = pb_z * tg_zzzzz_z_0[j] + wp_z[j] * tg_zzzzz_z_1[j] + 2.5 * fl1_fx * tg_zzzz_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_zzzz_z_1[j] + 0.5 * fl1_fxn * tg_zzzzz_0_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSK(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSPSK_0_54(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSPSK_54_108(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForSPSK_0_54(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,54)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {6, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fxn = osFactors.data(4 * idx);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_0_xxxxxxx_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx); 

                auto tg_0_xxxxxxy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 1); 

                auto tg_0_xxxxxxz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 2); 

                auto tg_0_xxxxxyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 3); 

                auto tg_0_xxxxxyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 4); 

                auto tg_0_xxxxxzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 5); 

                auto tg_0_xxxxyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 6); 

                auto tg_0_xxxxyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 7); 

                auto tg_0_xxxxyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 8); 

                auto tg_0_xxxxzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 9); 

                auto tg_0_xxxyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 10); 

                auto tg_0_xxxyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 11); 

                auto tg_0_xxxyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 12); 

                auto tg_0_xxxyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 13); 

                auto tg_0_xxxzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 14); 

                auto tg_0_xxyyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 15); 

                auto tg_0_xxyyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 16); 

                auto tg_0_xxyyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 17); 

                auto tg_0_xxyyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 18); 

                auto tg_0_xxyzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 19); 

                auto tg_0_xxzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 20); 

                auto tg_0_xyyyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 21); 

                auto tg_0_xyyyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 22); 

                auto tg_0_xyyyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 23); 

                auto tg_0_xyyyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 24); 

                auto tg_0_xyyzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 25); 

                auto tg_0_xyzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 26); 

                auto tg_0_xzzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 27); 

                auto tg_0_yyyyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 28); 

                auto tg_0_yyyyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 29); 

                auto tg_0_yyyyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 30); 

                auto tg_0_yyyyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 31); 

                auto tg_0_yyyzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 32); 

                auto tg_0_yyzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 33); 

                auto tg_0_yzzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 34); 

                auto tg_0_zzzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 35); 

                auto tg_0_xxxxxxx_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx); 

                auto tg_0_xxxxxxy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 1); 

                auto tg_0_xxxxxxz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 2); 

                auto tg_0_xxxxxyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 3); 

                auto tg_0_xxxxxyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 4); 

                auto tg_0_xxxxxzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 5); 

                auto tg_0_xxxxyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 6); 

                auto tg_0_xxxxyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 7); 

                auto tg_0_xxxxyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 8); 

                auto tg_0_xxxxzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 9); 

                auto tg_0_xxxyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 10); 

                auto tg_0_xxxyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 11); 

                auto tg_0_xxxyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 12); 

                auto tg_0_xxxyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 13); 

                auto tg_0_xxxzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 14); 

                auto tg_0_xxyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 15); 

                auto tg_0_xxyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 16); 

                auto tg_0_xxyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 17); 

                auto tg_0_xxyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 18); 

                auto tg_0_xxyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 19); 

                auto tg_0_xxzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 20); 

                auto tg_0_xyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 21); 

                auto tg_0_xyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 22); 

                auto tg_0_xyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 23); 

                auto tg_0_xyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 24); 

                auto tg_0_xyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 25); 

                auto tg_0_xyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 26); 

                auto tg_0_xzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 27); 

                auto tg_0_yyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 28); 

                auto tg_0_yyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 29); 

                auto tg_0_yyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 30); 

                auto tg_0_yyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 31); 

                auto tg_0_yyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 32); 

                auto tg_0_yyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 33); 

                auto tg_0_yzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 34); 

                auto tg_0_zzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 35); 

                auto tg_0_xxxxxx_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx); 

                auto tg_0_xxxxxy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 1); 

                auto tg_0_xxxxxz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 2); 

                auto tg_0_xxxxyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 3); 

                auto tg_0_xxxxyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 4); 

                auto tg_0_xxxxzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 5); 

                auto tg_0_xxxyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 6); 

                auto tg_0_xxxyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 7); 

                auto tg_0_xxxyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 8); 

                auto tg_0_xxxzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 9); 

                auto tg_0_xxyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 10); 

                auto tg_0_xxyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 11); 

                auto tg_0_xxyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 12); 

                auto tg_0_xxyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 13); 

                auto tg_0_xxzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 14); 

                auto tg_0_xyyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 15); 

                auto tg_0_xyyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 16); 

                auto tg_0_xyyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 17); 

                auto tg_0_xyyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 18); 

                auto tg_0_xyzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 19); 

                auto tg_0_xzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 20); 

                auto tg_0_yyyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 21); 

                auto tg_0_yyyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 22); 

                auto tg_0_yyyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 23); 

                auto tg_0_yyyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 24); 

                auto tg_0_yyzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 25); 

                auto tg_0_yzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 26); 

                auto tg_0_zzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 27); 

                // set up pointers to integrals

                auto tg_x_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx); 

                auto tg_x_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 1); 

                auto tg_x_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 2); 

                auto tg_x_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 3); 

                auto tg_x_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 4); 

                auto tg_x_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 5); 

                auto tg_x_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 6); 

                auto tg_x_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 7); 

                auto tg_x_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 8); 

                auto tg_x_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 9); 

                auto tg_x_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 10); 

                auto tg_x_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 11); 

                auto tg_x_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 12); 

                auto tg_x_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 13); 

                auto tg_x_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 14); 

                auto tg_x_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 15); 

                auto tg_x_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 16); 

                auto tg_x_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 17); 

                auto tg_x_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 18); 

                auto tg_x_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 19); 

                auto tg_x_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 20); 

                auto tg_x_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 21); 

                auto tg_x_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 22); 

                auto tg_x_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 23); 

                auto tg_x_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 24); 

                auto tg_x_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 25); 

                auto tg_x_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 26); 

                auto tg_x_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 27); 

                auto tg_x_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 28); 

                auto tg_x_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 29); 

                auto tg_x_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 30); 

                auto tg_x_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 31); 

                auto tg_x_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 32); 

                auto tg_x_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 33); 

                auto tg_x_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 34); 

                auto tg_x_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 35); 

                auto tg_y_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 36); 

                auto tg_y_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 37); 

                auto tg_y_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 38); 

                auto tg_y_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 39); 

                auto tg_y_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 40); 

                auto tg_y_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 41); 

                auto tg_y_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 42); 

                auto tg_y_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 43); 

                auto tg_y_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 44); 

                auto tg_y_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 45); 

                auto tg_y_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 46); 

                auto tg_y_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 47); 

                auto tg_y_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 48); 

                auto tg_y_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 49); 

                auto tg_y_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 50); 

                auto tg_y_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 51); 

                auto tg_y_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 52); 

                auto tg_y_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 53); 

                // Batch of Integrals (0,54)

                #pragma omp simd aligned(fxn, tg_0_xxxxxx_1, tg_0_xxxxxxx_0, tg_0_xxxxxxx_1, tg_0_xxxxxxy_0, \
                                         tg_0_xxxxxxy_1, tg_0_xxxxxxz_0, tg_0_xxxxxxz_1, tg_0_xxxxxy_1, tg_0_xxxxxyy_0, \
                                         tg_0_xxxxxyy_1, tg_0_xxxxxyz_0, tg_0_xxxxxyz_1, tg_0_xxxxxz_1, tg_0_xxxxxzz_0, \
                                         tg_0_xxxxxzz_1, tg_0_xxxxyy_1, tg_0_xxxxyyy_0, tg_0_xxxxyyy_1, tg_0_xxxxyyz_0, \
                                         tg_0_xxxxyyz_1, tg_0_xxxxyz_1, tg_0_xxxxyzz_0, tg_0_xxxxyzz_1, tg_0_xxxxzz_1, \
                                         tg_0_xxxxzzz_0, tg_0_xxxxzzz_1, tg_0_xxxyyy_1, tg_0_xxxyyyy_0, tg_0_xxxyyyy_1, \
                                         tg_0_xxxyyyz_0, tg_0_xxxyyyz_1, tg_0_xxxyyz_1, tg_0_xxxyyzz_0, tg_0_xxxyyzz_1, \
                                         tg_0_xxxyzz_1, tg_0_xxxyzzz_0, tg_0_xxxyzzz_1, tg_0_xxxzzz_1, tg_0_xxxzzzz_0, \
                                         tg_0_xxxzzzz_1, tg_0_xxyyyy_1, tg_0_xxyyyyy_0, tg_0_xxyyyyy_1, tg_0_xxyyyyz_0, \
                                         tg_0_xxyyyyz_1, tg_0_xxyyyz_1, tg_0_xxyyyzz_0, tg_0_xxyyyzz_1, tg_0_xxyyzz_1, \
                                         tg_0_xxyyzzz_0, tg_0_xxyyzzz_1, tg_0_xxyzzz_1, tg_0_xxyzzzz_0, tg_0_xxyzzzz_1, \
                                         tg_0_xxzzzz_1, tg_0_xxzzzzz_0, tg_0_xxzzzzz_1, tg_0_xyyyyy_1, tg_0_xyyyyyy_0, \
                                         tg_0_xyyyyyy_1, tg_0_xyyyyyz_0, tg_0_xyyyyyz_1, tg_0_xyyyyz_1, tg_0_xyyyyzz_0, \
                                         tg_0_xyyyyzz_1, tg_0_xyyyzz_1, tg_0_xyyyzzz_0, tg_0_xyyyzzz_1, tg_0_xyyzzz_1, \
                                         tg_0_xyyzzzz_0, tg_0_xyyzzzz_1, tg_0_xyzzzz_1, tg_0_xyzzzzz_0, tg_0_xyzzzzz_1, \
                                         tg_0_xzzzzz_1, tg_0_xzzzzzz_0, tg_0_xzzzzzz_1, tg_0_yyyyyy_1, tg_0_yyyyyyy_0, \
                                         tg_0_yyyyyyy_1, tg_0_yyyyyyz_0, tg_0_yyyyyyz_1, tg_0_yyyyyz_1, tg_0_yyyyyzz_0, \
                                         tg_0_yyyyyzz_1, tg_0_yyyyzz_1, tg_0_yyyyzzz_0, tg_0_yyyyzzz_1, tg_0_yyyzzz_1, \
                                         tg_0_yyyzzzz_0, tg_0_yyyzzzz_1, tg_0_yyzzzz_1, tg_0_yyzzzzz_0, tg_0_yyzzzzz_1, \
                                         tg_0_yzzzzz_1, tg_0_yzzzzzz_0, tg_0_yzzzzzz_1, tg_0_zzzzzz_1, tg_0_zzzzzzz_0, \
                                         tg_0_zzzzzzz_1, tg_x_xxxxxxx_0, tg_x_xxxxxxy_0, tg_x_xxxxxxz_0, tg_x_xxxxxyy_0, \
                                         tg_x_xxxxxyz_0, tg_x_xxxxxzz_0, tg_x_xxxxyyy_0, tg_x_xxxxyyz_0, tg_x_xxxxyzz_0, \
                                         tg_x_xxxxzzz_0, tg_x_xxxyyyy_0, tg_x_xxxyyyz_0, tg_x_xxxyyzz_0, tg_x_xxxyzzz_0, \
                                         tg_x_xxxzzzz_0, tg_x_xxyyyyy_0, tg_x_xxyyyyz_0, tg_x_xxyyyzz_0, tg_x_xxyyzzz_0, \
                                         tg_x_xxyzzzz_0, tg_x_xxzzzzz_0, tg_x_xyyyyyy_0, tg_x_xyyyyyz_0, tg_x_xyyyyzz_0, \
                                         tg_x_xyyyzzz_0, tg_x_xyyzzzz_0, tg_x_xyzzzzz_0, tg_x_xzzzzzz_0, tg_x_yyyyyyy_0, \
                                         tg_x_yyyyyyz_0, tg_x_yyyyyzz_0, tg_x_yyyyzzz_0, tg_x_yyyzzzz_0, tg_x_yyzzzzz_0, \
                                         tg_x_yzzzzzz_0, tg_x_zzzzzzz_0, tg_y_xxxxxxx_0, tg_y_xxxxxxy_0, tg_y_xxxxxxz_0, \
                                         tg_y_xxxxxyy_0, tg_y_xxxxxyz_0, tg_y_xxxxxzz_0, tg_y_xxxxyyy_0, tg_y_xxxxyyz_0, \
                                         tg_y_xxxxyzz_0, tg_y_xxxxzzz_0, tg_y_xxxyyyy_0, tg_y_xxxyyyz_0, tg_y_xxxyyzz_0, \
                                         tg_y_xxxyzzz_0, tg_y_xxxzzzz_0, tg_y_xxyyyyy_0, tg_y_xxyyyyz_0, tg_y_xxyyyzz_0, wp_x, \
                                         wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_xxxxxxx_0[j] = pb_x * tg_0_xxxxxxx_0[j] + wp_x[j] * tg_0_xxxxxxx_1[j] + 3.5 * fl1_fxn * tg_0_xxxxxx_1[j];

                    tg_x_xxxxxxy_0[j] = pb_x * tg_0_xxxxxxy_0[j] + wp_x[j] * tg_0_xxxxxxy_1[j] + 3.0 * fl1_fxn * tg_0_xxxxxy_1[j];

                    tg_x_xxxxxxz_0[j] = pb_x * tg_0_xxxxxxz_0[j] + wp_x[j] * tg_0_xxxxxxz_1[j] + 3.0 * fl1_fxn * tg_0_xxxxxz_1[j];

                    tg_x_xxxxxyy_0[j] = pb_x * tg_0_xxxxxyy_0[j] + wp_x[j] * tg_0_xxxxxyy_1[j] + 2.5 * fl1_fxn * tg_0_xxxxyy_1[j];

                    tg_x_xxxxxyz_0[j] = pb_x * tg_0_xxxxxyz_0[j] + wp_x[j] * tg_0_xxxxxyz_1[j] + 2.5 * fl1_fxn * tg_0_xxxxyz_1[j];

                    tg_x_xxxxxzz_0[j] = pb_x * tg_0_xxxxxzz_0[j] + wp_x[j] * tg_0_xxxxxzz_1[j] + 2.5 * fl1_fxn * tg_0_xxxxzz_1[j];

                    tg_x_xxxxyyy_0[j] = pb_x * tg_0_xxxxyyy_0[j] + wp_x[j] * tg_0_xxxxyyy_1[j] + 2.0 * fl1_fxn * tg_0_xxxyyy_1[j];

                    tg_x_xxxxyyz_0[j] = pb_x * tg_0_xxxxyyz_0[j] + wp_x[j] * tg_0_xxxxyyz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyyz_1[j];

                    tg_x_xxxxyzz_0[j] = pb_x * tg_0_xxxxyzz_0[j] + wp_x[j] * tg_0_xxxxyzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyzz_1[j];

                    tg_x_xxxxzzz_0[j] = pb_x * tg_0_xxxxzzz_0[j] + wp_x[j] * tg_0_xxxxzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxzzz_1[j];

                    tg_x_xxxyyyy_0[j] = pb_x * tg_0_xxxyyyy_0[j] + wp_x[j] * tg_0_xxxyyyy_1[j] + 1.5 * fl1_fxn * tg_0_xxyyyy_1[j];

                    tg_x_xxxyyyz_0[j] = pb_x * tg_0_xxxyyyz_0[j] + wp_x[j] * tg_0_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyyz_1[j];

                    tg_x_xxxyyzz_0[j] = pb_x * tg_0_xxxyyzz_0[j] + wp_x[j] * tg_0_xxxyyzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyzz_1[j];

                    tg_x_xxxyzzz_0[j] = pb_x * tg_0_xxxyzzz_0[j] + wp_x[j] * tg_0_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyzzz_1[j];

                    tg_x_xxxzzzz_0[j] = pb_x * tg_0_xxxzzzz_0[j] + wp_x[j] * tg_0_xxxzzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxzzzz_1[j];

                    tg_x_xxyyyyy_0[j] = pb_x * tg_0_xxyyyyy_0[j] + wp_x[j] * tg_0_xxyyyyy_1[j] + fl1_fxn * tg_0_xyyyyy_1[j];

                    tg_x_xxyyyyz_0[j] = pb_x * tg_0_xxyyyyz_0[j] + wp_x[j] * tg_0_xxyyyyz_1[j] + fl1_fxn * tg_0_xyyyyz_1[j];

                    tg_x_xxyyyzz_0[j] = pb_x * tg_0_xxyyyzz_0[j] + wp_x[j] * tg_0_xxyyyzz_1[j] + fl1_fxn * tg_0_xyyyzz_1[j];

                    tg_x_xxyyzzz_0[j] = pb_x * tg_0_xxyyzzz_0[j] + wp_x[j] * tg_0_xxyyzzz_1[j] + fl1_fxn * tg_0_xyyzzz_1[j];

                    tg_x_xxyzzzz_0[j] = pb_x * tg_0_xxyzzzz_0[j] + wp_x[j] * tg_0_xxyzzzz_1[j] + fl1_fxn * tg_0_xyzzzz_1[j];

                    tg_x_xxzzzzz_0[j] = pb_x * tg_0_xxzzzzz_0[j] + wp_x[j] * tg_0_xxzzzzz_1[j] + fl1_fxn * tg_0_xzzzzz_1[j];

                    tg_x_xyyyyyy_0[j] = pb_x * tg_0_xyyyyyy_0[j] + wp_x[j] * tg_0_xyyyyyy_1[j] + 0.5 * fl1_fxn * tg_0_yyyyyy_1[j];

                    tg_x_xyyyyyz_0[j] = pb_x * tg_0_xyyyyyz_0[j] + wp_x[j] * tg_0_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyyz_1[j];

                    tg_x_xyyyyzz_0[j] = pb_x * tg_0_xyyyyzz_0[j] + wp_x[j] * tg_0_xyyyyzz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyzz_1[j];

                    tg_x_xyyyzzz_0[j] = pb_x * tg_0_xyyyzzz_0[j] + wp_x[j] * tg_0_xyyyzzz_1[j] + 0.5 * fl1_fxn * tg_0_yyyzzz_1[j];

                    tg_x_xyyzzzz_0[j] = pb_x * tg_0_xyyzzzz_0[j] + wp_x[j] * tg_0_xyyzzzz_1[j] + 0.5 * fl1_fxn * tg_0_yyzzzz_1[j];

                    tg_x_xyzzzzz_0[j] = pb_x * tg_0_xyzzzzz_0[j] + wp_x[j] * tg_0_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_yzzzzz_1[j];

                    tg_x_xzzzzzz_0[j] = pb_x * tg_0_xzzzzzz_0[j] + wp_x[j] * tg_0_xzzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzzzz_1[j];

                    tg_x_yyyyyyy_0[j] = pb_x * tg_0_yyyyyyy_0[j] + wp_x[j] * tg_0_yyyyyyy_1[j];

                    tg_x_yyyyyyz_0[j] = pb_x * tg_0_yyyyyyz_0[j] + wp_x[j] * tg_0_yyyyyyz_1[j];

                    tg_x_yyyyyzz_0[j] = pb_x * tg_0_yyyyyzz_0[j] + wp_x[j] * tg_0_yyyyyzz_1[j];

                    tg_x_yyyyzzz_0[j] = pb_x * tg_0_yyyyzzz_0[j] + wp_x[j] * tg_0_yyyyzzz_1[j];

                    tg_x_yyyzzzz_0[j] = pb_x * tg_0_yyyzzzz_0[j] + wp_x[j] * tg_0_yyyzzzz_1[j];

                    tg_x_yyzzzzz_0[j] = pb_x * tg_0_yyzzzzz_0[j] + wp_x[j] * tg_0_yyzzzzz_1[j];

                    tg_x_yzzzzzz_0[j] = pb_x * tg_0_yzzzzzz_0[j] + wp_x[j] * tg_0_yzzzzzz_1[j];

                    tg_x_zzzzzzz_0[j] = pb_x * tg_0_zzzzzzz_0[j] + wp_x[j] * tg_0_zzzzzzz_1[j];

                    tg_y_xxxxxxx_0[j] = pb_y * tg_0_xxxxxxx_0[j] + wp_y[j] * tg_0_xxxxxxx_1[j];

                    tg_y_xxxxxxy_0[j] = pb_y * tg_0_xxxxxxy_0[j] + wp_y[j] * tg_0_xxxxxxy_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxx_1[j];

                    tg_y_xxxxxxz_0[j] = pb_y * tg_0_xxxxxxz_0[j] + wp_y[j] * tg_0_xxxxxxz_1[j];

                    tg_y_xxxxxyy_0[j] = pb_y * tg_0_xxxxxyy_0[j] + wp_y[j] * tg_0_xxxxxyy_1[j] + fl1_fxn * tg_0_xxxxxy_1[j];

                    tg_y_xxxxxyz_0[j] = pb_y * tg_0_xxxxxyz_0[j] + wp_y[j] * tg_0_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxz_1[j];

                    tg_y_xxxxxzz_0[j] = pb_y * tg_0_xxxxxzz_0[j] + wp_y[j] * tg_0_xxxxxzz_1[j];

                    tg_y_xxxxyyy_0[j] = pb_y * tg_0_xxxxyyy_0[j] + wp_y[j] * tg_0_xxxxyyy_1[j] + 1.5 * fl1_fxn * tg_0_xxxxyy_1[j];

                    tg_y_xxxxyyz_0[j] = pb_y * tg_0_xxxxyyz_0[j] + wp_y[j] * tg_0_xxxxyyz_1[j] + fl1_fxn * tg_0_xxxxyz_1[j];

                    tg_y_xxxxyzz_0[j] = pb_y * tg_0_xxxxyzz_0[j] + wp_y[j] * tg_0_xxxxyzz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxzz_1[j];

                    tg_y_xxxxzzz_0[j] = pb_y * tg_0_xxxxzzz_0[j] + wp_y[j] * tg_0_xxxxzzz_1[j];

                    tg_y_xxxyyyy_0[j] = pb_y * tg_0_xxxyyyy_0[j] + wp_y[j] * tg_0_xxxyyyy_1[j] + 2.0 * fl1_fxn * tg_0_xxxyyy_1[j];

                    tg_y_xxxyyyz_0[j] = pb_y * tg_0_xxxyyyz_0[j] + wp_y[j] * tg_0_xxxyyyz_1[j] + 1.5 * fl1_fxn * tg_0_xxxyyz_1[j];

                    tg_y_xxxyyzz_0[j] = pb_y * tg_0_xxxyyzz_0[j] + wp_y[j] * tg_0_xxxyyzz_1[j] + fl1_fxn * tg_0_xxxyzz_1[j];

                    tg_y_xxxyzzz_0[j] = pb_y * tg_0_xxxyzzz_0[j] + wp_y[j] * tg_0_xxxyzzz_1[j] + 0.5 * fl1_fxn * tg_0_xxxzzz_1[j];

                    tg_y_xxxzzzz_0[j] = pb_y * tg_0_xxxzzzz_0[j] + wp_y[j] * tg_0_xxxzzzz_1[j];

                    tg_y_xxyyyyy_0[j] = pb_y * tg_0_xxyyyyy_0[j] + wp_y[j] * tg_0_xxyyyyy_1[j] + 2.5 * fl1_fxn * tg_0_xxyyyy_1[j];

                    tg_y_xxyyyyz_0[j] = pb_y * tg_0_xxyyyyz_0[j] + wp_y[j] * tg_0_xxyyyyz_1[j] + 2.0 * fl1_fxn * tg_0_xxyyyz_1[j];

                    tg_y_xxyyyzz_0[j] = pb_y * tg_0_xxyyyzz_0[j] + wp_y[j] * tg_0_xxyyyzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSK_54_108(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (54,108)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_7_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_6_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {6, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fxn = osFactors.data(4 * idx);

                // set up distances R(PB) = P - B

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xxxxxxx_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx); 

                auto tg_0_xxxxxxy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 1); 

                auto tg_0_xxxxxxz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 2); 

                auto tg_0_xxxxxyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 3); 

                auto tg_0_xxxxxyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 4); 

                auto tg_0_xxxxxzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 5); 

                auto tg_0_xxxxyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 6); 

                auto tg_0_xxxxyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 7); 

                auto tg_0_xxxxyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 8); 

                auto tg_0_xxxxzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 9); 

                auto tg_0_xxxyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 10); 

                auto tg_0_xxxyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 11); 

                auto tg_0_xxxyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 12); 

                auto tg_0_xxxyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 13); 

                auto tg_0_xxxzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 14); 

                auto tg_0_xxyyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 15); 

                auto tg_0_xxyyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 16); 

                auto tg_0_xxyyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 17); 

                auto tg_0_xxyyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 18); 

                auto tg_0_xxyzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 19); 

                auto tg_0_xxzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 20); 

                auto tg_0_xyyyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 21); 

                auto tg_0_xyyyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 22); 

                auto tg_0_xyyyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 23); 

                auto tg_0_xyyyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 24); 

                auto tg_0_xyyzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 25); 

                auto tg_0_xyzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 26); 

                auto tg_0_xzzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 27); 

                auto tg_0_yyyyyyy_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 28); 

                auto tg_0_yyyyyyz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 29); 

                auto tg_0_yyyyyzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 30); 

                auto tg_0_yyyyzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 31); 

                auto tg_0_yyyzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 32); 

                auto tg_0_yyzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 33); 

                auto tg_0_yzzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 34); 

                auto tg_0_zzzzzzz_0 = primBuffer.data(pidx_g_0_7_m0 + 36 * idx + 35); 

                auto tg_0_xxxxxxx_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx); 

                auto tg_0_xxxxxxy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 1); 

                auto tg_0_xxxxxxz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 2); 

                auto tg_0_xxxxxyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 3); 

                auto tg_0_xxxxxyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 4); 

                auto tg_0_xxxxxzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 5); 

                auto tg_0_xxxxyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 6); 

                auto tg_0_xxxxyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 7); 

                auto tg_0_xxxxyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 8); 

                auto tg_0_xxxxzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 9); 

                auto tg_0_xxxyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 10); 

                auto tg_0_xxxyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 11); 

                auto tg_0_xxxyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 12); 

                auto tg_0_xxxyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 13); 

                auto tg_0_xxxzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 14); 

                auto tg_0_xxyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 15); 

                auto tg_0_xxyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 16); 

                auto tg_0_xxyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 17); 

                auto tg_0_xxyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 18); 

                auto tg_0_xxyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 19); 

                auto tg_0_xxzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 20); 

                auto tg_0_xyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 21); 

                auto tg_0_xyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 22); 

                auto tg_0_xyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 23); 

                auto tg_0_xyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 24); 

                auto tg_0_xyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 25); 

                auto tg_0_xyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 26); 

                auto tg_0_xzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 27); 

                auto tg_0_yyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 28); 

                auto tg_0_yyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 29); 

                auto tg_0_yyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 30); 

                auto tg_0_yyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 31); 

                auto tg_0_yyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 32); 

                auto tg_0_yyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 33); 

                auto tg_0_yzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 34); 

                auto tg_0_zzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 35); 

                auto tg_0_xxxxxx_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx); 

                auto tg_0_xxxxxy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 1); 

                auto tg_0_xxxxxz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 2); 

                auto tg_0_xxxxyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 3); 

                auto tg_0_xxxxyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 4); 

                auto tg_0_xxxxzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 5); 

                auto tg_0_xxxyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 6); 

                auto tg_0_xxxyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 7); 

                auto tg_0_xxxyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 8); 

                auto tg_0_xxxzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 9); 

                auto tg_0_xxyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 10); 

                auto tg_0_xxyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 11); 

                auto tg_0_xxyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 12); 

                auto tg_0_xxyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 13); 

                auto tg_0_xxzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 14); 

                auto tg_0_xyyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 15); 

                auto tg_0_xyyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 16); 

                auto tg_0_xyyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 17); 

                auto tg_0_xyyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 18); 

                auto tg_0_xyzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 19); 

                auto tg_0_xzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 20); 

                auto tg_0_yyyyyy_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 21); 

                auto tg_0_yyyyyz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 22); 

                auto tg_0_yyyyzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 23); 

                auto tg_0_yyyzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 24); 

                auto tg_0_yyzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 25); 

                auto tg_0_yzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 26); 

                auto tg_0_zzzzzz_1 = primBuffer.data(pidx_g_0_6_m1 + 28 * idx + 27); 

                // set up pointers to integrals

                auto tg_y_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 54); 

                auto tg_y_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 55); 

                auto tg_y_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 56); 

                auto tg_y_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 57); 

                auto tg_y_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 58); 

                auto tg_y_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 59); 

                auto tg_y_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 60); 

                auto tg_y_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 61); 

                auto tg_y_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 62); 

                auto tg_y_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 63); 

                auto tg_y_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 64); 

                auto tg_y_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 65); 

                auto tg_y_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 66); 

                auto tg_y_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 67); 

                auto tg_y_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 68); 

                auto tg_y_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 69); 

                auto tg_y_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 70); 

                auto tg_y_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 71); 

                auto tg_z_xxxxxxx_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 72); 

                auto tg_z_xxxxxxy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 73); 

                auto tg_z_xxxxxxz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 74); 

                auto tg_z_xxxxxyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 75); 

                auto tg_z_xxxxxyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 76); 

                auto tg_z_xxxxxzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 77); 

                auto tg_z_xxxxyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 78); 

                auto tg_z_xxxxyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 79); 

                auto tg_z_xxxxyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 80); 

                auto tg_z_xxxxzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 81); 

                auto tg_z_xxxyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 82); 

                auto tg_z_xxxyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 83); 

                auto tg_z_xxxyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 84); 

                auto tg_z_xxxyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 85); 

                auto tg_z_xxxzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 86); 

                auto tg_z_xxyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 87); 

                auto tg_z_xxyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 88); 

                auto tg_z_xxyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 89); 

                auto tg_z_xxyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 90); 

                auto tg_z_xxyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 91); 

                auto tg_z_xxzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 92); 

                auto tg_z_xyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 93); 

                auto tg_z_xyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 94); 

                auto tg_z_xyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 95); 

                auto tg_z_xyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 96); 

                auto tg_z_xyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 97); 

                auto tg_z_xyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 98); 

                auto tg_z_xzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 99); 

                auto tg_z_yyyyyyy_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 100); 

                auto tg_z_yyyyyyz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 101); 

                auto tg_z_yyyyyzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 102); 

                auto tg_z_yyyyzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 103); 

                auto tg_z_yyyzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 104); 

                auto tg_z_yyzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 105); 

                auto tg_z_yzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 106); 

                auto tg_z_zzzzzzz_0 = primBuffer.data(pidx_g_1_7_m0 + 108 * idx + 107); 

                // Batch of Integrals (54,108)

                #pragma omp simd aligned(fxn, tg_0_xxxxxx_1, tg_0_xxxxxxx_0, tg_0_xxxxxxx_1, tg_0_xxxxxxy_0, \
                                         tg_0_xxxxxxy_1, tg_0_xxxxxxz_0, tg_0_xxxxxxz_1, tg_0_xxxxxy_1, tg_0_xxxxxyy_0, \
                                         tg_0_xxxxxyy_1, tg_0_xxxxxyz_0, tg_0_xxxxxyz_1, tg_0_xxxxxz_1, tg_0_xxxxxzz_0, \
                                         tg_0_xxxxxzz_1, tg_0_xxxxyy_1, tg_0_xxxxyyy_0, tg_0_xxxxyyy_1, tg_0_xxxxyyz_0, \
                                         tg_0_xxxxyyz_1, tg_0_xxxxyz_1, tg_0_xxxxyzz_0, tg_0_xxxxyzz_1, tg_0_xxxxzz_1, \
                                         tg_0_xxxxzzz_0, tg_0_xxxxzzz_1, tg_0_xxxyyy_1, tg_0_xxxyyyy_0, tg_0_xxxyyyy_1, \
                                         tg_0_xxxyyyz_0, tg_0_xxxyyyz_1, tg_0_xxxyyz_1, tg_0_xxxyyzz_0, tg_0_xxxyyzz_1, \
                                         tg_0_xxxyzz_1, tg_0_xxxyzzz_0, tg_0_xxxyzzz_1, tg_0_xxxzzz_1, tg_0_xxxzzzz_0, \
                                         tg_0_xxxzzzz_1, tg_0_xxyyyy_1, tg_0_xxyyyyy_0, tg_0_xxyyyyy_1, tg_0_xxyyyyz_0, \
                                         tg_0_xxyyyyz_1, tg_0_xxyyyz_1, tg_0_xxyyyzz_0, tg_0_xxyyyzz_1, tg_0_xxyyzz_1, \
                                         tg_0_xxyyzzz_0, tg_0_xxyyzzz_1, tg_0_xxyzzz_1, tg_0_xxyzzzz_0, tg_0_xxyzzzz_1, \
                                         tg_0_xxzzzz_1, tg_0_xxzzzzz_0, tg_0_xxzzzzz_1, tg_0_xyyyyy_1, tg_0_xyyyyyy_0, \
                                         tg_0_xyyyyyy_1, tg_0_xyyyyyz_0, tg_0_xyyyyyz_1, tg_0_xyyyyz_1, tg_0_xyyyyzz_0, \
                                         tg_0_xyyyyzz_1, tg_0_xyyyzz_1, tg_0_xyyyzzz_0, tg_0_xyyyzzz_1, tg_0_xyyzzz_1, \
                                         tg_0_xyyzzzz_0, tg_0_xyyzzzz_1, tg_0_xyzzzz_1, tg_0_xyzzzzz_0, tg_0_xyzzzzz_1, \
                                         tg_0_xzzzzz_1, tg_0_xzzzzzz_0, tg_0_xzzzzzz_1, tg_0_yyyyyy_1, tg_0_yyyyyyy_0, \
                                         tg_0_yyyyyyy_1, tg_0_yyyyyyz_0, tg_0_yyyyyyz_1, tg_0_yyyyyz_1, tg_0_yyyyyzz_0, \
                                         tg_0_yyyyyzz_1, tg_0_yyyyzz_1, tg_0_yyyyzzz_0, tg_0_yyyyzzz_1, tg_0_yyyzzz_1, \
                                         tg_0_yyyzzzz_0, tg_0_yyyzzzz_1, tg_0_yyzzzz_1, tg_0_yyzzzzz_0, tg_0_yyzzzzz_1, \
                                         tg_0_yzzzzz_1, tg_0_yzzzzzz_0, tg_0_yzzzzzz_1, tg_0_zzzzzz_1, tg_0_zzzzzzz_0, \
                                         tg_0_zzzzzzz_1, tg_y_xxyyzzz_0, tg_y_xxyzzzz_0, tg_y_xxzzzzz_0, tg_y_xyyyyyy_0, \
                                         tg_y_xyyyyyz_0, tg_y_xyyyyzz_0, tg_y_xyyyzzz_0, tg_y_xyyzzzz_0, tg_y_xyzzzzz_0, \
                                         tg_y_xzzzzzz_0, tg_y_yyyyyyy_0, tg_y_yyyyyyz_0, tg_y_yyyyyzz_0, tg_y_yyyyzzz_0, \
                                         tg_y_yyyzzzz_0, tg_y_yyzzzzz_0, tg_y_yzzzzzz_0, tg_y_zzzzzzz_0, tg_z_xxxxxxx_0, \
                                         tg_z_xxxxxxy_0, tg_z_xxxxxxz_0, tg_z_xxxxxyy_0, tg_z_xxxxxyz_0, tg_z_xxxxxzz_0, \
                                         tg_z_xxxxyyy_0, tg_z_xxxxyyz_0, tg_z_xxxxyzz_0, tg_z_xxxxzzz_0, tg_z_xxxyyyy_0, \
                                         tg_z_xxxyyyz_0, tg_z_xxxyyzz_0, tg_z_xxxyzzz_0, tg_z_xxxzzzz_0, tg_z_xxyyyyy_0, \
                                         tg_z_xxyyyyz_0, tg_z_xxyyyzz_0, tg_z_xxyyzzz_0, tg_z_xxyzzzz_0, tg_z_xxzzzzz_0, \
                                         tg_z_xyyyyyy_0, tg_z_xyyyyyz_0, tg_z_xyyyyzz_0, tg_z_xyyyzzz_0, tg_z_xyyzzzz_0, \
                                         tg_z_xyzzzzz_0, tg_z_xzzzzzz_0, tg_z_yyyyyyy_0, tg_z_yyyyyyz_0, tg_z_yyyyyzz_0, \
                                         tg_z_yyyyzzz_0, tg_z_yyyzzzz_0, tg_z_yyzzzzz_0, tg_z_yzzzzzz_0, tg_z_zzzzzzz_0, wp_y, \
                                         wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_y_xxyyzzz_0[j] = pb_y * tg_0_xxyyzzz_0[j] + wp_y[j] * tg_0_xxyyzzz_1[j] + fl1_fxn * tg_0_xxyzzz_1[j];

                    tg_y_xxyzzzz_0[j] = pb_y * tg_0_xxyzzzz_0[j] + wp_y[j] * tg_0_xxyzzzz_1[j] + 0.5 * fl1_fxn * tg_0_xxzzzz_1[j];

                    tg_y_xxzzzzz_0[j] = pb_y * tg_0_xxzzzzz_0[j] + wp_y[j] * tg_0_xxzzzzz_1[j];

                    tg_y_xyyyyyy_0[j] = pb_y * tg_0_xyyyyyy_0[j] + wp_y[j] * tg_0_xyyyyyy_1[j] + 3.0 * fl1_fxn * tg_0_xyyyyy_1[j];

                    tg_y_xyyyyyz_0[j] = pb_y * tg_0_xyyyyyz_0[j] + wp_y[j] * tg_0_xyyyyyz_1[j] + 2.5 * fl1_fxn * tg_0_xyyyyz_1[j];

                    tg_y_xyyyyzz_0[j] = pb_y * tg_0_xyyyyzz_0[j] + wp_y[j] * tg_0_xyyyyzz_1[j] + 2.0 * fl1_fxn * tg_0_xyyyzz_1[j];

                    tg_y_xyyyzzz_0[j] = pb_y * tg_0_xyyyzzz_0[j] + wp_y[j] * tg_0_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xyyzzz_1[j];

                    tg_y_xyyzzzz_0[j] = pb_y * tg_0_xyyzzzz_0[j] + wp_y[j] * tg_0_xyyzzzz_1[j] + fl1_fxn * tg_0_xyzzzz_1[j];

                    tg_y_xyzzzzz_0[j] = pb_y * tg_0_xyzzzzz_0[j] + wp_y[j] * tg_0_xyzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_xzzzzz_1[j];

                    tg_y_xzzzzzz_0[j] = pb_y * tg_0_xzzzzzz_0[j] + wp_y[j] * tg_0_xzzzzzz_1[j];

                    tg_y_yyyyyyy_0[j] = pb_y * tg_0_yyyyyyy_0[j] + wp_y[j] * tg_0_yyyyyyy_1[j] + 3.5 * fl1_fxn * tg_0_yyyyyy_1[j];

                    tg_y_yyyyyyz_0[j] = pb_y * tg_0_yyyyyyz_0[j] + wp_y[j] * tg_0_yyyyyyz_1[j] + 3.0 * fl1_fxn * tg_0_yyyyyz_1[j];

                    tg_y_yyyyyzz_0[j] = pb_y * tg_0_yyyyyzz_0[j] + wp_y[j] * tg_0_yyyyyzz_1[j] + 2.5 * fl1_fxn * tg_0_yyyyzz_1[j];

                    tg_y_yyyyzzz_0[j] = pb_y * tg_0_yyyyzzz_0[j] + wp_y[j] * tg_0_yyyyzzz_1[j] + 2.0 * fl1_fxn * tg_0_yyyzzz_1[j];

                    tg_y_yyyzzzz_0[j] = pb_y * tg_0_yyyzzzz_0[j] + wp_y[j] * tg_0_yyyzzzz_1[j] + 1.5 * fl1_fxn * tg_0_yyzzzz_1[j];

                    tg_y_yyzzzzz_0[j] = pb_y * tg_0_yyzzzzz_0[j] + wp_y[j] * tg_0_yyzzzzz_1[j] + fl1_fxn * tg_0_yzzzzz_1[j];

                    tg_y_yzzzzzz_0[j] = pb_y * tg_0_yzzzzzz_0[j] + wp_y[j] * tg_0_yzzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzzzz_1[j];

                    tg_y_zzzzzzz_0[j] = pb_y * tg_0_zzzzzzz_0[j] + wp_y[j] * tg_0_zzzzzzz_1[j];

                    tg_z_xxxxxxx_0[j] = pb_z * tg_0_xxxxxxx_0[j] + wp_z[j] * tg_0_xxxxxxx_1[j];

                    tg_z_xxxxxxy_0[j] = pb_z * tg_0_xxxxxxy_0[j] + wp_z[j] * tg_0_xxxxxxy_1[j];

                    tg_z_xxxxxxz_0[j] = pb_z * tg_0_xxxxxxz_0[j] + wp_z[j] * tg_0_xxxxxxz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxx_1[j];

                    tg_z_xxxxxyy_0[j] = pb_z * tg_0_xxxxxyy_0[j] + wp_z[j] * tg_0_xxxxxyy_1[j];

                    tg_z_xxxxxyz_0[j] = pb_z * tg_0_xxxxxyz_0[j] + wp_z[j] * tg_0_xxxxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxy_1[j];

                    tg_z_xxxxxzz_0[j] = pb_z * tg_0_xxxxxzz_0[j] + wp_z[j] * tg_0_xxxxxzz_1[j] + fl1_fxn * tg_0_xxxxxz_1[j];

                    tg_z_xxxxyyy_0[j] = pb_z * tg_0_xxxxyyy_0[j] + wp_z[j] * tg_0_xxxxyyy_1[j];

                    tg_z_xxxxyyz_0[j] = pb_z * tg_0_xxxxyyz_0[j] + wp_z[j] * tg_0_xxxxyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxyy_1[j];

                    tg_z_xxxxyzz_0[j] = pb_z * tg_0_xxxxyzz_0[j] + wp_z[j] * tg_0_xxxxyzz_1[j] + fl1_fxn * tg_0_xxxxyz_1[j];

                    tg_z_xxxxzzz_0[j] = pb_z * tg_0_xxxxzzz_0[j] + wp_z[j] * tg_0_xxxxzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxxxzz_1[j];

                    tg_z_xxxyyyy_0[j] = pb_z * tg_0_xxxyyyy_0[j] + wp_z[j] * tg_0_xxxyyyy_1[j];

                    tg_z_xxxyyyz_0[j] = pb_z * tg_0_xxxyyyz_0[j] + wp_z[j] * tg_0_xxxyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxyyy_1[j];

                    tg_z_xxxyyzz_0[j] = pb_z * tg_0_xxxyyzz_0[j] + wp_z[j] * tg_0_xxxyyzz_1[j] + fl1_fxn * tg_0_xxxyyz_1[j];

                    tg_z_xxxyzzz_0[j] = pb_z * tg_0_xxxyzzz_0[j] + wp_z[j] * tg_0_xxxyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxxyzz_1[j];

                    tg_z_xxxzzzz_0[j] = pb_z * tg_0_xxxzzzz_0[j] + wp_z[j] * tg_0_xxxzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxzzz_1[j];

                    tg_z_xxyyyyy_0[j] = pb_z * tg_0_xxyyyyy_0[j] + wp_z[j] * tg_0_xxyyyyy_1[j];

                    tg_z_xxyyyyz_0[j] = pb_z * tg_0_xxyyyyz_0[j] + wp_z[j] * tg_0_xxyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxyyyy_1[j];

                    tg_z_xxyyyzz_0[j] = pb_z * tg_0_xxyyyzz_0[j] + wp_z[j] * tg_0_xxyyyzz_1[j] + fl1_fxn * tg_0_xxyyyz_1[j];

                    tg_z_xxyyzzz_0[j] = pb_z * tg_0_xxyyzzz_0[j] + wp_z[j] * tg_0_xxyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyzz_1[j];

                    tg_z_xxyzzzz_0[j] = pb_z * tg_0_xxyzzzz_0[j] + wp_z[j] * tg_0_xxyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxyzzz_1[j];

                    tg_z_xxzzzzz_0[j] = pb_z * tg_0_xxzzzzz_0[j] + wp_z[j] * tg_0_xxzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_xxzzzz_1[j];

                    tg_z_xyyyyyy_0[j] = pb_z * tg_0_xyyyyyy_0[j] + wp_z[j] * tg_0_xyyyyyy_1[j];

                    tg_z_xyyyyyz_0[j] = pb_z * tg_0_xyyyyyz_0[j] + wp_z[j] * tg_0_xyyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xyyyyy_1[j];

                    tg_z_xyyyyzz_0[j] = pb_z * tg_0_xyyyyzz_0[j] + wp_z[j] * tg_0_xyyyyzz_1[j] + fl1_fxn * tg_0_xyyyyz_1[j];

                    tg_z_xyyyzzz_0[j] = pb_z * tg_0_xyyyzzz_0[j] + wp_z[j] * tg_0_xyyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xyyyzz_1[j];

                    tg_z_xyyzzzz_0[j] = pb_z * tg_0_xyyzzzz_0[j] + wp_z[j] * tg_0_xyyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xyyzzz_1[j];

                    tg_z_xyzzzzz_0[j] = pb_z * tg_0_xyzzzzz_0[j] + wp_z[j] * tg_0_xyzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_xyzzzz_1[j];

                    tg_z_xzzzzzz_0[j] = pb_z * tg_0_xzzzzzz_0[j] + wp_z[j] * tg_0_xzzzzzz_1[j] + 3.0 * fl1_fxn * tg_0_xzzzzz_1[j];

                    tg_z_yyyyyyy_0[j] = pb_z * tg_0_yyyyyyy_0[j] + wp_z[j] * tg_0_yyyyyyy_1[j];

                    tg_z_yyyyyyz_0[j] = pb_z * tg_0_yyyyyyz_0[j] + wp_z[j] * tg_0_yyyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyyy_1[j];

                    tg_z_yyyyyzz_0[j] = pb_z * tg_0_yyyyyzz_0[j] + wp_z[j] * tg_0_yyyyyzz_1[j] + fl1_fxn * tg_0_yyyyyz_1[j];

                    tg_z_yyyyzzz_0[j] = pb_z * tg_0_yyyyzzz_0[j] + wp_z[j] * tg_0_yyyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_yyyyzz_1[j];

                    tg_z_yyyzzzz_0[j] = pb_z * tg_0_yyyzzzz_0[j] + wp_z[j] * tg_0_yyyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_yyyzzz_1[j];

                    tg_z_yyzzzzz_0[j] = pb_z * tg_0_yyzzzzz_0[j] + wp_z[j] * tg_0_yyzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_yyzzzz_1[j];

                    tg_z_yzzzzzz_0[j] = pb_z * tg_0_yzzzzzz_0[j] + wp_z[j] * tg_0_yzzzzzz_1[j] + 3.0 * fl1_fxn * tg_0_yzzzzz_1[j];

                    tg_z_zzzzzzz_0[j] = pb_z * tg_0_zzzzzzz_0[j] + wp_z[j] * tg_0_zzzzzzz_1[j] + 3.5 * fl1_fxn * tg_0_zzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSKSP_0_54(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSKSP_54_108(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForSKSP_0_54(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,54)

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
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xxxxxx_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx); 

                auto tg_xxxxxx_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 1); 

                auto tg_xxxxxx_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 2); 

                auto tg_xxxxxy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 3); 

                auto tg_xxxxxy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 4); 

                auto tg_xxxxxy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 5); 

                auto tg_xxxxxz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 6); 

                auto tg_xxxxxz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 7); 

                auto tg_xxxxxz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 8); 

                auto tg_xxxxyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 9); 

                auto tg_xxxxyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 10); 

                auto tg_xxxxyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 11); 

                auto tg_xxxxyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 12); 

                auto tg_xxxxyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 13); 

                auto tg_xxxxyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 14); 

                auto tg_xxxxzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 15); 

                auto tg_xxxxzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 16); 

                auto tg_xxxxzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 17); 

                auto tg_xxxyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 18); 

                auto tg_xxxyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 19); 

                auto tg_xxxyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 20); 

                auto tg_xxxyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 21); 

                auto tg_xxxyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 22); 

                auto tg_xxxyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 23); 

                auto tg_xxxyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 24); 

                auto tg_xxxyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 25); 

                auto tg_xxxyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 26); 

                auto tg_xxxzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 27); 

                auto tg_xxxzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 28); 

                auto tg_xxxzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 29); 

                auto tg_xxyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 30); 

                auto tg_xxyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 31); 

                auto tg_xxyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 32); 

                auto tg_xxyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 33); 

                auto tg_xxyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 34); 

                auto tg_xxyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 35); 

                auto tg_xxyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 36); 

                auto tg_xxyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 37); 

                auto tg_xxyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 38); 

                auto tg_xxyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 39); 

                auto tg_xxyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 40); 

                auto tg_xxyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 41); 

                auto tg_xxzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 42); 

                auto tg_xxzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 43); 

                auto tg_xxzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 44); 

                auto tg_xyyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 45); 

                auto tg_xyyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 46); 

                auto tg_xyyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 47); 

                auto tg_xyyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 48); 

                auto tg_xyyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 49); 

                auto tg_xyyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 50); 

                auto tg_xyyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 51); 

                auto tg_xyyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 52); 

                auto tg_xyyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 53); 

                auto tg_xxxxxx_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx); 

                auto tg_xxxxxx_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 1); 

                auto tg_xxxxxx_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 2); 

                auto tg_xxxxxy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 3); 

                auto tg_xxxxxy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 4); 

                auto tg_xxxxxy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 5); 

                auto tg_xxxxxz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 6); 

                auto tg_xxxxxz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 7); 

                auto tg_xxxxxz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 8); 

                auto tg_xxxxyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 9); 

                auto tg_xxxxyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 10); 

                auto tg_xxxxyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 11); 

                auto tg_xxxxyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 12); 

                auto tg_xxxxyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 13); 

                auto tg_xxxxyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 14); 

                auto tg_xxxxzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 15); 

                auto tg_xxxxzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 16); 

                auto tg_xxxxzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 17); 

                auto tg_xxxyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 18); 

                auto tg_xxxyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 19); 

                auto tg_xxxyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 20); 

                auto tg_xxxyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 21); 

                auto tg_xxxyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 22); 

                auto tg_xxxyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 23); 

                auto tg_xxxyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 24); 

                auto tg_xxxyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 25); 

                auto tg_xxxyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 26); 

                auto tg_xxxzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 27); 

                auto tg_xxxzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 28); 

                auto tg_xxxzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 29); 

                auto tg_xxyyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 30); 

                auto tg_xxyyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 31); 

                auto tg_xxyyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 32); 

                auto tg_xxyyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 33); 

                auto tg_xxyyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 34); 

                auto tg_xxyyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 35); 

                auto tg_xxyyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 36); 

                auto tg_xxyyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 37); 

                auto tg_xxyyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 38); 

                auto tg_xxyzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 39); 

                auto tg_xxyzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 40); 

                auto tg_xxyzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 41); 

                auto tg_xxzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 42); 

                auto tg_xxzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 43); 

                auto tg_xxzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 44); 

                auto tg_xyyyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 45); 

                auto tg_xyyyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 46); 

                auto tg_xyyyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 47); 

                auto tg_xyyyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 48); 

                auto tg_xyyyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 49); 

                auto tg_xyyyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 50); 

                auto tg_xyyyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 51); 

                auto tg_xyyyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 52); 

                auto tg_xyyyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 53); 

                auto tg_xxxxx_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx); 

                auto tg_xxxxx_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 1); 

                auto tg_xxxxx_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 2); 

                auto tg_xxxxy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 3); 

                auto tg_xxxxy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 4); 

                auto tg_xxxxy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 5); 

                auto tg_xxxxz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 6); 

                auto tg_xxxxz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 7); 

                auto tg_xxxxz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 8); 

                auto tg_xxxyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 9); 

                auto tg_xxxyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 10); 

                auto tg_xxxyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 11); 

                auto tg_xxxyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 12); 

                auto tg_xxxyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 13); 

                auto tg_xxxyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 14); 

                auto tg_xxxzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 15); 

                auto tg_xxxzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 16); 

                auto tg_xxxzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 17); 

                auto tg_xxyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 18); 

                auto tg_xxyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 19); 

                auto tg_xxyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 20); 

                auto tg_xxyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 21); 

                auto tg_xxyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 22); 

                auto tg_xxyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 23); 

                auto tg_xxyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 24); 

                auto tg_xxyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 25); 

                auto tg_xxyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 26); 

                auto tg_xxzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 27); 

                auto tg_xxzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 28); 

                auto tg_xxzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 29); 

                auto tg_xyyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 30); 

                auto tg_xyyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 31); 

                auto tg_xyyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 32); 

                auto tg_xyyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 33); 

                auto tg_xyyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 34); 

                auto tg_xyyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 35); 

                auto tg_xyyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 36); 

                auto tg_xyyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 37); 

                auto tg_xyyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 38); 

                auto tg_xyzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 39); 

                auto tg_xyzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 40); 

                auto tg_xyzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 41); 

                auto tg_xzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 42); 

                auto tg_xzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 43); 

                auto tg_xzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 44); 

                auto tg_yyyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 45); 

                auto tg_yyyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 46); 

                auto tg_yyyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 47); 

                auto tg_yyyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 48); 

                auto tg_yyyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 49); 

                auto tg_yyyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 50); 

                auto tg_yyyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 51); 

                auto tg_yyyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 52); 

                auto tg_yyyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 53); 

                auto tg_xxxxx_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx); 

                auto tg_xxxxx_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 1); 

                auto tg_xxxxx_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 2); 

                auto tg_xxxxy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 3); 

                auto tg_xxxxy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 4); 

                auto tg_xxxxy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 5); 

                auto tg_xxxxz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 6); 

                auto tg_xxxxz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 7); 

                auto tg_xxxxz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 8); 

                auto tg_xxxyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 9); 

                auto tg_xxxyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 10); 

                auto tg_xxxyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 11); 

                auto tg_xxxyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 12); 

                auto tg_xxxyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 13); 

                auto tg_xxxyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 14); 

                auto tg_xxxzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 15); 

                auto tg_xxxzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 16); 

                auto tg_xxxzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 17); 

                auto tg_xxyyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 18); 

                auto tg_xxyyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 19); 

                auto tg_xxyyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 20); 

                auto tg_xxyyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 21); 

                auto tg_xxyyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 22); 

                auto tg_xxyyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 23); 

                auto tg_xxyzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 24); 

                auto tg_xxyzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 25); 

                auto tg_xxyzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 26); 

                auto tg_xxzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 27); 

                auto tg_xxzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 28); 

                auto tg_xxzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 29); 

                auto tg_xyyyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 30); 

                auto tg_xyyyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 31); 

                auto tg_xyyyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 32); 

                auto tg_xyyyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 33); 

                auto tg_xyyyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 34); 

                auto tg_xyyyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 35); 

                auto tg_xyyzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 36); 

                auto tg_xyyzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 37); 

                auto tg_xyyzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 38); 

                auto tg_xyzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 39); 

                auto tg_xyzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 40); 

                auto tg_xyzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 41); 

                auto tg_xzzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 42); 

                auto tg_xzzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 43); 

                auto tg_xzzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 44); 

                auto tg_yyyyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 45); 

                auto tg_yyyyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 46); 

                auto tg_yyyyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 47); 

                auto tg_yyyyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 48); 

                auto tg_yyyyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 49); 

                auto tg_yyyyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 50); 

                auto tg_yyyzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 51); 

                auto tg_yyyzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 52); 

                auto tg_yyyzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 53); 

                auto tg_xxxxxx_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx); 

                auto tg_xxxxxy_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 1); 

                auto tg_xxxxxz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 2); 

                auto tg_xxxxyy_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 3); 

                auto tg_xxxxyz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 4); 

                auto tg_xxxxzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 5); 

                auto tg_xxxyyy_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 6); 

                auto tg_xxxyyz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 7); 

                auto tg_xxxyzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 8); 

                auto tg_xxxzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 9); 

                auto tg_xxyyyy_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 10); 

                auto tg_xxyyyz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 11); 

                auto tg_xxyyzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 12); 

                auto tg_xxyzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 13); 

                auto tg_xxzzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 14); 

                auto tg_xyyyyy_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 15); 

                auto tg_xyyyyz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 16); 

                auto tg_xyyyzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 17); 

                // set up pointers to integrals

                auto tg_xxxxxxx_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx); 

                auto tg_xxxxxxx_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 1); 

                auto tg_xxxxxxx_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 2); 

                auto tg_xxxxxxy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 3); 

                auto tg_xxxxxxy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 4); 

                auto tg_xxxxxxy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 5); 

                auto tg_xxxxxxz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 6); 

                auto tg_xxxxxxz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 7); 

                auto tg_xxxxxxz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 8); 

                auto tg_xxxxxyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 9); 

                auto tg_xxxxxyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 10); 

                auto tg_xxxxxyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 11); 

                auto tg_xxxxxyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 12); 

                auto tg_xxxxxyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 13); 

                auto tg_xxxxxyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 14); 

                auto tg_xxxxxzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 15); 

                auto tg_xxxxxzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 16); 

                auto tg_xxxxxzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 17); 

                auto tg_xxxxyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 18); 

                auto tg_xxxxyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 19); 

                auto tg_xxxxyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 20); 

                auto tg_xxxxyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 21); 

                auto tg_xxxxyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 22); 

                auto tg_xxxxyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 23); 

                auto tg_xxxxyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 24); 

                auto tg_xxxxyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 25); 

                auto tg_xxxxyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 26); 

                auto tg_xxxxzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 27); 

                auto tg_xxxxzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 28); 

                auto tg_xxxxzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 29); 

                auto tg_xxxyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 30); 

                auto tg_xxxyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 31); 

                auto tg_xxxyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 32); 

                auto tg_xxxyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 33); 

                auto tg_xxxyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 34); 

                auto tg_xxxyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 35); 

                auto tg_xxxyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 36); 

                auto tg_xxxyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 37); 

                auto tg_xxxyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 38); 

                auto tg_xxxyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 39); 

                auto tg_xxxyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 40); 

                auto tg_xxxyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 41); 

                auto tg_xxxzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 42); 

                auto tg_xxxzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 43); 

                auto tg_xxxzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 44); 

                auto tg_xxyyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 45); 

                auto tg_xxyyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 46); 

                auto tg_xxyyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 47); 

                auto tg_xxyyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 48); 

                auto tg_xxyyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 49); 

                auto tg_xxyyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 50); 

                auto tg_xxyyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 51); 

                auto tg_xxyyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 52); 

                auto tg_xxyyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 53); 

                // Batch of Integrals (0,54)

                #pragma omp simd aligned(fxn, fza, tg_xxxxx_x_0, tg_xxxxx_x_1, tg_xxxxx_y_0, tg_xxxxx_y_1, \
                                         tg_xxxxx_z_0, tg_xxxxx_z_1, tg_xxxxxx_0_1, tg_xxxxxx_x_0, tg_xxxxxx_x_1, \
                                         tg_xxxxxx_y_0, tg_xxxxxx_y_1, tg_xxxxxx_z_0, tg_xxxxxx_z_1, tg_xxxxxxx_x_0, \
                                         tg_xxxxxxx_y_0, tg_xxxxxxx_z_0, tg_xxxxxxy_x_0, tg_xxxxxxy_y_0, tg_xxxxxxy_z_0, \
                                         tg_xxxxxxz_x_0, tg_xxxxxxz_y_0, tg_xxxxxxz_z_0, tg_xxxxxy_0_1, tg_xxxxxy_x_0, \
                                         tg_xxxxxy_x_1, tg_xxxxxy_y_0, tg_xxxxxy_y_1, tg_xxxxxy_z_0, tg_xxxxxy_z_1, \
                                         tg_xxxxxyy_x_0, tg_xxxxxyy_y_0, tg_xxxxxyy_z_0, tg_xxxxxyz_x_0, tg_xxxxxyz_y_0, \
                                         tg_xxxxxyz_z_0, tg_xxxxxz_0_1, tg_xxxxxz_x_0, tg_xxxxxz_x_1, tg_xxxxxz_y_0, \
                                         tg_xxxxxz_y_1, tg_xxxxxz_z_0, tg_xxxxxz_z_1, tg_xxxxxzz_x_0, tg_xxxxxzz_y_0, \
                                         tg_xxxxxzz_z_0, tg_xxxxy_x_0, tg_xxxxy_x_1, tg_xxxxy_y_0, tg_xxxxy_y_1, tg_xxxxy_z_0, \
                                         tg_xxxxy_z_1, tg_xxxxyy_0_1, tg_xxxxyy_x_0, tg_xxxxyy_x_1, tg_xxxxyy_y_0, \
                                         tg_xxxxyy_y_1, tg_xxxxyy_z_0, tg_xxxxyy_z_1, tg_xxxxyyy_x_0, tg_xxxxyyy_y_0, \
                                         tg_xxxxyyy_z_0, tg_xxxxyyz_x_0, tg_xxxxyyz_y_0, tg_xxxxyyz_z_0, tg_xxxxyz_0_1, \
                                         tg_xxxxyz_x_0, tg_xxxxyz_x_1, tg_xxxxyz_y_0, tg_xxxxyz_y_1, tg_xxxxyz_z_0, \
                                         tg_xxxxyz_z_1, tg_xxxxyzz_x_0, tg_xxxxyzz_y_0, tg_xxxxyzz_z_0, tg_xxxxz_x_0, \
                                         tg_xxxxz_x_1, tg_xxxxz_y_0, tg_xxxxz_y_1, tg_xxxxz_z_0, tg_xxxxz_z_1, \
                                         tg_xxxxzz_0_1, tg_xxxxzz_x_0, tg_xxxxzz_x_1, tg_xxxxzz_y_0, tg_xxxxzz_y_1, \
                                         tg_xxxxzz_z_0, tg_xxxxzz_z_1, tg_xxxxzzz_x_0, tg_xxxxzzz_y_0, tg_xxxxzzz_z_0, \
                                         tg_xxxyy_x_0, tg_xxxyy_x_1, tg_xxxyy_y_0, tg_xxxyy_y_1, tg_xxxyy_z_0, tg_xxxyy_z_1, \
                                         tg_xxxyyy_0_1, tg_xxxyyy_x_0, tg_xxxyyy_x_1, tg_xxxyyy_y_0, tg_xxxyyy_y_1, \
                                         tg_xxxyyy_z_0, tg_xxxyyy_z_1, tg_xxxyyyy_x_0, tg_xxxyyyy_y_0, tg_xxxyyyy_z_0, \
                                         tg_xxxyyyz_x_0, tg_xxxyyyz_y_0, tg_xxxyyyz_z_0, tg_xxxyyz_0_1, tg_xxxyyz_x_0, \
                                         tg_xxxyyz_x_1, tg_xxxyyz_y_0, tg_xxxyyz_y_1, tg_xxxyyz_z_0, tg_xxxyyz_z_1, \
                                         tg_xxxyyzz_x_0, tg_xxxyyzz_y_0, tg_xxxyyzz_z_0, tg_xxxyz_x_0, tg_xxxyz_x_1, \
                                         tg_xxxyz_y_0, tg_xxxyz_y_1, tg_xxxyz_z_0, tg_xxxyz_z_1, tg_xxxyzz_0_1, \
                                         tg_xxxyzz_x_0, tg_xxxyzz_x_1, tg_xxxyzz_y_0, tg_xxxyzz_y_1, tg_xxxyzz_z_0, \
                                         tg_xxxyzz_z_1, tg_xxxyzzz_x_0, tg_xxxyzzz_y_0, tg_xxxyzzz_z_0, tg_xxxzz_x_0, \
                                         tg_xxxzz_x_1, tg_xxxzz_y_0, tg_xxxzz_y_1, tg_xxxzz_z_0, tg_xxxzz_z_1, \
                                         tg_xxxzzz_0_1, tg_xxxzzz_x_0, tg_xxxzzz_x_1, tg_xxxzzz_y_0, tg_xxxzzz_y_1, \
                                         tg_xxxzzz_z_0, tg_xxxzzz_z_1, tg_xxxzzzz_x_0, tg_xxxzzzz_y_0, tg_xxxzzzz_z_0, \
                                         tg_xxyyy_x_0, tg_xxyyy_x_1, tg_xxyyy_y_0, tg_xxyyy_y_1, tg_xxyyy_z_0, tg_xxyyy_z_1, \
                                         tg_xxyyyy_0_1, tg_xxyyyy_x_0, tg_xxyyyy_x_1, tg_xxyyyy_y_0, tg_xxyyyy_y_1, \
                                         tg_xxyyyy_z_0, tg_xxyyyy_z_1, tg_xxyyyyy_x_0, tg_xxyyyyy_y_0, tg_xxyyyyy_z_0, \
                                         tg_xxyyyyz_x_0, tg_xxyyyyz_y_0, tg_xxyyyyz_z_0, tg_xxyyyz_0_1, tg_xxyyyz_x_0, \
                                         tg_xxyyyz_x_1, tg_xxyyyz_y_0, tg_xxyyyz_y_1, tg_xxyyyz_z_0, tg_xxyyyz_z_1, \
                                         tg_xxyyyzz_x_0, tg_xxyyyzz_y_0, tg_xxyyyzz_z_0, tg_xxyyz_x_0, tg_xxyyz_x_1, \
                                         tg_xxyyz_y_0, tg_xxyyz_y_1, tg_xxyyz_z_0, tg_xxyyz_z_1, tg_xxyyzz_0_1, \
                                         tg_xxyyzz_x_0, tg_xxyyzz_x_1, tg_xxyyzz_y_0, tg_xxyyzz_y_1, tg_xxyyzz_z_0, \
                                         tg_xxyyzz_z_1, tg_xxyzz_x_0, tg_xxyzz_x_1, tg_xxyzz_y_0, tg_xxyzz_y_1, tg_xxyzz_z_0, \
                                         tg_xxyzz_z_1, tg_xxyzzz_0_1, tg_xxyzzz_x_0, tg_xxyzzz_x_1, tg_xxyzzz_y_0, \
                                         tg_xxyzzz_y_1, tg_xxyzzz_z_0, tg_xxyzzz_z_1, tg_xxzzz_x_0, tg_xxzzz_x_1, \
                                         tg_xxzzz_y_0, tg_xxzzz_y_1, tg_xxzzz_z_0, tg_xxzzz_z_1, tg_xxzzzz_0_1, \
                                         tg_xxzzzz_x_0, tg_xxzzzz_x_1, tg_xxzzzz_y_0, tg_xxzzzz_y_1, tg_xxzzzz_z_0, \
                                         tg_xxzzzz_z_1, tg_xyyyy_x_0, tg_xyyyy_x_1, tg_xyyyy_y_0, tg_xyyyy_y_1, tg_xyyyy_z_0, \
                                         tg_xyyyy_z_1, tg_xyyyyy_0_1, tg_xyyyyy_x_0, tg_xyyyyy_x_1, tg_xyyyyy_y_0, \
                                         tg_xyyyyy_y_1, tg_xyyyyy_z_0, tg_xyyyyy_z_1, tg_xyyyyz_0_1, tg_xyyyyz_x_0, \
                                         tg_xyyyyz_x_1, tg_xyyyyz_y_0, tg_xyyyyz_y_1, tg_xyyyyz_z_0, tg_xyyyyz_z_1, \
                                         tg_xyyyz_x_0, tg_xyyyz_x_1, tg_xyyyz_y_0, tg_xyyyz_y_1, tg_xyyyz_z_0, tg_xyyyz_z_1, \
                                         tg_xyyyzz_0_1, tg_xyyyzz_x_0, tg_xyyyzz_x_1, tg_xyyyzz_y_0, tg_xyyyzz_y_1, \
                                         tg_xyyyzz_z_0, tg_xyyyzz_z_1, tg_xyyzz_x_0, tg_xyyzz_x_1, tg_xyyzz_y_0, \
                                         tg_xyyzz_y_1, tg_xyyzz_z_0, tg_xyyzz_z_1, tg_xyzzz_x_0, tg_xyzzz_x_1, tg_xyzzz_y_0, \
                                         tg_xyzzz_y_1, tg_xyzzz_z_0, tg_xyzzz_z_1, tg_xzzzz_x_0, tg_xzzzz_x_1, tg_xzzzz_y_0, \
                                         tg_xzzzz_y_1, tg_xzzzz_z_0, tg_xzzzz_z_1, tg_yyyyy_x_0, tg_yyyyy_x_1, tg_yyyyy_y_0, \
                                         tg_yyyyy_y_1, tg_yyyyy_z_0, tg_yyyyy_z_1, tg_yyyyz_x_0, tg_yyyyz_x_1, tg_yyyyz_y_0, \
                                         tg_yyyyz_y_1, tg_yyyyz_z_0, tg_yyyyz_z_1, tg_yyyzz_x_0, tg_yyyzz_x_1, tg_yyyzz_y_0, \
                                         tg_yyyzz_y_1, tg_yyyzz_z_0, tg_yyyzz_z_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxxx_x_0[j] = pb_x * tg_xxxxxx_x_0[j] + wp_x[j] * tg_xxxxxx_x_1[j] + 3.0 * fl1_fx * tg_xxxxx_x_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxx_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxx_0_1[j];

                    tg_xxxxxxx_y_0[j] = pb_x * tg_xxxxxx_y_0[j] + wp_x[j] * tg_xxxxxx_y_1[j] + 3.0 * fl1_fx * tg_xxxxx_y_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxx_y_1[j];

                    tg_xxxxxxx_z_0[j] = pb_x * tg_xxxxxx_z_0[j] + wp_x[j] * tg_xxxxxx_z_1[j] + 3.0 * fl1_fx * tg_xxxxx_z_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxx_z_1[j];

                    tg_xxxxxxy_x_0[j] = pb_x * tg_xxxxxy_x_0[j] + wp_x[j] * tg_xxxxxy_x_1[j] + 2.5 * fl1_fx * tg_xxxxy_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxy_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxy_0_1[j];

                    tg_xxxxxxy_y_0[j] = pb_x * tg_xxxxxy_y_0[j] + wp_x[j] * tg_xxxxxy_y_1[j] + 2.5 * fl1_fx * tg_xxxxy_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxy_y_1[j];

                    tg_xxxxxxy_z_0[j] = pb_x * tg_xxxxxy_z_0[j] + wp_x[j] * tg_xxxxxy_z_1[j] + 2.5 * fl1_fx * tg_xxxxy_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxy_z_1[j];

                    tg_xxxxxxz_x_0[j] = pb_x * tg_xxxxxz_x_0[j] + wp_x[j] * tg_xxxxxz_x_1[j] + 2.5 * fl1_fx * tg_xxxxz_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxz_0_1[j];

                    tg_xxxxxxz_y_0[j] = pb_x * tg_xxxxxz_y_0[j] + wp_x[j] * tg_xxxxxz_y_1[j] + 2.5 * fl1_fx * tg_xxxxz_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxz_y_1[j];

                    tg_xxxxxxz_z_0[j] = pb_x * tg_xxxxxz_z_0[j] + wp_x[j] * tg_xxxxxz_z_1[j] + 2.5 * fl1_fx * tg_xxxxz_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxz_z_1[j];

                    tg_xxxxxyy_x_0[j] = pb_x * tg_xxxxyy_x_0[j] + wp_x[j] * tg_xxxxyy_x_1[j] + 2.0 * fl1_fx * tg_xxxyy_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyy_x_1[j] + 0.5 * fl1_fxn * tg_xxxxyy_0_1[j];

                    tg_xxxxxyy_y_0[j] = pb_x * tg_xxxxyy_y_0[j] + wp_x[j] * tg_xxxxyy_y_1[j] + 2.0 * fl1_fx * tg_xxxyy_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyy_y_1[j];

                    tg_xxxxxyy_z_0[j] = pb_x * tg_xxxxyy_z_0[j] + wp_x[j] * tg_xxxxyy_z_1[j] + 2.0 * fl1_fx * tg_xxxyy_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyy_z_1[j];

                    tg_xxxxxyz_x_0[j] = pb_x * tg_xxxxyz_x_0[j] + wp_x[j] * tg_xxxxyz_x_1[j] + 2.0 * fl1_fx * tg_xxxyz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxyz_0_1[j];

                    tg_xxxxxyz_y_0[j] = pb_x * tg_xxxxyz_y_0[j] + wp_x[j] * tg_xxxxyz_y_1[j] + 2.0 * fl1_fx * tg_xxxyz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyz_y_1[j];

                    tg_xxxxxyz_z_0[j] = pb_x * tg_xxxxyz_z_0[j] + wp_x[j] * tg_xxxxyz_z_1[j] + 2.0 * fl1_fx * tg_xxxyz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyz_z_1[j];

                    tg_xxxxxzz_x_0[j] = pb_x * tg_xxxxzz_x_0[j] + wp_x[j] * tg_xxxxzz_x_1[j] + 2.0 * fl1_fx * tg_xxxzz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxzz_0_1[j];

                    tg_xxxxxzz_y_0[j] = pb_x * tg_xxxxzz_y_0[j] + wp_x[j] * tg_xxxxzz_y_1[j] + 2.0 * fl1_fx * tg_xxxzz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzz_y_1[j];

                    tg_xxxxxzz_z_0[j] = pb_x * tg_xxxxzz_z_0[j] + wp_x[j] * tg_xxxxzz_z_1[j] + 2.0 * fl1_fx * tg_xxxzz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzz_z_1[j];

                    tg_xxxxyyy_x_0[j] = pb_x * tg_xxxyyy_x_0[j] + wp_x[j] * tg_xxxyyy_x_1[j] + 1.5 * fl1_fx * tg_xxyyy_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyy_x_1[j] + 0.5 * fl1_fxn * tg_xxxyyy_0_1[j];

                    tg_xxxxyyy_y_0[j] = pb_x * tg_xxxyyy_y_0[j] + wp_x[j] * tg_xxxyyy_y_1[j] + 1.5 * fl1_fx * tg_xxyyy_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyy_y_1[j];

                    tg_xxxxyyy_z_0[j] = pb_x * tg_xxxyyy_z_0[j] + wp_x[j] * tg_xxxyyy_z_1[j] + 1.5 * fl1_fx * tg_xxyyy_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyy_z_1[j];

                    tg_xxxxyyz_x_0[j] = pb_x * tg_xxxyyz_x_0[j] + wp_x[j] * tg_xxxyyz_x_1[j] + 1.5 * fl1_fx * tg_xxyyz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyz_x_1[j] + 0.5 * fl1_fxn * tg_xxxyyz_0_1[j];

                    tg_xxxxyyz_y_0[j] = pb_x * tg_xxxyyz_y_0[j] + wp_x[j] * tg_xxxyyz_y_1[j] + 1.5 * fl1_fx * tg_xxyyz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyz_y_1[j];

                    tg_xxxxyyz_z_0[j] = pb_x * tg_xxxyyz_z_0[j] + wp_x[j] * tg_xxxyyz_z_1[j] + 1.5 * fl1_fx * tg_xxyyz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyz_z_1[j];

                    tg_xxxxyzz_x_0[j] = pb_x * tg_xxxyzz_x_0[j] + wp_x[j] * tg_xxxyzz_x_1[j] + 1.5 * fl1_fx * tg_xxyzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxyzz_0_1[j];

                    tg_xxxxyzz_y_0[j] = pb_x * tg_xxxyzz_y_0[j] + wp_x[j] * tg_xxxyzz_y_1[j] + 1.5 * fl1_fx * tg_xxyzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzz_y_1[j];

                    tg_xxxxyzz_z_0[j] = pb_x * tg_xxxyzz_z_0[j] + wp_x[j] * tg_xxxyzz_z_1[j] + 1.5 * fl1_fx * tg_xxyzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzz_z_1[j];

                    tg_xxxxzzz_x_0[j] = pb_x * tg_xxxzzz_x_0[j] + wp_x[j] * tg_xxxzzz_x_1[j] + 1.5 * fl1_fx * tg_xxzzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxzzz_0_1[j];

                    tg_xxxxzzz_y_0[j] = pb_x * tg_xxxzzz_y_0[j] + wp_x[j] * tg_xxxzzz_y_1[j] + 1.5 * fl1_fx * tg_xxzzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzz_y_1[j];

                    tg_xxxxzzz_z_0[j] = pb_x * tg_xxxzzz_z_0[j] + wp_x[j] * tg_xxxzzz_z_1[j] + 1.5 * fl1_fx * tg_xxzzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzz_z_1[j];

                    tg_xxxyyyy_x_0[j] = pb_x * tg_xxyyyy_x_0[j] + wp_x[j] * tg_xxyyyy_x_1[j] + fl1_fx * tg_xyyyy_x_0[j] - fl1_fx * fl1_fza * tg_xyyyy_x_1[j] + 0.5 * fl1_fxn * tg_xxyyyy_0_1[j];

                    tg_xxxyyyy_y_0[j] = pb_x * tg_xxyyyy_y_0[j] + wp_x[j] * tg_xxyyyy_y_1[j] + fl1_fx * tg_xyyyy_y_0[j] - fl1_fx * fl1_fza * tg_xyyyy_y_1[j];

                    tg_xxxyyyy_z_0[j] = pb_x * tg_xxyyyy_z_0[j] + wp_x[j] * tg_xxyyyy_z_1[j] + fl1_fx * tg_xyyyy_z_0[j] - fl1_fx * fl1_fza * tg_xyyyy_z_1[j];

                    tg_xxxyyyz_x_0[j] = pb_x * tg_xxyyyz_x_0[j] + wp_x[j] * tg_xxyyyz_x_1[j] + fl1_fx * tg_xyyyz_x_0[j] - fl1_fx * fl1_fza * tg_xyyyz_x_1[j] + 0.5 * fl1_fxn * tg_xxyyyz_0_1[j];

                    tg_xxxyyyz_y_0[j] = pb_x * tg_xxyyyz_y_0[j] + wp_x[j] * tg_xxyyyz_y_1[j] + fl1_fx * tg_xyyyz_y_0[j] - fl1_fx * fl1_fza * tg_xyyyz_y_1[j];

                    tg_xxxyyyz_z_0[j] = pb_x * tg_xxyyyz_z_0[j] + wp_x[j] * tg_xxyyyz_z_1[j] + fl1_fx * tg_xyyyz_z_0[j] - fl1_fx * fl1_fza * tg_xyyyz_z_1[j];

                    tg_xxxyyzz_x_0[j] = pb_x * tg_xxyyzz_x_0[j] + wp_x[j] * tg_xxyyzz_x_1[j] + fl1_fx * tg_xyyzz_x_0[j] - fl1_fx * fl1_fza * tg_xyyzz_x_1[j] + 0.5 * fl1_fxn * tg_xxyyzz_0_1[j];

                    tg_xxxyyzz_y_0[j] = pb_x * tg_xxyyzz_y_0[j] + wp_x[j] * tg_xxyyzz_y_1[j] + fl1_fx * tg_xyyzz_y_0[j] - fl1_fx * fl1_fza * tg_xyyzz_y_1[j];

                    tg_xxxyyzz_z_0[j] = pb_x * tg_xxyyzz_z_0[j] + wp_x[j] * tg_xxyyzz_z_1[j] + fl1_fx * tg_xyyzz_z_0[j] - fl1_fx * fl1_fza * tg_xyyzz_z_1[j];

                    tg_xxxyzzz_x_0[j] = pb_x * tg_xxyzzz_x_0[j] + wp_x[j] * tg_xxyzzz_x_1[j] + fl1_fx * tg_xyzzz_x_0[j] - fl1_fx * fl1_fza * tg_xyzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxyzzz_0_1[j];

                    tg_xxxyzzz_y_0[j] = pb_x * tg_xxyzzz_y_0[j] + wp_x[j] * tg_xxyzzz_y_1[j] + fl1_fx * tg_xyzzz_y_0[j] - fl1_fx * fl1_fza * tg_xyzzz_y_1[j];

                    tg_xxxyzzz_z_0[j] = pb_x * tg_xxyzzz_z_0[j] + wp_x[j] * tg_xxyzzz_z_1[j] + fl1_fx * tg_xyzzz_z_0[j] - fl1_fx * fl1_fza * tg_xyzzz_z_1[j];

                    tg_xxxzzzz_x_0[j] = pb_x * tg_xxzzzz_x_0[j] + wp_x[j] * tg_xxzzzz_x_1[j] + fl1_fx * tg_xzzzz_x_0[j] - fl1_fx * fl1_fza * tg_xzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxzzzz_0_1[j];

                    tg_xxxzzzz_y_0[j] = pb_x * tg_xxzzzz_y_0[j] + wp_x[j] * tg_xxzzzz_y_1[j] + fl1_fx * tg_xzzzz_y_0[j] - fl1_fx * fl1_fza * tg_xzzzz_y_1[j];

                    tg_xxxzzzz_z_0[j] = pb_x * tg_xxzzzz_z_0[j] + wp_x[j] * tg_xxzzzz_z_1[j] + fl1_fx * tg_xzzzz_z_0[j] - fl1_fx * fl1_fza * tg_xzzzz_z_1[j];

                    tg_xxyyyyy_x_0[j] = pb_x * tg_xyyyyy_x_0[j] + wp_x[j] * tg_xyyyyy_x_1[j] + 0.5 * fl1_fx * tg_yyyyy_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyy_x_1[j] + 0.5 * fl1_fxn * tg_xyyyyy_0_1[j];

                    tg_xxyyyyy_y_0[j] = pb_x * tg_xyyyyy_y_0[j] + wp_x[j] * tg_xyyyyy_y_1[j] + 0.5 * fl1_fx * tg_yyyyy_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyy_y_1[j];

                    tg_xxyyyyy_z_0[j] = pb_x * tg_xyyyyy_z_0[j] + wp_x[j] * tg_xyyyyy_z_1[j] + 0.5 * fl1_fx * tg_yyyyy_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyy_z_1[j];

                    tg_xxyyyyz_x_0[j] = pb_x * tg_xyyyyz_x_0[j] + wp_x[j] * tg_xyyyyz_x_1[j] + 0.5 * fl1_fx * tg_yyyyz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyz_x_1[j] + 0.5 * fl1_fxn * tg_xyyyyz_0_1[j];

                    tg_xxyyyyz_y_0[j] = pb_x * tg_xyyyyz_y_0[j] + wp_x[j] * tg_xyyyyz_y_1[j] + 0.5 * fl1_fx * tg_yyyyz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyz_y_1[j];

                    tg_xxyyyyz_z_0[j] = pb_x * tg_xyyyyz_z_0[j] + wp_x[j] * tg_xyyyyz_z_1[j] + 0.5 * fl1_fx * tg_yyyyz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyz_z_1[j];

                    tg_xxyyyzz_x_0[j] = pb_x * tg_xyyyzz_x_0[j] + wp_x[j] * tg_xyyyzz_x_1[j] + 0.5 * fl1_fx * tg_yyyzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzz_x_1[j] + 0.5 * fl1_fxn * tg_xyyyzz_0_1[j];

                    tg_xxyyyzz_y_0[j] = pb_x * tg_xyyyzz_y_0[j] + wp_x[j] * tg_xyyyzz_y_1[j] + 0.5 * fl1_fx * tg_yyyzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzz_y_1[j];

                    tg_xxyyyzz_z_0[j] = pb_x * tg_xyyyzz_z_0[j] + wp_x[j] * tg_xyyyzz_z_1[j] + 0.5 * fl1_fx * tg_yyyzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzz_z_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSP_54_108(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (54,108)

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
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_5_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_5_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {5, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xyyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 54); 

                auto tg_xyyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 55); 

                auto tg_xyyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 56); 

                auto tg_xyzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 57); 

                auto tg_xyzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 58); 

                auto tg_xyzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 59); 

                auto tg_xzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 60); 

                auto tg_xzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 61); 

                auto tg_xzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 62); 

                auto tg_yyyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 63); 

                auto tg_yyyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 64); 

                auto tg_yyyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 65); 

                auto tg_yyyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 66); 

                auto tg_yyyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 67); 

                auto tg_yyyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 68); 

                auto tg_yyyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 69); 

                auto tg_yyyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 70); 

                auto tg_yyyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 71); 

                auto tg_yyyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 72); 

                auto tg_yyyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 73); 

                auto tg_yyyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 74); 

                auto tg_yyzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 75); 

                auto tg_yyzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 76); 

                auto tg_yyzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 77); 

                auto tg_yzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 78); 

                auto tg_yzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 79); 

                auto tg_yzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 80); 

                auto tg_zzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 81); 

                auto tg_zzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 82); 

                auto tg_zzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 83); 

                auto tg_xyyzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 54); 

                auto tg_xyyzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 55); 

                auto tg_xyyzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 56); 

                auto tg_xyzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 57); 

                auto tg_xyzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 58); 

                auto tg_xyzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 59); 

                auto tg_xzzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 60); 

                auto tg_xzzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 61); 

                auto tg_xzzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 62); 

                auto tg_yyyyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 63); 

                auto tg_yyyyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 64); 

                auto tg_yyyyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 65); 

                auto tg_yyyyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 66); 

                auto tg_yyyyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 67); 

                auto tg_yyyyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 68); 

                auto tg_yyyyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 69); 

                auto tg_yyyyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 70); 

                auto tg_yyyyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 71); 

                auto tg_yyyzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 72); 

                auto tg_yyyzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 73); 

                auto tg_yyyzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 74); 

                auto tg_yyzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 75); 

                auto tg_yyzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 76); 

                auto tg_yyzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 77); 

                auto tg_yzzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 78); 

                auto tg_yzzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 79); 

                auto tg_yzzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 80); 

                auto tg_zzzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 81); 

                auto tg_zzzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 82); 

                auto tg_zzzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 83); 

                auto tg_yyyyy_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 45); 

                auto tg_yyyyy_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 46); 

                auto tg_yyyyy_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 47); 

                auto tg_yyyyz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 48); 

                auto tg_yyyyz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 49); 

                auto tg_yyyyz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 50); 

                auto tg_yyyzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 51); 

                auto tg_yyyzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 52); 

                auto tg_yyyzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 53); 

                auto tg_yyzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 54); 

                auto tg_yyzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 55); 

                auto tg_yyzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 56); 

                auto tg_yzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 57); 

                auto tg_yzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 58); 

                auto tg_yzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 59); 

                auto tg_zzzzz_x_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 60); 

                auto tg_zzzzz_y_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 61); 

                auto tg_zzzzz_z_0 = primBuffer.data(pidx_g_5_1_m0 + 63 * idx + 62); 

                auto tg_yyyyy_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 45); 

                auto tg_yyyyy_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 46); 

                auto tg_yyyyy_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 47); 

                auto tg_yyyyz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 48); 

                auto tg_yyyyz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 49); 

                auto tg_yyyyz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 50); 

                auto tg_yyyzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 51); 

                auto tg_yyyzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 52); 

                auto tg_yyyzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 53); 

                auto tg_yyzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 54); 

                auto tg_yyzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 55); 

                auto tg_yyzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 56); 

                auto tg_yzzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 57); 

                auto tg_yzzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 58); 

                auto tg_yzzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 59); 

                auto tg_zzzzz_x_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 60); 

                auto tg_zzzzz_y_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 61); 

                auto tg_zzzzz_z_1 = primBuffer.data(pidx_g_5_1_m1 + 63 * idx + 62); 

                auto tg_xyyzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 18); 

                auto tg_xyzzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 19); 

                auto tg_xzzzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 20); 

                auto tg_yyyyyy_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 21); 

                auto tg_yyyyyz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 22); 

                auto tg_yyyyzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 23); 

                auto tg_yyyzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 24); 

                auto tg_yyzzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 25); 

                auto tg_yzzzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 26); 

                auto tg_zzzzzz_0_1 = primBuffer.data(pidx_g_6_0_m1 + 28 * idx + 27); 

                // set up pointers to integrals

                auto tg_xxyyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 54); 

                auto tg_xxyyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 55); 

                auto tg_xxyyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 56); 

                auto tg_xxyzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 57); 

                auto tg_xxyzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 58); 

                auto tg_xxyzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 59); 

                auto tg_xxzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 60); 

                auto tg_xxzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 61); 

                auto tg_xxzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 62); 

                auto tg_xyyyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 63); 

                auto tg_xyyyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 64); 

                auto tg_xyyyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 65); 

                auto tg_xyyyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 66); 

                auto tg_xyyyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 67); 

                auto tg_xyyyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 68); 

                auto tg_xyyyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 69); 

                auto tg_xyyyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 70); 

                auto tg_xyyyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 71); 

                auto tg_xyyyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 72); 

                auto tg_xyyyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 73); 

                auto tg_xyyyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 74); 

                auto tg_xyyzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 75); 

                auto tg_xyyzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 76); 

                auto tg_xyyzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 77); 

                auto tg_xyzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 78); 

                auto tg_xyzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 79); 

                auto tg_xyzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 80); 

                auto tg_xzzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 81); 

                auto tg_xzzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 82); 

                auto tg_xzzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 83); 

                auto tg_yyyyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 84); 

                auto tg_yyyyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 85); 

                auto tg_yyyyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 86); 

                auto tg_yyyyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 87); 

                auto tg_yyyyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 88); 

                auto tg_yyyyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 89); 

                auto tg_yyyyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 90); 

                auto tg_yyyyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 91); 

                auto tg_yyyyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 92); 

                auto tg_yyyyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 93); 

                auto tg_yyyyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 94); 

                auto tg_yyyyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 95); 

                auto tg_yyyzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 96); 

                auto tg_yyyzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 97); 

                auto tg_yyyzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 98); 

                auto tg_yyzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 99); 

                auto tg_yyzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 100); 

                auto tg_yyzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 101); 

                auto tg_yzzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 102); 

                auto tg_yzzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 103); 

                auto tg_yzzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 104); 

                auto tg_zzzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 105); 

                auto tg_zzzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 106); 

                auto tg_zzzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 107); 

                // Batch of Integrals (54,108)

                #pragma omp simd aligned(fxn, fza, tg_xxyyzzz_x_0, tg_xxyyzzz_y_0, tg_xxyyzzz_z_0, \
                                         tg_xxyzzzz_x_0, tg_xxyzzzz_y_0, tg_xxyzzzz_z_0, tg_xxzzzzz_x_0, tg_xxzzzzz_y_0, \
                                         tg_xxzzzzz_z_0, tg_xyyyyyy_x_0, tg_xyyyyyy_y_0, tg_xyyyyyy_z_0, tg_xyyyyyz_x_0, \
                                         tg_xyyyyyz_y_0, tg_xyyyyyz_z_0, tg_xyyyyzz_x_0, tg_xyyyyzz_y_0, tg_xyyyyzz_z_0, \
                                         tg_xyyyzzz_x_0, tg_xyyyzzz_y_0, tg_xyyyzzz_z_0, tg_xyyzzz_0_1, tg_xyyzzz_x_0, \
                                         tg_xyyzzz_x_1, tg_xyyzzz_y_0, tg_xyyzzz_y_1, tg_xyyzzz_z_0, tg_xyyzzz_z_1, \
                                         tg_xyyzzzz_x_0, tg_xyyzzzz_y_0, tg_xyyzzzz_z_0, tg_xyzzzz_0_1, tg_xyzzzz_x_0, \
                                         tg_xyzzzz_x_1, tg_xyzzzz_y_0, tg_xyzzzz_y_1, tg_xyzzzz_z_0, tg_xyzzzz_z_1, \
                                         tg_xyzzzzz_x_0, tg_xyzzzzz_y_0, tg_xyzzzzz_z_0, tg_xzzzzz_0_1, tg_xzzzzz_x_0, \
                                         tg_xzzzzz_x_1, tg_xzzzzz_y_0, tg_xzzzzz_y_1, tg_xzzzzz_z_0, tg_xzzzzz_z_1, \
                                         tg_xzzzzzz_x_0, tg_xzzzzzz_y_0, tg_xzzzzzz_z_0, tg_yyyyy_x_0, tg_yyyyy_x_1, \
                                         tg_yyyyy_y_0, tg_yyyyy_y_1, tg_yyyyy_z_0, tg_yyyyy_z_1, tg_yyyyyy_0_1, \
                                         tg_yyyyyy_x_0, tg_yyyyyy_x_1, tg_yyyyyy_y_0, tg_yyyyyy_y_1, tg_yyyyyy_z_0, \
                                         tg_yyyyyy_z_1, tg_yyyyyyy_x_0, tg_yyyyyyy_y_0, tg_yyyyyyy_z_0, tg_yyyyyyz_x_0, \
                                         tg_yyyyyyz_y_0, tg_yyyyyyz_z_0, tg_yyyyyz_0_1, tg_yyyyyz_x_0, tg_yyyyyz_x_1, \
                                         tg_yyyyyz_y_0, tg_yyyyyz_y_1, tg_yyyyyz_z_0, tg_yyyyyz_z_1, tg_yyyyyzz_x_0, \
                                         tg_yyyyyzz_y_0, tg_yyyyyzz_z_0, tg_yyyyz_x_0, tg_yyyyz_x_1, tg_yyyyz_y_0, \
                                         tg_yyyyz_y_1, tg_yyyyz_z_0, tg_yyyyz_z_1, tg_yyyyzz_0_1, tg_yyyyzz_x_0, \
                                         tg_yyyyzz_x_1, tg_yyyyzz_y_0, tg_yyyyzz_y_1, tg_yyyyzz_z_0, tg_yyyyzz_z_1, \
                                         tg_yyyyzzz_x_0, tg_yyyyzzz_y_0, tg_yyyyzzz_z_0, tg_yyyzz_x_0, tg_yyyzz_x_1, \
                                         tg_yyyzz_y_0, tg_yyyzz_y_1, tg_yyyzz_z_0, tg_yyyzz_z_1, tg_yyyzzz_0_1, \
                                         tg_yyyzzz_x_0, tg_yyyzzz_x_1, tg_yyyzzz_y_0, tg_yyyzzz_y_1, tg_yyyzzz_z_0, \
                                         tg_yyyzzz_z_1, tg_yyyzzzz_x_0, tg_yyyzzzz_y_0, tg_yyyzzzz_z_0, tg_yyzzz_x_0, \
                                         tg_yyzzz_x_1, tg_yyzzz_y_0, tg_yyzzz_y_1, tg_yyzzz_z_0, tg_yyzzz_z_1, \
                                         tg_yyzzzz_0_1, tg_yyzzzz_x_0, tg_yyzzzz_x_1, tg_yyzzzz_y_0, tg_yyzzzz_y_1, \
                                         tg_yyzzzz_z_0, tg_yyzzzz_z_1, tg_yyzzzzz_x_0, tg_yyzzzzz_y_0, tg_yyzzzzz_z_0, \
                                         tg_yzzzz_x_0, tg_yzzzz_x_1, tg_yzzzz_y_0, tg_yzzzz_y_1, tg_yzzzz_z_0, tg_yzzzz_z_1, \
                                         tg_yzzzzz_0_1, tg_yzzzzz_x_0, tg_yzzzzz_x_1, tg_yzzzzz_y_0, tg_yzzzzz_y_1, \
                                         tg_yzzzzz_z_0, tg_yzzzzz_z_1, tg_yzzzzzz_x_0, tg_yzzzzzz_y_0, tg_yzzzzzz_z_0, \
                                         tg_zzzzz_x_0, tg_zzzzz_x_1, tg_zzzzz_y_0, tg_zzzzz_y_1, tg_zzzzz_z_0, tg_zzzzz_z_1, \
                                         tg_zzzzzz_0_1, tg_zzzzzz_x_0, tg_zzzzzz_x_1, tg_zzzzzz_y_0, tg_zzzzzz_y_1, \
                                         tg_zzzzzz_z_0, tg_zzzzzz_z_1, tg_zzzzzzz_x_0, tg_zzzzzzz_y_0, tg_zzzzzzz_z_0, wp_x, \
                                         wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyyzzz_x_0[j] = pb_x * tg_xyyzzz_x_0[j] + wp_x[j] * tg_xyyzzz_x_1[j] + 0.5 * fl1_fx * tg_yyzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzz_x_1[j] + 0.5 * fl1_fxn * tg_xyyzzz_0_1[j];

                    tg_xxyyzzz_y_0[j] = pb_x * tg_xyyzzz_y_0[j] + wp_x[j] * tg_xyyzzz_y_1[j] + 0.5 * fl1_fx * tg_yyzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzz_y_1[j];

                    tg_xxyyzzz_z_0[j] = pb_x * tg_xyyzzz_z_0[j] + wp_x[j] * tg_xyyzzz_z_1[j] + 0.5 * fl1_fx * tg_yyzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzz_z_1[j];

                    tg_xxyzzzz_x_0[j] = pb_x * tg_xyzzzz_x_0[j] + wp_x[j] * tg_xyzzzz_x_1[j] + 0.5 * fl1_fx * tg_yzzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xyzzzz_0_1[j];

                    tg_xxyzzzz_y_0[j] = pb_x * tg_xyzzzz_y_0[j] + wp_x[j] * tg_xyzzzz_y_1[j] + 0.5 * fl1_fx * tg_yzzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzz_y_1[j];

                    tg_xxyzzzz_z_0[j] = pb_x * tg_xyzzzz_z_0[j] + wp_x[j] * tg_xyzzzz_z_1[j] + 0.5 * fl1_fx * tg_yzzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzz_z_1[j];

                    tg_xxzzzzz_x_0[j] = pb_x * tg_xzzzzz_x_0[j] + wp_x[j] * tg_xzzzzz_x_1[j] + 0.5 * fl1_fx * tg_zzzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xzzzzz_0_1[j];

                    tg_xxzzzzz_y_0[j] = pb_x * tg_xzzzzz_y_0[j] + wp_x[j] * tg_xzzzzz_y_1[j] + 0.5 * fl1_fx * tg_zzzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzz_y_1[j];

                    tg_xxzzzzz_z_0[j] = pb_x * tg_xzzzzz_z_0[j] + wp_x[j] * tg_xzzzzz_z_1[j] + 0.5 * fl1_fx * tg_zzzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzz_z_1[j];

                    tg_xyyyyyy_x_0[j] = pb_x * tg_yyyyyy_x_0[j] + wp_x[j] * tg_yyyyyy_x_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_0_1[j];

                    tg_xyyyyyy_y_0[j] = pb_x * tg_yyyyyy_y_0[j] + wp_x[j] * tg_yyyyyy_y_1[j];

                    tg_xyyyyyy_z_0[j] = pb_x * tg_yyyyyy_z_0[j] + wp_x[j] * tg_yyyyyy_z_1[j];

                    tg_xyyyyyz_x_0[j] = pb_x * tg_yyyyyz_x_0[j] + wp_x[j] * tg_yyyyyz_x_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_0_1[j];

                    tg_xyyyyyz_y_0[j] = pb_x * tg_yyyyyz_y_0[j] + wp_x[j] * tg_yyyyyz_y_1[j];

                    tg_xyyyyyz_z_0[j] = pb_x * tg_yyyyyz_z_0[j] + wp_x[j] * tg_yyyyyz_z_1[j];

                    tg_xyyyyzz_x_0[j] = pb_x * tg_yyyyzz_x_0[j] + wp_x[j] * tg_yyyyzz_x_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_0_1[j];

                    tg_xyyyyzz_y_0[j] = pb_x * tg_yyyyzz_y_0[j] + wp_x[j] * tg_yyyyzz_y_1[j];

                    tg_xyyyyzz_z_0[j] = pb_x * tg_yyyyzz_z_0[j] + wp_x[j] * tg_yyyyzz_z_1[j];

                    tg_xyyyzzz_x_0[j] = pb_x * tg_yyyzzz_x_0[j] + wp_x[j] * tg_yyyzzz_x_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_0_1[j];

                    tg_xyyyzzz_y_0[j] = pb_x * tg_yyyzzz_y_0[j] + wp_x[j] * tg_yyyzzz_y_1[j];

                    tg_xyyyzzz_z_0[j] = pb_x * tg_yyyzzz_z_0[j] + wp_x[j] * tg_yyyzzz_z_1[j];

                    tg_xyyzzzz_x_0[j] = pb_x * tg_yyzzzz_x_0[j] + wp_x[j] * tg_yyzzzz_x_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_0_1[j];

                    tg_xyyzzzz_y_0[j] = pb_x * tg_yyzzzz_y_0[j] + wp_x[j] * tg_yyzzzz_y_1[j];

                    tg_xyyzzzz_z_0[j] = pb_x * tg_yyzzzz_z_0[j] + wp_x[j] * tg_yyzzzz_z_1[j];

                    tg_xyzzzzz_x_0[j] = pb_x * tg_yzzzzz_x_0[j] + wp_x[j] * tg_yzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_0_1[j];

                    tg_xyzzzzz_y_0[j] = pb_x * tg_yzzzzz_y_0[j] + wp_x[j] * tg_yzzzzz_y_1[j];

                    tg_xyzzzzz_z_0[j] = pb_x * tg_yzzzzz_z_0[j] + wp_x[j] * tg_yzzzzz_z_1[j];

                    tg_xzzzzzz_x_0[j] = pb_x * tg_zzzzzz_x_0[j] + wp_x[j] * tg_zzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_0_1[j];

                    tg_xzzzzzz_y_0[j] = pb_x * tg_zzzzzz_y_0[j] + wp_x[j] * tg_zzzzzz_y_1[j];

                    tg_xzzzzzz_z_0[j] = pb_x * tg_zzzzzz_z_0[j] + wp_x[j] * tg_zzzzzz_z_1[j];

                    tg_yyyyyyy_x_0[j] = pb_y * tg_yyyyyy_x_0[j] + wp_y[j] * tg_yyyyyy_x_1[j] + 3.0 * fl1_fx * tg_yyyyy_x_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyy_x_1[j];

                    tg_yyyyyyy_y_0[j] = pb_y * tg_yyyyyy_y_0[j] + wp_y[j] * tg_yyyyyy_y_1[j] + 3.0 * fl1_fx * tg_yyyyy_y_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyy_y_1[j] + 0.5 * fl1_fxn * tg_yyyyyy_0_1[j];

                    tg_yyyyyyy_z_0[j] = pb_y * tg_yyyyyy_z_0[j] + wp_y[j] * tg_yyyyyy_z_1[j] + 3.0 * fl1_fx * tg_yyyyy_z_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyy_z_1[j];

                    tg_yyyyyyz_x_0[j] = pb_y * tg_yyyyyz_x_0[j] + wp_y[j] * tg_yyyyyz_x_1[j] + 2.5 * fl1_fx * tg_yyyyz_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyz_x_1[j];

                    tg_yyyyyyz_y_0[j] = pb_y * tg_yyyyyz_y_0[j] + wp_y[j] * tg_yyyyyz_y_1[j] + 2.5 * fl1_fx * tg_yyyyz_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyz_y_1[j] + 0.5 * fl1_fxn * tg_yyyyyz_0_1[j];

                    tg_yyyyyyz_z_0[j] = pb_y * tg_yyyyyz_z_0[j] + wp_y[j] * tg_yyyyyz_z_1[j] + 2.5 * fl1_fx * tg_yyyyz_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyz_z_1[j];

                    tg_yyyyyzz_x_0[j] = pb_y * tg_yyyyzz_x_0[j] + wp_y[j] * tg_yyyyzz_x_1[j] + 2.0 * fl1_fx * tg_yyyzz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzz_x_1[j];

                    tg_yyyyyzz_y_0[j] = pb_y * tg_yyyyzz_y_0[j] + wp_y[j] * tg_yyyyzz_y_1[j] + 2.0 * fl1_fx * tg_yyyzz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzz_y_1[j] + 0.5 * fl1_fxn * tg_yyyyzz_0_1[j];

                    tg_yyyyyzz_z_0[j] = pb_y * tg_yyyyzz_z_0[j] + wp_y[j] * tg_yyyyzz_z_1[j] + 2.0 * fl1_fx * tg_yyyzz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzz_z_1[j];

                    tg_yyyyzzz_x_0[j] = pb_y * tg_yyyzzz_x_0[j] + wp_y[j] * tg_yyyzzz_x_1[j] + 1.5 * fl1_fx * tg_yyzzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzz_x_1[j];

                    tg_yyyyzzz_y_0[j] = pb_y * tg_yyyzzz_y_0[j] + wp_y[j] * tg_yyyzzz_y_1[j] + 1.5 * fl1_fx * tg_yyzzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzz_y_1[j] + 0.5 * fl1_fxn * tg_yyyzzz_0_1[j];

                    tg_yyyyzzz_z_0[j] = pb_y * tg_yyyzzz_z_0[j] + wp_y[j] * tg_yyyzzz_z_1[j] + 1.5 * fl1_fx * tg_yyzzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzz_z_1[j];

                    tg_yyyzzzz_x_0[j] = pb_y * tg_yyzzzz_x_0[j] + wp_y[j] * tg_yyzzzz_x_1[j] + fl1_fx * tg_yzzzz_x_0[j] - fl1_fx * fl1_fza * tg_yzzzz_x_1[j];

                    tg_yyyzzzz_y_0[j] = pb_y * tg_yyzzzz_y_0[j] + wp_y[j] * tg_yyzzzz_y_1[j] + fl1_fx * tg_yzzzz_y_0[j] - fl1_fx * fl1_fza * tg_yzzzz_y_1[j] + 0.5 * fl1_fxn * tg_yyzzzz_0_1[j];

                    tg_yyyzzzz_z_0[j] = pb_y * tg_yyzzzz_z_0[j] + wp_y[j] * tg_yyzzzz_z_1[j] + fl1_fx * tg_yzzzz_z_0[j] - fl1_fx * fl1_fza * tg_yzzzz_z_1[j];

                    tg_yyzzzzz_x_0[j] = pb_y * tg_yzzzzz_x_0[j] + wp_y[j] * tg_yzzzzz_x_1[j] + 0.5 * fl1_fx * tg_zzzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzz_x_1[j];

                    tg_yyzzzzz_y_0[j] = pb_y * tg_yzzzzz_y_0[j] + wp_y[j] * tg_yzzzzz_y_1[j] + 0.5 * fl1_fx * tg_zzzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzz_y_1[j] + 0.5 * fl1_fxn * tg_yzzzzz_0_1[j];

                    tg_yyzzzzz_z_0[j] = pb_y * tg_yzzzzz_z_0[j] + wp_y[j] * tg_yzzzzz_z_1[j] + 0.5 * fl1_fx * tg_zzzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzz_z_1[j];

                    tg_yzzzzzz_x_0[j] = pb_y * tg_zzzzzz_x_0[j] + wp_y[j] * tg_zzzzzz_x_1[j];

                    tg_yzzzzzz_y_0[j] = pb_y * tg_zzzzzz_y_0[j] + wp_y[j] * tg_zzzzzz_y_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_0_1[j];

                    tg_yzzzzzz_z_0[j] = pb_y * tg_zzzzzz_z_0[j] + wp_y[j] * tg_zzzzzz_z_1[j];

                    tg_zzzzzzz_x_0[j] = pb_z * tg_zzzzzz_x_0[j] + wp_z[j] * tg_zzzzzz_x_1[j] + 3.0 * fl1_fx * tg_zzzzz_x_0[j] - 3.0 * fl1_fx * fl1_fza * tg_zzzzz_x_1[j];

                    tg_zzzzzzz_y_0[j] = pb_z * tg_zzzzzz_y_0[j] + wp_z[j] * tg_zzzzzz_y_1[j] + 3.0 * fl1_fx * tg_zzzzz_y_0[j] - 3.0 * fl1_fx * fl1_fza * tg_zzzzz_y_1[j];

                    tg_zzzzzzz_z_0[j] = pb_z * tg_zzzzzz_z_0[j] + wp_z[j] * tg_zzzzzz_z_1[j] + 3.0 * fl1_fx * tg_zzzzz_z_0[j] - 3.0 * fl1_fx * fl1_fza * tg_zzzzz_z_1[j] + 0.5 * fl1_fxn * tg_zzzzzz_0_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSL(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSPSL_0_68(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSPSL_68_135(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForSPSL_0_68(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,68)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_x = braGtoPairsBlock.getDistancesPBX();

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fxn = osFactors.data(4 * idx);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tg_0_xxxxxxxx_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx); 

                auto tg_0_xxxxxxxy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 1); 

                auto tg_0_xxxxxxxz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 2); 

                auto tg_0_xxxxxxyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 3); 

                auto tg_0_xxxxxxyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 4); 

                auto tg_0_xxxxxxzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 5); 

                auto tg_0_xxxxxyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 6); 

                auto tg_0_xxxxxyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 7); 

                auto tg_0_xxxxxyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 8); 

                auto tg_0_xxxxxzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 9); 

                auto tg_0_xxxxyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 10); 

                auto tg_0_xxxxyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 11); 

                auto tg_0_xxxxyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 12); 

                auto tg_0_xxxxyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 13); 

                auto tg_0_xxxxzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 14); 

                auto tg_0_xxxyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 15); 

                auto tg_0_xxxyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 16); 

                auto tg_0_xxxyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 17); 

                auto tg_0_xxxyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 18); 

                auto tg_0_xxxyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 19); 

                auto tg_0_xxxzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 20); 

                auto tg_0_xxyyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 21); 

                auto tg_0_xxyyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 22); 

                auto tg_0_xxyyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 23); 

                auto tg_0_xxyyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 24); 

                auto tg_0_xxyyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 25); 

                auto tg_0_xxyzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 26); 

                auto tg_0_xxzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 27); 

                auto tg_0_xyyyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 28); 

                auto tg_0_xyyyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 29); 

                auto tg_0_xyyyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 30); 

                auto tg_0_xyyyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 31); 

                auto tg_0_xyyyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 32); 

                auto tg_0_xyyzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 33); 

                auto tg_0_xyzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 34); 

                auto tg_0_xzzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 35); 

                auto tg_0_yyyyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 36); 

                auto tg_0_yyyyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 37); 

                auto tg_0_yyyyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 38); 

                auto tg_0_yyyyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 39); 

                auto tg_0_yyyyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 40); 

                auto tg_0_yyyzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 41); 

                auto tg_0_yyzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 42); 

                auto tg_0_yzzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 43); 

                auto tg_0_zzzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 44); 

                auto tg_0_xxxxxxxx_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx); 

                auto tg_0_xxxxxxxy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 1); 

                auto tg_0_xxxxxxxz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 2); 

                auto tg_0_xxxxxxyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 3); 

                auto tg_0_xxxxxxyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 4); 

                auto tg_0_xxxxxxzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 5); 

                auto tg_0_xxxxxyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 6); 

                auto tg_0_xxxxxyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 7); 

                auto tg_0_xxxxxyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 8); 

                auto tg_0_xxxxxzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 9); 

                auto tg_0_xxxxyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 10); 

                auto tg_0_xxxxyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 11); 

                auto tg_0_xxxxyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 12); 

                auto tg_0_xxxxyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 13); 

                auto tg_0_xxxxzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 14); 

                auto tg_0_xxxyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 15); 

                auto tg_0_xxxyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 16); 

                auto tg_0_xxxyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 17); 

                auto tg_0_xxxyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 18); 

                auto tg_0_xxxyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 19); 

                auto tg_0_xxxzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 20); 

                auto tg_0_xxyyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 21); 

                auto tg_0_xxyyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 22); 

                auto tg_0_xxyyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 23); 

                auto tg_0_xxyyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 24); 

                auto tg_0_xxyyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 25); 

                auto tg_0_xxyzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 26); 

                auto tg_0_xxzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 27); 

                auto tg_0_xyyyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 28); 

                auto tg_0_xyyyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 29); 

                auto tg_0_xyyyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 30); 

                auto tg_0_xyyyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 31); 

                auto tg_0_xyyyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 32); 

                auto tg_0_xyyzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 33); 

                auto tg_0_xyzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 34); 

                auto tg_0_xzzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 35); 

                auto tg_0_yyyyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 36); 

                auto tg_0_yyyyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 37); 

                auto tg_0_yyyyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 38); 

                auto tg_0_yyyyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 39); 

                auto tg_0_yyyyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 40); 

                auto tg_0_yyyzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 41); 

                auto tg_0_yyzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 42); 

                auto tg_0_yzzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 43); 

                auto tg_0_zzzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 44); 

                auto tg_0_xxxxxxx_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx); 

                auto tg_0_xxxxxxy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 1); 

                auto tg_0_xxxxxxz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 2); 

                auto tg_0_xxxxxyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 3); 

                auto tg_0_xxxxxyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 4); 

                auto tg_0_xxxxxzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 5); 

                auto tg_0_xxxxyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 6); 

                auto tg_0_xxxxyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 7); 

                auto tg_0_xxxxyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 8); 

                auto tg_0_xxxxzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 9); 

                auto tg_0_xxxyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 10); 

                auto tg_0_xxxyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 11); 

                auto tg_0_xxxyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 12); 

                auto tg_0_xxxyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 13); 

                auto tg_0_xxxzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 14); 

                auto tg_0_xxyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 15); 

                auto tg_0_xxyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 16); 

                auto tg_0_xxyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 17); 

                auto tg_0_xxyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 18); 

                auto tg_0_xxyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 19); 

                auto tg_0_xxzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 20); 

                auto tg_0_xyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 21); 

                auto tg_0_xyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 22); 

                auto tg_0_xyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 23); 

                auto tg_0_xyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 24); 

                auto tg_0_xyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 25); 

                auto tg_0_xyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 26); 

                auto tg_0_xzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 27); 

                auto tg_0_yyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 28); 

                auto tg_0_yyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 29); 

                auto tg_0_yyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 30); 

                auto tg_0_yyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 31); 

                auto tg_0_yyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 32); 

                auto tg_0_yyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 33); 

                auto tg_0_yzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 34); 

                auto tg_0_zzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 35); 

                // set up pointers to integrals

                auto tg_x_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx); 

                auto tg_x_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 1); 

                auto tg_x_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 2); 

                auto tg_x_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 3); 

                auto tg_x_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 4); 

                auto tg_x_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 5); 

                auto tg_x_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 6); 

                auto tg_x_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 7); 

                auto tg_x_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 8); 

                auto tg_x_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 9); 

                auto tg_x_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 10); 

                auto tg_x_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 11); 

                auto tg_x_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 12); 

                auto tg_x_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 13); 

                auto tg_x_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 14); 

                auto tg_x_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 15); 

                auto tg_x_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 16); 

                auto tg_x_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 17); 

                auto tg_x_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 18); 

                auto tg_x_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 19); 

                auto tg_x_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 20); 

                auto tg_x_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 21); 

                auto tg_x_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 22); 

                auto tg_x_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 23); 

                auto tg_x_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 24); 

                auto tg_x_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 25); 

                auto tg_x_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 26); 

                auto tg_x_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 27); 

                auto tg_x_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 28); 

                auto tg_x_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 29); 

                auto tg_x_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 30); 

                auto tg_x_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 31); 

                auto tg_x_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 32); 

                auto tg_x_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 33); 

                auto tg_x_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 34); 

                auto tg_x_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 35); 

                auto tg_x_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 36); 

                auto tg_x_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 37); 

                auto tg_x_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 38); 

                auto tg_x_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 39); 

                auto tg_x_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 40); 

                auto tg_x_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 41); 

                auto tg_x_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 42); 

                auto tg_x_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 43); 

                auto tg_x_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 44); 

                auto tg_y_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 45); 

                auto tg_y_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 46); 

                auto tg_y_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 47); 

                auto tg_y_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 48); 

                auto tg_y_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 49); 

                auto tg_y_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 50); 

                auto tg_y_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 51); 

                auto tg_y_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 52); 

                auto tg_y_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 53); 

                auto tg_y_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 54); 

                auto tg_y_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 55); 

                auto tg_y_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 56); 

                auto tg_y_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 57); 

                auto tg_y_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 58); 

                auto tg_y_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 59); 

                auto tg_y_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 60); 

                auto tg_y_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 61); 

                auto tg_y_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 62); 

                auto tg_y_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 63); 

                auto tg_y_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 64); 

                auto tg_y_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 65); 

                auto tg_y_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 66); 

                auto tg_y_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 67); 

                // Batch of Integrals (0,68)

                #pragma omp simd aligned(fxn, tg_0_xxxxxxx_1, tg_0_xxxxxxxx_0, tg_0_xxxxxxxx_1, \
                                         tg_0_xxxxxxxy_0, tg_0_xxxxxxxy_1, tg_0_xxxxxxxz_0, tg_0_xxxxxxxz_1, tg_0_xxxxxxy_1, \
                                         tg_0_xxxxxxyy_0, tg_0_xxxxxxyy_1, tg_0_xxxxxxyz_0, tg_0_xxxxxxyz_1, tg_0_xxxxxxz_1, \
                                         tg_0_xxxxxxzz_0, tg_0_xxxxxxzz_1, tg_0_xxxxxyy_1, tg_0_xxxxxyyy_0, tg_0_xxxxxyyy_1, \
                                         tg_0_xxxxxyyz_0, tg_0_xxxxxyyz_1, tg_0_xxxxxyz_1, tg_0_xxxxxyzz_0, tg_0_xxxxxyzz_1, \
                                         tg_0_xxxxxzz_1, tg_0_xxxxxzzz_0, tg_0_xxxxxzzz_1, tg_0_xxxxyyy_1, tg_0_xxxxyyyy_0, \
                                         tg_0_xxxxyyyy_1, tg_0_xxxxyyyz_0, tg_0_xxxxyyyz_1, tg_0_xxxxyyz_1, tg_0_xxxxyyzz_0, \
                                         tg_0_xxxxyyzz_1, tg_0_xxxxyzz_1, tg_0_xxxxyzzz_0, tg_0_xxxxyzzz_1, tg_0_xxxxzzz_1, \
                                         tg_0_xxxxzzzz_0, tg_0_xxxxzzzz_1, tg_0_xxxyyyy_1, tg_0_xxxyyyyy_0, tg_0_xxxyyyyy_1, \
                                         tg_0_xxxyyyyz_0, tg_0_xxxyyyyz_1, tg_0_xxxyyyz_1, tg_0_xxxyyyzz_0, tg_0_xxxyyyzz_1, \
                                         tg_0_xxxyyzz_1, tg_0_xxxyyzzz_0, tg_0_xxxyyzzz_1, tg_0_xxxyzzz_1, tg_0_xxxyzzzz_0, \
                                         tg_0_xxxyzzzz_1, tg_0_xxxzzzz_1, tg_0_xxxzzzzz_0, tg_0_xxxzzzzz_1, tg_0_xxyyyyy_1, \
                                         tg_0_xxyyyyyy_0, tg_0_xxyyyyyy_1, tg_0_xxyyyyyz_0, tg_0_xxyyyyyz_1, tg_0_xxyyyyz_1, \
                                         tg_0_xxyyyyzz_0, tg_0_xxyyyyzz_1, tg_0_xxyyyzz_1, tg_0_xxyyyzzz_0, tg_0_xxyyyzzz_1, \
                                         tg_0_xxyyzzz_1, tg_0_xxyyzzzz_0, tg_0_xxyyzzzz_1, tg_0_xxyzzzz_1, tg_0_xxyzzzzz_0, \
                                         tg_0_xxyzzzzz_1, tg_0_xxzzzzz_1, tg_0_xxzzzzzz_0, tg_0_xxzzzzzz_1, tg_0_xyyyyyy_1, \
                                         tg_0_xyyyyyyy_0, tg_0_xyyyyyyy_1, tg_0_xyyyyyyz_0, tg_0_xyyyyyyz_1, tg_0_xyyyyyz_1, \
                                         tg_0_xyyyyyzz_0, tg_0_xyyyyyzz_1, tg_0_xyyyyzz_1, tg_0_xyyyyzzz_0, tg_0_xyyyyzzz_1, \
                                         tg_0_xyyyzzz_1, tg_0_xyyyzzzz_0, tg_0_xyyyzzzz_1, tg_0_xyyzzzz_1, tg_0_xyyzzzzz_0, \
                                         tg_0_xyyzzzzz_1, tg_0_xyzzzzz_1, tg_0_xyzzzzzz_0, tg_0_xyzzzzzz_1, tg_0_xzzzzzz_1, \
                                         tg_0_xzzzzzzz_0, tg_0_xzzzzzzz_1, tg_0_yyyyyyy_1, tg_0_yyyyyyyy_0, tg_0_yyyyyyyy_1, \
                                         tg_0_yyyyyyyz_0, tg_0_yyyyyyyz_1, tg_0_yyyyyyz_1, tg_0_yyyyyyzz_0, tg_0_yyyyyyzz_1, \
                                         tg_0_yyyyyzz_1, tg_0_yyyyyzzz_0, tg_0_yyyyyzzz_1, tg_0_yyyyzzz_1, tg_0_yyyyzzzz_0, \
                                         tg_0_yyyyzzzz_1, tg_0_yyyzzzz_1, tg_0_yyyzzzzz_0, tg_0_yyyzzzzz_1, tg_0_yyzzzzz_1, \
                                         tg_0_yyzzzzzz_0, tg_0_yyzzzzzz_1, tg_0_yzzzzzz_1, tg_0_yzzzzzzz_0, tg_0_yzzzzzzz_1, \
                                         tg_0_zzzzzzz_1, tg_0_zzzzzzzz_0, tg_0_zzzzzzzz_1, tg_x_xxxxxxxx_0, tg_x_xxxxxxxy_0, \
                                         tg_x_xxxxxxxz_0, tg_x_xxxxxxyy_0, tg_x_xxxxxxyz_0, tg_x_xxxxxxzz_0, tg_x_xxxxxyyy_0, \
                                         tg_x_xxxxxyyz_0, tg_x_xxxxxyzz_0, tg_x_xxxxxzzz_0, tg_x_xxxxyyyy_0, tg_x_xxxxyyyz_0, \
                                         tg_x_xxxxyyzz_0, tg_x_xxxxyzzz_0, tg_x_xxxxzzzz_0, tg_x_xxxyyyyy_0, tg_x_xxxyyyyz_0, \
                                         tg_x_xxxyyyzz_0, tg_x_xxxyyzzz_0, tg_x_xxxyzzzz_0, tg_x_xxxzzzzz_0, tg_x_xxyyyyyy_0, \
                                         tg_x_xxyyyyyz_0, tg_x_xxyyyyzz_0, tg_x_xxyyyzzz_0, tg_x_xxyyzzzz_0, tg_x_xxyzzzzz_0, \
                                         tg_x_xxzzzzzz_0, tg_x_xyyyyyyy_0, tg_x_xyyyyyyz_0, tg_x_xyyyyyzz_0, tg_x_xyyyyzzz_0, \
                                         tg_x_xyyyzzzz_0, tg_x_xyyzzzzz_0, tg_x_xyzzzzzz_0, tg_x_xzzzzzzz_0, tg_x_yyyyyyyy_0, \
                                         tg_x_yyyyyyyz_0, tg_x_yyyyyyzz_0, tg_x_yyyyyzzz_0, tg_x_yyyyzzzz_0, tg_x_yyyzzzzz_0, \
                                         tg_x_yyzzzzzz_0, tg_x_yzzzzzzz_0, tg_x_zzzzzzzz_0, tg_y_xxxxxxxx_0, tg_y_xxxxxxxy_0, \
                                         tg_y_xxxxxxxz_0, tg_y_xxxxxxyy_0, tg_y_xxxxxxyz_0, tg_y_xxxxxxzz_0, tg_y_xxxxxyyy_0, \
                                         tg_y_xxxxxyyz_0, tg_y_xxxxxyzz_0, tg_y_xxxxxzzz_0, tg_y_xxxxyyyy_0, tg_y_xxxxyyyz_0, \
                                         tg_y_xxxxyyzz_0, tg_y_xxxxyzzz_0, tg_y_xxxxzzzz_0, tg_y_xxxyyyyy_0, tg_y_xxxyyyyz_0, \
                                         tg_y_xxxyyyzz_0, tg_y_xxxyyzzz_0, tg_y_xxxyzzzz_0, tg_y_xxxzzzzz_0, tg_y_xxyyyyyy_0, \
                                         tg_y_xxyyyyyz_0, wp_x, wp_y: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_x_xxxxxxxx_0[j] = pb_x * tg_0_xxxxxxxx_0[j] + wp_x[j] * tg_0_xxxxxxxx_1[j] + 4.0 * fl1_fxn * tg_0_xxxxxxx_1[j];

                    tg_x_xxxxxxxy_0[j] = pb_x * tg_0_xxxxxxxy_0[j] + wp_x[j] * tg_0_xxxxxxxy_1[j] + 3.5 * fl1_fxn * tg_0_xxxxxxy_1[j];

                    tg_x_xxxxxxxz_0[j] = pb_x * tg_0_xxxxxxxz_0[j] + wp_x[j] * tg_0_xxxxxxxz_1[j] + 3.5 * fl1_fxn * tg_0_xxxxxxz_1[j];

                    tg_x_xxxxxxyy_0[j] = pb_x * tg_0_xxxxxxyy_0[j] + wp_x[j] * tg_0_xxxxxxyy_1[j] + 3.0 * fl1_fxn * tg_0_xxxxxyy_1[j];

                    tg_x_xxxxxxyz_0[j] = pb_x * tg_0_xxxxxxyz_0[j] + wp_x[j] * tg_0_xxxxxxyz_1[j] + 3.0 * fl1_fxn * tg_0_xxxxxyz_1[j];

                    tg_x_xxxxxxzz_0[j] = pb_x * tg_0_xxxxxxzz_0[j] + wp_x[j] * tg_0_xxxxxxzz_1[j] + 3.0 * fl1_fxn * tg_0_xxxxxzz_1[j];

                    tg_x_xxxxxyyy_0[j] = pb_x * tg_0_xxxxxyyy_0[j] + wp_x[j] * tg_0_xxxxxyyy_1[j] + 2.5 * fl1_fxn * tg_0_xxxxyyy_1[j];

                    tg_x_xxxxxyyz_0[j] = pb_x * tg_0_xxxxxyyz_0[j] + wp_x[j] * tg_0_xxxxxyyz_1[j] + 2.5 * fl1_fxn * tg_0_xxxxyyz_1[j];

                    tg_x_xxxxxyzz_0[j] = pb_x * tg_0_xxxxxyzz_0[j] + wp_x[j] * tg_0_xxxxxyzz_1[j] + 2.5 * fl1_fxn * tg_0_xxxxyzz_1[j];

                    tg_x_xxxxxzzz_0[j] = pb_x * tg_0_xxxxxzzz_0[j] + wp_x[j] * tg_0_xxxxxzzz_1[j] + 2.5 * fl1_fxn * tg_0_xxxxzzz_1[j];

                    tg_x_xxxxyyyy_0[j] = pb_x * tg_0_xxxxyyyy_0[j] + wp_x[j] * tg_0_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_0_xxxyyyy_1[j];

                    tg_x_xxxxyyyz_0[j] = pb_x * tg_0_xxxxyyyz_0[j] + wp_x[j] * tg_0_xxxxyyyz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyyyz_1[j];

                    tg_x_xxxxyyzz_0[j] = pb_x * tg_0_xxxxyyzz_0[j] + wp_x[j] * tg_0_xxxxyyzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyyzz_1[j];

                    tg_x_xxxxyzzz_0[j] = pb_x * tg_0_xxxxyzzz_0[j] + wp_x[j] * tg_0_xxxxyzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyzzz_1[j];

                    tg_x_xxxxzzzz_0[j] = pb_x * tg_0_xxxxzzzz_0[j] + wp_x[j] * tg_0_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxzzzz_1[j];

                    tg_x_xxxyyyyy_0[j] = pb_x * tg_0_xxxyyyyy_0[j] + wp_x[j] * tg_0_xxxyyyyy_1[j] + 1.5 * fl1_fxn * tg_0_xxyyyyy_1[j];

                    tg_x_xxxyyyyz_0[j] = pb_x * tg_0_xxxyyyyz_0[j] + wp_x[j] * tg_0_xxxyyyyz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyyyz_1[j];

                    tg_x_xxxyyyzz_0[j] = pb_x * tg_0_xxxyyyzz_0[j] + wp_x[j] * tg_0_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyyzz_1[j];

                    tg_x_xxxyyzzz_0[j] = pb_x * tg_0_xxxyyzzz_0[j] + wp_x[j] * tg_0_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyzzz_1[j];

                    tg_x_xxxyzzzz_0[j] = pb_x * tg_0_xxxyzzzz_0[j] + wp_x[j] * tg_0_xxxyzzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyzzzz_1[j];

                    tg_x_xxxzzzzz_0[j] = pb_x * tg_0_xxxzzzzz_0[j] + wp_x[j] * tg_0_xxxzzzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxzzzzz_1[j];

                    tg_x_xxyyyyyy_0[j] = pb_x * tg_0_xxyyyyyy_0[j] + wp_x[j] * tg_0_xxyyyyyy_1[j] + fl1_fxn * tg_0_xyyyyyy_1[j];

                    tg_x_xxyyyyyz_0[j] = pb_x * tg_0_xxyyyyyz_0[j] + wp_x[j] * tg_0_xxyyyyyz_1[j] + fl1_fxn * tg_0_xyyyyyz_1[j];

                    tg_x_xxyyyyzz_0[j] = pb_x * tg_0_xxyyyyzz_0[j] + wp_x[j] * tg_0_xxyyyyzz_1[j] + fl1_fxn * tg_0_xyyyyzz_1[j];

                    tg_x_xxyyyzzz_0[j] = pb_x * tg_0_xxyyyzzz_0[j] + wp_x[j] * tg_0_xxyyyzzz_1[j] + fl1_fxn * tg_0_xyyyzzz_1[j];

                    tg_x_xxyyzzzz_0[j] = pb_x * tg_0_xxyyzzzz_0[j] + wp_x[j] * tg_0_xxyyzzzz_1[j] + fl1_fxn * tg_0_xyyzzzz_1[j];

                    tg_x_xxyzzzzz_0[j] = pb_x * tg_0_xxyzzzzz_0[j] + wp_x[j] * tg_0_xxyzzzzz_1[j] + fl1_fxn * tg_0_xyzzzzz_1[j];

                    tg_x_xxzzzzzz_0[j] = pb_x * tg_0_xxzzzzzz_0[j] + wp_x[j] * tg_0_xxzzzzzz_1[j] + fl1_fxn * tg_0_xzzzzzz_1[j];

                    tg_x_xyyyyyyy_0[j] = pb_x * tg_0_xyyyyyyy_0[j] + wp_x[j] * tg_0_xyyyyyyy_1[j] + 0.5 * fl1_fxn * tg_0_yyyyyyy_1[j];

                    tg_x_xyyyyyyz_0[j] = pb_x * tg_0_xyyyyyyz_0[j] + wp_x[j] * tg_0_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyyyz_1[j];

                    tg_x_xyyyyyzz_0[j] = pb_x * tg_0_xyyyyyzz_0[j] + wp_x[j] * tg_0_xyyyyyzz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyyzz_1[j];

                    tg_x_xyyyyzzz_0[j] = pb_x * tg_0_xyyyyzzz_0[j] + wp_x[j] * tg_0_xyyyyzzz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyzzz_1[j];

                    tg_x_xyyyzzzz_0[j] = pb_x * tg_0_xyyyzzzz_0[j] + wp_x[j] * tg_0_xyyyzzzz_1[j] + 0.5 * fl1_fxn * tg_0_yyyzzzz_1[j];

                    tg_x_xyyzzzzz_0[j] = pb_x * tg_0_xyyzzzzz_0[j] + wp_x[j] * tg_0_xyyzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_yyzzzzz_1[j];

                    tg_x_xyzzzzzz_0[j] = pb_x * tg_0_xyzzzzzz_0[j] + wp_x[j] * tg_0_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_yzzzzzz_1[j];

                    tg_x_xzzzzzzz_0[j] = pb_x * tg_0_xzzzzzzz_0[j] + wp_x[j] * tg_0_xzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzzzzz_1[j];

                    tg_x_yyyyyyyy_0[j] = pb_x * tg_0_yyyyyyyy_0[j] + wp_x[j] * tg_0_yyyyyyyy_1[j];

                    tg_x_yyyyyyyz_0[j] = pb_x * tg_0_yyyyyyyz_0[j] + wp_x[j] * tg_0_yyyyyyyz_1[j];

                    tg_x_yyyyyyzz_0[j] = pb_x * tg_0_yyyyyyzz_0[j] + wp_x[j] * tg_0_yyyyyyzz_1[j];

                    tg_x_yyyyyzzz_0[j] = pb_x * tg_0_yyyyyzzz_0[j] + wp_x[j] * tg_0_yyyyyzzz_1[j];

                    tg_x_yyyyzzzz_0[j] = pb_x * tg_0_yyyyzzzz_0[j] + wp_x[j] * tg_0_yyyyzzzz_1[j];

                    tg_x_yyyzzzzz_0[j] = pb_x * tg_0_yyyzzzzz_0[j] + wp_x[j] * tg_0_yyyzzzzz_1[j];

                    tg_x_yyzzzzzz_0[j] = pb_x * tg_0_yyzzzzzz_0[j] + wp_x[j] * tg_0_yyzzzzzz_1[j];

                    tg_x_yzzzzzzz_0[j] = pb_x * tg_0_yzzzzzzz_0[j] + wp_x[j] * tg_0_yzzzzzzz_1[j];

                    tg_x_zzzzzzzz_0[j] = pb_x * tg_0_zzzzzzzz_0[j] + wp_x[j] * tg_0_zzzzzzzz_1[j];

                    tg_y_xxxxxxxx_0[j] = pb_y * tg_0_xxxxxxxx_0[j] + wp_y[j] * tg_0_xxxxxxxx_1[j];

                    tg_y_xxxxxxxy_0[j] = pb_y * tg_0_xxxxxxxy_0[j] + wp_y[j] * tg_0_xxxxxxxy_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxxx_1[j];

                    tg_y_xxxxxxxz_0[j] = pb_y * tg_0_xxxxxxxz_0[j] + wp_y[j] * tg_0_xxxxxxxz_1[j];

                    tg_y_xxxxxxyy_0[j] = pb_y * tg_0_xxxxxxyy_0[j] + wp_y[j] * tg_0_xxxxxxyy_1[j] + fl1_fxn * tg_0_xxxxxxy_1[j];

                    tg_y_xxxxxxyz_0[j] = pb_y * tg_0_xxxxxxyz_0[j] + wp_y[j] * tg_0_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxxz_1[j];

                    tg_y_xxxxxxzz_0[j] = pb_y * tg_0_xxxxxxzz_0[j] + wp_y[j] * tg_0_xxxxxxzz_1[j];

                    tg_y_xxxxxyyy_0[j] = pb_y * tg_0_xxxxxyyy_0[j] + wp_y[j] * tg_0_xxxxxyyy_1[j] + 1.5 * fl1_fxn * tg_0_xxxxxyy_1[j];

                    tg_y_xxxxxyyz_0[j] = pb_y * tg_0_xxxxxyyz_0[j] + wp_y[j] * tg_0_xxxxxyyz_1[j] + fl1_fxn * tg_0_xxxxxyz_1[j];

                    tg_y_xxxxxyzz_0[j] = pb_y * tg_0_xxxxxyzz_0[j] + wp_y[j] * tg_0_xxxxxyzz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxzz_1[j];

                    tg_y_xxxxxzzz_0[j] = pb_y * tg_0_xxxxxzzz_0[j] + wp_y[j] * tg_0_xxxxxzzz_1[j];

                    tg_y_xxxxyyyy_0[j] = pb_y * tg_0_xxxxyyyy_0[j] + wp_y[j] * tg_0_xxxxyyyy_1[j] + 2.0 * fl1_fxn * tg_0_xxxxyyy_1[j];

                    tg_y_xxxxyyyz_0[j] = pb_y * tg_0_xxxxyyyz_0[j] + wp_y[j] * tg_0_xxxxyyyz_1[j] + 1.5 * fl1_fxn * tg_0_xxxxyyz_1[j];

                    tg_y_xxxxyyzz_0[j] = pb_y * tg_0_xxxxyyzz_0[j] + wp_y[j] * tg_0_xxxxyyzz_1[j] + fl1_fxn * tg_0_xxxxyzz_1[j];

                    tg_y_xxxxyzzz_0[j] = pb_y * tg_0_xxxxyzzz_0[j] + wp_y[j] * tg_0_xxxxyzzz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxzzz_1[j];

                    tg_y_xxxxzzzz_0[j] = pb_y * tg_0_xxxxzzzz_0[j] + wp_y[j] * tg_0_xxxxzzzz_1[j];

                    tg_y_xxxyyyyy_0[j] = pb_y * tg_0_xxxyyyyy_0[j] + wp_y[j] * tg_0_xxxyyyyy_1[j] + 2.5 * fl1_fxn * tg_0_xxxyyyy_1[j];

                    tg_y_xxxyyyyz_0[j] = pb_y * tg_0_xxxyyyyz_0[j] + wp_y[j] * tg_0_xxxyyyyz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyyyz_1[j];

                    tg_y_xxxyyyzz_0[j] = pb_y * tg_0_xxxyyyzz_0[j] + wp_y[j] * tg_0_xxxyyyzz_1[j] + 1.5 * fl1_fxn * tg_0_xxxyyzz_1[j];

                    tg_y_xxxyyzzz_0[j] = pb_y * tg_0_xxxyyzzz_0[j] + wp_y[j] * tg_0_xxxyyzzz_1[j] + fl1_fxn * tg_0_xxxyzzz_1[j];

                    tg_y_xxxyzzzz_0[j] = pb_y * tg_0_xxxyzzzz_0[j] + wp_y[j] * tg_0_xxxyzzzz_1[j] + 0.5 * fl1_fxn * tg_0_xxxzzzz_1[j];

                    tg_y_xxxzzzzz_0[j] = pb_y * tg_0_xxxzzzzz_0[j] + wp_y[j] * tg_0_xxxzzzzz_1[j];

                    tg_y_xxyyyyyy_0[j] = pb_y * tg_0_xxyyyyyy_0[j] + wp_y[j] * tg_0_xxyyyyyy_1[j] + 3.0 * fl1_fxn * tg_0_xxyyyyy_1[j];

                    tg_y_xxyyyyyz_0[j] = pb_y * tg_0_xxyyyyyz_0[j] + wp_y[j] * tg_0_xxyyyyyz_1[j] + 2.5 * fl1_fxn * tg_0_xxyyyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSL_68_135(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (68,135)

        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(PB) = P - B

        auto r_pb_y = braGtoPairsBlock.getDistancesPBY();

        auto r_pb_z = braGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {1, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_8_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_0_8_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {8, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_0_7_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {0, -1, -1, -1}, {7, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fxn = osFactors.data(4 * idx);

                // set up distances R(PB) = P - B

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xxxxxxxx_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx); 

                auto tg_0_xxxxxxxy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 1); 

                auto tg_0_xxxxxxxz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 2); 

                auto tg_0_xxxxxxyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 3); 

                auto tg_0_xxxxxxyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 4); 

                auto tg_0_xxxxxxzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 5); 

                auto tg_0_xxxxxyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 6); 

                auto tg_0_xxxxxyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 7); 

                auto tg_0_xxxxxyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 8); 

                auto tg_0_xxxxxzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 9); 

                auto tg_0_xxxxyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 10); 

                auto tg_0_xxxxyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 11); 

                auto tg_0_xxxxyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 12); 

                auto tg_0_xxxxyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 13); 

                auto tg_0_xxxxzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 14); 

                auto tg_0_xxxyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 15); 

                auto tg_0_xxxyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 16); 

                auto tg_0_xxxyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 17); 

                auto tg_0_xxxyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 18); 

                auto tg_0_xxxyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 19); 

                auto tg_0_xxxzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 20); 

                auto tg_0_xxyyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 21); 

                auto tg_0_xxyyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 22); 

                auto tg_0_xxyyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 23); 

                auto tg_0_xxyyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 24); 

                auto tg_0_xxyyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 25); 

                auto tg_0_xxyzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 26); 

                auto tg_0_xxzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 27); 

                auto tg_0_xyyyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 28); 

                auto tg_0_xyyyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 29); 

                auto tg_0_xyyyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 30); 

                auto tg_0_xyyyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 31); 

                auto tg_0_xyyyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 32); 

                auto tg_0_xyyzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 33); 

                auto tg_0_xyzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 34); 

                auto tg_0_xzzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 35); 

                auto tg_0_yyyyyyyy_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 36); 

                auto tg_0_yyyyyyyz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 37); 

                auto tg_0_yyyyyyzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 38); 

                auto tg_0_yyyyyzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 39); 

                auto tg_0_yyyyzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 40); 

                auto tg_0_yyyzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 41); 

                auto tg_0_yyzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 42); 

                auto tg_0_yzzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 43); 

                auto tg_0_zzzzzzzz_0 = primBuffer.data(pidx_g_0_8_m0 + 45 * idx + 44); 

                auto tg_0_xxxxxxxx_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx); 

                auto tg_0_xxxxxxxy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 1); 

                auto tg_0_xxxxxxxz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 2); 

                auto tg_0_xxxxxxyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 3); 

                auto tg_0_xxxxxxyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 4); 

                auto tg_0_xxxxxxzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 5); 

                auto tg_0_xxxxxyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 6); 

                auto tg_0_xxxxxyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 7); 

                auto tg_0_xxxxxyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 8); 

                auto tg_0_xxxxxzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 9); 

                auto tg_0_xxxxyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 10); 

                auto tg_0_xxxxyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 11); 

                auto tg_0_xxxxyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 12); 

                auto tg_0_xxxxyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 13); 

                auto tg_0_xxxxzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 14); 

                auto tg_0_xxxyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 15); 

                auto tg_0_xxxyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 16); 

                auto tg_0_xxxyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 17); 

                auto tg_0_xxxyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 18); 

                auto tg_0_xxxyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 19); 

                auto tg_0_xxxzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 20); 

                auto tg_0_xxyyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 21); 

                auto tg_0_xxyyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 22); 

                auto tg_0_xxyyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 23); 

                auto tg_0_xxyyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 24); 

                auto tg_0_xxyyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 25); 

                auto tg_0_xxyzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 26); 

                auto tg_0_xxzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 27); 

                auto tg_0_xyyyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 28); 

                auto tg_0_xyyyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 29); 

                auto tg_0_xyyyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 30); 

                auto tg_0_xyyyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 31); 

                auto tg_0_xyyyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 32); 

                auto tg_0_xyyzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 33); 

                auto tg_0_xyzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 34); 

                auto tg_0_xzzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 35); 

                auto tg_0_yyyyyyyy_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 36); 

                auto tg_0_yyyyyyyz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 37); 

                auto tg_0_yyyyyyzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 38); 

                auto tg_0_yyyyyzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 39); 

                auto tg_0_yyyyzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 40); 

                auto tg_0_yyyzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 41); 

                auto tg_0_yyzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 42); 

                auto tg_0_yzzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 43); 

                auto tg_0_zzzzzzzz_1 = primBuffer.data(pidx_g_0_8_m1 + 45 * idx + 44); 

                auto tg_0_xxxxxxx_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx); 

                auto tg_0_xxxxxxy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 1); 

                auto tg_0_xxxxxxz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 2); 

                auto tg_0_xxxxxyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 3); 

                auto tg_0_xxxxxyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 4); 

                auto tg_0_xxxxxzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 5); 

                auto tg_0_xxxxyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 6); 

                auto tg_0_xxxxyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 7); 

                auto tg_0_xxxxyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 8); 

                auto tg_0_xxxxzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 9); 

                auto tg_0_xxxyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 10); 

                auto tg_0_xxxyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 11); 

                auto tg_0_xxxyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 12); 

                auto tg_0_xxxyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 13); 

                auto tg_0_xxxzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 14); 

                auto tg_0_xxyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 15); 

                auto tg_0_xxyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 16); 

                auto tg_0_xxyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 17); 

                auto tg_0_xxyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 18); 

                auto tg_0_xxyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 19); 

                auto tg_0_xxzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 20); 

                auto tg_0_xyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 21); 

                auto tg_0_xyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 22); 

                auto tg_0_xyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 23); 

                auto tg_0_xyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 24); 

                auto tg_0_xyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 25); 

                auto tg_0_xyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 26); 

                auto tg_0_xzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 27); 

                auto tg_0_yyyyyyy_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 28); 

                auto tg_0_yyyyyyz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 29); 

                auto tg_0_yyyyyzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 30); 

                auto tg_0_yyyyzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 31); 

                auto tg_0_yyyzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 32); 

                auto tg_0_yyzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 33); 

                auto tg_0_yzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 34); 

                auto tg_0_zzzzzzz_1 = primBuffer.data(pidx_g_0_7_m1 + 36 * idx + 35); 

                // set up pointers to integrals

                auto tg_y_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 68); 

                auto tg_y_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 69); 

                auto tg_y_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 70); 

                auto tg_y_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 71); 

                auto tg_y_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 72); 

                auto tg_y_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 73); 

                auto tg_y_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 74); 

                auto tg_y_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 75); 

                auto tg_y_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 76); 

                auto tg_y_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 77); 

                auto tg_y_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 78); 

                auto tg_y_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 79); 

                auto tg_y_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 80); 

                auto tg_y_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 81); 

                auto tg_y_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 82); 

                auto tg_y_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 83); 

                auto tg_y_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 84); 

                auto tg_y_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 85); 

                auto tg_y_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 86); 

                auto tg_y_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 87); 

                auto tg_y_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 88); 

                auto tg_y_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 89); 

                auto tg_z_xxxxxxxx_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 90); 

                auto tg_z_xxxxxxxy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 91); 

                auto tg_z_xxxxxxxz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 92); 

                auto tg_z_xxxxxxyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 93); 

                auto tg_z_xxxxxxyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 94); 

                auto tg_z_xxxxxxzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 95); 

                auto tg_z_xxxxxyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 96); 

                auto tg_z_xxxxxyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 97); 

                auto tg_z_xxxxxyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 98); 

                auto tg_z_xxxxxzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 99); 

                auto tg_z_xxxxyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 100); 

                auto tg_z_xxxxyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 101); 

                auto tg_z_xxxxyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 102); 

                auto tg_z_xxxxyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 103); 

                auto tg_z_xxxxzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 104); 

                auto tg_z_xxxyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 105); 

                auto tg_z_xxxyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 106); 

                auto tg_z_xxxyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 107); 

                auto tg_z_xxxyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 108); 

                auto tg_z_xxxyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 109); 

                auto tg_z_xxxzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 110); 

                auto tg_z_xxyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 111); 

                auto tg_z_xxyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 112); 

                auto tg_z_xxyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 113); 

                auto tg_z_xxyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 114); 

                auto tg_z_xxyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 115); 

                auto tg_z_xxyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 116); 

                auto tg_z_xxzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 117); 

                auto tg_z_xyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 118); 

                auto tg_z_xyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 119); 

                auto tg_z_xyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 120); 

                auto tg_z_xyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 121); 

                auto tg_z_xyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 122); 

                auto tg_z_xyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 123); 

                auto tg_z_xyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 124); 

                auto tg_z_xzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 125); 

                auto tg_z_yyyyyyyy_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 126); 

                auto tg_z_yyyyyyyz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 127); 

                auto tg_z_yyyyyyzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 128); 

                auto tg_z_yyyyyzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 129); 

                auto tg_z_yyyyzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 130); 

                auto tg_z_yyyzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 131); 

                auto tg_z_yyzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 132); 

                auto tg_z_yzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 133); 

                auto tg_z_zzzzzzzz_0 = primBuffer.data(pidx_g_1_8_m0 + 135 * idx + 134); 

                // Batch of Integrals (68,135)

                #pragma omp simd aligned(fxn, tg_0_xxxxxxx_1, tg_0_xxxxxxxx_0, tg_0_xxxxxxxx_1, \
                                         tg_0_xxxxxxxy_0, tg_0_xxxxxxxy_1, tg_0_xxxxxxxz_0, tg_0_xxxxxxxz_1, tg_0_xxxxxxy_1, \
                                         tg_0_xxxxxxyy_0, tg_0_xxxxxxyy_1, tg_0_xxxxxxyz_0, tg_0_xxxxxxyz_1, tg_0_xxxxxxz_1, \
                                         tg_0_xxxxxxzz_0, tg_0_xxxxxxzz_1, tg_0_xxxxxyy_1, tg_0_xxxxxyyy_0, tg_0_xxxxxyyy_1, \
                                         tg_0_xxxxxyyz_0, tg_0_xxxxxyyz_1, tg_0_xxxxxyz_1, tg_0_xxxxxyzz_0, tg_0_xxxxxyzz_1, \
                                         tg_0_xxxxxzz_1, tg_0_xxxxxzzz_0, tg_0_xxxxxzzz_1, tg_0_xxxxyyy_1, tg_0_xxxxyyyy_0, \
                                         tg_0_xxxxyyyy_1, tg_0_xxxxyyyz_0, tg_0_xxxxyyyz_1, tg_0_xxxxyyz_1, tg_0_xxxxyyzz_0, \
                                         tg_0_xxxxyyzz_1, tg_0_xxxxyzz_1, tg_0_xxxxyzzz_0, tg_0_xxxxyzzz_1, tg_0_xxxxzzz_1, \
                                         tg_0_xxxxzzzz_0, tg_0_xxxxzzzz_1, tg_0_xxxyyyy_1, tg_0_xxxyyyyy_0, tg_0_xxxyyyyy_1, \
                                         tg_0_xxxyyyyz_0, tg_0_xxxyyyyz_1, tg_0_xxxyyyz_1, tg_0_xxxyyyzz_0, tg_0_xxxyyyzz_1, \
                                         tg_0_xxxyyzz_1, tg_0_xxxyyzzz_0, tg_0_xxxyyzzz_1, tg_0_xxxyzzz_1, tg_0_xxxyzzzz_0, \
                                         tg_0_xxxyzzzz_1, tg_0_xxxzzzz_1, tg_0_xxxzzzzz_0, tg_0_xxxzzzzz_1, tg_0_xxyyyyy_1, \
                                         tg_0_xxyyyyyy_0, tg_0_xxyyyyyy_1, tg_0_xxyyyyyz_0, tg_0_xxyyyyyz_1, tg_0_xxyyyyz_1, \
                                         tg_0_xxyyyyzz_0, tg_0_xxyyyyzz_1, tg_0_xxyyyzz_1, tg_0_xxyyyzzz_0, tg_0_xxyyyzzz_1, \
                                         tg_0_xxyyzzz_1, tg_0_xxyyzzzz_0, tg_0_xxyyzzzz_1, tg_0_xxyzzzz_1, tg_0_xxyzzzzz_0, \
                                         tg_0_xxyzzzzz_1, tg_0_xxzzzzz_1, tg_0_xxzzzzzz_0, tg_0_xxzzzzzz_1, tg_0_xyyyyyy_1, \
                                         tg_0_xyyyyyyy_0, tg_0_xyyyyyyy_1, tg_0_xyyyyyyz_0, tg_0_xyyyyyyz_1, tg_0_xyyyyyz_1, \
                                         tg_0_xyyyyyzz_0, tg_0_xyyyyyzz_1, tg_0_xyyyyzz_1, tg_0_xyyyyzzz_0, tg_0_xyyyyzzz_1, \
                                         tg_0_xyyyzzz_1, tg_0_xyyyzzzz_0, tg_0_xyyyzzzz_1, tg_0_xyyzzzz_1, tg_0_xyyzzzzz_0, \
                                         tg_0_xyyzzzzz_1, tg_0_xyzzzzz_1, tg_0_xyzzzzzz_0, tg_0_xyzzzzzz_1, tg_0_xzzzzzz_1, \
                                         tg_0_xzzzzzzz_0, tg_0_xzzzzzzz_1, tg_0_yyyyyyy_1, tg_0_yyyyyyyy_0, tg_0_yyyyyyyy_1, \
                                         tg_0_yyyyyyyz_0, tg_0_yyyyyyyz_1, tg_0_yyyyyyz_1, tg_0_yyyyyyzz_0, tg_0_yyyyyyzz_1, \
                                         tg_0_yyyyyzz_1, tg_0_yyyyyzzz_0, tg_0_yyyyyzzz_1, tg_0_yyyyzzz_1, tg_0_yyyyzzzz_0, \
                                         tg_0_yyyyzzzz_1, tg_0_yyyzzzz_1, tg_0_yyyzzzzz_0, tg_0_yyyzzzzz_1, tg_0_yyzzzzz_1, \
                                         tg_0_yyzzzzzz_0, tg_0_yyzzzzzz_1, tg_0_yzzzzzz_1, tg_0_yzzzzzzz_0, tg_0_yzzzzzzz_1, \
                                         tg_0_zzzzzzz_1, tg_0_zzzzzzzz_0, tg_0_zzzzzzzz_1, tg_y_xxyyyyzz_0, tg_y_xxyyyzzz_0, \
                                         tg_y_xxyyzzzz_0, tg_y_xxyzzzzz_0, tg_y_xxzzzzzz_0, tg_y_xyyyyyyy_0, tg_y_xyyyyyyz_0, \
                                         tg_y_xyyyyyzz_0, tg_y_xyyyyzzz_0, tg_y_xyyyzzzz_0, tg_y_xyyzzzzz_0, tg_y_xyzzzzzz_0, \
                                         tg_y_xzzzzzzz_0, tg_y_yyyyyyyy_0, tg_y_yyyyyyyz_0, tg_y_yyyyyyzz_0, tg_y_yyyyyzzz_0, \
                                         tg_y_yyyyzzzz_0, tg_y_yyyzzzzz_0, tg_y_yyzzzzzz_0, tg_y_yzzzzzzz_0, tg_y_zzzzzzzz_0, \
                                         tg_z_xxxxxxxx_0, tg_z_xxxxxxxy_0, tg_z_xxxxxxxz_0, tg_z_xxxxxxyy_0, tg_z_xxxxxxyz_0, \
                                         tg_z_xxxxxxzz_0, tg_z_xxxxxyyy_0, tg_z_xxxxxyyz_0, tg_z_xxxxxyzz_0, tg_z_xxxxxzzz_0, \
                                         tg_z_xxxxyyyy_0, tg_z_xxxxyyyz_0, tg_z_xxxxyyzz_0, tg_z_xxxxyzzz_0, tg_z_xxxxzzzz_0, \
                                         tg_z_xxxyyyyy_0, tg_z_xxxyyyyz_0, tg_z_xxxyyyzz_0, tg_z_xxxyyzzz_0, tg_z_xxxyzzzz_0, \
                                         tg_z_xxxzzzzz_0, tg_z_xxyyyyyy_0, tg_z_xxyyyyyz_0, tg_z_xxyyyyzz_0, tg_z_xxyyyzzz_0, \
                                         tg_z_xxyyzzzz_0, tg_z_xxyzzzzz_0, tg_z_xxzzzzzz_0, tg_z_xyyyyyyy_0, tg_z_xyyyyyyz_0, \
                                         tg_z_xyyyyyzz_0, tg_z_xyyyyzzz_0, tg_z_xyyyzzzz_0, tg_z_xyyzzzzz_0, tg_z_xyzzzzzz_0, \
                                         tg_z_xzzzzzzz_0, tg_z_yyyyyyyy_0, tg_z_yyyyyyyz_0, tg_z_yyyyyyzz_0, tg_z_yyyyyzzz_0, \
                                         tg_z_yyyyzzzz_0, tg_z_yyyzzzzz_0, tg_z_yyzzzzzz_0, tg_z_yzzzzzzz_0, tg_z_zzzzzzzz_0, \
                                         wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fxn = fxn[j];

                    tg_y_xxyyyyzz_0[j] = pb_y * tg_0_xxyyyyzz_0[j] + wp_y[j] * tg_0_xxyyyyzz_1[j] + 2.0 * fl1_fxn * tg_0_xxyyyzz_1[j];

                    tg_y_xxyyyzzz_0[j] = pb_y * tg_0_xxyyyzzz_0[j] + wp_y[j] * tg_0_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyzzz_1[j];

                    tg_y_xxyyzzzz_0[j] = pb_y * tg_0_xxyyzzzz_0[j] + wp_y[j] * tg_0_xxyyzzzz_1[j] + fl1_fxn * tg_0_xxyzzzz_1[j];

                    tg_y_xxyzzzzz_0[j] = pb_y * tg_0_xxyzzzzz_0[j] + wp_y[j] * tg_0_xxyzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_xxzzzzz_1[j];

                    tg_y_xxzzzzzz_0[j] = pb_y * tg_0_xxzzzzzz_0[j] + wp_y[j] * tg_0_xxzzzzzz_1[j];

                    tg_y_xyyyyyyy_0[j] = pb_y * tg_0_xyyyyyyy_0[j] + wp_y[j] * tg_0_xyyyyyyy_1[j] + 3.5 * fl1_fxn * tg_0_xyyyyyy_1[j];

                    tg_y_xyyyyyyz_0[j] = pb_y * tg_0_xyyyyyyz_0[j] + wp_y[j] * tg_0_xyyyyyyz_1[j] + 3.0 * fl1_fxn * tg_0_xyyyyyz_1[j];

                    tg_y_xyyyyyzz_0[j] = pb_y * tg_0_xyyyyyzz_0[j] + wp_y[j] * tg_0_xyyyyyzz_1[j] + 2.5 * fl1_fxn * tg_0_xyyyyzz_1[j];

                    tg_y_xyyyyzzz_0[j] = pb_y * tg_0_xyyyyzzz_0[j] + wp_y[j] * tg_0_xyyyyzzz_1[j] + 2.0 * fl1_fxn * tg_0_xyyyzzz_1[j];

                    tg_y_xyyyzzzz_0[j] = pb_y * tg_0_xyyyzzzz_0[j] + wp_y[j] * tg_0_xyyyzzzz_1[j] + 1.5 * fl1_fxn * tg_0_xyyzzzz_1[j];

                    tg_y_xyyzzzzz_0[j] = pb_y * tg_0_xyyzzzzz_0[j] + wp_y[j] * tg_0_xyyzzzzz_1[j] + fl1_fxn * tg_0_xyzzzzz_1[j];

                    tg_y_xyzzzzzz_0[j] = pb_y * tg_0_xyzzzzzz_0[j] + wp_y[j] * tg_0_xyzzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_xzzzzzz_1[j];

                    tg_y_xzzzzzzz_0[j] = pb_y * tg_0_xzzzzzzz_0[j] + wp_y[j] * tg_0_xzzzzzzz_1[j];

                    tg_y_yyyyyyyy_0[j] = pb_y * tg_0_yyyyyyyy_0[j] + wp_y[j] * tg_0_yyyyyyyy_1[j] + 4.0 * fl1_fxn * tg_0_yyyyyyy_1[j];

                    tg_y_yyyyyyyz_0[j] = pb_y * tg_0_yyyyyyyz_0[j] + wp_y[j] * tg_0_yyyyyyyz_1[j] + 3.5 * fl1_fxn * tg_0_yyyyyyz_1[j];

                    tg_y_yyyyyyzz_0[j] = pb_y * tg_0_yyyyyyzz_0[j] + wp_y[j] * tg_0_yyyyyyzz_1[j] + 3.0 * fl1_fxn * tg_0_yyyyyzz_1[j];

                    tg_y_yyyyyzzz_0[j] = pb_y * tg_0_yyyyyzzz_0[j] + wp_y[j] * tg_0_yyyyyzzz_1[j] + 2.5 * fl1_fxn * tg_0_yyyyzzz_1[j];

                    tg_y_yyyyzzzz_0[j] = pb_y * tg_0_yyyyzzzz_0[j] + wp_y[j] * tg_0_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_yyyzzzz_1[j];

                    tg_y_yyyzzzzz_0[j] = pb_y * tg_0_yyyzzzzz_0[j] + wp_y[j] * tg_0_yyyzzzzz_1[j] + 1.5 * fl1_fxn * tg_0_yyzzzzz_1[j];

                    tg_y_yyzzzzzz_0[j] = pb_y * tg_0_yyzzzzzz_0[j] + wp_y[j] * tg_0_yyzzzzzz_1[j] + fl1_fxn * tg_0_yzzzzzz_1[j];

                    tg_y_yzzzzzzz_0[j] = pb_y * tg_0_yzzzzzzz_0[j] + wp_y[j] * tg_0_yzzzzzzz_1[j] + 0.5 * fl1_fxn * tg_0_zzzzzzz_1[j];

                    tg_y_zzzzzzzz_0[j] = pb_y * tg_0_zzzzzzzz_0[j] + wp_y[j] * tg_0_zzzzzzzz_1[j];

                    tg_z_xxxxxxxx_0[j] = pb_z * tg_0_xxxxxxxx_0[j] + wp_z[j] * tg_0_xxxxxxxx_1[j];

                    tg_z_xxxxxxxy_0[j] = pb_z * tg_0_xxxxxxxy_0[j] + wp_z[j] * tg_0_xxxxxxxy_1[j];

                    tg_z_xxxxxxxz_0[j] = pb_z * tg_0_xxxxxxxz_0[j] + wp_z[j] * tg_0_xxxxxxxz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxxx_1[j];

                    tg_z_xxxxxxyy_0[j] = pb_z * tg_0_xxxxxxyy_0[j] + wp_z[j] * tg_0_xxxxxxyy_1[j];

                    tg_z_xxxxxxyz_0[j] = pb_z * tg_0_xxxxxxyz_0[j] + wp_z[j] * tg_0_xxxxxxyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxxy_1[j];

                    tg_z_xxxxxxzz_0[j] = pb_z * tg_0_xxxxxxzz_0[j] + wp_z[j] * tg_0_xxxxxxzz_1[j] + fl1_fxn * tg_0_xxxxxxz_1[j];

                    tg_z_xxxxxyyy_0[j] = pb_z * tg_0_xxxxxyyy_0[j] + wp_z[j] * tg_0_xxxxxyyy_1[j];

                    tg_z_xxxxxyyz_0[j] = pb_z * tg_0_xxxxxyyz_0[j] + wp_z[j] * tg_0_xxxxxyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxxyy_1[j];

                    tg_z_xxxxxyzz_0[j] = pb_z * tg_0_xxxxxyzz_0[j] + wp_z[j] * tg_0_xxxxxyzz_1[j] + fl1_fxn * tg_0_xxxxxyz_1[j];

                    tg_z_xxxxxzzz_0[j] = pb_z * tg_0_xxxxxzzz_0[j] + wp_z[j] * tg_0_xxxxxzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxxxxzz_1[j];

                    tg_z_xxxxyyyy_0[j] = pb_z * tg_0_xxxxyyyy_0[j] + wp_z[j] * tg_0_xxxxyyyy_1[j];

                    tg_z_xxxxyyyz_0[j] = pb_z * tg_0_xxxxyyyz_0[j] + wp_z[j] * tg_0_xxxxyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxxyyy_1[j];

                    tg_z_xxxxyyzz_0[j] = pb_z * tg_0_xxxxyyzz_0[j] + wp_z[j] * tg_0_xxxxyyzz_1[j] + fl1_fxn * tg_0_xxxxyyz_1[j];

                    tg_z_xxxxyzzz_0[j] = pb_z * tg_0_xxxxyzzz_0[j] + wp_z[j] * tg_0_xxxxyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxxxyzz_1[j];

                    tg_z_xxxxzzzz_0[j] = pb_z * tg_0_xxxxzzzz_0[j] + wp_z[j] * tg_0_xxxxzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxxzzz_1[j];

                    tg_z_xxxyyyyy_0[j] = pb_z * tg_0_xxxyyyyy_0[j] + wp_z[j] * tg_0_xxxyyyyy_1[j];

                    tg_z_xxxyyyyz_0[j] = pb_z * tg_0_xxxyyyyz_0[j] + wp_z[j] * tg_0_xxxyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxxyyyy_1[j];

                    tg_z_xxxyyyzz_0[j] = pb_z * tg_0_xxxyyyzz_0[j] + wp_z[j] * tg_0_xxxyyyzz_1[j] + fl1_fxn * tg_0_xxxyyyz_1[j];

                    tg_z_xxxyyzzz_0[j] = pb_z * tg_0_xxxyyzzz_0[j] + wp_z[j] * tg_0_xxxyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxxyyzz_1[j];

                    tg_z_xxxyzzzz_0[j] = pb_z * tg_0_xxxyzzzz_0[j] + wp_z[j] * tg_0_xxxyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxxyzzz_1[j];

                    tg_z_xxxzzzzz_0[j] = pb_z * tg_0_xxxzzzzz_0[j] + wp_z[j] * tg_0_xxxzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_xxxzzzz_1[j];

                    tg_z_xxyyyyyy_0[j] = pb_z * tg_0_xxyyyyyy_0[j] + wp_z[j] * tg_0_xxyyyyyy_1[j];

                    tg_z_xxyyyyyz_0[j] = pb_z * tg_0_xxyyyyyz_0[j] + wp_z[j] * tg_0_xxyyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xxyyyyy_1[j];

                    tg_z_xxyyyyzz_0[j] = pb_z * tg_0_xxyyyyzz_0[j] + wp_z[j] * tg_0_xxyyyyzz_1[j] + fl1_fxn * tg_0_xxyyyyz_1[j];

                    tg_z_xxyyyzzz_0[j] = pb_z * tg_0_xxyyyzzz_0[j] + wp_z[j] * tg_0_xxyyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xxyyyzz_1[j];

                    tg_z_xxyyzzzz_0[j] = pb_z * tg_0_xxyyzzzz_0[j] + wp_z[j] * tg_0_xxyyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xxyyzzz_1[j];

                    tg_z_xxyzzzzz_0[j] = pb_z * tg_0_xxyzzzzz_0[j] + wp_z[j] * tg_0_xxyzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_xxyzzzz_1[j];

                    tg_z_xxzzzzzz_0[j] = pb_z * tg_0_xxzzzzzz_0[j] + wp_z[j] * tg_0_xxzzzzzz_1[j] + 3.0 * fl1_fxn * tg_0_xxzzzzz_1[j];

                    tg_z_xyyyyyyy_0[j] = pb_z * tg_0_xyyyyyyy_0[j] + wp_z[j] * tg_0_xyyyyyyy_1[j];

                    tg_z_xyyyyyyz_0[j] = pb_z * tg_0_xyyyyyyz_0[j] + wp_z[j] * tg_0_xyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_xyyyyyy_1[j];

                    tg_z_xyyyyyzz_0[j] = pb_z * tg_0_xyyyyyzz_0[j] + wp_z[j] * tg_0_xyyyyyzz_1[j] + fl1_fxn * tg_0_xyyyyyz_1[j];

                    tg_z_xyyyyzzz_0[j] = pb_z * tg_0_xyyyyzzz_0[j] + wp_z[j] * tg_0_xyyyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_xyyyyzz_1[j];

                    tg_z_xyyyzzzz_0[j] = pb_z * tg_0_xyyyzzzz_0[j] + wp_z[j] * tg_0_xyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_xyyyzzz_1[j];

                    tg_z_xyyzzzzz_0[j] = pb_z * tg_0_xyyzzzzz_0[j] + wp_z[j] * tg_0_xyyzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_xyyzzzz_1[j];

                    tg_z_xyzzzzzz_0[j] = pb_z * tg_0_xyzzzzzz_0[j] + wp_z[j] * tg_0_xyzzzzzz_1[j] + 3.0 * fl1_fxn * tg_0_xyzzzzz_1[j];

                    tg_z_xzzzzzzz_0[j] = pb_z * tg_0_xzzzzzzz_0[j] + wp_z[j] * tg_0_xzzzzzzz_1[j] + 3.5 * fl1_fxn * tg_0_xzzzzzz_1[j];

                    tg_z_yyyyyyyy_0[j] = pb_z * tg_0_yyyyyyyy_0[j] + wp_z[j] * tg_0_yyyyyyyy_1[j];

                    tg_z_yyyyyyyz_0[j] = pb_z * tg_0_yyyyyyyz_0[j] + wp_z[j] * tg_0_yyyyyyyz_1[j] + 0.5 * fl1_fxn * tg_0_yyyyyyy_1[j];

                    tg_z_yyyyyyzz_0[j] = pb_z * tg_0_yyyyyyzz_0[j] + wp_z[j] * tg_0_yyyyyyzz_1[j] + fl1_fxn * tg_0_yyyyyyz_1[j];

                    tg_z_yyyyyzzz_0[j] = pb_z * tg_0_yyyyyzzz_0[j] + wp_z[j] * tg_0_yyyyyzzz_1[j] + 1.5 * fl1_fxn * tg_0_yyyyyzz_1[j];

                    tg_z_yyyyzzzz_0[j] = pb_z * tg_0_yyyyzzzz_0[j] + wp_z[j] * tg_0_yyyyzzzz_1[j] + 2.0 * fl1_fxn * tg_0_yyyyzzz_1[j];

                    tg_z_yyyzzzzz_0[j] = pb_z * tg_0_yyyzzzzz_0[j] + wp_z[j] * tg_0_yyyzzzzz_1[j] + 2.5 * fl1_fxn * tg_0_yyyzzzz_1[j];

                    tg_z_yyzzzzzz_0[j] = pb_z * tg_0_yyzzzzzz_0[j] + wp_z[j] * tg_0_yyzzzzzz_1[j] + 3.0 * fl1_fxn * tg_0_yyzzzzz_1[j];

                    tg_z_yzzzzzzz_0[j] = pb_z * tg_0_yzzzzzzz_0[j] + wp_z[j] * tg_0_yzzzzzzz_1[j] + 3.5 * fl1_fxn * tg_0_yzzzzzz_1[j];

                    tg_z_zzzzzzzz_0[j] = pb_z * tg_0_zzzzzzzz_0[j] + wp_z[j] * tg_0_zzzzzzzz_1[j] + 4.0 * fl1_fxn * tg_0_zzzzzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSP(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wpDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        erirecfunc::compElectronRepulsionForSLSP_0_68(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      wpDistances, 
                                                      braGtoPairsBlock,
                                                      ketGtoPairsBlock,
                                                      nKetPrimPairs,
                                                      iContrPair); 

        erirecfunc::compElectronRepulsionForSLSP_68_135(primBuffer,
                                                        recursionMap,
                                                        osFactors,
                                                        wpDistances, 
                                                        braGtoPairsBlock,
                                                        ketGtoPairsBlock,
                                                        nKetPrimPairs,
                                                        iContrPair); 
    }

    void
    compElectronRepulsionForSLSP_0_68(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& wpDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetPrimPairs,
                                      const int32_t              iContrPair)
    {
        // Batch of Integrals (0,68)

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
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xxxxxxx_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx); 

                auto tg_xxxxxxx_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 1); 

                auto tg_xxxxxxx_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 2); 

                auto tg_xxxxxxy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 3); 

                auto tg_xxxxxxy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 4); 

                auto tg_xxxxxxy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 5); 

                auto tg_xxxxxxz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 6); 

                auto tg_xxxxxxz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 7); 

                auto tg_xxxxxxz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 8); 

                auto tg_xxxxxyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 9); 

                auto tg_xxxxxyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 10); 

                auto tg_xxxxxyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 11); 

                auto tg_xxxxxyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 12); 

                auto tg_xxxxxyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 13); 

                auto tg_xxxxxyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 14); 

                auto tg_xxxxxzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 15); 

                auto tg_xxxxxzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 16); 

                auto tg_xxxxxzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 17); 

                auto tg_xxxxyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 18); 

                auto tg_xxxxyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 19); 

                auto tg_xxxxyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 20); 

                auto tg_xxxxyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 21); 

                auto tg_xxxxyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 22); 

                auto tg_xxxxyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 23); 

                auto tg_xxxxyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 24); 

                auto tg_xxxxyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 25); 

                auto tg_xxxxyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 26); 

                auto tg_xxxxzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 27); 

                auto tg_xxxxzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 28); 

                auto tg_xxxxzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 29); 

                auto tg_xxxyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 30); 

                auto tg_xxxyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 31); 

                auto tg_xxxyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 32); 

                auto tg_xxxyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 33); 

                auto tg_xxxyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 34); 

                auto tg_xxxyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 35); 

                auto tg_xxxyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 36); 

                auto tg_xxxyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 37); 

                auto tg_xxxyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 38); 

                auto tg_xxxyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 39); 

                auto tg_xxxyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 40); 

                auto tg_xxxyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 41); 

                auto tg_xxxzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 42); 

                auto tg_xxxzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 43); 

                auto tg_xxxzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 44); 

                auto tg_xxyyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 45); 

                auto tg_xxyyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 46); 

                auto tg_xxyyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 47); 

                auto tg_xxyyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 48); 

                auto tg_xxyyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 49); 

                auto tg_xxyyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 50); 

                auto tg_xxyyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 51); 

                auto tg_xxyyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 52); 

                auto tg_xxyyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 53); 

                auto tg_xxyyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 54); 

                auto tg_xxyyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 55); 

                auto tg_xxyyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 56); 

                auto tg_xxyzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 57); 

                auto tg_xxyzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 58); 

                auto tg_xxyzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 59); 

                auto tg_xxzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 60); 

                auto tg_xxzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 61); 

                auto tg_xxzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 62); 

                auto tg_xyyyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 63); 

                auto tg_xyyyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 64); 

                auto tg_xyyyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 65); 

                auto tg_xyyyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 66); 

                auto tg_xyyyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 67); 

                auto tg_xxxxxxx_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx); 

                auto tg_xxxxxxx_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 1); 

                auto tg_xxxxxxx_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 2); 

                auto tg_xxxxxxy_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 3); 

                auto tg_xxxxxxy_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 4); 

                auto tg_xxxxxxy_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 5); 

                auto tg_xxxxxxz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 6); 

                auto tg_xxxxxxz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 7); 

                auto tg_xxxxxxz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 8); 

                auto tg_xxxxxyy_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 9); 

                auto tg_xxxxxyy_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 10); 

                auto tg_xxxxxyy_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 11); 

                auto tg_xxxxxyz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 12); 

                auto tg_xxxxxyz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 13); 

                auto tg_xxxxxyz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 14); 

                auto tg_xxxxxzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 15); 

                auto tg_xxxxxzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 16); 

                auto tg_xxxxxzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 17); 

                auto tg_xxxxyyy_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 18); 

                auto tg_xxxxyyy_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 19); 

                auto tg_xxxxyyy_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 20); 

                auto tg_xxxxyyz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 21); 

                auto tg_xxxxyyz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 22); 

                auto tg_xxxxyyz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 23); 

                auto tg_xxxxyzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 24); 

                auto tg_xxxxyzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 25); 

                auto tg_xxxxyzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 26); 

                auto tg_xxxxzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 27); 

                auto tg_xxxxzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 28); 

                auto tg_xxxxzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 29); 

                auto tg_xxxyyyy_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 30); 

                auto tg_xxxyyyy_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 31); 

                auto tg_xxxyyyy_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 32); 

                auto tg_xxxyyyz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 33); 

                auto tg_xxxyyyz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 34); 

                auto tg_xxxyyyz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 35); 

                auto tg_xxxyyzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 36); 

                auto tg_xxxyyzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 37); 

                auto tg_xxxyyzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 38); 

                auto tg_xxxyzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 39); 

                auto tg_xxxyzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 40); 

                auto tg_xxxyzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 41); 

                auto tg_xxxzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 42); 

                auto tg_xxxzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 43); 

                auto tg_xxxzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 44); 

                auto tg_xxyyyyy_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 45); 

                auto tg_xxyyyyy_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 46); 

                auto tg_xxyyyyy_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 47); 

                auto tg_xxyyyyz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 48); 

                auto tg_xxyyyyz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 49); 

                auto tg_xxyyyyz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 50); 

                auto tg_xxyyyzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 51); 

                auto tg_xxyyyzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 52); 

                auto tg_xxyyyzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 53); 

                auto tg_xxyyzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 54); 

                auto tg_xxyyzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 55); 

                auto tg_xxyyzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 56); 

                auto tg_xxyzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 57); 

                auto tg_xxyzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 58); 

                auto tg_xxyzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 59); 

                auto tg_xxzzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 60); 

                auto tg_xxzzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 61); 

                auto tg_xxzzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 62); 

                auto tg_xyyyyyy_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 63); 

                auto tg_xyyyyyy_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 64); 

                auto tg_xyyyyyy_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 65); 

                auto tg_xyyyyyz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 66); 

                auto tg_xyyyyyz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 67); 

                auto tg_xxxxxx_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx); 

                auto tg_xxxxxx_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 1); 

                auto tg_xxxxxx_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 2); 

                auto tg_xxxxxy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 3); 

                auto tg_xxxxxy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 4); 

                auto tg_xxxxxy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 5); 

                auto tg_xxxxxz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 6); 

                auto tg_xxxxxz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 7); 

                auto tg_xxxxxz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 8); 

                auto tg_xxxxyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 9); 

                auto tg_xxxxyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 10); 

                auto tg_xxxxyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 11); 

                auto tg_xxxxyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 12); 

                auto tg_xxxxyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 13); 

                auto tg_xxxxyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 14); 

                auto tg_xxxxzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 15); 

                auto tg_xxxxzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 16); 

                auto tg_xxxxzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 17); 

                auto tg_xxxyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 18); 

                auto tg_xxxyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 19); 

                auto tg_xxxyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 20); 

                auto tg_xxxyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 21); 

                auto tg_xxxyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 22); 

                auto tg_xxxyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 23); 

                auto tg_xxxyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 24); 

                auto tg_xxxyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 25); 

                auto tg_xxxyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 26); 

                auto tg_xxxzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 27); 

                auto tg_xxxzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 28); 

                auto tg_xxxzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 29); 

                auto tg_xxyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 30); 

                auto tg_xxyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 31); 

                auto tg_xxyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 32); 

                auto tg_xxyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 33); 

                auto tg_xxyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 34); 

                auto tg_xxyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 35); 

                auto tg_xxyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 36); 

                auto tg_xxyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 37); 

                auto tg_xxyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 38); 

                auto tg_xxyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 39); 

                auto tg_xxyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 40); 

                auto tg_xxyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 41); 

                auto tg_xxzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 42); 

                auto tg_xxzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 43); 

                auto tg_xxzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 44); 

                auto tg_xyyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 45); 

                auto tg_xyyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 46); 

                auto tg_xyyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 47); 

                auto tg_xyyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 48); 

                auto tg_xyyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 49); 

                auto tg_xyyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 50); 

                auto tg_xyyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 51); 

                auto tg_xyyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 52); 

                auto tg_xyyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 53); 

                auto tg_xyyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 54); 

                auto tg_xyyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 55); 

                auto tg_xyyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 56); 

                auto tg_xyzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 57); 

                auto tg_xyzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 58); 

                auto tg_xyzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 59); 

                auto tg_xzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 60); 

                auto tg_xzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 61); 

                auto tg_xzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 62); 

                auto tg_yyyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 63); 

                auto tg_yyyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 64); 

                auto tg_yyyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 65); 

                auto tg_yyyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 66); 

                auto tg_yyyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 67); 

                auto tg_xxxxxx_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx); 

                auto tg_xxxxxx_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 1); 

                auto tg_xxxxxx_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 2); 

                auto tg_xxxxxy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 3); 

                auto tg_xxxxxy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 4); 

                auto tg_xxxxxy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 5); 

                auto tg_xxxxxz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 6); 

                auto tg_xxxxxz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 7); 

                auto tg_xxxxxz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 8); 

                auto tg_xxxxyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 9); 

                auto tg_xxxxyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 10); 

                auto tg_xxxxyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 11); 

                auto tg_xxxxyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 12); 

                auto tg_xxxxyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 13); 

                auto tg_xxxxyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 14); 

                auto tg_xxxxzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 15); 

                auto tg_xxxxzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 16); 

                auto tg_xxxxzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 17); 

                auto tg_xxxyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 18); 

                auto tg_xxxyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 19); 

                auto tg_xxxyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 20); 

                auto tg_xxxyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 21); 

                auto tg_xxxyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 22); 

                auto tg_xxxyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 23); 

                auto tg_xxxyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 24); 

                auto tg_xxxyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 25); 

                auto tg_xxxyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 26); 

                auto tg_xxxzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 27); 

                auto tg_xxxzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 28); 

                auto tg_xxxzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 29); 

                auto tg_xxyyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 30); 

                auto tg_xxyyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 31); 

                auto tg_xxyyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 32); 

                auto tg_xxyyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 33); 

                auto tg_xxyyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 34); 

                auto tg_xxyyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 35); 

                auto tg_xxyyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 36); 

                auto tg_xxyyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 37); 

                auto tg_xxyyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 38); 

                auto tg_xxyzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 39); 

                auto tg_xxyzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 40); 

                auto tg_xxyzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 41); 

                auto tg_xxzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 42); 

                auto tg_xxzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 43); 

                auto tg_xxzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 44); 

                auto tg_xyyyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 45); 

                auto tg_xyyyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 46); 

                auto tg_xyyyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 47); 

                auto tg_xyyyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 48); 

                auto tg_xyyyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 49); 

                auto tg_xyyyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 50); 

                auto tg_xyyyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 51); 

                auto tg_xyyyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 52); 

                auto tg_xyyyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 53); 

                auto tg_xyyzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 54); 

                auto tg_xyyzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 55); 

                auto tg_xyyzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 56); 

                auto tg_xyzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 57); 

                auto tg_xyzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 58); 

                auto tg_xyzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 59); 

                auto tg_xzzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 60); 

                auto tg_xzzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 61); 

                auto tg_xzzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 62); 

                auto tg_yyyyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 63); 

                auto tg_yyyyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 64); 

                auto tg_yyyyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 65); 

                auto tg_yyyyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 66); 

                auto tg_yyyyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 67); 

                auto tg_xxxxxxx_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx); 

                auto tg_xxxxxxy_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 1); 

                auto tg_xxxxxxz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 2); 

                auto tg_xxxxxyy_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 3); 

                auto tg_xxxxxyz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 4); 

                auto tg_xxxxxzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 5); 

                auto tg_xxxxyyy_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 6); 

                auto tg_xxxxyyz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 7); 

                auto tg_xxxxyzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 8); 

                auto tg_xxxxzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 9); 

                auto tg_xxxyyyy_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 10); 

                auto tg_xxxyyyz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 11); 

                auto tg_xxxyyzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 12); 

                auto tg_xxxyzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 13); 

                auto tg_xxxzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 14); 

                auto tg_xxyyyyy_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 15); 

                auto tg_xxyyyyz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 16); 

                auto tg_xxyyyzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 17); 

                auto tg_xxyyzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 18); 

                auto tg_xxyzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 19); 

                auto tg_xxzzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 20); 

                auto tg_xyyyyyy_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 21); 

                auto tg_xyyyyyz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 22); 

                // set up pointers to integrals

                auto tg_xxxxxxxx_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx); 

                auto tg_xxxxxxxx_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 1); 

                auto tg_xxxxxxxx_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 2); 

                auto tg_xxxxxxxy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 3); 

                auto tg_xxxxxxxy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 4); 

                auto tg_xxxxxxxy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 5); 

                auto tg_xxxxxxxz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 6); 

                auto tg_xxxxxxxz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 7); 

                auto tg_xxxxxxxz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 8); 

                auto tg_xxxxxxyy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 9); 

                auto tg_xxxxxxyy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 10); 

                auto tg_xxxxxxyy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 11); 

                auto tg_xxxxxxyz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 12); 

                auto tg_xxxxxxyz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 13); 

                auto tg_xxxxxxyz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 14); 

                auto tg_xxxxxxzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 15); 

                auto tg_xxxxxxzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 16); 

                auto tg_xxxxxxzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 17); 

                auto tg_xxxxxyyy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 18); 

                auto tg_xxxxxyyy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 19); 

                auto tg_xxxxxyyy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 20); 

                auto tg_xxxxxyyz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 21); 

                auto tg_xxxxxyyz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 22); 

                auto tg_xxxxxyyz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 23); 

                auto tg_xxxxxyzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 24); 

                auto tg_xxxxxyzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 25); 

                auto tg_xxxxxyzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 26); 

                auto tg_xxxxxzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 27); 

                auto tg_xxxxxzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 28); 

                auto tg_xxxxxzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 29); 

                auto tg_xxxxyyyy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 30); 

                auto tg_xxxxyyyy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 31); 

                auto tg_xxxxyyyy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 32); 

                auto tg_xxxxyyyz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 33); 

                auto tg_xxxxyyyz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 34); 

                auto tg_xxxxyyyz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 35); 

                auto tg_xxxxyyzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 36); 

                auto tg_xxxxyyzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 37); 

                auto tg_xxxxyyzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 38); 

                auto tg_xxxxyzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 39); 

                auto tg_xxxxyzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 40); 

                auto tg_xxxxyzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 41); 

                auto tg_xxxxzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 42); 

                auto tg_xxxxzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 43); 

                auto tg_xxxxzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 44); 

                auto tg_xxxyyyyy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 45); 

                auto tg_xxxyyyyy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 46); 

                auto tg_xxxyyyyy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 47); 

                auto tg_xxxyyyyz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 48); 

                auto tg_xxxyyyyz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 49); 

                auto tg_xxxyyyyz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 50); 

                auto tg_xxxyyyzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 51); 

                auto tg_xxxyyyzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 52); 

                auto tg_xxxyyyzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 53); 

                auto tg_xxxyyzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 54); 

                auto tg_xxxyyzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 55); 

                auto tg_xxxyyzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 56); 

                auto tg_xxxyzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 57); 

                auto tg_xxxyzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 58); 

                auto tg_xxxyzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 59); 

                auto tg_xxxzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 60); 

                auto tg_xxxzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 61); 

                auto tg_xxxzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 62); 

                auto tg_xxyyyyyy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 63); 

                auto tg_xxyyyyyy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 64); 

                auto tg_xxyyyyyy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 65); 

                auto tg_xxyyyyyz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 66); 

                auto tg_xxyyyyyz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 67); 

                // Batch of Integrals (0,68)

                #pragma omp simd aligned(fxn, fza, tg_xxxxxx_x_0, tg_xxxxxx_x_1, tg_xxxxxx_y_0, tg_xxxxxx_y_1, \
                                         tg_xxxxxx_z_0, tg_xxxxxx_z_1, tg_xxxxxxx_0_1, tg_xxxxxxx_x_0, tg_xxxxxxx_x_1, \
                                         tg_xxxxxxx_y_0, tg_xxxxxxx_y_1, tg_xxxxxxx_z_0, tg_xxxxxxx_z_1, tg_xxxxxxxx_x_0, \
                                         tg_xxxxxxxx_y_0, tg_xxxxxxxx_z_0, tg_xxxxxxxy_x_0, tg_xxxxxxxy_y_0, tg_xxxxxxxy_z_0, \
                                         tg_xxxxxxxz_x_0, tg_xxxxxxxz_y_0, tg_xxxxxxxz_z_0, tg_xxxxxxy_0_1, tg_xxxxxxy_x_0, \
                                         tg_xxxxxxy_x_1, tg_xxxxxxy_y_0, tg_xxxxxxy_y_1, tg_xxxxxxy_z_0, tg_xxxxxxy_z_1, \
                                         tg_xxxxxxyy_x_0, tg_xxxxxxyy_y_0, tg_xxxxxxyy_z_0, tg_xxxxxxyz_x_0, tg_xxxxxxyz_y_0, \
                                         tg_xxxxxxyz_z_0, tg_xxxxxxz_0_1, tg_xxxxxxz_x_0, tg_xxxxxxz_x_1, tg_xxxxxxz_y_0, \
                                         tg_xxxxxxz_y_1, tg_xxxxxxz_z_0, tg_xxxxxxz_z_1, tg_xxxxxxzz_x_0, tg_xxxxxxzz_y_0, \
                                         tg_xxxxxxzz_z_0, tg_xxxxxy_x_0, tg_xxxxxy_x_1, tg_xxxxxy_y_0, tg_xxxxxy_y_1, \
                                         tg_xxxxxy_z_0, tg_xxxxxy_z_1, tg_xxxxxyy_0_1, tg_xxxxxyy_x_0, tg_xxxxxyy_x_1, \
                                         tg_xxxxxyy_y_0, tg_xxxxxyy_y_1, tg_xxxxxyy_z_0, tg_xxxxxyy_z_1, tg_xxxxxyyy_x_0, \
                                         tg_xxxxxyyy_y_0, tg_xxxxxyyy_z_0, tg_xxxxxyyz_x_0, tg_xxxxxyyz_y_0, tg_xxxxxyyz_z_0, \
                                         tg_xxxxxyz_0_1, tg_xxxxxyz_x_0, tg_xxxxxyz_x_1, tg_xxxxxyz_y_0, tg_xxxxxyz_y_1, \
                                         tg_xxxxxyz_z_0, tg_xxxxxyz_z_1, tg_xxxxxyzz_x_0, tg_xxxxxyzz_y_0, tg_xxxxxyzz_z_0, \
                                         tg_xxxxxz_x_0, tg_xxxxxz_x_1, tg_xxxxxz_y_0, tg_xxxxxz_y_1, tg_xxxxxz_z_0, \
                                         tg_xxxxxz_z_1, tg_xxxxxzz_0_1, tg_xxxxxzz_x_0, tg_xxxxxzz_x_1, tg_xxxxxzz_y_0, \
                                         tg_xxxxxzz_y_1, tg_xxxxxzz_z_0, tg_xxxxxzz_z_1, tg_xxxxxzzz_x_0, tg_xxxxxzzz_y_0, \
                                         tg_xxxxxzzz_z_0, tg_xxxxyy_x_0, tg_xxxxyy_x_1, tg_xxxxyy_y_0, tg_xxxxyy_y_1, \
                                         tg_xxxxyy_z_0, tg_xxxxyy_z_1, tg_xxxxyyy_0_1, tg_xxxxyyy_x_0, tg_xxxxyyy_x_1, \
                                         tg_xxxxyyy_y_0, tg_xxxxyyy_y_1, tg_xxxxyyy_z_0, tg_xxxxyyy_z_1, tg_xxxxyyyy_x_0, \
                                         tg_xxxxyyyy_y_0, tg_xxxxyyyy_z_0, tg_xxxxyyyz_x_0, tg_xxxxyyyz_y_0, tg_xxxxyyyz_z_0, \
                                         tg_xxxxyyz_0_1, tg_xxxxyyz_x_0, tg_xxxxyyz_x_1, tg_xxxxyyz_y_0, tg_xxxxyyz_y_1, \
                                         tg_xxxxyyz_z_0, tg_xxxxyyz_z_1, tg_xxxxyyzz_x_0, tg_xxxxyyzz_y_0, tg_xxxxyyzz_z_0, \
                                         tg_xxxxyz_x_0, tg_xxxxyz_x_1, tg_xxxxyz_y_0, tg_xxxxyz_y_1, tg_xxxxyz_z_0, \
                                         tg_xxxxyz_z_1, tg_xxxxyzz_0_1, tg_xxxxyzz_x_0, tg_xxxxyzz_x_1, tg_xxxxyzz_y_0, \
                                         tg_xxxxyzz_y_1, tg_xxxxyzz_z_0, tg_xxxxyzz_z_1, tg_xxxxyzzz_x_0, tg_xxxxyzzz_y_0, \
                                         tg_xxxxyzzz_z_0, tg_xxxxzz_x_0, tg_xxxxzz_x_1, tg_xxxxzz_y_0, tg_xxxxzz_y_1, \
                                         tg_xxxxzz_z_0, tg_xxxxzz_z_1, tg_xxxxzzz_0_1, tg_xxxxzzz_x_0, tg_xxxxzzz_x_1, \
                                         tg_xxxxzzz_y_0, tg_xxxxzzz_y_1, tg_xxxxzzz_z_0, tg_xxxxzzz_z_1, tg_xxxxzzzz_x_0, \
                                         tg_xxxxzzzz_y_0, tg_xxxxzzzz_z_0, tg_xxxyyy_x_0, tg_xxxyyy_x_1, tg_xxxyyy_y_0, \
                                         tg_xxxyyy_y_1, tg_xxxyyy_z_0, tg_xxxyyy_z_1, tg_xxxyyyy_0_1, tg_xxxyyyy_x_0, \
                                         tg_xxxyyyy_x_1, tg_xxxyyyy_y_0, tg_xxxyyyy_y_1, tg_xxxyyyy_z_0, tg_xxxyyyy_z_1, \
                                         tg_xxxyyyyy_x_0, tg_xxxyyyyy_y_0, tg_xxxyyyyy_z_0, tg_xxxyyyyz_x_0, tg_xxxyyyyz_y_0, \
                                         tg_xxxyyyyz_z_0, tg_xxxyyyz_0_1, tg_xxxyyyz_x_0, tg_xxxyyyz_x_1, tg_xxxyyyz_y_0, \
                                         tg_xxxyyyz_y_1, tg_xxxyyyz_z_0, tg_xxxyyyz_z_1, tg_xxxyyyzz_x_0, tg_xxxyyyzz_y_0, \
                                         tg_xxxyyyzz_z_0, tg_xxxyyz_x_0, tg_xxxyyz_x_1, tg_xxxyyz_y_0, tg_xxxyyz_y_1, \
                                         tg_xxxyyz_z_0, tg_xxxyyz_z_1, tg_xxxyyzz_0_1, tg_xxxyyzz_x_0, tg_xxxyyzz_x_1, \
                                         tg_xxxyyzz_y_0, tg_xxxyyzz_y_1, tg_xxxyyzz_z_0, tg_xxxyyzz_z_1, tg_xxxyyzzz_x_0, \
                                         tg_xxxyyzzz_y_0, tg_xxxyyzzz_z_0, tg_xxxyzz_x_0, tg_xxxyzz_x_1, tg_xxxyzz_y_0, \
                                         tg_xxxyzz_y_1, tg_xxxyzz_z_0, tg_xxxyzz_z_1, tg_xxxyzzz_0_1, tg_xxxyzzz_x_0, \
                                         tg_xxxyzzz_x_1, tg_xxxyzzz_y_0, tg_xxxyzzz_y_1, tg_xxxyzzz_z_0, tg_xxxyzzz_z_1, \
                                         tg_xxxyzzzz_x_0, tg_xxxyzzzz_y_0, tg_xxxyzzzz_z_0, tg_xxxzzz_x_0, tg_xxxzzz_x_1, \
                                         tg_xxxzzz_y_0, tg_xxxzzz_y_1, tg_xxxzzz_z_0, tg_xxxzzz_z_1, tg_xxxzzzz_0_1, \
                                         tg_xxxzzzz_x_0, tg_xxxzzzz_x_1, tg_xxxzzzz_y_0, tg_xxxzzzz_y_1, tg_xxxzzzz_z_0, \
                                         tg_xxxzzzz_z_1, tg_xxxzzzzz_x_0, tg_xxxzzzzz_y_0, tg_xxxzzzzz_z_0, tg_xxyyyy_x_0, \
                                         tg_xxyyyy_x_1, tg_xxyyyy_y_0, tg_xxyyyy_y_1, tg_xxyyyy_z_0, tg_xxyyyy_z_1, \
                                         tg_xxyyyyy_0_1, tg_xxyyyyy_x_0, tg_xxyyyyy_x_1, tg_xxyyyyy_y_0, tg_xxyyyyy_y_1, \
                                         tg_xxyyyyy_z_0, tg_xxyyyyy_z_1, tg_xxyyyyyy_x_0, tg_xxyyyyyy_y_0, tg_xxyyyyyy_z_0, \
                                         tg_xxyyyyyz_x_0, tg_xxyyyyyz_y_0, tg_xxyyyyz_0_1, tg_xxyyyyz_x_0, tg_xxyyyyz_x_1, \
                                         tg_xxyyyyz_y_0, tg_xxyyyyz_y_1, tg_xxyyyyz_z_0, tg_xxyyyyz_z_1, tg_xxyyyz_x_0, \
                                         tg_xxyyyz_x_1, tg_xxyyyz_y_0, tg_xxyyyz_y_1, tg_xxyyyz_z_0, tg_xxyyyz_z_1, \
                                         tg_xxyyyzz_0_1, tg_xxyyyzz_x_0, tg_xxyyyzz_x_1, tg_xxyyyzz_y_0, tg_xxyyyzz_y_1, \
                                         tg_xxyyyzz_z_0, tg_xxyyyzz_z_1, tg_xxyyzz_x_0, tg_xxyyzz_x_1, tg_xxyyzz_y_0, \
                                         tg_xxyyzz_y_1, tg_xxyyzz_z_0, tg_xxyyzz_z_1, tg_xxyyzzz_0_1, tg_xxyyzzz_x_0, \
                                         tg_xxyyzzz_x_1, tg_xxyyzzz_y_0, tg_xxyyzzz_y_1, tg_xxyyzzz_z_0, tg_xxyyzzz_z_1, \
                                         tg_xxyzzz_x_0, tg_xxyzzz_x_1, tg_xxyzzz_y_0, tg_xxyzzz_y_1, tg_xxyzzz_z_0, \
                                         tg_xxyzzz_z_1, tg_xxyzzzz_0_1, tg_xxyzzzz_x_0, tg_xxyzzzz_x_1, tg_xxyzzzz_y_0, \
                                         tg_xxyzzzz_y_1, tg_xxyzzzz_z_0, tg_xxyzzzz_z_1, tg_xxzzzz_x_0, tg_xxzzzz_x_1, \
                                         tg_xxzzzz_y_0, tg_xxzzzz_y_1, tg_xxzzzz_z_0, tg_xxzzzz_z_1, tg_xxzzzzz_0_1, \
                                         tg_xxzzzzz_x_0, tg_xxzzzzz_x_1, tg_xxzzzzz_y_0, tg_xxzzzzz_y_1, tg_xxzzzzz_z_0, \
                                         tg_xxzzzzz_z_1, tg_xyyyyy_x_0, tg_xyyyyy_x_1, tg_xyyyyy_y_0, tg_xyyyyy_y_1, \
                                         tg_xyyyyy_z_0, tg_xyyyyy_z_1, tg_xyyyyyy_0_1, tg_xyyyyyy_x_0, tg_xyyyyyy_x_1, \
                                         tg_xyyyyyy_y_0, tg_xyyyyyy_y_1, tg_xyyyyyy_z_0, tg_xyyyyyy_z_1, tg_xyyyyyz_0_1, \
                                         tg_xyyyyyz_x_0, tg_xyyyyyz_x_1, tg_xyyyyyz_y_0, tg_xyyyyyz_y_1, tg_xyyyyz_x_0, \
                                         tg_xyyyyz_x_1, tg_xyyyyz_y_0, tg_xyyyyz_y_1, tg_xyyyyz_z_0, tg_xyyyyz_z_1, \
                                         tg_xyyyzz_x_0, tg_xyyyzz_x_1, tg_xyyyzz_y_0, tg_xyyyzz_y_1, tg_xyyyzz_z_0, \
                                         tg_xyyyzz_z_1, tg_xyyzzz_x_0, tg_xyyzzz_x_1, tg_xyyzzz_y_0, tg_xyyzzz_y_1, \
                                         tg_xyyzzz_z_0, tg_xyyzzz_z_1, tg_xyzzzz_x_0, tg_xyzzzz_x_1, tg_xyzzzz_y_0, \
                                         tg_xyzzzz_y_1, tg_xyzzzz_z_0, tg_xyzzzz_z_1, tg_xzzzzz_x_0, tg_xzzzzz_x_1, \
                                         tg_xzzzzz_y_0, tg_xzzzzz_y_1, tg_xzzzzz_z_0, tg_xzzzzz_z_1, tg_yyyyyy_x_0, \
                                         tg_yyyyyy_x_1, tg_yyyyyy_y_0, tg_yyyyyy_y_1, tg_yyyyyy_z_0, tg_yyyyyy_z_1, \
                                         tg_yyyyyz_x_0, tg_yyyyyz_x_1, tg_yyyyyz_y_0, tg_yyyyyz_y_1, wp_x: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxxxxxxx_x_0[j] = pb_x * tg_xxxxxxx_x_0[j] + wp_x[j] * tg_xxxxxxx_x_1[j] + 3.5 * fl1_fx * tg_xxxxxx_x_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxxx_0_1[j];

                    tg_xxxxxxxx_y_0[j] = pb_x * tg_xxxxxxx_y_0[j] + wp_x[j] * tg_xxxxxxx_y_1[j] + 3.5 * fl1_fx * tg_xxxxxx_y_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_y_1[j];

                    tg_xxxxxxxx_z_0[j] = pb_x * tg_xxxxxxx_z_0[j] + wp_x[j] * tg_xxxxxxx_z_1[j] + 3.5 * fl1_fx * tg_xxxxxx_z_0[j] - 3.5 * fl1_fx * fl1_fza * tg_xxxxxx_z_1[j];

                    tg_xxxxxxxy_x_0[j] = pb_x * tg_xxxxxxy_x_0[j] + wp_x[j] * tg_xxxxxxy_x_1[j] + 3.0 * fl1_fx * tg_xxxxxy_x_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxxy_0_1[j];

                    tg_xxxxxxxy_y_0[j] = pb_x * tg_xxxxxxy_y_0[j] + wp_x[j] * tg_xxxxxxy_y_1[j] + 3.0 * fl1_fx * tg_xxxxxy_y_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_y_1[j];

                    tg_xxxxxxxy_z_0[j] = pb_x * tg_xxxxxxy_z_0[j] + wp_x[j] * tg_xxxxxxy_z_1[j] + 3.0 * fl1_fx * tg_xxxxxy_z_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxy_z_1[j];

                    tg_xxxxxxxz_x_0[j] = pb_x * tg_xxxxxxz_x_0[j] + wp_x[j] * tg_xxxxxxz_x_1[j] + 3.0 * fl1_fx * tg_xxxxxz_x_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxxz_0_1[j];

                    tg_xxxxxxxz_y_0[j] = pb_x * tg_xxxxxxz_y_0[j] + wp_x[j] * tg_xxxxxxz_y_1[j] + 3.0 * fl1_fx * tg_xxxxxz_y_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_y_1[j];

                    tg_xxxxxxxz_z_0[j] = pb_x * tg_xxxxxxz_z_0[j] + wp_x[j] * tg_xxxxxxz_z_1[j] + 3.0 * fl1_fx * tg_xxxxxz_z_0[j] - 3.0 * fl1_fx * fl1_fza * tg_xxxxxz_z_1[j];

                    tg_xxxxxxyy_x_0[j] = pb_x * tg_xxxxxyy_x_0[j] + wp_x[j] * tg_xxxxxyy_x_1[j] + 2.5 * fl1_fx * tg_xxxxyy_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxyy_0_1[j];

                    tg_xxxxxxyy_y_0[j] = pb_x * tg_xxxxxyy_y_0[j] + wp_x[j] * tg_xxxxxyy_y_1[j] + 2.5 * fl1_fx * tg_xxxxyy_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_y_1[j];

                    tg_xxxxxxyy_z_0[j] = pb_x * tg_xxxxxyy_z_0[j] + wp_x[j] * tg_xxxxxyy_z_1[j] + 2.5 * fl1_fx * tg_xxxxyy_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyy_z_1[j];

                    tg_xxxxxxyz_x_0[j] = pb_x * tg_xxxxxyz_x_0[j] + wp_x[j] * tg_xxxxxyz_x_1[j] + 2.5 * fl1_fx * tg_xxxxyz_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxyz_0_1[j];

                    tg_xxxxxxyz_y_0[j] = pb_x * tg_xxxxxyz_y_0[j] + wp_x[j] * tg_xxxxxyz_y_1[j] + 2.5 * fl1_fx * tg_xxxxyz_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_y_1[j];

                    tg_xxxxxxyz_z_0[j] = pb_x * tg_xxxxxyz_z_0[j] + wp_x[j] * tg_xxxxxyz_z_1[j] + 2.5 * fl1_fx * tg_xxxxyz_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxyz_z_1[j];

                    tg_xxxxxxzz_x_0[j] = pb_x * tg_xxxxxzz_x_0[j] + wp_x[j] * tg_xxxxxzz_x_1[j] + 2.5 * fl1_fx * tg_xxxxzz_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxxzz_0_1[j];

                    tg_xxxxxxzz_y_0[j] = pb_x * tg_xxxxxzz_y_0[j] + wp_x[j] * tg_xxxxxzz_y_1[j] + 2.5 * fl1_fx * tg_xxxxzz_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_y_1[j];

                    tg_xxxxxxzz_z_0[j] = pb_x * tg_xxxxxzz_z_0[j] + wp_x[j] * tg_xxxxxzz_z_1[j] + 2.5 * fl1_fx * tg_xxxxzz_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_xxxxzz_z_1[j];

                    tg_xxxxxyyy_x_0[j] = pb_x * tg_xxxxyyy_x_0[j] + wp_x[j] * tg_xxxxyyy_x_1[j] + 2.0 * fl1_fx * tg_xxxyyy_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_x_1[j] + 0.5 * fl1_fxn * tg_xxxxyyy_0_1[j];

                    tg_xxxxxyyy_y_0[j] = pb_x * tg_xxxxyyy_y_0[j] + wp_x[j] * tg_xxxxyyy_y_1[j] + 2.0 * fl1_fx * tg_xxxyyy_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_y_1[j];

                    tg_xxxxxyyy_z_0[j] = pb_x * tg_xxxxyyy_z_0[j] + wp_x[j] * tg_xxxxyyy_z_1[j] + 2.0 * fl1_fx * tg_xxxyyy_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyy_z_1[j];

                    tg_xxxxxyyz_x_0[j] = pb_x * tg_xxxxyyz_x_0[j] + wp_x[j] * tg_xxxxyyz_x_1[j] + 2.0 * fl1_fx * tg_xxxyyz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxyyz_0_1[j];

                    tg_xxxxxyyz_y_0[j] = pb_x * tg_xxxxyyz_y_0[j] + wp_x[j] * tg_xxxxyyz_y_1[j] + 2.0 * fl1_fx * tg_xxxyyz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_y_1[j];

                    tg_xxxxxyyz_z_0[j] = pb_x * tg_xxxxyyz_z_0[j] + wp_x[j] * tg_xxxxyyz_z_1[j] + 2.0 * fl1_fx * tg_xxxyyz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyyz_z_1[j];

                    tg_xxxxxyzz_x_0[j] = pb_x * tg_xxxxyzz_x_0[j] + wp_x[j] * tg_xxxxyzz_x_1[j] + 2.0 * fl1_fx * tg_xxxyzz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxyzz_0_1[j];

                    tg_xxxxxyzz_y_0[j] = pb_x * tg_xxxxyzz_y_0[j] + wp_x[j] * tg_xxxxyzz_y_1[j] + 2.0 * fl1_fx * tg_xxxyzz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_y_1[j];

                    tg_xxxxxyzz_z_0[j] = pb_x * tg_xxxxyzz_z_0[j] + wp_x[j] * tg_xxxxyzz_z_1[j] + 2.0 * fl1_fx * tg_xxxyzz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxyzz_z_1[j];

                    tg_xxxxxzzz_x_0[j] = pb_x * tg_xxxxzzz_x_0[j] + wp_x[j] * tg_xxxxzzz_x_1[j] + 2.0 * fl1_fx * tg_xxxzzz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxxzzz_0_1[j];

                    tg_xxxxxzzz_y_0[j] = pb_x * tg_xxxxzzz_y_0[j] + wp_x[j] * tg_xxxxzzz_y_1[j] + 2.0 * fl1_fx * tg_xxxzzz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_y_1[j];

                    tg_xxxxxzzz_z_0[j] = pb_x * tg_xxxxzzz_z_0[j] + wp_x[j] * tg_xxxxzzz_z_1[j] + 2.0 * fl1_fx * tg_xxxzzz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_xxxzzz_z_1[j];

                    tg_xxxxyyyy_x_0[j] = pb_x * tg_xxxyyyy_x_0[j] + wp_x[j] * tg_xxxyyyy_x_1[j] + 1.5 * fl1_fx * tg_xxyyyy_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_x_1[j] + 0.5 * fl1_fxn * tg_xxxyyyy_0_1[j];

                    tg_xxxxyyyy_y_0[j] = pb_x * tg_xxxyyyy_y_0[j] + wp_x[j] * tg_xxxyyyy_y_1[j] + 1.5 * fl1_fx * tg_xxyyyy_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_y_1[j];

                    tg_xxxxyyyy_z_0[j] = pb_x * tg_xxxyyyy_z_0[j] + wp_x[j] * tg_xxxyyyy_z_1[j] + 1.5 * fl1_fx * tg_xxyyyy_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyy_z_1[j];

                    tg_xxxxyyyz_x_0[j] = pb_x * tg_xxxyyyz_x_0[j] + wp_x[j] * tg_xxxyyyz_x_1[j] + 1.5 * fl1_fx * tg_xxyyyz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_x_1[j] + 0.5 * fl1_fxn * tg_xxxyyyz_0_1[j];

                    tg_xxxxyyyz_y_0[j] = pb_x * tg_xxxyyyz_y_0[j] + wp_x[j] * tg_xxxyyyz_y_1[j] + 1.5 * fl1_fx * tg_xxyyyz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_y_1[j];

                    tg_xxxxyyyz_z_0[j] = pb_x * tg_xxxyyyz_z_0[j] + wp_x[j] * tg_xxxyyyz_z_1[j] + 1.5 * fl1_fx * tg_xxyyyz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyyz_z_1[j];

                    tg_xxxxyyzz_x_0[j] = pb_x * tg_xxxyyzz_x_0[j] + wp_x[j] * tg_xxxyyzz_x_1[j] + 1.5 * fl1_fx * tg_xxyyzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxyyzz_0_1[j];

                    tg_xxxxyyzz_y_0[j] = pb_x * tg_xxxyyzz_y_0[j] + wp_x[j] * tg_xxxyyzz_y_1[j] + 1.5 * fl1_fx * tg_xxyyzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_y_1[j];

                    tg_xxxxyyzz_z_0[j] = pb_x * tg_xxxyyzz_z_0[j] + wp_x[j] * tg_xxxyyzz_z_1[j] + 1.5 * fl1_fx * tg_xxyyzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyyzz_z_1[j];

                    tg_xxxxyzzz_x_0[j] = pb_x * tg_xxxyzzz_x_0[j] + wp_x[j] * tg_xxxyzzz_x_1[j] + 1.5 * fl1_fx * tg_xxyzzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxyzzz_0_1[j];

                    tg_xxxxyzzz_y_0[j] = pb_x * tg_xxxyzzz_y_0[j] + wp_x[j] * tg_xxxyzzz_y_1[j] + 1.5 * fl1_fx * tg_xxyzzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_y_1[j];

                    tg_xxxxyzzz_z_0[j] = pb_x * tg_xxxyzzz_z_0[j] + wp_x[j] * tg_xxxyzzz_z_1[j] + 1.5 * fl1_fx * tg_xxyzzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxyzzz_z_1[j];

                    tg_xxxxzzzz_x_0[j] = pb_x * tg_xxxzzzz_x_0[j] + wp_x[j] * tg_xxxzzzz_x_1[j] + 1.5 * fl1_fx * tg_xxzzzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxxzzzz_0_1[j];

                    tg_xxxxzzzz_y_0[j] = pb_x * tg_xxxzzzz_y_0[j] + wp_x[j] * tg_xxxzzzz_y_1[j] + 1.5 * fl1_fx * tg_xxzzzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_y_1[j];

                    tg_xxxxzzzz_z_0[j] = pb_x * tg_xxxzzzz_z_0[j] + wp_x[j] * tg_xxxzzzz_z_1[j] + 1.5 * fl1_fx * tg_xxzzzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_xxzzzz_z_1[j];

                    tg_xxxyyyyy_x_0[j] = pb_x * tg_xxyyyyy_x_0[j] + wp_x[j] * tg_xxyyyyy_x_1[j] + fl1_fx * tg_xyyyyy_x_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_x_1[j] + 0.5 * fl1_fxn * tg_xxyyyyy_0_1[j];

                    tg_xxxyyyyy_y_0[j] = pb_x * tg_xxyyyyy_y_0[j] + wp_x[j] * tg_xxyyyyy_y_1[j] + fl1_fx * tg_xyyyyy_y_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_y_1[j];

                    tg_xxxyyyyy_z_0[j] = pb_x * tg_xxyyyyy_z_0[j] + wp_x[j] * tg_xxyyyyy_z_1[j] + fl1_fx * tg_xyyyyy_z_0[j] - fl1_fx * fl1_fza * tg_xyyyyy_z_1[j];

                    tg_xxxyyyyz_x_0[j] = pb_x * tg_xxyyyyz_x_0[j] + wp_x[j] * tg_xxyyyyz_x_1[j] + fl1_fx * tg_xyyyyz_x_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_x_1[j] + 0.5 * fl1_fxn * tg_xxyyyyz_0_1[j];

                    tg_xxxyyyyz_y_0[j] = pb_x * tg_xxyyyyz_y_0[j] + wp_x[j] * tg_xxyyyyz_y_1[j] + fl1_fx * tg_xyyyyz_y_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_y_1[j];

                    tg_xxxyyyyz_z_0[j] = pb_x * tg_xxyyyyz_z_0[j] + wp_x[j] * tg_xxyyyyz_z_1[j] + fl1_fx * tg_xyyyyz_z_0[j] - fl1_fx * fl1_fza * tg_xyyyyz_z_1[j];

                    tg_xxxyyyzz_x_0[j] = pb_x * tg_xxyyyzz_x_0[j] + wp_x[j] * tg_xxyyyzz_x_1[j] + fl1_fx * tg_xyyyzz_x_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_x_1[j] + 0.5 * fl1_fxn * tg_xxyyyzz_0_1[j];

                    tg_xxxyyyzz_y_0[j] = pb_x * tg_xxyyyzz_y_0[j] + wp_x[j] * tg_xxyyyzz_y_1[j] + fl1_fx * tg_xyyyzz_y_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_y_1[j];

                    tg_xxxyyyzz_z_0[j] = pb_x * tg_xxyyyzz_z_0[j] + wp_x[j] * tg_xxyyyzz_z_1[j] + fl1_fx * tg_xyyyzz_z_0[j] - fl1_fx * fl1_fza * tg_xyyyzz_z_1[j];

                    tg_xxxyyzzz_x_0[j] = pb_x * tg_xxyyzzz_x_0[j] + wp_x[j] * tg_xxyyzzz_x_1[j] + fl1_fx * tg_xyyzzz_x_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxyyzzz_0_1[j];

                    tg_xxxyyzzz_y_0[j] = pb_x * tg_xxyyzzz_y_0[j] + wp_x[j] * tg_xxyyzzz_y_1[j] + fl1_fx * tg_xyyzzz_y_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_y_1[j];

                    tg_xxxyyzzz_z_0[j] = pb_x * tg_xxyyzzz_z_0[j] + wp_x[j] * tg_xxyyzzz_z_1[j] + fl1_fx * tg_xyyzzz_z_0[j] - fl1_fx * fl1_fza * tg_xyyzzz_z_1[j];

                    tg_xxxyzzzz_x_0[j] = pb_x * tg_xxyzzzz_x_0[j] + wp_x[j] * tg_xxyzzzz_x_1[j] + fl1_fx * tg_xyzzzz_x_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxyzzzz_0_1[j];

                    tg_xxxyzzzz_y_0[j] = pb_x * tg_xxyzzzz_y_0[j] + wp_x[j] * tg_xxyzzzz_y_1[j] + fl1_fx * tg_xyzzzz_y_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_y_1[j];

                    tg_xxxyzzzz_z_0[j] = pb_x * tg_xxyzzzz_z_0[j] + wp_x[j] * tg_xxyzzzz_z_1[j] + fl1_fx * tg_xyzzzz_z_0[j] - fl1_fx * fl1_fza * tg_xyzzzz_z_1[j];

                    tg_xxxzzzzz_x_0[j] = pb_x * tg_xxzzzzz_x_0[j] + wp_x[j] * tg_xxzzzzz_x_1[j] + fl1_fx * tg_xzzzzz_x_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xxzzzzz_0_1[j];

                    tg_xxxzzzzz_y_0[j] = pb_x * tg_xxzzzzz_y_0[j] + wp_x[j] * tg_xxzzzzz_y_1[j] + fl1_fx * tg_xzzzzz_y_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_y_1[j];

                    tg_xxxzzzzz_z_0[j] = pb_x * tg_xxzzzzz_z_0[j] + wp_x[j] * tg_xxzzzzz_z_1[j] + fl1_fx * tg_xzzzzz_z_0[j] - fl1_fx * fl1_fza * tg_xzzzzz_z_1[j];

                    tg_xxyyyyyy_x_0[j] = pb_x * tg_xyyyyyy_x_0[j] + wp_x[j] * tg_xyyyyyy_x_1[j] + 0.5 * fl1_fx * tg_yyyyyy_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_x_1[j] + 0.5 * fl1_fxn * tg_xyyyyyy_0_1[j];

                    tg_xxyyyyyy_y_0[j] = pb_x * tg_xyyyyyy_y_0[j] + wp_x[j] * tg_xyyyyyy_y_1[j] + 0.5 * fl1_fx * tg_yyyyyy_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_y_1[j];

                    tg_xxyyyyyy_z_0[j] = pb_x * tg_xyyyyyy_z_0[j] + wp_x[j] * tg_xyyyyyy_z_1[j] + 0.5 * fl1_fx * tg_yyyyyy_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyy_z_1[j];

                    tg_xxyyyyyz_x_0[j] = pb_x * tg_xyyyyyz_x_0[j] + wp_x[j] * tg_xyyyyyz_x_1[j] + 0.5 * fl1_fx * tg_yyyyyz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_x_1[j] + 0.5 * fl1_fxn * tg_xyyyyyz_0_1[j];

                    tg_xxyyyyyz_y_0[j] = pb_x * tg_xyyyyyz_y_0[j] + wp_x[j] * tg_xyyyyyz_y_1[j] + 0.5 * fl1_fx * tg_yyyyyz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_y_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSP_68_135(      CMemBlock2D<double>& primBuffer,
                                        const CRecursionMap&       recursionMap,
                                        const CMemBlock2D<double>& osFactors,
                                        const CMemBlock2D<double>& wpDistances,
                                        const CGtoPairsBlock&      braGtoPairsBlock,
                                        const CGtoPairsBlock&      ketGtoPairsBlock,
                                        const int32_t              nKetPrimPairs,
                                        const int32_t              iContrPair)
    {
        // Batch of Integrals (68,135)

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
                                             {8, -1, -1, -1},
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {8, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_7_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_6_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_g_6_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {6, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_g_7_0_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electron Repulsion"}, 0, true, 
                                                             {7, -1, -1, -1}, {0, -1, -1, -1}, 
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

                auto tg_xyyyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 68); 

                auto tg_xyyyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 69); 

                auto tg_xyyyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 70); 

                auto tg_xyyyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 71); 

                auto tg_xyyyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 72); 

                auto tg_xyyyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 73); 

                auto tg_xyyyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 74); 

                auto tg_xyyzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 75); 

                auto tg_xyyzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 76); 

                auto tg_xyyzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 77); 

                auto tg_xyzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 78); 

                auto tg_xyzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 79); 

                auto tg_xyzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 80); 

                auto tg_xzzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 81); 

                auto tg_xzzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 82); 

                auto tg_xzzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 83); 

                auto tg_yyyyyyy_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 84); 

                auto tg_yyyyyyy_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 85); 

                auto tg_yyyyyyy_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 86); 

                auto tg_yyyyyyz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 87); 

                auto tg_yyyyyyz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 88); 

                auto tg_yyyyyyz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 89); 

                auto tg_yyyyyzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 90); 

                auto tg_yyyyyzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 91); 

                auto tg_yyyyyzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 92); 

                auto tg_yyyyzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 93); 

                auto tg_yyyyzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 94); 

                auto tg_yyyyzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 95); 

                auto tg_yyyzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 96); 

                auto tg_yyyzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 97); 

                auto tg_yyyzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 98); 

                auto tg_yyzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 99); 

                auto tg_yyzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 100); 

                auto tg_yyzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 101); 

                auto tg_yzzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 102); 

                auto tg_yzzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 103); 

                auto tg_yzzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 104); 

                auto tg_zzzzzzz_x_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 105); 

                auto tg_zzzzzzz_y_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 106); 

                auto tg_zzzzzzz_z_0 = primBuffer.data(pidx_g_7_1_m0 + 108 * idx + 107); 

                auto tg_xyyyyyz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 68); 

                auto tg_xyyyyzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 69); 

                auto tg_xyyyyzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 70); 

                auto tg_xyyyyzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 71); 

                auto tg_xyyyzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 72); 

                auto tg_xyyyzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 73); 

                auto tg_xyyyzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 74); 

                auto tg_xyyzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 75); 

                auto tg_xyyzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 76); 

                auto tg_xyyzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 77); 

                auto tg_xyzzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 78); 

                auto tg_xyzzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 79); 

                auto tg_xyzzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 80); 

                auto tg_xzzzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 81); 

                auto tg_xzzzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 82); 

                auto tg_xzzzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 83); 

                auto tg_yyyyyyy_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 84); 

                auto tg_yyyyyyy_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 85); 

                auto tg_yyyyyyy_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 86); 

                auto tg_yyyyyyz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 87); 

                auto tg_yyyyyyz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 88); 

                auto tg_yyyyyyz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 89); 

                auto tg_yyyyyzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 90); 

                auto tg_yyyyyzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 91); 

                auto tg_yyyyyzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 92); 

                auto tg_yyyyzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 93); 

                auto tg_yyyyzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 94); 

                auto tg_yyyyzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 95); 

                auto tg_yyyzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 96); 

                auto tg_yyyzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 97); 

                auto tg_yyyzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 98); 

                auto tg_yyzzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 99); 

                auto tg_yyzzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 100); 

                auto tg_yyzzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 101); 

                auto tg_yzzzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 102); 

                auto tg_yzzzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 103); 

                auto tg_yzzzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 104); 

                auto tg_zzzzzzz_x_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 105); 

                auto tg_zzzzzzz_y_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 106); 

                auto tg_zzzzzzz_z_1 = primBuffer.data(pidx_g_7_1_m1 + 108 * idx + 107); 

                auto tg_yyyyyy_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 63); 

                auto tg_yyyyyy_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 64); 

                auto tg_yyyyyy_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 65); 

                auto tg_yyyyyz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 66); 

                auto tg_yyyyyz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 67); 

                auto tg_yyyyyz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 68); 

                auto tg_yyyyzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 69); 

                auto tg_yyyyzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 70); 

                auto tg_yyyyzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 71); 

                auto tg_yyyzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 72); 

                auto tg_yyyzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 73); 

                auto tg_yyyzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 74); 

                auto tg_yyzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 75); 

                auto tg_yyzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 76); 

                auto tg_yyzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 77); 

                auto tg_yzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 78); 

                auto tg_yzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 79); 

                auto tg_yzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 80); 

                auto tg_zzzzzz_x_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 81); 

                auto tg_zzzzzz_y_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 82); 

                auto tg_zzzzzz_z_0 = primBuffer.data(pidx_g_6_1_m0 + 84 * idx + 83); 

                auto tg_yyyyyy_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 63); 

                auto tg_yyyyyy_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 64); 

                auto tg_yyyyyy_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 65); 

                auto tg_yyyyyz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 66); 

                auto tg_yyyyyz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 67); 

                auto tg_yyyyyz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 68); 

                auto tg_yyyyzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 69); 

                auto tg_yyyyzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 70); 

                auto tg_yyyyzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 71); 

                auto tg_yyyzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 72); 

                auto tg_yyyzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 73); 

                auto tg_yyyzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 74); 

                auto tg_yyzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 75); 

                auto tg_yyzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 76); 

                auto tg_yyzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 77); 

                auto tg_yzzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 78); 

                auto tg_yzzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 79); 

                auto tg_yzzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 80); 

                auto tg_zzzzzz_x_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 81); 

                auto tg_zzzzzz_y_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 82); 

                auto tg_zzzzzz_z_1 = primBuffer.data(pidx_g_6_1_m1 + 84 * idx + 83); 

                auto tg_xyyyyzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 23); 

                auto tg_xyyyzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 24); 

                auto tg_xyyzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 25); 

                auto tg_xyzzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 26); 

                auto tg_xzzzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 27); 

                auto tg_yyyyyyy_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 28); 

                auto tg_yyyyyyz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 29); 

                auto tg_yyyyyzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 30); 

                auto tg_yyyyzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 31); 

                auto tg_yyyzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 32); 

                auto tg_yyzzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 33); 

                auto tg_yzzzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 34); 

                auto tg_zzzzzzz_0_1 = primBuffer.data(pidx_g_7_0_m1 + 36 * idx + 35); 

                // set up pointers to integrals

                auto tg_xxyyyyyz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 68); 

                auto tg_xxyyyyzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 69); 

                auto tg_xxyyyyzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 70); 

                auto tg_xxyyyyzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 71); 

                auto tg_xxyyyzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 72); 

                auto tg_xxyyyzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 73); 

                auto tg_xxyyyzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 74); 

                auto tg_xxyyzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 75); 

                auto tg_xxyyzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 76); 

                auto tg_xxyyzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 77); 

                auto tg_xxyzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 78); 

                auto tg_xxyzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 79); 

                auto tg_xxyzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 80); 

                auto tg_xxzzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 81); 

                auto tg_xxzzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 82); 

                auto tg_xxzzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 83); 

                auto tg_xyyyyyyy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 84); 

                auto tg_xyyyyyyy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 85); 

                auto tg_xyyyyyyy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 86); 

                auto tg_xyyyyyyz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 87); 

                auto tg_xyyyyyyz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 88); 

                auto tg_xyyyyyyz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 89); 

                auto tg_xyyyyyzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 90); 

                auto tg_xyyyyyzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 91); 

                auto tg_xyyyyyzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 92); 

                auto tg_xyyyyzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 93); 

                auto tg_xyyyyzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 94); 

                auto tg_xyyyyzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 95); 

                auto tg_xyyyzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 96); 

                auto tg_xyyyzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 97); 

                auto tg_xyyyzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 98); 

                auto tg_xyyzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 99); 

                auto tg_xyyzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 100); 

                auto tg_xyyzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 101); 

                auto tg_xyzzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 102); 

                auto tg_xyzzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 103); 

                auto tg_xyzzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 104); 

                auto tg_xzzzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 105); 

                auto tg_xzzzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 106); 

                auto tg_xzzzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 107); 

                auto tg_yyyyyyyy_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 108); 

                auto tg_yyyyyyyy_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 109); 

                auto tg_yyyyyyyy_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 110); 

                auto tg_yyyyyyyz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 111); 

                auto tg_yyyyyyyz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 112); 

                auto tg_yyyyyyyz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 113); 

                auto tg_yyyyyyzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 114); 

                auto tg_yyyyyyzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 115); 

                auto tg_yyyyyyzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 116); 

                auto tg_yyyyyzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 117); 

                auto tg_yyyyyzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 118); 

                auto tg_yyyyyzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 119); 

                auto tg_yyyyzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 120); 

                auto tg_yyyyzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 121); 

                auto tg_yyyyzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 122); 

                auto tg_yyyzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 123); 

                auto tg_yyyzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 124); 

                auto tg_yyyzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 125); 

                auto tg_yyzzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 126); 

                auto tg_yyzzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 127); 

                auto tg_yyzzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 128); 

                auto tg_yzzzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 129); 

                auto tg_yzzzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 130); 

                auto tg_yzzzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 131); 

                auto tg_zzzzzzzz_x_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 132); 

                auto tg_zzzzzzzz_y_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 133); 

                auto tg_zzzzzzzz_z_0 = primBuffer.data(pidx_g_8_1_m0 + 135 * idx + 134); 

                // Batch of Integrals (68,135)

                #pragma omp simd aligned(fxn, fza, tg_xxyyyyyz_z_0, tg_xxyyyyzz_x_0, tg_xxyyyyzz_y_0, \
                                         tg_xxyyyyzz_z_0, tg_xxyyyzzz_x_0, tg_xxyyyzzz_y_0, tg_xxyyyzzz_z_0, tg_xxyyzzzz_x_0, \
                                         tg_xxyyzzzz_y_0, tg_xxyyzzzz_z_0, tg_xxyzzzzz_x_0, tg_xxyzzzzz_y_0, tg_xxyzzzzz_z_0, \
                                         tg_xxzzzzzz_x_0, tg_xxzzzzzz_y_0, tg_xxzzzzzz_z_0, tg_xyyyyyyy_x_0, tg_xyyyyyyy_y_0, \
                                         tg_xyyyyyyy_z_0, tg_xyyyyyyz_x_0, tg_xyyyyyyz_y_0, tg_xyyyyyyz_z_0, tg_xyyyyyz_z_0, \
                                         tg_xyyyyyz_z_1, tg_xyyyyyzz_x_0, tg_xyyyyyzz_y_0, tg_xyyyyyzz_z_0, tg_xyyyyzz_0_1, \
                                         tg_xyyyyzz_x_0, tg_xyyyyzz_x_1, tg_xyyyyzz_y_0, tg_xyyyyzz_y_1, tg_xyyyyzz_z_0, \
                                         tg_xyyyyzz_z_1, tg_xyyyyzzz_x_0, tg_xyyyyzzz_y_0, tg_xyyyyzzz_z_0, tg_xyyyzzz_0_1, \
                                         tg_xyyyzzz_x_0, tg_xyyyzzz_x_1, tg_xyyyzzz_y_0, tg_xyyyzzz_y_1, tg_xyyyzzz_z_0, \
                                         tg_xyyyzzz_z_1, tg_xyyyzzzz_x_0, tg_xyyyzzzz_y_0, tg_xyyyzzzz_z_0, tg_xyyzzzz_0_1, \
                                         tg_xyyzzzz_x_0, tg_xyyzzzz_x_1, tg_xyyzzzz_y_0, tg_xyyzzzz_y_1, tg_xyyzzzz_z_0, \
                                         tg_xyyzzzz_z_1, tg_xyyzzzzz_x_0, tg_xyyzzzzz_y_0, tg_xyyzzzzz_z_0, tg_xyzzzzz_0_1, \
                                         tg_xyzzzzz_x_0, tg_xyzzzzz_x_1, tg_xyzzzzz_y_0, tg_xyzzzzz_y_1, tg_xyzzzzz_z_0, \
                                         tg_xyzzzzz_z_1, tg_xyzzzzzz_x_0, tg_xyzzzzzz_y_0, tg_xyzzzzzz_z_0, tg_xzzzzzz_0_1, \
                                         tg_xzzzzzz_x_0, tg_xzzzzzz_x_1, tg_xzzzzzz_y_0, tg_xzzzzzz_y_1, tg_xzzzzzz_z_0, \
                                         tg_xzzzzzz_z_1, tg_xzzzzzzz_x_0, tg_xzzzzzzz_y_0, tg_xzzzzzzz_z_0, tg_yyyyyy_x_0, \
                                         tg_yyyyyy_x_1, tg_yyyyyy_y_0, tg_yyyyyy_y_1, tg_yyyyyy_z_0, tg_yyyyyy_z_1, \
                                         tg_yyyyyyy_0_1, tg_yyyyyyy_x_0, tg_yyyyyyy_x_1, tg_yyyyyyy_y_0, tg_yyyyyyy_y_1, \
                                         tg_yyyyyyy_z_0, tg_yyyyyyy_z_1, tg_yyyyyyyy_x_0, tg_yyyyyyyy_y_0, tg_yyyyyyyy_z_0, \
                                         tg_yyyyyyyz_x_0, tg_yyyyyyyz_y_0, tg_yyyyyyyz_z_0, tg_yyyyyyz_0_1, tg_yyyyyyz_x_0, \
                                         tg_yyyyyyz_x_1, tg_yyyyyyz_y_0, tg_yyyyyyz_y_1, tg_yyyyyyz_z_0, tg_yyyyyyz_z_1, \
                                         tg_yyyyyyzz_x_0, tg_yyyyyyzz_y_0, tg_yyyyyyzz_z_0, tg_yyyyyz_x_0, tg_yyyyyz_x_1, \
                                         tg_yyyyyz_y_0, tg_yyyyyz_y_1, tg_yyyyyz_z_0, tg_yyyyyz_z_1, tg_yyyyyzz_0_1, \
                                         tg_yyyyyzz_x_0, tg_yyyyyzz_x_1, tg_yyyyyzz_y_0, tg_yyyyyzz_y_1, tg_yyyyyzz_z_0, \
                                         tg_yyyyyzz_z_1, tg_yyyyyzzz_x_0, tg_yyyyyzzz_y_0, tg_yyyyyzzz_z_0, tg_yyyyzz_x_0, \
                                         tg_yyyyzz_x_1, tg_yyyyzz_y_0, tg_yyyyzz_y_1, tg_yyyyzz_z_0, tg_yyyyzz_z_1, \
                                         tg_yyyyzzz_0_1, tg_yyyyzzz_x_0, tg_yyyyzzz_x_1, tg_yyyyzzz_y_0, tg_yyyyzzz_y_1, \
                                         tg_yyyyzzz_z_0, tg_yyyyzzz_z_1, tg_yyyyzzzz_x_0, tg_yyyyzzzz_y_0, tg_yyyyzzzz_z_0, \
                                         tg_yyyzzz_x_0, tg_yyyzzz_x_1, tg_yyyzzz_y_0, tg_yyyzzz_y_1, tg_yyyzzz_z_0, \
                                         tg_yyyzzz_z_1, tg_yyyzzzz_0_1, tg_yyyzzzz_x_0, tg_yyyzzzz_x_1, tg_yyyzzzz_y_0, \
                                         tg_yyyzzzz_y_1, tg_yyyzzzz_z_0, tg_yyyzzzz_z_1, tg_yyyzzzzz_x_0, tg_yyyzzzzz_y_0, \
                                         tg_yyyzzzzz_z_0, tg_yyzzzz_x_0, tg_yyzzzz_x_1, tg_yyzzzz_y_0, tg_yyzzzz_y_1, \
                                         tg_yyzzzz_z_0, tg_yyzzzz_z_1, tg_yyzzzzz_0_1, tg_yyzzzzz_x_0, tg_yyzzzzz_x_1, \
                                         tg_yyzzzzz_y_0, tg_yyzzzzz_y_1, tg_yyzzzzz_z_0, tg_yyzzzzz_z_1, tg_yyzzzzzz_x_0, \
                                         tg_yyzzzzzz_y_0, tg_yyzzzzzz_z_0, tg_yzzzzz_x_0, tg_yzzzzz_x_1, tg_yzzzzz_y_0, \
                                         tg_yzzzzz_y_1, tg_yzzzzz_z_0, tg_yzzzzz_z_1, tg_yzzzzzz_0_1, tg_yzzzzzz_x_0, \
                                         tg_yzzzzzz_x_1, tg_yzzzzzz_y_0, tg_yzzzzzz_y_1, tg_yzzzzzz_z_0, tg_yzzzzzz_z_1, \
                                         tg_yzzzzzzz_x_0, tg_yzzzzzzz_y_0, tg_yzzzzzzz_z_0, tg_zzzzzz_x_0, tg_zzzzzz_x_1, \
                                         tg_zzzzzz_y_0, tg_zzzzzz_y_1, tg_zzzzzz_z_0, tg_zzzzzz_z_1, tg_zzzzzzz_0_1, \
                                         tg_zzzzzzz_x_0, tg_zzzzzzz_x_1, tg_zzzzzzz_y_0, tg_zzzzzzz_y_1, tg_zzzzzzz_z_0, \
                                         tg_zzzzzzz_z_1, tg_zzzzzzzz_x_0, tg_zzzzzzzz_y_0, tg_zzzzzzzz_z_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fxn = fxn[j];

                    double fl1_fza = fza[j];

                    tg_xxyyyyyz_z_0[j] = pb_x * tg_xyyyyyz_z_0[j] + wp_x[j] * tg_xyyyyyz_z_1[j] + 0.5 * fl1_fx * tg_yyyyyz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyyz_z_1[j];

                    tg_xxyyyyzz_x_0[j] = pb_x * tg_xyyyyzz_x_0[j] + wp_x[j] * tg_xyyyyzz_x_1[j] + 0.5 * fl1_fx * tg_yyyyzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_x_1[j] + 0.5 * fl1_fxn * tg_xyyyyzz_0_1[j];

                    tg_xxyyyyzz_y_0[j] = pb_x * tg_xyyyyzz_y_0[j] + wp_x[j] * tg_xyyyyzz_y_1[j] + 0.5 * fl1_fx * tg_yyyyzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_y_1[j];

                    tg_xxyyyyzz_z_0[j] = pb_x * tg_xyyyyzz_z_0[j] + wp_x[j] * tg_xyyyyzz_z_1[j] + 0.5 * fl1_fx * tg_yyyyzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyyzz_z_1[j];

                    tg_xxyyyzzz_x_0[j] = pb_x * tg_xyyyzzz_x_0[j] + wp_x[j] * tg_xyyyzzz_x_1[j] + 0.5 * fl1_fx * tg_yyyzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_x_1[j] + 0.5 * fl1_fxn * tg_xyyyzzz_0_1[j];

                    tg_xxyyyzzz_y_0[j] = pb_x * tg_xyyyzzz_y_0[j] + wp_x[j] * tg_xyyyzzz_y_1[j] + 0.5 * fl1_fx * tg_yyyzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_y_1[j];

                    tg_xxyyyzzz_z_0[j] = pb_x * tg_xyyyzzz_z_0[j] + wp_x[j] * tg_xyyyzzz_z_1[j] + 0.5 * fl1_fx * tg_yyyzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyyzzz_z_1[j];

                    tg_xxyyzzzz_x_0[j] = pb_x * tg_xyyzzzz_x_0[j] + wp_x[j] * tg_xyyzzzz_x_1[j] + 0.5 * fl1_fx * tg_yyzzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xyyzzzz_0_1[j];

                    tg_xxyyzzzz_y_0[j] = pb_x * tg_xyyzzzz_y_0[j] + wp_x[j] * tg_xyyzzzz_y_1[j] + 0.5 * fl1_fx * tg_yyzzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_y_1[j];

                    tg_xxyyzzzz_z_0[j] = pb_x * tg_xyyzzzz_z_0[j] + wp_x[j] * tg_xyyzzzz_z_1[j] + 0.5 * fl1_fx * tg_yyzzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yyzzzz_z_1[j];

                    tg_xxyzzzzz_x_0[j] = pb_x * tg_xyzzzzz_x_0[j] + wp_x[j] * tg_xyzzzzz_x_1[j] + 0.5 * fl1_fx * tg_yzzzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xyzzzzz_0_1[j];

                    tg_xxyzzzzz_y_0[j] = pb_x * tg_xyzzzzz_y_0[j] + wp_x[j] * tg_xyzzzzz_y_1[j] + 0.5 * fl1_fx * tg_yzzzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_y_1[j];

                    tg_xxyzzzzz_z_0[j] = pb_x * tg_xyzzzzz_z_0[j] + wp_x[j] * tg_xyzzzzz_z_1[j] + 0.5 * fl1_fx * tg_yzzzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_yzzzzz_z_1[j];

                    tg_xxzzzzzz_x_0[j] = pb_x * tg_xzzzzzz_x_0[j] + wp_x[j] * tg_xzzzzzz_x_1[j] + 0.5 * fl1_fx * tg_zzzzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_xzzzzzz_0_1[j];

                    tg_xxzzzzzz_y_0[j] = pb_x * tg_xzzzzzz_y_0[j] + wp_x[j] * tg_xzzzzzz_y_1[j] + 0.5 * fl1_fx * tg_zzzzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_y_1[j];

                    tg_xxzzzzzz_z_0[j] = pb_x * tg_xzzzzzz_z_0[j] + wp_x[j] * tg_xzzzzzz_z_1[j] + 0.5 * fl1_fx * tg_zzzzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_z_1[j];

                    tg_xyyyyyyy_x_0[j] = pb_x * tg_yyyyyyy_x_0[j] + wp_x[j] * tg_yyyyyyy_x_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_0_1[j];

                    tg_xyyyyyyy_y_0[j] = pb_x * tg_yyyyyyy_y_0[j] + wp_x[j] * tg_yyyyyyy_y_1[j];

                    tg_xyyyyyyy_z_0[j] = pb_x * tg_yyyyyyy_z_0[j] + wp_x[j] * tg_yyyyyyy_z_1[j];

                    tg_xyyyyyyz_x_0[j] = pb_x * tg_yyyyyyz_x_0[j] + wp_x[j] * tg_yyyyyyz_x_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_0_1[j];

                    tg_xyyyyyyz_y_0[j] = pb_x * tg_yyyyyyz_y_0[j] + wp_x[j] * tg_yyyyyyz_y_1[j];

                    tg_xyyyyyyz_z_0[j] = pb_x * tg_yyyyyyz_z_0[j] + wp_x[j] * tg_yyyyyyz_z_1[j];

                    tg_xyyyyyzz_x_0[j] = pb_x * tg_yyyyyzz_x_0[j] + wp_x[j] * tg_yyyyyzz_x_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_0_1[j];

                    tg_xyyyyyzz_y_0[j] = pb_x * tg_yyyyyzz_y_0[j] + wp_x[j] * tg_yyyyyzz_y_1[j];

                    tg_xyyyyyzz_z_0[j] = pb_x * tg_yyyyyzz_z_0[j] + wp_x[j] * tg_yyyyyzz_z_1[j];

                    tg_xyyyyzzz_x_0[j] = pb_x * tg_yyyyzzz_x_0[j] + wp_x[j] * tg_yyyyzzz_x_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_0_1[j];

                    tg_xyyyyzzz_y_0[j] = pb_x * tg_yyyyzzz_y_0[j] + wp_x[j] * tg_yyyyzzz_y_1[j];

                    tg_xyyyyzzz_z_0[j] = pb_x * tg_yyyyzzz_z_0[j] + wp_x[j] * tg_yyyyzzz_z_1[j];

                    tg_xyyyzzzz_x_0[j] = pb_x * tg_yyyzzzz_x_0[j] + wp_x[j] * tg_yyyzzzz_x_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_0_1[j];

                    tg_xyyyzzzz_y_0[j] = pb_x * tg_yyyzzzz_y_0[j] + wp_x[j] * tg_yyyzzzz_y_1[j];

                    tg_xyyyzzzz_z_0[j] = pb_x * tg_yyyzzzz_z_0[j] + wp_x[j] * tg_yyyzzzz_z_1[j];

                    tg_xyyzzzzz_x_0[j] = pb_x * tg_yyzzzzz_x_0[j] + wp_x[j] * tg_yyzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_0_1[j];

                    tg_xyyzzzzz_y_0[j] = pb_x * tg_yyzzzzz_y_0[j] + wp_x[j] * tg_yyzzzzz_y_1[j];

                    tg_xyyzzzzz_z_0[j] = pb_x * tg_yyzzzzz_z_0[j] + wp_x[j] * tg_yyzzzzz_z_1[j];

                    tg_xyzzzzzz_x_0[j] = pb_x * tg_yzzzzzz_x_0[j] + wp_x[j] * tg_yzzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_0_1[j];

                    tg_xyzzzzzz_y_0[j] = pb_x * tg_yzzzzzz_y_0[j] + wp_x[j] * tg_yzzzzzz_y_1[j];

                    tg_xyzzzzzz_z_0[j] = pb_x * tg_yzzzzzz_z_0[j] + wp_x[j] * tg_yzzzzzz_z_1[j];

                    tg_xzzzzzzz_x_0[j] = pb_x * tg_zzzzzzz_x_0[j] + wp_x[j] * tg_zzzzzzz_x_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_0_1[j];

                    tg_xzzzzzzz_y_0[j] = pb_x * tg_zzzzzzz_y_0[j] + wp_x[j] * tg_zzzzzzz_y_1[j];

                    tg_xzzzzzzz_z_0[j] = pb_x * tg_zzzzzzz_z_0[j] + wp_x[j] * tg_zzzzzzz_z_1[j];

                    tg_yyyyyyyy_x_0[j] = pb_y * tg_yyyyyyy_x_0[j] + wp_y[j] * tg_yyyyyyy_x_1[j] + 3.5 * fl1_fx * tg_yyyyyy_x_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_x_1[j];

                    tg_yyyyyyyy_y_0[j] = pb_y * tg_yyyyyyy_y_0[j] + wp_y[j] * tg_yyyyyyy_y_1[j] + 3.5 * fl1_fx * tg_yyyyyy_y_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_y_1[j] + 0.5 * fl1_fxn * tg_yyyyyyy_0_1[j];

                    tg_yyyyyyyy_z_0[j] = pb_y * tg_yyyyyyy_z_0[j] + wp_y[j] * tg_yyyyyyy_z_1[j] + 3.5 * fl1_fx * tg_yyyyyy_z_0[j] - 3.5 * fl1_fx * fl1_fza * tg_yyyyyy_z_1[j];

                    tg_yyyyyyyz_x_0[j] = pb_y * tg_yyyyyyz_x_0[j] + wp_y[j] * tg_yyyyyyz_x_1[j] + 3.0 * fl1_fx * tg_yyyyyz_x_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_x_1[j];

                    tg_yyyyyyyz_y_0[j] = pb_y * tg_yyyyyyz_y_0[j] + wp_y[j] * tg_yyyyyyz_y_1[j] + 3.0 * fl1_fx * tg_yyyyyz_y_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_y_1[j] + 0.5 * fl1_fxn * tg_yyyyyyz_0_1[j];

                    tg_yyyyyyyz_z_0[j] = pb_y * tg_yyyyyyz_z_0[j] + wp_y[j] * tg_yyyyyyz_z_1[j] + 3.0 * fl1_fx * tg_yyyyyz_z_0[j] - 3.0 * fl1_fx * fl1_fza * tg_yyyyyz_z_1[j];

                    tg_yyyyyyzz_x_0[j] = pb_y * tg_yyyyyzz_x_0[j] + wp_y[j] * tg_yyyyyzz_x_1[j] + 2.5 * fl1_fx * tg_yyyyzz_x_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_x_1[j];

                    tg_yyyyyyzz_y_0[j] = pb_y * tg_yyyyyzz_y_0[j] + wp_y[j] * tg_yyyyyzz_y_1[j] + 2.5 * fl1_fx * tg_yyyyzz_y_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_y_1[j] + 0.5 * fl1_fxn * tg_yyyyyzz_0_1[j];

                    tg_yyyyyyzz_z_0[j] = pb_y * tg_yyyyyzz_z_0[j] + wp_y[j] * tg_yyyyyzz_z_1[j] + 2.5 * fl1_fx * tg_yyyyzz_z_0[j] - 2.5 * fl1_fx * fl1_fza * tg_yyyyzz_z_1[j];

                    tg_yyyyyzzz_x_0[j] = pb_y * tg_yyyyzzz_x_0[j] + wp_y[j] * tg_yyyyzzz_x_1[j] + 2.0 * fl1_fx * tg_yyyzzz_x_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_x_1[j];

                    tg_yyyyyzzz_y_0[j] = pb_y * tg_yyyyzzz_y_0[j] + wp_y[j] * tg_yyyyzzz_y_1[j] + 2.0 * fl1_fx * tg_yyyzzz_y_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_y_1[j] + 0.5 * fl1_fxn * tg_yyyyzzz_0_1[j];

                    tg_yyyyyzzz_z_0[j] = pb_y * tg_yyyyzzz_z_0[j] + wp_y[j] * tg_yyyyzzz_z_1[j] + 2.0 * fl1_fx * tg_yyyzzz_z_0[j] - 2.0 * fl1_fx * fl1_fza * tg_yyyzzz_z_1[j];

                    tg_yyyyzzzz_x_0[j] = pb_y * tg_yyyzzzz_x_0[j] + wp_y[j] * tg_yyyzzzz_x_1[j] + 1.5 * fl1_fx * tg_yyzzzz_x_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_x_1[j];

                    tg_yyyyzzzz_y_0[j] = pb_y * tg_yyyzzzz_y_0[j] + wp_y[j] * tg_yyyzzzz_y_1[j] + 1.5 * fl1_fx * tg_yyzzzz_y_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_y_1[j] + 0.5 * fl1_fxn * tg_yyyzzzz_0_1[j];

                    tg_yyyyzzzz_z_0[j] = pb_y * tg_yyyzzzz_z_0[j] + wp_y[j] * tg_yyyzzzz_z_1[j] + 1.5 * fl1_fx * tg_yyzzzz_z_0[j] - 1.5 * fl1_fx * fl1_fza * tg_yyzzzz_z_1[j];

                    tg_yyyzzzzz_x_0[j] = pb_y * tg_yyzzzzz_x_0[j] + wp_y[j] * tg_yyzzzzz_x_1[j] + fl1_fx * tg_yzzzzz_x_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_x_1[j];

                    tg_yyyzzzzz_y_0[j] = pb_y * tg_yyzzzzz_y_0[j] + wp_y[j] * tg_yyzzzzz_y_1[j] + fl1_fx * tg_yzzzzz_y_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_y_1[j] + 0.5 * fl1_fxn * tg_yyzzzzz_0_1[j];

                    tg_yyyzzzzz_z_0[j] = pb_y * tg_yyzzzzz_z_0[j] + wp_y[j] * tg_yyzzzzz_z_1[j] + fl1_fx * tg_yzzzzz_z_0[j] - fl1_fx * fl1_fza * tg_yzzzzz_z_1[j];

                    tg_yyzzzzzz_x_0[j] = pb_y * tg_yzzzzzz_x_0[j] + wp_y[j] * tg_yzzzzzz_x_1[j] + 0.5 * fl1_fx * tg_zzzzzz_x_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_x_1[j];

                    tg_yyzzzzzz_y_0[j] = pb_y * tg_yzzzzzz_y_0[j] + wp_y[j] * tg_yzzzzzz_y_1[j] + 0.5 * fl1_fx * tg_zzzzzz_y_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_y_1[j] + 0.5 * fl1_fxn * tg_yzzzzzz_0_1[j];

                    tg_yyzzzzzz_z_0[j] = pb_y * tg_yzzzzzz_z_0[j] + wp_y[j] * tg_yzzzzzz_z_1[j] + 0.5 * fl1_fx * tg_zzzzzz_z_0[j] - 0.5 * fl1_fx * fl1_fza * tg_zzzzzz_z_1[j];

                    tg_yzzzzzzz_x_0[j] = pb_y * tg_zzzzzzz_x_0[j] + wp_y[j] * tg_zzzzzzz_x_1[j];

                    tg_yzzzzzzz_y_0[j] = pb_y * tg_zzzzzzz_y_0[j] + wp_y[j] * tg_zzzzzzz_y_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_0_1[j];

                    tg_yzzzzzzz_z_0[j] = pb_y * tg_zzzzzzz_z_0[j] + wp_y[j] * tg_zzzzzzz_z_1[j];

                    tg_zzzzzzzz_x_0[j] = pb_z * tg_zzzzzzz_x_0[j] + wp_z[j] * tg_zzzzzzz_x_1[j] + 3.5 * fl1_fx * tg_zzzzzz_x_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_x_1[j];

                    tg_zzzzzzzz_y_0[j] = pb_z * tg_zzzzzzz_y_0[j] + wp_z[j] * tg_zzzzzzz_y_1[j] + 3.5 * fl1_fx * tg_zzzzzz_y_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_y_1[j];

                    tg_zzzzzzzz_z_0[j] = pb_z * tg_zzzzzzz_z_0[j] + wp_z[j] * tg_zzzzzzz_z_1[j] + 3.5 * fl1_fx * tg_zzzzzz_z_0[j] - 3.5 * fl1_fx * fl1_fza * tg_zzzzzz_z_1[j] + 0.5 * fl1_fxn * tg_zzzzzzz_0_1[j];
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

