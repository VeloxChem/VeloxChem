//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectronRepulsionRecFuncForSX.hpp"

#include "MathConst.hpp"

namespace erirecfunc { // erirecfunc namespace

    void
    compElectronRepulsionForSSSS(      CMemBlock2D<double>*  primBuffer,
                                 const CRecursionMap&        recursionMap,
                                 const CBoysFunction&        bfTable,
                                       CMemBlock<double>&    bfArguments,
                                       CMemBlock2D<double>&  bfValues,
                                 const int32_t               bfOrder,
                                 const CMemBlock2D<double>&  osFactors,
                                 const CMemBlock2D<double>&  pqDistances,
                                 const CGtoPairsBlock&       braGtoPairsBlock,
                                 const CGtoPairsBlock&       ketGtoPairsBlock,
                                 const int32_t               nKetPrimPairs,
                                 const int32_t               iContrPair)
    {
        // set up pointers to primitive pairs data on bra side
        
        auto bfss = braGtoPairsBlock.getOverlaps();
        
        auto spos = braGtoPairsBlock.getStartPositions();
        
        auto epos = braGtoPairsBlock.getEndPositions();
        
        // set up pointers to primitives data on ket side
        
        auto kfss = ketGtoPairsBlock.getOverlaps();
        
        // set up pi prefactor
        
        auto fpi = 2.0 / std::sqrt(mathconst::getPiValue());
        
        // loop over contracted GTO on bra side
        
        int32_t idx = 0;
        
        for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
        {
            // set up Obara-Saika prefactors
            
            auto fz = osFactors.data(4 * idx + 1);

            // set up R(PQ) distances
            
            auto pqx = pqDistances.data(3 * idx);
            
            auto pqy = pqDistances.data(3 * idx + 1);
            
            auto pqz = pqDistances.data(3 * idx + 2);
            
            // compute Boys function argument
            
            auto fargs = bfArguments.data();
            
            #pragma omp simd aligned(fargs, fz, pqx, pqy, pqz: VLX_ALIGN)
            for (int32_t j = 0; j < nKetPrimPairs; j++)
            {
                fargs[j] = fz[j] * (pqx[j] * pqx[j] + pqy[j] * pqy[j] +
                                    
                                    pqz[j] * pqz[j]);
            }
            
            // evaluate Boys function values
            
            bfTable.compute(bfValues, bfArguments, nKetPrimPairs, bfOrder);
            
            // set up pointers to Obara-Saika factors
            
            auto fss = bfss[i];
            
            // compute overlap scaling factor
            
            #pragma omp simd aligned(fz, kfss, fargs: VLX_ALIGN)
            for (int32_t j = 0; j < nKetPrimPairs; j++)
            {
                fargs[j] = fss * kfss[j] * fpi * std::sqrt(fz[j]);
            }
            
            // distribute (SS|g(r,r')|SS) integrals
            
            for (int32_t j = 0; j <= bfOrder; j++)
            {
                auto pidx = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                              {0, -1, -1, -1}, {0, -1, -1, -1},
                                                               1, 1, j));
                
                if (pidx != -1)
                {
                    auto g_00_00 = primBuffer[pidx].data(idx);
                    
                    auto bvals = bfValues.data(j);
                    
                    #pragma omp simd aligned(g_00_00, bvals, fargs: VLX_ALIGN)
                    for (int32_t k = 0; k < nKetPrimPairs; k++)
                    {
                        g_00_00[k] = bvals[k] * fargs[k];
                    }
                }
            }
            
            idx++;
        }
    }

    void
    compElectronRepulsionForSSSP(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {1, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_1_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {1, -1, -1, -1},
                                                                    1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_1_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_0_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_0_0 = primBuffer[pidx_g_0_0_m0].data(idx);

                auto tg_0_0_1 = primBuffer[pidx_g_0_0_m1].data(idx);

                // set up pointers to integrals

                auto tg_0_x_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx);

                auto tg_0_y_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx + 1);

                auto tg_0_z_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx + 2);

                #pragma omp simd aligned(qd_x, qd_y, qd_z, tg_0_0_0, tg_0_0_1, tg_0_x_0, tg_0_y_0, tg_0_z_0, wq_x, wq_y, \
                                         wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double tgval = tg_0_0_1[j];
                    
                    tg_0_x_0[j] = qd_x[j] * tg_0_0_0[j] + wq_x[j] * tgval;

                    tg_0_y_0[j] = qd_y[j] * tg_0_0_0[j] + wq_y[j] * tgval;

                    tg_0_z_0[j] = qd_z[j] * tg_0_0_0[j] + wq_z[j] * tgval;
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSPSS(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
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

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {1, -1, -1, -1},
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_1_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {1, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_1_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_0_0 = primBuffer[pidx_g_0_0_m0].data(idx);

                auto tg_0_0_1 = primBuffer[pidx_g_0_0_m1].data(idx);

                // set up pointers to integrals

                auto tg_x_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx);

                auto tg_y_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx + 1);

                auto tg_z_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx + 2);

                #pragma omp simd aligned(tg_0_0_0, tg_0_0_1, tg_x_0_0, tg_y_0_0, tg_z_0_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double tgval = tg_0_0_1[j];
                    
                    tg_x_0_0[j] = pb_x * tg_0_0_0[j] + wp_x[j] * tgval;

                    tg_y_0_0[j] = pb_y * tg_0_0_0[j] + wp_y[j] * tgval;

                    tg_z_0_0[j] = pb_z * tg_0_0_0[j] + wp_z[j] * tgval;
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSSSD(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto fn = ketGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_2_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {2, -1, -1, -1},
                                                                    1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_1_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {1, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_0_1_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {1, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fzb = osFactors.data(4 * idx + 3);

                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_x_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx);

                auto tg_0_y_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx + 1);

                auto tg_0_z_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx + 2);

                auto tg_0_x_1 = primBuffer[pidx_g_0_1_m1].data(3 * idx);

                auto tg_0_y_1 = primBuffer[pidx_g_0_1_m1].data(3 * idx + 1);

                auto tg_0_z_1 = primBuffer[pidx_g_0_1_m1].data(3 * idx + 2);

                auto tg_0_0_0 = primBuffer[pidx_g_0_0_m0].data(idx);

                auto tg_0_0_1 = primBuffer[pidx_g_0_0_m1].data(idx);

                // set up pointers to integrals

                auto tg_0_xx_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx);

                auto tg_0_xy_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 1);

                auto tg_0_xz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 2);

                auto tg_0_yy_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 3);

                auto tg_0_yz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 4);

                auto tg_0_zz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 5);

                #pragma omp simd aligned(fn, fzb, qd_x, qd_y, qd_z, tg_0_0_0, tg_0_0_1, tg_0_x_0, tg_0_x_1, tg_0_xx_0, \
                                         tg_0_xy_0, tg_0_xz_0, tg_0_y_0, tg_0_y_1, tg_0_yy_0, tg_0_yz_0, tg_0_z_0, tg_0_z_1, \
                                         tg_0_zz_0, wq_x, wq_y, wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double tgval =  0.5 * fn[j] * (tg_0_0_0[j] - fzb[j] * tg_0_0_1[j]);

                    double fr0 = qd_x[j]; double fr1 = wq_x[j];
                    
                    tg_0_xx_0[j] = fr0 * tg_0_x_0[j] + fr1 * tg_0_x_1[j] + tgval;

                    tg_0_xy_0[j] = fr0 * tg_0_y_0[j] + fr1 * tg_0_y_1[j];

                    tg_0_xz_0[j] = fr0 * tg_0_z_0[j] + fr1 * tg_0_z_1[j];
                    
                    fr0 = qd_y[j]; fr1 = wq_y[j];

                    tg_0_yy_0[j] = fr0 * tg_0_y_0[j] + fr1 * tg_0_y_1[j] + tgval;

                    tg_0_yz_0[j] = fr0 * tg_0_z_0[j] + fr1 * tg_0_z_1[j];

                    tg_0_zz_0[j] = qd_z[j] * tg_0_z_0[j] + wq_z[j] * tg_0_z_1[j] + tgval;
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSDSS(      CMemBlock2D<double>* primBuffer,
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
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_2_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                    {2, -1, -1, -1}, {0, -1, -1, -1},
                                                                     1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_2_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_1_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {1, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_1_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {1, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_x_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx);

                auto tg_y_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx + 1);

                auto tg_z_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx + 2);

                auto tg_x_0_1 = primBuffer[pidx_g_1_0_m1].data(3 * idx);

                auto tg_y_0_1 = primBuffer[pidx_g_1_0_m1].data(3 * idx + 1);

                auto tg_z_0_1 = primBuffer[pidx_g_1_0_m1].data(3 * idx + 2);

                auto tg_0_0_0 = primBuffer[pidx_g_0_0_m0].data(idx);

                auto tg_0_0_1 = primBuffer[pidx_g_0_0_m1].data(idx);

                // set up pointers to integrals

                auto tg_xx_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx);

                auto tg_xy_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 1);

                auto tg_xz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 2);

                auto tg_yy_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 3);

                auto tg_yz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 4);

                auto tg_zz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 5);

                #pragma omp simd aligned(fza, tg_0_0_0, tg_0_0_1, tg_x_0_0, tg_x_0_1, tg_xx_0_0, tg_xy_0_0, \
                                         tg_xz_0_0, tg_y_0_0, tg_y_0_1, tg_yy_0_0, tg_yz_0_0, tg_z_0_0, tg_z_0_1, tg_zz_0_0, \
                                         wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double tgval =  0.5 * fx * (tg_0_0_0[j] - fza[j] * tg_0_0_1[j]);
                    
                    double fr1 = wp_x[j];

                    tg_xx_0_0[j] = pb_x * tg_x_0_0[j] + fr1 * tg_x_0_1[j] + tgval;

                    tg_xy_0_0[j] = pb_x * tg_y_0_0[j] + fr1 * tg_y_0_1[j];

                    tg_xz_0_0[j] = pb_x * tg_z_0_0[j] + fr1 * tg_z_0_1[j];
                    
                    fr1 = wp_y[j];

                    tg_yy_0_0[j] = pb_y * tg_y_0_0[j] + fr1 * tg_y_0_1[j] + tgval;

                    tg_yz_0_0[j] = pb_y * tg_z_0_0[j] + fr1 * tg_z_0_1[j];

                    tg_zz_0_0[j] = pb_z * tg_z_0_0[j] + wp_z[j] * tg_z_0_1[j] + tgval;
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSSSF(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto fn = ketGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {3, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_2_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {2, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {2, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_1_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {1, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_1_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {1, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fzb = osFactors.data(4 * idx + 3);

                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xx_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx);

                auto tg_0_xy_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 1);

                auto tg_0_xz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 2);

                auto tg_0_yy_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 3);

                auto tg_0_yz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 4);

                auto tg_0_zz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 5);

                auto tg_0_xx_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx);

                auto tg_0_xy_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 1);

                auto tg_0_xz_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 2);

                auto tg_0_yy_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 3);

                auto tg_0_yz_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 4);

                auto tg_0_zz_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 5);

                auto tg_0_x_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx);

                auto tg_0_y_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx + 1);

                auto tg_0_z_0 = primBuffer[pidx_g_0_1_m0].data(3 * idx + 2);

                auto tg_0_x_1 = primBuffer[pidx_g_0_1_m1].data(3 * idx);

                auto tg_0_y_1 = primBuffer[pidx_g_0_1_m1].data(3 * idx + 1);

                auto tg_0_z_1 = primBuffer[pidx_g_0_1_m1].data(3 * idx + 2);

                // set up pointers to integrals

                auto tg_0_xxx_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx);

                auto tg_0_xxy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 1);

                auto tg_0_xxz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 2);

                auto tg_0_xyy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 3);

                auto tg_0_xyz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 4);

                auto tg_0_xzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 5);

                auto tg_0_yyy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 6);

                auto tg_0_yyz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 7);

                auto tg_0_yzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 8);

                auto tg_0_zzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 9);

                #pragma omp simd aligned(fn, fzb, qd_x, qd_y, qd_z, tg_0_x_0, tg_0_x_1, tg_0_xx_0, tg_0_xx_1, \
                                         tg_0_xxx_0, tg_0_xxy_0, tg_0_xxz_0, tg_0_xy_0, tg_0_xy_1, tg_0_xyy_0, tg_0_xyz_0, \
                                         tg_0_xz_0, tg_0_xz_1, tg_0_xzz_0, tg_0_y_0, tg_0_y_1, tg_0_yy_0, tg_0_yy_1, \
                                         tg_0_yyy_0, tg_0_yyz_0, tg_0_yz_0, tg_0_yz_1, tg_0_yzz_0, tg_0_z_0, tg_0_z_1, \
                                         tg_0_zz_0, tg_0_zz_1, tg_0_zzz_0, wq_x, wq_y, wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fn = fn[j];

                    double fl1_fzb = fzb[j];
                    
                    double fr0 = qd_x[j]; double fr1 = wq_x[j];

                    tg_0_xxx_0[j] = fr0 * tg_0_xx_0[j] + fr1 * tg_0_xx_1[j] + fl1_fn * (tg_0_x_0[j] - fl1_fzb * tg_0_x_1[j]);

                    tg_0_xxy_0[j] = fr0 * tg_0_xy_0[j] + fr1 * tg_0_xy_1[j] + 0.5 * fl1_fn * (tg_0_y_0[j] - fl1_fzb * tg_0_y_1[j]);

                    tg_0_xxz_0[j] = fr0 * tg_0_xz_0[j] + fr1 * tg_0_xz_1[j] + 0.5 * fl1_fn * (tg_0_z_0[j] - fl1_fzb * tg_0_z_1[j]);

                    tg_0_xyy_0[j] = fr0 * tg_0_yy_0[j] + fr1 * tg_0_yy_1[j];

                    tg_0_xyz_0[j] = fr0 * tg_0_yz_0[j] + fr1 * tg_0_yz_1[j];

                    tg_0_xzz_0[j] = fr0 * tg_0_zz_0[j] + fr1 * tg_0_zz_1[j];
                    
                    fr0 = qd_y[j]; fr1 = wq_y[j];

                    tg_0_yyy_0[j] = fr0 * tg_0_yy_0[j] + fr1 * tg_0_yy_1[j] + fl1_fn * (tg_0_y_0[j] - fl1_fzb * tg_0_y_1[j]);

                    tg_0_yyz_0[j] = fr0 * tg_0_yz_0[j] + fr1 * tg_0_yz_1[j] + 0.5 * fl1_fn * (tg_0_z_0[j] - fl1_fzb * tg_0_z_1[j]);

                    tg_0_yzz_0[j] = fr0 * tg_0_zz_0[j] + fr1 * tg_0_zz_1[j];

                    tg_0_zzz_0[j] = qd_z[j] * tg_0_zz_0[j] + wq_z[j] * tg_0_zz_1[j] + fl1_fn * (tg_0_z_0[j] - fl1_fzb * tg_0_z_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSFSS(      CMemBlock2D<double>* primBuffer,
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
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_3_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {3, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_3_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_2_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {2, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_2_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {2, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord + 1));

            auto pidx_g_1_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {1, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_1_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {1, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_xx_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx);

                auto tg_xy_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 1);

                auto tg_xz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 2);

                auto tg_yy_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 3);

                auto tg_yz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 4);

                auto tg_zz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 5);

                auto tg_xx_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx);

                auto tg_xy_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 1);

                auto tg_xz_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 2);

                auto tg_yy_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 3);

                auto tg_yz_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 4);

                auto tg_zz_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 5);

                auto tg_x_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx);

                auto tg_y_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx + 1);

                auto tg_z_0_0 = primBuffer[pidx_g_1_0_m0].data(3 * idx + 2);

                auto tg_x_0_1 = primBuffer[pidx_g_1_0_m1].data(3 * idx);

                auto tg_y_0_1 = primBuffer[pidx_g_1_0_m1].data(3 * idx + 1);

                auto tg_z_0_1 = primBuffer[pidx_g_1_0_m1].data(3 * idx + 2);

                // set up pointers to integrals

                auto tg_xxx_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx);

                auto tg_xxy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 1);

                auto tg_xxz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 2);

                auto tg_xyy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 3);

                auto tg_xyz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 4);

                auto tg_xzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 5);

                auto tg_yyy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 6);

                auto tg_yyz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 7);

                auto tg_yzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 8);

                auto tg_zzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 9);

                #pragma omp simd aligned(fza, tg_x_0_0, tg_x_0_1, tg_xx_0_0, tg_xx_0_1, tg_xxx_0_0, tg_xxy_0_0, \
                                         tg_xxz_0_0, tg_xy_0_0, tg_xy_0_1, tg_xyy_0_0, tg_xyz_0_0, tg_xz_0_0, tg_xz_0_1, \
                                         tg_xzz_0_0, tg_y_0_0, tg_y_0_1, tg_yy_0_0, tg_yy_0_1, tg_yyy_0_0, tg_yyz_0_0, \
                                         tg_yz_0_0, tg_yz_0_1, tg_yzz_0_0, tg_z_0_0, tg_z_0_1, tg_zz_0_0, tg_zz_0_1, \
                                         tg_zzz_0_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fza = fza[j];
                    
                    double fr1 = wp_x[j];

                    tg_xxx_0_0[j] = pb_x * tg_xx_0_0[j] + fr1 * tg_xx_0_1[j] + fl1_fx * (tg_x_0_0[j] - fl1_fza * tg_x_0_1[j]);

                    tg_xxy_0_0[j] = pb_x * tg_xy_0_0[j] + fr1 * tg_xy_0_1[j] + 0.5 * fl1_fx * (tg_y_0_0[j] - fl1_fza * tg_y_0_1[j]);

                    tg_xxz_0_0[j] = pb_x * tg_xz_0_0[j] + fr1 * tg_xz_0_1[j] + 0.5 * fl1_fx * (tg_z_0_0[j] - fl1_fza * tg_z_0_1[j]);

                    tg_xyy_0_0[j] = pb_x * tg_yy_0_0[j] + fr1 * tg_yy_0_1[j];

                    tg_xyz_0_0[j] = pb_x * tg_yz_0_0[j] + fr1 * tg_yz_0_1[j];

                    tg_xzz_0_0[j] = pb_x * tg_zz_0_0[j] + fr1 * tg_zz_0_1[j];
                    
                    fr1 = wp_y[j];

                    tg_yyy_0_0[j] = pb_y * tg_yy_0_0[j] + fr1 * tg_yy_0_1[j] + fl1_fx * (tg_y_0_0[j] - fl1_fza * tg_y_0_1[j]);

                    tg_yyz_0_0[j] = pb_y * tg_yz_0_0[j] + fr1 * tg_yz_0_1[j] + 0.5 * fl1_fx * (tg_z_0_0[j] - fl1_fza * tg_z_0_1[j]);

                    tg_yzz_0_0[j] = pb_y * tg_zz_0_0[j] + fr1 * tg_zz_0_1[j];

                    tg_zzz_0_0[j] = pb_z * tg_zz_0_0[j] + wp_z[j] * tg_zz_0_1[j] + fl1_fx * (tg_z_0_0[j] - fl1_fza * tg_z_0_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSSSG(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto fn = ketGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {4, -1, -1, -1},
                                                                    1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {3, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_0_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {3, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_2_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {2, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_0_2_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {2, -1, -1, -1},
                                                                    1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fzb = osFactors.data(4 * idx + 3);

                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xxx_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx);

                auto tg_0_xxy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 1);

                auto tg_0_xxz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 2);

                auto tg_0_xyy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 3);

                auto tg_0_xyz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 4);

                auto tg_0_xzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 5);

                auto tg_0_yyy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 6);

                auto tg_0_yyz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 7);

                auto tg_0_yzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 8);

                auto tg_0_zzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 9);

                auto tg_0_xxx_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx);

                auto tg_0_xxy_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 1);

                auto tg_0_xxz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 2);

                auto tg_0_xyy_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 3);

                auto tg_0_xyz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 4);

                auto tg_0_xzz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 5);

                auto tg_0_yyy_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 6);

                auto tg_0_yyz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 7);

                auto tg_0_yzz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 8);

                auto tg_0_zzz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 9);

                auto tg_0_xx_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx);

                auto tg_0_xy_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 1);

                auto tg_0_xz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 2);

                auto tg_0_yy_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 3);

                auto tg_0_yz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 4);

                auto tg_0_zz_0 = primBuffer[pidx_g_0_2_m0].data(6 * idx + 5);

                auto tg_0_xx_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx);

                auto tg_0_xy_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 1);

                auto tg_0_xz_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 2);

                auto tg_0_yy_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 3);

                auto tg_0_yz_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 4);

                auto tg_0_zz_1 = primBuffer[pidx_g_0_2_m1].data(6 * idx + 5);

                // set up pointers to integrals

                auto tg_0_xxxx_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx);

                auto tg_0_xxxy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 1);

                auto tg_0_xxxz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 2);

                auto tg_0_xxyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 3);

                auto tg_0_xxyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 4);

                auto tg_0_xxzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 5);

                auto tg_0_xyyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 6);

                auto tg_0_xyyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 7);

                auto tg_0_xyzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 8);

                auto tg_0_xzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 9);

                auto tg_0_yyyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 10);

                auto tg_0_yyyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 11);

                auto tg_0_yyzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 12);

                auto tg_0_yzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 13);

                auto tg_0_zzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 14);

                #pragma omp simd aligned(fn, fzb, qd_x, qd_y, qd_z, tg_0_xx_0, tg_0_xx_1, tg_0_xxx_0, tg_0_xxx_1, \
                                         tg_0_xxxx_0, tg_0_xxxy_0, tg_0_xxxz_0, tg_0_xxy_0, tg_0_xxy_1, tg_0_xxyy_0, \
                                         tg_0_xxyz_0, tg_0_xxz_0, tg_0_xxz_1, tg_0_xxzz_0, tg_0_xy_0, tg_0_xy_1, tg_0_xyy_0, \
                                         tg_0_xyy_1, tg_0_xyyy_0, tg_0_xyyz_0, tg_0_xyz_0, tg_0_xyz_1, tg_0_xyzz_0, \
                                         tg_0_xz_0, tg_0_xz_1, tg_0_xzz_0, tg_0_xzz_1, tg_0_xzzz_0, tg_0_yy_0, tg_0_yy_1, \
                                         tg_0_yyy_0, tg_0_yyy_1, tg_0_yyyy_0, tg_0_yyyz_0, tg_0_yyz_0, tg_0_yyz_1, \
                                         tg_0_yyzz_0, tg_0_yz_0, tg_0_yz_1, tg_0_yzz_0, tg_0_yzz_1, tg_0_yzzz_0, tg_0_zz_0, \
                                         tg_0_zz_1, tg_0_zzz_0, tg_0_zzz_1, tg_0_zzzz_0, wq_x, wq_y, wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fn = fn[j];

                    double fl1_fzb = fzb[j];
                    
                    double fr0 = qd_x[j]; double fr1 = wq_x[j];

                    tg_0_xxxx_0[j] = fr0 * tg_0_xxx_0[j] + fr1 * tg_0_xxx_1[j] + 1.5 * fl1_fn * (tg_0_xx_0[j] - fl1_fzb * tg_0_xx_1[j]);

                    tg_0_xxxy_0[j] = fr0 * tg_0_xxy_0[j] + fr1 * tg_0_xxy_1[j] + fl1_fn * (tg_0_xy_0[j] - fl1_fzb * tg_0_xy_1[j]);

                    tg_0_xxxz_0[j] = fr0 * tg_0_xxz_0[j] + fr1 * tg_0_xxz_1[j] + fl1_fn * (tg_0_xz_0[j] - fl1_fzb * tg_0_xz_1[j]);

                    tg_0_xxyy_0[j] = fr0 * tg_0_xyy_0[j] + fr1 * tg_0_xyy_1[j] + 0.5 * fl1_fn * (tg_0_yy_0[j] - fl1_fzb * tg_0_yy_1[j]);

                    tg_0_xxyz_0[j] = fr0 * tg_0_xyz_0[j] + fr1 * tg_0_xyz_1[j] + 0.5 * fl1_fn * (tg_0_yz_0[j] - fl1_fzb * tg_0_yz_1[j]);

                    tg_0_xxzz_0[j] = fr0 * tg_0_xzz_0[j] + fr1 * tg_0_xzz_1[j] + 0.5 * fl1_fn * (tg_0_zz_0[j] - fl1_fzb * tg_0_zz_1[j]);

                    tg_0_xyyy_0[j] = fr0 * tg_0_yyy_0[j] + fr1 * tg_0_yyy_1[j];

                    tg_0_xyyz_0[j] = fr0 * tg_0_yyz_0[j] + fr1 * tg_0_yyz_1[j];

                    tg_0_xyzz_0[j] = fr0 * tg_0_yzz_0[j] + fr1 * tg_0_yzz_1[j];

                    tg_0_xzzz_0[j] = fr0 * tg_0_zzz_0[j] + fr1 * tg_0_zzz_1[j];
                    
                    fr0 = qd_y[j]; fr1 = wq_y[j];

                    tg_0_yyyy_0[j] = fr0 * tg_0_yyy_0[j] + fr1 * tg_0_yyy_1[j] + 1.5 * fl1_fn * (tg_0_yy_0[j] - fl1_fzb * tg_0_yy_1[j]);

                    tg_0_yyyz_0[j] = fr0 * tg_0_yyz_0[j] + fr1 * tg_0_yyz_1[j] + fl1_fn * (tg_0_yz_0[j] - fl1_fzb * tg_0_yz_1[j]);

                    tg_0_yyzz_0[j] = fr0 * tg_0_yzz_0[j] + fr1 * tg_0_yzz_1[j] + 0.5 * fl1_fn * (tg_0_zz_0[j] - fl1_fzb * tg_0_zz_1[j]);

                    tg_0_yzzz_0[j] = fr0 * tg_0_zzz_0[j] + fr1 * tg_0_zzz_1[j];

                    tg_0_zzzz_0[j] = qd_z[j] * tg_0_zzz_0[j] + wq_z[j] * tg_0_zzz_1[j] + 1.5 * fl1_fn * (tg_0_zz_0[j] - fl1_fzb * tg_0_zz_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSGSS(      CMemBlock2D<double>* primBuffer,
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
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_4_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {4, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_4_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_3_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {3, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_3_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {3, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord + 1));

            auto pidx_g_2_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {2, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_2_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {2, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_xxx_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx);

                auto tg_xxy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 1);

                auto tg_xxz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 2);

                auto tg_xyy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 3);

                auto tg_xyz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 4);

                auto tg_xzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 5);

                auto tg_yyy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 6);

                auto tg_yyz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 7);

                auto tg_yzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 8);

                auto tg_zzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 9);

                auto tg_xxx_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx);

                auto tg_xxy_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 1);

                auto tg_xxz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 2);

                auto tg_xyy_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 3);

                auto tg_xyz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 4);

                auto tg_xzz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 5);

                auto tg_yyy_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 6);

                auto tg_yyz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 7);

                auto tg_yzz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 8);

                auto tg_zzz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 9);

                auto tg_xx_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx);

                auto tg_xy_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 1);

                auto tg_xz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 2);

                auto tg_yy_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 3);

                auto tg_yz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 4);

                auto tg_zz_0_0 = primBuffer[pidx_g_2_0_m0].data(6 * idx + 5);

                auto tg_xx_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx);

                auto tg_xy_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 1);

                auto tg_xz_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 2);

                auto tg_yy_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 3);

                auto tg_yz_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 4);

                auto tg_zz_0_1 = primBuffer[pidx_g_2_0_m1].data(6 * idx + 5);

                // set up pointers to integrals

                auto tg_xxxx_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx);

                auto tg_xxxy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 1);

                auto tg_xxxz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 2);

                auto tg_xxyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 3);

                auto tg_xxyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 4);

                auto tg_xxzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 5);

                auto tg_xyyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 6);

                auto tg_xyyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 7);

                auto tg_xyzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 8);

                auto tg_xzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 9);

                auto tg_yyyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 10);

                auto tg_yyyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 11);

                auto tg_yyzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 12);

                auto tg_yzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 13);

                auto tg_zzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 14);

                #pragma omp simd aligned(fza, tg_xx_0_0, tg_xx_0_1, tg_xxx_0_0, tg_xxx_0_1, tg_xxxx_0_0, \
                                         tg_xxxy_0_0, tg_xxxz_0_0, tg_xxy_0_0, tg_xxy_0_1, tg_xxyy_0_0, tg_xxyz_0_0, \
                                         tg_xxz_0_0, tg_xxz_0_1, tg_xxzz_0_0, tg_xy_0_0, tg_xy_0_1, tg_xyy_0_0, tg_xyy_0_1, \
                                         tg_xyyy_0_0, tg_xyyz_0_0, tg_xyz_0_0, tg_xyz_0_1, tg_xyzz_0_0, tg_xz_0_0, tg_xz_0_1, \
                                         tg_xzz_0_0, tg_xzz_0_1, tg_xzzz_0_0, tg_yy_0_0, tg_yy_0_1, tg_yyy_0_0, tg_yyy_0_1, \
                                         tg_yyyy_0_0, tg_yyyz_0_0, tg_yyz_0_0, tg_yyz_0_1, tg_yyzz_0_0, tg_yz_0_0, tg_yz_0_1, \
                                         tg_yzz_0_0, tg_yzz_0_1, tg_yzzz_0_0, tg_zz_0_0, tg_zz_0_1, tg_zzz_0_0, tg_zzz_0_1, \
                                         tg_zzzz_0_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fza = fza[j];
                    
                    double fr1 = wp_x[j];

                    tg_xxxx_0_0[j] = pb_x * tg_xxx_0_0[j] + fr1 * tg_xxx_0_1[j] + 1.5 * fl1_fx * (tg_xx_0_0[j] - fl1_fza * tg_xx_0_1[j]);

                    tg_xxxy_0_0[j] = pb_x * tg_xxy_0_0[j] + fr1 * tg_xxy_0_1[j] + fl1_fx * (tg_xy_0_0[j] - fl1_fza * tg_xy_0_1[j]);

                    tg_xxxz_0_0[j] = pb_x * tg_xxz_0_0[j] + fr1 * tg_xxz_0_1[j] + fl1_fx * (tg_xz_0_0[j] - fl1_fza * tg_xz_0_1[j]);

                    tg_xxyy_0_0[j] = pb_x * tg_xyy_0_0[j] + fr1 * tg_xyy_0_1[j] + 0.5 * fl1_fx * (tg_yy_0_0[j] - fl1_fza * tg_yy_0_1[j]);

                    tg_xxyz_0_0[j] = pb_x * tg_xyz_0_0[j] + fr1 * tg_xyz_0_1[j] + 0.5 * fl1_fx * (tg_yz_0_0[j] - fl1_fza * tg_yz_0_1[j]);

                    tg_xxzz_0_0[j] = pb_x * tg_xzz_0_0[j] + fr1 * tg_xzz_0_1[j] + 0.5 * fl1_fx * (tg_zz_0_0[j] - fl1_fza * tg_zz_0_1[j]);

                    tg_xyyy_0_0[j] = pb_x * tg_yyy_0_0[j] + fr1 * tg_yyy_0_1[j];

                    tg_xyyz_0_0[j] = pb_x * tg_yyz_0_0[j] + fr1 * tg_yyz_0_1[j];

                    tg_xyzz_0_0[j] = pb_x * tg_yzz_0_0[j] + fr1 * tg_yzz_0_1[j];

                    tg_xzzz_0_0[j] = pb_x * tg_zzz_0_0[j] + fr1 * tg_zzz_0_1[j];
                    
                    fr1 = wp_y[j];

                    tg_yyyy_0_0[j] = pb_y * tg_yyy_0_0[j] + fr1 * tg_yyy_0_1[j] + 1.5 * fl1_fx * (tg_yy_0_0[j] - fl1_fza * tg_yy_0_1[j]);

                    tg_yyyz_0_0[j] = pb_y * tg_yyz_0_0[j] + fr1 * tg_yyz_0_1[j] + fl1_fx * (tg_yz_0_0[j] - fl1_fza * tg_yz_0_1[j]);

                    tg_yyzz_0_0[j] = pb_y * tg_yzz_0_0[j] + fr1 * tg_yzz_0_1[j] + 0.5 * fl1_fx * (tg_zz_0_0[j] - fl1_fza * tg_zz_0_1[j]);

                    tg_yzzz_0_0[j] = pb_y * tg_zzz_0_0[j] + fr1 * tg_zzz_0_1[j];

                    tg_zzzz_0_0[j] = pb_z * tg_zzz_0_0[j] + wp_z[j] * tg_zzz_0_1[j] + 1.5 * fl1_fx * (tg_zz_0_0[j] - fl1_fza * tg_zz_0_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSSSH(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto fn = ketGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {5, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {5, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_5_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {4, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {4, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_3_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {3, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_3_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {3, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fzb = osFactors.data(4 * idx + 3);

                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xxxx_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx);

                auto tg_0_xxxy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 1);

                auto tg_0_xxxz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 2);

                auto tg_0_xxyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 3);

                auto tg_0_xxyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 4);

                auto tg_0_xxzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 5);

                auto tg_0_xyyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 6);

                auto tg_0_xyyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 7);

                auto tg_0_xyzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 8);

                auto tg_0_xzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 9);

                auto tg_0_yyyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 10);

                auto tg_0_yyyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 11);

                auto tg_0_yyzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 12);

                auto tg_0_yzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 13);

                auto tg_0_zzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 14);

                auto tg_0_xxxx_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx);

                auto tg_0_xxxy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 1);

                auto tg_0_xxxz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 2);

                auto tg_0_xxyy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 3);

                auto tg_0_xxyz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 4);

                auto tg_0_xxzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 5);

                auto tg_0_xyyy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 6);

                auto tg_0_xyyz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 7);

                auto tg_0_xyzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 8);

                auto tg_0_xzzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 9);

                auto tg_0_yyyy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 10);

                auto tg_0_yyyz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 11);

                auto tg_0_yyzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 12);

                auto tg_0_yzzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 13);

                auto tg_0_zzzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 14);

                auto tg_0_xxx_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx);

                auto tg_0_xxy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 1);

                auto tg_0_xxz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 2);

                auto tg_0_xyy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 3);

                auto tg_0_xyz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 4);

                auto tg_0_xzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 5);

                auto tg_0_yyy_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 6);

                auto tg_0_yyz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 7);

                auto tg_0_yzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 8);

                auto tg_0_zzz_0 = primBuffer[pidx_g_0_3_m0].data(10 * idx + 9);

                auto tg_0_xxx_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx);

                auto tg_0_xxy_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 1);

                auto tg_0_xxz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 2);

                auto tg_0_xyy_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 3);

                auto tg_0_xyz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 4);

                auto tg_0_xzz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 5);

                auto tg_0_yyy_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 6);

                auto tg_0_yyz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 7);

                auto tg_0_yzz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 8);

                auto tg_0_zzz_1 = primBuffer[pidx_g_0_3_m1].data(10 * idx + 9);

                // set up pointers to integrals

                auto tg_0_xxxxx_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx);

                auto tg_0_xxxxy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 1);

                auto tg_0_xxxxz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 2);

                auto tg_0_xxxyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 3);

                auto tg_0_xxxyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 4);

                auto tg_0_xxxzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 5);

                auto tg_0_xxyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 6);

                auto tg_0_xxyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 7);

                auto tg_0_xxyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 8);

                auto tg_0_xxzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 9);

                auto tg_0_xyyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 10);

                auto tg_0_xyyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 11);

                auto tg_0_xyyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 12);

                auto tg_0_xyzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 13);

                auto tg_0_xzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 14);

                auto tg_0_yyyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 15);

                auto tg_0_yyyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 16);

                auto tg_0_yyyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 17);

                auto tg_0_yyzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 18);

                auto tg_0_yzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 19);

                auto tg_0_zzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 20);

                #pragma omp simd aligned(fn, fzb, qd_x, qd_y, qd_z, tg_0_xxx_0, tg_0_xxx_1, tg_0_xxxx_0, tg_0_xxxx_1, \
                                         tg_0_xxxxx_0, tg_0_xxxxy_0, tg_0_xxxxz_0, tg_0_xxxy_0, tg_0_xxxy_1, tg_0_xxxyy_0, \
                                         tg_0_xxxyz_0, tg_0_xxxz_0, tg_0_xxxz_1, tg_0_xxxzz_0, tg_0_xxy_0, tg_0_xxy_1, \
                                         tg_0_xxyy_0, tg_0_xxyy_1, tg_0_xxyyy_0, tg_0_xxyyz_0, tg_0_xxyz_0, tg_0_xxyz_1, \
                                         tg_0_xxyzz_0, tg_0_xxz_0, tg_0_xxz_1, tg_0_xxzz_0, tg_0_xxzz_1, tg_0_xxzzz_0, \
                                         tg_0_xyy_0, tg_0_xyy_1, tg_0_xyyy_0, tg_0_xyyy_1, tg_0_xyyyy_0, tg_0_xyyyz_0, \
                                         tg_0_xyyz_0, tg_0_xyyz_1, tg_0_xyyzz_0, tg_0_xyz_0, tg_0_xyz_1, tg_0_xyzz_0, \
                                         tg_0_xyzz_1, tg_0_xyzzz_0, tg_0_xzz_0, tg_0_xzz_1, tg_0_xzzz_0, tg_0_xzzz_1, \
                                         tg_0_xzzzz_0, tg_0_yyy_0, tg_0_yyy_1, tg_0_yyyy_0, tg_0_yyyy_1, tg_0_yyyyy_0, \
                                         tg_0_yyyyz_0, tg_0_yyyz_0, tg_0_yyyz_1, tg_0_yyyzz_0, tg_0_yyz_0, tg_0_yyz_1, \
                                         tg_0_yyzz_0, tg_0_yyzz_1, tg_0_yyzzz_0, tg_0_yzz_0, tg_0_yzz_1, tg_0_yzzz_0, \
                                         tg_0_yzzz_1, tg_0_yzzzz_0, tg_0_zzz_0, tg_0_zzz_1, tg_0_zzzz_0, tg_0_zzzz_1, \
                                         tg_0_zzzzz_0, wq_x, wq_y, wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fn = fn[j];

                    double fl1_fzb = fzb[j];
                    
                    double fr0 = qd_x[j]; double fr1 = wq_x[j];

                    tg_0_xxxxx_0[j] = fr0 * tg_0_xxxx_0[j] + fr1 * tg_0_xxxx_1[j] + 2.0 * fl1_fn * (tg_0_xxx_0[j] - fl1_fzb * tg_0_xxx_1[j]);

                    tg_0_xxxxy_0[j] = fr0 * tg_0_xxxy_0[j] + fr1 * tg_0_xxxy_1[j] + 1.5 * fl1_fn * (tg_0_xxy_0[j] - fl1_fzb * tg_0_xxy_1[j]);

                    tg_0_xxxxz_0[j] = fr0 * tg_0_xxxz_0[j] + fr1 * tg_0_xxxz_1[j] + 1.5 * fl1_fn * (tg_0_xxz_0[j] - fl1_fzb * tg_0_xxz_1[j]);

                    tg_0_xxxyy_0[j] = fr0 * tg_0_xxyy_0[j] + fr1 * tg_0_xxyy_1[j] + fl1_fn * (tg_0_xyy_0[j] - fl1_fzb * tg_0_xyy_1[j]);

                    tg_0_xxxyz_0[j] = fr0 * tg_0_xxyz_0[j] + fr1 * tg_0_xxyz_1[j] + fl1_fn * (tg_0_xyz_0[j] - fl1_fzb * tg_0_xyz_1[j]);

                    tg_0_xxxzz_0[j] = fr0 * tg_0_xxzz_0[j] + fr1 * tg_0_xxzz_1[j] + fl1_fn * (tg_0_xzz_0[j] - fl1_fzb * tg_0_xzz_1[j]);

                    tg_0_xxyyy_0[j] = fr0 * tg_0_xyyy_0[j] + fr1 * tg_0_xyyy_1[j] + 0.5 * fl1_fn * (tg_0_yyy_0[j] - fl1_fzb * tg_0_yyy_1[j]);

                    tg_0_xxyyz_0[j] = fr0 * tg_0_xyyz_0[j] + fr1 * tg_0_xyyz_1[j] + 0.5 * fl1_fn * (tg_0_yyz_0[j] - fl1_fzb * tg_0_yyz_1[j]);

                    tg_0_xxyzz_0[j] = fr0 * tg_0_xyzz_0[j] + fr1 * tg_0_xyzz_1[j] + 0.5 * fl1_fn * (tg_0_yzz_0[j] - fl1_fzb * tg_0_yzz_1[j]);

                    tg_0_xxzzz_0[j] = fr0 * tg_0_xzzz_0[j] + fr1 * tg_0_xzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzz_0[j] - fl1_fzb * tg_0_zzz_1[j]);

                    tg_0_xyyyy_0[j] = fr0 * tg_0_yyyy_0[j] + fr1 * tg_0_yyyy_1[j];

                    tg_0_xyyyz_0[j] = fr0 * tg_0_yyyz_0[j] + fr1 * tg_0_yyyz_1[j];

                    tg_0_xyyzz_0[j] = fr0 * tg_0_yyzz_0[j] + fr1 * tg_0_yyzz_1[j];

                    tg_0_xyzzz_0[j] = fr0 * tg_0_yzzz_0[j] + fr1 * tg_0_yzzz_1[j];

                    tg_0_xzzzz_0[j] = fr0 * tg_0_zzzz_0[j] + fr1 * tg_0_zzzz_1[j];
                    
                    fr0 = qd_y[j]; fr1 = wq_y[j];

                    tg_0_yyyyy_0[j] = fr0 * tg_0_yyyy_0[j] + fr1 * tg_0_yyyy_1[j] + 2.0 * fl1_fn * (tg_0_yyy_0[j] - fl1_fzb * tg_0_yyy_1[j]);

                    tg_0_yyyyz_0[j] = fr0 * tg_0_yyyz_0[j] + fr1 * tg_0_yyyz_1[j] + 1.5 * fl1_fn * (tg_0_yyz_0[j] - fl1_fzb * tg_0_yyz_1[j]);

                    tg_0_yyyzz_0[j] = fr0 * tg_0_yyzz_0[j] + fr1 * tg_0_yyzz_1[j] + fl1_fn * (tg_0_yzz_0[j] - fl1_fzb * tg_0_yzz_1[j]);

                    tg_0_yyzzz_0[j] = fr0 * tg_0_yzzz_0[j] + fr1 * tg_0_yzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzz_0[j] - fl1_fzb * tg_0_zzz_1[j]);

                    tg_0_yzzzz_0[j] = fr0 * tg_0_zzzz_0[j] + fr1 * tg_0_zzzz_1[j];

                    tg_0_zzzzz_0[j] = qd_z[j] * tg_0_zzzz_0[j] + wq_z[j] * tg_0_zzzz_1[j] + 2.0 * fl1_fn * (tg_0_zzz_0[j] - fl1_fzb * tg_0_zzz_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSHSS(      CMemBlock2D<double>* primBuffer,
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
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_5_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {5, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_5_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_4_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {4, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_4_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {4, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord + 1));

            auto pidx_g_3_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {3, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_3_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {3, -1, -1, -1}, {0, -1, -1, -1},
                                                                    1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_xxxx_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx);

                auto tg_xxxy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 1);

                auto tg_xxxz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 2);

                auto tg_xxyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 3);

                auto tg_xxyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 4);

                auto tg_xxzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 5);

                auto tg_xyyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 6);

                auto tg_xyyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 7);

                auto tg_xyzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 8);

                auto tg_xzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 9);

                auto tg_yyyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 10);

                auto tg_yyyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 11);

                auto tg_yyzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 12);

                auto tg_yzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 13);

                auto tg_zzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 14);

                auto tg_xxxx_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx);

                auto tg_xxxy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 1);

                auto tg_xxxz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 2);

                auto tg_xxyy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 3);

                auto tg_xxyz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 4);

                auto tg_xxzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 5);

                auto tg_xyyy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 6);

                auto tg_xyyz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 7);

                auto tg_xyzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 8);

                auto tg_xzzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 9);

                auto tg_yyyy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 10);

                auto tg_yyyz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 11);

                auto tg_yyzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 12);

                auto tg_yzzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 13);

                auto tg_zzzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 14);

                auto tg_xxx_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx);

                auto tg_xxy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 1);

                auto tg_xxz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 2);

                auto tg_xyy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 3);

                auto tg_xyz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 4);

                auto tg_xzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 5);

                auto tg_yyy_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 6);

                auto tg_yyz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 7);

                auto tg_yzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 8);

                auto tg_zzz_0_0 = primBuffer[pidx_g_3_0_m0].data(10 * idx + 9);

                auto tg_xxx_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx);

                auto tg_xxy_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 1);

                auto tg_xxz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 2);

                auto tg_xyy_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 3);

                auto tg_xyz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 4);

                auto tg_xzz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 5);

                auto tg_yyy_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 6);

                auto tg_yyz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 7);

                auto tg_yzz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 8);

                auto tg_zzz_0_1 = primBuffer[pidx_g_3_0_m1].data(10 * idx + 9);

                // set up pointers to integrals

                auto tg_xxxxx_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx);

                auto tg_xxxxy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 1);

                auto tg_xxxxz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 2);

                auto tg_xxxyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 3);

                auto tg_xxxyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 4);

                auto tg_xxxzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 5);

                auto tg_xxyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 6);

                auto tg_xxyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 7);

                auto tg_xxyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 8);

                auto tg_xxzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 9);

                auto tg_xyyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 10);

                auto tg_xyyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 11);

                auto tg_xyyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 12);

                auto tg_xyzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 13);

                auto tg_xzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 14);

                auto tg_yyyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 15);

                auto tg_yyyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 16);

                auto tg_yyyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 17);

                auto tg_yyzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 18);

                auto tg_yzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 19);

                auto tg_zzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 20);

                #pragma omp simd aligned(fza, tg_xxx_0_0, tg_xxx_0_1, tg_xxxx_0_0, tg_xxxx_0_1, tg_xxxxx_0_0, \
                                         tg_xxxxy_0_0, tg_xxxxz_0_0, tg_xxxy_0_0, tg_xxxy_0_1, tg_xxxyy_0_0, tg_xxxyz_0_0, \
                                         tg_xxxz_0_0, tg_xxxz_0_1, tg_xxxzz_0_0, tg_xxy_0_0, tg_xxy_0_1, tg_xxyy_0_0, \
                                         tg_xxyy_0_1, tg_xxyyy_0_0, tg_xxyyz_0_0, tg_xxyz_0_0, tg_xxyz_0_1, tg_xxyzz_0_0, \
                                         tg_xxz_0_0, tg_xxz_0_1, tg_xxzz_0_0, tg_xxzz_0_1, tg_xxzzz_0_0, tg_xyy_0_0, \
                                         tg_xyy_0_1, tg_xyyy_0_0, tg_xyyy_0_1, tg_xyyyy_0_0, tg_xyyyz_0_0, tg_xyyz_0_0, \
                                         tg_xyyz_0_1, tg_xyyzz_0_0, tg_xyz_0_0, tg_xyz_0_1, tg_xyzz_0_0, tg_xyzz_0_1, \
                                         tg_xyzzz_0_0, tg_xzz_0_0, tg_xzz_0_1, tg_xzzz_0_0, tg_xzzz_0_1, tg_xzzzz_0_0, \
                                         tg_yyy_0_0, tg_yyy_0_1, tg_yyyy_0_0, tg_yyyy_0_1, tg_yyyyy_0_0, tg_yyyyz_0_0, \
                                         tg_yyyz_0_0, tg_yyyz_0_1, tg_yyyzz_0_0, tg_yyz_0_0, tg_yyz_0_1, tg_yyzz_0_0, \
                                         tg_yyzz_0_1, tg_yyzzz_0_0, tg_yzz_0_0, tg_yzz_0_1, tg_yzzz_0_0, tg_yzzz_0_1, \
                                         tg_yzzzz_0_0, tg_zzz_0_0, tg_zzz_0_1, tg_zzzz_0_0, tg_zzzz_0_1, tg_zzzzz_0_0, wp_x, \
                                         wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fza = fza[j];
                    
                    double fr1 = wp_x[j];

                    tg_xxxxx_0_0[j] = pb_x * tg_xxxx_0_0[j] + fr1 * tg_xxxx_0_1[j] + 2.0 * fl1_fx * (tg_xxx_0_0[j] -  fl1_fza * tg_xxx_0_1[j]);

                    tg_xxxxy_0_0[j] = pb_x * tg_xxxy_0_0[j] + fr1 * tg_xxxy_0_1[j] + 1.5 * fl1_fx * (tg_xxy_0_0[j] -  fl1_fza * tg_xxy_0_1[j]);

                    tg_xxxxz_0_0[j] = pb_x * tg_xxxz_0_0[j] + fr1 * tg_xxxz_0_1[j] + 1.5 * fl1_fx * (tg_xxz_0_0[j] -  fl1_fza * tg_xxz_0_1[j]);

                    tg_xxxyy_0_0[j] = pb_x * tg_xxyy_0_0[j] + fr1 * tg_xxyy_0_1[j] + fl1_fx * (tg_xyy_0_0[j] -  fl1_fza * tg_xyy_0_1[j]);

                    tg_xxxyz_0_0[j] = pb_x * tg_xxyz_0_0[j] + fr1 * tg_xxyz_0_1[j] + fl1_fx * (tg_xyz_0_0[j] -  fl1_fza * tg_xyz_0_1[j]);

                    tg_xxxzz_0_0[j] = pb_x * tg_xxzz_0_0[j] + fr1 * tg_xxzz_0_1[j] + fl1_fx * (tg_xzz_0_0[j] -  fl1_fza * tg_xzz_0_1[j]);

                    tg_xxyyy_0_0[j] = pb_x * tg_xyyy_0_0[j] + fr1 * tg_xyyy_0_1[j] + 0.5 * fl1_fx * (tg_yyy_0_0[j] -  fl1_fza * tg_yyy_0_1[j]);

                    tg_xxyyz_0_0[j] = pb_x * tg_xyyz_0_0[j] + fr1 * tg_xyyz_0_1[j] + 0.5 * fl1_fx * (tg_yyz_0_0[j] -  fl1_fza * tg_yyz_0_1[j]);

                    tg_xxyzz_0_0[j] = pb_x * tg_xyzz_0_0[j] + fr1 * tg_xyzz_0_1[j] + 0.5 * fl1_fx * (tg_yzz_0_0[j] -  fl1_fza * tg_yzz_0_1[j]);

                    tg_xxzzz_0_0[j] = pb_x * tg_xzzz_0_0[j] + fr1 * tg_xzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzz_0_0[j] -  fl1_fza * tg_zzz_0_1[j]);

                    tg_xyyyy_0_0[j] = pb_x * tg_yyyy_0_0[j] + fr1 * tg_yyyy_0_1[j];

                    tg_xyyyz_0_0[j] = pb_x * tg_yyyz_0_0[j] + fr1 * tg_yyyz_0_1[j];

                    tg_xyyzz_0_0[j] = pb_x * tg_yyzz_0_0[j] + fr1 * tg_yyzz_0_1[j];

                    tg_xyzzz_0_0[j] = pb_x * tg_yzzz_0_0[j] + fr1 * tg_yzzz_0_1[j];

                    tg_xzzzz_0_0[j] = pb_x * tg_zzzz_0_0[j] + fr1 * tg_zzzz_0_1[j];
                    
                    fr1 = wp_y[j];

                    tg_yyyyy_0_0[j] = pb_y * tg_yyyy_0_0[j] + fr1 * tg_yyyy_0_1[j] + 2.0 * fl1_fx * (tg_yyy_0_0[j] -  fl1_fza * tg_yyy_0_1[j]);

                    tg_yyyyz_0_0[j] = pb_y * tg_yyyz_0_0[j] + fr1 * tg_yyyz_0_1[j] + 1.5 * fl1_fx * (tg_yyz_0_0[j] -  fl1_fza * tg_yyz_0_1[j]);

                    tg_yyyzz_0_0[j] = pb_y * tg_yyzz_0_0[j] + fr1 * tg_yyzz_0_1[j] + fl1_fx * (tg_yzz_0_0[j] -  fl1_fza * tg_yzz_0_1[j]);

                    tg_yyzzz_0_0[j] = pb_y * tg_yzzz_0_0[j] + fr1 * tg_yzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzz_0_0[j] -  fl1_fza * tg_zzz_0_1[j]);

                    tg_yzzzz_0_0[j] = pb_y * tg_zzzz_0_0[j] + fr1 * tg_zzzz_0_1[j];

                    tg_zzzzz_0_0[j] = pb_z * tg_zzzz_0_0[j] + wp_z[j] * tg_zzzz_0_1[j] + 2.0 * fl1_fx * (tg_zzz_0_0[j] -  fl1_fza * tg_zzz_0_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSSSI(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto fn = ketGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {6, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {6, -1, -1, -1},
                                                                    1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_6_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {5, -1, -1, -1},
                                                                    1, 1, iord));

            auto pidx_g_0_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {5, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_4_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {4, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_4_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {4, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fzb = osFactors.data(4 * idx + 3);

                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xxxxx_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx);

                auto tg_0_xxxxy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 1);

                auto tg_0_xxxxz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 2);

                auto tg_0_xxxyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 3);

                auto tg_0_xxxyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 4);

                auto tg_0_xxxzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 5);

                auto tg_0_xxyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 6);

                auto tg_0_xxyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 7);

                auto tg_0_xxyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 8);

                auto tg_0_xxzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 9);

                auto tg_0_xyyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 10);

                auto tg_0_xyyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 11);

                auto tg_0_xyyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 12);

                auto tg_0_xyzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 13);

                auto tg_0_xzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 14);

                auto tg_0_yyyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 15);

                auto tg_0_yyyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 16);

                auto tg_0_yyyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 17);

                auto tg_0_yyzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 18);

                auto tg_0_yzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 19);

                auto tg_0_zzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 20);

                auto tg_0_xxxxx_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx);

                auto tg_0_xxxxy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 1);

                auto tg_0_xxxxz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 2);

                auto tg_0_xxxyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 3);

                auto tg_0_xxxyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 4);

                auto tg_0_xxxzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 5);

                auto tg_0_xxyyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 6);

                auto tg_0_xxyyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 7);

                auto tg_0_xxyzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 8);

                auto tg_0_xxzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 9);

                auto tg_0_xyyyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 10);

                auto tg_0_xyyyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 11);

                auto tg_0_xyyzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 12);

                auto tg_0_xyzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 13);

                auto tg_0_xzzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 14);

                auto tg_0_yyyyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 15);

                auto tg_0_yyyyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 16);

                auto tg_0_yyyzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 17);

                auto tg_0_yyzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 18);

                auto tg_0_yzzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 19);

                auto tg_0_zzzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 20);

                auto tg_0_xxxx_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx);

                auto tg_0_xxxy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 1);

                auto tg_0_xxxz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 2);

                auto tg_0_xxyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 3);

                auto tg_0_xxyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 4);

                auto tg_0_xxzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 5);

                auto tg_0_xyyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 6);

                auto tg_0_xyyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 7);

                auto tg_0_xyzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 8);

                auto tg_0_xzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 9);

                auto tg_0_yyyy_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 10);

                auto tg_0_yyyz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 11);

                auto tg_0_yyzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 12);

                auto tg_0_yzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 13);

                auto tg_0_zzzz_0 = primBuffer[pidx_g_0_4_m0].data(15 * idx + 14);

                auto tg_0_xxxx_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx);

                auto tg_0_xxxy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 1);

                auto tg_0_xxxz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 2);

                auto tg_0_xxyy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 3);

                auto tg_0_xxyz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 4);

                auto tg_0_xxzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 5);

                auto tg_0_xyyy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 6);

                auto tg_0_xyyz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 7);

                auto tg_0_xyzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 8);

                auto tg_0_xzzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 9);

                auto tg_0_yyyy_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 10);

                auto tg_0_yyyz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 11);

                auto tg_0_yyzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 12);

                auto tg_0_yzzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 13);

                auto tg_0_zzzz_1 = primBuffer[pidx_g_0_4_m1].data(15 * idx + 14);

                // set up pointers to integrals

                auto tg_0_xxxxxx_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx);

                auto tg_0_xxxxxy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 1);

                auto tg_0_xxxxxz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 2);

                auto tg_0_xxxxyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 3);

                auto tg_0_xxxxyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 4);

                auto tg_0_xxxxzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 5);

                auto tg_0_xxxyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 6);

                auto tg_0_xxxyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 7);

                auto tg_0_xxxyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 8);

                auto tg_0_xxxzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 9);

                auto tg_0_xxyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 10);

                auto tg_0_xxyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 11);

                auto tg_0_xxyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 12);

                auto tg_0_xxyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 13);

                auto tg_0_xxzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 14);

                auto tg_0_xyyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 15);

                auto tg_0_xyyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 16);

                auto tg_0_xyyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 17);

                auto tg_0_xyyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 18);

                auto tg_0_xyzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 19);

                auto tg_0_xzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 20);

                auto tg_0_yyyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 21);

                auto tg_0_yyyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 22);

                auto tg_0_yyyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 23);

                auto tg_0_yyyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 24);

                auto tg_0_yyzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 25);

                auto tg_0_yzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 26);

                auto tg_0_zzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 27);

                #pragma omp simd aligned(fn, fzb, qd_x, qd_y, qd_z, tg_0_xxxx_0, tg_0_xxxx_1, tg_0_xxxxx_0, \
                                         tg_0_xxxxx_1, tg_0_xxxxxx_0, tg_0_xxxxxy_0, tg_0_xxxxxz_0, tg_0_xxxxy_0, \
                                         tg_0_xxxxy_1, tg_0_xxxxyy_0, tg_0_xxxxyz_0, tg_0_xxxxz_0, tg_0_xxxxz_1, \
                                         tg_0_xxxxzz_0, tg_0_xxxy_0, tg_0_xxxy_1, tg_0_xxxyy_0, tg_0_xxxyy_1, tg_0_xxxyyy_0, \
                                         tg_0_xxxyyz_0, tg_0_xxxyz_0, tg_0_xxxyz_1, tg_0_xxxyzz_0, tg_0_xxxz_0, tg_0_xxxz_1, \
                                         tg_0_xxxzz_0, tg_0_xxxzz_1, tg_0_xxxzzz_0, tg_0_xxyy_0, tg_0_xxyy_1, tg_0_xxyyy_0, \
                                         tg_0_xxyyy_1, tg_0_xxyyyy_0, tg_0_xxyyyz_0, tg_0_xxyyz_0, tg_0_xxyyz_1, \
                                         tg_0_xxyyzz_0, tg_0_xxyz_0, tg_0_xxyz_1, tg_0_xxyzz_0, tg_0_xxyzz_1, tg_0_xxyzzz_0, \
                                         tg_0_xxzz_0, tg_0_xxzz_1, tg_0_xxzzz_0, tg_0_xxzzz_1, tg_0_xxzzzz_0, tg_0_xyyy_0, \
                                         tg_0_xyyy_1, tg_0_xyyyy_0, tg_0_xyyyy_1, tg_0_xyyyyy_0, tg_0_xyyyyz_0, \
                                         tg_0_xyyyz_0, tg_0_xyyyz_1, tg_0_xyyyzz_0, tg_0_xyyz_0, tg_0_xyyz_1, tg_0_xyyzz_0, \
                                         tg_0_xyyzz_1, tg_0_xyyzzz_0, tg_0_xyzz_0, tg_0_xyzz_1, tg_0_xyzzz_0, tg_0_xyzzz_1, \
                                         tg_0_xyzzzz_0, tg_0_xzzz_0, tg_0_xzzz_1, tg_0_xzzzz_0, tg_0_xzzzz_1, tg_0_xzzzzz_0, \
                                         tg_0_yyyy_0, tg_0_yyyy_1, tg_0_yyyyy_0, tg_0_yyyyy_1, tg_0_yyyyyy_0, \
                                         tg_0_yyyyyz_0, tg_0_yyyyz_0, tg_0_yyyyz_1, tg_0_yyyyzz_0, tg_0_yyyz_0, tg_0_yyyz_1, \
                                         tg_0_yyyzz_0, tg_0_yyyzz_1, tg_0_yyyzzz_0, tg_0_yyzz_0, tg_0_yyzz_1, tg_0_yyzzz_0, \
                                         tg_0_yyzzz_1, tg_0_yyzzzz_0, tg_0_yzzz_0, tg_0_yzzz_1, tg_0_yzzzz_0, tg_0_yzzzz_1, \
                                         tg_0_yzzzzz_0, tg_0_zzzz_0, tg_0_zzzz_1, tg_0_zzzzz_0, tg_0_zzzzz_1, tg_0_zzzzzz_0, \
                                         wq_x, wq_y, wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fn = fn[j];

                    double fl1_fzb = fzb[j];
                    
                    double fr0 = qd_x[j]; double fr1 = wq_x[j];

                    tg_0_xxxxxx_0[j] = fr0 * tg_0_xxxxx_0[j] + fr1 * tg_0_xxxxx_1[j] + 2.5 * fl1_fn * (tg_0_xxxx_0[j] - fl1_fzb * tg_0_xxxx_1[j]);

                    tg_0_xxxxxy_0[j] = fr0 * tg_0_xxxxy_0[j] + fr1 * tg_0_xxxxy_1[j] + 2.0 * fl1_fn * (tg_0_xxxy_0[j] - fl1_fzb * tg_0_xxxy_1[j]);

                    tg_0_xxxxxz_0[j] = fr0 * tg_0_xxxxz_0[j] + fr1 * tg_0_xxxxz_1[j] + 2.0 * fl1_fn * (tg_0_xxxz_0[j] - fl1_fzb * tg_0_xxxz_1[j]);

                    tg_0_xxxxyy_0[j] = fr0 * tg_0_xxxyy_0[j] + fr1 * tg_0_xxxyy_1[j] + 1.5 * fl1_fn * (tg_0_xxyy_0[j] - fl1_fzb * tg_0_xxyy_1[j]);

                    tg_0_xxxxyz_0[j] = fr0 * tg_0_xxxyz_0[j] + fr1 * tg_0_xxxyz_1[j] + 1.5 * fl1_fn * (tg_0_xxyz_0[j] - fl1_fzb * tg_0_xxyz_1[j]);

                    tg_0_xxxxzz_0[j] = fr0 * tg_0_xxxzz_0[j] + fr1 * tg_0_xxxzz_1[j] + 1.5 * fl1_fn * (tg_0_xxzz_0[j] - fl1_fzb * tg_0_xxzz_1[j]);

                    tg_0_xxxyyy_0[j] = fr0 * tg_0_xxyyy_0[j] + fr1 * tg_0_xxyyy_1[j] + fl1_fn * (tg_0_xyyy_0[j] - fl1_fzb * tg_0_xyyy_1[j]);

                    tg_0_xxxyyz_0[j] = fr0 * tg_0_xxyyz_0[j] + fr1 * tg_0_xxyyz_1[j] + fl1_fn * (tg_0_xyyz_0[j] - fl1_fzb * tg_0_xyyz_1[j]);

                    tg_0_xxxyzz_0[j] = fr0 * tg_0_xxyzz_0[j] + fr1 * tg_0_xxyzz_1[j] + fl1_fn * (tg_0_xyzz_0[j] - fl1_fzb * tg_0_xyzz_1[j]);

                    tg_0_xxxzzz_0[j] = fr0 * tg_0_xxzzz_0[j] + fr1 * tg_0_xxzzz_1[j] + fl1_fn * (tg_0_xzzz_0[j] - fl1_fzb * tg_0_xzzz_1[j]);

                    tg_0_xxyyyy_0[j] = fr0 * tg_0_xyyyy_0[j] + fr1 * tg_0_xyyyy_1[j] + 0.5 * fl1_fn * (tg_0_yyyy_0[j] - fl1_fzb * tg_0_yyyy_1[j]);

                    tg_0_xxyyyz_0[j] = fr0 * tg_0_xyyyz_0[j] + fr1 * tg_0_xyyyz_1[j] + 0.5 * fl1_fn * (tg_0_yyyz_0[j] - fl1_fzb * tg_0_yyyz_1[j]);

                    tg_0_xxyyzz_0[j] = fr0 * tg_0_xyyzz_0[j] + fr1 * tg_0_xyyzz_1[j] + 0.5 * fl1_fn * (tg_0_yyzz_0[j] - fl1_fzb * tg_0_yyzz_1[j]);

                    tg_0_xxyzzz_0[j] = fr0 * tg_0_xyzzz_0[j] + fr1 * tg_0_xyzzz_1[j] + 0.5 * fl1_fn * (tg_0_yzzz_0[j] - fl1_fzb * tg_0_yzzz_1[j]);

                    tg_0_xxzzzz_0[j] = fr0 * tg_0_xzzzz_0[j] + fr1 * tg_0_xzzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzzz_0[j] -  fl1_fzb * tg_0_zzzz_1[j]);

                    tg_0_xyyyyy_0[j] = fr0 * tg_0_yyyyy_0[j] + fr1 * tg_0_yyyyy_1[j];

                    tg_0_xyyyyz_0[j] = fr0 * tg_0_yyyyz_0[j] + fr1 * tg_0_yyyyz_1[j];

                    tg_0_xyyyzz_0[j] = fr0 * tg_0_yyyzz_0[j] + fr1 * tg_0_yyyzz_1[j];

                    tg_0_xyyzzz_0[j] = fr0 * tg_0_yyzzz_0[j] + fr1 * tg_0_yyzzz_1[j];

                    tg_0_xyzzzz_0[j] = fr0 * tg_0_yzzzz_0[j] + fr1 * tg_0_yzzzz_1[j];

                    tg_0_xzzzzz_0[j] = fr0 * tg_0_zzzzz_0[j] + fr1 * tg_0_zzzzz_1[j];
                    
                    fr0 = qd_y[j]; fr1 = wq_y[j];

                    tg_0_yyyyyy_0[j] = fr0 * tg_0_yyyyy_0[j] + fr1 * tg_0_yyyyy_1[j] + 2.5 * fl1_fn * (tg_0_yyyy_0[j] - fl1_fzb * tg_0_yyyy_1[j]);

                    tg_0_yyyyyz_0[j] = fr0 * tg_0_yyyyz_0[j] + fr1 * tg_0_yyyyz_1[j] + 2.0 * fl1_fn * (tg_0_yyyz_0[j] - fl1_fzb * tg_0_yyyz_1[j]);

                    tg_0_yyyyzz_0[j] = fr0 * tg_0_yyyzz_0[j] + fr1 * tg_0_yyyzz_1[j] + 1.5 * fl1_fn * (tg_0_yyzz_0[j] - fl1_fzb * tg_0_yyzz_1[j]);

                    tg_0_yyyzzz_0[j] = fr0 * tg_0_yyzzz_0[j] + fr1 * tg_0_yyzzz_1[j] + fl1_fn * (tg_0_yzzz_0[j] - fl1_fzb * tg_0_yzzz_1[j]);

                    tg_0_yyzzzz_0[j] = fr0 * tg_0_yzzzz_0[j] + fr1 * tg_0_yzzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzzz_0[j] - fl1_fzb * tg_0_zzzz_1[j]);

                    tg_0_yzzzzz_0[j] = fr0 * tg_0_zzzzz_0[j] + fr1 * tg_0_zzzzz_1[j];

                    tg_0_zzzzzz_0[j] = qd_z[j] * tg_0_zzzzz_0[j] + wq_z[j] * tg_0_zzzzz_1[j] + 2.5 * fl1_fn * (tg_0_zzzz_0[j] - fl1_fzb * tg_0_zzzz_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSISS(      CMemBlock2D<double>* primBuffer,
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
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_6_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {6, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_6_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_5_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {5, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_5_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {5, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_4_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {4, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_4_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {4, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_xxxxx_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx);

                auto tg_xxxxy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 1);

                auto tg_xxxxz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 2);

                auto tg_xxxyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 3);

                auto tg_xxxyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 4);

                auto tg_xxxzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 5);

                auto tg_xxyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 6);

                auto tg_xxyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 7);

                auto tg_xxyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 8);

                auto tg_xxzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 9);

                auto tg_xyyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 10);

                auto tg_xyyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 11);

                auto tg_xyyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 12);

                auto tg_xyzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 13);

                auto tg_xzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 14);

                auto tg_yyyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 15);

                auto tg_yyyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 16);

                auto tg_yyyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 17);

                auto tg_yyzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 18);

                auto tg_yzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 19);

                auto tg_zzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 20);

                auto tg_xxxxx_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx);

                auto tg_xxxxy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 1);

                auto tg_xxxxz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 2);

                auto tg_xxxyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 3);

                auto tg_xxxyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 4);

                auto tg_xxxzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 5);

                auto tg_xxyyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 6);

                auto tg_xxyyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 7);

                auto tg_xxyzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 8);

                auto tg_xxzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 9);

                auto tg_xyyyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 10);

                auto tg_xyyyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 11);

                auto tg_xyyzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 12);

                auto tg_xyzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 13);

                auto tg_xzzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 14);

                auto tg_yyyyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 15);

                auto tg_yyyyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 16);

                auto tg_yyyzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 17);

                auto tg_yyzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 18);

                auto tg_yzzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 19);

                auto tg_zzzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 20);

                auto tg_xxxx_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx);

                auto tg_xxxy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 1);

                auto tg_xxxz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 2);

                auto tg_xxyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 3);

                auto tg_xxyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 4);

                auto tg_xxzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 5);

                auto tg_xyyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 6);

                auto tg_xyyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 7);

                auto tg_xyzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 8);

                auto tg_xzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 9);

                auto tg_yyyy_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 10);

                auto tg_yyyz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 11);

                auto tg_yyzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 12);

                auto tg_yzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 13);

                auto tg_zzzz_0_0 = primBuffer[pidx_g_4_0_m0].data(15 * idx + 14);

                auto tg_xxxx_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx);

                auto tg_xxxy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 1);

                auto tg_xxxz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 2);

                auto tg_xxyy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 3);

                auto tg_xxyz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 4);

                auto tg_xxzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 5);

                auto tg_xyyy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 6);

                auto tg_xyyz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 7);

                auto tg_xyzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 8);

                auto tg_xzzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 9);

                auto tg_yyyy_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 10);

                auto tg_yyyz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 11);

                auto tg_yyzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 12);

                auto tg_yzzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 13);

                auto tg_zzzz_0_1 = primBuffer[pidx_g_4_0_m1].data(15 * idx + 14);

                // set up pointers to integrals

                auto tg_xxxxxx_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx);

                auto tg_xxxxxy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 1);

                auto tg_xxxxxz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 2);

                auto tg_xxxxyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 3);

                auto tg_xxxxyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 4);

                auto tg_xxxxzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 5);

                auto tg_xxxyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 6);

                auto tg_xxxyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 7);

                auto tg_xxxyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 8);

                auto tg_xxxzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 9);

                auto tg_xxyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 10);

                auto tg_xxyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 11);

                auto tg_xxyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 12);

                auto tg_xxyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 13);

                auto tg_xxzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 14);

                auto tg_xyyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 15);

                auto tg_xyyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 16);

                auto tg_xyyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 17);

                auto tg_xyyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 18);

                auto tg_xyzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 19);

                auto tg_xzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 20);

                auto tg_yyyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 21);

                auto tg_yyyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 22);

                auto tg_yyyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 23);

                auto tg_yyyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 24);

                auto tg_yyzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 25);

                auto tg_yzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 26);

                auto tg_zzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 27);

                #pragma omp simd aligned(fza, tg_xxxx_0_0, tg_xxxx_0_1, tg_xxxxx_0_0, tg_xxxxx_0_1, \
                                         tg_xxxxxx_0_0, tg_xxxxxy_0_0, tg_xxxxxz_0_0, tg_xxxxy_0_0, tg_xxxxy_0_1, \
                                         tg_xxxxyy_0_0, tg_xxxxyz_0_0, tg_xxxxz_0_0, tg_xxxxz_0_1, tg_xxxxzz_0_0, \
                                         tg_xxxy_0_0, tg_xxxy_0_1, tg_xxxyy_0_0, tg_xxxyy_0_1, tg_xxxyyy_0_0, \
                                         tg_xxxyyz_0_0, tg_xxxyz_0_0, tg_xxxyz_0_1, tg_xxxyzz_0_0, tg_xxxz_0_0, tg_xxxz_0_1, \
                                         tg_xxxzz_0_0, tg_xxxzz_0_1, tg_xxxzzz_0_0, tg_xxyy_0_0, tg_xxyy_0_1, tg_xxyyy_0_0, \
                                         tg_xxyyy_0_1, tg_xxyyyy_0_0, tg_xxyyyz_0_0, tg_xxyyz_0_0, tg_xxyyz_0_1, \
                                         tg_xxyyzz_0_0, tg_xxyz_0_0, tg_xxyz_0_1, tg_xxyzz_0_0, tg_xxyzz_0_1, tg_xxyzzz_0_0, \
                                         tg_xxzz_0_0, tg_xxzz_0_1, tg_xxzzz_0_0, tg_xxzzz_0_1, tg_xxzzzz_0_0, tg_xyyy_0_0, \
                                         tg_xyyy_0_1, tg_xyyyy_0_0, tg_xyyyy_0_1, tg_xyyyyy_0_0, tg_xyyyyz_0_0, \
                                         tg_xyyyz_0_0, tg_xyyyz_0_1, tg_xyyyzz_0_0, tg_xyyz_0_0, tg_xyyz_0_1, tg_xyyzz_0_0, \
                                         tg_xyyzz_0_1, tg_xyyzzz_0_0, tg_xyzz_0_0, tg_xyzz_0_1, tg_xyzzz_0_0, tg_xyzzz_0_1, \
                                         tg_xyzzzz_0_0, tg_xzzz_0_0, tg_xzzz_0_1, tg_xzzzz_0_0, tg_xzzzz_0_1, tg_xzzzzz_0_0, \
                                         tg_yyyy_0_0, tg_yyyy_0_1, tg_yyyyy_0_0, tg_yyyyy_0_1, tg_yyyyyy_0_0, \
                                         tg_yyyyyz_0_0, tg_yyyyz_0_0, tg_yyyyz_0_1, tg_yyyyzz_0_0, tg_yyyz_0_0, tg_yyyz_0_1, \
                                         tg_yyyzz_0_0, tg_yyyzz_0_1, tg_yyyzzz_0_0, tg_yyzz_0_0, tg_yyzz_0_1, tg_yyzzz_0_0, \
                                         tg_yyzzz_0_1, tg_yyzzzz_0_0, tg_yzzz_0_0, tg_yzzz_0_1, tg_yzzzz_0_0, tg_yzzzz_0_1, \
                                         tg_yzzzzz_0_0, tg_zzzz_0_0, tg_zzzz_0_1, tg_zzzzz_0_0, tg_zzzzz_0_1, tg_zzzzzz_0_0, \
                                         wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fza = fza[j];
                    
                    double fr1 = wp_x[j];

                    tg_xxxxxx_0_0[j] = pb_x * tg_xxxxx_0_0[j] + fr1 * tg_xxxxx_0_1[j] + 2.5 * fl1_fx * (tg_xxxx_0_0[j] -  fl1_fza * tg_xxxx_0_1[j]);

                    tg_xxxxxy_0_0[j] = pb_x * tg_xxxxy_0_0[j] + fr1 * tg_xxxxy_0_1[j] + 2.0 * fl1_fx * (tg_xxxy_0_0[j] -  fl1_fza * tg_xxxy_0_1[j]);

                    tg_xxxxxz_0_0[j] = pb_x * tg_xxxxz_0_0[j] + fr1 * tg_xxxxz_0_1[j] + 2.0 * fl1_fx * (tg_xxxz_0_0[j] -  fl1_fza * tg_xxxz_0_1[j]);

                    tg_xxxxyy_0_0[j] = pb_x * tg_xxxyy_0_0[j] + fr1 * tg_xxxyy_0_1[j] + 1.5 * fl1_fx * (tg_xxyy_0_0[j] -  fl1_fza * tg_xxyy_0_1[j]);

                    tg_xxxxyz_0_0[j] = pb_x * tg_xxxyz_0_0[j] + fr1 * tg_xxxyz_0_1[j] + 1.5 * fl1_fx * (tg_xxyz_0_0[j] -  fl1_fza * tg_xxyz_0_1[j]);

                    tg_xxxxzz_0_0[j] = pb_x * tg_xxxzz_0_0[j] + fr1 * tg_xxxzz_0_1[j] + 1.5 * fl1_fx * (tg_xxzz_0_0[j] -  fl1_fza * tg_xxzz_0_1[j]);

                    tg_xxxyyy_0_0[j] = pb_x * tg_xxyyy_0_0[j] + fr1 * tg_xxyyy_0_1[j] + fl1_fx * (tg_xyyy_0_0[j] -  fl1_fza * tg_xyyy_0_1[j]);

                    tg_xxxyyz_0_0[j] = pb_x * tg_xxyyz_0_0[j] + fr1 * tg_xxyyz_0_1[j] + fl1_fx * (tg_xyyz_0_0[j] -  fl1_fza * tg_xyyz_0_1[j]);

                    tg_xxxyzz_0_0[j] = pb_x * tg_xxyzz_0_0[j] + fr1 * tg_xxyzz_0_1[j] + fl1_fx * (tg_xyzz_0_0[j] -  fl1_fza * tg_xyzz_0_1[j]);

                    tg_xxxzzz_0_0[j] = pb_x * tg_xxzzz_0_0[j] + fr1 * tg_xxzzz_0_1[j] + fl1_fx * (tg_xzzz_0_0[j] -  fl1_fza * tg_xzzz_0_1[j]);

                    tg_xxyyyy_0_0[j] = pb_x * tg_xyyyy_0_0[j] + fr1 * tg_xyyyy_0_1[j] + 0.5 * fl1_fx * (tg_yyyy_0_0[j] -  fl1_fza * tg_yyyy_0_1[j]);

                    tg_xxyyyz_0_0[j] = pb_x * tg_xyyyz_0_0[j] + fr1 * tg_xyyyz_0_1[j] + 0.5 * fl1_fx * (tg_yyyz_0_0[j] -  fl1_fza * tg_yyyz_0_1[j]);

                    tg_xxyyzz_0_0[j] = pb_x * tg_xyyzz_0_0[j] + fr1 * tg_xyyzz_0_1[j] + 0.5 * fl1_fx * (tg_yyzz_0_0[j] -  fl1_fza * tg_yyzz_0_1[j]);

                    tg_xxyzzz_0_0[j] = pb_x * tg_xyzzz_0_0[j] + fr1 * tg_xyzzz_0_1[j] + 0.5 * fl1_fx * (tg_yzzz_0_0[j] -  fl1_fza * tg_yzzz_0_1[j]);

                    tg_xxzzzz_0_0[j] = pb_x * tg_xzzzz_0_0[j] + fr1 * tg_xzzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzzz_0_0[j] -  fl1_fza * tg_zzzz_0_1[j]);

                    tg_xyyyyy_0_0[j] = pb_x * tg_yyyyy_0_0[j] + fr1 * tg_yyyyy_0_1[j];

                    tg_xyyyyz_0_0[j] = pb_x * tg_yyyyz_0_0[j] + fr1 * tg_yyyyz_0_1[j];

                    tg_xyyyzz_0_0[j] = pb_x * tg_yyyzz_0_0[j] + fr1 * tg_yyyzz_0_1[j];

                    tg_xyyzzz_0_0[j] = pb_x * tg_yyzzz_0_0[j] + fr1 * tg_yyzzz_0_1[j];

                    tg_xyzzzz_0_0[j] = pb_x * tg_yzzzz_0_0[j] + fr1 * tg_yzzzz_0_1[j];

                    tg_xzzzzz_0_0[j] = pb_x * tg_zzzzz_0_0[j] + fr1 * tg_zzzzz_0_1[j];
                    
                    fr1 = wp_y[j];

                    tg_yyyyyy_0_0[j] = pb_y * tg_yyyyy_0_0[j] + fr1 * tg_yyyyy_0_1[j] + 2.5 * fl1_fx * (tg_yyyy_0_0[j] -  fl1_fza * tg_yyyy_0_1[j]);

                    tg_yyyyyz_0_0[j] = pb_y * tg_yyyyz_0_0[j] + fr1 * tg_yyyyz_0_1[j] + 2.0 * fl1_fx * (tg_yyyz_0_0[j] -  fl1_fza * tg_yyyz_0_1[j]);

                    tg_yyyyzz_0_0[j] = pb_y * tg_yyyzz_0_0[j] + fr1 * tg_yyyzz_0_1[j] + 1.5 * fl1_fx * (tg_yyzz_0_0[j] -  fl1_fza * tg_yyzz_0_1[j]);

                    tg_yyyzzz_0_0[j] = pb_y * tg_yyzzz_0_0[j] + fr1 * tg_yyzzz_0_1[j] + fl1_fx * (tg_yzzz_0_0[j] -  fl1_fza * tg_yzzz_0_1[j]);

                    tg_yyzzzz_0_0[j] = pb_y * tg_yzzzz_0_0[j] + fr1 * tg_yzzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzzz_0_0[j] -  fl1_fza * tg_zzzz_0_1[j]);

                    tg_yzzzzz_0_0[j] = pb_y * tg_zzzzz_0_0[j] + fr1 * tg_zzzzz_0_1[j];

                    tg_zzzzzz_0_0[j] = pb_z * tg_zzzzz_0_0[j] + wp_z[j] * tg_zzzzz_0_1[j] + 2.5 * fl1_fx * (tg_zzzz_0_0[j] -  fl1_fza * tg_zzzz_0_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSSSK(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto fn = ketGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {7, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {7, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_7_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {6, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {6, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_5_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {5, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_5_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {5, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fzb = osFactors.data(4 * idx + 3);

                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xxxxxx_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx);

                auto tg_0_xxxxxy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 1);

                auto tg_0_xxxxxz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 2);

                auto tg_0_xxxxyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 3);

                auto tg_0_xxxxyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 4);

                auto tg_0_xxxxzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 5);

                auto tg_0_xxxyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 6);

                auto tg_0_xxxyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 7);

                auto tg_0_xxxyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 8);

                auto tg_0_xxxzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 9);

                auto tg_0_xxyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 10);

                auto tg_0_xxyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 11);

                auto tg_0_xxyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 12);

                auto tg_0_xxyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 13);

                auto tg_0_xxzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 14);

                auto tg_0_xyyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 15);

                auto tg_0_xyyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 16);

                auto tg_0_xyyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 17);

                auto tg_0_xyyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 18);

                auto tg_0_xyzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 19);

                auto tg_0_xzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 20);

                auto tg_0_yyyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 21);

                auto tg_0_yyyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 22);

                auto tg_0_yyyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 23);

                auto tg_0_yyyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 24);

                auto tg_0_yyzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 25);

                auto tg_0_yzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 26);

                auto tg_0_zzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 27);

                auto tg_0_xxxxxx_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx);

                auto tg_0_xxxxxy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 1);

                auto tg_0_xxxxxz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 2);

                auto tg_0_xxxxyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 3);

                auto tg_0_xxxxyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 4);

                auto tg_0_xxxxzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 5);

                auto tg_0_xxxyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 6);

                auto tg_0_xxxyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 7);

                auto tg_0_xxxyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 8);

                auto tg_0_xxxzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 9);

                auto tg_0_xxyyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 10);

                auto tg_0_xxyyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 11);

                auto tg_0_xxyyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 12);

                auto tg_0_xxyzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 13);

                auto tg_0_xxzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 14);

                auto tg_0_xyyyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 15);

                auto tg_0_xyyyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 16);

                auto tg_0_xyyyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 17);

                auto tg_0_xyyzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 18);

                auto tg_0_xyzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 19);

                auto tg_0_xzzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 20);

                auto tg_0_yyyyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 21);

                auto tg_0_yyyyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 22);

                auto tg_0_yyyyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 23);

                auto tg_0_yyyzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 24);

                auto tg_0_yyzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 25);

                auto tg_0_yzzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 26);

                auto tg_0_zzzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 27);

                auto tg_0_xxxxx_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx);

                auto tg_0_xxxxy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 1);

                auto tg_0_xxxxz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 2);

                auto tg_0_xxxyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 3);

                auto tg_0_xxxyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 4);

                auto tg_0_xxxzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 5);

                auto tg_0_xxyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 6);

                auto tg_0_xxyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 7);

                auto tg_0_xxyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 8);

                auto tg_0_xxzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 9);

                auto tg_0_xyyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 10);

                auto tg_0_xyyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 11);

                auto tg_0_xyyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 12);

                auto tg_0_xyzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 13);

                auto tg_0_xzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 14);

                auto tg_0_yyyyy_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 15);

                auto tg_0_yyyyz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 16);

                auto tg_0_yyyzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 17);

                auto tg_0_yyzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 18);

                auto tg_0_yzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 19);

                auto tg_0_zzzzz_0 = primBuffer[pidx_g_0_5_m0].data(21 * idx + 20);

                auto tg_0_xxxxx_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx);

                auto tg_0_xxxxy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 1);

                auto tg_0_xxxxz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 2);

                auto tg_0_xxxyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 3);

                auto tg_0_xxxyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 4);

                auto tg_0_xxxzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 5);

                auto tg_0_xxyyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 6);

                auto tg_0_xxyyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 7);

                auto tg_0_xxyzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 8);

                auto tg_0_xxzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 9);

                auto tg_0_xyyyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 10);

                auto tg_0_xyyyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 11);

                auto tg_0_xyyzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 12);

                auto tg_0_xyzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 13);

                auto tg_0_xzzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 14);

                auto tg_0_yyyyy_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 15);

                auto tg_0_yyyyz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 16);

                auto tg_0_yyyzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 17);

                auto tg_0_yyzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 18);

                auto tg_0_yzzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 19);

                auto tg_0_zzzzz_1 = primBuffer[pidx_g_0_5_m1].data(21 * idx + 20);

                // set up pointers to integrals

                auto tg_0_xxxxxxx_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx);

                auto tg_0_xxxxxxy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 1);

                auto tg_0_xxxxxxz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 2);

                auto tg_0_xxxxxyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 3);

                auto tg_0_xxxxxyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 4);

                auto tg_0_xxxxxzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 5);

                auto tg_0_xxxxyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 6);

                auto tg_0_xxxxyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 7);

                auto tg_0_xxxxyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 8);

                auto tg_0_xxxxzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 9);

                auto tg_0_xxxyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 10);

                auto tg_0_xxxyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 11);

                auto tg_0_xxxyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 12);

                auto tg_0_xxxyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 13);

                auto tg_0_xxxzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 14);

                auto tg_0_xxyyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 15);

                auto tg_0_xxyyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 16);

                auto tg_0_xxyyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 17);

                auto tg_0_xxyyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 18);

                auto tg_0_xxyzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 19);

                auto tg_0_xxzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 20);

                auto tg_0_xyyyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 21);

                auto tg_0_xyyyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 22);

                auto tg_0_xyyyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 23);

                auto tg_0_xyyyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 24);

                auto tg_0_xyyzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 25);

                auto tg_0_xyzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 26);

                auto tg_0_xzzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 27);

                auto tg_0_yyyyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 28);

                auto tg_0_yyyyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 29);

                auto tg_0_yyyyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 30);

                auto tg_0_yyyyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 31);

                auto tg_0_yyyzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 32);

                auto tg_0_yyzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 33);

                auto tg_0_yzzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 34);

                auto tg_0_zzzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 35);

                #pragma omp simd aligned(fn, fzb, qd_x, qd_y, qd_z, tg_0_xxxxx_0, tg_0_xxxxx_1, tg_0_xxxxxx_0, \
                                         tg_0_xxxxxx_1, tg_0_xxxxxxx_0, tg_0_xxxxxxy_0, tg_0_xxxxxxz_0, tg_0_xxxxxy_0, \
                                         tg_0_xxxxxy_1, tg_0_xxxxxyy_0, tg_0_xxxxxyz_0, tg_0_xxxxxz_0, tg_0_xxxxxz_1, \
                                         tg_0_xxxxxzz_0, tg_0_xxxxy_0, tg_0_xxxxy_1, tg_0_xxxxyy_0, tg_0_xxxxyy_1, \
                                         tg_0_xxxxyyy_0, tg_0_xxxxyyz_0, tg_0_xxxxyz_0, tg_0_xxxxyz_1, tg_0_xxxxyzz_0, \
                                         tg_0_xxxxz_0, tg_0_xxxxz_1, tg_0_xxxxzz_0, tg_0_xxxxzz_1, tg_0_xxxxzzz_0, \
                                         tg_0_xxxyy_0, tg_0_xxxyy_1, tg_0_xxxyyy_0, tg_0_xxxyyy_1, tg_0_xxxyyyy_0, \
                                         tg_0_xxxyyyz_0, tg_0_xxxyyz_0, tg_0_xxxyyz_1, tg_0_xxxyyzz_0, tg_0_xxxyz_0, \
                                         tg_0_xxxyz_1, tg_0_xxxyzz_0, tg_0_xxxyzz_1, tg_0_xxxyzzz_0, tg_0_xxxzz_0, \
                                         tg_0_xxxzz_1, tg_0_xxxzzz_0, tg_0_xxxzzz_1, tg_0_xxxzzzz_0, tg_0_xxyyy_0, \
                                         tg_0_xxyyy_1, tg_0_xxyyyy_0, tg_0_xxyyyy_1, tg_0_xxyyyyy_0, tg_0_xxyyyyz_0, \
                                         tg_0_xxyyyz_0, tg_0_xxyyyz_1, tg_0_xxyyyzz_0, tg_0_xxyyz_0, tg_0_xxyyz_1, \
                                         tg_0_xxyyzz_0, tg_0_xxyyzz_1, tg_0_xxyyzzz_0, tg_0_xxyzz_0, tg_0_xxyzz_1, \
                                         tg_0_xxyzzz_0, tg_0_xxyzzz_1, tg_0_xxyzzzz_0, tg_0_xxzzz_0, tg_0_xxzzz_1, \
                                         tg_0_xxzzzz_0, tg_0_xxzzzz_1, tg_0_xxzzzzz_0, tg_0_xyyyy_0, tg_0_xyyyy_1, \
                                         tg_0_xyyyyy_0, tg_0_xyyyyy_1, tg_0_xyyyyyy_0, tg_0_xyyyyyz_0, tg_0_xyyyyz_0, \
                                         tg_0_xyyyyz_1, tg_0_xyyyyzz_0, tg_0_xyyyz_0, tg_0_xyyyz_1, tg_0_xyyyzz_0, \
                                         tg_0_xyyyzz_1, tg_0_xyyyzzz_0, tg_0_xyyzz_0, tg_0_xyyzz_1, tg_0_xyyzzz_0, \
                                         tg_0_xyyzzz_1, tg_0_xyyzzzz_0, tg_0_xyzzz_0, tg_0_xyzzz_1, tg_0_xyzzzz_0, \
                                         tg_0_xyzzzz_1, tg_0_xyzzzzz_0, tg_0_xzzzz_0, tg_0_xzzzz_1, tg_0_xzzzzz_0, \
                                         tg_0_xzzzzz_1, tg_0_xzzzzzz_0, tg_0_yyyyy_0, tg_0_yyyyy_1, tg_0_yyyyyy_0, \
                                         tg_0_yyyyyy_1, tg_0_yyyyyyy_0, tg_0_yyyyyyz_0, tg_0_yyyyyz_0, tg_0_yyyyyz_1, \
                                         tg_0_yyyyyzz_0, tg_0_yyyyz_0, tg_0_yyyyz_1, tg_0_yyyyzz_0, tg_0_yyyyzz_1, \
                                         tg_0_yyyyzzz_0, tg_0_yyyzz_0, tg_0_yyyzz_1, tg_0_yyyzzz_0, tg_0_yyyzzz_1, \
                                         tg_0_yyyzzzz_0, tg_0_yyzzz_0, tg_0_yyzzz_1, tg_0_yyzzzz_0, tg_0_yyzzzz_1, \
                                         tg_0_yyzzzzz_0, tg_0_yzzzz_0, tg_0_yzzzz_1, tg_0_yzzzzz_0, tg_0_yzzzzz_1, \
                                         tg_0_yzzzzzz_0, tg_0_zzzzz_0, tg_0_zzzzz_1, tg_0_zzzzzz_0, tg_0_zzzzzz_1, \
                                         tg_0_zzzzzzz_0, wq_x, wq_y, wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fn = fn[j];

                    double fl1_fzb = fzb[j];
                    
                    double fr0 = qd_x[j]; double fr1 = wq_x[j];

                    tg_0_xxxxxxx_0[j] = fr0 * tg_0_xxxxxx_0[j] + fr1 * tg_0_xxxxxx_1[j] + 3.0 * fl1_fn * (tg_0_xxxxx_0[j] - fl1_fzb * tg_0_xxxxx_1[j]);

                    tg_0_xxxxxxy_0[j] = fr0 * tg_0_xxxxxy_0[j] + fr1 * tg_0_xxxxxy_1[j] + 2.5 * fl1_fn * (tg_0_xxxxy_0[j] - fl1_fzb * tg_0_xxxxy_1[j]);

                    tg_0_xxxxxxz_0[j] = fr0 * tg_0_xxxxxz_0[j] + fr1 * tg_0_xxxxxz_1[j] + 2.5 * fl1_fn * (tg_0_xxxxz_0[j] - fl1_fzb * tg_0_xxxxz_1[j]);

                    tg_0_xxxxxyy_0[j] = fr0 * tg_0_xxxxyy_0[j] + fr1 * tg_0_xxxxyy_1[j] + 2.0 * fl1_fn * (tg_0_xxxyy_0[j] - fl1_fzb * tg_0_xxxyy_1[j]);

                    tg_0_xxxxxyz_0[j] = fr0 * tg_0_xxxxyz_0[j] + fr1 * tg_0_xxxxyz_1[j] + 2.0 * fl1_fn * (tg_0_xxxyz_0[j] - fl1_fzb * tg_0_xxxyz_1[j]);

                    tg_0_xxxxxzz_0[j] = fr0 * tg_0_xxxxzz_0[j] + fr1 * tg_0_xxxxzz_1[j] + 2.0 * fl1_fn * (tg_0_xxxzz_0[j] - fl1_fzb * tg_0_xxxzz_1[j]);

                    tg_0_xxxxyyy_0[j] = fr0 * tg_0_xxxyyy_0[j] + fr1 * tg_0_xxxyyy_1[j] + 1.5 * fl1_fn * (tg_0_xxyyy_0[j] - fl1_fzb * tg_0_xxyyy_1[j]);

                    tg_0_xxxxyyz_0[j] = fr0 * tg_0_xxxyyz_0[j] + fr1 * tg_0_xxxyyz_1[j] + 1.5 * fl1_fn * (tg_0_xxyyz_0[j] - fl1_fzb * tg_0_xxyyz_1[j]);

                    tg_0_xxxxyzz_0[j] = fr0 * tg_0_xxxyzz_0[j] + fr1 * tg_0_xxxyzz_1[j] + 1.5 * fl1_fn * (tg_0_xxyzz_0[j] - fl1_fzb * tg_0_xxyzz_1[j]);

                    tg_0_xxxxzzz_0[j] = fr0 * tg_0_xxxzzz_0[j] + fr1 * tg_0_xxxzzz_1[j] + 1.5 * fl1_fn * (tg_0_xxzzz_0[j] - fl1_fzb * tg_0_xxzzz_1[j]);

                    tg_0_xxxyyyy_0[j] = fr0 * tg_0_xxyyyy_0[j] + fr1 * tg_0_xxyyyy_1[j] + fl1_fn * (tg_0_xyyyy_0[j] - fl1_fzb * tg_0_xyyyy_1[j]);

                    tg_0_xxxyyyz_0[j] = fr0 * tg_0_xxyyyz_0[j] + fr1 * tg_0_xxyyyz_1[j] + fl1_fn * (tg_0_xyyyz_0[j] - fl1_fzb * tg_0_xyyyz_1[j]);

                    tg_0_xxxyyzz_0[j] = fr0 * tg_0_xxyyzz_0[j] + fr1 * tg_0_xxyyzz_1[j] + fl1_fn * (tg_0_xyyzz_0[j] - fl1_fzb * tg_0_xyyzz_1[j]);

                    tg_0_xxxyzzz_0[j] = fr0 * tg_0_xxyzzz_0[j] + fr1 * tg_0_xxyzzz_1[j] + fl1_fn * (tg_0_xyzzz_0[j] - fl1_fzb * tg_0_xyzzz_1[j]);

                    tg_0_xxxzzzz_0[j] = fr0 * tg_0_xxzzzz_0[j] + fr1 * tg_0_xxzzzz_1[j] + fl1_fn * (tg_0_xzzzz_0[j] - fl1_fzb * tg_0_xzzzz_1[j]);

                    tg_0_xxyyyyy_0[j] = fr0 * tg_0_xyyyyy_0[j] + fr1 * tg_0_xyyyyy_1[j] + 0.5 * fl1_fn * (tg_0_yyyyy_0[j] - fl1_fzb * tg_0_yyyyy_1[j]);

                    tg_0_xxyyyyz_0[j] = fr0 * tg_0_xyyyyz_0[j] + fr1 * tg_0_xyyyyz_1[j] + 0.5 * fl1_fn * (tg_0_yyyyz_0[j] - fl1_fzb * tg_0_yyyyz_1[j]);

                    tg_0_xxyyyzz_0[j] = fr0 * tg_0_xyyyzz_0[j] + fr1 * tg_0_xyyyzz_1[j] + 0.5 * fl1_fn * (tg_0_yyyzz_0[j] - fl1_fzb * tg_0_yyyzz_1[j]);

                    tg_0_xxyyzzz_0[j] = fr0 * tg_0_xyyzzz_0[j] + fr1 * tg_0_xyyzzz_1[j] + 0.5 * fl1_fn * (tg_0_yyzzz_0[j] - fl1_fzb * tg_0_yyzzz_1[j]);

                    tg_0_xxyzzzz_0[j] = fr0 * tg_0_xyzzzz_0[j] + fr1 * tg_0_xyzzzz_1[j] + 0.5 * fl1_fn * (tg_0_yzzzz_0[j] - fl1_fzb * tg_0_yzzzz_1[j]);

                    tg_0_xxzzzzz_0[j] = fr0 * tg_0_xzzzzz_0[j] + fr1 * tg_0_xzzzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzzzz_0[j] - fl1_fzb * tg_0_zzzzz_1[j]);

                    tg_0_xyyyyyy_0[j] = fr0 * tg_0_yyyyyy_0[j] + fr1 * tg_0_yyyyyy_1[j];

                    tg_0_xyyyyyz_0[j] = fr0 * tg_0_yyyyyz_0[j] + fr1 * tg_0_yyyyyz_1[j];

                    tg_0_xyyyyzz_0[j] = fr0 * tg_0_yyyyzz_0[j] + fr1 * tg_0_yyyyzz_1[j];

                    tg_0_xyyyzzz_0[j] = fr0 * tg_0_yyyzzz_0[j] + fr1 * tg_0_yyyzzz_1[j];

                    tg_0_xyyzzzz_0[j] = fr0 * tg_0_yyzzzz_0[j] + fr1 * tg_0_yyzzzz_1[j];

                    tg_0_xyzzzzz_0[j] = fr0 * tg_0_yzzzzz_0[j] + fr1 * tg_0_yzzzzz_1[j];

                    tg_0_xzzzzzz_0[j] = fr0 * tg_0_zzzzzz_0[j] + fr1 * tg_0_zzzzzz_1[j];
                    
                    fr0 = qd_y[j]; fr1 = wq_y[j];

                    tg_0_yyyyyyy_0[j] = fr0 * tg_0_yyyyyy_0[j] + fr1 * tg_0_yyyyyy_1[j] + 3.0 * fl1_fn * (tg_0_yyyyy_0[j] - fl1_fzb * tg_0_yyyyy_1[j]);

                    tg_0_yyyyyyz_0[j] = fr0 * tg_0_yyyyyz_0[j] + fr1 * tg_0_yyyyyz_1[j] + 2.5 * fl1_fn * (tg_0_yyyyz_0[j] - fl1_fzb * tg_0_yyyyz_1[j]);

                    tg_0_yyyyyzz_0[j] = fr0 * tg_0_yyyyzz_0[j] + fr1 * tg_0_yyyyzz_1[j] + 2.0 * fl1_fn * (tg_0_yyyzz_0[j] - fl1_fzb * tg_0_yyyzz_1[j]);

                    tg_0_yyyyzzz_0[j] = fr0 * tg_0_yyyzzz_0[j] + fr1 * tg_0_yyyzzz_1[j] + 1.5 * fl1_fn * (tg_0_yyzzz_0[j] - fl1_fzb * tg_0_yyzzz_1[j]);

                    tg_0_yyyzzzz_0[j] = fr0 * tg_0_yyzzzz_0[j] + fr1 * tg_0_yyzzzz_1[j] + fl1_fn * (tg_0_yzzzz_0[j] - fl1_fzb * tg_0_yzzzz_1[j]);

                    tg_0_yyzzzzz_0[j] = fr0 * tg_0_yzzzzz_0[j] + fr1 * tg_0_yzzzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzzzz_0[j] - fl1_fzb * tg_0_zzzzz_1[j]);

                    tg_0_yzzzzzz_0[j] = fr0 * tg_0_zzzzzz_0[j] + fr1 * tg_0_zzzzzz_1[j];

                    tg_0_zzzzzzz_0[j] = qd_z[j] * tg_0_zzzzzz_0[j] + wq_z[j] * tg_0_zzzzzz_1[j] + 3.0 * fl1_fn * (tg_0_zzzzz_0[j] - fl1_fzb * tg_0_zzzzz_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSKSS(      CMemBlock2D<double>* primBuffer,
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
                                             {7, -1, -1, -1},
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_7_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {7, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_7_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_6_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {6, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_6_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {6, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_5_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {5, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_5_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {5, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_xxxxxx_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx);

                auto tg_xxxxxy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 1);

                auto tg_xxxxxz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 2);

                auto tg_xxxxyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 3);

                auto tg_xxxxyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 4);

                auto tg_xxxxzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 5);

                auto tg_xxxyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 6);

                auto tg_xxxyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 7);

                auto tg_xxxyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 8);

                auto tg_xxxzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 9);

                auto tg_xxyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 10);

                auto tg_xxyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 11);

                auto tg_xxyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 12);

                auto tg_xxyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 13);

                auto tg_xxzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 14);

                auto tg_xyyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 15);

                auto tg_xyyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 16);

                auto tg_xyyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 17);

                auto tg_xyyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 18);

                auto tg_xyzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 19);

                auto tg_xzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 20);

                auto tg_yyyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 21);

                auto tg_yyyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 22);

                auto tg_yyyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 23);

                auto tg_yyyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 24);

                auto tg_yyzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 25);

                auto tg_yzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 26);

                auto tg_zzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 27);

                auto tg_xxxxxx_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx);

                auto tg_xxxxxy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 1);

                auto tg_xxxxxz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 2);

                auto tg_xxxxyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 3);

                auto tg_xxxxyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 4);

                auto tg_xxxxzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 5);

                auto tg_xxxyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 6);

                auto tg_xxxyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 7);

                auto tg_xxxyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 8);

                auto tg_xxxzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 9);

                auto tg_xxyyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 10);

                auto tg_xxyyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 11);

                auto tg_xxyyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 12);

                auto tg_xxyzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 13);

                auto tg_xxzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 14);

                auto tg_xyyyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 15);

                auto tg_xyyyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 16);

                auto tg_xyyyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 17);

                auto tg_xyyzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 18);

                auto tg_xyzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 19);

                auto tg_xzzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 20);

                auto tg_yyyyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 21);

                auto tg_yyyyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 22);

                auto tg_yyyyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 23);

                auto tg_yyyzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 24);

                auto tg_yyzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 25);

                auto tg_yzzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 26);

                auto tg_zzzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 27);

                auto tg_xxxxx_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx);

                auto tg_xxxxy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 1);

                auto tg_xxxxz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 2);

                auto tg_xxxyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 3);

                auto tg_xxxyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 4);

                auto tg_xxxzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 5);

                auto tg_xxyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 6);

                auto tg_xxyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 7);

                auto tg_xxyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 8);

                auto tg_xxzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 9);

                auto tg_xyyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 10);

                auto tg_xyyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 11);

                auto tg_xyyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 12);

                auto tg_xyzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 13);

                auto tg_xzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 14);

                auto tg_yyyyy_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 15);

                auto tg_yyyyz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 16);

                auto tg_yyyzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 17);

                auto tg_yyzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 18);

                auto tg_yzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 19);

                auto tg_zzzzz_0_0 = primBuffer[pidx_g_5_0_m0].data(21 * idx + 20);

                auto tg_xxxxx_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx);

                auto tg_xxxxy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 1);

                auto tg_xxxxz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 2);

                auto tg_xxxyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 3);

                auto tg_xxxyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 4);

                auto tg_xxxzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 5);

                auto tg_xxyyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 6);

                auto tg_xxyyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 7);

                auto tg_xxyzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 8);

                auto tg_xxzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 9);

                auto tg_xyyyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 10);

                auto tg_xyyyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 11);

                auto tg_xyyzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 12);

                auto tg_xyzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 13);

                auto tg_xzzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 14);

                auto tg_yyyyy_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 15);

                auto tg_yyyyz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 16);

                auto tg_yyyzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 17);

                auto tg_yyzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 18);

                auto tg_yzzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 19);

                auto tg_zzzzz_0_1 = primBuffer[pidx_g_5_0_m1].data(21 * idx + 20);

                // set up pointers to integrals

                auto tg_xxxxxxx_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx);

                auto tg_xxxxxxy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 1);

                auto tg_xxxxxxz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 2);

                auto tg_xxxxxyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 3);

                auto tg_xxxxxyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 4);

                auto tg_xxxxxzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 5);

                auto tg_xxxxyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 6);

                auto tg_xxxxyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 7);

                auto tg_xxxxyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 8);

                auto tg_xxxxzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 9);

                auto tg_xxxyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 10);

                auto tg_xxxyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 11);

                auto tg_xxxyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 12);

                auto tg_xxxyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 13);

                auto tg_xxxzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 14);

                auto tg_xxyyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 15);

                auto tg_xxyyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 16);

                auto tg_xxyyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 17);

                auto tg_xxyyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 18);

                auto tg_xxyzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 19);

                auto tg_xxzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 20);

                auto tg_xyyyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 21);

                auto tg_xyyyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 22);

                auto tg_xyyyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 23);

                auto tg_xyyyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 24);

                auto tg_xyyzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 25);

                auto tg_xyzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 26);

                auto tg_xzzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 27);

                auto tg_yyyyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 28);

                auto tg_yyyyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 29);

                auto tg_yyyyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 30);

                auto tg_yyyyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 31);

                auto tg_yyyzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 32);

                auto tg_yyzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 33);

                auto tg_yzzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 34);

                auto tg_zzzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 35);

                #pragma omp simd aligned(fza, tg_xxxxx_0_0, tg_xxxxx_0_1, tg_xxxxxx_0_0, tg_xxxxxx_0_1, \
                                         tg_xxxxxxx_0_0, tg_xxxxxxy_0_0, tg_xxxxxxz_0_0, tg_xxxxxy_0_0, tg_xxxxxy_0_1, \
                                         tg_xxxxxyy_0_0, tg_xxxxxyz_0_0, tg_xxxxxz_0_0, tg_xxxxxz_0_1, tg_xxxxxzz_0_0, \
                                         tg_xxxxy_0_0, tg_xxxxy_0_1, tg_xxxxyy_0_0, tg_xxxxyy_0_1, tg_xxxxyyy_0_0, \
                                         tg_xxxxyyz_0_0, tg_xxxxyz_0_0, tg_xxxxyz_0_1, tg_xxxxyzz_0_0, tg_xxxxz_0_0, \
                                         tg_xxxxz_0_1, tg_xxxxzz_0_0, tg_xxxxzz_0_1, tg_xxxxzzz_0_0, tg_xxxyy_0_0, \
                                         tg_xxxyy_0_1, tg_xxxyyy_0_0, tg_xxxyyy_0_1, tg_xxxyyyy_0_0, tg_xxxyyyz_0_0, \
                                         tg_xxxyyz_0_0, tg_xxxyyz_0_1, tg_xxxyyzz_0_0, tg_xxxyz_0_0, tg_xxxyz_0_1, \
                                         tg_xxxyzz_0_0, tg_xxxyzz_0_1, tg_xxxyzzz_0_0, tg_xxxzz_0_0, tg_xxxzz_0_1, \
                                         tg_xxxzzz_0_0, tg_xxxzzz_0_1, tg_xxxzzzz_0_0, tg_xxyyy_0_0, tg_xxyyy_0_1, \
                                         tg_xxyyyy_0_0, tg_xxyyyy_0_1, tg_xxyyyyy_0_0, tg_xxyyyyz_0_0, tg_xxyyyz_0_0, \
                                         tg_xxyyyz_0_1, tg_xxyyyzz_0_0, tg_xxyyz_0_0, tg_xxyyz_0_1, tg_xxyyzz_0_0, \
                                         tg_xxyyzz_0_1, tg_xxyyzzz_0_0, tg_xxyzz_0_0, tg_xxyzz_0_1, tg_xxyzzz_0_0, \
                                         tg_xxyzzz_0_1, tg_xxyzzzz_0_0, tg_xxzzz_0_0, tg_xxzzz_0_1, tg_xxzzzz_0_0, \
                                         tg_xxzzzz_0_1, tg_xxzzzzz_0_0, tg_xyyyy_0_0, tg_xyyyy_0_1, tg_xyyyyy_0_0, \
                                         tg_xyyyyy_0_1, tg_xyyyyyy_0_0, tg_xyyyyyz_0_0, tg_xyyyyz_0_0, tg_xyyyyz_0_1, \
                                         tg_xyyyyzz_0_0, tg_xyyyz_0_0, tg_xyyyz_0_1, tg_xyyyzz_0_0, tg_xyyyzz_0_1, \
                                         tg_xyyyzzz_0_0, tg_xyyzz_0_0, tg_xyyzz_0_1, tg_xyyzzz_0_0, tg_xyyzzz_0_1, \
                                         tg_xyyzzzz_0_0, tg_xyzzz_0_0, tg_xyzzz_0_1, tg_xyzzzz_0_0, tg_xyzzzz_0_1, \
                                         tg_xyzzzzz_0_0, tg_xzzzz_0_0, tg_xzzzz_0_1, tg_xzzzzz_0_0, tg_xzzzzz_0_1, \
                                         tg_xzzzzzz_0_0, tg_yyyyy_0_0, tg_yyyyy_0_1, tg_yyyyyy_0_0, tg_yyyyyy_0_1, \
                                         tg_yyyyyyy_0_0, tg_yyyyyyz_0_0, tg_yyyyyz_0_0, tg_yyyyyz_0_1, tg_yyyyyzz_0_0, \
                                         tg_yyyyz_0_0, tg_yyyyz_0_1, tg_yyyyzz_0_0, tg_yyyyzz_0_1, tg_yyyyzzz_0_0, \
                                         tg_yyyzz_0_0, tg_yyyzz_0_1, tg_yyyzzz_0_0, tg_yyyzzz_0_1, tg_yyyzzzz_0_0, \
                                         tg_yyzzz_0_0, tg_yyzzz_0_1, tg_yyzzzz_0_0, tg_yyzzzz_0_1, tg_yyzzzzz_0_0, \
                                         tg_yzzzz_0_0, tg_yzzzz_0_1, tg_yzzzzz_0_0, tg_yzzzzz_0_1, tg_yzzzzzz_0_0, \
                                         tg_zzzzz_0_0, tg_zzzzz_0_1, tg_zzzzzz_0_0, tg_zzzzzz_0_1, tg_zzzzzzz_0_0, wp_x, wp_y, \
                                         wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fza = fza[j];
                    
                    double fr1 = wp_x[j];

                    tg_xxxxxxx_0_0[j] = pb_x * tg_xxxxxx_0_0[j] + fr1 * tg_xxxxxx_0_1[j] + 3.0 * fl1_fx * (tg_xxxxx_0_0[j] -  fl1_fza * tg_xxxxx_0_1[j]);

                    tg_xxxxxxy_0_0[j] = pb_x * tg_xxxxxy_0_0[j] + fr1 * tg_xxxxxy_0_1[j] + 2.5 * fl1_fx * (tg_xxxxy_0_0[j] -  fl1_fza * tg_xxxxy_0_1[j]);

                    tg_xxxxxxz_0_0[j] = pb_x * tg_xxxxxz_0_0[j] + fr1 * tg_xxxxxz_0_1[j] + 2.5 * fl1_fx * (tg_xxxxz_0_0[j] -  fl1_fza * tg_xxxxz_0_1[j]);

                    tg_xxxxxyy_0_0[j] = pb_x * tg_xxxxyy_0_0[j] + fr1 * tg_xxxxyy_0_1[j] + 2.0 * fl1_fx * (tg_xxxyy_0_0[j] -  fl1_fza * tg_xxxyy_0_1[j]);

                    tg_xxxxxyz_0_0[j] = pb_x * tg_xxxxyz_0_0[j] + fr1 * tg_xxxxyz_0_1[j] + 2.0 * fl1_fx * (tg_xxxyz_0_0[j] -  fl1_fza * tg_xxxyz_0_1[j]);

                    tg_xxxxxzz_0_0[j] = pb_x * tg_xxxxzz_0_0[j] + fr1 * tg_xxxxzz_0_1[j] + 2.0 * fl1_fx * (tg_xxxzz_0_0[j] -  fl1_fza * tg_xxxzz_0_1[j]);

                    tg_xxxxyyy_0_0[j] = pb_x * tg_xxxyyy_0_0[j] + fr1 * tg_xxxyyy_0_1[j] + 1.5 * fl1_fx * (tg_xxyyy_0_0[j] -  fl1_fza * tg_xxyyy_0_1[j]);

                    tg_xxxxyyz_0_0[j] = pb_x * tg_xxxyyz_0_0[j] + fr1 * tg_xxxyyz_0_1[j] + 1.5 * fl1_fx * (tg_xxyyz_0_0[j] -  fl1_fza * tg_xxyyz_0_1[j]);

                    tg_xxxxyzz_0_0[j] = pb_x * tg_xxxyzz_0_0[j] + fr1 * tg_xxxyzz_0_1[j] + 1.5 * fl1_fx * (tg_xxyzz_0_0[j] -  fl1_fza * tg_xxyzz_0_1[j]);

                    tg_xxxxzzz_0_0[j] = pb_x * tg_xxxzzz_0_0[j] + fr1 * tg_xxxzzz_0_1[j] + 1.5 * fl1_fx * (tg_xxzzz_0_0[j] -  fl1_fza * tg_xxzzz_0_1[j]);

                    tg_xxxyyyy_0_0[j] = pb_x * tg_xxyyyy_0_0[j] + fr1 * tg_xxyyyy_0_1[j] + fl1_fx * (tg_xyyyy_0_0[j] -  fl1_fza * tg_xyyyy_0_1[j]);

                    tg_xxxyyyz_0_0[j] = pb_x * tg_xxyyyz_0_0[j] + fr1 * tg_xxyyyz_0_1[j] + fl1_fx * (tg_xyyyz_0_0[j] -  fl1_fza * tg_xyyyz_0_1[j]);

                    tg_xxxyyzz_0_0[j] = pb_x * tg_xxyyzz_0_0[j] + fr1 * tg_xxyyzz_0_1[j] + fl1_fx * (tg_xyyzz_0_0[j] -  fl1_fza * tg_xyyzz_0_1[j]);

                    tg_xxxyzzz_0_0[j] = pb_x * tg_xxyzzz_0_0[j] + fr1 * tg_xxyzzz_0_1[j] + fl1_fx * (tg_xyzzz_0_0[j] -  fl1_fza * tg_xyzzz_0_1[j]);

                    tg_xxxzzzz_0_0[j] = pb_x * tg_xxzzzz_0_0[j] + fr1 * tg_xxzzzz_0_1[j] + fl1_fx * (tg_xzzzz_0_0[j] -  fl1_fza * tg_xzzzz_0_1[j]);

                    tg_xxyyyyy_0_0[j] = pb_x * tg_xyyyyy_0_0[j] + fr1 * tg_xyyyyy_0_1[j] + 0.5 * fl1_fx * (tg_yyyyy_0_0[j] -  fl1_fza * tg_yyyyy_0_1[j]);

                    tg_xxyyyyz_0_0[j] = pb_x * tg_xyyyyz_0_0[j] + fr1 * tg_xyyyyz_0_1[j] + 0.5 * fl1_fx * (tg_yyyyz_0_0[j] -  fl1_fza * tg_yyyyz_0_1[j]);

                    tg_xxyyyzz_0_0[j] = pb_x * tg_xyyyzz_0_0[j] + fr1 * tg_xyyyzz_0_1[j] + 0.5 * fl1_fx * (tg_yyyzz_0_0[j] -  fl1_fza * tg_yyyzz_0_1[j]);

                    tg_xxyyzzz_0_0[j] = pb_x * tg_xyyzzz_0_0[j] + fr1 * tg_xyyzzz_0_1[j] + 0.5 * fl1_fx * (tg_yyzzz_0_0[j] -  fl1_fza * tg_yyzzz_0_1[j]);

                    tg_xxyzzzz_0_0[j] = pb_x * tg_xyzzzz_0_0[j] + fr1 * tg_xyzzzz_0_1[j] + 0.5 * fl1_fx * (tg_yzzzz_0_0[j] -  fl1_fza * tg_yzzzz_0_1[j]);

                    tg_xxzzzzz_0_0[j] = pb_x * tg_xzzzzz_0_0[j] + fr1 * tg_xzzzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzzzz_0_0[j] -  fl1_fza * tg_zzzzz_0_1[j]);

                    tg_xyyyyyy_0_0[j] = pb_x * tg_yyyyyy_0_0[j] + fr1 * tg_yyyyyy_0_1[j];

                    tg_xyyyyyz_0_0[j] = pb_x * tg_yyyyyz_0_0[j] + fr1 * tg_yyyyyz_0_1[j];

                    tg_xyyyyzz_0_0[j] = pb_x * tg_yyyyzz_0_0[j] + fr1 * tg_yyyyzz_0_1[j];

                    tg_xyyyzzz_0_0[j] = pb_x * tg_yyyzzz_0_0[j] + fr1 * tg_yyyzzz_0_1[j];

                    tg_xyyzzzz_0_0[j] = pb_x * tg_yyzzzz_0_0[j] + fr1 * tg_yyzzzz_0_1[j];

                    tg_xyzzzzz_0_0[j] = pb_x * tg_yzzzzz_0_0[j] + fr1 * tg_yzzzzz_0_1[j];

                    tg_xzzzzzz_0_0[j] = pb_x * tg_zzzzzz_0_0[j] + fr1 * tg_zzzzzz_0_1[j];
                                                                                                           
                    fr1 = wp_y[j];

                    tg_yyyyyyy_0_0[j] = pb_y * tg_yyyyyy_0_0[j] + fr1 * tg_yyyyyy_0_1[j] + 3.0 * fl1_fx * (tg_yyyyy_0_0[j] -  fl1_fza * tg_yyyyy_0_1[j]);

                    tg_yyyyyyz_0_0[j] = pb_y * tg_yyyyyz_0_0[j] + fr1 * tg_yyyyyz_0_1[j] + 2.5 * fl1_fx * (tg_yyyyz_0_0[j] -  fl1_fza * tg_yyyyz_0_1[j]);

                    tg_yyyyyzz_0_0[j] = pb_y * tg_yyyyzz_0_0[j] + fr1 * tg_yyyyzz_0_1[j] + 2.0 * fl1_fx * (tg_yyyzz_0_0[j] -  fl1_fza * tg_yyyzz_0_1[j]);

                    tg_yyyyzzz_0_0[j] = pb_y * tg_yyyzzz_0_0[j] + fr1 * tg_yyyzzz_0_1[j] + 1.5 * fl1_fx * (tg_yyzzz_0_0[j] -  fl1_fza * tg_yyzzz_0_1[j]);

                    tg_yyyzzzz_0_0[j] = pb_y * tg_yyzzzz_0_0[j] + fr1 * tg_yyzzzz_0_1[j] + fl1_fx * (tg_yzzzz_0_0[j] -  fl1_fza * tg_yzzzz_0_1[j]);

                    tg_yyzzzzz_0_0[j] = pb_y * tg_yzzzzz_0_0[j] + fr1 * tg_yzzzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzzzz_0_0[j] -  fl1_fza * tg_zzzzz_0_1[j]);

                    tg_yzzzzzz_0_0[j] = pb_y * tg_zzzzzz_0_0[j] + fr1 * tg_zzzzzz_0_1[j];

                    tg_zzzzzzz_0_0[j] = pb_z * tg_zzzzzz_0_0[j] + wp_z[j] * tg_zzzzzz_0_1[j] + 3.0 * fl1_fx * (tg_zzzzz_0_0[j] -  fl1_fza * tg_zzzzz_0_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSSSL(      CMemBlock2D<double>* primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& wqDistances,
                                 const CGtoPairsBlock&      braGtoPairsBlock,
                                 const CGtoPairsBlock&      ketGtoPairsBlock,
                                 const int32_t              nKetPrimPairs,
                                 const int32_t              iContrPair)
    {
        // set up pointers to primitives data on bra side

        auto spos = braGtoPairsBlock.getStartPositions();

        auto epos = braGtoPairsBlock.getEndPositions();

        // set up pointers to tensor of distance R(QD) = Q - D

        auto qd_x = ketGtoPairsBlock.getDistancesPBX();

        auto qd_y = ketGtoPairsBlock.getDistancesPBY();

        auto qd_z = ketGtoPairsBlock.getDistancesPBZ();

        // set up pointers to common Obara-Saika factors

        auto fn = ketGtoPairsBlock.getFactorsOneOverXi();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electron Repulsion"},
                                             {0, -1, -1, -1},
                                             {8, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_0_8_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {8, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_0_8_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_0_7_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {7, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_7_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {7, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_0_6_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {6, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_0_6_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {0, -1, -1, -1}, {6, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fzb = osFactors.data(4 * idx + 3);

                // set up pointers to tensors product of distances R(WQ) = W - Q

                auto wq_x = wqDistances.data(3 * idx);

                auto wq_y = wqDistances.data(3 * idx + 1);

                auto wq_z = wqDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_0_xxxxxxx_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx);

                auto tg_0_xxxxxxy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 1);

                auto tg_0_xxxxxxz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 2);

                auto tg_0_xxxxxyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 3);

                auto tg_0_xxxxxyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 4);

                auto tg_0_xxxxxzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 5);

                auto tg_0_xxxxyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 6);

                auto tg_0_xxxxyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 7);

                auto tg_0_xxxxyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 8);

                auto tg_0_xxxxzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 9);

                auto tg_0_xxxyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 10);

                auto tg_0_xxxyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 11);

                auto tg_0_xxxyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 12);

                auto tg_0_xxxyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 13);

                auto tg_0_xxxzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 14);

                auto tg_0_xxyyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 15);

                auto tg_0_xxyyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 16);

                auto tg_0_xxyyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 17);

                auto tg_0_xxyyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 18);

                auto tg_0_xxyzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 19);

                auto tg_0_xxzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 20);

                auto tg_0_xyyyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 21);

                auto tg_0_xyyyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 22);

                auto tg_0_xyyyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 23);

                auto tg_0_xyyyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 24);

                auto tg_0_xyyzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 25);

                auto tg_0_xyzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 26);

                auto tg_0_xzzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 27);

                auto tg_0_yyyyyyy_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 28);

                auto tg_0_yyyyyyz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 29);

                auto tg_0_yyyyyzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 30);

                auto tg_0_yyyyzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 31);

                auto tg_0_yyyzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 32);

                auto tg_0_yyzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 33);

                auto tg_0_yzzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 34);

                auto tg_0_zzzzzzz_0 = primBuffer[pidx_g_0_7_m0].data(36 * idx + 35);

                auto tg_0_xxxxxxx_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx);

                auto tg_0_xxxxxxy_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 1);

                auto tg_0_xxxxxxz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 2);

                auto tg_0_xxxxxyy_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 3);

                auto tg_0_xxxxxyz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 4);

                auto tg_0_xxxxxzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 5);

                auto tg_0_xxxxyyy_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 6);

                auto tg_0_xxxxyyz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 7);

                auto tg_0_xxxxyzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 8);

                auto tg_0_xxxxzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 9);

                auto tg_0_xxxyyyy_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 10);

                auto tg_0_xxxyyyz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 11);

                auto tg_0_xxxyyzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 12);

                auto tg_0_xxxyzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 13);

                auto tg_0_xxxzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 14);

                auto tg_0_xxyyyyy_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 15);

                auto tg_0_xxyyyyz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 16);

                auto tg_0_xxyyyzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 17);

                auto tg_0_xxyyzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 18);

                auto tg_0_xxyzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 19);

                auto tg_0_xxzzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 20);

                auto tg_0_xyyyyyy_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 21);

                auto tg_0_xyyyyyz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 22);

                auto tg_0_xyyyyzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 23);

                auto tg_0_xyyyzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 24);

                auto tg_0_xyyzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 25);

                auto tg_0_xyzzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 26);

                auto tg_0_xzzzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 27);

                auto tg_0_yyyyyyy_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 28);

                auto tg_0_yyyyyyz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 29);

                auto tg_0_yyyyyzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 30);

                auto tg_0_yyyyzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 31);

                auto tg_0_yyyzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 32);

                auto tg_0_yyzzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 33);

                auto tg_0_yzzzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 34);

                auto tg_0_zzzzzzz_1 = primBuffer[pidx_g_0_7_m1].data(36 * idx + 35);

                auto tg_0_xxxxxx_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx);

                auto tg_0_xxxxxy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 1);

                auto tg_0_xxxxxz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 2);

                auto tg_0_xxxxyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 3);

                auto tg_0_xxxxyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 4);

                auto tg_0_xxxxzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 5);

                auto tg_0_xxxyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 6);

                auto tg_0_xxxyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 7);

                auto tg_0_xxxyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 8);

                auto tg_0_xxxzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 9);

                auto tg_0_xxyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 10);

                auto tg_0_xxyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 11);

                auto tg_0_xxyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 12);

                auto tg_0_xxyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 13);

                auto tg_0_xxzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 14);

                auto tg_0_xyyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 15);

                auto tg_0_xyyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 16);

                auto tg_0_xyyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 17);

                auto tg_0_xyyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 18);

                auto tg_0_xyzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 19);

                auto tg_0_xzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 20);

                auto tg_0_yyyyyy_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 21);

                auto tg_0_yyyyyz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 22);

                auto tg_0_yyyyzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 23);

                auto tg_0_yyyzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 24);

                auto tg_0_yyzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 25);

                auto tg_0_yzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 26);

                auto tg_0_zzzzzz_0 = primBuffer[pidx_g_0_6_m0].data(28 * idx + 27);

                auto tg_0_xxxxxx_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx);

                auto tg_0_xxxxxy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 1);

                auto tg_0_xxxxxz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 2);

                auto tg_0_xxxxyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 3);

                auto tg_0_xxxxyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 4);

                auto tg_0_xxxxzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 5);

                auto tg_0_xxxyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 6);

                auto tg_0_xxxyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 7);

                auto tg_0_xxxyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 8);

                auto tg_0_xxxzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 9);

                auto tg_0_xxyyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 10);

                auto tg_0_xxyyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 11);

                auto tg_0_xxyyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 12);

                auto tg_0_xxyzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 13);

                auto tg_0_xxzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 14);

                auto tg_0_xyyyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 15);

                auto tg_0_xyyyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 16);

                auto tg_0_xyyyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 17);

                auto tg_0_xyyzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 18);

                auto tg_0_xyzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 19);

                auto tg_0_xzzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 20);

                auto tg_0_yyyyyy_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 21);

                auto tg_0_yyyyyz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 22);

                auto tg_0_yyyyzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 23);

                auto tg_0_yyyzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 24);

                auto tg_0_yyzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 25);

                auto tg_0_yzzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 26);

                auto tg_0_zzzzzz_1 = primBuffer[pidx_g_0_6_m1].data(28 * idx + 27);

                // set up pointers to integrals

                auto tg_0_xxxxxxxx_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx);

                auto tg_0_xxxxxxxy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 1);

                auto tg_0_xxxxxxxz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 2);

                auto tg_0_xxxxxxyy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 3);

                auto tg_0_xxxxxxyz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 4);

                auto tg_0_xxxxxxzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 5);

                auto tg_0_xxxxxyyy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 6);

                auto tg_0_xxxxxyyz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 7);

                auto tg_0_xxxxxyzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 8);

                auto tg_0_xxxxxzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 9);

                auto tg_0_xxxxyyyy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 10);

                auto tg_0_xxxxyyyz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 11);

                auto tg_0_xxxxyyzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 12);

                auto tg_0_xxxxyzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 13);

                auto tg_0_xxxxzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 14);

                auto tg_0_xxxyyyyy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 15);

                auto tg_0_xxxyyyyz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 16);

                auto tg_0_xxxyyyzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 17);

                auto tg_0_xxxyyzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 18);

                auto tg_0_xxxyzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 19);

                auto tg_0_xxxzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 20);

                auto tg_0_xxyyyyyy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 21);

                auto tg_0_xxyyyyyz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 22);

                auto tg_0_xxyyyyzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 23);

                auto tg_0_xxyyyzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 24);

                auto tg_0_xxyyzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 25);

                auto tg_0_xxyzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 26);

                auto tg_0_xxzzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 27);

                auto tg_0_xyyyyyyy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 28);

                auto tg_0_xyyyyyyz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 29);

                auto tg_0_xyyyyyzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 30);

                auto tg_0_xyyyyzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 31);

                auto tg_0_xyyyzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 32);

                auto tg_0_xyyzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 33);

                auto tg_0_xyzzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 34);

                auto tg_0_xzzzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 35);

                auto tg_0_yyyyyyyy_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 36);

                auto tg_0_yyyyyyyz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 37);

                auto tg_0_yyyyyyzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 38);

                auto tg_0_yyyyyzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 39);

                auto tg_0_yyyyzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 40);

                auto tg_0_yyyzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 41);

                auto tg_0_yyzzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 42);

                auto tg_0_yzzzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 43);

                auto tg_0_zzzzzzzz_0 = primBuffer[pidx_g_0_8_m0].data(45 * idx + 44);

                #pragma omp simd aligned(fn, fzb, qd_x, qd_y, qd_z, tg_0_xxxxxx_0, tg_0_xxxxxx_1, tg_0_xxxxxxx_0, \
                                         tg_0_xxxxxxx_1, tg_0_xxxxxxxx_0, tg_0_xxxxxxxy_0, tg_0_xxxxxxxz_0, tg_0_xxxxxxy_0, \
                                         tg_0_xxxxxxy_1, tg_0_xxxxxxyy_0, tg_0_xxxxxxyz_0, tg_0_xxxxxxz_0, tg_0_xxxxxxz_1, \
                                         tg_0_xxxxxxzz_0, tg_0_xxxxxy_0, tg_0_xxxxxy_1, tg_0_xxxxxyy_0, tg_0_xxxxxyy_1, \
                                         tg_0_xxxxxyyy_0, tg_0_xxxxxyyz_0, tg_0_xxxxxyz_0, tg_0_xxxxxyz_1, tg_0_xxxxxyzz_0, \
                                         tg_0_xxxxxz_0, tg_0_xxxxxz_1, tg_0_xxxxxzz_0, tg_0_xxxxxzz_1, tg_0_xxxxxzzz_0, \
                                         tg_0_xxxxyy_0, tg_0_xxxxyy_1, tg_0_xxxxyyy_0, tg_0_xxxxyyy_1, tg_0_xxxxyyyy_0, \
                                         tg_0_xxxxyyyz_0, tg_0_xxxxyyz_0, tg_0_xxxxyyz_1, tg_0_xxxxyyzz_0, tg_0_xxxxyz_0, \
                                         tg_0_xxxxyz_1, tg_0_xxxxyzz_0, tg_0_xxxxyzz_1, tg_0_xxxxyzzz_0, tg_0_xxxxzz_0, \
                                         tg_0_xxxxzz_1, tg_0_xxxxzzz_0, tg_0_xxxxzzz_1, tg_0_xxxxzzzz_0, tg_0_xxxyyy_0, \
                                         tg_0_xxxyyy_1, tg_0_xxxyyyy_0, tg_0_xxxyyyy_1, tg_0_xxxyyyyy_0, tg_0_xxxyyyyz_0, \
                                         tg_0_xxxyyyz_0, tg_0_xxxyyyz_1, tg_0_xxxyyyzz_0, tg_0_xxxyyz_0, tg_0_xxxyyz_1, \
                                         tg_0_xxxyyzz_0, tg_0_xxxyyzz_1, tg_0_xxxyyzzz_0, tg_0_xxxyzz_0, tg_0_xxxyzz_1, \
                                         tg_0_xxxyzzz_0, tg_0_xxxyzzz_1, tg_0_xxxyzzzz_0, tg_0_xxxzzz_0, tg_0_xxxzzz_1, \
                                         tg_0_xxxzzzz_0, tg_0_xxxzzzz_1, tg_0_xxxzzzzz_0, tg_0_xxyyyy_0, tg_0_xxyyyy_1, \
                                         tg_0_xxyyyyy_0, tg_0_xxyyyyy_1, tg_0_xxyyyyyy_0, tg_0_xxyyyyyz_0, tg_0_xxyyyyz_0, \
                                         tg_0_xxyyyyz_1, tg_0_xxyyyyzz_0, tg_0_xxyyyz_0, tg_0_xxyyyz_1, tg_0_xxyyyzz_0, \
                                         tg_0_xxyyyzz_1, tg_0_xxyyyzzz_0, tg_0_xxyyzz_0, tg_0_xxyyzz_1, tg_0_xxyyzzz_0, \
                                         tg_0_xxyyzzz_1, tg_0_xxyyzzzz_0, tg_0_xxyzzz_0, tg_0_xxyzzz_1, tg_0_xxyzzzz_0, \
                                         tg_0_xxyzzzz_1, tg_0_xxyzzzzz_0, tg_0_xxzzzz_0, tg_0_xxzzzz_1, tg_0_xxzzzzz_0, \
                                         tg_0_xxzzzzz_1, tg_0_xxzzzzzz_0, tg_0_xyyyyy_0, tg_0_xyyyyy_1, tg_0_xyyyyyy_0, \
                                         tg_0_xyyyyyy_1, tg_0_xyyyyyyy_0, tg_0_xyyyyyyz_0, tg_0_xyyyyyz_0, tg_0_xyyyyyz_1, \
                                         tg_0_xyyyyyzz_0, tg_0_xyyyyz_0, tg_0_xyyyyz_1, tg_0_xyyyyzz_0, tg_0_xyyyyzz_1, \
                                         tg_0_xyyyyzzz_0, tg_0_xyyyzz_0, tg_0_xyyyzz_1, tg_0_xyyyzzz_0, tg_0_xyyyzzz_1, \
                                         tg_0_xyyyzzzz_0, tg_0_xyyzzz_0, tg_0_xyyzzz_1, tg_0_xyyzzzz_0, tg_0_xyyzzzz_1, \
                                         tg_0_xyyzzzzz_0, tg_0_xyzzzz_0, tg_0_xyzzzz_1, tg_0_xyzzzzz_0, tg_0_xyzzzzz_1, \
                                         tg_0_xyzzzzzz_0, tg_0_xzzzzz_0, tg_0_xzzzzz_1, tg_0_xzzzzzz_0, tg_0_xzzzzzz_1, \
                                         tg_0_xzzzzzzz_0, tg_0_yyyyyy_0, tg_0_yyyyyy_1, tg_0_yyyyyyy_0, tg_0_yyyyyyy_1, \
                                         tg_0_yyyyyyyy_0, tg_0_yyyyyyyz_0, tg_0_yyyyyyz_0, tg_0_yyyyyyz_1, tg_0_yyyyyyzz_0, \
                                         tg_0_yyyyyz_0, tg_0_yyyyyz_1, tg_0_yyyyyzz_0, tg_0_yyyyyzz_1, tg_0_yyyyyzzz_0, \
                                         tg_0_yyyyzz_0, tg_0_yyyyzz_1, tg_0_yyyyzzz_0, tg_0_yyyyzzz_1, tg_0_yyyyzzzz_0, \
                                         tg_0_yyyzzz_0, tg_0_yyyzzz_1, tg_0_yyyzzzz_0, tg_0_yyyzzzz_1, tg_0_yyyzzzzz_0, \
                                         tg_0_yyzzzz_0, tg_0_yyzzzz_1, tg_0_yyzzzzz_0, tg_0_yyzzzzz_1, tg_0_yyzzzzzz_0, \
                                         tg_0_yzzzzz_0, tg_0_yzzzzz_1, tg_0_yzzzzzz_0, tg_0_yzzzzzz_1, tg_0_yzzzzzzz_0, \
                                         tg_0_zzzzzz_0, tg_0_zzzzzz_1, tg_0_zzzzzzz_0, tg_0_zzzzzzz_1, tg_0_zzzzzzzz_0, wq_x, \
                                         wq_y, wq_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fn = fn[j];

                    double fl1_fzb = fzb[j];
                    
                    double fr0 = qd_x[j]; double fr1 = wq_x[j];

                    tg_0_xxxxxxxx_0[j] = fr0 * tg_0_xxxxxxx_0[j] + fr1 * tg_0_xxxxxxx_1[j] + 3.5 * fl1_fn * (tg_0_xxxxxx_0[j] - fl1_fzb * tg_0_xxxxxx_1[j]);

                    tg_0_xxxxxxxy_0[j] = fr0 * tg_0_xxxxxxy_0[j] + fr1 * tg_0_xxxxxxy_1[j] + 3.0 * fl1_fn * (tg_0_xxxxxy_0[j] - fl1_fzb * tg_0_xxxxxy_1[j]);

                    tg_0_xxxxxxxz_0[j] = fr0 * tg_0_xxxxxxz_0[j] + fr1 * tg_0_xxxxxxz_1[j] + 3.0 * fl1_fn * (tg_0_xxxxxz_0[j] - fl1_fzb * tg_0_xxxxxz_1[j]);

                    tg_0_xxxxxxyy_0[j] = fr0 * tg_0_xxxxxyy_0[j] + fr1 * tg_0_xxxxxyy_1[j] + 2.5 * fl1_fn * (tg_0_xxxxyy_0[j] - fl1_fzb * tg_0_xxxxyy_1[j]);

                    tg_0_xxxxxxyz_0[j] = fr0 * tg_0_xxxxxyz_0[j] + fr1 * tg_0_xxxxxyz_1[j] + 2.5 * fl1_fn * (tg_0_xxxxyz_0[j] - fl1_fzb * tg_0_xxxxyz_1[j]);

                    tg_0_xxxxxxzz_0[j] = fr0 * tg_0_xxxxxzz_0[j] + fr1 * tg_0_xxxxxzz_1[j] + 2.5 * fl1_fn * (tg_0_xxxxzz_0[j] - fl1_fzb * tg_0_xxxxzz_1[j]);

                    tg_0_xxxxxyyy_0[j] = fr0 * tg_0_xxxxyyy_0[j] + fr1 * tg_0_xxxxyyy_1[j] + 2.0 * fl1_fn * (tg_0_xxxyyy_0[j] - fl1_fzb * tg_0_xxxyyy_1[j]);

                    tg_0_xxxxxyyz_0[j] = fr0 * tg_0_xxxxyyz_0[j] + fr1 * tg_0_xxxxyyz_1[j] + 2.0 * fl1_fn * (tg_0_xxxyyz_0[j] - fl1_fzb * tg_0_xxxyyz_1[j]);

                    tg_0_xxxxxyzz_0[j] = fr0 * tg_0_xxxxyzz_0[j] + fr1 * tg_0_xxxxyzz_1[j] + 2.0 * fl1_fn * (tg_0_xxxyzz_0[j] - fl1_fzb * tg_0_xxxyzz_1[j]);

                    tg_0_xxxxxzzz_0[j] = fr0 * tg_0_xxxxzzz_0[j] + fr1 * tg_0_xxxxzzz_1[j] + 2.0 * fl1_fn * (tg_0_xxxzzz_0[j] - fl1_fzb * tg_0_xxxzzz_1[j]);

                    tg_0_xxxxyyyy_0[j] = fr0 * tg_0_xxxyyyy_0[j] + fr1 * tg_0_xxxyyyy_1[j] + 1.5 * fl1_fn * (tg_0_xxyyyy_0[j] - fl1_fzb * tg_0_xxyyyy_1[j]);

                    tg_0_xxxxyyyz_0[j] = fr0 * tg_0_xxxyyyz_0[j] + fr1 * tg_0_xxxyyyz_1[j] + 1.5 * fl1_fn * (tg_0_xxyyyz_0[j] - fl1_fzb * tg_0_xxyyyz_1[j]);

                    tg_0_xxxxyyzz_0[j] = fr0 * tg_0_xxxyyzz_0[j] + fr1 * tg_0_xxxyyzz_1[j] + 1.5 * fl1_fn * (tg_0_xxyyzz_0[j] - fl1_fzb * tg_0_xxyyzz_1[j]);

                    tg_0_xxxxyzzz_0[j] = fr0 * tg_0_xxxyzzz_0[j] + fr1 * tg_0_xxxyzzz_1[j] + 1.5 * fl1_fn * (tg_0_xxyzzz_0[j] - fl1_fzb * tg_0_xxyzzz_1[j]);

                    tg_0_xxxxzzzz_0[j] = fr0 * tg_0_xxxzzzz_0[j] + fr1 * tg_0_xxxzzzz_1[j] + 1.5 * fl1_fn * (tg_0_xxzzzz_0[j] - fl1_fzb * tg_0_xxzzzz_1[j]);

                    tg_0_xxxyyyyy_0[j] = fr0 * tg_0_xxyyyyy_0[j] + fr1 * tg_0_xxyyyyy_1[j] + fl1_fn * (tg_0_xyyyyy_0[j] - fl1_fzb * tg_0_xyyyyy_1[j]);

                    tg_0_xxxyyyyz_0[j] = fr0 * tg_0_xxyyyyz_0[j] + fr1 * tg_0_xxyyyyz_1[j] + fl1_fn * (tg_0_xyyyyz_0[j] - fl1_fzb * tg_0_xyyyyz_1[j]);

                    tg_0_xxxyyyzz_0[j] = fr0 * tg_0_xxyyyzz_0[j] + fr1 * tg_0_xxyyyzz_1[j] + fl1_fn * (tg_0_xyyyzz_0[j] - fl1_fzb * tg_0_xyyyzz_1[j]);

                    tg_0_xxxyyzzz_0[j] = fr0 * tg_0_xxyyzzz_0[j] + fr1 * tg_0_xxyyzzz_1[j] + fl1_fn * (tg_0_xyyzzz_0[j] - fl1_fzb * tg_0_xyyzzz_1[j]);

                    tg_0_xxxyzzzz_0[j] = fr0 * tg_0_xxyzzzz_0[j] + fr1 * tg_0_xxyzzzz_1[j] + fl1_fn * (tg_0_xyzzzz_0[j] - fl1_fzb * tg_0_xyzzzz_1[j]);

                    tg_0_xxxzzzzz_0[j] = fr0 * tg_0_xxzzzzz_0[j] + fr1 * tg_0_xxzzzzz_1[j] + fl1_fn * (tg_0_xzzzzz_0[j] - fl1_fzb * tg_0_xzzzzz_1[j]);

                    tg_0_xxyyyyyy_0[j] = fr0 * tg_0_xyyyyyy_0[j] + fr1 * tg_0_xyyyyyy_1[j] + 0.5 * fl1_fn * (tg_0_yyyyyy_0[j] - fl1_fzb * tg_0_yyyyyy_1[j]);

                    tg_0_xxyyyyyz_0[j] = fr0 * tg_0_xyyyyyz_0[j] + fr1 * tg_0_xyyyyyz_1[j] + 0.5 * fl1_fn * (tg_0_yyyyyz_0[j] - fl1_fzb * tg_0_yyyyyz_1[j]);

                    tg_0_xxyyyyzz_0[j] = fr0 * tg_0_xyyyyzz_0[j] + fr1 * tg_0_xyyyyzz_1[j] + 0.5 * fl1_fn * (tg_0_yyyyzz_0[j] - fl1_fzb * tg_0_yyyyzz_1[j]);

                    tg_0_xxyyyzzz_0[j] = fr0 * tg_0_xyyyzzz_0[j] + fr1 * tg_0_xyyyzzz_1[j] + 0.5 * fl1_fn * (tg_0_yyyzzz_0[j] - fl1_fzb * tg_0_yyyzzz_1[j]);

                    tg_0_xxyyzzzz_0[j] = fr0 * tg_0_xyyzzzz_0[j] + fr1 * tg_0_xyyzzzz_1[j] + 0.5 * fl1_fn * (tg_0_yyzzzz_0[j] - fl1_fzb * tg_0_yyzzzz_1[j]);

                    tg_0_xxyzzzzz_0[j] = fr0 * tg_0_xyzzzzz_0[j] + fr1 * tg_0_xyzzzzz_1[j] + 0.5 * fl1_fn * (tg_0_yzzzzz_0[j] - fl1_fzb * tg_0_yzzzzz_1[j]);

                    tg_0_xxzzzzzz_0[j] = fr0 * tg_0_xzzzzzz_0[j] + fr1 * tg_0_xzzzzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzzzzz_0[j] - fl1_fzb * tg_0_zzzzzz_1[j]);

                    tg_0_xyyyyyyy_0[j] = fr0 * tg_0_yyyyyyy_0[j] + fr1 * tg_0_yyyyyyy_1[j];

                    tg_0_xyyyyyyz_0[j] = fr0 * tg_0_yyyyyyz_0[j] + fr1 * tg_0_yyyyyyz_1[j];

                    tg_0_xyyyyyzz_0[j] = fr0 * tg_0_yyyyyzz_0[j] + fr1 * tg_0_yyyyyzz_1[j];

                    tg_0_xyyyyzzz_0[j] = fr0 * tg_0_yyyyzzz_0[j] + fr1 * tg_0_yyyyzzz_1[j];

                    tg_0_xyyyzzzz_0[j] = fr0 * tg_0_yyyzzzz_0[j] + fr1 * tg_0_yyyzzzz_1[j];

                    tg_0_xyyzzzzz_0[j] = fr0 * tg_0_yyzzzzz_0[j] + fr1 * tg_0_yyzzzzz_1[j];

                    tg_0_xyzzzzzz_0[j] = fr0 * tg_0_yzzzzzz_0[j] + fr1 * tg_0_yzzzzzz_1[j];

                    tg_0_xzzzzzzz_0[j] = fr0 * tg_0_zzzzzzz_0[j] + fr1 * tg_0_zzzzzzz_1[j];
                                                                                                                 
                    fr0 = qd_y[j]; fr1 = wq_y[j];

                    tg_0_yyyyyyyy_0[j] = fr0 * tg_0_yyyyyyy_0[j] + fr1 * tg_0_yyyyyyy_1[j] + 3.5 * fl1_fn * (tg_0_yyyyyy_0[j] - fl1_fzb * tg_0_yyyyyy_1[j]);

                    tg_0_yyyyyyyz_0[j] = fr0 * tg_0_yyyyyyz_0[j] + fr1 * tg_0_yyyyyyz_1[j] + 3.0 * fl1_fn * (tg_0_yyyyyz_0[j] - fl1_fzb * tg_0_yyyyyz_1[j]);

                    tg_0_yyyyyyzz_0[j] = fr0 * tg_0_yyyyyzz_0[j] + fr1 * tg_0_yyyyyzz_1[j] + 2.5 * fl1_fn * (tg_0_yyyyzz_0[j] - fl1_fzb * tg_0_yyyyzz_1[j]);

                    tg_0_yyyyyzzz_0[j] = fr0 * tg_0_yyyyzzz_0[j] + fr1 * tg_0_yyyyzzz_1[j] + 2.0 * fl1_fn * (tg_0_yyyzzz_0[j] - fl1_fzb * tg_0_yyyzzz_1[j]);

                    tg_0_yyyyzzzz_0[j] = fr0 * tg_0_yyyzzzz_0[j] + fr1 * tg_0_yyyzzzz_1[j] + 1.5 * fl1_fn * (tg_0_yyzzzz_0[j] - fl1_fzb * tg_0_yyzzzz_1[j]);

                    tg_0_yyyzzzzz_0[j] = fr0 * tg_0_yyzzzzz_0[j] + fr1 * tg_0_yyzzzzz_1[j] + fl1_fn * (tg_0_yzzzzz_0[j] - fl1_fzb * tg_0_yzzzzz_1[j]);

                    tg_0_yyzzzzzz_0[j] = fr0 * tg_0_yzzzzzz_0[j] + fr1 * tg_0_yzzzzzz_1[j] + 0.5 * fl1_fn * (tg_0_zzzzzz_0[j] - fl1_fzb * tg_0_zzzzzz_1[j]);

                    tg_0_yzzzzzzz_0[j] = fr0 * tg_0_zzzzzzz_0[j] + fr1 * tg_0_zzzzzzz_1[j];

                    tg_0_zzzzzzzz_0[j] = qd_z[j] * tg_0_zzzzzzz_0[j] + wq_z[j] * tg_0_zzzzzzz_1[j] + 3.5 * fl1_fn * (tg_0_zzzzzz_0[j] - fl1_fzb * tg_0_zzzzzz_1[j]);
                }

                idx++;
            }
        }
    }

    void
    compElectronRepulsionForSLSS(      CMemBlock2D<double>* primBuffer,
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
                                             {8, -1, -1, -1},
                                             {0, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_g_8_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {8, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_g_8_0_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_g_7_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {7, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_7_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {7, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            auto pidx_g_6_0_m0 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {6, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord));

            auto pidx_g_6_0_m1 = recursionMap.index(CRecursionTerm({"Electron Repulsion"}, 0, true,
                                                                   {6, -1, -1, -1}, {0, -1, -1, -1},
                                                                   1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrPair]; i < epos[iContrPair]; i++)
            {
                // set up pointers to Obara-Saika factors

                double fx = b_fx[i];

                auto fza = osFactors.data(4 * idx + 2);

                // set up distances R(PB) = P - B

                auto pb_x = r_pb_x[i];

                auto pb_y = r_pb_y[i];

                auto pb_z = r_pb_z[i];

                // set up pointers to tensors product of distances R(WP) = W - P

                auto wp_x = wpDistances.data(3 * idx);

                auto wp_y = wpDistances.data(3 * idx + 1);

                auto wp_z = wpDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto tg_xxxxxxx_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx);

                auto tg_xxxxxxy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 1);

                auto tg_xxxxxxz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 2);

                auto tg_xxxxxyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 3);

                auto tg_xxxxxyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 4);

                auto tg_xxxxxzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 5);

                auto tg_xxxxyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 6);

                auto tg_xxxxyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 7);

                auto tg_xxxxyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 8);

                auto tg_xxxxzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 9);

                auto tg_xxxyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 10);

                auto tg_xxxyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 11);

                auto tg_xxxyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 12);

                auto tg_xxxyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 13);

                auto tg_xxxzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 14);

                auto tg_xxyyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 15);

                auto tg_xxyyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 16);

                auto tg_xxyyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 17);

                auto tg_xxyyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 18);

                auto tg_xxyzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 19);

                auto tg_xxzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 20);

                auto tg_xyyyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 21);

                auto tg_xyyyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 22);

                auto tg_xyyyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 23);

                auto tg_xyyyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 24);

                auto tg_xyyzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 25);

                auto tg_xyzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 26);

                auto tg_xzzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 27);

                auto tg_yyyyyyy_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 28);

                auto tg_yyyyyyz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 29);

                auto tg_yyyyyzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 30);

                auto tg_yyyyzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 31);

                auto tg_yyyzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 32);

                auto tg_yyzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 33);

                auto tg_yzzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 34);

                auto tg_zzzzzzz_0_0 = primBuffer[pidx_g_7_0_m0].data(36 * idx + 35);

                auto tg_xxxxxxx_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx);

                auto tg_xxxxxxy_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 1);

                auto tg_xxxxxxz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 2);

                auto tg_xxxxxyy_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 3);

                auto tg_xxxxxyz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 4);

                auto tg_xxxxxzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 5);

                auto tg_xxxxyyy_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 6);

                auto tg_xxxxyyz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 7);

                auto tg_xxxxyzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 8);

                auto tg_xxxxzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 9);

                auto tg_xxxyyyy_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 10);

                auto tg_xxxyyyz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 11);

                auto tg_xxxyyzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 12);

                auto tg_xxxyzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 13);

                auto tg_xxxzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 14);

                auto tg_xxyyyyy_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 15);

                auto tg_xxyyyyz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 16);

                auto tg_xxyyyzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 17);

                auto tg_xxyyzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 18);

                auto tg_xxyzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 19);

                auto tg_xxzzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 20);

                auto tg_xyyyyyy_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 21);

                auto tg_xyyyyyz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 22);

                auto tg_xyyyyzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 23);

                auto tg_xyyyzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 24);

                auto tg_xyyzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 25);

                auto tg_xyzzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 26);

                auto tg_xzzzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 27);

                auto tg_yyyyyyy_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 28);

                auto tg_yyyyyyz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 29);

                auto tg_yyyyyzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 30);

                auto tg_yyyyzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 31);

                auto tg_yyyzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 32);

                auto tg_yyzzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 33);

                auto tg_yzzzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 34);

                auto tg_zzzzzzz_0_1 = primBuffer[pidx_g_7_0_m1].data(36 * idx + 35);

                auto tg_xxxxxx_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx);

                auto tg_xxxxxy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 1);

                auto tg_xxxxxz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 2);

                auto tg_xxxxyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 3);

                auto tg_xxxxyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 4);

                auto tg_xxxxzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 5);

                auto tg_xxxyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 6);

                auto tg_xxxyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 7);

                auto tg_xxxyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 8);

                auto tg_xxxzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 9);

                auto tg_xxyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 10);

                auto tg_xxyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 11);

                auto tg_xxyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 12);

                auto tg_xxyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 13);

                auto tg_xxzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 14);

                auto tg_xyyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 15);

                auto tg_xyyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 16);

                auto tg_xyyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 17);

                auto tg_xyyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 18);

                auto tg_xyzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 19);

                auto tg_xzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 20);

                auto tg_yyyyyy_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 21);

                auto tg_yyyyyz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 22);

                auto tg_yyyyzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 23);

                auto tg_yyyzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 24);

                auto tg_yyzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 25);

                auto tg_yzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 26);

                auto tg_zzzzzz_0_0 = primBuffer[pidx_g_6_0_m0].data(28 * idx + 27);

                auto tg_xxxxxx_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx);

                auto tg_xxxxxy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 1);

                auto tg_xxxxxz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 2);

                auto tg_xxxxyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 3);

                auto tg_xxxxyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 4);

                auto tg_xxxxzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 5);

                auto tg_xxxyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 6);

                auto tg_xxxyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 7);

                auto tg_xxxyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 8);

                auto tg_xxxzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 9);

                auto tg_xxyyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 10);

                auto tg_xxyyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 11);

                auto tg_xxyyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 12);

                auto tg_xxyzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 13);

                auto tg_xxzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 14);

                auto tg_xyyyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 15);

                auto tg_xyyyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 16);

                auto tg_xyyyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 17);

                auto tg_xyyzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 18);

                auto tg_xyzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 19);

                auto tg_xzzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 20);

                auto tg_yyyyyy_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 21);

                auto tg_yyyyyz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 22);

                auto tg_yyyyzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 23);

                auto tg_yyyzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 24);

                auto tg_yyzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 25);

                auto tg_yzzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 26);

                auto tg_zzzzzz_0_1 = primBuffer[pidx_g_6_0_m1].data(28 * idx + 27);

                // set up pointers to integrals

                auto tg_xxxxxxxx_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx);

                auto tg_xxxxxxxy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 1);

                auto tg_xxxxxxxz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 2);

                auto tg_xxxxxxyy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 3);

                auto tg_xxxxxxyz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 4);

                auto tg_xxxxxxzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 5);

                auto tg_xxxxxyyy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 6);

                auto tg_xxxxxyyz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 7);

                auto tg_xxxxxyzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 8);

                auto tg_xxxxxzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 9);

                auto tg_xxxxyyyy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 10);

                auto tg_xxxxyyyz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 11);

                auto tg_xxxxyyzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 12);

                auto tg_xxxxyzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 13);

                auto tg_xxxxzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 14);

                auto tg_xxxyyyyy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 15);

                auto tg_xxxyyyyz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 16);

                auto tg_xxxyyyzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 17);

                auto tg_xxxyyzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 18);

                auto tg_xxxyzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 19);

                auto tg_xxxzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 20);

                auto tg_xxyyyyyy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 21);

                auto tg_xxyyyyyz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 22);

                auto tg_xxyyyyzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 23);

                auto tg_xxyyyzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 24);

                auto tg_xxyyzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 25);

                auto tg_xxyzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 26);

                auto tg_xxzzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 27);

                auto tg_xyyyyyyy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 28);

                auto tg_xyyyyyyz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 29);

                auto tg_xyyyyyzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 30);

                auto tg_xyyyyzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 31);

                auto tg_xyyyzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 32);

                auto tg_xyyzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 33);

                auto tg_xyzzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 34);

                auto tg_xzzzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 35);

                auto tg_yyyyyyyy_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 36);

                auto tg_yyyyyyyz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 37);

                auto tg_yyyyyyzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 38);

                auto tg_yyyyyzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 39);

                auto tg_yyyyzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 40);

                auto tg_yyyzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 41);

                auto tg_yyzzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 42);

                auto tg_yzzzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 43);

                auto tg_zzzzzzzz_0_0 = primBuffer[pidx_g_8_0_m0].data(45 * idx + 44); 

                #pragma omp simd aligned(fza, tg_xxxxxx_0_0, tg_xxxxxx_0_1, tg_xxxxxxx_0_0, tg_xxxxxxx_0_1, \
                                         tg_xxxxxxxx_0_0, tg_xxxxxxxy_0_0, tg_xxxxxxxz_0_0, tg_xxxxxxy_0_0, tg_xxxxxxy_0_1, \
                                         tg_xxxxxxyy_0_0, tg_xxxxxxyz_0_0, tg_xxxxxxz_0_0, tg_xxxxxxz_0_1, tg_xxxxxxzz_0_0, \
                                         tg_xxxxxy_0_0, tg_xxxxxy_0_1, tg_xxxxxyy_0_0, tg_xxxxxyy_0_1, tg_xxxxxyyy_0_0, \
                                         tg_xxxxxyyz_0_0, tg_xxxxxyz_0_0, tg_xxxxxyz_0_1, tg_xxxxxyzz_0_0, tg_xxxxxz_0_0, \
                                         tg_xxxxxz_0_1, tg_xxxxxzz_0_0, tg_xxxxxzz_0_1, tg_xxxxxzzz_0_0, tg_xxxxyy_0_0, \
                                         tg_xxxxyy_0_1, tg_xxxxyyy_0_0, tg_xxxxyyy_0_1, tg_xxxxyyyy_0_0, tg_xxxxyyyz_0_0, \
                                         tg_xxxxyyz_0_0, tg_xxxxyyz_0_1, tg_xxxxyyzz_0_0, tg_xxxxyz_0_0, tg_xxxxyz_0_1, \
                                         tg_xxxxyzz_0_0, tg_xxxxyzz_0_1, tg_xxxxyzzz_0_0, tg_xxxxzz_0_0, tg_xxxxzz_0_1, \
                                         tg_xxxxzzz_0_0, tg_xxxxzzz_0_1, tg_xxxxzzzz_0_0, tg_xxxyyy_0_0, tg_xxxyyy_0_1, \
                                         tg_xxxyyyy_0_0, tg_xxxyyyy_0_1, tg_xxxyyyyy_0_0, tg_xxxyyyyz_0_0, tg_xxxyyyz_0_0, \
                                         tg_xxxyyyz_0_1, tg_xxxyyyzz_0_0, tg_xxxyyz_0_0, tg_xxxyyz_0_1, tg_xxxyyzz_0_0, \
                                         tg_xxxyyzz_0_1, tg_xxxyyzzz_0_0, tg_xxxyzz_0_0, tg_xxxyzz_0_1, tg_xxxyzzz_0_0, \
                                         tg_xxxyzzz_0_1, tg_xxxyzzzz_0_0, tg_xxxzzz_0_0, tg_xxxzzz_0_1, tg_xxxzzzz_0_0, \
                                         tg_xxxzzzz_0_1, tg_xxxzzzzz_0_0, tg_xxyyyy_0_0, tg_xxyyyy_0_1, tg_xxyyyyy_0_0, \
                                         tg_xxyyyyy_0_1, tg_xxyyyyyy_0_0, tg_xxyyyyyz_0_0, tg_xxyyyyz_0_0, tg_xxyyyyz_0_1, \
                                         tg_xxyyyyzz_0_0, tg_xxyyyz_0_0, tg_xxyyyz_0_1, tg_xxyyyzz_0_0, tg_xxyyyzz_0_1, \
                                         tg_xxyyyzzz_0_0, tg_xxyyzz_0_0, tg_xxyyzz_0_1, tg_xxyyzzz_0_0, tg_xxyyzzz_0_1, \
                                         tg_xxyyzzzz_0_0, tg_xxyzzz_0_0, tg_xxyzzz_0_1, tg_xxyzzzz_0_0, tg_xxyzzzz_0_1, \
                                         tg_xxyzzzzz_0_0, tg_xxzzzz_0_0, tg_xxzzzz_0_1, tg_xxzzzzz_0_0, tg_xxzzzzz_0_1, \
                                         tg_xxzzzzzz_0_0, tg_xyyyyy_0_0, tg_xyyyyy_0_1, tg_xyyyyyy_0_0, tg_xyyyyyy_0_1, \
                                         tg_xyyyyyyy_0_0, tg_xyyyyyyz_0_0, tg_xyyyyyz_0_0, tg_xyyyyyz_0_1, tg_xyyyyyzz_0_0, \
                                         tg_xyyyyz_0_0, tg_xyyyyz_0_1, tg_xyyyyzz_0_0, tg_xyyyyzz_0_1, tg_xyyyyzzz_0_0, \
                                         tg_xyyyzz_0_0, tg_xyyyzz_0_1, tg_xyyyzzz_0_0, tg_xyyyzzz_0_1, tg_xyyyzzzz_0_0, \
                                         tg_xyyzzz_0_0, tg_xyyzzz_0_1, tg_xyyzzzz_0_0, tg_xyyzzzz_0_1, tg_xyyzzzzz_0_0, \
                                         tg_xyzzzz_0_0, tg_xyzzzz_0_1, tg_xyzzzzz_0_0, tg_xyzzzzz_0_1, tg_xyzzzzzz_0_0, \
                                         tg_xzzzzz_0_0, tg_xzzzzz_0_1, tg_xzzzzzz_0_0, tg_xzzzzzz_0_1, tg_xzzzzzzz_0_0, \
                                         tg_yyyyyy_0_0, tg_yyyyyy_0_1, tg_yyyyyyy_0_0, tg_yyyyyyy_0_1, tg_yyyyyyyy_0_0, \
                                         tg_yyyyyyyz_0_0, tg_yyyyyyz_0_0, tg_yyyyyyz_0_1, tg_yyyyyyzz_0_0, tg_yyyyyz_0_0, \
                                         tg_yyyyyz_0_1, tg_yyyyyzz_0_0, tg_yyyyyzz_0_1, tg_yyyyyzzz_0_0, tg_yyyyzz_0_0, \
                                         tg_yyyyzz_0_1, tg_yyyyzzz_0_0, tg_yyyyzzz_0_1, tg_yyyyzzzz_0_0, tg_yyyzzz_0_0, \
                                         tg_yyyzzz_0_1, tg_yyyzzzz_0_0, tg_yyyzzzz_0_1, tg_yyyzzzzz_0_0, tg_yyzzzz_0_0, \
                                         tg_yyzzzz_0_1, tg_yyzzzzz_0_0, tg_yyzzzzz_0_1, tg_yyzzzzzz_0_0, tg_yzzzzz_0_0, \
                                         tg_yzzzzz_0_1, tg_yzzzzzz_0_0, tg_yzzzzzz_0_1, tg_yzzzzzzz_0_0, tg_zzzzzz_0_0, \
                                         tg_zzzzzz_0_1, tg_zzzzzzz_0_0, tg_zzzzzzz_0_1, tg_zzzzzzzz_0_0, wp_x, wp_y, wp_z: VLX_ALIGN)
                for (int32_t j = 0; j < nKetPrimPairs; j++)
                {
                    double fl1_fx = fx;

                    double fl1_fza = fza[j];
                    
                    double fr1 = wp_x[j];

                    tg_xxxxxxxx_0_0[j] = pb_x * tg_xxxxxxx_0_0[j] + fr1 * tg_xxxxxxx_0_1[j] + 3.5 * fl1_fx * (tg_xxxxxx_0_0[j] -  fl1_fza * tg_xxxxxx_0_1[j]);

                    tg_xxxxxxxy_0_0[j] = pb_x * tg_xxxxxxy_0_0[j] + fr1 * tg_xxxxxxy_0_1[j] + 3.0 * fl1_fx * (tg_xxxxxy_0_0[j] -  fl1_fza * tg_xxxxxy_0_1[j]);

                    tg_xxxxxxxz_0_0[j] = pb_x * tg_xxxxxxz_0_0[j] + fr1 * tg_xxxxxxz_0_1[j] + 3.0 * fl1_fx * (tg_xxxxxz_0_0[j] -  fl1_fza * tg_xxxxxz_0_1[j]);

                    tg_xxxxxxyy_0_0[j] = pb_x * tg_xxxxxyy_0_0[j] + fr1 * tg_xxxxxyy_0_1[j] + 2.5 * fl1_fx * (tg_xxxxyy_0_0[j] -  fl1_fza * tg_xxxxyy_0_1[j]);

                    tg_xxxxxxyz_0_0[j] = pb_x * tg_xxxxxyz_0_0[j] + fr1 * tg_xxxxxyz_0_1[j] + 2.5 * fl1_fx * (tg_xxxxyz_0_0[j] -  fl1_fza * tg_xxxxyz_0_1[j]);

                    tg_xxxxxxzz_0_0[j] = pb_x * tg_xxxxxzz_0_0[j] + fr1 * tg_xxxxxzz_0_1[j] + 2.5 * fl1_fx * (tg_xxxxzz_0_0[j] -  fl1_fza * tg_xxxxzz_0_1[j]);

                    tg_xxxxxyyy_0_0[j] = pb_x * tg_xxxxyyy_0_0[j] + fr1 * tg_xxxxyyy_0_1[j] + 2.0 * fl1_fx * (tg_xxxyyy_0_0[j] -  fl1_fza * tg_xxxyyy_0_1[j]);

                    tg_xxxxxyyz_0_0[j] = pb_x * tg_xxxxyyz_0_0[j] + fr1 * tg_xxxxyyz_0_1[j] + 2.0 * fl1_fx * (tg_xxxyyz_0_0[j] -  fl1_fza * tg_xxxyyz_0_1[j]);

                    tg_xxxxxyzz_0_0[j] = pb_x * tg_xxxxyzz_0_0[j] + fr1 * tg_xxxxyzz_0_1[j] + 2.0 * fl1_fx * (tg_xxxyzz_0_0[j] -  fl1_fza * tg_xxxyzz_0_1[j]);

                    tg_xxxxxzzz_0_0[j] = pb_x * tg_xxxxzzz_0_0[j] + fr1 * tg_xxxxzzz_0_1[j] + 2.0 * fl1_fx * (tg_xxxzzz_0_0[j] -  fl1_fza * tg_xxxzzz_0_1[j]);

                    tg_xxxxyyyy_0_0[j] = pb_x * tg_xxxyyyy_0_0[j] + fr1 * tg_xxxyyyy_0_1[j] + 1.5 * fl1_fx * (tg_xxyyyy_0_0[j] -  fl1_fza * tg_xxyyyy_0_1[j]);

                    tg_xxxxyyyz_0_0[j] = pb_x * tg_xxxyyyz_0_0[j] + fr1 * tg_xxxyyyz_0_1[j] + 1.5 * fl1_fx * (tg_xxyyyz_0_0[j] -  fl1_fza * tg_xxyyyz_0_1[j]);

                    tg_xxxxyyzz_0_0[j] = pb_x * tg_xxxyyzz_0_0[j] + fr1 * tg_xxxyyzz_0_1[j] + 1.5 * fl1_fx * (tg_xxyyzz_0_0[j] -  fl1_fza * tg_xxyyzz_0_1[j]);

                    tg_xxxxyzzz_0_0[j] = pb_x * tg_xxxyzzz_0_0[j] + fr1 * tg_xxxyzzz_0_1[j] + 1.5 * fl1_fx * (tg_xxyzzz_0_0[j] -  fl1_fza * tg_xxyzzz_0_1[j]);

                    tg_xxxxzzzz_0_0[j] = pb_x * tg_xxxzzzz_0_0[j] + fr1 * tg_xxxzzzz_0_1[j] + 1.5 * fl1_fx * (tg_xxzzzz_0_0[j] -  fl1_fza * tg_xxzzzz_0_1[j]);

                    tg_xxxyyyyy_0_0[j] = pb_x * tg_xxyyyyy_0_0[j] + fr1 * tg_xxyyyyy_0_1[j] + fl1_fx * (tg_xyyyyy_0_0[j] -  fl1_fza * tg_xyyyyy_0_1[j]);

                    tg_xxxyyyyz_0_0[j] = pb_x * tg_xxyyyyz_0_0[j] + fr1 * tg_xxyyyyz_0_1[j] + fl1_fx * (tg_xyyyyz_0_0[j] -  fl1_fza * tg_xyyyyz_0_1[j]);

                    tg_xxxyyyzz_0_0[j] = pb_x * tg_xxyyyzz_0_0[j] + fr1 * tg_xxyyyzz_0_1[j] + fl1_fx * (tg_xyyyzz_0_0[j] -  fl1_fza * tg_xyyyzz_0_1[j]);

                    tg_xxxyyzzz_0_0[j] = pb_x * tg_xxyyzzz_0_0[j] + fr1 * tg_xxyyzzz_0_1[j] + fl1_fx * (tg_xyyzzz_0_0[j] -  fl1_fza * tg_xyyzzz_0_1[j]);

                    tg_xxxyzzzz_0_0[j] = pb_x * tg_xxyzzzz_0_0[j] + fr1 * tg_xxyzzzz_0_1[j] + fl1_fx * (tg_xyzzzz_0_0[j] -  fl1_fza * tg_xyzzzz_0_1[j]);

                    tg_xxxzzzzz_0_0[j] = pb_x * tg_xxzzzzz_0_0[j] + fr1 * tg_xxzzzzz_0_1[j] + fl1_fx * (tg_xzzzzz_0_0[j] -  fl1_fza * tg_xzzzzz_0_1[j]);

                    tg_xxyyyyyy_0_0[j] = pb_x * tg_xyyyyyy_0_0[j] + fr1 * tg_xyyyyyy_0_1[j] + 0.5 * fl1_fx * (tg_yyyyyy_0_0[j] -  fl1_fza * tg_yyyyyy_0_1[j]);

                    tg_xxyyyyyz_0_0[j] = pb_x * tg_xyyyyyz_0_0[j] + fr1 * tg_xyyyyyz_0_1[j] + 0.5 * fl1_fx * (tg_yyyyyz_0_0[j] -  fl1_fza * tg_yyyyyz_0_1[j]);

                    tg_xxyyyyzz_0_0[j] = pb_x * tg_xyyyyzz_0_0[j] + fr1 * tg_xyyyyzz_0_1[j] + 0.5 * fl1_fx * (tg_yyyyzz_0_0[j] -  fl1_fza * tg_yyyyzz_0_1[j]);

                    tg_xxyyyzzz_0_0[j] = pb_x * tg_xyyyzzz_0_0[j] + fr1 * tg_xyyyzzz_0_1[j] + 0.5 * fl1_fx * (tg_yyyzzz_0_0[j] -  fl1_fza * tg_yyyzzz_0_1[j]);

                    tg_xxyyzzzz_0_0[j] = pb_x * tg_xyyzzzz_0_0[j] + fr1 * tg_xyyzzzz_0_1[j] + 0.5 * fl1_fx * (tg_yyzzzz_0_0[j] -  fl1_fza * tg_yyzzzz_0_1[j]);

                    tg_xxyzzzzz_0_0[j] = pb_x * tg_xyzzzzz_0_0[j] + fr1 * tg_xyzzzzz_0_1[j] + 0.5 * fl1_fx * (tg_yzzzzz_0_0[j] -  fl1_fza * tg_yzzzzz_0_1[j]);

                    tg_xxzzzzzz_0_0[j] = pb_x * tg_xzzzzzz_0_0[j] + fr1 * tg_xzzzzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_0_0[j] -  fl1_fza * tg_zzzzzz_0_1[j]);

                    tg_xyyyyyyy_0_0[j] = pb_x * tg_yyyyyyy_0_0[j] + fr1 * tg_yyyyyyy_0_1[j];

                    tg_xyyyyyyz_0_0[j] = pb_x * tg_yyyyyyz_0_0[j] + fr1 * tg_yyyyyyz_0_1[j];

                    tg_xyyyyyzz_0_0[j] = pb_x * tg_yyyyyzz_0_0[j] + fr1 * tg_yyyyyzz_0_1[j];

                    tg_xyyyyzzz_0_0[j] = pb_x * tg_yyyyzzz_0_0[j] + fr1 * tg_yyyyzzz_0_1[j];

                    tg_xyyyzzzz_0_0[j] = pb_x * tg_yyyzzzz_0_0[j] + fr1 * tg_yyyzzzz_0_1[j];

                    tg_xyyzzzzz_0_0[j] = pb_x * tg_yyzzzzz_0_0[j] + fr1 * tg_yyzzzzz_0_1[j];

                    tg_xyzzzzzz_0_0[j] = pb_x * tg_yzzzzzz_0_0[j] + fr1 * tg_yzzzzzz_0_1[j];

                    tg_xzzzzzzz_0_0[j] = pb_x * tg_zzzzzzz_0_0[j] + fr1 * tg_zzzzzzz_0_1[j];
                                                                                                              
                    fr1 = wp_y[j];

                    tg_yyyyyyyy_0_0[j] = pb_y * tg_yyyyyyy_0_0[j] + fr1 * tg_yyyyyyy_0_1[j] + 3.5 * fl1_fx * (tg_yyyyyy_0_0[j] -  fl1_fza * tg_yyyyyy_0_1[j]);

                    tg_yyyyyyyz_0_0[j] = pb_y * tg_yyyyyyz_0_0[j] + fr1 * tg_yyyyyyz_0_1[j] + 3.0 * fl1_fx * (tg_yyyyyz_0_0[j] -  fl1_fza * tg_yyyyyz_0_1[j]);

                    tg_yyyyyyzz_0_0[j] = pb_y * tg_yyyyyzz_0_0[j] + fr1 * tg_yyyyyzz_0_1[j] + 2.5 * fl1_fx * (tg_yyyyzz_0_0[j] -  fl1_fza * tg_yyyyzz_0_1[j]);

                    tg_yyyyyzzz_0_0[j] = pb_y * tg_yyyyzzz_0_0[j] + fr1 * tg_yyyyzzz_0_1[j] + 2.0 * fl1_fx * (tg_yyyzzz_0_0[j] -  fl1_fza * tg_yyyzzz_0_1[j]);

                    tg_yyyyzzzz_0_0[j] = pb_y * tg_yyyzzzz_0_0[j] + fr1 * tg_yyyzzzz_0_1[j] + 1.5 * fl1_fx * (tg_yyzzzz_0_0[j] -  fl1_fza * tg_yyzzzz_0_1[j]);

                    tg_yyyzzzzz_0_0[j] = pb_y * tg_yyzzzzz_0_0[j] + fr1 * tg_yyzzzzz_0_1[j] + fl1_fx * (tg_yzzzzz_0_0[j] -  fl1_fza * tg_yzzzzz_0_1[j]);

                    tg_yyzzzzzz_0_0[j] = pb_y * tg_yzzzzzz_0_0[j] + fr1 * tg_yzzzzzz_0_1[j] + 0.5 * fl1_fx * (tg_zzzzzz_0_0[j] - fl1_fza * tg_zzzzzz_0_1[j]);

                    tg_yzzzzzzz_0_0[j] = pb_y * tg_zzzzzzz_0_0[j] + fr1 * tg_zzzzzzz_0_1[j];

                    tg_zzzzzzzz_0_0[j] = pb_z * tg_zzzzzzz_0_0[j] + wp_z[j] * tg_zzzzzzz_0_1[j] + 3.5 * fl1_fx * (tg_zzzzzz_0_0[j] - fl1_fza * tg_zzzzzz_0_1[j]);
                }

                idx++;
            }
        }
    }


} // erirecfunc namespace

