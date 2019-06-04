//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "NuclearPotentialRecFuncForGF.hpp"

namespace npotrecfunc { // npotrecfunc namespace

    void
    compNuclearPotentialForGF(      CMemBlock2D<double>& primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
    {
        npotrecfunc::compNuclearPotentialForGF_0_50(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    pcDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        npotrecfunc::compNuclearPotentialForGF_50_100(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        npotrecfunc::compNuclearPotentialForGF_100_150(primBuffer,
                                                       recursionMap,
                                                       osFactors,
                                                       paDistances, 
                                                       pcDistances, 
                                                       braGtoBlock,
                                                       ketGtoBlock,
                                                       iContrGto); 
    }

    void
    compNuclearPotentialForGF_0_50(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (0,50)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Nuclear Potential"},
                                             {4, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_a_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_a_4_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_a_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to tensors product of distances R(PA) = P - A

                auto pa_x = paDistances.data(3 * idx);

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_x = pcDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto ta_xxx_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx); 

                auto ta_xxx_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 1); 

                auto ta_xxx_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 2); 

                auto ta_xxx_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 3); 

                auto ta_xxx_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 4); 

                auto ta_xxx_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 5); 

                auto ta_xxx_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 6); 

                auto ta_xxx_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 7); 

                auto ta_xxx_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 8); 

                auto ta_xxx_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 9); 

                auto ta_xxy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 10); 

                auto ta_xxy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 11); 

                auto ta_xxy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 12); 

                auto ta_xxy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 13); 

                auto ta_xxy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 14); 

                auto ta_xxy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 15); 

                auto ta_xxy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 16); 

                auto ta_xxy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 17); 

                auto ta_xxy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 18); 

                auto ta_xxy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 19); 

                auto ta_xxz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 20); 

                auto ta_xxz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 21); 

                auto ta_xxz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 22); 

                auto ta_xxz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 23); 

                auto ta_xxz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 24); 

                auto ta_xxz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 25); 

                auto ta_xxz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 26); 

                auto ta_xxz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 27); 

                auto ta_xxz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 28); 

                auto ta_xxz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 29); 

                auto ta_xyy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 30); 

                auto ta_xyy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 31); 

                auto ta_xyy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 32); 

                auto ta_xyy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 33); 

                auto ta_xyy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 34); 

                auto ta_xyy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 35); 

                auto ta_xyy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 36); 

                auto ta_xyy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 37); 

                auto ta_xyy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 38); 

                auto ta_xyy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 39); 

                auto ta_xyz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 40); 

                auto ta_xyz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 41); 

                auto ta_xyz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 42); 

                auto ta_xyz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 43); 

                auto ta_xyz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 44); 

                auto ta_xyz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 45); 

                auto ta_xyz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 46); 

                auto ta_xyz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 47); 

                auto ta_xyz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 48); 

                auto ta_xyz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 49); 

                auto ta_xxx_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx); 

                auto ta_xxx_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 1); 

                auto ta_xxx_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 2); 

                auto ta_xxx_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 3); 

                auto ta_xxx_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 4); 

                auto ta_xxx_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 5); 

                auto ta_xxx_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 6); 

                auto ta_xxx_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 7); 

                auto ta_xxx_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 8); 

                auto ta_xxx_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 9); 

                auto ta_xxy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 10); 

                auto ta_xxy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 11); 

                auto ta_xxy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 12); 

                auto ta_xxy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 13); 

                auto ta_xxy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 14); 

                auto ta_xxy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 15); 

                auto ta_xxy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 16); 

                auto ta_xxy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 17); 

                auto ta_xxy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 18); 

                auto ta_xxy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 19); 

                auto ta_xxz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 20); 

                auto ta_xxz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 21); 

                auto ta_xxz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 22); 

                auto ta_xxz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 23); 

                auto ta_xxz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 24); 

                auto ta_xxz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 25); 

                auto ta_xxz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 26); 

                auto ta_xxz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 27); 

                auto ta_xxz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 28); 

                auto ta_xxz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 29); 

                auto ta_xyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 30); 

                auto ta_xyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 31); 

                auto ta_xyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 32); 

                auto ta_xyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 33); 

                auto ta_xyy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 34); 

                auto ta_xyy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 35); 

                auto ta_xyy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 36); 

                auto ta_xyy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 37); 

                auto ta_xyy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 38); 

                auto ta_xyy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 39); 

                auto ta_xyz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 40); 

                auto ta_xyz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 41); 

                auto ta_xyz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 42); 

                auto ta_xyz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 43); 

                auto ta_xyz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 44); 

                auto ta_xyz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 45); 

                auto ta_xyz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 46); 

                auto ta_xyz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 47); 

                auto ta_xyz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 48); 

                auto ta_xyz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 49); 

                auto ta_xx_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx); 

                auto ta_xx_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 1); 

                auto ta_xx_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 2); 

                auto ta_xx_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 3); 

                auto ta_xx_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 4); 

                auto ta_xx_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 5); 

                auto ta_xx_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 6); 

                auto ta_xx_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 7); 

                auto ta_xx_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 8); 

                auto ta_xx_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 9); 

                auto ta_xy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 10); 

                auto ta_xy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 11); 

                auto ta_xy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 12); 

                auto ta_xy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 13); 

                auto ta_xy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 14); 

                auto ta_xy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 15); 

                auto ta_xy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 16); 

                auto ta_xy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 17); 

                auto ta_xy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 18); 

                auto ta_xy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 19); 

                auto ta_xz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 20); 

                auto ta_xz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 21); 

                auto ta_xz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 22); 

                auto ta_xz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 23); 

                auto ta_xz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 24); 

                auto ta_xz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 25); 

                auto ta_xz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 26); 

                auto ta_xz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 27); 

                auto ta_xz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 28); 

                auto ta_xz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 29); 

                auto ta_yy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 30); 

                auto ta_yy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 31); 

                auto ta_yy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 32); 

                auto ta_yy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 33); 

                auto ta_yy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 34); 

                auto ta_yy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 35); 

                auto ta_yy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 36); 

                auto ta_yy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 37); 

                auto ta_yy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 38); 

                auto ta_yy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 39); 

                auto ta_yz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 40); 

                auto ta_yz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 41); 

                auto ta_yz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 42); 

                auto ta_yz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 43); 

                auto ta_yz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 44); 

                auto ta_yz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 45); 

                auto ta_yz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 46); 

                auto ta_yz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 47); 

                auto ta_yz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 48); 

                auto ta_yz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 49); 

                auto ta_xx_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx); 

                auto ta_xx_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 1); 

                auto ta_xx_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 2); 

                auto ta_xx_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 3); 

                auto ta_xx_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 4); 

                auto ta_xx_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 5); 

                auto ta_xx_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 6); 

                auto ta_xx_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 7); 

                auto ta_xx_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 8); 

                auto ta_xx_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 9); 

                auto ta_xy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 10); 

                auto ta_xy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 11); 

                auto ta_xy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 12); 

                auto ta_xy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 13); 

                auto ta_xy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 14); 

                auto ta_xy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 15); 

                auto ta_xy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 16); 

                auto ta_xy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 17); 

                auto ta_xy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 18); 

                auto ta_xy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 19); 

                auto ta_xz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 20); 

                auto ta_xz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 21); 

                auto ta_xz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 22); 

                auto ta_xz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 23); 

                auto ta_xz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 24); 

                auto ta_xz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 25); 

                auto ta_xz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 26); 

                auto ta_xz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 27); 

                auto ta_xz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 28); 

                auto ta_xz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 29); 

                auto ta_yy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 30); 

                auto ta_yy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 31); 

                auto ta_yy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 32); 

                auto ta_yy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 33); 

                auto ta_yy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 34); 

                auto ta_yy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 35); 

                auto ta_yy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 36); 

                auto ta_yy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 37); 

                auto ta_yy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 38); 

                auto ta_yy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 39); 

                auto ta_yz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 40); 

                auto ta_yz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 41); 

                auto ta_yz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 42); 

                auto ta_yz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 43); 

                auto ta_yz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 44); 

                auto ta_yz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 45); 

                auto ta_yz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 46); 

                auto ta_yz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 47); 

                auto ta_yz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 48); 

                auto ta_yz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 49); 

                auto ta_xxx_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx); 

                auto ta_xxx_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 1); 

                auto ta_xxx_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 2); 

                auto ta_xxx_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 3); 

                auto ta_xxx_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 4); 

                auto ta_xxx_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 5); 

                auto ta_xxy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 6); 

                auto ta_xxy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 7); 

                auto ta_xxy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 8); 

                auto ta_xxy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 9); 

                auto ta_xxy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 10); 

                auto ta_xxy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 11); 

                auto ta_xxz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 12); 

                auto ta_xxz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 13); 

                auto ta_xxz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 14); 

                auto ta_xxz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 15); 

                auto ta_xxz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 16); 

                auto ta_xxz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 17); 

                auto ta_xyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 18); 

                auto ta_xyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 19); 

                auto ta_xyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 20); 

                auto ta_xyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 21); 

                auto ta_xyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 22); 

                auto ta_xyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 23); 

                auto ta_xyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 24); 

                auto ta_xyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 25); 

                auto ta_xyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 26); 

                auto ta_xyz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 27); 

                auto ta_xyz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 28); 

                auto ta_xyz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 29); 

                auto ta_xxx_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx); 

                auto ta_xxx_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 1); 

                auto ta_xxx_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 2); 

                auto ta_xxx_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 3); 

                auto ta_xxx_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 4); 

                auto ta_xxx_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 5); 

                auto ta_xxy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 6); 

                auto ta_xxy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 7); 

                auto ta_xxy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 8); 

                auto ta_xxy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 9); 

                auto ta_xxy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 10); 

                auto ta_xxy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 11); 

                auto ta_xxz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 12); 

                auto ta_xxz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 13); 

                auto ta_xxz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 14); 

                auto ta_xxz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 15); 

                auto ta_xxz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 16); 

                auto ta_xxz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 17); 

                auto ta_xyy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 18); 

                auto ta_xyy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 19); 

                auto ta_xyy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 20); 

                auto ta_xyy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 21); 

                auto ta_xyy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 22); 

                auto ta_xyy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 23); 

                auto ta_xyz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 24); 

                auto ta_xyz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 25); 

                auto ta_xyz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 26); 

                auto ta_xyz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 27); 

                auto ta_xyz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 28); 

                auto ta_xyz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 29); 

                // set up pointers to integrals

                auto ta_xxxx_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx); 

                auto ta_xxxx_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 1); 

                auto ta_xxxx_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 2); 

                auto ta_xxxx_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 3); 

                auto ta_xxxx_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 4); 

                auto ta_xxxx_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 5); 

                auto ta_xxxx_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 6); 

                auto ta_xxxx_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 7); 

                auto ta_xxxx_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 8); 

                auto ta_xxxx_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 9); 

                auto ta_xxxy_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 10); 

                auto ta_xxxy_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 11); 

                auto ta_xxxy_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 12); 

                auto ta_xxxy_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 13); 

                auto ta_xxxy_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 14); 

                auto ta_xxxy_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 15); 

                auto ta_xxxy_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 16); 

                auto ta_xxxy_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 17); 

                auto ta_xxxy_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 18); 

                auto ta_xxxy_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 19); 

                auto ta_xxxz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 20); 

                auto ta_xxxz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 21); 

                auto ta_xxxz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 22); 

                auto ta_xxxz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 23); 

                auto ta_xxxz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 24); 

                auto ta_xxxz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 25); 

                auto ta_xxxz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 26); 

                auto ta_xxxz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 27); 

                auto ta_xxxz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 28); 

                auto ta_xxxz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 29); 

                auto ta_xxyy_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 30); 

                auto ta_xxyy_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 31); 

                auto ta_xxyy_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 32); 

                auto ta_xxyy_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 33); 

                auto ta_xxyy_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 34); 

                auto ta_xxyy_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 35); 

                auto ta_xxyy_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 36); 

                auto ta_xxyy_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 37); 

                auto ta_xxyy_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 38); 

                auto ta_xxyy_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 39); 

                auto ta_xxyz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 40); 

                auto ta_xxyz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 41); 

                auto ta_xxyz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 42); 

                auto ta_xxyz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 43); 

                auto ta_xxyz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 44); 

                auto ta_xxyz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 45); 

                auto ta_xxyz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 46); 

                auto ta_xxyz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 47); 

                auto ta_xxyz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 48); 

                auto ta_xxyz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 49); 

                // Batch of Integrals (0,50)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xx_xxx_0, ta_xx_xxx_1, ta_xx_xxy_0, ta_xx_xxy_1, \
                                         ta_xx_xxz_0, ta_xx_xxz_1, ta_xx_xyy_0, ta_xx_xyy_1, ta_xx_xyz_0, ta_xx_xyz_1, \
                                         ta_xx_xzz_0, ta_xx_xzz_1, ta_xx_yyy_0, ta_xx_yyy_1, ta_xx_yyz_0, ta_xx_yyz_1, \
                                         ta_xx_yzz_0, ta_xx_yzz_1, ta_xx_zzz_0, ta_xx_zzz_1, ta_xxx_xx_0, ta_xxx_xx_1, \
                                         ta_xxx_xxx_0, ta_xxx_xxx_1, ta_xxx_xxy_0, ta_xxx_xxy_1, ta_xxx_xxz_0, ta_xxx_xxz_1, \
                                         ta_xxx_xy_0, ta_xxx_xy_1, ta_xxx_xyy_0, ta_xxx_xyy_1, ta_xxx_xyz_0, ta_xxx_xyz_1, \
                                         ta_xxx_xz_0, ta_xxx_xz_1, ta_xxx_xzz_0, ta_xxx_xzz_1, ta_xxx_yy_0, ta_xxx_yy_1, \
                                         ta_xxx_yyy_0, ta_xxx_yyy_1, ta_xxx_yyz_0, ta_xxx_yyz_1, ta_xxx_yz_0, ta_xxx_yz_1, \
                                         ta_xxx_yzz_0, ta_xxx_yzz_1, ta_xxx_zz_0, ta_xxx_zz_1, ta_xxx_zzz_0, ta_xxx_zzz_1, \
                                         ta_xxxx_xxx_0, ta_xxxx_xxy_0, ta_xxxx_xxz_0, ta_xxxx_xyy_0, ta_xxxx_xyz_0, \
                                         ta_xxxx_xzz_0, ta_xxxx_yyy_0, ta_xxxx_yyz_0, ta_xxxx_yzz_0, ta_xxxx_zzz_0, \
                                         ta_xxxy_xxx_0, ta_xxxy_xxy_0, ta_xxxy_xxz_0, ta_xxxy_xyy_0, ta_xxxy_xyz_0, \
                                         ta_xxxy_xzz_0, ta_xxxy_yyy_0, ta_xxxy_yyz_0, ta_xxxy_yzz_0, ta_xxxy_zzz_0, \
                                         ta_xxxz_xxx_0, ta_xxxz_xxy_0, ta_xxxz_xxz_0, ta_xxxz_xyy_0, ta_xxxz_xyz_0, \
                                         ta_xxxz_xzz_0, ta_xxxz_yyy_0, ta_xxxz_yyz_0, ta_xxxz_yzz_0, ta_xxxz_zzz_0, \
                                         ta_xxy_xx_0, ta_xxy_xx_1, ta_xxy_xxx_0, ta_xxy_xxx_1, ta_xxy_xxy_0, ta_xxy_xxy_1, \
                                         ta_xxy_xxz_0, ta_xxy_xxz_1, ta_xxy_xy_0, ta_xxy_xy_1, ta_xxy_xyy_0, ta_xxy_xyy_1, \
                                         ta_xxy_xyz_0, ta_xxy_xyz_1, ta_xxy_xz_0, ta_xxy_xz_1, ta_xxy_xzz_0, ta_xxy_xzz_1, \
                                         ta_xxy_yy_0, ta_xxy_yy_1, ta_xxy_yyy_0, ta_xxy_yyy_1, ta_xxy_yyz_0, ta_xxy_yyz_1, \
                                         ta_xxy_yz_0, ta_xxy_yz_1, ta_xxy_yzz_0, ta_xxy_yzz_1, ta_xxy_zz_0, ta_xxy_zz_1, \
                                         ta_xxy_zzz_0, ta_xxy_zzz_1, ta_xxyy_xxx_0, ta_xxyy_xxy_0, ta_xxyy_xxz_0, \
                                         ta_xxyy_xyy_0, ta_xxyy_xyz_0, ta_xxyy_xzz_0, ta_xxyy_yyy_0, ta_xxyy_yyz_0, \
                                         ta_xxyy_yzz_0, ta_xxyy_zzz_0, ta_xxyz_xxx_0, ta_xxyz_xxy_0, ta_xxyz_xxz_0, \
                                         ta_xxyz_xyy_0, ta_xxyz_xyz_0, ta_xxyz_xzz_0, ta_xxyz_yyy_0, ta_xxyz_yyz_0, \
                                         ta_xxyz_yzz_0, ta_xxyz_zzz_0, ta_xxz_xx_0, ta_xxz_xx_1, ta_xxz_xxx_0, ta_xxz_xxx_1, \
                                         ta_xxz_xxy_0, ta_xxz_xxy_1, ta_xxz_xxz_0, ta_xxz_xxz_1, ta_xxz_xy_0, ta_xxz_xy_1, \
                                         ta_xxz_xyy_0, ta_xxz_xyy_1, ta_xxz_xyz_0, ta_xxz_xyz_1, ta_xxz_xz_0, ta_xxz_xz_1, \
                                         ta_xxz_xzz_0, ta_xxz_xzz_1, ta_xxz_yy_0, ta_xxz_yy_1, ta_xxz_yyy_0, ta_xxz_yyy_1, \
                                         ta_xxz_yyz_0, ta_xxz_yyz_1, ta_xxz_yz_0, ta_xxz_yz_1, ta_xxz_yzz_0, ta_xxz_yzz_1, \
                                         ta_xxz_zz_0, ta_xxz_zz_1, ta_xxz_zzz_0, ta_xxz_zzz_1, ta_xy_xxx_0, ta_xy_xxx_1, \
                                         ta_xy_xxy_0, ta_xy_xxy_1, ta_xy_xxz_0, ta_xy_xxz_1, ta_xy_xyy_0, ta_xy_xyy_1, \
                                         ta_xy_xyz_0, ta_xy_xyz_1, ta_xy_xzz_0, ta_xy_xzz_1, ta_xy_yyy_0, ta_xy_yyy_1, \
                                         ta_xy_yyz_0, ta_xy_yyz_1, ta_xy_yzz_0, ta_xy_yzz_1, ta_xy_zzz_0, ta_xy_zzz_1, \
                                         ta_xyy_xx_0, ta_xyy_xx_1, ta_xyy_xxx_0, ta_xyy_xxx_1, ta_xyy_xxy_0, ta_xyy_xxy_1, \
                                         ta_xyy_xxz_0, ta_xyy_xxz_1, ta_xyy_xy_0, ta_xyy_xy_1, ta_xyy_xyy_0, ta_xyy_xyy_1, \
                                         ta_xyy_xyz_0, ta_xyy_xyz_1, ta_xyy_xz_0, ta_xyy_xz_1, ta_xyy_xzz_0, ta_xyy_xzz_1, \
                                         ta_xyy_yy_0, ta_xyy_yy_1, ta_xyy_yyy_0, ta_xyy_yyy_1, ta_xyy_yyz_0, ta_xyy_yyz_1, \
                                         ta_xyy_yz_0, ta_xyy_yz_1, ta_xyy_yzz_0, ta_xyy_yzz_1, ta_xyy_zz_0, ta_xyy_zz_1, \
                                         ta_xyy_zzz_0, ta_xyy_zzz_1, ta_xyz_xx_0, ta_xyz_xx_1, ta_xyz_xxx_0, ta_xyz_xxx_1, \
                                         ta_xyz_xxy_0, ta_xyz_xxy_1, ta_xyz_xxz_0, ta_xyz_xxz_1, ta_xyz_xy_0, ta_xyz_xy_1, \
                                         ta_xyz_xyy_0, ta_xyz_xyy_1, ta_xyz_xyz_0, ta_xyz_xyz_1, ta_xyz_xz_0, ta_xyz_xz_1, \
                                         ta_xyz_xzz_0, ta_xyz_xzz_1, ta_xyz_yy_0, ta_xyz_yy_1, ta_xyz_yyy_0, ta_xyz_yyy_1, \
                                         ta_xyz_yyz_0, ta_xyz_yyz_1, ta_xyz_yz_0, ta_xyz_yz_1, ta_xyz_yzz_0, ta_xyz_yzz_1, \
                                         ta_xyz_zz_0, ta_xyz_zz_1, ta_xyz_zzz_0, ta_xyz_zzz_1, ta_xz_xxx_0, ta_xz_xxx_1, \
                                         ta_xz_xxy_0, ta_xz_xxy_1, ta_xz_xxz_0, ta_xz_xxz_1, ta_xz_xyy_0, ta_xz_xyy_1, \
                                         ta_xz_xyz_0, ta_xz_xyz_1, ta_xz_xzz_0, ta_xz_xzz_1, ta_xz_yyy_0, ta_xz_yyy_1, \
                                         ta_xz_yyz_0, ta_xz_yyz_1, ta_xz_yzz_0, ta_xz_yzz_1, ta_xz_zzz_0, ta_xz_zzz_1, \
                                         ta_yy_xxx_0, ta_yy_xxx_1, ta_yy_xxy_0, ta_yy_xxy_1, ta_yy_xxz_0, ta_yy_xxz_1, \
                                         ta_yy_xyy_0, ta_yy_xyy_1, ta_yy_xyz_0, ta_yy_xyz_1, ta_yy_xzz_0, ta_yy_xzz_1, \
                                         ta_yy_yyy_0, ta_yy_yyy_1, ta_yy_yyz_0, ta_yy_yyz_1, ta_yy_yzz_0, ta_yy_yzz_1, \
                                         ta_yy_zzz_0, ta_yy_zzz_1, ta_yz_xxx_0, ta_yz_xxx_1, ta_yz_xxy_0, ta_yz_xxy_1, \
                                         ta_yz_xxz_0, ta_yz_xxz_1, ta_yz_xyy_0, ta_yz_xyy_1, ta_yz_xyz_0, ta_yz_xyz_1, \
                                         ta_yz_xzz_0, ta_yz_xzz_1, ta_yz_yyy_0, ta_yz_yyy_1, ta_yz_yyz_0, ta_yz_yyz_1, \
                                         ta_yz_yzz_0, ta_yz_yzz_1, ta_yz_zzz_0, ta_yz_zzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    ta_xxxx_xxx_0[j] = pa_x[j] * ta_xxx_xxx_0[j] - pc_x[j] * ta_xxx_xxx_1[j] + 1.5 * fl1_fx * ta_xx_xxx_0[j] - 1.5 * fl1_fx * ta_xx_xxx_1[j] + 1.5 * fl1_fx * ta_xxx_xx_0[j] - 1.5 * fl1_fx * ta_xxx_xx_1[j];

                    ta_xxxx_xxy_0[j] = pa_x[j] * ta_xxx_xxy_0[j] - pc_x[j] * ta_xxx_xxy_1[j] + 1.5 * fl1_fx * ta_xx_xxy_0[j] - 1.5 * fl1_fx * ta_xx_xxy_1[j] + fl1_fx * ta_xxx_xy_0[j] - fl1_fx * ta_xxx_xy_1[j];

                    ta_xxxx_xxz_0[j] = pa_x[j] * ta_xxx_xxz_0[j] - pc_x[j] * ta_xxx_xxz_1[j] + 1.5 * fl1_fx * ta_xx_xxz_0[j] - 1.5 * fl1_fx * ta_xx_xxz_1[j] + fl1_fx * ta_xxx_xz_0[j] - fl1_fx * ta_xxx_xz_1[j];

                    ta_xxxx_xyy_0[j] = pa_x[j] * ta_xxx_xyy_0[j] - pc_x[j] * ta_xxx_xyy_1[j] + 1.5 * fl1_fx * ta_xx_xyy_0[j] - 1.5 * fl1_fx * ta_xx_xyy_1[j] + 0.5 * fl1_fx * ta_xxx_yy_0[j] - 0.5 * fl1_fx * ta_xxx_yy_1[j];

                    ta_xxxx_xyz_0[j] = pa_x[j] * ta_xxx_xyz_0[j] - pc_x[j] * ta_xxx_xyz_1[j] + 1.5 * fl1_fx * ta_xx_xyz_0[j] - 1.5 * fl1_fx * ta_xx_xyz_1[j] + 0.5 * fl1_fx * ta_xxx_yz_0[j] - 0.5 * fl1_fx * ta_xxx_yz_1[j];

                    ta_xxxx_xzz_0[j] = pa_x[j] * ta_xxx_xzz_0[j] - pc_x[j] * ta_xxx_xzz_1[j] + 1.5 * fl1_fx * ta_xx_xzz_0[j] - 1.5 * fl1_fx * ta_xx_xzz_1[j] + 0.5 * fl1_fx * ta_xxx_zz_0[j] - 0.5 * fl1_fx * ta_xxx_zz_1[j];

                    ta_xxxx_yyy_0[j] = pa_x[j] * ta_xxx_yyy_0[j] - pc_x[j] * ta_xxx_yyy_1[j] + 1.5 * fl1_fx * ta_xx_yyy_0[j] - 1.5 * fl1_fx * ta_xx_yyy_1[j];

                    ta_xxxx_yyz_0[j] = pa_x[j] * ta_xxx_yyz_0[j] - pc_x[j] * ta_xxx_yyz_1[j] + 1.5 * fl1_fx * ta_xx_yyz_0[j] - 1.5 * fl1_fx * ta_xx_yyz_1[j];

                    ta_xxxx_yzz_0[j] = pa_x[j] * ta_xxx_yzz_0[j] - pc_x[j] * ta_xxx_yzz_1[j] + 1.5 * fl1_fx * ta_xx_yzz_0[j] - 1.5 * fl1_fx * ta_xx_yzz_1[j];

                    ta_xxxx_zzz_0[j] = pa_x[j] * ta_xxx_zzz_0[j] - pc_x[j] * ta_xxx_zzz_1[j] + 1.5 * fl1_fx * ta_xx_zzz_0[j] - 1.5 * fl1_fx * ta_xx_zzz_1[j];

                    ta_xxxy_xxx_0[j] = pa_x[j] * ta_xxy_xxx_0[j] - pc_x[j] * ta_xxy_xxx_1[j] + fl1_fx * ta_xy_xxx_0[j] - fl1_fx * ta_xy_xxx_1[j] + 1.5 * fl1_fx * ta_xxy_xx_0[j] - 1.5 * fl1_fx * ta_xxy_xx_1[j];

                    ta_xxxy_xxy_0[j] = pa_x[j] * ta_xxy_xxy_0[j] - pc_x[j] * ta_xxy_xxy_1[j] + fl1_fx * ta_xy_xxy_0[j] - fl1_fx * ta_xy_xxy_1[j] + fl1_fx * ta_xxy_xy_0[j] - fl1_fx * ta_xxy_xy_1[j];

                    ta_xxxy_xxz_0[j] = pa_x[j] * ta_xxy_xxz_0[j] - pc_x[j] * ta_xxy_xxz_1[j] + fl1_fx * ta_xy_xxz_0[j] - fl1_fx * ta_xy_xxz_1[j] + fl1_fx * ta_xxy_xz_0[j] - fl1_fx * ta_xxy_xz_1[j];

                    ta_xxxy_xyy_0[j] = pa_x[j] * ta_xxy_xyy_0[j] - pc_x[j] * ta_xxy_xyy_1[j] + fl1_fx * ta_xy_xyy_0[j] - fl1_fx * ta_xy_xyy_1[j] + 0.5 * fl1_fx * ta_xxy_yy_0[j] - 0.5 * fl1_fx * ta_xxy_yy_1[j];

                    ta_xxxy_xyz_0[j] = pa_x[j] * ta_xxy_xyz_0[j] - pc_x[j] * ta_xxy_xyz_1[j] + fl1_fx * ta_xy_xyz_0[j] - fl1_fx * ta_xy_xyz_1[j] + 0.5 * fl1_fx * ta_xxy_yz_0[j] - 0.5 * fl1_fx * ta_xxy_yz_1[j];

                    ta_xxxy_xzz_0[j] = pa_x[j] * ta_xxy_xzz_0[j] - pc_x[j] * ta_xxy_xzz_1[j] + fl1_fx * ta_xy_xzz_0[j] - fl1_fx * ta_xy_xzz_1[j] + 0.5 * fl1_fx * ta_xxy_zz_0[j] - 0.5 * fl1_fx * ta_xxy_zz_1[j];

                    ta_xxxy_yyy_0[j] = pa_x[j] * ta_xxy_yyy_0[j] - pc_x[j] * ta_xxy_yyy_1[j] + fl1_fx * ta_xy_yyy_0[j] - fl1_fx * ta_xy_yyy_1[j];

                    ta_xxxy_yyz_0[j] = pa_x[j] * ta_xxy_yyz_0[j] - pc_x[j] * ta_xxy_yyz_1[j] + fl1_fx * ta_xy_yyz_0[j] - fl1_fx * ta_xy_yyz_1[j];

                    ta_xxxy_yzz_0[j] = pa_x[j] * ta_xxy_yzz_0[j] - pc_x[j] * ta_xxy_yzz_1[j] + fl1_fx * ta_xy_yzz_0[j] - fl1_fx * ta_xy_yzz_1[j];

                    ta_xxxy_zzz_0[j] = pa_x[j] * ta_xxy_zzz_0[j] - pc_x[j] * ta_xxy_zzz_1[j] + fl1_fx * ta_xy_zzz_0[j] - fl1_fx * ta_xy_zzz_1[j];

                    ta_xxxz_xxx_0[j] = pa_x[j] * ta_xxz_xxx_0[j] - pc_x[j] * ta_xxz_xxx_1[j] + fl1_fx * ta_xz_xxx_0[j] - fl1_fx * ta_xz_xxx_1[j] + 1.5 * fl1_fx * ta_xxz_xx_0[j] - 1.5 * fl1_fx * ta_xxz_xx_1[j];

                    ta_xxxz_xxy_0[j] = pa_x[j] * ta_xxz_xxy_0[j] - pc_x[j] * ta_xxz_xxy_1[j] + fl1_fx * ta_xz_xxy_0[j] - fl1_fx * ta_xz_xxy_1[j] + fl1_fx * ta_xxz_xy_0[j] - fl1_fx * ta_xxz_xy_1[j];

                    ta_xxxz_xxz_0[j] = pa_x[j] * ta_xxz_xxz_0[j] - pc_x[j] * ta_xxz_xxz_1[j] + fl1_fx * ta_xz_xxz_0[j] - fl1_fx * ta_xz_xxz_1[j] + fl1_fx * ta_xxz_xz_0[j] - fl1_fx * ta_xxz_xz_1[j];

                    ta_xxxz_xyy_0[j] = pa_x[j] * ta_xxz_xyy_0[j] - pc_x[j] * ta_xxz_xyy_1[j] + fl1_fx * ta_xz_xyy_0[j] - fl1_fx * ta_xz_xyy_1[j] + 0.5 * fl1_fx * ta_xxz_yy_0[j] - 0.5 * fl1_fx * ta_xxz_yy_1[j];

                    ta_xxxz_xyz_0[j] = pa_x[j] * ta_xxz_xyz_0[j] - pc_x[j] * ta_xxz_xyz_1[j] + fl1_fx * ta_xz_xyz_0[j] - fl1_fx * ta_xz_xyz_1[j] + 0.5 * fl1_fx * ta_xxz_yz_0[j] - 0.5 * fl1_fx * ta_xxz_yz_1[j];

                    ta_xxxz_xzz_0[j] = pa_x[j] * ta_xxz_xzz_0[j] - pc_x[j] * ta_xxz_xzz_1[j] + fl1_fx * ta_xz_xzz_0[j] - fl1_fx * ta_xz_xzz_1[j] + 0.5 * fl1_fx * ta_xxz_zz_0[j] - 0.5 * fl1_fx * ta_xxz_zz_1[j];

                    ta_xxxz_yyy_0[j] = pa_x[j] * ta_xxz_yyy_0[j] - pc_x[j] * ta_xxz_yyy_1[j] + fl1_fx * ta_xz_yyy_0[j] - fl1_fx * ta_xz_yyy_1[j];

                    ta_xxxz_yyz_0[j] = pa_x[j] * ta_xxz_yyz_0[j] - pc_x[j] * ta_xxz_yyz_1[j] + fl1_fx * ta_xz_yyz_0[j] - fl1_fx * ta_xz_yyz_1[j];

                    ta_xxxz_yzz_0[j] = pa_x[j] * ta_xxz_yzz_0[j] - pc_x[j] * ta_xxz_yzz_1[j] + fl1_fx * ta_xz_yzz_0[j] - fl1_fx * ta_xz_yzz_1[j];

                    ta_xxxz_zzz_0[j] = pa_x[j] * ta_xxz_zzz_0[j] - pc_x[j] * ta_xxz_zzz_1[j] + fl1_fx * ta_xz_zzz_0[j] - fl1_fx * ta_xz_zzz_1[j];

                    ta_xxyy_xxx_0[j] = pa_x[j] * ta_xyy_xxx_0[j] - pc_x[j] * ta_xyy_xxx_1[j] + 0.5 * fl1_fx * ta_yy_xxx_0[j] - 0.5 * fl1_fx * ta_yy_xxx_1[j] + 1.5 * fl1_fx * ta_xyy_xx_0[j] - 1.5 * fl1_fx * ta_xyy_xx_1[j];

                    ta_xxyy_xxy_0[j] = pa_x[j] * ta_xyy_xxy_0[j] - pc_x[j] * ta_xyy_xxy_1[j] + 0.5 * fl1_fx * ta_yy_xxy_0[j] - 0.5 * fl1_fx * ta_yy_xxy_1[j] + fl1_fx * ta_xyy_xy_0[j] - fl1_fx * ta_xyy_xy_1[j];

                    ta_xxyy_xxz_0[j] = pa_x[j] * ta_xyy_xxz_0[j] - pc_x[j] * ta_xyy_xxz_1[j] + 0.5 * fl1_fx * ta_yy_xxz_0[j] - 0.5 * fl1_fx * ta_yy_xxz_1[j] + fl1_fx * ta_xyy_xz_0[j] - fl1_fx * ta_xyy_xz_1[j];

                    ta_xxyy_xyy_0[j] = pa_x[j] * ta_xyy_xyy_0[j] - pc_x[j] * ta_xyy_xyy_1[j] + 0.5 * fl1_fx * ta_yy_xyy_0[j] - 0.5 * fl1_fx * ta_yy_xyy_1[j] + 0.5 * fl1_fx * ta_xyy_yy_0[j] - 0.5 * fl1_fx * ta_xyy_yy_1[j];

                    ta_xxyy_xyz_0[j] = pa_x[j] * ta_xyy_xyz_0[j] - pc_x[j] * ta_xyy_xyz_1[j] + 0.5 * fl1_fx * ta_yy_xyz_0[j] - 0.5 * fl1_fx * ta_yy_xyz_1[j] + 0.5 * fl1_fx * ta_xyy_yz_0[j] - 0.5 * fl1_fx * ta_xyy_yz_1[j];

                    ta_xxyy_xzz_0[j] = pa_x[j] * ta_xyy_xzz_0[j] - pc_x[j] * ta_xyy_xzz_1[j] + 0.5 * fl1_fx * ta_yy_xzz_0[j] - 0.5 * fl1_fx * ta_yy_xzz_1[j] + 0.5 * fl1_fx * ta_xyy_zz_0[j] - 0.5 * fl1_fx * ta_xyy_zz_1[j];

                    ta_xxyy_yyy_0[j] = pa_x[j] * ta_xyy_yyy_0[j] - pc_x[j] * ta_xyy_yyy_1[j] + 0.5 * fl1_fx * ta_yy_yyy_0[j] - 0.5 * fl1_fx * ta_yy_yyy_1[j];

                    ta_xxyy_yyz_0[j] = pa_x[j] * ta_xyy_yyz_0[j] - pc_x[j] * ta_xyy_yyz_1[j] + 0.5 * fl1_fx * ta_yy_yyz_0[j] - 0.5 * fl1_fx * ta_yy_yyz_1[j];

                    ta_xxyy_yzz_0[j] = pa_x[j] * ta_xyy_yzz_0[j] - pc_x[j] * ta_xyy_yzz_1[j] + 0.5 * fl1_fx * ta_yy_yzz_0[j] - 0.5 * fl1_fx * ta_yy_yzz_1[j];

                    ta_xxyy_zzz_0[j] = pa_x[j] * ta_xyy_zzz_0[j] - pc_x[j] * ta_xyy_zzz_1[j] + 0.5 * fl1_fx * ta_yy_zzz_0[j] - 0.5 * fl1_fx * ta_yy_zzz_1[j];

                    ta_xxyz_xxx_0[j] = pa_x[j] * ta_xyz_xxx_0[j] - pc_x[j] * ta_xyz_xxx_1[j] + 0.5 * fl1_fx * ta_yz_xxx_0[j] - 0.5 * fl1_fx * ta_yz_xxx_1[j] + 1.5 * fl1_fx * ta_xyz_xx_0[j] - 1.5 * fl1_fx * ta_xyz_xx_1[j];

                    ta_xxyz_xxy_0[j] = pa_x[j] * ta_xyz_xxy_0[j] - pc_x[j] * ta_xyz_xxy_1[j] + 0.5 * fl1_fx * ta_yz_xxy_0[j] - 0.5 * fl1_fx * ta_yz_xxy_1[j] + fl1_fx * ta_xyz_xy_0[j] - fl1_fx * ta_xyz_xy_1[j];

                    ta_xxyz_xxz_0[j] = pa_x[j] * ta_xyz_xxz_0[j] - pc_x[j] * ta_xyz_xxz_1[j] + 0.5 * fl1_fx * ta_yz_xxz_0[j] - 0.5 * fl1_fx * ta_yz_xxz_1[j] + fl1_fx * ta_xyz_xz_0[j] - fl1_fx * ta_xyz_xz_1[j];

                    ta_xxyz_xyy_0[j] = pa_x[j] * ta_xyz_xyy_0[j] - pc_x[j] * ta_xyz_xyy_1[j] + 0.5 * fl1_fx * ta_yz_xyy_0[j] - 0.5 * fl1_fx * ta_yz_xyy_1[j] + 0.5 * fl1_fx * ta_xyz_yy_0[j] - 0.5 * fl1_fx * ta_xyz_yy_1[j];

                    ta_xxyz_xyz_0[j] = pa_x[j] * ta_xyz_xyz_0[j] - pc_x[j] * ta_xyz_xyz_1[j] + 0.5 * fl1_fx * ta_yz_xyz_0[j] - 0.5 * fl1_fx * ta_yz_xyz_1[j] + 0.5 * fl1_fx * ta_xyz_yz_0[j] - 0.5 * fl1_fx * ta_xyz_yz_1[j];

                    ta_xxyz_xzz_0[j] = pa_x[j] * ta_xyz_xzz_0[j] - pc_x[j] * ta_xyz_xzz_1[j] + 0.5 * fl1_fx * ta_yz_xzz_0[j] - 0.5 * fl1_fx * ta_yz_xzz_1[j] + 0.5 * fl1_fx * ta_xyz_zz_0[j] - 0.5 * fl1_fx * ta_xyz_zz_1[j];

                    ta_xxyz_yyy_0[j] = pa_x[j] * ta_xyz_yyy_0[j] - pc_x[j] * ta_xyz_yyy_1[j] + 0.5 * fl1_fx * ta_yz_yyy_0[j] - 0.5 * fl1_fx * ta_yz_yyy_1[j];

                    ta_xxyz_yyz_0[j] = pa_x[j] * ta_xyz_yyz_0[j] - pc_x[j] * ta_xyz_yyz_1[j] + 0.5 * fl1_fx * ta_yz_yyz_0[j] - 0.5 * fl1_fx * ta_yz_yyz_1[j];

                    ta_xxyz_yzz_0[j] = pa_x[j] * ta_xyz_yzz_0[j] - pc_x[j] * ta_xyz_yzz_1[j] + 0.5 * fl1_fx * ta_yz_yzz_0[j] - 0.5 * fl1_fx * ta_yz_yzz_1[j];

                    ta_xxyz_zzz_0[j] = pa_x[j] * ta_xyz_zzz_0[j] - pc_x[j] * ta_xyz_zzz_1[j] + 0.5 * fl1_fx * ta_yz_zzz_0[j] - 0.5 * fl1_fx * ta_yz_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compNuclearPotentialForGF_50_100(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
                                     const CMemBlock2D<double>& pcDistances,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto)
    {
        // Batch of Integrals (50,100)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Nuclear Potential"},
                                             {4, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_a_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_a_4_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_a_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to tensors product of distances R(PA) = P - A

                auto pa_x = paDistances.data(3 * idx);

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_x = pcDistances.data(3 * idx);

                // set up pointers to auxilary integrals

                auto ta_xzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 50); 

                auto ta_xzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 51); 

                auto ta_xzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 52); 

                auto ta_xzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 53); 

                auto ta_xzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 54); 

                auto ta_xzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 55); 

                auto ta_xzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 56); 

                auto ta_xzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 57); 

                auto ta_xzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 58); 

                auto ta_xzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 59); 

                auto ta_yyy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 60); 

                auto ta_yyy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 61); 

                auto ta_yyy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 62); 

                auto ta_yyy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 63); 

                auto ta_yyy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 64); 

                auto ta_yyy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 65); 

                auto ta_yyy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 66); 

                auto ta_yyy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 67); 

                auto ta_yyy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 68); 

                auto ta_yyy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 69); 

                auto ta_yyz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 70); 

                auto ta_yyz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 71); 

                auto ta_yyz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 72); 

                auto ta_yyz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 73); 

                auto ta_yyz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 74); 

                auto ta_yyz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 75); 

                auto ta_yyz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 76); 

                auto ta_yyz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 77); 

                auto ta_yyz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 78); 

                auto ta_yyz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 79); 

                auto ta_yzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 80); 

                auto ta_yzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 81); 

                auto ta_yzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 82); 

                auto ta_yzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 83); 

                auto ta_yzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 84); 

                auto ta_yzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 85); 

                auto ta_yzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 86); 

                auto ta_yzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 87); 

                auto ta_yzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 88); 

                auto ta_yzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 89); 

                auto ta_zzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 90); 

                auto ta_zzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 91); 

                auto ta_zzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 92); 

                auto ta_zzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 93); 

                auto ta_zzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 94); 

                auto ta_zzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 95); 

                auto ta_zzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 96); 

                auto ta_zzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 97); 

                auto ta_zzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 98); 

                auto ta_zzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 99); 

                auto ta_xzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 50); 

                auto ta_xzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 51); 

                auto ta_xzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 52); 

                auto ta_xzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 53); 

                auto ta_xzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 54); 

                auto ta_xzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 55); 

                auto ta_xzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 56); 

                auto ta_xzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 57); 

                auto ta_xzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 58); 

                auto ta_xzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 59); 

                auto ta_yyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 60); 

                auto ta_yyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 61); 

                auto ta_yyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 62); 

                auto ta_yyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 63); 

                auto ta_yyy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 64); 

                auto ta_yyy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 65); 

                auto ta_yyy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 66); 

                auto ta_yyy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 67); 

                auto ta_yyy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 68); 

                auto ta_yyy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 69); 

                auto ta_yyz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 70); 

                auto ta_yyz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 71); 

                auto ta_yyz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 72); 

                auto ta_yyz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 73); 

                auto ta_yyz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 74); 

                auto ta_yyz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 75); 

                auto ta_yyz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 76); 

                auto ta_yyz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 77); 

                auto ta_yyz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 78); 

                auto ta_yyz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 79); 

                auto ta_yzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 80); 

                auto ta_yzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 81); 

                auto ta_yzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 82); 

                auto ta_yzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 83); 

                auto ta_yzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 84); 

                auto ta_yzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 85); 

                auto ta_yzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 86); 

                auto ta_yzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 87); 

                auto ta_yzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 88); 

                auto ta_yzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 89); 

                auto ta_zzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 90); 

                auto ta_zzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 91); 

                auto ta_zzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 92); 

                auto ta_zzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 93); 

                auto ta_zzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 94); 

                auto ta_zzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 95); 

                auto ta_zzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 96); 

                auto ta_zzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 97); 

                auto ta_zzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 98); 

                auto ta_zzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 99); 

                auto ta_zz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 50); 

                auto ta_zz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 51); 

                auto ta_zz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 52); 

                auto ta_zz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 53); 

                auto ta_zz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 54); 

                auto ta_zz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 55); 

                auto ta_zz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 56); 

                auto ta_zz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 57); 

                auto ta_zz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 58); 

                auto ta_zz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 59); 

                auto ta_zz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 50); 

                auto ta_zz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 51); 

                auto ta_zz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 52); 

                auto ta_zz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 53); 

                auto ta_zz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 54); 

                auto ta_zz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 55); 

                auto ta_zz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 56); 

                auto ta_zz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 57); 

                auto ta_zz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 58); 

                auto ta_zz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 59); 

                auto ta_xzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 30); 

                auto ta_xzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 31); 

                auto ta_xzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 32); 

                auto ta_xzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 33); 

                auto ta_xzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 34); 

                auto ta_xzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 35); 

                auto ta_yyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 36); 

                auto ta_yyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 37); 

                auto ta_yyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 38); 

                auto ta_yyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 39); 

                auto ta_yyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 40); 

                auto ta_yyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 41); 

                auto ta_yyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 42); 

                auto ta_yyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 43); 

                auto ta_yyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 44); 

                auto ta_yyz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 45); 

                auto ta_yyz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 46); 

                auto ta_yyz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 47); 

                auto ta_yzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 48); 

                auto ta_yzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 49); 

                auto ta_yzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 50); 

                auto ta_yzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 51); 

                auto ta_yzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 52); 

                auto ta_yzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 53); 

                auto ta_zzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 54); 

                auto ta_zzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 55); 

                auto ta_zzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 56); 

                auto ta_zzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 57); 

                auto ta_zzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 58); 

                auto ta_zzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 59); 

                auto ta_xzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 30); 

                auto ta_xzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 31); 

                auto ta_xzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 32); 

                auto ta_xzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 33); 

                auto ta_xzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 34); 

                auto ta_xzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 35); 

                auto ta_yyy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 36); 

                auto ta_yyy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 37); 

                auto ta_yyy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 38); 

                auto ta_yyy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 39); 

                auto ta_yyy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 40); 

                auto ta_yyy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 41); 

                auto ta_yyz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 42); 

                auto ta_yyz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 43); 

                auto ta_yyz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 44); 

                auto ta_yyz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 45); 

                auto ta_yyz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 46); 

                auto ta_yyz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 47); 

                auto ta_yzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 48); 

                auto ta_yzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 49); 

                auto ta_yzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 50); 

                auto ta_yzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 51); 

                auto ta_yzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 52); 

                auto ta_yzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 53); 

                auto ta_zzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 54); 

                auto ta_zzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 55); 

                auto ta_zzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 56); 

                auto ta_zzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 57); 

                auto ta_zzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 58); 

                auto ta_zzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 59); 

                // set up pointers to integrals

                auto ta_xxzz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 50); 

                auto ta_xxzz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 51); 

                auto ta_xxzz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 52); 

                auto ta_xxzz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 53); 

                auto ta_xxzz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 54); 

                auto ta_xxzz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 55); 

                auto ta_xxzz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 56); 

                auto ta_xxzz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 57); 

                auto ta_xxzz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 58); 

                auto ta_xxzz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 59); 

                auto ta_xyyy_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 60); 

                auto ta_xyyy_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 61); 

                auto ta_xyyy_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 62); 

                auto ta_xyyy_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 63); 

                auto ta_xyyy_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 64); 

                auto ta_xyyy_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 65); 

                auto ta_xyyy_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 66); 

                auto ta_xyyy_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 67); 

                auto ta_xyyy_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 68); 

                auto ta_xyyy_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 69); 

                auto ta_xyyz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 70); 

                auto ta_xyyz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 71); 

                auto ta_xyyz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 72); 

                auto ta_xyyz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 73); 

                auto ta_xyyz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 74); 

                auto ta_xyyz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 75); 

                auto ta_xyyz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 76); 

                auto ta_xyyz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 77); 

                auto ta_xyyz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 78); 

                auto ta_xyyz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 79); 

                auto ta_xyzz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 80); 

                auto ta_xyzz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 81); 

                auto ta_xyzz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 82); 

                auto ta_xyzz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 83); 

                auto ta_xyzz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 84); 

                auto ta_xyzz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 85); 

                auto ta_xyzz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 86); 

                auto ta_xyzz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 87); 

                auto ta_xyzz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 88); 

                auto ta_xyzz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 89); 

                auto ta_xzzz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 90); 

                auto ta_xzzz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 91); 

                auto ta_xzzz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 92); 

                auto ta_xzzz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 93); 

                auto ta_xzzz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 94); 

                auto ta_xzzz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 95); 

                auto ta_xzzz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 96); 

                auto ta_xzzz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 97); 

                auto ta_xzzz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 98); 

                auto ta_xzzz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 99); 

                // Batch of Integrals (50,100)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxzz_xxx_0, ta_xxzz_xxy_0, ta_xxzz_xxz_0, \
                                         ta_xxzz_xyy_0, ta_xxzz_xyz_0, ta_xxzz_xzz_0, ta_xxzz_yyy_0, ta_xxzz_yyz_0, \
                                         ta_xxzz_yzz_0, ta_xxzz_zzz_0, ta_xyyy_xxx_0, ta_xyyy_xxy_0, ta_xyyy_xxz_0, \
                                         ta_xyyy_xyy_0, ta_xyyy_xyz_0, ta_xyyy_xzz_0, ta_xyyy_yyy_0, ta_xyyy_yyz_0, \
                                         ta_xyyy_yzz_0, ta_xyyy_zzz_0, ta_xyyz_xxx_0, ta_xyyz_xxy_0, ta_xyyz_xxz_0, \
                                         ta_xyyz_xyy_0, ta_xyyz_xyz_0, ta_xyyz_xzz_0, ta_xyyz_yyy_0, ta_xyyz_yyz_0, \
                                         ta_xyyz_yzz_0, ta_xyyz_zzz_0, ta_xyzz_xxx_0, ta_xyzz_xxy_0, ta_xyzz_xxz_0, \
                                         ta_xyzz_xyy_0, ta_xyzz_xyz_0, ta_xyzz_xzz_0, ta_xyzz_yyy_0, ta_xyzz_yyz_0, \
                                         ta_xyzz_yzz_0, ta_xyzz_zzz_0, ta_xzz_xx_0, ta_xzz_xx_1, ta_xzz_xxx_0, ta_xzz_xxx_1, \
                                         ta_xzz_xxy_0, ta_xzz_xxy_1, ta_xzz_xxz_0, ta_xzz_xxz_1, ta_xzz_xy_0, ta_xzz_xy_1, \
                                         ta_xzz_xyy_0, ta_xzz_xyy_1, ta_xzz_xyz_0, ta_xzz_xyz_1, ta_xzz_xz_0, ta_xzz_xz_1, \
                                         ta_xzz_xzz_0, ta_xzz_xzz_1, ta_xzz_yy_0, ta_xzz_yy_1, ta_xzz_yyy_0, ta_xzz_yyy_1, \
                                         ta_xzz_yyz_0, ta_xzz_yyz_1, ta_xzz_yz_0, ta_xzz_yz_1, ta_xzz_yzz_0, ta_xzz_yzz_1, \
                                         ta_xzz_zz_0, ta_xzz_zz_1, ta_xzz_zzz_0, ta_xzz_zzz_1, ta_xzzz_xxx_0, \
                                         ta_xzzz_xxy_0, ta_xzzz_xxz_0, ta_xzzz_xyy_0, ta_xzzz_xyz_0, ta_xzzz_xzz_0, \
                                         ta_xzzz_yyy_0, ta_xzzz_yyz_0, ta_xzzz_yzz_0, ta_xzzz_zzz_0, ta_yyy_xx_0, \
                                         ta_yyy_xx_1, ta_yyy_xxx_0, ta_yyy_xxx_1, ta_yyy_xxy_0, ta_yyy_xxy_1, ta_yyy_xxz_0, \
                                         ta_yyy_xxz_1, ta_yyy_xy_0, ta_yyy_xy_1, ta_yyy_xyy_0, ta_yyy_xyy_1, ta_yyy_xyz_0, \
                                         ta_yyy_xyz_1, ta_yyy_xz_0, ta_yyy_xz_1, ta_yyy_xzz_0, ta_yyy_xzz_1, ta_yyy_yy_0, \
                                         ta_yyy_yy_1, ta_yyy_yyy_0, ta_yyy_yyy_1, ta_yyy_yyz_0, ta_yyy_yyz_1, ta_yyy_yz_0, \
                                         ta_yyy_yz_1, ta_yyy_yzz_0, ta_yyy_yzz_1, ta_yyy_zz_0, ta_yyy_zz_1, ta_yyy_zzz_0, \
                                         ta_yyy_zzz_1, ta_yyz_xx_0, ta_yyz_xx_1, ta_yyz_xxx_0, ta_yyz_xxx_1, ta_yyz_xxy_0, \
                                         ta_yyz_xxy_1, ta_yyz_xxz_0, ta_yyz_xxz_1, ta_yyz_xy_0, ta_yyz_xy_1, ta_yyz_xyy_0, \
                                         ta_yyz_xyy_1, ta_yyz_xyz_0, ta_yyz_xyz_1, ta_yyz_xz_0, ta_yyz_xz_1, ta_yyz_xzz_0, \
                                         ta_yyz_xzz_1, ta_yyz_yy_0, ta_yyz_yy_1, ta_yyz_yyy_0, ta_yyz_yyy_1, ta_yyz_yyz_0, \
                                         ta_yyz_yyz_1, ta_yyz_yz_0, ta_yyz_yz_1, ta_yyz_yzz_0, ta_yyz_yzz_1, ta_yyz_zz_0, \
                                         ta_yyz_zz_1, ta_yyz_zzz_0, ta_yyz_zzz_1, ta_yzz_xx_0, ta_yzz_xx_1, ta_yzz_xxx_0, \
                                         ta_yzz_xxx_1, ta_yzz_xxy_0, ta_yzz_xxy_1, ta_yzz_xxz_0, ta_yzz_xxz_1, ta_yzz_xy_0, \
                                         ta_yzz_xy_1, ta_yzz_xyy_0, ta_yzz_xyy_1, ta_yzz_xyz_0, ta_yzz_xyz_1, ta_yzz_xz_0, \
                                         ta_yzz_xz_1, ta_yzz_xzz_0, ta_yzz_xzz_1, ta_yzz_yy_0, ta_yzz_yy_1, ta_yzz_yyy_0, \
                                         ta_yzz_yyy_1, ta_yzz_yyz_0, ta_yzz_yyz_1, ta_yzz_yz_0, ta_yzz_yz_1, ta_yzz_yzz_0, \
                                         ta_yzz_yzz_1, ta_yzz_zz_0, ta_yzz_zz_1, ta_yzz_zzz_0, ta_yzz_zzz_1, ta_zz_xxx_0, \
                                         ta_zz_xxx_1, ta_zz_xxy_0, ta_zz_xxy_1, ta_zz_xxz_0, ta_zz_xxz_1, ta_zz_xyy_0, \
                                         ta_zz_xyy_1, ta_zz_xyz_0, ta_zz_xyz_1, ta_zz_xzz_0, ta_zz_xzz_1, ta_zz_yyy_0, \
                                         ta_zz_yyy_1, ta_zz_yyz_0, ta_zz_yyz_1, ta_zz_yzz_0, ta_zz_yzz_1, ta_zz_zzz_0, \
                                         ta_zz_zzz_1, ta_zzz_xx_0, ta_zzz_xx_1, ta_zzz_xxx_0, ta_zzz_xxx_1, ta_zzz_xxy_0, \
                                         ta_zzz_xxy_1, ta_zzz_xxz_0, ta_zzz_xxz_1, ta_zzz_xy_0, ta_zzz_xy_1, ta_zzz_xyy_0, \
                                         ta_zzz_xyy_1, ta_zzz_xyz_0, ta_zzz_xyz_1, ta_zzz_xz_0, ta_zzz_xz_1, ta_zzz_xzz_0, \
                                         ta_zzz_xzz_1, ta_zzz_yy_0, ta_zzz_yy_1, ta_zzz_yyy_0, ta_zzz_yyy_1, ta_zzz_yyz_0, \
                                         ta_zzz_yyz_1, ta_zzz_yz_0, ta_zzz_yz_1, ta_zzz_yzz_0, ta_zzz_yzz_1, ta_zzz_zz_0, \
                                         ta_zzz_zz_1, ta_zzz_zzz_0, ta_zzz_zzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    ta_xxzz_xxx_0[j] = pa_x[j] * ta_xzz_xxx_0[j] - pc_x[j] * ta_xzz_xxx_1[j] + 0.5 * fl1_fx * ta_zz_xxx_0[j] - 0.5 * fl1_fx * ta_zz_xxx_1[j] + 1.5 * fl1_fx * ta_xzz_xx_0[j] - 1.5 * fl1_fx * ta_xzz_xx_1[j];

                    ta_xxzz_xxy_0[j] = pa_x[j] * ta_xzz_xxy_0[j] - pc_x[j] * ta_xzz_xxy_1[j] + 0.5 * fl1_fx * ta_zz_xxy_0[j] - 0.5 * fl1_fx * ta_zz_xxy_1[j] + fl1_fx * ta_xzz_xy_0[j] - fl1_fx * ta_xzz_xy_1[j];

                    ta_xxzz_xxz_0[j] = pa_x[j] * ta_xzz_xxz_0[j] - pc_x[j] * ta_xzz_xxz_1[j] + 0.5 * fl1_fx * ta_zz_xxz_0[j] - 0.5 * fl1_fx * ta_zz_xxz_1[j] + fl1_fx * ta_xzz_xz_0[j] - fl1_fx * ta_xzz_xz_1[j];

                    ta_xxzz_xyy_0[j] = pa_x[j] * ta_xzz_xyy_0[j] - pc_x[j] * ta_xzz_xyy_1[j] + 0.5 * fl1_fx * ta_zz_xyy_0[j] - 0.5 * fl1_fx * ta_zz_xyy_1[j] + 0.5 * fl1_fx * ta_xzz_yy_0[j] - 0.5 * fl1_fx * ta_xzz_yy_1[j];

                    ta_xxzz_xyz_0[j] = pa_x[j] * ta_xzz_xyz_0[j] - pc_x[j] * ta_xzz_xyz_1[j] + 0.5 * fl1_fx * ta_zz_xyz_0[j] - 0.5 * fl1_fx * ta_zz_xyz_1[j] + 0.5 * fl1_fx * ta_xzz_yz_0[j] - 0.5 * fl1_fx * ta_xzz_yz_1[j];

                    ta_xxzz_xzz_0[j] = pa_x[j] * ta_xzz_xzz_0[j] - pc_x[j] * ta_xzz_xzz_1[j] + 0.5 * fl1_fx * ta_zz_xzz_0[j] - 0.5 * fl1_fx * ta_zz_xzz_1[j] + 0.5 * fl1_fx * ta_xzz_zz_0[j] - 0.5 * fl1_fx * ta_xzz_zz_1[j];

                    ta_xxzz_yyy_0[j] = pa_x[j] * ta_xzz_yyy_0[j] - pc_x[j] * ta_xzz_yyy_1[j] + 0.5 * fl1_fx * ta_zz_yyy_0[j] - 0.5 * fl1_fx * ta_zz_yyy_1[j];

                    ta_xxzz_yyz_0[j] = pa_x[j] * ta_xzz_yyz_0[j] - pc_x[j] * ta_xzz_yyz_1[j] + 0.5 * fl1_fx * ta_zz_yyz_0[j] - 0.5 * fl1_fx * ta_zz_yyz_1[j];

                    ta_xxzz_yzz_0[j] = pa_x[j] * ta_xzz_yzz_0[j] - pc_x[j] * ta_xzz_yzz_1[j] + 0.5 * fl1_fx * ta_zz_yzz_0[j] - 0.5 * fl1_fx * ta_zz_yzz_1[j];

                    ta_xxzz_zzz_0[j] = pa_x[j] * ta_xzz_zzz_0[j] - pc_x[j] * ta_xzz_zzz_1[j] + 0.5 * fl1_fx * ta_zz_zzz_0[j] - 0.5 * fl1_fx * ta_zz_zzz_1[j];

                    ta_xyyy_xxx_0[j] = pa_x[j] * ta_yyy_xxx_0[j] - pc_x[j] * ta_yyy_xxx_1[j] + 1.5 * fl1_fx * ta_yyy_xx_0[j] - 1.5 * fl1_fx * ta_yyy_xx_1[j];

                    ta_xyyy_xxy_0[j] = pa_x[j] * ta_yyy_xxy_0[j] - pc_x[j] * ta_yyy_xxy_1[j] + fl1_fx * ta_yyy_xy_0[j] - fl1_fx * ta_yyy_xy_1[j];

                    ta_xyyy_xxz_0[j] = pa_x[j] * ta_yyy_xxz_0[j] - pc_x[j] * ta_yyy_xxz_1[j] + fl1_fx * ta_yyy_xz_0[j] - fl1_fx * ta_yyy_xz_1[j];

                    ta_xyyy_xyy_0[j] = pa_x[j] * ta_yyy_xyy_0[j] - pc_x[j] * ta_yyy_xyy_1[j] + 0.5 * fl1_fx * ta_yyy_yy_0[j] - 0.5 * fl1_fx * ta_yyy_yy_1[j];

                    ta_xyyy_xyz_0[j] = pa_x[j] * ta_yyy_xyz_0[j] - pc_x[j] * ta_yyy_xyz_1[j] + 0.5 * fl1_fx * ta_yyy_yz_0[j] - 0.5 * fl1_fx * ta_yyy_yz_1[j];

                    ta_xyyy_xzz_0[j] = pa_x[j] * ta_yyy_xzz_0[j] - pc_x[j] * ta_yyy_xzz_1[j] + 0.5 * fl1_fx * ta_yyy_zz_0[j] - 0.5 * fl1_fx * ta_yyy_zz_1[j];

                    ta_xyyy_yyy_0[j] = pa_x[j] * ta_yyy_yyy_0[j] - pc_x[j] * ta_yyy_yyy_1[j];

                    ta_xyyy_yyz_0[j] = pa_x[j] * ta_yyy_yyz_0[j] - pc_x[j] * ta_yyy_yyz_1[j];

                    ta_xyyy_yzz_0[j] = pa_x[j] * ta_yyy_yzz_0[j] - pc_x[j] * ta_yyy_yzz_1[j];

                    ta_xyyy_zzz_0[j] = pa_x[j] * ta_yyy_zzz_0[j] - pc_x[j] * ta_yyy_zzz_1[j];

                    ta_xyyz_xxx_0[j] = pa_x[j] * ta_yyz_xxx_0[j] - pc_x[j] * ta_yyz_xxx_1[j] + 1.5 * fl1_fx * ta_yyz_xx_0[j] - 1.5 * fl1_fx * ta_yyz_xx_1[j];

                    ta_xyyz_xxy_0[j] = pa_x[j] * ta_yyz_xxy_0[j] - pc_x[j] * ta_yyz_xxy_1[j] + fl1_fx * ta_yyz_xy_0[j] - fl1_fx * ta_yyz_xy_1[j];

                    ta_xyyz_xxz_0[j] = pa_x[j] * ta_yyz_xxz_0[j] - pc_x[j] * ta_yyz_xxz_1[j] + fl1_fx * ta_yyz_xz_0[j] - fl1_fx * ta_yyz_xz_1[j];

                    ta_xyyz_xyy_0[j] = pa_x[j] * ta_yyz_xyy_0[j] - pc_x[j] * ta_yyz_xyy_1[j] + 0.5 * fl1_fx * ta_yyz_yy_0[j] - 0.5 * fl1_fx * ta_yyz_yy_1[j];

                    ta_xyyz_xyz_0[j] = pa_x[j] * ta_yyz_xyz_0[j] - pc_x[j] * ta_yyz_xyz_1[j] + 0.5 * fl1_fx * ta_yyz_yz_0[j] - 0.5 * fl1_fx * ta_yyz_yz_1[j];

                    ta_xyyz_xzz_0[j] = pa_x[j] * ta_yyz_xzz_0[j] - pc_x[j] * ta_yyz_xzz_1[j] + 0.5 * fl1_fx * ta_yyz_zz_0[j] - 0.5 * fl1_fx * ta_yyz_zz_1[j];

                    ta_xyyz_yyy_0[j] = pa_x[j] * ta_yyz_yyy_0[j] - pc_x[j] * ta_yyz_yyy_1[j];

                    ta_xyyz_yyz_0[j] = pa_x[j] * ta_yyz_yyz_0[j] - pc_x[j] * ta_yyz_yyz_1[j];

                    ta_xyyz_yzz_0[j] = pa_x[j] * ta_yyz_yzz_0[j] - pc_x[j] * ta_yyz_yzz_1[j];

                    ta_xyyz_zzz_0[j] = pa_x[j] * ta_yyz_zzz_0[j] - pc_x[j] * ta_yyz_zzz_1[j];

                    ta_xyzz_xxx_0[j] = pa_x[j] * ta_yzz_xxx_0[j] - pc_x[j] * ta_yzz_xxx_1[j] + 1.5 * fl1_fx * ta_yzz_xx_0[j] - 1.5 * fl1_fx * ta_yzz_xx_1[j];

                    ta_xyzz_xxy_0[j] = pa_x[j] * ta_yzz_xxy_0[j] - pc_x[j] * ta_yzz_xxy_1[j] + fl1_fx * ta_yzz_xy_0[j] - fl1_fx * ta_yzz_xy_1[j];

                    ta_xyzz_xxz_0[j] = pa_x[j] * ta_yzz_xxz_0[j] - pc_x[j] * ta_yzz_xxz_1[j] + fl1_fx * ta_yzz_xz_0[j] - fl1_fx * ta_yzz_xz_1[j];

                    ta_xyzz_xyy_0[j] = pa_x[j] * ta_yzz_xyy_0[j] - pc_x[j] * ta_yzz_xyy_1[j] + 0.5 * fl1_fx * ta_yzz_yy_0[j] - 0.5 * fl1_fx * ta_yzz_yy_1[j];

                    ta_xyzz_xyz_0[j] = pa_x[j] * ta_yzz_xyz_0[j] - pc_x[j] * ta_yzz_xyz_1[j] + 0.5 * fl1_fx * ta_yzz_yz_0[j] - 0.5 * fl1_fx * ta_yzz_yz_1[j];

                    ta_xyzz_xzz_0[j] = pa_x[j] * ta_yzz_xzz_0[j] - pc_x[j] * ta_yzz_xzz_1[j] + 0.5 * fl1_fx * ta_yzz_zz_0[j] - 0.5 * fl1_fx * ta_yzz_zz_1[j];

                    ta_xyzz_yyy_0[j] = pa_x[j] * ta_yzz_yyy_0[j] - pc_x[j] * ta_yzz_yyy_1[j];

                    ta_xyzz_yyz_0[j] = pa_x[j] * ta_yzz_yyz_0[j] - pc_x[j] * ta_yzz_yyz_1[j];

                    ta_xyzz_yzz_0[j] = pa_x[j] * ta_yzz_yzz_0[j] - pc_x[j] * ta_yzz_yzz_1[j];

                    ta_xyzz_zzz_0[j] = pa_x[j] * ta_yzz_zzz_0[j] - pc_x[j] * ta_yzz_zzz_1[j];

                    ta_xzzz_xxx_0[j] = pa_x[j] * ta_zzz_xxx_0[j] - pc_x[j] * ta_zzz_xxx_1[j] + 1.5 * fl1_fx * ta_zzz_xx_0[j] - 1.5 * fl1_fx * ta_zzz_xx_1[j];

                    ta_xzzz_xxy_0[j] = pa_x[j] * ta_zzz_xxy_0[j] - pc_x[j] * ta_zzz_xxy_1[j] + fl1_fx * ta_zzz_xy_0[j] - fl1_fx * ta_zzz_xy_1[j];

                    ta_xzzz_xxz_0[j] = pa_x[j] * ta_zzz_xxz_0[j] - pc_x[j] * ta_zzz_xxz_1[j] + fl1_fx * ta_zzz_xz_0[j] - fl1_fx * ta_zzz_xz_1[j];

                    ta_xzzz_xyy_0[j] = pa_x[j] * ta_zzz_xyy_0[j] - pc_x[j] * ta_zzz_xyy_1[j] + 0.5 * fl1_fx * ta_zzz_yy_0[j] - 0.5 * fl1_fx * ta_zzz_yy_1[j];

                    ta_xzzz_xyz_0[j] = pa_x[j] * ta_zzz_xyz_0[j] - pc_x[j] * ta_zzz_xyz_1[j] + 0.5 * fl1_fx * ta_zzz_yz_0[j] - 0.5 * fl1_fx * ta_zzz_yz_1[j];

                    ta_xzzz_xzz_0[j] = pa_x[j] * ta_zzz_xzz_0[j] - pc_x[j] * ta_zzz_xzz_1[j] + 0.5 * fl1_fx * ta_zzz_zz_0[j] - 0.5 * fl1_fx * ta_zzz_zz_1[j];

                    ta_xzzz_yyy_0[j] = pa_x[j] * ta_zzz_yyy_0[j] - pc_x[j] * ta_zzz_yyy_1[j];

                    ta_xzzz_yyz_0[j] = pa_x[j] * ta_zzz_yyz_0[j] - pc_x[j] * ta_zzz_yyz_1[j];

                    ta_xzzz_yzz_0[j] = pa_x[j] * ta_zzz_yzz_0[j] - pc_x[j] * ta_zzz_yzz_1[j];

                    ta_xzzz_zzz_0[j] = pa_x[j] * ta_zzz_zzz_0[j] - pc_x[j] * ta_zzz_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compNuclearPotentialForGF_100_150(      CMemBlock2D<double>& primBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& osFactors,
                                      const CMemBlock2D<double>& paDistances,
                                      const CMemBlock2D<double>& pcDistances,
                                      const CGtoBlock&           braGtoBlock,
                                      const CGtoBlock&           ketGtoBlock,
                                      const int32_t              iContrGto)
    {
        // Batch of Integrals (100,150)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Nuclear Potential"},
                                             {4, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_a_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {4, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_a_4_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_a_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_3_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_a_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to tensors product of distances R(PA) = P - A

                auto pa_y = paDistances.data(3 * idx + 1);

                auto pa_z = paDistances.data(3 * idx + 2);

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_y = pcDistances.data(3 * idx + 1);

                auto pc_z = pcDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

                auto ta_yyy_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 60); 

                auto ta_yyy_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 61); 

                auto ta_yyy_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 62); 

                auto ta_yyy_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 63); 

                auto ta_yyy_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 64); 

                auto ta_yyy_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 65); 

                auto ta_yyy_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 66); 

                auto ta_yyy_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 67); 

                auto ta_yyy_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 68); 

                auto ta_yyy_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 69); 

                auto ta_yyz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 70); 

                auto ta_yyz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 71); 

                auto ta_yyz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 72); 

                auto ta_yyz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 73); 

                auto ta_yyz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 74); 

                auto ta_yyz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 75); 

                auto ta_yyz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 76); 

                auto ta_yyz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 77); 

                auto ta_yyz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 78); 

                auto ta_yyz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 79); 

                auto ta_yzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 80); 

                auto ta_yzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 81); 

                auto ta_yzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 82); 

                auto ta_yzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 83); 

                auto ta_yzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 84); 

                auto ta_yzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 85); 

                auto ta_yzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 86); 

                auto ta_yzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 87); 

                auto ta_yzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 88); 

                auto ta_yzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 89); 

                auto ta_zzz_xxx_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 90); 

                auto ta_zzz_xxy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 91); 

                auto ta_zzz_xxz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 92); 

                auto ta_zzz_xyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 93); 

                auto ta_zzz_xyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 94); 

                auto ta_zzz_xzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 95); 

                auto ta_zzz_yyy_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 96); 

                auto ta_zzz_yyz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 97); 

                auto ta_zzz_yzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 98); 

                auto ta_zzz_zzz_0 = primBuffer.data(pidx_a_3_3_m0 + 100 * idx + 99); 

                auto ta_yyy_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 60); 

                auto ta_yyy_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 61); 

                auto ta_yyy_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 62); 

                auto ta_yyy_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 63); 

                auto ta_yyy_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 64); 

                auto ta_yyy_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 65); 

                auto ta_yyy_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 66); 

                auto ta_yyy_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 67); 

                auto ta_yyy_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 68); 

                auto ta_yyy_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 69); 

                auto ta_yyz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 70); 

                auto ta_yyz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 71); 

                auto ta_yyz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 72); 

                auto ta_yyz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 73); 

                auto ta_yyz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 74); 

                auto ta_yyz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 75); 

                auto ta_yyz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 76); 

                auto ta_yyz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 77); 

                auto ta_yyz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 78); 

                auto ta_yyz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 79); 

                auto ta_yzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 80); 

                auto ta_yzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 81); 

                auto ta_yzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 82); 

                auto ta_yzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 83); 

                auto ta_yzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 84); 

                auto ta_yzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 85); 

                auto ta_yzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 86); 

                auto ta_yzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 87); 

                auto ta_yzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 88); 

                auto ta_yzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 89); 

                auto ta_zzz_xxx_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 90); 

                auto ta_zzz_xxy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 91); 

                auto ta_zzz_xxz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 92); 

                auto ta_zzz_xyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 93); 

                auto ta_zzz_xyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 94); 

                auto ta_zzz_xzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 95); 

                auto ta_zzz_yyy_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 96); 

                auto ta_zzz_yyz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 97); 

                auto ta_zzz_yzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 98); 

                auto ta_zzz_zzz_1 = primBuffer.data(pidx_a_3_3_m1 + 100 * idx + 99); 

                auto ta_yy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 30); 

                auto ta_yy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 31); 

                auto ta_yy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 32); 

                auto ta_yy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 33); 

                auto ta_yy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 34); 

                auto ta_yy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 35); 

                auto ta_yy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 36); 

                auto ta_yy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 37); 

                auto ta_yy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 38); 

                auto ta_yy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 39); 

                auto ta_yz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 40); 

                auto ta_yz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 41); 

                auto ta_yz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 42); 

                auto ta_yz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 43); 

                auto ta_yz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 44); 

                auto ta_yz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 45); 

                auto ta_yz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 46); 

                auto ta_yz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 47); 

                auto ta_yz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 48); 

                auto ta_yz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 49); 

                auto ta_zz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 50); 

                auto ta_zz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 51); 

                auto ta_zz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 52); 

                auto ta_zz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 53); 

                auto ta_zz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 54); 

                auto ta_zz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 55); 

                auto ta_zz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 56); 

                auto ta_zz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 57); 

                auto ta_zz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 58); 

                auto ta_zz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 59); 

                auto ta_yy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 30); 

                auto ta_yy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 31); 

                auto ta_yy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 32); 

                auto ta_yy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 33); 

                auto ta_yy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 34); 

                auto ta_yy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 35); 

                auto ta_yy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 36); 

                auto ta_yy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 37); 

                auto ta_yy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 38); 

                auto ta_yy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 39); 

                auto ta_yz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 40); 

                auto ta_yz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 41); 

                auto ta_yz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 42); 

                auto ta_yz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 43); 

                auto ta_yz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 44); 

                auto ta_yz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 45); 

                auto ta_yz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 46); 

                auto ta_yz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 47); 

                auto ta_yz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 48); 

                auto ta_yz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 49); 

                auto ta_zz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 50); 

                auto ta_zz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 51); 

                auto ta_zz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 52); 

                auto ta_zz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 53); 

                auto ta_zz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 54); 

                auto ta_zz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 55); 

                auto ta_zz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 56); 

                auto ta_zz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 57); 

                auto ta_zz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 58); 

                auto ta_zz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 59); 

                auto ta_yyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 36); 

                auto ta_yyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 37); 

                auto ta_yyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 38); 

                auto ta_yyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 39); 

                auto ta_yyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 40); 

                auto ta_yyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 41); 

                auto ta_yyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 42); 

                auto ta_yyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 43); 

                auto ta_yyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 44); 

                auto ta_yyz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 45); 

                auto ta_yyz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 46); 

                auto ta_yyz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 47); 

                auto ta_yzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 48); 

                auto ta_yzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 49); 

                auto ta_yzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 50); 

                auto ta_yzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 51); 

                auto ta_yzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 52); 

                auto ta_yzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 53); 

                auto ta_zzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 54); 

                auto ta_zzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 55); 

                auto ta_zzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 56); 

                auto ta_zzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 57); 

                auto ta_zzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 58); 

                auto ta_zzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 59); 

                auto ta_yyy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 36); 

                auto ta_yyy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 37); 

                auto ta_yyy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 38); 

                auto ta_yyy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 39); 

                auto ta_yyy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 40); 

                auto ta_yyy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 41); 

                auto ta_yyz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 42); 

                auto ta_yyz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 43); 

                auto ta_yyz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 44); 

                auto ta_yyz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 45); 

                auto ta_yyz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 46); 

                auto ta_yyz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 47); 

                auto ta_yzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 48); 

                auto ta_yzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 49); 

                auto ta_yzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 50); 

                auto ta_yzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 51); 

                auto ta_yzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 52); 

                auto ta_yzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 53); 

                auto ta_zzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 54); 

                auto ta_zzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 55); 

                auto ta_zzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 56); 

                auto ta_zzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 57); 

                auto ta_zzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 58); 

                auto ta_zzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 59); 

                // set up pointers to integrals

                auto ta_yyyy_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 100); 

                auto ta_yyyy_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 101); 

                auto ta_yyyy_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 102); 

                auto ta_yyyy_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 103); 

                auto ta_yyyy_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 104); 

                auto ta_yyyy_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 105); 

                auto ta_yyyy_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 106); 

                auto ta_yyyy_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 107); 

                auto ta_yyyy_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 108); 

                auto ta_yyyy_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 109); 

                auto ta_yyyz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 110); 

                auto ta_yyyz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 111); 

                auto ta_yyyz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 112); 

                auto ta_yyyz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 113); 

                auto ta_yyyz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 114); 

                auto ta_yyyz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 115); 

                auto ta_yyyz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 116); 

                auto ta_yyyz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 117); 

                auto ta_yyyz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 118); 

                auto ta_yyyz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 119); 

                auto ta_yyzz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 120); 

                auto ta_yyzz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 121); 

                auto ta_yyzz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 122); 

                auto ta_yyzz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 123); 

                auto ta_yyzz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 124); 

                auto ta_yyzz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 125); 

                auto ta_yyzz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 126); 

                auto ta_yyzz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 127); 

                auto ta_yyzz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 128); 

                auto ta_yyzz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 129); 

                auto ta_yzzz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 130); 

                auto ta_yzzz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 131); 

                auto ta_yzzz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 132); 

                auto ta_yzzz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 133); 

                auto ta_yzzz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 134); 

                auto ta_yzzz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 135); 

                auto ta_yzzz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 136); 

                auto ta_yzzz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 137); 

                auto ta_yzzz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 138); 

                auto ta_yzzz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 139); 

                auto ta_zzzz_xxx_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 140); 

                auto ta_zzzz_xxy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 141); 

                auto ta_zzzz_xxz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 142); 

                auto ta_zzzz_xyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 143); 

                auto ta_zzzz_xyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 144); 

                auto ta_zzzz_xzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 145); 

                auto ta_zzzz_yyy_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 146); 

                auto ta_zzzz_yyz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 147); 

                auto ta_zzzz_yzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 148); 

                auto ta_zzzz_zzz_0 = primBuffer.data(pidx_a_4_3_m0 + 150 * idx + 149); 

                // Batch of Integrals (100,150)

                #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_yy_xxx_0, ta_yy_xxx_1, ta_yy_xxy_0, \
                                         ta_yy_xxy_1, ta_yy_xxz_0, ta_yy_xxz_1, ta_yy_xyy_0, ta_yy_xyy_1, ta_yy_xyz_0, \
                                         ta_yy_xyz_1, ta_yy_xzz_0, ta_yy_xzz_1, ta_yy_yyy_0, ta_yy_yyy_1, ta_yy_yyz_0, \
                                         ta_yy_yyz_1, ta_yy_yzz_0, ta_yy_yzz_1, ta_yy_zzz_0, ta_yy_zzz_1, ta_yyy_xx_0, \
                                         ta_yyy_xx_1, ta_yyy_xxx_0, ta_yyy_xxx_1, ta_yyy_xxy_0, ta_yyy_xxy_1, ta_yyy_xxz_0, \
                                         ta_yyy_xxz_1, ta_yyy_xy_0, ta_yyy_xy_1, ta_yyy_xyy_0, ta_yyy_xyy_1, ta_yyy_xyz_0, \
                                         ta_yyy_xyz_1, ta_yyy_xz_0, ta_yyy_xz_1, ta_yyy_xzz_0, ta_yyy_xzz_1, ta_yyy_yy_0, \
                                         ta_yyy_yy_1, ta_yyy_yyy_0, ta_yyy_yyy_1, ta_yyy_yyz_0, ta_yyy_yyz_1, ta_yyy_yz_0, \
                                         ta_yyy_yz_1, ta_yyy_yzz_0, ta_yyy_yzz_1, ta_yyy_zz_0, ta_yyy_zz_1, ta_yyy_zzz_0, \
                                         ta_yyy_zzz_1, ta_yyyy_xxx_0, ta_yyyy_xxy_0, ta_yyyy_xxz_0, ta_yyyy_xyy_0, \
                                         ta_yyyy_xyz_0, ta_yyyy_xzz_0, ta_yyyy_yyy_0, ta_yyyy_yyz_0, ta_yyyy_yzz_0, \
                                         ta_yyyy_zzz_0, ta_yyyz_xxx_0, ta_yyyz_xxy_0, ta_yyyz_xxz_0, ta_yyyz_xyy_0, \
                                         ta_yyyz_xyz_0, ta_yyyz_xzz_0, ta_yyyz_yyy_0, ta_yyyz_yyz_0, ta_yyyz_yzz_0, \
                                         ta_yyyz_zzz_0, ta_yyz_xx_0, ta_yyz_xx_1, ta_yyz_xxx_0, ta_yyz_xxx_1, ta_yyz_xxy_0, \
                                         ta_yyz_xxy_1, ta_yyz_xxz_0, ta_yyz_xxz_1, ta_yyz_xy_0, ta_yyz_xy_1, ta_yyz_xyy_0, \
                                         ta_yyz_xyy_1, ta_yyz_xyz_0, ta_yyz_xyz_1, ta_yyz_xz_0, ta_yyz_xz_1, ta_yyz_xzz_0, \
                                         ta_yyz_xzz_1, ta_yyz_yy_0, ta_yyz_yy_1, ta_yyz_yyy_0, ta_yyz_yyy_1, ta_yyz_yyz_0, \
                                         ta_yyz_yyz_1, ta_yyz_yz_0, ta_yyz_yz_1, ta_yyz_yzz_0, ta_yyz_yzz_1, ta_yyz_zz_0, \
                                         ta_yyz_zz_1, ta_yyz_zzz_0, ta_yyz_zzz_1, ta_yyzz_xxx_0, ta_yyzz_xxy_0, \
                                         ta_yyzz_xxz_0, ta_yyzz_xyy_0, ta_yyzz_xyz_0, ta_yyzz_xzz_0, ta_yyzz_yyy_0, \
                                         ta_yyzz_yyz_0, ta_yyzz_yzz_0, ta_yyzz_zzz_0, ta_yz_xxx_0, ta_yz_xxx_1, ta_yz_xxy_0, \
                                         ta_yz_xxy_1, ta_yz_xxz_0, ta_yz_xxz_1, ta_yz_xyy_0, ta_yz_xyy_1, ta_yz_xyz_0, \
                                         ta_yz_xyz_1, ta_yz_xzz_0, ta_yz_xzz_1, ta_yz_yyy_0, ta_yz_yyy_1, ta_yz_yyz_0, \
                                         ta_yz_yyz_1, ta_yz_yzz_0, ta_yz_yzz_1, ta_yz_zzz_0, ta_yz_zzz_1, ta_yzz_xx_0, \
                                         ta_yzz_xx_1, ta_yzz_xxx_0, ta_yzz_xxx_1, ta_yzz_xxy_0, ta_yzz_xxy_1, ta_yzz_xxz_0, \
                                         ta_yzz_xxz_1, ta_yzz_xy_0, ta_yzz_xy_1, ta_yzz_xyy_0, ta_yzz_xyy_1, ta_yzz_xyz_0, \
                                         ta_yzz_xyz_1, ta_yzz_xz_0, ta_yzz_xz_1, ta_yzz_xzz_0, ta_yzz_xzz_1, ta_yzz_yy_0, \
                                         ta_yzz_yy_1, ta_yzz_yyy_0, ta_yzz_yyy_1, ta_yzz_yyz_0, ta_yzz_yyz_1, ta_yzz_yz_0, \
                                         ta_yzz_yz_1, ta_yzz_yzz_0, ta_yzz_yzz_1, ta_yzz_zz_0, ta_yzz_zz_1, ta_yzz_zzz_0, \
                                         ta_yzz_zzz_1, ta_yzzz_xxx_0, ta_yzzz_xxy_0, ta_yzzz_xxz_0, ta_yzzz_xyy_0, \
                                         ta_yzzz_xyz_0, ta_yzzz_xzz_0, ta_yzzz_yyy_0, ta_yzzz_yyz_0, ta_yzzz_yzz_0, \
                                         ta_yzzz_zzz_0, ta_zz_xxx_0, ta_zz_xxx_1, ta_zz_xxy_0, ta_zz_xxy_1, ta_zz_xxz_0, \
                                         ta_zz_xxz_1, ta_zz_xyy_0, ta_zz_xyy_1, ta_zz_xyz_0, ta_zz_xyz_1, ta_zz_xzz_0, \
                                         ta_zz_xzz_1, ta_zz_yyy_0, ta_zz_yyy_1, ta_zz_yyz_0, ta_zz_yyz_1, ta_zz_yzz_0, \
                                         ta_zz_yzz_1, ta_zz_zzz_0, ta_zz_zzz_1, ta_zzz_xx_0, ta_zzz_xx_1, ta_zzz_xxx_0, \
                                         ta_zzz_xxx_1, ta_zzz_xxy_0, ta_zzz_xxy_1, ta_zzz_xxz_0, ta_zzz_xxz_1, ta_zzz_xy_0, \
                                         ta_zzz_xy_1, ta_zzz_xyy_0, ta_zzz_xyy_1, ta_zzz_xyz_0, ta_zzz_xyz_1, ta_zzz_xz_0, \
                                         ta_zzz_xz_1, ta_zzz_xzz_0, ta_zzz_xzz_1, ta_zzz_yy_0, ta_zzz_yy_1, ta_zzz_yyy_0, \
                                         ta_zzz_yyy_1, ta_zzz_yyz_0, ta_zzz_yyz_1, ta_zzz_yz_0, ta_zzz_yz_1, ta_zzz_yzz_0, \
                                         ta_zzz_yzz_1, ta_zzz_zz_0, ta_zzz_zz_1, ta_zzz_zzz_0, ta_zzz_zzz_1, ta_zzzz_xxx_0, \
                                         ta_zzzz_xxy_0, ta_zzzz_xxz_0, ta_zzzz_xyy_0, ta_zzzz_xyz_0, ta_zzzz_xzz_0, \
                                         ta_zzzz_yyy_0, ta_zzzz_yyz_0, ta_zzzz_yzz_0, ta_zzzz_zzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    ta_yyyy_xxx_0[j] = pa_y[j] * ta_yyy_xxx_0[j] - pc_y[j] * ta_yyy_xxx_1[j] + 1.5 * fl1_fx * ta_yy_xxx_0[j] - 1.5 * fl1_fx * ta_yy_xxx_1[j];

                    ta_yyyy_xxy_0[j] = pa_y[j] * ta_yyy_xxy_0[j] - pc_y[j] * ta_yyy_xxy_1[j] + 1.5 * fl1_fx * ta_yy_xxy_0[j] - 1.5 * fl1_fx * ta_yy_xxy_1[j] + 0.5 * fl1_fx * ta_yyy_xx_0[j] - 0.5 * fl1_fx * ta_yyy_xx_1[j];

                    ta_yyyy_xxz_0[j] = pa_y[j] * ta_yyy_xxz_0[j] - pc_y[j] * ta_yyy_xxz_1[j] + 1.5 * fl1_fx * ta_yy_xxz_0[j] - 1.5 * fl1_fx * ta_yy_xxz_1[j];

                    ta_yyyy_xyy_0[j] = pa_y[j] * ta_yyy_xyy_0[j] - pc_y[j] * ta_yyy_xyy_1[j] + 1.5 * fl1_fx * ta_yy_xyy_0[j] - 1.5 * fl1_fx * ta_yy_xyy_1[j] + fl1_fx * ta_yyy_xy_0[j] - fl1_fx * ta_yyy_xy_1[j];

                    ta_yyyy_xyz_0[j] = pa_y[j] * ta_yyy_xyz_0[j] - pc_y[j] * ta_yyy_xyz_1[j] + 1.5 * fl1_fx * ta_yy_xyz_0[j] - 1.5 * fl1_fx * ta_yy_xyz_1[j] + 0.5 * fl1_fx * ta_yyy_xz_0[j] - 0.5 * fl1_fx * ta_yyy_xz_1[j];

                    ta_yyyy_xzz_0[j] = pa_y[j] * ta_yyy_xzz_0[j] - pc_y[j] * ta_yyy_xzz_1[j] + 1.5 * fl1_fx * ta_yy_xzz_0[j] - 1.5 * fl1_fx * ta_yy_xzz_1[j];

                    ta_yyyy_yyy_0[j] = pa_y[j] * ta_yyy_yyy_0[j] - pc_y[j] * ta_yyy_yyy_1[j] + 1.5 * fl1_fx * ta_yy_yyy_0[j] - 1.5 * fl1_fx * ta_yy_yyy_1[j] + 1.5 * fl1_fx * ta_yyy_yy_0[j] - 1.5 * fl1_fx * ta_yyy_yy_1[j];

                    ta_yyyy_yyz_0[j] = pa_y[j] * ta_yyy_yyz_0[j] - pc_y[j] * ta_yyy_yyz_1[j] + 1.5 * fl1_fx * ta_yy_yyz_0[j] - 1.5 * fl1_fx * ta_yy_yyz_1[j] + fl1_fx * ta_yyy_yz_0[j] - fl1_fx * ta_yyy_yz_1[j];

                    ta_yyyy_yzz_0[j] = pa_y[j] * ta_yyy_yzz_0[j] - pc_y[j] * ta_yyy_yzz_1[j] + 1.5 * fl1_fx * ta_yy_yzz_0[j] - 1.5 * fl1_fx * ta_yy_yzz_1[j] + 0.5 * fl1_fx * ta_yyy_zz_0[j] - 0.5 * fl1_fx * ta_yyy_zz_1[j];

                    ta_yyyy_zzz_0[j] = pa_y[j] * ta_yyy_zzz_0[j] - pc_y[j] * ta_yyy_zzz_1[j] + 1.5 * fl1_fx * ta_yy_zzz_0[j] - 1.5 * fl1_fx * ta_yy_zzz_1[j];

                    ta_yyyz_xxx_0[j] = pa_y[j] * ta_yyz_xxx_0[j] - pc_y[j] * ta_yyz_xxx_1[j] + fl1_fx * ta_yz_xxx_0[j] - fl1_fx * ta_yz_xxx_1[j];

                    ta_yyyz_xxy_0[j] = pa_y[j] * ta_yyz_xxy_0[j] - pc_y[j] * ta_yyz_xxy_1[j] + fl1_fx * ta_yz_xxy_0[j] - fl1_fx * ta_yz_xxy_1[j] + 0.5 * fl1_fx * ta_yyz_xx_0[j] - 0.5 * fl1_fx * ta_yyz_xx_1[j];

                    ta_yyyz_xxz_0[j] = pa_y[j] * ta_yyz_xxz_0[j] - pc_y[j] * ta_yyz_xxz_1[j] + fl1_fx * ta_yz_xxz_0[j] - fl1_fx * ta_yz_xxz_1[j];

                    ta_yyyz_xyy_0[j] = pa_y[j] * ta_yyz_xyy_0[j] - pc_y[j] * ta_yyz_xyy_1[j] + fl1_fx * ta_yz_xyy_0[j] - fl1_fx * ta_yz_xyy_1[j] + fl1_fx * ta_yyz_xy_0[j] - fl1_fx * ta_yyz_xy_1[j];

                    ta_yyyz_xyz_0[j] = pa_y[j] * ta_yyz_xyz_0[j] - pc_y[j] * ta_yyz_xyz_1[j] + fl1_fx * ta_yz_xyz_0[j] - fl1_fx * ta_yz_xyz_1[j] + 0.5 * fl1_fx * ta_yyz_xz_0[j] - 0.5 * fl1_fx * ta_yyz_xz_1[j];

                    ta_yyyz_xzz_0[j] = pa_y[j] * ta_yyz_xzz_0[j] - pc_y[j] * ta_yyz_xzz_1[j] + fl1_fx * ta_yz_xzz_0[j] - fl1_fx * ta_yz_xzz_1[j];

                    ta_yyyz_yyy_0[j] = pa_y[j] * ta_yyz_yyy_0[j] - pc_y[j] * ta_yyz_yyy_1[j] + fl1_fx * ta_yz_yyy_0[j] - fl1_fx * ta_yz_yyy_1[j] + 1.5 * fl1_fx * ta_yyz_yy_0[j] - 1.5 * fl1_fx * ta_yyz_yy_1[j];

                    ta_yyyz_yyz_0[j] = pa_y[j] * ta_yyz_yyz_0[j] - pc_y[j] * ta_yyz_yyz_1[j] + fl1_fx * ta_yz_yyz_0[j] - fl1_fx * ta_yz_yyz_1[j] + fl1_fx * ta_yyz_yz_0[j] - fl1_fx * ta_yyz_yz_1[j];

                    ta_yyyz_yzz_0[j] = pa_y[j] * ta_yyz_yzz_0[j] - pc_y[j] * ta_yyz_yzz_1[j] + fl1_fx * ta_yz_yzz_0[j] - fl1_fx * ta_yz_yzz_1[j] + 0.5 * fl1_fx * ta_yyz_zz_0[j] - 0.5 * fl1_fx * ta_yyz_zz_1[j];

                    ta_yyyz_zzz_0[j] = pa_y[j] * ta_yyz_zzz_0[j] - pc_y[j] * ta_yyz_zzz_1[j] + fl1_fx * ta_yz_zzz_0[j] - fl1_fx * ta_yz_zzz_1[j];

                    ta_yyzz_xxx_0[j] = pa_y[j] * ta_yzz_xxx_0[j] - pc_y[j] * ta_yzz_xxx_1[j] + 0.5 * fl1_fx * ta_zz_xxx_0[j] - 0.5 * fl1_fx * ta_zz_xxx_1[j];

                    ta_yyzz_xxy_0[j] = pa_y[j] * ta_yzz_xxy_0[j] - pc_y[j] * ta_yzz_xxy_1[j] + 0.5 * fl1_fx * ta_zz_xxy_0[j] - 0.5 * fl1_fx * ta_zz_xxy_1[j] + 0.5 * fl1_fx * ta_yzz_xx_0[j] - 0.5 * fl1_fx * ta_yzz_xx_1[j];

                    ta_yyzz_xxz_0[j] = pa_y[j] * ta_yzz_xxz_0[j] - pc_y[j] * ta_yzz_xxz_1[j] + 0.5 * fl1_fx * ta_zz_xxz_0[j] - 0.5 * fl1_fx * ta_zz_xxz_1[j];

                    ta_yyzz_xyy_0[j] = pa_y[j] * ta_yzz_xyy_0[j] - pc_y[j] * ta_yzz_xyy_1[j] + 0.5 * fl1_fx * ta_zz_xyy_0[j] - 0.5 * fl1_fx * ta_zz_xyy_1[j] + fl1_fx * ta_yzz_xy_0[j] - fl1_fx * ta_yzz_xy_1[j];

                    ta_yyzz_xyz_0[j] = pa_y[j] * ta_yzz_xyz_0[j] - pc_y[j] * ta_yzz_xyz_1[j] + 0.5 * fl1_fx * ta_zz_xyz_0[j] - 0.5 * fl1_fx * ta_zz_xyz_1[j] + 0.5 * fl1_fx * ta_yzz_xz_0[j] - 0.5 * fl1_fx * ta_yzz_xz_1[j];

                    ta_yyzz_xzz_0[j] = pa_y[j] * ta_yzz_xzz_0[j] - pc_y[j] * ta_yzz_xzz_1[j] + 0.5 * fl1_fx * ta_zz_xzz_0[j] - 0.5 * fl1_fx * ta_zz_xzz_1[j];

                    ta_yyzz_yyy_0[j] = pa_y[j] * ta_yzz_yyy_0[j] - pc_y[j] * ta_yzz_yyy_1[j] + 0.5 * fl1_fx * ta_zz_yyy_0[j] - 0.5 * fl1_fx * ta_zz_yyy_1[j] + 1.5 * fl1_fx * ta_yzz_yy_0[j] - 1.5 * fl1_fx * ta_yzz_yy_1[j];

                    ta_yyzz_yyz_0[j] = pa_y[j] * ta_yzz_yyz_0[j] - pc_y[j] * ta_yzz_yyz_1[j] + 0.5 * fl1_fx * ta_zz_yyz_0[j] - 0.5 * fl1_fx * ta_zz_yyz_1[j] + fl1_fx * ta_yzz_yz_0[j] - fl1_fx * ta_yzz_yz_1[j];

                    ta_yyzz_yzz_0[j] = pa_y[j] * ta_yzz_yzz_0[j] - pc_y[j] * ta_yzz_yzz_1[j] + 0.5 * fl1_fx * ta_zz_yzz_0[j] - 0.5 * fl1_fx * ta_zz_yzz_1[j] + 0.5 * fl1_fx * ta_yzz_zz_0[j] - 0.5 * fl1_fx * ta_yzz_zz_1[j];

                    ta_yyzz_zzz_0[j] = pa_y[j] * ta_yzz_zzz_0[j] - pc_y[j] * ta_yzz_zzz_1[j] + 0.5 * fl1_fx * ta_zz_zzz_0[j] - 0.5 * fl1_fx * ta_zz_zzz_1[j];

                    ta_yzzz_xxx_0[j] = pa_y[j] * ta_zzz_xxx_0[j] - pc_y[j] * ta_zzz_xxx_1[j];

                    ta_yzzz_xxy_0[j] = pa_y[j] * ta_zzz_xxy_0[j] - pc_y[j] * ta_zzz_xxy_1[j] + 0.5 * fl1_fx * ta_zzz_xx_0[j] - 0.5 * fl1_fx * ta_zzz_xx_1[j];

                    ta_yzzz_xxz_0[j] = pa_y[j] * ta_zzz_xxz_0[j] - pc_y[j] * ta_zzz_xxz_1[j];

                    ta_yzzz_xyy_0[j] = pa_y[j] * ta_zzz_xyy_0[j] - pc_y[j] * ta_zzz_xyy_1[j] + fl1_fx * ta_zzz_xy_0[j] - fl1_fx * ta_zzz_xy_1[j];

                    ta_yzzz_xyz_0[j] = pa_y[j] * ta_zzz_xyz_0[j] - pc_y[j] * ta_zzz_xyz_1[j] + 0.5 * fl1_fx * ta_zzz_xz_0[j] - 0.5 * fl1_fx * ta_zzz_xz_1[j];

                    ta_yzzz_xzz_0[j] = pa_y[j] * ta_zzz_xzz_0[j] - pc_y[j] * ta_zzz_xzz_1[j];

                    ta_yzzz_yyy_0[j] = pa_y[j] * ta_zzz_yyy_0[j] - pc_y[j] * ta_zzz_yyy_1[j] + 1.5 * fl1_fx * ta_zzz_yy_0[j] - 1.5 * fl1_fx * ta_zzz_yy_1[j];

                    ta_yzzz_yyz_0[j] = pa_y[j] * ta_zzz_yyz_0[j] - pc_y[j] * ta_zzz_yyz_1[j] + fl1_fx * ta_zzz_yz_0[j] - fl1_fx * ta_zzz_yz_1[j];

                    ta_yzzz_yzz_0[j] = pa_y[j] * ta_zzz_yzz_0[j] - pc_y[j] * ta_zzz_yzz_1[j] + 0.5 * fl1_fx * ta_zzz_zz_0[j] - 0.5 * fl1_fx * ta_zzz_zz_1[j];

                    ta_yzzz_zzz_0[j] = pa_y[j] * ta_zzz_zzz_0[j] - pc_y[j] * ta_zzz_zzz_1[j];

                    ta_zzzz_xxx_0[j] = pa_z[j] * ta_zzz_xxx_0[j] - pc_z[j] * ta_zzz_xxx_1[j] + 1.5 * fl1_fx * ta_zz_xxx_0[j] - 1.5 * fl1_fx * ta_zz_xxx_1[j];

                    ta_zzzz_xxy_0[j] = pa_z[j] * ta_zzz_xxy_0[j] - pc_z[j] * ta_zzz_xxy_1[j] + 1.5 * fl1_fx * ta_zz_xxy_0[j] - 1.5 * fl1_fx * ta_zz_xxy_1[j];

                    ta_zzzz_xxz_0[j] = pa_z[j] * ta_zzz_xxz_0[j] - pc_z[j] * ta_zzz_xxz_1[j] + 1.5 * fl1_fx * ta_zz_xxz_0[j] - 1.5 * fl1_fx * ta_zz_xxz_1[j] + 0.5 * fl1_fx * ta_zzz_xx_0[j] - 0.5 * fl1_fx * ta_zzz_xx_1[j];

                    ta_zzzz_xyy_0[j] = pa_z[j] * ta_zzz_xyy_0[j] - pc_z[j] * ta_zzz_xyy_1[j] + 1.5 * fl1_fx * ta_zz_xyy_0[j] - 1.5 * fl1_fx * ta_zz_xyy_1[j];

                    ta_zzzz_xyz_0[j] = pa_z[j] * ta_zzz_xyz_0[j] - pc_z[j] * ta_zzz_xyz_1[j] + 1.5 * fl1_fx * ta_zz_xyz_0[j] - 1.5 * fl1_fx * ta_zz_xyz_1[j] + 0.5 * fl1_fx * ta_zzz_xy_0[j] - 0.5 * fl1_fx * ta_zzz_xy_1[j];

                    ta_zzzz_xzz_0[j] = pa_z[j] * ta_zzz_xzz_0[j] - pc_z[j] * ta_zzz_xzz_1[j] + 1.5 * fl1_fx * ta_zz_xzz_0[j] - 1.5 * fl1_fx * ta_zz_xzz_1[j] + fl1_fx * ta_zzz_xz_0[j] - fl1_fx * ta_zzz_xz_1[j];

                    ta_zzzz_yyy_0[j] = pa_z[j] * ta_zzz_yyy_0[j] - pc_z[j] * ta_zzz_yyy_1[j] + 1.5 * fl1_fx * ta_zz_yyy_0[j] - 1.5 * fl1_fx * ta_zz_yyy_1[j];

                    ta_zzzz_yyz_0[j] = pa_z[j] * ta_zzz_yyz_0[j] - pc_z[j] * ta_zzz_yyz_1[j] + 1.5 * fl1_fx * ta_zz_yyz_0[j] - 1.5 * fl1_fx * ta_zz_yyz_1[j] + 0.5 * fl1_fx * ta_zzz_yy_0[j] - 0.5 * fl1_fx * ta_zzz_yy_1[j];

                    ta_zzzz_yzz_0[j] = pa_z[j] * ta_zzz_yzz_0[j] - pc_z[j] * ta_zzz_yzz_1[j] + 1.5 * fl1_fx * ta_zz_yzz_0[j] - 1.5 * fl1_fx * ta_zz_yzz_1[j] + fl1_fx * ta_zzz_yz_0[j] - fl1_fx * ta_zzz_yz_1[j];

                    ta_zzzz_zzz_0[j] = pa_z[j] * ta_zzz_zzz_0[j] - pc_z[j] * ta_zzz_zzz_1[j] + 1.5 * fl1_fx * ta_zz_zzz_0[j] - 1.5 * fl1_fx * ta_zz_zzz_1[j] + 1.5 * fl1_fx * ta_zzz_zz_0[j] - 1.5 * fl1_fx * ta_zzz_zz_1[j];
                }

                idx++;
            }
        }
    }


} // npotrecfunc namespace

