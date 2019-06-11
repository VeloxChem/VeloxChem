//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldRecFuncForFG.hpp"

namespace efieldrecfunc { // efieldrecfunc namespace

    void
    compElectricFieldForFG(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pcDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        efieldrecfunc::compElectricFieldForFG_0_50(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   pcDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        efieldrecfunc::compElectricFieldForFG_50_100(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     pcDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        efieldrecfunc::compElectricFieldForFG_100_150(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForFG_150_200(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForFG_200_250(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForFG_250_300(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForFG_300_350(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForFG_350_400(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForFG_400_450(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compElectricFieldForFG_0_50(      CMemBlock2D<double>& primBuffer,
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

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tex_xx_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx); 

                auto tey_xx_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx); 

                auto tez_xx_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx); 

                auto tex_xx_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 1); 

                auto tey_xx_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 1); 

                auto tez_xx_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 1); 

                auto tex_xx_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 2); 

                auto tey_xx_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 2); 

                auto tez_xx_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 2); 

                auto tex_xx_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 3); 

                auto tey_xx_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 3); 

                auto tez_xx_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 3); 

                auto tex_xx_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 4); 

                auto tey_xx_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 4); 

                auto tez_xx_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 4); 

                auto tex_xx_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 5); 

                auto tey_xx_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 5); 

                auto tez_xx_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 5); 

                auto tex_xx_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 6); 

                auto tey_xx_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 6); 

                auto tez_xx_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 6); 

                auto tex_xx_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 7); 

                auto tey_xx_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 7); 

                auto tez_xx_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 7); 

                auto tex_xx_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 8); 

                auto tey_xx_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 8); 

                auto tez_xx_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 8); 

                auto tex_xx_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 9); 

                auto tey_xx_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 9); 

                auto tez_xx_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 9); 

                auto tex_xx_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 10); 

                auto tey_xx_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 10); 

                auto tez_xx_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 10); 

                auto tex_xx_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 11); 

                auto tey_xx_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 11); 

                auto tez_xx_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 11); 

                auto tex_xx_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 12); 

                auto tey_xx_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 12); 

                auto tez_xx_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 12); 

                auto tex_xx_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 13); 

                auto tey_xx_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 13); 

                auto tez_xx_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 13); 

                auto tex_xx_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 14); 

                auto tey_xx_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 14); 

                auto tez_xx_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 14); 

                auto tex_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 15); 

                auto tey_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 15); 

                auto tez_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 15); 

                auto tex_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 16); 

                auto tey_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 16); 

                auto tex_xx_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx); 

                auto tey_xx_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx); 

                auto tez_xx_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx); 

                auto tex_xx_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 1); 

                auto tey_xx_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 1); 

                auto tez_xx_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 1); 

                auto tex_xx_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 2); 

                auto tey_xx_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 2); 

                auto tez_xx_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 2); 

                auto tex_xx_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 3); 

                auto tey_xx_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 3); 

                auto tez_xx_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 3); 

                auto tex_xx_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 4); 

                auto tey_xx_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 4); 

                auto tez_xx_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 4); 

                auto tex_xx_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 5); 

                auto tey_xx_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 5); 

                auto tez_xx_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 5); 

                auto tex_xx_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 6); 

                auto tey_xx_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 6); 

                auto tez_xx_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 6); 

                auto tex_xx_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 7); 

                auto tey_xx_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 7); 

                auto tez_xx_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 7); 

                auto tex_xx_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 8); 

                auto tey_xx_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 8); 

                auto tez_xx_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 8); 

                auto tex_xx_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 9); 

                auto tey_xx_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 9); 

                auto tez_xx_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 9); 

                auto tex_xx_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 10); 

                auto tey_xx_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 10); 

                auto tez_xx_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 10); 

                auto tex_xx_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 11); 

                auto tey_xx_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 11); 

                auto tez_xx_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 11); 

                auto tex_xx_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 12); 

                auto tey_xx_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 12); 

                auto tez_xx_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 12); 

                auto tex_xx_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 13); 

                auto tey_xx_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 13); 

                auto tez_xx_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 13); 

                auto tex_xx_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 14); 

                auto tey_xx_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 14); 

                auto tez_xx_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 14); 

                auto tex_xy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 15); 

                auto tey_xy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 15); 

                auto tez_xy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 15); 

                auto tex_xy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 16); 

                auto tey_xy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 16); 

                auto tex_x_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx); 

                auto tey_x_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx); 

                auto tez_x_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx); 

                auto tex_x_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 1); 

                auto tey_x_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 1); 

                auto tez_x_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 1); 

                auto tex_x_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 2); 

                auto tey_x_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 2); 

                auto tez_x_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 2); 

                auto tex_x_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 3); 

                auto tey_x_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 3); 

                auto tez_x_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 3); 

                auto tex_x_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 4); 

                auto tey_x_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 4); 

                auto tez_x_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 4); 

                auto tex_x_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 5); 

                auto tey_x_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 5); 

                auto tez_x_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 5); 

                auto tex_x_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 6); 

                auto tey_x_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 6); 

                auto tez_x_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 6); 

                auto tex_x_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 7); 

                auto tey_x_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 7); 

                auto tez_x_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 7); 

                auto tex_x_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 8); 

                auto tey_x_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 8); 

                auto tez_x_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 8); 

                auto tex_x_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 9); 

                auto tey_x_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 9); 

                auto tez_x_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 9); 

                auto tex_x_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 10); 

                auto tey_x_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 10); 

                auto tez_x_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 10); 

                auto tex_x_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 11); 

                auto tey_x_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 11); 

                auto tez_x_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 11); 

                auto tex_x_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 12); 

                auto tey_x_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 12); 

                auto tez_x_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 12); 

                auto tex_x_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 13); 

                auto tey_x_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 13); 

                auto tez_x_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 13); 

                auto tex_x_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 14); 

                auto tey_x_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 14); 

                auto tez_x_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 14); 

                auto tex_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 15); 

                auto tey_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 15); 

                auto tez_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 15); 

                auto tex_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 16); 

                auto tey_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 16); 

                auto tex_x_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx); 

                auto tey_x_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx); 

                auto tez_x_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx); 

                auto tex_x_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 1); 

                auto tey_x_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 1); 

                auto tez_x_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 1); 

                auto tex_x_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 2); 

                auto tey_x_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 2); 

                auto tez_x_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 2); 

                auto tex_x_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 3); 

                auto tey_x_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 3); 

                auto tez_x_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 3); 

                auto tex_x_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 4); 

                auto tey_x_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 4); 

                auto tez_x_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 4); 

                auto tex_x_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 5); 

                auto tey_x_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 5); 

                auto tez_x_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 5); 

                auto tex_x_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 6); 

                auto tey_x_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 6); 

                auto tez_x_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 6); 

                auto tex_x_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 7); 

                auto tey_x_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 7); 

                auto tez_x_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 7); 

                auto tex_x_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 8); 

                auto tey_x_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 8); 

                auto tez_x_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 8); 

                auto tex_x_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 9); 

                auto tey_x_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 9); 

                auto tez_x_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 9); 

                auto tex_x_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 10); 

                auto tey_x_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 10); 

                auto tez_x_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 10); 

                auto tex_x_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 11); 

                auto tey_x_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 11); 

                auto tez_x_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 11); 

                auto tex_x_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 12); 

                auto tey_x_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 12); 

                auto tez_x_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 12); 

                auto tex_x_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 13); 

                auto tey_x_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 13); 

                auto tez_x_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 13); 

                auto tex_x_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 14); 

                auto tey_x_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 14); 

                auto tez_x_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 14); 

                auto tex_y_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 15); 

                auto tey_y_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 15); 

                auto tez_y_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 15); 

                auto tex_y_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 16); 

                auto tey_y_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 16); 

                auto tex_xx_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx); 

                auto tey_xx_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx); 

                auto tez_xx_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx); 

                auto tex_xx_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 1); 

                auto tey_xx_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 1); 

                auto tez_xx_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 1); 

                auto tex_xx_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 2); 

                auto tey_xx_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 2); 

                auto tez_xx_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 2); 

                auto tex_xx_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 3); 

                auto tey_xx_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 3); 

                auto tez_xx_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 3); 

                auto tex_xx_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 4); 

                auto tey_xx_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 4); 

                auto tez_xx_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 4); 

                auto tex_xx_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 5); 

                auto tey_xx_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 5); 

                auto tez_xx_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 5); 

                auto tex_xx_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 6); 

                auto tey_xx_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 6); 

                auto tez_xx_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 6); 

                auto tex_xx_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 7); 

                auto tey_xx_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 7); 

                auto tez_xx_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 7); 

                auto tex_xx_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 8); 

                auto tey_xx_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 8); 

                auto tez_xx_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 8); 

                auto tex_xx_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 9); 

                auto tey_xx_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 9); 

                auto tez_xx_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 9); 

                auto tex_xy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 10); 

                auto tey_xy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 10); 

                auto tez_xy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 10); 

                auto tex_xy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 11); 

                auto tey_xy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 11); 

                auto tex_xx_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx); 

                auto tey_xx_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx); 

                auto tez_xx_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx); 

                auto tex_xx_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 1); 

                auto tey_xx_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 1); 

                auto tez_xx_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 1); 

                auto tex_xx_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 2); 

                auto tey_xx_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 2); 

                auto tez_xx_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 2); 

                auto tex_xx_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 3); 

                auto tey_xx_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 3); 

                auto tez_xx_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 3); 

                auto tex_xx_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 4); 

                auto tey_xx_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 4); 

                auto tez_xx_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 4); 

                auto tex_xx_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 5); 

                auto tey_xx_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 5); 

                auto tez_xx_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 5); 

                auto tex_xx_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 6); 

                auto tey_xx_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 6); 

                auto tez_xx_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 6); 

                auto tex_xx_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 7); 

                auto tey_xx_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 7); 

                auto tez_xx_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 7); 

                auto tex_xx_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 8); 

                auto tey_xx_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 8); 

                auto tez_xx_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 8); 

                auto tex_xx_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 9); 

                auto tey_xx_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 9); 

                auto tez_xx_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 9); 

                auto tex_xy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 10); 

                auto tey_xy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 10); 

                auto tez_xy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 10); 

                auto tex_xy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 11); 

                auto tey_xy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 11); 

                auto ta_xx_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx); 

                auto ta_xx_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 1); 

                auto ta_xx_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 2); 

                auto ta_xx_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 3); 

                auto ta_xx_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 4); 

                auto ta_xx_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 5); 

                auto ta_xx_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 6); 

                auto ta_xx_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 7); 

                auto ta_xx_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 8); 

                auto ta_xx_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 9); 

                auto ta_xx_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 10); 

                auto ta_xx_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 11); 

                auto ta_xx_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 12); 

                auto ta_xx_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 13); 

                auto ta_xx_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 14); 

                auto ta_xy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 15); 

                auto ta_xy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 16); 

                // set up pointers to integrals

                auto tex_xxx_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx); 

                auto tey_xxx_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx); 

                auto tez_xxx_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx); 

                auto tex_xxx_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 1); 

                auto tey_xxx_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 1); 

                auto tez_xxx_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 1); 

                auto tex_xxx_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 2); 

                auto tey_xxx_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 2); 

                auto tez_xxx_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 2); 

                auto tex_xxx_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 3); 

                auto tey_xxx_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 3); 

                auto tez_xxx_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 3); 

                auto tex_xxx_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 4); 

                auto tey_xxx_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 4); 

                auto tez_xxx_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 4); 

                auto tex_xxx_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 5); 

                auto tey_xxx_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 5); 

                auto tez_xxx_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 5); 

                auto tex_xxx_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 6); 

                auto tey_xxx_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 6); 

                auto tez_xxx_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 6); 

                auto tex_xxx_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 7); 

                auto tey_xxx_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 7); 

                auto tez_xxx_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 7); 

                auto tex_xxx_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 8); 

                auto tey_xxx_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 8); 

                auto tez_xxx_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 8); 

                auto tex_xxx_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 9); 

                auto tey_xxx_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 9); 

                auto tez_xxx_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 9); 

                auto tex_xxx_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 10); 

                auto tey_xxx_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 10); 

                auto tez_xxx_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 10); 

                auto tex_xxx_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 11); 

                auto tey_xxx_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 11); 

                auto tez_xxx_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 11); 

                auto tex_xxx_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 12); 

                auto tey_xxx_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 12); 

                auto tez_xxx_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 12); 

                auto tex_xxx_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 13); 

                auto tey_xxx_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 13); 

                auto tez_xxx_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 13); 

                auto tex_xxx_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 14); 

                auto tey_xxx_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 14); 

                auto tez_xxx_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 14); 

                auto tex_xxy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 15); 

                auto tey_xxy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 15); 

                auto tez_xxy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 15); 

                auto tex_xxy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 16); 

                auto tey_xxy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 16); 

                // Batch of Integrals (0,50)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xx_xxxx_1, ta_xx_xxxy_1, ta_xx_xxxz_1, ta_xx_xxyy_1, \
                                         ta_xx_xxyz_1, ta_xx_xxzz_1, ta_xx_xyyy_1, ta_xx_xyyz_1, ta_xx_xyzz_1, ta_xx_xzzz_1, \
                                         ta_xx_yyyy_1, ta_xx_yyyz_1, ta_xx_yyzz_1, ta_xx_yzzz_1, ta_xx_zzzz_1, ta_xy_xxxx_1, \
                                         ta_xy_xxxy_1, tex_x_xxxx_0, tex_x_xxxx_1, tex_x_xxxy_0, tex_x_xxxy_1, tex_x_xxxz_0, \
                                         tex_x_xxxz_1, tex_x_xxyy_0, tex_x_xxyy_1, tex_x_xxyz_0, tex_x_xxyz_1, tex_x_xxzz_0, \
                                         tex_x_xxzz_1, tex_x_xyyy_0, tex_x_xyyy_1, tex_x_xyyz_0, tex_x_xyyz_1, tex_x_xyzz_0, \
                                         tex_x_xyzz_1, tex_x_xzzz_0, tex_x_xzzz_1, tex_x_yyyy_0, tex_x_yyyy_1, tex_x_yyyz_0, \
                                         tex_x_yyyz_1, tex_x_yyzz_0, tex_x_yyzz_1, tex_x_yzzz_0, tex_x_yzzz_1, tex_x_zzzz_0, \
                                         tex_x_zzzz_1, tex_xx_xxx_0, tex_xx_xxx_1, tex_xx_xxxx_0, tex_xx_xxxx_1, \
                                         tex_xx_xxxy_0, tex_xx_xxxy_1, tex_xx_xxxz_0, tex_xx_xxxz_1, tex_xx_xxy_0, \
                                         tex_xx_xxy_1, tex_xx_xxyy_0, tex_xx_xxyy_1, tex_xx_xxyz_0, tex_xx_xxyz_1, \
                                         tex_xx_xxz_0, tex_xx_xxz_1, tex_xx_xxzz_0, tex_xx_xxzz_1, tex_xx_xyy_0, \
                                         tex_xx_xyy_1, tex_xx_xyyy_0, tex_xx_xyyy_1, tex_xx_xyyz_0, tex_xx_xyyz_1, \
                                         tex_xx_xyz_0, tex_xx_xyz_1, tex_xx_xyzz_0, tex_xx_xyzz_1, tex_xx_xzz_0, \
                                         tex_xx_xzz_1, tex_xx_xzzz_0, tex_xx_xzzz_1, tex_xx_yyy_0, tex_xx_yyy_1, \
                                         tex_xx_yyyy_0, tex_xx_yyyy_1, tex_xx_yyyz_0, tex_xx_yyyz_1, tex_xx_yyz_0, \
                                         tex_xx_yyz_1, tex_xx_yyzz_0, tex_xx_yyzz_1, tex_xx_yzz_0, tex_xx_yzz_1, \
                                         tex_xx_yzzz_0, tex_xx_yzzz_1, tex_xx_zzz_0, tex_xx_zzz_1, tex_xx_zzzz_0, \
                                         tex_xx_zzzz_1, tex_xxx_xxxx_0, tex_xxx_xxxy_0, tex_xxx_xxxz_0, tex_xxx_xxyy_0, \
                                         tex_xxx_xxyz_0, tex_xxx_xxzz_0, tex_xxx_xyyy_0, tex_xxx_xyyz_0, tex_xxx_xyzz_0, \
                                         tex_xxx_xzzz_0, tex_xxx_yyyy_0, tex_xxx_yyyz_0, tex_xxx_yyzz_0, tex_xxx_yzzz_0, \
                                         tex_xxx_zzzz_0, tex_xxy_xxxx_0, tex_xxy_xxxy_0, tex_xy_xxx_0, tex_xy_xxx_1, \
                                         tex_xy_xxxx_0, tex_xy_xxxx_1, tex_xy_xxxy_0, tex_xy_xxxy_1, tex_xy_xxy_0, \
                                         tex_xy_xxy_1, tex_y_xxxx_0, tex_y_xxxx_1, tex_y_xxxy_0, tex_y_xxxy_1, tey_x_xxxx_0, \
                                         tey_x_xxxx_1, tey_x_xxxy_0, tey_x_xxxy_1, tey_x_xxxz_0, tey_x_xxxz_1, tey_x_xxyy_0, \
                                         tey_x_xxyy_1, tey_x_xxyz_0, tey_x_xxyz_1, tey_x_xxzz_0, tey_x_xxzz_1, tey_x_xyyy_0, \
                                         tey_x_xyyy_1, tey_x_xyyz_0, tey_x_xyyz_1, tey_x_xyzz_0, tey_x_xyzz_1, tey_x_xzzz_0, \
                                         tey_x_xzzz_1, tey_x_yyyy_0, tey_x_yyyy_1, tey_x_yyyz_0, tey_x_yyyz_1, tey_x_yyzz_0, \
                                         tey_x_yyzz_1, tey_x_yzzz_0, tey_x_yzzz_1, tey_x_zzzz_0, tey_x_zzzz_1, tey_xx_xxx_0, \
                                         tey_xx_xxx_1, tey_xx_xxxx_0, tey_xx_xxxx_1, tey_xx_xxxy_0, tey_xx_xxxy_1, \
                                         tey_xx_xxxz_0, tey_xx_xxxz_1, tey_xx_xxy_0, tey_xx_xxy_1, tey_xx_xxyy_0, \
                                         tey_xx_xxyy_1, tey_xx_xxyz_0, tey_xx_xxyz_1, tey_xx_xxz_0, tey_xx_xxz_1, \
                                         tey_xx_xxzz_0, tey_xx_xxzz_1, tey_xx_xyy_0, tey_xx_xyy_1, tey_xx_xyyy_0, \
                                         tey_xx_xyyy_1, tey_xx_xyyz_0, tey_xx_xyyz_1, tey_xx_xyz_0, tey_xx_xyz_1, \
                                         tey_xx_xyzz_0, tey_xx_xyzz_1, tey_xx_xzz_0, tey_xx_xzz_1, tey_xx_xzzz_0, \
                                         tey_xx_xzzz_1, tey_xx_yyy_0, tey_xx_yyy_1, tey_xx_yyyy_0, tey_xx_yyyy_1, \
                                         tey_xx_yyyz_0, tey_xx_yyyz_1, tey_xx_yyz_0, tey_xx_yyz_1, tey_xx_yyzz_0, \
                                         tey_xx_yyzz_1, tey_xx_yzz_0, tey_xx_yzz_1, tey_xx_yzzz_0, tey_xx_yzzz_1, \
                                         tey_xx_zzz_0, tey_xx_zzz_1, tey_xx_zzzz_0, tey_xx_zzzz_1, tey_xxx_xxxx_0, \
                                         tey_xxx_xxxy_0, tey_xxx_xxxz_0, tey_xxx_xxyy_0, tey_xxx_xxyz_0, tey_xxx_xxzz_0, \
                                         tey_xxx_xyyy_0, tey_xxx_xyyz_0, tey_xxx_xyzz_0, tey_xxx_xzzz_0, tey_xxx_yyyy_0, \
                                         tey_xxx_yyyz_0, tey_xxx_yyzz_0, tey_xxx_yzzz_0, tey_xxx_zzzz_0, tey_xxy_xxxx_0, \
                                         tey_xxy_xxxy_0, tey_xy_xxx_0, tey_xy_xxx_1, tey_xy_xxxx_0, tey_xy_xxxx_1, \
                                         tey_xy_xxxy_0, tey_xy_xxxy_1, tey_xy_xxy_0, tey_xy_xxy_1, tey_y_xxxx_0, \
                                         tey_y_xxxx_1, tey_y_xxxy_0, tey_y_xxxy_1, tez_x_xxxx_0, tez_x_xxxx_1, tez_x_xxxy_0, \
                                         tez_x_xxxy_1, tez_x_xxxz_0, tez_x_xxxz_1, tez_x_xxyy_0, tez_x_xxyy_1, tez_x_xxyz_0, \
                                         tez_x_xxyz_1, tez_x_xxzz_0, tez_x_xxzz_1, tez_x_xyyy_0, tez_x_xyyy_1, tez_x_xyyz_0, \
                                         tez_x_xyyz_1, tez_x_xyzz_0, tez_x_xyzz_1, tez_x_xzzz_0, tez_x_xzzz_1, tez_x_yyyy_0, \
                                         tez_x_yyyy_1, tez_x_yyyz_0, tez_x_yyyz_1, tez_x_yyzz_0, tez_x_yyzz_1, tez_x_yzzz_0, \
                                         tez_x_yzzz_1, tez_x_zzzz_0, tez_x_zzzz_1, tez_xx_xxx_0, tez_xx_xxx_1, \
                                         tez_xx_xxxx_0, tez_xx_xxxx_1, tez_xx_xxxy_0, tez_xx_xxxy_1, tez_xx_xxxz_0, \
                                         tez_xx_xxxz_1, tez_xx_xxy_0, tez_xx_xxy_1, tez_xx_xxyy_0, tez_xx_xxyy_1, \
                                         tez_xx_xxyz_0, tez_xx_xxyz_1, tez_xx_xxz_0, tez_xx_xxz_1, tez_xx_xxzz_0, \
                                         tez_xx_xxzz_1, tez_xx_xyy_0, tez_xx_xyy_1, tez_xx_xyyy_0, tez_xx_xyyy_1, \
                                         tez_xx_xyyz_0, tez_xx_xyyz_1, tez_xx_xyz_0, tez_xx_xyz_1, tez_xx_xyzz_0, \
                                         tez_xx_xyzz_1, tez_xx_xzz_0, tez_xx_xzz_1, tez_xx_xzzz_0, tez_xx_xzzz_1, \
                                         tez_xx_yyy_0, tez_xx_yyy_1, tez_xx_yyyy_0, tez_xx_yyyy_1, tez_xx_yyyz_0, \
                                         tez_xx_yyyz_1, tez_xx_yyz_0, tez_xx_yyz_1, tez_xx_yyzz_0, tez_xx_yyzz_1, \
                                         tez_xx_yzz_0, tez_xx_yzz_1, tez_xx_yzzz_0, tez_xx_yzzz_1, tez_xx_zzz_0, \
                                         tez_xx_zzz_1, tez_xx_zzzz_0, tez_xx_zzzz_1, tez_xxx_xxxx_0, tez_xxx_xxxy_0, \
                                         tez_xxx_xxxz_0, tez_xxx_xxyy_0, tez_xxx_xxyz_0, tez_xxx_xxzz_0, tez_xxx_xyyy_0, \
                                         tez_xxx_xyyz_0, tez_xxx_xyzz_0, tez_xxx_xzzz_0, tez_xxx_yyyy_0, tez_xxx_yyyz_0, \
                                         tez_xxx_yyzz_0, tez_xxx_yzzz_0, tez_xxx_zzzz_0, tez_xxy_xxxx_0, tez_xy_xxx_0, \
                                         tez_xy_xxx_1, tez_xy_xxxx_0, tez_xy_xxxx_1, tez_y_xxxx_0, tez_y_xxxx_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xxx_xxxx_0[j] = pa_x[j] * tex_xx_xxxx_0[j] - pc_x[j] * tex_xx_xxxx_1[j] + fl1_fx * tex_x_xxxx_0[j] - fl1_fx * tex_x_xxxx_1[j] + 2.0 * fl1_fx * tex_xx_xxx_0[j] - 2.0 * fl1_fx * tex_xx_xxx_1[j] + ta_xx_xxxx_1[j];

                    tey_xxx_xxxx_0[j] = pa_x[j] * tey_xx_xxxx_0[j] - pc_x[j] * tey_xx_xxxx_1[j] + fl1_fx * tey_x_xxxx_0[j] - fl1_fx * tey_x_xxxx_1[j] + 2.0 * fl1_fx * tey_xx_xxx_0[j] - 2.0 * fl1_fx * tey_xx_xxx_1[j];

                    tez_xxx_xxxx_0[j] = pa_x[j] * tez_xx_xxxx_0[j] - pc_x[j] * tez_xx_xxxx_1[j] + fl1_fx * tez_x_xxxx_0[j] - fl1_fx * tez_x_xxxx_1[j] + 2.0 * fl1_fx * tez_xx_xxx_0[j] - 2.0 * fl1_fx * tez_xx_xxx_1[j];

                    tex_xxx_xxxy_0[j] = pa_x[j] * tex_xx_xxxy_0[j] - pc_x[j] * tex_xx_xxxy_1[j] + fl1_fx * tex_x_xxxy_0[j] - fl1_fx * tex_x_xxxy_1[j] + 1.5 * fl1_fx * tex_xx_xxy_0[j] - 1.5 * fl1_fx * tex_xx_xxy_1[j] + ta_xx_xxxy_1[j];

                    tey_xxx_xxxy_0[j] = pa_x[j] * tey_xx_xxxy_0[j] - pc_x[j] * tey_xx_xxxy_1[j] + fl1_fx * tey_x_xxxy_0[j] - fl1_fx * tey_x_xxxy_1[j] + 1.5 * fl1_fx * tey_xx_xxy_0[j] - 1.5 * fl1_fx * tey_xx_xxy_1[j];

                    tez_xxx_xxxy_0[j] = pa_x[j] * tez_xx_xxxy_0[j] - pc_x[j] * tez_xx_xxxy_1[j] + fl1_fx * tez_x_xxxy_0[j] - fl1_fx * tez_x_xxxy_1[j] + 1.5 * fl1_fx * tez_xx_xxy_0[j] - 1.5 * fl1_fx * tez_xx_xxy_1[j];

                    tex_xxx_xxxz_0[j] = pa_x[j] * tex_xx_xxxz_0[j] - pc_x[j] * tex_xx_xxxz_1[j] + fl1_fx * tex_x_xxxz_0[j] - fl1_fx * tex_x_xxxz_1[j] + 1.5 * fl1_fx * tex_xx_xxz_0[j] - 1.5 * fl1_fx * tex_xx_xxz_1[j] + ta_xx_xxxz_1[j];

                    tey_xxx_xxxz_0[j] = pa_x[j] * tey_xx_xxxz_0[j] - pc_x[j] * tey_xx_xxxz_1[j] + fl1_fx * tey_x_xxxz_0[j] - fl1_fx * tey_x_xxxz_1[j] + 1.5 * fl1_fx * tey_xx_xxz_0[j] - 1.5 * fl1_fx * tey_xx_xxz_1[j];

                    tez_xxx_xxxz_0[j] = pa_x[j] * tez_xx_xxxz_0[j] - pc_x[j] * tez_xx_xxxz_1[j] + fl1_fx * tez_x_xxxz_0[j] - fl1_fx * tez_x_xxxz_1[j] + 1.5 * fl1_fx * tez_xx_xxz_0[j] - 1.5 * fl1_fx * tez_xx_xxz_1[j];

                    tex_xxx_xxyy_0[j] = pa_x[j] * tex_xx_xxyy_0[j] - pc_x[j] * tex_xx_xxyy_1[j] + fl1_fx * tex_x_xxyy_0[j] - fl1_fx * tex_x_xxyy_1[j] + fl1_fx * tex_xx_xyy_0[j] - fl1_fx * tex_xx_xyy_1[j] + ta_xx_xxyy_1[j];

                    tey_xxx_xxyy_0[j] = pa_x[j] * tey_xx_xxyy_0[j] - pc_x[j] * tey_xx_xxyy_1[j] + fl1_fx * tey_x_xxyy_0[j] - fl1_fx * tey_x_xxyy_1[j] + fl1_fx * tey_xx_xyy_0[j] - fl1_fx * tey_xx_xyy_1[j];

                    tez_xxx_xxyy_0[j] = pa_x[j] * tez_xx_xxyy_0[j] - pc_x[j] * tez_xx_xxyy_1[j] + fl1_fx * tez_x_xxyy_0[j] - fl1_fx * tez_x_xxyy_1[j] + fl1_fx * tez_xx_xyy_0[j] - fl1_fx * tez_xx_xyy_1[j];

                    tex_xxx_xxyz_0[j] = pa_x[j] * tex_xx_xxyz_0[j] - pc_x[j] * tex_xx_xxyz_1[j] + fl1_fx * tex_x_xxyz_0[j] - fl1_fx * tex_x_xxyz_1[j] + fl1_fx * tex_xx_xyz_0[j] - fl1_fx * tex_xx_xyz_1[j] + ta_xx_xxyz_1[j];

                    tey_xxx_xxyz_0[j] = pa_x[j] * tey_xx_xxyz_0[j] - pc_x[j] * tey_xx_xxyz_1[j] + fl1_fx * tey_x_xxyz_0[j] - fl1_fx * tey_x_xxyz_1[j] + fl1_fx * tey_xx_xyz_0[j] - fl1_fx * tey_xx_xyz_1[j];

                    tez_xxx_xxyz_0[j] = pa_x[j] * tez_xx_xxyz_0[j] - pc_x[j] * tez_xx_xxyz_1[j] + fl1_fx * tez_x_xxyz_0[j] - fl1_fx * tez_x_xxyz_1[j] + fl1_fx * tez_xx_xyz_0[j] - fl1_fx * tez_xx_xyz_1[j];

                    tex_xxx_xxzz_0[j] = pa_x[j] * tex_xx_xxzz_0[j] - pc_x[j] * tex_xx_xxzz_1[j] + fl1_fx * tex_x_xxzz_0[j] - fl1_fx * tex_x_xxzz_1[j] + fl1_fx * tex_xx_xzz_0[j] - fl1_fx * tex_xx_xzz_1[j] + ta_xx_xxzz_1[j];

                    tey_xxx_xxzz_0[j] = pa_x[j] * tey_xx_xxzz_0[j] - pc_x[j] * tey_xx_xxzz_1[j] + fl1_fx * tey_x_xxzz_0[j] - fl1_fx * tey_x_xxzz_1[j] + fl1_fx * tey_xx_xzz_0[j] - fl1_fx * tey_xx_xzz_1[j];

                    tez_xxx_xxzz_0[j] = pa_x[j] * tez_xx_xxzz_0[j] - pc_x[j] * tez_xx_xxzz_1[j] + fl1_fx * tez_x_xxzz_0[j] - fl1_fx * tez_x_xxzz_1[j] + fl1_fx * tez_xx_xzz_0[j] - fl1_fx * tez_xx_xzz_1[j];

                    tex_xxx_xyyy_0[j] = pa_x[j] * tex_xx_xyyy_0[j] - pc_x[j] * tex_xx_xyyy_1[j] + fl1_fx * tex_x_xyyy_0[j] - fl1_fx * tex_x_xyyy_1[j] + 0.5 * fl1_fx * tex_xx_yyy_0[j] - 0.5 * fl1_fx * tex_xx_yyy_1[j] + ta_xx_xyyy_1[j];

                    tey_xxx_xyyy_0[j] = pa_x[j] * tey_xx_xyyy_0[j] - pc_x[j] * tey_xx_xyyy_1[j] + fl1_fx * tey_x_xyyy_0[j] - fl1_fx * tey_x_xyyy_1[j] + 0.5 * fl1_fx * tey_xx_yyy_0[j] - 0.5 * fl1_fx * tey_xx_yyy_1[j];

                    tez_xxx_xyyy_0[j] = pa_x[j] * tez_xx_xyyy_0[j] - pc_x[j] * tez_xx_xyyy_1[j] + fl1_fx * tez_x_xyyy_0[j] - fl1_fx * tez_x_xyyy_1[j] + 0.5 * fl1_fx * tez_xx_yyy_0[j] - 0.5 * fl1_fx * tez_xx_yyy_1[j];

                    tex_xxx_xyyz_0[j] = pa_x[j] * tex_xx_xyyz_0[j] - pc_x[j] * tex_xx_xyyz_1[j] + fl1_fx * tex_x_xyyz_0[j] - fl1_fx * tex_x_xyyz_1[j] + 0.5 * fl1_fx * tex_xx_yyz_0[j] - 0.5 * fl1_fx * tex_xx_yyz_1[j] + ta_xx_xyyz_1[j];

                    tey_xxx_xyyz_0[j] = pa_x[j] * tey_xx_xyyz_0[j] - pc_x[j] * tey_xx_xyyz_1[j] + fl1_fx * tey_x_xyyz_0[j] - fl1_fx * tey_x_xyyz_1[j] + 0.5 * fl1_fx * tey_xx_yyz_0[j] - 0.5 * fl1_fx * tey_xx_yyz_1[j];

                    tez_xxx_xyyz_0[j] = pa_x[j] * tez_xx_xyyz_0[j] - pc_x[j] * tez_xx_xyyz_1[j] + fl1_fx * tez_x_xyyz_0[j] - fl1_fx * tez_x_xyyz_1[j] + 0.5 * fl1_fx * tez_xx_yyz_0[j] - 0.5 * fl1_fx * tez_xx_yyz_1[j];

                    tex_xxx_xyzz_0[j] = pa_x[j] * tex_xx_xyzz_0[j] - pc_x[j] * tex_xx_xyzz_1[j] + fl1_fx * tex_x_xyzz_0[j] - fl1_fx * tex_x_xyzz_1[j] + 0.5 * fl1_fx * tex_xx_yzz_0[j] - 0.5 * fl1_fx * tex_xx_yzz_1[j] + ta_xx_xyzz_1[j];

                    tey_xxx_xyzz_0[j] = pa_x[j] * tey_xx_xyzz_0[j] - pc_x[j] * tey_xx_xyzz_1[j] + fl1_fx * tey_x_xyzz_0[j] - fl1_fx * tey_x_xyzz_1[j] + 0.5 * fl1_fx * tey_xx_yzz_0[j] - 0.5 * fl1_fx * tey_xx_yzz_1[j];

                    tez_xxx_xyzz_0[j] = pa_x[j] * tez_xx_xyzz_0[j] - pc_x[j] * tez_xx_xyzz_1[j] + fl1_fx * tez_x_xyzz_0[j] - fl1_fx * tez_x_xyzz_1[j] + 0.5 * fl1_fx * tez_xx_yzz_0[j] - 0.5 * fl1_fx * tez_xx_yzz_1[j];

                    tex_xxx_xzzz_0[j] = pa_x[j] * tex_xx_xzzz_0[j] - pc_x[j] * tex_xx_xzzz_1[j] + fl1_fx * tex_x_xzzz_0[j] - fl1_fx * tex_x_xzzz_1[j] + 0.5 * fl1_fx * tex_xx_zzz_0[j] - 0.5 * fl1_fx * tex_xx_zzz_1[j] + ta_xx_xzzz_1[j];

                    tey_xxx_xzzz_0[j] = pa_x[j] * tey_xx_xzzz_0[j] - pc_x[j] * tey_xx_xzzz_1[j] + fl1_fx * tey_x_xzzz_0[j] - fl1_fx * tey_x_xzzz_1[j] + 0.5 * fl1_fx * tey_xx_zzz_0[j] - 0.5 * fl1_fx * tey_xx_zzz_1[j];

                    tez_xxx_xzzz_0[j] = pa_x[j] * tez_xx_xzzz_0[j] - pc_x[j] * tez_xx_xzzz_1[j] + fl1_fx * tez_x_xzzz_0[j] - fl1_fx * tez_x_xzzz_1[j] + 0.5 * fl1_fx * tez_xx_zzz_0[j] - 0.5 * fl1_fx * tez_xx_zzz_1[j];

                    tex_xxx_yyyy_0[j] = pa_x[j] * tex_xx_yyyy_0[j] - pc_x[j] * tex_xx_yyyy_1[j] + fl1_fx * tex_x_yyyy_0[j] - fl1_fx * tex_x_yyyy_1[j] + ta_xx_yyyy_1[j];

                    tey_xxx_yyyy_0[j] = pa_x[j] * tey_xx_yyyy_0[j] - pc_x[j] * tey_xx_yyyy_1[j] + fl1_fx * tey_x_yyyy_0[j] - fl1_fx * tey_x_yyyy_1[j];

                    tez_xxx_yyyy_0[j] = pa_x[j] * tez_xx_yyyy_0[j] - pc_x[j] * tez_xx_yyyy_1[j] + fl1_fx * tez_x_yyyy_0[j] - fl1_fx * tez_x_yyyy_1[j];

                    tex_xxx_yyyz_0[j] = pa_x[j] * tex_xx_yyyz_0[j] - pc_x[j] * tex_xx_yyyz_1[j] + fl1_fx * tex_x_yyyz_0[j] - fl1_fx * tex_x_yyyz_1[j] + ta_xx_yyyz_1[j];

                    tey_xxx_yyyz_0[j] = pa_x[j] * tey_xx_yyyz_0[j] - pc_x[j] * tey_xx_yyyz_1[j] + fl1_fx * tey_x_yyyz_0[j] - fl1_fx * tey_x_yyyz_1[j];

                    tez_xxx_yyyz_0[j] = pa_x[j] * tez_xx_yyyz_0[j] - pc_x[j] * tez_xx_yyyz_1[j] + fl1_fx * tez_x_yyyz_0[j] - fl1_fx * tez_x_yyyz_1[j];

                    tex_xxx_yyzz_0[j] = pa_x[j] * tex_xx_yyzz_0[j] - pc_x[j] * tex_xx_yyzz_1[j] + fl1_fx * tex_x_yyzz_0[j] - fl1_fx * tex_x_yyzz_1[j] + ta_xx_yyzz_1[j];

                    tey_xxx_yyzz_0[j] = pa_x[j] * tey_xx_yyzz_0[j] - pc_x[j] * tey_xx_yyzz_1[j] + fl1_fx * tey_x_yyzz_0[j] - fl1_fx * tey_x_yyzz_1[j];

                    tez_xxx_yyzz_0[j] = pa_x[j] * tez_xx_yyzz_0[j] - pc_x[j] * tez_xx_yyzz_1[j] + fl1_fx * tez_x_yyzz_0[j] - fl1_fx * tez_x_yyzz_1[j];

                    tex_xxx_yzzz_0[j] = pa_x[j] * tex_xx_yzzz_0[j] - pc_x[j] * tex_xx_yzzz_1[j] + fl1_fx * tex_x_yzzz_0[j] - fl1_fx * tex_x_yzzz_1[j] + ta_xx_yzzz_1[j];

                    tey_xxx_yzzz_0[j] = pa_x[j] * tey_xx_yzzz_0[j] - pc_x[j] * tey_xx_yzzz_1[j] + fl1_fx * tey_x_yzzz_0[j] - fl1_fx * tey_x_yzzz_1[j];

                    tez_xxx_yzzz_0[j] = pa_x[j] * tez_xx_yzzz_0[j] - pc_x[j] * tez_xx_yzzz_1[j] + fl1_fx * tez_x_yzzz_0[j] - fl1_fx * tez_x_yzzz_1[j];

                    tex_xxx_zzzz_0[j] = pa_x[j] * tex_xx_zzzz_0[j] - pc_x[j] * tex_xx_zzzz_1[j] + fl1_fx * tex_x_zzzz_0[j] - fl1_fx * tex_x_zzzz_1[j] + ta_xx_zzzz_1[j];

                    tey_xxx_zzzz_0[j] = pa_x[j] * tey_xx_zzzz_0[j] - pc_x[j] * tey_xx_zzzz_1[j] + fl1_fx * tey_x_zzzz_0[j] - fl1_fx * tey_x_zzzz_1[j];

                    tez_xxx_zzzz_0[j] = pa_x[j] * tez_xx_zzzz_0[j] - pc_x[j] * tez_xx_zzzz_1[j] + fl1_fx * tez_x_zzzz_0[j] - fl1_fx * tez_x_zzzz_1[j];

                    tex_xxy_xxxx_0[j] = pa_x[j] * tex_xy_xxxx_0[j] - pc_x[j] * tex_xy_xxxx_1[j] + 0.5 * fl1_fx * tex_y_xxxx_0[j] - 0.5 * fl1_fx * tex_y_xxxx_1[j] + 2.0 * fl1_fx * tex_xy_xxx_0[j] - 2.0 * fl1_fx * tex_xy_xxx_1[j] + ta_xy_xxxx_1[j];

                    tey_xxy_xxxx_0[j] = pa_x[j] * tey_xy_xxxx_0[j] - pc_x[j] * tey_xy_xxxx_1[j] + 0.5 * fl1_fx * tey_y_xxxx_0[j] - 0.5 * fl1_fx * tey_y_xxxx_1[j] + 2.0 * fl1_fx * tey_xy_xxx_0[j] - 2.0 * fl1_fx * tey_xy_xxx_1[j];

                    tez_xxy_xxxx_0[j] = pa_x[j] * tez_xy_xxxx_0[j] - pc_x[j] * tez_xy_xxxx_1[j] + 0.5 * fl1_fx * tez_y_xxxx_0[j] - 0.5 * fl1_fx * tez_y_xxxx_1[j] + 2.0 * fl1_fx * tez_xy_xxx_0[j] - 2.0 * fl1_fx * tez_xy_xxx_1[j];

                    tex_xxy_xxxy_0[j] = pa_x[j] * tex_xy_xxxy_0[j] - pc_x[j] * tex_xy_xxxy_1[j] + 0.5 * fl1_fx * tex_y_xxxy_0[j] - 0.5 * fl1_fx * tex_y_xxxy_1[j] + 1.5 * fl1_fx * tex_xy_xxy_0[j] - 1.5 * fl1_fx * tex_xy_xxy_1[j] + ta_xy_xxxy_1[j];

                    tey_xxy_xxxy_0[j] = pa_x[j] * tey_xy_xxxy_0[j] - pc_x[j] * tey_xy_xxxy_1[j] + 0.5 * fl1_fx * tey_y_xxxy_0[j] - 0.5 * fl1_fx * tey_y_xxxy_1[j] + 1.5 * fl1_fx * tey_xy_xxy_0[j] - 1.5 * fl1_fx * tey_xy_xxy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_50_100(      CMemBlock2D<double>& primBuffer,
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

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tez_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 16); 

                auto tex_xy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 17); 

                auto tey_xy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 17); 

                auto tez_xy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 17); 

                auto tex_xy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 18); 

                auto tey_xy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 18); 

                auto tez_xy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 18); 

                auto tex_xy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 19); 

                auto tey_xy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 19); 

                auto tez_xy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 19); 

                auto tex_xy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 20); 

                auto tey_xy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 20); 

                auto tez_xy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 20); 

                auto tex_xy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 21); 

                auto tey_xy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 21); 

                auto tez_xy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 21); 

                auto tex_xy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 22); 

                auto tey_xy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 22); 

                auto tez_xy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 22); 

                auto tex_xy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 23); 

                auto tey_xy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 23); 

                auto tez_xy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 23); 

                auto tex_xy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 24); 

                auto tey_xy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 24); 

                auto tez_xy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 24); 

                auto tex_xy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 25); 

                auto tey_xy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 25); 

                auto tez_xy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 25); 

                auto tex_xy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 26); 

                auto tey_xy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 26); 

                auto tez_xy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 26); 

                auto tex_xy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 27); 

                auto tey_xy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 27); 

                auto tez_xy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 27); 

                auto tex_xy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 28); 

                auto tey_xy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 28); 

                auto tez_xy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 28); 

                auto tex_xy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 29); 

                auto tey_xy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 29); 

                auto tez_xy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 29); 

                auto tex_xz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 30); 

                auto tey_xz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 30); 

                auto tez_xz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 30); 

                auto tex_xz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 31); 

                auto tey_xz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 31); 

                auto tez_xz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 31); 

                auto tex_xz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 32); 

                auto tey_xz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 32); 

                auto tez_xz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 32); 

                auto tex_xz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 33); 

                auto tez_xy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 16); 

                auto tex_xy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 17); 

                auto tey_xy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 17); 

                auto tez_xy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 17); 

                auto tex_xy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 18); 

                auto tey_xy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 18); 

                auto tez_xy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 18); 

                auto tex_xy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 19); 

                auto tey_xy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 19); 

                auto tez_xy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 19); 

                auto tex_xy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 20); 

                auto tey_xy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 20); 

                auto tez_xy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 20); 

                auto tex_xy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 21); 

                auto tey_xy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 21); 

                auto tez_xy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 21); 

                auto tex_xy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 22); 

                auto tey_xy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 22); 

                auto tez_xy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 22); 

                auto tex_xy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 23); 

                auto tey_xy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 23); 

                auto tez_xy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 23); 

                auto tex_xy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 24); 

                auto tey_xy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 24); 

                auto tez_xy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 24); 

                auto tex_xy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 25); 

                auto tey_xy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 25); 

                auto tez_xy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 25); 

                auto tex_xy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 26); 

                auto tey_xy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 26); 

                auto tez_xy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 26); 

                auto tex_xy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 27); 

                auto tey_xy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 27); 

                auto tez_xy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 27); 

                auto tex_xy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 28); 

                auto tey_xy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 28); 

                auto tez_xy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 28); 

                auto tex_xy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 29); 

                auto tey_xy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 29); 

                auto tez_xy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 29); 

                auto tex_xz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 30); 

                auto tey_xz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 30); 

                auto tez_xz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 30); 

                auto tex_xz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 31); 

                auto tey_xz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 31); 

                auto tez_xz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 31); 

                auto tex_xz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 32); 

                auto tey_xz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 32); 

                auto tez_xz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 32); 

                auto tex_xz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 33); 

                auto tez_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 16); 

                auto tex_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 17); 

                auto tey_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 17); 

                auto tez_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 17); 

                auto tex_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 18); 

                auto tey_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 18); 

                auto tez_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 18); 

                auto tex_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 19); 

                auto tey_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 19); 

                auto tez_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 19); 

                auto tex_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 20); 

                auto tey_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 20); 

                auto tez_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 20); 

                auto tex_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 21); 

                auto tey_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 21); 

                auto tez_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 21); 

                auto tex_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 22); 

                auto tey_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 22); 

                auto tez_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 22); 

                auto tex_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 23); 

                auto tey_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 23); 

                auto tez_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 23); 

                auto tex_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 24); 

                auto tey_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 24); 

                auto tez_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 24); 

                auto tex_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 25); 

                auto tey_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 25); 

                auto tez_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 25); 

                auto tex_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 26); 

                auto tey_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 26); 

                auto tez_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 26); 

                auto tex_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 27); 

                auto tey_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 27); 

                auto tez_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 27); 

                auto tex_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 28); 

                auto tey_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 28); 

                auto tez_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 28); 

                auto tex_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 29); 

                auto tey_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 29); 

                auto tez_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 29); 

                auto tex_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 30); 

                auto tey_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 30); 

                auto tez_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 30); 

                auto tex_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 31); 

                auto tey_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 31); 

                auto tez_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 31); 

                auto tex_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 32); 

                auto tey_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 32); 

                auto tez_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 32); 

                auto tex_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 33); 

                auto tez_y_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 16); 

                auto tex_y_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 17); 

                auto tey_y_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 17); 

                auto tez_y_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 17); 

                auto tex_y_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 18); 

                auto tey_y_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 18); 

                auto tez_y_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 18); 

                auto tex_y_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 19); 

                auto tey_y_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 19); 

                auto tez_y_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 19); 

                auto tex_y_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 20); 

                auto tey_y_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 20); 

                auto tez_y_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 20); 

                auto tex_y_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 21); 

                auto tey_y_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 21); 

                auto tez_y_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 21); 

                auto tex_y_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 22); 

                auto tey_y_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 22); 

                auto tez_y_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 22); 

                auto tex_y_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 23); 

                auto tey_y_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 23); 

                auto tez_y_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 23); 

                auto tex_y_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 24); 

                auto tey_y_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 24); 

                auto tez_y_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 24); 

                auto tex_y_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 25); 

                auto tey_y_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 25); 

                auto tez_y_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 25); 

                auto tex_y_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 26); 

                auto tey_y_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 26); 

                auto tez_y_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 26); 

                auto tex_y_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 27); 

                auto tey_y_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 27); 

                auto tez_y_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 27); 

                auto tex_y_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 28); 

                auto tey_y_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 28); 

                auto tez_y_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 28); 

                auto tex_y_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 29); 

                auto tey_y_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 29); 

                auto tez_y_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 29); 

                auto tex_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 30); 

                auto tey_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 30); 

                auto tez_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 30); 

                auto tex_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 31); 

                auto tey_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 31); 

                auto tez_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 31); 

                auto tex_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 32); 

                auto tey_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 32); 

                auto tez_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 32); 

                auto tex_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 33); 

                auto tez_xy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 11); 

                auto tex_xy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 12); 

                auto tey_xy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 12); 

                auto tez_xy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 12); 

                auto tex_xy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 13); 

                auto tey_xy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 13); 

                auto tez_xy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 13); 

                auto tex_xy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 14); 

                auto tey_xy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 14); 

                auto tez_xy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 14); 

                auto tex_xy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 15); 

                auto tey_xy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 15); 

                auto tez_xy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 15); 

                auto tex_xy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 16); 

                auto tey_xy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 16); 

                auto tez_xy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 16); 

                auto tex_xy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 17); 

                auto tey_xy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 17); 

                auto tez_xy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 17); 

                auto tex_xy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 18); 

                auto tey_xy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 18); 

                auto tez_xy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 18); 

                auto tex_xy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 19); 

                auto tey_xy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 19); 

                auto tez_xy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 19); 

                auto tex_xz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 20); 

                auto tey_xz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 20); 

                auto tez_xz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 20); 

                auto tex_xz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 21); 

                auto tey_xz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 21); 

                auto tez_xz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 21); 

                auto tex_xz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 22); 

                auto tey_xz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 22); 

                auto tez_xz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 22); 

                auto tex_xz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 23); 

                auto tez_xy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 11); 

                auto tex_xy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 12); 

                auto tey_xy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 12); 

                auto tez_xy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 12); 

                auto tex_xy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 13); 

                auto tey_xy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 13); 

                auto tez_xy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 13); 

                auto tex_xy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 14); 

                auto tey_xy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 14); 

                auto tez_xy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 14); 

                auto tex_xy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 15); 

                auto tey_xy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 15); 

                auto tez_xy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 15); 

                auto tex_xy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 16); 

                auto tey_xy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 16); 

                auto tez_xy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 16); 

                auto tex_xy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 17); 

                auto tey_xy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 17); 

                auto tez_xy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 17); 

                auto tex_xy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 18); 

                auto tey_xy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 18); 

                auto tez_xy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 18); 

                auto tex_xy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 19); 

                auto tey_xy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 19); 

                auto tez_xy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 19); 

                auto tex_xz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 20); 

                auto tey_xz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 20); 

                auto tez_xz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 20); 

                auto tex_xz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 21); 

                auto tey_xz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 21); 

                auto tez_xz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 21); 

                auto tex_xz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 22); 

                auto tey_xz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 22); 

                auto tez_xz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 22); 

                auto tex_xz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 23); 

                auto ta_xy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 17); 

                auto ta_xy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 18); 

                auto ta_xy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 19); 

                auto ta_xy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 20); 

                auto ta_xy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 21); 

                auto ta_xy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 22); 

                auto ta_xy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 23); 

                auto ta_xy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 24); 

                auto ta_xy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 25); 

                auto ta_xy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 26); 

                auto ta_xy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 27); 

                auto ta_xy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 28); 

                auto ta_xy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 29); 

                auto ta_xz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 30); 

                auto ta_xz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 31); 

                auto ta_xz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 32); 

                auto ta_xz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 33); 

                // set up pointers to integrals

                auto tez_xxy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 16); 

                auto tex_xxy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 17); 

                auto tey_xxy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 17); 

                auto tez_xxy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 17); 

                auto tex_xxy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 18); 

                auto tey_xxy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 18); 

                auto tez_xxy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 18); 

                auto tex_xxy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 19); 

                auto tey_xxy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 19); 

                auto tez_xxy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 19); 

                auto tex_xxy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 20); 

                auto tey_xxy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 20); 

                auto tez_xxy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 20); 

                auto tex_xxy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 21); 

                auto tey_xxy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 21); 

                auto tez_xxy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 21); 

                auto tex_xxy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 22); 

                auto tey_xxy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 22); 

                auto tez_xxy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 22); 

                auto tex_xxy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 23); 

                auto tey_xxy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 23); 

                auto tez_xxy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 23); 

                auto tex_xxy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 24); 

                auto tey_xxy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 24); 

                auto tez_xxy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 24); 

                auto tex_xxy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 25); 

                auto tey_xxy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 25); 

                auto tez_xxy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 25); 

                auto tex_xxy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 26); 

                auto tey_xxy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 26); 

                auto tez_xxy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 26); 

                auto tex_xxy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 27); 

                auto tey_xxy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 27); 

                auto tez_xxy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 27); 

                auto tex_xxy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 28); 

                auto tey_xxy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 28); 

                auto tez_xxy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 28); 

                auto tex_xxy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 29); 

                auto tey_xxy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 29); 

                auto tez_xxy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 29); 

                auto tex_xxz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 30); 

                auto tey_xxz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 30); 

                auto tez_xxz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 30); 

                auto tex_xxz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 31); 

                auto tey_xxz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 31); 

                auto tez_xxz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 31); 

                auto tex_xxz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 32); 

                auto tey_xxz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 32); 

                auto tez_xxz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 32); 

                auto tex_xxz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 33); 

                // Batch of Integrals (50,100)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xy_xxxz_1, ta_xy_xxyy_1, ta_xy_xxyz_1, ta_xy_xxzz_1, \
                                         ta_xy_xyyy_1, ta_xy_xyyz_1, ta_xy_xyzz_1, ta_xy_xzzz_1, ta_xy_yyyy_1, ta_xy_yyyz_1, \
                                         ta_xy_yyzz_1, ta_xy_yzzz_1, ta_xy_zzzz_1, ta_xz_xxxx_1, ta_xz_xxxy_1, ta_xz_xxxz_1, \
                                         ta_xz_xxyy_1, tex_xxy_xxxz_0, tex_xxy_xxyy_0, tex_xxy_xxyz_0, tex_xxy_xxzz_0, \
                                         tex_xxy_xyyy_0, tex_xxy_xyyz_0, tex_xxy_xyzz_0, tex_xxy_xzzz_0, tex_xxy_yyyy_0, \
                                         tex_xxy_yyyz_0, tex_xxy_yyzz_0, tex_xxy_yzzz_0, tex_xxy_zzzz_0, tex_xxz_xxxx_0, \
                                         tex_xxz_xxxy_0, tex_xxz_xxxz_0, tex_xxz_xxyy_0, tex_xy_xxxz_0, tex_xy_xxxz_1, \
                                         tex_xy_xxyy_0, tex_xy_xxyy_1, tex_xy_xxyz_0, tex_xy_xxyz_1, tex_xy_xxz_0, \
                                         tex_xy_xxz_1, tex_xy_xxzz_0, tex_xy_xxzz_1, tex_xy_xyy_0, tex_xy_xyy_1, \
                                         tex_xy_xyyy_0, tex_xy_xyyy_1, tex_xy_xyyz_0, tex_xy_xyyz_1, tex_xy_xyz_0, \
                                         tex_xy_xyz_1, tex_xy_xyzz_0, tex_xy_xyzz_1, tex_xy_xzz_0, tex_xy_xzz_1, \
                                         tex_xy_xzzz_0, tex_xy_xzzz_1, tex_xy_yyy_0, tex_xy_yyy_1, tex_xy_yyyy_0, \
                                         tex_xy_yyyy_1, tex_xy_yyyz_0, tex_xy_yyyz_1, tex_xy_yyz_0, tex_xy_yyz_1, \
                                         tex_xy_yyzz_0, tex_xy_yyzz_1, tex_xy_yzz_0, tex_xy_yzz_1, tex_xy_yzzz_0, \
                                         tex_xy_yzzz_1, tex_xy_zzz_0, tex_xy_zzz_1, tex_xy_zzzz_0, tex_xy_zzzz_1, \
                                         tex_xz_xxx_0, tex_xz_xxx_1, tex_xz_xxxx_0, tex_xz_xxxx_1, tex_xz_xxxy_0, \
                                         tex_xz_xxxy_1, tex_xz_xxxz_0, tex_xz_xxxz_1, tex_xz_xxy_0, tex_xz_xxy_1, \
                                         tex_xz_xxyy_0, tex_xz_xxyy_1, tex_xz_xxz_0, tex_xz_xxz_1, tex_xz_xyy_0, \
                                         tex_xz_xyy_1, tex_y_xxxz_0, tex_y_xxxz_1, tex_y_xxyy_0, tex_y_xxyy_1, tex_y_xxyz_0, \
                                         tex_y_xxyz_1, tex_y_xxzz_0, tex_y_xxzz_1, tex_y_xyyy_0, tex_y_xyyy_1, tex_y_xyyz_0, \
                                         tex_y_xyyz_1, tex_y_xyzz_0, tex_y_xyzz_1, tex_y_xzzz_0, tex_y_xzzz_1, tex_y_yyyy_0, \
                                         tex_y_yyyy_1, tex_y_yyyz_0, tex_y_yyyz_1, tex_y_yyzz_0, tex_y_yyzz_1, tex_y_yzzz_0, \
                                         tex_y_yzzz_1, tex_y_zzzz_0, tex_y_zzzz_1, tex_z_xxxx_0, tex_z_xxxx_1, tex_z_xxxy_0, \
                                         tex_z_xxxy_1, tex_z_xxxz_0, tex_z_xxxz_1, tex_z_xxyy_0, tex_z_xxyy_1, \
                                         tey_xxy_xxxz_0, tey_xxy_xxyy_0, tey_xxy_xxyz_0, tey_xxy_xxzz_0, tey_xxy_xyyy_0, \
                                         tey_xxy_xyyz_0, tey_xxy_xyzz_0, tey_xxy_xzzz_0, tey_xxy_yyyy_0, tey_xxy_yyyz_0, \
                                         tey_xxy_yyzz_0, tey_xxy_yzzz_0, tey_xxy_zzzz_0, tey_xxz_xxxx_0, tey_xxz_xxxy_0, \
                                         tey_xxz_xxxz_0, tey_xy_xxxz_0, tey_xy_xxxz_1, tey_xy_xxyy_0, tey_xy_xxyy_1, \
                                         tey_xy_xxyz_0, tey_xy_xxyz_1, tey_xy_xxz_0, tey_xy_xxz_1, tey_xy_xxzz_0, \
                                         tey_xy_xxzz_1, tey_xy_xyy_0, tey_xy_xyy_1, tey_xy_xyyy_0, tey_xy_xyyy_1, \
                                         tey_xy_xyyz_0, tey_xy_xyyz_1, tey_xy_xyz_0, tey_xy_xyz_1, tey_xy_xyzz_0, \
                                         tey_xy_xyzz_1, tey_xy_xzz_0, tey_xy_xzz_1, tey_xy_xzzz_0, tey_xy_xzzz_1, \
                                         tey_xy_yyy_0, tey_xy_yyy_1, tey_xy_yyyy_0, tey_xy_yyyy_1, tey_xy_yyyz_0, \
                                         tey_xy_yyyz_1, tey_xy_yyz_0, tey_xy_yyz_1, tey_xy_yyzz_0, tey_xy_yyzz_1, \
                                         tey_xy_yzz_0, tey_xy_yzz_1, tey_xy_yzzz_0, tey_xy_yzzz_1, tey_xy_zzz_0, \
                                         tey_xy_zzz_1, tey_xy_zzzz_0, tey_xy_zzzz_1, tey_xz_xxx_0, tey_xz_xxx_1, \
                                         tey_xz_xxxx_0, tey_xz_xxxx_1, tey_xz_xxxy_0, tey_xz_xxxy_1, tey_xz_xxxz_0, \
                                         tey_xz_xxxz_1, tey_xz_xxy_0, tey_xz_xxy_1, tey_xz_xxz_0, tey_xz_xxz_1, tey_y_xxxz_0, \
                                         tey_y_xxxz_1, tey_y_xxyy_0, tey_y_xxyy_1, tey_y_xxyz_0, tey_y_xxyz_1, tey_y_xxzz_0, \
                                         tey_y_xxzz_1, tey_y_xyyy_0, tey_y_xyyy_1, tey_y_xyyz_0, tey_y_xyyz_1, tey_y_xyzz_0, \
                                         tey_y_xyzz_1, tey_y_xzzz_0, tey_y_xzzz_1, tey_y_yyyy_0, tey_y_yyyy_1, tey_y_yyyz_0, \
                                         tey_y_yyyz_1, tey_y_yyzz_0, tey_y_yyzz_1, tey_y_yzzz_0, tey_y_yzzz_1, tey_y_zzzz_0, \
                                         tey_y_zzzz_1, tey_z_xxxx_0, tey_z_xxxx_1, tey_z_xxxy_0, tey_z_xxxy_1, tey_z_xxxz_0, \
                                         tey_z_xxxz_1, tez_xxy_xxxy_0, tez_xxy_xxxz_0, tez_xxy_xxyy_0, tez_xxy_xxyz_0, \
                                         tez_xxy_xxzz_0, tez_xxy_xyyy_0, tez_xxy_xyyz_0, tez_xxy_xyzz_0, tez_xxy_xzzz_0, \
                                         tez_xxy_yyyy_0, tez_xxy_yyyz_0, tez_xxy_yyzz_0, tez_xxy_yzzz_0, tez_xxy_zzzz_0, \
                                         tez_xxz_xxxx_0, tez_xxz_xxxy_0, tez_xxz_xxxz_0, tez_xy_xxxy_0, tez_xy_xxxy_1, \
                                         tez_xy_xxxz_0, tez_xy_xxxz_1, tez_xy_xxy_0, tez_xy_xxy_1, tez_xy_xxyy_0, \
                                         tez_xy_xxyy_1, tez_xy_xxyz_0, tez_xy_xxyz_1, tez_xy_xxz_0, tez_xy_xxz_1, \
                                         tez_xy_xxzz_0, tez_xy_xxzz_1, tez_xy_xyy_0, tez_xy_xyy_1, tez_xy_xyyy_0, \
                                         tez_xy_xyyy_1, tez_xy_xyyz_0, tez_xy_xyyz_1, tez_xy_xyz_0, tez_xy_xyz_1, \
                                         tez_xy_xyzz_0, tez_xy_xyzz_1, tez_xy_xzz_0, tez_xy_xzz_1, tez_xy_xzzz_0, \
                                         tez_xy_xzzz_1, tez_xy_yyy_0, tez_xy_yyy_1, tez_xy_yyyy_0, tez_xy_yyyy_1, \
                                         tez_xy_yyyz_0, tez_xy_yyyz_1, tez_xy_yyz_0, tez_xy_yyz_1, tez_xy_yyzz_0, \
                                         tez_xy_yyzz_1, tez_xy_yzz_0, tez_xy_yzz_1, tez_xy_yzzz_0, tez_xy_yzzz_1, \
                                         tez_xy_zzz_0, tez_xy_zzz_1, tez_xy_zzzz_0, tez_xy_zzzz_1, tez_xz_xxx_0, \
                                         tez_xz_xxx_1, tez_xz_xxxx_0, tez_xz_xxxx_1, tez_xz_xxxy_0, tez_xz_xxxy_1, \
                                         tez_xz_xxxz_0, tez_xz_xxxz_1, tez_xz_xxy_0, tez_xz_xxy_1, tez_xz_xxz_0, \
                                         tez_xz_xxz_1, tez_y_xxxy_0, tez_y_xxxy_1, tez_y_xxxz_0, tez_y_xxxz_1, tez_y_xxyy_0, \
                                         tez_y_xxyy_1, tez_y_xxyz_0, tez_y_xxyz_1, tez_y_xxzz_0, tez_y_xxzz_1, tez_y_xyyy_0, \
                                         tez_y_xyyy_1, tez_y_xyyz_0, tez_y_xyyz_1, tez_y_xyzz_0, tez_y_xyzz_1, tez_y_xzzz_0, \
                                         tez_y_xzzz_1, tez_y_yyyy_0, tez_y_yyyy_1, tez_y_yyyz_0, tez_y_yyyz_1, tez_y_yyzz_0, \
                                         tez_y_yyzz_1, tez_y_yzzz_0, tez_y_yzzz_1, tez_y_zzzz_0, tez_y_zzzz_1, tez_z_xxxx_0, \
                                         tez_z_xxxx_1, tez_z_xxxy_0, tez_z_xxxy_1, tez_z_xxxz_0, tez_z_xxxz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tez_xxy_xxxy_0[j] = pa_x[j] * tez_xy_xxxy_0[j] - pc_x[j] * tez_xy_xxxy_1[j] + 0.5 * fl1_fx * tez_y_xxxy_0[j] - 0.5 * fl1_fx * tez_y_xxxy_1[j] + 1.5 * fl1_fx * tez_xy_xxy_0[j] - 1.5 * fl1_fx * tez_xy_xxy_1[j];

                    tex_xxy_xxxz_0[j] = pa_x[j] * tex_xy_xxxz_0[j] - pc_x[j] * tex_xy_xxxz_1[j] + 0.5 * fl1_fx * tex_y_xxxz_0[j] - 0.5 * fl1_fx * tex_y_xxxz_1[j] + 1.5 * fl1_fx * tex_xy_xxz_0[j] - 1.5 * fl1_fx * tex_xy_xxz_1[j] + ta_xy_xxxz_1[j];

                    tey_xxy_xxxz_0[j] = pa_x[j] * tey_xy_xxxz_0[j] - pc_x[j] * tey_xy_xxxz_1[j] + 0.5 * fl1_fx * tey_y_xxxz_0[j] - 0.5 * fl1_fx * tey_y_xxxz_1[j] + 1.5 * fl1_fx * tey_xy_xxz_0[j] - 1.5 * fl1_fx * tey_xy_xxz_1[j];

                    tez_xxy_xxxz_0[j] = pa_x[j] * tez_xy_xxxz_0[j] - pc_x[j] * tez_xy_xxxz_1[j] + 0.5 * fl1_fx * tez_y_xxxz_0[j] - 0.5 * fl1_fx * tez_y_xxxz_1[j] + 1.5 * fl1_fx * tez_xy_xxz_0[j] - 1.5 * fl1_fx * tez_xy_xxz_1[j];

                    tex_xxy_xxyy_0[j] = pa_x[j] * tex_xy_xxyy_0[j] - pc_x[j] * tex_xy_xxyy_1[j] + 0.5 * fl1_fx * tex_y_xxyy_0[j] - 0.5 * fl1_fx * tex_y_xxyy_1[j] + fl1_fx * tex_xy_xyy_0[j] - fl1_fx * tex_xy_xyy_1[j] + ta_xy_xxyy_1[j];

                    tey_xxy_xxyy_0[j] = pa_x[j] * tey_xy_xxyy_0[j] - pc_x[j] * tey_xy_xxyy_1[j] + 0.5 * fl1_fx * tey_y_xxyy_0[j] - 0.5 * fl1_fx * tey_y_xxyy_1[j] + fl1_fx * tey_xy_xyy_0[j] - fl1_fx * tey_xy_xyy_1[j];

                    tez_xxy_xxyy_0[j] = pa_x[j] * tez_xy_xxyy_0[j] - pc_x[j] * tez_xy_xxyy_1[j] + 0.5 * fl1_fx * tez_y_xxyy_0[j] - 0.5 * fl1_fx * tez_y_xxyy_1[j] + fl1_fx * tez_xy_xyy_0[j] - fl1_fx * tez_xy_xyy_1[j];

                    tex_xxy_xxyz_0[j] = pa_x[j] * tex_xy_xxyz_0[j] - pc_x[j] * tex_xy_xxyz_1[j] + 0.5 * fl1_fx * tex_y_xxyz_0[j] - 0.5 * fl1_fx * tex_y_xxyz_1[j] + fl1_fx * tex_xy_xyz_0[j] - fl1_fx * tex_xy_xyz_1[j] + ta_xy_xxyz_1[j];

                    tey_xxy_xxyz_0[j] = pa_x[j] * tey_xy_xxyz_0[j] - pc_x[j] * tey_xy_xxyz_1[j] + 0.5 * fl1_fx * tey_y_xxyz_0[j] - 0.5 * fl1_fx * tey_y_xxyz_1[j] + fl1_fx * tey_xy_xyz_0[j] - fl1_fx * tey_xy_xyz_1[j];

                    tez_xxy_xxyz_0[j] = pa_x[j] * tez_xy_xxyz_0[j] - pc_x[j] * tez_xy_xxyz_1[j] + 0.5 * fl1_fx * tez_y_xxyz_0[j] - 0.5 * fl1_fx * tez_y_xxyz_1[j] + fl1_fx * tez_xy_xyz_0[j] - fl1_fx * tez_xy_xyz_1[j];

                    tex_xxy_xxzz_0[j] = pa_x[j] * tex_xy_xxzz_0[j] - pc_x[j] * tex_xy_xxzz_1[j] + 0.5 * fl1_fx * tex_y_xxzz_0[j] - 0.5 * fl1_fx * tex_y_xxzz_1[j] + fl1_fx * tex_xy_xzz_0[j] - fl1_fx * tex_xy_xzz_1[j] + ta_xy_xxzz_1[j];

                    tey_xxy_xxzz_0[j] = pa_x[j] * tey_xy_xxzz_0[j] - pc_x[j] * tey_xy_xxzz_1[j] + 0.5 * fl1_fx * tey_y_xxzz_0[j] - 0.5 * fl1_fx * tey_y_xxzz_1[j] + fl1_fx * tey_xy_xzz_0[j] - fl1_fx * tey_xy_xzz_1[j];

                    tez_xxy_xxzz_0[j] = pa_x[j] * tez_xy_xxzz_0[j] - pc_x[j] * tez_xy_xxzz_1[j] + 0.5 * fl1_fx * tez_y_xxzz_0[j] - 0.5 * fl1_fx * tez_y_xxzz_1[j] + fl1_fx * tez_xy_xzz_0[j] - fl1_fx * tez_xy_xzz_1[j];

                    tex_xxy_xyyy_0[j] = pa_x[j] * tex_xy_xyyy_0[j] - pc_x[j] * tex_xy_xyyy_1[j] + 0.5 * fl1_fx * tex_y_xyyy_0[j] - 0.5 * fl1_fx * tex_y_xyyy_1[j] + 0.5 * fl1_fx * tex_xy_yyy_0[j] - 0.5 * fl1_fx * tex_xy_yyy_1[j] + ta_xy_xyyy_1[j];

                    tey_xxy_xyyy_0[j] = pa_x[j] * tey_xy_xyyy_0[j] - pc_x[j] * tey_xy_xyyy_1[j] + 0.5 * fl1_fx * tey_y_xyyy_0[j] - 0.5 * fl1_fx * tey_y_xyyy_1[j] + 0.5 * fl1_fx * tey_xy_yyy_0[j] - 0.5 * fl1_fx * tey_xy_yyy_1[j];

                    tez_xxy_xyyy_0[j] = pa_x[j] * tez_xy_xyyy_0[j] - pc_x[j] * tez_xy_xyyy_1[j] + 0.5 * fl1_fx * tez_y_xyyy_0[j] - 0.5 * fl1_fx * tez_y_xyyy_1[j] + 0.5 * fl1_fx * tez_xy_yyy_0[j] - 0.5 * fl1_fx * tez_xy_yyy_1[j];

                    tex_xxy_xyyz_0[j] = pa_x[j] * tex_xy_xyyz_0[j] - pc_x[j] * tex_xy_xyyz_1[j] + 0.5 * fl1_fx * tex_y_xyyz_0[j] - 0.5 * fl1_fx * tex_y_xyyz_1[j] + 0.5 * fl1_fx * tex_xy_yyz_0[j] - 0.5 * fl1_fx * tex_xy_yyz_1[j] + ta_xy_xyyz_1[j];

                    tey_xxy_xyyz_0[j] = pa_x[j] * tey_xy_xyyz_0[j] - pc_x[j] * tey_xy_xyyz_1[j] + 0.5 * fl1_fx * tey_y_xyyz_0[j] - 0.5 * fl1_fx * tey_y_xyyz_1[j] + 0.5 * fl1_fx * tey_xy_yyz_0[j] - 0.5 * fl1_fx * tey_xy_yyz_1[j];

                    tez_xxy_xyyz_0[j] = pa_x[j] * tez_xy_xyyz_0[j] - pc_x[j] * tez_xy_xyyz_1[j] + 0.5 * fl1_fx * tez_y_xyyz_0[j] - 0.5 * fl1_fx * tez_y_xyyz_1[j] + 0.5 * fl1_fx * tez_xy_yyz_0[j] - 0.5 * fl1_fx * tez_xy_yyz_1[j];

                    tex_xxy_xyzz_0[j] = pa_x[j] * tex_xy_xyzz_0[j] - pc_x[j] * tex_xy_xyzz_1[j] + 0.5 * fl1_fx * tex_y_xyzz_0[j] - 0.5 * fl1_fx * tex_y_xyzz_1[j] + 0.5 * fl1_fx * tex_xy_yzz_0[j] - 0.5 * fl1_fx * tex_xy_yzz_1[j] + ta_xy_xyzz_1[j];

                    tey_xxy_xyzz_0[j] = pa_x[j] * tey_xy_xyzz_0[j] - pc_x[j] * tey_xy_xyzz_1[j] + 0.5 * fl1_fx * tey_y_xyzz_0[j] - 0.5 * fl1_fx * tey_y_xyzz_1[j] + 0.5 * fl1_fx * tey_xy_yzz_0[j] - 0.5 * fl1_fx * tey_xy_yzz_1[j];

                    tez_xxy_xyzz_0[j] = pa_x[j] * tez_xy_xyzz_0[j] - pc_x[j] * tez_xy_xyzz_1[j] + 0.5 * fl1_fx * tez_y_xyzz_0[j] - 0.5 * fl1_fx * tez_y_xyzz_1[j] + 0.5 * fl1_fx * tez_xy_yzz_0[j] - 0.5 * fl1_fx * tez_xy_yzz_1[j];

                    tex_xxy_xzzz_0[j] = pa_x[j] * tex_xy_xzzz_0[j] - pc_x[j] * tex_xy_xzzz_1[j] + 0.5 * fl1_fx * tex_y_xzzz_0[j] - 0.5 * fl1_fx * tex_y_xzzz_1[j] + 0.5 * fl1_fx * tex_xy_zzz_0[j] - 0.5 * fl1_fx * tex_xy_zzz_1[j] + ta_xy_xzzz_1[j];

                    tey_xxy_xzzz_0[j] = pa_x[j] * tey_xy_xzzz_0[j] - pc_x[j] * tey_xy_xzzz_1[j] + 0.5 * fl1_fx * tey_y_xzzz_0[j] - 0.5 * fl1_fx * tey_y_xzzz_1[j] + 0.5 * fl1_fx * tey_xy_zzz_0[j] - 0.5 * fl1_fx * tey_xy_zzz_1[j];

                    tez_xxy_xzzz_0[j] = pa_x[j] * tez_xy_xzzz_0[j] - pc_x[j] * tez_xy_xzzz_1[j] + 0.5 * fl1_fx * tez_y_xzzz_0[j] - 0.5 * fl1_fx * tez_y_xzzz_1[j] + 0.5 * fl1_fx * tez_xy_zzz_0[j] - 0.5 * fl1_fx * tez_xy_zzz_1[j];

                    tex_xxy_yyyy_0[j] = pa_x[j] * tex_xy_yyyy_0[j] - pc_x[j] * tex_xy_yyyy_1[j] + 0.5 * fl1_fx * tex_y_yyyy_0[j] - 0.5 * fl1_fx * tex_y_yyyy_1[j] + ta_xy_yyyy_1[j];

                    tey_xxy_yyyy_0[j] = pa_x[j] * tey_xy_yyyy_0[j] - pc_x[j] * tey_xy_yyyy_1[j] + 0.5 * fl1_fx * tey_y_yyyy_0[j] - 0.5 * fl1_fx * tey_y_yyyy_1[j];

                    tez_xxy_yyyy_0[j] = pa_x[j] * tez_xy_yyyy_0[j] - pc_x[j] * tez_xy_yyyy_1[j] + 0.5 * fl1_fx * tez_y_yyyy_0[j] - 0.5 * fl1_fx * tez_y_yyyy_1[j];

                    tex_xxy_yyyz_0[j] = pa_x[j] * tex_xy_yyyz_0[j] - pc_x[j] * tex_xy_yyyz_1[j] + 0.5 * fl1_fx * tex_y_yyyz_0[j] - 0.5 * fl1_fx * tex_y_yyyz_1[j] + ta_xy_yyyz_1[j];

                    tey_xxy_yyyz_0[j] = pa_x[j] * tey_xy_yyyz_0[j] - pc_x[j] * tey_xy_yyyz_1[j] + 0.5 * fl1_fx * tey_y_yyyz_0[j] - 0.5 * fl1_fx * tey_y_yyyz_1[j];

                    tez_xxy_yyyz_0[j] = pa_x[j] * tez_xy_yyyz_0[j] - pc_x[j] * tez_xy_yyyz_1[j] + 0.5 * fl1_fx * tez_y_yyyz_0[j] - 0.5 * fl1_fx * tez_y_yyyz_1[j];

                    tex_xxy_yyzz_0[j] = pa_x[j] * tex_xy_yyzz_0[j] - pc_x[j] * tex_xy_yyzz_1[j] + 0.5 * fl1_fx * tex_y_yyzz_0[j] - 0.5 * fl1_fx * tex_y_yyzz_1[j] + ta_xy_yyzz_1[j];

                    tey_xxy_yyzz_0[j] = pa_x[j] * tey_xy_yyzz_0[j] - pc_x[j] * tey_xy_yyzz_1[j] + 0.5 * fl1_fx * tey_y_yyzz_0[j] - 0.5 * fl1_fx * tey_y_yyzz_1[j];

                    tez_xxy_yyzz_0[j] = pa_x[j] * tez_xy_yyzz_0[j] - pc_x[j] * tez_xy_yyzz_1[j] + 0.5 * fl1_fx * tez_y_yyzz_0[j] - 0.5 * fl1_fx * tez_y_yyzz_1[j];

                    tex_xxy_yzzz_0[j] = pa_x[j] * tex_xy_yzzz_0[j] - pc_x[j] * tex_xy_yzzz_1[j] + 0.5 * fl1_fx * tex_y_yzzz_0[j] - 0.5 * fl1_fx * tex_y_yzzz_1[j] + ta_xy_yzzz_1[j];

                    tey_xxy_yzzz_0[j] = pa_x[j] * tey_xy_yzzz_0[j] - pc_x[j] * tey_xy_yzzz_1[j] + 0.5 * fl1_fx * tey_y_yzzz_0[j] - 0.5 * fl1_fx * tey_y_yzzz_1[j];

                    tez_xxy_yzzz_0[j] = pa_x[j] * tez_xy_yzzz_0[j] - pc_x[j] * tez_xy_yzzz_1[j] + 0.5 * fl1_fx * tez_y_yzzz_0[j] - 0.5 * fl1_fx * tez_y_yzzz_1[j];

                    tex_xxy_zzzz_0[j] = pa_x[j] * tex_xy_zzzz_0[j] - pc_x[j] * tex_xy_zzzz_1[j] + 0.5 * fl1_fx * tex_y_zzzz_0[j] - 0.5 * fl1_fx * tex_y_zzzz_1[j] + ta_xy_zzzz_1[j];

                    tey_xxy_zzzz_0[j] = pa_x[j] * tey_xy_zzzz_0[j] - pc_x[j] * tey_xy_zzzz_1[j] + 0.5 * fl1_fx * tey_y_zzzz_0[j] - 0.5 * fl1_fx * tey_y_zzzz_1[j];

                    tez_xxy_zzzz_0[j] = pa_x[j] * tez_xy_zzzz_0[j] - pc_x[j] * tez_xy_zzzz_1[j] + 0.5 * fl1_fx * tez_y_zzzz_0[j] - 0.5 * fl1_fx * tez_y_zzzz_1[j];

                    tex_xxz_xxxx_0[j] = pa_x[j] * tex_xz_xxxx_0[j] - pc_x[j] * tex_xz_xxxx_1[j] + 0.5 * fl1_fx * tex_z_xxxx_0[j] - 0.5 * fl1_fx * tex_z_xxxx_1[j] + 2.0 * fl1_fx * tex_xz_xxx_0[j] - 2.0 * fl1_fx * tex_xz_xxx_1[j] + ta_xz_xxxx_1[j];

                    tey_xxz_xxxx_0[j] = pa_x[j] * tey_xz_xxxx_0[j] - pc_x[j] * tey_xz_xxxx_1[j] + 0.5 * fl1_fx * tey_z_xxxx_0[j] - 0.5 * fl1_fx * tey_z_xxxx_1[j] + 2.0 * fl1_fx * tey_xz_xxx_0[j] - 2.0 * fl1_fx * tey_xz_xxx_1[j];

                    tez_xxz_xxxx_0[j] = pa_x[j] * tez_xz_xxxx_0[j] - pc_x[j] * tez_xz_xxxx_1[j] + 0.5 * fl1_fx * tez_z_xxxx_0[j] - 0.5 * fl1_fx * tez_z_xxxx_1[j] + 2.0 * fl1_fx * tez_xz_xxx_0[j] - 2.0 * fl1_fx * tez_xz_xxx_1[j];

                    tex_xxz_xxxy_0[j] = pa_x[j] * tex_xz_xxxy_0[j] - pc_x[j] * tex_xz_xxxy_1[j] + 0.5 * fl1_fx * tex_z_xxxy_0[j] - 0.5 * fl1_fx * tex_z_xxxy_1[j] + 1.5 * fl1_fx * tex_xz_xxy_0[j] - 1.5 * fl1_fx * tex_xz_xxy_1[j] + ta_xz_xxxy_1[j];

                    tey_xxz_xxxy_0[j] = pa_x[j] * tey_xz_xxxy_0[j] - pc_x[j] * tey_xz_xxxy_1[j] + 0.5 * fl1_fx * tey_z_xxxy_0[j] - 0.5 * fl1_fx * tey_z_xxxy_1[j] + 1.5 * fl1_fx * tey_xz_xxy_0[j] - 1.5 * fl1_fx * tey_xz_xxy_1[j];

                    tez_xxz_xxxy_0[j] = pa_x[j] * tez_xz_xxxy_0[j] - pc_x[j] * tez_xz_xxxy_1[j] + 0.5 * fl1_fx * tez_z_xxxy_0[j] - 0.5 * fl1_fx * tez_z_xxxy_1[j] + 1.5 * fl1_fx * tez_xz_xxy_0[j] - 1.5 * fl1_fx * tez_xz_xxy_1[j];

                    tex_xxz_xxxz_0[j] = pa_x[j] * tex_xz_xxxz_0[j] - pc_x[j] * tex_xz_xxxz_1[j] + 0.5 * fl1_fx * tex_z_xxxz_0[j] - 0.5 * fl1_fx * tex_z_xxxz_1[j] + 1.5 * fl1_fx * tex_xz_xxz_0[j] - 1.5 * fl1_fx * tex_xz_xxz_1[j] + ta_xz_xxxz_1[j];

                    tey_xxz_xxxz_0[j] = pa_x[j] * tey_xz_xxxz_0[j] - pc_x[j] * tey_xz_xxxz_1[j] + 0.5 * fl1_fx * tey_z_xxxz_0[j] - 0.5 * fl1_fx * tey_z_xxxz_1[j] + 1.5 * fl1_fx * tey_xz_xxz_0[j] - 1.5 * fl1_fx * tey_xz_xxz_1[j];

                    tez_xxz_xxxz_0[j] = pa_x[j] * tez_xz_xxxz_0[j] - pc_x[j] * tez_xz_xxxz_1[j] + 0.5 * fl1_fx * tez_z_xxxz_0[j] - 0.5 * fl1_fx * tez_z_xxxz_1[j] + 1.5 * fl1_fx * tez_xz_xxz_0[j] - 1.5 * fl1_fx * tez_xz_xxz_1[j];

                    tex_xxz_xxyy_0[j] = pa_x[j] * tex_xz_xxyy_0[j] - pc_x[j] * tex_xz_xxyy_1[j] + 0.5 * fl1_fx * tex_z_xxyy_0[j] - 0.5 * fl1_fx * tex_z_xxyy_1[j] + fl1_fx * tex_xz_xyy_0[j] - fl1_fx * tex_xz_xyy_1[j] + ta_xz_xxyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_100_150(      CMemBlock2D<double>& primBuffer,
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

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tey_xz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 33); 

                auto tez_xz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 33); 

                auto tex_xz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 34); 

                auto tey_xz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 34); 

                auto tez_xz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 34); 

                auto tex_xz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 35); 

                auto tey_xz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 35); 

                auto tez_xz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 35); 

                auto tex_xz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 36); 

                auto tey_xz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 36); 

                auto tez_xz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 36); 

                auto tex_xz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 37); 

                auto tey_xz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 37); 

                auto tez_xz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 37); 

                auto tex_xz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 38); 

                auto tey_xz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 38); 

                auto tez_xz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 38); 

                auto tex_xz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 39); 

                auto tey_xz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 39); 

                auto tez_xz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 39); 

                auto tex_xz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 40); 

                auto tey_xz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 40); 

                auto tez_xz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 40); 

                auto tex_xz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 41); 

                auto tey_xz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 41); 

                auto tez_xz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 41); 

                auto tex_xz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 42); 

                auto tey_xz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 42); 

                auto tez_xz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 42); 

                auto tex_xz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 43); 

                auto tey_xz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 43); 

                auto tez_xz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 43); 

                auto tex_xz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 44); 

                auto tey_xz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 44); 

                auto tez_xz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 44); 

                auto tex_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 45); 

                auto tey_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 45); 

                auto tez_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 45); 

                auto tex_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 46); 

                auto tey_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 46); 

                auto tez_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 46); 

                auto tex_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 47); 

                auto tey_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 47); 

                auto tez_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 47); 

                auto tex_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 48); 

                auto tey_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 48); 

                auto tez_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 48); 

                auto tex_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 49); 

                auto tey_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 49); 

                auto tez_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 49); 

                auto tey_xz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 33); 

                auto tez_xz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 33); 

                auto tex_xz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 34); 

                auto tey_xz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 34); 

                auto tez_xz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 34); 

                auto tex_xz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 35); 

                auto tey_xz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 35); 

                auto tez_xz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 35); 

                auto tex_xz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 36); 

                auto tey_xz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 36); 

                auto tez_xz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 36); 

                auto tex_xz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 37); 

                auto tey_xz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 37); 

                auto tez_xz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 37); 

                auto tex_xz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 38); 

                auto tey_xz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 38); 

                auto tez_xz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 38); 

                auto tex_xz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 39); 

                auto tey_xz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 39); 

                auto tez_xz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 39); 

                auto tex_xz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 40); 

                auto tey_xz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 40); 

                auto tez_xz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 40); 

                auto tex_xz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 41); 

                auto tey_xz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 41); 

                auto tez_xz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 41); 

                auto tex_xz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 42); 

                auto tey_xz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 42); 

                auto tez_xz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 42); 

                auto tex_xz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 43); 

                auto tey_xz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 43); 

                auto tez_xz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 43); 

                auto tex_xz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 44); 

                auto tey_xz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 44); 

                auto tez_xz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 44); 

                auto tex_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 45); 

                auto tey_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 45); 

                auto tez_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 45); 

                auto tex_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 46); 

                auto tey_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 46); 

                auto tez_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 46); 

                auto tex_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 47); 

                auto tey_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 47); 

                auto tez_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 47); 

                auto tex_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 48); 

                auto tey_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 48); 

                auto tez_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 48); 

                auto tex_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 49); 

                auto tey_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 49); 

                auto tez_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 49); 

                auto tey_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 33); 

                auto tez_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 33); 

                auto tex_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 34); 

                auto tey_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 34); 

                auto tez_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 34); 

                auto tex_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 35); 

                auto tey_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 35); 

                auto tez_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 35); 

                auto tex_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 36); 

                auto tey_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 36); 

                auto tez_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 36); 

                auto tex_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 37); 

                auto tey_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 37); 

                auto tez_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 37); 

                auto tex_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 38); 

                auto tey_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 38); 

                auto tez_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 38); 

                auto tex_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 39); 

                auto tey_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 39); 

                auto tez_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 39); 

                auto tex_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 40); 

                auto tey_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 40); 

                auto tez_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 40); 

                auto tex_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 41); 

                auto tey_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 41); 

                auto tez_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 41); 

                auto tex_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 42); 

                auto tey_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 42); 

                auto tez_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 42); 

                auto tex_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 43); 

                auto tey_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 43); 

                auto tez_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 43); 

                auto tex_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 44); 

                auto tey_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 44); 

                auto tez_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 44); 

                auto tey_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 33); 

                auto tez_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 33); 

                auto tex_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 34); 

                auto tey_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 34); 

                auto tez_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 34); 

                auto tex_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 35); 

                auto tey_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 35); 

                auto tez_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 35); 

                auto tex_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 36); 

                auto tey_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 36); 

                auto tez_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 36); 

                auto tex_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 37); 

                auto tey_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 37); 

                auto tez_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 37); 

                auto tex_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 38); 

                auto tey_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 38); 

                auto tez_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 38); 

                auto tex_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 39); 

                auto tey_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 39); 

                auto tez_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 39); 

                auto tex_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 40); 

                auto tey_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 40); 

                auto tez_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 40); 

                auto tex_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 41); 

                auto tey_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 41); 

                auto tez_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 41); 

                auto tex_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 42); 

                auto tey_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 42); 

                auto tez_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 42); 

                auto tex_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 43); 

                auto tey_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 43); 

                auto tez_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 43); 

                auto tex_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 44); 

                auto tey_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 44); 

                auto tez_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 44); 

                auto tey_xz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 23); 

                auto tez_xz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 23); 

                auto tex_xz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 24); 

                auto tey_xz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 24); 

                auto tez_xz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 24); 

                auto tex_xz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 25); 

                auto tey_xz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 25); 

                auto tez_xz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 25); 

                auto tex_xz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 26); 

                auto tey_xz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 26); 

                auto tez_xz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 26); 

                auto tex_xz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 27); 

                auto tey_xz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 27); 

                auto tez_xz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 27); 

                auto tex_xz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 28); 

                auto tey_xz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 28); 

                auto tez_xz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 28); 

                auto tex_xz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 29); 

                auto tey_xz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 29); 

                auto tez_xz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 29); 

                auto tex_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 30); 

                auto tey_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 30); 

                auto tez_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 30); 

                auto tex_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 31); 

                auto tey_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 31); 

                auto tez_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 31); 

                auto tex_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 32); 

                auto tey_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 32); 

                auto tez_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 32); 

                auto tex_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 33); 

                auto tey_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 33); 

                auto tez_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 33); 

                auto tex_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 34); 

                auto tey_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 34); 

                auto tez_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 34); 

                auto tey_xz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 23); 

                auto tez_xz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 23); 

                auto tex_xz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 24); 

                auto tey_xz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 24); 

                auto tez_xz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 24); 

                auto tex_xz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 25); 

                auto tey_xz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 25); 

                auto tez_xz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 25); 

                auto tex_xz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 26); 

                auto tey_xz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 26); 

                auto tez_xz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 26); 

                auto tex_xz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 27); 

                auto tey_xz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 27); 

                auto tez_xz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 27); 

                auto tex_xz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 28); 

                auto tey_xz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 28); 

                auto tez_xz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 28); 

                auto tex_xz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 29); 

                auto tey_xz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 29); 

                auto tez_xz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 29); 

                auto tex_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 30); 

                auto tey_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 30); 

                auto tez_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 30); 

                auto tex_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 31); 

                auto tey_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 31); 

                auto tez_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 31); 

                auto tex_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 32); 

                auto tey_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 32); 

                auto tez_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 32); 

                auto tex_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 33); 

                auto tey_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 33); 

                auto tez_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 33); 

                auto tex_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 34); 

                auto tey_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 34); 

                auto tez_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 34); 

                auto ta_xz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 34); 

                auto ta_xz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 35); 

                auto ta_xz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 36); 

                auto ta_xz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 37); 

                auto ta_xz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 38); 

                auto ta_xz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 39); 

                auto ta_xz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 40); 

                auto ta_xz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 41); 

                auto ta_xz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 42); 

                auto ta_xz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 43); 

                auto ta_xz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 44); 

                auto ta_yy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 45); 

                auto ta_yy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 46); 

                auto ta_yy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 47); 

                auto ta_yy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 48); 

                auto ta_yy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 49); 

                // set up pointers to integrals

                auto tey_xxz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 33); 

                auto tez_xxz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 33); 

                auto tex_xxz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 34); 

                auto tey_xxz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 34); 

                auto tez_xxz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 34); 

                auto tex_xxz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 35); 

                auto tey_xxz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 35); 

                auto tez_xxz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 35); 

                auto tex_xxz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 36); 

                auto tey_xxz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 36); 

                auto tez_xxz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 36); 

                auto tex_xxz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 37); 

                auto tey_xxz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 37); 

                auto tez_xxz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 37); 

                auto tex_xxz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 38); 

                auto tey_xxz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 38); 

                auto tez_xxz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 38); 

                auto tex_xxz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 39); 

                auto tey_xxz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 39); 

                auto tez_xxz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 39); 

                auto tex_xxz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 40); 

                auto tey_xxz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 40); 

                auto tez_xxz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 40); 

                auto tex_xxz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 41); 

                auto tey_xxz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 41); 

                auto tez_xxz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 41); 

                auto tex_xxz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 42); 

                auto tey_xxz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 42); 

                auto tez_xxz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 42); 

                auto tex_xxz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 43); 

                auto tey_xxz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 43); 

                auto tez_xxz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 43); 

                auto tex_xxz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 44); 

                auto tey_xxz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 44); 

                auto tez_xxz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 44); 

                auto tex_xyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 45); 

                auto tey_xyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 45); 

                auto tez_xyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 45); 

                auto tex_xyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 46); 

                auto tey_xyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 46); 

                auto tez_xyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 46); 

                auto tex_xyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 47); 

                auto tey_xyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 47); 

                auto tez_xyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 47); 

                auto tex_xyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 48); 

                auto tey_xyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 48); 

                auto tez_xyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 48); 

                auto tex_xyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 49); 

                auto tey_xyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 49); 

                auto tez_xyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 49); 

                // Batch of Integrals (100,150)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xz_xxyz_1, ta_xz_xxzz_1, ta_xz_xyyy_1, ta_xz_xyyz_1, \
                                         ta_xz_xyzz_1, ta_xz_xzzz_1, ta_xz_yyyy_1, ta_xz_yyyz_1, ta_xz_yyzz_1, ta_xz_yzzz_1, \
                                         ta_xz_zzzz_1, ta_yy_xxxx_1, ta_yy_xxxy_1, ta_yy_xxxz_1, ta_yy_xxyy_1, ta_yy_xxyz_1, \
                                         tex_xxz_xxyz_0, tex_xxz_xxzz_0, tex_xxz_xyyy_0, tex_xxz_xyyz_0, tex_xxz_xyzz_0, \
                                         tex_xxz_xzzz_0, tex_xxz_yyyy_0, tex_xxz_yyyz_0, tex_xxz_yyzz_0, tex_xxz_yzzz_0, \
                                         tex_xxz_zzzz_0, tex_xyy_xxxx_0, tex_xyy_xxxy_0, tex_xyy_xxxz_0, tex_xyy_xxyy_0, \
                                         tex_xyy_xxyz_0, tex_xz_xxyz_0, tex_xz_xxyz_1, tex_xz_xxzz_0, tex_xz_xxzz_1, \
                                         tex_xz_xyyy_0, tex_xz_xyyy_1, tex_xz_xyyz_0, tex_xz_xyyz_1, tex_xz_xyz_0, \
                                         tex_xz_xyz_1, tex_xz_xyzz_0, tex_xz_xyzz_1, tex_xz_xzz_0, tex_xz_xzz_1, \
                                         tex_xz_xzzz_0, tex_xz_xzzz_1, tex_xz_yyy_0, tex_xz_yyy_1, tex_xz_yyyy_0, \
                                         tex_xz_yyyy_1, tex_xz_yyyz_0, tex_xz_yyyz_1, tex_xz_yyz_0, tex_xz_yyz_1, \
                                         tex_xz_yyzz_0, tex_xz_yyzz_1, tex_xz_yzz_0, tex_xz_yzz_1, tex_xz_yzzz_0, \
                                         tex_xz_yzzz_1, tex_xz_zzz_0, tex_xz_zzz_1, tex_xz_zzzz_0, tex_xz_zzzz_1, \
                                         tex_yy_xxx_0, tex_yy_xxx_1, tex_yy_xxxx_0, tex_yy_xxxx_1, tex_yy_xxxy_0, \
                                         tex_yy_xxxy_1, tex_yy_xxxz_0, tex_yy_xxxz_1, tex_yy_xxy_0, tex_yy_xxy_1, \
                                         tex_yy_xxyy_0, tex_yy_xxyy_1, tex_yy_xxyz_0, tex_yy_xxyz_1, tex_yy_xxz_0, \
                                         tex_yy_xxz_1, tex_yy_xyy_0, tex_yy_xyy_1, tex_yy_xyz_0, tex_yy_xyz_1, tex_z_xxyz_0, \
                                         tex_z_xxyz_1, tex_z_xxzz_0, tex_z_xxzz_1, tex_z_xyyy_0, tex_z_xyyy_1, tex_z_xyyz_0, \
                                         tex_z_xyyz_1, tex_z_xyzz_0, tex_z_xyzz_1, tex_z_xzzz_0, tex_z_xzzz_1, tex_z_yyyy_0, \
                                         tex_z_yyyy_1, tex_z_yyyz_0, tex_z_yyyz_1, tex_z_yyzz_0, tex_z_yyzz_1, tex_z_yzzz_0, \
                                         tex_z_yzzz_1, tex_z_zzzz_0, tex_z_zzzz_1, tey_xxz_xxyy_0, tey_xxz_xxyz_0, \
                                         tey_xxz_xxzz_0, tey_xxz_xyyy_0, tey_xxz_xyyz_0, tey_xxz_xyzz_0, tey_xxz_xzzz_0, \
                                         tey_xxz_yyyy_0, tey_xxz_yyyz_0, tey_xxz_yyzz_0, tey_xxz_yzzz_0, tey_xxz_zzzz_0, \
                                         tey_xyy_xxxx_0, tey_xyy_xxxy_0, tey_xyy_xxxz_0, tey_xyy_xxyy_0, tey_xyy_xxyz_0, \
                                         tey_xz_xxyy_0, tey_xz_xxyy_1, tey_xz_xxyz_0, tey_xz_xxyz_1, tey_xz_xxzz_0, \
                                         tey_xz_xxzz_1, tey_xz_xyy_0, tey_xz_xyy_1, tey_xz_xyyy_0, tey_xz_xyyy_1, \
                                         tey_xz_xyyz_0, tey_xz_xyyz_1, tey_xz_xyz_0, tey_xz_xyz_1, tey_xz_xyzz_0, \
                                         tey_xz_xyzz_1, tey_xz_xzz_0, tey_xz_xzz_1, tey_xz_xzzz_0, tey_xz_xzzz_1, \
                                         tey_xz_yyy_0, tey_xz_yyy_1, tey_xz_yyyy_0, tey_xz_yyyy_1, tey_xz_yyyz_0, \
                                         tey_xz_yyyz_1, tey_xz_yyz_0, tey_xz_yyz_1, tey_xz_yyzz_0, tey_xz_yyzz_1, \
                                         tey_xz_yzz_0, tey_xz_yzz_1, tey_xz_yzzz_0, tey_xz_yzzz_1, tey_xz_zzz_0, \
                                         tey_xz_zzz_1, tey_xz_zzzz_0, tey_xz_zzzz_1, tey_yy_xxx_0, tey_yy_xxx_1, \
                                         tey_yy_xxxx_0, tey_yy_xxxx_1, tey_yy_xxxy_0, tey_yy_xxxy_1, tey_yy_xxxz_0, \
                                         tey_yy_xxxz_1, tey_yy_xxy_0, tey_yy_xxy_1, tey_yy_xxyy_0, tey_yy_xxyy_1, \
                                         tey_yy_xxyz_0, tey_yy_xxyz_1, tey_yy_xxz_0, tey_yy_xxz_1, tey_yy_xyy_0, \
                                         tey_yy_xyy_1, tey_yy_xyz_0, tey_yy_xyz_1, tey_z_xxyy_0, tey_z_xxyy_1, tey_z_xxyz_0, \
                                         tey_z_xxyz_1, tey_z_xxzz_0, tey_z_xxzz_1, tey_z_xyyy_0, tey_z_xyyy_1, tey_z_xyyz_0, \
                                         tey_z_xyyz_1, tey_z_xyzz_0, tey_z_xyzz_1, tey_z_xzzz_0, tey_z_xzzz_1, tey_z_yyyy_0, \
                                         tey_z_yyyy_1, tey_z_yyyz_0, tey_z_yyyz_1, tey_z_yyzz_0, tey_z_yyzz_1, tey_z_yzzz_0, \
                                         tey_z_yzzz_1, tey_z_zzzz_0, tey_z_zzzz_1, tez_xxz_xxyy_0, tez_xxz_xxyz_0, \
                                         tez_xxz_xxzz_0, tez_xxz_xyyy_0, tez_xxz_xyyz_0, tez_xxz_xyzz_0, tez_xxz_xzzz_0, \
                                         tez_xxz_yyyy_0, tez_xxz_yyyz_0, tez_xxz_yyzz_0, tez_xxz_yzzz_0, tez_xxz_zzzz_0, \
                                         tez_xyy_xxxx_0, tez_xyy_xxxy_0, tez_xyy_xxxz_0, tez_xyy_xxyy_0, tez_xyy_xxyz_0, \
                                         tez_xz_xxyy_0, tez_xz_xxyy_1, tez_xz_xxyz_0, tez_xz_xxyz_1, tez_xz_xxzz_0, \
                                         tez_xz_xxzz_1, tez_xz_xyy_0, tez_xz_xyy_1, tez_xz_xyyy_0, tez_xz_xyyy_1, \
                                         tez_xz_xyyz_0, tez_xz_xyyz_1, tez_xz_xyz_0, tez_xz_xyz_1, tez_xz_xyzz_0, \
                                         tez_xz_xyzz_1, tez_xz_xzz_0, tez_xz_xzz_1, tez_xz_xzzz_0, tez_xz_xzzz_1, \
                                         tez_xz_yyy_0, tez_xz_yyy_1, tez_xz_yyyy_0, tez_xz_yyyy_1, tez_xz_yyyz_0, \
                                         tez_xz_yyyz_1, tez_xz_yyz_0, tez_xz_yyz_1, tez_xz_yyzz_0, tez_xz_yyzz_1, \
                                         tez_xz_yzz_0, tez_xz_yzz_1, tez_xz_yzzz_0, tez_xz_yzzz_1, tez_xz_zzz_0, \
                                         tez_xz_zzz_1, tez_xz_zzzz_0, tez_xz_zzzz_1, tez_yy_xxx_0, tez_yy_xxx_1, \
                                         tez_yy_xxxx_0, tez_yy_xxxx_1, tez_yy_xxxy_0, tez_yy_xxxy_1, tez_yy_xxxz_0, \
                                         tez_yy_xxxz_1, tez_yy_xxy_0, tez_yy_xxy_1, tez_yy_xxyy_0, tez_yy_xxyy_1, \
                                         tez_yy_xxyz_0, tez_yy_xxyz_1, tez_yy_xxz_0, tez_yy_xxz_1, tez_yy_xyy_0, \
                                         tez_yy_xyy_1, tez_yy_xyz_0, tez_yy_xyz_1, tez_z_xxyy_0, tez_z_xxyy_1, tez_z_xxyz_0, \
                                         tez_z_xxyz_1, tez_z_xxzz_0, tez_z_xxzz_1, tez_z_xyyy_0, tez_z_xyyy_1, tez_z_xyyz_0, \
                                         tez_z_xyyz_1, tez_z_xyzz_0, tez_z_xyzz_1, tez_z_xzzz_0, tez_z_xzzz_1, tez_z_yyyy_0, \
                                         tez_z_yyyy_1, tez_z_yyyz_0, tez_z_yyyz_1, tez_z_yyzz_0, tez_z_yyzz_1, tez_z_yzzz_0, \
                                         tez_z_yzzz_1, tez_z_zzzz_0, tez_z_zzzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tey_xxz_xxyy_0[j] = pa_x[j] * tey_xz_xxyy_0[j] - pc_x[j] * tey_xz_xxyy_1[j] + 0.5 * fl1_fx * tey_z_xxyy_0[j] - 0.5 * fl1_fx * tey_z_xxyy_1[j] + fl1_fx * tey_xz_xyy_0[j] - fl1_fx * tey_xz_xyy_1[j];

                    tez_xxz_xxyy_0[j] = pa_x[j] * tez_xz_xxyy_0[j] - pc_x[j] * tez_xz_xxyy_1[j] + 0.5 * fl1_fx * tez_z_xxyy_0[j] - 0.5 * fl1_fx * tez_z_xxyy_1[j] + fl1_fx * tez_xz_xyy_0[j] - fl1_fx * tez_xz_xyy_1[j];

                    tex_xxz_xxyz_0[j] = pa_x[j] * tex_xz_xxyz_0[j] - pc_x[j] * tex_xz_xxyz_1[j] + 0.5 * fl1_fx * tex_z_xxyz_0[j] - 0.5 * fl1_fx * tex_z_xxyz_1[j] + fl1_fx * tex_xz_xyz_0[j] - fl1_fx * tex_xz_xyz_1[j] + ta_xz_xxyz_1[j];

                    tey_xxz_xxyz_0[j] = pa_x[j] * tey_xz_xxyz_0[j] - pc_x[j] * tey_xz_xxyz_1[j] + 0.5 * fl1_fx * tey_z_xxyz_0[j] - 0.5 * fl1_fx * tey_z_xxyz_1[j] + fl1_fx * tey_xz_xyz_0[j] - fl1_fx * tey_xz_xyz_1[j];

                    tez_xxz_xxyz_0[j] = pa_x[j] * tez_xz_xxyz_0[j] - pc_x[j] * tez_xz_xxyz_1[j] + 0.5 * fl1_fx * tez_z_xxyz_0[j] - 0.5 * fl1_fx * tez_z_xxyz_1[j] + fl1_fx * tez_xz_xyz_0[j] - fl1_fx * tez_xz_xyz_1[j];

                    tex_xxz_xxzz_0[j] = pa_x[j] * tex_xz_xxzz_0[j] - pc_x[j] * tex_xz_xxzz_1[j] + 0.5 * fl1_fx * tex_z_xxzz_0[j] - 0.5 * fl1_fx * tex_z_xxzz_1[j] + fl1_fx * tex_xz_xzz_0[j] - fl1_fx * tex_xz_xzz_1[j] + ta_xz_xxzz_1[j];

                    tey_xxz_xxzz_0[j] = pa_x[j] * tey_xz_xxzz_0[j] - pc_x[j] * tey_xz_xxzz_1[j] + 0.5 * fl1_fx * tey_z_xxzz_0[j] - 0.5 * fl1_fx * tey_z_xxzz_1[j] + fl1_fx * tey_xz_xzz_0[j] - fl1_fx * tey_xz_xzz_1[j];

                    tez_xxz_xxzz_0[j] = pa_x[j] * tez_xz_xxzz_0[j] - pc_x[j] * tez_xz_xxzz_1[j] + 0.5 * fl1_fx * tez_z_xxzz_0[j] - 0.5 * fl1_fx * tez_z_xxzz_1[j] + fl1_fx * tez_xz_xzz_0[j] - fl1_fx * tez_xz_xzz_1[j];

                    tex_xxz_xyyy_0[j] = pa_x[j] * tex_xz_xyyy_0[j] - pc_x[j] * tex_xz_xyyy_1[j] + 0.5 * fl1_fx * tex_z_xyyy_0[j] - 0.5 * fl1_fx * tex_z_xyyy_1[j] + 0.5 * fl1_fx * tex_xz_yyy_0[j] - 0.5 * fl1_fx * tex_xz_yyy_1[j] + ta_xz_xyyy_1[j];

                    tey_xxz_xyyy_0[j] = pa_x[j] * tey_xz_xyyy_0[j] - pc_x[j] * tey_xz_xyyy_1[j] + 0.5 * fl1_fx * tey_z_xyyy_0[j] - 0.5 * fl1_fx * tey_z_xyyy_1[j] + 0.5 * fl1_fx * tey_xz_yyy_0[j] - 0.5 * fl1_fx * tey_xz_yyy_1[j];

                    tez_xxz_xyyy_0[j] = pa_x[j] * tez_xz_xyyy_0[j] - pc_x[j] * tez_xz_xyyy_1[j] + 0.5 * fl1_fx * tez_z_xyyy_0[j] - 0.5 * fl1_fx * tez_z_xyyy_1[j] + 0.5 * fl1_fx * tez_xz_yyy_0[j] - 0.5 * fl1_fx * tez_xz_yyy_1[j];

                    tex_xxz_xyyz_0[j] = pa_x[j] * tex_xz_xyyz_0[j] - pc_x[j] * tex_xz_xyyz_1[j] + 0.5 * fl1_fx * tex_z_xyyz_0[j] - 0.5 * fl1_fx * tex_z_xyyz_1[j] + 0.5 * fl1_fx * tex_xz_yyz_0[j] - 0.5 * fl1_fx * tex_xz_yyz_1[j] + ta_xz_xyyz_1[j];

                    tey_xxz_xyyz_0[j] = pa_x[j] * tey_xz_xyyz_0[j] - pc_x[j] * tey_xz_xyyz_1[j] + 0.5 * fl1_fx * tey_z_xyyz_0[j] - 0.5 * fl1_fx * tey_z_xyyz_1[j] + 0.5 * fl1_fx * tey_xz_yyz_0[j] - 0.5 * fl1_fx * tey_xz_yyz_1[j];

                    tez_xxz_xyyz_0[j] = pa_x[j] * tez_xz_xyyz_0[j] - pc_x[j] * tez_xz_xyyz_1[j] + 0.5 * fl1_fx * tez_z_xyyz_0[j] - 0.5 * fl1_fx * tez_z_xyyz_1[j] + 0.5 * fl1_fx * tez_xz_yyz_0[j] - 0.5 * fl1_fx * tez_xz_yyz_1[j];

                    tex_xxz_xyzz_0[j] = pa_x[j] * tex_xz_xyzz_0[j] - pc_x[j] * tex_xz_xyzz_1[j] + 0.5 * fl1_fx * tex_z_xyzz_0[j] - 0.5 * fl1_fx * tex_z_xyzz_1[j] + 0.5 * fl1_fx * tex_xz_yzz_0[j] - 0.5 * fl1_fx * tex_xz_yzz_1[j] + ta_xz_xyzz_1[j];

                    tey_xxz_xyzz_0[j] = pa_x[j] * tey_xz_xyzz_0[j] - pc_x[j] * tey_xz_xyzz_1[j] + 0.5 * fl1_fx * tey_z_xyzz_0[j] - 0.5 * fl1_fx * tey_z_xyzz_1[j] + 0.5 * fl1_fx * tey_xz_yzz_0[j] - 0.5 * fl1_fx * tey_xz_yzz_1[j];

                    tez_xxz_xyzz_0[j] = pa_x[j] * tez_xz_xyzz_0[j] - pc_x[j] * tez_xz_xyzz_1[j] + 0.5 * fl1_fx * tez_z_xyzz_0[j] - 0.5 * fl1_fx * tez_z_xyzz_1[j] + 0.5 * fl1_fx * tez_xz_yzz_0[j] - 0.5 * fl1_fx * tez_xz_yzz_1[j];

                    tex_xxz_xzzz_0[j] = pa_x[j] * tex_xz_xzzz_0[j] - pc_x[j] * tex_xz_xzzz_1[j] + 0.5 * fl1_fx * tex_z_xzzz_0[j] - 0.5 * fl1_fx * tex_z_xzzz_1[j] + 0.5 * fl1_fx * tex_xz_zzz_0[j] - 0.5 * fl1_fx * tex_xz_zzz_1[j] + ta_xz_xzzz_1[j];

                    tey_xxz_xzzz_0[j] = pa_x[j] * tey_xz_xzzz_0[j] - pc_x[j] * tey_xz_xzzz_1[j] + 0.5 * fl1_fx * tey_z_xzzz_0[j] - 0.5 * fl1_fx * tey_z_xzzz_1[j] + 0.5 * fl1_fx * tey_xz_zzz_0[j] - 0.5 * fl1_fx * tey_xz_zzz_1[j];

                    tez_xxz_xzzz_0[j] = pa_x[j] * tez_xz_xzzz_0[j] - pc_x[j] * tez_xz_xzzz_1[j] + 0.5 * fl1_fx * tez_z_xzzz_0[j] - 0.5 * fl1_fx * tez_z_xzzz_1[j] + 0.5 * fl1_fx * tez_xz_zzz_0[j] - 0.5 * fl1_fx * tez_xz_zzz_1[j];

                    tex_xxz_yyyy_0[j] = pa_x[j] * tex_xz_yyyy_0[j] - pc_x[j] * tex_xz_yyyy_1[j] + 0.5 * fl1_fx * tex_z_yyyy_0[j] - 0.5 * fl1_fx * tex_z_yyyy_1[j] + ta_xz_yyyy_1[j];

                    tey_xxz_yyyy_0[j] = pa_x[j] * tey_xz_yyyy_0[j] - pc_x[j] * tey_xz_yyyy_1[j] + 0.5 * fl1_fx * tey_z_yyyy_0[j] - 0.5 * fl1_fx * tey_z_yyyy_1[j];

                    tez_xxz_yyyy_0[j] = pa_x[j] * tez_xz_yyyy_0[j] - pc_x[j] * tez_xz_yyyy_1[j] + 0.5 * fl1_fx * tez_z_yyyy_0[j] - 0.5 * fl1_fx * tez_z_yyyy_1[j];

                    tex_xxz_yyyz_0[j] = pa_x[j] * tex_xz_yyyz_0[j] - pc_x[j] * tex_xz_yyyz_1[j] + 0.5 * fl1_fx * tex_z_yyyz_0[j] - 0.5 * fl1_fx * tex_z_yyyz_1[j] + ta_xz_yyyz_1[j];

                    tey_xxz_yyyz_0[j] = pa_x[j] * tey_xz_yyyz_0[j] - pc_x[j] * tey_xz_yyyz_1[j] + 0.5 * fl1_fx * tey_z_yyyz_0[j] - 0.5 * fl1_fx * tey_z_yyyz_1[j];

                    tez_xxz_yyyz_0[j] = pa_x[j] * tez_xz_yyyz_0[j] - pc_x[j] * tez_xz_yyyz_1[j] + 0.5 * fl1_fx * tez_z_yyyz_0[j] - 0.5 * fl1_fx * tez_z_yyyz_1[j];

                    tex_xxz_yyzz_0[j] = pa_x[j] * tex_xz_yyzz_0[j] - pc_x[j] * tex_xz_yyzz_1[j] + 0.5 * fl1_fx * tex_z_yyzz_0[j] - 0.5 * fl1_fx * tex_z_yyzz_1[j] + ta_xz_yyzz_1[j];

                    tey_xxz_yyzz_0[j] = pa_x[j] * tey_xz_yyzz_0[j] - pc_x[j] * tey_xz_yyzz_1[j] + 0.5 * fl1_fx * tey_z_yyzz_0[j] - 0.5 * fl1_fx * tey_z_yyzz_1[j];

                    tez_xxz_yyzz_0[j] = pa_x[j] * tez_xz_yyzz_0[j] - pc_x[j] * tez_xz_yyzz_1[j] + 0.5 * fl1_fx * tez_z_yyzz_0[j] - 0.5 * fl1_fx * tez_z_yyzz_1[j];

                    tex_xxz_yzzz_0[j] = pa_x[j] * tex_xz_yzzz_0[j] - pc_x[j] * tex_xz_yzzz_1[j] + 0.5 * fl1_fx * tex_z_yzzz_0[j] - 0.5 * fl1_fx * tex_z_yzzz_1[j] + ta_xz_yzzz_1[j];

                    tey_xxz_yzzz_0[j] = pa_x[j] * tey_xz_yzzz_0[j] - pc_x[j] * tey_xz_yzzz_1[j] + 0.5 * fl1_fx * tey_z_yzzz_0[j] - 0.5 * fl1_fx * tey_z_yzzz_1[j];

                    tez_xxz_yzzz_0[j] = pa_x[j] * tez_xz_yzzz_0[j] - pc_x[j] * tez_xz_yzzz_1[j] + 0.5 * fl1_fx * tez_z_yzzz_0[j] - 0.5 * fl1_fx * tez_z_yzzz_1[j];

                    tex_xxz_zzzz_0[j] = pa_x[j] * tex_xz_zzzz_0[j] - pc_x[j] * tex_xz_zzzz_1[j] + 0.5 * fl1_fx * tex_z_zzzz_0[j] - 0.5 * fl1_fx * tex_z_zzzz_1[j] + ta_xz_zzzz_1[j];

                    tey_xxz_zzzz_0[j] = pa_x[j] * tey_xz_zzzz_0[j] - pc_x[j] * tey_xz_zzzz_1[j] + 0.5 * fl1_fx * tey_z_zzzz_0[j] - 0.5 * fl1_fx * tey_z_zzzz_1[j];

                    tez_xxz_zzzz_0[j] = pa_x[j] * tez_xz_zzzz_0[j] - pc_x[j] * tez_xz_zzzz_1[j] + 0.5 * fl1_fx * tez_z_zzzz_0[j] - 0.5 * fl1_fx * tez_z_zzzz_1[j];

                    tex_xyy_xxxx_0[j] = pa_x[j] * tex_yy_xxxx_0[j] - pc_x[j] * tex_yy_xxxx_1[j] + 2.0 * fl1_fx * tex_yy_xxx_0[j] - 2.0 * fl1_fx * tex_yy_xxx_1[j] + ta_yy_xxxx_1[j];

                    tey_xyy_xxxx_0[j] = pa_x[j] * tey_yy_xxxx_0[j] - pc_x[j] * tey_yy_xxxx_1[j] + 2.0 * fl1_fx * tey_yy_xxx_0[j] - 2.0 * fl1_fx * tey_yy_xxx_1[j];

                    tez_xyy_xxxx_0[j] = pa_x[j] * tez_yy_xxxx_0[j] - pc_x[j] * tez_yy_xxxx_1[j] + 2.0 * fl1_fx * tez_yy_xxx_0[j] - 2.0 * fl1_fx * tez_yy_xxx_1[j];

                    tex_xyy_xxxy_0[j] = pa_x[j] * tex_yy_xxxy_0[j] - pc_x[j] * tex_yy_xxxy_1[j] + 1.5 * fl1_fx * tex_yy_xxy_0[j] - 1.5 * fl1_fx * tex_yy_xxy_1[j] + ta_yy_xxxy_1[j];

                    tey_xyy_xxxy_0[j] = pa_x[j] * tey_yy_xxxy_0[j] - pc_x[j] * tey_yy_xxxy_1[j] + 1.5 * fl1_fx * tey_yy_xxy_0[j] - 1.5 * fl1_fx * tey_yy_xxy_1[j];

                    tez_xyy_xxxy_0[j] = pa_x[j] * tez_yy_xxxy_0[j] - pc_x[j] * tez_yy_xxxy_1[j] + 1.5 * fl1_fx * tez_yy_xxy_0[j] - 1.5 * fl1_fx * tez_yy_xxy_1[j];

                    tex_xyy_xxxz_0[j] = pa_x[j] * tex_yy_xxxz_0[j] - pc_x[j] * tex_yy_xxxz_1[j] + 1.5 * fl1_fx * tex_yy_xxz_0[j] - 1.5 * fl1_fx * tex_yy_xxz_1[j] + ta_yy_xxxz_1[j];

                    tey_xyy_xxxz_0[j] = pa_x[j] * tey_yy_xxxz_0[j] - pc_x[j] * tey_yy_xxxz_1[j] + 1.5 * fl1_fx * tey_yy_xxz_0[j] - 1.5 * fl1_fx * tey_yy_xxz_1[j];

                    tez_xyy_xxxz_0[j] = pa_x[j] * tez_yy_xxxz_0[j] - pc_x[j] * tez_yy_xxxz_1[j] + 1.5 * fl1_fx * tez_yy_xxz_0[j] - 1.5 * fl1_fx * tez_yy_xxz_1[j];

                    tex_xyy_xxyy_0[j] = pa_x[j] * tex_yy_xxyy_0[j] - pc_x[j] * tex_yy_xxyy_1[j] + fl1_fx * tex_yy_xyy_0[j] - fl1_fx * tex_yy_xyy_1[j] + ta_yy_xxyy_1[j];

                    tey_xyy_xxyy_0[j] = pa_x[j] * tey_yy_xxyy_0[j] - pc_x[j] * tey_yy_xxyy_1[j] + fl1_fx * tey_yy_xyy_0[j] - fl1_fx * tey_yy_xyy_1[j];

                    tez_xyy_xxyy_0[j] = pa_x[j] * tez_yy_xxyy_0[j] - pc_x[j] * tez_yy_xxyy_1[j] + fl1_fx * tez_yy_xyy_0[j] - fl1_fx * tez_yy_xyy_1[j];

                    tex_xyy_xxyz_0[j] = pa_x[j] * tex_yy_xxyz_0[j] - pc_x[j] * tex_yy_xxyz_1[j] + fl1_fx * tex_yy_xyz_0[j] - fl1_fx * tex_yy_xyz_1[j] + ta_yy_xxyz_1[j];

                    tey_xyy_xxyz_0[j] = pa_x[j] * tey_yy_xxyz_0[j] - pc_x[j] * tey_yy_xxyz_1[j] + fl1_fx * tey_yy_xyz_0[j] - fl1_fx * tey_yy_xyz_1[j];

                    tez_xyy_xxyz_0[j] = pa_x[j] * tez_yy_xxyz_0[j] - pc_x[j] * tez_yy_xxyz_1[j] + fl1_fx * tez_yy_xyz_0[j] - fl1_fx * tez_yy_xyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_150_200(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (150,200)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tex_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 50); 

                auto tey_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 50); 

                auto tez_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 50); 

                auto tex_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 51); 

                auto tey_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 51); 

                auto tez_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 51); 

                auto tex_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 52); 

                auto tey_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 52); 

                auto tez_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 52); 

                auto tex_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 53); 

                auto tey_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 53); 

                auto tez_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 53); 

                auto tex_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 54); 

                auto tey_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 54); 

                auto tez_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 54); 

                auto tex_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 55); 

                auto tey_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 55); 

                auto tez_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 55); 

                auto tex_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 56); 

                auto tey_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 56); 

                auto tez_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 56); 

                auto tex_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 57); 

                auto tey_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 57); 

                auto tez_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 57); 

                auto tex_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 58); 

                auto tey_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 58); 

                auto tez_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 58); 

                auto tex_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 59); 

                auto tey_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 59); 

                auto tez_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 59); 

                auto tex_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 60); 

                auto tey_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 60); 

                auto tez_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 60); 

                auto tex_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 61); 

                auto tey_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 61); 

                auto tez_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 61); 

                auto tex_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 62); 

                auto tey_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 62); 

                auto tez_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 62); 

                auto tex_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 63); 

                auto tey_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 63); 

                auto tez_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 63); 

                auto tex_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 64); 

                auto tey_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 64); 

                auto tez_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 64); 

                auto tex_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 65); 

                auto tey_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 65); 

                auto tez_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 65); 

                auto tex_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 66); 

                auto tey_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 66); 

                auto tex_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 50); 

                auto tey_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 50); 

                auto tez_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 50); 

                auto tex_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 51); 

                auto tey_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 51); 

                auto tez_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 51); 

                auto tex_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 52); 

                auto tey_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 52); 

                auto tez_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 52); 

                auto tex_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 53); 

                auto tey_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 53); 

                auto tez_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 53); 

                auto tex_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 54); 

                auto tey_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 54); 

                auto tez_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 54); 

                auto tex_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 55); 

                auto tey_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 55); 

                auto tez_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 55); 

                auto tex_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 56); 

                auto tey_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 56); 

                auto tez_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 56); 

                auto tex_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 57); 

                auto tey_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 57); 

                auto tez_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 57); 

                auto tex_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 58); 

                auto tey_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 58); 

                auto tez_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 58); 

                auto tex_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 59); 

                auto tey_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 59); 

                auto tez_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 59); 

                auto tex_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 60); 

                auto tey_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 60); 

                auto tez_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 60); 

                auto tex_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 61); 

                auto tey_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 61); 

                auto tez_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 61); 

                auto tex_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 62); 

                auto tey_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 62); 

                auto tez_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 62); 

                auto tex_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 63); 

                auto tey_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 63); 

                auto tez_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 63); 

                auto tex_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 64); 

                auto tey_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 64); 

                auto tez_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 64); 

                auto tex_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 65); 

                auto tey_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 65); 

                auto tez_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 65); 

                auto tex_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 66); 

                auto tey_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 66); 

                auto tex_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 35); 

                auto tey_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 35); 

                auto tez_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 35); 

                auto tex_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 36); 

                auto tey_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 36); 

                auto tez_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 36); 

                auto tex_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 37); 

                auto tey_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 37); 

                auto tez_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 37); 

                auto tex_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 38); 

                auto tey_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 38); 

                auto tez_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 38); 

                auto tex_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 39); 

                auto tey_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 39); 

                auto tez_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 39); 

                auto tex_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 40); 

                auto tey_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 40); 

                auto tez_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 40); 

                auto tex_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 41); 

                auto tey_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 41); 

                auto tez_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 41); 

                auto tex_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 42); 

                auto tey_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 42); 

                auto tez_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 42); 

                auto tex_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 43); 

                auto tey_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 43); 

                auto tez_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 43); 

                auto tex_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 44); 

                auto tey_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 44); 

                auto tez_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 44); 

                auto tex_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 45); 

                auto tey_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 45); 

                auto tez_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 45); 

                auto tex_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 46); 

                auto tey_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 46); 

                auto tex_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 35); 

                auto tey_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 35); 

                auto tez_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 35); 

                auto tex_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 36); 

                auto tey_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 36); 

                auto tez_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 36); 

                auto tex_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 37); 

                auto tey_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 37); 

                auto tez_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 37); 

                auto tex_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 38); 

                auto tey_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 38); 

                auto tez_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 38); 

                auto tex_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 39); 

                auto tey_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 39); 

                auto tez_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 39); 

                auto tex_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 40); 

                auto tey_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 40); 

                auto tez_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 40); 

                auto tex_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 41); 

                auto tey_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 41); 

                auto tez_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 41); 

                auto tex_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 42); 

                auto tey_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 42); 

                auto tez_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 42); 

                auto tex_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 43); 

                auto tey_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 43); 

                auto tez_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 43); 

                auto tex_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 44); 

                auto tey_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 44); 

                auto tez_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 44); 

                auto tex_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 45); 

                auto tey_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 45); 

                auto tez_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 45); 

                auto tex_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 46); 

                auto tey_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 46); 

                auto ta_yy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 50); 

                auto ta_yy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 51); 

                auto ta_yy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 52); 

                auto ta_yy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 53); 

                auto ta_yy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 54); 

                auto ta_yy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 55); 

                auto ta_yy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 56); 

                auto ta_yy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 57); 

                auto ta_yy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 58); 

                auto ta_yy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 59); 

                auto ta_yz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 60); 

                auto ta_yz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 61); 

                auto ta_yz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 62); 

                auto ta_yz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 63); 

                auto ta_yz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 64); 

                auto ta_yz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 65); 

                auto ta_yz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 66); 

                // set up pointers to integrals

                auto tex_xyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 50); 

                auto tey_xyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 50); 

                auto tez_xyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 50); 

                auto tex_xyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 51); 

                auto tey_xyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 51); 

                auto tez_xyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 51); 

                auto tex_xyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 52); 

                auto tey_xyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 52); 

                auto tez_xyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 52); 

                auto tex_xyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 53); 

                auto tey_xyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 53); 

                auto tez_xyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 53); 

                auto tex_xyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 54); 

                auto tey_xyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 54); 

                auto tez_xyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 54); 

                auto tex_xyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 55); 

                auto tey_xyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 55); 

                auto tez_xyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 55); 

                auto tex_xyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 56); 

                auto tey_xyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 56); 

                auto tez_xyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 56); 

                auto tex_xyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 57); 

                auto tey_xyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 57); 

                auto tez_xyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 57); 

                auto tex_xyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 58); 

                auto tey_xyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 58); 

                auto tez_xyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 58); 

                auto tex_xyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 59); 

                auto tey_xyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 59); 

                auto tez_xyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 59); 

                auto tex_xyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 60); 

                auto tey_xyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 60); 

                auto tez_xyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 60); 

                auto tex_xyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 61); 

                auto tey_xyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 61); 

                auto tez_xyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 61); 

                auto tex_xyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 62); 

                auto tey_xyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 62); 

                auto tez_xyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 62); 

                auto tex_xyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 63); 

                auto tey_xyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 63); 

                auto tez_xyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 63); 

                auto tex_xyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 64); 

                auto tey_xyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 64); 

                auto tez_xyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 64); 

                auto tex_xyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 65); 

                auto tey_xyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 65); 

                auto tez_xyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 65); 

                auto tex_xyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 66); 

                auto tey_xyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 66); 

                // Batch of Integrals (150,200)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_yy_xxzz_1, ta_yy_xyyy_1, ta_yy_xyyz_1, ta_yy_xyzz_1, \
                                         ta_yy_xzzz_1, ta_yy_yyyy_1, ta_yy_yyyz_1, ta_yy_yyzz_1, ta_yy_yzzz_1, ta_yy_zzzz_1, \
                                         ta_yz_xxxx_1, ta_yz_xxxy_1, ta_yz_xxxz_1, ta_yz_xxyy_1, ta_yz_xxyz_1, ta_yz_xxzz_1, \
                                         ta_yz_xyyy_1, tex_xyy_xxzz_0, tex_xyy_xyyy_0, tex_xyy_xyyz_0, tex_xyy_xyzz_0, \
                                         tex_xyy_xzzz_0, tex_xyy_yyyy_0, tex_xyy_yyyz_0, tex_xyy_yyzz_0, tex_xyy_yzzz_0, \
                                         tex_xyy_zzzz_0, tex_xyz_xxxx_0, tex_xyz_xxxy_0, tex_xyz_xxxz_0, tex_xyz_xxyy_0, \
                                         tex_xyz_xxyz_0, tex_xyz_xxzz_0, tex_xyz_xyyy_0, tex_yy_xxzz_0, tex_yy_xxzz_1, \
                                         tex_yy_xyyy_0, tex_yy_xyyy_1, tex_yy_xyyz_0, tex_yy_xyyz_1, tex_yy_xyzz_0, \
                                         tex_yy_xyzz_1, tex_yy_xzz_0, tex_yy_xzz_1, tex_yy_xzzz_0, tex_yy_xzzz_1, \
                                         tex_yy_yyy_0, tex_yy_yyy_1, tex_yy_yyyy_0, tex_yy_yyyy_1, tex_yy_yyyz_0, \
                                         tex_yy_yyyz_1, tex_yy_yyz_0, tex_yy_yyz_1, tex_yy_yyzz_0, tex_yy_yyzz_1, \
                                         tex_yy_yzz_0, tex_yy_yzz_1, tex_yy_yzzz_0, tex_yy_yzzz_1, tex_yy_zzz_0, \
                                         tex_yy_zzz_1, tex_yy_zzzz_0, tex_yy_zzzz_1, tex_yz_xxx_0, tex_yz_xxx_1, \
                                         tex_yz_xxxx_0, tex_yz_xxxx_1, tex_yz_xxxy_0, tex_yz_xxxy_1, tex_yz_xxxz_0, \
                                         tex_yz_xxxz_1, tex_yz_xxy_0, tex_yz_xxy_1, tex_yz_xxyy_0, tex_yz_xxyy_1, \
                                         tex_yz_xxyz_0, tex_yz_xxyz_1, tex_yz_xxz_0, tex_yz_xxz_1, tex_yz_xxzz_0, \
                                         tex_yz_xxzz_1, tex_yz_xyy_0, tex_yz_xyy_1, tex_yz_xyyy_0, tex_yz_xyyy_1, \
                                         tex_yz_xyz_0, tex_yz_xyz_1, tex_yz_xzz_0, tex_yz_xzz_1, tex_yz_yyy_0, tex_yz_yyy_1, \
                                         tey_xyy_xxzz_0, tey_xyy_xyyy_0, tey_xyy_xyyz_0, tey_xyy_xyzz_0, tey_xyy_xzzz_0, \
                                         tey_xyy_yyyy_0, tey_xyy_yyyz_0, tey_xyy_yyzz_0, tey_xyy_yzzz_0, tey_xyy_zzzz_0, \
                                         tey_xyz_xxxx_0, tey_xyz_xxxy_0, tey_xyz_xxxz_0, tey_xyz_xxyy_0, tey_xyz_xxyz_0, \
                                         tey_xyz_xxzz_0, tey_xyz_xyyy_0, tey_yy_xxzz_0, tey_yy_xxzz_1, tey_yy_xyyy_0, \
                                         tey_yy_xyyy_1, tey_yy_xyyz_0, tey_yy_xyyz_1, tey_yy_xyzz_0, tey_yy_xyzz_1, \
                                         tey_yy_xzz_0, tey_yy_xzz_1, tey_yy_xzzz_0, tey_yy_xzzz_1, tey_yy_yyy_0, \
                                         tey_yy_yyy_1, tey_yy_yyyy_0, tey_yy_yyyy_1, tey_yy_yyyz_0, tey_yy_yyyz_1, \
                                         tey_yy_yyz_0, tey_yy_yyz_1, tey_yy_yyzz_0, tey_yy_yyzz_1, tey_yy_yzz_0, \
                                         tey_yy_yzz_1, tey_yy_yzzz_0, tey_yy_yzzz_1, tey_yy_zzz_0, tey_yy_zzz_1, \
                                         tey_yy_zzzz_0, tey_yy_zzzz_1, tey_yz_xxx_0, tey_yz_xxx_1, tey_yz_xxxx_0, \
                                         tey_yz_xxxx_1, tey_yz_xxxy_0, tey_yz_xxxy_1, tey_yz_xxxz_0, tey_yz_xxxz_1, \
                                         tey_yz_xxy_0, tey_yz_xxy_1, tey_yz_xxyy_0, tey_yz_xxyy_1, tey_yz_xxyz_0, \
                                         tey_yz_xxyz_1, tey_yz_xxz_0, tey_yz_xxz_1, tey_yz_xxzz_0, tey_yz_xxzz_1, \
                                         tey_yz_xyy_0, tey_yz_xyy_1, tey_yz_xyyy_0, tey_yz_xyyy_1, tey_yz_xyz_0, \
                                         tey_yz_xyz_1, tey_yz_xzz_0, tey_yz_xzz_1, tey_yz_yyy_0, tey_yz_yyy_1, \
                                         tez_xyy_xxzz_0, tez_xyy_xyyy_0, tez_xyy_xyyz_0, tez_xyy_xyzz_0, tez_xyy_xzzz_0, \
                                         tez_xyy_yyyy_0, tez_xyy_yyyz_0, tez_xyy_yyzz_0, tez_xyy_yzzz_0, tez_xyy_zzzz_0, \
                                         tez_xyz_xxxx_0, tez_xyz_xxxy_0, tez_xyz_xxxz_0, tez_xyz_xxyy_0, tez_xyz_xxyz_0, \
                                         tez_xyz_xxzz_0, tez_yy_xxzz_0, tez_yy_xxzz_1, tez_yy_xyyy_0, tez_yy_xyyy_1, \
                                         tez_yy_xyyz_0, tez_yy_xyyz_1, tez_yy_xyzz_0, tez_yy_xyzz_1, tez_yy_xzz_0, \
                                         tez_yy_xzz_1, tez_yy_xzzz_0, tez_yy_xzzz_1, tez_yy_yyy_0, tez_yy_yyy_1, \
                                         tez_yy_yyyy_0, tez_yy_yyyy_1, tez_yy_yyyz_0, tez_yy_yyyz_1, tez_yy_yyz_0, \
                                         tez_yy_yyz_1, tez_yy_yyzz_0, tez_yy_yyzz_1, tez_yy_yzz_0, tez_yy_yzz_1, \
                                         tez_yy_yzzz_0, tez_yy_yzzz_1, tez_yy_zzz_0, tez_yy_zzz_1, tez_yy_zzzz_0, \
                                         tez_yy_zzzz_1, tez_yz_xxx_0, tez_yz_xxx_1, tez_yz_xxxx_0, tez_yz_xxxx_1, \
                                         tez_yz_xxxy_0, tez_yz_xxxy_1, tez_yz_xxxz_0, tez_yz_xxxz_1, tez_yz_xxy_0, \
                                         tez_yz_xxy_1, tez_yz_xxyy_0, tez_yz_xxyy_1, tez_yz_xxyz_0, tez_yz_xxyz_1, \
                                         tez_yz_xxz_0, tez_yz_xxz_1, tez_yz_xxzz_0, tez_yz_xxzz_1, tez_yz_xyy_0, \
                                         tez_yz_xyy_1, tez_yz_xyz_0, tez_yz_xyz_1, tez_yz_xzz_0, tez_yz_xzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xyy_xxzz_0[j] = pa_x[j] * tex_yy_xxzz_0[j] - pc_x[j] * tex_yy_xxzz_1[j] + fl1_fx * tex_yy_xzz_0[j] - fl1_fx * tex_yy_xzz_1[j] + ta_yy_xxzz_1[j];

                    tey_xyy_xxzz_0[j] = pa_x[j] * tey_yy_xxzz_0[j] - pc_x[j] * tey_yy_xxzz_1[j] + fl1_fx * tey_yy_xzz_0[j] - fl1_fx * tey_yy_xzz_1[j];

                    tez_xyy_xxzz_0[j] = pa_x[j] * tez_yy_xxzz_0[j] - pc_x[j] * tez_yy_xxzz_1[j] + fl1_fx * tez_yy_xzz_0[j] - fl1_fx * tez_yy_xzz_1[j];

                    tex_xyy_xyyy_0[j] = pa_x[j] * tex_yy_xyyy_0[j] - pc_x[j] * tex_yy_xyyy_1[j] + 0.5 * fl1_fx * tex_yy_yyy_0[j] - 0.5 * fl1_fx * tex_yy_yyy_1[j] + ta_yy_xyyy_1[j];

                    tey_xyy_xyyy_0[j] = pa_x[j] * tey_yy_xyyy_0[j] - pc_x[j] * tey_yy_xyyy_1[j] + 0.5 * fl1_fx * tey_yy_yyy_0[j] - 0.5 * fl1_fx * tey_yy_yyy_1[j];

                    tez_xyy_xyyy_0[j] = pa_x[j] * tez_yy_xyyy_0[j] - pc_x[j] * tez_yy_xyyy_1[j] + 0.5 * fl1_fx * tez_yy_yyy_0[j] - 0.5 * fl1_fx * tez_yy_yyy_1[j];

                    tex_xyy_xyyz_0[j] = pa_x[j] * tex_yy_xyyz_0[j] - pc_x[j] * tex_yy_xyyz_1[j] + 0.5 * fl1_fx * tex_yy_yyz_0[j] - 0.5 * fl1_fx * tex_yy_yyz_1[j] + ta_yy_xyyz_1[j];

                    tey_xyy_xyyz_0[j] = pa_x[j] * tey_yy_xyyz_0[j] - pc_x[j] * tey_yy_xyyz_1[j] + 0.5 * fl1_fx * tey_yy_yyz_0[j] - 0.5 * fl1_fx * tey_yy_yyz_1[j];

                    tez_xyy_xyyz_0[j] = pa_x[j] * tez_yy_xyyz_0[j] - pc_x[j] * tez_yy_xyyz_1[j] + 0.5 * fl1_fx * tez_yy_yyz_0[j] - 0.5 * fl1_fx * tez_yy_yyz_1[j];

                    tex_xyy_xyzz_0[j] = pa_x[j] * tex_yy_xyzz_0[j] - pc_x[j] * tex_yy_xyzz_1[j] + 0.5 * fl1_fx * tex_yy_yzz_0[j] - 0.5 * fl1_fx * tex_yy_yzz_1[j] + ta_yy_xyzz_1[j];

                    tey_xyy_xyzz_0[j] = pa_x[j] * tey_yy_xyzz_0[j] - pc_x[j] * tey_yy_xyzz_1[j] + 0.5 * fl1_fx * tey_yy_yzz_0[j] - 0.5 * fl1_fx * tey_yy_yzz_1[j];

                    tez_xyy_xyzz_0[j] = pa_x[j] * tez_yy_xyzz_0[j] - pc_x[j] * tez_yy_xyzz_1[j] + 0.5 * fl1_fx * tez_yy_yzz_0[j] - 0.5 * fl1_fx * tez_yy_yzz_1[j];

                    tex_xyy_xzzz_0[j] = pa_x[j] * tex_yy_xzzz_0[j] - pc_x[j] * tex_yy_xzzz_1[j] + 0.5 * fl1_fx * tex_yy_zzz_0[j] - 0.5 * fl1_fx * tex_yy_zzz_1[j] + ta_yy_xzzz_1[j];

                    tey_xyy_xzzz_0[j] = pa_x[j] * tey_yy_xzzz_0[j] - pc_x[j] * tey_yy_xzzz_1[j] + 0.5 * fl1_fx * tey_yy_zzz_0[j] - 0.5 * fl1_fx * tey_yy_zzz_1[j];

                    tez_xyy_xzzz_0[j] = pa_x[j] * tez_yy_xzzz_0[j] - pc_x[j] * tez_yy_xzzz_1[j] + 0.5 * fl1_fx * tez_yy_zzz_0[j] - 0.5 * fl1_fx * tez_yy_zzz_1[j];

                    tex_xyy_yyyy_0[j] = pa_x[j] * tex_yy_yyyy_0[j] - pc_x[j] * tex_yy_yyyy_1[j] + ta_yy_yyyy_1[j];

                    tey_xyy_yyyy_0[j] = pa_x[j] * tey_yy_yyyy_0[j] - pc_x[j] * tey_yy_yyyy_1[j];

                    tez_xyy_yyyy_0[j] = pa_x[j] * tez_yy_yyyy_0[j] - pc_x[j] * tez_yy_yyyy_1[j];

                    tex_xyy_yyyz_0[j] = pa_x[j] * tex_yy_yyyz_0[j] - pc_x[j] * tex_yy_yyyz_1[j] + ta_yy_yyyz_1[j];

                    tey_xyy_yyyz_0[j] = pa_x[j] * tey_yy_yyyz_0[j] - pc_x[j] * tey_yy_yyyz_1[j];

                    tez_xyy_yyyz_0[j] = pa_x[j] * tez_yy_yyyz_0[j] - pc_x[j] * tez_yy_yyyz_1[j];

                    tex_xyy_yyzz_0[j] = pa_x[j] * tex_yy_yyzz_0[j] - pc_x[j] * tex_yy_yyzz_1[j] + ta_yy_yyzz_1[j];

                    tey_xyy_yyzz_0[j] = pa_x[j] * tey_yy_yyzz_0[j] - pc_x[j] * tey_yy_yyzz_1[j];

                    tez_xyy_yyzz_0[j] = pa_x[j] * tez_yy_yyzz_0[j] - pc_x[j] * tez_yy_yyzz_1[j];

                    tex_xyy_yzzz_0[j] = pa_x[j] * tex_yy_yzzz_0[j] - pc_x[j] * tex_yy_yzzz_1[j] + ta_yy_yzzz_1[j];

                    tey_xyy_yzzz_0[j] = pa_x[j] * tey_yy_yzzz_0[j] - pc_x[j] * tey_yy_yzzz_1[j];

                    tez_xyy_yzzz_0[j] = pa_x[j] * tez_yy_yzzz_0[j] - pc_x[j] * tez_yy_yzzz_1[j];

                    tex_xyy_zzzz_0[j] = pa_x[j] * tex_yy_zzzz_0[j] - pc_x[j] * tex_yy_zzzz_1[j] + ta_yy_zzzz_1[j];

                    tey_xyy_zzzz_0[j] = pa_x[j] * tey_yy_zzzz_0[j] - pc_x[j] * tey_yy_zzzz_1[j];

                    tez_xyy_zzzz_0[j] = pa_x[j] * tez_yy_zzzz_0[j] - pc_x[j] * tez_yy_zzzz_1[j];

                    tex_xyz_xxxx_0[j] = pa_x[j] * tex_yz_xxxx_0[j] - pc_x[j] * tex_yz_xxxx_1[j] + 2.0 * fl1_fx * tex_yz_xxx_0[j] - 2.0 * fl1_fx * tex_yz_xxx_1[j] + ta_yz_xxxx_1[j];

                    tey_xyz_xxxx_0[j] = pa_x[j] * tey_yz_xxxx_0[j] - pc_x[j] * tey_yz_xxxx_1[j] + 2.0 * fl1_fx * tey_yz_xxx_0[j] - 2.0 * fl1_fx * tey_yz_xxx_1[j];

                    tez_xyz_xxxx_0[j] = pa_x[j] * tez_yz_xxxx_0[j] - pc_x[j] * tez_yz_xxxx_1[j] + 2.0 * fl1_fx * tez_yz_xxx_0[j] - 2.0 * fl1_fx * tez_yz_xxx_1[j];

                    tex_xyz_xxxy_0[j] = pa_x[j] * tex_yz_xxxy_0[j] - pc_x[j] * tex_yz_xxxy_1[j] + 1.5 * fl1_fx * tex_yz_xxy_0[j] - 1.5 * fl1_fx * tex_yz_xxy_1[j] + ta_yz_xxxy_1[j];

                    tey_xyz_xxxy_0[j] = pa_x[j] * tey_yz_xxxy_0[j] - pc_x[j] * tey_yz_xxxy_1[j] + 1.5 * fl1_fx * tey_yz_xxy_0[j] - 1.5 * fl1_fx * tey_yz_xxy_1[j];

                    tez_xyz_xxxy_0[j] = pa_x[j] * tez_yz_xxxy_0[j] - pc_x[j] * tez_yz_xxxy_1[j] + 1.5 * fl1_fx * tez_yz_xxy_0[j] - 1.5 * fl1_fx * tez_yz_xxy_1[j];

                    tex_xyz_xxxz_0[j] = pa_x[j] * tex_yz_xxxz_0[j] - pc_x[j] * tex_yz_xxxz_1[j] + 1.5 * fl1_fx * tex_yz_xxz_0[j] - 1.5 * fl1_fx * tex_yz_xxz_1[j] + ta_yz_xxxz_1[j];

                    tey_xyz_xxxz_0[j] = pa_x[j] * tey_yz_xxxz_0[j] - pc_x[j] * tey_yz_xxxz_1[j] + 1.5 * fl1_fx * tey_yz_xxz_0[j] - 1.5 * fl1_fx * tey_yz_xxz_1[j];

                    tez_xyz_xxxz_0[j] = pa_x[j] * tez_yz_xxxz_0[j] - pc_x[j] * tez_yz_xxxz_1[j] + 1.5 * fl1_fx * tez_yz_xxz_0[j] - 1.5 * fl1_fx * tez_yz_xxz_1[j];

                    tex_xyz_xxyy_0[j] = pa_x[j] * tex_yz_xxyy_0[j] - pc_x[j] * tex_yz_xxyy_1[j] + fl1_fx * tex_yz_xyy_0[j] - fl1_fx * tex_yz_xyy_1[j] + ta_yz_xxyy_1[j];

                    tey_xyz_xxyy_0[j] = pa_x[j] * tey_yz_xxyy_0[j] - pc_x[j] * tey_yz_xxyy_1[j] + fl1_fx * tey_yz_xyy_0[j] - fl1_fx * tey_yz_xyy_1[j];

                    tez_xyz_xxyy_0[j] = pa_x[j] * tez_yz_xxyy_0[j] - pc_x[j] * tez_yz_xxyy_1[j] + fl1_fx * tez_yz_xyy_0[j] - fl1_fx * tez_yz_xyy_1[j];

                    tex_xyz_xxyz_0[j] = pa_x[j] * tex_yz_xxyz_0[j] - pc_x[j] * tex_yz_xxyz_1[j] + fl1_fx * tex_yz_xyz_0[j] - fl1_fx * tex_yz_xyz_1[j] + ta_yz_xxyz_1[j];

                    tey_xyz_xxyz_0[j] = pa_x[j] * tey_yz_xxyz_0[j] - pc_x[j] * tey_yz_xxyz_1[j] + fl1_fx * tey_yz_xyz_0[j] - fl1_fx * tey_yz_xyz_1[j];

                    tez_xyz_xxyz_0[j] = pa_x[j] * tez_yz_xxyz_0[j] - pc_x[j] * tez_yz_xxyz_1[j] + fl1_fx * tez_yz_xyz_0[j] - fl1_fx * tez_yz_xyz_1[j];

                    tex_xyz_xxzz_0[j] = pa_x[j] * tex_yz_xxzz_0[j] - pc_x[j] * tex_yz_xxzz_1[j] + fl1_fx * tex_yz_xzz_0[j] - fl1_fx * tex_yz_xzz_1[j] + ta_yz_xxzz_1[j];

                    tey_xyz_xxzz_0[j] = pa_x[j] * tey_yz_xxzz_0[j] - pc_x[j] * tey_yz_xxzz_1[j] + fl1_fx * tey_yz_xzz_0[j] - fl1_fx * tey_yz_xzz_1[j];

                    tez_xyz_xxzz_0[j] = pa_x[j] * tez_yz_xxzz_0[j] - pc_x[j] * tez_yz_xxzz_1[j] + fl1_fx * tez_yz_xzz_0[j] - fl1_fx * tez_yz_xzz_1[j];

                    tex_xyz_xyyy_0[j] = pa_x[j] * tex_yz_xyyy_0[j] - pc_x[j] * tex_yz_xyyy_1[j] + 0.5 * fl1_fx * tex_yz_yyy_0[j] - 0.5 * fl1_fx * tex_yz_yyy_1[j] + ta_yz_xyyy_1[j];

                    tey_xyz_xyyy_0[j] = pa_x[j] * tey_yz_xyyy_0[j] - pc_x[j] * tey_yz_xyyy_1[j] + 0.5 * fl1_fx * tey_yz_yyy_0[j] - 0.5 * fl1_fx * tey_yz_yyy_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_200_250(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (200,250)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tez_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 66); 

                auto tex_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 67); 

                auto tey_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 67); 

                auto tez_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 67); 

                auto tex_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 68); 

                auto tey_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 68); 

                auto tez_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 68); 

                auto tex_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 69); 

                auto tey_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 69); 

                auto tez_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 69); 

                auto tex_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 70); 

                auto tey_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 70); 

                auto tez_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 70); 

                auto tex_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 71); 

                auto tey_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 71); 

                auto tez_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 71); 

                auto tex_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 72); 

                auto tey_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 72); 

                auto tez_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 72); 

                auto tex_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 73); 

                auto tey_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 73); 

                auto tez_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 73); 

                auto tex_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 74); 

                auto tey_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 74); 

                auto tez_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 74); 

                auto tex_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 75); 

                auto tey_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 75); 

                auto tez_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 75); 

                auto tex_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 76); 

                auto tey_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 76); 

                auto tez_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 76); 

                auto tex_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 77); 

                auto tey_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 77); 

                auto tez_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 77); 

                auto tex_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 78); 

                auto tey_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 78); 

                auto tez_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 78); 

                auto tex_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 79); 

                auto tey_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 79); 

                auto tez_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 79); 

                auto tex_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 80); 

                auto tey_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 80); 

                auto tez_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 80); 

                auto tex_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 81); 

                auto tey_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 81); 

                auto tez_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 81); 

                auto tex_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 82); 

                auto tey_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 82); 

                auto tez_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 82); 

                auto tex_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 83); 

                auto tez_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 66); 

                auto tex_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 67); 

                auto tey_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 67); 

                auto tez_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 67); 

                auto tex_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 68); 

                auto tey_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 68); 

                auto tez_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 68); 

                auto tex_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 69); 

                auto tey_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 69); 

                auto tez_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 69); 

                auto tex_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 70); 

                auto tey_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 70); 

                auto tez_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 70); 

                auto tex_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 71); 

                auto tey_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 71); 

                auto tez_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 71); 

                auto tex_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 72); 

                auto tey_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 72); 

                auto tez_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 72); 

                auto tex_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 73); 

                auto tey_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 73); 

                auto tez_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 73); 

                auto tex_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 74); 

                auto tey_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 74); 

                auto tez_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 74); 

                auto tex_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 75); 

                auto tey_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 75); 

                auto tez_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 75); 

                auto tex_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 76); 

                auto tey_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 76); 

                auto tez_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 76); 

                auto tex_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 77); 

                auto tey_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 77); 

                auto tez_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 77); 

                auto tex_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 78); 

                auto tey_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 78); 

                auto tez_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 78); 

                auto tex_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 79); 

                auto tey_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 79); 

                auto tez_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 79); 

                auto tex_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 80); 

                auto tey_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 80); 

                auto tez_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 80); 

                auto tex_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 81); 

                auto tey_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 81); 

                auto tez_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 81); 

                auto tex_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 82); 

                auto tey_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 82); 

                auto tez_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 82); 

                auto tex_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 83); 

                auto tez_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 46); 

                auto tex_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 47); 

                auto tey_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 47); 

                auto tez_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 47); 

                auto tex_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 48); 

                auto tey_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 48); 

                auto tez_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 48); 

                auto tex_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 49); 

                auto tey_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 49); 

                auto tez_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 49); 

                auto tex_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 50); 

                auto tey_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 50); 

                auto tez_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 50); 

                auto tex_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 51); 

                auto tey_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 51); 

                auto tez_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 51); 

                auto tex_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 52); 

                auto tey_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 52); 

                auto tez_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 52); 

                auto tex_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 53); 

                auto tey_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 53); 

                auto tez_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 53); 

                auto tex_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 54); 

                auto tey_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 54); 

                auto tez_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 54); 

                auto tex_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 55); 

                auto tey_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 55); 

                auto tez_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 55); 

                auto tex_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 56); 

                auto tey_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 56); 

                auto tez_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 56); 

                auto tex_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 57); 

                auto tey_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 57); 

                auto tez_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 57); 

                auto tex_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 58); 

                auto tez_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 46); 

                auto tex_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 47); 

                auto tey_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 47); 

                auto tez_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 47); 

                auto tex_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 48); 

                auto tey_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 48); 

                auto tez_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 48); 

                auto tex_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 49); 

                auto tey_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 49); 

                auto tez_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 49); 

                auto tex_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 50); 

                auto tey_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 50); 

                auto tez_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 50); 

                auto tex_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 51); 

                auto tey_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 51); 

                auto tez_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 51); 

                auto tex_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 52); 

                auto tey_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 52); 

                auto tez_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 52); 

                auto tex_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 53); 

                auto tey_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 53); 

                auto tez_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 53); 

                auto tex_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 54); 

                auto tey_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 54); 

                auto tez_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 54); 

                auto tex_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 55); 

                auto tey_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 55); 

                auto tez_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 55); 

                auto tex_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 56); 

                auto tey_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 56); 

                auto tez_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 56); 

                auto tex_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 57); 

                auto tey_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 57); 

                auto tez_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 57); 

                auto tex_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 58); 

                auto ta_yz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 67); 

                auto ta_yz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 68); 

                auto ta_yz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 69); 

                auto ta_yz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 70); 

                auto ta_yz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 71); 

                auto ta_yz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 72); 

                auto ta_yz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 73); 

                auto ta_yz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 74); 

                auto ta_zz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 75); 

                auto ta_zz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 76); 

                auto ta_zz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 77); 

                auto ta_zz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 78); 

                auto ta_zz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 79); 

                auto ta_zz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 80); 

                auto ta_zz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 81); 

                auto ta_zz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 82); 

                auto ta_zz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 83); 

                // set up pointers to integrals

                auto tez_xyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 66); 

                auto tex_xyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 67); 

                auto tey_xyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 67); 

                auto tez_xyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 67); 

                auto tex_xyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 68); 

                auto tey_xyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 68); 

                auto tez_xyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 68); 

                auto tex_xyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 69); 

                auto tey_xyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 69); 

                auto tez_xyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 69); 

                auto tex_xyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 70); 

                auto tey_xyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 70); 

                auto tez_xyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 70); 

                auto tex_xyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 71); 

                auto tey_xyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 71); 

                auto tez_xyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 71); 

                auto tex_xyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 72); 

                auto tey_xyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 72); 

                auto tez_xyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 72); 

                auto tex_xyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 73); 

                auto tey_xyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 73); 

                auto tez_xyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 73); 

                auto tex_xyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 74); 

                auto tey_xyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 74); 

                auto tez_xyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 74); 

                auto tex_xzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 75); 

                auto tey_xzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 75); 

                auto tez_xzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 75); 

                auto tex_xzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 76); 

                auto tey_xzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 76); 

                auto tez_xzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 76); 

                auto tex_xzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 77); 

                auto tey_xzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 77); 

                auto tez_xzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 77); 

                auto tex_xzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 78); 

                auto tey_xzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 78); 

                auto tez_xzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 78); 

                auto tex_xzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 79); 

                auto tey_xzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 79); 

                auto tez_xzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 79); 

                auto tex_xzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 80); 

                auto tey_xzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 80); 

                auto tez_xzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 80); 

                auto tex_xzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 81); 

                auto tey_xzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 81); 

                auto tez_xzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 81); 

                auto tex_xzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 82); 

                auto tey_xzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 82); 

                auto tez_xzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 82); 

                auto tex_xzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 83); 

                // Batch of Integrals (200,250)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_yz_xyyz_1, ta_yz_xyzz_1, ta_yz_xzzz_1, ta_yz_yyyy_1, \
                                         ta_yz_yyyz_1, ta_yz_yyzz_1, ta_yz_yzzz_1, ta_yz_zzzz_1, ta_zz_xxxx_1, ta_zz_xxxy_1, \
                                         ta_zz_xxxz_1, ta_zz_xxyy_1, ta_zz_xxyz_1, ta_zz_xxzz_1, ta_zz_xyyy_1, ta_zz_xyyz_1, \
                                         ta_zz_xyzz_1, tex_xyz_xyyz_0, tex_xyz_xyzz_0, tex_xyz_xzzz_0, tex_xyz_yyyy_0, \
                                         tex_xyz_yyyz_0, tex_xyz_yyzz_0, tex_xyz_yzzz_0, tex_xyz_zzzz_0, tex_xzz_xxxx_0, \
                                         tex_xzz_xxxy_0, tex_xzz_xxxz_0, tex_xzz_xxyy_0, tex_xzz_xxyz_0, tex_xzz_xxzz_0, \
                                         tex_xzz_xyyy_0, tex_xzz_xyyz_0, tex_xzz_xyzz_0, tex_yz_xyyz_0, tex_yz_xyyz_1, \
                                         tex_yz_xyzz_0, tex_yz_xyzz_1, tex_yz_xzzz_0, tex_yz_xzzz_1, tex_yz_yyyy_0, \
                                         tex_yz_yyyy_1, tex_yz_yyyz_0, tex_yz_yyyz_1, tex_yz_yyz_0, tex_yz_yyz_1, \
                                         tex_yz_yyzz_0, tex_yz_yyzz_1, tex_yz_yzz_0, tex_yz_yzz_1, tex_yz_yzzz_0, \
                                         tex_yz_yzzz_1, tex_yz_zzz_0, tex_yz_zzz_1, tex_yz_zzzz_0, tex_yz_zzzz_1, \
                                         tex_zz_xxx_0, tex_zz_xxx_1, tex_zz_xxxx_0, tex_zz_xxxx_1, tex_zz_xxxy_0, \
                                         tex_zz_xxxy_1, tex_zz_xxxz_0, tex_zz_xxxz_1, tex_zz_xxy_0, tex_zz_xxy_1, \
                                         tex_zz_xxyy_0, tex_zz_xxyy_1, tex_zz_xxyz_0, tex_zz_xxyz_1, tex_zz_xxz_0, \
                                         tex_zz_xxz_1, tex_zz_xxzz_0, tex_zz_xxzz_1, tex_zz_xyy_0, tex_zz_xyy_1, \
                                         tex_zz_xyyy_0, tex_zz_xyyy_1, tex_zz_xyyz_0, tex_zz_xyyz_1, tex_zz_xyz_0, \
                                         tex_zz_xyz_1, tex_zz_xyzz_0, tex_zz_xyzz_1, tex_zz_xzz_0, tex_zz_xzz_1, \
                                         tex_zz_yyy_0, tex_zz_yyy_1, tex_zz_yyz_0, tex_zz_yyz_1, tex_zz_yzz_0, tex_zz_yzz_1, \
                                         tey_xyz_xyyz_0, tey_xyz_xyzz_0, tey_xyz_xzzz_0, tey_xyz_yyyy_0, tey_xyz_yyyz_0, \
                                         tey_xyz_yyzz_0, tey_xyz_yzzz_0, tey_xyz_zzzz_0, tey_xzz_xxxx_0, tey_xzz_xxxy_0, \
                                         tey_xzz_xxxz_0, tey_xzz_xxyy_0, tey_xzz_xxyz_0, tey_xzz_xxzz_0, tey_xzz_xyyy_0, \
                                         tey_xzz_xyyz_0, tey_yz_xyyz_0, tey_yz_xyyz_1, tey_yz_xyzz_0, tey_yz_xyzz_1, \
                                         tey_yz_xzzz_0, tey_yz_xzzz_1, tey_yz_yyyy_0, tey_yz_yyyy_1, tey_yz_yyyz_0, \
                                         tey_yz_yyyz_1, tey_yz_yyz_0, tey_yz_yyz_1, tey_yz_yyzz_0, tey_yz_yyzz_1, \
                                         tey_yz_yzz_0, tey_yz_yzz_1, tey_yz_yzzz_0, tey_yz_yzzz_1, tey_yz_zzz_0, \
                                         tey_yz_zzz_1, tey_yz_zzzz_0, tey_yz_zzzz_1, tey_zz_xxx_0, tey_zz_xxx_1, \
                                         tey_zz_xxxx_0, tey_zz_xxxx_1, tey_zz_xxxy_0, tey_zz_xxxy_1, tey_zz_xxxz_0, \
                                         tey_zz_xxxz_1, tey_zz_xxy_0, tey_zz_xxy_1, tey_zz_xxyy_0, tey_zz_xxyy_1, \
                                         tey_zz_xxyz_0, tey_zz_xxyz_1, tey_zz_xxz_0, tey_zz_xxz_1, tey_zz_xxzz_0, \
                                         tey_zz_xxzz_1, tey_zz_xyy_0, tey_zz_xyy_1, tey_zz_xyyy_0, tey_zz_xyyy_1, \
                                         tey_zz_xyyz_0, tey_zz_xyyz_1, tey_zz_xyz_0, tey_zz_xyz_1, tey_zz_xzz_0, \
                                         tey_zz_xzz_1, tey_zz_yyy_0, tey_zz_yyy_1, tey_zz_yyz_0, tey_zz_yyz_1, \
                                         tez_xyz_xyyy_0, tez_xyz_xyyz_0, tez_xyz_xyzz_0, tez_xyz_xzzz_0, tez_xyz_yyyy_0, \
                                         tez_xyz_yyyz_0, tez_xyz_yyzz_0, tez_xyz_yzzz_0, tez_xyz_zzzz_0, tez_xzz_xxxx_0, \
                                         tez_xzz_xxxy_0, tez_xzz_xxxz_0, tez_xzz_xxyy_0, tez_xzz_xxyz_0, tez_xzz_xxzz_0, \
                                         tez_xzz_xyyy_0, tez_xzz_xyyz_0, tez_yz_xyyy_0, tez_yz_xyyy_1, tez_yz_xyyz_0, \
                                         tez_yz_xyyz_1, tez_yz_xyzz_0, tez_yz_xyzz_1, tez_yz_xzzz_0, tez_yz_xzzz_1, \
                                         tez_yz_yyy_0, tez_yz_yyy_1, tez_yz_yyyy_0, tez_yz_yyyy_1, tez_yz_yyyz_0, \
                                         tez_yz_yyyz_1, tez_yz_yyz_0, tez_yz_yyz_1, tez_yz_yyzz_0, tez_yz_yyzz_1, \
                                         tez_yz_yzz_0, tez_yz_yzz_1, tez_yz_yzzz_0, tez_yz_yzzz_1, tez_yz_zzz_0, \
                                         tez_yz_zzz_1, tez_yz_zzzz_0, tez_yz_zzzz_1, tez_zz_xxx_0, tez_zz_xxx_1, \
                                         tez_zz_xxxx_0, tez_zz_xxxx_1, tez_zz_xxxy_0, tez_zz_xxxy_1, tez_zz_xxxz_0, \
                                         tez_zz_xxxz_1, tez_zz_xxy_0, tez_zz_xxy_1, tez_zz_xxyy_0, tez_zz_xxyy_1, \
                                         tez_zz_xxyz_0, tez_zz_xxyz_1, tez_zz_xxz_0, tez_zz_xxz_1, tez_zz_xxzz_0, \
                                         tez_zz_xxzz_1, tez_zz_xyy_0, tez_zz_xyy_1, tez_zz_xyyy_0, tez_zz_xyyy_1, \
                                         tez_zz_xyyz_0, tez_zz_xyyz_1, tez_zz_xyz_0, tez_zz_xyz_1, tez_zz_xzz_0, \
                                         tez_zz_xzz_1, tez_zz_yyy_0, tez_zz_yyy_1, tez_zz_yyz_0, tez_zz_yyz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tez_xyz_xyyy_0[j] = pa_x[j] * tez_yz_xyyy_0[j] - pc_x[j] * tez_yz_xyyy_1[j] + 0.5 * fl1_fx * tez_yz_yyy_0[j] - 0.5 * fl1_fx * tez_yz_yyy_1[j];

                    tex_xyz_xyyz_0[j] = pa_x[j] * tex_yz_xyyz_0[j] - pc_x[j] * tex_yz_xyyz_1[j] + 0.5 * fl1_fx * tex_yz_yyz_0[j] - 0.5 * fl1_fx * tex_yz_yyz_1[j] + ta_yz_xyyz_1[j];

                    tey_xyz_xyyz_0[j] = pa_x[j] * tey_yz_xyyz_0[j] - pc_x[j] * tey_yz_xyyz_1[j] + 0.5 * fl1_fx * tey_yz_yyz_0[j] - 0.5 * fl1_fx * tey_yz_yyz_1[j];

                    tez_xyz_xyyz_0[j] = pa_x[j] * tez_yz_xyyz_0[j] - pc_x[j] * tez_yz_xyyz_1[j] + 0.5 * fl1_fx * tez_yz_yyz_0[j] - 0.5 * fl1_fx * tez_yz_yyz_1[j];

                    tex_xyz_xyzz_0[j] = pa_x[j] * tex_yz_xyzz_0[j] - pc_x[j] * tex_yz_xyzz_1[j] + 0.5 * fl1_fx * tex_yz_yzz_0[j] - 0.5 * fl1_fx * tex_yz_yzz_1[j] + ta_yz_xyzz_1[j];

                    tey_xyz_xyzz_0[j] = pa_x[j] * tey_yz_xyzz_0[j] - pc_x[j] * tey_yz_xyzz_1[j] + 0.5 * fl1_fx * tey_yz_yzz_0[j] - 0.5 * fl1_fx * tey_yz_yzz_1[j];

                    tez_xyz_xyzz_0[j] = pa_x[j] * tez_yz_xyzz_0[j] - pc_x[j] * tez_yz_xyzz_1[j] + 0.5 * fl1_fx * tez_yz_yzz_0[j] - 0.5 * fl1_fx * tez_yz_yzz_1[j];

                    tex_xyz_xzzz_0[j] = pa_x[j] * tex_yz_xzzz_0[j] - pc_x[j] * tex_yz_xzzz_1[j] + 0.5 * fl1_fx * tex_yz_zzz_0[j] - 0.5 * fl1_fx * tex_yz_zzz_1[j] + ta_yz_xzzz_1[j];

                    tey_xyz_xzzz_0[j] = pa_x[j] * tey_yz_xzzz_0[j] - pc_x[j] * tey_yz_xzzz_1[j] + 0.5 * fl1_fx * tey_yz_zzz_0[j] - 0.5 * fl1_fx * tey_yz_zzz_1[j];

                    tez_xyz_xzzz_0[j] = pa_x[j] * tez_yz_xzzz_0[j] - pc_x[j] * tez_yz_xzzz_1[j] + 0.5 * fl1_fx * tez_yz_zzz_0[j] - 0.5 * fl1_fx * tez_yz_zzz_1[j];

                    tex_xyz_yyyy_0[j] = pa_x[j] * tex_yz_yyyy_0[j] - pc_x[j] * tex_yz_yyyy_1[j] + ta_yz_yyyy_1[j];

                    tey_xyz_yyyy_0[j] = pa_x[j] * tey_yz_yyyy_0[j] - pc_x[j] * tey_yz_yyyy_1[j];

                    tez_xyz_yyyy_0[j] = pa_x[j] * tez_yz_yyyy_0[j] - pc_x[j] * tez_yz_yyyy_1[j];

                    tex_xyz_yyyz_0[j] = pa_x[j] * tex_yz_yyyz_0[j] - pc_x[j] * tex_yz_yyyz_1[j] + ta_yz_yyyz_1[j];

                    tey_xyz_yyyz_0[j] = pa_x[j] * tey_yz_yyyz_0[j] - pc_x[j] * tey_yz_yyyz_1[j];

                    tez_xyz_yyyz_0[j] = pa_x[j] * tez_yz_yyyz_0[j] - pc_x[j] * tez_yz_yyyz_1[j];

                    tex_xyz_yyzz_0[j] = pa_x[j] * tex_yz_yyzz_0[j] - pc_x[j] * tex_yz_yyzz_1[j] + ta_yz_yyzz_1[j];

                    tey_xyz_yyzz_0[j] = pa_x[j] * tey_yz_yyzz_0[j] - pc_x[j] * tey_yz_yyzz_1[j];

                    tez_xyz_yyzz_0[j] = pa_x[j] * tez_yz_yyzz_0[j] - pc_x[j] * tez_yz_yyzz_1[j];

                    tex_xyz_yzzz_0[j] = pa_x[j] * tex_yz_yzzz_0[j] - pc_x[j] * tex_yz_yzzz_1[j] + ta_yz_yzzz_1[j];

                    tey_xyz_yzzz_0[j] = pa_x[j] * tey_yz_yzzz_0[j] - pc_x[j] * tey_yz_yzzz_1[j];

                    tez_xyz_yzzz_0[j] = pa_x[j] * tez_yz_yzzz_0[j] - pc_x[j] * tez_yz_yzzz_1[j];

                    tex_xyz_zzzz_0[j] = pa_x[j] * tex_yz_zzzz_0[j] - pc_x[j] * tex_yz_zzzz_1[j] + ta_yz_zzzz_1[j];

                    tey_xyz_zzzz_0[j] = pa_x[j] * tey_yz_zzzz_0[j] - pc_x[j] * tey_yz_zzzz_1[j];

                    tez_xyz_zzzz_0[j] = pa_x[j] * tez_yz_zzzz_0[j] - pc_x[j] * tez_yz_zzzz_1[j];

                    tex_xzz_xxxx_0[j] = pa_x[j] * tex_zz_xxxx_0[j] - pc_x[j] * tex_zz_xxxx_1[j] + 2.0 * fl1_fx * tex_zz_xxx_0[j] - 2.0 * fl1_fx * tex_zz_xxx_1[j] + ta_zz_xxxx_1[j];

                    tey_xzz_xxxx_0[j] = pa_x[j] * tey_zz_xxxx_0[j] - pc_x[j] * tey_zz_xxxx_1[j] + 2.0 * fl1_fx * tey_zz_xxx_0[j] - 2.0 * fl1_fx * tey_zz_xxx_1[j];

                    tez_xzz_xxxx_0[j] = pa_x[j] * tez_zz_xxxx_0[j] - pc_x[j] * tez_zz_xxxx_1[j] + 2.0 * fl1_fx * tez_zz_xxx_0[j] - 2.0 * fl1_fx * tez_zz_xxx_1[j];

                    tex_xzz_xxxy_0[j] = pa_x[j] * tex_zz_xxxy_0[j] - pc_x[j] * tex_zz_xxxy_1[j] + 1.5 * fl1_fx * tex_zz_xxy_0[j] - 1.5 * fl1_fx * tex_zz_xxy_1[j] + ta_zz_xxxy_1[j];

                    tey_xzz_xxxy_0[j] = pa_x[j] * tey_zz_xxxy_0[j] - pc_x[j] * tey_zz_xxxy_1[j] + 1.5 * fl1_fx * tey_zz_xxy_0[j] - 1.5 * fl1_fx * tey_zz_xxy_1[j];

                    tez_xzz_xxxy_0[j] = pa_x[j] * tez_zz_xxxy_0[j] - pc_x[j] * tez_zz_xxxy_1[j] + 1.5 * fl1_fx * tez_zz_xxy_0[j] - 1.5 * fl1_fx * tez_zz_xxy_1[j];

                    tex_xzz_xxxz_0[j] = pa_x[j] * tex_zz_xxxz_0[j] - pc_x[j] * tex_zz_xxxz_1[j] + 1.5 * fl1_fx * tex_zz_xxz_0[j] - 1.5 * fl1_fx * tex_zz_xxz_1[j] + ta_zz_xxxz_1[j];

                    tey_xzz_xxxz_0[j] = pa_x[j] * tey_zz_xxxz_0[j] - pc_x[j] * tey_zz_xxxz_1[j] + 1.5 * fl1_fx * tey_zz_xxz_0[j] - 1.5 * fl1_fx * tey_zz_xxz_1[j];

                    tez_xzz_xxxz_0[j] = pa_x[j] * tez_zz_xxxz_0[j] - pc_x[j] * tez_zz_xxxz_1[j] + 1.5 * fl1_fx * tez_zz_xxz_0[j] - 1.5 * fl1_fx * tez_zz_xxz_1[j];

                    tex_xzz_xxyy_0[j] = pa_x[j] * tex_zz_xxyy_0[j] - pc_x[j] * tex_zz_xxyy_1[j] + fl1_fx * tex_zz_xyy_0[j] - fl1_fx * tex_zz_xyy_1[j] + ta_zz_xxyy_1[j];

                    tey_xzz_xxyy_0[j] = pa_x[j] * tey_zz_xxyy_0[j] - pc_x[j] * tey_zz_xxyy_1[j] + fl1_fx * tey_zz_xyy_0[j] - fl1_fx * tey_zz_xyy_1[j];

                    tez_xzz_xxyy_0[j] = pa_x[j] * tez_zz_xxyy_0[j] - pc_x[j] * tez_zz_xxyy_1[j] + fl1_fx * tez_zz_xyy_0[j] - fl1_fx * tez_zz_xyy_1[j];

                    tex_xzz_xxyz_0[j] = pa_x[j] * tex_zz_xxyz_0[j] - pc_x[j] * tex_zz_xxyz_1[j] + fl1_fx * tex_zz_xyz_0[j] - fl1_fx * tex_zz_xyz_1[j] + ta_zz_xxyz_1[j];

                    tey_xzz_xxyz_0[j] = pa_x[j] * tey_zz_xxyz_0[j] - pc_x[j] * tey_zz_xxyz_1[j] + fl1_fx * tey_zz_xyz_0[j] - fl1_fx * tey_zz_xyz_1[j];

                    tez_xzz_xxyz_0[j] = pa_x[j] * tez_zz_xxyz_0[j] - pc_x[j] * tez_zz_xxyz_1[j] + fl1_fx * tez_zz_xyz_0[j] - fl1_fx * tez_zz_xyz_1[j];

                    tex_xzz_xxzz_0[j] = pa_x[j] * tex_zz_xxzz_0[j] - pc_x[j] * tex_zz_xxzz_1[j] + fl1_fx * tex_zz_xzz_0[j] - fl1_fx * tex_zz_xzz_1[j] + ta_zz_xxzz_1[j];

                    tey_xzz_xxzz_0[j] = pa_x[j] * tey_zz_xxzz_0[j] - pc_x[j] * tey_zz_xxzz_1[j] + fl1_fx * tey_zz_xzz_0[j] - fl1_fx * tey_zz_xzz_1[j];

                    tez_xzz_xxzz_0[j] = pa_x[j] * tez_zz_xxzz_0[j] - pc_x[j] * tez_zz_xxzz_1[j] + fl1_fx * tez_zz_xzz_0[j] - fl1_fx * tez_zz_xzz_1[j];

                    tex_xzz_xyyy_0[j] = pa_x[j] * tex_zz_xyyy_0[j] - pc_x[j] * tex_zz_xyyy_1[j] + 0.5 * fl1_fx * tex_zz_yyy_0[j] - 0.5 * fl1_fx * tex_zz_yyy_1[j] + ta_zz_xyyy_1[j];

                    tey_xzz_xyyy_0[j] = pa_x[j] * tey_zz_xyyy_0[j] - pc_x[j] * tey_zz_xyyy_1[j] + 0.5 * fl1_fx * tey_zz_yyy_0[j] - 0.5 * fl1_fx * tey_zz_yyy_1[j];

                    tez_xzz_xyyy_0[j] = pa_x[j] * tez_zz_xyyy_0[j] - pc_x[j] * tez_zz_xyyy_1[j] + 0.5 * fl1_fx * tez_zz_yyy_0[j] - 0.5 * fl1_fx * tez_zz_yyy_1[j];

                    tex_xzz_xyyz_0[j] = pa_x[j] * tex_zz_xyyz_0[j] - pc_x[j] * tex_zz_xyyz_1[j] + 0.5 * fl1_fx * tex_zz_yyz_0[j] - 0.5 * fl1_fx * tex_zz_yyz_1[j] + ta_zz_xyyz_1[j];

                    tey_xzz_xyyz_0[j] = pa_x[j] * tey_zz_xyyz_0[j] - pc_x[j] * tey_zz_xyyz_1[j] + 0.5 * fl1_fx * tey_zz_yyz_0[j] - 0.5 * fl1_fx * tey_zz_yyz_1[j];

                    tez_xzz_xyyz_0[j] = pa_x[j] * tez_zz_xyyz_0[j] - pc_x[j] * tez_zz_xyyz_1[j] + 0.5 * fl1_fx * tez_zz_yyz_0[j] - 0.5 * fl1_fx * tez_zz_yyz_1[j];

                    tex_xzz_xyzz_0[j] = pa_x[j] * tex_zz_xyzz_0[j] - pc_x[j] * tex_zz_xyzz_1[j] + 0.5 * fl1_fx * tex_zz_yzz_0[j] - 0.5 * fl1_fx * tex_zz_yzz_1[j] + ta_zz_xyzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_250_300(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (250,300)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to tensors product of distances R(PA) = P - A

                auto pa_x = paDistances.data(3 * idx);

                auto pa_y = paDistances.data(3 * idx + 1);

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_x = pcDistances.data(3 * idx);

                auto pc_y = pcDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tex_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 45); 

                auto tey_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 45); 

                auto tez_yy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 45); 

                auto tex_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 46); 

                auto tey_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 46); 

                auto tez_yy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 46); 

                auto tex_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 47); 

                auto tey_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 47); 

                auto tez_yy_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 47); 

                auto tex_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 48); 

                auto tey_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 48); 

                auto tez_yy_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 48); 

                auto tex_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 49); 

                auto tey_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 49); 

                auto tez_yy_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 49); 

                auto tex_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 50); 

                auto tey_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 50); 

                auto tez_yy_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 50); 

                auto tex_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 51); 

                auto tey_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 51); 

                auto tez_yy_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 51); 

                auto tex_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 52); 

                auto tey_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 52); 

                auto tez_yy_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 52); 

                auto tex_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 53); 

                auto tey_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 53); 

                auto tez_yy_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 53); 

                auto tex_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 54); 

                auto tey_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 54); 

                auto tez_yy_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 54); 

                auto tey_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 83); 

                auto tez_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 83); 

                auto tex_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 84); 

                auto tey_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 84); 

                auto tez_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 84); 

                auto tex_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 85); 

                auto tey_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 85); 

                auto tez_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 85); 

                auto tex_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 86); 

                auto tey_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 86); 

                auto tez_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 86); 

                auto tex_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 87); 

                auto tey_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 87); 

                auto tez_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 87); 

                auto tex_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 88); 

                auto tey_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 88); 

                auto tez_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 88); 

                auto tex_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 89); 

                auto tey_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 89); 

                auto tez_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 89); 

                auto tex_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 45); 

                auto tey_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 45); 

                auto tez_yy_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 45); 

                auto tex_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 46); 

                auto tey_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 46); 

                auto tez_yy_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 46); 

                auto tex_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 47); 

                auto tey_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 47); 

                auto tez_yy_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 47); 

                auto tex_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 48); 

                auto tey_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 48); 

                auto tez_yy_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 48); 

                auto tex_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 49); 

                auto tey_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 49); 

                auto tez_yy_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 49); 

                auto tex_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 50); 

                auto tey_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 50); 

                auto tez_yy_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 50); 

                auto tex_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 51); 

                auto tey_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 51); 

                auto tez_yy_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 51); 

                auto tex_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 52); 

                auto tey_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 52); 

                auto tez_yy_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 52); 

                auto tex_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 53); 

                auto tey_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 53); 

                auto tez_yy_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 53); 

                auto tex_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 54); 

                auto tey_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 54); 

                auto tez_yy_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 54); 

                auto tey_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 83); 

                auto tez_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 83); 

                auto tex_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 84); 

                auto tey_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 84); 

                auto tez_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 84); 

                auto tex_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 85); 

                auto tey_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 85); 

                auto tez_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 85); 

                auto tex_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 86); 

                auto tey_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 86); 

                auto tez_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 86); 

                auto tex_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 87); 

                auto tey_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 87); 

                auto tez_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 87); 

                auto tex_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 88); 

                auto tey_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 88); 

                auto tez_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 88); 

                auto tex_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 89); 

                auto tey_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 89); 

                auto tez_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 89); 

                auto tex_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 15); 

                auto tey_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 15); 

                auto tez_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 15); 

                auto tex_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 16); 

                auto tey_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 16); 

                auto tez_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 16); 

                auto tex_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 17); 

                auto tey_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 17); 

                auto tez_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 17); 

                auto tex_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 18); 

                auto tey_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 18); 

                auto tez_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 18); 

                auto tex_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 19); 

                auto tey_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 19); 

                auto tez_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 19); 

                auto tex_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 20); 

                auto tey_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 20); 

                auto tez_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 20); 

                auto tex_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 21); 

                auto tey_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 21); 

                auto tez_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 21); 

                auto tex_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 22); 

                auto tey_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 22); 

                auto tez_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 22); 

                auto tex_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 23); 

                auto tey_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 23); 

                auto tez_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 23); 

                auto tex_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 24); 

                auto tey_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 24); 

                auto tez_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 24); 

                auto tex_y_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 15); 

                auto tey_y_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 15); 

                auto tez_y_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 15); 

                auto tex_y_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 16); 

                auto tey_y_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 16); 

                auto tez_y_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 16); 

                auto tex_y_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 17); 

                auto tey_y_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 17); 

                auto tez_y_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 17); 

                auto tex_y_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 18); 

                auto tey_y_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 18); 

                auto tez_y_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 18); 

                auto tex_y_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 19); 

                auto tey_y_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 19); 

                auto tez_y_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 19); 

                auto tex_y_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 20); 

                auto tey_y_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 20); 

                auto tez_y_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 20); 

                auto tex_y_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 21); 

                auto tey_y_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 21); 

                auto tez_y_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 21); 

                auto tex_y_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 22); 

                auto tey_y_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 22); 

                auto tez_y_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 22); 

                auto tex_y_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 23); 

                auto tey_y_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 23); 

                auto tez_y_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 23); 

                auto tex_y_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 24); 

                auto tey_y_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 24); 

                auto tez_y_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 24); 

                auto tex_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 30); 

                auto tey_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 30); 

                auto tez_yy_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 30); 

                auto tex_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 31); 

                auto tey_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 31); 

                auto tez_yy_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 31); 

                auto tex_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 32); 

                auto tey_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 32); 

                auto tez_yy_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 32); 

                auto tex_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 33); 

                auto tey_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 33); 

                auto tez_yy_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 33); 

                auto tex_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 34); 

                auto tey_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 34); 

                auto tez_yy_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 34); 

                auto tex_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 35); 

                auto tey_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 35); 

                auto tez_yy_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 35); 

                auto tey_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 58); 

                auto tez_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 58); 

                auto tex_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 59); 

                auto tey_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 59); 

                auto tez_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 59); 

                auto tex_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 30); 

                auto tey_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 30); 

                auto tez_yy_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 30); 

                auto tex_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 31); 

                auto tey_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 31); 

                auto tez_yy_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 31); 

                auto tex_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 32); 

                auto tey_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 32); 

                auto tez_yy_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 32); 

                auto tex_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 33); 

                auto tey_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 33); 

                auto tez_yy_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 33); 

                auto tex_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 34); 

                auto tey_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 34); 

                auto tez_yy_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 34); 

                auto tex_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 35); 

                auto tey_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 35); 

                auto tez_yy_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 35); 

                auto tey_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 58); 

                auto tez_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 58); 

                auto tex_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 59); 

                auto tey_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 59); 

                auto tez_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 59); 

                auto ta_yy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 45); 

                auto ta_yy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 46); 

                auto ta_yy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 47); 

                auto ta_yy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 48); 

                auto ta_yy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 49); 

                auto ta_yy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 50); 

                auto ta_yy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 51); 

                auto ta_yy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 52); 

                auto ta_yy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 53); 

                auto ta_yy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 54); 

                auto ta_zz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 84); 

                auto ta_zz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 85); 

                auto ta_zz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 86); 

                auto ta_zz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 87); 

                auto ta_zz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 88); 

                auto ta_zz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 89); 

                // set up pointers to integrals

                auto tey_xzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 83); 

                auto tez_xzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 83); 

                auto tex_xzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 84); 

                auto tey_xzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 84); 

                auto tez_xzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 84); 

                auto tex_xzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 85); 

                auto tey_xzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 85); 

                auto tez_xzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 85); 

                auto tex_xzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 86); 

                auto tey_xzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 86); 

                auto tez_xzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 86); 

                auto tex_xzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 87); 

                auto tey_xzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 87); 

                auto tez_xzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 87); 

                auto tex_xzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 88); 

                auto tey_xzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 88); 

                auto tez_xzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 88); 

                auto tex_xzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 89); 

                auto tey_xzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 89); 

                auto tez_xzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 89); 

                auto tex_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 90); 

                auto tey_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 90); 

                auto tez_yyy_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 90); 

                auto tex_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 91); 

                auto tey_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 91); 

                auto tez_yyy_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 91); 

                auto tex_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 92); 

                auto tey_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 92); 

                auto tez_yyy_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 92); 

                auto tex_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 93); 

                auto tey_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 93); 

                auto tez_yyy_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 93); 

                auto tex_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 94); 

                auto tey_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 94); 

                auto tez_yyy_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 94); 

                auto tex_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 95); 

                auto tey_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 95); 

                auto tez_yyy_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 95); 

                auto tex_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 96); 

                auto tey_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 96); 

                auto tez_yyy_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 96); 

                auto tex_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 97); 

                auto tey_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 97); 

                auto tez_yyy_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 97); 

                auto tex_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 98); 

                auto tey_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 98); 

                auto tez_yyy_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 98); 

                auto tex_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 99); 

                auto tey_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 99); 

                auto tez_yyy_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 99); 

                // Batch of Integrals (250,300)

                #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_yy_xxxx_1, ta_yy_xxxy_1, ta_yy_xxxz_1, \
                                         ta_yy_xxyy_1, ta_yy_xxyz_1, ta_yy_xxzz_1, ta_yy_xyyy_1, ta_yy_xyyz_1, ta_yy_xyzz_1, \
                                         ta_yy_xzzz_1, ta_zz_xzzz_1, ta_zz_yyyy_1, ta_zz_yyyz_1, ta_zz_yyzz_1, ta_zz_yzzz_1, \
                                         ta_zz_zzzz_1, tex_xzz_xzzz_0, tex_xzz_yyyy_0, tex_xzz_yyyz_0, tex_xzz_yyzz_0, \
                                         tex_xzz_yzzz_0, tex_xzz_zzzz_0, tex_y_xxxx_0, tex_y_xxxx_1, tex_y_xxxy_0, \
                                         tex_y_xxxy_1, tex_y_xxxz_0, tex_y_xxxz_1, tex_y_xxyy_0, tex_y_xxyy_1, tex_y_xxyz_0, \
                                         tex_y_xxyz_1, tex_y_xxzz_0, tex_y_xxzz_1, tex_y_xyyy_0, tex_y_xyyy_1, tex_y_xyyz_0, \
                                         tex_y_xyyz_1, tex_y_xyzz_0, tex_y_xyzz_1, tex_y_xzzz_0, tex_y_xzzz_1, tex_yy_xxx_0, \
                                         tex_yy_xxx_1, tex_yy_xxxx_0, tex_yy_xxxx_1, tex_yy_xxxy_0, tex_yy_xxxy_1, \
                                         tex_yy_xxxz_0, tex_yy_xxxz_1, tex_yy_xxy_0, tex_yy_xxy_1, tex_yy_xxyy_0, \
                                         tex_yy_xxyy_1, tex_yy_xxyz_0, tex_yy_xxyz_1, tex_yy_xxz_0, tex_yy_xxz_1, \
                                         tex_yy_xxzz_0, tex_yy_xxzz_1, tex_yy_xyy_0, tex_yy_xyy_1, tex_yy_xyyy_0, \
                                         tex_yy_xyyy_1, tex_yy_xyyz_0, tex_yy_xyyz_1, tex_yy_xyz_0, tex_yy_xyz_1, \
                                         tex_yy_xyzz_0, tex_yy_xyzz_1, tex_yy_xzz_0, tex_yy_xzz_1, tex_yy_xzzz_0, \
                                         tex_yy_xzzz_1, tex_yyy_xxxx_0, tex_yyy_xxxy_0, tex_yyy_xxxz_0, tex_yyy_xxyy_0, \
                                         tex_yyy_xxyz_0, tex_yyy_xxzz_0, tex_yyy_xyyy_0, tex_yyy_xyyz_0, tex_yyy_xyzz_0, \
                                         tex_yyy_xzzz_0, tex_zz_xzzz_0, tex_zz_xzzz_1, tex_zz_yyyy_0, tex_zz_yyyy_1, \
                                         tex_zz_yyyz_0, tex_zz_yyyz_1, tex_zz_yyzz_0, tex_zz_yyzz_1, tex_zz_yzzz_0, \
                                         tex_zz_yzzz_1, tex_zz_zzz_0, tex_zz_zzz_1, tex_zz_zzzz_0, tex_zz_zzzz_1, \
                                         tey_xzz_xyzz_0, tey_xzz_xzzz_0, tey_xzz_yyyy_0, tey_xzz_yyyz_0, tey_xzz_yyzz_0, \
                                         tey_xzz_yzzz_0, tey_xzz_zzzz_0, tey_y_xxxx_0, tey_y_xxxx_1, tey_y_xxxy_0, \
                                         tey_y_xxxy_1, tey_y_xxxz_0, tey_y_xxxz_1, tey_y_xxyy_0, tey_y_xxyy_1, tey_y_xxyz_0, \
                                         tey_y_xxyz_1, tey_y_xxzz_0, tey_y_xxzz_1, tey_y_xyyy_0, tey_y_xyyy_1, tey_y_xyyz_0, \
                                         tey_y_xyyz_1, tey_y_xyzz_0, tey_y_xyzz_1, tey_y_xzzz_0, tey_y_xzzz_1, tey_yy_xxx_0, \
                                         tey_yy_xxx_1, tey_yy_xxxx_0, tey_yy_xxxx_1, tey_yy_xxxy_0, tey_yy_xxxy_1, \
                                         tey_yy_xxxz_0, tey_yy_xxxz_1, tey_yy_xxy_0, tey_yy_xxy_1, tey_yy_xxyy_0, \
                                         tey_yy_xxyy_1, tey_yy_xxyz_0, tey_yy_xxyz_1, tey_yy_xxz_0, tey_yy_xxz_1, \
                                         tey_yy_xxzz_0, tey_yy_xxzz_1, tey_yy_xyy_0, tey_yy_xyy_1, tey_yy_xyyy_0, \
                                         tey_yy_xyyy_1, tey_yy_xyyz_0, tey_yy_xyyz_1, tey_yy_xyz_0, tey_yy_xyz_1, \
                                         tey_yy_xyzz_0, tey_yy_xyzz_1, tey_yy_xzz_0, tey_yy_xzz_1, tey_yy_xzzz_0, \
                                         tey_yy_xzzz_1, tey_yyy_xxxx_0, tey_yyy_xxxy_0, tey_yyy_xxxz_0, tey_yyy_xxyy_0, \
                                         tey_yyy_xxyz_0, tey_yyy_xxzz_0, tey_yyy_xyyy_0, tey_yyy_xyyz_0, tey_yyy_xyzz_0, \
                                         tey_yyy_xzzz_0, tey_zz_xyzz_0, tey_zz_xyzz_1, tey_zz_xzzz_0, tey_zz_xzzz_1, \
                                         tey_zz_yyyy_0, tey_zz_yyyy_1, tey_zz_yyyz_0, tey_zz_yyyz_1, tey_zz_yyzz_0, \
                                         tey_zz_yyzz_1, tey_zz_yzz_0, tey_zz_yzz_1, tey_zz_yzzz_0, tey_zz_yzzz_1, \
                                         tey_zz_zzz_0, tey_zz_zzz_1, tey_zz_zzzz_0, tey_zz_zzzz_1, tez_xzz_xyzz_0, \
                                         tez_xzz_xzzz_0, tez_xzz_yyyy_0, tez_xzz_yyyz_0, tez_xzz_yyzz_0, tez_xzz_yzzz_0, \
                                         tez_xzz_zzzz_0, tez_y_xxxx_0, tez_y_xxxx_1, tez_y_xxxy_0, tez_y_xxxy_1, tez_y_xxxz_0, \
                                         tez_y_xxxz_1, tez_y_xxyy_0, tez_y_xxyy_1, tez_y_xxyz_0, tez_y_xxyz_1, tez_y_xxzz_0, \
                                         tez_y_xxzz_1, tez_y_xyyy_0, tez_y_xyyy_1, tez_y_xyyz_0, tez_y_xyyz_1, tez_y_xyzz_0, \
                                         tez_y_xyzz_1, tez_y_xzzz_0, tez_y_xzzz_1, tez_yy_xxx_0, tez_yy_xxx_1, \
                                         tez_yy_xxxx_0, tez_yy_xxxx_1, tez_yy_xxxy_0, tez_yy_xxxy_1, tez_yy_xxxz_0, \
                                         tez_yy_xxxz_1, tez_yy_xxy_0, tez_yy_xxy_1, tez_yy_xxyy_0, tez_yy_xxyy_1, \
                                         tez_yy_xxyz_0, tez_yy_xxyz_1, tez_yy_xxz_0, tez_yy_xxz_1, tez_yy_xxzz_0, \
                                         tez_yy_xxzz_1, tez_yy_xyy_0, tez_yy_xyy_1, tez_yy_xyyy_0, tez_yy_xyyy_1, \
                                         tez_yy_xyyz_0, tez_yy_xyyz_1, tez_yy_xyz_0, tez_yy_xyz_1, tez_yy_xyzz_0, \
                                         tez_yy_xyzz_1, tez_yy_xzz_0, tez_yy_xzz_1, tez_yy_xzzz_0, tez_yy_xzzz_1, \
                                         tez_yyy_xxxx_0, tez_yyy_xxxy_0, tez_yyy_xxxz_0, tez_yyy_xxyy_0, tez_yyy_xxyz_0, \
                                         tez_yyy_xxzz_0, tez_yyy_xyyy_0, tez_yyy_xyyz_0, tez_yyy_xyzz_0, tez_yyy_xzzz_0, \
                                         tez_zz_xyzz_0, tez_zz_xyzz_1, tez_zz_xzzz_0, tez_zz_xzzz_1, tez_zz_yyyy_0, \
                                         tez_zz_yyyy_1, tez_zz_yyyz_0, tez_zz_yyyz_1, tez_zz_yyzz_0, tez_zz_yyzz_1, \
                                         tez_zz_yzz_0, tez_zz_yzz_1, tez_zz_yzzz_0, tez_zz_yzzz_1, tez_zz_zzz_0, \
                                         tez_zz_zzz_1, tez_zz_zzzz_0, tez_zz_zzzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tey_xzz_xyzz_0[j] = pa_x[j] * tey_zz_xyzz_0[j] - pc_x[j] * tey_zz_xyzz_1[j] + 0.5 * fl1_fx * tey_zz_yzz_0[j] - 0.5 * fl1_fx * tey_zz_yzz_1[j];

                    tez_xzz_xyzz_0[j] = pa_x[j] * tez_zz_xyzz_0[j] - pc_x[j] * tez_zz_xyzz_1[j] + 0.5 * fl1_fx * tez_zz_yzz_0[j] - 0.5 * fl1_fx * tez_zz_yzz_1[j];

                    tex_xzz_xzzz_0[j] = pa_x[j] * tex_zz_xzzz_0[j] - pc_x[j] * tex_zz_xzzz_1[j] + 0.5 * fl1_fx * tex_zz_zzz_0[j] - 0.5 * fl1_fx * tex_zz_zzz_1[j] + ta_zz_xzzz_1[j];

                    tey_xzz_xzzz_0[j] = pa_x[j] * tey_zz_xzzz_0[j] - pc_x[j] * tey_zz_xzzz_1[j] + 0.5 * fl1_fx * tey_zz_zzz_0[j] - 0.5 * fl1_fx * tey_zz_zzz_1[j];

                    tez_xzz_xzzz_0[j] = pa_x[j] * tez_zz_xzzz_0[j] - pc_x[j] * tez_zz_xzzz_1[j] + 0.5 * fl1_fx * tez_zz_zzz_0[j] - 0.5 * fl1_fx * tez_zz_zzz_1[j];

                    tex_xzz_yyyy_0[j] = pa_x[j] * tex_zz_yyyy_0[j] - pc_x[j] * tex_zz_yyyy_1[j] + ta_zz_yyyy_1[j];

                    tey_xzz_yyyy_0[j] = pa_x[j] * tey_zz_yyyy_0[j] - pc_x[j] * tey_zz_yyyy_1[j];

                    tez_xzz_yyyy_0[j] = pa_x[j] * tez_zz_yyyy_0[j] - pc_x[j] * tez_zz_yyyy_1[j];

                    tex_xzz_yyyz_0[j] = pa_x[j] * tex_zz_yyyz_0[j] - pc_x[j] * tex_zz_yyyz_1[j] + ta_zz_yyyz_1[j];

                    tey_xzz_yyyz_0[j] = pa_x[j] * tey_zz_yyyz_0[j] - pc_x[j] * tey_zz_yyyz_1[j];

                    tez_xzz_yyyz_0[j] = pa_x[j] * tez_zz_yyyz_0[j] - pc_x[j] * tez_zz_yyyz_1[j];

                    tex_xzz_yyzz_0[j] = pa_x[j] * tex_zz_yyzz_0[j] - pc_x[j] * tex_zz_yyzz_1[j] + ta_zz_yyzz_1[j];

                    tey_xzz_yyzz_0[j] = pa_x[j] * tey_zz_yyzz_0[j] - pc_x[j] * tey_zz_yyzz_1[j];

                    tez_xzz_yyzz_0[j] = pa_x[j] * tez_zz_yyzz_0[j] - pc_x[j] * tez_zz_yyzz_1[j];

                    tex_xzz_yzzz_0[j] = pa_x[j] * tex_zz_yzzz_0[j] - pc_x[j] * tex_zz_yzzz_1[j] + ta_zz_yzzz_1[j];

                    tey_xzz_yzzz_0[j] = pa_x[j] * tey_zz_yzzz_0[j] - pc_x[j] * tey_zz_yzzz_1[j];

                    tez_xzz_yzzz_0[j] = pa_x[j] * tez_zz_yzzz_0[j] - pc_x[j] * tez_zz_yzzz_1[j];

                    tex_xzz_zzzz_0[j] = pa_x[j] * tex_zz_zzzz_0[j] - pc_x[j] * tex_zz_zzzz_1[j] + ta_zz_zzzz_1[j];

                    tey_xzz_zzzz_0[j] = pa_x[j] * tey_zz_zzzz_0[j] - pc_x[j] * tey_zz_zzzz_1[j];

                    tez_xzz_zzzz_0[j] = pa_x[j] * tez_zz_zzzz_0[j] - pc_x[j] * tez_zz_zzzz_1[j];

                    tex_yyy_xxxx_0[j] = pa_y[j] * tex_yy_xxxx_0[j] - pc_y[j] * tex_yy_xxxx_1[j] + fl1_fx * tex_y_xxxx_0[j] - fl1_fx * tex_y_xxxx_1[j];

                    tey_yyy_xxxx_0[j] = pa_y[j] * tey_yy_xxxx_0[j] - pc_y[j] * tey_yy_xxxx_1[j] + fl1_fx * tey_y_xxxx_0[j] - fl1_fx * tey_y_xxxx_1[j] + ta_yy_xxxx_1[j];

                    tez_yyy_xxxx_0[j] = pa_y[j] * tez_yy_xxxx_0[j] - pc_y[j] * tez_yy_xxxx_1[j] + fl1_fx * tez_y_xxxx_0[j] - fl1_fx * tez_y_xxxx_1[j];

                    tex_yyy_xxxy_0[j] = pa_y[j] * tex_yy_xxxy_0[j] - pc_y[j] * tex_yy_xxxy_1[j] + fl1_fx * tex_y_xxxy_0[j] - fl1_fx * tex_y_xxxy_1[j] + 0.5 * fl1_fx * tex_yy_xxx_0[j] - 0.5 * fl1_fx * tex_yy_xxx_1[j];

                    tey_yyy_xxxy_0[j] = pa_y[j] * tey_yy_xxxy_0[j] - pc_y[j] * tey_yy_xxxy_1[j] + fl1_fx * tey_y_xxxy_0[j] - fl1_fx * tey_y_xxxy_1[j] + 0.5 * fl1_fx * tey_yy_xxx_0[j] - 0.5 * fl1_fx * tey_yy_xxx_1[j] + ta_yy_xxxy_1[j];

                    tez_yyy_xxxy_0[j] = pa_y[j] * tez_yy_xxxy_0[j] - pc_y[j] * tez_yy_xxxy_1[j] + fl1_fx * tez_y_xxxy_0[j] - fl1_fx * tez_y_xxxy_1[j] + 0.5 * fl1_fx * tez_yy_xxx_0[j] - 0.5 * fl1_fx * tez_yy_xxx_1[j];

                    tex_yyy_xxxz_0[j] = pa_y[j] * tex_yy_xxxz_0[j] - pc_y[j] * tex_yy_xxxz_1[j] + fl1_fx * tex_y_xxxz_0[j] - fl1_fx * tex_y_xxxz_1[j];

                    tey_yyy_xxxz_0[j] = pa_y[j] * tey_yy_xxxz_0[j] - pc_y[j] * tey_yy_xxxz_1[j] + fl1_fx * tey_y_xxxz_0[j] - fl1_fx * tey_y_xxxz_1[j] + ta_yy_xxxz_1[j];

                    tez_yyy_xxxz_0[j] = pa_y[j] * tez_yy_xxxz_0[j] - pc_y[j] * tez_yy_xxxz_1[j] + fl1_fx * tez_y_xxxz_0[j] - fl1_fx * tez_y_xxxz_1[j];

                    tex_yyy_xxyy_0[j] = pa_y[j] * tex_yy_xxyy_0[j] - pc_y[j] * tex_yy_xxyy_1[j] + fl1_fx * tex_y_xxyy_0[j] - fl1_fx * tex_y_xxyy_1[j] + fl1_fx * tex_yy_xxy_0[j] - fl1_fx * tex_yy_xxy_1[j];

                    tey_yyy_xxyy_0[j] = pa_y[j] * tey_yy_xxyy_0[j] - pc_y[j] * tey_yy_xxyy_1[j] + fl1_fx * tey_y_xxyy_0[j] - fl1_fx * tey_y_xxyy_1[j] + fl1_fx * tey_yy_xxy_0[j] - fl1_fx * tey_yy_xxy_1[j] + ta_yy_xxyy_1[j];

                    tez_yyy_xxyy_0[j] = pa_y[j] * tez_yy_xxyy_0[j] - pc_y[j] * tez_yy_xxyy_1[j] + fl1_fx * tez_y_xxyy_0[j] - fl1_fx * tez_y_xxyy_1[j] + fl1_fx * tez_yy_xxy_0[j] - fl1_fx * tez_yy_xxy_1[j];

                    tex_yyy_xxyz_0[j] = pa_y[j] * tex_yy_xxyz_0[j] - pc_y[j] * tex_yy_xxyz_1[j] + fl1_fx * tex_y_xxyz_0[j] - fl1_fx * tex_y_xxyz_1[j] + 0.5 * fl1_fx * tex_yy_xxz_0[j] - 0.5 * fl1_fx * tex_yy_xxz_1[j];

                    tey_yyy_xxyz_0[j] = pa_y[j] * tey_yy_xxyz_0[j] - pc_y[j] * tey_yy_xxyz_1[j] + fl1_fx * tey_y_xxyz_0[j] - fl1_fx * tey_y_xxyz_1[j] + 0.5 * fl1_fx * tey_yy_xxz_0[j] - 0.5 * fl1_fx * tey_yy_xxz_1[j] + ta_yy_xxyz_1[j];

                    tez_yyy_xxyz_0[j] = pa_y[j] * tez_yy_xxyz_0[j] - pc_y[j] * tez_yy_xxyz_1[j] + fl1_fx * tez_y_xxyz_0[j] - fl1_fx * tez_y_xxyz_1[j] + 0.5 * fl1_fx * tez_yy_xxz_0[j] - 0.5 * fl1_fx * tez_yy_xxz_1[j];

                    tex_yyy_xxzz_0[j] = pa_y[j] * tex_yy_xxzz_0[j] - pc_y[j] * tex_yy_xxzz_1[j] + fl1_fx * tex_y_xxzz_0[j] - fl1_fx * tex_y_xxzz_1[j];

                    tey_yyy_xxzz_0[j] = pa_y[j] * tey_yy_xxzz_0[j] - pc_y[j] * tey_yy_xxzz_1[j] + fl1_fx * tey_y_xxzz_0[j] - fl1_fx * tey_y_xxzz_1[j] + ta_yy_xxzz_1[j];

                    tez_yyy_xxzz_0[j] = pa_y[j] * tez_yy_xxzz_0[j] - pc_y[j] * tez_yy_xxzz_1[j] + fl1_fx * tez_y_xxzz_0[j] - fl1_fx * tez_y_xxzz_1[j];

                    tex_yyy_xyyy_0[j] = pa_y[j] * tex_yy_xyyy_0[j] - pc_y[j] * tex_yy_xyyy_1[j] + fl1_fx * tex_y_xyyy_0[j] - fl1_fx * tex_y_xyyy_1[j] + 1.5 * fl1_fx * tex_yy_xyy_0[j] - 1.5 * fl1_fx * tex_yy_xyy_1[j];

                    tey_yyy_xyyy_0[j] = pa_y[j] * tey_yy_xyyy_0[j] - pc_y[j] * tey_yy_xyyy_1[j] + fl1_fx * tey_y_xyyy_0[j] - fl1_fx * tey_y_xyyy_1[j] + 1.5 * fl1_fx * tey_yy_xyy_0[j] - 1.5 * fl1_fx * tey_yy_xyy_1[j] + ta_yy_xyyy_1[j];

                    tez_yyy_xyyy_0[j] = pa_y[j] * tez_yy_xyyy_0[j] - pc_y[j] * tez_yy_xyyy_1[j] + fl1_fx * tez_y_xyyy_0[j] - fl1_fx * tez_y_xyyy_1[j] + 1.5 * fl1_fx * tez_yy_xyy_0[j] - 1.5 * fl1_fx * tez_yy_xyy_1[j];

                    tex_yyy_xyyz_0[j] = pa_y[j] * tex_yy_xyyz_0[j] - pc_y[j] * tex_yy_xyyz_1[j] + fl1_fx * tex_y_xyyz_0[j] - fl1_fx * tex_y_xyyz_1[j] + fl1_fx * tex_yy_xyz_0[j] - fl1_fx * tex_yy_xyz_1[j];

                    tey_yyy_xyyz_0[j] = pa_y[j] * tey_yy_xyyz_0[j] - pc_y[j] * tey_yy_xyyz_1[j] + fl1_fx * tey_y_xyyz_0[j] - fl1_fx * tey_y_xyyz_1[j] + fl1_fx * tey_yy_xyz_0[j] - fl1_fx * tey_yy_xyz_1[j] + ta_yy_xyyz_1[j];

                    tez_yyy_xyyz_0[j] = pa_y[j] * tez_yy_xyyz_0[j] - pc_y[j] * tez_yy_xyyz_1[j] + fl1_fx * tez_y_xyyz_0[j] - fl1_fx * tez_y_xyyz_1[j] + fl1_fx * tez_yy_xyz_0[j] - fl1_fx * tez_yy_xyz_1[j];

                    tex_yyy_xyzz_0[j] = pa_y[j] * tex_yy_xyzz_0[j] - pc_y[j] * tex_yy_xyzz_1[j] + fl1_fx * tex_y_xyzz_0[j] - fl1_fx * tex_y_xyzz_1[j] + 0.5 * fl1_fx * tex_yy_xzz_0[j] - 0.5 * fl1_fx * tex_yy_xzz_1[j];

                    tey_yyy_xyzz_0[j] = pa_y[j] * tey_yy_xyzz_0[j] - pc_y[j] * tey_yy_xyzz_1[j] + fl1_fx * tey_y_xyzz_0[j] - fl1_fx * tey_y_xyzz_1[j] + 0.5 * fl1_fx * tey_yy_xzz_0[j] - 0.5 * fl1_fx * tey_yy_xzz_1[j] + ta_yy_xyzz_1[j];

                    tez_yyy_xyzz_0[j] = pa_y[j] * tez_yy_xyzz_0[j] - pc_y[j] * tez_yy_xyzz_1[j] + fl1_fx * tez_y_xyzz_0[j] - fl1_fx * tez_y_xyzz_1[j] + 0.5 * fl1_fx * tez_yy_xzz_0[j] - 0.5 * fl1_fx * tez_yy_xzz_1[j];

                    tex_yyy_xzzz_0[j] = pa_y[j] * tex_yy_xzzz_0[j] - pc_y[j] * tex_yy_xzzz_1[j] + fl1_fx * tex_y_xzzz_0[j] - fl1_fx * tex_y_xzzz_1[j];

                    tey_yyy_xzzz_0[j] = pa_y[j] * tey_yy_xzzz_0[j] - pc_y[j] * tey_yy_xzzz_1[j] + fl1_fx * tey_y_xzzz_0[j] - fl1_fx * tey_y_xzzz_1[j] + ta_yy_xzzz_1[j];

                    tez_yyy_xzzz_0[j] = pa_y[j] * tez_yy_xzzz_0[j] - pc_y[j] * tez_yy_xzzz_1[j] + fl1_fx * tez_y_xzzz_0[j] - fl1_fx * tez_y_xzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_300_350(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (300,350)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to tensors product of distances R(PA) = P - A

                auto pa_y = paDistances.data(3 * idx + 1);

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_y = pcDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tex_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 55); 

                auto tey_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 55); 

                auto tez_yy_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 55); 

                auto tex_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 56); 

                auto tey_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 56); 

                auto tez_yy_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 56); 

                auto tex_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 57); 

                auto tey_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 57); 

                auto tez_yy_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 57); 

                auto tex_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 58); 

                auto tey_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 58); 

                auto tez_yy_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 58); 

                auto tex_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 59); 

                auto tey_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 59); 

                auto tez_yy_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 59); 

                auto tex_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 60); 

                auto tey_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 60); 

                auto tez_yz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 60); 

                auto tex_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 61); 

                auto tey_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 61); 

                auto tez_yz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 61); 

                auto tex_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 62); 

                auto tey_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 62); 

                auto tez_yz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 62); 

                auto tex_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 63); 

                auto tey_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 63); 

                auto tez_yz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 63); 

                auto tex_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 64); 

                auto tey_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 64); 

                auto tez_yz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 64); 

                auto tex_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 65); 

                auto tey_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 65); 

                auto tez_yz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 65); 

                auto tex_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 66); 

                auto tey_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 66); 

                auto tez_yz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 66); 

                auto tex_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 67); 

                auto tey_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 67); 

                auto tez_yz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 67); 

                auto tex_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 68); 

                auto tey_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 68); 

                auto tez_yz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 68); 

                auto tex_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 69); 

                auto tey_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 69); 

                auto tez_yz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 69); 

                auto tex_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 70); 

                auto tey_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 70); 

                auto tez_yz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 70); 

                auto tex_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 71); 

                auto tey_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 71); 

                auto tex_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 55); 

                auto tey_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 55); 

                auto tez_yy_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 55); 

                auto tex_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 56); 

                auto tey_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 56); 

                auto tez_yy_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 56); 

                auto tex_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 57); 

                auto tey_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 57); 

                auto tez_yy_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 57); 

                auto tex_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 58); 

                auto tey_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 58); 

                auto tez_yy_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 58); 

                auto tex_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 59); 

                auto tey_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 59); 

                auto tez_yy_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 59); 

                auto tex_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 60); 

                auto tey_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 60); 

                auto tez_yz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 60); 

                auto tex_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 61); 

                auto tey_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 61); 

                auto tez_yz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 61); 

                auto tex_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 62); 

                auto tey_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 62); 

                auto tez_yz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 62); 

                auto tex_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 63); 

                auto tey_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 63); 

                auto tez_yz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 63); 

                auto tex_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 64); 

                auto tey_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 64); 

                auto tez_yz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 64); 

                auto tex_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 65); 

                auto tey_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 65); 

                auto tez_yz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 65); 

                auto tex_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 66); 

                auto tey_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 66); 

                auto tez_yz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 66); 

                auto tex_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 67); 

                auto tey_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 67); 

                auto tez_yz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 67); 

                auto tex_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 68); 

                auto tey_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 68); 

                auto tez_yz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 68); 

                auto tex_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 69); 

                auto tey_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 69); 

                auto tez_yz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 69); 

                auto tex_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 70); 

                auto tey_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 70); 

                auto tez_yz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 70); 

                auto tex_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 71); 

                auto tey_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 71); 

                auto tex_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 25); 

                auto tey_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 25); 

                auto tez_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 25); 

                auto tex_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 26); 

                auto tey_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 26); 

                auto tez_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 26); 

                auto tex_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 27); 

                auto tey_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 27); 

                auto tez_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 27); 

                auto tex_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 28); 

                auto tey_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 28); 

                auto tez_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 28); 

                auto tex_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 29); 

                auto tey_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 29); 

                auto tez_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 29); 

                auto tex_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 30); 

                auto tey_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 30); 

                auto tez_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 30); 

                auto tex_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 31); 

                auto tey_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 31); 

                auto tez_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 31); 

                auto tex_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 32); 

                auto tey_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 32); 

                auto tez_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 32); 

                auto tex_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 33); 

                auto tey_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 33); 

                auto tez_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 33); 

                auto tex_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 34); 

                auto tey_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 34); 

                auto tez_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 34); 

                auto tex_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 35); 

                auto tey_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 35); 

                auto tez_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 35); 

                auto tex_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 36); 

                auto tey_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 36); 

                auto tez_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 36); 

                auto tex_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 37); 

                auto tey_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 37); 

                auto tez_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 37); 

                auto tex_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 38); 

                auto tey_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 38); 

                auto tez_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 38); 

                auto tex_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 39); 

                auto tey_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 39); 

                auto tez_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 39); 

                auto tex_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 40); 

                auto tey_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 40); 

                auto tez_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 40); 

                auto tex_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 41); 

                auto tey_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 41); 

                auto tex_y_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 25); 

                auto tey_y_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 25); 

                auto tez_y_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 25); 

                auto tex_y_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 26); 

                auto tey_y_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 26); 

                auto tez_y_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 26); 

                auto tex_y_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 27); 

                auto tey_y_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 27); 

                auto tez_y_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 27); 

                auto tex_y_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 28); 

                auto tey_y_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 28); 

                auto tez_y_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 28); 

                auto tex_y_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 29); 

                auto tey_y_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 29); 

                auto tez_y_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 29); 

                auto tex_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 30); 

                auto tey_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 30); 

                auto tez_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 30); 

                auto tex_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 31); 

                auto tey_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 31); 

                auto tez_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 31); 

                auto tex_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 32); 

                auto tey_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 32); 

                auto tez_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 32); 

                auto tex_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 33); 

                auto tey_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 33); 

                auto tez_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 33); 

                auto tex_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 34); 

                auto tey_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 34); 

                auto tez_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 34); 

                auto tex_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 35); 

                auto tey_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 35); 

                auto tez_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 35); 

                auto tex_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 36); 

                auto tey_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 36); 

                auto tez_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 36); 

                auto tex_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 37); 

                auto tey_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 37); 

                auto tez_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 37); 

                auto tex_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 38); 

                auto tey_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 38); 

                auto tez_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 38); 

                auto tex_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 39); 

                auto tey_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 39); 

                auto tez_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 39); 

                auto tex_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 40); 

                auto tey_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 40); 

                auto tez_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 40); 

                auto tex_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 41); 

                auto tey_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 41); 

                auto tex_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 36); 

                auto tey_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 36); 

                auto tez_yy_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 36); 

                auto tex_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 37); 

                auto tey_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 37); 

                auto tez_yy_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 37); 

                auto tex_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 38); 

                auto tey_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 38); 

                auto tez_yy_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 38); 

                auto tex_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 39); 

                auto tey_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 39); 

                auto tez_yy_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 39); 

                auto tex_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 40); 

                auto tey_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 40); 

                auto tez_yz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 40); 

                auto tex_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 41); 

                auto tey_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 41); 

                auto tez_yz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 41); 

                auto tex_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 42); 

                auto tey_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 42); 

                auto tez_yz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 42); 

                auto tex_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 43); 

                auto tey_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 43); 

                auto tez_yz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 43); 

                auto tex_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 44); 

                auto tey_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 44); 

                auto tez_yz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 44); 

                auto tex_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 45); 

                auto tey_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 45); 

                auto tez_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 45); 

                auto tex_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 46); 

                auto tey_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 46); 

                auto tez_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 46); 

                auto tex_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 47); 

                auto tey_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 47); 

                auto tex_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 36); 

                auto tey_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 36); 

                auto tez_yy_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 36); 

                auto tex_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 37); 

                auto tey_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 37); 

                auto tez_yy_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 37); 

                auto tex_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 38); 

                auto tey_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 38); 

                auto tez_yy_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 38); 

                auto tex_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 39); 

                auto tey_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 39); 

                auto tez_yy_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 39); 

                auto tex_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 40); 

                auto tey_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 40); 

                auto tez_yz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 40); 

                auto tex_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 41); 

                auto tey_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 41); 

                auto tez_yz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 41); 

                auto tex_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 42); 

                auto tey_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 42); 

                auto tez_yz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 42); 

                auto tex_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 43); 

                auto tey_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 43); 

                auto tez_yz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 43); 

                auto tex_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 44); 

                auto tey_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 44); 

                auto tez_yz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 44); 

                auto tex_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 45); 

                auto tey_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 45); 

                auto tez_yz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 45); 

                auto tex_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 46); 

                auto tey_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 46); 

                auto tez_yz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 46); 

                auto tex_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 47); 

                auto tey_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 47); 

                auto ta_yy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 55); 

                auto ta_yy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 56); 

                auto ta_yy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 57); 

                auto ta_yy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 58); 

                auto ta_yy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 59); 

                auto ta_yz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 60); 

                auto ta_yz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 61); 

                auto ta_yz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 62); 

                auto ta_yz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 63); 

                auto ta_yz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 64); 

                auto ta_yz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 65); 

                auto ta_yz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 66); 

                auto ta_yz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 67); 

                auto ta_yz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 68); 

                auto ta_yz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 69); 

                auto ta_yz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 70); 

                auto ta_yz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 71); 

                // set up pointers to integrals

                auto tex_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 100); 

                auto tey_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 100); 

                auto tez_yyy_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 100); 

                auto tex_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 101); 

                auto tey_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 101); 

                auto tez_yyy_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 101); 

                auto tex_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 102); 

                auto tey_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 102); 

                auto tez_yyy_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 102); 

                auto tex_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 103); 

                auto tey_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 103); 

                auto tez_yyy_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 103); 

                auto tex_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 104); 

                auto tey_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 104); 

                auto tez_yyy_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 104); 

                auto tex_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 105); 

                auto tey_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 105); 

                auto tez_yyz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 105); 

                auto tex_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 106); 

                auto tey_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 106); 

                auto tez_yyz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 106); 

                auto tex_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 107); 

                auto tey_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 107); 

                auto tez_yyz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 107); 

                auto tex_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 108); 

                auto tey_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 108); 

                auto tez_yyz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 108); 

                auto tex_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 109); 

                auto tey_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 109); 

                auto tez_yyz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 109); 

                auto tex_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 110); 

                auto tey_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 110); 

                auto tez_yyz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 110); 

                auto tex_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 111); 

                auto tey_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 111); 

                auto tez_yyz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 111); 

                auto tex_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 112); 

                auto tey_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 112); 

                auto tez_yyz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 112); 

                auto tex_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 113); 

                auto tey_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 113); 

                auto tez_yyz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 113); 

                auto tex_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 114); 

                auto tey_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 114); 

                auto tez_yyz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 114); 

                auto tex_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 115); 

                auto tey_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 115); 

                auto tez_yyz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 115); 

                auto tex_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 116); 

                auto tey_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 116); 

                // Batch of Integrals (300,350)

                #pragma omp simd aligned(fx, pa_y, pc_y, ta_yy_yyyy_1, ta_yy_yyyz_1, ta_yy_yyzz_1, ta_yy_yzzz_1, \
                                         ta_yy_zzzz_1, ta_yz_xxxx_1, ta_yz_xxxy_1, ta_yz_xxxz_1, ta_yz_xxyy_1, ta_yz_xxyz_1, \
                                         ta_yz_xxzz_1, ta_yz_xyyy_1, ta_yz_xyyz_1, ta_yz_xyzz_1, ta_yz_xzzz_1, ta_yz_yyyy_1, \
                                         ta_yz_yyyz_1, tex_y_yyyy_0, tex_y_yyyy_1, tex_y_yyyz_0, tex_y_yyyz_1, tex_y_yyzz_0, \
                                         tex_y_yyzz_1, tex_y_yzzz_0, tex_y_yzzz_1, tex_y_zzzz_0, tex_y_zzzz_1, tex_yy_yyy_0, \
                                         tex_yy_yyy_1, tex_yy_yyyy_0, tex_yy_yyyy_1, tex_yy_yyyz_0, tex_yy_yyyz_1, \
                                         tex_yy_yyz_0, tex_yy_yyz_1, tex_yy_yyzz_0, tex_yy_yyzz_1, tex_yy_yzz_0, \
                                         tex_yy_yzz_1, tex_yy_yzzz_0, tex_yy_yzzz_1, tex_yy_zzz_0, tex_yy_zzz_1, \
                                         tex_yy_zzzz_0, tex_yy_zzzz_1, tex_yyy_yyyy_0, tex_yyy_yyyz_0, tex_yyy_yyzz_0, \
                                         tex_yyy_yzzz_0, tex_yyy_zzzz_0, tex_yyz_xxxx_0, tex_yyz_xxxy_0, tex_yyz_xxxz_0, \
                                         tex_yyz_xxyy_0, tex_yyz_xxyz_0, tex_yyz_xxzz_0, tex_yyz_xyyy_0, tex_yyz_xyyz_0, \
                                         tex_yyz_xyzz_0, tex_yyz_xzzz_0, tex_yyz_yyyy_0, tex_yyz_yyyz_0, tex_yz_xxx_0, \
                                         tex_yz_xxx_1, tex_yz_xxxx_0, tex_yz_xxxx_1, tex_yz_xxxy_0, tex_yz_xxxy_1, \
                                         tex_yz_xxxz_0, tex_yz_xxxz_1, tex_yz_xxy_0, tex_yz_xxy_1, tex_yz_xxyy_0, \
                                         tex_yz_xxyy_1, tex_yz_xxyz_0, tex_yz_xxyz_1, tex_yz_xxz_0, tex_yz_xxz_1, \
                                         tex_yz_xxzz_0, tex_yz_xxzz_1, tex_yz_xyy_0, tex_yz_xyy_1, tex_yz_xyyy_0, \
                                         tex_yz_xyyy_1, tex_yz_xyyz_0, tex_yz_xyyz_1, tex_yz_xyz_0, tex_yz_xyz_1, \
                                         tex_yz_xyzz_0, tex_yz_xyzz_1, tex_yz_xzz_0, tex_yz_xzz_1, tex_yz_xzzz_0, \
                                         tex_yz_xzzz_1, tex_yz_yyy_0, tex_yz_yyy_1, tex_yz_yyyy_0, tex_yz_yyyy_1, \
                                         tex_yz_yyyz_0, tex_yz_yyyz_1, tex_yz_yyz_0, tex_yz_yyz_1, tex_z_xxxx_0, \
                                         tex_z_xxxx_1, tex_z_xxxy_0, tex_z_xxxy_1, tex_z_xxxz_0, tex_z_xxxz_1, tex_z_xxyy_0, \
                                         tex_z_xxyy_1, tex_z_xxyz_0, tex_z_xxyz_1, tex_z_xxzz_0, tex_z_xxzz_1, tex_z_xyyy_0, \
                                         tex_z_xyyy_1, tex_z_xyyz_0, tex_z_xyyz_1, tex_z_xyzz_0, tex_z_xyzz_1, tex_z_xzzz_0, \
                                         tex_z_xzzz_1, tex_z_yyyy_0, tex_z_yyyy_1, tex_z_yyyz_0, tex_z_yyyz_1, tey_y_yyyy_0, \
                                         tey_y_yyyy_1, tey_y_yyyz_0, tey_y_yyyz_1, tey_y_yyzz_0, tey_y_yyzz_1, tey_y_yzzz_0, \
                                         tey_y_yzzz_1, tey_y_zzzz_0, tey_y_zzzz_1, tey_yy_yyy_0, tey_yy_yyy_1, \
                                         tey_yy_yyyy_0, tey_yy_yyyy_1, tey_yy_yyyz_0, tey_yy_yyyz_1, tey_yy_yyz_0, \
                                         tey_yy_yyz_1, tey_yy_yyzz_0, tey_yy_yyzz_1, tey_yy_yzz_0, tey_yy_yzz_1, \
                                         tey_yy_yzzz_0, tey_yy_yzzz_1, tey_yy_zzz_0, tey_yy_zzz_1, tey_yy_zzzz_0, \
                                         tey_yy_zzzz_1, tey_yyy_yyyy_0, tey_yyy_yyyz_0, tey_yyy_yyzz_0, tey_yyy_yzzz_0, \
                                         tey_yyy_zzzz_0, tey_yyz_xxxx_0, tey_yyz_xxxy_0, tey_yyz_xxxz_0, tey_yyz_xxyy_0, \
                                         tey_yyz_xxyz_0, tey_yyz_xxzz_0, tey_yyz_xyyy_0, tey_yyz_xyyz_0, tey_yyz_xyzz_0, \
                                         tey_yyz_xzzz_0, tey_yyz_yyyy_0, tey_yyz_yyyz_0, tey_yz_xxx_0, tey_yz_xxx_1, \
                                         tey_yz_xxxx_0, tey_yz_xxxx_1, tey_yz_xxxy_0, tey_yz_xxxy_1, tey_yz_xxxz_0, \
                                         tey_yz_xxxz_1, tey_yz_xxy_0, tey_yz_xxy_1, tey_yz_xxyy_0, tey_yz_xxyy_1, \
                                         tey_yz_xxyz_0, tey_yz_xxyz_1, tey_yz_xxz_0, tey_yz_xxz_1, tey_yz_xxzz_0, \
                                         tey_yz_xxzz_1, tey_yz_xyy_0, tey_yz_xyy_1, tey_yz_xyyy_0, tey_yz_xyyy_1, \
                                         tey_yz_xyyz_0, tey_yz_xyyz_1, tey_yz_xyz_0, tey_yz_xyz_1, tey_yz_xyzz_0, \
                                         tey_yz_xyzz_1, tey_yz_xzz_0, tey_yz_xzz_1, tey_yz_xzzz_0, tey_yz_xzzz_1, \
                                         tey_yz_yyy_0, tey_yz_yyy_1, tey_yz_yyyy_0, tey_yz_yyyy_1, tey_yz_yyyz_0, \
                                         tey_yz_yyyz_1, tey_yz_yyz_0, tey_yz_yyz_1, tey_z_xxxx_0, tey_z_xxxx_1, tey_z_xxxy_0, \
                                         tey_z_xxxy_1, tey_z_xxxz_0, tey_z_xxxz_1, tey_z_xxyy_0, tey_z_xxyy_1, tey_z_xxyz_0, \
                                         tey_z_xxyz_1, tey_z_xxzz_0, tey_z_xxzz_1, tey_z_xyyy_0, tey_z_xyyy_1, tey_z_xyyz_0, \
                                         tey_z_xyyz_1, tey_z_xyzz_0, tey_z_xyzz_1, tey_z_xzzz_0, tey_z_xzzz_1, tey_z_yyyy_0, \
                                         tey_z_yyyy_1, tey_z_yyyz_0, tey_z_yyyz_1, tez_y_yyyy_0, tez_y_yyyy_1, tez_y_yyyz_0, \
                                         tez_y_yyyz_1, tez_y_yyzz_0, tez_y_yyzz_1, tez_y_yzzz_0, tez_y_yzzz_1, tez_y_zzzz_0, \
                                         tez_y_zzzz_1, tez_yy_yyy_0, tez_yy_yyy_1, tez_yy_yyyy_0, tez_yy_yyyy_1, \
                                         tez_yy_yyyz_0, tez_yy_yyyz_1, tez_yy_yyz_0, tez_yy_yyz_1, tez_yy_yyzz_0, \
                                         tez_yy_yyzz_1, tez_yy_yzz_0, tez_yy_yzz_1, tez_yy_yzzz_0, tez_yy_yzzz_1, \
                                         tez_yy_zzz_0, tez_yy_zzz_1, tez_yy_zzzz_0, tez_yy_zzzz_1, tez_yyy_yyyy_0, \
                                         tez_yyy_yyyz_0, tez_yyy_yyzz_0, tez_yyy_yzzz_0, tez_yyy_zzzz_0, tez_yyz_xxxx_0, \
                                         tez_yyz_xxxy_0, tez_yyz_xxxz_0, tez_yyz_xxyy_0, tez_yyz_xxyz_0, tez_yyz_xxzz_0, \
                                         tez_yyz_xyyy_0, tez_yyz_xyyz_0, tez_yyz_xyzz_0, tez_yyz_xzzz_0, tez_yyz_yyyy_0, \
                                         tez_yz_xxx_0, tez_yz_xxx_1, tez_yz_xxxx_0, tez_yz_xxxx_1, tez_yz_xxxy_0, \
                                         tez_yz_xxxy_1, tez_yz_xxxz_0, tez_yz_xxxz_1, tez_yz_xxy_0, tez_yz_xxy_1, \
                                         tez_yz_xxyy_0, tez_yz_xxyy_1, tez_yz_xxyz_0, tez_yz_xxyz_1, tez_yz_xxz_0, \
                                         tez_yz_xxz_1, tez_yz_xxzz_0, tez_yz_xxzz_1, tez_yz_xyy_0, tez_yz_xyy_1, \
                                         tez_yz_xyyy_0, tez_yz_xyyy_1, tez_yz_xyyz_0, tez_yz_xyyz_1, tez_yz_xyz_0, \
                                         tez_yz_xyz_1, tez_yz_xyzz_0, tez_yz_xyzz_1, tez_yz_xzz_0, tez_yz_xzz_1, \
                                         tez_yz_xzzz_0, tez_yz_xzzz_1, tez_yz_yyy_0, tez_yz_yyy_1, tez_yz_yyyy_0, \
                                         tez_yz_yyyy_1, tez_z_xxxx_0, tez_z_xxxx_1, tez_z_xxxy_0, tez_z_xxxy_1, tez_z_xxxz_0, \
                                         tez_z_xxxz_1, tez_z_xxyy_0, tez_z_xxyy_1, tez_z_xxyz_0, tez_z_xxyz_1, tez_z_xxzz_0, \
                                         tez_z_xxzz_1, tez_z_xyyy_0, tez_z_xyyy_1, tez_z_xyyz_0, tez_z_xyyz_1, tez_z_xyzz_0, \
                                         tez_z_xyzz_1, tez_z_xzzz_0, tez_z_xzzz_1, tez_z_yyyy_0, tez_z_yyyy_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yyy_yyyy_0[j] = pa_y[j] * tex_yy_yyyy_0[j] - pc_y[j] * tex_yy_yyyy_1[j] + fl1_fx * tex_y_yyyy_0[j] - fl1_fx * tex_y_yyyy_1[j] + 2.0 * fl1_fx * tex_yy_yyy_0[j] - 2.0 * fl1_fx * tex_yy_yyy_1[j];

                    tey_yyy_yyyy_0[j] = pa_y[j] * tey_yy_yyyy_0[j] - pc_y[j] * tey_yy_yyyy_1[j] + fl1_fx * tey_y_yyyy_0[j] - fl1_fx * tey_y_yyyy_1[j] + 2.0 * fl1_fx * tey_yy_yyy_0[j] - 2.0 * fl1_fx * tey_yy_yyy_1[j] + ta_yy_yyyy_1[j];

                    tez_yyy_yyyy_0[j] = pa_y[j] * tez_yy_yyyy_0[j] - pc_y[j] * tez_yy_yyyy_1[j] + fl1_fx * tez_y_yyyy_0[j] - fl1_fx * tez_y_yyyy_1[j] + 2.0 * fl1_fx * tez_yy_yyy_0[j] - 2.0 * fl1_fx * tez_yy_yyy_1[j];

                    tex_yyy_yyyz_0[j] = pa_y[j] * tex_yy_yyyz_0[j] - pc_y[j] * tex_yy_yyyz_1[j] + fl1_fx * tex_y_yyyz_0[j] - fl1_fx * tex_y_yyyz_1[j] + 1.5 * fl1_fx * tex_yy_yyz_0[j] - 1.5 * fl1_fx * tex_yy_yyz_1[j];

                    tey_yyy_yyyz_0[j] = pa_y[j] * tey_yy_yyyz_0[j] - pc_y[j] * tey_yy_yyyz_1[j] + fl1_fx * tey_y_yyyz_0[j] - fl1_fx * tey_y_yyyz_1[j] + 1.5 * fl1_fx * tey_yy_yyz_0[j] - 1.5 * fl1_fx * tey_yy_yyz_1[j] + ta_yy_yyyz_1[j];

                    tez_yyy_yyyz_0[j] = pa_y[j] * tez_yy_yyyz_0[j] - pc_y[j] * tez_yy_yyyz_1[j] + fl1_fx * tez_y_yyyz_0[j] - fl1_fx * tez_y_yyyz_1[j] + 1.5 * fl1_fx * tez_yy_yyz_0[j] - 1.5 * fl1_fx * tez_yy_yyz_1[j];

                    tex_yyy_yyzz_0[j] = pa_y[j] * tex_yy_yyzz_0[j] - pc_y[j] * tex_yy_yyzz_1[j] + fl1_fx * tex_y_yyzz_0[j] - fl1_fx * tex_y_yyzz_1[j] + fl1_fx * tex_yy_yzz_0[j] - fl1_fx * tex_yy_yzz_1[j];

                    tey_yyy_yyzz_0[j] = pa_y[j] * tey_yy_yyzz_0[j] - pc_y[j] * tey_yy_yyzz_1[j] + fl1_fx * tey_y_yyzz_0[j] - fl1_fx * tey_y_yyzz_1[j] + fl1_fx * tey_yy_yzz_0[j] - fl1_fx * tey_yy_yzz_1[j] + ta_yy_yyzz_1[j];

                    tez_yyy_yyzz_0[j] = pa_y[j] * tez_yy_yyzz_0[j] - pc_y[j] * tez_yy_yyzz_1[j] + fl1_fx * tez_y_yyzz_0[j] - fl1_fx * tez_y_yyzz_1[j] + fl1_fx * tez_yy_yzz_0[j] - fl1_fx * tez_yy_yzz_1[j];

                    tex_yyy_yzzz_0[j] = pa_y[j] * tex_yy_yzzz_0[j] - pc_y[j] * tex_yy_yzzz_1[j] + fl1_fx * tex_y_yzzz_0[j] - fl1_fx * tex_y_yzzz_1[j] + 0.5 * fl1_fx * tex_yy_zzz_0[j] - 0.5 * fl1_fx * tex_yy_zzz_1[j];

                    tey_yyy_yzzz_0[j] = pa_y[j] * tey_yy_yzzz_0[j] - pc_y[j] * tey_yy_yzzz_1[j] + fl1_fx * tey_y_yzzz_0[j] - fl1_fx * tey_y_yzzz_1[j] + 0.5 * fl1_fx * tey_yy_zzz_0[j] - 0.5 * fl1_fx * tey_yy_zzz_1[j] + ta_yy_yzzz_1[j];

                    tez_yyy_yzzz_0[j] = pa_y[j] * tez_yy_yzzz_0[j] - pc_y[j] * tez_yy_yzzz_1[j] + fl1_fx * tez_y_yzzz_0[j] - fl1_fx * tez_y_yzzz_1[j] + 0.5 * fl1_fx * tez_yy_zzz_0[j] - 0.5 * fl1_fx * tez_yy_zzz_1[j];

                    tex_yyy_zzzz_0[j] = pa_y[j] * tex_yy_zzzz_0[j] - pc_y[j] * tex_yy_zzzz_1[j] + fl1_fx * tex_y_zzzz_0[j] - fl1_fx * tex_y_zzzz_1[j];

                    tey_yyy_zzzz_0[j] = pa_y[j] * tey_yy_zzzz_0[j] - pc_y[j] * tey_yy_zzzz_1[j] + fl1_fx * tey_y_zzzz_0[j] - fl1_fx * tey_y_zzzz_1[j] + ta_yy_zzzz_1[j];

                    tez_yyy_zzzz_0[j] = pa_y[j] * tez_yy_zzzz_0[j] - pc_y[j] * tez_yy_zzzz_1[j] + fl1_fx * tez_y_zzzz_0[j] - fl1_fx * tez_y_zzzz_1[j];

                    tex_yyz_xxxx_0[j] = pa_y[j] * tex_yz_xxxx_0[j] - pc_y[j] * tex_yz_xxxx_1[j] + 0.5 * fl1_fx * tex_z_xxxx_0[j] - 0.5 * fl1_fx * tex_z_xxxx_1[j];

                    tey_yyz_xxxx_0[j] = pa_y[j] * tey_yz_xxxx_0[j] - pc_y[j] * tey_yz_xxxx_1[j] + 0.5 * fl1_fx * tey_z_xxxx_0[j] - 0.5 * fl1_fx * tey_z_xxxx_1[j] + ta_yz_xxxx_1[j];

                    tez_yyz_xxxx_0[j] = pa_y[j] * tez_yz_xxxx_0[j] - pc_y[j] * tez_yz_xxxx_1[j] + 0.5 * fl1_fx * tez_z_xxxx_0[j] - 0.5 * fl1_fx * tez_z_xxxx_1[j];

                    tex_yyz_xxxy_0[j] = pa_y[j] * tex_yz_xxxy_0[j] - pc_y[j] * tex_yz_xxxy_1[j] + 0.5 * fl1_fx * tex_z_xxxy_0[j] - 0.5 * fl1_fx * tex_z_xxxy_1[j] + 0.5 * fl1_fx * tex_yz_xxx_0[j] - 0.5 * fl1_fx * tex_yz_xxx_1[j];

                    tey_yyz_xxxy_0[j] = pa_y[j] * tey_yz_xxxy_0[j] - pc_y[j] * tey_yz_xxxy_1[j] + 0.5 * fl1_fx * tey_z_xxxy_0[j] - 0.5 * fl1_fx * tey_z_xxxy_1[j] + 0.5 * fl1_fx * tey_yz_xxx_0[j] - 0.5 * fl1_fx * tey_yz_xxx_1[j] + ta_yz_xxxy_1[j];

                    tez_yyz_xxxy_0[j] = pa_y[j] * tez_yz_xxxy_0[j] - pc_y[j] * tez_yz_xxxy_1[j] + 0.5 * fl1_fx * tez_z_xxxy_0[j] - 0.5 * fl1_fx * tez_z_xxxy_1[j] + 0.5 * fl1_fx * tez_yz_xxx_0[j] - 0.5 * fl1_fx * tez_yz_xxx_1[j];

                    tex_yyz_xxxz_0[j] = pa_y[j] * tex_yz_xxxz_0[j] - pc_y[j] * tex_yz_xxxz_1[j] + 0.5 * fl1_fx * tex_z_xxxz_0[j] - 0.5 * fl1_fx * tex_z_xxxz_1[j];

                    tey_yyz_xxxz_0[j] = pa_y[j] * tey_yz_xxxz_0[j] - pc_y[j] * tey_yz_xxxz_1[j] + 0.5 * fl1_fx * tey_z_xxxz_0[j] - 0.5 * fl1_fx * tey_z_xxxz_1[j] + ta_yz_xxxz_1[j];

                    tez_yyz_xxxz_0[j] = pa_y[j] * tez_yz_xxxz_0[j] - pc_y[j] * tez_yz_xxxz_1[j] + 0.5 * fl1_fx * tez_z_xxxz_0[j] - 0.5 * fl1_fx * tez_z_xxxz_1[j];

                    tex_yyz_xxyy_0[j] = pa_y[j] * tex_yz_xxyy_0[j] - pc_y[j] * tex_yz_xxyy_1[j] + 0.5 * fl1_fx * tex_z_xxyy_0[j] - 0.5 * fl1_fx * tex_z_xxyy_1[j] + fl1_fx * tex_yz_xxy_0[j] - fl1_fx * tex_yz_xxy_1[j];

                    tey_yyz_xxyy_0[j] = pa_y[j] * tey_yz_xxyy_0[j] - pc_y[j] * tey_yz_xxyy_1[j] + 0.5 * fl1_fx * tey_z_xxyy_0[j] - 0.5 * fl1_fx * tey_z_xxyy_1[j] + fl1_fx * tey_yz_xxy_0[j] - fl1_fx * tey_yz_xxy_1[j] + ta_yz_xxyy_1[j];

                    tez_yyz_xxyy_0[j] = pa_y[j] * tez_yz_xxyy_0[j] - pc_y[j] * tez_yz_xxyy_1[j] + 0.5 * fl1_fx * tez_z_xxyy_0[j] - 0.5 * fl1_fx * tez_z_xxyy_1[j] + fl1_fx * tez_yz_xxy_0[j] - fl1_fx * tez_yz_xxy_1[j];

                    tex_yyz_xxyz_0[j] = pa_y[j] * tex_yz_xxyz_0[j] - pc_y[j] * tex_yz_xxyz_1[j] + 0.5 * fl1_fx * tex_z_xxyz_0[j] - 0.5 * fl1_fx * tex_z_xxyz_1[j] + 0.5 * fl1_fx * tex_yz_xxz_0[j] - 0.5 * fl1_fx * tex_yz_xxz_1[j];

                    tey_yyz_xxyz_0[j] = pa_y[j] * tey_yz_xxyz_0[j] - pc_y[j] * tey_yz_xxyz_1[j] + 0.5 * fl1_fx * tey_z_xxyz_0[j] - 0.5 * fl1_fx * tey_z_xxyz_1[j] + 0.5 * fl1_fx * tey_yz_xxz_0[j] - 0.5 * fl1_fx * tey_yz_xxz_1[j] + ta_yz_xxyz_1[j];

                    tez_yyz_xxyz_0[j] = pa_y[j] * tez_yz_xxyz_0[j] - pc_y[j] * tez_yz_xxyz_1[j] + 0.5 * fl1_fx * tez_z_xxyz_0[j] - 0.5 * fl1_fx * tez_z_xxyz_1[j] + 0.5 * fl1_fx * tez_yz_xxz_0[j] - 0.5 * fl1_fx * tez_yz_xxz_1[j];

                    tex_yyz_xxzz_0[j] = pa_y[j] * tex_yz_xxzz_0[j] - pc_y[j] * tex_yz_xxzz_1[j] + 0.5 * fl1_fx * tex_z_xxzz_0[j] - 0.5 * fl1_fx * tex_z_xxzz_1[j];

                    tey_yyz_xxzz_0[j] = pa_y[j] * tey_yz_xxzz_0[j] - pc_y[j] * tey_yz_xxzz_1[j] + 0.5 * fl1_fx * tey_z_xxzz_0[j] - 0.5 * fl1_fx * tey_z_xxzz_1[j] + ta_yz_xxzz_1[j];

                    tez_yyz_xxzz_0[j] = pa_y[j] * tez_yz_xxzz_0[j] - pc_y[j] * tez_yz_xxzz_1[j] + 0.5 * fl1_fx * tez_z_xxzz_0[j] - 0.5 * fl1_fx * tez_z_xxzz_1[j];

                    tex_yyz_xyyy_0[j] = pa_y[j] * tex_yz_xyyy_0[j] - pc_y[j] * tex_yz_xyyy_1[j] + 0.5 * fl1_fx * tex_z_xyyy_0[j] - 0.5 * fl1_fx * tex_z_xyyy_1[j] + 1.5 * fl1_fx * tex_yz_xyy_0[j] - 1.5 * fl1_fx * tex_yz_xyy_1[j];

                    tey_yyz_xyyy_0[j] = pa_y[j] * tey_yz_xyyy_0[j] - pc_y[j] * tey_yz_xyyy_1[j] + 0.5 * fl1_fx * tey_z_xyyy_0[j] - 0.5 * fl1_fx * tey_z_xyyy_1[j] + 1.5 * fl1_fx * tey_yz_xyy_0[j] - 1.5 * fl1_fx * tey_yz_xyy_1[j] + ta_yz_xyyy_1[j];

                    tez_yyz_xyyy_0[j] = pa_y[j] * tez_yz_xyyy_0[j] - pc_y[j] * tez_yz_xyyy_1[j] + 0.5 * fl1_fx * tez_z_xyyy_0[j] - 0.5 * fl1_fx * tez_z_xyyy_1[j] + 1.5 * fl1_fx * tez_yz_xyy_0[j] - 1.5 * fl1_fx * tez_yz_xyy_1[j];

                    tex_yyz_xyyz_0[j] = pa_y[j] * tex_yz_xyyz_0[j] - pc_y[j] * tex_yz_xyyz_1[j] + 0.5 * fl1_fx * tex_z_xyyz_0[j] - 0.5 * fl1_fx * tex_z_xyyz_1[j] + fl1_fx * tex_yz_xyz_0[j] - fl1_fx * tex_yz_xyz_1[j];

                    tey_yyz_xyyz_0[j] = pa_y[j] * tey_yz_xyyz_0[j] - pc_y[j] * tey_yz_xyyz_1[j] + 0.5 * fl1_fx * tey_z_xyyz_0[j] - 0.5 * fl1_fx * tey_z_xyyz_1[j] + fl1_fx * tey_yz_xyz_0[j] - fl1_fx * tey_yz_xyz_1[j] + ta_yz_xyyz_1[j];

                    tez_yyz_xyyz_0[j] = pa_y[j] * tez_yz_xyyz_0[j] - pc_y[j] * tez_yz_xyyz_1[j] + 0.5 * fl1_fx * tez_z_xyyz_0[j] - 0.5 * fl1_fx * tez_z_xyyz_1[j] + fl1_fx * tez_yz_xyz_0[j] - fl1_fx * tez_yz_xyz_1[j];

                    tex_yyz_xyzz_0[j] = pa_y[j] * tex_yz_xyzz_0[j] - pc_y[j] * tex_yz_xyzz_1[j] + 0.5 * fl1_fx * tex_z_xyzz_0[j] - 0.5 * fl1_fx * tex_z_xyzz_1[j] + 0.5 * fl1_fx * tex_yz_xzz_0[j] - 0.5 * fl1_fx * tex_yz_xzz_1[j];

                    tey_yyz_xyzz_0[j] = pa_y[j] * tey_yz_xyzz_0[j] - pc_y[j] * tey_yz_xyzz_1[j] + 0.5 * fl1_fx * tey_z_xyzz_0[j] - 0.5 * fl1_fx * tey_z_xyzz_1[j] + 0.5 * fl1_fx * tey_yz_xzz_0[j] - 0.5 * fl1_fx * tey_yz_xzz_1[j] + ta_yz_xyzz_1[j];

                    tez_yyz_xyzz_0[j] = pa_y[j] * tez_yz_xyzz_0[j] - pc_y[j] * tez_yz_xyzz_1[j] + 0.5 * fl1_fx * tez_z_xyzz_0[j] - 0.5 * fl1_fx * tez_z_xyzz_1[j] + 0.5 * fl1_fx * tez_yz_xzz_0[j] - 0.5 * fl1_fx * tez_yz_xzz_1[j];

                    tex_yyz_xzzz_0[j] = pa_y[j] * tex_yz_xzzz_0[j] - pc_y[j] * tex_yz_xzzz_1[j] + 0.5 * fl1_fx * tex_z_xzzz_0[j] - 0.5 * fl1_fx * tex_z_xzzz_1[j];

                    tey_yyz_xzzz_0[j] = pa_y[j] * tey_yz_xzzz_0[j] - pc_y[j] * tey_yz_xzzz_1[j] + 0.5 * fl1_fx * tey_z_xzzz_0[j] - 0.5 * fl1_fx * tey_z_xzzz_1[j] + ta_yz_xzzz_1[j];

                    tez_yyz_xzzz_0[j] = pa_y[j] * tez_yz_xzzz_0[j] - pc_y[j] * tez_yz_xzzz_1[j] + 0.5 * fl1_fx * tez_z_xzzz_0[j] - 0.5 * fl1_fx * tez_z_xzzz_1[j];

                    tex_yyz_yyyy_0[j] = pa_y[j] * tex_yz_yyyy_0[j] - pc_y[j] * tex_yz_yyyy_1[j] + 0.5 * fl1_fx * tex_z_yyyy_0[j] - 0.5 * fl1_fx * tex_z_yyyy_1[j] + 2.0 * fl1_fx * tex_yz_yyy_0[j] - 2.0 * fl1_fx * tex_yz_yyy_1[j];

                    tey_yyz_yyyy_0[j] = pa_y[j] * tey_yz_yyyy_0[j] - pc_y[j] * tey_yz_yyyy_1[j] + 0.5 * fl1_fx * tey_z_yyyy_0[j] - 0.5 * fl1_fx * tey_z_yyyy_1[j] + 2.0 * fl1_fx * tey_yz_yyy_0[j] - 2.0 * fl1_fx * tey_yz_yyy_1[j] + ta_yz_yyyy_1[j];

                    tez_yyz_yyyy_0[j] = pa_y[j] * tez_yz_yyyy_0[j] - pc_y[j] * tez_yz_yyyy_1[j] + 0.5 * fl1_fx * tez_z_yyyy_0[j] - 0.5 * fl1_fx * tez_z_yyyy_1[j] + 2.0 * fl1_fx * tez_yz_yyy_0[j] - 2.0 * fl1_fx * tez_yz_yyy_1[j];

                    tex_yyz_yyyz_0[j] = pa_y[j] * tex_yz_yyyz_0[j] - pc_y[j] * tex_yz_yyyz_1[j] + 0.5 * fl1_fx * tex_z_yyyz_0[j] - 0.5 * fl1_fx * tex_z_yyyz_1[j] + 1.5 * fl1_fx * tex_yz_yyz_0[j] - 1.5 * fl1_fx * tex_yz_yyz_1[j];

                    tey_yyz_yyyz_0[j] = pa_y[j] * tey_yz_yyyz_0[j] - pc_y[j] * tey_yz_yyyz_1[j] + 0.5 * fl1_fx * tey_z_yyyz_0[j] - 0.5 * fl1_fx * tey_z_yyyz_1[j] + 1.5 * fl1_fx * tey_yz_yyz_0[j] - 1.5 * fl1_fx * tey_yz_yyz_1[j] + ta_yz_yyyz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_350_400(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (350,400)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to tensors product of distances R(PA) = P - A

                auto pa_y = paDistances.data(3 * idx + 1);

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_y = pcDistances.data(3 * idx + 1);

                // set up pointers to auxilary integrals

                auto tez_yz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 71); 

                auto tex_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 72); 

                auto tey_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 72); 

                auto tez_yz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 72); 

                auto tex_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 73); 

                auto tey_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 73); 

                auto tez_yz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 73); 

                auto tex_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 74); 

                auto tey_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 74); 

                auto tez_yz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 74); 

                auto tex_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 75); 

                auto tey_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 75); 

                auto tez_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 75); 

                auto tex_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 76); 

                auto tey_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 76); 

                auto tez_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 76); 

                auto tex_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 77); 

                auto tey_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 77); 

                auto tez_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 77); 

                auto tex_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 78); 

                auto tey_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 78); 

                auto tez_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 78); 

                auto tex_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 79); 

                auto tey_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 79); 

                auto tez_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 79); 

                auto tex_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 80); 

                auto tey_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 80); 

                auto tez_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 80); 

                auto tex_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 81); 

                auto tey_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 81); 

                auto tez_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 81); 

                auto tex_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 82); 

                auto tey_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 82); 

                auto tez_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 82); 

                auto tex_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 83); 

                auto tey_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 83); 

                auto tez_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 83); 

                auto tex_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 84); 

                auto tey_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 84); 

                auto tez_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 84); 

                auto tex_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 85); 

                auto tey_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 85); 

                auto tez_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 85); 

                auto tex_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 86); 

                auto tey_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 86); 

                auto tez_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 86); 

                auto tex_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 87); 

                auto tey_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 87); 

                auto tez_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 87); 

                auto tex_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 88); 

                auto tez_yz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 71); 

                auto tex_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 72); 

                auto tey_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 72); 

                auto tez_yz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 72); 

                auto tex_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 73); 

                auto tey_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 73); 

                auto tez_yz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 73); 

                auto tex_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 74); 

                auto tey_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 74); 

                auto tez_yz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 74); 

                auto tex_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 75); 

                auto tey_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 75); 

                auto tez_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 75); 

                auto tex_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 76); 

                auto tey_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 76); 

                auto tez_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 76); 

                auto tex_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 77); 

                auto tey_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 77); 

                auto tez_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 77); 

                auto tex_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 78); 

                auto tey_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 78); 

                auto tez_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 78); 

                auto tex_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 79); 

                auto tey_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 79); 

                auto tez_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 79); 

                auto tex_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 80); 

                auto tey_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 80); 

                auto tez_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 80); 

                auto tex_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 81); 

                auto tey_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 81); 

                auto tez_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 81); 

                auto tex_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 82); 

                auto tey_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 82); 

                auto tez_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 82); 

                auto tex_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 83); 

                auto tey_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 83); 

                auto tez_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 83); 

                auto tex_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 84); 

                auto tey_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 84); 

                auto tez_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 84); 

                auto tex_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 85); 

                auto tey_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 85); 

                auto tez_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 85); 

                auto tex_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 86); 

                auto tey_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 86); 

                auto tez_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 86); 

                auto tex_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 87); 

                auto tey_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 87); 

                auto tez_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 87); 

                auto tex_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 88); 

                auto tez_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 41); 

                auto tex_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 42); 

                auto tey_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 42); 

                auto tez_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 42); 

                auto tex_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 43); 

                auto tey_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 43); 

                auto tez_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 43); 

                auto tex_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 44); 

                auto tey_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 44); 

                auto tez_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 44); 

                auto tez_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 41); 

                auto tex_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 42); 

                auto tey_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 42); 

                auto tez_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 42); 

                auto tex_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 43); 

                auto tey_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 43); 

                auto tez_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 43); 

                auto tex_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 44); 

                auto tey_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 44); 

                auto tez_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 44); 

                auto tez_yz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 47); 

                auto tex_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 48); 

                auto tey_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 48); 

                auto tez_yz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 48); 

                auto tex_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 49); 

                auto tey_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 49); 

                auto tez_yz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 49); 

                auto tex_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 50); 

                auto tey_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 50); 

                auto tez_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 50); 

                auto tex_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 51); 

                auto tey_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 51); 

                auto tez_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 51); 

                auto tex_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 52); 

                auto tey_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 52); 

                auto tez_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 52); 

                auto tex_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 53); 

                auto tey_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 53); 

                auto tez_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 53); 

                auto tex_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 54); 

                auto tey_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 54); 

                auto tez_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 54); 

                auto tex_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 55); 

                auto tey_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 55); 

                auto tez_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 55); 

                auto tex_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 56); 

                auto tey_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 56); 

                auto tez_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 56); 

                auto tex_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 57); 

                auto tey_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 57); 

                auto tez_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 57); 

                auto tex_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 58); 

                auto tey_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 58); 

                auto tez_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 58); 

                auto tex_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 59); 

                auto tez_yz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 47); 

                auto tex_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 48); 

                auto tey_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 48); 

                auto tez_yz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 48); 

                auto tex_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 49); 

                auto tey_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 49); 

                auto tez_yz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 49); 

                auto tex_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 50); 

                auto tey_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 50); 

                auto tez_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 50); 

                auto tex_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 51); 

                auto tey_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 51); 

                auto tez_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 51); 

                auto tex_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 52); 

                auto tey_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 52); 

                auto tez_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 52); 

                auto tex_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 53); 

                auto tey_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 53); 

                auto tez_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 53); 

                auto tex_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 54); 

                auto tey_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 54); 

                auto tez_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 54); 

                auto tex_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 55); 

                auto tey_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 55); 

                auto tez_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 55); 

                auto tex_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 56); 

                auto tey_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 56); 

                auto tez_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 56); 

                auto tex_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 57); 

                auto tey_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 57); 

                auto tez_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 57); 

                auto tex_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 58); 

                auto tey_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 58); 

                auto tez_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 58); 

                auto tex_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 59); 

                auto ta_yz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 72); 

                auto ta_yz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 73); 

                auto ta_yz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 74); 

                auto ta_zz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 75); 

                auto ta_zz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 76); 

                auto ta_zz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 77); 

                auto ta_zz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 78); 

                auto ta_zz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 79); 

                auto ta_zz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 80); 

                auto ta_zz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 81); 

                auto ta_zz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 82); 

                auto ta_zz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 83); 

                auto ta_zz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 84); 

                auto ta_zz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 85); 

                auto ta_zz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 86); 

                auto ta_zz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 87); 

                // set up pointers to integrals

                auto tez_yyz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 116); 

                auto tex_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 117); 

                auto tey_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 117); 

                auto tez_yyz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 117); 

                auto tex_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 118); 

                auto tey_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 118); 

                auto tez_yyz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 118); 

                auto tex_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 119); 

                auto tey_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 119); 

                auto tez_yyz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 119); 

                auto tex_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 120); 

                auto tey_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 120); 

                auto tez_yzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 120); 

                auto tex_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 121); 

                auto tey_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 121); 

                auto tez_yzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 121); 

                auto tex_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 122); 

                auto tey_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 122); 

                auto tez_yzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 122); 

                auto tex_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 123); 

                auto tey_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 123); 

                auto tez_yzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 123); 

                auto tex_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 124); 

                auto tey_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 124); 

                auto tez_yzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 124); 

                auto tex_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 125); 

                auto tey_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 125); 

                auto tez_yzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 125); 

                auto tex_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 126); 

                auto tey_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 126); 

                auto tez_yzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 126); 

                auto tex_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 127); 

                auto tey_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 127); 

                auto tez_yzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 127); 

                auto tex_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 128); 

                auto tey_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 128); 

                auto tez_yzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 128); 

                auto tex_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 129); 

                auto tey_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 129); 

                auto tez_yzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 129); 

                auto tex_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 130); 

                auto tey_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 130); 

                auto tez_yzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 130); 

                auto tex_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 131); 

                auto tey_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 131); 

                auto tez_yzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 131); 

                auto tex_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 132); 

                auto tey_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 132); 

                auto tez_yzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 132); 

                auto tex_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 133); 

                // Batch of Integrals (350,400)

                #pragma omp simd aligned(fx, pa_y, pc_y, ta_yz_yyzz_1, ta_yz_yzzz_1, ta_yz_zzzz_1, ta_zz_xxxx_1, \
                                         ta_zz_xxxy_1, ta_zz_xxxz_1, ta_zz_xxyy_1, ta_zz_xxyz_1, ta_zz_xxzz_1, ta_zz_xyyy_1, \
                                         ta_zz_xyyz_1, ta_zz_xyzz_1, ta_zz_xzzz_1, ta_zz_yyyy_1, ta_zz_yyyz_1, ta_zz_yyzz_1, \
                                         tex_yyz_yyzz_0, tex_yyz_yzzz_0, tex_yyz_zzzz_0, tex_yz_yyzz_0, tex_yz_yyzz_1, \
                                         tex_yz_yzz_0, tex_yz_yzz_1, tex_yz_yzzz_0, tex_yz_yzzz_1, tex_yz_zzz_0, \
                                         tex_yz_zzz_1, tex_yz_zzzz_0, tex_yz_zzzz_1, tex_yzz_xxxx_0, tex_yzz_xxxy_0, \
                                         tex_yzz_xxxz_0, tex_yzz_xxyy_0, tex_yzz_xxyz_0, tex_yzz_xxzz_0, tex_yzz_xyyy_0, \
                                         tex_yzz_xyyz_0, tex_yzz_xyzz_0, tex_yzz_xzzz_0, tex_yzz_yyyy_0, tex_yzz_yyyz_0, \
                                         tex_yzz_yyzz_0, tex_yzz_yzzz_0, tex_z_yyzz_0, tex_z_yyzz_1, tex_z_yzzz_0, \
                                         tex_z_yzzz_1, tex_z_zzzz_0, tex_z_zzzz_1, tex_zz_xxx_0, tex_zz_xxx_1, \
                                         tex_zz_xxxx_0, tex_zz_xxxx_1, tex_zz_xxxy_0, tex_zz_xxxy_1, tex_zz_xxxz_0, \
                                         tex_zz_xxxz_1, tex_zz_xxy_0, tex_zz_xxy_1, tex_zz_xxyy_0, tex_zz_xxyy_1, \
                                         tex_zz_xxyz_0, tex_zz_xxyz_1, tex_zz_xxz_0, tex_zz_xxz_1, tex_zz_xxzz_0, \
                                         tex_zz_xxzz_1, tex_zz_xyy_0, tex_zz_xyy_1, tex_zz_xyyy_0, tex_zz_xyyy_1, \
                                         tex_zz_xyyz_0, tex_zz_xyyz_1, tex_zz_xyz_0, tex_zz_xyz_1, tex_zz_xyzz_0, \
                                         tex_zz_xyzz_1, tex_zz_xzz_0, tex_zz_xzz_1, tex_zz_xzzz_0, tex_zz_xzzz_1, \
                                         tex_zz_yyy_0, tex_zz_yyy_1, tex_zz_yyyy_0, tex_zz_yyyy_1, tex_zz_yyyz_0, \
                                         tex_zz_yyyz_1, tex_zz_yyz_0, tex_zz_yyz_1, tex_zz_yyzz_0, tex_zz_yyzz_1, \
                                         tex_zz_yzz_0, tex_zz_yzz_1, tex_zz_yzzz_0, tex_zz_yzzz_1, tex_zz_zzz_0, \
                                         tex_zz_zzz_1, tey_yyz_yyzz_0, tey_yyz_yzzz_0, tey_yyz_zzzz_0, tey_yz_yyzz_0, \
                                         tey_yz_yyzz_1, tey_yz_yzz_0, tey_yz_yzz_1, tey_yz_yzzz_0, tey_yz_yzzz_1, \
                                         tey_yz_zzz_0, tey_yz_zzz_1, tey_yz_zzzz_0, tey_yz_zzzz_1, tey_yzz_xxxx_0, \
                                         tey_yzz_xxxy_0, tey_yzz_xxxz_0, tey_yzz_xxyy_0, tey_yzz_xxyz_0, tey_yzz_xxzz_0, \
                                         tey_yzz_xyyy_0, tey_yzz_xyyz_0, tey_yzz_xyzz_0, tey_yzz_xzzz_0, tey_yzz_yyyy_0, \
                                         tey_yzz_yyyz_0, tey_yzz_yyzz_0, tey_z_yyzz_0, tey_z_yyzz_1, tey_z_yzzz_0, \
                                         tey_z_yzzz_1, tey_z_zzzz_0, tey_z_zzzz_1, tey_zz_xxx_0, tey_zz_xxx_1, \
                                         tey_zz_xxxx_0, tey_zz_xxxx_1, tey_zz_xxxy_0, tey_zz_xxxy_1, tey_zz_xxxz_0, \
                                         tey_zz_xxxz_1, tey_zz_xxy_0, tey_zz_xxy_1, tey_zz_xxyy_0, tey_zz_xxyy_1, \
                                         tey_zz_xxyz_0, tey_zz_xxyz_1, tey_zz_xxz_0, tey_zz_xxz_1, tey_zz_xxzz_0, \
                                         tey_zz_xxzz_1, tey_zz_xyy_0, tey_zz_xyy_1, tey_zz_xyyy_0, tey_zz_xyyy_1, \
                                         tey_zz_xyyz_0, tey_zz_xyyz_1, tey_zz_xyz_0, tey_zz_xyz_1, tey_zz_xyzz_0, \
                                         tey_zz_xyzz_1, tey_zz_xzz_0, tey_zz_xzz_1, tey_zz_xzzz_0, tey_zz_xzzz_1, \
                                         tey_zz_yyy_0, tey_zz_yyy_1, tey_zz_yyyy_0, tey_zz_yyyy_1, tey_zz_yyyz_0, \
                                         tey_zz_yyyz_1, tey_zz_yyz_0, tey_zz_yyz_1, tey_zz_yyzz_0, tey_zz_yyzz_1, \
                                         tey_zz_yzz_0, tey_zz_yzz_1, tez_yyz_yyyz_0, tez_yyz_yyzz_0, tez_yyz_yzzz_0, \
                                         tez_yyz_zzzz_0, tez_yz_yyyz_0, tez_yz_yyyz_1, tez_yz_yyz_0, tez_yz_yyz_1, \
                                         tez_yz_yyzz_0, tez_yz_yyzz_1, tez_yz_yzz_0, tez_yz_yzz_1, tez_yz_yzzz_0, \
                                         tez_yz_yzzz_1, tez_yz_zzz_0, tez_yz_zzz_1, tez_yz_zzzz_0, tez_yz_zzzz_1, \
                                         tez_yzz_xxxx_0, tez_yzz_xxxy_0, tez_yzz_xxxz_0, tez_yzz_xxyy_0, tez_yzz_xxyz_0, \
                                         tez_yzz_xxzz_0, tez_yzz_xyyy_0, tez_yzz_xyyz_0, tez_yzz_xyzz_0, tez_yzz_xzzz_0, \
                                         tez_yzz_yyyy_0, tez_yzz_yyyz_0, tez_yzz_yyzz_0, tez_z_yyyz_0, tez_z_yyyz_1, \
                                         tez_z_yyzz_0, tez_z_yyzz_1, tez_z_yzzz_0, tez_z_yzzz_1, tez_z_zzzz_0, tez_z_zzzz_1, \
                                         tez_zz_xxx_0, tez_zz_xxx_1, tez_zz_xxxx_0, tez_zz_xxxx_1, tez_zz_xxxy_0, \
                                         tez_zz_xxxy_1, tez_zz_xxxz_0, tez_zz_xxxz_1, tez_zz_xxy_0, tez_zz_xxy_1, \
                                         tez_zz_xxyy_0, tez_zz_xxyy_1, tez_zz_xxyz_0, tez_zz_xxyz_1, tez_zz_xxz_0, \
                                         tez_zz_xxz_1, tez_zz_xxzz_0, tez_zz_xxzz_1, tez_zz_xyy_0, tez_zz_xyy_1, \
                                         tez_zz_xyyy_0, tez_zz_xyyy_1, tez_zz_xyyz_0, tez_zz_xyyz_1, tez_zz_xyz_0, \
                                         tez_zz_xyz_1, tez_zz_xyzz_0, tez_zz_xyzz_1, tez_zz_xzz_0, tez_zz_xzz_1, \
                                         tez_zz_xzzz_0, tez_zz_xzzz_1, tez_zz_yyy_0, tez_zz_yyy_1, tez_zz_yyyy_0, \
                                         tez_zz_yyyy_1, tez_zz_yyyz_0, tez_zz_yyyz_1, tez_zz_yyz_0, tez_zz_yyz_1, \
                                         tez_zz_yyzz_0, tez_zz_yyzz_1, tez_zz_yzz_0, tez_zz_yzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tez_yyz_yyyz_0[j] = pa_y[j] * tez_yz_yyyz_0[j] - pc_y[j] * tez_yz_yyyz_1[j] + 0.5 * fl1_fx * tez_z_yyyz_0[j] - 0.5 * fl1_fx * tez_z_yyyz_1[j] + 1.5 * fl1_fx * tez_yz_yyz_0[j] - 1.5 * fl1_fx * tez_yz_yyz_1[j];

                    tex_yyz_yyzz_0[j] = pa_y[j] * tex_yz_yyzz_0[j] - pc_y[j] * tex_yz_yyzz_1[j] + 0.5 * fl1_fx * tex_z_yyzz_0[j] - 0.5 * fl1_fx * tex_z_yyzz_1[j] + fl1_fx * tex_yz_yzz_0[j] - fl1_fx * tex_yz_yzz_1[j];

                    tey_yyz_yyzz_0[j] = pa_y[j] * tey_yz_yyzz_0[j] - pc_y[j] * tey_yz_yyzz_1[j] + 0.5 * fl1_fx * tey_z_yyzz_0[j] - 0.5 * fl1_fx * tey_z_yyzz_1[j] + fl1_fx * tey_yz_yzz_0[j] - fl1_fx * tey_yz_yzz_1[j] + ta_yz_yyzz_1[j];

                    tez_yyz_yyzz_0[j] = pa_y[j] * tez_yz_yyzz_0[j] - pc_y[j] * tez_yz_yyzz_1[j] + 0.5 * fl1_fx * tez_z_yyzz_0[j] - 0.5 * fl1_fx * tez_z_yyzz_1[j] + fl1_fx * tez_yz_yzz_0[j] - fl1_fx * tez_yz_yzz_1[j];

                    tex_yyz_yzzz_0[j] = pa_y[j] * tex_yz_yzzz_0[j] - pc_y[j] * tex_yz_yzzz_1[j] + 0.5 * fl1_fx * tex_z_yzzz_0[j] - 0.5 * fl1_fx * tex_z_yzzz_1[j] + 0.5 * fl1_fx * tex_yz_zzz_0[j] - 0.5 * fl1_fx * tex_yz_zzz_1[j];

                    tey_yyz_yzzz_0[j] = pa_y[j] * tey_yz_yzzz_0[j] - pc_y[j] * tey_yz_yzzz_1[j] + 0.5 * fl1_fx * tey_z_yzzz_0[j] - 0.5 * fl1_fx * tey_z_yzzz_1[j] + 0.5 * fl1_fx * tey_yz_zzz_0[j] - 0.5 * fl1_fx * tey_yz_zzz_1[j] + ta_yz_yzzz_1[j];

                    tez_yyz_yzzz_0[j] = pa_y[j] * tez_yz_yzzz_0[j] - pc_y[j] * tez_yz_yzzz_1[j] + 0.5 * fl1_fx * tez_z_yzzz_0[j] - 0.5 * fl1_fx * tez_z_yzzz_1[j] + 0.5 * fl1_fx * tez_yz_zzz_0[j] - 0.5 * fl1_fx * tez_yz_zzz_1[j];

                    tex_yyz_zzzz_0[j] = pa_y[j] * tex_yz_zzzz_0[j] - pc_y[j] * tex_yz_zzzz_1[j] + 0.5 * fl1_fx * tex_z_zzzz_0[j] - 0.5 * fl1_fx * tex_z_zzzz_1[j];

                    tey_yyz_zzzz_0[j] = pa_y[j] * tey_yz_zzzz_0[j] - pc_y[j] * tey_yz_zzzz_1[j] + 0.5 * fl1_fx * tey_z_zzzz_0[j] - 0.5 * fl1_fx * tey_z_zzzz_1[j] + ta_yz_zzzz_1[j];

                    tez_yyz_zzzz_0[j] = pa_y[j] * tez_yz_zzzz_0[j] - pc_y[j] * tez_yz_zzzz_1[j] + 0.5 * fl1_fx * tez_z_zzzz_0[j] - 0.5 * fl1_fx * tez_z_zzzz_1[j];

                    tex_yzz_xxxx_0[j] = pa_y[j] * tex_zz_xxxx_0[j] - pc_y[j] * tex_zz_xxxx_1[j];

                    tey_yzz_xxxx_0[j] = pa_y[j] * tey_zz_xxxx_0[j] - pc_y[j] * tey_zz_xxxx_1[j] + ta_zz_xxxx_1[j];

                    tez_yzz_xxxx_0[j] = pa_y[j] * tez_zz_xxxx_0[j] - pc_y[j] * tez_zz_xxxx_1[j];

                    tex_yzz_xxxy_0[j] = pa_y[j] * tex_zz_xxxy_0[j] - pc_y[j] * tex_zz_xxxy_1[j] + 0.5 * fl1_fx * tex_zz_xxx_0[j] - 0.5 * fl1_fx * tex_zz_xxx_1[j];

                    tey_yzz_xxxy_0[j] = pa_y[j] * tey_zz_xxxy_0[j] - pc_y[j] * tey_zz_xxxy_1[j] + 0.5 * fl1_fx * tey_zz_xxx_0[j] - 0.5 * fl1_fx * tey_zz_xxx_1[j] + ta_zz_xxxy_1[j];

                    tez_yzz_xxxy_0[j] = pa_y[j] * tez_zz_xxxy_0[j] - pc_y[j] * tez_zz_xxxy_1[j] + 0.5 * fl1_fx * tez_zz_xxx_0[j] - 0.5 * fl1_fx * tez_zz_xxx_1[j];

                    tex_yzz_xxxz_0[j] = pa_y[j] * tex_zz_xxxz_0[j] - pc_y[j] * tex_zz_xxxz_1[j];

                    tey_yzz_xxxz_0[j] = pa_y[j] * tey_zz_xxxz_0[j] - pc_y[j] * tey_zz_xxxz_1[j] + ta_zz_xxxz_1[j];

                    tez_yzz_xxxz_0[j] = pa_y[j] * tez_zz_xxxz_0[j] - pc_y[j] * tez_zz_xxxz_1[j];

                    tex_yzz_xxyy_0[j] = pa_y[j] * tex_zz_xxyy_0[j] - pc_y[j] * tex_zz_xxyy_1[j] + fl1_fx * tex_zz_xxy_0[j] - fl1_fx * tex_zz_xxy_1[j];

                    tey_yzz_xxyy_0[j] = pa_y[j] * tey_zz_xxyy_0[j] - pc_y[j] * tey_zz_xxyy_1[j] + fl1_fx * tey_zz_xxy_0[j] - fl1_fx * tey_zz_xxy_1[j] + ta_zz_xxyy_1[j];

                    tez_yzz_xxyy_0[j] = pa_y[j] * tez_zz_xxyy_0[j] - pc_y[j] * tez_zz_xxyy_1[j] + fl1_fx * tez_zz_xxy_0[j] - fl1_fx * tez_zz_xxy_1[j];

                    tex_yzz_xxyz_0[j] = pa_y[j] * tex_zz_xxyz_0[j] - pc_y[j] * tex_zz_xxyz_1[j] + 0.5 * fl1_fx * tex_zz_xxz_0[j] - 0.5 * fl1_fx * tex_zz_xxz_1[j];

                    tey_yzz_xxyz_0[j] = pa_y[j] * tey_zz_xxyz_0[j] - pc_y[j] * tey_zz_xxyz_1[j] + 0.5 * fl1_fx * tey_zz_xxz_0[j] - 0.5 * fl1_fx * tey_zz_xxz_1[j] + ta_zz_xxyz_1[j];

                    tez_yzz_xxyz_0[j] = pa_y[j] * tez_zz_xxyz_0[j] - pc_y[j] * tez_zz_xxyz_1[j] + 0.5 * fl1_fx * tez_zz_xxz_0[j] - 0.5 * fl1_fx * tez_zz_xxz_1[j];

                    tex_yzz_xxzz_0[j] = pa_y[j] * tex_zz_xxzz_0[j] - pc_y[j] * tex_zz_xxzz_1[j];

                    tey_yzz_xxzz_0[j] = pa_y[j] * tey_zz_xxzz_0[j] - pc_y[j] * tey_zz_xxzz_1[j] + ta_zz_xxzz_1[j];

                    tez_yzz_xxzz_0[j] = pa_y[j] * tez_zz_xxzz_0[j] - pc_y[j] * tez_zz_xxzz_1[j];

                    tex_yzz_xyyy_0[j] = pa_y[j] * tex_zz_xyyy_0[j] - pc_y[j] * tex_zz_xyyy_1[j] + 1.5 * fl1_fx * tex_zz_xyy_0[j] - 1.5 * fl1_fx * tex_zz_xyy_1[j];

                    tey_yzz_xyyy_0[j] = pa_y[j] * tey_zz_xyyy_0[j] - pc_y[j] * tey_zz_xyyy_1[j] + 1.5 * fl1_fx * tey_zz_xyy_0[j] - 1.5 * fl1_fx * tey_zz_xyy_1[j] + ta_zz_xyyy_1[j];

                    tez_yzz_xyyy_0[j] = pa_y[j] * tez_zz_xyyy_0[j] - pc_y[j] * tez_zz_xyyy_1[j] + 1.5 * fl1_fx * tez_zz_xyy_0[j] - 1.5 * fl1_fx * tez_zz_xyy_1[j];

                    tex_yzz_xyyz_0[j] = pa_y[j] * tex_zz_xyyz_0[j] - pc_y[j] * tex_zz_xyyz_1[j] + fl1_fx * tex_zz_xyz_0[j] - fl1_fx * tex_zz_xyz_1[j];

                    tey_yzz_xyyz_0[j] = pa_y[j] * tey_zz_xyyz_0[j] - pc_y[j] * tey_zz_xyyz_1[j] + fl1_fx * tey_zz_xyz_0[j] - fl1_fx * tey_zz_xyz_1[j] + ta_zz_xyyz_1[j];

                    tez_yzz_xyyz_0[j] = pa_y[j] * tez_zz_xyyz_0[j] - pc_y[j] * tez_zz_xyyz_1[j] + fl1_fx * tez_zz_xyz_0[j] - fl1_fx * tez_zz_xyz_1[j];

                    tex_yzz_xyzz_0[j] = pa_y[j] * tex_zz_xyzz_0[j] - pc_y[j] * tex_zz_xyzz_1[j] + 0.5 * fl1_fx * tex_zz_xzz_0[j] - 0.5 * fl1_fx * tex_zz_xzz_1[j];

                    tey_yzz_xyzz_0[j] = pa_y[j] * tey_zz_xyzz_0[j] - pc_y[j] * tey_zz_xyzz_1[j] + 0.5 * fl1_fx * tey_zz_xzz_0[j] - 0.5 * fl1_fx * tey_zz_xzz_1[j] + ta_zz_xyzz_1[j];

                    tez_yzz_xyzz_0[j] = pa_y[j] * tez_zz_xyzz_0[j] - pc_y[j] * tez_zz_xyzz_1[j] + 0.5 * fl1_fx * tez_zz_xzz_0[j] - 0.5 * fl1_fx * tez_zz_xzz_1[j];

                    tex_yzz_xzzz_0[j] = pa_y[j] * tex_zz_xzzz_0[j] - pc_y[j] * tex_zz_xzzz_1[j];

                    tey_yzz_xzzz_0[j] = pa_y[j] * tey_zz_xzzz_0[j] - pc_y[j] * tey_zz_xzzz_1[j] + ta_zz_xzzz_1[j];

                    tez_yzz_xzzz_0[j] = pa_y[j] * tez_zz_xzzz_0[j] - pc_y[j] * tez_zz_xzzz_1[j];

                    tex_yzz_yyyy_0[j] = pa_y[j] * tex_zz_yyyy_0[j] - pc_y[j] * tex_zz_yyyy_1[j] + 2.0 * fl1_fx * tex_zz_yyy_0[j] - 2.0 * fl1_fx * tex_zz_yyy_1[j];

                    tey_yzz_yyyy_0[j] = pa_y[j] * tey_zz_yyyy_0[j] - pc_y[j] * tey_zz_yyyy_1[j] + 2.0 * fl1_fx * tey_zz_yyy_0[j] - 2.0 * fl1_fx * tey_zz_yyy_1[j] + ta_zz_yyyy_1[j];

                    tez_yzz_yyyy_0[j] = pa_y[j] * tez_zz_yyyy_0[j] - pc_y[j] * tez_zz_yyyy_1[j] + 2.0 * fl1_fx * tez_zz_yyy_0[j] - 2.0 * fl1_fx * tez_zz_yyy_1[j];

                    tex_yzz_yyyz_0[j] = pa_y[j] * tex_zz_yyyz_0[j] - pc_y[j] * tex_zz_yyyz_1[j] + 1.5 * fl1_fx * tex_zz_yyz_0[j] - 1.5 * fl1_fx * tex_zz_yyz_1[j];

                    tey_yzz_yyyz_0[j] = pa_y[j] * tey_zz_yyyz_0[j] - pc_y[j] * tey_zz_yyyz_1[j] + 1.5 * fl1_fx * tey_zz_yyz_0[j] - 1.5 * fl1_fx * tey_zz_yyz_1[j] + ta_zz_yyyz_1[j];

                    tez_yzz_yyyz_0[j] = pa_y[j] * tez_zz_yyyz_0[j] - pc_y[j] * tez_zz_yyyz_1[j] + 1.5 * fl1_fx * tez_zz_yyz_0[j] - 1.5 * fl1_fx * tez_zz_yyz_1[j];

                    tex_yzz_yyzz_0[j] = pa_y[j] * tex_zz_yyzz_0[j] - pc_y[j] * tex_zz_yyzz_1[j] + fl1_fx * tex_zz_yzz_0[j] - fl1_fx * tex_zz_yzz_1[j];

                    tey_yzz_yyzz_0[j] = pa_y[j] * tey_zz_yyzz_0[j] - pc_y[j] * tey_zz_yyzz_1[j] + fl1_fx * tey_zz_yzz_0[j] - fl1_fx * tey_zz_yzz_1[j] + ta_zz_yyzz_1[j];

                    tez_yzz_yyzz_0[j] = pa_y[j] * tez_zz_yyzz_0[j] - pc_y[j] * tez_zz_yyzz_1[j] + fl1_fx * tez_zz_yzz_0[j] - fl1_fx * tez_zz_yzz_1[j];

                    tex_yzz_yzzz_0[j] = pa_y[j] * tex_zz_yzzz_0[j] - pc_y[j] * tex_zz_yzzz_1[j] + 0.5 * fl1_fx * tex_zz_zzz_0[j] - 0.5 * fl1_fx * tex_zz_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFG_400_450(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (400,450)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
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

                auto tex_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 75); 

                auto tey_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 75); 

                auto tez_zz_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 75); 

                auto tex_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 76); 

                auto tey_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 76); 

                auto tez_zz_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 76); 

                auto tex_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 77); 

                auto tey_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 77); 

                auto tez_zz_xxxz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 77); 

                auto tex_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 78); 

                auto tey_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 78); 

                auto tez_zz_xxyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 78); 

                auto tex_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 79); 

                auto tey_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 79); 

                auto tez_zz_xxyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 79); 

                auto tex_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 80); 

                auto tey_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 80); 

                auto tez_zz_xxzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 80); 

                auto tex_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 81); 

                auto tey_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 81); 

                auto tez_zz_xyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 81); 

                auto tex_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 82); 

                auto tey_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 82); 

                auto tez_zz_xyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 82); 

                auto tex_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 83); 

                auto tey_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 83); 

                auto tez_zz_xyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 83); 

                auto tex_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 84); 

                auto tey_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 84); 

                auto tez_zz_xzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 84); 

                auto tex_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 85); 

                auto tey_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 85); 

                auto tez_zz_yyyy_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 85); 

                auto tex_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 86); 

                auto tey_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 86); 

                auto tez_zz_yyyz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 86); 

                auto tex_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 87); 

                auto tey_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 87); 

                auto tez_zz_yyzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 87); 

                auto tex_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 88); 

                auto tey_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 88); 

                auto tez_zz_yzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 88); 

                auto tex_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 89); 

                auto tey_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 89); 

                auto tez_zz_zzzz_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 89); 

                auto tex_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 75); 

                auto tey_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 75); 

                auto tez_zz_xxxx_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 75); 

                auto tex_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 76); 

                auto tey_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 76); 

                auto tez_zz_xxxy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 76); 

                auto tex_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 77); 

                auto tey_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 77); 

                auto tez_zz_xxxz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 77); 

                auto tex_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 78); 

                auto tey_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 78); 

                auto tez_zz_xxyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 78); 

                auto tex_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 79); 

                auto tey_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 79); 

                auto tez_zz_xxyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 79); 

                auto tex_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 80); 

                auto tey_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 80); 

                auto tez_zz_xxzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 80); 

                auto tex_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 81); 

                auto tey_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 81); 

                auto tez_zz_xyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 81); 

                auto tex_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 82); 

                auto tey_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 82); 

                auto tez_zz_xyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 82); 

                auto tex_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 83); 

                auto tey_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 83); 

                auto tez_zz_xyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 83); 

                auto tex_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 84); 

                auto tey_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 84); 

                auto tez_zz_xzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 84); 

                auto tex_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 85); 

                auto tey_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 85); 

                auto tez_zz_yyyy_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 85); 

                auto tex_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 86); 

                auto tey_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 86); 

                auto tez_zz_yyyz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 86); 

                auto tex_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 87); 

                auto tey_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 87); 

                auto tez_zz_yyzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 87); 

                auto tex_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 88); 

                auto tey_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 88); 

                auto tez_zz_yzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 88); 

                auto tex_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * idx + 89); 

                auto tey_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 90 * bdim + 90 * idx + 89); 

                auto tez_zz_zzzz_1 = primBuffer.data(pidx_e_2_4_m1 + 180 * bdim + 90 * idx + 89); 

                auto tex_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 30); 

                auto tey_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 30); 

                auto tez_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 30); 

                auto tex_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 31); 

                auto tey_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 31); 

                auto tez_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 31); 

                auto tex_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 32); 

                auto tey_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 32); 

                auto tez_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 32); 

                auto tex_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 33); 

                auto tey_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 33); 

                auto tez_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 33); 

                auto tex_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 34); 

                auto tey_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 34); 

                auto tez_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 34); 

                auto tex_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 35); 

                auto tey_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 35); 

                auto tez_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 35); 

                auto tex_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 36); 

                auto tey_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 36); 

                auto tez_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 36); 

                auto tex_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 37); 

                auto tey_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 37); 

                auto tez_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 37); 

                auto tex_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 38); 

                auto tey_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 38); 

                auto tez_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 38); 

                auto tex_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 39); 

                auto tey_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 39); 

                auto tez_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 39); 

                auto tex_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 40); 

                auto tey_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 40); 

                auto tez_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 40); 

                auto tex_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 41); 

                auto tey_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 41); 

                auto tez_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 41); 

                auto tex_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 42); 

                auto tey_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 42); 

                auto tez_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 42); 

                auto tex_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 43); 

                auto tey_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 43); 

                auto tez_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 43); 

                auto tex_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 44); 

                auto tey_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 44); 

                auto tez_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 44); 

                auto tex_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 30); 

                auto tey_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 30); 

                auto tez_z_xxxx_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 30); 

                auto tex_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 31); 

                auto tey_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 31); 

                auto tez_z_xxxy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 31); 

                auto tex_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 32); 

                auto tey_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 32); 

                auto tez_z_xxxz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 32); 

                auto tex_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 33); 

                auto tey_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 33); 

                auto tez_z_xxyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 33); 

                auto tex_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 34); 

                auto tey_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 34); 

                auto tez_z_xxyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 34); 

                auto tex_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 35); 

                auto tey_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 35); 

                auto tez_z_xxzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 35); 

                auto tex_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 36); 

                auto tey_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 36); 

                auto tez_z_xyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 36); 

                auto tex_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 37); 

                auto tey_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 37); 

                auto tez_z_xyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 37); 

                auto tex_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 38); 

                auto tey_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 38); 

                auto tez_z_xyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 38); 

                auto tex_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 39); 

                auto tey_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 39); 

                auto tez_z_xzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 39); 

                auto tex_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 40); 

                auto tey_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 40); 

                auto tez_z_yyyy_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 40); 

                auto tex_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 41); 

                auto tey_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 41); 

                auto tez_z_yyyz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 41); 

                auto tex_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 42); 

                auto tey_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 42); 

                auto tez_z_yyzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 42); 

                auto tex_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 43); 

                auto tey_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 43); 

                auto tez_z_yzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 43); 

                auto tex_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * idx + 44); 

                auto tey_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 45 * bdim + 45 * idx + 44); 

                auto tez_z_zzzz_1 = primBuffer.data(pidx_e_1_4_m1 + 90 * bdim + 45 * idx + 44); 

                auto tex_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 50); 

                auto tey_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 50); 

                auto tez_zz_xxx_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 50); 

                auto tex_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 51); 

                auto tey_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 51); 

                auto tez_zz_xxy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 51); 

                auto tex_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 52); 

                auto tey_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 52); 

                auto tez_zz_xxz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 52); 

                auto tex_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 53); 

                auto tey_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 53); 

                auto tez_zz_xyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 53); 

                auto tex_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 54); 

                auto tey_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 54); 

                auto tez_zz_xyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 54); 

                auto tex_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 55); 

                auto tey_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 55); 

                auto tez_zz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 55); 

                auto tex_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 56); 

                auto tey_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 56); 

                auto tez_zz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 56); 

                auto tex_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 57); 

                auto tey_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 57); 

                auto tez_zz_yyz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 57); 

                auto tex_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 58); 

                auto tey_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 58); 

                auto tez_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 58); 

                auto tex_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 59); 

                auto tey_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 59); 

                auto tez_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 59); 

                auto tex_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 50); 

                auto tey_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 50); 

                auto tez_zz_xxx_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 50); 

                auto tex_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 51); 

                auto tey_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 51); 

                auto tez_zz_xxy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 51); 

                auto tex_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 52); 

                auto tey_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 52); 

                auto tez_zz_xxz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 52); 

                auto tex_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 53); 

                auto tey_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 53); 

                auto tez_zz_xyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 53); 

                auto tex_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 54); 

                auto tey_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 54); 

                auto tez_zz_xyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 54); 

                auto tex_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 55); 

                auto tey_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 55); 

                auto tez_zz_xzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 55); 

                auto tex_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 56); 

                auto tey_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 56); 

                auto tez_zz_yyy_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 56); 

                auto tex_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 57); 

                auto tey_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 57); 

                auto tez_zz_yyz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 57); 

                auto tex_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 58); 

                auto tey_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 58); 

                auto tez_zz_yzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 58); 

                auto tex_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * idx + 59); 

                auto tey_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 60 * bdim + 60 * idx + 59); 

                auto tez_zz_zzz_1 = primBuffer.data(pidx_e_2_3_m1 + 120 * bdim + 60 * idx + 59); 

                auto ta_zz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 75); 

                auto ta_zz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 76); 

                auto ta_zz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 77); 

                auto ta_zz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 78); 

                auto ta_zz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 79); 

                auto ta_zz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 80); 

                auto ta_zz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 81); 

                auto ta_zz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 82); 

                auto ta_zz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 83); 

                auto ta_zz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 84); 

                auto ta_zz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 85); 

                auto ta_zz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 86); 

                auto ta_zz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 87); 

                auto ta_zz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 88); 

                auto ta_zz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 89); 

                // set up pointers to integrals

                auto tey_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 133); 

                auto tez_yzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 133); 

                auto tex_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 134); 

                auto tey_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 134); 

                auto tez_yzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 134); 

                auto tex_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 135); 

                auto tey_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 135); 

                auto tez_zzz_xxxx_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 135); 

                auto tex_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 136); 

                auto tey_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 136); 

                auto tez_zzz_xxxy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 136); 

                auto tex_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 137); 

                auto tey_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 137); 

                auto tez_zzz_xxxz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 137); 

                auto tex_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 138); 

                auto tey_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 138); 

                auto tez_zzz_xxyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 138); 

                auto tex_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 139); 

                auto tey_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 139); 

                auto tez_zzz_xxyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 139); 

                auto tex_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 140); 

                auto tey_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 140); 

                auto tez_zzz_xxzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 140); 

                auto tex_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 141); 

                auto tey_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 141); 

                auto tez_zzz_xyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 141); 

                auto tex_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 142); 

                auto tey_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 142); 

                auto tez_zzz_xyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 142); 

                auto tex_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 143); 

                auto tey_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 143); 

                auto tez_zzz_xyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 143); 

                auto tex_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 144); 

                auto tey_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 144); 

                auto tez_zzz_xzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 144); 

                auto tex_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 145); 

                auto tey_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 145); 

                auto tez_zzz_yyyy_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 145); 

                auto tex_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 146); 

                auto tey_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 146); 

                auto tez_zzz_yyyz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 146); 

                auto tex_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 147); 

                auto tey_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 147); 

                auto tez_zzz_yyzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 147); 

                auto tex_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 148); 

                auto tey_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 148); 

                auto tez_zzz_yzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 148); 

                auto tex_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * idx + 149); 

                auto tey_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 150 * bdim + 150 * idx + 149); 

                auto tez_zzz_zzzz_0 = primBuffer.data(pidx_e_3_4_m0 + 300 * bdim + 150 * idx + 149); 

                // Batch of Integrals (400,450)

                #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_zz_xxxx_1, ta_zz_xxxy_1, ta_zz_xxxz_1, \
                                         ta_zz_xxyy_1, ta_zz_xxyz_1, ta_zz_xxzz_1, ta_zz_xyyy_1, ta_zz_xyyz_1, ta_zz_xyzz_1, \
                                         ta_zz_xzzz_1, ta_zz_yyyy_1, ta_zz_yyyz_1, ta_zz_yyzz_1, ta_zz_yzzz_1, ta_zz_zzzz_1, \
                                         tex_yzz_zzzz_0, tex_z_xxxx_0, tex_z_xxxx_1, tex_z_xxxy_0, tex_z_xxxy_1, tex_z_xxxz_0, \
                                         tex_z_xxxz_1, tex_z_xxyy_0, tex_z_xxyy_1, tex_z_xxyz_0, tex_z_xxyz_1, tex_z_xxzz_0, \
                                         tex_z_xxzz_1, tex_z_xyyy_0, tex_z_xyyy_1, tex_z_xyyz_0, tex_z_xyyz_1, tex_z_xyzz_0, \
                                         tex_z_xyzz_1, tex_z_xzzz_0, tex_z_xzzz_1, tex_z_yyyy_0, tex_z_yyyy_1, tex_z_yyyz_0, \
                                         tex_z_yyyz_1, tex_z_yyzz_0, tex_z_yyzz_1, tex_z_yzzz_0, tex_z_yzzz_1, tex_z_zzzz_0, \
                                         tex_z_zzzz_1, tex_zz_xxx_0, tex_zz_xxx_1, tex_zz_xxxx_0, tex_zz_xxxx_1, \
                                         tex_zz_xxxy_0, tex_zz_xxxy_1, tex_zz_xxxz_0, tex_zz_xxxz_1, tex_zz_xxy_0, \
                                         tex_zz_xxy_1, tex_zz_xxyy_0, tex_zz_xxyy_1, tex_zz_xxyz_0, tex_zz_xxyz_1, \
                                         tex_zz_xxz_0, tex_zz_xxz_1, tex_zz_xxzz_0, tex_zz_xxzz_1, tex_zz_xyy_0, \
                                         tex_zz_xyy_1, tex_zz_xyyy_0, tex_zz_xyyy_1, tex_zz_xyyz_0, tex_zz_xyyz_1, \
                                         tex_zz_xyz_0, tex_zz_xyz_1, tex_zz_xyzz_0, tex_zz_xyzz_1, tex_zz_xzz_0, \
                                         tex_zz_xzz_1, tex_zz_xzzz_0, tex_zz_xzzz_1, tex_zz_yyy_0, tex_zz_yyy_1, \
                                         tex_zz_yyyy_0, tex_zz_yyyy_1, tex_zz_yyyz_0, tex_zz_yyyz_1, tex_zz_yyz_0, \
                                         tex_zz_yyz_1, tex_zz_yyzz_0, tex_zz_yyzz_1, tex_zz_yzz_0, tex_zz_yzz_1, \
                                         tex_zz_yzzz_0, tex_zz_yzzz_1, tex_zz_zzz_0, tex_zz_zzz_1, tex_zz_zzzz_0, \
                                         tex_zz_zzzz_1, tex_zzz_xxxx_0, tex_zzz_xxxy_0, tex_zzz_xxxz_0, tex_zzz_xxyy_0, \
                                         tex_zzz_xxyz_0, tex_zzz_xxzz_0, tex_zzz_xyyy_0, tex_zzz_xyyz_0, tex_zzz_xyzz_0, \
                                         tex_zzz_xzzz_0, tex_zzz_yyyy_0, tex_zzz_yyyz_0, tex_zzz_yyzz_0, tex_zzz_yzzz_0, \
                                         tex_zzz_zzzz_0, tey_yzz_yzzz_0, tey_yzz_zzzz_0, tey_z_xxxx_0, tey_z_xxxx_1, \
                                         tey_z_xxxy_0, tey_z_xxxy_1, tey_z_xxxz_0, tey_z_xxxz_1, tey_z_xxyy_0, tey_z_xxyy_1, \
                                         tey_z_xxyz_0, tey_z_xxyz_1, tey_z_xxzz_0, tey_z_xxzz_1, tey_z_xyyy_0, tey_z_xyyy_1, \
                                         tey_z_xyyz_0, tey_z_xyyz_1, tey_z_xyzz_0, tey_z_xyzz_1, tey_z_xzzz_0, tey_z_xzzz_1, \
                                         tey_z_yyyy_0, tey_z_yyyy_1, tey_z_yyyz_0, tey_z_yyyz_1, tey_z_yyzz_0, tey_z_yyzz_1, \
                                         tey_z_yzzz_0, tey_z_yzzz_1, tey_z_zzzz_0, tey_z_zzzz_1, tey_zz_xxx_0, tey_zz_xxx_1, \
                                         tey_zz_xxxx_0, tey_zz_xxxx_1, tey_zz_xxxy_0, tey_zz_xxxy_1, tey_zz_xxxz_0, \
                                         tey_zz_xxxz_1, tey_zz_xxy_0, tey_zz_xxy_1, tey_zz_xxyy_0, tey_zz_xxyy_1, \
                                         tey_zz_xxyz_0, tey_zz_xxyz_1, tey_zz_xxz_0, tey_zz_xxz_1, tey_zz_xxzz_0, \
                                         tey_zz_xxzz_1, tey_zz_xyy_0, tey_zz_xyy_1, tey_zz_xyyy_0, tey_zz_xyyy_1, \
                                         tey_zz_xyyz_0, tey_zz_xyyz_1, tey_zz_xyz_0, tey_zz_xyz_1, tey_zz_xyzz_0, \
                                         tey_zz_xyzz_1, tey_zz_xzz_0, tey_zz_xzz_1, tey_zz_xzzz_0, tey_zz_xzzz_1, \
                                         tey_zz_yyy_0, tey_zz_yyy_1, tey_zz_yyyy_0, tey_zz_yyyy_1, tey_zz_yyyz_0, \
                                         tey_zz_yyyz_1, tey_zz_yyz_0, tey_zz_yyz_1, tey_zz_yyzz_0, tey_zz_yyzz_1, \
                                         tey_zz_yzz_0, tey_zz_yzz_1, tey_zz_yzzz_0, tey_zz_yzzz_1, tey_zz_zzz_0, \
                                         tey_zz_zzz_1, tey_zz_zzzz_0, tey_zz_zzzz_1, tey_zzz_xxxx_0, tey_zzz_xxxy_0, \
                                         tey_zzz_xxxz_0, tey_zzz_xxyy_0, tey_zzz_xxyz_0, tey_zzz_xxzz_0, tey_zzz_xyyy_0, \
                                         tey_zzz_xyyz_0, tey_zzz_xyzz_0, tey_zzz_xzzz_0, tey_zzz_yyyy_0, tey_zzz_yyyz_0, \
                                         tey_zzz_yyzz_0, tey_zzz_yzzz_0, tey_zzz_zzzz_0, tez_yzz_yzzz_0, tez_yzz_zzzz_0, \
                                         tez_z_xxxx_0, tez_z_xxxx_1, tez_z_xxxy_0, tez_z_xxxy_1, tez_z_xxxz_0, tez_z_xxxz_1, \
                                         tez_z_xxyy_0, tez_z_xxyy_1, tez_z_xxyz_0, tez_z_xxyz_1, tez_z_xxzz_0, tez_z_xxzz_1, \
                                         tez_z_xyyy_0, tez_z_xyyy_1, tez_z_xyyz_0, tez_z_xyyz_1, tez_z_xyzz_0, tez_z_xyzz_1, \
                                         tez_z_xzzz_0, tez_z_xzzz_1, tez_z_yyyy_0, tez_z_yyyy_1, tez_z_yyyz_0, tez_z_yyyz_1, \
                                         tez_z_yyzz_0, tez_z_yyzz_1, tez_z_yzzz_0, tez_z_yzzz_1, tez_z_zzzz_0, tez_z_zzzz_1, \
                                         tez_zz_xxx_0, tez_zz_xxx_1, tez_zz_xxxx_0, tez_zz_xxxx_1, tez_zz_xxxy_0, \
                                         tez_zz_xxxy_1, tez_zz_xxxz_0, tez_zz_xxxz_1, tez_zz_xxy_0, tez_zz_xxy_1, \
                                         tez_zz_xxyy_0, tez_zz_xxyy_1, tez_zz_xxyz_0, tez_zz_xxyz_1, tez_zz_xxz_0, \
                                         tez_zz_xxz_1, tez_zz_xxzz_0, tez_zz_xxzz_1, tez_zz_xyy_0, tez_zz_xyy_1, \
                                         tez_zz_xyyy_0, tez_zz_xyyy_1, tez_zz_xyyz_0, tez_zz_xyyz_1, tez_zz_xyz_0, \
                                         tez_zz_xyz_1, tez_zz_xyzz_0, tez_zz_xyzz_1, tez_zz_xzz_0, tez_zz_xzz_1, \
                                         tez_zz_xzzz_0, tez_zz_xzzz_1, tez_zz_yyy_0, tez_zz_yyy_1, tez_zz_yyyy_0, \
                                         tez_zz_yyyy_1, tez_zz_yyyz_0, tez_zz_yyyz_1, tez_zz_yyz_0, tez_zz_yyz_1, \
                                         tez_zz_yyzz_0, tez_zz_yyzz_1, tez_zz_yzz_0, tez_zz_yzz_1, tez_zz_yzzz_0, \
                                         tez_zz_yzzz_1, tez_zz_zzz_0, tez_zz_zzz_1, tez_zz_zzzz_0, tez_zz_zzzz_1, \
                                         tez_zzz_xxxx_0, tez_zzz_xxxy_0, tez_zzz_xxxz_0, tez_zzz_xxyy_0, tez_zzz_xxyz_0, \
                                         tez_zzz_xxzz_0, tez_zzz_xyyy_0, tez_zzz_xyyz_0, tez_zzz_xyzz_0, tez_zzz_xzzz_0, \
                                         tez_zzz_yyyy_0, tez_zzz_yyyz_0, tez_zzz_yyzz_0, tez_zzz_yzzz_0, tez_zzz_zzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tey_yzz_yzzz_0[j] = pa_y[j] * tey_zz_yzzz_0[j] - pc_y[j] * tey_zz_yzzz_1[j] + 0.5 * fl1_fx * tey_zz_zzz_0[j] - 0.5 * fl1_fx * tey_zz_zzz_1[j] + ta_zz_yzzz_1[j];

                    tez_yzz_yzzz_0[j] = pa_y[j] * tez_zz_yzzz_0[j] - pc_y[j] * tez_zz_yzzz_1[j] + 0.5 * fl1_fx * tez_zz_zzz_0[j] - 0.5 * fl1_fx * tez_zz_zzz_1[j];

                    tex_yzz_zzzz_0[j] = pa_y[j] * tex_zz_zzzz_0[j] - pc_y[j] * tex_zz_zzzz_1[j];

                    tey_yzz_zzzz_0[j] = pa_y[j] * tey_zz_zzzz_0[j] - pc_y[j] * tey_zz_zzzz_1[j] + ta_zz_zzzz_1[j];

                    tez_yzz_zzzz_0[j] = pa_y[j] * tez_zz_zzzz_0[j] - pc_y[j] * tez_zz_zzzz_1[j];

                    tex_zzz_xxxx_0[j] = pa_z[j] * tex_zz_xxxx_0[j] - pc_z[j] * tex_zz_xxxx_1[j] + fl1_fx * tex_z_xxxx_0[j] - fl1_fx * tex_z_xxxx_1[j];

                    tey_zzz_xxxx_0[j] = pa_z[j] * tey_zz_xxxx_0[j] - pc_z[j] * tey_zz_xxxx_1[j] + fl1_fx * tey_z_xxxx_0[j] - fl1_fx * tey_z_xxxx_1[j];

                    tez_zzz_xxxx_0[j] = pa_z[j] * tez_zz_xxxx_0[j] - pc_z[j] * tez_zz_xxxx_1[j] + fl1_fx * tez_z_xxxx_0[j] - fl1_fx * tez_z_xxxx_1[j] + ta_zz_xxxx_1[j];

                    tex_zzz_xxxy_0[j] = pa_z[j] * tex_zz_xxxy_0[j] - pc_z[j] * tex_zz_xxxy_1[j] + fl1_fx * tex_z_xxxy_0[j] - fl1_fx * tex_z_xxxy_1[j];

                    tey_zzz_xxxy_0[j] = pa_z[j] * tey_zz_xxxy_0[j] - pc_z[j] * tey_zz_xxxy_1[j] + fl1_fx * tey_z_xxxy_0[j] - fl1_fx * tey_z_xxxy_1[j];

                    tez_zzz_xxxy_0[j] = pa_z[j] * tez_zz_xxxy_0[j] - pc_z[j] * tez_zz_xxxy_1[j] + fl1_fx * tez_z_xxxy_0[j] - fl1_fx * tez_z_xxxy_1[j] + ta_zz_xxxy_1[j];

                    tex_zzz_xxxz_0[j] = pa_z[j] * tex_zz_xxxz_0[j] - pc_z[j] * tex_zz_xxxz_1[j] + fl1_fx * tex_z_xxxz_0[j] - fl1_fx * tex_z_xxxz_1[j] + 0.5 * fl1_fx * tex_zz_xxx_0[j] - 0.5 * fl1_fx * tex_zz_xxx_1[j];

                    tey_zzz_xxxz_0[j] = pa_z[j] * tey_zz_xxxz_0[j] - pc_z[j] * tey_zz_xxxz_1[j] + fl1_fx * tey_z_xxxz_0[j] - fl1_fx * tey_z_xxxz_1[j] + 0.5 * fl1_fx * tey_zz_xxx_0[j] - 0.5 * fl1_fx * tey_zz_xxx_1[j];

                    tez_zzz_xxxz_0[j] = pa_z[j] * tez_zz_xxxz_0[j] - pc_z[j] * tez_zz_xxxz_1[j] + fl1_fx * tez_z_xxxz_0[j] - fl1_fx * tez_z_xxxz_1[j] + 0.5 * fl1_fx * tez_zz_xxx_0[j] - 0.5 * fl1_fx * tez_zz_xxx_1[j] + ta_zz_xxxz_1[j];

                    tex_zzz_xxyy_0[j] = pa_z[j] * tex_zz_xxyy_0[j] - pc_z[j] * tex_zz_xxyy_1[j] + fl1_fx * tex_z_xxyy_0[j] - fl1_fx * tex_z_xxyy_1[j];

                    tey_zzz_xxyy_0[j] = pa_z[j] * tey_zz_xxyy_0[j] - pc_z[j] * tey_zz_xxyy_1[j] + fl1_fx * tey_z_xxyy_0[j] - fl1_fx * tey_z_xxyy_1[j];

                    tez_zzz_xxyy_0[j] = pa_z[j] * tez_zz_xxyy_0[j] - pc_z[j] * tez_zz_xxyy_1[j] + fl1_fx * tez_z_xxyy_0[j] - fl1_fx * tez_z_xxyy_1[j] + ta_zz_xxyy_1[j];

                    tex_zzz_xxyz_0[j] = pa_z[j] * tex_zz_xxyz_0[j] - pc_z[j] * tex_zz_xxyz_1[j] + fl1_fx * tex_z_xxyz_0[j] - fl1_fx * tex_z_xxyz_1[j] + 0.5 * fl1_fx * tex_zz_xxy_0[j] - 0.5 * fl1_fx * tex_zz_xxy_1[j];

                    tey_zzz_xxyz_0[j] = pa_z[j] * tey_zz_xxyz_0[j] - pc_z[j] * tey_zz_xxyz_1[j] + fl1_fx * tey_z_xxyz_0[j] - fl1_fx * tey_z_xxyz_1[j] + 0.5 * fl1_fx * tey_zz_xxy_0[j] - 0.5 * fl1_fx * tey_zz_xxy_1[j];

                    tez_zzz_xxyz_0[j] = pa_z[j] * tez_zz_xxyz_0[j] - pc_z[j] * tez_zz_xxyz_1[j] + fl1_fx * tez_z_xxyz_0[j] - fl1_fx * tez_z_xxyz_1[j] + 0.5 * fl1_fx * tez_zz_xxy_0[j] - 0.5 * fl1_fx * tez_zz_xxy_1[j] + ta_zz_xxyz_1[j];

                    tex_zzz_xxzz_0[j] = pa_z[j] * tex_zz_xxzz_0[j] - pc_z[j] * tex_zz_xxzz_1[j] + fl1_fx * tex_z_xxzz_0[j] - fl1_fx * tex_z_xxzz_1[j] + fl1_fx * tex_zz_xxz_0[j] - fl1_fx * tex_zz_xxz_1[j];

                    tey_zzz_xxzz_0[j] = pa_z[j] * tey_zz_xxzz_0[j] - pc_z[j] * tey_zz_xxzz_1[j] + fl1_fx * tey_z_xxzz_0[j] - fl1_fx * tey_z_xxzz_1[j] + fl1_fx * tey_zz_xxz_0[j] - fl1_fx * tey_zz_xxz_1[j];

                    tez_zzz_xxzz_0[j] = pa_z[j] * tez_zz_xxzz_0[j] - pc_z[j] * tez_zz_xxzz_1[j] + fl1_fx * tez_z_xxzz_0[j] - fl1_fx * tez_z_xxzz_1[j] + fl1_fx * tez_zz_xxz_0[j] - fl1_fx * tez_zz_xxz_1[j] + ta_zz_xxzz_1[j];

                    tex_zzz_xyyy_0[j] = pa_z[j] * tex_zz_xyyy_0[j] - pc_z[j] * tex_zz_xyyy_1[j] + fl1_fx * tex_z_xyyy_0[j] - fl1_fx * tex_z_xyyy_1[j];

                    tey_zzz_xyyy_0[j] = pa_z[j] * tey_zz_xyyy_0[j] - pc_z[j] * tey_zz_xyyy_1[j] + fl1_fx * tey_z_xyyy_0[j] - fl1_fx * tey_z_xyyy_1[j];

                    tez_zzz_xyyy_0[j] = pa_z[j] * tez_zz_xyyy_0[j] - pc_z[j] * tez_zz_xyyy_1[j] + fl1_fx * tez_z_xyyy_0[j] - fl1_fx * tez_z_xyyy_1[j] + ta_zz_xyyy_1[j];

                    tex_zzz_xyyz_0[j] = pa_z[j] * tex_zz_xyyz_0[j] - pc_z[j] * tex_zz_xyyz_1[j] + fl1_fx * tex_z_xyyz_0[j] - fl1_fx * tex_z_xyyz_1[j] + 0.5 * fl1_fx * tex_zz_xyy_0[j] - 0.5 * fl1_fx * tex_zz_xyy_1[j];

                    tey_zzz_xyyz_0[j] = pa_z[j] * tey_zz_xyyz_0[j] - pc_z[j] * tey_zz_xyyz_1[j] + fl1_fx * tey_z_xyyz_0[j] - fl1_fx * tey_z_xyyz_1[j] + 0.5 * fl1_fx * tey_zz_xyy_0[j] - 0.5 * fl1_fx * tey_zz_xyy_1[j];

                    tez_zzz_xyyz_0[j] = pa_z[j] * tez_zz_xyyz_0[j] - pc_z[j] * tez_zz_xyyz_1[j] + fl1_fx * tez_z_xyyz_0[j] - fl1_fx * tez_z_xyyz_1[j] + 0.5 * fl1_fx * tez_zz_xyy_0[j] - 0.5 * fl1_fx * tez_zz_xyy_1[j] + ta_zz_xyyz_1[j];

                    tex_zzz_xyzz_0[j] = pa_z[j] * tex_zz_xyzz_0[j] - pc_z[j] * tex_zz_xyzz_1[j] + fl1_fx * tex_z_xyzz_0[j] - fl1_fx * tex_z_xyzz_1[j] + fl1_fx * tex_zz_xyz_0[j] - fl1_fx * tex_zz_xyz_1[j];

                    tey_zzz_xyzz_0[j] = pa_z[j] * tey_zz_xyzz_0[j] - pc_z[j] * tey_zz_xyzz_1[j] + fl1_fx * tey_z_xyzz_0[j] - fl1_fx * tey_z_xyzz_1[j] + fl1_fx * tey_zz_xyz_0[j] - fl1_fx * tey_zz_xyz_1[j];

                    tez_zzz_xyzz_0[j] = pa_z[j] * tez_zz_xyzz_0[j] - pc_z[j] * tez_zz_xyzz_1[j] + fl1_fx * tez_z_xyzz_0[j] - fl1_fx * tez_z_xyzz_1[j] + fl1_fx * tez_zz_xyz_0[j] - fl1_fx * tez_zz_xyz_1[j] + ta_zz_xyzz_1[j];

                    tex_zzz_xzzz_0[j] = pa_z[j] * tex_zz_xzzz_0[j] - pc_z[j] * tex_zz_xzzz_1[j] + fl1_fx * tex_z_xzzz_0[j] - fl1_fx * tex_z_xzzz_1[j] + 1.5 * fl1_fx * tex_zz_xzz_0[j] - 1.5 * fl1_fx * tex_zz_xzz_1[j];

                    tey_zzz_xzzz_0[j] = pa_z[j] * tey_zz_xzzz_0[j] - pc_z[j] * tey_zz_xzzz_1[j] + fl1_fx * tey_z_xzzz_0[j] - fl1_fx * tey_z_xzzz_1[j] + 1.5 * fl1_fx * tey_zz_xzz_0[j] - 1.5 * fl1_fx * tey_zz_xzz_1[j];

                    tez_zzz_xzzz_0[j] = pa_z[j] * tez_zz_xzzz_0[j] - pc_z[j] * tez_zz_xzzz_1[j] + fl1_fx * tez_z_xzzz_0[j] - fl1_fx * tez_z_xzzz_1[j] + 1.5 * fl1_fx * tez_zz_xzz_0[j] - 1.5 * fl1_fx * tez_zz_xzz_1[j] + ta_zz_xzzz_1[j];

                    tex_zzz_yyyy_0[j] = pa_z[j] * tex_zz_yyyy_0[j] - pc_z[j] * tex_zz_yyyy_1[j] + fl1_fx * tex_z_yyyy_0[j] - fl1_fx * tex_z_yyyy_1[j];

                    tey_zzz_yyyy_0[j] = pa_z[j] * tey_zz_yyyy_0[j] - pc_z[j] * tey_zz_yyyy_1[j] + fl1_fx * tey_z_yyyy_0[j] - fl1_fx * tey_z_yyyy_1[j];

                    tez_zzz_yyyy_0[j] = pa_z[j] * tez_zz_yyyy_0[j] - pc_z[j] * tez_zz_yyyy_1[j] + fl1_fx * tez_z_yyyy_0[j] - fl1_fx * tez_z_yyyy_1[j] + ta_zz_yyyy_1[j];

                    tex_zzz_yyyz_0[j] = pa_z[j] * tex_zz_yyyz_0[j] - pc_z[j] * tex_zz_yyyz_1[j] + fl1_fx * tex_z_yyyz_0[j] - fl1_fx * tex_z_yyyz_1[j] + 0.5 * fl1_fx * tex_zz_yyy_0[j] - 0.5 * fl1_fx * tex_zz_yyy_1[j];

                    tey_zzz_yyyz_0[j] = pa_z[j] * tey_zz_yyyz_0[j] - pc_z[j] * tey_zz_yyyz_1[j] + fl1_fx * tey_z_yyyz_0[j] - fl1_fx * tey_z_yyyz_1[j] + 0.5 * fl1_fx * tey_zz_yyy_0[j] - 0.5 * fl1_fx * tey_zz_yyy_1[j];

                    tez_zzz_yyyz_0[j] = pa_z[j] * tez_zz_yyyz_0[j] - pc_z[j] * tez_zz_yyyz_1[j] + fl1_fx * tez_z_yyyz_0[j] - fl1_fx * tez_z_yyyz_1[j] + 0.5 * fl1_fx * tez_zz_yyy_0[j] - 0.5 * fl1_fx * tez_zz_yyy_1[j] + ta_zz_yyyz_1[j];

                    tex_zzz_yyzz_0[j] = pa_z[j] * tex_zz_yyzz_0[j] - pc_z[j] * tex_zz_yyzz_1[j] + fl1_fx * tex_z_yyzz_0[j] - fl1_fx * tex_z_yyzz_1[j] + fl1_fx * tex_zz_yyz_0[j] - fl1_fx * tex_zz_yyz_1[j];

                    tey_zzz_yyzz_0[j] = pa_z[j] * tey_zz_yyzz_0[j] - pc_z[j] * tey_zz_yyzz_1[j] + fl1_fx * tey_z_yyzz_0[j] - fl1_fx * tey_z_yyzz_1[j] + fl1_fx * tey_zz_yyz_0[j] - fl1_fx * tey_zz_yyz_1[j];

                    tez_zzz_yyzz_0[j] = pa_z[j] * tez_zz_yyzz_0[j] - pc_z[j] * tez_zz_yyzz_1[j] + fl1_fx * tez_z_yyzz_0[j] - fl1_fx * tez_z_yyzz_1[j] + fl1_fx * tez_zz_yyz_0[j] - fl1_fx * tez_zz_yyz_1[j] + ta_zz_yyzz_1[j];

                    tex_zzz_yzzz_0[j] = pa_z[j] * tex_zz_yzzz_0[j] - pc_z[j] * tex_zz_yzzz_1[j] + fl1_fx * tex_z_yzzz_0[j] - fl1_fx * tex_z_yzzz_1[j] + 1.5 * fl1_fx * tex_zz_yzz_0[j] - 1.5 * fl1_fx * tex_zz_yzz_1[j];

                    tey_zzz_yzzz_0[j] = pa_z[j] * tey_zz_yzzz_0[j] - pc_z[j] * tey_zz_yzzz_1[j] + fl1_fx * tey_z_yzzz_0[j] - fl1_fx * tey_z_yzzz_1[j] + 1.5 * fl1_fx * tey_zz_yzz_0[j] - 1.5 * fl1_fx * tey_zz_yzz_1[j];

                    tez_zzz_yzzz_0[j] = pa_z[j] * tez_zz_yzzz_0[j] - pc_z[j] * tez_zz_yzzz_1[j] + fl1_fx * tez_z_yzzz_0[j] - fl1_fx * tez_z_yzzz_1[j] + 1.5 * fl1_fx * tez_zz_yzz_0[j] - 1.5 * fl1_fx * tez_zz_yzz_1[j] + ta_zz_yzzz_1[j];

                    tex_zzz_zzzz_0[j] = pa_z[j] * tex_zz_zzzz_0[j] - pc_z[j] * tex_zz_zzzz_1[j] + fl1_fx * tex_z_zzzz_0[j] - fl1_fx * tex_z_zzzz_1[j] + 2.0 * fl1_fx * tex_zz_zzz_0[j] - 2.0 * fl1_fx * tex_zz_zzz_1[j];

                    tey_zzz_zzzz_0[j] = pa_z[j] * tey_zz_zzzz_0[j] - pc_z[j] * tey_zz_zzzz_1[j] + fl1_fx * tey_z_zzzz_0[j] - fl1_fx * tey_z_zzzz_1[j] + 2.0 * fl1_fx * tey_zz_zzz_0[j] - 2.0 * fl1_fx * tey_zz_zzz_1[j];

                    tez_zzz_zzzz_0[j] = pa_z[j] * tez_zz_zzzz_0[j] - pc_z[j] * tez_zz_zzzz_1[j] + fl1_fx * tez_z_zzzz_0[j] - fl1_fx * tez_z_zzzz_1[j] + 2.0 * fl1_fx * tez_zz_zzz_0[j] - 2.0 * fl1_fx * tez_zz_zzz_1[j] + ta_zz_zzzz_1[j];
                }

                idx++;
            }
        }
    }


} // efieldrecfunc namespace

