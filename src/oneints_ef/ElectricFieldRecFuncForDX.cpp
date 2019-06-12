//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldRecFuncForDX.hpp"

namespace efieldrecfunc { // efieldrecfunc namespace

    void
    compElectricFieldForDD(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pcDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        efieldrecfunc::compElectricFieldForDD_0_36(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   pcDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        efieldrecfunc::compElectricFieldForDD_36_72(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    pcDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        efieldrecfunc::compElectricFieldForDD_72_108(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     pcDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 
    }

    void
    compElectricFieldForDD_0_36(      CMemBlock2D<double>& primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
    {
        // Batch of Integrals (0,36)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tex_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx); 

                auto tey_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx); 

                auto tez_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx); 

                auto tex_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 1); 

                auto tey_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 1); 

                auto tez_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 1); 

                auto tex_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 2); 

                auto tey_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 2); 

                auto tez_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 2); 

                auto tex_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 3); 

                auto tey_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 3); 

                auto tez_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 3); 

                auto tex_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 4); 

                auto tey_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 4); 

                auto tez_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 4); 

                auto tex_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 5); 

                auto tey_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 5); 

                auto tez_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 5); 

                auto tex_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 6); 

                auto tey_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 7); 

                auto tey_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 8); 

                auto tey_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 9); 

                auto tey_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 10); 

                auto tey_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 11); 

                auto tey_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx); 

                auto tey_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx); 

                auto tez_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx); 

                auto tex_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 1); 

                auto tey_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 1); 

                auto tez_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 1); 

                auto tex_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 2); 

                auto tey_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 2); 

                auto tez_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 2); 

                auto tex_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 3); 

                auto tey_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 3); 

                auto tez_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 3); 

                auto tex_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 4); 

                auto tey_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 4); 

                auto tez_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 4); 

                auto tex_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 5); 

                auto tey_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 5); 

                auto tez_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 5); 

                auto tex_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 6); 

                auto tey_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 7); 

                auto tey_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 8); 

                auto tey_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 9); 

                auto tey_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 10); 

                auto tey_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 11); 

                auto tey_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx); 

                auto tey_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx); 

                auto tez_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx); 

                auto tex_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 1); 

                auto tey_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 1); 

                auto tez_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 1); 

                auto tex_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 2); 

                auto tey_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 2); 

                auto tez_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 2); 

                auto tex_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 3); 

                auto tey_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 3); 

                auto tez_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 3); 

                auto tex_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 4); 

                auto tey_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 4); 

                auto tez_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 4); 

                auto tex_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 5); 

                auto tey_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 5); 

                auto tez_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 5); 

                auto tex_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx); 

                auto tey_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx); 

                auto tez_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx); 

                auto tex_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 1); 

                auto tey_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 1); 

                auto tez_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 1); 

                auto tex_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 2); 

                auto tey_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 2); 

                auto tez_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 2); 

                auto tex_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 3); 

                auto tey_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 3); 

                auto tez_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 3); 

                auto tex_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 4); 

                auto tey_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 4); 

                auto tez_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 4); 

                auto tex_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 5); 

                auto tey_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 5); 

                auto tez_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 5); 

                auto tex_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx); 

                auto tey_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx); 

                auto tez_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx); 

                auto tex_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 1); 

                auto tey_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 1); 

                auto tez_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 1); 

                auto tex_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 2); 

                auto tey_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 2); 

                auto tez_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 2); 

                auto tex_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 3); 

                auto tey_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 3); 

                auto tez_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 3); 

                auto tex_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 4); 

                auto tey_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 4); 

                auto tez_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 4); 

                auto tex_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 5); 

                auto tey_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 5); 

                auto tez_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 5); 

                auto tex_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx); 

                auto tey_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx); 

                auto tez_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx); 

                auto tex_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 1); 

                auto tey_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 1); 

                auto tez_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 1); 

                auto tex_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 2); 

                auto tey_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 2); 

                auto tez_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 2); 

                auto tex_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 3); 

                auto tey_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 3); 

                auto tez_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 3); 

                auto tex_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 4); 

                auto tey_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 4); 

                auto tez_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 4); 

                auto tex_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 5); 

                auto tey_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 5); 

                auto tez_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 5); 

                auto ta_x_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx); 

                auto ta_x_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 1); 

                auto ta_x_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 2); 

                auto ta_x_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 3); 

                auto ta_x_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 4); 

                auto ta_x_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 5); 

                auto ta_y_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 6); 

                auto ta_y_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 7); 

                auto ta_y_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 8); 

                auto ta_y_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 9); 

                auto ta_y_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 10); 

                auto ta_y_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 11); 

                // set up pointers to integrals

                auto tex_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx); 

                auto tey_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx); 

                auto tez_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx); 

                auto tex_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 1); 

                auto tey_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 1); 

                auto tez_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 1); 

                auto tex_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 2); 

                auto tey_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 2); 

                auto tez_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 2); 

                auto tex_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 3); 

                auto tey_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 3); 

                auto tez_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 3); 

                auto tex_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 4); 

                auto tey_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 4); 

                auto tez_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 4); 

                auto tex_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 5); 

                auto tey_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 5); 

                auto tez_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 5); 

                auto tex_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 6); 

                auto tey_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 6); 

                auto tez_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 6); 

                auto tex_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 7); 

                auto tey_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 7); 

                auto tez_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 7); 

                auto tex_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 8); 

                auto tey_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 8); 

                auto tez_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 8); 

                auto tex_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 9); 

                auto tey_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 9); 

                auto tez_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 9); 

                auto tex_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 10); 

                auto tey_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 10); 

                auto tez_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 10); 

                auto tex_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 11); 

                auto tey_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 11); 

                auto tez_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 11); 

                // Batch of Integrals (0,36)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_x_xx_1, ta_x_xy_1, ta_x_xz_1, ta_x_yy_1, ta_x_yz_1, \
                                         ta_x_zz_1, ta_y_xx_1, ta_y_xy_1, ta_y_xz_1, ta_y_yy_1, ta_y_yz_1, ta_y_zz_1, \
                                         tex_0_xx_0, tex_0_xx_1, tex_0_xy_0, tex_0_xy_1, tex_0_xz_0, tex_0_xz_1, tex_0_yy_0, \
                                         tex_0_yy_1, tex_0_yz_0, tex_0_yz_1, tex_0_zz_0, tex_0_zz_1, tex_x_x_0, tex_x_x_1, \
                                         tex_x_xx_0, tex_x_xx_1, tex_x_xy_0, tex_x_xy_1, tex_x_xz_0, tex_x_xz_1, tex_x_y_0, \
                                         tex_x_y_1, tex_x_yy_0, tex_x_yy_1, tex_x_yz_0, tex_x_yz_1, tex_x_z_0, tex_x_z_1, \
                                         tex_x_zz_0, tex_x_zz_1, tex_xx_xx_0, tex_xx_xy_0, tex_xx_xz_0, tex_xx_yy_0, \
                                         tex_xx_yz_0, tex_xx_zz_0, tex_xy_xx_0, tex_xy_xy_0, tex_xy_xz_0, tex_xy_yy_0, \
                                         tex_xy_yz_0, tex_xy_zz_0, tex_y_x_0, tex_y_x_1, tex_y_xx_0, tex_y_xx_1, tex_y_xy_0, \
                                         tex_y_xy_1, tex_y_xz_0, tex_y_xz_1, tex_y_y_0, tex_y_y_1, tex_y_yy_0, tex_y_yy_1, \
                                         tex_y_yz_0, tex_y_yz_1, tex_y_z_0, tex_y_z_1, tex_y_zz_0, tex_y_zz_1, tey_0_xx_0, \
                                         tey_0_xx_1, tey_0_xy_0, tey_0_xy_1, tey_0_xz_0, tey_0_xz_1, tey_0_yy_0, tey_0_yy_1, \
                                         tey_0_yz_0, tey_0_yz_1, tey_0_zz_0, tey_0_zz_1, tey_x_x_0, tey_x_x_1, tey_x_xx_0, \
                                         tey_x_xx_1, tey_x_xy_0, tey_x_xy_1, tey_x_xz_0, tey_x_xz_1, tey_x_y_0, tey_x_y_1, \
                                         tey_x_yy_0, tey_x_yy_1, tey_x_yz_0, tey_x_yz_1, tey_x_z_0, tey_x_z_1, tey_x_zz_0, \
                                         tey_x_zz_1, tey_xx_xx_0, tey_xx_xy_0, tey_xx_xz_0, tey_xx_yy_0, tey_xx_yz_0, \
                                         tey_xx_zz_0, tey_xy_xx_0, tey_xy_xy_0, tey_xy_xz_0, tey_xy_yy_0, tey_xy_yz_0, \
                                         tey_xy_zz_0, tey_y_x_0, tey_y_x_1, tey_y_xx_0, tey_y_xx_1, tey_y_xy_0, tey_y_xy_1, \
                                         tey_y_xz_0, tey_y_xz_1, tey_y_y_0, tey_y_y_1, tey_y_yy_0, tey_y_yy_1, tey_y_yz_0, \
                                         tey_y_yz_1, tey_y_z_0, tey_y_z_1, tey_y_zz_0, tey_y_zz_1, tez_0_xx_0, tez_0_xx_1, \
                                         tez_0_xy_0, tez_0_xy_1, tez_0_xz_0, tez_0_xz_1, tez_0_yy_0, tez_0_yy_1, tez_0_yz_0, \
                                         tez_0_yz_1, tez_0_zz_0, tez_0_zz_1, tez_x_x_0, tez_x_x_1, tez_x_xx_0, tez_x_xx_1, \
                                         tez_x_xy_0, tez_x_xy_1, tez_x_xz_0, tez_x_xz_1, tez_x_y_0, tez_x_y_1, tez_x_yy_0, \
                                         tez_x_yy_1, tez_x_yz_0, tez_x_yz_1, tez_x_z_0, tez_x_z_1, tez_x_zz_0, tez_x_zz_1, \
                                         tez_xx_xx_0, tez_xx_xy_0, tez_xx_xz_0, tez_xx_yy_0, tez_xx_yz_0, tez_xx_zz_0, \
                                         tez_xy_xx_0, tez_xy_xy_0, tez_xy_xz_0, tez_xy_yy_0, tez_xy_yz_0, tez_xy_zz_0, \
                                         tez_y_x_0, tez_y_x_1, tez_y_xx_0, tez_y_xx_1, tez_y_xy_0, tez_y_xy_1, tez_y_xz_0, \
                                         tez_y_xz_1, tez_y_y_0, tez_y_y_1, tez_y_yy_0, tez_y_yy_1, tez_y_yz_0, tez_y_yz_1, \
                                         tez_y_z_0, tez_y_z_1, tez_y_zz_0, tez_y_zz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xx_xx_0[j] = pa_x[j] * tex_x_xx_0[j] - pc_x[j] * tex_x_xx_1[j] + 0.5 * fl1_fx * tex_0_xx_0[j] - 0.5 * fl1_fx * tex_0_xx_1[j] + fl1_fx * tex_x_x_0[j] - fl1_fx * tex_x_x_1[j] + ta_x_xx_1[j];

                    tey_xx_xx_0[j] = pa_x[j] * tey_x_xx_0[j] - pc_x[j] * tey_x_xx_1[j] + 0.5 * fl1_fx * tey_0_xx_0[j] - 0.5 * fl1_fx * tey_0_xx_1[j] + fl1_fx * tey_x_x_0[j] - fl1_fx * tey_x_x_1[j];

                    tez_xx_xx_0[j] = pa_x[j] * tez_x_xx_0[j] - pc_x[j] * tez_x_xx_1[j] + 0.5 * fl1_fx * tez_0_xx_0[j] - 0.5 * fl1_fx * tez_0_xx_1[j] + fl1_fx * tez_x_x_0[j] - fl1_fx * tez_x_x_1[j];

                    tex_xx_xy_0[j] = pa_x[j] * tex_x_xy_0[j] - pc_x[j] * tex_x_xy_1[j] + 0.5 * fl1_fx * tex_0_xy_0[j] - 0.5 * fl1_fx * tex_0_xy_1[j] + 0.5 * fl1_fx * tex_x_y_0[j] - 0.5 * fl1_fx * tex_x_y_1[j] + ta_x_xy_1[j];

                    tey_xx_xy_0[j] = pa_x[j] * tey_x_xy_0[j] - pc_x[j] * tey_x_xy_1[j] + 0.5 * fl1_fx * tey_0_xy_0[j] - 0.5 * fl1_fx * tey_0_xy_1[j] + 0.5 * fl1_fx * tey_x_y_0[j] - 0.5 * fl1_fx * tey_x_y_1[j];

                    tez_xx_xy_0[j] = pa_x[j] * tez_x_xy_0[j] - pc_x[j] * tez_x_xy_1[j] + 0.5 * fl1_fx * tez_0_xy_0[j] - 0.5 * fl1_fx * tez_0_xy_1[j] + 0.5 * fl1_fx * tez_x_y_0[j] - 0.5 * fl1_fx * tez_x_y_1[j];

                    tex_xx_xz_0[j] = pa_x[j] * tex_x_xz_0[j] - pc_x[j] * tex_x_xz_1[j] + 0.5 * fl1_fx * tex_0_xz_0[j] - 0.5 * fl1_fx * tex_0_xz_1[j] + 0.5 * fl1_fx * tex_x_z_0[j] - 0.5 * fl1_fx * tex_x_z_1[j] + ta_x_xz_1[j];

                    tey_xx_xz_0[j] = pa_x[j] * tey_x_xz_0[j] - pc_x[j] * tey_x_xz_1[j] + 0.5 * fl1_fx * tey_0_xz_0[j] - 0.5 * fl1_fx * tey_0_xz_1[j] + 0.5 * fl1_fx * tey_x_z_0[j] - 0.5 * fl1_fx * tey_x_z_1[j];

                    tez_xx_xz_0[j] = pa_x[j] * tez_x_xz_0[j] - pc_x[j] * tez_x_xz_1[j] + 0.5 * fl1_fx * tez_0_xz_0[j] - 0.5 * fl1_fx * tez_0_xz_1[j] + 0.5 * fl1_fx * tez_x_z_0[j] - 0.5 * fl1_fx * tez_x_z_1[j];

                    tex_xx_yy_0[j] = pa_x[j] * tex_x_yy_0[j] - pc_x[j] * tex_x_yy_1[j] + 0.5 * fl1_fx * tex_0_yy_0[j] - 0.5 * fl1_fx * tex_0_yy_1[j] + ta_x_yy_1[j];

                    tey_xx_yy_0[j] = pa_x[j] * tey_x_yy_0[j] - pc_x[j] * tey_x_yy_1[j] + 0.5 * fl1_fx * tey_0_yy_0[j] - 0.5 * fl1_fx * tey_0_yy_1[j];

                    tez_xx_yy_0[j] = pa_x[j] * tez_x_yy_0[j] - pc_x[j] * tez_x_yy_1[j] + 0.5 * fl1_fx * tez_0_yy_0[j] - 0.5 * fl1_fx * tez_0_yy_1[j];

                    tex_xx_yz_0[j] = pa_x[j] * tex_x_yz_0[j] - pc_x[j] * tex_x_yz_1[j] + 0.5 * fl1_fx * tex_0_yz_0[j] - 0.5 * fl1_fx * tex_0_yz_1[j] + ta_x_yz_1[j];

                    tey_xx_yz_0[j] = pa_x[j] * tey_x_yz_0[j] - pc_x[j] * tey_x_yz_1[j] + 0.5 * fl1_fx * tey_0_yz_0[j] - 0.5 * fl1_fx * tey_0_yz_1[j];

                    tez_xx_yz_0[j] = pa_x[j] * tez_x_yz_0[j] - pc_x[j] * tez_x_yz_1[j] + 0.5 * fl1_fx * tez_0_yz_0[j] - 0.5 * fl1_fx * tez_0_yz_1[j];

                    tex_xx_zz_0[j] = pa_x[j] * tex_x_zz_0[j] - pc_x[j] * tex_x_zz_1[j] + 0.5 * fl1_fx * tex_0_zz_0[j] - 0.5 * fl1_fx * tex_0_zz_1[j] + ta_x_zz_1[j];

                    tey_xx_zz_0[j] = pa_x[j] * tey_x_zz_0[j] - pc_x[j] * tey_x_zz_1[j] + 0.5 * fl1_fx * tey_0_zz_0[j] - 0.5 * fl1_fx * tey_0_zz_1[j];

                    tez_xx_zz_0[j] = pa_x[j] * tez_x_zz_0[j] - pc_x[j] * tez_x_zz_1[j] + 0.5 * fl1_fx * tez_0_zz_0[j] - 0.5 * fl1_fx * tez_0_zz_1[j];

                    tex_xy_xx_0[j] = pa_x[j] * tex_y_xx_0[j] - pc_x[j] * tex_y_xx_1[j] + fl1_fx * tex_y_x_0[j] - fl1_fx * tex_y_x_1[j] + ta_y_xx_1[j];

                    tey_xy_xx_0[j] = pa_x[j] * tey_y_xx_0[j] - pc_x[j] * tey_y_xx_1[j] + fl1_fx * tey_y_x_0[j] - fl1_fx * tey_y_x_1[j];

                    tez_xy_xx_0[j] = pa_x[j] * tez_y_xx_0[j] - pc_x[j] * tez_y_xx_1[j] + fl1_fx * tez_y_x_0[j] - fl1_fx * tez_y_x_1[j];

                    tex_xy_xy_0[j] = pa_x[j] * tex_y_xy_0[j] - pc_x[j] * tex_y_xy_1[j] + 0.5 * fl1_fx * tex_y_y_0[j] - 0.5 * fl1_fx * tex_y_y_1[j] + ta_y_xy_1[j];

                    tey_xy_xy_0[j] = pa_x[j] * tey_y_xy_0[j] - pc_x[j] * tey_y_xy_1[j] + 0.5 * fl1_fx * tey_y_y_0[j] - 0.5 * fl1_fx * tey_y_y_1[j];

                    tez_xy_xy_0[j] = pa_x[j] * tez_y_xy_0[j] - pc_x[j] * tez_y_xy_1[j] + 0.5 * fl1_fx * tez_y_y_0[j] - 0.5 * fl1_fx * tez_y_y_1[j];

                    tex_xy_xz_0[j] = pa_x[j] * tex_y_xz_0[j] - pc_x[j] * tex_y_xz_1[j] + 0.5 * fl1_fx * tex_y_z_0[j] - 0.5 * fl1_fx * tex_y_z_1[j] + ta_y_xz_1[j];

                    tey_xy_xz_0[j] = pa_x[j] * tey_y_xz_0[j] - pc_x[j] * tey_y_xz_1[j] + 0.5 * fl1_fx * tey_y_z_0[j] - 0.5 * fl1_fx * tey_y_z_1[j];

                    tez_xy_xz_0[j] = pa_x[j] * tez_y_xz_0[j] - pc_x[j] * tez_y_xz_1[j] + 0.5 * fl1_fx * tez_y_z_0[j] - 0.5 * fl1_fx * tez_y_z_1[j];

                    tex_xy_yy_0[j] = pa_x[j] * tex_y_yy_0[j] - pc_x[j] * tex_y_yy_1[j] + ta_y_yy_1[j];

                    tey_xy_yy_0[j] = pa_x[j] * tey_y_yy_0[j] - pc_x[j] * tey_y_yy_1[j];

                    tez_xy_yy_0[j] = pa_x[j] * tez_y_yy_0[j] - pc_x[j] * tez_y_yy_1[j];

                    tex_xy_yz_0[j] = pa_x[j] * tex_y_yz_0[j] - pc_x[j] * tex_y_yz_1[j] + ta_y_yz_1[j];

                    tey_xy_yz_0[j] = pa_x[j] * tey_y_yz_0[j] - pc_x[j] * tey_y_yz_1[j];

                    tez_xy_yz_0[j] = pa_x[j] * tez_y_yz_0[j] - pc_x[j] * tez_y_yz_1[j];

                    tex_xy_zz_0[j] = pa_x[j] * tex_y_zz_0[j] - pc_x[j] * tex_y_zz_1[j] + ta_y_zz_1[j];

                    tey_xy_zz_0[j] = pa_x[j] * tey_y_zz_0[j] - pc_x[j] * tey_y_zz_1[j];

                    tez_xy_zz_0[j] = pa_x[j] * tez_y_zz_0[j] - pc_x[j] * tez_y_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDD_36_72(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (36,72)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tex_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 6); 

                auto tey_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 7); 

                auto tey_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 8); 

                auto tey_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 9); 

                auto tey_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 10); 

                auto tey_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 11); 

                auto tey_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 15); 

                auto tey_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 16); 

                auto tey_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 17); 

                auto tey_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 6); 

                auto tey_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 7); 

                auto tey_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 8); 

                auto tey_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 9); 

                auto tey_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 10); 

                auto tey_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 11); 

                auto tey_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 15); 

                auto tey_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 16); 

                auto tey_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 17); 

                auto tey_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 17); 

                auto tex_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx); 

                auto tey_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx); 

                auto tez_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx); 

                auto tex_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 1); 

                auto tey_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 1); 

                auto tez_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 1); 

                auto tex_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 2); 

                auto tey_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 2); 

                auto tez_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 2); 

                auto tex_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 3); 

                auto tey_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 3); 

                auto tez_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 3); 

                auto tex_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 4); 

                auto tey_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 4); 

                auto tez_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 4); 

                auto tex_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 5); 

                auto tey_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 5); 

                auto tez_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 5); 

                auto tex_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx); 

                auto tey_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx); 

                auto tez_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx); 

                auto tex_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 1); 

                auto tey_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 1); 

                auto tez_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 1); 

                auto tex_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 2); 

                auto tey_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 2); 

                auto tez_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 2); 

                auto tex_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 3); 

                auto tey_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 3); 

                auto tez_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 3); 

                auto tex_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 4); 

                auto tey_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 4); 

                auto tez_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 4); 

                auto tex_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 5); 

                auto tey_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 5); 

                auto tez_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 5); 

                auto tex_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 3); 

                auto tey_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 3); 

                auto tez_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 3); 

                auto tex_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 4); 

                auto tey_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 4); 

                auto tez_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 4); 

                auto tex_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 5); 

                auto tey_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 5); 

                auto tez_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 5); 

                auto tex_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 6); 

                auto tey_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 6); 

                auto tez_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 6); 

                auto tex_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 7); 

                auto tey_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 7); 

                auto tez_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 7); 

                auto tex_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 8); 

                auto tey_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 8); 

                auto tez_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 8); 

                auto tex_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 3); 

                auto tey_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 3); 

                auto tez_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 3); 

                auto tex_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 4); 

                auto tey_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 4); 

                auto tez_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 4); 

                auto tex_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 5); 

                auto tey_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 5); 

                auto tez_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 5); 

                auto tex_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 6); 

                auto tey_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 6); 

                auto tez_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 6); 

                auto tex_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 7); 

                auto tey_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 7); 

                auto tez_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 7); 

                auto tex_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 8); 

                auto tey_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 8); 

                auto tez_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 8); 

                auto ta_y_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 6); 

                auto ta_y_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 7); 

                auto ta_y_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 8); 

                auto ta_y_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 9); 

                auto ta_y_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 10); 

                auto ta_y_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 11); 

                auto ta_z_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 12); 

                auto ta_z_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 13); 

                auto ta_z_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 14); 

                auto ta_z_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 15); 

                auto ta_z_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 16); 

                auto ta_z_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 17); 

                // set up pointers to integrals

                auto tex_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 12); 

                auto tey_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 12); 

                auto tez_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 12); 

                auto tex_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 13); 

                auto tey_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 13); 

                auto tez_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 13); 

                auto tex_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 14); 

                auto tey_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 14); 

                auto tez_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 14); 

                auto tex_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 15); 

                auto tey_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 15); 

                auto tez_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 15); 

                auto tex_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 16); 

                auto tey_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 16); 

                auto tez_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 16); 

                auto tex_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 17); 

                auto tey_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 17); 

                auto tez_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 17); 

                auto tex_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 18); 

                auto tey_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 19); 

                auto tey_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 20); 

                auto tey_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 21); 

                auto tey_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 22); 

                auto tey_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 23); 

                auto tey_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 23); 

                // Batch of Integrals (36,72)

                #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_y_xx_1, ta_y_xy_1, ta_y_xz_1, ta_y_yy_1, \
                                         ta_y_yz_1, ta_y_zz_1, ta_z_xx_1, ta_z_xy_1, ta_z_xz_1, ta_z_yy_1, ta_z_yz_1, \
                                         ta_z_zz_1, tex_0_xx_0, tex_0_xx_1, tex_0_xy_0, tex_0_xy_1, tex_0_xz_0, tex_0_xz_1, \
                                         tex_0_yy_0, tex_0_yy_1, tex_0_yz_0, tex_0_yz_1, tex_0_zz_0, tex_0_zz_1, \
                                         tex_xz_xx_0, tex_xz_xy_0, tex_xz_xz_0, tex_xz_yy_0, tex_xz_yz_0, tex_xz_zz_0, \
                                         tex_y_x_0, tex_y_x_1, tex_y_xx_0, tex_y_xx_1, tex_y_xy_0, tex_y_xy_1, tex_y_xz_0, \
                                         tex_y_xz_1, tex_y_y_0, tex_y_y_1, tex_y_yy_0, tex_y_yy_1, tex_y_yz_0, tex_y_yz_1, \
                                         tex_y_z_0, tex_y_z_1, tex_y_zz_0, tex_y_zz_1, tex_yy_xx_0, tex_yy_xy_0, \
                                         tex_yy_xz_0, tex_yy_yy_0, tex_yy_yz_0, tex_yy_zz_0, tex_z_x_0, tex_z_x_1, \
                                         tex_z_xx_0, tex_z_xx_1, tex_z_xy_0, tex_z_xy_1, tex_z_xz_0, tex_z_xz_1, tex_z_y_0, \
                                         tex_z_y_1, tex_z_yy_0, tex_z_yy_1, tex_z_yz_0, tex_z_yz_1, tex_z_z_0, tex_z_z_1, \
                                         tex_z_zz_0, tex_z_zz_1, tey_0_xx_0, tey_0_xx_1, tey_0_xy_0, tey_0_xy_1, tey_0_xz_0, \
                                         tey_0_xz_1, tey_0_yy_0, tey_0_yy_1, tey_0_yz_0, tey_0_yz_1, tey_0_zz_0, tey_0_zz_1, \
                                         tey_xz_xx_0, tey_xz_xy_0, tey_xz_xz_0, tey_xz_yy_0, tey_xz_yz_0, tey_xz_zz_0, \
                                         tey_y_x_0, tey_y_x_1, tey_y_xx_0, tey_y_xx_1, tey_y_xy_0, tey_y_xy_1, tey_y_xz_0, \
                                         tey_y_xz_1, tey_y_y_0, tey_y_y_1, tey_y_yy_0, tey_y_yy_1, tey_y_yz_0, tey_y_yz_1, \
                                         tey_y_z_0, tey_y_z_1, tey_y_zz_0, tey_y_zz_1, tey_yy_xx_0, tey_yy_xy_0, \
                                         tey_yy_xz_0, tey_yy_yy_0, tey_yy_yz_0, tey_yy_zz_0, tey_z_x_0, tey_z_x_1, \
                                         tey_z_xx_0, tey_z_xx_1, tey_z_xy_0, tey_z_xy_1, tey_z_xz_0, tey_z_xz_1, tey_z_y_0, \
                                         tey_z_y_1, tey_z_yy_0, tey_z_yy_1, tey_z_yz_0, tey_z_yz_1, tey_z_z_0, tey_z_z_1, \
                                         tey_z_zz_0, tey_z_zz_1, tez_0_xx_0, tez_0_xx_1, tez_0_xy_0, tez_0_xy_1, tez_0_xz_0, \
                                         tez_0_xz_1, tez_0_yy_0, tez_0_yy_1, tez_0_yz_0, tez_0_yz_1, tez_0_zz_0, tez_0_zz_1, \
                                         tez_xz_xx_0, tez_xz_xy_0, tez_xz_xz_0, tez_xz_yy_0, tez_xz_yz_0, tez_xz_zz_0, \
                                         tez_y_x_0, tez_y_x_1, tez_y_xx_0, tez_y_xx_1, tez_y_xy_0, tez_y_xy_1, tez_y_xz_0, \
                                         tez_y_xz_1, tez_y_y_0, tez_y_y_1, tez_y_yy_0, tez_y_yy_1, tez_y_yz_0, tez_y_yz_1, \
                                         tez_y_z_0, tez_y_z_1, tez_y_zz_0, tez_y_zz_1, tez_yy_xx_0, tez_yy_xy_0, \
                                         tez_yy_xz_0, tez_yy_yy_0, tez_yy_yz_0, tez_yy_zz_0, tez_z_x_0, tez_z_x_1, \
                                         tez_z_xx_0, tez_z_xx_1, tez_z_xy_0, tez_z_xy_1, tez_z_xz_0, tez_z_xz_1, tez_z_y_0, \
                                         tez_z_y_1, tez_z_yy_0, tez_z_yy_1, tez_z_yz_0, tez_z_yz_1, tez_z_z_0, tez_z_z_1, \
                                         tez_z_zz_0, tez_z_zz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xz_xx_0[j] = pa_x[j] * tex_z_xx_0[j] - pc_x[j] * tex_z_xx_1[j] + fl1_fx * tex_z_x_0[j] - fl1_fx * tex_z_x_1[j] + ta_z_xx_1[j];

                    tey_xz_xx_0[j] = pa_x[j] * tey_z_xx_0[j] - pc_x[j] * tey_z_xx_1[j] + fl1_fx * tey_z_x_0[j] - fl1_fx * tey_z_x_1[j];

                    tez_xz_xx_0[j] = pa_x[j] * tez_z_xx_0[j] - pc_x[j] * tez_z_xx_1[j] + fl1_fx * tez_z_x_0[j] - fl1_fx * tez_z_x_1[j];

                    tex_xz_xy_0[j] = pa_x[j] * tex_z_xy_0[j] - pc_x[j] * tex_z_xy_1[j] + 0.5 * fl1_fx * tex_z_y_0[j] - 0.5 * fl1_fx * tex_z_y_1[j] + ta_z_xy_1[j];

                    tey_xz_xy_0[j] = pa_x[j] * tey_z_xy_0[j] - pc_x[j] * tey_z_xy_1[j] + 0.5 * fl1_fx * tey_z_y_0[j] - 0.5 * fl1_fx * tey_z_y_1[j];

                    tez_xz_xy_0[j] = pa_x[j] * tez_z_xy_0[j] - pc_x[j] * tez_z_xy_1[j] + 0.5 * fl1_fx * tez_z_y_0[j] - 0.5 * fl1_fx * tez_z_y_1[j];

                    tex_xz_xz_0[j] = pa_x[j] * tex_z_xz_0[j] - pc_x[j] * tex_z_xz_1[j] + 0.5 * fl1_fx * tex_z_z_0[j] - 0.5 * fl1_fx * tex_z_z_1[j] + ta_z_xz_1[j];

                    tey_xz_xz_0[j] = pa_x[j] * tey_z_xz_0[j] - pc_x[j] * tey_z_xz_1[j] + 0.5 * fl1_fx * tey_z_z_0[j] - 0.5 * fl1_fx * tey_z_z_1[j];

                    tez_xz_xz_0[j] = pa_x[j] * tez_z_xz_0[j] - pc_x[j] * tez_z_xz_1[j] + 0.5 * fl1_fx * tez_z_z_0[j] - 0.5 * fl1_fx * tez_z_z_1[j];

                    tex_xz_yy_0[j] = pa_x[j] * tex_z_yy_0[j] - pc_x[j] * tex_z_yy_1[j] + ta_z_yy_1[j];

                    tey_xz_yy_0[j] = pa_x[j] * tey_z_yy_0[j] - pc_x[j] * tey_z_yy_1[j];

                    tez_xz_yy_0[j] = pa_x[j] * tez_z_yy_0[j] - pc_x[j] * tez_z_yy_1[j];

                    tex_xz_yz_0[j] = pa_x[j] * tex_z_yz_0[j] - pc_x[j] * tex_z_yz_1[j] + ta_z_yz_1[j];

                    tey_xz_yz_0[j] = pa_x[j] * tey_z_yz_0[j] - pc_x[j] * tey_z_yz_1[j];

                    tez_xz_yz_0[j] = pa_x[j] * tez_z_yz_0[j] - pc_x[j] * tez_z_yz_1[j];

                    tex_xz_zz_0[j] = pa_x[j] * tex_z_zz_0[j] - pc_x[j] * tex_z_zz_1[j] + ta_z_zz_1[j];

                    tey_xz_zz_0[j] = pa_x[j] * tey_z_zz_0[j] - pc_x[j] * tey_z_zz_1[j];

                    tez_xz_zz_0[j] = pa_x[j] * tez_z_zz_0[j] - pc_x[j] * tez_z_zz_1[j];

                    tex_yy_xx_0[j] = pa_y[j] * tex_y_xx_0[j] - pc_y[j] * tex_y_xx_1[j] + 0.5 * fl1_fx * tex_0_xx_0[j] - 0.5 * fl1_fx * tex_0_xx_1[j];

                    tey_yy_xx_0[j] = pa_y[j] * tey_y_xx_0[j] - pc_y[j] * tey_y_xx_1[j] + 0.5 * fl1_fx * tey_0_xx_0[j] - 0.5 * fl1_fx * tey_0_xx_1[j] + ta_y_xx_1[j];

                    tez_yy_xx_0[j] = pa_y[j] * tez_y_xx_0[j] - pc_y[j] * tez_y_xx_1[j] + 0.5 * fl1_fx * tez_0_xx_0[j] - 0.5 * fl1_fx * tez_0_xx_1[j];

                    tex_yy_xy_0[j] = pa_y[j] * tex_y_xy_0[j] - pc_y[j] * tex_y_xy_1[j] + 0.5 * fl1_fx * tex_0_xy_0[j] - 0.5 * fl1_fx * tex_0_xy_1[j] + 0.5 * fl1_fx * tex_y_x_0[j] - 0.5 * fl1_fx * tex_y_x_1[j];

                    tey_yy_xy_0[j] = pa_y[j] * tey_y_xy_0[j] - pc_y[j] * tey_y_xy_1[j] + 0.5 * fl1_fx * tey_0_xy_0[j] - 0.5 * fl1_fx * tey_0_xy_1[j] + 0.5 * fl1_fx * tey_y_x_0[j] - 0.5 * fl1_fx * tey_y_x_1[j] + ta_y_xy_1[j];

                    tez_yy_xy_0[j] = pa_y[j] * tez_y_xy_0[j] - pc_y[j] * tez_y_xy_1[j] + 0.5 * fl1_fx * tez_0_xy_0[j] - 0.5 * fl1_fx * tez_0_xy_1[j] + 0.5 * fl1_fx * tez_y_x_0[j] - 0.5 * fl1_fx * tez_y_x_1[j];

                    tex_yy_xz_0[j] = pa_y[j] * tex_y_xz_0[j] - pc_y[j] * tex_y_xz_1[j] + 0.5 * fl1_fx * tex_0_xz_0[j] - 0.5 * fl1_fx * tex_0_xz_1[j];

                    tey_yy_xz_0[j] = pa_y[j] * tey_y_xz_0[j] - pc_y[j] * tey_y_xz_1[j] + 0.5 * fl1_fx * tey_0_xz_0[j] - 0.5 * fl1_fx * tey_0_xz_1[j] + ta_y_xz_1[j];

                    tez_yy_xz_0[j] = pa_y[j] * tez_y_xz_0[j] - pc_y[j] * tez_y_xz_1[j] + 0.5 * fl1_fx * tez_0_xz_0[j] - 0.5 * fl1_fx * tez_0_xz_1[j];

                    tex_yy_yy_0[j] = pa_y[j] * tex_y_yy_0[j] - pc_y[j] * tex_y_yy_1[j] + 0.5 * fl1_fx * tex_0_yy_0[j] - 0.5 * fl1_fx * tex_0_yy_1[j] + fl1_fx * tex_y_y_0[j] - fl1_fx * tex_y_y_1[j];

                    tey_yy_yy_0[j] = pa_y[j] * tey_y_yy_0[j] - pc_y[j] * tey_y_yy_1[j] + 0.5 * fl1_fx * tey_0_yy_0[j] - 0.5 * fl1_fx * tey_0_yy_1[j] + fl1_fx * tey_y_y_0[j] - fl1_fx * tey_y_y_1[j] + ta_y_yy_1[j];

                    tez_yy_yy_0[j] = pa_y[j] * tez_y_yy_0[j] - pc_y[j] * tez_y_yy_1[j] + 0.5 * fl1_fx * tez_0_yy_0[j] - 0.5 * fl1_fx * tez_0_yy_1[j] + fl1_fx * tez_y_y_0[j] - fl1_fx * tez_y_y_1[j];

                    tex_yy_yz_0[j] = pa_y[j] * tex_y_yz_0[j] - pc_y[j] * tex_y_yz_1[j] + 0.5 * fl1_fx * tex_0_yz_0[j] - 0.5 * fl1_fx * tex_0_yz_1[j] + 0.5 * fl1_fx * tex_y_z_0[j] - 0.5 * fl1_fx * tex_y_z_1[j];

                    tey_yy_yz_0[j] = pa_y[j] * tey_y_yz_0[j] - pc_y[j] * tey_y_yz_1[j] + 0.5 * fl1_fx * tey_0_yz_0[j] - 0.5 * fl1_fx * tey_0_yz_1[j] + 0.5 * fl1_fx * tey_y_z_0[j] - 0.5 * fl1_fx * tey_y_z_1[j] + ta_y_yz_1[j];

                    tez_yy_yz_0[j] = pa_y[j] * tez_y_yz_0[j] - pc_y[j] * tez_y_yz_1[j] + 0.5 * fl1_fx * tez_0_yz_0[j] - 0.5 * fl1_fx * tez_0_yz_1[j] + 0.5 * fl1_fx * tez_y_z_0[j] - 0.5 * fl1_fx * tez_y_z_1[j];

                    tex_yy_zz_0[j] = pa_y[j] * tex_y_zz_0[j] - pc_y[j] * tex_y_zz_1[j] + 0.5 * fl1_fx * tex_0_zz_0[j] - 0.5 * fl1_fx * tex_0_zz_1[j];

                    tey_yy_zz_0[j] = pa_y[j] * tey_y_zz_0[j] - pc_y[j] * tey_y_zz_1[j] + 0.5 * fl1_fx * tey_0_zz_0[j] - 0.5 * fl1_fx * tey_0_zz_1[j] + ta_y_zz_1[j];

                    tez_yy_zz_0[j] = pa_y[j] * tez_y_zz_0[j] - pc_y[j] * tez_y_zz_1[j] + 0.5 * fl1_fx * tez_0_zz_0[j] - 0.5 * fl1_fx * tez_0_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDD_72_108(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (72,108)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 15); 

                auto tey_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 16); 

                auto tey_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 17); 

                auto tey_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 15); 

                auto tey_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 16); 

                auto tey_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 17); 

                auto tey_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 17); 

                auto tex_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx); 

                auto tey_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx); 

                auto tez_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx); 

                auto tex_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 1); 

                auto tey_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 1); 

                auto tez_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 1); 

                auto tex_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 2); 

                auto tey_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 2); 

                auto tez_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 2); 

                auto tex_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 3); 

                auto tey_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 3); 

                auto tez_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 3); 

                auto tex_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 4); 

                auto tey_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 4); 

                auto tez_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 4); 

                auto tex_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 5); 

                auto tey_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 5); 

                auto tez_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 5); 

                auto tex_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx); 

                auto tey_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx); 

                auto tez_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx); 

                auto tex_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 1); 

                auto tey_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 1); 

                auto tez_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 1); 

                auto tex_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 2); 

                auto tey_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 2); 

                auto tez_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 2); 

                auto tex_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 3); 

                auto tey_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 3); 

                auto tez_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 3); 

                auto tex_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 4); 

                auto tey_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 4); 

                auto tez_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 4); 

                auto tex_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 5); 

                auto tey_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 5); 

                auto tez_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 5); 

                auto tex_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 6); 

                auto tey_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 6); 

                auto tez_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 6); 

                auto tex_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 7); 

                auto tey_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 7); 

                auto tez_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 7); 

                auto tex_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 8); 

                auto tey_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 8); 

                auto tez_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 8); 

                auto tex_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 6); 

                auto tey_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 6); 

                auto tez_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 6); 

                auto tex_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 7); 

                auto tey_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 7); 

                auto tez_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 7); 

                auto tex_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 8); 

                auto tey_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 8); 

                auto tez_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 8); 

                auto ta_z_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 12); 

                auto ta_z_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 13); 

                auto ta_z_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 14); 

                auto ta_z_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 15); 

                auto ta_z_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 16); 

                auto ta_z_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 17); 

                // set up pointers to integrals

                auto tex_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 24); 

                auto tey_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 25); 

                auto tey_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 26); 

                auto tey_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 26); 

                auto tex_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 27); 

                auto tey_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 28); 

                auto tey_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 29); 

                auto tey_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 29); 

                auto tex_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 30); 

                auto tey_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 31); 

                auto tey_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 32); 

                auto tey_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 33); 

                auto tey_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 34); 

                auto tey_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 35); 

                auto tey_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 35); 

                // Batch of Integrals (72,108)

                #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_z_xx_1, ta_z_xy_1, ta_z_xz_1, ta_z_yy_1, \
                                         ta_z_yz_1, ta_z_zz_1, tex_0_xx_0, tex_0_xx_1, tex_0_xy_0, tex_0_xy_1, tex_0_xz_0, \
                                         tex_0_xz_1, tex_0_yy_0, tex_0_yy_1, tex_0_yz_0, tex_0_yz_1, tex_0_zz_0, tex_0_zz_1, \
                                         tex_yz_xx_0, tex_yz_xy_0, tex_yz_xz_0, tex_yz_yy_0, tex_yz_yz_0, tex_yz_zz_0, \
                                         tex_z_x_0, tex_z_x_1, tex_z_xx_0, tex_z_xx_1, tex_z_xy_0, tex_z_xy_1, tex_z_xz_0, \
                                         tex_z_xz_1, tex_z_y_0, tex_z_y_1, tex_z_yy_0, tex_z_yy_1, tex_z_yz_0, tex_z_yz_1, \
                                         tex_z_z_0, tex_z_z_1, tex_z_zz_0, tex_z_zz_1, tex_zz_xx_0, tex_zz_xy_0, \
                                         tex_zz_xz_0, tex_zz_yy_0, tex_zz_yz_0, tex_zz_zz_0, tey_0_xx_0, tey_0_xx_1, \
                                         tey_0_xy_0, tey_0_xy_1, tey_0_xz_0, tey_0_xz_1, tey_0_yy_0, tey_0_yy_1, tey_0_yz_0, \
                                         tey_0_yz_1, tey_0_zz_0, tey_0_zz_1, tey_yz_xx_0, tey_yz_xy_0, tey_yz_xz_0, \
                                         tey_yz_yy_0, tey_yz_yz_0, tey_yz_zz_0, tey_z_x_0, tey_z_x_1, tey_z_xx_0, tey_z_xx_1, \
                                         tey_z_xy_0, tey_z_xy_1, tey_z_xz_0, tey_z_xz_1, tey_z_y_0, tey_z_y_1, tey_z_yy_0, \
                                         tey_z_yy_1, tey_z_yz_0, tey_z_yz_1, tey_z_z_0, tey_z_z_1, tey_z_zz_0, tey_z_zz_1, \
                                         tey_zz_xx_0, tey_zz_xy_0, tey_zz_xz_0, tey_zz_yy_0, tey_zz_yz_0, tey_zz_zz_0, \
                                         tez_0_xx_0, tez_0_xx_1, tez_0_xy_0, tez_0_xy_1, tez_0_xz_0, tez_0_xz_1, tez_0_yy_0, \
                                         tez_0_yy_1, tez_0_yz_0, tez_0_yz_1, tez_0_zz_0, tez_0_zz_1, tez_yz_xx_0, \
                                         tez_yz_xy_0, tez_yz_xz_0, tez_yz_yy_0, tez_yz_yz_0, tez_yz_zz_0, tez_z_x_0, \
                                         tez_z_x_1, tez_z_xx_0, tez_z_xx_1, tez_z_xy_0, tez_z_xy_1, tez_z_xz_0, tez_z_xz_1, \
                                         tez_z_y_0, tez_z_y_1, tez_z_yy_0, tez_z_yy_1, tez_z_yz_0, tez_z_yz_1, tez_z_z_0, \
                                         tez_z_z_1, tez_z_zz_0, tez_z_zz_1, tez_zz_xx_0, tez_zz_xy_0, tez_zz_xz_0, \
                                         tez_zz_yy_0, tez_zz_yz_0, tez_zz_zz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yz_xx_0[j] = pa_y[j] * tex_z_xx_0[j] - pc_y[j] * tex_z_xx_1[j];

                    tey_yz_xx_0[j] = pa_y[j] * tey_z_xx_0[j] - pc_y[j] * tey_z_xx_1[j] + ta_z_xx_1[j];

                    tez_yz_xx_0[j] = pa_y[j] * tez_z_xx_0[j] - pc_y[j] * tez_z_xx_1[j];

                    tex_yz_xy_0[j] = pa_y[j] * tex_z_xy_0[j] - pc_y[j] * tex_z_xy_1[j] + 0.5 * fl1_fx * tex_z_x_0[j] - 0.5 * fl1_fx * tex_z_x_1[j];

                    tey_yz_xy_0[j] = pa_y[j] * tey_z_xy_0[j] - pc_y[j] * tey_z_xy_1[j] + 0.5 * fl1_fx * tey_z_x_0[j] - 0.5 * fl1_fx * tey_z_x_1[j] + ta_z_xy_1[j];

                    tez_yz_xy_0[j] = pa_y[j] * tez_z_xy_0[j] - pc_y[j] * tez_z_xy_1[j] + 0.5 * fl1_fx * tez_z_x_0[j] - 0.5 * fl1_fx * tez_z_x_1[j];

                    tex_yz_xz_0[j] = pa_y[j] * tex_z_xz_0[j] - pc_y[j] * tex_z_xz_1[j];

                    tey_yz_xz_0[j] = pa_y[j] * tey_z_xz_0[j] - pc_y[j] * tey_z_xz_1[j] + ta_z_xz_1[j];

                    tez_yz_xz_0[j] = pa_y[j] * tez_z_xz_0[j] - pc_y[j] * tez_z_xz_1[j];

                    tex_yz_yy_0[j] = pa_y[j] * tex_z_yy_0[j] - pc_y[j] * tex_z_yy_1[j] + fl1_fx * tex_z_y_0[j] - fl1_fx * tex_z_y_1[j];

                    tey_yz_yy_0[j] = pa_y[j] * tey_z_yy_0[j] - pc_y[j] * tey_z_yy_1[j] + fl1_fx * tey_z_y_0[j] - fl1_fx * tey_z_y_1[j] + ta_z_yy_1[j];

                    tez_yz_yy_0[j] = pa_y[j] * tez_z_yy_0[j] - pc_y[j] * tez_z_yy_1[j] + fl1_fx * tez_z_y_0[j] - fl1_fx * tez_z_y_1[j];

                    tex_yz_yz_0[j] = pa_y[j] * tex_z_yz_0[j] - pc_y[j] * tex_z_yz_1[j] + 0.5 * fl1_fx * tex_z_z_0[j] - 0.5 * fl1_fx * tex_z_z_1[j];

                    tey_yz_yz_0[j] = pa_y[j] * tey_z_yz_0[j] - pc_y[j] * tey_z_yz_1[j] + 0.5 * fl1_fx * tey_z_z_0[j] - 0.5 * fl1_fx * tey_z_z_1[j] + ta_z_yz_1[j];

                    tez_yz_yz_0[j] = pa_y[j] * tez_z_yz_0[j] - pc_y[j] * tez_z_yz_1[j] + 0.5 * fl1_fx * tez_z_z_0[j] - 0.5 * fl1_fx * tez_z_z_1[j];

                    tex_yz_zz_0[j] = pa_y[j] * tex_z_zz_0[j] - pc_y[j] * tex_z_zz_1[j];

                    tey_yz_zz_0[j] = pa_y[j] * tey_z_zz_0[j] - pc_y[j] * tey_z_zz_1[j] + ta_z_zz_1[j];

                    tez_yz_zz_0[j] = pa_y[j] * tez_z_zz_0[j] - pc_y[j] * tez_z_zz_1[j];

                    tex_zz_xx_0[j] = pa_z[j] * tex_z_xx_0[j] - pc_z[j] * tex_z_xx_1[j] + 0.5 * fl1_fx * tex_0_xx_0[j] - 0.5 * fl1_fx * tex_0_xx_1[j];

                    tey_zz_xx_0[j] = pa_z[j] * tey_z_xx_0[j] - pc_z[j] * tey_z_xx_1[j] + 0.5 * fl1_fx * tey_0_xx_0[j] - 0.5 * fl1_fx * tey_0_xx_1[j];

                    tez_zz_xx_0[j] = pa_z[j] * tez_z_xx_0[j] - pc_z[j] * tez_z_xx_1[j] + 0.5 * fl1_fx * tez_0_xx_0[j] - 0.5 * fl1_fx * tez_0_xx_1[j] + ta_z_xx_1[j];

                    tex_zz_xy_0[j] = pa_z[j] * tex_z_xy_0[j] - pc_z[j] * tex_z_xy_1[j] + 0.5 * fl1_fx * tex_0_xy_0[j] - 0.5 * fl1_fx * tex_0_xy_1[j];

                    tey_zz_xy_0[j] = pa_z[j] * tey_z_xy_0[j] - pc_z[j] * tey_z_xy_1[j] + 0.5 * fl1_fx * tey_0_xy_0[j] - 0.5 * fl1_fx * tey_0_xy_1[j];

                    tez_zz_xy_0[j] = pa_z[j] * tez_z_xy_0[j] - pc_z[j] * tez_z_xy_1[j] + 0.5 * fl1_fx * tez_0_xy_0[j] - 0.5 * fl1_fx * tez_0_xy_1[j] + ta_z_xy_1[j];

                    tex_zz_xz_0[j] = pa_z[j] * tex_z_xz_0[j] - pc_z[j] * tex_z_xz_1[j] + 0.5 * fl1_fx * tex_0_xz_0[j] - 0.5 * fl1_fx * tex_0_xz_1[j] + 0.5 * fl1_fx * tex_z_x_0[j] - 0.5 * fl1_fx * tex_z_x_1[j];

                    tey_zz_xz_0[j] = pa_z[j] * tey_z_xz_0[j] - pc_z[j] * tey_z_xz_1[j] + 0.5 * fl1_fx * tey_0_xz_0[j] - 0.5 * fl1_fx * tey_0_xz_1[j] + 0.5 * fl1_fx * tey_z_x_0[j] - 0.5 * fl1_fx * tey_z_x_1[j];

                    tez_zz_xz_0[j] = pa_z[j] * tez_z_xz_0[j] - pc_z[j] * tez_z_xz_1[j] + 0.5 * fl1_fx * tez_0_xz_0[j] - 0.5 * fl1_fx * tez_0_xz_1[j] + 0.5 * fl1_fx * tez_z_x_0[j] - 0.5 * fl1_fx * tez_z_x_1[j] + ta_z_xz_1[j];

                    tex_zz_yy_0[j] = pa_z[j] * tex_z_yy_0[j] - pc_z[j] * tex_z_yy_1[j] + 0.5 * fl1_fx * tex_0_yy_0[j] - 0.5 * fl1_fx * tex_0_yy_1[j];

                    tey_zz_yy_0[j] = pa_z[j] * tey_z_yy_0[j] - pc_z[j] * tey_z_yy_1[j] + 0.5 * fl1_fx * tey_0_yy_0[j] - 0.5 * fl1_fx * tey_0_yy_1[j];

                    tez_zz_yy_0[j] = pa_z[j] * tez_z_yy_0[j] - pc_z[j] * tez_z_yy_1[j] + 0.5 * fl1_fx * tez_0_yy_0[j] - 0.5 * fl1_fx * tez_0_yy_1[j] + ta_z_yy_1[j];

                    tex_zz_yz_0[j] = pa_z[j] * tex_z_yz_0[j] - pc_z[j] * tex_z_yz_1[j] + 0.5 * fl1_fx * tex_0_yz_0[j] - 0.5 * fl1_fx * tex_0_yz_1[j] + 0.5 * fl1_fx * tex_z_y_0[j] - 0.5 * fl1_fx * tex_z_y_1[j];

                    tey_zz_yz_0[j] = pa_z[j] * tey_z_yz_0[j] - pc_z[j] * tey_z_yz_1[j] + 0.5 * fl1_fx * tey_0_yz_0[j] - 0.5 * fl1_fx * tey_0_yz_1[j] + 0.5 * fl1_fx * tey_z_y_0[j] - 0.5 * fl1_fx * tey_z_y_1[j];

                    tez_zz_yz_0[j] = pa_z[j] * tez_z_yz_0[j] - pc_z[j] * tez_z_yz_1[j] + 0.5 * fl1_fx * tez_0_yz_0[j] - 0.5 * fl1_fx * tez_0_yz_1[j] + 0.5 * fl1_fx * tez_z_y_0[j] - 0.5 * fl1_fx * tez_z_y_1[j] + ta_z_yz_1[j];

                    tex_zz_zz_0[j] = pa_z[j] * tex_z_zz_0[j] - pc_z[j] * tex_z_zz_1[j] + 0.5 * fl1_fx * tex_0_zz_0[j] - 0.5 * fl1_fx * tex_0_zz_1[j] + fl1_fx * tex_z_z_0[j] - fl1_fx * tex_z_z_1[j];

                    tey_zz_zz_0[j] = pa_z[j] * tey_z_zz_0[j] - pc_z[j] * tey_z_zz_1[j] + 0.5 * fl1_fx * tey_0_zz_0[j] - 0.5 * fl1_fx * tey_0_zz_1[j] + fl1_fx * tey_z_z_0[j] - fl1_fx * tey_z_z_1[j];

                    tez_zz_zz_0[j] = pa_z[j] * tez_z_zz_0[j] - pc_z[j] * tez_z_zz_1[j] + 0.5 * fl1_fx * tez_0_zz_0[j] - 0.5 * fl1_fx * tez_0_zz_1[j] + fl1_fx * tez_z_z_0[j] - fl1_fx * tez_z_z_1[j] + ta_z_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDF(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pcDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        efieldrecfunc::compElectricFieldForDF_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   pcDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        efieldrecfunc::compElectricFieldForDF_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    pcDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        efieldrecfunc::compElectricFieldForDF_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     pcDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        efieldrecfunc::compElectricFieldForDF_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compElectricFieldForDF_0_45(      CMemBlock2D<double>& primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
    {
        // Batch of Integrals (0,45)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tex_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx); 

                auto tey_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx); 

                auto tez_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx); 

                auto tex_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 1); 

                auto tey_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 1); 

                auto tez_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 1); 

                auto tex_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 2); 

                auto tey_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 2); 

                auto tez_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 2); 

                auto tex_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 3); 

                auto tey_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 3); 

                auto tez_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 3); 

                auto tex_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 4); 

                auto tey_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 4); 

                auto tez_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 4); 

                auto tex_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 5); 

                auto tey_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 5); 

                auto tez_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 5); 

                auto tex_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 6); 

                auto tey_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 6); 

                auto tez_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 6); 

                auto tex_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 7); 

                auto tey_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 7); 

                auto tez_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 7); 

                auto tex_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 8); 

                auto tey_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 8); 

                auto tez_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 8); 

                auto tex_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 9); 

                auto tey_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 9); 

                auto tez_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 9); 

                auto tex_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 10); 

                auto tey_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 11); 

                auto tey_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 12); 

                auto tey_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 13); 

                auto tey_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 14); 

                auto tey_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 14); 

                auto tex_x_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx); 

                auto tey_x_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx); 

                auto tez_x_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx); 

                auto tex_x_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 1); 

                auto tey_x_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 1); 

                auto tez_x_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 1); 

                auto tex_x_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 2); 

                auto tey_x_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 2); 

                auto tez_x_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 2); 

                auto tex_x_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 3); 

                auto tey_x_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 3); 

                auto tez_x_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 3); 

                auto tex_x_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 4); 

                auto tey_x_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 4); 

                auto tez_x_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 4); 

                auto tex_x_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 5); 

                auto tey_x_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 5); 

                auto tez_x_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 5); 

                auto tex_x_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 6); 

                auto tey_x_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 6); 

                auto tez_x_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 6); 

                auto tex_x_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 7); 

                auto tey_x_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 7); 

                auto tez_x_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 7); 

                auto tex_x_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 8); 

                auto tey_x_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 8); 

                auto tez_x_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 8); 

                auto tex_x_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 9); 

                auto tey_x_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 9); 

                auto tez_x_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 9); 

                auto tex_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 10); 

                auto tey_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 11); 

                auto tey_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 12); 

                auto tey_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 13); 

                auto tey_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 14); 

                auto tey_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 14); 

                auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx); 

                auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx); 

                auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx); 

                auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1); 

                auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1); 

                auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1); 

                auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2); 

                auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2); 

                auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2); 

                auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3); 

                auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3); 

                auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3); 

                auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4); 

                auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4); 

                auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4); 

                auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5); 

                auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5); 

                auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5); 

                auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6); 

                auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6); 

                auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6); 

                auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7); 

                auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7); 

                auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7); 

                auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8); 

                auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8); 

                auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8); 

                auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9); 

                auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9); 

                auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9); 

                auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx); 

                auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx); 

                auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx); 

                auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1); 

                auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1); 

                auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1); 

                auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2); 

                auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2); 

                auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2); 

                auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3); 

                auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3); 

                auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3); 

                auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4); 

                auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4); 

                auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4); 

                auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5); 

                auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5); 

                auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5); 

                auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6); 

                auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6); 

                auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6); 

                auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7); 

                auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7); 

                auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7); 

                auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8); 

                auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8); 

                auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8); 

                auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9); 

                auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9); 

                auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9); 

                auto tex_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx); 

                auto tey_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx); 

                auto tez_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx); 

                auto tex_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 1); 

                auto tey_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 1); 

                auto tez_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 1); 

                auto tex_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 2); 

                auto tey_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 2); 

                auto tez_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 2); 

                auto tex_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 3); 

                auto tey_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 3); 

                auto tez_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 3); 

                auto tex_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 4); 

                auto tey_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 4); 

                auto tez_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 4); 

                auto tex_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 5); 

                auto tey_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 5); 

                auto tez_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 5); 

                auto tex_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 6); 

                auto tey_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 7); 

                auto tey_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 8); 

                auto tey_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 9); 

                auto tey_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 10); 

                auto tey_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx); 

                auto tey_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx); 

                auto tez_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx); 

                auto tex_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 1); 

                auto tey_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 1); 

                auto tez_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 1); 

                auto tex_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 2); 

                auto tey_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 2); 

                auto tez_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 2); 

                auto tex_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 3); 

                auto tey_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 3); 

                auto tez_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 3); 

                auto tex_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 4); 

                auto tey_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 4); 

                auto tez_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 4); 

                auto tex_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 5); 

                auto tey_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 5); 

                auto tez_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 5); 

                auto tex_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 6); 

                auto tey_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 7); 

                auto tey_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 8); 

                auto tey_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 9); 

                auto tey_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 10); 

                auto tey_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 10); 

                auto ta_x_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx); 

                auto ta_x_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 1); 

                auto ta_x_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 2); 

                auto ta_x_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 3); 

                auto ta_x_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 4); 

                auto ta_x_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 5); 

                auto ta_x_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 6); 

                auto ta_x_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 7); 

                auto ta_x_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 8); 

                auto ta_x_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 9); 

                auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10); 

                auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11); 

                auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12); 

                auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13); 

                auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14); 

                // set up pointers to integrals

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

                // Batch of Integrals (0,45)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_x_xxx_1, ta_x_xxy_1, ta_x_xxz_1, ta_x_xyy_1, ta_x_xyz_1, \
                                         ta_x_xzz_1, ta_x_yyy_1, ta_x_yyz_1, ta_x_yzz_1, ta_x_zzz_1, ta_y_xxx_1, ta_y_xxy_1, \
                                         ta_y_xxz_1, ta_y_xyy_1, ta_y_xyz_1, tex_0_xxx_0, tex_0_xxx_1, tex_0_xxy_0, \
                                         tex_0_xxy_1, tex_0_xxz_0, tex_0_xxz_1, tex_0_xyy_0, tex_0_xyy_1, tex_0_xyz_0, \
                                         tex_0_xyz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_yyy_0, tex_0_yyy_1, tex_0_yyz_0, \
                                         tex_0_yyz_1, tex_0_yzz_0, tex_0_yzz_1, tex_0_zzz_0, tex_0_zzz_1, tex_x_xx_0, \
                                         tex_x_xx_1, tex_x_xxx_0, tex_x_xxx_1, tex_x_xxy_0, tex_x_xxy_1, tex_x_xxz_0, \
                                         tex_x_xxz_1, tex_x_xy_0, tex_x_xy_1, tex_x_xyy_0, tex_x_xyy_1, tex_x_xyz_0, \
                                         tex_x_xyz_1, tex_x_xz_0, tex_x_xz_1, tex_x_xzz_0, tex_x_xzz_1, tex_x_yy_0, \
                                         tex_x_yy_1, tex_x_yyy_0, tex_x_yyy_1, tex_x_yyz_0, tex_x_yyz_1, tex_x_yz_0, \
                                         tex_x_yz_1, tex_x_yzz_0, tex_x_yzz_1, tex_x_zz_0, tex_x_zz_1, tex_x_zzz_0, \
                                         tex_x_zzz_1, tex_xx_xxx_0, tex_xx_xxy_0, tex_xx_xxz_0, tex_xx_xyy_0, tex_xx_xyz_0, \
                                         tex_xx_xzz_0, tex_xx_yyy_0, tex_xx_yyz_0, tex_xx_yzz_0, tex_xx_zzz_0, tex_xy_xxx_0, \
                                         tex_xy_xxy_0, tex_xy_xxz_0, tex_xy_xyy_0, tex_xy_xyz_0, tex_y_xx_0, tex_y_xx_1, \
                                         tex_y_xxx_0, tex_y_xxx_1, tex_y_xxy_0, tex_y_xxy_1, tex_y_xxz_0, tex_y_xxz_1, \
                                         tex_y_xy_0, tex_y_xy_1, tex_y_xyy_0, tex_y_xyy_1, tex_y_xyz_0, tex_y_xyz_1, \
                                         tex_y_xz_0, tex_y_xz_1, tex_y_yy_0, tex_y_yy_1, tex_y_yz_0, tex_y_yz_1, \
                                         tey_0_xxx_0, tey_0_xxx_1, tey_0_xxy_0, tey_0_xxy_1, tey_0_xxz_0, tey_0_xxz_1, \
                                         tey_0_xyy_0, tey_0_xyy_1, tey_0_xyz_0, tey_0_xyz_1, tey_0_xzz_0, tey_0_xzz_1, \
                                         tey_0_yyy_0, tey_0_yyy_1, tey_0_yyz_0, tey_0_yyz_1, tey_0_yzz_0, tey_0_yzz_1, \
                                         tey_0_zzz_0, tey_0_zzz_1, tey_x_xx_0, tey_x_xx_1, tey_x_xxx_0, tey_x_xxx_1, \
                                         tey_x_xxy_0, tey_x_xxy_1, tey_x_xxz_0, tey_x_xxz_1, tey_x_xy_0, tey_x_xy_1, \
                                         tey_x_xyy_0, tey_x_xyy_1, tey_x_xyz_0, tey_x_xyz_1, tey_x_xz_0, tey_x_xz_1, \
                                         tey_x_xzz_0, tey_x_xzz_1, tey_x_yy_0, tey_x_yy_1, tey_x_yyy_0, tey_x_yyy_1, \
                                         tey_x_yyz_0, tey_x_yyz_1, tey_x_yz_0, tey_x_yz_1, tey_x_yzz_0, tey_x_yzz_1, \
                                         tey_x_zz_0, tey_x_zz_1, tey_x_zzz_0, tey_x_zzz_1, tey_xx_xxx_0, tey_xx_xxy_0, \
                                         tey_xx_xxz_0, tey_xx_xyy_0, tey_xx_xyz_0, tey_xx_xzz_0, tey_xx_yyy_0, tey_xx_yyz_0, \
                                         tey_xx_yzz_0, tey_xx_zzz_0, tey_xy_xxx_0, tey_xy_xxy_0, tey_xy_xxz_0, tey_xy_xyy_0, \
                                         tey_xy_xyz_0, tey_y_xx_0, tey_y_xx_1, tey_y_xxx_0, tey_y_xxx_1, tey_y_xxy_0, \
                                         tey_y_xxy_1, tey_y_xxz_0, tey_y_xxz_1, tey_y_xy_0, tey_y_xy_1, tey_y_xyy_0, \
                                         tey_y_xyy_1, tey_y_xyz_0, tey_y_xyz_1, tey_y_xz_0, tey_y_xz_1, tey_y_yy_0, \
                                         tey_y_yy_1, tey_y_yz_0, tey_y_yz_1, tez_0_xxx_0, tez_0_xxx_1, tez_0_xxy_0, \
                                         tez_0_xxy_1, tez_0_xxz_0, tez_0_xxz_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyz_0, \
                                         tez_0_xyz_1, tez_0_xzz_0, tez_0_xzz_1, tez_0_yyy_0, tez_0_yyy_1, tez_0_yyz_0, \
                                         tez_0_yyz_1, tez_0_yzz_0, tez_0_yzz_1, tez_0_zzz_0, tez_0_zzz_1, tez_x_xx_0, \
                                         tez_x_xx_1, tez_x_xxx_0, tez_x_xxx_1, tez_x_xxy_0, tez_x_xxy_1, tez_x_xxz_0, \
                                         tez_x_xxz_1, tez_x_xy_0, tez_x_xy_1, tez_x_xyy_0, tez_x_xyy_1, tez_x_xyz_0, \
                                         tez_x_xyz_1, tez_x_xz_0, tez_x_xz_1, tez_x_xzz_0, tez_x_xzz_1, tez_x_yy_0, \
                                         tez_x_yy_1, tez_x_yyy_0, tez_x_yyy_1, tez_x_yyz_0, tez_x_yyz_1, tez_x_yz_0, \
                                         tez_x_yz_1, tez_x_yzz_0, tez_x_yzz_1, tez_x_zz_0, tez_x_zz_1, tez_x_zzz_0, \
                                         tez_x_zzz_1, tez_xx_xxx_0, tez_xx_xxy_0, tez_xx_xxz_0, tez_xx_xyy_0, tez_xx_xyz_0, \
                                         tez_xx_xzz_0, tez_xx_yyy_0, tez_xx_yyz_0, tez_xx_yzz_0, tez_xx_zzz_0, tez_xy_xxx_0, \
                                         tez_xy_xxy_0, tez_xy_xxz_0, tez_xy_xyy_0, tez_xy_xyz_0, tez_y_xx_0, tez_y_xx_1, \
                                         tez_y_xxx_0, tez_y_xxx_1, tez_y_xxy_0, tez_y_xxy_1, tez_y_xxz_0, tez_y_xxz_1, \
                                         tez_y_xy_0, tez_y_xy_1, tez_y_xyy_0, tez_y_xyy_1, tez_y_xyz_0, tez_y_xyz_1, \
                                         tez_y_xz_0, tez_y_xz_1, tez_y_yy_0, tez_y_yy_1, tez_y_yz_0, tez_y_yz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xx_xxx_0[j] = pa_x[j] * tex_x_xxx_0[j] - pc_x[j] * tex_x_xxx_1[j] + 0.5 * fl1_fx * tex_0_xxx_0[j] - 0.5 * fl1_fx * tex_0_xxx_1[j] + 1.5 * fl1_fx * tex_x_xx_0[j] - 1.5 * fl1_fx * tex_x_xx_1[j] + ta_x_xxx_1[j];

                    tey_xx_xxx_0[j] = pa_x[j] * tey_x_xxx_0[j] - pc_x[j] * tey_x_xxx_1[j] + 0.5 * fl1_fx * tey_0_xxx_0[j] - 0.5 * fl1_fx * tey_0_xxx_1[j] + 1.5 * fl1_fx * tey_x_xx_0[j] - 1.5 * fl1_fx * tey_x_xx_1[j];

                    tez_xx_xxx_0[j] = pa_x[j] * tez_x_xxx_0[j] - pc_x[j] * tez_x_xxx_1[j] + 0.5 * fl1_fx * tez_0_xxx_0[j] - 0.5 * fl1_fx * tez_0_xxx_1[j] + 1.5 * fl1_fx * tez_x_xx_0[j] - 1.5 * fl1_fx * tez_x_xx_1[j];

                    tex_xx_xxy_0[j] = pa_x[j] * tex_x_xxy_0[j] - pc_x[j] * tex_x_xxy_1[j] + 0.5 * fl1_fx * tex_0_xxy_0[j] - 0.5 * fl1_fx * tex_0_xxy_1[j] + fl1_fx * tex_x_xy_0[j] - fl1_fx * tex_x_xy_1[j] + ta_x_xxy_1[j];

                    tey_xx_xxy_0[j] = pa_x[j] * tey_x_xxy_0[j] - pc_x[j] * tey_x_xxy_1[j] + 0.5 * fl1_fx * tey_0_xxy_0[j] - 0.5 * fl1_fx * tey_0_xxy_1[j] + fl1_fx * tey_x_xy_0[j] - fl1_fx * tey_x_xy_1[j];

                    tez_xx_xxy_0[j] = pa_x[j] * tez_x_xxy_0[j] - pc_x[j] * tez_x_xxy_1[j] + 0.5 * fl1_fx * tez_0_xxy_0[j] - 0.5 * fl1_fx * tez_0_xxy_1[j] + fl1_fx * tez_x_xy_0[j] - fl1_fx * tez_x_xy_1[j];

                    tex_xx_xxz_0[j] = pa_x[j] * tex_x_xxz_0[j] - pc_x[j] * tex_x_xxz_1[j] + 0.5 * fl1_fx * tex_0_xxz_0[j] - 0.5 * fl1_fx * tex_0_xxz_1[j] + fl1_fx * tex_x_xz_0[j] - fl1_fx * tex_x_xz_1[j] + ta_x_xxz_1[j];

                    tey_xx_xxz_0[j] = pa_x[j] * tey_x_xxz_0[j] - pc_x[j] * tey_x_xxz_1[j] + 0.5 * fl1_fx * tey_0_xxz_0[j] - 0.5 * fl1_fx * tey_0_xxz_1[j] + fl1_fx * tey_x_xz_0[j] - fl1_fx * tey_x_xz_1[j];

                    tez_xx_xxz_0[j] = pa_x[j] * tez_x_xxz_0[j] - pc_x[j] * tez_x_xxz_1[j] + 0.5 * fl1_fx * tez_0_xxz_0[j] - 0.5 * fl1_fx * tez_0_xxz_1[j] + fl1_fx * tez_x_xz_0[j] - fl1_fx * tez_x_xz_1[j];

                    tex_xx_xyy_0[j] = pa_x[j] * tex_x_xyy_0[j] - pc_x[j] * tex_x_xyy_1[j] + 0.5 * fl1_fx * tex_0_xyy_0[j] - 0.5 * fl1_fx * tex_0_xyy_1[j] + 0.5 * fl1_fx * tex_x_yy_0[j] - 0.5 * fl1_fx * tex_x_yy_1[j] + ta_x_xyy_1[j];

                    tey_xx_xyy_0[j] = pa_x[j] * tey_x_xyy_0[j] - pc_x[j] * tey_x_xyy_1[j] + 0.5 * fl1_fx * tey_0_xyy_0[j] - 0.5 * fl1_fx * tey_0_xyy_1[j] + 0.5 * fl1_fx * tey_x_yy_0[j] - 0.5 * fl1_fx * tey_x_yy_1[j];

                    tez_xx_xyy_0[j] = pa_x[j] * tez_x_xyy_0[j] - pc_x[j] * tez_x_xyy_1[j] + 0.5 * fl1_fx * tez_0_xyy_0[j] - 0.5 * fl1_fx * tez_0_xyy_1[j] + 0.5 * fl1_fx * tez_x_yy_0[j] - 0.5 * fl1_fx * tez_x_yy_1[j];

                    tex_xx_xyz_0[j] = pa_x[j] * tex_x_xyz_0[j] - pc_x[j] * tex_x_xyz_1[j] + 0.5 * fl1_fx * tex_0_xyz_0[j] - 0.5 * fl1_fx * tex_0_xyz_1[j] + 0.5 * fl1_fx * tex_x_yz_0[j] - 0.5 * fl1_fx * tex_x_yz_1[j] + ta_x_xyz_1[j];

                    tey_xx_xyz_0[j] = pa_x[j] * tey_x_xyz_0[j] - pc_x[j] * tey_x_xyz_1[j] + 0.5 * fl1_fx * tey_0_xyz_0[j] - 0.5 * fl1_fx * tey_0_xyz_1[j] + 0.5 * fl1_fx * tey_x_yz_0[j] - 0.5 * fl1_fx * tey_x_yz_1[j];

                    tez_xx_xyz_0[j] = pa_x[j] * tez_x_xyz_0[j] - pc_x[j] * tez_x_xyz_1[j] + 0.5 * fl1_fx * tez_0_xyz_0[j] - 0.5 * fl1_fx * tez_0_xyz_1[j] + 0.5 * fl1_fx * tez_x_yz_0[j] - 0.5 * fl1_fx * tez_x_yz_1[j];

                    tex_xx_xzz_0[j] = pa_x[j] * tex_x_xzz_0[j] - pc_x[j] * tex_x_xzz_1[j] + 0.5 * fl1_fx * tex_0_xzz_0[j] - 0.5 * fl1_fx * tex_0_xzz_1[j] + 0.5 * fl1_fx * tex_x_zz_0[j] - 0.5 * fl1_fx * tex_x_zz_1[j] + ta_x_xzz_1[j];

                    tey_xx_xzz_0[j] = pa_x[j] * tey_x_xzz_0[j] - pc_x[j] * tey_x_xzz_1[j] + 0.5 * fl1_fx * tey_0_xzz_0[j] - 0.5 * fl1_fx * tey_0_xzz_1[j] + 0.5 * fl1_fx * tey_x_zz_0[j] - 0.5 * fl1_fx * tey_x_zz_1[j];

                    tez_xx_xzz_0[j] = pa_x[j] * tez_x_xzz_0[j] - pc_x[j] * tez_x_xzz_1[j] + 0.5 * fl1_fx * tez_0_xzz_0[j] - 0.5 * fl1_fx * tez_0_xzz_1[j] + 0.5 * fl1_fx * tez_x_zz_0[j] - 0.5 * fl1_fx * tez_x_zz_1[j];

                    tex_xx_yyy_0[j] = pa_x[j] * tex_x_yyy_0[j] - pc_x[j] * tex_x_yyy_1[j] + 0.5 * fl1_fx * tex_0_yyy_0[j] - 0.5 * fl1_fx * tex_0_yyy_1[j] + ta_x_yyy_1[j];

                    tey_xx_yyy_0[j] = pa_x[j] * tey_x_yyy_0[j] - pc_x[j] * tey_x_yyy_1[j] + 0.5 * fl1_fx * tey_0_yyy_0[j] - 0.5 * fl1_fx * tey_0_yyy_1[j];

                    tez_xx_yyy_0[j] = pa_x[j] * tez_x_yyy_0[j] - pc_x[j] * tez_x_yyy_1[j] + 0.5 * fl1_fx * tez_0_yyy_0[j] - 0.5 * fl1_fx * tez_0_yyy_1[j];

                    tex_xx_yyz_0[j] = pa_x[j] * tex_x_yyz_0[j] - pc_x[j] * tex_x_yyz_1[j] + 0.5 * fl1_fx * tex_0_yyz_0[j] - 0.5 * fl1_fx * tex_0_yyz_1[j] + ta_x_yyz_1[j];

                    tey_xx_yyz_0[j] = pa_x[j] * tey_x_yyz_0[j] - pc_x[j] * tey_x_yyz_1[j] + 0.5 * fl1_fx * tey_0_yyz_0[j] - 0.5 * fl1_fx * tey_0_yyz_1[j];

                    tez_xx_yyz_0[j] = pa_x[j] * tez_x_yyz_0[j] - pc_x[j] * tez_x_yyz_1[j] + 0.5 * fl1_fx * tez_0_yyz_0[j] - 0.5 * fl1_fx * tez_0_yyz_1[j];

                    tex_xx_yzz_0[j] = pa_x[j] * tex_x_yzz_0[j] - pc_x[j] * tex_x_yzz_1[j] + 0.5 * fl1_fx * tex_0_yzz_0[j] - 0.5 * fl1_fx * tex_0_yzz_1[j] + ta_x_yzz_1[j];

                    tey_xx_yzz_0[j] = pa_x[j] * tey_x_yzz_0[j] - pc_x[j] * tey_x_yzz_1[j] + 0.5 * fl1_fx * tey_0_yzz_0[j] - 0.5 * fl1_fx * tey_0_yzz_1[j];

                    tez_xx_yzz_0[j] = pa_x[j] * tez_x_yzz_0[j] - pc_x[j] * tez_x_yzz_1[j] + 0.5 * fl1_fx * tez_0_yzz_0[j] - 0.5 * fl1_fx * tez_0_yzz_1[j];

                    tex_xx_zzz_0[j] = pa_x[j] * tex_x_zzz_0[j] - pc_x[j] * tex_x_zzz_1[j] + 0.5 * fl1_fx * tex_0_zzz_0[j] - 0.5 * fl1_fx * tex_0_zzz_1[j] + ta_x_zzz_1[j];

                    tey_xx_zzz_0[j] = pa_x[j] * tey_x_zzz_0[j] - pc_x[j] * tey_x_zzz_1[j] + 0.5 * fl1_fx * tey_0_zzz_0[j] - 0.5 * fl1_fx * tey_0_zzz_1[j];

                    tez_xx_zzz_0[j] = pa_x[j] * tez_x_zzz_0[j] - pc_x[j] * tez_x_zzz_1[j] + 0.5 * fl1_fx * tez_0_zzz_0[j] - 0.5 * fl1_fx * tez_0_zzz_1[j];

                    tex_xy_xxx_0[j] = pa_x[j] * tex_y_xxx_0[j] - pc_x[j] * tex_y_xxx_1[j] + 1.5 * fl1_fx * tex_y_xx_0[j] - 1.5 * fl1_fx * tex_y_xx_1[j] + ta_y_xxx_1[j];

                    tey_xy_xxx_0[j] = pa_x[j] * tey_y_xxx_0[j] - pc_x[j] * tey_y_xxx_1[j] + 1.5 * fl1_fx * tey_y_xx_0[j] - 1.5 * fl1_fx * tey_y_xx_1[j];

                    tez_xy_xxx_0[j] = pa_x[j] * tez_y_xxx_0[j] - pc_x[j] * tez_y_xxx_1[j] + 1.5 * fl1_fx * tez_y_xx_0[j] - 1.5 * fl1_fx * tez_y_xx_1[j];

                    tex_xy_xxy_0[j] = pa_x[j] * tex_y_xxy_0[j] - pc_x[j] * tex_y_xxy_1[j] + fl1_fx * tex_y_xy_0[j] - fl1_fx * tex_y_xy_1[j] + ta_y_xxy_1[j];

                    tey_xy_xxy_0[j] = pa_x[j] * tey_y_xxy_0[j] - pc_x[j] * tey_y_xxy_1[j] + fl1_fx * tey_y_xy_0[j] - fl1_fx * tey_y_xy_1[j];

                    tez_xy_xxy_0[j] = pa_x[j] * tez_y_xxy_0[j] - pc_x[j] * tez_y_xxy_1[j] + fl1_fx * tez_y_xy_0[j] - fl1_fx * tez_y_xy_1[j];

                    tex_xy_xxz_0[j] = pa_x[j] * tex_y_xxz_0[j] - pc_x[j] * tex_y_xxz_1[j] + fl1_fx * tex_y_xz_0[j] - fl1_fx * tex_y_xz_1[j] + ta_y_xxz_1[j];

                    tey_xy_xxz_0[j] = pa_x[j] * tey_y_xxz_0[j] - pc_x[j] * tey_y_xxz_1[j] + fl1_fx * tey_y_xz_0[j] - fl1_fx * tey_y_xz_1[j];

                    tez_xy_xxz_0[j] = pa_x[j] * tez_y_xxz_0[j] - pc_x[j] * tez_y_xxz_1[j] + fl1_fx * tez_y_xz_0[j] - fl1_fx * tez_y_xz_1[j];

                    tex_xy_xyy_0[j] = pa_x[j] * tex_y_xyy_0[j] - pc_x[j] * tex_y_xyy_1[j] + 0.5 * fl1_fx * tex_y_yy_0[j] - 0.5 * fl1_fx * tex_y_yy_1[j] + ta_y_xyy_1[j];

                    tey_xy_xyy_0[j] = pa_x[j] * tey_y_xyy_0[j] - pc_x[j] * tey_y_xyy_1[j] + 0.5 * fl1_fx * tey_y_yy_0[j] - 0.5 * fl1_fx * tey_y_yy_1[j];

                    tez_xy_xyy_0[j] = pa_x[j] * tez_y_xyy_0[j] - pc_x[j] * tez_y_xyy_1[j] + 0.5 * fl1_fx * tez_y_yy_0[j] - 0.5 * fl1_fx * tez_y_yy_1[j];

                    tex_xy_xyz_0[j] = pa_x[j] * tex_y_xyz_0[j] - pc_x[j] * tex_y_xyz_1[j] + 0.5 * fl1_fx * tex_y_yz_0[j] - 0.5 * fl1_fx * tex_y_yz_1[j] + ta_y_xyz_1[j];

                    tey_xy_xyz_0[j] = pa_x[j] * tey_y_xyz_0[j] - pc_x[j] * tey_y_xyz_1[j] + 0.5 * fl1_fx * tey_y_yz_0[j] - 0.5 * fl1_fx * tey_y_yz_1[j];

                    tez_xy_xyz_0[j] = pa_x[j] * tez_y_xyz_0[j] - pc_x[j] * tez_y_xyz_1[j] + 0.5 * fl1_fx * tez_y_yz_0[j] - 0.5 * fl1_fx * tez_y_yz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDF_45_90(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (45,90)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tex_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 15); 

                auto tey_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 16); 

                auto tey_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 17); 

                auto tey_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 18); 

                auto tey_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 19); 

                auto tey_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 19); 

                auto tex_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 20); 

                auto tey_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 21); 

                auto tey_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 22); 

                auto tey_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 23); 

                auto tey_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 24); 

                auto tey_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 25); 

                auto tey_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 26); 

                auto tey_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 27); 

                auto tey_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 28); 

                auto tey_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 29); 

                auto tey_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 29); 

                auto tex_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 15); 

                auto tey_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 16); 

                auto tey_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 17); 

                auto tey_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 18); 

                auto tey_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 19); 

                auto tey_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 19); 

                auto tex_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 20); 

                auto tey_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 21); 

                auto tey_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 22); 

                auto tey_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 23); 

                auto tey_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 24); 

                auto tey_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 25); 

                auto tey_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 26); 

                auto tey_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 27); 

                auto tey_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 28); 

                auto tey_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 29); 

                auto tey_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 29); 

                auto tex_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 11); 

                auto tey_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 15); 

                auto tey_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 16); 

                auto tey_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 17); 

                auto tey_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 11); 

                auto tey_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 15); 

                auto tey_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 16); 

                auto tey_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 17); 

                auto tey_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 17); 

                auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15); 

                auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16); 

                auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17); 

                auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18); 

                auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19); 

                auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20); 

                auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21); 

                auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22); 

                auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23); 

                auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24); 

                auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25); 

                auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26); 

                auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27); 

                auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28); 

                auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29); 

                // set up pointers to integrals

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

                // Batch of Integrals (45,90)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_y_xzz_1, ta_y_yyy_1, ta_y_yyz_1, ta_y_yzz_1, ta_y_zzz_1, \
                                         ta_z_xxx_1, ta_z_xxy_1, ta_z_xxz_1, ta_z_xyy_1, ta_z_xyz_1, ta_z_xzz_1, ta_z_yyy_1, \
                                         ta_z_yyz_1, ta_z_yzz_1, ta_z_zzz_1, tex_xy_xzz_0, tex_xy_yyy_0, tex_xy_yyz_0, \
                                         tex_xy_yzz_0, tex_xy_zzz_0, tex_xz_xxx_0, tex_xz_xxy_0, tex_xz_xxz_0, tex_xz_xyy_0, \
                                         tex_xz_xyz_0, tex_xz_xzz_0, tex_xz_yyy_0, tex_xz_yyz_0, tex_xz_yzz_0, tex_xz_zzz_0, \
                                         tex_y_xzz_0, tex_y_xzz_1, tex_y_yyy_0, tex_y_yyy_1, tex_y_yyz_0, tex_y_yyz_1, \
                                         tex_y_yzz_0, tex_y_yzz_1, tex_y_zz_0, tex_y_zz_1, tex_y_zzz_0, tex_y_zzz_1, \
                                         tex_z_xx_0, tex_z_xx_1, tex_z_xxx_0, tex_z_xxx_1, tex_z_xxy_0, tex_z_xxy_1, \
                                         tex_z_xxz_0, tex_z_xxz_1, tex_z_xy_0, tex_z_xy_1, tex_z_xyy_0, tex_z_xyy_1, \
                                         tex_z_xyz_0, tex_z_xyz_1, tex_z_xz_0, tex_z_xz_1, tex_z_xzz_0, tex_z_xzz_1, \
                                         tex_z_yy_0, tex_z_yy_1, tex_z_yyy_0, tex_z_yyy_1, tex_z_yyz_0, tex_z_yyz_1, \
                                         tex_z_yz_0, tex_z_yz_1, tex_z_yzz_0, tex_z_yzz_1, tex_z_zz_0, tex_z_zz_1, \
                                         tex_z_zzz_0, tex_z_zzz_1, tey_xy_xzz_0, tey_xy_yyy_0, tey_xy_yyz_0, tey_xy_yzz_0, \
                                         tey_xy_zzz_0, tey_xz_xxx_0, tey_xz_xxy_0, tey_xz_xxz_0, tey_xz_xyy_0, tey_xz_xyz_0, \
                                         tey_xz_xzz_0, tey_xz_yyy_0, tey_xz_yyz_0, tey_xz_yzz_0, tey_xz_zzz_0, tey_y_xzz_0, \
                                         tey_y_xzz_1, tey_y_yyy_0, tey_y_yyy_1, tey_y_yyz_0, tey_y_yyz_1, tey_y_yzz_0, \
                                         tey_y_yzz_1, tey_y_zz_0, tey_y_zz_1, tey_y_zzz_0, tey_y_zzz_1, tey_z_xx_0, \
                                         tey_z_xx_1, tey_z_xxx_0, tey_z_xxx_1, tey_z_xxy_0, tey_z_xxy_1, tey_z_xxz_0, \
                                         tey_z_xxz_1, tey_z_xy_0, tey_z_xy_1, tey_z_xyy_0, tey_z_xyy_1, tey_z_xyz_0, \
                                         tey_z_xyz_1, tey_z_xz_0, tey_z_xz_1, tey_z_xzz_0, tey_z_xzz_1, tey_z_yy_0, \
                                         tey_z_yy_1, tey_z_yyy_0, tey_z_yyy_1, tey_z_yyz_0, tey_z_yyz_1, tey_z_yz_0, \
                                         tey_z_yz_1, tey_z_yzz_0, tey_z_yzz_1, tey_z_zz_0, tey_z_zz_1, tey_z_zzz_0, \
                                         tey_z_zzz_1, tez_xy_xzz_0, tez_xy_yyy_0, tez_xy_yyz_0, tez_xy_yzz_0, tez_xy_zzz_0, \
                                         tez_xz_xxx_0, tez_xz_xxy_0, tez_xz_xxz_0, tez_xz_xyy_0, tez_xz_xyz_0, tez_xz_xzz_0, \
                                         tez_xz_yyy_0, tez_xz_yyz_0, tez_xz_yzz_0, tez_xz_zzz_0, tez_y_xzz_0, tez_y_xzz_1, \
                                         tez_y_yyy_0, tez_y_yyy_1, tez_y_yyz_0, tez_y_yyz_1, tez_y_yzz_0, tez_y_yzz_1, \
                                         tez_y_zz_0, tez_y_zz_1, tez_y_zzz_0, tez_y_zzz_1, tez_z_xx_0, tez_z_xx_1, \
                                         tez_z_xxx_0, tez_z_xxx_1, tez_z_xxy_0, tez_z_xxy_1, tez_z_xxz_0, tez_z_xxz_1, \
                                         tez_z_xy_0, tez_z_xy_1, tez_z_xyy_0, tez_z_xyy_1, tez_z_xyz_0, tez_z_xyz_1, \
                                         tez_z_xz_0, tez_z_xz_1, tez_z_xzz_0, tez_z_xzz_1, tez_z_yy_0, tez_z_yy_1, \
                                         tez_z_yyy_0, tez_z_yyy_1, tez_z_yyz_0, tez_z_yyz_1, tez_z_yz_0, tez_z_yz_1, \
                                         tez_z_yzz_0, tez_z_yzz_1, tez_z_zz_0, tez_z_zz_1, tez_z_zzz_0, tez_z_zzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xy_xzz_0[j] = pa_x[j] * tex_y_xzz_0[j] - pc_x[j] * tex_y_xzz_1[j] + 0.5 * fl1_fx * tex_y_zz_0[j] - 0.5 * fl1_fx * tex_y_zz_1[j] + ta_y_xzz_1[j];

                    tey_xy_xzz_0[j] = pa_x[j] * tey_y_xzz_0[j] - pc_x[j] * tey_y_xzz_1[j] + 0.5 * fl1_fx * tey_y_zz_0[j] - 0.5 * fl1_fx * tey_y_zz_1[j];

                    tez_xy_xzz_0[j] = pa_x[j] * tez_y_xzz_0[j] - pc_x[j] * tez_y_xzz_1[j] + 0.5 * fl1_fx * tez_y_zz_0[j] - 0.5 * fl1_fx * tez_y_zz_1[j];

                    tex_xy_yyy_0[j] = pa_x[j] * tex_y_yyy_0[j] - pc_x[j] * tex_y_yyy_1[j] + ta_y_yyy_1[j];

                    tey_xy_yyy_0[j] = pa_x[j] * tey_y_yyy_0[j] - pc_x[j] * tey_y_yyy_1[j];

                    tez_xy_yyy_0[j] = pa_x[j] * tez_y_yyy_0[j] - pc_x[j] * tez_y_yyy_1[j];

                    tex_xy_yyz_0[j] = pa_x[j] * tex_y_yyz_0[j] - pc_x[j] * tex_y_yyz_1[j] + ta_y_yyz_1[j];

                    tey_xy_yyz_0[j] = pa_x[j] * tey_y_yyz_0[j] - pc_x[j] * tey_y_yyz_1[j];

                    tez_xy_yyz_0[j] = pa_x[j] * tez_y_yyz_0[j] - pc_x[j] * tez_y_yyz_1[j];

                    tex_xy_yzz_0[j] = pa_x[j] * tex_y_yzz_0[j] - pc_x[j] * tex_y_yzz_1[j] + ta_y_yzz_1[j];

                    tey_xy_yzz_0[j] = pa_x[j] * tey_y_yzz_0[j] - pc_x[j] * tey_y_yzz_1[j];

                    tez_xy_yzz_0[j] = pa_x[j] * tez_y_yzz_0[j] - pc_x[j] * tez_y_yzz_1[j];

                    tex_xy_zzz_0[j] = pa_x[j] * tex_y_zzz_0[j] - pc_x[j] * tex_y_zzz_1[j] + ta_y_zzz_1[j];

                    tey_xy_zzz_0[j] = pa_x[j] * tey_y_zzz_0[j] - pc_x[j] * tey_y_zzz_1[j];

                    tez_xy_zzz_0[j] = pa_x[j] * tez_y_zzz_0[j] - pc_x[j] * tez_y_zzz_1[j];

                    tex_xz_xxx_0[j] = pa_x[j] * tex_z_xxx_0[j] - pc_x[j] * tex_z_xxx_1[j] + 1.5 * fl1_fx * tex_z_xx_0[j] - 1.5 * fl1_fx * tex_z_xx_1[j] + ta_z_xxx_1[j];

                    tey_xz_xxx_0[j] = pa_x[j] * tey_z_xxx_0[j] - pc_x[j] * tey_z_xxx_1[j] + 1.5 * fl1_fx * tey_z_xx_0[j] - 1.5 * fl1_fx * tey_z_xx_1[j];

                    tez_xz_xxx_0[j] = pa_x[j] * tez_z_xxx_0[j] - pc_x[j] * tez_z_xxx_1[j] + 1.5 * fl1_fx * tez_z_xx_0[j] - 1.5 * fl1_fx * tez_z_xx_1[j];

                    tex_xz_xxy_0[j] = pa_x[j] * tex_z_xxy_0[j] - pc_x[j] * tex_z_xxy_1[j] + fl1_fx * tex_z_xy_0[j] - fl1_fx * tex_z_xy_1[j] + ta_z_xxy_1[j];

                    tey_xz_xxy_0[j] = pa_x[j] * tey_z_xxy_0[j] - pc_x[j] * tey_z_xxy_1[j] + fl1_fx * tey_z_xy_0[j] - fl1_fx * tey_z_xy_1[j];

                    tez_xz_xxy_0[j] = pa_x[j] * tez_z_xxy_0[j] - pc_x[j] * tez_z_xxy_1[j] + fl1_fx * tez_z_xy_0[j] - fl1_fx * tez_z_xy_1[j];

                    tex_xz_xxz_0[j] = pa_x[j] * tex_z_xxz_0[j] - pc_x[j] * tex_z_xxz_1[j] + fl1_fx * tex_z_xz_0[j] - fl1_fx * tex_z_xz_1[j] + ta_z_xxz_1[j];

                    tey_xz_xxz_0[j] = pa_x[j] * tey_z_xxz_0[j] - pc_x[j] * tey_z_xxz_1[j] + fl1_fx * tey_z_xz_0[j] - fl1_fx * tey_z_xz_1[j];

                    tez_xz_xxz_0[j] = pa_x[j] * tez_z_xxz_0[j] - pc_x[j] * tez_z_xxz_1[j] + fl1_fx * tez_z_xz_0[j] - fl1_fx * tez_z_xz_1[j];

                    tex_xz_xyy_0[j] = pa_x[j] * tex_z_xyy_0[j] - pc_x[j] * tex_z_xyy_1[j] + 0.5 * fl1_fx * tex_z_yy_0[j] - 0.5 * fl1_fx * tex_z_yy_1[j] + ta_z_xyy_1[j];

                    tey_xz_xyy_0[j] = pa_x[j] * tey_z_xyy_0[j] - pc_x[j] * tey_z_xyy_1[j] + 0.5 * fl1_fx * tey_z_yy_0[j] - 0.5 * fl1_fx * tey_z_yy_1[j];

                    tez_xz_xyy_0[j] = pa_x[j] * tez_z_xyy_0[j] - pc_x[j] * tez_z_xyy_1[j] + 0.5 * fl1_fx * tez_z_yy_0[j] - 0.5 * fl1_fx * tez_z_yy_1[j];

                    tex_xz_xyz_0[j] = pa_x[j] * tex_z_xyz_0[j] - pc_x[j] * tex_z_xyz_1[j] + 0.5 * fl1_fx * tex_z_yz_0[j] - 0.5 * fl1_fx * tex_z_yz_1[j] + ta_z_xyz_1[j];

                    tey_xz_xyz_0[j] = pa_x[j] * tey_z_xyz_0[j] - pc_x[j] * tey_z_xyz_1[j] + 0.5 * fl1_fx * tey_z_yz_0[j] - 0.5 * fl1_fx * tey_z_yz_1[j];

                    tez_xz_xyz_0[j] = pa_x[j] * tez_z_xyz_0[j] - pc_x[j] * tez_z_xyz_1[j] + 0.5 * fl1_fx * tez_z_yz_0[j] - 0.5 * fl1_fx * tez_z_yz_1[j];

                    tex_xz_xzz_0[j] = pa_x[j] * tex_z_xzz_0[j] - pc_x[j] * tex_z_xzz_1[j] + 0.5 * fl1_fx * tex_z_zz_0[j] - 0.5 * fl1_fx * tex_z_zz_1[j] + ta_z_xzz_1[j];

                    tey_xz_xzz_0[j] = pa_x[j] * tey_z_xzz_0[j] - pc_x[j] * tey_z_xzz_1[j] + 0.5 * fl1_fx * tey_z_zz_0[j] - 0.5 * fl1_fx * tey_z_zz_1[j];

                    tez_xz_xzz_0[j] = pa_x[j] * tez_z_xzz_0[j] - pc_x[j] * tez_z_xzz_1[j] + 0.5 * fl1_fx * tez_z_zz_0[j] - 0.5 * fl1_fx * tez_z_zz_1[j];

                    tex_xz_yyy_0[j] = pa_x[j] * tex_z_yyy_0[j] - pc_x[j] * tex_z_yyy_1[j] + ta_z_yyy_1[j];

                    tey_xz_yyy_0[j] = pa_x[j] * tey_z_yyy_0[j] - pc_x[j] * tey_z_yyy_1[j];

                    tez_xz_yyy_0[j] = pa_x[j] * tez_z_yyy_0[j] - pc_x[j] * tez_z_yyy_1[j];

                    tex_xz_yyz_0[j] = pa_x[j] * tex_z_yyz_0[j] - pc_x[j] * tex_z_yyz_1[j] + ta_z_yyz_1[j];

                    tey_xz_yyz_0[j] = pa_x[j] * tey_z_yyz_0[j] - pc_x[j] * tey_z_yyz_1[j];

                    tez_xz_yyz_0[j] = pa_x[j] * tez_z_yyz_0[j] - pc_x[j] * tez_z_yyz_1[j];

                    tex_xz_yzz_0[j] = pa_x[j] * tex_z_yzz_0[j] - pc_x[j] * tex_z_yzz_1[j] + ta_z_yzz_1[j];

                    tey_xz_yzz_0[j] = pa_x[j] * tey_z_yzz_0[j] - pc_x[j] * tey_z_yzz_1[j];

                    tez_xz_yzz_0[j] = pa_x[j] * tez_z_yzz_0[j] - pc_x[j] * tez_z_yzz_1[j];

                    tex_xz_zzz_0[j] = pa_x[j] * tex_z_zzz_0[j] - pc_x[j] * tex_z_zzz_1[j] + ta_z_zzz_1[j];

                    tey_xz_zzz_0[j] = pa_x[j] * tey_z_zzz_0[j] - pc_x[j] * tey_z_zzz_1[j];

                    tez_xz_zzz_0[j] = pa_x[j] * tez_z_zzz_0[j] - pc_x[j] * tez_z_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDF_90_135(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (90,135)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tex_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 10); 

                auto tey_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 11); 

                auto tey_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 12); 

                auto tey_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 13); 

                auto tey_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 14); 

                auto tey_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 14); 

                auto tex_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 15); 

                auto tey_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 16); 

                auto tey_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 17); 

                auto tey_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 18); 

                auto tey_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 19); 

                auto tey_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 19); 

                auto tex_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 20); 

                auto tey_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 21); 

                auto tey_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 22); 

                auto tey_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 23); 

                auto tey_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 24); 

                auto tey_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 10); 

                auto tey_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 11); 

                auto tey_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 12); 

                auto tey_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 13); 

                auto tey_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 14); 

                auto tey_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 14); 

                auto tex_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 15); 

                auto tey_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 16); 

                auto tey_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 17); 

                auto tey_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 18); 

                auto tey_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 19); 

                auto tey_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 19); 

                auto tex_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 20); 

                auto tey_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 21); 

                auto tey_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 22); 

                auto tey_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 23); 

                auto tey_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 24); 

                auto tey_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 24); 

                auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx); 

                auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx); 

                auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx); 

                auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1); 

                auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1); 

                auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1); 

                auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2); 

                auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2); 

                auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2); 

                auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3); 

                auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3); 

                auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3); 

                auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4); 

                auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4); 

                auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4); 

                auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5); 

                auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5); 

                auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5); 

                auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6); 

                auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6); 

                auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6); 

                auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7); 

                auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7); 

                auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7); 

                auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8); 

                auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8); 

                auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8); 

                auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9); 

                auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9); 

                auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9); 

                auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx); 

                auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx); 

                auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx); 

                auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1); 

                auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1); 

                auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1); 

                auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2); 

                auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2); 

                auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2); 

                auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3); 

                auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3); 

                auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3); 

                auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4); 

                auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4); 

                auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4); 

                auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5); 

                auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5); 

                auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5); 

                auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6); 

                auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6); 

                auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6); 

                auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7); 

                auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7); 

                auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7); 

                auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8); 

                auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8); 

                auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8); 

                auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9); 

                auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9); 

                auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9); 

                auto tex_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 6); 

                auto tey_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 7); 

                auto tey_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 8); 

                auto tey_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 9); 

                auto tey_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 10); 

                auto tey_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 11); 

                auto tey_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 6); 

                auto tey_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 7); 

                auto tey_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 8); 

                auto tey_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 9); 

                auto tey_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 10); 

                auto tey_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 11); 

                auto tey_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10); 

                auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11); 

                auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12); 

                auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13); 

                auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14); 

                auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15); 

                auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16); 

                auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17); 

                auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18); 

                auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19); 

                auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20); 

                auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21); 

                auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22); 

                auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23); 

                auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24); 

                // set up pointers to integrals

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

                // Batch of Integrals (90,135)

                #pragma omp simd aligned(fx, pa_y, pc_y, ta_y_xxx_1, ta_y_xxy_1, ta_y_xxz_1, ta_y_xyy_1, ta_y_xyz_1, \
                                         ta_y_xzz_1, ta_y_yyy_1, ta_y_yyz_1, ta_y_yzz_1, ta_y_zzz_1, ta_z_xxx_1, ta_z_xxy_1, \
                                         ta_z_xxz_1, ta_z_xyy_1, ta_z_xyz_1, tex_0_xxx_0, tex_0_xxx_1, tex_0_xxy_0, \
                                         tex_0_xxy_1, tex_0_xxz_0, tex_0_xxz_1, tex_0_xyy_0, tex_0_xyy_1, tex_0_xyz_0, \
                                         tex_0_xyz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_yyy_0, tex_0_yyy_1, tex_0_yyz_0, \
                                         tex_0_yyz_1, tex_0_yzz_0, tex_0_yzz_1, tex_0_zzz_0, tex_0_zzz_1, tex_y_xx_0, \
                                         tex_y_xx_1, tex_y_xxx_0, tex_y_xxx_1, tex_y_xxy_0, tex_y_xxy_1, tex_y_xxz_0, \
                                         tex_y_xxz_1, tex_y_xy_0, tex_y_xy_1, tex_y_xyy_0, tex_y_xyy_1, tex_y_xyz_0, \
                                         tex_y_xyz_1, tex_y_xz_0, tex_y_xz_1, tex_y_xzz_0, tex_y_xzz_1, tex_y_yy_0, \
                                         tex_y_yy_1, tex_y_yyy_0, tex_y_yyy_1, tex_y_yyz_0, tex_y_yyz_1, tex_y_yz_0, \
                                         tex_y_yz_1, tex_y_yzz_0, tex_y_yzz_1, tex_y_zz_0, tex_y_zz_1, tex_y_zzz_0, \
                                         tex_y_zzz_1, tex_yy_xxx_0, tex_yy_xxy_0, tex_yy_xxz_0, tex_yy_xyy_0, tex_yy_xyz_0, \
                                         tex_yy_xzz_0, tex_yy_yyy_0, tex_yy_yyz_0, tex_yy_yzz_0, tex_yy_zzz_0, tex_yz_xxx_0, \
                                         tex_yz_xxy_0, tex_yz_xxz_0, tex_yz_xyy_0, tex_yz_xyz_0, tex_z_xx_0, tex_z_xx_1, \
                                         tex_z_xxx_0, tex_z_xxx_1, tex_z_xxy_0, tex_z_xxy_1, tex_z_xxz_0, tex_z_xxz_1, \
                                         tex_z_xy_0, tex_z_xy_1, tex_z_xyy_0, tex_z_xyy_1, tex_z_xyz_0, tex_z_xyz_1, \
                                         tex_z_xz_0, tex_z_xz_1, tey_0_xxx_0, tey_0_xxx_1, tey_0_xxy_0, tey_0_xxy_1, \
                                         tey_0_xxz_0, tey_0_xxz_1, tey_0_xyy_0, tey_0_xyy_1, tey_0_xyz_0, tey_0_xyz_1, \
                                         tey_0_xzz_0, tey_0_xzz_1, tey_0_yyy_0, tey_0_yyy_1, tey_0_yyz_0, tey_0_yyz_1, \
                                         tey_0_yzz_0, tey_0_yzz_1, tey_0_zzz_0, tey_0_zzz_1, tey_y_xx_0, tey_y_xx_1, \
                                         tey_y_xxx_0, tey_y_xxx_1, tey_y_xxy_0, tey_y_xxy_1, tey_y_xxz_0, tey_y_xxz_1, \
                                         tey_y_xy_0, tey_y_xy_1, tey_y_xyy_0, tey_y_xyy_1, tey_y_xyz_0, tey_y_xyz_1, \
                                         tey_y_xz_0, tey_y_xz_1, tey_y_xzz_0, tey_y_xzz_1, tey_y_yy_0, tey_y_yy_1, \
                                         tey_y_yyy_0, tey_y_yyy_1, tey_y_yyz_0, tey_y_yyz_1, tey_y_yz_0, tey_y_yz_1, \
                                         tey_y_yzz_0, tey_y_yzz_1, tey_y_zz_0, tey_y_zz_1, tey_y_zzz_0, tey_y_zzz_1, \
                                         tey_yy_xxx_0, tey_yy_xxy_0, tey_yy_xxz_0, tey_yy_xyy_0, tey_yy_xyz_0, tey_yy_xzz_0, \
                                         tey_yy_yyy_0, tey_yy_yyz_0, tey_yy_yzz_0, tey_yy_zzz_0, tey_yz_xxx_0, tey_yz_xxy_0, \
                                         tey_yz_xxz_0, tey_yz_xyy_0, tey_yz_xyz_0, tey_z_xx_0, tey_z_xx_1, tey_z_xxx_0, \
                                         tey_z_xxx_1, tey_z_xxy_0, tey_z_xxy_1, tey_z_xxz_0, tey_z_xxz_1, tey_z_xy_0, \
                                         tey_z_xy_1, tey_z_xyy_0, tey_z_xyy_1, tey_z_xyz_0, tey_z_xyz_1, tey_z_xz_0, \
                                         tey_z_xz_1, tez_0_xxx_0, tez_0_xxx_1, tez_0_xxy_0, tez_0_xxy_1, tez_0_xxz_0, \
                                         tez_0_xxz_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyz_0, tez_0_xyz_1, tez_0_xzz_0, \
                                         tez_0_xzz_1, tez_0_yyy_0, tez_0_yyy_1, tez_0_yyz_0, tez_0_yyz_1, tez_0_yzz_0, \
                                         tez_0_yzz_1, tez_0_zzz_0, tez_0_zzz_1, tez_y_xx_0, tez_y_xx_1, tez_y_xxx_0, \
                                         tez_y_xxx_1, tez_y_xxy_0, tez_y_xxy_1, tez_y_xxz_0, tez_y_xxz_1, tez_y_xy_0, \
                                         tez_y_xy_1, tez_y_xyy_0, tez_y_xyy_1, tez_y_xyz_0, tez_y_xyz_1, tez_y_xz_0, \
                                         tez_y_xz_1, tez_y_xzz_0, tez_y_xzz_1, tez_y_yy_0, tez_y_yy_1, tez_y_yyy_0, \
                                         tez_y_yyy_1, tez_y_yyz_0, tez_y_yyz_1, tez_y_yz_0, tez_y_yz_1, tez_y_yzz_0, \
                                         tez_y_yzz_1, tez_y_zz_0, tez_y_zz_1, tez_y_zzz_0, tez_y_zzz_1, tez_yy_xxx_0, \
                                         tez_yy_xxy_0, tez_yy_xxz_0, tez_yy_xyy_0, tez_yy_xyz_0, tez_yy_xzz_0, tez_yy_yyy_0, \
                                         tez_yy_yyz_0, tez_yy_yzz_0, tez_yy_zzz_0, tez_yz_xxx_0, tez_yz_xxy_0, tez_yz_xxz_0, \
                                         tez_yz_xyy_0, tez_yz_xyz_0, tez_z_xx_0, tez_z_xx_1, tez_z_xxx_0, tez_z_xxx_1, \
                                         tez_z_xxy_0, tez_z_xxy_1, tez_z_xxz_0, tez_z_xxz_1, tez_z_xy_0, tez_z_xy_1, \
                                         tez_z_xyy_0, tez_z_xyy_1, tez_z_xyz_0, tez_z_xyz_1, tez_z_xz_0, tez_z_xz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yy_xxx_0[j] = pa_y[j] * tex_y_xxx_0[j] - pc_y[j] * tex_y_xxx_1[j] + 0.5 * fl1_fx * tex_0_xxx_0[j] - 0.5 * fl1_fx * tex_0_xxx_1[j];

                    tey_yy_xxx_0[j] = pa_y[j] * tey_y_xxx_0[j] - pc_y[j] * tey_y_xxx_1[j] + 0.5 * fl1_fx * tey_0_xxx_0[j] - 0.5 * fl1_fx * tey_0_xxx_1[j] + ta_y_xxx_1[j];

                    tez_yy_xxx_0[j] = pa_y[j] * tez_y_xxx_0[j] - pc_y[j] * tez_y_xxx_1[j] + 0.5 * fl1_fx * tez_0_xxx_0[j] - 0.5 * fl1_fx * tez_0_xxx_1[j];

                    tex_yy_xxy_0[j] = pa_y[j] * tex_y_xxy_0[j] - pc_y[j] * tex_y_xxy_1[j] + 0.5 * fl1_fx * tex_0_xxy_0[j] - 0.5 * fl1_fx * tex_0_xxy_1[j] + 0.5 * fl1_fx * tex_y_xx_0[j] - 0.5 * fl1_fx * tex_y_xx_1[j];

                    tey_yy_xxy_0[j] = pa_y[j] * tey_y_xxy_0[j] - pc_y[j] * tey_y_xxy_1[j] + 0.5 * fl1_fx * tey_0_xxy_0[j] - 0.5 * fl1_fx * tey_0_xxy_1[j] + 0.5 * fl1_fx * tey_y_xx_0[j] - 0.5 * fl1_fx * tey_y_xx_1[j] + ta_y_xxy_1[j];

                    tez_yy_xxy_0[j] = pa_y[j] * tez_y_xxy_0[j] - pc_y[j] * tez_y_xxy_1[j] + 0.5 * fl1_fx * tez_0_xxy_0[j] - 0.5 * fl1_fx * tez_0_xxy_1[j] + 0.5 * fl1_fx * tez_y_xx_0[j] - 0.5 * fl1_fx * tez_y_xx_1[j];

                    tex_yy_xxz_0[j] = pa_y[j] * tex_y_xxz_0[j] - pc_y[j] * tex_y_xxz_1[j] + 0.5 * fl1_fx * tex_0_xxz_0[j] - 0.5 * fl1_fx * tex_0_xxz_1[j];

                    tey_yy_xxz_0[j] = pa_y[j] * tey_y_xxz_0[j] - pc_y[j] * tey_y_xxz_1[j] + 0.5 * fl1_fx * tey_0_xxz_0[j] - 0.5 * fl1_fx * tey_0_xxz_1[j] + ta_y_xxz_1[j];

                    tez_yy_xxz_0[j] = pa_y[j] * tez_y_xxz_0[j] - pc_y[j] * tez_y_xxz_1[j] + 0.5 * fl1_fx * tez_0_xxz_0[j] - 0.5 * fl1_fx * tez_0_xxz_1[j];

                    tex_yy_xyy_0[j] = pa_y[j] * tex_y_xyy_0[j] - pc_y[j] * tex_y_xyy_1[j] + 0.5 * fl1_fx * tex_0_xyy_0[j] - 0.5 * fl1_fx * tex_0_xyy_1[j] + fl1_fx * tex_y_xy_0[j] - fl1_fx * tex_y_xy_1[j];

                    tey_yy_xyy_0[j] = pa_y[j] * tey_y_xyy_0[j] - pc_y[j] * tey_y_xyy_1[j] + 0.5 * fl1_fx * tey_0_xyy_0[j] - 0.5 * fl1_fx * tey_0_xyy_1[j] + fl1_fx * tey_y_xy_0[j] - fl1_fx * tey_y_xy_1[j] + ta_y_xyy_1[j];

                    tez_yy_xyy_0[j] = pa_y[j] * tez_y_xyy_0[j] - pc_y[j] * tez_y_xyy_1[j] + 0.5 * fl1_fx * tez_0_xyy_0[j] - 0.5 * fl1_fx * tez_0_xyy_1[j] + fl1_fx * tez_y_xy_0[j] - fl1_fx * tez_y_xy_1[j];

                    tex_yy_xyz_0[j] = pa_y[j] * tex_y_xyz_0[j] - pc_y[j] * tex_y_xyz_1[j] + 0.5 * fl1_fx * tex_0_xyz_0[j] - 0.5 * fl1_fx * tex_0_xyz_1[j] + 0.5 * fl1_fx * tex_y_xz_0[j] - 0.5 * fl1_fx * tex_y_xz_1[j];

                    tey_yy_xyz_0[j] = pa_y[j] * tey_y_xyz_0[j] - pc_y[j] * tey_y_xyz_1[j] + 0.5 * fl1_fx * tey_0_xyz_0[j] - 0.5 * fl1_fx * tey_0_xyz_1[j] + 0.5 * fl1_fx * tey_y_xz_0[j] - 0.5 * fl1_fx * tey_y_xz_1[j] + ta_y_xyz_1[j];

                    tez_yy_xyz_0[j] = pa_y[j] * tez_y_xyz_0[j] - pc_y[j] * tez_y_xyz_1[j] + 0.5 * fl1_fx * tez_0_xyz_0[j] - 0.5 * fl1_fx * tez_0_xyz_1[j] + 0.5 * fl1_fx * tez_y_xz_0[j] - 0.5 * fl1_fx * tez_y_xz_1[j];

                    tex_yy_xzz_0[j] = pa_y[j] * tex_y_xzz_0[j] - pc_y[j] * tex_y_xzz_1[j] + 0.5 * fl1_fx * tex_0_xzz_0[j] - 0.5 * fl1_fx * tex_0_xzz_1[j];

                    tey_yy_xzz_0[j] = pa_y[j] * tey_y_xzz_0[j] - pc_y[j] * tey_y_xzz_1[j] + 0.5 * fl1_fx * tey_0_xzz_0[j] - 0.5 * fl1_fx * tey_0_xzz_1[j] + ta_y_xzz_1[j];

                    tez_yy_xzz_0[j] = pa_y[j] * tez_y_xzz_0[j] - pc_y[j] * tez_y_xzz_1[j] + 0.5 * fl1_fx * tez_0_xzz_0[j] - 0.5 * fl1_fx * tez_0_xzz_1[j];

                    tex_yy_yyy_0[j] = pa_y[j] * tex_y_yyy_0[j] - pc_y[j] * tex_y_yyy_1[j] + 0.5 * fl1_fx * tex_0_yyy_0[j] - 0.5 * fl1_fx * tex_0_yyy_1[j] + 1.5 * fl1_fx * tex_y_yy_0[j] - 1.5 * fl1_fx * tex_y_yy_1[j];

                    tey_yy_yyy_0[j] = pa_y[j] * tey_y_yyy_0[j] - pc_y[j] * tey_y_yyy_1[j] + 0.5 * fl1_fx * tey_0_yyy_0[j] - 0.5 * fl1_fx * tey_0_yyy_1[j] + 1.5 * fl1_fx * tey_y_yy_0[j] - 1.5 * fl1_fx * tey_y_yy_1[j] + ta_y_yyy_1[j];

                    tez_yy_yyy_0[j] = pa_y[j] * tez_y_yyy_0[j] - pc_y[j] * tez_y_yyy_1[j] + 0.5 * fl1_fx * tez_0_yyy_0[j] - 0.5 * fl1_fx * tez_0_yyy_1[j] + 1.5 * fl1_fx * tez_y_yy_0[j] - 1.5 * fl1_fx * tez_y_yy_1[j];

                    tex_yy_yyz_0[j] = pa_y[j] * tex_y_yyz_0[j] - pc_y[j] * tex_y_yyz_1[j] + 0.5 * fl1_fx * tex_0_yyz_0[j] - 0.5 * fl1_fx * tex_0_yyz_1[j] + fl1_fx * tex_y_yz_0[j] - fl1_fx * tex_y_yz_1[j];

                    tey_yy_yyz_0[j] = pa_y[j] * tey_y_yyz_0[j] - pc_y[j] * tey_y_yyz_1[j] + 0.5 * fl1_fx * tey_0_yyz_0[j] - 0.5 * fl1_fx * tey_0_yyz_1[j] + fl1_fx * tey_y_yz_0[j] - fl1_fx * tey_y_yz_1[j] + ta_y_yyz_1[j];

                    tez_yy_yyz_0[j] = pa_y[j] * tez_y_yyz_0[j] - pc_y[j] * tez_y_yyz_1[j] + 0.5 * fl1_fx * tez_0_yyz_0[j] - 0.5 * fl1_fx * tez_0_yyz_1[j] + fl1_fx * tez_y_yz_0[j] - fl1_fx * tez_y_yz_1[j];

                    tex_yy_yzz_0[j] = pa_y[j] * tex_y_yzz_0[j] - pc_y[j] * tex_y_yzz_1[j] + 0.5 * fl1_fx * tex_0_yzz_0[j] - 0.5 * fl1_fx * tex_0_yzz_1[j] + 0.5 * fl1_fx * tex_y_zz_0[j] - 0.5 * fl1_fx * tex_y_zz_1[j];

                    tey_yy_yzz_0[j] = pa_y[j] * tey_y_yzz_0[j] - pc_y[j] * tey_y_yzz_1[j] + 0.5 * fl1_fx * tey_0_yzz_0[j] - 0.5 * fl1_fx * tey_0_yzz_1[j] + 0.5 * fl1_fx * tey_y_zz_0[j] - 0.5 * fl1_fx * tey_y_zz_1[j] + ta_y_yzz_1[j];

                    tez_yy_yzz_0[j] = pa_y[j] * tez_y_yzz_0[j] - pc_y[j] * tez_y_yzz_1[j] + 0.5 * fl1_fx * tez_0_yzz_0[j] - 0.5 * fl1_fx * tez_0_yzz_1[j] + 0.5 * fl1_fx * tez_y_zz_0[j] - 0.5 * fl1_fx * tez_y_zz_1[j];

                    tex_yy_zzz_0[j] = pa_y[j] * tex_y_zzz_0[j] - pc_y[j] * tex_y_zzz_1[j] + 0.5 * fl1_fx * tex_0_zzz_0[j] - 0.5 * fl1_fx * tex_0_zzz_1[j];

                    tey_yy_zzz_0[j] = pa_y[j] * tey_y_zzz_0[j] - pc_y[j] * tey_y_zzz_1[j] + 0.5 * fl1_fx * tey_0_zzz_0[j] - 0.5 * fl1_fx * tey_0_zzz_1[j] + ta_y_zzz_1[j];

                    tez_yy_zzz_0[j] = pa_y[j] * tez_y_zzz_0[j] - pc_y[j] * tez_y_zzz_1[j] + 0.5 * fl1_fx * tez_0_zzz_0[j] - 0.5 * fl1_fx * tez_0_zzz_1[j];

                    tex_yz_xxx_0[j] = pa_y[j] * tex_z_xxx_0[j] - pc_y[j] * tex_z_xxx_1[j];

                    tey_yz_xxx_0[j] = pa_y[j] * tey_z_xxx_0[j] - pc_y[j] * tey_z_xxx_1[j] + ta_z_xxx_1[j];

                    tez_yz_xxx_0[j] = pa_y[j] * tez_z_xxx_0[j] - pc_y[j] * tez_z_xxx_1[j];

                    tex_yz_xxy_0[j] = pa_y[j] * tex_z_xxy_0[j] - pc_y[j] * tex_z_xxy_1[j] + 0.5 * fl1_fx * tex_z_xx_0[j] - 0.5 * fl1_fx * tex_z_xx_1[j];

                    tey_yz_xxy_0[j] = pa_y[j] * tey_z_xxy_0[j] - pc_y[j] * tey_z_xxy_1[j] + 0.5 * fl1_fx * tey_z_xx_0[j] - 0.5 * fl1_fx * tey_z_xx_1[j] + ta_z_xxy_1[j];

                    tez_yz_xxy_0[j] = pa_y[j] * tez_z_xxy_0[j] - pc_y[j] * tez_z_xxy_1[j] + 0.5 * fl1_fx * tez_z_xx_0[j] - 0.5 * fl1_fx * tez_z_xx_1[j];

                    tex_yz_xxz_0[j] = pa_y[j] * tex_z_xxz_0[j] - pc_y[j] * tex_z_xxz_1[j];

                    tey_yz_xxz_0[j] = pa_y[j] * tey_z_xxz_0[j] - pc_y[j] * tey_z_xxz_1[j] + ta_z_xxz_1[j];

                    tez_yz_xxz_0[j] = pa_y[j] * tez_z_xxz_0[j] - pc_y[j] * tez_z_xxz_1[j];

                    tex_yz_xyy_0[j] = pa_y[j] * tex_z_xyy_0[j] - pc_y[j] * tex_z_xyy_1[j] + fl1_fx * tex_z_xy_0[j] - fl1_fx * tex_z_xy_1[j];

                    tey_yz_xyy_0[j] = pa_y[j] * tey_z_xyy_0[j] - pc_y[j] * tey_z_xyy_1[j] + fl1_fx * tey_z_xy_0[j] - fl1_fx * tey_z_xy_1[j] + ta_z_xyy_1[j];

                    tez_yz_xyy_0[j] = pa_y[j] * tez_z_xyy_0[j] - pc_y[j] * tez_z_xyy_1[j] + fl1_fx * tez_z_xy_0[j] - fl1_fx * tez_z_xy_1[j];

                    tex_yz_xyz_0[j] = pa_y[j] * tex_z_xyz_0[j] - pc_y[j] * tex_z_xyz_1[j] + 0.5 * fl1_fx * tex_z_xz_0[j] - 0.5 * fl1_fx * tex_z_xz_1[j];

                    tey_yz_xyz_0[j] = pa_y[j] * tey_z_xyz_0[j] - pc_y[j] * tey_z_xyz_1[j] + 0.5 * fl1_fx * tey_z_xz_0[j] - 0.5 * fl1_fx * tey_z_xz_1[j] + ta_z_xyz_1[j];

                    tez_yz_xyz_0[j] = pa_y[j] * tez_z_xyz_0[j] - pc_y[j] * tez_z_xyz_1[j] + 0.5 * fl1_fx * tez_z_xz_0[j] - 0.5 * fl1_fx * tez_z_xz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDF_135_180(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (135,180)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {3, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_3_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tex_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 20); 

                auto tey_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 21); 

                auto tey_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 22); 

                auto tey_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 23); 

                auto tey_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 24); 

                auto tey_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 25); 

                auto tey_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 26); 

                auto tey_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 27); 

                auto tey_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 28); 

                auto tey_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 29); 

                auto tey_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 29); 

                auto tex_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 20); 

                auto tey_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 21); 

                auto tey_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 22); 

                auto tey_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 23); 

                auto tey_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 24); 

                auto tey_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 25); 

                auto tey_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 26); 

                auto tey_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 27); 

                auto tey_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 28); 

                auto tey_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 29); 

                auto tey_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 29); 

                auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx); 

                auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx); 

                auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx); 

                auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1); 

                auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1); 

                auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1); 

                auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2); 

                auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2); 

                auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2); 

                auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3); 

                auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3); 

                auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3); 

                auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4); 

                auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4); 

                auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4); 

                auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5); 

                auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5); 

                auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5); 

                auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6); 

                auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6); 

                auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6); 

                auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7); 

                auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7); 

                auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7); 

                auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8); 

                auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8); 

                auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8); 

                auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9); 

                auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9); 

                auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9); 

                auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx); 

                auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx); 

                auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx); 

                auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1); 

                auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1); 

                auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1); 

                auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2); 

                auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2); 

                auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2); 

                auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3); 

                auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3); 

                auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3); 

                auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4); 

                auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4); 

                auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4); 

                auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5); 

                auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5); 

                auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5); 

                auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6); 

                auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6); 

                auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6); 

                auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7); 

                auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7); 

                auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7); 

                auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8); 

                auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8); 

                auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8); 

                auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9); 

                auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9); 

                auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9); 

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 15); 

                auto tey_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 16); 

                auto tey_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 17); 

                auto tey_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 15); 

                auto tey_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 16); 

                auto tey_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 17); 

                auto tey_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 17); 

                auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20); 

                auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21); 

                auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22); 

                auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23); 

                auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24); 

                auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25); 

                auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26); 

                auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27); 

                auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28); 

                auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29); 

                // set up pointers to integrals

                auto tex_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 45); 

                auto tey_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 45); 

                auto tez_yz_xzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 45); 

                auto tex_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 46); 

                auto tey_yz_yyy_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 46); 

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

                auto tey_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 58); 

                auto tez_zz_yzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 58); 

                auto tex_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * idx + 59); 

                auto tey_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 60 * bdim + 60 * idx + 59); 

                auto tez_zz_zzz_0 = primBuffer.data(pidx_e_2_3_m0 + 120 * bdim + 60 * idx + 59); 

                // Batch of Integrals (135,180)

                #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_z_xxx_1, ta_z_xxy_1, ta_z_xxz_1, ta_z_xyy_1, \
                                         ta_z_xyz_1, ta_z_xzz_1, ta_z_yyy_1, ta_z_yyz_1, ta_z_yzz_1, ta_z_zzz_1, \
                                         tex_0_xxx_0, tex_0_xxx_1, tex_0_xxy_0, tex_0_xxy_1, tex_0_xxz_0, tex_0_xxz_1, \
                                         tex_0_xyy_0, tex_0_xyy_1, tex_0_xyz_0, tex_0_xyz_1, tex_0_xzz_0, tex_0_xzz_1, \
                                         tex_0_yyy_0, tex_0_yyy_1, tex_0_yyz_0, tex_0_yyz_1, tex_0_yzz_0, tex_0_yzz_1, \
                                         tex_0_zzz_0, tex_0_zzz_1, tex_yz_xzz_0, tex_yz_yyy_0, tex_yz_yyz_0, tex_yz_yzz_0, \
                                         tex_yz_zzz_0, tex_z_xx_0, tex_z_xx_1, tex_z_xxx_0, tex_z_xxx_1, tex_z_xxy_0, \
                                         tex_z_xxy_1, tex_z_xxz_0, tex_z_xxz_1, tex_z_xy_0, tex_z_xy_1, tex_z_xyy_0, \
                                         tex_z_xyy_1, tex_z_xyz_0, tex_z_xyz_1, tex_z_xz_0, tex_z_xz_1, tex_z_xzz_0, \
                                         tex_z_xzz_1, tex_z_yy_0, tex_z_yy_1, tex_z_yyy_0, tex_z_yyy_1, tex_z_yyz_0, \
                                         tex_z_yyz_1, tex_z_yz_0, tex_z_yz_1, tex_z_yzz_0, tex_z_yzz_1, tex_z_zz_0, \
                                         tex_z_zz_1, tex_z_zzz_0, tex_z_zzz_1, tex_zz_xxx_0, tex_zz_xxy_0, tex_zz_xxz_0, \
                                         tex_zz_xyy_0, tex_zz_xyz_0, tex_zz_xzz_0, tex_zz_yyy_0, tex_zz_yyz_0, tex_zz_yzz_0, \
                                         tex_zz_zzz_0, tey_0_xxx_0, tey_0_xxx_1, tey_0_xxy_0, tey_0_xxy_1, tey_0_xxz_0, \
                                         tey_0_xxz_1, tey_0_xyy_0, tey_0_xyy_1, tey_0_xyz_0, tey_0_xyz_1, tey_0_xzz_0, \
                                         tey_0_xzz_1, tey_0_yyy_0, tey_0_yyy_1, tey_0_yyz_0, tey_0_yyz_1, tey_0_yzz_0, \
                                         tey_0_yzz_1, tey_0_zzz_0, tey_0_zzz_1, tey_yz_xzz_0, tey_yz_yyy_0, tey_yz_yyz_0, \
                                         tey_yz_yzz_0, tey_yz_zzz_0, tey_z_xx_0, tey_z_xx_1, tey_z_xxx_0, tey_z_xxx_1, \
                                         tey_z_xxy_0, tey_z_xxy_1, tey_z_xxz_0, tey_z_xxz_1, tey_z_xy_0, tey_z_xy_1, \
                                         tey_z_xyy_0, tey_z_xyy_1, tey_z_xyz_0, tey_z_xyz_1, tey_z_xz_0, tey_z_xz_1, \
                                         tey_z_xzz_0, tey_z_xzz_1, tey_z_yy_0, tey_z_yy_1, tey_z_yyy_0, tey_z_yyy_1, \
                                         tey_z_yyz_0, tey_z_yyz_1, tey_z_yz_0, tey_z_yz_1, tey_z_yzz_0, tey_z_yzz_1, \
                                         tey_z_zz_0, tey_z_zz_1, tey_z_zzz_0, tey_z_zzz_1, tey_zz_xxx_0, tey_zz_xxy_0, \
                                         tey_zz_xxz_0, tey_zz_xyy_0, tey_zz_xyz_0, tey_zz_xzz_0, tey_zz_yyy_0, tey_zz_yyz_0, \
                                         tey_zz_yzz_0, tey_zz_zzz_0, tez_0_xxx_0, tez_0_xxx_1, tez_0_xxy_0, tez_0_xxy_1, \
                                         tez_0_xxz_0, tez_0_xxz_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyz_0, tez_0_xyz_1, \
                                         tez_0_xzz_0, tez_0_xzz_1, tez_0_yyy_0, tez_0_yyy_1, tez_0_yyz_0, tez_0_yyz_1, \
                                         tez_0_yzz_0, tez_0_yzz_1, tez_0_zzz_0, tez_0_zzz_1, tez_yz_xzz_0, tez_yz_yyy_0, \
                                         tez_yz_yyz_0, tez_yz_yzz_0, tez_yz_zzz_0, tez_z_xx_0, tez_z_xx_1, tez_z_xxx_0, \
                                         tez_z_xxx_1, tez_z_xxy_0, tez_z_xxy_1, tez_z_xxz_0, tez_z_xxz_1, tez_z_xy_0, \
                                         tez_z_xy_1, tez_z_xyy_0, tez_z_xyy_1, tez_z_xyz_0, tez_z_xyz_1, tez_z_xz_0, \
                                         tez_z_xz_1, tez_z_xzz_0, tez_z_xzz_1, tez_z_yy_0, tez_z_yy_1, tez_z_yyy_0, \
                                         tez_z_yyy_1, tez_z_yyz_0, tez_z_yyz_1, tez_z_yz_0, tez_z_yz_1, tez_z_yzz_0, \
                                         tez_z_yzz_1, tez_z_zz_0, tez_z_zz_1, tez_z_zzz_0, tez_z_zzz_1, tez_zz_xxx_0, \
                                         tez_zz_xxy_0, tez_zz_xxz_0, tez_zz_xyy_0, tez_zz_xyz_0, tez_zz_xzz_0, tez_zz_yyy_0, \
                                         tez_zz_yyz_0, tez_zz_yzz_0, tez_zz_zzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yz_xzz_0[j] = pa_y[j] * tex_z_xzz_0[j] - pc_y[j] * tex_z_xzz_1[j];

                    tey_yz_xzz_0[j] = pa_y[j] * tey_z_xzz_0[j] - pc_y[j] * tey_z_xzz_1[j] + ta_z_xzz_1[j];

                    tez_yz_xzz_0[j] = pa_y[j] * tez_z_xzz_0[j] - pc_y[j] * tez_z_xzz_1[j];

                    tex_yz_yyy_0[j] = pa_y[j] * tex_z_yyy_0[j] - pc_y[j] * tex_z_yyy_1[j] + 1.5 * fl1_fx * tex_z_yy_0[j] - 1.5 * fl1_fx * tex_z_yy_1[j];

                    tey_yz_yyy_0[j] = pa_y[j] * tey_z_yyy_0[j] - pc_y[j] * tey_z_yyy_1[j] + 1.5 * fl1_fx * tey_z_yy_0[j] - 1.5 * fl1_fx * tey_z_yy_1[j] + ta_z_yyy_1[j];

                    tez_yz_yyy_0[j] = pa_y[j] * tez_z_yyy_0[j] - pc_y[j] * tez_z_yyy_1[j] + 1.5 * fl1_fx * tez_z_yy_0[j] - 1.5 * fl1_fx * tez_z_yy_1[j];

                    tex_yz_yyz_0[j] = pa_y[j] * tex_z_yyz_0[j] - pc_y[j] * tex_z_yyz_1[j] + fl1_fx * tex_z_yz_0[j] - fl1_fx * tex_z_yz_1[j];

                    tey_yz_yyz_0[j] = pa_y[j] * tey_z_yyz_0[j] - pc_y[j] * tey_z_yyz_1[j] + fl1_fx * tey_z_yz_0[j] - fl1_fx * tey_z_yz_1[j] + ta_z_yyz_1[j];

                    tez_yz_yyz_0[j] = pa_y[j] * tez_z_yyz_0[j] - pc_y[j] * tez_z_yyz_1[j] + fl1_fx * tez_z_yz_0[j] - fl1_fx * tez_z_yz_1[j];

                    tex_yz_yzz_0[j] = pa_y[j] * tex_z_yzz_0[j] - pc_y[j] * tex_z_yzz_1[j] + 0.5 * fl1_fx * tex_z_zz_0[j] - 0.5 * fl1_fx * tex_z_zz_1[j];

                    tey_yz_yzz_0[j] = pa_y[j] * tey_z_yzz_0[j] - pc_y[j] * tey_z_yzz_1[j] + 0.5 * fl1_fx * tey_z_zz_0[j] - 0.5 * fl1_fx * tey_z_zz_1[j] + ta_z_yzz_1[j];

                    tez_yz_yzz_0[j] = pa_y[j] * tez_z_yzz_0[j] - pc_y[j] * tez_z_yzz_1[j] + 0.5 * fl1_fx * tez_z_zz_0[j] - 0.5 * fl1_fx * tez_z_zz_1[j];

                    tex_yz_zzz_0[j] = pa_y[j] * tex_z_zzz_0[j] - pc_y[j] * tex_z_zzz_1[j];

                    tey_yz_zzz_0[j] = pa_y[j] * tey_z_zzz_0[j] - pc_y[j] * tey_z_zzz_1[j] + ta_z_zzz_1[j];

                    tez_yz_zzz_0[j] = pa_y[j] * tez_z_zzz_0[j] - pc_y[j] * tez_z_zzz_1[j];

                    tex_zz_xxx_0[j] = pa_z[j] * tex_z_xxx_0[j] - pc_z[j] * tex_z_xxx_1[j] + 0.5 * fl1_fx * tex_0_xxx_0[j] - 0.5 * fl1_fx * tex_0_xxx_1[j];

                    tey_zz_xxx_0[j] = pa_z[j] * tey_z_xxx_0[j] - pc_z[j] * tey_z_xxx_1[j] + 0.5 * fl1_fx * tey_0_xxx_0[j] - 0.5 * fl1_fx * tey_0_xxx_1[j];

                    tez_zz_xxx_0[j] = pa_z[j] * tez_z_xxx_0[j] - pc_z[j] * tez_z_xxx_1[j] + 0.5 * fl1_fx * tez_0_xxx_0[j] - 0.5 * fl1_fx * tez_0_xxx_1[j] + ta_z_xxx_1[j];

                    tex_zz_xxy_0[j] = pa_z[j] * tex_z_xxy_0[j] - pc_z[j] * tex_z_xxy_1[j] + 0.5 * fl1_fx * tex_0_xxy_0[j] - 0.5 * fl1_fx * tex_0_xxy_1[j];

                    tey_zz_xxy_0[j] = pa_z[j] * tey_z_xxy_0[j] - pc_z[j] * tey_z_xxy_1[j] + 0.5 * fl1_fx * tey_0_xxy_0[j] - 0.5 * fl1_fx * tey_0_xxy_1[j];

                    tez_zz_xxy_0[j] = pa_z[j] * tez_z_xxy_0[j] - pc_z[j] * tez_z_xxy_1[j] + 0.5 * fl1_fx * tez_0_xxy_0[j] - 0.5 * fl1_fx * tez_0_xxy_1[j] + ta_z_xxy_1[j];

                    tex_zz_xxz_0[j] = pa_z[j] * tex_z_xxz_0[j] - pc_z[j] * tex_z_xxz_1[j] + 0.5 * fl1_fx * tex_0_xxz_0[j] - 0.5 * fl1_fx * tex_0_xxz_1[j] + 0.5 * fl1_fx * tex_z_xx_0[j] - 0.5 * fl1_fx * tex_z_xx_1[j];

                    tey_zz_xxz_0[j] = pa_z[j] * tey_z_xxz_0[j] - pc_z[j] * tey_z_xxz_1[j] + 0.5 * fl1_fx * tey_0_xxz_0[j] - 0.5 * fl1_fx * tey_0_xxz_1[j] + 0.5 * fl1_fx * tey_z_xx_0[j] - 0.5 * fl1_fx * tey_z_xx_1[j];

                    tez_zz_xxz_0[j] = pa_z[j] * tez_z_xxz_0[j] - pc_z[j] * tez_z_xxz_1[j] + 0.5 * fl1_fx * tez_0_xxz_0[j] - 0.5 * fl1_fx * tez_0_xxz_1[j] + 0.5 * fl1_fx * tez_z_xx_0[j] - 0.5 * fl1_fx * tez_z_xx_1[j] + ta_z_xxz_1[j];

                    tex_zz_xyy_0[j] = pa_z[j] * tex_z_xyy_0[j] - pc_z[j] * tex_z_xyy_1[j] + 0.5 * fl1_fx * tex_0_xyy_0[j] - 0.5 * fl1_fx * tex_0_xyy_1[j];

                    tey_zz_xyy_0[j] = pa_z[j] * tey_z_xyy_0[j] - pc_z[j] * tey_z_xyy_1[j] + 0.5 * fl1_fx * tey_0_xyy_0[j] - 0.5 * fl1_fx * tey_0_xyy_1[j];

                    tez_zz_xyy_0[j] = pa_z[j] * tez_z_xyy_0[j] - pc_z[j] * tez_z_xyy_1[j] + 0.5 * fl1_fx * tez_0_xyy_0[j] - 0.5 * fl1_fx * tez_0_xyy_1[j] + ta_z_xyy_1[j];

                    tex_zz_xyz_0[j] = pa_z[j] * tex_z_xyz_0[j] - pc_z[j] * tex_z_xyz_1[j] + 0.5 * fl1_fx * tex_0_xyz_0[j] - 0.5 * fl1_fx * tex_0_xyz_1[j] + 0.5 * fl1_fx * tex_z_xy_0[j] - 0.5 * fl1_fx * tex_z_xy_1[j];

                    tey_zz_xyz_0[j] = pa_z[j] * tey_z_xyz_0[j] - pc_z[j] * tey_z_xyz_1[j] + 0.5 * fl1_fx * tey_0_xyz_0[j] - 0.5 * fl1_fx * tey_0_xyz_1[j] + 0.5 * fl1_fx * tey_z_xy_0[j] - 0.5 * fl1_fx * tey_z_xy_1[j];

                    tez_zz_xyz_0[j] = pa_z[j] * tez_z_xyz_0[j] - pc_z[j] * tez_z_xyz_1[j] + 0.5 * fl1_fx * tez_0_xyz_0[j] - 0.5 * fl1_fx * tez_0_xyz_1[j] + 0.5 * fl1_fx * tez_z_xy_0[j] - 0.5 * fl1_fx * tez_z_xy_1[j] + ta_z_xyz_1[j];

                    tex_zz_xzz_0[j] = pa_z[j] * tex_z_xzz_0[j] - pc_z[j] * tex_z_xzz_1[j] + 0.5 * fl1_fx * tex_0_xzz_0[j] - 0.5 * fl1_fx * tex_0_xzz_1[j] + fl1_fx * tex_z_xz_0[j] - fl1_fx * tex_z_xz_1[j];

                    tey_zz_xzz_0[j] = pa_z[j] * tey_z_xzz_0[j] - pc_z[j] * tey_z_xzz_1[j] + 0.5 * fl1_fx * tey_0_xzz_0[j] - 0.5 * fl1_fx * tey_0_xzz_1[j] + fl1_fx * tey_z_xz_0[j] - fl1_fx * tey_z_xz_1[j];

                    tez_zz_xzz_0[j] = pa_z[j] * tez_z_xzz_0[j] - pc_z[j] * tez_z_xzz_1[j] + 0.5 * fl1_fx * tez_0_xzz_0[j] - 0.5 * fl1_fx * tez_0_xzz_1[j] + fl1_fx * tez_z_xz_0[j] - fl1_fx * tez_z_xz_1[j] + ta_z_xzz_1[j];

                    tex_zz_yyy_0[j] = pa_z[j] * tex_z_yyy_0[j] - pc_z[j] * tex_z_yyy_1[j] + 0.5 * fl1_fx * tex_0_yyy_0[j] - 0.5 * fl1_fx * tex_0_yyy_1[j];

                    tey_zz_yyy_0[j] = pa_z[j] * tey_z_yyy_0[j] - pc_z[j] * tey_z_yyy_1[j] + 0.5 * fl1_fx * tey_0_yyy_0[j] - 0.5 * fl1_fx * tey_0_yyy_1[j];

                    tez_zz_yyy_0[j] = pa_z[j] * tez_z_yyy_0[j] - pc_z[j] * tez_z_yyy_1[j] + 0.5 * fl1_fx * tez_0_yyy_0[j] - 0.5 * fl1_fx * tez_0_yyy_1[j] + ta_z_yyy_1[j];

                    tex_zz_yyz_0[j] = pa_z[j] * tex_z_yyz_0[j] - pc_z[j] * tex_z_yyz_1[j] + 0.5 * fl1_fx * tex_0_yyz_0[j] - 0.5 * fl1_fx * tex_0_yyz_1[j] + 0.5 * fl1_fx * tex_z_yy_0[j] - 0.5 * fl1_fx * tex_z_yy_1[j];

                    tey_zz_yyz_0[j] = pa_z[j] * tey_z_yyz_0[j] - pc_z[j] * tey_z_yyz_1[j] + 0.5 * fl1_fx * tey_0_yyz_0[j] - 0.5 * fl1_fx * tey_0_yyz_1[j] + 0.5 * fl1_fx * tey_z_yy_0[j] - 0.5 * fl1_fx * tey_z_yy_1[j];

                    tez_zz_yyz_0[j] = pa_z[j] * tez_z_yyz_0[j] - pc_z[j] * tez_z_yyz_1[j] + 0.5 * fl1_fx * tez_0_yyz_0[j] - 0.5 * fl1_fx * tez_0_yyz_1[j] + 0.5 * fl1_fx * tez_z_yy_0[j] - 0.5 * fl1_fx * tez_z_yy_1[j] + ta_z_yyz_1[j];

                    tex_zz_yzz_0[j] = pa_z[j] * tex_z_yzz_0[j] - pc_z[j] * tex_z_yzz_1[j] + 0.5 * fl1_fx * tex_0_yzz_0[j] - 0.5 * fl1_fx * tex_0_yzz_1[j] + fl1_fx * tex_z_yz_0[j] - fl1_fx * tex_z_yz_1[j];

                    tey_zz_yzz_0[j] = pa_z[j] * tey_z_yzz_0[j] - pc_z[j] * tey_z_yzz_1[j] + 0.5 * fl1_fx * tey_0_yzz_0[j] - 0.5 * fl1_fx * tey_0_yzz_1[j] + fl1_fx * tey_z_yz_0[j] - fl1_fx * tey_z_yz_1[j];

                    tez_zz_yzz_0[j] = pa_z[j] * tez_z_yzz_0[j] - pc_z[j] * tez_z_yzz_1[j] + 0.5 * fl1_fx * tez_0_yzz_0[j] - 0.5 * fl1_fx * tez_0_yzz_1[j] + fl1_fx * tez_z_yz_0[j] - fl1_fx * tez_z_yz_1[j] + ta_z_yzz_1[j];

                    tex_zz_zzz_0[j] = pa_z[j] * tex_z_zzz_0[j] - pc_z[j] * tex_z_zzz_1[j] + 0.5 * fl1_fx * tex_0_zzz_0[j] - 0.5 * fl1_fx * tex_0_zzz_1[j] + 1.5 * fl1_fx * tex_z_zz_0[j] - 1.5 * fl1_fx * tex_z_zz_1[j];

                    tey_zz_zzz_0[j] = pa_z[j] * tey_z_zzz_0[j] - pc_z[j] * tey_z_zzz_1[j] + 0.5 * fl1_fx * tey_0_zzz_0[j] - 0.5 * fl1_fx * tey_0_zzz_1[j] + 1.5 * fl1_fx * tey_z_zz_0[j] - 1.5 * fl1_fx * tey_z_zz_1[j];

                    tez_zz_zzz_0[j] = pa_z[j] * tez_z_zzz_0[j] - pc_z[j] * tez_z_zzz_1[j] + 0.5 * fl1_fx * tez_0_zzz_0[j] - 0.5 * fl1_fx * tez_0_zzz_1[j] + 1.5 * fl1_fx * tez_z_zz_0[j] - 1.5 * fl1_fx * tez_z_zz_1[j] + ta_z_zzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFD(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pcDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        efieldrecfunc::compElectricFieldForFD_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   pcDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        efieldrecfunc::compElectricFieldForFD_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    pcDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        efieldrecfunc::compElectricFieldForFD_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     pcDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        efieldrecfunc::compElectricFieldForFD_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compElectricFieldForFD_0_45(      CMemBlock2D<double>& primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
    {
        // Batch of Integrals (0,45)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tex_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx); 

                auto tey_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx); 

                auto tez_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx); 

                auto tex_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 1); 

                auto tey_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 1); 

                auto tez_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 1); 

                auto tex_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 2); 

                auto tey_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 2); 

                auto tez_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 2); 

                auto tex_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 3); 

                auto tey_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 3); 

                auto tez_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 3); 

                auto tex_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 4); 

                auto tey_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 4); 

                auto tez_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 4); 

                auto tex_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 5); 

                auto tey_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 5); 

                auto tez_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 5); 

                auto tex_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 6); 

                auto tey_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 6); 

                auto tez_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 6); 

                auto tex_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 7); 

                auto tey_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 7); 

                auto tez_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 7); 

                auto tex_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 8); 

                auto tey_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 8); 

                auto tez_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 8); 

                auto tex_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 9); 

                auto tey_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 9); 

                auto tez_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 9); 

                auto tex_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 10); 

                auto tey_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 10); 

                auto tez_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 10); 

                auto tex_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 11); 

                auto tey_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 11); 

                auto tez_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 11); 

                auto tex_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 12); 

                auto tey_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 12); 

                auto tez_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 12); 

                auto tex_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 13); 

                auto tey_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 13); 

                auto tez_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 13); 

                auto tex_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 14); 

                auto tey_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 14); 

                auto tez_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 14); 

                auto tex_xx_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx); 

                auto tey_xx_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx); 

                auto tez_xx_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx); 

                auto tex_xx_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 1); 

                auto tey_xx_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 1); 

                auto tez_xx_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 1); 

                auto tex_xx_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 2); 

                auto tey_xx_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 2); 

                auto tez_xx_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 2); 

                auto tex_xx_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 3); 

                auto tey_xx_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 3); 

                auto tez_xx_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 3); 

                auto tex_xx_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 4); 

                auto tey_xx_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 4); 

                auto tez_xx_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 4); 

                auto tex_xx_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 5); 

                auto tey_xx_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 5); 

                auto tez_xx_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 5); 

                auto tex_xy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 6); 

                auto tey_xy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 6); 

                auto tez_xy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 6); 

                auto tex_xy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 7); 

                auto tey_xy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 7); 

                auto tez_xy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 7); 

                auto tex_xy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 8); 

                auto tey_xy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 8); 

                auto tez_xy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 8); 

                auto tex_xy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 9); 

                auto tey_xy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 9); 

                auto tez_xy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 9); 

                auto tex_xy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 10); 

                auto tey_xy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 10); 

                auto tez_xy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 10); 

                auto tex_xy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 11); 

                auto tey_xy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 11); 

                auto tez_xy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 11); 

                auto tex_xz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 12); 

                auto tey_xz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 12); 

                auto tez_xz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 12); 

                auto tex_xz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 13); 

                auto tey_xz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 13); 

                auto tez_xz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 13); 

                auto tex_xz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 14); 

                auto tey_xz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 14); 

                auto tez_xz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 14); 

                auto tex_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx); 

                auto tey_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx); 

                auto tez_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx); 

                auto tex_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 1); 

                auto tey_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 1); 

                auto tez_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 1); 

                auto tex_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 2); 

                auto tey_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 2); 

                auto tez_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 2); 

                auto tex_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 3); 

                auto tey_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 3); 

                auto tez_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 3); 

                auto tex_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 4); 

                auto tey_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 4); 

                auto tez_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 4); 

                auto tex_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 5); 

                auto tey_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 5); 

                auto tez_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 5); 

                auto tex_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 6); 

                auto tey_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 7); 

                auto tey_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 8); 

                auto tey_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 9); 

                auto tey_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 10); 

                auto tey_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 11); 

                auto tey_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx); 

                auto tey_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx); 

                auto tez_x_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx); 

                auto tex_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 1); 

                auto tey_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 1); 

                auto tez_x_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 1); 

                auto tex_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 2); 

                auto tey_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 2); 

                auto tez_x_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 2); 

                auto tex_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 3); 

                auto tey_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 3); 

                auto tez_x_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 3); 

                auto tex_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 4); 

                auto tey_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 4); 

                auto tez_x_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 4); 

                auto tex_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 5); 

                auto tey_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 5); 

                auto tez_x_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 5); 

                auto tex_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 6); 

                auto tey_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 7); 

                auto tey_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 8); 

                auto tey_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 9); 

                auto tey_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 10); 

                auto tey_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 11); 

                auto tey_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx); 

                auto tey_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx); 

                auto tez_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx); 

                auto tex_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 1); 

                auto tey_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 1); 

                auto tez_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 1); 

                auto tex_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 2); 

                auto tey_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 2); 

                auto tez_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 2); 

                auto tex_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 3); 

                auto tey_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 3); 

                auto tez_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 3); 

                auto tex_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 4); 

                auto tey_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 4); 

                auto tez_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 4); 

                auto tex_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 5); 

                auto tey_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 5); 

                auto tez_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 5); 

                auto tex_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 6); 

                auto tey_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 6); 

                auto tez_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 6); 

                auto tex_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 7); 

                auto tey_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 7); 

                auto tez_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 7); 

                auto tex_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 8); 

                auto tey_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 8); 

                auto tez_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 8); 

                auto tex_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx); 

                auto tey_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx); 

                auto tez_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx); 

                auto tex_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 1); 

                auto tey_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 1); 

                auto tez_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 1); 

                auto tex_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 2); 

                auto tey_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 2); 

                auto tez_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 2); 

                auto tex_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 3); 

                auto tey_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 3); 

                auto tez_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 3); 

                auto tex_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 4); 

                auto tey_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 4); 

                auto tez_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 4); 

                auto tex_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 5); 

                auto tey_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 5); 

                auto tez_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 5); 

                auto tex_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 6); 

                auto tey_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 6); 

                auto tez_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 6); 

                auto tex_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 7); 

                auto tey_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 7); 

                auto tez_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 7); 

                auto tex_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 8); 

                auto tey_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 8); 

                auto tez_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 8); 

                auto ta_xx_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx); 

                auto ta_xx_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 1); 

                auto ta_xx_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 2); 

                auto ta_xx_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 3); 

                auto ta_xx_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 4); 

                auto ta_xx_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 5); 

                auto ta_xy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 6); 

                auto ta_xy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 7); 

                auto ta_xy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 8); 

                auto ta_xy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 9); 

                auto ta_xy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 10); 

                auto ta_xy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 11); 

                auto ta_xz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 12); 

                auto ta_xz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 13); 

                auto ta_xz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 14); 

                // set up pointers to integrals

                auto tex_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx); 

                auto tey_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx); 

                auto tez_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx); 

                auto tex_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 1); 

                auto tey_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 1); 

                auto tez_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 1); 

                auto tex_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 2); 

                auto tey_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 2); 

                auto tez_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 2); 

                auto tex_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 3); 

                auto tey_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 3); 

                auto tez_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 3); 

                auto tex_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 4); 

                auto tey_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 4); 

                auto tez_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 4); 

                auto tex_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 5); 

                auto tey_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 5); 

                auto tez_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 5); 

                auto tex_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 6); 

                auto tey_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 6); 

                auto tez_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 6); 

                auto tex_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 7); 

                auto tey_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 7); 

                auto tez_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 7); 

                auto tex_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 8); 

                auto tey_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 8); 

                auto tez_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 8); 

                auto tex_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 9); 

                auto tey_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 9); 

                auto tez_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 9); 

                auto tex_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 10); 

                auto tey_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 10); 

                auto tez_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 10); 

                auto tex_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 11); 

                auto tey_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 11); 

                auto tez_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 11); 

                auto tex_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 12); 

                auto tey_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 12); 

                auto tez_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 12); 

                auto tex_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 13); 

                auto tey_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 13); 

                auto tez_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 13); 

                auto tex_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 14); 

                auto tey_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 14); 

                auto tez_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 14); 

                // Batch of Integrals (0,45)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xx_xx_1, ta_xx_xy_1, ta_xx_xz_1, ta_xx_yy_1, ta_xx_yz_1, \
                                         ta_xx_zz_1, ta_xy_xx_1, ta_xy_xy_1, ta_xy_xz_1, ta_xy_yy_1, ta_xy_yz_1, ta_xy_zz_1, \
                                         ta_xz_xx_1, ta_xz_xy_1, ta_xz_xz_1, tex_x_xx_0, tex_x_xx_1, tex_x_xy_0, tex_x_xy_1, \
                                         tex_x_xz_0, tex_x_xz_1, tex_x_yy_0, tex_x_yy_1, tex_x_yz_0, tex_x_yz_1, tex_x_zz_0, \
                                         tex_x_zz_1, tex_xx_x_0, tex_xx_x_1, tex_xx_xx_0, tex_xx_xx_1, tex_xx_xy_0, \
                                         tex_xx_xy_1, tex_xx_xz_0, tex_xx_xz_1, tex_xx_y_0, tex_xx_y_1, tex_xx_yy_0, \
                                         tex_xx_yy_1, tex_xx_yz_0, tex_xx_yz_1, tex_xx_z_0, tex_xx_z_1, tex_xx_zz_0, \
                                         tex_xx_zz_1, tex_xxx_xx_0, tex_xxx_xy_0, tex_xxx_xz_0, tex_xxx_yy_0, tex_xxx_yz_0, \
                                         tex_xxx_zz_0, tex_xxy_xx_0, tex_xxy_xy_0, tex_xxy_xz_0, tex_xxy_yy_0, tex_xxy_yz_0, \
                                         tex_xxy_zz_0, tex_xxz_xx_0, tex_xxz_xy_0, tex_xxz_xz_0, tex_xy_x_0, tex_xy_x_1, \
                                         tex_xy_xx_0, tex_xy_xx_1, tex_xy_xy_0, tex_xy_xy_1, tex_xy_xz_0, tex_xy_xz_1, \
                                         tex_xy_y_0, tex_xy_y_1, tex_xy_yy_0, tex_xy_yy_1, tex_xy_yz_0, tex_xy_yz_1, \
                                         tex_xy_z_0, tex_xy_z_1, tex_xy_zz_0, tex_xy_zz_1, tex_xz_x_0, tex_xz_x_1, \
                                         tex_xz_xx_0, tex_xz_xx_1, tex_xz_xy_0, tex_xz_xy_1, tex_xz_xz_0, tex_xz_xz_1, \
                                         tex_xz_y_0, tex_xz_y_1, tex_xz_z_0, tex_xz_z_1, tex_y_xx_0, tex_y_xx_1, tex_y_xy_0, \
                                         tex_y_xy_1, tex_y_xz_0, tex_y_xz_1, tex_y_yy_0, tex_y_yy_1, tex_y_yz_0, tex_y_yz_1, \
                                         tex_y_zz_0, tex_y_zz_1, tex_z_xx_0, tex_z_xx_1, tex_z_xy_0, tex_z_xy_1, tex_z_xz_0, \
                                         tex_z_xz_1, tey_x_xx_0, tey_x_xx_1, tey_x_xy_0, tey_x_xy_1, tey_x_xz_0, tey_x_xz_1, \
                                         tey_x_yy_0, tey_x_yy_1, tey_x_yz_0, tey_x_yz_1, tey_x_zz_0, tey_x_zz_1, tey_xx_x_0, \
                                         tey_xx_x_1, tey_xx_xx_0, tey_xx_xx_1, tey_xx_xy_0, tey_xx_xy_1, tey_xx_xz_0, \
                                         tey_xx_xz_1, tey_xx_y_0, tey_xx_y_1, tey_xx_yy_0, tey_xx_yy_1, tey_xx_yz_0, \
                                         tey_xx_yz_1, tey_xx_z_0, tey_xx_z_1, tey_xx_zz_0, tey_xx_zz_1, tey_xxx_xx_0, \
                                         tey_xxx_xy_0, tey_xxx_xz_0, tey_xxx_yy_0, tey_xxx_yz_0, tey_xxx_zz_0, tey_xxy_xx_0, \
                                         tey_xxy_xy_0, tey_xxy_xz_0, tey_xxy_yy_0, tey_xxy_yz_0, tey_xxy_zz_0, tey_xxz_xx_0, \
                                         tey_xxz_xy_0, tey_xxz_xz_0, tey_xy_x_0, tey_xy_x_1, tey_xy_xx_0, tey_xy_xx_1, \
                                         tey_xy_xy_0, tey_xy_xy_1, tey_xy_xz_0, tey_xy_xz_1, tey_xy_y_0, tey_xy_y_1, \
                                         tey_xy_yy_0, tey_xy_yy_1, tey_xy_yz_0, tey_xy_yz_1, tey_xy_z_0, tey_xy_z_1, \
                                         tey_xy_zz_0, tey_xy_zz_1, tey_xz_x_0, tey_xz_x_1, tey_xz_xx_0, tey_xz_xx_1, \
                                         tey_xz_xy_0, tey_xz_xy_1, tey_xz_xz_0, tey_xz_xz_1, tey_xz_y_0, tey_xz_y_1, \
                                         tey_xz_z_0, tey_xz_z_1, tey_y_xx_0, tey_y_xx_1, tey_y_xy_0, tey_y_xy_1, tey_y_xz_0, \
                                         tey_y_xz_1, tey_y_yy_0, tey_y_yy_1, tey_y_yz_0, tey_y_yz_1, tey_y_zz_0, tey_y_zz_1, \
                                         tey_z_xx_0, tey_z_xx_1, tey_z_xy_0, tey_z_xy_1, tey_z_xz_0, tey_z_xz_1, tez_x_xx_0, \
                                         tez_x_xx_1, tez_x_xy_0, tez_x_xy_1, tez_x_xz_0, tez_x_xz_1, tez_x_yy_0, tez_x_yy_1, \
                                         tez_x_yz_0, tez_x_yz_1, tez_x_zz_0, tez_x_zz_1, tez_xx_x_0, tez_xx_x_1, \
                                         tez_xx_xx_0, tez_xx_xx_1, tez_xx_xy_0, tez_xx_xy_1, tez_xx_xz_0, tez_xx_xz_1, \
                                         tez_xx_y_0, tez_xx_y_1, tez_xx_yy_0, tez_xx_yy_1, tez_xx_yz_0, tez_xx_yz_1, \
                                         tez_xx_z_0, tez_xx_z_1, tez_xx_zz_0, tez_xx_zz_1, tez_xxx_xx_0, tez_xxx_xy_0, \
                                         tez_xxx_xz_0, tez_xxx_yy_0, tez_xxx_yz_0, tez_xxx_zz_0, tez_xxy_xx_0, tez_xxy_xy_0, \
                                         tez_xxy_xz_0, tez_xxy_yy_0, tez_xxy_yz_0, tez_xxy_zz_0, tez_xxz_xx_0, tez_xxz_xy_0, \
                                         tez_xxz_xz_0, tez_xy_x_0, tez_xy_x_1, tez_xy_xx_0, tez_xy_xx_1, tez_xy_xy_0, \
                                         tez_xy_xy_1, tez_xy_xz_0, tez_xy_xz_1, tez_xy_y_0, tez_xy_y_1, tez_xy_yy_0, \
                                         tez_xy_yy_1, tez_xy_yz_0, tez_xy_yz_1, tez_xy_z_0, tez_xy_z_1, tez_xy_zz_0, \
                                         tez_xy_zz_1, tez_xz_x_0, tez_xz_x_1, tez_xz_xx_0, tez_xz_xx_1, tez_xz_xy_0, \
                                         tez_xz_xy_1, tez_xz_xz_0, tez_xz_xz_1, tez_xz_y_0, tez_xz_y_1, tez_xz_z_0, \
                                         tez_xz_z_1, tez_y_xx_0, tez_y_xx_1, tez_y_xy_0, tez_y_xy_1, tez_y_xz_0, tez_y_xz_1, \
                                         tez_y_yy_0, tez_y_yy_1, tez_y_yz_0, tez_y_yz_1, tez_y_zz_0, tez_y_zz_1, tez_z_xx_0, \
                                         tez_z_xx_1, tez_z_xy_0, tez_z_xy_1, tez_z_xz_0, tez_z_xz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xxx_xx_0[j] = pa_x[j] * tex_xx_xx_0[j] - pc_x[j] * tex_xx_xx_1[j] + fl1_fx * tex_x_xx_0[j] - fl1_fx * tex_x_xx_1[j] + fl1_fx * tex_xx_x_0[j] - fl1_fx * tex_xx_x_1[j] + ta_xx_xx_1[j];

                    tey_xxx_xx_0[j] = pa_x[j] * tey_xx_xx_0[j] - pc_x[j] * tey_xx_xx_1[j] + fl1_fx * tey_x_xx_0[j] - fl1_fx * tey_x_xx_1[j] + fl1_fx * tey_xx_x_0[j] - fl1_fx * tey_xx_x_1[j];

                    tez_xxx_xx_0[j] = pa_x[j] * tez_xx_xx_0[j] - pc_x[j] * tez_xx_xx_1[j] + fl1_fx * tez_x_xx_0[j] - fl1_fx * tez_x_xx_1[j] + fl1_fx * tez_xx_x_0[j] - fl1_fx * tez_xx_x_1[j];

                    tex_xxx_xy_0[j] = pa_x[j] * tex_xx_xy_0[j] - pc_x[j] * tex_xx_xy_1[j] + fl1_fx * tex_x_xy_0[j] - fl1_fx * tex_x_xy_1[j] + 0.5 * fl1_fx * tex_xx_y_0[j] - 0.5 * fl1_fx * tex_xx_y_1[j] + ta_xx_xy_1[j];

                    tey_xxx_xy_0[j] = pa_x[j] * tey_xx_xy_0[j] - pc_x[j] * tey_xx_xy_1[j] + fl1_fx * tey_x_xy_0[j] - fl1_fx * tey_x_xy_1[j] + 0.5 * fl1_fx * tey_xx_y_0[j] - 0.5 * fl1_fx * tey_xx_y_1[j];

                    tez_xxx_xy_0[j] = pa_x[j] * tez_xx_xy_0[j] - pc_x[j] * tez_xx_xy_1[j] + fl1_fx * tez_x_xy_0[j] - fl1_fx * tez_x_xy_1[j] + 0.5 * fl1_fx * tez_xx_y_0[j] - 0.5 * fl1_fx * tez_xx_y_1[j];

                    tex_xxx_xz_0[j] = pa_x[j] * tex_xx_xz_0[j] - pc_x[j] * tex_xx_xz_1[j] + fl1_fx * tex_x_xz_0[j] - fl1_fx * tex_x_xz_1[j] + 0.5 * fl1_fx * tex_xx_z_0[j] - 0.5 * fl1_fx * tex_xx_z_1[j] + ta_xx_xz_1[j];

                    tey_xxx_xz_0[j] = pa_x[j] * tey_xx_xz_0[j] - pc_x[j] * tey_xx_xz_1[j] + fl1_fx * tey_x_xz_0[j] - fl1_fx * tey_x_xz_1[j] + 0.5 * fl1_fx * tey_xx_z_0[j] - 0.5 * fl1_fx * tey_xx_z_1[j];

                    tez_xxx_xz_0[j] = pa_x[j] * tez_xx_xz_0[j] - pc_x[j] * tez_xx_xz_1[j] + fl1_fx * tez_x_xz_0[j] - fl1_fx * tez_x_xz_1[j] + 0.5 * fl1_fx * tez_xx_z_0[j] - 0.5 * fl1_fx * tez_xx_z_1[j];

                    tex_xxx_yy_0[j] = pa_x[j] * tex_xx_yy_0[j] - pc_x[j] * tex_xx_yy_1[j] + fl1_fx * tex_x_yy_0[j] - fl1_fx * tex_x_yy_1[j] + ta_xx_yy_1[j];

                    tey_xxx_yy_0[j] = pa_x[j] * tey_xx_yy_0[j] - pc_x[j] * tey_xx_yy_1[j] + fl1_fx * tey_x_yy_0[j] - fl1_fx * tey_x_yy_1[j];

                    tez_xxx_yy_0[j] = pa_x[j] * tez_xx_yy_0[j] - pc_x[j] * tez_xx_yy_1[j] + fl1_fx * tez_x_yy_0[j] - fl1_fx * tez_x_yy_1[j];

                    tex_xxx_yz_0[j] = pa_x[j] * tex_xx_yz_0[j] - pc_x[j] * tex_xx_yz_1[j] + fl1_fx * tex_x_yz_0[j] - fl1_fx * tex_x_yz_1[j] + ta_xx_yz_1[j];

                    tey_xxx_yz_0[j] = pa_x[j] * tey_xx_yz_0[j] - pc_x[j] * tey_xx_yz_1[j] + fl1_fx * tey_x_yz_0[j] - fl1_fx * tey_x_yz_1[j];

                    tez_xxx_yz_0[j] = pa_x[j] * tez_xx_yz_0[j] - pc_x[j] * tez_xx_yz_1[j] + fl1_fx * tez_x_yz_0[j] - fl1_fx * tez_x_yz_1[j];

                    tex_xxx_zz_0[j] = pa_x[j] * tex_xx_zz_0[j] - pc_x[j] * tex_xx_zz_1[j] + fl1_fx * tex_x_zz_0[j] - fl1_fx * tex_x_zz_1[j] + ta_xx_zz_1[j];

                    tey_xxx_zz_0[j] = pa_x[j] * tey_xx_zz_0[j] - pc_x[j] * tey_xx_zz_1[j] + fl1_fx * tey_x_zz_0[j] - fl1_fx * tey_x_zz_1[j];

                    tez_xxx_zz_0[j] = pa_x[j] * tez_xx_zz_0[j] - pc_x[j] * tez_xx_zz_1[j] + fl1_fx * tez_x_zz_0[j] - fl1_fx * tez_x_zz_1[j];

                    tex_xxy_xx_0[j] = pa_x[j] * tex_xy_xx_0[j] - pc_x[j] * tex_xy_xx_1[j] + 0.5 * fl1_fx * tex_y_xx_0[j] - 0.5 * fl1_fx * tex_y_xx_1[j] + fl1_fx * tex_xy_x_0[j] - fl1_fx * tex_xy_x_1[j] + ta_xy_xx_1[j];

                    tey_xxy_xx_0[j] = pa_x[j] * tey_xy_xx_0[j] - pc_x[j] * tey_xy_xx_1[j] + 0.5 * fl1_fx * tey_y_xx_0[j] - 0.5 * fl1_fx * tey_y_xx_1[j] + fl1_fx * tey_xy_x_0[j] - fl1_fx * tey_xy_x_1[j];

                    tez_xxy_xx_0[j] = pa_x[j] * tez_xy_xx_0[j] - pc_x[j] * tez_xy_xx_1[j] + 0.5 * fl1_fx * tez_y_xx_0[j] - 0.5 * fl1_fx * tez_y_xx_1[j] + fl1_fx * tez_xy_x_0[j] - fl1_fx * tez_xy_x_1[j];

                    tex_xxy_xy_0[j] = pa_x[j] * tex_xy_xy_0[j] - pc_x[j] * tex_xy_xy_1[j] + 0.5 * fl1_fx * tex_y_xy_0[j] - 0.5 * fl1_fx * tex_y_xy_1[j] + 0.5 * fl1_fx * tex_xy_y_0[j] - 0.5 * fl1_fx * tex_xy_y_1[j] + ta_xy_xy_1[j];

                    tey_xxy_xy_0[j] = pa_x[j] * tey_xy_xy_0[j] - pc_x[j] * tey_xy_xy_1[j] + 0.5 * fl1_fx * tey_y_xy_0[j] - 0.5 * fl1_fx * tey_y_xy_1[j] + 0.5 * fl1_fx * tey_xy_y_0[j] - 0.5 * fl1_fx * tey_xy_y_1[j];

                    tez_xxy_xy_0[j] = pa_x[j] * tez_xy_xy_0[j] - pc_x[j] * tez_xy_xy_1[j] + 0.5 * fl1_fx * tez_y_xy_0[j] - 0.5 * fl1_fx * tez_y_xy_1[j] + 0.5 * fl1_fx * tez_xy_y_0[j] - 0.5 * fl1_fx * tez_xy_y_1[j];

                    tex_xxy_xz_0[j] = pa_x[j] * tex_xy_xz_0[j] - pc_x[j] * tex_xy_xz_1[j] + 0.5 * fl1_fx * tex_y_xz_0[j] - 0.5 * fl1_fx * tex_y_xz_1[j] + 0.5 * fl1_fx * tex_xy_z_0[j] - 0.5 * fl1_fx * tex_xy_z_1[j] + ta_xy_xz_1[j];

                    tey_xxy_xz_0[j] = pa_x[j] * tey_xy_xz_0[j] - pc_x[j] * tey_xy_xz_1[j] + 0.5 * fl1_fx * tey_y_xz_0[j] - 0.5 * fl1_fx * tey_y_xz_1[j] + 0.5 * fl1_fx * tey_xy_z_0[j] - 0.5 * fl1_fx * tey_xy_z_1[j];

                    tez_xxy_xz_0[j] = pa_x[j] * tez_xy_xz_0[j] - pc_x[j] * tez_xy_xz_1[j] + 0.5 * fl1_fx * tez_y_xz_0[j] - 0.5 * fl1_fx * tez_y_xz_1[j] + 0.5 * fl1_fx * tez_xy_z_0[j] - 0.5 * fl1_fx * tez_xy_z_1[j];

                    tex_xxy_yy_0[j] = pa_x[j] * tex_xy_yy_0[j] - pc_x[j] * tex_xy_yy_1[j] + 0.5 * fl1_fx * tex_y_yy_0[j] - 0.5 * fl1_fx * tex_y_yy_1[j] + ta_xy_yy_1[j];

                    tey_xxy_yy_0[j] = pa_x[j] * tey_xy_yy_0[j] - pc_x[j] * tey_xy_yy_1[j] + 0.5 * fl1_fx * tey_y_yy_0[j] - 0.5 * fl1_fx * tey_y_yy_1[j];

                    tez_xxy_yy_0[j] = pa_x[j] * tez_xy_yy_0[j] - pc_x[j] * tez_xy_yy_1[j] + 0.5 * fl1_fx * tez_y_yy_0[j] - 0.5 * fl1_fx * tez_y_yy_1[j];

                    tex_xxy_yz_0[j] = pa_x[j] * tex_xy_yz_0[j] - pc_x[j] * tex_xy_yz_1[j] + 0.5 * fl1_fx * tex_y_yz_0[j] - 0.5 * fl1_fx * tex_y_yz_1[j] + ta_xy_yz_1[j];

                    tey_xxy_yz_0[j] = pa_x[j] * tey_xy_yz_0[j] - pc_x[j] * tey_xy_yz_1[j] + 0.5 * fl1_fx * tey_y_yz_0[j] - 0.5 * fl1_fx * tey_y_yz_1[j];

                    tez_xxy_yz_0[j] = pa_x[j] * tez_xy_yz_0[j] - pc_x[j] * tez_xy_yz_1[j] + 0.5 * fl1_fx * tez_y_yz_0[j] - 0.5 * fl1_fx * tez_y_yz_1[j];

                    tex_xxy_zz_0[j] = pa_x[j] * tex_xy_zz_0[j] - pc_x[j] * tex_xy_zz_1[j] + 0.5 * fl1_fx * tex_y_zz_0[j] - 0.5 * fl1_fx * tex_y_zz_1[j] + ta_xy_zz_1[j];

                    tey_xxy_zz_0[j] = pa_x[j] * tey_xy_zz_0[j] - pc_x[j] * tey_xy_zz_1[j] + 0.5 * fl1_fx * tey_y_zz_0[j] - 0.5 * fl1_fx * tey_y_zz_1[j];

                    tez_xxy_zz_0[j] = pa_x[j] * tez_xy_zz_0[j] - pc_x[j] * tez_xy_zz_1[j] + 0.5 * fl1_fx * tez_y_zz_0[j] - 0.5 * fl1_fx * tez_y_zz_1[j];

                    tex_xxz_xx_0[j] = pa_x[j] * tex_xz_xx_0[j] - pc_x[j] * tex_xz_xx_1[j] + 0.5 * fl1_fx * tex_z_xx_0[j] - 0.5 * fl1_fx * tex_z_xx_1[j] + fl1_fx * tex_xz_x_0[j] - fl1_fx * tex_xz_x_1[j] + ta_xz_xx_1[j];

                    tey_xxz_xx_0[j] = pa_x[j] * tey_xz_xx_0[j] - pc_x[j] * tey_xz_xx_1[j] + 0.5 * fl1_fx * tey_z_xx_0[j] - 0.5 * fl1_fx * tey_z_xx_1[j] + fl1_fx * tey_xz_x_0[j] - fl1_fx * tey_xz_x_1[j];

                    tez_xxz_xx_0[j] = pa_x[j] * tez_xz_xx_0[j] - pc_x[j] * tez_xz_xx_1[j] + 0.5 * fl1_fx * tez_z_xx_0[j] - 0.5 * fl1_fx * tez_z_xx_1[j] + fl1_fx * tez_xz_x_0[j] - fl1_fx * tez_xz_x_1[j];

                    tex_xxz_xy_0[j] = pa_x[j] * tex_xz_xy_0[j] - pc_x[j] * tex_xz_xy_1[j] + 0.5 * fl1_fx * tex_z_xy_0[j] - 0.5 * fl1_fx * tex_z_xy_1[j] + 0.5 * fl1_fx * tex_xz_y_0[j] - 0.5 * fl1_fx * tex_xz_y_1[j] + ta_xz_xy_1[j];

                    tey_xxz_xy_0[j] = pa_x[j] * tey_xz_xy_0[j] - pc_x[j] * tey_xz_xy_1[j] + 0.5 * fl1_fx * tey_z_xy_0[j] - 0.5 * fl1_fx * tey_z_xy_1[j] + 0.5 * fl1_fx * tey_xz_y_0[j] - 0.5 * fl1_fx * tey_xz_y_1[j];

                    tez_xxz_xy_0[j] = pa_x[j] * tez_xz_xy_0[j] - pc_x[j] * tez_xz_xy_1[j] + 0.5 * fl1_fx * tez_z_xy_0[j] - 0.5 * fl1_fx * tez_z_xy_1[j] + 0.5 * fl1_fx * tez_xz_y_0[j] - 0.5 * fl1_fx * tez_xz_y_1[j];

                    tex_xxz_xz_0[j] = pa_x[j] * tex_xz_xz_0[j] - pc_x[j] * tex_xz_xz_1[j] + 0.5 * fl1_fx * tex_z_xz_0[j] - 0.5 * fl1_fx * tex_z_xz_1[j] + 0.5 * fl1_fx * tex_xz_z_0[j] - 0.5 * fl1_fx * tex_xz_z_1[j] + ta_xz_xz_1[j];

                    tey_xxz_xz_0[j] = pa_x[j] * tey_xz_xz_0[j] - pc_x[j] * tey_xz_xz_1[j] + 0.5 * fl1_fx * tey_z_xz_0[j] - 0.5 * fl1_fx * tey_z_xz_1[j] + 0.5 * fl1_fx * tey_xz_z_0[j] - 0.5 * fl1_fx * tey_xz_z_1[j];

                    tez_xxz_xz_0[j] = pa_x[j] * tez_xz_xz_0[j] - pc_x[j] * tez_xz_xz_1[j] + 0.5 * fl1_fx * tez_z_xz_0[j] - 0.5 * fl1_fx * tez_z_xz_1[j] + 0.5 * fl1_fx * tez_xz_z_0[j] - 0.5 * fl1_fx * tez_xz_z_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFD_45_90(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (45,90)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
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

                auto tex_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 15); 

                auto tey_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 15); 

                auto tez_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 15); 

                auto tex_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 16); 

                auto tey_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 16); 

                auto tez_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 16); 

                auto tex_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 17); 

                auto tey_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 17); 

                auto tez_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 17); 

                auto tex_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 18); 

                auto tey_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 19); 

                auto tey_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 20); 

                auto tey_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 21); 

                auto tey_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 22); 

                auto tey_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 23); 

                auto tey_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 24); 

                auto tey_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 25); 

                auto tey_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 26); 

                auto tey_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 26); 

                auto tex_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 27); 

                auto tey_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 28); 

                auto tey_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 29); 

                auto tey_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 29); 

                auto tex_xz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 15); 

                auto tey_xz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 15); 

                auto tez_xz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 15); 

                auto tex_xz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 16); 

                auto tey_xz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 16); 

                auto tez_xz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 16); 

                auto tex_xz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 17); 

                auto tey_xz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 17); 

                auto tez_xz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 17); 

                auto tex_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 18); 

                auto tey_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 19); 

                auto tey_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 20); 

                auto tey_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 21); 

                auto tey_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 22); 

                auto tey_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 23); 

                auto tey_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 24); 

                auto tey_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 25); 

                auto tey_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 26); 

                auto tey_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 26); 

                auto tex_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 27); 

                auto tey_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 28); 

                auto tey_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 29); 

                auto tey_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 29); 

                auto tex_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 15); 

                auto tey_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 16); 

                auto tey_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 17); 

                auto tey_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 15); 

                auto tey_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 16); 

                auto tey_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 17); 

                auto tey_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 17); 

                auto tex_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 9); 

                auto tey_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 10); 

                auto tey_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 11); 

                auto tey_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 12); 

                auto tey_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 13); 

                auto tey_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 14); 

                auto tey_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 9); 

                auto tey_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 10); 

                auto tey_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 10); 

                auto tex_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 11); 

                auto tey_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 12); 

                auto tey_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 13); 

                auto tey_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 14); 

                auto tey_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 14); 

                auto ta_xz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 15); 

                auto ta_xz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 16); 

                auto ta_xz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 17); 

                auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18); 

                auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19); 

                auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20); 

                auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21); 

                auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22); 

                auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23); 

                auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24); 

                auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25); 

                auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26); 

                auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27); 

                auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28); 

                auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29); 

                // set up pointers to integrals

                auto tex_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 15); 

                auto tey_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 15); 

                auto tez_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 15); 

                auto tex_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 16); 

                auto tey_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 16); 

                auto tez_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 16); 

                auto tex_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 17); 

                auto tey_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 17); 

                auto tez_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 17); 

                auto tex_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 18); 

                auto tey_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 18); 

                auto tez_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 18); 

                auto tex_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 19); 

                auto tey_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 19); 

                auto tez_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 19); 

                auto tex_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 20); 

                auto tey_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 20); 

                auto tez_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 20); 

                auto tex_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 21); 

                auto tey_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 21); 

                auto tez_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 21); 

                auto tex_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 22); 

                auto tey_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 22); 

                auto tez_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 22); 

                auto tex_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 23); 

                auto tey_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 23); 

                auto tez_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 23); 

                auto tex_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 24); 

                auto tey_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 24); 

                auto tez_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 24); 

                auto tex_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 25); 

                auto tey_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 25); 

                auto tez_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 25); 

                auto tex_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 26); 

                auto tey_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 26); 

                auto tez_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 26); 

                auto tex_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 27); 

                auto tey_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 27); 

                auto tez_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 27); 

                auto tex_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 28); 

                auto tey_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 28); 

                auto tez_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 28); 

                auto tex_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 29); 

                auto tey_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 29); 

                auto tez_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 29); 

                // Batch of Integrals (45,90)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xz_yy_1, ta_xz_yz_1, ta_xz_zz_1, ta_yy_xx_1, ta_yy_xy_1, \
                                         ta_yy_xz_1, ta_yy_yy_1, ta_yy_yz_1, ta_yy_zz_1, ta_yz_xx_1, ta_yz_xy_1, ta_yz_xz_1, \
                                         ta_yz_yy_1, ta_yz_yz_1, ta_yz_zz_1, tex_xxz_yy_0, tex_xxz_yz_0, tex_xxz_zz_0, \
                                         tex_xyy_xx_0, tex_xyy_xy_0, tex_xyy_xz_0, tex_xyy_yy_0, tex_xyy_yz_0, tex_xyy_zz_0, \
                                         tex_xyz_xx_0, tex_xyz_xy_0, tex_xyz_xz_0, tex_xyz_yy_0, tex_xyz_yz_0, tex_xyz_zz_0, \
                                         tex_xz_yy_0, tex_xz_yy_1, tex_xz_yz_0, tex_xz_yz_1, tex_xz_zz_0, tex_xz_zz_1, \
                                         tex_yy_x_0, tex_yy_x_1, tex_yy_xx_0, tex_yy_xx_1, tex_yy_xy_0, tex_yy_xy_1, \
                                         tex_yy_xz_0, tex_yy_xz_1, tex_yy_y_0, tex_yy_y_1, tex_yy_yy_0, tex_yy_yy_1, \
                                         tex_yy_yz_0, tex_yy_yz_1, tex_yy_z_0, tex_yy_z_1, tex_yy_zz_0, tex_yy_zz_1, \
                                         tex_yz_x_0, tex_yz_x_1, tex_yz_xx_0, tex_yz_xx_1, tex_yz_xy_0, tex_yz_xy_1, \
                                         tex_yz_xz_0, tex_yz_xz_1, tex_yz_y_0, tex_yz_y_1, tex_yz_yy_0, tex_yz_yy_1, \
                                         tex_yz_yz_0, tex_yz_yz_1, tex_yz_z_0, tex_yz_z_1, tex_yz_zz_0, tex_yz_zz_1, \
                                         tex_z_yy_0, tex_z_yy_1, tex_z_yz_0, tex_z_yz_1, tex_z_zz_0, tex_z_zz_1, \
                                         tey_xxz_yy_0, tey_xxz_yz_0, tey_xxz_zz_0, tey_xyy_xx_0, tey_xyy_xy_0, tey_xyy_xz_0, \
                                         tey_xyy_yy_0, tey_xyy_yz_0, tey_xyy_zz_0, tey_xyz_xx_0, tey_xyz_xy_0, tey_xyz_xz_0, \
                                         tey_xyz_yy_0, tey_xyz_yz_0, tey_xyz_zz_0, tey_xz_yy_0, tey_xz_yy_1, tey_xz_yz_0, \
                                         tey_xz_yz_1, tey_xz_zz_0, tey_xz_zz_1, tey_yy_x_0, tey_yy_x_1, tey_yy_xx_0, \
                                         tey_yy_xx_1, tey_yy_xy_0, tey_yy_xy_1, tey_yy_xz_0, tey_yy_xz_1, tey_yy_y_0, \
                                         tey_yy_y_1, tey_yy_yy_0, tey_yy_yy_1, tey_yy_yz_0, tey_yy_yz_1, tey_yy_z_0, \
                                         tey_yy_z_1, tey_yy_zz_0, tey_yy_zz_1, tey_yz_x_0, tey_yz_x_1, tey_yz_xx_0, \
                                         tey_yz_xx_1, tey_yz_xy_0, tey_yz_xy_1, tey_yz_xz_0, tey_yz_xz_1, tey_yz_y_0, \
                                         tey_yz_y_1, tey_yz_yy_0, tey_yz_yy_1, tey_yz_yz_0, tey_yz_yz_1, tey_yz_z_0, \
                                         tey_yz_z_1, tey_yz_zz_0, tey_yz_zz_1, tey_z_yy_0, tey_z_yy_1, tey_z_yz_0, \
                                         tey_z_yz_1, tey_z_zz_0, tey_z_zz_1, tez_xxz_yy_0, tez_xxz_yz_0, tez_xxz_zz_0, \
                                         tez_xyy_xx_0, tez_xyy_xy_0, tez_xyy_xz_0, tez_xyy_yy_0, tez_xyy_yz_0, tez_xyy_zz_0, \
                                         tez_xyz_xx_0, tez_xyz_xy_0, tez_xyz_xz_0, tez_xyz_yy_0, tez_xyz_yz_0, tez_xyz_zz_0, \
                                         tez_xz_yy_0, tez_xz_yy_1, tez_xz_yz_0, tez_xz_yz_1, tez_xz_zz_0, tez_xz_zz_1, \
                                         tez_yy_x_0, tez_yy_x_1, tez_yy_xx_0, tez_yy_xx_1, tez_yy_xy_0, tez_yy_xy_1, \
                                         tez_yy_xz_0, tez_yy_xz_1, tez_yy_y_0, tez_yy_y_1, tez_yy_yy_0, tez_yy_yy_1, \
                                         tez_yy_yz_0, tez_yy_yz_1, tez_yy_z_0, tez_yy_z_1, tez_yy_zz_0, tez_yy_zz_1, \
                                         tez_yz_x_0, tez_yz_x_1, tez_yz_xx_0, tez_yz_xx_1, tez_yz_xy_0, tez_yz_xy_1, \
                                         tez_yz_xz_0, tez_yz_xz_1, tez_yz_y_0, tez_yz_y_1, tez_yz_yy_0, tez_yz_yy_1, \
                                         tez_yz_yz_0, tez_yz_yz_1, tez_yz_z_0, tez_yz_z_1, tez_yz_zz_0, tez_yz_zz_1, \
                                         tez_z_yy_0, tez_z_yy_1, tez_z_yz_0, tez_z_yz_1, tez_z_zz_0, tez_z_zz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xxz_yy_0[j] = pa_x[j] * tex_xz_yy_0[j] - pc_x[j] * tex_xz_yy_1[j] + 0.5 * fl1_fx * tex_z_yy_0[j] - 0.5 * fl1_fx * tex_z_yy_1[j] + ta_xz_yy_1[j];

                    tey_xxz_yy_0[j] = pa_x[j] * tey_xz_yy_0[j] - pc_x[j] * tey_xz_yy_1[j] + 0.5 * fl1_fx * tey_z_yy_0[j] - 0.5 * fl1_fx * tey_z_yy_1[j];

                    tez_xxz_yy_0[j] = pa_x[j] * tez_xz_yy_0[j] - pc_x[j] * tez_xz_yy_1[j] + 0.5 * fl1_fx * tez_z_yy_0[j] - 0.5 * fl1_fx * tez_z_yy_1[j];

                    tex_xxz_yz_0[j] = pa_x[j] * tex_xz_yz_0[j] - pc_x[j] * tex_xz_yz_1[j] + 0.5 * fl1_fx * tex_z_yz_0[j] - 0.5 * fl1_fx * tex_z_yz_1[j] + ta_xz_yz_1[j];

                    tey_xxz_yz_0[j] = pa_x[j] * tey_xz_yz_0[j] - pc_x[j] * tey_xz_yz_1[j] + 0.5 * fl1_fx * tey_z_yz_0[j] - 0.5 * fl1_fx * tey_z_yz_1[j];

                    tez_xxz_yz_0[j] = pa_x[j] * tez_xz_yz_0[j] - pc_x[j] * tez_xz_yz_1[j] + 0.5 * fl1_fx * tez_z_yz_0[j] - 0.5 * fl1_fx * tez_z_yz_1[j];

                    tex_xxz_zz_0[j] = pa_x[j] * tex_xz_zz_0[j] - pc_x[j] * tex_xz_zz_1[j] + 0.5 * fl1_fx * tex_z_zz_0[j] - 0.5 * fl1_fx * tex_z_zz_1[j] + ta_xz_zz_1[j];

                    tey_xxz_zz_0[j] = pa_x[j] * tey_xz_zz_0[j] - pc_x[j] * tey_xz_zz_1[j] + 0.5 * fl1_fx * tey_z_zz_0[j] - 0.5 * fl1_fx * tey_z_zz_1[j];

                    tez_xxz_zz_0[j] = pa_x[j] * tez_xz_zz_0[j] - pc_x[j] * tez_xz_zz_1[j] + 0.5 * fl1_fx * tez_z_zz_0[j] - 0.5 * fl1_fx * tez_z_zz_1[j];

                    tex_xyy_xx_0[j] = pa_x[j] * tex_yy_xx_0[j] - pc_x[j] * tex_yy_xx_1[j] + fl1_fx * tex_yy_x_0[j] - fl1_fx * tex_yy_x_1[j] + ta_yy_xx_1[j];

                    tey_xyy_xx_0[j] = pa_x[j] * tey_yy_xx_0[j] - pc_x[j] * tey_yy_xx_1[j] + fl1_fx * tey_yy_x_0[j] - fl1_fx * tey_yy_x_1[j];

                    tez_xyy_xx_0[j] = pa_x[j] * tez_yy_xx_0[j] - pc_x[j] * tez_yy_xx_1[j] + fl1_fx * tez_yy_x_0[j] - fl1_fx * tez_yy_x_1[j];

                    tex_xyy_xy_0[j] = pa_x[j] * tex_yy_xy_0[j] - pc_x[j] * tex_yy_xy_1[j] + 0.5 * fl1_fx * tex_yy_y_0[j] - 0.5 * fl1_fx * tex_yy_y_1[j] + ta_yy_xy_1[j];

                    tey_xyy_xy_0[j] = pa_x[j] * tey_yy_xy_0[j] - pc_x[j] * tey_yy_xy_1[j] + 0.5 * fl1_fx * tey_yy_y_0[j] - 0.5 * fl1_fx * tey_yy_y_1[j];

                    tez_xyy_xy_0[j] = pa_x[j] * tez_yy_xy_0[j] - pc_x[j] * tez_yy_xy_1[j] + 0.5 * fl1_fx * tez_yy_y_0[j] - 0.5 * fl1_fx * tez_yy_y_1[j];

                    tex_xyy_xz_0[j] = pa_x[j] * tex_yy_xz_0[j] - pc_x[j] * tex_yy_xz_1[j] + 0.5 * fl1_fx * tex_yy_z_0[j] - 0.5 * fl1_fx * tex_yy_z_1[j] + ta_yy_xz_1[j];

                    tey_xyy_xz_0[j] = pa_x[j] * tey_yy_xz_0[j] - pc_x[j] * tey_yy_xz_1[j] + 0.5 * fl1_fx * tey_yy_z_0[j] - 0.5 * fl1_fx * tey_yy_z_1[j];

                    tez_xyy_xz_0[j] = pa_x[j] * tez_yy_xz_0[j] - pc_x[j] * tez_yy_xz_1[j] + 0.5 * fl1_fx * tez_yy_z_0[j] - 0.5 * fl1_fx * tez_yy_z_1[j];

                    tex_xyy_yy_0[j] = pa_x[j] * tex_yy_yy_0[j] - pc_x[j] * tex_yy_yy_1[j] + ta_yy_yy_1[j];

                    tey_xyy_yy_0[j] = pa_x[j] * tey_yy_yy_0[j] - pc_x[j] * tey_yy_yy_1[j];

                    tez_xyy_yy_0[j] = pa_x[j] * tez_yy_yy_0[j] - pc_x[j] * tez_yy_yy_1[j];

                    tex_xyy_yz_0[j] = pa_x[j] * tex_yy_yz_0[j] - pc_x[j] * tex_yy_yz_1[j] + ta_yy_yz_1[j];

                    tey_xyy_yz_0[j] = pa_x[j] * tey_yy_yz_0[j] - pc_x[j] * tey_yy_yz_1[j];

                    tez_xyy_yz_0[j] = pa_x[j] * tez_yy_yz_0[j] - pc_x[j] * tez_yy_yz_1[j];

                    tex_xyy_zz_0[j] = pa_x[j] * tex_yy_zz_0[j] - pc_x[j] * tex_yy_zz_1[j] + ta_yy_zz_1[j];

                    tey_xyy_zz_0[j] = pa_x[j] * tey_yy_zz_0[j] - pc_x[j] * tey_yy_zz_1[j];

                    tez_xyy_zz_0[j] = pa_x[j] * tez_yy_zz_0[j] - pc_x[j] * tez_yy_zz_1[j];

                    tex_xyz_xx_0[j] = pa_x[j] * tex_yz_xx_0[j] - pc_x[j] * tex_yz_xx_1[j] + fl1_fx * tex_yz_x_0[j] - fl1_fx * tex_yz_x_1[j] + ta_yz_xx_1[j];

                    tey_xyz_xx_0[j] = pa_x[j] * tey_yz_xx_0[j] - pc_x[j] * tey_yz_xx_1[j] + fl1_fx * tey_yz_x_0[j] - fl1_fx * tey_yz_x_1[j];

                    tez_xyz_xx_0[j] = pa_x[j] * tez_yz_xx_0[j] - pc_x[j] * tez_yz_xx_1[j] + fl1_fx * tez_yz_x_0[j] - fl1_fx * tez_yz_x_1[j];

                    tex_xyz_xy_0[j] = pa_x[j] * tex_yz_xy_0[j] - pc_x[j] * tex_yz_xy_1[j] + 0.5 * fl1_fx * tex_yz_y_0[j] - 0.5 * fl1_fx * tex_yz_y_1[j] + ta_yz_xy_1[j];

                    tey_xyz_xy_0[j] = pa_x[j] * tey_yz_xy_0[j] - pc_x[j] * tey_yz_xy_1[j] + 0.5 * fl1_fx * tey_yz_y_0[j] - 0.5 * fl1_fx * tey_yz_y_1[j];

                    tez_xyz_xy_0[j] = pa_x[j] * tez_yz_xy_0[j] - pc_x[j] * tez_yz_xy_1[j] + 0.5 * fl1_fx * tez_yz_y_0[j] - 0.5 * fl1_fx * tez_yz_y_1[j];

                    tex_xyz_xz_0[j] = pa_x[j] * tex_yz_xz_0[j] - pc_x[j] * tex_yz_xz_1[j] + 0.5 * fl1_fx * tex_yz_z_0[j] - 0.5 * fl1_fx * tex_yz_z_1[j] + ta_yz_xz_1[j];

                    tey_xyz_xz_0[j] = pa_x[j] * tey_yz_xz_0[j] - pc_x[j] * tey_yz_xz_1[j] + 0.5 * fl1_fx * tey_yz_z_0[j] - 0.5 * fl1_fx * tey_yz_z_1[j];

                    tez_xyz_xz_0[j] = pa_x[j] * tez_yz_xz_0[j] - pc_x[j] * tez_yz_xz_1[j] + 0.5 * fl1_fx * tez_yz_z_0[j] - 0.5 * fl1_fx * tez_yz_z_1[j];

                    tex_xyz_yy_0[j] = pa_x[j] * tex_yz_yy_0[j] - pc_x[j] * tex_yz_yy_1[j] + ta_yz_yy_1[j];

                    tey_xyz_yy_0[j] = pa_x[j] * tey_yz_yy_0[j] - pc_x[j] * tey_yz_yy_1[j];

                    tez_xyz_yy_0[j] = pa_x[j] * tez_yz_yy_0[j] - pc_x[j] * tez_yz_yy_1[j];

                    tex_xyz_yz_0[j] = pa_x[j] * tex_yz_yz_0[j] - pc_x[j] * tex_yz_yz_1[j] + ta_yz_yz_1[j];

                    tey_xyz_yz_0[j] = pa_x[j] * tey_yz_yz_0[j] - pc_x[j] * tey_yz_yz_1[j];

                    tez_xyz_yz_0[j] = pa_x[j] * tez_yz_yz_0[j] - pc_x[j] * tez_yz_yz_1[j];

                    tex_xyz_zz_0[j] = pa_x[j] * tex_yz_zz_0[j] - pc_x[j] * tex_yz_zz_1[j] + ta_yz_zz_1[j];

                    tey_xyz_zz_0[j] = pa_x[j] * tey_yz_zz_0[j] - pc_x[j] * tey_yz_zz_1[j];

                    tez_xyz_zz_0[j] = pa_x[j] * tez_yz_zz_0[j] - pc_x[j] * tez_yz_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFD_90_135(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (90,135)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tex_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 18); 

                auto tey_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 19); 

                auto tey_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 20); 

                auto tey_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 21); 

                auto tey_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 22); 

                auto tey_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 23); 

                auto tey_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 24); 

                auto tey_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 25); 

                auto tey_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 26); 

                auto tey_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 26); 

                auto tex_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 30); 

                auto tey_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 31); 

                auto tey_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 32); 

                auto tey_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 33); 

                auto tey_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 34); 

                auto tey_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 35); 

                auto tey_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 35); 

                auto tex_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 18); 

                auto tey_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 19); 

                auto tey_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 20); 

                auto tey_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 21); 

                auto tey_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 22); 

                auto tey_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 23); 

                auto tey_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 24); 

                auto tey_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 25); 

                auto tey_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 26); 

                auto tey_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 26); 

                auto tex_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 30); 

                auto tey_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 31); 

                auto tey_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 32); 

                auto tey_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 33); 

                auto tey_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 34); 

                auto tey_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 35); 

                auto tey_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 35); 

                auto tex_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 6); 

                auto tey_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 7); 

                auto tey_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 8); 

                auto tey_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 9); 

                auto tey_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 10); 

                auto tey_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 11); 

                auto tey_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 6); 

                auto tey_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 6); 

                auto tez_y_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 6); 

                auto tex_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 7); 

                auto tey_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 7); 

                auto tez_y_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 7); 

                auto tex_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 8); 

                auto tey_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 8); 

                auto tez_y_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 8); 

                auto tex_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 9); 

                auto tey_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_y_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 10); 

                auto tey_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_y_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 10); 

                auto tex_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 11); 

                auto tey_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_y_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 9); 

                auto tey_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 9); 

                auto tez_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 9); 

                auto tex_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 10); 

                auto tey_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 10); 

                auto tez_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 10); 

                auto tex_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 11); 

                auto tey_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 11); 

                auto tez_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 11); 

                auto tex_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 12); 

                auto tey_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 15); 

                auto tey_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 16); 

                auto tey_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 17); 

                auto tey_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 9); 

                auto tey_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 9); 

                auto tez_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 9); 

                auto tex_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 10); 

                auto tey_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 10); 

                auto tez_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 10); 

                auto tex_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 11); 

                auto tey_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 11); 

                auto tez_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 11); 

                auto tex_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 12); 

                auto tey_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 15); 

                auto tey_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 16); 

                auto tey_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 17); 

                auto tey_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 17); 

                auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18); 

                auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19); 

                auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20); 

                auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21); 

                auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22); 

                auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23); 

                auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24); 

                auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25); 

                auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26); 

                auto ta_zz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 30); 

                auto ta_zz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 31); 

                auto ta_zz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 32); 

                auto ta_zz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 33); 

                auto ta_zz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 34); 

                auto ta_zz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 35); 

                // set up pointers to integrals

                auto tex_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 30); 

                auto tey_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 30); 

                auto tez_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 30); 

                auto tex_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 31); 

                auto tey_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 31); 

                auto tez_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 31); 

                auto tex_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 32); 

                auto tey_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 32); 

                auto tez_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 32); 

                auto tex_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 33); 

                auto tey_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 33); 

                auto tez_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 33); 

                auto tex_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 34); 

                auto tey_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 34); 

                auto tez_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 34); 

                auto tex_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 35); 

                auto tey_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 35); 

                auto tez_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 35); 

                auto tex_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 36); 

                auto tey_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 36); 

                auto tez_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 36); 

                auto tex_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 37); 

                auto tey_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 37); 

                auto tez_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 37); 

                auto tex_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 38); 

                auto tey_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 38); 

                auto tez_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 38); 

                auto tex_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 39); 

                auto tey_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 39); 

                auto tez_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 39); 

                auto tex_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 40); 

                auto tey_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 40); 

                auto tez_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 40); 

                auto tex_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 41); 

                auto tey_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 41); 

                auto tez_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 41); 

                auto tex_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 42); 

                auto tey_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 42); 

                auto tez_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 42); 

                auto tex_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 43); 

                auto tey_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 43); 

                auto tez_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 43); 

                auto tex_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 44); 

                auto tey_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 44); 

                auto tez_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 44); 

                // Batch of Integrals (90,135)

                #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_yy_xx_1, ta_yy_xy_1, ta_yy_xz_1, ta_yy_yy_1, \
                                         ta_yy_yz_1, ta_yy_zz_1, ta_yz_xx_1, ta_yz_xy_1, ta_yz_xz_1, ta_zz_xx_1, ta_zz_xy_1, \
                                         ta_zz_xz_1, ta_zz_yy_1, ta_zz_yz_1, ta_zz_zz_1, tex_xzz_xx_0, tex_xzz_xy_0, \
                                         tex_xzz_xz_0, tex_xzz_yy_0, tex_xzz_yz_0, tex_xzz_zz_0, tex_y_xx_0, tex_y_xx_1, \
                                         tex_y_xy_0, tex_y_xy_1, tex_y_xz_0, tex_y_xz_1, tex_y_yy_0, tex_y_yy_1, tex_y_yz_0, \
                                         tex_y_yz_1, tex_y_zz_0, tex_y_zz_1, tex_yy_x_0, tex_yy_x_1, tex_yy_xx_0, \
                                         tex_yy_xx_1, tex_yy_xy_0, tex_yy_xy_1, tex_yy_xz_0, tex_yy_xz_1, tex_yy_y_0, \
                                         tex_yy_y_1, tex_yy_yy_0, tex_yy_yy_1, tex_yy_yz_0, tex_yy_yz_1, tex_yy_z_0, \
                                         tex_yy_z_1, tex_yy_zz_0, tex_yy_zz_1, tex_yyy_xx_0, tex_yyy_xy_0, tex_yyy_xz_0, \
                                         tex_yyy_yy_0, tex_yyy_yz_0, tex_yyy_zz_0, tex_yyz_xx_0, tex_yyz_xy_0, tex_yyz_xz_0, \
                                         tex_yz_x_0, tex_yz_x_1, tex_yz_xx_0, tex_yz_xx_1, tex_yz_xy_0, tex_yz_xy_1, \
                                         tex_yz_xz_0, tex_yz_xz_1, tex_z_xx_0, tex_z_xx_1, tex_z_xy_0, tex_z_xy_1, \
                                         tex_z_xz_0, tex_z_xz_1, tex_zz_x_0, tex_zz_x_1, tex_zz_xx_0, tex_zz_xx_1, \
                                         tex_zz_xy_0, tex_zz_xy_1, tex_zz_xz_0, tex_zz_xz_1, tex_zz_y_0, tex_zz_y_1, \
                                         tex_zz_yy_0, tex_zz_yy_1, tex_zz_yz_0, tex_zz_yz_1, tex_zz_z_0, tex_zz_z_1, \
                                         tex_zz_zz_0, tex_zz_zz_1, tey_xzz_xx_0, tey_xzz_xy_0, tey_xzz_xz_0, tey_xzz_yy_0, \
                                         tey_xzz_yz_0, tey_xzz_zz_0, tey_y_xx_0, tey_y_xx_1, tey_y_xy_0, tey_y_xy_1, \
                                         tey_y_xz_0, tey_y_xz_1, tey_y_yy_0, tey_y_yy_1, tey_y_yz_0, tey_y_yz_1, tey_y_zz_0, \
                                         tey_y_zz_1, tey_yy_x_0, tey_yy_x_1, tey_yy_xx_0, tey_yy_xx_1, tey_yy_xy_0, \
                                         tey_yy_xy_1, tey_yy_xz_0, tey_yy_xz_1, tey_yy_y_0, tey_yy_y_1, tey_yy_yy_0, \
                                         tey_yy_yy_1, tey_yy_yz_0, tey_yy_yz_1, tey_yy_z_0, tey_yy_z_1, tey_yy_zz_0, \
                                         tey_yy_zz_1, tey_yyy_xx_0, tey_yyy_xy_0, tey_yyy_xz_0, tey_yyy_yy_0, tey_yyy_yz_0, \
                                         tey_yyy_zz_0, tey_yyz_xx_0, tey_yyz_xy_0, tey_yyz_xz_0, tey_yz_x_0, tey_yz_x_1, \
                                         tey_yz_xx_0, tey_yz_xx_1, tey_yz_xy_0, tey_yz_xy_1, tey_yz_xz_0, tey_yz_xz_1, \
                                         tey_z_xx_0, tey_z_xx_1, tey_z_xy_0, tey_z_xy_1, tey_z_xz_0, tey_z_xz_1, tey_zz_x_0, \
                                         tey_zz_x_1, tey_zz_xx_0, tey_zz_xx_1, tey_zz_xy_0, tey_zz_xy_1, tey_zz_xz_0, \
                                         tey_zz_xz_1, tey_zz_y_0, tey_zz_y_1, tey_zz_yy_0, tey_zz_yy_1, tey_zz_yz_0, \
                                         tey_zz_yz_1, tey_zz_z_0, tey_zz_z_1, tey_zz_zz_0, tey_zz_zz_1, tez_xzz_xx_0, \
                                         tez_xzz_xy_0, tez_xzz_xz_0, tez_xzz_yy_0, tez_xzz_yz_0, tez_xzz_zz_0, tez_y_xx_0, \
                                         tez_y_xx_1, tez_y_xy_0, tez_y_xy_1, tez_y_xz_0, tez_y_xz_1, tez_y_yy_0, tez_y_yy_1, \
                                         tez_y_yz_0, tez_y_yz_1, tez_y_zz_0, tez_y_zz_1, tez_yy_x_0, tez_yy_x_1, \
                                         tez_yy_xx_0, tez_yy_xx_1, tez_yy_xy_0, tez_yy_xy_1, tez_yy_xz_0, tez_yy_xz_1, \
                                         tez_yy_y_0, tez_yy_y_1, tez_yy_yy_0, tez_yy_yy_1, tez_yy_yz_0, tez_yy_yz_1, \
                                         tez_yy_z_0, tez_yy_z_1, tez_yy_zz_0, tez_yy_zz_1, tez_yyy_xx_0, tez_yyy_xy_0, \
                                         tez_yyy_xz_0, tez_yyy_yy_0, tez_yyy_yz_0, tez_yyy_zz_0, tez_yyz_xx_0, tez_yyz_xy_0, \
                                         tez_yyz_xz_0, tez_yz_x_0, tez_yz_x_1, tez_yz_xx_0, tez_yz_xx_1, tez_yz_xy_0, \
                                         tez_yz_xy_1, tez_yz_xz_0, tez_yz_xz_1, tez_z_xx_0, tez_z_xx_1, tez_z_xy_0, \
                                         tez_z_xy_1, tez_z_xz_0, tez_z_xz_1, tez_zz_x_0, tez_zz_x_1, tez_zz_xx_0, \
                                         tez_zz_xx_1, tez_zz_xy_0, tez_zz_xy_1, tez_zz_xz_0, tez_zz_xz_1, tez_zz_y_0, \
                                         tez_zz_y_1, tez_zz_yy_0, tez_zz_yy_1, tez_zz_yz_0, tez_zz_yz_1, tez_zz_z_0, \
                                         tez_zz_z_1, tez_zz_zz_0, tez_zz_zz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xzz_xx_0[j] = pa_x[j] * tex_zz_xx_0[j] - pc_x[j] * tex_zz_xx_1[j] + fl1_fx * tex_zz_x_0[j] - fl1_fx * tex_zz_x_1[j] + ta_zz_xx_1[j];

                    tey_xzz_xx_0[j] = pa_x[j] * tey_zz_xx_0[j] - pc_x[j] * tey_zz_xx_1[j] + fl1_fx * tey_zz_x_0[j] - fl1_fx * tey_zz_x_1[j];

                    tez_xzz_xx_0[j] = pa_x[j] * tez_zz_xx_0[j] - pc_x[j] * tez_zz_xx_1[j] + fl1_fx * tez_zz_x_0[j] - fl1_fx * tez_zz_x_1[j];

                    tex_xzz_xy_0[j] = pa_x[j] * tex_zz_xy_0[j] - pc_x[j] * tex_zz_xy_1[j] + 0.5 * fl1_fx * tex_zz_y_0[j] - 0.5 * fl1_fx * tex_zz_y_1[j] + ta_zz_xy_1[j];

                    tey_xzz_xy_0[j] = pa_x[j] * tey_zz_xy_0[j] - pc_x[j] * tey_zz_xy_1[j] + 0.5 * fl1_fx * tey_zz_y_0[j] - 0.5 * fl1_fx * tey_zz_y_1[j];

                    tez_xzz_xy_0[j] = pa_x[j] * tez_zz_xy_0[j] - pc_x[j] * tez_zz_xy_1[j] + 0.5 * fl1_fx * tez_zz_y_0[j] - 0.5 * fl1_fx * tez_zz_y_1[j];

                    tex_xzz_xz_0[j] = pa_x[j] * tex_zz_xz_0[j] - pc_x[j] * tex_zz_xz_1[j] + 0.5 * fl1_fx * tex_zz_z_0[j] - 0.5 * fl1_fx * tex_zz_z_1[j] + ta_zz_xz_1[j];

                    tey_xzz_xz_0[j] = pa_x[j] * tey_zz_xz_0[j] - pc_x[j] * tey_zz_xz_1[j] + 0.5 * fl1_fx * tey_zz_z_0[j] - 0.5 * fl1_fx * tey_zz_z_1[j];

                    tez_xzz_xz_0[j] = pa_x[j] * tez_zz_xz_0[j] - pc_x[j] * tez_zz_xz_1[j] + 0.5 * fl1_fx * tez_zz_z_0[j] - 0.5 * fl1_fx * tez_zz_z_1[j];

                    tex_xzz_yy_0[j] = pa_x[j] * tex_zz_yy_0[j] - pc_x[j] * tex_zz_yy_1[j] + ta_zz_yy_1[j];

                    tey_xzz_yy_0[j] = pa_x[j] * tey_zz_yy_0[j] - pc_x[j] * tey_zz_yy_1[j];

                    tez_xzz_yy_0[j] = pa_x[j] * tez_zz_yy_0[j] - pc_x[j] * tez_zz_yy_1[j];

                    tex_xzz_yz_0[j] = pa_x[j] * tex_zz_yz_0[j] - pc_x[j] * tex_zz_yz_1[j] + ta_zz_yz_1[j];

                    tey_xzz_yz_0[j] = pa_x[j] * tey_zz_yz_0[j] - pc_x[j] * tey_zz_yz_1[j];

                    tez_xzz_yz_0[j] = pa_x[j] * tez_zz_yz_0[j] - pc_x[j] * tez_zz_yz_1[j];

                    tex_xzz_zz_0[j] = pa_x[j] * tex_zz_zz_0[j] - pc_x[j] * tex_zz_zz_1[j] + ta_zz_zz_1[j];

                    tey_xzz_zz_0[j] = pa_x[j] * tey_zz_zz_0[j] - pc_x[j] * tey_zz_zz_1[j];

                    tez_xzz_zz_0[j] = pa_x[j] * tez_zz_zz_0[j] - pc_x[j] * tez_zz_zz_1[j];

                    tex_yyy_xx_0[j] = pa_y[j] * tex_yy_xx_0[j] - pc_y[j] * tex_yy_xx_1[j] + fl1_fx * tex_y_xx_0[j] - fl1_fx * tex_y_xx_1[j];

                    tey_yyy_xx_0[j] = pa_y[j] * tey_yy_xx_0[j] - pc_y[j] * tey_yy_xx_1[j] + fl1_fx * tey_y_xx_0[j] - fl1_fx * tey_y_xx_1[j] + ta_yy_xx_1[j];

                    tez_yyy_xx_0[j] = pa_y[j] * tez_yy_xx_0[j] - pc_y[j] * tez_yy_xx_1[j] + fl1_fx * tez_y_xx_0[j] - fl1_fx * tez_y_xx_1[j];

                    tex_yyy_xy_0[j] = pa_y[j] * tex_yy_xy_0[j] - pc_y[j] * tex_yy_xy_1[j] + fl1_fx * tex_y_xy_0[j] - fl1_fx * tex_y_xy_1[j] + 0.5 * fl1_fx * tex_yy_x_0[j] - 0.5 * fl1_fx * tex_yy_x_1[j];

                    tey_yyy_xy_0[j] = pa_y[j] * tey_yy_xy_0[j] - pc_y[j] * tey_yy_xy_1[j] + fl1_fx * tey_y_xy_0[j] - fl1_fx * tey_y_xy_1[j] + 0.5 * fl1_fx * tey_yy_x_0[j] - 0.5 * fl1_fx * tey_yy_x_1[j] + ta_yy_xy_1[j];

                    tez_yyy_xy_0[j] = pa_y[j] * tez_yy_xy_0[j] - pc_y[j] * tez_yy_xy_1[j] + fl1_fx * tez_y_xy_0[j] - fl1_fx * tez_y_xy_1[j] + 0.5 * fl1_fx * tez_yy_x_0[j] - 0.5 * fl1_fx * tez_yy_x_1[j];

                    tex_yyy_xz_0[j] = pa_y[j] * tex_yy_xz_0[j] - pc_y[j] * tex_yy_xz_1[j] + fl1_fx * tex_y_xz_0[j] - fl1_fx * tex_y_xz_1[j];

                    tey_yyy_xz_0[j] = pa_y[j] * tey_yy_xz_0[j] - pc_y[j] * tey_yy_xz_1[j] + fl1_fx * tey_y_xz_0[j] - fl1_fx * tey_y_xz_1[j] + ta_yy_xz_1[j];

                    tez_yyy_xz_0[j] = pa_y[j] * tez_yy_xz_0[j] - pc_y[j] * tez_yy_xz_1[j] + fl1_fx * tez_y_xz_0[j] - fl1_fx * tez_y_xz_1[j];

                    tex_yyy_yy_0[j] = pa_y[j] * tex_yy_yy_0[j] - pc_y[j] * tex_yy_yy_1[j] + fl1_fx * tex_y_yy_0[j] - fl1_fx * tex_y_yy_1[j] + fl1_fx * tex_yy_y_0[j] - fl1_fx * tex_yy_y_1[j];

                    tey_yyy_yy_0[j] = pa_y[j] * tey_yy_yy_0[j] - pc_y[j] * tey_yy_yy_1[j] + fl1_fx * tey_y_yy_0[j] - fl1_fx * tey_y_yy_1[j] + fl1_fx * tey_yy_y_0[j] - fl1_fx * tey_yy_y_1[j] + ta_yy_yy_1[j];

                    tez_yyy_yy_0[j] = pa_y[j] * tez_yy_yy_0[j] - pc_y[j] * tez_yy_yy_1[j] + fl1_fx * tez_y_yy_0[j] - fl1_fx * tez_y_yy_1[j] + fl1_fx * tez_yy_y_0[j] - fl1_fx * tez_yy_y_1[j];

                    tex_yyy_yz_0[j] = pa_y[j] * tex_yy_yz_0[j] - pc_y[j] * tex_yy_yz_1[j] + fl1_fx * tex_y_yz_0[j] - fl1_fx * tex_y_yz_1[j] + 0.5 * fl1_fx * tex_yy_z_0[j] - 0.5 * fl1_fx * tex_yy_z_1[j];

                    tey_yyy_yz_0[j] = pa_y[j] * tey_yy_yz_0[j] - pc_y[j] * tey_yy_yz_1[j] + fl1_fx * tey_y_yz_0[j] - fl1_fx * tey_y_yz_1[j] + 0.5 * fl1_fx * tey_yy_z_0[j] - 0.5 * fl1_fx * tey_yy_z_1[j] + ta_yy_yz_1[j];

                    tez_yyy_yz_0[j] = pa_y[j] * tez_yy_yz_0[j] - pc_y[j] * tez_yy_yz_1[j] + fl1_fx * tez_y_yz_0[j] - fl1_fx * tez_y_yz_1[j] + 0.5 * fl1_fx * tez_yy_z_0[j] - 0.5 * fl1_fx * tez_yy_z_1[j];

                    tex_yyy_zz_0[j] = pa_y[j] * tex_yy_zz_0[j] - pc_y[j] * tex_yy_zz_1[j] + fl1_fx * tex_y_zz_0[j] - fl1_fx * tex_y_zz_1[j];

                    tey_yyy_zz_0[j] = pa_y[j] * tey_yy_zz_0[j] - pc_y[j] * tey_yy_zz_1[j] + fl1_fx * tey_y_zz_0[j] - fl1_fx * tey_y_zz_1[j] + ta_yy_zz_1[j];

                    tez_yyy_zz_0[j] = pa_y[j] * tez_yy_zz_0[j] - pc_y[j] * tez_yy_zz_1[j] + fl1_fx * tez_y_zz_0[j] - fl1_fx * tez_y_zz_1[j];

                    tex_yyz_xx_0[j] = pa_y[j] * tex_yz_xx_0[j] - pc_y[j] * tex_yz_xx_1[j] + 0.5 * fl1_fx * tex_z_xx_0[j] - 0.5 * fl1_fx * tex_z_xx_1[j];

                    tey_yyz_xx_0[j] = pa_y[j] * tey_yz_xx_0[j] - pc_y[j] * tey_yz_xx_1[j] + 0.5 * fl1_fx * tey_z_xx_0[j] - 0.5 * fl1_fx * tey_z_xx_1[j] + ta_yz_xx_1[j];

                    tez_yyz_xx_0[j] = pa_y[j] * tez_yz_xx_0[j] - pc_y[j] * tez_yz_xx_1[j] + 0.5 * fl1_fx * tez_z_xx_0[j] - 0.5 * fl1_fx * tez_z_xx_1[j];

                    tex_yyz_xy_0[j] = pa_y[j] * tex_yz_xy_0[j] - pc_y[j] * tex_yz_xy_1[j] + 0.5 * fl1_fx * tex_z_xy_0[j] - 0.5 * fl1_fx * tex_z_xy_1[j] + 0.5 * fl1_fx * tex_yz_x_0[j] - 0.5 * fl1_fx * tex_yz_x_1[j];

                    tey_yyz_xy_0[j] = pa_y[j] * tey_yz_xy_0[j] - pc_y[j] * tey_yz_xy_1[j] + 0.5 * fl1_fx * tey_z_xy_0[j] - 0.5 * fl1_fx * tey_z_xy_1[j] + 0.5 * fl1_fx * tey_yz_x_0[j] - 0.5 * fl1_fx * tey_yz_x_1[j] + ta_yz_xy_1[j];

                    tez_yyz_xy_0[j] = pa_y[j] * tez_yz_xy_0[j] - pc_y[j] * tez_yz_xy_1[j] + 0.5 * fl1_fx * tez_z_xy_0[j] - 0.5 * fl1_fx * tez_z_xy_1[j] + 0.5 * fl1_fx * tez_yz_x_0[j] - 0.5 * fl1_fx * tez_yz_x_1[j];

                    tex_yyz_xz_0[j] = pa_y[j] * tex_yz_xz_0[j] - pc_y[j] * tex_yz_xz_1[j] + 0.5 * fl1_fx * tex_z_xz_0[j] - 0.5 * fl1_fx * tex_z_xz_1[j];

                    tey_yyz_xz_0[j] = pa_y[j] * tey_yz_xz_0[j] - pc_y[j] * tey_yz_xz_1[j] + 0.5 * fl1_fx * tey_z_xz_0[j] - 0.5 * fl1_fx * tey_z_xz_1[j] + ta_yz_xz_1[j];

                    tez_yyz_xz_0[j] = pa_y[j] * tez_yz_xz_0[j] - pc_y[j] * tez_yz_xz_1[j] + 0.5 * fl1_fx * tez_z_xz_0[j] - 0.5 * fl1_fx * tez_z_xz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForFD_135_180(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (135,180)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {3, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_3_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
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

                auto tex_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 27); 

                auto tey_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 28); 

                auto tey_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 29); 

                auto tey_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 29); 

                auto tex_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 30); 

                auto tey_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 31); 

                auto tey_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 32); 

                auto tey_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 33); 

                auto tey_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 34); 

                auto tey_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 35); 

                auto tey_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 35); 

                auto tex_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 27); 

                auto tey_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 28); 

                auto tey_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 29); 

                auto tey_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 29); 

                auto tex_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 30); 

                auto tey_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 31); 

                auto tey_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 32); 

                auto tey_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 33); 

                auto tey_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 34); 

                auto tey_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 35); 

                auto tey_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 35); 

                auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12); 

                auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13); 

                auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14); 

                auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 15); 

                auto tey_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 16); 

                auto tey_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 17); 

                auto tey_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 12); 

                auto tey_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 12); 

                auto tez_z_xx_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 12); 

                auto tex_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 13); 

                auto tey_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_z_xy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 14); 

                auto tey_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_z_xz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 15); 

                auto tey_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_z_yy_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 16); 

                auto tey_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_z_yz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * idx + 17); 

                auto tey_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_z_zz_1 = primBuffer.data(pidx_e_1_2_m1 + 36 * bdim + 18 * idx + 17); 

                auto tex_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 13); 

                auto tey_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 13); 

                auto tez_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 13); 

                auto tex_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 14); 

                auto tey_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 14); 

                auto tez_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 14); 

                auto tex_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 15); 

                auto tey_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 15); 

                auto tez_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 15); 

                auto tex_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 16); 

                auto tey_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 16); 

                auto tez_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 16); 

                auto tex_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 17); 

                auto tey_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 17); 

                auto tez_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 17); 

                auto tex_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 13); 

                auto tey_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 13); 

                auto tez_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 13); 

                auto tex_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 14); 

                auto tey_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 14); 

                auto tez_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 14); 

                auto tex_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 15); 

                auto tey_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 15); 

                auto tez_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 15); 

                auto tex_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 16); 

                auto tey_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 16); 

                auto tez_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 16); 

                auto tex_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 17); 

                auto tey_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 17); 

                auto tez_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 17); 

                auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27); 

                auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28); 

                auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29); 

                auto ta_zz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 30); 

                auto ta_zz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 31); 

                auto ta_zz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 32); 

                auto ta_zz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 33); 

                auto ta_zz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 34); 

                auto ta_zz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 35); 

                // set up pointers to integrals

                auto tex_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 45); 

                auto tey_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 45); 

                auto tez_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 45); 

                auto tex_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 46); 

                auto tey_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 46); 

                auto tez_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 46); 

                auto tex_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 47); 

                auto tey_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 47); 

                auto tez_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 47); 

                auto tex_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 48); 

                auto tey_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 48); 

                auto tez_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 48); 

                auto tex_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 49); 

                auto tey_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 49); 

                auto tez_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 49); 

                auto tex_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 50); 

                auto tey_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 50); 

                auto tez_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 50); 

                auto tex_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 51); 

                auto tey_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 51); 

                auto tez_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 51); 

                auto tex_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 52); 

                auto tey_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 52); 

                auto tez_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 52); 

                auto tex_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 53); 

                auto tey_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 53); 

                auto tez_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 53); 

                auto tex_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 54); 

                auto tey_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 54); 

                auto tez_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 54); 

                auto tex_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 55); 

                auto tey_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 55); 

                auto tez_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 55); 

                auto tex_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 56); 

                auto tey_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 56); 

                auto tez_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 56); 

                auto tex_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 57); 

                auto tey_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 57); 

                auto tez_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 57); 

                auto tex_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 58); 

                auto tey_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 58); 

                auto tez_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 58); 

                auto tex_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 59); 

                auto tey_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 59); 

                auto tez_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 59); 

                // Batch of Integrals (135,180)

                #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_yz_yy_1, ta_yz_yz_1, ta_yz_zz_1, ta_zz_xx_1, \
                                         ta_zz_xy_1, ta_zz_xz_1, ta_zz_yy_1, ta_zz_yz_1, ta_zz_zz_1, tex_yyz_yy_0, \
                                         tex_yyz_yz_0, tex_yyz_zz_0, tex_yz_y_0, tex_yz_y_1, tex_yz_yy_0, tex_yz_yy_1, \
                                         tex_yz_yz_0, tex_yz_yz_1, tex_yz_z_0, tex_yz_z_1, tex_yz_zz_0, tex_yz_zz_1, \
                                         tex_yzz_xx_0, tex_yzz_xy_0, tex_yzz_xz_0, tex_yzz_yy_0, tex_yzz_yz_0, tex_yzz_zz_0, \
                                         tex_z_xx_0, tex_z_xx_1, tex_z_xy_0, tex_z_xy_1, tex_z_xz_0, tex_z_xz_1, tex_z_yy_0, \
                                         tex_z_yy_1, tex_z_yz_0, tex_z_yz_1, tex_z_zz_0, tex_z_zz_1, tex_zz_x_0, tex_zz_x_1, \
                                         tex_zz_xx_0, tex_zz_xx_1, tex_zz_xy_0, tex_zz_xy_1, tex_zz_xz_0, tex_zz_xz_1, \
                                         tex_zz_y_0, tex_zz_y_1, tex_zz_yy_0, tex_zz_yy_1, tex_zz_yz_0, tex_zz_yz_1, \
                                         tex_zz_z_0, tex_zz_z_1, tex_zz_zz_0, tex_zz_zz_1, tex_zzz_xx_0, tex_zzz_xy_0, \
                                         tex_zzz_xz_0, tex_zzz_yy_0, tex_zzz_yz_0, tex_zzz_zz_0, tey_yyz_yy_0, tey_yyz_yz_0, \
                                         tey_yyz_zz_0, tey_yz_y_0, tey_yz_y_1, tey_yz_yy_0, tey_yz_yy_1, tey_yz_yz_0, \
                                         tey_yz_yz_1, tey_yz_z_0, tey_yz_z_1, tey_yz_zz_0, tey_yz_zz_1, tey_yzz_xx_0, \
                                         tey_yzz_xy_0, tey_yzz_xz_0, tey_yzz_yy_0, tey_yzz_yz_0, tey_yzz_zz_0, tey_z_xx_0, \
                                         tey_z_xx_1, tey_z_xy_0, tey_z_xy_1, tey_z_xz_0, tey_z_xz_1, tey_z_yy_0, tey_z_yy_1, \
                                         tey_z_yz_0, tey_z_yz_1, tey_z_zz_0, tey_z_zz_1, tey_zz_x_0, tey_zz_x_1, \
                                         tey_zz_xx_0, tey_zz_xx_1, tey_zz_xy_0, tey_zz_xy_1, tey_zz_xz_0, tey_zz_xz_1, \
                                         tey_zz_y_0, tey_zz_y_1, tey_zz_yy_0, tey_zz_yy_1, tey_zz_yz_0, tey_zz_yz_1, \
                                         tey_zz_z_0, tey_zz_z_1, tey_zz_zz_0, tey_zz_zz_1, tey_zzz_xx_0, tey_zzz_xy_0, \
                                         tey_zzz_xz_0, tey_zzz_yy_0, tey_zzz_yz_0, tey_zzz_zz_0, tez_yyz_yy_0, tez_yyz_yz_0, \
                                         tez_yyz_zz_0, tez_yz_y_0, tez_yz_y_1, tez_yz_yy_0, tez_yz_yy_1, tez_yz_yz_0, \
                                         tez_yz_yz_1, tez_yz_z_0, tez_yz_z_1, tez_yz_zz_0, tez_yz_zz_1, tez_yzz_xx_0, \
                                         tez_yzz_xy_0, tez_yzz_xz_0, tez_yzz_yy_0, tez_yzz_yz_0, tez_yzz_zz_0, tez_z_xx_0, \
                                         tez_z_xx_1, tez_z_xy_0, tez_z_xy_1, tez_z_xz_0, tez_z_xz_1, tez_z_yy_0, tez_z_yy_1, \
                                         tez_z_yz_0, tez_z_yz_1, tez_z_zz_0, tez_z_zz_1, tez_zz_x_0, tez_zz_x_1, \
                                         tez_zz_xx_0, tez_zz_xx_1, tez_zz_xy_0, tez_zz_xy_1, tez_zz_xz_0, tez_zz_xz_1, \
                                         tez_zz_y_0, tez_zz_y_1, tez_zz_yy_0, tez_zz_yy_1, tez_zz_yz_0, tez_zz_yz_1, \
                                         tez_zz_z_0, tez_zz_z_1, tez_zz_zz_0, tez_zz_zz_1, tez_zzz_xx_0, tez_zzz_xy_0, \
                                         tez_zzz_xz_0, tez_zzz_yy_0, tez_zzz_yz_0, tez_zzz_zz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yyz_yy_0[j] = pa_y[j] * tex_yz_yy_0[j] - pc_y[j] * tex_yz_yy_1[j] + 0.5 * fl1_fx * tex_z_yy_0[j] - 0.5 * fl1_fx * tex_z_yy_1[j] + fl1_fx * tex_yz_y_0[j] - fl1_fx * tex_yz_y_1[j];

                    tey_yyz_yy_0[j] = pa_y[j] * tey_yz_yy_0[j] - pc_y[j] * tey_yz_yy_1[j] + 0.5 * fl1_fx * tey_z_yy_0[j] - 0.5 * fl1_fx * tey_z_yy_1[j] + fl1_fx * tey_yz_y_0[j] - fl1_fx * tey_yz_y_1[j] + ta_yz_yy_1[j];

                    tez_yyz_yy_0[j] = pa_y[j] * tez_yz_yy_0[j] - pc_y[j] * tez_yz_yy_1[j] + 0.5 * fl1_fx * tez_z_yy_0[j] - 0.5 * fl1_fx * tez_z_yy_1[j] + fl1_fx * tez_yz_y_0[j] - fl1_fx * tez_yz_y_1[j];

                    tex_yyz_yz_0[j] = pa_y[j] * tex_yz_yz_0[j] - pc_y[j] * tex_yz_yz_1[j] + 0.5 * fl1_fx * tex_z_yz_0[j] - 0.5 * fl1_fx * tex_z_yz_1[j] + 0.5 * fl1_fx * tex_yz_z_0[j] - 0.5 * fl1_fx * tex_yz_z_1[j];

                    tey_yyz_yz_0[j] = pa_y[j] * tey_yz_yz_0[j] - pc_y[j] * tey_yz_yz_1[j] + 0.5 * fl1_fx * tey_z_yz_0[j] - 0.5 * fl1_fx * tey_z_yz_1[j] + 0.5 * fl1_fx * tey_yz_z_0[j] - 0.5 * fl1_fx * tey_yz_z_1[j] + ta_yz_yz_1[j];

                    tez_yyz_yz_0[j] = pa_y[j] * tez_yz_yz_0[j] - pc_y[j] * tez_yz_yz_1[j] + 0.5 * fl1_fx * tez_z_yz_0[j] - 0.5 * fl1_fx * tez_z_yz_1[j] + 0.5 * fl1_fx * tez_yz_z_0[j] - 0.5 * fl1_fx * tez_yz_z_1[j];

                    tex_yyz_zz_0[j] = pa_y[j] * tex_yz_zz_0[j] - pc_y[j] * tex_yz_zz_1[j] + 0.5 * fl1_fx * tex_z_zz_0[j] - 0.5 * fl1_fx * tex_z_zz_1[j];

                    tey_yyz_zz_0[j] = pa_y[j] * tey_yz_zz_0[j] - pc_y[j] * tey_yz_zz_1[j] + 0.5 * fl1_fx * tey_z_zz_0[j] - 0.5 * fl1_fx * tey_z_zz_1[j] + ta_yz_zz_1[j];

                    tez_yyz_zz_0[j] = pa_y[j] * tez_yz_zz_0[j] - pc_y[j] * tez_yz_zz_1[j] + 0.5 * fl1_fx * tez_z_zz_0[j] - 0.5 * fl1_fx * tez_z_zz_1[j];

                    tex_yzz_xx_0[j] = pa_y[j] * tex_zz_xx_0[j] - pc_y[j] * tex_zz_xx_1[j];

                    tey_yzz_xx_0[j] = pa_y[j] * tey_zz_xx_0[j] - pc_y[j] * tey_zz_xx_1[j] + ta_zz_xx_1[j];

                    tez_yzz_xx_0[j] = pa_y[j] * tez_zz_xx_0[j] - pc_y[j] * tez_zz_xx_1[j];

                    tex_yzz_xy_0[j] = pa_y[j] * tex_zz_xy_0[j] - pc_y[j] * tex_zz_xy_1[j] + 0.5 * fl1_fx * tex_zz_x_0[j] - 0.5 * fl1_fx * tex_zz_x_1[j];

                    tey_yzz_xy_0[j] = pa_y[j] * tey_zz_xy_0[j] - pc_y[j] * tey_zz_xy_1[j] + 0.5 * fl1_fx * tey_zz_x_0[j] - 0.5 * fl1_fx * tey_zz_x_1[j] + ta_zz_xy_1[j];

                    tez_yzz_xy_0[j] = pa_y[j] * tez_zz_xy_0[j] - pc_y[j] * tez_zz_xy_1[j] + 0.5 * fl1_fx * tez_zz_x_0[j] - 0.5 * fl1_fx * tez_zz_x_1[j];

                    tex_yzz_xz_0[j] = pa_y[j] * tex_zz_xz_0[j] - pc_y[j] * tex_zz_xz_1[j];

                    tey_yzz_xz_0[j] = pa_y[j] * tey_zz_xz_0[j] - pc_y[j] * tey_zz_xz_1[j] + ta_zz_xz_1[j];

                    tez_yzz_xz_0[j] = pa_y[j] * tez_zz_xz_0[j] - pc_y[j] * tez_zz_xz_1[j];

                    tex_yzz_yy_0[j] = pa_y[j] * tex_zz_yy_0[j] - pc_y[j] * tex_zz_yy_1[j] + fl1_fx * tex_zz_y_0[j] - fl1_fx * tex_zz_y_1[j];

                    tey_yzz_yy_0[j] = pa_y[j] * tey_zz_yy_0[j] - pc_y[j] * tey_zz_yy_1[j] + fl1_fx * tey_zz_y_0[j] - fl1_fx * tey_zz_y_1[j] + ta_zz_yy_1[j];

                    tez_yzz_yy_0[j] = pa_y[j] * tez_zz_yy_0[j] - pc_y[j] * tez_zz_yy_1[j] + fl1_fx * tez_zz_y_0[j] - fl1_fx * tez_zz_y_1[j];

                    tex_yzz_yz_0[j] = pa_y[j] * tex_zz_yz_0[j] - pc_y[j] * tex_zz_yz_1[j] + 0.5 * fl1_fx * tex_zz_z_0[j] - 0.5 * fl1_fx * tex_zz_z_1[j];

                    tey_yzz_yz_0[j] = pa_y[j] * tey_zz_yz_0[j] - pc_y[j] * tey_zz_yz_1[j] + 0.5 * fl1_fx * tey_zz_z_0[j] - 0.5 * fl1_fx * tey_zz_z_1[j] + ta_zz_yz_1[j];

                    tez_yzz_yz_0[j] = pa_y[j] * tez_zz_yz_0[j] - pc_y[j] * tez_zz_yz_1[j] + 0.5 * fl1_fx * tez_zz_z_0[j] - 0.5 * fl1_fx * tez_zz_z_1[j];

                    tex_yzz_zz_0[j] = pa_y[j] * tex_zz_zz_0[j] - pc_y[j] * tex_zz_zz_1[j];

                    tey_yzz_zz_0[j] = pa_y[j] * tey_zz_zz_0[j] - pc_y[j] * tey_zz_zz_1[j] + ta_zz_zz_1[j];

                    tez_yzz_zz_0[j] = pa_y[j] * tez_zz_zz_0[j] - pc_y[j] * tez_zz_zz_1[j];

                    tex_zzz_xx_0[j] = pa_z[j] * tex_zz_xx_0[j] - pc_z[j] * tex_zz_xx_1[j] + fl1_fx * tex_z_xx_0[j] - fl1_fx * tex_z_xx_1[j];

                    tey_zzz_xx_0[j] = pa_z[j] * tey_zz_xx_0[j] - pc_z[j] * tey_zz_xx_1[j] + fl1_fx * tey_z_xx_0[j] - fl1_fx * tey_z_xx_1[j];

                    tez_zzz_xx_0[j] = pa_z[j] * tez_zz_xx_0[j] - pc_z[j] * tez_zz_xx_1[j] + fl1_fx * tez_z_xx_0[j] - fl1_fx * tez_z_xx_1[j] + ta_zz_xx_1[j];

                    tex_zzz_xy_0[j] = pa_z[j] * tex_zz_xy_0[j] - pc_z[j] * tex_zz_xy_1[j] + fl1_fx * tex_z_xy_0[j] - fl1_fx * tex_z_xy_1[j];

                    tey_zzz_xy_0[j] = pa_z[j] * tey_zz_xy_0[j] - pc_z[j] * tey_zz_xy_1[j] + fl1_fx * tey_z_xy_0[j] - fl1_fx * tey_z_xy_1[j];

                    tez_zzz_xy_0[j] = pa_z[j] * tez_zz_xy_0[j] - pc_z[j] * tez_zz_xy_1[j] + fl1_fx * tez_z_xy_0[j] - fl1_fx * tez_z_xy_1[j] + ta_zz_xy_1[j];

                    tex_zzz_xz_0[j] = pa_z[j] * tex_zz_xz_0[j] - pc_z[j] * tex_zz_xz_1[j] + fl1_fx * tex_z_xz_0[j] - fl1_fx * tex_z_xz_1[j] + 0.5 * fl1_fx * tex_zz_x_0[j] - 0.5 * fl1_fx * tex_zz_x_1[j];

                    tey_zzz_xz_0[j] = pa_z[j] * tey_zz_xz_0[j] - pc_z[j] * tey_zz_xz_1[j] + fl1_fx * tey_z_xz_0[j] - fl1_fx * tey_z_xz_1[j] + 0.5 * fl1_fx * tey_zz_x_0[j] - 0.5 * fl1_fx * tey_zz_x_1[j];

                    tez_zzz_xz_0[j] = pa_z[j] * tez_zz_xz_0[j] - pc_z[j] * tez_zz_xz_1[j] + fl1_fx * tez_z_xz_0[j] - fl1_fx * tez_z_xz_1[j] + 0.5 * fl1_fx * tez_zz_x_0[j] - 0.5 * fl1_fx * tez_zz_x_1[j] + ta_zz_xz_1[j];

                    tex_zzz_yy_0[j] = pa_z[j] * tex_zz_yy_0[j] - pc_z[j] * tex_zz_yy_1[j] + fl1_fx * tex_z_yy_0[j] - fl1_fx * tex_z_yy_1[j];

                    tey_zzz_yy_0[j] = pa_z[j] * tey_zz_yy_0[j] - pc_z[j] * tey_zz_yy_1[j] + fl1_fx * tey_z_yy_0[j] - fl1_fx * tey_z_yy_1[j];

                    tez_zzz_yy_0[j] = pa_z[j] * tez_zz_yy_0[j] - pc_z[j] * tez_zz_yy_1[j] + fl1_fx * tez_z_yy_0[j] - fl1_fx * tez_z_yy_1[j] + ta_zz_yy_1[j];

                    tex_zzz_yz_0[j] = pa_z[j] * tex_zz_yz_0[j] - pc_z[j] * tex_zz_yz_1[j] + fl1_fx * tex_z_yz_0[j] - fl1_fx * tex_z_yz_1[j] + 0.5 * fl1_fx * tex_zz_y_0[j] - 0.5 * fl1_fx * tex_zz_y_1[j];

                    tey_zzz_yz_0[j] = pa_z[j] * tey_zz_yz_0[j] - pc_z[j] * tey_zz_yz_1[j] + fl1_fx * tey_z_yz_0[j] - fl1_fx * tey_z_yz_1[j] + 0.5 * fl1_fx * tey_zz_y_0[j] - 0.5 * fl1_fx * tey_zz_y_1[j];

                    tez_zzz_yz_0[j] = pa_z[j] * tez_zz_yz_0[j] - pc_z[j] * tez_zz_yz_1[j] + fl1_fx * tez_z_yz_0[j] - fl1_fx * tez_z_yz_1[j] + 0.5 * fl1_fx * tez_zz_y_0[j] - 0.5 * fl1_fx * tez_zz_y_1[j] + ta_zz_yz_1[j];

                    tex_zzz_zz_0[j] = pa_z[j] * tex_zz_zz_0[j] - pc_z[j] * tex_zz_zz_1[j] + fl1_fx * tex_z_zz_0[j] - fl1_fx * tex_z_zz_1[j] + fl1_fx * tex_zz_z_0[j] - fl1_fx * tex_zz_z_1[j];

                    tey_zzz_zz_0[j] = pa_z[j] * tey_zz_zz_0[j] - pc_z[j] * tey_zz_zz_1[j] + fl1_fx * tey_z_zz_0[j] - fl1_fx * tey_z_zz_1[j] + fl1_fx * tey_zz_z_0[j] - fl1_fx * tey_zz_z_1[j];

                    tez_zzz_zz_0[j] = pa_z[j] * tez_zz_zz_0[j] - pc_z[j] * tez_zz_zz_1[j] + fl1_fx * tez_z_zz_0[j] - fl1_fx * tez_z_zz_1[j] + fl1_fx * tez_zz_z_0[j] - fl1_fx * tez_zz_z_1[j] + ta_zz_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDG(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pcDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        efieldrecfunc::compElectricFieldForDG_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   pcDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        efieldrecfunc::compElectricFieldForDG_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    pcDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        efieldrecfunc::compElectricFieldForDG_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     pcDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        efieldrecfunc::compElectricFieldForDG_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForDG_180_225(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForDG_225_270(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compElectricFieldForDG_0_45(      CMemBlock2D<double>& primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
    {
        // Batch of Integrals (0,45)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
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

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_x = pcDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tex_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx); 

                auto tey_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx); 

                auto tez_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx); 

                auto tex_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 1); 

                auto tey_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 1); 

                auto tez_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 1); 

                auto tex_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 2); 

                auto tey_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 2); 

                auto tez_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 2); 

                auto tex_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 3); 

                auto tey_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 3); 

                auto tez_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 3); 

                auto tex_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 4); 

                auto tey_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 4); 

                auto tez_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 4); 

                auto tex_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 5); 

                auto tey_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 5); 

                auto tez_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 5); 

                auto tex_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 6); 

                auto tey_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 6); 

                auto tez_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 6); 

                auto tex_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 7); 

                auto tey_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 7); 

                auto tez_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 7); 

                auto tex_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 8); 

                auto tey_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 8); 

                auto tez_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 8); 

                auto tex_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 9); 

                auto tey_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 9); 

                auto tez_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 9); 

                auto tex_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 10); 

                auto tey_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 10); 

                auto tez_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 10); 

                auto tex_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 11); 

                auto tey_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 11); 

                auto tez_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 11); 

                auto tex_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 12); 

                auto tey_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 12); 

                auto tez_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 12); 

                auto tex_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 13); 

                auto tey_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 13); 

                auto tez_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 13); 

                auto tex_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 14); 

                auto tey_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 14); 

                auto tez_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 14); 

                auto tex_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx); 

                auto tey_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx); 

                auto tez_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx); 

                auto tex_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 1); 

                auto tey_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 1); 

                auto tez_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 1); 

                auto tex_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 2); 

                auto tey_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 2); 

                auto tez_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 2); 

                auto tex_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 3); 

                auto tey_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 3); 

                auto tez_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 3); 

                auto tex_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 4); 

                auto tey_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 4); 

                auto tez_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 4); 

                auto tex_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 5); 

                auto tey_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 5); 

                auto tez_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 5); 

                auto tex_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 6); 

                auto tey_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 6); 

                auto tez_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 6); 

                auto tex_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 7); 

                auto tey_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 7); 

                auto tez_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 7); 

                auto tex_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 8); 

                auto tey_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 8); 

                auto tez_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 8); 

                auto tex_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 9); 

                auto tey_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 9); 

                auto tez_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 9); 

                auto tex_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 10); 

                auto tey_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 10); 

                auto tez_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 10); 

                auto tex_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 11); 

                auto tey_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 11); 

                auto tez_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 11); 

                auto tex_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 12); 

                auto tey_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 12); 

                auto tez_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 12); 

                auto tex_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 13); 

                auto tey_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 13); 

                auto tez_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 13); 

                auto tex_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 14); 

                auto tey_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 14); 

                auto tez_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 14); 

                auto tex_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx); 

                auto tey_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx); 

                auto tez_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx); 

                auto tex_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 1); 

                auto tey_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 1); 

                auto tez_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 1); 

                auto tex_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 2); 

                auto tey_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 2); 

                auto tez_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 2); 

                auto tex_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 3); 

                auto tey_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 3); 

                auto tez_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 3); 

                auto tex_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 4); 

                auto tey_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 4); 

                auto tez_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 4); 

                auto tex_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 5); 

                auto tey_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 5); 

                auto tez_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 5); 

                auto tex_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 6); 

                auto tey_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 6); 

                auto tez_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 6); 

                auto tex_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 7); 

                auto tey_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 7); 

                auto tez_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 7); 

                auto tex_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 8); 

                auto tey_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 8); 

                auto tez_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 8); 

                auto tex_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 9); 

                auto tey_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 9); 

                auto tez_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 9); 

                auto tex_x_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx); 

                auto tey_x_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx); 

                auto tez_x_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx); 

                auto tex_x_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 1); 

                auto tey_x_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 1); 

                auto tez_x_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 1); 

                auto tex_x_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 2); 

                auto tey_x_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 2); 

                auto tez_x_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 2); 

                auto tex_x_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 3); 

                auto tey_x_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 3); 

                auto tez_x_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 3); 

                auto tex_x_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 4); 

                auto tey_x_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 4); 

                auto tez_x_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 4); 

                auto tex_x_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 5); 

                auto tey_x_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 5); 

                auto tez_x_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 5); 

                auto tex_x_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 6); 

                auto tey_x_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 6); 

                auto tez_x_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 6); 

                auto tex_x_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 7); 

                auto tey_x_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 7); 

                auto tez_x_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 7); 

                auto tex_x_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 8); 

                auto tey_x_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 8); 

                auto tez_x_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 8); 

                auto tex_x_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 9); 

                auto tey_x_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 9); 

                auto tez_x_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 9); 

                auto ta_x_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx); 

                auto ta_x_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 1); 

                auto ta_x_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 2); 

                auto ta_x_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 3); 

                auto ta_x_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 4); 

                auto ta_x_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 5); 

                auto ta_x_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 6); 

                auto ta_x_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 7); 

                auto ta_x_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 8); 

                auto ta_x_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 9); 

                auto ta_x_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 10); 

                auto ta_x_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 11); 

                auto ta_x_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 12); 

                auto ta_x_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 13); 

                auto ta_x_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 14); 

                // set up pointers to integrals

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

                // Batch of Integrals (0,45)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_x_xxxx_1, ta_x_xxxy_1, ta_x_xxxz_1, ta_x_xxyy_1, \
                                         ta_x_xxyz_1, ta_x_xxzz_1, ta_x_xyyy_1, ta_x_xyyz_1, ta_x_xyzz_1, ta_x_xzzz_1, \
                                         ta_x_yyyy_1, ta_x_yyyz_1, ta_x_yyzz_1, ta_x_yzzz_1, ta_x_zzzz_1, tex_0_xxxx_0, \
                                         tex_0_xxxx_1, tex_0_xxxy_0, tex_0_xxxy_1, tex_0_xxxz_0, tex_0_xxxz_1, tex_0_xxyy_0, \
                                         tex_0_xxyy_1, tex_0_xxyz_0, tex_0_xxyz_1, tex_0_xxzz_0, tex_0_xxzz_1, tex_0_xyyy_0, \
                                         tex_0_xyyy_1, tex_0_xyyz_0, tex_0_xyyz_1, tex_0_xyzz_0, tex_0_xyzz_1, tex_0_xzzz_0, \
                                         tex_0_xzzz_1, tex_0_yyyy_0, tex_0_yyyy_1, tex_0_yyyz_0, tex_0_yyyz_1, tex_0_yyzz_0, \
                                         tex_0_yyzz_1, tex_0_yzzz_0, tex_0_yzzz_1, tex_0_zzzz_0, tex_0_zzzz_1, tex_x_xxx_0, \
                                         tex_x_xxx_1, tex_x_xxxx_0, tex_x_xxxx_1, tex_x_xxxy_0, tex_x_xxxy_1, tex_x_xxxz_0, \
                                         tex_x_xxxz_1, tex_x_xxy_0, tex_x_xxy_1, tex_x_xxyy_0, tex_x_xxyy_1, tex_x_xxyz_0, \
                                         tex_x_xxyz_1, tex_x_xxz_0, tex_x_xxz_1, tex_x_xxzz_0, tex_x_xxzz_1, tex_x_xyy_0, \
                                         tex_x_xyy_1, tex_x_xyyy_0, tex_x_xyyy_1, tex_x_xyyz_0, tex_x_xyyz_1, tex_x_xyz_0, \
                                         tex_x_xyz_1, tex_x_xyzz_0, tex_x_xyzz_1, tex_x_xzz_0, tex_x_xzz_1, tex_x_xzzz_0, \
                                         tex_x_xzzz_1, tex_x_yyy_0, tex_x_yyy_1, tex_x_yyyy_0, tex_x_yyyy_1, tex_x_yyyz_0, \
                                         tex_x_yyyz_1, tex_x_yyz_0, tex_x_yyz_1, tex_x_yyzz_0, tex_x_yyzz_1, tex_x_yzz_0, \
                                         tex_x_yzz_1, tex_x_yzzz_0, tex_x_yzzz_1, tex_x_zzz_0, tex_x_zzz_1, tex_x_zzzz_0, \
                                         tex_x_zzzz_1, tex_xx_xxxx_0, tex_xx_xxxy_0, tex_xx_xxxz_0, tex_xx_xxyy_0, \
                                         tex_xx_xxyz_0, tex_xx_xxzz_0, tex_xx_xyyy_0, tex_xx_xyyz_0, tex_xx_xyzz_0, \
                                         tex_xx_xzzz_0, tex_xx_yyyy_0, tex_xx_yyyz_0, tex_xx_yyzz_0, tex_xx_yzzz_0, \
                                         tex_xx_zzzz_0, tey_0_xxxx_0, tey_0_xxxx_1, tey_0_xxxy_0, tey_0_xxxy_1, tey_0_xxxz_0, \
                                         tey_0_xxxz_1, tey_0_xxyy_0, tey_0_xxyy_1, tey_0_xxyz_0, tey_0_xxyz_1, tey_0_xxzz_0, \
                                         tey_0_xxzz_1, tey_0_xyyy_0, tey_0_xyyy_1, tey_0_xyyz_0, tey_0_xyyz_1, tey_0_xyzz_0, \
                                         tey_0_xyzz_1, tey_0_xzzz_0, tey_0_xzzz_1, tey_0_yyyy_0, tey_0_yyyy_1, tey_0_yyyz_0, \
                                         tey_0_yyyz_1, tey_0_yyzz_0, tey_0_yyzz_1, tey_0_yzzz_0, tey_0_yzzz_1, tey_0_zzzz_0, \
                                         tey_0_zzzz_1, tey_x_xxx_0, tey_x_xxx_1, tey_x_xxxx_0, tey_x_xxxx_1, tey_x_xxxy_0, \
                                         tey_x_xxxy_1, tey_x_xxxz_0, tey_x_xxxz_1, tey_x_xxy_0, tey_x_xxy_1, tey_x_xxyy_0, \
                                         tey_x_xxyy_1, tey_x_xxyz_0, tey_x_xxyz_1, tey_x_xxz_0, tey_x_xxz_1, tey_x_xxzz_0, \
                                         tey_x_xxzz_1, tey_x_xyy_0, tey_x_xyy_1, tey_x_xyyy_0, tey_x_xyyy_1, tey_x_xyyz_0, \
                                         tey_x_xyyz_1, tey_x_xyz_0, tey_x_xyz_1, tey_x_xyzz_0, tey_x_xyzz_1, tey_x_xzz_0, \
                                         tey_x_xzz_1, tey_x_xzzz_0, tey_x_xzzz_1, tey_x_yyy_0, tey_x_yyy_1, tey_x_yyyy_0, \
                                         tey_x_yyyy_1, tey_x_yyyz_0, tey_x_yyyz_1, tey_x_yyz_0, tey_x_yyz_1, tey_x_yyzz_0, \
                                         tey_x_yyzz_1, tey_x_yzz_0, tey_x_yzz_1, tey_x_yzzz_0, tey_x_yzzz_1, tey_x_zzz_0, \
                                         tey_x_zzz_1, tey_x_zzzz_0, tey_x_zzzz_1, tey_xx_xxxx_0, tey_xx_xxxy_0, \
                                         tey_xx_xxxz_0, tey_xx_xxyy_0, tey_xx_xxyz_0, tey_xx_xxzz_0, tey_xx_xyyy_0, \
                                         tey_xx_xyyz_0, tey_xx_xyzz_0, tey_xx_xzzz_0, tey_xx_yyyy_0, tey_xx_yyyz_0, \
                                         tey_xx_yyzz_0, tey_xx_yzzz_0, tey_xx_zzzz_0, tez_0_xxxx_0, tez_0_xxxx_1, \
                                         tez_0_xxxy_0, tez_0_xxxy_1, tez_0_xxxz_0, tez_0_xxxz_1, tez_0_xxyy_0, tez_0_xxyy_1, \
                                         tez_0_xxyz_0, tez_0_xxyz_1, tez_0_xxzz_0, tez_0_xxzz_1, tez_0_xyyy_0, tez_0_xyyy_1, \
                                         tez_0_xyyz_0, tez_0_xyyz_1, tez_0_xyzz_0, tez_0_xyzz_1, tez_0_xzzz_0, tez_0_xzzz_1, \
                                         tez_0_yyyy_0, tez_0_yyyy_1, tez_0_yyyz_0, tez_0_yyyz_1, tez_0_yyzz_0, tez_0_yyzz_1, \
                                         tez_0_yzzz_0, tez_0_yzzz_1, tez_0_zzzz_0, tez_0_zzzz_1, tez_x_xxx_0, tez_x_xxx_1, \
                                         tez_x_xxxx_0, tez_x_xxxx_1, tez_x_xxxy_0, tez_x_xxxy_1, tez_x_xxxz_0, tez_x_xxxz_1, \
                                         tez_x_xxy_0, tez_x_xxy_1, tez_x_xxyy_0, tez_x_xxyy_1, tez_x_xxyz_0, tez_x_xxyz_1, \
                                         tez_x_xxz_0, tez_x_xxz_1, tez_x_xxzz_0, tez_x_xxzz_1, tez_x_xyy_0, tez_x_xyy_1, \
                                         tez_x_xyyy_0, tez_x_xyyy_1, tez_x_xyyz_0, tez_x_xyyz_1, tez_x_xyz_0, tez_x_xyz_1, \
                                         tez_x_xyzz_0, tez_x_xyzz_1, tez_x_xzz_0, tez_x_xzz_1, tez_x_xzzz_0, tez_x_xzzz_1, \
                                         tez_x_yyy_0, tez_x_yyy_1, tez_x_yyyy_0, tez_x_yyyy_1, tez_x_yyyz_0, tez_x_yyyz_1, \
                                         tez_x_yyz_0, tez_x_yyz_1, tez_x_yyzz_0, tez_x_yyzz_1, tez_x_yzz_0, tez_x_yzz_1, \
                                         tez_x_yzzz_0, tez_x_yzzz_1, tez_x_zzz_0, tez_x_zzz_1, tez_x_zzzz_0, tez_x_zzzz_1, \
                                         tez_xx_xxxx_0, tez_xx_xxxy_0, tez_xx_xxxz_0, tez_xx_xxyy_0, tez_xx_xxyz_0, \
                                         tez_xx_xxzz_0, tez_xx_xyyy_0, tez_xx_xyyz_0, tez_xx_xyzz_0, tez_xx_xzzz_0, \
                                         tez_xx_yyyy_0, tez_xx_yyyz_0, tez_xx_yyzz_0, tez_xx_yzzz_0, tez_xx_zzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xx_xxxx_0[j] = pa_x[j] * tex_x_xxxx_0[j] - pc_x[j] * tex_x_xxxx_1[j] + 0.5 * fl1_fx * tex_0_xxxx_0[j] - 0.5 * fl1_fx * tex_0_xxxx_1[j] + 2.0 * fl1_fx * tex_x_xxx_0[j] - 2.0 * fl1_fx * tex_x_xxx_1[j] + ta_x_xxxx_1[j];

                    tey_xx_xxxx_0[j] = pa_x[j] * tey_x_xxxx_0[j] - pc_x[j] * tey_x_xxxx_1[j] + 0.5 * fl1_fx * tey_0_xxxx_0[j] - 0.5 * fl1_fx * tey_0_xxxx_1[j] + 2.0 * fl1_fx * tey_x_xxx_0[j] - 2.0 * fl1_fx * tey_x_xxx_1[j];

                    tez_xx_xxxx_0[j] = pa_x[j] * tez_x_xxxx_0[j] - pc_x[j] * tez_x_xxxx_1[j] + 0.5 * fl1_fx * tez_0_xxxx_0[j] - 0.5 * fl1_fx * tez_0_xxxx_1[j] + 2.0 * fl1_fx * tez_x_xxx_0[j] - 2.0 * fl1_fx * tez_x_xxx_1[j];

                    tex_xx_xxxy_0[j] = pa_x[j] * tex_x_xxxy_0[j] - pc_x[j] * tex_x_xxxy_1[j] + 0.5 * fl1_fx * tex_0_xxxy_0[j] - 0.5 * fl1_fx * tex_0_xxxy_1[j] + 1.5 * fl1_fx * tex_x_xxy_0[j] - 1.5 * fl1_fx * tex_x_xxy_1[j] + ta_x_xxxy_1[j];

                    tey_xx_xxxy_0[j] = pa_x[j] * tey_x_xxxy_0[j] - pc_x[j] * tey_x_xxxy_1[j] + 0.5 * fl1_fx * tey_0_xxxy_0[j] - 0.5 * fl1_fx * tey_0_xxxy_1[j] + 1.5 * fl1_fx * tey_x_xxy_0[j] - 1.5 * fl1_fx * tey_x_xxy_1[j];

                    tez_xx_xxxy_0[j] = pa_x[j] * tez_x_xxxy_0[j] - pc_x[j] * tez_x_xxxy_1[j] + 0.5 * fl1_fx * tez_0_xxxy_0[j] - 0.5 * fl1_fx * tez_0_xxxy_1[j] + 1.5 * fl1_fx * tez_x_xxy_0[j] - 1.5 * fl1_fx * tez_x_xxy_1[j];

                    tex_xx_xxxz_0[j] = pa_x[j] * tex_x_xxxz_0[j] - pc_x[j] * tex_x_xxxz_1[j] + 0.5 * fl1_fx * tex_0_xxxz_0[j] - 0.5 * fl1_fx * tex_0_xxxz_1[j] + 1.5 * fl1_fx * tex_x_xxz_0[j] - 1.5 * fl1_fx * tex_x_xxz_1[j] + ta_x_xxxz_1[j];

                    tey_xx_xxxz_0[j] = pa_x[j] * tey_x_xxxz_0[j] - pc_x[j] * tey_x_xxxz_1[j] + 0.5 * fl1_fx * tey_0_xxxz_0[j] - 0.5 * fl1_fx * tey_0_xxxz_1[j] + 1.5 * fl1_fx * tey_x_xxz_0[j] - 1.5 * fl1_fx * tey_x_xxz_1[j];

                    tez_xx_xxxz_0[j] = pa_x[j] * tez_x_xxxz_0[j] - pc_x[j] * tez_x_xxxz_1[j] + 0.5 * fl1_fx * tez_0_xxxz_0[j] - 0.5 * fl1_fx * tez_0_xxxz_1[j] + 1.5 * fl1_fx * tez_x_xxz_0[j] - 1.5 * fl1_fx * tez_x_xxz_1[j];

                    tex_xx_xxyy_0[j] = pa_x[j] * tex_x_xxyy_0[j] - pc_x[j] * tex_x_xxyy_1[j] + 0.5 * fl1_fx * tex_0_xxyy_0[j] - 0.5 * fl1_fx * tex_0_xxyy_1[j] + fl1_fx * tex_x_xyy_0[j] - fl1_fx * tex_x_xyy_1[j] + ta_x_xxyy_1[j];

                    tey_xx_xxyy_0[j] = pa_x[j] * tey_x_xxyy_0[j] - pc_x[j] * tey_x_xxyy_1[j] + 0.5 * fl1_fx * tey_0_xxyy_0[j] - 0.5 * fl1_fx * tey_0_xxyy_1[j] + fl1_fx * tey_x_xyy_0[j] - fl1_fx * tey_x_xyy_1[j];

                    tez_xx_xxyy_0[j] = pa_x[j] * tez_x_xxyy_0[j] - pc_x[j] * tez_x_xxyy_1[j] + 0.5 * fl1_fx * tez_0_xxyy_0[j] - 0.5 * fl1_fx * tez_0_xxyy_1[j] + fl1_fx * tez_x_xyy_0[j] - fl1_fx * tez_x_xyy_1[j];

                    tex_xx_xxyz_0[j] = pa_x[j] * tex_x_xxyz_0[j] - pc_x[j] * tex_x_xxyz_1[j] + 0.5 * fl1_fx * tex_0_xxyz_0[j] - 0.5 * fl1_fx * tex_0_xxyz_1[j] + fl1_fx * tex_x_xyz_0[j] - fl1_fx * tex_x_xyz_1[j] + ta_x_xxyz_1[j];

                    tey_xx_xxyz_0[j] = pa_x[j] * tey_x_xxyz_0[j] - pc_x[j] * tey_x_xxyz_1[j] + 0.5 * fl1_fx * tey_0_xxyz_0[j] - 0.5 * fl1_fx * tey_0_xxyz_1[j] + fl1_fx * tey_x_xyz_0[j] - fl1_fx * tey_x_xyz_1[j];

                    tez_xx_xxyz_0[j] = pa_x[j] * tez_x_xxyz_0[j] - pc_x[j] * tez_x_xxyz_1[j] + 0.5 * fl1_fx * tez_0_xxyz_0[j] - 0.5 * fl1_fx * tez_0_xxyz_1[j] + fl1_fx * tez_x_xyz_0[j] - fl1_fx * tez_x_xyz_1[j];

                    tex_xx_xxzz_0[j] = pa_x[j] * tex_x_xxzz_0[j] - pc_x[j] * tex_x_xxzz_1[j] + 0.5 * fl1_fx * tex_0_xxzz_0[j] - 0.5 * fl1_fx * tex_0_xxzz_1[j] + fl1_fx * tex_x_xzz_0[j] - fl1_fx * tex_x_xzz_1[j] + ta_x_xxzz_1[j];

                    tey_xx_xxzz_0[j] = pa_x[j] * tey_x_xxzz_0[j] - pc_x[j] * tey_x_xxzz_1[j] + 0.5 * fl1_fx * tey_0_xxzz_0[j] - 0.5 * fl1_fx * tey_0_xxzz_1[j] + fl1_fx * tey_x_xzz_0[j] - fl1_fx * tey_x_xzz_1[j];

                    tez_xx_xxzz_0[j] = pa_x[j] * tez_x_xxzz_0[j] - pc_x[j] * tez_x_xxzz_1[j] + 0.5 * fl1_fx * tez_0_xxzz_0[j] - 0.5 * fl1_fx * tez_0_xxzz_1[j] + fl1_fx * tez_x_xzz_0[j] - fl1_fx * tez_x_xzz_1[j];

                    tex_xx_xyyy_0[j] = pa_x[j] * tex_x_xyyy_0[j] - pc_x[j] * tex_x_xyyy_1[j] + 0.5 * fl1_fx * tex_0_xyyy_0[j] - 0.5 * fl1_fx * tex_0_xyyy_1[j] + 0.5 * fl1_fx * tex_x_yyy_0[j] - 0.5 * fl1_fx * tex_x_yyy_1[j] + ta_x_xyyy_1[j];

                    tey_xx_xyyy_0[j] = pa_x[j] * tey_x_xyyy_0[j] - pc_x[j] * tey_x_xyyy_1[j] + 0.5 * fl1_fx * tey_0_xyyy_0[j] - 0.5 * fl1_fx * tey_0_xyyy_1[j] + 0.5 * fl1_fx * tey_x_yyy_0[j] - 0.5 * fl1_fx * tey_x_yyy_1[j];

                    tez_xx_xyyy_0[j] = pa_x[j] * tez_x_xyyy_0[j] - pc_x[j] * tez_x_xyyy_1[j] + 0.5 * fl1_fx * tez_0_xyyy_0[j] - 0.5 * fl1_fx * tez_0_xyyy_1[j] + 0.5 * fl1_fx * tez_x_yyy_0[j] - 0.5 * fl1_fx * tez_x_yyy_1[j];

                    tex_xx_xyyz_0[j] = pa_x[j] * tex_x_xyyz_0[j] - pc_x[j] * tex_x_xyyz_1[j] + 0.5 * fl1_fx * tex_0_xyyz_0[j] - 0.5 * fl1_fx * tex_0_xyyz_1[j] + 0.5 * fl1_fx * tex_x_yyz_0[j] - 0.5 * fl1_fx * tex_x_yyz_1[j] + ta_x_xyyz_1[j];

                    tey_xx_xyyz_0[j] = pa_x[j] * tey_x_xyyz_0[j] - pc_x[j] * tey_x_xyyz_1[j] + 0.5 * fl1_fx * tey_0_xyyz_0[j] - 0.5 * fl1_fx * tey_0_xyyz_1[j] + 0.5 * fl1_fx * tey_x_yyz_0[j] - 0.5 * fl1_fx * tey_x_yyz_1[j];

                    tez_xx_xyyz_0[j] = pa_x[j] * tez_x_xyyz_0[j] - pc_x[j] * tez_x_xyyz_1[j] + 0.5 * fl1_fx * tez_0_xyyz_0[j] - 0.5 * fl1_fx * tez_0_xyyz_1[j] + 0.5 * fl1_fx * tez_x_yyz_0[j] - 0.5 * fl1_fx * tez_x_yyz_1[j];

                    tex_xx_xyzz_0[j] = pa_x[j] * tex_x_xyzz_0[j] - pc_x[j] * tex_x_xyzz_1[j] + 0.5 * fl1_fx * tex_0_xyzz_0[j] - 0.5 * fl1_fx * tex_0_xyzz_1[j] + 0.5 * fl1_fx * tex_x_yzz_0[j] - 0.5 * fl1_fx * tex_x_yzz_1[j] + ta_x_xyzz_1[j];

                    tey_xx_xyzz_0[j] = pa_x[j] * tey_x_xyzz_0[j] - pc_x[j] * tey_x_xyzz_1[j] + 0.5 * fl1_fx * tey_0_xyzz_0[j] - 0.5 * fl1_fx * tey_0_xyzz_1[j] + 0.5 * fl1_fx * tey_x_yzz_0[j] - 0.5 * fl1_fx * tey_x_yzz_1[j];

                    tez_xx_xyzz_0[j] = pa_x[j] * tez_x_xyzz_0[j] - pc_x[j] * tez_x_xyzz_1[j] + 0.5 * fl1_fx * tez_0_xyzz_0[j] - 0.5 * fl1_fx * tez_0_xyzz_1[j] + 0.5 * fl1_fx * tez_x_yzz_0[j] - 0.5 * fl1_fx * tez_x_yzz_1[j];

                    tex_xx_xzzz_0[j] = pa_x[j] * tex_x_xzzz_0[j] - pc_x[j] * tex_x_xzzz_1[j] + 0.5 * fl1_fx * tex_0_xzzz_0[j] - 0.5 * fl1_fx * tex_0_xzzz_1[j] + 0.5 * fl1_fx * tex_x_zzz_0[j] - 0.5 * fl1_fx * tex_x_zzz_1[j] + ta_x_xzzz_1[j];

                    tey_xx_xzzz_0[j] = pa_x[j] * tey_x_xzzz_0[j] - pc_x[j] * tey_x_xzzz_1[j] + 0.5 * fl1_fx * tey_0_xzzz_0[j] - 0.5 * fl1_fx * tey_0_xzzz_1[j] + 0.5 * fl1_fx * tey_x_zzz_0[j] - 0.5 * fl1_fx * tey_x_zzz_1[j];

                    tez_xx_xzzz_0[j] = pa_x[j] * tez_x_xzzz_0[j] - pc_x[j] * tez_x_xzzz_1[j] + 0.5 * fl1_fx * tez_0_xzzz_0[j] - 0.5 * fl1_fx * tez_0_xzzz_1[j] + 0.5 * fl1_fx * tez_x_zzz_0[j] - 0.5 * fl1_fx * tez_x_zzz_1[j];

                    tex_xx_yyyy_0[j] = pa_x[j] * tex_x_yyyy_0[j] - pc_x[j] * tex_x_yyyy_1[j] + 0.5 * fl1_fx * tex_0_yyyy_0[j] - 0.5 * fl1_fx * tex_0_yyyy_1[j] + ta_x_yyyy_1[j];

                    tey_xx_yyyy_0[j] = pa_x[j] * tey_x_yyyy_0[j] - pc_x[j] * tey_x_yyyy_1[j] + 0.5 * fl1_fx * tey_0_yyyy_0[j] - 0.5 * fl1_fx * tey_0_yyyy_1[j];

                    tez_xx_yyyy_0[j] = pa_x[j] * tez_x_yyyy_0[j] - pc_x[j] * tez_x_yyyy_1[j] + 0.5 * fl1_fx * tez_0_yyyy_0[j] - 0.5 * fl1_fx * tez_0_yyyy_1[j];

                    tex_xx_yyyz_0[j] = pa_x[j] * tex_x_yyyz_0[j] - pc_x[j] * tex_x_yyyz_1[j] + 0.5 * fl1_fx * tex_0_yyyz_0[j] - 0.5 * fl1_fx * tex_0_yyyz_1[j] + ta_x_yyyz_1[j];

                    tey_xx_yyyz_0[j] = pa_x[j] * tey_x_yyyz_0[j] - pc_x[j] * tey_x_yyyz_1[j] + 0.5 * fl1_fx * tey_0_yyyz_0[j] - 0.5 * fl1_fx * tey_0_yyyz_1[j];

                    tez_xx_yyyz_0[j] = pa_x[j] * tez_x_yyyz_0[j] - pc_x[j] * tez_x_yyyz_1[j] + 0.5 * fl1_fx * tez_0_yyyz_0[j] - 0.5 * fl1_fx * tez_0_yyyz_1[j];

                    tex_xx_yyzz_0[j] = pa_x[j] * tex_x_yyzz_0[j] - pc_x[j] * tex_x_yyzz_1[j] + 0.5 * fl1_fx * tex_0_yyzz_0[j] - 0.5 * fl1_fx * tex_0_yyzz_1[j] + ta_x_yyzz_1[j];

                    tey_xx_yyzz_0[j] = pa_x[j] * tey_x_yyzz_0[j] - pc_x[j] * tey_x_yyzz_1[j] + 0.5 * fl1_fx * tey_0_yyzz_0[j] - 0.5 * fl1_fx * tey_0_yyzz_1[j];

                    tez_xx_yyzz_0[j] = pa_x[j] * tez_x_yyzz_0[j] - pc_x[j] * tez_x_yyzz_1[j] + 0.5 * fl1_fx * tez_0_yyzz_0[j] - 0.5 * fl1_fx * tez_0_yyzz_1[j];

                    tex_xx_yzzz_0[j] = pa_x[j] * tex_x_yzzz_0[j] - pc_x[j] * tex_x_yzzz_1[j] + 0.5 * fl1_fx * tex_0_yzzz_0[j] - 0.5 * fl1_fx * tex_0_yzzz_1[j] + ta_x_yzzz_1[j];

                    tey_xx_yzzz_0[j] = pa_x[j] * tey_x_yzzz_0[j] - pc_x[j] * tey_x_yzzz_1[j] + 0.5 * fl1_fx * tey_0_yzzz_0[j] - 0.5 * fl1_fx * tey_0_yzzz_1[j];

                    tez_xx_yzzz_0[j] = pa_x[j] * tez_x_yzzz_0[j] - pc_x[j] * tez_x_yzzz_1[j] + 0.5 * fl1_fx * tez_0_yzzz_0[j] - 0.5 * fl1_fx * tez_0_yzzz_1[j];

                    tex_xx_zzzz_0[j] = pa_x[j] * tex_x_zzzz_0[j] - pc_x[j] * tex_x_zzzz_1[j] + 0.5 * fl1_fx * tex_0_zzzz_0[j] - 0.5 * fl1_fx * tex_0_zzzz_1[j] + ta_x_zzzz_1[j];

                    tey_xx_zzzz_0[j] = pa_x[j] * tey_x_zzzz_0[j] - pc_x[j] * tey_x_zzzz_1[j] + 0.5 * fl1_fx * tey_0_zzzz_0[j] - 0.5 * fl1_fx * tey_0_zzzz_1[j];

                    tez_xx_zzzz_0[j] = pa_x[j] * tez_x_zzzz_0[j] - pc_x[j] * tez_x_zzzz_1[j] + 0.5 * fl1_fx * tez_0_zzzz_0[j] - 0.5 * fl1_fx * tez_0_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDG_45_90(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (45,90)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
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

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_x = pcDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tex_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 10); 

                auto tey_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 11); 

                auto tey_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 12); 

                auto tey_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 13); 

                auto tey_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 14); 

                auto tey_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 14); 

                auto tex_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 15); 

                auto tey_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 16); 

                auto tey_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 17); 

                auto tey_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 18); 

                auto tey_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 19); 

                auto tey_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 19); 

                auto tex_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 10); 

                auto tey_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 11); 

                auto tey_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 12); 

                auto tey_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 13); 

                auto tey_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 14); 

                auto tey_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 14); 

                auto tex_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 15); 

                auto tey_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 16); 

                auto tey_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 17); 

                auto tey_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 18); 

                auto tey_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 19); 

                auto tey_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 19); 

                auto ta_y_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 15); 

                auto ta_y_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 16); 

                auto ta_y_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 17); 

                auto ta_y_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 18); 

                auto ta_y_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 19); 

                auto ta_y_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 20); 

                auto ta_y_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 21); 

                auto ta_y_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 22); 

                auto ta_y_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 23); 

                auto ta_y_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 24); 

                auto ta_y_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 25); 

                auto ta_y_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 26); 

                auto ta_y_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 27); 

                auto ta_y_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 28); 

                auto ta_y_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 29); 

                // set up pointers to integrals

                auto tex_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 15); 

                auto tey_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 15); 

                auto tez_xy_xxxx_0 = primBuffer.data(pidx_e_2_4_m0 + 180 * bdim + 90 * idx + 15); 

                auto tex_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * idx + 16); 

                auto tey_xy_xxxy_0 = primBuffer.data(pidx_e_2_4_m0 + 90 * bdim + 90 * idx + 16); 

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

                // Batch of Integrals (45,90)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_y_xxxx_1, ta_y_xxxy_1, ta_y_xxxz_1, ta_y_xxyy_1, \
                                         ta_y_xxyz_1, ta_y_xxzz_1, ta_y_xyyy_1, ta_y_xyyz_1, ta_y_xyzz_1, ta_y_xzzz_1, \
                                         ta_y_yyyy_1, ta_y_yyyz_1, ta_y_yyzz_1, ta_y_yzzz_1, ta_y_zzzz_1, tex_xy_xxxx_0, \
                                         tex_xy_xxxy_0, tex_xy_xxxz_0, tex_xy_xxyy_0, tex_xy_xxyz_0, tex_xy_xxzz_0, \
                                         tex_xy_xyyy_0, tex_xy_xyyz_0, tex_xy_xyzz_0, tex_xy_xzzz_0, tex_xy_yyyy_0, \
                                         tex_xy_yyyz_0, tex_xy_yyzz_0, tex_xy_yzzz_0, tex_xy_zzzz_0, tex_y_xxx_0, \
                                         tex_y_xxx_1, tex_y_xxxx_0, tex_y_xxxx_1, tex_y_xxxy_0, tex_y_xxxy_1, tex_y_xxxz_0, \
                                         tex_y_xxxz_1, tex_y_xxy_0, tex_y_xxy_1, tex_y_xxyy_0, tex_y_xxyy_1, tex_y_xxyz_0, \
                                         tex_y_xxyz_1, tex_y_xxz_0, tex_y_xxz_1, tex_y_xxzz_0, tex_y_xxzz_1, tex_y_xyy_0, \
                                         tex_y_xyy_1, tex_y_xyyy_0, tex_y_xyyy_1, tex_y_xyyz_0, tex_y_xyyz_1, tex_y_xyz_0, \
                                         tex_y_xyz_1, tex_y_xyzz_0, tex_y_xyzz_1, tex_y_xzz_0, tex_y_xzz_1, tex_y_xzzz_0, \
                                         tex_y_xzzz_1, tex_y_yyy_0, tex_y_yyy_1, tex_y_yyyy_0, tex_y_yyyy_1, tex_y_yyyz_0, \
                                         tex_y_yyyz_1, tex_y_yyz_0, tex_y_yyz_1, tex_y_yyzz_0, tex_y_yyzz_1, tex_y_yzz_0, \
                                         tex_y_yzz_1, tex_y_yzzz_0, tex_y_yzzz_1, tex_y_zzz_0, tex_y_zzz_1, tex_y_zzzz_0, \
                                         tex_y_zzzz_1, tey_xy_xxxx_0, tey_xy_xxxy_0, tey_xy_xxxz_0, tey_xy_xxyy_0, \
                                         tey_xy_xxyz_0, tey_xy_xxzz_0, tey_xy_xyyy_0, tey_xy_xyyz_0, tey_xy_xyzz_0, \
                                         tey_xy_xzzz_0, tey_xy_yyyy_0, tey_xy_yyyz_0, tey_xy_yyzz_0, tey_xy_yzzz_0, \
                                         tey_xy_zzzz_0, tey_y_xxx_0, tey_y_xxx_1, tey_y_xxxx_0, tey_y_xxxx_1, tey_y_xxxy_0, \
                                         tey_y_xxxy_1, tey_y_xxxz_0, tey_y_xxxz_1, tey_y_xxy_0, tey_y_xxy_1, tey_y_xxyy_0, \
                                         tey_y_xxyy_1, tey_y_xxyz_0, tey_y_xxyz_1, tey_y_xxz_0, tey_y_xxz_1, tey_y_xxzz_0, \
                                         tey_y_xxzz_1, tey_y_xyy_0, tey_y_xyy_1, tey_y_xyyy_0, tey_y_xyyy_1, tey_y_xyyz_0, \
                                         tey_y_xyyz_1, tey_y_xyz_0, tey_y_xyz_1, tey_y_xyzz_0, tey_y_xyzz_1, tey_y_xzz_0, \
                                         tey_y_xzz_1, tey_y_xzzz_0, tey_y_xzzz_1, tey_y_yyy_0, tey_y_yyy_1, tey_y_yyyy_0, \
                                         tey_y_yyyy_1, tey_y_yyyz_0, tey_y_yyyz_1, tey_y_yyz_0, tey_y_yyz_1, tey_y_yyzz_0, \
                                         tey_y_yyzz_1, tey_y_yzz_0, tey_y_yzz_1, tey_y_yzzz_0, tey_y_yzzz_1, tey_y_zzz_0, \
                                         tey_y_zzz_1, tey_y_zzzz_0, tey_y_zzzz_1, tez_xy_xxxx_0, tez_xy_xxxy_0, \
                                         tez_xy_xxxz_0, tez_xy_xxyy_0, tez_xy_xxyz_0, tez_xy_xxzz_0, tez_xy_xyyy_0, \
                                         tez_xy_xyyz_0, tez_xy_xyzz_0, tez_xy_xzzz_0, tez_xy_yyyy_0, tez_xy_yyyz_0, \
                                         tez_xy_yyzz_0, tez_xy_yzzz_0, tez_xy_zzzz_0, tez_y_xxx_0, tez_y_xxx_1, tez_y_xxxx_0, \
                                         tez_y_xxxx_1, tez_y_xxxy_0, tez_y_xxxy_1, tez_y_xxxz_0, tez_y_xxxz_1, tez_y_xxy_0, \
                                         tez_y_xxy_1, tez_y_xxyy_0, tez_y_xxyy_1, tez_y_xxyz_0, tez_y_xxyz_1, tez_y_xxz_0, \
                                         tez_y_xxz_1, tez_y_xxzz_0, tez_y_xxzz_1, tez_y_xyy_0, tez_y_xyy_1, tez_y_xyyy_0, \
                                         tez_y_xyyy_1, tez_y_xyyz_0, tez_y_xyyz_1, tez_y_xyz_0, tez_y_xyz_1, tez_y_xyzz_0, \
                                         tez_y_xyzz_1, tez_y_xzz_0, tez_y_xzz_1, tez_y_xzzz_0, tez_y_xzzz_1, tez_y_yyy_0, \
                                         tez_y_yyy_1, tez_y_yyyy_0, tez_y_yyyy_1, tez_y_yyyz_0, tez_y_yyyz_1, tez_y_yyz_0, \
                                         tez_y_yyz_1, tez_y_yyzz_0, tez_y_yyzz_1, tez_y_yzz_0, tez_y_yzz_1, tez_y_yzzz_0, \
                                         tez_y_yzzz_1, tez_y_zzz_0, tez_y_zzz_1, tez_y_zzzz_0, tez_y_zzzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xy_xxxx_0[j] = pa_x[j] * tex_y_xxxx_0[j] - pc_x[j] * tex_y_xxxx_1[j] + 2.0 * fl1_fx * tex_y_xxx_0[j] - 2.0 * fl1_fx * tex_y_xxx_1[j] + ta_y_xxxx_1[j];

                    tey_xy_xxxx_0[j] = pa_x[j] * tey_y_xxxx_0[j] - pc_x[j] * tey_y_xxxx_1[j] + 2.0 * fl1_fx * tey_y_xxx_0[j] - 2.0 * fl1_fx * tey_y_xxx_1[j];

                    tez_xy_xxxx_0[j] = pa_x[j] * tez_y_xxxx_0[j] - pc_x[j] * tez_y_xxxx_1[j] + 2.0 * fl1_fx * tez_y_xxx_0[j] - 2.0 * fl1_fx * tez_y_xxx_1[j];

                    tex_xy_xxxy_0[j] = pa_x[j] * tex_y_xxxy_0[j] - pc_x[j] * tex_y_xxxy_1[j] + 1.5 * fl1_fx * tex_y_xxy_0[j] - 1.5 * fl1_fx * tex_y_xxy_1[j] + ta_y_xxxy_1[j];

                    tey_xy_xxxy_0[j] = pa_x[j] * tey_y_xxxy_0[j] - pc_x[j] * tey_y_xxxy_1[j] + 1.5 * fl1_fx * tey_y_xxy_0[j] - 1.5 * fl1_fx * tey_y_xxy_1[j];

                    tez_xy_xxxy_0[j] = pa_x[j] * tez_y_xxxy_0[j] - pc_x[j] * tez_y_xxxy_1[j] + 1.5 * fl1_fx * tez_y_xxy_0[j] - 1.5 * fl1_fx * tez_y_xxy_1[j];

                    tex_xy_xxxz_0[j] = pa_x[j] * tex_y_xxxz_0[j] - pc_x[j] * tex_y_xxxz_1[j] + 1.5 * fl1_fx * tex_y_xxz_0[j] - 1.5 * fl1_fx * tex_y_xxz_1[j] + ta_y_xxxz_1[j];

                    tey_xy_xxxz_0[j] = pa_x[j] * tey_y_xxxz_0[j] - pc_x[j] * tey_y_xxxz_1[j] + 1.5 * fl1_fx * tey_y_xxz_0[j] - 1.5 * fl1_fx * tey_y_xxz_1[j];

                    tez_xy_xxxz_0[j] = pa_x[j] * tez_y_xxxz_0[j] - pc_x[j] * tez_y_xxxz_1[j] + 1.5 * fl1_fx * tez_y_xxz_0[j] - 1.5 * fl1_fx * tez_y_xxz_1[j];

                    tex_xy_xxyy_0[j] = pa_x[j] * tex_y_xxyy_0[j] - pc_x[j] * tex_y_xxyy_1[j] + fl1_fx * tex_y_xyy_0[j] - fl1_fx * tex_y_xyy_1[j] + ta_y_xxyy_1[j];

                    tey_xy_xxyy_0[j] = pa_x[j] * tey_y_xxyy_0[j] - pc_x[j] * tey_y_xxyy_1[j] + fl1_fx * tey_y_xyy_0[j] - fl1_fx * tey_y_xyy_1[j];

                    tez_xy_xxyy_0[j] = pa_x[j] * tez_y_xxyy_0[j] - pc_x[j] * tez_y_xxyy_1[j] + fl1_fx * tez_y_xyy_0[j] - fl1_fx * tez_y_xyy_1[j];

                    tex_xy_xxyz_0[j] = pa_x[j] * tex_y_xxyz_0[j] - pc_x[j] * tex_y_xxyz_1[j] + fl1_fx * tex_y_xyz_0[j] - fl1_fx * tex_y_xyz_1[j] + ta_y_xxyz_1[j];

                    tey_xy_xxyz_0[j] = pa_x[j] * tey_y_xxyz_0[j] - pc_x[j] * tey_y_xxyz_1[j] + fl1_fx * tey_y_xyz_0[j] - fl1_fx * tey_y_xyz_1[j];

                    tez_xy_xxyz_0[j] = pa_x[j] * tez_y_xxyz_0[j] - pc_x[j] * tez_y_xxyz_1[j] + fl1_fx * tez_y_xyz_0[j] - fl1_fx * tez_y_xyz_1[j];

                    tex_xy_xxzz_0[j] = pa_x[j] * tex_y_xxzz_0[j] - pc_x[j] * tex_y_xxzz_1[j] + fl1_fx * tex_y_xzz_0[j] - fl1_fx * tex_y_xzz_1[j] + ta_y_xxzz_1[j];

                    tey_xy_xxzz_0[j] = pa_x[j] * tey_y_xxzz_0[j] - pc_x[j] * tey_y_xxzz_1[j] + fl1_fx * tey_y_xzz_0[j] - fl1_fx * tey_y_xzz_1[j];

                    tez_xy_xxzz_0[j] = pa_x[j] * tez_y_xxzz_0[j] - pc_x[j] * tez_y_xxzz_1[j] + fl1_fx * tez_y_xzz_0[j] - fl1_fx * tez_y_xzz_1[j];

                    tex_xy_xyyy_0[j] = pa_x[j] * tex_y_xyyy_0[j] - pc_x[j] * tex_y_xyyy_1[j] + 0.5 * fl1_fx * tex_y_yyy_0[j] - 0.5 * fl1_fx * tex_y_yyy_1[j] + ta_y_xyyy_1[j];

                    tey_xy_xyyy_0[j] = pa_x[j] * tey_y_xyyy_0[j] - pc_x[j] * tey_y_xyyy_1[j] + 0.5 * fl1_fx * tey_y_yyy_0[j] - 0.5 * fl1_fx * tey_y_yyy_1[j];

                    tez_xy_xyyy_0[j] = pa_x[j] * tez_y_xyyy_0[j] - pc_x[j] * tez_y_xyyy_1[j] + 0.5 * fl1_fx * tez_y_yyy_0[j] - 0.5 * fl1_fx * tez_y_yyy_1[j];

                    tex_xy_xyyz_0[j] = pa_x[j] * tex_y_xyyz_0[j] - pc_x[j] * tex_y_xyyz_1[j] + 0.5 * fl1_fx * tex_y_yyz_0[j] - 0.5 * fl1_fx * tex_y_yyz_1[j] + ta_y_xyyz_1[j];

                    tey_xy_xyyz_0[j] = pa_x[j] * tey_y_xyyz_0[j] - pc_x[j] * tey_y_xyyz_1[j] + 0.5 * fl1_fx * tey_y_yyz_0[j] - 0.5 * fl1_fx * tey_y_yyz_1[j];

                    tez_xy_xyyz_0[j] = pa_x[j] * tez_y_xyyz_0[j] - pc_x[j] * tez_y_xyyz_1[j] + 0.5 * fl1_fx * tez_y_yyz_0[j] - 0.5 * fl1_fx * tez_y_yyz_1[j];

                    tex_xy_xyzz_0[j] = pa_x[j] * tex_y_xyzz_0[j] - pc_x[j] * tex_y_xyzz_1[j] + 0.5 * fl1_fx * tex_y_yzz_0[j] - 0.5 * fl1_fx * tex_y_yzz_1[j] + ta_y_xyzz_1[j];

                    tey_xy_xyzz_0[j] = pa_x[j] * tey_y_xyzz_0[j] - pc_x[j] * tey_y_xyzz_1[j] + 0.5 * fl1_fx * tey_y_yzz_0[j] - 0.5 * fl1_fx * tey_y_yzz_1[j];

                    tez_xy_xyzz_0[j] = pa_x[j] * tez_y_xyzz_0[j] - pc_x[j] * tez_y_xyzz_1[j] + 0.5 * fl1_fx * tez_y_yzz_0[j] - 0.5 * fl1_fx * tez_y_yzz_1[j];

                    tex_xy_xzzz_0[j] = pa_x[j] * tex_y_xzzz_0[j] - pc_x[j] * tex_y_xzzz_1[j] + 0.5 * fl1_fx * tex_y_zzz_0[j] - 0.5 * fl1_fx * tex_y_zzz_1[j] + ta_y_xzzz_1[j];

                    tey_xy_xzzz_0[j] = pa_x[j] * tey_y_xzzz_0[j] - pc_x[j] * tey_y_xzzz_1[j] + 0.5 * fl1_fx * tey_y_zzz_0[j] - 0.5 * fl1_fx * tey_y_zzz_1[j];

                    tez_xy_xzzz_0[j] = pa_x[j] * tez_y_xzzz_0[j] - pc_x[j] * tez_y_xzzz_1[j] + 0.5 * fl1_fx * tez_y_zzz_0[j] - 0.5 * fl1_fx * tez_y_zzz_1[j];

                    tex_xy_yyyy_0[j] = pa_x[j] * tex_y_yyyy_0[j] - pc_x[j] * tex_y_yyyy_1[j] + ta_y_yyyy_1[j];

                    tey_xy_yyyy_0[j] = pa_x[j] * tey_y_yyyy_0[j] - pc_x[j] * tey_y_yyyy_1[j];

                    tez_xy_yyyy_0[j] = pa_x[j] * tez_y_yyyy_0[j] - pc_x[j] * tez_y_yyyy_1[j];

                    tex_xy_yyyz_0[j] = pa_x[j] * tex_y_yyyz_0[j] - pc_x[j] * tex_y_yyyz_1[j] + ta_y_yyyz_1[j];

                    tey_xy_yyyz_0[j] = pa_x[j] * tey_y_yyyz_0[j] - pc_x[j] * tey_y_yyyz_1[j];

                    tez_xy_yyyz_0[j] = pa_x[j] * tez_y_yyyz_0[j] - pc_x[j] * tez_y_yyyz_1[j];

                    tex_xy_yyzz_0[j] = pa_x[j] * tex_y_yyzz_0[j] - pc_x[j] * tex_y_yyzz_1[j] + ta_y_yyzz_1[j];

                    tey_xy_yyzz_0[j] = pa_x[j] * tey_y_yyzz_0[j] - pc_x[j] * tey_y_yyzz_1[j];

                    tez_xy_yyzz_0[j] = pa_x[j] * tez_y_yyzz_0[j] - pc_x[j] * tez_y_yyzz_1[j];

                    tex_xy_yzzz_0[j] = pa_x[j] * tex_y_yzzz_0[j] - pc_x[j] * tex_y_yzzz_1[j] + ta_y_yzzz_1[j];

                    tey_xy_yzzz_0[j] = pa_x[j] * tey_y_yzzz_0[j] - pc_x[j] * tey_y_yzzz_1[j];

                    tez_xy_yzzz_0[j] = pa_x[j] * tez_y_yzzz_0[j] - pc_x[j] * tez_y_yzzz_1[j];

                    tex_xy_zzzz_0[j] = pa_x[j] * tex_y_zzzz_0[j] - pc_x[j] * tex_y_zzzz_1[j] + ta_y_zzzz_1[j];

                    tey_xy_zzzz_0[j] = pa_x[j] * tey_y_zzzz_0[j] - pc_x[j] * tey_y_zzzz_1[j];

                    tez_xy_zzzz_0[j] = pa_x[j] * tez_y_zzzz_0[j] - pc_x[j] * tez_y_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDG_90_135(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (90,135)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
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

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_x = pcDistances.data(3 * idx);

                // set up pointers to auxilary integrals

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

                auto tex_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 20); 

                auto tey_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 21); 

                auto tey_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 22); 

                auto tey_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 23); 

                auto tey_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 24); 

                auto tey_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 25); 

                auto tey_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 26); 

                auto tey_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 27); 

                auto tey_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 28); 

                auto tey_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 29); 

                auto tey_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 29); 

                auto tex_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 20); 

                auto tey_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 21); 

                auto tey_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 22); 

                auto tey_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 23); 

                auto tey_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 24); 

                auto tey_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 25); 

                auto tey_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 26); 

                auto tey_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 27); 

                auto tey_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 28); 

                auto tey_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 29); 

                auto tey_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 29); 

                auto ta_z_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 30); 

                auto ta_z_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 31); 

                auto ta_z_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 32); 

                auto ta_z_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 33); 

                auto ta_z_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 34); 

                auto ta_z_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 35); 

                auto ta_z_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 36); 

                auto ta_z_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 37); 

                auto ta_z_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 38); 

                auto ta_z_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 39); 

                auto ta_z_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 40); 

                auto ta_z_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 41); 

                auto ta_z_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 42); 

                auto ta_z_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 43); 

                auto ta_z_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 44); 

                // set up pointers to integrals

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

                // Batch of Integrals (90,135)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_z_xxxx_1, ta_z_xxxy_1, ta_z_xxxz_1, ta_z_xxyy_1, \
                                         ta_z_xxyz_1, ta_z_xxzz_1, ta_z_xyyy_1, ta_z_xyyz_1, ta_z_xyzz_1, ta_z_xzzz_1, \
                                         ta_z_yyyy_1, ta_z_yyyz_1, ta_z_yyzz_1, ta_z_yzzz_1, ta_z_zzzz_1, tex_xz_xxxx_0, \
                                         tex_xz_xxxy_0, tex_xz_xxxz_0, tex_xz_xxyy_0, tex_xz_xxyz_0, tex_xz_xxzz_0, \
                                         tex_xz_xyyy_0, tex_xz_xyyz_0, tex_xz_xyzz_0, tex_xz_xzzz_0, tex_xz_yyyy_0, \
                                         tex_xz_yyyz_0, tex_xz_yyzz_0, tex_xz_yzzz_0, tex_xz_zzzz_0, tex_z_xxx_0, \
                                         tex_z_xxx_1, tex_z_xxxx_0, tex_z_xxxx_1, tex_z_xxxy_0, tex_z_xxxy_1, tex_z_xxxz_0, \
                                         tex_z_xxxz_1, tex_z_xxy_0, tex_z_xxy_1, tex_z_xxyy_0, tex_z_xxyy_1, tex_z_xxyz_0, \
                                         tex_z_xxyz_1, tex_z_xxz_0, tex_z_xxz_1, tex_z_xxzz_0, tex_z_xxzz_1, tex_z_xyy_0, \
                                         tex_z_xyy_1, tex_z_xyyy_0, tex_z_xyyy_1, tex_z_xyyz_0, tex_z_xyyz_1, tex_z_xyz_0, \
                                         tex_z_xyz_1, tex_z_xyzz_0, tex_z_xyzz_1, tex_z_xzz_0, tex_z_xzz_1, tex_z_xzzz_0, \
                                         tex_z_xzzz_1, tex_z_yyy_0, tex_z_yyy_1, tex_z_yyyy_0, tex_z_yyyy_1, tex_z_yyyz_0, \
                                         tex_z_yyyz_1, tex_z_yyz_0, tex_z_yyz_1, tex_z_yyzz_0, tex_z_yyzz_1, tex_z_yzz_0, \
                                         tex_z_yzz_1, tex_z_yzzz_0, tex_z_yzzz_1, tex_z_zzz_0, tex_z_zzz_1, tex_z_zzzz_0, \
                                         tex_z_zzzz_1, tey_xz_xxxx_0, tey_xz_xxxy_0, tey_xz_xxxz_0, tey_xz_xxyy_0, \
                                         tey_xz_xxyz_0, tey_xz_xxzz_0, tey_xz_xyyy_0, tey_xz_xyyz_0, tey_xz_xyzz_0, \
                                         tey_xz_xzzz_0, tey_xz_yyyy_0, tey_xz_yyyz_0, tey_xz_yyzz_0, tey_xz_yzzz_0, \
                                         tey_xz_zzzz_0, tey_z_xxx_0, tey_z_xxx_1, tey_z_xxxx_0, tey_z_xxxx_1, tey_z_xxxy_0, \
                                         tey_z_xxxy_1, tey_z_xxxz_0, tey_z_xxxz_1, tey_z_xxy_0, tey_z_xxy_1, tey_z_xxyy_0, \
                                         tey_z_xxyy_1, tey_z_xxyz_0, tey_z_xxyz_1, tey_z_xxz_0, tey_z_xxz_1, tey_z_xxzz_0, \
                                         tey_z_xxzz_1, tey_z_xyy_0, tey_z_xyy_1, tey_z_xyyy_0, tey_z_xyyy_1, tey_z_xyyz_0, \
                                         tey_z_xyyz_1, tey_z_xyz_0, tey_z_xyz_1, tey_z_xyzz_0, tey_z_xyzz_1, tey_z_xzz_0, \
                                         tey_z_xzz_1, tey_z_xzzz_0, tey_z_xzzz_1, tey_z_yyy_0, tey_z_yyy_1, tey_z_yyyy_0, \
                                         tey_z_yyyy_1, tey_z_yyyz_0, tey_z_yyyz_1, tey_z_yyz_0, tey_z_yyz_1, tey_z_yyzz_0, \
                                         tey_z_yyzz_1, tey_z_yzz_0, tey_z_yzz_1, tey_z_yzzz_0, tey_z_yzzz_1, tey_z_zzz_0, \
                                         tey_z_zzz_1, tey_z_zzzz_0, tey_z_zzzz_1, tez_xz_xxxx_0, tez_xz_xxxy_0, \
                                         tez_xz_xxxz_0, tez_xz_xxyy_0, tez_xz_xxyz_0, tez_xz_xxzz_0, tez_xz_xyyy_0, \
                                         tez_xz_xyyz_0, tez_xz_xyzz_0, tez_xz_xzzz_0, tez_xz_yyyy_0, tez_xz_yyyz_0, \
                                         tez_xz_yyzz_0, tez_xz_yzzz_0, tez_xz_zzzz_0, tez_z_xxx_0, tez_z_xxx_1, tez_z_xxxx_0, \
                                         tez_z_xxxx_1, tez_z_xxxy_0, tez_z_xxxy_1, tez_z_xxxz_0, tez_z_xxxz_1, tez_z_xxy_0, \
                                         tez_z_xxy_1, tez_z_xxyy_0, tez_z_xxyy_1, tez_z_xxyz_0, tez_z_xxyz_1, tez_z_xxz_0, \
                                         tez_z_xxz_1, tez_z_xxzz_0, tez_z_xxzz_1, tez_z_xyy_0, tez_z_xyy_1, tez_z_xyyy_0, \
                                         tez_z_xyyy_1, tez_z_xyyz_0, tez_z_xyyz_1, tez_z_xyz_0, tez_z_xyz_1, tez_z_xyzz_0, \
                                         tez_z_xyzz_1, tez_z_xzz_0, tez_z_xzz_1, tez_z_xzzz_0, tez_z_xzzz_1, tez_z_yyy_0, \
                                         tez_z_yyy_1, tez_z_yyyy_0, tez_z_yyyy_1, tez_z_yyyz_0, tez_z_yyyz_1, tez_z_yyz_0, \
                                         tez_z_yyz_1, tez_z_yyzz_0, tez_z_yyzz_1, tez_z_yzz_0, tez_z_yzz_1, tez_z_yzzz_0, \
                                         tez_z_yzzz_1, tez_z_zzz_0, tez_z_zzz_1, tez_z_zzzz_0, tez_z_zzzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xz_xxxx_0[j] = pa_x[j] * tex_z_xxxx_0[j] - pc_x[j] * tex_z_xxxx_1[j] + 2.0 * fl1_fx * tex_z_xxx_0[j] - 2.0 * fl1_fx * tex_z_xxx_1[j] + ta_z_xxxx_1[j];

                    tey_xz_xxxx_0[j] = pa_x[j] * tey_z_xxxx_0[j] - pc_x[j] * tey_z_xxxx_1[j] + 2.0 * fl1_fx * tey_z_xxx_0[j] - 2.0 * fl1_fx * tey_z_xxx_1[j];

                    tez_xz_xxxx_0[j] = pa_x[j] * tez_z_xxxx_0[j] - pc_x[j] * tez_z_xxxx_1[j] + 2.0 * fl1_fx * tez_z_xxx_0[j] - 2.0 * fl1_fx * tez_z_xxx_1[j];

                    tex_xz_xxxy_0[j] = pa_x[j] * tex_z_xxxy_0[j] - pc_x[j] * tex_z_xxxy_1[j] + 1.5 * fl1_fx * tex_z_xxy_0[j] - 1.5 * fl1_fx * tex_z_xxy_1[j] + ta_z_xxxy_1[j];

                    tey_xz_xxxy_0[j] = pa_x[j] * tey_z_xxxy_0[j] - pc_x[j] * tey_z_xxxy_1[j] + 1.5 * fl1_fx * tey_z_xxy_0[j] - 1.5 * fl1_fx * tey_z_xxy_1[j];

                    tez_xz_xxxy_0[j] = pa_x[j] * tez_z_xxxy_0[j] - pc_x[j] * tez_z_xxxy_1[j] + 1.5 * fl1_fx * tez_z_xxy_0[j] - 1.5 * fl1_fx * tez_z_xxy_1[j];

                    tex_xz_xxxz_0[j] = pa_x[j] * tex_z_xxxz_0[j] - pc_x[j] * tex_z_xxxz_1[j] + 1.5 * fl1_fx * tex_z_xxz_0[j] - 1.5 * fl1_fx * tex_z_xxz_1[j] + ta_z_xxxz_1[j];

                    tey_xz_xxxz_0[j] = pa_x[j] * tey_z_xxxz_0[j] - pc_x[j] * tey_z_xxxz_1[j] + 1.5 * fl1_fx * tey_z_xxz_0[j] - 1.5 * fl1_fx * tey_z_xxz_1[j];

                    tez_xz_xxxz_0[j] = pa_x[j] * tez_z_xxxz_0[j] - pc_x[j] * tez_z_xxxz_1[j] + 1.5 * fl1_fx * tez_z_xxz_0[j] - 1.5 * fl1_fx * tez_z_xxz_1[j];

                    tex_xz_xxyy_0[j] = pa_x[j] * tex_z_xxyy_0[j] - pc_x[j] * tex_z_xxyy_1[j] + fl1_fx * tex_z_xyy_0[j] - fl1_fx * tex_z_xyy_1[j] + ta_z_xxyy_1[j];

                    tey_xz_xxyy_0[j] = pa_x[j] * tey_z_xxyy_0[j] - pc_x[j] * tey_z_xxyy_1[j] + fl1_fx * tey_z_xyy_0[j] - fl1_fx * tey_z_xyy_1[j];

                    tez_xz_xxyy_0[j] = pa_x[j] * tez_z_xxyy_0[j] - pc_x[j] * tez_z_xxyy_1[j] + fl1_fx * tez_z_xyy_0[j] - fl1_fx * tez_z_xyy_1[j];

                    tex_xz_xxyz_0[j] = pa_x[j] * tex_z_xxyz_0[j] - pc_x[j] * tex_z_xxyz_1[j] + fl1_fx * tex_z_xyz_0[j] - fl1_fx * tex_z_xyz_1[j] + ta_z_xxyz_1[j];

                    tey_xz_xxyz_0[j] = pa_x[j] * tey_z_xxyz_0[j] - pc_x[j] * tey_z_xxyz_1[j] + fl1_fx * tey_z_xyz_0[j] - fl1_fx * tey_z_xyz_1[j];

                    tez_xz_xxyz_0[j] = pa_x[j] * tez_z_xxyz_0[j] - pc_x[j] * tez_z_xxyz_1[j] + fl1_fx * tez_z_xyz_0[j] - fl1_fx * tez_z_xyz_1[j];

                    tex_xz_xxzz_0[j] = pa_x[j] * tex_z_xxzz_0[j] - pc_x[j] * tex_z_xxzz_1[j] + fl1_fx * tex_z_xzz_0[j] - fl1_fx * tex_z_xzz_1[j] + ta_z_xxzz_1[j];

                    tey_xz_xxzz_0[j] = pa_x[j] * tey_z_xxzz_0[j] - pc_x[j] * tey_z_xxzz_1[j] + fl1_fx * tey_z_xzz_0[j] - fl1_fx * tey_z_xzz_1[j];

                    tez_xz_xxzz_0[j] = pa_x[j] * tez_z_xxzz_0[j] - pc_x[j] * tez_z_xxzz_1[j] + fl1_fx * tez_z_xzz_0[j] - fl1_fx * tez_z_xzz_1[j];

                    tex_xz_xyyy_0[j] = pa_x[j] * tex_z_xyyy_0[j] - pc_x[j] * tex_z_xyyy_1[j] + 0.5 * fl1_fx * tex_z_yyy_0[j] - 0.5 * fl1_fx * tex_z_yyy_1[j] + ta_z_xyyy_1[j];

                    tey_xz_xyyy_0[j] = pa_x[j] * tey_z_xyyy_0[j] - pc_x[j] * tey_z_xyyy_1[j] + 0.5 * fl1_fx * tey_z_yyy_0[j] - 0.5 * fl1_fx * tey_z_yyy_1[j];

                    tez_xz_xyyy_0[j] = pa_x[j] * tez_z_xyyy_0[j] - pc_x[j] * tez_z_xyyy_1[j] + 0.5 * fl1_fx * tez_z_yyy_0[j] - 0.5 * fl1_fx * tez_z_yyy_1[j];

                    tex_xz_xyyz_0[j] = pa_x[j] * tex_z_xyyz_0[j] - pc_x[j] * tex_z_xyyz_1[j] + 0.5 * fl1_fx * tex_z_yyz_0[j] - 0.5 * fl1_fx * tex_z_yyz_1[j] + ta_z_xyyz_1[j];

                    tey_xz_xyyz_0[j] = pa_x[j] * tey_z_xyyz_0[j] - pc_x[j] * tey_z_xyyz_1[j] + 0.5 * fl1_fx * tey_z_yyz_0[j] - 0.5 * fl1_fx * tey_z_yyz_1[j];

                    tez_xz_xyyz_0[j] = pa_x[j] * tez_z_xyyz_0[j] - pc_x[j] * tez_z_xyyz_1[j] + 0.5 * fl1_fx * tez_z_yyz_0[j] - 0.5 * fl1_fx * tez_z_yyz_1[j];

                    tex_xz_xyzz_0[j] = pa_x[j] * tex_z_xyzz_0[j] - pc_x[j] * tex_z_xyzz_1[j] + 0.5 * fl1_fx * tex_z_yzz_0[j] - 0.5 * fl1_fx * tex_z_yzz_1[j] + ta_z_xyzz_1[j];

                    tey_xz_xyzz_0[j] = pa_x[j] * tey_z_xyzz_0[j] - pc_x[j] * tey_z_xyzz_1[j] + 0.5 * fl1_fx * tey_z_yzz_0[j] - 0.5 * fl1_fx * tey_z_yzz_1[j];

                    tez_xz_xyzz_0[j] = pa_x[j] * tez_z_xyzz_0[j] - pc_x[j] * tez_z_xyzz_1[j] + 0.5 * fl1_fx * tez_z_yzz_0[j] - 0.5 * fl1_fx * tez_z_yzz_1[j];

                    tex_xz_xzzz_0[j] = pa_x[j] * tex_z_xzzz_0[j] - pc_x[j] * tex_z_xzzz_1[j] + 0.5 * fl1_fx * tex_z_zzz_0[j] - 0.5 * fl1_fx * tex_z_zzz_1[j] + ta_z_xzzz_1[j];

                    tey_xz_xzzz_0[j] = pa_x[j] * tey_z_xzzz_0[j] - pc_x[j] * tey_z_xzzz_1[j] + 0.5 * fl1_fx * tey_z_zzz_0[j] - 0.5 * fl1_fx * tey_z_zzz_1[j];

                    tez_xz_xzzz_0[j] = pa_x[j] * tez_z_xzzz_0[j] - pc_x[j] * tez_z_xzzz_1[j] + 0.5 * fl1_fx * tez_z_zzz_0[j] - 0.5 * fl1_fx * tez_z_zzz_1[j];

                    tex_xz_yyyy_0[j] = pa_x[j] * tex_z_yyyy_0[j] - pc_x[j] * tex_z_yyyy_1[j] + ta_z_yyyy_1[j];

                    tey_xz_yyyy_0[j] = pa_x[j] * tey_z_yyyy_0[j] - pc_x[j] * tey_z_yyyy_1[j];

                    tez_xz_yyyy_0[j] = pa_x[j] * tez_z_yyyy_0[j] - pc_x[j] * tez_z_yyyy_1[j];

                    tex_xz_yyyz_0[j] = pa_x[j] * tex_z_yyyz_0[j] - pc_x[j] * tex_z_yyyz_1[j] + ta_z_yyyz_1[j];

                    tey_xz_yyyz_0[j] = pa_x[j] * tey_z_yyyz_0[j] - pc_x[j] * tey_z_yyyz_1[j];

                    tez_xz_yyyz_0[j] = pa_x[j] * tez_z_yyyz_0[j] - pc_x[j] * tez_z_yyyz_1[j];

                    tex_xz_yyzz_0[j] = pa_x[j] * tex_z_yyzz_0[j] - pc_x[j] * tex_z_yyzz_1[j] + ta_z_yyzz_1[j];

                    tey_xz_yyzz_0[j] = pa_x[j] * tey_z_yyzz_0[j] - pc_x[j] * tey_z_yyzz_1[j];

                    tez_xz_yyzz_0[j] = pa_x[j] * tez_z_yyzz_0[j] - pc_x[j] * tez_z_yyzz_1[j];

                    tex_xz_yzzz_0[j] = pa_x[j] * tex_z_yzzz_0[j] - pc_x[j] * tex_z_yzzz_1[j] + ta_z_yzzz_1[j];

                    tey_xz_yzzz_0[j] = pa_x[j] * tey_z_yzzz_0[j] - pc_x[j] * tey_z_yzzz_1[j];

                    tez_xz_yzzz_0[j] = pa_x[j] * tez_z_yzzz_0[j] - pc_x[j] * tez_z_yzzz_1[j];

                    tex_xz_zzzz_0[j] = pa_x[j] * tex_z_zzzz_0[j] - pc_x[j] * tex_z_zzzz_1[j] + ta_z_zzzz_1[j];

                    tey_xz_zzzz_0[j] = pa_x[j] * tey_z_zzzz_0[j] - pc_x[j] * tey_z_zzzz_1[j];

                    tez_xz_zzzz_0[j] = pa_x[j] * tez_z_zzzz_0[j] - pc_x[j] * tez_z_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDG_135_180(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (135,180)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tex_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx); 

                auto tey_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx); 

                auto tez_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx); 

                auto tex_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 1); 

                auto tey_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 1); 

                auto tez_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 1); 

                auto tex_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 2); 

                auto tey_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 2); 

                auto tez_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 2); 

                auto tex_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 3); 

                auto tey_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 3); 

                auto tez_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 3); 

                auto tex_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 4); 

                auto tey_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 4); 

                auto tez_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 4); 

                auto tex_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 5); 

                auto tey_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 5); 

                auto tez_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 5); 

                auto tex_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 6); 

                auto tey_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 6); 

                auto tez_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 6); 

                auto tex_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 7); 

                auto tey_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 7); 

                auto tez_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 7); 

                auto tex_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 8); 

                auto tey_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 8); 

                auto tez_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 8); 

                auto tex_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 9); 

                auto tey_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 9); 

                auto tez_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 9); 

                auto tex_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 10); 

                auto tey_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 10); 

                auto tez_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 10); 

                auto tex_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 11); 

                auto tey_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 11); 

                auto tez_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 11); 

                auto tex_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 12); 

                auto tey_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 12); 

                auto tez_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 12); 

                auto tex_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 13); 

                auto tey_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 13); 

                auto tez_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 13); 

                auto tex_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 14); 

                auto tey_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 14); 

                auto tez_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 14); 

                auto tex_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx); 

                auto tey_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx); 

                auto tez_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx); 

                auto tex_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 1); 

                auto tey_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 1); 

                auto tez_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 1); 

                auto tex_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 2); 

                auto tey_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 2); 

                auto tez_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 2); 

                auto tex_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 3); 

                auto tey_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 3); 

                auto tez_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 3); 

                auto tex_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 4); 

                auto tey_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 4); 

                auto tez_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 4); 

                auto tex_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 5); 

                auto tey_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 5); 

                auto tez_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 5); 

                auto tex_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 6); 

                auto tey_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 6); 

                auto tez_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 6); 

                auto tex_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 7); 

                auto tey_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 7); 

                auto tez_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 7); 

                auto tex_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 8); 

                auto tey_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 8); 

                auto tez_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 8); 

                auto tex_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 9); 

                auto tey_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 9); 

                auto tez_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 9); 

                auto tex_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 10); 

                auto tey_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 10); 

                auto tez_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 10); 

                auto tex_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 11); 

                auto tey_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 11); 

                auto tez_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 11); 

                auto tex_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 12); 

                auto tey_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 12); 

                auto tez_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 12); 

                auto tex_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 13); 

                auto tey_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 13); 

                auto tez_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 13); 

                auto tex_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 14); 

                auto tey_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 14); 

                auto tez_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 14); 

                auto tex_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 10); 

                auto tey_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 11); 

                auto tey_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 12); 

                auto tey_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 13); 

                auto tey_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 14); 

                auto tey_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 14); 

                auto tex_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 15); 

                auto tey_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 16); 

                auto tey_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 17); 

                auto tey_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 18); 

                auto tey_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 19); 

                auto tey_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 19); 

                auto tex_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 10); 

                auto tey_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 10); 

                auto tez_y_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 10); 

                auto tex_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 11); 

                auto tey_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 11); 

                auto tez_y_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 11); 

                auto tex_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 12); 

                auto tey_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 12); 

                auto tez_y_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 12); 

                auto tex_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 13); 

                auto tey_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 13); 

                auto tez_y_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 13); 

                auto tex_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 14); 

                auto tey_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 14); 

                auto tez_y_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 14); 

                auto tex_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 15); 

                auto tey_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 15); 

                auto tez_y_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 15); 

                auto tex_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 16); 

                auto tey_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 16); 

                auto tez_y_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 16); 

                auto tex_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 17); 

                auto tey_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 17); 

                auto tez_y_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 17); 

                auto tex_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 18); 

                auto tey_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 18); 

                auto tez_y_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 18); 

                auto tex_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 19); 

                auto tey_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 19); 

                auto tez_y_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 19); 

                auto ta_y_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 15); 

                auto ta_y_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 16); 

                auto ta_y_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 17); 

                auto ta_y_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 18); 

                auto ta_y_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 19); 

                auto ta_y_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 20); 

                auto ta_y_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 21); 

                auto ta_y_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 22); 

                auto ta_y_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 23); 

                auto ta_y_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 24); 

                auto ta_y_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 25); 

                auto ta_y_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 26); 

                auto ta_y_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 27); 

                auto ta_y_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 28); 

                auto ta_y_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 29); 

                // set up pointers to integrals

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

                // Batch of Integrals (135,180)

                #pragma omp simd aligned(fx, pa_y, pc_y, ta_y_xxxx_1, ta_y_xxxy_1, ta_y_xxxz_1, ta_y_xxyy_1, \
                                         ta_y_xxyz_1, ta_y_xxzz_1, ta_y_xyyy_1, ta_y_xyyz_1, ta_y_xyzz_1, ta_y_xzzz_1, \
                                         ta_y_yyyy_1, ta_y_yyyz_1, ta_y_yyzz_1, ta_y_yzzz_1, ta_y_zzzz_1, tex_0_xxxx_0, \
                                         tex_0_xxxx_1, tex_0_xxxy_0, tex_0_xxxy_1, tex_0_xxxz_0, tex_0_xxxz_1, tex_0_xxyy_0, \
                                         tex_0_xxyy_1, tex_0_xxyz_0, tex_0_xxyz_1, tex_0_xxzz_0, tex_0_xxzz_1, tex_0_xyyy_0, \
                                         tex_0_xyyy_1, tex_0_xyyz_0, tex_0_xyyz_1, tex_0_xyzz_0, tex_0_xyzz_1, tex_0_xzzz_0, \
                                         tex_0_xzzz_1, tex_0_yyyy_0, tex_0_yyyy_1, tex_0_yyyz_0, tex_0_yyyz_1, tex_0_yyzz_0, \
                                         tex_0_yyzz_1, tex_0_yzzz_0, tex_0_yzzz_1, tex_0_zzzz_0, tex_0_zzzz_1, tex_y_xxx_0, \
                                         tex_y_xxx_1, tex_y_xxxx_0, tex_y_xxxx_1, tex_y_xxxy_0, tex_y_xxxy_1, tex_y_xxxz_0, \
                                         tex_y_xxxz_1, tex_y_xxy_0, tex_y_xxy_1, tex_y_xxyy_0, tex_y_xxyy_1, tex_y_xxyz_0, \
                                         tex_y_xxyz_1, tex_y_xxz_0, tex_y_xxz_1, tex_y_xxzz_0, tex_y_xxzz_1, tex_y_xyy_0, \
                                         tex_y_xyy_1, tex_y_xyyy_0, tex_y_xyyy_1, tex_y_xyyz_0, tex_y_xyyz_1, tex_y_xyz_0, \
                                         tex_y_xyz_1, tex_y_xyzz_0, tex_y_xyzz_1, tex_y_xzz_0, tex_y_xzz_1, tex_y_xzzz_0, \
                                         tex_y_xzzz_1, tex_y_yyy_0, tex_y_yyy_1, tex_y_yyyy_0, tex_y_yyyy_1, tex_y_yyyz_0, \
                                         tex_y_yyyz_1, tex_y_yyz_0, tex_y_yyz_1, tex_y_yyzz_0, tex_y_yyzz_1, tex_y_yzz_0, \
                                         tex_y_yzz_1, tex_y_yzzz_0, tex_y_yzzz_1, tex_y_zzz_0, tex_y_zzz_1, tex_y_zzzz_0, \
                                         tex_y_zzzz_1, tex_yy_xxxx_0, tex_yy_xxxy_0, tex_yy_xxxz_0, tex_yy_xxyy_0, \
                                         tex_yy_xxyz_0, tex_yy_xxzz_0, tex_yy_xyyy_0, tex_yy_xyyz_0, tex_yy_xyzz_0, \
                                         tex_yy_xzzz_0, tex_yy_yyyy_0, tex_yy_yyyz_0, tex_yy_yyzz_0, tex_yy_yzzz_0, \
                                         tex_yy_zzzz_0, tey_0_xxxx_0, tey_0_xxxx_1, tey_0_xxxy_0, tey_0_xxxy_1, tey_0_xxxz_0, \
                                         tey_0_xxxz_1, tey_0_xxyy_0, tey_0_xxyy_1, tey_0_xxyz_0, tey_0_xxyz_1, tey_0_xxzz_0, \
                                         tey_0_xxzz_1, tey_0_xyyy_0, tey_0_xyyy_1, tey_0_xyyz_0, tey_0_xyyz_1, tey_0_xyzz_0, \
                                         tey_0_xyzz_1, tey_0_xzzz_0, tey_0_xzzz_1, tey_0_yyyy_0, tey_0_yyyy_1, tey_0_yyyz_0, \
                                         tey_0_yyyz_1, tey_0_yyzz_0, tey_0_yyzz_1, tey_0_yzzz_0, tey_0_yzzz_1, tey_0_zzzz_0, \
                                         tey_0_zzzz_1, tey_y_xxx_0, tey_y_xxx_1, tey_y_xxxx_0, tey_y_xxxx_1, tey_y_xxxy_0, \
                                         tey_y_xxxy_1, tey_y_xxxz_0, tey_y_xxxz_1, tey_y_xxy_0, tey_y_xxy_1, tey_y_xxyy_0, \
                                         tey_y_xxyy_1, tey_y_xxyz_0, tey_y_xxyz_1, tey_y_xxz_0, tey_y_xxz_1, tey_y_xxzz_0, \
                                         tey_y_xxzz_1, tey_y_xyy_0, tey_y_xyy_1, tey_y_xyyy_0, tey_y_xyyy_1, tey_y_xyyz_0, \
                                         tey_y_xyyz_1, tey_y_xyz_0, tey_y_xyz_1, tey_y_xyzz_0, tey_y_xyzz_1, tey_y_xzz_0, \
                                         tey_y_xzz_1, tey_y_xzzz_0, tey_y_xzzz_1, tey_y_yyy_0, tey_y_yyy_1, tey_y_yyyy_0, \
                                         tey_y_yyyy_1, tey_y_yyyz_0, tey_y_yyyz_1, tey_y_yyz_0, tey_y_yyz_1, tey_y_yyzz_0, \
                                         tey_y_yyzz_1, tey_y_yzz_0, tey_y_yzz_1, tey_y_yzzz_0, tey_y_yzzz_1, tey_y_zzz_0, \
                                         tey_y_zzz_1, tey_y_zzzz_0, tey_y_zzzz_1, tey_yy_xxxx_0, tey_yy_xxxy_0, \
                                         tey_yy_xxxz_0, tey_yy_xxyy_0, tey_yy_xxyz_0, tey_yy_xxzz_0, tey_yy_xyyy_0, \
                                         tey_yy_xyyz_0, tey_yy_xyzz_0, tey_yy_xzzz_0, tey_yy_yyyy_0, tey_yy_yyyz_0, \
                                         tey_yy_yyzz_0, tey_yy_yzzz_0, tey_yy_zzzz_0, tez_0_xxxx_0, tez_0_xxxx_1, \
                                         tez_0_xxxy_0, tez_0_xxxy_1, tez_0_xxxz_0, tez_0_xxxz_1, tez_0_xxyy_0, tez_0_xxyy_1, \
                                         tez_0_xxyz_0, tez_0_xxyz_1, tez_0_xxzz_0, tez_0_xxzz_1, tez_0_xyyy_0, tez_0_xyyy_1, \
                                         tez_0_xyyz_0, tez_0_xyyz_1, tez_0_xyzz_0, tez_0_xyzz_1, tez_0_xzzz_0, tez_0_xzzz_1, \
                                         tez_0_yyyy_0, tez_0_yyyy_1, tez_0_yyyz_0, tez_0_yyyz_1, tez_0_yyzz_0, tez_0_yyzz_1, \
                                         tez_0_yzzz_0, tez_0_yzzz_1, tez_0_zzzz_0, tez_0_zzzz_1, tez_y_xxx_0, tez_y_xxx_1, \
                                         tez_y_xxxx_0, tez_y_xxxx_1, tez_y_xxxy_0, tez_y_xxxy_1, tez_y_xxxz_0, tez_y_xxxz_1, \
                                         tez_y_xxy_0, tez_y_xxy_1, tez_y_xxyy_0, tez_y_xxyy_1, tez_y_xxyz_0, tez_y_xxyz_1, \
                                         tez_y_xxz_0, tez_y_xxz_1, tez_y_xxzz_0, tez_y_xxzz_1, tez_y_xyy_0, tez_y_xyy_1, \
                                         tez_y_xyyy_0, tez_y_xyyy_1, tez_y_xyyz_0, tez_y_xyyz_1, tez_y_xyz_0, tez_y_xyz_1, \
                                         tez_y_xyzz_0, tez_y_xyzz_1, tez_y_xzz_0, tez_y_xzz_1, tez_y_xzzz_0, tez_y_xzzz_1, \
                                         tez_y_yyy_0, tez_y_yyy_1, tez_y_yyyy_0, tez_y_yyyy_1, tez_y_yyyz_0, tez_y_yyyz_1, \
                                         tez_y_yyz_0, tez_y_yyz_1, tez_y_yyzz_0, tez_y_yyzz_1, tez_y_yzz_0, tez_y_yzz_1, \
                                         tez_y_yzzz_0, tez_y_yzzz_1, tez_y_zzz_0, tez_y_zzz_1, tez_y_zzzz_0, tez_y_zzzz_1, \
                                         tez_yy_xxxx_0, tez_yy_xxxy_0, tez_yy_xxxz_0, tez_yy_xxyy_0, tez_yy_xxyz_0, \
                                         tez_yy_xxzz_0, tez_yy_xyyy_0, tez_yy_xyyz_0, tez_yy_xyzz_0, tez_yy_xzzz_0, \
                                         tez_yy_yyyy_0, tez_yy_yyyz_0, tez_yy_yyzz_0, tez_yy_yzzz_0, tez_yy_zzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yy_xxxx_0[j] = pa_y[j] * tex_y_xxxx_0[j] - pc_y[j] * tex_y_xxxx_1[j] + 0.5 * fl1_fx * tex_0_xxxx_0[j] - 0.5 * fl1_fx * tex_0_xxxx_1[j];

                    tey_yy_xxxx_0[j] = pa_y[j] * tey_y_xxxx_0[j] - pc_y[j] * tey_y_xxxx_1[j] + 0.5 * fl1_fx * tey_0_xxxx_0[j] - 0.5 * fl1_fx * tey_0_xxxx_1[j] + ta_y_xxxx_1[j];

                    tez_yy_xxxx_0[j] = pa_y[j] * tez_y_xxxx_0[j] - pc_y[j] * tez_y_xxxx_1[j] + 0.5 * fl1_fx * tez_0_xxxx_0[j] - 0.5 * fl1_fx * tez_0_xxxx_1[j];

                    tex_yy_xxxy_0[j] = pa_y[j] * tex_y_xxxy_0[j] - pc_y[j] * tex_y_xxxy_1[j] + 0.5 * fl1_fx * tex_0_xxxy_0[j] - 0.5 * fl1_fx * tex_0_xxxy_1[j] + 0.5 * fl1_fx * tex_y_xxx_0[j] - 0.5 * fl1_fx * tex_y_xxx_1[j];

                    tey_yy_xxxy_0[j] = pa_y[j] * tey_y_xxxy_0[j] - pc_y[j] * tey_y_xxxy_1[j] + 0.5 * fl1_fx * tey_0_xxxy_0[j] - 0.5 * fl1_fx * tey_0_xxxy_1[j] + 0.5 * fl1_fx * tey_y_xxx_0[j] - 0.5 * fl1_fx * tey_y_xxx_1[j] + ta_y_xxxy_1[j];

                    tez_yy_xxxy_0[j] = pa_y[j] * tez_y_xxxy_0[j] - pc_y[j] * tez_y_xxxy_1[j] + 0.5 * fl1_fx * tez_0_xxxy_0[j] - 0.5 * fl1_fx * tez_0_xxxy_1[j] + 0.5 * fl1_fx * tez_y_xxx_0[j] - 0.5 * fl1_fx * tez_y_xxx_1[j];

                    tex_yy_xxxz_0[j] = pa_y[j] * tex_y_xxxz_0[j] - pc_y[j] * tex_y_xxxz_1[j] + 0.5 * fl1_fx * tex_0_xxxz_0[j] - 0.5 * fl1_fx * tex_0_xxxz_1[j];

                    tey_yy_xxxz_0[j] = pa_y[j] * tey_y_xxxz_0[j] - pc_y[j] * tey_y_xxxz_1[j] + 0.5 * fl1_fx * tey_0_xxxz_0[j] - 0.5 * fl1_fx * tey_0_xxxz_1[j] + ta_y_xxxz_1[j];

                    tez_yy_xxxz_0[j] = pa_y[j] * tez_y_xxxz_0[j] - pc_y[j] * tez_y_xxxz_1[j] + 0.5 * fl1_fx * tez_0_xxxz_0[j] - 0.5 * fl1_fx * tez_0_xxxz_1[j];

                    tex_yy_xxyy_0[j] = pa_y[j] * tex_y_xxyy_0[j] - pc_y[j] * tex_y_xxyy_1[j] + 0.5 * fl1_fx * tex_0_xxyy_0[j] - 0.5 * fl1_fx * tex_0_xxyy_1[j] + fl1_fx * tex_y_xxy_0[j] - fl1_fx * tex_y_xxy_1[j];

                    tey_yy_xxyy_0[j] = pa_y[j] * tey_y_xxyy_0[j] - pc_y[j] * tey_y_xxyy_1[j] + 0.5 * fl1_fx * tey_0_xxyy_0[j] - 0.5 * fl1_fx * tey_0_xxyy_1[j] + fl1_fx * tey_y_xxy_0[j] - fl1_fx * tey_y_xxy_1[j] + ta_y_xxyy_1[j];

                    tez_yy_xxyy_0[j] = pa_y[j] * tez_y_xxyy_0[j] - pc_y[j] * tez_y_xxyy_1[j] + 0.5 * fl1_fx * tez_0_xxyy_0[j] - 0.5 * fl1_fx * tez_0_xxyy_1[j] + fl1_fx * tez_y_xxy_0[j] - fl1_fx * tez_y_xxy_1[j];

                    tex_yy_xxyz_0[j] = pa_y[j] * tex_y_xxyz_0[j] - pc_y[j] * tex_y_xxyz_1[j] + 0.5 * fl1_fx * tex_0_xxyz_0[j] - 0.5 * fl1_fx * tex_0_xxyz_1[j] + 0.5 * fl1_fx * tex_y_xxz_0[j] - 0.5 * fl1_fx * tex_y_xxz_1[j];

                    tey_yy_xxyz_0[j] = pa_y[j] * tey_y_xxyz_0[j] - pc_y[j] * tey_y_xxyz_1[j] + 0.5 * fl1_fx * tey_0_xxyz_0[j] - 0.5 * fl1_fx * tey_0_xxyz_1[j] + 0.5 * fl1_fx * tey_y_xxz_0[j] - 0.5 * fl1_fx * tey_y_xxz_1[j] + ta_y_xxyz_1[j];

                    tez_yy_xxyz_0[j] = pa_y[j] * tez_y_xxyz_0[j] - pc_y[j] * tez_y_xxyz_1[j] + 0.5 * fl1_fx * tez_0_xxyz_0[j] - 0.5 * fl1_fx * tez_0_xxyz_1[j] + 0.5 * fl1_fx * tez_y_xxz_0[j] - 0.5 * fl1_fx * tez_y_xxz_1[j];

                    tex_yy_xxzz_0[j] = pa_y[j] * tex_y_xxzz_0[j] - pc_y[j] * tex_y_xxzz_1[j] + 0.5 * fl1_fx * tex_0_xxzz_0[j] - 0.5 * fl1_fx * tex_0_xxzz_1[j];

                    tey_yy_xxzz_0[j] = pa_y[j] * tey_y_xxzz_0[j] - pc_y[j] * tey_y_xxzz_1[j] + 0.5 * fl1_fx * tey_0_xxzz_0[j] - 0.5 * fl1_fx * tey_0_xxzz_1[j] + ta_y_xxzz_1[j];

                    tez_yy_xxzz_0[j] = pa_y[j] * tez_y_xxzz_0[j] - pc_y[j] * tez_y_xxzz_1[j] + 0.5 * fl1_fx * tez_0_xxzz_0[j] - 0.5 * fl1_fx * tez_0_xxzz_1[j];

                    tex_yy_xyyy_0[j] = pa_y[j] * tex_y_xyyy_0[j] - pc_y[j] * tex_y_xyyy_1[j] + 0.5 * fl1_fx * tex_0_xyyy_0[j] - 0.5 * fl1_fx * tex_0_xyyy_1[j] + 1.5 * fl1_fx * tex_y_xyy_0[j] - 1.5 * fl1_fx * tex_y_xyy_1[j];

                    tey_yy_xyyy_0[j] = pa_y[j] * tey_y_xyyy_0[j] - pc_y[j] * tey_y_xyyy_1[j] + 0.5 * fl1_fx * tey_0_xyyy_0[j] - 0.5 * fl1_fx * tey_0_xyyy_1[j] + 1.5 * fl1_fx * tey_y_xyy_0[j] - 1.5 * fl1_fx * tey_y_xyy_1[j] + ta_y_xyyy_1[j];

                    tez_yy_xyyy_0[j] = pa_y[j] * tez_y_xyyy_0[j] - pc_y[j] * tez_y_xyyy_1[j] + 0.5 * fl1_fx * tez_0_xyyy_0[j] - 0.5 * fl1_fx * tez_0_xyyy_1[j] + 1.5 * fl1_fx * tez_y_xyy_0[j] - 1.5 * fl1_fx * tez_y_xyy_1[j];

                    tex_yy_xyyz_0[j] = pa_y[j] * tex_y_xyyz_0[j] - pc_y[j] * tex_y_xyyz_1[j] + 0.5 * fl1_fx * tex_0_xyyz_0[j] - 0.5 * fl1_fx * tex_0_xyyz_1[j] + fl1_fx * tex_y_xyz_0[j] - fl1_fx * tex_y_xyz_1[j];

                    tey_yy_xyyz_0[j] = pa_y[j] * tey_y_xyyz_0[j] - pc_y[j] * tey_y_xyyz_1[j] + 0.5 * fl1_fx * tey_0_xyyz_0[j] - 0.5 * fl1_fx * tey_0_xyyz_1[j] + fl1_fx * tey_y_xyz_0[j] - fl1_fx * tey_y_xyz_1[j] + ta_y_xyyz_1[j];

                    tez_yy_xyyz_0[j] = pa_y[j] * tez_y_xyyz_0[j] - pc_y[j] * tez_y_xyyz_1[j] + 0.5 * fl1_fx * tez_0_xyyz_0[j] - 0.5 * fl1_fx * tez_0_xyyz_1[j] + fl1_fx * tez_y_xyz_0[j] - fl1_fx * tez_y_xyz_1[j];

                    tex_yy_xyzz_0[j] = pa_y[j] * tex_y_xyzz_0[j] - pc_y[j] * tex_y_xyzz_1[j] + 0.5 * fl1_fx * tex_0_xyzz_0[j] - 0.5 * fl1_fx * tex_0_xyzz_1[j] + 0.5 * fl1_fx * tex_y_xzz_0[j] - 0.5 * fl1_fx * tex_y_xzz_1[j];

                    tey_yy_xyzz_0[j] = pa_y[j] * tey_y_xyzz_0[j] - pc_y[j] * tey_y_xyzz_1[j] + 0.5 * fl1_fx * tey_0_xyzz_0[j] - 0.5 * fl1_fx * tey_0_xyzz_1[j] + 0.5 * fl1_fx * tey_y_xzz_0[j] - 0.5 * fl1_fx * tey_y_xzz_1[j] + ta_y_xyzz_1[j];

                    tez_yy_xyzz_0[j] = pa_y[j] * tez_y_xyzz_0[j] - pc_y[j] * tez_y_xyzz_1[j] + 0.5 * fl1_fx * tez_0_xyzz_0[j] - 0.5 * fl1_fx * tez_0_xyzz_1[j] + 0.5 * fl1_fx * tez_y_xzz_0[j] - 0.5 * fl1_fx * tez_y_xzz_1[j];

                    tex_yy_xzzz_0[j] = pa_y[j] * tex_y_xzzz_0[j] - pc_y[j] * tex_y_xzzz_1[j] + 0.5 * fl1_fx * tex_0_xzzz_0[j] - 0.5 * fl1_fx * tex_0_xzzz_1[j];

                    tey_yy_xzzz_0[j] = pa_y[j] * tey_y_xzzz_0[j] - pc_y[j] * tey_y_xzzz_1[j] + 0.5 * fl1_fx * tey_0_xzzz_0[j] - 0.5 * fl1_fx * tey_0_xzzz_1[j] + ta_y_xzzz_1[j];

                    tez_yy_xzzz_0[j] = pa_y[j] * tez_y_xzzz_0[j] - pc_y[j] * tez_y_xzzz_1[j] + 0.5 * fl1_fx * tez_0_xzzz_0[j] - 0.5 * fl1_fx * tez_0_xzzz_1[j];

                    tex_yy_yyyy_0[j] = pa_y[j] * tex_y_yyyy_0[j] - pc_y[j] * tex_y_yyyy_1[j] + 0.5 * fl1_fx * tex_0_yyyy_0[j] - 0.5 * fl1_fx * tex_0_yyyy_1[j] + 2.0 * fl1_fx * tex_y_yyy_0[j] - 2.0 * fl1_fx * tex_y_yyy_1[j];

                    tey_yy_yyyy_0[j] = pa_y[j] * tey_y_yyyy_0[j] - pc_y[j] * tey_y_yyyy_1[j] + 0.5 * fl1_fx * tey_0_yyyy_0[j] - 0.5 * fl1_fx * tey_0_yyyy_1[j] + 2.0 * fl1_fx * tey_y_yyy_0[j] - 2.0 * fl1_fx * tey_y_yyy_1[j] + ta_y_yyyy_1[j];

                    tez_yy_yyyy_0[j] = pa_y[j] * tez_y_yyyy_0[j] - pc_y[j] * tez_y_yyyy_1[j] + 0.5 * fl1_fx * tez_0_yyyy_0[j] - 0.5 * fl1_fx * tez_0_yyyy_1[j] + 2.0 * fl1_fx * tez_y_yyy_0[j] - 2.0 * fl1_fx * tez_y_yyy_1[j];

                    tex_yy_yyyz_0[j] = pa_y[j] * tex_y_yyyz_0[j] - pc_y[j] * tex_y_yyyz_1[j] + 0.5 * fl1_fx * tex_0_yyyz_0[j] - 0.5 * fl1_fx * tex_0_yyyz_1[j] + 1.5 * fl1_fx * tex_y_yyz_0[j] - 1.5 * fl1_fx * tex_y_yyz_1[j];

                    tey_yy_yyyz_0[j] = pa_y[j] * tey_y_yyyz_0[j] - pc_y[j] * tey_y_yyyz_1[j] + 0.5 * fl1_fx * tey_0_yyyz_0[j] - 0.5 * fl1_fx * tey_0_yyyz_1[j] + 1.5 * fl1_fx * tey_y_yyz_0[j] - 1.5 * fl1_fx * tey_y_yyz_1[j] + ta_y_yyyz_1[j];

                    tez_yy_yyyz_0[j] = pa_y[j] * tez_y_yyyz_0[j] - pc_y[j] * tez_y_yyyz_1[j] + 0.5 * fl1_fx * tez_0_yyyz_0[j] - 0.5 * fl1_fx * tez_0_yyyz_1[j] + 1.5 * fl1_fx * tez_y_yyz_0[j] - 1.5 * fl1_fx * tez_y_yyz_1[j];

                    tex_yy_yyzz_0[j] = pa_y[j] * tex_y_yyzz_0[j] - pc_y[j] * tex_y_yyzz_1[j] + 0.5 * fl1_fx * tex_0_yyzz_0[j] - 0.5 * fl1_fx * tex_0_yyzz_1[j] + fl1_fx * tex_y_yzz_0[j] - fl1_fx * tex_y_yzz_1[j];

                    tey_yy_yyzz_0[j] = pa_y[j] * tey_y_yyzz_0[j] - pc_y[j] * tey_y_yyzz_1[j] + 0.5 * fl1_fx * tey_0_yyzz_0[j] - 0.5 * fl1_fx * tey_0_yyzz_1[j] + fl1_fx * tey_y_yzz_0[j] - fl1_fx * tey_y_yzz_1[j] + ta_y_yyzz_1[j];

                    tez_yy_yyzz_0[j] = pa_y[j] * tez_y_yyzz_0[j] - pc_y[j] * tez_y_yyzz_1[j] + 0.5 * fl1_fx * tez_0_yyzz_0[j] - 0.5 * fl1_fx * tez_0_yyzz_1[j] + fl1_fx * tez_y_yzz_0[j] - fl1_fx * tez_y_yzz_1[j];

                    tex_yy_yzzz_0[j] = pa_y[j] * tex_y_yzzz_0[j] - pc_y[j] * tex_y_yzzz_1[j] + 0.5 * fl1_fx * tex_0_yzzz_0[j] - 0.5 * fl1_fx * tex_0_yzzz_1[j] + 0.5 * fl1_fx * tex_y_zzz_0[j] - 0.5 * fl1_fx * tex_y_zzz_1[j];

                    tey_yy_yzzz_0[j] = pa_y[j] * tey_y_yzzz_0[j] - pc_y[j] * tey_y_yzzz_1[j] + 0.5 * fl1_fx * tey_0_yzzz_0[j] - 0.5 * fl1_fx * tey_0_yzzz_1[j] + 0.5 * fl1_fx * tey_y_zzz_0[j] - 0.5 * fl1_fx * tey_y_zzz_1[j] + ta_y_yzzz_1[j];

                    tez_yy_yzzz_0[j] = pa_y[j] * tez_y_yzzz_0[j] - pc_y[j] * tez_y_yzzz_1[j] + 0.5 * fl1_fx * tez_0_yzzz_0[j] - 0.5 * fl1_fx * tez_0_yzzz_1[j] + 0.5 * fl1_fx * tez_y_zzz_0[j] - 0.5 * fl1_fx * tez_y_zzz_1[j];

                    tex_yy_zzzz_0[j] = pa_y[j] * tex_y_zzzz_0[j] - pc_y[j] * tex_y_zzzz_1[j] + 0.5 * fl1_fx * tex_0_zzzz_0[j] - 0.5 * fl1_fx * tex_0_zzzz_1[j];

                    tey_yy_zzzz_0[j] = pa_y[j] * tey_y_zzzz_0[j] - pc_y[j] * tey_y_zzzz_1[j] + 0.5 * fl1_fx * tey_0_zzzz_0[j] - 0.5 * fl1_fx * tey_0_zzzz_1[j] + ta_y_zzzz_1[j];

                    tez_yy_zzzz_0[j] = pa_y[j] * tez_y_zzzz_0[j] - pc_y[j] * tez_y_zzzz_1[j] + 0.5 * fl1_fx * tez_0_zzzz_0[j] - 0.5 * fl1_fx * tez_0_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDG_180_225(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (180,225)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
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

                auto tex_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 20); 

                auto tey_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 21); 

                auto tey_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 22); 

                auto tey_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 23); 

                auto tey_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 24); 

                auto tey_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 25); 

                auto tey_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 26); 

                auto tey_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 27); 

                auto tey_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 28); 

                auto tey_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 29); 

                auto tey_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 29); 

                auto tex_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 20); 

                auto tey_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 21); 

                auto tey_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 22); 

                auto tey_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 23); 

                auto tey_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 24); 

                auto tey_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 25); 

                auto tey_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 26); 

                auto tey_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 27); 

                auto tey_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 28); 

                auto tey_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 29); 

                auto tey_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 29); 

                auto ta_z_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 30); 

                auto ta_z_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 31); 

                auto ta_z_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 32); 

                auto ta_z_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 33); 

                auto ta_z_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 34); 

                auto ta_z_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 35); 

                auto ta_z_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 36); 

                auto ta_z_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 37); 

                auto ta_z_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 38); 

                auto ta_z_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 39); 

                auto ta_z_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 40); 

                auto ta_z_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 41); 

                auto ta_z_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 42); 

                auto ta_z_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 43); 

                auto ta_z_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 44); 

                // set up pointers to integrals

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

                // Batch of Integrals (180,225)

                #pragma omp simd aligned(fx, pa_y, pc_y, ta_z_xxxx_1, ta_z_xxxy_1, ta_z_xxxz_1, ta_z_xxyy_1, \
                                         ta_z_xxyz_1, ta_z_xxzz_1, ta_z_xyyy_1, ta_z_xyyz_1, ta_z_xyzz_1, ta_z_xzzz_1, \
                                         ta_z_yyyy_1, ta_z_yyyz_1, ta_z_yyzz_1, ta_z_yzzz_1, ta_z_zzzz_1, tex_yz_xxxx_0, \
                                         tex_yz_xxxy_0, tex_yz_xxxz_0, tex_yz_xxyy_0, tex_yz_xxyz_0, tex_yz_xxzz_0, \
                                         tex_yz_xyyy_0, tex_yz_xyyz_0, tex_yz_xyzz_0, tex_yz_xzzz_0, tex_yz_yyyy_0, \
                                         tex_yz_yyyz_0, tex_yz_yyzz_0, tex_yz_yzzz_0, tex_yz_zzzz_0, tex_z_xxx_0, \
                                         tex_z_xxx_1, tex_z_xxxx_0, tex_z_xxxx_1, tex_z_xxxy_0, tex_z_xxxy_1, tex_z_xxxz_0, \
                                         tex_z_xxxz_1, tex_z_xxy_0, tex_z_xxy_1, tex_z_xxyy_0, tex_z_xxyy_1, tex_z_xxyz_0, \
                                         tex_z_xxyz_1, tex_z_xxz_0, tex_z_xxz_1, tex_z_xxzz_0, tex_z_xxzz_1, tex_z_xyy_0, \
                                         tex_z_xyy_1, tex_z_xyyy_0, tex_z_xyyy_1, tex_z_xyyz_0, tex_z_xyyz_1, tex_z_xyz_0, \
                                         tex_z_xyz_1, tex_z_xyzz_0, tex_z_xyzz_1, tex_z_xzz_0, tex_z_xzz_1, tex_z_xzzz_0, \
                                         tex_z_xzzz_1, tex_z_yyy_0, tex_z_yyy_1, tex_z_yyyy_0, tex_z_yyyy_1, tex_z_yyyz_0, \
                                         tex_z_yyyz_1, tex_z_yyz_0, tex_z_yyz_1, tex_z_yyzz_0, tex_z_yyzz_1, tex_z_yzz_0, \
                                         tex_z_yzz_1, tex_z_yzzz_0, tex_z_yzzz_1, tex_z_zzz_0, tex_z_zzz_1, tex_z_zzzz_0, \
                                         tex_z_zzzz_1, tey_yz_xxxx_0, tey_yz_xxxy_0, tey_yz_xxxz_0, tey_yz_xxyy_0, \
                                         tey_yz_xxyz_0, tey_yz_xxzz_0, tey_yz_xyyy_0, tey_yz_xyyz_0, tey_yz_xyzz_0, \
                                         tey_yz_xzzz_0, tey_yz_yyyy_0, tey_yz_yyyz_0, tey_yz_yyzz_0, tey_yz_yzzz_0, \
                                         tey_yz_zzzz_0, tey_z_xxx_0, tey_z_xxx_1, tey_z_xxxx_0, tey_z_xxxx_1, tey_z_xxxy_0, \
                                         tey_z_xxxy_1, tey_z_xxxz_0, tey_z_xxxz_1, tey_z_xxy_0, tey_z_xxy_1, tey_z_xxyy_0, \
                                         tey_z_xxyy_1, tey_z_xxyz_0, tey_z_xxyz_1, tey_z_xxz_0, tey_z_xxz_1, tey_z_xxzz_0, \
                                         tey_z_xxzz_1, tey_z_xyy_0, tey_z_xyy_1, tey_z_xyyy_0, tey_z_xyyy_1, tey_z_xyyz_0, \
                                         tey_z_xyyz_1, tey_z_xyz_0, tey_z_xyz_1, tey_z_xyzz_0, tey_z_xyzz_1, tey_z_xzz_0, \
                                         tey_z_xzz_1, tey_z_xzzz_0, tey_z_xzzz_1, tey_z_yyy_0, tey_z_yyy_1, tey_z_yyyy_0, \
                                         tey_z_yyyy_1, tey_z_yyyz_0, tey_z_yyyz_1, tey_z_yyz_0, tey_z_yyz_1, tey_z_yyzz_0, \
                                         tey_z_yyzz_1, tey_z_yzz_0, tey_z_yzz_1, tey_z_yzzz_0, tey_z_yzzz_1, tey_z_zzz_0, \
                                         tey_z_zzz_1, tey_z_zzzz_0, tey_z_zzzz_1, tez_yz_xxxx_0, tez_yz_xxxy_0, \
                                         tez_yz_xxxz_0, tez_yz_xxyy_0, tez_yz_xxyz_0, tez_yz_xxzz_0, tez_yz_xyyy_0, \
                                         tez_yz_xyyz_0, tez_yz_xyzz_0, tez_yz_xzzz_0, tez_yz_yyyy_0, tez_yz_yyyz_0, \
                                         tez_yz_yyzz_0, tez_yz_yzzz_0, tez_yz_zzzz_0, tez_z_xxx_0, tez_z_xxx_1, tez_z_xxxx_0, \
                                         tez_z_xxxx_1, tez_z_xxxy_0, tez_z_xxxy_1, tez_z_xxxz_0, tez_z_xxxz_1, tez_z_xxy_0, \
                                         tez_z_xxy_1, tez_z_xxyy_0, tez_z_xxyy_1, tez_z_xxyz_0, tez_z_xxyz_1, tez_z_xxz_0, \
                                         tez_z_xxz_1, tez_z_xxzz_0, tez_z_xxzz_1, tez_z_xyy_0, tez_z_xyy_1, tez_z_xyyy_0, \
                                         tez_z_xyyy_1, tez_z_xyyz_0, tez_z_xyyz_1, tez_z_xyz_0, tez_z_xyz_1, tez_z_xyzz_0, \
                                         tez_z_xyzz_1, tez_z_xzz_0, tez_z_xzz_1, tez_z_xzzz_0, tez_z_xzzz_1, tez_z_yyy_0, \
                                         tez_z_yyy_1, tez_z_yyyy_0, tez_z_yyyy_1, tez_z_yyyz_0, tez_z_yyyz_1, tez_z_yyz_0, \
                                         tez_z_yyz_1, tez_z_yyzz_0, tez_z_yyzz_1, tez_z_yzz_0, tez_z_yzz_1, tez_z_yzzz_0, \
                                         tez_z_yzzz_1, tez_z_zzz_0, tez_z_zzz_1, tez_z_zzzz_0, tez_z_zzzz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yz_xxxx_0[j] = pa_y[j] * tex_z_xxxx_0[j] - pc_y[j] * tex_z_xxxx_1[j];

                    tey_yz_xxxx_0[j] = pa_y[j] * tey_z_xxxx_0[j] - pc_y[j] * tey_z_xxxx_1[j] + ta_z_xxxx_1[j];

                    tez_yz_xxxx_0[j] = pa_y[j] * tez_z_xxxx_0[j] - pc_y[j] * tez_z_xxxx_1[j];

                    tex_yz_xxxy_0[j] = pa_y[j] * tex_z_xxxy_0[j] - pc_y[j] * tex_z_xxxy_1[j] + 0.5 * fl1_fx * tex_z_xxx_0[j] - 0.5 * fl1_fx * tex_z_xxx_1[j];

                    tey_yz_xxxy_0[j] = pa_y[j] * tey_z_xxxy_0[j] - pc_y[j] * tey_z_xxxy_1[j] + 0.5 * fl1_fx * tey_z_xxx_0[j] - 0.5 * fl1_fx * tey_z_xxx_1[j] + ta_z_xxxy_1[j];

                    tez_yz_xxxy_0[j] = pa_y[j] * tez_z_xxxy_0[j] - pc_y[j] * tez_z_xxxy_1[j] + 0.5 * fl1_fx * tez_z_xxx_0[j] - 0.5 * fl1_fx * tez_z_xxx_1[j];

                    tex_yz_xxxz_0[j] = pa_y[j] * tex_z_xxxz_0[j] - pc_y[j] * tex_z_xxxz_1[j];

                    tey_yz_xxxz_0[j] = pa_y[j] * tey_z_xxxz_0[j] - pc_y[j] * tey_z_xxxz_1[j] + ta_z_xxxz_1[j];

                    tez_yz_xxxz_0[j] = pa_y[j] * tez_z_xxxz_0[j] - pc_y[j] * tez_z_xxxz_1[j];

                    tex_yz_xxyy_0[j] = pa_y[j] * tex_z_xxyy_0[j] - pc_y[j] * tex_z_xxyy_1[j] + fl1_fx * tex_z_xxy_0[j] - fl1_fx * tex_z_xxy_1[j];

                    tey_yz_xxyy_0[j] = pa_y[j] * tey_z_xxyy_0[j] - pc_y[j] * tey_z_xxyy_1[j] + fl1_fx * tey_z_xxy_0[j] - fl1_fx * tey_z_xxy_1[j] + ta_z_xxyy_1[j];

                    tez_yz_xxyy_0[j] = pa_y[j] * tez_z_xxyy_0[j] - pc_y[j] * tez_z_xxyy_1[j] + fl1_fx * tez_z_xxy_0[j] - fl1_fx * tez_z_xxy_1[j];

                    tex_yz_xxyz_0[j] = pa_y[j] * tex_z_xxyz_0[j] - pc_y[j] * tex_z_xxyz_1[j] + 0.5 * fl1_fx * tex_z_xxz_0[j] - 0.5 * fl1_fx * tex_z_xxz_1[j];

                    tey_yz_xxyz_0[j] = pa_y[j] * tey_z_xxyz_0[j] - pc_y[j] * tey_z_xxyz_1[j] + 0.5 * fl1_fx * tey_z_xxz_0[j] - 0.5 * fl1_fx * tey_z_xxz_1[j] + ta_z_xxyz_1[j];

                    tez_yz_xxyz_0[j] = pa_y[j] * tez_z_xxyz_0[j] - pc_y[j] * tez_z_xxyz_1[j] + 0.5 * fl1_fx * tez_z_xxz_0[j] - 0.5 * fl1_fx * tez_z_xxz_1[j];

                    tex_yz_xxzz_0[j] = pa_y[j] * tex_z_xxzz_0[j] - pc_y[j] * tex_z_xxzz_1[j];

                    tey_yz_xxzz_0[j] = pa_y[j] * tey_z_xxzz_0[j] - pc_y[j] * tey_z_xxzz_1[j] + ta_z_xxzz_1[j];

                    tez_yz_xxzz_0[j] = pa_y[j] * tez_z_xxzz_0[j] - pc_y[j] * tez_z_xxzz_1[j];

                    tex_yz_xyyy_0[j] = pa_y[j] * tex_z_xyyy_0[j] - pc_y[j] * tex_z_xyyy_1[j] + 1.5 * fl1_fx * tex_z_xyy_0[j] - 1.5 * fl1_fx * tex_z_xyy_1[j];

                    tey_yz_xyyy_0[j] = pa_y[j] * tey_z_xyyy_0[j] - pc_y[j] * tey_z_xyyy_1[j] + 1.5 * fl1_fx * tey_z_xyy_0[j] - 1.5 * fl1_fx * tey_z_xyy_1[j] + ta_z_xyyy_1[j];

                    tez_yz_xyyy_0[j] = pa_y[j] * tez_z_xyyy_0[j] - pc_y[j] * tez_z_xyyy_1[j] + 1.5 * fl1_fx * tez_z_xyy_0[j] - 1.5 * fl1_fx * tez_z_xyy_1[j];

                    tex_yz_xyyz_0[j] = pa_y[j] * tex_z_xyyz_0[j] - pc_y[j] * tex_z_xyyz_1[j] + fl1_fx * tex_z_xyz_0[j] - fl1_fx * tex_z_xyz_1[j];

                    tey_yz_xyyz_0[j] = pa_y[j] * tey_z_xyyz_0[j] - pc_y[j] * tey_z_xyyz_1[j] + fl1_fx * tey_z_xyz_0[j] - fl1_fx * tey_z_xyz_1[j] + ta_z_xyyz_1[j];

                    tez_yz_xyyz_0[j] = pa_y[j] * tez_z_xyyz_0[j] - pc_y[j] * tez_z_xyyz_1[j] + fl1_fx * tez_z_xyz_0[j] - fl1_fx * tez_z_xyz_1[j];

                    tex_yz_xyzz_0[j] = pa_y[j] * tex_z_xyzz_0[j] - pc_y[j] * tex_z_xyzz_1[j] + 0.5 * fl1_fx * tex_z_xzz_0[j] - 0.5 * fl1_fx * tex_z_xzz_1[j];

                    tey_yz_xyzz_0[j] = pa_y[j] * tey_z_xyzz_0[j] - pc_y[j] * tey_z_xyzz_1[j] + 0.5 * fl1_fx * tey_z_xzz_0[j] - 0.5 * fl1_fx * tey_z_xzz_1[j] + ta_z_xyzz_1[j];

                    tez_yz_xyzz_0[j] = pa_y[j] * tez_z_xyzz_0[j] - pc_y[j] * tez_z_xyzz_1[j] + 0.5 * fl1_fx * tez_z_xzz_0[j] - 0.5 * fl1_fx * tez_z_xzz_1[j];

                    tex_yz_xzzz_0[j] = pa_y[j] * tex_z_xzzz_0[j] - pc_y[j] * tex_z_xzzz_1[j];

                    tey_yz_xzzz_0[j] = pa_y[j] * tey_z_xzzz_0[j] - pc_y[j] * tey_z_xzzz_1[j] + ta_z_xzzz_1[j];

                    tez_yz_xzzz_0[j] = pa_y[j] * tez_z_xzzz_0[j] - pc_y[j] * tez_z_xzzz_1[j];

                    tex_yz_yyyy_0[j] = pa_y[j] * tex_z_yyyy_0[j] - pc_y[j] * tex_z_yyyy_1[j] + 2.0 * fl1_fx * tex_z_yyy_0[j] - 2.0 * fl1_fx * tex_z_yyy_1[j];

                    tey_yz_yyyy_0[j] = pa_y[j] * tey_z_yyyy_0[j] - pc_y[j] * tey_z_yyyy_1[j] + 2.0 * fl1_fx * tey_z_yyy_0[j] - 2.0 * fl1_fx * tey_z_yyy_1[j] + ta_z_yyyy_1[j];

                    tez_yz_yyyy_0[j] = pa_y[j] * tez_z_yyyy_0[j] - pc_y[j] * tez_z_yyyy_1[j] + 2.0 * fl1_fx * tez_z_yyy_0[j] - 2.0 * fl1_fx * tez_z_yyy_1[j];

                    tex_yz_yyyz_0[j] = pa_y[j] * tex_z_yyyz_0[j] - pc_y[j] * tex_z_yyyz_1[j] + 1.5 * fl1_fx * tex_z_yyz_0[j] - 1.5 * fl1_fx * tex_z_yyz_1[j];

                    tey_yz_yyyz_0[j] = pa_y[j] * tey_z_yyyz_0[j] - pc_y[j] * tey_z_yyyz_1[j] + 1.5 * fl1_fx * tey_z_yyz_0[j] - 1.5 * fl1_fx * tey_z_yyz_1[j] + ta_z_yyyz_1[j];

                    tez_yz_yyyz_0[j] = pa_y[j] * tez_z_yyyz_0[j] - pc_y[j] * tez_z_yyyz_1[j] + 1.5 * fl1_fx * tez_z_yyz_0[j] - 1.5 * fl1_fx * tez_z_yyz_1[j];

                    tex_yz_yyzz_0[j] = pa_y[j] * tex_z_yyzz_0[j] - pc_y[j] * tex_z_yyzz_1[j] + fl1_fx * tex_z_yzz_0[j] - fl1_fx * tex_z_yzz_1[j];

                    tey_yz_yyzz_0[j] = pa_y[j] * tey_z_yyzz_0[j] - pc_y[j] * tey_z_yyzz_1[j] + fl1_fx * tey_z_yzz_0[j] - fl1_fx * tey_z_yzz_1[j] + ta_z_yyzz_1[j];

                    tez_yz_yyzz_0[j] = pa_y[j] * tez_z_yyzz_0[j] - pc_y[j] * tez_z_yyzz_1[j] + fl1_fx * tez_z_yzz_0[j] - fl1_fx * tez_z_yzz_1[j];

                    tex_yz_yzzz_0[j] = pa_y[j] * tex_z_yzzz_0[j] - pc_y[j] * tex_z_yzzz_1[j] + 0.5 * fl1_fx * tex_z_zzz_0[j] - 0.5 * fl1_fx * tex_z_zzz_1[j];

                    tey_yz_yzzz_0[j] = pa_y[j] * tey_z_yzzz_0[j] - pc_y[j] * tey_z_yzzz_1[j] + 0.5 * fl1_fx * tey_z_zzz_0[j] - 0.5 * fl1_fx * tey_z_zzz_1[j] + ta_z_yzzz_1[j];

                    tez_yz_yzzz_0[j] = pa_y[j] * tez_z_yzzz_0[j] - pc_y[j] * tez_z_yzzz_1[j] + 0.5 * fl1_fx * tez_z_zzz_0[j] - 0.5 * fl1_fx * tez_z_zzz_1[j];

                    tex_yz_zzzz_0[j] = pa_y[j] * tex_z_zzzz_0[j] - pc_y[j] * tex_z_zzzz_1[j];

                    tey_yz_zzzz_0[j] = pa_y[j] * tey_z_zzzz_0[j] - pc_y[j] * tey_z_zzzz_1[j] + ta_z_zzzz_1[j];

                    tez_yz_zzzz_0[j] = pa_y[j] * tez_z_zzzz_0[j] - pc_y[j] * tez_z_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForDG_225_270(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (225,270)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {2, -1, -1, -1},
                                             {4, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_2_4_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_0_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_1_4_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_1_3_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            // loop over contracted GTO on bra side

            int32_t idx = 0;

            for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
            {
                // set up pointers to Obara-Saika factors

                auto fx = osFactors.data(3 * idx);

                // set up pointers to tensors product of distances R(PA) = P - A

                auto pa_z = paDistances.data(3 * idx + 2);

                // set up pointers to tensors product of distances R(PC) = P - C

                auto pc_z = pcDistances.data(3 * idx + 2);

                // set up pointers to auxilary integrals

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

                auto tex_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx); 

                auto tey_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx); 

                auto tez_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx); 

                auto tex_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 1); 

                auto tey_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 1); 

                auto tez_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 1); 

                auto tex_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 2); 

                auto tey_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 2); 

                auto tez_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 2); 

                auto tex_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 3); 

                auto tey_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 3); 

                auto tez_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 3); 

                auto tex_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 4); 

                auto tey_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 4); 

                auto tez_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 4); 

                auto tex_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 5); 

                auto tey_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 5); 

                auto tez_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 5); 

                auto tex_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 6); 

                auto tey_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 6); 

                auto tez_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 6); 

                auto tex_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 7); 

                auto tey_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 7); 

                auto tez_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 7); 

                auto tex_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 8); 

                auto tey_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 8); 

                auto tez_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 8); 

                auto tex_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 9); 

                auto tey_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 9); 

                auto tez_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 9); 

                auto tex_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 10); 

                auto tey_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 10); 

                auto tez_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 10); 

                auto tex_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 11); 

                auto tey_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 11); 

                auto tez_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 11); 

                auto tex_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 12); 

                auto tey_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 12); 

                auto tez_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 12); 

                auto tex_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 13); 

                auto tey_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 13); 

                auto tez_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 13); 

                auto tex_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 14); 

                auto tey_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 14); 

                auto tez_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 14); 

                auto tex_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx); 

                auto tey_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx); 

                auto tez_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx); 

                auto tex_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 1); 

                auto tey_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 1); 

                auto tez_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 1); 

                auto tex_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 2); 

                auto tey_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 2); 

                auto tez_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 2); 

                auto tex_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 3); 

                auto tey_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 3); 

                auto tez_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 3); 

                auto tex_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 4); 

                auto tey_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 4); 

                auto tez_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 4); 

                auto tex_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 5); 

                auto tey_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 5); 

                auto tez_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 5); 

                auto tex_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 6); 

                auto tey_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 6); 

                auto tez_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 6); 

                auto tex_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 7); 

                auto tey_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 7); 

                auto tez_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 7); 

                auto tex_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 8); 

                auto tey_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 8); 

                auto tez_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 8); 

                auto tex_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 9); 

                auto tey_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 9); 

                auto tez_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 9); 

                auto tex_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 10); 

                auto tey_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 10); 

                auto tez_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 10); 

                auto tex_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 11); 

                auto tey_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 11); 

                auto tez_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 11); 

                auto tex_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 12); 

                auto tey_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 12); 

                auto tez_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 12); 

                auto tex_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 13); 

                auto tey_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 13); 

                auto tez_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 13); 

                auto tex_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 14); 

                auto tey_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 14); 

                auto tez_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 14); 

                auto tex_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 20); 

                auto tey_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 21); 

                auto tey_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 22); 

                auto tey_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 23); 

                auto tey_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 24); 

                auto tey_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 25); 

                auto tey_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 26); 

                auto tey_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 27); 

                auto tey_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 28); 

                auto tey_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 29); 

                auto tey_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 29); 

                auto tex_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 20); 

                auto tey_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_z_xxx_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 21); 

                auto tey_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_z_xxy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 22); 

                auto tey_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_z_xxz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 23); 

                auto tey_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_z_xyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 23); 

                auto tex_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 24); 

                auto tey_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_z_xyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 24); 

                auto tex_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 25); 

                auto tey_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 25); 

                auto tez_z_xzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 25); 

                auto tex_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 26); 

                auto tey_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 26); 

                auto tez_z_yyy_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 26); 

                auto tex_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 27); 

                auto tey_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 27); 

                auto tez_z_yyz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 27); 

                auto tex_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 28); 

                auto tey_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 28); 

                auto tez_z_yzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 28); 

                auto tex_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * idx + 29); 

                auto tey_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 30 * bdim + 30 * idx + 29); 

                auto tez_z_zzz_1 = primBuffer.data(pidx_e_1_3_m1 + 60 * bdim + 30 * idx + 29); 

                auto ta_z_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 30); 

                auto ta_z_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 31); 

                auto ta_z_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 32); 

                auto ta_z_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 33); 

                auto ta_z_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 34); 

                auto ta_z_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 35); 

                auto ta_z_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 36); 

                auto ta_z_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 37); 

                auto ta_z_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 38); 

                auto ta_z_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 39); 

                auto ta_z_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 40); 

                auto ta_z_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 41); 

                auto ta_z_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 42); 

                auto ta_z_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 43); 

                auto ta_z_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 44); 

                // set up pointers to integrals

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

                // Batch of Integrals (225,270)

                #pragma omp simd aligned(fx, pa_z, pc_z, ta_z_xxxx_1, ta_z_xxxy_1, ta_z_xxxz_1, ta_z_xxyy_1, \
                                         ta_z_xxyz_1, ta_z_xxzz_1, ta_z_xyyy_1, ta_z_xyyz_1, ta_z_xyzz_1, ta_z_xzzz_1, \
                                         ta_z_yyyy_1, ta_z_yyyz_1, ta_z_yyzz_1, ta_z_yzzz_1, ta_z_zzzz_1, tex_0_xxxx_0, \
                                         tex_0_xxxx_1, tex_0_xxxy_0, tex_0_xxxy_1, tex_0_xxxz_0, tex_0_xxxz_1, tex_0_xxyy_0, \
                                         tex_0_xxyy_1, tex_0_xxyz_0, tex_0_xxyz_1, tex_0_xxzz_0, tex_0_xxzz_1, tex_0_xyyy_0, \
                                         tex_0_xyyy_1, tex_0_xyyz_0, tex_0_xyyz_1, tex_0_xyzz_0, tex_0_xyzz_1, tex_0_xzzz_0, \
                                         tex_0_xzzz_1, tex_0_yyyy_0, tex_0_yyyy_1, tex_0_yyyz_0, tex_0_yyyz_1, tex_0_yyzz_0, \
                                         tex_0_yyzz_1, tex_0_yzzz_0, tex_0_yzzz_1, tex_0_zzzz_0, tex_0_zzzz_1, tex_z_xxx_0, \
                                         tex_z_xxx_1, tex_z_xxxx_0, tex_z_xxxx_1, tex_z_xxxy_0, tex_z_xxxy_1, tex_z_xxxz_0, \
                                         tex_z_xxxz_1, tex_z_xxy_0, tex_z_xxy_1, tex_z_xxyy_0, tex_z_xxyy_1, tex_z_xxyz_0, \
                                         tex_z_xxyz_1, tex_z_xxz_0, tex_z_xxz_1, tex_z_xxzz_0, tex_z_xxzz_1, tex_z_xyy_0, \
                                         tex_z_xyy_1, tex_z_xyyy_0, tex_z_xyyy_1, tex_z_xyyz_0, tex_z_xyyz_1, tex_z_xyz_0, \
                                         tex_z_xyz_1, tex_z_xyzz_0, tex_z_xyzz_1, tex_z_xzz_0, tex_z_xzz_1, tex_z_xzzz_0, \
                                         tex_z_xzzz_1, tex_z_yyy_0, tex_z_yyy_1, tex_z_yyyy_0, tex_z_yyyy_1, tex_z_yyyz_0, \
                                         tex_z_yyyz_1, tex_z_yyz_0, tex_z_yyz_1, tex_z_yyzz_0, tex_z_yyzz_1, tex_z_yzz_0, \
                                         tex_z_yzz_1, tex_z_yzzz_0, tex_z_yzzz_1, tex_z_zzz_0, tex_z_zzz_1, tex_z_zzzz_0, \
                                         tex_z_zzzz_1, tex_zz_xxxx_0, tex_zz_xxxy_0, tex_zz_xxxz_0, tex_zz_xxyy_0, \
                                         tex_zz_xxyz_0, tex_zz_xxzz_0, tex_zz_xyyy_0, tex_zz_xyyz_0, tex_zz_xyzz_0, \
                                         tex_zz_xzzz_0, tex_zz_yyyy_0, tex_zz_yyyz_0, tex_zz_yyzz_0, tex_zz_yzzz_0, \
                                         tex_zz_zzzz_0, tey_0_xxxx_0, tey_0_xxxx_1, tey_0_xxxy_0, tey_0_xxxy_1, tey_0_xxxz_0, \
                                         tey_0_xxxz_1, tey_0_xxyy_0, tey_0_xxyy_1, tey_0_xxyz_0, tey_0_xxyz_1, tey_0_xxzz_0, \
                                         tey_0_xxzz_1, tey_0_xyyy_0, tey_0_xyyy_1, tey_0_xyyz_0, tey_0_xyyz_1, tey_0_xyzz_0, \
                                         tey_0_xyzz_1, tey_0_xzzz_0, tey_0_xzzz_1, tey_0_yyyy_0, tey_0_yyyy_1, tey_0_yyyz_0, \
                                         tey_0_yyyz_1, tey_0_yyzz_0, tey_0_yyzz_1, tey_0_yzzz_0, tey_0_yzzz_1, tey_0_zzzz_0, \
                                         tey_0_zzzz_1, tey_z_xxx_0, tey_z_xxx_1, tey_z_xxxx_0, tey_z_xxxx_1, tey_z_xxxy_0, \
                                         tey_z_xxxy_1, tey_z_xxxz_0, tey_z_xxxz_1, tey_z_xxy_0, tey_z_xxy_1, tey_z_xxyy_0, \
                                         tey_z_xxyy_1, tey_z_xxyz_0, tey_z_xxyz_1, tey_z_xxz_0, tey_z_xxz_1, tey_z_xxzz_0, \
                                         tey_z_xxzz_1, tey_z_xyy_0, tey_z_xyy_1, tey_z_xyyy_0, tey_z_xyyy_1, tey_z_xyyz_0, \
                                         tey_z_xyyz_1, tey_z_xyz_0, tey_z_xyz_1, tey_z_xyzz_0, tey_z_xyzz_1, tey_z_xzz_0, \
                                         tey_z_xzz_1, tey_z_xzzz_0, tey_z_xzzz_1, tey_z_yyy_0, tey_z_yyy_1, tey_z_yyyy_0, \
                                         tey_z_yyyy_1, tey_z_yyyz_0, tey_z_yyyz_1, tey_z_yyz_0, tey_z_yyz_1, tey_z_yyzz_0, \
                                         tey_z_yyzz_1, tey_z_yzz_0, tey_z_yzz_1, tey_z_yzzz_0, tey_z_yzzz_1, tey_z_zzz_0, \
                                         tey_z_zzz_1, tey_z_zzzz_0, tey_z_zzzz_1, tey_zz_xxxx_0, tey_zz_xxxy_0, \
                                         tey_zz_xxxz_0, tey_zz_xxyy_0, tey_zz_xxyz_0, tey_zz_xxzz_0, tey_zz_xyyy_0, \
                                         tey_zz_xyyz_0, tey_zz_xyzz_0, tey_zz_xzzz_0, tey_zz_yyyy_0, tey_zz_yyyz_0, \
                                         tey_zz_yyzz_0, tey_zz_yzzz_0, tey_zz_zzzz_0, tez_0_xxxx_0, tez_0_xxxx_1, \
                                         tez_0_xxxy_0, tez_0_xxxy_1, tez_0_xxxz_0, tez_0_xxxz_1, tez_0_xxyy_0, tez_0_xxyy_1, \
                                         tez_0_xxyz_0, tez_0_xxyz_1, tez_0_xxzz_0, tez_0_xxzz_1, tez_0_xyyy_0, tez_0_xyyy_1, \
                                         tez_0_xyyz_0, tez_0_xyyz_1, tez_0_xyzz_0, tez_0_xyzz_1, tez_0_xzzz_0, tez_0_xzzz_1, \
                                         tez_0_yyyy_0, tez_0_yyyy_1, tez_0_yyyz_0, tez_0_yyyz_1, tez_0_yyzz_0, tez_0_yyzz_1, \
                                         tez_0_yzzz_0, tez_0_yzzz_1, tez_0_zzzz_0, tez_0_zzzz_1, tez_z_xxx_0, tez_z_xxx_1, \
                                         tez_z_xxxx_0, tez_z_xxxx_1, tez_z_xxxy_0, tez_z_xxxy_1, tez_z_xxxz_0, tez_z_xxxz_1, \
                                         tez_z_xxy_0, tez_z_xxy_1, tez_z_xxyy_0, tez_z_xxyy_1, tez_z_xxyz_0, tez_z_xxyz_1, \
                                         tez_z_xxz_0, tez_z_xxz_1, tez_z_xxzz_0, tez_z_xxzz_1, tez_z_xyy_0, tez_z_xyy_1, \
                                         tez_z_xyyy_0, tez_z_xyyy_1, tez_z_xyyz_0, tez_z_xyyz_1, tez_z_xyz_0, tez_z_xyz_1, \
                                         tez_z_xyzz_0, tez_z_xyzz_1, tez_z_xzz_0, tez_z_xzz_1, tez_z_xzzz_0, tez_z_xzzz_1, \
                                         tez_z_yyy_0, tez_z_yyy_1, tez_z_yyyy_0, tez_z_yyyy_1, tez_z_yyyz_0, tez_z_yyyz_1, \
                                         tez_z_yyz_0, tez_z_yyz_1, tez_z_yyzz_0, tez_z_yyzz_1, tez_z_yzz_0, tez_z_yzz_1, \
                                         tez_z_yzzz_0, tez_z_yzzz_1, tez_z_zzz_0, tez_z_zzz_1, tez_z_zzzz_0, tez_z_zzzz_1, \
                                         tez_zz_xxxx_0, tez_zz_xxxy_0, tez_zz_xxxz_0, tez_zz_xxyy_0, tez_zz_xxyz_0, \
                                         tez_zz_xxzz_0, tez_zz_xyyy_0, tez_zz_xyyz_0, tez_zz_xyzz_0, tez_zz_xzzz_0, \
                                         tez_zz_yyyy_0, tez_zz_yyyz_0, tez_zz_yyzz_0, tez_zz_yzzz_0, tez_zz_zzzz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_zz_xxxx_0[j] = pa_z[j] * tex_z_xxxx_0[j] - pc_z[j] * tex_z_xxxx_1[j] + 0.5 * fl1_fx * tex_0_xxxx_0[j] - 0.5 * fl1_fx * tex_0_xxxx_1[j];

                    tey_zz_xxxx_0[j] = pa_z[j] * tey_z_xxxx_0[j] - pc_z[j] * tey_z_xxxx_1[j] + 0.5 * fl1_fx * tey_0_xxxx_0[j] - 0.5 * fl1_fx * tey_0_xxxx_1[j];

                    tez_zz_xxxx_0[j] = pa_z[j] * tez_z_xxxx_0[j] - pc_z[j] * tez_z_xxxx_1[j] + 0.5 * fl1_fx * tez_0_xxxx_0[j] - 0.5 * fl1_fx * tez_0_xxxx_1[j] + ta_z_xxxx_1[j];

                    tex_zz_xxxy_0[j] = pa_z[j] * tex_z_xxxy_0[j] - pc_z[j] * tex_z_xxxy_1[j] + 0.5 * fl1_fx * tex_0_xxxy_0[j] - 0.5 * fl1_fx * tex_0_xxxy_1[j];

                    tey_zz_xxxy_0[j] = pa_z[j] * tey_z_xxxy_0[j] - pc_z[j] * tey_z_xxxy_1[j] + 0.5 * fl1_fx * tey_0_xxxy_0[j] - 0.5 * fl1_fx * tey_0_xxxy_1[j];

                    tez_zz_xxxy_0[j] = pa_z[j] * tez_z_xxxy_0[j] - pc_z[j] * tez_z_xxxy_1[j] + 0.5 * fl1_fx * tez_0_xxxy_0[j] - 0.5 * fl1_fx * tez_0_xxxy_1[j] + ta_z_xxxy_1[j];

                    tex_zz_xxxz_0[j] = pa_z[j] * tex_z_xxxz_0[j] - pc_z[j] * tex_z_xxxz_1[j] + 0.5 * fl1_fx * tex_0_xxxz_0[j] - 0.5 * fl1_fx * tex_0_xxxz_1[j] + 0.5 * fl1_fx * tex_z_xxx_0[j] - 0.5 * fl1_fx * tex_z_xxx_1[j];

                    tey_zz_xxxz_0[j] = pa_z[j] * tey_z_xxxz_0[j] - pc_z[j] * tey_z_xxxz_1[j] + 0.5 * fl1_fx * tey_0_xxxz_0[j] - 0.5 * fl1_fx * tey_0_xxxz_1[j] + 0.5 * fl1_fx * tey_z_xxx_0[j] - 0.5 * fl1_fx * tey_z_xxx_1[j];

                    tez_zz_xxxz_0[j] = pa_z[j] * tez_z_xxxz_0[j] - pc_z[j] * tez_z_xxxz_1[j] + 0.5 * fl1_fx * tez_0_xxxz_0[j] - 0.5 * fl1_fx * tez_0_xxxz_1[j] + 0.5 * fl1_fx * tez_z_xxx_0[j] - 0.5 * fl1_fx * tez_z_xxx_1[j] + ta_z_xxxz_1[j];

                    tex_zz_xxyy_0[j] = pa_z[j] * tex_z_xxyy_0[j] - pc_z[j] * tex_z_xxyy_1[j] + 0.5 * fl1_fx * tex_0_xxyy_0[j] - 0.5 * fl1_fx * tex_0_xxyy_1[j];

                    tey_zz_xxyy_0[j] = pa_z[j] * tey_z_xxyy_0[j] - pc_z[j] * tey_z_xxyy_1[j] + 0.5 * fl1_fx * tey_0_xxyy_0[j] - 0.5 * fl1_fx * tey_0_xxyy_1[j];

                    tez_zz_xxyy_0[j] = pa_z[j] * tez_z_xxyy_0[j] - pc_z[j] * tez_z_xxyy_1[j] + 0.5 * fl1_fx * tez_0_xxyy_0[j] - 0.5 * fl1_fx * tez_0_xxyy_1[j] + ta_z_xxyy_1[j];

                    tex_zz_xxyz_0[j] = pa_z[j] * tex_z_xxyz_0[j] - pc_z[j] * tex_z_xxyz_1[j] + 0.5 * fl1_fx * tex_0_xxyz_0[j] - 0.5 * fl1_fx * tex_0_xxyz_1[j] + 0.5 * fl1_fx * tex_z_xxy_0[j] - 0.5 * fl1_fx * tex_z_xxy_1[j];

                    tey_zz_xxyz_0[j] = pa_z[j] * tey_z_xxyz_0[j] - pc_z[j] * tey_z_xxyz_1[j] + 0.5 * fl1_fx * tey_0_xxyz_0[j] - 0.5 * fl1_fx * tey_0_xxyz_1[j] + 0.5 * fl1_fx * tey_z_xxy_0[j] - 0.5 * fl1_fx * tey_z_xxy_1[j];

                    tez_zz_xxyz_0[j] = pa_z[j] * tez_z_xxyz_0[j] - pc_z[j] * tez_z_xxyz_1[j] + 0.5 * fl1_fx * tez_0_xxyz_0[j] - 0.5 * fl1_fx * tez_0_xxyz_1[j] + 0.5 * fl1_fx * tez_z_xxy_0[j] - 0.5 * fl1_fx * tez_z_xxy_1[j] + ta_z_xxyz_1[j];

                    tex_zz_xxzz_0[j] = pa_z[j] * tex_z_xxzz_0[j] - pc_z[j] * tex_z_xxzz_1[j] + 0.5 * fl1_fx * tex_0_xxzz_0[j] - 0.5 * fl1_fx * tex_0_xxzz_1[j] + fl1_fx * tex_z_xxz_0[j] - fl1_fx * tex_z_xxz_1[j];

                    tey_zz_xxzz_0[j] = pa_z[j] * tey_z_xxzz_0[j] - pc_z[j] * tey_z_xxzz_1[j] + 0.5 * fl1_fx * tey_0_xxzz_0[j] - 0.5 * fl1_fx * tey_0_xxzz_1[j] + fl1_fx * tey_z_xxz_0[j] - fl1_fx * tey_z_xxz_1[j];

                    tez_zz_xxzz_0[j] = pa_z[j] * tez_z_xxzz_0[j] - pc_z[j] * tez_z_xxzz_1[j] + 0.5 * fl1_fx * tez_0_xxzz_0[j] - 0.5 * fl1_fx * tez_0_xxzz_1[j] + fl1_fx * tez_z_xxz_0[j] - fl1_fx * tez_z_xxz_1[j] + ta_z_xxzz_1[j];

                    tex_zz_xyyy_0[j] = pa_z[j] * tex_z_xyyy_0[j] - pc_z[j] * tex_z_xyyy_1[j] + 0.5 * fl1_fx * tex_0_xyyy_0[j] - 0.5 * fl1_fx * tex_0_xyyy_1[j];

                    tey_zz_xyyy_0[j] = pa_z[j] * tey_z_xyyy_0[j] - pc_z[j] * tey_z_xyyy_1[j] + 0.5 * fl1_fx * tey_0_xyyy_0[j] - 0.5 * fl1_fx * tey_0_xyyy_1[j];

                    tez_zz_xyyy_0[j] = pa_z[j] * tez_z_xyyy_0[j] - pc_z[j] * tez_z_xyyy_1[j] + 0.5 * fl1_fx * tez_0_xyyy_0[j] - 0.5 * fl1_fx * tez_0_xyyy_1[j] + ta_z_xyyy_1[j];

                    tex_zz_xyyz_0[j] = pa_z[j] * tex_z_xyyz_0[j] - pc_z[j] * tex_z_xyyz_1[j] + 0.5 * fl1_fx * tex_0_xyyz_0[j] - 0.5 * fl1_fx * tex_0_xyyz_1[j] + 0.5 * fl1_fx * tex_z_xyy_0[j] - 0.5 * fl1_fx * tex_z_xyy_1[j];

                    tey_zz_xyyz_0[j] = pa_z[j] * tey_z_xyyz_0[j] - pc_z[j] * tey_z_xyyz_1[j] + 0.5 * fl1_fx * tey_0_xyyz_0[j] - 0.5 * fl1_fx * tey_0_xyyz_1[j] + 0.5 * fl1_fx * tey_z_xyy_0[j] - 0.5 * fl1_fx * tey_z_xyy_1[j];

                    tez_zz_xyyz_0[j] = pa_z[j] * tez_z_xyyz_0[j] - pc_z[j] * tez_z_xyyz_1[j] + 0.5 * fl1_fx * tez_0_xyyz_0[j] - 0.5 * fl1_fx * tez_0_xyyz_1[j] + 0.5 * fl1_fx * tez_z_xyy_0[j] - 0.5 * fl1_fx * tez_z_xyy_1[j] + ta_z_xyyz_1[j];

                    tex_zz_xyzz_0[j] = pa_z[j] * tex_z_xyzz_0[j] - pc_z[j] * tex_z_xyzz_1[j] + 0.5 * fl1_fx * tex_0_xyzz_0[j] - 0.5 * fl1_fx * tex_0_xyzz_1[j] + fl1_fx * tex_z_xyz_0[j] - fl1_fx * tex_z_xyz_1[j];

                    tey_zz_xyzz_0[j] = pa_z[j] * tey_z_xyzz_0[j] - pc_z[j] * tey_z_xyzz_1[j] + 0.5 * fl1_fx * tey_0_xyzz_0[j] - 0.5 * fl1_fx * tey_0_xyzz_1[j] + fl1_fx * tey_z_xyz_0[j] - fl1_fx * tey_z_xyz_1[j];

                    tez_zz_xyzz_0[j] = pa_z[j] * tez_z_xyzz_0[j] - pc_z[j] * tez_z_xyzz_1[j] + 0.5 * fl1_fx * tez_0_xyzz_0[j] - 0.5 * fl1_fx * tez_0_xyzz_1[j] + fl1_fx * tez_z_xyz_0[j] - fl1_fx * tez_z_xyz_1[j] + ta_z_xyzz_1[j];

                    tex_zz_xzzz_0[j] = pa_z[j] * tex_z_xzzz_0[j] - pc_z[j] * tex_z_xzzz_1[j] + 0.5 * fl1_fx * tex_0_xzzz_0[j] - 0.5 * fl1_fx * tex_0_xzzz_1[j] + 1.5 * fl1_fx * tex_z_xzz_0[j] - 1.5 * fl1_fx * tex_z_xzz_1[j];

                    tey_zz_xzzz_0[j] = pa_z[j] * tey_z_xzzz_0[j] - pc_z[j] * tey_z_xzzz_1[j] + 0.5 * fl1_fx * tey_0_xzzz_0[j] - 0.5 * fl1_fx * tey_0_xzzz_1[j] + 1.5 * fl1_fx * tey_z_xzz_0[j] - 1.5 * fl1_fx * tey_z_xzz_1[j];

                    tez_zz_xzzz_0[j] = pa_z[j] * tez_z_xzzz_0[j] - pc_z[j] * tez_z_xzzz_1[j] + 0.5 * fl1_fx * tez_0_xzzz_0[j] - 0.5 * fl1_fx * tez_0_xzzz_1[j] + 1.5 * fl1_fx * tez_z_xzz_0[j] - 1.5 * fl1_fx * tez_z_xzz_1[j] + ta_z_xzzz_1[j];

                    tex_zz_yyyy_0[j] = pa_z[j] * tex_z_yyyy_0[j] - pc_z[j] * tex_z_yyyy_1[j] + 0.5 * fl1_fx * tex_0_yyyy_0[j] - 0.5 * fl1_fx * tex_0_yyyy_1[j];

                    tey_zz_yyyy_0[j] = pa_z[j] * tey_z_yyyy_0[j] - pc_z[j] * tey_z_yyyy_1[j] + 0.5 * fl1_fx * tey_0_yyyy_0[j] - 0.5 * fl1_fx * tey_0_yyyy_1[j];

                    tez_zz_yyyy_0[j] = pa_z[j] * tez_z_yyyy_0[j] - pc_z[j] * tez_z_yyyy_1[j] + 0.5 * fl1_fx * tez_0_yyyy_0[j] - 0.5 * fl1_fx * tez_0_yyyy_1[j] + ta_z_yyyy_1[j];

                    tex_zz_yyyz_0[j] = pa_z[j] * tex_z_yyyz_0[j] - pc_z[j] * tex_z_yyyz_1[j] + 0.5 * fl1_fx * tex_0_yyyz_0[j] - 0.5 * fl1_fx * tex_0_yyyz_1[j] + 0.5 * fl1_fx * tex_z_yyy_0[j] - 0.5 * fl1_fx * tex_z_yyy_1[j];

                    tey_zz_yyyz_0[j] = pa_z[j] * tey_z_yyyz_0[j] - pc_z[j] * tey_z_yyyz_1[j] + 0.5 * fl1_fx * tey_0_yyyz_0[j] - 0.5 * fl1_fx * tey_0_yyyz_1[j] + 0.5 * fl1_fx * tey_z_yyy_0[j] - 0.5 * fl1_fx * tey_z_yyy_1[j];

                    tez_zz_yyyz_0[j] = pa_z[j] * tez_z_yyyz_0[j] - pc_z[j] * tez_z_yyyz_1[j] + 0.5 * fl1_fx * tez_0_yyyz_0[j] - 0.5 * fl1_fx * tez_0_yyyz_1[j] + 0.5 * fl1_fx * tez_z_yyy_0[j] - 0.5 * fl1_fx * tez_z_yyy_1[j] + ta_z_yyyz_1[j];

                    tex_zz_yyzz_0[j] = pa_z[j] * tex_z_yyzz_0[j] - pc_z[j] * tex_z_yyzz_1[j] + 0.5 * fl1_fx * tex_0_yyzz_0[j] - 0.5 * fl1_fx * tex_0_yyzz_1[j] + fl1_fx * tex_z_yyz_0[j] - fl1_fx * tex_z_yyz_1[j];

                    tey_zz_yyzz_0[j] = pa_z[j] * tey_z_yyzz_0[j] - pc_z[j] * tey_z_yyzz_1[j] + 0.5 * fl1_fx * tey_0_yyzz_0[j] - 0.5 * fl1_fx * tey_0_yyzz_1[j] + fl1_fx * tey_z_yyz_0[j] - fl1_fx * tey_z_yyz_1[j];

                    tez_zz_yyzz_0[j] = pa_z[j] * tez_z_yyzz_0[j] - pc_z[j] * tez_z_yyzz_1[j] + 0.5 * fl1_fx * tez_0_yyzz_0[j] - 0.5 * fl1_fx * tez_0_yyzz_1[j] + fl1_fx * tez_z_yyz_0[j] - fl1_fx * tez_z_yyz_1[j] + ta_z_yyzz_1[j];

                    tex_zz_yzzz_0[j] = pa_z[j] * tex_z_yzzz_0[j] - pc_z[j] * tex_z_yzzz_1[j] + 0.5 * fl1_fx * tex_0_yzzz_0[j] - 0.5 * fl1_fx * tex_0_yzzz_1[j] + 1.5 * fl1_fx * tex_z_yzz_0[j] - 1.5 * fl1_fx * tex_z_yzz_1[j];

                    tey_zz_yzzz_0[j] = pa_z[j] * tey_z_yzzz_0[j] - pc_z[j] * tey_z_yzzz_1[j] + 0.5 * fl1_fx * tey_0_yzzz_0[j] - 0.5 * fl1_fx * tey_0_yzzz_1[j] + 1.5 * fl1_fx * tey_z_yzz_0[j] - 1.5 * fl1_fx * tey_z_yzz_1[j];

                    tez_zz_yzzz_0[j] = pa_z[j] * tez_z_yzzz_0[j] - pc_z[j] * tez_z_yzzz_1[j] + 0.5 * fl1_fx * tez_0_yzzz_0[j] - 0.5 * fl1_fx * tez_0_yzzz_1[j] + 1.5 * fl1_fx * tez_z_yzz_0[j] - 1.5 * fl1_fx * tez_z_yzz_1[j] + ta_z_yzzz_1[j];

                    tex_zz_zzzz_0[j] = pa_z[j] * tex_z_zzzz_0[j] - pc_z[j] * tex_z_zzzz_1[j] + 0.5 * fl1_fx * tex_0_zzzz_0[j] - 0.5 * fl1_fx * tex_0_zzzz_1[j] + 2.0 * fl1_fx * tex_z_zzz_0[j] - 2.0 * fl1_fx * tex_z_zzz_1[j];

                    tey_zz_zzzz_0[j] = pa_z[j] * tey_z_zzzz_0[j] - pc_z[j] * tey_z_zzzz_1[j] + 0.5 * fl1_fx * tey_0_zzzz_0[j] - 0.5 * fl1_fx * tey_0_zzzz_1[j] + 2.0 * fl1_fx * tey_z_zzz_0[j] - 2.0 * fl1_fx * tey_z_zzz_1[j];

                    tez_zz_zzzz_0[j] = pa_z[j] * tez_z_zzzz_0[j] - pc_z[j] * tez_z_zzzz_1[j] + 0.5 * fl1_fx * tez_0_zzzz_0[j] - 0.5 * fl1_fx * tez_0_zzzz_1[j] + 2.0 * fl1_fx * tez_z_zzz_0[j] - 2.0 * fl1_fx * tez_z_zzz_1[j] + ta_z_zzzz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForGD(      CMemBlock2D<double>& primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CMemBlock2D<double>& pcDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto)
    {
        efieldrecfunc::compElectricFieldForGD_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   pcDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        efieldrecfunc::compElectricFieldForGD_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    pcDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        efieldrecfunc::compElectricFieldForGD_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     pcDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        efieldrecfunc::compElectricFieldForGD_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForGD_180_225(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        efieldrecfunc::compElectricFieldForGD_225_270(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      pcDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compElectricFieldForGD_0_45(      CMemBlock2D<double>& primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
    {
        // Batch of Integrals (0,45)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {4, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_4_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                auto tex_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx); 

                auto tey_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx); 

                auto tez_xxx_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx); 

                auto tex_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 1); 

                auto tey_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 1); 

                auto tez_xxx_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 1); 

                auto tex_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 2); 

                auto tey_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 2); 

                auto tez_xxx_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 2); 

                auto tex_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 3); 

                auto tey_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 3); 

                auto tez_xxx_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 3); 

                auto tex_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 4); 

                auto tey_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 4); 

                auto tez_xxx_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 4); 

                auto tex_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 5); 

                auto tey_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 5); 

                auto tez_xxx_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 5); 

                auto tex_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 6); 

                auto tey_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 6); 

                auto tez_xxy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 6); 

                auto tex_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 7); 

                auto tey_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 7); 

                auto tez_xxy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 7); 

                auto tex_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 8); 

                auto tey_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 8); 

                auto tez_xxy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 8); 

                auto tex_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 9); 

                auto tey_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 9); 

                auto tez_xxy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 9); 

                auto tex_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 10); 

                auto tey_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 10); 

                auto tez_xxy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 10); 

                auto tex_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 11); 

                auto tey_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 11); 

                auto tez_xxy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 11); 

                auto tex_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 12); 

                auto tey_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 12); 

                auto tez_xxz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 12); 

                auto tex_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 13); 

                auto tey_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 13); 

                auto tez_xxz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 13); 

                auto tex_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 14); 

                auto tey_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 14); 

                auto tez_xxz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 14); 

                auto tex_xxx_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx); 

                auto tey_xxx_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx); 

                auto tez_xxx_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx); 

                auto tex_xxx_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 1); 

                auto tey_xxx_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 1); 

                auto tez_xxx_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 1); 

                auto tex_xxx_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 2); 

                auto tey_xxx_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 2); 

                auto tez_xxx_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 2); 

                auto tex_xxx_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 3); 

                auto tey_xxx_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 3); 

                auto tez_xxx_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 3); 

                auto tex_xxx_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 4); 

                auto tey_xxx_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 4); 

                auto tez_xxx_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 4); 

                auto tex_xxx_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 5); 

                auto tey_xxx_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 5); 

                auto tez_xxx_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 5); 

                auto tex_xxy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 6); 

                auto tey_xxy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 6); 

                auto tez_xxy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 6); 

                auto tex_xxy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 7); 

                auto tey_xxy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 7); 

                auto tez_xxy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 7); 

                auto tex_xxy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 8); 

                auto tey_xxy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 8); 

                auto tez_xxy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 8); 

                auto tex_xxy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 9); 

                auto tey_xxy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 9); 

                auto tez_xxy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 9); 

                auto tex_xxy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 10); 

                auto tey_xxy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 10); 

                auto tez_xxy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 10); 

                auto tex_xxy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 11); 

                auto tey_xxy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 11); 

                auto tez_xxy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 11); 

                auto tex_xxz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 12); 

                auto tey_xxz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 12); 

                auto tez_xxz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 12); 

                auto tex_xxz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 13); 

                auto tey_xxz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 13); 

                auto tez_xxz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 13); 

                auto tex_xxz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 14); 

                auto tey_xxz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 14); 

                auto tez_xxz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 14); 

                auto tex_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx); 

                auto tey_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx); 

                auto tez_xx_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx); 

                auto tex_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 1); 

                auto tey_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 1); 

                auto tez_xx_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 1); 

                auto tex_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 2); 

                auto tey_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 2); 

                auto tez_xx_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 2); 

                auto tex_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 3); 

                auto tey_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 3); 

                auto tez_xx_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 3); 

                auto tex_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 4); 

                auto tey_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 4); 

                auto tez_xx_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 4); 

                auto tex_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 5); 

                auto tey_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 5); 

                auto tez_xx_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 5); 

                auto tex_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 6); 

                auto tey_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 6); 

                auto tez_xy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 6); 

                auto tex_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 7); 

                auto tey_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 7); 

                auto tez_xy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 7); 

                auto tex_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 8); 

                auto tey_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 8); 

                auto tez_xy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 8); 

                auto tex_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 9); 

                auto tey_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 9); 

                auto tez_xy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 9); 

                auto tex_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 10); 

                auto tey_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 10); 

                auto tez_xy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 10); 

                auto tex_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 11); 

                auto tey_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 11); 

                auto tez_xy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 11); 

                auto tex_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 12); 

                auto tey_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 12); 

                auto tez_xz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 12); 

                auto tex_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 13); 

                auto tey_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 13); 

                auto tez_xz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 13); 

                auto tex_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 14); 

                auto tey_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 14); 

                auto tez_xz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 14); 

                auto tex_xx_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx); 

                auto tey_xx_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx); 

                auto tez_xx_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx); 

                auto tex_xx_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 1); 

                auto tey_xx_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 1); 

                auto tez_xx_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 1); 

                auto tex_xx_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 2); 

                auto tey_xx_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 2); 

                auto tez_xx_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 2); 

                auto tex_xx_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 3); 

                auto tey_xx_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 3); 

                auto tez_xx_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 3); 

                auto tex_xx_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 4); 

                auto tey_xx_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 4); 

                auto tez_xx_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 4); 

                auto tex_xx_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 5); 

                auto tey_xx_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 5); 

                auto tez_xx_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 5); 

                auto tex_xy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 6); 

                auto tey_xy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 6); 

                auto tez_xy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 6); 

                auto tex_xy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 7); 

                auto tey_xy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 7); 

                auto tez_xy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 7); 

                auto tex_xy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 8); 

                auto tey_xy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 8); 

                auto tez_xy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 8); 

                auto tex_xy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 9); 

                auto tey_xy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 9); 

                auto tez_xy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 9); 

                auto tex_xy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 10); 

                auto tey_xy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 10); 

                auto tez_xy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 10); 

                auto tex_xy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 11); 

                auto tey_xy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 11); 

                auto tez_xy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 11); 

                auto tex_xz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 12); 

                auto tey_xz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 12); 

                auto tez_xz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 12); 

                auto tex_xz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 13); 

                auto tey_xz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 13); 

                auto tez_xz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 13); 

                auto tex_xz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 14); 

                auto tey_xz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 14); 

                auto tez_xz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 14); 

                auto tex_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx); 

                auto tey_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx); 

                auto tez_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx); 

                auto tex_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 1); 

                auto tey_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 1); 

                auto tez_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 1); 

                auto tex_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 2); 

                auto tey_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 2); 

                auto tez_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 2); 

                auto tex_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 3); 

                auto tey_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 3); 

                auto tez_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 3); 

                auto tex_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 4); 

                auto tey_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 4); 

                auto tez_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 4); 

                auto tex_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 5); 

                auto tey_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 5); 

                auto tez_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 5); 

                auto tex_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 6); 

                auto tey_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 6); 

                auto tez_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 6); 

                auto tex_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 7); 

                auto tey_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 7); 

                auto tez_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 7); 

                auto tex_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 8); 

                auto tey_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 8); 

                auto tez_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 8); 

                auto tex_xxx_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx); 

                auto tey_xxx_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx); 

                auto tez_xxx_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx); 

                auto tex_xxx_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 1); 

                auto tey_xxx_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 1); 

                auto tez_xxx_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 1); 

                auto tex_xxx_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 2); 

                auto tey_xxx_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 2); 

                auto tez_xxx_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 2); 

                auto tex_xxy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 3); 

                auto tey_xxy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 3); 

                auto tez_xxy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 3); 

                auto tex_xxy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 4); 

                auto tey_xxy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 4); 

                auto tez_xxy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 4); 

                auto tex_xxy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 5); 

                auto tey_xxy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 5); 

                auto tez_xxy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 5); 

                auto tex_xxz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 6); 

                auto tey_xxz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 6); 

                auto tez_xxz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 6); 

                auto tex_xxz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 7); 

                auto tey_xxz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 7); 

                auto tez_xxz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 7); 

                auto tex_xxz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 8); 

                auto tey_xxz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 8); 

                auto tez_xxz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 8); 

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

                // set up pointers to integrals

                auto tex_xxxx_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx); 

                auto tey_xxxx_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx); 

                auto tez_xxxx_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx); 

                auto tex_xxxx_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 1); 

                auto tey_xxxx_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 1); 

                auto tez_xxxx_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 1); 

                auto tex_xxxx_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 2); 

                auto tey_xxxx_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 2); 

                auto tez_xxxx_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 2); 

                auto tex_xxxx_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 3); 

                auto tey_xxxx_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 3); 

                auto tez_xxxx_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 3); 

                auto tex_xxxx_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 4); 

                auto tey_xxxx_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 4); 

                auto tez_xxxx_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 4); 

                auto tex_xxxx_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 5); 

                auto tey_xxxx_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 5); 

                auto tez_xxxx_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 5); 

                auto tex_xxxy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 6); 

                auto tey_xxxy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 6); 

                auto tez_xxxy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 6); 

                auto tex_xxxy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 7); 

                auto tey_xxxy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 7); 

                auto tez_xxxy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 7); 

                auto tex_xxxy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 8); 

                auto tey_xxxy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 8); 

                auto tez_xxxy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 8); 

                auto tex_xxxy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 9); 

                auto tey_xxxy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 9); 

                auto tez_xxxy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 9); 

                auto tex_xxxy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 10); 

                auto tey_xxxy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 10); 

                auto tez_xxxy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 10); 

                auto tex_xxxy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 11); 

                auto tey_xxxy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 11); 

                auto tez_xxxy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 11); 

                auto tex_xxxz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 12); 

                auto tey_xxxz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 12); 

                auto tez_xxxz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 12); 

                auto tex_xxxz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 13); 

                auto tey_xxxz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 13); 

                auto tez_xxxz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 13); 

                auto tex_xxxz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 14); 

                auto tey_xxxz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 14); 

                auto tez_xxxz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 14); 

                // Batch of Integrals (0,45)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxx_xx_1, ta_xxx_xy_1, ta_xxx_xz_1, ta_xxx_yy_1, \
                                         ta_xxx_yz_1, ta_xxx_zz_1, ta_xxy_xx_1, ta_xxy_xy_1, ta_xxy_xz_1, ta_xxy_yy_1, \
                                         ta_xxy_yz_1, ta_xxy_zz_1, ta_xxz_xx_1, ta_xxz_xy_1, ta_xxz_xz_1, tex_xx_xx_0, \
                                         tex_xx_xx_1, tex_xx_xy_0, tex_xx_xy_1, tex_xx_xz_0, tex_xx_xz_1, tex_xx_yy_0, \
                                         tex_xx_yy_1, tex_xx_yz_0, tex_xx_yz_1, tex_xx_zz_0, tex_xx_zz_1, tex_xxx_x_0, \
                                         tex_xxx_x_1, tex_xxx_xx_0, tex_xxx_xx_1, tex_xxx_xy_0, tex_xxx_xy_1, tex_xxx_xz_0, \
                                         tex_xxx_xz_1, tex_xxx_y_0, tex_xxx_y_1, tex_xxx_yy_0, tex_xxx_yy_1, tex_xxx_yz_0, \
                                         tex_xxx_yz_1, tex_xxx_z_0, tex_xxx_z_1, tex_xxx_zz_0, tex_xxx_zz_1, tex_xxxx_xx_0, \
                                         tex_xxxx_xy_0, tex_xxxx_xz_0, tex_xxxx_yy_0, tex_xxxx_yz_0, tex_xxxx_zz_0, \
                                         tex_xxxy_xx_0, tex_xxxy_xy_0, tex_xxxy_xz_0, tex_xxxy_yy_0, tex_xxxy_yz_0, \
                                         tex_xxxy_zz_0, tex_xxxz_xx_0, tex_xxxz_xy_0, tex_xxxz_xz_0, tex_xxy_x_0, \
                                         tex_xxy_x_1, tex_xxy_xx_0, tex_xxy_xx_1, tex_xxy_xy_0, tex_xxy_xy_1, tex_xxy_xz_0, \
                                         tex_xxy_xz_1, tex_xxy_y_0, tex_xxy_y_1, tex_xxy_yy_0, tex_xxy_yy_1, tex_xxy_yz_0, \
                                         tex_xxy_yz_1, tex_xxy_z_0, tex_xxy_z_1, tex_xxy_zz_0, tex_xxy_zz_1, tex_xxz_x_0, \
                                         tex_xxz_x_1, tex_xxz_xx_0, tex_xxz_xx_1, tex_xxz_xy_0, tex_xxz_xy_1, tex_xxz_xz_0, \
                                         tex_xxz_xz_1, tex_xxz_y_0, tex_xxz_y_1, tex_xxz_z_0, tex_xxz_z_1, tex_xy_xx_0, \
                                         tex_xy_xx_1, tex_xy_xy_0, tex_xy_xy_1, tex_xy_xz_0, tex_xy_xz_1, tex_xy_yy_0, \
                                         tex_xy_yy_1, tex_xy_yz_0, tex_xy_yz_1, tex_xy_zz_0, tex_xy_zz_1, tex_xz_xx_0, \
                                         tex_xz_xx_1, tex_xz_xy_0, tex_xz_xy_1, tex_xz_xz_0, tex_xz_xz_1, tey_xx_xx_0, \
                                         tey_xx_xx_1, tey_xx_xy_0, tey_xx_xy_1, tey_xx_xz_0, tey_xx_xz_1, tey_xx_yy_0, \
                                         tey_xx_yy_1, tey_xx_yz_0, tey_xx_yz_1, tey_xx_zz_0, tey_xx_zz_1, tey_xxx_x_0, \
                                         tey_xxx_x_1, tey_xxx_xx_0, tey_xxx_xx_1, tey_xxx_xy_0, tey_xxx_xy_1, tey_xxx_xz_0, \
                                         tey_xxx_xz_1, tey_xxx_y_0, tey_xxx_y_1, tey_xxx_yy_0, tey_xxx_yy_1, tey_xxx_yz_0, \
                                         tey_xxx_yz_1, tey_xxx_z_0, tey_xxx_z_1, tey_xxx_zz_0, tey_xxx_zz_1, tey_xxxx_xx_0, \
                                         tey_xxxx_xy_0, tey_xxxx_xz_0, tey_xxxx_yy_0, tey_xxxx_yz_0, tey_xxxx_zz_0, \
                                         tey_xxxy_xx_0, tey_xxxy_xy_0, tey_xxxy_xz_0, tey_xxxy_yy_0, tey_xxxy_yz_0, \
                                         tey_xxxy_zz_0, tey_xxxz_xx_0, tey_xxxz_xy_0, tey_xxxz_xz_0, tey_xxy_x_0, \
                                         tey_xxy_x_1, tey_xxy_xx_0, tey_xxy_xx_1, tey_xxy_xy_0, tey_xxy_xy_1, tey_xxy_xz_0, \
                                         tey_xxy_xz_1, tey_xxy_y_0, tey_xxy_y_1, tey_xxy_yy_0, tey_xxy_yy_1, tey_xxy_yz_0, \
                                         tey_xxy_yz_1, tey_xxy_z_0, tey_xxy_z_1, tey_xxy_zz_0, tey_xxy_zz_1, tey_xxz_x_0, \
                                         tey_xxz_x_1, tey_xxz_xx_0, tey_xxz_xx_1, tey_xxz_xy_0, tey_xxz_xy_1, tey_xxz_xz_0, \
                                         tey_xxz_xz_1, tey_xxz_y_0, tey_xxz_y_1, tey_xxz_z_0, tey_xxz_z_1, tey_xy_xx_0, \
                                         tey_xy_xx_1, tey_xy_xy_0, tey_xy_xy_1, tey_xy_xz_0, tey_xy_xz_1, tey_xy_yy_0, \
                                         tey_xy_yy_1, tey_xy_yz_0, tey_xy_yz_1, tey_xy_zz_0, tey_xy_zz_1, tey_xz_xx_0, \
                                         tey_xz_xx_1, tey_xz_xy_0, tey_xz_xy_1, tey_xz_xz_0, tey_xz_xz_1, tez_xx_xx_0, \
                                         tez_xx_xx_1, tez_xx_xy_0, tez_xx_xy_1, tez_xx_xz_0, tez_xx_xz_1, tez_xx_yy_0, \
                                         tez_xx_yy_1, tez_xx_yz_0, tez_xx_yz_1, tez_xx_zz_0, tez_xx_zz_1, tez_xxx_x_0, \
                                         tez_xxx_x_1, tez_xxx_xx_0, tez_xxx_xx_1, tez_xxx_xy_0, tez_xxx_xy_1, tez_xxx_xz_0, \
                                         tez_xxx_xz_1, tez_xxx_y_0, tez_xxx_y_1, tez_xxx_yy_0, tez_xxx_yy_1, tez_xxx_yz_0, \
                                         tez_xxx_yz_1, tez_xxx_z_0, tez_xxx_z_1, tez_xxx_zz_0, tez_xxx_zz_1, tez_xxxx_xx_0, \
                                         tez_xxxx_xy_0, tez_xxxx_xz_0, tez_xxxx_yy_0, tez_xxxx_yz_0, tez_xxxx_zz_0, \
                                         tez_xxxy_xx_0, tez_xxxy_xy_0, tez_xxxy_xz_0, tez_xxxy_yy_0, tez_xxxy_yz_0, \
                                         tez_xxxy_zz_0, tez_xxxz_xx_0, tez_xxxz_xy_0, tez_xxxz_xz_0, tez_xxy_x_0, \
                                         tez_xxy_x_1, tez_xxy_xx_0, tez_xxy_xx_1, tez_xxy_xy_0, tez_xxy_xy_1, tez_xxy_xz_0, \
                                         tez_xxy_xz_1, tez_xxy_y_0, tez_xxy_y_1, tez_xxy_yy_0, tez_xxy_yy_1, tez_xxy_yz_0, \
                                         tez_xxy_yz_1, tez_xxy_z_0, tez_xxy_z_1, tez_xxy_zz_0, tez_xxy_zz_1, tez_xxz_x_0, \
                                         tez_xxz_x_1, tez_xxz_xx_0, tez_xxz_xx_1, tez_xxz_xy_0, tez_xxz_xy_1, tez_xxz_xz_0, \
                                         tez_xxz_xz_1, tez_xxz_y_0, tez_xxz_y_1, tez_xxz_z_0, tez_xxz_z_1, tez_xy_xx_0, \
                                         tez_xy_xx_1, tez_xy_xy_0, tez_xy_xy_1, tez_xy_xz_0, tez_xy_xz_1, tez_xy_yy_0, \
                                         tez_xy_yy_1, tez_xy_yz_0, tez_xy_yz_1, tez_xy_zz_0, tez_xy_zz_1, tez_xz_xx_0, \
                                         tez_xz_xx_1, tez_xz_xy_0, tez_xz_xy_1, tez_xz_xz_0, tez_xz_xz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xxxx_xx_0[j] = pa_x[j] * tex_xxx_xx_0[j] - pc_x[j] * tex_xxx_xx_1[j] + 1.5 * fl1_fx * tex_xx_xx_0[j] - 1.5 * fl1_fx * tex_xx_xx_1[j] + fl1_fx * tex_xxx_x_0[j] - fl1_fx * tex_xxx_x_1[j] + ta_xxx_xx_1[j];

                    tey_xxxx_xx_0[j] = pa_x[j] * tey_xxx_xx_0[j] - pc_x[j] * tey_xxx_xx_1[j] + 1.5 * fl1_fx * tey_xx_xx_0[j] - 1.5 * fl1_fx * tey_xx_xx_1[j] + fl1_fx * tey_xxx_x_0[j] - fl1_fx * tey_xxx_x_1[j];

                    tez_xxxx_xx_0[j] = pa_x[j] * tez_xxx_xx_0[j] - pc_x[j] * tez_xxx_xx_1[j] + 1.5 * fl1_fx * tez_xx_xx_0[j] - 1.5 * fl1_fx * tez_xx_xx_1[j] + fl1_fx * tez_xxx_x_0[j] - fl1_fx * tez_xxx_x_1[j];

                    tex_xxxx_xy_0[j] = pa_x[j] * tex_xxx_xy_0[j] - pc_x[j] * tex_xxx_xy_1[j] + 1.5 * fl1_fx * tex_xx_xy_0[j] - 1.5 * fl1_fx * tex_xx_xy_1[j] + 0.5 * fl1_fx * tex_xxx_y_0[j] - 0.5 * fl1_fx * tex_xxx_y_1[j] + ta_xxx_xy_1[j];

                    tey_xxxx_xy_0[j] = pa_x[j] * tey_xxx_xy_0[j] - pc_x[j] * tey_xxx_xy_1[j] + 1.5 * fl1_fx * tey_xx_xy_0[j] - 1.5 * fl1_fx * tey_xx_xy_1[j] + 0.5 * fl1_fx * tey_xxx_y_0[j] - 0.5 * fl1_fx * tey_xxx_y_1[j];

                    tez_xxxx_xy_0[j] = pa_x[j] * tez_xxx_xy_0[j] - pc_x[j] * tez_xxx_xy_1[j] + 1.5 * fl1_fx * tez_xx_xy_0[j] - 1.5 * fl1_fx * tez_xx_xy_1[j] + 0.5 * fl1_fx * tez_xxx_y_0[j] - 0.5 * fl1_fx * tez_xxx_y_1[j];

                    tex_xxxx_xz_0[j] = pa_x[j] * tex_xxx_xz_0[j] - pc_x[j] * tex_xxx_xz_1[j] + 1.5 * fl1_fx * tex_xx_xz_0[j] - 1.5 * fl1_fx * tex_xx_xz_1[j] + 0.5 * fl1_fx * tex_xxx_z_0[j] - 0.5 * fl1_fx * tex_xxx_z_1[j] + ta_xxx_xz_1[j];

                    tey_xxxx_xz_0[j] = pa_x[j] * tey_xxx_xz_0[j] - pc_x[j] * tey_xxx_xz_1[j] + 1.5 * fl1_fx * tey_xx_xz_0[j] - 1.5 * fl1_fx * tey_xx_xz_1[j] + 0.5 * fl1_fx * tey_xxx_z_0[j] - 0.5 * fl1_fx * tey_xxx_z_1[j];

                    tez_xxxx_xz_0[j] = pa_x[j] * tez_xxx_xz_0[j] - pc_x[j] * tez_xxx_xz_1[j] + 1.5 * fl1_fx * tez_xx_xz_0[j] - 1.5 * fl1_fx * tez_xx_xz_1[j] + 0.5 * fl1_fx * tez_xxx_z_0[j] - 0.5 * fl1_fx * tez_xxx_z_1[j];

                    tex_xxxx_yy_0[j] = pa_x[j] * tex_xxx_yy_0[j] - pc_x[j] * tex_xxx_yy_1[j] + 1.5 * fl1_fx * tex_xx_yy_0[j] - 1.5 * fl1_fx * tex_xx_yy_1[j] + ta_xxx_yy_1[j];

                    tey_xxxx_yy_0[j] = pa_x[j] * tey_xxx_yy_0[j] - pc_x[j] * tey_xxx_yy_1[j] + 1.5 * fl1_fx * tey_xx_yy_0[j] - 1.5 * fl1_fx * tey_xx_yy_1[j];

                    tez_xxxx_yy_0[j] = pa_x[j] * tez_xxx_yy_0[j] - pc_x[j] * tez_xxx_yy_1[j] + 1.5 * fl1_fx * tez_xx_yy_0[j] - 1.5 * fl1_fx * tez_xx_yy_1[j];

                    tex_xxxx_yz_0[j] = pa_x[j] * tex_xxx_yz_0[j] - pc_x[j] * tex_xxx_yz_1[j] + 1.5 * fl1_fx * tex_xx_yz_0[j] - 1.5 * fl1_fx * tex_xx_yz_1[j] + ta_xxx_yz_1[j];

                    tey_xxxx_yz_0[j] = pa_x[j] * tey_xxx_yz_0[j] - pc_x[j] * tey_xxx_yz_1[j] + 1.5 * fl1_fx * tey_xx_yz_0[j] - 1.5 * fl1_fx * tey_xx_yz_1[j];

                    tez_xxxx_yz_0[j] = pa_x[j] * tez_xxx_yz_0[j] - pc_x[j] * tez_xxx_yz_1[j] + 1.5 * fl1_fx * tez_xx_yz_0[j] - 1.5 * fl1_fx * tez_xx_yz_1[j];

                    tex_xxxx_zz_0[j] = pa_x[j] * tex_xxx_zz_0[j] - pc_x[j] * tex_xxx_zz_1[j] + 1.5 * fl1_fx * tex_xx_zz_0[j] - 1.5 * fl1_fx * tex_xx_zz_1[j] + ta_xxx_zz_1[j];

                    tey_xxxx_zz_0[j] = pa_x[j] * tey_xxx_zz_0[j] - pc_x[j] * tey_xxx_zz_1[j] + 1.5 * fl1_fx * tey_xx_zz_0[j] - 1.5 * fl1_fx * tey_xx_zz_1[j];

                    tez_xxxx_zz_0[j] = pa_x[j] * tez_xxx_zz_0[j] - pc_x[j] * tez_xxx_zz_1[j] + 1.5 * fl1_fx * tez_xx_zz_0[j] - 1.5 * fl1_fx * tez_xx_zz_1[j];

                    tex_xxxy_xx_0[j] = pa_x[j] * tex_xxy_xx_0[j] - pc_x[j] * tex_xxy_xx_1[j] + fl1_fx * tex_xy_xx_0[j] - fl1_fx * tex_xy_xx_1[j] + fl1_fx * tex_xxy_x_0[j] - fl1_fx * tex_xxy_x_1[j] + ta_xxy_xx_1[j];

                    tey_xxxy_xx_0[j] = pa_x[j] * tey_xxy_xx_0[j] - pc_x[j] * tey_xxy_xx_1[j] + fl1_fx * tey_xy_xx_0[j] - fl1_fx * tey_xy_xx_1[j] + fl1_fx * tey_xxy_x_0[j] - fl1_fx * tey_xxy_x_1[j];

                    tez_xxxy_xx_0[j] = pa_x[j] * tez_xxy_xx_0[j] - pc_x[j] * tez_xxy_xx_1[j] + fl1_fx * tez_xy_xx_0[j] - fl1_fx * tez_xy_xx_1[j] + fl1_fx * tez_xxy_x_0[j] - fl1_fx * tez_xxy_x_1[j];

                    tex_xxxy_xy_0[j] = pa_x[j] * tex_xxy_xy_0[j] - pc_x[j] * tex_xxy_xy_1[j] + fl1_fx * tex_xy_xy_0[j] - fl1_fx * tex_xy_xy_1[j] + 0.5 * fl1_fx * tex_xxy_y_0[j] - 0.5 * fl1_fx * tex_xxy_y_1[j] + ta_xxy_xy_1[j];

                    tey_xxxy_xy_0[j] = pa_x[j] * tey_xxy_xy_0[j] - pc_x[j] * tey_xxy_xy_1[j] + fl1_fx * tey_xy_xy_0[j] - fl1_fx * tey_xy_xy_1[j] + 0.5 * fl1_fx * tey_xxy_y_0[j] - 0.5 * fl1_fx * tey_xxy_y_1[j];

                    tez_xxxy_xy_0[j] = pa_x[j] * tez_xxy_xy_0[j] - pc_x[j] * tez_xxy_xy_1[j] + fl1_fx * tez_xy_xy_0[j] - fl1_fx * tez_xy_xy_1[j] + 0.5 * fl1_fx * tez_xxy_y_0[j] - 0.5 * fl1_fx * tez_xxy_y_1[j];

                    tex_xxxy_xz_0[j] = pa_x[j] * tex_xxy_xz_0[j] - pc_x[j] * tex_xxy_xz_1[j] + fl1_fx * tex_xy_xz_0[j] - fl1_fx * tex_xy_xz_1[j] + 0.5 * fl1_fx * tex_xxy_z_0[j] - 0.5 * fl1_fx * tex_xxy_z_1[j] + ta_xxy_xz_1[j];

                    tey_xxxy_xz_0[j] = pa_x[j] * tey_xxy_xz_0[j] - pc_x[j] * tey_xxy_xz_1[j] + fl1_fx * tey_xy_xz_0[j] - fl1_fx * tey_xy_xz_1[j] + 0.5 * fl1_fx * tey_xxy_z_0[j] - 0.5 * fl1_fx * tey_xxy_z_1[j];

                    tez_xxxy_xz_0[j] = pa_x[j] * tez_xxy_xz_0[j] - pc_x[j] * tez_xxy_xz_1[j] + fl1_fx * tez_xy_xz_0[j] - fl1_fx * tez_xy_xz_1[j] + 0.5 * fl1_fx * tez_xxy_z_0[j] - 0.5 * fl1_fx * tez_xxy_z_1[j];

                    tex_xxxy_yy_0[j] = pa_x[j] * tex_xxy_yy_0[j] - pc_x[j] * tex_xxy_yy_1[j] + fl1_fx * tex_xy_yy_0[j] - fl1_fx * tex_xy_yy_1[j] + ta_xxy_yy_1[j];

                    tey_xxxy_yy_0[j] = pa_x[j] * tey_xxy_yy_0[j] - pc_x[j] * tey_xxy_yy_1[j] + fl1_fx * tey_xy_yy_0[j] - fl1_fx * tey_xy_yy_1[j];

                    tez_xxxy_yy_0[j] = pa_x[j] * tez_xxy_yy_0[j] - pc_x[j] * tez_xxy_yy_1[j] + fl1_fx * tez_xy_yy_0[j] - fl1_fx * tez_xy_yy_1[j];

                    tex_xxxy_yz_0[j] = pa_x[j] * tex_xxy_yz_0[j] - pc_x[j] * tex_xxy_yz_1[j] + fl1_fx * tex_xy_yz_0[j] - fl1_fx * tex_xy_yz_1[j] + ta_xxy_yz_1[j];

                    tey_xxxy_yz_0[j] = pa_x[j] * tey_xxy_yz_0[j] - pc_x[j] * tey_xxy_yz_1[j] + fl1_fx * tey_xy_yz_0[j] - fl1_fx * tey_xy_yz_1[j];

                    tez_xxxy_yz_0[j] = pa_x[j] * tez_xxy_yz_0[j] - pc_x[j] * tez_xxy_yz_1[j] + fl1_fx * tez_xy_yz_0[j] - fl1_fx * tez_xy_yz_1[j];

                    tex_xxxy_zz_0[j] = pa_x[j] * tex_xxy_zz_0[j] - pc_x[j] * tex_xxy_zz_1[j] + fl1_fx * tex_xy_zz_0[j] - fl1_fx * tex_xy_zz_1[j] + ta_xxy_zz_1[j];

                    tey_xxxy_zz_0[j] = pa_x[j] * tey_xxy_zz_0[j] - pc_x[j] * tey_xxy_zz_1[j] + fl1_fx * tey_xy_zz_0[j] - fl1_fx * tey_xy_zz_1[j];

                    tez_xxxy_zz_0[j] = pa_x[j] * tez_xxy_zz_0[j] - pc_x[j] * tez_xxy_zz_1[j] + fl1_fx * tez_xy_zz_0[j] - fl1_fx * tez_xy_zz_1[j];

                    tex_xxxz_xx_0[j] = pa_x[j] * tex_xxz_xx_0[j] - pc_x[j] * tex_xxz_xx_1[j] + fl1_fx * tex_xz_xx_0[j] - fl1_fx * tex_xz_xx_1[j] + fl1_fx * tex_xxz_x_0[j] - fl1_fx * tex_xxz_x_1[j] + ta_xxz_xx_1[j];

                    tey_xxxz_xx_0[j] = pa_x[j] * tey_xxz_xx_0[j] - pc_x[j] * tey_xxz_xx_1[j] + fl1_fx * tey_xz_xx_0[j] - fl1_fx * tey_xz_xx_1[j] + fl1_fx * tey_xxz_x_0[j] - fl1_fx * tey_xxz_x_1[j];

                    tez_xxxz_xx_0[j] = pa_x[j] * tez_xxz_xx_0[j] - pc_x[j] * tez_xxz_xx_1[j] + fl1_fx * tez_xz_xx_0[j] - fl1_fx * tez_xz_xx_1[j] + fl1_fx * tez_xxz_x_0[j] - fl1_fx * tez_xxz_x_1[j];

                    tex_xxxz_xy_0[j] = pa_x[j] * tex_xxz_xy_0[j] - pc_x[j] * tex_xxz_xy_1[j] + fl1_fx * tex_xz_xy_0[j] - fl1_fx * tex_xz_xy_1[j] + 0.5 * fl1_fx * tex_xxz_y_0[j] - 0.5 * fl1_fx * tex_xxz_y_1[j] + ta_xxz_xy_1[j];

                    tey_xxxz_xy_0[j] = pa_x[j] * tey_xxz_xy_0[j] - pc_x[j] * tey_xxz_xy_1[j] + fl1_fx * tey_xz_xy_0[j] - fl1_fx * tey_xz_xy_1[j] + 0.5 * fl1_fx * tey_xxz_y_0[j] - 0.5 * fl1_fx * tey_xxz_y_1[j];

                    tez_xxxz_xy_0[j] = pa_x[j] * tez_xxz_xy_0[j] - pc_x[j] * tez_xxz_xy_1[j] + fl1_fx * tez_xz_xy_0[j] - fl1_fx * tez_xz_xy_1[j] + 0.5 * fl1_fx * tez_xxz_y_0[j] - 0.5 * fl1_fx * tez_xxz_y_1[j];

                    tex_xxxz_xz_0[j] = pa_x[j] * tex_xxz_xz_0[j] - pc_x[j] * tex_xxz_xz_1[j] + fl1_fx * tex_xz_xz_0[j] - fl1_fx * tex_xz_xz_1[j] + 0.5 * fl1_fx * tex_xxz_z_0[j] - 0.5 * fl1_fx * tex_xxz_z_1[j] + ta_xxz_xz_1[j];

                    tey_xxxz_xz_0[j] = pa_x[j] * tey_xxz_xz_0[j] - pc_x[j] * tey_xxz_xz_1[j] + fl1_fx * tey_xz_xz_0[j] - fl1_fx * tey_xz_xz_1[j] + 0.5 * fl1_fx * tey_xxz_z_0[j] - 0.5 * fl1_fx * tey_xxz_z_1[j];

                    tez_xxxz_xz_0[j] = pa_x[j] * tez_xxz_xz_0[j] - pc_x[j] * tez_xxz_xz_1[j] + fl1_fx * tez_xz_xz_0[j] - fl1_fx * tez_xz_xz_1[j] + 0.5 * fl1_fx * tez_xxz_z_0[j] - 0.5 * fl1_fx * tez_xxz_z_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForGD_45_90(      CMemBlock2D<double>& primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto)
    {
        // Batch of Integrals (45,90)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {4, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_4_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
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

                auto tex_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 15); 

                auto tey_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 15); 

                auto tez_xxz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 15); 

                auto tex_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 16); 

                auto tey_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 16); 

                auto tez_xxz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 16); 

                auto tex_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 17); 

                auto tey_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 17); 

                auto tez_xxz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 17); 

                auto tex_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 18); 

                auto tey_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 18); 

                auto tez_xyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 18); 

                auto tex_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 19); 

                auto tey_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 19); 

                auto tez_xyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 19); 

                auto tex_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 20); 

                auto tey_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 20); 

                auto tez_xyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 20); 

                auto tex_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 21); 

                auto tey_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 21); 

                auto tez_xyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 21); 

                auto tex_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 22); 

                auto tey_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 22); 

                auto tez_xyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 22); 

                auto tex_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 23); 

                auto tey_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 23); 

                auto tez_xyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 23); 

                auto tex_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 24); 

                auto tey_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 24); 

                auto tez_xyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 24); 

                auto tex_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 25); 

                auto tey_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 25); 

                auto tez_xyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 25); 

                auto tex_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 26); 

                auto tey_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 26); 

                auto tez_xyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 26); 

                auto tex_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 27); 

                auto tey_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 27); 

                auto tez_xyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 27); 

                auto tex_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 28); 

                auto tey_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 28); 

                auto tez_xyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 28); 

                auto tex_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 29); 

                auto tey_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 29); 

                auto tez_xyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 29); 

                auto tex_xxz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 15); 

                auto tey_xxz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 15); 

                auto tez_xxz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 15); 

                auto tex_xxz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 16); 

                auto tey_xxz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 16); 

                auto tez_xxz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 16); 

                auto tex_xxz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 17); 

                auto tey_xxz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 17); 

                auto tez_xxz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 17); 

                auto tex_xyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 18); 

                auto tey_xyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 18); 

                auto tez_xyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 18); 

                auto tex_xyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 19); 

                auto tey_xyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 19); 

                auto tez_xyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 19); 

                auto tex_xyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 20); 

                auto tey_xyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 20); 

                auto tez_xyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 20); 

                auto tex_xyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 21); 

                auto tey_xyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 21); 

                auto tez_xyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 21); 

                auto tex_xyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 22); 

                auto tey_xyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 22); 

                auto tez_xyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 22); 

                auto tex_xyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 23); 

                auto tey_xyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 23); 

                auto tez_xyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 23); 

                auto tex_xyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 24); 

                auto tey_xyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 24); 

                auto tez_xyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 24); 

                auto tex_xyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 25); 

                auto tey_xyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 25); 

                auto tez_xyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 25); 

                auto tex_xyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 26); 

                auto tey_xyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 26); 

                auto tez_xyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 26); 

                auto tex_xyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 27); 

                auto tey_xyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 27); 

                auto tez_xyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 27); 

                auto tex_xyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 28); 

                auto tey_xyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 28); 

                auto tez_xyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 28); 

                auto tex_xyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 29); 

                auto tey_xyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 29); 

                auto tez_xyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 29); 

                auto tex_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 15); 

                auto tey_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 15); 

                auto tez_xz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 15); 

                auto tex_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 16); 

                auto tey_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 16); 

                auto tez_xz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 16); 

                auto tex_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 17); 

                auto tey_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 17); 

                auto tez_xz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 17); 

                auto tex_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 18); 

                auto tey_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 19); 

                auto tey_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 20); 

                auto tey_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 21); 

                auto tey_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 22); 

                auto tey_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 23); 

                auto tey_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 24); 

                auto tey_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 25); 

                auto tey_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 26); 

                auto tey_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 26); 

                auto tex_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 27); 

                auto tey_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 28); 

                auto tey_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 29); 

                auto tey_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 29); 

                auto tex_xz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 15); 

                auto tey_xz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 15); 

                auto tez_xz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 15); 

                auto tex_xz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 16); 

                auto tey_xz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 16); 

                auto tez_xz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 16); 

                auto tex_xz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 17); 

                auto tey_xz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 17); 

                auto tez_xz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 17); 

                auto tex_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 18); 

                auto tey_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 19); 

                auto tey_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 20); 

                auto tey_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 21); 

                auto tey_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 22); 

                auto tey_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 23); 

                auto tey_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 24); 

                auto tey_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 25); 

                auto tey_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 26); 

                auto tey_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 26); 

                auto tex_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 27); 

                auto tey_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 28); 

                auto tey_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 29); 

                auto tey_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 29); 

                auto tex_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 9); 

                auto tey_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 9); 

                auto tez_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 9); 

                auto tex_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 10); 

                auto tey_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 10); 

                auto tez_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 10); 

                auto tex_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 11); 

                auto tey_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 11); 

                auto tez_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 11); 

                auto tex_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 12); 

                auto tey_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 12); 

                auto tez_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 12); 

                auto tex_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 13); 

                auto tey_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 13); 

                auto tez_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 13); 

                auto tex_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 14); 

                auto tey_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 14); 

                auto tez_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 14); 

                auto tex_xyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 9); 

                auto tey_xyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 9); 

                auto tez_xyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 9); 

                auto tex_xyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 10); 

                auto tey_xyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 10); 

                auto tez_xyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 10); 

                auto tex_xyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 11); 

                auto tey_xyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 11); 

                auto tez_xyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 11); 

                auto tex_xyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 12); 

                auto tey_xyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 12); 

                auto tez_xyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 12); 

                auto tex_xyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 13); 

                auto tey_xyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 13); 

                auto tez_xyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 13); 

                auto tex_xyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 14); 

                auto tey_xyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 14); 

                auto tez_xyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 14); 

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

                auto tex_xxxz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 15); 

                auto tey_xxxz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 15); 

                auto tez_xxxz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 15); 

                auto tex_xxxz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 16); 

                auto tey_xxxz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 16); 

                auto tez_xxxz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 16); 

                auto tex_xxxz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 17); 

                auto tey_xxxz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 17); 

                auto tez_xxxz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 17); 

                auto tex_xxyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 18); 

                auto tey_xxyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 18); 

                auto tez_xxyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 18); 

                auto tex_xxyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 19); 

                auto tey_xxyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 19); 

                auto tez_xxyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 19); 

                auto tex_xxyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 20); 

                auto tey_xxyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 20); 

                auto tez_xxyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 20); 

                auto tex_xxyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 21); 

                auto tey_xxyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 21); 

                auto tez_xxyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 21); 

                auto tex_xxyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 22); 

                auto tey_xxyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 22); 

                auto tez_xxyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 22); 

                auto tex_xxyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 23); 

                auto tey_xxyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 23); 

                auto tez_xxyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 23); 

                auto tex_xxyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 24); 

                auto tey_xxyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 24); 

                auto tez_xxyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 24); 

                auto tex_xxyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 25); 

                auto tey_xxyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 25); 

                auto tez_xxyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 25); 

                auto tex_xxyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 26); 

                auto tey_xxyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 26); 

                auto tez_xxyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 26); 

                auto tex_xxyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 27); 

                auto tey_xxyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 27); 

                auto tez_xxyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 27); 

                auto tex_xxyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 28); 

                auto tey_xxyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 28); 

                auto tez_xxyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 28); 

                auto tex_xxyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 29); 

                auto tey_xxyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 29); 

                auto tez_xxyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 29); 

                // Batch of Integrals (45,90)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxz_yy_1, ta_xxz_yz_1, ta_xxz_zz_1, ta_xyy_xx_1, \
                                         ta_xyy_xy_1, ta_xyy_xz_1, ta_xyy_yy_1, ta_xyy_yz_1, ta_xyy_zz_1, ta_xyz_xx_1, \
                                         ta_xyz_xy_1, ta_xyz_xz_1, ta_xyz_yy_1, ta_xyz_yz_1, ta_xyz_zz_1, tex_xxxz_yy_0, \
                                         tex_xxxz_yz_0, tex_xxxz_zz_0, tex_xxyy_xx_0, tex_xxyy_xy_0, tex_xxyy_xz_0, \
                                         tex_xxyy_yy_0, tex_xxyy_yz_0, tex_xxyy_zz_0, tex_xxyz_xx_0, tex_xxyz_xy_0, \
                                         tex_xxyz_xz_0, tex_xxyz_yy_0, tex_xxyz_yz_0, tex_xxyz_zz_0, tex_xxz_yy_0, \
                                         tex_xxz_yy_1, tex_xxz_yz_0, tex_xxz_yz_1, tex_xxz_zz_0, tex_xxz_zz_1, tex_xyy_x_0, \
                                         tex_xyy_x_1, tex_xyy_xx_0, tex_xyy_xx_1, tex_xyy_xy_0, tex_xyy_xy_1, tex_xyy_xz_0, \
                                         tex_xyy_xz_1, tex_xyy_y_0, tex_xyy_y_1, tex_xyy_yy_0, tex_xyy_yy_1, tex_xyy_yz_0, \
                                         tex_xyy_yz_1, tex_xyy_z_0, tex_xyy_z_1, tex_xyy_zz_0, tex_xyy_zz_1, tex_xyz_x_0, \
                                         tex_xyz_x_1, tex_xyz_xx_0, tex_xyz_xx_1, tex_xyz_xy_0, tex_xyz_xy_1, tex_xyz_xz_0, \
                                         tex_xyz_xz_1, tex_xyz_y_0, tex_xyz_y_1, tex_xyz_yy_0, tex_xyz_yy_1, tex_xyz_yz_0, \
                                         tex_xyz_yz_1, tex_xyz_z_0, tex_xyz_z_1, tex_xyz_zz_0, tex_xyz_zz_1, tex_xz_yy_0, \
                                         tex_xz_yy_1, tex_xz_yz_0, tex_xz_yz_1, tex_xz_zz_0, tex_xz_zz_1, tex_yy_xx_0, \
                                         tex_yy_xx_1, tex_yy_xy_0, tex_yy_xy_1, tex_yy_xz_0, tex_yy_xz_1, tex_yy_yy_0, \
                                         tex_yy_yy_1, tex_yy_yz_0, tex_yy_yz_1, tex_yy_zz_0, tex_yy_zz_1, tex_yz_xx_0, \
                                         tex_yz_xx_1, tex_yz_xy_0, tex_yz_xy_1, tex_yz_xz_0, tex_yz_xz_1, tex_yz_yy_0, \
                                         tex_yz_yy_1, tex_yz_yz_0, tex_yz_yz_1, tex_yz_zz_0, tex_yz_zz_1, tey_xxxz_yy_0, \
                                         tey_xxxz_yz_0, tey_xxxz_zz_0, tey_xxyy_xx_0, tey_xxyy_xy_0, tey_xxyy_xz_0, \
                                         tey_xxyy_yy_0, tey_xxyy_yz_0, tey_xxyy_zz_0, tey_xxyz_xx_0, tey_xxyz_xy_0, \
                                         tey_xxyz_xz_0, tey_xxyz_yy_0, tey_xxyz_yz_0, tey_xxyz_zz_0, tey_xxz_yy_0, \
                                         tey_xxz_yy_1, tey_xxz_yz_0, tey_xxz_yz_1, tey_xxz_zz_0, tey_xxz_zz_1, tey_xyy_x_0, \
                                         tey_xyy_x_1, tey_xyy_xx_0, tey_xyy_xx_1, tey_xyy_xy_0, tey_xyy_xy_1, tey_xyy_xz_0, \
                                         tey_xyy_xz_1, tey_xyy_y_0, tey_xyy_y_1, tey_xyy_yy_0, tey_xyy_yy_1, tey_xyy_yz_0, \
                                         tey_xyy_yz_1, tey_xyy_z_0, tey_xyy_z_1, tey_xyy_zz_0, tey_xyy_zz_1, tey_xyz_x_0, \
                                         tey_xyz_x_1, tey_xyz_xx_0, tey_xyz_xx_1, tey_xyz_xy_0, tey_xyz_xy_1, tey_xyz_xz_0, \
                                         tey_xyz_xz_1, tey_xyz_y_0, tey_xyz_y_1, tey_xyz_yy_0, tey_xyz_yy_1, tey_xyz_yz_0, \
                                         tey_xyz_yz_1, tey_xyz_z_0, tey_xyz_z_1, tey_xyz_zz_0, tey_xyz_zz_1, tey_xz_yy_0, \
                                         tey_xz_yy_1, tey_xz_yz_0, tey_xz_yz_1, tey_xz_zz_0, tey_xz_zz_1, tey_yy_xx_0, \
                                         tey_yy_xx_1, tey_yy_xy_0, tey_yy_xy_1, tey_yy_xz_0, tey_yy_xz_1, tey_yy_yy_0, \
                                         tey_yy_yy_1, tey_yy_yz_0, tey_yy_yz_1, tey_yy_zz_0, tey_yy_zz_1, tey_yz_xx_0, \
                                         tey_yz_xx_1, tey_yz_xy_0, tey_yz_xy_1, tey_yz_xz_0, tey_yz_xz_1, tey_yz_yy_0, \
                                         tey_yz_yy_1, tey_yz_yz_0, tey_yz_yz_1, tey_yz_zz_0, tey_yz_zz_1, tez_xxxz_yy_0, \
                                         tez_xxxz_yz_0, tez_xxxz_zz_0, tez_xxyy_xx_0, tez_xxyy_xy_0, tez_xxyy_xz_0, \
                                         tez_xxyy_yy_0, tez_xxyy_yz_0, tez_xxyy_zz_0, tez_xxyz_xx_0, tez_xxyz_xy_0, \
                                         tez_xxyz_xz_0, tez_xxyz_yy_0, tez_xxyz_yz_0, tez_xxyz_zz_0, tez_xxz_yy_0, \
                                         tez_xxz_yy_1, tez_xxz_yz_0, tez_xxz_yz_1, tez_xxz_zz_0, tez_xxz_zz_1, tez_xyy_x_0, \
                                         tez_xyy_x_1, tez_xyy_xx_0, tez_xyy_xx_1, tez_xyy_xy_0, tez_xyy_xy_1, tez_xyy_xz_0, \
                                         tez_xyy_xz_1, tez_xyy_y_0, tez_xyy_y_1, tez_xyy_yy_0, tez_xyy_yy_1, tez_xyy_yz_0, \
                                         tez_xyy_yz_1, tez_xyy_z_0, tez_xyy_z_1, tez_xyy_zz_0, tez_xyy_zz_1, tez_xyz_x_0, \
                                         tez_xyz_x_1, tez_xyz_xx_0, tez_xyz_xx_1, tez_xyz_xy_0, tez_xyz_xy_1, tez_xyz_xz_0, \
                                         tez_xyz_xz_1, tez_xyz_y_0, tez_xyz_y_1, tez_xyz_yy_0, tez_xyz_yy_1, tez_xyz_yz_0, \
                                         tez_xyz_yz_1, tez_xyz_z_0, tez_xyz_z_1, tez_xyz_zz_0, tez_xyz_zz_1, tez_xz_yy_0, \
                                         tez_xz_yy_1, tez_xz_yz_0, tez_xz_yz_1, tez_xz_zz_0, tez_xz_zz_1, tez_yy_xx_0, \
                                         tez_yy_xx_1, tez_yy_xy_0, tez_yy_xy_1, tez_yy_xz_0, tez_yy_xz_1, tez_yy_yy_0, \
                                         tez_yy_yy_1, tez_yy_yz_0, tez_yy_yz_1, tez_yy_zz_0, tez_yy_zz_1, tez_yz_xx_0, \
                                         tez_yz_xx_1, tez_yz_xy_0, tez_yz_xy_1, tez_yz_xz_0, tez_yz_xz_1, tez_yz_yy_0, \
                                         tez_yz_yy_1, tez_yz_yz_0, tez_yz_yz_1, tez_yz_zz_0, tez_yz_zz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xxxz_yy_0[j] = pa_x[j] * tex_xxz_yy_0[j] - pc_x[j] * tex_xxz_yy_1[j] + fl1_fx * tex_xz_yy_0[j] - fl1_fx * tex_xz_yy_1[j] + ta_xxz_yy_1[j];

                    tey_xxxz_yy_0[j] = pa_x[j] * tey_xxz_yy_0[j] - pc_x[j] * tey_xxz_yy_1[j] + fl1_fx * tey_xz_yy_0[j] - fl1_fx * tey_xz_yy_1[j];

                    tez_xxxz_yy_0[j] = pa_x[j] * tez_xxz_yy_0[j] - pc_x[j] * tez_xxz_yy_1[j] + fl1_fx * tez_xz_yy_0[j] - fl1_fx * tez_xz_yy_1[j];

                    tex_xxxz_yz_0[j] = pa_x[j] * tex_xxz_yz_0[j] - pc_x[j] * tex_xxz_yz_1[j] + fl1_fx * tex_xz_yz_0[j] - fl1_fx * tex_xz_yz_1[j] + ta_xxz_yz_1[j];

                    tey_xxxz_yz_0[j] = pa_x[j] * tey_xxz_yz_0[j] - pc_x[j] * tey_xxz_yz_1[j] + fl1_fx * tey_xz_yz_0[j] - fl1_fx * tey_xz_yz_1[j];

                    tez_xxxz_yz_0[j] = pa_x[j] * tez_xxz_yz_0[j] - pc_x[j] * tez_xxz_yz_1[j] + fl1_fx * tez_xz_yz_0[j] - fl1_fx * tez_xz_yz_1[j];

                    tex_xxxz_zz_0[j] = pa_x[j] * tex_xxz_zz_0[j] - pc_x[j] * tex_xxz_zz_1[j] + fl1_fx * tex_xz_zz_0[j] - fl1_fx * tex_xz_zz_1[j] + ta_xxz_zz_1[j];

                    tey_xxxz_zz_0[j] = pa_x[j] * tey_xxz_zz_0[j] - pc_x[j] * tey_xxz_zz_1[j] + fl1_fx * tey_xz_zz_0[j] - fl1_fx * tey_xz_zz_1[j];

                    tez_xxxz_zz_0[j] = pa_x[j] * tez_xxz_zz_0[j] - pc_x[j] * tez_xxz_zz_1[j] + fl1_fx * tez_xz_zz_0[j] - fl1_fx * tez_xz_zz_1[j];

                    tex_xxyy_xx_0[j] = pa_x[j] * tex_xyy_xx_0[j] - pc_x[j] * tex_xyy_xx_1[j] + 0.5 * fl1_fx * tex_yy_xx_0[j] - 0.5 * fl1_fx * tex_yy_xx_1[j] + fl1_fx * tex_xyy_x_0[j] - fl1_fx * tex_xyy_x_1[j] + ta_xyy_xx_1[j];

                    tey_xxyy_xx_0[j] = pa_x[j] * tey_xyy_xx_0[j] - pc_x[j] * tey_xyy_xx_1[j] + 0.5 * fl1_fx * tey_yy_xx_0[j] - 0.5 * fl1_fx * tey_yy_xx_1[j] + fl1_fx * tey_xyy_x_0[j] - fl1_fx * tey_xyy_x_1[j];

                    tez_xxyy_xx_0[j] = pa_x[j] * tez_xyy_xx_0[j] - pc_x[j] * tez_xyy_xx_1[j] + 0.5 * fl1_fx * tez_yy_xx_0[j] - 0.5 * fl1_fx * tez_yy_xx_1[j] + fl1_fx * tez_xyy_x_0[j] - fl1_fx * tez_xyy_x_1[j];

                    tex_xxyy_xy_0[j] = pa_x[j] * tex_xyy_xy_0[j] - pc_x[j] * tex_xyy_xy_1[j] + 0.5 * fl1_fx * tex_yy_xy_0[j] - 0.5 * fl1_fx * tex_yy_xy_1[j] + 0.5 * fl1_fx * tex_xyy_y_0[j] - 0.5 * fl1_fx * tex_xyy_y_1[j] + ta_xyy_xy_1[j];

                    tey_xxyy_xy_0[j] = pa_x[j] * tey_xyy_xy_0[j] - pc_x[j] * tey_xyy_xy_1[j] + 0.5 * fl1_fx * tey_yy_xy_0[j] - 0.5 * fl1_fx * tey_yy_xy_1[j] + 0.5 * fl1_fx * tey_xyy_y_0[j] - 0.5 * fl1_fx * tey_xyy_y_1[j];

                    tez_xxyy_xy_0[j] = pa_x[j] * tez_xyy_xy_0[j] - pc_x[j] * tez_xyy_xy_1[j] + 0.5 * fl1_fx * tez_yy_xy_0[j] - 0.5 * fl1_fx * tez_yy_xy_1[j] + 0.5 * fl1_fx * tez_xyy_y_0[j] - 0.5 * fl1_fx * tez_xyy_y_1[j];

                    tex_xxyy_xz_0[j] = pa_x[j] * tex_xyy_xz_0[j] - pc_x[j] * tex_xyy_xz_1[j] + 0.5 * fl1_fx * tex_yy_xz_0[j] - 0.5 * fl1_fx * tex_yy_xz_1[j] + 0.5 * fl1_fx * tex_xyy_z_0[j] - 0.5 * fl1_fx * tex_xyy_z_1[j] + ta_xyy_xz_1[j];

                    tey_xxyy_xz_0[j] = pa_x[j] * tey_xyy_xz_0[j] - pc_x[j] * tey_xyy_xz_1[j] + 0.5 * fl1_fx * tey_yy_xz_0[j] - 0.5 * fl1_fx * tey_yy_xz_1[j] + 0.5 * fl1_fx * tey_xyy_z_0[j] - 0.5 * fl1_fx * tey_xyy_z_1[j];

                    tez_xxyy_xz_0[j] = pa_x[j] * tez_xyy_xz_0[j] - pc_x[j] * tez_xyy_xz_1[j] + 0.5 * fl1_fx * tez_yy_xz_0[j] - 0.5 * fl1_fx * tez_yy_xz_1[j] + 0.5 * fl1_fx * tez_xyy_z_0[j] - 0.5 * fl1_fx * tez_xyy_z_1[j];

                    tex_xxyy_yy_0[j] = pa_x[j] * tex_xyy_yy_0[j] - pc_x[j] * tex_xyy_yy_1[j] + 0.5 * fl1_fx * tex_yy_yy_0[j] - 0.5 * fl1_fx * tex_yy_yy_1[j] + ta_xyy_yy_1[j];

                    tey_xxyy_yy_0[j] = pa_x[j] * tey_xyy_yy_0[j] - pc_x[j] * tey_xyy_yy_1[j] + 0.5 * fl1_fx * tey_yy_yy_0[j] - 0.5 * fl1_fx * tey_yy_yy_1[j];

                    tez_xxyy_yy_0[j] = pa_x[j] * tez_xyy_yy_0[j] - pc_x[j] * tez_xyy_yy_1[j] + 0.5 * fl1_fx * tez_yy_yy_0[j] - 0.5 * fl1_fx * tez_yy_yy_1[j];

                    tex_xxyy_yz_0[j] = pa_x[j] * tex_xyy_yz_0[j] - pc_x[j] * tex_xyy_yz_1[j] + 0.5 * fl1_fx * tex_yy_yz_0[j] - 0.5 * fl1_fx * tex_yy_yz_1[j] + ta_xyy_yz_1[j];

                    tey_xxyy_yz_0[j] = pa_x[j] * tey_xyy_yz_0[j] - pc_x[j] * tey_xyy_yz_1[j] + 0.5 * fl1_fx * tey_yy_yz_0[j] - 0.5 * fl1_fx * tey_yy_yz_1[j];

                    tez_xxyy_yz_0[j] = pa_x[j] * tez_xyy_yz_0[j] - pc_x[j] * tez_xyy_yz_1[j] + 0.5 * fl1_fx * tez_yy_yz_0[j] - 0.5 * fl1_fx * tez_yy_yz_1[j];

                    tex_xxyy_zz_0[j] = pa_x[j] * tex_xyy_zz_0[j] - pc_x[j] * tex_xyy_zz_1[j] + 0.5 * fl1_fx * tex_yy_zz_0[j] - 0.5 * fl1_fx * tex_yy_zz_1[j] + ta_xyy_zz_1[j];

                    tey_xxyy_zz_0[j] = pa_x[j] * tey_xyy_zz_0[j] - pc_x[j] * tey_xyy_zz_1[j] + 0.5 * fl1_fx * tey_yy_zz_0[j] - 0.5 * fl1_fx * tey_yy_zz_1[j];

                    tez_xxyy_zz_0[j] = pa_x[j] * tez_xyy_zz_0[j] - pc_x[j] * tez_xyy_zz_1[j] + 0.5 * fl1_fx * tez_yy_zz_0[j] - 0.5 * fl1_fx * tez_yy_zz_1[j];

                    tex_xxyz_xx_0[j] = pa_x[j] * tex_xyz_xx_0[j] - pc_x[j] * tex_xyz_xx_1[j] + 0.5 * fl1_fx * tex_yz_xx_0[j] - 0.5 * fl1_fx * tex_yz_xx_1[j] + fl1_fx * tex_xyz_x_0[j] - fl1_fx * tex_xyz_x_1[j] + ta_xyz_xx_1[j];

                    tey_xxyz_xx_0[j] = pa_x[j] * tey_xyz_xx_0[j] - pc_x[j] * tey_xyz_xx_1[j] + 0.5 * fl1_fx * tey_yz_xx_0[j] - 0.5 * fl1_fx * tey_yz_xx_1[j] + fl1_fx * tey_xyz_x_0[j] - fl1_fx * tey_xyz_x_1[j];

                    tez_xxyz_xx_0[j] = pa_x[j] * tez_xyz_xx_0[j] - pc_x[j] * tez_xyz_xx_1[j] + 0.5 * fl1_fx * tez_yz_xx_0[j] - 0.5 * fl1_fx * tez_yz_xx_1[j] + fl1_fx * tez_xyz_x_0[j] - fl1_fx * tez_xyz_x_1[j];

                    tex_xxyz_xy_0[j] = pa_x[j] * tex_xyz_xy_0[j] - pc_x[j] * tex_xyz_xy_1[j] + 0.5 * fl1_fx * tex_yz_xy_0[j] - 0.5 * fl1_fx * tex_yz_xy_1[j] + 0.5 * fl1_fx * tex_xyz_y_0[j] - 0.5 * fl1_fx * tex_xyz_y_1[j] + ta_xyz_xy_1[j];

                    tey_xxyz_xy_0[j] = pa_x[j] * tey_xyz_xy_0[j] - pc_x[j] * tey_xyz_xy_1[j] + 0.5 * fl1_fx * tey_yz_xy_0[j] - 0.5 * fl1_fx * tey_yz_xy_1[j] + 0.5 * fl1_fx * tey_xyz_y_0[j] - 0.5 * fl1_fx * tey_xyz_y_1[j];

                    tez_xxyz_xy_0[j] = pa_x[j] * tez_xyz_xy_0[j] - pc_x[j] * tez_xyz_xy_1[j] + 0.5 * fl1_fx * tez_yz_xy_0[j] - 0.5 * fl1_fx * tez_yz_xy_1[j] + 0.5 * fl1_fx * tez_xyz_y_0[j] - 0.5 * fl1_fx * tez_xyz_y_1[j];

                    tex_xxyz_xz_0[j] = pa_x[j] * tex_xyz_xz_0[j] - pc_x[j] * tex_xyz_xz_1[j] + 0.5 * fl1_fx * tex_yz_xz_0[j] - 0.5 * fl1_fx * tex_yz_xz_1[j] + 0.5 * fl1_fx * tex_xyz_z_0[j] - 0.5 * fl1_fx * tex_xyz_z_1[j] + ta_xyz_xz_1[j];

                    tey_xxyz_xz_0[j] = pa_x[j] * tey_xyz_xz_0[j] - pc_x[j] * tey_xyz_xz_1[j] + 0.5 * fl1_fx * tey_yz_xz_0[j] - 0.5 * fl1_fx * tey_yz_xz_1[j] + 0.5 * fl1_fx * tey_xyz_z_0[j] - 0.5 * fl1_fx * tey_xyz_z_1[j];

                    tez_xxyz_xz_0[j] = pa_x[j] * tez_xyz_xz_0[j] - pc_x[j] * tez_xyz_xz_1[j] + 0.5 * fl1_fx * tez_yz_xz_0[j] - 0.5 * fl1_fx * tez_yz_xz_1[j] + 0.5 * fl1_fx * tez_xyz_z_0[j] - 0.5 * fl1_fx * tez_xyz_z_1[j];

                    tex_xxyz_yy_0[j] = pa_x[j] * tex_xyz_yy_0[j] - pc_x[j] * tex_xyz_yy_1[j] + 0.5 * fl1_fx * tex_yz_yy_0[j] - 0.5 * fl1_fx * tex_yz_yy_1[j] + ta_xyz_yy_1[j];

                    tey_xxyz_yy_0[j] = pa_x[j] * tey_xyz_yy_0[j] - pc_x[j] * tey_xyz_yy_1[j] + 0.5 * fl1_fx * tey_yz_yy_0[j] - 0.5 * fl1_fx * tey_yz_yy_1[j];

                    tez_xxyz_yy_0[j] = pa_x[j] * tez_xyz_yy_0[j] - pc_x[j] * tez_xyz_yy_1[j] + 0.5 * fl1_fx * tez_yz_yy_0[j] - 0.5 * fl1_fx * tez_yz_yy_1[j];

                    tex_xxyz_yz_0[j] = pa_x[j] * tex_xyz_yz_0[j] - pc_x[j] * tex_xyz_yz_1[j] + 0.5 * fl1_fx * tex_yz_yz_0[j] - 0.5 * fl1_fx * tex_yz_yz_1[j] + ta_xyz_yz_1[j];

                    tey_xxyz_yz_0[j] = pa_x[j] * tey_xyz_yz_0[j] - pc_x[j] * tey_xyz_yz_1[j] + 0.5 * fl1_fx * tey_yz_yz_0[j] - 0.5 * fl1_fx * tey_yz_yz_1[j];

                    tez_xxyz_yz_0[j] = pa_x[j] * tez_xyz_yz_0[j] - pc_x[j] * tez_xyz_yz_1[j] + 0.5 * fl1_fx * tez_yz_yz_0[j] - 0.5 * fl1_fx * tez_yz_yz_1[j];

                    tex_xxyz_zz_0[j] = pa_x[j] * tex_xyz_zz_0[j] - pc_x[j] * tex_xyz_zz_1[j] + 0.5 * fl1_fx * tex_yz_zz_0[j] - 0.5 * fl1_fx * tex_yz_zz_1[j] + ta_xyz_zz_1[j];

                    tey_xxyz_zz_0[j] = pa_x[j] * tey_xyz_zz_0[j] - pc_x[j] * tey_xyz_zz_1[j] + 0.5 * fl1_fx * tey_yz_zz_0[j] - 0.5 * fl1_fx * tey_yz_zz_1[j];

                    tez_xxyz_zz_0[j] = pa_x[j] * tez_xyz_zz_0[j] - pc_x[j] * tez_xyz_zz_1[j] + 0.5 * fl1_fx * tez_yz_zz_0[j] - 0.5 * fl1_fx * tez_yz_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForGD_90_135(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
                                  const CGtoBlock&           braGtoBlock,
                                  const CGtoBlock&           ketGtoBlock,
                                  const int32_t              iContrGto)
    {
        // Batch of Integrals (90,135)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {4, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_4_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                auto tex_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 30); 

                auto tey_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 30); 

                auto tez_xzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 30); 

                auto tex_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 31); 

                auto tey_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 31); 

                auto tez_xzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 31); 

                auto tex_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 32); 

                auto tey_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 32); 

                auto tez_xzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 32); 

                auto tex_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 33); 

                auto tey_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 33); 

                auto tez_xzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 33); 

                auto tex_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 34); 

                auto tey_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 34); 

                auto tez_xzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 34); 

                auto tex_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 35); 

                auto tey_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 35); 

                auto tez_xzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 35); 

                auto tex_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 36); 

                auto tey_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 36); 

                auto tez_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 36); 

                auto tex_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 37); 

                auto tey_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 37); 

                auto tez_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 37); 

                auto tex_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 38); 

                auto tey_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 38); 

                auto tez_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 38); 

                auto tex_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 39); 

                auto tey_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 39); 

                auto tez_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 39); 

                auto tex_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 40); 

                auto tey_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 40); 

                auto tez_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 40); 

                auto tex_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 41); 

                auto tey_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 41); 

                auto tez_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 41); 

                auto tex_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 42); 

                auto tey_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 42); 

                auto tez_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 42); 

                auto tex_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 43); 

                auto tey_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 43); 

                auto tez_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 43); 

                auto tex_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 44); 

                auto tey_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 44); 

                auto tez_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 44); 

                auto tex_xzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 30); 

                auto tey_xzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 30); 

                auto tez_xzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 30); 

                auto tex_xzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 31); 

                auto tey_xzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 31); 

                auto tez_xzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 31); 

                auto tex_xzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 32); 

                auto tey_xzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 32); 

                auto tez_xzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 32); 

                auto tex_xzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 33); 

                auto tey_xzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 33); 

                auto tez_xzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 33); 

                auto tex_xzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 34); 

                auto tey_xzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 34); 

                auto tez_xzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 34); 

                auto tex_xzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 35); 

                auto tey_xzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 35); 

                auto tez_xzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 35); 

                auto tex_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 36); 

                auto tey_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 36); 

                auto tez_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 36); 

                auto tex_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 37); 

                auto tey_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 37); 

                auto tez_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 37); 

                auto tex_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 38); 

                auto tey_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 38); 

                auto tez_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 38); 

                auto tex_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 39); 

                auto tey_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 39); 

                auto tez_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 39); 

                auto tex_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 40); 

                auto tey_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 40); 

                auto tez_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 40); 

                auto tex_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 41); 

                auto tey_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 41); 

                auto tez_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 41); 

                auto tex_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 42); 

                auto tey_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 42); 

                auto tez_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 42); 

                auto tex_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 43); 

                auto tey_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 43); 

                auto tez_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 43); 

                auto tex_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 44); 

                auto tey_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 44); 

                auto tez_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 44); 

                auto tex_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 30); 

                auto tey_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 31); 

                auto tey_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 32); 

                auto tey_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 33); 

                auto tey_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 34); 

                auto tey_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 35); 

                auto tey_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 35); 

                auto tex_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 30); 

                auto tey_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 31); 

                auto tey_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 32); 

                auto tey_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 33); 

                auto tey_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 34); 

                auto tey_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 35); 

                auto tey_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 35); 

                auto tex_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 15); 

                auto tey_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 15); 

                auto tez_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 15); 

                auto tex_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 16); 

                auto tey_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 16); 

                auto tez_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 16); 

                auto tex_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 17); 

                auto tey_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 17); 

                auto tez_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 17); 

                auto tex_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 18); 

                auto tey_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 18); 

                auto tez_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 18); 

                auto tex_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 19); 

                auto tey_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 19); 

                auto tez_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 19); 

                auto tex_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 20); 

                auto tey_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 21); 

                auto tey_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 22); 

                auto tey_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 23); 

                auto tey_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_xzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 15); 

                auto tey_xzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 15); 

                auto tez_xzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 15); 

                auto tex_xzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 16); 

                auto tey_xzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 16); 

                auto tez_xzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 16); 

                auto tex_xzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 17); 

                auto tey_xzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 17); 

                auto tez_xzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 17); 

                auto tex_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 18); 

                auto tey_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 18); 

                auto tez_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 18); 

                auto tex_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 19); 

                auto tey_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 19); 

                auto tez_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 19); 

                auto tex_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 20); 

                auto tey_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 21); 

                auto tey_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 22); 

                auto tey_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 23); 

                auto tey_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 23); 

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

                // set up pointers to integrals

                auto tex_xxzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 30); 

                auto tey_xxzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 30); 

                auto tez_xxzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 30); 

                auto tex_xxzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 31); 

                auto tey_xxzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 31); 

                auto tez_xxzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 31); 

                auto tex_xxzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 32); 

                auto tey_xxzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 32); 

                auto tez_xxzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 32); 

                auto tex_xxzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 33); 

                auto tey_xxzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 33); 

                auto tez_xxzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 33); 

                auto tex_xxzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 34); 

                auto tey_xxzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 34); 

                auto tez_xxzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 34); 

                auto tex_xxzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 35); 

                auto tey_xxzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 35); 

                auto tez_xxzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 35); 

                auto tex_xyyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 36); 

                auto tey_xyyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 36); 

                auto tez_xyyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 36); 

                auto tex_xyyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 37); 

                auto tey_xyyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 37); 

                auto tez_xyyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 37); 

                auto tex_xyyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 38); 

                auto tey_xyyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 38); 

                auto tez_xyyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 38); 

                auto tex_xyyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 39); 

                auto tey_xyyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 39); 

                auto tez_xyyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 39); 

                auto tex_xyyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 40); 

                auto tey_xyyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 40); 

                auto tez_xyyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 40); 

                auto tex_xyyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 41); 

                auto tey_xyyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 41); 

                auto tez_xyyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 41); 

                auto tex_xyyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 42); 

                auto tey_xyyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 42); 

                auto tez_xyyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 42); 

                auto tex_xyyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 43); 

                auto tey_xyyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 43); 

                auto tez_xyyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 43); 

                auto tex_xyyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 44); 

                auto tey_xyyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 44); 

                auto tez_xyyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 44); 

                // Batch of Integrals (90,135)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_xzz_xx_1, ta_xzz_xy_1, ta_xzz_xz_1, ta_xzz_yy_1, \
                                         ta_xzz_yz_1, ta_xzz_zz_1, ta_yyy_xx_1, ta_yyy_xy_1, ta_yyy_xz_1, ta_yyy_yy_1, \
                                         ta_yyy_yz_1, ta_yyy_zz_1, ta_yyz_xx_1, ta_yyz_xy_1, ta_yyz_xz_1, tex_xxzz_xx_0, \
                                         tex_xxzz_xy_0, tex_xxzz_xz_0, tex_xxzz_yy_0, tex_xxzz_yz_0, tex_xxzz_zz_0, \
                                         tex_xyyy_xx_0, tex_xyyy_xy_0, tex_xyyy_xz_0, tex_xyyy_yy_0, tex_xyyy_yz_0, \
                                         tex_xyyy_zz_0, tex_xyyz_xx_0, tex_xyyz_xy_0, tex_xyyz_xz_0, tex_xzz_x_0, \
                                         tex_xzz_x_1, tex_xzz_xx_0, tex_xzz_xx_1, tex_xzz_xy_0, tex_xzz_xy_1, tex_xzz_xz_0, \
                                         tex_xzz_xz_1, tex_xzz_y_0, tex_xzz_y_1, tex_xzz_yy_0, tex_xzz_yy_1, tex_xzz_yz_0, \
                                         tex_xzz_yz_1, tex_xzz_z_0, tex_xzz_z_1, tex_xzz_zz_0, tex_xzz_zz_1, tex_yyy_x_0, \
                                         tex_yyy_x_1, tex_yyy_xx_0, tex_yyy_xx_1, tex_yyy_xy_0, tex_yyy_xy_1, tex_yyy_xz_0, \
                                         tex_yyy_xz_1, tex_yyy_y_0, tex_yyy_y_1, tex_yyy_yy_0, tex_yyy_yy_1, tex_yyy_yz_0, \
                                         tex_yyy_yz_1, tex_yyy_z_0, tex_yyy_z_1, tex_yyy_zz_0, tex_yyy_zz_1, tex_yyz_x_0, \
                                         tex_yyz_x_1, tex_yyz_xx_0, tex_yyz_xx_1, tex_yyz_xy_0, tex_yyz_xy_1, tex_yyz_xz_0, \
                                         tex_yyz_xz_1, tex_yyz_y_0, tex_yyz_y_1, tex_yyz_z_0, tex_yyz_z_1, tex_zz_xx_0, \
                                         tex_zz_xx_1, tex_zz_xy_0, tex_zz_xy_1, tex_zz_xz_0, tex_zz_xz_1, tex_zz_yy_0, \
                                         tex_zz_yy_1, tex_zz_yz_0, tex_zz_yz_1, tex_zz_zz_0, tex_zz_zz_1, tey_xxzz_xx_0, \
                                         tey_xxzz_xy_0, tey_xxzz_xz_0, tey_xxzz_yy_0, tey_xxzz_yz_0, tey_xxzz_zz_0, \
                                         tey_xyyy_xx_0, tey_xyyy_xy_0, tey_xyyy_xz_0, tey_xyyy_yy_0, tey_xyyy_yz_0, \
                                         tey_xyyy_zz_0, tey_xyyz_xx_0, tey_xyyz_xy_0, tey_xyyz_xz_0, tey_xzz_x_0, \
                                         tey_xzz_x_1, tey_xzz_xx_0, tey_xzz_xx_1, tey_xzz_xy_0, tey_xzz_xy_1, tey_xzz_xz_0, \
                                         tey_xzz_xz_1, tey_xzz_y_0, tey_xzz_y_1, tey_xzz_yy_0, tey_xzz_yy_1, tey_xzz_yz_0, \
                                         tey_xzz_yz_1, tey_xzz_z_0, tey_xzz_z_1, tey_xzz_zz_0, tey_xzz_zz_1, tey_yyy_x_0, \
                                         tey_yyy_x_1, tey_yyy_xx_0, tey_yyy_xx_1, tey_yyy_xy_0, tey_yyy_xy_1, tey_yyy_xz_0, \
                                         tey_yyy_xz_1, tey_yyy_y_0, tey_yyy_y_1, tey_yyy_yy_0, tey_yyy_yy_1, tey_yyy_yz_0, \
                                         tey_yyy_yz_1, tey_yyy_z_0, tey_yyy_z_1, tey_yyy_zz_0, tey_yyy_zz_1, tey_yyz_x_0, \
                                         tey_yyz_x_1, tey_yyz_xx_0, tey_yyz_xx_1, tey_yyz_xy_0, tey_yyz_xy_1, tey_yyz_xz_0, \
                                         tey_yyz_xz_1, tey_yyz_y_0, tey_yyz_y_1, tey_yyz_z_0, tey_yyz_z_1, tey_zz_xx_0, \
                                         tey_zz_xx_1, tey_zz_xy_0, tey_zz_xy_1, tey_zz_xz_0, tey_zz_xz_1, tey_zz_yy_0, \
                                         tey_zz_yy_1, tey_zz_yz_0, tey_zz_yz_1, tey_zz_zz_0, tey_zz_zz_1, tez_xxzz_xx_0, \
                                         tez_xxzz_xy_0, tez_xxzz_xz_0, tez_xxzz_yy_0, tez_xxzz_yz_0, tez_xxzz_zz_0, \
                                         tez_xyyy_xx_0, tez_xyyy_xy_0, tez_xyyy_xz_0, tez_xyyy_yy_0, tez_xyyy_yz_0, \
                                         tez_xyyy_zz_0, tez_xyyz_xx_0, tez_xyyz_xy_0, tez_xyyz_xz_0, tez_xzz_x_0, \
                                         tez_xzz_x_1, tez_xzz_xx_0, tez_xzz_xx_1, tez_xzz_xy_0, tez_xzz_xy_1, tez_xzz_xz_0, \
                                         tez_xzz_xz_1, tez_xzz_y_0, tez_xzz_y_1, tez_xzz_yy_0, tez_xzz_yy_1, tez_xzz_yz_0, \
                                         tez_xzz_yz_1, tez_xzz_z_0, tez_xzz_z_1, tez_xzz_zz_0, tez_xzz_zz_1, tez_yyy_x_0, \
                                         tez_yyy_x_1, tez_yyy_xx_0, tez_yyy_xx_1, tez_yyy_xy_0, tez_yyy_xy_1, tez_yyy_xz_0, \
                                         tez_yyy_xz_1, tez_yyy_y_0, tez_yyy_y_1, tez_yyy_yy_0, tez_yyy_yy_1, tez_yyy_yz_0, \
                                         tez_yyy_yz_1, tez_yyy_z_0, tez_yyy_z_1, tez_yyy_zz_0, tez_yyy_zz_1, tez_yyz_x_0, \
                                         tez_yyz_x_1, tez_yyz_xx_0, tez_yyz_xx_1, tez_yyz_xy_0, tez_yyz_xy_1, tez_yyz_xz_0, \
                                         tez_yyz_xz_1, tez_yyz_y_0, tez_yyz_y_1, tez_yyz_z_0, tez_yyz_z_1, tez_zz_xx_0, \
                                         tez_zz_xx_1, tez_zz_xy_0, tez_zz_xy_1, tez_zz_xz_0, tez_zz_xz_1, tez_zz_yy_0, \
                                         tez_zz_yy_1, tez_zz_yz_0, tez_zz_yz_1, tez_zz_zz_0, tez_zz_zz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xxzz_xx_0[j] = pa_x[j] * tex_xzz_xx_0[j] - pc_x[j] * tex_xzz_xx_1[j] + 0.5 * fl1_fx * tex_zz_xx_0[j] - 0.5 * fl1_fx * tex_zz_xx_1[j] + fl1_fx * tex_xzz_x_0[j] - fl1_fx * tex_xzz_x_1[j] + ta_xzz_xx_1[j];

                    tey_xxzz_xx_0[j] = pa_x[j] * tey_xzz_xx_0[j] - pc_x[j] * tey_xzz_xx_1[j] + 0.5 * fl1_fx * tey_zz_xx_0[j] - 0.5 * fl1_fx * tey_zz_xx_1[j] + fl1_fx * tey_xzz_x_0[j] - fl1_fx * tey_xzz_x_1[j];

                    tez_xxzz_xx_0[j] = pa_x[j] * tez_xzz_xx_0[j] - pc_x[j] * tez_xzz_xx_1[j] + 0.5 * fl1_fx * tez_zz_xx_0[j] - 0.5 * fl1_fx * tez_zz_xx_1[j] + fl1_fx * tez_xzz_x_0[j] - fl1_fx * tez_xzz_x_1[j];

                    tex_xxzz_xy_0[j] = pa_x[j] * tex_xzz_xy_0[j] - pc_x[j] * tex_xzz_xy_1[j] + 0.5 * fl1_fx * tex_zz_xy_0[j] - 0.5 * fl1_fx * tex_zz_xy_1[j] + 0.5 * fl1_fx * tex_xzz_y_0[j] - 0.5 * fl1_fx * tex_xzz_y_1[j] + ta_xzz_xy_1[j];

                    tey_xxzz_xy_0[j] = pa_x[j] * tey_xzz_xy_0[j] - pc_x[j] * tey_xzz_xy_1[j] + 0.5 * fl1_fx * tey_zz_xy_0[j] - 0.5 * fl1_fx * tey_zz_xy_1[j] + 0.5 * fl1_fx * tey_xzz_y_0[j] - 0.5 * fl1_fx * tey_xzz_y_1[j];

                    tez_xxzz_xy_0[j] = pa_x[j] * tez_xzz_xy_0[j] - pc_x[j] * tez_xzz_xy_1[j] + 0.5 * fl1_fx * tez_zz_xy_0[j] - 0.5 * fl1_fx * tez_zz_xy_1[j] + 0.5 * fl1_fx * tez_xzz_y_0[j] - 0.5 * fl1_fx * tez_xzz_y_1[j];

                    tex_xxzz_xz_0[j] = pa_x[j] * tex_xzz_xz_0[j] - pc_x[j] * tex_xzz_xz_1[j] + 0.5 * fl1_fx * tex_zz_xz_0[j] - 0.5 * fl1_fx * tex_zz_xz_1[j] + 0.5 * fl1_fx * tex_xzz_z_0[j] - 0.5 * fl1_fx * tex_xzz_z_1[j] + ta_xzz_xz_1[j];

                    tey_xxzz_xz_0[j] = pa_x[j] * tey_xzz_xz_0[j] - pc_x[j] * tey_xzz_xz_1[j] + 0.5 * fl1_fx * tey_zz_xz_0[j] - 0.5 * fl1_fx * tey_zz_xz_1[j] + 0.5 * fl1_fx * tey_xzz_z_0[j] - 0.5 * fl1_fx * tey_xzz_z_1[j];

                    tez_xxzz_xz_0[j] = pa_x[j] * tez_xzz_xz_0[j] - pc_x[j] * tez_xzz_xz_1[j] + 0.5 * fl1_fx * tez_zz_xz_0[j] - 0.5 * fl1_fx * tez_zz_xz_1[j] + 0.5 * fl1_fx * tez_xzz_z_0[j] - 0.5 * fl1_fx * tez_xzz_z_1[j];

                    tex_xxzz_yy_0[j] = pa_x[j] * tex_xzz_yy_0[j] - pc_x[j] * tex_xzz_yy_1[j] + 0.5 * fl1_fx * tex_zz_yy_0[j] - 0.5 * fl1_fx * tex_zz_yy_1[j] + ta_xzz_yy_1[j];

                    tey_xxzz_yy_0[j] = pa_x[j] * tey_xzz_yy_0[j] - pc_x[j] * tey_xzz_yy_1[j] + 0.5 * fl1_fx * tey_zz_yy_0[j] - 0.5 * fl1_fx * tey_zz_yy_1[j];

                    tez_xxzz_yy_0[j] = pa_x[j] * tez_xzz_yy_0[j] - pc_x[j] * tez_xzz_yy_1[j] + 0.5 * fl1_fx * tez_zz_yy_0[j] - 0.5 * fl1_fx * tez_zz_yy_1[j];

                    tex_xxzz_yz_0[j] = pa_x[j] * tex_xzz_yz_0[j] - pc_x[j] * tex_xzz_yz_1[j] + 0.5 * fl1_fx * tex_zz_yz_0[j] - 0.5 * fl1_fx * tex_zz_yz_1[j] + ta_xzz_yz_1[j];

                    tey_xxzz_yz_0[j] = pa_x[j] * tey_xzz_yz_0[j] - pc_x[j] * tey_xzz_yz_1[j] + 0.5 * fl1_fx * tey_zz_yz_0[j] - 0.5 * fl1_fx * tey_zz_yz_1[j];

                    tez_xxzz_yz_0[j] = pa_x[j] * tez_xzz_yz_0[j] - pc_x[j] * tez_xzz_yz_1[j] + 0.5 * fl1_fx * tez_zz_yz_0[j] - 0.5 * fl1_fx * tez_zz_yz_1[j];

                    tex_xxzz_zz_0[j] = pa_x[j] * tex_xzz_zz_0[j] - pc_x[j] * tex_xzz_zz_1[j] + 0.5 * fl1_fx * tex_zz_zz_0[j] - 0.5 * fl1_fx * tex_zz_zz_1[j] + ta_xzz_zz_1[j];

                    tey_xxzz_zz_0[j] = pa_x[j] * tey_xzz_zz_0[j] - pc_x[j] * tey_xzz_zz_1[j] + 0.5 * fl1_fx * tey_zz_zz_0[j] - 0.5 * fl1_fx * tey_zz_zz_1[j];

                    tez_xxzz_zz_0[j] = pa_x[j] * tez_xzz_zz_0[j] - pc_x[j] * tez_xzz_zz_1[j] + 0.5 * fl1_fx * tez_zz_zz_0[j] - 0.5 * fl1_fx * tez_zz_zz_1[j];

                    tex_xyyy_xx_0[j] = pa_x[j] * tex_yyy_xx_0[j] - pc_x[j] * tex_yyy_xx_1[j] + fl1_fx * tex_yyy_x_0[j] - fl1_fx * tex_yyy_x_1[j] + ta_yyy_xx_1[j];

                    tey_xyyy_xx_0[j] = pa_x[j] * tey_yyy_xx_0[j] - pc_x[j] * tey_yyy_xx_1[j] + fl1_fx * tey_yyy_x_0[j] - fl1_fx * tey_yyy_x_1[j];

                    tez_xyyy_xx_0[j] = pa_x[j] * tez_yyy_xx_0[j] - pc_x[j] * tez_yyy_xx_1[j] + fl1_fx * tez_yyy_x_0[j] - fl1_fx * tez_yyy_x_1[j];

                    tex_xyyy_xy_0[j] = pa_x[j] * tex_yyy_xy_0[j] - pc_x[j] * tex_yyy_xy_1[j] + 0.5 * fl1_fx * tex_yyy_y_0[j] - 0.5 * fl1_fx * tex_yyy_y_1[j] + ta_yyy_xy_1[j];

                    tey_xyyy_xy_0[j] = pa_x[j] * tey_yyy_xy_0[j] - pc_x[j] * tey_yyy_xy_1[j] + 0.5 * fl1_fx * tey_yyy_y_0[j] - 0.5 * fl1_fx * tey_yyy_y_1[j];

                    tez_xyyy_xy_0[j] = pa_x[j] * tez_yyy_xy_0[j] - pc_x[j] * tez_yyy_xy_1[j] + 0.5 * fl1_fx * tez_yyy_y_0[j] - 0.5 * fl1_fx * tez_yyy_y_1[j];

                    tex_xyyy_xz_0[j] = pa_x[j] * tex_yyy_xz_0[j] - pc_x[j] * tex_yyy_xz_1[j] + 0.5 * fl1_fx * tex_yyy_z_0[j] - 0.5 * fl1_fx * tex_yyy_z_1[j] + ta_yyy_xz_1[j];

                    tey_xyyy_xz_0[j] = pa_x[j] * tey_yyy_xz_0[j] - pc_x[j] * tey_yyy_xz_1[j] + 0.5 * fl1_fx * tey_yyy_z_0[j] - 0.5 * fl1_fx * tey_yyy_z_1[j];

                    tez_xyyy_xz_0[j] = pa_x[j] * tez_yyy_xz_0[j] - pc_x[j] * tez_yyy_xz_1[j] + 0.5 * fl1_fx * tez_yyy_z_0[j] - 0.5 * fl1_fx * tez_yyy_z_1[j];

                    tex_xyyy_yy_0[j] = pa_x[j] * tex_yyy_yy_0[j] - pc_x[j] * tex_yyy_yy_1[j] + ta_yyy_yy_1[j];

                    tey_xyyy_yy_0[j] = pa_x[j] * tey_yyy_yy_0[j] - pc_x[j] * tey_yyy_yy_1[j];

                    tez_xyyy_yy_0[j] = pa_x[j] * tez_yyy_yy_0[j] - pc_x[j] * tez_yyy_yy_1[j];

                    tex_xyyy_yz_0[j] = pa_x[j] * tex_yyy_yz_0[j] - pc_x[j] * tex_yyy_yz_1[j] + ta_yyy_yz_1[j];

                    tey_xyyy_yz_0[j] = pa_x[j] * tey_yyy_yz_0[j] - pc_x[j] * tey_yyy_yz_1[j];

                    tez_xyyy_yz_0[j] = pa_x[j] * tez_yyy_yz_0[j] - pc_x[j] * tez_yyy_yz_1[j];

                    tex_xyyy_zz_0[j] = pa_x[j] * tex_yyy_zz_0[j] - pc_x[j] * tex_yyy_zz_1[j] + ta_yyy_zz_1[j];

                    tey_xyyy_zz_0[j] = pa_x[j] * tey_yyy_zz_0[j] - pc_x[j] * tey_yyy_zz_1[j];

                    tez_xyyy_zz_0[j] = pa_x[j] * tez_yyy_zz_0[j] - pc_x[j] * tez_yyy_zz_1[j];

                    tex_xyyz_xx_0[j] = pa_x[j] * tex_yyz_xx_0[j] - pc_x[j] * tex_yyz_xx_1[j] + fl1_fx * tex_yyz_x_0[j] - fl1_fx * tex_yyz_x_1[j] + ta_yyz_xx_1[j];

                    tey_xyyz_xx_0[j] = pa_x[j] * tey_yyz_xx_0[j] - pc_x[j] * tey_yyz_xx_1[j] + fl1_fx * tey_yyz_x_0[j] - fl1_fx * tey_yyz_x_1[j];

                    tez_xyyz_xx_0[j] = pa_x[j] * tez_yyz_xx_0[j] - pc_x[j] * tez_yyz_xx_1[j] + fl1_fx * tez_yyz_x_0[j] - fl1_fx * tez_yyz_x_1[j];

                    tex_xyyz_xy_0[j] = pa_x[j] * tex_yyz_xy_0[j] - pc_x[j] * tex_yyz_xy_1[j] + 0.5 * fl1_fx * tex_yyz_y_0[j] - 0.5 * fl1_fx * tex_yyz_y_1[j] + ta_yyz_xy_1[j];

                    tey_xyyz_xy_0[j] = pa_x[j] * tey_yyz_xy_0[j] - pc_x[j] * tey_yyz_xy_1[j] + 0.5 * fl1_fx * tey_yyz_y_0[j] - 0.5 * fl1_fx * tey_yyz_y_1[j];

                    tez_xyyz_xy_0[j] = pa_x[j] * tez_yyz_xy_0[j] - pc_x[j] * tez_yyz_xy_1[j] + 0.5 * fl1_fx * tez_yyz_y_0[j] - 0.5 * fl1_fx * tez_yyz_y_1[j];

                    tex_xyyz_xz_0[j] = pa_x[j] * tex_yyz_xz_0[j] - pc_x[j] * tex_yyz_xz_1[j] + 0.5 * fl1_fx * tex_yyz_z_0[j] - 0.5 * fl1_fx * tex_yyz_z_1[j] + ta_yyz_xz_1[j];

                    tey_xyyz_xz_0[j] = pa_x[j] * tey_yyz_xz_0[j] - pc_x[j] * tey_yyz_xz_1[j] + 0.5 * fl1_fx * tey_yyz_z_0[j] - 0.5 * fl1_fx * tey_yyz_z_1[j];

                    tez_xyyz_xz_0[j] = pa_x[j] * tez_yyz_xz_0[j] - pc_x[j] * tez_yyz_xz_1[j] + 0.5 * fl1_fx * tez_yyz_z_0[j] - 0.5 * fl1_fx * tez_yyz_z_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForGD_135_180(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (135,180)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {4, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_4_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
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

                auto tex_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 45); 

                auto tey_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 45); 

                auto tez_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 45); 

                auto tex_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 46); 

                auto tey_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 46); 

                auto tez_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 46); 

                auto tex_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 47); 

                auto tey_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 47); 

                auto tez_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 47); 

                auto tex_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 48); 

                auto tey_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 48); 

                auto tez_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 48); 

                auto tex_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 49); 

                auto tey_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 49); 

                auto tez_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 49); 

                auto tex_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 50); 

                auto tey_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 50); 

                auto tez_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 50); 

                auto tex_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 51); 

                auto tey_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 51); 

                auto tez_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 51); 

                auto tex_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 52); 

                auto tey_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 52); 

                auto tez_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 52); 

                auto tex_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 53); 

                auto tey_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 53); 

                auto tez_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 53); 

                auto tex_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 54); 

                auto tey_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 54); 

                auto tez_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 54); 

                auto tex_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 55); 

                auto tey_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 55); 

                auto tez_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 55); 

                auto tex_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 56); 

                auto tey_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 56); 

                auto tez_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 56); 

                auto tex_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 57); 

                auto tey_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 57); 

                auto tez_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 57); 

                auto tex_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 58); 

                auto tey_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 58); 

                auto tez_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 58); 

                auto tex_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 59); 

                auto tey_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 59); 

                auto tez_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 59); 

                auto tex_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 45); 

                auto tey_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 45); 

                auto tez_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 45); 

                auto tex_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 46); 

                auto tey_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 46); 

                auto tez_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 46); 

                auto tex_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 47); 

                auto tey_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 47); 

                auto tez_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 47); 

                auto tex_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 48); 

                auto tey_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 48); 

                auto tez_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 48); 

                auto tex_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 49); 

                auto tey_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 49); 

                auto tez_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 49); 

                auto tex_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 50); 

                auto tey_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 50); 

                auto tez_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 50); 

                auto tex_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 51); 

                auto tey_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 51); 

                auto tez_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 51); 

                auto tex_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 52); 

                auto tey_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 52); 

                auto tez_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 52); 

                auto tex_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 53); 

                auto tey_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 53); 

                auto tez_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 53); 

                auto tex_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 54); 

                auto tey_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 54); 

                auto tez_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 54); 

                auto tex_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 55); 

                auto tey_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 55); 

                auto tez_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 55); 

                auto tex_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 56); 

                auto tey_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 56); 

                auto tez_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 56); 

                auto tex_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 57); 

                auto tey_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 57); 

                auto tez_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 57); 

                auto tex_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 58); 

                auto tey_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 58); 

                auto tez_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 58); 

                auto tex_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 59); 

                auto tey_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 59); 

                auto tez_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 59); 

                auto tex_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 24); 

                auto tey_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 25); 

                auto tey_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 25); 

                auto tez_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 25); 

                auto tex_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 26); 

                auto tey_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 26); 

                auto tez_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 26); 

                auto tex_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 27); 

                auto tey_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 27); 

                auto tez_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 27); 

                auto tex_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 28); 

                auto tey_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 28); 

                auto tez_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 28); 

                auto tex_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 29); 

                auto tey_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 29); 

                auto tez_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 29); 

                auto tex_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 24); 

                auto tey_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 24); 

                auto tex_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 25); 

                auto tey_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 25); 

                auto tez_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 25); 

                auto tex_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 26); 

                auto tey_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 26); 

                auto tez_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 26); 

                auto tex_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 27); 

                auto tey_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 27); 

                auto tez_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 27); 

                auto tex_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 28); 

                auto tey_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 28); 

                auto tez_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 28); 

                auto tex_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 29); 

                auto tey_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 29); 

                auto tez_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 29); 

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

                auto tex_xyyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 45); 

                auto tey_xyyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 45); 

                auto tez_xyyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 45); 

                auto tex_xyyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 46); 

                auto tey_xyyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 46); 

                auto tez_xyyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 46); 

                auto tex_xyyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 47); 

                auto tey_xyyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 47); 

                auto tez_xyyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 47); 

                auto tex_xyzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 48); 

                auto tey_xyzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 48); 

                auto tez_xyzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 48); 

                auto tex_xyzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 49); 

                auto tey_xyzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 49); 

                auto tez_xyzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 49); 

                auto tex_xyzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 50); 

                auto tey_xyzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 50); 

                auto tez_xyzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 50); 

                auto tex_xyzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 51); 

                auto tey_xyzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 51); 

                auto tez_xyzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 51); 

                auto tex_xyzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 52); 

                auto tey_xyzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 52); 

                auto tez_xyzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 52); 

                auto tex_xyzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 53); 

                auto tey_xyzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 53); 

                auto tez_xyzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 53); 

                auto tex_xzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 54); 

                auto tey_xzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 54); 

                auto tez_xzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 54); 

                auto tex_xzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 55); 

                auto tey_xzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 55); 

                auto tez_xzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 55); 

                auto tex_xzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 56); 

                auto tey_xzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 56); 

                auto tez_xzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 56); 

                auto tex_xzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 57); 

                auto tey_xzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 57); 

                auto tez_xzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 57); 

                auto tex_xzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 58); 

                auto tey_xzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 58); 

                auto tez_xzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 58); 

                auto tex_xzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 59); 

                auto tey_xzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 59); 

                auto tez_xzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 59); 

                // Batch of Integrals (135,180)

                #pragma omp simd aligned(fx, pa_x, pc_x, ta_yyz_yy_1, ta_yyz_yz_1, ta_yyz_zz_1, ta_yzz_xx_1, \
                                         ta_yzz_xy_1, ta_yzz_xz_1, ta_yzz_yy_1, ta_yzz_yz_1, ta_yzz_zz_1, ta_zzz_xx_1, \
                                         ta_zzz_xy_1, ta_zzz_xz_1, ta_zzz_yy_1, ta_zzz_yz_1, ta_zzz_zz_1, tex_xyyz_yy_0, \
                                         tex_xyyz_yz_0, tex_xyyz_zz_0, tex_xyzz_xx_0, tex_xyzz_xy_0, tex_xyzz_xz_0, \
                                         tex_xyzz_yy_0, tex_xyzz_yz_0, tex_xyzz_zz_0, tex_xzzz_xx_0, tex_xzzz_xy_0, \
                                         tex_xzzz_xz_0, tex_xzzz_yy_0, tex_xzzz_yz_0, tex_xzzz_zz_0, tex_yyz_yy_0, \
                                         tex_yyz_yy_1, tex_yyz_yz_0, tex_yyz_yz_1, tex_yyz_zz_0, tex_yyz_zz_1, tex_yzz_x_0, \
                                         tex_yzz_x_1, tex_yzz_xx_0, tex_yzz_xx_1, tex_yzz_xy_0, tex_yzz_xy_1, tex_yzz_xz_0, \
                                         tex_yzz_xz_1, tex_yzz_y_0, tex_yzz_y_1, tex_yzz_yy_0, tex_yzz_yy_1, tex_yzz_yz_0, \
                                         tex_yzz_yz_1, tex_yzz_z_0, tex_yzz_z_1, tex_yzz_zz_0, tex_yzz_zz_1, tex_zzz_x_0, \
                                         tex_zzz_x_1, tex_zzz_xx_0, tex_zzz_xx_1, tex_zzz_xy_0, tex_zzz_xy_1, tex_zzz_xz_0, \
                                         tex_zzz_xz_1, tex_zzz_y_0, tex_zzz_y_1, tex_zzz_yy_0, tex_zzz_yy_1, tex_zzz_yz_0, \
                                         tex_zzz_yz_1, tex_zzz_z_0, tex_zzz_z_1, tex_zzz_zz_0, tex_zzz_zz_1, tey_xyyz_yy_0, \
                                         tey_xyyz_yz_0, tey_xyyz_zz_0, tey_xyzz_xx_0, tey_xyzz_xy_0, tey_xyzz_xz_0, \
                                         tey_xyzz_yy_0, tey_xyzz_yz_0, tey_xyzz_zz_0, tey_xzzz_xx_0, tey_xzzz_xy_0, \
                                         tey_xzzz_xz_0, tey_xzzz_yy_0, tey_xzzz_yz_0, tey_xzzz_zz_0, tey_yyz_yy_0, \
                                         tey_yyz_yy_1, tey_yyz_yz_0, tey_yyz_yz_1, tey_yyz_zz_0, tey_yyz_zz_1, tey_yzz_x_0, \
                                         tey_yzz_x_1, tey_yzz_xx_0, tey_yzz_xx_1, tey_yzz_xy_0, tey_yzz_xy_1, tey_yzz_xz_0, \
                                         tey_yzz_xz_1, tey_yzz_y_0, tey_yzz_y_1, tey_yzz_yy_0, tey_yzz_yy_1, tey_yzz_yz_0, \
                                         tey_yzz_yz_1, tey_yzz_z_0, tey_yzz_z_1, tey_yzz_zz_0, tey_yzz_zz_1, tey_zzz_x_0, \
                                         tey_zzz_x_1, tey_zzz_xx_0, tey_zzz_xx_1, tey_zzz_xy_0, tey_zzz_xy_1, tey_zzz_xz_0, \
                                         tey_zzz_xz_1, tey_zzz_y_0, tey_zzz_y_1, tey_zzz_yy_0, tey_zzz_yy_1, tey_zzz_yz_0, \
                                         tey_zzz_yz_1, tey_zzz_z_0, tey_zzz_z_1, tey_zzz_zz_0, tey_zzz_zz_1, tez_xyyz_yy_0, \
                                         tez_xyyz_yz_0, tez_xyyz_zz_0, tez_xyzz_xx_0, tez_xyzz_xy_0, tez_xyzz_xz_0, \
                                         tez_xyzz_yy_0, tez_xyzz_yz_0, tez_xyzz_zz_0, tez_xzzz_xx_0, tez_xzzz_xy_0, \
                                         tez_xzzz_xz_0, tez_xzzz_yy_0, tez_xzzz_yz_0, tez_xzzz_zz_0, tez_yyz_yy_0, \
                                         tez_yyz_yy_1, tez_yyz_yz_0, tez_yyz_yz_1, tez_yyz_zz_0, tez_yyz_zz_1, tez_yzz_x_0, \
                                         tez_yzz_x_1, tez_yzz_xx_0, tez_yzz_xx_1, tez_yzz_xy_0, tez_yzz_xy_1, tez_yzz_xz_0, \
                                         tez_yzz_xz_1, tez_yzz_y_0, tez_yzz_y_1, tez_yzz_yy_0, tez_yzz_yy_1, tez_yzz_yz_0, \
                                         tez_yzz_yz_1, tez_yzz_z_0, tez_yzz_z_1, tez_yzz_zz_0, tez_yzz_zz_1, tez_zzz_x_0, \
                                         tez_zzz_x_1, tez_zzz_xx_0, tez_zzz_xx_1, tez_zzz_xy_0, tez_zzz_xy_1, tez_zzz_xz_0, \
                                         tez_zzz_xz_1, tez_zzz_y_0, tez_zzz_y_1, tez_zzz_yy_0, tez_zzz_yy_1, tez_zzz_yz_0, \
                                         tez_zzz_yz_1, tez_zzz_z_0, tez_zzz_z_1, tez_zzz_zz_0, tez_zzz_zz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_xyyz_yy_0[j] = pa_x[j] * tex_yyz_yy_0[j] - pc_x[j] * tex_yyz_yy_1[j] + ta_yyz_yy_1[j];

                    tey_xyyz_yy_0[j] = pa_x[j] * tey_yyz_yy_0[j] - pc_x[j] * tey_yyz_yy_1[j];

                    tez_xyyz_yy_0[j] = pa_x[j] * tez_yyz_yy_0[j] - pc_x[j] * tez_yyz_yy_1[j];

                    tex_xyyz_yz_0[j] = pa_x[j] * tex_yyz_yz_0[j] - pc_x[j] * tex_yyz_yz_1[j] + ta_yyz_yz_1[j];

                    tey_xyyz_yz_0[j] = pa_x[j] * tey_yyz_yz_0[j] - pc_x[j] * tey_yyz_yz_1[j];

                    tez_xyyz_yz_0[j] = pa_x[j] * tez_yyz_yz_0[j] - pc_x[j] * tez_yyz_yz_1[j];

                    tex_xyyz_zz_0[j] = pa_x[j] * tex_yyz_zz_0[j] - pc_x[j] * tex_yyz_zz_1[j] + ta_yyz_zz_1[j];

                    tey_xyyz_zz_0[j] = pa_x[j] * tey_yyz_zz_0[j] - pc_x[j] * tey_yyz_zz_1[j];

                    tez_xyyz_zz_0[j] = pa_x[j] * tez_yyz_zz_0[j] - pc_x[j] * tez_yyz_zz_1[j];

                    tex_xyzz_xx_0[j] = pa_x[j] * tex_yzz_xx_0[j] - pc_x[j] * tex_yzz_xx_1[j] + fl1_fx * tex_yzz_x_0[j] - fl1_fx * tex_yzz_x_1[j] + ta_yzz_xx_1[j];

                    tey_xyzz_xx_0[j] = pa_x[j] * tey_yzz_xx_0[j] - pc_x[j] * tey_yzz_xx_1[j] + fl1_fx * tey_yzz_x_0[j] - fl1_fx * tey_yzz_x_1[j];

                    tez_xyzz_xx_0[j] = pa_x[j] * tez_yzz_xx_0[j] - pc_x[j] * tez_yzz_xx_1[j] + fl1_fx * tez_yzz_x_0[j] - fl1_fx * tez_yzz_x_1[j];

                    tex_xyzz_xy_0[j] = pa_x[j] * tex_yzz_xy_0[j] - pc_x[j] * tex_yzz_xy_1[j] + 0.5 * fl1_fx * tex_yzz_y_0[j] - 0.5 * fl1_fx * tex_yzz_y_1[j] + ta_yzz_xy_1[j];

                    tey_xyzz_xy_0[j] = pa_x[j] * tey_yzz_xy_0[j] - pc_x[j] * tey_yzz_xy_1[j] + 0.5 * fl1_fx * tey_yzz_y_0[j] - 0.5 * fl1_fx * tey_yzz_y_1[j];

                    tez_xyzz_xy_0[j] = pa_x[j] * tez_yzz_xy_0[j] - pc_x[j] * tez_yzz_xy_1[j] + 0.5 * fl1_fx * tez_yzz_y_0[j] - 0.5 * fl1_fx * tez_yzz_y_1[j];

                    tex_xyzz_xz_0[j] = pa_x[j] * tex_yzz_xz_0[j] - pc_x[j] * tex_yzz_xz_1[j] + 0.5 * fl1_fx * tex_yzz_z_0[j] - 0.5 * fl1_fx * tex_yzz_z_1[j] + ta_yzz_xz_1[j];

                    tey_xyzz_xz_0[j] = pa_x[j] * tey_yzz_xz_0[j] - pc_x[j] * tey_yzz_xz_1[j] + 0.5 * fl1_fx * tey_yzz_z_0[j] - 0.5 * fl1_fx * tey_yzz_z_1[j];

                    tez_xyzz_xz_0[j] = pa_x[j] * tez_yzz_xz_0[j] - pc_x[j] * tez_yzz_xz_1[j] + 0.5 * fl1_fx * tez_yzz_z_0[j] - 0.5 * fl1_fx * tez_yzz_z_1[j];

                    tex_xyzz_yy_0[j] = pa_x[j] * tex_yzz_yy_0[j] - pc_x[j] * tex_yzz_yy_1[j] + ta_yzz_yy_1[j];

                    tey_xyzz_yy_0[j] = pa_x[j] * tey_yzz_yy_0[j] - pc_x[j] * tey_yzz_yy_1[j];

                    tez_xyzz_yy_0[j] = pa_x[j] * tez_yzz_yy_0[j] - pc_x[j] * tez_yzz_yy_1[j];

                    tex_xyzz_yz_0[j] = pa_x[j] * tex_yzz_yz_0[j] - pc_x[j] * tex_yzz_yz_1[j] + ta_yzz_yz_1[j];

                    tey_xyzz_yz_0[j] = pa_x[j] * tey_yzz_yz_0[j] - pc_x[j] * tey_yzz_yz_1[j];

                    tez_xyzz_yz_0[j] = pa_x[j] * tez_yzz_yz_0[j] - pc_x[j] * tez_yzz_yz_1[j];

                    tex_xyzz_zz_0[j] = pa_x[j] * tex_yzz_zz_0[j] - pc_x[j] * tex_yzz_zz_1[j] + ta_yzz_zz_1[j];

                    tey_xyzz_zz_0[j] = pa_x[j] * tey_yzz_zz_0[j] - pc_x[j] * tey_yzz_zz_1[j];

                    tez_xyzz_zz_0[j] = pa_x[j] * tez_yzz_zz_0[j] - pc_x[j] * tez_yzz_zz_1[j];

                    tex_xzzz_xx_0[j] = pa_x[j] * tex_zzz_xx_0[j] - pc_x[j] * tex_zzz_xx_1[j] + fl1_fx * tex_zzz_x_0[j] - fl1_fx * tex_zzz_x_1[j] + ta_zzz_xx_1[j];

                    tey_xzzz_xx_0[j] = pa_x[j] * tey_zzz_xx_0[j] - pc_x[j] * tey_zzz_xx_1[j] + fl1_fx * tey_zzz_x_0[j] - fl1_fx * tey_zzz_x_1[j];

                    tez_xzzz_xx_0[j] = pa_x[j] * tez_zzz_xx_0[j] - pc_x[j] * tez_zzz_xx_1[j] + fl1_fx * tez_zzz_x_0[j] - fl1_fx * tez_zzz_x_1[j];

                    tex_xzzz_xy_0[j] = pa_x[j] * tex_zzz_xy_0[j] - pc_x[j] * tex_zzz_xy_1[j] + 0.5 * fl1_fx * tex_zzz_y_0[j] - 0.5 * fl1_fx * tex_zzz_y_1[j] + ta_zzz_xy_1[j];

                    tey_xzzz_xy_0[j] = pa_x[j] * tey_zzz_xy_0[j] - pc_x[j] * tey_zzz_xy_1[j] + 0.5 * fl1_fx * tey_zzz_y_0[j] - 0.5 * fl1_fx * tey_zzz_y_1[j];

                    tez_xzzz_xy_0[j] = pa_x[j] * tez_zzz_xy_0[j] - pc_x[j] * tez_zzz_xy_1[j] + 0.5 * fl1_fx * tez_zzz_y_0[j] - 0.5 * fl1_fx * tez_zzz_y_1[j];

                    tex_xzzz_xz_0[j] = pa_x[j] * tex_zzz_xz_0[j] - pc_x[j] * tex_zzz_xz_1[j] + 0.5 * fl1_fx * tex_zzz_z_0[j] - 0.5 * fl1_fx * tex_zzz_z_1[j] + ta_zzz_xz_1[j];

                    tey_xzzz_xz_0[j] = pa_x[j] * tey_zzz_xz_0[j] - pc_x[j] * tey_zzz_xz_1[j] + 0.5 * fl1_fx * tey_zzz_z_0[j] - 0.5 * fl1_fx * tey_zzz_z_1[j];

                    tez_xzzz_xz_0[j] = pa_x[j] * tez_zzz_xz_0[j] - pc_x[j] * tez_zzz_xz_1[j] + 0.5 * fl1_fx * tez_zzz_z_0[j] - 0.5 * fl1_fx * tez_zzz_z_1[j];

                    tex_xzzz_yy_0[j] = pa_x[j] * tex_zzz_yy_0[j] - pc_x[j] * tex_zzz_yy_1[j] + ta_zzz_yy_1[j];

                    tey_xzzz_yy_0[j] = pa_x[j] * tey_zzz_yy_0[j] - pc_x[j] * tey_zzz_yy_1[j];

                    tez_xzzz_yy_0[j] = pa_x[j] * tez_zzz_yy_0[j] - pc_x[j] * tez_zzz_yy_1[j];

                    tex_xzzz_yz_0[j] = pa_x[j] * tex_zzz_yz_0[j] - pc_x[j] * tex_zzz_yz_1[j] + ta_zzz_yz_1[j];

                    tey_xzzz_yz_0[j] = pa_x[j] * tey_zzz_yz_0[j] - pc_x[j] * tey_zzz_yz_1[j];

                    tez_xzzz_yz_0[j] = pa_x[j] * tez_zzz_yz_0[j] - pc_x[j] * tez_zzz_yz_1[j];

                    tex_xzzz_zz_0[j] = pa_x[j] * tex_zzz_zz_0[j] - pc_x[j] * tex_zzz_zz_1[j] + ta_zzz_zz_1[j];

                    tey_xzzz_zz_0[j] = pa_x[j] * tey_zzz_zz_0[j] - pc_x[j] * tey_zzz_zz_1[j];

                    tez_xzzz_zz_0[j] = pa_x[j] * tez_zzz_zz_0[j] - pc_x[j] * tez_zzz_zz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForGD_180_225(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (180,225)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {4, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_4_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_a_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
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

                auto tex_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 36); 

                auto tey_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 36); 

                auto tez_yyy_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 36); 

                auto tex_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 37); 

                auto tey_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 37); 

                auto tez_yyy_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 37); 

                auto tex_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 38); 

                auto tey_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 38); 

                auto tez_yyy_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 38); 

                auto tex_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 39); 

                auto tey_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 39); 

                auto tez_yyy_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 39); 

                auto tex_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 40); 

                auto tey_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 40); 

                auto tez_yyy_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 40); 

                auto tex_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 41); 

                auto tey_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 41); 

                auto tez_yyy_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 41); 

                auto tex_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 42); 

                auto tey_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 42); 

                auto tez_yyz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 42); 

                auto tex_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 43); 

                auto tey_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 43); 

                auto tez_yyz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 43); 

                auto tex_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 44); 

                auto tey_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 44); 

                auto tez_yyz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 44); 

                auto tex_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 45); 

                auto tey_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 45); 

                auto tez_yyz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 45); 

                auto tex_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 46); 

                auto tey_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 46); 

                auto tez_yyz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 46); 

                auto tex_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 47); 

                auto tey_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 47); 

                auto tez_yyz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 47); 

                auto tex_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 48); 

                auto tey_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 48); 

                auto tez_yzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 48); 

                auto tex_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 49); 

                auto tey_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 49); 

                auto tez_yzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 49); 

                auto tex_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 50); 

                auto tey_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 50); 

                auto tez_yzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 50); 

                auto tex_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 36); 

                auto tey_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 36); 

                auto tez_yyy_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 36); 

                auto tex_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 37); 

                auto tey_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 37); 

                auto tez_yyy_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 37); 

                auto tex_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 38); 

                auto tey_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 38); 

                auto tez_yyy_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 38); 

                auto tex_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 39); 

                auto tey_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 39); 

                auto tez_yyy_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 39); 

                auto tex_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 40); 

                auto tey_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 40); 

                auto tez_yyy_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 40); 

                auto tex_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 41); 

                auto tey_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 41); 

                auto tez_yyy_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 41); 

                auto tex_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 42); 

                auto tey_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 42); 

                auto tez_yyz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 42); 

                auto tex_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 43); 

                auto tey_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 43); 

                auto tez_yyz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 43); 

                auto tex_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 44); 

                auto tey_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 44); 

                auto tez_yyz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 44); 

                auto tex_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 45); 

                auto tey_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 45); 

                auto tez_yyz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 45); 

                auto tex_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 46); 

                auto tey_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 46); 

                auto tez_yyz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 46); 

                auto tex_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 47); 

                auto tey_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 47); 

                auto tez_yyz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 47); 

                auto tex_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 48); 

                auto tey_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 48); 

                auto tez_yzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 48); 

                auto tex_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 49); 

                auto tey_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 49); 

                auto tez_yzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 49); 

                auto tex_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 50); 

                auto tey_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 50); 

                auto tez_yzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 50); 

                auto tex_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 18); 

                auto tey_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 19); 

                auto tey_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 20); 

                auto tey_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 21); 

                auto tey_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 22); 

                auto tey_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 23); 

                auto tey_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 24); 

                auto tey_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 25); 

                auto tey_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 26); 

                auto tey_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 26); 

                auto tex_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 27); 

                auto tey_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 28); 

                auto tey_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 29); 

                auto tey_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 29); 

                auto tex_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 30); 

                auto tey_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 31); 

                auto tey_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 32); 

                auto tey_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 32); 

                auto tex_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 18); 

                auto tey_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 18); 

                auto tez_yy_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 18); 

                auto tex_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 19); 

                auto tey_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 19); 

                auto tez_yy_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 19); 

                auto tex_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 20); 

                auto tey_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 20); 

                auto tez_yy_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 20); 

                auto tex_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 21); 

                auto tey_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 21); 

                auto tez_yy_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 21); 

                auto tex_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 22); 

                auto tey_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 22); 

                auto tez_yy_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 22); 

                auto tex_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 23); 

                auto tey_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 23); 

                auto tez_yy_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 23); 

                auto tex_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 24); 

                auto tey_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 24); 

                auto tez_yz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 24); 

                auto tex_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 25); 

                auto tey_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 25); 

                auto tez_yz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 25); 

                auto tex_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 26); 

                auto tey_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 26); 

                auto tez_yz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 26); 

                auto tex_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 27); 

                auto tey_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 27); 

                auto tez_yz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 27); 

                auto tex_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 28); 

                auto tey_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 28); 

                auto tez_yz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 28); 

                auto tex_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 29); 

                auto tey_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 29); 

                auto tez_yz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 29); 

                auto tex_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 30); 

                auto tey_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 31); 

                auto tey_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 32); 

                auto tey_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 32); 

                auto tex_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 18); 

                auto tey_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 18); 

                auto tez_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 18); 

                auto tex_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 19); 

                auto tey_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 19); 

                auto tez_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 19); 

                auto tex_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 20); 

                auto tey_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 20); 

                auto tez_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 20); 

                auto tex_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 21); 

                auto tey_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 21); 

                auto tez_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 21); 

                auto tex_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 22); 

                auto tey_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 22); 

                auto tez_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 22); 

                auto tex_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 23); 

                auto tey_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 23); 

                auto tez_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 23); 

                auto tex_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 24); 

                auto tey_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 24); 

                auto tez_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 24); 

                auto tex_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 18); 

                auto tey_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 18); 

                auto tez_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 18); 

                auto tex_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 19); 

                auto tey_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 19); 

                auto tez_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 19); 

                auto tex_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 20); 

                auto tey_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 20); 

                auto tez_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 20); 

                auto tex_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 21); 

                auto tey_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 21); 

                auto tez_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 21); 

                auto tex_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 22); 

                auto tey_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 22); 

                auto tez_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 22); 

                auto tex_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 23); 

                auto tey_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 23); 

                auto tez_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 23); 

                auto tex_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 24); 

                auto tey_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 24); 

                auto tez_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 24); 

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

                // set up pointers to integrals

                auto tex_yyyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 60); 

                auto tey_yyyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 60); 

                auto tez_yyyy_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 60); 

                auto tex_yyyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 61); 

                auto tey_yyyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 61); 

                auto tez_yyyy_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 61); 

                auto tex_yyyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 62); 

                auto tey_yyyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 62); 

                auto tez_yyyy_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 62); 

                auto tex_yyyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 63); 

                auto tey_yyyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 63); 

                auto tez_yyyy_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 63); 

                auto tex_yyyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 64); 

                auto tey_yyyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 64); 

                auto tez_yyyy_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 64); 

                auto tex_yyyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 65); 

                auto tey_yyyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 65); 

                auto tez_yyyy_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 65); 

                auto tex_yyyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 66); 

                auto tey_yyyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 66); 

                auto tez_yyyz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 66); 

                auto tex_yyyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 67); 

                auto tey_yyyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 67); 

                auto tez_yyyz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 67); 

                auto tex_yyyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 68); 

                auto tey_yyyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 68); 

                auto tez_yyyz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 68); 

                auto tex_yyyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 69); 

                auto tey_yyyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 69); 

                auto tez_yyyz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 69); 

                auto tex_yyyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 70); 

                auto tey_yyyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 70); 

                auto tez_yyyz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 70); 

                auto tex_yyyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 71); 

                auto tey_yyyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 71); 

                auto tez_yyyz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 71); 

                auto tex_yyzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 72); 

                auto tey_yyzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 72); 

                auto tez_yyzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 72); 

                auto tex_yyzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 73); 

                auto tey_yyzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 73); 

                auto tez_yyzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 73); 

                auto tex_yyzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 74); 

                auto tey_yyzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 74); 

                auto tez_yyzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 74); 

                // Batch of Integrals (180,225)

                #pragma omp simd aligned(fx, pa_y, pc_y, ta_yyy_xx_1, ta_yyy_xy_1, ta_yyy_xz_1, ta_yyy_yy_1, \
                                         ta_yyy_yz_1, ta_yyy_zz_1, ta_yyz_xx_1, ta_yyz_xy_1, ta_yyz_xz_1, ta_yyz_yy_1, \
                                         ta_yyz_yz_1, ta_yyz_zz_1, ta_yzz_xx_1, ta_yzz_xy_1, ta_yzz_xz_1, tex_yy_xx_0, \
                                         tex_yy_xx_1, tex_yy_xy_0, tex_yy_xy_1, tex_yy_xz_0, tex_yy_xz_1, tex_yy_yy_0, \
                                         tex_yy_yy_1, tex_yy_yz_0, tex_yy_yz_1, tex_yy_zz_0, tex_yy_zz_1, tex_yyy_x_0, \
                                         tex_yyy_x_1, tex_yyy_xx_0, tex_yyy_xx_1, tex_yyy_xy_0, tex_yyy_xy_1, tex_yyy_xz_0, \
                                         tex_yyy_xz_1, tex_yyy_y_0, tex_yyy_y_1, tex_yyy_yy_0, tex_yyy_yy_1, tex_yyy_yz_0, \
                                         tex_yyy_yz_1, tex_yyy_z_0, tex_yyy_z_1, tex_yyy_zz_0, tex_yyy_zz_1, tex_yyyy_xx_0, \
                                         tex_yyyy_xy_0, tex_yyyy_xz_0, tex_yyyy_yy_0, tex_yyyy_yz_0, tex_yyyy_zz_0, \
                                         tex_yyyz_xx_0, tex_yyyz_xy_0, tex_yyyz_xz_0, tex_yyyz_yy_0, tex_yyyz_yz_0, \
                                         tex_yyyz_zz_0, tex_yyz_x_0, tex_yyz_x_1, tex_yyz_xx_0, tex_yyz_xx_1, tex_yyz_xy_0, \
                                         tex_yyz_xy_1, tex_yyz_xz_0, tex_yyz_xz_1, tex_yyz_y_0, tex_yyz_y_1, tex_yyz_yy_0, \
                                         tex_yyz_yy_1, tex_yyz_yz_0, tex_yyz_yz_1, tex_yyz_z_0, tex_yyz_z_1, tex_yyz_zz_0, \
                                         tex_yyz_zz_1, tex_yyzz_xx_0, tex_yyzz_xy_0, tex_yyzz_xz_0, tex_yz_xx_0, \
                                         tex_yz_xx_1, tex_yz_xy_0, tex_yz_xy_1, tex_yz_xz_0, tex_yz_xz_1, tex_yz_yy_0, \
                                         tex_yz_yy_1, tex_yz_yz_0, tex_yz_yz_1, tex_yz_zz_0, tex_yz_zz_1, tex_yzz_x_0, \
                                         tex_yzz_x_1, tex_yzz_xx_0, tex_yzz_xx_1, tex_yzz_xy_0, tex_yzz_xy_1, tex_yzz_xz_0, \
                                         tex_yzz_xz_1, tex_zz_xx_0, tex_zz_xx_1, tex_zz_xy_0, tex_zz_xy_1, tex_zz_xz_0, \
                                         tex_zz_xz_1, tey_yy_xx_0, tey_yy_xx_1, tey_yy_xy_0, tey_yy_xy_1, tey_yy_xz_0, \
                                         tey_yy_xz_1, tey_yy_yy_0, tey_yy_yy_1, tey_yy_yz_0, tey_yy_yz_1, tey_yy_zz_0, \
                                         tey_yy_zz_1, tey_yyy_x_0, tey_yyy_x_1, tey_yyy_xx_0, tey_yyy_xx_1, tey_yyy_xy_0, \
                                         tey_yyy_xy_1, tey_yyy_xz_0, tey_yyy_xz_1, tey_yyy_y_0, tey_yyy_y_1, tey_yyy_yy_0, \
                                         tey_yyy_yy_1, tey_yyy_yz_0, tey_yyy_yz_1, tey_yyy_z_0, tey_yyy_z_1, tey_yyy_zz_0, \
                                         tey_yyy_zz_1, tey_yyyy_xx_0, tey_yyyy_xy_0, tey_yyyy_xz_0, tey_yyyy_yy_0, \
                                         tey_yyyy_yz_0, tey_yyyy_zz_0, tey_yyyz_xx_0, tey_yyyz_xy_0, tey_yyyz_xz_0, \
                                         tey_yyyz_yy_0, tey_yyyz_yz_0, tey_yyyz_zz_0, tey_yyz_x_0, tey_yyz_x_1, tey_yyz_xx_0, \
                                         tey_yyz_xx_1, tey_yyz_xy_0, tey_yyz_xy_1, tey_yyz_xz_0, tey_yyz_xz_1, tey_yyz_y_0, \
                                         tey_yyz_y_1, tey_yyz_yy_0, tey_yyz_yy_1, tey_yyz_yz_0, tey_yyz_yz_1, tey_yyz_z_0, \
                                         tey_yyz_z_1, tey_yyz_zz_0, tey_yyz_zz_1, tey_yyzz_xx_0, tey_yyzz_xy_0, \
                                         tey_yyzz_xz_0, tey_yz_xx_0, tey_yz_xx_1, tey_yz_xy_0, tey_yz_xy_1, tey_yz_xz_0, \
                                         tey_yz_xz_1, tey_yz_yy_0, tey_yz_yy_1, tey_yz_yz_0, tey_yz_yz_1, tey_yz_zz_0, \
                                         tey_yz_zz_1, tey_yzz_x_0, tey_yzz_x_1, tey_yzz_xx_0, tey_yzz_xx_1, tey_yzz_xy_0, \
                                         tey_yzz_xy_1, tey_yzz_xz_0, tey_yzz_xz_1, tey_zz_xx_0, tey_zz_xx_1, tey_zz_xy_0, \
                                         tey_zz_xy_1, tey_zz_xz_0, tey_zz_xz_1, tez_yy_xx_0, tez_yy_xx_1, tez_yy_xy_0, \
                                         tez_yy_xy_1, tez_yy_xz_0, tez_yy_xz_1, tez_yy_yy_0, tez_yy_yy_1, tez_yy_yz_0, \
                                         tez_yy_yz_1, tez_yy_zz_0, tez_yy_zz_1, tez_yyy_x_0, tez_yyy_x_1, tez_yyy_xx_0, \
                                         tez_yyy_xx_1, tez_yyy_xy_0, tez_yyy_xy_1, tez_yyy_xz_0, tez_yyy_xz_1, tez_yyy_y_0, \
                                         tez_yyy_y_1, tez_yyy_yy_0, tez_yyy_yy_1, tez_yyy_yz_0, tez_yyy_yz_1, tez_yyy_z_0, \
                                         tez_yyy_z_1, tez_yyy_zz_0, tez_yyy_zz_1, tez_yyyy_xx_0, tez_yyyy_xy_0, \
                                         tez_yyyy_xz_0, tez_yyyy_yy_0, tez_yyyy_yz_0, tez_yyyy_zz_0, tez_yyyz_xx_0, \
                                         tez_yyyz_xy_0, tez_yyyz_xz_0, tez_yyyz_yy_0, tez_yyyz_yz_0, tez_yyyz_zz_0, \
                                         tez_yyz_x_0, tez_yyz_x_1, tez_yyz_xx_0, tez_yyz_xx_1, tez_yyz_xy_0, tez_yyz_xy_1, \
                                         tez_yyz_xz_0, tez_yyz_xz_1, tez_yyz_y_0, tez_yyz_y_1, tez_yyz_yy_0, tez_yyz_yy_1, \
                                         tez_yyz_yz_0, tez_yyz_yz_1, tez_yyz_z_0, tez_yyz_z_1, tez_yyz_zz_0, tez_yyz_zz_1, \
                                         tez_yyzz_xx_0, tez_yyzz_xy_0, tez_yyzz_xz_0, tez_yz_xx_0, tez_yz_xx_1, tez_yz_xy_0, \
                                         tez_yz_xy_1, tez_yz_xz_0, tez_yz_xz_1, tez_yz_yy_0, tez_yz_yy_1, tez_yz_yz_0, \
                                         tez_yz_yz_1, tez_yz_zz_0, tez_yz_zz_1, tez_yzz_x_0, tez_yzz_x_1, tez_yzz_xx_0, \
                                         tez_yzz_xx_1, tez_yzz_xy_0, tez_yzz_xy_1, tez_yzz_xz_0, tez_yzz_xz_1, tez_zz_xx_0, \
                                         tez_zz_xx_1, tez_zz_xy_0, tez_zz_xy_1, tez_zz_xz_0, tez_zz_xz_1: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yyyy_xx_0[j] = pa_y[j] * tex_yyy_xx_0[j] - pc_y[j] * tex_yyy_xx_1[j] + 1.5 * fl1_fx * tex_yy_xx_0[j] - 1.5 * fl1_fx * tex_yy_xx_1[j];

                    tey_yyyy_xx_0[j] = pa_y[j] * tey_yyy_xx_0[j] - pc_y[j] * tey_yyy_xx_1[j] + 1.5 * fl1_fx * tey_yy_xx_0[j] - 1.5 * fl1_fx * tey_yy_xx_1[j] + ta_yyy_xx_1[j];

                    tez_yyyy_xx_0[j] = pa_y[j] * tez_yyy_xx_0[j] - pc_y[j] * tez_yyy_xx_1[j] + 1.5 * fl1_fx * tez_yy_xx_0[j] - 1.5 * fl1_fx * tez_yy_xx_1[j];

                    tex_yyyy_xy_0[j] = pa_y[j] * tex_yyy_xy_0[j] - pc_y[j] * tex_yyy_xy_1[j] + 1.5 * fl1_fx * tex_yy_xy_0[j] - 1.5 * fl1_fx * tex_yy_xy_1[j] + 0.5 * fl1_fx * tex_yyy_x_0[j] - 0.5 * fl1_fx * tex_yyy_x_1[j];

                    tey_yyyy_xy_0[j] = pa_y[j] * tey_yyy_xy_0[j] - pc_y[j] * tey_yyy_xy_1[j] + 1.5 * fl1_fx * tey_yy_xy_0[j] - 1.5 * fl1_fx * tey_yy_xy_1[j] + 0.5 * fl1_fx * tey_yyy_x_0[j] - 0.5 * fl1_fx * tey_yyy_x_1[j] + ta_yyy_xy_1[j];

                    tez_yyyy_xy_0[j] = pa_y[j] * tez_yyy_xy_0[j] - pc_y[j] * tez_yyy_xy_1[j] + 1.5 * fl1_fx * tez_yy_xy_0[j] - 1.5 * fl1_fx * tez_yy_xy_1[j] + 0.5 * fl1_fx * tez_yyy_x_0[j] - 0.5 * fl1_fx * tez_yyy_x_1[j];

                    tex_yyyy_xz_0[j] = pa_y[j] * tex_yyy_xz_0[j] - pc_y[j] * tex_yyy_xz_1[j] + 1.5 * fl1_fx * tex_yy_xz_0[j] - 1.5 * fl1_fx * tex_yy_xz_1[j];

                    tey_yyyy_xz_0[j] = pa_y[j] * tey_yyy_xz_0[j] - pc_y[j] * tey_yyy_xz_1[j] + 1.5 * fl1_fx * tey_yy_xz_0[j] - 1.5 * fl1_fx * tey_yy_xz_1[j] + ta_yyy_xz_1[j];

                    tez_yyyy_xz_0[j] = pa_y[j] * tez_yyy_xz_0[j] - pc_y[j] * tez_yyy_xz_1[j] + 1.5 * fl1_fx * tez_yy_xz_0[j] - 1.5 * fl1_fx * tez_yy_xz_1[j];

                    tex_yyyy_yy_0[j] = pa_y[j] * tex_yyy_yy_0[j] - pc_y[j] * tex_yyy_yy_1[j] + 1.5 * fl1_fx * tex_yy_yy_0[j] - 1.5 * fl1_fx * tex_yy_yy_1[j] + fl1_fx * tex_yyy_y_0[j] - fl1_fx * tex_yyy_y_1[j];

                    tey_yyyy_yy_0[j] = pa_y[j] * tey_yyy_yy_0[j] - pc_y[j] * tey_yyy_yy_1[j] + 1.5 * fl1_fx * tey_yy_yy_0[j] - 1.5 * fl1_fx * tey_yy_yy_1[j] + fl1_fx * tey_yyy_y_0[j] - fl1_fx * tey_yyy_y_1[j] + ta_yyy_yy_1[j];

                    tez_yyyy_yy_0[j] = pa_y[j] * tez_yyy_yy_0[j] - pc_y[j] * tez_yyy_yy_1[j] + 1.5 * fl1_fx * tez_yy_yy_0[j] - 1.5 * fl1_fx * tez_yy_yy_1[j] + fl1_fx * tez_yyy_y_0[j] - fl1_fx * tez_yyy_y_1[j];

                    tex_yyyy_yz_0[j] = pa_y[j] * tex_yyy_yz_0[j] - pc_y[j] * tex_yyy_yz_1[j] + 1.5 * fl1_fx * tex_yy_yz_0[j] - 1.5 * fl1_fx * tex_yy_yz_1[j] + 0.5 * fl1_fx * tex_yyy_z_0[j] - 0.5 * fl1_fx * tex_yyy_z_1[j];

                    tey_yyyy_yz_0[j] = pa_y[j] * tey_yyy_yz_0[j] - pc_y[j] * tey_yyy_yz_1[j] + 1.5 * fl1_fx * tey_yy_yz_0[j] - 1.5 * fl1_fx * tey_yy_yz_1[j] + 0.5 * fl1_fx * tey_yyy_z_0[j] - 0.5 * fl1_fx * tey_yyy_z_1[j] + ta_yyy_yz_1[j];

                    tez_yyyy_yz_0[j] = pa_y[j] * tez_yyy_yz_0[j] - pc_y[j] * tez_yyy_yz_1[j] + 1.5 * fl1_fx * tez_yy_yz_0[j] - 1.5 * fl1_fx * tez_yy_yz_1[j] + 0.5 * fl1_fx * tez_yyy_z_0[j] - 0.5 * fl1_fx * tez_yyy_z_1[j];

                    tex_yyyy_zz_0[j] = pa_y[j] * tex_yyy_zz_0[j] - pc_y[j] * tex_yyy_zz_1[j] + 1.5 * fl1_fx * tex_yy_zz_0[j] - 1.5 * fl1_fx * tex_yy_zz_1[j];

                    tey_yyyy_zz_0[j] = pa_y[j] * tey_yyy_zz_0[j] - pc_y[j] * tey_yyy_zz_1[j] + 1.5 * fl1_fx * tey_yy_zz_0[j] - 1.5 * fl1_fx * tey_yy_zz_1[j] + ta_yyy_zz_1[j];

                    tez_yyyy_zz_0[j] = pa_y[j] * tez_yyy_zz_0[j] - pc_y[j] * tez_yyy_zz_1[j] + 1.5 * fl1_fx * tez_yy_zz_0[j] - 1.5 * fl1_fx * tez_yy_zz_1[j];

                    tex_yyyz_xx_0[j] = pa_y[j] * tex_yyz_xx_0[j] - pc_y[j] * tex_yyz_xx_1[j] + fl1_fx * tex_yz_xx_0[j] - fl1_fx * tex_yz_xx_1[j];

                    tey_yyyz_xx_0[j] = pa_y[j] * tey_yyz_xx_0[j] - pc_y[j] * tey_yyz_xx_1[j] + fl1_fx * tey_yz_xx_0[j] - fl1_fx * tey_yz_xx_1[j] + ta_yyz_xx_1[j];

                    tez_yyyz_xx_0[j] = pa_y[j] * tez_yyz_xx_0[j] - pc_y[j] * tez_yyz_xx_1[j] + fl1_fx * tez_yz_xx_0[j] - fl1_fx * tez_yz_xx_1[j];

                    tex_yyyz_xy_0[j] = pa_y[j] * tex_yyz_xy_0[j] - pc_y[j] * tex_yyz_xy_1[j] + fl1_fx * tex_yz_xy_0[j] - fl1_fx * tex_yz_xy_1[j] + 0.5 * fl1_fx * tex_yyz_x_0[j] - 0.5 * fl1_fx * tex_yyz_x_1[j];

                    tey_yyyz_xy_0[j] = pa_y[j] * tey_yyz_xy_0[j] - pc_y[j] * tey_yyz_xy_1[j] + fl1_fx * tey_yz_xy_0[j] - fl1_fx * tey_yz_xy_1[j] + 0.5 * fl1_fx * tey_yyz_x_0[j] - 0.5 * fl1_fx * tey_yyz_x_1[j] + ta_yyz_xy_1[j];

                    tez_yyyz_xy_0[j] = pa_y[j] * tez_yyz_xy_0[j] - pc_y[j] * tez_yyz_xy_1[j] + fl1_fx * tez_yz_xy_0[j] - fl1_fx * tez_yz_xy_1[j] + 0.5 * fl1_fx * tez_yyz_x_0[j] - 0.5 * fl1_fx * tez_yyz_x_1[j];

                    tex_yyyz_xz_0[j] = pa_y[j] * tex_yyz_xz_0[j] - pc_y[j] * tex_yyz_xz_1[j] + fl1_fx * tex_yz_xz_0[j] - fl1_fx * tex_yz_xz_1[j];

                    tey_yyyz_xz_0[j] = pa_y[j] * tey_yyz_xz_0[j] - pc_y[j] * tey_yyz_xz_1[j] + fl1_fx * tey_yz_xz_0[j] - fl1_fx * tey_yz_xz_1[j] + ta_yyz_xz_1[j];

                    tez_yyyz_xz_0[j] = pa_y[j] * tez_yyz_xz_0[j] - pc_y[j] * tez_yyz_xz_1[j] + fl1_fx * tez_yz_xz_0[j] - fl1_fx * tez_yz_xz_1[j];

                    tex_yyyz_yy_0[j] = pa_y[j] * tex_yyz_yy_0[j] - pc_y[j] * tex_yyz_yy_1[j] + fl1_fx * tex_yz_yy_0[j] - fl1_fx * tex_yz_yy_1[j] + fl1_fx * tex_yyz_y_0[j] - fl1_fx * tex_yyz_y_1[j];

                    tey_yyyz_yy_0[j] = pa_y[j] * tey_yyz_yy_0[j] - pc_y[j] * tey_yyz_yy_1[j] + fl1_fx * tey_yz_yy_0[j] - fl1_fx * tey_yz_yy_1[j] + fl1_fx * tey_yyz_y_0[j] - fl1_fx * tey_yyz_y_1[j] + ta_yyz_yy_1[j];

                    tez_yyyz_yy_0[j] = pa_y[j] * tez_yyz_yy_0[j] - pc_y[j] * tez_yyz_yy_1[j] + fl1_fx * tez_yz_yy_0[j] - fl1_fx * tez_yz_yy_1[j] + fl1_fx * tez_yyz_y_0[j] - fl1_fx * tez_yyz_y_1[j];

                    tex_yyyz_yz_0[j] = pa_y[j] * tex_yyz_yz_0[j] - pc_y[j] * tex_yyz_yz_1[j] + fl1_fx * tex_yz_yz_0[j] - fl1_fx * tex_yz_yz_1[j] + 0.5 * fl1_fx * tex_yyz_z_0[j] - 0.5 * fl1_fx * tex_yyz_z_1[j];

                    tey_yyyz_yz_0[j] = pa_y[j] * tey_yyz_yz_0[j] - pc_y[j] * tey_yyz_yz_1[j] + fl1_fx * tey_yz_yz_0[j] - fl1_fx * tey_yz_yz_1[j] + 0.5 * fl1_fx * tey_yyz_z_0[j] - 0.5 * fl1_fx * tey_yyz_z_1[j] + ta_yyz_yz_1[j];

                    tez_yyyz_yz_0[j] = pa_y[j] * tez_yyz_yz_0[j] - pc_y[j] * tez_yyz_yz_1[j] + fl1_fx * tez_yz_yz_0[j] - fl1_fx * tez_yz_yz_1[j] + 0.5 * fl1_fx * tez_yyz_z_0[j] - 0.5 * fl1_fx * tez_yyz_z_1[j];

                    tex_yyyz_zz_0[j] = pa_y[j] * tex_yyz_zz_0[j] - pc_y[j] * tex_yyz_zz_1[j] + fl1_fx * tex_yz_zz_0[j] - fl1_fx * tex_yz_zz_1[j];

                    tey_yyyz_zz_0[j] = pa_y[j] * tey_yyz_zz_0[j] - pc_y[j] * tey_yyz_zz_1[j] + fl1_fx * tey_yz_zz_0[j] - fl1_fx * tey_yz_zz_1[j] + ta_yyz_zz_1[j];

                    tez_yyyz_zz_0[j] = pa_y[j] * tez_yyz_zz_0[j] - pc_y[j] * tez_yyz_zz_1[j] + fl1_fx * tez_yz_zz_0[j] - fl1_fx * tez_yz_zz_1[j];

                    tex_yyzz_xx_0[j] = pa_y[j] * tex_yzz_xx_0[j] - pc_y[j] * tex_yzz_xx_1[j] + 0.5 * fl1_fx * tex_zz_xx_0[j] - 0.5 * fl1_fx * tex_zz_xx_1[j];

                    tey_yyzz_xx_0[j] = pa_y[j] * tey_yzz_xx_0[j] - pc_y[j] * tey_yzz_xx_1[j] + 0.5 * fl1_fx * tey_zz_xx_0[j] - 0.5 * fl1_fx * tey_zz_xx_1[j] + ta_yzz_xx_1[j];

                    tez_yyzz_xx_0[j] = pa_y[j] * tez_yzz_xx_0[j] - pc_y[j] * tez_yzz_xx_1[j] + 0.5 * fl1_fx * tez_zz_xx_0[j] - 0.5 * fl1_fx * tez_zz_xx_1[j];

                    tex_yyzz_xy_0[j] = pa_y[j] * tex_yzz_xy_0[j] - pc_y[j] * tex_yzz_xy_1[j] + 0.5 * fl1_fx * tex_zz_xy_0[j] - 0.5 * fl1_fx * tex_zz_xy_1[j] + 0.5 * fl1_fx * tex_yzz_x_0[j] - 0.5 * fl1_fx * tex_yzz_x_1[j];

                    tey_yyzz_xy_0[j] = pa_y[j] * tey_yzz_xy_0[j] - pc_y[j] * tey_yzz_xy_1[j] + 0.5 * fl1_fx * tey_zz_xy_0[j] - 0.5 * fl1_fx * tey_zz_xy_1[j] + 0.5 * fl1_fx * tey_yzz_x_0[j] - 0.5 * fl1_fx * tey_yzz_x_1[j] + ta_yzz_xy_1[j];

                    tez_yyzz_xy_0[j] = pa_y[j] * tez_yzz_xy_0[j] - pc_y[j] * tez_yzz_xy_1[j] + 0.5 * fl1_fx * tez_zz_xy_0[j] - 0.5 * fl1_fx * tez_zz_xy_1[j] + 0.5 * fl1_fx * tez_yzz_x_0[j] - 0.5 * fl1_fx * tez_yzz_x_1[j];

                    tex_yyzz_xz_0[j] = pa_y[j] * tex_yzz_xz_0[j] - pc_y[j] * tex_yzz_xz_1[j] + 0.5 * fl1_fx * tex_zz_xz_0[j] - 0.5 * fl1_fx * tex_zz_xz_1[j];

                    tey_yyzz_xz_0[j] = pa_y[j] * tey_yzz_xz_0[j] - pc_y[j] * tey_yzz_xz_1[j] + 0.5 * fl1_fx * tey_zz_xz_0[j] - 0.5 * fl1_fx * tey_zz_xz_1[j] + ta_yzz_xz_1[j];

                    tez_yyzz_xz_0[j] = pa_y[j] * tez_yzz_xz_0[j] - pc_y[j] * tez_yzz_xz_1[j] + 0.5 * fl1_fx * tez_zz_xz_0[j] - 0.5 * fl1_fx * tez_zz_xz_1[j];
                }

                idx++;
            }
        }
    }

    void
    compElectricFieldForGD_225_270(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto)
    {
        // Batch of Integrals (225,270)

        // set up pointers to primitives data on bra side

        auto spos = braGtoBlock.getStartPositions();

        auto epos = braGtoBlock.getEndPositions();

        auto bdim = epos[iContrGto] - spos[iContrGto];

        // set up pointers to primitives data on ket side

        auto nprim = ketGtoBlock.getNumberOfPrimGtos();

        // set up maximum order of integral

        auto mord = recursionMap.getMaxOrder({"Electric Field"},
                                             {4, -1, -1, -1},
                                             {2, -1, -1, -1},
                                             1, 1);

        for (int32_t iord = 0; iord <= mord; iord++)
        {
            // set up index of integral

            auto pidx_e_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            // check if integral is needed in recursion expansion

            if (pidx_e_4_2_m0 == -1) continue;

            // set up indexes of auxilary integral

            auto pidx_e_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_2_2_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                             1, 1, iord + 1));

            auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord));

            auto pidx_e_3_1_m1 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, 
                                                             {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                             1, 1, iord + 1));

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

                auto tex_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 51); 

                auto tey_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 51); 

                auto tez_yzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 51); 

                auto tex_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 52); 

                auto tey_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 52); 

                auto tez_yzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 52); 

                auto tex_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 53); 

                auto tey_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 53); 

                auto tez_yzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 53); 

                auto tex_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 54); 

                auto tey_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 54); 

                auto tez_zzz_xx_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 54); 

                auto tex_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 55); 

                auto tey_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 55); 

                auto tez_zzz_xy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 55); 

                auto tex_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 56); 

                auto tey_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 56); 

                auto tez_zzz_xz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 56); 

                auto tex_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 57); 

                auto tey_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 57); 

                auto tez_zzz_yy_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 57); 

                auto tex_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 58); 

                auto tey_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 58); 

                auto tez_zzz_yz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 58); 

                auto tex_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * idx + 59); 

                auto tey_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 60 * bdim + 60 * idx + 59); 

                auto tez_zzz_zz_0 = primBuffer.data(pidx_e_3_2_m0 + 120 * bdim + 60 * idx + 59); 

                auto tex_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 51); 

                auto tey_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 51); 

                auto tez_yzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 51); 

                auto tex_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 52); 

                auto tey_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 52); 

                auto tez_yzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 52); 

                auto tex_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 53); 

                auto tey_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 53); 

                auto tez_yzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 53); 

                auto tex_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 54); 

                auto tey_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 54); 

                auto tez_zzz_xx_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 54); 

                auto tex_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 55); 

                auto tey_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 55); 

                auto tez_zzz_xy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 55); 

                auto tex_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 56); 

                auto tey_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 56); 

                auto tez_zzz_xz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 56); 

                auto tex_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 57); 

                auto tey_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 57); 

                auto tez_zzz_yy_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 57); 

                auto tex_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 58); 

                auto tey_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 58); 

                auto tez_zzz_yz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 58); 

                auto tex_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * idx + 59); 

                auto tey_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 60 * bdim + 60 * idx + 59); 

                auto tez_zzz_zz_1 = primBuffer.data(pidx_e_3_2_m1 + 120 * bdim + 60 * idx + 59); 

                auto tex_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 30); 

                auto tey_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 31); 

                auto tey_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 32); 

                auto tey_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 33); 

                auto tey_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 34); 

                auto tey_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * idx + 35); 

                auto tey_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_0 = primBuffer.data(pidx_e_2_2_m0 + 72 * bdim + 36 * idx + 35); 

                auto tex_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 30); 

                auto tey_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 30); 

                auto tez_zz_xx_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 30); 

                auto tex_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 31); 

                auto tey_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 31); 

                auto tez_zz_xy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 31); 

                auto tex_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 32); 

                auto tey_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 32); 

                auto tez_zz_xz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 32); 

                auto tex_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 33); 

                auto tey_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 33); 

                auto tez_zz_yy_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 33); 

                auto tex_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 34); 

                auto tey_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 34); 

                auto tez_zz_yz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 34); 

                auto tex_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * idx + 35); 

                auto tey_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 36 * bdim + 36 * idx + 35); 

                auto tez_zz_zz_1 = primBuffer.data(pidx_e_2_2_m1 + 72 * bdim + 36 * idx + 35); 

                auto tex_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 25); 

                auto tey_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 25); 

                auto tez_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 25); 

                auto tex_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 26); 

                auto tey_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 26); 

                auto tez_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 26); 

                auto tex_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 27); 

                auto tey_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 27); 

                auto tez_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 27); 

                auto tex_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 28); 

                auto tey_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 28); 

                auto tez_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 28); 

                auto tex_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 29); 

                auto tey_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 29); 

                auto tez_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 29); 

                auto tex_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 25); 

                auto tey_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 25); 

                auto tez_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 25); 

                auto tex_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 26); 

                auto tey_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 26); 

                auto tez_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 26); 

                auto tex_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 27); 

                auto tey_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 27); 

                auto tez_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 27); 

                auto tex_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 28); 

                auto tey_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 28); 

                auto tez_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 28); 

                auto tex_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 29); 

                auto tey_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 29); 

                auto tez_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 29); 

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

                auto tex_yyzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 75); 

                auto tey_yyzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 75); 

                auto tez_yyzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 75); 

                auto tex_yyzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 76); 

                auto tey_yyzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 76); 

                auto tez_yyzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 76); 

                auto tex_yyzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 77); 

                auto tey_yyzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 77); 

                auto tez_yyzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 77); 

                auto tex_yzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 78); 

                auto tey_yzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 78); 

                auto tez_yzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 78); 

                auto tex_yzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 79); 

                auto tey_yzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 79); 

                auto tez_yzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 79); 

                auto tex_yzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 80); 

                auto tey_yzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 80); 

                auto tez_yzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 80); 

                auto tex_yzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 81); 

                auto tey_yzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 81); 

                auto tez_yzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 81); 

                auto tex_yzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 82); 

                auto tey_yzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 82); 

                auto tez_yzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 82); 

                auto tex_yzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 83); 

                auto tey_yzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 83); 

                auto tez_yzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 83); 

                auto tex_zzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 84); 

                auto tey_zzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 84); 

                auto tez_zzzz_xx_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 84); 

                auto tex_zzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 85); 

                auto tey_zzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 85); 

                auto tez_zzzz_xy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 85); 

                auto tex_zzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 86); 

                auto tey_zzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 86); 

                auto tez_zzzz_xz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 86); 

                auto tex_zzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 87); 

                auto tey_zzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 87); 

                auto tez_zzzz_yy_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 87); 

                auto tex_zzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 88); 

                auto tey_zzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 88); 

                auto tez_zzzz_yz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 88); 

                auto tex_zzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * idx + 89); 

                auto tey_zzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 90 * bdim + 90 * idx + 89); 

                auto tez_zzzz_zz_0 = primBuffer.data(pidx_e_4_2_m0 + 180 * bdim + 90 * idx + 89); 

                // Batch of Integrals (225,270)

                #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_yzz_yy_1, ta_yzz_yz_1, ta_yzz_zz_1, \
                                         ta_zzz_xx_1, ta_zzz_xy_1, ta_zzz_xz_1, ta_zzz_yy_1, ta_zzz_yz_1, ta_zzz_zz_1, \
                                         tex_yyzz_yy_0, tex_yyzz_yz_0, tex_yyzz_zz_0, tex_yzz_y_0, tex_yzz_y_1, tex_yzz_yy_0, \
                                         tex_yzz_yy_1, tex_yzz_yz_0, tex_yzz_yz_1, tex_yzz_z_0, tex_yzz_z_1, tex_yzz_zz_0, \
                                         tex_yzz_zz_1, tex_yzzz_xx_0, tex_yzzz_xy_0, tex_yzzz_xz_0, tex_yzzz_yy_0, \
                                         tex_yzzz_yz_0, tex_yzzz_zz_0, tex_zz_xx_0, tex_zz_xx_1, tex_zz_xy_0, tex_zz_xy_1, \
                                         tex_zz_xz_0, tex_zz_xz_1, tex_zz_yy_0, tex_zz_yy_1, tex_zz_yz_0, tex_zz_yz_1, \
                                         tex_zz_zz_0, tex_zz_zz_1, tex_zzz_x_0, tex_zzz_x_1, tex_zzz_xx_0, tex_zzz_xx_1, \
                                         tex_zzz_xy_0, tex_zzz_xy_1, tex_zzz_xz_0, tex_zzz_xz_1, tex_zzz_y_0, tex_zzz_y_1, \
                                         tex_zzz_yy_0, tex_zzz_yy_1, tex_zzz_yz_0, tex_zzz_yz_1, tex_zzz_z_0, tex_zzz_z_1, \
                                         tex_zzz_zz_0, tex_zzz_zz_1, tex_zzzz_xx_0, tex_zzzz_xy_0, tex_zzzz_xz_0, \
                                         tex_zzzz_yy_0, tex_zzzz_yz_0, tex_zzzz_zz_0, tey_yyzz_yy_0, tey_yyzz_yz_0, \
                                         tey_yyzz_zz_0, tey_yzz_y_0, tey_yzz_y_1, tey_yzz_yy_0, tey_yzz_yy_1, tey_yzz_yz_0, \
                                         tey_yzz_yz_1, tey_yzz_z_0, tey_yzz_z_1, tey_yzz_zz_0, tey_yzz_zz_1, tey_yzzz_xx_0, \
                                         tey_yzzz_xy_0, tey_yzzz_xz_0, tey_yzzz_yy_0, tey_yzzz_yz_0, tey_yzzz_zz_0, \
                                         tey_zz_xx_0, tey_zz_xx_1, tey_zz_xy_0, tey_zz_xy_1, tey_zz_xz_0, tey_zz_xz_1, \
                                         tey_zz_yy_0, tey_zz_yy_1, tey_zz_yz_0, tey_zz_yz_1, tey_zz_zz_0, tey_zz_zz_1, \
                                         tey_zzz_x_0, tey_zzz_x_1, tey_zzz_xx_0, tey_zzz_xx_1, tey_zzz_xy_0, tey_zzz_xy_1, \
                                         tey_zzz_xz_0, tey_zzz_xz_1, tey_zzz_y_0, tey_zzz_y_1, tey_zzz_yy_0, tey_zzz_yy_1, \
                                         tey_zzz_yz_0, tey_zzz_yz_1, tey_zzz_z_0, tey_zzz_z_1, tey_zzz_zz_0, tey_zzz_zz_1, \
                                         tey_zzzz_xx_0, tey_zzzz_xy_0, tey_zzzz_xz_0, tey_zzzz_yy_0, tey_zzzz_yz_0, \
                                         tey_zzzz_zz_0, tez_yyzz_yy_0, tez_yyzz_yz_0, tez_yyzz_zz_0, tez_yzz_y_0, \
                                         tez_yzz_y_1, tez_yzz_yy_0, tez_yzz_yy_1, tez_yzz_yz_0, tez_yzz_yz_1, tez_yzz_z_0, \
                                         tez_yzz_z_1, tez_yzz_zz_0, tez_yzz_zz_1, tez_yzzz_xx_0, tez_yzzz_xy_0, \
                                         tez_yzzz_xz_0, tez_yzzz_yy_0, tez_yzzz_yz_0, tez_yzzz_zz_0, tez_zz_xx_0, \
                                         tez_zz_xx_1, tez_zz_xy_0, tez_zz_xy_1, tez_zz_xz_0, tez_zz_xz_1, tez_zz_yy_0, \
                                         tez_zz_yy_1, tez_zz_yz_0, tez_zz_yz_1, tez_zz_zz_0, tez_zz_zz_1, tez_zzz_x_0, \
                                         tez_zzz_x_1, tez_zzz_xx_0, tez_zzz_xx_1, tez_zzz_xy_0, tez_zzz_xy_1, tez_zzz_xz_0, \
                                         tez_zzz_xz_1, tez_zzz_y_0, tez_zzz_y_1, tez_zzz_yy_0, tez_zzz_yy_1, tez_zzz_yz_0, \
                                         tez_zzz_yz_1, tez_zzz_z_0, tez_zzz_z_1, tez_zzz_zz_0, tez_zzz_zz_1, tez_zzzz_xx_0, \
                                         tez_zzzz_xy_0, tez_zzzz_xz_0, tez_zzzz_yy_0, tez_zzzz_yz_0, tez_zzzz_zz_0: VLX_ALIGN)
                for (int32_t j = 0; j < nprim; j++)
                {
                    double fl1_fx = fx[j];

                    tex_yyzz_yy_0[j] = pa_y[j] * tex_yzz_yy_0[j] - pc_y[j] * tex_yzz_yy_1[j] + 0.5 * fl1_fx * tex_zz_yy_0[j] - 0.5 * fl1_fx * tex_zz_yy_1[j] + fl1_fx * tex_yzz_y_0[j] - fl1_fx * tex_yzz_y_1[j];

                    tey_yyzz_yy_0[j] = pa_y[j] * tey_yzz_yy_0[j] - pc_y[j] * tey_yzz_yy_1[j] + 0.5 * fl1_fx * tey_zz_yy_0[j] - 0.5 * fl1_fx * tey_zz_yy_1[j] + fl1_fx * tey_yzz_y_0[j] - fl1_fx * tey_yzz_y_1[j] + ta_yzz_yy_1[j];

                    tez_yyzz_yy_0[j] = pa_y[j] * tez_yzz_yy_0[j] - pc_y[j] * tez_yzz_yy_1[j] + 0.5 * fl1_fx * tez_zz_yy_0[j] - 0.5 * fl1_fx * tez_zz_yy_1[j] + fl1_fx * tez_yzz_y_0[j] - fl1_fx * tez_yzz_y_1[j];

                    tex_yyzz_yz_0[j] = pa_y[j] * tex_yzz_yz_0[j] - pc_y[j] * tex_yzz_yz_1[j] + 0.5 * fl1_fx * tex_zz_yz_0[j] - 0.5 * fl1_fx * tex_zz_yz_1[j] + 0.5 * fl1_fx * tex_yzz_z_0[j] - 0.5 * fl1_fx * tex_yzz_z_1[j];

                    tey_yyzz_yz_0[j] = pa_y[j] * tey_yzz_yz_0[j] - pc_y[j] * tey_yzz_yz_1[j] + 0.5 * fl1_fx * tey_zz_yz_0[j] - 0.5 * fl1_fx * tey_zz_yz_1[j] + 0.5 * fl1_fx * tey_yzz_z_0[j] - 0.5 * fl1_fx * tey_yzz_z_1[j] + ta_yzz_yz_1[j];

                    tez_yyzz_yz_0[j] = pa_y[j] * tez_yzz_yz_0[j] - pc_y[j] * tez_yzz_yz_1[j] + 0.5 * fl1_fx * tez_zz_yz_0[j] - 0.5 * fl1_fx * tez_zz_yz_1[j] + 0.5 * fl1_fx * tez_yzz_z_0[j] - 0.5 * fl1_fx * tez_yzz_z_1[j];

                    tex_yyzz_zz_0[j] = pa_y[j] * tex_yzz_zz_0[j] - pc_y[j] * tex_yzz_zz_1[j] + 0.5 * fl1_fx * tex_zz_zz_0[j] - 0.5 * fl1_fx * tex_zz_zz_1[j];

                    tey_yyzz_zz_0[j] = pa_y[j] * tey_yzz_zz_0[j] - pc_y[j] * tey_yzz_zz_1[j] + 0.5 * fl1_fx * tey_zz_zz_0[j] - 0.5 * fl1_fx * tey_zz_zz_1[j] + ta_yzz_zz_1[j];

                    tez_yyzz_zz_0[j] = pa_y[j] * tez_yzz_zz_0[j] - pc_y[j] * tez_yzz_zz_1[j] + 0.5 * fl1_fx * tez_zz_zz_0[j] - 0.5 * fl1_fx * tez_zz_zz_1[j];

                    tex_yzzz_xx_0[j] = pa_y[j] * tex_zzz_xx_0[j] - pc_y[j] * tex_zzz_xx_1[j];

                    tey_yzzz_xx_0[j] = pa_y[j] * tey_zzz_xx_0[j] - pc_y[j] * tey_zzz_xx_1[j] + ta_zzz_xx_1[j];

                    tez_yzzz_xx_0[j] = pa_y[j] * tez_zzz_xx_0[j] - pc_y[j] * tez_zzz_xx_1[j];

                    tex_yzzz_xy_0[j] = pa_y[j] * tex_zzz_xy_0[j] - pc_y[j] * tex_zzz_xy_1[j] + 0.5 * fl1_fx * tex_zzz_x_0[j] - 0.5 * fl1_fx * tex_zzz_x_1[j];

                    tey_yzzz_xy_0[j] = pa_y[j] * tey_zzz_xy_0[j] - pc_y[j] * tey_zzz_xy_1[j] + 0.5 * fl1_fx * tey_zzz_x_0[j] - 0.5 * fl1_fx * tey_zzz_x_1[j] + ta_zzz_xy_1[j];

                    tez_yzzz_xy_0[j] = pa_y[j] * tez_zzz_xy_0[j] - pc_y[j] * tez_zzz_xy_1[j] + 0.5 * fl1_fx * tez_zzz_x_0[j] - 0.5 * fl1_fx * tez_zzz_x_1[j];

                    tex_yzzz_xz_0[j] = pa_y[j] * tex_zzz_xz_0[j] - pc_y[j] * tex_zzz_xz_1[j];

                    tey_yzzz_xz_0[j] = pa_y[j] * tey_zzz_xz_0[j] - pc_y[j] * tey_zzz_xz_1[j] + ta_zzz_xz_1[j];

                    tez_yzzz_xz_0[j] = pa_y[j] * tez_zzz_xz_0[j] - pc_y[j] * tez_zzz_xz_1[j];

                    tex_yzzz_yy_0[j] = pa_y[j] * tex_zzz_yy_0[j] - pc_y[j] * tex_zzz_yy_1[j] + fl1_fx * tex_zzz_y_0[j] - fl1_fx * tex_zzz_y_1[j];

                    tey_yzzz_yy_0[j] = pa_y[j] * tey_zzz_yy_0[j] - pc_y[j] * tey_zzz_yy_1[j] + fl1_fx * tey_zzz_y_0[j] - fl1_fx * tey_zzz_y_1[j] + ta_zzz_yy_1[j];

                    tez_yzzz_yy_0[j] = pa_y[j] * tez_zzz_yy_0[j] - pc_y[j] * tez_zzz_yy_1[j] + fl1_fx * tez_zzz_y_0[j] - fl1_fx * tez_zzz_y_1[j];

                    tex_yzzz_yz_0[j] = pa_y[j] * tex_zzz_yz_0[j] - pc_y[j] * tex_zzz_yz_1[j] + 0.5 * fl1_fx * tex_zzz_z_0[j] - 0.5 * fl1_fx * tex_zzz_z_1[j];

                    tey_yzzz_yz_0[j] = pa_y[j] * tey_zzz_yz_0[j] - pc_y[j] * tey_zzz_yz_1[j] + 0.5 * fl1_fx * tey_zzz_z_0[j] - 0.5 * fl1_fx * tey_zzz_z_1[j] + ta_zzz_yz_1[j];

                    tez_yzzz_yz_0[j] = pa_y[j] * tez_zzz_yz_0[j] - pc_y[j] * tez_zzz_yz_1[j] + 0.5 * fl1_fx * tez_zzz_z_0[j] - 0.5 * fl1_fx * tez_zzz_z_1[j];

                    tex_yzzz_zz_0[j] = pa_y[j] * tex_zzz_zz_0[j] - pc_y[j] * tex_zzz_zz_1[j];

                    tey_yzzz_zz_0[j] = pa_y[j] * tey_zzz_zz_0[j] - pc_y[j] * tey_zzz_zz_1[j] + ta_zzz_zz_1[j];

                    tez_yzzz_zz_0[j] = pa_y[j] * tez_zzz_zz_0[j] - pc_y[j] * tez_zzz_zz_1[j];

                    tex_zzzz_xx_0[j] = pa_z[j] * tex_zzz_xx_0[j] - pc_z[j] * tex_zzz_xx_1[j] + 1.5 * fl1_fx * tex_zz_xx_0[j] - 1.5 * fl1_fx * tex_zz_xx_1[j];

                    tey_zzzz_xx_0[j] = pa_z[j] * tey_zzz_xx_0[j] - pc_z[j] * tey_zzz_xx_1[j] + 1.5 * fl1_fx * tey_zz_xx_0[j] - 1.5 * fl1_fx * tey_zz_xx_1[j];

                    tez_zzzz_xx_0[j] = pa_z[j] * tez_zzz_xx_0[j] - pc_z[j] * tez_zzz_xx_1[j] + 1.5 * fl1_fx * tez_zz_xx_0[j] - 1.5 * fl1_fx * tez_zz_xx_1[j] + ta_zzz_xx_1[j];

                    tex_zzzz_xy_0[j] = pa_z[j] * tex_zzz_xy_0[j] - pc_z[j] * tex_zzz_xy_1[j] + 1.5 * fl1_fx * tex_zz_xy_0[j] - 1.5 * fl1_fx * tex_zz_xy_1[j];

                    tey_zzzz_xy_0[j] = pa_z[j] * tey_zzz_xy_0[j] - pc_z[j] * tey_zzz_xy_1[j] + 1.5 * fl1_fx * tey_zz_xy_0[j] - 1.5 * fl1_fx * tey_zz_xy_1[j];

                    tez_zzzz_xy_0[j] = pa_z[j] * tez_zzz_xy_0[j] - pc_z[j] * tez_zzz_xy_1[j] + 1.5 * fl1_fx * tez_zz_xy_0[j] - 1.5 * fl1_fx * tez_zz_xy_1[j] + ta_zzz_xy_1[j];

                    tex_zzzz_xz_0[j] = pa_z[j] * tex_zzz_xz_0[j] - pc_z[j] * tex_zzz_xz_1[j] + 1.5 * fl1_fx * tex_zz_xz_0[j] - 1.5 * fl1_fx * tex_zz_xz_1[j] + 0.5 * fl1_fx * tex_zzz_x_0[j] - 0.5 * fl1_fx * tex_zzz_x_1[j];

                    tey_zzzz_xz_0[j] = pa_z[j] * tey_zzz_xz_0[j] - pc_z[j] * tey_zzz_xz_1[j] + 1.5 * fl1_fx * tey_zz_xz_0[j] - 1.5 * fl1_fx * tey_zz_xz_1[j] + 0.5 * fl1_fx * tey_zzz_x_0[j] - 0.5 * fl1_fx * tey_zzz_x_1[j];

                    tez_zzzz_xz_0[j] = pa_z[j] * tez_zzz_xz_0[j] - pc_z[j] * tez_zzz_xz_1[j] + 1.5 * fl1_fx * tez_zz_xz_0[j] - 1.5 * fl1_fx * tez_zz_xz_1[j] + 0.5 * fl1_fx * tez_zzz_x_0[j] - 0.5 * fl1_fx * tez_zzz_x_1[j] + ta_zzz_xz_1[j];

                    tex_zzzz_yy_0[j] = pa_z[j] * tex_zzz_yy_0[j] - pc_z[j] * tex_zzz_yy_1[j] + 1.5 * fl1_fx * tex_zz_yy_0[j] - 1.5 * fl1_fx * tex_zz_yy_1[j];

                    tey_zzzz_yy_0[j] = pa_z[j] * tey_zzz_yy_0[j] - pc_z[j] * tey_zzz_yy_1[j] + 1.5 * fl1_fx * tey_zz_yy_0[j] - 1.5 * fl1_fx * tey_zz_yy_1[j];

                    tez_zzzz_yy_0[j] = pa_z[j] * tez_zzz_yy_0[j] - pc_z[j] * tez_zzz_yy_1[j] + 1.5 * fl1_fx * tez_zz_yy_0[j] - 1.5 * fl1_fx * tez_zz_yy_1[j] + ta_zzz_yy_1[j];

                    tex_zzzz_yz_0[j] = pa_z[j] * tex_zzz_yz_0[j] - pc_z[j] * tex_zzz_yz_1[j] + 1.5 * fl1_fx * tex_zz_yz_0[j] - 1.5 * fl1_fx * tex_zz_yz_1[j] + 0.5 * fl1_fx * tex_zzz_y_0[j] - 0.5 * fl1_fx * tex_zzz_y_1[j];

                    tey_zzzz_yz_0[j] = pa_z[j] * tey_zzz_yz_0[j] - pc_z[j] * tey_zzz_yz_1[j] + 1.5 * fl1_fx * tey_zz_yz_0[j] - 1.5 * fl1_fx * tey_zz_yz_1[j] + 0.5 * fl1_fx * tey_zzz_y_0[j] - 0.5 * fl1_fx * tey_zzz_y_1[j];

                    tez_zzzz_yz_0[j] = pa_z[j] * tez_zzz_yz_0[j] - pc_z[j] * tez_zzz_yz_1[j] + 1.5 * fl1_fx * tez_zz_yz_0[j] - 1.5 * fl1_fx * tez_zz_yz_1[j] + 0.5 * fl1_fx * tez_zzz_y_0[j] - 0.5 * fl1_fx * tez_zzz_y_1[j] + ta_zzz_yz_1[j];

                    tex_zzzz_zz_0[j] = pa_z[j] * tex_zzz_zz_0[j] - pc_z[j] * tex_zzz_zz_1[j] + 1.5 * fl1_fx * tex_zz_zz_0[j] - 1.5 * fl1_fx * tex_zz_zz_1[j] + fl1_fx * tex_zzz_z_0[j] - fl1_fx * tex_zzz_z_1[j];

                    tey_zzzz_zz_0[j] = pa_z[j] * tey_zzz_zz_0[j] - pc_z[j] * tey_zzz_zz_1[j] + 1.5 * fl1_fx * tey_zz_zz_0[j] - 1.5 * fl1_fx * tey_zz_zz_1[j] + fl1_fx * tey_zzz_z_0[j] - fl1_fx * tey_zzz_z_1[j];

                    tez_zzzz_zz_0[j] = pa_z[j] * tez_zzz_zz_0[j] - pc_z[j] * tez_zzz_zz_1[j] + 1.5 * fl1_fx * tez_zz_zz_0[j] - 1.5 * fl1_fx * tez_zz_zz_1[j] + fl1_fx * tez_zzz_z_0[j] - fl1_fx * tez_zzz_z_1[j] + ta_zzz_zz_1[j];
                }

                idx++;
            }
        }
    }


} // efieldrecfunc namespace

