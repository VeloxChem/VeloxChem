//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AngularMomentumRecFuncForDX.hpp"

namespace amomrecfunc { // amomrecfunc namespace

    void
    compAngularMomentumForDD(      CMemBlock2D<double>& primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        amomrecfunc::compAngularMomentumForDD_0_36(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        amomrecfunc::compAngularMomentumForDD_36_72(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        amomrecfunc::compAngularMomentumForDD_72_108(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 
    }

    void
    compAngularMomentumForDD_0_36(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx); 

            auto tly_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tlz_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tlx_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 1); 

            auto tly_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tlz_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tlx_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 2); 

            auto tly_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tlz_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tlx_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 3); 

            auto tly_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tlz_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tlx_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 4); 

            auto tly_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tlz_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tlx_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 5); 

            auto tly_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tlz_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tlx_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 6); 

            auto tly_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tlz_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tlx_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 7); 

            auto tly_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tlz_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tlx_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 8); 

            auto tly_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tlz_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tlx_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 9); 

            auto tly_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 10); 

            auto tly_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tlx_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 11); 

            auto tly_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx); 

            auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx); 

            auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx); 

            auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1); 

            auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1); 

            auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1); 

            auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2); 

            auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2); 

            auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2); 

            auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3); 

            auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3); 

            auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3); 

            auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4); 

            auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4); 

            auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4); 

            auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5); 

            auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5); 

            auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5); 

            auto tlx_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx); 

            auto tly_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx); 

            auto tlz_x_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx); 

            auto tlx_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 1); 

            auto tly_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 1); 

            auto tlz_x_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 1); 

            auto tlx_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 2); 

            auto tly_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 2); 

            auto tlz_x_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 2); 

            auto tlx_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 3); 

            auto tly_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 3); 

            auto tlz_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 3); 

            auto tlx_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 4); 

            auto tly_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 4); 

            auto tlz_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 4); 

            auto tlx_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 5); 

            auto tly_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 5); 

            auto tlz_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 5); 

            auto tpy_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tpz_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tpy_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tpz_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tpy_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tpz_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tpy_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tpz_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tpy_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tpz_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tpy_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tpz_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdy_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tdz_x_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tdy_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tdz_x_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tdy_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tdz_x_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tdy_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tdz_x_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tdy_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tdz_x_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tdy_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tdz_x_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tdy_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdy_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdy_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdy_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdy_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdy_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            // set up pointers to integrals

            auto tlx_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx); 

            auto tly_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tlz_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tlx_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 1); 

            auto tly_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tlz_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tlx_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 2); 

            auto tly_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tlz_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tlx_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 3); 

            auto tly_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tlz_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tlx_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 4); 

            auto tly_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tlz_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tlx_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 5); 

            auto tly_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tlz_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tlx_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 6); 

            auto tly_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tlz_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tlx_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 7); 

            auto tly_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tlz_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tlx_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 8); 

            auto tly_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tlz_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tlx_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 9); 

            auto tly_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tlz_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tlx_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 10); 

            auto tly_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tlz_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tlx_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 11); 

            auto tly_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tlz_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            // Batch of Integrals (0,36)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_x_xx_0, tdy_x_xy_0, tdy_x_xz_0, tdy_x_yy_0, tdy_x_yz_0, \
                                     tdy_x_zz_0, tdy_y_xx_0, tdy_y_xy_0, tdy_y_xz_0, tdy_y_yy_0, tdy_y_yz_0, tdy_y_zz_0, \
                                     tdz_x_xx_0, tdz_x_xy_0, tdz_x_xz_0, tdz_x_yy_0, tdz_x_yz_0, tdz_x_zz_0, tdz_y_xx_0, \
                                     tdz_y_xy_0, tdz_y_xz_0, tdz_y_yy_0, tdz_y_yz_0, tdz_y_zz_0, tlx_0_xx_0, tlx_0_xy_0, \
                                     tlx_0_xz_0, tlx_0_yy_0, tlx_0_yz_0, tlx_0_zz_0, tlx_x_x_0, tlx_x_xx_0, tlx_x_xy_0, \
                                     tlx_x_xz_0, tlx_x_y_0, tlx_x_yy_0, tlx_x_yz_0, tlx_x_z_0, tlx_x_zz_0, tlx_xx_xx_0, \
                                     tlx_xx_xy_0, tlx_xx_xz_0, tlx_xx_yy_0, tlx_xx_yz_0, tlx_xx_zz_0, tlx_xy_xx_0, \
                                     tlx_xy_xy_0, tlx_xy_xz_0, tlx_xy_yy_0, tlx_xy_yz_0, tlx_xy_zz_0, tlx_y_x_0, \
                                     tlx_y_xx_0, tlx_y_xy_0, tlx_y_xz_0, tlx_y_y_0, tlx_y_yy_0, tlx_y_yz_0, tlx_y_z_0, \
                                     tlx_y_zz_0, tly_0_xx_0, tly_0_xy_0, tly_0_xz_0, tly_0_yy_0, tly_0_yz_0, tly_0_zz_0, \
                                     tly_x_x_0, tly_x_xx_0, tly_x_xy_0, tly_x_xz_0, tly_x_y_0, tly_x_yy_0, tly_x_yz_0, \
                                     tly_x_z_0, tly_x_zz_0, tly_xx_xx_0, tly_xx_xy_0, tly_xx_xz_0, tly_xx_yy_0, \
                                     tly_xx_yz_0, tly_xx_zz_0, tly_xy_xx_0, tly_xy_xy_0, tly_xy_xz_0, tly_xy_yy_0, \
                                     tly_xy_yz_0, tly_xy_zz_0, tly_y_x_0, tly_y_xx_0, tly_y_xy_0, tly_y_xz_0, tly_y_y_0, \
                                     tly_y_yy_0, tly_y_yz_0, tly_y_z_0, tly_y_zz_0, tlz_0_xx_0, tlz_0_xy_0, tlz_0_xz_0, \
                                     tlz_0_yy_0, tlz_0_yz_0, tlz_0_zz_0, tlz_x_x_0, tlz_x_xx_0, tlz_x_xy_0, tlz_x_xz_0, \
                                     tlz_x_y_0, tlz_x_yy_0, tlz_x_yz_0, tlz_x_z_0, tlz_x_zz_0, tlz_xx_xx_0, \
                                     tlz_xx_xy_0, tlz_xx_xz_0, tlz_xx_yy_0, tlz_xx_yz_0, tlz_xx_zz_0, tlz_xy_xx_0, \
                                     tlz_xy_xy_0, tlz_xy_xz_0, tlz_xy_yy_0, tlz_xy_yz_0, tlz_xy_zz_0, tlz_y_x_0, \
                                     tlz_y_xx_0, tlz_y_xy_0, tlz_y_xz_0, tlz_y_y_0, tlz_y_yy_0, tlz_y_yz_0, tlz_y_z_0, \
                                     tlz_y_zz_0, tpy_x_xx_0, tpy_x_xy_0, tpy_x_xz_0, tpy_x_yy_0, tpy_x_yz_0, tpy_x_zz_0, \
                                     tpy_y_xx_0, tpy_y_xy_0, tpy_y_xz_0, tpy_y_yy_0, tpy_y_yz_0, tpy_y_zz_0, tpz_x_xx_0, \
                                     tpz_x_xy_0, tpz_x_xz_0, tpz_x_yy_0, tpz_x_yz_0, tpz_x_zz_0, tpz_y_xx_0, tpz_y_xy_0, \
                                     tpz_y_xz_0, tpz_y_yy_0, tpz_y_yz_0, tpz_y_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xx_xx_0[j] = pa_x[j] * tlx_x_xx_0[j] + 0.5 * fl1_fx * tlx_0_xx_0[j] + fl1_fx * tlx_x_x_0[j];

                tly_xx_xx_0[j] = pa_x[j] * tly_x_xx_0[j] + 0.5 * fl1_fx * tly_0_xx_0[j] + fl1_fx * tly_x_x_0[j] + 0.5 * fl1_fx * tpz_x_xx_0[j] + fl1_fx * fl1_fgb * tdz_x_xx_0[j];

                tlz_xx_xx_0[j] = pa_x[j] * tlz_x_xx_0[j] + 0.5 * fl1_fx * tlz_0_xx_0[j] + fl1_fx * tlz_x_x_0[j] - 0.5 * fl1_fx * tpy_x_xx_0[j] - fl1_fx * fl1_fgb * tdy_x_xx_0[j];

                tlx_xx_xy_0[j] = pa_x[j] * tlx_x_xy_0[j] + 0.5 * fl1_fx * tlx_0_xy_0[j] + 0.5 * fl1_fx * tlx_x_y_0[j];

                tly_xx_xy_0[j] = pa_x[j] * tly_x_xy_0[j] + 0.5 * fl1_fx * tly_0_xy_0[j] + 0.5 * fl1_fx * tly_x_y_0[j] + 0.5 * fl1_fx * tpz_x_xy_0[j] + fl1_fx * fl1_fgb * tdz_x_xy_0[j];

                tlz_xx_xy_0[j] = pa_x[j] * tlz_x_xy_0[j] + 0.5 * fl1_fx * tlz_0_xy_0[j] + 0.5 * fl1_fx * tlz_x_y_0[j] - 0.5 * fl1_fx * tpy_x_xy_0[j] - fl1_fx * fl1_fgb * tdy_x_xy_0[j];

                tlx_xx_xz_0[j] = pa_x[j] * tlx_x_xz_0[j] + 0.5 * fl1_fx * tlx_0_xz_0[j] + 0.5 * fl1_fx * tlx_x_z_0[j];

                tly_xx_xz_0[j] = pa_x[j] * tly_x_xz_0[j] + 0.5 * fl1_fx * tly_0_xz_0[j] + 0.5 * fl1_fx * tly_x_z_0[j] + 0.5 * fl1_fx * tpz_x_xz_0[j] + fl1_fx * fl1_fgb * tdz_x_xz_0[j];

                tlz_xx_xz_0[j] = pa_x[j] * tlz_x_xz_0[j] + 0.5 * fl1_fx * tlz_0_xz_0[j] + 0.5 * fl1_fx * tlz_x_z_0[j] - 0.5 * fl1_fx * tpy_x_xz_0[j] - fl1_fx * fl1_fgb * tdy_x_xz_0[j];

                tlx_xx_yy_0[j] = pa_x[j] * tlx_x_yy_0[j] + 0.5 * fl1_fx * tlx_0_yy_0[j];

                tly_xx_yy_0[j] = pa_x[j] * tly_x_yy_0[j] + 0.5 * fl1_fx * tly_0_yy_0[j] + 0.5 * fl1_fx * tpz_x_yy_0[j] + fl1_fx * fl1_fgb * tdz_x_yy_0[j];

                tlz_xx_yy_0[j] = pa_x[j] * tlz_x_yy_0[j] + 0.5 * fl1_fx * tlz_0_yy_0[j] - 0.5 * fl1_fx * tpy_x_yy_0[j] - fl1_fx * fl1_fgb * tdy_x_yy_0[j];

                tlx_xx_yz_0[j] = pa_x[j] * tlx_x_yz_0[j] + 0.5 * fl1_fx * tlx_0_yz_0[j];

                tly_xx_yz_0[j] = pa_x[j] * tly_x_yz_0[j] + 0.5 * fl1_fx * tly_0_yz_0[j] + 0.5 * fl1_fx * tpz_x_yz_0[j] + fl1_fx * fl1_fgb * tdz_x_yz_0[j];

                tlz_xx_yz_0[j] = pa_x[j] * tlz_x_yz_0[j] + 0.5 * fl1_fx * tlz_0_yz_0[j] - 0.5 * fl1_fx * tpy_x_yz_0[j] - fl1_fx * fl1_fgb * tdy_x_yz_0[j];

                tlx_xx_zz_0[j] = pa_x[j] * tlx_x_zz_0[j] + 0.5 * fl1_fx * tlx_0_zz_0[j];

                tly_xx_zz_0[j] = pa_x[j] * tly_x_zz_0[j] + 0.5 * fl1_fx * tly_0_zz_0[j] + 0.5 * fl1_fx * tpz_x_zz_0[j] + fl1_fx * fl1_fgb * tdz_x_zz_0[j];

                tlz_xx_zz_0[j] = pa_x[j] * tlz_x_zz_0[j] + 0.5 * fl1_fx * tlz_0_zz_0[j] - 0.5 * fl1_fx * tpy_x_zz_0[j] - fl1_fx * fl1_fgb * tdy_x_zz_0[j];

                tlx_xy_xx_0[j] = pa_x[j] * tlx_y_xx_0[j] + fl1_fx * tlx_y_x_0[j];

                tly_xy_xx_0[j] = pa_x[j] * tly_y_xx_0[j] + fl1_fx * tly_y_x_0[j] + 0.5 * fl1_fx * tpz_y_xx_0[j] + fl1_fx * fl1_fgb * tdz_y_xx_0[j];

                tlz_xy_xx_0[j] = pa_x[j] * tlz_y_xx_0[j] + fl1_fx * tlz_y_x_0[j] - 0.5 * fl1_fx * tpy_y_xx_0[j] - fl1_fx * fl1_fgb * tdy_y_xx_0[j];

                tlx_xy_xy_0[j] = pa_x[j] * tlx_y_xy_0[j] + 0.5 * fl1_fx * tlx_y_y_0[j];

                tly_xy_xy_0[j] = pa_x[j] * tly_y_xy_0[j] + 0.5 * fl1_fx * tly_y_y_0[j] + 0.5 * fl1_fx * tpz_y_xy_0[j] + fl1_fx * fl1_fgb * tdz_y_xy_0[j];

                tlz_xy_xy_0[j] = pa_x[j] * tlz_y_xy_0[j] + 0.5 * fl1_fx * tlz_y_y_0[j] - 0.5 * fl1_fx * tpy_y_xy_0[j] - fl1_fx * fl1_fgb * tdy_y_xy_0[j];

                tlx_xy_xz_0[j] = pa_x[j] * tlx_y_xz_0[j] + 0.5 * fl1_fx * tlx_y_z_0[j];

                tly_xy_xz_0[j] = pa_x[j] * tly_y_xz_0[j] + 0.5 * fl1_fx * tly_y_z_0[j] + 0.5 * fl1_fx * tpz_y_xz_0[j] + fl1_fx * fl1_fgb * tdz_y_xz_0[j];

                tlz_xy_xz_0[j] = pa_x[j] * tlz_y_xz_0[j] + 0.5 * fl1_fx * tlz_y_z_0[j] - 0.5 * fl1_fx * tpy_y_xz_0[j] - fl1_fx * fl1_fgb * tdy_y_xz_0[j];

                tlx_xy_yy_0[j] = pa_x[j] * tlx_y_yy_0[j];

                tly_xy_yy_0[j] = pa_x[j] * tly_y_yy_0[j] + 0.5 * fl1_fx * tpz_y_yy_0[j] + fl1_fx * fl1_fgb * tdz_y_yy_0[j];

                tlz_xy_yy_0[j] = pa_x[j] * tlz_y_yy_0[j] - 0.5 * fl1_fx * tpy_y_yy_0[j] - fl1_fx * fl1_fgb * tdy_y_yy_0[j];

                tlx_xy_yz_0[j] = pa_x[j] * tlx_y_yz_0[j];

                tly_xy_yz_0[j] = pa_x[j] * tly_y_yz_0[j] + 0.5 * fl1_fx * tpz_y_yz_0[j] + fl1_fx * fl1_fgb * tdz_y_yz_0[j];

                tlz_xy_yz_0[j] = pa_x[j] * tlz_y_yz_0[j] - 0.5 * fl1_fx * tpy_y_yz_0[j] - fl1_fx * fl1_fgb * tdy_y_yz_0[j];

                tlx_xy_zz_0[j] = pa_x[j] * tlx_y_zz_0[j];

                tly_xy_zz_0[j] = pa_x[j] * tly_y_zz_0[j] + 0.5 * fl1_fx * tpz_y_zz_0[j] + fl1_fx * fl1_fgb * tdz_y_zz_0[j];

                tlz_xy_zz_0[j] = pa_x[j] * tlz_y_zz_0[j] - 0.5 * fl1_fx * tpy_y_zz_0[j] - fl1_fx * fl1_fgb * tdy_y_zz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDD_36_72(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tlx_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 6); 

            auto tly_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tlz_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tlx_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 7); 

            auto tly_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tlz_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tlx_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 8); 

            auto tly_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tlz_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tlx_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 9); 

            auto tly_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 10); 

            auto tly_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tlx_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 11); 

            auto tly_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 15); 

            auto tly_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 16); 

            auto tly_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 17); 

            auto tly_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx); 

            auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx); 

            auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx); 

            auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1); 

            auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1); 

            auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1); 

            auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2); 

            auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2); 

            auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2); 

            auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3); 

            auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3); 

            auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3); 

            auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4); 

            auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4); 

            auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4); 

            auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5); 

            auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5); 

            auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5); 

            auto tlx_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 3); 

            auto tly_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 3); 

            auto tlz_y_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 3); 

            auto tlx_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 4); 

            auto tly_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 4); 

            auto tlz_y_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 4); 

            auto tlx_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 5); 

            auto tly_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 5); 

            auto tlz_y_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 5); 

            auto tlx_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 6); 

            auto tly_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 6); 

            auto tlz_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 6); 

            auto tlx_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 7); 

            auto tly_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 7); 

            auto tlz_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 7); 

            auto tlx_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 8); 

            auto tly_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 8); 

            auto tlz_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 8); 

            auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6); 

            auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7); 

            auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8); 

            auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9); 

            auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10); 

            auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11); 

            auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tdx_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 6); 

            auto tdz_y_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tdx_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 7); 

            auto tdz_y_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tdx_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 8); 

            auto tdz_y_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tdx_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 9); 

            auto tdz_y_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tdx_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 10); 

            auto tdz_y_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tdx_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 11); 

            auto tdz_y_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tdy_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tdz_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tdy_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tdz_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tdy_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tdz_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tdy_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdy_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdy_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            // set up pointers to integrals

            auto tlx_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 12); 

            auto tly_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tlz_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tlx_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 13); 

            auto tly_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tlz_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tlx_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 14); 

            auto tly_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tlz_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            auto tlx_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 15); 

            auto tly_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tlz_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tlx_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 16); 

            auto tly_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tlz_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tlx_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 17); 

            auto tly_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tlz_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tlx_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 18); 

            auto tly_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tlz_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tlx_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 19); 

            auto tly_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tlz_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tlx_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 20); 

            auto tly_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tlz_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tlx_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 21); 

            auto tly_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tlz_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tlx_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 22); 

            auto tly_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tlz_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tlx_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 23); 

            auto tly_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tlz_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            // Batch of Integrals (36,72)

            #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tdx_y_xx_0, tdx_y_xy_0, tdx_y_xz_0, tdx_y_yy_0, \
                                     tdx_y_yz_0, tdx_y_zz_0, tdy_z_xx_0, tdy_z_xy_0, tdy_z_xz_0, tdy_z_yy_0, tdy_z_yz_0, \
                                     tdy_z_zz_0, tdz_y_xx_0, tdz_y_xy_0, tdz_y_xz_0, tdz_y_yy_0, tdz_y_yz_0, tdz_y_zz_0, \
                                     tdz_z_xx_0, tdz_z_xy_0, tdz_z_xz_0, tdz_z_yy_0, tdz_z_yz_0, tdz_z_zz_0, tlx_0_xx_0, \
                                     tlx_0_xy_0, tlx_0_xz_0, tlx_0_yy_0, tlx_0_yz_0, tlx_0_zz_0, tlx_xz_xx_0, \
                                     tlx_xz_xy_0, tlx_xz_xz_0, tlx_xz_yy_0, tlx_xz_yz_0, tlx_xz_zz_0, tlx_y_x_0, \
                                     tlx_y_xx_0, tlx_y_xy_0, tlx_y_xz_0, tlx_y_y_0, tlx_y_yy_0, tlx_y_yz_0, tlx_y_z_0, \
                                     tlx_y_zz_0, tlx_yy_xx_0, tlx_yy_xy_0, tlx_yy_xz_0, tlx_yy_yy_0, tlx_yy_yz_0, \
                                     tlx_yy_zz_0, tlx_z_x_0, tlx_z_xx_0, tlx_z_xy_0, tlx_z_xz_0, tlx_z_y_0, tlx_z_yy_0, \
                                     tlx_z_yz_0, tlx_z_z_0, tlx_z_zz_0, tly_0_xx_0, tly_0_xy_0, tly_0_xz_0, tly_0_yy_0, \
                                     tly_0_yz_0, tly_0_zz_0, tly_xz_xx_0, tly_xz_xy_0, tly_xz_xz_0, tly_xz_yy_0, \
                                     tly_xz_yz_0, tly_xz_zz_0, tly_y_x_0, tly_y_xx_0, tly_y_xy_0, tly_y_xz_0, tly_y_y_0, \
                                     tly_y_yy_0, tly_y_yz_0, tly_y_z_0, tly_y_zz_0, tly_yy_xx_0, tly_yy_xy_0, \
                                     tly_yy_xz_0, tly_yy_yy_0, tly_yy_yz_0, tly_yy_zz_0, tly_z_x_0, tly_z_xx_0, \
                                     tly_z_xy_0, tly_z_xz_0, tly_z_y_0, tly_z_yy_0, tly_z_yz_0, tly_z_z_0, tly_z_zz_0, \
                                     tlz_0_xx_0, tlz_0_xy_0, tlz_0_xz_0, tlz_0_yy_0, tlz_0_yz_0, tlz_0_zz_0, \
                                     tlz_xz_xx_0, tlz_xz_xy_0, tlz_xz_xz_0, tlz_xz_yy_0, tlz_xz_yz_0, tlz_xz_zz_0, \
                                     tlz_y_x_0, tlz_y_xx_0, tlz_y_xy_0, tlz_y_xz_0, tlz_y_y_0, tlz_y_yy_0, tlz_y_yz_0, \
                                     tlz_y_z_0, tlz_y_zz_0, tlz_yy_xx_0, tlz_yy_xy_0, tlz_yy_xz_0, tlz_yy_yy_0, \
                                     tlz_yy_yz_0, tlz_yy_zz_0, tlz_z_x_0, tlz_z_xx_0, tlz_z_xy_0, tlz_z_xz_0, tlz_z_y_0, \
                                     tlz_z_yy_0, tlz_z_yz_0, tlz_z_z_0, tlz_z_zz_0, tpx_y_xx_0, tpx_y_xy_0, tpx_y_xz_0, \
                                     tpx_y_yy_0, tpx_y_yz_0, tpx_y_zz_0, tpy_z_xx_0, tpy_z_xy_0, tpy_z_xz_0, tpy_z_yy_0, \
                                     tpy_z_yz_0, tpy_z_zz_0, tpz_y_xx_0, tpz_y_xy_0, tpz_y_xz_0, tpz_y_yy_0, tpz_y_yz_0, \
                                     tpz_y_zz_0, tpz_z_xx_0, tpz_z_xy_0, tpz_z_xz_0, tpz_z_yy_0, tpz_z_yz_0, tpz_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xz_xx_0[j] = pa_x[j] * tlx_z_xx_0[j] + fl1_fx * tlx_z_x_0[j];

                tly_xz_xx_0[j] = pa_x[j] * tly_z_xx_0[j] + fl1_fx * tly_z_x_0[j] + 0.5 * fl1_fx * tpz_z_xx_0[j] + fl1_fx * fl1_fgb * tdz_z_xx_0[j];

                tlz_xz_xx_0[j] = pa_x[j] * tlz_z_xx_0[j] + fl1_fx * tlz_z_x_0[j] - 0.5 * fl1_fx * tpy_z_xx_0[j] - fl1_fx * fl1_fgb * tdy_z_xx_0[j];

                tlx_xz_xy_0[j] = pa_x[j] * tlx_z_xy_0[j] + 0.5 * fl1_fx * tlx_z_y_0[j];

                tly_xz_xy_0[j] = pa_x[j] * tly_z_xy_0[j] + 0.5 * fl1_fx * tly_z_y_0[j] + 0.5 * fl1_fx * tpz_z_xy_0[j] + fl1_fx * fl1_fgb * tdz_z_xy_0[j];

                tlz_xz_xy_0[j] = pa_x[j] * tlz_z_xy_0[j] + 0.5 * fl1_fx * tlz_z_y_0[j] - 0.5 * fl1_fx * tpy_z_xy_0[j] - fl1_fx * fl1_fgb * tdy_z_xy_0[j];

                tlx_xz_xz_0[j] = pa_x[j] * tlx_z_xz_0[j] + 0.5 * fl1_fx * tlx_z_z_0[j];

                tly_xz_xz_0[j] = pa_x[j] * tly_z_xz_0[j] + 0.5 * fl1_fx * tly_z_z_0[j] + 0.5 * fl1_fx * tpz_z_xz_0[j] + fl1_fx * fl1_fgb * tdz_z_xz_0[j];

                tlz_xz_xz_0[j] = pa_x[j] * tlz_z_xz_0[j] + 0.5 * fl1_fx * tlz_z_z_0[j] - 0.5 * fl1_fx * tpy_z_xz_0[j] - fl1_fx * fl1_fgb * tdy_z_xz_0[j];

                tlx_xz_yy_0[j] = pa_x[j] * tlx_z_yy_0[j];

                tly_xz_yy_0[j] = pa_x[j] * tly_z_yy_0[j] + 0.5 * fl1_fx * tpz_z_yy_0[j] + fl1_fx * fl1_fgb * tdz_z_yy_0[j];

                tlz_xz_yy_0[j] = pa_x[j] * tlz_z_yy_0[j] - 0.5 * fl1_fx * tpy_z_yy_0[j] - fl1_fx * fl1_fgb * tdy_z_yy_0[j];

                tlx_xz_yz_0[j] = pa_x[j] * tlx_z_yz_0[j];

                tly_xz_yz_0[j] = pa_x[j] * tly_z_yz_0[j] + 0.5 * fl1_fx * tpz_z_yz_0[j] + fl1_fx * fl1_fgb * tdz_z_yz_0[j];

                tlz_xz_yz_0[j] = pa_x[j] * tlz_z_yz_0[j] - 0.5 * fl1_fx * tpy_z_yz_0[j] - fl1_fx * fl1_fgb * tdy_z_yz_0[j];

                tlx_xz_zz_0[j] = pa_x[j] * tlx_z_zz_0[j];

                tly_xz_zz_0[j] = pa_x[j] * tly_z_zz_0[j] + 0.5 * fl1_fx * tpz_z_zz_0[j] + fl1_fx * fl1_fgb * tdz_z_zz_0[j];

                tlz_xz_zz_0[j] = pa_x[j] * tlz_z_zz_0[j] - 0.5 * fl1_fx * tpy_z_zz_0[j] - fl1_fx * fl1_fgb * tdy_z_zz_0[j];

                tlx_yy_xx_0[j] = pa_y[j] * tlx_y_xx_0[j] + 0.5 * fl1_fx * tlx_0_xx_0[j] - 0.5 * fl1_fx * tpz_y_xx_0[j] - fl1_fx * fl1_fgb * tdz_y_xx_0[j];

                tly_yy_xx_0[j] = pa_y[j] * tly_y_xx_0[j] + 0.5 * fl1_fx * tly_0_xx_0[j];

                tlz_yy_xx_0[j] = pa_y[j] * tlz_y_xx_0[j] + 0.5 * fl1_fx * tlz_0_xx_0[j] + 0.5 * fl1_fx * tpx_y_xx_0[j] + fl1_fx * fl1_fgb * tdx_y_xx_0[j];

                tlx_yy_xy_0[j] = pa_y[j] * tlx_y_xy_0[j] + 0.5 * fl1_fx * tlx_0_xy_0[j] + 0.5 * fl1_fx * tlx_y_x_0[j] - 0.5 * fl1_fx * tpz_y_xy_0[j] - fl1_fx * fl1_fgb * tdz_y_xy_0[j];

                tly_yy_xy_0[j] = pa_y[j] * tly_y_xy_0[j] + 0.5 * fl1_fx * tly_0_xy_0[j] + 0.5 * fl1_fx * tly_y_x_0[j];

                tlz_yy_xy_0[j] = pa_y[j] * tlz_y_xy_0[j] + 0.5 * fl1_fx * tlz_0_xy_0[j] + 0.5 * fl1_fx * tlz_y_x_0[j] + 0.5 * fl1_fx * tpx_y_xy_0[j] + fl1_fx * fl1_fgb * tdx_y_xy_0[j];

                tlx_yy_xz_0[j] = pa_y[j] * tlx_y_xz_0[j] + 0.5 * fl1_fx * tlx_0_xz_0[j] - 0.5 * fl1_fx * tpz_y_xz_0[j] - fl1_fx * fl1_fgb * tdz_y_xz_0[j];

                tly_yy_xz_0[j] = pa_y[j] * tly_y_xz_0[j] + 0.5 * fl1_fx * tly_0_xz_0[j];

                tlz_yy_xz_0[j] = pa_y[j] * tlz_y_xz_0[j] + 0.5 * fl1_fx * tlz_0_xz_0[j] + 0.5 * fl1_fx * tpx_y_xz_0[j] + fl1_fx * fl1_fgb * tdx_y_xz_0[j];

                tlx_yy_yy_0[j] = pa_y[j] * tlx_y_yy_0[j] + 0.5 * fl1_fx * tlx_0_yy_0[j] + fl1_fx * tlx_y_y_0[j] - 0.5 * fl1_fx * tpz_y_yy_0[j] - fl1_fx * fl1_fgb * tdz_y_yy_0[j];

                tly_yy_yy_0[j] = pa_y[j] * tly_y_yy_0[j] + 0.5 * fl1_fx * tly_0_yy_0[j] + fl1_fx * tly_y_y_0[j];

                tlz_yy_yy_0[j] = pa_y[j] * tlz_y_yy_0[j] + 0.5 * fl1_fx * tlz_0_yy_0[j] + fl1_fx * tlz_y_y_0[j] + 0.5 * fl1_fx * tpx_y_yy_0[j] + fl1_fx * fl1_fgb * tdx_y_yy_0[j];

                tlx_yy_yz_0[j] = pa_y[j] * tlx_y_yz_0[j] + 0.5 * fl1_fx * tlx_0_yz_0[j] + 0.5 * fl1_fx * tlx_y_z_0[j] - 0.5 * fl1_fx * tpz_y_yz_0[j] - fl1_fx * fl1_fgb * tdz_y_yz_0[j];

                tly_yy_yz_0[j] = pa_y[j] * tly_y_yz_0[j] + 0.5 * fl1_fx * tly_0_yz_0[j] + 0.5 * fl1_fx * tly_y_z_0[j];

                tlz_yy_yz_0[j] = pa_y[j] * tlz_y_yz_0[j] + 0.5 * fl1_fx * tlz_0_yz_0[j] + 0.5 * fl1_fx * tlz_y_z_0[j] + 0.5 * fl1_fx * tpx_y_yz_0[j] + fl1_fx * fl1_fgb * tdx_y_yz_0[j];

                tlx_yy_zz_0[j] = pa_y[j] * tlx_y_zz_0[j] + 0.5 * fl1_fx * tlx_0_zz_0[j] - 0.5 * fl1_fx * tpz_y_zz_0[j] - fl1_fx * fl1_fgb * tdz_y_zz_0[j];

                tly_yy_zz_0[j] = pa_y[j] * tly_y_zz_0[j] + 0.5 * fl1_fx * tly_0_zz_0[j];

                tlz_yy_zz_0[j] = pa_y[j] * tlz_y_zz_0[j] + 0.5 * fl1_fx * tlz_0_zz_0[j] + 0.5 * fl1_fx * tpx_y_zz_0[j] + fl1_fx * fl1_fgb * tdx_y_zz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDD_72_108(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 15); 

            auto tly_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 16); 

            auto tly_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 17); 

            auto tly_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx); 

            auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx); 

            auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx); 

            auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1); 

            auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1); 

            auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1); 

            auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2); 

            auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2); 

            auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2); 

            auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3); 

            auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3); 

            auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3); 

            auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4); 

            auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4); 

            auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4); 

            auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5); 

            auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5); 

            auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5); 

            auto tlx_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 6); 

            auto tly_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 6); 

            auto tlz_z_x_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 6); 

            auto tlx_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 7); 

            auto tly_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 7); 

            auto tlz_z_y_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 7); 

            auto tlx_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * idx + 8); 

            auto tly_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 9 * bdim + 9 * idx + 8); 

            auto tlz_z_z_0 = primBuffer.data(pidx_l_1_1_m0 + 18 * bdim + 9 * idx + 8); 

            auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12); 

            auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13); 

            auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14); 

            auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15); 

            auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16); 

            auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17); 

            auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tdx_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 12); 

            auto tdy_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tdz_z_xx_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tdx_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 13); 

            auto tdy_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tdz_z_xy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tdx_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 14); 

            auto tdy_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tdz_z_xz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tdx_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 15); 

            auto tdy_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tdz_z_yy_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tdx_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 16); 

            auto tdy_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tdz_z_yz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tdx_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * idx + 17); 

            auto tdy_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tdz_z_zz_0 = primBuffer.data(pidx_d_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            // set up pointers to integrals

            auto tlx_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 24); 

            auto tly_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tlz_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tlx_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 25); 

            auto tly_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tlz_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tlx_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 26); 

            auto tly_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tlz_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tlx_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 27); 

            auto tly_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tlz_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tlx_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 28); 

            auto tly_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tlz_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tlx_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 29); 

            auto tly_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tlz_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30); 

            auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31); 

            auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32); 

            auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tlx_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 33); 

            auto tly_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tlz_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tlx_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 34); 

            auto tly_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tlz_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tlx_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 35); 

            auto tly_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tlz_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            // Batch of Integrals (72,108)

            #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_z_xx_0, tdx_z_xy_0, tdx_z_xz_0, tdx_z_yy_0, \
                                     tdx_z_yz_0, tdx_z_zz_0, tdy_z_xx_0, tdy_z_xy_0, tdy_z_xz_0, tdy_z_yy_0, tdy_z_yz_0, \
                                     tdy_z_zz_0, tdz_z_xx_0, tdz_z_xy_0, tdz_z_xz_0, tdz_z_yy_0, tdz_z_yz_0, tdz_z_zz_0, \
                                     tlx_0_xx_0, tlx_0_xy_0, tlx_0_xz_0, tlx_0_yy_0, tlx_0_yz_0, tlx_0_zz_0, \
                                     tlx_yz_xx_0, tlx_yz_xy_0, tlx_yz_xz_0, tlx_yz_yy_0, tlx_yz_yz_0, tlx_yz_zz_0, \
                                     tlx_z_x_0, tlx_z_xx_0, tlx_z_xy_0, tlx_z_xz_0, tlx_z_y_0, tlx_z_yy_0, tlx_z_yz_0, \
                                     tlx_z_z_0, tlx_z_zz_0, tlx_zz_xx_0, tlx_zz_xy_0, tlx_zz_xz_0, tlx_zz_yy_0, \
                                     tlx_zz_yz_0, tlx_zz_zz_0, tly_0_xx_0, tly_0_xy_0, tly_0_xz_0, tly_0_yy_0, \
                                     tly_0_yz_0, tly_0_zz_0, tly_yz_xx_0, tly_yz_xy_0, tly_yz_xz_0, tly_yz_yy_0, \
                                     tly_yz_yz_0, tly_yz_zz_0, tly_z_x_0, tly_z_xx_0, tly_z_xy_0, tly_z_xz_0, tly_z_y_0, \
                                     tly_z_yy_0, tly_z_yz_0, tly_z_z_0, tly_z_zz_0, tly_zz_xx_0, tly_zz_xy_0, \
                                     tly_zz_xz_0, tly_zz_yy_0, tly_zz_yz_0, tly_zz_zz_0, tlz_0_xx_0, tlz_0_xy_0, \
                                     tlz_0_xz_0, tlz_0_yy_0, tlz_0_yz_0, tlz_0_zz_0, tlz_yz_xx_0, tlz_yz_xy_0, \
                                     tlz_yz_xz_0, tlz_yz_yy_0, tlz_yz_yz_0, tlz_yz_zz_0, tlz_z_x_0, tlz_z_xx_0, \
                                     tlz_z_xy_0, tlz_z_xz_0, tlz_z_y_0, tlz_z_yy_0, tlz_z_yz_0, tlz_z_z_0, tlz_z_zz_0, \
                                     tlz_zz_xx_0, tlz_zz_xy_0, tlz_zz_xz_0, tlz_zz_yy_0, tlz_zz_yz_0, tlz_zz_zz_0, \
                                     tpx_z_xx_0, tpx_z_xy_0, tpx_z_xz_0, tpx_z_yy_0, tpx_z_yz_0, tpx_z_zz_0, tpy_z_xx_0, \
                                     tpy_z_xy_0, tpy_z_xz_0, tpy_z_yy_0, tpy_z_yz_0, tpy_z_zz_0, tpz_z_xx_0, tpz_z_xy_0, \
                                     tpz_z_xz_0, tpz_z_yy_0, tpz_z_yz_0, tpz_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yz_xx_0[j] = pa_y[j] * tlx_z_xx_0[j] - 0.5 * fl1_fx * tpz_z_xx_0[j] - fl1_fx * fl1_fgb * tdz_z_xx_0[j];

                tly_yz_xx_0[j] = pa_y[j] * tly_z_xx_0[j];

                tlz_yz_xx_0[j] = pa_y[j] * tlz_z_xx_0[j] + 0.5 * fl1_fx * tpx_z_xx_0[j] + fl1_fx * fl1_fgb * tdx_z_xx_0[j];

                tlx_yz_xy_0[j] = pa_y[j] * tlx_z_xy_0[j] + 0.5 * fl1_fx * tlx_z_x_0[j] - 0.5 * fl1_fx * tpz_z_xy_0[j] - fl1_fx * fl1_fgb * tdz_z_xy_0[j];

                tly_yz_xy_0[j] = pa_y[j] * tly_z_xy_0[j] + 0.5 * fl1_fx * tly_z_x_0[j];

                tlz_yz_xy_0[j] = pa_y[j] * tlz_z_xy_0[j] + 0.5 * fl1_fx * tlz_z_x_0[j] + 0.5 * fl1_fx * tpx_z_xy_0[j] + fl1_fx * fl1_fgb * tdx_z_xy_0[j];

                tlx_yz_xz_0[j] = pa_y[j] * tlx_z_xz_0[j] - 0.5 * fl1_fx * tpz_z_xz_0[j] - fl1_fx * fl1_fgb * tdz_z_xz_0[j];

                tly_yz_xz_0[j] = pa_y[j] * tly_z_xz_0[j];

                tlz_yz_xz_0[j] = pa_y[j] * tlz_z_xz_0[j] + 0.5 * fl1_fx * tpx_z_xz_0[j] + fl1_fx * fl1_fgb * tdx_z_xz_0[j];

                tlx_yz_yy_0[j] = pa_y[j] * tlx_z_yy_0[j] + fl1_fx * tlx_z_y_0[j] - 0.5 * fl1_fx * tpz_z_yy_0[j] - fl1_fx * fl1_fgb * tdz_z_yy_0[j];

                tly_yz_yy_0[j] = pa_y[j] * tly_z_yy_0[j] + fl1_fx * tly_z_y_0[j];

                tlz_yz_yy_0[j] = pa_y[j] * tlz_z_yy_0[j] + fl1_fx * tlz_z_y_0[j] + 0.5 * fl1_fx * tpx_z_yy_0[j] + fl1_fx * fl1_fgb * tdx_z_yy_0[j];

                tlx_yz_yz_0[j] = pa_y[j] * tlx_z_yz_0[j] + 0.5 * fl1_fx * tlx_z_z_0[j] - 0.5 * fl1_fx * tpz_z_yz_0[j] - fl1_fx * fl1_fgb * tdz_z_yz_0[j];

                tly_yz_yz_0[j] = pa_y[j] * tly_z_yz_0[j] + 0.5 * fl1_fx * tly_z_z_0[j];

                tlz_yz_yz_0[j] = pa_y[j] * tlz_z_yz_0[j] + 0.5 * fl1_fx * tlz_z_z_0[j] + 0.5 * fl1_fx * tpx_z_yz_0[j] + fl1_fx * fl1_fgb * tdx_z_yz_0[j];

                tlx_yz_zz_0[j] = pa_y[j] * tlx_z_zz_0[j] - 0.5 * fl1_fx * tpz_z_zz_0[j] - fl1_fx * fl1_fgb * tdz_z_zz_0[j];

                tly_yz_zz_0[j] = pa_y[j] * tly_z_zz_0[j];

                tlz_yz_zz_0[j] = pa_y[j] * tlz_z_zz_0[j] + 0.5 * fl1_fx * tpx_z_zz_0[j] + fl1_fx * fl1_fgb * tdx_z_zz_0[j];

                tlx_zz_xx_0[j] = pa_z[j] * tlx_z_xx_0[j] + 0.5 * fl1_fx * tlx_0_xx_0[j] + 0.5 * fl1_fx * tpy_z_xx_0[j] + fl1_fx * fl1_fgb * tdy_z_xx_0[j];

                tly_zz_xx_0[j] = pa_z[j] * tly_z_xx_0[j] + 0.5 * fl1_fx * tly_0_xx_0[j] - 0.5 * fl1_fx * tpx_z_xx_0[j] - fl1_fx * fl1_fgb * tdx_z_xx_0[j];

                tlz_zz_xx_0[j] = pa_z[j] * tlz_z_xx_0[j] + 0.5 * fl1_fx * tlz_0_xx_0[j];

                tlx_zz_xy_0[j] = pa_z[j] * tlx_z_xy_0[j] + 0.5 * fl1_fx * tlx_0_xy_0[j] + 0.5 * fl1_fx * tpy_z_xy_0[j] + fl1_fx * fl1_fgb * tdy_z_xy_0[j];

                tly_zz_xy_0[j] = pa_z[j] * tly_z_xy_0[j] + 0.5 * fl1_fx * tly_0_xy_0[j] - 0.5 * fl1_fx * tpx_z_xy_0[j] - fl1_fx * fl1_fgb * tdx_z_xy_0[j];

                tlz_zz_xy_0[j] = pa_z[j] * tlz_z_xy_0[j] + 0.5 * fl1_fx * tlz_0_xy_0[j];

                tlx_zz_xz_0[j] = pa_z[j] * tlx_z_xz_0[j] + 0.5 * fl1_fx * tlx_0_xz_0[j] + 0.5 * fl1_fx * tlx_z_x_0[j] + 0.5 * fl1_fx * tpy_z_xz_0[j] + fl1_fx * fl1_fgb * tdy_z_xz_0[j];

                tly_zz_xz_0[j] = pa_z[j] * tly_z_xz_0[j] + 0.5 * fl1_fx * tly_0_xz_0[j] + 0.5 * fl1_fx * tly_z_x_0[j] - 0.5 * fl1_fx * tpx_z_xz_0[j] - fl1_fx * fl1_fgb * tdx_z_xz_0[j];

                tlz_zz_xz_0[j] = pa_z[j] * tlz_z_xz_0[j] + 0.5 * fl1_fx * tlz_0_xz_0[j] + 0.5 * fl1_fx * tlz_z_x_0[j];

                tlx_zz_yy_0[j] = pa_z[j] * tlx_z_yy_0[j] + 0.5 * fl1_fx * tlx_0_yy_0[j] + 0.5 * fl1_fx * tpy_z_yy_0[j] + fl1_fx * fl1_fgb * tdy_z_yy_0[j];

                tly_zz_yy_0[j] = pa_z[j] * tly_z_yy_0[j] + 0.5 * fl1_fx * tly_0_yy_0[j] - 0.5 * fl1_fx * tpx_z_yy_0[j] - fl1_fx * fl1_fgb * tdx_z_yy_0[j];

                tlz_zz_yy_0[j] = pa_z[j] * tlz_z_yy_0[j] + 0.5 * fl1_fx * tlz_0_yy_0[j];

                tlx_zz_yz_0[j] = pa_z[j] * tlx_z_yz_0[j] + 0.5 * fl1_fx * tlx_0_yz_0[j] + 0.5 * fl1_fx * tlx_z_y_0[j] + 0.5 * fl1_fx * tpy_z_yz_0[j] + fl1_fx * fl1_fgb * tdy_z_yz_0[j];

                tly_zz_yz_0[j] = pa_z[j] * tly_z_yz_0[j] + 0.5 * fl1_fx * tly_0_yz_0[j] + 0.5 * fl1_fx * tly_z_y_0[j] - 0.5 * fl1_fx * tpx_z_yz_0[j] - fl1_fx * fl1_fgb * tdx_z_yz_0[j];

                tlz_zz_yz_0[j] = pa_z[j] * tlz_z_yz_0[j] + 0.5 * fl1_fx * tlz_0_yz_0[j] + 0.5 * fl1_fx * tlz_z_y_0[j];

                tlx_zz_zz_0[j] = pa_z[j] * tlx_z_zz_0[j] + 0.5 * fl1_fx * tlx_0_zz_0[j] + fl1_fx * tlx_z_z_0[j] + 0.5 * fl1_fx * tpy_z_zz_0[j] + fl1_fx * fl1_fgb * tdy_z_zz_0[j];

                tly_zz_zz_0[j] = pa_z[j] * tly_z_zz_0[j] + 0.5 * fl1_fx * tly_0_zz_0[j] + fl1_fx * tly_z_z_0[j] - 0.5 * fl1_fx * tpx_z_zz_0[j] - fl1_fx * fl1_fgb * tdx_z_zz_0[j];

                tlz_zz_zz_0[j] = pa_z[j] * tlz_z_zz_0[j] + 0.5 * fl1_fx * tlz_0_zz_0[j] + fl1_fx * tlz_z_z_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDF(      CMemBlock2D<double>& primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        amomrecfunc::compAngularMomentumForDF_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        amomrecfunc::compAngularMomentumForDF_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        amomrecfunc::compAngularMomentumForDF_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        amomrecfunc::compAngularMomentumForDF_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compAngularMomentumForDF_0_45(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx); 

            auto tly_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx); 

            auto tlz_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx); 

            auto tlx_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 1); 

            auto tly_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 1); 

            auto tlz_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 1); 

            auto tlx_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 2); 

            auto tly_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 2); 

            auto tlz_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 2); 

            auto tlx_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 3); 

            auto tly_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 3); 

            auto tlz_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 3); 

            auto tlx_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 4); 

            auto tly_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 4); 

            auto tlz_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 4); 

            auto tlx_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 5); 

            auto tly_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 5); 

            auto tlz_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 5); 

            auto tlx_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 6); 

            auto tly_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 6); 

            auto tlz_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 6); 

            auto tlx_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 7); 

            auto tly_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 7); 

            auto tlz_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 7); 

            auto tlx_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 8); 

            auto tly_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 8); 

            auto tlz_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 8); 

            auto tlx_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 9); 

            auto tly_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 9); 

            auto tlz_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 9); 

            auto tlx_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 10); 

            auto tly_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tlz_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tlx_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 11); 

            auto tly_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tlz_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tlx_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 12); 

            auto tly_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tlz_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tlx_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 13); 

            auto tly_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tlz_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tlx_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 14); 

            auto tly_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tlz_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx); 

            auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx); 

            auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx); 

            auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1); 

            auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1); 

            auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1); 

            auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2); 

            auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2); 

            auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2); 

            auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3); 

            auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3); 

            auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3); 

            auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4); 

            auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4); 

            auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4); 

            auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5); 

            auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5); 

            auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5); 

            auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6); 

            auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6); 

            auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6); 

            auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7); 

            auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7); 

            auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7); 

            auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8); 

            auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8); 

            auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8); 

            auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9); 

            auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9); 

            auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9); 

            auto tlx_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx); 

            auto tly_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tlz_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tlx_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 1); 

            auto tly_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tlz_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tlx_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 2); 

            auto tly_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tlz_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tlx_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 3); 

            auto tly_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tlz_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tlx_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 4); 

            auto tly_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tlz_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tlx_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 5); 

            auto tly_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tlz_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tlx_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 6); 

            auto tly_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tlz_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tlx_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 7); 

            auto tly_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tlz_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tlx_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 8); 

            auto tly_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tlz_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tlx_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 9); 

            auto tly_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 10); 

            auto tly_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tpy_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx); 

            auto tpz_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx); 

            auto tpy_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 1); 

            auto tpz_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 1); 

            auto tpy_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 2); 

            auto tpz_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 2); 

            auto tpy_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 3); 

            auto tpz_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 3); 

            auto tpy_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 4); 

            auto tpz_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 4); 

            auto tpy_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 5); 

            auto tpz_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 5); 

            auto tpy_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 6); 

            auto tpz_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 6); 

            auto tpy_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 7); 

            auto tpz_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 7); 

            auto tpy_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 8); 

            auto tpz_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 8); 

            auto tpy_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 9); 

            auto tpz_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 9); 

            auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tdy_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx); 

            auto tdz_x_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx); 

            auto tdy_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 1); 

            auto tdz_x_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 1); 

            auto tdy_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 2); 

            auto tdz_x_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 2); 

            auto tdy_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 3); 

            auto tdz_x_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 3); 

            auto tdy_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 4); 

            auto tdz_x_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 4); 

            auto tdy_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 5); 

            auto tdz_x_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 5); 

            auto tdy_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 6); 

            auto tdz_x_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 6); 

            auto tdy_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 7); 

            auto tdz_x_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 7); 

            auto tdy_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 8); 

            auto tdz_x_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 8); 

            auto tdy_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 9); 

            auto tdz_x_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 9); 

            auto tdy_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tdy_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tdy_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tdy_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tdy_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            // set up pointers to integrals

            auto tlx_xx_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx); 

            auto tly_xx_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx); 

            auto tlz_xx_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx); 

            auto tlx_xx_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 1); 

            auto tly_xx_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 1); 

            auto tlz_xx_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 1); 

            auto tlx_xx_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 2); 

            auto tly_xx_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 2); 

            auto tlz_xx_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 2); 

            auto tlx_xx_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 3); 

            auto tly_xx_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 3); 

            auto tlz_xx_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 3); 

            auto tlx_xx_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 4); 

            auto tly_xx_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 4); 

            auto tlz_xx_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 4); 

            auto tlx_xx_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 5); 

            auto tly_xx_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 5); 

            auto tlz_xx_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 5); 

            auto tlx_xx_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 6); 

            auto tly_xx_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 6); 

            auto tlz_xx_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 6); 

            auto tlx_xx_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 7); 

            auto tly_xx_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 7); 

            auto tlz_xx_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 7); 

            auto tlx_xx_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 8); 

            auto tly_xx_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 8); 

            auto tlz_xx_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 8); 

            auto tlx_xx_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 9); 

            auto tly_xx_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 9); 

            auto tlz_xx_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 9); 

            auto tlx_xy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 10); 

            auto tly_xy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 10); 

            auto tlz_xy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 10); 

            auto tlx_xy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 11); 

            auto tly_xy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 11); 

            auto tlz_xy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 11); 

            auto tlx_xy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 12); 

            auto tly_xy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 12); 

            auto tlz_xy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 12); 

            auto tlx_xy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 13); 

            auto tly_xy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 13); 

            auto tlz_xy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 13); 

            auto tlx_xy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 14); 

            auto tly_xy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 14); 

            auto tlz_xy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_x_xxx_0, tdy_x_xxy_0, tdy_x_xxz_0, tdy_x_xyy_0, \
                                     tdy_x_xyz_0, tdy_x_xzz_0, tdy_x_yyy_0, tdy_x_yyz_0, tdy_x_yzz_0, tdy_x_zzz_0, \
                                     tdy_y_xxx_0, tdy_y_xxy_0, tdy_y_xxz_0, tdy_y_xyy_0, tdy_y_xyz_0, tdz_x_xxx_0, \
                                     tdz_x_xxy_0, tdz_x_xxz_0, tdz_x_xyy_0, tdz_x_xyz_0, tdz_x_xzz_0, tdz_x_yyy_0, \
                                     tdz_x_yyz_0, tdz_x_yzz_0, tdz_x_zzz_0, tdz_y_xxx_0, tdz_y_xxy_0, tdz_y_xxz_0, \
                                     tdz_y_xyy_0, tdz_y_xyz_0, tlx_0_xxx_0, tlx_0_xxy_0, tlx_0_xxz_0, tlx_0_xyy_0, \
                                     tlx_0_xyz_0, tlx_0_xzz_0, tlx_0_yyy_0, tlx_0_yyz_0, tlx_0_yzz_0, tlx_0_zzz_0, \
                                     tlx_x_xx_0, tlx_x_xxx_0, tlx_x_xxy_0, tlx_x_xxz_0, tlx_x_xy_0, tlx_x_xyy_0, \
                                     tlx_x_xyz_0, tlx_x_xz_0, tlx_x_xzz_0, tlx_x_yy_0, tlx_x_yyy_0, tlx_x_yyz_0, \
                                     tlx_x_yz_0, tlx_x_yzz_0, tlx_x_zz_0, tlx_x_zzz_0, tlx_xx_xxx_0, tlx_xx_xxy_0, \
                                     tlx_xx_xxz_0, tlx_xx_xyy_0, tlx_xx_xyz_0, tlx_xx_xzz_0, tlx_xx_yyy_0, tlx_xx_yyz_0, \
                                     tlx_xx_yzz_0, tlx_xx_zzz_0, tlx_xy_xxx_0, tlx_xy_xxy_0, tlx_xy_xxz_0, tlx_xy_xyy_0, \
                                     tlx_xy_xyz_0, tlx_y_xx_0, tlx_y_xxx_0, tlx_y_xxy_0, tlx_y_xxz_0, tlx_y_xy_0, \
                                     tlx_y_xyy_0, tlx_y_xyz_0, tlx_y_xz_0, tlx_y_yy_0, tlx_y_yz_0, tly_0_xxx_0, \
                                     tly_0_xxy_0, tly_0_xxz_0, tly_0_xyy_0, tly_0_xyz_0, tly_0_xzz_0, tly_0_yyy_0, \
                                     tly_0_yyz_0, tly_0_yzz_0, tly_0_zzz_0, tly_x_xx_0, tly_x_xxx_0, tly_x_xxy_0, \
                                     tly_x_xxz_0, tly_x_xy_0, tly_x_xyy_0, tly_x_xyz_0, tly_x_xz_0, tly_x_xzz_0, \
                                     tly_x_yy_0, tly_x_yyy_0, tly_x_yyz_0, tly_x_yz_0, tly_x_yzz_0, tly_x_zz_0, \
                                     tly_x_zzz_0, tly_xx_xxx_0, tly_xx_xxy_0, tly_xx_xxz_0, tly_xx_xyy_0, tly_xx_xyz_0, \
                                     tly_xx_xzz_0, tly_xx_yyy_0, tly_xx_yyz_0, tly_xx_yzz_0, tly_xx_zzz_0, tly_xy_xxx_0, \
                                     tly_xy_xxy_0, tly_xy_xxz_0, tly_xy_xyy_0, tly_xy_xyz_0, tly_y_xx_0, tly_y_xxx_0, \
                                     tly_y_xxy_0, tly_y_xxz_0, tly_y_xy_0, tly_y_xyy_0, tly_y_xyz_0, tly_y_xz_0, \
                                     tly_y_yy_0, tly_y_yz_0, tlz_0_xxx_0, tlz_0_xxy_0, tlz_0_xxz_0, tlz_0_xyy_0, \
                                     tlz_0_xyz_0, tlz_0_xzz_0, tlz_0_yyy_0, tlz_0_yyz_0, tlz_0_yzz_0, tlz_0_zzz_0, \
                                     tlz_x_xx_0, tlz_x_xxx_0, tlz_x_xxy_0, tlz_x_xxz_0, tlz_x_xy_0, tlz_x_xyy_0, \
                                     tlz_x_xyz_0, tlz_x_xz_0, tlz_x_xzz_0, tlz_x_yy_0, tlz_x_yyy_0, tlz_x_yyz_0, \
                                     tlz_x_yz_0, tlz_x_yzz_0, tlz_x_zz_0, tlz_x_zzz_0, tlz_xx_xxx_0, tlz_xx_xxy_0, \
                                     tlz_xx_xxz_0, tlz_xx_xyy_0, tlz_xx_xyz_0, tlz_xx_xzz_0, tlz_xx_yyy_0, tlz_xx_yyz_0, \
                                     tlz_xx_yzz_0, tlz_xx_zzz_0, tlz_xy_xxx_0, tlz_xy_xxy_0, tlz_xy_xxz_0, tlz_xy_xyy_0, \
                                     tlz_xy_xyz_0, tlz_y_xx_0, tlz_y_xxx_0, tlz_y_xxy_0, tlz_y_xxz_0, tlz_y_xy_0, \
                                     tlz_y_xyy_0, tlz_y_xyz_0, tlz_y_xz_0, tlz_y_yy_0, tlz_y_yz_0, tpy_x_xxx_0, \
                                     tpy_x_xxy_0, tpy_x_xxz_0, tpy_x_xyy_0, tpy_x_xyz_0, tpy_x_xzz_0, tpy_x_yyy_0, \
                                     tpy_x_yyz_0, tpy_x_yzz_0, tpy_x_zzz_0, tpy_y_xxx_0, tpy_y_xxy_0, tpy_y_xxz_0, \
                                     tpy_y_xyy_0, tpy_y_xyz_0, tpz_x_xxx_0, tpz_x_xxy_0, tpz_x_xxz_0, tpz_x_xyy_0, \
                                     tpz_x_xyz_0, tpz_x_xzz_0, tpz_x_yyy_0, tpz_x_yyz_0, tpz_x_yzz_0, tpz_x_zzz_0, \
                                     tpz_y_xxx_0, tpz_y_xxy_0, tpz_y_xxz_0, tpz_y_xyy_0, tpz_y_xyz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xx_xxx_0[j] = pa_x[j] * tlx_x_xxx_0[j] + 0.5 * fl1_fx * tlx_0_xxx_0[j] + 1.5 * fl1_fx * tlx_x_xx_0[j];

                tly_xx_xxx_0[j] = pa_x[j] * tly_x_xxx_0[j] + 0.5 * fl1_fx * tly_0_xxx_0[j] + 1.5 * fl1_fx * tly_x_xx_0[j] + 0.5 * fl1_fx * tpz_x_xxx_0[j] + fl1_fx * fl1_fgb * tdz_x_xxx_0[j];

                tlz_xx_xxx_0[j] = pa_x[j] * tlz_x_xxx_0[j] + 0.5 * fl1_fx * tlz_0_xxx_0[j] + 1.5 * fl1_fx * tlz_x_xx_0[j] - 0.5 * fl1_fx * tpy_x_xxx_0[j] - fl1_fx * fl1_fgb * tdy_x_xxx_0[j];

                tlx_xx_xxy_0[j] = pa_x[j] * tlx_x_xxy_0[j] + 0.5 * fl1_fx * tlx_0_xxy_0[j] + fl1_fx * tlx_x_xy_0[j];

                tly_xx_xxy_0[j] = pa_x[j] * tly_x_xxy_0[j] + 0.5 * fl1_fx * tly_0_xxy_0[j] + fl1_fx * tly_x_xy_0[j] + 0.5 * fl1_fx * tpz_x_xxy_0[j] + fl1_fx * fl1_fgb * tdz_x_xxy_0[j];

                tlz_xx_xxy_0[j] = pa_x[j] * tlz_x_xxy_0[j] + 0.5 * fl1_fx * tlz_0_xxy_0[j] + fl1_fx * tlz_x_xy_0[j] - 0.5 * fl1_fx * tpy_x_xxy_0[j] - fl1_fx * fl1_fgb * tdy_x_xxy_0[j];

                tlx_xx_xxz_0[j] = pa_x[j] * tlx_x_xxz_0[j] + 0.5 * fl1_fx * tlx_0_xxz_0[j] + fl1_fx * tlx_x_xz_0[j];

                tly_xx_xxz_0[j] = pa_x[j] * tly_x_xxz_0[j] + 0.5 * fl1_fx * tly_0_xxz_0[j] + fl1_fx * tly_x_xz_0[j] + 0.5 * fl1_fx * tpz_x_xxz_0[j] + fl1_fx * fl1_fgb * tdz_x_xxz_0[j];

                tlz_xx_xxz_0[j] = pa_x[j] * tlz_x_xxz_0[j] + 0.5 * fl1_fx * tlz_0_xxz_0[j] + fl1_fx * tlz_x_xz_0[j] - 0.5 * fl1_fx * tpy_x_xxz_0[j] - fl1_fx * fl1_fgb * tdy_x_xxz_0[j];

                tlx_xx_xyy_0[j] = pa_x[j] * tlx_x_xyy_0[j] + 0.5 * fl1_fx * tlx_0_xyy_0[j] + 0.5 * fl1_fx * tlx_x_yy_0[j];

                tly_xx_xyy_0[j] = pa_x[j] * tly_x_xyy_0[j] + 0.5 * fl1_fx * tly_0_xyy_0[j] + 0.5 * fl1_fx * tly_x_yy_0[j] + 0.5 * fl1_fx * tpz_x_xyy_0[j] + fl1_fx * fl1_fgb * tdz_x_xyy_0[j];

                tlz_xx_xyy_0[j] = pa_x[j] * tlz_x_xyy_0[j] + 0.5 * fl1_fx * tlz_0_xyy_0[j] + 0.5 * fl1_fx * tlz_x_yy_0[j] - 0.5 * fl1_fx * tpy_x_xyy_0[j] - fl1_fx * fl1_fgb * tdy_x_xyy_0[j];

                tlx_xx_xyz_0[j] = pa_x[j] * tlx_x_xyz_0[j] + 0.5 * fl1_fx * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tlx_x_yz_0[j];

                tly_xx_xyz_0[j] = pa_x[j] * tly_x_xyz_0[j] + 0.5 * fl1_fx * tly_0_xyz_0[j] + 0.5 * fl1_fx * tly_x_yz_0[j] + 0.5 * fl1_fx * tpz_x_xyz_0[j] + fl1_fx * fl1_fgb * tdz_x_xyz_0[j];

                tlz_xx_xyz_0[j] = pa_x[j] * tlz_x_xyz_0[j] + 0.5 * fl1_fx * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tlz_x_yz_0[j] - 0.5 * fl1_fx * tpy_x_xyz_0[j] - fl1_fx * fl1_fgb * tdy_x_xyz_0[j];

                tlx_xx_xzz_0[j] = pa_x[j] * tlx_x_xzz_0[j] + 0.5 * fl1_fx * tlx_0_xzz_0[j] + 0.5 * fl1_fx * tlx_x_zz_0[j];

                tly_xx_xzz_0[j] = pa_x[j] * tly_x_xzz_0[j] + 0.5 * fl1_fx * tly_0_xzz_0[j] + 0.5 * fl1_fx * tly_x_zz_0[j] + 0.5 * fl1_fx * tpz_x_xzz_0[j] + fl1_fx * fl1_fgb * tdz_x_xzz_0[j];

                tlz_xx_xzz_0[j] = pa_x[j] * tlz_x_xzz_0[j] + 0.5 * fl1_fx * tlz_0_xzz_0[j] + 0.5 * fl1_fx * tlz_x_zz_0[j] - 0.5 * fl1_fx * tpy_x_xzz_0[j] - fl1_fx * fl1_fgb * tdy_x_xzz_0[j];

                tlx_xx_yyy_0[j] = pa_x[j] * tlx_x_yyy_0[j] + 0.5 * fl1_fx * tlx_0_yyy_0[j];

                tly_xx_yyy_0[j] = pa_x[j] * tly_x_yyy_0[j] + 0.5 * fl1_fx * tly_0_yyy_0[j] + 0.5 * fl1_fx * tpz_x_yyy_0[j] + fl1_fx * fl1_fgb * tdz_x_yyy_0[j];

                tlz_xx_yyy_0[j] = pa_x[j] * tlz_x_yyy_0[j] + 0.5 * fl1_fx * tlz_0_yyy_0[j] - 0.5 * fl1_fx * tpy_x_yyy_0[j] - fl1_fx * fl1_fgb * tdy_x_yyy_0[j];

                tlx_xx_yyz_0[j] = pa_x[j] * tlx_x_yyz_0[j] + 0.5 * fl1_fx * tlx_0_yyz_0[j];

                tly_xx_yyz_0[j] = pa_x[j] * tly_x_yyz_0[j] + 0.5 * fl1_fx * tly_0_yyz_0[j] + 0.5 * fl1_fx * tpz_x_yyz_0[j] + fl1_fx * fl1_fgb * tdz_x_yyz_0[j];

                tlz_xx_yyz_0[j] = pa_x[j] * tlz_x_yyz_0[j] + 0.5 * fl1_fx * tlz_0_yyz_0[j] - 0.5 * fl1_fx * tpy_x_yyz_0[j] - fl1_fx * fl1_fgb * tdy_x_yyz_0[j];

                tlx_xx_yzz_0[j] = pa_x[j] * tlx_x_yzz_0[j] + 0.5 * fl1_fx * tlx_0_yzz_0[j];

                tly_xx_yzz_0[j] = pa_x[j] * tly_x_yzz_0[j] + 0.5 * fl1_fx * tly_0_yzz_0[j] + 0.5 * fl1_fx * tpz_x_yzz_0[j] + fl1_fx * fl1_fgb * tdz_x_yzz_0[j];

                tlz_xx_yzz_0[j] = pa_x[j] * tlz_x_yzz_0[j] + 0.5 * fl1_fx * tlz_0_yzz_0[j] - 0.5 * fl1_fx * tpy_x_yzz_0[j] - fl1_fx * fl1_fgb * tdy_x_yzz_0[j];

                tlx_xx_zzz_0[j] = pa_x[j] * tlx_x_zzz_0[j] + 0.5 * fl1_fx * tlx_0_zzz_0[j];

                tly_xx_zzz_0[j] = pa_x[j] * tly_x_zzz_0[j] + 0.5 * fl1_fx * tly_0_zzz_0[j] + 0.5 * fl1_fx * tpz_x_zzz_0[j] + fl1_fx * fl1_fgb * tdz_x_zzz_0[j];

                tlz_xx_zzz_0[j] = pa_x[j] * tlz_x_zzz_0[j] + 0.5 * fl1_fx * tlz_0_zzz_0[j] - 0.5 * fl1_fx * tpy_x_zzz_0[j] - fl1_fx * fl1_fgb * tdy_x_zzz_0[j];

                tlx_xy_xxx_0[j] = pa_x[j] * tlx_y_xxx_0[j] + 1.5 * fl1_fx * tlx_y_xx_0[j];

                tly_xy_xxx_0[j] = pa_x[j] * tly_y_xxx_0[j] + 1.5 * fl1_fx * tly_y_xx_0[j] + 0.5 * fl1_fx * tpz_y_xxx_0[j] + fl1_fx * fl1_fgb * tdz_y_xxx_0[j];

                tlz_xy_xxx_0[j] = pa_x[j] * tlz_y_xxx_0[j] + 1.5 * fl1_fx * tlz_y_xx_0[j] - 0.5 * fl1_fx * tpy_y_xxx_0[j] - fl1_fx * fl1_fgb * tdy_y_xxx_0[j];

                tlx_xy_xxy_0[j] = pa_x[j] * tlx_y_xxy_0[j] + fl1_fx * tlx_y_xy_0[j];

                tly_xy_xxy_0[j] = pa_x[j] * tly_y_xxy_0[j] + fl1_fx * tly_y_xy_0[j] + 0.5 * fl1_fx * tpz_y_xxy_0[j] + fl1_fx * fl1_fgb * tdz_y_xxy_0[j];

                tlz_xy_xxy_0[j] = pa_x[j] * tlz_y_xxy_0[j] + fl1_fx * tlz_y_xy_0[j] - 0.5 * fl1_fx * tpy_y_xxy_0[j] - fl1_fx * fl1_fgb * tdy_y_xxy_0[j];

                tlx_xy_xxz_0[j] = pa_x[j] * tlx_y_xxz_0[j] + fl1_fx * tlx_y_xz_0[j];

                tly_xy_xxz_0[j] = pa_x[j] * tly_y_xxz_0[j] + fl1_fx * tly_y_xz_0[j] + 0.5 * fl1_fx * tpz_y_xxz_0[j] + fl1_fx * fl1_fgb * tdz_y_xxz_0[j];

                tlz_xy_xxz_0[j] = pa_x[j] * tlz_y_xxz_0[j] + fl1_fx * tlz_y_xz_0[j] - 0.5 * fl1_fx * tpy_y_xxz_0[j] - fl1_fx * fl1_fgb * tdy_y_xxz_0[j];

                tlx_xy_xyy_0[j] = pa_x[j] * tlx_y_xyy_0[j] + 0.5 * fl1_fx * tlx_y_yy_0[j];

                tly_xy_xyy_0[j] = pa_x[j] * tly_y_xyy_0[j] + 0.5 * fl1_fx * tly_y_yy_0[j] + 0.5 * fl1_fx * tpz_y_xyy_0[j] + fl1_fx * fl1_fgb * tdz_y_xyy_0[j];

                tlz_xy_xyy_0[j] = pa_x[j] * tlz_y_xyy_0[j] + 0.5 * fl1_fx * tlz_y_yy_0[j] - 0.5 * fl1_fx * tpy_y_xyy_0[j] - fl1_fx * fl1_fgb * tdy_y_xyy_0[j];

                tlx_xy_xyz_0[j] = pa_x[j] * tlx_y_xyz_0[j] + 0.5 * fl1_fx * tlx_y_yz_0[j];

                tly_xy_xyz_0[j] = pa_x[j] * tly_y_xyz_0[j] + 0.5 * fl1_fx * tly_y_yz_0[j] + 0.5 * fl1_fx * tpz_y_xyz_0[j] + fl1_fx * fl1_fgb * tdz_y_xyz_0[j];

                tlz_xy_xyz_0[j] = pa_x[j] * tlz_y_xyz_0[j] + 0.5 * fl1_fx * tlz_y_yz_0[j] - 0.5 * fl1_fx * tpy_y_xyz_0[j] - fl1_fx * fl1_fgb * tdy_y_xyz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDF_45_90(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 15); 

            auto tly_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tlz_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tlx_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 16); 

            auto tly_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tlz_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tlx_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 17); 

            auto tly_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tlz_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tlx_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 18); 

            auto tly_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tlz_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tlx_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 19); 

            auto tly_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tlz_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tlx_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 20); 

            auto tly_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 21); 

            auto tly_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 22); 

            auto tly_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 23); 

            auto tly_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tlx_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 24); 

            auto tly_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tlx_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 25); 

            auto tly_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tlz_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tlx_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 26); 

            auto tly_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tlz_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tlx_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 27); 

            auto tly_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tlz_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tlx_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 28); 

            auto tly_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tlz_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tlx_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 29); 

            auto tly_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tlz_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            auto tlx_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 11); 

            auto tly_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 15); 

            auto tly_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 16); 

            auto tly_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 17); 

            auto tly_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            auto tdy_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tdy_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tdy_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tdy_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdy_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tdy_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tdz_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tdy_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tdz_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tdy_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tdz_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tdy_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tdz_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tdy_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tdz_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            // set up pointers to integrals

            auto tlx_xy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 15); 

            auto tly_xy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 15); 

            auto tlz_xy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 15); 

            auto tlx_xy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 16); 

            auto tly_xy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 16); 

            auto tlz_xy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 16); 

            auto tlx_xy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 17); 

            auto tly_xy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 17); 

            auto tlz_xy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 17); 

            auto tlx_xy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 18); 

            auto tly_xy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 18); 

            auto tlz_xy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 18); 

            auto tlx_xy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 19); 

            auto tly_xy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 19); 

            auto tlz_xy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 19); 

            auto tlx_xz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 20); 

            auto tly_xz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 20); 

            auto tlz_xz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 20); 

            auto tlx_xz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 21); 

            auto tly_xz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 21); 

            auto tlz_xz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 21); 

            auto tlx_xz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 22); 

            auto tly_xz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 22); 

            auto tlz_xz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 22); 

            auto tlx_xz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 23); 

            auto tly_xz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 23); 

            auto tlz_xz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 23); 

            auto tlx_xz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 24); 

            auto tly_xz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 24); 

            auto tlz_xz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 24); 

            auto tlx_xz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 25); 

            auto tly_xz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 25); 

            auto tlz_xz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 25); 

            auto tlx_xz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 26); 

            auto tly_xz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 26); 

            auto tlz_xz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 26); 

            auto tlx_xz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 27); 

            auto tly_xz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 27); 

            auto tlz_xz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 27); 

            auto tlx_xz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 28); 

            auto tly_xz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 28); 

            auto tlz_xz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 28); 

            auto tlx_xz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 29); 

            auto tly_xz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 29); 

            auto tlz_xz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_y_xzz_0, tdy_y_yyy_0, tdy_y_yyz_0, tdy_y_yzz_0, \
                                     tdy_y_zzz_0, tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, tdy_z_xyy_0, tdy_z_xyz_0, \
                                     tdy_z_xzz_0, tdy_z_yyy_0, tdy_z_yyz_0, tdy_z_yzz_0, tdy_z_zzz_0, tdz_y_xzz_0, \
                                     tdz_y_yyy_0, tdz_y_yyz_0, tdz_y_yzz_0, tdz_y_zzz_0, tdz_z_xxx_0, tdz_z_xxy_0, \
                                     tdz_z_xxz_0, tdz_z_xyy_0, tdz_z_xyz_0, tdz_z_xzz_0, tdz_z_yyy_0, tdz_z_yyz_0, \
                                     tdz_z_yzz_0, tdz_z_zzz_0, tlx_xy_xzz_0, tlx_xy_yyy_0, tlx_xy_yyz_0, tlx_xy_yzz_0, \
                                     tlx_xy_zzz_0, tlx_xz_xxx_0, tlx_xz_xxy_0, tlx_xz_xxz_0, tlx_xz_xyy_0, tlx_xz_xyz_0, \
                                     tlx_xz_xzz_0, tlx_xz_yyy_0, tlx_xz_yyz_0, tlx_xz_yzz_0, tlx_xz_zzz_0, tlx_y_xzz_0, \
                                     tlx_y_yyy_0, tlx_y_yyz_0, tlx_y_yzz_0, tlx_y_zz_0, tlx_y_zzz_0, tlx_z_xx_0, \
                                     tlx_z_xxx_0, tlx_z_xxy_0, tlx_z_xxz_0, tlx_z_xy_0, tlx_z_xyy_0, tlx_z_xyz_0, \
                                     tlx_z_xz_0, tlx_z_xzz_0, tlx_z_yy_0, tlx_z_yyy_0, tlx_z_yyz_0, tlx_z_yz_0, \
                                     tlx_z_yzz_0, tlx_z_zz_0, tlx_z_zzz_0, tly_xy_xzz_0, tly_xy_yyy_0, tly_xy_yyz_0, \
                                     tly_xy_yzz_0, tly_xy_zzz_0, tly_xz_xxx_0, tly_xz_xxy_0, tly_xz_xxz_0, tly_xz_xyy_0, \
                                     tly_xz_xyz_0, tly_xz_xzz_0, tly_xz_yyy_0, tly_xz_yyz_0, tly_xz_yzz_0, tly_xz_zzz_0, \
                                     tly_y_xzz_0, tly_y_yyy_0, tly_y_yyz_0, tly_y_yzz_0, tly_y_zz_0, tly_y_zzz_0, \
                                     tly_z_xx_0, tly_z_xxx_0, tly_z_xxy_0, tly_z_xxz_0, tly_z_xy_0, tly_z_xyy_0, \
                                     tly_z_xyz_0, tly_z_xz_0, tly_z_xzz_0, tly_z_yy_0, tly_z_yyy_0, tly_z_yyz_0, \
                                     tly_z_yz_0, tly_z_yzz_0, tly_z_zz_0, tly_z_zzz_0, tlz_xy_xzz_0, tlz_xy_yyy_0, \
                                     tlz_xy_yyz_0, tlz_xy_yzz_0, tlz_xy_zzz_0, tlz_xz_xxx_0, tlz_xz_xxy_0, tlz_xz_xxz_0, \
                                     tlz_xz_xyy_0, tlz_xz_xyz_0, tlz_xz_xzz_0, tlz_xz_yyy_0, tlz_xz_yyz_0, tlz_xz_yzz_0, \
                                     tlz_xz_zzz_0, tlz_y_xzz_0, tlz_y_yyy_0, tlz_y_yyz_0, tlz_y_yzz_0, tlz_y_zz_0, \
                                     tlz_y_zzz_0, tlz_z_xx_0, tlz_z_xxx_0, tlz_z_xxy_0, tlz_z_xxz_0, tlz_z_xy_0, \
                                     tlz_z_xyy_0, tlz_z_xyz_0, tlz_z_xz_0, tlz_z_xzz_0, tlz_z_yy_0, tlz_z_yyy_0, \
                                     tlz_z_yyz_0, tlz_z_yz_0, tlz_z_yzz_0, tlz_z_zz_0, tlz_z_zzz_0, tpy_y_xzz_0, \
                                     tpy_y_yyy_0, tpy_y_yyz_0, tpy_y_yzz_0, tpy_y_zzz_0, tpy_z_xxx_0, tpy_z_xxy_0, \
                                     tpy_z_xxz_0, tpy_z_xyy_0, tpy_z_xyz_0, tpy_z_xzz_0, tpy_z_yyy_0, tpy_z_yyz_0, \
                                     tpy_z_yzz_0, tpy_z_zzz_0, tpz_y_xzz_0, tpz_y_yyy_0, tpz_y_yyz_0, tpz_y_yzz_0, \
                                     tpz_y_zzz_0, tpz_z_xxx_0, tpz_z_xxy_0, tpz_z_xxz_0, tpz_z_xyy_0, tpz_z_xyz_0, \
                                     tpz_z_xzz_0, tpz_z_yyy_0, tpz_z_yyz_0, tpz_z_yzz_0, tpz_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xy_xzz_0[j] = pa_x[j] * tlx_y_xzz_0[j] + 0.5 * fl1_fx * tlx_y_zz_0[j];

                tly_xy_xzz_0[j] = pa_x[j] * tly_y_xzz_0[j] + 0.5 * fl1_fx * tly_y_zz_0[j] + 0.5 * fl1_fx * tpz_y_xzz_0[j] + fl1_fx * fl1_fgb * tdz_y_xzz_0[j];

                tlz_xy_xzz_0[j] = pa_x[j] * tlz_y_xzz_0[j] + 0.5 * fl1_fx * tlz_y_zz_0[j] - 0.5 * fl1_fx * tpy_y_xzz_0[j] - fl1_fx * fl1_fgb * tdy_y_xzz_0[j];

                tlx_xy_yyy_0[j] = pa_x[j] * tlx_y_yyy_0[j];

                tly_xy_yyy_0[j] = pa_x[j] * tly_y_yyy_0[j] + 0.5 * fl1_fx * tpz_y_yyy_0[j] + fl1_fx * fl1_fgb * tdz_y_yyy_0[j];

                tlz_xy_yyy_0[j] = pa_x[j] * tlz_y_yyy_0[j] - 0.5 * fl1_fx * tpy_y_yyy_0[j] - fl1_fx * fl1_fgb * tdy_y_yyy_0[j];

                tlx_xy_yyz_0[j] = pa_x[j] * tlx_y_yyz_0[j];

                tly_xy_yyz_0[j] = pa_x[j] * tly_y_yyz_0[j] + 0.5 * fl1_fx * tpz_y_yyz_0[j] + fl1_fx * fl1_fgb * tdz_y_yyz_0[j];

                tlz_xy_yyz_0[j] = pa_x[j] * tlz_y_yyz_0[j] - 0.5 * fl1_fx * tpy_y_yyz_0[j] - fl1_fx * fl1_fgb * tdy_y_yyz_0[j];

                tlx_xy_yzz_0[j] = pa_x[j] * tlx_y_yzz_0[j];

                tly_xy_yzz_0[j] = pa_x[j] * tly_y_yzz_0[j] + 0.5 * fl1_fx * tpz_y_yzz_0[j] + fl1_fx * fl1_fgb * tdz_y_yzz_0[j];

                tlz_xy_yzz_0[j] = pa_x[j] * tlz_y_yzz_0[j] - 0.5 * fl1_fx * tpy_y_yzz_0[j] - fl1_fx * fl1_fgb * tdy_y_yzz_0[j];

                tlx_xy_zzz_0[j] = pa_x[j] * tlx_y_zzz_0[j];

                tly_xy_zzz_0[j] = pa_x[j] * tly_y_zzz_0[j] + 0.5 * fl1_fx * tpz_y_zzz_0[j] + fl1_fx * fl1_fgb * tdz_y_zzz_0[j];

                tlz_xy_zzz_0[j] = pa_x[j] * tlz_y_zzz_0[j] - 0.5 * fl1_fx * tpy_y_zzz_0[j] - fl1_fx * fl1_fgb * tdy_y_zzz_0[j];

                tlx_xz_xxx_0[j] = pa_x[j] * tlx_z_xxx_0[j] + 1.5 * fl1_fx * tlx_z_xx_0[j];

                tly_xz_xxx_0[j] = pa_x[j] * tly_z_xxx_0[j] + 1.5 * fl1_fx * tly_z_xx_0[j] + 0.5 * fl1_fx * tpz_z_xxx_0[j] + fl1_fx * fl1_fgb * tdz_z_xxx_0[j];

                tlz_xz_xxx_0[j] = pa_x[j] * tlz_z_xxx_0[j] + 1.5 * fl1_fx * tlz_z_xx_0[j] - 0.5 * fl1_fx * tpy_z_xxx_0[j] - fl1_fx * fl1_fgb * tdy_z_xxx_0[j];

                tlx_xz_xxy_0[j] = pa_x[j] * tlx_z_xxy_0[j] + fl1_fx * tlx_z_xy_0[j];

                tly_xz_xxy_0[j] = pa_x[j] * tly_z_xxy_0[j] + fl1_fx * tly_z_xy_0[j] + 0.5 * fl1_fx * tpz_z_xxy_0[j] + fl1_fx * fl1_fgb * tdz_z_xxy_0[j];

                tlz_xz_xxy_0[j] = pa_x[j] * tlz_z_xxy_0[j] + fl1_fx * tlz_z_xy_0[j] - 0.5 * fl1_fx * tpy_z_xxy_0[j] - fl1_fx * fl1_fgb * tdy_z_xxy_0[j];

                tlx_xz_xxz_0[j] = pa_x[j] * tlx_z_xxz_0[j] + fl1_fx * tlx_z_xz_0[j];

                tly_xz_xxz_0[j] = pa_x[j] * tly_z_xxz_0[j] + fl1_fx * tly_z_xz_0[j] + 0.5 * fl1_fx * tpz_z_xxz_0[j] + fl1_fx * fl1_fgb * tdz_z_xxz_0[j];

                tlz_xz_xxz_0[j] = pa_x[j] * tlz_z_xxz_0[j] + fl1_fx * tlz_z_xz_0[j] - 0.5 * fl1_fx * tpy_z_xxz_0[j] - fl1_fx * fl1_fgb * tdy_z_xxz_0[j];

                tlx_xz_xyy_0[j] = pa_x[j] * tlx_z_xyy_0[j] + 0.5 * fl1_fx * tlx_z_yy_0[j];

                tly_xz_xyy_0[j] = pa_x[j] * tly_z_xyy_0[j] + 0.5 * fl1_fx * tly_z_yy_0[j] + 0.5 * fl1_fx * tpz_z_xyy_0[j] + fl1_fx * fl1_fgb * tdz_z_xyy_0[j];

                tlz_xz_xyy_0[j] = pa_x[j] * tlz_z_xyy_0[j] + 0.5 * fl1_fx * tlz_z_yy_0[j] - 0.5 * fl1_fx * tpy_z_xyy_0[j] - fl1_fx * fl1_fgb * tdy_z_xyy_0[j];

                tlx_xz_xyz_0[j] = pa_x[j] * tlx_z_xyz_0[j] + 0.5 * fl1_fx * tlx_z_yz_0[j];

                tly_xz_xyz_0[j] = pa_x[j] * tly_z_xyz_0[j] + 0.5 * fl1_fx * tly_z_yz_0[j] + 0.5 * fl1_fx * tpz_z_xyz_0[j] + fl1_fx * fl1_fgb * tdz_z_xyz_0[j];

                tlz_xz_xyz_0[j] = pa_x[j] * tlz_z_xyz_0[j] + 0.5 * fl1_fx * tlz_z_yz_0[j] - 0.5 * fl1_fx * tpy_z_xyz_0[j] - fl1_fx * fl1_fgb * tdy_z_xyz_0[j];

                tlx_xz_xzz_0[j] = pa_x[j] * tlx_z_xzz_0[j] + 0.5 * fl1_fx * tlx_z_zz_0[j];

                tly_xz_xzz_0[j] = pa_x[j] * tly_z_xzz_0[j] + 0.5 * fl1_fx * tly_z_zz_0[j] + 0.5 * fl1_fx * tpz_z_xzz_0[j] + fl1_fx * fl1_fgb * tdz_z_xzz_0[j];

                tlz_xz_xzz_0[j] = pa_x[j] * tlz_z_xzz_0[j] + 0.5 * fl1_fx * tlz_z_zz_0[j] - 0.5 * fl1_fx * tpy_z_xzz_0[j] - fl1_fx * fl1_fgb * tdy_z_xzz_0[j];

                tlx_xz_yyy_0[j] = pa_x[j] * tlx_z_yyy_0[j];

                tly_xz_yyy_0[j] = pa_x[j] * tly_z_yyy_0[j] + 0.5 * fl1_fx * tpz_z_yyy_0[j] + fl1_fx * fl1_fgb * tdz_z_yyy_0[j];

                tlz_xz_yyy_0[j] = pa_x[j] * tlz_z_yyy_0[j] - 0.5 * fl1_fx * tpy_z_yyy_0[j] - fl1_fx * fl1_fgb * tdy_z_yyy_0[j];

                tlx_xz_yyz_0[j] = pa_x[j] * tlx_z_yyz_0[j];

                tly_xz_yyz_0[j] = pa_x[j] * tly_z_yyz_0[j] + 0.5 * fl1_fx * tpz_z_yyz_0[j] + fl1_fx * fl1_fgb * tdz_z_yyz_0[j];

                tlz_xz_yyz_0[j] = pa_x[j] * tlz_z_yyz_0[j] - 0.5 * fl1_fx * tpy_z_yyz_0[j] - fl1_fx * fl1_fgb * tdy_z_yyz_0[j];

                tlx_xz_yzz_0[j] = pa_x[j] * tlx_z_yzz_0[j];

                tly_xz_yzz_0[j] = pa_x[j] * tly_z_yzz_0[j] + 0.5 * fl1_fx * tpz_z_yzz_0[j] + fl1_fx * fl1_fgb * tdz_z_yzz_0[j];

                tlz_xz_yzz_0[j] = pa_x[j] * tlz_z_yzz_0[j] - 0.5 * fl1_fx * tpy_z_yzz_0[j] - fl1_fx * fl1_fgb * tdy_z_yzz_0[j];

                tlx_xz_zzz_0[j] = pa_x[j] * tlx_z_zzz_0[j];

                tly_xz_zzz_0[j] = pa_x[j] * tly_z_zzz_0[j] + 0.5 * fl1_fx * tpz_z_zzz_0[j] + fl1_fx * fl1_fgb * tdz_z_zzz_0[j];

                tlz_xz_zzz_0[j] = pa_x[j] * tlz_z_zzz_0[j] - 0.5 * fl1_fx * tpy_z_zzz_0[j] - fl1_fx * fl1_fgb * tdy_z_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDF_90_135(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tlx_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 10); 

            auto tly_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tlz_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tlx_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 11); 

            auto tly_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tlz_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tlx_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 12); 

            auto tly_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tlz_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tlx_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 13); 

            auto tly_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tlz_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tlx_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 14); 

            auto tly_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tlz_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tlx_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 15); 

            auto tly_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tlz_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tlx_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 16); 

            auto tly_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tlz_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tlx_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 17); 

            auto tly_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tlz_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tlx_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 18); 

            auto tly_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tlz_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tlx_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 19); 

            auto tly_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tlz_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tlx_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 20); 

            auto tly_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 21); 

            auto tly_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 22); 

            auto tly_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 23); 

            auto tly_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tlx_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 24); 

            auto tly_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx); 

            auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx); 

            auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx); 

            auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1); 

            auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1); 

            auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1); 

            auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2); 

            auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2); 

            auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2); 

            auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3); 

            auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3); 

            auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3); 

            auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4); 

            auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4); 

            auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4); 

            auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5); 

            auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5); 

            auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5); 

            auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6); 

            auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6); 

            auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6); 

            auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7); 

            auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7); 

            auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7); 

            auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8); 

            auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8); 

            auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8); 

            auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9); 

            auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9); 

            auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9); 

            auto tlx_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 6); 

            auto tly_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tlz_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tlx_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 7); 

            auto tly_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tlz_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tlx_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 8); 

            auto tly_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tlz_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tlx_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 9); 

            auto tly_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 10); 

            auto tly_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tlx_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 11); 

            auto tly_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10); 

            auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11); 

            auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12); 

            auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13); 

            auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14); 

            auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15); 

            auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16); 

            auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17); 

            auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18); 

            auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19); 

            auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20); 

            auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21); 

            auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22); 

            auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23); 

            auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24); 

            auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tdx_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 10); 

            auto tdz_y_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tdx_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 11); 

            auto tdz_y_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tdx_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 12); 

            auto tdz_y_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tdx_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 13); 

            auto tdz_y_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tdx_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 14); 

            auto tdz_y_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tdx_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 15); 

            auto tdz_y_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tdx_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 16); 

            auto tdz_y_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tdx_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 17); 

            auto tdz_y_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tdx_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 18); 

            auto tdz_y_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tdx_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 19); 

            auto tdz_y_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdz_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdz_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdz_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdz_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdz_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            // set up pointers to integrals

            auto tlx_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 30); 

            auto tly_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 30); 

            auto tlz_yy_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 30); 

            auto tlx_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 31); 

            auto tly_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 31); 

            auto tlz_yy_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 31); 

            auto tlx_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 32); 

            auto tly_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 32); 

            auto tlz_yy_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 32); 

            auto tlx_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 33); 

            auto tly_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 33); 

            auto tlz_yy_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 33); 

            auto tlx_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 34); 

            auto tly_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 34); 

            auto tlz_yy_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 34); 

            auto tlx_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 35); 

            auto tly_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 35); 

            auto tlz_yy_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 35); 

            auto tlx_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 36); 

            auto tly_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 36); 

            auto tlz_yy_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 36); 

            auto tlx_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 37); 

            auto tly_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 37); 

            auto tlz_yy_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 37); 

            auto tlx_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 38); 

            auto tly_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 38); 

            auto tlz_yy_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 38); 

            auto tlx_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 39); 

            auto tly_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 39); 

            auto tlz_yy_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 39); 

            auto tlx_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 40); 

            auto tly_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 40); 

            auto tlz_yz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 40); 

            auto tlx_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 41); 

            auto tly_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 41); 

            auto tlz_yz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 41); 

            auto tlx_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 42); 

            auto tly_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 42); 

            auto tlz_yz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 42); 

            auto tlx_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 43); 

            auto tly_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 43); 

            auto tlz_yz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 43); 

            auto tlx_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 44); 

            auto tly_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 44); 

            auto tlz_yz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fgb, fx, pa_y, tdx_y_xxx_0, tdx_y_xxy_0, tdx_y_xxz_0, tdx_y_xyy_0, \
                                     tdx_y_xyz_0, tdx_y_xzz_0, tdx_y_yyy_0, tdx_y_yyz_0, tdx_y_yzz_0, tdx_y_zzz_0, \
                                     tdx_z_xxx_0, tdx_z_xxy_0, tdx_z_xxz_0, tdx_z_xyy_0, tdx_z_xyz_0, tdz_y_xxx_0, \
                                     tdz_y_xxy_0, tdz_y_xxz_0, tdz_y_xyy_0, tdz_y_xyz_0, tdz_y_xzz_0, tdz_y_yyy_0, \
                                     tdz_y_yyz_0, tdz_y_yzz_0, tdz_y_zzz_0, tdz_z_xxx_0, tdz_z_xxy_0, tdz_z_xxz_0, \
                                     tdz_z_xyy_0, tdz_z_xyz_0, tlx_0_xxx_0, tlx_0_xxy_0, tlx_0_xxz_0, tlx_0_xyy_0, \
                                     tlx_0_xyz_0, tlx_0_xzz_0, tlx_0_yyy_0, tlx_0_yyz_0, tlx_0_yzz_0, tlx_0_zzz_0, \
                                     tlx_y_xx_0, tlx_y_xxx_0, tlx_y_xxy_0, tlx_y_xxz_0, tlx_y_xy_0, tlx_y_xyy_0, \
                                     tlx_y_xyz_0, tlx_y_xz_0, tlx_y_xzz_0, tlx_y_yy_0, tlx_y_yyy_0, tlx_y_yyz_0, \
                                     tlx_y_yz_0, tlx_y_yzz_0, tlx_y_zz_0, tlx_y_zzz_0, tlx_yy_xxx_0, tlx_yy_xxy_0, \
                                     tlx_yy_xxz_0, tlx_yy_xyy_0, tlx_yy_xyz_0, tlx_yy_xzz_0, tlx_yy_yyy_0, tlx_yy_yyz_0, \
                                     tlx_yy_yzz_0, tlx_yy_zzz_0, tlx_yz_xxx_0, tlx_yz_xxy_0, tlx_yz_xxz_0, tlx_yz_xyy_0, \
                                     tlx_yz_xyz_0, tlx_z_xx_0, tlx_z_xxx_0, tlx_z_xxy_0, tlx_z_xxz_0, tlx_z_xy_0, \
                                     tlx_z_xyy_0, tlx_z_xyz_0, tlx_z_xz_0, tly_0_xxx_0, tly_0_xxy_0, tly_0_xxz_0, \
                                     tly_0_xyy_0, tly_0_xyz_0, tly_0_xzz_0, tly_0_yyy_0, tly_0_yyz_0, tly_0_yzz_0, \
                                     tly_0_zzz_0, tly_y_xx_0, tly_y_xxx_0, tly_y_xxy_0, tly_y_xxz_0, tly_y_xy_0, \
                                     tly_y_xyy_0, tly_y_xyz_0, tly_y_xz_0, tly_y_xzz_0, tly_y_yy_0, tly_y_yyy_0, \
                                     tly_y_yyz_0, tly_y_yz_0, tly_y_yzz_0, tly_y_zz_0, tly_y_zzz_0, tly_yy_xxx_0, \
                                     tly_yy_xxy_0, tly_yy_xxz_0, tly_yy_xyy_0, tly_yy_xyz_0, tly_yy_xzz_0, tly_yy_yyy_0, \
                                     tly_yy_yyz_0, tly_yy_yzz_0, tly_yy_zzz_0, tly_yz_xxx_0, tly_yz_xxy_0, tly_yz_xxz_0, \
                                     tly_yz_xyy_0, tly_yz_xyz_0, tly_z_xx_0, tly_z_xxx_0, tly_z_xxy_0, tly_z_xxz_0, \
                                     tly_z_xy_0, tly_z_xyy_0, tly_z_xyz_0, tly_z_xz_0, tlz_0_xxx_0, tlz_0_xxy_0, \
                                     tlz_0_xxz_0, tlz_0_xyy_0, tlz_0_xyz_0, tlz_0_xzz_0, tlz_0_yyy_0, tlz_0_yyz_0, \
                                     tlz_0_yzz_0, tlz_0_zzz_0, tlz_y_xx_0, tlz_y_xxx_0, tlz_y_xxy_0, tlz_y_xxz_0, \
                                     tlz_y_xy_0, tlz_y_xyy_0, tlz_y_xyz_0, tlz_y_xz_0, tlz_y_xzz_0, tlz_y_yy_0, \
                                     tlz_y_yyy_0, tlz_y_yyz_0, tlz_y_yz_0, tlz_y_yzz_0, tlz_y_zz_0, tlz_y_zzz_0, \
                                     tlz_yy_xxx_0, tlz_yy_xxy_0, tlz_yy_xxz_0, tlz_yy_xyy_0, tlz_yy_xyz_0, tlz_yy_xzz_0, \
                                     tlz_yy_yyy_0, tlz_yy_yyz_0, tlz_yy_yzz_0, tlz_yy_zzz_0, tlz_yz_xxx_0, tlz_yz_xxy_0, \
                                     tlz_yz_xxz_0, tlz_yz_xyy_0, tlz_yz_xyz_0, tlz_z_xx_0, tlz_z_xxx_0, tlz_z_xxy_0, \
                                     tlz_z_xxz_0, tlz_z_xy_0, tlz_z_xyy_0, tlz_z_xyz_0, tlz_z_xz_0, tpx_y_xxx_0, \
                                     tpx_y_xxy_0, tpx_y_xxz_0, tpx_y_xyy_0, tpx_y_xyz_0, tpx_y_xzz_0, tpx_y_yyy_0, \
                                     tpx_y_yyz_0, tpx_y_yzz_0, tpx_y_zzz_0, tpx_z_xxx_0, tpx_z_xxy_0, tpx_z_xxz_0, \
                                     tpx_z_xyy_0, tpx_z_xyz_0, tpz_y_xxx_0, tpz_y_xxy_0, tpz_y_xxz_0, tpz_y_xyy_0, \
                                     tpz_y_xyz_0, tpz_y_xzz_0, tpz_y_yyy_0, tpz_y_yyz_0, tpz_y_yzz_0, tpz_y_zzz_0, \
                                     tpz_z_xxx_0, tpz_z_xxy_0, tpz_z_xxz_0, tpz_z_xyy_0, tpz_z_xyz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yy_xxx_0[j] = pa_y[j] * tlx_y_xxx_0[j] + 0.5 * fl1_fx * tlx_0_xxx_0[j] - 0.5 * fl1_fx * tpz_y_xxx_0[j] - fl1_fx * fl1_fgb * tdz_y_xxx_0[j];

                tly_yy_xxx_0[j] = pa_y[j] * tly_y_xxx_0[j] + 0.5 * fl1_fx * tly_0_xxx_0[j];

                tlz_yy_xxx_0[j] = pa_y[j] * tlz_y_xxx_0[j] + 0.5 * fl1_fx * tlz_0_xxx_0[j] + 0.5 * fl1_fx * tpx_y_xxx_0[j] + fl1_fx * fl1_fgb * tdx_y_xxx_0[j];

                tlx_yy_xxy_0[j] = pa_y[j] * tlx_y_xxy_0[j] + 0.5 * fl1_fx * tlx_0_xxy_0[j] + 0.5 * fl1_fx * tlx_y_xx_0[j] - 0.5 * fl1_fx * tpz_y_xxy_0[j] - fl1_fx * fl1_fgb * tdz_y_xxy_0[j];

                tly_yy_xxy_0[j] = pa_y[j] * tly_y_xxy_0[j] + 0.5 * fl1_fx * tly_0_xxy_0[j] + 0.5 * fl1_fx * tly_y_xx_0[j];

                tlz_yy_xxy_0[j] = pa_y[j] * tlz_y_xxy_0[j] + 0.5 * fl1_fx * tlz_0_xxy_0[j] + 0.5 * fl1_fx * tlz_y_xx_0[j] + 0.5 * fl1_fx * tpx_y_xxy_0[j] + fl1_fx * fl1_fgb * tdx_y_xxy_0[j];

                tlx_yy_xxz_0[j] = pa_y[j] * tlx_y_xxz_0[j] + 0.5 * fl1_fx * tlx_0_xxz_0[j] - 0.5 * fl1_fx * tpz_y_xxz_0[j] - fl1_fx * fl1_fgb * tdz_y_xxz_0[j];

                tly_yy_xxz_0[j] = pa_y[j] * tly_y_xxz_0[j] + 0.5 * fl1_fx * tly_0_xxz_0[j];

                tlz_yy_xxz_0[j] = pa_y[j] * tlz_y_xxz_0[j] + 0.5 * fl1_fx * tlz_0_xxz_0[j] + 0.5 * fl1_fx * tpx_y_xxz_0[j] + fl1_fx * fl1_fgb * tdx_y_xxz_0[j];

                tlx_yy_xyy_0[j] = pa_y[j] * tlx_y_xyy_0[j] + 0.5 * fl1_fx * tlx_0_xyy_0[j] + fl1_fx * tlx_y_xy_0[j] - 0.5 * fl1_fx * tpz_y_xyy_0[j] - fl1_fx * fl1_fgb * tdz_y_xyy_0[j];

                tly_yy_xyy_0[j] = pa_y[j] * tly_y_xyy_0[j] + 0.5 * fl1_fx * tly_0_xyy_0[j] + fl1_fx * tly_y_xy_0[j];

                tlz_yy_xyy_0[j] = pa_y[j] * tlz_y_xyy_0[j] + 0.5 * fl1_fx * tlz_0_xyy_0[j] + fl1_fx * tlz_y_xy_0[j] + 0.5 * fl1_fx * tpx_y_xyy_0[j] + fl1_fx * fl1_fgb * tdx_y_xyy_0[j];

                tlx_yy_xyz_0[j] = pa_y[j] * tlx_y_xyz_0[j] + 0.5 * fl1_fx * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tlx_y_xz_0[j] - 0.5 * fl1_fx * tpz_y_xyz_0[j] - fl1_fx * fl1_fgb * tdz_y_xyz_0[j];

                tly_yy_xyz_0[j] = pa_y[j] * tly_y_xyz_0[j] + 0.5 * fl1_fx * tly_0_xyz_0[j] + 0.5 * fl1_fx * tly_y_xz_0[j];

                tlz_yy_xyz_0[j] = pa_y[j] * tlz_y_xyz_0[j] + 0.5 * fl1_fx * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tlz_y_xz_0[j] + 0.5 * fl1_fx * tpx_y_xyz_0[j] + fl1_fx * fl1_fgb * tdx_y_xyz_0[j];

                tlx_yy_xzz_0[j] = pa_y[j] * tlx_y_xzz_0[j] + 0.5 * fl1_fx * tlx_0_xzz_0[j] - 0.5 * fl1_fx * tpz_y_xzz_0[j] - fl1_fx * fl1_fgb * tdz_y_xzz_0[j];

                tly_yy_xzz_0[j] = pa_y[j] * tly_y_xzz_0[j] + 0.5 * fl1_fx * tly_0_xzz_0[j];

                tlz_yy_xzz_0[j] = pa_y[j] * tlz_y_xzz_0[j] + 0.5 * fl1_fx * tlz_0_xzz_0[j] + 0.5 * fl1_fx * tpx_y_xzz_0[j] + fl1_fx * fl1_fgb * tdx_y_xzz_0[j];

                tlx_yy_yyy_0[j] = pa_y[j] * tlx_y_yyy_0[j] + 0.5 * fl1_fx * tlx_0_yyy_0[j] + 1.5 * fl1_fx * tlx_y_yy_0[j] - 0.5 * fl1_fx * tpz_y_yyy_0[j] - fl1_fx * fl1_fgb * tdz_y_yyy_0[j];

                tly_yy_yyy_0[j] = pa_y[j] * tly_y_yyy_0[j] + 0.5 * fl1_fx * tly_0_yyy_0[j] + 1.5 * fl1_fx * tly_y_yy_0[j];

                tlz_yy_yyy_0[j] = pa_y[j] * tlz_y_yyy_0[j] + 0.5 * fl1_fx * tlz_0_yyy_0[j] + 1.5 * fl1_fx * tlz_y_yy_0[j] + 0.5 * fl1_fx * tpx_y_yyy_0[j] + fl1_fx * fl1_fgb * tdx_y_yyy_0[j];

                tlx_yy_yyz_0[j] = pa_y[j] * tlx_y_yyz_0[j] + 0.5 * fl1_fx * tlx_0_yyz_0[j] + fl1_fx * tlx_y_yz_0[j] - 0.5 * fl1_fx * tpz_y_yyz_0[j] - fl1_fx * fl1_fgb * tdz_y_yyz_0[j];

                tly_yy_yyz_0[j] = pa_y[j] * tly_y_yyz_0[j] + 0.5 * fl1_fx * tly_0_yyz_0[j] + fl1_fx * tly_y_yz_0[j];

                tlz_yy_yyz_0[j] = pa_y[j] * tlz_y_yyz_0[j] + 0.5 * fl1_fx * tlz_0_yyz_0[j] + fl1_fx * tlz_y_yz_0[j] + 0.5 * fl1_fx * tpx_y_yyz_0[j] + fl1_fx * fl1_fgb * tdx_y_yyz_0[j];

                tlx_yy_yzz_0[j] = pa_y[j] * tlx_y_yzz_0[j] + 0.5 * fl1_fx * tlx_0_yzz_0[j] + 0.5 * fl1_fx * tlx_y_zz_0[j] - 0.5 * fl1_fx * tpz_y_yzz_0[j] - fl1_fx * fl1_fgb * tdz_y_yzz_0[j];

                tly_yy_yzz_0[j] = pa_y[j] * tly_y_yzz_0[j] + 0.5 * fl1_fx * tly_0_yzz_0[j] + 0.5 * fl1_fx * tly_y_zz_0[j];

                tlz_yy_yzz_0[j] = pa_y[j] * tlz_y_yzz_0[j] + 0.5 * fl1_fx * tlz_0_yzz_0[j] + 0.5 * fl1_fx * tlz_y_zz_0[j] + 0.5 * fl1_fx * tpx_y_yzz_0[j] + fl1_fx * fl1_fgb * tdx_y_yzz_0[j];

                tlx_yy_zzz_0[j] = pa_y[j] * tlx_y_zzz_0[j] + 0.5 * fl1_fx * tlx_0_zzz_0[j] - 0.5 * fl1_fx * tpz_y_zzz_0[j] - fl1_fx * fl1_fgb * tdz_y_zzz_0[j];

                tly_yy_zzz_0[j] = pa_y[j] * tly_y_zzz_0[j] + 0.5 * fl1_fx * tly_0_zzz_0[j];

                tlz_yy_zzz_0[j] = pa_y[j] * tlz_y_zzz_0[j] + 0.5 * fl1_fx * tlz_0_zzz_0[j] + 0.5 * fl1_fx * tpx_y_zzz_0[j] + fl1_fx * fl1_fgb * tdx_y_zzz_0[j];

                tlx_yz_xxx_0[j] = pa_y[j] * tlx_z_xxx_0[j] - 0.5 * fl1_fx * tpz_z_xxx_0[j] - fl1_fx * fl1_fgb * tdz_z_xxx_0[j];

                tly_yz_xxx_0[j] = pa_y[j] * tly_z_xxx_0[j];

                tlz_yz_xxx_0[j] = pa_y[j] * tlz_z_xxx_0[j] + 0.5 * fl1_fx * tpx_z_xxx_0[j] + fl1_fx * fl1_fgb * tdx_z_xxx_0[j];

                tlx_yz_xxy_0[j] = pa_y[j] * tlx_z_xxy_0[j] + 0.5 * fl1_fx * tlx_z_xx_0[j] - 0.5 * fl1_fx * tpz_z_xxy_0[j] - fl1_fx * fl1_fgb * tdz_z_xxy_0[j];

                tly_yz_xxy_0[j] = pa_y[j] * tly_z_xxy_0[j] + 0.5 * fl1_fx * tly_z_xx_0[j];

                tlz_yz_xxy_0[j] = pa_y[j] * tlz_z_xxy_0[j] + 0.5 * fl1_fx * tlz_z_xx_0[j] + 0.5 * fl1_fx * tpx_z_xxy_0[j] + fl1_fx * fl1_fgb * tdx_z_xxy_0[j];

                tlx_yz_xxz_0[j] = pa_y[j] * tlx_z_xxz_0[j] - 0.5 * fl1_fx * tpz_z_xxz_0[j] - fl1_fx * fl1_fgb * tdz_z_xxz_0[j];

                tly_yz_xxz_0[j] = pa_y[j] * tly_z_xxz_0[j];

                tlz_yz_xxz_0[j] = pa_y[j] * tlz_z_xxz_0[j] + 0.5 * fl1_fx * tpx_z_xxz_0[j] + fl1_fx * fl1_fgb * tdx_z_xxz_0[j];

                tlx_yz_xyy_0[j] = pa_y[j] * tlx_z_xyy_0[j] + fl1_fx * tlx_z_xy_0[j] - 0.5 * fl1_fx * tpz_z_xyy_0[j] - fl1_fx * fl1_fgb * tdz_z_xyy_0[j];

                tly_yz_xyy_0[j] = pa_y[j] * tly_z_xyy_0[j] + fl1_fx * tly_z_xy_0[j];

                tlz_yz_xyy_0[j] = pa_y[j] * tlz_z_xyy_0[j] + fl1_fx * tlz_z_xy_0[j] + 0.5 * fl1_fx * tpx_z_xyy_0[j] + fl1_fx * fl1_fgb * tdx_z_xyy_0[j];

                tlx_yz_xyz_0[j] = pa_y[j] * tlx_z_xyz_0[j] + 0.5 * fl1_fx * tlx_z_xz_0[j] - 0.5 * fl1_fx * tpz_z_xyz_0[j] - fl1_fx * fl1_fgb * tdz_z_xyz_0[j];

                tly_yz_xyz_0[j] = pa_y[j] * tly_z_xyz_0[j] + 0.5 * fl1_fx * tly_z_xz_0[j];

                tlz_yz_xyz_0[j] = pa_y[j] * tlz_z_xyz_0[j] + 0.5 * fl1_fx * tlz_z_xz_0[j] + 0.5 * fl1_fx * tpx_z_xyz_0[j] + fl1_fx * fl1_fgb * tdx_z_xyz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDF_135_180(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_3_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tlx_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 20); 

            auto tly_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 21); 

            auto tly_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 22); 

            auto tly_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 23); 

            auto tly_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tlx_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 24); 

            auto tly_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tlx_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 25); 

            auto tly_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tlz_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tlx_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 26); 

            auto tly_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tlz_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tlx_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 27); 

            auto tly_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tlz_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tlx_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 28); 

            auto tly_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tlz_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tlx_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 29); 

            auto tly_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tlz_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx); 

            auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx); 

            auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx); 

            auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1); 

            auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1); 

            auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1); 

            auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2); 

            auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2); 

            auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2); 

            auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3); 

            auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3); 

            auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3); 

            auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4); 

            auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4); 

            auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4); 

            auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5); 

            auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5); 

            auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5); 

            auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6); 

            auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6); 

            auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6); 

            auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7); 

            auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7); 

            auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7); 

            auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8); 

            auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8); 

            auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8); 

            auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9); 

            auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9); 

            auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9); 

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 15); 

            auto tly_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 16); 

            auto tly_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 17); 

            auto tly_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20); 

            auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21); 

            auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22); 

            auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23); 

            auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24); 

            auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25); 

            auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26); 

            auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27); 

            auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28); 

            auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29); 

            auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            auto tdx_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 20); 

            auto tdy_z_xxx_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tdx_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 21); 

            auto tdy_z_xxy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tdx_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 22); 

            auto tdy_z_xxz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tdx_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 23); 

            auto tdy_z_xyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tdx_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 24); 

            auto tdy_z_xyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tdx_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 25); 

            auto tdy_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tdz_z_xzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tdx_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 26); 

            auto tdy_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tdz_z_yyy_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tdx_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 27); 

            auto tdy_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tdz_z_yyz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tdx_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 28); 

            auto tdy_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tdz_z_yzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tdx_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * idx + 29); 

            auto tdy_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tdz_z_zzz_0 = primBuffer.data(pidx_d_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            // set up pointers to integrals

            auto tlx_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 45); 

            auto tly_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 45); 

            auto tlz_yz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 45); 

            auto tlx_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 46); 

            auto tly_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 46); 

            auto tlz_yz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 46); 

            auto tlx_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 47); 

            auto tly_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 47); 

            auto tlz_yz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 47); 

            auto tlx_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 48); 

            auto tly_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 48); 

            auto tlz_yz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 48); 

            auto tlx_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 49); 

            auto tly_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 49); 

            auto tlz_yz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 49); 

            auto tlx_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 50); 

            auto tly_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 50); 

            auto tlz_zz_xxx_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 50); 

            auto tlx_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 51); 

            auto tly_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 51); 

            auto tlz_zz_xxy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 51); 

            auto tlx_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 52); 

            auto tly_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 52); 

            auto tlz_zz_xxz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 52); 

            auto tlx_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 53); 

            auto tly_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 53); 

            auto tlz_zz_xyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 53); 

            auto tlx_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 54); 

            auto tly_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 54); 

            auto tlz_zz_xyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 54); 

            auto tlx_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 55); 

            auto tly_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 55); 

            auto tlz_zz_xzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 55); 

            auto tlx_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 56); 

            auto tly_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 56); 

            auto tlz_zz_yyy_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 56); 

            auto tlx_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 57); 

            auto tly_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 57); 

            auto tlz_zz_yyz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 57); 

            auto tlx_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 58); 

            auto tly_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 58); 

            auto tlz_zz_yzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 58); 

            auto tlx_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * idx + 59); 

            auto tly_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 60 * bdim + 60 * idx + 59); 

            auto tlz_zz_zzz_0 = primBuffer.data(pidx_l_2_3_m0 + 120 * bdim + 60 * idx + 59); 

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_z_xxx_0, tdx_z_xxy_0, tdx_z_xxz_0, tdx_z_xyy_0, \
                                     tdx_z_xyz_0, tdx_z_xzz_0, tdx_z_yyy_0, tdx_z_yyz_0, tdx_z_yzz_0, tdx_z_zzz_0, \
                                     tdy_z_xxx_0, tdy_z_xxy_0, tdy_z_xxz_0, tdy_z_xyy_0, tdy_z_xyz_0, tdy_z_xzz_0, \
                                     tdy_z_yyy_0, tdy_z_yyz_0, tdy_z_yzz_0, tdy_z_zzz_0, tdz_z_xzz_0, tdz_z_yyy_0, \
                                     tdz_z_yyz_0, tdz_z_yzz_0, tdz_z_zzz_0, tlx_0_xxx_0, tlx_0_xxy_0, tlx_0_xxz_0, \
                                     tlx_0_xyy_0, tlx_0_xyz_0, tlx_0_xzz_0, tlx_0_yyy_0, tlx_0_yyz_0, tlx_0_yzz_0, \
                                     tlx_0_zzz_0, tlx_yz_xzz_0, tlx_yz_yyy_0, tlx_yz_yyz_0, tlx_yz_yzz_0, tlx_yz_zzz_0, \
                                     tlx_z_xx_0, tlx_z_xxx_0, tlx_z_xxy_0, tlx_z_xxz_0, tlx_z_xy_0, tlx_z_xyy_0, \
                                     tlx_z_xyz_0, tlx_z_xz_0, tlx_z_xzz_0, tlx_z_yy_0, tlx_z_yyy_0, tlx_z_yyz_0, \
                                     tlx_z_yz_0, tlx_z_yzz_0, tlx_z_zz_0, tlx_z_zzz_0, tlx_zz_xxx_0, tlx_zz_xxy_0, \
                                     tlx_zz_xxz_0, tlx_zz_xyy_0, tlx_zz_xyz_0, tlx_zz_xzz_0, tlx_zz_yyy_0, tlx_zz_yyz_0, \
                                     tlx_zz_yzz_0, tlx_zz_zzz_0, tly_0_xxx_0, tly_0_xxy_0, tly_0_xxz_0, tly_0_xyy_0, \
                                     tly_0_xyz_0, tly_0_xzz_0, tly_0_yyy_0, tly_0_yyz_0, tly_0_yzz_0, tly_0_zzz_0, \
                                     tly_yz_xzz_0, tly_yz_yyy_0, tly_yz_yyz_0, tly_yz_yzz_0, tly_yz_zzz_0, tly_z_xx_0, \
                                     tly_z_xxx_0, tly_z_xxy_0, tly_z_xxz_0, tly_z_xy_0, tly_z_xyy_0, tly_z_xyz_0, \
                                     tly_z_xz_0, tly_z_xzz_0, tly_z_yy_0, tly_z_yyy_0, tly_z_yyz_0, tly_z_yz_0, \
                                     tly_z_yzz_0, tly_z_zz_0, tly_z_zzz_0, tly_zz_xxx_0, tly_zz_xxy_0, tly_zz_xxz_0, \
                                     tly_zz_xyy_0, tly_zz_xyz_0, tly_zz_xzz_0, tly_zz_yyy_0, tly_zz_yyz_0, tly_zz_yzz_0, \
                                     tly_zz_zzz_0, tlz_0_xxx_0, tlz_0_xxy_0, tlz_0_xxz_0, tlz_0_xyy_0, tlz_0_xyz_0, \
                                     tlz_0_xzz_0, tlz_0_yyy_0, tlz_0_yyz_0, tlz_0_yzz_0, tlz_0_zzz_0, tlz_yz_xzz_0, \
                                     tlz_yz_yyy_0, tlz_yz_yyz_0, tlz_yz_yzz_0, tlz_yz_zzz_0, tlz_z_xx_0, tlz_z_xxx_0, \
                                     tlz_z_xxy_0, tlz_z_xxz_0, tlz_z_xy_0, tlz_z_xyy_0, tlz_z_xyz_0, tlz_z_xz_0, \
                                     tlz_z_xzz_0, tlz_z_yy_0, tlz_z_yyy_0, tlz_z_yyz_0, tlz_z_yz_0, tlz_z_yzz_0, \
                                     tlz_z_zz_0, tlz_z_zzz_0, tlz_zz_xxx_0, tlz_zz_xxy_0, tlz_zz_xxz_0, tlz_zz_xyy_0, \
                                     tlz_zz_xyz_0, tlz_zz_xzz_0, tlz_zz_yyy_0, tlz_zz_yyz_0, tlz_zz_yzz_0, tlz_zz_zzz_0, \
                                     tpx_z_xxx_0, tpx_z_xxy_0, tpx_z_xxz_0, tpx_z_xyy_0, tpx_z_xyz_0, tpx_z_xzz_0, \
                                     tpx_z_yyy_0, tpx_z_yyz_0, tpx_z_yzz_0, tpx_z_zzz_0, tpy_z_xxx_0, tpy_z_xxy_0, \
                                     tpy_z_xxz_0, tpy_z_xyy_0, tpy_z_xyz_0, tpy_z_xzz_0, tpy_z_yyy_0, tpy_z_yyz_0, \
                                     tpy_z_yzz_0, tpy_z_zzz_0, tpz_z_xzz_0, tpz_z_yyy_0, tpz_z_yyz_0, tpz_z_yzz_0, \
                                     tpz_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yz_xzz_0[j] = pa_y[j] * tlx_z_xzz_0[j] - 0.5 * fl1_fx * tpz_z_xzz_0[j] - fl1_fx * fl1_fgb * tdz_z_xzz_0[j];

                tly_yz_xzz_0[j] = pa_y[j] * tly_z_xzz_0[j];

                tlz_yz_xzz_0[j] = pa_y[j] * tlz_z_xzz_0[j] + 0.5 * fl1_fx * tpx_z_xzz_0[j] + fl1_fx * fl1_fgb * tdx_z_xzz_0[j];

                tlx_yz_yyy_0[j] = pa_y[j] * tlx_z_yyy_0[j] + 1.5 * fl1_fx * tlx_z_yy_0[j] - 0.5 * fl1_fx * tpz_z_yyy_0[j] - fl1_fx * fl1_fgb * tdz_z_yyy_0[j];

                tly_yz_yyy_0[j] = pa_y[j] * tly_z_yyy_0[j] + 1.5 * fl1_fx * tly_z_yy_0[j];

                tlz_yz_yyy_0[j] = pa_y[j] * tlz_z_yyy_0[j] + 1.5 * fl1_fx * tlz_z_yy_0[j] + 0.5 * fl1_fx * tpx_z_yyy_0[j] + fl1_fx * fl1_fgb * tdx_z_yyy_0[j];

                tlx_yz_yyz_0[j] = pa_y[j] * tlx_z_yyz_0[j] + fl1_fx * tlx_z_yz_0[j] - 0.5 * fl1_fx * tpz_z_yyz_0[j] - fl1_fx * fl1_fgb * tdz_z_yyz_0[j];

                tly_yz_yyz_0[j] = pa_y[j] * tly_z_yyz_0[j] + fl1_fx * tly_z_yz_0[j];

                tlz_yz_yyz_0[j] = pa_y[j] * tlz_z_yyz_0[j] + fl1_fx * tlz_z_yz_0[j] + 0.5 * fl1_fx * tpx_z_yyz_0[j] + fl1_fx * fl1_fgb * tdx_z_yyz_0[j];

                tlx_yz_yzz_0[j] = pa_y[j] * tlx_z_yzz_0[j] + 0.5 * fl1_fx * tlx_z_zz_0[j] - 0.5 * fl1_fx * tpz_z_yzz_0[j] - fl1_fx * fl1_fgb * tdz_z_yzz_0[j];

                tly_yz_yzz_0[j] = pa_y[j] * tly_z_yzz_0[j] + 0.5 * fl1_fx * tly_z_zz_0[j];

                tlz_yz_yzz_0[j] = pa_y[j] * tlz_z_yzz_0[j] + 0.5 * fl1_fx * tlz_z_zz_0[j] + 0.5 * fl1_fx * tpx_z_yzz_0[j] + fl1_fx * fl1_fgb * tdx_z_yzz_0[j];

                tlx_yz_zzz_0[j] = pa_y[j] * tlx_z_zzz_0[j] - 0.5 * fl1_fx * tpz_z_zzz_0[j] - fl1_fx * fl1_fgb * tdz_z_zzz_0[j];

                tly_yz_zzz_0[j] = pa_y[j] * tly_z_zzz_0[j];

                tlz_yz_zzz_0[j] = pa_y[j] * tlz_z_zzz_0[j] + 0.5 * fl1_fx * tpx_z_zzz_0[j] + fl1_fx * fl1_fgb * tdx_z_zzz_0[j];

                tlx_zz_xxx_0[j] = pa_z[j] * tlx_z_xxx_0[j] + 0.5 * fl1_fx * tlx_0_xxx_0[j] + 0.5 * fl1_fx * tpy_z_xxx_0[j] + fl1_fx * fl1_fgb * tdy_z_xxx_0[j];

                tly_zz_xxx_0[j] = pa_z[j] * tly_z_xxx_0[j] + 0.5 * fl1_fx * tly_0_xxx_0[j] - 0.5 * fl1_fx * tpx_z_xxx_0[j] - fl1_fx * fl1_fgb * tdx_z_xxx_0[j];

                tlz_zz_xxx_0[j] = pa_z[j] * tlz_z_xxx_0[j] + 0.5 * fl1_fx * tlz_0_xxx_0[j];

                tlx_zz_xxy_0[j] = pa_z[j] * tlx_z_xxy_0[j] + 0.5 * fl1_fx * tlx_0_xxy_0[j] + 0.5 * fl1_fx * tpy_z_xxy_0[j] + fl1_fx * fl1_fgb * tdy_z_xxy_0[j];

                tly_zz_xxy_0[j] = pa_z[j] * tly_z_xxy_0[j] + 0.5 * fl1_fx * tly_0_xxy_0[j] - 0.5 * fl1_fx * tpx_z_xxy_0[j] - fl1_fx * fl1_fgb * tdx_z_xxy_0[j];

                tlz_zz_xxy_0[j] = pa_z[j] * tlz_z_xxy_0[j] + 0.5 * fl1_fx * tlz_0_xxy_0[j];

                tlx_zz_xxz_0[j] = pa_z[j] * tlx_z_xxz_0[j] + 0.5 * fl1_fx * tlx_0_xxz_0[j] + 0.5 * fl1_fx * tlx_z_xx_0[j] + 0.5 * fl1_fx * tpy_z_xxz_0[j] + fl1_fx * fl1_fgb * tdy_z_xxz_0[j];

                tly_zz_xxz_0[j] = pa_z[j] * tly_z_xxz_0[j] + 0.5 * fl1_fx * tly_0_xxz_0[j] + 0.5 * fl1_fx * tly_z_xx_0[j] - 0.5 * fl1_fx * tpx_z_xxz_0[j] - fl1_fx * fl1_fgb * tdx_z_xxz_0[j];

                tlz_zz_xxz_0[j] = pa_z[j] * tlz_z_xxz_0[j] + 0.5 * fl1_fx * tlz_0_xxz_0[j] + 0.5 * fl1_fx * tlz_z_xx_0[j];

                tlx_zz_xyy_0[j] = pa_z[j] * tlx_z_xyy_0[j] + 0.5 * fl1_fx * tlx_0_xyy_0[j] + 0.5 * fl1_fx * tpy_z_xyy_0[j] + fl1_fx * fl1_fgb * tdy_z_xyy_0[j];

                tly_zz_xyy_0[j] = pa_z[j] * tly_z_xyy_0[j] + 0.5 * fl1_fx * tly_0_xyy_0[j] - 0.5 * fl1_fx * tpx_z_xyy_0[j] - fl1_fx * fl1_fgb * tdx_z_xyy_0[j];

                tlz_zz_xyy_0[j] = pa_z[j] * tlz_z_xyy_0[j] + 0.5 * fl1_fx * tlz_0_xyy_0[j];

                tlx_zz_xyz_0[j] = pa_z[j] * tlx_z_xyz_0[j] + 0.5 * fl1_fx * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tlx_z_xy_0[j] + 0.5 * fl1_fx * tpy_z_xyz_0[j] + fl1_fx * fl1_fgb * tdy_z_xyz_0[j];

                tly_zz_xyz_0[j] = pa_z[j] * tly_z_xyz_0[j] + 0.5 * fl1_fx * tly_0_xyz_0[j] + 0.5 * fl1_fx * tly_z_xy_0[j] - 0.5 * fl1_fx * tpx_z_xyz_0[j] - fl1_fx * fl1_fgb * tdx_z_xyz_0[j];

                tlz_zz_xyz_0[j] = pa_z[j] * tlz_z_xyz_0[j] + 0.5 * fl1_fx * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tlz_z_xy_0[j];

                tlx_zz_xzz_0[j] = pa_z[j] * tlx_z_xzz_0[j] + 0.5 * fl1_fx * tlx_0_xzz_0[j] + fl1_fx * tlx_z_xz_0[j] + 0.5 * fl1_fx * tpy_z_xzz_0[j] + fl1_fx * fl1_fgb * tdy_z_xzz_0[j];

                tly_zz_xzz_0[j] = pa_z[j] * tly_z_xzz_0[j] + 0.5 * fl1_fx * tly_0_xzz_0[j] + fl1_fx * tly_z_xz_0[j] - 0.5 * fl1_fx * tpx_z_xzz_0[j] - fl1_fx * fl1_fgb * tdx_z_xzz_0[j];

                tlz_zz_xzz_0[j] = pa_z[j] * tlz_z_xzz_0[j] + 0.5 * fl1_fx * tlz_0_xzz_0[j] + fl1_fx * tlz_z_xz_0[j];

                tlx_zz_yyy_0[j] = pa_z[j] * tlx_z_yyy_0[j] + 0.5 * fl1_fx * tlx_0_yyy_0[j] + 0.5 * fl1_fx * tpy_z_yyy_0[j] + fl1_fx * fl1_fgb * tdy_z_yyy_0[j];

                tly_zz_yyy_0[j] = pa_z[j] * tly_z_yyy_0[j] + 0.5 * fl1_fx * tly_0_yyy_0[j] - 0.5 * fl1_fx * tpx_z_yyy_0[j] - fl1_fx * fl1_fgb * tdx_z_yyy_0[j];

                tlz_zz_yyy_0[j] = pa_z[j] * tlz_z_yyy_0[j] + 0.5 * fl1_fx * tlz_0_yyy_0[j];

                tlx_zz_yyz_0[j] = pa_z[j] * tlx_z_yyz_0[j] + 0.5 * fl1_fx * tlx_0_yyz_0[j] + 0.5 * fl1_fx * tlx_z_yy_0[j] + 0.5 * fl1_fx * tpy_z_yyz_0[j] + fl1_fx * fl1_fgb * tdy_z_yyz_0[j];

                tly_zz_yyz_0[j] = pa_z[j] * tly_z_yyz_0[j] + 0.5 * fl1_fx * tly_0_yyz_0[j] + 0.5 * fl1_fx * tly_z_yy_0[j] - 0.5 * fl1_fx * tpx_z_yyz_0[j] - fl1_fx * fl1_fgb * tdx_z_yyz_0[j];

                tlz_zz_yyz_0[j] = pa_z[j] * tlz_z_yyz_0[j] + 0.5 * fl1_fx * tlz_0_yyz_0[j] + 0.5 * fl1_fx * tlz_z_yy_0[j];

                tlx_zz_yzz_0[j] = pa_z[j] * tlx_z_yzz_0[j] + 0.5 * fl1_fx * tlx_0_yzz_0[j] + fl1_fx * tlx_z_yz_0[j] + 0.5 * fl1_fx * tpy_z_yzz_0[j] + fl1_fx * fl1_fgb * tdy_z_yzz_0[j];

                tly_zz_yzz_0[j] = pa_z[j] * tly_z_yzz_0[j] + 0.5 * fl1_fx * tly_0_yzz_0[j] + fl1_fx * tly_z_yz_0[j] - 0.5 * fl1_fx * tpx_z_yzz_0[j] - fl1_fx * fl1_fgb * tdx_z_yzz_0[j];

                tlz_zz_yzz_0[j] = pa_z[j] * tlz_z_yzz_0[j] + 0.5 * fl1_fx * tlz_0_yzz_0[j] + fl1_fx * tlz_z_yz_0[j];

                tlx_zz_zzz_0[j] = pa_z[j] * tlx_z_zzz_0[j] + 0.5 * fl1_fx * tlx_0_zzz_0[j] + 1.5 * fl1_fx * tlx_z_zz_0[j] + 0.5 * fl1_fx * tpy_z_zzz_0[j] + fl1_fx * fl1_fgb * tdy_z_zzz_0[j];

                tly_zz_zzz_0[j] = pa_z[j] * tly_z_zzz_0[j] + 0.5 * fl1_fx * tly_0_zzz_0[j] + 1.5 * fl1_fx * tly_z_zz_0[j] - 0.5 * fl1_fx * tpx_z_zzz_0[j] - fl1_fx * fl1_fgb * tdx_z_zzz_0[j];

                tlz_zz_zzz_0[j] = pa_z[j] * tlz_z_zzz_0[j] + 0.5 * fl1_fx * tlz_0_zzz_0[j] + 1.5 * fl1_fx * tlz_z_zz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForFD(      CMemBlock2D<double>& primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        amomrecfunc::compAngularMomentumForFD_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        amomrecfunc::compAngularMomentumForFD_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        amomrecfunc::compAngularMomentumForFD_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        amomrecfunc::compAngularMomentumForFD_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compAngularMomentumForFD_0_45(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx); 

            auto tly_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tlz_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tlx_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 1); 

            auto tly_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tlz_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tlx_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 2); 

            auto tly_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tlz_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tlx_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 3); 

            auto tly_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tlz_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tlx_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 4); 

            auto tly_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tlz_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tlx_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 5); 

            auto tly_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tlz_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tlx_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 6); 

            auto tly_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tlz_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tlx_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 7); 

            auto tly_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tlz_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tlx_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 8); 

            auto tly_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tlz_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tlx_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 9); 

            auto tly_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tlz_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tlx_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 10); 

            auto tly_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tlz_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tlx_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 11); 

            auto tly_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tlz_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            auto tlx_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 12); 

            auto tly_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tlz_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tlx_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 13); 

            auto tly_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tlz_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tlx_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 14); 

            auto tly_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tlz_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            auto tlx_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx); 

            auto tly_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx); 

            auto tlz_x_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx); 

            auto tlx_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 1); 

            auto tly_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 1); 

            auto tlz_x_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 1); 

            auto tlx_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 2); 

            auto tly_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 2); 

            auto tlz_x_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 2); 

            auto tlx_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 3); 

            auto tly_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 3); 

            auto tlz_x_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 3); 

            auto tlx_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 4); 

            auto tly_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 4); 

            auto tlz_x_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 4); 

            auto tlx_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 5); 

            auto tly_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 5); 

            auto tlz_x_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 5); 

            auto tlx_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 6); 

            auto tly_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tlz_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tlx_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 7); 

            auto tly_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tlz_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tlx_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 8); 

            auto tly_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tlz_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tlx_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 9); 

            auto tly_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 10); 

            auto tly_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tlx_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 11); 

            auto tly_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx); 

            auto tly_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx); 

            auto tlz_xx_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx); 

            auto tlx_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 1); 

            auto tly_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 1); 

            auto tlz_xx_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 1); 

            auto tlx_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 2); 

            auto tly_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 2); 

            auto tlz_xx_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 2); 

            auto tlx_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 3); 

            auto tly_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 3); 

            auto tlz_xy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 3); 

            auto tlx_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 4); 

            auto tly_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 4); 

            auto tlz_xy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 4); 

            auto tlx_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 5); 

            auto tly_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 5); 

            auto tlz_xy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 5); 

            auto tlx_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 6); 

            auto tly_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 6); 

            auto tlz_xz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 6); 

            auto tlx_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 7); 

            auto tly_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 7); 

            auto tlz_xz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 7); 

            auto tlx_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 8); 

            auto tly_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 8); 

            auto tlz_xz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 8); 

            auto tpy_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tpz_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tpy_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tpz_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tpy_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tpz_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tpy_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tpz_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tpy_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tpz_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tpy_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tpz_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tpy_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tpz_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tpy_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tpz_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tpy_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tpz_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tpy_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tpz_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tpy_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tpz_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tpy_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tpz_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            auto tpy_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tpz_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tpy_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tpz_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tpy_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tpz_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            auto tdy_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tdz_xx_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tdy_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tdz_xx_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tdy_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tdz_xx_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tdy_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tdz_xx_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tdy_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tdz_xx_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tdy_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tdz_xx_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tdy_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tdz_xy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tdy_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tdz_xy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tdy_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tdz_xy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tdy_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tdz_xy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tdy_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tdz_xy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tdy_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tdz_xy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            auto tdy_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tdz_xz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tdy_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tdz_xz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tdy_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tdz_xz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            // set up pointers to integrals

            auto tlx_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx); 

            auto tly_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx); 

            auto tlz_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx); 

            auto tlx_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 1); 

            auto tly_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 1); 

            auto tlz_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 1); 

            auto tlx_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 2); 

            auto tly_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 2); 

            auto tlz_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 2); 

            auto tlx_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 3); 

            auto tly_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 3); 

            auto tlz_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 3); 

            auto tlx_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 4); 

            auto tly_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 4); 

            auto tlz_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 4); 

            auto tlx_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 5); 

            auto tly_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 5); 

            auto tlz_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 5); 

            auto tlx_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 6); 

            auto tly_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 6); 

            auto tlz_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 6); 

            auto tlx_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 7); 

            auto tly_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 7); 

            auto tlz_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 7); 

            auto tlx_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 8); 

            auto tly_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 8); 

            auto tlz_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 8); 

            auto tlx_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 9); 

            auto tly_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 9); 

            auto tlz_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 9); 

            auto tlx_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 10); 

            auto tly_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 10); 

            auto tlz_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 10); 

            auto tlx_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 11); 

            auto tly_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 11); 

            auto tlz_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 11); 

            auto tlx_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 12); 

            auto tly_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 12); 

            auto tlz_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 12); 

            auto tlx_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 13); 

            auto tly_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 13); 

            auto tlz_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 13); 

            auto tlx_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 14); 

            auto tly_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 14); 

            auto tlz_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_xx_xx_0, tdy_xx_xy_0, tdy_xx_xz_0, tdy_xx_yy_0, \
                                     tdy_xx_yz_0, tdy_xx_zz_0, tdy_xy_xx_0, tdy_xy_xy_0, tdy_xy_xz_0, tdy_xy_yy_0, \
                                     tdy_xy_yz_0, tdy_xy_zz_0, tdy_xz_xx_0, tdy_xz_xy_0, tdy_xz_xz_0, tdz_xx_xx_0, \
                                     tdz_xx_xy_0, tdz_xx_xz_0, tdz_xx_yy_0, tdz_xx_yz_0, tdz_xx_zz_0, tdz_xy_xx_0, \
                                     tdz_xy_xy_0, tdz_xy_xz_0, tdz_xy_yy_0, tdz_xy_yz_0, tdz_xy_zz_0, tdz_xz_xx_0, \
                                     tdz_xz_xy_0, tdz_xz_xz_0, tlx_x_xx_0, tlx_x_xy_0, tlx_x_xz_0, tlx_x_yy_0, \
                                     tlx_x_yz_0, tlx_x_zz_0, tlx_xx_x_0, tlx_xx_xx_0, tlx_xx_xy_0, tlx_xx_xz_0, \
                                     tlx_xx_y_0, tlx_xx_yy_0, tlx_xx_yz_0, tlx_xx_z_0, tlx_xx_zz_0, tlx_xxx_xx_0, \
                                     tlx_xxx_xy_0, tlx_xxx_xz_0, tlx_xxx_yy_0, tlx_xxx_yz_0, tlx_xxx_zz_0, tlx_xxy_xx_0, \
                                     tlx_xxy_xy_0, tlx_xxy_xz_0, tlx_xxy_yy_0, tlx_xxy_yz_0, tlx_xxy_zz_0, tlx_xxz_xx_0, \
                                     tlx_xxz_xy_0, tlx_xxz_xz_0, tlx_xy_x_0, tlx_xy_xx_0, tlx_xy_xy_0, tlx_xy_xz_0, \
                                     tlx_xy_y_0, tlx_xy_yy_0, tlx_xy_yz_0, tlx_xy_z_0, tlx_xy_zz_0, tlx_xz_x_0, \
                                     tlx_xz_xx_0, tlx_xz_xy_0, tlx_xz_xz_0, tlx_xz_y_0, tlx_xz_z_0, tlx_y_xx_0, \
                                     tlx_y_xy_0, tlx_y_xz_0, tlx_y_yy_0, tlx_y_yz_0, tlx_y_zz_0, tlx_z_xx_0, tlx_z_xy_0, \
                                     tlx_z_xz_0, tly_x_xx_0, tly_x_xy_0, tly_x_xz_0, tly_x_yy_0, tly_x_yz_0, tly_x_zz_0, \
                                     tly_xx_x_0, tly_xx_xx_0, tly_xx_xy_0, tly_xx_xz_0, tly_xx_y_0, tly_xx_yy_0, \
                                     tly_xx_yz_0, tly_xx_z_0, tly_xx_zz_0, tly_xxx_xx_0, tly_xxx_xy_0, tly_xxx_xz_0, \
                                     tly_xxx_yy_0, tly_xxx_yz_0, tly_xxx_zz_0, tly_xxy_xx_0, tly_xxy_xy_0, tly_xxy_xz_0, \
                                     tly_xxy_yy_0, tly_xxy_yz_0, tly_xxy_zz_0, tly_xxz_xx_0, tly_xxz_xy_0, tly_xxz_xz_0, \
                                     tly_xy_x_0, tly_xy_xx_0, tly_xy_xy_0, tly_xy_xz_0, tly_xy_y_0, tly_xy_yy_0, \
                                     tly_xy_yz_0, tly_xy_z_0, tly_xy_zz_0, tly_xz_x_0, tly_xz_xx_0, tly_xz_xy_0, \
                                     tly_xz_xz_0, tly_xz_y_0, tly_xz_z_0, tly_y_xx_0, tly_y_xy_0, tly_y_xz_0, tly_y_yy_0, \
                                     tly_y_yz_0, tly_y_zz_0, tly_z_xx_0, tly_z_xy_0, tly_z_xz_0, tlz_x_xx_0, tlz_x_xy_0, \
                                     tlz_x_xz_0, tlz_x_yy_0, tlz_x_yz_0, tlz_x_zz_0, tlz_xx_x_0, tlz_xx_xx_0, \
                                     tlz_xx_xy_0, tlz_xx_xz_0, tlz_xx_y_0, tlz_xx_yy_0, tlz_xx_yz_0, tlz_xx_z_0, \
                                     tlz_xx_zz_0, tlz_xxx_xx_0, tlz_xxx_xy_0, tlz_xxx_xz_0, tlz_xxx_yy_0, tlz_xxx_yz_0, \
                                     tlz_xxx_zz_0, tlz_xxy_xx_0, tlz_xxy_xy_0, tlz_xxy_xz_0, tlz_xxy_yy_0, tlz_xxy_yz_0, \
                                     tlz_xxy_zz_0, tlz_xxz_xx_0, tlz_xxz_xy_0, tlz_xxz_xz_0, tlz_xy_x_0, tlz_xy_xx_0, \
                                     tlz_xy_xy_0, tlz_xy_xz_0, tlz_xy_y_0, tlz_xy_yy_0, tlz_xy_yz_0, tlz_xy_z_0, \
                                     tlz_xy_zz_0, tlz_xz_x_0, tlz_xz_xx_0, tlz_xz_xy_0, tlz_xz_xz_0, tlz_xz_y_0, \
                                     tlz_xz_z_0, tlz_y_xx_0, tlz_y_xy_0, tlz_y_xz_0, tlz_y_yy_0, tlz_y_yz_0, tlz_y_zz_0, \
                                     tlz_z_xx_0, tlz_z_xy_0, tlz_z_xz_0, tpy_xx_xx_0, tpy_xx_xy_0, tpy_xx_xz_0, \
                                     tpy_xx_yy_0, tpy_xx_yz_0, tpy_xx_zz_0, tpy_xy_xx_0, tpy_xy_xy_0, tpy_xy_xz_0, \
                                     tpy_xy_yy_0, tpy_xy_yz_0, tpy_xy_zz_0, tpy_xz_xx_0, tpy_xz_xy_0, tpy_xz_xz_0, \
                                     tpz_xx_xx_0, tpz_xx_xy_0, tpz_xx_xz_0, tpz_xx_yy_0, tpz_xx_yz_0, tpz_xx_zz_0, \
                                     tpz_xy_xx_0, tpz_xy_xy_0, tpz_xy_xz_0, tpz_xy_yy_0, tpz_xy_yz_0, tpz_xy_zz_0, \
                                     tpz_xz_xx_0, tpz_xz_xy_0, tpz_xz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xxx_xx_0[j] = pa_x[j] * tlx_xx_xx_0[j] + fl1_fx * tlx_x_xx_0[j] + fl1_fx * tlx_xx_x_0[j];

                tly_xxx_xx_0[j] = pa_x[j] * tly_xx_xx_0[j] + fl1_fx * tly_x_xx_0[j] + fl1_fx * tly_xx_x_0[j] + 0.5 * fl1_fx * tpz_xx_xx_0[j] + fl1_fx * fl1_fgb * tdz_xx_xx_0[j];

                tlz_xxx_xx_0[j] = pa_x[j] * tlz_xx_xx_0[j] + fl1_fx * tlz_x_xx_0[j] + fl1_fx * tlz_xx_x_0[j] - 0.5 * fl1_fx * tpy_xx_xx_0[j] - fl1_fx * fl1_fgb * tdy_xx_xx_0[j];

                tlx_xxx_xy_0[j] = pa_x[j] * tlx_xx_xy_0[j] + fl1_fx * tlx_x_xy_0[j] + 0.5 * fl1_fx * tlx_xx_y_0[j];

                tly_xxx_xy_0[j] = pa_x[j] * tly_xx_xy_0[j] + fl1_fx * tly_x_xy_0[j] + 0.5 * fl1_fx * tly_xx_y_0[j] + 0.5 * fl1_fx * tpz_xx_xy_0[j] + fl1_fx * fl1_fgb * tdz_xx_xy_0[j];

                tlz_xxx_xy_0[j] = pa_x[j] * tlz_xx_xy_0[j] + fl1_fx * tlz_x_xy_0[j] + 0.5 * fl1_fx * tlz_xx_y_0[j] - 0.5 * fl1_fx * tpy_xx_xy_0[j] - fl1_fx * fl1_fgb * tdy_xx_xy_0[j];

                tlx_xxx_xz_0[j] = pa_x[j] * tlx_xx_xz_0[j] + fl1_fx * tlx_x_xz_0[j] + 0.5 * fl1_fx * tlx_xx_z_0[j];

                tly_xxx_xz_0[j] = pa_x[j] * tly_xx_xz_0[j] + fl1_fx * tly_x_xz_0[j] + 0.5 * fl1_fx * tly_xx_z_0[j] + 0.5 * fl1_fx * tpz_xx_xz_0[j] + fl1_fx * fl1_fgb * tdz_xx_xz_0[j];

                tlz_xxx_xz_0[j] = pa_x[j] * tlz_xx_xz_0[j] + fl1_fx * tlz_x_xz_0[j] + 0.5 * fl1_fx * tlz_xx_z_0[j] - 0.5 * fl1_fx * tpy_xx_xz_0[j] - fl1_fx * fl1_fgb * tdy_xx_xz_0[j];

                tlx_xxx_yy_0[j] = pa_x[j] * tlx_xx_yy_0[j] + fl1_fx * tlx_x_yy_0[j];

                tly_xxx_yy_0[j] = pa_x[j] * tly_xx_yy_0[j] + fl1_fx * tly_x_yy_0[j] + 0.5 * fl1_fx * tpz_xx_yy_0[j] + fl1_fx * fl1_fgb * tdz_xx_yy_0[j];

                tlz_xxx_yy_0[j] = pa_x[j] * tlz_xx_yy_0[j] + fl1_fx * tlz_x_yy_0[j] - 0.5 * fl1_fx * tpy_xx_yy_0[j] - fl1_fx * fl1_fgb * tdy_xx_yy_0[j];

                tlx_xxx_yz_0[j] = pa_x[j] * tlx_xx_yz_0[j] + fl1_fx * tlx_x_yz_0[j];

                tly_xxx_yz_0[j] = pa_x[j] * tly_xx_yz_0[j] + fl1_fx * tly_x_yz_0[j] + 0.5 * fl1_fx * tpz_xx_yz_0[j] + fl1_fx * fl1_fgb * tdz_xx_yz_0[j];

                tlz_xxx_yz_0[j] = pa_x[j] * tlz_xx_yz_0[j] + fl1_fx * tlz_x_yz_0[j] - 0.5 * fl1_fx * tpy_xx_yz_0[j] - fl1_fx * fl1_fgb * tdy_xx_yz_0[j];

                tlx_xxx_zz_0[j] = pa_x[j] * tlx_xx_zz_0[j] + fl1_fx * tlx_x_zz_0[j];

                tly_xxx_zz_0[j] = pa_x[j] * tly_xx_zz_0[j] + fl1_fx * tly_x_zz_0[j] + 0.5 * fl1_fx * tpz_xx_zz_0[j] + fl1_fx * fl1_fgb * tdz_xx_zz_0[j];

                tlz_xxx_zz_0[j] = pa_x[j] * tlz_xx_zz_0[j] + fl1_fx * tlz_x_zz_0[j] - 0.5 * fl1_fx * tpy_xx_zz_0[j] - fl1_fx * fl1_fgb * tdy_xx_zz_0[j];

                tlx_xxy_xx_0[j] = pa_x[j] * tlx_xy_xx_0[j] + 0.5 * fl1_fx * tlx_y_xx_0[j] + fl1_fx * tlx_xy_x_0[j];

                tly_xxy_xx_0[j] = pa_x[j] * tly_xy_xx_0[j] + 0.5 * fl1_fx * tly_y_xx_0[j] + fl1_fx * tly_xy_x_0[j] + 0.5 * fl1_fx * tpz_xy_xx_0[j] + fl1_fx * fl1_fgb * tdz_xy_xx_0[j];

                tlz_xxy_xx_0[j] = pa_x[j] * tlz_xy_xx_0[j] + 0.5 * fl1_fx * tlz_y_xx_0[j] + fl1_fx * tlz_xy_x_0[j] - 0.5 * fl1_fx * tpy_xy_xx_0[j] - fl1_fx * fl1_fgb * tdy_xy_xx_0[j];

                tlx_xxy_xy_0[j] = pa_x[j] * tlx_xy_xy_0[j] + 0.5 * fl1_fx * tlx_y_xy_0[j] + 0.5 * fl1_fx * tlx_xy_y_0[j];

                tly_xxy_xy_0[j] = pa_x[j] * tly_xy_xy_0[j] + 0.5 * fl1_fx * tly_y_xy_0[j] + 0.5 * fl1_fx * tly_xy_y_0[j] + 0.5 * fl1_fx * tpz_xy_xy_0[j] + fl1_fx * fl1_fgb * tdz_xy_xy_0[j];

                tlz_xxy_xy_0[j] = pa_x[j] * tlz_xy_xy_0[j] + 0.5 * fl1_fx * tlz_y_xy_0[j] + 0.5 * fl1_fx * tlz_xy_y_0[j] - 0.5 * fl1_fx * tpy_xy_xy_0[j] - fl1_fx * fl1_fgb * tdy_xy_xy_0[j];

                tlx_xxy_xz_0[j] = pa_x[j] * tlx_xy_xz_0[j] + 0.5 * fl1_fx * tlx_y_xz_0[j] + 0.5 * fl1_fx * tlx_xy_z_0[j];

                tly_xxy_xz_0[j] = pa_x[j] * tly_xy_xz_0[j] + 0.5 * fl1_fx * tly_y_xz_0[j] + 0.5 * fl1_fx * tly_xy_z_0[j] + 0.5 * fl1_fx * tpz_xy_xz_0[j] + fl1_fx * fl1_fgb * tdz_xy_xz_0[j];

                tlz_xxy_xz_0[j] = pa_x[j] * tlz_xy_xz_0[j] + 0.5 * fl1_fx * tlz_y_xz_0[j] + 0.5 * fl1_fx * tlz_xy_z_0[j] - 0.5 * fl1_fx * tpy_xy_xz_0[j] - fl1_fx * fl1_fgb * tdy_xy_xz_0[j];

                tlx_xxy_yy_0[j] = pa_x[j] * tlx_xy_yy_0[j] + 0.5 * fl1_fx * tlx_y_yy_0[j];

                tly_xxy_yy_0[j] = pa_x[j] * tly_xy_yy_0[j] + 0.5 * fl1_fx * tly_y_yy_0[j] + 0.5 * fl1_fx * tpz_xy_yy_0[j] + fl1_fx * fl1_fgb * tdz_xy_yy_0[j];

                tlz_xxy_yy_0[j] = pa_x[j] * tlz_xy_yy_0[j] + 0.5 * fl1_fx * tlz_y_yy_0[j] - 0.5 * fl1_fx * tpy_xy_yy_0[j] - fl1_fx * fl1_fgb * tdy_xy_yy_0[j];

                tlx_xxy_yz_0[j] = pa_x[j] * tlx_xy_yz_0[j] + 0.5 * fl1_fx * tlx_y_yz_0[j];

                tly_xxy_yz_0[j] = pa_x[j] * tly_xy_yz_0[j] + 0.5 * fl1_fx * tly_y_yz_0[j] + 0.5 * fl1_fx * tpz_xy_yz_0[j] + fl1_fx * fl1_fgb * tdz_xy_yz_0[j];

                tlz_xxy_yz_0[j] = pa_x[j] * tlz_xy_yz_0[j] + 0.5 * fl1_fx * tlz_y_yz_0[j] - 0.5 * fl1_fx * tpy_xy_yz_0[j] - fl1_fx * fl1_fgb * tdy_xy_yz_0[j];

                tlx_xxy_zz_0[j] = pa_x[j] * tlx_xy_zz_0[j] + 0.5 * fl1_fx * tlx_y_zz_0[j];

                tly_xxy_zz_0[j] = pa_x[j] * tly_xy_zz_0[j] + 0.5 * fl1_fx * tly_y_zz_0[j] + 0.5 * fl1_fx * tpz_xy_zz_0[j] + fl1_fx * fl1_fgb * tdz_xy_zz_0[j];

                tlz_xxy_zz_0[j] = pa_x[j] * tlz_xy_zz_0[j] + 0.5 * fl1_fx * tlz_y_zz_0[j] - 0.5 * fl1_fx * tpy_xy_zz_0[j] - fl1_fx * fl1_fgb * tdy_xy_zz_0[j];

                tlx_xxz_xx_0[j] = pa_x[j] * tlx_xz_xx_0[j] + 0.5 * fl1_fx * tlx_z_xx_0[j] + fl1_fx * tlx_xz_x_0[j];

                tly_xxz_xx_0[j] = pa_x[j] * tly_xz_xx_0[j] + 0.5 * fl1_fx * tly_z_xx_0[j] + fl1_fx * tly_xz_x_0[j] + 0.5 * fl1_fx * tpz_xz_xx_0[j] + fl1_fx * fl1_fgb * tdz_xz_xx_0[j];

                tlz_xxz_xx_0[j] = pa_x[j] * tlz_xz_xx_0[j] + 0.5 * fl1_fx * tlz_z_xx_0[j] + fl1_fx * tlz_xz_x_0[j] - 0.5 * fl1_fx * tpy_xz_xx_0[j] - fl1_fx * fl1_fgb * tdy_xz_xx_0[j];

                tlx_xxz_xy_0[j] = pa_x[j] * tlx_xz_xy_0[j] + 0.5 * fl1_fx * tlx_z_xy_0[j] + 0.5 * fl1_fx * tlx_xz_y_0[j];

                tly_xxz_xy_0[j] = pa_x[j] * tly_xz_xy_0[j] + 0.5 * fl1_fx * tly_z_xy_0[j] + 0.5 * fl1_fx * tly_xz_y_0[j] + 0.5 * fl1_fx * tpz_xz_xy_0[j] + fl1_fx * fl1_fgb * tdz_xz_xy_0[j];

                tlz_xxz_xy_0[j] = pa_x[j] * tlz_xz_xy_0[j] + 0.5 * fl1_fx * tlz_z_xy_0[j] + 0.5 * fl1_fx * tlz_xz_y_0[j] - 0.5 * fl1_fx * tpy_xz_xy_0[j] - fl1_fx * fl1_fgb * tdy_xz_xy_0[j];

                tlx_xxz_xz_0[j] = pa_x[j] * tlx_xz_xz_0[j] + 0.5 * fl1_fx * tlx_z_xz_0[j] + 0.5 * fl1_fx * tlx_xz_z_0[j];

                tly_xxz_xz_0[j] = pa_x[j] * tly_xz_xz_0[j] + 0.5 * fl1_fx * tly_z_xz_0[j] + 0.5 * fl1_fx * tly_xz_z_0[j] + 0.5 * fl1_fx * tpz_xz_xz_0[j] + fl1_fx * fl1_fgb * tdz_xz_xz_0[j];

                tlz_xxz_xz_0[j] = pa_x[j] * tlz_xz_xz_0[j] + 0.5 * fl1_fx * tlz_z_xz_0[j] + 0.5 * fl1_fx * tlz_xz_z_0[j] - 0.5 * fl1_fx * tpy_xz_xz_0[j] - fl1_fx * fl1_fgb * tdy_xz_xz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForFD_45_90(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 15); 

            auto tly_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tlz_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tlx_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 16); 

            auto tly_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tlz_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tlx_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 17); 

            auto tly_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tlz_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tlx_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 18); 

            auto tly_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tlz_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tlx_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 19); 

            auto tly_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tlz_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tlx_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 20); 

            auto tly_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tlz_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tlx_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 21); 

            auto tly_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tlz_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tlx_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 22); 

            auto tly_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tlz_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tlx_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 23); 

            auto tly_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tlz_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tlx_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 24); 

            auto tly_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tlz_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tlx_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 25); 

            auto tly_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tlz_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tlx_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 26); 

            auto tly_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tlz_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tlx_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 27); 

            auto tly_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tlz_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tlx_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 28); 

            auto tly_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tlz_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tlx_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 29); 

            auto tly_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tlz_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tlx_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 15); 

            auto tly_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 16); 

            auto tly_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 17); 

            auto tly_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tlx_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 9); 

            auto tly_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 10); 

            auto tly_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 10); 

            auto tlx_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 11); 

            auto tly_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 12); 

            auto tly_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 13); 

            auto tly_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 14); 

            auto tly_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 14); 

            auto tpy_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tpz_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tpy_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tpz_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tpy_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tpz_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tdy_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tdz_xz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tdy_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tdz_xz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tdy_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tdz_xz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tdy_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tdy_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tdy_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tdy_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tdy_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tdy_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tdy_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tdy_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tdy_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tdy_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tdy_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tdy_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tdz_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            // set up pointers to integrals

            auto tlx_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 15); 

            auto tly_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 15); 

            auto tlz_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 15); 

            auto tlx_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 16); 

            auto tly_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 16); 

            auto tlz_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 16); 

            auto tlx_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 17); 

            auto tly_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 17); 

            auto tlz_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 17); 

            auto tlx_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 18); 

            auto tly_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 18); 

            auto tlz_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 18); 

            auto tlx_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 19); 

            auto tly_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 19); 

            auto tlz_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 19); 

            auto tlx_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 20); 

            auto tly_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 20); 

            auto tlz_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 20); 

            auto tlx_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 21); 

            auto tly_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 21); 

            auto tlz_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 21); 

            auto tlx_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 22); 

            auto tly_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 22); 

            auto tlz_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 22); 

            auto tlx_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 23); 

            auto tly_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 23); 

            auto tlz_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 23); 

            auto tlx_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 24); 

            auto tly_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 24); 

            auto tlz_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 24); 

            auto tlx_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 25); 

            auto tly_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 25); 

            auto tlz_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 25); 

            auto tlx_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 26); 

            auto tly_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 26); 

            auto tlz_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 26); 

            auto tlx_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 27); 

            auto tly_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 27); 

            auto tlz_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 27); 

            auto tlx_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 28); 

            auto tly_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 28); 

            auto tlz_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 28); 

            auto tlx_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 29); 

            auto tly_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 29); 

            auto tlz_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_xz_yy_0, tdy_xz_yz_0, tdy_xz_zz_0, tdy_yy_xx_0, \
                                     tdy_yy_xy_0, tdy_yy_xz_0, tdy_yy_yy_0, tdy_yy_yz_0, tdy_yy_zz_0, tdy_yz_xx_0, \
                                     tdy_yz_xy_0, tdy_yz_xz_0, tdy_yz_yy_0, tdy_yz_yz_0, tdy_yz_zz_0, tdz_xz_yy_0, \
                                     tdz_xz_yz_0, tdz_xz_zz_0, tdz_yy_xx_0, tdz_yy_xy_0, tdz_yy_xz_0, tdz_yy_yy_0, \
                                     tdz_yy_yz_0, tdz_yy_zz_0, tdz_yz_xx_0, tdz_yz_xy_0, tdz_yz_xz_0, tdz_yz_yy_0, \
                                     tdz_yz_yz_0, tdz_yz_zz_0, tlx_xxz_yy_0, tlx_xxz_yz_0, tlx_xxz_zz_0, tlx_xyy_xx_0, \
                                     tlx_xyy_xy_0, tlx_xyy_xz_0, tlx_xyy_yy_0, tlx_xyy_yz_0, tlx_xyy_zz_0, tlx_xyz_xx_0, \
                                     tlx_xyz_xy_0, tlx_xyz_xz_0, tlx_xyz_yy_0, tlx_xyz_yz_0, tlx_xyz_zz_0, tlx_xz_yy_0, \
                                     tlx_xz_yz_0, tlx_xz_zz_0, tlx_yy_x_0, tlx_yy_xx_0, tlx_yy_xy_0, tlx_yy_xz_0, \
                                     tlx_yy_y_0, tlx_yy_yy_0, tlx_yy_yz_0, tlx_yy_z_0, tlx_yy_zz_0, tlx_yz_x_0, \
                                     tlx_yz_xx_0, tlx_yz_xy_0, tlx_yz_xz_0, tlx_yz_y_0, tlx_yz_yy_0, tlx_yz_yz_0, \
                                     tlx_yz_z_0, tlx_yz_zz_0, tlx_z_yy_0, tlx_z_yz_0, tlx_z_zz_0, tly_xxz_yy_0, \
                                     tly_xxz_yz_0, tly_xxz_zz_0, tly_xyy_xx_0, tly_xyy_xy_0, tly_xyy_xz_0, tly_xyy_yy_0, \
                                     tly_xyy_yz_0, tly_xyy_zz_0, tly_xyz_xx_0, tly_xyz_xy_0, tly_xyz_xz_0, tly_xyz_yy_0, \
                                     tly_xyz_yz_0, tly_xyz_zz_0, tly_xz_yy_0, tly_xz_yz_0, tly_xz_zz_0, tly_yy_x_0, \
                                     tly_yy_xx_0, tly_yy_xy_0, tly_yy_xz_0, tly_yy_y_0, tly_yy_yy_0, tly_yy_yz_0, \
                                     tly_yy_z_0, tly_yy_zz_0, tly_yz_x_0, tly_yz_xx_0, tly_yz_xy_0, tly_yz_xz_0, \
                                     tly_yz_y_0, tly_yz_yy_0, tly_yz_yz_0, tly_yz_z_0, tly_yz_zz_0, tly_z_yy_0, \
                                     tly_z_yz_0, tly_z_zz_0, tlz_xxz_yy_0, tlz_xxz_yz_0, tlz_xxz_zz_0, tlz_xyy_xx_0, \
                                     tlz_xyy_xy_0, tlz_xyy_xz_0, tlz_xyy_yy_0, tlz_xyy_yz_0, tlz_xyy_zz_0, tlz_xyz_xx_0, \
                                     tlz_xyz_xy_0, tlz_xyz_xz_0, tlz_xyz_yy_0, tlz_xyz_yz_0, tlz_xyz_zz_0, tlz_xz_yy_0, \
                                     tlz_xz_yz_0, tlz_xz_zz_0, tlz_yy_x_0, tlz_yy_xx_0, tlz_yy_xy_0, tlz_yy_xz_0, \
                                     tlz_yy_y_0, tlz_yy_yy_0, tlz_yy_yz_0, tlz_yy_z_0, tlz_yy_zz_0, tlz_yz_x_0, \
                                     tlz_yz_xx_0, tlz_yz_xy_0, tlz_yz_xz_0, tlz_yz_y_0, tlz_yz_yy_0, tlz_yz_yz_0, \
                                     tlz_yz_z_0, tlz_yz_zz_0, tlz_z_yy_0, tlz_z_yz_0, tlz_z_zz_0, tpy_xz_yy_0, \
                                     tpy_xz_yz_0, tpy_xz_zz_0, tpy_yy_xx_0, tpy_yy_xy_0, tpy_yy_xz_0, tpy_yy_yy_0, \
                                     tpy_yy_yz_0, tpy_yy_zz_0, tpy_yz_xx_0, tpy_yz_xy_0, tpy_yz_xz_0, tpy_yz_yy_0, \
                                     tpy_yz_yz_0, tpy_yz_zz_0, tpz_xz_yy_0, tpz_xz_yz_0, tpz_xz_zz_0, tpz_yy_xx_0, \
                                     tpz_yy_xy_0, tpz_yy_xz_0, tpz_yy_yy_0, tpz_yy_yz_0, tpz_yy_zz_0, tpz_yz_xx_0, \
                                     tpz_yz_xy_0, tpz_yz_xz_0, tpz_yz_yy_0, tpz_yz_yz_0, tpz_yz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xxz_yy_0[j] = pa_x[j] * tlx_xz_yy_0[j] + 0.5 * fl1_fx * tlx_z_yy_0[j];

                tly_xxz_yy_0[j] = pa_x[j] * tly_xz_yy_0[j] + 0.5 * fl1_fx * tly_z_yy_0[j] + 0.5 * fl1_fx * tpz_xz_yy_0[j] + fl1_fx * fl1_fgb * tdz_xz_yy_0[j];

                tlz_xxz_yy_0[j] = pa_x[j] * tlz_xz_yy_0[j] + 0.5 * fl1_fx * tlz_z_yy_0[j] - 0.5 * fl1_fx * tpy_xz_yy_0[j] - fl1_fx * fl1_fgb * tdy_xz_yy_0[j];

                tlx_xxz_yz_0[j] = pa_x[j] * tlx_xz_yz_0[j] + 0.5 * fl1_fx * tlx_z_yz_0[j];

                tly_xxz_yz_0[j] = pa_x[j] * tly_xz_yz_0[j] + 0.5 * fl1_fx * tly_z_yz_0[j] + 0.5 * fl1_fx * tpz_xz_yz_0[j] + fl1_fx * fl1_fgb * tdz_xz_yz_0[j];

                tlz_xxz_yz_0[j] = pa_x[j] * tlz_xz_yz_0[j] + 0.5 * fl1_fx * tlz_z_yz_0[j] - 0.5 * fl1_fx * tpy_xz_yz_0[j] - fl1_fx * fl1_fgb * tdy_xz_yz_0[j];

                tlx_xxz_zz_0[j] = pa_x[j] * tlx_xz_zz_0[j] + 0.5 * fl1_fx * tlx_z_zz_0[j];

                tly_xxz_zz_0[j] = pa_x[j] * tly_xz_zz_0[j] + 0.5 * fl1_fx * tly_z_zz_0[j] + 0.5 * fl1_fx * tpz_xz_zz_0[j] + fl1_fx * fl1_fgb * tdz_xz_zz_0[j];

                tlz_xxz_zz_0[j] = pa_x[j] * tlz_xz_zz_0[j] + 0.5 * fl1_fx * tlz_z_zz_0[j] - 0.5 * fl1_fx * tpy_xz_zz_0[j] - fl1_fx * fl1_fgb * tdy_xz_zz_0[j];

                tlx_xyy_xx_0[j] = pa_x[j] * tlx_yy_xx_0[j] + fl1_fx * tlx_yy_x_0[j];

                tly_xyy_xx_0[j] = pa_x[j] * tly_yy_xx_0[j] + fl1_fx * tly_yy_x_0[j] + 0.5 * fl1_fx * tpz_yy_xx_0[j] + fl1_fx * fl1_fgb * tdz_yy_xx_0[j];

                tlz_xyy_xx_0[j] = pa_x[j] * tlz_yy_xx_0[j] + fl1_fx * tlz_yy_x_0[j] - 0.5 * fl1_fx * tpy_yy_xx_0[j] - fl1_fx * fl1_fgb * tdy_yy_xx_0[j];

                tlx_xyy_xy_0[j] = pa_x[j] * tlx_yy_xy_0[j] + 0.5 * fl1_fx * tlx_yy_y_0[j];

                tly_xyy_xy_0[j] = pa_x[j] * tly_yy_xy_0[j] + 0.5 * fl1_fx * tly_yy_y_0[j] + 0.5 * fl1_fx * tpz_yy_xy_0[j] + fl1_fx * fl1_fgb * tdz_yy_xy_0[j];

                tlz_xyy_xy_0[j] = pa_x[j] * tlz_yy_xy_0[j] + 0.5 * fl1_fx * tlz_yy_y_0[j] - 0.5 * fl1_fx * tpy_yy_xy_0[j] - fl1_fx * fl1_fgb * tdy_yy_xy_0[j];

                tlx_xyy_xz_0[j] = pa_x[j] * tlx_yy_xz_0[j] + 0.5 * fl1_fx * tlx_yy_z_0[j];

                tly_xyy_xz_0[j] = pa_x[j] * tly_yy_xz_0[j] + 0.5 * fl1_fx * tly_yy_z_0[j] + 0.5 * fl1_fx * tpz_yy_xz_0[j] + fl1_fx * fl1_fgb * tdz_yy_xz_0[j];

                tlz_xyy_xz_0[j] = pa_x[j] * tlz_yy_xz_0[j] + 0.5 * fl1_fx * tlz_yy_z_0[j] - 0.5 * fl1_fx * tpy_yy_xz_0[j] - fl1_fx * fl1_fgb * tdy_yy_xz_0[j];

                tlx_xyy_yy_0[j] = pa_x[j] * tlx_yy_yy_0[j];

                tly_xyy_yy_0[j] = pa_x[j] * tly_yy_yy_0[j] + 0.5 * fl1_fx * tpz_yy_yy_0[j] + fl1_fx * fl1_fgb * tdz_yy_yy_0[j];

                tlz_xyy_yy_0[j] = pa_x[j] * tlz_yy_yy_0[j] - 0.5 * fl1_fx * tpy_yy_yy_0[j] - fl1_fx * fl1_fgb * tdy_yy_yy_0[j];

                tlx_xyy_yz_0[j] = pa_x[j] * tlx_yy_yz_0[j];

                tly_xyy_yz_0[j] = pa_x[j] * tly_yy_yz_0[j] + 0.5 * fl1_fx * tpz_yy_yz_0[j] + fl1_fx * fl1_fgb * tdz_yy_yz_0[j];

                tlz_xyy_yz_0[j] = pa_x[j] * tlz_yy_yz_0[j] - 0.5 * fl1_fx * tpy_yy_yz_0[j] - fl1_fx * fl1_fgb * tdy_yy_yz_0[j];

                tlx_xyy_zz_0[j] = pa_x[j] * tlx_yy_zz_0[j];

                tly_xyy_zz_0[j] = pa_x[j] * tly_yy_zz_0[j] + 0.5 * fl1_fx * tpz_yy_zz_0[j] + fl1_fx * fl1_fgb * tdz_yy_zz_0[j];

                tlz_xyy_zz_0[j] = pa_x[j] * tlz_yy_zz_0[j] - 0.5 * fl1_fx * tpy_yy_zz_0[j] - fl1_fx * fl1_fgb * tdy_yy_zz_0[j];

                tlx_xyz_xx_0[j] = pa_x[j] * tlx_yz_xx_0[j] + fl1_fx * tlx_yz_x_0[j];

                tly_xyz_xx_0[j] = pa_x[j] * tly_yz_xx_0[j] + fl1_fx * tly_yz_x_0[j] + 0.5 * fl1_fx * tpz_yz_xx_0[j] + fl1_fx * fl1_fgb * tdz_yz_xx_0[j];

                tlz_xyz_xx_0[j] = pa_x[j] * tlz_yz_xx_0[j] + fl1_fx * tlz_yz_x_0[j] - 0.5 * fl1_fx * tpy_yz_xx_0[j] - fl1_fx * fl1_fgb * tdy_yz_xx_0[j];

                tlx_xyz_xy_0[j] = pa_x[j] * tlx_yz_xy_0[j] + 0.5 * fl1_fx * tlx_yz_y_0[j];

                tly_xyz_xy_0[j] = pa_x[j] * tly_yz_xy_0[j] + 0.5 * fl1_fx * tly_yz_y_0[j] + 0.5 * fl1_fx * tpz_yz_xy_0[j] + fl1_fx * fl1_fgb * tdz_yz_xy_0[j];

                tlz_xyz_xy_0[j] = pa_x[j] * tlz_yz_xy_0[j] + 0.5 * fl1_fx * tlz_yz_y_0[j] - 0.5 * fl1_fx * tpy_yz_xy_0[j] - fl1_fx * fl1_fgb * tdy_yz_xy_0[j];

                tlx_xyz_xz_0[j] = pa_x[j] * tlx_yz_xz_0[j] + 0.5 * fl1_fx * tlx_yz_z_0[j];

                tly_xyz_xz_0[j] = pa_x[j] * tly_yz_xz_0[j] + 0.5 * fl1_fx * tly_yz_z_0[j] + 0.5 * fl1_fx * tpz_yz_xz_0[j] + fl1_fx * fl1_fgb * tdz_yz_xz_0[j];

                tlz_xyz_xz_0[j] = pa_x[j] * tlz_yz_xz_0[j] + 0.5 * fl1_fx * tlz_yz_z_0[j] - 0.5 * fl1_fx * tpy_yz_xz_0[j] - fl1_fx * fl1_fgb * tdy_yz_xz_0[j];

                tlx_xyz_yy_0[j] = pa_x[j] * tlx_yz_yy_0[j];

                tly_xyz_yy_0[j] = pa_x[j] * tly_yz_yy_0[j] + 0.5 * fl1_fx * tpz_yz_yy_0[j] + fl1_fx * fl1_fgb * tdz_yz_yy_0[j];

                tlz_xyz_yy_0[j] = pa_x[j] * tlz_yz_yy_0[j] - 0.5 * fl1_fx * tpy_yz_yy_0[j] - fl1_fx * fl1_fgb * tdy_yz_yy_0[j];

                tlx_xyz_yz_0[j] = pa_x[j] * tlx_yz_yz_0[j];

                tly_xyz_yz_0[j] = pa_x[j] * tly_yz_yz_0[j] + 0.5 * fl1_fx * tpz_yz_yz_0[j] + fl1_fx * fl1_fgb * tdz_yz_yz_0[j];

                tlz_xyz_yz_0[j] = pa_x[j] * tlz_yz_yz_0[j] - 0.5 * fl1_fx * tpy_yz_yz_0[j] - fl1_fx * fl1_fgb * tdy_yz_yz_0[j];

                tlx_xyz_zz_0[j] = pa_x[j] * tlx_yz_zz_0[j];

                tly_xyz_zz_0[j] = pa_x[j] * tly_yz_zz_0[j] + 0.5 * fl1_fx * tpz_yz_zz_0[j] + fl1_fx * fl1_fgb * tdz_yz_zz_0[j];

                tlz_xyz_zz_0[j] = pa_x[j] * tlz_yz_zz_0[j] - 0.5 * fl1_fx * tpy_yz_zz_0[j] - fl1_fx * fl1_fgb * tdy_yz_zz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForFD_90_135(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tlx_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 18); 

            auto tly_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tlz_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tlx_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 19); 

            auto tly_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tlz_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tlx_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 20); 

            auto tly_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tlz_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tlx_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 21); 

            auto tly_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tlz_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tlx_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 22); 

            auto tly_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tlz_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tlx_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 23); 

            auto tly_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tlz_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tlx_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 24); 

            auto tly_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tlz_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tlx_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 25); 

            auto tly_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tlz_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tlx_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 26); 

            auto tly_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tlz_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30); 

            auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31); 

            auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32); 

            auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tlx_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 33); 

            auto tly_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tlz_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tlx_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 34); 

            auto tly_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tlz_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tlx_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 35); 

            auto tly_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tlz_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            auto tlx_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 6); 

            auto tly_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 6); 

            auto tlz_y_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 6); 

            auto tlx_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 7); 

            auto tly_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 7); 

            auto tlz_y_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 7); 

            auto tlx_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 8); 

            auto tly_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 8); 

            auto tlz_y_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 8); 

            auto tlx_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 9); 

            auto tly_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_y_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 10); 

            auto tly_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_y_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 10); 

            auto tlx_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 11); 

            auto tly_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_y_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 9); 

            auto tly_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 9); 

            auto tlz_yy_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 9); 

            auto tlx_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 10); 

            auto tly_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 10); 

            auto tlz_yy_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 10); 

            auto tlx_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 11); 

            auto tly_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 11); 

            auto tlz_yy_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 11); 

            auto tlx_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 12); 

            auto tly_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_yz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 15); 

            auto tly_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 16); 

            auto tly_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 17); 

            auto tly_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 17); 

            auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18); 

            auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19); 

            auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20); 

            auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21); 

            auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22); 

            auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23); 

            auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24); 

            auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25); 

            auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26); 

            auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            auto tdx_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 18); 

            auto tdz_yy_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tdx_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 19); 

            auto tdz_yy_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tdx_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 20); 

            auto tdz_yy_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tdx_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 21); 

            auto tdz_yy_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tdx_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 22); 

            auto tdz_yy_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tdx_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 23); 

            auto tdz_yy_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tdx_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 24); 

            auto tdz_yz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tdx_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 25); 

            auto tdz_yz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tdx_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 26); 

            auto tdz_yz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tdy_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tdz_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tdy_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tdz_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tdy_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tdz_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tdy_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tdz_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tdy_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tdz_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tdy_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tdz_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            // set up pointers to integrals

            auto tlx_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 30); 

            auto tly_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 30); 

            auto tlz_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 30); 

            auto tlx_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 31); 

            auto tly_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 31); 

            auto tlz_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 31); 

            auto tlx_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 32); 

            auto tly_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 32); 

            auto tlz_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 32); 

            auto tlx_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 33); 

            auto tly_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 33); 

            auto tlz_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 33); 

            auto tlx_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 34); 

            auto tly_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 34); 

            auto tlz_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 34); 

            auto tlx_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 35); 

            auto tly_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 35); 

            auto tlz_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 35); 

            auto tlx_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 36); 

            auto tly_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tlz_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tlx_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 37); 

            auto tly_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tlz_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tlx_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 38); 

            auto tly_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tlz_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tlx_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 39); 

            auto tly_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tlz_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tlx_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 40); 

            auto tly_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tlz_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tlx_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 41); 

            auto tly_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tlz_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tlx_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 42); 

            auto tly_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tlz_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tlx_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 43); 

            auto tly_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tlz_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tlx_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 44); 

            auto tly_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tlz_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tdx_yy_xx_0, tdx_yy_xy_0, tdx_yy_xz_0, tdx_yy_yy_0, \
                                     tdx_yy_yz_0, tdx_yy_zz_0, tdx_yz_xx_0, tdx_yz_xy_0, tdx_yz_xz_0, tdy_zz_xx_0, \
                                     tdy_zz_xy_0, tdy_zz_xz_0, tdy_zz_yy_0, tdy_zz_yz_0, tdy_zz_zz_0, tdz_yy_xx_0, \
                                     tdz_yy_xy_0, tdz_yy_xz_0, tdz_yy_yy_0, tdz_yy_yz_0, tdz_yy_zz_0, tdz_yz_xx_0, \
                                     tdz_yz_xy_0, tdz_yz_xz_0, tdz_zz_xx_0, tdz_zz_xy_0, tdz_zz_xz_0, tdz_zz_yy_0, \
                                     tdz_zz_yz_0, tdz_zz_zz_0, tlx_xzz_xx_0, tlx_xzz_xy_0, tlx_xzz_xz_0, tlx_xzz_yy_0, \
                                     tlx_xzz_yz_0, tlx_xzz_zz_0, tlx_y_xx_0, tlx_y_xy_0, tlx_y_xz_0, tlx_y_yy_0, \
                                     tlx_y_yz_0, tlx_y_zz_0, tlx_yy_x_0, tlx_yy_xx_0, tlx_yy_xy_0, tlx_yy_xz_0, \
                                     tlx_yy_y_0, tlx_yy_yy_0, tlx_yy_yz_0, tlx_yy_z_0, tlx_yy_zz_0, tlx_yyy_xx_0, \
                                     tlx_yyy_xy_0, tlx_yyy_xz_0, tlx_yyy_yy_0, tlx_yyy_yz_0, tlx_yyy_zz_0, tlx_yyz_xx_0, \
                                     tlx_yyz_xy_0, tlx_yyz_xz_0, tlx_yz_x_0, tlx_yz_xx_0, tlx_yz_xy_0, tlx_yz_xz_0, \
                                     tlx_z_xx_0, tlx_z_xy_0, tlx_z_xz_0, tlx_zz_x_0, tlx_zz_xx_0, tlx_zz_xy_0, \
                                     tlx_zz_xz_0, tlx_zz_y_0, tlx_zz_yy_0, tlx_zz_yz_0, tlx_zz_z_0, tlx_zz_zz_0, \
                                     tly_xzz_xx_0, tly_xzz_xy_0, tly_xzz_xz_0, tly_xzz_yy_0, tly_xzz_yz_0, tly_xzz_zz_0, \
                                     tly_y_xx_0, tly_y_xy_0, tly_y_xz_0, tly_y_yy_0, tly_y_yz_0, tly_y_zz_0, tly_yy_x_0, \
                                     tly_yy_xx_0, tly_yy_xy_0, tly_yy_xz_0, tly_yy_y_0, tly_yy_yy_0, tly_yy_yz_0, \
                                     tly_yy_z_0, tly_yy_zz_0, tly_yyy_xx_0, tly_yyy_xy_0, tly_yyy_xz_0, tly_yyy_yy_0, \
                                     tly_yyy_yz_0, tly_yyy_zz_0, tly_yyz_xx_0, tly_yyz_xy_0, tly_yyz_xz_0, tly_yz_x_0, \
                                     tly_yz_xx_0, tly_yz_xy_0, tly_yz_xz_0, tly_z_xx_0, tly_z_xy_0, tly_z_xz_0, \
                                     tly_zz_x_0, tly_zz_xx_0, tly_zz_xy_0, tly_zz_xz_0, tly_zz_y_0, tly_zz_yy_0, \
                                     tly_zz_yz_0, tly_zz_z_0, tly_zz_zz_0, tlz_xzz_xx_0, tlz_xzz_xy_0, tlz_xzz_xz_0, \
                                     tlz_xzz_yy_0, tlz_xzz_yz_0, tlz_xzz_zz_0, tlz_y_xx_0, tlz_y_xy_0, tlz_y_xz_0, \
                                     tlz_y_yy_0, tlz_y_yz_0, tlz_y_zz_0, tlz_yy_x_0, tlz_yy_xx_0, tlz_yy_xy_0, \
                                     tlz_yy_xz_0, tlz_yy_y_0, tlz_yy_yy_0, tlz_yy_yz_0, tlz_yy_z_0, tlz_yy_zz_0, \
                                     tlz_yyy_xx_0, tlz_yyy_xy_0, tlz_yyy_xz_0, tlz_yyy_yy_0, tlz_yyy_yz_0, tlz_yyy_zz_0, \
                                     tlz_yyz_xx_0, tlz_yyz_xy_0, tlz_yyz_xz_0, tlz_yz_x_0, tlz_yz_xx_0, tlz_yz_xy_0, \
                                     tlz_yz_xz_0, tlz_z_xx_0, tlz_z_xy_0, tlz_z_xz_0, tlz_zz_x_0, tlz_zz_xx_0, \
                                     tlz_zz_xy_0, tlz_zz_xz_0, tlz_zz_y_0, tlz_zz_yy_0, tlz_zz_yz_0, tlz_zz_z_0, \
                                     tlz_zz_zz_0, tpx_yy_xx_0, tpx_yy_xy_0, tpx_yy_xz_0, tpx_yy_yy_0, tpx_yy_yz_0, \
                                     tpx_yy_zz_0, tpx_yz_xx_0, tpx_yz_xy_0, tpx_yz_xz_0, tpy_zz_xx_0, tpy_zz_xy_0, \
                                     tpy_zz_xz_0, tpy_zz_yy_0, tpy_zz_yz_0, tpy_zz_zz_0, tpz_yy_xx_0, tpz_yy_xy_0, \
                                     tpz_yy_xz_0, tpz_yy_yy_0, tpz_yy_yz_0, tpz_yy_zz_0, tpz_yz_xx_0, tpz_yz_xy_0, \
                                     tpz_yz_xz_0, tpz_zz_xx_0, tpz_zz_xy_0, tpz_zz_xz_0, tpz_zz_yy_0, tpz_zz_yz_0, \
                                     tpz_zz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xzz_xx_0[j] = pa_x[j] * tlx_zz_xx_0[j] + fl1_fx * tlx_zz_x_0[j];

                tly_xzz_xx_0[j] = pa_x[j] * tly_zz_xx_0[j] + fl1_fx * tly_zz_x_0[j] + 0.5 * fl1_fx * tpz_zz_xx_0[j] + fl1_fx * fl1_fgb * tdz_zz_xx_0[j];

                tlz_xzz_xx_0[j] = pa_x[j] * tlz_zz_xx_0[j] + fl1_fx * tlz_zz_x_0[j] - 0.5 * fl1_fx * tpy_zz_xx_0[j] - fl1_fx * fl1_fgb * tdy_zz_xx_0[j];

                tlx_xzz_xy_0[j] = pa_x[j] * tlx_zz_xy_0[j] + 0.5 * fl1_fx * tlx_zz_y_0[j];

                tly_xzz_xy_0[j] = pa_x[j] * tly_zz_xy_0[j] + 0.5 * fl1_fx * tly_zz_y_0[j] + 0.5 * fl1_fx * tpz_zz_xy_0[j] + fl1_fx * fl1_fgb * tdz_zz_xy_0[j];

                tlz_xzz_xy_0[j] = pa_x[j] * tlz_zz_xy_0[j] + 0.5 * fl1_fx * tlz_zz_y_0[j] - 0.5 * fl1_fx * tpy_zz_xy_0[j] - fl1_fx * fl1_fgb * tdy_zz_xy_0[j];

                tlx_xzz_xz_0[j] = pa_x[j] * tlx_zz_xz_0[j] + 0.5 * fl1_fx * tlx_zz_z_0[j];

                tly_xzz_xz_0[j] = pa_x[j] * tly_zz_xz_0[j] + 0.5 * fl1_fx * tly_zz_z_0[j] + 0.5 * fl1_fx * tpz_zz_xz_0[j] + fl1_fx * fl1_fgb * tdz_zz_xz_0[j];

                tlz_xzz_xz_0[j] = pa_x[j] * tlz_zz_xz_0[j] + 0.5 * fl1_fx * tlz_zz_z_0[j] - 0.5 * fl1_fx * tpy_zz_xz_0[j] - fl1_fx * fl1_fgb * tdy_zz_xz_0[j];

                tlx_xzz_yy_0[j] = pa_x[j] * tlx_zz_yy_0[j];

                tly_xzz_yy_0[j] = pa_x[j] * tly_zz_yy_0[j] + 0.5 * fl1_fx * tpz_zz_yy_0[j] + fl1_fx * fl1_fgb * tdz_zz_yy_0[j];

                tlz_xzz_yy_0[j] = pa_x[j] * tlz_zz_yy_0[j] - 0.5 * fl1_fx * tpy_zz_yy_0[j] - fl1_fx * fl1_fgb * tdy_zz_yy_0[j];

                tlx_xzz_yz_0[j] = pa_x[j] * tlx_zz_yz_0[j];

                tly_xzz_yz_0[j] = pa_x[j] * tly_zz_yz_0[j] + 0.5 * fl1_fx * tpz_zz_yz_0[j] + fl1_fx * fl1_fgb * tdz_zz_yz_0[j];

                tlz_xzz_yz_0[j] = pa_x[j] * tlz_zz_yz_0[j] - 0.5 * fl1_fx * tpy_zz_yz_0[j] - fl1_fx * fl1_fgb * tdy_zz_yz_0[j];

                tlx_xzz_zz_0[j] = pa_x[j] * tlx_zz_zz_0[j];

                tly_xzz_zz_0[j] = pa_x[j] * tly_zz_zz_0[j] + 0.5 * fl1_fx * tpz_zz_zz_0[j] + fl1_fx * fl1_fgb * tdz_zz_zz_0[j];

                tlz_xzz_zz_0[j] = pa_x[j] * tlz_zz_zz_0[j] - 0.5 * fl1_fx * tpy_zz_zz_0[j] - fl1_fx * fl1_fgb * tdy_zz_zz_0[j];

                tlx_yyy_xx_0[j] = pa_y[j] * tlx_yy_xx_0[j] + fl1_fx * tlx_y_xx_0[j] - 0.5 * fl1_fx * tpz_yy_xx_0[j] - fl1_fx * fl1_fgb * tdz_yy_xx_0[j];

                tly_yyy_xx_0[j] = pa_y[j] * tly_yy_xx_0[j] + fl1_fx * tly_y_xx_0[j];

                tlz_yyy_xx_0[j] = pa_y[j] * tlz_yy_xx_0[j] + fl1_fx * tlz_y_xx_0[j] + 0.5 * fl1_fx * tpx_yy_xx_0[j] + fl1_fx * fl1_fgb * tdx_yy_xx_0[j];

                tlx_yyy_xy_0[j] = pa_y[j] * tlx_yy_xy_0[j] + fl1_fx * tlx_y_xy_0[j] + 0.5 * fl1_fx * tlx_yy_x_0[j] - 0.5 * fl1_fx * tpz_yy_xy_0[j] - fl1_fx * fl1_fgb * tdz_yy_xy_0[j];

                tly_yyy_xy_0[j] = pa_y[j] * tly_yy_xy_0[j] + fl1_fx * tly_y_xy_0[j] + 0.5 * fl1_fx * tly_yy_x_0[j];

                tlz_yyy_xy_0[j] = pa_y[j] * tlz_yy_xy_0[j] + fl1_fx * tlz_y_xy_0[j] + 0.5 * fl1_fx * tlz_yy_x_0[j] + 0.5 * fl1_fx * tpx_yy_xy_0[j] + fl1_fx * fl1_fgb * tdx_yy_xy_0[j];

                tlx_yyy_xz_0[j] = pa_y[j] * tlx_yy_xz_0[j] + fl1_fx * tlx_y_xz_0[j] - 0.5 * fl1_fx * tpz_yy_xz_0[j] - fl1_fx * fl1_fgb * tdz_yy_xz_0[j];

                tly_yyy_xz_0[j] = pa_y[j] * tly_yy_xz_0[j] + fl1_fx * tly_y_xz_0[j];

                tlz_yyy_xz_0[j] = pa_y[j] * tlz_yy_xz_0[j] + fl1_fx * tlz_y_xz_0[j] + 0.5 * fl1_fx * tpx_yy_xz_0[j] + fl1_fx * fl1_fgb * tdx_yy_xz_0[j];

                tlx_yyy_yy_0[j] = pa_y[j] * tlx_yy_yy_0[j] + fl1_fx * tlx_y_yy_0[j] + fl1_fx * tlx_yy_y_0[j] - 0.5 * fl1_fx * tpz_yy_yy_0[j] - fl1_fx * fl1_fgb * tdz_yy_yy_0[j];

                tly_yyy_yy_0[j] = pa_y[j] * tly_yy_yy_0[j] + fl1_fx * tly_y_yy_0[j] + fl1_fx * tly_yy_y_0[j];

                tlz_yyy_yy_0[j] = pa_y[j] * tlz_yy_yy_0[j] + fl1_fx * tlz_y_yy_0[j] + fl1_fx * tlz_yy_y_0[j] + 0.5 * fl1_fx * tpx_yy_yy_0[j] + fl1_fx * fl1_fgb * tdx_yy_yy_0[j];

                tlx_yyy_yz_0[j] = pa_y[j] * tlx_yy_yz_0[j] + fl1_fx * tlx_y_yz_0[j] + 0.5 * fl1_fx * tlx_yy_z_0[j] - 0.5 * fl1_fx * tpz_yy_yz_0[j] - fl1_fx * fl1_fgb * tdz_yy_yz_0[j];

                tly_yyy_yz_0[j] = pa_y[j] * tly_yy_yz_0[j] + fl1_fx * tly_y_yz_0[j] + 0.5 * fl1_fx * tly_yy_z_0[j];

                tlz_yyy_yz_0[j] = pa_y[j] * tlz_yy_yz_0[j] + fl1_fx * tlz_y_yz_0[j] + 0.5 * fl1_fx * tlz_yy_z_0[j] + 0.5 * fl1_fx * tpx_yy_yz_0[j] + fl1_fx * fl1_fgb * tdx_yy_yz_0[j];

                tlx_yyy_zz_0[j] = pa_y[j] * tlx_yy_zz_0[j] + fl1_fx * tlx_y_zz_0[j] - 0.5 * fl1_fx * tpz_yy_zz_0[j] - fl1_fx * fl1_fgb * tdz_yy_zz_0[j];

                tly_yyy_zz_0[j] = pa_y[j] * tly_yy_zz_0[j] + fl1_fx * tly_y_zz_0[j];

                tlz_yyy_zz_0[j] = pa_y[j] * tlz_yy_zz_0[j] + fl1_fx * tlz_y_zz_0[j] + 0.5 * fl1_fx * tpx_yy_zz_0[j] + fl1_fx * fl1_fgb * tdx_yy_zz_0[j];

                tlx_yyz_xx_0[j] = pa_y[j] * tlx_yz_xx_0[j] + 0.5 * fl1_fx * tlx_z_xx_0[j] - 0.5 * fl1_fx * tpz_yz_xx_0[j] - fl1_fx * fl1_fgb * tdz_yz_xx_0[j];

                tly_yyz_xx_0[j] = pa_y[j] * tly_yz_xx_0[j] + 0.5 * fl1_fx * tly_z_xx_0[j];

                tlz_yyz_xx_0[j] = pa_y[j] * tlz_yz_xx_0[j] + 0.5 * fl1_fx * tlz_z_xx_0[j] + 0.5 * fl1_fx * tpx_yz_xx_0[j] + fl1_fx * fl1_fgb * tdx_yz_xx_0[j];

                tlx_yyz_xy_0[j] = pa_y[j] * tlx_yz_xy_0[j] + 0.5 * fl1_fx * tlx_z_xy_0[j] + 0.5 * fl1_fx * tlx_yz_x_0[j] - 0.5 * fl1_fx * tpz_yz_xy_0[j] - fl1_fx * fl1_fgb * tdz_yz_xy_0[j];

                tly_yyz_xy_0[j] = pa_y[j] * tly_yz_xy_0[j] + 0.5 * fl1_fx * tly_z_xy_0[j] + 0.5 * fl1_fx * tly_yz_x_0[j];

                tlz_yyz_xy_0[j] = pa_y[j] * tlz_yz_xy_0[j] + 0.5 * fl1_fx * tlz_z_xy_0[j] + 0.5 * fl1_fx * tlz_yz_x_0[j] + 0.5 * fl1_fx * tpx_yz_xy_0[j] + fl1_fx * fl1_fgb * tdx_yz_xy_0[j];

                tlx_yyz_xz_0[j] = pa_y[j] * tlx_yz_xz_0[j] + 0.5 * fl1_fx * tlx_z_xz_0[j] - 0.5 * fl1_fx * tpz_yz_xz_0[j] - fl1_fx * fl1_fgb * tdz_yz_xz_0[j];

                tly_yyz_xz_0[j] = pa_y[j] * tly_yz_xz_0[j] + 0.5 * fl1_fx * tly_z_xz_0[j];

                tlz_yyz_xz_0[j] = pa_y[j] * tlz_yz_xz_0[j] + 0.5 * fl1_fx * tlz_z_xz_0[j] + 0.5 * fl1_fx * tpx_yz_xz_0[j] + fl1_fx * fl1_fgb * tdx_yz_xz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForFD_135_180(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_3_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tlx_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 27); 

            auto tly_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tlz_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tlx_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 28); 

            auto tly_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tlz_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tlx_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 29); 

            auto tly_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tlz_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30); 

            auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31); 

            auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32); 

            auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tlx_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 33); 

            auto tly_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tlz_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tlx_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 34); 

            auto tly_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tlz_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tlx_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 35); 

            auto tly_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tlz_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            auto tlx_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 12); 

            auto tly_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 12); 

            auto tlz_z_xx_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 12); 

            auto tlx_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 13); 

            auto tly_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_z_xy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 14); 

            auto tly_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_z_xz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 15); 

            auto tly_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_z_yy_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 16); 

            auto tly_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_z_yz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * idx + 17); 

            auto tly_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_z_zz_0 = primBuffer.data(pidx_l_1_2_m0 + 36 * bdim + 18 * idx + 17); 

            auto tlx_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 13); 

            auto tly_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 13); 

            auto tlz_yz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 13); 

            auto tlx_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 14); 

            auto tly_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 14); 

            auto tlz_yz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 14); 

            auto tlx_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 15); 

            auto tly_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 15); 

            auto tlz_zz_x_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 15); 

            auto tlx_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 16); 

            auto tly_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 16); 

            auto tlz_zz_y_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 16); 

            auto tlx_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * idx + 17); 

            auto tly_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 18 * bdim + 18 * idx + 17); 

            auto tlz_zz_z_0 = primBuffer.data(pidx_l_2_1_m0 + 36 * bdim + 18 * idx + 17); 

            auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27); 

            auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28); 

            auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29); 

            auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30); 

            auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31); 

            auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32); 

            auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33); 

            auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34); 

            auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35); 

            auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            auto tdx_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 27); 

            auto tdz_yz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tdx_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 28); 

            auto tdz_yz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tdx_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 29); 

            auto tdz_yz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tdx_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 30); 

            auto tdy_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tdz_zz_xx_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tdx_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 31); 

            auto tdy_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tdz_zz_xy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tdx_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 32); 

            auto tdy_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tdz_zz_xz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tdx_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 33); 

            auto tdy_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tdz_zz_yy_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tdx_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 34); 

            auto tdy_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tdz_zz_yz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tdx_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * idx + 35); 

            auto tdy_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tdz_zz_zz_0 = primBuffer.data(pidx_d_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            // set up pointers to integrals

            auto tlx_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 45); 

            auto tly_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tlz_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tlx_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 46); 

            auto tly_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tlz_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tlx_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 47); 

            auto tly_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tlz_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tlx_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 48); 

            auto tly_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tlz_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tlx_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 49); 

            auto tly_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tlz_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tlx_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 50); 

            auto tly_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tlz_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tlx_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 51); 

            auto tly_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tlz_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tlx_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 52); 

            auto tly_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tlz_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tlx_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 53); 

            auto tly_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 53); 

            auto tlz_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tlx_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 54); 

            auto tly_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tlz_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tlx_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 55); 

            auto tly_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tlz_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tlx_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 56); 

            auto tly_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tlz_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tlx_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 57); 

            auto tly_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tlz_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tlx_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 58); 

            auto tly_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tlz_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tlx_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 59); 

            auto tly_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tlz_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_yz_yy_0, tdx_yz_yz_0, tdx_yz_zz_0, tdx_zz_xx_0, \
                                     tdx_zz_xy_0, tdx_zz_xz_0, tdx_zz_yy_0, tdx_zz_yz_0, tdx_zz_zz_0, tdy_zz_xx_0, \
                                     tdy_zz_xy_0, tdy_zz_xz_0, tdy_zz_yy_0, tdy_zz_yz_0, tdy_zz_zz_0, tdz_yz_yy_0, \
                                     tdz_yz_yz_0, tdz_yz_zz_0, tdz_zz_xx_0, tdz_zz_xy_0, tdz_zz_xz_0, tdz_zz_yy_0, \
                                     tdz_zz_yz_0, tdz_zz_zz_0, tlx_yyz_yy_0, tlx_yyz_yz_0, tlx_yyz_zz_0, tlx_yz_y_0, \
                                     tlx_yz_yy_0, tlx_yz_yz_0, tlx_yz_z_0, tlx_yz_zz_0, tlx_yzz_xx_0, tlx_yzz_xy_0, \
                                     tlx_yzz_xz_0, tlx_yzz_yy_0, tlx_yzz_yz_0, tlx_yzz_zz_0, tlx_z_xx_0, tlx_z_xy_0, \
                                     tlx_z_xz_0, tlx_z_yy_0, tlx_z_yz_0, tlx_z_zz_0, tlx_zz_x_0, tlx_zz_xx_0, \
                                     tlx_zz_xy_0, tlx_zz_xz_0, tlx_zz_y_0, tlx_zz_yy_0, tlx_zz_yz_0, tlx_zz_z_0, \
                                     tlx_zz_zz_0, tlx_zzz_xx_0, tlx_zzz_xy_0, tlx_zzz_xz_0, tlx_zzz_yy_0, tlx_zzz_yz_0, \
                                     tlx_zzz_zz_0, tly_yyz_yy_0, tly_yyz_yz_0, tly_yyz_zz_0, tly_yz_y_0, tly_yz_yy_0, \
                                     tly_yz_yz_0, tly_yz_z_0, tly_yz_zz_0, tly_yzz_xx_0, tly_yzz_xy_0, tly_yzz_xz_0, \
                                     tly_yzz_yy_0, tly_yzz_yz_0, tly_yzz_zz_0, tly_z_xx_0, tly_z_xy_0, tly_z_xz_0, \
                                     tly_z_yy_0, tly_z_yz_0, tly_z_zz_0, tly_zz_x_0, tly_zz_xx_0, tly_zz_xy_0, \
                                     tly_zz_xz_0, tly_zz_y_0, tly_zz_yy_0, tly_zz_yz_0, tly_zz_z_0, tly_zz_zz_0, \
                                     tly_zzz_xx_0, tly_zzz_xy_0, tly_zzz_xz_0, tly_zzz_yy_0, tly_zzz_yz_0, tly_zzz_zz_0, \
                                     tlz_yyz_yy_0, tlz_yyz_yz_0, tlz_yyz_zz_0, tlz_yz_y_0, tlz_yz_yy_0, tlz_yz_yz_0, \
                                     tlz_yz_z_0, tlz_yz_zz_0, tlz_yzz_xx_0, tlz_yzz_xy_0, tlz_yzz_xz_0, tlz_yzz_yy_0, \
                                     tlz_yzz_yz_0, tlz_yzz_zz_0, tlz_z_xx_0, tlz_z_xy_0, tlz_z_xz_0, tlz_z_yy_0, \
                                     tlz_z_yz_0, tlz_z_zz_0, tlz_zz_x_0, tlz_zz_xx_0, tlz_zz_xy_0, tlz_zz_xz_0, \
                                     tlz_zz_y_0, tlz_zz_yy_0, tlz_zz_yz_0, tlz_zz_z_0, tlz_zz_zz_0, tlz_zzz_xx_0, \
                                     tlz_zzz_xy_0, tlz_zzz_xz_0, tlz_zzz_yy_0, tlz_zzz_yz_0, tlz_zzz_zz_0, tpx_yz_yy_0, \
                                     tpx_yz_yz_0, tpx_yz_zz_0, tpx_zz_xx_0, tpx_zz_xy_0, tpx_zz_xz_0, tpx_zz_yy_0, \
                                     tpx_zz_yz_0, tpx_zz_zz_0, tpy_zz_xx_0, tpy_zz_xy_0, tpy_zz_xz_0, tpy_zz_yy_0, \
                                     tpy_zz_yz_0, tpy_zz_zz_0, tpz_yz_yy_0, tpz_yz_yz_0, tpz_yz_zz_0, tpz_zz_xx_0, \
                                     tpz_zz_xy_0, tpz_zz_xz_0, tpz_zz_yy_0, tpz_zz_yz_0, tpz_zz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yyz_yy_0[j] = pa_y[j] * tlx_yz_yy_0[j] + 0.5 * fl1_fx * tlx_z_yy_0[j] + fl1_fx * tlx_yz_y_0[j] - 0.5 * fl1_fx * tpz_yz_yy_0[j] - fl1_fx * fl1_fgb * tdz_yz_yy_0[j];

                tly_yyz_yy_0[j] = pa_y[j] * tly_yz_yy_0[j] + 0.5 * fl1_fx * tly_z_yy_0[j] + fl1_fx * tly_yz_y_0[j];

                tlz_yyz_yy_0[j] = pa_y[j] * tlz_yz_yy_0[j] + 0.5 * fl1_fx * tlz_z_yy_0[j] + fl1_fx * tlz_yz_y_0[j] + 0.5 * fl1_fx * tpx_yz_yy_0[j] + fl1_fx * fl1_fgb * tdx_yz_yy_0[j];

                tlx_yyz_yz_0[j] = pa_y[j] * tlx_yz_yz_0[j] + 0.5 * fl1_fx * tlx_z_yz_0[j] + 0.5 * fl1_fx * tlx_yz_z_0[j] - 0.5 * fl1_fx * tpz_yz_yz_0[j] - fl1_fx * fl1_fgb * tdz_yz_yz_0[j];

                tly_yyz_yz_0[j] = pa_y[j] * tly_yz_yz_0[j] + 0.5 * fl1_fx * tly_z_yz_0[j] + 0.5 * fl1_fx * tly_yz_z_0[j];

                tlz_yyz_yz_0[j] = pa_y[j] * tlz_yz_yz_0[j] + 0.5 * fl1_fx * tlz_z_yz_0[j] + 0.5 * fl1_fx * tlz_yz_z_0[j] + 0.5 * fl1_fx * tpx_yz_yz_0[j] + fl1_fx * fl1_fgb * tdx_yz_yz_0[j];

                tlx_yyz_zz_0[j] = pa_y[j] * tlx_yz_zz_0[j] + 0.5 * fl1_fx * tlx_z_zz_0[j] - 0.5 * fl1_fx * tpz_yz_zz_0[j] - fl1_fx * fl1_fgb * tdz_yz_zz_0[j];

                tly_yyz_zz_0[j] = pa_y[j] * tly_yz_zz_0[j] + 0.5 * fl1_fx * tly_z_zz_0[j];

                tlz_yyz_zz_0[j] = pa_y[j] * tlz_yz_zz_0[j] + 0.5 * fl1_fx * tlz_z_zz_0[j] + 0.5 * fl1_fx * tpx_yz_zz_0[j] + fl1_fx * fl1_fgb * tdx_yz_zz_0[j];

                tlx_yzz_xx_0[j] = pa_y[j] * tlx_zz_xx_0[j] - 0.5 * fl1_fx * tpz_zz_xx_0[j] - fl1_fx * fl1_fgb * tdz_zz_xx_0[j];

                tly_yzz_xx_0[j] = pa_y[j] * tly_zz_xx_0[j];

                tlz_yzz_xx_0[j] = pa_y[j] * tlz_zz_xx_0[j] + 0.5 * fl1_fx * tpx_zz_xx_0[j] + fl1_fx * fl1_fgb * tdx_zz_xx_0[j];

                tlx_yzz_xy_0[j] = pa_y[j] * tlx_zz_xy_0[j] + 0.5 * fl1_fx * tlx_zz_x_0[j] - 0.5 * fl1_fx * tpz_zz_xy_0[j] - fl1_fx * fl1_fgb * tdz_zz_xy_0[j];

                tly_yzz_xy_0[j] = pa_y[j] * tly_zz_xy_0[j] + 0.5 * fl1_fx * tly_zz_x_0[j];

                tlz_yzz_xy_0[j] = pa_y[j] * tlz_zz_xy_0[j] + 0.5 * fl1_fx * tlz_zz_x_0[j] + 0.5 * fl1_fx * tpx_zz_xy_0[j] + fl1_fx * fl1_fgb * tdx_zz_xy_0[j];

                tlx_yzz_xz_0[j] = pa_y[j] * tlx_zz_xz_0[j] - 0.5 * fl1_fx * tpz_zz_xz_0[j] - fl1_fx * fl1_fgb * tdz_zz_xz_0[j];

                tly_yzz_xz_0[j] = pa_y[j] * tly_zz_xz_0[j];

                tlz_yzz_xz_0[j] = pa_y[j] * tlz_zz_xz_0[j] + 0.5 * fl1_fx * tpx_zz_xz_0[j] + fl1_fx * fl1_fgb * tdx_zz_xz_0[j];

                tlx_yzz_yy_0[j] = pa_y[j] * tlx_zz_yy_0[j] + fl1_fx * tlx_zz_y_0[j] - 0.5 * fl1_fx * tpz_zz_yy_0[j] - fl1_fx * fl1_fgb * tdz_zz_yy_0[j];

                tly_yzz_yy_0[j] = pa_y[j] * tly_zz_yy_0[j] + fl1_fx * tly_zz_y_0[j];

                tlz_yzz_yy_0[j] = pa_y[j] * tlz_zz_yy_0[j] + fl1_fx * tlz_zz_y_0[j] + 0.5 * fl1_fx * tpx_zz_yy_0[j] + fl1_fx * fl1_fgb * tdx_zz_yy_0[j];

                tlx_yzz_yz_0[j] = pa_y[j] * tlx_zz_yz_0[j] + 0.5 * fl1_fx * tlx_zz_z_0[j] - 0.5 * fl1_fx * tpz_zz_yz_0[j] - fl1_fx * fl1_fgb * tdz_zz_yz_0[j];

                tly_yzz_yz_0[j] = pa_y[j] * tly_zz_yz_0[j] + 0.5 * fl1_fx * tly_zz_z_0[j];

                tlz_yzz_yz_0[j] = pa_y[j] * tlz_zz_yz_0[j] + 0.5 * fl1_fx * tlz_zz_z_0[j] + 0.5 * fl1_fx * tpx_zz_yz_0[j] + fl1_fx * fl1_fgb * tdx_zz_yz_0[j];

                tlx_yzz_zz_0[j] = pa_y[j] * tlx_zz_zz_0[j] - 0.5 * fl1_fx * tpz_zz_zz_0[j] - fl1_fx * fl1_fgb * tdz_zz_zz_0[j];

                tly_yzz_zz_0[j] = pa_y[j] * tly_zz_zz_0[j];

                tlz_yzz_zz_0[j] = pa_y[j] * tlz_zz_zz_0[j] + 0.5 * fl1_fx * tpx_zz_zz_0[j] + fl1_fx * fl1_fgb * tdx_zz_zz_0[j];

                tlx_zzz_xx_0[j] = pa_z[j] * tlx_zz_xx_0[j] + fl1_fx * tlx_z_xx_0[j] + 0.5 * fl1_fx * tpy_zz_xx_0[j] + fl1_fx * fl1_fgb * tdy_zz_xx_0[j];

                tly_zzz_xx_0[j] = pa_z[j] * tly_zz_xx_0[j] + fl1_fx * tly_z_xx_0[j] - 0.5 * fl1_fx * tpx_zz_xx_0[j] - fl1_fx * fl1_fgb * tdx_zz_xx_0[j];

                tlz_zzz_xx_0[j] = pa_z[j] * tlz_zz_xx_0[j] + fl1_fx * tlz_z_xx_0[j];

                tlx_zzz_xy_0[j] = pa_z[j] * tlx_zz_xy_0[j] + fl1_fx * tlx_z_xy_0[j] + 0.5 * fl1_fx * tpy_zz_xy_0[j] + fl1_fx * fl1_fgb * tdy_zz_xy_0[j];

                tly_zzz_xy_0[j] = pa_z[j] * tly_zz_xy_0[j] + fl1_fx * tly_z_xy_0[j] - 0.5 * fl1_fx * tpx_zz_xy_0[j] - fl1_fx * fl1_fgb * tdx_zz_xy_0[j];

                tlz_zzz_xy_0[j] = pa_z[j] * tlz_zz_xy_0[j] + fl1_fx * tlz_z_xy_0[j];

                tlx_zzz_xz_0[j] = pa_z[j] * tlx_zz_xz_0[j] + fl1_fx * tlx_z_xz_0[j] + 0.5 * fl1_fx * tlx_zz_x_0[j] + 0.5 * fl1_fx * tpy_zz_xz_0[j] + fl1_fx * fl1_fgb * tdy_zz_xz_0[j];

                tly_zzz_xz_0[j] = pa_z[j] * tly_zz_xz_0[j] + fl1_fx * tly_z_xz_0[j] + 0.5 * fl1_fx * tly_zz_x_0[j] - 0.5 * fl1_fx * tpx_zz_xz_0[j] - fl1_fx * fl1_fgb * tdx_zz_xz_0[j];

                tlz_zzz_xz_0[j] = pa_z[j] * tlz_zz_xz_0[j] + fl1_fx * tlz_z_xz_0[j] + 0.5 * fl1_fx * tlz_zz_x_0[j];

                tlx_zzz_yy_0[j] = pa_z[j] * tlx_zz_yy_0[j] + fl1_fx * tlx_z_yy_0[j] + 0.5 * fl1_fx * tpy_zz_yy_0[j] + fl1_fx * fl1_fgb * tdy_zz_yy_0[j];

                tly_zzz_yy_0[j] = pa_z[j] * tly_zz_yy_0[j] + fl1_fx * tly_z_yy_0[j] - 0.5 * fl1_fx * tpx_zz_yy_0[j] - fl1_fx * fl1_fgb * tdx_zz_yy_0[j];

                tlz_zzz_yy_0[j] = pa_z[j] * tlz_zz_yy_0[j] + fl1_fx * tlz_z_yy_0[j];

                tlx_zzz_yz_0[j] = pa_z[j] * tlx_zz_yz_0[j] + fl1_fx * tlx_z_yz_0[j] + 0.5 * fl1_fx * tlx_zz_y_0[j] + 0.5 * fl1_fx * tpy_zz_yz_0[j] + fl1_fx * fl1_fgb * tdy_zz_yz_0[j];

                tly_zzz_yz_0[j] = pa_z[j] * tly_zz_yz_0[j] + fl1_fx * tly_z_yz_0[j] + 0.5 * fl1_fx * tly_zz_y_0[j] - 0.5 * fl1_fx * tpx_zz_yz_0[j] - fl1_fx * fl1_fgb * tdx_zz_yz_0[j];

                tlz_zzz_yz_0[j] = pa_z[j] * tlz_zz_yz_0[j] + fl1_fx * tlz_z_yz_0[j] + 0.5 * fl1_fx * tlz_zz_y_0[j];

                tlx_zzz_zz_0[j] = pa_z[j] * tlx_zz_zz_0[j] + fl1_fx * tlx_z_zz_0[j] + fl1_fx * tlx_zz_z_0[j] + 0.5 * fl1_fx * tpy_zz_zz_0[j] + fl1_fx * fl1_fgb * tdy_zz_zz_0[j];

                tly_zzz_zz_0[j] = pa_z[j] * tly_zz_zz_0[j] + fl1_fx * tly_z_zz_0[j] + fl1_fx * tly_zz_z_0[j] - 0.5 * fl1_fx * tpx_zz_zz_0[j] - fl1_fx * fl1_fgb * tdx_zz_zz_0[j];

                tlz_zzz_zz_0[j] = pa_z[j] * tlz_zz_zz_0[j] + fl1_fx * tlz_z_zz_0[j] + fl1_fx * tlz_zz_z_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDG(      CMemBlock2D<double>& primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        amomrecfunc::compAngularMomentumForDG_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        amomrecfunc::compAngularMomentumForDG_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        amomrecfunc::compAngularMomentumForDG_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        amomrecfunc::compAngularMomentumForDG_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        amomrecfunc::compAngularMomentumForDG_180_225(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        amomrecfunc::compAngularMomentumForDG_225_270(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compAngularMomentumForDG_0_45(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_x_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx); 

            auto tly_x_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx); 

            auto tlz_x_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx); 

            auto tlx_x_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 1); 

            auto tly_x_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 1); 

            auto tlz_x_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 1); 

            auto tlx_x_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 2); 

            auto tly_x_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 2); 

            auto tlz_x_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 2); 

            auto tlx_x_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 3); 

            auto tly_x_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 3); 

            auto tlz_x_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 3); 

            auto tlx_x_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 4); 

            auto tly_x_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 4); 

            auto tlz_x_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 4); 

            auto tlx_x_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 5); 

            auto tly_x_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 5); 

            auto tlz_x_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 5); 

            auto tlx_x_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 6); 

            auto tly_x_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 6); 

            auto tlz_x_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 6); 

            auto tlx_x_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 7); 

            auto tly_x_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 7); 

            auto tlz_x_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 7); 

            auto tlx_x_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 8); 

            auto tly_x_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 8); 

            auto tlz_x_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 8); 

            auto tlx_x_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 9); 

            auto tly_x_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 9); 

            auto tlz_x_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 9); 

            auto tlx_x_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 10); 

            auto tly_x_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 10); 

            auto tlz_x_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 10); 

            auto tlx_x_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 11); 

            auto tly_x_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 11); 

            auto tlz_x_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 11); 

            auto tlx_x_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 12); 

            auto tly_x_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 12); 

            auto tlz_x_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 12); 

            auto tlx_x_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 13); 

            auto tly_x_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 13); 

            auto tlz_x_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 13); 

            auto tlx_x_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 14); 

            auto tly_x_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 14); 

            auto tlz_x_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 14); 

            auto tlx_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx); 

            auto tly_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx); 

            auto tlz_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx); 

            auto tlx_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 1); 

            auto tly_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 1); 

            auto tlz_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 1); 

            auto tlx_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 2); 

            auto tly_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 2); 

            auto tlz_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 2); 

            auto tlx_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 3); 

            auto tly_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 3); 

            auto tlz_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 3); 

            auto tlx_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 4); 

            auto tly_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 4); 

            auto tlz_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 4); 

            auto tlx_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 5); 

            auto tly_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 5); 

            auto tlz_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 5); 

            auto tlx_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 6); 

            auto tly_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 6); 

            auto tlz_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 6); 

            auto tlx_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 7); 

            auto tly_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 7); 

            auto tlz_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 7); 

            auto tlx_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 8); 

            auto tly_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 8); 

            auto tlz_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 8); 

            auto tlx_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 9); 

            auto tly_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 9); 

            auto tlz_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 9); 

            auto tlx_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 10); 

            auto tly_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 10); 

            auto tlz_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 10); 

            auto tlx_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 11); 

            auto tly_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 11); 

            auto tlz_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 11); 

            auto tlx_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 12); 

            auto tly_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 12); 

            auto tlz_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 12); 

            auto tlx_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 13); 

            auto tly_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 13); 

            auto tlz_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 13); 

            auto tlx_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 14); 

            auto tly_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 14); 

            auto tlz_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 14); 

            auto tlx_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx); 

            auto tly_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx); 

            auto tlz_x_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx); 

            auto tlx_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 1); 

            auto tly_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 1); 

            auto tlz_x_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 1); 

            auto tlx_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 2); 

            auto tly_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 2); 

            auto tlz_x_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 2); 

            auto tlx_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 3); 

            auto tly_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 3); 

            auto tlz_x_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 3); 

            auto tlx_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 4); 

            auto tly_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 4); 

            auto tlz_x_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 4); 

            auto tlx_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 5); 

            auto tly_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 5); 

            auto tlz_x_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 5); 

            auto tlx_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 6); 

            auto tly_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 6); 

            auto tlz_x_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 6); 

            auto tlx_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 7); 

            auto tly_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 7); 

            auto tlz_x_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 7); 

            auto tlx_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 8); 

            auto tly_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 8); 

            auto tlz_x_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 8); 

            auto tlx_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 9); 

            auto tly_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 9); 

            auto tlz_x_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 9); 

            auto tpy_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx); 

            auto tpz_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx); 

            auto tpy_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 1); 

            auto tpz_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 1); 

            auto tpy_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 2); 

            auto tpz_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 2); 

            auto tpy_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 3); 

            auto tpz_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 3); 

            auto tpy_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 4); 

            auto tpz_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 4); 

            auto tpy_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 5); 

            auto tpz_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 5); 

            auto tpy_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 6); 

            auto tpz_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 6); 

            auto tpy_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 7); 

            auto tpz_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 7); 

            auto tpy_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 8); 

            auto tpz_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 8); 

            auto tpy_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 9); 

            auto tpz_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 9); 

            auto tpy_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 10); 

            auto tpz_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 10); 

            auto tpy_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 11); 

            auto tpz_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 11); 

            auto tpy_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 12); 

            auto tpz_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 12); 

            auto tpy_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 13); 

            auto tpz_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 13); 

            auto tpy_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 14); 

            auto tpz_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 14); 

            auto tdy_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx); 

            auto tdz_x_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx); 

            auto tdy_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 1); 

            auto tdz_x_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 1); 

            auto tdy_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 2); 

            auto tdz_x_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 2); 

            auto tdy_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 3); 

            auto tdz_x_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 3); 

            auto tdy_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 4); 

            auto tdz_x_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 4); 

            auto tdy_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 5); 

            auto tdz_x_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 5); 

            auto tdy_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 6); 

            auto tdz_x_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 6); 

            auto tdy_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 7); 

            auto tdz_x_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 7); 

            auto tdy_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 8); 

            auto tdz_x_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 8); 

            auto tdy_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 9); 

            auto tdz_x_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 9); 

            auto tdy_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 10); 

            auto tdz_x_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 10); 

            auto tdy_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 11); 

            auto tdz_x_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 11); 

            auto tdy_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 12); 

            auto tdz_x_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 12); 

            auto tdy_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 13); 

            auto tdz_x_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 13); 

            auto tdy_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 14); 

            auto tdz_x_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 14); 

            // set up pointers to integrals

            auto tlx_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx); 

            auto tly_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx); 

            auto tlz_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx); 

            auto tlx_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 1); 

            auto tly_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 1); 

            auto tlz_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 1); 

            auto tlx_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 2); 

            auto tly_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 2); 

            auto tlz_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 2); 

            auto tlx_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 3); 

            auto tly_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 3); 

            auto tlz_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 3); 

            auto tlx_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 4); 

            auto tly_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 4); 

            auto tlz_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 4); 

            auto tlx_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 5); 

            auto tly_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 5); 

            auto tlz_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 5); 

            auto tlx_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 6); 

            auto tly_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 6); 

            auto tlz_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 6); 

            auto tlx_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 7); 

            auto tly_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 7); 

            auto tlz_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 7); 

            auto tlx_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 8); 

            auto tly_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 8); 

            auto tlz_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 8); 

            auto tlx_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 9); 

            auto tly_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 9); 

            auto tlz_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 9); 

            auto tlx_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 10); 

            auto tly_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 10); 

            auto tlz_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 10); 

            auto tlx_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 11); 

            auto tly_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 11); 

            auto tlz_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 11); 

            auto tlx_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 12); 

            auto tly_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 12); 

            auto tlz_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 12); 

            auto tlx_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 13); 

            auto tly_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 13); 

            auto tlz_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 13); 

            auto tlx_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 14); 

            auto tly_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 14); 

            auto tlz_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_x_xxxx_0, tdy_x_xxxy_0, tdy_x_xxxz_0, tdy_x_xxyy_0, \
                                     tdy_x_xxyz_0, tdy_x_xxzz_0, tdy_x_xyyy_0, tdy_x_xyyz_0, tdy_x_xyzz_0, tdy_x_xzzz_0, \
                                     tdy_x_yyyy_0, tdy_x_yyyz_0, tdy_x_yyzz_0, tdy_x_yzzz_0, tdy_x_zzzz_0, tdz_x_xxxx_0, \
                                     tdz_x_xxxy_0, tdz_x_xxxz_0, tdz_x_xxyy_0, tdz_x_xxyz_0, tdz_x_xxzz_0, tdz_x_xyyy_0, \
                                     tdz_x_xyyz_0, tdz_x_xyzz_0, tdz_x_xzzz_0, tdz_x_yyyy_0, tdz_x_yyyz_0, tdz_x_yyzz_0, \
                                     tdz_x_yzzz_0, tdz_x_zzzz_0, tlx_0_xxxx_0, tlx_0_xxxy_0, tlx_0_xxxz_0, tlx_0_xxyy_0, \
                                     tlx_0_xxyz_0, tlx_0_xxzz_0, tlx_0_xyyy_0, tlx_0_xyyz_0, tlx_0_xyzz_0, tlx_0_xzzz_0, \
                                     tlx_0_yyyy_0, tlx_0_yyyz_0, tlx_0_yyzz_0, tlx_0_yzzz_0, tlx_0_zzzz_0, tlx_x_xxx_0, \
                                     tlx_x_xxxx_0, tlx_x_xxxy_0, tlx_x_xxxz_0, tlx_x_xxy_0, tlx_x_xxyy_0, tlx_x_xxyz_0, \
                                     tlx_x_xxz_0, tlx_x_xxzz_0, tlx_x_xyy_0, tlx_x_xyyy_0, tlx_x_xyyz_0, tlx_x_xyz_0, \
                                     tlx_x_xyzz_0, tlx_x_xzz_0, tlx_x_xzzz_0, tlx_x_yyy_0, tlx_x_yyyy_0, tlx_x_yyyz_0, \
                                     tlx_x_yyz_0, tlx_x_yyzz_0, tlx_x_yzz_0, tlx_x_yzzz_0, tlx_x_zzz_0, tlx_x_zzzz_0, \
                                     tlx_xx_xxxx_0, tlx_xx_xxxy_0, tlx_xx_xxxz_0, tlx_xx_xxyy_0, tlx_xx_xxyz_0, \
                                     tlx_xx_xxzz_0, tlx_xx_xyyy_0, tlx_xx_xyyz_0, tlx_xx_xyzz_0, tlx_xx_xzzz_0, \
                                     tlx_xx_yyyy_0, tlx_xx_yyyz_0, tlx_xx_yyzz_0, tlx_xx_yzzz_0, tlx_xx_zzzz_0, \
                                     tly_0_xxxx_0, tly_0_xxxy_0, tly_0_xxxz_0, tly_0_xxyy_0, tly_0_xxyz_0, tly_0_xxzz_0, \
                                     tly_0_xyyy_0, tly_0_xyyz_0, tly_0_xyzz_0, tly_0_xzzz_0, tly_0_yyyy_0, tly_0_yyyz_0, \
                                     tly_0_yyzz_0, tly_0_yzzz_0, tly_0_zzzz_0, tly_x_xxx_0, tly_x_xxxx_0, tly_x_xxxy_0, \
                                     tly_x_xxxz_0, tly_x_xxy_0, tly_x_xxyy_0, tly_x_xxyz_0, tly_x_xxz_0, tly_x_xxzz_0, \
                                     tly_x_xyy_0, tly_x_xyyy_0, tly_x_xyyz_0, tly_x_xyz_0, tly_x_xyzz_0, tly_x_xzz_0, \
                                     tly_x_xzzz_0, tly_x_yyy_0, tly_x_yyyy_0, tly_x_yyyz_0, tly_x_yyz_0, tly_x_yyzz_0, \
                                     tly_x_yzz_0, tly_x_yzzz_0, tly_x_zzz_0, tly_x_zzzz_0, tly_xx_xxxx_0, \
                                     tly_xx_xxxy_0, tly_xx_xxxz_0, tly_xx_xxyy_0, tly_xx_xxyz_0, tly_xx_xxzz_0, \
                                     tly_xx_xyyy_0, tly_xx_xyyz_0, tly_xx_xyzz_0, tly_xx_xzzz_0, tly_xx_yyyy_0, \
                                     tly_xx_yyyz_0, tly_xx_yyzz_0, tly_xx_yzzz_0, tly_xx_zzzz_0, tlz_0_xxxx_0, \
                                     tlz_0_xxxy_0, tlz_0_xxxz_0, tlz_0_xxyy_0, tlz_0_xxyz_0, tlz_0_xxzz_0, tlz_0_xyyy_0, \
                                     tlz_0_xyyz_0, tlz_0_xyzz_0, tlz_0_xzzz_0, tlz_0_yyyy_0, tlz_0_yyyz_0, tlz_0_yyzz_0, \
                                     tlz_0_yzzz_0, tlz_0_zzzz_0, tlz_x_xxx_0, tlz_x_xxxx_0, tlz_x_xxxy_0, tlz_x_xxxz_0, \
                                     tlz_x_xxy_0, tlz_x_xxyy_0, tlz_x_xxyz_0, tlz_x_xxz_0, tlz_x_xxzz_0, tlz_x_xyy_0, \
                                     tlz_x_xyyy_0, tlz_x_xyyz_0, tlz_x_xyz_0, tlz_x_xyzz_0, tlz_x_xzz_0, tlz_x_xzzz_0, \
                                     tlz_x_yyy_0, tlz_x_yyyy_0, tlz_x_yyyz_0, tlz_x_yyz_0, tlz_x_yyzz_0, tlz_x_yzz_0, \
                                     tlz_x_yzzz_0, tlz_x_zzz_0, tlz_x_zzzz_0, tlz_xx_xxxx_0, tlz_xx_xxxy_0, \
                                     tlz_xx_xxxz_0, tlz_xx_xxyy_0, tlz_xx_xxyz_0, tlz_xx_xxzz_0, tlz_xx_xyyy_0, \
                                     tlz_xx_xyyz_0, tlz_xx_xyzz_0, tlz_xx_xzzz_0, tlz_xx_yyyy_0, tlz_xx_yyyz_0, \
                                     tlz_xx_yyzz_0, tlz_xx_yzzz_0, tlz_xx_zzzz_0, tpy_x_xxxx_0, tpy_x_xxxy_0, \
                                     tpy_x_xxxz_0, tpy_x_xxyy_0, tpy_x_xxyz_0, tpy_x_xxzz_0, tpy_x_xyyy_0, tpy_x_xyyz_0, \
                                     tpy_x_xyzz_0, tpy_x_xzzz_0, tpy_x_yyyy_0, tpy_x_yyyz_0, tpy_x_yyzz_0, tpy_x_yzzz_0, \
                                     tpy_x_zzzz_0, tpz_x_xxxx_0, tpz_x_xxxy_0, tpz_x_xxxz_0, tpz_x_xxyy_0, tpz_x_xxyz_0, \
                                     tpz_x_xxzz_0, tpz_x_xyyy_0, tpz_x_xyyz_0, tpz_x_xyzz_0, tpz_x_xzzz_0, tpz_x_yyyy_0, \
                                     tpz_x_yyyz_0, tpz_x_yyzz_0, tpz_x_yzzz_0, tpz_x_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xx_xxxx_0[j] = pa_x[j] * tlx_x_xxxx_0[j] + 0.5 * fl1_fx * tlx_0_xxxx_0[j] + 2.0 * fl1_fx * tlx_x_xxx_0[j];

                tly_xx_xxxx_0[j] = pa_x[j] * tly_x_xxxx_0[j] + 0.5 * fl1_fx * tly_0_xxxx_0[j] + 2.0 * fl1_fx * tly_x_xxx_0[j] + 0.5 * fl1_fx * tpz_x_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_x_xxxx_0[j];

                tlz_xx_xxxx_0[j] = pa_x[j] * tlz_x_xxxx_0[j] + 0.5 * fl1_fx * tlz_0_xxxx_0[j] + 2.0 * fl1_fx * tlz_x_xxx_0[j] - 0.5 * fl1_fx * tpy_x_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_x_xxxx_0[j];

                tlx_xx_xxxy_0[j] = pa_x[j] * tlx_x_xxxy_0[j] + 0.5 * fl1_fx * tlx_0_xxxy_0[j] + 1.5 * fl1_fx * tlx_x_xxy_0[j];

                tly_xx_xxxy_0[j] = pa_x[j] * tly_x_xxxy_0[j] + 0.5 * fl1_fx * tly_0_xxxy_0[j] + 1.5 * fl1_fx * tly_x_xxy_0[j] + 0.5 * fl1_fx * tpz_x_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_x_xxxy_0[j];

                tlz_xx_xxxy_0[j] = pa_x[j] * tlz_x_xxxy_0[j] + 0.5 * fl1_fx * tlz_0_xxxy_0[j] + 1.5 * fl1_fx * tlz_x_xxy_0[j] - 0.5 * fl1_fx * tpy_x_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_x_xxxy_0[j];

                tlx_xx_xxxz_0[j] = pa_x[j] * tlx_x_xxxz_0[j] + 0.5 * fl1_fx * tlx_0_xxxz_0[j] + 1.5 * fl1_fx * tlx_x_xxz_0[j];

                tly_xx_xxxz_0[j] = pa_x[j] * tly_x_xxxz_0[j] + 0.5 * fl1_fx * tly_0_xxxz_0[j] + 1.5 * fl1_fx * tly_x_xxz_0[j] + 0.5 * fl1_fx * tpz_x_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_x_xxxz_0[j];

                tlz_xx_xxxz_0[j] = pa_x[j] * tlz_x_xxxz_0[j] + 0.5 * fl1_fx * tlz_0_xxxz_0[j] + 1.5 * fl1_fx * tlz_x_xxz_0[j] - 0.5 * fl1_fx * tpy_x_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_x_xxxz_0[j];

                tlx_xx_xxyy_0[j] = pa_x[j] * tlx_x_xxyy_0[j] + 0.5 * fl1_fx * tlx_0_xxyy_0[j] + fl1_fx * tlx_x_xyy_0[j];

                tly_xx_xxyy_0[j] = pa_x[j] * tly_x_xxyy_0[j] + 0.5 * fl1_fx * tly_0_xxyy_0[j] + fl1_fx * tly_x_xyy_0[j] + 0.5 * fl1_fx * tpz_x_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_x_xxyy_0[j];

                tlz_xx_xxyy_0[j] = pa_x[j] * tlz_x_xxyy_0[j] + 0.5 * fl1_fx * tlz_0_xxyy_0[j] + fl1_fx * tlz_x_xyy_0[j] - 0.5 * fl1_fx * tpy_x_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_x_xxyy_0[j];

                tlx_xx_xxyz_0[j] = pa_x[j] * tlx_x_xxyz_0[j] + 0.5 * fl1_fx * tlx_0_xxyz_0[j] + fl1_fx * tlx_x_xyz_0[j];

                tly_xx_xxyz_0[j] = pa_x[j] * tly_x_xxyz_0[j] + 0.5 * fl1_fx * tly_0_xxyz_0[j] + fl1_fx * tly_x_xyz_0[j] + 0.5 * fl1_fx * tpz_x_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_x_xxyz_0[j];

                tlz_xx_xxyz_0[j] = pa_x[j] * tlz_x_xxyz_0[j] + 0.5 * fl1_fx * tlz_0_xxyz_0[j] + fl1_fx * tlz_x_xyz_0[j] - 0.5 * fl1_fx * tpy_x_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_x_xxyz_0[j];

                tlx_xx_xxzz_0[j] = pa_x[j] * tlx_x_xxzz_0[j] + 0.5 * fl1_fx * tlx_0_xxzz_0[j] + fl1_fx * tlx_x_xzz_0[j];

                tly_xx_xxzz_0[j] = pa_x[j] * tly_x_xxzz_0[j] + 0.5 * fl1_fx * tly_0_xxzz_0[j] + fl1_fx * tly_x_xzz_0[j] + 0.5 * fl1_fx * tpz_x_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_x_xxzz_0[j];

                tlz_xx_xxzz_0[j] = pa_x[j] * tlz_x_xxzz_0[j] + 0.5 * fl1_fx * tlz_0_xxzz_0[j] + fl1_fx * tlz_x_xzz_0[j] - 0.5 * fl1_fx * tpy_x_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_x_xxzz_0[j];

                tlx_xx_xyyy_0[j] = pa_x[j] * tlx_x_xyyy_0[j] + 0.5 * fl1_fx * tlx_0_xyyy_0[j] + 0.5 * fl1_fx * tlx_x_yyy_0[j];

                tly_xx_xyyy_0[j] = pa_x[j] * tly_x_xyyy_0[j] + 0.5 * fl1_fx * tly_0_xyyy_0[j] + 0.5 * fl1_fx * tly_x_yyy_0[j] + 0.5 * fl1_fx * tpz_x_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_x_xyyy_0[j];

                tlz_xx_xyyy_0[j] = pa_x[j] * tlz_x_xyyy_0[j] + 0.5 * fl1_fx * tlz_0_xyyy_0[j] + 0.5 * fl1_fx * tlz_x_yyy_0[j] - 0.5 * fl1_fx * tpy_x_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_x_xyyy_0[j];

                tlx_xx_xyyz_0[j] = pa_x[j] * tlx_x_xyyz_0[j] + 0.5 * fl1_fx * tlx_0_xyyz_0[j] + 0.5 * fl1_fx * tlx_x_yyz_0[j];

                tly_xx_xyyz_0[j] = pa_x[j] * tly_x_xyyz_0[j] + 0.5 * fl1_fx * tly_0_xyyz_0[j] + 0.5 * fl1_fx * tly_x_yyz_0[j] + 0.5 * fl1_fx * tpz_x_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_x_xyyz_0[j];

                tlz_xx_xyyz_0[j] = pa_x[j] * tlz_x_xyyz_0[j] + 0.5 * fl1_fx * tlz_0_xyyz_0[j] + 0.5 * fl1_fx * tlz_x_yyz_0[j] - 0.5 * fl1_fx * tpy_x_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_x_xyyz_0[j];

                tlx_xx_xyzz_0[j] = pa_x[j] * tlx_x_xyzz_0[j] + 0.5 * fl1_fx * tlx_0_xyzz_0[j] + 0.5 * fl1_fx * tlx_x_yzz_0[j];

                tly_xx_xyzz_0[j] = pa_x[j] * tly_x_xyzz_0[j] + 0.5 * fl1_fx * tly_0_xyzz_0[j] + 0.5 * fl1_fx * tly_x_yzz_0[j] + 0.5 * fl1_fx * tpz_x_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_x_xyzz_0[j];

                tlz_xx_xyzz_0[j] = pa_x[j] * tlz_x_xyzz_0[j] + 0.5 * fl1_fx * tlz_0_xyzz_0[j] + 0.5 * fl1_fx * tlz_x_yzz_0[j] - 0.5 * fl1_fx * tpy_x_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_x_xyzz_0[j];

                tlx_xx_xzzz_0[j] = pa_x[j] * tlx_x_xzzz_0[j] + 0.5 * fl1_fx * tlx_0_xzzz_0[j] + 0.5 * fl1_fx * tlx_x_zzz_0[j];

                tly_xx_xzzz_0[j] = pa_x[j] * tly_x_xzzz_0[j] + 0.5 * fl1_fx * tly_0_xzzz_0[j] + 0.5 * fl1_fx * tly_x_zzz_0[j] + 0.5 * fl1_fx * tpz_x_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_x_xzzz_0[j];

                tlz_xx_xzzz_0[j] = pa_x[j] * tlz_x_xzzz_0[j] + 0.5 * fl1_fx * tlz_0_xzzz_0[j] + 0.5 * fl1_fx * tlz_x_zzz_0[j] - 0.5 * fl1_fx * tpy_x_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_x_xzzz_0[j];

                tlx_xx_yyyy_0[j] = pa_x[j] * tlx_x_yyyy_0[j] + 0.5 * fl1_fx * tlx_0_yyyy_0[j];

                tly_xx_yyyy_0[j] = pa_x[j] * tly_x_yyyy_0[j] + 0.5 * fl1_fx * tly_0_yyyy_0[j] + 0.5 * fl1_fx * tpz_x_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_x_yyyy_0[j];

                tlz_xx_yyyy_0[j] = pa_x[j] * tlz_x_yyyy_0[j] + 0.5 * fl1_fx * tlz_0_yyyy_0[j] - 0.5 * fl1_fx * tpy_x_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_x_yyyy_0[j];

                tlx_xx_yyyz_0[j] = pa_x[j] * tlx_x_yyyz_0[j] + 0.5 * fl1_fx * tlx_0_yyyz_0[j];

                tly_xx_yyyz_0[j] = pa_x[j] * tly_x_yyyz_0[j] + 0.5 * fl1_fx * tly_0_yyyz_0[j] + 0.5 * fl1_fx * tpz_x_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_x_yyyz_0[j];

                tlz_xx_yyyz_0[j] = pa_x[j] * tlz_x_yyyz_0[j] + 0.5 * fl1_fx * tlz_0_yyyz_0[j] - 0.5 * fl1_fx * tpy_x_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_x_yyyz_0[j];

                tlx_xx_yyzz_0[j] = pa_x[j] * tlx_x_yyzz_0[j] + 0.5 * fl1_fx * tlx_0_yyzz_0[j];

                tly_xx_yyzz_0[j] = pa_x[j] * tly_x_yyzz_0[j] + 0.5 * fl1_fx * tly_0_yyzz_0[j] + 0.5 * fl1_fx * tpz_x_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_x_yyzz_0[j];

                tlz_xx_yyzz_0[j] = pa_x[j] * tlz_x_yyzz_0[j] + 0.5 * fl1_fx * tlz_0_yyzz_0[j] - 0.5 * fl1_fx * tpy_x_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_x_yyzz_0[j];

                tlx_xx_yzzz_0[j] = pa_x[j] * tlx_x_yzzz_0[j] + 0.5 * fl1_fx * tlx_0_yzzz_0[j];

                tly_xx_yzzz_0[j] = pa_x[j] * tly_x_yzzz_0[j] + 0.5 * fl1_fx * tly_0_yzzz_0[j] + 0.5 * fl1_fx * tpz_x_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_x_yzzz_0[j];

                tlz_xx_yzzz_0[j] = pa_x[j] * tlz_x_yzzz_0[j] + 0.5 * fl1_fx * tlz_0_yzzz_0[j] - 0.5 * fl1_fx * tpy_x_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_x_yzzz_0[j];

                tlx_xx_zzzz_0[j] = pa_x[j] * tlx_x_zzzz_0[j] + 0.5 * fl1_fx * tlx_0_zzzz_0[j];

                tly_xx_zzzz_0[j] = pa_x[j] * tly_x_zzzz_0[j] + 0.5 * fl1_fx * tly_0_zzzz_0[j] + 0.5 * fl1_fx * tpz_x_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_x_zzzz_0[j];

                tlz_xx_zzzz_0[j] = pa_x[j] * tlz_x_zzzz_0[j] + 0.5 * fl1_fx * tlz_0_zzzz_0[j] - 0.5 * fl1_fx * tpy_x_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_x_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDG_45_90(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 15); 

            auto tly_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 15); 

            auto tlz_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tlx_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 16); 

            auto tly_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 16); 

            auto tlz_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tlx_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 17); 

            auto tly_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 17); 

            auto tlz_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tlx_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 18); 

            auto tly_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 18); 

            auto tlz_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tlx_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 19); 

            auto tly_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 19); 

            auto tlz_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tlx_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 20); 

            auto tly_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 20); 

            auto tlz_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tlx_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 21); 

            auto tly_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 21); 

            auto tlz_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tlx_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 22); 

            auto tly_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 22); 

            auto tlz_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tlx_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 23); 

            auto tly_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 23); 

            auto tlz_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tlx_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 24); 

            auto tly_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 24); 

            auto tlz_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tlx_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 25); 

            auto tly_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 25); 

            auto tlz_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tlx_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 26); 

            auto tly_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 26); 

            auto tlz_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tlx_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 27); 

            auto tly_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 27); 

            auto tlz_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tlx_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 28); 

            auto tly_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 28); 

            auto tlz_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tlx_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 29); 

            auto tly_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 29); 

            auto tlz_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            auto tlx_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 10); 

            auto tly_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tlz_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tlx_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 11); 

            auto tly_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tlz_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tlx_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 12); 

            auto tly_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tlz_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tlx_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 13); 

            auto tly_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tlz_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tlx_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 14); 

            auto tly_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tlz_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tlx_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 15); 

            auto tly_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tlz_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tlx_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 16); 

            auto tly_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tlz_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tlx_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 17); 

            auto tly_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tlz_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tlx_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 18); 

            auto tly_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tlz_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tlx_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 19); 

            auto tly_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tlz_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tpy_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 15); 

            auto tpz_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tpy_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 16); 

            auto tpz_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tpy_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 17); 

            auto tpz_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tpy_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 18); 

            auto tpz_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tpy_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 19); 

            auto tpz_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tpy_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 20); 

            auto tpz_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tpy_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 21); 

            auto tpz_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tpy_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 22); 

            auto tpz_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tpy_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 23); 

            auto tpz_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tpy_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 24); 

            auto tpz_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tpy_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 25); 

            auto tpz_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tpy_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 26); 

            auto tpz_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tpy_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 27); 

            auto tpz_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tpy_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 28); 

            auto tpz_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tpy_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 29); 

            auto tpz_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            auto tdy_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 15); 

            auto tdz_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tdy_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 16); 

            auto tdz_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tdy_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 17); 

            auto tdz_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tdy_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 18); 

            auto tdz_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tdy_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 19); 

            auto tdz_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tdy_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 20); 

            auto tdz_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tdy_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 21); 

            auto tdz_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tdy_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 22); 

            auto tdz_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tdy_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 23); 

            auto tdz_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tdy_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 24); 

            auto tdz_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tdy_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 25); 

            auto tdz_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tdy_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 26); 

            auto tdz_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tdy_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 27); 

            auto tdz_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tdy_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 28); 

            auto tdz_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tdy_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 29); 

            auto tdz_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            // set up pointers to integrals

            auto tlx_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 15); 

            auto tly_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 15); 

            auto tlz_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 15); 

            auto tlx_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 16); 

            auto tly_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 16); 

            auto tlz_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 16); 

            auto tlx_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 17); 

            auto tly_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 17); 

            auto tlz_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 17); 

            auto tlx_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 18); 

            auto tly_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 18); 

            auto tlz_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 18); 

            auto tlx_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 19); 

            auto tly_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 19); 

            auto tlz_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 19); 

            auto tlx_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 20); 

            auto tly_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 20); 

            auto tlz_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 20); 

            auto tlx_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 21); 

            auto tly_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 21); 

            auto tlz_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 21); 

            auto tlx_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 22); 

            auto tly_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 22); 

            auto tlz_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 22); 

            auto tlx_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 23); 

            auto tly_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 23); 

            auto tlz_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 23); 

            auto tlx_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 24); 

            auto tly_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 24); 

            auto tlz_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 24); 

            auto tlx_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 25); 

            auto tly_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 25); 

            auto tlz_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 25); 

            auto tlx_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 26); 

            auto tly_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 26); 

            auto tlz_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 26); 

            auto tlx_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 27); 

            auto tly_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 27); 

            auto tlz_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 27); 

            auto tlx_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 28); 

            auto tly_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 28); 

            auto tlz_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 28); 

            auto tlx_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 29); 

            auto tly_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 29); 

            auto tlz_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_y_xxxx_0, tdy_y_xxxy_0, tdy_y_xxxz_0, tdy_y_xxyy_0, \
                                     tdy_y_xxyz_0, tdy_y_xxzz_0, tdy_y_xyyy_0, tdy_y_xyyz_0, tdy_y_xyzz_0, tdy_y_xzzz_0, \
                                     tdy_y_yyyy_0, tdy_y_yyyz_0, tdy_y_yyzz_0, tdy_y_yzzz_0, tdy_y_zzzz_0, tdz_y_xxxx_0, \
                                     tdz_y_xxxy_0, tdz_y_xxxz_0, tdz_y_xxyy_0, tdz_y_xxyz_0, tdz_y_xxzz_0, tdz_y_xyyy_0, \
                                     tdz_y_xyyz_0, tdz_y_xyzz_0, tdz_y_xzzz_0, tdz_y_yyyy_0, tdz_y_yyyz_0, tdz_y_yyzz_0, \
                                     tdz_y_yzzz_0, tdz_y_zzzz_0, tlx_xy_xxxx_0, tlx_xy_xxxy_0, tlx_xy_xxxz_0, \
                                     tlx_xy_xxyy_0, tlx_xy_xxyz_0, tlx_xy_xxzz_0, tlx_xy_xyyy_0, tlx_xy_xyyz_0, \
                                     tlx_xy_xyzz_0, tlx_xy_xzzz_0, tlx_xy_yyyy_0, tlx_xy_yyyz_0, tlx_xy_yyzz_0, \
                                     tlx_xy_yzzz_0, tlx_xy_zzzz_0, tlx_y_xxx_0, tlx_y_xxxx_0, tlx_y_xxxy_0, tlx_y_xxxz_0, \
                                     tlx_y_xxy_0, tlx_y_xxyy_0, tlx_y_xxyz_0, tlx_y_xxz_0, tlx_y_xxzz_0, tlx_y_xyy_0, \
                                     tlx_y_xyyy_0, tlx_y_xyyz_0, tlx_y_xyz_0, tlx_y_xyzz_0, tlx_y_xzz_0, tlx_y_xzzz_0, \
                                     tlx_y_yyy_0, tlx_y_yyyy_0, tlx_y_yyyz_0, tlx_y_yyz_0, tlx_y_yyzz_0, tlx_y_yzz_0, \
                                     tlx_y_yzzz_0, tlx_y_zzz_0, tlx_y_zzzz_0, tly_xy_xxxx_0, tly_xy_xxxy_0, \
                                     tly_xy_xxxz_0, tly_xy_xxyy_0, tly_xy_xxyz_0, tly_xy_xxzz_0, tly_xy_xyyy_0, \
                                     tly_xy_xyyz_0, tly_xy_xyzz_0, tly_xy_xzzz_0, tly_xy_yyyy_0, tly_xy_yyyz_0, \
                                     tly_xy_yyzz_0, tly_xy_yzzz_0, tly_xy_zzzz_0, tly_y_xxx_0, tly_y_xxxx_0, \
                                     tly_y_xxxy_0, tly_y_xxxz_0, tly_y_xxy_0, tly_y_xxyy_0, tly_y_xxyz_0, tly_y_xxz_0, \
                                     tly_y_xxzz_0, tly_y_xyy_0, tly_y_xyyy_0, tly_y_xyyz_0, tly_y_xyz_0, tly_y_xyzz_0, \
                                     tly_y_xzz_0, tly_y_xzzz_0, tly_y_yyy_0, tly_y_yyyy_0, tly_y_yyyz_0, tly_y_yyz_0, \
                                     tly_y_yyzz_0, tly_y_yzz_0, tly_y_yzzz_0, tly_y_zzz_0, tly_y_zzzz_0, tlz_xy_xxxx_0, \
                                     tlz_xy_xxxy_0, tlz_xy_xxxz_0, tlz_xy_xxyy_0, tlz_xy_xxyz_0, tlz_xy_xxzz_0, \
                                     tlz_xy_xyyy_0, tlz_xy_xyyz_0, tlz_xy_xyzz_0, tlz_xy_xzzz_0, tlz_xy_yyyy_0, \
                                     tlz_xy_yyyz_0, tlz_xy_yyzz_0, tlz_xy_yzzz_0, tlz_xy_zzzz_0, tlz_y_xxx_0, \
                                     tlz_y_xxxx_0, tlz_y_xxxy_0, tlz_y_xxxz_0, tlz_y_xxy_0, tlz_y_xxyy_0, tlz_y_xxyz_0, \
                                     tlz_y_xxz_0, tlz_y_xxzz_0, tlz_y_xyy_0, tlz_y_xyyy_0, tlz_y_xyyz_0, tlz_y_xyz_0, \
                                     tlz_y_xyzz_0, tlz_y_xzz_0, tlz_y_xzzz_0, tlz_y_yyy_0, tlz_y_yyyy_0, tlz_y_yyyz_0, \
                                     tlz_y_yyz_0, tlz_y_yyzz_0, tlz_y_yzz_0, tlz_y_yzzz_0, tlz_y_zzz_0, tlz_y_zzzz_0, \
                                     tpy_y_xxxx_0, tpy_y_xxxy_0, tpy_y_xxxz_0, tpy_y_xxyy_0, tpy_y_xxyz_0, tpy_y_xxzz_0, \
                                     tpy_y_xyyy_0, tpy_y_xyyz_0, tpy_y_xyzz_0, tpy_y_xzzz_0, tpy_y_yyyy_0, tpy_y_yyyz_0, \
                                     tpy_y_yyzz_0, tpy_y_yzzz_0, tpy_y_zzzz_0, tpz_y_xxxx_0, tpz_y_xxxy_0, tpz_y_xxxz_0, \
                                     tpz_y_xxyy_0, tpz_y_xxyz_0, tpz_y_xxzz_0, tpz_y_xyyy_0, tpz_y_xyyz_0, tpz_y_xyzz_0, \
                                     tpz_y_xzzz_0, tpz_y_yyyy_0, tpz_y_yyyz_0, tpz_y_yyzz_0, tpz_y_yzzz_0, tpz_y_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xy_xxxx_0[j] = pa_x[j] * tlx_y_xxxx_0[j] + 2.0 * fl1_fx * tlx_y_xxx_0[j];

                tly_xy_xxxx_0[j] = pa_x[j] * tly_y_xxxx_0[j] + 2.0 * fl1_fx * tly_y_xxx_0[j] + 0.5 * fl1_fx * tpz_y_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_y_xxxx_0[j];

                tlz_xy_xxxx_0[j] = pa_x[j] * tlz_y_xxxx_0[j] + 2.0 * fl1_fx * tlz_y_xxx_0[j] - 0.5 * fl1_fx * tpy_y_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_y_xxxx_0[j];

                tlx_xy_xxxy_0[j] = pa_x[j] * tlx_y_xxxy_0[j] + 1.5 * fl1_fx * tlx_y_xxy_0[j];

                tly_xy_xxxy_0[j] = pa_x[j] * tly_y_xxxy_0[j] + 1.5 * fl1_fx * tly_y_xxy_0[j] + 0.5 * fl1_fx * tpz_y_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_y_xxxy_0[j];

                tlz_xy_xxxy_0[j] = pa_x[j] * tlz_y_xxxy_0[j] + 1.5 * fl1_fx * tlz_y_xxy_0[j] - 0.5 * fl1_fx * tpy_y_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_y_xxxy_0[j];

                tlx_xy_xxxz_0[j] = pa_x[j] * tlx_y_xxxz_0[j] + 1.5 * fl1_fx * tlx_y_xxz_0[j];

                tly_xy_xxxz_0[j] = pa_x[j] * tly_y_xxxz_0[j] + 1.5 * fl1_fx * tly_y_xxz_0[j] + 0.5 * fl1_fx * tpz_y_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_y_xxxz_0[j];

                tlz_xy_xxxz_0[j] = pa_x[j] * tlz_y_xxxz_0[j] + 1.5 * fl1_fx * tlz_y_xxz_0[j] - 0.5 * fl1_fx * tpy_y_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_y_xxxz_0[j];

                tlx_xy_xxyy_0[j] = pa_x[j] * tlx_y_xxyy_0[j] + fl1_fx * tlx_y_xyy_0[j];

                tly_xy_xxyy_0[j] = pa_x[j] * tly_y_xxyy_0[j] + fl1_fx * tly_y_xyy_0[j] + 0.5 * fl1_fx * tpz_y_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_y_xxyy_0[j];

                tlz_xy_xxyy_0[j] = pa_x[j] * tlz_y_xxyy_0[j] + fl1_fx * tlz_y_xyy_0[j] - 0.5 * fl1_fx * tpy_y_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_y_xxyy_0[j];

                tlx_xy_xxyz_0[j] = pa_x[j] * tlx_y_xxyz_0[j] + fl1_fx * tlx_y_xyz_0[j];

                tly_xy_xxyz_0[j] = pa_x[j] * tly_y_xxyz_0[j] + fl1_fx * tly_y_xyz_0[j] + 0.5 * fl1_fx * tpz_y_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_y_xxyz_0[j];

                tlz_xy_xxyz_0[j] = pa_x[j] * tlz_y_xxyz_0[j] + fl1_fx * tlz_y_xyz_0[j] - 0.5 * fl1_fx * tpy_y_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_y_xxyz_0[j];

                tlx_xy_xxzz_0[j] = pa_x[j] * tlx_y_xxzz_0[j] + fl1_fx * tlx_y_xzz_0[j];

                tly_xy_xxzz_0[j] = pa_x[j] * tly_y_xxzz_0[j] + fl1_fx * tly_y_xzz_0[j] + 0.5 * fl1_fx * tpz_y_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_y_xxzz_0[j];

                tlz_xy_xxzz_0[j] = pa_x[j] * tlz_y_xxzz_0[j] + fl1_fx * tlz_y_xzz_0[j] - 0.5 * fl1_fx * tpy_y_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_y_xxzz_0[j];

                tlx_xy_xyyy_0[j] = pa_x[j] * tlx_y_xyyy_0[j] + 0.5 * fl1_fx * tlx_y_yyy_0[j];

                tly_xy_xyyy_0[j] = pa_x[j] * tly_y_xyyy_0[j] + 0.5 * fl1_fx * tly_y_yyy_0[j] + 0.5 * fl1_fx * tpz_y_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_y_xyyy_0[j];

                tlz_xy_xyyy_0[j] = pa_x[j] * tlz_y_xyyy_0[j] + 0.5 * fl1_fx * tlz_y_yyy_0[j] - 0.5 * fl1_fx * tpy_y_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_y_xyyy_0[j];

                tlx_xy_xyyz_0[j] = pa_x[j] * tlx_y_xyyz_0[j] + 0.5 * fl1_fx * tlx_y_yyz_0[j];

                tly_xy_xyyz_0[j] = pa_x[j] * tly_y_xyyz_0[j] + 0.5 * fl1_fx * tly_y_yyz_0[j] + 0.5 * fl1_fx * tpz_y_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_y_xyyz_0[j];

                tlz_xy_xyyz_0[j] = pa_x[j] * tlz_y_xyyz_0[j] + 0.5 * fl1_fx * tlz_y_yyz_0[j] - 0.5 * fl1_fx * tpy_y_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_y_xyyz_0[j];

                tlx_xy_xyzz_0[j] = pa_x[j] * tlx_y_xyzz_0[j] + 0.5 * fl1_fx * tlx_y_yzz_0[j];

                tly_xy_xyzz_0[j] = pa_x[j] * tly_y_xyzz_0[j] + 0.5 * fl1_fx * tly_y_yzz_0[j] + 0.5 * fl1_fx * tpz_y_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_y_xyzz_0[j];

                tlz_xy_xyzz_0[j] = pa_x[j] * tlz_y_xyzz_0[j] + 0.5 * fl1_fx * tlz_y_yzz_0[j] - 0.5 * fl1_fx * tpy_y_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_y_xyzz_0[j];

                tlx_xy_xzzz_0[j] = pa_x[j] * tlx_y_xzzz_0[j] + 0.5 * fl1_fx * tlx_y_zzz_0[j];

                tly_xy_xzzz_0[j] = pa_x[j] * tly_y_xzzz_0[j] + 0.5 * fl1_fx * tly_y_zzz_0[j] + 0.5 * fl1_fx * tpz_y_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_y_xzzz_0[j];

                tlz_xy_xzzz_0[j] = pa_x[j] * tlz_y_xzzz_0[j] + 0.5 * fl1_fx * tlz_y_zzz_0[j] - 0.5 * fl1_fx * tpy_y_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_y_xzzz_0[j];

                tlx_xy_yyyy_0[j] = pa_x[j] * tlx_y_yyyy_0[j];

                tly_xy_yyyy_0[j] = pa_x[j] * tly_y_yyyy_0[j] + 0.5 * fl1_fx * tpz_y_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_y_yyyy_0[j];

                tlz_xy_yyyy_0[j] = pa_x[j] * tlz_y_yyyy_0[j] - 0.5 * fl1_fx * tpy_y_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_y_yyyy_0[j];

                tlx_xy_yyyz_0[j] = pa_x[j] * tlx_y_yyyz_0[j];

                tly_xy_yyyz_0[j] = pa_x[j] * tly_y_yyyz_0[j] + 0.5 * fl1_fx * tpz_y_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_y_yyyz_0[j];

                tlz_xy_yyyz_0[j] = pa_x[j] * tlz_y_yyyz_0[j] - 0.5 * fl1_fx * tpy_y_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_y_yyyz_0[j];

                tlx_xy_yyzz_0[j] = pa_x[j] * tlx_y_yyzz_0[j];

                tly_xy_yyzz_0[j] = pa_x[j] * tly_y_yyzz_0[j] + 0.5 * fl1_fx * tpz_y_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_y_yyzz_0[j];

                tlz_xy_yyzz_0[j] = pa_x[j] * tlz_y_yyzz_0[j] - 0.5 * fl1_fx * tpy_y_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_y_yyzz_0[j];

                tlx_xy_yzzz_0[j] = pa_x[j] * tlx_y_yzzz_0[j];

                tly_xy_yzzz_0[j] = pa_x[j] * tly_y_yzzz_0[j] + 0.5 * fl1_fx * tpz_y_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_y_yzzz_0[j];

                tlz_xy_yzzz_0[j] = pa_x[j] * tlz_y_yzzz_0[j] - 0.5 * fl1_fx * tpy_y_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_y_yzzz_0[j];

                tlx_xy_zzzz_0[j] = pa_x[j] * tlx_y_zzzz_0[j];

                tly_xy_zzzz_0[j] = pa_x[j] * tly_y_zzzz_0[j] + 0.5 * fl1_fx * tpz_y_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_y_zzzz_0[j];

                tlz_xy_zzzz_0[j] = pa_x[j] * tlz_y_zzzz_0[j] - 0.5 * fl1_fx * tpy_y_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_y_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDG_90_135(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 30); 

            auto tly_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tlz_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tlx_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 31); 

            auto tly_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tlz_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tlx_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 32); 

            auto tly_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tlz_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tlx_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 33); 

            auto tly_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tlz_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tlx_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 34); 

            auto tly_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tlz_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tlx_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 35); 

            auto tly_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tlz_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tlx_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 36); 

            auto tly_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tlz_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tlx_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 37); 

            auto tly_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tlz_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tlx_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 38); 

            auto tly_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tlz_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tlx_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 39); 

            auto tly_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tlz_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tlx_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 40); 

            auto tly_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tlz_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tlx_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 41); 

            auto tly_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tlz_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tlx_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 42); 

            auto tly_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tlz_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tlx_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 43); 

            auto tly_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tlz_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tlx_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 44); 

            auto tly_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tlz_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tlx_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 20); 

            auto tly_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 21); 

            auto tly_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 22); 

            auto tly_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 23); 

            auto tly_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tlx_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 24); 

            auto tly_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tlx_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 25); 

            auto tly_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tlz_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tlx_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 26); 

            auto tly_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tlz_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tlx_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 27); 

            auto tly_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tlz_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tlx_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 28); 

            auto tly_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tlz_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tlx_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 29); 

            auto tly_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tlz_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            // set up pointers to integrals

            auto tlx_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 30); 

            auto tly_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 30); 

            auto tlz_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 30); 

            auto tlx_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 31); 

            auto tly_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 31); 

            auto tlz_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 31); 

            auto tlx_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 32); 

            auto tly_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 32); 

            auto tlz_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 32); 

            auto tlx_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 33); 

            auto tly_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 33); 

            auto tlz_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 33); 

            auto tlx_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 34); 

            auto tly_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 34); 

            auto tlz_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 34); 

            auto tlx_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 35); 

            auto tly_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 35); 

            auto tlz_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 35); 

            auto tlx_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 36); 

            auto tly_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 36); 

            auto tlz_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 36); 

            auto tlx_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 37); 

            auto tly_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 37); 

            auto tlz_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 37); 

            auto tlx_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 38); 

            auto tly_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 38); 

            auto tlz_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 38); 

            auto tlx_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 39); 

            auto tly_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 39); 

            auto tlz_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 39); 

            auto tlx_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 40); 

            auto tly_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 40); 

            auto tlz_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 40); 

            auto tlx_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 41); 

            auto tly_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 41); 

            auto tlz_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 41); 

            auto tlx_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 42); 

            auto tly_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 42); 

            auto tlz_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 42); 

            auto tlx_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 43); 

            auto tly_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 43); 

            auto tlz_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 43); 

            auto tlx_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 44); 

            auto tly_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 44); 

            auto tlz_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_z_xxxx_0, tdy_z_xxxy_0, tdy_z_xxxz_0, tdy_z_xxyy_0, \
                                     tdy_z_xxyz_0, tdy_z_xxzz_0, tdy_z_xyyy_0, tdy_z_xyyz_0, tdy_z_xyzz_0, tdy_z_xzzz_0, \
                                     tdy_z_yyyy_0, tdy_z_yyyz_0, tdy_z_yyzz_0, tdy_z_yzzz_0, tdy_z_zzzz_0, tdz_z_xxxx_0, \
                                     tdz_z_xxxy_0, tdz_z_xxxz_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxzz_0, tdz_z_xyyy_0, \
                                     tdz_z_xyyz_0, tdz_z_xyzz_0, tdz_z_xzzz_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyzz_0, \
                                     tdz_z_yzzz_0, tdz_z_zzzz_0, tlx_xz_xxxx_0, tlx_xz_xxxy_0, tlx_xz_xxxz_0, \
                                     tlx_xz_xxyy_0, tlx_xz_xxyz_0, tlx_xz_xxzz_0, tlx_xz_xyyy_0, tlx_xz_xyyz_0, \
                                     tlx_xz_xyzz_0, tlx_xz_xzzz_0, tlx_xz_yyyy_0, tlx_xz_yyyz_0, tlx_xz_yyzz_0, \
                                     tlx_xz_yzzz_0, tlx_xz_zzzz_0, tlx_z_xxx_0, tlx_z_xxxx_0, tlx_z_xxxy_0, tlx_z_xxxz_0, \
                                     tlx_z_xxy_0, tlx_z_xxyy_0, tlx_z_xxyz_0, tlx_z_xxz_0, tlx_z_xxzz_0, tlx_z_xyy_0, \
                                     tlx_z_xyyy_0, tlx_z_xyyz_0, tlx_z_xyz_0, tlx_z_xyzz_0, tlx_z_xzz_0, tlx_z_xzzz_0, \
                                     tlx_z_yyy_0, tlx_z_yyyy_0, tlx_z_yyyz_0, tlx_z_yyz_0, tlx_z_yyzz_0, tlx_z_yzz_0, \
                                     tlx_z_yzzz_0, tlx_z_zzz_0, tlx_z_zzzz_0, tly_xz_xxxx_0, tly_xz_xxxy_0, \
                                     tly_xz_xxxz_0, tly_xz_xxyy_0, tly_xz_xxyz_0, tly_xz_xxzz_0, tly_xz_xyyy_0, \
                                     tly_xz_xyyz_0, tly_xz_xyzz_0, tly_xz_xzzz_0, tly_xz_yyyy_0, tly_xz_yyyz_0, \
                                     tly_xz_yyzz_0, tly_xz_yzzz_0, tly_xz_zzzz_0, tly_z_xxx_0, tly_z_xxxx_0, \
                                     tly_z_xxxy_0, tly_z_xxxz_0, tly_z_xxy_0, tly_z_xxyy_0, tly_z_xxyz_0, tly_z_xxz_0, \
                                     tly_z_xxzz_0, tly_z_xyy_0, tly_z_xyyy_0, tly_z_xyyz_0, tly_z_xyz_0, tly_z_xyzz_0, \
                                     tly_z_xzz_0, tly_z_xzzz_0, tly_z_yyy_0, tly_z_yyyy_0, tly_z_yyyz_0, tly_z_yyz_0, \
                                     tly_z_yyzz_0, tly_z_yzz_0, tly_z_yzzz_0, tly_z_zzz_0, tly_z_zzzz_0, tlz_xz_xxxx_0, \
                                     tlz_xz_xxxy_0, tlz_xz_xxxz_0, tlz_xz_xxyy_0, tlz_xz_xxyz_0, tlz_xz_xxzz_0, \
                                     tlz_xz_xyyy_0, tlz_xz_xyyz_0, tlz_xz_xyzz_0, tlz_xz_xzzz_0, tlz_xz_yyyy_0, \
                                     tlz_xz_yyyz_0, tlz_xz_yyzz_0, tlz_xz_yzzz_0, tlz_xz_zzzz_0, tlz_z_xxx_0, \
                                     tlz_z_xxxx_0, tlz_z_xxxy_0, tlz_z_xxxz_0, tlz_z_xxy_0, tlz_z_xxyy_0, tlz_z_xxyz_0, \
                                     tlz_z_xxz_0, tlz_z_xxzz_0, tlz_z_xyy_0, tlz_z_xyyy_0, tlz_z_xyyz_0, tlz_z_xyz_0, \
                                     tlz_z_xyzz_0, tlz_z_xzz_0, tlz_z_xzzz_0, tlz_z_yyy_0, tlz_z_yyyy_0, tlz_z_yyyz_0, \
                                     tlz_z_yyz_0, tlz_z_yyzz_0, tlz_z_yzz_0, tlz_z_yzzz_0, tlz_z_zzz_0, tlz_z_zzzz_0, \
                                     tpy_z_xxxx_0, tpy_z_xxxy_0, tpy_z_xxxz_0, tpy_z_xxyy_0, tpy_z_xxyz_0, tpy_z_xxzz_0, \
                                     tpy_z_xyyy_0, tpy_z_xyyz_0, tpy_z_xyzz_0, tpy_z_xzzz_0, tpy_z_yyyy_0, tpy_z_yyyz_0, \
                                     tpy_z_yyzz_0, tpy_z_yzzz_0, tpy_z_zzzz_0, tpz_z_xxxx_0, tpz_z_xxxy_0, tpz_z_xxxz_0, \
                                     tpz_z_xxyy_0, tpz_z_xxyz_0, tpz_z_xxzz_0, tpz_z_xyyy_0, tpz_z_xyyz_0, tpz_z_xyzz_0, \
                                     tpz_z_xzzz_0, tpz_z_yyyy_0, tpz_z_yyyz_0, tpz_z_yyzz_0, tpz_z_yzzz_0, tpz_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xz_xxxx_0[j] = pa_x[j] * tlx_z_xxxx_0[j] + 2.0 * fl1_fx * tlx_z_xxx_0[j];

                tly_xz_xxxx_0[j] = pa_x[j] * tly_z_xxxx_0[j] + 2.0 * fl1_fx * tly_z_xxx_0[j] + 0.5 * fl1_fx * tpz_z_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_z_xxxx_0[j];

                tlz_xz_xxxx_0[j] = pa_x[j] * tlz_z_xxxx_0[j] + 2.0 * fl1_fx * tlz_z_xxx_0[j] - 0.5 * fl1_fx * tpy_z_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_z_xxxx_0[j];

                tlx_xz_xxxy_0[j] = pa_x[j] * tlx_z_xxxy_0[j] + 1.5 * fl1_fx * tlx_z_xxy_0[j];

                tly_xz_xxxy_0[j] = pa_x[j] * tly_z_xxxy_0[j] + 1.5 * fl1_fx * tly_z_xxy_0[j] + 0.5 * fl1_fx * tpz_z_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_z_xxxy_0[j];

                tlz_xz_xxxy_0[j] = pa_x[j] * tlz_z_xxxy_0[j] + 1.5 * fl1_fx * tlz_z_xxy_0[j] - 0.5 * fl1_fx * tpy_z_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_z_xxxy_0[j];

                tlx_xz_xxxz_0[j] = pa_x[j] * tlx_z_xxxz_0[j] + 1.5 * fl1_fx * tlx_z_xxz_0[j];

                tly_xz_xxxz_0[j] = pa_x[j] * tly_z_xxxz_0[j] + 1.5 * fl1_fx * tly_z_xxz_0[j] + 0.5 * fl1_fx * tpz_z_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_z_xxxz_0[j];

                tlz_xz_xxxz_0[j] = pa_x[j] * tlz_z_xxxz_0[j] + 1.5 * fl1_fx * tlz_z_xxz_0[j] - 0.5 * fl1_fx * tpy_z_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_z_xxxz_0[j];

                tlx_xz_xxyy_0[j] = pa_x[j] * tlx_z_xxyy_0[j] + fl1_fx * tlx_z_xyy_0[j];

                tly_xz_xxyy_0[j] = pa_x[j] * tly_z_xxyy_0[j] + fl1_fx * tly_z_xyy_0[j] + 0.5 * fl1_fx * tpz_z_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_z_xxyy_0[j];

                tlz_xz_xxyy_0[j] = pa_x[j] * tlz_z_xxyy_0[j] + fl1_fx * tlz_z_xyy_0[j] - 0.5 * fl1_fx * tpy_z_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_z_xxyy_0[j];

                tlx_xz_xxyz_0[j] = pa_x[j] * tlx_z_xxyz_0[j] + fl1_fx * tlx_z_xyz_0[j];

                tly_xz_xxyz_0[j] = pa_x[j] * tly_z_xxyz_0[j] + fl1_fx * tly_z_xyz_0[j] + 0.5 * fl1_fx * tpz_z_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_z_xxyz_0[j];

                tlz_xz_xxyz_0[j] = pa_x[j] * tlz_z_xxyz_0[j] + fl1_fx * tlz_z_xyz_0[j] - 0.5 * fl1_fx * tpy_z_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_z_xxyz_0[j];

                tlx_xz_xxzz_0[j] = pa_x[j] * tlx_z_xxzz_0[j] + fl1_fx * tlx_z_xzz_0[j];

                tly_xz_xxzz_0[j] = pa_x[j] * tly_z_xxzz_0[j] + fl1_fx * tly_z_xzz_0[j] + 0.5 * fl1_fx * tpz_z_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_z_xxzz_0[j];

                tlz_xz_xxzz_0[j] = pa_x[j] * tlz_z_xxzz_0[j] + fl1_fx * tlz_z_xzz_0[j] - 0.5 * fl1_fx * tpy_z_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_z_xxzz_0[j];

                tlx_xz_xyyy_0[j] = pa_x[j] * tlx_z_xyyy_0[j] + 0.5 * fl1_fx * tlx_z_yyy_0[j];

                tly_xz_xyyy_0[j] = pa_x[j] * tly_z_xyyy_0[j] + 0.5 * fl1_fx * tly_z_yyy_0[j] + 0.5 * fl1_fx * tpz_z_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_z_xyyy_0[j];

                tlz_xz_xyyy_0[j] = pa_x[j] * tlz_z_xyyy_0[j] + 0.5 * fl1_fx * tlz_z_yyy_0[j] - 0.5 * fl1_fx * tpy_z_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_z_xyyy_0[j];

                tlx_xz_xyyz_0[j] = pa_x[j] * tlx_z_xyyz_0[j] + 0.5 * fl1_fx * tlx_z_yyz_0[j];

                tly_xz_xyyz_0[j] = pa_x[j] * tly_z_xyyz_0[j] + 0.5 * fl1_fx * tly_z_yyz_0[j] + 0.5 * fl1_fx * tpz_z_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_z_xyyz_0[j];

                tlz_xz_xyyz_0[j] = pa_x[j] * tlz_z_xyyz_0[j] + 0.5 * fl1_fx * tlz_z_yyz_0[j] - 0.5 * fl1_fx * tpy_z_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_z_xyyz_0[j];

                tlx_xz_xyzz_0[j] = pa_x[j] * tlx_z_xyzz_0[j] + 0.5 * fl1_fx * tlx_z_yzz_0[j];

                tly_xz_xyzz_0[j] = pa_x[j] * tly_z_xyzz_0[j] + 0.5 * fl1_fx * tly_z_yzz_0[j] + 0.5 * fl1_fx * tpz_z_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_z_xyzz_0[j];

                tlz_xz_xyzz_0[j] = pa_x[j] * tlz_z_xyzz_0[j] + 0.5 * fl1_fx * tlz_z_yzz_0[j] - 0.5 * fl1_fx * tpy_z_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_z_xyzz_0[j];

                tlx_xz_xzzz_0[j] = pa_x[j] * tlx_z_xzzz_0[j] + 0.5 * fl1_fx * tlx_z_zzz_0[j];

                tly_xz_xzzz_0[j] = pa_x[j] * tly_z_xzzz_0[j] + 0.5 * fl1_fx * tly_z_zzz_0[j] + 0.5 * fl1_fx * tpz_z_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_z_xzzz_0[j];

                tlz_xz_xzzz_0[j] = pa_x[j] * tlz_z_xzzz_0[j] + 0.5 * fl1_fx * tlz_z_zzz_0[j] - 0.5 * fl1_fx * tpy_z_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_z_xzzz_0[j];

                tlx_xz_yyyy_0[j] = pa_x[j] * tlx_z_yyyy_0[j];

                tly_xz_yyyy_0[j] = pa_x[j] * tly_z_yyyy_0[j] + 0.5 * fl1_fx * tpz_z_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_z_yyyy_0[j];

                tlz_xz_yyyy_0[j] = pa_x[j] * tlz_z_yyyy_0[j] - 0.5 * fl1_fx * tpy_z_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_z_yyyy_0[j];

                tlx_xz_yyyz_0[j] = pa_x[j] * tlx_z_yyyz_0[j];

                tly_xz_yyyz_0[j] = pa_x[j] * tly_z_yyyz_0[j] + 0.5 * fl1_fx * tpz_z_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_z_yyyz_0[j];

                tlz_xz_yyyz_0[j] = pa_x[j] * tlz_z_yyyz_0[j] - 0.5 * fl1_fx * tpy_z_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_z_yyyz_0[j];

                tlx_xz_yyzz_0[j] = pa_x[j] * tlx_z_yyzz_0[j];

                tly_xz_yyzz_0[j] = pa_x[j] * tly_z_yyzz_0[j] + 0.5 * fl1_fx * tpz_z_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_z_yyzz_0[j];

                tlz_xz_yyzz_0[j] = pa_x[j] * tlz_z_yyzz_0[j] - 0.5 * fl1_fx * tpy_z_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_z_yyzz_0[j];

                tlx_xz_yzzz_0[j] = pa_x[j] * tlx_z_yzzz_0[j];

                tly_xz_yzzz_0[j] = pa_x[j] * tly_z_yzzz_0[j] + 0.5 * fl1_fx * tpz_z_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_z_yzzz_0[j];

                tlz_xz_yzzz_0[j] = pa_x[j] * tlz_z_yzzz_0[j] - 0.5 * fl1_fx * tpy_z_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_z_yzzz_0[j];

                tlx_xz_zzzz_0[j] = pa_x[j] * tlx_z_zzzz_0[j];

                tly_xz_zzzz_0[j] = pa_x[j] * tly_z_zzzz_0[j] + 0.5 * fl1_fx * tpz_z_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_z_zzzz_0[j];

                tlz_xz_zzzz_0[j] = pa_x[j] * tlz_z_zzzz_0[j] - 0.5 * fl1_fx * tpy_z_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_z_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDG_135_180(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tlx_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 15); 

            auto tly_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 15); 

            auto tlz_y_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tlx_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 16); 

            auto tly_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 16); 

            auto tlz_y_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tlx_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 17); 

            auto tly_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 17); 

            auto tlz_y_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tlx_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 18); 

            auto tly_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 18); 

            auto tlz_y_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tlx_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 19); 

            auto tly_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 19); 

            auto tlz_y_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tlx_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 20); 

            auto tly_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 20); 

            auto tlz_y_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tlx_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 21); 

            auto tly_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 21); 

            auto tlz_y_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tlx_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 22); 

            auto tly_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 22); 

            auto tlz_y_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tlx_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 23); 

            auto tly_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 23); 

            auto tlz_y_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tlx_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 24); 

            auto tly_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 24); 

            auto tlz_y_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tlx_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 25); 

            auto tly_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 25); 

            auto tlz_y_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tlx_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 26); 

            auto tly_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 26); 

            auto tlz_y_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tlx_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 27); 

            auto tly_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 27); 

            auto tlz_y_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tlx_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 28); 

            auto tly_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 28); 

            auto tlz_y_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tlx_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 29); 

            auto tly_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 29); 

            auto tlz_y_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            auto tlx_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx); 

            auto tly_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx); 

            auto tlz_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx); 

            auto tlx_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 1); 

            auto tly_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 1); 

            auto tlz_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 1); 

            auto tlx_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 2); 

            auto tly_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 2); 

            auto tlz_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 2); 

            auto tlx_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 3); 

            auto tly_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 3); 

            auto tlz_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 3); 

            auto tlx_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 4); 

            auto tly_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 4); 

            auto tlz_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 4); 

            auto tlx_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 5); 

            auto tly_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 5); 

            auto tlz_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 5); 

            auto tlx_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 6); 

            auto tly_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 6); 

            auto tlz_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 6); 

            auto tlx_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 7); 

            auto tly_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 7); 

            auto tlz_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 7); 

            auto tlx_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 8); 

            auto tly_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 8); 

            auto tlz_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 8); 

            auto tlx_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 9); 

            auto tly_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 9); 

            auto tlz_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 9); 

            auto tlx_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 10); 

            auto tly_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 10); 

            auto tlz_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 10); 

            auto tlx_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 11); 

            auto tly_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 11); 

            auto tlz_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 11); 

            auto tlx_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 12); 

            auto tly_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 12); 

            auto tlz_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 12); 

            auto tlx_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 13); 

            auto tly_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 13); 

            auto tlz_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 13); 

            auto tlx_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 14); 

            auto tly_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 14); 

            auto tlz_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 14); 

            auto tlx_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 10); 

            auto tly_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 10); 

            auto tlz_y_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 10); 

            auto tlx_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 11); 

            auto tly_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 11); 

            auto tlz_y_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 11); 

            auto tlx_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 12); 

            auto tly_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 12); 

            auto tlz_y_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 12); 

            auto tlx_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 13); 

            auto tly_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 13); 

            auto tlz_y_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 13); 

            auto tlx_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 14); 

            auto tly_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 14); 

            auto tlz_y_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 14); 

            auto tlx_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 15); 

            auto tly_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 15); 

            auto tlz_y_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 15); 

            auto tlx_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 16); 

            auto tly_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 16); 

            auto tlz_y_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 16); 

            auto tlx_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 17); 

            auto tly_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 17); 

            auto tlz_y_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 17); 

            auto tlx_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 18); 

            auto tly_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 18); 

            auto tlz_y_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 18); 

            auto tlx_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 19); 

            auto tly_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 19); 

            auto tlz_y_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 19); 

            auto tpx_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 15); 

            auto tpz_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tpx_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 16); 

            auto tpz_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tpx_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 17); 

            auto tpz_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tpx_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 18); 

            auto tpz_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tpx_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 19); 

            auto tpz_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tpx_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 20); 

            auto tpz_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tpx_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 21); 

            auto tpz_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tpx_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 22); 

            auto tpz_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tpx_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 23); 

            auto tpz_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tpx_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 24); 

            auto tpz_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tpx_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 25); 

            auto tpz_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tpx_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 26); 

            auto tpz_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tpx_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 27); 

            auto tpz_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tpx_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 28); 

            auto tpz_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tpx_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 29); 

            auto tpz_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            auto tdx_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 15); 

            auto tdz_y_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 15); 

            auto tdx_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 16); 

            auto tdz_y_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 16); 

            auto tdx_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 17); 

            auto tdz_y_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 17); 

            auto tdx_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 18); 

            auto tdz_y_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 18); 

            auto tdx_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 19); 

            auto tdz_y_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 19); 

            auto tdx_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 20); 

            auto tdz_y_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 20); 

            auto tdx_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 21); 

            auto tdz_y_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 21); 

            auto tdx_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 22); 

            auto tdz_y_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 22); 

            auto tdx_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 23); 

            auto tdz_y_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 23); 

            auto tdx_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 24); 

            auto tdz_y_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 24); 

            auto tdx_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 25); 

            auto tdz_y_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 25); 

            auto tdx_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 26); 

            auto tdz_y_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 26); 

            auto tdx_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 27); 

            auto tdz_y_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 27); 

            auto tdx_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 28); 

            auto tdz_y_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 28); 

            auto tdx_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 29); 

            auto tdz_y_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 29); 

            // set up pointers to integrals

            auto tlx_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 45); 

            auto tly_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 45); 

            auto tlz_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 45); 

            auto tlx_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 46); 

            auto tly_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 46); 

            auto tlz_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 46); 

            auto tlx_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 47); 

            auto tly_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 47); 

            auto tlz_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 47); 

            auto tlx_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 48); 

            auto tly_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 48); 

            auto tlz_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 48); 

            auto tlx_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 49); 

            auto tly_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 49); 

            auto tlz_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 49); 

            auto tlx_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 50); 

            auto tly_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 50); 

            auto tlz_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 50); 

            auto tlx_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 51); 

            auto tly_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 51); 

            auto tlz_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 51); 

            auto tlx_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 52); 

            auto tly_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 52); 

            auto tlz_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 52); 

            auto tlx_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 53); 

            auto tly_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 53); 

            auto tlz_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 53); 

            auto tlx_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 54); 

            auto tly_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 54); 

            auto tlz_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 54); 

            auto tlx_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 55); 

            auto tly_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 55); 

            auto tlz_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 55); 

            auto tlx_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 56); 

            auto tly_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 56); 

            auto tlz_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 56); 

            auto tlx_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 57); 

            auto tly_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 57); 

            auto tlz_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 57); 

            auto tlx_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 58); 

            auto tly_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 58); 

            auto tlz_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 58); 

            auto tlx_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 59); 

            auto tly_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 59); 

            auto tlz_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 59); 

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fgb, fx, pa_y, tdx_y_xxxx_0, tdx_y_xxxy_0, tdx_y_xxxz_0, tdx_y_xxyy_0, \
                                     tdx_y_xxyz_0, tdx_y_xxzz_0, tdx_y_xyyy_0, tdx_y_xyyz_0, tdx_y_xyzz_0, tdx_y_xzzz_0, \
                                     tdx_y_yyyy_0, tdx_y_yyyz_0, tdx_y_yyzz_0, tdx_y_yzzz_0, tdx_y_zzzz_0, tdz_y_xxxx_0, \
                                     tdz_y_xxxy_0, tdz_y_xxxz_0, tdz_y_xxyy_0, tdz_y_xxyz_0, tdz_y_xxzz_0, tdz_y_xyyy_0, \
                                     tdz_y_xyyz_0, tdz_y_xyzz_0, tdz_y_xzzz_0, tdz_y_yyyy_0, tdz_y_yyyz_0, tdz_y_yyzz_0, \
                                     tdz_y_yzzz_0, tdz_y_zzzz_0, tlx_0_xxxx_0, tlx_0_xxxy_0, tlx_0_xxxz_0, tlx_0_xxyy_0, \
                                     tlx_0_xxyz_0, tlx_0_xxzz_0, tlx_0_xyyy_0, tlx_0_xyyz_0, tlx_0_xyzz_0, tlx_0_xzzz_0, \
                                     tlx_0_yyyy_0, tlx_0_yyyz_0, tlx_0_yyzz_0, tlx_0_yzzz_0, tlx_0_zzzz_0, tlx_y_xxx_0, \
                                     tlx_y_xxxx_0, tlx_y_xxxy_0, tlx_y_xxxz_0, tlx_y_xxy_0, tlx_y_xxyy_0, tlx_y_xxyz_0, \
                                     tlx_y_xxz_0, tlx_y_xxzz_0, tlx_y_xyy_0, tlx_y_xyyy_0, tlx_y_xyyz_0, tlx_y_xyz_0, \
                                     tlx_y_xyzz_0, tlx_y_xzz_0, tlx_y_xzzz_0, tlx_y_yyy_0, tlx_y_yyyy_0, tlx_y_yyyz_0, \
                                     tlx_y_yyz_0, tlx_y_yyzz_0, tlx_y_yzz_0, tlx_y_yzzz_0, tlx_y_zzz_0, tlx_y_zzzz_0, \
                                     tlx_yy_xxxx_0, tlx_yy_xxxy_0, tlx_yy_xxxz_0, tlx_yy_xxyy_0, tlx_yy_xxyz_0, \
                                     tlx_yy_xxzz_0, tlx_yy_xyyy_0, tlx_yy_xyyz_0, tlx_yy_xyzz_0, tlx_yy_xzzz_0, \
                                     tlx_yy_yyyy_0, tlx_yy_yyyz_0, tlx_yy_yyzz_0, tlx_yy_yzzz_0, tlx_yy_zzzz_0, \
                                     tly_0_xxxx_0, tly_0_xxxy_0, tly_0_xxxz_0, tly_0_xxyy_0, tly_0_xxyz_0, tly_0_xxzz_0, \
                                     tly_0_xyyy_0, tly_0_xyyz_0, tly_0_xyzz_0, tly_0_xzzz_0, tly_0_yyyy_0, tly_0_yyyz_0, \
                                     tly_0_yyzz_0, tly_0_yzzz_0, tly_0_zzzz_0, tly_y_xxx_0, tly_y_xxxx_0, tly_y_xxxy_0, \
                                     tly_y_xxxz_0, tly_y_xxy_0, tly_y_xxyy_0, tly_y_xxyz_0, tly_y_xxz_0, tly_y_xxzz_0, \
                                     tly_y_xyy_0, tly_y_xyyy_0, tly_y_xyyz_0, tly_y_xyz_0, tly_y_xyzz_0, tly_y_xzz_0, \
                                     tly_y_xzzz_0, tly_y_yyy_0, tly_y_yyyy_0, tly_y_yyyz_0, tly_y_yyz_0, tly_y_yyzz_0, \
                                     tly_y_yzz_0, tly_y_yzzz_0, tly_y_zzz_0, tly_y_zzzz_0, tly_yy_xxxx_0, \
                                     tly_yy_xxxy_0, tly_yy_xxxz_0, tly_yy_xxyy_0, tly_yy_xxyz_0, tly_yy_xxzz_0, \
                                     tly_yy_xyyy_0, tly_yy_xyyz_0, tly_yy_xyzz_0, tly_yy_xzzz_0, tly_yy_yyyy_0, \
                                     tly_yy_yyyz_0, tly_yy_yyzz_0, tly_yy_yzzz_0, tly_yy_zzzz_0, tlz_0_xxxx_0, \
                                     tlz_0_xxxy_0, tlz_0_xxxz_0, tlz_0_xxyy_0, tlz_0_xxyz_0, tlz_0_xxzz_0, tlz_0_xyyy_0, \
                                     tlz_0_xyyz_0, tlz_0_xyzz_0, tlz_0_xzzz_0, tlz_0_yyyy_0, tlz_0_yyyz_0, tlz_0_yyzz_0, \
                                     tlz_0_yzzz_0, tlz_0_zzzz_0, tlz_y_xxx_0, tlz_y_xxxx_0, tlz_y_xxxy_0, tlz_y_xxxz_0, \
                                     tlz_y_xxy_0, tlz_y_xxyy_0, tlz_y_xxyz_0, tlz_y_xxz_0, tlz_y_xxzz_0, tlz_y_xyy_0, \
                                     tlz_y_xyyy_0, tlz_y_xyyz_0, tlz_y_xyz_0, tlz_y_xyzz_0, tlz_y_xzz_0, tlz_y_xzzz_0, \
                                     tlz_y_yyy_0, tlz_y_yyyy_0, tlz_y_yyyz_0, tlz_y_yyz_0, tlz_y_yyzz_0, tlz_y_yzz_0, \
                                     tlz_y_yzzz_0, tlz_y_zzz_0, tlz_y_zzzz_0, tlz_yy_xxxx_0, tlz_yy_xxxy_0, \
                                     tlz_yy_xxxz_0, tlz_yy_xxyy_0, tlz_yy_xxyz_0, tlz_yy_xxzz_0, tlz_yy_xyyy_0, \
                                     tlz_yy_xyyz_0, tlz_yy_xyzz_0, tlz_yy_xzzz_0, tlz_yy_yyyy_0, tlz_yy_yyyz_0, \
                                     tlz_yy_yyzz_0, tlz_yy_yzzz_0, tlz_yy_zzzz_0, tpx_y_xxxx_0, tpx_y_xxxy_0, \
                                     tpx_y_xxxz_0, tpx_y_xxyy_0, tpx_y_xxyz_0, tpx_y_xxzz_0, tpx_y_xyyy_0, tpx_y_xyyz_0, \
                                     tpx_y_xyzz_0, tpx_y_xzzz_0, tpx_y_yyyy_0, tpx_y_yyyz_0, tpx_y_yyzz_0, tpx_y_yzzz_0, \
                                     tpx_y_zzzz_0, tpz_y_xxxx_0, tpz_y_xxxy_0, tpz_y_xxxz_0, tpz_y_xxyy_0, tpz_y_xxyz_0, \
                                     tpz_y_xxzz_0, tpz_y_xyyy_0, tpz_y_xyyz_0, tpz_y_xyzz_0, tpz_y_xzzz_0, tpz_y_yyyy_0, \
                                     tpz_y_yyyz_0, tpz_y_yyzz_0, tpz_y_yzzz_0, tpz_y_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yy_xxxx_0[j] = pa_y[j] * tlx_y_xxxx_0[j] + 0.5 * fl1_fx * tlx_0_xxxx_0[j] - 0.5 * fl1_fx * tpz_y_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_y_xxxx_0[j];

                tly_yy_xxxx_0[j] = pa_y[j] * tly_y_xxxx_0[j] + 0.5 * fl1_fx * tly_0_xxxx_0[j];

                tlz_yy_xxxx_0[j] = pa_y[j] * tlz_y_xxxx_0[j] + 0.5 * fl1_fx * tlz_0_xxxx_0[j] + 0.5 * fl1_fx * tpx_y_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_y_xxxx_0[j];

                tlx_yy_xxxy_0[j] = pa_y[j] * tlx_y_xxxy_0[j] + 0.5 * fl1_fx * tlx_0_xxxy_0[j] + 0.5 * fl1_fx * tlx_y_xxx_0[j] - 0.5 * fl1_fx * tpz_y_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_y_xxxy_0[j];

                tly_yy_xxxy_0[j] = pa_y[j] * tly_y_xxxy_0[j] + 0.5 * fl1_fx * tly_0_xxxy_0[j] + 0.5 * fl1_fx * tly_y_xxx_0[j];

                tlz_yy_xxxy_0[j] = pa_y[j] * tlz_y_xxxy_0[j] + 0.5 * fl1_fx * tlz_0_xxxy_0[j] + 0.5 * fl1_fx * tlz_y_xxx_0[j] + 0.5 * fl1_fx * tpx_y_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_y_xxxy_0[j];

                tlx_yy_xxxz_0[j] = pa_y[j] * tlx_y_xxxz_0[j] + 0.5 * fl1_fx * tlx_0_xxxz_0[j] - 0.5 * fl1_fx * tpz_y_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_y_xxxz_0[j];

                tly_yy_xxxz_0[j] = pa_y[j] * tly_y_xxxz_0[j] + 0.5 * fl1_fx * tly_0_xxxz_0[j];

                tlz_yy_xxxz_0[j] = pa_y[j] * tlz_y_xxxz_0[j] + 0.5 * fl1_fx * tlz_0_xxxz_0[j] + 0.5 * fl1_fx * tpx_y_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_y_xxxz_0[j];

                tlx_yy_xxyy_0[j] = pa_y[j] * tlx_y_xxyy_0[j] + 0.5 * fl1_fx * tlx_0_xxyy_0[j] + fl1_fx * tlx_y_xxy_0[j] - 0.5 * fl1_fx * tpz_y_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_y_xxyy_0[j];

                tly_yy_xxyy_0[j] = pa_y[j] * tly_y_xxyy_0[j] + 0.5 * fl1_fx * tly_0_xxyy_0[j] + fl1_fx * tly_y_xxy_0[j];

                tlz_yy_xxyy_0[j] = pa_y[j] * tlz_y_xxyy_0[j] + 0.5 * fl1_fx * tlz_0_xxyy_0[j] + fl1_fx * tlz_y_xxy_0[j] + 0.5 * fl1_fx * tpx_y_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_y_xxyy_0[j];

                tlx_yy_xxyz_0[j] = pa_y[j] * tlx_y_xxyz_0[j] + 0.5 * fl1_fx * tlx_0_xxyz_0[j] + 0.5 * fl1_fx * tlx_y_xxz_0[j] - 0.5 * fl1_fx * tpz_y_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_y_xxyz_0[j];

                tly_yy_xxyz_0[j] = pa_y[j] * tly_y_xxyz_0[j] + 0.5 * fl1_fx * tly_0_xxyz_0[j] + 0.5 * fl1_fx * tly_y_xxz_0[j];

                tlz_yy_xxyz_0[j] = pa_y[j] * tlz_y_xxyz_0[j] + 0.5 * fl1_fx * tlz_0_xxyz_0[j] + 0.5 * fl1_fx * tlz_y_xxz_0[j] + 0.5 * fl1_fx * tpx_y_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_y_xxyz_0[j];

                tlx_yy_xxzz_0[j] = pa_y[j] * tlx_y_xxzz_0[j] + 0.5 * fl1_fx * tlx_0_xxzz_0[j] - 0.5 * fl1_fx * tpz_y_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_y_xxzz_0[j];

                tly_yy_xxzz_0[j] = pa_y[j] * tly_y_xxzz_0[j] + 0.5 * fl1_fx * tly_0_xxzz_0[j];

                tlz_yy_xxzz_0[j] = pa_y[j] * tlz_y_xxzz_0[j] + 0.5 * fl1_fx * tlz_0_xxzz_0[j] + 0.5 * fl1_fx * tpx_y_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_y_xxzz_0[j];

                tlx_yy_xyyy_0[j] = pa_y[j] * tlx_y_xyyy_0[j] + 0.5 * fl1_fx * tlx_0_xyyy_0[j] + 1.5 * fl1_fx * tlx_y_xyy_0[j] - 0.5 * fl1_fx * tpz_y_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_y_xyyy_0[j];

                tly_yy_xyyy_0[j] = pa_y[j] * tly_y_xyyy_0[j] + 0.5 * fl1_fx * tly_0_xyyy_0[j] + 1.5 * fl1_fx * tly_y_xyy_0[j];

                tlz_yy_xyyy_0[j] = pa_y[j] * tlz_y_xyyy_0[j] + 0.5 * fl1_fx * tlz_0_xyyy_0[j] + 1.5 * fl1_fx * tlz_y_xyy_0[j] + 0.5 * fl1_fx * tpx_y_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_y_xyyy_0[j];

                tlx_yy_xyyz_0[j] = pa_y[j] * tlx_y_xyyz_0[j] + 0.5 * fl1_fx * tlx_0_xyyz_0[j] + fl1_fx * tlx_y_xyz_0[j] - 0.5 * fl1_fx * tpz_y_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_y_xyyz_0[j];

                tly_yy_xyyz_0[j] = pa_y[j] * tly_y_xyyz_0[j] + 0.5 * fl1_fx * tly_0_xyyz_0[j] + fl1_fx * tly_y_xyz_0[j];

                tlz_yy_xyyz_0[j] = pa_y[j] * tlz_y_xyyz_0[j] + 0.5 * fl1_fx * tlz_0_xyyz_0[j] + fl1_fx * tlz_y_xyz_0[j] + 0.5 * fl1_fx * tpx_y_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_y_xyyz_0[j];

                tlx_yy_xyzz_0[j] = pa_y[j] * tlx_y_xyzz_0[j] + 0.5 * fl1_fx * tlx_0_xyzz_0[j] + 0.5 * fl1_fx * tlx_y_xzz_0[j] - 0.5 * fl1_fx * tpz_y_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_y_xyzz_0[j];

                tly_yy_xyzz_0[j] = pa_y[j] * tly_y_xyzz_0[j] + 0.5 * fl1_fx * tly_0_xyzz_0[j] + 0.5 * fl1_fx * tly_y_xzz_0[j];

                tlz_yy_xyzz_0[j] = pa_y[j] * tlz_y_xyzz_0[j] + 0.5 * fl1_fx * tlz_0_xyzz_0[j] + 0.5 * fl1_fx * tlz_y_xzz_0[j] + 0.5 * fl1_fx * tpx_y_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_y_xyzz_0[j];

                tlx_yy_xzzz_0[j] = pa_y[j] * tlx_y_xzzz_0[j] + 0.5 * fl1_fx * tlx_0_xzzz_0[j] - 0.5 * fl1_fx * tpz_y_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_y_xzzz_0[j];

                tly_yy_xzzz_0[j] = pa_y[j] * tly_y_xzzz_0[j] + 0.5 * fl1_fx * tly_0_xzzz_0[j];

                tlz_yy_xzzz_0[j] = pa_y[j] * tlz_y_xzzz_0[j] + 0.5 * fl1_fx * tlz_0_xzzz_0[j] + 0.5 * fl1_fx * tpx_y_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_y_xzzz_0[j];

                tlx_yy_yyyy_0[j] = pa_y[j] * tlx_y_yyyy_0[j] + 0.5 * fl1_fx * tlx_0_yyyy_0[j] + 2.0 * fl1_fx * tlx_y_yyy_0[j] - 0.5 * fl1_fx * tpz_y_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_y_yyyy_0[j];

                tly_yy_yyyy_0[j] = pa_y[j] * tly_y_yyyy_0[j] + 0.5 * fl1_fx * tly_0_yyyy_0[j] + 2.0 * fl1_fx * tly_y_yyy_0[j];

                tlz_yy_yyyy_0[j] = pa_y[j] * tlz_y_yyyy_0[j] + 0.5 * fl1_fx * tlz_0_yyyy_0[j] + 2.0 * fl1_fx * tlz_y_yyy_0[j] + 0.5 * fl1_fx * tpx_y_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_y_yyyy_0[j];

                tlx_yy_yyyz_0[j] = pa_y[j] * tlx_y_yyyz_0[j] + 0.5 * fl1_fx * tlx_0_yyyz_0[j] + 1.5 * fl1_fx * tlx_y_yyz_0[j] - 0.5 * fl1_fx * tpz_y_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_y_yyyz_0[j];

                tly_yy_yyyz_0[j] = pa_y[j] * tly_y_yyyz_0[j] + 0.5 * fl1_fx * tly_0_yyyz_0[j] + 1.5 * fl1_fx * tly_y_yyz_0[j];

                tlz_yy_yyyz_0[j] = pa_y[j] * tlz_y_yyyz_0[j] + 0.5 * fl1_fx * tlz_0_yyyz_0[j] + 1.5 * fl1_fx * tlz_y_yyz_0[j] + 0.5 * fl1_fx * tpx_y_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_y_yyyz_0[j];

                tlx_yy_yyzz_0[j] = pa_y[j] * tlx_y_yyzz_0[j] + 0.5 * fl1_fx * tlx_0_yyzz_0[j] + fl1_fx * tlx_y_yzz_0[j] - 0.5 * fl1_fx * tpz_y_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_y_yyzz_0[j];

                tly_yy_yyzz_0[j] = pa_y[j] * tly_y_yyzz_0[j] + 0.5 * fl1_fx * tly_0_yyzz_0[j] + fl1_fx * tly_y_yzz_0[j];

                tlz_yy_yyzz_0[j] = pa_y[j] * tlz_y_yyzz_0[j] + 0.5 * fl1_fx * tlz_0_yyzz_0[j] + fl1_fx * tlz_y_yzz_0[j] + 0.5 * fl1_fx * tpx_y_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_y_yyzz_0[j];

                tlx_yy_yzzz_0[j] = pa_y[j] * tlx_y_yzzz_0[j] + 0.5 * fl1_fx * tlx_0_yzzz_0[j] + 0.5 * fl1_fx * tlx_y_zzz_0[j] - 0.5 * fl1_fx * tpz_y_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_y_yzzz_0[j];

                tly_yy_yzzz_0[j] = pa_y[j] * tly_y_yzzz_0[j] + 0.5 * fl1_fx * tly_0_yzzz_0[j] + 0.5 * fl1_fx * tly_y_zzz_0[j];

                tlz_yy_yzzz_0[j] = pa_y[j] * tlz_y_yzzz_0[j] + 0.5 * fl1_fx * tlz_0_yzzz_0[j] + 0.5 * fl1_fx * tlz_y_zzz_0[j] + 0.5 * fl1_fx * tpx_y_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_y_yzzz_0[j];

                tlx_yy_zzzz_0[j] = pa_y[j] * tlx_y_zzzz_0[j] + 0.5 * fl1_fx * tlx_0_zzzz_0[j] - 0.5 * fl1_fx * tpz_y_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_y_zzzz_0[j];

                tly_yy_zzzz_0[j] = pa_y[j] * tly_y_zzzz_0[j] + 0.5 * fl1_fx * tly_0_zzzz_0[j];

                tlz_yy_zzzz_0[j] = pa_y[j] * tlz_y_zzzz_0[j] + 0.5 * fl1_fx * tlz_0_zzzz_0[j] + 0.5 * fl1_fx * tpx_y_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_y_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDG_180_225(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tlx_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 30); 

            auto tly_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tlz_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tlx_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 31); 

            auto tly_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tlz_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tlx_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 32); 

            auto tly_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tlz_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tlx_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 33); 

            auto tly_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tlz_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tlx_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 34); 

            auto tly_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tlz_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tlx_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 35); 

            auto tly_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tlz_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tlx_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 36); 

            auto tly_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tlz_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tlx_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 37); 

            auto tly_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tlz_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tlx_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 38); 

            auto tly_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tlz_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tlx_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 39); 

            auto tly_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tlz_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tlx_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 40); 

            auto tly_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tlz_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tlx_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 41); 

            auto tly_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tlz_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tlx_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 42); 

            auto tly_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tlz_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tlx_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 43); 

            auto tly_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tlz_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tlx_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 44); 

            auto tly_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tlz_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tlx_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 20); 

            auto tly_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 21); 

            auto tly_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 22); 

            auto tly_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 23); 

            auto tly_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tlx_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 24); 

            auto tly_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tlx_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 25); 

            auto tly_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tlz_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tlx_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 26); 

            auto tly_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tlz_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tlx_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 27); 

            auto tly_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tlz_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tlx_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 28); 

            auto tly_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tlz_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tlx_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 29); 

            auto tly_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tlz_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30); 

            auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31); 

            auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32); 

            auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33); 

            auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34); 

            auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35); 

            auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36); 

            auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37); 

            auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38); 

            auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39); 

            auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40); 

            auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41); 

            auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42); 

            auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43); 

            auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44); 

            auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30); 

            auto tdz_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31); 

            auto tdz_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32); 

            auto tdz_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33); 

            auto tdz_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34); 

            auto tdz_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35); 

            auto tdz_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36); 

            auto tdz_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37); 

            auto tdz_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38); 

            auto tdz_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39); 

            auto tdz_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40); 

            auto tdz_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41); 

            auto tdz_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42); 

            auto tdz_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43); 

            auto tdz_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44); 

            auto tdz_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            // set up pointers to integrals

            auto tlx_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 60); 

            auto tly_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 60); 

            auto tlz_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 60); 

            auto tlx_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 61); 

            auto tly_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 61); 

            auto tlz_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 61); 

            auto tlx_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 62); 

            auto tly_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 62); 

            auto tlz_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 62); 

            auto tlx_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 63); 

            auto tly_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 63); 

            auto tlz_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 63); 

            auto tlx_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 64); 

            auto tly_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 64); 

            auto tlz_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 64); 

            auto tlx_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 65); 

            auto tly_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 65); 

            auto tlz_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 65); 

            auto tlx_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 66); 

            auto tly_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 66); 

            auto tlz_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 66); 

            auto tlx_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 67); 

            auto tly_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 67); 

            auto tlz_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 67); 

            auto tlx_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 68); 

            auto tly_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 68); 

            auto tlz_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 68); 

            auto tlx_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 69); 

            auto tly_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 69); 

            auto tlz_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 69); 

            auto tlx_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 70); 

            auto tly_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 70); 

            auto tlz_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 70); 

            auto tlx_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 71); 

            auto tly_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 71); 

            auto tlz_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 71); 

            auto tlx_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 72); 

            auto tly_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 72); 

            auto tlz_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 72); 

            auto tlx_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 73); 

            auto tly_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 73); 

            auto tlz_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 73); 

            auto tlx_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 74); 

            auto tly_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 74); 

            auto tlz_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 74); 

            // Batch of Integrals (180,225)

            #pragma omp simd aligned(fgb, fx, pa_y, tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, tdx_z_xxyy_0, \
                                     tdx_z_xxyz_0, tdx_z_xxzz_0, tdx_z_xyyy_0, tdx_z_xyyz_0, tdx_z_xyzz_0, tdx_z_xzzz_0, \
                                     tdx_z_yyyy_0, tdx_z_yyyz_0, tdx_z_yyzz_0, tdx_z_yzzz_0, tdx_z_zzzz_0, tdz_z_xxxx_0, \
                                     tdz_z_xxxy_0, tdz_z_xxxz_0, tdz_z_xxyy_0, tdz_z_xxyz_0, tdz_z_xxzz_0, tdz_z_xyyy_0, \
                                     tdz_z_xyyz_0, tdz_z_xyzz_0, tdz_z_xzzz_0, tdz_z_yyyy_0, tdz_z_yyyz_0, tdz_z_yyzz_0, \
                                     tdz_z_yzzz_0, tdz_z_zzzz_0, tlx_yz_xxxx_0, tlx_yz_xxxy_0, tlx_yz_xxxz_0, \
                                     tlx_yz_xxyy_0, tlx_yz_xxyz_0, tlx_yz_xxzz_0, tlx_yz_xyyy_0, tlx_yz_xyyz_0, \
                                     tlx_yz_xyzz_0, tlx_yz_xzzz_0, tlx_yz_yyyy_0, tlx_yz_yyyz_0, tlx_yz_yyzz_0, \
                                     tlx_yz_yzzz_0, tlx_yz_zzzz_0, tlx_z_xxx_0, tlx_z_xxxx_0, tlx_z_xxxy_0, tlx_z_xxxz_0, \
                                     tlx_z_xxy_0, tlx_z_xxyy_0, tlx_z_xxyz_0, tlx_z_xxz_0, tlx_z_xxzz_0, tlx_z_xyy_0, \
                                     tlx_z_xyyy_0, tlx_z_xyyz_0, tlx_z_xyz_0, tlx_z_xyzz_0, tlx_z_xzz_0, tlx_z_xzzz_0, \
                                     tlx_z_yyy_0, tlx_z_yyyy_0, tlx_z_yyyz_0, tlx_z_yyz_0, tlx_z_yyzz_0, tlx_z_yzz_0, \
                                     tlx_z_yzzz_0, tlx_z_zzz_0, tlx_z_zzzz_0, tly_yz_xxxx_0, tly_yz_xxxy_0, \
                                     tly_yz_xxxz_0, tly_yz_xxyy_0, tly_yz_xxyz_0, tly_yz_xxzz_0, tly_yz_xyyy_0, \
                                     tly_yz_xyyz_0, tly_yz_xyzz_0, tly_yz_xzzz_0, tly_yz_yyyy_0, tly_yz_yyyz_0, \
                                     tly_yz_yyzz_0, tly_yz_yzzz_0, tly_yz_zzzz_0, tly_z_xxx_0, tly_z_xxxx_0, \
                                     tly_z_xxxy_0, tly_z_xxxz_0, tly_z_xxy_0, tly_z_xxyy_0, tly_z_xxyz_0, tly_z_xxz_0, \
                                     tly_z_xxzz_0, tly_z_xyy_0, tly_z_xyyy_0, tly_z_xyyz_0, tly_z_xyz_0, tly_z_xyzz_0, \
                                     tly_z_xzz_0, tly_z_xzzz_0, tly_z_yyy_0, tly_z_yyyy_0, tly_z_yyyz_0, tly_z_yyz_0, \
                                     tly_z_yyzz_0, tly_z_yzz_0, tly_z_yzzz_0, tly_z_zzz_0, tly_z_zzzz_0, tlz_yz_xxxx_0, \
                                     tlz_yz_xxxy_0, tlz_yz_xxxz_0, tlz_yz_xxyy_0, tlz_yz_xxyz_0, tlz_yz_xxzz_0, \
                                     tlz_yz_xyyy_0, tlz_yz_xyyz_0, tlz_yz_xyzz_0, tlz_yz_xzzz_0, tlz_yz_yyyy_0, \
                                     tlz_yz_yyyz_0, tlz_yz_yyzz_0, tlz_yz_yzzz_0, tlz_yz_zzzz_0, tlz_z_xxx_0, \
                                     tlz_z_xxxx_0, tlz_z_xxxy_0, tlz_z_xxxz_0, tlz_z_xxy_0, tlz_z_xxyy_0, tlz_z_xxyz_0, \
                                     tlz_z_xxz_0, tlz_z_xxzz_0, tlz_z_xyy_0, tlz_z_xyyy_0, tlz_z_xyyz_0, tlz_z_xyz_0, \
                                     tlz_z_xyzz_0, tlz_z_xzz_0, tlz_z_xzzz_0, tlz_z_yyy_0, tlz_z_yyyy_0, tlz_z_yyyz_0, \
                                     tlz_z_yyz_0, tlz_z_yyzz_0, tlz_z_yzz_0, tlz_z_yzzz_0, tlz_z_zzz_0, tlz_z_zzzz_0, \
                                     tpx_z_xxxx_0, tpx_z_xxxy_0, tpx_z_xxxz_0, tpx_z_xxyy_0, tpx_z_xxyz_0, tpx_z_xxzz_0, \
                                     tpx_z_xyyy_0, tpx_z_xyyz_0, tpx_z_xyzz_0, tpx_z_xzzz_0, tpx_z_yyyy_0, tpx_z_yyyz_0, \
                                     tpx_z_yyzz_0, tpx_z_yzzz_0, tpx_z_zzzz_0, tpz_z_xxxx_0, tpz_z_xxxy_0, tpz_z_xxxz_0, \
                                     tpz_z_xxyy_0, tpz_z_xxyz_0, tpz_z_xxzz_0, tpz_z_xyyy_0, tpz_z_xyyz_0, tpz_z_xyzz_0, \
                                     tpz_z_xzzz_0, tpz_z_yyyy_0, tpz_z_yyyz_0, tpz_z_yyzz_0, tpz_z_yzzz_0, tpz_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yz_xxxx_0[j] = pa_y[j] * tlx_z_xxxx_0[j] - 0.5 * fl1_fx * tpz_z_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_z_xxxx_0[j];

                tly_yz_xxxx_0[j] = pa_y[j] * tly_z_xxxx_0[j];

                tlz_yz_xxxx_0[j] = pa_y[j] * tlz_z_xxxx_0[j] + 0.5 * fl1_fx * tpx_z_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_z_xxxx_0[j];

                tlx_yz_xxxy_0[j] = pa_y[j] * tlx_z_xxxy_0[j] + 0.5 * fl1_fx * tlx_z_xxx_0[j] - 0.5 * fl1_fx * tpz_z_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_z_xxxy_0[j];

                tly_yz_xxxy_0[j] = pa_y[j] * tly_z_xxxy_0[j] + 0.5 * fl1_fx * tly_z_xxx_0[j];

                tlz_yz_xxxy_0[j] = pa_y[j] * tlz_z_xxxy_0[j] + 0.5 * fl1_fx * tlz_z_xxx_0[j] + 0.5 * fl1_fx * tpx_z_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_z_xxxy_0[j];

                tlx_yz_xxxz_0[j] = pa_y[j] * tlx_z_xxxz_0[j] - 0.5 * fl1_fx * tpz_z_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_z_xxxz_0[j];

                tly_yz_xxxz_0[j] = pa_y[j] * tly_z_xxxz_0[j];

                tlz_yz_xxxz_0[j] = pa_y[j] * tlz_z_xxxz_0[j] + 0.5 * fl1_fx * tpx_z_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_z_xxxz_0[j];

                tlx_yz_xxyy_0[j] = pa_y[j] * tlx_z_xxyy_0[j] + fl1_fx * tlx_z_xxy_0[j] - 0.5 * fl1_fx * tpz_z_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_z_xxyy_0[j];

                tly_yz_xxyy_0[j] = pa_y[j] * tly_z_xxyy_0[j] + fl1_fx * tly_z_xxy_0[j];

                tlz_yz_xxyy_0[j] = pa_y[j] * tlz_z_xxyy_0[j] + fl1_fx * tlz_z_xxy_0[j] + 0.5 * fl1_fx * tpx_z_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_z_xxyy_0[j];

                tlx_yz_xxyz_0[j] = pa_y[j] * tlx_z_xxyz_0[j] + 0.5 * fl1_fx * tlx_z_xxz_0[j] - 0.5 * fl1_fx * tpz_z_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_z_xxyz_0[j];

                tly_yz_xxyz_0[j] = pa_y[j] * tly_z_xxyz_0[j] + 0.5 * fl1_fx * tly_z_xxz_0[j];

                tlz_yz_xxyz_0[j] = pa_y[j] * tlz_z_xxyz_0[j] + 0.5 * fl1_fx * tlz_z_xxz_0[j] + 0.5 * fl1_fx * tpx_z_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_z_xxyz_0[j];

                tlx_yz_xxzz_0[j] = pa_y[j] * tlx_z_xxzz_0[j] - 0.5 * fl1_fx * tpz_z_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_z_xxzz_0[j];

                tly_yz_xxzz_0[j] = pa_y[j] * tly_z_xxzz_0[j];

                tlz_yz_xxzz_0[j] = pa_y[j] * tlz_z_xxzz_0[j] + 0.5 * fl1_fx * tpx_z_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_z_xxzz_0[j];

                tlx_yz_xyyy_0[j] = pa_y[j] * tlx_z_xyyy_0[j] + 1.5 * fl1_fx * tlx_z_xyy_0[j] - 0.5 * fl1_fx * tpz_z_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_z_xyyy_0[j];

                tly_yz_xyyy_0[j] = pa_y[j] * tly_z_xyyy_0[j] + 1.5 * fl1_fx * tly_z_xyy_0[j];

                tlz_yz_xyyy_0[j] = pa_y[j] * tlz_z_xyyy_0[j] + 1.5 * fl1_fx * tlz_z_xyy_0[j] + 0.5 * fl1_fx * tpx_z_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_z_xyyy_0[j];

                tlx_yz_xyyz_0[j] = pa_y[j] * tlx_z_xyyz_0[j] + fl1_fx * tlx_z_xyz_0[j] - 0.5 * fl1_fx * tpz_z_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_z_xyyz_0[j];

                tly_yz_xyyz_0[j] = pa_y[j] * tly_z_xyyz_0[j] + fl1_fx * tly_z_xyz_0[j];

                tlz_yz_xyyz_0[j] = pa_y[j] * tlz_z_xyyz_0[j] + fl1_fx * tlz_z_xyz_0[j] + 0.5 * fl1_fx * tpx_z_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_z_xyyz_0[j];

                tlx_yz_xyzz_0[j] = pa_y[j] * tlx_z_xyzz_0[j] + 0.5 * fl1_fx * tlx_z_xzz_0[j] - 0.5 * fl1_fx * tpz_z_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_z_xyzz_0[j];

                tly_yz_xyzz_0[j] = pa_y[j] * tly_z_xyzz_0[j] + 0.5 * fl1_fx * tly_z_xzz_0[j];

                tlz_yz_xyzz_0[j] = pa_y[j] * tlz_z_xyzz_0[j] + 0.5 * fl1_fx * tlz_z_xzz_0[j] + 0.5 * fl1_fx * tpx_z_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_z_xyzz_0[j];

                tlx_yz_xzzz_0[j] = pa_y[j] * tlx_z_xzzz_0[j] - 0.5 * fl1_fx * tpz_z_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_z_xzzz_0[j];

                tly_yz_xzzz_0[j] = pa_y[j] * tly_z_xzzz_0[j];

                tlz_yz_xzzz_0[j] = pa_y[j] * tlz_z_xzzz_0[j] + 0.5 * fl1_fx * tpx_z_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_z_xzzz_0[j];

                tlx_yz_yyyy_0[j] = pa_y[j] * tlx_z_yyyy_0[j] + 2.0 * fl1_fx * tlx_z_yyy_0[j] - 0.5 * fl1_fx * tpz_z_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_z_yyyy_0[j];

                tly_yz_yyyy_0[j] = pa_y[j] * tly_z_yyyy_0[j] + 2.0 * fl1_fx * tly_z_yyy_0[j];

                tlz_yz_yyyy_0[j] = pa_y[j] * tlz_z_yyyy_0[j] + 2.0 * fl1_fx * tlz_z_yyy_0[j] + 0.5 * fl1_fx * tpx_z_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_z_yyyy_0[j];

                tlx_yz_yyyz_0[j] = pa_y[j] * tlx_z_yyyz_0[j] + 1.5 * fl1_fx * tlx_z_yyz_0[j] - 0.5 * fl1_fx * tpz_z_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_z_yyyz_0[j];

                tly_yz_yyyz_0[j] = pa_y[j] * tly_z_yyyz_0[j] + 1.5 * fl1_fx * tly_z_yyz_0[j];

                tlz_yz_yyyz_0[j] = pa_y[j] * tlz_z_yyyz_0[j] + 1.5 * fl1_fx * tlz_z_yyz_0[j] + 0.5 * fl1_fx * tpx_z_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_z_yyyz_0[j];

                tlx_yz_yyzz_0[j] = pa_y[j] * tlx_z_yyzz_0[j] + fl1_fx * tlx_z_yzz_0[j] - 0.5 * fl1_fx * tpz_z_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_z_yyzz_0[j];

                tly_yz_yyzz_0[j] = pa_y[j] * tly_z_yyzz_0[j] + fl1_fx * tly_z_yzz_0[j];

                tlz_yz_yyzz_0[j] = pa_y[j] * tlz_z_yyzz_0[j] + fl1_fx * tlz_z_yzz_0[j] + 0.5 * fl1_fx * tpx_z_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_z_yyzz_0[j];

                tlx_yz_yzzz_0[j] = pa_y[j] * tlx_z_yzzz_0[j] + 0.5 * fl1_fx * tlx_z_zzz_0[j] - 0.5 * fl1_fx * tpz_z_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_z_yzzz_0[j];

                tly_yz_yzzz_0[j] = pa_y[j] * tly_z_yzzz_0[j] + 0.5 * fl1_fx * tly_z_zzz_0[j];

                tlz_yz_yzzz_0[j] = pa_y[j] * tlz_z_yzzz_0[j] + 0.5 * fl1_fx * tlz_z_zzz_0[j] + 0.5 * fl1_fx * tpx_z_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_z_yzzz_0[j];

                tlx_yz_zzzz_0[j] = pa_y[j] * tlx_z_zzzz_0[j] - 0.5 * fl1_fx * tpz_z_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_z_zzzz_0[j];

                tly_yz_zzzz_0[j] = pa_y[j] * tly_z_zzzz_0[j];

                tlz_yz_zzzz_0[j] = pa_y[j] * tlz_z_zzzz_0[j] + 0.5 * fl1_fx * tpx_z_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_z_zzzz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForDG_225_270(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_2_4_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {0, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {1, -1, -1, -1}, {4, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {1, -1, -1, -1}, {3, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tlx_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 30); 

            auto tly_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tlz_z_xxxx_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 30); 

            auto tlx_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 31); 

            auto tly_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tlz_z_xxxy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 31); 

            auto tlx_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 32); 

            auto tly_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tlz_z_xxxz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 32); 

            auto tlx_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 33); 

            auto tly_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tlz_z_xxyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 33); 

            auto tlx_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 34); 

            auto tly_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tlz_z_xxyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 34); 

            auto tlx_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 35); 

            auto tly_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tlz_z_xxzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 35); 

            auto tlx_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 36); 

            auto tly_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tlz_z_xyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 36); 

            auto tlx_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 37); 

            auto tly_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tlz_z_xyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 37); 

            auto tlx_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 38); 

            auto tly_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tlz_z_xyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 38); 

            auto tlx_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 39); 

            auto tly_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tlz_z_xzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 39); 

            auto tlx_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 40); 

            auto tly_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tlz_z_yyyy_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 40); 

            auto tlx_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 41); 

            auto tly_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tlz_z_yyyz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 41); 

            auto tlx_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 42); 

            auto tly_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tlz_z_yyzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 42); 

            auto tlx_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 43); 

            auto tly_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tlz_z_yzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 43); 

            auto tlx_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * idx + 44); 

            auto tly_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tlz_z_zzzz_0 = primBuffer.data(pidx_l_1_4_m0 + 90 * bdim + 45 * idx + 44); 

            auto tlx_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx); 

            auto tly_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx); 

            auto tlz_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx); 

            auto tlx_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 1); 

            auto tly_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 1); 

            auto tlz_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 1); 

            auto tlx_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 2); 

            auto tly_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 2); 

            auto tlz_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 2); 

            auto tlx_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 3); 

            auto tly_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 3); 

            auto tlz_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 3); 

            auto tlx_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 4); 

            auto tly_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 4); 

            auto tlz_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 4); 

            auto tlx_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 5); 

            auto tly_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 5); 

            auto tlz_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 5); 

            auto tlx_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 6); 

            auto tly_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 6); 

            auto tlz_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 6); 

            auto tlx_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 7); 

            auto tly_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 7); 

            auto tlz_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 7); 

            auto tlx_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 8); 

            auto tly_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 8); 

            auto tlz_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 8); 

            auto tlx_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 9); 

            auto tly_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 9); 

            auto tlz_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 9); 

            auto tlx_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 10); 

            auto tly_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 10); 

            auto tlz_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 10); 

            auto tlx_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 11); 

            auto tly_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 11); 

            auto tlz_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 11); 

            auto tlx_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 12); 

            auto tly_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 12); 

            auto tlz_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 12); 

            auto tlx_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 13); 

            auto tly_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 13); 

            auto tlz_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 13); 

            auto tlx_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 14); 

            auto tly_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 14); 

            auto tlz_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 14); 

            auto tlx_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 20); 

            auto tly_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_z_xxx_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 21); 

            auto tly_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_z_xxy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 22); 

            auto tly_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_z_xxz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 23); 

            auto tly_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_z_xyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 23); 

            auto tlx_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 24); 

            auto tly_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_z_xyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 24); 

            auto tlx_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 25); 

            auto tly_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 25); 

            auto tlz_z_xzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 25); 

            auto tlx_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 26); 

            auto tly_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 26); 

            auto tlz_z_yyy_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 26); 

            auto tlx_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 27); 

            auto tly_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 27); 

            auto tlz_z_yyz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 27); 

            auto tlx_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 28); 

            auto tly_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 28); 

            auto tlz_z_yzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 28); 

            auto tlx_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * idx + 29); 

            auto tly_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 30 * bdim + 30 * idx + 29); 

            auto tlz_z_zzz_0 = primBuffer.data(pidx_l_1_3_m0 + 60 * bdim + 30 * idx + 29); 

            auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30); 

            auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31); 

            auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32); 

            auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33); 

            auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34); 

            auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35); 

            auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36); 

            auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37); 

            auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38); 

            auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39); 

            auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40); 

            auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41); 

            auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42); 

            auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43); 

            auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44); 

            auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            auto tdx_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 30); 

            auto tdy_z_xxxx_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 30); 

            auto tdx_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 31); 

            auto tdy_z_xxxy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 31); 

            auto tdx_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 32); 

            auto tdy_z_xxxz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 32); 

            auto tdx_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 33); 

            auto tdy_z_xxyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 33); 

            auto tdx_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 34); 

            auto tdy_z_xxyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 34); 

            auto tdx_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 35); 

            auto tdy_z_xxzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 35); 

            auto tdx_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 36); 

            auto tdy_z_xyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 36); 

            auto tdx_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 37); 

            auto tdy_z_xyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 37); 

            auto tdx_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 38); 

            auto tdy_z_xyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 38); 

            auto tdx_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 39); 

            auto tdy_z_xzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 39); 

            auto tdx_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 40); 

            auto tdy_z_yyyy_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 40); 

            auto tdx_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 41); 

            auto tdy_z_yyyz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 41); 

            auto tdx_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 42); 

            auto tdy_z_yyzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 42); 

            auto tdx_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 43); 

            auto tdy_z_yzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 43); 

            auto tdx_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * idx + 44); 

            auto tdy_z_zzzz_0 = primBuffer.data(pidx_d_1_4_m0 + 45 * bdim + 45 * idx + 44); 

            // set up pointers to integrals

            auto tlx_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 75); 

            auto tly_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 75); 

            auto tlz_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 75); 

            auto tlx_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 76); 

            auto tly_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 76); 

            auto tlz_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 76); 

            auto tlx_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 77); 

            auto tly_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 77); 

            auto tlz_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 77); 

            auto tlx_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 78); 

            auto tly_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 78); 

            auto tlz_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 78); 

            auto tlx_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 79); 

            auto tly_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 79); 

            auto tlz_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 79); 

            auto tlx_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 80); 

            auto tly_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 80); 

            auto tlz_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 80); 

            auto tlx_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 81); 

            auto tly_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 81); 

            auto tlz_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 81); 

            auto tlx_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 82); 

            auto tly_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 82); 

            auto tlz_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 82); 

            auto tlx_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 83); 

            auto tly_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 83); 

            auto tlz_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 83); 

            auto tlx_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 84); 

            auto tly_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 84); 

            auto tlz_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 84); 

            auto tlx_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 85); 

            auto tly_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 85); 

            auto tlz_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 85); 

            auto tlx_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 86); 

            auto tly_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 86); 

            auto tlz_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 86); 

            auto tlx_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 87); 

            auto tly_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 87); 

            auto tlz_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 87); 

            auto tlx_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 88); 

            auto tly_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 88); 

            auto tlz_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 88); 

            auto tlx_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 89); 

            auto tly_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 89); 

            auto tlz_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 89); 

            // Batch of Integrals (225,270)

            #pragma omp simd aligned(fgb, fx, pa_z, tdx_z_xxxx_0, tdx_z_xxxy_0, tdx_z_xxxz_0, tdx_z_xxyy_0, \
                                     tdx_z_xxyz_0, tdx_z_xxzz_0, tdx_z_xyyy_0, tdx_z_xyyz_0, tdx_z_xyzz_0, tdx_z_xzzz_0, \
                                     tdx_z_yyyy_0, tdx_z_yyyz_0, tdx_z_yyzz_0, tdx_z_yzzz_0, tdx_z_zzzz_0, tdy_z_xxxx_0, \
                                     tdy_z_xxxy_0, tdy_z_xxxz_0, tdy_z_xxyy_0, tdy_z_xxyz_0, tdy_z_xxzz_0, tdy_z_xyyy_0, \
                                     tdy_z_xyyz_0, tdy_z_xyzz_0, tdy_z_xzzz_0, tdy_z_yyyy_0, tdy_z_yyyz_0, tdy_z_yyzz_0, \
                                     tdy_z_yzzz_0, tdy_z_zzzz_0, tlx_0_xxxx_0, tlx_0_xxxy_0, tlx_0_xxxz_0, tlx_0_xxyy_0, \
                                     tlx_0_xxyz_0, tlx_0_xxzz_0, tlx_0_xyyy_0, tlx_0_xyyz_0, tlx_0_xyzz_0, tlx_0_xzzz_0, \
                                     tlx_0_yyyy_0, tlx_0_yyyz_0, tlx_0_yyzz_0, tlx_0_yzzz_0, tlx_0_zzzz_0, tlx_z_xxx_0, \
                                     tlx_z_xxxx_0, tlx_z_xxxy_0, tlx_z_xxxz_0, tlx_z_xxy_0, tlx_z_xxyy_0, tlx_z_xxyz_0, \
                                     tlx_z_xxz_0, tlx_z_xxzz_0, tlx_z_xyy_0, tlx_z_xyyy_0, tlx_z_xyyz_0, tlx_z_xyz_0, \
                                     tlx_z_xyzz_0, tlx_z_xzz_0, tlx_z_xzzz_0, tlx_z_yyy_0, tlx_z_yyyy_0, tlx_z_yyyz_0, \
                                     tlx_z_yyz_0, tlx_z_yyzz_0, tlx_z_yzz_0, tlx_z_yzzz_0, tlx_z_zzz_0, tlx_z_zzzz_0, \
                                     tlx_zz_xxxx_0, tlx_zz_xxxy_0, tlx_zz_xxxz_0, tlx_zz_xxyy_0, tlx_zz_xxyz_0, \
                                     tlx_zz_xxzz_0, tlx_zz_xyyy_0, tlx_zz_xyyz_0, tlx_zz_xyzz_0, tlx_zz_xzzz_0, \
                                     tlx_zz_yyyy_0, tlx_zz_yyyz_0, tlx_zz_yyzz_0, tlx_zz_yzzz_0, tlx_zz_zzzz_0, \
                                     tly_0_xxxx_0, tly_0_xxxy_0, tly_0_xxxz_0, tly_0_xxyy_0, tly_0_xxyz_0, tly_0_xxzz_0, \
                                     tly_0_xyyy_0, tly_0_xyyz_0, tly_0_xyzz_0, tly_0_xzzz_0, tly_0_yyyy_0, tly_0_yyyz_0, \
                                     tly_0_yyzz_0, tly_0_yzzz_0, tly_0_zzzz_0, tly_z_xxx_0, tly_z_xxxx_0, tly_z_xxxy_0, \
                                     tly_z_xxxz_0, tly_z_xxy_0, tly_z_xxyy_0, tly_z_xxyz_0, tly_z_xxz_0, tly_z_xxzz_0, \
                                     tly_z_xyy_0, tly_z_xyyy_0, tly_z_xyyz_0, tly_z_xyz_0, tly_z_xyzz_0, tly_z_xzz_0, \
                                     tly_z_xzzz_0, tly_z_yyy_0, tly_z_yyyy_0, tly_z_yyyz_0, tly_z_yyz_0, tly_z_yyzz_0, \
                                     tly_z_yzz_0, tly_z_yzzz_0, tly_z_zzz_0, tly_z_zzzz_0, tly_zz_xxxx_0, \
                                     tly_zz_xxxy_0, tly_zz_xxxz_0, tly_zz_xxyy_0, tly_zz_xxyz_0, tly_zz_xxzz_0, \
                                     tly_zz_xyyy_0, tly_zz_xyyz_0, tly_zz_xyzz_0, tly_zz_xzzz_0, tly_zz_yyyy_0, \
                                     tly_zz_yyyz_0, tly_zz_yyzz_0, tly_zz_yzzz_0, tly_zz_zzzz_0, tlz_0_xxxx_0, \
                                     tlz_0_xxxy_0, tlz_0_xxxz_0, tlz_0_xxyy_0, tlz_0_xxyz_0, tlz_0_xxzz_0, tlz_0_xyyy_0, \
                                     tlz_0_xyyz_0, tlz_0_xyzz_0, tlz_0_xzzz_0, tlz_0_yyyy_0, tlz_0_yyyz_0, tlz_0_yyzz_0, \
                                     tlz_0_yzzz_0, tlz_0_zzzz_0, tlz_z_xxx_0, tlz_z_xxxx_0, tlz_z_xxxy_0, tlz_z_xxxz_0, \
                                     tlz_z_xxy_0, tlz_z_xxyy_0, tlz_z_xxyz_0, tlz_z_xxz_0, tlz_z_xxzz_0, tlz_z_xyy_0, \
                                     tlz_z_xyyy_0, tlz_z_xyyz_0, tlz_z_xyz_0, tlz_z_xyzz_0, tlz_z_xzz_0, tlz_z_xzzz_0, \
                                     tlz_z_yyy_0, tlz_z_yyyy_0, tlz_z_yyyz_0, tlz_z_yyz_0, tlz_z_yyzz_0, tlz_z_yzz_0, \
                                     tlz_z_yzzz_0, tlz_z_zzz_0, tlz_z_zzzz_0, tlz_zz_xxxx_0, tlz_zz_xxxy_0, \
                                     tlz_zz_xxxz_0, tlz_zz_xxyy_0, tlz_zz_xxyz_0, tlz_zz_xxzz_0, tlz_zz_xyyy_0, \
                                     tlz_zz_xyyz_0, tlz_zz_xyzz_0, tlz_zz_xzzz_0, tlz_zz_yyyy_0, tlz_zz_yyyz_0, \
                                     tlz_zz_yyzz_0, tlz_zz_yzzz_0, tlz_zz_zzzz_0, tpx_z_xxxx_0, tpx_z_xxxy_0, \
                                     tpx_z_xxxz_0, tpx_z_xxyy_0, tpx_z_xxyz_0, tpx_z_xxzz_0, tpx_z_xyyy_0, tpx_z_xyyz_0, \
                                     tpx_z_xyzz_0, tpx_z_xzzz_0, tpx_z_yyyy_0, tpx_z_yyyz_0, tpx_z_yyzz_0, tpx_z_yzzz_0, \
                                     tpx_z_zzzz_0, tpy_z_xxxx_0, tpy_z_xxxy_0, tpy_z_xxxz_0, tpy_z_xxyy_0, tpy_z_xxyz_0, \
                                     tpy_z_xxzz_0, tpy_z_xyyy_0, tpy_z_xyyz_0, tpy_z_xyzz_0, tpy_z_xzzz_0, tpy_z_yyyy_0, \
                                     tpy_z_yyyz_0, tpy_z_yyzz_0, tpy_z_yzzz_0, tpy_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_zz_xxxx_0[j] = pa_z[j] * tlx_z_xxxx_0[j] + 0.5 * fl1_fx * tlx_0_xxxx_0[j] + 0.5 * fl1_fx * tpy_z_xxxx_0[j] + fl1_fx * fl1_fgb * tdy_z_xxxx_0[j];

                tly_zz_xxxx_0[j] = pa_z[j] * tly_z_xxxx_0[j] + 0.5 * fl1_fx * tly_0_xxxx_0[j] - 0.5 * fl1_fx * tpx_z_xxxx_0[j] - fl1_fx * fl1_fgb * tdx_z_xxxx_0[j];

                tlz_zz_xxxx_0[j] = pa_z[j] * tlz_z_xxxx_0[j] + 0.5 * fl1_fx * tlz_0_xxxx_0[j];

                tlx_zz_xxxy_0[j] = pa_z[j] * tlx_z_xxxy_0[j] + 0.5 * fl1_fx * tlx_0_xxxy_0[j] + 0.5 * fl1_fx * tpy_z_xxxy_0[j] + fl1_fx * fl1_fgb * tdy_z_xxxy_0[j];

                tly_zz_xxxy_0[j] = pa_z[j] * tly_z_xxxy_0[j] + 0.5 * fl1_fx * tly_0_xxxy_0[j] - 0.5 * fl1_fx * tpx_z_xxxy_0[j] - fl1_fx * fl1_fgb * tdx_z_xxxy_0[j];

                tlz_zz_xxxy_0[j] = pa_z[j] * tlz_z_xxxy_0[j] + 0.5 * fl1_fx * tlz_0_xxxy_0[j];

                tlx_zz_xxxz_0[j] = pa_z[j] * tlx_z_xxxz_0[j] + 0.5 * fl1_fx * tlx_0_xxxz_0[j] + 0.5 * fl1_fx * tlx_z_xxx_0[j] + 0.5 * fl1_fx * tpy_z_xxxz_0[j] + fl1_fx * fl1_fgb * tdy_z_xxxz_0[j];

                tly_zz_xxxz_0[j] = pa_z[j] * tly_z_xxxz_0[j] + 0.5 * fl1_fx * tly_0_xxxz_0[j] + 0.5 * fl1_fx * tly_z_xxx_0[j] - 0.5 * fl1_fx * tpx_z_xxxz_0[j] - fl1_fx * fl1_fgb * tdx_z_xxxz_0[j];

                tlz_zz_xxxz_0[j] = pa_z[j] * tlz_z_xxxz_0[j] + 0.5 * fl1_fx * tlz_0_xxxz_0[j] + 0.5 * fl1_fx * tlz_z_xxx_0[j];

                tlx_zz_xxyy_0[j] = pa_z[j] * tlx_z_xxyy_0[j] + 0.5 * fl1_fx * tlx_0_xxyy_0[j] + 0.5 * fl1_fx * tpy_z_xxyy_0[j] + fl1_fx * fl1_fgb * tdy_z_xxyy_0[j];

                tly_zz_xxyy_0[j] = pa_z[j] * tly_z_xxyy_0[j] + 0.5 * fl1_fx * tly_0_xxyy_0[j] - 0.5 * fl1_fx * tpx_z_xxyy_0[j] - fl1_fx * fl1_fgb * tdx_z_xxyy_0[j];

                tlz_zz_xxyy_0[j] = pa_z[j] * tlz_z_xxyy_0[j] + 0.5 * fl1_fx * tlz_0_xxyy_0[j];

                tlx_zz_xxyz_0[j] = pa_z[j] * tlx_z_xxyz_0[j] + 0.5 * fl1_fx * tlx_0_xxyz_0[j] + 0.5 * fl1_fx * tlx_z_xxy_0[j] + 0.5 * fl1_fx * tpy_z_xxyz_0[j] + fl1_fx * fl1_fgb * tdy_z_xxyz_0[j];

                tly_zz_xxyz_0[j] = pa_z[j] * tly_z_xxyz_0[j] + 0.5 * fl1_fx * tly_0_xxyz_0[j] + 0.5 * fl1_fx * tly_z_xxy_0[j] - 0.5 * fl1_fx * tpx_z_xxyz_0[j] - fl1_fx * fl1_fgb * tdx_z_xxyz_0[j];

                tlz_zz_xxyz_0[j] = pa_z[j] * tlz_z_xxyz_0[j] + 0.5 * fl1_fx * tlz_0_xxyz_0[j] + 0.5 * fl1_fx * tlz_z_xxy_0[j];

                tlx_zz_xxzz_0[j] = pa_z[j] * tlx_z_xxzz_0[j] + 0.5 * fl1_fx * tlx_0_xxzz_0[j] + fl1_fx * tlx_z_xxz_0[j] + 0.5 * fl1_fx * tpy_z_xxzz_0[j] + fl1_fx * fl1_fgb * tdy_z_xxzz_0[j];

                tly_zz_xxzz_0[j] = pa_z[j] * tly_z_xxzz_0[j] + 0.5 * fl1_fx * tly_0_xxzz_0[j] + fl1_fx * tly_z_xxz_0[j] - 0.5 * fl1_fx * tpx_z_xxzz_0[j] - fl1_fx * fl1_fgb * tdx_z_xxzz_0[j];

                tlz_zz_xxzz_0[j] = pa_z[j] * tlz_z_xxzz_0[j] + 0.5 * fl1_fx * tlz_0_xxzz_0[j] + fl1_fx * tlz_z_xxz_0[j];

                tlx_zz_xyyy_0[j] = pa_z[j] * tlx_z_xyyy_0[j] + 0.5 * fl1_fx * tlx_0_xyyy_0[j] + 0.5 * fl1_fx * tpy_z_xyyy_0[j] + fl1_fx * fl1_fgb * tdy_z_xyyy_0[j];

                tly_zz_xyyy_0[j] = pa_z[j] * tly_z_xyyy_0[j] + 0.5 * fl1_fx * tly_0_xyyy_0[j] - 0.5 * fl1_fx * tpx_z_xyyy_0[j] - fl1_fx * fl1_fgb * tdx_z_xyyy_0[j];

                tlz_zz_xyyy_0[j] = pa_z[j] * tlz_z_xyyy_0[j] + 0.5 * fl1_fx * tlz_0_xyyy_0[j];

                tlx_zz_xyyz_0[j] = pa_z[j] * tlx_z_xyyz_0[j] + 0.5 * fl1_fx * tlx_0_xyyz_0[j] + 0.5 * fl1_fx * tlx_z_xyy_0[j] + 0.5 * fl1_fx * tpy_z_xyyz_0[j] + fl1_fx * fl1_fgb * tdy_z_xyyz_0[j];

                tly_zz_xyyz_0[j] = pa_z[j] * tly_z_xyyz_0[j] + 0.5 * fl1_fx * tly_0_xyyz_0[j] + 0.5 * fl1_fx * tly_z_xyy_0[j] - 0.5 * fl1_fx * tpx_z_xyyz_0[j] - fl1_fx * fl1_fgb * tdx_z_xyyz_0[j];

                tlz_zz_xyyz_0[j] = pa_z[j] * tlz_z_xyyz_0[j] + 0.5 * fl1_fx * tlz_0_xyyz_0[j] + 0.5 * fl1_fx * tlz_z_xyy_0[j];

                tlx_zz_xyzz_0[j] = pa_z[j] * tlx_z_xyzz_0[j] + 0.5 * fl1_fx * tlx_0_xyzz_0[j] + fl1_fx * tlx_z_xyz_0[j] + 0.5 * fl1_fx * tpy_z_xyzz_0[j] + fl1_fx * fl1_fgb * tdy_z_xyzz_0[j];

                tly_zz_xyzz_0[j] = pa_z[j] * tly_z_xyzz_0[j] + 0.5 * fl1_fx * tly_0_xyzz_0[j] + fl1_fx * tly_z_xyz_0[j] - 0.5 * fl1_fx * tpx_z_xyzz_0[j] - fl1_fx * fl1_fgb * tdx_z_xyzz_0[j];

                tlz_zz_xyzz_0[j] = pa_z[j] * tlz_z_xyzz_0[j] + 0.5 * fl1_fx * tlz_0_xyzz_0[j] + fl1_fx * tlz_z_xyz_0[j];

                tlx_zz_xzzz_0[j] = pa_z[j] * tlx_z_xzzz_0[j] + 0.5 * fl1_fx * tlx_0_xzzz_0[j] + 1.5 * fl1_fx * tlx_z_xzz_0[j] + 0.5 * fl1_fx * tpy_z_xzzz_0[j] + fl1_fx * fl1_fgb * tdy_z_xzzz_0[j];

                tly_zz_xzzz_0[j] = pa_z[j] * tly_z_xzzz_0[j] + 0.5 * fl1_fx * tly_0_xzzz_0[j] + 1.5 * fl1_fx * tly_z_xzz_0[j] - 0.5 * fl1_fx * tpx_z_xzzz_0[j] - fl1_fx * fl1_fgb * tdx_z_xzzz_0[j];

                tlz_zz_xzzz_0[j] = pa_z[j] * tlz_z_xzzz_0[j] + 0.5 * fl1_fx * tlz_0_xzzz_0[j] + 1.5 * fl1_fx * tlz_z_xzz_0[j];

                tlx_zz_yyyy_0[j] = pa_z[j] * tlx_z_yyyy_0[j] + 0.5 * fl1_fx * tlx_0_yyyy_0[j] + 0.5 * fl1_fx * tpy_z_yyyy_0[j] + fl1_fx * fl1_fgb * tdy_z_yyyy_0[j];

                tly_zz_yyyy_0[j] = pa_z[j] * tly_z_yyyy_0[j] + 0.5 * fl1_fx * tly_0_yyyy_0[j] - 0.5 * fl1_fx * tpx_z_yyyy_0[j] - fl1_fx * fl1_fgb * tdx_z_yyyy_0[j];

                tlz_zz_yyyy_0[j] = pa_z[j] * tlz_z_yyyy_0[j] + 0.5 * fl1_fx * tlz_0_yyyy_0[j];

                tlx_zz_yyyz_0[j] = pa_z[j] * tlx_z_yyyz_0[j] + 0.5 * fl1_fx * tlx_0_yyyz_0[j] + 0.5 * fl1_fx * tlx_z_yyy_0[j] + 0.5 * fl1_fx * tpy_z_yyyz_0[j] + fl1_fx * fl1_fgb * tdy_z_yyyz_0[j];

                tly_zz_yyyz_0[j] = pa_z[j] * tly_z_yyyz_0[j] + 0.5 * fl1_fx * tly_0_yyyz_0[j] + 0.5 * fl1_fx * tly_z_yyy_0[j] - 0.5 * fl1_fx * tpx_z_yyyz_0[j] - fl1_fx * fl1_fgb * tdx_z_yyyz_0[j];

                tlz_zz_yyyz_0[j] = pa_z[j] * tlz_z_yyyz_0[j] + 0.5 * fl1_fx * tlz_0_yyyz_0[j] + 0.5 * fl1_fx * tlz_z_yyy_0[j];

                tlx_zz_yyzz_0[j] = pa_z[j] * tlx_z_yyzz_0[j] + 0.5 * fl1_fx * tlx_0_yyzz_0[j] + fl1_fx * tlx_z_yyz_0[j] + 0.5 * fl1_fx * tpy_z_yyzz_0[j] + fl1_fx * fl1_fgb * tdy_z_yyzz_0[j];

                tly_zz_yyzz_0[j] = pa_z[j] * tly_z_yyzz_0[j] + 0.5 * fl1_fx * tly_0_yyzz_0[j] + fl1_fx * tly_z_yyz_0[j] - 0.5 * fl1_fx * tpx_z_yyzz_0[j] - fl1_fx * fl1_fgb * tdx_z_yyzz_0[j];

                tlz_zz_yyzz_0[j] = pa_z[j] * tlz_z_yyzz_0[j] + 0.5 * fl1_fx * tlz_0_yyzz_0[j] + fl1_fx * tlz_z_yyz_0[j];

                tlx_zz_yzzz_0[j] = pa_z[j] * tlx_z_yzzz_0[j] + 0.5 * fl1_fx * tlx_0_yzzz_0[j] + 1.5 * fl1_fx * tlx_z_yzz_0[j] + 0.5 * fl1_fx * tpy_z_yzzz_0[j] + fl1_fx * fl1_fgb * tdy_z_yzzz_0[j];

                tly_zz_yzzz_0[j] = pa_z[j] * tly_z_yzzz_0[j] + 0.5 * fl1_fx * tly_0_yzzz_0[j] + 1.5 * fl1_fx * tly_z_yzz_0[j] - 0.5 * fl1_fx * tpx_z_yzzz_0[j] - fl1_fx * fl1_fgb * tdx_z_yzzz_0[j];

                tlz_zz_yzzz_0[j] = pa_z[j] * tlz_z_yzzz_0[j] + 0.5 * fl1_fx * tlz_0_yzzz_0[j] + 1.5 * fl1_fx * tlz_z_yzz_0[j];

                tlx_zz_zzzz_0[j] = pa_z[j] * tlx_z_zzzz_0[j] + 0.5 * fl1_fx * tlx_0_zzzz_0[j] + 2.0 * fl1_fx * tlx_z_zzz_0[j] + 0.5 * fl1_fx * tpy_z_zzzz_0[j] + fl1_fx * fl1_fgb * tdy_z_zzzz_0[j];

                tly_zz_zzzz_0[j] = pa_z[j] * tly_z_zzzz_0[j] + 0.5 * fl1_fx * tly_0_zzzz_0[j] + 2.0 * fl1_fx * tly_z_zzz_0[j] - 0.5 * fl1_fx * tpx_z_zzzz_0[j] - fl1_fx * fl1_fgb * tdx_z_zzzz_0[j];

                tlz_zz_zzzz_0[j] = pa_z[j] * tlz_z_zzzz_0[j] + 0.5 * fl1_fx * tlz_0_zzzz_0[j] + 2.0 * fl1_fx * tlz_z_zzz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForGD(      CMemBlock2D<double>& primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
    {
        amomrecfunc::compAngularMomentumForGD_0_45(primBuffer,
                                                   recursionMap,
                                                   osFactors,
                                                   paDistances, 
                                                   braGtoBlock,
                                                   ketGtoBlock,
                                                   iContrGto); 

        amomrecfunc::compAngularMomentumForGD_45_90(primBuffer,
                                                    recursionMap,
                                                    osFactors,
                                                    paDistances, 
                                                    braGtoBlock,
                                                    ketGtoBlock,
                                                    iContrGto); 

        amomrecfunc::compAngularMomentumForGD_90_135(primBuffer,
                                                     recursionMap,
                                                     osFactors,
                                                     paDistances, 
                                                     braGtoBlock,
                                                     ketGtoBlock,
                                                     iContrGto); 

        amomrecfunc::compAngularMomentumForGD_135_180(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        amomrecfunc::compAngularMomentumForGD_180_225(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 

        amomrecfunc::compAngularMomentumForGD_225_270(primBuffer,
                                                      recursionMap,
                                                      osFactors,
                                                      paDistances, 
                                                      braGtoBlock,
                                                      ketGtoBlock,
                                                      iContrGto); 
    }

    void
    compAngularMomentumForGD_0_45(      CMemBlock2D<double>& primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx); 

            auto tly_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx); 

            auto tlz_xxx_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx); 

            auto tlx_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 1); 

            auto tly_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 1); 

            auto tlz_xxx_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 1); 

            auto tlx_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 2); 

            auto tly_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 2); 

            auto tlz_xxx_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 2); 

            auto tlx_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 3); 

            auto tly_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 3); 

            auto tlz_xxx_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 3); 

            auto tlx_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 4); 

            auto tly_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 4); 

            auto tlz_xxx_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 4); 

            auto tlx_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 5); 

            auto tly_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 5); 

            auto tlz_xxx_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 5); 

            auto tlx_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 6); 

            auto tly_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 6); 

            auto tlz_xxy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 6); 

            auto tlx_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 7); 

            auto tly_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 7); 

            auto tlz_xxy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 7); 

            auto tlx_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 8); 

            auto tly_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 8); 

            auto tlz_xxy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 8); 

            auto tlx_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 9); 

            auto tly_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 9); 

            auto tlz_xxy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 9); 

            auto tlx_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 10); 

            auto tly_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 10); 

            auto tlz_xxy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 10); 

            auto tlx_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 11); 

            auto tly_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 11); 

            auto tlz_xxy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 11); 

            auto tlx_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 12); 

            auto tly_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 12); 

            auto tlz_xxz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 12); 

            auto tlx_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 13); 

            auto tly_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 13); 

            auto tlz_xxz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 13); 

            auto tlx_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 14); 

            auto tly_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 14); 

            auto tlz_xxz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 14); 

            auto tlx_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx); 

            auto tly_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx); 

            auto tlz_xx_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx); 

            auto tlx_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 1); 

            auto tly_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 1); 

            auto tlz_xx_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 1); 

            auto tlx_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 2); 

            auto tly_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 2); 

            auto tlz_xx_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 2); 

            auto tlx_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 3); 

            auto tly_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 3); 

            auto tlz_xx_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 3); 

            auto tlx_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 4); 

            auto tly_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 4); 

            auto tlz_xx_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 4); 

            auto tlx_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 5); 

            auto tly_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 5); 

            auto tlz_xx_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 5); 

            auto tlx_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 6); 

            auto tly_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 6); 

            auto tlz_xy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 6); 

            auto tlx_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 7); 

            auto tly_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 7); 

            auto tlz_xy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 7); 

            auto tlx_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 8); 

            auto tly_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 8); 

            auto tlz_xy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 8); 

            auto tlx_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 9); 

            auto tly_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 9); 

            auto tlz_xy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 9); 

            auto tlx_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 10); 

            auto tly_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 10); 

            auto tlz_xy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 10); 

            auto tlx_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 11); 

            auto tly_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 11); 

            auto tlz_xy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 11); 

            auto tlx_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 12); 

            auto tly_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 12); 

            auto tlz_xz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 12); 

            auto tlx_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 13); 

            auto tly_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 13); 

            auto tlz_xz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 13); 

            auto tlx_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 14); 

            auto tly_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 14); 

            auto tlz_xz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 14); 

            auto tlx_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx); 

            auto tly_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx); 

            auto tlz_xxx_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx); 

            auto tlx_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 1); 

            auto tly_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 1); 

            auto tlz_xxx_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 1); 

            auto tlx_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 2); 

            auto tly_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 2); 

            auto tlz_xxx_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 2); 

            auto tlx_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 3); 

            auto tly_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 3); 

            auto tlz_xxy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 3); 

            auto tlx_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 4); 

            auto tly_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 4); 

            auto tlz_xxy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 4); 

            auto tlx_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 5); 

            auto tly_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 5); 

            auto tlz_xxy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 5); 

            auto tlx_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 6); 

            auto tly_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 6); 

            auto tlz_xxz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 6); 

            auto tlx_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 7); 

            auto tly_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 7); 

            auto tlz_xxz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 7); 

            auto tlx_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 8); 

            auto tly_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 8); 

            auto tlz_xxz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 8); 

            auto tpy_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx); 

            auto tpz_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx); 

            auto tpy_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 1); 

            auto tpz_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 1); 

            auto tpy_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 2); 

            auto tpz_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 2); 

            auto tpy_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 3); 

            auto tpz_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 3); 

            auto tpy_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 4); 

            auto tpz_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 4); 

            auto tpy_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 5); 

            auto tpz_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 5); 

            auto tpy_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 6); 

            auto tpz_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 6); 

            auto tpy_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 7); 

            auto tpz_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 7); 

            auto tpy_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 8); 

            auto tpz_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 8); 

            auto tpy_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 9); 

            auto tpz_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 9); 

            auto tpy_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 10); 

            auto tpz_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 10); 

            auto tpy_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 11); 

            auto tpz_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 11); 

            auto tpy_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 12); 

            auto tpz_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 12); 

            auto tpy_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 13); 

            auto tpz_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 13); 

            auto tpy_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 14); 

            auto tpz_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 14); 

            auto tdy_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx); 

            auto tdz_xxx_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx); 

            auto tdy_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 1); 

            auto tdz_xxx_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 1); 

            auto tdy_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 2); 

            auto tdz_xxx_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 2); 

            auto tdy_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 3); 

            auto tdz_xxx_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 3); 

            auto tdy_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 4); 

            auto tdz_xxx_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 4); 

            auto tdy_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 5); 

            auto tdz_xxx_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 5); 

            auto tdy_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 6); 

            auto tdz_xxy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 6); 

            auto tdy_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 7); 

            auto tdz_xxy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 7); 

            auto tdy_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 8); 

            auto tdz_xxy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 8); 

            auto tdy_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 9); 

            auto tdz_xxy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 9); 

            auto tdy_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 10); 

            auto tdz_xxy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 10); 

            auto tdy_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 11); 

            auto tdz_xxy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 11); 

            auto tdy_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 12); 

            auto tdz_xxz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 12); 

            auto tdy_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 13); 

            auto tdz_xxz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 13); 

            auto tdy_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 14); 

            auto tdz_xxz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 14); 

            // set up pointers to integrals

            auto tlx_xxxx_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx); 

            auto tly_xxxx_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx); 

            auto tlz_xxxx_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx); 

            auto tlx_xxxx_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 1); 

            auto tly_xxxx_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 1); 

            auto tlz_xxxx_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 1); 

            auto tlx_xxxx_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 2); 

            auto tly_xxxx_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 2); 

            auto tlz_xxxx_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 2); 

            auto tlx_xxxx_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 3); 

            auto tly_xxxx_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 3); 

            auto tlz_xxxx_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 3); 

            auto tlx_xxxx_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 4); 

            auto tly_xxxx_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 4); 

            auto tlz_xxxx_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 4); 

            auto tlx_xxxx_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 5); 

            auto tly_xxxx_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 5); 

            auto tlz_xxxx_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 5); 

            auto tlx_xxxy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 6); 

            auto tly_xxxy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 6); 

            auto tlz_xxxy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 6); 

            auto tlx_xxxy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 7); 

            auto tly_xxxy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 7); 

            auto tlz_xxxy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 7); 

            auto tlx_xxxy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 8); 

            auto tly_xxxy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 8); 

            auto tlz_xxxy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 8); 

            auto tlx_xxxy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 9); 

            auto tly_xxxy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 9); 

            auto tlz_xxxy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 9); 

            auto tlx_xxxy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 10); 

            auto tly_xxxy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 10); 

            auto tlz_xxxy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 10); 

            auto tlx_xxxy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 11); 

            auto tly_xxxy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 11); 

            auto tlz_xxxy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 11); 

            auto tlx_xxxz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 12); 

            auto tly_xxxz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 12); 

            auto tlz_xxxz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 12); 

            auto tlx_xxxz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 13); 

            auto tly_xxxz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 13); 

            auto tlz_xxxz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 13); 

            auto tlx_xxxz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 14); 

            auto tly_xxxz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 14); 

            auto tlz_xxxz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 14); 

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxx_xx_0, tdy_xxx_xy_0, tdy_xxx_xz_0, tdy_xxx_yy_0, \
                                     tdy_xxx_yz_0, tdy_xxx_zz_0, tdy_xxy_xx_0, tdy_xxy_xy_0, tdy_xxy_xz_0, tdy_xxy_yy_0, \
                                     tdy_xxy_yz_0, tdy_xxy_zz_0, tdy_xxz_xx_0, tdy_xxz_xy_0, tdy_xxz_xz_0, tdz_xxx_xx_0, \
                                     tdz_xxx_xy_0, tdz_xxx_xz_0, tdz_xxx_yy_0, tdz_xxx_yz_0, tdz_xxx_zz_0, tdz_xxy_xx_0, \
                                     tdz_xxy_xy_0, tdz_xxy_xz_0, tdz_xxy_yy_0, tdz_xxy_yz_0, tdz_xxy_zz_0, tdz_xxz_xx_0, \
                                     tdz_xxz_xy_0, tdz_xxz_xz_0, tlx_xx_xx_0, tlx_xx_xy_0, tlx_xx_xz_0, tlx_xx_yy_0, \
                                     tlx_xx_yz_0, tlx_xx_zz_0, tlx_xxx_x_0, tlx_xxx_xx_0, tlx_xxx_xy_0, tlx_xxx_xz_0, \
                                     tlx_xxx_y_0, tlx_xxx_yy_0, tlx_xxx_yz_0, tlx_xxx_z_0, tlx_xxx_zz_0, tlx_xxxx_xx_0, \
                                     tlx_xxxx_xy_0, tlx_xxxx_xz_0, tlx_xxxx_yy_0, tlx_xxxx_yz_0, tlx_xxxx_zz_0, \
                                     tlx_xxxy_xx_0, tlx_xxxy_xy_0, tlx_xxxy_xz_0, tlx_xxxy_yy_0, tlx_xxxy_yz_0, \
                                     tlx_xxxy_zz_0, tlx_xxxz_xx_0, tlx_xxxz_xy_0, tlx_xxxz_xz_0, tlx_xxy_x_0, \
                                     tlx_xxy_xx_0, tlx_xxy_xy_0, tlx_xxy_xz_0, tlx_xxy_y_0, tlx_xxy_yy_0, tlx_xxy_yz_0, \
                                     tlx_xxy_z_0, tlx_xxy_zz_0, tlx_xxz_x_0, tlx_xxz_xx_0, tlx_xxz_xy_0, tlx_xxz_xz_0, \
                                     tlx_xxz_y_0, tlx_xxz_z_0, tlx_xy_xx_0, tlx_xy_xy_0, tlx_xy_xz_0, tlx_xy_yy_0, \
                                     tlx_xy_yz_0, tlx_xy_zz_0, tlx_xz_xx_0, tlx_xz_xy_0, tlx_xz_xz_0, tly_xx_xx_0, \
                                     tly_xx_xy_0, tly_xx_xz_0, tly_xx_yy_0, tly_xx_yz_0, tly_xx_zz_0, tly_xxx_x_0, \
                                     tly_xxx_xx_0, tly_xxx_xy_0, tly_xxx_xz_0, tly_xxx_y_0, tly_xxx_yy_0, tly_xxx_yz_0, \
                                     tly_xxx_z_0, tly_xxx_zz_0, tly_xxxx_xx_0, tly_xxxx_xy_0, tly_xxxx_xz_0, \
                                     tly_xxxx_yy_0, tly_xxxx_yz_0, tly_xxxx_zz_0, tly_xxxy_xx_0, tly_xxxy_xy_0, \
                                     tly_xxxy_xz_0, tly_xxxy_yy_0, tly_xxxy_yz_0, tly_xxxy_zz_0, tly_xxxz_xx_0, \
                                     tly_xxxz_xy_0, tly_xxxz_xz_0, tly_xxy_x_0, tly_xxy_xx_0, tly_xxy_xy_0, tly_xxy_xz_0, \
                                     tly_xxy_y_0, tly_xxy_yy_0, tly_xxy_yz_0, tly_xxy_z_0, tly_xxy_zz_0, tly_xxz_x_0, \
                                     tly_xxz_xx_0, tly_xxz_xy_0, tly_xxz_xz_0, tly_xxz_y_0, tly_xxz_z_0, tly_xy_xx_0, \
                                     tly_xy_xy_0, tly_xy_xz_0, tly_xy_yy_0, tly_xy_yz_0, tly_xy_zz_0, tly_xz_xx_0, \
                                     tly_xz_xy_0, tly_xz_xz_0, tlz_xx_xx_0, tlz_xx_xy_0, tlz_xx_xz_0, tlz_xx_yy_0, \
                                     tlz_xx_yz_0, tlz_xx_zz_0, tlz_xxx_x_0, tlz_xxx_xx_0, tlz_xxx_xy_0, tlz_xxx_xz_0, \
                                     tlz_xxx_y_0, tlz_xxx_yy_0, tlz_xxx_yz_0, tlz_xxx_z_0, tlz_xxx_zz_0, tlz_xxxx_xx_0, \
                                     tlz_xxxx_xy_0, tlz_xxxx_xz_0, tlz_xxxx_yy_0, tlz_xxxx_yz_0, tlz_xxxx_zz_0, \
                                     tlz_xxxy_xx_0, tlz_xxxy_xy_0, tlz_xxxy_xz_0, tlz_xxxy_yy_0, tlz_xxxy_yz_0, \
                                     tlz_xxxy_zz_0, tlz_xxxz_xx_0, tlz_xxxz_xy_0, tlz_xxxz_xz_0, tlz_xxy_x_0, \
                                     tlz_xxy_xx_0, tlz_xxy_xy_0, tlz_xxy_xz_0, tlz_xxy_y_0, tlz_xxy_yy_0, tlz_xxy_yz_0, \
                                     tlz_xxy_z_0, tlz_xxy_zz_0, tlz_xxz_x_0, tlz_xxz_xx_0, tlz_xxz_xy_0, tlz_xxz_xz_0, \
                                     tlz_xxz_y_0, tlz_xxz_z_0, tlz_xy_xx_0, tlz_xy_xy_0, tlz_xy_xz_0, tlz_xy_yy_0, \
                                     tlz_xy_yz_0, tlz_xy_zz_0, tlz_xz_xx_0, tlz_xz_xy_0, tlz_xz_xz_0, tpy_xxx_xx_0, \
                                     tpy_xxx_xy_0, tpy_xxx_xz_0, tpy_xxx_yy_0, tpy_xxx_yz_0, tpy_xxx_zz_0, tpy_xxy_xx_0, \
                                     tpy_xxy_xy_0, tpy_xxy_xz_0, tpy_xxy_yy_0, tpy_xxy_yz_0, tpy_xxy_zz_0, tpy_xxz_xx_0, \
                                     tpy_xxz_xy_0, tpy_xxz_xz_0, tpz_xxx_xx_0, tpz_xxx_xy_0, tpz_xxx_xz_0, tpz_xxx_yy_0, \
                                     tpz_xxx_yz_0, tpz_xxx_zz_0, tpz_xxy_xx_0, tpz_xxy_xy_0, tpz_xxy_xz_0, tpz_xxy_yy_0, \
                                     tpz_xxy_yz_0, tpz_xxy_zz_0, tpz_xxz_xx_0, tpz_xxz_xy_0, tpz_xxz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xxxx_xx_0[j] = pa_x[j] * tlx_xxx_xx_0[j] + 1.5 * fl1_fx * tlx_xx_xx_0[j] + fl1_fx * tlx_xxx_x_0[j];

                tly_xxxx_xx_0[j] = pa_x[j] * tly_xxx_xx_0[j] + 1.5 * fl1_fx * tly_xx_xx_0[j] + fl1_fx * tly_xxx_x_0[j] + 0.5 * fl1_fx * tpz_xxx_xx_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xx_0[j];

                tlz_xxxx_xx_0[j] = pa_x[j] * tlz_xxx_xx_0[j] + 1.5 * fl1_fx * tlz_xx_xx_0[j] + fl1_fx * tlz_xxx_x_0[j] - 0.5 * fl1_fx * tpy_xxx_xx_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xx_0[j];

                tlx_xxxx_xy_0[j] = pa_x[j] * tlx_xxx_xy_0[j] + 1.5 * fl1_fx * tlx_xx_xy_0[j] + 0.5 * fl1_fx * tlx_xxx_y_0[j];

                tly_xxxx_xy_0[j] = pa_x[j] * tly_xxx_xy_0[j] + 1.5 * fl1_fx * tly_xx_xy_0[j] + 0.5 * fl1_fx * tly_xxx_y_0[j] + 0.5 * fl1_fx * tpz_xxx_xy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xy_0[j];

                tlz_xxxx_xy_0[j] = pa_x[j] * tlz_xxx_xy_0[j] + 1.5 * fl1_fx * tlz_xx_xy_0[j] + 0.5 * fl1_fx * tlz_xxx_y_0[j] - 0.5 * fl1_fx * tpy_xxx_xy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xy_0[j];

                tlx_xxxx_xz_0[j] = pa_x[j] * tlx_xxx_xz_0[j] + 1.5 * fl1_fx * tlx_xx_xz_0[j] + 0.5 * fl1_fx * tlx_xxx_z_0[j];

                tly_xxxx_xz_0[j] = pa_x[j] * tly_xxx_xz_0[j] + 1.5 * fl1_fx * tly_xx_xz_0[j] + 0.5 * fl1_fx * tly_xxx_z_0[j] + 0.5 * fl1_fx * tpz_xxx_xz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xz_0[j];

                tlz_xxxx_xz_0[j] = pa_x[j] * tlz_xxx_xz_0[j] + 1.5 * fl1_fx * tlz_xx_xz_0[j] + 0.5 * fl1_fx * tlz_xxx_z_0[j] - 0.5 * fl1_fx * tpy_xxx_xz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xz_0[j];

                tlx_xxxx_yy_0[j] = pa_x[j] * tlx_xxx_yy_0[j] + 1.5 * fl1_fx * tlx_xx_yy_0[j];

                tly_xxxx_yy_0[j] = pa_x[j] * tly_xxx_yy_0[j] + 1.5 * fl1_fx * tly_xx_yy_0[j] + 0.5 * fl1_fx * tpz_xxx_yy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_yy_0[j];

                tlz_xxxx_yy_0[j] = pa_x[j] * tlz_xxx_yy_0[j] + 1.5 * fl1_fx * tlz_xx_yy_0[j] - 0.5 * fl1_fx * tpy_xxx_yy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_yy_0[j];

                tlx_xxxx_yz_0[j] = pa_x[j] * tlx_xxx_yz_0[j] + 1.5 * fl1_fx * tlx_xx_yz_0[j];

                tly_xxxx_yz_0[j] = pa_x[j] * tly_xxx_yz_0[j] + 1.5 * fl1_fx * tly_xx_yz_0[j] + 0.5 * fl1_fx * tpz_xxx_yz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_yz_0[j];

                tlz_xxxx_yz_0[j] = pa_x[j] * tlz_xxx_yz_0[j] + 1.5 * fl1_fx * tlz_xx_yz_0[j] - 0.5 * fl1_fx * tpy_xxx_yz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_yz_0[j];

                tlx_xxxx_zz_0[j] = pa_x[j] * tlx_xxx_zz_0[j] + 1.5 * fl1_fx * tlx_xx_zz_0[j];

                tly_xxxx_zz_0[j] = pa_x[j] * tly_xxx_zz_0[j] + 1.5 * fl1_fx * tly_xx_zz_0[j] + 0.5 * fl1_fx * tpz_xxx_zz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_zz_0[j];

                tlz_xxxx_zz_0[j] = pa_x[j] * tlz_xxx_zz_0[j] + 1.5 * fl1_fx * tlz_xx_zz_0[j] - 0.5 * fl1_fx * tpy_xxx_zz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_zz_0[j];

                tlx_xxxy_xx_0[j] = pa_x[j] * tlx_xxy_xx_0[j] + fl1_fx * tlx_xy_xx_0[j] + fl1_fx * tlx_xxy_x_0[j];

                tly_xxxy_xx_0[j] = pa_x[j] * tly_xxy_xx_0[j] + fl1_fx * tly_xy_xx_0[j] + fl1_fx * tly_xxy_x_0[j] + 0.5 * fl1_fx * tpz_xxy_xx_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xx_0[j];

                tlz_xxxy_xx_0[j] = pa_x[j] * tlz_xxy_xx_0[j] + fl1_fx * tlz_xy_xx_0[j] + fl1_fx * tlz_xxy_x_0[j] - 0.5 * fl1_fx * tpy_xxy_xx_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xx_0[j];

                tlx_xxxy_xy_0[j] = pa_x[j] * tlx_xxy_xy_0[j] + fl1_fx * tlx_xy_xy_0[j] + 0.5 * fl1_fx * tlx_xxy_y_0[j];

                tly_xxxy_xy_0[j] = pa_x[j] * tly_xxy_xy_0[j] + fl1_fx * tly_xy_xy_0[j] + 0.5 * fl1_fx * tly_xxy_y_0[j] + 0.5 * fl1_fx * tpz_xxy_xy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xy_0[j];

                tlz_xxxy_xy_0[j] = pa_x[j] * tlz_xxy_xy_0[j] + fl1_fx * tlz_xy_xy_0[j] + 0.5 * fl1_fx * tlz_xxy_y_0[j] - 0.5 * fl1_fx * tpy_xxy_xy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xy_0[j];

                tlx_xxxy_xz_0[j] = pa_x[j] * tlx_xxy_xz_0[j] + fl1_fx * tlx_xy_xz_0[j] + 0.5 * fl1_fx * tlx_xxy_z_0[j];

                tly_xxxy_xz_0[j] = pa_x[j] * tly_xxy_xz_0[j] + fl1_fx * tly_xy_xz_0[j] + 0.5 * fl1_fx * tly_xxy_z_0[j] + 0.5 * fl1_fx * tpz_xxy_xz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xz_0[j];

                tlz_xxxy_xz_0[j] = pa_x[j] * tlz_xxy_xz_0[j] + fl1_fx * tlz_xy_xz_0[j] + 0.5 * fl1_fx * tlz_xxy_z_0[j] - 0.5 * fl1_fx * tpy_xxy_xz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xz_0[j];

                tlx_xxxy_yy_0[j] = pa_x[j] * tlx_xxy_yy_0[j] + fl1_fx * tlx_xy_yy_0[j];

                tly_xxxy_yy_0[j] = pa_x[j] * tly_xxy_yy_0[j] + fl1_fx * tly_xy_yy_0[j] + 0.5 * fl1_fx * tpz_xxy_yy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yy_0[j];

                tlz_xxxy_yy_0[j] = pa_x[j] * tlz_xxy_yy_0[j] + fl1_fx * tlz_xy_yy_0[j] - 0.5 * fl1_fx * tpy_xxy_yy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yy_0[j];

                tlx_xxxy_yz_0[j] = pa_x[j] * tlx_xxy_yz_0[j] + fl1_fx * tlx_xy_yz_0[j];

                tly_xxxy_yz_0[j] = pa_x[j] * tly_xxy_yz_0[j] + fl1_fx * tly_xy_yz_0[j] + 0.5 * fl1_fx * tpz_xxy_yz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yz_0[j];

                tlz_xxxy_yz_0[j] = pa_x[j] * tlz_xxy_yz_0[j] + fl1_fx * tlz_xy_yz_0[j] - 0.5 * fl1_fx * tpy_xxy_yz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yz_0[j];

                tlx_xxxy_zz_0[j] = pa_x[j] * tlx_xxy_zz_0[j] + fl1_fx * tlx_xy_zz_0[j];

                tly_xxxy_zz_0[j] = pa_x[j] * tly_xxy_zz_0[j] + fl1_fx * tly_xy_zz_0[j] + 0.5 * fl1_fx * tpz_xxy_zz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_zz_0[j];

                tlz_xxxy_zz_0[j] = pa_x[j] * tlz_xxy_zz_0[j] + fl1_fx * tlz_xy_zz_0[j] - 0.5 * fl1_fx * tpy_xxy_zz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_zz_0[j];

                tlx_xxxz_xx_0[j] = pa_x[j] * tlx_xxz_xx_0[j] + fl1_fx * tlx_xz_xx_0[j] + fl1_fx * tlx_xxz_x_0[j];

                tly_xxxz_xx_0[j] = pa_x[j] * tly_xxz_xx_0[j] + fl1_fx * tly_xz_xx_0[j] + fl1_fx * tly_xxz_x_0[j] + 0.5 * fl1_fx * tpz_xxz_xx_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xx_0[j];

                tlz_xxxz_xx_0[j] = pa_x[j] * tlz_xxz_xx_0[j] + fl1_fx * tlz_xz_xx_0[j] + fl1_fx * tlz_xxz_x_0[j] - 0.5 * fl1_fx * tpy_xxz_xx_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xx_0[j];

                tlx_xxxz_xy_0[j] = pa_x[j] * tlx_xxz_xy_0[j] + fl1_fx * tlx_xz_xy_0[j] + 0.5 * fl1_fx * tlx_xxz_y_0[j];

                tly_xxxz_xy_0[j] = pa_x[j] * tly_xxz_xy_0[j] + fl1_fx * tly_xz_xy_0[j] + 0.5 * fl1_fx * tly_xxz_y_0[j] + 0.5 * fl1_fx * tpz_xxz_xy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xy_0[j];

                tlz_xxxz_xy_0[j] = pa_x[j] * tlz_xxz_xy_0[j] + fl1_fx * tlz_xz_xy_0[j] + 0.5 * fl1_fx * tlz_xxz_y_0[j] - 0.5 * fl1_fx * tpy_xxz_xy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xy_0[j];

                tlx_xxxz_xz_0[j] = pa_x[j] * tlx_xxz_xz_0[j] + fl1_fx * tlx_xz_xz_0[j] + 0.5 * fl1_fx * tlx_xxz_z_0[j];

                tly_xxxz_xz_0[j] = pa_x[j] * tly_xxz_xz_0[j] + fl1_fx * tly_xz_xz_0[j] + 0.5 * fl1_fx * tly_xxz_z_0[j] + 0.5 * fl1_fx * tpz_xxz_xz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xz_0[j];

                tlz_xxxz_xz_0[j] = pa_x[j] * tlz_xxz_xz_0[j] + fl1_fx * tlz_xz_xz_0[j] + 0.5 * fl1_fx * tlz_xxz_z_0[j] - 0.5 * fl1_fx * tpy_xxz_xz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForGD_45_90(      CMemBlock2D<double>& primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 15); 

            auto tly_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 15); 

            auto tlz_xxz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 15); 

            auto tlx_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 16); 

            auto tly_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 16); 

            auto tlz_xxz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 16); 

            auto tlx_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 17); 

            auto tly_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 17); 

            auto tlz_xxz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 17); 

            auto tlx_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 18); 

            auto tly_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 18); 

            auto tlz_xyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 18); 

            auto tlx_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 19); 

            auto tly_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 19); 

            auto tlz_xyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 19); 

            auto tlx_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 20); 

            auto tly_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 20); 

            auto tlz_xyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 20); 

            auto tlx_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 21); 

            auto tly_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 21); 

            auto tlz_xyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 21); 

            auto tlx_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 22); 

            auto tly_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 22); 

            auto tlz_xyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 22); 

            auto tlx_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 23); 

            auto tly_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 23); 

            auto tlz_xyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 23); 

            auto tlx_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 24); 

            auto tly_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 24); 

            auto tlz_xyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 24); 

            auto tlx_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 25); 

            auto tly_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 25); 

            auto tlz_xyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 25); 

            auto tlx_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 26); 

            auto tly_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 26); 

            auto tlz_xyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 26); 

            auto tlx_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 27); 

            auto tly_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 27); 

            auto tlz_xyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 27); 

            auto tlx_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 28); 

            auto tly_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 28); 

            auto tlz_xyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 28); 

            auto tlx_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 29); 

            auto tly_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 29); 

            auto tlz_xyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 29); 

            auto tlx_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 15); 

            auto tly_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 15); 

            auto tlz_xz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 15); 

            auto tlx_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 16); 

            auto tly_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 16); 

            auto tlz_xz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 16); 

            auto tlx_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 17); 

            auto tly_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 17); 

            auto tlz_xz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 17); 

            auto tlx_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 18); 

            auto tly_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tlz_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tlx_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 19); 

            auto tly_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tlz_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tlx_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 20); 

            auto tly_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tlz_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tlx_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 21); 

            auto tly_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tlz_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tlx_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 22); 

            auto tly_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tlz_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tlx_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 23); 

            auto tly_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tlz_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tlx_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 24); 

            auto tly_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tlz_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tlx_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 25); 

            auto tly_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tlz_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tlx_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 26); 

            auto tly_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tlz_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tlx_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 27); 

            auto tly_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tlz_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tlx_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 28); 

            auto tly_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tlz_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tlx_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 29); 

            auto tly_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tlz_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tlx_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 9); 

            auto tly_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 9); 

            auto tlz_xyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 9); 

            auto tlx_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 10); 

            auto tly_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 10); 

            auto tlz_xyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 10); 

            auto tlx_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 11); 

            auto tly_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 11); 

            auto tlz_xyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 11); 

            auto tlx_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 12); 

            auto tly_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 12); 

            auto tlz_xyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 12); 

            auto tlx_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 13); 

            auto tly_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 13); 

            auto tlz_xyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 13); 

            auto tlx_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 14); 

            auto tly_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 14); 

            auto tlz_xyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 14); 

            auto tpy_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 15); 

            auto tpz_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 15); 

            auto tpy_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 16); 

            auto tpz_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 16); 

            auto tpy_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 17); 

            auto tpz_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 17); 

            auto tpy_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 18); 

            auto tpz_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 18); 

            auto tpy_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 19); 

            auto tpz_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 19); 

            auto tpy_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 20); 

            auto tpz_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 20); 

            auto tpy_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 21); 

            auto tpz_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 21); 

            auto tpy_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 22); 

            auto tpz_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 22); 

            auto tpy_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 23); 

            auto tpz_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 23); 

            auto tpy_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 24); 

            auto tpz_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 24); 

            auto tpy_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 25); 

            auto tpz_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 25); 

            auto tpy_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 26); 

            auto tpz_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 26); 

            auto tpy_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 27); 

            auto tpz_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 27); 

            auto tpy_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 28); 

            auto tpz_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 28); 

            auto tpy_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 29); 

            auto tpz_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 29); 

            auto tdy_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 15); 

            auto tdz_xxz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 15); 

            auto tdy_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 16); 

            auto tdz_xxz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 16); 

            auto tdy_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 17); 

            auto tdz_xxz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 17); 

            auto tdy_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 18); 

            auto tdz_xyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 18); 

            auto tdy_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 19); 

            auto tdz_xyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 19); 

            auto tdy_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 20); 

            auto tdz_xyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 20); 

            auto tdy_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 21); 

            auto tdz_xyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 21); 

            auto tdy_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 22); 

            auto tdz_xyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 22); 

            auto tdy_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 23); 

            auto tdz_xyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 23); 

            auto tdy_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 24); 

            auto tdz_xyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 24); 

            auto tdy_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 25); 

            auto tdz_xyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 25); 

            auto tdy_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 26); 

            auto tdz_xyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 26); 

            auto tdy_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 27); 

            auto tdz_xyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 27); 

            auto tdy_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 28); 

            auto tdz_xyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 28); 

            auto tdy_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 29); 

            auto tdz_xyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 29); 

            // set up pointers to integrals

            auto tlx_xxxz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 15); 

            auto tly_xxxz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 15); 

            auto tlz_xxxz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 15); 

            auto tlx_xxxz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 16); 

            auto tly_xxxz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 16); 

            auto tlz_xxxz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 16); 

            auto tlx_xxxz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 17); 

            auto tly_xxxz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 17); 

            auto tlz_xxxz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 17); 

            auto tlx_xxyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 18); 

            auto tly_xxyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 18); 

            auto tlz_xxyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 18); 

            auto tlx_xxyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 19); 

            auto tly_xxyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 19); 

            auto tlz_xxyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 19); 

            auto tlx_xxyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 20); 

            auto tly_xxyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 20); 

            auto tlz_xxyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 20); 

            auto tlx_xxyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 21); 

            auto tly_xxyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 21); 

            auto tlz_xxyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 21); 

            auto tlx_xxyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 22); 

            auto tly_xxyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 22); 

            auto tlz_xxyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 22); 

            auto tlx_xxyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 23); 

            auto tly_xxyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 23); 

            auto tlz_xxyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 23); 

            auto tlx_xxyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 24); 

            auto tly_xxyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 24); 

            auto tlz_xxyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 24); 

            auto tlx_xxyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 25); 

            auto tly_xxyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 25); 

            auto tlz_xxyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 25); 

            auto tlx_xxyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 26); 

            auto tly_xxyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 26); 

            auto tlz_xxyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 26); 

            auto tlx_xxyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 27); 

            auto tly_xxyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 27); 

            auto tlz_xxyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 27); 

            auto tlx_xxyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 28); 

            auto tly_xxyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 28); 

            auto tlz_xxyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 28); 

            auto tlx_xxyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 29); 

            auto tly_xxyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 29); 

            auto tlz_xxyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 29); 

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxz_yy_0, tdy_xxz_yz_0, tdy_xxz_zz_0, tdy_xyy_xx_0, \
                                     tdy_xyy_xy_0, tdy_xyy_xz_0, tdy_xyy_yy_0, tdy_xyy_yz_0, tdy_xyy_zz_0, tdy_xyz_xx_0, \
                                     tdy_xyz_xy_0, tdy_xyz_xz_0, tdy_xyz_yy_0, tdy_xyz_yz_0, tdy_xyz_zz_0, tdz_xxz_yy_0, \
                                     tdz_xxz_yz_0, tdz_xxz_zz_0, tdz_xyy_xx_0, tdz_xyy_xy_0, tdz_xyy_xz_0, tdz_xyy_yy_0, \
                                     tdz_xyy_yz_0, tdz_xyy_zz_0, tdz_xyz_xx_0, tdz_xyz_xy_0, tdz_xyz_xz_0, tdz_xyz_yy_0, \
                                     tdz_xyz_yz_0, tdz_xyz_zz_0, tlx_xxxz_yy_0, tlx_xxxz_yz_0, tlx_xxxz_zz_0, \
                                     tlx_xxyy_xx_0, tlx_xxyy_xy_0, tlx_xxyy_xz_0, tlx_xxyy_yy_0, tlx_xxyy_yz_0, \
                                     tlx_xxyy_zz_0, tlx_xxyz_xx_0, tlx_xxyz_xy_0, tlx_xxyz_xz_0, tlx_xxyz_yy_0, \
                                     tlx_xxyz_yz_0, tlx_xxyz_zz_0, tlx_xxz_yy_0, tlx_xxz_yz_0, tlx_xxz_zz_0, tlx_xyy_x_0, \
                                     tlx_xyy_xx_0, tlx_xyy_xy_0, tlx_xyy_xz_0, tlx_xyy_y_0, tlx_xyy_yy_0, tlx_xyy_yz_0, \
                                     tlx_xyy_z_0, tlx_xyy_zz_0, tlx_xyz_x_0, tlx_xyz_xx_0, tlx_xyz_xy_0, tlx_xyz_xz_0, \
                                     tlx_xyz_y_0, tlx_xyz_yy_0, tlx_xyz_yz_0, tlx_xyz_z_0, tlx_xyz_zz_0, tlx_xz_yy_0, \
                                     tlx_xz_yz_0, tlx_xz_zz_0, tlx_yy_xx_0, tlx_yy_xy_0, tlx_yy_xz_0, tlx_yy_yy_0, \
                                     tlx_yy_yz_0, tlx_yy_zz_0, tlx_yz_xx_0, tlx_yz_xy_0, tlx_yz_xz_0, tlx_yz_yy_0, \
                                     tlx_yz_yz_0, tlx_yz_zz_0, tly_xxxz_yy_0, tly_xxxz_yz_0, tly_xxxz_zz_0, \
                                     tly_xxyy_xx_0, tly_xxyy_xy_0, tly_xxyy_xz_0, tly_xxyy_yy_0, tly_xxyy_yz_0, \
                                     tly_xxyy_zz_0, tly_xxyz_xx_0, tly_xxyz_xy_0, tly_xxyz_xz_0, tly_xxyz_yy_0, \
                                     tly_xxyz_yz_0, tly_xxyz_zz_0, tly_xxz_yy_0, tly_xxz_yz_0, tly_xxz_zz_0, tly_xyy_x_0, \
                                     tly_xyy_xx_0, tly_xyy_xy_0, tly_xyy_xz_0, tly_xyy_y_0, tly_xyy_yy_0, tly_xyy_yz_0, \
                                     tly_xyy_z_0, tly_xyy_zz_0, tly_xyz_x_0, tly_xyz_xx_0, tly_xyz_xy_0, tly_xyz_xz_0, \
                                     tly_xyz_y_0, tly_xyz_yy_0, tly_xyz_yz_0, tly_xyz_z_0, tly_xyz_zz_0, tly_xz_yy_0, \
                                     tly_xz_yz_0, tly_xz_zz_0, tly_yy_xx_0, tly_yy_xy_0, tly_yy_xz_0, tly_yy_yy_0, \
                                     tly_yy_yz_0, tly_yy_zz_0, tly_yz_xx_0, tly_yz_xy_0, tly_yz_xz_0, tly_yz_yy_0, \
                                     tly_yz_yz_0, tly_yz_zz_0, tlz_xxxz_yy_0, tlz_xxxz_yz_0, tlz_xxxz_zz_0, \
                                     tlz_xxyy_xx_0, tlz_xxyy_xy_0, tlz_xxyy_xz_0, tlz_xxyy_yy_0, tlz_xxyy_yz_0, \
                                     tlz_xxyy_zz_0, tlz_xxyz_xx_0, tlz_xxyz_xy_0, tlz_xxyz_xz_0, tlz_xxyz_yy_0, \
                                     tlz_xxyz_yz_0, tlz_xxyz_zz_0, tlz_xxz_yy_0, tlz_xxz_yz_0, tlz_xxz_zz_0, tlz_xyy_x_0, \
                                     tlz_xyy_xx_0, tlz_xyy_xy_0, tlz_xyy_xz_0, tlz_xyy_y_0, tlz_xyy_yy_0, tlz_xyy_yz_0, \
                                     tlz_xyy_z_0, tlz_xyy_zz_0, tlz_xyz_x_0, tlz_xyz_xx_0, tlz_xyz_xy_0, tlz_xyz_xz_0, \
                                     tlz_xyz_y_0, tlz_xyz_yy_0, tlz_xyz_yz_0, tlz_xyz_z_0, tlz_xyz_zz_0, tlz_xz_yy_0, \
                                     tlz_xz_yz_0, tlz_xz_zz_0, tlz_yy_xx_0, tlz_yy_xy_0, tlz_yy_xz_0, tlz_yy_yy_0, \
                                     tlz_yy_yz_0, tlz_yy_zz_0, tlz_yz_xx_0, tlz_yz_xy_0, tlz_yz_xz_0, tlz_yz_yy_0, \
                                     tlz_yz_yz_0, tlz_yz_zz_0, tpy_xxz_yy_0, tpy_xxz_yz_0, tpy_xxz_zz_0, tpy_xyy_xx_0, \
                                     tpy_xyy_xy_0, tpy_xyy_xz_0, tpy_xyy_yy_0, tpy_xyy_yz_0, tpy_xyy_zz_0, tpy_xyz_xx_0, \
                                     tpy_xyz_xy_0, tpy_xyz_xz_0, tpy_xyz_yy_0, tpy_xyz_yz_0, tpy_xyz_zz_0, tpz_xxz_yy_0, \
                                     tpz_xxz_yz_0, tpz_xxz_zz_0, tpz_xyy_xx_0, tpz_xyy_xy_0, tpz_xyy_xz_0, tpz_xyy_yy_0, \
                                     tpz_xyy_yz_0, tpz_xyy_zz_0, tpz_xyz_xx_0, tpz_xyz_xy_0, tpz_xyz_xz_0, tpz_xyz_yy_0, \
                                     tpz_xyz_yz_0, tpz_xyz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xxxz_yy_0[j] = pa_x[j] * tlx_xxz_yy_0[j] + fl1_fx * tlx_xz_yy_0[j];

                tly_xxxz_yy_0[j] = pa_x[j] * tly_xxz_yy_0[j] + fl1_fx * tly_xz_yy_0[j] + 0.5 * fl1_fx * tpz_xxz_yy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yy_0[j];

                tlz_xxxz_yy_0[j] = pa_x[j] * tlz_xxz_yy_0[j] + fl1_fx * tlz_xz_yy_0[j] - 0.5 * fl1_fx * tpy_xxz_yy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yy_0[j];

                tlx_xxxz_yz_0[j] = pa_x[j] * tlx_xxz_yz_0[j] + fl1_fx * tlx_xz_yz_0[j];

                tly_xxxz_yz_0[j] = pa_x[j] * tly_xxz_yz_0[j] + fl1_fx * tly_xz_yz_0[j] + 0.5 * fl1_fx * tpz_xxz_yz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yz_0[j];

                tlz_xxxz_yz_0[j] = pa_x[j] * tlz_xxz_yz_0[j] + fl1_fx * tlz_xz_yz_0[j] - 0.5 * fl1_fx * tpy_xxz_yz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yz_0[j];

                tlx_xxxz_zz_0[j] = pa_x[j] * tlx_xxz_zz_0[j] + fl1_fx * tlx_xz_zz_0[j];

                tly_xxxz_zz_0[j] = pa_x[j] * tly_xxz_zz_0[j] + fl1_fx * tly_xz_zz_0[j] + 0.5 * fl1_fx * tpz_xxz_zz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_zz_0[j];

                tlz_xxxz_zz_0[j] = pa_x[j] * tlz_xxz_zz_0[j] + fl1_fx * tlz_xz_zz_0[j] - 0.5 * fl1_fx * tpy_xxz_zz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_zz_0[j];

                tlx_xxyy_xx_0[j] = pa_x[j] * tlx_xyy_xx_0[j] + 0.5 * fl1_fx * tlx_yy_xx_0[j] + fl1_fx * tlx_xyy_x_0[j];

                tly_xxyy_xx_0[j] = pa_x[j] * tly_xyy_xx_0[j] + 0.5 * fl1_fx * tly_yy_xx_0[j] + fl1_fx * tly_xyy_x_0[j] + 0.5 * fl1_fx * tpz_xyy_xx_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xx_0[j];

                tlz_xxyy_xx_0[j] = pa_x[j] * tlz_xyy_xx_0[j] + 0.5 * fl1_fx * tlz_yy_xx_0[j] + fl1_fx * tlz_xyy_x_0[j] - 0.5 * fl1_fx * tpy_xyy_xx_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xx_0[j];

                tlx_xxyy_xy_0[j] = pa_x[j] * tlx_xyy_xy_0[j] + 0.5 * fl1_fx * tlx_yy_xy_0[j] + 0.5 * fl1_fx * tlx_xyy_y_0[j];

                tly_xxyy_xy_0[j] = pa_x[j] * tly_xyy_xy_0[j] + 0.5 * fl1_fx * tly_yy_xy_0[j] + 0.5 * fl1_fx * tly_xyy_y_0[j] + 0.5 * fl1_fx * tpz_xyy_xy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xy_0[j];

                tlz_xxyy_xy_0[j] = pa_x[j] * tlz_xyy_xy_0[j] + 0.5 * fl1_fx * tlz_yy_xy_0[j] + 0.5 * fl1_fx * tlz_xyy_y_0[j] - 0.5 * fl1_fx * tpy_xyy_xy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xy_0[j];

                tlx_xxyy_xz_0[j] = pa_x[j] * tlx_xyy_xz_0[j] + 0.5 * fl1_fx * tlx_yy_xz_0[j] + 0.5 * fl1_fx * tlx_xyy_z_0[j];

                tly_xxyy_xz_0[j] = pa_x[j] * tly_xyy_xz_0[j] + 0.5 * fl1_fx * tly_yy_xz_0[j] + 0.5 * fl1_fx * tly_xyy_z_0[j] + 0.5 * fl1_fx * tpz_xyy_xz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xz_0[j];

                tlz_xxyy_xz_0[j] = pa_x[j] * tlz_xyy_xz_0[j] + 0.5 * fl1_fx * tlz_yy_xz_0[j] + 0.5 * fl1_fx * tlz_xyy_z_0[j] - 0.5 * fl1_fx * tpy_xyy_xz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xz_0[j];

                tlx_xxyy_yy_0[j] = pa_x[j] * tlx_xyy_yy_0[j] + 0.5 * fl1_fx * tlx_yy_yy_0[j];

                tly_xxyy_yy_0[j] = pa_x[j] * tly_xyy_yy_0[j] + 0.5 * fl1_fx * tly_yy_yy_0[j] + 0.5 * fl1_fx * tpz_xyy_yy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_yy_0[j];

                tlz_xxyy_yy_0[j] = pa_x[j] * tlz_xyy_yy_0[j] + 0.5 * fl1_fx * tlz_yy_yy_0[j] - 0.5 * fl1_fx * tpy_xyy_yy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_yy_0[j];

                tlx_xxyy_yz_0[j] = pa_x[j] * tlx_xyy_yz_0[j] + 0.5 * fl1_fx * tlx_yy_yz_0[j];

                tly_xxyy_yz_0[j] = pa_x[j] * tly_xyy_yz_0[j] + 0.5 * fl1_fx * tly_yy_yz_0[j] + 0.5 * fl1_fx * tpz_xyy_yz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_yz_0[j];

                tlz_xxyy_yz_0[j] = pa_x[j] * tlz_xyy_yz_0[j] + 0.5 * fl1_fx * tlz_yy_yz_0[j] - 0.5 * fl1_fx * tpy_xyy_yz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_yz_0[j];

                tlx_xxyy_zz_0[j] = pa_x[j] * tlx_xyy_zz_0[j] + 0.5 * fl1_fx * tlx_yy_zz_0[j];

                tly_xxyy_zz_0[j] = pa_x[j] * tly_xyy_zz_0[j] + 0.5 * fl1_fx * tly_yy_zz_0[j] + 0.5 * fl1_fx * tpz_xyy_zz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_zz_0[j];

                tlz_xxyy_zz_0[j] = pa_x[j] * tlz_xyy_zz_0[j] + 0.5 * fl1_fx * tlz_yy_zz_0[j] - 0.5 * fl1_fx * tpy_xyy_zz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_zz_0[j];

                tlx_xxyz_xx_0[j] = pa_x[j] * tlx_xyz_xx_0[j] + 0.5 * fl1_fx * tlx_yz_xx_0[j] + fl1_fx * tlx_xyz_x_0[j];

                tly_xxyz_xx_0[j] = pa_x[j] * tly_xyz_xx_0[j] + 0.5 * fl1_fx * tly_yz_xx_0[j] + fl1_fx * tly_xyz_x_0[j] + 0.5 * fl1_fx * tpz_xyz_xx_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xx_0[j];

                tlz_xxyz_xx_0[j] = pa_x[j] * tlz_xyz_xx_0[j] + 0.5 * fl1_fx * tlz_yz_xx_0[j] + fl1_fx * tlz_xyz_x_0[j] - 0.5 * fl1_fx * tpy_xyz_xx_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xx_0[j];

                tlx_xxyz_xy_0[j] = pa_x[j] * tlx_xyz_xy_0[j] + 0.5 * fl1_fx * tlx_yz_xy_0[j] + 0.5 * fl1_fx * tlx_xyz_y_0[j];

                tly_xxyz_xy_0[j] = pa_x[j] * tly_xyz_xy_0[j] + 0.5 * fl1_fx * tly_yz_xy_0[j] + 0.5 * fl1_fx * tly_xyz_y_0[j] + 0.5 * fl1_fx * tpz_xyz_xy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xy_0[j];

                tlz_xxyz_xy_0[j] = pa_x[j] * tlz_xyz_xy_0[j] + 0.5 * fl1_fx * tlz_yz_xy_0[j] + 0.5 * fl1_fx * tlz_xyz_y_0[j] - 0.5 * fl1_fx * tpy_xyz_xy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xy_0[j];

                tlx_xxyz_xz_0[j] = pa_x[j] * tlx_xyz_xz_0[j] + 0.5 * fl1_fx * tlx_yz_xz_0[j] + 0.5 * fl1_fx * tlx_xyz_z_0[j];

                tly_xxyz_xz_0[j] = pa_x[j] * tly_xyz_xz_0[j] + 0.5 * fl1_fx * tly_yz_xz_0[j] + 0.5 * fl1_fx * tly_xyz_z_0[j] + 0.5 * fl1_fx * tpz_xyz_xz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xz_0[j];

                tlz_xxyz_xz_0[j] = pa_x[j] * tlz_xyz_xz_0[j] + 0.5 * fl1_fx * tlz_yz_xz_0[j] + 0.5 * fl1_fx * tlz_xyz_z_0[j] - 0.5 * fl1_fx * tpy_xyz_xz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xz_0[j];

                tlx_xxyz_yy_0[j] = pa_x[j] * tlx_xyz_yy_0[j] + 0.5 * fl1_fx * tlx_yz_yy_0[j];

                tly_xxyz_yy_0[j] = pa_x[j] * tly_xyz_yy_0[j] + 0.5 * fl1_fx * tly_yz_yy_0[j] + 0.5 * fl1_fx * tpz_xyz_yy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_yy_0[j];

                tlz_xxyz_yy_0[j] = pa_x[j] * tlz_xyz_yy_0[j] + 0.5 * fl1_fx * tlz_yz_yy_0[j] - 0.5 * fl1_fx * tpy_xyz_yy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_yy_0[j];

                tlx_xxyz_yz_0[j] = pa_x[j] * tlx_xyz_yz_0[j] + 0.5 * fl1_fx * tlx_yz_yz_0[j];

                tly_xxyz_yz_0[j] = pa_x[j] * tly_xyz_yz_0[j] + 0.5 * fl1_fx * tly_yz_yz_0[j] + 0.5 * fl1_fx * tpz_xyz_yz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_yz_0[j];

                tlz_xxyz_yz_0[j] = pa_x[j] * tlz_xyz_yz_0[j] + 0.5 * fl1_fx * tlz_yz_yz_0[j] - 0.5 * fl1_fx * tpy_xyz_yz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_yz_0[j];

                tlx_xxyz_zz_0[j] = pa_x[j] * tlx_xyz_zz_0[j] + 0.5 * fl1_fx * tlx_yz_zz_0[j];

                tly_xxyz_zz_0[j] = pa_x[j] * tly_xyz_zz_0[j] + 0.5 * fl1_fx * tly_yz_zz_0[j] + 0.5 * fl1_fx * tpz_xyz_zz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_zz_0[j];

                tlz_xxyz_zz_0[j] = pa_x[j] * tlz_xyz_zz_0[j] + 0.5 * fl1_fx * tlz_yz_zz_0[j] - 0.5 * fl1_fx * tpy_xyz_zz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_zz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForGD_90_135(      CMemBlock2D<double>& primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 30); 

            auto tly_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 30); 

            auto tlz_xzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 30); 

            auto tlx_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 31); 

            auto tly_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 31); 

            auto tlz_xzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 31); 

            auto tlx_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 32); 

            auto tly_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 32); 

            auto tlz_xzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 32); 

            auto tlx_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 33); 

            auto tly_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 33); 

            auto tlz_xzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 33); 

            auto tlx_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 34); 

            auto tly_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 34); 

            auto tlz_xzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 34); 

            auto tlx_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 35); 

            auto tly_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 35); 

            auto tlz_xzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 35); 

            auto tlx_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 36); 

            auto tly_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tlz_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tlx_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 37); 

            auto tly_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tlz_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tlx_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 38); 

            auto tly_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tlz_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tlx_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 39); 

            auto tly_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tlz_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tlx_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 40); 

            auto tly_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tlz_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tlx_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 41); 

            auto tly_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tlz_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tlx_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 42); 

            auto tly_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tlz_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tlx_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 43); 

            auto tly_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tlz_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tlx_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 44); 

            auto tly_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tlz_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30); 

            auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31); 

            auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32); 

            auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tlx_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 33); 

            auto tly_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tlz_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tlx_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 34); 

            auto tly_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tlz_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tlx_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 35); 

            auto tly_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tlz_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            auto tlx_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 15); 

            auto tly_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 15); 

            auto tlz_xzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 15); 

            auto tlx_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 16); 

            auto tly_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 16); 

            auto tlz_xzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 16); 

            auto tlx_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 17); 

            auto tly_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 17); 

            auto tlz_xzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 17); 

            auto tlx_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 18); 

            auto tly_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 18); 

            auto tlz_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 18); 

            auto tlx_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 19); 

            auto tly_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 19); 

            auto tlz_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 19); 

            auto tlx_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 20); 

            auto tly_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 21); 

            auto tly_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 22); 

            auto tly_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 23); 

            auto tly_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 23); 

            auto tpy_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 30); 

            auto tpz_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 30); 

            auto tpy_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 31); 

            auto tpz_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 31); 

            auto tpy_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 32); 

            auto tpz_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 32); 

            auto tpy_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 33); 

            auto tpz_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 33); 

            auto tpy_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 34); 

            auto tpz_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 34); 

            auto tpy_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 35); 

            auto tpz_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 35); 

            auto tpy_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tpz_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tpy_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tpz_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tpy_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tpz_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tpy_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tpz_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tpy_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tpz_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tpy_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tpz_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tpy_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tpz_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tpy_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tpz_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tpy_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tpz_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tdy_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 30); 

            auto tdz_xzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 30); 

            auto tdy_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 31); 

            auto tdz_xzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 31); 

            auto tdy_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 32); 

            auto tdz_xzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 32); 

            auto tdy_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 33); 

            auto tdz_xzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 33); 

            auto tdy_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 34); 

            auto tdz_xzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 34); 

            auto tdy_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 35); 

            auto tdz_xzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 35); 

            auto tdy_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tdz_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdy_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tdz_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdy_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tdz_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdy_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tdz_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdy_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tdz_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdy_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tdz_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdy_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tdz_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdy_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tdz_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdy_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tdz_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            // set up pointers to integrals

            auto tlx_xxzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 30); 

            auto tly_xxzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 30); 

            auto tlz_xxzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 30); 

            auto tlx_xxzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 31); 

            auto tly_xxzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 31); 

            auto tlz_xxzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 31); 

            auto tlx_xxzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 32); 

            auto tly_xxzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 32); 

            auto tlz_xxzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 32); 

            auto tlx_xxzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 33); 

            auto tly_xxzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 33); 

            auto tlz_xxzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 33); 

            auto tlx_xxzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 34); 

            auto tly_xxzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 34); 

            auto tlz_xxzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 34); 

            auto tlx_xxzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 35); 

            auto tly_xxzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 35); 

            auto tlz_xxzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 35); 

            auto tlx_xyyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 36); 

            auto tly_xyyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 36); 

            auto tlz_xyyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 36); 

            auto tlx_xyyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 37); 

            auto tly_xyyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 37); 

            auto tlz_xyyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 37); 

            auto tlx_xyyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 38); 

            auto tly_xyyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 38); 

            auto tlz_xyyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 38); 

            auto tlx_xyyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 39); 

            auto tly_xyyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 39); 

            auto tlz_xyyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 39); 

            auto tlx_xyyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 40); 

            auto tly_xyyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 40); 

            auto tlz_xyyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 40); 

            auto tlx_xyyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 41); 

            auto tly_xyyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 41); 

            auto tlz_xyyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 41); 

            auto tlx_xyyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 42); 

            auto tly_xyyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 42); 

            auto tlz_xyyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 42); 

            auto tlx_xyyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 43); 

            auto tly_xyyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 43); 

            auto tlz_xyyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 43); 

            auto tlx_xyyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 44); 

            auto tly_xyyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 44); 

            auto tlz_xyyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 44); 

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_xzz_xx_0, tdy_xzz_xy_0, tdy_xzz_xz_0, tdy_xzz_yy_0, \
                                     tdy_xzz_yz_0, tdy_xzz_zz_0, tdy_yyy_xx_0, tdy_yyy_xy_0, tdy_yyy_xz_0, tdy_yyy_yy_0, \
                                     tdy_yyy_yz_0, tdy_yyy_zz_0, tdy_yyz_xx_0, tdy_yyz_xy_0, tdy_yyz_xz_0, tdz_xzz_xx_0, \
                                     tdz_xzz_xy_0, tdz_xzz_xz_0, tdz_xzz_yy_0, tdz_xzz_yz_0, tdz_xzz_zz_0, tdz_yyy_xx_0, \
                                     tdz_yyy_xy_0, tdz_yyy_xz_0, tdz_yyy_yy_0, tdz_yyy_yz_0, tdz_yyy_zz_0, tdz_yyz_xx_0, \
                                     tdz_yyz_xy_0, tdz_yyz_xz_0, tlx_xxzz_xx_0, tlx_xxzz_xy_0, tlx_xxzz_xz_0, \
                                     tlx_xxzz_yy_0, tlx_xxzz_yz_0, tlx_xxzz_zz_0, tlx_xyyy_xx_0, tlx_xyyy_xy_0, \
                                     tlx_xyyy_xz_0, tlx_xyyy_yy_0, tlx_xyyy_yz_0, tlx_xyyy_zz_0, tlx_xyyz_xx_0, \
                                     tlx_xyyz_xy_0, tlx_xyyz_xz_0, tlx_xzz_x_0, tlx_xzz_xx_0, tlx_xzz_xy_0, tlx_xzz_xz_0, \
                                     tlx_xzz_y_0, tlx_xzz_yy_0, tlx_xzz_yz_0, tlx_xzz_z_0, tlx_xzz_zz_0, tlx_yyy_x_0, \
                                     tlx_yyy_xx_0, tlx_yyy_xy_0, tlx_yyy_xz_0, tlx_yyy_y_0, tlx_yyy_yy_0, tlx_yyy_yz_0, \
                                     tlx_yyy_z_0, tlx_yyy_zz_0, tlx_yyz_x_0, tlx_yyz_xx_0, tlx_yyz_xy_0, tlx_yyz_xz_0, \
                                     tlx_yyz_y_0, tlx_yyz_z_0, tlx_zz_xx_0, tlx_zz_xy_0, tlx_zz_xz_0, tlx_zz_yy_0, \
                                     tlx_zz_yz_0, tlx_zz_zz_0, tly_xxzz_xx_0, tly_xxzz_xy_0, tly_xxzz_xz_0, \
                                     tly_xxzz_yy_0, tly_xxzz_yz_0, tly_xxzz_zz_0, tly_xyyy_xx_0, tly_xyyy_xy_0, \
                                     tly_xyyy_xz_0, tly_xyyy_yy_0, tly_xyyy_yz_0, tly_xyyy_zz_0, tly_xyyz_xx_0, \
                                     tly_xyyz_xy_0, tly_xyyz_xz_0, tly_xzz_x_0, tly_xzz_xx_0, tly_xzz_xy_0, tly_xzz_xz_0, \
                                     tly_xzz_y_0, tly_xzz_yy_0, tly_xzz_yz_0, tly_xzz_z_0, tly_xzz_zz_0, tly_yyy_x_0, \
                                     tly_yyy_xx_0, tly_yyy_xy_0, tly_yyy_xz_0, tly_yyy_y_0, tly_yyy_yy_0, tly_yyy_yz_0, \
                                     tly_yyy_z_0, tly_yyy_zz_0, tly_yyz_x_0, tly_yyz_xx_0, tly_yyz_xy_0, tly_yyz_xz_0, \
                                     tly_yyz_y_0, tly_yyz_z_0, tly_zz_xx_0, tly_zz_xy_0, tly_zz_xz_0, tly_zz_yy_0, \
                                     tly_zz_yz_0, tly_zz_zz_0, tlz_xxzz_xx_0, tlz_xxzz_xy_0, tlz_xxzz_xz_0, \
                                     tlz_xxzz_yy_0, tlz_xxzz_yz_0, tlz_xxzz_zz_0, tlz_xyyy_xx_0, tlz_xyyy_xy_0, \
                                     tlz_xyyy_xz_0, tlz_xyyy_yy_0, tlz_xyyy_yz_0, tlz_xyyy_zz_0, tlz_xyyz_xx_0, \
                                     tlz_xyyz_xy_0, tlz_xyyz_xz_0, tlz_xzz_x_0, tlz_xzz_xx_0, tlz_xzz_xy_0, tlz_xzz_xz_0, \
                                     tlz_xzz_y_0, tlz_xzz_yy_0, tlz_xzz_yz_0, tlz_xzz_z_0, tlz_xzz_zz_0, tlz_yyy_x_0, \
                                     tlz_yyy_xx_0, tlz_yyy_xy_0, tlz_yyy_xz_0, tlz_yyy_y_0, tlz_yyy_yy_0, tlz_yyy_yz_0, \
                                     tlz_yyy_z_0, tlz_yyy_zz_0, tlz_yyz_x_0, tlz_yyz_xx_0, tlz_yyz_xy_0, tlz_yyz_xz_0, \
                                     tlz_yyz_y_0, tlz_yyz_z_0, tlz_zz_xx_0, tlz_zz_xy_0, tlz_zz_xz_0, tlz_zz_yy_0, \
                                     tlz_zz_yz_0, tlz_zz_zz_0, tpy_xzz_xx_0, tpy_xzz_xy_0, tpy_xzz_xz_0, tpy_xzz_yy_0, \
                                     tpy_xzz_yz_0, tpy_xzz_zz_0, tpy_yyy_xx_0, tpy_yyy_xy_0, tpy_yyy_xz_0, tpy_yyy_yy_0, \
                                     tpy_yyy_yz_0, tpy_yyy_zz_0, tpy_yyz_xx_0, tpy_yyz_xy_0, tpy_yyz_xz_0, tpz_xzz_xx_0, \
                                     tpz_xzz_xy_0, tpz_xzz_xz_0, tpz_xzz_yy_0, tpz_xzz_yz_0, tpz_xzz_zz_0, tpz_yyy_xx_0, \
                                     tpz_yyy_xy_0, tpz_yyy_xz_0, tpz_yyy_yy_0, tpz_yyy_yz_0, tpz_yyy_zz_0, tpz_yyz_xx_0, \
                                     tpz_yyz_xy_0, tpz_yyz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xxzz_xx_0[j] = pa_x[j] * tlx_xzz_xx_0[j] + 0.5 * fl1_fx * tlx_zz_xx_0[j] + fl1_fx * tlx_xzz_x_0[j];

                tly_xxzz_xx_0[j] = pa_x[j] * tly_xzz_xx_0[j] + 0.5 * fl1_fx * tly_zz_xx_0[j] + fl1_fx * tly_xzz_x_0[j] + 0.5 * fl1_fx * tpz_xzz_xx_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xx_0[j];

                tlz_xxzz_xx_0[j] = pa_x[j] * tlz_xzz_xx_0[j] + 0.5 * fl1_fx * tlz_zz_xx_0[j] + fl1_fx * tlz_xzz_x_0[j] - 0.5 * fl1_fx * tpy_xzz_xx_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xx_0[j];

                tlx_xxzz_xy_0[j] = pa_x[j] * tlx_xzz_xy_0[j] + 0.5 * fl1_fx * tlx_zz_xy_0[j] + 0.5 * fl1_fx * tlx_xzz_y_0[j];

                tly_xxzz_xy_0[j] = pa_x[j] * tly_xzz_xy_0[j] + 0.5 * fl1_fx * tly_zz_xy_0[j] + 0.5 * fl1_fx * tly_xzz_y_0[j] + 0.5 * fl1_fx * tpz_xzz_xy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xy_0[j];

                tlz_xxzz_xy_0[j] = pa_x[j] * tlz_xzz_xy_0[j] + 0.5 * fl1_fx * tlz_zz_xy_0[j] + 0.5 * fl1_fx * tlz_xzz_y_0[j] - 0.5 * fl1_fx * tpy_xzz_xy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xy_0[j];

                tlx_xxzz_xz_0[j] = pa_x[j] * tlx_xzz_xz_0[j] + 0.5 * fl1_fx * tlx_zz_xz_0[j] + 0.5 * fl1_fx * tlx_xzz_z_0[j];

                tly_xxzz_xz_0[j] = pa_x[j] * tly_xzz_xz_0[j] + 0.5 * fl1_fx * tly_zz_xz_0[j] + 0.5 * fl1_fx * tly_xzz_z_0[j] + 0.5 * fl1_fx * tpz_xzz_xz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xz_0[j];

                tlz_xxzz_xz_0[j] = pa_x[j] * tlz_xzz_xz_0[j] + 0.5 * fl1_fx * tlz_zz_xz_0[j] + 0.5 * fl1_fx * tlz_xzz_z_0[j] - 0.5 * fl1_fx * tpy_xzz_xz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xz_0[j];

                tlx_xxzz_yy_0[j] = pa_x[j] * tlx_xzz_yy_0[j] + 0.5 * fl1_fx * tlx_zz_yy_0[j];

                tly_xxzz_yy_0[j] = pa_x[j] * tly_xzz_yy_0[j] + 0.5 * fl1_fx * tly_zz_yy_0[j] + 0.5 * fl1_fx * tpz_xzz_yy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_yy_0[j];

                tlz_xxzz_yy_0[j] = pa_x[j] * tlz_xzz_yy_0[j] + 0.5 * fl1_fx * tlz_zz_yy_0[j] - 0.5 * fl1_fx * tpy_xzz_yy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_yy_0[j];

                tlx_xxzz_yz_0[j] = pa_x[j] * tlx_xzz_yz_0[j] + 0.5 * fl1_fx * tlx_zz_yz_0[j];

                tly_xxzz_yz_0[j] = pa_x[j] * tly_xzz_yz_0[j] + 0.5 * fl1_fx * tly_zz_yz_0[j] + 0.5 * fl1_fx * tpz_xzz_yz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_yz_0[j];

                tlz_xxzz_yz_0[j] = pa_x[j] * tlz_xzz_yz_0[j] + 0.5 * fl1_fx * tlz_zz_yz_0[j] - 0.5 * fl1_fx * tpy_xzz_yz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_yz_0[j];

                tlx_xxzz_zz_0[j] = pa_x[j] * tlx_xzz_zz_0[j] + 0.5 * fl1_fx * tlx_zz_zz_0[j];

                tly_xxzz_zz_0[j] = pa_x[j] * tly_xzz_zz_0[j] + 0.5 * fl1_fx * tly_zz_zz_0[j] + 0.5 * fl1_fx * tpz_xzz_zz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_zz_0[j];

                tlz_xxzz_zz_0[j] = pa_x[j] * tlz_xzz_zz_0[j] + 0.5 * fl1_fx * tlz_zz_zz_0[j] - 0.5 * fl1_fx * tpy_xzz_zz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_zz_0[j];

                tlx_xyyy_xx_0[j] = pa_x[j] * tlx_yyy_xx_0[j] + fl1_fx * tlx_yyy_x_0[j];

                tly_xyyy_xx_0[j] = pa_x[j] * tly_yyy_xx_0[j] + fl1_fx * tly_yyy_x_0[j] + 0.5 * fl1_fx * tpz_yyy_xx_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xx_0[j];

                tlz_xyyy_xx_0[j] = pa_x[j] * tlz_yyy_xx_0[j] + fl1_fx * tlz_yyy_x_0[j] - 0.5 * fl1_fx * tpy_yyy_xx_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xx_0[j];

                tlx_xyyy_xy_0[j] = pa_x[j] * tlx_yyy_xy_0[j] + 0.5 * fl1_fx * tlx_yyy_y_0[j];

                tly_xyyy_xy_0[j] = pa_x[j] * tly_yyy_xy_0[j] + 0.5 * fl1_fx * tly_yyy_y_0[j] + 0.5 * fl1_fx * tpz_yyy_xy_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xy_0[j];

                tlz_xyyy_xy_0[j] = pa_x[j] * tlz_yyy_xy_0[j] + 0.5 * fl1_fx * tlz_yyy_y_0[j] - 0.5 * fl1_fx * tpy_yyy_xy_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xy_0[j];

                tlx_xyyy_xz_0[j] = pa_x[j] * tlx_yyy_xz_0[j] + 0.5 * fl1_fx * tlx_yyy_z_0[j];

                tly_xyyy_xz_0[j] = pa_x[j] * tly_yyy_xz_0[j] + 0.5 * fl1_fx * tly_yyy_z_0[j] + 0.5 * fl1_fx * tpz_yyy_xz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xz_0[j];

                tlz_xyyy_xz_0[j] = pa_x[j] * tlz_yyy_xz_0[j] + 0.5 * fl1_fx * tlz_yyy_z_0[j] - 0.5 * fl1_fx * tpy_yyy_xz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xz_0[j];

                tlx_xyyy_yy_0[j] = pa_x[j] * tlx_yyy_yy_0[j];

                tly_xyyy_yy_0[j] = pa_x[j] * tly_yyy_yy_0[j] + 0.5 * fl1_fx * tpz_yyy_yy_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yy_0[j];

                tlz_xyyy_yy_0[j] = pa_x[j] * tlz_yyy_yy_0[j] - 0.5 * fl1_fx * tpy_yyy_yy_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yy_0[j];

                tlx_xyyy_yz_0[j] = pa_x[j] * tlx_yyy_yz_0[j];

                tly_xyyy_yz_0[j] = pa_x[j] * tly_yyy_yz_0[j] + 0.5 * fl1_fx * tpz_yyy_yz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yz_0[j];

                tlz_xyyy_yz_0[j] = pa_x[j] * tlz_yyy_yz_0[j] - 0.5 * fl1_fx * tpy_yyy_yz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yz_0[j];

                tlx_xyyy_zz_0[j] = pa_x[j] * tlx_yyy_zz_0[j];

                tly_xyyy_zz_0[j] = pa_x[j] * tly_yyy_zz_0[j] + 0.5 * fl1_fx * tpz_yyy_zz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_zz_0[j];

                tlz_xyyy_zz_0[j] = pa_x[j] * tlz_yyy_zz_0[j] - 0.5 * fl1_fx * tpy_yyy_zz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_zz_0[j];

                tlx_xyyz_xx_0[j] = pa_x[j] * tlx_yyz_xx_0[j] + fl1_fx * tlx_yyz_x_0[j];

                tly_xyyz_xx_0[j] = pa_x[j] * tly_yyz_xx_0[j] + fl1_fx * tly_yyz_x_0[j] + 0.5 * fl1_fx * tpz_yyz_xx_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xx_0[j];

                tlz_xyyz_xx_0[j] = pa_x[j] * tlz_yyz_xx_0[j] + fl1_fx * tlz_yyz_x_0[j] - 0.5 * fl1_fx * tpy_yyz_xx_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xx_0[j];

                tlx_xyyz_xy_0[j] = pa_x[j] * tlx_yyz_xy_0[j] + 0.5 * fl1_fx * tlx_yyz_y_0[j];

                tly_xyyz_xy_0[j] = pa_x[j] * tly_yyz_xy_0[j] + 0.5 * fl1_fx * tly_yyz_y_0[j] + 0.5 * fl1_fx * tpz_yyz_xy_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xy_0[j];

                tlz_xyyz_xy_0[j] = pa_x[j] * tlz_yyz_xy_0[j] + 0.5 * fl1_fx * tlz_yyz_y_0[j] - 0.5 * fl1_fx * tpy_yyz_xy_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xy_0[j];

                tlx_xyyz_xz_0[j] = pa_x[j] * tlx_yyz_xz_0[j] + 0.5 * fl1_fx * tlx_yyz_z_0[j];

                tly_xyyz_xz_0[j] = pa_x[j] * tly_yyz_xz_0[j] + 0.5 * fl1_fx * tly_yyz_z_0[j] + 0.5 * fl1_fx * tpz_yyz_xz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xz_0[j];

                tlz_xyyz_xz_0[j] = pa_x[j] * tlz_yyz_xz_0[j] + 0.5 * fl1_fx * tlz_yyz_z_0[j] - 0.5 * fl1_fx * tpy_yyz_xz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForGD_135_180(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tlx_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 45); 

            auto tly_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tlz_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tlx_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 46); 

            auto tly_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tlz_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tlx_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 47); 

            auto tly_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tlz_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tlx_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 48); 

            auto tly_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tlz_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tlx_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 49); 

            auto tly_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tlz_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tlx_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 50); 

            auto tly_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tlz_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tlx_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 51); 

            auto tly_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tlz_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tlx_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 52); 

            auto tly_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tlz_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tlx_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 53); 

            auto tly_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 53); 

            auto tlz_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tlx_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 54); 

            auto tly_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tlz_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tlx_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 55); 

            auto tly_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tlz_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tlx_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 56); 

            auto tly_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tlz_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tlx_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 57); 

            auto tly_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tlz_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tlx_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 58); 

            auto tly_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tlz_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tlx_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 59); 

            auto tly_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tlz_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            auto tlx_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 24); 

            auto tly_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 24); 

            auto tlx_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 25); 

            auto tly_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 25); 

            auto tlz_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 25); 

            auto tlx_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 26); 

            auto tly_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 26); 

            auto tlz_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 26); 

            auto tlx_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 27); 

            auto tly_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 27); 

            auto tlz_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 27); 

            auto tlx_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 28); 

            auto tly_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 28); 

            auto tlz_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 28); 

            auto tlx_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 29); 

            auto tly_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 29); 

            auto tlz_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 29); 

            auto tpy_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tpz_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tpy_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tpz_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tpy_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tpz_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tpy_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tpz_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tpy_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tpz_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tpy_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tpz_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tpy_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tpz_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tpy_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tpz_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tpy_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 53); 

            auto tpz_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tpy_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tpz_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tpy_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tpz_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tpy_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tpz_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tpy_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tpz_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tpy_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tpz_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            auto tdy_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tdz_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdy_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tdz_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdy_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tdz_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdy_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tdz_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdy_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tdz_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdy_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tdz_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdy_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tdz_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdy_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tdz_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdy_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 53); 

            auto tdz_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdy_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdy_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdy_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdy_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdy_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdy_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            // set up pointers to integrals

            auto tlx_xyyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 45); 

            auto tly_xyyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 45); 

            auto tlz_xyyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 45); 

            auto tlx_xyyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 46); 

            auto tly_xyyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 46); 

            auto tlz_xyyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 46); 

            auto tlx_xyyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 47); 

            auto tly_xyyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 47); 

            auto tlz_xyyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 47); 

            auto tlx_xyzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 48); 

            auto tly_xyzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 48); 

            auto tlz_xyzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 48); 

            auto tlx_xyzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 49); 

            auto tly_xyzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 49); 

            auto tlz_xyzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 49); 

            auto tlx_xyzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 50); 

            auto tly_xyzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 50); 

            auto tlz_xyzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 50); 

            auto tlx_xyzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 51); 

            auto tly_xyzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 51); 

            auto tlz_xyzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 51); 

            auto tlx_xyzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 52); 

            auto tly_xyzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 52); 

            auto tlz_xyzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 52); 

            auto tlx_xyzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 53); 

            auto tly_xyzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 53); 

            auto tlz_xyzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 53); 

            auto tlx_xzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 54); 

            auto tly_xzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 54); 

            auto tlz_xzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 54); 

            auto tlx_xzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 55); 

            auto tly_xzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 55); 

            auto tlz_xzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 55); 

            auto tlx_xzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 56); 

            auto tly_xzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 56); 

            auto tlz_xzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 56); 

            auto tlx_xzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 57); 

            auto tly_xzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 57); 

            auto tlz_xzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 57); 

            auto tlx_xzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 58); 

            auto tly_xzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 58); 

            auto tlz_xzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 58); 

            auto tlx_xzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 59); 

            auto tly_xzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 59); 

            auto tlz_xzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 59); 

            // Batch of Integrals (135,180)

            #pragma omp simd aligned(fgb, fx, pa_x, tdy_yyz_yy_0, tdy_yyz_yz_0, tdy_yyz_zz_0, tdy_yzz_xx_0, \
                                     tdy_yzz_xy_0, tdy_yzz_xz_0, tdy_yzz_yy_0, tdy_yzz_yz_0, tdy_yzz_zz_0, tdy_zzz_xx_0, \
                                     tdy_zzz_xy_0, tdy_zzz_xz_0, tdy_zzz_yy_0, tdy_zzz_yz_0, tdy_zzz_zz_0, tdz_yyz_yy_0, \
                                     tdz_yyz_yz_0, tdz_yyz_zz_0, tdz_yzz_xx_0, tdz_yzz_xy_0, tdz_yzz_xz_0, tdz_yzz_yy_0, \
                                     tdz_yzz_yz_0, tdz_yzz_zz_0, tdz_zzz_xx_0, tdz_zzz_xy_0, tdz_zzz_xz_0, tdz_zzz_yy_0, \
                                     tdz_zzz_yz_0, tdz_zzz_zz_0, tlx_xyyz_yy_0, tlx_xyyz_yz_0, tlx_xyyz_zz_0, \
                                     tlx_xyzz_xx_0, tlx_xyzz_xy_0, tlx_xyzz_xz_0, tlx_xyzz_yy_0, tlx_xyzz_yz_0, \
                                     tlx_xyzz_zz_0, tlx_xzzz_xx_0, tlx_xzzz_xy_0, tlx_xzzz_xz_0, tlx_xzzz_yy_0, \
                                     tlx_xzzz_yz_0, tlx_xzzz_zz_0, tlx_yyz_yy_0, tlx_yyz_yz_0, tlx_yyz_zz_0, tlx_yzz_x_0, \
                                     tlx_yzz_xx_0, tlx_yzz_xy_0, tlx_yzz_xz_0, tlx_yzz_y_0, tlx_yzz_yy_0, tlx_yzz_yz_0, \
                                     tlx_yzz_z_0, tlx_yzz_zz_0, tlx_zzz_x_0, tlx_zzz_xx_0, tlx_zzz_xy_0, tlx_zzz_xz_0, \
                                     tlx_zzz_y_0, tlx_zzz_yy_0, tlx_zzz_yz_0, tlx_zzz_z_0, tlx_zzz_zz_0, tly_xyyz_yy_0, \
                                     tly_xyyz_yz_0, tly_xyyz_zz_0, tly_xyzz_xx_0, tly_xyzz_xy_0, tly_xyzz_xz_0, \
                                     tly_xyzz_yy_0, tly_xyzz_yz_0, tly_xyzz_zz_0, tly_xzzz_xx_0, tly_xzzz_xy_0, \
                                     tly_xzzz_xz_0, tly_xzzz_yy_0, tly_xzzz_yz_0, tly_xzzz_zz_0, tly_yyz_yy_0, \
                                     tly_yyz_yz_0, tly_yyz_zz_0, tly_yzz_x_0, tly_yzz_xx_0, tly_yzz_xy_0, tly_yzz_xz_0, \
                                     tly_yzz_y_0, tly_yzz_yy_0, tly_yzz_yz_0, tly_yzz_z_0, tly_yzz_zz_0, tly_zzz_x_0, \
                                     tly_zzz_xx_0, tly_zzz_xy_0, tly_zzz_xz_0, tly_zzz_y_0, tly_zzz_yy_0, tly_zzz_yz_0, \
                                     tly_zzz_z_0, tly_zzz_zz_0, tlz_xyyz_yy_0, tlz_xyyz_yz_0, tlz_xyyz_zz_0, \
                                     tlz_xyzz_xx_0, tlz_xyzz_xy_0, tlz_xyzz_xz_0, tlz_xyzz_yy_0, tlz_xyzz_yz_0, \
                                     tlz_xyzz_zz_0, tlz_xzzz_xx_0, tlz_xzzz_xy_0, tlz_xzzz_xz_0, tlz_xzzz_yy_0, \
                                     tlz_xzzz_yz_0, tlz_xzzz_zz_0, tlz_yyz_yy_0, tlz_yyz_yz_0, tlz_yyz_zz_0, tlz_yzz_x_0, \
                                     tlz_yzz_xx_0, tlz_yzz_xy_0, tlz_yzz_xz_0, tlz_yzz_y_0, tlz_yzz_yy_0, tlz_yzz_yz_0, \
                                     tlz_yzz_z_0, tlz_yzz_zz_0, tlz_zzz_x_0, tlz_zzz_xx_0, tlz_zzz_xy_0, tlz_zzz_xz_0, \
                                     tlz_zzz_y_0, tlz_zzz_yy_0, tlz_zzz_yz_0, tlz_zzz_z_0, tlz_zzz_zz_0, tpy_yyz_yy_0, \
                                     tpy_yyz_yz_0, tpy_yyz_zz_0, tpy_yzz_xx_0, tpy_yzz_xy_0, tpy_yzz_xz_0, tpy_yzz_yy_0, \
                                     tpy_yzz_yz_0, tpy_yzz_zz_0, tpy_zzz_xx_0, tpy_zzz_xy_0, tpy_zzz_xz_0, tpy_zzz_yy_0, \
                                     tpy_zzz_yz_0, tpy_zzz_zz_0, tpz_yyz_yy_0, tpz_yyz_yz_0, tpz_yyz_zz_0, tpz_yzz_xx_0, \
                                     tpz_yzz_xy_0, tpz_yzz_xz_0, tpz_yzz_yy_0, tpz_yzz_yz_0, tpz_yzz_zz_0, tpz_zzz_xx_0, \
                                     tpz_zzz_xy_0, tpz_zzz_xz_0, tpz_zzz_yy_0, tpz_zzz_yz_0, tpz_zzz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_xyyz_yy_0[j] = pa_x[j] * tlx_yyz_yy_0[j];

                tly_xyyz_yy_0[j] = pa_x[j] * tly_yyz_yy_0[j] + 0.5 * fl1_fx * tpz_yyz_yy_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yy_0[j];

                tlz_xyyz_yy_0[j] = pa_x[j] * tlz_yyz_yy_0[j] - 0.5 * fl1_fx * tpy_yyz_yy_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yy_0[j];

                tlx_xyyz_yz_0[j] = pa_x[j] * tlx_yyz_yz_0[j];

                tly_xyyz_yz_0[j] = pa_x[j] * tly_yyz_yz_0[j] + 0.5 * fl1_fx * tpz_yyz_yz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yz_0[j];

                tlz_xyyz_yz_0[j] = pa_x[j] * tlz_yyz_yz_0[j] - 0.5 * fl1_fx * tpy_yyz_yz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yz_0[j];

                tlx_xyyz_zz_0[j] = pa_x[j] * tlx_yyz_zz_0[j];

                tly_xyyz_zz_0[j] = pa_x[j] * tly_yyz_zz_0[j] + 0.5 * fl1_fx * tpz_yyz_zz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_zz_0[j];

                tlz_xyyz_zz_0[j] = pa_x[j] * tlz_yyz_zz_0[j] - 0.5 * fl1_fx * tpy_yyz_zz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_zz_0[j];

                tlx_xyzz_xx_0[j] = pa_x[j] * tlx_yzz_xx_0[j] + fl1_fx * tlx_yzz_x_0[j];

                tly_xyzz_xx_0[j] = pa_x[j] * tly_yzz_xx_0[j] + fl1_fx * tly_yzz_x_0[j] + 0.5 * fl1_fx * tpz_yzz_xx_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xx_0[j];

                tlz_xyzz_xx_0[j] = pa_x[j] * tlz_yzz_xx_0[j] + fl1_fx * tlz_yzz_x_0[j] - 0.5 * fl1_fx * tpy_yzz_xx_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xx_0[j];

                tlx_xyzz_xy_0[j] = pa_x[j] * tlx_yzz_xy_0[j] + 0.5 * fl1_fx * tlx_yzz_y_0[j];

                tly_xyzz_xy_0[j] = pa_x[j] * tly_yzz_xy_0[j] + 0.5 * fl1_fx * tly_yzz_y_0[j] + 0.5 * fl1_fx * tpz_yzz_xy_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xy_0[j];

                tlz_xyzz_xy_0[j] = pa_x[j] * tlz_yzz_xy_0[j] + 0.5 * fl1_fx * tlz_yzz_y_0[j] - 0.5 * fl1_fx * tpy_yzz_xy_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xy_0[j];

                tlx_xyzz_xz_0[j] = pa_x[j] * tlx_yzz_xz_0[j] + 0.5 * fl1_fx * tlx_yzz_z_0[j];

                tly_xyzz_xz_0[j] = pa_x[j] * tly_yzz_xz_0[j] + 0.5 * fl1_fx * tly_yzz_z_0[j] + 0.5 * fl1_fx * tpz_yzz_xz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xz_0[j];

                tlz_xyzz_xz_0[j] = pa_x[j] * tlz_yzz_xz_0[j] + 0.5 * fl1_fx * tlz_yzz_z_0[j] - 0.5 * fl1_fx * tpy_yzz_xz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xz_0[j];

                tlx_xyzz_yy_0[j] = pa_x[j] * tlx_yzz_yy_0[j];

                tly_xyzz_yy_0[j] = pa_x[j] * tly_yzz_yy_0[j] + 0.5 * fl1_fx * tpz_yzz_yy_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yy_0[j];

                tlz_xyzz_yy_0[j] = pa_x[j] * tlz_yzz_yy_0[j] - 0.5 * fl1_fx * tpy_yzz_yy_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yy_0[j];

                tlx_xyzz_yz_0[j] = pa_x[j] * tlx_yzz_yz_0[j];

                tly_xyzz_yz_0[j] = pa_x[j] * tly_yzz_yz_0[j] + 0.5 * fl1_fx * tpz_yzz_yz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yz_0[j];

                tlz_xyzz_yz_0[j] = pa_x[j] * tlz_yzz_yz_0[j] - 0.5 * fl1_fx * tpy_yzz_yz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yz_0[j];

                tlx_xyzz_zz_0[j] = pa_x[j] * tlx_yzz_zz_0[j];

                tly_xyzz_zz_0[j] = pa_x[j] * tly_yzz_zz_0[j] + 0.5 * fl1_fx * tpz_yzz_zz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_zz_0[j];

                tlz_xyzz_zz_0[j] = pa_x[j] * tlz_yzz_zz_0[j] - 0.5 * fl1_fx * tpy_yzz_zz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_zz_0[j];

                tlx_xzzz_xx_0[j] = pa_x[j] * tlx_zzz_xx_0[j] + fl1_fx * tlx_zzz_x_0[j];

                tly_xzzz_xx_0[j] = pa_x[j] * tly_zzz_xx_0[j] + fl1_fx * tly_zzz_x_0[j] + 0.5 * fl1_fx * tpz_zzz_xx_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xx_0[j];

                tlz_xzzz_xx_0[j] = pa_x[j] * tlz_zzz_xx_0[j] + fl1_fx * tlz_zzz_x_0[j] - 0.5 * fl1_fx * tpy_zzz_xx_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xx_0[j];

                tlx_xzzz_xy_0[j] = pa_x[j] * tlx_zzz_xy_0[j] + 0.5 * fl1_fx * tlx_zzz_y_0[j];

                tly_xzzz_xy_0[j] = pa_x[j] * tly_zzz_xy_0[j] + 0.5 * fl1_fx * tly_zzz_y_0[j] + 0.5 * fl1_fx * tpz_zzz_xy_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xy_0[j];

                tlz_xzzz_xy_0[j] = pa_x[j] * tlz_zzz_xy_0[j] + 0.5 * fl1_fx * tlz_zzz_y_0[j] - 0.5 * fl1_fx * tpy_zzz_xy_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xy_0[j];

                tlx_xzzz_xz_0[j] = pa_x[j] * tlx_zzz_xz_0[j] + 0.5 * fl1_fx * tlx_zzz_z_0[j];

                tly_xzzz_xz_0[j] = pa_x[j] * tly_zzz_xz_0[j] + 0.5 * fl1_fx * tly_zzz_z_0[j] + 0.5 * fl1_fx * tpz_zzz_xz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xz_0[j];

                tlz_xzzz_xz_0[j] = pa_x[j] * tlz_zzz_xz_0[j] + 0.5 * fl1_fx * tlz_zzz_z_0[j] - 0.5 * fl1_fx * tpy_zzz_xz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xz_0[j];

                tlx_xzzz_yy_0[j] = pa_x[j] * tlx_zzz_yy_0[j];

                tly_xzzz_yy_0[j] = pa_x[j] * tly_zzz_yy_0[j] + 0.5 * fl1_fx * tpz_zzz_yy_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yy_0[j];

                tlz_xzzz_yy_0[j] = pa_x[j] * tlz_zzz_yy_0[j] - 0.5 * fl1_fx * tpy_zzz_yy_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yy_0[j];

                tlx_xzzz_yz_0[j] = pa_x[j] * tlx_zzz_yz_0[j];

                tly_xzzz_yz_0[j] = pa_x[j] * tly_zzz_yz_0[j] + 0.5 * fl1_fx * tpz_zzz_yz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yz_0[j];

                tlz_xzzz_yz_0[j] = pa_x[j] * tlz_zzz_yz_0[j] - 0.5 * fl1_fx * tpy_zzz_yz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yz_0[j];

                tlx_xzzz_zz_0[j] = pa_x[j] * tlx_zzz_zz_0[j];

                tly_xzzz_zz_0[j] = pa_x[j] * tly_zzz_zz_0[j] + 0.5 * fl1_fx * tpz_zzz_zz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_zz_0[j];

                tlz_xzzz_zz_0[j] = pa_x[j] * tlz_zzz_zz_0[j] - 0.5 * fl1_fx * tpy_zzz_zz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_zz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForGD_180_225(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tlx_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 36); 

            auto tly_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 36); 

            auto tlz_yyy_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tlx_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 37); 

            auto tly_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 37); 

            auto tlz_yyy_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tlx_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 38); 

            auto tly_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 38); 

            auto tlz_yyy_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tlx_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 39); 

            auto tly_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 39); 

            auto tlz_yyy_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tlx_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 40); 

            auto tly_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 40); 

            auto tlz_yyy_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tlx_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 41); 

            auto tly_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 41); 

            auto tlz_yyy_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tlx_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 42); 

            auto tly_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 42); 

            auto tlz_yyz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tlx_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 43); 

            auto tly_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 43); 

            auto tlz_yyz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tlx_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 44); 

            auto tly_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 44); 

            auto tlz_yyz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tlx_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 45); 

            auto tly_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 45); 

            auto tlz_yyz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tlx_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 46); 

            auto tly_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 46); 

            auto tlz_yyz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tlx_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 47); 

            auto tly_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 47); 

            auto tlz_yyz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tlx_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 48); 

            auto tly_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 48); 

            auto tlz_yzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tlx_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 49); 

            auto tly_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 49); 

            auto tlz_yzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tlx_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 50); 

            auto tly_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 50); 

            auto tlz_yzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tlx_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 18); 

            auto tly_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 18); 

            auto tlz_yy_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 18); 

            auto tlx_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 19); 

            auto tly_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 19); 

            auto tlz_yy_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 19); 

            auto tlx_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 20); 

            auto tly_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 20); 

            auto tlz_yy_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 20); 

            auto tlx_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 21); 

            auto tly_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 21); 

            auto tlz_yy_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 21); 

            auto tlx_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 22); 

            auto tly_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 22); 

            auto tlz_yy_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 22); 

            auto tlx_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 23); 

            auto tly_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 23); 

            auto tlz_yy_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 23); 

            auto tlx_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 24); 

            auto tly_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 24); 

            auto tlz_yz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 24); 

            auto tlx_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 25); 

            auto tly_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 25); 

            auto tlz_yz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 25); 

            auto tlx_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 26); 

            auto tly_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 26); 

            auto tlz_yz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 26); 

            auto tlx_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 27); 

            auto tly_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 27); 

            auto tlz_yz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 27); 

            auto tlx_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 28); 

            auto tly_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 28); 

            auto tlz_yz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 28); 

            auto tlx_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 29); 

            auto tly_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 29); 

            auto tlz_yz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 29); 

            auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30); 

            auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31); 

            auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32); 

            auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tlx_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 18); 

            auto tly_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 18); 

            auto tlz_yyy_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 18); 

            auto tlx_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 19); 

            auto tly_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 19); 

            auto tlz_yyy_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 19); 

            auto tlx_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 20); 

            auto tly_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 20); 

            auto tlz_yyy_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 20); 

            auto tlx_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 21); 

            auto tly_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 21); 

            auto tlz_yyz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 21); 

            auto tlx_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 22); 

            auto tly_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 22); 

            auto tlz_yyz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 22); 

            auto tlx_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 23); 

            auto tly_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 23); 

            auto tlz_yyz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 23); 

            auto tlx_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 24); 

            auto tly_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 24); 

            auto tlz_yzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 24); 

            auto tpx_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 36); 

            auto tpz_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tpx_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 37); 

            auto tpz_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tpx_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 38); 

            auto tpz_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tpx_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 39); 

            auto tpz_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tpx_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 40); 

            auto tpz_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tpx_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 41); 

            auto tpz_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tpx_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 42); 

            auto tpz_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tpx_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 43); 

            auto tpz_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tpx_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 44); 

            auto tpz_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tpx_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 45); 

            auto tpz_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tpx_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 46); 

            auto tpz_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tpx_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 47); 

            auto tpz_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tpx_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 48); 

            auto tpz_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tpx_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 49); 

            auto tpz_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tpx_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 50); 

            auto tpz_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            auto tdx_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 36); 

            auto tdz_yyy_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 36); 

            auto tdx_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 37); 

            auto tdz_yyy_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 37); 

            auto tdx_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 38); 

            auto tdz_yyy_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 38); 

            auto tdx_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 39); 

            auto tdz_yyy_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 39); 

            auto tdx_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 40); 

            auto tdz_yyy_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 40); 

            auto tdx_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 41); 

            auto tdz_yyy_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 41); 

            auto tdx_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 42); 

            auto tdz_yyz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 42); 

            auto tdx_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 43); 

            auto tdz_yyz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 43); 

            auto tdx_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 44); 

            auto tdz_yyz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 44); 

            auto tdx_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 45); 

            auto tdz_yyz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 45); 

            auto tdx_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 46); 

            auto tdz_yyz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 46); 

            auto tdx_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 47); 

            auto tdz_yyz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 47); 

            auto tdx_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 48); 

            auto tdz_yzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 48); 

            auto tdx_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 49); 

            auto tdz_yzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 49); 

            auto tdx_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 50); 

            auto tdz_yzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 50); 

            // set up pointers to integrals

            auto tlx_yyyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 60); 

            auto tly_yyyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 60); 

            auto tlz_yyyy_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 60); 

            auto tlx_yyyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 61); 

            auto tly_yyyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 61); 

            auto tlz_yyyy_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 61); 

            auto tlx_yyyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 62); 

            auto tly_yyyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 62); 

            auto tlz_yyyy_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 62); 

            auto tlx_yyyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 63); 

            auto tly_yyyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 63); 

            auto tlz_yyyy_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 63); 

            auto tlx_yyyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 64); 

            auto tly_yyyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 64); 

            auto tlz_yyyy_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 64); 

            auto tlx_yyyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 65); 

            auto tly_yyyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 65); 

            auto tlz_yyyy_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 65); 

            auto tlx_yyyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 66); 

            auto tly_yyyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 66); 

            auto tlz_yyyz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 66); 

            auto tlx_yyyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 67); 

            auto tly_yyyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 67); 

            auto tlz_yyyz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 67); 

            auto tlx_yyyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 68); 

            auto tly_yyyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 68); 

            auto tlz_yyyz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 68); 

            auto tlx_yyyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 69); 

            auto tly_yyyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 69); 

            auto tlz_yyyz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 69); 

            auto tlx_yyyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 70); 

            auto tly_yyyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 70); 

            auto tlz_yyyz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 70); 

            auto tlx_yyyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 71); 

            auto tly_yyyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 71); 

            auto tlz_yyyz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 71); 

            auto tlx_yyzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 72); 

            auto tly_yyzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 72); 

            auto tlz_yyzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 72); 

            auto tlx_yyzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 73); 

            auto tly_yyzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 73); 

            auto tlz_yyzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 73); 

            auto tlx_yyzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 74); 

            auto tly_yyzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 74); 

            auto tlz_yyzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 74); 

            // Batch of Integrals (180,225)

            #pragma omp simd aligned(fgb, fx, pa_y, tdx_yyy_xx_0, tdx_yyy_xy_0, tdx_yyy_xz_0, tdx_yyy_yy_0, \
                                     tdx_yyy_yz_0, tdx_yyy_zz_0, tdx_yyz_xx_0, tdx_yyz_xy_0, tdx_yyz_xz_0, tdx_yyz_yy_0, \
                                     tdx_yyz_yz_0, tdx_yyz_zz_0, tdx_yzz_xx_0, tdx_yzz_xy_0, tdx_yzz_xz_0, tdz_yyy_xx_0, \
                                     tdz_yyy_xy_0, tdz_yyy_xz_0, tdz_yyy_yy_0, tdz_yyy_yz_0, tdz_yyy_zz_0, tdz_yyz_xx_0, \
                                     tdz_yyz_xy_0, tdz_yyz_xz_0, tdz_yyz_yy_0, tdz_yyz_yz_0, tdz_yyz_zz_0, tdz_yzz_xx_0, \
                                     tdz_yzz_xy_0, tdz_yzz_xz_0, tlx_yy_xx_0, tlx_yy_xy_0, tlx_yy_xz_0, tlx_yy_yy_0, \
                                     tlx_yy_yz_0, tlx_yy_zz_0, tlx_yyy_x_0, tlx_yyy_xx_0, tlx_yyy_xy_0, tlx_yyy_xz_0, \
                                     tlx_yyy_y_0, tlx_yyy_yy_0, tlx_yyy_yz_0, tlx_yyy_z_0, tlx_yyy_zz_0, tlx_yyyy_xx_0, \
                                     tlx_yyyy_xy_0, tlx_yyyy_xz_0, tlx_yyyy_yy_0, tlx_yyyy_yz_0, tlx_yyyy_zz_0, \
                                     tlx_yyyz_xx_0, tlx_yyyz_xy_0, tlx_yyyz_xz_0, tlx_yyyz_yy_0, tlx_yyyz_yz_0, \
                                     tlx_yyyz_zz_0, tlx_yyz_x_0, tlx_yyz_xx_0, tlx_yyz_xy_0, tlx_yyz_xz_0, tlx_yyz_y_0, \
                                     tlx_yyz_yy_0, tlx_yyz_yz_0, tlx_yyz_z_0, tlx_yyz_zz_0, tlx_yyzz_xx_0, \
                                     tlx_yyzz_xy_0, tlx_yyzz_xz_0, tlx_yz_xx_0, tlx_yz_xy_0, tlx_yz_xz_0, tlx_yz_yy_0, \
                                     tlx_yz_yz_0, tlx_yz_zz_0, tlx_yzz_x_0, tlx_yzz_xx_0, tlx_yzz_xy_0, tlx_yzz_xz_0, \
                                     tlx_zz_xx_0, tlx_zz_xy_0, tlx_zz_xz_0, tly_yy_xx_0, tly_yy_xy_0, tly_yy_xz_0, \
                                     tly_yy_yy_0, tly_yy_yz_0, tly_yy_zz_0, tly_yyy_x_0, tly_yyy_xx_0, tly_yyy_xy_0, \
                                     tly_yyy_xz_0, tly_yyy_y_0, tly_yyy_yy_0, tly_yyy_yz_0, tly_yyy_z_0, tly_yyy_zz_0, \
                                     tly_yyyy_xx_0, tly_yyyy_xy_0, tly_yyyy_xz_0, tly_yyyy_yy_0, tly_yyyy_yz_0, \
                                     tly_yyyy_zz_0, tly_yyyz_xx_0, tly_yyyz_xy_0, tly_yyyz_xz_0, tly_yyyz_yy_0, \
                                     tly_yyyz_yz_0, tly_yyyz_zz_0, tly_yyz_x_0, tly_yyz_xx_0, tly_yyz_xy_0, tly_yyz_xz_0, \
                                     tly_yyz_y_0, tly_yyz_yy_0, tly_yyz_yz_0, tly_yyz_z_0, tly_yyz_zz_0, tly_yyzz_xx_0, \
                                     tly_yyzz_xy_0, tly_yyzz_xz_0, tly_yz_xx_0, tly_yz_xy_0, tly_yz_xz_0, tly_yz_yy_0, \
                                     tly_yz_yz_0, tly_yz_zz_0, tly_yzz_x_0, tly_yzz_xx_0, tly_yzz_xy_0, tly_yzz_xz_0, \
                                     tly_zz_xx_0, tly_zz_xy_0, tly_zz_xz_0, tlz_yy_xx_0, tlz_yy_xy_0, tlz_yy_xz_0, \
                                     tlz_yy_yy_0, tlz_yy_yz_0, tlz_yy_zz_0, tlz_yyy_x_0, tlz_yyy_xx_0, tlz_yyy_xy_0, \
                                     tlz_yyy_xz_0, tlz_yyy_y_0, tlz_yyy_yy_0, tlz_yyy_yz_0, tlz_yyy_z_0, tlz_yyy_zz_0, \
                                     tlz_yyyy_xx_0, tlz_yyyy_xy_0, tlz_yyyy_xz_0, tlz_yyyy_yy_0, tlz_yyyy_yz_0, \
                                     tlz_yyyy_zz_0, tlz_yyyz_xx_0, tlz_yyyz_xy_0, tlz_yyyz_xz_0, tlz_yyyz_yy_0, \
                                     tlz_yyyz_yz_0, tlz_yyyz_zz_0, tlz_yyz_x_0, tlz_yyz_xx_0, tlz_yyz_xy_0, tlz_yyz_xz_0, \
                                     tlz_yyz_y_0, tlz_yyz_yy_0, tlz_yyz_yz_0, tlz_yyz_z_0, tlz_yyz_zz_0, tlz_yyzz_xx_0, \
                                     tlz_yyzz_xy_0, tlz_yyzz_xz_0, tlz_yz_xx_0, tlz_yz_xy_0, tlz_yz_xz_0, tlz_yz_yy_0, \
                                     tlz_yz_yz_0, tlz_yz_zz_0, tlz_yzz_x_0, tlz_yzz_xx_0, tlz_yzz_xy_0, tlz_yzz_xz_0, \
                                     tlz_zz_xx_0, tlz_zz_xy_0, tlz_zz_xz_0, tpx_yyy_xx_0, tpx_yyy_xy_0, tpx_yyy_xz_0, \
                                     tpx_yyy_yy_0, tpx_yyy_yz_0, tpx_yyy_zz_0, tpx_yyz_xx_0, tpx_yyz_xy_0, tpx_yyz_xz_0, \
                                     tpx_yyz_yy_0, tpx_yyz_yz_0, tpx_yyz_zz_0, tpx_yzz_xx_0, tpx_yzz_xy_0, tpx_yzz_xz_0, \
                                     tpz_yyy_xx_0, tpz_yyy_xy_0, tpz_yyy_xz_0, tpz_yyy_yy_0, tpz_yyy_yz_0, tpz_yyy_zz_0, \
                                     tpz_yyz_xx_0, tpz_yyz_xy_0, tpz_yyz_xz_0, tpz_yyz_yy_0, tpz_yyz_yz_0, tpz_yyz_zz_0, \
                                     tpz_yzz_xx_0, tpz_yzz_xy_0, tpz_yzz_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yyyy_xx_0[j] = pa_y[j] * tlx_yyy_xx_0[j] + 1.5 * fl1_fx * tlx_yy_xx_0[j] - 0.5 * fl1_fx * tpz_yyy_xx_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xx_0[j];

                tly_yyyy_xx_0[j] = pa_y[j] * tly_yyy_xx_0[j] + 1.5 * fl1_fx * tly_yy_xx_0[j];

                tlz_yyyy_xx_0[j] = pa_y[j] * tlz_yyy_xx_0[j] + 1.5 * fl1_fx * tlz_yy_xx_0[j] + 0.5 * fl1_fx * tpx_yyy_xx_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xx_0[j];

                tlx_yyyy_xy_0[j] = pa_y[j] * tlx_yyy_xy_0[j] + 1.5 * fl1_fx * tlx_yy_xy_0[j] + 0.5 * fl1_fx * tlx_yyy_x_0[j] - 0.5 * fl1_fx * tpz_yyy_xy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xy_0[j];

                tly_yyyy_xy_0[j] = pa_y[j] * tly_yyy_xy_0[j] + 1.5 * fl1_fx * tly_yy_xy_0[j] + 0.5 * fl1_fx * tly_yyy_x_0[j];

                tlz_yyyy_xy_0[j] = pa_y[j] * tlz_yyy_xy_0[j] + 1.5 * fl1_fx * tlz_yy_xy_0[j] + 0.5 * fl1_fx * tlz_yyy_x_0[j] + 0.5 * fl1_fx * tpx_yyy_xy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xy_0[j];

                tlx_yyyy_xz_0[j] = pa_y[j] * tlx_yyy_xz_0[j] + 1.5 * fl1_fx * tlx_yy_xz_0[j] - 0.5 * fl1_fx * tpz_yyy_xz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xz_0[j];

                tly_yyyy_xz_0[j] = pa_y[j] * tly_yyy_xz_0[j] + 1.5 * fl1_fx * tly_yy_xz_0[j];

                tlz_yyyy_xz_0[j] = pa_y[j] * tlz_yyy_xz_0[j] + 1.5 * fl1_fx * tlz_yy_xz_0[j] + 0.5 * fl1_fx * tpx_yyy_xz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xz_0[j];

                tlx_yyyy_yy_0[j] = pa_y[j] * tlx_yyy_yy_0[j] + 1.5 * fl1_fx * tlx_yy_yy_0[j] + fl1_fx * tlx_yyy_y_0[j] - 0.5 * fl1_fx * tpz_yyy_yy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yy_0[j];

                tly_yyyy_yy_0[j] = pa_y[j] * tly_yyy_yy_0[j] + 1.5 * fl1_fx * tly_yy_yy_0[j] + fl1_fx * tly_yyy_y_0[j];

                tlz_yyyy_yy_0[j] = pa_y[j] * tlz_yyy_yy_0[j] + 1.5 * fl1_fx * tlz_yy_yy_0[j] + fl1_fx * tlz_yyy_y_0[j] + 0.5 * fl1_fx * tpx_yyy_yy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yy_0[j];

                tlx_yyyy_yz_0[j] = pa_y[j] * tlx_yyy_yz_0[j] + 1.5 * fl1_fx * tlx_yy_yz_0[j] + 0.5 * fl1_fx * tlx_yyy_z_0[j] - 0.5 * fl1_fx * tpz_yyy_yz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yz_0[j];

                tly_yyyy_yz_0[j] = pa_y[j] * tly_yyy_yz_0[j] + 1.5 * fl1_fx * tly_yy_yz_0[j] + 0.5 * fl1_fx * tly_yyy_z_0[j];

                tlz_yyyy_yz_0[j] = pa_y[j] * tlz_yyy_yz_0[j] + 1.5 * fl1_fx * tlz_yy_yz_0[j] + 0.5 * fl1_fx * tlz_yyy_z_0[j] + 0.5 * fl1_fx * tpx_yyy_yz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yz_0[j];

                tlx_yyyy_zz_0[j] = pa_y[j] * tlx_yyy_zz_0[j] + 1.5 * fl1_fx * tlx_yy_zz_0[j] - 0.5 * fl1_fx * tpz_yyy_zz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_zz_0[j];

                tly_yyyy_zz_0[j] = pa_y[j] * tly_yyy_zz_0[j] + 1.5 * fl1_fx * tly_yy_zz_0[j];

                tlz_yyyy_zz_0[j] = pa_y[j] * tlz_yyy_zz_0[j] + 1.5 * fl1_fx * tlz_yy_zz_0[j] + 0.5 * fl1_fx * tpx_yyy_zz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_zz_0[j];

                tlx_yyyz_xx_0[j] = pa_y[j] * tlx_yyz_xx_0[j] + fl1_fx * tlx_yz_xx_0[j] - 0.5 * fl1_fx * tpz_yyz_xx_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xx_0[j];

                tly_yyyz_xx_0[j] = pa_y[j] * tly_yyz_xx_0[j] + fl1_fx * tly_yz_xx_0[j];

                tlz_yyyz_xx_0[j] = pa_y[j] * tlz_yyz_xx_0[j] + fl1_fx * tlz_yz_xx_0[j] + 0.5 * fl1_fx * tpx_yyz_xx_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xx_0[j];

                tlx_yyyz_xy_0[j] = pa_y[j] * tlx_yyz_xy_0[j] + fl1_fx * tlx_yz_xy_0[j] + 0.5 * fl1_fx * tlx_yyz_x_0[j] - 0.5 * fl1_fx * tpz_yyz_xy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xy_0[j];

                tly_yyyz_xy_0[j] = pa_y[j] * tly_yyz_xy_0[j] + fl1_fx * tly_yz_xy_0[j] + 0.5 * fl1_fx * tly_yyz_x_0[j];

                tlz_yyyz_xy_0[j] = pa_y[j] * tlz_yyz_xy_0[j] + fl1_fx * tlz_yz_xy_0[j] + 0.5 * fl1_fx * tlz_yyz_x_0[j] + 0.5 * fl1_fx * tpx_yyz_xy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xy_0[j];

                tlx_yyyz_xz_0[j] = pa_y[j] * tlx_yyz_xz_0[j] + fl1_fx * tlx_yz_xz_0[j] - 0.5 * fl1_fx * tpz_yyz_xz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xz_0[j];

                tly_yyyz_xz_0[j] = pa_y[j] * tly_yyz_xz_0[j] + fl1_fx * tly_yz_xz_0[j];

                tlz_yyyz_xz_0[j] = pa_y[j] * tlz_yyz_xz_0[j] + fl1_fx * tlz_yz_xz_0[j] + 0.5 * fl1_fx * tpx_yyz_xz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xz_0[j];

                tlx_yyyz_yy_0[j] = pa_y[j] * tlx_yyz_yy_0[j] + fl1_fx * tlx_yz_yy_0[j] + fl1_fx * tlx_yyz_y_0[j] - 0.5 * fl1_fx * tpz_yyz_yy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yy_0[j];

                tly_yyyz_yy_0[j] = pa_y[j] * tly_yyz_yy_0[j] + fl1_fx * tly_yz_yy_0[j] + fl1_fx * tly_yyz_y_0[j];

                tlz_yyyz_yy_0[j] = pa_y[j] * tlz_yyz_yy_0[j] + fl1_fx * tlz_yz_yy_0[j] + fl1_fx * tlz_yyz_y_0[j] + 0.5 * fl1_fx * tpx_yyz_yy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yy_0[j];

                tlx_yyyz_yz_0[j] = pa_y[j] * tlx_yyz_yz_0[j] + fl1_fx * tlx_yz_yz_0[j] + 0.5 * fl1_fx * tlx_yyz_z_0[j] - 0.5 * fl1_fx * tpz_yyz_yz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yz_0[j];

                tly_yyyz_yz_0[j] = pa_y[j] * tly_yyz_yz_0[j] + fl1_fx * tly_yz_yz_0[j] + 0.5 * fl1_fx * tly_yyz_z_0[j];

                tlz_yyyz_yz_0[j] = pa_y[j] * tlz_yyz_yz_0[j] + fl1_fx * tlz_yz_yz_0[j] + 0.5 * fl1_fx * tlz_yyz_z_0[j] + 0.5 * fl1_fx * tpx_yyz_yz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yz_0[j];

                tlx_yyyz_zz_0[j] = pa_y[j] * tlx_yyz_zz_0[j] + fl1_fx * tlx_yz_zz_0[j] - 0.5 * fl1_fx * tpz_yyz_zz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_zz_0[j];

                tly_yyyz_zz_0[j] = pa_y[j] * tly_yyz_zz_0[j] + fl1_fx * tly_yz_zz_0[j];

                tlz_yyyz_zz_0[j] = pa_y[j] * tlz_yyz_zz_0[j] + fl1_fx * tlz_yz_zz_0[j] + 0.5 * fl1_fx * tpx_yyz_zz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_zz_0[j];

                tlx_yyzz_xx_0[j] = pa_y[j] * tlx_yzz_xx_0[j] + 0.5 * fl1_fx * tlx_zz_xx_0[j] - 0.5 * fl1_fx * tpz_yzz_xx_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xx_0[j];

                tly_yyzz_xx_0[j] = pa_y[j] * tly_yzz_xx_0[j] + 0.5 * fl1_fx * tly_zz_xx_0[j];

                tlz_yyzz_xx_0[j] = pa_y[j] * tlz_yzz_xx_0[j] + 0.5 * fl1_fx * tlz_zz_xx_0[j] + 0.5 * fl1_fx * tpx_yzz_xx_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xx_0[j];

                tlx_yyzz_xy_0[j] = pa_y[j] * tlx_yzz_xy_0[j] + 0.5 * fl1_fx * tlx_zz_xy_0[j] + 0.5 * fl1_fx * tlx_yzz_x_0[j] - 0.5 * fl1_fx * tpz_yzz_xy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xy_0[j];

                tly_yyzz_xy_0[j] = pa_y[j] * tly_yzz_xy_0[j] + 0.5 * fl1_fx * tly_zz_xy_0[j] + 0.5 * fl1_fx * tly_yzz_x_0[j];

                tlz_yyzz_xy_0[j] = pa_y[j] * tlz_yzz_xy_0[j] + 0.5 * fl1_fx * tlz_zz_xy_0[j] + 0.5 * fl1_fx * tlz_yzz_x_0[j] + 0.5 * fl1_fx * tpx_yzz_xy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xy_0[j];

                tlx_yyzz_xz_0[j] = pa_y[j] * tlx_yzz_xz_0[j] + 0.5 * fl1_fx * tlx_zz_xz_0[j] - 0.5 * fl1_fx * tpz_yzz_xz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xz_0[j];

                tly_yyzz_xz_0[j] = pa_y[j] * tly_yzz_xz_0[j] + 0.5 * fl1_fx * tly_zz_xz_0[j];

                tlz_yyzz_xz_0[j] = pa_y[j] * tlz_yzz_xz_0[j] + 0.5 * fl1_fx * tlz_zz_xz_0[j] + 0.5 * fl1_fx * tpx_yzz_xz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xz_0[j];
            }

            idx++;
        }
    }

    void
    compAngularMomentumForGD_225_270(      CMemBlock2D<double>& primBuffer,
                                     const CRecursionMap&       recursionMap,
                                     const CMemBlock2D<double>& osFactors,
                                     const CMemBlock2D<double>& paDistances,
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

        // set up index of integral

        auto pidx_l_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {4, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // check if integral is needed in recursion expansion

        if (pidx_l_4_2_m0 == -1) return;

        // set up indexes of auxilary integral

        auto pidx_l_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {2, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_l_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {1, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        auto pidx_d_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, 
                                                         {3, -1, -1, -1}, {2, -1, -1, -1}, 
                                                         1, 1, 0));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(4 * idx);

            auto fgb = osFactors.data(4 * idx + 3);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tlx_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 51); 

            auto tly_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 51); 

            auto tlz_yzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tlx_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 52); 

            auto tly_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 52); 

            auto tlz_yzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tlx_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 53); 

            auto tly_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 53); 

            auto tlz_yzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tlx_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 54); 

            auto tly_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tlz_zzz_xx_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tlx_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 55); 

            auto tly_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tlz_zzz_xy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tlx_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 56); 

            auto tly_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tlz_zzz_xz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tlx_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 57); 

            auto tly_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tlz_zzz_yy_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tlx_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 58); 

            auto tly_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tlz_zzz_yz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tlx_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * idx + 59); 

            auto tly_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tlz_zzz_zz_0 = primBuffer.data(pidx_l_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            auto tlx_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 30); 

            auto tly_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 30); 

            auto tlz_zz_xx_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 30); 

            auto tlx_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 31); 

            auto tly_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 31); 

            auto tlz_zz_xy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 31); 

            auto tlx_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 32); 

            auto tly_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 32); 

            auto tlz_zz_xz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 32); 

            auto tlx_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 33); 

            auto tly_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 33); 

            auto tlz_zz_yy_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 33); 

            auto tlx_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 34); 

            auto tly_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 34); 

            auto tlz_zz_yz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 34); 

            auto tlx_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * idx + 35); 

            auto tly_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 36 * bdim + 36 * idx + 35); 

            auto tlz_zz_zz_0 = primBuffer.data(pidx_l_2_2_m0 + 72 * bdim + 36 * idx + 35); 

            auto tlx_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 25); 

            auto tly_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 25); 

            auto tlz_yzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 25); 

            auto tlx_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 26); 

            auto tly_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 26); 

            auto tlz_yzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 26); 

            auto tlx_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 27); 

            auto tly_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 27); 

            auto tlz_zzz_x_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 27); 

            auto tlx_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 28); 

            auto tly_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 28); 

            auto tlz_zzz_y_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 28); 

            auto tlx_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * idx + 29); 

            auto tly_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 30 * bdim + 30 * idx + 29); 

            auto tlz_zzz_z_0 = primBuffer.data(pidx_l_3_1_m0 + 60 * bdim + 30 * idx + 29); 

            auto tpx_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 51); 

            auto tpz_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tpx_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 52); 

            auto tpz_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tpx_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 53); 

            auto tpz_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tpx_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 54); 

            auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tpx_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 55); 

            auto tpy_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tpz_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tpx_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 56); 

            auto tpy_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tpz_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tpx_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 57); 

            auto tpy_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tpz_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tpx_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 58); 

            auto tpy_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tpz_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tpx_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 59); 

            auto tpy_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tpz_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            auto tdx_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 51); 

            auto tdz_yzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 51); 

            auto tdx_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 52); 

            auto tdz_yzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 52); 

            auto tdx_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 53); 

            auto tdz_yzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 53); 

            auto tdx_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 54); 

            auto tdy_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 54); 

            auto tdz_zzz_xx_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 54); 

            auto tdx_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 55); 

            auto tdy_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 55); 

            auto tdz_zzz_xy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 55); 

            auto tdx_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 56); 

            auto tdy_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 56); 

            auto tdz_zzz_xz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 56); 

            auto tdx_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 57); 

            auto tdy_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 57); 

            auto tdz_zzz_yy_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 57); 

            auto tdx_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 58); 

            auto tdy_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 58); 

            auto tdz_zzz_yz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 58); 

            auto tdx_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * idx + 59); 

            auto tdy_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 60 * bdim + 60 * idx + 59); 

            auto tdz_zzz_zz_0 = primBuffer.data(pidx_d_3_2_m0 + 120 * bdim + 60 * idx + 59); 

            // set up pointers to integrals

            auto tlx_yyzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 75); 

            auto tly_yyzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 75); 

            auto tlz_yyzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 75); 

            auto tlx_yyzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 76); 

            auto tly_yyzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 76); 

            auto tlz_yyzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 76); 

            auto tlx_yyzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 77); 

            auto tly_yyzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 77); 

            auto tlz_yyzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 77); 

            auto tlx_yzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 78); 

            auto tly_yzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 78); 

            auto tlz_yzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 78); 

            auto tlx_yzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 79); 

            auto tly_yzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 79); 

            auto tlz_yzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 79); 

            auto tlx_yzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 80); 

            auto tly_yzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 80); 

            auto tlz_yzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 80); 

            auto tlx_yzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 81); 

            auto tly_yzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 81); 

            auto tlz_yzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 81); 

            auto tlx_yzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 82); 

            auto tly_yzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 82); 

            auto tlz_yzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 82); 

            auto tlx_yzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 83); 

            auto tly_yzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 83); 

            auto tlz_yzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 83); 

            auto tlx_zzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 84); 

            auto tly_zzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 84); 

            auto tlz_zzzz_xx_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 84); 

            auto tlx_zzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 85); 

            auto tly_zzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 85); 

            auto tlz_zzzz_xy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 85); 

            auto tlx_zzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 86); 

            auto tly_zzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 86); 

            auto tlz_zzzz_xz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 86); 

            auto tlx_zzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 87); 

            auto tly_zzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 87); 

            auto tlz_zzzz_yy_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 87); 

            auto tlx_zzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 88); 

            auto tly_zzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 88); 

            auto tlz_zzzz_yz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 88); 

            auto tlx_zzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * idx + 89); 

            auto tly_zzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 90 * bdim + 90 * idx + 89); 

            auto tlz_zzzz_zz_0 = primBuffer.data(pidx_l_4_2_m0 + 180 * bdim + 90 * idx + 89); 

            // Batch of Integrals (225,270)

            #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_yzz_yy_0, tdx_yzz_yz_0, tdx_yzz_zz_0, \
                                     tdx_zzz_xx_0, tdx_zzz_xy_0, tdx_zzz_xz_0, tdx_zzz_yy_0, tdx_zzz_yz_0, tdx_zzz_zz_0, \
                                     tdy_zzz_xx_0, tdy_zzz_xy_0, tdy_zzz_xz_0, tdy_zzz_yy_0, tdy_zzz_yz_0, tdy_zzz_zz_0, \
                                     tdz_yzz_yy_0, tdz_yzz_yz_0, tdz_yzz_zz_0, tdz_zzz_xx_0, tdz_zzz_xy_0, tdz_zzz_xz_0, \
                                     tdz_zzz_yy_0, tdz_zzz_yz_0, tdz_zzz_zz_0, tlx_yyzz_yy_0, tlx_yyzz_yz_0, \
                                     tlx_yyzz_zz_0, tlx_yzz_y_0, tlx_yzz_yy_0, tlx_yzz_yz_0, tlx_yzz_z_0, tlx_yzz_zz_0, \
                                     tlx_yzzz_xx_0, tlx_yzzz_xy_0, tlx_yzzz_xz_0, tlx_yzzz_yy_0, tlx_yzzz_yz_0, \
                                     tlx_yzzz_zz_0, tlx_zz_xx_0, tlx_zz_xy_0, tlx_zz_xz_0, tlx_zz_yy_0, tlx_zz_yz_0, \
                                     tlx_zz_zz_0, tlx_zzz_x_0, tlx_zzz_xx_0, tlx_zzz_xy_0, tlx_zzz_xz_0, tlx_zzz_y_0, \
                                     tlx_zzz_yy_0, tlx_zzz_yz_0, tlx_zzz_z_0, tlx_zzz_zz_0, tlx_zzzz_xx_0, \
                                     tlx_zzzz_xy_0, tlx_zzzz_xz_0, tlx_zzzz_yy_0, tlx_zzzz_yz_0, tlx_zzzz_zz_0, \
                                     tly_yyzz_yy_0, tly_yyzz_yz_0, tly_yyzz_zz_0, tly_yzz_y_0, tly_yzz_yy_0, \
                                     tly_yzz_yz_0, tly_yzz_z_0, tly_yzz_zz_0, tly_yzzz_xx_0, tly_yzzz_xy_0, \
                                     tly_yzzz_xz_0, tly_yzzz_yy_0, tly_yzzz_yz_0, tly_yzzz_zz_0, tly_zz_xx_0, \
                                     tly_zz_xy_0, tly_zz_xz_0, tly_zz_yy_0, tly_zz_yz_0, tly_zz_zz_0, tly_zzz_x_0, \
                                     tly_zzz_xx_0, tly_zzz_xy_0, tly_zzz_xz_0, tly_zzz_y_0, tly_zzz_yy_0, tly_zzz_yz_0, \
                                     tly_zzz_z_0, tly_zzz_zz_0, tly_zzzz_xx_0, tly_zzzz_xy_0, tly_zzzz_xz_0, \
                                     tly_zzzz_yy_0, tly_zzzz_yz_0, tly_zzzz_zz_0, tlz_yyzz_yy_0, tlz_yyzz_yz_0, \
                                     tlz_yyzz_zz_0, tlz_yzz_y_0, tlz_yzz_yy_0, tlz_yzz_yz_0, tlz_yzz_z_0, tlz_yzz_zz_0, \
                                     tlz_yzzz_xx_0, tlz_yzzz_xy_0, tlz_yzzz_xz_0, tlz_yzzz_yy_0, tlz_yzzz_yz_0, \
                                     tlz_yzzz_zz_0, tlz_zz_xx_0, tlz_zz_xy_0, tlz_zz_xz_0, tlz_zz_yy_0, tlz_zz_yz_0, \
                                     tlz_zz_zz_0, tlz_zzz_x_0, tlz_zzz_xx_0, tlz_zzz_xy_0, tlz_zzz_xz_0, tlz_zzz_y_0, \
                                     tlz_zzz_yy_0, tlz_zzz_yz_0, tlz_zzz_z_0, tlz_zzz_zz_0, tlz_zzzz_xx_0, \
                                     tlz_zzzz_xy_0, tlz_zzzz_xz_0, tlz_zzzz_yy_0, tlz_zzzz_yz_0, tlz_zzzz_zz_0, \
                                     tpx_yzz_yy_0, tpx_yzz_yz_0, tpx_yzz_zz_0, tpx_zzz_xx_0, tpx_zzz_xy_0, tpx_zzz_xz_0, \
                                     tpx_zzz_yy_0, tpx_zzz_yz_0, tpx_zzz_zz_0, tpy_zzz_xx_0, tpy_zzz_xy_0, tpy_zzz_xz_0, \
                                     tpy_zzz_yy_0, tpy_zzz_yz_0, tpy_zzz_zz_0, tpz_yzz_yy_0, tpz_yzz_yz_0, tpz_yzz_zz_0, \
                                     tpz_zzz_xx_0, tpz_zzz_xy_0, tpz_zzz_xz_0, tpz_zzz_yy_0, tpz_zzz_yz_0, tpz_zzz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fgb = fgb[j];

                double fl1_fx = fx[j];

                tlx_yyzz_yy_0[j] = pa_y[j] * tlx_yzz_yy_0[j] + 0.5 * fl1_fx * tlx_zz_yy_0[j] + fl1_fx * tlx_yzz_y_0[j] - 0.5 * fl1_fx * tpz_yzz_yy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yy_0[j];

                tly_yyzz_yy_0[j] = pa_y[j] * tly_yzz_yy_0[j] + 0.5 * fl1_fx * tly_zz_yy_0[j] + fl1_fx * tly_yzz_y_0[j];

                tlz_yyzz_yy_0[j] = pa_y[j] * tlz_yzz_yy_0[j] + 0.5 * fl1_fx * tlz_zz_yy_0[j] + fl1_fx * tlz_yzz_y_0[j] + 0.5 * fl1_fx * tpx_yzz_yy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yy_0[j];

                tlx_yyzz_yz_0[j] = pa_y[j] * tlx_yzz_yz_0[j] + 0.5 * fl1_fx * tlx_zz_yz_0[j] + 0.5 * fl1_fx * tlx_yzz_z_0[j] - 0.5 * fl1_fx * tpz_yzz_yz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yz_0[j];

                tly_yyzz_yz_0[j] = pa_y[j] * tly_yzz_yz_0[j] + 0.5 * fl1_fx * tly_zz_yz_0[j] + 0.5 * fl1_fx * tly_yzz_z_0[j];

                tlz_yyzz_yz_0[j] = pa_y[j] * tlz_yzz_yz_0[j] + 0.5 * fl1_fx * tlz_zz_yz_0[j] + 0.5 * fl1_fx * tlz_yzz_z_0[j] + 0.5 * fl1_fx * tpx_yzz_yz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yz_0[j];

                tlx_yyzz_zz_0[j] = pa_y[j] * tlx_yzz_zz_0[j] + 0.5 * fl1_fx * tlx_zz_zz_0[j] - 0.5 * fl1_fx * tpz_yzz_zz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_zz_0[j];

                tly_yyzz_zz_0[j] = pa_y[j] * tly_yzz_zz_0[j] + 0.5 * fl1_fx * tly_zz_zz_0[j];

                tlz_yyzz_zz_0[j] = pa_y[j] * tlz_yzz_zz_0[j] + 0.5 * fl1_fx * tlz_zz_zz_0[j] + 0.5 * fl1_fx * tpx_yzz_zz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_zz_0[j];

                tlx_yzzz_xx_0[j] = pa_y[j] * tlx_zzz_xx_0[j] - 0.5 * fl1_fx * tpz_zzz_xx_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xx_0[j];

                tly_yzzz_xx_0[j] = pa_y[j] * tly_zzz_xx_0[j];

                tlz_yzzz_xx_0[j] = pa_y[j] * tlz_zzz_xx_0[j] + 0.5 * fl1_fx * tpx_zzz_xx_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xx_0[j];

                tlx_yzzz_xy_0[j] = pa_y[j] * tlx_zzz_xy_0[j] + 0.5 * fl1_fx * tlx_zzz_x_0[j] - 0.5 * fl1_fx * tpz_zzz_xy_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xy_0[j];

                tly_yzzz_xy_0[j] = pa_y[j] * tly_zzz_xy_0[j] + 0.5 * fl1_fx * tly_zzz_x_0[j];

                tlz_yzzz_xy_0[j] = pa_y[j] * tlz_zzz_xy_0[j] + 0.5 * fl1_fx * tlz_zzz_x_0[j] + 0.5 * fl1_fx * tpx_zzz_xy_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xy_0[j];

                tlx_yzzz_xz_0[j] = pa_y[j] * tlx_zzz_xz_0[j] - 0.5 * fl1_fx * tpz_zzz_xz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xz_0[j];

                tly_yzzz_xz_0[j] = pa_y[j] * tly_zzz_xz_0[j];

                tlz_yzzz_xz_0[j] = pa_y[j] * tlz_zzz_xz_0[j] + 0.5 * fl1_fx * tpx_zzz_xz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xz_0[j];

                tlx_yzzz_yy_0[j] = pa_y[j] * tlx_zzz_yy_0[j] + fl1_fx * tlx_zzz_y_0[j] - 0.5 * fl1_fx * tpz_zzz_yy_0[j] - fl1_fx * fl1_fgb * tdz_zzz_yy_0[j];

                tly_yzzz_yy_0[j] = pa_y[j] * tly_zzz_yy_0[j] + fl1_fx * tly_zzz_y_0[j];

                tlz_yzzz_yy_0[j] = pa_y[j] * tlz_zzz_yy_0[j] + fl1_fx * tlz_zzz_y_0[j] + 0.5 * fl1_fx * tpx_zzz_yy_0[j] + fl1_fx * fl1_fgb * tdx_zzz_yy_0[j];

                tlx_yzzz_yz_0[j] = pa_y[j] * tlx_zzz_yz_0[j] + 0.5 * fl1_fx * tlx_zzz_z_0[j] - 0.5 * fl1_fx * tpz_zzz_yz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_yz_0[j];

                tly_yzzz_yz_0[j] = pa_y[j] * tly_zzz_yz_0[j] + 0.5 * fl1_fx * tly_zzz_z_0[j];

                tlz_yzzz_yz_0[j] = pa_y[j] * tlz_zzz_yz_0[j] + 0.5 * fl1_fx * tlz_zzz_z_0[j] + 0.5 * fl1_fx * tpx_zzz_yz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_yz_0[j];

                tlx_yzzz_zz_0[j] = pa_y[j] * tlx_zzz_zz_0[j] - 0.5 * fl1_fx * tpz_zzz_zz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_zz_0[j];

                tly_yzzz_zz_0[j] = pa_y[j] * tly_zzz_zz_0[j];

                tlz_yzzz_zz_0[j] = pa_y[j] * tlz_zzz_zz_0[j] + 0.5 * fl1_fx * tpx_zzz_zz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_zz_0[j];

                tlx_zzzz_xx_0[j] = pa_z[j] * tlx_zzz_xx_0[j] + 1.5 * fl1_fx * tlx_zz_xx_0[j] + 0.5 * fl1_fx * tpy_zzz_xx_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xx_0[j];

                tly_zzzz_xx_0[j] = pa_z[j] * tly_zzz_xx_0[j] + 1.5 * fl1_fx * tly_zz_xx_0[j] - 0.5 * fl1_fx * tpx_zzz_xx_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xx_0[j];

                tlz_zzzz_xx_0[j] = pa_z[j] * tlz_zzz_xx_0[j] + 1.5 * fl1_fx * tlz_zz_xx_0[j];

                tlx_zzzz_xy_0[j] = pa_z[j] * tlx_zzz_xy_0[j] + 1.5 * fl1_fx * tlx_zz_xy_0[j] + 0.5 * fl1_fx * tpy_zzz_xy_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xy_0[j];

                tly_zzzz_xy_0[j] = pa_z[j] * tly_zzz_xy_0[j] + 1.5 * fl1_fx * tly_zz_xy_0[j] - 0.5 * fl1_fx * tpx_zzz_xy_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xy_0[j];

                tlz_zzzz_xy_0[j] = pa_z[j] * tlz_zzz_xy_0[j] + 1.5 * fl1_fx * tlz_zz_xy_0[j];

                tlx_zzzz_xz_0[j] = pa_z[j] * tlx_zzz_xz_0[j] + 1.5 * fl1_fx * tlx_zz_xz_0[j] + 0.5 * fl1_fx * tlx_zzz_x_0[j] + 0.5 * fl1_fx * tpy_zzz_xz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xz_0[j];

                tly_zzzz_xz_0[j] = pa_z[j] * tly_zzz_xz_0[j] + 1.5 * fl1_fx * tly_zz_xz_0[j] + 0.5 * fl1_fx * tly_zzz_x_0[j] - 0.5 * fl1_fx * tpx_zzz_xz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xz_0[j];

                tlz_zzzz_xz_0[j] = pa_z[j] * tlz_zzz_xz_0[j] + 1.5 * fl1_fx * tlz_zz_xz_0[j] + 0.5 * fl1_fx * tlz_zzz_x_0[j];

                tlx_zzzz_yy_0[j] = pa_z[j] * tlx_zzz_yy_0[j] + 1.5 * fl1_fx * tlx_zz_yy_0[j] + 0.5 * fl1_fx * tpy_zzz_yy_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yy_0[j];

                tly_zzzz_yy_0[j] = pa_z[j] * tly_zzz_yy_0[j] + 1.5 * fl1_fx * tly_zz_yy_0[j] - 0.5 * fl1_fx * tpx_zzz_yy_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yy_0[j];

                tlz_zzzz_yy_0[j] = pa_z[j] * tlz_zzz_yy_0[j] + 1.5 * fl1_fx * tlz_zz_yy_0[j];

                tlx_zzzz_yz_0[j] = pa_z[j] * tlx_zzz_yz_0[j] + 1.5 * fl1_fx * tlx_zz_yz_0[j] + 0.5 * fl1_fx * tlx_zzz_y_0[j] + 0.5 * fl1_fx * tpy_zzz_yz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yz_0[j];

                tly_zzzz_yz_0[j] = pa_z[j] * tly_zzz_yz_0[j] + 1.5 * fl1_fx * tly_zz_yz_0[j] + 0.5 * fl1_fx * tly_zzz_y_0[j] - 0.5 * fl1_fx * tpx_zzz_yz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yz_0[j];

                tlz_zzzz_yz_0[j] = pa_z[j] * tlz_zzz_yz_0[j] + 1.5 * fl1_fx * tlz_zz_yz_0[j] + 0.5 * fl1_fx * tlz_zzz_y_0[j];

                tlx_zzzz_zz_0[j] = pa_z[j] * tlx_zzz_zz_0[j] + 1.5 * fl1_fx * tlx_zz_zz_0[j] + fl1_fx * tlx_zzz_z_0[j] + 0.5 * fl1_fx * tpy_zzz_zz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_zz_0[j];

                tly_zzzz_zz_0[j] = pa_z[j] * tly_zzz_zz_0[j] + 1.5 * fl1_fx * tly_zz_zz_0[j] + fl1_fx * tly_zzz_z_0[j] - 0.5 * fl1_fx * tpx_zzz_zz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_zz_0[j];

                tlz_zzzz_zz_0[j] = pa_z[j] * tlz_zzz_zz_0[j] + 1.5 * fl1_fx * tlz_zz_zz_0[j] + fl1_fx * tlz_zzz_z_0[j];
            }

            idx++;
        }
    }


} // amomrecfunc namespace

